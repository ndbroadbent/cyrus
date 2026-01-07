//! CYTools-compatible intersection number computation.
//!
//! This module implements the exact algorithm from CYTools `_construct_intnum_equations_4d`.
//! It builds a sparse linear system and solves for intersection numbers.

mod distinct;
mod equations;
mod extract;
mod solver;
mod variables;

use crate::Point;
use crate::error::Result;
use crate::intersection::Intersection;
use crate::triangulation::Triangulation;

use distinct::compute_distinct_intnums;
use equations::{build_c_dict, build_eqn_structures};
use extract::extract_intersection_numbers;
use solver::{build_linear_system, solve_sparse_system};
use variables::build_variable_array;

/// Compute intersection numbers using the CYTools algorithm.
///
/// # Arguments
/// * `tri` - Triangulation with simplices
/// * `points` - All points including origin at index 0
/// * `linear_relations_no_origin` - Linear relations matrix WITHOUT origin column
///   Shape: (dim, n_points-1) where column i corresponds to point i+1
///
/// # Returns
/// Intersection tensor indexed by point indices (including origin).
///
/// # Errors
/// Returns an error if the sparse linear system cannot be solved.
pub fn compute_intersection_cytools(
    tri: &Triangulation,
    points: &[Point],
    linear_relations_no_origin: &[Vec<i64>],
) -> Result<Intersection> {
    let n_points = points.len();
    let _dim = if linear_relations_no_origin.is_empty() {
        points[0].dim()
    } else {
        linear_relations_no_origin.len()
    };

    // Step 1: Build pts_ext (points augmented with 1s column) for determinant computation
    // Step 2: Compute distinct intersection numbers from simplices
    let distintnum_array = compute_distinct_intnums(tri, points);

    // Step 3: Build variable array (4-tuples with repeated indices)
    let (variable_array, variable_dict) = build_variable_array(tri, n_points);

    // Step 4: Build c_dict (known 4-form contributions for each 3-tuple)
    let c_dict = build_c_dict(tri, &distintnum_array);

    // Step 5: Build equation array and eqn_dict
    let (eqn_array, eqn_dict) =
        build_eqn_structures(tri, n_points, &variable_array, &variable_dict);

    // Step 6: Build sparse linear system M*x + C = 0
    let (m_triplets, c_vec) = build_linear_system(
        &eqn_array,
        &eqn_dict,
        &c_dict,
        linear_relations_no_origin,
        variable_array.len(),
    );

    // Step 7: Solve the system
    let solution = solve_sparse_system(&m_triplets, &c_vec, variable_array.len())?;

    // Step 8: Extract intersection numbers
    // The variables are 4-form values kappa^V_{ijkl}
    // We need to sum them to get 3-form values kappa_{ijk}
    let kappa =
        extract_intersection_numbers(&variable_array, &solution, &distintnum_array, n_points);

    Ok(kappa)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_build_variable_array_simple() {
        // Simple triangulation with one simplex [0, 1, 2, 3, 4]
        let tri = crate::triangulation::Triangulation::new(vec![vec![0, 1, 2, 3, 4]]);

        let (vars, dict) = build_variable_array(&tri, 5);

        // Should have variables for:
        // - 3-tuples: (1,2,3), (1,2,4), (1,3,4), (2,3,4) -> 3 vars each = 12
        // - 2-tuples: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4) -> 3 vars each = 18
        // - self-intersections: (1,1,1,1), (2,2,2,2), (3,3,3,3), (4,4,4,4) = 4
        // Total: some with overlap

        assert!(!vars.is_empty());
        assert_eq!(vars.len(), dict.len());

        // All variables should have at least one repeated index
        for v in &vars {
            let unique: HashSet<usize> = v.iter().copied().collect();
            assert!(
                unique.len() < 4,
                "Variable {v:?} should have repeated indices"
            );
        }
    }
}
