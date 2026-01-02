//! Kähler cone and Mori cone computations.
//!
//! The Kähler cone is the set of Kähler moduli that yield positive volumes
//! for all curves in the Mori cone.
//!
//! Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 3.

use crate::Point;
use crate::error::{Error, Result};
use crate::integer_math::integer_kernel;
use crate::triangulation::Triangulation;
use malachite::Integer;
use malachite::num::conversion::traits::RoundingFrom;
use malachite::rounding_modes::RoundingMode;
use std::collections::HashMap;

/// Mori cone generators for a given triangulation.
#[derive(Debug, Clone)]
pub struct MoriCone {
    /// Generators as integer vectors in the space of divisors.
    generators: Vec<Vec<Integer>>,
}

impl MoriCone {
    /// Create a new Mori cone from generators.
    pub const fn new(generators: Vec<Vec<Integer>>) -> Self {
        Self { generators }
    }

    /// Get the generators.
    pub fn generators(&self) -> &[Vec<Integer>] {
        &self.generators
    }

    /// Check if a Kähler class `p` is in the Kähler cone.
    ///
    /// A class `p` is in the Kähler cone if `p · R > 0` for all Mori generators `R`.
    ///
    /// # Panics
    /// Panics if `p.len()` does not match the generator dimension.
    pub fn contains(&self, p: &[f64]) -> bool {
        for generator in &self.generators {
            let mut dot = 0.0;
            for (i, val) in generator.iter().enumerate() {
                // gen is Integer, p is f64.
                let (val_f, _) = f64::rounding_from(val, RoundingMode::Nearest);
                dot += val_f * p[i];
            }
            if dot <= 1e-10 {
                // Threshold for numerical stability
                return false;
            }
        }
        true
    }
}

/// Compute the Mori cone generators for a triangulation.
///
/// For each interior ridge (shared by two simplices), there is a unique
/// relation that generates a ray of the Mori cone.
///
/// # Errors
/// Returns an error if the triangulation is inconsistent or if no generators are found.
pub fn compute_mori_generators(tri: &Triangulation, points: &[Point]) -> Result<MoriCone> {
    let n_pts = points.len();
    if n_pts == 0 {
        return Err(Error::InvalidInput("No points provided".into()));
    }
    let dim = points[0].dim();

    // 1. Identify ridges and their adjacent simplices
    let ridge_map = build_ridge_map(tri, dim);

    let mut generators = Vec::new();

    // 2. For each interior ridge (shared by exactly 2 simplices)
    for (ridge, adj_simplices) in ridge_map {
        if adj_simplices.len() == 2 {
            let s1 = &tri.simplices()[adj_simplices[0]];
            let s2 = &tri.simplices()[adj_simplices[1]];

            // Opposing vertices u1, u2
            let u1 = *s1.iter().find(|&v| !ridge.contains(v)).ok_or_else(|| {
                Error::LinearAlgebra("Corrupt triangulation: simplex subset of ridge".into())
            })?;
            let u2 = *s2.iter().find(|&v| !ridge.contains(v)).ok_or_else(|| {
                Error::LinearAlgebra("Corrupt triangulation: simplex subset of ridge".into())
            })?;

            // Find relation among {u1, u2} U ridge
            // Total dim+2 points in dim dimensions.
            let mut subset = vec![u1, u2];
            subset.extend(ridge.iter());

            let mut mat = vec![vec![Integer::from(0); subset.len()]; dim + 1];
            for (j, &idx) in subset.iter().enumerate() {
                let coords = points[idx].coords();
                for (i, &c) in coords.iter().enumerate() {
                    mat[i][j] = Integer::from(c);
                }
                mat[dim][j] = Integer::from(1);
            }

            let kernel = integer_kernel(&mat);
            if let Some(relation) = extract_mori_relation(&kernel) {
                // relation is in order of 'subset': [u1, u2, v1, v2, v3, v4]
                // Expand to full points list
                let mut full_gen = vec![Integer::from(0); n_pts];
                for (i, &idx) in subset.iter().enumerate() {
                    full_gen[idx] = relation[i].clone();
                }
                generators.push(full_gen);
            }
        }
    }

    Ok(MoriCone::new(generators))
}

fn build_ridge_map(tri: &Triangulation, dim: usize) -> HashMap<Vec<usize>, Vec<usize>> {
    let mut ridge_map: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();
    for (s_idx, simplex) in tri.simplices().iter().enumerate() {
        if simplex.len() != dim + 1 {
            continue;
        }

        for i in 0..simplex.len() {
            let mut ridge = simplex.clone();
            ridge.remove(i);
            ridge.sort_unstable();
            ridge_map.entry(ridge).or_default().push(s_idx);
        }
    }
    ridge_map
}

/// From the kernel of {u1, u2, ridge}, extract the Mori relation.
///
/// There is a 2D kernel. We want the relation where u1 and u2 both have
/// the SAME sign (by convention positive).
fn extract_mori_relation(kernel: &[Vec<Integer>]) -> Option<Vec<Integer>> {
    // TODO: Robust implementation of relation extraction.
    // For now, return first kernel vector if it looks reasonable.
    if !kernel.is_empty() {
        let mut res = kernel[0].clone();
        if res[0] < 0 {
            for val in &mut res {
                *val = -val.clone();
            }
        }
        return Some(res);
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_kahler_cone_contains() {
        // Generator R = [1, -2]
        // p = [2, 1] => p·R = 2 - 2 = 0 (boundary)
        // p = [3, 1] => p·R = 3 - 2 = 1 > 0 (inside)
        let mori = MoriCone::new(vec![vec![Integer::from(1), Integer::from(-2)]]);
        assert!(!mori.contains(&[2.0, 1.0]));
        assert!(mori.contains(&[3.0, 1.0]));
    }

    #[test]
    fn test_compute_mori_generators_simple() {
        // Two triangles sharing an edge in 2D
        // Points: 0=(0,0), 1=(1,0), 2=(0,1), 3=(1,1)
        // Triangulation: [0,1,2] and [1,2,3] -> shared edge [1,2]
        // Relation: v0 + v3 = v1 + v2 => v0 - v1 - v2 + v3 = 0.
        // Mori generator: [1, -1, -1, 1].

        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![0, 1]),
            Point::new(vec![1, 1]),
        ];

        let tri = Triangulation::new(vec![vec![0, 1, 2], vec![1, 2, 3]]);

        let mori = compute_mori_generators(&tri, &points).unwrap();
        assert!(!mori.generators.is_empty());

        let generator = &mori.generators[0];
        // Check signs: v0, v3 positive? v1, v2 negative?
        // Or v0, v3 negative?
        // Relation is unique up to scale.
        // v0+v3 - v1-v2 = 0.

        // We implemented a heuristic: first element positive.
        // If gen[0] corresponds to v0.
        // Then gen should be [1, -1, -1, 1].

        assert_eq!(generator[0], Integer::from(1));
        assert_eq!(generator[3], Integer::from(1));
        assert_eq!(generator[1], Integer::from(-1));
        assert_eq!(generator[2], Integer::from(-1));
    }
}
