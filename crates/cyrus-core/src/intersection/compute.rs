//! Intersection number computation algorithm.
//!
//! Computes κ_{ijk} = (D_i · D_j · D_k)_X for a CY hypersurface X.
//!
//! Algorithm (from clean room docs):
//! 1. Compute known 4-form values: distinct tuples have κ^V_{ijkl} = 1/|det(pts)|
//! 2. Enumerate unknown 4-form variables (self-intersection, at least one repeated index)
//! 3. Build GLSM constraints using 3-face probes: Σ_m (Q_m - Q_0) κ^V_{m,face} = 0
//! 4. Solve for unknown 4-forms
//! 5. Sum ALL 4-forms to get 3-forms: κ_{ijk} = Σ_l κ^V_{ijkl}
//! 6. Derive origin-indexed values via linear equivalence
//!
//! Reference: [[project_docs/clean_room/INTERSECTION_NUMBERS.md]]
//! Reference: [[project_docs/clean_room/INTERSECTION_TRANSFORMATION.md]]

use super::Intersection;
use super::four_form::{
    collect_3face_probes, compute_known_4form_values, enumerate_4form_variables,
};
use super::helpers::solve_normal_equations;
use super::linear_system::{
    build_4form_glsm_system_no_origin, build_4form_system_with_linear_relations,
};
use super::three_form::{compute_3form_from_4form, derive_origin_values, kappa_to_intersection};
use crate::Point;
use crate::error::Result;
use crate::triangulation::Triangulation;
use malachite::Integer;

/// Compute the intersection numbers for a Calabi-Yau hypersurface.
pub fn compute_intersection_numbers(
    tri: &Triangulation,
    points: &[Point],
    glsm: &[Vec<Integer>],
) -> Result<Intersection> {
    compute_intersection_numbers_with_offset(tri, points, glsm, None)
}

/// Compute intersection numbers with explicit linear relations matrix.
///
/// This is the preferred method when you have the CYTools-style reduced linear_relations
/// matrix (with -I on the left and basis coefficients on the right).
///
/// # Arguments
/// * `tri` - Triangulation
/// * `points` - Polytope points (including origin at index 0)
/// * `linear_relations` - Reduced h11 x (n_points-1) matrix from CYTools
///
/// The linear_relations matrix should be in CYTools form:
/// - Columns correspond to non-origin points (indices 1 to n_points-1)
/// - Each row expresses one non-basis divisor in terms of basis divisors
/// - Has -1 in the pivot columns (non-basis) and coefficients in basis columns
pub fn compute_intersection_numbers_with_linear_relations(
    tri: &Triangulation,
    points: &[Point],
    linear_relations: &[Vec<Integer>],
) -> Result<Intersection> {
    let n_pts = points.len();
    let origin_idx = points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0));

    eprintln!("[DEBUG] Origin at index: {origin_idx:?}");

    // Step 1: Compute known 4-form values
    let known_4form = compute_known_4form_values(tri, points, origin_idx);
    eprintln!("[DEBUG] Known 4-form values: {}", known_4form.len());

    // Step 2: Enumerate unknown variables
    let (var_list, var_map) = enumerate_4form_variables(tri, origin_idx, &known_4form);
    eprintln!("[DEBUG] 4-form variables: {}", var_list.len());

    // Step 3: Collect probes
    let probe_3faces = collect_3face_probes(tri, origin_idx);
    eprintln!("[DEBUG] 3-face probes: {}", probe_3faces.len());

    // Step 4: Build and solve using the reduced linear_relations
    let solved_4form = if var_list.is_empty() {
        known_4form
    } else {
        let (triplets, c_vec) = build_4form_system_with_linear_relations(
            &probe_3faces,
            &known_4form,
            &var_map,
            n_pts,
            origin_idx,
            linear_relations,
        );

        let n_equations = c_vec.len();
        let n_vars = var_list.len();
        eprintln!("[DEBUG] Equations: {n_equations}, Variables: {n_vars}");

        if n_equations == 0 {
            known_4form
        } else {
            let sol = solve_normal_equations(&triplets, &c_vec, n_equations, n_vars)?;

            let mut all_4form = known_4form;
            for (i, key) in var_list.iter().enumerate() {
                all_4form.insert(key.clone(), sol[i]);
            }
            all_4form
        }
    };

    // Debug some solved 4-forms
    for key in [vec![3, 3, 4, 5], vec![3, 4, 4, 5], vec![3, 4, 5, 5]] {
        if let Some(&val) = solved_4form.get(&key) {
            eprintln!("[DEBUG] κ^V{key:?} = {val:.6}");
        }
    }

    // Step 5: Sum 4-forms to get 3-forms
    let kappa_3form = compute_3form_from_4form(tri, origin_idx, &solved_4form, n_pts);

    for key in [vec![3, 4, 5], vec![1, 2, 3], vec![1, 1, 1]] {
        if let Some(&val) = kappa_3form.get(&key) {
            eprintln!("[DEBUG] κ^X{key:?} = {val:.6}");
        }
    }

    // Step 6: Derive origin values
    let kappa_full = derive_origin_values(&kappa_3form, n_pts, origin_idx);

    Ok(kappa_to_intersection(n_pts, &kappa_full))
}

/// Compute intersection numbers with explicit GLSM index offset.
///
/// Algorithm: 4-form approach with origin exclusion
/// 1. Compute known 4-form values: distinct tuples (all indices ≠ 0) have κ^V_{ijkl} = 1/|det(pts)|
/// 2. Enumerate unknown 4-form variables (self-intersection, all indices ≠ 0)
/// 3. Build GLSM constraints using 3-face probes (all indices ≠ 0), GLSM without origin column
/// 4. Solve for unknown 4-forms via least squares
/// 5. Sum ALL 4-forms to get 3-forms: κ_{ijk} = Σ_l κ^V_{ijkl}
/// 6. Derive origin-indexed values via linear equivalence D_0 ~ -Σ_{i>0} D_i
pub fn compute_intersection_numbers_with_offset(
    tri: &Triangulation,
    points: &[Point],
    glsm: &[Vec<Integer>],
    glsm_point_offset: Option<i32>,
) -> Result<Intersection> {
    let n_pts = points.len();

    // Find origin index
    let origin_idx = points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0));

    eprintln!("[DEBUG] Origin at index: {origin_idx:?}");

    // Step 1: Compute known 4-form values for DISTINCT non-origin 4-tuples
    let known_4form = compute_known_4form_values(tri, points, origin_idx);
    eprintln!(
        "[DEBUG] Known 4-form values (distinct): {}",
        known_4form.len()
    );

    // Step 2: Enumerate unknown 4-form variables (self-intersection tuples, all indices ≠ 0)
    let (var_list, var_map) = enumerate_4form_variables(tri, origin_idx, &known_4form);
    eprintln!(
        "[DEBUG] 4-form variables (self-intersection): {}",
        var_list.len()
    );

    // Step 3: Collect 3-face probes (all indices ≠ 0)
    let probe_3faces = collect_3face_probes(tri, origin_idx);
    eprintln!("[DEBUG] 3-face probes: {}", probe_3faces.len());

    // Step 4: Build and solve GLSM system (excluding origin entirely)
    let solved_4form = if var_list.is_empty() {
        known_4form
    } else {
        let (triplets, c_vec) = build_4form_glsm_system_no_origin(
            &probe_3faces,
            &known_4form,
            &var_map,
            n_pts,
            origin_idx,
            glsm,
            glsm_point_offset,
        );

        let n_equations = c_vec.len();
        let n_vars = var_list.len();
        eprintln!("[DEBUG] GLSM equations: {n_equations}, Variables: {n_vars}");

        if n_equations == 0 {
            known_4form
        } else {
            let sol = solve_normal_equations(&triplets, &c_vec, n_equations, n_vars)?;

            // Combine known and solved values
            let mut all_4form = known_4form;
            for (i, key) in var_list.iter().enumerate() {
                all_4form.insert(key.clone(), sol[i]);
            }
            all_4form
        }
    };

    // Debug some solved 4-forms
    for key in [vec![3, 3, 4, 5], vec![3, 4, 4, 5], vec![3, 4, 5, 5]] {
        if let Some(&val) = solved_4form.get(&key) {
            eprintln!("[DEBUG] κ^V{key:?} = {val:.6}");
        }
    }

    // Step 5: Sum 4-forms to get 3-forms: κ_{ijk} = Σ_l κ^V_{ijkl} (l ≠ 0)
    let kappa_3form = compute_3form_from_4form(tri, origin_idx, &solved_4form, n_pts);

    // Debug: Check specific 3-form values
    for key in [vec![3, 4, 5], vec![1, 2, 3], vec![1, 1, 1]] {
        if let Some(&val) = kappa_3form.get(&key) {
            eprintln!("[DEBUG] κ^X{key:?} = {val:.6}");
        } else {
            eprintln!("[DEBUG] κ^X{key:?} = not found");
        }
    }

    // Step 6: Derive origin-indexed values via D_0 ~ -Σ_{i>0} D_i
    let kappa_full = derive_origin_values(&kappa_3form, n_pts, origin_idx);

    // Step 7: Convert to Intersection struct
    Ok(kappa_to_intersection(n_pts, &kappa_full))
}
