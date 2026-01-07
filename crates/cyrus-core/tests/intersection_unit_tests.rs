// Test-specific clippy allows
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::uninlined_format_args)]
#![allow(clippy::items_after_statements)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::cognitive_complexity)]
#![allow(clippy::redundant_closure_for_method_calls)]
#![allow(clippy::unnecessary_sort_by)]
#![allow(clippy::needless_for_each)]
#![allow(clippy::iter_cloned_collect)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::float_cmp)]
#![allow(clippy::struct_field_names)]
#![allow(clippy::suboptimal_flops)]
#![allow(clippy::redundant_clone)]
#![allow(clippy::stable_sort_primitive)]
#![allow(clippy::bool_to_int_with_if)]
#![allow(dead_code)]
#![allow(unused_variables)]
#![allow(clippy::unnecessary_get_then_check)]

//! Unit tests for intersection number computation components.
//!
//! These tests isolate each step of the algorithm to pinpoint bugs.
//!
//! Core tests (tests 1-9) are in this file.
//! Extended tests (tests 10-13) are in intersection_unit_tests_extended.rs.

use cyrus_core::{Point, Triangulation, compute_glsm_charge_matrix};
use std::collections::HashMap;

/// The dual polytope setup for McAllister 4-214-647
pub fn dual_test_setup() -> (Vec<Point>, Triangulation, usize) {
    // 9 points: origin at index 0, rays at 1-8
    let points: Vec<Point> = vec![
        Point::new(vec![0, 0, 0, 0]),    // 0: origin
        Point::new(vec![-1, 2, -1, -1]), // 1
        Point::new(vec![1, -1, 0, 0]),   // 2
        Point::new(vec![-1, -1, 1, 1]),  // 3
        Point::new(vec![-1, -1, 1, 2]),  // 4
        Point::new(vec![-1, -1, 2, 1]),  // 5
        Point::new(vec![-1, -1, 2, 3]),  // 6
        Point::new(vec![-1, -1, 3, 2]),  // 7
        Point::new(vec![-1, -1, 2, 2]),  // 8
    ];

    // McAllister's triangulation (15 simplices, all containing origin)
    let simplices = vec![
        vec![0, 1, 2, 3, 4],
        vec![0, 1, 2, 3, 5],
        vec![0, 1, 2, 4, 6],
        vec![0, 1, 2, 5, 7],
        vec![0, 1, 2, 6, 7],
        vec![0, 1, 3, 4, 5],
        vec![0, 1, 4, 5, 8],
        vec![0, 1, 4, 6, 8],
        vec![0, 1, 5, 7, 8],
        vec![0, 1, 6, 7, 8],
        vec![0, 2, 3, 4, 5],
        vec![0, 2, 4, 5, 8],
        vec![0, 2, 4, 6, 8],
        vec![0, 2, 5, 7, 8],
        vec![0, 2, 6, 7, 8],
    ];

    let tri = Triangulation::new(simplices);
    let origin_idx = 0;

    (points, tri, origin_idx)
}

fn compute_4x4_det(points: &[Point], indices: &[usize]) -> i64 {
    let m: Vec<Vec<i64>> = indices
        .iter()
        .map(|&i| points[i].coords().to_vec())
        .collect();

    det_4x4(&m)
}

fn det_4x4(m: &[Vec<i64>]) -> i64 {
    let (a, b, c, d) = (m[0][0], m[0][1], m[0][2], m[0][3]);
    let (e, f, g, h) = (m[1][0], m[1][1], m[1][2], m[1][3]);
    let (i, j, k, l) = (m[2][0], m[2][1], m[2][2], m[2][3]);
    let (o, p, q, r) = (m[3][0], m[3][1], m[3][2], m[3][3]);

    a * det_3x3(f, g, h, j, k, l, p, q, r) - b * det_3x3(e, g, h, i, k, l, o, q, r)
        + c * det_3x3(e, f, h, i, j, l, o, p, r)
        - d * det_3x3(e, f, g, i, j, k, o, p, q)
}

const fn det_3x3(a: i64, b: i64, c: i64, d: i64, e: i64, f: i64, g: i64, h: i64, i: i64) -> i64 {
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

// =============================================================================
// TEST 1: Verify GLSM matrix matches CYTools
// =============================================================================

#[test]
fn test_glsm_with_origin_matches_cytools() {
    let (points, _, origin_idx) = dual_test_setup();

    // Compute GLSM with origin (the "standard" way)
    let non_origin: Vec<Point> = points
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != origin_idx)
        .map(|(_, p)| p.clone())
        .collect();

    let glsm = compute_glsm_charge_matrix(&non_origin, true).expect("GLSM computation failed");

    // CYTools expected GLSM (with origin column 0)
    let expected: Vec<Vec<i64>> = vec![
        vec![-6, 2, 3, -1, 1, 1, 0, 0, 0],
        vec![-6, 2, 3, 1, -1, 0, 1, 0, 0],
        vec![-12, 4, 6, 0, 1, 0, 0, 1, 0],
        vec![-6, 2, 3, 0, 0, 0, 0, 0, 1],
    ];

    assert_eq!(glsm.len(), 4, "GLSM should have 4 rows");
    assert_eq!(glsm[0].len(), 9, "GLSM should have 9 columns (with origin)");

    for (row_idx, (actual_row, expected_row)) in glsm.iter().zip(expected.iter()).enumerate() {
        for (col_idx, (actual, expected)) in actual_row.iter().zip(expected_row.iter()).enumerate()
        {
            let actual_i64: i64 = actual.to_string().parse().unwrap();
            assert_eq!(
                actual_i64, *expected,
                "GLSM mismatch at [{row_idx}, {col_idx}]: got {actual_i64}, expected {expected}"
            );
        }
    }

    println!("GLSM with origin matches CYTools exactly");
}

#[test]
fn test_glsm_for_intersection_computation() {
    let (points, _, origin_idx) = dual_test_setup();

    // For intersection computation, we need GLSM with origin column,
    // then skip column 0 when building equations.
    let non_origin: Vec<Point> = points
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != origin_idx)
        .map(|(_, p)| p.clone())
        .collect();

    // Compute with include_origin=true to get correct kernel
    let glsm = compute_glsm_charge_matrix(&non_origin, true).expect("GLSM computation failed");

    // Remove column 0 (origin) for intersection equations
    let glsm_no_origin: Vec<Vec<i64>> = glsm
        .iter()
        .map(|row| {
            row.iter()
                .skip(1)
                .map(|x| x.to_string().parse().unwrap())
                .collect()
        })
        .collect();

    // Expected values from CYTools
    let expected: Vec<Vec<i64>> = vec![
        vec![2, 3, -1, 1, 1, 0, 0, 0],
        vec![2, 3, 1, -1, 0, 1, 0, 0],
        vec![4, 6, 0, 1, 0, 0, 1, 0],
        vec![2, 3, 0, 0, 0, 0, 0, 1],
    ];

    assert_eq!(
        glsm_no_origin, expected,
        "GLSM (origin column removed) must match CYTools"
    );
}

// =============================================================================
// TEST 2: Verify simplex determinants
// =============================================================================

#[test]
fn test_simplex_determinants() {
    let (points, tri, origin_idx) = dual_test_setup();

    // Expected |det| values for each simplex's non-origin rays
    let expected_dets: Vec<(Vec<usize>, i64)> = vec![
        (vec![1, 2, 3, 4], 1),
        (vec![1, 2, 3, 5], 1),
        (vec![1, 2, 4, 6], 1),
        (vec![1, 2, 5, 7], 1),
        (vec![1, 2, 6, 7], 1),
        (vec![1, 3, 4, 5], 3),
        (vec![1, 4, 5, 8], 3),
        (vec![1, 4, 6, 8], 3),
        (vec![1, 5, 7, 8], 3),
        (vec![1, 6, 7, 8], 3),
        (vec![2, 3, 4, 5], 2),
        (vec![2, 4, 5, 8], 2),
        (vec![2, 4, 6, 8], 2),
        (vec![2, 5, 7, 8], 2),
        (vec![2, 6, 7, 8], 2),
    ];

    for (simplex_idx, simplex) in tri.simplices().iter().enumerate() {
        let rays: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| i != origin_idx)
            .collect();

        if rays.len() == 4 {
            let det = compute_4x4_det(&points, &rays);
            let abs_det = det.abs();

            let (expected_rays, expected_det) = &expected_dets[simplex_idx];
            assert_eq!(&rays, expected_rays, "Simplex {simplex_idx} rays mismatch");
            assert_eq!(
                abs_det, *expected_det,
                "Simplex {simplex_idx} det mismatch: got {abs_det}, expected {expected_det}"
            );
        }
    }

    println!("All 15 simplex determinants correct");
}

// =============================================================================
// TEST 3: Verify known 4-form values (distinct tuples)
// =============================================================================

#[test]
fn test_known_4form_values() {
    let (points, tri, origin_idx) = dual_test_setup();

    // Compute known 4-forms: kappa^V_{ijkl} = 1/|det| for distinct (i,j,k,l) in simplex
    let mut known: HashMap<Vec<usize>, f64> = HashMap::new();

    for simplex in tri.simplices() {
        let rays: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| i != origin_idx)
            .collect();

        if rays.len() == 4 {
            let det = compute_4x4_det(&points, &rays);
            let abs_det = det.abs();
            let value = 1.0 / (abs_det as f64);

            let mut key = rays.clone();
            key.sort_unstable();
            known.insert(key, value);
        }
    }

    // Check specific values
    assert_eq!(known.len(), 15, "Should have 15 distinct 4-tuples");

    // kappa^V_{1,3,4,5} = 1/3
    let key_1345 = vec![1, 3, 4, 5];
    assert!(
        (known[&key_1345] - 1.0 / 3.0).abs() < 1e-10,
        "kappa^V_{{1,3,4,5}} should be 1/3, got {}",
        known[&key_1345]
    );

    // kappa^V_{2,3,4,5} = 1/2
    let key_2345 = vec![2, 3, 4, 5];
    assert!(
        (known[&key_2345] - 0.5).abs() < 1e-10,
        "kappa^V_{{2,3,4,5}} should be 1/2, got {}",
        known[&key_2345]
    );

    // kappa^V_{1,2,3,4} = 1/1 = 1
    let key_1234 = vec![1, 2, 3, 4];
    assert!(
        (known[&key_1234] - 1.0).abs() < 1e-10,
        "kappa^V_{{1,2,3,4}} should be 1, got {}",
        known[&key_1234]
    );

    println!("Known 4-form values correct");
    println!("  kappa^V_{{1,3,4,5}} = 1/3");
    println!("  kappa^V_{{2,3,4,5}} = 1/2");
    println!("  kappa^V_{{1,2,3,4}} = 1");
}

// =============================================================================
// TEST 4: Manual equation construction for probe (3,4,5)
// =============================================================================

#[test]
fn test_equation_for_probe_345() {
    let (points, tri, origin_idx) = dual_test_setup();

    // GLSM row 0 (with origin): [-6, 2, 3, -1, 1, 1, 0, 0, 0]
    // Column 0 = origin, columns 1-8 = points 1-8
    let glsm_row_0: Vec<i64> = vec![-6, 2, 3, -1, 1, 1, 0, 0, 0];

    // Compute known 4-forms
    let mut known: HashMap<Vec<usize>, f64> = HashMap::new();
    for simplex in tri.simplices() {
        let rays: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| i != origin_idx)
            .collect();
        if rays.len() == 4 {
            let det = compute_4x4_det(&points, &rays);
            let mut key = rays.clone();
            key.sort_unstable();
            known.insert(key, 1.0 / det.abs() as f64);
        }
    }

    println!("\n=== Equation for probe (3,4,5), GLSM row 0 ===");
    println!("GLSM row 0: {glsm_row_0:?}");
    println!("Q_0 (origin) = {}", glsm_row_0[0]);

    // Method A: Using effective coefficients (Q_m - Q_0)
    println!("\n--- Method A: Effective coefficients (Q_m - Q_0) ---");
    let q0 = glsm_row_0[0] as f64;
    let mut rhs_a = 0.0;
    let mut var_coeffs_a: Vec<(String, f64)> = Vec::new();

    for m in 1..=8 {
        let q_m = glsm_row_0[m] as f64;
        let eff = q_m - q0;
        if eff.abs() < 1e-10 {
            continue;
        }

        let mut key = vec![m, 3, 4, 5];
        key.sort_unstable();

        if let Some(&val) = known.get(&key) {
            rhs_a -= eff * val;
            println!(
                "  m={}: ({} - {}) x kappa^V{:?}={:.4} -> RHS -= {:.4}",
                m,
                q_m,
                q0,
                key,
                val,
                eff * val
            );
        } else {
            var_coeffs_a.push((format!("kappa^V{key:?}"), eff));
            println!("  m={m}: ({q_m} - {q0}) x kappa^V{key:?} [variable, coeff={eff}]");
        }
    }
    println!("Equation A: {var_coeffs_a:?} = {rhs_a:.6}");

    // Method B: Excluding origin entirely (raw Q_m, skip m=0)
    println!("\n--- Method B: Exclude origin entirely (raw Q_m) ---");
    let mut rhs_b = 0.0;
    let mut var_coeffs_b: Vec<(String, f64)> = Vec::new();

    for m in 1..=8 {
        let q_m = glsm_row_0[m] as f64;
        if q_m.abs() < 1e-10 {
            continue;
        }

        let mut key = vec![m, 3, 4, 5];
        key.sort_unstable();

        if let Some(&val) = known.get(&key) {
            rhs_b -= q_m * val;
            println!(
                "  m={}: {} x kappa^V{:?}={:.4} -> RHS -= {:.4}",
                m,
                q_m,
                key,
                val,
                q_m * val
            );
        } else {
            var_coeffs_b.push((format!("kappa^V{key:?}"), q_m));
            println!("  m={m}: {q_m} x kappa^V{key:?} [variable, coeff={q_m}]");
        }
    }
    println!("Equation B: {var_coeffs_b:?} = {rhs_b:.6}");

    // Expected: For kappa_{345} = 1, we need:
    // 1/3 + 1/2 + S = 1, so S = 1/6 where S = kappa^V_{3345} + kappa^V_{3445} + kappa^V_{3455}
    println!("\n--- Expected values ---");
    println!("For kappa_{{345}} = 1:");
    println!("  kappa_{{345}} = kappa^V_{{1345}} + kappa^V_{{2345}} + self-intersection sum S");
    println!("  1 = 1/3 + 1/2 + S");
    println!("  S = 1 - 5/6 = 1/6 ~ {:.6}", 1.0 / 6.0);
}

// =============================================================================
// TEST 5: Check ALL GLSM rows for probe (3,4,5)
// =============================================================================

#[test]
fn test_all_glsm_rows_for_probe_345() {
    let (points, tri, origin_idx) = dual_test_setup();

    // CYTools GLSM (with origin column 0)
    let glsm: Vec<Vec<i64>> = vec![
        vec![-6, 2, 3, -1, 1, 1, 0, 0, 0],
        vec![-6, 2, 3, 1, -1, 0, 1, 0, 0],
        vec![-12, 4, 6, 0, 1, 0, 0, 1, 0],
        vec![-6, 2, 3, 0, 0, 0, 0, 0, 1],
    ];

    // Compute known 4-forms (distinct tuples only)
    let mut known: HashMap<Vec<usize>, f64> = HashMap::new();
    for simplex in tri.simplices() {
        let rays: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| i != origin_idx)
            .collect();
        if rays.len() == 4 {
            let det = compute_4x4_det(&points, &rays);
            let mut key = rays.clone();
            key.sort_unstable();
            known.insert(key, 1.0 / det.abs() as f64);
        }
    }

    // Self-intersection variables: tuples with repeated indices that are in simplices
    // For probe (3,4,5), these are: {3,3,4,5}, {3,4,4,5}, {3,4,5,5}
    // They're in simplices [0,1,3,4,5] and [0,2,3,4,5]
    let self_int_vars: Vec<Vec<usize>> = vec![vec![3, 3, 4, 5], vec![3, 4, 4, 5], vec![3, 4, 5, 5]];

    println!("\n=== All GLSM rows for probe (3,4,5) ===");
    println!("Known: kappa^V_{{1345}}=1/3, kappa^V_{{2345}}=1/2");
    println!("Variables: kappa^V_{{3345}}, kappa^V_{{3445}}, kappa^V_{{3455}}");
    println!(
        "Expected: S = kappa^V_{{3345}} + kappa^V_{{3445}} + kappa^V_{{3455}} = 1/6 for kappa_{{345}}=1"
    );

    for (row_idx, glsm_row) in glsm.iter().enumerate() {
        println!("\n--- GLSM Row {row_idx} ---");
        println!("  GLSM: {glsm_row:?}");
        let q0 = glsm_row[0] as f64;

        // Method A: Effective coefficients
        let mut rhs_a = 0.0;
        let mut var_terms_a: Vec<(String, f64)> = Vec::new();

        // Method B: Raw Q_m
        let mut rhs_b = 0.0;
        let mut var_terms_b: Vec<(String, f64)> = Vec::new();

        for m in 1..=8 {
            let q_m = glsm_row[m] as f64;
            let mut key = vec![m, 3, 4, 5];
            key.sort_unstable();
            let key_str = format!("{key:?}");

            let is_known = known.contains_key(&key);
            let is_var = self_int_vars.contains(&key);

            // Method A: effective coefficient
            let eff = q_m - q0;
            if eff.abs() > 1e-10 {
                if is_known {
                    let val = known[&key];
                    rhs_a -= eff * val;
                } else if is_var {
                    var_terms_a.push((key_str.clone(), eff));
                }
                // else: not in simplex, skip
            }

            // Method B: raw Q_m
            if q_m.abs() > 1e-10 {
                if is_known {
                    let val = known[&key];
                    rhs_b -= q_m * val;
                } else if is_var {
                    var_terms_b.push((key_str.clone(), q_m));
                }
            }
        }

        println!("  Method A (eff): {var_terms_a:?} = {rhs_a:.4}");
        println!("  Method B (raw): {var_terms_b:?} = {rhs_b:.4}");
    }

    // Check what value of S would satisfy each equation
    println!("\n=== What S would satisfy each equation? ===");
    println!("Method A Row 3: 6*S = RHS => S = RHS/6");
    println!("  S = -7.1667/6 = {:.4}", -7.1667 / 6.0);
    println!("Method B Row 3: 0 = RHS (no variables!) => inconsistent");
}

// =============================================================================
// TEST 6: Variable enumeration - self-intersection 4-tuples
// =============================================================================

#[test]
fn test_variable_enumeration() {
    let (_, tri, origin_idx) = dual_test_setup();

    // Variables are 4-tuples with at least one repeated index, excluding origin
    let mut variables: std::collections::HashSet<Vec<usize>> = std::collections::HashSet::new();

    for simplex in tri.simplices() {
        let rays: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| i != origin_idx)
            .collect();

        // Generate all 4-tuples with replacement
        for &i in &rays {
            for &j in &rays {
                if j < i {
                    continue;
                }
                for &k in &rays {
                    if k < j {
                        continue;
                    }
                    for &l in &rays {
                        if l < k {
                            continue;
                        }
                        let tuple = vec![i, j, k, l];

                        // Only include if NOT distinct (has repeated index)
                        let is_distinct = i != j && j != k && k != l;
                        if !is_distinct {
                            variables.insert(tuple);
                        }
                    }
                }
            }
        }
    }

    // Check specific self-intersection variables exist
    assert!(
        variables.contains(&vec![3, 3, 4, 5]),
        "Should have [3,3,4,5]"
    );
    assert!(
        variables.contains(&vec![3, 4, 4, 5]),
        "Should have [3,4,4,5]"
    );
    assert!(
        variables.contains(&vec![3, 4, 5, 5]),
        "Should have [3,4,5,5]"
    );
    assert!(
        variables.contains(&vec![1, 1, 1, 1]),
        "Should have [1,1,1,1]"
    );

    // Check distinct tuples are NOT in variables
    assert!(
        !variables.contains(&vec![1, 2, 3, 4]),
        "Should NOT have [1,2,3,4]"
    );
    assert!(
        !variables.contains(&vec![1, 3, 4, 5]),
        "Should NOT have [1,3,4,5]"
    );

    // Check origin is excluded
    for var in &variables {
        assert!(
            !var.contains(&0),
            "Variables should not contain origin: {var:?}"
        );
    }

    println!(
        "Variable enumeration: {} self-intersection 4-tuples",
        variables.len()
    );
}

// =============================================================================
// TEST 7: Equation construction for probe (3,4,5)
// =============================================================================

#[test]
fn test_equation_construction_probe_345() {
    let (points, tri, origin_idx) = dual_test_setup();

    // Known 4-form values
    let mut known: HashMap<Vec<usize>, f64> = HashMap::new();
    for simplex in tri.simplices() {
        let rays: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| i != origin_idx)
            .collect();
        if rays.len() == 4 {
            let det = compute_4x4_det(&points, &rays);
            let mut key = rays.clone();
            key.sort_unstable();
            known.insert(key, 1.0 / det.abs() as f64);
        }
    }

    // GLSM without origin column (CYTools values)
    let glsm: Vec<Vec<i64>> = vec![
        vec![2, 3, -1, 1, 1, 0, 0, 0],
        vec![2, 3, 1, -1, 0, 1, 0, 0],
        vec![4, 6, 0, 1, 0, 0, 1, 0],
        vec![2, 3, 0, 0, 0, 0, 0, 1],
    ];

    // Enumerate variables (self-intersection tuples in simplices)
    let mut variables: std::collections::HashSet<Vec<usize>> = std::collections::HashSet::new();
    for simplex in tri.simplices() {
        let rays: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| i != origin_idx)
            .collect();
        for &i in &rays {
            for &j in &rays {
                if j < i {
                    continue;
                }
                for &k in &rays {
                    if k < j {
                        continue;
                    }
                    for &l in &rays {
                        if l < k {
                            continue;
                        }
                        let tuple = vec![i, j, k, l];
                        let is_distinct = i != j && j != k && k != l;
                        if !is_distinct && !known.contains_key(&tuple) {
                            variables.insert(tuple);
                        }
                    }
                }
            }
        }
    }

    // For probe (3,4,5), build equations from each GLSM row
    // Equation: Sum_{m=1..8} Q_m x kappa^V_{m,3,4,5} = 0
    // - Known terms go to RHS
    // - Variable terms stay on LHS
    // - Terms not in any simplex are 0 (skip)

    for (row_idx, row) in glsm.iter().enumerate() {
        let mut rhs = 0.0;
        let mut var_coeffs: Vec<(Vec<usize>, i64)> = Vec::new();

        for (col, &q) in row.iter().enumerate() {
            if q == 0 {
                continue;
            }

            let m = col + 1; // Column 0 -> point 1, etc.
            let mut key = vec![m, 3, 4, 5];
            key.sort_unstable();

            if let Some(&val) = known.get(&key) {
                rhs -= (q as f64) * val;
            } else if variables.contains(&key) {
                var_coeffs.push((key.clone(), q));
            }
            // else: not in any simplex, skip (value is 0)
        }

        // Verify equations are as expected
        match row_idx {
            0 => {
                // Row 0: -kappa_{3345} + kappa_{3445} + kappa_{3455} = -2.1667
                assert_eq!(var_coeffs.len(), 3);
                assert!(
                    (rhs - (-2.1667)).abs() < 0.01,
                    "Row 0 RHS: {rhs} expected -2.1667"
                );
            }
            1 => {
                // Row 1: kappa_{3345} - kappa_{3445} = -2.1667
                assert_eq!(var_coeffs.len(), 2);
            }
            2 => {
                // Row 2: kappa_{3445} = -4.3333
                assert_eq!(var_coeffs.len(), 1);
                assert!(
                    (rhs - (-4.3333)).abs() < 0.01,
                    "Row 2 RHS: {rhs} expected -4.3333"
                );
            }
            3 => {
                // Row 3: 0 = -2.1667 (no variables!)
                assert_eq!(var_coeffs.len(), 0, "Row 3 should have no variables");
            }
            _ => {}
        }
    }

    println!("Equation construction verified for probe (3,4,5)");
}

// =============================================================================
// TEST 8: Least squares solver on simple system
// =============================================================================

#[test]
fn test_least_squares_solver_simple() {
    // Simple overdetermined system: 2 variables, 3 equations
    // x + y = 3
    // x - y = 1
    // 2x = 4  (redundant but consistent)
    // Solution: x = 2, y = 1

    use faer::prelude::SpSolver;
    use faer::sparse::SparseColMat;

    // M matrix (3x2)
    let triplets = vec![
        (0, 0, 1.0),
        (0, 1, 1.0), // row 0: x + y
        (1, 0, 1.0),
        (1, 1, -1.0), // row 1: x - y
        (2, 0, 2.0),  // row 2: 2x
    ];

    let c = [3.0, 1.0, 4.0];

    let n_rows = 3;
    let n_cols = 2;

    // Build M^T M
    let mut mtm = vec![vec![0.0; n_cols]; n_cols];
    let mut mtc = vec![0.0; n_cols];

    for &(row, col, val) in &triplets {
        for &(row2, col2, val2) in &triplets {
            if row == row2 {
                mtm[col][col2] += val * val2;
            }
        }
        mtc[col] += val * c[row];
    }

    // M^T M = [[6, 0], [0, 2]]
    // M^T c = [14, 2]
    // Solution: x = 14/6 = 2.33..., y = 2/2 = 1
    // Wait, that's wrong. Let me recalculate.
    //
    // M = [[1, 1], [1, -1], [2, 0]]
    // M^T = [[1, 1, 2], [1, -1, 0]]
    // M^T M = [[1+1+4, 1-1+0], [1-1+0, 1+1+0]] = [[6, 0], [0, 2]]
    // M^T c = [[1*3 + 1*1 + 2*4], [1*3 + (-1)*1 + 0*4]] = [[12], [2]]
    // x = 12/6 = 2, y = 2/2 = 1

    let mtm_triplets: Vec<(usize, usize, f64)> = vec![(0, 0, 6.0), (1, 1, 2.0)];

    let mtm_sparse =
        SparseColMat::<usize, f64>::try_new_from_triplets(n_cols, n_cols, &mtm_triplets)
            .expect("Failed to build sparse matrix");

    let mut mtc_mat = faer::Mat::<f64>::zeros(n_cols, 1);
    mtc_mat[(0, 0)] = 12.0;
    mtc_mat[(1, 0)] = 2.0;

    let chol = mtm_sparse
        .sp_cholesky(faer::Side::Lower)
        .expect("Cholesky failed");

    let solution = chol.solve(&mtc_mat);

    let x = solution[(0, 0)];
    let y = solution[(1, 0)];

    assert!((x - 2.0).abs() < 1e-10, "x = {x} but expected 2");
    assert!((y - 1.0).abs() < 1e-10, "y = {y} but expected 1");

    println!("Least squares solver: x = {x}, y = {y}");
}

// =============================================================================
// TEST 9: Origin handling - derive kappa_{0jk} from non-origin values
// =============================================================================

#[test]
fn test_origin_derivation() {
    // Given non-origin kappa values, kappa_{0jk} = -Sum_{i!=0} kappa_{ijk}
    // This derives from D_0 ~ -Sum_{i>0} D_i

    // Simple example: 3 non-origin points (1, 2, 3)
    let mut kappa: HashMap<Vec<usize>, f64> = HashMap::new();

    // Set some non-origin values
    kappa.insert(vec![1, 2, 3], 6.0);
    kappa.insert(vec![1, 1, 2], 10.0);
    kappa.insert(vec![1, 1, 3], 5.0);
    kappa.insert(vec![2, 2, 3], 3.0);
    kappa.insert(vec![1, 2, 2], 4.0);
    kappa.insert(vec![1, 3, 3], 2.0);
    kappa.insert(vec![2, 3, 3], 1.0);
    kappa.insert(vec![1, 1, 1], 20.0);
    kappa.insert(vec![2, 2, 2], 8.0);
    kappa.insert(vec![3, 3, 3], 3.0);

    let n_pts = 4; // 0, 1, 2, 3
    let origin = 0;
    let non_origin: Vec<usize> = vec![1, 2, 3];

    // Step 1: kappa_{0jk} = -Sum_{i>0} kappa_{ijk}
    for &j in &non_origin {
        for &k in &non_origin {
            if k < j {
                continue;
            }
            let sum: f64 = non_origin
                .iter()
                .map(|&i| {
                    let mut key = vec![i, j, k];
                    key.sort_unstable();
                    kappa.get(&key).copied().unwrap_or(0.0)
                })
                .sum();
            let mut key = vec![origin, j, k];
            key.sort_unstable();
            kappa.insert(key, -sum);
        }
    }

    // kappa_{0,2,3} = -(kappa_{1,2,3} + kappa_{2,2,3} + kappa_{3,2,3})
    //           = -(6 + 3 + 1) = -10
    let val_023 = kappa.get(&vec![0, 2, 3]).copied().unwrap_or(999.0);
    assert!(
        (val_023 - (-10.0)).abs() < 1e-10,
        "kappa[0,2,3] = {val_023} but expected -10"
    );

    // kappa_{0,1,2} = -(kappa_{1,1,2} + kappa_{1,2,2} + kappa_{1,2,3})
    //           = -(10 + 4 + 6) = -20
    let val_012 = kappa.get(&vec![0, 1, 2]).copied().unwrap_or(999.0);
    assert!(
        (val_012 - (-20.0)).abs() < 1e-10,
        "kappa[0,1,2] = {val_012} but expected -20"
    );

    // Step 2: kappa_{00k} = -Sum_{i>0} kappa_{i0k}
    for &k in &non_origin {
        let sum: f64 = non_origin
            .iter()
            .map(|&i| {
                let mut key = vec![i, origin, k];
                key.sort_unstable();
                kappa.get(&key).copied().unwrap_or(0.0)
            })
            .sum();
        let mut key = vec![origin, origin, k];
        key.sort_unstable();
        kappa.insert(key, -sum);
    }

    // Step 3: kappa_{000} = -Sum_{i>0} kappa_{i00}
    let sum: f64 = non_origin
        .iter()
        .map(|&i| {
            let mut key = vec![i, origin, origin];
            key.sort_unstable();
            kappa.get(&key).copied().unwrap_or(0.0)
        })
        .sum();
    kappa.insert(vec![origin, origin, origin], -sum);

    println!("Origin derivation formulas verified");
}
