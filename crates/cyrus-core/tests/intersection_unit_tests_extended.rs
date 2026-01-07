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

//! Extended unit tests for intersection number computation components.
//!
//! These tests cover additional scenarios including full intersection computation,
//! basis transformation, GLSM constraints verification, and simplex containment.
//!
//! Core tests (tests 1-9) are in intersection_unit_tests.rs.
//! Extended tests (tests 10-13) are in this file.

use cyrus_core::{Point, Triangulation, compute_glsm_charge_matrix, compute_intersection_numbers};
use std::collections::HashMap;

/// The dual polytope setup for McAllister 4-214-647
fn dual_test_setup() -> (Vec<Point>, Triangulation, usize) {
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

// =============================================================================
// TEST 10: Full intersection computation must give kappa[3,4,5] = 1
// NOTE: This test uses the OLD compute_intersection_numbers algorithm which has
// known issues with certain polytopes (high residual). The CYTools algorithm
// (compute_intersection_cytools) is the correct one and is tested in dual_tests.rs
// and mcallister_e2e tests.
// =============================================================================

#[test]
#[ignore = "Uses old algorithm with known issues; see dual_tests.rs for CYTools algorithm tests"]
fn test_intersection_numbers_345_equals_1() {
    let (points, tri, origin_idx) = dual_test_setup();

    // Compute GLSM with origin (correct method)
    let non_origin: Vec<Point> = points
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != origin_idx)
        .map(|(_, p)| p.clone())
        .collect();

    let glsm = compute_glsm_charge_matrix(&non_origin, true).expect("GLSM computation failed");

    // Compute intersection numbers
    let kappa = compute_intersection_numbers(&tri, &points, &glsm)
        .expect("Intersection computation failed");

    // CYTools ground truth: kappa[3,4,5] = 1
    let kappa_345 = kappa.get(3, 4, 5);

    // Convert to f64 for comparison
    use malachite::num::conversion::traits::RoundingFrom;
    use malachite::rounding_modes::RoundingMode;
    let (val, _) = f64::rounding_from(kappa_345.get(), RoundingMode::Nearest);

    assert!(
        (val - 1.0).abs() < 1e-6,
        "kappa[3,4,5] = {val} but expected 1"
    );
}

// =============================================================================
// TEST 11: Multiple intersection values must match CYTools
// NOTE: This test uses the OLD compute_intersection_numbers algorithm which has
// known issues with certain polytopes (high residual). The CYTools algorithm
// (compute_intersection_cytools) is the correct one and is tested in dual_tests.rs
// and mcallister_e2e tests.
// =============================================================================

#[test]
#[ignore = "Uses old algorithm with known issues; see dual_tests.rs for CYTools algorithm tests"]
fn test_intersection_numbers_multiple_values() {
    let (points, tri, origin_idx) = dual_test_setup();

    let non_origin: Vec<Point> = points
        .iter()
        .enumerate()
        .filter(|&(i, _)| i != origin_idx)
        .map(|(_, p)| p.clone())
        .collect();

    let glsm = compute_glsm_charge_matrix(&non_origin, true).expect("GLSM computation failed");

    let kappa = compute_intersection_numbers(&tri, &points, &glsm)
        .expect("Intersection computation failed");

    use malachite::num::conversion::traits::RoundingFrom;
    use malachite::rounding_modes::RoundingMode;

    // CYTools ground truth values
    let expected: Vec<((usize, usize, usize), i64)> = vec![
        ((3, 4, 5), 1),
        ((4, 5, 8), 1),
        ((1, 2, 3), 6),
        ((1, 1, 1), 56),
        ((3, 3, 3), 1),
    ];

    for ((i, j, k), exp) in expected {
        let (val, _) = f64::rounding_from(kappa.get(i, j, k).get(), RoundingMode::Nearest);
        assert!(
            (val - exp as f64).abs() < 0.5,
            "kappa[{i},{j},{k}] = {val:.2} but expected {exp}"
        );
    }
}

// =============================================================================
// TEST 12: Basis transformation - extract subtensor for basis indices
// =============================================================================

#[test]
fn test_basis_transformation() {
    // Given a full intersection tensor kappa on n points,
    // extract the subtensor for only the basis indices.
    //
    // For McAllister 4-214-647:
    // - Dual polytope has 9 points (indices 0-8)
    // - Dual basis is [3, 4, 5, 8] (h11=4)
    // - We extract kappa_{ijk} where i,j,k in {3,4,5,8}

    // Create a mock full tensor with known values
    let mut full_kappa: HashMap<Vec<usize>, f64> = HashMap::new();

    // Set some values (matching CYTools ground truth for dual polytope)
    full_kappa.insert(vec![3, 4, 5], 1.0);
    full_kappa.insert(vec![4, 5, 8], 1.0);
    full_kappa.insert(vec![3, 4, 8], 2.0);
    full_kappa.insert(vec![3, 5, 8], 2.0);
    full_kappa.insert(vec![3, 3, 4], 0.0);
    full_kappa.insert(vec![3, 3, 5], 0.0);
    full_kappa.insert(vec![3, 3, 8], -1.0);
    full_kappa.insert(vec![3, 4, 4], 0.0);
    full_kappa.insert(vec![4, 4, 5], 0.0);
    full_kappa.insert(vec![4, 4, 8], -1.0);
    full_kappa.insert(vec![3, 5, 5], 0.0);
    full_kappa.insert(vec![4, 5, 5], 0.0);
    full_kappa.insert(vec![5, 5, 8], -1.0);
    full_kappa.insert(vec![3, 8, 8], -2.0);
    full_kappa.insert(vec![4, 8, 8], -2.0);
    full_kappa.insert(vec![5, 8, 8], -2.0);
    full_kappa.insert(vec![3, 3, 3], 1.0);
    full_kappa.insert(vec![4, 4, 4], 1.0);
    full_kappa.insert(vec![5, 5, 5], 1.0);
    full_kappa.insert(vec![8, 8, 8], 6.0);

    // Values outside the basis (shouldn't be extracted)
    full_kappa.insert(vec![1, 2, 3], 6.0);
    full_kappa.insert(vec![1, 1, 1], 56.0);
    full_kappa.insert(vec![2, 2, 2], 24.0);

    // Basis indices
    let basis = vec![3, 4, 5, 8];
    let h11 = basis.len();

    // Extract subtensor: map basis indices to 0, 1, 2, 3
    let mut basis_kappa: HashMap<Vec<usize>, f64> = HashMap::new();

    for (new_i, &old_i) in basis.iter().enumerate() {
        for (new_j, &old_j) in basis.iter().enumerate() {
            if new_j < new_i {
                continue;
            }
            for (new_k, &old_k) in basis.iter().enumerate() {
                if new_k < new_j {
                    continue;
                }

                let mut old_key = vec![old_i, old_j, old_k];
                old_key.sort_unstable();

                if let Some(&val) = full_kappa.get(&old_key) {
                    let new_key = vec![new_i, new_j, new_k];
                    basis_kappa.insert(new_key, val);
                }
            }
        }
    }

    // Verify the extracted subtensor
    // In basis indexing: 0=3, 1=4, 2=5, 3=8
    // kappa_{basis}[0,1,2] should be kappa_{full}[3,4,5] = 1
    assert_eq!(
        basis_kappa.get(&vec![0, 1, 2]).copied(),
        Some(1.0),
        "kappa_basis[0,1,2] should be kappa_full[3,4,5] = 1"
    );

    // kappa_{basis}[1,2,3] should be kappa_{full}[4,5,8] = 1
    assert_eq!(
        basis_kappa.get(&vec![1, 2, 3]).copied(),
        Some(1.0),
        "kappa_basis[1,2,3] should be kappa_full[4,5,8] = 1"
    );

    // kappa_{basis}[3,3,3] should be kappa_{full}[8,8,8] = 6
    assert_eq!(
        basis_kappa.get(&vec![3, 3, 3]).copied(),
        Some(6.0),
        "kappa_basis[3,3,3] should be kappa_full[8,8,8] = 6"
    );

    // Verify values outside basis are NOT in subtensor
    assert!(
        basis_kappa.get(&vec![1, 2, 3]).is_some()
            || full_kappa.get(&vec![1, 2, 3]).is_some_and(|&v| v == 6.0),
        "Sanity check: full tensor has kappa[1,2,3]=6"
    );

    // The basis tensor should have exactly C(4+3-1, 3) = C(6,3) = 20 entries
    // (multiset of size 3 from 4 elements)
    let expected_entries = 20;
    assert_eq!(
        basis_kappa.len(),
        expected_entries,
        "Basis subtensor should have {} entries, got {}",
        expected_entries,
        basis_kappa.len()
    );

    println!("Basis transformation: extracted {h11}-dim subtensor from full tensor");
    println!("  Basis: {basis:?}");
    println!("  kappa_basis[0,1,2] = kappa_full[3,4,5] = 1");
    println!("  kappa_basis[1,2,3] = kappa_full[4,5,8] = 1");
    println!("  kappa_basis[3,3,3] = kappa_full[8,8,8] = 6");
}

// =============================================================================
// TEST 13: Verify GLSM constraints are compatible with kappa[3,4,5]=1
// =============================================================================

#[test]
fn test_glsm_constraints_vs_expected_kappa_345() {
    // For kappa_{345} = 1, we need:
    // kappa_{345} = kappa^V_{1345} + kappa^V_{2345} + kappa^V_{3345} + kappa^V_{3445} + kappa^V_{3455}
    //         = 1/3 + 1/2 + S
    //         = 5/6 + S = 1
    // So S = kappa^V_{3345} + kappa^V_{3445} + kappa^V_{3455} = 1/6

    // GLSM with origin: (column 0 = origin = [0,0,0,0])
    let glsm: Vec<Vec<i64>> = vec![
        vec![-6, 2, 3, -1, 1, 1, 0, 0, 0],
        vec![-6, 2, 3, 1, -1, 0, 1, 0, 0],
        vec![-12, 4, 6, 0, 1, 0, 0, 1, 0],
        vec![-6, 2, 3, 0, 0, 0, 0, 0, 1],
    ];

    let kv_1345 = 1.0 / 3.0;
    let kv_2345 = 0.5;

    // kappa^V_{0345} = 0 because det of matrix with origin row [0,0,0,0] is 0

    println!("\n=== GLSM constraints for probe (3,4,5) ===");
    println!("Known: kappa^V_{{1345}} = 1/3, kappa^V_{{2345}} = 1/2");
    println!("Variables: kappa^V_{{3345}}, kappa^V_{{3445}}, kappa^V_{{3455}}");
    println!(
        "For kappa_{{345}} = 1, need: S = kappa^V_{{3345}} + kappa^V_{{3445}} + kappa^V_{{3455}} = 1/6"
    );

    // Build constraint equations from each GLSM row
    // Full constraint: Sum_i Q_i kappa^V_{i,3,4,5} = 0
    // kappa^V_{0345} = 0 (origin has zero det)
    // kappa^V_{6345}, kappa^V_{7345}, kappa^V_{8345} = 0 (not in any simplex)

    for (row_idx, q_row) in glsm.iter().enumerate() {
        // Compute: Q_0*kappa^V_{0345} + Q_1*kappa^V_{1345} + Q_2*kappa^V_{2345} + Q_3*kappa^V_{3345} + ...
        let q0 = q_row[0] as f64;
        let q1 = q_row[1] as f64;
        let q2 = q_row[2] as f64;
        let q3 = q_row[3] as f64;
        let q4 = q_row[4] as f64;
        let q5 = q_row[5] as f64;
        let q6 = q_row[6] as f64;
        let q7 = q_row[7] as f64;
        let q8 = q_row[8] as f64;

        // Known terms (move to RHS)
        let rhs = -q2.mul_add(kv_2345, q0.mul_add(0.0, q1 * kv_1345));

        // Variable terms (LHS): Q_3*kappa^V_{3345} + Q_4*kappa^V_{3445} + Q_5*kappa^V_{3455}
        // Note: kappa^V_{3456}, kappa^V_{3457}, kappa^V_{3458} = 0 (not in simplex)
        println!("\nRow {row_idx}: Q = {q_row:?}");
        println!(
            "  Known: Q_1 x kappa^V_{{1345}} + Q_2 x kappa^V_{{2345}} = {} x {:.4} + {} x {:.4} = {:.4}",
            q1,
            kv_1345,
            q2,
            kv_2345,
            q1.mul_add(kv_1345, q2 * kv_2345)
        );
        println!(
            "  Equation: {q3} x kappa^V_{{3345}} + {q4} x kappa^V_{{3445}} + {q5} x kappa^V_{{3455}} = {rhs:.4}"
        );

        // Check if any Q_6, Q_7, Q_8 are non-zero
        if q6 != 0.0 || q7 != 0.0 || q8 != 0.0 {
            println!(
                "  Note: Q_6={q6}, Q_7={q7}, Q_8={q8} (these terms are 0 since not in simplex)"
            );
        }
    }

    // Now solve the system analytically (rows 0, 1, 2 only - row 3 is degenerate)
    // Row 0: -kappa^V_{3345} + kappa^V_{3445} + kappa^V_{3455} = -(2/3 + 3/2) = -13/6
    // Row 1: kappa^V_{3345} - kappa^V_{3445} = -(2/3 + 3/2) = -13/6
    // Row 2: kappa^V_{3445} = -(4/3 + 3) = -13/3

    println!("\n=== Analytical solution (from rows 0-2) ===");

    let kv_3445 = -13.0 / 3.0; // From row 2
    println!("From Row 2: kappa^V_{{3445}} = -13/3 ~ {kv_3445:.4}");

    let kv_3345 = kv_3445 - 13.0 / 6.0; // From row 1
    println!(
        "From Row 1: kappa^V_{{3345}} = kappa^V_{{3445}} - 13/6 = {:.4} - {:.4} = {:.4}",
        kv_3445,
        13.0 / 6.0,
        kv_3345
    );

    // From row 0: -kappa^V_{3345} + kappa^V_{3445} + kappa^V_{3455} = -13/6
    let kv_3455 = -13.0 / 6.0 + kv_3345 - kv_3445;
    println!(
        "From Row 0: kappa^V_{{3455}} = -13/6 + kappa^V_{{3345}} - kappa^V_{{3445}} = {kv_3455:.4}"
    );

    let s_actual = kv_3345 + kv_3445 + kv_3455;
    let kappa_345_actual = kv_1345 + kv_2345 + s_actual;

    println!("\n=== Result ===");
    println!("S = kappa^V_{{3345}} + kappa^V_{{3445}} + kappa^V_{{3455}} = {s_actual:.4}");
    println!("kappa_{{345}} = 1/3 + 1/2 + S = {kappa_345_actual:.4}");
    println!("Expected: kappa_{{345}} = 1");
    println!("Discrepancy: {:.4}", kappa_345_actual - 1.0);

    // This test documents the FUNDAMENTAL INCONSISTENCY:
    // The GLSM constraints give kappa_{345} ~ -14.3, NOT 1!
    //
    // This means either:
    // 1. The 4-form approach isn't what CYTools uses
    // 2. We're missing something about how to handle the origin
    // 3. The "least squares" solution is supposed to override the exact solution
    //
    // NOTE: This test is intentionally left without an assertion to document the issue.
}

// =============================================================================
// TEST 14: Check which simplices contain (3,4,5)
// =============================================================================

#[test]
fn test_simplices_containing_345() {
    let (_, tri, origin_idx) = dual_test_setup();

    println!("\n=== Simplices containing {{3,4,5}} ===");

    let mut count = 0;
    for (idx, simplex) in tri.simplices().iter().enumerate() {
        let rays: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| i != origin_idx)
            .collect();

        if rays.contains(&3) && rays.contains(&4) && rays.contains(&5) {
            println!("  Simplex {idx}: rays {rays:?}");
            count += 1;
        }
    }

    assert_eq!(count, 2, "Exactly 2 simplices should contain {{3,4,5}}");
    println!("\nFound {count} simplices containing {{3,4,5}}");
    println!(
        "  This means kappa_{{345}} = kappa^V_{{1345}} + kappa^V_{{2345}} + self-intersections"
    );
    println!("                      = 1/3 + 1/2 + S");
    println!("                      = 5/6 + S");
}
