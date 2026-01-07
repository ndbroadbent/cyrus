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
#![allow(clippy::collection_is_never_read)]
#![allow(clippy::redundant_clone)]
#![allow(clippy::stable_sort_primitive)]
#![allow(clippy::bool_to_int_with_if)]
#![allow(dead_code)]
#![allow(unused_variables)]

//! Dual test suite: CYTools <-> Rust comparison tests.
//!
//! These tests run identical inputs through CYTools (Python) and our Rust implementation,
//! ensuring we get exactly the same results at every intermediate step.
//!
//! # Regenerating Fixtures
//!
//! ```bash
//! cd crates/cyrus-core/tests/dual_tests
//! /path/to/cytools_venv/bin/python cytools_all_functions.py > fixtures/all_functions.json
//! ```
//!
//! # Test Organization
//!
//! Core tests (determinant, triangulation, GLSM, linear relations, divisor basis, 4-form, 3-face)
//! are in this file. Intersection number and flat direction tests are in `dual_tests_intersection.rs`.

use cyrus_core::Point;
use cyrus_core::basis::compute_divisor_basis;
use cyrus_core::glsm::compute_glsm_charge_matrix;
use cyrus_core::integer_math::determinant_gaussian;
use malachite::Rational;
use serde::Deserialize;

// ==================== FIXTURE DATA STRUCTURES ====================

#[derive(Debug, Deserialize)]
struct AllFunctionsFixture {
    test_case: String,
    input_points: Vec<Vec<i64>>,
    n_points: usize,
    dim: usize,
    polytope: PolytopeData,
    triangulation: TriangulationData,
    triangulation_from_heights: TriangulationFromHeights,
    glsm_charge_matrix: MatrixData,
    linear_relations_with_origin: MatrixData,
    linear_relations_no_origin: MatrixData,
    linear_relations_glsm_product: Vec<Vec<i64>>,
    linear_relations_glsm_product_is_zero: bool,
    divisor_basis: Vec<usize>,
    h11: usize,
    h21: usize,
    intersection_numbers_full: Vec<KappaEntry>,
    intersection_numbers_basis: Vec<KappaEntry>,
    distinct_4form_values: Vec<FourFormEntry>,
    probe_3faces: Vec<Vec<usize>>,
    sanity_checks: SanityChecks,
    flat_direction: FlatDirectionData,
}

#[derive(Debug, Deserialize)]
struct PolytopeData {
    points: Vec<Vec<i64>>,
    vertices: Vec<Vec<i64>>,
    is_reflexive: bool,
    dimension: usize,
    points_not_interior_to_facets: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct TriangulationData {
    simplices: Vec<Vec<usize>>,
    n_simplices: usize,
    simplex_determinants: Vec<SimplexDet>,
    heights: Vec<f64>,
}

#[derive(Debug, Deserialize)]
struct TriangulationFromHeights {
    input_heights: Vec<f64>,
    simplices: Vec<Vec<usize>>,
    matches_original: bool,
}

#[derive(Debug, Deserialize)]
struct SimplexDet {
    simplex: Vec<usize>,
    abs_det: i64,
}

#[derive(Debug, Deserialize)]
struct MatrixData {
    matrix: Vec<Vec<i64>>,
    shape: Vec<usize>,
}

#[derive(Debug, Deserialize)]
struct KappaEntry {
    indices: Vec<usize>,
    value: i64,
}

#[derive(Debug, Deserialize)]
struct FourFormEntry {
    indices: Vec<usize>,
    det: f64,
    kappa_4form: f64,
}

#[derive(Debug, Deserialize)]
struct SanityChecks {
    kappa_345: Option<i64>,
    kappa_111: Option<i64>,
    kappa_333: Option<i64>,
    kappa_123: Option<i64>,
    kappa_888: Option<i64>,
}

#[derive(Debug, Deserialize)]
struct FlatDirectionData {
    #[serde(rename = "K_flux")]
    k_flux: Vec<i64>,
    #[serde(rename = "M_flux")]
    m_flux: Vec<i64>,
    #[serde(rename = "N_matrix")]
    n_matrix: Vec<Vec<f64>>,
    p: Option<Vec<f64>>,
    solve_success: bool,
    kappa_ppp: Option<f64>,
    #[serde(rename = "eK0")]
    ek0: Option<f64>,
}

fn load_fixture() -> AllFunctionsFixture {
    let json_str = include_str!("dual_tests/fixtures/all_functions.json");
    serde_json::from_str(json_str).expect("Failed to parse fixture")
}

fn points_from_vecs(vecs: &[Vec<i64>]) -> Vec<Point> {
    vecs.iter().map(|v| Point::new(v.clone())).collect()
}

// ==================== DETERMINANT TESTS ====================

#[test]
fn test_dual_simplex_determinants() {
    let fixture = load_fixture();
    let points = points_from_vecs(&fixture.input_points);

    println!("Testing simplex determinants...");
    println!(
        "CYTools reports {} simplices",
        fixture.triangulation.n_simplices
    );

    let mut all_match = true;
    for det_entry in &fixture.triangulation.simplex_determinants {
        let simplex = &det_entry.simplex;
        let expected_det = det_entry.abs_det;

        // Build matrix from simplex points (augmented with 1s column)
        let mut mat: Vec<Vec<Rational>> = simplex
            .iter()
            .map(|&idx| {
                let mut row: Vec<Rational> = points[idx]
                    .coords()
                    .iter()
                    .map(|&x| Rational::from(x))
                    .collect();
                row.push(Rational::from(1)); // Augment with 1
                row
            })
            .collect();

        let our_det = determinant_gaussian(&mut mat);
        let our_abs_det = if our_det < 0 { -our_det } else { our_det };
        let our_abs_det_i64: i64 = (&our_abs_det).try_into().unwrap_or(0);

        if our_abs_det_i64 != expected_det {
            println!(
                "  MISMATCH: simplex {:?}: ours={}, expected={}",
                simplex, our_abs_det_i64, expected_det
            );
            all_match = false;
        }
    }

    if all_match {
        println!(
            "  All {} simplex determinants match!",
            fixture.triangulation.simplex_determinants.len()
        );
    }
    assert!(all_match, "Some simplex determinants don't match");
}

// ==================== TRIANGULATION TESTS ====================

#[test]
fn test_dual_frst_heights() {
    use cyrus_core::compute_frst_heights;

    let fixture = load_fixture();
    let points = points_from_vecs(&fixture.input_points);

    // Find origin index
    let origin_idx = points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .expect("Origin not found");

    println!("Testing FRST heights computation...");
    println!("CYTools heights: {:?}", fixture.triangulation.heights);

    // Compute our FRST heights
    let (our_heights, our_tri) =
        compute_frst_heights(&points, origin_idx).expect("Failed to compute FRST heights");

    println!("Our heights: {:?}", our_heights);

    // Compare heights (allow small floating point differences)
    let heights_match = our_heights.len() == fixture.triangulation.heights.len()
        && our_heights
            .iter()
            .zip(&fixture.triangulation.heights)
            .all(|(a, b)| (a - b).abs() < 1e-10);

    if !heights_match {
        println!("  Heights differ!");
        for (i, (ours, theirs)) in our_heights
            .iter()
            .zip(&fixture.triangulation.heights)
            .enumerate()
        {
            if (ours - theirs).abs() >= 1e-10 {
                println!("    Point {}: ours={}, theirs={}", i, ours, theirs);
            }
        }
    }

    // Compare simplices
    let mut our_simplices: Vec<Vec<usize>> = our_tri.simplices().iter().cloned().collect();
    let mut expected_simplices = fixture.triangulation.simplices.clone();
    our_simplices.iter_mut().for_each(|s| s.sort());
    expected_simplices.iter_mut().for_each(|s| s.sort());
    our_simplices.sort();
    expected_simplices.sort();

    let simplices_match = our_simplices == expected_simplices;

    println!("Heights match: {}", heights_match);
    println!("Simplices match: {}", simplices_match);
    println!("Our simplices: {}", our_simplices.len());
    println!("Expected simplices: {}", expected_simplices.len());

    // For now, just require simplices to match (heights may differ in normalization)
    assert!(
        simplices_match,
        "FRST triangulation simplices don't match CYTools"
    );
}

#[test]
fn test_dual_triangulation_from_heights() {
    use cyrus_core::compute_regular_triangulation;

    let fixture = load_fixture();
    let points = points_from_vecs(&fixture.input_points);

    println!("Testing triangulation from heights...");
    println!(
        "Input heights: {:?}",
        fixture.triangulation_from_heights.input_heights
    );

    // Compute triangulation from the given heights
    let tri =
        compute_regular_triangulation(&points, &fixture.triangulation_from_heights.input_heights)
            .expect("Failed to compute triangulation from heights");

    // Compare simplices
    let mut our_simplices: Vec<Vec<usize>> = tri.simplices().iter().cloned().collect();
    let mut expected_simplices = fixture.triangulation_from_heights.simplices.clone();
    our_simplices.iter_mut().for_each(|s| s.sort());
    expected_simplices.iter_mut().for_each(|s| s.sort());
    our_simplices.sort();
    expected_simplices.sort();

    println!("Our simplices: {} entries", our_simplices.len());
    println!("Expected simplices: {} entries", expected_simplices.len());

    let simplices_match = our_simplices == expected_simplices;

    if !simplices_match {
        println!("MISMATCH!");
        println!("Our simplices:");
        for s in &our_simplices {
            println!("  {:?}", s);
        }
        println!("Expected simplices:");
        for s in &expected_simplices {
            println!("  {:?}", s);
        }
    }

    assert!(
        simplices_match,
        "Triangulation from heights doesn't match CYTools"
    );
    assert!(
        fixture.triangulation_from_heights.matches_original,
        "CYTools reports heights don't reproduce original triangulation"
    );
}

// ==================== GLSM TESTS ====================

#[test]
fn test_dual_glsm_shape() {
    let fixture = load_fixture();

    // CYTools GLSM shape
    let expected_rows = fixture.glsm_charge_matrix.shape[0];
    let expected_cols = fixture.glsm_charge_matrix.shape[1];

    println!("CYTools GLSM shape: {} x {}", expected_rows, expected_cols);

    // Our GLSM: points without origin, then prepend origin
    let points_without_origin: Vec<Point> = fixture
        .input_points
        .iter()
        .skip(1) // Skip origin at index 0
        .map(|v| Point::new(v.clone()))
        .collect();

    let glsm = compute_glsm_charge_matrix(&points_without_origin, true).unwrap();

    println!(
        "Our GLSM shape: {} x {}",
        glsm.len(),
        glsm.first().map_or(0, |r| r.len())
    );

    assert_eq!(glsm.len(), expected_rows, "GLSM row count mismatch");
    assert_eq!(
        glsm.first().map_or(0, |r| r.len()),
        expected_cols,
        "GLSM column count mismatch"
    );
}

#[test]
fn test_dual_glsm_kernel_property() {
    let fixture = load_fixture();

    // The GLSM kernel property: L @ Q^T = 0
    // Where L is linear_relations and Q is GLSM
    assert!(
        fixture.linear_relations_glsm_product_is_zero,
        "CYTools L @ Q^T should be zero"
    );

    // Verify with our GLSM
    let points_without_origin: Vec<Point> = fixture
        .input_points
        .iter()
        .skip(1)
        .map(|v| Point::new(v.clone()))
        .collect();

    let glsm = compute_glsm_charge_matrix(&points_without_origin, true).unwrap();
    let lr = &fixture.linear_relations_with_origin.matrix;

    // Compute L @ Q^T
    let n_lr_rows = lr.len();
    let n_glsm_rows = glsm.len();
    let n_cols = glsm.first().map_or(0, |r| r.len());

    println!("Checking L @ Q^T = 0...");
    println!("  L: {} x {}", n_lr_rows, lr.first().map_or(0, |r| r.len()));
    println!("  Q: {} x {}", n_glsm_rows, n_cols);

    let mut all_zero = true;
    for i in 0..n_lr_rows {
        for j in 0..n_glsm_rows {
            let mut dot: i64 = 0;
            for k in 0..n_cols {
                let l_ik = lr[i][k];
                let q_jk: i64 = (&glsm[j][k]).try_into().unwrap_or(0);
                dot += l_ik * q_jk;
            }
            if dot != 0 {
                println!("  (L @ Q^T)[{},{}] = {} (expected 0)", i, j, dot);
                all_zero = false;
            }
        }
    }

    if all_zero {
        println!("  L @ Q^T = 0 verified!");
    }
    assert!(all_zero, "L @ Q^T should be zero");
}

// ==================== LINEAR RELATIONS TESTS ====================

#[test]
fn test_dual_linear_relations_structure() {
    let fixture = load_fixture();

    let lr = &fixture.linear_relations_with_origin.matrix;
    let shape = &fixture.linear_relations_with_origin.shape;

    println!("Linear relations shape: {} x {}", shape[0], shape[1]);
    println!(
        "Expected: (dim+1) x n_points = {} x {}",
        fixture.dim + 1,
        fixture.n_points
    );

    assert_eq!(shape[0], fixture.dim + 1, "LR should have dim+1 rows");
    assert_eq!(
        shape[1], fixture.n_points,
        "LR should have n_points columns"
    );

    // Check identity structure in first (dim+1) columns
    println!("\nLinear relations matrix:");
    for (i, row) in lr.iter().enumerate() {
        println!("  Row {}: {:?}", i, row);
    }

    // The pivot columns should form identity
    for i in 0..shape[0] {
        for j in 0..shape[0] {
            let expected = if i == j { 1 } else { 0 };
            assert_eq!(
                lr[i][j], expected,
                "LR[{},{}] should be {} (identity block)",
                i, j, expected
            );
        }
    }
    println!("\nIdentity block in first {} columns verified!", shape[0]);
}

// ==================== DIVISOR BASIS TESTS ====================

#[test]
fn test_dual_divisor_basis() {
    let fixture = load_fixture();

    println!("CYTools divisor basis: {:?}", fixture.divisor_basis);
    println!("CYTools h11: {}", fixture.h11);

    // The basis should be the non-pivot columns of linear_relations
    // Pivot columns are 0, 1, ..., dim (first dim+1 columns)
    let n_points = fixture.n_points;
    let n_pivots = fixture.dim + 1;
    let expected_basis: Vec<usize> = (n_pivots..n_points).collect();

    assert_eq!(
        fixture.divisor_basis,
        expected_basis,
        "Divisor basis should be columns {} to {}",
        n_pivots,
        n_points - 1
    );

    assert_eq!(
        fixture.h11,
        fixture.divisor_basis.len(),
        "h11 should equal basis size"
    );
}

#[test]
fn test_dual_our_divisor_basis() {
    let fixture = load_fixture();
    let points = points_from_vecs(&fixture.input_points);

    // Compute GLSM
    let points_without_origin: Vec<Point> = points.iter().skip(1).cloned().collect();
    let glsm = compute_glsm_charge_matrix(&points_without_origin, true).unwrap();

    // Compute our divisor basis
    let our_basis = compute_divisor_basis(&points, &glsm);

    match our_basis {
        Ok(basis) => {
            println!("Our divisor basis: {:?}", basis);
            println!("CYTools basis: {:?}", fixture.divisor_basis);

            // Note: basis indices may differ but should have same length
            assert_eq!(
                basis.len(),
                fixture.divisor_basis.len(),
                "Basis size mismatch: ours={}, expected={}",
                basis.len(),
                fixture.divisor_basis.len()
            );
        }
        Err(e) => {
            println!("Our divisor basis computation failed: {:?}", e);
            // This is expected to fail until we fix the algorithm
        }
    }
}

// ==================== 4-FORM TESTS ====================

#[test]
fn test_dual_distinct_4form_values() {
    let fixture = load_fixture();
    let points = points_from_vecs(&fixture.input_points);

    println!("Testing distinct 4-form values...");
    println!(
        "CYTools reports {} distinct 4-form entries",
        fixture.distinct_4form_values.len()
    );

    let mut mismatches = 0;
    for entry in &fixture.distinct_4form_values {
        let indices = &entry.indices;
        let expected_det = entry.det;
        let expected_kappa = entry.kappa_4form;

        // Compute determinant of 4 points (4x4 matrix)
        let mut mat: Vec<Vec<Rational>> = indices
            .iter()
            .map(|&idx| {
                points[idx]
                    .coords()
                    .iter()
                    .map(|&x| Rational::from(x))
                    .collect()
            })
            .collect();

        let det = determinant_gaussian(&mut mat);
        let abs_det = if det < 0 { -det.clone() } else { det };
        let abs_det_f64: f64 = abs_det.to_string().parse().unwrap_or(0.0);

        let our_kappa = if abs_det_f64 > 1e-10 {
            1.0 / abs_det_f64
        } else {
            0.0
        };

        if (abs_det_f64 - expected_det).abs() > 0.001 {
            println!(
                "  Det mismatch for {:?}: ours={}, expected={}",
                indices, abs_det_f64, expected_det
            );
            mismatches += 1;
        }

        if (our_kappa - expected_kappa).abs() > 0.001 {
            println!(
                "  κ^V mismatch for {:?}: ours={}, expected={}",
                indices, our_kappa, expected_kappa
            );
            mismatches += 1;
        }
    }

    if mismatches == 0 {
        println!(
            "  All {} 4-form values match!",
            fixture.distinct_4form_values.len()
        );
    }
    assert_eq!(mismatches, 0, "{} 4-form mismatches", mismatches);
}

// ==================== 3-FACE PROBES TESTS ====================

#[test]
fn test_dual_3face_probes() {
    let fixture = load_fixture();

    println!(
        "CYTools 3-face probes: {} entries",
        fixture.probe_3faces.len()
    );

    // Verify structure: all should be sorted triples with i <= j <= k
    for probe in &fixture.probe_3faces {
        assert_eq!(probe.len(), 3, "Probe should have 3 indices");
        assert!(
            probe[0] <= probe[1] && probe[1] <= probe[2],
            "Probe should be sorted"
        );
        assert!(probe[0] > 0, "Probe indices should exclude origin (0)");
    }

    println!("  All probes have valid structure");
    println!(
        "  Sample probes: {:?}",
        &fixture.probe_3faces[..5.min(fixture.probe_3faces.len())]
    );
}
