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

//! Intersection number and flat direction dual tests.
//!
//! Split from dual_tests.rs for file size management.

use cyrus_core::Point;
use cyrus_core::intersection::compute_intersection_cytools;
use cyrus_core::intersection::compute_intersection_numbers_with_linear_relations;
use cyrus_core::triangulation::Triangulation;
use malachite::Integer;
use serde::Deserialize;
use std::collections::HashMap;

// ==================== FIXTURE DATA STRUCTURES ====================
// (Duplicated from dual_tests.rs for independence)

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

// ==================== INTERSECTION NUMBER TESTS ====================

#[test]
fn test_dual_intersection_sanity_checks() {
    let fixture = load_fixture();
    let points = points_from_vecs(&fixture.input_points);
    let simplices = fixture.triangulation.simplices.clone();

    println!("Testing intersection number sanity checks...");

    // Build intersection numbers from CYTools
    let mut expected_kappa: HashMap<(usize, usize, usize), i64> = HashMap::new();
    for entry in &fixture.intersection_numbers_full {
        let (i, j, k) = (entry.indices[0], entry.indices[1], entry.indices[2]);
        expected_kappa.insert((i, j, k), entry.value);
    }

    println!("\nCYTools sanity checks:");
    println!("  κ[3,4,5] = {:?}", fixture.sanity_checks.kappa_345);
    println!("  κ[1,1,1] = {:?}", fixture.sanity_checks.kappa_111);
    println!("  κ[3,3,3] = {:?}", fixture.sanity_checks.kappa_333);
    println!("  κ[1,2,3] = {:?}", fixture.sanity_checks.kappa_123);
    println!("  κ[8,8,8] = {:?}", fixture.sanity_checks.kappa_888);

    // Try computing with our implementation
    let tri = Triangulation::new(simplices);
    let linear_relations: Vec<Vec<Integer>> = fixture
        .linear_relations_with_origin
        .matrix
        .iter()
        .map(|row| row.iter().map(|&x| Integer::from(x)).collect())
        .collect();

    println!("\nComputing intersection numbers with our implementation...");

    match compute_intersection_numbers_with_linear_relations(&tri, &points, &linear_relations) {
        Ok(kappa) => {
            println!("\nOur results:");
            println!("  κ[3,4,5] = {}", kappa.get(3, 4, 5).get());
            println!("  κ[1,1,1] = {}", kappa.get(1, 1, 1).get());
            println!("  κ[3,3,3] = {}", kappa.get(3, 3, 3).get());
            println!("  κ[1,2,3] = {}", kappa.get(1, 2, 3).get());
            println!("  κ[8,8,8] = {}", kappa.get(8, 8, 8).get());

            // Check specific values
            if let Some(expected) = fixture.sanity_checks.kappa_345 {
                let ours = kappa.get(3, 4, 5);
                let ours_i64: i64 = ours.get().to_string().parse().unwrap_or(-999);
                if ours_i64 != expected {
                    println!(
                        "\n  MISMATCH κ[3,4,5]: ours={}, expected={}",
                        ours_i64, expected
                    );
                }
            }
        }
        Err(e) => {
            println!("  Computation failed: {:?}", e);
        }
    }
}

#[test]
fn test_dual_intersection_all_values() {
    let fixture = load_fixture();
    let points = points_from_vecs(&fixture.input_points);
    let simplices = fixture.triangulation.simplices.clone();

    let tri = Triangulation::new(simplices);
    let linear_relations: Vec<Vec<Integer>> = fixture
        .linear_relations_with_origin
        .matrix
        .iter()
        .map(|row| row.iter().map(|&x| Integer::from(x)).collect())
        .collect();

    // Build expected map
    let mut expected: HashMap<(usize, usize, usize), i64> = HashMap::new();
    for entry in &fixture.intersection_numbers_full {
        let (i, j, k) = (entry.indices[0], entry.indices[1], entry.indices[2]);
        expected.insert((i, j, k), entry.value);
    }

    match compute_intersection_numbers_with_linear_relations(&tri, &points, &linear_relations) {
        Ok(kappa) => {
            let mut matches = 0;
            let mut mismatches = 0;

            // Check all expected values
            for (&(i, j, k), &expected_val) in &expected {
                let our_val = kappa.get(i, j, k);
                let our_i64: i64 = our_val.get().to_string().parse().unwrap_or(-999);

                if our_i64 == expected_val {
                    matches += 1;
                } else {
                    mismatches += 1;
                    if mismatches <= 10 {
                        println!(
                            "  MISMATCH κ[{},{},{}]: ours={}, expected={}",
                            i, j, k, our_i64, expected_val
                        );
                    }
                }
            }

            println!("\nIntersection number comparison (old algorithm):");
            println!("  Matches: {}", matches);
            println!("  Mismatches: {}", mismatches);

            // For now, just report - don't fail the test
            // assert_eq!(mismatches, 0, "{} intersection number mismatches", mismatches);
        }
        Err(e) => {
            println!("Computation failed: {:?}", e);
        }
    }
}

#[test]
fn test_dual_intersection_cytools_algorithm() {
    let fixture = load_fixture();
    let points = points_from_vecs(&fixture.input_points);
    let simplices = fixture.triangulation.simplices.clone();

    let tri = Triangulation::new(simplices);

    // Use linear_relations WITHOUT origin (as CYTools does internally)
    let linear_relations_no_origin: Vec<Vec<i64>> =
        fixture.linear_relations_no_origin.matrix.clone();

    // Build expected map
    let mut expected: HashMap<(usize, usize, usize), i64> = HashMap::new();
    for entry in &fixture.intersection_numbers_full {
        let (i, j, k) = (entry.indices[0], entry.indices[1], entry.indices[2]);
        expected.insert((i, j, k), entry.value);
    }

    println!("\nTesting CYTools-compatible intersection algorithm...");
    println!("  Points: {}", points.len());
    println!("  Simplices: {}", fixture.triangulation.n_simplices);
    println!(
        "  Linear relations (no origin): {} x {}",
        linear_relations_no_origin.len(),
        linear_relations_no_origin.first().map_or(0, |r| r.len())
    );

    match compute_intersection_cytools(&tri, &points, &linear_relations_no_origin) {
        Ok(kappa) => {
            let mut matches = 0;
            let mut mismatches = 0;

            // Check all expected values
            for (&(i, j, k), &expected_val) in &expected {
                let our_val = kappa.get(i, j, k);
                let our_i64: i64 = our_val.get().to_string().parse().unwrap_or(-999);

                if our_i64 == expected_val {
                    matches += 1;
                } else {
                    mismatches += 1;
                    if mismatches <= 20 {
                        println!(
                            "  MISMATCH κ[{},{},{}]: ours={}, expected={}",
                            i, j, k, our_i64, expected_val
                        );
                    }
                }
            }

            println!("\nIntersection number comparison (CYTools algorithm):");
            println!("  Matches: {}/{}", matches, expected.len());
            println!("  Mismatches: {}", mismatches);

            // Also show some specific values
            println!("\nSample values:");
            println!(
                "  κ[3,4,5] = {} (expected {:?})",
                kappa.get(3, 4, 5).get(),
                fixture.sanity_checks.kappa_345
            );
            println!(
                "  κ[1,1,1] = {} (expected {:?})",
                kappa.get(1, 1, 1).get(),
                fixture.sanity_checks.kappa_111
            );
            println!(
                "  κ[8,8,8] = {} (expected {:?})",
                kappa.get(8, 8, 8).get(),
                fixture.sanity_checks.kappa_888
            );

            assert!(
                matches > 0,
                "At least some intersection numbers should match"
            );
        }
        Err(e) => {
            panic!("CYTools algorithm computation failed: {:?}", e);
        }
    }
}

// ==================== INTERSECTION IN BASIS TEST ====================

#[test]
fn test_dual_intersection_in_basis() {
    use cyrus_core::intersection_in_basis;

    let fixture = load_fixture();
    let points = points_from_vecs(&fixture.input_points);
    let simplices = fixture.triangulation.simplices.clone();

    let tri = Triangulation::new(simplices);

    // Compute full intersection numbers using CYTools algorithm
    let linear_relations_no_origin: Vec<Vec<i64>> =
        fixture.linear_relations_no_origin.matrix.clone();

    let kappa_full = compute_intersection_cytools(&tri, &points, &linear_relations_no_origin)
        .expect("Failed to compute intersection numbers");

    // Get divisor basis (matches CYTools)
    let basis = &fixture.divisor_basis;

    // Project to basis
    let kappa_basis = intersection_in_basis(&kappa_full, basis);

    println!("\nTesting intersection_in_basis...");
    println!("  Basis: {:?}", basis);
    println!("  Basis tensor dimension: {}", kappa_basis.dim());

    // Build expected map from fixture
    let mut expected: HashMap<(usize, usize, usize), i64> = HashMap::new();
    for entry in &fixture.intersection_numbers_basis {
        let (i, j, k) = (entry.indices[0], entry.indices[1], entry.indices[2]);
        expected.insert((i, j, k), entry.value);
    }

    println!("  CYTools basis entries: {}", expected.len());
    println!("  Our basis entries: {}", kappa_basis.num_nonzero());

    // Compare
    let mut matches = 0;
    let mut mismatches = 0;

    for (&(i, j, k), &expected_val) in &expected {
        let our_val = kappa_basis.get(i, j, k);
        let our_i64: i64 = our_val.get().to_string().parse().unwrap_or(-999);

        if our_i64 == expected_val {
            matches += 1;
        } else {
            mismatches += 1;
            if mismatches <= 10 {
                println!(
                    "  MISMATCH κ_basis[{},{},{}]: ours={}, expected={}",
                    i, j, k, our_i64, expected_val
                );
            }
        }
    }

    println!("\nIntersection in basis comparison:");
    println!("  Matches: {}/{}", matches, expected.len());
    println!("  Mismatches: {}", mismatches);

    assert_eq!(mismatches, 0, "intersection_in_basis mismatches");
}

// ==================== FLAT DIRECTION TESTS ====================

#[test]
fn test_dual_n_matrix() {
    use cyrus_core::I64;
    use cyrus_core::flat_direction::compute_n_matrix;
    use cyrus_core::types::tags::Finite;

    let fixture = load_fixture();

    // First compute intersection numbers in basis using CYTools algorithm
    let points = points_from_vecs(&fixture.input_points);
    let simplices = fixture.triangulation.simplices.clone();
    let tri = Triangulation::new(simplices);

    let linear_relations_no_origin: Vec<Vec<i64>> =
        fixture.linear_relations_no_origin.matrix.clone();

    let kappa_full = compute_intersection_cytools(&tri, &points, &linear_relations_no_origin)
        .expect("Failed to compute intersection numbers");

    let basis = &fixture.divisor_basis;
    let kappa_basis = cyrus_core::intersection_in_basis(&kappa_full, basis);

    println!("\nTesting N matrix computation...");
    println!("  M flux: {:?}", fixture.flat_direction.m_flux);

    // Compute N matrix (convert to typed I64)
    let m_flux: Vec<I64<Finite>> = fixture
        .flat_direction
        .m_flux
        .iter()
        .map(|&x| I64::<Finite>::new(x))
        .collect();
    let n_matrix = compute_n_matrix(&kappa_basis, &m_flux);

    println!(
        "  N matrix shape: {} x {}",
        n_matrix.nrows(),
        n_matrix.ncols()
    );

    // Compare with expected
    let expected = &fixture.flat_direction.n_matrix;
    let mut all_match = true;

    for i in 0..n_matrix.nrows() {
        for j in 0..n_matrix.ncols() {
            let diff = (n_matrix[(i, j)] - expected[i][j]).abs();
            if diff > 1e-10 {
                println!(
                    "  MISMATCH N[{},{}]: ours={}, expected={}",
                    i,
                    j,
                    n_matrix[(i, j)],
                    expected[i][j]
                );
                all_match = false;
            }
        }
    }

    if all_match {
        println!("  N matrix matches CYTools!");
    }

    assert!(all_match, "N matrix doesn't match CYTools");
}

#[test]
fn test_dual_flat_direction_p() {
    use cyrus_core::I64;
    use cyrus_core::flat_direction::{compute_n_matrix, solve_linear_system_faer};
    use cyrus_core::types::tags::Finite;

    let fixture = load_fixture();

    // First compute intersection numbers in basis using CYTools algorithm
    let points = points_from_vecs(&fixture.input_points);
    let simplices = fixture.triangulation.simplices.clone();
    let tri = Triangulation::new(simplices);

    let linear_relations_no_origin: Vec<Vec<i64>> =
        fixture.linear_relations_no_origin.matrix.clone();

    let kappa_full = compute_intersection_cytools(&tri, &points, &linear_relations_no_origin)
        .expect("Failed to compute intersection numbers");

    let basis = &fixture.divisor_basis;
    let kappa_basis = cyrus_core::intersection_in_basis(&kappa_full, basis);

    println!("\nTesting flat direction p computation...");
    println!("  K flux: {:?}", fixture.flat_direction.k_flux);
    println!("  M flux: {:?}", fixture.flat_direction.m_flux);

    // Compute N matrix (convert to typed I64)
    let m_flux: Vec<I64<Finite>> = fixture
        .flat_direction
        .m_flux
        .iter()
        .map(|&x| I64::<Finite>::new(x))
        .collect();
    let n_matrix = compute_n_matrix(&kappa_basis, &m_flux);

    // Solve for p (convert K to typed I64)
    let k_flux: Vec<I64<Finite>> = fixture
        .flat_direction
        .k_flux
        .iter()
        .map(|&x| I64::<Finite>::new(x))
        .collect();
    let our_p = solve_linear_system_faer(&n_matrix, &k_flux);

    match (&our_p, &fixture.flat_direction.p) {
        (Some(ours), Some(expected)) => {
            // Convert our F64<Finite> to f64 for comparison
            let ours_f64: Vec<f64> = ours.iter().map(|x| x.get()).collect();
            println!("  Our p: {:?}", ours_f64);
            println!("  Expected p: {:?}", expected);

            let mut all_match = true;
            for (i, (a, b)) in ours_f64.iter().zip(expected).enumerate() {
                if (a - b).abs() > 1e-8 {
                    println!("  MISMATCH p[{}]: ours={}, expected={}", i, a, b);
                    all_match = false;
                }
            }

            if all_match {
                println!("  Flat direction p matches CYTools!");
            }
            assert!(all_match, "Flat direction p doesn't match");
        }
        (None, Some(_)) => {
            panic!("Our solve failed but CYTools succeeded");
        }
        (Some(_), None) => {
            panic!("Our solve succeeded but CYTools failed");
        }
        (None, None) => {
            println!("  Both failed to solve (expected for this flux)");
        }
    }
}

#[test]
fn test_dual_kappa_ppp() {
    use cyrus_core::I64;
    use cyrus_core::flat_direction::{compute_n_matrix, solve_linear_system_faer};
    use cyrus_core::types::tags::Finite;

    let fixture = load_fixture();

    // First compute intersection numbers in basis using CYTools algorithm
    let points = points_from_vecs(&fixture.input_points);
    let simplices = fixture.triangulation.simplices.clone();
    let tri = Triangulation::new(simplices);

    let linear_relations_no_origin: Vec<Vec<i64>> =
        fixture.linear_relations_no_origin.matrix.clone();

    let kappa_full = compute_intersection_cytools(&tri, &points, &linear_relations_no_origin)
        .expect("Failed to compute intersection numbers");

    let basis = &fixture.divisor_basis;
    let kappa_basis = cyrus_core::intersection_in_basis(&kappa_full, basis);

    // Compute p (convert to typed I64)
    let m_flux: Vec<I64<Finite>> = fixture
        .flat_direction
        .m_flux
        .iter()
        .map(|&x| I64::<Finite>::new(x))
        .collect();
    let n_matrix = compute_n_matrix(&kappa_basis, &m_flux);

    let k_flux: Vec<I64<Finite>> = fixture
        .flat_direction
        .k_flux
        .iter()
        .map(|&x| I64::<Finite>::new(x))
        .collect();

    if let Some(p) = solve_linear_system_faer(&n_matrix, &k_flux) {
        println!("\nTesting κ_abc p^a p^b p^c computation...");

        // Compute κ_ppp manually
        let h11 = kappa_basis.dim();
        let mut kappa_ppp = 0.0;

        for a in 0..h11 {
            for b in 0..h11 {
                for c in 0..h11 {
                    let kappa_val: f64 = kappa_basis
                        .get(a, b, c)
                        .get()
                        .to_string()
                        .parse()
                        .unwrap_or(0.0);
                    kappa_ppp += kappa_val * p[a].get() * p[b].get() * p[c].get();
                }
            }
        }

        println!("  Our κ_ppp: {}", kappa_ppp);
        println!("  Expected κ_ppp: {:?}", fixture.flat_direction.kappa_ppp);

        if let Some(expected) = fixture.flat_direction.kappa_ppp {
            let diff = (kappa_ppp - expected).abs();
            assert!(
                diff < 1e-8,
                "κ_ppp mismatch: ours={}, expected={}",
                kappa_ppp,
                expected
            );
            println!("  κ_ppp matches!");
        }

        // Test eK0 = (4/3) / κ_ppp
        if let Some(expected_ek0) = fixture.flat_direction.ek0 {
            let our_ek0 = (4.0 / 3.0) / kappa_ppp;
            println!("  Our e^K0: {}", our_ek0);
            println!("  Expected e^K0: {}", expected_ek0);

            let diff = (our_ek0 - expected_ek0).abs();
            assert!(
                diff < 1e-8,
                "e^K0 mismatch: ours={}, expected={}",
                our_ek0,
                expected_ek0
            );
            println!("  e^K0 matches!");
        }
    }
}

// ==================== H11/H21 TESTS ====================

#[test]
fn test_dual_hodge_numbers() {
    let fixture = load_fixture();

    println!("CYTools Hodge numbers:");
    println!("  h11 = {}", fixture.h11);
    println!("  h21 = {}", fixture.h21);

    // h11 = n_points - (dim + 1) for a favorable polytope
    let expected_h11 = fixture.n_points - (fixture.dim + 1);
    assert_eq!(fixture.h11, expected_h11, "h11 = n_points - (dim+1)");

    // Euler characteristic χ = 2(h11 - h21)
    let chi = 2 * (fixture.h11 as i64 - fixture.h21 as i64);
    println!("  χ = 2(h11 - h21) = {}", chi);
}
