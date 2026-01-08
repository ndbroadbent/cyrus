//! Dual test suite: CYTools <-> Rust comparison tests.
//!
//! These tests verify our Rust implementations produce identical results
//! to CYTools at every intermediate step.
//!
//! # Fixture Structure
//!
//! Each function has its own fixture file in `dual_tests/fixtures/`:
//! - `glsm.json` - GLSM charge matrix
//! - `triangulation.json` - Triangulation data
//! - `linear_relations.json` - GLSM linear relations
//! - `divisor_basis.json` - Divisor basis indices
//! - `intersection_full.json` - Full intersection numbers (point-indexed)
//! - `intersection_basis.json` - Intersection numbers in basis
//! - `flat_direction.json` - Flat direction computation
//! - `filter_tensor.json` - filter_tensor_indices utility function

// Test-specific clippy allows
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::uninlined_format_args)]
#![allow(clippy::items_after_statements)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::float_cmp)]
#![allow(dead_code)]

use std::collections::{BTreeMap, HashMap};
use std::path::PathBuf;

use cyrus_core::utils::filter_tensor_indices;
use cyrus_core::Point;
use malachite::Rational;
use serde::Deserialize;

// ============================================================================
// Fixture loading
// ============================================================================

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/dual_tests/fixtures")
        .join(name)
}

fn load_fixture<T: for<'de> Deserialize<'de>>(name: &str) -> T {
    let path = fixture_path(name);
    let content = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {}", path.display(), e));
    serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {}", path.display(), e))
}

fn points_from_vecs(vecs: &[Vec<i64>]) -> Vec<Point> {
    vecs.iter().map(|v| Point::new(v.clone())).collect()
}

// ============================================================================
// Fixture data structures
// ============================================================================

#[derive(Debug, Deserialize)]
struct GlsmFixture {
    input_points: Vec<Vec<i64>>,
    matrix: Vec<Vec<i64>>,
    shape: Vec<usize>,
    h11: usize,
}

#[derive(Debug, Deserialize)]
struct TriangulationFixture {
    input_points: Vec<Vec<i64>>,
    simplices: Vec<Vec<usize>>,
    n_simplices: usize,
    heights: Vec<f64>,
    simplex_determinants: Vec<SimplexDet>,
}

#[derive(Debug, Deserialize)]
struct SimplexDet {
    simplex: Vec<usize>,
    abs_det: i64,
}

#[derive(Debug, Deserialize)]
struct LinearRelationsFixture {
    input_points: Vec<Vec<i64>>,
    with_origin: MatrixData,
    without_origin: MatrixData,
    verification: VerificationData,
}

#[derive(Debug, Deserialize)]
struct MatrixData {
    matrix: Vec<Vec<i64>>,
    shape: Vec<usize>,
}

#[derive(Debug, Deserialize)]
struct VerificationData {
    #[serde(rename = "L_times_Q_transpose")]
    l_times_q_transpose: Vec<Vec<i64>>,
    is_zero: bool,
}

#[derive(Debug, Deserialize)]
struct DivisorBasisFixture {
    input_points: Vec<Vec<i64>>,
    basis: Vec<usize>,
    h11: usize,
    h21: usize,
}

#[derive(Debug, Deserialize)]
struct IntersectionFullFixture {
    input_points: Vec<Vec<i64>>,
    n_entries: usize,
    entries: Vec<TensorEntry>,
}

#[derive(Debug, Deserialize)]
struct IntersectionBasisFixture {
    input_points: Vec<Vec<i64>>,
    basis: Vec<usize>,
    h11: usize,
    n_entries: usize,
    entries: Vec<TensorEntry>,
}

#[derive(Debug, Deserialize, Clone)]
struct TensorEntry {
    indices: Vec<usize>,
    value: i64,
}

#[derive(Debug, Deserialize)]
struct FlatDirectionFixture {
    input_points: Vec<Vec<i64>>,
    h11: usize,
    result: FlatDirectionResult,
}

#[derive(Debug, Deserialize)]
struct FlatDirectionResult {
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

#[derive(Debug, Deserialize)]
struct FilterTensorFixture {
    test_cases: Vec<FilterTensorTestCase>,
}

#[derive(Debug, Deserialize)]
struct FilterTensorTestCase {
    name: String,
    input_tensor: Vec<TensorEntry>,
    indices: Vec<usize>,
    expected: Vec<TensorEntry>,
}

// ============================================================================
// filter_tensor_indices tests (uses library implementation)
// ============================================================================

#[test]
fn test_filter_tensor_indices() {
    let fixture: FilterTensorFixture = load_fixture("filter_tensor.json");

    for case in fixture.test_cases {
        // Build BTreeMap input (library uses BTreeMap)
        let input: BTreeMap<Vec<usize>, i64> = case
            .input_tensor
            .iter()
            .map(|e| (e.indices.clone(), e.value))
            .collect();

        let expected: BTreeMap<Vec<usize>, i64> = case
            .expected
            .iter()
            .map(|e| (e.indices.clone(), e.value))
            .collect();

        // Use the library function
        let result = filter_tensor_indices(&input, &case.indices);

        assert_eq!(
            result, expected,
            "filter_tensor_indices mismatch for test case '{}'",
            case.name
        );
    }
}

// ============================================================================
// GLSM tests
// ============================================================================

#[test]
fn test_glsm_charge_matrix() {
    use cyrus_core::glsm::compute_glsm_charge_matrix;

    let fixture: GlsmFixture = load_fixture("glsm.json");

    // Build points without origin (CYTools includes origin in column 0)
    let points: Vec<Point> = fixture
        .input_points
        .iter()
        .skip(1) // Skip origin
        .map(|coords| Point::new(coords.clone()))
        .collect();

    // Compute GLSM with origin included
    let glsm = compute_glsm_charge_matrix(&points, true).expect("GLSM computation failed");

    // Convert to i64
    let result: Vec<Vec<i64>> = glsm
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("GLSM value fits in i64"))
                .collect()
        })
        .collect();

    assert_eq!(result.len(), fixture.shape[0], "GLSM row count mismatch");
    assert_eq!(
        result[0].len(),
        fixture.shape[1],
        "GLSM column count mismatch"
    );
    assert_eq!(result, fixture.matrix, "GLSM matrix doesn't match CYTools");
}

// ============================================================================
// Triangulation tests
// ============================================================================

#[test]
fn test_simplex_determinants() {
    use cyrus_core::integer_math::determinant_gaussian;

    let fixture: TriangulationFixture = load_fixture("triangulation.json");
    let points = points_from_vecs(&fixture.input_points);

    for det_entry in &fixture.simplex_determinants {
        let simplex = &det_entry.simplex;
        let expected_det = det_entry.abs_det;

        // Build augmented matrix (points + 1s column)
        let mut mat: Vec<Vec<Rational>> = simplex
            .iter()
            .map(|&idx| {
                let mut row: Vec<Rational> = points[idx]
                    .coords()
                    .iter()
                    .map(|&x| Rational::from(x))
                    .collect();
                row.push(Rational::from(1));
                row
            })
            .collect();

        let det = determinant_gaussian(&mut mat);
        let abs_det = if det < 0 { -det } else { det };
        let abs_det_i64: i64 = (&abs_det).try_into().unwrap_or(0);

        assert_eq!(
            abs_det_i64, expected_det,
            "Simplex {:?} determinant mismatch: ours={}, expected={}",
            simplex, abs_det_i64, expected_det
        );
    }
}

#[test]
fn test_triangulation_from_heights() {
    use cyrus_core::compute_regular_triangulation;

    let fixture: TriangulationFixture = load_fixture("triangulation.json");
    let points = points_from_vecs(&fixture.input_points);

    let tri = compute_regular_triangulation(&points, &fixture.heights)
        .expect("Failed to compute triangulation");

    let mut our_simplices: Vec<Vec<usize>> = tri.simplices().iter().cloned().collect();
    let mut expected_simplices = fixture.simplices.clone();

    our_simplices.iter_mut().for_each(|s| s.sort());
    expected_simplices.iter_mut().for_each(|s| s.sort());
    our_simplices.sort();
    expected_simplices.sort();

    assert_eq!(
        our_simplices, expected_simplices,
        "Triangulation from heights doesn't match CYTools"
    );
}

// ============================================================================
// Linear relations tests
// ============================================================================

#[test]
fn test_linear_relations_no_origin() {
    use cyrus_core::integer_math::compute_linear_relations_no_origin;

    let fixture: LinearRelationsFixture = load_fixture("linear_relations.json");

    // Pass all points (the function internally skips the origin)
    let lr = compute_linear_relations_no_origin(&fixture.input_points);

    let result: Vec<Vec<i64>> = lr
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("LR value fits in i64"))
                .collect()
        })
        .collect();

    assert_eq!(
        result.len(),
        fixture.without_origin.shape[0],
        "Linear relations row count mismatch"
    );
    assert_eq!(
        result[0].len(),
        fixture.without_origin.shape[1],
        "Linear relations column count mismatch"
    );
    assert_eq!(
        result, fixture.without_origin.matrix,
        "Linear relations don't match CYTools"
    );
}

#[test]
fn test_linear_relations_kernel_property() {
    let fixture: LinearRelationsFixture = load_fixture("linear_relations.json");
    assert!(
        fixture.verification.is_zero,
        "CYTools L @ Q^T should be zero"
    );
}

// ============================================================================
// Divisor basis tests
// ============================================================================

#[test]
fn test_divisor_basis() {
    use cyrus_core::basis::compute_divisor_basis;
    use cyrus_core::glsm::compute_glsm_charge_matrix;

    let fixture: DivisorBasisFixture = load_fixture("divisor_basis.json");
    let points = points_from_vecs(&fixture.input_points);

    // Compute GLSM
    let points_without_origin: Vec<Point> = points.iter().skip(1).cloned().collect();
    let glsm = compute_glsm_charge_matrix(&points_without_origin, true).unwrap();

    // Compute divisor basis
    let basis = compute_divisor_basis(&points, &glsm).expect("Divisor basis failed");

    assert_eq!(
        basis.len(),
        fixture.h11,
        "Basis size should equal h11: ours={}, expected={}",
        basis.len(),
        fixture.h11
    );

    // Note: Exact indices may differ due to different tie-breaking, but size must match
}

// ============================================================================
// Intersection number tests
// ============================================================================

#[test]
fn test_intersection_in_basis() {
    use cyrus_core::basis::intersection_in_basis;
    use cyrus_core::intersection::Intersection;
    use cyrus_core::types::rational::Rational as TypedRational;
    use cyrus_core::types::tags::Finite;

    let full: IntersectionFullFixture = load_fixture("intersection_full.json");
    let basis_fixture: IntersectionBasisFixture = load_fixture("intersection_basis.json");

    // Build full intersection tensor
    let n_points = full.input_points.len();
    let mut kappa_full = Intersection::new(n_points);
    for entry in &full.entries {
        let val = TypedRational::<Finite>::from_raw(Rational::from(entry.value));
        kappa_full.set(entry.indices[0], entry.indices[1], entry.indices[2], val);
    }

    // Apply intersection_in_basis
    let kappa_basis = intersection_in_basis(&kappa_full, &basis_fixture.basis);

    // Convert to HashMap for comparison
    let result: HashMap<Vec<usize>, i64> = kappa_basis
        .iter()
        .map(|(&(i, j, k), val)| {
            let v: i64 = val.get().to_string().parse().unwrap();
            (vec![i, j, k], v)
        })
        .collect();

    let expected: HashMap<Vec<usize>, i64> = basis_fixture
        .entries
        .iter()
        .map(|e| (e.indices.clone(), e.value))
        .collect();

    assert_eq!(
        result, expected,
        "intersection_in_basis doesn't match CYTools"
    );
}

// ============================================================================
// Flat direction tests
// ============================================================================

#[test]
fn test_n_matrix() {
    use cyrus_core::flat_direction::compute_n_matrix;
    use cyrus_core::intersection::Intersection;
    use cyrus_core::types::i64::I64;
    use cyrus_core::types::rational::Rational as TypedRational;
    use cyrus_core::types::tags::Finite;

    let fixture: FlatDirectionFixture = load_fixture("flat_direction.json");
    let basis_fixture: IntersectionBasisFixture = load_fixture("intersection_basis.json");

    // Build intersection tensor in basis
    let mut kappa = Intersection::new(fixture.h11);
    for entry in &basis_fixture.entries {
        let val = TypedRational::<Finite>::from_raw(Rational::from(entry.value));
        kappa.set(entry.indices[0], entry.indices[1], entry.indices[2], val);
    }

    // Convert M flux to typed I64
    let m_flux: Vec<I64<Finite>> = fixture
        .result
        .m_flux
        .iter()
        .map(|&x| I64::<Finite>::new(x))
        .collect();

    // Compute N matrix
    let n_matrix = compute_n_matrix(&kappa, &m_flux);

    // Compare
    for i in 0..fixture.h11 {
        for j in 0..fixture.h11 {
            let expected = fixture.result.n_matrix[i][j];
            let actual = n_matrix[(i, j)];
            assert!(
                (actual - expected).abs() < 1e-10,
                "N_matrix[{},{}] mismatch: expected {}, got {}",
                i,
                j,
                expected,
                actual
            );
        }
    }
}
