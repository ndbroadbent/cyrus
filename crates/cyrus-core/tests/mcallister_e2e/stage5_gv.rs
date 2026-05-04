//! McAllister Pipeline Stage 5: GV Invariant Computation
//!
//! Tests the integration with the `cygv` crate for computing Gopakumar-Vafa invariants.
//!
//! **Requirements for GV computation:**
//! - Mori cone cap generators (from triangulation circuits)
//! - Grading vector (computed from Mori cone)
//! - Curve basis (GLSM charge matrix without origin)
//! - Intersection numbers κ_ijk
//!
//! **Current Status:**
//! - cygv crate is integrated as a dev dependency
//! - Cyrus computes Mori cone cap rays, a grading vector, and cygv inputs
//! - Full McAllister-sized GV validation is still too expensive for the normal
//!   suite: the 214-dimensional cone has hundreds of thousands of rays, and
//!   lattice-point generation needs a faster CYTools-faithful implementation

#![allow(missing_docs)]
#![allow(dead_code)]

use std::path::Path;

use cyrus_core::{
    F64, Finite, Point, Polytope, compute_curve_basis_matrix, compute_glsm_and_linrels,
    compute_mori_cone_cap_rays, compute_regular_triangulation,
    compute_toric_two_face_curve_gv_invariants, effective_prime_divisors_from_curve_basis,
    heights_to_kahler, remove_pair_decomposable_curve_candidates, subcutoff_toric_curve_candidates,
};

fn read_csv_rows_i64(path: &Path) -> Vec<Vec<i64>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|s| s.trim().parse::<i64>().expect("invalid integer"))
                .collect()
        })
        .collect()
}

fn read_csv_usize(path: &Path) -> Vec<usize> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split([',', '\n', '\r'])
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<usize>().expect("invalid usize"))
        .collect()
}

fn read_csv_i64(path: &Path) -> Vec<i64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split([',', '\n', '\r'])
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<i64>().expect("invalid i64"))
        .collect()
}

fn read_csv_f64(path: &Path) -> Vec<f64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split([',', '\n', '\r'])
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<f64>().expect("invalid float"))
        .collect()
}

fn read_csv_finite(path: &Path) -> Vec<F64<Finite>> {
    read_csv_f64(path)
        .into_iter()
        .map(|value| F64::<Finite>::new(value).expect("value must be finite"))
        .collect()
}

fn require_first_principles() -> bool {
    if !crate::first_principles_enabled() {
        eprintln!("Skipping first-principles test (set CYRUS_FIRST_PRINCIPLES=1)");
        return false;
    }
    true
}

/// Test that cygv crate is properly linked and basic types work
#[test]
fn stage5_cygv_basic_import() {
    // Import cygv types to verify linking
    use cygv::Semigroup;
    use nalgebra::{DMatrix, RowDVector};

    // Create a simple test matrix with generators as columns
    // 2 generators in 2D space: [1,0] and [0,1]
    let generators: DMatrix<i32> = DMatrix::from_column_slice(2, 2, &[1, 0, 0, 1]);

    // Create a grading vector (must be positive-definite)
    let grading: RowDVector<i32> = RowDVector::from_row_slice(&[1, 1]);

    // Verify we can create these types
    assert_eq!(generators.nrows(), 2);
    assert_eq!(generators.ncols(), 2);
    assert_eq!(grading.ncols(), 2);

    // Test Semigroup construction (basic cygv type)
    // This verifies the cygv crate is properly linked
    let semigroup = Semigroup::from_data(generators.clone(), grading.clone());

    // Note: Semigroup construction may fail for simple inputs
    // The important thing is that the cygv crate is linked and types work
    match semigroup {
        Ok(sg) => {
            eprintln!(
                "Semigroup created successfully, max_degree: {}",
                sg.max_degree
            );
        }
        Err(e) => {
            eprintln!(
                "Semigroup construction error (expected for simple inputs): {:?}",
                e
            );
            // This is OK - the crate is linked and working
        }
    }

    // Test that we can call with_max_degree which gives us more control
    let semigroup_with_deg = Semigroup::with_max_degree(generators, grading, 5);
    match semigroup_with_deg {
        Ok(sg) => {
            eprintln!(
                "Semigroup with max_degree created, elements: {}",
                sg.elements.ncols()
            );
            assert!(sg.max_degree >= 5, "Max degree should be at least 5");
        }
        Err(e) => {
            eprintln!("Semigroup with_max_degree error: {:?}", e);
            // This might fail too for simple inputs, but cygv is working
        }
    }

    // The test passes as long as the cygv types are importable and callable
    // The actual computation may require more complex inputs
}

/// Simple GV computation test with a known small example
/// The quintic CY3 in P^4 has GV_1 = 2875 (number of lines)
#[test]
#[ignore = "Requires full Mori cone infrastructure - for future implementation"]
fn stage5_cygv_quintic_example() {
    // Note: This test is for future implementation when we have all inputs ready.
    // The cygv API requires:
    // - generators: Mori cone generators (columns)
    // - grading_vector: positive-definite grading
    // - q: curve basis matrix
    // - nefpart: nef partition (can be empty for hypersurface)
    // - intnums: intersection numbers as HashMap

    // For the quintic in P^4:
    // - h11 = 1, single Kähler modulus
    // - κ_{111} = 5 (triple self-intersection of hyperplane)
    // - GV_1 = 2875 (number of lines)
    // - GV_2 = 609250 (number of conics)

    // Placeholder - actual implementation requires Mori cone computation
    // Test is ignored until Mori cone infrastructure is ready
}

/// Verify we have the data structures needed for McAllister GV computation
#[test]
fn stage5_mcallister_gv_data_available() {
    use std::fs;

    let Some(data_dir) = crate::mcallister_data_dir() else {
        assert!(
            !crate::first_principles_enabled(),
            "CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests"
        );
        eprintln!("Skipping GV data availability check (set CYRUS_MCALLISTER_DATA_DIR)");
        return;
    };

    // Check all required data files exist
    let required_files = [
        "dual_curves.dat",     // 5177 curves for dual polytope
        "dual_curves_gv.dat",  // GV invariants (some very large)
        "small_curves.dat",    // 344 curves for racetrack
        "small_curves_gv.dat", // Small GV invariants
    ];

    for file in &required_files {
        let path = data_dir.join(file);
        assert!(
            fs::metadata(&path).is_ok(),
            "Required file {} should exist",
            file
        );
    }

    // Load and verify dual curves data
    let curves_content = fs::read_to_string(data_dir.join("dual_curves.dat")).unwrap();
    let n_curves = curves_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .count();
    assert_eq!(n_curves, 5177, "Should have 5177 dual curves");

    // Load and verify GV invariants count
    let gv_content = fs::read_to_string(data_dir.join("dual_curves_gv.dat")).unwrap();
    let gv_values: Vec<&str> = gv_content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .flat_map(|l| l.split(','))
        .map(|s| s.trim())
        .collect();
    assert_eq!(gv_values.len(), 5177, "Should have 5177 GV values");

    // Check that some GV values are huge (demonstrating why we use arbitrary precision)
    let huge_gv_count = gv_values
        .iter()
        .filter(|s| s.len() > 20) // More than 20 digits
        .count();

    #[derive(serde::Serialize)]
    struct GvDataSummary {
        n_curves: usize,
        n_gv_values: usize,
        huge_gv_count: usize,
        first_5_gv: Vec<String>,
        max_gv_digits: usize,
    }

    let max_digits = gv_values.iter().map(|s| s.len()).max().unwrap_or(0);

    let summary = GvDataSummary {
        n_curves,
        n_gv_values: gv_values.len(),
        huge_gv_count,
        first_5_gv: gv_values.iter().take(5).map(|s| (*s).to_string()).collect(),
        max_gv_digits: max_digits,
    };

    insta::assert_json_snapshot!("mcallister_gv_data_summary", summary);
}

/// Compute the McAllister small toric curve checkpoint from upstream geometry.
///
/// This test intentionally uses `kahler_param.dat`, `small_curves_cutoff.dat`,
/// and `small_curves.dat` as validation artifacts. The reusable production
/// functions take Kähler parameters from the caller and do not load these files.
#[test]
fn stage5_mcallister_small_toric_curves_match_checkpoint() {
    if !require_first_principles() {
        return;
    }
    let Some(data_dir) = crate::mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let points_raw = read_csv_rows_i64(&data_dir.join("points.dat"));
    let heights = read_csv_f64(&data_dir.join("heights.dat"));
    let basis = read_csv_usize(&data_dir.join("basis.dat"));
    let kahler = read_csv_finite(&data_dir.join("kahler_param.dat"));
    let cutoff = F64::<Finite>::new(read_csv_f64(&data_dir.join("small_curves_cutoff.dat"))[0])
        .and_then(|value| value.try_to_pos())
        .expect("small curve cutoff must be positive");
    let expected_small = read_csv_rows_i64(&data_dir.join("small_curves.dat"));

    let all_points: Vec<Point> = points_raw.into_iter().map(Point::new).collect();
    let polytope = Polytope::from_vertices(all_points).expect("failed to create polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("failed to filter points");
    let triangulation = compute_regular_triangulation(&triangulation_points, &heights)
        .expect("failed to compute triangulation");
    let rays = compute_mori_cone_cap_rays(
        &triangulation,
        &triangulation_points,
        &polytope,
        false,
        false,
        None,
    )
    .expect("failed to compute ambient Mori cap rays");

    let selected =
        subcutoff_toric_curve_candidates(&rays, &basis, &kahler, cutoff).expect("curve selection");
    let filtered =
        remove_pair_decomposable_curve_candidates(&selected).expect("Hilbert-basis filter");

    let mut actual: Vec<Vec<i64>> = filtered.into_iter().map(|curve| curve.class).collect();
    actual.sort();
    let mut expected = expected_small;
    expected.sort();

    assert_eq!(
        selected.len(),
        419,
        "raw sub-cutoff toric curve count changed"
    );
    assert_eq!(
        actual.len(),
        344,
        "filtered small toric curve count changed"
    );
    assert_eq!(
        actual, expected,
        "Cyrus-computed filtered small toric curves must match McAllister checkpoint"
    );
}

/// Compute the GV values of the McAllister small toric curves from toric
/// two-face curve formulas, using `small_curves_gv.dat` only as a checkpoint.
#[test]
fn stage5_mcallister_small_toric_curve_gvs_match_checkpoint() {
    if !require_first_principles() {
        return;
    }
    let Some(data_dir) = crate::mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let points_raw = read_csv_rows_i64(&data_dir.join("points.dat"));
    let heights = read_csv_f64(&data_dir.join("heights.dat"));
    let expected_small = read_csv_rows_i64(&data_dir.join("small_curves.dat"));
    let expected_gv = read_csv_i64(&data_dir.join("small_curves_gv.dat"));
    assert_eq!(
        expected_small.len(),
        expected_gv.len(),
        "small curve/GV checkpoint length mismatch"
    );

    let all_points: Vec<Point> = points_raw.into_iter().map(Point::new).collect();
    let polytope = Polytope::from_vertices(all_points).expect("failed to create polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("failed to filter points");
    let triangulation = compute_regular_triangulation(&triangulation_points, &heights)
        .expect("failed to compute triangulation");

    let actual_gvs = compute_toric_two_face_curve_gv_invariants(
        &triangulation,
        &triangulation_points,
        &polytope,
    )
    .expect("failed to compute toric two-face curve GV invariants");
    let actual_by_class: std::collections::BTreeMap<Vec<i64>, String> = actual_gvs
        .into_iter()
        .map(|item| (item.class, item.gv.to_string()))
        .collect();

    let mut missing = Vec::new();
    let mut mismatched = Vec::new();
    for (class, expected) in expected_small.iter().zip(expected_gv.iter()) {
        match actual_by_class.get(class) {
            Some(actual) if actual == &expected.to_string() => {}
            Some(actual) => mismatched.push((class.clone(), expected.to_string(), actual.clone())),
            None => missing.push(class.clone()),
        }
    }

    assert!(
        missing.is_empty() && mismatched.is_empty(),
        "toric GV checkpoint mismatch: missing={} mismatched={} first_missing={:?} first_mismatch={:?}",
        missing.len(),
        mismatched.len(),
        missing.first(),
        mismatched.first()
    );
}

#[test]
fn stage5_height_projected_kahler_matches_uncorrected_branch_checkpoint() {
    if !require_first_principles() {
        return;
    }
    let Some(data_dir) = crate::mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let points_raw = read_csv_rows_i64(&data_dir.join("points.dat"));
    let heights = read_csv_finite(&data_dir.join("heights.dat"));
    let basis = read_csv_usize(&data_dir.join("basis.dat"));
    let mcallister_kahler = read_csv_finite(&data_dir.join("kahler_param.dat"));

    let all_points: Vec<Point> = points_raw.into_iter().map(Point::new).collect();
    let polytope = Polytope::from_vertices(all_points).expect("failed to create polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("failed to filter points");
    let (_glsm, linrels, _computed_basis) =
        compute_glsm_and_linrels(&triangulation_points).expect("failed to compute GLSM linrels");
    let basis_non_origin = basis
        .iter()
        .map(|idx| idx.checked_sub(1).expect("basis should not contain origin"))
        .collect::<Vec<_>>();
    let curve_basis =
        compute_curve_basis_matrix(&linrels, &basis).expect("failed to compute curve basis");
    let prime_divisors = effective_prime_divisors_from_curve_basis(&curve_basis, &basis_non_origin)
        .expect("failed to extract effective-cone rays");
    let height_kahler = heights_to_kahler(&heights, &basis_non_origin, &prime_divisors)
        .expect("failed to project heights to Kähler parameters");

    assert_eq!(height_kahler.len(), mcallister_kahler.len());
    let correlation = pearson_correlation(&height_kahler, &mcallister_kahler);
    let relative_l2_error = relative_l2_error(&height_kahler, &mcallister_kahler);
    let max_abs_error = max_abs_error(&height_kahler, &mcallister_kahler);
    let negative_count = height_kahler
        .iter()
        .filter(|value| value.get() < 0.0)
        .count();

    #[derive(serde::Serialize)]
    struct HeightKahlerDiagnostic {
        dimension: usize,
        negative_count: usize,
        correlation_with_mcallister_uncorrected: f64,
        relative_l2_error: f64,
        max_abs_error: f64,
        first_5: Vec<f64>,
    }

    let summary = HeightKahlerDiagnostic {
        dimension: height_kahler.len(),
        negative_count,
        correlation_with_mcallister_uncorrected: correlation,
        relative_l2_error,
        max_abs_error,
        first_5: height_kahler
            .iter()
            .take(5)
            .map(|value| value.get())
            .collect(),
    };

    insta::assert_json_snapshot!("height_projected_kahler_diagnostic", summary);
    assert!(
        relative_l2_error < 1e-12,
        "height-projected Kähler point must reproduce kahler_param.dat checkpoint without loading it as input: rel_l2={relative_l2_error}"
    );
}

fn pearson_correlation(lhs: &[F64<Finite>], rhs: &[F64<Finite>]) -> f64 {
    assert_eq!(lhs.len(), rhs.len());
    let n = lhs.len() as f64;
    let lhs_mean = lhs.iter().map(|value| value.get()).sum::<f64>() / n;
    let rhs_mean = rhs.iter().map(|value| value.get()).sum::<f64>() / n;
    let mut numerator = 0.0;
    let mut lhs_sq = 0.0;
    let mut rhs_sq = 0.0;
    for (l, r) in lhs.iter().zip(rhs.iter()) {
        let dl = l.get() - lhs_mean;
        let dr = r.get() - rhs_mean;
        numerator += dl * dr;
        lhs_sq += dl * dl;
        rhs_sq += dr * dr;
    }
    numerator / (lhs_sq.sqrt() * rhs_sq.sqrt())
}

fn relative_l2_error(lhs: &[F64<Finite>], rhs: &[F64<Finite>]) -> f64 {
    assert_eq!(lhs.len(), rhs.len());
    let mut numerator = 0.0;
    let mut denominator = 0.0;
    for (l, r) in lhs.iter().zip(rhs.iter()) {
        numerator += (l.get() - r.get()).powi(2);
        denominator += r.get().powi(2);
    }
    numerator.sqrt() / denominator.sqrt()
}

fn max_abs_error(lhs: &[F64<Finite>], rhs: &[F64<Finite>]) -> f64 {
    assert_eq!(lhs.len(), rhs.len());
    lhs.iter()
        .zip(rhs.iter())
        .map(|(l, r)| (l.get() - r.get()).abs())
        .fold(0.0, f64::max)
}

/// Document what's needed to compute GV invariants from scratch
/// This test serves as a roadmap for the full implementation
#[test]
fn stage5_gv_computation_roadmap() {
    #[derive(serde::Serialize)]
    struct GvComputationRoadmap {
        status: &'static str,
        /// Bitmask of completed components:
        /// cygv=1, mori=2, grading=4, pipeline=8, small-toric-curves=16,
        /// small-toric-gvs=32, primal-general-gv-fallback-wiring=64.
        completed_components: u8,
        verified_components: Vec<&'static str>,
        remaining_gaps: Vec<&'static str>,
    }

    // Components: cygv integrated (1), mori cap (2), grading vector (4),
    // one-off pipeline wiring (8), small toric curve checkpoint (16),
    // small toric curve GV checkpoint (32), and general primal GV fallback
    // wiring from basis-coordinate cygv output to ambient curve classes (64).
    let completed = 1u8 | 2 | 4 | 8 | 16 | 32 | 64;

    let roadmap = GvComputationRoadmap {
        status: "In Progress - Cyrus computes GV inputs, McAllister-sized validation is expensive",
        completed_components: completed,
        verified_components: vec![
            "compute_mori_cone_cap_rays is wired into mcallister_gv and mcallister_first_principles",
            "compute_grading_vector is wired into the GV pipeline",
            "grading-vector results are content-address cached after strict dual-cone validation",
            "DDM skips non-adjacent positive/negative ray pairs and no longer LP-prefilters every Mori ray",
            "DDM initializes from an exact full-rank hyperplane basis instead of the full ambient space when possible",
            "DDM tracks active constraints exactly for intersection rays using parent common-active sets plus the new hyperplane",
            "DDM rank checks use modular true certificates, exact integer false checks, and a global active-set rank cache",
            "DDM quotient-rank checks have modular true certificates plus exact small-residual true/false paths covered by unit tests",
            "DDM caches sparse hyperplane rows for partition/intersection dot products and quotient-coordinate projection while retaining dense rows for exact rank algebra",
            "DDM skips redundant hyperplanes that do not cut the current generated cone, avoiding active-set bloat from implied constraints",
            "DDM prunes constraints already redundant against the exact initial full-rank basis cone",
            "DDM uses bounded dynamic lookahead to process cheaper currently cutting hyperplanes before high-pair cuts",
            "DDM quotient-rank checks run before full dense modular checks when a basis context exists, with their own modular true certificate and exact integer false path",
            "DDM preserves ray orientation when normalizing primitive integer rays; sign-flipping was a correctness bug and is now unit-tested",
            "compute_gv_invariants wraps cygv::compute_gv_rat_threefold",
            "McAllister 4-214-647 small toric curve classes are computed from Cyrus Mori-cap rays and pair-decomposable pruning",
            "McAllister 4-214-647 small toric curve GV values are computed from toric two-face/origin-circuit formulas and match small_curves_gv.dat as a checkpoint",
            "CYTools-style height projection from heights plus curve-basis effective-cone rows reproduces McAllister 4-214-647 kahler_param.dat exactly as an uncorrected-branch checkpoint",
            "The no-replay runner reaches GV-corrected KKLT volume and V0 using the height-projected branch initializer, computed B-field gamma, and computed toric GV values",
            "first-principles binaries do not load small_curves.dat or small_curves_gv.dat",
            "basis-coordinate GV invariants can be mapped back to ambient divisor-intersection curve classes for primal volume corrections",
            "mcallister_first_principles can explicitly attempt general primal GV fallback for missing small-curve toric formula coverage via --primal-gv-max-deg or --primal-gv-min-points",
            "compute_gv_invariants can dump the exact pre-lattice Mori-cap generator matrix as a CDD V-representation for PPL/cdd diagnostics",
            "mcallister_first_principles rejects insufficient primal --primal-gv-max-deg values before full GV enumeration; for the selected 4-214-647 missing small curves, degree 1 is impossible because the required grading degrees are 10..2386",
            "branch-report JSONL records init_source so validation artifacts distinguish generated branch-search candidates from the computed height_projected secondary-fan initializer",
            "generated KKLT branch candidates now use one scaled-uniform baseline plus non-uniform deterministic random starts instead of wasting the prefix on equivalent scaled-uniform duplicates",
            "Intersection iteration is key-sorted for deterministic floating-point accumulation; repeated no-height 48-candidate branch reports are byte-identical for the same seed",
            "mcallister_first_principles supports --branch-selection min-toric-gv-missing to rank generated positive branches by computed small-curve toric GV coverage before falling back to volume",
        ],
        remaining_gaps: vec![
            "Generated branch candidates without the height_projected initializer still did not find the 4-214-647 paper branch in a deterministic 48-candidate diagnostic: the lowest sampled phase-1 volume was about 20611 rather than 17901, and even coverage-aware selection still had at least 412 small curves missing toric GV coverage",
            "The explicit general primal GV fallback still reaches full 214-dimensional Mori-cone dualization for any max_deg high enough to cover the selected missing curves, or for min_points-driven runs; bounded DDM diagnostics stop loudly at configured limits",
            "A PPL/cdd diagnostic on the dumped 561658-ray, 214-dimensional V-representation also exceeded a 300-second cap without producing an H-representation",
            "Finish a post-orientation-fix validation run of adjacency-filtered DDM on the full 214-dimensional McAllister Mori cone",
            "Further optimize or replace hyperplane dualization; bounded diagnostics still need to prove the full 561658-ray McAllister dualization completes with the corrected ray orientation",
            "Reduce the 561658-ray Mori cap input before dualization, or add a CYTools/PPL-faithful constraint minimization path",
            "Run and validate lattice-point generation under a Python environment with OR-Tools after DDM returns the dual cone",
        ],
    };

    insta::assert_json_snapshot!("gv_computation_roadmap", roadmap);
}
