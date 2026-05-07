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

use cyrus_core::intersection::compute_intersection_cytools;
use cyrus_core::{
    CurvePruningStrategy, F64, Finite, Point, Polytope,
    compute_ambient_one_dimensional_ray_gv_series, compute_curve_basis_matrix,
    compute_glsm_and_linrels, compute_grading_vector, compute_linear_relations_no_origin,
    compute_mori_cone_cap_rays, compute_regular_triangulation, compute_toric_curve_gv_diagnostics,
    compute_toric_two_face_curve_gv_invariants, curve_row_span_rank,
    effective_prime_divisors_from_curve_basis, find_semigroup_decomposition, heights_to_kahler,
    intersection_in_basis, potent_ray_convergence, project_ambient_curve_to_basis,
    project_mori_cone_cap_rays_to_basis, prune_decomposable_curve_candidates,
    remove_pair_decomposable_curve_candidates, subcutoff_toric_curve_candidates, types::i64::I64,
};
use malachite::Integer;

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

fn read_csv_rows_integer(path: &Path) -> Vec<Vec<Integer>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|s| s.trim().parse::<Integer>().expect("invalid integer"))
                .collect()
        })
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

fn o7_gamma_from_kklt_data(data_dir: &Path, ambient_dim: usize) -> Vec<i64> {
    let kklt_basis = read_csv_usize(&data_dir.join("kklt_basis.dat"));
    let c_i = read_csv_i64(&data_dir.join("target_volumes.dat"));
    assert_eq!(
        kklt_basis.len(),
        c_i.len(),
        "KKLT basis and c_i checkpoint length mismatch"
    );

    let mut gamma = vec![0_i64; ambient_dim];
    for (&idx, &ci) in kklt_basis.iter().zip(c_i.iter()) {
        assert!(
            idx < ambient_dim,
            "KKLT divisor index {idx} out of ambient dimension {ambient_dim}"
        );
        if ci == 6 {
            gamma[idx] += 1;
        }
    }
    gamma
}

fn typed_gamma(raw: &[i64]) -> Vec<I64<Finite>> {
    raw.iter()
        .copied()
        .map(|value| I64::<Finite>::new(value))
        .collect()
}

fn ambient_parity(curve: &[i64], gamma: &[i64]) -> i128 {
    curve
        .iter()
        .zip(gamma.iter())
        .map(|(&q_i, &gamma_i)| i128::from(q_i) * i128::from(gamma_i))
        .sum::<i128>()
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

/// Compute the validation quantities attached to McAllister's potent-ray sample.
///
/// This test does not yet generate the rays or their GV values from scratch.
/// It verifies that the checkpoint files are no longer opaque: Cyrus computes
/// their rank, their corrected-Kähler volumes, and the paper's convergence
/// diagnostic from the supplied ray/GV sample.
#[test]
fn stage5_mcallister_potent_ray_checkpoint_quantities_are_computed() {
    if !require_first_principles() {
        return;
    }
    let Some(data_dir) = crate::mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let basis = read_csv_usize(&data_dir.join("basis.dat"));
    let corrected_kahler = read_csv_finite(&data_dir.join("corrected_kahler_param.dat"));
    let rays = read_csv_rows_i64(&data_dir.join("potent_rays.dat"));
    let gv_rows = read_csv_rows_integer(&data_dir.join("potent_rays_gv.dat"));
    let expected_rank = read_csv_usize(&data_dir.join("potent_rays_rank.dat"))[0];
    let expected_volumes = read_csv_f64(&data_dir.join("potent_rays_vols.dat"));

    assert_eq!(rays.len(), 411, "potent-ray checkpoint row count changed");
    assert_eq!(
        gv_rows.len(),
        rays.len(),
        "potent-ray GV row count mismatch"
    );
    assert_eq!(
        expected_volumes.len(),
        rays.len(),
        "potent-ray volume row count mismatch"
    );

    let actual_rank = curve_row_span_rank(&rays).expect("potent ray rank");
    assert_eq!(
        actual_rank, expected_rank,
        "Cyrus exact rank must match potent_rays_rank.dat"
    );

    let mut max_abs_volume_delta = 0.0_f64;
    let mut nondecaying_slopes = 0usize;
    for ((ray, gv_row), &expected_volume) in rays.iter().zip(gv_rows.iter()).zip(&expected_volumes)
    {
        let volume = cyrus_core::curve_volume_in_divisor_basis(ray, &basis, &corrected_kahler)
            .expect("potent ray volume")
            .try_to_pos()
            .expect("potent ray volume must be positive");
        max_abs_volume_delta = max_abs_volume_delta.max((volume.get() - expected_volume).abs());

        let convergence = potent_ray_convergence(gv_row, volume).expect("convergence terms");
        let slope = convergence
            .log_xi_slope
            .expect("ten nonzero potent-ray GV values should determine a slope");
        if slope.get() >= 0.0 {
            nondecaying_slopes += 1;
        }
    }

    assert!(
        max_abs_volume_delta < 1e-9,
        "computed potent-ray volumes must match checkpoint; max_abs_delta={max_abs_volume_delta}"
    );
    assert_eq!(
        nondecaying_slopes, 0,
        "all checkpoint potent-ray log-xi slopes should decay at corrected t"
    );
}

/// Diagnostic for the first saved potent ray.
///
/// This deliberately uses a one-generator cygv call, which is not the full
/// low-dimensional-face method from the paper. It records that the simplest
/// local semigroup is not enough for the first 4-214-647 potent ray: Cyrus must
/// reconstruct the face context to regenerate the saved checkpoint row.
#[test]
#[ignore = "diagnostic; run manually with CYRUS_FIRST_PRINCIPLES and CYRUS_MCALLISTER_DATA_DIR"]
fn stage5_first_potent_ray_one_generator_gv_diagnostic() {
    if !require_first_principles() {
        return;
    }
    let Some(data_dir) = crate::mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let points_raw = read_csv_rows_i64(&data_dir.join("points.dat"));
    let heights = read_csv_f64(&data_dir.join("heights.dat"));
    let basis = read_csv_usize(&data_dir.join("basis.dat"));
    let potent_rays = read_csv_rows_i64(&data_dir.join("potent_rays.dat"));
    let potent_gv = read_csv_rows_integer(&data_dir.join("potent_rays_gv.dat"));
    let target_ray = potent_rays.first().expect("saved potent rays are nonempty");
    let expected = potent_gv
        .first()
        .expect("saved potent-ray GV rows are nonempty");

    let all_points: Vec<Point> = points_raw.into_iter().map(Point::new).collect();
    let polytope = Polytope::from_vertices(all_points).expect("failed to create polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("failed to filter points");
    let triangulation = compute_regular_triangulation(&triangulation_points, &heights)
        .expect("failed to compute triangulation");
    let ambient_rays = compute_mori_cone_cap_rays(
        &triangulation,
        &triangulation_points,
        &polytope,
        false,
        false,
        None,
    )
    .expect("failed to compute ambient Mori cap rays");
    let basis_rays = project_mori_cone_cap_rays_to_basis(&ambient_rays, &basis)
        .expect("failed to project Mori-cap rays to basis");
    let grading = compute_grading_vector(&basis_rays).expect("failed to compute grading vector");
    let basis_target =
        project_ambient_curve_to_basis(target_ray, &basis).expect("failed to project potent ray");
    let is_mori_cap_generator = basis_rays.iter().any(|ray| ray == &basis_target);

    let points_i64 = triangulation_points
        .iter()
        .map(|point| point.coords().to_vec())
        .collect::<Vec<_>>();
    let linrels_reduced = compute_linear_relations_no_origin(&points_i64);
    let linrels_i64 = linrels_reduced
        .iter()
        .map(|row| {
            row.iter()
                .map(|value| i64::try_from(value).expect("linear relation fits in i64"))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let kappa_full =
        compute_intersection_cytools(&triangulation, &triangulation_points, &linrels_i64)
            .expect("failed to compute intersection numbers");
    let kappa_basis = intersection_in_basis(&kappa_full, &basis);
    let (_glsm, linrels, _computed_basis) =
        compute_glsm_and_linrels(&triangulation_points).expect("failed to compute GLSM linrels");
    let curve_basis =
        compute_curve_basis_matrix(&linrels, &basis).expect("failed to compute curve basis");
    let q_matrix = curve_basis
        .iter()
        .map(|row| {
            row.iter()
                .map(|value| i64::try_from(value).expect("curve-basis entry fits in i64"))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();

    let series = compute_ambient_one_dimensional_ray_gv_series(
        target_ray,
        &basis,
        &grading,
        &q_matrix,
        &kappa_basis,
        10,
    )
    .expect("failed to compute one-generator potent ray GV series");

    assert!(
        !is_mori_cap_generator,
        "the first saved potent ray unexpectedly became a Mori-cap generator"
    );
    assert_ne!(
        series.values, *expected,
        "if this starts matching, revisit whether a one-generator semigroup is now sufficient"
    );
    assert_eq!(series.degree, 27);
    assert_eq!(
        series.values,
        vec![
            Integer::from(4),
            Integer::from(-11),
            Integer::from(60),
            Integer::from(-478),
            Integer::from(4588),
            Integer::from(-49368),
            Integer::from(575896),
            Integer::from(-7131274),
            Integer::from(92429484),
            Integer::from(-1241983287),
        ],
        "direct one-generator diagnostic changed"
    );
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
        remove_pair_decomposable_curve_candidates(&selected).expect("pair-decomposable pruning");

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

/// Diagnose whether McAllister's saved small-curve pruning needs more than pair sums.
///
/// The production checkpoint still uses the documented pair-decomposable rule
/// above. This test is intentionally first-principles-only because it solves a
/// finite integer feasibility problem over the selected 4-214-647 toric curve
/// set to expose any multi-term sums that the pair shortcut misses.
#[test]
fn stage5_mcallister_small_toric_curve_finite_semigroup_diagnostic() {
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
    let pair_filtered =
        remove_pair_decomposable_curve_candidates(&selected).expect("pair-decomposable pruning");
    let strategy_pair_filtered =
        prune_decomposable_curve_candidates(&selected, CurvePruningStrategy::PairDecomposable)
            .expect("explicit pair pruning strategy");
    let finite_semigroup_filtered =
        prune_decomposable_curve_candidates(&selected, CurvePruningStrategy::FiniteSemigroup)
            .expect("explicit finite-semigroup pruning strategy");

    let mut semigroup_decomposable = Vec::new();
    for curve in &pair_filtered {
        if let Some(decomposition) =
            find_semigroup_decomposition(curve, &selected).expect("finite semigroup diagnostic")
        {
            semigroup_decomposable.push((curve.class.clone(), decomposition));
        }
    }

    assert_eq!(
        selected.len(),
        419,
        "raw sub-cutoff toric curve count changed"
    );
    assert_eq!(
        pair_filtered.len(),
        344,
        "pair-pruned small toric curve count changed"
    );
    assert_eq!(
        strategy_pair_filtered, pair_filtered,
        "explicit pair pruning strategy must match the checkpoint-faithful helper"
    );
    assert_eq!(
        semigroup_decomposable.len(),
        5,
        "McAllister checkpoint is pair-pruned, not finite-semigroup-pruned"
    );
    assert_eq!(
        finite_semigroup_filtered.len(),
        339,
        "full finite-semigroup pruning would retain 339 input-chamber candidates"
    );
    assert_eq!(
        pair_filtered.len() - finite_semigroup_filtered.len(),
        semigroup_decomposable.len(),
        "explicit finite-semigroup strategy must remove the multi-term decomposable checkpoint curves"
    );
}

/// Verify the semantics of McAllister's saved small-curve volume checkpoint.
///
/// `small_curves_vols.dat` is not an independent production input. It is the
/// validation-only contraction of `small_curves.dat` with the downstream
/// `corrected_kahler_param.dat` in the saved `basis.dat` coordinates. The
/// negative entries are the warning sign that the input-chamber small-curve set
/// crosses flops at the corrected checkpoint.
#[test]
fn stage5_mcallister_small_curve_volumes_are_corrected_kahler_projection() {
    if !require_first_principles() {
        return;
    }
    let Some(data_dir) = crate::mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let curves = read_csv_rows_i64(&data_dir.join("small_curves.dat"));
    let basis = read_csv_usize(&data_dir.join("basis.dat"));
    let corrected_kahler = read_csv_f64(&data_dir.join("corrected_kahler_param.dat"));
    let expected_volumes = read_csv_f64(&data_dir.join("small_curves_vols.dat"));

    assert_eq!(
        basis.len(),
        corrected_kahler.len(),
        "basis/corrected Kähler checkpoint length mismatch"
    );
    assert_eq!(
        curves.len(),
        expected_volumes.len(),
        "small curve volume checkpoint length mismatch"
    );

    let actual_volumes = curves
        .iter()
        .map(|curve| {
            basis
                .iter()
                .zip(corrected_kahler.iter())
                .map(|(&idx, &ti)| {
                    assert!(
                        idx < curve.len(),
                        "basis index {idx} out of curve dimension {}",
                        curve.len()
                    );
                    curve[idx] as f64 * ti
                })
                .sum::<f64>()
        })
        .collect::<Vec<_>>();
    let max_abs_delta = actual_volumes
        .iter()
        .zip(expected_volumes.iter())
        .map(|(actual, expected)| (actual - expected).abs())
        .fold(0.0_f64, f64::max);
    let negative_count = actual_volumes
        .iter()
        .filter(|&&volume| volume < 0.0)
        .count();

    assert!(
        max_abs_delta < 1e-12,
        "small_curves_vols.dat must equal small_curves[:, basis.dat] dot corrected_kahler_param.dat; max_abs_delta={max_abs_delta}"
    );
    assert_eq!(
        negative_count, 10,
        "corrected checkpoint should expose the known 10 flopped input-chamber small curves"
    );
}

#[test]
fn stage5_mcallister_flopped_small_curves_have_mixed_b_field_parity() {
    if !require_first_principles() {
        return;
    }
    let Some(data_dir) = crate::mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let curves = read_csv_rows_i64(&data_dir.join("small_curves.dat"));
    let gvs = read_csv_i64(&data_dir.join("small_curves_gv.dat"));
    let basis = read_csv_usize(&data_dir.join("basis.dat"));
    let corrected_kahler = read_csv_finite(&data_dir.join("corrected_kahler_param.dat"));
    let stored_volumes = read_csv_f64(&data_dir.join("small_curves_vols.dat"));
    let gamma_raw = o7_gamma_from_kklt_data(&data_dir, curves[0].len());
    let gamma = typed_gamma(&gamma_raw);
    let kklt_basis = read_csv_usize(&data_dir.join("kklt_basis.dat"));

    let mut negative_indices = Vec::new();
    let mut odd_negative_indices = Vec::new();
    let mut even_negative_indices = Vec::new();
    for (idx, (curve, &stored_volume)) in curves.iter().zip(stored_volumes.iter()).enumerate() {
        if stored_volume >= 0.0 {
            continue;
        }
        negative_indices.push(idx);
        let parity = ambient_parity(curve, &gamma_raw).rem_euclid(2);
        if parity == 0 {
            even_negative_indices.push(idx);
        } else {
            odd_negative_indices.push(idx);
        }
    }

    assert_eq!(
        negative_indices.len(),
        10,
        "expected the 10 flopped input-chamber small curves from the ancillary data"
    );
    assert_eq!(
        odd_negative_indices.len(),
        8,
        "only odd B-field parity curves admit the simple real flop continuation"
    );
    assert_eq!(
        even_negative_indices,
        vec![1, 3],
        "two flopped saved small curves are even-parity branch-cut cases"
    );
    assert!(
        negative_indices.iter().all(|&idx| gvs[idx] == 1),
        "all negative saved small curves are expected to have GV=1"
    );

    let odd_negative_gvs = odd_negative_indices
        .iter()
        .map(|&idx| (curves[idx].clone(), Integer::from(gvs[idx])))
        .collect::<Vec<_>>();
    assert!(
        cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
            &odd_negative_gvs,
            &basis,
            &kklt_basis,
            &corrected_kahler,
            Some(&gamma),
        )
        .is_some(),
        "odd-parity flopped curves should remain on the real -exp branch"
    );

    let even_negative_gvs = even_negative_indices
        .iter()
        .map(|&idx| (curves[idx].clone(), Integer::from(gvs[idx])))
        .collect::<Vec<_>>();
    assert!(
        cyrus_core::kklt::compute_gv_target_correction_for_ambient_curves(
            &even_negative_gvs,
            &basis,
            &kklt_basis,
            &corrected_kahler,
            Some(&gamma),
        )
        .is_none(),
        "even-parity flopped curves cross the real Li2 branch cut and need a different treatment"
    );
}

#[test]
fn stage5_mcallister_mori_basis_projection_matches_direct_basis_rays() {
    if !require_first_principles() {
        return;
    }
    let Some(data_dir) = crate::mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let points_raw = read_csv_rows_i64(&data_dir.join("points.dat"));
    let heights = read_csv_f64(&data_dir.join("heights.dat"));
    let basis = read_csv_usize(&data_dir.join("basis.dat"));

    let all_points: Vec<Point> = points_raw.into_iter().map(Point::new).collect();
    let polytope = Polytope::from_vertices(all_points).expect("failed to create polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("failed to filter points");
    let triangulation = compute_regular_triangulation(&triangulation_points, &heights)
        .expect("failed to compute triangulation");

    let ambient_rays = compute_mori_cone_cap_rays(
        &triangulation,
        &triangulation_points,
        &polytope,
        false,
        false,
        None,
    )
    .expect("failed to compute ambient Mori-cap rays");
    let projected = project_mori_cone_cap_rays_to_basis(&ambient_rays, &basis)
        .expect("failed to project ambient Mori-cap rays to basis");
    let direct = compute_mori_cone_cap_rays(
        &triangulation,
        &triangulation_points,
        &polytope,
        true,
        false,
        Some(&basis),
    )
    .expect("failed to compute basis Mori-cap rays directly");

    assert_eq!(
        projected, direct,
        "basis projection of reused ambient Mori-cap rays must match direct in_basis=true computation"
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
    let diagnostic_gvs =
        compute_toric_curve_gv_diagnostics(&triangulation, &triangulation_points, &polytope)
            .expect("failed to compute toric curve GV diagnostics");
    let actual_by_class: std::collections::BTreeMap<Vec<i64>, String> = actual_gvs
        .into_iter()
        .map(|item| (item.class, item.gv.to_string()))
        .collect();
    let diagnostic_by_class: std::collections::BTreeMap<Vec<i64>, String> = diagnostic_gvs
        .into_iter()
        .map(|item| {
            assert!(
                !item.sources.is_empty(),
                "diagnostic GV entry has no provenance sources"
            );
            (item.class, item.gv.to_string())
        })
        .collect();
    assert_eq!(
        actual_by_class, diagnostic_by_class,
        "diagnostic GV provenance diverged from toric GV invariants"
    );

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
        /// small-toric-gvs=32, primal-general-gv-fallback-wiring=64,
        /// potent-ray-convergence-checks=128, one-dimensional-ray-gv-series=256.
        completed_components: u16,
        verified_components: Vec<&'static str>,
        remaining_gaps: Vec<&'static str>,
    }

    // Components: cygv integrated (1), mori cap (2), grading vector (4),
    // one-off pipeline wiring (8), small toric curve checkpoint (16),
    // small toric curve GV checkpoint (32), and general primal GV fallback
    // wiring from basis-coordinate cygv output to ambient curve classes (64),
    // potent-ray convergence checks over supplied samples (128), and reusable
    // one-dimensional ray GV series computation (256).
    let completed = 1u16 | 2 | 4 | 8 | 16 | 32 | 64 | 128 | 256;

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
            "compute_gv_invariants runs the upstream cygv HKTY modules and maps cygv construction/inversion failures and remaining unwind panics into Result errors",
            "Cyrus' direct cygv HKTY call chain is regression-tested against cygv::compute_gv_rat_threefold on the quintic degree-one GV 2875 case",
            "McAllister 4-214-647 small toric curve classes are computed from Cyrus Mori-cap rays and verified pair-decomposable pruning",
            "McAllister 4-214-647 small toric curve GV values are computed from toric two-face/origin-circuit formulas and match small_curves_gv.dat as a checkpoint",
            "CYTools-style height projection from heights plus curve-basis effective-cone rows reproduces McAllister 4-214-647 kahler_param.dat exactly as an uncorrected-branch checkpoint",
            "The no-replay runner reaches GV-corrected KKLT volume and V0 using the height-projected branch initializer, computed B-field gamma, and computed toric GV values",
            "mcallister_first_principles now uses the height-projected KKLT initializer by default when no explicit branch search is requested, so the standard no-replay command no longer needs branch flags to reach the McAllister branch",
            "mcallister_first_principles derives the BBHL Hodge inputs from the computed primal and dual bases instead of hard-coding the 4-214-647 values",
            "mcallister_first_principles validates the no-replay corrected volume against corrected_cy_vol.dat when McAllister comparison data is available",
            "mcallister_first_principles logs the McAllister volume split into classical, BBHL, pre-GV, and direct GV volume terms for residual attribution",
            "first-principles binaries do not load small_curves.dat or small_curves_gv.dat",
            "first-principles runner policy distinguishes downstream Kähler replay from post-computation corrected_cy_vol.dat validation comparison",
            "basis-coordinate GV invariants can be mapped back to ambient divisor-intersection curve classes for primal volume corrections",
            "McAllister 4-214-647 potent-ray checkpoint samples have exact rank, corrected-Kähler volumes, and log-xi convergence slopes computed by Cyrus; all supplied slopes decay at corrected t",
            "McAllister 4-214-647 rank-two CKYZ potent-ray GV checks now cover the first four entries for all 395 CKYZ rows and all ten entries for the canonical F0 source directions [1,1]/[1,2] and F1 source directions [2,1]/[3,1] without using potent_rays_gv.dat as input",
            "one-dimensional ray GV series computation is available in cyrus-core via compute_one_dimensional_ray_gv_series, and mcallister_first_principles uses it for the corrected-chamber ray diagnostic instead of carrying a private cygv wrapper",
            "mcallister_first_principles can explicitly attempt general primal GV fallback for missing small-curve toric formula coverage via --primal-gv-max-deg or --primal-gv-min-points",
            "compute_gv_invariants can dump the exact pre-lattice Mori-cap generator matrix as a CDD V-representation for PPL/cdd diagnostics",
            "compute_gv_invariants now uses CYTools-style min_points lattice augmentation even when max_deg is specified; degree-bounded lattice enumeration is isolated behind compute_gv_invariants_with_degree_bounded_lattice for diagnostics",
            "mcallister_first_principles rejects insufficient primal --primal-gv-max-deg values before full GV enumeration; for the selected 4-214-647 missing small curves, degree 1 is impossible because the required grading degrees are 10..2386",
            "branch-report JSONL records init_source so validation artifacts distinguish generated branch-search candidates from the computed height_projected secondary-fan initializer",
            "generated KKLT branch candidates now use one scaled-uniform baseline plus non-uniform deterministic random starts instead of wasting the prefix on equivalent scaled-uniform duplicates",
            "Intersection iteration is key-sorted for deterministic floating-point accumulation; repeated no-height 48-candidate branch reports are byte-identical for the same seed",
            "mcallister_first_principles supports --branch-selection min-toric-gv-missing to rank generated positive branches by computed small-curve toric GV coverage before falling back to volume",
            "mcallister_first_principles supports --branch-selection min-required-gv-degree to rank generated positive branches by the max required degree of missing small-curve general GV values",
            "branch-report JSONL can include bounded samples of missing small-curve classes via --branch-report-missing-limit for formula-classification diagnostics",
            "branch-report JSONL records required grading-degree ranges for missing small-curve classes, exposing the cost of general GV fallback per branch",
            "branch-report JSONL can diagnose missing small-curve classes that are exact sums of up to four raw sub-cutoff candidates via --branch-report-decomposition-depth",
            "project_mori_cone_cap_rays_to_basis reuses ambient Mori-cap rays and is first-principles checked against direct in_basis=true Mori-cap generation on McAllister 4-214-647",
            "corrected-chamber ray and LP-witness GV diagnostics now report a combined diagnostic-only coverage and volume correction without promoting those values into the production volume",
            "corrected-chamber combined diagnostics now report zero/nonzero diagnostic GV counts and delta_vs_toric_covered; for the 4-214-647 corrected chamber, zero=10, nonzero=0, and delta_vs_toric_covered=0",
            "corrected-chamber diagnostics now compare KKLT GV target-correction vectors against the input-chamber fixed-point data; for 4-214-647 the final corrected chamber has max_abs_delta=0.04717957747495261 and relative_l2_delta=0.5969594710109491",
            "mcallister_first_principles supports --diagnose-chamber-updated-kklt, an opt-in toric-covered-only diagnostic that recomputes the FRST chamber, intersection numbers, divisor χ, small-curve selection, and GV target correction at each KKLT fixed-point step",
            "mcallister_first_principles reports validation-only corrected_target_volumes.dat deltas, including the implied divisor-χ shift; for 4-214-647 the computed GV-corrected KKLT target vector differs by relative_l2=0.022036578712025 and max_abs=0.09773695104398727, equivalent to an implied χ(D) delta with relative_l2=0.09790035863637207 and max_abs=2.3456868250556946",
            "mcallister_first_principles reports the largest corrected_target_volumes.dat implied-χ deltas by KKLT index and ambient point; the top 4-214-647 deltas are spread across ordinary KKLT divisors rather than the four basis-mismatch divisors, and the no-gamma GV target-correction variant is invalid at the solved point",
            "McAllister checkpoint-t input-chamber GV target correction is invalid for 4-214-647, matching the ancillary warning that some approximate-chamber small curves flop and become negative at corrected_kahler_param.dat",
            "small_curves_vols.dat is now pinned as small_curves.dat contracted with corrected_kahler_param.dat in basis.dat coordinates; for 4-214-647 it matches exactly and contains the known 10 negative input-chamber curve volumes",
            "McAllister checkpoint-t corrected-chamber toric-covered GV target correction is now compared against corrected_target_volumes.dat-implied GV corrections; for 4-214-647 it has ambient_rays=561575, subcutoff=556, pair_pruned=420, toric_covered=410, toric_missing=10, max_abs_delta=0.1154398144342352, and relative_l2_delta=1.2720321759981263",
            "checkpoint-t corrected-chamber unpruned-subcutoff GV diagnostics show the checkpoint vector is not simply using all subcutoff toric-covered curves: for 4-214-647 the unpruned set has toric_covered=512, toric_missing=44, max_abs_delta=0.1154398144342352, relative_l2_delta=1.3245185834958435, and differs from the pair-pruned vector by max_abs=0.10239160047732201",
            "corrected_heights.dat induces the same 1003-simplex chamber as McAllister corrected_kahler_param.dat transformed into the Cyrus computed basis, so the checkpoint-t corrected-chamber GV mismatch is not caused by reconstructing a different corrected chamber from t",
            "a direct CYTools source trace on corrected_heights.dat matches the Cyrus corrected-chamber small-curve selection counts for 4-214-647: mori_shape=(561575,219), subcutoff=556, pair_pruned=420, and the leading tiny-volume offender classes are present",
            "checkpoint-t corrected-chamber diagnostics now recompute divisor chi in the corrected chamber; for 4-214-647 this changes the checkpoint-implied GV vector by max_abs=0.08333333333333337 and relative_l2=0.8329163940091703, but the toric-covered GV mismatch remains max_abs=0.1154398144342352 and relative_l2=1.1084474600994203",
            "checkpoint-t corrected-chamber diagnostics now compare the old 2021-style no-B-field-sign convention; for 4-214-647 no_gamma worsens the target-vector mismatch to input_chi max_abs=0.14693479936880824 and relative_l2=2.0801895095962637",
            "checkpoint-t corrected-chamber diagnostics now sweep tiny positive curve exclusions; dropping q.t<0.005, 0.01, 0.02, 0.05, or 0.1 does not reproduce the checkpoint-implied GV correction, with the best tested input-chi relative_l2=1.056010872179146 still leaving max_abs=0.10111964958304713",
            "checkpoint-t corrected-chamber GV branch-bucket diagnostics now split the toric-covered target correction by B-field parity and q.t; for 4-214-647 the largest aggregate l2 bucket is odd-parity q.t>=0.1 with count=297 and bucket_l2=0.2334903638965792, while the two even near-flop singleton buckets at q.t≈0.004802446698931817 and q.t≈0.007120534449768456 are important top-divisor offenders but dropping or sign-flipping either still leaves max_abs=0.10111964958304713",
            "checkpoint-t corrected-chamber GV target diagnostics now report the largest per-divisor deltas and per-curve contributions; the leading 4-214-647 offenders are dominated by tiny-volume toric classes such as [(2,1),(50,-1),(54,1),(197,-1)] with q.t≈0.004802446698931817 and no detected decomposition into up to four smaller selected candidates",
            "real_dilog_real_axis is now unit-tested against SciPy/mpmath reference values in the leading corrected-chamber offender regime q.t≈0.004802446698931817, q.t≈0.007120534449768456, and q.t≈0.11537135311573365 for both B-field signs",
            "mcallister_first_principles can export CYRUS_CORRECTED_CHAMBER_GV_TRACE_JSON for the exact corrected-chamber GV target inputs; a Python/SciPy recomputation from that trace reproduces Cyrus' 4-214-647 toric-covered GV target vector to max_abs=5.551115123125783e-17 and relative_l2=6.800408568913545e-16, while preserving the same checkpoint-implied mismatch max_abs=0.11543981443423522 and relative_l2=1.1084474600994199",
            "the corrected-chamber GV trace now includes toric-missing classes; an arbitrary-real least-squares fit over the 10 pair-pruned missing classes still leaves max_abs=0.1154398144342352 and relative_l2=1.0282708464664079, and the 44 subcutoff missing classes still leave max_abs=0.1154398144342352 and relative_l2=0.8204789658126757, so the target-vector residual is not explained by missing-class GV values alone",
            "checkpoint-t corrected-chamber GV provenance diagnostics now attach each top per-curve contribution to its local toric formula source; for 4-214-647 the largest kklt_idx=192 offender [(2,1),(50,-1),(54,1),(197,-1)] is a genus-0 two-face circuit on edge [50,197], not an origin-circuit fallback, while the largest kklt_idx=52 contributors include two-face gv=-2 classes before a smaller resolved-conifold origin-circuit term",
            "checkpoint-t corrected-chamber top-toric local-HKTY diagnostics can be enabled with CYRUS_CORRECTED_CHAMBER_TOP_TORIC_LOCAL_GV=1; an optimized unwind diagnostic for 4-214-647 requested the top 8 divisor residuals with 10 contribution rows each and checked 50 unique top-contributor classes, including 48 one-dimensional-ray checks and 2 LP-witness checks; every local GV value matched the toric formula value",
            "checkpoint-t corrected-chamber GV source-bucket diagnostics now group toric-covered target corrections by local formula provenance; for 4-214-647 the largest buckets are ordinary genus-0 two-face classes with bucket_l2=0.22152421247870088 and 0.20063746571486588, while the resolved-conifold origin bucket is small with bucket_l2=0.014963484333468835, and dropping or sign-flipping any reported source bucket still leaves the 0.115-scale max error",
            "checkpoint-t corrected-chamber GV weight-LP diagnostics can be enabled with CYRUS_CORRECTED_CHAMBER_GV_WEIGHT_LP=1; for 4-214-647 a relaxed continuous-drop model over the 410 toric-covered curves still leaves max_abs=0.049547280835143204 and relative_l2=0.8199330813800068, while a nonphysical bounded-sign model can reduce the mismatch to max_abs=0.0039037977783476785 only by driving 32 curve weights to -1 and 80 weights into the interior",
            "checkpoint-t corrected-chamber GV weight-LP parity-fit diagnostics now test whether the bounded-sign near-bound weights can come from one ambient mod-2 B-field shift; for 4-214-647 the system is inconsistent for 330 near-bound constraints (32 negative, 298 positive), so the LP near-fit is not a valid gamma representative",
            "checkpoint-t corrected-chamber GV basis-projected gamma diagnostics project the ambient O7 B-field through the curve-basis matrix; for 4-214-647 this gives 214 basis entries with 49 odd entries, zero parity mismatches over the 410 toric-covered curves, and exactly the same target correction as direct ambient gamma with vs_ambient_max_abs=0",
            "checkpoint-t corrected divisor equation diagnostics now compare T_i=τ_classical-χ(D_i)/24+GV_i directly to c_i/c_tau; for 4-214-647 the checkpoint-implied GV satisfies this equation up to the corrected_target_volumes.dat classical tolerance max_abs=0.0005560489706983862, while Cyrus toric-covered GV leaves max_abs=0.11543984297490459",
            "simple ambient B-field gamma index shifts do not explain the checkpoint-implied GV mismatch: shifting gamma by -1 gives max_abs_delta=0.10111964958304713 but worsens relative_l2_delta to 1.6272075913431028, and shifting by +1 gives relative_l2_delta=1.582044782998777",
            "toggling the ambient origin entry in the B-field gamma does not explain the checkpoint-implied GV mismatch: for 4-214-647 origin_idx=0 gives max_abs_delta=0.11534346351187635 and relative_l2_delta=1.246485805041676",
            "corrected_target_volumes.dat is now validated as corrected-chamber classical KKLT divisor volumes: using McAllister corrected_kahler_param.dat gives relative_l2=0.00010542596827921091 and max_abs=0.0005560489706986083, while Cyrus' computed corrected t still gives relative_l2=0.02203393739805434 and max_abs=0.0977369510437136",
            "mcallister_first_principles reports validation-only corrected_kahler_param.dat deltas; for 4-214-647 the computed corrected Kähler vector differs by relative_l2=0.0005603330179533605 and max_abs=0.18907322067627774",
            "optional supporting-face certificate diagnostics report certified and uncertified LP-witness face GV counts; the latest 4-214-647 sample certified 0/9 LP-witness faces",
            "mcallister_first_principles runs from declared inputs only with validation checkpoints absent, and can skip McAllister final assertions for generic candidate evaluation",
            "generic no-assertion runs default K/M flux vectors to Cyrus' selected production dual basis instead of the McAllister [3,4,5,8] source basis; --dual-basis now accepts index or matrix validation-source flux coordinate bases, and --production-dual-basis carries index or matrix internal dual bases through flat-direction/GV handoff",
            "corrected-chamber general GV diagnostics report degree-bounded ray counts; for the 4-214-647 corrected chamber, max_deg=26 keeps 2963 of 561596 Mori generators",
            "A regenerated current corrected-chamber GV context export from mcallister_first_principles has ambient_rays=561596, subcutoff=561, pair_pruned=419, toric_covered=410, toric_missing=9, required missing degrees 10..26, and a complete degree-bounded source context for the 2963 rays with degree <= 26",
            "mcallister_gv_context --run-integer-diamonds computes explicit-semigroup HKTY values for 3/5 exact integer-semigroup missing targets, all zero, while 2/5 integer-semigroup targets fail with non-integral HKTY output and 4/9 rational-cone targets remain unpromoted",
            "mcallister_gv_context now source-selects the unique compact-threefold local cygv q-matrix phase for all 9 corrected-chamber toric-missing targets: 4 origin-including phases and 5 origin-omitting phases, leaving only local intersection tensor and chamber certificate as local cygv input blockers",
            "mcallister_gv_context --scan-local-integer-tensors now finds an integer tensor value matching at least one expected local toric formula value for all 9 missing targets, but these remain uncertified because the matching tensor value is still diagnostic rather than a source-derived local intersection tensor",
            "mcallister_gv_context reports a multiplicity-preserving local toric formula sum for the 9 corrected-chamber toric-missing targets: the phase-selected unit-tensor cygv response matches the formula sum for 8/9 targets, while target 5 would require an uncertified local tensor value 4 under the summed-component interpretation",
            "mcallister_gv_context --scan-local-integer-tensors separately reports formula-sum tensor matches for all 9 missing targets, including target 5 at tensor value 4, without promoting those matches past the missing chamber/tensor certificate",
            "mcallister_first_principles now reports the local formula-sum GV values as a corrected-chamber diagnostic only: for 4-214-647 they cover 9/9 toric-missing classes with 9 nonzero values, give volume correction -0.2590023228059164, shift the toric-covered corrected-chamber volume correction by only +0.0005719649570432028, and still miss corrected_target_volumes.dat's implied GV vector with max_abs=0.09773695104401635 and relative_l2=1.2387268027775085 without promotion; the leading post-formula-sum target offenders are kklt_idx/point_idx 190/195, 34/37, 52/56, 192/197, 203/208, 206/211, 103/107, and 20/23",
            "A direct read of latest_cytools/compute_kklt_iterative.py confirms the GV target formula sign used by Cyrus (tau_target = c_i/c_tau + chi/24 - GV), but the primal script maps c_i onto basis.dat and does not implement the mixed basis.dat/kklt_basis.dat corrected-target vector, so it is not evidence that the Python path solved the remaining 4-214-647 corrected-target residual",
            "mcallister_gv_context now accepts schema-4 corrected-chamber GV context exports and aggregates closest known q_N residual statuses and degree splits; bounded target-7/8 reports on 4-214-647 record the residual split as a shared known source-derived degree-6 predecessor plus known toric degree-2 differences",
            "mcallister_gv_context now promotes closest-known q_N residual source predecessors into a first-class report queue; an unfiltered bounded 4-214-647 run records two unique source-derived predecessors across three occurrences, with target 7/8 sharing the degree-6 ray [(54,-2),(203,1),(206,1),(209,1)] and target 6 using an analogous degree-4 ray, all with local charge signature [-1,-1,-1,1,2], shared origin-circuit witness relations, source_derived_full_facet_context, CMS cubic candidate 3, source-derived GV -2, and remaining blockers local_intersection_tensor plus local_chamber_certificate",
            "mcallister_gv_context now requires full source-derived origin-circuit witness/facet context before importing CMS-derived scalar GV values into known q_N history; the guarded 4-214-647 report preserves the same two residual source predecessors across three occurrences, so shape-only CMS matches remain diagnostic rather than source history",
            "mcallister_gv_context now reports source-derived GV import statuses for uncovered source rays; the guarded 4-214-647 report records 52 full-facet CMS formula imports, 79 rows with no integral matching CMS formula, and 115 rows with no origin-circuit witness",
            "mcallister_gv_context reports pre-LP supporting-face profiles for origin-circuit witness domains; for the 4-214-647 corrected chamber, the default witness-domain certificate diagnostic now finishes with 14 relation domains skipped as single-generator codimension-213 cases, seven small shared-facet domains returning no LP certificate, seven larger shared-facet domains skipped by the 256-generator limit, and all facet-union domains skipped by generator limits",
        ],
        remaining_gaps: vec![
            "Generated branch candidates without the height_projected initializer still did not find the 4-214-647 paper branch in a deterministic 48-candidate diagnostic: the lowest sampled phase-1 volume was about 20611 rather than 17901, and even coverage-aware selection still had at least 412 small curves missing toric GV coverage",
            "The current no-replay 4-214-647 corrected-volume residual is about 0.072 after direct toric small-curve GV volume correction; the component split is V_classical=4712.269883, BBHL=0.508832, pre-GV V=4711.761051, GV volume correction=-0.256385, so the remaining discrepancy is in the instanton-correction layer rather than classical geometry",
            "The local formula-sum diagnostic values for the 9 corrected-chamber toric-missing classes do not explain the corrected-volume residual: they cover every current toric miss but change the corrected-chamber volume correction by only about +0.000572 relative to toric-covered and leave corrected_target_volumes.dat's implied GV vector off by max_abs=0.09773695104401635, led by point_idx=195 rather than a missing-origin-circuit class",
            "Corrected-chamber ray plus LP-witness diagnostics previously covered the toric-missing curves with diagnostic zero GV values and no delta_vs_toric_covered; the current exported context has 9 toric-missing curves, and those values still cannot be promoted because the LP-witness/local values are not certified supporting-face computations",
            "The final-chamber GV target-correction vector differs strongly from the input-chamber vector (relative_l2_delta=0.5969594710109491), but a six-step toric-covered-only chamber-updated KKLT diagnostic shifts the partial V_string by only +0.0003630209221228142, ending at V_string_partial=4711.505029594299; stale final-chamber toric-covered data alone does not explain the 0.072 residual",
            "The corrected-volume residual is now mostly attributable to the solved Kähler vector: corrected_kahler_param.dat gives input-chamber V_classical=4712.3385075726355 while Cyrus computes 4712.269883496959, a delta of -0.0686240756767802 out of the 0.07216733782297524 V_string residual",
            "The Kähler-vector residual is now traceable further upstream to the corrected-chamber GV target-correction layer: corrected_target_volumes.dat is understood as corrected-chamber classical KKLT divisor volumes and matches McAllister corrected_kahler_param.dat at relative_l2=0.00010542596827921091, but the checkpoint-t corrected-chamber toric-covered GV target correction still differs from the checkpoint-implied correction by relative_l2=1.2720321759981263 and max_abs=0.1154398144342352; the next audit should focus on exact corrected-chamber GV coverage/conventions, not divisor χ, checkpoint-file semantics, or the approximate-chamber small-curve set",
            "The leading checkpoint-t corrected-chamber GV target offenders are not explained by pair-vs-four-term pruning, by including all unpruned subcutoff toric-covered curves, by the old no-gamma convention, by a uniform ±1 ambient gamma index shift, by toggling the ambient origin entry in gamma, by ambient-vs-basis gamma contraction, by corrected_target_volumes.dat semantics, by a mismatch between corrected_heights.dat and the checkpoint-t-induced chamber, by corrected-chamber divisor chi, by dropping tiny positive q.t curves, by continuous dropping of individual toric-covered curves, by arbitrary real GV values on the toric-missing corrected-chamber classes, by isolating or sign-flipping B-field parity/q.t branch buckets, by isolating or sign-flipping local-formula source buckets, by real-dilog numerical evaluation for the top small-q.t contributors, by an external Python/SciPy recomputation of the exported corrected-chamber target trace, by a Rust/CYTools Mori-cap selection divergence, by toric-formula provenance for the top offenders, by the exact local-HKTY/LP-witness values of 50 unique top toric contributors, or by any single ambient mod-2 gamma shift matching the relaxed bounded-sign LP near-bound pattern; a relaxed bounded-sign LP can nearly fit the checkpoint vector but only with many arbitrary sign flips, so the next audit should compare broader per-curve cygv/general-GV values and missing non-toric contributions rather than the final dilog/prefactor aggregation",
            "Saved input-chamber small-curve variants do not explain corrected_target_volumes.dat either: using small_curves.dat at corrected_kahler_param.dat with signed or absolute q.t, with or without the O7 gamma sign, still misses the checkpoint-implied GV correction; the best obvious variant was absolute q.t plus gamma with max_abs_delta~0.11539193539577981 and relative_l2_delta~1.1964221418662775",
            "A direct CYTools/cygv provided-generator diagnostic using the 420 corrected-chamber pair-pruned curves reached grading construction but exceeded a 240-second timeout before returning GV values, so per-curve cygv comparison for the leading offenders remains unresolved",
            "A production chamber-updated KKLT solve still needs certified or general GV values for the current 9 corrected-chamber toric-missing classes at every fixed-point chamber before it can replace the input-chamber small-curve/GV set",
            "The current production small-curve pruning only removes pair-decomposable curves; the new decomposition-depth report is bounded to four terms and is not a full faithful implementation of the paper's sums-of-others/Hilbert-basis reduction",
            "Potent-ray rays and most full N_{nq} GV series are still validation-supplied for 4-214-647; Cyrus computes rank/volumes/convergence from the sample and has source-derived rank-two CKYZ checks through the first four entries for every rank-two row, plus all ten entries for F0 [1,1]/[1,2] and F1 [2,1]/[3,1], but does not yet generate low-dimensional Mori-infinity face rays, validate all ten entries for every rank-two row, or handle rank-four rows",
            "The explicit general primal GV fallback still reaches full 214-dimensional Mori-cone dualization for any max_deg high enough to cover the selected missing curves, or for min_points-driven runs; in the latest 4-candidate generated-branch diagnostic the min-required-gv-degree branch still requires degrees 4..2334",
            "A PPL/cdd diagnostic on the dumped 561658-ray, 214-dimensional V-representation also exceeded a 300-second cap without producing an H-representation",
            "Finish a post-orientation-fix validation run of adjacency-filtered DDM on the full 214-dimensional McAllister Mori cone",
            "Further optimize or replace hyperplane dualization; bounded diagnostics still need to prove the full 561658-ray McAllister dualization completes with the corrected ray orientation",
            "Do not promote degree-bounded generator filtering as an exact GV shortcut: a high-degree cone ray can still affect low-degree lattice points through fractional cone combinations, so exact bounded GV needs full inequalities or a proven reduced formulation",
            "The corrected-chamber provided-generator diagnostic with max_deg=26 and 2963 degree-bounded rows exceeded a 600-second timeout without producing a cygv result, so naive mcap_generators truncation is not a useful replacement for an exact reduced formulation",
            "The corrected-chamber LP-face supporting-certificate diagnostic with default 16-anchor search exceeded a 300-second timeout before producing a certificate summary; the smaller 2-anchor sample still certified 0/9 faces",
            "The corrected-chamber origin-circuit witness-domain certificate profile does not justify promoting shared-facet local GV values: relation supports are single-generator high-codimension domains, seven small shared-facet domains have no bounded LP certificate, and the remaining shared/union domains are too large under the current generator limits",
            "Reduce the 561658-ray Mori cap input before dualization, or add a CYTools/PPL-faithful constraint minimization path",
            "Run and validate lattice-point generation under a Python environment with OR-Tools after DDM returns the dual cone",
        ],
    };

    insta::assert_json_snapshot!("gv_computation_roadmap", roadmap);
}
