//! Stage 11 (in-process): drive `verify_kklt_vacuum` directly.
//!
//! The heavy `stage10_volume::stage11_first_principles_runner_reaches_*` test
//! runs the whole `mcallister_first_principles` binary as a subprocess. This
//! test instead exercises the library orchestrator IN-PROCESS: it computes the
//! primal geometry and intersection numbers from first principles (the same
//! production code the binary's triangulation/intersection stages call), then
//! calls `cyrus_core::kklt_vacuum::verify_kklt_vacuum` and checks it reproduces
//! McAllister 4-214's corrected volume and `V0`.
//!
//! The orchestrator's *inputs* `g_s`, `|W0|`, `e^{K0}`, `c_i` and the KKLT
//! basis are NOT computed by the orchestrator - they are produced by the
//! upstream flat-direction/racetrack stages (validated by stage4/stage5) and
//! are passed in here from the published checkpoint. The orchestrator's own
//! scope - the corrected target, branch solve, small-curve GV coverage, the
//! instanton-corrected solve, and the BBHL/GV-corrected volume down to `V0` -
//! is computed live, so a regression in any of those fails this test.

use std::path::{Path, PathBuf};

use cyrus_core::kklt_vacuum::{
    BranchSelection, PrimalGeom, PrimalIntersection, StabilizationInputs, VacuumConfig,
    computed_primal_basis, verify_kklt_vacuum,
};
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::{Finite, Pos};
use cyrus_core::{
    CurvePruningStrategy, Point, Polytope, compute_glsm_and_linrels, compute_intersection_cytools,
    compute_linear_relations_no_origin, compute_regular_triangulation, intersection_in_basis,
};

/// Published `e^{K0}` for 4-214-647 (arXiv:2107.09064 SS6.4). The checkpoint
/// ships no `e^{K0}` scalar file; `e^{K0}` is an orchestrator input from the
/// flat-direction stage (stage4) and only enters `V0`, not the volume.
const MCALLISTER_4_214_EK0: f64 = 0.234_393;
/// `h^{2,1}` for 4-214-647 (so `χ = 2(h11 - h21) = 2(214 - 4) = 420`).
const MCALLISTER_4_214_H21: usize = 4;

fn require_first_principles() -> bool {
    if !crate::first_principles_enabled() {
        eprintln!("Skipping first-principles test (set CYRUS_FIRST_PRINCIPLES=1)");
        return false;
    }
    true
}

fn require_runner_heavy() -> bool {
    if std::env::var_os("CYRUS_MCALLISTER_RUNNER_HEAVY").is_none() {
        eprintln!("Skipping heavy orchestrator test (set CYRUS_MCALLISTER_RUNNER_HEAVY=1)");
        return false;
    }
    true
}

fn require_data_dir() -> Option<PathBuf> {
    let Some(dir) = crate::mcallister_data_dir() else {
        assert!(
            !crate::first_principles_enabled(),
            "CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests"
        );
        return None;
    };
    Some(dir)
}

fn read_floats_csv(path: &Path) -> Vec<f64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split([',', '\n', '\r'])
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<f64>().expect("invalid float"))
        .collect()
}

fn read_ints_csv(path: &Path) -> Vec<i64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split([',', '\n', '\r'])
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<i64>().expect("invalid integer"))
        .collect()
}

fn read_usize_csv(path: &Path) -> Vec<usize> {
    read_ints_csv(path)
        .into_iter()
        .map(|v| usize::try_from(v).expect("index must be non-negative"))
        .collect()
}

fn read_scalar_f64(path: &Path) -> f64 {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()))
        .trim()
        .parse::<f64>()
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

fn read_points(path: &Path) -> Vec<Vec<i64>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|s| s.trim().parse::<i64>().expect("invalid point entry"))
                .collect()
        })
        .collect()
}

/// Build the primal geometry from `points.dat` + `heights.dat` (production
/// triangulation code).
fn build_geom(data_dir: &Path) -> PrimalGeom {
    let points_raw = read_points(&data_dir.join("points.dat"));
    let heights_raw = read_floats_csv(&data_dir.join("heights.dat"));
    let heights: Vec<F64<Finite>> = heights_raw
        .iter()
        .map(|&h| F64::<Finite>::new(h).expect("height must be finite"))
        .collect();
    let points: Vec<Point> = points_raw.iter().map(|p| Point::new(p.clone())).collect();
    let polytope = Polytope::from_vertices(points).expect("failed to build polytope");
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("failed to filter triangulation points");
    let triangulation = compute_regular_triangulation(&triangulation_points, &heights_raw)
        .expect("failed to compute triangulation");
    PrimalGeom {
        polytope,
        heights,
        triangulation_points,
        triangulation,
    }
}

/// Build the primal intersection data (GLSM/linrels/basis + intersection
/// tensors) from the geometry (production intersection code).
fn build_intersection(geom: &PrimalGeom) -> PrimalIntersection {
    let (glsm, linrels, basis) =
        compute_glsm_and_linrels(&geom.triangulation_points).expect("failed GLSM/linrels");
    let points_i64: Vec<Vec<i64>> = geom
        .triangulation_points
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();
    let linrels_i64: Vec<Vec<i64>> = compute_linear_relations_no_origin(&points_i64)
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("linrel fits in i64"))
                .collect()
        })
        .collect();
    let kappa_full = compute_intersection_cytools(
        &geom.triangulation,
        &geom.triangulation_points,
        &linrels_i64,
    )
    .expect("failed intersection numbers");
    let kappa_basis = intersection_in_basis(&kappa_full, &basis);
    PrimalIntersection {
        glsm,
        linrels,
        basis,
        kappa_full,
        kappa_basis,
    }
}

#[test]
fn verify_kklt_vacuum_reaches_corrected_volume_and_v0() {
    if !require_first_principles() || !require_runner_heavy() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };

    let geom = build_geom(&data_dir);
    let intersection = build_intersection(&geom);
    let production_primal_basis = computed_primal_basis(&intersection);

    let c_i: Vec<I64<Pos>> = read_ints_csv(&data_dir.join("target_volumes.dat"))
        .into_iter()
        .map(|v| I64::<Pos>::new(v).expect("target_volumes c_i must be positive"))
        .collect();
    let kklt_basis: Vec<usize> = read_usize_csv(&data_dir.join("kklt_basis.dat"));

    let g_s = F64::<Pos>::new(read_scalar_f64(&data_dir.join("g_s.dat"))).expect("g_s positive");
    let w0 = F64::<Pos>::new(read_scalar_f64(&data_dir.join("W_0.dat"))).expect("W0 positive");
    let ek0 = F64::<Pos>::new(MCALLISTER_4_214_EK0).expect("ek0 positive");
    let small_curve_cutoff =
        F64::<Pos>::new(read_scalar_f64(&data_dir.join("small_curves_cutoff.dat")))
            .expect("small-curve cutoff positive");

    let inputs = StabilizationInputs {
        geom: &geom,
        intersection: &intersection,
        production_primal_basis: &production_primal_basis,
        c_i: &c_i,
        kklt_basis: &kklt_basis,
        g_s,
        w0,
        ek0,
        h21: MCALLISTER_4_214_H21,
    };
    // Mirrors the validated runner invocation:
    // --branch-candidates 1 --branch-height-init --branch-selection first-positive --kklt-steps 64
    let config = VacuumConfig {
        kklt_steps: 64,
        branch_candidates: 1,
        branch_seed: 42,
        branch_height_init: true,
        branch_selection: BranchSelection::FirstPositive,
        small_curve_cutoff,
        small_curve_pruning: CurvePruningStrategy::PairDecomposable,
        primal_gv_min_points: None,
        primal_gv_max_deg: None,
        max_missing_gv_impact: 0.0,
        missing_gv_abs_bound: 10,
    };

    let verdict = verify_kklt_vacuum(&inputs, &config)
        .unwrap_or_else(|e| panic!("verify_kklt_vacuum failed: {e}"));

    assert_eq!(
        verdict.small_curve_gv_count, 344,
        "expected 344 selected small-curve GVs, got {}",
        verdict.small_curve_gv_count
    );
    assert!(
        verdict.gv_controlled && verdict.deferred_missing_gv_count == 0,
        "4-214 should be fully GV-controlled with no deferred curves"
    );

    let corrected_checkpoint = read_scalar_f64(&data_dir.join("corrected_cy_vol.dat"));
    let uncorrected_checkpoint = read_scalar_f64(&data_dir.join("cy_vol.dat"));
    let v_string = verdict.v_string.get();
    assert!(
        (v_string - corrected_checkpoint).abs() < 0.1,
        "orchestrator V_string {v_string} should match corrected_cy_vol.dat {corrected_checkpoint}"
    );
    assert!(
        (v_string - corrected_checkpoint).abs() < (v_string - uncorrected_checkpoint).abs(),
        "orchestrator V_string should be closer to corrected than uncorrected checkpoint"
    );

    let v0_log10_abs = verdict.v0.get().abs().log10();
    assert!(
        (-203.0..-201.0).contains(&v0_log10_abs),
        "orchestrator log10(|V0|) should be near -202, got {v0_log10_abs}"
    );
}
