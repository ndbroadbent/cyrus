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

/// De-risking the GA deep-verify integration: how much of the orientifold
/// model that `verify_kklt_vacuum` needs (the involution, `c_i`, `kklt_basis`)
/// can be DERIVED from the geometry rather than read from McAllister's data
/// files?
///
/// Established here for 4-214-647:
/// 1. `enumerate_involutions` finds the correct involution `sigma=[1,0,0,0]`
///    with 51 O7 divisors (matching the runner's own gamma derivation).
/// 2. The `c_i` rule - `6` on O7-parity divisors, else the irreducible-
///    component count - reproduces `target_volumes.dat` EXACTLY on
///    McAllister's `kklt_basis`.
/// 3. The remaining gap is basis SELECTION: the rigid prime toric divisors
///    number 216, and McAllister's `kklt_basis` is those minus exactly
///    `{46, 130}` (both rigid; 2 is correctly excluded as non-rigid). So
///    `kklt_basis` is a 214-element H^{1,1} basis chosen from the rigid prime
///    divisors; the deterministic drop rule (and whether V0 is invariant to
///    it) is the open piece tracked for the GA integration.
#[test]
fn derived_orientifold_model_vs_mcallister_data() {
    if !require_first_principles() || !require_runner_heavy() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    let geom = build_geom(&data_dir);
    let intersection = build_intersection(&geom);

    let h11_extra =
        cyrus_core::divisor::batyrev_h11_extra_classes(&geom.polytope, &geom.triangulation_points)
            .expect("batyrev h11 extra");
    let h11 = intersection.basis.len() + h11_extra;

    let models = cyrus_core::orientifold::enumerate_involutions(
        &geom.polytope,
        &geom.triangulation_points,
        &intersection.kappa_full,
        &intersection.basis,
        h11,
        MCALLISTER_4_214_H21,
    )
    .expect("enumerate involutions");

    // (1) The correct involution is found.
    let model = models
        .iter()
        .find(|m| m.sigma == vec![1, 0, 0, 0])
        .expect("sigma=[1,0,0,0] model present");
    assert_eq!(model.o7_points.len(), 51, "expected 51 O7 divisors");

    let data_kklt_basis: Vec<usize> = read_usize_csv(&data_dir.join("kklt_basis.dat"));
    let data_c_i: Vec<i64> = read_ints_csv(&data_dir.join("target_volumes.dat"));

    let c_i_rule = |idx: usize| -> i64 {
        if model.gamma[idx] == 1 {
            6
        } else {
            i64::try_from(model.components[idx]).unwrap()
        }
    };

    // (2) The c_i rule reproduces target_volumes.dat on McAllister's basis.
    let c_i_on_data: Vec<i64> = data_kklt_basis.iter().map(|&idx| c_i_rule(idx)).collect();
    assert_eq!(
        c_i_on_data, data_c_i,
        "derived c_i rule must reproduce target_volumes.dat on McAllister's kklt_basis"
    );

    // (3) Basis-selection gap: rigid prime divisors minus McAllister's basis
    // is exactly {46, 130}, and every McAllister divisor is rigid.
    let mut rigid: Vec<usize> = model
        .divisor_classes
        .iter()
        .filter(|(_, class)| class.is_rigid())
        .map(|(idx, _)| *idx)
        .collect();
    rigid.sort_unstable();
    let dropped: Vec<usize> = rigid
        .iter()
        .filter(|i| !data_kklt_basis.contains(i))
        .copied()
        .collect();
    let extra: Vec<usize> = data_kklt_basis
        .iter()
        .filter(|i| !rigid.contains(i))
        .copied()
        .collect();
    assert_eq!(dropped, vec![46, 130], "rigid-minus-McAllister drop set");
    assert!(extra.is_empty(), "every McAllister kklt divisor is rigid");
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
