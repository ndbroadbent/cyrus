//! BENCHMARK (not a correctness gate): time the first N admissible-chamber
//! coverage probes and split each probe's cost between the phase-1 branch
//! homotopy and the small-curve GV coverage cascade. This is the measurement
//! harness for the scan-speedup work: it lets each optimization (convex-dual
//! initializer for the homotopy, dangerous-ray-only coverage for the cascade)
//! be evaluated on the first ~10 chambers (~1 min) instead of the full 47-min
//! scan.
//!
//! Gated like the heavy scan tests
//! (`CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_RUNNER_HEAVY=1 CYRUS_GENERAL_GV=1`)
//! and restricted to the 4-214-647 target. Env knobs: `CYRUS_BENCH_LIMIT`
//! (default 10) candidates probed; `CYRUS_PROBE_STEPS` (default 40, the
//! production `COVERAGE_PROBE_STEPS`).

use std::time::{Duration, Instant};

use cyrus_core::CurvePruningStrategy;
use cyrus_core::kklt_vacuum::{
    BranchSelection, StabilizationInputs, VacuumConfig, compute_small_curve_geometry,
    computed_primal_basis, primal_gv_grading, probe_chamber_coverage_timed,
};
use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::Pos;

use crate::stage11_verify_orchestrator::{
    MCALLISTER_4_214_EK0, MCALLISTER_4_214_H21, build_geom, build_intersection, build_o7_model,
    c_i_for_basis, read_scalar_f64, require_data_dir, require_first_principles, require_general_gv,
    require_runner_heavy,
};

fn env_usize(key: &str, default: usize) -> usize {
    std::env::var(key)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(default)
}

#[test]
fn bench_first_n_chamber_probes() {
    if !require_first_principles() || !require_runner_heavy() || !require_general_gv() {
        return;
    }
    let Some(data_dir) = require_data_dir() else {
        return;
    };
    if data_dir.file_name().unwrap().to_string_lossy() != "4-214-647" {
        return;
    }
    let limit = env_usize("CYRUS_BENCH_LIMIT", 10);
    let probe_steps = env_usize("CYRUS_PROBE_STEPS", 40);

    let setup_start = Instant::now();
    let geom = build_geom(&data_dir);
    let intersection = build_intersection(&geom);
    let production_primal_basis = computed_primal_basis(&intersection);
    let scgeom = compute_small_curve_geometry(&geom, &intersection).expect("small-curve geometry");
    let grading = primal_gv_grading(&scgeom, &intersection).expect("gv grading");
    let bases = cyrus_core::orientifold::admissible_kklt_bases(
        &geom.polytope,
        &geom.triangulation_points,
        0,
    )
    .expect("admissible bases");
    let o7_model = build_o7_model(&geom, &intersection);
    let g_s = F64::<Pos>::new(read_scalar_f64(&data_dir.join("g_s.dat"))).expect("g_s");
    let w0 = F64::<Pos>::new(read_scalar_f64(&data_dir.join("W_0.dat"))).expect("W0");
    let ek0 = F64::<Pos>::new(MCALLISTER_4_214_EK0).expect("ek0");
    let config = VacuumConfig {
        kklt_steps: probe_steps,
        branch_candidates: 1,
        branch_seed: 42,
        branch_height_init: true,
        branch_selection: BranchSelection::FirstPositive,
        small_curve_cutoff: F64::<Pos>::new(1.0).unwrap(),
        small_curve_pruning: CurvePruningStrategy::PairDecomposable,
        primal_gv_min_points: None,
        primal_gv_max_deg: None,
        max_missing_gv_impact: 0.0,
        missing_gv_abs_bound: 10,
    };
    eprintln!(
        "[BENCH] setup {:.1}s; probing first {limit} of {} bases at {probe_steps} steps",
        setup_start.elapsed().as_secs_f64(),
        bases.len()
    );

    let (mut total_homotopy, mut total_coverage) = (Duration::ZERO, Duration::ZERO);
    let (mut indeterminate, mut covered) = (0usize, 0usize);
    let scan_start = Instant::now();
    for (rank, basis) in bases.iter().take(limit).enumerate() {
        let c_i = c_i_for_basis(&o7_model, basis);
        let inputs = StabilizationInputs {
            geom: &geom,
            intersection: &intersection,
            production_primal_basis: &production_primal_basis,
            c_i: &c_i,
            kklt_basis: basis,
            g_s,
            w0,
            ek0,
            h21: MCALLISTER_4_214_H21,
        };
        let (report, timings) =
            probe_chamber_coverage_timed(&scgeom, &inputs, &config, &grading, probe_steps);
        total_homotopy += timings.homotopy;
        total_coverage += timings.coverage;
        if report.indeterminate {
            indeterminate += 1;
        } else if report.uncovered_count == 0 {
            covered += 1;
        }
        eprintln!(
            "[BENCH] rank {rank}: total={:>5}ms homotopy={:>5}ms coverage={:>5}ms uncovered={} indeterminate={}",
            (timings.homotopy + timings.coverage).as_millis(),
            timings.homotopy.as_millis(),
            timings.coverage.as_millis(),
            report.uncovered_count,
            report.indeterminate,
        );
    }
    let wall = scan_start.elapsed();
    let wall_ms = wall.as_millis();
    let pct = |d: Duration| {
        if wall_ms == 0 {
            0
        } else {
            100 * d.as_millis() / wall_ms
        }
    };
    let mean = wall / u32::try_from(limit).unwrap_or(1).max(1);
    eprintln!(
        "[BENCH] first {limit}: wall={:.1}s homotopy={:.1}s ({}%) coverage={:.1}s ({}%) mean={}ms/probe covered={covered} indeterminate={indeterminate}",
        wall.as_secs_f64(),
        total_homotopy.as_secs_f64(),
        pct(total_homotopy),
        total_coverage.as_secs_f64(),
        pct(total_coverage),
        mean.as_millis(),
    );
}
