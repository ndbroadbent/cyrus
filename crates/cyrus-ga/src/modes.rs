//! Run modes beyond the single-geometry loop: the multi-polytope landscape
//! bandit and the orientifold verification emitter.

use std::io::Write as _;
use std::path::PathBuf;

use rayon::prelude::*;

use cyrus_ga::fitness::{FitnessConfig, evaluate_fitness};
use cyrus_ga::geometry::{DEFAULT_GV_MIN_POINTS, GaGeometry};
use cyrus_ga::multi::{PolytopeStats, load_pool, prepare_or_mark_dead, select_next};
use cyrus_ga::population::{GaParams, GaState};

use crate::{inject_isotropic, parse_arg_value};

/// Shared GA hyperparameter / fitness-config parsing for both modes.
pub fn parse_params_and_fitness() -> (GaParams, FitnessConfig) {
    let mut params = GaParams::default();
    if let Some(v) = parse_arg_value("--population") {
        params.population_size = v;
    }
    if let Some(v) = parse_arg_value("--frontier") {
        params.frontier_size = v;
    }
    if let Some(v) = parse_arg_value("--flux-range") {
        params.flux_range = v;
    }
    if let Some(v) = parse_arg_value("--crossover-rate") {
        params.crossover_rate = v;
    }
    if let Some(v) = parse_arg_value("--asteroid-threshold") {
        params.asteroid_threshold = v;
    }
    let mut fitness_cfg = FitnessConfig::default();
    if let Some(v) = parse_arg_value("--q-max") {
        fitness_cfg.q_max = v;
    }
    if let Some(v) = parse_arg_value("--target-log10-v0") {
        fitness_cfg.target_log10_v0 = v;
    }
    if let Some(v) = parse_arg_value("--desi-w0") {
        fitness_cfg.desi_w0 = v;
    }
    if let Some(v) = parse_arg_value("--desi-wa") {
        fitness_cfg.desi_wa = v;
    }
    if let Some(v) = parse_arg_value("--decay-constant") {
        fitness_cfg.decay_constant = v;
    }
    (params, fitness_cfg)
}

/// Multi-polytope landscape mode: a UCB bandit allocates inner-GA rounds
/// across the polytope pool; per-polytope state persists under
/// `<run-dir>/polytopes/<name>/`, the bandit summary in `summary.json`.
#[allow(clippy::too_many_lines)] // linear orchestration
pub fn run_multi(pool_path: &std::path::Path) {
    let run_dir: PathBuf = parse_arg_value::<String>("--run-dir")
        .map_or_else(|| PathBuf::from("runs/landscape"), PathBuf::from);
    let rounds: u64 = parse_arg_value("--rounds").unwrap_or(u64::MAX);
    let gens_per_round: u64 = parse_arg_value("--gens-per-round").unwrap_or(10);
    let seed: u64 = parse_arg_value("--seed").unwrap_or(20260611);
    let gv_min_points: u32 = parse_arg_value("--gv-min-points").unwrap_or(DEFAULT_GV_MIN_POINTS);
    let (params, fitness_cfg) = parse_params_and_fitness();

    let pool = load_pool(pool_path).unwrap_or_else(|e| {
        eprintln!("[ERROR] {e}");
        std::process::exit(2);
    });
    std::fs::create_dir_all(run_dir.join("polytopes")).expect("create run dir");
    let summary_path = run_dir.join("summary.json");
    let improvements_path = run_dir.join("improvements.jsonl");

    let mut stats: Vec<PolytopeStats> = if summary_path.exists() {
        let text = std::fs::read_to_string(&summary_path).expect("read summary");
        let loaded: Vec<PolytopeStats> = serde_json::from_str(&text).expect("parse summary");
        eprintln!(
            "[INFO] resumed landscape run: {} polytopes, {} rounds so far",
            loaded.len(),
            loaded.iter().map(|s| s.rounds).sum::<u64>()
        );
        loaded
    } else {
        pool.iter()
            .map(|p| PolytopeStats::new(p.name.clone()))
            .collect()
    };
    assert_eq!(stats.len(), pool.len(), "summary does not match pool file");

    let mut geometries: std::collections::HashMap<usize, GaGeometry> =
        std::collections::HashMap::new();
    let mut global_best = stats
        .iter()
        .map(|s| s.best_fitness)
        .fold(f64::NEG_INFINITY, f64::max);
    let mut total_rounds: u64 = stats.iter().map(|s| s.rounds).sum();
    let target_rounds = total_rounds.saturating_add(rounds);

    while total_rounds < target_rounds {
        let Some(idx) = select_next(&stats, total_rounds) else {
            eprintln!("[ERROR] all polytopes are dead; nothing to schedule");
            std::process::exit(2);
        };
        let record = &pool[idx];
        if let std::collections::hash_map::Entry::Vacant(slot) = geometries.entry(idx) {
            let t0 = std::time::Instant::now();
            let prep_timeout = std::time::Duration::from_secs(
                parse_arg_value("--prep-timeout-secs").unwrap_or(120),
            );
            match prepare_or_mark_dead(record, gv_min_points, prep_timeout, &mut stats[idx]) {
                Some(geom) => {
                    eprintln!(
                        "[INFO] prepared {} in {:.1?} (h21={}, {} GV curves)",
                        record.name,
                        t0.elapsed(),
                        geom.h21_primal,
                        geom.gv.len()
                    );
                    slot.insert(geom);
                }
                None => {
                    std::fs::write(&summary_path, serde_json::to_string_pretty(&stats).unwrap())
                        .expect("write summary");
                    continue;
                }
            }
        }
        let geom = &geometries[&idx];

        let poly_dir = run_dir.join("polytopes").join(&record.name);
        std::fs::create_dir_all(&poly_dir).expect("create polytope dir");
        let state_path = poly_dir.join("state.json");
        let mut state = if state_path.exists() {
            serde_json::from_str(&std::fs::read_to_string(&state_path).expect("read state"))
                .expect("parse state")
        } else {
            GaState::new(params.clone(), geom.h21_primal, seed ^ (idx as u64) << 1)
        };

        let mut round_valid = 0usize;
        for _ in 0..gens_per_round {
            let reports: Vec<_> = state
                .population
                .par_iter()
                .map(|individual| evaluate_fitness(geom, &fitness_cfg, &individual.genome))
                .collect();
            round_valid += reports.iter().filter(|r| r.tier == "valid").count();
            state.absorb_reports(reports);
            inject_isotropic(&mut state, geom);
            state.next_generation();
        }
        stats[idx].rounds += 1;
        stats[idx].evaluations = state.evaluations;
        stats[idx].valid_seen += round_valid as u64;
        if state.best_fitness > stats[idx].best_fitness {
            stats[idx].best_fitness = state.best_fitness;
        }
        if state.best_fitness > global_best {
            global_best = state.best_fitness;
            if let Some(best) = state.hall_of_fame.first() {
                let record_json = serde_json::json!({
                    "polytope": record.name,
                    "round": total_rounds,
                    "fitness": state.best_fitness,
                    "genome": best.genome,
                    "report": best.report,
                });
                let mut file = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&improvements_path)
                    .expect("open improvements log");
                writeln!(file, "{record_json}").expect("append improvement");
            }
        }
        total_rounds += 1;

        let tmp = poly_dir.join("state.json.tmp");
        std::fs::write(
            &tmp,
            serde_json::to_string(&state).expect("serialize state"),
        )
        .expect("write state");
        std::fs::rename(&tmp, &state_path).expect("replace state");
        std::fs::write(&summary_path, serde_json::to_string_pretty(&stats).unwrap())
            .expect("write summary");

        eprintln!(
            "[ROUND {:>5}] {} gens+={} evals={} valid+={} best_here={:.2} global_best={:.2} live={}/{}",
            total_rounds,
            record.name,
            gens_per_round,
            state.evaluations,
            round_valid,
            state.best_fitness,
            global_best,
            stats.iter().filter(|s| !s.dead).count(),
            stats.len()
        );
    }
    eprintln!(
        "[INFO] reached round target; resume any time with the same --run-dir ({})",
        run_dir.display()
    );
}

/// Emit a runner-compatible data directory for one GA candidate:
/// `--emit-verification-dir <dir> --polytope-file <jsonl> --polytope-name
/// <name> --k a,b,... --m a,b,...`. Constructs the automated orientifold
/// (involution, c_i, KKLT divisor set, tadpole) and writes the declared
/// inputs the McAllister first-principles runner needs for a full
/// instanton-corrected KKLT verification of the candidate.
#[allow(clippy::too_many_lines)] // linear orchestration
pub fn emit_verification_dir(out_dir: &std::path::Path) {
    let pool_path: String = parse_arg_value("--polytope-file").unwrap_or_else(|| {
        eprintln!("[ERROR] --polytope-file is required");
        std::process::exit(2);
    });
    let name: String = parse_arg_value("--polytope-name").unwrap_or_else(|| {
        eprintln!("[ERROR] --polytope-name is required");
        std::process::exit(2);
    });
    let parse_flux = |flag: &str| -> Vec<i64> {
        parse_arg_value::<String>(flag)
            .unwrap_or_else(|| {
                eprintln!("[ERROR] {flag} a,b,c,... is required");
                std::process::exit(2);
            })
            .split(',')
            .map(|x| x.trim().parse().expect("flux integer"))
            .collect()
    };
    let k = parse_flux("--k");
    let m = parse_flux("--m");

    let pool = load_pool(&PathBuf::from(pool_path)).unwrap_or_else(|e| {
        eprintln!("[ERROR] {e}");
        std::process::exit(2);
    });
    let record = pool.iter().find(|p| p.name == name).unwrap_or_else(|| {
        eprintln!("[ERROR] polytope {name} not found in pool");
        std::process::exit(2);
    });

    // Primal-side data the runner needs: points, default-chamber heights,
    // GLSM basis, intersection numbers, true Hodge numbers, orientifold.
    use cyrus_core::{Point, Polytope};
    let point_objs: Vec<Point> = record
        .points
        .iter()
        .map(|p| Point::new(p.clone()))
        .collect();
    let polytope = Polytope::from_vertices(point_objs).expect("polytope");
    let tri_points = polytope
        .points_not_interior_to_facets()
        .expect("triangulation points");
    let origin_idx = tri_points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .expect("origin");
    let (heights, triangulation) =
        cyrus_core::compute_frst_heights(&tri_points, origin_idx).expect("FRST");
    let (_, _, basis) = cyrus_core::compute_glsm_and_linrels(&tri_points).expect("glsm");
    let points_i64: Vec<Vec<i64>> = tri_points.iter().map(|p| p.coords().to_vec()).collect();
    let linrels = cyrus_core::compute_linear_relations_no_origin(&points_i64);
    let linrels_i64: Vec<Vec<i64>> = linrels
        .iter()
        .map(|r| r.iter().map(|x| i64::try_from(x).expect("i64")).collect())
        .collect();
    let kappa_full =
        cyrus_core::compute_intersection_cytools(&triangulation, &tri_points, &linrels_i64)
            .expect("intersection numbers");
    let h11_extra =
        cyrus_core::divisor::batyrev_h11_extra_classes(&polytope, &tri_points).expect("batyrev");
    let h11 = basis.len() + h11_extra;
    let h21 = record.h21.unwrap_or_else(|| {
        eprintln!("[ERROR] pool record lacks h21");
        std::process::exit(2);
    });

    let models = cyrus_core::orientifold::enumerate_involutions(
        &polytope,
        &tri_points,
        &kappa_full,
        &basis,
        h11,
        h21,
    )
    .expect("orientifold enumeration");
    let Some(model) = cyrus_core::orientifold::select_orientifold(&models) else {
        eprintln!(
            "[ERROR] no viable orientifold (15 involutions enumerated): candidate cannot be verified in the McAllister model class"
        );
        std::process::exit(2);
    };
    eprintln!(
        "[INFO] selected involution sigma={:?}: {} O7-planes, tadpole Q_D3={}, purity UNVERIFIED",
        model.sigma,
        model.o7_points.len(),
        model.q_d3
    );
    let q_flux = -0.5
        * k.iter()
            .zip(m.iter())
            .map(|(&ki, &mi)| (ki * mi) as f64)
            .sum::<f64>();
    if q_flux < 0.0 || q_flux > model.q_d3 {
        eprintln!(
            "[ERROR] candidate violates the automated tadpole bound: -M.K/2 = {q_flux} not in [0, {}]",
            model.q_d3
        );
        std::process::exit(2);
    }

    std::fs::create_dir_all(out_dir).expect("create out dir");
    let csv_i64 = |v: &[i64]| {
        v.iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(",")
    };
    let csv_usize = |v: &[usize]| {
        v.iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(",")
    };
    let points_csv: String = record
        .points
        .iter()
        .map(|p| csv_i64(p))
        .collect::<Vec<_>>()
        .join("\n");
    std::fs::write(out_dir.join("points.dat"), points_csv).expect("points");
    std::fs::write(
        out_dir.join("heights.dat"),
        heights
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(","),
    )
    .expect("heights");
    std::fs::write(out_dir.join("K_vec.dat"), csv_i64(&k)).expect("K");
    std::fs::write(out_dir.join("M_vec.dat"), csv_i64(&m)).expect("M");
    std::fs::write(out_dir.join("kklt_basis.dat"), csv_usize(&model.kklt_basis)).expect("kklt");
    std::fs::write(out_dir.join("target_volumes.dat"), csv_i64(&model.c_i)).expect("c_i");

    // The GA generated K/M in the COMPUTED dual divisor basis; declare that
    // explicitly so the runner cannot misinterpret the coordinates.
    let geom = GaGeometry::prepare_from_points(&record.points, DEFAULT_GV_MIN_POINTS)
        .expect("geometry for dual basis");
    std::fs::write(
        out_dir.join("flux_basis.json"),
        serde_json::to_string(&serde_json::json!({"indices": geom.dual_basis})).expect("json"),
    )
    .expect("flux basis");

    eprintln!(
        "[INFO] verification inputs written to {}",
        out_dir.display()
    );
    eprintln!("[INFO] run the full first-principles verification with:");
    eprintln!(
        "  CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR={d} ./target/release/mcallister_first_principles --data-dir {d} --skip-mcallister-assertions --dual-basis {d}/flux_basis.json",
        d = out_dir.display()
    );
}
