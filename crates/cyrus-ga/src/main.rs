//! GA landscape search runner with persistent, resumable run directories.
//!
//! Start (or resume) a run:
//! ```text
//! cargo run --release -p cyrus-ga -- \
//!     --data-dir <mcallister-example-dir> --run-dir runs/my-search \
//!     [--generations N] [--population N] [--seed N] [--q-max Q] ...
//! ```
//! State is checkpointed to `<run-dir>/state.json` after every generation,
//! so the process can be stopped (Ctrl-C, kill, laptop sleep) at any time
//! and resumed by re-running the same command: if `state.json` exists it is
//! loaded and the search continues exactly where it left off (the RNG state
//! is part of the checkpoint). Improvements append to
//! `<run-dir>/improvements.jsonl`; the hall of fame lives in
//! `<run-dir>/hall_of_fame.json`.

use std::io::Write as _;
use std::path::PathBuf;

use rayon::prelude::*;

use cyrus_ga::fitness::{FitnessConfig, evaluate_fitness};
use cyrus_ga::genome::Genome;
use cyrus_ga::geometry::{DEFAULT_GV_MIN_POINTS, GaGeometry};
use cyrus_ga::multi::{PolytopeStats, load_pool, prepare_or_mark_dead, select_next};
use cyrus_ga::pfv::sample_isotropic_genome;
use cyrus_ga::population::{GaParams, GaState};

fn parse_arg_value<T: std::str::FromStr>(name: &str) -> Option<T> {
    let args: Vec<String> = std::env::args().collect();
    args.iter()
        .position(|arg| arg == name)
        .and_then(|idx| args.get(idx + 1))
        .and_then(|raw| raw.parse::<T>().ok())
}

#[derive(serde::Serialize)]
struct ImprovementRecord<'a> {
    generation: u64,
    evaluations: u64,
    fitness: f64,
    genome: &'a Genome,
    report: &'a cyrus_ga::fitness::FitnessReport,
}

#[allow(clippy::too_many_lines)] // linear CLI orchestration
fn main() {
    if let Some(pool_path) = parse_arg_value::<String>("--polytope-file") {
        run_multi(&PathBuf::from(pool_path));
        return;
    }
    let data_dir: String = parse_arg_value("--data-dir").unwrap_or_else(|| {
        eprintln!(
            "[ERROR] --data-dir <dir with points.dat> (or --polytope-file <jsonl>) is required"
        );
        std::process::exit(2);
    });
    let run_dir: PathBuf = parse_arg_value::<String>("--run-dir")
        .map_or_else(|| PathBuf::from("runs/default"), PathBuf::from);
    let generations: u64 = parse_arg_value("--generations").unwrap_or(u64::MAX);
    let seed: u64 = parse_arg_value("--seed").unwrap_or(20260611);
    let gv_min_points: u32 = parse_arg_value("--gv-min-points").unwrap_or(DEFAULT_GV_MIN_POINTS);

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

    std::fs::create_dir_all(&run_dir).unwrap_or_else(|e| {
        eprintln!("[ERROR] cannot create run dir {}: {e}", run_dir.display());
        std::process::exit(2);
    });
    let state_path = run_dir.join("state.json");
    let config_path = run_dir.join("config.json");
    let improvements_path = run_dir.join("improvements.jsonl");
    let hall_path = run_dir.join("hall_of_fame.json");

    eprintln!("[INFO] preparing geometry from {data_dir} (one-time cost)...");
    let t0 = std::time::Instant::now();
    let geom = GaGeometry::prepare_from_data_dir(std::path::Path::new(&data_dir), gv_min_points)
        .unwrap_or_else(|e| {
            eprintln!("[ERROR] geometry preparation failed: {e}");
            std::process::exit(2);
        });
    eprintln!(
        "[INFO] geometry ready in {:.1?}: h21={} (flux vectors), {} mirror GV curves, dual basis {:?}",
        t0.elapsed(),
        geom.h21_primal,
        geom.gv.len(),
        geom.dual_basis
    );

    // Resume if a checkpoint exists; otherwise start fresh.
    let mut state = if state_path.exists() {
        let text = std::fs::read_to_string(&state_path).unwrap_or_else(|e| {
            eprintln!("[ERROR] cannot read {}: {e}", state_path.display());
            std::process::exit(2);
        });
        let state: GaState = serde_json::from_str(&text).unwrap_or_else(|e| {
            eprintln!("[ERROR] cannot parse {}: {e}", state_path.display());
            std::process::exit(2);
        });
        eprintln!(
            "[INFO] resumed from {} at generation {} ({} evaluations, best fitness {:.3})",
            state_path.display(),
            state.generation,
            state.evaluations,
            state.best_fitness
        );
        state
    } else {
        eprintln!("[INFO] starting fresh run in {}", run_dir.display());
        let mut state = GaState::new(params, geom.h21_primal, seed);
        // Optional: inject the data-dir's published flux pair (transformed
        // from --seed-flux-basis coordinates) so known vacua anchor the
        // search and validate the pipeline end to end.
        if let Some(basis_csv) = parse_arg_value::<String>("--seed-flux-basis") {
            let source: Vec<usize> = basis_csv
                .split(',')
                .map(|x| x.trim().parse().expect("basis index"))
                .collect();
            let read_flux = |name: &str| -> Vec<i64> {
                std::fs::read_to_string(std::path::Path::new(&data_dir).join(name))
                    .unwrap_or_else(|e| panic!("read {name}: {e}"))
                    .split(',')
                    .map(|cell| cell.trim().parse().expect("flux integer"))
                    .collect()
            };
            let (k, m) = geom
                .flux_pair_from_index_basis(
                    &source,
                    &read_flux("K_vec.dat"),
                    &read_flux("M_vec.dat"),
                )
                .unwrap_or_else(|e| {
                    eprintln!("[ERROR] seed flux transform failed: {e}");
                    std::process::exit(2);
                });
            eprintln!("[INFO] seeding known flux pair K={k:?} M={m:?}");
            state.population[0].genome = Genome { k, m };
        }
        state
    };
    std::fs::write(
        &config_path,
        serde_json::to_string_pretty(&serde_json::json!({
            "data_dir": data_dir,
            "seed": seed,
            "params": state.params,
            "fitness": fitness_cfg,
        }))
        .expect("config serializes"),
    )
    .expect("write config");

    let mut best_seen = state.best_fitness;
    let target_generation = state.generation + generations;
    while state.generation < target_generation {
        let gen_start = std::time::Instant::now();
        let reports: Vec<_> = state
            .population
            .par_iter()
            .map(|individual| evaluate_fitness(&geom, &fitness_cfg, &individual.genome))
            .collect();
        let valid = reports.iter().filter(|r| r.tier == "valid").count();
        state.absorb_reports(reports);
        inject_isotropic(&mut state, &geom);

        if state.best_fitness > best_seen {
            best_seen = state.best_fitness;
            if let Some(best) = state.hall_of_fame.first() {
                let record = ImprovementRecord {
                    generation: state.generation,
                    evaluations: state.evaluations,
                    fitness: state.best_fitness,
                    genome: &best.genome,
                    report: best.report.as_ref().expect("hall of fame is evaluated"),
                };
                let mut file = std::fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(&improvements_path)
                    .expect("open improvements log");
                writeln!(
                    file,
                    "{}",
                    serde_json::to_string(&record).expect("record serializes")
                )
                .expect("append improvement");
            }
        }

        let best_report = state.hall_of_fame.first().and_then(|b| b.report.as_ref());
        eprintln!(
            "[GEN {:>5}] evals={} valid={}/{} best={:.3} g_s={} log10|V0|={} cpl(w0,wa)=({},{}) stagnation={} [{:.2?}]",
            state.generation,
            state.evaluations,
            valid,
            state.population.len(),
            state.best_fitness,
            best_report
                .and_then(|r| r.g_s)
                .map_or("-".into(), |v| format!("{v:.5}")),
            best_report
                .and_then(|r| r.log10_abs_v0)
                .map_or("-".into(), |v| format!("{v:.1}")),
            best_report
                .and_then(|r| r.cpl_w0)
                .map_or("-".into(), |v| format!("{v:.2}")),
            best_report
                .and_then(|r| r.cpl_wa)
                .map_or("-".into(), |v| format!("{v:.2}")),
            state.global_stagnation,
            gen_start.elapsed()
        );

        state.next_generation();

        // Checkpoint after every generation: stop/kill at any time is safe.
        let serialized = serde_json::to_string(&state).expect("state serializes");
        let tmp_path = run_dir.join("state.json.tmp");
        std::fs::write(&tmp_path, &serialized).expect("write state tmp");
        std::fs::rename(&tmp_path, &state_path).expect("atomic state replace");
        std::fs::write(
            &hall_path,
            serde_json::to_string_pretty(&state.hall_of_fame).expect("hall serializes"),
        )
        .expect("write hall of fame");
    }
    eprintln!(
        "[INFO] reached generation target; resume any time with the same --run-dir ({})",
        run_dir.display()
    );
}

/// Shared GA hyperparameter / fitness-config parsing for both modes.
fn parse_params_and_fitness() -> (GaParams, FitnessConfig) {
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
fn run_multi(pool_path: &std::path::Path) {
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
    let target_rounds = total_rounds + rounds;

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

/// Constructive PFV injection: replace the worst slice of the population
/// with genomes that satisfy the orthogonality condition K.p = 0 EXACTLY
/// (isotropic vectors of N(M)^{-1}, tested in rational arithmetic). Random
/// mutation almost never hits this measure-zero set, so the population is
/// continuously re-seeded on the constraint manifold and evolution explores
/// along it.
fn inject_isotropic(state: &mut GaState, geom: &GaGeometry) {
    const INJECT_PER_GEN: usize = 32;
    let flux_range = state.params.flux_range;
    use rand::Rng as _;
    use rand::SeedableRng as _;
    let mut sampler_rng = rand_chacha::ChaCha8Rng::seed_from_u64(state.rng.r#gen());
    let genomes: Vec<Genome> = (0..INJECT_PER_GEN)
        .filter_map(|_| sample_isotropic_genome(geom, &mut sampler_rng, flux_range, 16, 64))
        .collect();
    state.inject_genomes(genomes);
}
