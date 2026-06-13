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
use cyrus_ga::pfv::sample_isotropic_genome;

mod modes;
use cyrus_ga::population::{GaParams, GaState};
use modes::{emit_verification_dir, run_multi};

pub(crate) fn parse_arg_value<T: std::str::FromStr>(name: &str) -> Option<T> {
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
    // Hidden worker mode: prepare one polytope's geometry and exit 0/1.
    // The landscape scheduler runs this as a KILLABLE subprocess probe so a
    // pathological geometry cannot leak a runaway thread into the service
    // (GV/lattice caches persist to disk, so the parent's re-prep is warm).
    if let Some(name) = parse_arg_value::<String>("--prep-probe") {
        let pool_path: String = parse_arg_value("--polytope-file").unwrap_or_else(|| {
            eprintln!("[ERROR] --prep-probe requires --polytope-file");
            std::process::exit(2);
        });
        let gv_min_points: u32 =
            parse_arg_value("--gv-min-points").unwrap_or(DEFAULT_GV_MIN_POINTS);
        let pool =
            cyrus_ga::multi::load_pool(std::path::Path::new(&pool_path)).unwrap_or_else(|e| {
                eprintln!("[ERROR] {e}");
                std::process::exit(2);
            });
        let Some(record) = pool.iter().find(|p| p.name == name) else {
            eprintln!("[ERROR] polytope {name} not in pool");
            std::process::exit(2);
        };
        let chamber: usize = parse_arg_value("--chamber").unwrap_or(record.chamber);
        match GaGeometry::prepare_from_points_in_chamber(
            &record.points,
            gv_min_points,
            chamber,
            false,
        ) {
            Ok(_) => std::process::exit(0),
            Err(reason) => {
                eprintln!("[PROBE] {name}: {reason}");
                std::process::exit(1);
            }
        }
    }
    if let Some(out_dir) = parse_arg_value::<String>("--emit-verification-dir") {
        emit_verification_dir(&PathBuf::from(out_dir));
        return;
    }
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
    let target_generation = state.generation.saturating_add(generations);
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

/// Lamarckian-style constructive injection shared by both run modes.
pub(crate) fn inject_isotropic(state: &mut GaState, geom: &GaGeometry) {
    const INJECT_PER_GEN: usize = 32;
    let flux_range = state.params.flux_range;
    use rand::Rng as _;
    use rand::SeedableRng as _;
    let mut sampler_rng = rand_chacha::ChaCha8Rng::seed_from_u64(state.rng.r#gen());
    // Draw from the geometry's precomputed exact-isotropic seeds. The seeds
    // were filtered to be cone-INTERIOR PFVs (q.p > 0 for every mirror GV
    // curve), and only POSITIVE integer scaling of K preserves both the
    // isotropy (K.N(M)^{-1}.K stays 0) AND the cone-interiority (q.p stays
    // positive). Sign-flipping K or M sends p -> -p, flipping every q.p
    // negative - it threw the carefully-found PFVs straight back outside
    // the cone. Diversity comes from the distinct seeds, not from sign
    // variants. (Random K sampling provably fails where the form is
    // anisotropic for most M; the seed scan solved that at prep time.)
    let genomes: Vec<Genome> = if geom.pfv_seeds.is_empty() {
        (0..INJECT_PER_GEN)
            .filter_map(|_| sample_isotropic_genome(geom, &mut sampler_rng, flux_range, 16, 64))
            .collect()
    } else {
        (0..INJECT_PER_GEN)
            .map(|_| {
                let seed = &geom.pfv_seeds[sampler_rng.gen_range(0..geom.pfv_seeds.len())];
                // n = 1 mostly (the raw PFV), with occasional 2x (tadpole
                // grows quadratically, so larger scalings rarely fit).
                let scale: i64 = if sampler_rng.gen_ratio(1, 4) { 2 } else { 1 };
                Genome {
                    k: seed.k.iter().map(|&x| x * scale).collect(),
                    m: seed.m.clone(),
                }
            })
            .collect()
    };
    state.inject_genomes(genomes);
}
