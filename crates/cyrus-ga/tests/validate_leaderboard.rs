//! Validate GA leaderboard candidates end-to-end.
//!
//! For each candidate this runs the SAME pipeline the GA uses
//! (`GaGeometry::prepare_from_points` -> `evaluate_fitness` with deep-verify
//! forced on) and reports how far it gets: the cheap mirror proxy `V0`, and
//! then the full primal KKLT `deep_verify` - either a real `log10|V0|` and a
//! `controlled`/`kklt-stabilized` tier, or the exact stage it fails at
//! (recorded in `deep_tier`, no panic).
//!
//! Driven by a JSON file of candidate records via `CYRUS_VALIDATE_FILE`
//! (skipped if unset). Each record is
//! `{name, h21, k, m, points, [g_s], [w0], [log10_v0]}`, where `points` is the
//! polytope vertex/point list and `k`/`m` are the flux genome in the computed
//! dual divisor basis (exactly as stored on the leaderboard).

use cyrus_ga::fitness::{FitnessConfig, evaluate_fitness};
use cyrus_ga::genome::Genome;
use cyrus_ga::geometry::GaGeometry;
use serde::Deserialize;

#[derive(Deserialize)]
struct Candidate {
    name: String,
    h21: usize,
    k: Vec<i64>,
    m: Vec<i64>,
    points: Vec<Vec<i64>>,
    #[serde(default)]
    g_s: Option<f64>,
    #[serde(default)]
    w0: Option<f64>,
    #[serde(default)]
    log10_v0: Option<f64>,
}

#[test]
fn validate_leaderboard_candidates() {
    let Ok(file) = std::env::var("CYRUS_VALIDATE_FILE") else {
        eprintln!("CYRUS_VALIDATE_FILE not set; skipping leaderboard validation");
        return;
    };
    let raw = std::fs::read_to_string(&file).unwrap_or_else(|e| panic!("read {file}: {e}"));
    let candidates: Vec<Candidate> =
        serde_json::from_str(&raw).unwrap_or_else(|e| panic!("parse {file}: {e}"));

    // Force the primal deep-verify to run on every successful candidate,
    // regardless of its valid-band fitness.
    let cfg = FitnessConfig {
        deep_verify_threshold: Some(f64::NEG_INFINITY),
        ..FitnessConfig::default()
    };

    for c in &candidates {
        eprintln!(
            "\n[VALIDATE] === {} (h21={}) k={:?} m={:?} ===",
            c.name, c.h21, c.k, c.m
        );
        let geom = match GaGeometry::prepare_from_points(&c.points, 20_000) {
            Ok(geom) => geom,
            Err(e) => {
                eprintln!("[VALIDATE] {} prepare_from_points FAILED: {e}", c.name);
                continue;
            }
        };
        let genome = Genome {
            k: c.k.clone(),
            m: c.m.clone(),
        };
        let report = evaluate_fitness(&geom, &cfg, &genome);
        eprintln!(
            "[VALIDATE] {} tier={} fitness={:.3} q_flux={}",
            c.name, report.tier, report.fitness, report.q_flux
        );
        eprintln!(
            "[VALIDATE] {} PROXY: g_s={:?} (table {:?})  w0={:?} (table {:?})  log10|V0|={:?} (table {:?})",
            c.name, report.g_s, c.g_s, report.w0, c.w0, report.log10_abs_v0, c.log10_v0
        );
        eprintln!(
            "[VALIDATE] {} DEEP:  deep_tier={:?}  deep_log10|V0|={:?}",
            c.name, report.deep_tier, report.deep_log10_v0
        );
    }
}
