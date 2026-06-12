//! Where do exact-isotropic seeds land in the evaluation chain? Gated on
//! CYRUS_RACETRACK_DIAG = "<pool>|<name>".

use std::collections::BTreeMap;

use cyrus_ga::fitness::{FitnessConfig, evaluate_fitness};
use cyrus_ga::geometry::{DEFAULT_GV_MIN_POINTS, GaGeometry};
use cyrus_ga::multi::load_pool;

#[test]
fn seed_tier_histogram() {
    let Ok(spec) = std::env::var("CYRUS_RACETRACK_DIAG") else {
        eprintln!("skipping");
        return;
    };
    let parts: Vec<&str> = spec.split('|').collect();
    let pool = load_pool(std::path::Path::new(parts[0])).expect("pool");
    let record = pool.iter().find(|p| p.name == parts[1]).expect("polytope");
    let geom =
        GaGeometry::prepare_from_points(&record.points, DEFAULT_GV_MIN_POINTS).expect("geometry");
    println!("{}: {} seeds", parts[1], geom.pfv_seeds.len());
    let cfg = FitnessConfig::default();
    let mut hist: BTreeMap<String, usize> = BTreeMap::new();
    let mut best: Option<(f64, String)> = None;
    for seed in &geom.pfv_seeds {
        let report = evaluate_fitness(&geom, &cfg, seed);
        let key: String = report.tier.chars().take(40).collect();
        *hist.entry(key).or_default() += 1;
        if best.as_ref().is_none_or(|(f, _)| report.fitness > *f) {
            best = Some((
                report.fitness,
                format!(
                    "K={:?} M={:?} tier={} logV0={:?}",
                    seed.k, seed.m, report.tier, report.log10_abs_v0
                ),
            ));
        }
    }
    for (tier, n) in &hist {
        println!("  {n:>3} x {tier}");
    }
    if let Some((f, desc)) = best {
        println!("best seed: fitness={f:.2} {desc}");
    }
}
