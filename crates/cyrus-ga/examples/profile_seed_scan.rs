//! Profile where the PFV seed-scan time goes on a barren high-h21 geometry.
//!
//! Splits the per-M scan into three buckets - adjugate construction,
//! isotropic-K enumeration, and the keep() gate - so we can see which one
//! dominates before optimizing. Run with a polytope JSONL record path:
//!
//!   cargo run --release -p cyrus-ga --example profile_seed_scan -- /tmp/barren_7.jsonl
//!
//! Replicates the production scan parameters exactly (m_range=15, k_box=4,
//! max_m_tries=1500, max_seeds=64, seed 0xC1B05).

use std::time::Instant;

use rand::Rng as _;
use rand::SeedableRng as _;
use rand_chacha::ChaCha8Rng;

use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::Finite;
use cyrus_ga::genome::Genome;
use cyrus_ga::geometry::{DEFAULT_GV_MIN_POINTS, GaGeometry};
use cyrus_ga::isotropic_enum::enumerate_isotropic_in_box;
use cyrus_ga::pfv::integer_adjugate_for_m;

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: profile_seed_scan <record.jsonl>");
    let line = std::fs::read_to_string(&path).expect("read record");
    let record: serde_json::Value = serde_json::from_str(&line).expect("parse json");
    let name = record["name"].as_str().unwrap_or("?");
    let points: Vec<Vec<i64>> = record["points"]
        .as_array()
        .expect("points array")
        .iter()
        .map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap())
                .collect()
        })
        .collect();

    eprintln!("[prep] building geometry for {name} (gv_min={DEFAULT_GV_MIN_POINTS})...");
    let t_prep = Instant::now();
    let geom = GaGeometry::prepare_from_points_in_chamber(&points, DEFAULT_GV_MIN_POINTS, 0, false)
        .expect("prepare geometry");
    eprintln!(
        "[prep] done in {:.2}s | h21_primal={} | gv_curves={}",
        t_prep.elapsed().as_secs_f64(),
        geom.h21_primal,
        geom.gv.len()
    );

    // Reconstruct the production keep() gate from public fields.
    let kappa = &geom.kappa_basis;
    let gv = &geom.gv;
    let q_d3 = geom.q_d3;
    let keep = |genome: &Genome| -> bool {
        let q: f64 = -0.5
            * genome
                .k
                .iter()
                .zip(genome.m.iter())
                .map(|(&ki, &mi)| (ki * mi) as f64)
                .sum::<f64>();
        if !(0.0..=q_d3).contains(&q) {
            return false;
        }
        let m_typed: Vec<I64<Finite>> = genome.m.iter().map(|&x| I64::<Finite>::new(x)).collect();
        let k_typed: Vec<I64<Finite>> = genome.k.iter().map(|&x| I64::<Finite>::new(x)).collect();
        let n_mat = cyrus_core::flat_direction::compute_n_matrix(kappa, &m_typed);
        let Some(p) = cyrus_core::flat_direction::solve_linear_system_faer(&n_mat, &k_typed) else {
            return false;
        };
        gv.iter().all(|inv| {
            let a: f64 = inv
                .curve
                .iter()
                .zip(p.iter())
                .map(|(qi, pi)| qi.to_f64().get() * pi.get())
                .sum();
            a > 1e-8
        })
    };

    // Production scan parameters (geometry.rs:286-295).
    let dim = geom.h21_primal;
    let m_range = 15i64;
    let k_box = 4i64;
    let max_m_tries = 1500usize;
    let max_seeds = 64usize;
    const PER_M_CAP: usize = 4_000_000;
    let mut rng = ChaCha8Rng::seed_from_u64(0x000C_1B05);

    let mut t_adj = std::time::Duration::ZERO;
    let mut t_enum = std::time::Duration::ZERO;
    let mut t_keep = std::time::Duration::ZERO;
    let mut n_m_used = 0usize;
    let mut n_singular = 0usize;
    let mut total_iso = 0u64;
    let mut total_kept = 0u64;

    let t_all = Instant::now();
    let mut seeds: Vec<Genome> = Vec::new();
    'outer: for _ in 0..max_m_tries {
        if seeds.len() >= max_seeds {
            break;
        }
        let m: Vec<i64> = (0..dim)
            .map(|_| rng.gen_range(-m_range..=m_range))
            .collect();
        if m.iter().all(|&x| x == 0) {
            continue;
        }
        n_m_used += 1;

        let t0 = Instant::now();
        let adj = integer_adjugate_for_m(kappa, &m);
        t_adj += t0.elapsed();
        let Some(adj) = adj else {
            n_singular += 1;
            continue;
        };

        let t1 = Instant::now();
        let ks = enumerate_isotropic_in_box(&adj, dim, k_box, PER_M_CAP);
        t_enum += t1.elapsed();
        total_iso += ks.len() as u64;

        for k in ks {
            let genome = Genome { k, m: m.clone() };
            let t2 = Instant::now();
            let kept = keep(&genome);
            t_keep += t2.elapsed();
            if kept {
                total_kept += 1;
                seeds.push(genome);
                if seeds.len() >= max_seeds {
                    break 'outer;
                }
            }
        }
    }
    let wall = t_all.elapsed();

    let pct = |d: std::time::Duration| 100.0 * d.as_secs_f64() / wall.as_secs_f64();
    println!("\n=== seed-scan profile: {name} (dim={dim}) ===");
    println!("M tries used: {n_m_used} (singular skipped: {n_singular})");
    println!("isotropic K enumerated: {total_iso} | kept (passed keep): {total_kept}");
    println!("seeds found: {} / {max_seeds}", seeds.len());
    println!("wall: {:.3}s", wall.as_secs_f64());
    println!(
        "  adjugate:  {:>7.3}s  ({:>5.1}%)",
        t_adj.as_secs_f64(),
        pct(t_adj)
    );
    println!(
        "  enumerate: {:>7.3}s  ({:>5.1}%)",
        t_enum.as_secs_f64(),
        pct(t_enum)
    );
    println!(
        "  keep:      {:>7.3}s  ({:>5.1}%)",
        t_keep.as_secs_f64(),
        pct(t_keep)
    );
    if n_m_used > 0 {
        println!(
            "per-M: enumerate {:.3}ms | adjugate {:.3}ms",
            1000.0 * t_enum.as_secs_f64() / n_m_used as f64,
            1000.0 * t_adj.as_secs_f64() / n_m_used as f64
        );
    }
}
