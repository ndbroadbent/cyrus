//! How dense are EXACT isotropic flux pairs on a given geometry? Counts
//! isotropic K in the full box for a few random M, and reports the kappa
//! scale (the float K.p of exact solutions, i.e. the solve noise floor).
//! Gated on CYRUS_RACETRACK_DIAG = "<pool>|<name>|<flux_range>|<n_m>".

use cyrus_ga::geometry::{DEFAULT_GV_MIN_POINTS, GaGeometry};
use cyrus_ga::multi::load_pool;
use cyrus_ga::pfv::{is_isotropic, n_inverse_for_m};
use rand::{Rng, SeedableRng};

#[test]
fn isotropic_density() {
    let Ok(spec) = std::env::var("CYRUS_RACETRACK_DIAG") else {
        eprintln!("skipping");
        return;
    };
    let parts: Vec<&str> = spec.split('|').collect();
    let range: i64 = parts[2].parse().unwrap();
    let n_m: usize = parts[3].parse().unwrap();
    let pool = load_pool(std::path::Path::new(parts[0])).expect("pool");
    let record = pool.iter().find(|p| p.name == parts[1]).expect("polytope");
    let geom =
        GaGeometry::prepare_from_points(&record.points, DEFAULT_GV_MIN_POINTS).expect("geometry");

    // kappa scale
    let mut max_kappa: f64 = 0.0;
    for a in 0..4 {
        for b in 0..4 {
            for c in 0..4 {
                let v: f64 = geom
                    .kappa_basis
                    .get(a, b, c)
                    .to_string()
                    .parse()
                    .unwrap_or(0.0);
                max_kappa = max_kappa.max(v.abs());
            }
        }
    }
    println!("max |kappa| entry: {max_kappa}");

    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(42);
    let side = 2 * range + 1;
    let total = (side as u64).pow(4);
    for trial in 0..n_m {
        let m: Vec<i64> = (0..4).map(|_| rng.gen_range(-range..=range)).collect();
        if m.iter().all(|&x| x == 0) {
            continue;
        }
        let Some(n_inv) = n_inverse_for_m(&geom.kappa_basis, &m) else {
            println!("M={m:?}: N singular");
            continue;
        };
        let mut found = 0u64;
        let mut example: Option<Vec<i64>> = None;
        for k0 in -range..=range {
            for k1 in -range..=range {
                for k2 in -range..=range {
                    for k3 in -range..=range {
                        let k = [k0, k1, k2, k3];
                        if k.iter().all(|&x| x == 0) {
                            continue;
                        }
                        if is_isotropic(&n_inv, &k) {
                            found += 1;
                            if example.is_none() {
                                example = Some(k.to_vec());
                            }
                        }
                    }
                }
            }
        }
        println!(
            "M={m:?}: {found}/{total} isotropic K in box (example: {example:?}), trial {trial}"
        );
    }
}
