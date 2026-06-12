//! Is a float-gate "orthogonality violation" actually an EXACT isotropic
//! candidate corrupted by solve conditioning? Gated on CYRUS_RACETRACK_DIAG.

use cyrus_ga::geometry::{DEFAULT_GV_MIN_POINTS, GaGeometry};
use cyrus_ga::multi::load_pool;
use cyrus_ga::pfv::{is_isotropic, n_inverse_for_m};

#[test]
fn exact_check() {
    let Ok(spec) = std::env::var("CYRUS_RACETRACK_DIAG") else {
        eprintln!("skipping");
        return;
    };
    let parts: Vec<&str> = spec.split('|').collect();
    let parse = |s: &str| -> Vec<i64> { s.split(',').map(|x| x.parse().unwrap()).collect() };
    let (k, m) = (parse(parts[2]), parse(parts[3]));
    let pool = load_pool(std::path::Path::new(parts[0])).expect("pool");
    let record = pool.iter().find(|p| p.name == parts[1]).expect("polytope");
    let geom =
        GaGeometry::prepare_from_points(&record.points, DEFAULT_GV_MIN_POINTS).expect("geometry");
    let n_inv = n_inverse_for_m(&geom.kappa_basis, &m).expect("N invertible");
    println!(
        "EXACT isotropy of K={k:?} M={m:?}: {}",
        is_isotropic(&n_inv, &k)
    );
}
