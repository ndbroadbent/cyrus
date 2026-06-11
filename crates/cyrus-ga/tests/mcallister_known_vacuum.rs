//! Integration test: McAllister's published 4-214-647 flux pair must come
//! out of the GA fitness pipeline as a valid vacuum with the right physics.
//!
//! Requires `CYRUS_MCALLISTER_DATA_DIR` pointing at the 4-214-647 paper-data
//! directory (skipped otherwise, like the cyrus-core e2e suite). For this
//! geometry Cyrus's computed dual divisor basis coincides with the basis
//! the published `K_vec.dat`/`M_vec.dat` are written in, so the files can
//! be used as genomes directly.

use cyrus_ga::fitness::{FitnessConfig, evaluate_fitness};
use cyrus_ga::genome::Genome;
use cyrus_ga::geometry::GaGeometry;

fn read_flux(path: &std::path::Path) -> Vec<i64> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()))
        .split(',')
        .map(|cell| cell.trim().parse::<i64>().expect("integer flux"))
        .collect()
}

#[test]
fn mcallister_4_214_flux_pair_is_a_valid_high_fitness_vacuum() {
    let Ok(data_dir) = std::env::var("CYRUS_MCALLISTER_DATA_DIR") else {
        eprintln!("CYRUS_MCALLISTER_DATA_DIR not set; skipping");
        return;
    };
    let data_dir = std::path::PathBuf::from(data_dir);
    if !data_dir.join("points.dat").exists() {
        eprintln!("points.dat not found; skipping");
        return;
    }

    let geom = GaGeometry::prepare_from_data_dir(&data_dir, 20_000).expect("geometry");
    // The published flux files are written in the CYTools-2021 basis
    // [3, 4, 5, 8]; transform into the computed basis like the runner.
    let (k, m) = geom
        .flux_pair_from_index_basis(
            &[3, 4, 5, 8],
            &read_flux(&data_dir.join("K_vec.dat")),
            &read_flux(&data_dir.join("M_vec.dat")),
        )
        .expect("flux basis transform");
    let genome = Genome { k, m };
    let cfg = FitnessConfig::default();
    let report = evaluate_fitness(&geom, &cfg, &genome);

    assert_eq!(report.tier, "valid", "report: {report:?}");
    let g_s = report.g_s.expect("g_s");
    assert!(
        (g_s - 0.00911134).abs() / 0.00911134 < 1e-3,
        "g_s = {g_s}, expected ~0.00911134"
    );
    let w0 = report.w0.expect("w0");
    assert!(
        (w0.log10() - (2.3e-90f64).log10()).abs() < 0.1,
        "W0 = {w0:e}, expected ~2.3e-90"
    );
    let log_v0 = report.log10_abs_v0.expect("V0");
    assert!(
        (-215.0..=-195.0).contains(&log_v0),
        "log10|V0| = {log_v0}, expected near -203 (proxy-volume chain)"
    );
    // Deep AdS scale: the quintessence integration must NOT reward this
    // candidate's slope (it freezes at w = -1 or fails to dominate).
    if let (Some(w0_cpl), Some(_)) = (report.cpl_w0, report.cpl_wa) {
        assert!(
            (w0_cpl - -1.0).abs() < 0.1,
            "deep vacuum should look like a frozen field, got w0 = {w0_cpl}"
        );
    }
}
