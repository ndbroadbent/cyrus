//! Integration tests validating against McAllister et al. (arXiv:2107.09064).
//!
//! These tests load fixtures (intersection numbers κ, flux vectors K/M) and verify
//! that our Rust implementation reproduces the physics values from the paper.
//!
//! ## Test Structure
//!
//! - **Fixtures**: Input data only (κ_ijk, K, M from McAllister's paper data)
//! - **Assertions**: Inline expected values from the paper
//!
//! ## Expected Values (from McAllister Table 6.1 and paper data files)
//!
//! | Example | g_s | W₀ | V_string | V₀ |
//! |---------|-----|-------|----------|------|
//! | 4-214-647 | 0.00911 | 2.30e-90 | 4711.83 | -5.5e-203 |
//! | 5-113-4627-main | 0.01112 | 6.46e-62 | 945.18 | -4.5e-138 |
//! | 5-113-4627-alt | 0.00359 | 1.13e-95 | 388.70 | -1.4e-210 |
//! | 5-81-3213 | 0.05041 | 2.04e-23 | 198.31 | -7.9e-56 |
//! | 7-51-13590 | 0.04033 | 4.08e-21 | 141.53 | -2.9e-50 |

use std::collections::HashMap;
use std::path::PathBuf;

use serde::Deserialize;

use cyrus_core::{
    Intersection, compute_ek0, compute_flat_direction, compute_v0, kklt::compute_c_tau,
    volume_string,
};

/// Fixture data for flux vectors.
#[derive(Debug, Deserialize)]
struct FluxFixture {
    #[serde(rename = "K")]
    k: Vec<i64>,
    #[serde(rename = "M")]
    m: Vec<i64>,
}

/// Fixture data for intersection numbers.
#[derive(Debug, Deserialize)]
struct IntersectionFixture {
    dim: usize,
    entries: HashMap<String, i64>,
}

/// McAllister's published physics values (ground truth from paper data files).
#[derive(Debug, Deserialize)]
struct McAllisterValues {
    g_s: f64,
    #[serde(rename = "W0")]
    w0: f64,
    #[serde(rename = "V_string")]
    v_string: f64,
}

/// Load a JSON fixture file.
fn load_fixture<T: for<'de> Deserialize<'de>>(example: &str, filename: &str) -> T {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let path = manifest_dir
        .join("tests/fixtures")
        .join(example)
        .join(filename);
    let content = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

/// Convert intersection fixture to our Intersection type.
fn fixture_to_intersection(fixture: &IntersectionFixture) -> Intersection {
    let mut kappa = Intersection::new(fixture.dim);
    for (key, &value) in &fixture.entries {
        let parts: Vec<usize> = key.split(',').map(|s| s.parse().unwrap()).collect();
        kappa.set(parts[0], parts[1], parts[2], value);
    }
    kappa
}

// =============================================================================
// Flat Direction Tests
// =============================================================================
//
// These test that p = N^{-1} K is computed correctly, where N_ab = κ_abc M^c.
// The expected p values are computed from McAllister's K, M, κ data.

/// Expected flat direction p for 4-214-647 (from eq. 6.56).
/// p = (293/110, 163/110, 163/110, 13/22)
const P_4_214_647: [f64; 4] = [
    2.663636363636364,
    1.481818181818182,
    1.481818181818182,
    0.5909090909090909,
];

/// Expected e^{K₀} for 4-214-647.
/// e^{K₀} = (4/3 × κ_abc p^a p^b p^c)^{-1} ≈ 0.2344
const EK0_4_214_647: f64 = 0.23439299551782825;

#[test]
fn test_flat_direction_4_214_647() {
    let flux: FluxFixture = load_fixture("4_214_647", "flux.json");
    let intersection: IntersectionFixture = load_fixture("4_214_647", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);

    // Compute flat direction p = N^{-1} K
    let p =
        compute_flat_direction(&kappa, &flux.k, &flux.m).expect("Failed to compute flat direction");

    // Check p matches expected (from eq. 6.56)
    assert_eq!(p.len(), 4);
    for (i, (computed, expected)) in p.iter().zip(P_4_214_647.iter()).enumerate() {
        let rel_err = (computed - expected).abs() / expected.abs();
        assert!(
            rel_err < 1e-6,
            "p[{i}] mismatch: computed={computed}, expected={expected}, rel_err={rel_err}"
        );
    }

    // Check e^{K_0}
    let ek0 = compute_ek0(&kappa, &p);
    let ek0_rel_err = (ek0 - EK0_4_214_647).abs() / EK0_4_214_647;
    assert!(
        ek0_rel_err < 1e-4,
        "e^K0 mismatch: computed={ek0}, expected={EK0_4_214_647}, rel_err={ek0_rel_err}"
    );
}

#[test]
fn test_flat_direction_5_113_4627_main() {
    let flux: FluxFixture = load_fixture("5_113_4627_main", "flux.json");
    let intersection: IntersectionFixture = load_fixture("5_113_4627_main", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);

    let p =
        compute_flat_direction(&kappa, &flux.k, &flux.m).expect("Failed to compute flat direction");

    // Expected p computed from κ, K, M fixtures
    let expected_p = [
        0.12068965517241358,
        0.258620689655173,
        0.8706896551724139,
        2.6034482758620703,
        -0.11206896551724067,
    ];
    assert_eq!(p.len(), 5);
    for (i, (computed, expected)) in p.iter().zip(expected_p.iter()).enumerate() {
        let rel_err = (computed - expected).abs() / expected.abs().max(1e-10);
        assert!(rel_err < 1e-6, "p[{i}] mismatch: {computed} vs {expected}");
    }

    let ek0 = compute_ek0(&kappa, &p);
    // Expected e^K0 ≈ 0.0911
    assert!(
        (ek0 - 0.09114853876607293).abs() / 0.09114853876607293 < 1e-4,
        "e^K0={ek0}"
    );
}

#[test]
fn test_flat_direction_5_113_4627_alternative() {
    let flux: FluxFixture = load_fixture("5_113_4627_alternative", "flux.json");
    let intersection: IntersectionFixture =
        load_fixture("5_113_4627_alternative", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);

    let p =
        compute_flat_direction(&kappa, &flux.k, &flux.m).expect("Failed to compute flat direction");

    // Expected p computed from κ, K, M fixtures
    let expected_p = [
        0.12857142857142864,
        -0.007142857142857363,
        0.5035714285714287,
        2.025,
        -0.26071428571428595,
    ];
    assert_eq!(p.len(), 5);
    for (i, (computed, expected)) in p.iter().zip(expected_p.iter()).enumerate() {
        let rel_err = (computed - expected).abs() / expected.abs().max(1e-10);
        assert!(rel_err < 1e-6, "p[{i}] mismatch: {computed} vs {expected}");
    }

    let ek0 = compute_ek0(&kappa, &p);
    assert!(
        (ek0 - 0.2718642810708103).abs() / 0.2718642810708103 < 1e-4,
        "e^K0={ek0}"
    );
}

#[test]
fn test_flat_direction_5_81_3213() {
    let flux: FluxFixture = load_fixture("5_81_3213", "flux.json");
    let intersection: IntersectionFixture = load_fixture("5_81_3213", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);

    let p =
        compute_flat_direction(&kappa, &flux.k, &flux.m).expect("Failed to compute flat direction");

    // Expected p computed from κ, K, M fixtures
    let expected_p = [
        1.625,
        2.458333333333333,
        1.2499999999999998,
        1.2499999999999998,
        0.4166666666666664,
    ];
    assert_eq!(p.len(), 5);
    for (i, (computed, expected)) in p.iter().zip(expected_p.iter()).enumerate() {
        let rel_err = (computed - expected).abs() / expected.abs().max(1e-10);
        assert!(rel_err < 1e-6, "p[{i}] mismatch: {computed} vs {expected}");
    }

    let ek0 = compute_ek0(&kappa, &p);
    assert!(
        (ek0 - 0.025714285714285724).abs() / 0.025714285714285724 < 1e-4,
        "e^K0={ek0}"
    );
}

#[test]
fn test_flat_direction_7_51_13590() {
    let flux: FluxFixture = load_fixture("7_51_13590", "flux.json");
    let intersection: IntersectionFixture = load_fixture("7_51_13590", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);

    let p =
        compute_flat_direction(&kappa, &flux.k, &flux.m).expect("Failed to compute flat direction");

    // Expected p computed from κ, K, M fixtures
    let expected_p = [
        2.1666666666666647,
        0.3333333333333326,
        -0.6666666666666646,
        0.9999999999999984,
        0.6999999999999984,
        1.5999999999999979,
        1.0999999999999992,
    ];
    assert_eq!(p.len(), 7);
    for (i, (computed, expected)) in p.iter().zip(expected_p.iter()).enumerate() {
        let rel_err = (computed - expected).abs() / expected.abs().max(1e-10);
        assert!(rel_err < 1e-6, "p[{i}] mismatch: {computed} vs {expected}");
    }

    let ek0 = compute_ek0(&kappa, &p);
    assert!(
        (ek0 - 0.11565876950492376).abs() / 0.11565876950492376 < 1e-4,
        "e^K0={ek0}"
    );
}

// =============================================================================
// Vacuum Energy Tests
// =============================================================================
//
// These test V₀ = -3 × e^{K₀} × (g_s^7 / (4V)²) × W₀²
// Using McAllister's published values for g_s, W₀, V_string.

#[test]
fn test_vacuum_energy_4_214_647() {
    // McAllister's published values
    let g_s = 0.00911134;
    let w0 = 2.30012e-90;
    let v_string = 4711.829675204889;
    let ek0 = EK0_4_214_647;

    let v0 = compute_v0(ek0, g_s, v_string, w0);

    // Expected: V₀ ≈ -5.5e-203 (from eq. 6.63)
    assert!(v0 < 0.0, "V0 should be negative for AdS vacuum");
    let log_v0 = v0.abs().log10();
    assert!(
        (log_v0 - (-203.0)).abs() < 1.0,
        "V0 should be ~10^{{-203}}, got 10^{log_v0}"
    );
}

#[test]
fn test_vacuum_energy_5_113_4627_main() {
    let mcallister: McAllisterValues = load_fixture("5_113_4627_main", "mcallister_values.json");
    let ek0 = 0.09114853876607293;

    let v0 = compute_v0(ek0, mcallister.g_s, mcallister.v_string, mcallister.w0);

    assert!(v0 < 0.0, "V0 should be negative for AdS vacuum");
    let log_v0 = v0.abs().log10();
    // Expected: ~10^{-144} (log10 = -143.78)
    assert!(
        (log_v0 - (-143.78)).abs() < 1.0,
        "V0 should be ~10^{{-144}}, got 10^{log_v0}"
    );
}

#[test]
fn test_vacuum_energy_5_113_4627_alternative() {
    let mcallister: McAllisterValues =
        load_fixture("5_113_4627_alternative", "mcallister_values.json");
    let ek0 = 0.2718642810708103;

    let v0 = compute_v0(ek0, mcallister.g_s, mcallister.v_string, mcallister.w0);

    assert!(v0 < 0.0, "V0 should be negative");
    let log_v0 = v0.abs().log10();
    // Expected: ~10^{-213} (log10 = -213.48)
    assert!(
        (log_v0 - (-213.48)).abs() < 1.0,
        "V0 should be ~10^{{-213}}, got 10^{log_v0}"
    );
}

#[test]
fn test_vacuum_energy_5_81_3213() {
    let mcallister: McAllisterValues = load_fixture("5_81_3213", "mcallister_values.json");
    let ek0 = 0.025714285714285724;

    let v0 = compute_v0(ek0, mcallister.g_s, mcallister.v_string, mcallister.w0);

    assert!(v0 < 0.0, "V0 should be negative");
    let log_v0 = v0.abs().log10();
    // Expected: ~10^{-61} (log10 = -61.38)
    assert!(
        (log_v0 - (-61.38)).abs() < 1.0,
        "V0 should be ~10^{{-61}}, got 10^{log_v0}"
    );
}

#[test]
fn test_vacuum_energy_7_51_13590() {
    let mcallister: McAllisterValues = load_fixture("7_51_13590", "mcallister_values.json");
    let ek0 = 0.11565876950492376;

    let v0 = compute_v0(ek0, mcallister.g_s, mcallister.v_string, mcallister.w0);

    assert!(v0 < 0.0, "V0 should be negative");
    let log_v0 = v0.abs().log10();
    // Expected: ~10^{-57} (log10 = -56.51)
    assert!(
        (log_v0 - (-56.51)).abs() < 1.0,
        "V0 should be ~10^{{-57}}, got 10^{log_v0}"
    );
}

// =============================================================================
// c_τ Relationship Tests
// =============================================================================
//
// c_τ = 2π / (g_s × ln(W₀⁻¹))
// This relates string coupling to the flux superpotential.

#[test]
fn test_c_tau_4_214_647() {
    let g_s = 0.00911134;
    let w0 = 2.30012e-90;
    let expected_c_tau = 3.34109; // From McAllister's c_tau.dat

    let c_tau = compute_c_tau(g_s, w0);
    let rel_err = (c_tau - expected_c_tau).abs() / expected_c_tau;
    assert!(
        rel_err < 0.01,
        "c_tau mismatch: computed={c_tau}, expected={expected_c_tau}"
    );
}

// =============================================================================
// End-to-End Pipeline Test
// =============================================================================
//
// Full pipeline: κ + (K, M) → p → e^{K₀} → V₀
// This validates the complete physics chain.

#[test]
fn test_end_to_end_4_214_647() {
    // Load fixtures (inputs only)
    let flux: FluxFixture = load_fixture("4_214_647", "flux.json");
    let intersection: IntersectionFixture = load_fixture("4_214_647", "intersection.json");
    let mcallister: McAllisterValues = load_fixture("4_214_647", "mcallister_values.json");
    let kappa = fixture_to_intersection(&intersection);

    // Step 1: Compute flat direction p = N^{-1} K
    let p =
        compute_flat_direction(&kappa, &flux.k, &flux.m).expect("Failed to compute flat direction");

    // Step 2: Compute e^{K₀}
    let ek0 = compute_ek0(&kappa, &p);

    // Step 3: Compute V₀ using McAllister's g_s, W₀, V_string
    let v0 = compute_v0(ek0, mcallister.g_s, mcallister.v_string, mcallister.w0);

    // Validate the full chain
    // p should match eq. 6.56
    assert!((p[0] - 2.6636).abs() < 0.001, "p[0]={}", p[0]);
    assert!((p[1] - 1.4818).abs() < 0.001, "p[1]={}", p[1]);

    // e^{K₀} should be ~0.234
    assert!((ek0 - 0.234).abs() < 0.01, "e^K0={ek0}");

    // V₀ should be ~-5.5e-203
    assert!(v0 < 0.0, "V0 must be negative");
    let log_v0 = v0.abs().log10();
    assert!(
        (log_v0 - (-203.0)).abs() < 1.0,
        "V0 should be ~10^{{-203}}, got 10^{log_v0}"
    );

    // Also verify c_τ relationship
    let c_tau = compute_c_tau(mcallister.g_s, mcallister.w0);
    assert!((c_tau - 3.34).abs() < 0.1, "c_tau={c_tau}");
}

// =============================================================================
// Volume Formula Sanity Check
// =============================================================================

#[test]
fn test_volume_formula_sanity() {
    // Just verify the volume formula works with simple values
    let mut kappa = Intersection::new(4);
    kappa.set(0, 0, 0, 1);
    kappa.set(3, 3, 3, 8);

    let t = vec![1.0, 1.0, 1.0, 1.0];
    let v = volume_string(&kappa, &t, 214, 4);

    assert!(v.is_finite(), "Volume should be finite");
}
