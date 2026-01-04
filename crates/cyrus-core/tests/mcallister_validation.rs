#![allow(missing_docs)]
//! Integration tests against McAllister et al. (arXiv:2107.09064).
//!
//! These tests compare our computed values against published reference data.

use serde::Deserialize;
use std::collections::HashMap;
use std::path::PathBuf;

use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::Finite;
use cyrus_core::{
    EvaluationRequest, GvInvariant, H11, H21, Intersection, MoriCone, build_racetrack_terms,
    compute_flat_direction, compute_w0, evaluate_vacuum, solve_racetrack,
};

// ============================================================================
// Fixture Types
// ============================================================================

/// Flux vectors from fixture.
#[derive(Debug, Deserialize)]
struct FluxFixture {
    #[serde(rename = "K")]
    k: Vec<i64>,
    #[serde(rename = "M")]
    m: Vec<i64>,
}

/// Intersection numbers from fixture.
#[derive(Debug, Deserialize)]
struct IntersectionFixture {
    dim: usize,
    entries: HashMap<String, i64>,
}

/// GV Invariant entry from fixture.
#[derive(Debug, Deserialize)]
struct GvFixtureEntry {
    curve: Vec<i64>,
    value: f64,
}

/// Published reference values from McAllister et al.
///
/// These are the target values we're trying to reproduce.
#[derive(Debug, Deserialize)]
struct McAllisterValues {
    /// String coupling g_s
    g_s: f64,
    /// Superpotential W₀
    #[serde(rename = "W0")]
    w0: f64,
    /// String frame volume V_string
    #[serde(rename = "V_string")]
    v_string: f64,
}

// ============================================================================
// Fixture Loaders
// ============================================================================

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

/// Load GV invariants from fixture.
fn load_gv_invariants(example: &str) -> Vec<GvInvariant> {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let basis_path = manifest_dir
        .join("tests/fixtures")
        .join(example)
        .join("gv_basis_old.json");

    let convert_entry = |e: GvFixtureEntry| GvInvariant {
        curve: e.curve.into_iter().map(I64::<Finite>::new).collect(),
        value: F64::<Finite>::new(e.value).expect("GV value must be finite"),
    };

    if basis_path.exists() {
        let content = std::fs::read_to_string(&basis_path).unwrap();
        let entries: Vec<GvFixtureEntry> = serde_json::from_str(&content).unwrap();
        return entries.into_iter().map(convert_entry).collect();
    }

    panic!("GV invariants file not found: {}", basis_path.display());
}

use cyrus_core::types::rational::Rational as TypedRational;
use malachite::Rational;

/// Convert intersection fixture to our Intersection type.
fn fixture_to_intersection(fixture: &IntersectionFixture) -> Intersection {
    let mut kappa = Intersection::new(fixture.dim);
    for (key, &value) in &fixture.entries {
        let parts: Vec<usize> = key.split(',').map(|s| s.parse().unwrap()).collect();
        let val = TypedRational::<Finite>::from_raw(Rational::from(value));
        kappa.set(parts[0], parts[1], parts[2], val);
    }
    kappa
}

// ============================================================================
// Validation Tests Against Published Values
// ============================================================================

/// Test g_s against McAllister published value.
#[test]
fn test_gs_against_mcallister_4_214_647() {
    let reference: McAllisterValues = load_fixture("4_214_647", "mcallister_values.json");
    let flux: FluxFixture = load_fixture("4_214_647", "flux.json");
    let intersection: IntersectionFixture = load_fixture("4_214_647", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);

    let k_typed: Vec<I64<Finite>> = flux.k.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let m_typed: Vec<I64<Finite>> = flux.m.iter().map(|&x| I64::<Finite>::new(x)).collect();

    let p = compute_flat_direction(&kappa, &k_typed, &m_typed)
        .expect("Failed to compute flat direction");

    let gv = load_gv_invariants("4_214_647");
    let terms = build_racetrack_terms(&gv, &m_typed, &p);
    let result = solve_racetrack(&terms).expect("Racetrack should stabilize");

    let computed_gs = result.g_s.get();
    let expected_gs = reference.g_s;
    let relative_error = (computed_gs - expected_gs).abs() / expected_gs;

    eprintln!("g_s comparison:");
    eprintln!("  Expected (McAllister): {:.10}", expected_gs);
    eprintln!("  Computed:              {:.10}", computed_gs);
    eprintln!("  Relative error:        {:.2e}", relative_error);

    // g_s should match to 0.1% (this currently passes)
    assert!(
        relative_error < 0.001,
        "g_s mismatch: expected {}, got {} (error: {:.2}%)",
        expected_gs,
        computed_gs,
        relative_error * 100.0
    );
}

/// Test W₀ against McAllister published value.
#[test]
fn test_w0_against_mcallister_4_214_647() {
    let reference: McAllisterValues = load_fixture("4_214_647", "mcallister_values.json");
    let flux: FluxFixture = load_fixture("4_214_647", "flux.json");
    let intersection: IntersectionFixture = load_fixture("4_214_647", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);

    let k_typed: Vec<I64<Finite>> = flux.k.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let m_typed: Vec<I64<Finite>> = flux.m.iter().map(|&x| I64::<Finite>::new(x)).collect();

    let p = compute_flat_direction(&kappa, &k_typed, &m_typed)
        .expect("Failed to compute flat direction");

    let gv = load_gv_invariants("4_214_647");
    let terms = build_racetrack_terms(&gv, &m_typed, &p);
    let result = solve_racetrack(&terms).expect("Racetrack should stabilize");

    assert!(terms.len() >= 2, "Need at least 2 racetrack terms");
    let w0 = compute_w0(&result, &terms[0], &terms[1]);

    let computed_w0 = w0.get().abs();
    let expected_w0 = reference.w0.abs();

    let computed_log = computed_w0.log10();
    let expected_log = expected_w0.log10();
    let log_difference = (computed_log - expected_log).abs();

    eprintln!("W₀ comparison:");
    eprintln!("  Expected (McAllister): {:.6e} (log10: {:.2})", expected_w0, expected_log);
    eprintln!("  Computed:              {:.6e} (log10: {:.2})", computed_w0, computed_log);
    eprintln!("  Log10 difference:      {:.2} orders of magnitude", log_difference);

    // W₀ should match to within 1 order of magnitude
    // Currently fails: computed is ~10^-88, expected is ~10^-90
    assert!(
        log_difference < 1.0,
        "W₀ mismatch: expected {:.2e}, got {:.2e} (off by 10^{:.1})",
        expected_w0,
        computed_w0,
        log_difference
    );
}

/// Diagnostic test that runs without ignoring - shows current discrepancy.
#[test]
fn test_w0_diagnostic_4_214_647() {
    let reference: McAllisterValues = load_fixture("4_214_647", "mcallister_values.json");
    let flux: FluxFixture = load_fixture("4_214_647", "flux.json");
    let intersection: IntersectionFixture = load_fixture("4_214_647", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);

    let k_typed: Vec<I64<Finite>> = flux.k.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let m_typed: Vec<I64<Finite>> = flux.m.iter().map(|&x| I64::<Finite>::new(x)).collect();

    let p = compute_flat_direction(&kappa, &k_typed, &m_typed)
        .expect("Failed to compute flat direction");

    let gv = load_gv_invariants("4_214_647");
    eprintln!("Loaded {} GV invariants", gv.len());

    let terms = build_racetrack_terms(&gv, &m_typed, &p);
    eprintln!("Built {} racetrack terms", terms.len());

    // Show first few terms for debugging
    for (i, term) in terms.iter().take(5).enumerate() {
        eprintln!(
            "  Term {}: coeff={:.6e}, exp={:.6}",
            i,
            term.coefficient.get(),
            term.exponent.get()
        );
    }

    let result = solve_racetrack(&terms).expect("Racetrack should stabilize");
    eprintln!("\nRacetrack solution:");
    eprintln!("  g_s = {:.10}", result.g_s.get());
    eprintln!("  im_tau = {:.6}", result.im_tau.get());

    assert!(terms.len() >= 2, "Need at least 2 racetrack terms");
    let w0 = compute_w0(&result, &terms[0], &terms[1]);

    let computed_w0 = w0.get();
    let expected_w0 = reference.w0;

    eprintln!("\nW₀ comparison:");
    eprintln!("  Expected (McAllister): {:.6e}", expected_w0);
    eprintln!("  Computed:              {:.6e}", computed_w0);
    eprintln!(
        "  Ratio (computed/expected): {:.2e}",
        computed_w0.abs() / expected_w0.abs()
    );

    // This test always passes - it's just for diagnostics
    // The actual validation is in test_w0_against_mcallister_4_214_647
}

/// Test V_string against McAllister published value.
#[test]
fn test_v_string_against_mcallister_4_214_647() {
    let reference: McAllisterValues = load_fixture("4_214_647", "mcallister_values.json");
    let flux: FluxFixture = load_fixture("4_214_647", "flux.json");
    let intersection: IntersectionFixture = load_fixture("4_214_647", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);
    let gv = load_gv_invariants("4_214_647");

    let k_typed: Vec<I64<Finite>> = flux.k.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let m_typed: Vec<I64<Finite>> = flux.m.iter().map(|&x| I64::<Finite>::new(x)).collect();

    let p = compute_flat_direction(&kappa, &k_typed, &m_typed).expect("Flat direction failed");

    // Create Mori cone containing the flat direction
    let pos_idx = p.iter().position(|x| x.get() > 0.0).unwrap_or(0);
    let mut generator: Vec<I64<Finite>> = vec![I64::<Finite>::new(0); p.len()];
    generator[pos_idx] = I64::<Finite>::new(1);
    let mori = MoriCone::new(vec![generator]);

    let req = EvaluationRequest {
        kappa: &kappa,
        mori: &mori,
        gv: &gv,
        h11: H11::new(214).unwrap(),
        h21: H21::new(4).unwrap(),
        q_max: 200.0,
    };

    let result = evaluate_vacuum(&req, &flux.k, &flux.m).expect("Pipeline failed");
    assert!(result.success, "Pipeline should succeed");

    let vac = result.vacuum.unwrap();
    let computed_v_string = vac.v_string.get();
    let expected_v_string = reference.v_string;
    let relative_error = (computed_v_string - expected_v_string).abs() / expected_v_string;

    eprintln!("V_string comparison:");
    eprintln!("  Expected (McAllister): {:.6}", expected_v_string);
    eprintln!("  Computed:              {:.6}", computed_v_string);
    eprintln!("  Relative error:        {:.2e}", relative_error);

    assert!(
        relative_error < 0.01,
        "V_string mismatch: expected {}, got {} (error: {:.2}%)",
        expected_v_string,
        computed_v_string,
        relative_error * 100.0
    );
}

// ============================================================================
// Full Pipeline Tests
// ============================================================================

/// Full pipeline test.
#[test]
fn test_vacuum_pipeline_4_214_647() {
    let flux: FluxFixture = load_fixture("4_214_647", "flux.json");
    let intersection: IntersectionFixture = load_fixture("4_214_647", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);
    let gv = load_gv_invariants("4_214_647");

    let k_typed: Vec<I64<Finite>> = flux.k.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let m_typed: Vec<I64<Finite>> = flux.m.iter().map(|&x| I64::<Finite>::new(x)).collect();

    // Compute flat direction and create a Mori cone that contains it
    let p = compute_flat_direction(&kappa, &k_typed, &m_typed).expect("Flat direction failed");

    // Find the first positive component to use as the generator direction
    // NOTE: This is a simplification - real Mori cone would come from triangulation
    let pos_idx = p.iter().position(|x| x.get() > 0.0).unwrap_or(0);
    let mut generator: Vec<I64<Finite>> = vec![I64::<Finite>::new(0); p.len()];
    generator[pos_idx] = I64::<Finite>::new(1);
    let mori = MoriCone::new(vec![generator]);

    let req = EvaluationRequest {
        kappa: &kappa,
        mori: &mori,
        gv: &gv,
        h11: H11::new(214).unwrap(),
        h21: H21::new(4).unwrap(),
        q_max: 200.0,
    };

    let result = evaluate_vacuum(&req, &flux.k, &flux.m).expect("Pipeline failed");

    if !result.success {
        eprintln!("Vacuum evaluation failed: {:?}", result.reason);
    }
    assert!(result.success, "Pipeline should succeed");
    assert!(result.vacuum.is_some(), "Should produce vacuum result");

    let vac = result.vacuum.unwrap();
    let log_v0 = vac.v0.get().abs().log10();

    // V₀ should be extremely small (cosmological constant scale)
    // Note: This is a loose bound - we need to tighten once W₀ is fixed
    assert!(
        log_v0 < -100.0,
        "V₀ should be very small, got 10^{:.1}",
        log_v0
    );
}
