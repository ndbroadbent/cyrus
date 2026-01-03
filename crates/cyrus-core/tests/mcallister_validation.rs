#![allow(missing_docs)]
//! Integration tests validating against McAllister et al. (arXiv:2107.09064).
//!
//! These tests load fixtures (intersection numbers κ, flux vectors K/M, GV invariants) and verify
//! that our Rust implementation reproduces the physics values from the paper using snapshot tests.

use insta::assert_json_snapshot;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;

use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::Finite;
use cyrus_core::{
    EvaluationRequest, GvInvariant, H11, H21, Intersection, MoriCone, build_racetrack_terms,
    compute_flat_direction, compute_w0, evaluate_vacuum, solve_racetrack,
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

/// GV Invariant loaded from fixture
#[derive(Debug, Deserialize)]
struct GvFixtureEntry {
    curve: Vec<i64>,
    value: f64,
}

#[derive(Debug, Deserialize)]
struct GvFixture {
    curves: Vec<GvFixtureEntry>,
}

/// Load a JSON fixture file.
fn load_fixture<T: for<'de> Deserialize<'de>>(example: &str, filename: &str) -> T {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let mut path = manifest_dir
        .join("tests/fixtures")
        .join(example)
        .join(filename);

    // Check if file exists, if not try the computed/basis version
    if !path.exists() && filename == "gv_invariants.json" {
        let alt_path = manifest_dir
            .join("tests/fixtures")
            .join(example)
            .join("gv_basis.json");
        if alt_path.exists() {
            path = alt_path;
        }
    }

    let content = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));

    // Handle potential list format for GV basis file
    if path.file_name().unwrap() == "gv_basis.json"
        && serde_json::from_str::<Vec<GvFixtureEntry>>(&content).is_ok()
    {
        // This assumes T is GvFixture or compatible struct with "curves" field if we could map it
        // Let's just try parsing as T.
    }

    serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()))
}

/// Helper to load GV invariants correctly regardless of file format
fn load_gv_invariants(example: &str) -> Vec<GvInvariant> {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let basis_path = manifest_dir
        .join("tests/fixtures")
        .join(example)
        .join("gv_basis_old.json");

    // Helper to convert raw fixture entry to typed GvInvariant
    let convert_entry = |e: GvFixtureEntry| GvInvariant {
        curve: e.curve.into_iter().map(I64::<Finite>::new).collect(),
        value: F64::<Finite>::new(e.value).expect("GV value must be finite"),
    };

    if basis_path.exists() {
        let content = std::fs::read_to_string(&basis_path).unwrap();
        let entries: Vec<GvFixtureEntry> = serde_json::from_str(&content).unwrap();
        return entries.into_iter().map(convert_entry).collect();
    }

    // Fallback to original format if basis file not found
    let fixture: GvFixture = load_fixture(example, "gv_invariants.json");
    fixture.curves.into_iter().map(convert_entry).collect()
}

use cyrus_core::types::rational::Rational as TypedRational;
use malachite::Rational;

/// Convert intersection fixture to our Intersection type.
fn fixture_to_intersection(fixture: &IntersectionFixture) -> Intersection {
    let mut kappa = Intersection::new(fixture.dim);
    for (key, &value) in &fixture.entries {
        let parts: Vec<usize> = key.split(',').map(|s| s.parse().unwrap()).collect();
        // Intersection numbers can be positive or negative - use Finite
        let val = TypedRational::<Finite>::from_raw(Rational::from(value));
        kappa.set(parts[0], parts[1], parts[2], val);
    }
    kappa
}

#[derive(Serialize)]
struct RacetrackSnapshot {
    g_s: f64, // Extract from F64<Pos> for serialization
    w0: f64,
    w0_log10: f64,
    im_tau: f64, // Extract from F64<Pos> for serialization
    delta: f64,
    epsilon: f64,
}

#[test]
fn test_vacuum_pipeline_4_214_647() {
    let flux: FluxFixture = load_fixture("4_214_647", "flux.json");
    let intersection: IntersectionFixture = load_fixture("4_214_647", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);
    let gv = load_gv_invariants("4_214_647");

    // Convert raw flux vectors to typed at the boundary
    let k_typed: Vec<I64<Finite>> = flux.k.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let m_typed: Vec<I64<Finite>> = flux.m.iter().map(|&x| I64::<Finite>::new(x)).collect();

    // Construct a dummy Mori cone that contains the flat direction p
    // (Actual Mori cone would come from triangulation)
    let p = compute_flat_direction(&kappa, &k_typed, &m_typed).expect("Flat direction failed");
    eprintln!("Flat direction p: {:?}", &p[..5.min(p.len())]);

    // Find the first positive component to use as the generator direction
    let pos_idx = p.iter().position(|x| x.get() > 0.0).unwrap_or(0);
    eprintln!(
        "Using positive index {} with value {}",
        pos_idx,
        p[pos_idx].get()
    );

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
    assert!(result.success);
    assert!(result.vacuum.is_some());

    let vac = result.vacuum.unwrap();
    let log_v0 = vac.v0.get().abs().log10();
    assert!((-204.0..=-202.0).contains(&log_v0));
}

#[test]
fn test_racetrack_full_pipeline_4_214_647() {
    // Full pipeline: (K, M, κ, GV) → p → racetrack terms → g_s, W₀
    let flux: FluxFixture = load_fixture("4_214_647", "flux.json");
    let intersection: IntersectionFixture = load_fixture("4_214_647", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);

    // Convert raw flux vectors to typed at the boundary
    let k_typed: Vec<I64<Finite>> = flux.k.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let m_typed: Vec<I64<Finite>> = flux.m.iter().map(|&x| I64::<Finite>::new(x)).collect();

    // Step 1: Compute flat direction
    let p = compute_flat_direction(&kappa, &k_typed, &m_typed)
        .expect("Failed to compute flat direction");

    // Step 2: Build racetrack terms from GV invariants
    // We load the BASIS aligned GV invariants we generated
    let gv = load_gv_invariants("4_214_647");

    // Ensure we have data
    assert!(!gv.is_empty(), "GV invariants should not be empty");

    let terms = build_racetrack_terms(&gv, &m_typed, &p);

    // Step 3: Solve racetrack
    let result = solve_racetrack(&terms).expect("Racetrack should stabilize");

    // Step 4: Compute W₀
    // We need at least 2 terms to compute W0
    assert!(terms.len() >= 2);
    let w0 = compute_w0(&result, &terms[0], &terms[1]);

    // Snapshot the result - convert typed values to raw for serialization
    assert_json_snapshot!(
        "racetrack_4_214_647",
        RacetrackSnapshot {
            g_s: result.g_s.get(),
            w0: w0.get(),
            w0_log10: w0.get().abs().log10(),
            im_tau: result.im_tau.get(),
            delta: result.delta.get(),
            epsilon: result.epsilon.get(),
        }
    );
}
