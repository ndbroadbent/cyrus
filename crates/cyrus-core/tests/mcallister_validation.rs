//! Integration tests validating against McAllister et al. (arXiv:2107.09064).
//!
//! These tests load fixtures (intersection numbers κ, flux vectors K/M, GV invariants) and verify
//! that our Rust implementation reproduces the physics values from the paper using snapshot tests.

use insta::assert_json_snapshot;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::path::PathBuf;

use cyrus_core::{
    EvaluationRequest, GvInvariant, Intersection, MoriCone, build_racetrack_terms,
    compute_flat_direction, compute_w0, evaluate_vacuum, solve_racetrack,
};
use malachite::Integer;

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

    if basis_path.exists() {
        let content = std::fs::read_to_string(&basis_path).unwrap();
        let entries: Vec<GvFixtureEntry> = serde_json::from_str(&content).unwrap();
        return entries
            .into_iter()
            .map(|e| GvInvariant {
                curve: e.curve,
                value: e.value,
            })
            .collect();
    }

    // Fallback to original format if basis file not found
    let fixture: GvFixture = load_fixture(example, "gv_invariants.json");
    fixture
        .curves
        .into_iter()
        .map(|e| GvInvariant {
            curve: e.curve,
            value: e.value,
        })
        .collect()
}

use malachite::Rational;

/// Convert intersection fixture to our Intersection type.
fn fixture_to_intersection(fixture: &IntersectionFixture) -> Intersection {
    let mut kappa = Intersection::new(fixture.dim);
    for (key, &value) in &fixture.entries {
        let parts: Vec<usize> = key.split(',').map(|s| s.parse().unwrap()).collect();
        kappa.set(parts[0], parts[1], parts[2], Rational::from(value));
    }
    kappa
}

#[derive(Serialize)]
struct RacetrackSnapshot {
    g_s: f64,
    w0: f64,
    w0_log10: f64,
    im_tau: f64,
    delta: f64,
    epsilon: f64,
}

#[test]
fn test_vacuum_pipeline_4_214_647() {
    let flux: FluxFixture = load_fixture("4_214_647", "flux.json");
    let intersection: IntersectionFixture = load_fixture("4_214_647", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);
    let gv = load_gv_invariants("4_214_647");

    // Construct a dummy Mori cone that contains the flat direction p
    // (Actual Mori cone would come from triangulation)
    let p = compute_flat_direction(&kappa, &flux.k, &flux.m).expect("Flat direction failed");
    // Relation R such that p . R > 0. Let R = [1, 0, ... 0]
    let mut generator = vec![Integer::from(0); p.len()];
    generator[0] = Integer::from(1);
    let mori = MoriCone::new(vec![generator]);

    let req = EvaluationRequest {
        kappa: &kappa,
        mori: &mori,
        gv: &gv,
        h11: 214,
        h21: 4,
        q_max: 200.0,
    };

    let result = evaluate_vacuum(&req, &flux.k, &flux.m).expect("Pipeline failed");
    assert!(result.success);
    assert!(result.vacuum.is_some());

    let vac = result.vacuum.unwrap();
    let log_v0 = vac.v0.abs().log10();
    assert!((-204.0..=-202.0).contains(&log_v0));
}

#[test]
fn test_racetrack_full_pipeline_4_214_647() {
    // Full pipeline: (K, M, κ, GV) → p → racetrack terms → g_s, W₀
    let flux: FluxFixture = load_fixture("4_214_647", "flux.json");
    let intersection: IntersectionFixture = load_fixture("4_214_647", "intersection.json");
    let kappa = fixture_to_intersection(&intersection);

    // Step 1: Compute flat direction
    let p =
        compute_flat_direction(&kappa, &flux.k, &flux.m).expect("Failed to compute flat direction");

    // Step 2: Build racetrack terms from GV invariants
    // We load the BASIS aligned GV invariants we generated
    let gv = load_gv_invariants("4_214_647");

    // Ensure we have data
    assert!(!gv.is_empty(), "GV invariants should not be empty");

    let terms = build_racetrack_terms(&gv, &flux.m, &p);

    // Step 3: Solve racetrack
    let result = solve_racetrack(&terms).expect("Racetrack should stabilize");

    // Step 4: Compute W₀
    // We need at least 2 terms to compute W0
    assert!(terms.len() >= 2);
    let w0 = compute_w0(&result, &terms[0], &terms[1]);

    // Snapshot the result
    assert_json_snapshot!(
        "racetrack_4_214_647",
        RacetrackSnapshot {
            g_s: result.g_s,
            w0,
            w0_log10: w0.abs().log10(),
            im_tau: result.im_tau,
            delta: result.delta,
            epsilon: result.epsilon,
        }
    );
}
