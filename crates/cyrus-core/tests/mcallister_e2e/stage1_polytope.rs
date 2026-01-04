//! McAllister Pipeline Stage 1: Polytope Points
//!
//! The raw primitive input - lattice points defining the reflexive polytope.
//! Everything else is COMPUTED from this.
//!
//! Polytope 4-214-647 from Kreuzer-Skarke database:
//! - Primal: 294 lattice points, h11 = 214
//! - Dual (mirror): 12 lattice points, h21 = 4
//!
//! Reference: arXiv:2107.09064

#![allow(missing_docs)]

use serde::Deserialize;
use std::collections::HashSet;
use std::path::PathBuf;

use cyrus_core::{Point, Polytope};

/// Input format: only primal points (we compute the dual)
#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

/// Expected dual points for assertion
#[derive(Debug, Deserialize)]
struct DualAssertion {
    points: Vec<Vec<i64>>,
}

/// Loaded fixture with computed dual
struct PolytopeFixture {
    primal_points: Vec<Vec<i64>>,
    dual_points: Vec<Vec<i64>>,
}

fn load_polytope() -> PolytopeFixture {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    // Load primal points from inputs
    let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
    let content = std::fs::read_to_string(&input_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
    let input: PolytopeInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()));

    // Compute the dual polytope
    let primal_verts: Vec<Point> = input
        .points
        .iter()
        .filter(|p| !p.iter().all(|&x| x == 0)) // Exclude origin for vertex set
        .map(|coords| Point::new(coords.clone()))
        .collect();

    let polytope = Polytope::from_vertices(primal_verts)
        .expect("Failed to create polytope from primal vertices");

    let dual = polytope
        .compute_dual()
        .expect("Failed to compute dual polytope");

    // Dual now includes all lattice points (including origin)
    let dual_points: Vec<Vec<i64>> = dual.vertices().iter().map(|p| p.coords().to_vec()).collect();

    PolytopeFixture {
        primal_points: input.points,
        dual_points,
    }
}

fn to_points(raw: &[Vec<i64>]) -> Vec<Point> {
    raw.iter()
        .map(|coords| Point::new(coords.clone()))
        .collect()
}

#[test]
fn stage1_primal_point_count() {
    let fixture = load_polytope();
    insta::assert_snapshot!(
        "primal_point_count",
        format!("{}", fixture.primal_points.len())
    );
}

#[test]
fn stage1_dual_point_count() {
    let fixture = load_polytope();
    insta::assert_snapshot!(
        "dual_point_count",
        format!("{}", fixture.dual_points.len())
    );
}

#[test]
fn stage1_primal_dimension() {
    let fixture = load_polytope();
    let points = to_points(&fixture.primal_points);

    // All points should be 4-dimensional (CY3 in 4D ambient)
    assert!(!points.is_empty(), "No primal points");
    let dim = points[0].dim();

    insta::assert_snapshot!("primal_dimension", format!("{}", dim));

    for (i, p) in points.iter().enumerate() {
        assert_eq!(
            p.dim(),
            dim,
            "Point {} has dimension {} but expected {}",
            i,
            p.dim(),
            dim
        );
    }
}

#[test]
fn stage1_primal_contains_origin() {
    let fixture = load_polytope();
    let points = to_points(&fixture.primal_points);

    let origin_count = points
        .iter()
        .filter(|p| p.coords().iter().all(|&x| x == 0))
        .count();

    // Reflexive polytope must contain exactly one interior point (the origin)
    assert_eq!(
        origin_count, 1,
        "Reflexive polytope must have exactly one origin, found {}",
        origin_count
    );
}

#[test]
fn stage1_dual_contains_origin() {
    let fixture = load_polytope();
    let points = to_points(&fixture.dual_points);

    let origin_count = points
        .iter()
        .filter(|p| p.coords().iter().all(|&x| x == 0))
        .count();

    assert_eq!(
        origin_count, 1,
        "Dual polytope must have exactly one origin, found {}",
        origin_count
    );
}

#[test]
fn stage1_primal_sample_points() {
    let fixture = load_polytope();

    // Snapshot first few points for verification
    let mut output = String::new();
    for (i, pt) in fixture.primal_points.iter().take(10).enumerate() {
        output.push_str(&format!("p[{}] = {:?}\n", i, pt));
    }

    insta::assert_snapshot!("primal_sample_points", output);
}

#[test]
fn stage1_dual_sample_points() {
    let fixture = load_polytope();

    let mut output = String::new();
    for (i, pt) in fixture.dual_points.iter().enumerate() {
        output.push_str(&format!("d[{}] = {:?}\n", i, pt));
    }

    insta::assert_snapshot!("dual_sample_points", output);
}

/// Assert that our computed dual matches McAllister's expected dual points.
/// This is an EXACT assertion - the sets must match exactly.
#[test]
fn stage1_dual_matches_expected() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    // Load expected dual from assertions
    let assertion_path = manifest_dir.join("tests/mcallister_e2e/assertions/dual_points.json");
    let content = std::fs::read_to_string(&assertion_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", assertion_path.display()));
    let expected: DualAssertion = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", assertion_path.display()));

    // Load computed dual
    let fixture = load_polytope();

    // Compare as sets (order doesn't matter)
    let computed_set: HashSet<Vec<i64>> = fixture.dual_points.into_iter().collect();
    let expected_set: HashSet<Vec<i64>> = expected.points.into_iter().collect();

    assert_eq!(
        computed_set.len(),
        expected_set.len(),
        "Dual point count mismatch: computed {} vs expected {}",
        computed_set.len(),
        expected_set.len()
    );

    // Check for missing points
    let missing: Vec<_> = expected_set.difference(&computed_set).collect();
    let extra: Vec<_> = computed_set.difference(&expected_set).collect();

    assert!(
        missing.is_empty() && extra.is_empty(),
        "Dual polytope mismatch!\nMissing (in expected but not computed): {:?}\nExtra (in computed but not expected): {:?}",
        missing,
        extra
    );
}
