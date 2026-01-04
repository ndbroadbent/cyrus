//! McAllister Pipeline Stage 2: Triangulation + Heights
//!
//! Computes regular triangulation from polytope points + heights.
//!
//! - "ours" branch: Our computed FRST heights (|p|² with star adjustment)
//! - "theirs" branch: McAllister's heights from inputs/heights.json

#![allow(missing_docs)]

use serde::Deserialize;
use std::path::PathBuf;

use cyrus_core::{Point, Polytope, compute_frst_heights, compute_regular_triangulation};

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

struct Stage2Fixture {
    /// Points not interior to facets (used for triangulation)
    triangulation_points: Vec<Point>,
    /// Index of origin in triangulation_points
    origin_idx: usize,
}

fn load_stage2_fixture() -> Stage2Fixture {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    // Load primal points from inputs
    let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
    let content = std::fs::read_to_string(&input_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
    let input: PolytopeInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()));

    // Create polytope from all primal points
    let all_points: Vec<Point> = input
        .points
        .iter()
        .map(|coords| Point::new(coords.clone()))
        .collect();

    let polytope = Polytope::from_vertices(all_points)
        .expect("Failed to create polytope");

    // Filter to points not interior to facets
    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter points");

    // Find origin index
    let origin_idx = triangulation_points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .expect("Origin not found in triangulation points");

    Stage2Fixture {
        triangulation_points,
        origin_idx,
    }
}

// =============================================================================
// Common assertions (apply to both branches)
// =============================================================================

#[test]
fn stage2_triangulation_point_count() {
    let fixture = load_stage2_fixture();
    assert_eq!(
        fixture.triangulation_points.len(),
        219,
        "Expected 219 triangulation points (points not interior to facets)"
    );
}

#[test]
fn stage2_origin_index() {
    let fixture = load_stage2_fixture();
    assert_eq!(fixture.origin_idx, 0, "Origin should be at index 0");
}

// =============================================================================
// "ours" branch: Our computed FRST heights (|p|² with star adjustment)
// =============================================================================

#[test]
#[cfg(feature = "slow-tests")]
fn stage2_ours_frst_heights() {
    let fixture = load_stage2_fixture();

    let (heights, triangulation) = compute_frst_heights(
        &fixture.triangulation_points,
        fixture.origin_idx,
    ).expect("Failed to compute FRST heights");

    // Verify star property
    assert!(
        triangulation.is_star(fixture.origin_idx),
        "Triangulation must have star property"
    );

    insta::assert_json_snapshot!("heights_ours", heights);
}

#[test]
#[cfg(feature = "slow-tests")]
fn stage2_ours_triangulation_simplex_count() {
    let fixture = load_stage2_fixture();

    let (_heights, triangulation) = compute_frst_heights(
        &fixture.triangulation_points,
        fixture.origin_idx,
    ).expect("Failed to compute FRST heights");

    assert_eq!(
        triangulation.simplices().len(),
        1107,
        "Expected 1107 simplices in FRST triangulation"
    );
}

#[test]
#[cfg(feature = "slow-tests")]
fn stage2_ours_triangulation_simplices() {
    let fixture = load_stage2_fixture();

    let (_heights, triangulation) = compute_frst_heights(
        &fixture.triangulation_points,
        fixture.origin_idx,
    ).expect("Failed to compute FRST heights");

    // Sort simplices for deterministic snapshot
    let mut simplices: Vec<_> = triangulation.simplices().to_vec();
    simplices.sort();

    insta::assert_json_snapshot!("triangulation_simplices_ours", simplices);
}

#[test]
#[cfg(feature = "slow-tests")]
fn stage2_ours_all_simplices_contain_origin() {
    let fixture = load_stage2_fixture();

    let (_heights, triangulation) = compute_frst_heights(
        &fixture.triangulation_points,
        fixture.origin_idx,
    ).expect("Failed to compute FRST heights");

    // Verify star property explicitly
    for (i, simplex) in triangulation.simplices().iter().enumerate() {
        assert!(
            simplex.contains(&fixture.origin_idx),
            "Simplex {} does not contain origin (idx {}): {:?}",
            i, fixture.origin_idx, simplex
        );
    }
}

// =============================================================================
// "theirs" branch: McAllister's heights from inputs/heights.json
// =============================================================================

#[derive(Debug, Deserialize)]
struct HeightsInput {
    values: Vec<f64>,
}

fn load_mcallister_heights() -> Vec<f64> {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let path = manifest_dir.join("tests/mcallister_e2e/inputs/heights.json");
    let content = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    let input: HeightsInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", path.display()));
    input.values
}

#[test]
fn stage2_theirs_heights_count() {
    let heights = load_mcallister_heights();
    assert_eq!(heights.len(), 219, "McAllister heights should have 219 values");
}

#[test]
#[cfg(feature = "slow-tests")]
fn stage2_theirs_triangulation() {
    let fixture = load_stage2_fixture();
    let heights = load_mcallister_heights();

    let triangulation = compute_regular_triangulation(
        &fixture.triangulation_points,
        &heights,
    ).expect("Failed to compute triangulation from McAllister heights");

    // Sort simplices for deterministic snapshot
    let mut simplices: Vec<_> = triangulation.simplices().to_vec();
    simplices.sort();

    insta::assert_json_snapshot!("triangulation_simplices_theirs", simplices);
}

#[test]
#[cfg(feature = "slow-tests")]
fn stage2_theirs_is_star() {
    let fixture = load_stage2_fixture();
    let heights = load_mcallister_heights();

    let triangulation = compute_regular_triangulation(
        &fixture.triangulation_points,
        &heights,
    ).expect("Failed to compute triangulation from McAllister heights");

    assert!(
        triangulation.is_star(fixture.origin_idx),
        "McAllister's triangulation must have star property"
    );
}
