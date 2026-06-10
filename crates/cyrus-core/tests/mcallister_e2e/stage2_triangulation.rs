//! McAllister Pipeline Stage 2: Triangulation + Heights
//!
//! Computes regular triangulation from polytope points + heights.
//!
//! - "ours" branch: Our computed FRST heights (|p|² with star adjustment)
//! - "theirs" branch: McAllister's heights from data files

#![allow(missing_docs)]

use serde::Deserialize;
use std::path::PathBuf;

use cyrus_core::{Point, Polytope};
use cyrus_core::{compute_frst_heights, compute_regular_triangulation};

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

fn read_csv_rows_i64(path: &PathBuf) -> Vec<Vec<i64>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|s| s.trim().parse::<i64>().expect("invalid integer"))
                .collect::<Vec<_>>()
        })
        .collect()
}

fn read_csv_f64(path: &PathBuf) -> Vec<f64> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .split([',', '\n', '\r'])
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse::<f64>().expect("invalid float"))
        .collect()
}

fn require_first_principles() -> bool {
    if !crate::first_principles_enabled() {
        eprintln!("Skipping first-principles test (set CYRUS_FIRST_PRINCIPLES=1)");
        return false;
    }
    true
}

fn load_stage2_fixture() -> Stage2Fixture {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let data_dir = crate::mcallister_data_dir();
    assert!(
        !(crate::first_principles_enabled() && data_dir.is_none()),
        "CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests"
    );

    let points_raw = if let Some(dir) = data_dir {
        read_csv_rows_i64(&dir.join("points.dat"))
    } else {
        assert!(
            crate::fixtures_enabled(),
            "Set CYRUS_ALLOW_FIXTURES=1 to use JSON fixtures"
        );
        let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
        let content = std::fs::read_to_string(&input_path)
            .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
        let input: PolytopeInput = serde_json::from_str(&content)
            .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()));
        input.points
    };

    // Create polytope from all primal points
    let all_points: Vec<Point> = points_raw
        .iter()
        .map(|coords| Point::new(coords.clone()))
        .collect();

    let polytope = Polytope::from_vertices(all_points).expect("Failed to create polytope");

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

fn load_primal_points_raw() -> Vec<Vec<i64>> {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let data_dir = crate::mcallister_data_dir();
    assert!(
        !(crate::first_principles_enabled() && data_dir.is_none()),
        "CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests"
    );

    if let Some(dir) = data_dir {
        return read_csv_rows_i64(&dir.join("points.dat"));
    }

    assert!(
        crate::fixtures_enabled(),
        "Set CYRUS_ALLOW_FIXTURES=1 to use JSON fixtures"
    );
    let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
    let content = std::fs::read_to_string(&input_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
    let input: PolytopeInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()));
    input.points
}

fn read_csv_rows_usize(path: &PathBuf) -> Vec<Vec<usize>> {
    let content = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", path.display()));
    content
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|s| s.trim().parse::<usize>().expect("invalid integer"))
                .collect::<Vec<_>>()
        })
        .collect()
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
fn stage2_ours_frst_heights() {
    if !require_first_principles() {
        return;
    }
    let fixture = load_stage2_fixture();

    let (heights, triangulation) =
        compute_frst_heights(&fixture.triangulation_points, fixture.origin_idx)
            .expect("Failed to compute FRST heights");

    // Verify star property
    assert!(
        triangulation.is_star(fixture.origin_idx),
        "Triangulation must have star property"
    );

    insta::assert_json_snapshot!("heights_ours", heights);
}

#[test]
fn stage2_ours_triangulation_simplex_count() {
    if !require_first_principles() {
        return;
    }
    let fixture = load_stage2_fixture();

    let (_heights, triangulation) =
        compute_frst_heights(&fixture.triangulation_points, fixture.origin_idx)
            .expect("Failed to compute FRST heights");

    assert_eq!(
        triangulation.simplices().len(),
        1107,
        "Expected 1107 simplices in FRST triangulation"
    );
}

#[test]
fn stage2_ours_triangulation_simplices() {
    if !require_first_principles() {
        return;
    }
    let fixture = load_stage2_fixture();

    let (_heights, triangulation) =
        compute_frst_heights(&fixture.triangulation_points, fixture.origin_idx)
            .expect("Failed to compute FRST heights");

    // Sort simplices for deterministic snapshot
    let mut simplices: Vec<_> = triangulation.simplices().to_vec();
    simplices.sort();

    insta::assert_json_snapshot!("triangulation_simplices_ours", simplices);
}

#[test]
fn stage2_ours_all_simplices_contain_origin() {
    if !require_first_principles() {
        return;
    }
    let fixture = load_stage2_fixture();

    let (_heights, triangulation) =
        compute_frst_heights(&fixture.triangulation_points, fixture.origin_idx)
            .expect("Failed to compute FRST heights");

    // Verify star property explicitly
    for (i, simplex) in triangulation.simplices().iter().enumerate() {
        assert!(
            simplex.contains(&fixture.origin_idx),
            "Simplex {} does not contain origin (idx {}): {:?}",
            i,
            fixture.origin_idx,
            simplex
        );
    }
}

// =============================================================================
// "theirs" branch: McAllister's heights from data files
// =============================================================================

#[derive(Debug, Deserialize)]
struct HeightsInput {
    values: Vec<f64>,
}

fn load_mcallister_heights() -> Vec<f64> {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let data_dir = crate::mcallister_data_dir();
    assert!(
        !(crate::first_principles_enabled() && data_dir.is_none()),
        "CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests"
    );
    if let Some(dir) = data_dir {
        return read_csv_f64(&dir.join("heights.dat"));
    }
    assert!(
        crate::fixtures_enabled(),
        "Set CYRUS_ALLOW_FIXTURES=1 to use JSON fixtures"
    );
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
    assert_eq!(
        heights.len(),
        219,
        "McAllister heights should have 219 values"
    );
}

#[test]
fn stage2_theirs_triangulation() {
    if !require_first_principles() {
        return;
    }
    let fixture = load_stage2_fixture();
    let heights = load_mcallister_heights();

    let triangulation = compute_regular_triangulation(&fixture.triangulation_points, &heights)
        .expect("Failed to compute triangulation from McAllister heights");

    // Sort simplices for deterministic snapshot
    let mut simplices: Vec<_> = triangulation.simplices().to_vec();
    simplices.sort();

    insta::assert_json_snapshot!("triangulation_simplices_theirs", simplices);
}

#[test]
fn stage2_theirs_is_star() {
    if !require_first_principles() {
        return;
    }
    let fixture = load_stage2_fixture();
    let heights = load_mcallister_heights();

    let triangulation = compute_regular_triangulation(&fixture.triangulation_points, &heights)
        .expect("Failed to compute triangulation from McAllister heights");

    assert!(
        triangulation.is_star(fixture.origin_idx),
        "McAllister's triangulation must have star property"
    );
}

#[test]
fn stage2_dual_frst_matches_mcallister_dual_simplices_checkpoint() {
    if !require_first_principles() {
        return;
    }
    let Some(data_dir) = crate::mcallister_data_dir() else {
        panic!("CYRUS_MCALLISTER_DATA_DIR must be set for first-principles tests");
    };

    let primal_points = load_primal_points_raw();
    let all_points: Vec<Point> = primal_points
        .iter()
        .map(|coords| Point::new(coords.clone()))
        .collect();
    let primal = Polytope::from_vertices(all_points).expect("failed to create primal polytope");
    let dual = primal
        .compute_dual()
        .expect("failed to compute dual polytope");
    let dual_triangulation_points = dual
        .points_not_interior_to_facets()
        .expect("failed to filter dual triangulation points");
    let origin_idx = dual_triangulation_points
        .iter()
        .position(|point| point.coords().iter().all(|&x| x == 0))
        .expect("dual origin not found");

    let (_heights, triangulation) = compute_frst_heights(&dual_triangulation_points, origin_idx)
        .expect("failed to compute dual FRST triangulation");
    let mut actual = triangulation.simplices().to_vec();
    actual.sort();

    let mut expected = read_csv_rows_usize(&data_dir.join("dual_simplices.dat"));
    expected.sort();

    assert_eq!(
        dual_triangulation_points.len(),
        9,
        "dual triangulation should use origin plus non-facet-interior points"
    );
    assert_eq!(
        actual, expected,
        "Cyrus-computed dual FRST should match the McAllister dual_simplices.dat checkpoint"
    );
}
