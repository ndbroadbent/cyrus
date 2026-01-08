//! McAllister Pipeline Stage 3: GLSM Charges + Intersection Numbers
//!
//! Computes:
//! - GLSM charge matrix from triangulation points
//! - Intersection numbers κ_ijk from triangulation + GLSM
//!
//! - "ours" branch: Using our FRST triangulation
//! - "theirs" branch: Using McAllister's triangulation (from their heights)

#![allow(missing_docs)]

use serde::Deserialize;
use std::path::PathBuf;

use cyrus_core::{
    Point, Polytope, compute_glsm_charge_matrix,
};

#[derive(Debug, Deserialize)]
struct PolytopeInput {
    points: Vec<Vec<i64>>,
}

#[derive(Debug, Deserialize)]
struct HeightsInput {
    values: Vec<f64>,
}

struct Stage3Fixture {
    triangulation_points: Vec<Point>,
    origin_idx: usize,
}

fn load_fixture() -> Stage3Fixture {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let input_path = manifest_dir.join("tests/mcallister_e2e/inputs/polytope.json");
    let content = std::fs::read_to_string(&input_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
    let input: PolytopeInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()));

    let all_points: Vec<Point> = input
        .points
        .iter()
        .map(|coords| Point::new(coords.clone()))
        .collect();

    let polytope = Polytope::from_vertices(all_points).expect("Failed to create polytope");

    let triangulation_points = polytope
        .points_not_interior_to_facets()
        .expect("Failed to filter points");

    let origin_idx = triangulation_points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .expect("Origin not found");

    Stage3Fixture {
        triangulation_points,
        origin_idx,
    }
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

// =============================================================================
// GLSM charge matrix (same for both branches - depends only on points)
// =============================================================================

#[test]
fn stage3_glsm_charge_matrix() {
    let fixture = load_fixture();

    let glsm = compute_glsm_charge_matrix(&fixture.triangulation_points, true)
        .expect("Failed to compute GLSM charge matrix");

    // GLSM is the kernel of the point matrix
    // With include_origin=true, we have 219 + 1 = 220 points
    // Kernel has (n_points - dim - 1) = 220 - 4 - 1 = 215 rows
    // Each row has 220 columns
    assert!(!glsm.is_empty(), "GLSM should not be empty");
    assert_eq!(
        glsm.len(),
        215,
        "GLSM should have 215 kernel vectors (h11+1)"
    );
    assert_eq!(
        glsm[0].len(),
        220,
        "GLSM columns should equal n_points+origin = 220"
    );

    // Snapshot the GLSM matrix - convert Integer to i64
    let glsm_i64: Vec<Vec<i64>> = glsm
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).expect("GLSM value fits in i64"))
                .collect()
        })
        .collect();

    insta::assert_json_snapshot!("glsm_charge_matrix", glsm_i64);
}

// =============================================================================
// "ours" branch: Intersection numbers from our FRST triangulation
// =============================================================================

#[test]
#[cfg(feature = "slow-tests")]
fn stage3_ours_intersection_numbers() {
    let fixture = load_fixture();

    // Compute our triangulation
    let (_heights, triangulation) =
        compute_frst_heights(&fixture.triangulation_points, fixture.origin_idx)
            .expect("Failed to compute FRST heights");

    // Compute GLSM
    let glsm = compute_glsm_charge_matrix(&fixture.triangulation_points, true)
        .expect("Failed to compute GLSM");

    // Compute intersection numbers
    let kappa = compute_intersection_numbers(&triangulation, &fixture.triangulation_points, &glsm)
        .expect("Failed to compute intersection numbers");

    // Convert to serializable format: sorted list of ((i,j,k), value)
    let mut entries: Vec<((usize, usize, usize), String)> = kappa
        .iter()
        .map(|(&key, val)| (key, val.get().to_string()))
        .collect();
    entries.sort_by_key(|(key, _)| *key);

    #[derive(Serialize)]
    struct IntersectionSnapshot {
        dim: usize,
        num_nonzero: usize,
        entries: Vec<((usize, usize, usize), String)>,
    }

    let snapshot = IntersectionSnapshot {
        dim: kappa.dim(),
        num_nonzero: kappa.num_nonzero(),
        entries,
    };

    insta::assert_json_snapshot!("intersection_ours", snapshot);
}

// =============================================================================
// "theirs" branch: Intersection numbers from McAllister's triangulation
// =============================================================================

#[test]
#[cfg(feature = "slow-tests")]
fn stage3_theirs_intersection_numbers() {
    let fixture = load_fixture();
    let heights = load_mcallister_heights();

    // Compute triangulation from their heights
    let triangulation = compute_regular_triangulation(&fixture.triangulation_points, &heights)
        .expect("Failed to compute triangulation");

    // Compute GLSM
    let glsm = compute_glsm_charge_matrix(&fixture.triangulation_points, true)
        .expect("Failed to compute GLSM");

    // Compute intersection numbers
    let kappa = compute_intersection_numbers(&triangulation, &fixture.triangulation_points, &glsm)
        .expect("Failed to compute intersection numbers");

    // Convert to serializable format: sorted list of ((i,j,k), value)
    let mut entries: Vec<((usize, usize, usize), String)> = kappa
        .iter()
        .map(|(&key, val)| (key, val.get().to_string()))
        .collect();
    entries.sort_by_key(|(key, _)| *key);

    #[derive(Serialize)]
    struct IntersectionSnapshot {
        dim: usize,
        num_nonzero: usize,
        entries: Vec<((usize, usize, usize), String)>,
    }

    let snapshot = IntersectionSnapshot {
        dim: kappa.dim(),
        num_nonzero: kappa.num_nonzero(),
        entries,
    };

    insta::assert_json_snapshot!("intersection_theirs", snapshot);
}
