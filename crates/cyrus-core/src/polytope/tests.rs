//! Tests for polytope module.

use super::*;
use crate::lattice::Point;
use std::collections::HashSet;

fn square() -> Polytope {
    Polytope::from_vertices(vec![
        Point::new(vec![-1, -1]),
        Point::new(vec![1, -1]),
        Point::new(vec![1, 1]),
        Point::new(vec![-1, 1]),
    ])
    .unwrap()
}

#[test]
fn test_polytope_creation() {
    let p = square();
    assert_eq!(p.dim(), 2);
    assert_eq!(p.num_vertices(), 4);
}

#[test]
fn test_origin_interior() {
    let p = square();
    assert!(p.origin_is_interior());
}

#[test]
fn test_origin_not_interior() {
    // Triangle not containing origin
    let p = Polytope::from_vertices(vec![
        Point::new(vec![1, 0]),
        Point::new(vec![2, 0]),
        Point::new(vec![1, 1]),
    ])
    .unwrap();
    assert!(!p.origin_is_interior());
}

#[test]
fn test_empty_vertices_error() {
    let result = Polytope::from_vertices(vec![]);
    assert!(result.is_err());
}

#[test]
fn test_dimension_mismatch_error() {
    let result = Polytope::from_vertices(vec![Point::new(vec![1, 2]), Point::new(vec![1, 2, 3])]);
    assert!(result.is_err());
}

#[test]
fn test_compute_dual_square() {
    // The dual of the unit square [-1,1]^2 with vertices (+/-1, +/-1)
    // is the diamond with lattice points: origin + vertices (+/-1, 0), (0, +/-1)
    let p = square();
    let dual = p.compute_dual().unwrap();

    assert_eq!(dual.dim(), 2);

    let dual_coords: HashSet<Vec<i64>> = dual
        .vertices()
        .iter()
        .map(|v| v.coords().to_vec())
        .collect();

    // Expected: all lattice points in dual (origin + 4 vertices)
    let expected: HashSet<Vec<i64>> = [
        vec![0, 0], // origin (interior)
        vec![1, 0],
        vec![-1, 0],
        vec![0, 1],
        vec![0, -1],
    ]
    .into_iter()
    .collect();

    assert_eq!(dual_coords, expected);
}

#[test]
fn test_compute_dual_not_reflexive() {
    // Triangle not containing origin - should fail
    let p = Polytope::from_vertices(vec![
        Point::new(vec![1, 0]),
        Point::new(vec![2, 0]),
        Point::new(vec![1, 1]),
    ])
    .unwrap();
    assert!(p.compute_dual().is_err());
}

#[test]
fn test_compute_dual_3d_cube() {
    // Unit cube [-1,1]^3 with 8 vertices
    let p = Polytope::from_vertices(vec![
        Point::new(vec![-1, -1, -1]),
        Point::new(vec![1, -1, -1]),
        Point::new(vec![-1, 1, -1]),
        Point::new(vec![1, 1, -1]),
        Point::new(vec![-1, -1, 1]),
        Point::new(vec![1, -1, 1]),
        Point::new(vec![-1, 1, 1]),
        Point::new(vec![1, 1, 1]),
    ])
    .unwrap();

    let dual = p.compute_dual().unwrap();

    // Dual of cube is octahedron with lattice points: origin + 6 vertices
    assert_eq!(dual.dim(), 3);
    assert_eq!(dual.num_vertices(), 7); // 6 vertices + 1 origin

    let dual_coords: HashSet<Vec<i64>> = dual
        .vertices()
        .iter()
        .map(|v| v.coords().to_vec())
        .collect();

    let expected: HashSet<Vec<i64>> = [
        vec![0, 0, 0], // origin (interior)
        vec![1, 0, 0],
        vec![-1, 0, 0],
        vec![0, 1, 0],
        vec![0, -1, 0],
        vec![0, 0, 1],
        vec![0, 0, -1],
    ]
    .into_iter()
    .collect();

    assert_eq!(dual_coords, expected);
}

#[test]
fn test_points_not_interior_to_facets_square() {
    // Square with 4 vertices - all are on 2 facets (edges)
    // No points interior to facets
    let p = square();
    let filtered = p.points_not_interior_to_facets().unwrap();

    // All 4 vertices should be kept (each on 2 edges)
    assert_eq!(filtered.len(), 4);
}

#[test]
fn test_points_not_interior_to_facets_with_edge_point() {
    // Square with an extra point on an edge
    // The edge point saturates exactly 1 facet, so it should be excluded
    let p = Polytope::from_vertices(vec![
        Point::new(vec![-1, -1]),
        Point::new(vec![1, -1]),
        Point::new(vec![1, 1]),
        Point::new(vec![-1, 1]),
        Point::new(vec![0, -1]), // Point on bottom edge - saturates 1 facet
    ])
    .unwrap();

    let filtered = p.points_not_interior_to_facets().unwrap();

    // Should have 4 vertices, edge point excluded
    assert_eq!(filtered.len(), 4);

    // Verify edge point is not in result
    let coords: Vec<Vec<i64>> = filtered.iter().map(|p| p.coords().to_vec()).collect();
    assert!(!coords.contains(&vec![0, -1]));
}

#[test]
fn test_points_not_interior_to_facets_with_interior() {
    // Square with origin inside
    let p = Polytope::from_vertices(vec![
        Point::new(vec![-1, -1]),
        Point::new(vec![1, -1]),
        Point::new(vec![1, 1]),
        Point::new(vec![-1, 1]),
        Point::new(vec![0, 0]), // Origin - interior, saturates 0 facets
    ])
    .unwrap();

    let filtered = p.points_not_interior_to_facets().unwrap();

    // Should have 5 points: 4 vertices + 1 interior
    assert_eq!(filtered.len(), 5);

    // Interior point should be FIRST (CYTools ordering)
    assert_eq!(filtered[0].coords(), &[0, 0]);
}

#[test]
fn test_points_not_interior_to_facets_ordering() {
    // Test that ordering is correct:
    // 1. Interior points first
    // 2. Then by decreasing saturation count
    let p = Polytope::from_vertices(vec![
        Point::new(vec![-1, -1]),
        Point::new(vec![1, -1]),
        Point::new(vec![1, 1]),
        Point::new(vec![-1, 1]),
        Point::new(vec![0, 0]), // Interior (0 saturations)
    ])
    .unwrap();

    let filtered = p.points_not_interior_to_facets().unwrap();

    // First point should be interior (0 saturations)
    assert_eq!(filtered[0].coords(), &[0, 0]);

    // Remaining points should all have saturation count >= 2
    // (they're on 2 edges each for a square)
}

#[test]
fn test_points_not_interior_to_facets_mcallister() {
    // Load McAllister polytope 4-214-647
    // Expected: 294 total points -> 219 after filtering (75 removed)
    let manifest_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let input_path = manifest_dir.join("tests/fixtures/4_214_647/polytope.json");

    #[derive(serde::Deserialize)]
    struct PolytopeInput {
        primal_points: Vec<Vec<i64>>,
    }

    let content = std::fs::read_to_string(&input_path)
        .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
    let input: PolytopeInput = serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()));

    let points: Vec<Point> = input
        .primal_points
        .iter()
        .map(|coords| Point::new(coords.clone()))
        .collect();

    assert_eq!(
        points.len(),
        294,
        "McAllister polytope should have 294 points"
    );

    let polytope = Polytope::from_vertices(points).unwrap();
    let filtered = polytope.points_not_interior_to_facets().unwrap();

    // McAllister uses 219 points for triangulation (75 are interior to facets)
    assert_eq!(
        filtered.len(),
        219,
        "Expected 219 points not interior to facets, got {}",
        filtered.len()
    );

    // Origin should be first (interior point with 0 saturations)
    assert!(
        filtered[0].coords().iter().all(|&x| x == 0),
        "First point should be the origin"
    );
}

#[test]
fn test_faces_4d_for_projective_simplex_points() {
    let points = vec![
        Point::new(vec![1, 0, 0, 0]),
        Point::new(vec![0, 1, 0, 0]),
        Point::new(vec![0, 0, 1, 0]),
        Point::new(vec![0, 0, 0, 1]),
        Point::new(vec![-1, -1, -1, -1]),
    ];
    let polytope = Polytope::from_vertices(points.clone()).unwrap();

    let faces = polytope.faces_4d_for_points(&points).unwrap();

    assert_eq!(faces.facets.len(), 5);
    assert!(faces.facets.iter().all(|facet| facet.len() == 4));
    assert_eq!(faces.twofaces.len(), 10);
    assert!(faces.twofaces.iter().all(|twoface| twoface.len() == 3));
}
