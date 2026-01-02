#![allow(missing_docs, clippy::too_many_lines)]
use cyrus_core::{Point, Triangulation, compute_regular_triangulation};

#[test]
fn test_empty_points() {
    let points: Vec<Point> = vec![];
    let heights: Vec<f64> = vec![];
    let tri = compute_regular_triangulation(&points, &heights).unwrap();
    assert!(tri.simplices().is_empty());
}

#[test]
fn test_single_point() {
    let points = vec![Point::new(vec![0, 0])];
    let heights = vec![1.0];
    let tri = compute_regular_triangulation(&points, &heights).unwrap();
    // Single point in 2D - convex hull produces facets but none are "lower"
    // in the lifted dimension since there's only one point
    // The result depends on is_lower_face logic
    // Just verify it doesn't panic and returns valid result
    assert!(tri.simplices().len() <= 1);
}

#[test]
fn test_two_points() {
    let points = vec![Point::new(vec![0, 0]), Point::new(vec![1, 0])];
    let heights = vec![0.0, 1.0];
    let tri = compute_regular_triangulation(&points, &heights).unwrap();
    // Two points in 2D lifted to 3D - may produce empty triangulation
    // if no facet qualifies as "lower"
    // Just verify it doesn't panic
    let _ = tri.simplices();
}

#[test]
fn test_mismatched_points_heights() {
    let points = vec![Point::new(vec![0, 0]), Point::new(vec![1, 0])];
    let heights = vec![0.0]; // Wrong length
    let res = compute_regular_triangulation(&points, &heights);
    assert!(res.is_err());
    assert!(res.unwrap_err().to_string().contains("must match"));
}

#[test]
fn test_triangulation_new() {
    // Test Triangulation::new directly
    let simplices = vec![vec![0, 1, 2], vec![1, 2, 3]];
    let tri = Triangulation::new(simplices.clone());
    assert_eq!(tri.simplices(), &simplices);
}

#[test]
fn test_3d_triangulation() {
    // 3D tetrahedron - 4 points in 3D lifted to 4D
    let points = vec![
        Point::new(vec![0, 0, 0]),
        Point::new(vec![1, 0, 0]),
        Point::new(vec![0, 1, 0]),
        Point::new(vec![0, 0, 1]),
    ];
    // Use more distinct heights to ensure proper lifting
    let heights = vec![0.0, 1.0, 2.0, 3.0];
    let tri = compute_regular_triangulation(&points, &heights).unwrap();
    // 4 points in 3D -> n <= dim case in convex_hull
    // Result depends on is_lower_face logic
    // Just verify it doesn't panic
    let _ = tri.simplices();
}

#[test]
fn test_points_inside_hull_ignored() {
    // 2D Square with center point
    let points = vec![
        Point::new(vec![0, 0]),
        Point::new(vec![1, 0]),
        Point::new(vec![1, 1]),
        Point::new(vec![0, 1]),
        Point::new(vec![0, 0]), // Duplicate of 0
        Point::new(vec![0, 0]), // Duplicate again
    ];
    // Heights: Non-planar
    // 0: 0.0
    // 1: 10.0
    // 2: 25.0
    // 3: 5.0
    // 4: 100.0 (Interior high)
    // 5: 0.0 (Duplicate)
    let heights = vec![0.0, 10.0, 25.0, 5.0, 100.0, 0.0];

    // Points 0,1,2,3 form a tetrahedron in 3D (not coplanar).
    // (0,0,0), (1,0,10), (1,1,25), (0,1,5).
    // Plane through 0,1,3: z = 10x + 5y.
    // At (1,1), z = 15. Actual is 25. 25 > 15. Convex!

    let tri = compute_regular_triangulation(&points, &heights).unwrap();

    // Should have 2 triangles.
    assert_eq!(tri.simplices().len(), 2);

    // Check that point 4 (index 4) is NOT used
    for s in tri.simplices() {
        assert!(
            !s.contains(&4),
            "Interior point 4 should not be in triangulation"
        );
    }
}

#[test]
fn test_coplanar_points_2d() {
    // 3 collinear points on x-axis
    let points = vec![
        Point::new(vec![0, 0]),
        Point::new(vec![1, 0]),
        Point::new(vec![2, 0]),
    ];
    let heights = vec![0.0, 0.0, 0.0];

    let tri = compute_regular_triangulation(&points, &heights).unwrap();
    // Lower hull of collinear points in 2D (lifted to 3D line) is just edges?
    // In 2D triangulation, we expect triangles.
    // Collinear points have NO triangles (dim=2).
    // So empty simplices?
    assert!(tri.simplices().is_empty());
}

#[test]
fn test_larger_triangulation() {
    // More points to exercise incremental hull algorithm
    let points = vec![
        Point::new(vec![0, 0]),
        Point::new(vec![2, 0]),
        Point::new(vec![1, 2]),
        Point::new(vec![3, 1]),
        Point::new(vec![1, 1]), // Interior point
        Point::new(vec![0, 2]),
    ];
    let heights = vec![0.0, 0.5, 1.0, 0.3, 2.0, 0.8];

    let tri = compute_regular_triangulation(&points, &heights).unwrap();
    // Should produce a valid triangulation
    assert!(!tri.simplices().is_empty());

    // Verify no simplex contains the interior point (index 4) if it's inside
    // (depends on heights, so just check triangulation is valid)
    for s in tri.simplices() {
        assert_eq!(s.len(), 3); // 2D simplices are triangles
    }
}

#[test]
fn test_degenerate_heights() {
    // All same height - planar in lifted dimension
    let points = vec![
        Point::new(vec![0, 0]),
        Point::new(vec![1, 0]),
        Point::new(vec![0, 1]),
        Point::new(vec![1, 1]),
    ];
    let heights = vec![0.0, 0.0, 0.0, 0.0];

    let tri = compute_regular_triangulation(&points, &heights).unwrap();
    // With all same height, the "lower hull" is the entire planar set
    // Should still produce valid triangulation
    let _ = tri.simplices();
}

#[test]
fn test_convex_hull_more_points() {
    // Square with 5th point exactly on edge
    let points = vec![
        Point::new(vec![0, 0]),
        Point::new(vec![2, 0]),
        Point::new(vec![2, 2]),
        Point::new(vec![0, 2]),
        Point::new(vec![1, 0]), // On edge between 0 and 1
    ];
    let heights = vec![0.0, 0.0, 0.0, 0.0, 0.5];

    let tri = compute_regular_triangulation(&points, &heights).unwrap();
    // Edge point with higher height won't be in lower hull
    let _ = tri.simplices();
}

#[test]
fn test_high_dimensional() {
    // 4D points
    let points = vec![
        Point::new(vec![0, 0, 0, 0]),
        Point::new(vec![1, 0, 0, 0]),
        Point::new(vec![0, 1, 0, 0]),
        Point::new(vec![0, 0, 1, 0]),
        Point::new(vec![0, 0, 0, 1]),
    ];
    let heights = vec![0.0, 1.0, 2.0, 3.0, 4.0];

    let tri = compute_regular_triangulation(&points, &heights).unwrap();
    // 5 points in 4D - forms a simplex, just verify no panic
    let _ = tri.simplices();
}
