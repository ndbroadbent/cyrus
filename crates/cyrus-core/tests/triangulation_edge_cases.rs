#![allow(missing_docs)]
use cyrus_core::{Point, compute_regular_triangulation};

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
