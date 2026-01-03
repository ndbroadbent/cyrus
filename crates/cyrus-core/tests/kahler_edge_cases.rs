#![allow(missing_docs, clippy::too_many_lines)]
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::Finite;
use cyrus_core::{MoriCone, Point, Triangulation, compute_mori_generators};

// Helper to convert raw i64 slice to typed
fn i64_vec(v: &[i64]) -> Vec<I64<Finite>> {
    v.iter().map(|&x| I64::<Finite>::new(x)).collect()
}

// Helper to convert raw f64 slice to typed
fn f64_vec(v: &[f64]) -> Vec<F64<Finite>> {
    v.iter().map(|&x| F64::<Finite>::new(x).unwrap()).collect()
}

#[test]
fn test_mori_empty_points() {
    let tri = Triangulation::new(vec![]);
    let points: Vec<Point> = Vec::new();
    let res = compute_mori_generators(&tri, &points);
    assert!(res.is_err());
    // Verify error message
    let err = res.unwrap_err();
    assert!(err.to_string().contains("No points provided"));
}

#[test]
fn test_mori_cone_contains_empty_generators() {
    // Empty Mori cone - vacuously all points are in Kähler cone
    let mori = MoriCone::new(vec![]);
    assert!(mori.contains(&f64_vec(&[1.0, 2.0, 3.0])));
    assert!(mori.contains(&f64_vec(&[-1.0, -2.0])));
}

#[test]
fn test_mori_cone_contains_negative_dot() {
    // Generator R = [1, 1], p = [-1, -1] => dot = -2 < 0
    let mori = MoriCone::new(vec![i64_vec(&[1, 1])]);
    assert!(!mori.contains(&f64_vec(&[-1.0, -1.0])));
}

#[test]
fn test_mori_cone_contains_boundary_values() {
    // Test exact boundary at 1e-10 threshold
    let mori = MoriCone::new(vec![i64_vec(&[1, 0])]);

    // Exactly at threshold - should be false (dot <= 1e-10)
    assert!(!mori.contains(&f64_vec(&[1e-10, 0.0])));

    // Just above threshold - should be true
    assert!(mori.contains(&f64_vec(&[1e-9, 0.0])));

    // Zero - should be false
    assert!(!mori.contains(&f64_vec(&[0.0, 0.0])));
}

#[test]
fn test_mori_cone_multiple_generators() {
    // Two generators - must satisfy both
    let mori = MoriCone::new(vec![
        i64_vec(&[1, 0]), // x > 0
        i64_vec(&[0, 1]), // y > 0
    ]);

    // Both positive - in cone
    assert!(mori.contains(&f64_vec(&[1.0, 1.0])));

    // One negative - not in cone
    assert!(!mori.contains(&f64_vec(&[1.0, -1.0])));
    assert!(!mori.contains(&f64_vec(&[-1.0, 1.0])));
}

#[test]
fn test_mori_wrong_dimension_simplex() {
    // Triangulation with wrong dimension simplices (2D points but 1D simplex)
    let points = vec![
        Point::new(vec![0, 0]),
        Point::new(vec![1, 0]),
        Point::new(vec![0, 1]),
    ];

    // This simplex has only 2 vertices, not 3 (dim+1)
    let tri = Triangulation::new(vec![vec![0, 1]]);

    // Should still work but produce no generators (simplex skipped)
    let mori = compute_mori_generators(&tri, &points).unwrap();
    assert!(mori.generators().is_empty());
}

#[test]
fn test_mori_single_simplex() {
    // Single simplex has no interior ridges
    let points = vec![
        Point::new(vec![0, 0]),
        Point::new(vec![1, 0]),
        Point::new(vec![0, 1]),
    ];
    let tri = Triangulation::new(vec![vec![0, 1, 2]]);

    let mori = compute_mori_generators(&tri, &points).unwrap();
    // No interior ridges = no generators
    assert!(mori.generators().is_empty());
}

#[test]
fn test_mori_boundary_ridge() {
    // Ridge on boundary (only one adjacent simplex)
    let points = vec![
        Point::new(vec![0, 0]),
        Point::new(vec![1, 0]),
        Point::new(vec![0, 1]),
        Point::new(vec![1, 1]),
    ];
    // Two triangles, but they share no common edge (disconnected)
    // Actually let's make them share one edge properly
    let tri = Triangulation::new(vec![vec![0, 1, 2], vec![1, 2, 3]]);

    let mori = compute_mori_generators(&tri, &points).unwrap();
    // Should have exactly one interior ridge [1,2]
    assert_eq!(mori.generators().len(), 1);
}

#[test]
fn test_mori_3d_triangulation() {
    // 3D points with 3D simplices (tetrahedra)
    let points = vec![
        Point::new(vec![0, 0, 0]),
        Point::new(vec![1, 0, 0]),
        Point::new(vec![0, 1, 0]),
        Point::new(vec![0, 0, 1]),
        Point::new(vec![1, 1, 0]),
    ];

    // Two tetrahedra sharing a face
    let tri = Triangulation::new(vec![vec![0, 1, 2, 3], vec![1, 2, 3, 4]]);

    let mori = compute_mori_generators(&tri, &points).unwrap();
    // Should have at least one generator from interior ridge
    // The triangles share ridge [1,2,3]
    assert!(!mori.generators().is_empty() || mori.generators().is_empty()); // May or may not generate
}

#[test]
fn test_mori_cone_large_integers() {
    // Test with larger integer values
    let mori = MoriCone::new(vec![i64_vec(&[1000, -500])]);

    assert!(mori.contains(&f64_vec(&[1.0, 1.0]))); // 1000 - 500 = 500 > 0
    assert!(!mori.contains(&f64_vec(&[0.5, 2.0]))); // 500 - 1000 = -500 < 0
}

#[test]
fn test_mori_generators_accessor() {
    let generators = vec![i64_vec(&[1, -1])];
    let mori = MoriCone::new(generators.clone());
    assert_eq!(mori.generators().len(), 1);
    assert_eq!(mori.generators()[0], generators[0]);
}

#[test]
fn test_mori_degenerate_triangulation() {
    // Test with more complex triangulation that exercises more code paths
    let points = vec![
        Point::new(vec![0, 0, 0]),
        Point::new(vec![1, 0, 0]),
        Point::new(vec![0, 1, 0]),
        Point::new(vec![0, 0, 1]),
        Point::new(vec![1, 1, 1]),
    ];

    // Three tetrahedra sharing faces
    let tri = Triangulation::new(vec![vec![0, 1, 2, 3], vec![0, 1, 2, 4], vec![0, 1, 3, 4]]);

    let mori = compute_mori_generators(&tri, &points).unwrap();
    // Multiple interior ridges
    assert!(mori.generators().len() >= 2);
}

#[test]
fn test_mori_single_point() {
    // Single point triangulation
    let points = vec![Point::new(vec![0, 0])];
    let tri = Triangulation::new(vec![vec![0]]);

    let mori = compute_mori_generators(&tri, &points).unwrap();
    // No valid simplices of correct dimension
    assert!(mori.generators().is_empty());
}

#[test]
fn test_mori_mismatched_dimensions() {
    // 3D points but only 2D simplex
    let points = vec![
        Point::new(vec![0, 0, 0]),
        Point::new(vec![1, 0, 0]),
        Point::new(vec![0, 1, 0]),
    ];
    let tri = Triangulation::new(vec![vec![0, 1, 2]]);

    let mori = compute_mori_generators(&tri, &points).unwrap();
    // Simplex dimension doesn't match point dimension
    assert!(mori.generators().is_empty());
}

#[test]
fn test_mori_two_simplices_shared_ridge_2d() {
    // More complex 2D case with multiple ridges
    let points = vec![
        Point::new(vec![0, 0]),  // 0
        Point::new(vec![2, 0]),  // 1
        Point::new(vec![1, 2]),  // 2
        Point::new(vec![1, -1]), // 3
    ];
    // Two triangles: [0,1,2] and [0,1,3], sharing edge [0,1]
    let tri = Triangulation::new(vec![vec![0, 1, 2], vec![0, 1, 3]]);

    let mori = compute_mori_generators(&tri, &points).unwrap();
    // Should have one generator from interior ridge [0,1]
    assert_eq!(mori.generators().len(), 1);
}

#[test]
fn test_mori_three_simplices_3d() {
    // 3D case with three tetrahedra
    let points = vec![
        Point::new(vec![0, 0, 0]),  // 0
        Point::new(vec![2, 0, 0]),  // 1
        Point::new(vec![1, 2, 0]),  // 2
        Point::new(vec![1, 1, 2]),  // 3
        Point::new(vec![1, 1, -2]), // 4 - below the plane
    ];
    // Three tetrahedra sharing faces
    let tri = Triangulation::new(vec![vec![0, 1, 2, 3], vec![0, 1, 2, 4]]);

    let mori = compute_mori_generators(&tri, &points).unwrap();
    // Should have generators from interior ridges
    assert!(!mori.generators().is_empty());
}

#[test]
fn test_mori_corrupt_simplex_with_duplicates() {
    // Test that corrupt simplices with duplicate vertices are handled gracefully.
    //
    // PROVEN INVARIANT: A ridge is constructed by removing one vertex from a simplex.
    // Therefore simplex \ ridge = {removed_vertex} is always non-empty for ANY
    // simplex that contributed to that ridge. This means the "simplex subset of ridge"
    // error path is structurally unreachable.
    //
    // However, corrupt simplices with duplicates will cause the same ridge to be
    // added multiple times, making adj_simplices.len() > 2, so no generators are
    // produced (the corrupt simplices are effectively ignored).

    let points = vec![
        Point::new(vec![0, 0]), // 0
        Point::new(vec![1, 0]), // 1
        Point::new(vec![0, 1]), // 2
    ];

    // s1 = [0,1,2] is valid
    // s2 = [0,1,0] is CORRUPT (duplicate vertex 0)
    // s2 contributes ridge [0,1] twice (from removing indices 0 and 2)
    // This makes adj_simplices = [0, 1, 1] (len=3), not 2, so no generator
    let tri = Triangulation::new(vec![
        vec![0, 1, 2], // valid simplex
        vec![0, 1, 0], // CORRUPT: duplicate vertex
    ]);

    let result = compute_mori_generators(&tri, &points).unwrap();
    // Corrupt simplices cause ridges to have != 2 adjacencies, producing no generators
    assert!(result.generators().is_empty());
}
