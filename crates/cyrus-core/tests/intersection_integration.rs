#![allow(missing_docs, clippy::too_many_lines)]
use cyrus_core::{
    Point, Triangulation, compute_intersection_numbers, compute_regular_triangulation,
};
use malachite::{Integer, Rational};

#[test]
fn test_p4_intersection_numbers() {
    // Construct P4 vertices: standard simplex in 4D
    // v0=(0,0,0,0), v1=(-1,-1,-1,-1), v2=(1,0,0,0), v3=(0,1,0,0), v4=(0,0,1,0), v5=(0,0,0,1)
    let points = vec![
        Point::new(vec![0, 0, 0, 0]),     // 0
        Point::new(vec![-1, -1, -1, -1]), // 1
        Point::new(vec![1, 0, 0, 0]),     // 2
        Point::new(vec![0, 1, 0, 0]),     // 3
        Point::new(vec![0, 0, 1, 0]),     // 4
        Point::new(vec![0, 0, 0, 1]),     // 5
    ];

    let glsm = vec![
        vec![
            Integer::from(0),
            Integer::from(-1),
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
        ],
        vec![
            Integer::from(0),
            Integer::from(-1),
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
        ],
        vec![
            Integer::from(0),
            Integer::from(-1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
        ],
        vec![
            Integer::from(0),
            Integer::from(-1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ],
    ];

    // Star triangulation: 5 simplices, each containing 0 and 4 vertices
    let tri = Triangulation::new(vec![
        vec![0, 1, 2, 3, 4],
        vec![0, 1, 2, 3, 5],
        vec![0, 1, 2, 4, 5],
        vec![0, 1, 3, 4, 5],
        vec![0, 2, 3, 4, 5],
    ]);

    // Compute intersection numbers
    let kappa = compute_intersection_numbers(&tri, &points, &glsm)
        .expect("Intersection computation failed");

    // Check kappa_234 (rays v2, v3, v4)
    // Sum over l: kappa_234l.
    // Term l=1 (v1): {1,2,3,4} is a facet.
    // Term l=5 (v5): {2,3,4,5} is a facet.
    // Others might be non-zero due to relations.
    let val = kappa.get(2, 3, 4);
    println!("kappa_234 = {val}");
    assert_eq!(val, Rational::from(5));
}

#[test]
fn test_intersection_with_triangulation() {
    // Same P4 but using compute_regular_triangulation to ensure integration
    let points = vec![
        Point::new(vec![0, 0, 0, 0]),
        Point::new(vec![-1, -1, -1, -1]),
        Point::new(vec![1, 0, 0, 0]),
        Point::new(vec![0, 1, 0, 0]),
        Point::new(vec![0, 0, 1, 0]),
        Point::new(vec![0, 0, 0, 1]),
    ];
    // Use linear relations again
    let glsm = vec![
        vec![
            Integer::from(0),
            Integer::from(-1),
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
        ],
        vec![
            Integer::from(0),
            Integer::from(-1),
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
        ],
        vec![
            Integer::from(0),
            Integer::from(-1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
        ],
        vec![
            Integer::from(0),
            Integer::from(-1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
        ],
    ];

    // For a star triangulation, the origin must be part of all simplices.
    // We lift the origin very low, others high.
    let heights = vec![-10.0, 1.0, 1.1, 1.2, 1.3, 1.4]; // Slight perturbation to be regular

    let tri = compute_regular_triangulation(&points, &heights).expect("Triangulation failed");
    // Should have 5 simplices, each containing 0.
    assert!(tri.simplices().len() >= 5);
    for s in tri.simplices() {
        assert!(s.contains(&0));
    }

    let kappa = compute_intersection_numbers(&tri, &points, &glsm).expect("Intersection failed");

    // Check D2.D3.D4 (indices 2, 3, 4)
    let val = kappa.get(2, 3, 4);
    println!("kappa_234 = {val}");
    assert_eq!(val, Rational::from(5));
}
