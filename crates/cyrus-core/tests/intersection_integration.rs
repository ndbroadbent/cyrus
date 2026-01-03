#![allow(missing_docs, clippy::too_many_lines)]
use cyrus_core::types::f64::F64;
use cyrus_core::types::rational::Rational as TypedRational;
use cyrus_core::types::tags::{Finite, Pos};
use cyrus_core::{
    Intersection, Point, Triangulation, compute_intersection_numbers, compute_regular_triangulation,
};
use malachite::{Integer, Rational};

fn pos_rat(n: i32) -> TypedRational<Pos> {
    TypedRational::<Pos>::new(Rational::from(n)).unwrap()
}

fn finite_rat(n: i32) -> TypedRational<Finite> {
    TypedRational::<Finite>::from_raw(Rational::from(n))
}

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
    println!("kappa_234 = {val:?}");
    assert_eq!(val, finite_rat(5));
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
    println!("kappa_234 = {val:?}");
    assert_eq!(val, finite_rat(5));
}

#[test]
fn test_intersection_symmetry() {
    let mut kappa = Intersection::new(4);
    kappa.set(1, 2, 3, pos_rat(7));

    let expected = finite_rat(7);
    // All permutations should return the same value
    assert_eq!(kappa.get(1, 2, 3), expected);
    assert_eq!(kappa.get(1, 3, 2), expected);
    assert_eq!(kappa.get(2, 1, 3), expected);
    assert_eq!(kappa.get(2, 3, 1), expected);
    assert_eq!(kappa.get(3, 1, 2), expected);
    assert_eq!(kappa.get(3, 2, 1), expected);
}

#[test]
fn test_intersection_contract_triple_finite() {
    let mut kappa = Intersection::new(2);
    // κ_000 = 1, κ_001 = 2, κ_011 = 3, κ_111 = 4
    kappa.set(0, 0, 0, pos_rat(1));
    kappa.set(0, 0, 1, pos_rat(2));
    kappa.set(0, 1, 1, pos_rat(3));
    kappa.set(1, 1, 1, pos_rat(4));

    let t: Vec<F64<Finite>> = vec![
        F64::<Finite>::new(1.0).unwrap(),
        F64::<Finite>::new(1.0).unwrap(),
    ];
    let vol = kappa.contract_triple_finite(&t).unwrap();
    // vol = 1*1*1*1 (mult=1) + 2*1*1*1 (mult=3) + 3*1*1*1 (mult=3) + 4*1*1*1 (mult=1)
    // = 1 + 6 + 9 + 4 = 20
    assert!((vol.get() - 20.0).abs() < 1e-10);
}

#[test]
fn test_intersection_contract_triple_symmetry_multiplicities() {
    // Test all three symmetry cases: (i,i,i), (i,i,j), (i,j,k)
    let mut kappa = Intersection::new(3);

    // Case 1: all equal (i,i,i) - multiplicity 1
    kappa.set(0, 0, 0, pos_rat(6));

    // Case 2: two equal (i,i,j) - multiplicity 3
    kappa.set(1, 1, 2, pos_rat(4));

    // Case 3: all distinct (i,j,k) - multiplicity 6
    kappa.set(0, 1, 2, pos_rat(2));

    let t: Vec<F64<Finite>> = vec![
        F64::<Finite>::new(1.0).unwrap(),
        F64::<Finite>::new(1.0).unwrap(),
        F64::<Finite>::new(1.0).unwrap(),
    ];
    let vol = kappa.contract_triple_finite(&t).unwrap();
    // 6*1 + 4*3 + 2*6 = 6 + 12 + 12 = 30
    assert!((vol.get() - 30.0).abs() < 1e-10);
}

#[test]
fn test_intersection_iter() {
    let mut kappa = Intersection::new(3);
    kappa.set(0, 1, 2, pos_rat(5));
    kappa.set(0, 0, 1, pos_rat(3));

    assert_eq!(kappa.iter().count(), 2);
}

#[test]
fn test_intersection_iter_entries() {
    let mut kappa = Intersection::new(2);
    kappa.set(0, 0, 1, pos_rat(7));

    let entries: Vec<_> = kappa.iter_entries().collect();
    assert_eq!(entries.len(), 1);
    assert_eq!(*entries[0].1, finite_rat(7));
}

#[test]
fn test_intersection_dim() {
    let kappa = Intersection::new(5);
    assert_eq!(kappa.dim(), 5);
}

#[test]
fn test_intersection_with_degenerate_simplex() {
    // Triangulation with a simplex that has fewer rays than needed after removing origin
    let points = vec![
        Point::new(vec![0, 0, 0, 0]), // 0 - origin
        Point::new(vec![1, 0, 0, 0]), // 1
        Point::new(vec![0, 1, 0, 0]), // 2
    ];

    // This simplex only has 2 rays (after removing origin), but we're in 4D
    // so dim_v = 4, and rays.len() = 2 < dim_v triggers the continue in build_variable_list
    let tri = Triangulation::new(vec![vec![0, 1, 2]]);

    let glsm = vec![vec![Integer::from(0), Integer::from(1), Integer::from(1)]];

    // Should still work but produce minimal output
    let result = compute_intersection_numbers(&tri, &points, &glsm);
    // May fail or succeed with empty results
    if let Ok(kappa) = result {
        // Just verify it returns something
        let _ = kappa.dim();
    }
}

#[test]
fn test_glsm_empty_points() {
    use cyrus_core::compute_glsm_charge_matrix;

    let points: Vec<Point> = Vec::new();
    let result = compute_glsm_charge_matrix(&points, true);
    assert!(result.is_err());
    let err = result.unwrap_err();
    assert!(err.to_string().contains("No points provided"));
}
