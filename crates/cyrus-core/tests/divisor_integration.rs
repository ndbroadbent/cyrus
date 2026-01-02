#![allow(missing_docs, clippy::too_many_lines)]
use cyrus_core::divisor::compute_divisor_jacobian;
use cyrus_core::{
    Point, compute_divisor_volumes, compute_intersection_numbers, compute_regular_triangulation,
};
use malachite::Integer;

#[test]
fn test_p4_divisor_volumes() {
    // P4 Quintic setup
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

    // Use star triangulation heights
    let heights = vec![-10.0, 1.0, 1.1, 1.2, 1.3, 1.4];
    let tri = compute_regular_triangulation(&points, &heights).unwrap();
    let kappa = compute_intersection_numbers(&tri, &points, &glsm).unwrap();

    // Kahler moduli t_i = 1.0 for all i (in basis)
    // Divisor volumes tau_i = 0.5 * int_X J^2 D_i
    // For P4 Quintic with t_i=1, J=5H. tau_i = 62.5.

    let t = vec![0.0, 1.0, 1.0, 1.0, 1.0, 1.0];

    let vols = compute_divisor_volumes(&kappa, &t);

    // Check non-origin volumes
    for (i, vol) in vols.iter().enumerate().skip(1) {
        println!("Vol D{i}: {vol}");
        assert!((vol - 62.5).abs() < 1e-10);
    }
}

#[test]
fn test_divisor_jacobian() {
    // P4 setup again
    let points = vec![
        Point::new(vec![0, 0, 0, 0]),
        Point::new(vec![-1, -1, -1, -1]),
        Point::new(vec![1, 0, 0, 0]),
        Point::new(vec![0, 1, 0, 0]),
        Point::new(vec![0, 0, 1, 0]),
        Point::new(vec![0, 0, 0, 1]),
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
    let heights = vec![-10.0, 1.0, 1.1, 1.2, 1.3, 1.4];
    let tri = compute_regular_triangulation(&points, &heights).unwrap();
    let kappa = compute_intersection_numbers(&tri, &points, &glsm).unwrap();
    let t = vec![0.0, 1.0, 1.0, 1.0, 1.0, 1.0];

    // Jacobian entries should be 25.

    let jac = compute_divisor_jacobian(&kappa, &t);

    for row in jac.iter().skip(1) {
        for val in row.iter().skip(1) {
            assert!((val - 25.0).abs() < 1e-10);
        }
    }
}
