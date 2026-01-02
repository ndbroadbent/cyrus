#![allow(missing_docs)]
use cyrus_core::{
    Intersection, compute_flat_direction_full, compute_n_matrix, solve_linear_system,
};
use malachite::Rational;

#[test]
fn test_solve_singular_system() {
    // 2x2 singular matrix [1, 1; 1, 1]
    let n = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
    let k = vec![1, 2]; // Inconsistent

    let result = solve_linear_system(&n, &k);
    assert!(result.is_none(), "Should return None for singular matrix");
}

#[test]
fn test_compute_n_matrix_permutations() {
    let mut kappa = Intersection::new(4);
    kappa.set(1, 2, 3, Rational::from(1));

    let m = vec![0, 0, 0, 1]; // Indices 0..3. 3 is the one.

    let n = compute_n_matrix(&kappa, &m);

    assert!((n[(1, 2)] - 1.0).abs() < 1e-10);
    assert!((n[(2, 1)] - 1.0).abs() < 1e-10);
    assert!((n[(1, 3)] - 0.0).abs() < 1e-10);

    // Also check self-intersections
    let mut kappa2 = Intersection::new(1);
    kappa2.set(0, 0, 0, Rational::from(1));
    let m2 = vec![1];
    let n2 = compute_n_matrix(&kappa2, &m2);
    assert!((n2[(0, 0)] - 1.0).abs() < 1e-10);
}

#[test]
fn test_compute_flat_direction_full() {
    let mut kappa = Intersection::new(1);
    kappa.set(0, 0, 0, Rational::from(6));
    let k = vec![3];
    let m = vec![1];

    let res = compute_flat_direction_full(&kappa, &k, &m).unwrap();
    assert!((res.p[0] - 0.5).abs() < 1e-10);
    // kappa(p,p,p) = 6 * 0.5^3 = 0.75.
    // ek0 = 1 / (4/3 * 0.75) = 1 / 1 = 1.0.
    assert!((res.ek0 - 1.0).abs() < 1e-10);
    assert_eq!(res.n_matrix.len(), 1);
}

#[test]
fn test_compute_n_matrix_parallel() {
    // Trigger parallel path (>100 nonzero entries)
    let dim = 10;
    let mut kappa = Intersection::new(dim);
    // Fill enough distinct entries
    // Total distinct entries for dim=10 is binom(10+3-1, 3) = binom(12,3) = 220.
    for i in 0..dim {
        for j in i..dim {
            for k in j..dim {
                kappa.set(i, j, k, Rational::from(1));
            }
        }
    }

    assert!(kappa.num_nonzero() > 100);

    let m = vec![1; dim];
    let n = compute_n_matrix(&kappa, &m);

    assert_eq!(n.nrows(), dim);
}
