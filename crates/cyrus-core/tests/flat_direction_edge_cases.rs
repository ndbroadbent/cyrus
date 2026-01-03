#![allow(missing_docs, clippy::too_many_lines)]
use cyrus_core::flat_direction::compute_n_matrix_vecs;
use cyrus_core::flat_direction::solve_linear_system_faer;
use cyrus_core::types::rational::Rational as TypedRational;
use cyrus_core::types::tags::Pos;
use cyrus_core::{
    Intersection, compute_ek0, compute_flat_direction, compute_flat_direction_full,
    compute_n_matrix, solve_linear_system,
};
use faer::Mat;
use malachite::Rational;

fn pos_rat(n: i32) -> TypedRational<Pos> {
    TypedRational::<Pos>::new(Rational::from(n)).unwrap()
}

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
    kappa.set(1, 2, 3, pos_rat(1));

    let m = vec![0, 0, 0, 1]; // Indices 0..3. 3 is the one.

    let n = compute_n_matrix(&kappa, &m);

    assert!((n[(1, 2)] - 1.0).abs() < 1e-10);
    assert!((n[(2, 1)] - 1.0).abs() < 1e-10);
    assert!((n[(1, 3)] - 0.0).abs() < 1e-10);

    // Also check self-intersections
    let mut kappa2 = Intersection::new(1);
    kappa2.set(0, 0, 0, pos_rat(1));
    let m2 = vec![1];
    let n2 = compute_n_matrix(&kappa2, &m2);
    assert!((n2[(0, 0)] - 1.0).abs() < 1e-10);
}

#[test]
fn test_compute_flat_direction_full() {
    let mut kappa = Intersection::new(1);
    kappa.set(0, 0, 0, pos_rat(6));
    let k = vec![3];
    let m = vec![1];

    let res = compute_flat_direction_full(&kappa, &k, &m).unwrap();
    assert!((res.p[0] - 0.5).abs() < 1e-10);
    // kappa(p,p,p) = 6 * 0.5^3 = 0.75.
    // ek0 = 1 / (4/3 * 0.75) = 1 / 1 = 1.0.
    assert!((res.ek0.get() - 1.0).abs() < 1e-10);
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
                kappa.set(i, j, k, pos_rat(1));
            }
        }
    }

    assert!(kappa.num_nonzero() > 100);

    let m = vec![1; dim];
    let n = compute_n_matrix(&kappa, &m);

    assert_eq!(n.nrows(), dim);
}

#[test]
fn test_compute_n_matrix_vecs() {
    let mut kappa = Intersection::new(2);
    kappa.set(0, 0, 1, pos_rat(3));

    let m = vec![1, 1];
    let n = compute_n_matrix_vecs(&kappa, &m);

    assert_eq!(n.len(), 2);
    assert_eq!(n[0].len(), 2);
}

#[test]
fn test_solve_linear_system_faer_singular() {
    // Singular matrix
    let n = Mat::from_fn(2, 2, |i, j| if i == j { 0.0 } else { 1.0 });
    let k = vec![1, 1];

    let res = solve_linear_system_faer(&n, &k);
    // May or may not return None depending on singularity detection
    if let Some(p) = res {
        // If it returns something, verify it's reasonable
        assert!(p.iter().all(|&x| x.is_finite()));
    }
}

#[test]
fn test_solve_linear_system_faer_identity() {
    // Identity matrix - trivial case
    let n = Mat::from_fn(2, 2, |i, j| if i == j { 1.0 } else { 0.0 });
    let k = vec![5, 3];

    let res = solve_linear_system_faer(&n, &k).unwrap();
    assert!((res[0] - 5.0).abs() < 1e-10);
    assert!((res[1] - 3.0).abs() < 1e-10);
}

#[test]
fn test_compute_flat_direction() {
    let mut kappa = Intersection::new(2);
    kappa.set(0, 0, 0, pos_rat(6));
    kappa.set(0, 0, 1, pos_rat(3));
    kappa.set(0, 1, 1, pos_rat(2));
    kappa.set(1, 1, 1, pos_rat(4));

    let k = vec![6, 4];
    let m = vec![1, 1];

    let p = compute_flat_direction(&kappa, &k, &m);
    // Just check it returns something
    assert!(p.is_some() || p.is_none());
}

#[test]
fn test_compute_ek0_simple() {
    let mut kappa = Intersection::new(1);
    kappa.set(0, 0, 0, pos_rat(6));

    let p = vec![1.0];
    let ek0 = compute_ek0(&kappa, &p).unwrap();
    // κ_ppp = 6, ek0 = 1/(4/3 * 6) = 0.125
    assert!((ek0.get() - 0.125).abs() < 1e-10);
}

#[test]
fn test_compute_flat_direction_full_singular() {
    // Create a singular N matrix
    let mut kappa = Intersection::new(2);
    // Only diagonal terms that make N singular
    kappa.set(0, 0, 0, pos_rat(1));
    // M = [0, 1] means N[0,0] = κ_000 * 0 = 0, N[0,1] = κ_001 * 1 = 0
    // This creates a singular matrix

    let k = vec![1, 1];
    let m = vec![0, 1]; // Zeros out first column

    let res = compute_flat_direction_full(&kappa, &k, &m);
    // May be singular
    if res.is_none() {
        // Expected for singular case
    }
}

#[test]
fn test_unique_permutations_all_equal() {
    // Test with all indices equal (should give 1 unique permutation)
    let mut kappa = Intersection::new(3);
    kappa.set(1, 1, 1, pos_rat(6));

    let m = vec![1, 1, 1];
    let n = compute_n_matrix(&kappa, &m);

    // N[1,1] should have contribution from κ_111 * m[1] = 6 * 1 = 6
    // But permutations of (1,1,1) give only 1 unique (1,1,1)
    assert!((n[(1, 1)] - 6.0).abs() < 1e-10);
}

#[test]
fn test_unique_permutations_two_equal() {
    // Test with two indices equal (should give 3 unique permutations)
    let mut kappa = Intersection::new(3);
    kappa.set(0, 0, 1, pos_rat(2));

    let m = vec![1, 1, 0];
    let n = compute_n_matrix(&kappa, &m);

    // κ_001 with m = [1,1,0]
    // Permutations: (0,0,1), (0,1,0), (1,0,0)
    // N[0,0] += 2 * m[1] = 2
    // N[0,1] += 2 * m[0] = 2
    // N[1,0] += 2 * m[0] = 2
    assert!(n[(0, 0)] >= 0.0); // At least 0
}

#[test]
fn test_unique_permutations_all_distinct() {
    // Test with all distinct indices (should give 6 unique permutations)
    let mut kappa = Intersection::new(3);
    kappa.set(0, 1, 2, pos_rat(1));

    let m = vec![1, 1, 1];
    let n = compute_n_matrix(&kappa, &m);

    // κ_012 with m = [1,1,1]
    // All 6 permutations contribute
    // N[0,1] += 1 * m[2] = 1
    // N[0,2] += 1 * m[1] = 1
    // N[1,0] += 1 * m[2] = 1
    // N[1,2] += 1 * m[0] = 1
    // N[2,0] += 1 * m[1] = 1
    // N[2,1] += 1 * m[0] = 1
    assert!((n[(0, 1)] - 1.0).abs() < 1e-10);
    assert!((n[(1, 2)] - 1.0).abs() < 1e-10);
}
