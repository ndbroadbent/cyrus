#![allow(missing_docs)]
use cyrus_core::integer_math::{integer_kernel, orientation, solve_linear_system_rational};
use malachite::{Integer, Rational};

fn i(n: i64) -> Integer {
    Integer::from(n)
}

fn r(n: i64) -> Rational {
    Rational::from(n)
}

#[test]
fn test_integer_kernel_empty() {
    let empty: Vec<Vec<Integer>> = Vec::new();
    let kernel = integer_kernel(&empty);
    assert!(kernel.is_empty());
}

#[test]
fn test_integer_kernel_zero() {
    // A = [0, 0, 0] (1x3)
    // Kernel is whole Z^3?
    // Algorithm augments with identity.
    let a = vec![vec![i(0), i(0), i(0)]];
    let kernel = integer_kernel(&a);
    // Should return basis for Z^3 (3 vectors)
    assert_eq!(kernel.len(), 3);
}

#[test]
fn test_orientation_zero_dim() {
    // Just 1 point: dim=0.
    let points = vec![vec![r(0)]];
    let vol = orientation(&points);
    assert_eq!(vol, r(0));
}

#[test]
fn test_solve_rational_singular() {
    // Singular system [1, 1; 1, 1] x = [1, 2]
    let m = vec![vec![r(1), r(1)], vec![r(1), r(1)]];
    let c = vec![r(1), r(2)];

    let res = solve_linear_system_rational(&m, &c);
    assert!(res.is_none());
}

#[test]
fn test_solve_rational_empty() {
    let m: Vec<Vec<Rational>> = Vec::new();
    let c: Vec<Rational> = Vec::new();
    let res = solve_linear_system_rational(&m, &c);
    assert_eq!(res, Some(vec![]));
}
