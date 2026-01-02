#![allow(missing_docs, clippy::too_many_lines)]
use cyrus_core::integer_math::{
    determinant_gaussian, integer_kernel, orientation, solve_linear_system_rational,
};
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

#[test]
fn test_determinant_empty_matrix() {
    let mut mat: Vec<Vec<Rational>> = Vec::new();
    let det = determinant_gaussian(&mut mat);
    assert_eq!(det, r(1)); // Convention: det of 0x0 is 1
}

#[test]
fn test_determinant_singular_zero_column() {
    // Matrix with zero column - singular
    let mut mat = vec![vec![r(1), r(0)], vec![r(2), r(0)]];
    let det = determinant_gaussian(&mut mat);
    assert_eq!(det, r(0));
}

#[test]
fn test_determinant_requires_pivot_swap() {
    // Matrix that needs row swap to find pivot
    let mut mat = vec![vec![r(0), r(1)], vec![r(1), r(0)]];
    let det = determinant_gaussian(&mut mat);
    // After swap: [[1,0],[0,1]] -> det = -1 (due to swap)
    assert_eq!(det, r(-1));
}

#[test]
fn test_determinant_3x3() {
    // 3x3 matrix
    let mut mat = vec![
        vec![r(1), r(2), r(3)],
        vec![r(4), r(5), r(6)],
        vec![r(7), r(8), r(9)],
    ];
    let det = determinant_gaussian(&mut mat);
    // This matrix is singular (rows are linearly dependent)
    assert_eq!(det, r(0));
}

#[test]
fn test_determinant_3x3_nonsingular() {
    // Non-singular 3x3 matrix
    let mut mat = vec![
        vec![r(1), r(0), r(0)],
        vec![r(0), r(2), r(0)],
        vec![r(0), r(0), r(3)],
    ];
    let det = determinant_gaussian(&mut mat);
    assert_eq!(det, r(6)); // 1 * 2 * 3 = 6
}

#[test]
fn test_orientation_2d() {
    // 3 points in 2D - counter-clockwise triangle
    let points = vec![vec![r(0), r(0)], vec![r(1), r(0)], vec![r(0), r(1)]];
    let vol = orientation(&points);
    // Should be positive (counter-clockwise)
    assert!(vol > r(0));
}

#[test]
fn test_orientation_2d_clockwise() {
    // 3 points in 2D - clockwise triangle
    let points = vec![vec![r(0), r(0)], vec![r(0), r(1)], vec![r(1), r(0)]];
    let vol = orientation(&points);
    // Should be negative (clockwise)
    assert!(vol < r(0));
}

#[test]
fn test_orientation_collinear() {
    // 3 collinear points in 2D
    let points = vec![vec![r(0), r(0)], vec![r(1), r(1)], vec![r(2), r(2)]];
    let vol = orientation(&points);
    // Should be zero (collinear)
    assert_eq!(vol, r(0));
}

#[test]
fn test_integer_kernel_full_rank() {
    // Full rank square matrix - kernel should be empty
    let a = vec![
        vec![i(1), i(0), i(0)],
        vec![i(0), i(1), i(0)],
        vec![i(0), i(0), i(1)],
    ];
    let kernel = integer_kernel(&a);
    assert!(kernel.is_empty());
}

#[test]
fn test_integer_kernel_rank_deficient() {
    // Rank-1 matrix (all rows same)
    let a = vec![vec![i(1), i(2), i(3)], vec![i(2), i(4), i(6)]];
    let kernel = integer_kernel(&a);
    // Kernel should be 2-dimensional (3 vars, rank 1)
    assert_eq!(kernel.len(), 2);
}

#[test]
fn test_integer_kernel_with_negative() {
    // Matrix with negative entries
    let a = vec![vec![i(1), i(-1)]];
    let kernel = integer_kernel(&a);
    // Kernel: x - y = 0 => x = y => [1, 1] is in kernel
    assert_eq!(kernel.len(), 1);
    // Check kernel vector: A * v = 0 => 1*v[0] + (-1)*v[1] = 0
    let v = &kernel[0];
    assert_eq!(&v[0] - &v[1], i(0));
}

#[test]
fn test_solve_rational_with_pivot_swap() {
    // Matrix where first column has zero in first row
    let m = vec![vec![r(0), r(1)], vec![r(1), r(0)]];
    let c = vec![r(2), r(3)];
    let res = solve_linear_system_rational(&m, &c).unwrap();
    // x = 3, y = 2
    assert_eq!(res[0], r(3));
    assert_eq!(res[1], r(2));
}

#[test]
fn test_solve_rational_3x3() {
    // 3x3 system
    let m = vec![
        vec![r(2), r(1), r(1)],
        vec![r(1), r(3), r(2)],
        vec![r(1), r(0), r(2)],
    ];
    let c = vec![r(4), r(6), r(3)];
    let res = solve_linear_system_rational(&m, &c);
    assert!(res.is_some());
    // Verify solution satisfies equations
    let x = res.unwrap();
    for (row, expected) in m.iter().zip(c.iter()) {
        let mut sum = r(0);
        for (coef, xi) in row.iter().zip(x.iter()) {
            sum += coef * xi;
        }
        assert_eq!(sum, *expected);
    }
}

#[test]
fn test_integer_kernel_wide_matrix() {
    // Matrix with more columns than rows (1x5)
    // This triggers the pivot_row >= n break (line 110)
    let a = vec![vec![i(1), i(2), i(3), i(4), i(5)]];
    let kernel = integer_kernel(&a);
    // Should have 4 kernel vectors (5 columns - 1 row = 4 dim kernel)
    assert_eq!(kernel.len(), 4);
    // Verify each kernel vector satisfies A*v = 0
    for v in &kernel {
        let mut sum = i(0);
        for (coef, vi) in a[0].iter().zip(v.iter()) {
            sum += coef * vi;
        }
        assert_eq!(sum, i(0));
    }
}

#[test]
fn test_integer_kernel_tall_matrix() {
    // Matrix with more rows than columns (3x2)
    let a = vec![
        vec![i(1), i(2)],
        vec![i(2), i(4)], // Multiple of first row
        vec![i(3), i(6)], // Also multiple
    ];
    let kernel = integer_kernel(&a);
    // All rows are multiples of [1, 2], so kernel is 1-dimensional
    assert!(!kernel.is_empty());
}
