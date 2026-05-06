//! Integer kernel (nullspace) computation.
//!
//! Computes an integer basis for the nullspace by performing exact
//! rational row-reduction and then clearing denominators.

use malachite::num::arithmetic::traits::Abs;
use malachite::num::conversion::traits::RoundingFrom;
use malachite::rounding_modes::RoundingMode;
use malachite::{Integer, Rational};

use super::matrix_utils::{gcd_integer, lcm_integer};

/// Compute the integer kernel (nullspace) of a matrix.
///
/// Given an m x n matrix A, find a basis for the set of vectors x in Z^n
/// such that A x = 0.
///
/// This performs exact row-reduction over Q to obtain a rational nullspace
/// basis, then scales each basis vector to integers and reduces by GCD.
pub fn integer_kernel(matrix: &[Vec<Integer>]) -> Vec<Vec<Integer>> {
    if matrix.is_empty() {
        return Vec::new();
    }
    let m = matrix.len();
    let n = matrix[0].len();
    if n == 0 {
        return Vec::new();
    }

    let mut mat = to_rational_matrix(matrix, m, n);
    let pivot_cols = rref(&mut mat, m, n);

    if std::env::var_os("CYRUS_DEBUG_KERNEL").is_some() {
        eprintln!("rref: {mat:?}");
        eprintln!("pivot_cols: {pivot_cols:?}");
    }

    let free_cols = free_columns(n, &pivot_cols);
    if free_cols.is_empty() {
        return Vec::new();
    }

    build_kernel_basis(&mat, &pivot_cols, &free_cols, n)
}

fn to_rational_matrix(matrix: &[Vec<Integer>], m: usize, n: usize) -> Vec<Vec<Rational>> {
    let mut mat: Vec<Vec<Rational>> = vec![vec![Rational::from(0); n]; m];
    for i in 0..m {
        for j in 0..n {
            mat[i][j] = Rational::from(&matrix[i][j]);
        }
    }
    mat
}

fn rref(mat: &mut [Vec<Rational>], m: usize, n: usize) -> Vec<usize> {
    let mut pivot_cols: Vec<usize> = Vec::new();
    let mut row = 0usize;
    for col in 0..n {
        if row >= m {
            break;
        }
        let Some(pivot_row) = find_pivot_row(mat, row, col, m) else {
            continue;
        };
        if pivot_row != row {
            mat.swap(pivot_row, row);
        }
        normalize_pivot_row(mat, row, col, n);
        eliminate_column(mat, row, col, m, n);
        pivot_cols.push(col);
        row += 1;
    }
    pivot_cols
}

fn find_pivot_row(mat: &[Vec<Rational>], start_row: usize, col: usize, m: usize) -> Option<usize> {
    (start_row..m).find(|&r| mat[r][col] != 0)
}

fn normalize_pivot_row(mat: &mut [Vec<Rational>], row: usize, col: usize, n: usize) {
    let pivot_val = mat[row][col].clone();
    for c in col..n {
        mat[row][c] /= pivot_val.clone();
    }
}

fn eliminate_column(mat: &mut [Vec<Rational>], row: usize, col: usize, m: usize, n: usize) {
    let pivot_row_vals = mat[row].clone();
    for r in 0..m {
        if r == row {
            continue;
        }
        if mat[r][col] != 0 {
            let factor = mat[r][col].clone();
            for c in col..n {
                mat[r][c] -= factor.clone() * pivot_row_vals[c].clone();
            }
        }
    }
}

fn free_columns(n: usize, pivot_cols: &[usize]) -> Vec<usize> {
    (0..n).filter(|col| !pivot_cols.contains(col)).collect()
}

fn build_kernel_basis(
    mat: &[Vec<Rational>],
    pivot_cols: &[usize],
    free_cols: &[usize],
    n: usize,
) -> Vec<Vec<Integer>> {
    let mut basis: Vec<Vec<Integer>> = Vec::new();
    for &free in free_cols {
        let vec = build_rational_kernel_vector(mat, pivot_cols, free, n);
        if std::env::var_os("CYRUS_DEBUG_KERNEL").is_some() {
            eprintln!("free={free} raw_vec={vec:?}");
        }
        let mut int_vec = scale_vector_to_integers(&vec);
        reduce_vector_gcd(&mut int_vec);
        basis.push(int_vec);
    }
    basis
}

fn build_rational_kernel_vector(
    mat: &[Vec<Rational>],
    pivot_cols: &[usize],
    free: usize,
    n: usize,
) -> Vec<Rational> {
    let mut vec = vec![Rational::from(0); n];
    vec[free] = Rational::from(1);
    for (r, &pcol) in pivot_cols.iter().enumerate() {
        let coeff = mat[r][free].clone();
        vec[pcol] = -coeff;
    }
    vec
}

fn scale_vector_to_integers(vec: &[Rational]) -> Vec<Integer> {
    let mut lcm = Integer::from(1);
    for v in vec {
        let (_, denom) = v.clone().into_numerator_and_denominator();
        lcm = lcm_integer(&lcm, &Integer::from(denom));
    }
    let mut int_vec: Vec<Integer> = Vec::with_capacity(vec.len());
    for v in vec {
        let scaled = v * Rational::from(&lcm);
        let int_val = Integer::rounding_from(&scaled, RoundingMode::Floor).0;
        debug_assert!(
            Rational::from(&int_val) == scaled,
            "Non-integer after scaling: {scaled}"
        );
        int_vec.push(int_val);
    }
    int_vec
}

fn reduce_vector_gcd(int_vec: &mut [Integer]) {
    let mut g = Integer::from(0);
    for v in int_vec.iter() {
        let abs = v.clone().abs();
        if g == 0 {
            g = abs;
        } else if abs != 0 {
            g = gcd_integer(&g, &abs);
        }
    }
    if g > 1 {
        for v in int_vec.iter_mut() {
            *v /= &g;
        }
    }
}
