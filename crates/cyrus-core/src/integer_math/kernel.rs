//! Integer kernel (nullspace) computation.
//!
//! Computes an integer basis for the nullspace by performing exact
//! rational row-reduction and then clearing denominators.

use malachite::{Integer, Rational};
use malachite::num::arithmetic::traits::Abs;
use malachite::num::conversion::traits::RoundingFrom;
use malachite::rounding_modes::RoundingMode;

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

    // Convert to rational matrix for exact elimination.
    let mut mat: Vec<Vec<Rational>> = vec![vec![Rational::from(0); n]; m];
    for i in 0..m {
        for j in 0..n {
            mat[i][j] = Rational::from(&matrix[i][j]);
        }
    }

    // Row-reduced echelon form.
    let mut pivot_cols: Vec<usize> = Vec::new();
    let mut row = 0usize;
    for col in 0..n {
        if row >= m {
            break;
        }

        // Find pivot row.
        let mut pivot_row = None;
        for r in row..m {
            if mat[r][col] != 0 {
                pivot_row = Some(r);
                break;
            }
        }
        let Some(pr) = pivot_row else {
            continue;
        };
        if pr != row {
            mat.swap(pr, row);
        }

        // Normalize pivot row.
        let pivot_val = mat[row][col].clone();
        for c in col..n {
            mat[row][c] /= pivot_val.clone();
        }

        // Eliminate this column from other rows.
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

        pivot_cols.push(col);
        row += 1;
    }

    if std::env::var_os("CYRUS_DEBUG_KERNEL").is_some() {
        eprintln!("rref: {mat:?}");
        eprintln!("pivot_cols: {pivot_cols:?}");
    }

    // Identify free columns.
    let mut free_cols: Vec<usize> = Vec::new();
    for col in 0..n {
        if !pivot_cols.contains(&col) {
            free_cols.push(col);
        }
    }
    if free_cols.is_empty() {
        return Vec::new();
    }

    // Build basis vectors for each free variable.
    let mut basis: Vec<Vec<Integer>> = Vec::new();
    for &free in &free_cols {
        let mut vec = vec![Rational::from(0); n];
        vec[free] = Rational::from(1);
        for (r, &pcol) in pivot_cols.iter().enumerate() {
            let coeff = mat[r][free].clone();
            vec[pcol] = -coeff;
        }
        if std::env::var_os("CYRUS_DEBUG_KERNEL").is_some() {
            eprintln!("free={free} raw_vec={vec:?}");
        }

        // Scale to integers.
        let mut lcm = Integer::from(1);
        for v in &vec {
            let (_, denom) = v.clone().into_numerator_and_denominator();
            lcm = lcm_integer(&lcm, &Integer::from(denom));
        }

        let mut int_vec: Vec<Integer> = Vec::with_capacity(n);
        for v in &vec {
            let scaled = v * Rational::from(&lcm);
            let int_val = Integer::rounding_from(&scaled, RoundingMode::Floor).0;
            debug_assert!(
                Rational::from(&int_val) == scaled,
                "Non-integer after scaling: {scaled}"
            );
            int_vec.push(int_val);
        }

        // Reduce by GCD for a primitive vector.
        let mut g = Integer::from(0);
        for v in &int_vec {
            let abs = v.clone().abs();
            if g == 0 {
                g = abs;
            } else if abs != 0 {
                g = gcd_integer(&g, &abs);
            }
        }
        if g > 1 {
            for v in &mut int_vec {
                *v /= &g;
            }
        }

        basis.push(int_vec);
    }

    basis
}
