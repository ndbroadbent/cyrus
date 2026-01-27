//! Hermite normal form utilities for integer matrices.
//!
//! Implements a basic row Hermite normal form (HNF) computation using
//! integer row operations. This is sufficient for extracting pivot
//! columns and sublattice indices in the CYTools GLSM algorithms.

use malachite::Integer;
use malachite::num::arithmetic::traits::{Abs, ExtendedGcd};

/// Compute the row Hermite normal form of an integer matrix.
///
/// The result has the following properties:
/// - Each non-zero row has a leftmost non-zero pivot.
/// - Pivot columns are strictly increasing down the rows.
/// - Pivot entries are positive.
/// - Entries above each pivot are reduced into the range [0, pivot).
pub fn hermite_normal_form(matrix: &[Vec<Integer>]) -> Vec<Vec<Integer>> {
    if matrix.is_empty() {
        return Vec::new();
    }

    let m = matrix.len();
    let n = matrix[0].len();
    let mut h = matrix.to_vec();

    let mut row = 0usize;
    for col in 0..n {
        if row >= m {
            break;
        }

        // Find a pivot row with minimal absolute value in this column.
        let mut pivot_row: Option<usize> = None;
        let mut pivot_abs: Option<Integer> = None;
        for r in row..m {
            if h[r][col] != 0 {
                let abs = h[r][col].clone().abs();
                let better = match &pivot_abs {
                    Some(best) => abs < *best,
                    None => true,
                };
                if better {
                    pivot_abs = Some(abs);
                    pivot_row = Some(r);
                }
            }
        }

        let Some(pr) = pivot_row else {
            continue;
        };
        if pr != row {
            h.swap(pr, row);
        }

        // Clear entries *below* the pivot using extended gcd row operations.
        // We must not disturb previously established pivots above `row`.
        for r in (row + 1)..m {
            if h[r][col] == 0 {
                continue;
            }
            let (g, u, v) = extended_gcd_int(&h[row][col], &h[r][col]);
            let row_i = h[row].clone();
            let row_r = h[r].clone();
            let ai = h[row][col].clone();
            let ar = h[r][col].clone();
            let g_i = g.clone();
            let factor_i = ai / &g_i;
            let factor_r = ar / &g_i;

            for c in 0..n {
                h[row][c] = &u * &row_i[c] + &v * &row_r[c];
                h[r][c] = -&factor_r * &row_i[c] + &factor_i * &row_r[c];
            }
        }

        if h[row][col] == 0 {
            continue;
        }

        // Ensure pivot is positive.
        if h[row][col] < 0 {
            for c in 0..n {
                h[row][c] = -h[row][c].clone();
            }
        }

        // Reduce entries above the pivot into [0, pivot).
        for r in 0..row {
            if h[r][col] != 0 {
                let q = div_floor(&h[r][col], &h[row][col]);
                let pivot_row = h[row].clone();
                for c in 0..n {
                    h[r][c] -= &q * &pivot_row[c];
                }
            }
        }

        row += 1;
    }

    h
}

/// Validate that a matrix is in row Hermite normal form.
pub fn is_row_hnf(matrix: &[Vec<Integer>]) -> bool {
    if matrix.is_empty() {
        return true;
    }
    let n_cols = matrix[0].len();
    let mut last_pivot: Option<usize> = None;

    for (r, row) in matrix.iter().enumerate() {
        if row.len() != n_cols {
            return false;
        }
        let pivot_col = row.iter().position(|v| *v != 0);
        let Some(pivot_col) = pivot_col else {
            continue;
        };
        if let Some(prev) = last_pivot {
            if pivot_col <= prev {
                return false;
            }
        }
        let pivot = &row[pivot_col];
        if *pivot <= 0 {
            return false;
        }
        if row[..pivot_col].iter().any(|v| *v != 0) {
            return false;
        }
        for rr in 0..r {
            let above = &matrix[rr][pivot_col];
            if *above < 0 || *above >= *pivot {
                return false;
            }
        }
        last_pivot = Some(pivot_col);
    }

    true
}

/// Compute the column Hermite normal form of an integer matrix.
///
/// This is equivalent to computing the row HNF of the transpose and
/// transposing the result back.
pub fn hermite_normal_form_columns(matrix: &[Vec<Integer>]) -> Vec<Vec<Integer>> {
    if matrix.is_empty() {
        return Vec::new();
    }

    let rows = matrix.len();
    let cols = matrix[0].len();
    let mut transposed = vec![vec![Integer::from(0); rows]; cols];
    for r in 0..rows {
        for c in 0..cols {
            transposed[c][r] = matrix[r][c].clone();
        }
    }

    let hnf_t = hermite_normal_form(&transposed);

    let mut out = vec![vec![Integer::from(0); cols]; rows];
    for r in 0..rows {
        for c in 0..cols {
            out[r][c] = hnf_t[c][r].clone();
        }
    }

    out
}

/// Compute the product of the pivot entries in a row-echelon integer matrix.
pub fn pivot_product(matrix: &[Vec<Integer>]) -> Integer {
    let mut prod = Integer::from(1);
    for row in matrix {
        if let Some(val) = row.iter().find(|v| **v != 0) {
            prod *= val.clone().abs();
        }
    }
    prod
}

/// Compute the sublattice index using the HNF pivot product.
pub fn sublattice_index(matrix: &[Vec<Integer>]) -> Integer {
    let hnf = hermite_normal_form(matrix);
    pivot_product(&hnf)
}

fn div_floor(a: &Integer, b: &Integer) -> Integer {
    let q = a / b;
    let r = a - &q * b;
    if r != 0 && r < 0 {
        q - Integer::from(1)
    } else {
        q
    }
}

fn extended_gcd_int(a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
    let (g, x, y) = a.extended_gcd(b);
    (Integer::from(g), x, y)
}
