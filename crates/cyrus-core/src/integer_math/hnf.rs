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

    let row_count = matrix.len();
    let col_count = matrix[0].len();
    let mut hnf_rows = matrix.to_vec();

    let mut row = 0usize;
    for col in 0..col_count {
        if row >= row_count {
            break;
        }

        let Some(pr) = find_pivot_row(&hnf_rows, row, col, row_count) else {
            continue;
        };
        if pr != row {
            hnf_rows.swap(pr, row);
        }

        clear_below_pivot(&mut hnf_rows, row, col, row_count, col_count);

        if hnf_rows[row][col] == 0 {
            continue;
        }

        ensure_positive_pivot(&mut hnf_rows, row, col, col_count);

        reduce_above_pivot(&mut hnf_rows, row, col, col_count);

        row += 1;
    }

    hnf_rows
}

fn find_pivot_row(
    hnf_rows: &[Vec<Integer>],
    start_row: usize,
    col: usize,
    row_count: usize,
) -> Option<usize> {
    let mut pivot_row: Option<usize> = None;
    let mut pivot_abs: Option<Integer> = None;
    for r in start_row..row_count {
        if hnf_rows[r][col] != 0 {
            let abs = hnf_rows[r][col].clone().abs();
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
    pivot_row
}

fn clear_below_pivot(
    hnf_rows: &mut [Vec<Integer>],
    row: usize,
    col: usize,
    row_count: usize,
    col_count: usize,
) {
    for r in (row + 1)..row_count {
        if hnf_rows[r][col] == 0 {
            continue;
        }
        let (gcd, coeff_u, coeff_v) = extended_gcd_int(&hnf_rows[row][col], &hnf_rows[r][col]);
        let row_i = hnf_rows[row].clone();
        let row_r = hnf_rows[r].clone();
        let pivot_val = hnf_rows[row][col].clone();
        let row_val = hnf_rows[r][col].clone();
        let gcd_val = gcd.clone();
        let factor_i = pivot_val / &gcd_val;
        let factor_r = row_val / &gcd_val;

        for c in 0..col_count {
            hnf_rows[row][c] = &coeff_u * &row_i[c] + &coeff_v * &row_r[c];
            hnf_rows[r][c] = -&factor_r * &row_i[c] + &factor_i * &row_r[c];
        }
    }
}

fn ensure_positive_pivot(hnf_rows: &mut [Vec<Integer>], row: usize, col: usize, col_count: usize) {
    if hnf_rows[row][col] < 0 {
        for c in 0..col_count {
            hnf_rows[row][c] = -hnf_rows[row][c].clone();
        }
    }
}

fn reduce_above_pivot(hnf_rows: &mut [Vec<Integer>], row: usize, col: usize, col_count: usize) {
    let pivot_row = hnf_rows[row].clone();
    for r in 0..row {
        if hnf_rows[r][col] != 0 {
            let q = div_floor(&hnf_rows[r][col], &pivot_row[col]);
            for c in 0..col_count {
                hnf_rows[r][c] -= &q * &pivot_row[c];
            }
        }
    }
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
        if let Some(prev) = last_pivot
            && pivot_col <= prev
        {
            return false;
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

fn extended_gcd_int(lhs: &Integer, rhs: &Integer) -> (Integer, Integer, Integer) {
    let (gcd, coeff_x, coeff_y) = lhs.extended_gcd(rhs);
    (Integer::from(gcd), coeff_x, coeff_y)
}
