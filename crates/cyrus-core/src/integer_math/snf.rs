//! Smith normal form utilities for integer matrices.
//!
//! This implementation targets small matrices and is used to compute
//! sublattice indices (product of invariant factors).

use malachite::Integer;
use malachite::num::arithmetic::traits::{Abs, ExtendedGcd};

use super::hnf::{hermite_normal_form, pivot_product};

/// Compute the diagonal entries of the Smith normal form.
pub fn smith_normal_form_diag(matrix: &[Vec<Integer>]) -> Vec<Integer> {
    if matrix.is_empty() {
        return Vec::new();
    }
    let row_count = matrix.len();
    let col_count = matrix[0].len();
    let mut work = matrix.to_vec();

    let mut row_idx = 0usize;
    let mut col_idx = 0usize;
    while row_idx < row_count && col_idx < col_count {
        let Some((pr, pc)) = find_pivot(&work, row_idx, col_idx, row_count, col_count) else {
            break;
        };
        if pr != row_idx {
            work.swap(pr, row_idx);
        }
        if pc != col_idx {
            swap_cols(&mut work, pc, col_idx);
        }

        loop {
            clear_column(&mut work, row_idx, col_idx, row_count, col_count);
            clear_row(&mut work, row_idx, col_idx, row_count, col_count);
            ensure_positive_pivot(&mut work, row_idx, col_idx, col_count);

            if let Some(bad_row) = find_non_divisible(&work, row_idx, col_idx, row_count, col_count)
            {
                adjust_pivot_with_row(&mut work, row_idx, bad_row, col_count);
                continue;
            }
            break;
        }

        row_idx += 1;
        col_idx += 1;
    }

    let diag_len = row_idx.min(col_idx);
    let mut diag = Vec::with_capacity(diag_len);
    for k in 0..diag_len {
        diag.push(work[k][k].clone().abs());
    }
    diag
}

fn find_pivot(
    work: &[Vec<Integer>],
    row_idx: usize,
    col_idx: usize,
    row_count: usize,
    col_count: usize,
) -> Option<(usize, usize)> {
    let mut pivot: Option<(usize, usize, Integer)> = None;
    for r in row_idx..row_count {
        for c in col_idx..col_count {
            if work[r][c] != 0 {
                let abs = work[r][c].clone().abs();
                let better = match &pivot {
                    Some((_, _, best)) => abs < *best,
                    None => true,
                };
                if better {
                    pivot = Some((r, c, abs));
                }
            }
        }
    }
    pivot.map(|(r, c, _)| (r, c))
}

fn clear_column(
    work: &mut [Vec<Integer>],
    row_idx: usize,
    col_idx: usize,
    row_count: usize,
    col_count: usize,
) {
    for r in 0..row_count {
        if r == row_idx || work[r][col_idx] == 0 {
            continue;
        }
        let (gcd, coeff_u, coeff_v) = extended_gcd_int(&work[row_idx][col_idx], &work[r][col_idx]);
        let row_i = work[row_idx].clone();
        let row_r = work[r].clone();
        let pivot_val = work[row_idx][col_idx].clone();
        let row_val = work[r][col_idx].clone();
        let gcd_val = gcd.clone();
        let factor_i = pivot_val / &gcd_val;
        let factor_r = row_val / &gcd_val;

        for c in 0..col_count {
            work[row_idx][c] = &coeff_u * &row_i[c] + &coeff_v * &row_r[c];
            work[r][c] = -&factor_r * &row_i[c] + &factor_i * &row_r[c];
        }
    }
}

fn clear_row(
    work: &mut [Vec<Integer>],
    row_idx: usize,
    col_idx: usize,
    row_count: usize,
    col_count: usize,
) {
    for c in 0..col_count {
        if c == col_idx || work[row_idx][c] == 0 {
            continue;
        }
        let (gcd, coeff_u, coeff_v) = extended_gcd_int(&work[row_idx][col_idx], &work[row_idx][c]);
        let col_pivot = get_col(work, col_idx);
        let col_target = get_col(work, c);
        let pivot_val = work[row_idx][col_idx].clone();
        let target_val = work[row_idx][c].clone();
        let gcd_val = gcd.clone();
        let factor_i = pivot_val / &gcd_val;
        let factor_c = target_val / &gcd_val;

        for r in 0..row_count {
            work[r][col_idx] = &coeff_u * &col_pivot[r] + &coeff_v * &col_target[r];
            work[r][c] = -&factor_c * &col_pivot[r] + &factor_i * &col_target[r];
        }
    }
}

fn ensure_positive_pivot(
    work: &mut [Vec<Integer>],
    row_idx: usize,
    col_idx: usize,
    col_count: usize,
) {
    if work[row_idx][col_idx] < 0 {
        for c in 0..col_count {
            work[row_idx][c] = -work[row_idx][c].clone();
        }
    }
}

fn find_non_divisible(
    work: &[Vec<Integer>],
    row_idx: usize,
    col_idx: usize,
    row_count: usize,
    col_count: usize,
) -> Option<usize> {
    for r in row_idx..row_count {
        for c in col_idx..col_count {
            if r == row_idx && c == col_idx {
                continue;
            }
            if &work[r][c] % &work[row_idx][col_idx] != 0 {
                return Some(r);
            }
        }
    }
    None
}

fn adjust_pivot_with_row(
    work: &mut [Vec<Integer>],
    row_idx: usize,
    bad_row: usize,
    col_count: usize,
) {
    let row_r = work[bad_row].clone();
    for col in 0..col_count {
        work[row_idx][col] += row_r[col].clone();
    }
}

/// Compute the sublattice index from Smith normal form diagonal entries.
pub fn sublattice_index_snf(matrix: &[Vec<Integer>]) -> Integer {
    if matrix.is_empty() {
        return Integer::from(1);
    }

    // The index of the column lattice equals the index of the row lattice of
    // the transpose. A row HNF provides a lattice basis whose pivot product
    // equals the sublattice index. This matches the SNF invariant product
    // used by CYTools without needing the full SNF.
    let rows = matrix.len();
    let cols = matrix[0].len();
    let mut transposed = vec![vec![Integer::from(0); rows]; cols];
    for r in 0..rows {
        for c in 0..cols {
            transposed[c][r] = matrix[r][c].clone();
        }
    }

    let hnf = hermite_normal_form(&transposed);
    pivot_product(&hnf)
}

fn extended_gcd_int(lhs: &Integer, rhs: &Integer) -> (Integer, Integer, Integer) {
    let (gcd, coeff_x, coeff_y) = lhs.extended_gcd(rhs);
    let gcd_val = Integer::from(gcd);
    (gcd_val, coeff_x, coeff_y)
}

fn swap_cols(matrix: &mut [Vec<Integer>], c1: usize, c2: usize) {
    for row in matrix {
        row.swap(c1, c2);
    }
}

fn get_col(matrix: &[Vec<Integer>], column: usize) -> Vec<Integer> {
    matrix.iter().map(|row| row[column].clone()).collect()
}
