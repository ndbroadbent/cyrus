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
    let m = matrix.len();
    let n = matrix[0].len();
    let mut a = matrix.to_vec();

    let mut i = 0usize;
    let mut j = 0usize;
    while i < m && j < n {
        // Find smallest non-zero entry in submatrix.
        let mut pivot: Option<(usize, usize, Integer)> = None;
        for r in i..m {
            for c in j..n {
                if a[r][c] != 0 {
                    let abs = a[r][c].clone().abs();
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
        let Some((pr, pc, _)) = pivot else {
            break;
        };
        if pr != i {
            a.swap(pr, i);
        }
        if pc != j {
            swap_cols(&mut a, pc, j);
        }

        loop {
            // Clear column j.
            for r in 0..m {
                if r == i || a[r][j] == 0 {
                    continue;
                }
                let (g, u, v) = extended_gcd_int(&a[i][j], &a[r][j]);
                let row_i = a[i].clone();
                let row_r = a[r].clone();
                let ai = a[i][j].clone();
                let ar = a[r][j].clone();
                let g_i = g.clone();
                let factor_i = ai / &g_i;
                let factor_r = ar / &g_i;

                for c in 0..n {
                    a[i][c] = &u * &row_i[c] + &v * &row_r[c];
                    a[r][c] = -&factor_r * &row_i[c] + &factor_i * &row_r[c];
                }
            }

            // Clear row i.
            for c in 0..n {
                if c == j || a[i][c] == 0 {
                    continue;
                }
                let (g, u, v) = extended_gcd_int(&a[i][j], &a[i][c]);
                let col_j = get_col(&a, j);
                let col_c = get_col(&a, c);
                let ai = a[i][j].clone();
                let ac = a[i][c].clone();
                let g_i = g.clone();
                let factor_i = ai / &g_i;
                let factor_c = ac / &g_i;

                for r in 0..m {
                    a[r][j] = &u * &col_j[r] + &v * &col_c[r];
                    a[r][c] = -&factor_c * &col_j[r] + &factor_i * &col_c[r];
                }
            }

            if a[i][j] < 0 {
                for c in 0..n {
                    a[i][c] = -a[i][c].clone();
                }
            }

            // Check divisibility of submatrix by pivot.
            let mut bad: Option<(usize, usize)> = None;
            for r in i..m {
                for c in j..n {
                    if r == i && c == j {
                        continue;
                    }
                    if &a[r][c] % &a[i][j] != 0 {
                        bad = Some((r, c));
                        break;
                    }
                }
                if bad.is_some() {
                    break;
                }
            }
            if let Some((r, _c)) = bad {
                // Adjust pivot by adding offending row.
                let row_r = a[r].clone();
                for col in 0..n {
                    a[i][col] += row_r[col].clone();
                }
                continue;
            }

            break;
        }

        i += 1;
        j += 1;
    }

    let diag_len = i.min(j);
    let mut diag = Vec::with_capacity(diag_len);
    for k in 0..diag_len {
        diag.push(a[k][k].clone().abs());
    }
    diag
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

fn extended_gcd_int(a: &Integer, b: &Integer) -> (Integer, Integer, Integer) {
    let (g, x, y) = a.extended_gcd(b);
    let g_i = Integer::from(g);
    (g_i, x, y)
}

fn swap_cols(m: &mut [Vec<Integer>], c1: usize, c2: usize) {
    for row in m {
        row.swap(c1, c2);
    }
}

fn get_col(m: &[Vec<Integer>], c: usize) -> Vec<Integer> {
    m.iter().map(|row| row[c].clone()).collect()
}
