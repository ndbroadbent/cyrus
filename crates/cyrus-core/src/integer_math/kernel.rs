//! Integer kernel (nullspace) computation.
//!
//! Provides algorithms for computing the integer kernel of matrices using
//! Hermite Normal Form (HNF) approach.

use malachite::Integer;
use malachite::num::arithmetic::traits::Abs;

/// Compute the integer kernel (nullspace) of a matrix.
///
/// Given an m x n matrix A, find a basis for the set of vectors x in Z^n
/// such that A x = 0.
///
/// # Algorithm
/// Uses the Hermite Normal Form (HNF) approach.
/// We compute the HNF of the augmented matrix [A^T | I].
///
/// Reference: Cohen, "A Course in Computational Algebraic Number Theory", Section 2.4.
///
/// # Panics
/// Panics if the internal state becomes inconsistent during pivot selection.
pub fn integer_kernel(matrix: &[Vec<Integer>]) -> Vec<Vec<Integer>> {
    if matrix.is_empty() {
        return Vec::new();
    }

    let m = matrix.len();
    let n = matrix[0].len();

    // 1. Construct the matrix M = [A^T | I] of size n x (m + n)
    let mut m_mat = vec![vec![Integer::from(0); m + n]; n];
    for (i, row) in m_mat.iter_mut().enumerate().take(n) {
        for (j, col) in matrix.iter().enumerate().take(m) {
            row[j] = col[i].clone();
        }
        row[m + i] = Integer::from(1);
    }

    // 2. Perform row operations to zero out the A^T part (first m columns)
    let mut pivot_row = 0;
    for col in 0..m {
        if pivot_row >= n {
            break;
        }

        // Find a row with a non-zero entry in this column
        let mut best_row: Option<usize> = None;
        for r in pivot_row..n {
            let val = &m_mat[r][col];
            if *val == 0 {
                continue;
            }

            if let Some(curr_best) = best_row {
                if val.clone().abs() < m_mat[curr_best][col].clone().abs() {
                    best_row = Some(r);
                }
            } else {
                best_row = Some(r);
            }
        }

        if let Some(r) = best_row {
            m_mat.swap(pivot_row, r);
            reduce_column(&mut m_mat, pivot_row, col, n, m);
            pivot_row += 1;
        }
    }

    // 3. The rows of M where the first m columns are all zero form the kernel
    collect_kernel(&m_mat, n, m)
}

fn reduce_column(m_mat: &mut [Vec<Integer>], pivot_row: usize, col: usize, n: usize, m: usize) {
    loop {
        let mut finished = true;
        for r in 0..n {
            if r == pivot_row || m_mat[r][col] == 0 {
                continue;
            }

            let q = &m_mat[r][col] / &m_mat[pivot_row][col];
            if q != 0 {
                for c in col..m + n {
                    let sub = &q * &m_mat[pivot_row][c];
                    m_mat[r][c] -= sub;
                }
            }

            if m_mat[r][col] != 0 {
                m_mat.swap(pivot_row, r);
                finished = false;
                break;
            }
        }
        if finished {
            break;
        }
    }
}

fn collect_kernel(m_mat: &[Vec<Integer>], n: usize, m: usize) -> Vec<Vec<Integer>> {
    let mut kernel = Vec::new();
    for row in m_mat.iter().take(n) {
        let is_zero = row.iter().take(m).all(|x| *x == 0);
        if is_zero {
            let vec: Vec<Integer> = row[m..m + n].to_vec();
            if vec.iter().any(|x| *x != 0) {
                kernel.push(vec);
            }
        }
    }
    kernel
}
