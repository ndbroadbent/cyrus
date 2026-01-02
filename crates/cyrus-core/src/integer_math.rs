//! Integer linear algebra primitives.
//!
//! Provides algorithms for working with integer matrices, such as
//! Hermite Normal Form (HNF) and integer kernel computation.
//!
//! These are essential for finding the integral basis of the GLSM charge matrix.

use malachite::Integer;
use malachite::Rational;
use malachite::num::arithmetic::traits::Abs;

/// Compute the orientation of d+1 points in d-dimensional space.
///
/// Returns the signed volume (times d!).
/// Positive if points are in counter-clockwise order.
pub fn orientation(points: &[Vec<Rational>]) -> Rational {
    let dim = points.len() - 1;
    if dim == 0 {
        return Rational::from(0);
    }

    // Construct matrix of differences
    let mut mat = vec![vec![Rational::from(0); dim]; dim];
    for (i, row) in mat.iter_mut().enumerate().take(dim) {
        for (j, item) in row.iter_mut().enumerate().take(dim) {
            *item = &points[i + 1][j] - &points[0][j];
        }
    }

    determinant_gaussian(&mut mat)
}

/// Compute determinant of a square matrix of Rationals using Gaussian elimination.
pub fn determinant_gaussian(mat: &mut [Vec<Rational>]) -> Rational {
    let n = mat.len();
    if n == 0 {
        return Rational::from(1);
    }

    let mut det = Rational::from(1);

    for i in 0..n {
        // Find pivot
        let mut pivot = i;
        while pivot < n && mat[pivot][i] == 0 {
            pivot += 1;
        }

        if pivot == n {
            return Rational::from(0);
        }

        if pivot != i {
            mat.swap(i, pivot);
            det = -det;
        }

        let diag = mat[i][i].clone();
        det *= &diag;

        // Eliminate rows below
        for j in (i + 1)..n {
            if mat[j][i] != 0 {
                let factor = &mat[j][i] / &diag;
                for k in i..n {
                    let val = &factor * &mat[i][k];
                    mat[j][k] -= val;
                }
            }
        }
    }

    det
}

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
            if *val != 0
                && (best_row.is_none()
                    || val.clone().abs() < m_mat[best_row.unwrap()][col].clone().abs())
            {
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

/// Solve a linear system Mx = C where M is square and non-singular.
///
/// Uses Gaussian elimination with pivoting.
/// Returns None if singular.
///
/// # Panics
/// Panics if matrix M is not square or if vector C dimension does not match M.
pub fn solve_linear_system_rational(m: &[Vec<Rational>], c: &[Rational]) -> Option<Vec<Rational>> {
    let n = m.len();
    if n == 0 {
        return Some(Vec::new());
    }
    assert_eq!(m[0].len(), n, "Matrix M must be square");
    assert_eq!(c.len(), n, "Vector C must have dimension n");

    // Augmented matrix [M | C]
    let mut mat = vec![vec![Rational::from(0); n + 1]; n];
    for (i, row) in mat.iter_mut().enumerate().take(n) {
        for (j, item) in row.iter_mut().enumerate().take(n) {
            *item = m[i][j].clone();
        }
        row[n] = c[i].clone();
    }

    // Gaussian elimination
    for i in 0..n {
        // Pivot
        let mut pivot = i;
        while pivot < n && mat[pivot][i] == 0 {
            pivot += 1;
        }
        if pivot == n {
            return None;
        } // Singular

        if pivot != i {
            mat.swap(i, pivot);
        }

        // Normalize pivot row
        let pivot_val = mat[i][i].clone();
        for j in i..=n {
            mat[i][j] /= &pivot_val;
        }

        // Eliminate other rows
        for r in 0..n {
            if r != i {
                let factor = mat[r][i].clone();
                if factor != 0 {
                    for j in i..=n {
                        let sub = &factor * &mat[i][j];
                        mat[r][j] -= sub;
                    }
                }
            }
        }
    }

    // Extract solution
    let mut x = Vec::with_capacity(n);
    for row in mat.iter().take(n) {
        x.push(row[n].clone());
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use malachite::Integer;
    use malachite::Rational;

    fn i(n: i64) -> Integer {
        Integer::from(n)
    }

    fn r(n: i64) -> Rational {
        Rational::from(n)
    }

    #[test]
    fn test_solve_linear_system_rational_simple() {
        // 2x + y = 5
        // x + 3y = 5
        // Sol: x=2, y=1
        let m = vec![vec![r(2), r(1)], vec![r(1), r(3)]];
        let c = vec![r(5), r(5)];

        let x = solve_linear_system_rational(&m, &c).unwrap();
        assert_eq!(x[0], r(2));
        assert_eq!(x[1], r(1));
    }

    #[test]
    fn test_integer_kernel_simple() {
        // A = [1, 1, 1] (1x3 matrix)
        let a = vec![vec![i(1), i(1), i(1)]];
        let kernel = integer_kernel(&a);

        assert_eq!(kernel.len(), 2);
        for v in &kernel {
            let mut sum = i(0);
            for val in v {
                sum += val;
            }
            assert_eq!(sum, i(0));
        }
    }

    #[test]
    fn test_integer_kernel_independent() {
        // A = [1, 0, 0]
        //     [0, 1, 0]
        let a = vec![vec![i(1), i(0), i(0)], vec![i(0), i(1), i(0)]];
        let kernel = integer_kernel(&a);
        assert_eq!(kernel.len(), 1);
        assert_eq!(kernel[0], vec![i(0), i(0), i(1)]);
    }
}
