//! Linear system solving for rational matrices.
//!
//! Provides Gaussian elimination solver for systems of linear equations
//! over the rationals.

use malachite::Rational;

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
