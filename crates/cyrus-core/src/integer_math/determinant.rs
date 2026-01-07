//! Determinant and orientation computations.
//!
//! Provides algorithms for computing determinants and orientations of point sets.

use malachite::Rational;

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
