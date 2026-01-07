//! Integer linear algebra for exact orientation tests.
//!
//! Uses Integer arithmetic to avoid floating-point precision issues
//! in convex hull and triangulation computations.

use malachite::Integer;

/// Compute orientation (signed volume) of a simplex using Integer determinant.
/// Returns positive, negative, or zero.
pub fn orientation_integer(points: &[Vec<Integer>]) -> Integer {
    let n = points.len();
    if n == 0 {
        return Integer::from(0);
    }
    let dim = n - 1;
    if dim == 0 {
        return Integer::from(0);
    }

    // Build matrix of differences: M[i][j] = points[i+1][j] - points[0][j]
    let mut mat: Vec<Vec<Integer>> = vec![vec![Integer::from(0); dim]; dim];
    for i in 0..dim {
        for j in 0..dim {
            mat[i][j] = &points[i + 1][j] - &points[0][j];
        }
    }

    determinant_integer(&mut mat)
}

/// Compute determinant of a square Integer matrix using Bareiss algorithm.
/// This is fraction-free and avoids the expression swell of Gaussian elimination.
pub fn determinant_integer(mat: &mut [Vec<Integer>]) -> Integer {
    let n = mat.len();
    if n == 0 {
        return Integer::from(1);
    }
    if n == 1 {
        return mat[0][0].clone();
    }

    let mut sign = Integer::from(1);
    let mut prev_pivot = Integer::from(1);

    for k in 0..n {
        // Find pivot
        let mut pivot_row = None;
        for i in k..n {
            if mat[i][k] != 0 {
                pivot_row = Some(i);
                break;
            }
        }

        let Some(pivot_row) = pivot_row else {
            return Integer::from(0); // Singular matrix
        };

        if pivot_row != k {
            mat.swap(k, pivot_row);
            sign = -sign;
        }

        let pivot = mat[k][k].clone();

        // Bareiss algorithm: M[i][j] = (M[i][j] * pivot - M[i][k] * M[k][j]) / prev_pivot
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = &mat[i][j] * &pivot - &mat[i][k] * &mat[k][j];
                mat[i][j] = numerator / &prev_pivot;
            }
        }

        prev_pivot = pivot;
    }

    sign * &mat[n - 1][n - 1]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn int(n: i64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn test_determinant_integer_2x2() {
        // det([[1, 2], [3, 4]]) = 1*4 - 2*3 = -2
        let mut mat = vec![vec![int(1), int(2)], vec![int(3), int(4)]];
        assert_eq!(determinant_integer(&mut mat), int(-2));
    }

    #[test]
    fn test_determinant_integer_3x3() {
        // det([[1, 2, 3], [4, 5, 6], [7, 8, 10]]) = -3
        let mut mat = vec![
            vec![int(1), int(2), int(3)],
            vec![int(4), int(5), int(6)],
            vec![int(7), int(8), int(10)],
        ];
        assert_eq!(determinant_integer(&mut mat), int(-3));
    }

    #[test]
    fn test_determinant_integer_singular() {
        // Singular matrix
        let mut mat = vec![vec![int(1), int(2)], vec![int(2), int(4)]];
        assert_eq!(determinant_integer(&mut mat), int(0));
    }
}
