//! Orientation tests and exact determinant computation.
//!
//! Uses malachite::Integer for exact arithmetic at decision boundaries.

use malachite::Integer;
use std::cmp::Ordering;

/// Result of an orientation test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    /// Point is above the hyperplane (positive side).
    Positive,
    /// Point is below the hyperplane (negative side).
    Negative,
    /// Point is exactly on the hyperplane.
    Zero,
}

impl Orientation {
    pub(crate) fn from_sign(det: &Integer) -> Self {
        match det.cmp(&Integer::from(0)) {
            Ordering::Greater => Self::Positive,
            Ordering::Less => Self::Negative,
            Ordering::Equal => Self::Zero,
        }
    }
}

/// Compute the orientation of point p relative to hyperplane through v1..v_d.
///
/// Returns the sign of the determinant of the (d+1)×(d+1) matrix:
/// ```text
/// | 1  v1[0]  v1[1]  ...  v1[d-1] |
/// | 1  v2[0]  v2[1]  ...  v2[d-1] |
/// | ...                           |
/// | 1  vd[0]  vd[1]  ...  vd[d-1] |
/// | 1  p[0]   p[1]   ...  p[d-1]  |
/// ```
///
/// Uses exact integer arithmetic (malachite::Integer).
pub fn orientation(hyperplane_points: &[&[i64]], query_point: &[i64]) -> Orientation {
    let d = query_point.len();
    assert_eq!(
        hyperplane_points.len(),
        d,
        "Need exactly d points for hyperplane in d dimensions"
    );
    for p in hyperplane_points {
        assert_eq!(p.len(), d, "All points must have same dimension");
    }

    // Build the (d+1) × (d+1) matrix
    let n = d + 1;
    let mut matrix: Vec<Vec<Integer>> = Vec::with_capacity(n);

    // First d rows from hyperplane points
    for hp in hyperplane_points {
        let mut row = Vec::with_capacity(n);
        row.push(Integer::from(1));
        for &coord in *hp {
            row.push(Integer::from(coord));
        }
        matrix.push(row);
    }

    // Last row from query point
    let mut last_row = Vec::with_capacity(n);
    last_row.push(Integer::from(1));
    for &coord in query_point {
        last_row.push(Integer::from(coord));
    }
    matrix.push(last_row);

    let det = determinant_exact(&matrix);
    Orientation::from_sign(&det)
}

/// Compute exact determinant of an n×n matrix using Bareiss algorithm.
///
/// Bareiss algorithm is fraction-free Gaussian elimination that maintains
/// exact integer arithmetic throughout.
pub fn determinant_exact(matrix: &[Vec<Integer>]) -> Integer {
    let n = matrix.len();
    if n == 0 {
        return Integer::from(1);
    }

    // Clone matrix for in-place modification
    let mut m: Vec<Vec<Integer>> = matrix.to_vec();

    let mut sign = Integer::from(1);
    let mut prev_pivot = Integer::from(1);

    for k in 0..n {
        // Find pivot
        let mut pivot_row = None;
        for i in k..n {
            if m[i][k] != 0 {
                pivot_row = Some(i);
                break;
            }
        }

        let pivot_row = match pivot_row {
            Some(r) => r,
            None => return Integer::from(0), // Singular matrix
        };

        // Swap rows if needed
        if pivot_row != k {
            m.swap(k, pivot_row);
            sign = -sign;
        }

        let pivot = m[k][k].clone();

        // Bareiss elimination
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                // m[i][j] = (m[k][k] * m[i][j] - m[k][j] * m[i][k]) / prev_pivot
                let numerator = &pivot * &m[i][j] - &m[k][j] * &m[i][k];
                m[i][j] = numerator / &prev_pivot;
            }
            m[i][k] = Integer::from(0);
        }

        prev_pivot = pivot;
    }

    sign * &m[n - 1][n - 1]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orientation_2d() {
        // Triangle with vertices (0,0), (1,0), (0,1)
        // Point (0.5, 0.5) should be... let's compute
        // Actually for integer points, test with (1,1) vs the line from (0,0) to (2,0)
        let hp = [&[0i64, 0][..], &[2, 0][..]];

        // Point above line
        assert_eq!(orientation(&hp, &[1, 1]), Orientation::Positive);
        // Point below line
        assert_eq!(orientation(&hp, &[1, -1]), Orientation::Negative);
        // Point on line
        assert_eq!(orientation(&hp, &[1, 0]), Orientation::Zero);
    }

    #[test]
    fn test_determinant_exact() {
        // 2x2 identity
        let m = vec![
            vec![Integer::from(1), Integer::from(0)],
            vec![Integer::from(0), Integer::from(1)],
        ];
        assert_eq!(determinant_exact(&m), Integer::from(1));

        // 2x2 with det = ad - bc = 1*4 - 2*3 = -2
        let m = vec![
            vec![Integer::from(1), Integer::from(2)],
            vec![Integer::from(3), Integer::from(4)],
        ];
        assert_eq!(determinant_exact(&m), Integer::from(-2));
    }
}
