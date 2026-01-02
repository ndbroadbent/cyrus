//! Geometric Langlands Super-Multiplet (GLSM) charge matrix computation.
//!
//! The GLSM charge matrix encodes the linear relations among divisors in a toric variety.
//! For Calabi-Yau manifolds, these relations are constrained by the sum of charges
//! being zero (vanishing anti-canonical class).
//!
//! Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 2.

use crate::Point;
use crate::error::{Error, Result};
use crate::integer_math::integer_kernel;
use malachite::Integer;

/// Compute the GLSM charge matrix for a set of lattice points.
///
/// The charge matrix Q is the integer kernel of the point matrix
/// augmented with a row of 1s (enforcing the Calabi-Yau condition).
///
/// If `include_origin` is true, the origin point (0,...,0) is prepended to the
/// points list, matching the convention used in CYTools.
///
/// # Arguments
/// * `points` - Lattice points v_i defining the toric variety rays.
/// * `include_origin` - Whether to automatically prepend the origin point.
///
/// # Returns
/// A matrix (represented as a vector of rows) where each row is an integer
/// charge vector Q^a satisfying Σ_i Q_i^a v_i = 0 and Σ_i Q_i^a = 0.
///
/// # Errors
/// Returns an error if the input points are empty or inconsistent.
///
/// Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 2.2
pub fn compute_glsm_charge_matrix(
    points: &[Point],
    include_origin: bool,
) -> Result<Vec<Vec<Integer>>> {
    if points.is_empty() {
        return Err(Error::InvalidInput(
            "No points provided for GLSM computation".into(),
        ));
    }

    let dim = points[0].dim();

    let mut augmented_pts = Vec::new();
    if include_origin {
        augmented_pts.push(Point::origin(dim));
    }
    augmented_pts.extend_from_slice(points);

    let n_pts = augmented_pts.len();

    // Construct the matrix A ((dim + 1) x n_pts) using Integer for exactness
    let mut a_mat = vec![vec![Integer::from(0); n_pts]; dim + 1];
    for (j, pt) in augmented_pts.iter().enumerate() {
        let coords = pt.coords();
        for (i, &coord) in coords.iter().enumerate() {
            a_mat[i][j] = Integer::from(coord);
        }
        a_mat[dim][j] = Integer::from(1);
    }

    // 2. Compute the integer kernel
    let kernel = integer_kernel(&a_mat);

    Ok(kernel)
}

#[cfg(test)]
mod tests {
    use super::*;
    use malachite::Integer;

    fn i(n: i64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn test_glsm_p2() {
        // P^2 has vertices (1,0), (0,1), (-1,-1)
        let points = vec![
            Point::new(vec![1, 0]),
            Point::new(vec![0, 1]),
            Point::new(vec![-1, -1]),
        ];

        let q = compute_glsm_charge_matrix(&points, true).unwrap();

        // Expected: 1 relation for P^2
        assert_eq!(q.len(), 1);

        let row = &q[0];
        assert_eq!(row.len(), 4);

        // Verify relation Σ Qi v_i = 0
        let mut x = i(0);
        let mut y = i(0);
        let mut sum_q = i(0);

        // Origin
        sum_q += &row[0];
        // Vertices
        x += &row[1] * i(1);
        y += &row[1] * i(0);
        sum_q += &row[1];

        x += &row[2] * i(0);
        y += &row[2] * i(1);
        sum_q += &row[2];

        x += &row[3] * i(-1);
        y += &row[3] * i(-1);
        sum_q += &row[3];

        assert_eq!(x, i(0));
        assert_eq!(y, i(0));
        assert_eq!(sum_q, i(0));
    }

    #[test]
    fn test_glsm_p1xp1() {
        // P^1 x P^1 has 4 rays: (1,0), (-1,0), (0,1), (0,-1)
        let points = vec![
            Point::new(vec![1, 0]),
            Point::new(vec![-1, 0]),
            Point::new(vec![0, 1]),
            Point::new(vec![0, -1]),
        ];

        let q = compute_glsm_charge_matrix(&points, true).unwrap();

        // Expected: 2 relations for P1xP1
        assert_eq!(q.len(), 2);

        for row in &q {
            assert_eq!(row.len(), 5);
            let mut x = i(0);
            let mut y = i(0);
            let mut sum_q = i(0);

            // Origin
            sum_q += &row[0];
            // Vertices
            x += &row[1] * i(1);
            sum_q += &row[1];
            x += &row[2] * i(-1);
            sum_q += &row[2];
            y += &row[3] * i(1);
            sum_q += &row[3];
            y += &row[4] * i(-1);
            sum_q += &row[4];

            assert_eq!(x, i(0));
            assert_eq!(y, i(0));
            assert_eq!(sum_q, i(0));
        }
    }
}
