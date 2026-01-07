//! Affine geometry utilities.
//!
//! Direct port of CYTools `utils.py` affine geometry functions.

use malachite::Integer;

use crate::integer_math::integer_kernel;

/// Find new points that are affinely independent to the input list of points.
///
/// This is useful when one wants to turn a polytope that is not
/// full-dimensional into one that is, without affecting the structure of the
/// triangulations.
///
/// Direct port of CYTools `utils.py` `find_new_affinely_independent_points`.
///
/// # Arguments
///
/// * `pts` - A list of points (as rows)
///
/// # Returns
///
/// A list of affinely independent points with respect to the ones inputted.
///
/// # Panics
///
/// Panics if the input is empty.
///
/// # Example
///
/// ```rust
/// use cyrus_core::utils::find_new_affinely_independent_points;
///
/// let pts = vec![
///     vec![1, 0, 1],
///     vec![0, 0, 1],
///     vec![0, 1, 1],
/// ];
/// let new_pts = find_new_affinely_independent_points(&pts);
/// // Returns points that together with input span the full space
/// ```
pub fn find_new_affinely_independent_points(pts: &[Vec<i64>]) -> Vec<Vec<i64>> {
    assert!(!pts.is_empty(), "List of points cannot be empty.");

    let n_pts = pts.len();
    let dim = pts[0].len();

    // Translate so first point is at origin
    let translation = pts[0].clone();
    let translated: Vec<Vec<i64>> = pts
        .iter()
        .map(|p| p.iter().zip(&translation).map(|(a, b)| a - b).collect())
        .collect();

    // Handle single point case
    let mut working_pts = translated;
    if n_pts == 1 {
        let mut unit = vec![0i64; dim];
        if dim > 0 {
            unit[0] = 1;
        }
        working_pts.push(unit);
    }

    // Compute rank using floating point (for efficiency)
    let rank = compute_rank(&working_pts);

    // Build basis of linearly independent points
    let mut basis: Vec<Vec<i64>> = Vec::new();
    let mut basis_rank = 0;

    for pt in &working_pts {
        basis.push(pt.clone());
        let new_rank = compute_rank(&basis);
        if new_rank > basis_rank {
            basis_rank = new_rank;
        } else {
            basis.pop();
        }
        if basis_rank == rank {
            break;
        }
    }

    // Find kernel (nullspace) of basis matrix
    // Convert to Integer matrix for exact computation
    let int_basis: Vec<Vec<Integer>> = basis
        .iter()
        .map(|row| row.iter().map(|&x| Integer::from(x)).collect())
        .collect();

    let kernel = integer_kernel(&int_basis);

    // Convert kernel to i64 and transpose (kernel rows -> new points)
    let mut new_pts: Vec<Vec<i64>> = kernel
        .iter()
        .map(|row| row.iter().map(|x| i64::try_from(x).unwrap_or(0)).collect())
        .collect();

    // Handle single point case: add unit vector
    if n_pts == 1 && dim > 0 {
        let mut unit = vec![0i64; dim];
        unit[0] = 1;
        new_pts.push(unit);
    }

    // Translate back
    new_pts
        .iter()
        .map(|p| p.iter().zip(&translation).map(|(a, b)| a + b).collect())
        .collect()
}

/// Compute the rank of a matrix using floating point arithmetic.
fn compute_rank(matrix: &[Vec<i64>]) -> usize {
    if matrix.is_empty() {
        return 0;
    }

    let rows = matrix.len();
    let cols = matrix[0].len();

    // Convert to f64 for SVD-based rank computation
    let mut mat: Vec<Vec<f64>> = matrix
        .iter()
        .map(|row| row.iter().map(|&x| x as f64).collect())
        .collect();

    // Gaussian elimination to count pivots
    let mut rank = 0;
    let mut pivot_col = 0;

    for row in 0..rows {
        if pivot_col >= cols {
            break;
        }

        // Find pivot
        let mut max_row = row;
        let mut max_val = mat[row][pivot_col].abs();
        for r in (row + 1)..rows {
            if mat[r][pivot_col].abs() > max_val {
                max_val = mat[r][pivot_col].abs();
                max_row = r;
            }
        }

        if max_val < 1e-10 {
            pivot_col += 1;
            continue;
        }

        mat.swap(row, max_row);
        rank += 1;

        // Eliminate
        for r in (row + 1)..rows {
            if mat[r][pivot_col].abs() > 1e-15 {
                let factor = mat[r][pivot_col] / mat[row][pivot_col];
                for c in pivot_col..cols {
                    mat[r][c] -= factor * mat[row][c];
                }
            }
        }

        pivot_col += 1;
    }

    rank
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_affinely_independent_simple() {
        // Points in a plane in 3D
        let pts = vec![vec![1, 0, 1], vec![0, 0, 1], vec![0, 1, 1]];
        let new_pts = find_new_affinely_independent_points(&pts);

        // Should find at least one new point
        // The new point(s) should make the span full-dimensional
        assert!(!new_pts.is_empty());
    }

    #[test]
    fn test_compute_rank() {
        // Full rank 2x2
        let m1 = vec![vec![1, 0], vec![0, 1]];
        assert_eq!(compute_rank(&m1), 2);

        // Rank 1 (second row is multiple of first)
        let m2 = vec![vec![1, 2], vec![2, 4]];
        assert_eq!(compute_rank(&m2), 1);

        // Rank 2 for 3x2 matrix
        let m3 = vec![vec![1, 0], vec![0, 1], vec![1, 1]];
        assert_eq!(compute_rank(&m3), 2);
    }

    #[test]
    fn test_single_point() {
        let pts = vec![vec![1, 2, 3]];
        let new_pts = find_new_affinely_independent_points(&pts);
        // Should return points to span the space
        assert!(!new_pts.is_empty());
    }

    #[test]
    #[should_panic(expected = "List of points cannot be empty")]
    fn test_empty_points_panics() {
        let pts: Vec<Vec<i64>> = Vec::new();
        find_new_affinely_independent_points(&pts);
    }
}
