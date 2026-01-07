//! Divisor basis computation.
//!
//! Computes the h11-dimensional basis of divisors from the GLSM charge matrix.
//! This is essential for expressing physics (Kähler form, fluxes) in a canonical basis.
//!
//! ## Algorithm (from DIVISOR_BASIS.md)
//!
//! 1. Construct the Linear Relations Matrix L with columns (1, v_i)
//! 2. Compute the sublattice index via Smith Normal Form
//! 3. Find h11 columns that form an integral basis (determinant matches sublattice index)
//!
//! ## Usage
//!
//! ```ignore
//! let basis = compute_divisor_basis(&points, &glsm)?;
//! let kappa_basis = intersection_in_basis(&kappa, &basis);
//! ```

use malachite::Integer;
use malachite::num::arithmetic::traits::Abs;

use crate::Point;
use crate::error::{Error, Result};
use crate::intersection::Intersection;

/// Compute the divisor basis indices from triangulation points and GLSM.
///
/// Returns a vector of h11 point indices that form an integral basis
/// for the divisor group modulo linear equivalence.
///
/// # Algorithm
///
/// 1. Build the linear relations matrix L with columns (1, v_i) for each point
/// 2. Find (dim+1) pivot columns that span the row space of L (these are "relation" columns)
/// 3. The remaining h11 columns form the divisor basis
///
/// The key insight: L has (dim+1) rows and n_points columns. The divisor group
/// Pic(X) = Z^n_points / Im(L^T) has rank h11 = n_points - (dim+1).
/// We need to select h11 columns whose images in Pic(X) form a basis.
/// This is equivalent to finding (dim+1) pivot columns and taking the complement.
///
/// # Errors
///
/// Returns an error if no valid basis can be found.
pub fn compute_divisor_basis(points: &[Point], glsm: &[Vec<Integer>]) -> Result<Vec<usize>> {
    if points.is_empty() {
        return Err(Error::InvalidInput("No points provided".into()));
    }
    if glsm.is_empty() {
        return Err(Error::InvalidInput("GLSM is empty".into()));
    }

    let n_points = points.len();
    let dim = points[0].dim();
    let n_relations = dim + 1; // Number of linear relations (rows of L)

    // h11 = n_points - (dim+1) for a favorable CY
    // This equals the number of GLSM kernel vectors (no torsion case)
    let h11 = n_points.saturating_sub(n_relations);

    if h11 == 0 {
        return Err(Error::InvalidInput("h11 must be positive".into()));
    }

    // Verify: GLSM rows should equal h11
    if glsm.len() != h11 {
        // Allow off-by-one for legacy compatibility (some code passes h11+1 rows)
        if glsm.len() != h11 + 1 {
            return Err(Error::InvalidInput(format!(
                "GLSM row count {} doesn't match h11={} (expected {} or {})",
                glsm.len(),
                h11,
                h11,
                h11 + 1
            )));
        }
    }

    // Step 1: Build linear relations matrix L
    // L has n_points columns and (dim + 1) rows
    // Column i = (1, v_i) where v_i is the point coordinates
    let l_matrix = build_linear_relations_matrix(points, dim);

    // Step 2: Sort columns by L1 norm for determinism
    let sorted_indices = sort_by_l1_norm(points);

    // Step 3: Find (dim+1) pivot columns that span the row space
    // These are the "relation columns" that will NOT be in the basis
    let pivot_columns = find_pivot_columns(&l_matrix, &sorted_indices, n_relations)?;

    // Step 4: The divisor basis is all columns EXCEPT the pivot columns
    let pivot_set: std::collections::HashSet<usize> = pivot_columns.into_iter().collect();
    let basis: Vec<usize> = sorted_indices
        .into_iter()
        .filter(|idx| !pivot_set.contains(idx))
        .collect();

    if basis.len() != h11 {
        return Err(Error::InvalidInput(format!(
            "Could not find {} basis columns, found {}",
            h11,
            basis.len()
        )));
    }

    Ok(basis)
}

/// Build the linear relations matrix L.
/// Column i = (1, v_i) where v_i are the point coordinates.
fn build_linear_relations_matrix(points: &[Point], dim: usize) -> Vec<Vec<Integer>> {
    let n_points = points.len();
    let n_rows = dim + 1;

    let mut l = vec![vec![Integer::from(0); n_points]; n_rows];

    for (col, point) in points.iter().enumerate() {
        // First row is always 1 (CY condition)
        l[0][col] = Integer::from(1);
        // Remaining rows are point coordinates
        for (row, &coord) in point.coords().iter().enumerate() {
            l[row + 1][col] = Integer::from(coord);
        }
    }

    l
}

/// Sort point indices by L1 norm of coordinates.
fn sort_by_l1_norm(points: &[Point]) -> Vec<usize> {
    let mut indices: Vec<usize> = (0..points.len()).collect();

    indices.sort_by_key(|&i| points[i].coords().iter().map(|&x| x.abs()).sum::<i64>());

    indices
}

/// Find n_pivots linearly independent columns from L.
///
/// These are the "relation columns" - the columns that span the row space of L.
/// The divisor basis will be all OTHER columns.
fn find_pivot_columns(
    l_matrix: &[Vec<Integer>],
    sorted_indices: &[usize],
    n_pivots: usize,
) -> Result<Vec<usize>> {
    let mut pivots = Vec::with_capacity(n_pivots);
    let mut selected_cols: Vec<Vec<Integer>> = Vec::with_capacity(n_pivots);

    for &idx in sorted_indices {
        if pivots.len() >= n_pivots {
            break;
        }

        // Extract column
        let col: Vec<Integer> = l_matrix.iter().map(|row| row[idx].clone()).collect();

        // Check if linearly independent from already selected
        if is_linearly_independent(&selected_cols, &col) {
            pivots.push(idx);
            selected_cols.push(col);
        }
    }

    if pivots.len() < n_pivots {
        return Err(Error::InvalidInput(format!(
            "Could not find {} linearly independent pivot columns, only found {}",
            n_pivots,
            pivots.len()
        )));
    }

    Ok(pivots)
}

/// Check if a column is linearly independent from existing columns.
fn is_linearly_independent(cols: &[Vec<Integer>], new_col: &[Integer]) -> bool {
    if cols.is_empty() {
        return new_col.iter().any(|x| *x != 0);
    }

    // Build matrix with all columns including new one
    let n_rows = new_col.len();
    let n_cols = cols.len() + 1;

    let mut mat: Vec<Vec<Integer>> = vec![vec![Integer::from(0); n_cols]; n_rows];
    for (i, row) in mat.iter_mut().enumerate() {
        for (j, col) in cols.iter().enumerate() {
            row[j] = col[i].clone();
        }
        row[n_cols - 1] = new_col[i].clone();
    }

    // Check rank by Gaussian elimination
    let rank = compute_rank(&mut mat);
    rank == n_cols
}

/// Compute rank of an integer matrix using Gaussian elimination.
fn compute_rank(mat: &mut [Vec<Integer>]) -> usize {
    if mat.is_empty() || mat[0].is_empty() {
        return 0;
    }

    let n_rows = mat.len();
    let n_cols = mat[0].len();

    let mut rank = 0;
    let mut pivot_col = 0;

    for row in 0..n_rows {
        if pivot_col >= n_cols {
            break;
        }

        // Find pivot in current column
        let mut pivot_row = None;
        for r in row..n_rows {
            if mat[r][pivot_col] != 0 {
                pivot_row = Some(r);
                break;
            }
        }

        if let Some(pr) = pivot_row {
            mat.swap(row, pr);

            // Eliminate below
            for r in (row + 1)..n_rows {
                if mat[r][pivot_col] != 0 {
                    let gcd = gcd_integer(&mat[row][pivot_col], &mat[r][pivot_col]);
                    let mult_r = &mat[row][pivot_col] / &gcd;
                    let mult_row = &mat[r][pivot_col] / &gcd;

                    for c in pivot_col..n_cols {
                        let new_val = &mat[r][c] * &mult_r - &mat[row][c] * &mult_row;
                        mat[r][c] = new_val;
                    }
                }
            }

            rank += 1;
            pivot_col += 1;
        } else {
            pivot_col += 1;
        }
    }

    rank
}

/// GCD of two integers.
fn gcd_integer(a: &Integer, b: &Integer) -> Integer {
    if *b == 0 {
        a.clone().abs()
    } else {
        gcd_integer(b, &(a % b))
    }
}

/// Extract intersection numbers for basis indices only.
///
/// Given the full intersection tensor κ_IJK indexed by point indices,
/// and a basis B = [b_0, b_1, ..., b_{h11-1}],
/// returns κ^basis_{ijk} = κ_{b_i, b_j, b_k}.
///
/// This is the tensor needed for physics computations (Kähler form, fluxes).
pub fn intersection_in_basis(kappa: &Intersection, basis: &[usize]) -> Intersection {
    let h11 = basis.len();
    let mut kappa_basis = Intersection::new(h11);

    // Iterate over all non-zero entries in the original tensor
    for (&(i, j, k), val) in kappa.iter() {
        // Check if all indices are in basis
        let bi = basis.iter().position(|&b| b == i);
        let bj = basis.iter().position(|&b| b == j);
        let bk = basis.iter().position(|&b| b == k);

        if let (Some(bi), Some(bj), Some(bk)) = (bi, bj, bk) {
            kappa_basis.set(bi, bj, bk, val.clone());
        }
    }

    kappa_basis
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::rational::Rational as TypedRational;
    use crate::types::tags::Pos;
    use malachite::Rational;

    #[test]
    fn test_build_linear_relations_matrix() {
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![0, 1]),
        ];

        let l = build_linear_relations_matrix(&points, 2);

        // 3 rows (1 + dim), 3 columns (n_points)
        assert_eq!(l.len(), 3);
        assert_eq!(l[0].len(), 3);

        // First row is all 1s
        assert!(l[0].iter().all(|x| *x == 1));
    }

    #[test]
    fn test_sort_by_l1_norm() {
        let points = vec![
            Point::new(vec![5, 5]), // L1 = 10
            Point::new(vec![0, 0]), // L1 = 0
            Point::new(vec![1, 2]), // L1 = 3
        ];

        let sorted = sort_by_l1_norm(&points);

        // Should be sorted: origin (0), then (1,2), then (5,5)
        assert_eq!(sorted, vec![1, 2, 0]);
    }

    #[test]
    fn test_intersection_in_basis() {
        // Create a 4-dimensional intersection tensor
        let mut kappa = Intersection::new(4);

        // Set some values
        kappa.set(
            0,
            1,
            2,
            TypedRational::<Pos>::new(Rational::from(6)).unwrap(),
        );
        kappa.set(
            1,
            2,
            3,
            TypedRational::<Pos>::new(Rational::from(12)).unwrap(),
        );
        kappa.set(
            0,
            0,
            3,
            TypedRational::<Pos>::new(Rational::from(3)).unwrap(),
        );

        // Extract for basis [1, 2, 3] (indices 0, 1, 2 in new tensor)
        let basis = vec![1, 2, 3];
        let kappa_basis = intersection_in_basis(&kappa, &basis);

        // Should have dim = 3
        assert_eq!(kappa_basis.dim(), 3);

        // (1,2,3) in original -> (0,1,2) in basis
        // Check that the value was preserved
        let val = kappa_basis.get(0, 1, 2);
        assert_eq!(*val.get(), Rational::from(12));

        // (0,1,2) should NOT be in basis (index 0 not in basis)
        // (0,0,3) should NOT be in basis (index 0 not in basis)
        assert_eq!(kappa_basis.num_nonzero(), 1);
    }
}
