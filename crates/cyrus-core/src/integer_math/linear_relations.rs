//! GLSM charge matrix and linear relations computation.
//!
//! These algorithms follow CYTools' approach for computing the GLSM charge
//! matrix and linear relations from polytope points.

use malachite::Integer;
use malachite::Rational;

use super::matrix_utils::{invert_matrix, lcm_integer, matrix_rank, rational_matrix_to_integer};

/// Compute the GLSM charge matrix and linear relations from polytope points.
///
/// This is the correct algorithm used by CYTools:
/// 1. Build augmented point matrix L = [1s; points.T]
/// 2. Compute row echelon form (HNF-style) with specific pivot selection
/// 3. Extract both GLSM and linear_relations from the result
///
/// # Arguments
/// * `points` - The lattice points (n_points x dim), INCLUDING the origin at index 0
///
/// # Returns
/// A tuple (glsm, linear_relations, basis_indices) where:
/// - glsm: h11 x n_points charge matrix with identity on basis columns
/// - linear_relations: h11 x n_points matrix with -I on non-basis columns
/// - basis_indices: indices of the basis columns
///
/// # Note
/// The linear_relations matrix (with origin excluded) is what's used for
/// intersection number computation.
pub fn compute_glsm_and_linear_relations(
    points: &[Vec<i64>],
) -> (Vec<Vec<Integer>>, Vec<Vec<Integer>>, Vec<usize>) {
    if points.is_empty() {
        return (Vec::new(), Vec::new(), Vec::new());
    }

    let n_pts = points.len();
    let dim = points[0].len();

    // Build augmented point matrix: L = [1s; points.T]
    // L has shape (dim+1) x n_pts
    let mut linrel: Vec<Vec<Rational>> = vec![vec![Rational::from(0); n_pts]; dim + 1];

    // First row: all 1s (homogeneous coordinate)
    for j in 0..n_pts {
        linrel[0][j] = Rational::from(1);
    }

    // Remaining rows: point coordinates transposed
    for (j, pt) in points.iter().enumerate() {
        for (i, &coord) in pt.iter().enumerate() {
            linrel[i + 1][j] = Rational::from(coord);
        }
    }

    // Compute row echelon form, tracking pivot columns
    let n_rows = dim + 1;
    let mut pivot_cols: Vec<usize> = Vec::new();
    let mut current_row = 0;

    for col in 0..n_pts {
        if current_row >= n_rows {
            break;
        }

        // Find pivot in this column
        let mut pivot_row = None;
        for r in current_row..n_rows {
            if linrel[r][col] != 0 {
                pivot_row = Some(r);
                break;
            }
        }

        if let Some(pr) = pivot_row {
            // Swap rows
            linrel.swap(current_row, pr);

            // Normalize pivot row
            let pivot_val = linrel[current_row][col].clone();
            for c in 0..n_pts {
                linrel[current_row][c] /= &pivot_val;
            }

            // Eliminate in all other rows
            for r in 0..n_rows {
                if r != current_row && linrel[r][col] != 0 {
                    let factor = linrel[r][col].clone();
                    for c in 0..n_pts {
                        let sub = &factor * &linrel[current_row][c];
                        linrel[r][c] -= sub;
                    }
                }
            }

            pivot_cols.push(col);
            current_row += 1;
        }
    }

    // The pivot columns are the "excluded" (non-basis) columns
    // The non-pivot columns are the basis columns
    let basis_ind: Vec<usize> = (0..n_pts).filter(|c| !pivot_cols.contains(c)).collect();
    let h11 = n_pts - (dim + 1); // h11 = n_points - (dim + 1)

    // Build GLSM: h11 x n_pts matrix
    // Identity on basis columns, derive non-basis columns from linrel
    let mut glsm: Vec<Vec<Rational>> = vec![vec![Rational::from(0); n_pts]; h11];

    // Set identity on basis columns
    for (i, &b) in basis_ind.iter().enumerate() {
        if i < h11 {
            glsm[i][b] = Rational::from(1);
        }
    }

    // Fill non-basis columns using the linear relations from REF
    // For each pivot column (non-basis), express it in terms of basis columns
    for &nb in pivot_cols.iter().rev() {
        // Find which row has the pivot in this column
        let pivot_row =
            (0..n_rows).find(|&r| linrel[r][nb] == 1 && (0..nb).all(|c| linrel[r][c] == 0));

        if let Some(pr) = pivot_row {
            // glsm[:, nb] = -glsm @ linrel[pr]
            for i in 0..h11 {
                let mut sum = Rational::from(0);
                for c in 0..n_pts {
                    sum += &glsm[i][c] * &linrel[pr][c];
                }
                glsm[i][nb] = -sum;
            }
        }
    }

    // Build linear_relations: negate and reorder
    // linear_relations has -1 on pivot columns (non-basis)
    let mut lr: Vec<Vec<Rational>> = linrel.clone();

    // Negate all entries
    for row in &mut lr {
        for val in row {
            *val = -val.clone();
        }
    }

    // Convert to integers
    let glsm_int = rational_matrix_to_integer(&glsm);
    let lr_int = rational_matrix_to_integer(&lr);

    (glsm_int, lr_int, basis_ind)
}

/// Compute linear relations from points, excluding the origin.
///
/// This is the main entry point for intersection number computation.
/// Returns the linear_relations matrix with origin column removed.
///
/// The algorithm follows CYTools' approach:
/// 1. Find dim points that form a lattice basis (excluding origin)
/// 2. The remaining points form the GLSM basis
/// 3. Express lattice basis points in terms of GLSM basis points
/// 4. Format as linear_relations with -1 on lattice basis columns
///
/// # Arguments
/// * `points` - The lattice points (n_points x dim), with origin at index 0
///
/// # Returns
/// The linear_relations matrix with shape (h11 x (n_points - 1)),
/// origin column excluded. Has -1 on first h11 columns (lattice basis),
/// coefficients on remaining columns (GLSM basis).
pub fn compute_linear_relations_no_origin(points: &[Vec<i64>]) -> Vec<Vec<Integer>> {
    if points.is_empty() || points.len() < 2 {
        return Vec::new();
    }

    let _n_pts = points.len();
    let dim = points[0].len();

    // Points without origin (indices 1..n_pts in original, 0..n_pts-1 here)
    let pts_no_origin: Vec<Vec<Rational>> = points[1..]
        .iter()
        .map(|p| p.iter().map(|&x| Rational::from(x)).collect())
        .collect();

    let n_pts_no_origin = pts_no_origin.len(); // n_pts - 1
    let h11 = n_pts_no_origin.saturating_sub(dim);

    if h11 == 0 || n_pts_no_origin < dim {
        return Vec::new();
    }

    // Find dim points that form a lattice basis (points that span R^dim)
    // Use greedy selection based on rank
    let mut lattice_basis_indices: Vec<usize> = Vec::new();

    for i in 0..n_pts_no_origin {
        let mut test_pts: Vec<Vec<Rational>> = lattice_basis_indices
            .iter()
            .map(|&j| pts_no_origin[j].clone())
            .collect();
        test_pts.push(pts_no_origin[i].clone());

        let rank = matrix_rank(&test_pts);
        if rank > lattice_basis_indices.len() {
            lattice_basis_indices.push(i);
            if lattice_basis_indices.len() == dim {
                break;
            }
        }
    }

    if lattice_basis_indices.len() != dim {
        // Couldn't find full rank lattice basis
        return Vec::new();
    }

    // GLSM basis = remaining points (not in lattice basis)
    let glsm_basis_indices: Vec<usize> = (0..n_pts_no_origin)
        .filter(|i| !lattice_basis_indices.contains(i))
        .collect();

    // Build matrix M from lattice basis points (dim x dim)
    let mut m_mat: Vec<Vec<Rational>> = vec![vec![Rational::from(0); dim]; dim];
    for (i, &lb_idx) in lattice_basis_indices.iter().enumerate() {
        for (j, val) in pts_no_origin[lb_idx].iter().enumerate() {
            m_mat[j][i] = val.clone(); // Transpose: columns are points
        }
    }

    // Compute M^(-1)
    let m_inv = match invert_matrix(&m_mat) {
        Some(inv) => inv,
        None => return Vec::new(),
    };

    // Build matrix of GLSM basis points (dim x h11)
    let mut glsm_pts: Vec<Vec<Rational>> =
        vec![vec![Rational::from(0); glsm_basis_indices.len()]; dim];
    for (j, &gb_idx) in glsm_basis_indices.iter().enumerate() {
        for (i, val) in pts_no_origin[gb_idx].iter().enumerate() {
            glsm_pts[i][j] = val.clone();
        }
    }

    // extra_pts = -M^(-1) @ glsm_pts
    // This expresses GLSM basis points as linear combos of lattice basis points
    // Shape: dim x h11
    let mut _extra_pts: Vec<Vec<Rational>> =
        vec![vec![Rational::from(0); glsm_basis_indices.len()]; dim];
    for i in 0..dim {
        for j in 0..glsm_basis_indices.len() {
            let mut sum = Rational::from(0);
            for k in 0..dim {
                sum += &m_inv[i][k] * &glsm_pts[k][j];
            }
            _extra_pts[i][j] = -sum;
        }
    }

    // Build linear_relations using RREF approach which is cleaner
    build_linear_relations_from_rref(points)
}

/// Build linear_relations using RREF of augmented point matrix.
fn build_linear_relations_from_rref(points: &[Vec<i64>]) -> Vec<Vec<Integer>> {
    let n_pts = points.len();
    let dim = points[0].len();
    let h11 = n_pts - (dim + 1);

    if h11 == 0 {
        return Vec::new();
    }

    // Build augmented point matrix: L = [1s; points.T]
    let mut mat: Vec<Vec<Rational>> = vec![vec![Rational::from(0); n_pts]; dim + 1];
    for j in 0..n_pts {
        mat[0][j] = Rational::from(1);
    }
    for (j, pt) in points.iter().enumerate() {
        for (i, &coord) in pt.iter().enumerate() {
            mat[i + 1][j] = Rational::from(coord);
        }
    }

    // Compute RREF
    let n_rows = dim + 1;
    let mut pivot_cols: Vec<usize> = Vec::new();
    let mut current_row = 0;

    for col in 0..n_pts {
        if current_row >= n_rows {
            break;
        }

        let mut pivot_row = None;
        for r in current_row..n_rows {
            if mat[r][col] != 0 {
                pivot_row = Some(r);
                break;
            }
        }

        if let Some(pr) = pivot_row {
            mat.swap(current_row, pr);
            let pivot_val = mat[current_row][col].clone();
            for c in 0..n_pts {
                mat[current_row][c] /= &pivot_val;
            }
            for r in 0..n_rows {
                if r != current_row && mat[r][col] != 0 {
                    let factor = mat[r][col].clone();
                    for c in 0..n_pts {
                        let sub = &factor * &mat[current_row][c];
                        mat[r][c] -= sub;
                    }
                }
            }
            pivot_cols.push(col);
            current_row += 1;
        }
    }

    // Free columns (GLSM basis) are those not in pivot_cols
    let _free_cols: Vec<usize> = (0..n_pts).filter(|c| !pivot_cols.contains(c)).collect();

    // Build linear_relations from RREF
    // For each non-origin pivot column, create a row in linear_relations
    // The row has -1 at the pivot column, and the negated RREF entries at free columns

    let mut lr: Vec<Vec<Integer>> = Vec::with_capacity(h11);

    for (row_idx, &pivot_col) in pivot_cols.iter().enumerate() {
        if pivot_col == 0 {
            continue; // Skip origin column
        }
        if lr.len() >= h11 {
            break;
        }

        let mut row: Vec<Integer> = vec![Integer::from(0); n_pts - 1]; // Exclude origin

        // Find LCM of denominators in this RREF row
        let mut lcm = Integer::from(1);
        for c in 0..n_pts {
            if c != 0 {
                // Skip origin
                let (_, denom) = mat[row_idx][c].clone().into_numerator_and_denominator();
                lcm = lcm_integer(&lcm, &Integer::from(denom));
            }
        }

        // Fill the row with scaled integer values
        // RREF row: [0, ..., 1, ..., a, b, c, d] where 1 is at pivot position
        // This encodes: D_pivot + a*D_4 + b*D_5 + c*D_6 + d*D_7 = 0
        //
        // CYTools linear_relations uses the RREF row directly (no negation):
        // [1, 0, 0, 0, a, b, c, d]
        //
        // The relation is: sum(coeff_i * D_i) = 0

        for c in 1..n_pts {
            // Skip origin (column 0)
            let val = &mat[row_idx][c] * Rational::from(&lcm);

            // After scaling by LCM, this should be an integer. Use try_from which
            // preserves sign (unlike into_numerator_and_denominator which returns Natural).
            let int_val =
                Integer::try_from(val).expect("Value should be integer after LCM scaling");

            // Adjusted column index (c-1 because we skip origin)
            let col_idx = c - 1;

            row[col_idx] = int_val;
        }

        lr.push(row);
    }

    lr
}
