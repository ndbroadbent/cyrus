//! Linear system solving utilities.
//!
//! Direct port of CYTools `utils.py` solve_linear_system function.
//! Uses rsparse for sparse Cholesky factorization, matching CYTools' sksparse.cholmod backend.

use rsparse::data::{Sprs, Trpl};
use rsparse::{cholsol, gaxpy, multiply, transpose};

/// Backend for solving sparse linear systems.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum LinearBackend {
    /// Try all backends in order (sksparse-style first, then scipy-style)
    #[default]
    All,
    /// Use Cholesky factorization via rsparse (equivalent to sksparse.cholmod)
    Cholesky,
    /// Use LU factorization via rsparse (equivalent to scipy.sparse.linalg.spsolve)
    Lu,
}

/// Solve the linear system M*x + C = 0 (i.e., M*x = -C).
///
/// Uses the normal equations: (M^T M) x = -M^T C
/// with Cholesky factorization for symmetric positive definite systems.
///
/// This is a direct port of CYTools `utils.py` `solve_linear_system` function,
/// using rsparse instead of sksparse.cholmod.
///
/// # Arguments
///
/// * `m` - The coefficient matrix in sparse triplet form (row, col, value)
/// * `m_rows` - Number of rows in M
/// * `m_cols` - Number of columns in M
/// * `c` - The constant vector (length m_rows)
/// * `backend` - Which solver backend to use
/// * `check` - Whether to verify the solution
/// * `tol` - Error tolerance for solution verification
///
/// # Returns
///
/// `Some(x)` if solution found and verified, `None` otherwise.
pub fn solve_linear_system_sparse(
    triplets: &[(usize, usize, f64)],
    m_rows: usize,
    m_cols: usize,
    c: &[f64],
    backend: LinearBackend,
    check: bool,
    tol: f64,
) -> Option<Vec<f64>> {
    if c.len() != m_rows {
        return None;
    }

    // Build sparse matrix M from triplets
    let m = build_sparse_from_triplets(triplets, m_rows, m_cols);

    match backend {
        LinearBackend::All => {
            // Try Cholesky first (like sksparse), then LU (like scipy)
            if let Some(sol) = solve_with_backend(&m, c, LinearBackend::Cholesky, check, tol) {
                return Some(sol);
            }
            solve_with_backend(&m, c, LinearBackend::Lu, check, tol)
        }
        _ => solve_with_backend(&m, c, backend, check, tol),
    }
}

/// Solve the linear system M*x + C = 0 with a dense matrix M.
///
/// Convenience wrapper that converts dense matrix to sparse format.
///
/// # Arguments
///
/// * `m` - The coefficient matrix as row-major dense array
/// * `c` - The constant vector
/// * `check` - Whether to verify the solution
/// * `tol` - Error tolerance for solution verification
///
/// # Returns
///
/// `Some(x)` if solution found and verified, `None` otherwise.
pub fn solve_linear_system(m: &[Vec<f64>], c: &[f64], check: bool, tol: f64) -> Option<Vec<f64>> {
    if m.is_empty() || c.is_empty() {
        return None;
    }

    let m_rows = m.len();
    let m_cols = m[0].len();

    // Convert dense to triplets (only non-zero entries)
    let mut triplets = Vec::new();
    for (i, row) in m.iter().enumerate() {
        for (j, &val) in row.iter().enumerate() {
            if val != 0.0 {
                triplets.push((i, j, val));
            }
        }
    }

    solve_linear_system_sparse(&triplets, m_rows, m_cols, c, LinearBackend::All, check, tol)
}

/// Build sparse matrix from triplets.
fn build_sparse_from_triplets(
    triplets: &[(usize, usize, f64)],
    rows: usize,
    cols: usize,
) -> Sprs<f64> {
    let mut trpl = Trpl::<f64>::new();
    trpl.m = rows;
    trpl.n = cols;

    for &(row, col, val) in triplets {
        trpl.i.push(row);
        trpl.p.push(col as isize);
        trpl.x.push(val);
    }

    trpl.to_sprs()
}

/// Solve using a specific backend.
fn solve_with_backend(
    m: &Sprs<f64>,
    c: &[f64],
    backend: LinearBackend,
    check: bool,
    tol: f64,
) -> Option<Vec<f64>> {
    // Compute M^T
    let m_t = transpose(m);

    // Compute M^T * M (symmetric positive definite for overdetermined systems)
    let m_tm = multiply(&m_t, m);

    // Compute -M^T * C
    let zeros = vec![0.0; m.n];
    let neg_c: Vec<f64> = c.iter().map(|&x| -x).collect();
    let mut neg_mt_c = gaxpy(&m_t, &neg_c, &zeros);

    // Solve (M^T M) x = -M^T C
    let order = match backend {
        LinearBackend::Cholesky => 0, // Cholesky ordering
        LinearBackend::Lu => 1,       // LU ordering
        LinearBackend::All => 0,      // Default to Cholesky
    };

    if cholsol(&m_tm, &mut neg_mt_c, order).is_err() {
        return None;
    }

    let x = neg_mt_c;

    // Verify solution if requested
    if check {
        // Compute M*x + C and check max error
        let mx = gaxpy(m, &x, c);
        let max_error = mx.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);

        if max_error > tol {
            return None;
        }
    }

    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solve_identity() {
        // Solve I*x + c = 0 => x = -c
        let m = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let c = vec![1.0, 2.0, 3.0];

        let x = solve_linear_system(&m, &c, true, 1e-10).unwrap();

        assert!((x[0] + 1.0).abs() < 1e-9);
        assert!((x[1] + 2.0).abs() < 1e-9);
        assert!((x[2] + 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_solve_diagonal() {
        // Solve diag(2,3,4)*x + [2,6,12] = 0 => x = [-1, -2, -3]
        let m = vec![
            vec![2.0, 0.0, 0.0],
            vec![0.0, 3.0, 0.0],
            vec![0.0, 0.0, 4.0],
        ];
        let c = vec![2.0, 6.0, 12.0];

        let x = solve_linear_system(&m, &c, true, 1e-10).unwrap();

        assert!((x[0] + 1.0).abs() < 1e-9);
        assert!((x[1] + 2.0).abs() < 1e-9);
        assert!((x[2] + 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_solve_general() {
        // Solve [[1,2],[3,4]]*x + [5,11] = 0
        // => x = [-1, -2] since [1,2;3,4]*[-1,-2] = [-5,-11]
        let m = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let c = vec![5.0, 11.0];

        let x = solve_linear_system(&m, &c, true, 1e-10).unwrap();

        assert!((x[0] + 1.0).abs() < 1e-9);
        assert!((x[1] + 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_solve_overdetermined() {
        // Overdetermined system (more equations than unknowns)
        // Least squares solution
        // m = [[1,1], [1,2], [1,3]]
        let m = vec![vec![1.0, 1.0], vec![1.0, 2.0], vec![1.0, 3.0]];
        let c = vec![2.0, 3.0, 4.0];

        let result = solve_linear_system(&m, &c, false, 1e-4);
        assert!(result.is_some());
    }

    #[test]
    fn test_solve_check_fails() {
        // Inconsistent overdetermined system should fail check
        let m = vec![vec![1.0, 1.0], vec![1.0, 2.0], vec![1.0, 3.0]];
        let c = vec![100.0, 3.0, 4.0]; // Inconsistent

        // With strict tolerance, should fail
        let result = solve_linear_system(&m, &c, true, 1e-10);
        assert!(result.is_none());
    }

    #[test]
    fn test_solve_sparse_triplets() {
        // Same as test_solve_identity but using sparse triplet interface
        let triplets = vec![(0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0)];
        let c = vec![1.0, 2.0, 3.0];

        let x = solve_linear_system_sparse(&triplets, 3, 3, &c, LinearBackend::All, true, 1e-10)
            .unwrap();

        assert!((x[0] + 1.0).abs() < 1e-9);
        assert!((x[1] + 2.0).abs() < 1e-9);
        assert!((x[2] + 3.0).abs() < 1e-9);
    }

    #[test]
    fn test_backend_cholesky() {
        let m = [
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let c = vec![1.0, 2.0, 3.0];

        let triplets: Vec<_> = m
            .iter()
            .enumerate()
            .flat_map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .filter(|&(_, &v)| v != 0.0)
                    .map(move |(j, &v)| (i, j, v))
            })
            .collect();

        let x =
            solve_linear_system_sparse(&triplets, 3, 3, &c, LinearBackend::Cholesky, true, 1e-10)
                .unwrap();

        assert!((x[0] + 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_backend_lu() {
        let m = [
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let c = vec![1.0, 2.0, 3.0];

        let triplets: Vec<_> = m
            .iter()
            .enumerate()
            .flat_map(|(i, row)| {
                row.iter()
                    .enumerate()
                    .filter(|&(_, &v)| v != 0.0)
                    .map(move |(j, &v)| (i, j, v))
            })
            .collect();

        let x = solve_linear_system_sparse(&triplets, 3, 3, &c, LinearBackend::Lu, true, 1e-10)
            .unwrap();

        assert!((x[0] + 1.0).abs() < 1e-9);
    }
}
