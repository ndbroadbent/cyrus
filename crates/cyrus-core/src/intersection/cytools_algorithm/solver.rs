//! Sparse linear system construction and solving for the CYTools algorithm.

use crate::error::{Error, Result};
use faer::prelude::SpSolver;
use faer::sparse::SparseColMat;
use std::collections::HashMap;

/// Build the sparse linear system M*x + C = 0.
///
/// For each equation (3-tuple) and each linear_relations row:
/// - C component: sum of lin[point-1] * inv_det for known 4-forms
/// - M coefficients: lin[point-1] for variable contributions
pub fn build_linear_system(
    eqn_array: &[[usize; 3]],
    eqn_dict: &HashMap<[usize; 3], Vec<(usize, usize)>>,
    c_dict: &HashMap<[usize; 3], Vec<(usize, f64)>>,
    linear_relations: &[Vec<i64>],
    _n_vars: usize,
) -> (Vec<(usize, usize, f64)>, Vec<f64>) {
    let n_lin_rows = linear_relations.len();
    let num_rows = n_lin_rows * eqn_array.len();

    let mut m_triplets: Vec<(usize, usize, f64)> = Vec::new();
    let mut c_vec = vec![0.0; num_rows];

    let mut row_ctr = 0;
    for &eqn in eqn_array {
        for lin in linear_relations {
            // Compute C component (known 4-form contributions)
            // Only for distinct 3-tuples (eqn[0] != eqn[1] and eqn[1] != eqn[2])
            if eqn[0] != eqn[1]
                && eqn[1] != eqn[2]
                && let Some(c_entries) = c_dict.get(&eqn)
            {
                for &(point_idx, inv_det) in c_entries {
                    // lin[point_idx - 1] because linear_relations excludes origin
                    if point_idx > 0 && point_idx <= lin.len() {
                        let coeff = lin[point_idx - 1] as f64;
                        c_vec[row_ctr] += coeff * inv_det;
                    }
                }
            }

            // Compute M coefficients (variable contributions)
            if let Some(eqn_entries) = eqn_dict.get(&eqn) {
                for &(point_idx, var_idx) in eqn_entries {
                    if point_idx > 0 && point_idx <= lin.len() {
                        let coeff = lin[point_idx - 1] as f64;
                        if coeff.abs() > 1e-12 {
                            m_triplets.push((row_ctr, var_idx, coeff));
                        }
                    }
                }
            }

            row_ctr += 1;
        }
    }

    (m_triplets, c_vec)
}

/// Solve the sparse system M*x = -C using least squares.
pub fn solve_sparse_system(
    m_triplets: &[(usize, usize, f64)],
    c_vec: &[f64],
    n_vars: usize,
) -> Result<Vec<f64>> {
    if m_triplets.is_empty() || n_vars == 0 {
        return Ok(vec![0.0; n_vars]);
    }

    let n_rows = c_vec.len();

    let rows_data = group_triplets_by_row(m_triplets);
    let mtm_triplets = build_mtm_triplets(&rows_data, n_vars);
    let mtm_sparse = build_mtm_sparse(n_vars, &mtm_triplets)?;
    let mtc = build_mtc(m_triplets, c_vec, n_vars);

    let sol = solve_cholesky(&mtm_sparse, &mtc, n_vars)?;
    log_sparse_residual(&rows_data, &sol, c_vec, n_rows, n_vars);
    Ok(sol)
}

fn group_triplets_by_row(m_triplets: &[(usize, usize, f64)]) -> HashMap<usize, Vec<(usize, f64)>> {
    let mut rows_data: HashMap<usize, Vec<(usize, f64)>> = HashMap::new();
    for &(row, col, val) in m_triplets {
        rows_data.entry(row).or_default().push((col, val));
    }
    rows_data
}

fn build_mtm_triplets(
    rows_data: &HashMap<usize, Vec<(usize, f64)>>,
    _n_vars: usize,
) -> Vec<(usize, usize, f64)> {
    let mut mtm_map: HashMap<(usize, usize), f64> = HashMap::new();
    for entries in rows_data.values() {
        for &(col_i, val_i) in entries {
            for &(col_j, val_j) in entries {
                *mtm_map.entry((col_i, col_j)).or_insert(0.0) += val_i * val_j;
            }
        }
    }
    mtm_map.into_iter().map(|((i, j), v)| (i, j, v)).collect()
}

fn build_mtm_sparse(
    n_vars: usize,
    mtm_triplets: &[(usize, usize, f64)],
) -> Result<SparseColMat<usize, f64>> {
    SparseColMat::<usize, f64>::try_new_from_triplets(n_vars, n_vars, mtm_triplets)
        .map_err(|e| Error::SingularMatrix(format!("Failed to build M^T M: {e:?}")))
}

fn build_mtc(m_triplets: &[(usize, usize, f64)], c_vec: &[f64], n_vars: usize) -> faer::Mat<f64> {
    let mut mtc = faer::Mat::<f64>::zeros(n_vars, 1);
    for &(row, col, val) in m_triplets {
        mtc[(col, 0)] -= val * c_vec[row];
    }
    mtc
}

fn solve_cholesky(
    mtm_sparse: &SparseColMat<usize, f64>,
    mtc: &faer::Mat<f64>,
    n_vars: usize,
) -> Result<Vec<f64>> {
    let chol = mtm_sparse
        .sp_cholesky(faer::Side::Lower)
        .map_err(|_| Error::SingularMatrix("Cholesky failed".into()))?;
    let solution = chol.solve(&mtc);
    Ok((0..n_vars).map(|i| solution[(i, 0)]).collect())
}

fn log_sparse_residual(
    rows_data: &HashMap<usize, Vec<(usize, f64)>>,
    sol: &[f64],
    c_vec: &[f64],
    n_rows: usize,
    n_vars: usize,
) {
    let mut residual_sq = 0.0;
    let mut rhs_sq = 0.0;
    for row in 0..n_rows {
        let mut mx_row = 0.0;
        if let Some(entries) = rows_data.get(&row) {
            for &(col, val) in entries {
                mx_row += val * sol[col];
            }
        }
        let diff = mx_row + c_vec[row];
        residual_sq += diff * diff;
        rhs_sq += c_vec[row] * c_vec[row];
    }
    let rel_residual = if rhs_sq > 1e-12 {
        (residual_sq / rhs_sq).sqrt()
    } else {
        residual_sq.sqrt()
    };
    eprintln!("[CYTools algo] Relative residual: {rel_residual:.6} ({n_vars} vars, {n_rows} eqns)");
}
