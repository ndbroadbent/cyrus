//! Helper functions for intersection number computation.
//!
//! Utility functions including combinatorics, linear algebra solving,
//! and rational number conversion.

use crate::error::{Error, Result};
use faer::prelude::SpSolver;
use faer::sparse::SparseColMat;
use malachite::Rational;
use std::collections::HashMap;

/// Generate all k-combinations of a slice.
pub fn combinations(items: &[usize], k: usize) -> Vec<Vec<usize>> {
    if k == 0 {
        return vec![vec![]];
    }
    if items.len() < k {
        return vec![];
    }

    let mut result = Vec::new();
    for (i, &item) in items.iter().enumerate() {
        let rest = &items[i + 1..];
        for mut sub in combinations(rest, k - 1) {
            let mut combo = vec![item];
            combo.append(&mut sub);
            result.push(combo);
        }
    }
    result
}

/// Solve using normal equations: M^T M x = M^T C with sparse Cholesky.
pub fn solve_normal_equations(
    triplets: &[(usize, usize, f64)],
    c_vec: &[f64],
    n_rows: usize,
    n_cols: usize,
) -> Result<Vec<f64>> {
    // Group triplets by row for efficient M^T M computation
    let mut rows_data: HashMap<usize, Vec<(usize, f64)>> = HashMap::new();
    for &(row, col, val) in triplets {
        rows_data.entry(row).or_default().push((col, val));
    }

    // Build M^T M as sparse triplets
    let mut mtm_map: HashMap<(usize, usize), f64> = HashMap::new();
    for entries in rows_data.values() {
        for &(col_i, val_i) in entries {
            for &(col_j, val_j) in entries {
                *mtm_map.entry((col_i, col_j)).or_insert(0.0) += val_i * val_j;
            }
        }
    }

    // Add regularization to diagonal
    for i in 0..n_cols {
        *mtm_map.entry((i, i)).or_insert(0.0) += 1e-12;
    }

    // Convert to triplets for sparse matrix construction
    let mtm_triplets: Vec<(usize, usize, f64)> =
        mtm_map.into_iter().map(|((i, j), v)| (i, j, v)).collect();

    // Build sparse M^T M matrix
    let mtm_sparse =
        SparseColMat::<usize, f64>::try_new_from_triplets(n_cols, n_cols, &mtm_triplets)
            .map_err(|e| Error::SingularMatrix(format!("Failed to build sparse M^T M: {e:?}")))?;

    // Build M^T C as dense vector
    let mut mtc = faer::Mat::<f64>::zeros(n_cols, 1);
    for &(row, col, val) in triplets {
        mtc[(col, 0)] += val * c_vec[row];
    }

    // Solve using sparse Cholesky
    let chol = mtm_sparse
        .sp_cholesky(faer::Side::Lower)
        .map_err(|_| Error::SingularMatrix("Sparse Cholesky factorization failed".into()))?;

    let solution = chol.solve(&mtc);

    // Extract solution
    let sol: Vec<f64> = (0..n_cols).map(|i| solution[(i, 0)]).collect();

    if sol.iter().any(|&x| !x.is_finite()) {
        return Err(Error::SingularMatrix(
            "Solution contains non-finite values".into(),
        ));
    }

    // Compute residual: ||Mx - c||
    let mut residual_sq = 0.0;
    for row in 0..n_rows {
        let mut mx_row = 0.0;
        if let Some(entries) = rows_data.get(&row) {
            for &(col, val) in entries {
                mx_row += val * sol[col];
            }
        }
        let diff = mx_row - c_vec[row];
        residual_sq += diff * diff;
    }
    let residual = residual_sq.sqrt();
    let rhs_norm: f64 = c_vec.iter().map(|x| x * x).sum::<f64>().sqrt();
    let rel_residual = if rhs_norm > 1e-12 {
        residual / rhs_norm
    } else {
        residual
    };
    eprintln!(
        "[DEBUG] Residual: {residual:.6}, RHS norm: {rhs_norm:.6}, Relative: {rel_residual:.6}"
    );

    Ok(sol)
}

/// Convert f64 to a rational number.
pub fn float_to_rational(x: f64) -> Rational {
    // First try: is it close to an integer?
    let rounded = x.round();
    if (x - rounded).abs() < 1e-9 {
        #[allow(clippy::cast_possible_truncation)]
        return Rational::from(rounded as i64);
    }

    // Try common denominators
    for denom in [2, 3, 4, 5, 6, 8, 10, 12] {
        let scaled = x * f64::from(denom);
        let rounded_scaled = scaled.round();
        if (scaled - rounded_scaled).abs() < 1e-9 {
            #[allow(clippy::cast_possible_truncation)]
            return Rational::from(rounded_scaled as i64) / Rational::from(denom);
        }
    }

    // Continued fraction approximation
    continued_fraction_approx(x, 1000)
}

/// Approximate a float using continued fractions.
fn continued_fraction_approx(x: f64, max_denom: i64) -> Rational {
    let sign = x.signum();
    let x = x.abs();

    let mut a = x.floor();
    let mut h1 = 1_i64;
    #[allow(clippy::cast_possible_truncation)]
    let mut h0 = a as i64;
    let mut k1 = 0_i64;
    let mut k0 = 1_i64;

    let mut remainder = x - a;

    while remainder > 1e-12 && k0.abs() < max_denom {
        let inv = 1.0 / remainder;
        a = inv.floor();

        #[allow(clippy::cast_possible_truncation)]
        let a_i64 = a as i64;
        let h2 = a_i64 * h0 + h1;
        let k2 = a_i64 * k0 + k1;

        if k2.abs() > max_denom {
            break;
        }

        h1 = h0;
        h0 = h2;
        k1 = k0;
        k0 = k2;

        remainder = inv - a;
    }

    if sign < 0.0 {
        Rational::from(-h0) / Rational::from(k0)
    } else {
        Rational::from(h0) / Rational::from(k0)
    }
}
