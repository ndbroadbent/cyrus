//! Curve basis computation from GLSM linear relations.
//!
//! Port of CYTools `set_divisor_basis` curve basis construction.

use crate::error::{Error, Result};
use crate::integer_math::{hermite_normal_form, sublattice_index_snf};
use malachite::Integer;

/// Compute the curve basis matrix given GLSM linear relations and a divisor basis.
///
/// - `linrels` is the GLSM linear relations matrix (including the origin column).
/// - `basis` is the list of divisor basis indices (columns).
pub fn compute_curve_basis_matrix(
    linrels: &[Vec<Integer>],
    basis: &[usize],
) -> Result<Vec<Vec<Integer>>> {
    if linrels.is_empty() {
        return Err(Error::InvalidInput("linrels matrix is empty".into()));
    }

    let n_cols = linrels[0].len();
    let mut nobasis: Vec<usize> = (0..n_cols).filter(|i| !basis.contains(i)).collect();
    nobasis.sort_unstable();

    // Reorder columns: nobasis first, then basis.
    let mut linrels_tmp: Vec<Vec<Integer>> = Vec::with_capacity(linrels.len());
    for row in linrels {
        let mut new_row = Vec::with_capacity(n_cols);
        for &idx in &nobasis {
            new_row.push(row[idx].clone());
        }
        for &idx in basis {
            new_row.push(row[idx].clone());
        }
        linrels_tmp.push(new_row);
    }

    // Compute HNF in the reordered basis.
    let linrels_hnf = hermite_normal_form(&linrels_tmp);

    // Reassemble columns back to original order.
    let mut linrels_new = vec![vec![Integer::from(0); n_cols]; linrels.len()];
    for (r, row) in linrels_hnf.iter().enumerate() {
        for (c, &idx) in nobasis.iter().enumerate() {
            linrels_new[r][idx] = row[c].clone();
        }
        for (c, &idx) in basis.iter().enumerate() {
            linrels_new[r][idx] = row[c + nobasis.len()].clone();
        }
    }

    let h11 = n_cols.saturating_sub(linrels.len());
    let mut curve_basis = vec![vec![Integer::from(0); n_cols]; h11];
    for (i, &b) in basis.iter().enumerate() {
        if i < h11 {
            curve_basis[i][b] = Integer::from(1);
        }
    }

    let sublat_ind = sublattice_index_snf(linrels);

    for &nb in nobasis.iter().rev() {
        let mut pivot: Option<(usize, Integer)> = None;
        for (k, row) in linrels_new.iter().enumerate() {
            let coeff = &row[nb];
            if *coeff != 0 {
                pivot = Some((k, coeff.clone()));
            }
        }
        let Some((i, ii)) = pivot else {
            continue;
        };
        if (&sublat_ind % &ii) != 0 {
            return Err(Error::InvalidInput("Problem with linear relations".into()));
        }
        for r in 0..h11 {
            let mut sum = Integer::from(0);
            for c in 0..n_cols {
                sum += &curve_basis[r][c] * &linrels_new[i][c];
            }
            curve_basis[r][nb] = -sum / &ii;
        }
    }

    Ok(curve_basis)
}
