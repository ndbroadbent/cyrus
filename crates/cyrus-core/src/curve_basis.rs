//! Curve basis computation from GLSM linear relations.
//!
//! Port of CYTools `set_divisor_basis` curve basis construction.

use crate::error::{Error, Result};
use crate::integer_math::{hermite_normal_form, invert_matrix, sublattice_index_snf};
use crate::types::f64::F64;
use crate::types::tags::Finite;
use malachite::Integer;
use malachite::Rational;

/// Divisor-basis representation accepted by CYTools-style curve-basis setup.
///
/// `Indices` is the ordinary CYTools path where the basis is a vector of ambient
/// divisor column indices. `Matrix` is the generic divisor-basis path where each
/// row is an ambient divisor linear combination and `standard_basis` is the
/// ordinary GLSM basis used to reduce and invert the matrix block.
#[derive(Clone, Copy, Debug)]
pub enum DivisorBasis<'a> {
    /// Ambient divisor column indices, matching CYTools' vector basis path.
    Indices(&'a [usize]),
    /// Ambient divisor linear-combination rows, matching CYTools' matrix path.
    Matrix {
        /// Standard GLSM column basis used to reduce the matrix block.
        standard_basis: &'a [usize],
        /// Divisor-basis rows in ambient divisor coordinates, including origin.
        basis_matrix: &'a [Vec<Integer>],
    },
}

/// Compute the GLSM-coordinate columns of a divisor basis.
///
/// The returned matrix has one column per basis divisor. For an index basis,
/// this is `glsm[:, basis]`. For a matrix basis whose rows are ambient divisor
/// combinations, column `a` is `glsm * basis_matrix[a]^T`.
///
/// # Errors
/// Returns an error if the GLSM matrix or basis representation is malformed.
pub fn divisor_basis_glsm_coordinate_matrix(
    glsm: &[Vec<Integer>],
    basis: DivisorBasis<'_>,
) -> Result<Vec<Vec<Integer>>> {
    let (h11, n_cols) = validate_glsm_matrix(glsm)?;
    match basis {
        DivisorBasis::Indices(indices) => {
            validate_standard_basis(indices, h11, n_cols)?;
            let mut coords = vec![vec![Integer::from(0); h11]; h11];
            for (basis_col, &ambient_col) in indices.iter().enumerate() {
                for row in 0..h11 {
                    coords[row][basis_col] = glsm[row][ambient_col].clone();
                }
            }
            Ok(coords)
        }
        DivisorBasis::Matrix {
            standard_basis,
            basis_matrix,
        } => {
            validate_standard_basis(standard_basis, h11, n_cols)?;
            validate_matrix_basis_shape(basis_matrix, h11, n_cols)?;
            let mut coords = vec![vec![Integer::from(0); h11]; h11];
            for (basis_col, basis_row) in basis_matrix.iter().enumerate() {
                for glsm_row in 0..h11 {
                    let mut acc = Integer::from(0);
                    for (ambient_col, coefficient) in basis_row.iter().enumerate() {
                        if *coefficient != 0 {
                            acc += &glsm[glsm_row][ambient_col] * coefficient;
                        }
                    }
                    coords[glsm_row][basis_col] = acc;
                }
            }
            Ok(coords)
        }
    }
}

/// Compute the integer change-of-basis matrix between divisor-basis shapes.
///
/// Returns `T` such that `to_coords = from_coords * T`, where coordinates are
/// GLSM-coordinate columns from [`divisor_basis_glsm_coordinate_matrix`].
///
/// # Errors
/// Returns an error if either basis is malformed, singular, or the transform is
/// not integral.
pub fn divisor_basis_change_matrix(
    glsm: &[Vec<Integer>],
    from_basis: DivisorBasis<'_>,
    to_basis: DivisorBasis<'_>,
) -> Result<Vec<Vec<Integer>>> {
    let from_coords = divisor_basis_glsm_coordinate_matrix(glsm, from_basis)?;
    let to_coords = divisor_basis_glsm_coordinate_matrix(glsm, to_basis)?;
    basis_change_matrix_from_coordinate_matrices(&from_coords, &to_coords)
}

/// Apply an integer basis transform to integer coordinates.
///
/// Given a transform `T` and coordinate vector `v`, returns `T v`. This is the
/// convention used for M-flux and Kähler-coordinate style transforms in the
/// McAllister pipeline.
///
/// # Errors
/// Returns an error if the transform is not square, the vector length is
/// inconsistent, or an output coordinate does not fit in `i64`.
pub fn apply_integer_basis_transform(
    transform: &[Vec<Integer>],
    values: &[i64],
    context: &str,
) -> Result<Vec<i64>> {
    validate_basis_transform_shape(transform, values.len(), context)?;
    transform
        .iter()
        .map(|row| {
            let mut acc = Integer::from(0);
            for (coeff, &value) in row.iter().zip(values.iter()) {
                acc += coeff * Integer::from(value);
            }
            i64::try_from(&acc).map_err(|_| {
                Error::InvalidInput(format!(
                    "{context} transformed coordinate does not fit in i64"
                ))
            })
        })
        .collect()
}

/// Apply the transpose of an integer basis transform to integer coordinates.
///
/// Given a transform `T` and coordinate vector `v`, returns `T^T v`. This is
/// the convention used for K-flux/source-coordinate transforms.
///
/// # Errors
/// Returns an error if the transform is not square, the vector length is
/// inconsistent, or an output coordinate does not fit in `i64`.
pub fn apply_integer_basis_transform_transpose(
    transform: &[Vec<Integer>],
    values: &[i64],
    context: &str,
) -> Result<Vec<i64>> {
    let len = validate_basis_transform_shape(transform, values.len(), context)?;
    (0..len)
        .map(|col| {
            let mut acc = Integer::from(0);
            for (row, &value) in transform.iter().zip(values.iter()) {
                acc += &row[col] * Integer::from(value);
            }
            i64::try_from(&acc).map_err(|_| {
                Error::InvalidInput(format!(
                    "{context} transformed coordinate does not fit in i64"
                ))
            })
        })
        .collect()
}

/// Apply an integer basis transform to finite floating-point coordinates.
///
/// Given a transform `T` and coordinate vector `v`, returns `T v`, preserving
/// the `F64<Finite>` boundary for Kähler-style coordinates.
///
/// # Errors
/// Returns an error if the transform shape is inconsistent, a transform
/// coefficient cannot be represented as finite `f64`, or an output coordinate
/// is not finite.
pub fn apply_finite_f64_basis_transform(
    transform: &[Vec<Integer>],
    values: &[F64<Finite>],
    context: &str,
) -> Result<Vec<F64<Finite>>> {
    validate_basis_transform_shape(transform, values.len(), context)?;
    transform
        .iter()
        .map(|row| {
            let mut out = F64::<Finite>::ZERO;
            for (coeff, value) in row.iter().zip(values.iter()) {
                let raw_coeff = coeff.to_string().parse::<f64>().map_err(|_| {
                    Error::InvalidInput(format!(
                        "{context} basis transform coefficient does not fit in f64"
                    ))
                })?;
                let coeff_f = F64::<Finite>::new(raw_coeff).ok_or_else(|| {
                    Error::InvalidInput(format!(
                        "{context} basis transform coefficient is not finite"
                    ))
                })?;
                out = out + coeff_f * *value;
            }
            F64::<Finite>::new(out.get()).ok_or_else(|| {
                Error::InvalidInput(format!("{context} transformed coordinate is not finite"))
            })
        })
        .collect()
}

/// Compute the curve basis matrix given GLSM linear relations and a divisor basis.
///
/// - `linrels` is the GLSM linear relations matrix (including the origin column).
/// - `basis` is the list of divisor basis indices (columns).
pub fn compute_curve_basis_matrix(
    linrels: &[Vec<Integer>],
    basis: &[usize],
) -> Result<Vec<Vec<Integer>>> {
    let n_cols = validate_linrels(linrels)?;
    let h11 = n_cols.saturating_sub(linrels.len());
    validate_standard_basis(basis, h11, n_cols)?;
    let nobasis = compute_nobasis(n_cols, basis);
    let linrels_tmp = reorder_linrels(linrels, &nobasis, basis, n_cols);
    let linrels_hnf = hermite_normal_form(&linrels_tmp);
    let linrels_new = restore_linrels(&linrels_hnf, &nobasis, basis, n_cols, linrels.len());
    let mut curve_basis = initial_curve_basis(n_cols, h11, basis);
    let sublat_ind = sublattice_index_snf(linrels);
    apply_nobasis_constraints(
        &mut curve_basis,
        &linrels_new,
        &nobasis,
        n_cols,
        h11,
        &sublat_ind,
    )?;
    Ok(curve_basis)
}

/// Compute the dual curve-basis matrix for either supported divisor-basis shape.
///
/// This is the typed dispatch point for callers that may receive either a
/// CYTools vector basis or a generic matrix basis. It avoids accidentally
/// routing a matrix basis through index-column selection.
///
/// # Errors
/// Returns an error if the selected basis representation is malformed or does
/// not define an integral divisor basis.
pub fn compute_curve_basis_matrix_for_divisor_basis(
    linrels: &[Vec<Integer>],
    basis: DivisorBasis<'_>,
) -> Result<Vec<Vec<Integer>>> {
    match basis {
        DivisorBasis::Indices(indices) => compute_curve_basis_matrix(linrels, indices),
        DivisorBasis::Matrix {
            standard_basis,
            basis_matrix,
        } => compute_curve_basis_matrix_from_divisor_basis_matrix(
            linrels,
            standard_basis,
            basis_matrix,
        ),
    }
}

/// Compute the dual curve-basis matrix for a matrix divisor basis.
///
/// This ports the matrix-basis branch of CYTools `set_divisor_basis`. The
/// `basis_matrix` rows are divisor-basis vectors in ambient divisor coordinates,
/// including the origin column. `standard_basis` is the ordinary CYTools GLSM
/// basis used to reduce those rows against `linrels`.
///
/// # Errors
/// Returns an error if the matrix dimensions are inconsistent, the divisor
/// matrix does not define an integral basis in the standard GLSM basis, or the
/// non-basis columns cannot be reconstructed from the linear relations.
pub fn compute_curve_basis_matrix_from_divisor_basis_matrix(
    linrels: &[Vec<Integer>],
    standard_basis: &[usize],
    basis_matrix: &[Vec<Integer>],
) -> Result<Vec<Vec<Integer>>> {
    let n_cols = validate_linrels(linrels)?;
    let h11 = n_cols.saturating_sub(linrels.len());
    validate_matrix_basis_shape(basis_matrix, h11, n_cols)?;
    validate_standard_basis(standard_basis, h11, n_cols)?;

    let mut divisor_basis_mat = basis_matrix.to_vec();
    let nobasis = compute_nobasis(n_cols, standard_basis);
    let sublat_ind = sublattice_index_snf(linrels);
    reduce_divisor_basis_matrix(
        &mut divisor_basis_mat,
        linrels,
        &nobasis,
        n_cols,
        &sublat_ind,
    )?;

    let standard_block = extract_standard_block(&divisor_basis_mat, standard_basis);
    let standard_block_inv = invert_integer_matrix(&standard_block).ok_or_else(|| {
        Error::InvalidInput("matrix divisor basis does not form an integral basis".into())
    })?;

    let mut curve_basis = vec![vec![Integer::from(0); n_cols]; h11];
    for row in 0..h11 {
        for (col, &standard_col) in standard_basis.iter().enumerate() {
            curve_basis[row][standard_col] = standard_block_inv[col][row].clone();
        }
    }
    apply_nobasis_constraints(
        &mut curve_basis,
        linrels,
        &nobasis,
        n_cols,
        h11,
        &sublat_ind,
    )?;
    Ok(curve_basis)
}

/// Compute the no-origin `q` matrix passed to `cygv` for either divisor-basis
/// shape.
///
/// # Errors
/// Returns an error if curve-basis construction fails or if the no-origin matrix
/// cannot be represented as `i64`.
pub fn curve_basis_q_matrix_for_divisor_basis_i64(
    linrels: &[Vec<Integer>],
    basis: DivisorBasis<'_>,
) -> Result<Vec<Vec<i64>>> {
    let curve_basis = compute_curve_basis_matrix_for_divisor_basis(linrels, basis)?;
    curve_basis_matrix_without_origin_i64(&curve_basis)
}

/// Return a CYTools-style no-origin `q` matrix as `i64` rows.
///
/// CYTools passes `curve_basis(include_origin=False, as_matrix=True)` to
/// `cygv`, which means the origin/canonical divisor column is omitted while the
/// rows remain Kähler-basis curve rows. This helper centralizes that boundary
/// so callers do not hand-roll `skip(1)` and accidentally keep the origin
/// column in GV inputs.
///
/// # Errors
/// Returns an error if the curve-basis matrix is empty, rows have inconsistent
/// widths, the origin column is absent, or an entry does not fit in `i64`.
pub fn curve_basis_matrix_without_origin_i64(
    curve_basis: &[Vec<Integer>],
) -> Result<Vec<Vec<i64>>> {
    if curve_basis.is_empty() {
        return Err(Error::InvalidInput("curve basis matrix is empty".into()));
    }
    let n_cols = curve_basis[0].len();
    if n_cols < 2 {
        return Err(Error::InvalidInput(
            "curve basis matrix must include an origin column".into(),
        ));
    }
    for row in curve_basis {
        if row.len() != n_cols {
            return Err(Error::InvalidInput(
                "curve basis matrix rows have inconsistent length".into(),
            ));
        }
    }

    curve_basis
        .iter()
        .map(|row| {
            row.iter()
                .skip(1)
                .map(|value| {
                    i64::try_from(value).map_err(|_| {
                        Error::InvalidInput("curve basis entry does not fit in i64".into())
                    })
                })
                .collect()
        })
        .collect()
}

fn validate_linrels(linrels: &[Vec<Integer>]) -> Result<usize> {
    if linrels.is_empty() {
        return Err(Error::InvalidInput("linrels matrix is empty".into()));
    }
    let n_cols = linrels[0].len();
    if n_cols == 0 {
        return Err(Error::InvalidInput("linrels matrix has no columns".into()));
    }
    for row in linrels {
        if row.len() != n_cols {
            return Err(Error::InvalidInput(
                "linrels matrix rows have inconsistent length".into(),
            ));
        }
    }
    Ok(n_cols)
}

fn validate_glsm_matrix(glsm: &[Vec<Integer>]) -> Result<(usize, usize)> {
    if glsm.is_empty() {
        return Err(Error::InvalidInput("GLSM matrix is empty".into()));
    }
    let n_cols = glsm[0].len();
    if n_cols == 0 {
        return Err(Error::InvalidInput("GLSM matrix has no columns".into()));
    }
    for row in glsm {
        if row.len() != n_cols {
            return Err(Error::InvalidInput(
                "GLSM matrix rows have inconsistent length".into(),
            ));
        }
    }
    Ok((glsm.len(), n_cols))
}

fn validate_basis_transform_shape(
    transform: &[Vec<Integer>],
    value_len: usize,
    context: &str,
) -> Result<usize> {
    if transform.is_empty() {
        return Err(Error::InvalidInput(format!(
            "{context} basis transform is empty"
        )));
    }
    let len = transform.len();
    if value_len != len {
        return Err(Error::InvalidInput(format!(
            "{context} basis transform shape does not match vector length"
        )));
    }
    for row in transform {
        if row.len() != len {
            return Err(Error::InvalidInput(format!(
                "{context} basis transform is not square"
            )));
        }
    }
    Ok(len)
}

fn validate_matrix_basis_shape(
    basis_matrix: &[Vec<Integer>],
    h11: usize,
    n_cols: usize,
) -> Result<()> {
    if basis_matrix.len() != h11 {
        return Err(Error::InvalidInput(format!(
            "matrix divisor basis row count {} does not match h11={h11}",
            basis_matrix.len()
        )));
    }
    for row in basis_matrix {
        if row.len() != n_cols {
            return Err(Error::InvalidInput(format!(
                "matrix divisor basis row width {} does not match ambient divisor count {n_cols}",
                row.len()
            )));
        }
    }
    Ok(())
}

fn validate_standard_basis(standard_basis: &[usize], h11: usize, n_cols: usize) -> Result<()> {
    if standard_basis.len() != h11 {
        return Err(Error::InvalidInput(format!(
            "standard basis length {} does not match h11={h11}",
            standard_basis.len()
        )));
    }
    let mut seen = std::collections::BTreeSet::new();
    for &idx in standard_basis {
        if idx >= n_cols {
            return Err(Error::InvalidInput(format!(
                "standard basis index {idx} is out of bounds for {n_cols} columns"
            )));
        }
        if !seen.insert(idx) {
            return Err(Error::InvalidInput(
                "standard basis contains duplicate indices".into(),
            ));
        }
    }
    Ok(())
}

fn reduce_divisor_basis_matrix(
    divisor_basis_mat: &mut [Vec<Integer>],
    linrels: &[Vec<Integer>],
    nobasis: &[usize],
    n_cols: usize,
    sublat_ind: &Integer,
) -> Result<()> {
    for &nb in nobasis.iter().rev() {
        let Some((row_idx, coeff)) = find_pivot(linrels, nb) else {
            continue;
        };
        if (sublat_ind % &coeff) != 0 {
            return Err(Error::InvalidInput("Problem with linear relations".into()));
        }
        for row in divisor_basis_mat.iter_mut() {
            let scale = row[nb].clone();
            if scale == 0 {
                continue;
            }
            for col in 0..n_cols {
                row[col] -= &scale * &linrels[row_idx][col];
            }
        }
    }
    Ok(())
}

fn extract_standard_block(
    divisor_basis_mat: &[Vec<Integer>],
    standard_basis: &[usize],
) -> Vec<Vec<Integer>> {
    divisor_basis_mat
        .iter()
        .map(|row| standard_basis.iter().map(|&idx| row[idx].clone()).collect())
        .collect()
}

fn invert_integer_matrix(matrix: &[Vec<Integer>]) -> Option<Vec<Vec<Integer>>> {
    let rational_matrix = matrix
        .iter()
        .map(|row| row.iter().map(Rational::from).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let inverse = invert_matrix(&rational_matrix)?;
    rational_matrix_to_integer(&inverse)
}

fn rational_matrix_to_integer(matrix: &[Vec<Rational>]) -> Option<Vec<Vec<Integer>>> {
    let mut out = Vec::with_capacity(matrix.len());
    for row in matrix {
        let mut out_row = Vec::with_capacity(row.len());
        for value in row {
            let integer = Integer::try_from(value.clone()).ok()?;
            out_row.push(integer);
        }
        out.push(out_row);
    }
    Some(out)
}

fn multiply_square(left: &[Vec<Rational>], right: &[Vec<Rational>]) -> Vec<Vec<Rational>> {
    let n = left.len();
    let mut out = vec![vec![Rational::from(0); n]; n];
    for i in 0..n {
        for j in 0..n {
            let mut acc = Rational::from(0);
            for k in 0..n {
                acc += &left[i][k] * &right[k][j];
            }
            out[i][j] = acc;
        }
    }
    out
}

fn to_rational_matrix(matrix: &[Vec<Integer>]) -> Vec<Vec<Rational>> {
    matrix
        .iter()
        .map(|row| row.iter().map(Rational::from).collect())
        .collect()
}

fn basis_change_matrix_from_coordinate_matrices(
    from_coords: &[Vec<Integer>],
    to_coords: &[Vec<Integer>],
) -> Result<Vec<Vec<Integer>>> {
    if from_coords.is_empty() || to_coords.is_empty() {
        return Err(Error::InvalidInput(
            "basis coordinate matrix is empty".into(),
        ));
    }
    let h11 = from_coords.len();
    if from_coords.iter().any(|row| row.len() != h11)
        || to_coords.len() != h11
        || to_coords.iter().any(|row| row.len() != h11)
    {
        return Err(Error::InvalidInput(
            "basis coordinate matrices must be square with matching dimensions".into(),
        ));
    }

    let from_r = to_rational_matrix(from_coords);
    let to_r = to_rational_matrix(to_coords);
    let from_inv = invert_matrix(&from_r).ok_or_else(|| {
        Error::InvalidInput("failed to invert source basis coordinate matrix".into())
    })?;
    let transform_r = multiply_square(&from_inv, &to_r);
    rational_matrix_to_integer(&transform_r)
        .ok_or_else(|| Error::InvalidInput("basis change matrix is not integral".into()))
}

fn compute_nobasis(n_cols: usize, basis: &[usize]) -> Vec<usize> {
    let mut nobasis: Vec<usize> = (0..n_cols).filter(|i| !basis.contains(i)).collect();
    nobasis.sort_unstable();
    nobasis
}

fn reorder_linrels(
    linrels: &[Vec<Integer>],
    nobasis: &[usize],
    basis: &[usize],
    n_cols: usize,
) -> Vec<Vec<Integer>> {
    let mut linrels_tmp: Vec<Vec<Integer>> = Vec::with_capacity(linrels.len());
    for row in linrels {
        let mut new_row = Vec::with_capacity(n_cols);
        for &idx in nobasis {
            new_row.push(row[idx].clone());
        }
        for &idx in basis {
            new_row.push(row[idx].clone());
        }
        linrels_tmp.push(new_row);
    }
    linrels_tmp
}

fn restore_linrels(
    linrels_hnf: &[Vec<Integer>],
    nobasis: &[usize],
    basis: &[usize],
    n_cols: usize,
    n_rows: usize,
) -> Vec<Vec<Integer>> {
    let mut linrels_new = vec![vec![Integer::from(0); n_cols]; n_rows];
    for (r, row) in linrels_hnf.iter().enumerate() {
        for (c, &idx) in nobasis.iter().enumerate() {
            linrels_new[r][idx] = row[c].clone();
        }
        for (c, &idx) in basis.iter().enumerate() {
            linrels_new[r][idx] = row[c + nobasis.len()].clone();
        }
    }
    linrels_new
}

fn initial_curve_basis(n_cols: usize, h11: usize, basis: &[usize]) -> Vec<Vec<Integer>> {
    let mut curve_basis = vec![vec![Integer::from(0); n_cols]; h11];
    for (i, &b) in basis.iter().enumerate() {
        if i < h11 {
            curve_basis[i][b] = Integer::from(1);
        }
    }
    curve_basis
}

fn apply_nobasis_constraints(
    curve_basis: &mut [Vec<Integer>],
    linrels_new: &[Vec<Integer>],
    nobasis: &[usize],
    n_cols: usize,
    h11: usize,
    sublat_ind: &Integer,
) -> Result<()> {
    for &nb in nobasis.iter().rev() {
        let pivot = find_pivot(linrels_new, nb);
        let Some((row_idx, coeff)) = pivot else {
            continue;
        };
        if (sublat_ind % &coeff) != 0 {
            return Err(Error::InvalidInput("Problem with linear relations".into()));
        }
        for r in 0..h11 {
            let mut sum = Integer::from(0);
            for c in 0..n_cols {
                sum += &curve_basis[r][c] * &linrels_new[row_idx][c];
            }
            curve_basis[r][nb] = -sum / &coeff;
        }
    }
    Ok(())
}

fn find_pivot(linrels_new: &[Vec<Integer>], col: usize) -> Option<(usize, Integer)> {
    let mut pivot: Option<(usize, Integer)> = None;
    for (k, row) in linrels_new.iter().enumerate() {
        let coeff = &row[col];
        if *coeff != 0 {
            pivot = Some((k, coeff.clone()));
        }
    }
    pivot
}

#[cfg(test)]
mod tests {
    use super::*;

    fn int_matrix(rows: &[&[i64]]) -> Vec<Vec<Integer>> {
        rows.iter()
            .map(|row| row.iter().map(|&entry| Integer::from(entry)).collect())
            .collect()
    }

    #[test]
    fn divisor_basis_dispatch_matches_vector_constructor() {
        let linrels = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let basis = vec![2, 3];

        let direct = compute_curve_basis_matrix(&linrels, &basis).unwrap();
        let dispatched =
            compute_curve_basis_matrix_for_divisor_basis(&linrels, DivisorBasis::Indices(&basis))
                .unwrap();

        assert_eq!(dispatched, direct);
    }

    #[test]
    fn divisor_basis_dispatch_matches_matrix_constructor() {
        let linrels = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let standard_basis = vec![2, 3];
        let divisor_basis_matrix = int_matrix(&[&[0, 0, 1, 1], &[0, 0, 0, 1]]);

        let direct = compute_curve_basis_matrix_from_divisor_basis_matrix(
            &linrels,
            &standard_basis,
            &divisor_basis_matrix,
        )
        .unwrap();
        let dispatched = compute_curve_basis_matrix_for_divisor_basis(
            &linrels,
            DivisorBasis::Matrix {
                standard_basis: &standard_basis,
                basis_matrix: &divisor_basis_matrix,
            },
        )
        .unwrap();

        assert_eq!(dispatched, direct);
    }

    #[test]
    fn divisor_basis_dispatch_builds_no_origin_q_matrix() {
        let linrels = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let standard_basis = vec![2, 3];
        let divisor_basis_matrix = int_matrix(&[&[0, 0, 1, 1], &[0, 0, 0, 1]]);

        let q_matrix = curve_basis_q_matrix_for_divisor_basis_i64(
            &linrels,
            DivisorBasis::Matrix {
                standard_basis: &standard_basis,
                basis_matrix: &divisor_basis_matrix,
            },
        )
        .unwrap();

        assert_eq!(q_matrix, vec![vec![2, 1, 0], vec![1, -1, 1]]);
    }

    #[test]
    fn divisor_basis_coordinate_matrix_dispatches_indices_and_matrix_basis() {
        let glsm = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let standard_basis = vec![2, 3];
        let divisor_basis_matrix = int_matrix(&[&[0, 0, 1, 1], &[0, 0, 0, 1]]);

        let vector_coords =
            divisor_basis_glsm_coordinate_matrix(&glsm, DivisorBasis::Indices(&standard_basis))
                .unwrap();
        let matrix_coords = divisor_basis_glsm_coordinate_matrix(
            &glsm,
            DivisorBasis::Matrix {
                standard_basis: &standard_basis,
                basis_matrix: &divisor_basis_matrix,
            },
        )
        .unwrap();

        assert_eq!(vector_coords, int_matrix(&[&[-1, -1], &[-2, -3]]));
        assert_eq!(matrix_coords, int_matrix(&[&[-2, -1], &[-5, -3]]));
    }

    #[test]
    fn divisor_basis_change_matrix_supports_matrix_basis() {
        let glsm = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let standard_basis = vec![2, 3];
        let divisor_basis_matrix = int_matrix(&[&[0, 0, 1, 1], &[0, 0, 0, 1]]);

        let computed_to_matrix = divisor_basis_change_matrix(
            &glsm,
            DivisorBasis::Indices(&standard_basis),
            DivisorBasis::Matrix {
                standard_basis: &standard_basis,
                basis_matrix: &divisor_basis_matrix,
            },
        )
        .unwrap();
        let matrix_to_computed = divisor_basis_change_matrix(
            &glsm,
            DivisorBasis::Matrix {
                standard_basis: &standard_basis,
                basis_matrix: &divisor_basis_matrix,
            },
            DivisorBasis::Indices(&standard_basis),
        )
        .unwrap();

        assert_eq!(computed_to_matrix, int_matrix(&[&[1, 0], &[1, 1]]));
        assert_eq!(matrix_to_computed, int_matrix(&[&[1, 0], &[-1, 1]]));
    }

    #[test]
    fn basis_transform_helpers_apply_integer_and_finite_coordinates() {
        let transform = int_matrix(&[&[1, 0], &[1, 1]]);
        let values = vec![5, 7];
        let finite_values = vec![
            F64::<Finite>::new(5.0).unwrap(),
            F64::<Finite>::new(7.0).unwrap(),
        ];

        let integer = apply_integer_basis_transform(&transform, &values, "unit").unwrap();
        let transpose =
            apply_integer_basis_transform_transpose(&transform, &values, "unit").unwrap();
        let finite = apply_finite_f64_basis_transform(&transform, &finite_values, "unit").unwrap();

        assert_eq!(integer, vec![5, 12]);
        assert_eq!(transpose, vec![12, 7]);
        assert_eq!(
            finite.iter().map(|value| value.get()).collect::<Vec<_>>(),
            vec![5.0, 12.0]
        );
    }

    #[test]
    fn basis_transform_helpers_reject_bad_shapes() {
        let transform = int_matrix(&[&[1, 0, 0], &[0, 1, 0]]);
        let err = apply_integer_basis_transform(&transform, &[1, 2], "unit")
            .expect_err("non-square transform should fail");

        assert!(err.to_string().contains("not square"));
    }

    #[test]
    fn vector_divisor_basis_rejects_wrong_length() {
        let linrels = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);

        let err =
            compute_curve_basis_matrix(&linrels, &[2]).expect_err("basis length must equal h11");

        assert!(err.to_string().contains("standard basis length"));
    }

    #[test]
    fn vector_divisor_basis_rejects_duplicate_indices() {
        let linrels = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);

        let err = compute_curve_basis_matrix(&linrels, &[2, 2])
            .expect_err("basis indices must be unique");

        assert!(err.to_string().contains("duplicate indices"));
    }

    #[test]
    fn vector_divisor_basis_rejects_out_of_bounds_index() {
        let linrels = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);

        let err = compute_curve_basis_matrix(&linrels, &[2, 4])
            .expect_err("basis index must be in range");

        assert!(err.to_string().contains("out of bounds"));
    }

    #[test]
    fn matrix_divisor_basis_identity_matches_vector_basis_curve_basis() {
        let linrels = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let standard_basis = vec![2, 3];
        let divisor_basis_matrix = int_matrix(&[&[0, 0, 1, 0], &[0, 0, 0, 1]]);

        let vector_curve_basis = compute_curve_basis_matrix(&linrels, &standard_basis).unwrap();
        let matrix_curve_basis = compute_curve_basis_matrix_from_divisor_basis_matrix(
            &linrels,
            &standard_basis,
            &divisor_basis_matrix,
        )
        .unwrap();

        assert_eq!(matrix_curve_basis, vector_curve_basis);
    }

    #[test]
    fn matrix_divisor_basis_constructs_dual_curve_basis_for_unimodular_block() {
        let linrels = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let standard_basis = vec![2, 3];
        let divisor_basis_matrix = int_matrix(&[&[0, 0, 1, 1], &[0, 0, 0, 1]]);

        let curve_basis = compute_curve_basis_matrix_from_divisor_basis_matrix(
            &linrels,
            &standard_basis,
            &divisor_basis_matrix,
        )
        .unwrap();

        assert_eq!(curve_basis, int_matrix(&[&[1, 2, 1, 0], &[0, 1, -1, 1]]));
    }

    #[test]
    fn matrix_divisor_basis_rejects_nonintegral_dual_basis() {
        let linrels = int_matrix(&[&[1, 0, -1, -1], &[0, 1, -2, -3]]);
        let standard_basis = vec![2, 3];
        let divisor_basis_matrix = int_matrix(&[&[0, 0, 2, 0], &[0, 0, 0, 1]]);

        let err = compute_curve_basis_matrix_from_divisor_basis_matrix(
            &linrels,
            &standard_basis,
            &divisor_basis_matrix,
        )
        .expect_err("non-unimodular matrix basis should be rejected");

        assert!(err.to_string().contains("integral basis"));
    }

    #[test]
    fn curve_basis_without_origin_drops_origin_column() {
        let curve_basis = int_matrix(&[&[-3, 1, 2], &[4, -5, 6]]);

        let q_matrix = curve_basis_matrix_without_origin_i64(&curve_basis).unwrap();

        assert_eq!(q_matrix, vec![vec![1, 2], vec![-5, 6]]);
    }

    #[test]
    fn curve_basis_without_origin_rejects_missing_origin_column() {
        let curve_basis = int_matrix(&[&[1], &[2]]);

        let err = curve_basis_matrix_without_origin_i64(&curve_basis)
            .expect_err("single-column curve basis should be rejected");

        assert!(err.to_string().contains("origin column"));
    }
}
