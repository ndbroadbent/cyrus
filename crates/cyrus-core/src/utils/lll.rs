//! LLL (Lenstra-Lenstra-Lovász) lattice basis reduction.
//!
//! Direct port of the LLL algorithm used by CYTools via flint.
//! Uses malachite for arbitrary precision integer arithmetic.

use malachite::num::arithmetic::traits::Abs;
use malachite::{Integer, Rational};

/// Result of LLL reduction with optional transformation matrix.
#[derive(Debug, Clone)]
pub struct LllResult {
    /// The reduced basis (rows are basis vectors)
    pub reduced: Vec<Vec<i64>>,
    /// Transformation matrix A such that reduced^T = A @ input^T
    /// Only present if transform=true was requested
    pub transform: Option<Vec<Vec<i64>>>,
    /// Inverse of transformation matrix
    /// Only present if transform=true was requested
    pub transform_inv: Option<Vec<Vec<i64>>>,
}

/// Apply LLL reduction to a set of points (rows).
///
/// Direct port of CYTools `utils.py` `lll_reduce` function.
///
/// Given input points as rows, computes an LLL-reduced basis.
/// Optionally returns the transformation matrix A and its inverse A_inv
/// such that reduced^T = A @ input^T.
///
/// # Arguments
///
/// * `pts` - Input points as rows of a matrix
/// * `transform` - Whether to compute and return the transformation matrix
///
/// # Returns
///
/// LllResult containing the reduced points and optionally the transformation.
///
/// # Example
///
/// ```rust
/// use cyrus_core::utils::lll_reduce;
///
/// let pts = vec![
///     vec![1, 0, 0],
///     vec![0, 1, 0],
///     vec![0, 0, 1],
/// ];
/// let result = lll_reduce(&pts, false);
/// // For an already-reduced basis, output should be similar
/// ```
pub fn lll_reduce(pts: &[Vec<i64>], transform: bool) -> LllResult {
    if pts.is_empty() {
        return LllResult {
            reduced: Vec::new(),
            transform: if transform { Some(Vec::new()) } else { None },
            transform_inv: if transform { Some(Vec::new()) } else { None },
        };
    }

    let n_rows = pts.len();
    let n_cols = pts[0].len();

    // Transpose: LLL operates on column vectors
    // pts[i][j] -> basis[j][i]
    let mut basis: Vec<Vec<Integer>> = (0..n_cols)
        .map(|j| pts.iter().map(|row| Integer::from(row[j])).collect())
        .collect();

    // Initialize transformation matrix as identity
    let mut trans: Vec<Vec<Integer>> = (0..n_cols)
        .map(|i| {
            (0..n_cols)
                .map(|j| {
                    if i == j {
                        Integer::from(1)
                    } else {
                        Integer::from(0)
                    }
                })
                .collect()
        })
        .collect();

    // Run LLL algorithm
    lll_algorithm(&mut basis, &mut trans);

    // Transpose back to row format
    let reduced: Vec<Vec<i64>> = (0..n_rows)
        .map(|i| {
            basis
                .iter()
                .map(|col| i64::try_from(&col[i]).unwrap_or(0))
                .collect()
        })
        .collect();

    if transform {
        // Convert transformation matrix
        let a: Vec<Vec<i64>> = trans
            .iter()
            .map(|row| row.iter().map(|x| i64::try_from(x).unwrap_or(0)).collect())
            .collect();

        // Compute inverse
        let a_inv = invert_integer_matrix(&a);

        LllResult {
            reduced,
            transform: Some(a),
            transform_inv: a_inv,
        }
    } else {
        LllResult {
            reduced,
            transform: None,
            transform_inv: None,
        }
    }
}

/// Core LLL algorithm implementation.
///
/// Reduces basis in-place, tracking transformation in trans.
fn lll_algorithm(basis: &mut [Vec<Integer>], trans: &mut [Vec<Integer>]) {
    let n = basis.len();
    if n == 0 {
        return;
    }
    let m = basis[0].len();

    // Gram-Schmidt coefficients (stored as rationals for exactness)
    let mut mu: Vec<Vec<Rational>> = vec![vec![Rational::from(0); n]; n];
    let mut b_star: Vec<Rational> = vec![Rational::from(0); n];

    // Initial Gram-Schmidt orthogonalization
    compute_gram_schmidt(basis, &mut mu, &mut b_star);

    // LLL parameter (typically 3/4)
    let delta = Rational::from_signeds(3, 4);

    let mut k = 1;
    while k < n {
        // Size reduction
        for j in (0..k).rev() {
            if mu[k][j].clone().abs() > Rational::from_signeds(1, 2) {
                // Round mu[k][j] to nearest integer
                let r = round_rational(&mu[k][j]);

                // Update basis[k] -= r * basis[j]
                for i in 0..m {
                    let sub = &r * &basis[j][i];
                    basis[k][i] -= sub;
                }

                // Update transformation matrix
                for i in 0..n {
                    let sub = &r * &trans[j][i];
                    trans[k][i] -= sub;
                }

                // Update Gram-Schmidt coefficients
                for i in 0..=j {
                    let sub = Rational::from(&r) * &mu[j][i];
                    mu[k][i] -= sub;
                }
            }
        }

        // Lovász condition
        let lhs = &b_star[k];
        let rhs = (&delta - &mu[k][k - 1] * &mu[k][k - 1]) * &b_star[k - 1];

        if *lhs >= rhs {
            k += 1;
        } else {
            // Swap basis[k] and basis[k-1]
            basis.swap(k, k - 1);
            trans.swap(k, k - 1);

            // Recompute Gram-Schmidt from position k-1
            recompute_gram_schmidt_from(basis, &mut mu, &mut b_star, k - 1);

            k = k.saturating_sub(1).max(1);
        }
    }
}

/// Compute Gram-Schmidt orthogonalization coefficients.
fn compute_gram_schmidt(basis: &[Vec<Integer>], mu: &mut [Vec<Rational>], b_star: &mut [Rational]) {
    let n = basis.len();

    for i in 0..n {
        // b*_i = b_i - sum_{j<i} mu[i][j] * b*_j
        // ||b*_i||^2 = <b_i, b*_i> = <b_i, b_i> - sum_{j<i} mu[i][j]^2 * ||b*_j||^2

        let mut b_star_i = dot_product(&basis[i], &basis[i]);

        for j in 0..i {
            let dot_ij = dot_product(&basis[i], &basis[j]);
            let mut sum = Rational::from(0);
            for l in 0..j {
                sum += &mu[j][l] * &mu[i][l] * &b_star[l];
            }
            if b_star[j] != 0 {
                mu[i][j] = (dot_ij - sum) / &b_star[j];
            } else {
                mu[i][j] = Rational::from(0);
            }

            b_star_i -= &mu[i][j] * &mu[i][j] * &b_star[j];
        }

        b_star[i] = b_star_i;
        mu[i][i] = Rational::from(1);
    }
}

/// Recompute Gram-Schmidt from a specific index after a swap.
fn recompute_gram_schmidt_from(
    basis: &[Vec<Integer>],
    mu: &mut [Vec<Rational>],
    b_star: &mut [Rational],
    from: usize,
) {
    let n = basis.len();

    for i in from..n {
        let mut b_star_i = dot_product(&basis[i], &basis[i]);

        for j in 0..i {
            let dot_ij = dot_product(&basis[i], &basis[j]);
            let mut sum = Rational::from(0);
            for l in 0..j {
                sum += &mu[j][l] * &mu[i][l] * &b_star[l];
            }
            if b_star[j] != 0 {
                mu[i][j] = (dot_ij - sum) / &b_star[j];
            } else {
                mu[i][j] = Rational::from(0);
            }

            b_star_i -= &mu[i][j] * &mu[i][j] * &b_star[j];
        }

        b_star[i] = b_star_i;
        mu[i][i] = Rational::from(1);
    }
}

/// Compute dot product of two integer vectors, returning Rational.
fn dot_product(a: &[Integer], b: &[Integer]) -> Rational {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| Rational::from(x * y))
        .sum()
}

/// Round a rational to the nearest integer.
fn round_rational(r: &Rational) -> Integer {
    // Add 1/2, then floor
    let half = Rational::from_signeds(1, 2);
    let shifted = if *r >= Rational::from(0) {
        r + &half
    } else {
        r - &half
    };

    // Extract floor
    use malachite::num::conversion::traits::RoundingFrom;
    use malachite::rounding_modes::RoundingMode;
    Integer::rounding_from(&shifted, RoundingMode::Floor).0
}

/// Invert an integer matrix (must be unimodular for integer inverse to exist).
fn invert_integer_matrix(a: &[Vec<i64>]) -> Option<Vec<Vec<i64>>> {
    let n = a.len();
    if n == 0 {
        return Some(Vec::new());
    }

    // Convert to rationals and perform Gaussian elimination
    let mut aug: Vec<Vec<Rational>> = a
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut r: Vec<Rational> = row.iter().map(|&x| Rational::from(x)).collect();
            // Append identity
            for j in 0..n {
                r.push(if i == j {
                    Rational::from(1)
                } else {
                    Rational::from(0)
                });
            }
            r
        })
        .collect();

    // Gaussian elimination with partial pivoting
    for col in 0..n {
        // Find pivot
        let mut pivot_row = None;
        for row in col..n {
            if aug[row][col] != 0 {
                pivot_row = Some(row);
                break;
            }
        }

        let pivot_row = pivot_row?;
        aug.swap(col, pivot_row);

        // Scale pivot row
        let pivot = aug[col][col].clone();
        for j in 0..2 * n {
            aug[col][j] /= &pivot;
        }

        // Eliminate column
        for row in 0..n {
            if row != col {
                let factor = aug[row][col].clone();
                for j in 0..2 * n {
                    let sub = &factor * &aug[col][j];
                    aug[row][j] -= sub;
                }
            }
        }
    }

    // Extract inverse and convert to integers
    let mut inv: Vec<Vec<i64>> = Vec::with_capacity(n);
    for row in &aug {
        let mut inv_row = Vec::with_capacity(n);
        for j in n..2 * n {
            // Check if it's an integer
            use malachite::num::conversion::traits::RoundingFrom;
            use malachite::rounding_modes::RoundingMode;
            let val = &row[j];
            let int_val = Integer::rounding_from(val, RoundingMode::Floor).0;
            if Rational::from(&int_val) != *val {
                return None; // Not an integer matrix
            }
            inv_row.push(i64::try_from(&int_val).ok()?);
        }
        inv.push(inv_row);
    }

    // Verify: A * A_inv should be identity
    // (Handle sign ambiguity as in CYTools)
    let mut product_is_id = true;
    let mut product_is_neg_id = true;

    for i in 0..n {
        for j in 0..n {
            let mut sum: i64 = 0;
            for k in 0..n {
                sum += a[i][k] * inv[k][j];
            }
            let expected = i64::from(i == j);
            if sum != expected {
                product_is_id = false;
            }
            if sum != -expected {
                product_is_neg_id = false;
            }
        }
    }

    if product_is_id {
        Some(inv)
    } else if product_is_neg_id {
        // Negate the inverse
        Some(
            inv.iter()
                .map(|row| row.iter().map(|&x| -x).collect())
                .collect(),
        )
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lll_identity() {
        // Identity-like basis should remain similar
        let pts = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
        let result = lll_reduce(&pts, false);
        assert_eq!(result.reduced.len(), 3);
    }

    #[test]
    fn test_lll_with_transform() {
        let pts = vec![vec![1, 0], vec![0, 1]];
        let result = lll_reduce(&pts, true);
        assert!(result.transform.is_some());
        assert!(result.transform_inv.is_some());
    }

    #[test]
    fn test_lll_reduces_basis() {
        // A basis that needs reduction
        let pts = vec![vec![1, 0], vec![1, 2]];
        let result = lll_reduce(&pts, false);

        // Check that result is a valid basis (non-empty)
        assert!(!result.reduced.is_empty());

        // Reduced basis should have shorter vectors
        let orig_norm: i64 = pts
            .iter()
            .map(|v| v.iter().map(|x| x * x).sum::<i64>())
            .sum();
        let red_norm: i64 = result
            .reduced
            .iter()
            .map(|v| v.iter().map(|x| x * x).sum::<i64>())
            .sum();
        assert!(red_norm <= orig_norm);
    }

    #[test]
    fn test_lll_transform_inverse() {
        let pts = vec![vec![1, 2], vec![3, 4]];
        let result = lll_reduce(&pts, true);

        let a = result.transform.unwrap();
        let a_inv = result.transform_inv.unwrap();

        // A * A_inv should be identity (or negative identity)
        let n = a.len();
        for i in 0..n {
            for j in 0..n {
                let mut sum: i64 = 0;
                for k in 0..n {
                    sum += a[i][k] * a_inv[k][j];
                }
                let expected = i64::from(i == j);
                assert!(sum == expected || sum == -expected);
            }
        }
    }

    #[test]
    fn test_lll_empty() {
        let pts: Vec<Vec<i64>> = Vec::new();
        let result = lll_reduce(&pts, true);
        assert!(result.reduced.is_empty());
    }

    #[test]
    fn test_round_rational() {
        assert_eq!(
            round_rational(&Rational::from_signeds(3, 2)),
            Integer::from(2)
        );
        assert_eq!(
            round_rational(&Rational::from_signeds(5, 2)),
            Integer::from(3)
        );
        assert_eq!(
            round_rational(&Rational::from_signeds(-3, 2)),
            Integer::from(-2)
        );
        assert_eq!(
            round_rational(&Rational::from_signeds(1, 3)),
            Integer::from(0)
        );
        assert_eq!(
            round_rational(&Rational::from_signeds(2, 3)),
            Integer::from(1)
        );
    }

    #[test]
    fn test_invert_integer_matrix() {
        // Simple 2x2 unimodular matrix
        let a = vec![vec![1, 1], vec![0, 1]];
        let inv = invert_integer_matrix(&a).unwrap();
        assert_eq!(inv, vec![vec![1, -1], vec![0, 1]]);
    }
}
