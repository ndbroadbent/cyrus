//! Flat direction computation for flux compactifications.
//!
//! Given flux vectors K and M and intersection numbers `κ_ijk`,
//! compute the flat direction p and Kähler potential factor `e^{K₀}`.
//!
//! The flat direction satisfies:
//! ```text
//! N_ab p^b = K_a  where  N_ab = κ_abc M^c
//! ```
//!
//! The Kähler potential factor is:
//! ```text
//! e^{K₀} = (4/3 × κ_abc p^a p^b p^c)^{-1}
//! ```
//!
//! ## Performance
//!
//! - Uses `faer` for optimized linear algebra (BLAS-like performance)
//! - Parallel N matrix construction via rayon
//! - GPU-ready: faer matrices can be ported to GPU backends
//!
//! Reference: arXiv:2107.09064, Section 6

use faer::{Mat, prelude::SpSolver};
use malachite::num::conversion::traits::RoundingFrom;
use malachite::rounding_modes::RoundingMode;
use rayon::prelude::*;

use crate::error::{Error, Result};
use crate::f64_pos;
use crate::intersection::Intersection;
use crate::types::f64::F64;
use crate::types::i64::I64;
use crate::types::range::CheckedRange;
use crate::types::tags::{Finite, Pos};

/// Compute the N matrix: `N_ab = κ_abc M^c`.
///
/// Contracts the intersection tensor with flux vector M.
/// Returns a faer `Mat<f64>` for efficient linear algebra.
///
/// # Panics
/// Panics if `m.len() != kappa.dim()`.
#[must_use]
#[allow(clippy::cast_precision_loss)] // Physics values fit in f64 mantissa
pub fn compute_n_matrix(kappa: &Intersection, m: &[I64<Finite>]) -> Mat<f64> {
    let dim = kappa.dim();
    assert_eq!(m.len(), dim, "dimension mismatch");

    // For large tensors, parallelize the accumulation
    let contributions: Vec<_> = if kappa.num_nonzero() > 100 {
        kappa
            .par_iter()
            .flat_map(|((i, j, k), val)| {
                let (val_f, _) = f64::rounding_from(val.get(), RoundingMode::Nearest);
                unique_permutations(*i, *j, *k)
                    .map(move |(a, b, c)| (a, b, val_f * m[c].get() as f64))
                    .collect::<Vec<_>>()
            })
            .collect()
    } else {
        kappa
            .iter()
            .flat_map(|((i, j, k), val)| {
                let (val_f, _) = f64::rounding_from(val.get(), RoundingMode::Nearest);
                unique_permutations(*i, *j, *k)
                    .map(move |(a, b, c)| (a, b, val_f * m[c].get() as f64))
                    .collect::<Vec<_>>()
            })
            .collect()
    };

    // Accumulate into matrix
    let mut n = Mat::zeros(dim, dim);
    for (a, b, contrib) in contributions {
        n[(a, b)] += contrib;
    }

    n
}

/// Compute the N matrix as `Vec<Vec<F64<Finite>>>`.
#[must_use]
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub fn compute_n_matrix_vecs(kappa: &Intersection, m: &[I64<Finite>]) -> Vec<Vec<F64<Finite>>> {
    let n = compute_n_matrix(kappa, m);
    let dim = n.nrows();
    let rows = CheckedRange::new(0, dim);
    let cols = CheckedRange::new(0, dim);

    rows.iter_non_neg()
        .map(|i| {
            let i_idx = i.get() as usize;
            cols.iter_non_neg()
                .map(|j| {
                    let j_idx = j.get() as usize;
                    // Boundary: faer uses raw f64, convert to typed
                    F64::from_raw(n[(i_idx, j_idx)])
                })
                .collect()
        })
        .collect()
}

/// Generate all unique permutations of three indices.
///
/// Returns an iterator over `(a, b, c)` tuples representing unique permutations.
fn unique_permutations(
    i: usize,
    j: usize,
    k: usize,
) -> impl Iterator<Item = (usize, usize, usize)> {
    let perms = [
        (i, j, k),
        (i, k, j),
        (j, i, k),
        (j, k, i),
        (k, i, j),
        (k, j, i),
    ];

    // Deduplicate based on which indices are equal
    let mut seen = [(false, (0, 0, 0)); 6];
    let mut count = 0;

    for perm in perms {
        let is_dup = seen[..count].iter().any(|(_, p)| *p == perm);
        if !is_dup {
            seen[count] = (true, perm);
            count += 1;
        }
    }

    seen.into_iter()
        .take(count)
        .filter(|(used, _)| *used)
        .map(|(_, perm)| perm)
}

/// Solve `N p = K` for p using faer's LU decomposition.
///
/// Uses full pivoting LU for numerical stability.
/// Returns `None` if the matrix is singular.
///
/// # Panics
/// Panics if dimensions don't match.
#[must_use]
#[allow(
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss
)]
pub fn solve_linear_system_faer(n: &Mat<f64>, k: &[I64<Finite>]) -> Option<Vec<F64<Finite>>> {
    let dim = n.nrows();
    assert_eq!(n.ncols(), dim, "N must be square");
    assert_eq!(k.len(), dim, "K dimension mismatch");

    let indices = CheckedRange::new(0, dim);

    // Convert K to faer column vector (boundary: typed → raw)
    let k_vec = Mat::from_fn(dim, 1, |i, _| k[i].get() as f64);

    // LU decomposition with full pivoting (most stable)
    let lu = n.full_piv_lu();

    // Solve using LU
    let p_mat = lu.solve(&k_vec);

    // Check for singularity by verifying finite solution
    let singularity_threshold = f64_pos!(1e15);
    for i in indices.iter_non_neg() {
        let val = p_mat[(i.get() as usize, 0)];
        if !val.is_finite() || val.abs() > singularity_threshold.get() {
            return None;
        }
    }

    // Extract solution (boundary: raw → typed)
    Some(
        CheckedRange::new(0, dim)
            .iter_non_neg()
            .map(|i| F64::from_raw(p_mat[(i.get() as usize, 0)]))
            .collect(),
    )
}

/// Solve `N p = K` for p (interface with typed Vec<Vec>).
///
/// Returns `None` if the matrix is singular.
///
/// # Panics
/// Panics if dimensions don't match.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn solve_linear_system(n: &[Vec<F64<Finite>>], k: &[I64<Finite>]) -> Option<Vec<F64<Finite>>> {
    let dim = n.len();
    assert!(n.iter().all(|row| row.len() == dim), "N must be square");
    assert_eq!(k.len(), dim, "K dimension mismatch");

    // Convert to faer matrix (boundary: typed → raw)
    let n_mat = Mat::from_fn(dim, dim, |i, j| n[i][j].get());

    solve_linear_system_faer(&n_mat, k)
}

/// Compute flat direction p = N⁻¹ K where `N_ab = κ_abc M^c`.
///
/// Returns `None` if N is singular.
///
/// # Panics
/// Panics if dimensions don't match.
#[must_use]
pub fn compute_flat_direction(
    kappa: &Intersection,
    k: &[I64<Finite>],
    m: &[I64<Finite>],
) -> Option<Vec<F64<Finite>>> {
    let n = compute_n_matrix(kappa, m);
    solve_linear_system_faer(&n, k)
}

/// Compute `e^{K₀} = (4/3 × κ_abc p^a p^b p^c)^{-1}`.
///
/// Returns `F64<Pos>` since `e^x` is always positive.
///
/// # Errors
/// Returns an error if the contraction is non-positive (invalid flat direction).
pub fn compute_ek0(kappa: &Intersection, p: &[F64<Finite>]) -> Result<F64<Pos>> {
    // Boundary: contraction must be positive for valid physics
    let kappa_ppp = kappa.contract_triple_finite(p).ok_or_else(|| {
        Error::InvalidInput("κ_abc p^a p^b p^c must be positive for valid flat direction".into())
    })?;

    // Type algebra: Pos / Pos = Pos
    Ok(f64_pos!(1.0) / (f64_pos!(4.0 / 3.0) * kappa_ppp))
}

/// Result of flat direction computation.
#[derive(Debug, Clone)]
pub struct FlatDirectionResult {
    /// N matrix: `N_ab = κ_abc M^c`.
    pub n_matrix: Vec<Vec<F64<Finite>>>,
    /// Flat direction: p = N⁻¹ K.
    pub p: Vec<F64<Finite>>,
    /// Kähler potential factor: `e^{K₀} = (4/3 × κ_abc p^a p^b p^c)^{-1}`.
    /// Always positive (it's an exponential).
    pub ek0: F64<Pos>,
}

/// Compute flat direction with full breakdown.
///
/// Returns `None` if N is singular or if the flat direction is invalid.
///
/// # Panics
/// Panics if dimensions don't match.
#[must_use]
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
pub fn compute_flat_direction_full(
    kappa: &Intersection,
    k: &[I64<Finite>],
    m: &[I64<Finite>],
) -> Option<FlatDirectionResult> {
    let n_mat = compute_n_matrix(kappa, m);
    let p = solve_linear_system_faer(&n_mat, k)?;
    let ek0 = compute_ek0(kappa, &p).ok()?;

    // Convert to Vec<Vec<F64<Finite>>>
    let dim = n_mat.nrows();
    let indices = CheckedRange::new(0, dim);
    let n_matrix = indices
        .iter_non_neg()
        .map(|i| {
            let i_idx = i.get() as usize;
            CheckedRange::new(0, dim)
                .iter_non_neg()
                .map(|j| F64::from_raw(n_mat[(i_idx, j.get() as usize)]))
                .collect()
        })
        .collect();

    Some(FlatDirectionResult { n_matrix, p, ek0 })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::rational::Rational as TypedRational;
    use malachite::Rational;

    fn finite_i64(v: i64) -> I64<Finite> {
        I64::<Finite>::new(v)
    }

    fn finite_f64(v: f64) -> F64<Finite> {
        F64::<Finite>::new(v).unwrap()
    }

    #[test]
    fn test_compute_n_matrix_simple() {
        // κ_000 = 6 means N_00 = 6 * M^0
        let mut kappa = Intersection::new(1);
        kappa.set(
            0,
            0,
            0,
            TypedRational::<Pos>::new(Rational::from(6)).unwrap(),
        );

        let m = vec![finite_i64(2)];
        let n = compute_n_matrix(&kappa, &m);

        assert_eq!(n.nrows(), 1);
        assert!((n[(0, 0)] - 12.0).abs() < 1e-10); // 6 * 2 = 12
    }

    #[test]
    fn test_solve_linear_system() {
        // Simple 2x2: [2, 1; 1, 3] x = [5, 5]
        // Solution: x = [2, 1]
        let n = vec![
            vec![finite_f64(2.0), finite_f64(1.0)],
            vec![finite_f64(1.0), finite_f64(3.0)],
        ];
        let k = vec![finite_i64(5), finite_i64(5)];

        let p = solve_linear_system(&n, &k).unwrap();

        assert!((p[0].get() - 2.0).abs() < 1e-10);
        assert!((p[1].get() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_flat_direction() {
        // With κ_000 = 6, M = [1], K = [3]:
        // N_00 = 6, so p = 6⁻¹ × 3 = 0.5
        let mut kappa = Intersection::new(1);
        kappa.set(
            0,
            0,
            0,
            TypedRational::<Pos>::new(Rational::from(6)).unwrap(),
        );

        let k = vec![finite_i64(3)];
        let m = vec![finite_i64(1)];

        let p = compute_flat_direction(&kappa, &k, &m).unwrap();
        assert!((p[0].get() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_compute_ek0() {
        // κ_000 = 6, p = [1]
        // κ_ppp = 6 × 1 = 6 (note: multiplicity 1 for iii)
        // e^K₀ = 1 / (4/3 × 6) = 1/8 = 0.125
        let mut kappa = Intersection::new(1);
        kappa.set(
            0,
            0,
            0,
            TypedRational::<Pos>::new(Rational::from(6)).unwrap(),
        );

        let p = vec![finite_f64(1.0)];
        let ek0 = compute_ek0(&kappa, &p).unwrap();

        assert!((ek0.get() - 0.125).abs() < 1e-10);
    }
}
