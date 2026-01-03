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
pub fn compute_n_matrix(kappa: &Intersection, m: &[i64]) -> Mat<f64> {
    let dim = kappa.dim();
    assert_eq!(m.len(), dim, "dimension mismatch");

    // For large tensors, parallelize the accumulation
    let contributions: Vec<_> = if kappa.num_nonzero() > 100 {
        kappa
            .par_iter()
            .flat_map(|((i, j, k), val)| {
                let (val_f, _) = f64::rounding_from(val.get(), RoundingMode::Nearest);
                unique_permutations(*i, *j, *k)
                    .map(move |(a, b, c)| (a, b, val_f * m[c] as f64))
                    .collect::<Vec<_>>()
            })
            .collect()
    } else {
        kappa
            .iter()
            .flat_map(|((i, j, k), val)| {
                let (val_f, _) = f64::rounding_from(val.get(), RoundingMode::Nearest);
                unique_permutations(*i, *j, *k)
                    .map(move |(a, b, c)| (a, b, val_f * m[c] as f64))
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

/// Compute the N matrix as `Vec<Vec<f64>>` (for backwards compatibility).
#[must_use]
pub fn compute_n_matrix_vecs(kappa: &Intersection, m: &[i64]) -> Vec<Vec<f64>> {
    let n = compute_n_matrix(kappa, m);
    let dim = n.nrows();
    (0..dim)
        .map(|i| (0..dim).map(|j| n[(i, j)]).collect())
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
#[allow(clippy::cast_precision_loss)] // Physics values fit in f64 mantissa
pub fn solve_linear_system_faer(n: &Mat<f64>, k: &[i64]) -> Option<Vec<f64>> {
    let dim = n.nrows();
    assert_eq!(n.ncols(), dim, "N must be square");
    assert_eq!(k.len(), dim, "K dimension mismatch");

    // Convert K to faer column vector
    let k_vec = Mat::from_fn(dim, 1, |i, _| k[i] as f64);

    // LU decomposition with full pivoting (most stable)
    let lu = n.full_piv_lu();

    // Solve using LU
    let p_mat = lu.solve(&k_vec);

    // Check for singularity by verifying finite solution
    for i in 0..dim {
        let val = p_mat[(i, 0)];
        if !val.is_finite() || val.abs() > 1e15 {
            // Heuristic for singular
            return None;
        }
    }

    // Extract solution
    Some((0..dim).map(|i| p_mat[(i, 0)]).collect())
}

/// Solve `N p = K` for p (legacy interface with Vec<Vec>).
///
/// Returns `None` if the matrix is singular.
///
/// # Panics
/// Panics if dimensions don't match.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn solve_linear_system(n: &[Vec<f64>], k: &[i64]) -> Option<Vec<f64>> {
    let dim = n.len();
    assert!(n.iter().all(|row| row.len() == dim), "N must be square");
    assert_eq!(k.len(), dim, "K dimension mismatch");

    // Convert to faer matrix
    let n_mat = Mat::from_fn(dim, dim, |i, j| n[i][j]);

    solve_linear_system_faer(&n_mat, k)
}

/// Compute flat direction p = N⁻¹ K where `N_ab = κ_abc M^c`.
///
/// Returns `None` if N is singular.
///
/// # Panics
/// Panics if dimensions don't match.
#[must_use]
pub fn compute_flat_direction(kappa: &Intersection, k: &[i64], m: &[i64]) -> Option<Vec<f64>> {
    let n = compute_n_matrix(kappa, m);
    solve_linear_system_faer(&n, k)
}

/// Compute `e^{K₀} = (4/3 × κ_abc p^a p^b p^c)^{-1}`.
///
/// Returns `F64<Pos>` since `e^x` is always positive.
///
/// # Errors
/// Returns an error if the contraction is non-positive (invalid flat direction).
pub fn compute_ek0(kappa: &Intersection, p: &[f64]) -> Result<F64<Pos>> {
    // Convert raw values to typed (flat direction can have any sign)
    let p_typed: Vec<F64<Finite>> = p
        .iter()
        .map(|&x| {
            F64::<Finite>::new(x).expect("flat direction components must be finite")
        })
        .collect();

    // Boundary: contraction must be positive for valid physics
    let kappa_ppp = kappa.contract_triple_finite(&p_typed).ok_or_else(|| {
        Error::InvalidInput(
            "κ_abc p^a p^b p^c must be positive for valid flat direction".into(),
        )
    })?;

    // Type algebra: Pos / Pos = Pos
    Ok(f64_pos!(1.0) / (f64_pos!(4.0 / 3.0) * kappa_ppp))
}

/// Result of flat direction computation.
#[derive(Debug, Clone)]
pub struct FlatDirectionResult {
    /// N matrix as Vec<Vec> for compatibility: `N_ab = κ_abc M^c`.
    pub n_matrix: Vec<Vec<f64>>,
    /// Flat direction: p = N⁻¹ K.
    pub p: Vec<f64>,
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
pub fn compute_flat_direction_full(
    kappa: &Intersection,
    k: &[i64],
    m: &[i64],
) -> Option<FlatDirectionResult> {
    let n_mat = compute_n_matrix(kappa, m);
    let p = solve_linear_system_faer(&n_mat, k)?;
    let ek0 = compute_ek0(kappa, &p).ok()?;

    // Convert to Vec<Vec> for result
    let dim = n_mat.nrows();
    let n_matrix = (0..dim)
        .map(|i| (0..dim).map(|j| n_mat[(i, j)]).collect())
        .collect();

    Some(FlatDirectionResult { n_matrix, p, ek0 })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::rational::Rational as TypedRational;
use crate::types::tags::Pos;
    use malachite::Rational;

    #[test]
    fn test_compute_n_matrix_simple() {
        // κ_000 = 6 means N_00 = 6 * M^0
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, TypedRational::<Pos>::new(Rational::from(6)).unwrap());

        let m = vec![2];
        let n = compute_n_matrix(&kappa, &m);

        assert_eq!(n.nrows(), 1);
        assert!((n[(0, 0)] - 12.0).abs() < 1e-10); // 6 * 2 = 12
    }

    #[test]
    fn test_solve_linear_system() {
        // Simple 2x2: [2, 1; 1, 3] x = [5, 5]
        // Solution: x = [2, 1]
        let n = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let k = vec![5, 5];

        let p = solve_linear_system(&n, &k).unwrap();

        assert!((p[0] - 2.0).abs() < 1e-10);
        assert!((p[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_flat_direction() {
        // With κ_000 = 6, M = [1], K = [3]:
        // N_00 = 6, so p = 6⁻¹ × 3 = 0.5
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, TypedRational::<Pos>::new(Rational::from(6)).unwrap());

        let k = vec![3];
        let m = vec![1];

        let p = compute_flat_direction(&kappa, &k, &m).unwrap();
        assert!((p[0] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_compute_ek0() {
        // κ_000 = 6, p = [1]
        // κ_ppp = 6 × 1 = 6 (note: multiplicity 1 for iii)
        // e^K₀ = 1 / (4/3 × 6) = 1/8 = 0.125
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, TypedRational::<Pos>::new(Rational::from(6)).unwrap());

        let p = vec![1.0];
        let ek0 = compute_ek0(&kappa, &p).unwrap();

        assert!((ek0.get() - 0.125).abs() < 1e-10);
    }
}
