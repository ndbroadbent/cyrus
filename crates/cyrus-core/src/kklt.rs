//! KKLT moduli stabilization for Type IIB string compactifications.
//!
//! The KKLT mechanism stabilizes Kähler moduli by balancing flux and
//! non-perturbative terms in the superpotential.
//!
//! ## Target Divisor Volumes
//!
//! The KKLT stabilization targets specific divisor volumes (eq. 5.7):
//! ```text
//! τ_i = (c_i / 2π) × ln(W₀⁻¹)
//! ```
//!
//! Or equivalently:
//! ```text
//! τ_i = c_i / c_τ
//! ```
//!
//! Where:
//! - `c_i` = dual Coxeter numbers (1 for D3-instanton, 6 for O7-plane with so(8))
//! - `c_τ` = 2π / (`g_s` × ln(W₀⁻¹))
//!
//! Reference: arXiv:2107.09064, Section 5

use std::f64::consts::PI;

use crate::f64_pos;
use crate::intersection::Intersection;
use crate::types::f64::F64;
use crate::types::i64::I64;
use crate::types::physics::{CTau, DivisorVolume, RelativeError, StringCoupling};
use crate::types::range::CheckedRange;
use crate::types::tags::{Finite, NonNeg, Pos};

/// 2π as a typed positive constant.
const TWO_PI: F64<Pos> = f64_pos!(2.0 * PI);

/// Compute target divisor volumes `τ_i` = `c_i` / `c_τ`.
///
/// # Arguments
/// * `c_i` - Dual Coxeter numbers (positive integers: 1 for D3, 6 for O7)
/// * `c_tau` - Parameter relating `g_s` to W₀ (positive)
#[must_use]
pub fn compute_target_tau(c_i: &[I64<Pos>], c_tau: CTau) -> Vec<DivisorVolume> {
    c_i.iter()
        .map(|ci| {
            // Pos / Pos = Pos
            ci.to_f64() / c_tau
        })
        .collect()
}

/// Compute `c_τ` = 2π / (`g_s` × ln(W₀⁻¹)).
///
/// This relates the string coupling to the flux superpotential.
///
/// # Arguments
/// * `g_s` - String coupling (positive)
/// * `w0` - Flux superpotential magnitude |W₀| (positive)
#[must_use]
pub fn compute_c_tau(g_s: StringCoupling, w0: F64<Pos>) -> CTau {
    // ln(1/w0) = -ln(w0) = ln(w0.recip())
    // For KKLT we need W0 < 1, so 1/w0 > 1, so ln(1/w0) > 0
    let ln_w0_inv = w0.recip().ln();

    // 2π / (g_s × ln(1/w0))
    let denominator = g_s * ln_w0_inv;

    // This should be positive for valid KKLT
    let result = TWO_PI / denominator;
    result.try_to_pos().expect("c_tau must be positive for valid KKLT")
}

/// Compute divisor volumes `τ_i` = (1/2) `κ_ijk` t^j t^k.
///
/// This is the classical formula for divisor volumes in terms of Kähler moduli.
///
/// # Panics
/// Panics if `t.len() != kappa.dim()`.
#[must_use]
pub fn compute_divisor_volumes(kappa: &Intersection, t: &[F64<Finite>]) -> Vec<F64<Finite>> {
    let dim = kappa.dim();
    assert_eq!(t.len(), dim, "dimension mismatch");

    let mut tau = vec![F64::<Finite>::ZERO; dim];

    // For each divisor i, compute τ_i = (1/2) Σ_{j,k} κ_{ijk} t^j t^k
    for ((i_idx, j_idx, k_idx), val) in kappa.iter() {
        let val_f = val.to_f64();

        // Generate all unique permutations of (i,j,k)
        for (a, b, c) in unique_permutations(*i_idx, *j_idx, *k_idx) {
            // Finite * Finite * Finite = Finite
            tau[a] = tau[a] + val_f * t[b] * t[c];
        }
    }

    // Apply the 1/2 factor
    let half = F64::<Finite>::new(0.5).expect("0.5 is finite");
    tau.into_iter().map(|t| half * t).collect()
}

/// Compute the Jacobian `J_ik` = ∂`τ_i/∂t^k` = `κ_ijk` t^j.
///
/// # Panics
/// Panics if `t.len() != kappa.dim()`.
#[must_use]
pub fn compute_jacobian(kappa: &Intersection, t: &[F64<Finite>]) -> Vec<Vec<F64<Finite>>> {
    let dim = kappa.dim();
    assert_eq!(t.len(), dim, "dimension mismatch");

    let mut j = vec![vec![F64::<Finite>::ZERO; dim]; dim];

    for ((i_idx, j_idx, k_idx), val) in kappa.iter() {
        let val_f = val.to_f64();

        for (a, b, c) in unique_permutations(*i_idx, *j_idx, *k_idx) {
            // J[a][b] += κ * t[c], J[a][c] += κ * t[b]
            j[a][b] = j[a][b] + val_f * t[c];
            if b != c {
                j[a][c] = j[a][c] + val_f * t[b];
            }
        }
    }

    j
}

/// Generate all unique permutations of three indices.
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

/// Result of KKLT path-following.
#[derive(Debug, Clone)]
pub struct KkltResult {
    /// Kähler moduli t (can be any finite values during optimization).
    pub t: Vec<F64<Finite>>,
    /// Achieved divisor volumes.
    pub tau: Vec<F64<Finite>>,
    /// Target divisor volumes.
    pub tau_target: Vec<DivisorVolume>,
    /// Whether the solver converged.
    pub converged: bool,
    /// Relative error in tau (non-negative).
    pub relative_error: RelativeError,
}

/// Solve KKLT using path-following.
///
/// This implements the predictor-corrector algorithm from `McAllister` Section 5.2.
///
/// # Arguments
/// * `kappa` - Intersection numbers
/// * `tau_target` - Target divisor volumes (positive)
/// * `t_init` - Initial Kähler moduli guess
/// * `steps` - Checked range for interpolation steps (e.g., `range!(0..100)`)
///
/// # Returns
/// `Some(result)` if converged, `None` if diverged.
pub fn solve_path_following(
    kappa: &Intersection,
    tau_target: &[DivisorVolume],
    t_init: &[F64<Finite>],
    steps: CheckedRange<usize>,
) -> Option<KkltResult> {
    let dim = kappa.dim();
    assert_eq!(t_init.len(), dim);
    assert_eq!(tau_target.len(), dim);

    let mut t = t_init.to_vec();
    let tau_init = compute_divisor_volumes(kappa, &t);

    let n_steps = I64::<Pos>::new(steps.end as i64)?;

    for m in steps.iter_pos() {
        // m is I64<Pos>, n_steps is I64<Pos>, division yields F64<Pos>
        let alpha = m.to_f64() / n_steps.to_f64();
        let one_minus_alpha = f64_pos!(1.0) - alpha;

        // Interpolate target
        let tau_step: Vec<F64<Finite>> = tau_init
            .iter()
            .zip(tau_target.iter())
            .map(|(ti, tt)| one_minus_alpha * *ti + alpha * tt.to_finite())
            .collect();

        let tau_current = compute_divisor_volumes(kappa, &t);

        // Residual
        let delta_tau: Vec<F64<Finite>> = tau_step
            .iter()
            .zip(tau_current.iter())
            .map(|(ts, tc)| *ts - *tc)
            .collect();

        // Solve J @ epsilon = delta_tau using least squares
        let j = compute_jacobian(kappa, &t);
        let epsilon = solve_least_squares(&j, &delta_tau)?;

        // Update t
        for (ti, ei) in t.iter_mut().zip(epsilon.iter()) {
            *ti = *ti + *ei;
        }

        // Check for divergence
        let divergence_threshold = f64_pos!(1e6);
        if t.iter().any(|ti| ti.abs() > divergence_threshold) {
            return None;
        }
    }

    let tau = compute_divisor_volumes(kappa, &t);

    // Compute error: sum of squared differences
    let error_sq: F64<NonNeg> = tau
        .iter()
        .zip(tau_target.iter())
        .map(|(ta, tt)| (*ta - tt.to_finite()).square())
        .fold(F64::<NonNeg>::ZERO, |acc, x| acc + x);

    let error = error_sq.sqrt();

    // Mean of target tau values (all positive, so sum and mean are positive)
    let n_targets = I64::<Pos>::new(tau_target.len() as i64)?;
    let sum_target: F64<Pos> = tau_target
        .iter()
        .copied()
        .reduce(|acc, x| acc + x)?;
    let mean_target = sum_target / n_targets.to_f64();

    // Relative error = error / mean_target (both NonNeg/Pos, result is NonNeg)
    let relative_error = error / mean_target;

    let converged = relative_error < f64_pos!(0.001);

    Some(KkltResult {
        t,
        tau,
        tau_target: tau_target.to_vec(),
        converged,
        relative_error,
    })
}

/// Solve Ax = b using least squares (normal equations).
fn solve_least_squares(a: &[Vec<F64<Finite>>], b: &[F64<Finite>]) -> Option<Vec<F64<Finite>>> {
    let m = a.len();
    if m == 0 {
        return None;
    }
    let n = a[0].len();
    if n == 0 || b.len() != m {
        return None;
    }

    // Compute A^T A
    let mut ata = vec![vec![F64::<Finite>::ZERO; n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..m {
                ata[i][j] = ata[i][j] + a[k][i] * a[k][j];
            }
        }
    }

    // Compute A^T b
    let mut atb = vec![F64::<Finite>::ZERO; n];
    for i in 0..n {
        for k in 0..m {
            atb[i] = atb[i] + a[k][i] * b[k];
        }
    }

    // Solve (A^T A) x = A^T b using Gaussian elimination
    solve_linear_system(&ata, &atb)
}

/// Solve Ax = b using Gaussian elimination with partial pivoting.
fn solve_linear_system(a: &[Vec<F64<Finite>>], b: &[F64<Finite>]) -> Option<Vec<F64<Finite>>> {
    let n = a.len();
    if n == 0 || b.len() != n {
        return None;
    }

    // Create augmented matrix
    let mut aug: Vec<Vec<F64<Finite>>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, bi)| {
            let mut r = row.clone();
            r.push(*bi);
            r
        })
        .collect();

    // Gaussian elimination with partial pivoting
    let singular_threshold = f64_pos!(1e-14);
    for col in 0..n {
        // Find pivot (largest absolute value in column)
        let (max_row, max_abs) = aug
            .iter()
            .enumerate()
            .skip(col)
            .map(|(row, r)| (row, r[col].abs()))
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())?;

        if max_abs < singular_threshold {
            return None; // Singular matrix
        }

        // Swap rows
        aug.swap(col, max_row);

        // Eliminate column
        let pivot = aug[col][col];
        for row in (col + 1)..n {
            let factor = aug[row][col] / pivot;
            for c in col..=n {
                aug[row][c] = aug[row][c] - factor * aug[col][c];
            }
        }
    }

    // Back substitution
    let mut x = vec![F64::<Finite>::ZERO; n];
    for i in (0..n).rev() {
        let mut sum = aug[i][n];
        for (j, xj) in x.iter().enumerate().skip(i + 1) {
            sum = sum - aug[i][j] * *xj;
        }
        x[i] = sum / aug[i][i];
    }

    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::rational::Rational as TypedRational;
    use malachite::Rational;

    fn pos_i64(v: i64) -> I64<Pos> {
        I64::<Pos>::new(v).unwrap()
    }

    fn finite_f64(v: f64) -> F64<Finite> {
        F64::<Finite>::new(v).unwrap()
    }

    #[test]
    fn test_compute_c_tau() {
        // For McAllister 4-214-647:
        // g_s ≈ 0.00911, W₀ ≈ 2.3e-90, expected c_τ ≈ 3.34
        let g_s = f64_pos!(0.00911134);
        let w0 = f64_pos!(2.3e-90);
        let c_tau = compute_c_tau(g_s, w0);
        assert!((c_tau.get() - 3.34).abs() < 0.01, "c_tau = {}", c_tau.get());
    }

    #[test]
    fn test_compute_target_tau() {
        let c_i = vec![pos_i64(6), pos_i64(1), pos_i64(6), pos_i64(1)];
        let c_tau = f64_pos!(3.0);
        let tau = compute_target_tau(&c_i, c_tau);
        assert_eq!(tau.len(), 4);
        assert!((tau[0].get() - 2.0).abs() < 1e-10);
        assert!((tau[1].get() - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_divisor_volumes_simple() {
        // κ_000 = 6, t = [2]
        // τ_0 = (1/2) × 6 × 2 × 2 = 12
        let mut kappa = Intersection::new(1);
        kappa.set(
            0, 0, 0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );

        let t = vec![finite_f64(2.0)];
        let tau = compute_divisor_volumes(&kappa, &t);
        assert!((tau[0].get() - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_jacobian_simple() {
        // κ_000 = 6, t = [2]
        // ∂τ_0/∂t = 6t = 12
        let mut kappa = Intersection::new(1);
        kappa.set(
            0, 0, 0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );

        let t = vec![finite_f64(2.0)];
        let j = compute_jacobian(&kappa, &t);
        assert!((j[0][0].get() - 12.0).abs() < 1e-10, "J[0][0] = {}", j[0][0].get());
    }
}
