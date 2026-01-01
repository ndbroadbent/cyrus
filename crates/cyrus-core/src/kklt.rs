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
//! ## Path-Following Algorithm
//!
//! Given target τ values, we solve for Kähler moduli t using path-following:
//! 1. Start with `t_init` (e.g., uniform or from triangulation)
//! 2. Interpolate from `τ_init` to `τ_target` in N steps
//! 3. At each step, solve the linear system: J @ ε = Δτ
//! 4. Update t = t + ε
//!
//! The Jacobian is: `J_ik` = ∂`τ_i/∂t^k` = `κ_ijk` t^j
//!
//! Reference: arXiv:2107.09064, Section 5

use crate::intersection::Intersection;

/// Compute target divisor volumes `τ_i` = `c_i` / `c_τ`.
///
/// # Arguments
/// * `c_i` - Dual Coxeter numbers (1 for D3, 6 for O7)
/// * `c_tau` - Parameter relating `g_s` to W₀
#[must_use]
#[allow(clippy::cast_possible_truncation)]
pub fn compute_target_tau(c_i: &[i64], c_tau: f64) -> Vec<f64> {
    c_i.iter().map(|&ci| f64::from(ci as i32) / c_tau).collect()
}

/// Compute `c_τ` = 2π / (`g_s` × ln(W₀⁻¹)).
///
/// This relates the string coupling to the flux superpotential.
#[must_use]
pub fn compute_c_tau(g_s: f64, w0: f64) -> f64 {
    use std::f64::consts::PI;
    2.0 * PI / (g_s * (1.0 / w0).ln())
}

/// Compute divisor volumes `τ_i` = (1/2) `κ_ijk` t^j t^k.
///
/// This is the classical formula for divisor volumes in terms of Kähler moduli.
///
/// # Panics
/// Panics if `t.len() != kappa.dim()`.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn compute_divisor_volumes(kappa: &Intersection, t: &[f64]) -> Vec<f64> {
    let dim = kappa.dim();
    assert_eq!(t.len(), dim, "dimension mismatch");

    let mut tau = vec![0.0; dim];

    // For each divisor i, compute τ_i = (1/2) Σ_{j,k} κ_{ijk} t^j t^k
    // Using the fact that κ is symmetric, we iterate over canonical entries
    // and account for all permutations.
    for ((i_idx, j_idx, k_idx), val) in kappa.iter() {
        let val_f = val as f64;

        // Generate all unique permutations of (i,j,k)
        // For each permutation (a,b,c), we add val * t[b] * t[c] to τ[a]
        for (a, b, c) in unique_permutations(i_idx, j_idx, k_idx) {
            tau[a] += val_f * t[b] * t[c];
        }
    }

    // Apply the 1/2 factor
    for tau_i in &mut tau {
        *tau_i *= 0.5;
    }

    tau
}

/// Compute the Jacobian `J_ik` = ∂`τ_i/∂t^k` = `κ_ijk` t^j.
///
/// # Panics
/// Panics if `t.len() != kappa.dim()`.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn compute_jacobian(kappa: &Intersection, t: &[f64]) -> Vec<Vec<f64>> {
    let dim = kappa.dim();
    assert_eq!(t.len(), dim, "dimension mismatch");

    let mut j = vec![vec![0.0; dim]; dim];

    // For each κ entry, compute contribution to Jacobian
    // J_ik = ∂τ_i/∂t^k = Σ_j κ_ijk t^j
    for ((i_idx, j_idx, k_idx), val) in kappa.iter() {
        let val_f = val as f64;

        // For each permutation (a,b,c):
        // τ_a = (1/2) κ_{abc} t^b t^c
        // ∂τ_a/∂t^k = κ_{abc} × (δ_{bk} t^c + t^b δ_{ck}) / 2
        //           = κ_{abc} t^c when b=k, plus κ_{abc} t^b when c=k
        // Due to symmetry in b,c: ∂τ_a/∂t^k = κ_{akc} t^c

        for (a, b, c) in unique_permutations(i_idx, j_idx, k_idx) {
            // When computing ∂τ_a/∂t^k, we need contributions where
            // one of b,c equals k and we sum over the other.
            // Due to symmetry: J_ak = κ_abk t^b + κ_akc t^c
            // But κ_abk = κ_akb so this is 2 × κ_akc t^c... wait, that's wrong.

            // Let's be more careful. τ_a = (1/2) Σ_{b,c} κ_{abc} t^b t^c
            // ∂τ_a/∂t^k = (1/2) Σ_c κ_{akc} t^c + (1/2) Σ_b κ_{abk} t^b
            //           = (1/2) Σ_c κ_{akc} t^c + (1/2) Σ_b κ_{abk} t^b
            // Since κ is symmetric: κ_{abk} = κ_{akb}
            // So this is: (1/2) Σ_c κ_{akc} t^c + (1/2) Σ_c κ_{akc} t^c = Σ_c κ_{akc} t^c

            // So J_ak = Σ_c κ_{akc} t^c (no factor of 1/2 after applying chain rule)

            // For permutation (a,b,c), this contributes to J[a][b] += κ * t[c]
            // and to J[a][c] += κ * t[b]
            j[a][b] += val_f * t[c];
            if b != c {
                j[a][c] += val_f * t[b];
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
    /// Kähler moduli t.
    pub t: Vec<f64>,
    /// Achieved divisor volumes.
    pub tau: Vec<f64>,
    /// Target divisor volumes.
    pub tau_target: Vec<f64>,
    /// Whether the solver converged.
    pub converged: bool,
    /// Relative error in tau.
    pub relative_error: f64,
}

/// Solve KKLT using path-following.
///
/// This implements the predictor-corrector algorithm from `McAllister` Section 5.2.
///
/// # Arguments
/// * `kappa` - Intersection numbers
/// * `tau_target` - Target divisor volumes
/// * `t_init` - Initial Kähler moduli guess
/// * `n_steps` - Number of interpolation steps
///
/// # Returns
/// `Some(result)` if converged, `None` if diverged.
///
/// # Panics
/// Panics if `t_init` or `tau_target` length doesn't match `kappa.dim()`.
#[allow(clippy::cast_precision_loss)]
pub fn solve_path_following(
    kappa: &Intersection,
    tau_target: &[f64],
    t_init: &[f64],
    n_steps: usize,
) -> Option<KkltResult> {
    let dim = kappa.dim();
    assert_eq!(t_init.len(), dim);
    assert_eq!(tau_target.len(), dim);

    let mut t = t_init.to_vec();
    let tau_init = compute_divisor_volumes(kappa, &t);

    for m in 0..n_steps {
        let alpha = (m + 1) as f64 / n_steps as f64;

        // Interpolate target
        let tau_step: Vec<f64> = tau_init
            .iter()
            .zip(tau_target.iter())
            .map(|(&ti, &tt)| (1.0 - alpha).mul_add(ti, alpha * tt))
            .collect();

        let tau_current = compute_divisor_volumes(kappa, &t);

        // Residual
        let delta_tau: Vec<f64> = tau_step
            .iter()
            .zip(tau_current.iter())
            .map(|(&ts, &tc)| ts - tc)
            .collect();

        // Solve J @ epsilon = delta_tau using least squares
        let j = compute_jacobian(kappa, &t);
        let epsilon = solve_least_squares(&j, &delta_tau)?;

        // Update t
        for (ti, ei) in t.iter_mut().zip(epsilon.iter()) {
            *ti += ei;
        }

        // Check for divergence
        if t.iter().any(|ti| ti.abs() > 1e6) {
            return None;
        }
    }

    let tau = compute_divisor_volumes(kappa, &t);

    // Compute error
    let error: f64 = tau
        .iter()
        .zip(tau_target.iter())
        .map(|(&ta, &tt)| (ta - tt).powi(2))
        .sum::<f64>()
        .sqrt();
    let mean_target: f64 = tau_target.iter().sum::<f64>() / tau_target.len() as f64;
    let relative_error = error / mean_target;
    let converged = relative_error < 0.001;

    Some(KkltResult {
        t,
        tau,
        tau_target: tau_target.to_vec(),
        converged,
        relative_error,
    })
}

/// Solve Ax = b using least squares (truncated SVD).
#[allow(clippy::needless_range_loop)]
fn solve_least_squares(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let m = a.len();
    if m == 0 {
        return None;
    }
    let n = a[0].len();
    if n == 0 || b.len() != m {
        return None;
    }

    // Simple least squares using normal equations: (A^T A) x = A^T b
    // This is less numerically stable than SVD but simpler to implement.
    // TODO: Use SVD or QR for better numerical stability.

    // Compute A^T A
    let mut ata = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            for k in 0..m {
                ata[i][j] += a[k][i] * a[k][j];
            }
        }
    }

    // Compute A^T b
    let mut atb = vec![0.0; n];
    for i in 0..n {
        for k in 0..m {
            atb[i] += a[k][i] * b[k];
        }
    }

    // Solve (A^T A) x = A^T b using Gaussian elimination
    solve_linear_system(&ata, &atb)
}

/// Solve Ax = b using Gaussian elimination with partial pivoting.
fn solve_linear_system(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = a.len();
    if n == 0 || b.len() != n {
        return None;
    }

    // Create augmented matrix
    let mut aug: Vec<Vec<f64>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, &bi)| {
            let mut r = row.clone();
            r.push(bi);
            r
        })
        .collect();

    // Gaussian elimination with partial pivoting
    for col in 0..n {
        // Find pivot
        let (max_row, max_val) = aug
            .iter()
            .enumerate()
            .skip(col)
            .map(|(row, r)| (row, r[col].abs()))
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())?;

        if max_val < 1e-14 {
            return None; // Singular matrix
        }

        // Swap rows
        aug.swap(col, max_row);

        // Eliminate column
        let pivot = aug[col][col];
        for row in (col + 1)..n {
            let factor = aug[row][col] / pivot;
            for c in col..=n {
                aug[row][c] -= factor * aug[col][c];
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = aug[i][n];
        for (j, &xj) in x.iter().enumerate().skip(i + 1) {
            sum -= aug[i][j] * xj;
        }
        x[i] = sum / aug[i][i];
    }

    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_c_tau() {
        // For McAllister 4-214-647:
        // g_s ≈ 0.00911, W₀ ≈ 2.3e-90, expected c_τ ≈ 3.34
        let g_s = 0.00911134;
        let w0 = 2.3e-90;
        let c_tau = compute_c_tau(g_s, w0);
        assert!((c_tau - 3.34).abs() < 0.01, "c_tau = {c_tau}");
    }

    #[test]
    fn test_compute_target_tau() {
        let c_i = vec![6, 1, 6, 1];
        let c_tau = 3.0;
        let tau = compute_target_tau(&c_i, c_tau);
        assert_eq!(tau.len(), 4);
        assert!((tau[0] - 2.0).abs() < 1e-10);
        assert!((tau[1] - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_divisor_volumes_simple() {
        // κ_000 = 6, t = [2]
        // τ_0 = (1/2) × 6 × 2 × 2 = 12
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, 6);

        let t = vec![2.0];
        let tau = compute_divisor_volumes(&kappa, &t);
        assert!((tau[0] - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_jacobian_simple() {
        // κ_000 = 6, t = [2]
        // τ_0 = (1/2) × 6 × t × t = 3t²
        // ∂τ_0/∂t = 6t = 12
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, 6);

        let t = vec![2.0];
        let j = compute_jacobian(&kappa, &t);
        assert!((j[0][0] - 12.0).abs() < 1e-10, "J[0][0] = {}", j[0][0]);
    }
}
