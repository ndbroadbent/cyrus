//! Divisor volume computations.
//!
//! The volume of a divisor `D_i` in a Calabi-Yau threefold is:
//! ```text
//! τ_i = (1/2) κ_ijk t^j t^k
//! ```
//!
//! The Jacobian `∂τ_i/∂t^m` is:
//! ```text
//! J_{im} = κ_ijm t^j
//! ```
//!
//! Reference: arXiv:2107.09064, Section 5

use crate::intersection::Intersection;

/// Compute divisor volumes `τ_i` = (1/2) `κ_ijk` t^j t^k.
///
/// # Panics
/// Panics if `t.len() != kappa.dim()`.
#[must_use]
pub fn compute_divisor_volumes(kappa: &Intersection, t: &[f64]) -> Vec<f64> {
    let dim = kappa.dim();
    assert_eq!(t.len(), dim, "dimension mismatch");

    // Build full tensor then contract
    // This is simpler and more correct than trying to handle symmetry manually
    let mut tau = vec![0.0; dim];

    // Use contract_triple to verify: V = (1/6) κ t³
    // τ_i = (1/2) κ_{ijk} t^j t^k = ∂V/∂t^i

    // Build the contractions manually
    for ((a, b, c), val) in kappa.iter() {
        #[allow(clippy::cast_precision_loss)]
        let val_f = val as f64;

        // Contribution to τ_a from κ_{abc} t^b t^c
        // κ_{abc} contributes to τ_a, τ_b, τ_c based on symmetry

        if a == b && b == c {
            // All same: κ_{aaa} t^a t^a - multiplicity 1
            tau[a] += 0.5 * val_f * t[a] * t[a];
        } else if a == b {
            // (a,a,c) with a < c - multiplicity 3
            // τ_a gets contribution from: κ_{aac}t^at^c (×2 from ac,ca)
            // τ_c gets contribution from: κ_{aac}t^at^a (×1)
            tau[a] += 0.5 * val_f * 2.0 * t[a] * t[c];
            tau[c] += 0.5 * val_f * t[a] * t[a];
        } else if b == c {
            // (a,b,b) with a < b - multiplicity 3
            // τ_a gets contribution from: κ_{abb}t^bt^b (×1)
            // τ_b gets contribution from: κ_{abb}t^at^b (×2 from ab,ba)
            tau[a] += 0.5 * val_f * t[b] * t[b];
            tau[b] += 0.5 * val_f * 2.0 * t[a] * t[b];
        } else {
            // All different (a,b,c) with a < b < c - multiplicity 6
            // τ_a gets κ_{abc}t^bt^c (×2 from bc,cb)
            // τ_b gets κ_{abc}t^at^c (×2 from ac,ca)
            // τ_c gets κ_{abc}t^at^b (×2 from ab,ba)
            tau[a] += 0.5 * val_f * 2.0 * t[b] * t[c];
            tau[b] += 0.5 * val_f * 2.0 * t[a] * t[c];
            tau[c] += 0.5 * val_f * 2.0 * t[a] * t[b];
        }
    }

    tau
}

/// Compute Jacobian of divisor volumes: `J_{im}` = `∂τ_i/∂t^m` = `κ_ijm` t^j.
///
/// Returns a `dim × dim` matrix where `J[i][m] = ∂τ_i/∂t^m`.
///
/// # Panics
/// Panics if `t.len() != kappa.dim()`.
#[must_use]
pub fn compute_divisor_jacobian(kappa: &Intersection, t: &[f64]) -> Vec<Vec<f64>> {
    let dim = kappa.dim();
    assert_eq!(t.len(), dim, "dimension mismatch");

    let mut jac = vec![vec![0.0; dim]; dim];

    // J[i][m] = κ_{ijm} t^j (summing over j)
    // From canonical (a,b,c) with a ≤ b ≤ c, we need to fill J entries

    for ((a, b, c), val) in kappa.iter() {
        #[allow(clippy::cast_precision_loss)]
        let val_f = val as f64;

        if a == b && b == c {
            // (a,a,a): J[a][a] += κ_{aaa} t^a
            jac[a][a] += val_f * t[a];
        } else if a == b {
            // (a,a,c) with a < c
            // J[a][a] += κ_{aac} t^c
            // J[a][c] += κ_{aac} t^a (from j=a)
            // J[c][a] += κ_{caa} t^a = κ_{aac} t^a (symmetric, from j=a)
            jac[a][a] += val_f * t[c];
            jac[a][c] += val_f * t[a];
            jac[c][a] += val_f * t[a];
        } else if b == c {
            // (a,b,b) with a < b
            // J[a][b] += κ_{abb} t^b (from j=b)
            // J[b][a] += κ_{bba} t^b = κ_{abb} t^b (symmetric)
            // J[b][b] += κ_{bba} t^a = κ_{abb} t^a
            jac[a][b] += val_f * t[b];
            jac[b][a] += val_f * t[b];
            jac[b][b] += val_f * t[a];
        } else {
            // (a,b,c) all different
            // J[a][b] += κ_{abc} t^c (from κ_{ajb} with j=c)
            // J[a][c] += κ_{abc} t^b (from κ_{ajc} with j=b)
            // J[b][a] += κ_{bac} t^c = κ_{abc} t^c
            // J[b][c] += κ_{bac} t^a = κ_{abc} t^a
            // J[c][a] += κ_{cab} t^b = κ_{abc} t^b
            // J[c][b] += κ_{cab} t^a = κ_{abc} t^a
            jac[a][b] += val_f * t[c];
            jac[a][c] += val_f * t[b];
            jac[b][a] += val_f * t[c];
            jac[b][c] += val_f * t[a];
            jac[c][a] += val_f * t[b];
            jac[c][b] += val_f * t[a];
        }
    }

    jac
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_divisor_volume_simple() {
        // κ_000 = 6, t = [2]
        // τ_0 = (1/2) × 6 × 2² = 12
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, 6);

        let t = vec![2.0];
        let tau = compute_divisor_volumes(&kappa, &t);

        assert_eq!(tau.len(), 1);
        assert!((tau[0] - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_divisor_volume_mixed() {
        // 2D case: κ_001 = 3, t = [1, 2]
        // τ_0 = (1/2) × 3 × 2 × 1 × 2 = 6 (κ_{001}t^0t^1 + κ_{010}t^1t^0)
        // τ_1 = (1/2) × 3 × 1² = 1.5 (κ_{100}t^0t^0)
        let mut kappa = Intersection::new(2);
        kappa.set(0, 0, 1, 3);

        let t = vec![1.0, 2.0];
        let tau = compute_divisor_volumes(&kappa, &t);

        assert_eq!(tau.len(), 2);
        assert!((tau[0] - 6.0).abs() < 1e-10);
        assert!((tau[1] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_jacobian_simple() {
        // κ_000 = 6, t = [2]
        // τ_0 = (1/2) × 6 × t² = 3t²
        // J_{00} = dτ_0/dt^0 = 6t = 12
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, 6);

        let t = vec![2.0];
        let jac = compute_divisor_jacobian(&kappa, &t);

        assert_eq!(jac.len(), 1);
        assert!((jac[0][0] - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_jacobian_consistency() {
        // Verify J numerically by finite differences
        let mut kappa = Intersection::new(2);
        kappa.set(0, 0, 1, 3);
        kappa.set(1, 1, 1, 2);

        let t = vec![1.5, 2.5];
        let jac = compute_divisor_jacobian(&kappa, &t);

        // Finite difference check
        let eps = 1e-6;
        for m in 0..2 {
            let mut t_plus = t.clone();
            let mut t_minus = t.clone();
            t_plus[m] += eps;
            t_minus[m] -= eps;

            let tau_plus = compute_divisor_volumes(&kappa, &t_plus);
            let tau_minus = compute_divisor_volumes(&kappa, &t_minus);

            for i in 0..2 {
                let j_numerical = (tau_plus[i] - tau_minus[i]) / (2.0 * eps);
                assert!(
                    (jac[i][m] - j_numerical).abs() < 1e-5,
                    "J[{i}][{m}] mismatch: {} vs {}",
                    jac[i][m],
                    j_numerical
                );
            }
        }
    }
}
