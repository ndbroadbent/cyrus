//! Racetrack moduli stabilization for flux compactifications.
//!
//! The racetrack mechanism uses multiple non-perturbative terms to stabilize
//! the axio-dilaton τ and compute the flux superpotential W₀.
//!
//! ## The Racetrack Superpotential
//!
//! From `McAllister` eq. 2.22:
//! ```text
//! W = -ζ Σ_q (M·q) N_q Li₂(e^{2πiτ(q·p)})
//! ```
//!
//! Where:
//! - ζ = 1 / (2^{3/2} π^{5/2}) ≈ 0.0114
//! - M = flux vector
//! - `N_q` = Gopakumar-Vafa invariants
//! - q = curve classes
//! - p = flat direction vector
//! - τ = axio-dilaton
//!
//! ## F-term Stabilization
//!
//! The F-term ∂W/∂τ = 0 is solved using the two leading terms:
//! ```text
//! (M·q₁)(p·q₁)N_{q₁} e^{2πiτα} + (M·q₂)(p·q₂)N_{q₂} e^{2πiτβ} = 0
//! ```
//!
//! This gives δ (eq. 2.25):
//! ```text
//! δ = -[(M·q₁)(p·q₁)N_{q₁}] / [(M·q₂)(p·q₂)N_{q₂}]
//! ```
//!
//! And the solution (eq. 2.26):
//! ```text
//! Im(τ) = ln(1/|δ|) / (2πε)  where ε = β - α
//! g_s = 1 / Im(τ)
//! ```
//!
//! Reference: arXiv:2107.09064, Section 2.4

use std::f64::consts::PI;

/// Constant ζ = 1 / (2^{3/2} π^{5/2}) ≈ 0.02021.
pub const ZETA: f64 = 0.020_210_652_027_623_285;

/// A racetrack term representing one contribution to the superpotential.
///
/// Each term has an exponent (q·p) and coefficient ((M·q) × `N_q`).
#[derive(Debug, Clone, Copy)]
pub struct RacetrackTerm {
    /// Exponent: q·p where q is curve class and p is flat direction.
    pub exponent: f64,
    /// Coefficient: (M·q) × `N_q` where M is flux vector and `N_q` is GV invariant.
    pub coefficient: i64,
}

/// Result of racetrack stabilization.
#[derive(Debug, Clone)]
pub struct RacetrackResult {
    /// String coupling `g_s` = 1/Im(τ).
    pub g_s: f64,
    /// Real part of τ at minimum (0 for opposite-sign coefficients).
    pub re_tau: f64,
    /// Imaginary part of τ at minimum.
    pub im_tau: f64,
    /// Parameter δ from eq. 2.25.
    pub delta: f64,
    /// Exponent gap ε = β - α.
    pub epsilon: f64,
}

/// Solve the racetrack for `g_s` using the two leading terms.
///
/// The terms should be sorted by exponent (ascending). Uses the first two
/// terms with non-zero coefficients to solve for τ.
///
/// # Returns
/// Returns `Some(result)` if stabilization succeeded, `None` if:
/// - Fewer than 2 terms with non-zero coefficients
/// - Invalid exponent ordering (β ≤ α)
/// - Im(τ) ≤ 0 (no valid stabilization)
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn solve_racetrack(terms: &[RacetrackTerm]) -> Option<RacetrackResult> {
    // Filter to terms with non-zero coefficients
    let valid_terms: Vec<_> = terms.iter().filter(|t| t.coefficient != 0).collect();

    if valid_terms.len() < 2 {
        return None;
    }

    let term1 = valid_terms[0];
    let term2 = valid_terms[1];

    let alpha = term1.exponent;
    let beta = term2.exponent;

    // Exponent gap
    let epsilon = beta - alpha;
    if epsilon <= 0.0 {
        return None;
    }

    // Derivative coefficients: (M·q)(p·q)N_q = coeff × exponent
    let d1 = (term1.coefficient as f64) * alpha;
    let d2 = (term2.coefficient as f64) * beta;

    if d2.abs() < 1e-10 {
        return None;
    }

    // δ = -d1/d2
    let delta = -d1 / d2;

    if delta == 0.0 {
        return None;
    }

    // Im(τ) = ln(1/|δ|) / (2πε)
    let im_tau = (1.0 / delta.abs()).ln() / (2.0 * PI * epsilon);

    if im_tau <= 0.0 {
        return None;
    }

    // Re(τ) = 1/(2ε) if δ < 0 (same-sign coefficients), else 0
    let re_tau = if delta < 0.0 {
        1.0 / (2.0 * epsilon)
    } else {
        0.0
    };

    let g_s = 1.0 / im_tau;

    Some(RacetrackResult {
        g_s,
        re_tau,
        im_tau,
        delta,
        epsilon,
    })
}

/// Compute the relationship `c_τ` = 2π / (`g_s` × ln(1/W₀)).
///
/// This is eq. 2.29 from `McAllister`, relating the string coupling to W₀.
#[must_use]
pub fn compute_c_tau(g_s: f64, w0_log_inv: f64) -> f64 {
    2.0 * PI / (g_s * w0_log_inv)
}

/// Compute `g_s` from `c_τ` and ln(1/W₀).
#[must_use]
pub fn compute_g_s_from_c_tau(c_tau: f64, w0_log_inv: f64) -> f64 {
    2.0 * PI / (c_tau * w0_log_inv)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[allow(clippy::suboptimal_flops)]
    fn test_zeta_constant() {
        // ζ = 1 / (2^{3/2} π^{5/2}) ≈ 0.02021
        let expected = 1.0 / (2.0_f64.powf(1.5) * PI.powf(2.5));
        assert!((ZETA - expected).abs() < 1e-15);
        assert!((ZETA - 0.02021).abs() < 0.0001);
    }

    #[test]
    fn test_solve_racetrack_edge_case() {
        // When δ = 1, ln(1/|δ|) = 0, so Im(τ) = 0, which means no stabilization
        let terms = vec![
            RacetrackTerm {
                exponent: 1.0,
                coefficient: 10,
            },
            RacetrackTerm {
                exponent: 2.0,
                coefficient: -5,
            },
        ];

        // δ = -(10×1) / (-5×2) = 1, so this should return None
        let result = solve_racetrack(&terms);
        assert!(
            result.is_none(),
            "Should fail when δ = 1 (no stabilization)"
        );
    }

    #[test]
    fn test_solve_racetrack_valid() {
        // Create two terms where |δ| < 1
        let terms = vec![
            RacetrackTerm {
                exponent: 1.0,
                coefficient: 2,
            },
            RacetrackTerm {
                exponent: 2.0,
                coefficient: -10,
            },
        ];

        let result = solve_racetrack(&terms).unwrap();

        // δ = -(2×1) / (-10×2) = -2 / -20 = 0.1
        assert!((result.delta - 0.1).abs() < 1e-10);

        // ε = 2 - 1 = 1
        assert!((result.epsilon - 1.0).abs() < 1e-10);

        // Im(τ) = ln(10) / (2π) ≈ 0.366
        let expected_im_tau = 10.0_f64.ln() / (2.0 * PI);
        assert!((result.im_tau - expected_im_tau).abs() < 1e-6);

        // g_s = 1/Im(τ)
        assert!((result.g_s - 1.0 / expected_im_tau).abs() < 1e-6);
    }

    #[test]
    fn test_c_tau_relationship() {
        // For McAllister 4-214-647:
        // g_s ≈ 0.00911, W₀ ≈ 2.3e-90, c_τ ≈ 3.34
        let g_s = 0.00911134;
        let w0 = 2.3e-90_f64;
        let w0_log_inv = (1.0 / w0).ln();

        let c_tau = compute_c_tau(g_s, w0_log_inv);
        assert!((c_tau - 3.34).abs() < 0.1);
    }
}
