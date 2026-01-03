//! Vacuum energy (cosmological constant) computation.
//!
//! The KKLT vacuum energy formula is:
//! ```text
//! V₀ = -3 × e^{K₀} × (g_s⁷ / (4V_string)²) × W₀²
//! ```
//!
//! Where:
//! - `e^{K₀}`: Kähler potential factor from flat direction
//! - `g_s`: String coupling
//! - `V_string`: String frame volume with BBHL correction
//! - `W₀`: Flux superpotential
//!
//! Since all inputs are positive, the output is **guaranteed negative**
//! by the type system.
//!
//! Reference: arXiv:2107.09064, eq. 6.63

use crate::types::f64::F64;
use crate::types::tags::{Neg, Pos};

/// Compute the flux tadpole contribution: `Q_flux = -1/2 * K · M`.
///
/// Reference: arXiv:2107.09064, eq. 6.13
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn compute_tadpole(k: &[i64], m: &[i64]) -> f64 {
    let dot: i64 = k.iter().zip(m.iter()).map(|(&ki, &mi)| ki * mi).sum();
    -0.5 * dot as f64
}

/// Check if the flux tadpole is canceled given a bound `q_max`.
///
/// `Q_flux <= Q_max`
#[must_use]
pub fn is_tadpole_canceled(k: &[i64], m: &[i64], q_max: f64) -> bool {
    compute_tadpole(k, m) <= q_max
}

/// Compute vacuum energy V₀ = -3 × `e^{K₀}` × (`g_s`⁷ / (4V)²) × W₀².
///
/// Returns the cosmological constant in Planck units (Mpl⁴).
///
/// # Type Safety
/// All inputs must be positive. The return type `F64<Neg>` guarantees
/// the result is always negative - no runtime assertion needed.
///
/// # Arguments
/// * `ek0` - Kähler potential factor `e^{K₀}` (positive)
/// * `g_s` - String coupling (positive)
/// * `v_string` - String frame volume with BBHL correction (positive)
/// * `w0` - Flux superpotential magnitude |W₀| (positive)
#[must_use]
pub fn compute_v0(ek0: F64<Pos>, g_s: F64<Pos>, v_string: F64<Pos>, w0: F64<Pos>) -> F64<Neg> {
    // All inputs positive → product is positive → multiply by -3 → negative
    let value =
        -3.0 * ek0.get() * (g_s.get().powi(7) / (4.0 * v_string.get()).powi(2)) * w0.get().powi(2);
    // SAFETY: formula guarantees negative result when all inputs are positive
    F64::<Neg>::new(value).expect("V₀ formula with positive inputs always yields negative result")
}

/// Result of vacuum energy computation with breakdown.
#[derive(Debug, Clone)]
pub struct VacuumResult {
    /// Kähler potential factor `e^{K₀}`.
    pub ek0: F64<Pos>,
    /// String coupling `g_s`.
    pub g_s: F64<Pos>,
    /// String frame volume.
    pub v_string: F64<Pos>,
    /// Flux superpotential |W₀|.
    pub w0: F64<Pos>,
    /// Vacuum energy V₀ (cosmological constant) - always negative.
    pub v0: F64<Neg>,
}

/// Compute vacuum energy with detailed breakdown.
#[must_use]
pub fn compute_vacuum(
    ek0: F64<Pos>,
    g_s: F64<Pos>,
    v_string: F64<Pos>,
    w0: F64<Pos>,
) -> VacuumResult {
    let v0 = compute_v0(ek0, g_s, v_string, w0);
    VacuumResult {
        ek0,
        g_s,
        v_string,
        w0,
        v0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper to create F64<Pos> in tests
    fn pos(v: f64) -> F64<Pos> {
        F64::<Pos>::new(v).unwrap()
    }

    #[test]
    fn test_compute_v0_mcallister_4_214_647() {
        // From McAllister arXiv:2107.09064 example 4-214-647
        // e^K₀ ≈ 0.234393
        // g_s ≈ 0.00911134
        // V_string ≈ 4711.83
        // W₀ ≈ 2.3e-90
        // Expected V₀ ≈ -5.5e-203

        let ek0 = pos(0.234_393);
        let g_s = pos(0.009_111_34);
        let v_string = pos(4711.83);
        let w0 = pos(2.3e-90);

        let v0 = compute_v0(ek0, g_s, v_string, w0);

        // No need to assert v0 < 0 - the type F64<Neg> guarantees it!

        // Check order of magnitude: should be around -5.5e-203
        // Due to floating point limits with such small numbers, we check the exponent
        let log_v0 = v0.get().abs().log10();
        assert!(
            (-204.0..=-202.0).contains(&log_v0),
            "log10(|V₀|) = {log_v0}"
        );
    }

    #[test]
    fn test_compute_v0_simple() {
        // Simple case with unit values
        // V₀ = -3 × 1 × (1⁷ / (4×1)²) × 1² = -3/16 = -0.1875
        let v0 = compute_v0(pos(1.0), pos(1.0), pos(1.0), pos(1.0));
        assert!((v0.get() - (-0.1875)).abs() < 1e-10);
    }

    #[test]
    fn test_vacuum_result() {
        let result = compute_vacuum(pos(0.25), pos(0.01), pos(1000.0), pos(1e-50));
        assert!((result.ek0.get() - 0.25).abs() < 1e-14);
        assert!((result.g_s.get() - 0.01).abs() < 1e-14);
        assert!((result.v_string.get() - 1000.0).abs() < 1e-10);
        assert!((result.w0.get() - 1e-50).abs() < 1e-60);
        // No need to assert v0 < 0 - type guarantees it!
        assert!(result.v0.get() < 0.0); // We can still check if we want
    }

    #[test]
    fn test_compute_tadpole() {
        // K=[-3, -5], M=[10, 11]
        // K.M = -30 - 55 = -85
        // Q = -0.5 * -85 = 42.5
        let k = vec![-3, -5];
        let m = vec![10, 11];
        let q = compute_tadpole(&k, &m);
        assert!((q - 42.5).abs() < 1e-10);
    }

    #[test]
    fn test_is_tadpole_canceled() {
        let k = vec![1];
        let m = vec![-2]; // dot = -2 -> Q = 1.0
        assert!(is_tadpole_canceled(&k, &m, 1.0));
        assert!(!is_tadpole_canceled(&k, &m, 0.9));
    }
}
