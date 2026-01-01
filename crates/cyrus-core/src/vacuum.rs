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
//! Reference: arXiv:2107.09064, eq. 6.63

/// Compute vacuum energy V₀ = -3 × `e^{K₀}` × (`g_s`⁷ / (4V)²) × W₀².
///
/// Returns the cosmological constant in Planck units (Mpl⁴).
///
/// # Arguments
/// * `ek0` - Kähler potential factor `e^{K₀}`
/// * `g_s` - String coupling
/// * `v_string` - String frame volume with BBHL correction
/// * `w0` - Flux superpotential magnitude |W₀|
#[must_use]
pub fn compute_v0(ek0: f64, g_s: f64, v_string: f64, w0: f64) -> f64 {
    -3.0 * ek0 * (g_s.powi(7) / (4.0 * v_string).powi(2)) * w0.powi(2)
}

/// Result of vacuum energy computation with breakdown.
#[derive(Debug, Clone)]
pub struct VacuumResult {
    /// Kähler potential factor `e^{K₀}`.
    pub ek0: f64,
    /// String coupling `g_s`.
    pub g_s: f64,
    /// String frame volume.
    pub v_string: f64,
    /// Flux superpotential |W₀|.
    pub w0: f64,
    /// Vacuum energy V₀ (cosmological constant).
    pub v0: f64,
}

/// Compute vacuum energy with detailed breakdown.
#[must_use]
pub fn compute_vacuum(ek0: f64, g_s: f64, v_string: f64, w0: f64) -> VacuumResult {
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

    #[test]
    fn test_compute_v0_mcallister_4_214_647() {
        // From McAllister arXiv:2107.09064 example 4-214-647
        // e^K₀ ≈ 0.234393
        // g_s ≈ 0.00911134
        // V_string ≈ 4711.83
        // W₀ ≈ 2.3e-90
        // Expected V₀ ≈ -5.5e-203

        let ek0 = 0.234_393;
        let g_s = 0.009_111_34;
        let v_string = 4711.83;
        let w0 = 2.3e-90;

        let v0 = compute_v0(ek0, g_s, v_string, w0);

        // Check sign is negative
        assert!(v0 < 0.0);

        // Check order of magnitude: should be around -5.5e-203
        // Due to floating point limits with such small numbers, we check the exponent
        let log_v0 = v0.abs().log10();
        assert!(
            (-204.0..=-202.0).contains(&log_v0),
            "log10(|V₀|) = {log_v0}"
        );
    }

    #[test]
    fn test_compute_v0_simple() {
        // Simple case with unit values
        // V₀ = -3 × 1 × (1⁷ / (4×1)²) × 1² = -3/16 = -0.1875
        let v0 = compute_v0(1.0, 1.0, 1.0, 1.0);
        assert!((v0 - (-0.1875)).abs() < 1e-10);
    }

    #[test]
    fn test_vacuum_result() {
        let result = compute_vacuum(0.25, 0.01, 1000.0, 1e-50);
        assert!((result.ek0 - 0.25).abs() < 1e-14);
        assert!((result.g_s - 0.01).abs() < 1e-14);
        assert!((result.v_string - 1000.0).abs() < 1e-10);
        assert!((result.w0 - 1e-50).abs() < 1e-60);
        assert!(result.v0 < 0.0);
    }
}
