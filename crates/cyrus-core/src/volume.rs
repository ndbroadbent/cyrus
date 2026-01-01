//! Calabi-Yau volume computations.
//!
//! The string frame volume is:
//! ```text
//! V_string = (1/6) κ_ijk t^i t^j t^k - BBHL
//! ```
//!
//! Where:
//! - `κ_ijk`: intersection numbers
//! - `t^i`: Kähler moduli
//! - BBHL = ζ(3) χ(X) / (4(2π)³) is the α' correction
//! - χ(X) = 2(h¹¹ - h²¹) is the Euler characteristic
//!
//! Reference: arXiv:2107.09064

use crate::intersection::Intersection;

/// Riemann zeta function at 3: ζ(3) ≈ 1.202056903159594.
const ZETA_3: f64 = 1.202_056_903_159_594;

/// Compute the BBHL α' correction.
///
/// ```text
/// BBHL = ζ(3) χ(X) / (4(2π)³)
/// ```
///
/// where χ = 2(h11 - h21) is the Euler characteristic.
pub fn bbhl_correction(h11: i32, h21: i32) -> f64 {
    let chi = 2 * (h11 - h21);
    let two_pi_cubed = (2.0 * std::f64::consts::PI).powi(3);
    ZETA_3 * f64::from(chi) / (4.0 * two_pi_cubed)
}

/// Compute classical CY volume: V = (1/6) `κ_ijk` t^i t^j t^k.
pub fn volume_classical(kappa: &Intersection, t: &[f64]) -> f64 {
    kappa.contract_triple(t) / 6.0
}

/// Compute string frame volume with BBHL correction.
///
/// ```text
/// V_string = (1/6) κ_ijk t^i t^j t^k - BBHL
/// ```
pub fn volume_string(kappa: &Intersection, t: &[f64], h11: i32, h21: i32) -> f64 {
    volume_classical(kappa, t) - bbhl_correction(h11, h21)
}

/// Result of volume computation with breakdown.
#[derive(Debug, Clone)]
pub struct VolumeResult {
    /// Classical volume (1/6) κ t³.
    pub classical: f64,
    /// BBHL α' correction.
    pub bbhl: f64,
    /// String frame volume (classical - BBHL).
    pub string_frame: f64,
}

/// Compute volume with detailed breakdown.
pub fn compute_volume(kappa: &Intersection, t: &[f64], h11: i32, h21: i32) -> VolumeResult {
    let classical = volume_classical(kappa, t);
    let bbhl = bbhl_correction(h11, h21);
    VolumeResult {
        classical,
        bbhl,
        string_frame: classical - bbhl,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bbhl_correction() {
        // For h11=214, h21=4: χ = 2(214-4) = 420
        // BBHL = ζ(3) * 420 / (4 * (2π)³)
        //      = 1.202... * 420 / (4 * 248.05...)
        //      ≈ 0.509
        let bbhl = bbhl_correction(214, 4);
        assert!((bbhl - 0.509).abs() < 0.001);
    }

    #[test]
    fn test_volume_classical() {
        // Simple case: κ_000 = 6, t = [2]
        // V = (1/6) * 6 * 8 = 8
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, 6);

        let v = volume_classical(&kappa, &[2.0]);
        assert!((v - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_volume_string() {
        // V_string = V_classical - BBHL
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, 6);

        let v_string = volume_string(&kappa, &[2.0], 5, 3);
        let v_classical = 8.0;
        let bbhl = bbhl_correction(5, 3);

        assert!((v_string - (v_classical - bbhl)).abs() < 1e-10);
    }
}
