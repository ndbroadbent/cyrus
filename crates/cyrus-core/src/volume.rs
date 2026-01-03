//! Calabi-Yau volume computations.
//!
//! The string frame volume is:
//! ```text
//! V_string = (1/6) κ_ijk t^i t^j t^k - BBHL
//! ```
//!
//! Where:
//! - `κ_ijk`: intersection numbers (positive)
//! - `t^i`: Kähler moduli (positive)
//! - BBHL = ζ(3) χ(X) / (4(2π)³) is the α' correction
//! - χ(X) = 2(h¹¹ - h²¹) is the Euler characteristic
//!
//! Volume is always positive for valid physics - this is enforced by the type system.
//!
//! Reference: arXiv:2107.09064

use crate::f64_pos;
use crate::intersection::{Intersection, NonEmptyIntersection};
use crate::policy::{Strict, VolumePolicy};
use crate::types::branded::Moduli;
use crate::types::f64::F64;
use crate::types::i32::I32;
use crate::types::physics::{BBHLCorrection, H11, H21, Volume};
use crate::types::tags::{Finite, Pos, Two};

/// Riemann zeta function at 3: ζ(3) ≈ 1.202056903159594 (positive).
const ZETA_3: F64<Pos> = f64_pos!(1.202_056_903_159_594);

/// 4(2π)³ - the denominator in BBHL correction (positive).
/// 4 * (2π)³ = 4 * 8π³ = 32π³ ≈ 992.17
const FOUR_TWO_PI_CUBED: F64<Pos> = f64_pos!(32.0 * 31.006_276_680_299_816); // π³ ≈ 31.006...

/// Compute the BBHL α' correction.
///
/// ```text
/// BBHL = ζ(3) χ(X) / (4(2π)³)
/// ```
///
/// where χ = 2(h11 - h21) is the Euler characteristic.
///
/// # Arguments
/// - `h11`: Hodge number h¹¹ (≥ 1 for CY3)
/// - `h21`: Hodge number h²¹ (≥ 0)
///
/// # Returns
/// BBHL correction (can be positive, negative, or zero depending on χ)
pub fn bbhl_correction(h11: H11, h21: H21) -> BBHLCorrection {
    // χ = 2(h11 - h21) can be positive, negative, or zero
    let diff = h11.get() - h21.get(); // Finite (can be pos, neg, or zero)
    let chi = I32::<Two>::TWO.get() * diff;

    // BBHL = ζ(3) × χ / (4(2π)³)
    // ζ(3) > 0, 4(2π)³ > 0, so sign of BBHL = sign of χ
    let result = (ZETA_3.get() * f64::from(chi)) / FOUR_TWO_PI_CUBED.get();
    F64::<Finite>::new(result).expect("BBHL is always finite")
}

/// Compute classical CY volume: V = (1/6) `κ_ijk` t^i t^j t^k.
///
/// Pure function: positive inputs, positive output.
#[must_use]
pub fn volume_classical(kappa: &NonEmptyIntersection, t: &Moduli<'_>) -> F64<Pos> {
    kappa.contract_triple(t) / f64_pos!(6.0)
}

/// Compute string frame volume (internal, returns Finite).
///
/// This is the raw computation that can produce negative values.
/// Use `volume_string` or `volume_string_with_policy` for the public API.
#[inline(always)]
fn volume_string_raw(
    kappa: &NonEmptyIntersection,
    t: &Moduli<'_>,
    h11: H11,
    h21: H21,
) -> F64<Finite> {
    let classical = volume_classical(kappa, t);
    let bbhl = bbhl_correction(h11, h21);
    // Pos - Finite = Finite (safe unwrap)
    F64::<Finite>::new(classical.get() - bbhl.get()).unwrap()
}

/// Compute string frame volume with policy-determined handling.
///
/// ```text
/// V_string = (1/6) κ_ijk t^i t^j t^k - BBHL
/// ```
///
/// # Policies
///
/// - `Strict`: Returns `Option<F64<Pos>>` - None if non-positive
/// - `ForGA`: Returns `F64<Finite>` - keeps negative values for fitness
/// - `Abort`: Returns `F64<Pos>` - panics if non-positive
///
/// # Example
///
/// ```rust,ignore
/// use cyrus_core::{Strict, ForGA, volume_string_with_policy};
///
/// // Strict: None if invalid
/// let v: Option<F64<Pos>> = volume_string_with_policy::<Strict>(&kappa, &t, h11, h21);
///
/// // GA: keeps negative for fitness
/// let v: F64<Finite> = volume_string_with_policy::<ForGA>(&kappa, &t, h11, h21);
/// ```
#[must_use]
#[inline(always)]
pub fn volume_string_with_policy<P: VolumePolicy>(
    kappa: &NonEmptyIntersection,
    t: &Moduli<'_>,
    h11: H11,
    h21: H21,
) -> P::Output {
    P::apply(volume_string_raw(kappa, t, h11, h21))
}

/// Compute string frame volume with BBHL correction (strict mode).
///
/// ```text
/// V_string = (1/6) κ_ijk t^i t^j t^k - BBHL
/// ```
///
/// Returns `None` if string frame volume is non-positive (invalid physics).
///
/// This is equivalent to `volume_string_with_policy::<Strict>(...)`.
#[must_use]
pub fn volume_string(
    kappa: &NonEmptyIntersection,
    t: &Moduli<'_>,
    h11: H11,
    h21: H21,
) -> Option<Volume> {
    volume_string_with_policy::<Strict>(kappa, t, h11, h21)
}

/// Result of volume computation with breakdown.
#[derive(Debug, Clone)]
pub struct VolumeResult {
    /// Classical volume (1/6) κ t³.
    pub classical: F64<Pos>,
    /// BBHL α' correction (can be positive, negative, or zero).
    pub bbhl: F64<Finite>,
    /// String frame volume (classical - BBHL).
    pub string_frame: F64<Pos>,
}

/// Compute volume with detailed breakdown.
///
/// Returns `None` if string frame volume is non-positive (invalid physics).
#[must_use]
pub fn compute_volume(
    kappa: &NonEmptyIntersection,
    t: &Moduli<'_>,
    h11: H11,
    h21: H21,
) -> Option<VolumeResult> {
    let classical = volume_classical(kappa, t);
    let bbhl = bbhl_correction(h11, h21);
    let string_frame = F64::<Pos>::new(classical.get() - bbhl.get())?;

    Some(VolumeResult {
        classical,
        bbhl,
        string_frame,
    })
}

/// Raw volume computation result for pipeline use.
#[derive(Debug, Clone)]
pub struct VolumeResultRaw {
    /// Classical volume (1/6) κ t³ (positive).
    pub classical: F64<Pos>,
    /// BBHL α' correction (can be positive, negative, or zero).
    pub bbhl: F64<Finite>,
    /// String frame volume (classical - BBHL).
    pub string_frame: F64<Finite>,
}

/// Compute volume with raw types (for pipeline compatibility).
///
/// Uses raw `&Intersection` and `&[f64]` for cases where moduli
/// come from flat direction solving rather than GA generation.
///
/// Returns `None` if classical volume is non-positive.
#[must_use]
pub fn compute_volume_raw(
    kappa: &Intersection,
    t: &[f64],
    h11: H11,
    h21: H21,
) -> Option<VolumeResultRaw> {
    // Convert raw values to typed (flat direction moduli can have any sign)
    let t_typed: Vec<F64<Finite>> = t
        .iter()
        .map(|&x| F64::<Finite>::new(x).expect("moduli must be finite"))
        .collect();

    // Contraction and division by 6
    let contraction = kappa.contract_triple_finite(&t_typed)?;
    let classical = contraction / f64_pos!(6.0);

    let bbhl = bbhl_correction(h11, h21);
    // Type algebra: Pos - Finite = Finite
    let string_frame = classical - bbhl;

    Some(VolumeResultRaw {
        classical,
        bbhl,
        string_frame,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::intersection::Intersection;
    use crate::types::rational::Rational as TypedRational;
    use crate::types::tags::{GTEOne, NonNeg, Pos};
    use malachite::Rational;

    #[test]
    fn test_bbhl_correction() {
        // For h11=214, h21=4: χ = 2(214-4) = 420
        // BBHL = ζ(3) * 420 / (4 * (2π)³)
        //      = 1.202... * 420 / (4 * 248.05...)
        //      ≈ 0.509
        let h11 = I32::<GTEOne>::new(214).unwrap();
        let h21 = I32::<NonNeg>::new(4).unwrap();
        let bbhl = bbhl_correction(h11, h21);
        assert!((bbhl.get() - 0.509).abs() < 0.001);
    }

    #[test]
    fn test_volume_classical() {
        // Simple case: κ_000 = 6, t = [2]
        // V = (1/6) * 6 * 8 = 8
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, TypedRational::<Pos>::new(Rational::from(6)).unwrap());

        let kappa = kappa.into_non_empty().unwrap();
        let t = kappa.dim().validate(&[2.0]).unwrap();

        let v = volume_classical(&kappa, &t);
        assert!((v.get() - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_volume_classical_negative_moduli_rejected() {
        // Negative moduli are rejected at validation
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, TypedRational::<Pos>::new(Rational::from(6)).unwrap());

        let kappa = kappa.into_non_empty().unwrap();
        let result = kappa.dim().validate(&[-2.0]);
        assert!(result.is_none());
    }

    #[test]
    fn test_volume_string() {
        // V_string = V_classical - BBHL
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, TypedRational::<Pos>::new(Rational::from(6)).unwrap());

        let kappa = kappa.into_non_empty().unwrap();
        let t = kappa.dim().validate(&[2.0]).unwrap();

        let h11 = I32::<GTEOne>::new(5).unwrap();
        let h21 = I32::<NonNeg>::new(3).unwrap();
        let v_string = volume_string(&kappa, &t, h11, h21).unwrap();
        let v_classical = 8.0;
        let bbhl = bbhl_correction(h11, h21);

        assert!((v_string.get() - (v_classical - bbhl.get())).abs() < 1e-10);
    }

    #[test]
    fn test_compute_volume() {
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, TypedRational::<Pos>::new(Rational::from(6)).unwrap());

        let kappa = kappa.into_non_empty().unwrap();
        let t = kappa.dim().validate(&[2.0]).unwrap();

        let h11 = I32::<GTEOne>::new(5).unwrap();
        let h21 = I32::<NonNeg>::new(3).unwrap();
        let result = compute_volume(&kappa, &t, h11, h21).unwrap();
        assert!((result.classical.get() - 8.0).abs() < 1e-10);
        assert!((result.string_frame.get() - (8.0 - result.bbhl.get())).abs() < 1e-10);
    }
}
