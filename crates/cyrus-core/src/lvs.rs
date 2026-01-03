//! Large Volume Scenario (LVS) potential computation.
//!
//! LVS stabilizes the volume of the Calabi-Yau at exponentially large values
//! using perturbative alpha' corrections and non-perturbative instanton effects.
//!
//! Reference: arXiv:2407.03405 (Cicoli et al. 2024)

use std::f64::consts::PI;

use crate::f64_pos;
use crate::types::f64::F64;
use crate::types::i32::I32;
use crate::types::physics::{
    LvsPotential, SmallCycleModulus, StringCoupling, Volume, XiCorrection,
};
use crate::types::tags::{Neg, Pos};

/// 2π as a typed positive constant.
const TWO_PI: F64<Pos> = f64_pos!(2.0 * PI);

/// ζ(3) ≈ 1.202056903159594 (positive).
const ZETA_3: F64<Pos> = f64_pos!(1.202_056_903_159_594);

/// Parameters for LVS potential.
#[derive(Debug, Clone)]
pub struct LvsParams {
    /// Flux superpotential magnitude |W0| (positive).
    pub w0: F64<Pos>,
    /// String coupling g_s (positive).
    pub g_s: StringCoupling,
    /// Euler characteristic of the CY (negative for LVS).
    pub chi: I32<Neg>,
    /// Non-perturbative prefactor A_s (positive).
    pub a_s: F64<Pos>,
    /// Non-perturbative exponent coefficient λ_s = 2π/N (positive).
    pub lambda_s: F64<Pos>,
}

/// Result of LVS potential evaluation.
#[derive(Debug, Clone)]
pub struct LvsResult {
    /// Volume of the CY in string units (positive).
    pub volume: Volume,
    /// Stabilized small cycle volume τ_s (positive).
    pub tau_s: SmallCycleModulus,
    /// Scalar potential V_LVS (can be any sign during optimization).
    pub potential: LvsPotential,
    /// Alpha' correction parameter ξ (positive for LVS).
    pub xi: XiCorrection,
}

/// Compute the scalar potential V(vol, τ_s).
///
/// V_LVS = (8/3) * (a_s² * √τ_s / vol) * exp(-2λ_s τ_s)
///       - 4 * (W0 * a_s * τ_s / vol²) * exp(-λ_s τ_s)
///       + 3 * ξ * W0² / (4 * g_s^(3/2) * vol³)
pub fn compute_lvs_potential(
    params: &LvsParams,
    vol: Volume,
    tau_s: SmallCycleModulus,
    xi: XiCorrection,
) -> LvsPotential {
    // All inputs are Pos, so we can compute safely

    // term1 = (8/3) * a_s² * √τ_s * exp(-2λ_s τ_s) / vol
    // Pos * Pos = Pos, Pos / Pos = Pos
    let a_s_sq = params.a_s * params.a_s;
    let sqrt_tau = F64::<Pos>::new(tau_s.get().sqrt()).expect("sqrt of Pos is Pos");
    let exp_arg1 = f64_pos!(2.0) * params.lambda_s * tau_s;
    let exp1 = F64::<Pos>::new((-exp_arg1.get()).exp()).expect("exp is always positive");
    let term1 = f64_pos!(8.0 / 3.0) * a_s_sq * sqrt_tau * exp1 / vol;

    // term2 = -4 * W0 * a_s * τ_s * exp(-λ_s τ_s) / vol²
    // This is negative (the minus sign)
    let exp_arg2 = params.lambda_s * tau_s;
    let exp2 = F64::<Pos>::new((-exp_arg2.get()).exp()).expect("exp is always positive");
    let vol_sq = vol * vol;
    let term2_magnitude = f64_pos!(4.0) * params.w0 * params.a_s * tau_s * exp2 / vol_sq;
    // Negate: -Pos = Neg
    let term2 = -term2_magnitude;

    // term3 = 3 * ξ * W0² / (4 * g_s^(3/2) * vol³)
    // All positive
    let w0_sq = params.w0 * params.w0;
    let g_s_32 = F64::<Pos>::new(params.g_s.get().powf(1.5)).expect("Pos^1.5 is Pos");
    let vol_cubed = vol * vol * vol;
    let term3 = f64_pos!(3.0) * xi * w0_sq / (f64_pos!(4.0) * g_s_32 * vol_cubed);

    // Sum: Pos + Neg + Pos = Finite (algebra handles it)
    term1 + term2 + term3
}

/// Compute the LVS potential minimum.
///
/// Returns the stabilized volume and potential energy.
pub fn compute_lvs_vacuum(params: &LvsParams) -> Option<LvsResult> {
    // ξ = -χ × ζ(3) / (2 × (2π)³)
    // χ is Neg, so -χ is Pos, making ξ positive
    let neg_chi = -params.chi; // Neg → Pos after negation... actually I32<Neg> negated is I32<Pos>
    let chi_f64 = F64::<Pos>::new(f64::from(neg_chi.get())).expect("negated Neg chi is positive");

    // (2π)³ = 8π³
    let two_pi_cubed = TWO_PI * TWO_PI * TWO_PI;
    let denominator = f64_pos!(2.0) * two_pi_cubed;

    // ξ = chi_f64 × ζ(3) / denominator - all Pos
    let xi = chi_f64 * ZETA_3 / denominator;

    // Initial guess for τ_s
    let tau_s_init = f64_pos!(2.0) / params.lambda_s;

    // Initial guess for volume
    let exp_lambda_tau = F64::<Pos>::new((params.lambda_s * tau_s_init).get().exp())
        .expect("exp is always positive");
    let vol_init = params.w0 * exp_lambda_tau;

    let mut tau_s = tau_s_init;
    let mut vol = vol_init;

    // Gradient descent loop
    let lr = f64_pos!(0.01);
    let eps = f64_pos!(1e-5);

    for _ in 0..100 {
        let v_curr = compute_lvs_potential(params, vol, tau_s, xi);

        // Numerical gradients
        let vol_plus = vol + eps; // Pos + Pos = Pos
        let v_vol = compute_lvs_potential(params, vol_plus, tau_s, xi);
        let grad_vol = (v_vol - v_curr) / eps; // Finite / Pos = Finite

        let tau_plus = tau_s + eps; // Pos + Pos = Pos
        let v_tau = compute_lvs_potential(params, vol, tau_plus, xi);
        let grad_tau = (v_tau - v_curr) / eps; // Finite / Pos = Finite

        // Update - if values go non-positive, optimization diverged
        let vol_update = vol - lr * f64_pos!(1000.0) * grad_vol;
        vol = vol_update.try_to_pos()?;

        let tau_update = tau_s - lr * grad_tau;
        tau_s = tau_update.try_to_pos()?;

        // Check convergence (compare raw values since abs() returns NonNeg)
        if grad_vol.abs().get() < 1e-6 && grad_tau.abs().get() < 1e-6 {
            break;
        }
    }

    Some(LvsResult {
        volume: vol,
        tau_s,
        potential: compute_lvs_potential(params, vol, tau_s, xi),
        xi,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lvs_potential_minimum() {
        let params = LvsParams {
            w0: f64_pos!(1.0),
            g_s: f64_pos!(0.1),
            chi: I32::<Neg>::new(-100).unwrap(),
            a_s: f64_pos!(1.0),
            lambda_s: TWO_PI,
        };

        let result = compute_lvs_vacuum(&params).unwrap();
        assert!(result.volume.get() > 10.0);
    }

    #[test]
    fn test_lvs_result_fields() {
        let params = LvsParams {
            w0: f64_pos!(10.0),
            g_s: f64_pos!(0.05),
            chi: I32::<Neg>::new(-200).unwrap(),
            a_s: f64_pos!(1.0),
            lambda_s: TWO_PI,
        };

        let result = compute_lvs_vacuum(&params).unwrap();
        assert!(result.xi.get() > 0.0);
        assert!(result.tau_s.get() > 0.0);
        assert!(result.volume.get().is_finite());
        assert!(result.potential.get().is_finite());
    }
}
