//! Quintessence model helpers for reproducing published string examples.
//!
//! This module implements the explicit formulas used in
//! "From Inflation to Quintessence: a History of the Universe in String Theory"
//! (Cicoli et al., arXiv:2407.03405), Section 4.5 and eqs. (CYvol), (sols),
//! (LamI), (DecayConstants), (Cond3), (Cond4).
//!
//! We keep this numerics-only and focused on reproducing the paper's
//! published parameter point without shortcuts.

use std::f64::consts::PI;

use crate::f64_pos;
use crate::types::f64::F64;
use crate::types::i64::I64;
use crate::types::physics::{
    DivisorVolume, SmallCycleModulus, StringCoupling, Superpotential, Volume, XiCorrection,
};
use crate::types::tags::{Finite, Pos};

const TWO_PI: F64<Pos> = f64_pos!(2.0 * PI);
const SQRT_TWO: F64<Pos> = f64_pos!(std::f64::consts::SQRT_2);

/// Input parameters for the Cicoli 2407.03405 numerical example.
#[derive(Debug, Clone)]
pub struct Cicoli2407Params {
    /// k in the fibred volume form.
    pub k: F64<Pos>,
    /// \hat{k} in the small-cycle term.
    pub hat_k: F64<Pos>,
    /// String coupling g_s.
    pub g_s: StringCoupling,
    /// \xi in the LVS \alpha'^3 correction.
    pub xi: XiCorrection,
    /// |W_0|.
    pub w0: Superpotential,
    /// A_1 in non-perturbative term.
    pub a1_prefactor: F64<Pos>,
    /// A_2 in non-perturbative term.
    pub a2_prefactor: F64<Pos>,
    /// N_1 (instantons/gaugino condensate).
    pub n1: I64<Pos>,
    /// N_2 (instantons/gaugino condensate).
    pub n2: I64<Pos>,
}

/// Derived quantities for the Cicoli 2407.03405 example.
#[derive(Debug, Clone)]
pub struct Cicoli2407Derived {
    /// Stabilized small cycle volume <tau_s>.
    pub tau_s: SmallCycleModulus,
    /// Stabilized K3/T4 divisor volume <tau_1>.
    pub tau_1: DivisorVolume,
    /// Stabilized large divisor volume <tau_2>.
    pub tau_2: DivisorVolume,
    /// Calabi-Yau volume at the late-time minimum.
    pub volume: Volume,
    /// Axion decay constant f1.
    pub f1: F64<Pos>,
    /// Axion decay constant f2.
    pub f2: F64<Pos>,
    /// Lambda_1^4 scale in the late-time potential.
    pub lambda_1_4: F64<Pos>,
    /// Lambda_2^4 scale in the late-time potential.
    pub lambda_2_4: F64<Pos>,
}

fn pow_pos(x: F64<Pos>, exp: f64, label: &str) -> F64<Pos> {
    F64::<Pos>::new(x.get().powf(exp))
        .unwrap_or_else(|| panic!("powf produced non-positive value for {label}"))
}

/// Compute <tau_s> from eq. (sols):
/// <tau_s> = (hat_k^{1/3} * (3 xi)^{2/3}) / (2 g_s)
pub fn tau_s_from_lvs(hat_k: F64<Pos>, xi: XiCorrection, g_s: StringCoupling) -> SmallCycleModulus {
    let hat_k_third = pow_pos(hat_k, 1.0 / 3.0, "hat_k^(1/3)");
    let three_xi = f64_pos!(3.0) * xi;
    let three_xi_twothirds = pow_pos(three_xi, 2.0 / 3.0, "(3 xi)^(2/3)");
    let denominator = f64_pos!(2.0) * g_s;
    hat_k_third * three_xi_twothirds / denominator
}

/// Fibred CY volume form (eq. CYvol):
/// V = (1/sqrt(2k)) * sqrt(tau1) * tau2 - (1/3) * sqrt(2/hat_k) * tau_s^{3/2}
pub fn fibred_volume(
    k: F64<Pos>,
    hat_k: F64<Pos>,
    tau_1: DivisorVolume,
    tau_2: DivisorVolume,
    tau_s: SmallCycleModulus,
) -> Volume {
    let inv_sqrt_two_k = (f64_pos!(1.0) / (f64_pos!(2.0) * k)).sqrt();
    let first = inv_sqrt_two_k * tau_1.sqrt() * tau_2;

    let sqrt_two_over_hat_k = (f64_pos!(2.0) / hat_k).sqrt();
    let tau_s_3_2 = pow_pos(tau_s, 1.5, "tau_s^(3/2)");
    let second = f64_pos!(1.0 / 3.0) * sqrt_two_over_hat_k * tau_s_3_2;

    (first - second)
        .try_to_pos()
        .expect("fibred volume should be positive for valid parameters")
}

/// Decay constants (eq. DecayConstants).
pub fn decay_constants(
    n1: I64<Pos>,
    n2: I64<Pos>,
    tau_1: DivisorVolume,
    tau_2: DivisorVolume,
) -> (F64<Pos>, F64<Pos>) {
    let denom1 = f64_pos!(2.0) * SQRT_TWO * f64_pos!(PI) * tau_1;
    let denom2 = TWO_PI * tau_2;
    let f1 = n1.to_f64() / denom1;
    let f2 = n2.to_f64() / denom2;
    (f1, f2)
}

/// a_i = 2π / N_i for gaugino condensation or E3 instantons.
pub fn instanton_action(n_i: I64<Pos>) -> F64<Pos> {
    TWO_PI / n_i.to_f64()
}

/// Lambda_2^4 from eq. (LamI).
pub fn lambda2_4(
    w0: Superpotential,
    a2_prefactor: F64<Pos>,
    a2: F64<Pos>,
    tau_2: DivisorVolume,
    volume: Volume,
) -> F64<Pos> {
    let volume_sq = volume * volume;
    let exp_arg = -a2 * tau_2;
    let exp_term = F64::<Pos>::new(exp_arg.get().exp()).expect("exp is always positive");
    f64_pos!(4.0) * w0 * a2_prefactor * a2 * tau_2 * exp_term / volume_sq
}

/// Lambda_1^4 from eq. (LamI).
pub fn lambda1_4(
    lambda2_4: F64<Pos>,
    a1_prefactor: F64<Pos>,
    a1: F64<Pos>,
    tau_1: DivisorVolume,
    a2: F64<Pos>,
    tau_2: DivisorVolume,
) -> F64<Pos> {
    let ratio = (a1 * tau_1) / (a2 * tau_2);
    let factor = f64_pos!(1.0) + ratio;
    let exp_arg = -a1 * tau_1;
    let exp_term = F64::<Pos>::new(exp_arg.get().exp()).expect("exp is always positive");
    lambda2_4 * factor * a1_prefactor * exp_term
}

/// Axion masses at the late-time minimum (approximate):
/// m_i ≃ Lambda_i^2 / f_i.
pub fn axion_masses(
    lambda_1_4: F64<Pos>,
    lambda_2_4: F64<Pos>,
    f1: F64<Pos>,
    f2: F64<Pos>,
) -> (F64<Pos>, F64<Pos>) {
    let lambda_1_2 = lambda_1_4.sqrt();
    let lambda_2_2 = lambda_2_4.sqrt();
    let m1 = lambda_1_2 / f1;
    let m2 = lambda_2_2 / f2;
    (m1, m2)
}

/// Single-field axion potential: V = Lambda^4 (1 - cos(phi / f)).
#[derive(Debug, Clone)]
pub struct AxionPotential1D {
    /// Lambda^4 scale for the axion potential.
    pub lambda_4: F64<Pos>,
    /// Axion decay constant f.
    pub f: F64<Pos>,
}

impl AxionPotential1D {
    /// Typed potential value V(phi) in Planck units.
    #[must_use]
    pub fn value_typed(&self, phi: F64<Finite>) -> F64<Finite> {
        let phase = phi / self.f;
        let raw = 1.0 - phase.get().cos();
        let factor = F64::<Finite>::new(raw).expect("cos is finite");
        self.lambda_4 * factor
    }

    /// Typed derivative dV/dphi in Planck units.
    #[must_use]
    pub fn deriv_typed(&self, phi: F64<Finite>) -> F64<Finite> {
        let phase = phi / self.f;
        let raw = phase.get().sin();
        let factor = F64::<Finite>::new(raw).expect("sin is finite");
        (self.lambda_4 / self.f) * factor
    }
}

impl crate::cosmology::Potential for AxionPotential1D {
    fn value(&self, phi: f64) -> f64 {
        self.value_typed(F64::<Finite>::new(phi).expect("phi must be finite"))
            .get()
    }

    fn deriv(&self, phi: f64) -> f64 {
        self.deriv_typed(F64::<Finite>::new(phi).expect("phi must be finite"))
            .get()
    }
}

/// Late-time two-axion potential:
/// V = Lambda2^4 [1 - cos(phi2 / f2)] + Lambda1^4 [1 - cos(phi1 / f1 + phi2 / f2)].
pub fn late_time_potential_2d(
    phi_1: F64<Finite>,
    phi_2: F64<Finite>,
    f1: F64<Pos>,
    f2: F64<Pos>,
    lambda_1_4: F64<Pos>,
    lambda_2_4: F64<Pos>,
) -> F64<Finite> {
    let arg_2 = (phi_2 / f2).get();
    let arg_12 = (phi_1 / f1 + phi_2 / f2).get();
    let term2 = F64::<Finite>::new(1.0 - arg_2.cos()).expect("cos is finite");
    let term1 = F64::<Finite>::new(1.0 - arg_12.cos()).expect("cos is finite");
    lambda_2_4 * term2 + lambda_1_4 * term1
}

/// Gradient of the late-time two-axion potential (dV/dphi1, dV/dphi2).
pub fn late_time_gradient_2d(
    phi_1: F64<Finite>,
    phi_2: F64<Finite>,
    f1: F64<Pos>,
    f2: F64<Pos>,
    lambda_1_4: F64<Pos>,
    lambda_2_4: F64<Pos>,
) -> (F64<Finite>, F64<Finite>) {
    let arg_2 = (phi_2 / f2).get();
    let arg_12 = (phi_1 / f1 + phi_2 / f2).get();
    let sin_2 = F64::<Finite>::new(arg_2.sin()).expect("sin is finite");
    let sin_12 = F64::<Finite>::new(arg_12.sin()).expect("sin is finite");
    let dphi1 = (lambda_1_4 / f1) * sin_12;
    let dphi2 = (lambda_2_4 / f2) * sin_2 + (lambda_1_4 / f2) * sin_12;
    (dphi1, dphi2)
}

/// Derive the Cicoli 2407.03405 numerical example from inputs.
pub fn derive_cicoli_2407(
    params: &Cicoli2407Params,
    f1_target: F64<Pos>,
    f2_target: F64<Pos>,
) -> Cicoli2407Derived {
    // From Cond3 and Cond4 (rearranged):
    let tau_1 = params.n1.to_f64() / (f64_pos!(2.0) * SQRT_TWO * f64_pos!(PI) * f1_target);
    let tau_2 = params.n2.to_f64() / (TWO_PI * f2_target);

    let tau_s = tau_s_from_lvs(params.hat_k, params.xi, params.g_s);
    let volume = fibred_volume(params.k, params.hat_k, tau_1, tau_2, tau_s);

    let (f1, f2) = decay_constants(params.n1, params.n2, tau_1, tau_2);

    let a1 = instanton_action(params.n1);
    let a2 = instanton_action(params.n2);
    let lambda_2_4 = lambda2_4(params.w0, params.a2_prefactor, a2, tau_2, volume);
    let lambda_1_4 = lambda1_4(lambda_2_4, params.a1_prefactor, a1, tau_1, a2, tau_2);

    Cicoli2407Derived {
        tau_s,
        tau_1,
        tau_2,
        volume,
        f1,
        f2,
        lambda_1_4,
        lambda_2_4,
    }
}
