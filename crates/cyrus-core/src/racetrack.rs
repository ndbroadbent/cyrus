//! Solver for the racetrack Kähler moduli stabilization.
//!
//! Given GV invariants and flat direction p, builds the racetrack potential
//! and solves for the stabilized Kähler moduli.
//!
//! Reference: arXiv:2107.09064, Section 6.4

use std::collections::HashMap;
use std::f64::consts::PI;

use crate::f64_pos;
use crate::types::f64::F64;
use crate::types::i64::I64;
use crate::types::physics::{
    GvValue, ImTau, RacetrackCoefficient, RacetrackExponent, ReTau, ResidualError, StringCoupling,
    Superpotential,
};
use crate::types::tags::{Finite, NonNeg, Pos};

/// Riemann zeta function at 3: ζ(3) ≈ 1.202056903159594 (positive).
pub const ZETA_3: F64<Pos> = f64_pos!(1.202_056_903_159_594);

/// Superpotential normalization from eq. 2.22: ζ = 1/(2^{3/2} π^{5/2}).
///
/// This factor multiplies the flux superpotential sum.
/// Value: 1 / (2.828... × 17.49...) ≈ 0.02024
const ZETA_NORM: F64<Pos> = f64_pos!(1.0 / (2.828_427_124_746_19 * 17.493_418_327_624_86));

/// 2π as a typed positive constant.
const TWO_PI: F64<Pos> = f64_pos!(2.0 * PI);
const DILOG_TOL: f64 = 1e-16;
const DILOG_MAX_TERMS: usize = 100_000;

/// A racetrack term: `coefficient * exp(-2π q·p / g_s)`.
#[derive(Debug, Clone)]
pub struct RacetrackTerm {
    /// Curve indices (q) - can be any integers.
    pub curve: Vec<I64<Finite>>,
    /// Pre-exponential factor: `(M·q) * N_q` (can be any sign).
    pub coefficient: RacetrackCoefficient,
    /// Effective exponent: `q·p` (can be any sign).
    pub exponent: RacetrackExponent,
}

/// Result of racetrack stabilization.
#[derive(Debug, Clone)]
pub struct RacetrackResult {
    /// String coupling g_s (positive, bounded 0 < g_s < 1).
    pub g_s: StringCoupling,
    /// Real part of the stabilized axio-dilaton.
    pub re_tau: ReTau,
    /// Stabilized Kähler moduli value: `Im(τ) = (p / g_s)` (positive).
    pub im_tau: ImTau,
    /// Numerical error delta (F-term residual, non-negative).
    pub delta: ResidualError,
    /// Numerical error epsilon (distance from target, non-negative).
    pub epsilon: ResidualError,
}

/// A Gopakumar-Vafa invariant.
#[derive(Debug, Clone)]
pub struct GvInvariant {
    /// Curve class vector - integers that can be positive, negative, or zero.
    pub curve: Vec<I64<Finite>>,
    /// GV invariant value N_q (can be any finite value).
    pub value: GvValue,
}

/// Build racetrack terms from GV invariants.
///
/// Group curves with identical action (q·p) and sum their coefficients.
pub fn build_racetrack_terms(
    gv_invariants: &[GvInvariant],
    m: &[I64<Finite>],
    p: &[F64<Finite>],
    cutoff: F64<Pos>,
) -> Vec<RacetrackTerm> {
    let mut raw_terms = Vec::new();

    for gv in gv_invariants {
        // q·p (curve dot flat direction) - Finite * Finite = Finite, sum of Finite = Finite
        let q_dot_p: F64<Finite> = gv
            .curve
            .iter()
            .zip(p.iter())
            .map(|(qi, pi)| qi.to_f64() * *pi)
            .fold(F64::<Finite>::ZERO, |acc, x| acc + x);

        // Match McAllister: keep only 0 < q·p < cutoff.
        if q_dot_p.get() <= 0.0 || q_dot_p.get() >= cutoff.get() {
            continue;
        }

        // M·q (flux dot curve) - Finite * Finite = Finite, sum of Finite = Finite
        let m_dot_q: I64<Finite> = gv
            .curve
            .iter()
            .zip(m.iter())
            .map(|(qi, mi)| *qi * *mi)
            .fold(I64::<Finite>::ZERO, |acc, x| acc + x);

        // coefficient = (M·q) × N_q - Finite * Finite = Finite
        let coefficient: RacetrackCoefficient = m_dot_q.to_f64() * gv.value;

        raw_terms.push(RacetrackTerm {
            curve: gv.curve.clone(),
            coefficient,
            exponent: q_dot_p,
        });
    }

    // Group terms by exponent (q·p), matching CYTools rounding.
    let mut grouped: HashMap<String, RacetrackTerm> = HashMap::new();
    for term in raw_terms {
        // CYTools uses round(q·p, 6)
        let key = format!("{:.6}", term.exponent.get());
        let entry = grouped.entry(key).or_insert(RacetrackTerm {
            curve: Vec::new(),
            coefficient: F64::<Finite>::ZERO,
            exponent: term.exponent,
        });
        // Finite + Finite = Finite
        entry.coefficient = entry.coefficient + term.coefficient;
    }

    // Sort by exponent (smallest first)
    let mut terms: Vec<RacetrackTerm> = grouped.into_values().collect();
    // Terms whose grouped coefficient vanishes (M·q = 0, or exact cancellation
    // within an exponent group) contribute nothing to W or its derivative;
    // keeping them would let a zero-coefficient leading exponent shadow the
    // actual racetrack pair (5-113-4627-main has one at q·p ≈ 0.12).
    terms.retain(|term| term.coefficient.get() != 0.0);
    terms.sort_by(|a, b| a.exponent.get().total_cmp(&b.exponent.get()));

    terms
}

/// Solve the racetrack equations for stabilization.
///
/// Implements the 2-term racetrack stabilization:
/// `d/dτ (A1 exp(-2π q1·p τ) + A2 exp(-2π q2·p τ)) = 0`
///
/// Returns `None` if no stable solution is found.
#[must_use]
pub fn solve_racetrack(terms: &[RacetrackTerm]) -> Option<RacetrackResult> {
    if terms.len() < 2 {
        return None;
    }

    // Use the two dominant terms (smallest exponents)
    let t1 = &terms[0];
    let t2 = &terms[1];

    // McAllister eqs. 2.25-2.26:
    // δ = -[(M·q1)(p·q1)N_q1] / [(M·q2)(p·q2)N_q2]
    // ε = (p·q2) - (p·q1)
    // Im(τ) = ln(1/|δ|)/(2πε), with Re(τ)=1/(2ε) for δ < 0.
    let epsilon = (t2.exponent - t1.exponent).try_to_pos()?;
    let derivative_coeff_1 = t1.exponent * t1.coefficient;
    let derivative_coeff_2 = t2.exponent * t2.coefficient;

    // Need NonZero to divide - check denominator isn't zero
    let derivative_coeff_2_nz = derivative_coeff_2.try_to_non_zero()?;
    let delta = -(derivative_coeff_1 / derivative_coeff_2_nz);
    let abs_delta = delta.abs().try_to_pos()?;
    let ln_inv_abs_delta = abs_delta.recip().ln().try_to_pos()?;
    let im_tau = ln_inv_abs_delta / (TWO_PI * epsilon);
    let g_s = im_tau.recip();

    // Validate g_s is in physical range (0, 1)
    if g_s.get() > 1.0 {
        return None;
    }

    let re_tau = if delta.get() < 0.0 {
        let branch = (f64_pos!(2.0) * epsilon).recip();
        F64::<NonNeg>::new(branch.get()).expect("1/(2 epsilon) is non-negative")
    } else {
        F64::<NonNeg>::ZERO
    };

    Some(RacetrackResult {
        g_s,
        re_tau,
        im_tau,
        delta: F64::<NonNeg>::ZERO,
        epsilon: F64::<NonNeg>::ZERO,
    })
}

fn complex_dilog_unit_disk(z_re: f64, z_im: f64) -> Option<(f64, f64)> {
    let z_abs_sq = z_re.mul_add(z_re, z_im * z_im);
    if !z_abs_sq.is_finite() || z_abs_sq >= 1.0 {
        return None;
    }
    if z_abs_sq == 0.0 {
        return Some((0.0, 0.0));
    }

    let mut sum_re = 0.0;
    let mut sum_im = 0.0;
    let mut pow_re = z_re;
    let mut pow_im = z_im;

    for k in 1..=DILOG_MAX_TERMS {
        let k_f = k as f64;
        let denom = k_f * k_f;
        let term_re = pow_re / denom;
        let term_im = pow_im / denom;
        sum_re += term_re;
        sum_im += term_im;

        let term_abs = term_re.hypot(term_im);
        let sum_abs = sum_re.hypot(sum_im).max(1.0);
        if term_abs <= DILOG_TOL * sum_abs {
            return Some((sum_re, sum_im));
        }

        let next_re = pow_re * z_re - pow_im * z_im;
        let next_im = pow_re * z_im + pow_im * z_re;
        if !next_re.is_finite() || !next_im.is_finite() {
            return None;
        }
        pow_re = next_re;
        pow_im = next_im;
    }

    None
}

fn racetrack_term_value(result: &RacetrackResult, term: &RacetrackTerm) -> Option<(f64, f64)> {
    // exp(2πiτ q·p) = exp(-2π Im(τ) q·p) × exp(i 2π Re(τ) q·p)
    let magnitude_arg = -TWO_PI * term.exponent * result.im_tau;
    let magnitude = magnitude_arg.get().exp();
    if !magnitude.is_finite() {
        return None;
    }
    if magnitude == 0.0 {
        return Some((0.0, 0.0));
    }
    let phase = (TWO_PI * term.exponent * result.re_tau).get();
    let z_re = magnitude * phase.cos();
    let z_im = magnitude * phase.sin();
    let (li2_re, li2_im) = complex_dilog_unit_disk(z_re, z_im)?;
    let coefficient = term.coefficient.get();
    Some((coefficient * li2_re, coefficient * li2_im))
}

/// Compute W₀ from all available racetrack terms.
///
/// From eq. 2.22:
/// ```text
/// W_flux = -ζ Σ_q (M·q) N_q Li₂(e^{2πiτ(q·p)})
/// ```
///
/// where ζ = 1/(2^{3/2} π^{5/2}) ≈ 0.02024 is the normalization factor.
///
/// Reference: arXiv:2107.09064, Eq. 2.22
#[must_use]
pub fn compute_w0_from_terms(
    result: &RacetrackResult,
    terms: &[RacetrackTerm],
) -> Option<Superpotential> {
    if terms.is_empty() {
        return None;
    }

    // W₀ = ζ × |sum of terms|, keeping the phase needed by same-sign racetracks.
    // The ζ normalization factor was missing - this caused ~50x error
    let mut sum_re = 0.0;
    let mut sum_im = 0.0;
    for term in terms {
        let (term_re, term_im) = racetrack_term_value(result, term)?;
        sum_re += term_re;
        sum_im += term_im;
    }
    let sum_abs = (sum_re.mul_add(sum_re, sum_im * sum_im)).sqrt();
    let sum_pos = F64::<Pos>::new(sum_abs)?;
    Some(ZETA_NORM * sum_pos)
}

/// Compute W₀ from the two leading stabilized racetrack terms.
///
/// Prefer [`compute_w0_from_terms`] when all computed terms are available.
///
/// # Panics
/// Panics if the two terms cancel exactly.
#[must_use]
pub fn compute_w0(
    result: &RacetrackResult,
    term1: &RacetrackTerm,
    term2: &RacetrackTerm,
) -> Superpotential {
    compute_w0_from_terms(result, &[term1.clone(), term2.clone()])
        .expect("racetrack terms never exactly cancel")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn finite_i64(v: i64) -> I64<Finite> {
        I64::<Finite>::new(v)
    }

    fn finite_f64(v: f64) -> F64<Finite> {
        F64::<Finite>::new(v).unwrap()
    }

    #[test]
    fn test_complex_dilog_known_real_value() {
        let (re, im) = complex_dilog_unit_disk(0.5, 0.0).unwrap();
        assert!((re - 0.582_240_526_465_012_5).abs() < 1e-14);
        assert!(im.abs() < 1e-14);
    }

    #[test]
    fn test_compute_w0_from_terms_uses_dilog_not_linear_approximation() {
        let result = RacetrackResult {
            g_s: f64_pos!(1.0),
            re_tau: F64::<NonNeg>::ZERO,
            im_tau: f64_pos!(std::f64::consts::LN_2 / (2.0 * PI)),
            delta: F64::<NonNeg>::ZERO,
            epsilon: F64::<NonNeg>::ZERO,
        };
        let terms = vec![RacetrackTerm {
            curve: vec![finite_i64(1)],
            coefficient: finite_f64(1.0),
            exponent: finite_f64(1.0),
        }];

        let w0 = compute_w0_from_terms(&result, &terms).unwrap();
        let linearized = ZETA_NORM.get() * 0.5;
        assert!(
            (w0.get() - linearized).abs() > 0.01 * ZETA_NORM.get(),
            "Li2(1/2) should be visibly different from the linear term"
        );
    }

    #[test]
    fn test_solve_racetrack_valid() {
        let terms = vec![
            RacetrackTerm {
                curve: vec![finite_i64(1)],
                coefficient: finite_f64(1.0),
                exponent: finite_f64(1.0),
            },
            RacetrackTerm {
                curve: vec![finite_i64(2)],
                coefficient: finite_f64(-1000.0),
                exponent: finite_f64(2.0),
            },
        ];

        let res = solve_racetrack(&terms).unwrap();
        assert!(res.g_s.get() > 0.0);
        assert!(res.g_s.get() < 1.0);
        assert_eq!(res.re_tau, F64::<NonNeg>::ZERO);
    }

    #[test]
    fn test_solve_racetrack_edge_case() {
        let terms = vec![RacetrackTerm {
            curve: vec![finite_i64(1)],
            coefficient: finite_f64(100.0),
            exponent: finite_f64(1.0),
        }];
        assert!(solve_racetrack(&terms).is_none());
    }

    #[test]
    fn test_zeta_constant() {
        assert!((ZETA_3.get() - 1.202_f64).abs() < 0.001);
    }
}
