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
    GvValue, ImTau, RacetrackCoefficient, RacetrackExponent, ResidualError, StringCoupling,
    Superpotential,
};
use crate::types::tags::{Finite, NonNeg, Pos};

/// Riemann zeta function at 3: ζ(3) ≈ 1.202056903159594 (positive).
pub const ZETA: F64<Pos> = f64_pos!(1.202_056_903_159_594);

/// 2π as a typed positive constant.
const TWO_PI: F64<Pos> = f64_pos!(2.0 * PI);

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

    // Group terms by exponent (q·p)
    let mut grouped: HashMap<String, RacetrackTerm> = HashMap::new();
    for term in raw_terms {
        // Use precision-aware key for grouping
        let key = format!("{:.10}", term.exponent.get());
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

    // Coefficients must have opposite signs for cancellation
    // Finite * Finite = Finite
    let product = t1.coefficient * t2.coefficient;
    if product.get() >= 0.0 {
        return None;
    }

    // Condition for dW/dτ = 0:
    // e^(-2π (q1-q2)·p / g_s) = - (q2 A2) / (q1 A1)
    // ratio = -(q1 A1) / (q2 A2)
    let numerator = t1.exponent * t1.coefficient;
    let denominator = t2.exponent * t2.coefficient;

    // Need NonZero to divide - check denominator isn't zero
    let denominator_nz = denominator.try_to_non_zero()?;
    // Finite / NonZero = Finite
    let ratio = -(numerator / denominator_nz);
    if ratio.get() <= 0.0 {
        return None;
    }

    // g_s = 2π × (q1·p - q2·p) / ln(ratio)
    // Finite - Finite = Finite
    let diff = t1.exponent - t2.exponent;

    // ln(ratio) - ratio is positive so ln is finite (and non-zero for ratio != 1)
    let ln_ratio = F64::<Finite>::new(ratio.get().ln())?.try_to_non_zero()?;

    // TWO_PI is Pos, diff is Finite, ln_ratio is NonZero
    // Pos * Finite = Finite, Finite / NonZero = Finite
    let g_s_raw = TWO_PI * diff / ln_ratio;

    // Validate g_s is in physical range (0, 1)
    if g_s_raw.get() <= 0.0 || g_s_raw.get() > 1.0 {
        return None;
    }

    // im_tau = q1·p / g_s
    // We already validated g_s_raw > 0, but type system doesn't know that yet
    // Since we narrow to Pos below, we can use try_to_pos here too
    let g_s_nz = g_s_raw.try_to_non_zero()?;
    // Finite / NonZero = Finite
    let im_tau_raw = t1.exponent / g_s_nz;

    // Narrow to positive types at the boundary
    let g_s = g_s_raw.try_to_pos()?;
    let im_tau = im_tau_raw.try_to_pos()?;

    Some(RacetrackResult {
        g_s,
        im_tau,
        delta: F64::<NonNeg>::ZERO,
        epsilon: F64::<NonNeg>::ZERO,
    })
}

/// Compute W₀ from stabilized racetrack.
///
/// `W_0 = - (A1 exp(-2π q1·p/gs) + A2 exp(-2π q2·p/gs))`
///
/// Reference: arXiv:2107.09064, Eq. 6.4
#[must_use]
pub fn compute_w0(
    result: &RacetrackResult,
    term1: &RacetrackTerm,
    term2: &RacetrackTerm,
) -> Superpotential {
    // -2π × exponent / g_s
    // -Pos = Neg, Neg * Finite = Finite, Finite / Pos = Finite
    let arg1 = -TWO_PI * term1.exponent / result.g_s;
    let arg2 = -TWO_PI * term2.exponent / result.g_s;

    // exp() of Finite is Pos (always positive)
    let exp1 = F64::<Pos>::new(arg1.get().exp()).expect("exp is always positive and finite");
    let exp2 = F64::<Pos>::new(arg2.get().exp()).expect("exp is always positive and finite");

    // Finite * Pos = Finite (algebra handles cross-type automatically)
    let term1_val = term1.coefficient * exp1;
    let term2_val = term2.coefficient * exp2;

    // -(Finite + Finite) = Finite
    -(term1_val + term2_val)
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
        assert!((ZETA.get() - 1.202_f64).abs() < 0.001);
    }
}
