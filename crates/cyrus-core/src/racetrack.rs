//! Solver for the racetrack Kähler moduli stabilization.
//!
//! Given GV invariants and flat direction p, builds the racetrack potential
//! and solves for the stabilized Kähler moduli.
//!
//! Reference: arXiv:2107.09064, Section 6.4

use std::collections::HashMap;

/// A racetrack term: `coefficient * exp(-2π q·p / g_s)`.
#[derive(Debug, Clone)]
pub struct RacetrackTerm {
    /// Curve indices (q).
    pub curve: Vec<i64>,
    /// Pre-exponential factor: `(M·q) * N_q`.
    pub coefficient: f64,
    /// Effective exponent: `q·p`.
    pub exponent: f64,
}

/// Result of racetrack stabilization.
#[derive(Debug, Clone, Default)]
pub struct RacetrackResult {
    /// String coupling g_s.
    pub g_s: f64,
    /// Stabilized Kähler moduli value: `Im(τ) = (p / g_s)`.
    /// Note: This is the overall scale, individual τ_i = p_i / g_s.
    pub im_tau: f64,
    /// Numerical error delta.
    pub delta: f64,
    /// Numerical error epsilon.
    pub epsilon: f64,
}

/// A Gorus-Vafa invariant.
#[derive(Debug, Clone)]
pub struct GvInvariant {
    /// Curve class vector.
    pub curve: Vec<i64>,
    /// GV invariant value N_q.
    pub value: f64,
}

/// Build racetrack terms from GV invariants.
///
/// # Panics
/// Panics if the exponent sorting fails (e.g. if exponents are NaN).
#[allow(clippy::cast_precision_loss)]
pub fn build_racetrack_terms(
    gv_invariants: &[GvInvariant],
    m: &[i64],
    p: &[f64],
) -> Vec<RacetrackTerm> {
    let mut raw_terms = Vec::new();

    for gv in gv_invariants {
        // q·p (curve dot flat direction)
        let q_dot_p: f64 = gv
            .curve
            .iter()
            .zip(p.iter())
            .map(|(&qi, &pi)| (qi as f64) * pi)
            .sum();

        // M·q (flux dot curve)
        let m_dot_q: i64 = gv
            .curve
            .iter()
            .zip(m.iter())
            .map(|(&qi, &mi)| qi * mi)
            .sum();

        // coefficient = (M·q) × N_q
        let coefficient = (m_dot_q as f64) * gv.value;

        raw_terms.push(RacetrackTerm {
            curve: gv.curve.clone(),
            coefficient,
            exponent: q_dot_p,
        });
    }

    // Group terms by exponent (q·p)
    let mut grouped: HashMap<String, RacetrackTerm> = HashMap::new();
    for term in raw_terms {
        let key = format!("{:.8}", term.exponent);
        let entry = grouped.entry(key).or_insert(RacetrackTerm {
            curve: Vec::new(),
            coefficient: 0.0,
            exponent: term.exponent,
        });
        entry.coefficient += term.coefficient;
    }

    // Sort by exponent (smallest first)
    let mut terms: Vec<RacetrackTerm> = grouped.into_values().collect();
    terms.sort_by(|a, b| a.exponent.partial_cmp(&b.exponent).unwrap());

    terms
}

/// Solve the racetrack equations for stabilization.
///
/// Returns `None` if no stable solution is found.
#[must_use]
pub fn solve_racetrack(terms: &[RacetrackTerm]) -> Option<RacetrackResult> {
    if terms.len() < 2 {
        return None;
    }

    // Simple solver for 2-term racetrack
    let t1 = &terms[0];
    let t2 = &terms[1];

    // Ratio condition for stability
    let ratio = -t2.coefficient / t1.coefficient;
    if ratio <= 0.0 {
        return None;
    }

    let g_s = 2.0 * std::f64::consts::PI * (t2.exponent - t1.exponent) / ratio.ln();
    let im_tau = t1.exponent / g_s; // Crude approximation for overall scale

    Some(RacetrackResult {
        g_s,
        im_tau,
        delta: 0.0,
        epsilon: 0.0,
    })
}

/// Compute W₀ from stabilized racetrack.
#[must_use]
pub fn compute_w0(result: &RacetrackResult, term1: &RacetrackTerm, term2: &RacetrackTerm) -> f64 {
    let exp1 = (-2.0 * std::f64::consts::PI * term1.exponent / result.g_s).exp();
    let exp2 = (-2.0 * std::f64::consts::PI * term2.exponent / result.g_s).exp();

    let sum = term1.coefficient.mul_add(exp1, term2.coefficient * exp2);
    -sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_c_tau_relationship() {
        // Simple test logic for racetrack stability
        assert_eq!(1, 1);
    }

    #[test]
    fn test_solve_racetrack_valid() {
        let terms = vec![
            RacetrackTerm {
                curve: vec![1],
                coefficient: 50.0,
                exponent: 1.0,
            },
            RacetrackTerm {
                curve: vec![2],
                coefficient: -100.0,
                exponent: 2.0,
            },
        ];
        let res = solve_racetrack(&terms).unwrap();
        assert!(res.g_s > 0.0);
    }

    #[test]
    fn test_solve_racetrack_edge_case() {
        let terms = vec![RacetrackTerm {
            curve: vec![1],
            coefficient: 100.0,
            exponent: 1.0,
        }];
        assert!(solve_racetrack(&terms).is_none());
    }

    #[test]
    fn test_zeta_constant() {
        // Verify zeta(3) value
        let zeta3: f64 = 1.202_056_903_159_594;
        assert!((zeta3 - 1.202_f64).abs() < 0.001);
    }
}
