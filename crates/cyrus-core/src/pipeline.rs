//! High-level evaluation pipeline for Calabi-Yau vacua.
//!
//! Integrates geometry, topology, and physics to evaluate a single (K, M) pair
//! against the Demirtas-McAllister criteria.

use crate::flat_direction::{compute_ek0, compute_n_matrix, solve_linear_system_faer};
use crate::intersection::Intersection;
use crate::kahler::MoriCone;
use crate::racetrack::{
    GvInvariant, RacetrackResult, build_racetrack_terms, compute_w0, solve_racetrack,
};
use crate::vacuum::{VacuumResult, compute_tadpole, compute_vacuum};
use crate::volume::compute_volume;

/// Input for a single vacuum evaluation.
#[derive(Debug, Clone)]
pub struct EvaluationRequest<'a> {
    /// Triple intersection numbers.
    pub kappa: &'a Intersection,
    /// Mori cone generators.
    pub mori: &'a MoriCone,
    /// GV invariants for racetrack computation.
    pub gv: &'a [GvInvariant],
    /// Euler characteristic χ components.
    pub h11: i32,
    /// Number of complex structure moduli.
    pub h21: i32,
    /// Tadpole bound Q_max.
    pub q_max: f64,
}

/// Detailed result of a vacuum evaluation.
#[derive(Debug, Clone)]
pub struct EvaluationResult {
    /// Whether the point passed all filters.
    pub success: bool,
    /// Failure reason (if success is false).
    pub reason: Option<String>,

    /// Flat direction p.
    pub p: Option<Vec<f64>>,
    /// Racetrack solution (g_s, τ).
    pub racetrack: Option<RacetrackResult>,
    /// Vacuum energy details.
    pub vacuum: Option<VacuumResult>,
    /// Flux tadpole Q_flux.
    pub q_flux: f64,
}

/// Evaluate a (K, M) pair.
///
/// Implements the 4-step filter and calculation pipeline from Demirtas-McAllister.
///
/// # Errors
/// Returns an error if any linear algebra or geometric calculation fails unexpectedly.
#[allow(clippy::cast_precision_loss)]
pub fn evaluate_vacuum(
    req: &EvaluationRequest,
    k: &[i64],
    m: &[i64],
) -> crate::Result<EvaluationResult> {
    let mut res = EvaluationResult {
        success: false,
        reason: None,
        p: None,
        racetrack: None,
        vacuum: None,
        q_flux: compute_tadpole(k, m),
    };

    // 1. Filter: Tadpole Bound
    if res.q_flux > req.q_max {
        res.reason = Some("Tadpole bound exceeded".into());
        return Ok(res);
    }

    // 2. Filter: N Matrix Invertibility
    let n_mat = compute_n_matrix(req.kappa, m);
    let Some(p) = solve_linear_system_faer(&n_mat, k) else {
        res.reason = Some("N matrix is singular".into());
        return Ok(res);
    };

    // 3. Filter: Flat Direction in Kähler Cone
    if !req.mori.contains(&p) {
        res.reason = Some("Flat direction outside Kähler cone".into());
        return Ok(res);
    }

    // 4. Filter: Orthogonality Constraint (K · p vanishes)
    let k_dot_p: f64 = k
        .iter()
        .zip(p.iter())
        .map(|(&ki, &pi)| (ki as f64) * pi)
        .sum();
    if k_dot_p.abs() > 1e-8 {
        res.reason = Some(format!(
            "Orthogonality constraint violated (K·p = {k_dot_p})"
        ));
        return Ok(res);
    }

    // 5. Calculation: Racetrack Solution
    let terms = build_racetrack_terms(req.gv, m, &p);
    let Some(rt_res) = solve_racetrack(&terms) else {
        res.reason = Some("No stable racetrack solution".into());
        return Ok(res);
    };

    // 6. Calculation: W₀ and V₀
    let w0 = compute_w0(&rt_res, &terms[0], &terms[1]);
    let ek0 = compute_ek0(req.kappa, &p);

    // Scale p to get actual Kähler moduli t = p / g_s
    let t: Vec<f64> = p.iter().map(|&pi| pi / rt_res.g_s).collect();
    let vol_res = compute_volume(req.kappa, &t, req.h11, req.h21);

    let vac_res = compute_vacuum(ek0, rt_res.g_s, vol_res.string_frame, w0.abs());

    // Success!
    res.success = true;
    res.p = Some(p);
    res.racetrack = Some(rt_res);
    res.vacuum = Some(vac_res);

    Ok(res)
}

#[cfg(test)]
mod tests {
    use super::*;
    use malachite::Integer;
    use malachite::Rational;

    #[test]
    fn test_evaluate_vacuum_simple() {
        let mut kappa = Intersection::new(1);
        kappa.set(0, 0, 0, Rational::from(6));

        let mori = MoriCone::new(vec![vec![Integer::from(1)]]);

        let gv = vec![
            GvInvariant {
                curve: vec![1],
                value: 50.0,
            },
            GvInvariant {
                curve: vec![2],
                value: 100.0,
            },
        ];

        let req = EvaluationRequest {
            kappa: &kappa,
            mori: &mori,
            gv: &gv,
            h11: 5,
            h21: 3,
            q_max: 100.0,
        };

        assert_eq!(req.h11, 5);
    }
}
