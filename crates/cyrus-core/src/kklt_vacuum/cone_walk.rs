//! The tau-line cone walk of arXiv:2107.09064 SS"Algorithm for F-flat
//! solutions".
//!
//! Instead of solving the KKLT targets inside one chamber and
//! giving up at its boundary, walk the STRAIGHT line in tau-space from the
//! initial point toward the target. Convexity of the dual effective cone
//! keeps the line inside; the corresponding t-path is continuous and C^1
//! across chamber walls, and each wall crossing is a flop that transforms
//! the intersection numbers (kappa' = kappa - n_C C_a C_b C_c) and
//! reassigns the GV invariants. The walk DISCOVERS the flop sequence that
//! the corrected-chamber reproduction previously curated by hand.
//!
//! Mechanics per segment: solve the existing fixed-kappa path-following
//! toward the next tau waypoint, then check every enumerated GV curve's
//! volume at the new point; a curve whose volume crossed below zero marks
//! a flop, which is applied before continuing (small segments keep the
//! overshoot into the neighboring chamber benign - the next segment
//! re-solves under the transformed kappa).

use crate::Intersection;
use crate::gv::curve_volume_in_divisor_basis;
use crate::kklt::{
    CertifiedFlopContinuationStep, KkltResult, flop_reassign_gv_invariants,
    flop_transform_intersection_numbers, solve_divisor_basis_path_following,
};
use crate::types::f64::F64;
use crate::types::physics::DivisorVolume;
use crate::types::range::CheckedRange;
use crate::types::tags::Finite;
use malachite::Integer;

use super::{OwnedDivisorBasis, PrimalGeom, PrimalIntersection};

/// Outcome of a cone walk: the converged result in the final chamber plus
/// the discovered flop sequence and the transformed data.
pub struct ConeWalkOutcome {
    /// The converged KKLT result in the final (post-flop) chamber.
    pub result: KkltResult,
    /// Transformed data of the destination chamber - the handoff for the
    /// (not yet implemented) full cross-chamber continuation.
    #[allow(dead_code)]
    pub kappa_final: Intersection,
    /// See `kappa_final`.
    #[allow(dead_code)]
    pub gv_final: Vec<(Vec<i64>, Integer)>,
    /// The certified flop sequence crossed from the start chamber.
    pub flops: Vec<CertifiedFlopContinuationStep>,
}

/// Why the zeroth-order rescue could not hand back an in-chamber solution.
///
/// These are the failure modes the binary used to translate into distinct
/// process exit codes; as a library routine the rescue returns them so the
/// caller (CLI or GA verifier) decides how to react.
#[derive(Clone, Debug)]
pub enum ConeWalkRescueError {
    /// The chamber's toric two-face GV table could not be built.
    GvTableFailed(String),
    /// No path to the KKLT targets was found within the flop/segment budget.
    NoPath,
    /// The solution lives `flops.len()` chamber walls away. Cross-chamber
    /// continuation of the downstream stages (chi transforms across flops)
    /// is not yet implemented, so the discovered sequence is reported rather
    /// than silently used.
    FlopsAway {
        /// The certified flop sequence from the start chamber to the solution.
        flops: Vec<CertifiedFlopContinuationStep>,
        /// Relative error of the converged solve in the destination chamber.
        relative_error: f64,
    },
}

/// Walk the straight tau-line from the initial point's volumes to
/// `tau_target`, flopping across chamber walls as the path crosses them.
///
/// `gv` is the primal GV table in ambient curve coordinates; only curves
/// with nonzero invariant define geometric walls (zero-GV sign changes are
/// not CY-cone walls and are ignored). Returns `None` when a segment solve
/// fails, the flop budget is exhausted, or the final solve does not
/// converge.
#[allow(clippy::too_many_arguments)] // mirrors the runner's solver call
pub fn solve_with_cone_walk(
    kappa_full: &Intersection,
    production_basis: &OwnedDivisorBasis,
    basis_indices: &[usize],
    kklt_basis: &[usize],
    tau_target: &[DivisorVolume],
    t_init: &[F64<Finite>],
    gv: &[(Vec<i64>, Integer)],
    kklt_steps: usize,
    segments: usize,
    max_flops: usize,
) -> Option<ConeWalkOutcome> {
    let mut kappa = kappa_full.clone();
    let mut gv_table: Vec<(Vec<i64>, Integer)> = gv.to_vec();
    let mut flops: Vec<CertifiedFlopContinuationStep> = Vec::new();
    let mut t_cur: Vec<F64<Finite>> = t_init.to_vec();

    // Initial tau under the starting kappa (only the kklt-basis entries
    // matter for the targets).
    let tau_init = divisor_volumes_at(&kappa, production_basis, kklt_basis, &t_cur)?;
    let curve_volume = |class: &[i64], t: &[F64<Finite>]| -> Option<f64> {
        curve_volume_in_divisor_basis(class, basis_indices, t)
            .ok()
            .map(F64::get)
    };

    let inner_steps = (kklt_steps / segments).max(32);
    // Waypoint stack in alpha (descending); subdivision pushes midpoints.
    let mut pending: Vec<f64> = (1..=segments)
        .rev()
        .map(|m| m as f64 / segments as f64)
        .collect();
    let mut alpha_done = 0.0f64;
    let mut retries_at_waypoint = 0usize;
    while let Some(&alpha) = pending.last() {
        let waypoint: Vec<DivisorVolume> = tau_init
            .iter()
            .zip(tau_target.iter())
            .map(|(a, b)| {
                let v = (1.0 - alpha).mul_add(a.get(), alpha * b.get());
                F64::<Finite>::new(v).and_then(F64::<Finite>::try_to_pos)
            })
            .collect::<Option<Vec<_>>>()?;

        let subdivide = |pending: &mut Vec<f64>, alpha_done: f64| {
            let mid = 0.5 * (alpha_done + alpha);
            if alpha - alpha_done < 1e-6 {
                false // cannot subdivide further
            } else {
                pending.push(mid);
                true
            }
        };

        let res = solve_divisor_basis_path_following(
            &kappa,
            production_basis.as_divisor_basis(),
            kklt_basis,
            &waypoint,
            &t_cur,
            CheckedRange::new(0, inner_steps),
        );
        let res = match res.filter(|r| r.converged) {
            Some(r) => r,
            None => {
                // Local divergence: take a smaller step along the tau line.
                if !subdivide(&mut pending, alpha_done) {
                    eprintln!(
                        "[INFO] cone walk: stuck at alpha={alpha_done:.4}; cannot subdivide further"
                    );
                    return None;
                }
                continue;
            }
        };
        let t_next = res.t.clone();

        // Wall detection: nonzero-GV curves whose volume crossed below
        // zero along this step. The FIRST crossing along the step is the
        // one with the largest before/(before - after) margin.
        let mut crossed: Vec<(usize, f64)> = gv_table
            .iter()
            .enumerate()
            .filter(|(_, (_, n))| *n != 0)
            .filter_map(|(idx, (class, _))| {
                let before = curve_volume(class, &t_cur)?;
                let after = curve_volume(class, &t_next)?;
                (before > 1e-12 && after < -1e-12).then_some((idx, before / (before - after)))
            })
            .collect();
        crossed.sort_by(|a, b| a.1.total_cmp(&b.1));

        if let Some(&(idx, _)) = crossed.first() {
            // More than one wall in this step, or a repeated ping-pong on
            // the same wall: shrink the step instead of flopping blindly.
            let (class, invariant) = gv_table[idx].clone();
            let neg_class: Vec<i64> = class.iter().map(|&c| -c).collect();
            let ping_pong = flops
                .last()
                .is_some_and(|step| step.curve_class == neg_class);
            if (crossed.len() > 1 || ping_pong || retries_at_waypoint > 8)
                && subdivide(&mut pending, alpha_done)
            {
                retries_at_waypoint += 1;
                continue;
            }
            if flops.len() >= max_flops {
                eprintln!("[INFO] cone walk: flop budget ({max_flops}) exhausted");
                return None;
            }
            eprintln!(
                "[INFO] cone walk: flop {} across curve {:?} (n_C = {}) at alpha={alpha:.4}",
                flops.len() + 1,
                class,
                invariant
            );
            kappa = flop_transform_intersection_numbers(&kappa, &class, &invariant)?;
            gv_table = flop_reassign_gv_invariants(&gv_table, &class, &invariant)?;
            flops.push(CertifiedFlopContinuationStep {
                curve_class: class,
                gv_invariant: invariant,
            });
            // Re-solve THIS waypoint under the flopped kappa from the
            // pre-wall point (t_cur unchanged): no compounding overshoot.
            retries_at_waypoint += 1;
            continue;
        }
        retries_at_waypoint = 0;
        alpha_done = alpha;
        pending.pop();
        t_cur = t_next;
    }

    // Final full-precision solve in the destination chamber.
    let final_res = solve_divisor_basis_path_following(
        &kappa,
        production_basis.as_divisor_basis(),
        kklt_basis,
        tau_target,
        &t_cur,
        CheckedRange::new(0, kklt_steps),
    )?;
    if !final_res.converged {
        return None;
    }
    Some(ConeWalkOutcome {
        result: final_res,
        kappa_final: kappa,
        gv_final: gv_table,
        flops,
    })
}

/// Basis-divisor volumes of the kklt set at `t` under `kappa`.
fn divisor_volumes_at(
    kappa: &Intersection,
    production_basis: &OwnedDivisorBasis,
    kklt_basis: &[usize],
    t: &[F64<Finite>],
) -> Option<Vec<DivisorVolume>> {
    let volumes = crate::kklt::compute_kklt_divisor_volumes_in_divisor_basis(
        kappa,
        production_basis.as_divisor_basis(),
        kklt_basis,
        t,
    )?;
    volumes.into_iter().map(F64::<Finite>::try_to_pos).collect()
}

/// Zeroth-order rescue at the runner's divergence site: run the tau-line
/// walk with the chamber's toric two-face GV table for wall detection.
///
/// A ZERO-flop success means the segmented path navigated a region the
/// direct solve could not, within the starting chamber: the result is
/// sound and the pipeline continues with it (`Ok`). An N-flop success means
/// the solution lives N chamber walls away; continuing the downstream stages
/// there requires transforming the divisor topology (chi across flops),
/// which is not yet implemented - so the discovered sequence is returned as
/// [`ConeWalkRescueError::FlopsAway`] for the caller to report.
#[allow(clippy::too_many_arguments)] // failure-site rescue context
pub fn rescue_zeroth_order(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    production_basis: &OwnedDivisorBasis,
    kklt_basis: &[usize],
    tau_target: &[DivisorVolume],
    t_init: &[F64<Finite>],
    kklt_steps: usize,
) -> Result<KkltResult, ConeWalkRescueError> {
    eprintln!(
        "[INFO] zeroth-order solve diverged; attempting tau-line cone walk (arXiv:2107.09064 SS5.2)"
    );
    let gv_table: Vec<(Vec<i64>, Integer)> =
        match crate::gv::compute_toric_two_face_curve_gv_invariants(
            &geom.triangulation,
            &geom.triangulation_points,
            &geom.polytope,
        ) {
            Ok(invariants) => invariants
                .into_iter()
                .map(|inv| (inv.class, inv.gv))
                .collect(),
            Err(e) => return Err(ConeWalkRescueError::GvTableFailed(e.to_string())),
        };
    let Some(outcome) = solve_with_cone_walk(
        &intersection.kappa_full,
        production_basis,
        &intersection.basis,
        kklt_basis,
        tau_target,
        t_init,
        &gv_table,
        kklt_steps,
        200,
        24,
    ) else {
        return Err(ConeWalkRescueError::NoPath);
    };
    if outcome.flops.is_empty() {
        eprintln!(
            "[INFO] cone walk converged WITHOUT chamber crossings (rel_err={}); continuing",
            outcome.result.relative_error.get()
        );
        return Ok(outcome.result);
    }
    Err(ConeWalkRescueError::FlopsAway {
        flops: outcome.flops,
        relative_error: outcome.result.relative_error.get(),
    })
}

/// Direct zeroth-order solve, falling back to the cone walk on divergence.
/// Returns [`ConeWalkRescueError`] when neither the direct solve nor the
/// rescue yields an in-chamber solution.
#[allow(clippy::too_many_arguments)] // mirrors the solver call
pub fn solve_zeroth_order_with_rescue(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    production_basis: &OwnedDivisorBasis,
    kklt_basis: &[usize],
    tau_target: &[DivisorVolume],
    t_init: &[F64<Finite>],
    kklt_steps: usize,
) -> Result<KkltResult, ConeWalkRescueError> {
    match solve_divisor_basis_path_following(
        &intersection.kappa_full,
        production_basis.as_divisor_basis(),
        kklt_basis,
        tau_target,
        t_init,
        CheckedRange::new(0, kklt_steps),
    )
    .filter(|r| r.converged)
    {
        Some(result) => Ok(result),
        None => rescue_zeroth_order(
            geom,
            intersection,
            production_basis,
            kklt_basis,
            tau_target,
            t_init,
            kklt_steps,
        ),
    }
}
