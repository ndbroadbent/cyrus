//! Diagnostics for the admissible-basis scan: per-candidate coverage plus the
//! GV grading degrees of any uncovered curves.
//!
//! The deep-verify scan rejects a candidate basis on the FIRST uncovered small
//! curve, telling us nothing about how many uncovered curves a chamber has or
//! how expensive they would be to compute. This module measures exactly that,
//! so we can evaluate the tradeoff: "accept an early chamber with a few
//! low/moderate-degree uncovered curves and just compute their GV invariants"
//! vs "scan further to a fully toric-covered chamber." The general (HKTY) GV
//! cost grows steeply with a curve's grading degree, so the MAX degree of a
//! chamber's uncovered curves is what decides whether it is cheaply completable.

use crate::compute_grading_vector;

use super::PrimalIntersection;
use super::branch_solve::{BranchSolveConfig, solve_kklt_phase1_branch_search};
use super::gv_basis::vector_gv_basis_data;
use super::gv_coverage::{SmallCurveGeometry, collect_small_curve_gvs};
use super::verify::{StabilizationInputs, VacuumConfig};

/// Per-chamber coverage report at the (low-step) phase-1 point.
#[derive(Debug, Clone)]
pub struct ChamberCoverageReport {
    /// Phase-1 did not reach a positive-volume branch at the probe step count;
    /// the coverage verdict is unknown (this chamber would fall through to a
    /// full solve in the real scan).
    pub indeterminate: bool,
    /// Number of selected small curves with no cheap (toric / minimal-degree)
    /// GV invariant.
    pub uncovered_count: usize,
    /// GV grading degrees of the uncovered curves, ascending. The general-GV
    /// (HKTY) cost grows steeply with degree, so `uncovered_degrees.last()`
    /// (the max) decides whether the chamber is cheaply completable.
    pub uncovered_degrees: Vec<i128>,
}

/// Compute the primal GV grading vector once (basis/geometry-fixed), to be
/// reused across [`probe_chamber_coverage`] calls over many candidate bases.
///
/// # Errors
/// Returns an error if the GV basis data or grading vector cannot be built.
pub fn primal_gv_grading(
    scgeom: &SmallCurveGeometry,
    intersection: &PrimalIntersection,
) -> Result<Vec<i64>, String> {
    let gv_basis_data = vector_gv_basis_data(
        &scgeom.ambient_rays,
        &intersection.linrels,
        &intersection.basis,
        "primal-diagnostic",
    )?;
    compute_grading_vector(&gv_basis_data.mori_rays)
        .ok_or_else(|| "failed to compute primal GV grading vector".to_string())
}

/// Grading degree of an ambient curve class: the basis-restricted dot product
/// of the class with the GV grading vector (mirrors the general-GV path).
fn ambient_curve_grading_degree(class: &[i64], basis: &[usize], grading: &[i64]) -> i128 {
    basis
        .iter()
        .zip(grading.iter())
        .map(|(&idx, &weight)| {
            i128::from(class.get(idx).copied().unwrap_or(0)) * i128::from(weight)
        })
        .sum()
}

/// Measure one candidate basis's chamber.
///
/// Reports how many small curves the cheap GV layers leave uncovered at the
/// low-step phase-1 point, and their grading degrees. Runs only the phase-1
/// branch search (no corrected solve) at `probe_steps`, then the cheap coverage
/// collection (general GV OFF), then the degree projection. Reuses the shared,
/// basis-independent `grading`.
#[must_use]
pub fn probe_chamber_coverage(
    scgeom: &SmallCurveGeometry,
    inputs: &StabilizationInputs<'_>,
    config: &VacuumConfig,
    grading: &[i64],
    probe_steps: usize,
) -> ChamberCoverageReport {
    let probe_branch = BranchSolveConfig {
        kklt_steps: probe_steps,
        branch_candidates: config.branch_candidates,
        branch_seed: config.branch_seed,
        branch_height_init: config.branch_height_init,
        branch_selection: config.branch_selection,
    };
    let indeterminate = ChamberCoverageReport {
        indeterminate: true,
        uncovered_count: 0,
        uncovered_degrees: Vec::new(),
    };
    let Ok(phase1) = solve_kklt_phase1_branch_search(
        inputs.geom,
        inputs.intersection,
        inputs.production_primal_basis,
        inputs.kklt_basis,
        inputs.c_i,
        &probe_branch,
        true,
    ) else {
        return indeterminate;
    };
    let Ok(selection) = collect_small_curve_gvs(
        scgeom,
        inputs.geom,
        inputs.intersection,
        inputs.kklt_basis,
        &phase1.small_curve_selection_t,
        &config.small_curve_gv_config(),
    ) else {
        return indeterminate;
    };
    let basis = &inputs.intersection.basis;
    let mut uncovered_degrees: Vec<i128> = selection
        .uncovered_missing
        .iter()
        .map(|class| ambient_curve_grading_degree(class, basis, grading))
        .collect();
    uncovered_degrees.sort_unstable();
    ChamberCoverageReport {
        indeterminate: false,
        uncovered_count: selection.uncovered_missing.len(),
        uncovered_degrees,
    }
}
