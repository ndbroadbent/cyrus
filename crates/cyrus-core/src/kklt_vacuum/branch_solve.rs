//! Phase-1 KKLT branch solving: finding the Kähler point that seeds the
//! worldsheet-instanton-corrected solve.
//!
//! Two routes mirror the runner. With `branch_candidates == 0` a single
//! height-projected initialization is path-followed to the phase-1 targets and
//! rescued by the tau-line cone walk on divergence. With `branch_candidates > 0`
//! a fan of scaled initializations (optionally prepended with the
//! height-projected one) is searched in parallel and one positive-volume branch
//! is selected by the requested policy. Either way the result is the
//! zeroth-order solve at the corrected target plus the two Kähler points the
//! downstream stages need.

use crate::kklt::{
    KkltResult, generate_scaled_divisor_basis_branch_initializations,
    solve_divisor_basis_path_following, solve_divisor_basis_path_following_branch_candidates,
};
use crate::types::f64::F64;
use crate::types::i64::I64;
use crate::types::physics::DivisorVolume;
use crate::types::range::CheckedRange;
use crate::types::tags::{Finite, Pos};

use super::basis_transform::{
    height_projected_branch_initialization, transform_production_primal_kahler_to_computed,
};
use super::cone_walk::{ConeWalkRescueError, solve_zeroth_order_with_rescue};
use super::{OwnedDivisorBasis, PrimalGeom, PrimalIntersection};

/// Policy for choosing among positive-volume phase-1 branches.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BranchSelection {
    /// Smallest classical volume.
    MinVolume,
    /// Largest classical volume.
    MaxVolume,
    /// Lowest initialization index (the height-projected seed, when present).
    FirstPositive,
    /// Best-conditioned final KKLT Jacobian.
    MinCondition,
}

impl BranchSelection {
    /// Stable string name for diagnostics.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::MinVolume => "min-volume",
            Self::MaxVolume => "max-volume",
            Self::FirstPositive => "first-positive",
            Self::MinCondition => "min-condition",
        }
    }
}

/// The phase-1 solve output: the zeroth-order solve at the corrected target,
/// plus the Kähler points used to select small curves and to seed the
/// GV-corrected iteration.
pub struct KkltBranchSolveOutput {
    /// Converged zeroth-order solve at the corrected (Euler-shifted) target.
    pub zeroth_order: KkltResult,
    /// Phase-1 Kähler point in Cyrus's computed primal basis, used to select
    /// the small toric curves.
    pub small_curve_selection_t: Vec<F64<Finite>>,
    /// Phase-1 Kähler point in the production basis, used to seed the
    /// GV-corrected solve iteration.
    pub gv_iteration_source_t: Vec<F64<Finite>>,
}

/// Render a cone-walk rescue failure as a one-line orchestrator error.
fn cone_walk_error_to_string(error: &ConeWalkRescueError) -> String {
    match error {
        ConeWalkRescueError::GvTableFailed(msg) => {
            format!("cone walk: two-face GV table failed: {msg}")
        }
        ConeWalkRescueError::NoPath => {
            "cone walk failed: no path to the KKLT targets found".to_string()
        }
        ConeWalkRescueError::FlopsAway {
            flops,
            relative_error,
        } => format!(
            "cone walk found the solution {} flops away (rel_err={relative_error}); cross-chamber continuation is not implemented",
            flops.len()
        ),
    }
}

/// Parameters controlling the phase-1 branch solve.
pub struct BranchSolveConfig {
    /// Path-following interpolation steps.
    pub kklt_steps: usize,
    /// Number of scaled branch initializations to search (0 = single height
    /// initialization with cone-walk rescue).
    pub branch_candidates: usize,
    /// Seed for the scaled branch initialization generator.
    pub branch_seed: u64,
    /// Prepend the height-projected initialization at index 0.
    pub branch_height_init: bool,
    /// Policy for selecting among positive-volume branches.
    pub branch_selection: BranchSelection,
}

/// Solve for the zeroth-order KKLT Kähler point feeding the corrected solve.
pub fn solve_kklt_branches(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    production_primal_basis: &OwnedDivisorBasis,
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
    tau_target: &[DivisorVolume],
    config: &BranchSolveConfig,
) -> Result<KkltBranchSolveOutput, String> {
    let tau_phase1: Vec<F64<Pos>> = c_i.iter().map(|ci| ci.to_f64()).collect();
    let steps = CheckedRange::new(0, config.kklt_steps);

    if config.branch_candidates == 0 {
        return solve_single_height_branch(
            geom,
            intersection,
            production_primal_basis,
            kklt_basis,
            &tau_phase1,
            tau_target,
            config.kklt_steps,
        );
    }

    let mut t_initializations = generate_scaled_divisor_basis_branch_initializations(
        &intersection.kappa_full,
        production_primal_basis.as_divisor_basis(),
        kklt_basis,
        &tau_phase1,
        config.branch_candidates,
        config.branch_seed,
    )
    .ok_or_else(|| "failed to generate KKLT branch initializations".to_string())?;
    if config.branch_height_init {
        let height_init = height_projected_branch_initialization(
            geom,
            intersection,
            production_primal_basis,
            kklt_basis,
            &tau_phase1,
        )?;
        t_initializations.insert(0, height_init);
    }

    let branch_search = solve_divisor_basis_path_following_branch_candidates(
        &intersection.kappa_full,
        production_primal_basis.as_divisor_basis(),
        kklt_basis,
        &tau_phase1,
        &t_initializations,
        steps,
    );
    eprintln!(
        "[INFO] KKLT branch search: attempted={} solved={} non_converged={} non_positive_volume={} positive_volume={}",
        branch_search.attempted,
        branch_search.solved,
        branch_search.non_converged,
        branch_search.non_positive_volume,
        branch_search.positive_volume.len()
    );
    let mut positive_branches = branch_search.positive_volume;
    positive_branches.sort_by(|a, b| {
        a.classical_volume
            .get()
            .total_cmp(&b.classical_volume.get())
    });
    if positive_branches.is_empty() {
        return Err("KKLT branch search found no positive-volume phase-1 branch".to_string());
    }

    let selected_rank = select_branch_rank(&positive_branches, config.branch_selection)?;
    let best_branch = positive_branches[selected_rank].clone();
    eprintln!(
        "[INFO] KKLT branch search selected policy={} rank_by_volume={} init={} phase1_volume={} rel_err={}",
        config.branch_selection.as_str(),
        selected_rank,
        best_branch.init_index,
        best_branch.classical_volume.get(),
        best_branch.result.relative_error.get()
    );

    let small_curve_selection_t = transform_production_primal_kahler_to_computed(
        intersection,
        production_primal_basis,
        &best_branch.result.t,
        "selected branch Kähler point",
    )?;

    let zeroth_order = solve_divisor_basis_path_following(
        &intersection.kappa_full,
        production_primal_basis.as_divisor_basis(),
        kklt_basis,
        tau_target,
        &best_branch.result.t,
        steps,
    )
    .ok_or_else(|| "corrected divisor-basis KKLT solve failed after branch search".to_string())?;
    if !zeroth_order.converged {
        return Err(format!(
            "corrected mixed-basis KKLT solve after branch search did not converge: rel_err={}",
            zeroth_order.relative_error.get()
        ));
    }

    Ok(KkltBranchSolveOutput {
        zeroth_order,
        small_curve_selection_t,
        gv_iteration_source_t: best_branch.result.t,
    })
}

/// The `branch_candidates == 0` route: a single height-projected
/// initialization, phase-1 path-followed and cone-walk rescued at the target.
fn solve_single_height_branch(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    production_primal_basis: &OwnedDivisorBasis,
    kklt_basis: &[usize],
    tau_phase1: &[F64<Pos>],
    tau_target: &[DivisorVolume],
    kklt_steps: usize,
) -> Result<KkltBranchSolveOutput, String> {
    let height_init = height_projected_branch_initialization(
        geom,
        intersection,
        production_primal_basis,
        kklt_basis,
        tau_phase1,
    )?;
    eprintln!("[INFO] using height-projected KKLT branch initialization");
    let phase1 = solve_divisor_basis_path_following(
        &intersection.kappa_full,
        production_primal_basis.as_divisor_basis(),
        kklt_basis,
        tau_phase1,
        &height_init,
        CheckedRange::new(0, kklt_steps),
    )
    .ok_or_else(|| "phase-1 divisor-basis KKLT path-following failed".to_string())?;
    if !phase1.converged {
        return Err(format!(
            "phase-1 mixed-basis KKLT path-following did not converge: rel_err={}",
            phase1.relative_error.get()
        ));
    }

    let zeroth_order = solve_zeroth_order_with_rescue(
        geom,
        intersection,
        production_primal_basis,
        kklt_basis,
        tau_target,
        &phase1.t,
        kklt_steps,
    )
    .map_err(|e| cone_walk_error_to_string(&e))?;

    let small_curve_selection_t = transform_production_primal_kahler_to_computed(
        intersection,
        production_primal_basis,
        &phase1.t,
        "phase-1 Kähler point",
    )?;

    Ok(KkltBranchSolveOutput {
        zeroth_order,
        small_curve_selection_t,
        gv_iteration_source_t: phase1.t,
    })
}

/// Choose the rank (index into the volume-sorted positive branches) per policy.
fn select_branch_rank(
    positive_branches: &[crate::kklt::KkltBranchSolution],
    branch_selection: BranchSelection,
) -> Result<usize, String> {
    match branch_selection {
        BranchSelection::MinVolume => Ok(0),
        BranchSelection::MaxVolume => Ok(positive_branches.len() - 1),
        BranchSelection::FirstPositive => positive_branches
            .iter()
            .enumerate()
            .min_by_key(|(_, branch)| branch.init_index)
            .map(|(rank, _)| rank)
            .ok_or_else(|| "positive branches unexpectedly empty".to_string()),
        BranchSelection::MinCondition => positive_branches
            .iter()
            .enumerate()
            .filter_map(|(rank, branch)| {
                branch
                    .jacobian_diagnostics
                    .condition_number
                    .map(|condition| (rank, condition.get()))
            })
            .min_by(|(_, a), (_, b)| a.total_cmp(b))
            .map(|(rank, _)| rank)
            .ok_or_else(|| {
                "KKLT branch search found no full-rank positive-volume phase-1 branch for min-condition selection"
                    .to_string()
            }),
    }
}
