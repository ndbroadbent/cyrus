//! `verify_kklt_vacuum`: the one-shot KKLT vacuum verification entry point.
//!
//! Given a primal geometry, its intersection data, the orientifold/flux model
//! (`c_i`, KKLT basis) and the racetrack-derived scalars (`g_s`, `W0`, `e^{K0}`),
//! this assembles the full arXiv:2107.09064 SS6 stabilization: the corrected
//! divisor-volume target, the phase-1 branch solve, the small-curve GV
//! coverage, the worldsheet-instanton-corrected solve, and the BBHL/GV-corrected
//! string-frame volume, ending at the vacuum energy `V0`.
//!
//! It is the library product the `mcallister_first_principles` binary and the
//! GA verifier both call in-process; the binary's `stage_volume`/`stage_vacuum`
//! are the intermediate scaffolding it replaces.

use crate::divisor::{batyrev_h11_extra_classes, compute_kklt_divisor_chi};
use crate::kklt::{
    compute_c_tau, compute_corrected_target_tau, compute_gv_volume_correction_for_ambient_curves,
};
use crate::types::f64::F64;
use crate::types::i32::I32;
use crate::types::i64::I64;
use crate::types::tags::{Finite, GTEOne, Neg, NonNeg, Pos};
use crate::vacuum::compute_vacuum;
use crate::volume::bbhl_correction;
use crate::{CurvePruningStrategy, Intersection};

use super::branch_solve::{BranchSelection, BranchSolveConfig, solve_kklt_branches};
use super::gamma::compute_b_field_gamma_for_o7_divisors;
use super::gv_corrected::solve_gv_corrected_kklt;
use super::gv_coverage::{
    SmallCurveGeometry, SmallCurveGvConfig, compute_small_curve_geometry, select_small_curve_gvs,
};
use super::missing_gv::verify_deferred_missing_gv_bounds;
use super::{OwnedDivisorBasis, PrimalGeom, PrimalIntersection};

/// Primal-geometry + model inputs to the KKLT vacuum verification.
pub struct StabilizationInputs<'a> {
    /// Primal toric geometry (polytope, heights, triangulation).
    pub geom: &'a PrimalGeom,
    /// Primal intersection data (GLSM, linrels, basis, intersection tensors).
    pub intersection: &'a PrimalIntersection,
    /// Production divisor basis for the KKLT solve (e.g. McAllister's
    /// `basis.dat`, or Cyrus's computed basis - see
    /// [`computed_primal_basis`](super::computed_primal_basis)).
    pub production_primal_basis: &'a OwnedDivisorBasis,
    /// Dual Coxeter numbers `c_i` of the gaugino-condensate stacks.
    pub c_i: &'a [I64<Pos>],
    /// Divisor indices contributing to the superpotential (the KKLT basis).
    pub kklt_basis: &'a [usize],
    /// String coupling `g_s` from the racetrack.
    pub g_s: F64<Pos>,
    /// Flux superpotential magnitude `|W0|` from the racetrack.
    pub w0: F64<Pos>,
    /// Kähler-potential factor `e^{K0}` from the flat direction.
    pub ek0: F64<Pos>,
    /// Hodge number `h^{2,1}` (for the BBHL Euler characteristic).
    pub h21: usize,
}

/// Knobs controlling the KKLT solve and small-curve GV coverage.
pub struct VacuumConfig {
    /// Path-following interpolation steps.
    pub kklt_steps: usize,
    /// Number of scaled branch initializations (0 = single height init).
    pub branch_candidates: usize,
    /// Seed for the scaled branch initialization generator.
    pub branch_seed: u64,
    /// Prepend the height-projected initialization at index 0.
    pub branch_height_init: bool,
    /// Policy for selecting among positive-volume branches.
    pub branch_selection: BranchSelection,
    /// Volume cutoff for selecting small toric curves.
    pub small_curve_cutoff: F64<Pos>,
    /// Strategy for pruning decomposable curve classes.
    pub small_curve_pruning: CurvePruningStrategy,
    /// Optional minimum lattice-point count for the general GV fallback.
    pub primal_gv_min_points: Option<u32>,
    /// Optional maximum grading degree for the general GV fallback.
    pub primal_gv_max_deg: Option<u32>,
    /// Maximum permitted impact for deferring an uncovered curve (0 disables).
    pub max_missing_gv_impact: f64,
    /// Assumed `|GV|` bound when bounding a deferred curve's impact.
    pub missing_gv_abs_bound: u32,
}

impl VacuumConfig {
    const fn small_curve_gv_config(&self) -> SmallCurveGvConfig {
        SmallCurveGvConfig {
            small_curve_cutoff: self.small_curve_cutoff,
            small_curve_pruning: self.small_curve_pruning,
            primal_gv_min_points: self.primal_gv_min_points,
            primal_gv_max_deg: self.primal_gv_max_deg,
            max_missing_gv_impact: self.max_missing_gv_impact,
            missing_gv_abs_bound: self.missing_gv_abs_bound,
        }
    }
}

/// The verdict of a KKLT vacuum verification: the stabilized moduli, the
/// string-frame volume breakdown, and the vacuum energy `V0`.
pub struct VacuumVerdict {
    /// Vacuum energy `V0 = -3 e^{K0} (g_s^7/(4V)^2) |W0|^2` (always negative).
    pub v0: F64<Neg>,
    /// String-frame Calabi-Yau volume (classical - BBHL + GV correction).
    pub v_string: F64<Pos>,
    /// String coupling `g_s`.
    pub g_s: F64<Pos>,
    /// Flux superpotential `|W0|`.
    pub w0: F64<Pos>,
    /// Kähler-potential factor `e^{K0}`.
    pub ek0: F64<Pos>,
    /// Stabilized Kähler moduli `t` in Cyrus's computed primal basis.
    pub t_moduli: Vec<F64<Finite>>,
    /// Classical volume `(1/6) κ_ijk t^i t^j t^k`.
    pub classical_volume: f64,
    /// BBHL `α'` correction subtracted from the classical volume.
    pub bbhl_correction: f64,
    /// Worldsheet-instanton (GV) volume correction.
    pub gv_volume_correction: f64,
    /// Number of small toric curves whose GV invariants entered the correction.
    pub small_curve_gv_count: usize,
    /// Whether every selected small curve was covered (no curve deferred): the
    /// instanton expansion is fully under control.
    pub gv_controlled: bool,
    /// Number of small curves deferred under the bounded-impact policy.
    pub deferred_missing_gv_count: usize,
    /// The KKLT divisor basis that produced this verdict. For a single-basis
    /// verification this is the input basis; for the auto-basis scan it is the
    /// winning admissible basis (the first whose chamber the cheap GV layers
    /// fully covered). The stabilized vacuum - hence `V0` - is the same for
    /// every admissible basis; this records WHICH chamber the homotopy took.
    pub selected_kklt_basis: Vec<usize>,
}

/// Verify a KKLT vacuum end-to-end: stabilize the Kähler moduli with the
/// worldsheet-instanton-corrected solve and compute the vacuum energy.
///
/// Returns `Err` (rather than a meaningless vacuum) on any failure of the
/// chain: target construction, branch solve, GV coverage, the corrected solve,
/// the deferred-bound re-verification, the volume correction, or a `V0` that
/// underflows `f64`. No physics step silently falls back.
pub fn verify_kklt_vacuum(
    inputs: &StabilizationInputs<'_>,
    config: &VacuumConfig,
) -> Result<VacuumVerdict, String> {
    let scgeom = compute_small_curve_geometry(inputs.geom)?;
    verify_one_basis(&scgeom, inputs, config)
}

/// Verify a KKLT vacuum for ONE candidate basis, reusing the precomputed
/// basis-independent [`SmallCurveGeometry`]. Errors (before the expensive
/// corrected solve) when the basis lands in a chamber whose small curves the
/// cheap GV methods cannot cover - so a basis search can try the next candidate
/// without paying for the corrected solve.
fn verify_one_basis(
    scgeom: &SmallCurveGeometry,
    inputs: &StabilizationInputs<'_>,
    config: &VacuumConfig,
) -> Result<VacuumVerdict, String> {
    let StabilizationInputs {
        geom,
        intersection,
        production_primal_basis,
        c_i,
        kklt_basis,
        g_s,
        w0,
        ek0,
        h21,
    } = *inputs;

    let c_tau = compute_c_tau(g_s, w0);
    let chi_divisor = compute_kklt_divisor_chi(
        &geom.polytope,
        &geom.triangulation_points,
        &intersection.kappa_full,
        kklt_basis,
    )
    .map_err(|e| format!("failed to compute KKLT divisor chi: {e}"))?;
    let gamma = compute_b_field_gamma_for_o7_divisors(&geom.triangulation_points, kklt_basis, c_i)?;
    let tau_target = compute_corrected_target_tau(c_i, &chi_divisor, c_tau)
        .ok_or_else(|| "corrected KKLT target construction failed".to_string())?;

    let branch = solve_kklt_branches(
        geom,
        intersection,
        production_primal_basis,
        kklt_basis,
        c_i,
        &tau_target,
        &BranchSolveConfig {
            kklt_steps: config.kklt_steps,
            branch_candidates: config.branch_candidates,
            branch_seed: config.branch_seed,
            branch_height_init: config.branch_height_init,
            branch_selection: config.branch_selection,
        },
    )?;
    eprintln!(
        "[INFO] zeroth-order mixed-basis KKLT converged={} rel_err={}",
        branch.zeroth_order.converged,
        branch.zeroth_order.relative_error.get()
    );

    let gv_selection = select_small_curve_gvs(
        scgeom,
        geom,
        intersection,
        kklt_basis,
        &branch.small_curve_selection_t,
        &config.small_curve_gv_config(),
    )?;

    let correction_source_t = solve_gv_corrected_kklt(
        intersection,
        production_primal_basis,
        kklt_basis,
        c_i,
        &chi_divisor,
        c_tau,
        &gv_selection.small_curve_gvs,
        &gamma,
        &branch.gv_iteration_source_t,
        config.kklt_steps,
    )?;

    verify_deferred_missing_gv_bounds(
        &gv_selection.deferred_missing,
        intersection,
        kklt_basis,
        &correction_source_t,
        config.max_missing_gv_impact,
        config.missing_gv_abs_bound,
    )?;

    let gv_volume_correction = compute_gv_volume_correction_for_ambient_curves(
        &gv_selection.small_curve_gvs,
        &intersection.basis,
        &correction_source_t,
        Some(&gamma),
    )
    .ok_or_else(|| "failed to compute primal ambient GV volume correction".to_string())?;
    eprintln!(
        "[INFO] GV volume correction = {}",
        gv_volume_correction.get()
    );

    let t = correction_source_t;
    let classical_volume = classical_volume_from_t(&intersection.kappa_basis, &t);
    let bbhl = bbhl_correction_for_geometry(geom, intersection, h21)?;
    let v_string_value = classical_volume - bbhl + gv_volume_correction.get();
    eprintln!("[INFO] V_classical = {classical_volume}");
    eprintln!("[INFO] BBHL correction = {bbhl}");
    eprintln!("[INFO] V_string = {v_string_value}");
    let v_string = F64::<Pos>::new(v_string_value)
        .ok_or_else(|| format!("V_string must be positive, got {v_string_value}"))?;

    let vacuum = compute_vacuum(ek0, g_s, v_string, w0)
        .ok_or_else(|| "V0 underflowed f64; |W0| too small to represent".to_string())?;

    Ok(VacuumVerdict {
        v0: vacuum.v0,
        v_string,
        g_s,
        w0,
        ek0,
        t_moduli: t,
        classical_volume,
        bbhl_correction: bbhl,
        gv_volume_correction: gv_volume_correction.get(),
        small_curve_gv_count: gv_selection.small_curve_gvs.len(),
        gv_controlled: gv_selection.deferred_missing.is_empty(),
        deferred_missing_gv_count: gv_selection.deferred_missing.len(),
        selected_kklt_basis: kklt_basis.to_vec(),
    })
}

/// Verify a KKLT vacuum, choosing the KKLT basis automatically from a list of
/// admissible `(kklt_basis, c_i)` candidates.
///
/// All admissible bases describe the same physical vacuum, but they place the
/// phase-1 selection point in different chambers, and only some land where the
/// cheap GV methods cover every selected small curve. This tries the candidates
/// in order, paying the dominant Mori-cap / toric-GV cost ONCE (shared
/// [`SmallCurveGeometry`]) and only the light per-basis solve+coverage for each,
/// and returns the first candidate that verifies cleanly. The `c_i`/`kklt_basis`
/// fields of `base` are ignored (each candidate supplies its own); the other
/// fields (geometry, racetrack scalars, `h21`) are shared.
///
/// # Errors
/// Returns an error if the shared geometry cannot be built, or if no candidate
/// verifies (the last candidate's error is included).
pub fn verify_kklt_vacuum_auto_basis(
    base: &StabilizationInputs<'_>,
    candidates: &[(Vec<usize>, Vec<I64<Pos>>)],
    config: &VacuumConfig,
) -> Result<VacuumVerdict, String> {
    if candidates.is_empty() {
        return Err("no candidate KKLT bases supplied".to_string());
    }
    let scgeom = compute_small_curve_geometry(base.geom)?;
    let mut last_err = String::new();
    for (i, (kklt_basis, c_i)) in candidates.iter().enumerate() {
        let inputs = StabilizationInputs {
            geom: base.geom,
            intersection: base.intersection,
            production_primal_basis: base.production_primal_basis,
            c_i,
            kklt_basis,
            g_s: base.g_s,
            w0: base.w0,
            ek0: base.ek0,
            h21: base.h21,
        };
        match verify_one_basis(&scgeom, &inputs, config) {
            Ok(verdict) => {
                eprintln!(
                    "[INFO] auto-basis: candidate {}/{} verified cleanly",
                    i + 1,
                    candidates.len()
                );
                return Ok(verdict);
            }
            Err(e) => {
                eprintln!("[INFO] auto-basis: candidate {} rejected: {e}", i + 1);
                last_err = e;
            }
        }
    }
    Err(format!(
        "no admissible KKLT basis among {} candidates verified cleanly; last error: {last_err}",
        candidates.len()
    ))
}

/// Classical Calabi-Yau volume `(1/6) κ_ijk t^i t^j t^k` in the divisor basis.
fn classical_volume_from_t(kappa_basis: &Intersection, t: &[F64<Finite>]) -> f64 {
    let mut volume_sum = 0.0f64;
    for (&(i, j, k), val) in kappa_basis.iter() {
        let kappa_val = val.to_f64().get();
        let t_product = t[i].get() * t[j].get() * t[k].get();
        let mult = if i == j && j == k {
            1.0
        } else if i == j || j == k || i == k {
            3.0
        } else {
            6.0
        };
        volume_sum += mult * kappa_val * t_product;
    }
    volume_sum / 6.0
}

/// BBHL `α'` correction `ζ(3) χ / (4(2π)^3)`, accounting for non-favorable
/// split-divisor classes when counting `h^{1,1}`.
fn bbhl_correction_for_geometry(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    h21: usize,
) -> Result<f64, String> {
    let h11_extra = batyrev_h11_extra_classes(&geom.polytope, &geom.triangulation_points)
        .map_err(|e| format!("failed to compute Batyrev h11 correction: {e}"))?;
    let h11_raw = i32::try_from(intersection.basis.len() + h11_extra)
        .map_err(|_| "h11 does not fit in i32".to_string())?;
    let h21_raw = i32::try_from(h21).map_err(|_| "h21 does not fit in i32".to_string())?;
    let h11 = I32::<GTEOne>::new(h11_raw).ok_or_else(|| "computed h11 must be >= 1".to_string())?;
    let h21 = I32::<NonNeg>::new(h21_raw)
        .ok_or_else(|| "computed h21 must be non-negative".to_string())?;
    Ok(bbhl_correction(h11, h21).get())
}
