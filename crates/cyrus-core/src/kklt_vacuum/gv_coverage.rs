//! Selecting and covering the GV invariants of the small toric curves that
//! enter the worldsheet-instanton volume correction.
//!
//! The selection cascade mirrors arXiv:2107.09064 SS6: enumerate the
//! sub-cutoff toric curves at the chosen Kähler point, prune decomposable
//! classes, then cover each surviving class with a GV invariant - first from
//! the toric two-face formulas, then a minimal-degree HKTY pre-pass, then the
//! optional general-degree fallback. Curves whose maximal possible
//! contribution is provably negligible may be deferred; anything still
//! uncovered is a hard failure (no silent omission of an instanton term).

use std::collections::HashMap;

use malachite::Integer;

use crate::{
    CurvePruningStrategy, ToricCurveCandidate, compute_mori_cone_cap_rays,
    compute_toric_two_face_curve_gv_invariants, prune_decomposable_curve_candidates,
    subcutoff_toric_curve_candidates,
};

use crate::types::f64::F64;
use crate::types::tags::{Finite, Pos};

use super::general_gv::compute_primal_general_gv_by_ambient_class;
use super::missing_gv::{
    compute_minimal_degree_gv_by_ambient_class, defer_bounded_impact_missing_gvs,
};
use super::{PrimalGeom, PrimalIntersection};

/// The GV invariants selected for the small toric curves, plus the curves
/// whose negligible contribution was deferred (and must be re-verified at the
/// solved Kähler point by the caller).
pub struct SmallCurveGvSelection {
    /// Selected `(ambient curve class, GV invariant)` pairs feeding the volume
    /// correction.
    pub small_curve_gvs: Vec<(Vec<i64>, Integer)>,
    /// Curves deferred under the bounded-impact policy; excluded from the
    /// correction but re-verified at the solved point.
    pub deferred_missing: Vec<Vec<i64>>,
}

/// Parameters controlling the small-curve GV selection cascade.
pub struct SmallCurveGvConfig {
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

/// Collect GV invariants for the small curves from `gv_by_class`, returning
/// the covered `(class, gv)` pairs and the still-missing classes.
fn partition_covered(
    small_curves: &[ToricCurveCandidate],
    gv_by_class: &HashMap<Vec<i64>, Integer>,
) -> (Vec<(Vec<i64>, Integer)>, Vec<Vec<i64>>) {
    let mut covered = Vec::with_capacity(small_curves.len());
    let mut missing = Vec::new();
    for curve in small_curves {
        match gv_by_class.get(&curve.class) {
            Some(gv) => covered.push((curve.class.clone(), gv.clone())),
            None => missing.push(curve.class.clone()),
        }
    }
    (covered, missing)
}

/// Merge newly computed `(class, gv)` values into `gv_by_class`, failing on a
/// conflicting value for an already-known class.
fn merge_gv_values(
    gv_by_class: &mut HashMap<Vec<i64>, Integer>,
    new_values: HashMap<Vec<i64>, Integer>,
    context: &str,
) -> Result<(), String> {
    for (class, gv) in new_values {
        match gv_by_class.entry(class) {
            std::collections::hash_map::Entry::Occupied(existing) => {
                if existing.get() != &gv {
                    return Err(format!(
                        "toric/{context} GV conflict for selected primal curve: {} vs {gv}",
                        existing.get()
                    ));
                }
            }
            std::collections::hash_map::Entry::Vacant(slot) => {
                slot.insert(gv);
            }
        }
    }
    Ok(())
}

/// Run the full small-curve GV selection cascade at `small_curve_selection_t`
/// (a Kähler point in Cyrus's computed primal basis).
///
/// Returns the selected GV invariants and the deferred classes. Fails when any
/// selected small curve remains uncovered after the toric, minimal-degree, and
/// (optional) general-degree layers and is not eligible for deferral.
pub fn select_small_curve_gvs(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    kklt_basis: &[usize],
    small_curve_selection_t: &[F64<Finite>],
    config: &SmallCurveGvConfig,
) -> Result<SmallCurveGvSelection, String> {
    let ambient_rays = compute_mori_cone_cap_rays(
        &geom.triangulation,
        &geom.triangulation_points,
        &geom.polytope,
        false,
        false,
        None,
    )
    .map_err(|e| format!("failed to compute primal ambient Mori-cap rays: {e}"))?;
    let small_curve_candidates = subcutoff_toric_curve_candidates(
        &ambient_rays,
        &intersection.basis,
        small_curve_selection_t,
        config.small_curve_cutoff,
    )
    .map_err(|e| format!("failed to select small toric curve candidates: {e}"))?;
    let small_curves =
        prune_decomposable_curve_candidates(&small_curve_candidates, config.small_curve_pruning)
            .map_err(|e| {
                format!(
                    "failed to prune primal small curves with {} strategy: {e}",
                    config.small_curve_pruning.as_str()
                )
            })?;
    eprintln!(
        "[INFO] primal small toric curves: ambient_rays={} subcutoff={} pruned={} pruning={} cutoff={}",
        ambient_rays.len(),
        small_curve_candidates.len(),
        small_curves.len(),
        config.small_curve_pruning.as_str(),
        config.small_curve_cutoff.get()
    );

    let toric_gvs = compute_toric_two_face_curve_gv_invariants(
        &geom.triangulation,
        &geom.triangulation_points,
        &geom.polytope,
    )
    .map_err(|e| format!("failed to compute primal toric curve GV values: {e}"))?;
    let mut gv_by_class: HashMap<Vec<i64>, Integer> = toric_gvs
        .into_iter()
        .map(|item| (item.class, item.gv))
        .collect();

    let (mut small_curve_gvs, mut missing_gv_classes) =
        partition_covered(&small_curves, &gv_by_class);

    // Minimal-degree HKTY pre-pass for classes the toric formulas missed.
    if !missing_gv_classes.is_empty() {
        let minimal_degree_gvs = compute_minimal_degree_gv_by_ambient_class(
            geom,
            intersection,
            &missing_gv_classes,
            &ambient_rays,
        )
        .map_err(|e| format!("minimal-degree GV pre-pass failed: {e}"))?;
        if !minimal_degree_gvs.is_empty() {
            merge_gv_values(&mut gv_by_class, minimal_degree_gvs, "minimal-degree")?;
            (small_curve_gvs, missing_gv_classes) = partition_covered(&small_curves, &gv_by_class);
        }
    }

    // Optional general-degree fallback (opt-in via min-points / max-deg).
    if !missing_gv_classes.is_empty()
        && (config.primal_gv_min_points.is_some() || config.primal_gv_max_deg.is_some())
    {
        eprintln!(
            "[INFO] toric formulas missed {} selected primal small curves; computing primal general GV fallback with min_points={:?} max_deg={:?}",
            missing_gv_classes.len(),
            config.primal_gv_min_points,
            config.primal_gv_max_deg
        );
        let general_gvs = compute_primal_general_gv_by_ambient_class(
            geom,
            intersection,
            &missing_gv_classes,
            Some(&ambient_rays),
            config.primal_gv_min_points,
            config.primal_gv_max_deg,
        )
        .map_err(|e| format!("failed to compute primal general GV fallback: {e}"))?;
        merge_gv_values(&mut gv_by_class, general_gvs, "general")?;
        (small_curve_gvs, missing_gv_classes) = partition_covered(&small_curves, &gv_by_class);
    }

    let deferred_missing = defer_bounded_impact_missing_gvs(
        &mut missing_gv_classes,
        intersection,
        kklt_basis,
        small_curve_selection_t,
        config.max_missing_gv_impact,
        config.missing_gv_abs_bound,
    );

    if !missing_gv_classes.is_empty() {
        return Err(format!(
            "{} selected primal small curve(s) have no GV invariant and exceed the deferral bound; first uncovered class={:?}",
            missing_gv_classes.len(),
            missing_gv_classes.first()
        ));
    }

    eprintln!(
        "[INFO] primal small-curve GVs selected={} deferred={}",
        small_curve_gvs.len(),
        deferred_missing.len()
    );

    Ok(SmallCurveGvSelection {
        small_curve_gvs,
        deferred_missing,
    })
}
