//! In-process deep verification of a flux candidate's KKLT vacuum.
//!
//! This runs the full primal-side stabilization the mirror scan only proxies,
//! and only on the rare candidates worth it (see below).
//!
//! The mirror-side scan (`cyrus_core::evaluate_vacuum`) is a fast proxy: it
//! racetrack-stabilizes `g_s`/`|W0|`, computes `e^{K0}` from the mirror flat
//! direction, and reports a proxy `V0` from a mirror-side volume. It is cheap
//! enough to run on every genome. What it does NOT do is the full
//! arXiv:2107.09064 SS6 stabilization on the *primal* Calabi-Yau X: derive the
//! orientifold, pick a KKLT divisor basis, solve the worldsheet-instanton-
//! corrected F-flat conditions in a covered chamber, and read off the real
//! `V0`. That is `cyrus_core::kklt_vacuum::verify_kklt_vacuum`, ~tens of
//! minutes, and is the toil this module automates: it is run ONLY on genomes
//! that already pass the scan and clear a fitness threshold (genuine vacua are
//! rare - on the order of one or two a day - so rarity, not throttling, is the
//! gate).
//!
//! KEY SUBTLETY (the reason this is a scan, not a single solve): the paper's
//! SS2 effective-cone criterion admits THOUSANDS of KKLT bases, and they are
//! NOT interchangeable. Each places the phase-1 selection point in a different
//! chamber; most chambers select small curves the cheap toric GV formulas
//! cannot cover, and whose general (HKTY-series) GV is intractable (grading
//! degree in the hundreds). Only McAllister-like chambers keep every selected
//! small curve low-degree/toric-coverable. So `deep_verify` walks the
//! admissible bases most-interior first and takes the first whose chamber is
//! fully covered by the cheap GV layers - reproducing the published `V0`
//! without any expensive general-GV computation.

use std::sync::OnceLock;

use cyrus_core::kklt_vacuum::{
    BranchSelection, OwnedDivisorBasis, PrimalGeom, PrimalIntersection, StabilizationInputs,
    VacuumConfig, VacuumVerdict, computed_primal_basis, verify_kklt_vacuum,
    verify_kklt_vacuum_auto_basis,
};
use cyrus_core::orientifold::OrientifoldModel;
use cyrus_core::types::f64::F64;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::{Finite, Pos};
use cyrus_core::{
    CurvePruningStrategy, Point, Polytope, compute_frst_heights, compute_glsm_and_linrels,
    compute_intersection_cytools, compute_linear_relations_no_origin, intersection_in_basis,
};

/// Mirror the runner's default small-curve volume cutoff (no `.dat` in the GA).
const SMALL_CURVE_CUTOFF: f64 = 1.0;
/// Path-following steps for the corrected F-flat solve (runner default for the
/// validated McAllister invocation).
const KKLT_STEPS: usize = 64;

/// Everything needed to deep-verify any flux candidate on a fixed primal
/// geometry, built once (lazily, on the first deep verification) and reused.
pub struct DeepVerifyContext {
    primal_geom: PrimalGeom,
    primal_intersection: PrimalIntersection,
    production_primal_basis: OwnedDivisorBasis,
    /// Admissible KKLT bases (paper SS2), most-interior first, each with its
    /// involution-derived `c_i`.
    candidates: Vec<(Vec<usize>, Vec<I64<Pos>>)>,
    /// Primal `h^{2,1}` (= flux length / mirror `h^{1,1}`), for the BBHL term.
    primal_h21: usize,
    /// Genome-independent coverage determination, computed once: `Ok(basis)` is
    /// the covered KKLT chamber, `Err` is the reason no chamber is coverable.
    ///
    /// The small-curve coverage that selects (or rejects) a chamber is driven by
    /// the phase-1 target `c_i` - a property of the geometry/orientifold, NOT
    /// the flux genome (which only scales `g_s`/`W0`/`e^{K0}` afterwards). So the
    /// expensive admissible-basis scan is the same for every genome on this
    /// polytope and is paid exactly once: without this, the GA (which
    /// re-evaluates the whole population every generation) would re-scan
    /// thousands of bases for every surviving elite and crawl.
    coverage: OnceLock<Result<Vec<usize>, String>>,
}

impl DeepVerifyContext {
    /// Build the primal-side deep-verify context: a primal FRST, its
    /// intersection numbers, the orientifold involution, and the admissible
    /// KKLT bases with their derived `c_i`.
    ///
    /// # Errors
    /// Returns a description of the first stage that fails.
    pub fn build(
        primal: &Polytope,
        primal_tri_points: &[Point],
        primal_h21: usize,
    ) -> Result<Self, String> {
        let origin_idx = primal_tri_points
            .iter()
            .position(|p| p.coords().iter().all(|&x| x == 0))
            .ok_or("primal origin not found")?;
        let (heights_raw, triangulation) = compute_frst_heights(primal_tri_points, origin_idx)
            .map_err(|e| format!("primal FRST: {e}"))?;
        let primal_geom = PrimalGeom {
            polytope: primal.clone(),
            heights: heights_raw
                .iter()
                .map(|&h| F64::<Finite>::new(h).ok_or("primal height not finite"))
                .collect::<Result<_, _>>()?,
            triangulation_points: primal_tri_points.to_vec(),
            triangulation,
        };

        let primal_intersection = build_primal_intersection(&primal_geom)?;
        let production_primal_basis = computed_primal_basis(&primal_intersection);

        let model = build_o7_model(&primal_geom, &primal_intersection, primal_h21)?;
        let bases = cyrus_core::orientifold::admissible_kklt_bases(primal, primal_tri_points, 0)
            .map_err(|e| format!("admissible KKLT bases: {e}"))?;
        if bases.is_empty() {
            return Err("no admissible KKLT basis (tau_* leaves the dual effective cone)".into());
        }
        let candidates: Vec<(Vec<usize>, Vec<I64<Pos>>)> = bases
            .into_iter()
            .map(|basis| {
                let c_i = c_i_for_basis(&model, &basis);
                (basis, c_i)
            })
            .collect();

        Ok(Self {
            primal_geom,
            primal_intersection,
            production_primal_basis,
            candidates,
            primal_h21,
            coverage: OnceLock::new(),
        })
    }

    /// Verify the primal KKLT vacuum for racetrack scalars produced by the
    /// mirror scan: scan the admissible bases (most-interior first) and return
    /// the first chamber the cheap GV layers fully cover.
    ///
    /// `ek0`, `g_s`, `w0` come straight from the candidate's
    /// `EvaluationResult` (`vacuum.ek0`, `racetrack.g_s`, `vacuum.w0`).
    ///
    /// # Errors
    /// Returns an error if no admissible basis verifies cleanly (the chain is
    /// fail-loud: no silent fallback to a proxy or an uncovered chamber).
    pub fn verify(
        &self,
        ek0: F64<Pos>,
        g_s: F64<Pos>,
        w0: F64<Pos>,
    ) -> Result<VacuumVerdict, String> {
        let config = Self::vacuum_config()?;

        // Determine the covered chamber ONCE (genome-independent - see the
        // `coverage` field): the first genome on this polytope pays the full
        // admissible-basis scan; every later genome reuses the verdict, so the
        // GA's per-generation re-evaluation never re-scans thousands of bases.
        let coverage = self.coverage.get_or_init(|| {
            let scan_inputs = self.base_inputs(ek0, g_s, w0);
            verify_kklt_vacuum_auto_basis(&scan_inputs, &self.candidates, &config)
                .map(|verdict| verdict.selected_kklt_basis)
        });
        let covered_basis = match coverage {
            Ok(basis) => basis,
            Err(reason) => return Err(reason.clone()),
        };

        // Verify in the cached covered chamber for THIS genome's scalars: the
        // chamber is genome-independent, only the resulting V0 is not.
        let candidate = self
            .candidates
            .iter()
            .find(|(basis, _)| basis == covered_basis)
            .ok_or("cached covered KKLT basis is not among the admissible candidates")?;
        let inputs = StabilizationInputs {
            geom: &self.primal_geom,
            intersection: &self.primal_intersection,
            production_primal_basis: &self.production_primal_basis,
            c_i: &candidate.1,
            kklt_basis: &candidate.0,
            g_s,
            w0,
            ek0,
            h21: self.primal_h21,
        };
        verify_kklt_vacuum(&inputs, &config)
    }

    /// Genome-independent stabilization config, shared by the one-time scan and
    /// the per-genome single-chamber verify.
    fn vacuum_config() -> Result<VacuumConfig, String> {
        Ok(VacuumConfig {
            kklt_steps: KKLT_STEPS,
            branch_candidates: 1,
            branch_seed: 42,
            branch_height_init: true,
            branch_selection: BranchSelection::FirstPositive,
            small_curve_cutoff: F64::<Pos>::new(SMALL_CURVE_CUTOFF)
                .ok_or("small-curve cutoff must be positive")?,
            small_curve_pruning: CurvePruningStrategy::PairDecomposable,
            // Cheap coverage only: the scan finds a chamber the toric +
            // minimal-degree layers fully cover, so general GV is never run.
            primal_gv_min_points: None,
            primal_gv_max_deg: None,
            max_missing_gv_impact: 0.0,
            missing_gv_abs_bound: 10,
        })
    }

    /// Inputs for the one-time auto-basis scan. The `c_i`/`kklt_basis` fields are
    /// ignored by the auto-basis entry point (each candidate supplies its own);
    /// they are present only to satisfy the struct.
    fn base_inputs(&self, ek0: F64<Pos>, g_s: F64<Pos>, w0: F64<Pos>) -> StabilizationInputs<'_> {
        StabilizationInputs {
            geom: &self.primal_geom,
            intersection: &self.primal_intersection,
            production_primal_basis: &self.production_primal_basis,
            c_i: &self.candidates[0].1,
            kklt_basis: &self.candidates[0].0,
            g_s,
            w0,
            ek0,
            h21: self.primal_h21,
        }
    }
}

/// Build the primal intersection data (GLSM/linrels/basis + intersection
/// tensors) from a primal geometry - the production intersection path.
fn build_primal_intersection(geom: &PrimalGeom) -> Result<PrimalIntersection, String> {
    let (glsm, linrels, basis) = compute_glsm_and_linrels(&geom.triangulation_points)
        .map_err(|e| format!("primal GLSM/linrels: {e}"))?;
    let points_i64: Vec<Vec<i64>> = geom
        .triangulation_points
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();
    let linrels_i64: Vec<Vec<i64>> = compute_linear_relations_no_origin(&points_i64)
        .iter()
        .map(|row| {
            row.iter()
                .map(|x| i64::try_from(x).map_err(|e| format!("primal linrel overflow: {e:?}")))
                .collect()
        })
        .collect::<Result<_, _>>()?;
    let kappa_full = compute_intersection_cytools(
        &geom.triangulation,
        &geom.triangulation_points,
        &linrels_i64,
    )
    .map_err(|e| format!("primal intersection numbers: {e}"))?;
    let kappa_basis = intersection_in_basis(&kappa_full, &basis);
    Ok(PrimalIntersection {
        glsm,
        linrels,
        basis,
        kappa_full,
        kappa_basis,
    })
}

/// The involution defining the O7 set: most O7-planes among the
/// pairwise-disjoint involutions (strongest racetrack hierarchy). Viability
/// against the GLSM basis is irrelevant - the KKLT basis is all-rigid by
/// construction; only gamma/components are read for `c_i`.
fn build_o7_model(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    primal_h21: usize,
) -> Result<OrientifoldModel, String> {
    let h11_extra =
        cyrus_core::divisor::batyrev_h11_extra_classes(&geom.polytope, &geom.triangulation_points)
            .map_err(|e| format!("primal Batyrev correction: {e}"))?;
    let h11 = intersection.basis.len() + h11_extra;
    let models = cyrus_core::orientifold::enumerate_involutions(
        &geom.polytope,
        &geom.triangulation_points,
        &intersection.kappa_full,
        &intersection.basis,
        h11,
        primal_h21,
    )
    .map_err(|e| format!("enumerate involutions: {e}"))?;
    models
        .into_iter()
        .filter(|m| m.o7_pairwise_disjoint && !m.o7_points.is_empty())
        .max_by_key(|m| m.o7_points.len())
        .ok_or_else(|| "no disjoint-O7 involution exists".to_string())
}

/// Derive `c_i` for a KKLT basis from a prebuilt involution model: 6 on
/// O7-parity divisors (so(8) gaugino condensation), else the irreducible-
/// component count.
fn c_i_for_basis(model: &OrientifoldModel, kklt_basis: &[usize]) -> Vec<I64<Pos>> {
    kklt_basis
        .iter()
        .map(|&idx| {
            let c = if model.gamma[idx] == 1 {
                6
            } else {
                i64::try_from(model.components[idx]).expect("component count fits i64")
            };
            I64::<Pos>::new(c).expect("c_i positive")
        })
        .collect()
}

/// Lazily-built deep-verify context attached to a geometry.
///
/// `OnceLock` makes the (expensive, one-time) primal build happen on the first
/// verification and be shared thereafter, safely across the parallel fitness
/// evaluations.
#[derive(Default)]
pub struct DeepVerifyCell {
    cell: OnceLock<Result<DeepVerifyContext, String>>,
}

impl DeepVerifyCell {
    /// Get (building on first use) the deep-verify context.
    pub fn get_or_build(
        &self,
        primal: &Polytope,
        primal_tri_points: &[Point],
        primal_h21: usize,
    ) -> Result<&DeepVerifyContext, String> {
        self.cell
            .get_or_init(|| DeepVerifyContext::build(primal, primal_tri_points, primal_h21))
            .as_ref()
            .map_err(Clone::clone)
    }
}
