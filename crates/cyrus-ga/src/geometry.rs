//! One-time geometry preparation for a GA run.
//!
//! Given a primal polytope (CYTools-ordered points file), this prepares the
//! mirror-side data that `cyrus_core::evaluate_vacuum` consumes per flux
//! genome: dual polytope, dual FRST (Cyrus default chamber - the GA defines
//! its own conventions), dual intersection numbers in the computed divisor
//! basis, the Kahler-cone test rays, and the mirror GV invariants for the
//! racetrack. These are the same cyrus-core entry points the validated
//! McAllister runner drives, so every quantity here is backed by the
//! five-example reproduction and the landscape smoke test.

use cyrus_core::racetrack::GvInvariant;
use cyrus_core::types::f64::F64;
use cyrus_core::types::i32::I32;
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::{Finite, GTEOne, NonNeg};
use cyrus_core::{
    DivisorBasis, Intersection, MoriCone, Point, Polytope, compute_frst_heights,
    compute_glsm_and_linrels, compute_grading_vector, compute_intersection_cytools,
    compute_linear_relations_no_origin, compute_mori_cone_cap_rays, gv_divisor_basis_data,
};

/// Mirror GV computation default, matching the McAllister runner.
pub const DEFAULT_GV_MIN_POINTS: u32 = 20_000;

/// Prepared geometry shared by every fitness evaluation in a run.
pub struct GaGeometry {
    /// Dual (mirror) intersection numbers in the computed divisor basis.
    pub kappa_basis: Intersection,
    /// Kahler-cone membership test rays (projected dual Mori-cap rays).
    pub mori: MoriCone,
    /// Mirror GV invariants in basis coordinates for the racetrack.
    pub gv: Vec<GvInvariant>,
    /// Mirror h11 (= h21 of the primal CY) - the flux vector length.
    pub h21_primal: usize,
    /// Mirror Hodge numbers for the BBHL term of the proxy volume.
    pub mirror_h11: I32<GTEOne>,
    /// Mirror h21 (= true h11 of the primal CY, Batyrev-corrected).
    pub mirror_h21: I32<NonNeg>,
    /// Computed dual divisor basis indices (flux coordinate convention).
    pub dual_basis: Vec<usize>,
    /// Dual GLSM charge matrix (for flux basis transforms).
    pub dual_glsm: Vec<Vec<malachite::Integer>>,
}

fn read_points_csv(path: &std::path::Path) -> Result<Vec<Vec<i64>>, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("failed to read {}: {e}", path.display()))?;
    text.lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|cell| {
                    cell.trim()
                        .parse::<i64>()
                        .map_err(|e| format!("bad integer {cell:?} in {}: {e}", path.display()))
                })
                .collect()
        })
        .collect()
}

impl GaGeometry {
    /// Prepare from a McAllister-style data directory (`points.dat`).
    ///
    /// # Errors
    /// Returns a description of the first stage that fails; the GA must not
    /// start on a geometry it cannot evaluate faithfully.
    pub fn prepare_from_data_dir(
        data_dir: &std::path::Path,
        gv_min_points: u32,
    ) -> Result<Self, String> {
        let points = read_points_csv(&data_dir.join("points.dat"))?;
        Self::prepare_from_points(&points, gv_min_points)
    }

    /// Prepare from a full CYTools-ordered primal point list.
    ///
    /// # Errors
    /// Returns a description of the first stage that fails.
    #[allow(clippy::too_many_lines)] // linear stage-by-stage orchestration
    pub fn prepare_from_points(
        primal_points: &[Vec<i64>],
        gv_min_points: u32,
    ) -> Result<Self, String> {
        // Primal polytope and its true h11 (Batyrev-corrected), which is the
        // mirror's h21 for the BBHL term.
        let primal_point_objs: Vec<Point> = primal_points
            .iter()
            .map(|p| Point::new(p.clone()))
            .collect();
        let primal = Polytope::from_vertices(primal_point_objs)
            .map_err(|e| format!("primal polytope: {e}"))?;
        let primal_tri_points = primal
            .points_not_interior_to_facets()
            .map_err(|e| format!("primal points: {e}"))?;
        let (_, _, primal_basis) = compute_glsm_and_linrels(&primal_tri_points)
            .map_err(|e| format!("primal GLSM: {e}"))?;
        let primal_h11_extra =
            cyrus_core::divisor::batyrev_h11_extra_classes(&primal, &primal_tri_points)
                .map_err(|e| format!("primal Batyrev correction: {e}"))?;
        let primal_h11 = primal_basis.len() + primal_h11_extra;

        // Dual (mirror) side: default-chamber FRST and intersection numbers.
        let dual = primal
            .compute_dual()
            .map_err(|e| format!("dual polytope: {e}"))?;
        let dual_points = dual
            .points_not_interior_to_facets()
            .map_err(|e| format!("dual points: {e}"))?;
        let origin_idx = dual_points
            .iter()
            .position(|p| p.coords().iter().all(|&x| x == 0))
            .ok_or("dual origin not found")?;
        let (_, dual_tri) = compute_frst_heights(&dual_points, origin_idx)
            .map_err(|e| format!("dual FRST: {e}"))?;

        let (dual_glsm, linrels, dual_basis) =
            compute_glsm_and_linrels(&dual_points).map_err(|e| format!("dual GLSM: {e}"))?;
        let dual_points_i64: Vec<Vec<i64>> =
            dual_points.iter().map(|p| p.coords().to_vec()).collect();
        let linrels_reduced = compute_linear_relations_no_origin(&dual_points_i64);
        let linrels_i64: Vec<Vec<i64>> = linrels_reduced
            .iter()
            .map(|row| {
                row.iter()
                    .map(|x| i64::try_from(x).map_err(|e| format!("linrel overflow: {e:?}")))
                    .collect()
            })
            .collect::<Result<_, _>>()?;
        let kappa_full = compute_intersection_cytools(&dual_tri, &dual_points, &linrels_i64)
            .map_err(|e| format!("dual intersection numbers: {e}"))?;
        let kappa_basis = cyrus_core::intersection_in_basis(&kappa_full, &dual_basis);

        // Mori-cap rays projected to the basis: both the Kahler-cone test
        // and the GV computation use them, exactly like the runner.
        let cap_rays =
            compute_mori_cone_cap_rays(&dual_tri, &dual_points, &dual, false, false, None)
                .map_err(|e| format!("dual Mori cap: {e}"))?;
        let basis_data =
            gv_divisor_basis_data(&cap_rays, &linrels, DivisorBasis::Indices(&dual_basis))
                .map_err(|e| format!("GV basis data: {e}"))?;
        let grading = compute_grading_vector(&basis_data.mori_rays)
            .ok_or("no grading vector for dual Mori rays")?;

        let invariants = cyrus_core::compute_gv_invariants(
            &basis_data.mori_rays,
            &grading,
            &basis_data.q_matrix,
            &kappa_basis,
            Some(gv_min_points),
            None,
        )
        .map_err(|e| format!("mirror GV invariants: {e}"))?;
        let gv: Vec<GvInvariant> = invariants
            .into_iter()
            .map(|(curve, value)| {
                let value_f64 = value
                    .to_string()
                    .parse::<f64>()
                    .map_err(|e| format!("GV value parse: {e}"))?;
                if !value_f64.is_finite() {
                    return Err("GV value overflowed f64".to_string());
                }
                Ok(GvInvariant {
                    curve: curve
                        .into_iter()
                        .map(|x| I64::<Finite>::new(i64::from(x)))
                        .collect(),
                    value: F64::<Finite>::new(value_f64).ok_or("GV value not finite")?,
                })
            })
            .collect::<Result<_, String>>()?;

        let mori = MoriCone::new(
            basis_data
                .mori_rays
                .iter()
                .map(|ray| ray.iter().map(|&x| I64::<Finite>::new(x)).collect())
                .collect(),
        );

        let mirror_h11 = I32::<GTEOne>::new(
            i32::try_from(dual_basis.len()).map_err(|e| format!("mirror h11 overflow: {e}"))?,
        )
        .ok_or("mirror h11 must be >= 1")?;
        let mirror_h21 = I32::<NonNeg>::new(
            i32::try_from(primal_h11).map_err(|e| format!("mirror h21 overflow: {e}"))?,
        )
        .ok_or("mirror h21 must be >= 0")?;

        Ok(Self {
            kappa_basis,
            mori,
            gv,
            h21_primal: dual_basis.len(),
            mirror_h11,
            mirror_h21,
            dual_basis,
            dual_glsm,
        })
    }

    /// Transform a flux pair written in another index basis (e.g. the
    /// published McAllister `K_vec.dat`/`M_vec.dat`, which use the
    /// CYTools-2021 basis) into this geometry's computed basis, so known
    /// vacua can be injected into a run. M transforms with the basis-change
    /// matrix, K contravariantly (transpose), matching the validated
    /// McAllister runner.
    ///
    /// # Errors
    /// Returns an error if the bases are incompatible or non-unimodular.
    pub fn flux_pair_from_index_basis(
        &self,
        source_basis: &[usize],
        k: &[i64],
        m: &[i64],
    ) -> Result<(Vec<i64>, Vec<i64>), String> {
        if source_basis == self.dual_basis.as_slice() {
            return Ok((k.to_vec(), m.to_vec()));
        }
        let t_m = cyrus_core::basis_change_matrix(&self.dual_glsm, &self.dual_basis, source_basis)
            .map_err(|e| format!("M basis transform: {e}"))?;
        if !cyrus_core::is_unimodular(&t_m) {
            return Err("M basis transform is not unimodular".into());
        }
        let m_out = cyrus_core::apply_integer_basis_transform(&t_m, m, "M flux")
            .map_err(|e| format!("M transform apply: {e}"))?;
        let t_k = cyrus_core::basis_change_matrix(&self.dual_glsm, source_basis, &self.dual_basis)
            .map_err(|e| format!("K basis transform: {e}"))?;
        let k_out = cyrus_core::apply_integer_basis_transform_transpose(&t_k, k, "K flux")
            .map_err(|e| format!("K transform apply: {e}"))?;
        Ok((k_out, m_out))
    }
}
