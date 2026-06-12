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
    /// Simplices of the selected dual chamber (for verification emission).
    pub dual_simplices: Vec<Vec<usize>>,
    /// Dual GLSM charge matrix (for flux basis transforms).
    pub dual_glsm: Vec<Vec<malachite::Integer>>,
    /// Exact isotropic flux seeds (small-box scan at preparation time).
    /// Empty means PFV-barren within the scan budget: the quadratic form
    /// K.N(M)^{-1}.K is anisotropic (no integer solutions) for every
    /// sampled M - measured reality on some KS geometries, where random
    /// flux sampling provably cannot find vacua.
    pub pfv_seeds: Vec<crate::genome::Genome>,
    /// Automated D3 tadpole bound for this geometry's orientifold class:
    /// chi_f/4 = (h11 + h21)/2 + 1 (single-coordinate involutions with
    /// h11_- = h21_+ = 0, validated against all five published examples).
    pub q_d3: f64,
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
        Self::prepare_from_points_in_chamber(primal_points, gv_min_points, 0)
    }

    /// Prepare the geometry with the mirror (dual) side in the `chamber`-th
    /// dual FRST: chamber 0 is the default; chamber k > 0 is reached by k
    /// deterministic wall flips from the default. Different mirror chambers
    /// have different intersection numbers and Kahler cones, hence
    /// different PFV solution sets - the published scans search "in any of
    /// the LCS cones" of the dual polytope, and restricting to one chamber
    /// samples only a slice of each polytope's flux vacua.
    ///
    /// # Errors
    /// Returns an error string when any preparation stage fails, including
    /// when no k-th distinct valid chamber exists.
    pub fn prepare_from_points_in_chamber(
        primal_points: &[Vec<i64>],
        gv_min_points: u32,
        chamber: usize,
    ) -> Result<Self, String> {
        // Primal polytope and its true h11 (Batyrev-corrected), which is the
        // mirror's h21 for the BBHL term.
        let primal_point_objs: Vec<Point> = primal_points
            .iter()
            .map(|p| Point::new(p.clone()))
            .collect();
        // Canonicalize to the FULL lattice point set via double dualization
        // (compute_dual enumerates lattice points; the input may be
        // vertices-only, e.g. pools built from the raw KS database, and
        // points_not_interior_to_facets only filters what it is given).
        let primal = Polytope::from_vertices(primal_point_objs)
            .map_err(|e| format!("primal polytope: {e}"))?
            .compute_dual()
            .map_err(|e| format!("primal dual: {e}"))?
            .compute_dual()
            .map_err(|e| format!("primal double-dual: {e}"))?;
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
        let dual_tri = if chamber == 0 {
            dual_tri
        } else {
            flip_walk_chamber(&dual_points, dual_tri.simplices(), origin_idx, chamber)?
        };

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

        let dual_basis_len = dual_basis.len();
        let mirror_h11 = I32::<GTEOne>::new(
            i32::try_from(dual_basis.len()).map_err(|e| format!("mirror h11 overflow: {e}"))?,
        )
        .ok_or("mirror h11 must be >= 1")?;
        let mirror_h21 = I32::<NonNeg>::new(
            i32::try_from(primal_h11).map_err(|e| format!("mirror h21 overflow: {e}"))?,
        )
        .ok_or("mirror h21 must be >= 0")?;

        // Exact PFV seed scan (deterministic; integer arithmetic), keeping
        // only seeds that also clear the tadpole window and have a flat
        // direction strictly interior to the mirror Kahler cone (tested
        // against the enumerated GV curves) - the full set of cheap PFV
        // conditions. Viable seeds are RARE (the published examples came
        // from scanning thousands of geometries), so the M budget is large;
        // the integer-exact test costs ~ns per K.
        use rand::SeedableRng as _;
        let q_d3 = (primal_h11 + dual_basis.len()) as f64 / 2.0 + 1.0;
        let mut seed_rng = rand_chacha::ChaCha8Rng::seed_from_u64(0xC1B05);
        let keep = |genome: &crate::genome::Genome| -> bool {
            let q: f64 = -0.5
                * genome
                    .k
                    .iter()
                    .zip(genome.m.iter())
                    .map(|(&ki, &mi)| (ki * mi) as f64)
                    .sum::<f64>();
            if !(0.0..=q_d3).contains(&q) {
                return false;
            }
            // Flat direction p = N^{-1}K strictly interior: q.p > 0 for
            // every enumerated GV curve (mirror Kahler cone test).
            let m_typed: Vec<I64<Finite>> =
                genome.m.iter().map(|&x| I64::<Finite>::new(x)).collect();
            let k_typed: Vec<I64<Finite>> =
                genome.k.iter().map(|&x| I64::<Finite>::new(x)).collect();
            let n_mat = cyrus_core::flat_direction::compute_n_matrix(&kappa_basis, &m_typed);
            let Some(p) = cyrus_core::flat_direction::solve_linear_system_faer(&n_mat, &k_typed)
            else {
                return false;
            };
            gv.iter().all(|inv| {
                let a: f64 = inv
                    .curve
                    .iter()
                    .zip(p.iter())
                    .map(|(qi, pi)| qi.to_f64().get() * pi.get())
                    .sum();
                a > 1e-8
            })
        };
        let pfv_seeds = crate::pfv::find_isotropic_seeds(
            &kappa_basis,
            dual_basis.len(),
            &mut seed_rng,
            15,
            4,
            1500,
            64,
            keep,
        );

        Ok(Self {
            kappa_basis,
            mori,
            gv,
            h21_primal: dual_basis.len(),
            mirror_h11,
            mirror_h21,
            dual_basis,
            dual_glsm,
            dual_simplices: dual_tri.simplices().to_vec(),
            pfv_seeds,
            q_d3: (primal_h11 + dual_basis_len) as f64 / 2.0 + 1.0,
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

/// Walk `k` deterministic wall flips from the given dual FRST, keeping only
/// fine, star, regular neighbors. Errors when no k-th distinct chamber is
/// reachable within the attempt budget.
fn flip_walk_chamber(
    points: &[cyrus_core::Point],
    start: &[Vec<usize>],
    origin_idx: usize,
    k: usize,
) -> Result<cyrus_core::Triangulation, String> {
    use cyrus_core::triangulation::{
        flip_circuit_in_triangulation, secondary_cone_hyperplanes_native,
        triangulation_heights_from_secondary_cone,
    };
    use rand::Rng as _;
    use rand::SeedableRng as _;
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(0xD0A1 + k as u64);
    let mut simplices = start.to_vec();
    let mut done = 0usize;
    for _attempt in 0..(k * 24).max(24) {
        if done == k {
            break;
        }
        let tri = cyrus_core::Triangulation::new(simplices.clone());
        let circuits = secondary_cone_hyperplanes_native(points, &tri)
            .map_err(|e| format!("dual chamber walls: {e}"))?;
        if circuits.is_empty() {
            return Err("dual secondary cone has no walls".into());
        }
        let pick = &circuits[rng.gen_range(0..circuits.len())];
        let sparse: Vec<(usize, i64)> = pick
            .iter()
            .enumerate()
            .filter(|&(_, &c)| c != 0)
            .map(|(i, &c)| (i, c))
            .collect();
        let Ok(flip) = flip_circuit_in_triangulation(&simplices, &sparse) else {
            continue;
        };
        let neighbor = flip.simplices;
        let fine_star = neighbor.iter().all(|s| s.contains(&origin_idx)) && {
            let mut used = vec![false; points.len()];
            for s in &neighbor {
                for &i in s {
                    if i >= points.len() {
                        return Err("flip produced out-of-range index".into());
                    }
                    used[i] = true;
                }
            }
            used.iter().all(|&u| u)
        };
        if !fine_star {
            continue;
        }
        let neighbor_tri = cyrus_core::Triangulation::new(neighbor.clone());
        match triangulation_heights_from_secondary_cone(points, &neighbor_tri) {
            Ok(Some(_)) => {}
            _ => continue,
        }
        simplices = neighbor;
        done += 1;
    }
    if done < k {
        return Err(format!(
            "no distinct dual chamber {k} reachable ({done} flips found)"
        ));
    }
    Ok(cyrus_core::Triangulation::new(simplices))
}
