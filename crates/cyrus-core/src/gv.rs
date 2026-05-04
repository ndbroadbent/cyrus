//! Gopakumar-Vafa invariant computation utilities.
//!
//! Implements the CYTools mori_cone_cap algorithm and grading vector selection,
//! then wires up the `cygv` crate for actual GV computation.
//!
//! Reference: CYTools `calabiyau.py` (mori_cone_cap) and `cone.py` (grading vector).

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;
use std::hash::{Hash, Hasher};
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use good_lp::{Expression, ProblemVariables, Solution, SolverModel, default_solver, variable};
use malachite::Integer;
use malachite::Rational;
use malachite::num::arithmetic::traits::Abs;
use nalgebra::{DMatrix, DVector, RowDVector};

use crate::cone::Cone;
use crate::error::{Error, Result};
use crate::geometry::ConvexHull;
use crate::integer_math::matrix_rank;
use crate::integer_math::{gcd_integer, integer_kernel};
use crate::intersection::Intersection;
use crate::lattice::Point;
use crate::polytope::Polytope;
use crate::triangulation::Triangulation;
use crate::types::{F64, Finite, I64, Pos};

const GRADING_CACHE_VERSION: &str = "grading-vector-cytools-lp-v1";
const LATTICE_CACHE_VERSION: &str = "lattice-points-v2";

/// Compute the Mori cone cap generators (rays) using the CYTools algorithm.
///
/// Returns a matrix where each row is a generator (ray) expressed in the
/// divisor/curve coordinates corresponding to `points`.
///
/// This follows `CalabiYau.mori_cone_cap` in CYTools.
pub fn compute_mori_cone_cap_rays(
    tri: &Triangulation,
    points: &[Point],
    polytope: &Polytope,
    in_basis: bool,
    exclude_origin: bool,
    basis: Option<&[usize]>,
) -> Result<Vec<Vec<i64>>> {
    if points.is_empty() {
        return Err(Error::InvalidInput("No points provided".into()));
    }
    if polytope.dim() != 4 {
        return Err(Error::InvalidInput(
            "mori_cone_cap is only implemented for 4D polytopes".into(),
        ));
    }

    let origin_idx = points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .ok_or_else(|| Error::InvalidInput("Origin not found in points".into()))?;

    // pts_ext: append 1 to each point (for affine embedding)
    let pts_ext: Vec<Vec<i64>> = points
        .iter()
        .map(|p| {
            let mut v = p.coords().to_vec();
            v.push(1);
            v
        })
        .collect();

    let (facets, twofaces) = compute_faces_4d(points, polytope)?;

    // Collect 2D simplices and circuits
    let mut mori_cap_rays: HashSet<Vec<(usize, i64)>> = HashSet::new();
    let mut simp_2d_all: HashSet<Vec<usize>> = HashSet::new();

    for f in &twofaces {
        if f.len() < 4 {
            let mut f_sorted = f.clone();
            f_sorted.sort_unstable();
            simp_2d_all.insert(f_sorted);
            continue;
        }

        let face_pts: HashSet<usize> = f.iter().copied().collect();
        let mut simp_2d: HashSet<Vec<usize>> = HashSet::new();

        for simplex in tri.simplices() {
            let inter: Vec<usize> = simplex
                .iter()
                .filter(|idx| face_pts.contains(idx))
                .copied()
                .collect();
            if inter.len() == 3 {
                let mut inter_sorted = inter;
                inter_sorted.sort_unstable();
                simp_2d.insert(inter_sorted.clone());
                simp_2d_all.insert(inter_sorted);
            }
        }

        let simps: Vec<Vec<usize>> = simp_2d.into_iter().collect();
        for i in 0..simps.len() {
            for j in i..simps.len() {
                let s1: HashSet<usize> = simps[i].iter().copied().collect();
                let s2: HashSet<usize> = simps[j].iter().copied().collect();
                let comm: Vec<usize> = s1.intersection(&s2).copied().collect();
                if comm.len() == 2 {
                    let diff: Vec<usize> = s1.symmetric_difference(&s2).copied().collect();
                    if diff.len() != 2 {
                        continue;
                    }

                    if let Some(v) = nullspace_vector(&pts_ext, &diff, &comm, false) {
                        let full_v = build_full_v(&diff, &comm, &v);
                        mori_cap_rays.insert(full_v);
                    }
                }
            }
        }
    }

    // Origin circuits
    for s2d in &simp_2d_all {
        let s2d_set: HashSet<usize> = s2d.iter().copied().collect();
        let mut f1: Option<Vec<usize>> = None;
        let mut f2: Option<Vec<usize>> = None;
        for facet in &facets {
            let facet_set: HashSet<usize> = facet.iter().copied().collect();
            if s2d_set.is_subset(&facet_set) {
                if f1.is_none() {
                    f1 = Some(facet.clone());
                } else {
                    f2 = Some(facet.clone());
                    break;
                }
            }
        }
        let (Some(f1), Some(f2)) = (f1, f2) else {
            continue;
        };

        let f1_set: HashSet<usize> = f1.iter().copied().collect();
        let f2_set: HashSet<usize> = f2.iter().copied().collect();
        let pts_f1: Vec<usize> = f1_set.difference(&f2_set).copied().collect();
        let pts_f2: Vec<usize> = f2_set.difference(&f1_set).copied().collect();

        for p1 in &pts_f1 {
            for p2 in &pts_f2 {
                let diff = vec![*p1, *p2];
                let mut comm = s2d.clone();
                comm.push(origin_idx);
                let Some(v) = nullspace_vector(&pts_ext, &diff, &comm, true) else {
                    continue;
                };
                let full_v = build_full_v(&diff, &comm, &v);

                // Skip if origin coefficient is zero or positive (CYTools behavior).
                let origin_coeff = full_v
                    .iter()
                    .find(|(idx, _)| *idx == origin_idx)
                    .map_or(0, |(_, coeff)| *coeff);
                if origin_coeff == 0 {
                    continue;
                }
                if origin_coeff > 0 {
                    continue;
                }
                mori_cap_rays.insert(full_v);
            }
        }
    }

    let mut rays: Vec<Vec<i64>> = Vec::new();
    let n_pts = pts_ext.len();
    for r in mori_cap_rays {
        let mut row = vec![0i64; n_pts];
        for (idx, coeff) in r {
            row[idx] = coeff;
        }
        rays.push(row);
    }

    if exclude_origin && !in_basis {
        rays = rays.into_iter().map(|r| r[1..].to_vec()).collect();
    } else if in_basis {
        let basis = basis.ok_or_else(|| {
            Error::InvalidInput("basis indices required for in_basis=true".into())
        })?;
        rays = rays
            .into_iter()
            .map(|r| basis.iter().map(|&idx| r[idx]).collect())
            .collect();
    }

    // Normalize by GCD and drop zero rows (CYTools Cone normalization).
    let mut normalized: Vec<Vec<i64>> = Vec::with_capacity(rays.len());
    for row in rays {
        let mut g = 0i64;
        for &x in &row {
            g = gcd_i64(g, x.abs());
        }
        if g == 0 {
            continue;
        }
        normalized.push(row.into_iter().map(|x| x / g).collect());
    }

    // Deduplicate rows after normalization.
    let mut uniq: HashSet<Vec<i64>> = HashSet::new();
    let mut deduped = Vec::with_capacity(normalized.len());
    for row in normalized {
        if uniq.insert(row.clone()) {
            deduped.push(row);
        }
    }

    deduped.sort();
    Ok(deduped)
}

/// A toric curve candidate with its volume at a specific point in Kähler moduli space.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ToricCurveCandidate {
    /// Curve class in ambient divisor-intersection coordinates.
    pub class: Vec<i64>,
    /// Positive curve volume at the selected Kähler point.
    pub volume: F64<Pos>,
}

/// A toric curve class with its genus-zero GV invariant.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ToricCurveGvInvariant {
    /// Curve class in ambient divisor-intersection coordinates.
    pub class: Vec<i64>,
    /// Genus-zero GV invariant for this toric curve class.
    pub gv: Integer,
}

/// Compute the volume of an ambient curve class from Kähler parameters in a divisor basis.
pub fn curve_volume_in_divisor_basis(
    curve: &[i64],
    basis: &[usize],
    kahler_parameters: &[F64<Finite>],
) -> Result<F64<Finite>> {
    if basis.len() != kahler_parameters.len() {
        return Err(Error::InvalidInput(format!(
            "basis length {} does not match Kähler parameter length {}",
            basis.len(),
            kahler_parameters.len()
        )));
    }
    let mut volume = F64::<Finite>::ZERO;
    for (&idx, &ti) in basis.iter().zip(kahler_parameters.iter()) {
        let Some(&coefficient) = curve.get(idx) else {
            return Err(Error::InvalidInput(format!(
                "basis index {idx} is out of bounds for curve dimension {}",
                curve.len()
            )));
        };
        volume = volume + I64::<Finite>::new(coefficient).to_f64() * ti;
    }
    Ok(volume)
}

/// Select toric curve candidates with positive volume below a cutoff.
pub fn subcutoff_toric_curve_candidates(
    rays: &[Vec<i64>],
    basis: &[usize],
    kahler_parameters: &[F64<Finite>],
    cutoff: F64<Pos>,
) -> Result<Vec<ToricCurveCandidate>> {
    let mut out = Vec::new();
    for ray in rays {
        let volume = curve_volume_in_divisor_basis(ray, basis, kahler_parameters)?;
        if let Some(volume) = volume.try_to_pos()
            && volume < cutoff
        {
            out.push(ToricCurveCandidate {
                class: ray.clone(),
                volume,
            });
        }
    }
    Ok(out)
}

/// Find a decomposition of a curve as a sum of two selected toric curve candidates.
pub fn find_pair_decomposition(
    target: &ToricCurveCandidate,
    candidates: &[ToricCurveCandidate],
) -> Result<Option<(Vec<i64>, Vec<i64>)>> {
    let dim = target.class.len();
    let selected_set: HashSet<Vec<i64>> = candidates
        .iter()
        .map(|candidate| {
            if candidate.class.len() != dim {
                return Err(Error::InvalidInput(
                    "candidate curve dimensions are inconsistent".into(),
                ));
            }
            Ok(candidate.class.clone())
        })
        .collect::<Result<_>>()?;

    Ok(find_pair_decomposition_with_set(
        target,
        candidates,
        &selected_set,
    ))
}

/// Remove curve candidates that are sums of two other selected candidates.
///
/// The McAllister paper removes small toric curves that can be written as sums
/// of others before computing GV invariants. This helper implements the
/// currently verified pair-sum part of that pruning and matches the published
/// 4-214-647 small-curve checkpoint. It is not a full Hilbert-basis or
/// arbitrary-length semigroup reduction.
pub fn remove_pair_decomposable_curve_candidates(
    candidates: &[ToricCurveCandidate],
) -> Result<Vec<ToricCurveCandidate>> {
    if candidates.is_empty() {
        return Ok(Vec::new());
    }
    let dim = candidates[0].class.len();
    let mut selected_set: HashSet<Vec<i64>> = HashSet::with_capacity(candidates.len());
    for candidate in candidates {
        if candidate.class.len() != dim {
            return Err(Error::InvalidInput(
                "candidate curve dimensions are inconsistent".into(),
            ));
        }
        selected_set.insert(candidate.class.clone());
    }

    Ok(candidates
        .iter()
        .filter(|candidate| {
            find_pair_decomposition_with_set(candidate, candidates, &selected_set).is_none()
        })
        .cloned()
        .collect())
}

/// Exact bounded semigroup-decomposition diagnostic for toric curve candidates.
///
/// This is intentionally diagnostic-only. It can prove that a curve is a sum
/// of up to three other selected candidates, which is useful for investigating
/// the remaining McAllister "sums of others" pruning gap. It does not replace
/// a full Hilbert-basis computation.
pub struct BoundedCurveDecompositionIndex<'a> {
    candidates: &'a [ToricCurveCandidate],
    class_set: HashSet<Vec<i64>>,
    pair_sums: HashMap<Vec<i64>, Vec<(usize, usize)>>,
    dim: usize,
}

impl<'a> BoundedCurveDecompositionIndex<'a> {
    /// Build a bounded decomposition index over selected curve candidates.
    pub fn new(candidates: &'a [ToricCurveCandidate]) -> Result<Self> {
        let Some(first) = candidates.first() else {
            return Ok(Self {
                candidates,
                class_set: HashSet::new(),
                pair_sums: HashMap::new(),
                dim: 0,
            });
        };
        let dim = first.class.len();
        let mut class_set = HashSet::with_capacity(candidates.len());
        for candidate in candidates {
            if candidate.class.len() != dim {
                return Err(Error::InvalidInput(
                    "candidate curve dimensions are inconsistent".into(),
                ));
            }
            class_set.insert(candidate.class.clone());
        }

        let mut pair_sums: HashMap<Vec<i64>, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..candidates.len() {
            for j in i..candidates.len() {
                let sum = add_curve(&candidates[i].class, &candidates[j].class);
                pair_sums.entry(sum).or_default().push((i, j));
            }
        }

        Ok(Self {
            candidates,
            class_set,
            pair_sums,
            dim,
        })
    }

    /// Find a decomposition of `target` as a sum of up to `max_terms` indexed candidates.
    ///
    /// Currently supports `max_terms <= 3`. Terms may repeat, matching
    /// semigroup membership rather than a set partition.
    pub fn find_decomposition(
        &self,
        target: &ToricCurveCandidate,
        max_terms: usize,
    ) -> Result<Option<Vec<Vec<i64>>>> {
        if max_terms < 2 || self.candidates.is_empty() {
            return Ok(None);
        }
        if target.class.len() != self.dim {
            return Err(Error::InvalidInput(
                "target curve dimension does not match decomposition index".into(),
            ));
        }
        if max_terms > 3 {
            return Err(Error::InvalidInput(
                "bounded curve decomposition currently supports max_terms <= 3".into(),
            ));
        }

        if let Some((first, second)) =
            find_pair_decomposition_with_set(target, self.candidates, &self.class_set)
        {
            return Ok(Some(vec![first, second]));
        }
        if max_terms < 3 {
            return Ok(None);
        }

        for summand in self
            .candidates
            .iter()
            .filter(|candidate| candidate.volume < target.volume)
        {
            let remainder = subtract_curve(&target.class, &summand.class);
            let Some(pairs) = self.pair_sums.get(&remainder) else {
                continue;
            };
            for &(i, j) in pairs {
                let first = &self.candidates[i];
                let second = &self.candidates[j];
                if first.volume < target.volume && second.volume < target.volume {
                    let mut decomposition = vec![
                        first.class.clone(),
                        second.class.clone(),
                        summand.class.clone(),
                    ];
                    decomposition.sort();
                    return Ok(Some(decomposition));
                }
            }
        }

        Ok(None)
    }
}

/// Compute genus-zero GV invariants for simple toric curve classes.
///
/// This implements the analytic formulas from Demirtas-Kim-McAllister-Moritz-
/// Rios-Tascon, "Computational Mirror Symmetry", section 6. The supported
/// classes are complete intersections `C_IJ = D_I cap D_J` where the edge
/// `(I,J)` lies in a two-face of the reflexive polytope, plus origin circuits
/// that realize the isolated resolved-conifold normal bundle
/// `O(-1) + O(-1)`. The returned classes use the same ambient
/// divisor-intersection coordinates as [`compute_mori_cone_cap_rays`] with
/// `in_basis=false`.
pub fn compute_toric_two_face_curve_gv_invariants(
    tri: &Triangulation,
    points: &[Point],
    polytope: &Polytope,
) -> Result<Vec<ToricCurveGvInvariant>> {
    if points.is_empty() {
        return Err(Error::InvalidInput("No points provided".into()));
    }
    if polytope.dim() != 4 {
        return Err(Error::InvalidInput(
            "toric two-face GV formulas are only implemented for 4D polytopes".into(),
        ));
    }

    let pts_ext: Vec<Vec<i64>> = points
        .iter()
        .map(|p| {
            let mut v = p.coords().to_vec();
            v.push(1);
            v
        })
        .collect();
    let face_data = compute_two_face_data_4d(points, polytope)?;
    let point_face_dims = classify_primal_point_face_dimensions(points, polytope)?;
    let one_face_genera = compute_primal_one_face_genera(points, polytope, &point_face_dims)?;
    let (facets, _) = compute_faces_4d(points, polytope)?;
    let origin_idx = points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .ok_or_else(|| Error::InvalidInput("Origin not found in points".into()))?;

    let mut gv_by_class: HashMap<Vec<i64>, Integer> = HashMap::new();
    let mut simp_2d_all: HashSet<Vec<usize>> = HashSet::new();
    for face in &face_data {
        if face.points.len() < 4 {
            continue;
        }
        let face_pts: HashSet<usize> = face.points.iter().copied().collect();
        let mut simp_2d: HashSet<Vec<usize>> = HashSet::new();

        for simplex in tri.simplices() {
            let inter: Vec<usize> = simplex
                .iter()
                .filter(|idx| face_pts.contains(idx))
                .copied()
                .collect();
            if inter.len() == 3 {
                let mut inter_sorted = inter;
                inter_sorted.sort_unstable();
                simp_2d.insert(inter_sorted.clone());
                simp_2d_all.insert(inter_sorted);
            }
        }

        let simps: Vec<Vec<usize>> = simp_2d.into_iter().collect();
        for i in 0..simps.len() {
            for j in i..simps.len() {
                let s1: HashSet<usize> = simps[i].iter().copied().collect();
                let s2: HashSet<usize> = simps[j].iter().copied().collect();
                let comm: Vec<usize> = s1.intersection(&s2).copied().collect();
                if comm.len() != 2 {
                    continue;
                }
                let diff: Vec<usize> = s1.symmetric_difference(&s2).copied().collect();
                if diff.len() != 2 {
                    continue;
                }

                let Some(v) = nullspace_vector(&pts_ext, &diff, &comm, false) else {
                    continue;
                };
                let full_v = build_full_v(&diff, &comm, &v);
                let class = normalized_row_from_sparse_relation(points.len(), full_v);
                if class.iter().all(|&x| x == 0) {
                    continue;
                }
                let gv = toric_two_face_curve_gv(
                    &class,
                    &comm,
                    face.genus,
                    &point_face_dims,
                    &one_face_genera,
                )?;
                insert_toric_gv(&mut gv_by_class, class, gv)?;
            }
        }
    }

    for class in compute_origin_circuit_curve_classes(&pts_ext, &facets, &simp_2d_all, origin_idx) {
        if let Some(gv) = resolved_conifold_origin_circuit_gv(&class, origin_idx) {
            insert_toric_gv(&mut gv_by_class, class, gv)?;
        }
    }

    let mut out: Vec<ToricCurveGvInvariant> = gv_by_class
        .into_iter()
        .map(|(class, gv)| ToricCurveGvInvariant { class, gv })
        .collect();
    out.sort_by(|a, b| a.class.cmp(&b.class));
    Ok(out)
}

fn insert_toric_gv(
    gv_by_class: &mut HashMap<Vec<i64>, Integer>,
    class: Vec<i64>,
    gv: Integer,
) -> Result<()> {
    match gv_by_class.entry(class) {
        std::collections::hash_map::Entry::Occupied(existing) => {
            if existing.get() != &gv {
                return Err(Error::InvalidInput(format!(
                    "conflicting toric GV values for duplicate curve class: {} vs {gv}",
                    existing.get()
                )));
            }
        }
        std::collections::hash_map::Entry::Vacant(slot) => {
            slot.insert(gv);
        }
    }
    Ok(())
}

fn find_pair_decomposition_with_set(
    target: &ToricCurveCandidate,
    candidates: &[ToricCurveCandidate],
    selected_set: &HashSet<Vec<i64>>,
) -> Option<(Vec<i64>, Vec<i64>)> {
    for summand in candidates {
        if summand.volume >= target.volume {
            continue;
        }
        let remainder = subtract_curve(&target.class, &summand.class);
        if selected_set.contains(&remainder) {
            return Some((summand.class.clone(), remainder));
        }
    }
    None
}

fn subtract_curve(lhs: &[i64], rhs: &[i64]) -> Vec<i64> {
    debug_assert_eq!(lhs.len(), rhs.len());
    lhs.iter().zip(rhs.iter()).map(|(&a, &b)| a - b).collect()
}

fn add_curve(lhs: &[i64], rhs: &[i64]) -> Vec<i64> {
    debug_assert_eq!(lhs.len(), rhs.len());
    lhs.iter().zip(rhs.iter()).map(|(&a, &b)| a + b).collect()
}

fn toric_two_face_curve_gv(
    class: &[i64],
    edge: &[usize],
    two_face_genus: usize,
    point_face_dims: &[usize],
    one_face_genera: &[Option<usize>],
) -> Result<Integer> {
    if edge.len() != 2 {
        return Err(Error::InvalidInput(
            "two-face toric curve edge must have two endpoints".into(),
        ));
    }
    let i = edge[0];
    let j = edge[1];
    let Some(&di) = class.get(i) else {
        return Err(Error::InvalidInput(
            "edge endpoint out of class bounds".into(),
        ));
    };
    let Some(&dj) = class.get(j) else {
        return Err(Error::InvalidInput(
            "edge endpoint out of class bounds".into(),
        ));
    };
    let (m, n_idx, n) = if di >= dj { (di, j, dj) } else { (dj, i, di) };
    if m + n != -2 {
        return Err(Error::InvalidInput(format!(
            "two-face curve normal degrees must sum to -2, got {m}+{n} for edge ({i},{j})"
        )));
    }

    let i_is_two_face = point_face_dims.get(i).copied() == Some(2);
    let j_is_two_face = point_face_dims.get(j).copied() == Some(2);

    if m == 0 && point_face_dims.get(n_idx).copied() == Some(1) {
        let Some(g_prime) = one_face_genera.get(n_idx).copied().flatten() else {
            return Err(Error::InvalidInput(format!(
                "missing dual two-face genus for one-face divisor {n_idx}"
            )));
        };
        return Ok(Integer::from(2i64 * g_prime as i64 - 2));
    }

    if i_is_two_face || j_is_two_face || (m, n) != (-1, -1) {
        let magnitude = m + 2;
        if magnitude < 0 {
            return Err(Error::InvalidInput(format!(
                "unsupported toric curve normal degree m={m} for edge ({i},{j})"
            )));
        }
        let sign = if (m + 1).rem_euclid(2) == 0 { 1 } else { -1 };
        return Ok(Integer::from(sign * magnitude));
    }

    Ok(Integer::from(two_face_genus as i64 + 1))
}

fn compute_origin_circuit_curve_classes(
    pts_ext: &[Vec<i64>],
    facets: &[Vec<usize>],
    simp_2d_all: &HashSet<Vec<usize>>,
    origin_idx: usize,
) -> Vec<Vec<i64>> {
    let mut out = Vec::new();
    for s2d in simp_2d_all {
        let s2d_set: HashSet<usize> = s2d.iter().copied().collect();
        let mut f1: Option<Vec<usize>> = None;
        let mut f2: Option<Vec<usize>> = None;
        for facet in facets {
            let facet_set: HashSet<usize> = facet.iter().copied().collect();
            if s2d_set.is_subset(&facet_set) {
                if f1.is_none() {
                    f1 = Some(facet.clone());
                } else {
                    f2 = Some(facet.clone());
                    break;
                }
            }
        }
        let (Some(f1), Some(f2)) = (f1, f2) else {
            continue;
        };

        let f1_set: HashSet<usize> = f1.iter().copied().collect();
        let f2_set: HashSet<usize> = f2.iter().copied().collect();
        let pts_f1: Vec<usize> = f1_set.difference(&f2_set).copied().collect();
        let pts_f2: Vec<usize> = f2_set.difference(&f1_set).copied().collect();

        for p1 in &pts_f1 {
            for p2 in &pts_f2 {
                let diff = vec![*p1, *p2];
                let mut comm = s2d.clone();
                comm.push(origin_idx);
                let Some(v) = nullspace_vector(pts_ext, &diff, &comm, true) else {
                    continue;
                };
                let full_v = build_full_v(&diff, &comm, &v);
                let origin_coeff = full_v
                    .iter()
                    .find(|(idx, _)| *idx == origin_idx)
                    .map_or(0, |(_, coeff)| *coeff);
                if origin_coeff >= 0 {
                    continue;
                }
                out.push(normalized_row_from_sparse_relation(pts_ext.len(), full_v));
            }
        }
    }
    out
}

fn resolved_conifold_origin_circuit_gv(class: &[i64], origin_idx: usize) -> Option<Integer> {
    if class.get(origin_idx).copied()? != -1 {
        return None;
    }

    let mut non_origin_neg_ones = 0usize;
    let mut pos_ones = 0usize;
    for (idx, &coeff) in class.iter().enumerate() {
        match (idx == origin_idx, coeff) {
            (_, 0) => {}
            (true, -1) => {}
            (false, -1) => non_origin_neg_ones += 1,
            (false, 1) => pos_ones += 1,
            _ => return None,
        }
    }

    if non_origin_neg_ones == 2 && pos_ones == 3 {
        // An isolated resolved conifold curve has normal bundle
        // O(-1) + O(-1), so section 6 gives GV^0 = 1.
        Some(Integer::from(1))
    } else {
        None
    }
}

/// Compute a grading vector for the Mori cone cap.
///
/// The grading vector must be strictly interior to the dual cone.
pub fn compute_grading_vector(rays: &[Vec<i64>]) -> Option<Vec<i64>> {
    if rays.is_empty() {
        return None;
    }
    let (cache_dir, cache_path) = grading_cache_paths(rays);
    if let Some(cached) = load_grading_cache(&cache_path, rays) {
        eprintln!("[DEBUG] grading vector: cache hit");
        return Some(cached);
    }

    let interior = cytools_grading_lp_solution(rays).or_else(|| {
        let (zero_rays, opposite_pairs) = analyze_rays(rays);
        eprintln!(
            "[WARN] grading vector search failed: rays={}, zero_rays={}, opposite_pairs={}",
            rays.len(),
            zero_rays,
            opposite_pairs
        );
        None
    })?;

    // Match CYTools Cone.find_interior_point(integral=True): scale the LP
    // solution by successive positive integers and round until every defining
    // ray has strictly positive pairing with the candidate.
    for scale in 1..1000 {
        let candidate: Vec<i64> = interior
            .iter()
            .map(|x| (x * f64::from(scale)).round() as i64)
            .collect();
        if candidate.iter().all(|&x| x == 0) {
            continue;
        }
        if is_strictly_dual(rays, &candidate) {
            write_grading_cache(&cache_dir, &cache_path, &candidate);
            return Some(candidate);
        }
    }

    None
}

fn cytools_grading_lp_solution(rays: &[Vec<i64>]) -> Option<Vec<f64>> {
    let dim = rays.first()?.len();
    if dim == 0 || rays.iter().any(|row| row.len() != dim) {
        return None;
    }

    let mut vars = ProblemVariables::new();
    let bound = 1.0e9;
    let x: Vec<_> = (0..dim)
        .map(|_| vars.add(variable().min(-bound).max(bound)))
        .collect();

    // CYTools feasibility(..., backend="glop") minimizes the average
    // hyperplane normal while imposing H x >= 1.
    let mut objective = Expression::from(0.0);
    for i in 0..dim {
        let coeff: f64 = rays.iter().map(|row| row[i] as f64).sum::<f64>() / rays.len() as f64;
        objective.add_mul(coeff, x[i]);
    }

    let mut model = vars.minimise(objective).using(default_solver);
    for row in rays {
        let mut expr = Expression::from(0.0);
        for (i, &coeff) in row.iter().enumerate() {
            if coeff != 0 {
                expr.add_mul(coeff as f64, x[i]);
            }
        }
        model = model.with(expr.geq(1.0));
    }

    let solution = match model.solve() {
        Ok(solution) => solution,
        Err(err) => {
            eprintln!("[WARN] grading vector LP failed: {err}");
            return None;
        }
    };
    let candidate: Vec<f64> = x.iter().map(|var| solution.value(*var)).collect();
    if candidate.iter().all(|v| v.is_finite()) {
        Some(candidate)
    } else {
        None
    }
}

fn grading_cache_paths(rays: &[Vec<i64>]) -> (PathBuf, PathBuf) {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    GRADING_CACHE_VERSION.hash(&mut hasher);
    let mut key_rays: Vec<&[i64]> = rays.iter().map(Vec::as_slice).collect();
    key_rays.sort();
    for row in &key_rays {
        row.hash(&mut hasher);
    }
    let key = hasher.finish();
    let cache_dir = env::var("CYRUS_CACHE_DIR")
        .map_or_else(|_| PathBuf::from("target/cyrus-cache"), PathBuf::from);
    let cache_path = cache_dir.join(format!("grading_vector_{key:x}.json"));
    (cache_dir, cache_path)
}

fn load_grading_cache(path: &Path, rays: &[Vec<i64>]) -> Option<Vec<i64>> {
    if !path.exists() {
        return None;
    }
    let data = match fs::read_to_string(path) {
        Ok(data) => data,
        Err(err) => {
            eprintln!(
                "[WARN] failed to read grading cache {}: {err}",
                path.display()
            );
            return None;
        }
    };
    let candidate: Vec<i64> = match serde_json::from_str(&data) {
        Ok(candidate) => candidate,
        Err(err) => {
            eprintln!(
                "[WARN] failed to parse grading cache {}: {err}",
                path.display()
            );
            return None;
        }
    };
    let dim = rays.first()?.len();
    if candidate.len() != dim || !is_strictly_dual(rays, &candidate) {
        eprintln!("[WARN] ignoring invalid grading cache {}", path.display());
        return None;
    }
    Some(candidate)
}

fn write_grading_cache(cache_dir: &Path, path: &Path, candidate: &[i64]) {
    if let Err(err) = fs::create_dir_all(cache_dir) {
        eprintln!(
            "[WARN] failed to create grading cache dir {}: {err}",
            cache_dir.display()
        );
        return;
    }
    match fs::File::create(path) {
        Ok(file) => {
            let writer = BufWriter::new(file);
            if let Err(err) = serde_json::to_writer(writer, candidate) {
                eprintln!(
                    "[WARN] failed to serialize grading cache {}: {err}",
                    path.display()
                );
            }
        }
        Err(err) => {
            eprintln!(
                "[WARN] failed to create grading cache {}: {err}",
                path.display()
            );
        }
    }
}

struct LatticeCacheControls {
    max_coord: i64,
    deg_window: i64,
    cacheable: bool,
}

impl LatticeCacheControls {
    fn from_env(default_max_coord: i64, default_deg_window: i64) -> Self {
        let max_coord = env::var("CYRUS_LATTICE_MAX_COORD")
            .ok()
            .and_then(|value| value.parse::<i64>().ok())
            .unwrap_or(default_max_coord);
        let deg_window = env::var("CYRUS_LATTICE_DEG_WINDOW")
            .ok()
            .and_then(|value| value.parse::<i64>().ok())
            .unwrap_or(default_deg_window);

        let strict = env::var("CYRUS_LATTICE_STRICT")
            .map(|value| value != "0")
            .unwrap_or(true);
        let has_bounded_search = env::var("CYRUS_LATTICE_MAX_TIME_SEC").is_ok()
            || env::var("CYRUS_LATTICE_MAX_SOLUTIONS").is_ok()
            || env::var("CYRUS_LATTICE_MAX_DEG").is_ok();

        Self {
            max_coord,
            deg_window,
            cacheable: strict && !has_bounded_search,
        }
    }
}

fn analyze_rays(rays: &[Vec<i64>]) -> (usize, usize) {
    let mut zero_rays = 0usize;
    let mut opposite_pairs = 0usize;
    let mut seen: HashSet<Vec<i64>> = HashSet::new();

    for r in rays {
        if r.iter().all(|&x| x == 0) {
            zero_rays += 1;
            continue;
        }
        let norm = normalize_ray(r);
        let neg: Vec<i64> = norm.iter().map(|&x| -x).collect();
        if seen.contains(&neg) {
            opposite_pairs += 1;
        }
        seen.insert(norm);
    }

    (zero_rays, opposite_pairs)
}

fn normalize_ray(r: &[i64]) -> Vec<i64> {
    let mut g = 0i64;
    for &x in r {
        g = gcd_i64(g, x.abs());
    }
    let out: Vec<i64> = if g == 0 {
        r.to_vec()
    } else {
        r.iter().map(|&x| x / g).collect()
    };
    for &x in &out {
        if x > 0 {
            return out;
        } else if x < 0 {
            return out.iter().map(|&v| -v).collect();
        }
    }
    out
}

fn gcd_i64(a: i64, b: i64) -> i64 {
    if b == 0 { a } else { gcd_i64(b, a % b) }
}

/// Compute GV invariants using cygv.
///
/// This requires:
/// - Mori cone cap generators (rays)
/// - Grading vector (interior to dual cone)
/// - Curve basis matrix (q)
/// - Intersection numbers in basis (intnums)
pub fn compute_gv_invariants(
    rays: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<Vec<(Vec<i32>, Integer)>> {
    let t0 = std::time::Instant::now();
    eprintln!(
        "[DEBUG] gv start: rays={}, dim={}, max_deg={:?}, min_points={:?}",
        rays.len(),
        rays.first().map_or(0, Vec::len),
        max_deg,
        min_points
    );
    if min_points.is_none() && max_deg.is_none() {
        return Err(Error::InvalidInput(
            "Either min_points or max_deg must be specified".into(),
        ));
    }

    let dim = rays[0].len();

    // Report generator degree range for diagnostics.
    if let Some(d) = max_deg {
        let mut min_deg: Option<i128> = None;
        let mut max_deg_seen: Option<i128> = None;
        for r in rays {
            let deg: i128 = r
                .iter()
                .zip(grading_vector.iter())
                .map(|(&x, &g)| i128::from(x) * i128::from(g))
                .sum();
            min_deg = Some(min_deg.map_or(deg, |m| m.min(deg)));
            max_deg_seen = Some(max_deg_seen.map_or(deg, |m| m.max(deg)));
        }
        eprintln!(
            "[DEBUG] gv generators: total={}, degree_range={:?}-{:?}, max_deg={}",
            rays.len(),
            min_deg,
            max_deg_seen,
            d
        );
    }

    let n_rays = rays.len();
    eprintln!("[DEBUG] gv generators used: {n_rays}");
    if let Ok(path) = env::var("CYRUS_GV_DUMP_MORI_RAYS_CDD") {
        dump_mori_rays_cdd(Path::new(&path), rays)?;
    }

    let grading_vec_i32: Vec<i32> = grading_vector
        .iter()
        .map(|&v| i32::try_from(v))
        .collect::<std::result::Result<_, _>>()
        .map_err(|_| Error::InvalidInput("grading vector does not fit in i32".into()))?;
    let grading = RowDVector::from_row_slice(&grading_vec_i32);
    let mut min_g = i32::MAX;
    let mut max_g = i32::MIN;
    let mut neg_g = 0usize;
    for &v in &grading_vec_i32 {
        min_g = min_g.min(v);
        max_g = max_g.max(v);
        if v < 0 {
            neg_g += 1;
        }
    }
    eprintln!(
        "[DEBUG] gv grading_vec abs max={}, min={}, max={}, neg_count={}, elapsed={:.2?}",
        grading_vec_i32.iter().map(|v| v.abs()).max().unwrap_or(0),
        min_g,
        max_g,
        neg_g,
        t0.elapsed()
    );

    // Augment generators with lattice points, matching CYTools:
    // lattice_pts = mori.find_lattice_points(min_points=100*h11)
    let factor = env::var("CYRUS_LATTICE_MIN_POINTS_FACTOR")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(100);
    let gen_min_points = factor * dim;
    eprintln!("[DEBUG] gv generator min_points: {gen_min_points}");

    let lattice_pts = {
        let lattice_cache = LatticeCacheControls::from_env(1000, 0);
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        LATTICE_CACHE_VERSION.hash(&mut hasher);
        let mut key_rays: Vec<&[i64]> = rays.iter().map(Vec::as_slice).collect();
        key_rays.sort();
        for row in &key_rays {
            row.hash(&mut hasher);
        }
        grading_vector.hash(&mut hasher);
        gen_min_points.hash(&mut hasher);
        lattice_cache.max_coord.hash(&mut hasher);
        lattice_cache.deg_window.hash(&mut hasher);
        let key = hasher.finish();

        let cache_dir = env::var("CYRUS_CACHE_DIR")
            .map_or_else(|_| PathBuf::from("target/cyrus-cache"), PathBuf::from);
        let cache_path = cache_dir.join(format!("lattice_points_{key:x}.json"));

        eprintln!("[DEBUG] lattice cache path: {}", cache_path.display());
        if lattice_cache.cacheable && cache_path.exists() {
            let data = fs::read_to_string(&cache_path).map_err(|e| {
                Error::InvalidInput(format!(
                    "Failed to read lattice point cache {}: {e}",
                    cache_path.display()
                ))
            })?;
            let pts: Vec<Vec<i64>> = serde_json::from_str(&data).map_err(|e| {
                Error::InvalidInput(format!(
                    "Failed to parse lattice point cache {}: {e}",
                    cache_path.display()
                ))
            })?;
            eprintln!("[DEBUG] gv lattice points: {} (cache hit)", pts.len());
            pts
        } else {
            let rays_i128: Vec<Vec<i128>> = rays
                .iter()
                .map(|r| r.iter().map(|&x| i128::from(x)).collect())
                .collect();
            let mut cone = Cone::from_rays(rays_i128);
            let pts = cone.find_lattice_points_ortools(
                Some(gen_min_points),
                None,
                grading_vector,
                lattice_cache.max_coord,
                lattice_cache.deg_window,
            )?;
            eprintln!("[DEBUG] gv lattice points: {} (cache miss)", pts.len());
            if !lattice_cache.cacheable {
                eprintln!("[DEBUG] gv lattice points: cache disabled by search overrides");
            } else if let Err(e) = fs::create_dir_all(&cache_dir) {
                eprintln!(
                    "[WARN] failed to create lattice cache dir {}: {}",
                    cache_dir.display(),
                    e
                );
            } else {
                match fs::File::create(&cache_path) {
                    Ok(file) => {
                        let writer = BufWriter::new(file);
                        if let Err(e) = serde_json::to_writer(writer, &pts) {
                            eprintln!(
                                "[WARN] failed to serialize lattice cache {}: {}",
                                cache_path.display(),
                                e
                            );
                        }
                    }
                    Err(e) => {
                        eprintln!(
                            "[WARN] failed to create lattice cache {}: {}",
                            cache_path.display(),
                            e
                        );
                    }
                }
            }
            pts
        }
    };

    let mut all_generators: Vec<Vec<i64>> = Vec::new();
    for r in rays {
        all_generators.push(r.clone());
    }
    for p in lattice_pts {
        all_generators.push(p);
    }

    let mut uniq: HashSet<Vec<i64>> = HashSet::new();
    let mut uniq_vecs: Vec<Vec<i64>> = Vec::new();
    for v in all_generators {
        if uniq.insert(v.clone()) {
            uniq_vecs.push(v);
        }
    }

    let mut filtered_vecs = uniq_vecs;
    if let Some(d) = max_deg {
        let d_i128 = i128::from(d);
        let before = filtered_vecs.len();
        filtered_vecs.retain(|v| {
            let deg: i128 = v
                .iter()
                .zip(grading_vector.iter())
                .map(|(&x, &g)| i128::from(x) * i128::from(g))
                .sum();
            deg <= d_i128
        });
        eprintln!(
            "[DEBUG] gv generators filtered by max_deg: {} -> {}",
            before,
            filtered_vecs.len()
        );
        if filtered_vecs.is_empty() {
            return Err(Error::InvalidInput(
                "No generators remain after max_deg filtering".into(),
            ));
        }
    }

    let n_gen = filtered_vecs.len();
    let mut gen_data: Vec<i32> = Vec::with_capacity(dim * n_gen);
    for col in 0..n_gen {
        for row in 0..dim {
            let val = i32::try_from(filtered_vecs[col][row]).map_err(|_| {
                Error::InvalidInput("Mori cone generator does not fit in i32".into())
            })?;
            gen_data.push(val);
        }
    }
    let generators = DMatrix::from_column_slice(dim, n_gen, &gen_data);

    let generator_hash = {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        let mut key_vecs: Vec<&[i64]> = filtered_vecs.iter().map(Vec::as_slice).collect();
        key_vecs.sort();
        for row in &key_vecs {
            row.hash(&mut hasher);
        }
        hasher.finish()
    };

    // q matrix (curve basis), rows = h11, cols = n_pts
    let q_rows = q_matrix.len();
    if q_rows == 0 {
        return Err(Error::InvalidInput("q_matrix is empty".into()));
    }
    if q_rows != dim {
        return Err(Error::InvalidInput(
            "q_matrix row count must match generator dimension".into(),
        ));
    }
    let q_cols = q_matrix[0].len();
    let mut q_data: Vec<i32> = Vec::with_capacity(q_rows * q_cols);
    for row in q_matrix {
        if row.len() != q_cols {
            return Err(Error::InvalidInput(
                "q_matrix rows have inconsistent length".into(),
            ));
        }
        for &v in row {
            let v_i32 = i32::try_from(v)
                .map_err(|_| Error::InvalidInput("q_matrix entry does not fit in i32".into()))?;
            q_data.push(v_i32);
        }
    }
    let q = DMatrix::from_row_slice(q_rows, q_cols, &q_data);
    // cygv expects q with shape (n_divisors, h11), i.e. transpose of curve basis.
    let q = q.transpose();
    eprintln!(
        "[DEBUG] gv q shape: {}x{}, elapsed={:.2?}",
        q.nrows(),
        q.ncols(),
        t0.elapsed()
    );

    // Intersection numbers (dok format)
    let mut intnums_map: HashMap<(usize, usize, usize), i32> = HashMap::new();
    for (&(i, j, k), val) in intnums.iter() {
        if val.get().denominator_ref() != &1u32 {
            return Err(Error::InvalidInput(
                "intersection number is not integral".into(),
            ));
        }
        let signed = Integer::try_from(val.get().clone())
            .map_err(|_| Error::InvalidInput("intersection number is not integral".into()))?;
        let v_i64: i64 = i64::try_from(&signed)
            .map_err(|_| Error::InvalidInput("intersection number does not fit in i64".into()))?;
        let v_i32: i32 = v_i64
            .try_into()
            .map_err(|_| Error::InvalidInput("intersection number does not fit in i32".into()))?;
        intnums_map.insert((i, j, k), v_i32);
    }

    if let Ok(path) = env::var("CYRUS_GV_DUMP_INPUTS") {
        dump_gv_inputs(
            Path::new(&path),
            &filtered_vecs,
            &grading_vec_i32,
            &q,
            &intnums_map,
            min_points,
            max_deg,
        )?;
    }

    let cache_key = {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        grading_vec_i32.hash(&mut hasher);
        max_deg.hash(&mut hasher);
        min_points.hash(&mut hasher);
        generator_hash.hash(&mut hasher);
        q.iter().for_each(|v| v.hash(&mut hasher));
        let mut intnums_items: Vec<_> = intnums_map.iter().collect();
        intnums_items.sort_unstable_by_key(|(k, _)| **k);
        for (k, v) in intnums_items {
            k.hash(&mut hasher);
            v.hash(&mut hasher);
        }
        hasher.finish()
    };
    let cache_dir = env::var("CYRUS_CACHE_DIR")
        .map_or_else(|_| PathBuf::from("target/cyrus-cache"), PathBuf::from);
    let gv_cache_path = cache_dir.join(format!("gv_invariants_{cache_key:x}.json"));
    eprintln!("[DEBUG] gv cache path: {}", gv_cache_path.display());
    if gv_cache_path.exists() {
        let data = fs::read_to_string(&gv_cache_path).map_err(|e| {
            Error::InvalidInput(format!(
                "Failed to read GV cache {}: {e}",
                gv_cache_path.display()
            ))
        })?;
        #[derive(serde::Deserialize)]
        struct CachedGv {
            charge: Vec<i32>,
            value: String,
        }
        let items: Vec<CachedGv> = serde_json::from_str(&data).map_err(|e| {
            Error::InvalidInput(format!(
                "Failed to parse GV cache {}: {e}",
                gv_cache_path.display()
            ))
        })?;
        let mut out = Vec::with_capacity(items.len());
        for item in items {
            let gv_int = item
                .value
                .parse::<Integer>()
                .map_err(|()| Error::InvalidInput("GV cache integer parse failed".into()))?;
            out.push((item.charge, gv_int));
        }
        eprintln!("[DEBUG] gv invariants: cache hit ({} entries)", out.len());
        return Ok(out);
    }

    // No nef-partition for hypersurfaces
    let nefpart: Vec<DVector<i32>> = Vec::new();

    let n_threads = env::var("CYRUS_GV_THREADS")
        .ok()
        .and_then(|v| v.parse::<u32>().ok());
    let pool_size = env::var("CYRUS_GV_POOL_SIZE")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(1000);

    let invariants = cygv::compute_gv_rat_threefold(
        generators,
        grading,
        max_deg,
        min_points,
        q,
        nefpart,
        intnums_map,
        n_threads,
        pool_size,
    );

    let mut out = Vec::with_capacity(invariants.len());
    for (v, gv) in invariants {
        let gv_str = gv.to_string();
        let gv_int = gv_str
            .parse::<Integer>()
            .map_err(|()| Error::InvalidInput("GV integer conversion failed".into()))?;
        out.push((v.as_slice().to_vec(), gv_int));
    }

    if let Err(e) = fs::create_dir_all(&cache_dir) {
        eprintln!(
            "[WARN] failed to create GV cache dir {}: {}",
            cache_dir.display(),
            e
        );
    } else {
        #[derive(serde::Serialize)]
        struct CachedGv<'a> {
            charge: &'a [i32],
            value: String,
        }
        let payload: Vec<CachedGv<'_>> = out
            .iter()
            .map(|(c, v)| CachedGv {
                charge: c,
                value: v.to_string(),
            })
            .collect();
        match fs::File::create(&gv_cache_path) {
            Ok(file) => {
                let writer = BufWriter::new(file);
                if let Err(e) = serde_json::to_writer(writer, &payload) {
                    eprintln!(
                        "[WARN] failed to serialize GV cache {}: {}",
                        gv_cache_path.display(),
                        e
                    );
                }
            }
            Err(e) => {
                eprintln!(
                    "[WARN] failed to create GV cache {}: {}",
                    gv_cache_path.display(),
                    e
                );
            }
        }
    }

    Ok(out)
}

/// Map GV invariants from Kähler-basis curve coordinates to ambient divisor coordinates.
///
/// `curve_basis_matrix` is the matrix returned by
/// [`crate::compute_curve_basis_matrix`]. Its rows express each Kähler-basis
/// curve coordinate as an ambient divisor-intersection row, including the
/// origin column. The returned curve classes can therefore be compared against
/// `compute_mori_cone_cap_rays(..., in_basis=false, exclude_origin=false, ...)`
/// rows.
pub fn map_basis_gv_invariants_to_ambient(
    gv_invariants: &[(Vec<i32>, Integer)],
    curve_basis_matrix: &[Vec<Integer>],
) -> Result<Vec<(Vec<i64>, Integer)>> {
    if curve_basis_matrix.is_empty() {
        return Err(Error::InvalidInput("curve basis matrix is empty".into()));
    }
    let ambient_dim = curve_basis_matrix[0].len();
    if ambient_dim == 0 {
        return Err(Error::InvalidInput(
            "curve basis matrix has empty ambient rows".into(),
        ));
    }
    if curve_basis_matrix
        .iter()
        .any(|row| row.len() != ambient_dim)
    {
        return Err(Error::InvalidInput(
            "curve basis matrix rows have inconsistent length".into(),
        ));
    }

    let basis_dim = curve_basis_matrix.len();
    let mut out = Vec::with_capacity(gv_invariants.len());
    for (curve, gv) in gv_invariants {
        if curve.len() != basis_dim {
            return Err(Error::InvalidInput(format!(
                "GV curve dimension {} does not match curve basis dimension {}",
                curve.len(),
                basis_dim
            )));
        }

        let mut ambient = vec![Integer::from(0); ambient_dim];
        for (&coeff, basis_row) in curve.iter().zip(curve_basis_matrix.iter()) {
            if coeff == 0 {
                continue;
            }
            let coeff = Integer::from(coeff);
            for (entry, basis_coeff) in ambient.iter_mut().zip(basis_row.iter()) {
                *entry += &coeff * basis_coeff;
            }
        }

        let ambient = ambient
            .iter()
            .map(|entry| {
                i64::try_from(entry).map_err(|_| {
                    Error::InvalidInput("ambient GV curve coordinate does not fit in i64".into())
                })
            })
            .collect::<Result<Vec<_>>>()?;
        out.push((ambient, gv.clone()));
    }
    Ok(out)
}

fn dump_gv_inputs(
    path: &Path,
    generators: &[Vec<i64>],
    grading_vector: &[i32],
    q: &DMatrix<i32>,
    intnums: &HashMap<(usize, usize, usize), i32>,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<()> {
    #[derive(serde::Serialize)]
    struct IntNumDump {
        i: usize,
        j: usize,
        k: usize,
        value: i32,
    }

    #[derive(serde::Serialize)]
    struct GvInputDump<'a> {
        generators: &'a [Vec<i64>],
        grading_vector: &'a [i32],
        q: Vec<Vec<i32>>,
        intnums: Vec<IntNumDump>,
        min_points: Option<u32>,
        max_deg: Option<u32>,
    }

    let q_rows = (0..q.nrows())
        .map(|r| (0..q.ncols()).map(|c| q[(r, c)]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut intnum_rows = intnums
        .iter()
        .map(|(&(i, j, k), &value)| IntNumDump { i, j, k, value })
        .collect::<Vec<_>>();
    intnum_rows.sort_unstable_by_key(|row| (row.i, row.j, row.k));
    let payload = GvInputDump {
        generators,
        grading_vector,
        q: q_rows,
        intnums: intnum_rows,
        min_points,
        max_deg,
    };
    let json = serde_json::to_string_pretty(&payload)
        .map_err(|err| Error::InvalidInput(format!("GV input dump encode failed: {err}")))?;
    fs::write(path, json)
        .map_err(|err| Error::InvalidInput(format!("GV input dump write failed: {err}")))?;
    eprintln!("[DEBUG] wrote GV input dump: {}", path.display());
    Ok(())
}

fn dump_mori_rays_cdd(path: &Path, rays: &[Vec<i64>]) -> Result<()> {
    let dim = rays.first().map_or(0, Vec::len);
    for row in rays {
        if row.len() != dim {
            return Err(Error::InvalidInput(
                "cannot dump CDD rays with inconsistent dimensions".into(),
            ));
        }
    }

    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent).map_err(|err| {
            Error::InvalidInput(format!(
                "failed to create CDD dump directory {}: {err}",
                parent.display()
            ))
        })?;
    }

    let mut out = String::new();
    out.push_str("V-representation\nbegin\n");
    out.push_str(&format!("{} {} integer\n", rays.len(), dim + 1));
    for row in rays {
        out.push('0');
        for value in row {
            out.push(' ');
            out.push_str(&value.to_string());
        }
        out.push('\n');
    }
    out.push_str("end\n");
    fs::write(path, out).map_err(|err| {
        Error::InvalidInput(format!(
            "failed to write CDD dump {}: {err}",
            path.display()
        ))
    })?;
    eprintln!("[DEBUG] wrote Mori ray CDD dump: {}", path.display());
    Ok(())
}

#[derive(Clone, Debug)]
struct TwoFaceData {
    points: Vec<usize>,
    genus: usize,
}

fn compute_two_face_data_4d(points: &[Point], polytope: &Polytope) -> Result<Vec<TwoFaceData>> {
    let dual_vertices = polytope.dual_vertices()?;
    if dual_vertices.is_empty() {
        return Err(Error::InvalidInput("no dual vertices found".into()));
    }

    let hull_vertices = hull_vertex_coords(polytope)?;
    let mut facet_vertex_sets: Vec<HashSet<usize>> = Vec::with_capacity(dual_vertices.len());
    for dv in &dual_vertices {
        let mut vert_set: HashSet<usize> = HashSet::new();
        for (idx, vtx) in hull_vertices.iter().enumerate() {
            let dot: i64 = vtx
                .iter()
                .zip(dv.coords().iter())
                .map(|(&a, &b)| a * b)
                .sum();
            if dot == -1 {
                vert_set.insert(idx);
            }
        }
        facet_vertex_sets.push(vert_set);
    }

    let mut out = Vec::new();
    for i in 0..facet_vertex_sets.len() {
        for j in (i + 1)..facet_vertex_sets.len() {
            let inter_vertices = facet_vertex_sets[i]
                .intersection(&facet_vertex_sets[j])
                .count();
            if inter_vertices < 3 {
                continue;
            }

            let mut face_pts: Vec<usize> = Vec::new();
            for (idx, pt) in points.iter().enumerate() {
                let dot_i: i64 = pt
                    .coords()
                    .iter()
                    .zip(dual_vertices[i].coords().iter())
                    .map(|(&a, &b)| a * b)
                    .sum();
                if dot_i != -1 {
                    continue;
                }
                let dot_j: i64 = pt
                    .coords()
                    .iter()
                    .zip(dual_vertices[j].coords().iter())
                    .map(|(&a, &b)| a * b)
                    .sum();
                if dot_j == -1 {
                    face_pts.push(idx);
                }
            }
            face_pts.sort_unstable();
            out.push(TwoFaceData {
                points: face_pts,
                genus: lattice_segment_interior_points(
                    dual_vertices[i].coords(),
                    dual_vertices[j].coords(),
                ),
            });
        }
    }
    Ok(out)
}

fn classify_primal_point_face_dimensions(
    points: &[Point],
    polytope: &Polytope,
) -> Result<Vec<usize>> {
    let dual_vertices = polytope.dual_vertices()?;
    Ok(points
        .iter()
        .map(|pt| {
            let active: Vec<Vec<i64>> = dual_vertices
                .iter()
                .filter(|dv| {
                    pt.coords()
                        .iter()
                        .zip(dv.coords().iter())
                        .map(|(&a, &b)| a * b)
                        .sum::<i64>()
                        == -1
                })
                .map(|dv| dv.coords().to_vec())
                .collect();
            4usize.saturating_sub(integer_matrix_rank(&active))
        })
        .collect())
}

fn compute_primal_one_face_genera(
    points: &[Point],
    polytope: &Polytope,
    point_face_dims: &[usize],
) -> Result<Vec<Option<usize>>> {
    let dual_polytope = polytope.compute_dual()?;
    let dual_points = dual_polytope.vertices();
    let primal_hull_vertices = hull_vertex_coords(polytope)?;

    let dual_face_dims: Vec<usize> = dual_points
        .iter()
        .map(|dual_pt| {
            let active: Vec<Vec<i64>> = primal_hull_vertices
                .iter()
                .filter(|vtx| {
                    vtx.iter()
                        .zip(dual_pt.coords().iter())
                        .map(|(&a, &b)| a * b)
                        .sum::<i64>()
                        == -1
                })
                .cloned()
                .collect();
            4usize.saturating_sub(integer_matrix_rank(&active))
        })
        .collect();

    Ok(points
        .iter()
        .enumerate()
        .map(|(idx, pt)| {
            if point_face_dims.get(idx).copied() != Some(1) {
                return None;
            }
            Some(
                dual_points
                    .iter()
                    .zip(dual_face_dims.iter())
                    .filter(|(dual_pt, dual_face_dim)| {
                        **dual_face_dim == 2
                            && pt
                                .coords()
                                .iter()
                                .zip(dual_pt.coords().iter())
                                .map(|(&a, &b)| a * b)
                                .sum::<i64>()
                                == -1
                    })
                    .count(),
            )
        })
        .collect())
}

fn hull_vertex_coords(polytope: &Polytope) -> Result<Vec<Vec<i64>>> {
    let all_points: Vec<Vec<i64>> = polytope
        .vertices()
        .iter()
        .map(|p| p.coords().to_vec())
        .collect();
    let hull = ConvexHull::compute(&all_points)
        .ok_or_else(|| Error::InvalidInput("failed to compute convex hull".into()))?;
    Ok(hull
        .vertex_indices
        .iter()
        .map(|&idx| all_points[idx].clone())
        .collect())
}

fn lattice_segment_interior_points(a: &[i64], b: &[i64]) -> usize {
    let mut g = 0i64;
    for (&ai, &bi) in a.iter().zip(b.iter()) {
        g = gcd_i64(g, (bi - ai).abs());
    }
    g.saturating_sub(1) as usize
}

fn integer_matrix_rank(rows: &[Vec<i64>]) -> usize {
    let rational_rows: Vec<Vec<Rational>> = rows
        .iter()
        .map(|row| row.iter().map(|&value| Rational::from(value)).collect())
        .collect();
    matrix_rank(&rational_rows)
}

fn normalized_row_from_sparse_relation(n_points: usize, relation: Vec<(usize, i64)>) -> Vec<i64> {
    let mut row = vec![0i64; n_points];
    for (idx, coeff) in relation {
        row[idx] = coeff;
    }
    let mut g = 0i64;
    for &x in &row {
        g = gcd_i64(g, x.abs());
    }
    if g == 0 {
        return row;
    }
    row.into_iter().map(|x| x / g).collect()
}

fn compute_faces_4d(
    points: &[Point],
    polytope: &Polytope,
) -> Result<(Vec<Vec<usize>>, Vec<Vec<usize>>)> {
    if polytope.dim() != 4 {
        return Err(Error::InvalidInput(
            "faces_4d only defined for 4D polytopes".into(),
        ));
    }

    let dual_vertices = polytope.dual_vertices()?;
    if dual_vertices.is_empty() {
        return Err(Error::InvalidInput("no dual vertices found".into()));
    }

    // Build facet vertex sets (polytope vertices only) to define 2-faces.
    let poly_vertices = polytope.vertices();
    let mut facet_vertex_sets: Vec<HashSet<usize>> = Vec::with_capacity(dual_vertices.len());
    let mut facets: Vec<Vec<usize>> = Vec::with_capacity(dual_vertices.len());
    for dv in &dual_vertices {
        let mut vert_set: HashSet<usize> = HashSet::new();
        for (idx, vtx) in poly_vertices.iter().enumerate() {
            let dot: i64 = vtx
                .coords()
                .iter()
                .zip(dv.coords().iter())
                .map(|(&a, &b)| a * b)
                .sum();
            if dot == -1 {
                vert_set.insert(idx);
            }
        }
        facet_vertex_sets.push(vert_set);

        // Facet points (triangulation points only, matching CYTools boundary_points()).
        let mut facet_pts: Vec<usize> = Vec::new();
        for (idx, pt) in points.iter().enumerate() {
            let dot: i64 = pt
                .coords()
                .iter()
                .zip(dv.coords().iter())
                .map(|(&a, &b)| a * b)
                .sum();
            if dot == -1 {
                facet_pts.push(idx);
            }
        }
        facet_pts.sort_unstable();
        facets.push(facet_pts);
    }

    // Build 2-faces via intersections of facet vertex sets (CYTools _faces4d).
    let mut twofaces: Vec<Vec<usize>> = Vec::new();
    for i in 0..facet_vertex_sets.len() {
        for j in (i + 1)..facet_vertex_sets.len() {
            let inter_vertices = facet_vertex_sets[i]
                .intersection(&facet_vertex_sets[j])
                .count();
            if inter_vertices >= 3 {
                // Collect triangulation points lying on both facets.
                let mut face_pts: Vec<usize> = Vec::new();
                for (idx, pt) in points.iter().enumerate() {
                    let dot_i: i64 = pt
                        .coords()
                        .iter()
                        .zip(dual_vertices[i].coords().iter())
                        .map(|(&a, &b)| a * b)
                        .sum();
                    if dot_i != -1 {
                        continue;
                    }
                    let dot_j: i64 = pt
                        .coords()
                        .iter()
                        .zip(dual_vertices[j].coords().iter())
                        .map(|(&a, &b)| a * b)
                        .sum();
                    if dot_j == -1 {
                        face_pts.push(idx);
                    }
                }
                face_pts.sort_unstable();
                twofaces.push(face_pts);
            }
        }
    }

    Ok((facets, twofaces))
}

fn nullspace_vector(
    pts_ext: &[Vec<i64>],
    diff_pts: &[usize],
    comm_pts: &[usize],
    require_unique: bool,
) -> Option<Vec<Integer>> {
    let rows = diff_pts.len() + comm_pts.len();
    let cols = pts_ext[0].len();
    let mut m: Vec<Vec<Integer>> = vec![vec![Integer::from(0); cols]; rows];

    for (r, &idx) in diff_pts.iter().enumerate() {
        for c in 0..cols {
            m[r][c] = Integer::from(pts_ext[idx][c]);
        }
    }
    for (r, &idx) in comm_pts.iter().enumerate() {
        for c in 0..cols {
            m[r + diff_pts.len()][c] = Integer::from(pts_ext[idx][c]);
        }
    }

    // Compute kernel of m^T
    let m_t = transpose(&m);
    let kernel = integer_kernel(&m_t);
    if kernel.is_empty() {
        return None;
    }
    if require_unique && kernel.len() != 1 {
        return None;
    }

    let mut v = kernel[0].clone();
    if v[0] < 0 {
        for val in &mut v {
            *val = -val.clone();
        }
    }

    let g = gcd_list(&v);
    if g != 0 {
        for val in &mut v {
            *val /= &g;
        }
    }

    Some(v)
}

fn build_full_v(diff_pts: &[usize], comm_pts: &[usize], v: &[Integer]) -> Vec<(usize, i64)> {
    let mut full_v: Vec<(usize, i64)> = Vec::new();
    for k in 0..diff_pts.len() {
        full_v.push((
            diff_pts[k],
            i64::try_from(&v[k]).expect("mori coeff fits in i64"),
        ));
    }
    for k in 0..comm_pts.len() {
        if v[k + diff_pts.len()] != 0 {
            full_v.push((
                comm_pts[k],
                i64::try_from(&v[k + diff_pts.len()]).expect("mori coeff fits in i64"),
            ));
        }
    }
    full_v.sort_by_key(|(idx, _)| *idx);
    full_v
}

fn transpose(m: &[Vec<Integer>]) -> Vec<Vec<Integer>> {
    if m.is_empty() {
        return Vec::new();
    }
    let rows = m.len();
    let cols = m[0].len();
    let mut t = vec![vec![Integer::from(0); rows]; cols];
    for r in 0..rows {
        for c in 0..cols {
            t[c][r] = m[r][c].clone();
        }
    }
    t
}

fn gcd_list(vals: &[Integer]) -> Integer {
    let mut g = Integer::from(0);
    for v in vals {
        let abs = v.clone().abs();
        if g == 0 {
            g = abs;
        } else if abs != 0 {
            g = gcd_integer(&g, &abs);
        }
    }
    g
}

fn is_strictly_dual(rays: &[Vec<i64>], v: &[i64]) -> bool {
    for r in rays {
        if r.len() != v.len() {
            return false;
        }
        let dot: i128 = r
            .iter()
            .zip(v.iter())
            .map(|(&a, &b)| i128::from(a) * i128::from(b))
            .sum();
        if dot <= 0 {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::{
        BoundedCurveDecompositionIndex, ToricCurveCandidate, compute_grading_vector,
        curve_volume_in_divisor_basis, dump_mori_rays_cdd, find_pair_decomposition,
        load_grading_cache, map_basis_gv_invariants_to_ambient,
        remove_pair_decomposable_curve_candidates, subcutoff_toric_curve_candidates,
        write_grading_cache,
    };
    use crate::{f64_finite, f64_pos};
    use malachite::Integer;

    #[test]
    fn basis_gv_invariants_map_to_ambient_curve_classes() {
        let curve_basis = vec![
            vec![Integer::from(0), Integer::from(1), Integer::from(-2)],
            vec![Integer::from(0), Integer::from(3), Integer::from(4)],
        ];
        let invariants = vec![(vec![2, -1], Integer::from(7))];

        let mapped = map_basis_gv_invariants_to_ambient(&invariants, &curve_basis).unwrap();

        assert_eq!(mapped, vec![(vec![0, -1, -8], Integer::from(7))]);
    }

    #[test]
    fn basis_gv_invariant_mapping_rejects_dimension_mismatch() {
        let curve_basis = vec![vec![Integer::from(1), Integer::from(0)]];
        let invariants = vec![(vec![1, 2], Integer::from(1))];

        assert!(map_basis_gv_invariants_to_ambient(&invariants, &curve_basis).is_err());
    }

    #[test]
    fn mori_ray_cdd_dump_uses_cdd_v_representation_for_cone_rays() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time is after unix epoch")
            .as_nanos();
        let cache_dir = PathBuf::from("target/cyrus-test-cache")
            .join(format!("cdd-{}-{nonce}", std::process::id()));
        let path = cache_dir.join("rays.ext");

        dump_mori_rays_cdd(&path, &[vec![1, 0], vec![0, -2]]).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert_eq!(
            content,
            "V-representation\nbegin\n2 3 integer\n0 1 0\n0 0 -2\nend\n"
        );

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_dir_all(&cache_dir);
    }

    #[test]
    fn grading_cache_validates_candidate_against_rays() {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system time is after unix epoch")
            .as_nanos();
        let cache_dir = PathBuf::from("target/cyrus-test-cache")
            .join(format!("grading-{}-{nonce}", std::process::id()));
        let path = cache_dir.join("grading.json");
        let rays = vec![vec![1, 0], vec![0, 1]];

        write_grading_cache(&cache_dir, &path, &[2, 4]);
        assert_eq!(load_grading_cache(&path, &rays), Some(vec![2, 4]));

        write_grading_cache(&cache_dir, &path, &[1, 0]);
        assert_eq!(load_grading_cache(&path, &rays), None);

        write_grading_cache(&cache_dir, &path, &[1, 1, 1]);
        assert_eq!(load_grading_cache(&path, &rays), None);

        let _ = std::fs::remove_file(&path);
        let _ = std::fs::remove_dir_all(&cache_dir);
    }

    #[test]
    fn grading_vector_matches_cytools_mcallister_dual() {
        let rays = vec![
            vec![0, 0, 1, -2],
            vec![0, 0, 0, 1],
            vec![-1, 0, 0, 1],
            vec![1, 1, 0, 0],
            vec![1, 0, 0, 0],
            vec![0, 1, 0, 0],
            vec![0, 0, 1, 0],
            vec![0, 1, 1, 0],
            vec![1, 1, 0, -2],
            vec![-1, 0, 1, 0],
        ];

        assert_eq!(compute_grading_vector(&rays), Some(vec![1, 4, 5, 2]));
    }

    #[test]
    fn curve_volume_in_divisor_basis_allows_mixed_basis_coordinates() {
        let curve = vec![1, -2, 3];
        let basis = vec![0, 1, 2];
        let kahler = vec![f64_finite!(2.0), f64_finite!(-1.0), f64_finite!(0.5)];

        let volume = curve_volume_in_divisor_basis(&curve, &basis, &kahler).unwrap();

        assert_eq!(volume, f64_finite!(5.5));
    }

    #[test]
    fn subcutoff_toric_curve_candidates_keeps_only_positive_below_cutoff() {
        let rays = vec![vec![1, 0], vec![-1, 0], vec![3, 0], vec![0, 1]];
        let basis = vec![0];
        let kahler = vec![f64_finite!(0.4)];

        let selected =
            subcutoff_toric_curve_candidates(&rays, &basis, &kahler, f64_pos!(1.0)).unwrap();

        assert_eq!(
            selected,
            vec![ToricCurveCandidate {
                class: vec![1, 0],
                volume: f64_pos!(0.4)
            }]
        );
    }

    #[test]
    fn pair_decomposable_curve_candidates_are_removed() {
        let candidates = vec![
            ToricCurveCandidate {
                class: vec![1, 0],
                volume: f64_pos!(0.2),
            },
            ToricCurveCandidate {
                class: vec![0, 1],
                volume: f64_pos!(0.3),
            },
            ToricCurveCandidate {
                class: vec![1, 1],
                volume: f64_pos!(0.5),
            },
            ToricCurveCandidate {
                class: vec![2, 0],
                volume: f64_pos!(0.4),
            },
        ];

        let decomposition = find_pair_decomposition(&candidates[2], &candidates).unwrap();
        assert_eq!(decomposition, Some((vec![1, 0], vec![0, 1])));

        let filtered = remove_pair_decomposable_curve_candidates(&candidates).unwrap();
        assert_eq!(
            filtered,
            vec![
                ToricCurveCandidate {
                    class: vec![1, 0],
                    volume: f64_pos!(0.2),
                },
                ToricCurveCandidate {
                    class: vec![0, 1],
                    volume: f64_pos!(0.3),
                },
            ]
        );
    }

    #[test]
    fn pair_pruning_does_not_claim_full_semigroup_reduction() {
        let candidates = vec![
            ToricCurveCandidate {
                class: vec![1, 0, 0],
                volume: f64_pos!(0.2),
            },
            ToricCurveCandidate {
                class: vec![0, 1, 0],
                volume: f64_pos!(0.3),
            },
            ToricCurveCandidate {
                class: vec![0, 0, 1],
                volume: f64_pos!(0.4),
            },
            ToricCurveCandidate {
                class: vec![1, 1, 1],
                volume: f64_pos!(0.9),
            },
        ];

        let decomposition = find_pair_decomposition(&candidates[3], &candidates).unwrap();
        assert_eq!(decomposition, None);

        let filtered = remove_pair_decomposable_curve_candidates(&candidates).unwrap();
        assert_eq!(filtered, candidates);
    }

    #[test]
    fn bounded_decomposition_index_finds_three_term_sums() {
        let candidates = vec![
            ToricCurveCandidate {
                class: vec![1, 0, 0],
                volume: f64_pos!(0.2),
            },
            ToricCurveCandidate {
                class: vec![0, 1, 0],
                volume: f64_pos!(0.3),
            },
            ToricCurveCandidate {
                class: vec![0, 0, 1],
                volume: f64_pos!(0.4),
            },
            ToricCurveCandidate {
                class: vec![1, 1, 1],
                volume: f64_pos!(0.9),
            },
        ];
        let index = BoundedCurveDecompositionIndex::new(&candidates).unwrap();

        assert_eq!(index.find_decomposition(&candidates[3], 2).unwrap(), None);
        assert_eq!(
            index.find_decomposition(&candidates[3], 3).unwrap(),
            Some(vec![vec![0, 0, 1], vec![0, 1, 0], vec![1, 0, 0]])
        );
    }

    #[test]
    fn curve_volume_rejects_basis_dimension_mismatch() {
        let err = curve_volume_in_divisor_basis(&[1, 2], &[0, 1], &[f64_finite!(1.0)])
            .expect_err("basis and Kähler coordinate lengths must match");

        assert!(format!("{err}").contains("basis length"));
    }
}
