//! Gopakumar-Vafa invariant computation utilities.
//!
//! Implements the CYTools mori_cone_cap algorithm and grading vector selection,
//! then wires up the `cygv` crate for actual GV computation.
//! Compact hypersurface GV computation must stay on the upstream `cygv` crate;
//! the CKYZ helpers in this module are local/noncompact diagnostics and are not
//! a replacement for `cygv`'s HKTY implementation.
//!
//! Reference: CYTools `calabiyau.py` (mori_cone_cap) and `cone.py` (grading vector).

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet, VecDeque};
use std::env;
use std::f64::consts::{LN_10, PI};
use std::fs;
use std::hash::{Hash, Hasher};
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use good_lp::{
    Expression, ProblemVariables, ResolutionError, Solution, SolverModel, Variable, default_solver,
    variable,
};
use malachite::Integer;
use malachite::Rational;
use malachite::num::arithmetic::traits::Abs;
use nalgebra::{DMatrix, DVector, RowDVector};
use rug::Rational as RugRational;

use crate::cone::Cone;
use crate::curve_basis::{
    DivisorBasis, compute_curve_basis_matrix_for_divisor_basis,
    curve_basis_matrix_without_origin_i64,
};
use crate::error::{Error, Result};
use crate::geometry::ConvexHull;
use crate::integer_math::hermite_normal_form;
use crate::integer_math::matrix_rank;
use crate::integer_math::{gcd_integer, integer_kernel, solve_linear_system_rational};
use crate::intersection::Intersection;
use crate::lattice::Point;
use crate::polytope::Polytope;
use crate::triangulation::Triangulation;
use crate::types::{F64, Finite, I64, Pos};
use crate::utils::lll_reduce;

const GRADING_CACHE_VERSION: &str = "grading-vector-cytools-lp-v1";
const LATTICE_CACHE_VERSION: &str = "lattice-points-v2";
const GV_CACHE_VERSION: &str = "gv-invariants-cygv-0.1.2-v1";
const CKYZ_ADDITION_TABLE_MAX_ENTRIES: usize = 5_000_000;
const CKYZ_ABSENT_ADDITION_INDEX: usize = usize::MAX;
const CKYZ_DENSE_DEGREE_INDEX_MAX_ENTRIES: usize = 5_000_000;
const CKYZ_ABSENT_DEGREE_INDEX: usize = usize::MAX;

/// One term in a compact `q_N` polynomial materialized by cygv.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CygvQnTraceTerm {
    /// Index of the monomial in cygv's finite semigroup.
    pub monomial_index: usize,
    /// Semigroup exponent vector for this monomial.
    pub exponent: Vec<i32>,
    /// Exact coefficient, serialized with rug's rational display format.
    pub coefficient: String,
}

/// Compact `q_N` polynomial history exported from cygv series inversion.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CygvQnTracePolynomial {
    /// Index of the curve element whose `q_N` polynomial was computed.
    pub element_index: usize,
    /// Grading degree of the curve element.
    pub degree: u32,
    /// Semigroup exponent vector of the curve element.
    pub element: Vec<i32>,
    /// Nonzero terms of the `q_N` polynomial.
    pub terms: Vec<CygvQnTraceTerm>,
    /// Nonzero terms of `Li2(q_N)` after cygv's finite monomial-map truncation.
    pub li2_terms: Vec<CygvQnTraceTerm>,
}

/// Inverse-series GV candidate read from cygv's mutable instanton polynomial.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CygvGvCoefficientTrace {
    /// Index of the curve element whose coefficient was read.
    pub element_index: usize,
    /// Grading degree of the curve element.
    pub degree: u32,
    /// Semigroup exponent vector of the curve element.
    pub element: Vec<i32>,
    /// Instanton polynomial coordinate used by cygv to read the coefficient.
    pub insertion_index: usize,
    /// First nonzero component of `element`, used as the divisor.
    pub pivot_component: i32,
    /// Exact mutable instanton coefficient at this element, if present.
    pub instanton_coefficient: Option<String>,
    /// Exact GV candidate before rounding/filtering, if present.
    pub gv_candidate: Option<String>,
    /// Rounded GV candidate used by cygv's integrality check, if present.
    pub rounded_gv_candidate: Option<String>,
    /// cygv decision status for this candidate.
    pub status: String,
}

/// GV output together with the compact `q_N` polynomials cygv materialized.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GvInvariantsWithQnTrace {
    /// Integral GV invariants returned by cygv.
    pub invariants: Vec<(Vec<i32>, Integer)>,
    /// Compact `q_N` polynomials materialized while computing those invariants.
    pub qn_trace: Vec<CygvQnTracePolynomial>,
    /// GV coefficient candidates read before `q_N` materialization.
    pub gv_coefficient_trace: Vec<CygvGvCoefficientTrace>,
}

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

/// Project ambient Mori-cap rays to a divisor basis, then normalize and deduplicate.
///
/// This matches the final `in_basis=true` projection/normalization step of
/// [`compute_mori_cone_cap_rays`] when the caller already has ambient rays.
pub fn project_mori_cone_cap_rays_to_basis(
    ambient_rays: &[Vec<i64>],
    basis: &[usize],
) -> Result<Vec<Vec<i64>>> {
    let mut projected: Vec<Vec<i64>> = Vec::with_capacity(ambient_rays.len());
    for ray in ambient_rays {
        let mut row = Vec::with_capacity(basis.len());
        for &idx in basis {
            let Some(&value) = ray.get(idx) else {
                return Err(Error::InvalidInput(format!(
                    "basis index {idx} is out of bounds for Mori ray dimension {}",
                    ray.len()
                )));
            };
            row.push(value);
        }

        let mut g = 0i64;
        for &x in &row {
            g = gcd_i64(g, x.abs());
        }
        if g == 0 {
            continue;
        }
        projected.push(row.into_iter().map(|x| x / g).collect());
    }

    let mut uniq: HashSet<Vec<i64>> = HashSet::new();
    let mut deduped = Vec::with_capacity(projected.len());
    for row in projected {
        if uniq.insert(row.clone()) {
            deduped.push(row);
        }
    }
    deduped.sort();
    Ok(deduped)
}

/// Project ambient Mori-cap rays to a matrix divisor basis.
///
/// This matches CYTools' generic-basis path in `mori_cone_cap(in_basis=True)`:
/// if `basis` is a matrix, CYTools computes `mori_cap_matrix.dot(basis.T)`.
/// Rows of `basis_matrix` are divisor-basis vectors in ambient divisor
/// coordinates, including the origin column when the ambient rays include it.
pub fn project_mori_cone_cap_rays_to_basis_matrix(
    ambient_rays: &[Vec<i64>],
    basis_matrix: &[Vec<i64>],
) -> Result<Vec<Vec<i64>>> {
    let ambient_dim = ambient_rays.first().map_or_else(
        || {
            basis_matrix.first().map_or_else(
                || Err(Error::InvalidInput("basis matrix is empty".into())),
                |row| Ok(row.len()),
            )
        },
        |ray| Ok(ray.len()),
    )?;
    validate_basis_matrix(basis_matrix, ambient_dim)?;

    let mut projected: Vec<Vec<i64>> = Vec::with_capacity(ambient_rays.len());
    for ray in ambient_rays {
        if ray.len() != ambient_dim {
            return Err(Error::InvalidInput(
                "ambient Mori rays have inconsistent dimensions".into(),
            ));
        }
        let row = project_ambient_curve_to_basis_matrix(ray, basis_matrix)?;

        let mut g = 0i64;
        for &x in &row {
            g = gcd_i64(g, x.abs());
        }
        if g == 0 {
            continue;
        }
        projected.push(row.into_iter().map(|x| x / g).collect());
    }

    let mut uniq: HashSet<Vec<i64>> = HashSet::new();
    let mut deduped = Vec::with_capacity(projected.len());
    for row in projected {
        if uniq.insert(row.clone()) {
            deduped.push(row);
        }
    }
    deduped.sort();
    Ok(deduped)
}

/// CYTools-style divisor-basis data needed by direct `cygv` calls.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GvDivisorBasisData {
    /// Mori-cap rays projected to the selected divisor basis.
    pub mori_rays: Vec<Vec<i64>>,
    /// Dual curve-basis matrix in ambient coordinates, including the origin.
    pub curve_basis_matrix: Vec<Vec<Integer>>,
    /// No-origin `q` matrix passed to `cygv`.
    pub q_matrix: Vec<Vec<i64>>,
}

/// Compute in-basis intersection numbers for either CYTools divisor-basis shape.
///
/// Vector bases match CYTools' `filter_tensor_indices` path. Matrix bases match
/// CYTools' `symmetric_sparse_to_dense(tensor, basis)` contraction, where each
/// matrix row is a divisor-basis vector in ambient divisor coordinates.
///
/// # Errors
/// Returns an error if the selected basis is malformed for the ambient
/// intersection tensor.
pub fn intersection_in_divisor_basis(
    kappa: &Intersection,
    basis: DivisorBasis<'_>,
) -> Result<Intersection> {
    match basis {
        DivisorBasis::Indices(indices) => {
            validate_intersection_index_basis(kappa, indices)?;
            Ok(crate::basis::intersection_in_basis(kappa, indices))
        }
        DivisorBasis::Matrix {
            standard_basis,
            basis_matrix,
        } => {
            if standard_basis.len() != basis_matrix.len() {
                return Err(Error::InvalidInput(format!(
                    "standard basis length {} does not match matrix divisor basis row count {}",
                    standard_basis.len(),
                    basis_matrix.len()
                )));
            }
            intersection_in_matrix_divisor_basis(kappa, basis_matrix)
        }
    }
}

/// Transform ambient intersection numbers to a matrix divisor basis.
///
/// This is the rectangular version of the exact tensor pullback:
///
/// ```text
/// kappa'_{abc} = kappa_{ijk} B^a_i B^b_j B^c_k
/// ```
///
/// where `B` is the divisor-basis row matrix. The output dimension is the
/// number of rows of `basis_matrix`.
///
/// # Errors
/// Returns an error if the basis matrix is empty, ragged, or has a width that
/// does not match `kappa.dim()`.
pub fn intersection_in_matrix_divisor_basis(
    kappa: &Intersection,
    basis_matrix: &[Vec<Integer>],
) -> Result<Intersection> {
    validate_intersection_basis_matrix(kappa, basis_matrix)?;
    let column_supports = intersection_basis_column_supports(basis_matrix, kappa.dim());

    let mut entries: HashMap<(usize, usize, usize), Rational> = HashMap::new();
    for (&(i, j, k), value) in kappa.iter() {
        for [ambient_a, ambient_b, ambient_c] in unique_triple_permutations(i, j, k) {
            for (basis_a, coeff_ai) in &column_supports[ambient_a] {
                for (basis_b, coeff_bj) in &column_supports[ambient_b] {
                    if basis_a > basis_b {
                        continue;
                    }
                    for (basis_c, coeff_ck) in &column_supports[ambient_c] {
                        if basis_b > basis_c {
                            continue;
                        }
                        let mut term = value.get().clone();
                        term *= Rational::from(coeff_ai);
                        term *= Rational::from(coeff_bj);
                        term *= Rational::from(coeff_ck);
                        if term != 0 {
                            *entries
                                .entry((*basis_a, *basis_b, *basis_c))
                                .or_insert_with(|| Rational::from(0)) += term;
                        }
                    }
                }
            }
        }
    }

    let mut transformed = Intersection::new(basis_matrix.len());
    for ((a, b, c), value) in entries {
        if value != 0 {
            transformed.set(
                a,
                b,
                c,
                crate::types::rational::Rational::<Finite>::from_raw(value),
            );
        }
    }
    Ok(transformed)
}

fn intersection_basis_column_supports(
    basis_matrix: &[Vec<Integer>],
    ambient_dim: usize,
) -> Vec<Vec<(usize, Integer)>> {
    let mut supports = vec![Vec::new(); ambient_dim];
    for (basis_row, row) in basis_matrix.iter().enumerate() {
        for (ambient_col, value) in row.iter().enumerate() {
            if *value != 0 {
                supports[ambient_col].push((basis_row, value.clone()));
            }
        }
    }
    supports
}

fn validate_intersection_index_basis(kappa: &Intersection, basis: &[usize]) -> Result<()> {
    let mut seen = BTreeSet::new();
    for &idx in basis {
        if idx >= kappa.dim() {
            return Err(Error::InvalidInput(format!(
                "basis index {idx} is out of bounds for intersection dimension {}",
                kappa.dim()
            )));
        }
        if !seen.insert(idx) {
            return Err(Error::InvalidInput(
                "basis indices contain duplicates".into(),
            ));
        }
    }
    Ok(())
}

fn validate_intersection_basis_matrix(
    kappa: &Intersection,
    basis_matrix: &[Vec<Integer>],
) -> Result<()> {
    if basis_matrix.is_empty() {
        return Err(Error::InvalidInput("matrix divisor basis is empty".into()));
    }
    for row in basis_matrix {
        if row.len() != kappa.dim() {
            return Err(Error::InvalidInput(format!(
                "matrix divisor basis row width {} does not match intersection dimension {}",
                row.len(),
                kappa.dim()
            )));
        }
    }
    Ok(())
}

fn unique_triple_permutations(i: usize, j: usize, k: usize) -> Vec<[usize; 3]> {
    let mut permutations = vec![
        [i, j, k],
        [i, k, j],
        [j, i, k],
        [j, k, i],
        [k, i, j],
        [k, j, i],
    ];
    permutations.sort_unstable();
    permutations.dedup();
    permutations
}

/// Project ambient Mori-cap rays through either CYTools divisor-basis shape.
///
/// Vector bases use CYTools' column-selection path. Matrix bases use
/// `mori_cap_matrix.dot(basis.T)`.
///
/// # Errors
/// Returns an error if the selected basis is malformed, out of range, or cannot
/// be represented as `i64` for Mori projection.
pub fn project_mori_cone_cap_rays_for_divisor_basis(
    ambient_rays: &[Vec<i64>],
    basis: DivisorBasis<'_>,
) -> Result<Vec<Vec<i64>>> {
    match basis {
        DivisorBasis::Indices(indices) => {
            project_mori_cone_cap_rays_to_basis(ambient_rays, indices)
        }
        DivisorBasis::Matrix { basis_matrix, .. } => {
            let basis_matrix = integer_matrix_to_i64(basis_matrix, "divisor basis matrix")?;
            project_mori_cone_cap_rays_to_basis_matrix(ambient_rays, &basis_matrix)
        }
    }
}

/// Build projected Mori rays, ambient curve basis, and `cygv` q matrix together.
///
/// This is the source-shaped handoff for GA callers that may use either a
/// CYTools vector basis or a generic matrix divisor basis.
///
/// # Errors
/// Returns an error if Mori projection or dual curve-basis construction fails.
pub fn gv_divisor_basis_data(
    ambient_mori_rays: &[Vec<i64>],
    linrels: &[Vec<Integer>],
    basis: DivisorBasis<'_>,
) -> Result<GvDivisorBasisData> {
    let mori_rays = project_mori_cone_cap_rays_for_divisor_basis(ambient_mori_rays, basis)?;
    let curve_basis_matrix = compute_curve_basis_matrix_for_divisor_basis(linrels, basis)?;
    let q_matrix = curve_basis_matrix_without_origin_i64(&curve_basis_matrix)?;
    Ok(GvDivisorBasisData {
        mori_rays,
        curve_basis_matrix,
        q_matrix,
    })
}

fn integer_matrix_to_i64(matrix: &[Vec<Integer>], context: &str) -> Result<Vec<Vec<i64>>> {
    matrix
        .iter()
        .map(|row| integer_vector_to_i64(row, &format!("{context} entry")))
        .collect()
}

fn integer_vector_to_i64(vector: &[Integer], context: &str) -> Result<Vec<i64>> {
    vector
        .iter()
        .map(|value| {
            i64::try_from(value)
                .map_err(|_| Error::InvalidInput(format!("{context} does not fit in i64")))
        })
        .collect()
}

fn checked_neg_i64_vector(vector: &[i64], context: &str) -> Result<Vec<i64>> {
    vector
        .iter()
        .map(|&value| {
            value
                .checked_neg()
                .ok_or_else(|| Error::InvalidInput(format!("{context} negation overflowed")))
        })
        .collect()
}

/// A toric curve candidate with its volume at a specific point in Kähler moduli space.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ToricCurveCandidate {
    /// Curve class in ambient divisor-intersection coordinates.
    pub class: Vec<i64>,
    /// Positive curve volume at the selected Kähler point.
    pub volume: F64<Pos>,
}

/// Convergence data for the instanton series along one potent ray.
#[derive(Clone, Debug, PartialEq)]
pub struct PotentRayConvergence {
    /// `log(|N_{nq}| exp(-2π n q.t))` for `n = 1..`.
    ///
    /// A zero GV invariant has no finite log-magnitude and is represented by
    /// `None`.
    pub log_xi_terms: Vec<Option<F64<Finite>>>,
    /// Least-squares slope of the finite log-xi terms against `n`.
    pub log_xi_slope: Option<F64<Finite>>,
}

/// Genus-zero GV values along a one-dimensional ray.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OneDimensionalRayGvSeries {
    /// Primitive ray in Kähler-basis curve coordinates.
    pub ray: Vec<i64>,
    /// Positive grading degree of `ray`.
    pub degree: i128,
    /// GV values for `ray`, `2*ray`, ..., `max_multiple*ray`.
    pub values: Vec<Integer>,
}

/// Finite-cutoff evidence that a primitive ray is apparently nilpotent.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NilpotentRayCandidate {
    /// Primitive co-prime curve ray being tested.
    pub primitive_ray: Vec<i64>,
    /// First positive multiple with zero genus-zero GV invariant.
    pub first_vanishing_multiple: i128,
    /// Grading degree of `first_vanishing_multiple * primitive_ray`.
    pub first_vanishing_degree: i128,
    /// `sum_{k=1}^{k*-1} k^2 n^0_{kC}` at the first vanishing multiple.
    pub weighted_lower_gv_sum: Integer,
}

/// Finite-cutoff partition of nonzero GV charges by apparent nilpotence.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FiniteCutoffGvChargePartition {
    /// Primitive rays that pass the finite-cutoff nilpotence test.
    pub nilpotent_rays: Vec<NilpotentRayCandidate>,
    /// Nonzero charges up to the cutoff whose primitive ray is not in
    /// `nilpotent_rays`.
    pub potent_charges: Vec<(Vec<i64>, Integer)>,
}

/// Degree-slice origin for the finite-cutoff nop divergence test.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NilpotentRayDegreeSlice {
    /// Primitive candidate nilpotent ray.
    pub primitive_ray: Vec<i64>,
    /// Numerator of the cutoff fraction used for this slice.
    pub cutoff_numerator: i128,
    /// Denominator of the cutoff fraction used for this slice.
    pub cutoff_denominator: i128,
    /// Largest positive integer `k` satisfying
    /// `k * primitive_ray.degree <= cutoff_degree * numerator / denominator`.
    pub slice_multiple: i128,
    /// Grading degree of `slice_origin`.
    pub slice_degree: i128,
    /// Curve class `slice_multiple * primitive_ray`.
    pub slice_origin: Vec<i64>,
}

/// Integer point where a comparison ray intersects a nilpotent degree slice.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NilpotentRaySliceComparisonPoint {
    /// Primitive comparison ray.
    pub primitive_ray: Vec<i64>,
    /// Positive grading degree of `primitive_ray`.
    pub primitive_degree: i128,
    /// Positive multiple that places this ray on the slice.
    pub slice_multiple: i128,
    /// Curve class on the slice.
    pub slice_point: Vec<i64>,
    /// `slice_point - slice_origin`.
    pub offset_from_origin: Vec<i64>,
}

/// LLL-reduced distance from a nilpotent ray to comparison rays on one degree
/// slice.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NilpotentRaySliceDistance {
    /// Degree slice used for the comparison.
    pub slice: NilpotentRayDegreeSlice,
    /// Deduplicated nonzero offsets `C' - slice_origin` used to define the
    /// affine slice lattice before LLL reduction.
    pub lattice_offsets: Vec<Vec<i64>>,
    /// Unimodular transform `A` returned by the CYTools-style LLL convention:
    /// `reduced_offsets.T = A * lattice_offsets.T`.
    pub lll_transform: Vec<Vec<i64>>,
    /// LLL-reduced lattice offsets, in the same row order as
    /// `lattice_offsets`.
    pub reduced_lattice_offsets: Vec<Vec<i64>>,
    /// Comparison-ray intersections with the slice.
    pub comparison_points: Vec<NilpotentRaySliceComparisonPoint>,
    /// Each comparison offset transformed by `lll_transform`.
    pub transformed_comparison_offsets: Vec<Vec<i64>>,
    /// Smallest infinity norm among transformed comparison offsets, or `None`
    /// if no comparison ray intersects the slice at an integer point.
    pub minimum_infinity_norm: Option<i64>,
}

/// Half/full cutoff comparison for the finite-cutoff nop divergence test.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NilpotentRayDivergenceCheck {
    /// Distance on the half-cutoff slice.
    pub half_cutoff: NilpotentRaySliceDistance,
    /// Distance on the full-cutoff slice.
    pub full_cutoff: NilpotentRaySliceDistance,
    /// `Some(true)` exactly when the full-cutoff distance is strictly larger
    /// than the half-cutoff distance. `None` means one of the slices had no
    /// comparison point, so the finite test is inconclusive at this layer.
    pub appears_divergent: Option<bool>,
}

/// Result of the paper's two-pass finite-cutoff nop classification once the
/// required divergence checks have been computed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NilpotentRayTwoPassNopClassification {
    /// Provisional nop rays `F0` from the first pass against `C \ N`.
    pub initial_candidate_nop_rays: Vec<NilpotentRayCandidate>,
    /// Final nop rays `F` that still diverge in the second pass against
    /// `C \ F0`.
    pub nop_rays: Vec<NilpotentRayCandidate>,
    /// Nilpotent rays whose first-pass distance comparison had no definite
    /// `d' > d` answer.
    pub first_pass_inconclusive_rays: Vec<NilpotentRayCandidate>,
    /// Provisional nop rays whose second-pass distance comparison had no
    /// definite `d' > d` answer.
    pub second_pass_inconclusive_rays: Vec<NilpotentRayCandidate>,
}

/// Full finite-GV-table result for the two-pass nop classification.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FiniteGvTableNopClassification {
    /// Initial partition into apparently nilpotent rays and `C \ N`.
    pub partition: FiniteCutoffGvChargePartition,
    /// First-pass checks against `C \ N`.
    pub first_pass_checks: Vec<(Vec<i64>, NilpotentRayDivergenceCheck)>,
    /// Finite-table comparison charges for the second pass, i.e. `C \ F0`.
    pub second_pass_comparison_charges: Vec<(Vec<i64>, Integer)>,
    /// Second-pass checks against `C \ F0`.
    pub second_pass_checks: Vec<(Vec<i64>, NilpotentRayDivergenceCheck)>,
    /// Provisional `F0`, final `F`, and inconclusive rays.
    pub classification: NilpotentRayTwoPassNopClassification,
}

/// Exact certificate that an integer normal cuts out a supporting Mori face.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SupportingMoriFaceCertificate {
    /// Integer normal with non-negative pairing on every Mori generator.
    pub normal: Vec<i64>,
    /// Number of Mori generators with zero pairing against `normal`.
    pub zero_generator_count: usize,
    /// Number of Mori generators with positive pairing against `normal`.
    pub positive_generator_count: usize,
}

/// Search controls for LP-assisted supporting Mori face certificates.
///
/// The LP only proposes candidate normals. Any returned certificate is verified
/// exactly with integer arithmetic by [`check_supporting_mori_face_normal`].
/// A failed search is inconclusive: it does not prove that no supporting normal
/// exists.
#[derive(Clone, Debug, PartialEq)]
pub struct SupportingMoriFaceLpSearchOptions {
    /// Maximum number of anchor rays to try after the aggregate-normal pass.
    pub anchor_attempts: usize,
    /// Maximum cutting-plane rounds used to add violated Mori inequalities.
    pub cutting_rounds: usize,
    /// Largest integer scale used when rounding an LP normal.
    pub scale_limit: i64,
    /// Symmetric bound for each LP normal coordinate.
    pub variable_bound: f64,
}

impl Default for SupportingMoriFaceLpSearchOptions {
    fn default() -> Self {
        Self {
            anchor_attempts: 16,
            cutting_rounds: 64,
            scale_limit: 100_000,
            variable_bound: 1.0e9,
        }
    }
}

/// Diagnostic trace for an LP-assisted supporting Mori face search.
///
/// The LP phases only propose real normals. A populated `certificate` has
/// still passed exact integer verification through
/// [`check_supporting_mori_face_normal`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SupportingMoriFaceLpSearchDiagnostic {
    /// Rank of the supplied face generators.
    pub face_rank: usize,
    /// Ambient curve-coordinate dimension.
    pub dim: usize,
    /// Final search status.
    pub status: String,
    /// Outcome of the exact-kernel certificate attempt.
    pub exact_kernel_status: String,
    /// Outcome of the full LP feasibility attempt with all Mori inequalities.
    pub full_status: String,
    /// Outcome of the aggregate LP attempt.
    pub aggregate_status: String,
    /// Number of anchor rays tested after the aggregate LP attempt.
    pub anchor_attempt_count: usize,
    /// Number of anchor LPs that returned a finite real normal.
    pub anchor_lp_solution_count: usize,
    /// Counts of per-anchor outcomes.
    pub anchor_status_counts: BTreeMap<String, usize>,
    /// Exact certificate, when any search phase found one.
    pub certificate: Option<SupportingMoriFaceCertificate>,
}

/// Exact certificate that a target curve spans an extremal Mori ray.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExtremalMoriRayCertificate {
    /// Integer separator normal. It pairs negatively with the target curve and
    /// non-negatively with every Mori generator not on the same positive ray.
    pub separator_normal: Vec<i64>,
    /// Number of Mori generators on the same positive rational ray as the
    /// target curve.
    pub same_ray_generator_count: usize,
    /// Number of non-same-ray Mori generators with zero separator pairing.
    pub zero_other_generator_count: usize,
    /// Number of non-same-ray Mori generators with positive separator pairing.
    pub positive_other_generator_count: usize,
}

/// Diagnostic for LP-assisted extremal Mori ray separator search.
///
/// The LP phase only proposes a real normal. `certificate` is populated only
/// when a rounded integer normal passes exact separator verification.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExtremalMoriRayLpSearchDiagnostic {
    /// Search outcome.
    pub status: String,
    /// Whether the LP solver returned a finite real normal.
    pub lp_solution_found: bool,
    /// Number of distinct rounded integer normals tested exactly.
    pub exact_normal_candidate_count: usize,
    /// Exact certificate, if one was found.
    pub certificate: Option<ExtremalMoriRayCertificate>,
}

/// Mori generators lying on an exactly certified supporting face.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SupportingMoriFace {
    /// Exact certificate for the supporting normal.
    pub certificate: SupportingMoriFaceCertificate,
    /// Mori generators with zero pairing against the supporting normal.
    pub generators: Vec<Vec<i64>>,
}

/// One term in a finite semigroup decomposition of a curve class.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CurveDecompositionTerm {
    /// Summand curve class.
    pub class: Vec<i64>,
    /// Non-negative integer multiplicity of this summand.
    pub multiplicity: u64,
}

/// Pruning rule for selected toric curve candidates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CurvePruningStrategy {
    /// Remove curves that are sums of two selected candidates.
    ///
    /// This is the rule that reproduces the McAllister `small_curves.dat`
    /// checkpoint for 4-214-647.
    PairDecomposable,
    /// Remove curves that are finite-semigroup sums of selected candidates.
    ///
    /// This is stricter than the McAllister checkpoint rule and is intended for
    /// GA/search runs that want a finite selected-set Hilbert-basis proxy.
    FiniteSemigroup,
}

impl CurvePruningStrategy {
    /// Stable CLI/log label for this pruning rule.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::PairDecomposable => "pair",
            Self::FiniteSemigroup => "finite-semigroup",
        }
    }
}

/// A toric curve class with its genus-zero GV invariant.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ToricCurveGvInvariant {
    /// Curve class in ambient divisor-intersection coordinates.
    pub class: Vec<i64>,
    /// Genus-zero GV invariant for this toric curve class.
    pub gv: Integer,
}

/// A nonzero point in an affine toric circuit relation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AffineCircuitRelationPoint {
    /// Index of the triangulation point.
    pub point_index: usize,
    /// Coefficient of the point in the affine relation.
    pub coefficient: i64,
    /// Lattice coordinates of the point.
    pub coordinates: Vec<i64>,
}

/// Two-dimensional coordinates for a point in a reconstructed local toric model.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalToricCoordinate2D {
    /// Index of the triangulation point in the ambient point list.
    pub point_index: usize,
    /// Coordinates in a deterministic rank-two local lattice basis.
    pub coordinates: [i64; 2],
}

/// Coordinates for a point in a reconstructed local affine toric lattice.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LocalToricCoordinate {
    /// Index of the triangulation point in the ambient point list.
    pub point_index: usize,
    /// Coordinates in a deterministic local lattice basis of the affine support.
    pub coordinates: Vec<i64>,
}

/// One point in a normalized rank-two local support signature.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RankTwoLocalSupportSignatureEntry {
    /// Affine relation coefficient at this local point.
    pub coefficient: i64,
    /// Translated local coordinate. The signature translation places the
    /// lexicographically smallest local point at the origin.
    pub coordinates: [i64; 2],
}

/// Point-index-free signature of a rank-two local toric support.
///
/// This deliberately does not assign any GV value. It is a grouping and audit
/// object for reconstructing the local toric model from upstream lattice
/// points before any mirror-symmetry computation is attempted.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RankTwoLocalSupportSignature {
    /// Coefficient-coordinate entries, sorted after translation and relation
    /// orientation normalization.
    pub entries: Vec<RankTwoLocalSupportSignatureEntry>,
}

/// One point in the canonical local charge model for a rank-two support.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RankTwoLocalChargeModelPoint {
    /// Coefficient of this point in the target potent-ray affine relation.
    pub relation_coefficient: i64,
    /// Local integral coordinate of this support point.
    pub coordinates: [i64; 2],
}

/// Point-index-free local toric charge model for a rank-two support.
///
/// This is the source-derived input object before any local mirror/HKTY series
/// is assigned: the canonical support points, the target relation in that point
/// order, and an integer basis for the affine charge lattice of the local
/// toric diagram.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RankTwoLocalChargeModel {
    /// Canonical support points. The order is the same as
    /// [`target_relation`](Self::target_relation) and the columns of
    /// [`charge_basis`](Self::charge_basis).
    pub points: Vec<RankTwoLocalChargeModelPoint>,
    /// The target affine relation coefficients in canonical point order.
    pub target_relation: Vec<i64>,
    /// Integer coordinates of [`target_relation`](Self::target_relation) in
    /// [`charge_basis`](Self::charge_basis).
    pub target_relation_in_charge_basis: Vec<i64>,
    /// Integer basis for the kernel of `[1; x; y]`, in canonical point order.
    pub charge_basis: Vec<Vec<i64>>,
}

/// CKYZ local toric surface source model matched to a reconstructed support.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CkyzLocalSurfaceKind {
    /// Canonical bundle over `P^2`.
    LocalP2,
    /// Canonical bundle over `F0 = P^1 x P^1`.
    HirzebruchF0,
    /// Canonical bundle over the first Hirzebruch surface.
    HirzebruchF1,
    /// CKYZ two-dimensional reflexive polygon 5.
    Polygon5,
}

/// One term in CKYZ's local base intersection expression.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CkyzLocalIntersectionTerm {
    /// Zero-based first divisor index.
    pub first: usize,
    /// Zero-based second divisor index.
    pub second: usize,
    /// Integral coefficient of `J_first J_second`.
    pub coefficient: i64,
}

/// Exact identification of a local charge model with CKYZ source relation data.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CkyzLocalSurfaceIdentification {
    /// CKYZ source surface.
    pub kind: CkyzLocalSurfaceKind,
    /// Column permutation putting Cyrus' canonical point order into CKYZ order.
    pub point_permutation: Vec<usize>,
    /// Unimodular row transform mapping the permuted Cyrus charge basis to CKYZ
    /// source relation rows.
    pub row_transform: Vec<Vec<i64>>,
    /// CKYZ source relation rows.
    pub source_relations: Vec<Vec<i64>>,
    /// Target potent-ray direction in the CKYZ source relation basis.
    pub source_target_direction: Vec<i64>,
    /// CKYZ `c1(B)` expression coefficients in source divisor coordinates.
    pub c1_coefficients: Vec<i64>,
    /// CKYZ local base intersection expression in source coordinates.
    pub local_intersection_terms: Vec<CkyzLocalIntersectionTerm>,
}

/// Finite local semigroup context used by CKYZ causal-domain extraction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CkyzLocalCausalDomainSpec {
    /// Effective generators in CKYZ source relation coordinates.
    pub generators: Vec<Vec<usize>>,
    /// Positive grading vector used to order and truncate the generated
    /// semigroup.
    pub grading_vector: Vec<usize>,
}

/// Size profile for CKYZ local finite monomial domains.
///
/// This is a source-auditing object: it does not assign GV values. It compares
/// the domains used by the local CKYZ extraction so we can see whether a
/// targeted computation is spending time in the componentwise target downset,
/// the support-predicted mirror-map/potential history, or a generated causal
/// semigroup domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CkyzLocalDomainProfile {
    /// Number of CKYZ flat coordinates.
    pub rank: usize,
    /// Number of requested degrees after adding primitive cover divisors.
    pub extraction_degree_count: usize,
    /// Largest total degree among the originally requested target degrees.
    pub max_target_total_degree: usize,
    /// Number of monomials in the cover-closed componentwise downset.
    pub target_downset_degree_count: usize,
    /// Whether the target downset built the fast addition table.
    pub target_downset_has_addition_table: bool,
    /// Number of monomials retained by the support-predicted path-history
    /// domain.
    pub predicted_support_degree_count: usize,
    /// Largest total degree retained by the support-predicted domain.
    pub predicted_support_max_total_degree: usize,
    /// Whether the support-predicted domain built the fast addition table.
    pub predicted_support_has_addition_table: bool,
    /// Number of monomials in an optional generated causal semigroup domain.
    pub causal_semigroup_degree_count: Option<usize>,
    /// Largest total degree retained by an optional causal semigroup domain.
    pub causal_semigroup_max_total_degree: Option<usize>,
    /// Whether an optional causal semigroup domain built the fast addition
    /// table.
    pub causal_semigroup_has_addition_table: Option<bool>,
}

/// Coefficient-work profile for CKYZ z-residual extraction.
///
/// This is a source-auditing object. It does not compute GV values; it counts
/// the degree/scale states that the z-domain residual subtraction would need to
/// inspect before rational coefficient arithmetic.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CkyzZResidualCoefficientWorkProfile {
    /// Number of CKYZ flat coordinates.
    pub rank: usize,
    /// Size of the finite monomial domain used for alpha/potential data.
    pub domain_degree_count: usize,
    /// Number of selected residual-history degrees.
    pub path_history_degree_count: usize,
    /// Number of ordered lower-grading residual degree pairs inspected by
    /// extraction.
    pub residual_pair_count: usize,
    /// Number of same-grading index pairs skipped by the degree-ordered batch order.
    pub same_grading_pair_skip_count: usize,
    /// Number of residual pairs that pass the componentwise divisibility gate.
    pub componentwise_pair_count: usize,
    /// Number of multiple-cover delta terms considered across those pairs.
    pub li2_delta_term_count: usize,
    /// Number of residual pairs whose Li2 support can affect the target.
    pub support_pair_count: usize,
    /// Number of multiple-cover delta terms that survive Li2 support pruning.
    pub support_li2_delta_term_count: usize,
    /// Number of distinct scale degrees `m*d` appearing in those terms.
    pub unique_scale_count: usize,
    /// Number of distinct delta degrees `target - m*d` appearing in those terms.
    pub unique_delta_count: usize,
    /// Number of distinct `(scale, delta)` exponential coefficient states.
    pub unique_exp_state_count: usize,
    /// Number of distinct `(scale, delta)` states that survive support pruning.
    pub support_unique_exp_state_count: usize,
    /// Supported `(scale, delta)` counts grouped by scale degree, sorted largest first.
    pub support_exp_state_counts_by_scale: Vec<(Vec<usize>, usize)>,
    /// Number of rolling history levels cygv would retain for this rank.
    pub qn_history_level_count: usize,
    /// Number of selected history degrees whose closest source-shaped `q_N`
    /// start can reuse a previous candidate.
    pub qn_history_candidate_hit_count: usize,
    /// Number of selected history degrees that would still start from the
    /// direct monomial in the source-shaped `q_N` history heuristic.
    pub qn_history_candidate_miss_count: usize,
    /// Number of distinct nonzero delta degrees selected by the source-shaped
    /// `q_N` history heuristic.
    pub qn_history_unique_delta_count: usize,
    /// Restricted-domain sizes for the distinct nonzero `q_delta` degrees
    /// selected by the source-shaped `q_N` history heuristic.
    pub q_delta_domain_profiles: Vec<CkyzQDeltaDomainProfile>,
}

/// Size profile for a restricted `q_delta = z^delta exp(delta.alpha)` domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CkyzQDeltaDomainProfile {
    /// Nonzero delta degree selected by the previous-`q_N` heuristic.
    pub delta_degree: Vec<usize>,
    /// Number of exponential monomials `e` for which `delta + e` is present in
    /// the finite CKYZ domain.
    pub shiftable_exp_degree_count: usize,
    /// Number of monomials needed after closing those shiftable output
    /// exponents under alpha-term predecessors.
    pub predecessor_closure_degree_count: usize,
    /// Whether that closure is small enough to build Cyrus' dense addition
    /// table for indexed finite-polynomial products.
    pub predecessor_closure_has_addition_table: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CkyzQnHistoryReuseProfile {
    level_count: usize,
    candidate_hit_count: usize,
    candidate_miss_count: usize,
    unique_delta_count: usize,
    unique_delta_degrees: Vec<Vec<usize>>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CkyzLocalSurfaceSource {
    kind: CkyzLocalSurfaceKind,
    relations: Vec<Vec<i64>>,
    c1_coefficients: Vec<i64>,
    local_intersection_terms: Vec<CkyzLocalIntersectionTerm>,
}

/// Recognized local toric circuit shape.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LocalToricCircuitKind {
    /// Local `P^2` / `O(-3) -> P^2` triangle relation.
    LocalP2Triangle {
        /// The interior point of the toric triangle.
        interior_point: usize,
        /// The three triangle vertices.
        vertex_points: Vec<usize>,
        /// Coefficient of the interior point in the supplied orientation.
        interior_coefficient: i64,
        /// Coefficient of each vertex in the supplied orientation.
        vertex_coefficient: i64,
        /// Local rank-two coordinates for the triangle support.
        local_coordinates: Vec<LocalToricCoordinate2D>,
    },
}

/// Diagnostic for an ambient curve row that is an affine toric circuit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AffineToricCircuitDiagnostic {
    /// Nonzero relation points in ambient index order.
    pub relation_points: Vec<AffineCircuitRelationPoint>,
    /// Affine rank of the support points.
    pub affine_rank: usize,
    /// Counts of nonzero coefficients in the affine relation.
    pub coefficient_counts: BTreeMap<i64, usize>,
    /// Sum of the relation coefficients. This is zero for an affine relation.
    pub coefficient_sum: i128,
    /// Weighted sum of lattice coordinates. This is zero for an affine relation.
    pub coordinate_sum: Vec<i128>,
    /// Integer basis for affine charge relations among the support points.
    ///
    /// Rows are expressed in the same support order as `relation_points`. This
    /// is the local toric charge context that must be used before assigning any
    /// local mirror-symmetry GV values.
    pub local_charge_basis: Vec<Vec<i64>>,
    /// Integral coordinates for the local affine support in its own lattice.
    ///
    /// The coordinates are computed from the affine lattice generated by the
    /// support points, not from any coefficient-pattern lookup. They are
    /// available for every affine rank.
    pub local_coordinates: Vec<LocalToricCoordinate>,
    /// Integral coordinates for rank-two local toric supports.
    ///
    /// The coordinates are computed from the affine lattice generated by the
    /// support points, not from any coefficient-pattern lookup. They are
    /// present only when the support has affine rank two.
    pub local_coordinates_2d: Option<Vec<LocalToricCoordinate2D>>,
    /// Recognized local toric shape, when currently known.
    pub kind: Option<LocalToricCircuitKind>,
}

/// A nonzero point in an origin-circuit affine relation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OriginCircuitRelationPoint {
    /// Index of the triangulation point.
    pub point_index: usize,
    /// Coefficient of the point in the affine relation.
    pub coefficient: i64,
    /// Lattice coordinates of the point.
    pub coordinates: Vec<i64>,
    /// Dimension of the smallest primal polytope face containing the point.
    pub face_dimension: Option<usize>,
}

/// A triangulation witness for an origin-circuit Mori-cap curve class.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OriginCircuitCurveWitness {
    /// Curve class in ambient divisor-intersection coordinates.
    pub class: Vec<i64>,
    /// First point outside the shared facet intersection used in the circuit.
    pub first_facet_exclusive_point: usize,
    /// Second point outside the shared facet intersection used in the circuit.
    pub second_facet_exclusive_point: usize,
    /// Two-dimensional simplex shared by the two facets.
    pub shared_two_simplex: Vec<usize>,
    /// First facet containing `shared_two_simplex`.
    pub first_facet: Vec<usize>,
    /// Second facet containing `shared_two_simplex`.
    pub second_facet: Vec<usize>,
    /// Sparse affine relation before normalization.
    pub sparse_relation: Vec<(usize, i64)>,
    /// Nonzero relation points with coordinates and face dimensions.
    pub relation_points: Vec<OriginCircuitRelationPoint>,
}

/// An origin-circuit Mori-cap curve class, with enough shape data to diagnose
/// whether a local toric GV formula is currently known.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OriginCircuitCurveDiagnostic {
    /// Curve class in ambient divisor-intersection coordinates.
    pub class: Vec<i64>,
    /// Coefficient of the origin in `class`.
    pub origin_coefficient: i64,
    /// Counts of negative non-origin coefficients by coefficient value.
    pub negative_coefficient_counts: BTreeMap<i64, usize>,
    /// Counts of positive non-origin coefficients by coefficient value.
    pub positive_coefficient_counts: BTreeMap<i64, usize>,
    /// Whether this is the isolated resolved-conifold pattern currently covered
    /// by [`compute_toric_two_face_curve_gv_invariants`].
    pub is_resolved_conifold_pattern: bool,
    /// Triangulation witnesses that produce this origin-circuit class.
    pub witnesses: Vec<OriginCircuitCurveWitness>,
}

/// The local formula source used for a toric curve GV value.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ToricCurveGvSource {
    /// Curve from a circuit in a primal two-face.
    TwoFace {
        /// Edge shared by the two adjacent two-face simplices.
        edge: Vec<usize>,
        /// All triangulation-point indices in the two-face.
        two_face_points: Vec<usize>,
        /// Genus of the dual one-face.
        two_face_genus: usize,
        /// Curve coefficients on the shared edge endpoints.
        edge_coefficients: (i64, i64),
        /// Smallest primal face dimensions for the shared edge endpoints.
        edge_face_dimensions: (Option<usize>, Option<usize>),
        /// Dual two-face genera for edge endpoints that are primal one-face points.
        edge_one_face_genera: (Option<usize>, Option<usize>),
    },
    /// Origin-circuit curve matching the isolated resolved-conifold pattern.
    ResolvedConifoldOriginCircuit {
        /// Index of the origin point in ambient curve coordinates.
        origin_index: usize,
        /// Triangulation witness that produced this curve class.
        witness: OriginCircuitCurveWitness,
    },
}

/// A toric curve class with its GV value and local formula provenance.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ToricCurveGvDiagnostic {
    /// Curve class in ambient divisor-intersection coordinates.
    pub class: Vec<i64>,
    /// Genus-zero GV invariant for this toric curve class.
    pub gv: Integer,
    /// Local formulas/witnesses that produced the value.
    pub sources: Vec<ToricCurveGvSource>,
}

/// Diagnose whether an ambient curve row is an affine toric circuit.
///
/// The input curve uses the same ambient coordinate convention as
/// McAllister's `potent_rays.dat`: one coefficient per triangulation point.
/// This function only recognizes upstream toric geometry in the row; it does
/// not assign any GV values.
///
/// # Errors
/// Returns an error if the row length does not match the point list or if the
/// points do not all live in the same lattice dimension.
pub fn diagnose_affine_toric_circuit(
    ambient_curve: &[i64],
    points: &[Point],
) -> Result<Option<AffineToricCircuitDiagnostic>> {
    if ambient_curve.len() != points.len() {
        return Err(Error::InvalidInput(format!(
            "ambient curve dimension {} does not match point count {}",
            ambient_curve.len(),
            points.len()
        )));
    }
    let Some(first_point) = points.first() else {
        return Ok(None);
    };
    let dim = first_point.coords().len();
    if points.iter().any(|point| point.coords().len() != dim) {
        return Err(Error::InvalidInput(
            "affine circuit points must have consistent lattice dimension".into(),
        ));
    }

    let mut relation_points = Vec::new();
    let mut coefficient_counts = BTreeMap::new();
    let mut coefficient_sum = 0i128;
    let mut coordinate_sum = vec![0i128; dim];
    for (point_index, (&coefficient, point)) in ambient_curve.iter().zip(points.iter()).enumerate()
    {
        if coefficient == 0 {
            continue;
        }
        *coefficient_counts.entry(coefficient).or_insert(0) += 1;
        coefficient_sum += i128::from(coefficient);
        for (acc, &coord) in coordinate_sum.iter_mut().zip(point.coords().iter()) {
            *acc += i128::from(coefficient) * i128::from(coord);
        }
        relation_points.push(AffineCircuitRelationPoint {
            point_index,
            coefficient,
            coordinates: point.coords().to_vec(),
        });
    }

    if relation_points.is_empty() {
        return Ok(None);
    }
    if coefficient_sum != 0 || coordinate_sum.iter().any(|&value| value != 0) {
        return Ok(None);
    }

    let affine_rank = affine_relation_rank(&relation_points);
    let local_charge_basis = local_affine_charge_basis(&relation_points)?;
    let local_coordinates = local_affine_coordinates(&relation_points, affine_rank)?;
    let local_coordinates_2d = local_rank_two_coordinates(&relation_points, affine_rank)?;
    let kind = classify_local_toric_circuit(&relation_points, affine_rank);
    Ok(Some(AffineToricCircuitDiagnostic {
        relation_points,
        affine_rank,
        coefficient_counts,
        coefficient_sum,
        coordinate_sum,
        local_charge_basis,
        local_coordinates,
        local_coordinates_2d,
        kind,
    }))
}

/// Build a stable point-index-free signature for a rank-two affine support.
///
/// The signature uses the reconstructed rank-two local coordinates, translates
/// the local diagram so the lexicographically smallest coordinate is at the
/// origin, sorts coefficient-coordinate entries, and normalizes the overall
/// relation sign by choosing the lexicographically smaller of `q` and `-q`.
#[must_use]
pub fn rank_two_local_support_signature(
    diagnostic: &AffineToricCircuitDiagnostic,
) -> Option<RankTwoLocalSupportSignature> {
    let coordinates = diagnostic.local_coordinates_2d.as_ref()?;
    let coordinates_by_point: BTreeMap<usize, [i64; 2]> = coordinates
        .iter()
        .map(|point| (point.point_index, point.coordinates))
        .collect();
    let entries = diagnostic
        .relation_points
        .iter()
        .map(|point| {
            let coordinates = coordinates_by_point.get(&point.point_index).copied()?;
            Some(RankTwoLocalSupportSignatureEntry {
                coefficient: point.coefficient,
                coordinates,
            })
        })
        .collect::<Option<Vec<_>>>()?;

    let entries = canonical_rank_two_signature_entries(&entries)?;
    Some(RankTwoLocalSupportSignature { entries })
}

/// Build the local toric charge model for a normalized rank-two support.
///
/// This reconstructs the integer kernel of the canonical matrix with rows
/// `[1; x; y]`. It also verifies that the target potent-ray relation belongs
/// to that reconstructed affine charge lattice. No GV value or local mirror
/// series is assigned here.
pub fn rank_two_local_charge_model(
    signature: &RankTwoLocalSupportSignature,
) -> Result<RankTwoLocalChargeModel> {
    if signature.entries.is_empty() {
        return Err(Error::InvalidInput(
            "rank-two local charge model requires at least one support point".into(),
        ));
    }

    let mut points = signature
        .entries
        .iter()
        .map(|entry| RankTwoLocalChargeModelPoint {
            relation_coefficient: entry.coefficient,
            coordinates: entry.coordinates,
        })
        .collect::<Vec<_>>();
    points.sort_unstable();

    let support_len = points.len();
    let mut matrix = vec![vec![Integer::from(0); support_len]; 3];
    for (col, point) in points.iter().enumerate() {
        matrix[0][col] = Integer::from(1);
        matrix[1][col] = Integer::from(point.coordinates[0]);
        matrix[2][col] = Integer::from(point.coordinates[1]);
    }

    let mut charge_basis = integer_kernel(&matrix)
        .into_iter()
        .map(|row| {
            let mut converted = row
                .iter()
                .map(|value| {
                    i64::try_from(value).map_err(|_| {
                        Error::InvalidInput(
                            "rank-two local charge-model basis entry does not fit in i64".into(),
                        )
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            normalize_relation_orientation(&mut converted);
            Ok(converted)
        })
        .collect::<Result<Vec<_>>>()?;
    charge_basis.sort_unstable();

    let target_relation = points
        .iter()
        .map(|point| point.relation_coefficient)
        .collect::<Vec<_>>();
    if !curve_in_rational_row_span(&target_relation, &charge_basis)? {
        return Err(Error::InvalidInput(
            "rank-two local target relation is not in the reconstructed charge lattice".into(),
        ));
    }
    let target_relation_in_charge_basis =
        relation_coordinates_in_row_basis(&target_relation, &charge_basis)?;

    Ok(RankTwoLocalChargeModel {
        points,
        target_relation,
        target_relation_in_charge_basis,
        charge_basis,
    })
}

/// Identify a rank-two local charge model with CKYZ local mirror source data.
///
/// The match is exact and pre-GV: Cyrus searches point permutations and
/// unimodular row transformations until the reconstructed local charge basis
/// equals one of the CKYZ relation systems. The returned target direction is
/// expressed in that CKYZ source basis. This function does not read or assign
/// any GV sequence.
pub fn identify_ckyz_local_surface(
    model: &RankTwoLocalChargeModel,
) -> Result<Option<CkyzLocalSurfaceIdentification>> {
    for source in ckyz_local_surface_sources() {
        if source.relations.is_empty() || source.relations[0].len() != model.points.len() {
            continue;
        }
        if source.relations.len() != model.charge_basis.len() {
            continue;
        }

        for point_permutation in permutations(model.points.len()) {
            let permuted_charge_basis =
                permute_matrix_columns(&model.charge_basis, &point_permutation);
            let Some(row_transform) =
                integer_row_transform_between_bases(&source.relations, &permuted_charge_basis)?
            else {
                continue;
            };
            let permuted_target = permute_vector_values(&model.target_relation, &point_permutation);
            let Ok(source_target_direction) =
                relation_coordinates_in_row_basis(&permuted_target, &source.relations)
            else {
                continue;
            };
            return Ok(Some(CkyzLocalSurfaceIdentification {
                kind: source.kind,
                point_permutation,
                row_transform,
                source_relations: source.relations,
                source_target_direction,
                c1_coefficients: source.c1_coefficients,
                local_intersection_terms: source.local_intersection_terms,
            }));
        }
    }

    Ok(None)
}

/// Return the finite-limit cover-weight coefficients for a CKYZ local surface.
///
/// These are the coefficients `x_i` used by
/// [`extract_ckyz_local_gv_invariants_from_potential`] in
/// `w(d) = -sum_i x_i d_i`. They are not always the printed CKYZ `C1`
/// coefficients: for example `F1` uses `[2, 1]`, and polygon 5 uses
/// `[1, 1, 1]`.
#[must_use]
pub fn ckyz_local_surface_cover_weight_coefficients(kind: &CkyzLocalSurfaceKind) -> &'static [i64] {
    match kind {
        CkyzLocalSurfaceKind::LocalP2 => &[3],
        CkyzLocalSurfaceKind::HirzebruchF0 => &[2, 2],
        CkyzLocalSurfaceKind::HirzebruchF1 => &[2, 1],
        CkyzLocalSurfaceKind::Polygon5 => &[1, 1, 1],
    }
}

/// Build the CKYZ source-coordinate semigroup context for a matched local
/// surface.
///
/// The generators are the CKYZ Mori-coordinate unit classes for the matched
/// relation basis. The grading uses the same finite-limit cover weights as GV
/// extraction, so the causal-domain truncation is tied to the source formula
/// rather than to a caller-invented total-degree box.
///
/// # Errors
/// Returns an error if the matched source has inconsistent rank data or if a
/// finite-limit cover weight is not strictly positive.
pub fn ckyz_local_surface_causal_domain_spec(
    identification: &CkyzLocalSurfaceIdentification,
) -> Result<CkyzLocalCausalDomainSpec> {
    let rank = identification.source_relations.len();
    if rank == 0 {
        return Err(Error::InvalidInput(
            "CKYZ local causal domain requires at least one source relation".into(),
        ));
    }

    let cover_weights = ckyz_local_surface_cover_weight_coefficients(&identification.kind);
    if cover_weights.len() != rank {
        return Err(Error::InvalidInput(
            "CKYZ local causal-domain cover-weight rank does not match source rank".into(),
        ));
    }

    let mut grading_vector = Vec::with_capacity(rank);
    for &weight in cover_weights {
        if weight <= 0 {
            return Err(Error::InvalidInput(
                "CKYZ local causal-domain grading weights must be positive".into(),
            ));
        }
        grading_vector.push(usize::try_from(weight).map_err(|_| {
            Error::InvalidInput("CKYZ local causal-domain grading weight does not fit usize".into())
        })?);
    }

    let generators = (0..rank)
        .map(|index| {
            let mut generator = vec![0; rank];
            generator[index] = 1;
            generator
        })
        .collect();

    Ok(CkyzLocalCausalDomainSpec {
        generators,
        grading_vector,
    })
}

/// Return the positive multiples of a matched CKYZ potent-ray target direction.
///
/// The target direction is expressed in the CKYZ source relation basis by
/// [`identify_ckyz_local_surface`].
///
/// # Errors
/// Returns an error if `multiples` is zero, if the target direction has the
/// wrong rank, if any coordinate is negative, or if scaling overflows.
pub fn ckyz_local_surface_target_degrees(
    identification: &CkyzLocalSurfaceIdentification,
    multiples: usize,
) -> Result<Vec<Vec<usize>>> {
    if multiples == 0 {
        return Err(Error::InvalidInput(
            "CKYZ local target-degree multiple count must be positive".into(),
        ));
    }
    let rank = identification.source_relations.len();
    if identification.source_target_direction.len() != rank {
        return Err(Error::InvalidInput(
            "CKYZ local target direction rank does not match source rank".into(),
        ));
    }

    let direction = identification
        .source_target_direction
        .iter()
        .map(|&entry| {
            if entry < 0 {
                return Err(Error::InvalidInput(
                    "CKYZ local target direction must be nonnegative".into(),
                ));
            }
            usize::try_from(entry).map_err(|_| {
                Error::InvalidInput("CKYZ local target direction entry does not fit usize".into())
            })
        })
        .collect::<Result<Vec<_>>>()?;
    if direction.iter().all(|&entry| entry == 0) {
        return Err(Error::InvalidInput(
            "CKYZ local target direction must be nonzero".into(),
        ));
    }

    (1..=multiples)
        .map(|multiple| {
            direction
                .iter()
                .map(|&entry| {
                    entry.checked_mul(multiple).ok_or_else(|| {
                        Error::InvalidInput("CKYZ local target-degree scaling overflowed".into())
                    })
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect()
}

/// Profile the CKYZ finite monomial domains for requested local degrees.
///
/// This is intended for source audits and performance work. It does not compute
/// any GV invariant; it only compares the componentwise cover-closed downset
/// with the support-predicted mirror-map/potential history domain.
///
/// # Errors
/// Returns an error for invalid CKYZ relation rows, invalid local intersection
/// terms, invalid target degrees, or domain-construction overflow.
pub fn ckyz_local_domain_profile_for_degrees(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    target_degrees: &[Vec<usize>],
) -> Result<CkyzLocalDomainProfile> {
    ckyz_local_domain_profile_for_degrees_impl(
        relations,
        local_intersection_terms,
        target_degrees,
        None,
    )
}

/// Profile CKYZ finite monomial domains, including an explicit causal
/// semigroup domain.
///
/// This mirrors [`compute_ckyz_local_gv_invariants_for_degrees_with_causal_domain`]
/// at the domain-construction level, while also reporting the componentwise and
/// support-predicted domains for comparison.
///
/// # Errors
/// Returns an error for invalid CKYZ source data, invalid causal generators,
/// non-generated target degrees, or domain-construction overflow.
pub fn ckyz_local_domain_profile_for_degrees_with_causal_domain(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    target_degrees: &[Vec<usize>],
    causal_generators: &[Vec<usize>],
    grading_vector: &[usize],
) -> Result<CkyzLocalDomainProfile> {
    ckyz_local_domain_profile_for_degrees_impl(
        relations,
        local_intersection_terms,
        target_degrees,
        Some((causal_generators, grading_vector)),
    )
}

/// Profile CKYZ finite monomial domains for multiples of a matched local
/// surface target direction.
///
/// The causal semigroup domain is built from
/// [`ckyz_local_surface_causal_domain_spec`], so the report is tied to the same
/// source-derived grading used by local-surface causal GV extraction.
///
/// # Errors
/// Returns an error for invalid matched source data, invalid target multiples,
/// or domain-construction failures.
pub fn ckyz_local_surface_domain_profile_for_multiples(
    identification: &CkyzLocalSurfaceIdentification,
    multiples: usize,
) -> Result<CkyzLocalDomainProfile> {
    let target_degrees = ckyz_local_surface_target_degrees(identification, multiples)?;
    let domain_spec = ckyz_local_surface_causal_domain_spec(identification)?;
    ckyz_local_domain_profile_for_degrees_with_causal_domain(
        &identification.source_relations,
        &identification.local_intersection_terms,
        &target_degrees,
        &domain_spec.generators,
        &domain_spec.grading_vector,
    )
}

fn ckyz_local_domain_profile_for_degrees_impl(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    target_degrees: &[Vec<usize>],
    causal_domain: Option<(&[Vec<usize>], &[usize])>,
) -> Result<CkyzLocalDomainProfile> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    validate_ckyz_target_degrees(target_degrees, rank)?;
    for term in local_intersection_terms {
        if term.first >= rank || term.second >= rank {
            return Err(Error::InvalidInput(
                "CKYZ local intersection term index is outside the relation rank".into(),
            ));
        }
    }

    let extraction_degrees = ckyz_cover_closed_target_degrees(target_degrees)?;
    let max_target_total_degree = target_degrees
        .iter()
        .map(|degree| ckyz_total_degree(degree))
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .max()
        .expect("validated target degrees are nonempty");
    let target_downset = CkyzMonomialDomain::target_downset(&extraction_degrees, rank)?;
    let predicted_support = ckyz_predicted_support_domain_for_degrees(
        relations,
        local_intersection_terms,
        target_degrees,
    )?;
    let causal_semigroup = causal_domain
        .map(|(generators, grading_vector)| {
            ckyz_causal_monomial_domain(rank, generators, grading_vector, &extraction_degrees)
        })
        .transpose()?;

    Ok(CkyzLocalDomainProfile {
        rank,
        extraction_degree_count: extraction_degrees.len(),
        max_target_total_degree,
        target_downset_degree_count: target_downset.degrees.len(),
        target_downset_has_addition_table: target_downset.addition_indices.is_some(),
        predicted_support_degree_count: predicted_support.degrees.len(),
        predicted_support_max_total_degree: predicted_support.max_total_degree,
        predicted_support_has_addition_table: predicted_support.addition_indices.is_some(),
        causal_semigroup_degree_count: causal_semigroup.as_ref().map(|domain| domain.degrees.len()),
        causal_semigroup_max_total_degree: causal_semigroup
            .as_ref()
            .map(|domain| domain.max_total_degree),
        causal_semigroup_has_addition_table: causal_semigroup
            .as_ref()
            .map(|domain| domain.addition_indices.is_some()),
    })
}

/// Profile the coefficient states implied by CKYZ z-residual extraction.
///
/// This computes the support-predicted finite domain and selected z-residual
/// history, then counts the formal `Li2(q^d exp(d.alpha))` scale/delta states
/// that may need coefficient extraction. It intentionally stops before rational
/// coefficient arithmetic and does not read or compare any GV table.
///
/// # Errors
/// Returns an error for invalid CKYZ source data, invalid target degrees, or
/// domain-construction overflow.
pub fn ckyz_z_residual_coefficient_work_profile_for_degrees(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    cover_weight_coefficients: &[i64],
    target_degrees: &[Vec<usize>],
) -> Result<CkyzZResidualCoefficientWorkProfile> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    validate_ckyz_target_degrees(target_degrees, rank)?;
    if cover_weight_coefficients.len() != rank {
        return Err(Error::InvalidInput(
            "CKYZ coefficient profile cover-weight rank does not match relation rank".into(),
        ));
    }
    for term in local_intersection_terms {
        if term.first >= rank || term.second >= rank {
            return Err(Error::InvalidInput(
                "CKYZ local intersection term index is outside the relation rank".into(),
            ));
        }
    }

    let extraction_degrees = ckyz_cover_closed_target_degrees(target_degrees)?;
    let domain = ckyz_predicted_support_domain_for_degrees(
        relations,
        local_intersection_terms,
        &extraction_degrees,
    )?;
    let alpha = compute_ckyz_log_period_corrections_domain(relations, &domain)?;
    let path_history_degrees =
        ckyz_z_residual_dependency_degrees(&alpha, &extraction_degrees, &domain)?;
    ckyz_z_residual_coefficient_work_profile(
        &path_history_degrees,
        &alpha,
        cover_weight_coefficients,
        &domain,
    )
}

fn ckyz_z_residual_coefficient_work_profile(
    extraction_degrees: &[Vec<usize>],
    alpha: &[BTreeMap<Vec<usize>, Rational>],
    cover_weight_coefficients: &[i64],
    domain: &CkyzMonomialDomain,
) -> Result<CkyzZResidualCoefficientWorkProfile> {
    validate_ckyz_target_degrees(extraction_degrees, domain.rank)?;
    if alpha.len() != domain.rank {
        return Err(Error::InvalidInput(
            "CKYZ residual coefficient profile alpha rank mismatch".into(),
        ));
    }
    if cover_weight_coefficients.len() != domain.rank {
        return Err(Error::InvalidInput(
            "CKYZ residual coefficient profile cover-weight rank mismatch".into(),
        ));
    }
    let grading_vector = ckyz_grading_vector_from_cover_weights(cover_weight_coefficients)?;
    let mut extraction_degrees = extraction_degrees.to_vec();
    ckyz_sort_degrees_for_extraction_with_grading(&mut extraction_degrees, &grading_vector)?;
    extraction_degrees.dedup();
    let extraction_gradings = extraction_degrees
        .iter()
        .map(|degree| ckyz_grading_degree(degree, &grading_vector))
        .collect::<Result<Vec<_>>>()?;
    let qn_history_profile =
        ckyz_qn_history_reuse_profile(&extraction_degrees, &extraction_gradings, domain)?;
    let mut residual_pair_count = 0usize;
    let mut same_grading_pair_skip_count = 0usize;
    let mut componentwise_pair_count = 0usize;
    let mut li2_delta_term_count = 0usize;
    let mut support_pair_count = 0usize;
    let mut support_li2_delta_term_count = 0usize;
    let mut unique_scales = BTreeSet::<Vec<usize>>::new();
    let mut unique_deltas = BTreeSet::<Vec<usize>>::new();
    let mut unique_exp_states = BTreeSet::<(Vec<usize>, Vec<usize>)>::new();
    let mut support_unique_exp_states = BTreeSet::<(Vec<usize>, Vec<usize>)>::new();
    let mut support_exp_states_by_scale = BTreeMap::<Vec<usize>, BTreeSet<Vec<usize>>>::new();
    let alpha_supports = alpha
        .iter()
        .map(|series| ckyz_series_support_indices(series, domain))
        .collect::<Vec<_>>();
    let mut exp_support_cache = HashMap::<Vec<usize>, BTreeSet<usize>>::new();

    for (degree_index, degree) in extraction_degrees.iter().enumerate() {
        let degree_grading = extraction_gradings[degree_index];
        let first_later_grading_index =
            extraction_gradings.partition_point(|&grading| grading <= degree_grading);
        same_grading_pair_skip_count = same_grading_pair_skip_count
            .checked_add(first_later_grading_index.saturating_sub(degree_index + 1))
            .ok_or_else(|| Error::InvalidInput("CKYZ same-grading pair count overflowed".into()))?;
        for target in extraction_degrees.iter().skip(first_later_grading_index) {
            residual_pair_count = residual_pair_count
                .checked_add(1)
                .ok_or_else(|| Error::InvalidInput("CKYZ residual pair count overflowed".into()))?;
            if !ckyz_componentwise_le(degree, target) {
                continue;
            }
            componentwise_pair_count =
                componentwise_pair_count.checked_add(1).ok_or_else(|| {
                    Error::InvalidInput("CKYZ componentwise pair count overflowed".into())
                })?;
            let max_multiple = degree
                .iter()
                .zip(target.iter())
                .filter_map(|(&degree_entry, &target_entry)| {
                    (degree_entry != 0).then(|| target_entry / degree_entry)
                })
                .min()
                .ok_or_else(|| {
                    Error::InvalidInput("CKYZ residual profile degree must be nonzero".into())
                })?;
            let coordinate_key = ckyz_q_degree_nonzero_coordinate_key(degree);
            if !exp_support_cache.contains_key(&coordinate_key) {
                let exp_support = ckyz_q_degree_exp_support_for_coordinate_key(
                    &coordinate_key,
                    degree,
                    &alpha_supports,
                    domain,
                )?;
                exp_support_cache.insert(coordinate_key.clone(), exp_support);
            }
            let exp_support = exp_support_cache
                .get(&coordinate_key)
                .expect("exponential support was inserted above");
            let mut pair_has_support = false;
            for multiple in 1..=max_multiple {
                let Some(delta) = ckyz_subtract_degree_multiple(target, degree, multiple) else {
                    continue;
                };
                let Some(delta_index) = domain.index_of(&delta) else {
                    continue;
                };
                let scale = degree
                    .iter()
                    .map(|entry| {
                        entry.checked_mul(multiple).ok_or_else(|| {
                            Error::InvalidInput(
                                "CKYZ residual profile scale degree overflowed".into(),
                            )
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                li2_delta_term_count = li2_delta_term_count.checked_add(1).ok_or_else(|| {
                    Error::InvalidInput("CKYZ Li2 delta term count overflowed".into())
                })?;
                unique_scales.insert(scale.clone());
                unique_deltas.insert(delta.clone());
                if exp_support.contains(&delta_index) {
                    support_li2_delta_term_count =
                        support_li2_delta_term_count.checked_add(1).ok_or_else(|| {
                            Error::InvalidInput(
                                "CKYZ support Li2 delta term count overflowed".into(),
                            )
                        })?;
                    support_unique_exp_states.insert((scale.clone(), delta.clone()));
                    support_exp_states_by_scale
                        .entry(scale.clone())
                        .or_default()
                        .insert(delta.clone());
                    pair_has_support = true;
                }
                unique_exp_states.insert((scale, delta));
            }
            if pair_has_support {
                support_pair_count = support_pair_count.checked_add(1).ok_or_else(|| {
                    Error::InvalidInput("CKYZ support pair count overflowed".into())
                })?;
            }
        }
    }
    let mut support_exp_state_counts_by_scale = support_exp_states_by_scale
        .into_iter()
        .map(|(scale, deltas)| (scale, deltas.len()))
        .collect::<Vec<_>>();
    support_exp_state_counts_by_scale.sort_by(|lhs, rhs| {
        rhs.1
            .cmp(&lhs.1)
            .then_with(|| {
                ckyz_total_degree(&lhs.0)
                    .expect("validated scale degree")
                    .cmp(&ckyz_total_degree(&rhs.0).expect("validated scale degree"))
            })
            .then_with(|| lhs.0.cmp(&rhs.0))
    });

    Ok(CkyzZResidualCoefficientWorkProfile {
        rank: domain.rank,
        domain_degree_count: domain.degrees.len(),
        path_history_degree_count: extraction_degrees.len(),
        residual_pair_count,
        same_grading_pair_skip_count,
        componentwise_pair_count,
        li2_delta_term_count,
        support_pair_count,
        support_li2_delta_term_count,
        unique_scale_count: unique_scales.len(),
        unique_delta_count: unique_deltas.len(),
        unique_exp_state_count: unique_exp_states.len(),
        support_unique_exp_state_count: support_unique_exp_states.len(),
        support_exp_state_counts_by_scale,
        qn_history_level_count: qn_history_profile.level_count,
        qn_history_candidate_hit_count: qn_history_profile.candidate_hit_count,
        qn_history_candidate_miss_count: qn_history_profile.candidate_miss_count,
        qn_history_unique_delta_count: qn_history_profile.unique_delta_count,
        q_delta_domain_profiles: ckyz_q_delta_domain_profiles(
            &qn_history_profile.unique_delta_degrees,
            alpha,
            domain,
        )?,
    })
}

fn ckyz_q_delta_domain_profiles(
    delta_degrees: &[Vec<usize>],
    alpha: &[BTreeMap<Vec<usize>, Rational>],
    domain: &CkyzMonomialDomain,
) -> Result<Vec<CkyzQDeltaDomainProfile>> {
    if alpha.len() != domain.rank {
        return Err(Error::InvalidInput(
            "CKYZ q-delta domain profile alpha rank mismatch".into(),
        ));
    }
    let alpha_support_indices = alpha
        .iter()
        .flat_map(|series| ckyz_series_support_indices(series, domain))
        .collect::<BTreeSet<_>>();
    let alpha_support_degrees = alpha_support_indices
        .into_iter()
        .map(|index| {
            domain.degrees.get(index).cloned().ok_or_else(|| {
                Error::InvalidInput("CKYZ q-delta alpha support index is outside the domain".into())
            })
        })
        .collect::<Result<Vec<_>>>()?;

    delta_degrees
        .iter()
        .map(|delta_degree| {
            ckyz_q_delta_domain_profile(delta_degree, &alpha_support_degrees, domain)
        })
        .collect()
}

fn ckyz_q_delta_domain_profile(
    delta_degree: &[usize],
    alpha_support_degrees: &[Vec<usize>],
    domain: &CkyzMonomialDomain,
) -> Result<CkyzQDeltaDomainProfile> {
    if delta_degree.len() != domain.rank {
        return Err(Error::InvalidInput(
            "CKYZ q-delta domain profile rank mismatch".into(),
        ));
    }
    let Some(delta_index) = domain.index_of(delta_degree) else {
        return Err(Error::InvalidInput(format!(
            "CKYZ q-delta domain profile degree {delta_degree:?} is outside the monomial domain"
        )));
    };

    let mut output_exp_indices = BTreeSet::new();
    for exp_index in 0..domain.degrees.len() {
        if domain.sum_index(delta_index, exp_index)?.is_some() {
            output_exp_indices.insert(exp_index);
        }
    }
    let shiftable_exp_degree_count = output_exp_indices.len();

    let mut closure_indices = output_exp_indices;
    let mut queue = VecDeque::from_iter(closure_indices.iter().copied());
    while let Some(exp_index) = queue.pop_front() {
        let exp_degree = &domain.degrees[exp_index];
        for alpha_degree in alpha_support_degrees {
            if !ckyz_componentwise_le(alpha_degree, exp_degree) {
                continue;
            }
            let Some(remainder) = ckyz_subtract_degree_multiple(exp_degree, alpha_degree, 1) else {
                continue;
            };
            let Some(remainder_index) = domain.index_of(&remainder) else {
                continue;
            };
            if closure_indices.insert(remainder_index) {
                queue.push_back(remainder_index);
            }
        }
    }
    let predecessor_closure_degree_count = closure_indices.len();
    let predecessor_closure_has_addition_table = predecessor_closure_degree_count
        .saturating_mul(predecessor_closure_degree_count)
        <= CKYZ_ADDITION_TABLE_MAX_ENTRIES;

    Ok(CkyzQDeltaDomainProfile {
        delta_degree: delta_degree.to_vec(),
        shiftable_exp_degree_count,
        predecessor_closure_degree_count,
        predecessor_closure_has_addition_table,
    })
}

/// Compute CKYZ local GV invariants for multiples of the matched target
/// direction using the source-derived causal semigroup context.
///
/// # Errors
/// Returns an error for invalid matched source data, invalid target multiples,
/// causal-domain construction failures, or non-integral GV extraction.
pub fn compute_ckyz_local_surface_gv_invariants_for_multiples_with_causal_domain(
    identification: &CkyzLocalSurfaceIdentification,
    multiples: usize,
) -> Result<BTreeMap<Vec<usize>, Integer>> {
    let target_degrees = ckyz_local_surface_target_degrees(identification, multiples)?;
    let cover_weights = ckyz_local_surface_cover_weight_coefficients(&identification.kind);
    let domain_spec = ckyz_local_surface_causal_domain_spec(identification)?;
    compute_ckyz_local_gv_invariants_for_degrees_with_causal_domain(
        &identification.source_relations,
        &identification.local_intersection_terms,
        cover_weights,
        &target_degrees,
        &domain_spec.generators,
        &domain_spec.grading_vector,
    )
}

/// Compute CKYZ logarithmic-period corrections from local relation rows.
///
/// CKYZ writes the local hypergeometric coefficient as
/// `1 / product_i Gamma(sum_a l_i^(a)(n_a + rho_a) + 1)`. This computes the
/// non-logarithmic correction terms in `partial_{rho_a} Pi_0 |_{rho=0}` for
/// all multi-degrees of total degree at most `max_total_degree`.
///
/// The returned vector is indexed by relation row / flat coordinate. Each map
/// is keyed by the multi-degree `n`.
///
/// # Errors
/// Returns an error if the relation matrix is empty, ragged, not Calabi-Yau
/// (`sum_i l_i^(a) = 0`), or if a degree conversion overflows.
pub fn compute_ckyz_log_period_corrections(
    relations: &[Vec<i64>],
    max_total_degree: usize,
) -> Result<Vec<BTreeMap<Vec<usize>, Rational>>> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    let mut corrections = vec![BTreeMap::new(); rank];
    for degree in ckyz_multi_degrees(rank, max_total_degree) {
        let point_pairings = ckyz_point_pairings(relations, &degree)?;
        let values = ckyz_log_period_coefficients_for_degree(relations, &point_pairings)?;
        for (coordinate_index, value) in values.into_iter().enumerate() {
            if value != 0 {
                corrections[coordinate_index].insert(degree.clone(), value);
            }
        }
    }
    Ok(corrections)
}

/// Compute CKYZ local double-log/prepotential period corrections.
///
/// This applies CKYZ's local intersection expression to the second
/// `rho`-derivatives of the local hypergeometric coefficient. The output is
/// still in B-model `z` coordinates; mirror-map substitution and
/// multiple-cover inversion are intentionally separate later steps.
///
/// # Errors
/// Returns an error for invalid CKYZ relation rows, out-of-range intersection
/// terms, or degree conversion overflow.
pub fn compute_ckyz_local_prepotential_period_corrections(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    max_total_degree: usize,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    for term in local_intersection_terms {
        if term.first >= rank || term.second >= rank {
            return Err(Error::InvalidInput(
                "CKYZ local intersection term index is outside the relation rank".into(),
            ));
        }
    }

    let mut corrections = BTreeMap::new();
    for degree in ckyz_multi_degrees(rank, max_total_degree) {
        let point_pairings = ckyz_point_pairings(relations, &degree)?;
        let second_derivatives =
            ckyz_second_log_period_coefficients_for_degree(relations, &point_pairings)?;
        let mut value = Rational::from(0);
        for term in local_intersection_terms {
            value += Rational::from(term.coefficient)
                * second_derivatives[term.first][term.second].clone();
        }
        if value != 0 {
            corrections.insert(degree, value);
        }
    }
    Ok(corrections)
}

/// Compute the CKYZ inverse mirror map `z(q)` from logarithmic period
/// corrections.
///
/// If the logarithmic periods have the form `t_i = log(z_i) + alpha_i(z)`,
/// then `q_i = z_i exp(alpha_i(z))`. This routine solves
/// `z_i(q) = q_i exp(-alpha_i(z(q)))` as a truncated multivariable formal
/// power series.
///
/// # Errors
/// Returns an error if the correction series are empty, have inconsistent
/// monomial ranks, contain nonzero constant terms, or overflow degree
/// arithmetic.
pub fn compute_ckyz_inverse_mirror_map(
    log_period_corrections: &[BTreeMap<Vec<usize>, Rational>],
    max_total_degree: usize,
) -> Result<Vec<BTreeMap<Vec<usize>, Rational>>> {
    let rank = log_period_corrections.len();
    if rank == 0 {
        return Err(Error::InvalidInput(
            "CKYZ inverse mirror map requires at least one coordinate".into(),
        ));
    }
    for correction in log_period_corrections {
        validate_ckyz_series(correction, rank, "CKYZ log-period correction")?;
        validate_ckyz_series_has_zero_constant(correction, rank, "CKYZ log-period correction")?;
    }

    let mut z_of_q = (0..rank)
        .map(|coordinate| ckyz_coordinate_series(rank, coordinate, max_total_degree))
        .collect::<Vec<_>>();
    for _ in 0..max_total_degree {
        let mut next = Vec::with_capacity(rank);
        for (coordinate, correction) in log_period_corrections.iter().enumerate() {
            let correction_at_z = ckyz_series_compose(correction, &z_of_q, max_total_degree)?;
            let negative_correction = ckyz_series_scale(&correction_at_z, Rational::from(-1));
            let exp_negative_correction =
                ckyz_series_exp(&negative_correction, rank, max_total_degree)?;
            let q_coordinate = ckyz_coordinate_series(rank, coordinate, max_total_degree);
            next.push(ckyz_series_mul(
                &q_coordinate,
                &exp_negative_correction,
                rank,
                max_total_degree,
            )?);
        }
        if next == z_of_q {
            break;
        }
        z_of_q = next;
    }
    Ok(z_of_q)
}

/// Substitute a CKYZ B-model `z`-series into flat `q` coordinates.
///
/// The `z_of_q` argument is normally produced by
/// [`compute_ckyz_inverse_mirror_map`].
///
/// # Errors
/// Returns an error if the series ranks are inconsistent or degree arithmetic
/// overflows.
pub fn substitute_ckyz_series_in_flat_coordinates(
    series_z: &BTreeMap<Vec<usize>, Rational>,
    z_of_q: &[BTreeMap<Vec<usize>, Rational>],
    max_total_degree: usize,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let rank = z_of_q.len();
    if rank == 0 {
        return Err(Error::InvalidInput(
            "CKYZ flat-coordinate substitution requires at least one coordinate".into(),
        ));
    }
    validate_ckyz_series(series_z, rank, "CKYZ B-model series")?;
    for argument in z_of_q {
        validate_ckyz_series(argument, rank, "CKYZ inverse mirror-map coordinate")?;
        validate_ckyz_series_has_zero_constant(
            argument,
            rank,
            "CKYZ inverse mirror-map coordinate",
        )?;
    }
    ckyz_series_compose(series_z, z_of_q, max_total_degree)
}

/// Compute CKYZ local prepotential-period corrections in flat coordinates.
///
/// This composes the B-model double-log/prepotential-period series with the
/// inverse mirror map. It still does not apply multiple-cover inversion, so its
/// output must not be interpreted as GV invariants.
///
/// # Errors
/// Returns an error for invalid CKYZ relation rows, invalid local intersection
/// terms, or formal-series degree/rank inconsistencies.
pub fn compute_ckyz_flat_prepotential_period_corrections(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    max_total_degree: usize,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let log_period_corrections = compute_ckyz_log_period_corrections(relations, max_total_degree)?;
    let z_of_q = compute_ckyz_inverse_mirror_map(&log_period_corrections, max_total_degree)?;
    let prepotential_z = compute_ckyz_local_prepotential_period_corrections(
        relations,
        local_intersection_terms,
        max_total_degree,
    )?;
    substitute_ckyz_series_in_flat_coordinates(&prepotential_z, &z_of_q, max_total_degree)
}

/// Compute CKYZ local instanton-potential corrections in flat coordinates.
///
/// This mirrors the cygv instanton-data layer for the local CKYZ source: it
/// forms `beta_ij - alpha_i alpha_j` from second and first `rho`-derivative
/// coefficients, contracts with the local base-intersection expression, and
/// substitutes the inverse mirror map. Diagonal intersection terms carry the
/// same `1/2` symmetry factor used by cygv when contracting symmetric
/// intersection data.
///
/// The returned coefficients are still the integrated local metric potential.
/// GV extraction requires the multiple-cover inversion implemented by
/// [`extract_ckyz_local_gv_invariants_from_potential`].
///
/// # Errors
/// Returns an error for invalid CKYZ relation rows, invalid local intersection
/// terms, or formal-series degree/rank inconsistencies.
pub fn compute_ckyz_local_instanton_potential_corrections(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    max_total_degree: usize,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    for term in local_intersection_terms {
        if term.first >= rank || term.second >= rank {
            return Err(Error::InvalidInput(
                "CKYZ local intersection term index is outside the relation rank".into(),
            ));
        }
    }

    let alpha = compute_ckyz_log_period_corrections(relations, max_total_degree)?;
    let z_of_q = compute_ckyz_inverse_mirror_map(&alpha, max_total_degree)?;
    let mut contracted = BTreeMap::new();

    for term in local_intersection_terms {
        let beta = ckyz_second_log_period_series_for_pair(
            relations,
            term.first,
            term.second,
            max_total_degree,
        )?;
        let alpha_product = ckyz_series_mul(
            &alpha[term.first],
            &alpha[term.second],
            rank,
            max_total_degree,
        )?;
        let mut f_pair = beta;
        ckyz_series_add_scaled_assign(&mut f_pair, &alpha_product, Rational::from(-1));

        let mut term_coefficient = Rational::from(term.coefficient);
        if term.first == term.second {
            term_coefficient /= Rational::from(2);
        }
        ckyz_series_add_scaled_assign(&mut contracted, &f_pair, term_coefficient);
    }

    substitute_ckyz_series_in_flat_coordinates(&contracted, &z_of_q, max_total_degree)
}

/// Extract local genus-zero GV invariants from CKYZ flat potential coefficients.
///
/// The input potential coefficients are interpreted as
/// `A_m = sum_{k d = m} w(d) N_d / k^2`, where
/// `w(d) = -sum_i cover_weight_coefficients_i d_i`. This is the
/// twice-integrated form of CKYZ's `instbase` expansion after the local
/// `beta - alpha alpha` conversion. The returned map stores nonzero invariants
/// only.
///
/// # Errors
/// Returns an error if ranks are inconsistent, weight arithmetic overflows, or
/// a recovered invariant is non-integral.
pub fn extract_ckyz_local_gv_invariants_from_potential(
    potential_coefficients: &BTreeMap<Vec<usize>, Rational>,
    cover_weight_coefficients: &[i64],
    max_total_degree: usize,
) -> Result<BTreeMap<Vec<usize>, Integer>> {
    let rank = cover_weight_coefficients.len();
    if rank == 0 {
        return Err(Error::InvalidInput(
            "CKYZ local GV extraction requires at least one coordinate".into(),
        ));
    }
    validate_ckyz_series(
        potential_coefficients,
        rank,
        "CKYZ local instanton potential",
    )?;
    let grading_vector = ckyz_grading_vector_from_cover_weights(cover_weight_coefficients)?;

    let mut degrees = ckyz_multi_degrees(rank, max_total_degree);
    ckyz_sort_degrees_for_extraction_with_grading(&mut degrees, &grading_vector)?;

    let mut invariants = BTreeMap::new();
    for degree in degrees {
        let mut residual = potential_coefficients
            .get(&degree)
            .cloned()
            .unwrap_or_else(|| Rational::from(0));
        for cover in ckyz_nontrivial_covers(&degree) {
            let primitive = ckyz_divide_degree(&degree, cover)?;
            let Some(gv) = invariants.get(&primitive) else {
                continue;
            };
            let weight = ckyz_cover_weight(cover_weight_coefficients, &primitive)?;
            let cover_squared = cover.checked_mul(cover).ok_or_else(|| {
                Error::InvalidInput("CKYZ cover multiplicity overflowed usize".into())
            })?;
            residual -= Rational::from(weight) * Rational::from(gv)
                / Rational::from(Integer::from(cover_squared));
        }

        if residual == 0 {
            continue;
        }
        let weight = ckyz_cover_weight(cover_weight_coefficients, &degree)?;
        if weight == 0 {
            return Err(Error::InvalidInput(
                "CKYZ local GV extraction encountered a nonzero coefficient with zero cover weight"
                    .into(),
            ));
        }
        let gv_rational = residual / Rational::from(weight);
        if gv_rational.denominator_ref() != &1u32 {
            return Err(Error::InvalidInput(format!(
                "CKYZ local GV invariant at degree {degree:?} is not integral: {gv_rational}"
            )));
        }
        let gv = Integer::try_from(gv_rational).map_err(|_| {
            Error::InvalidInput(format!(
                "CKYZ local GV invariant at degree {degree:?} is not integral"
            ))
        })?;
        if gv != 0 {
            invariants.insert(degree, gv);
        }
    }

    Ok(invariants)
}

/// Compute local CKYZ genus-zero GV invariants through a total-degree cutoff.
///
/// This combines the source-derived period coefficients, mirror-map
/// substitution, the HKTY-shaped `beta - alpha alpha` conversion, and the
/// `instbase` multiple-cover inversion. The `cover_weight_coefficients` are the
/// coefficients `x_i` in CKYZ's finite-limit expression
/// `sum_i x_i (J_i J_j) = c1(B) J_j`.
///
/// # Errors
/// Returns an error for invalid source data or non-integral extracted
/// invariants.
pub fn compute_ckyz_local_gv_invariants(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    cover_weight_coefficients: &[i64],
    max_total_degree: usize,
) -> Result<BTreeMap<Vec<usize>, Integer>> {
    let potential = compute_ckyz_local_instanton_potential_corrections(
        relations,
        local_intersection_terms,
        max_total_degree,
    )?;
    extract_ckyz_local_gv_invariants_from_potential(
        &potential,
        cover_weight_coefficients,
        max_total_degree,
    )
}

/// Compute local CKYZ genus-zero GV invariants only for requested degrees.
///
/// This follows the same source-derived path as
/// [`compute_ckyz_local_gv_invariants`], but truncates all formal series to the
/// cover-closed union of componentwise past downsets for `target_degrees`. It is
/// intended for ray checks such as McAllister potent-ray validation, where a
/// total-degree cutoff would compute many irrelevant monomials.
///
/// Cover divisors of the requested degrees are included automatically because
/// multiple-cover subtraction for a target degree depends on its primitive
/// divisors.
///
/// # Errors
/// Returns an error for invalid source data, invalid target degrees, or
/// non-integral extracted invariants.
pub fn compute_ckyz_local_gv_invariants_for_degrees(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    cover_weight_coefficients: &[i64],
    target_degrees: &[Vec<usize>],
) -> Result<BTreeMap<Vec<usize>, Integer>> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    validate_ckyz_target_degrees(target_degrees, rank)?;
    if cover_weight_coefficients.len() != rank {
        return Err(Error::InvalidInput(
            "CKYZ cover-weight rank does not match relation rank".into(),
        ));
    }

    let extraction_degrees = ckyz_cover_closed_target_degrees(target_degrees)?;
    let domain = CkyzMonomialDomain::target_downset(&extraction_degrees, rank)?;
    let potential = compute_ckyz_local_instanton_potential_corrections_domain(
        relations,
        local_intersection_terms,
        &domain,
    )?;
    let mut invariants = extract_ckyz_local_gv_invariants_from_potential_for_degrees(
        &potential,
        cover_weight_coefficients,
        &extraction_degrees,
    )?;
    invariants.retain(|degree, _| target_degrees.iter().any(|target| target == degree));
    Ok(invariants)
}

/// Compute CKYZ local GV invariants using an explicit causal semigroup domain.
///
/// This follows the same local CKYZ extraction as
/// [`compute_ckyz_local_gv_invariants_for_degrees`], but it builds the finite
/// monomial domain by closing `causal_generators` under addition up to the
/// largest requested grading degree. This mirrors the cygv polynomial-domain
/// contract more closely than a componentwise formal box: products whose summed
/// degree is absent from the generated semigroup are intentionally outside the
/// truncation domain.
///
/// The returned map contains only the requested target degrees. Cover divisors
/// are still included internally for multiple-cover subtraction.
///
/// # Errors
/// Returns an error for invalid CKYZ source data, invalid causal generators,
/// non-generated target degrees, or non-integral extracted invariants.
pub fn compute_ckyz_local_gv_invariants_for_degrees_with_causal_domain(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    cover_weight_coefficients: &[i64],
    target_degrees: &[Vec<usize>],
    causal_generators: &[Vec<usize>],
    grading_vector: &[usize],
) -> Result<BTreeMap<Vec<usize>, Integer>> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    validate_ckyz_target_degrees(target_degrees, rank)?;
    if cover_weight_coefficients.len() != rank {
        return Err(Error::InvalidInput(
            "CKYZ cover-weight rank does not match relation rank".into(),
        ));
    }

    let extraction_degrees = ckyz_cover_closed_target_degrees(target_degrees)?;
    let domain =
        ckyz_causal_monomial_domain(rank, causal_generators, grading_vector, &extraction_degrees)?;
    let potential = compute_ckyz_local_instanton_potential_corrections_domain(
        relations,
        local_intersection_terms,
        &domain,
    )?;
    let mut invariants = extract_ckyz_local_gv_invariants_from_potential_for_degrees(
        &potential,
        cover_weight_coefficients,
        &extraction_degrees,
    )?;
    invariants.retain(|degree, _| target_degrees.iter().any(|target| target == degree));
    Ok(invariants)
}

/// Compute CKYZ local GV invariants using a support-predicted finite domain.
///
/// This follows the same local CKYZ extraction as
/// [`compute_ckyz_local_gv_invariants_for_degrees`], but first predicts the
/// finite monomial domain needed by the inverse mirror map and
/// flat-coordinate potential composition at the support level. The support
/// prediction is exact with respect to the source-derived broad target downset:
/// it still enumerates alpha/beta supports on that downset, but it avoids
/// evaluating the full rational mirror-map/potential series before constructing
/// the smaller computation domain.
///
/// The returned map contains only the requested target degrees. Cover divisors
/// are included internally for multiple-cover subtraction.
///
/// # Errors
/// Returns an error for invalid CKYZ source data, invalid target degrees, or
/// non-integral extracted invariants.
pub fn compute_ckyz_local_gv_invariants_for_degrees_with_predicted_support_domain(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    cover_weight_coefficients: &[i64],
    target_degrees: &[Vec<usize>],
) -> Result<BTreeMap<Vec<usize>, Integer>> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    validate_ckyz_target_degrees(target_degrees, rank)?;
    if cover_weight_coefficients.len() != rank {
        return Err(Error::InvalidInput(
            "CKYZ cover-weight rank does not match relation rank".into(),
        ));
    }

    let extraction_degrees = ckyz_cover_closed_target_degrees(target_degrees)?;
    let domain = ckyz_predicted_support_domain_for_degrees(
        relations,
        local_intersection_terms,
        &extraction_degrees,
    )?;
    let (alpha, potential_z) = compute_ckyz_local_instanton_potential_z_domain(
        relations,
        local_intersection_terms,
        &domain,
    )?;
    let path_history_degrees =
        ckyz_z_residual_dependency_degrees(&alpha, &extraction_degrees, &domain)?;
    let mut invariants = extract_ckyz_local_gv_invariants_from_z_potential_for_degrees(
        &potential_z,
        &alpha,
        cover_weight_coefficients,
        &path_history_degrees,
        &domain,
    )?;
    invariants.retain(|degree, _| target_degrees.iter().any(|target| target == degree));
    Ok(invariants)
}

fn validate_ckyz_relations(relations: &[Vec<i64>]) -> Result<()> {
    let Some(first_relation) = relations.first() else {
        return Err(Error::InvalidInput(
            "CKYZ log-period corrections require at least one relation row".into(),
        ));
    };
    if first_relation.is_empty() {
        return Err(Error::InvalidInput(
            "CKYZ relation rows must contain at least one point".into(),
        ));
    }
    if relations
        .iter()
        .any(|relation| relation.len() != first_relation.len())
    {
        return Err(Error::InvalidInput(
            "CKYZ relation rows must all have the same length".into(),
        ));
    }
    for relation in relations {
        let row_sum = relation.iter().try_fold(0i64, |acc, &entry| {
            acc.checked_add(entry)
                .ok_or_else(|| Error::InvalidInput("CKYZ relation row sum overflowed i64".into()))
        })?;
        if row_sum != 0 {
            return Err(Error::InvalidInput(
                "CKYZ local relation rows must satisfy the Calabi-Yau sum-zero condition".into(),
            ));
        }
    }
    Ok(())
}

fn ckyz_local_surface_sources() -> Vec<CkyzLocalSurfaceSource> {
    vec![
        CkyzLocalSurfaceSource {
            kind: CkyzLocalSurfaceKind::LocalP2,
            relations: vec![vec![-3, 1, 1, 1]],
            c1_coefficients: vec![3],
            local_intersection_terms: vec![CkyzLocalIntersectionTerm {
                first: 0,
                second: 0,
                coefficient: 1,
            }],
        },
        CkyzLocalSurfaceSource {
            kind: CkyzLocalSurfaceKind::HirzebruchF0,
            relations: vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]],
            c1_coefficients: vec![2, 2],
            local_intersection_terms: vec![CkyzLocalIntersectionTerm {
                first: 0,
                second: 1,
                coefficient: 1,
            }],
        },
        CkyzLocalSurfaceSource {
            kind: CkyzLocalSurfaceKind::HirzebruchF1,
            relations: vec![vec![-2, 1, 0, 1, 0], vec![-1, 0, 1, -1, 1]],
            c1_coefficients: vec![3, 2],
            local_intersection_terms: vec![
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 1,
                    coefficient: 1,
                },
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 0,
                    coefficient: 1,
                },
            ],
        },
        CkyzLocalSurfaceSource {
            kind: CkyzLocalSurfaceKind::Polygon5,
            relations: vec![
                vec![-1, 1, -1, 1, 0, 0],
                vec![-1, -1, 1, 0, 0, 1],
                vec![-1, 0, 1, -1, 1, 0],
            ],
            c1_coefficients: vec![3, 2, 2],
            local_intersection_terms: vec![
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 0,
                    coefficient: 1,
                },
                CkyzLocalIntersectionTerm {
                    first: 1,
                    second: 0,
                    coefficient: 1,
                },
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 2,
                    coefficient: 1,
                },
                CkyzLocalIntersectionTerm {
                    first: 1,
                    second: 2,
                    coefficient: 1,
                },
            ],
        },
    ]
}

fn ckyz_multi_degrees(rank: usize, max_total_degree: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut current = vec![0; rank];
    push_ckyz_multi_degrees(rank, max_total_degree, 0, &mut current, &mut out);
    out.retain(|degree| degree.iter().any(|&entry| entry != 0));
    out
}

fn ckyz_multi_degrees_box(max_degrees: &[usize]) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut current = vec![0; max_degrees.len()];
    push_ckyz_multi_degrees_box(max_degrees, 0, &mut current, &mut out);
    out.retain(|degree| degree.iter().any(|&entry| entry != 0));
    out
}

fn push_ckyz_multi_degrees_box(
    max_degrees: &[usize],
    coordinate_index: usize,
    current: &mut [usize],
    out: &mut Vec<Vec<usize>>,
) {
    if coordinate_index == max_degrees.len() {
        out.push(current.to_vec());
        return;
    }
    for value in 0..=max_degrees[coordinate_index] {
        current[coordinate_index] = value;
        push_ckyz_multi_degrees_box(max_degrees, coordinate_index + 1, current, out);
    }
}

fn push_ckyz_multi_degrees(
    rank: usize,
    remaining_degree: usize,
    coordinate_index: usize,
    current: &mut [usize],
    out: &mut Vec<Vec<usize>>,
) {
    if coordinate_index == rank {
        out.push(current.to_vec());
        return;
    }
    for value in 0..=remaining_degree {
        current[coordinate_index] = value;
        push_ckyz_multi_degrees(
            rank,
            remaining_degree - value,
            coordinate_index + 1,
            current,
            out,
        );
    }
}

fn ckyz_point_pairings(relations: &[Vec<i64>], degree: &[usize]) -> Result<Vec<i64>> {
    let n_points = relations[0].len();
    let mut pairings = vec![0i64; n_points];
    for (relation, &degree_entry) in relations.iter().zip(degree.iter()) {
        if degree_entry == 0 {
            continue;
        }
        let degree_entry = i64::try_from(degree_entry)
            .map_err(|_| Error::InvalidInput("CKYZ multi-degree entry does not fit i64".into()))?;
        for (pairing, &relation_entry) in pairings.iter_mut().zip(relation.iter()) {
            let contribution = relation_entry.checked_mul(degree_entry).ok_or_else(|| {
                Error::InvalidInput("CKYZ point pairing multiplication overflowed i64".into())
            })?;
            *pairing = pairing.checked_add(contribution).ok_or_else(|| {
                Error::InvalidInput("CKYZ point pairing addition overflowed i64".into())
            })?;
        }
    }
    Ok(pairings)
}

fn ckyz_log_period_coefficients_for_degree(
    relations: &[Vec<i64>],
    point_pairings: &[i64],
) -> Result<Vec<Rational>> {
    let rank = relations.len();
    let negative_points = point_pairings
        .iter()
        .enumerate()
        .filter_map(|(index, &pairing)| (pairing < 0).then_some(index))
        .collect::<Vec<_>>();
    let mut values = vec![Rational::from(0); rank];

    match negative_points.as_slice() {
        [negative_point] => {
            let coefficient = ckyz_one_negative_log_period_base(point_pairings, *negative_point)?;
            for coordinate_index in 0..rank {
                let relation_entry = relations[coordinate_index][*negative_point];
                if relation_entry != 0 {
                    values[coordinate_index] = coefficient.clone() * Rational::from(relation_entry);
                }
            }
        }
        [] => {
            let coefficient = ckyz_zero_negative_hypergeometric_coeff(point_pairings)?;
            let regular_terms = ckyz_regular_harmonic_terms(relations, point_pairings)?;
            for coordinate_index in 0..rank {
                values[coordinate_index] =
                    coefficient.clone() * regular_terms[coordinate_index].clone();
            }
        }
        _ => {}
    }
    Ok(values)
}

fn ckyz_second_log_period_coefficients_for_degree(
    relations: &[Vec<i64>],
    point_pairings: &[i64],
) -> Result<Vec<Vec<Rational>>> {
    let rank = relations.len();
    let negative_points = point_pairings
        .iter()
        .enumerate()
        .filter_map(|(index, &pairing)| (pairing < 0).then_some(index))
        .collect::<Vec<_>>();
    let mut values = vec![vec![Rational::from(0); rank]; rank];

    match negative_points.as_slice() {
        [] => {
            let coefficient = ckyz_zero_negative_hypergeometric_coeff(point_pairings)?;
            let regular_terms = ckyz_regular_harmonic_terms(relations, point_pairings)?;
            for first in 0..rank {
                for second in 0..rank {
                    let mut term = regular_terms[first].clone() * regular_terms[second].clone();
                    for point_index in 0..point_pairings.len() {
                        let arg = ckyz_harmonic_argument(point_pairings[point_index])?;
                        term += Rational::from(relations[first][point_index])
                            * Rational::from(relations[second][point_index])
                            * harmonic_number_order_two(arg);
                    }
                    values[first][second] = coefficient.clone() * term;
                }
            }
        }
        [negative_point] => {
            let coefficient = ckyz_one_negative_log_period_base(point_pairings, *negative_point)?;
            let regular_terms = ckyz_regular_harmonic_terms(relations, point_pairings)?;
            for first in 0..rank {
                for second in 0..rank {
                    values[first][second] = coefficient.clone()
                        * (regular_terms[first].clone()
                            * Rational::from(relations[second][*negative_point])
                            + regular_terms[second].clone()
                                * Rational::from(relations[first][*negative_point]));
                }
            }
        }
        [first_negative, second_negative] => {
            let coefficient = ckyz_two_negative_double_log_base(
                point_pairings,
                *first_negative,
                *second_negative,
            )?;
            for first in 0..rank {
                for second in 0..rank {
                    values[first][second] = coefficient.clone()
                        * (Rational::from(relations[first][*first_negative])
                            * Rational::from(relations[second][*second_negative])
                            + Rational::from(relations[first][*second_negative])
                                * Rational::from(relations[second][*first_negative]));
                }
            }
        }
        _ => {}
    }
    Ok(values)
}

fn ckyz_second_log_period_series_for_pair(
    relations: &[Vec<i64>],
    first: usize,
    second: usize,
    max_total_degree: usize,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let rank = relations.len();
    if first >= rank || second >= rank {
        return Err(Error::InvalidInput(
            "CKYZ second-log period pair index is outside the relation rank".into(),
        ));
    }

    let mut out = BTreeMap::new();
    for degree in ckyz_multi_degrees(rank, max_total_degree) {
        let point_pairings = ckyz_point_pairings(relations, &degree)?;
        let values = ckyz_second_log_period_coefficients_for_degree(relations, &point_pairings)?;
        let value = values[first][second].clone();
        if value != 0 {
            out.insert(degree, value);
        }
    }
    Ok(out)
}

fn ckyz_second_log_period_series_for_pair_domain(
    relations: &[Vec<i64>],
    first: usize,
    second: usize,
    domain: &CkyzMonomialDomain,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let rank = relations.len();
    if first >= rank || second >= rank {
        return Err(Error::InvalidInput(
            "CKYZ second-log period pair index is outside the relation rank".into(),
        ));
    }
    if domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ componentwise cutoff rank does not match relation rank".into(),
        ));
    }

    let mut out = BTreeMap::new();
    for degree in domain.nonzero_degrees() {
        let point_pairings = ckyz_point_pairings(relations, degree)?;
        let values = ckyz_second_log_period_coefficients_for_degree(relations, &point_pairings)?;
        let value = values[first][second].clone();
        if value != 0 {
            out.insert(degree.clone(), value);
        }
    }
    Ok(out)
}

fn ckyz_regular_harmonic_terms(
    relations: &[Vec<i64>],
    point_pairings: &[i64],
) -> Result<Vec<Rational>> {
    let rank = relations.len();
    let mut terms = vec![Rational::from(0); rank];
    for coordinate_index in 0..rank {
        for (&relation_entry, &pairing) in relations[coordinate_index]
            .iter()
            .zip(point_pairings.iter())
        {
            if relation_entry == 0 {
                continue;
            }
            let harmonic_arg = ckyz_harmonic_argument(pairing)?;
            if harmonic_arg == 0 {
                continue;
            }
            terms[coordinate_index] -=
                Rational::from(relation_entry) * harmonic_number(harmonic_arg);
        }
    }
    Ok(terms)
}

fn ckyz_harmonic_argument(pairing: i64) -> Result<usize> {
    let arg = if pairing < 0 { -pairing - 1 } else { pairing };
    usize::try_from(arg)
        .map_err(|_| Error::InvalidInput("CKYZ harmonic argument does not fit usize".into()))
}

fn ckyz_one_negative_log_period_base(
    point_pairings: &[i64],
    negative_point: usize,
) -> Result<Rational> {
    let negative_pairing = point_pairings[negative_point];
    if negative_pairing >= 0 {
        return Err(Error::InvalidInput(
            "CKYZ one-negative coefficient received a nonnegative point".into(),
        ));
    }
    let pole_order = usize::try_from(-negative_pairing - 1).map_err(|_| {
        Error::InvalidInput("CKYZ negative point pairing does not fit usize".into())
    })?;
    let mut coefficient = Rational::from(factorial_integer(pole_order));
    for (point_index, &pairing) in point_pairings.iter().enumerate() {
        if point_index == negative_point {
            continue;
        }
        if pairing < 0 {
            return Err(Error::InvalidInput(
                "CKYZ one-negative coefficient received multiple negative points".into(),
            ));
        }
        let pairing = usize::try_from(pairing).map_err(|_| {
            Error::InvalidInput("CKYZ nonnegative point pairing does not fit usize".into())
        })?;
        coefficient /= Rational::from(factorial_integer(pairing));
    }
    if negative_pairing % 2 == 0 {
        coefficient = -coefficient;
    }
    Ok(coefficient)
}

fn ckyz_two_negative_double_log_base(
    point_pairings: &[i64],
    first_negative: usize,
    second_negative: usize,
) -> Result<Rational> {
    let mut coefficient = Rational::from(1);
    for (point_index, &pairing) in point_pairings.iter().enumerate() {
        if point_index == first_negative || point_index == second_negative {
            if pairing >= 0 {
                return Err(Error::InvalidInput(
                    "CKYZ two-negative coefficient received a nonnegative point".into(),
                ));
            }
            let pole_order = usize::try_from(-pairing - 1).map_err(|_| {
                Error::InvalidInput("CKYZ negative point pairing does not fit usize".into())
            })?;
            coefficient *= Rational::from(factorial_integer(pole_order));
        } else {
            if pairing < 0 {
                return Err(Error::InvalidInput(
                    "CKYZ two-negative coefficient received more than two negative points".into(),
                ));
            }
            let pairing = usize::try_from(pairing).map_err(|_| {
                Error::InvalidInput("CKYZ nonnegative point pairing does not fit usize".into())
            })?;
            coefficient /= Rational::from(factorial_integer(pairing));
        }
    }
    if (point_pairings[first_negative] + point_pairings[second_negative]) % 2 != 0 {
        coefficient = -coefficient;
    }
    Ok(coefficient)
}

fn ckyz_zero_negative_hypergeometric_coeff(point_pairings: &[i64]) -> Result<Rational> {
    let mut coefficient = Rational::from(1);
    for &pairing in point_pairings {
        if pairing < 0 {
            return Err(Error::InvalidInput(
                "CKYZ zero-negative coefficient received a negative point".into(),
            ));
        }
        let pairing = usize::try_from(pairing).map_err(|_| {
            Error::InvalidInput("CKYZ nonnegative point pairing does not fit usize".into())
        })?;
        coefficient /= Rational::from(factorial_integer(pairing));
    }
    Ok(coefficient)
}

fn harmonic_number(n: usize) -> Rational {
    let mut out = Rational::from(0);
    for denominator in 1..=n {
        out += Rational::from_signeds(1i64, i64::try_from(denominator).expect("usize fits i64"));
    }
    out
}

fn harmonic_number_order_two(n: usize) -> Rational {
    let mut out = Rational::from(0);
    for denominator in 1..=n {
        let denominator = i64::try_from(denominator).expect("usize fits i64");
        out += Rational::from_signeds(1i64, denominator * denominator);
    }
    out
}

fn validate_ckyz_series(
    series: &BTreeMap<Vec<usize>, Rational>,
    rank: usize,
    context: &str,
) -> Result<()> {
    for degree in series.keys() {
        if degree.len() != rank {
            return Err(Error::InvalidInput(format!(
                "{context} monomial rank does not match coordinate rank"
            )));
        }
        ckyz_total_degree(degree)?;
    }
    Ok(())
}

fn validate_ckyz_series_has_zero_constant(
    series: &BTreeMap<Vec<usize>, Rational>,
    rank: usize,
    context: &str,
) -> Result<()> {
    let zero_degree = vec![0; rank];
    if series
        .get(&zero_degree)
        .is_some_and(|coefficient| *coefficient != 0)
    {
        return Err(Error::InvalidInput(format!(
            "{context} must have zero constant term"
        )));
    }
    Ok(())
}

fn ckyz_total_degree(degree: &[usize]) -> Result<usize> {
    degree.iter().try_fold(0usize, |sum, &entry| {
        sum.checked_add(entry)
            .ok_or_else(|| Error::InvalidInput("CKYZ total degree overflowed usize".into()))
    })
}

#[derive(Clone, Debug)]
struct CkyzMonomialDomain {
    rank: usize,
    degrees: Vec<Vec<usize>>,
    degree_indices: HashMap<Vec<usize>, usize>,
    max_coordinate_degrees: Vec<usize>,
    dense_degree_strides: Option<Vec<usize>>,
    dense_degree_indices: Option<Vec<usize>>,
    addition_indices: Option<Vec<usize>>,
    addition_pairs_by_lhs: Option<Vec<Vec<(usize, usize)>>>,
    max_total_degree: usize,
}

impl CkyzMonomialDomain {
    #[cfg(test)]
    fn componentwise_box(max_degrees: &[usize]) -> Result<Self> {
        if max_degrees.is_empty() {
            return Err(Error::InvalidInput(
                "CKYZ monomial domain requires at least one coordinate".into(),
            ));
        }
        let rank = max_degrees.len();
        Self::from_degrees(rank, ckyz_multi_degrees_box(max_degrees))
    }

    fn target_downset(target_degrees: &[Vec<usize>], rank: usize) -> Result<Self> {
        validate_ckyz_target_degrees(target_degrees, rank)?;
        let mut degrees = Vec::new();
        for target in target_degrees {
            degrees.extend(ckyz_multi_degrees_box(target));
        }
        Self::from_degrees(rank, degrees)
    }

    fn from_degrees<I>(rank: usize, degrees: I) -> Result<Self>
    where
        I: IntoIterator<Item = Vec<usize>>,
    {
        if rank == 0 {
            return Err(Error::InvalidInput(
                "CKYZ monomial domain requires at least one coordinate".into(),
            ));
        }
        let mut degree_set = BTreeSet::new();
        degree_set.insert(vec![0; rank]);
        let mut max_total_degree = 0usize;
        for degree in degrees {
            if degree.len() != rank {
                return Err(Error::InvalidInput(
                    "CKYZ monomial domain degree rank mismatch".into(),
                ));
            }
            max_total_degree = max_total_degree.max(ckyz_total_degree(&degree)?);
            degree_set.insert(degree);
        }
        let degrees = degree_set.into_iter().collect::<Vec<_>>();
        let degree_indices = degrees
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, degree)| (degree, index))
            .collect::<HashMap<_, _>>();
        let mut max_coordinate_degrees = vec![0; rank];
        for degree in &degrees {
            for (coordinate, &entry) in degree.iter().enumerate() {
                max_coordinate_degrees[coordinate] = max_coordinate_degrees[coordinate].max(entry);
            }
        }
        let mut dense_entries = Some(1usize);
        for &max_degree in &max_coordinate_degrees {
            dense_entries =
                dense_entries.and_then(|entries| entries.checked_mul(max_degree.checked_add(1)?));
        }
        let (dense_degree_strides, dense_degree_indices) = dense_entries
            .filter(|&entries| entries <= CKYZ_DENSE_DEGREE_INDEX_MAX_ENTRIES)
            .map_or((None, None), |entries| {
                let mut strides = Vec::with_capacity(rank);
                let mut stride = 1usize;
                for &max_degree in &max_coordinate_degrees {
                    strides.push(stride);
                    stride = stride
                        .checked_mul(max_degree + 1)
                        .expect("dense entries were checked above");
                }
                let mut dense_indices = vec![CKYZ_ABSENT_DEGREE_INDEX; entries];
                for (index, degree) in degrees.iter().enumerate() {
                    let dense_index =
                        ckyz_dense_degree_index(degree, rank, &max_coordinate_degrees, &strides)
                            .expect("domain degree is inside its own dense index bounds")
                            .expect("domain degree has a dense index");
                    dense_indices[dense_index] = index;
                }
                (Some(strides), Some(dense_indices))
            });
        let addition_entries = degrees.len().saturating_mul(degrees.len());
        let addition_indices = if addition_entries <= CKYZ_ADDITION_TABLE_MAX_ENTRIES {
            let mut addition_indices = vec![CKYZ_ABSENT_ADDITION_INDEX; addition_entries];
            let mut addition_pairs_by_lhs = vec![Vec::new(); degrees.len()];
            let addition_stride = degrees.len();
            for (lhs_index, lhs_degree) in degrees.iter().enumerate() {
                for (rhs_index, rhs_degree) in degrees.iter().enumerate() {
                    let sum_index =
                        ckyz_sum_degree_index(lhs_degree, rhs_degree, rank, &degree_indices)?;
                    addition_indices[lhs_index * addition_stride + rhs_index] =
                        sum_index.unwrap_or(CKYZ_ABSENT_ADDITION_INDEX);
                    if let Some(sum_index) = sum_index {
                        addition_pairs_by_lhs[lhs_index].push((rhs_index, sum_index));
                    }
                }
            }
            Some((addition_indices, addition_pairs_by_lhs))
        } else {
            None
        };
        let (addition_indices, addition_pairs_by_lhs) = addition_indices
            .map_or((None, None), |(indices, pairs)| {
                (Some(indices), Some(pairs))
            });
        Ok(Self {
            rank,
            degrees,
            degree_indices,
            max_coordinate_degrees,
            dense_degree_strides,
            dense_degree_indices,
            addition_indices,
            addition_pairs_by_lhs,
            max_total_degree,
        })
    }

    fn contains(&self, degree: &[usize]) -> bool {
        self.index_of(degree).is_some()
    }

    fn nonzero_degrees(&self) -> impl Iterator<Item = &Vec<usize>> {
        self.degrees
            .iter()
            .filter(|degree| degree.iter().any(|&entry| entry != 0))
    }

    fn index_of(&self, degree: &[usize]) -> Option<usize> {
        if let (Some(strides), Some(dense_indices)) =
            (&self.dense_degree_strides, &self.dense_degree_indices)
        {
            let Some(dense_index) =
                ckyz_dense_degree_index(degree, self.rank, &self.max_coordinate_degrees, strides)
                    .ok()
                    .flatten()
            else {
                return None;
            };
            let index = dense_indices[dense_index];
            return (index != CKYZ_ABSENT_DEGREE_INDEX).then_some(index);
        }
        self.degree_indices.get(degree).copied()
    }

    fn sum_index(&self, lhs_index: usize, rhs_index: usize) -> Result<Option<usize>> {
        if lhs_index >= self.degrees.len() || rhs_index >= self.degrees.len() {
            return Err(Error::InvalidInput(
                "CKYZ monomial addition index is outside the domain".into(),
            ));
        }
        if let Some(addition_indices) = &self.addition_indices {
            let table_index = lhs_index * self.degrees.len() + rhs_index;
            let sum_index = addition_indices[table_index];
            return Ok((sum_index != CKYZ_ABSENT_ADDITION_INDEX).then_some(sum_index));
        }
        let lhs_degree = &self.degrees[lhs_index];
        let rhs_degree = &self.degrees[rhs_index];
        if let (Some(strides), Some(dense_indices)) =
            (&self.dense_degree_strides, &self.dense_degree_indices)
        {
            let Some(dense_index) = ckyz_dense_domain_degree_sum_index(
                lhs_degree,
                rhs_degree,
                self.rank,
                &self.max_coordinate_degrees,
                strides,
            ) else {
                return Ok(None);
            };
            let sum_index = dense_indices[dense_index];
            return Ok((sum_index != CKYZ_ABSENT_DEGREE_INDEX).then_some(sum_index));
        }
        ckyz_sum_degree_index(lhs_degree, rhs_degree, self.rank, &self.degree_indices)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CkyzIndexedSeries {
    terms: Vec<(usize, Rational)>,
}

impl CkyzIndexedSeries {
    fn one(domain: &CkyzMonomialDomain) -> Result<Self> {
        let identity_index = ckyz_support_identity_index(domain)?;
        Ok(Self {
            terms: vec![(identity_index, Rational::from(1))],
        })
    }

    fn coordinate(rank: usize, coordinate: usize, domain: &CkyzMonomialDomain) -> Self {
        if coordinate >= rank || rank != domain.rank {
            return Self { terms: Vec::new() };
        }
        let mut degree = vec![0; rank];
        degree[coordinate] = 1;
        let terms = domain
            .index_of(&degree)
            .map_or_else(Vec::new, |index| vec![(index, Rational::from(1))]);
        Self { terms }
    }

    #[cfg(test)]
    fn monomial(
        degree: &[usize],
        coefficient: Rational,
        domain: &CkyzMonomialDomain,
        context: &str,
    ) -> Result<Self> {
        if degree.len() != domain.rank {
            return Err(Error::InvalidInput(format!(
                "{context} monomial rank does not match coordinate rank"
            )));
        }
        if coefficient == 0 {
            return Ok(Self { terms: Vec::new() });
        }
        let terms = domain
            .index_of(degree)
            .map_or_else(Vec::new, |index| vec![(index, coefficient)]);
        Ok(Self { terms })
    }

    fn from_btree(
        series: &BTreeMap<Vec<usize>, Rational>,
        domain: &CkyzMonomialDomain,
        context: &str,
    ) -> Result<Self> {
        let mut terms = Vec::with_capacity(series.len());
        for (degree, coefficient) in series {
            if *coefficient == 0 {
                continue;
            }
            if degree.len() != domain.rank {
                return Err(Error::InvalidInput(format!(
                    "{context} monomial rank does not match coordinate rank"
                )));
            }
            if let Some(index) = domain.index_of(degree) {
                terms.push((index, coefficient.clone()));
            }
        }
        terms.sort_unstable_by_key(|(index, _)| *index);
        Ok(Self { terms })
    }

    fn from_domain_coefficients(coefficients: Vec<Option<Rational>>) -> Self {
        let terms = coefficients
            .into_iter()
            .enumerate()
            .filter_map(|(index, coefficient)| {
                coefficient
                    .and_then(|coefficient| (coefficient != 0).then_some((index, coefficient)))
            })
            .collect();
        Self { terms }
    }

    fn to_btree(&self, domain: &CkyzMonomialDomain) -> BTreeMap<Vec<usize>, Rational> {
        self.terms
            .iter()
            .filter_map(|(index, coefficient)| {
                (coefficient != &Rational::from(0))
                    .then(|| (domain.degrees[*index].clone(), coefficient.clone()))
            })
            .collect()
    }

    fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    fn min_total_degree(
        &self,
        domain: &CkyzMonomialDomain,
        context: &str,
    ) -> Result<Option<usize>> {
        let mut min_degree = None;
        for &(index, _) in &self.terms {
            let Some(degree) = domain.degrees.get(index) else {
                return Err(Error::InvalidInput(format!(
                    "{context} monomial index is outside the domain"
                )));
            };
            let total_degree = ckyz_total_degree(degree)?;
            if total_degree == 0 {
                return Err(Error::InvalidInput(format!(
                    "{context} must have zero constant term"
                )));
            }
            min_degree =
                Some(min_degree.map_or(total_degree, |current: usize| current.min(total_degree)));
        }
        Ok(min_degree)
    }

    fn add_scaled_assign(&mut self, rhs: &Self, scalar: Rational) {
        if scalar == 0 || rhs.terms.is_empty() {
            return;
        }

        let lhs = std::mem::take(&mut self.terms);
        let mut merged = Vec::with_capacity(lhs.len().max(rhs.terms.len()));
        let mut lhs_iter = lhs.into_iter().peekable();
        let mut rhs_iter = rhs.terms.iter().peekable();

        loop {
            match (lhs_iter.peek(), rhs_iter.peek()) {
                (Some((lhs_index, _)), Some((rhs_index, _))) => {
                    if lhs_index < rhs_index {
                        merged.push(lhs_iter.next().expect("lhs term was peeked"));
                    } else if rhs_index < lhs_index {
                        let (index, coefficient) = rhs_iter.next().expect("rhs term was peeked");
                        let scaled = coefficient.clone() * scalar.clone();
                        if scaled != 0 {
                            merged.push((*index, scaled));
                        }
                    } else {
                        let (index, lhs_coefficient) =
                            lhs_iter.next().expect("lhs term was peeked");
                        let (_, rhs_coefficient) = rhs_iter.next().expect("rhs term was peeked");
                        let coefficient =
                            lhs_coefficient + rhs_coefficient.clone() * scalar.clone();
                        if coefficient != 0 {
                            merged.push((index, coefficient));
                        }
                    }
                }
                (Some(_), None) => {
                    merged.extend(lhs_iter);
                    break;
                }
                (None, Some(_)) => {
                    for (index, coefficient) in rhs_iter {
                        let scaled = coefficient.clone() * scalar.clone();
                        if scaled != 0 {
                            merged.push((*index, scaled));
                        }
                    }
                    break;
                }
                (None, None) => break,
            }
        }

        self.terms = merged;
    }

    fn scaled(&self, scalar: Rational) -> Self {
        if scalar == 0 {
            return Self { terms: Vec::new() };
        }
        let terms = self
            .terms
            .iter()
            .filter_map(|(index, coefficient)| {
                let scaled = coefficient.clone() * scalar.clone();
                (scaled != 0).then_some((*index, scaled))
            })
            .collect();
        Self { terms }
    }

    fn mul(&self, rhs: &Self, domain: &CkyzMonomialDomain) -> Result<Self> {
        let mut out_by_index = vec![None::<Rational>; domain.degrees.len()];
        if let Some(addition_pairs_by_lhs) = &domain.addition_pairs_by_lhs {
            let (outer_terms, inner_terms) = if self.terms.len() <= rhs.terms.len() {
                (&self.terms, &rhs.terms)
            } else {
                (&rhs.terms, &self.terms)
            };
            let indexed_candidate_count = outer_terms
                .iter()
                .map(|(outer_index, _)| addition_pairs_by_lhs[*outer_index].len())
                .sum::<usize>();
            let direct_candidate_count = outer_terms
                .len()
                .checked_mul(inner_terms.len())
                .ok_or_else(|| {
                    Error::InvalidInput(
                        "CKYZ sparse product candidate count overflowed usize".into(),
                    )
                })?;
            if direct_candidate_count <= indexed_candidate_count {
                for (outer_index, outer_coefficient) in outer_terms {
                    for (inner_index, inner_coefficient) in inner_terms {
                        if let Some(product_index) = domain.sum_index(*outer_index, *inner_index)? {
                            let entry = out_by_index[product_index]
                                .get_or_insert_with(|| Rational::from(0));
                            *entry += outer_coefficient.clone() * inner_coefficient.clone();
                        }
                    }
                }
                return Ok(Self::from_domain_coefficients(out_by_index));
            }

            let mut inner_by_index = vec![None::<&Rational>; domain.degrees.len()];
            for (inner_index, inner_coefficient) in inner_terms {
                inner_by_index[*inner_index] = Some(inner_coefficient);
            }
            for (outer_index, outer_coefficient) in outer_terms {
                for &(inner_index, product_index) in &addition_pairs_by_lhs[*outer_index] {
                    let Some(inner_coefficient) = inner_by_index[inner_index] else {
                        continue;
                    };
                    let entry =
                        out_by_index[product_index].get_or_insert_with(|| Rational::from(0));
                    *entry += outer_coefficient.clone() * inner_coefficient.clone();
                }
            }
            return Ok(Self::from_domain_coefficients(out_by_index));
        }

        for (lhs_index, lhs_coefficient) in &self.terms {
            for (rhs_index, rhs_coefficient) in &rhs.terms {
                if let Some(product_index) = domain.sum_index(*lhs_index, *rhs_index)? {
                    let entry =
                        out_by_index[product_index].get_or_insert_with(|| Rational::from(0));
                    *entry += lhs_coefficient.clone() * rhs_coefficient.clone();
                }
            }
        }
        Ok(Self::from_domain_coefficients(out_by_index))
    }

    fn power_cache(&self, max_exponent: usize, domain: &CkyzMonomialDomain) -> Result<Vec<Self>> {
        let mut powers = Vec::with_capacity(max_exponent + 1);
        powers.push(Self::one(domain)?);
        for exponent in 1..=max_exponent {
            let next = powers[exponent - 1].mul(self, domain)?;
            if next.is_empty() {
                powers.push(next);
                powers.resize_with(max_exponent + 1, || Self { terms: Vec::new() });
                break;
            }
            powers.push(next);
        }
        Ok(powers)
    }

    fn compose(&self, arguments: &[Self], domain: &CkyzMonomialDomain) -> Result<Self> {
        let rank = arguments.len();
        if rank != domain.rank {
            return Err(Error::InvalidInput(
                "CKYZ series composition rank mismatch".into(),
            ));
        }

        let mut max_exponents = vec![0usize; rank];
        for (degree_index, _) in &self.terms {
            let degree = &domain.degrees[*degree_index];
            for (coordinate, &exponent) in degree.iter().enumerate() {
                max_exponents[coordinate] = max_exponents[coordinate].max(exponent);
            }
        }
        let power_caches = arguments
            .iter()
            .zip(max_exponents.iter())
            .map(|(argument, &max_exponent)| argument.power_cache(max_exponent, domain))
            .collect::<Result<Vec<_>>>()?;

        let mut out = Self { terms: Vec::new() };
        for (degree_index, coefficient) in &self.terms {
            let degree = &domain.degrees[*degree_index];
            let mut monomial = Self::one(domain)?;
            for (coordinate, &exponent) in degree.iter().enumerate() {
                if exponent == 0 {
                    continue;
                }
                monomial = monomial.mul(&power_caches[coordinate][exponent], domain)?;
                if monomial.is_empty() {
                    break;
                }
            }
            out.add_scaled_assign(&monomial, coefficient.clone());
        }
        Ok(out)
    }

    fn exp(&self, domain: &CkyzMonomialDomain) -> Result<Self> {
        let max_exponent = self
            .min_total_degree(domain, "CKYZ exponential input")?
            .map_or(0, |min_degree| domain.max_total_degree / min_degree);

        let mut out = Self::one(domain)?;
        let mut power = out.clone();
        let mut factorial = Integer::from(1);
        for exponent in 1..=max_exponent {
            power = power.mul(self, domain)?;
            if power.is_empty() {
                break;
            }
            factorial *= Integer::from(exponent);
            out.add_scaled_assign(
                &power,
                Rational::from(1) / Rational::from(factorial.clone()),
            );
        }
        Ok(out)
    }

    #[cfg(test)]
    fn li2(&self, domain: &CkyzMonomialDomain) -> Result<Self> {
        let max_exponent = self
            .min_total_degree(domain, "CKYZ Li2 input")?
            .map_or(0, |min_degree| domain.max_total_degree / min_degree);
        if max_exponent == 0 {
            return Ok(Self { terms: Vec::new() });
        }

        let mut out = self.clone();
        let mut power = self.clone();
        for exponent in 2..=max_exponent {
            power = power.mul(self, domain)?;
            if power.is_empty() {
                break;
            }
            let exponent_squared = exponent
                .checked_mul(exponent)
                .ok_or_else(|| Error::InvalidInput("CKYZ Li2 exponent overflowed usize".into()))?;
            out.add_scaled_assign(
                &power,
                Rational::from(1) / Rational::from(Integer::from(exponent_squared)),
            );
        }
        Ok(out)
    }
}

fn ckyz_causal_monomial_domain(
    rank: usize,
    generators: &[Vec<usize>],
    grading_vector: &[usize],
    target_degrees: &[Vec<usize>],
) -> Result<CkyzMonomialDomain> {
    if rank == 0 {
        return Err(Error::InvalidInput(
            "CKYZ causal domain requires at least one coordinate".into(),
        ));
    }
    if generators.is_empty() {
        return Err(Error::InvalidInput(
            "CKYZ causal domain requires at least one generator".into(),
        ));
    }
    if grading_vector.len() != rank {
        return Err(Error::InvalidInput(
            "CKYZ causal domain grading rank mismatch".into(),
        ));
    }
    validate_ckyz_target_degrees(target_degrees, rank)?;

    let mut normalized_generators = Vec::with_capacity(generators.len());
    for generator in generators {
        if generator.len() != rank {
            return Err(Error::InvalidInput(
                "CKYZ causal generator rank mismatch".into(),
            ));
        }
        if generator.iter().all(|&entry| entry == 0) {
            return Err(Error::InvalidInput(
                "CKYZ causal generators must be nonzero".into(),
            ));
        }
        let degree = ckyz_grading_degree(generator, grading_vector)?;
        if degree == 0 {
            return Err(Error::InvalidInput(
                "CKYZ causal generator has zero grading degree".into(),
            ));
        }
        normalized_generators.push(generator.clone());
    }

    let max_grading = target_degrees
        .iter()
        .map(|degree| ckyz_grading_degree(degree, grading_vector))
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .max()
        .expect("validated target degrees are nonempty");

    let zero = vec![0; rank];
    let mut degree_set = BTreeSet::from([zero.clone()]);
    let mut frontier = VecDeque::from([zero]);
    while let Some(current) = frontier.pop_front() {
        for generator in &normalized_generators {
            let next = ckyz_add_degrees_unbounded(&current, generator, rank)?;
            if ckyz_grading_degree(&next, grading_vector)? > max_grading {
                continue;
            }
            if degree_set.insert(next.clone()) {
                frontier.push_back(next);
            }
        }
    }

    for target in target_degrees {
        if !degree_set.contains(target) {
            return Err(Error::InvalidInput(format!(
                "CKYZ causal generators do not generate target degree {target:?}"
            )));
        }
    }

    CkyzMonomialDomain::from_degrees(rank, degree_set)
}

#[cfg(test)]
fn ckyz_record_series_support(
    degree_set: &mut BTreeSet<Vec<usize>>,
    series: &BTreeMap<Vec<usize>, Rational>,
) {
    degree_set.extend(
        series
            .iter()
            .filter(|(_, coefficient)| **coefficient != 0)
            .map(|(degree, _)| degree.clone()),
    );
}

#[cfg(test)]
fn ckyz_observed_support_domain_for_degrees(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    target_degrees: &[Vec<usize>],
) -> Result<CkyzMonomialDomain> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    validate_ckyz_target_degrees(target_degrees, rank)?;
    for term in local_intersection_terms {
        if term.first >= rank || term.second >= rank {
            return Err(Error::InvalidInput(
                "CKYZ local intersection term index is outside the relation rank".into(),
            ));
        }
    }

    let extraction_degrees = ckyz_cover_closed_target_degrees(target_degrees)?;
    let broad_domain = CkyzMonomialDomain::target_downset(&extraction_degrees, rank)?;
    let mut degree_set = BTreeSet::from([vec![0; rank]]);
    for degree in &extraction_degrees {
        degree_set.insert(degree.clone());
    }
    for coordinate in 0..rank {
        let mut degree = vec![0; rank];
        degree[coordinate] = 1;
        if broad_domain.contains(&degree) {
            degree_set.insert(degree);
        }
    }

    let alpha = compute_ckyz_log_period_corrections_domain(relations, &broad_domain)?;
    for series in &alpha {
        ckyz_record_series_support(&mut degree_set, series);
    }
    let z_of_q = compute_ckyz_inverse_mirror_map_domain(&alpha, &broad_domain)?;
    for series in &z_of_q {
        ckyz_record_series_support(&mut degree_set, series);
    }

    let mut contracted = BTreeMap::new();
    for term in local_intersection_terms {
        let beta = ckyz_second_log_period_series_for_pair_domain(
            relations,
            term.first,
            term.second,
            &broad_domain,
        )?;
        ckyz_record_series_support(&mut degree_set, &beta);

        let alpha_product =
            ckyz_series_mul_domain(&alpha[term.first], &alpha[term.second], &broad_domain)?;
        ckyz_record_series_support(&mut degree_set, &alpha_product);

        let mut f_pair = beta;
        ckyz_series_add_scaled_assign(&mut f_pair, &alpha_product, Rational::from(-1));
        ckyz_record_series_support(&mut degree_set, &f_pair);

        let mut term_coefficient = Rational::from(term.coefficient);
        if term.first == term.second {
            term_coefficient /= Rational::from(2);
        }
        ckyz_series_add_scaled_assign(&mut contracted, &f_pair, term_coefficient);
    }
    ckyz_record_series_support(&mut degree_set, &contracted);

    let potential =
        substitute_ckyz_series_in_flat_coordinates_domain(&contracted, &z_of_q, &broad_domain)?;
    ckyz_record_series_support(&mut degree_set, &potential);

    CkyzMonomialDomain::from_degrees(rank, degree_set)
}

fn ckyz_series_support_indices(
    series: &BTreeMap<Vec<usize>, Rational>,
    domain: &CkyzMonomialDomain,
) -> BTreeSet<usize> {
    series
        .iter()
        .filter(|(_, coefficient)| **coefficient != 0)
        .filter_map(|(degree, _)| domain.index_of(degree))
        .collect()
}

fn ckyz_log_period_support_indices_domain(
    relations: &[Vec<i64>],
    domain: &CkyzMonomialDomain,
) -> Result<Vec<BTreeSet<usize>>> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    if domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ monomial domain rank does not match relation rank".into(),
        ));
    }

    let mut supports = vec![BTreeSet::new(); rank];
    for (degree_index, degree) in domain.degrees.iter().enumerate() {
        if degree.iter().all(|&entry| entry == 0) {
            continue;
        }
        let point_pairings = ckyz_point_pairings(relations, degree)?;
        let negative_points = point_pairings
            .iter()
            .enumerate()
            .filter_map(|(index, &pairing)| (pairing < 0).then_some(index))
            .collect::<Vec<_>>();

        match negative_points.as_slice() {
            [negative_point] => {
                for coordinate_index in 0..rank {
                    if relations[coordinate_index][*negative_point] != 0 {
                        supports[coordinate_index].insert(degree_index);
                    }
                }
            }
            [] => {
                let regular_terms = ckyz_regular_harmonic_terms(relations, &point_pairings)?;
                for (coordinate_index, term) in regular_terms.iter().enumerate() {
                    if *term != 0 {
                        supports[coordinate_index].insert(degree_index);
                    }
                }
            }
            _ => {}
        }
    }

    Ok(supports)
}

fn ckyz_second_log_period_support_indices_for_pair_domain(
    relations: &[Vec<i64>],
    first: usize,
    second: usize,
    domain: &CkyzMonomialDomain,
) -> Result<BTreeSet<usize>> {
    let rank = relations.len();
    if first >= rank || second >= rank {
        return Err(Error::InvalidInput(
            "CKYZ second-log period pair index is outside the relation rank".into(),
        ));
    }
    if domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ componentwise cutoff rank does not match relation rank".into(),
        ));
    }

    let mut support = BTreeSet::new();
    for (degree_index, degree) in domain.degrees.iter().enumerate() {
        if degree.iter().all(|&entry| entry == 0) {
            continue;
        }
        let point_pairings = ckyz_point_pairings(relations, degree)?;
        let negative_points = point_pairings
            .iter()
            .enumerate()
            .filter_map(|(index, &pairing)| (pairing < 0).then_some(index))
            .collect::<Vec<_>>();

        let supported = match negative_points.as_slice() {
            [] => {
                let regular_terms = ckyz_regular_harmonic_terms(relations, &point_pairings)?;
                let mut term = regular_terms[first].clone() * regular_terms[second].clone();
                for point_index in 0..point_pairings.len() {
                    let arg = ckyz_harmonic_argument(point_pairings[point_index])?;
                    term += Rational::from(relations[first][point_index])
                        * Rational::from(relations[second][point_index])
                        * harmonic_number_order_two(arg);
                }
                term != 0
            }
            [negative_point] => {
                let regular_terms = ckyz_regular_harmonic_terms(relations, &point_pairings)?;
                regular_terms[first].clone() * Rational::from(relations[second][*negative_point])
                    + regular_terms[second].clone()
                        * Rational::from(relations[first][*negative_point])
                    != 0
            }
            [first_negative, second_negative] => {
                Rational::from(relations[first][*first_negative])
                    * Rational::from(relations[second][*second_negative])
                    + Rational::from(relations[first][*second_negative])
                        * Rational::from(relations[second][*first_negative])
                    != 0
            }
            _ => false,
        };
        if supported {
            support.insert(degree_index);
        }
    }

    Ok(support)
}

fn ckyz_support_identity_index(domain: &CkyzMonomialDomain) -> Result<usize> {
    domain.index_of(&vec![0; domain.rank]).ok_or_else(|| {
        Error::InvalidInput("CKYZ support domain is missing the identity monomial".into())
    })
}

fn ckyz_support_coordinate_series_indices(
    rank: usize,
    coordinate: usize,
    domain: &CkyzMonomialDomain,
) -> BTreeSet<usize> {
    let mut out = BTreeSet::new();
    if coordinate < rank && rank == domain.rank {
        let mut degree = vec![0; rank];
        degree[coordinate] = 1;
        if let Some(index) = domain.index_of(&degree) {
            out.insert(index);
        }
    }
    out
}

fn ckyz_support_min_total_degree(
    support: &BTreeSet<usize>,
    domain: &CkyzMonomialDomain,
    context: &str,
) -> Result<Option<usize>> {
    let mut min_degree = None;
    for &degree_index in support {
        let Some(degree) = domain.degrees.get(degree_index) else {
            return Err(Error::InvalidInput(format!(
                "{context} monomial index is outside the support domain"
            )));
        };
        let total_degree = ckyz_total_degree(degree)?;
        if total_degree == 0 {
            return Err(Error::InvalidInput(format!(
                "{context} must have zero constant term"
            )));
        }
        min_degree =
            Some(min_degree.map_or(total_degree, |current: usize| current.min(total_degree)));
    }
    Ok(min_degree)
}

fn ckyz_support_mul_domain(
    lhs: &BTreeSet<usize>,
    rhs: &BTreeSet<usize>,
    domain: &CkyzMonomialDomain,
) -> Result<BTreeSet<usize>> {
    let mut out = BTreeSet::new();
    for lhs_index in lhs {
        for rhs_index in rhs {
            if let Some(product_index) = domain.sum_index(*lhs_index, *rhs_index)? {
                out.insert(product_index);
            }
        }
    }
    Ok(out)
}

fn ckyz_support_power_cache_domain(
    support: &BTreeSet<usize>,
    max_exponent: usize,
    domain: &CkyzMonomialDomain,
) -> Result<Vec<BTreeSet<usize>>> {
    let mut powers = Vec::with_capacity(max_exponent + 1);
    powers.push(BTreeSet::from([ckyz_support_identity_index(domain)?]));
    for exponent in 1..=max_exponent {
        let next = ckyz_support_mul_domain(&powers[exponent - 1], support, domain)?;
        if next.is_empty() {
            powers.push(next);
            powers.resize_with(max_exponent + 1, BTreeSet::new);
            break;
        }
        powers.push(next);
    }
    Ok(powers)
}

fn ckyz_support_additive_closure_domain(
    support: &BTreeSet<usize>,
    domain: &CkyzMonomialDomain,
) -> Result<BTreeSet<usize>> {
    let identity = ckyz_support_identity_index(domain)?;
    let mut closure = BTreeSet::from([identity]);
    if support.is_empty() {
        return Ok(closure);
    }

    let mut generators = support.iter().copied().collect::<Vec<_>>();
    generators.sort_by(|lhs, rhs| {
        ckyz_total_degree(&domain.degrees[*lhs])
            .expect("domain degrees are validated")
            .cmp(&ckyz_total_degree(&domain.degrees[*rhs]).expect("domain degrees are validated"))
            .then_with(|| domain.degrees[*lhs].cmp(&domain.degrees[*rhs]))
    });

    for generator in generators {
        if generator >= domain.degrees.len() {
            return Err(Error::InvalidInput(
                "CKYZ support closure generator is outside the domain".into(),
            ));
        }
        if generator == identity || closure.contains(&generator) {
            continue;
        }

        let base_closure = closure.iter().copied().collect::<Vec<_>>();
        let mut frontier = VecDeque::new();
        for base in base_closure {
            if let Some(sum_index) = domain.sum_index(base, generator)? {
                if closure.insert(sum_index) {
                    frontier.push_back(sum_index);
                }
            }
        }
        while let Some(current) = frontier.pop_front() {
            if let Some(sum_index) = domain.sum_index(current, generator)? {
                if closure.insert(sum_index) {
                    frontier.push_back(sum_index);
                }
            }
        }
    }

    Ok(closure)
}

#[cfg(test)]
fn ckyz_support_exp_domain_by_powers(
    support: &BTreeSet<usize>,
    domain: &CkyzMonomialDomain,
) -> Result<BTreeSet<usize>> {
    let max_exponent =
        ckyz_support_min_total_degree(support, domain, "CKYZ support exponential input")?
            .map_or(0, |min_degree| domain.max_total_degree / min_degree);
    let mut out = BTreeSet::from([ckyz_support_identity_index(domain)?]);
    let mut power = out.clone();
    for _ in 1..=max_exponent {
        power = ckyz_support_mul_domain(&power, support, domain)?;
        if power.is_empty() {
            break;
        }
        out.extend(power.iter().copied());
    }
    Ok(out)
}

fn ckyz_support_exp_domain(
    support: &BTreeSet<usize>,
    domain: &CkyzMonomialDomain,
) -> Result<BTreeSet<usize>> {
    ckyz_support_min_total_degree(support, domain, "CKYZ support exponential input")?;
    ckyz_support_additive_closure_domain(support, domain)
}

fn ckyz_support_compose_domain(
    series_support: &BTreeSet<usize>,
    argument_supports: &[BTreeSet<usize>],
    domain: &CkyzMonomialDomain,
) -> Result<BTreeSet<usize>> {
    let rank = argument_supports.len();
    if rank != domain.rank {
        return Err(Error::InvalidInput(
            "CKYZ support composition rank mismatch".into(),
        ));
    }
    let mut max_exponents = vec![0usize; rank];
    for &degree_index in series_support {
        let Some(degree) = domain.degrees.get(degree_index) else {
            return Err(Error::InvalidInput(
                "CKYZ support composition index is outside the domain".into(),
            ));
        };
        for (coordinate, &exponent) in degree.iter().enumerate() {
            max_exponents[coordinate] = max_exponents[coordinate].max(exponent);
        }
    }
    let power_caches = argument_supports
        .iter()
        .zip(max_exponents.iter())
        .map(|(argument, &max_exponent)| {
            ckyz_support_power_cache_domain(argument, max_exponent, domain)
        })
        .collect::<Result<Vec<_>>>()?;

    let mut out = BTreeSet::new();
    for &degree_index in series_support {
        let Some(degree) = domain.degrees.get(degree_index) else {
            return Err(Error::InvalidInput(
                "CKYZ support composition index is outside the domain".into(),
            ));
        };
        let mut monomial = BTreeSet::from([ckyz_support_identity_index(domain)?]);
        for (coordinate, &exponent) in degree.iter().enumerate() {
            if exponent == 0 {
                continue;
            }
            monomial =
                ckyz_support_mul_domain(&monomial, &power_caches[coordinate][exponent], domain)?;
            if monomial.is_empty() {
                break;
            }
        }
        out.extend(monomial);
    }
    Ok(out)
}

fn ckyz_inverse_mirror_map_support_domain(
    alpha_supports: &[BTreeSet<usize>],
    domain: &CkyzMonomialDomain,
) -> Result<Vec<BTreeSet<usize>>> {
    let rank = alpha_supports.len();
    if rank == 0 || rank != domain.rank {
        return Err(Error::InvalidInput(
            "CKYZ support inverse mirror map rank mismatch".into(),
        ));
    }
    let trace_timing = env::var_os("CYRUS_TRACE_CKYZ_SUPPORT_DOMAIN").is_some();
    let first_iteration_start = trace_timing.then(std::time::Instant::now);
    let exp_negative_corrections = alpha_supports
        .iter()
        .map(|alpha_support| ckyz_support_exp_domain(alpha_support, domain))
        .collect::<Result<Vec<_>>>()?;
    let mut z_of_q = Vec::with_capacity(rank);
    for (coordinate, exp_negative_correction) in exp_negative_corrections.iter().enumerate() {
        let q_coordinate = ckyz_support_coordinate_series_indices(rank, coordinate, domain);
        z_of_q.push(ckyz_support_mul_domain(
            &q_coordinate,
            exp_negative_correction,
            domain,
        )?);
    }
    if let Some(start) = first_iteration_start {
        let sizes = z_of_q.iter().map(BTreeSet::len).collect::<Vec<_>>();
        eprintln!(
            "[CKYZ_SUPPORT_Z] iteration=1 sizes={sizes:?} elapsed={:?}",
            start.elapsed(),
        );
    }
    // If every alpha coordinate generates the same support semigroup, substituting
    // z_i = q_i * S into any alpha_j cannot leave S and still contains alpha_j.
    if exp_negative_corrections
        .windows(2)
        .all(|window| window[0] == window[1])
    {
        return Ok(z_of_q);
    }
    for iteration in 1..domain.max_total_degree {
        let iteration_start = trace_timing.then(std::time::Instant::now);
        let mut next = Vec::with_capacity(rank);
        for (coordinate, alpha_support) in alpha_supports.iter().enumerate() {
            let correction_at_z = ckyz_support_compose_domain(alpha_support, &z_of_q, domain)?;
            let exp_negative_correction = ckyz_support_exp_domain(&correction_at_z, domain)?;
            let q_coordinate = ckyz_support_coordinate_series_indices(rank, coordinate, domain);
            next.push(ckyz_support_mul_domain(
                &q_coordinate,
                &exp_negative_correction,
                domain,
            )?);
        }
        if let Some(start) = iteration_start {
            let sizes = next.iter().map(BTreeSet::len).collect::<Vec<_>>();
            eprintln!(
                "[CKYZ_SUPPORT_Z] iteration={} sizes={sizes:?} elapsed={:?}",
                iteration + 1,
                start.elapsed(),
            );
        }
        if next == z_of_q {
            break;
        }
        z_of_q = next;
    }
    Ok(z_of_q)
}

fn ckyz_predicted_support_domain_for_degrees(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    target_degrees: &[Vec<usize>],
) -> Result<CkyzMonomialDomain> {
    let trace_timing = env::var_os("CYRUS_TRACE_CKYZ_SUPPORT_DOMAIN").is_some();
    let trace_start = trace_timing.then(std::time::Instant::now);
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    validate_ckyz_target_degrees(target_degrees, rank)?;
    for term in local_intersection_terms {
        if term.first >= rank || term.second >= rank {
            return Err(Error::InvalidInput(
                "CKYZ local intersection term index is outside the relation rank".into(),
            ));
        }
    }

    let extraction_degrees = ckyz_cover_closed_target_degrees(target_degrees)?;
    let broad_domain = CkyzMonomialDomain::target_downset(&extraction_degrees, rank)?;
    let trace_after_domain = trace_timing.then(std::time::Instant::now);
    let mut degree_set = BTreeSet::from([vec![0; rank]]);
    for degree in &extraction_degrees {
        degree_set.insert(degree.clone());
    }
    for coordinate in 0..rank {
        let mut degree = vec![0; rank];
        degree[coordinate] = 1;
        if broad_domain.contains(&degree) {
            degree_set.insert(degree);
        }
    }

    let alpha_supports = ckyz_log_period_support_indices_domain(relations, &broad_domain)?;
    let trace_after_alpha = trace_timing.then(std::time::Instant::now);
    for support in &alpha_supports {
        degree_set.extend(
            support
                .iter()
                .map(|&index| broad_domain.degrees[index].clone()),
        );
    }
    let z_of_q_supports = ckyz_inverse_mirror_map_support_domain(&alpha_supports, &broad_domain)?;
    let trace_after_z = trace_timing.then(std::time::Instant::now);
    for support in &z_of_q_supports {
        degree_set.extend(
            support
                .iter()
                .map(|&index| broad_domain.degrees[index].clone()),
        );
    }

    let mut alpha_support_classes = Vec::with_capacity(rank);
    let mut unique_alpha_supports = Vec::<&BTreeSet<usize>>::new();
    for support in &alpha_supports {
        let class = unique_alpha_supports
            .iter()
            .position(|candidate| *candidate == support)
            .unwrap_or_else(|| {
                unique_alpha_supports.push(support);
                unique_alpha_supports.len() - 1
            });
        alpha_support_classes.push(class);
    }
    let mut alpha_product_cache = BTreeMap::<(usize, usize), BTreeSet<usize>>::new();
    let mut contracted_support = BTreeSet::new();
    for term in local_intersection_terms {
        let beta_support = ckyz_second_log_period_support_indices_for_pair_domain(
            relations,
            term.first,
            term.second,
            &broad_domain,
        )?;
        degree_set.extend(
            beta_support
                .iter()
                .map(|&index| broad_domain.degrees[index].clone()),
        );

        let mut alpha_product_key = (
            alpha_support_classes[term.first],
            alpha_support_classes[term.second],
        );
        if alpha_product_key.0 > alpha_product_key.1 {
            alpha_product_key = (alpha_product_key.1, alpha_product_key.0);
        }
        let alpha_product_support =
            if let Some(cached_support) = alpha_product_cache.get(&alpha_product_key) {
                cached_support.clone()
            } else {
                let support = ckyz_support_mul_domain(
                    unique_alpha_supports[alpha_product_key.0],
                    unique_alpha_supports[alpha_product_key.1],
                    &broad_domain,
                )?;
                alpha_product_cache.insert(alpha_product_key, support.clone());
                support
            };
        degree_set.extend(
            alpha_product_support
                .iter()
                .map(|&index| broad_domain.degrees[index].clone()),
        );

        contracted_support.extend(beta_support);
        contracted_support.extend(alpha_product_support);
    }
    let trace_after_contracted = trace_timing.then(std::time::Instant::now);
    degree_set.extend(
        contracted_support
            .iter()
            .map(|&index| broad_domain.degrees[index].clone()),
    );

    let trace_after_potential = trace_timing.then(std::time::Instant::now);
    if let (
        Some(start),
        Some(after_domain),
        Some(after_alpha),
        Some(after_z),
        Some(after_contracted),
        Some(after_potential),
    ) = (
        trace_start,
        trace_after_domain,
        trace_after_alpha,
        trace_after_z,
        trace_after_contracted,
        trace_after_potential,
    ) {
        let alpha_sizes = alpha_supports.iter().map(BTreeSet::len).collect::<Vec<_>>();
        let z_sizes = z_of_q_supports
            .iter()
            .map(BTreeSet::len)
            .collect::<Vec<_>>();
        eprintln!(
            "[CKYZ_SUPPORT_DOMAIN] rank={rank} broad={} selected={} alpha_sizes={alpha_sizes:?} z_sizes={z_sizes:?} contracted_size={} domain={:?} alpha={:?} z={:?} contracted={:?} potential=skipped({:?})",
            broad_domain.degrees.len(),
            degree_set.len(),
            contracted_support.len(),
            after_domain.duration_since(start),
            after_alpha.duration_since(after_domain),
            after_z.duration_since(after_alpha),
            after_contracted.duration_since(after_z),
            after_potential.duration_since(after_contracted),
        );
    }

    CkyzMonomialDomain::from_degrees(rank, degree_set)
}

fn ckyz_grading_degree(degree: &[usize], grading_vector: &[usize]) -> Result<usize> {
    if degree.len() != grading_vector.len() {
        return Err(Error::InvalidInput(
            "CKYZ grading degree rank mismatch".into(),
        ));
    }
    degree
        .iter()
        .zip(grading_vector.iter())
        .try_fold(0usize, |sum, (&entry, &weight)| {
            let term = entry.checked_mul(weight).ok_or_else(|| {
                Error::InvalidInput("CKYZ grading degree multiplication overflowed usize".into())
            })?;
            sum.checked_add(term).ok_or_else(|| {
                Error::InvalidInput("CKYZ grading degree addition overflowed usize".into())
            })
        })
}

fn ckyz_sum_degree_index(
    lhs_degree: &[usize],
    rhs_degree: &[usize],
    rank: usize,
    degree_indices: &HashMap<Vec<usize>, usize>,
) -> Result<Option<usize>> {
    if lhs_degree.len() != rank || rhs_degree.len() != rank {
        return Err(Error::InvalidInput(
            "CKYZ series multiplication rank mismatch".into(),
        ));
    }
    let mut degree = Vec::with_capacity(rank);
    for (&lhs_entry, &rhs_entry) in lhs_degree.iter().zip(rhs_degree.iter()) {
        degree.push(lhs_entry.checked_add(rhs_entry).ok_or_else(|| {
            Error::InvalidInput("CKYZ multidegree addition overflowed usize".into())
        })?);
    }
    Ok(degree_indices.get(&degree).copied())
}

fn ckyz_dense_degree_index(
    degree: &[usize],
    rank: usize,
    max_coordinate_degrees: &[usize],
    strides: &[usize],
) -> Result<Option<usize>> {
    if degree.len() != rank || max_coordinate_degrees.len() != rank || strides.len() != rank {
        return Err(Error::InvalidInput(
            "CKYZ dense degree index rank mismatch".into(),
        ));
    }
    let mut index = 0usize;
    for ((&entry, &max_degree), &stride) in degree
        .iter()
        .zip(max_coordinate_degrees.iter())
        .zip(strides.iter())
    {
        if entry > max_degree {
            return Ok(None);
        }
        index = index
            .checked_add(entry.checked_mul(stride).ok_or_else(|| {
                Error::InvalidInput("CKYZ dense degree index overflowed usize".into())
            })?)
            .ok_or_else(|| {
                Error::InvalidInput("CKYZ dense degree index overflowed usize".into())
            })?;
    }
    Ok(Some(index))
}

fn ckyz_dense_domain_degree_sum_index(
    lhs_degree: &[usize],
    rhs_degree: &[usize],
    rank: usize,
    max_coordinate_degrees: &[usize],
    strides: &[usize],
) -> Option<usize> {
    debug_assert_eq!(lhs_degree.len(), rank);
    debug_assert_eq!(rhs_degree.len(), rank);
    debug_assert_eq!(max_coordinate_degrees.len(), rank);
    debug_assert_eq!(strides.len(), rank);
    match rank {
        1 => {
            let degree0 = lhs_degree[0] + rhs_degree[0];
            (degree0 <= max_coordinate_degrees[0]).then_some(degree0 * strides[0])
        }
        2 => {
            let degree0 = lhs_degree[0] + rhs_degree[0];
            if degree0 > max_coordinate_degrees[0] {
                return None;
            }
            let degree1 = lhs_degree[1] + rhs_degree[1];
            if degree1 > max_coordinate_degrees[1] {
                return None;
            }
            Some(degree0 * strides[0] + degree1 * strides[1])
        }
        3 => {
            let degree0 = lhs_degree[0] + rhs_degree[0];
            if degree0 > max_coordinate_degrees[0] {
                return None;
            }
            let degree1 = lhs_degree[1] + rhs_degree[1];
            if degree1 > max_coordinate_degrees[1] {
                return None;
            }
            let degree2 = lhs_degree[2] + rhs_degree[2];
            if degree2 > max_coordinate_degrees[2] {
                return None;
            }
            Some(degree0 * strides[0] + degree1 * strides[1] + degree2 * strides[2])
        }
        _ => ckyz_dense_domain_degree_sum_index_general(
            lhs_degree,
            rhs_degree,
            rank,
            max_coordinate_degrees,
            strides,
        ),
    }
}

fn ckyz_dense_domain_degree_sum_index_general(
    lhs_degree: &[usize],
    rhs_degree: &[usize],
    rank: usize,
    max_coordinate_degrees: &[usize],
    strides: &[usize],
) -> Option<usize> {
    let mut index = 0usize;
    for coordinate in 0..rank {
        let entry = lhs_degree[coordinate] + rhs_degree[coordinate];
        if entry > max_coordinate_degrees[coordinate] {
            return None;
        }
        index += entry * strides[coordinate];
    }
    Some(index)
}

fn ckyz_add_degrees_unbounded(
    lhs_degree: &[usize],
    rhs_degree: &[usize],
    rank: usize,
) -> Result<Vec<usize>> {
    if lhs_degree.len() != rank || rhs_degree.len() != rank {
        return Err(Error::InvalidInput(
            "CKYZ series multiplication rank mismatch".into(),
        ));
    }
    lhs_degree
        .iter()
        .zip(rhs_degree.iter())
        .map(|(&lhs_entry, &rhs_entry)| {
            lhs_entry.checked_add(rhs_entry).ok_or_else(|| {
                Error::InvalidInput("CKYZ multidegree addition overflowed usize".into())
            })
        })
        .collect()
}

fn ckyz_coordinate_series(
    rank: usize,
    coordinate: usize,
    max_total_degree: usize,
) -> BTreeMap<Vec<usize>, Rational> {
    let mut out = BTreeMap::new();
    if max_total_degree > 0 {
        let mut degree = vec![0; rank];
        degree[coordinate] = 1;
        out.insert(degree, Rational::from(1));
    }
    out
}

fn ckyz_series_scale(
    series: &BTreeMap<Vec<usize>, Rational>,
    scalar: Rational,
) -> BTreeMap<Vec<usize>, Rational> {
    if scalar == 0 {
        return BTreeMap::new();
    }
    series
        .iter()
        .filter_map(|(degree, coefficient)| {
            let scaled = coefficient.clone() * scalar.clone();
            (scaled != 0).then(|| (degree.clone(), scaled))
        })
        .collect()
}

fn ckyz_series_add_scaled_assign(
    out: &mut BTreeMap<Vec<usize>, Rational>,
    series: &BTreeMap<Vec<usize>, Rational>,
    scalar: Rational,
) {
    if scalar == 0 {
        return;
    }
    let mut zero_degrees = Vec::new();
    for (degree, coefficient) in series {
        let entry = out
            .entry(degree.clone())
            .or_insert_with(|| Rational::from(0));
        *entry += coefficient.clone() * scalar.clone();
        if *entry == 0 {
            zero_degrees.push(degree.clone());
        }
    }
    for degree in zero_degrees {
        out.remove(&degree);
    }
}

fn ckyz_series_min_total_degree(
    series: &BTreeMap<Vec<usize>, Rational>,
    rank: usize,
    context: &str,
) -> Result<Option<usize>> {
    let mut min_degree = None;
    for (degree, coefficient) in series {
        if *coefficient == 0 {
            continue;
        }
        if degree.len() != rank {
            return Err(Error::InvalidInput(format!(
                "{context} monomial rank does not match coordinate rank"
            )));
        }
        let total_degree = ckyz_total_degree(degree)?;
        if total_degree == 0 {
            return Err(Error::InvalidInput(format!(
                "{context} must have zero constant term"
            )));
        }
        min_degree =
            Some(min_degree.map_or(total_degree, |current: usize| current.min(total_degree)));
    }
    Ok(min_degree)
}

fn ckyz_series_mul(
    lhs: &BTreeMap<Vec<usize>, Rational>,
    rhs: &BTreeMap<Vec<usize>, Rational>,
    rank: usize,
    max_total_degree: usize,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let mut out = BTreeMap::new();
    for (lhs_degree, lhs_coefficient) in lhs {
        if *lhs_coefficient == 0 {
            continue;
        }
        for (rhs_degree, rhs_coefficient) in rhs {
            if *rhs_coefficient == 0 {
                continue;
            }
            let product_degree = ckyz_add_degrees(lhs_degree, rhs_degree, rank, max_total_degree)?;
            let Some(product_degree) = product_degree else {
                continue;
            };
            let entry = out
                .entry(product_degree)
                .or_insert_with(|| Rational::from(0));
            *entry += lhs_coefficient.clone() * rhs_coefficient.clone();
        }
    }
    out.retain(|_, coefficient| *coefficient != 0);
    Ok(out)
}

fn ckyz_series_mul_domain(
    lhs: &BTreeMap<Vec<usize>, Rational>,
    rhs: &BTreeMap<Vec<usize>, Rational>,
    domain: &CkyzMonomialDomain,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let lhs = CkyzIndexedSeries::from_btree(lhs, domain, "CKYZ series multiplication")?;
    let rhs = CkyzIndexedSeries::from_btree(rhs, domain, "CKYZ series multiplication")?;
    Ok(lhs.mul(&rhs, domain)?.to_btree(domain))
}

fn ckyz_add_degrees(
    lhs: &[usize],
    rhs: &[usize],
    rank: usize,
    max_total_degree: usize,
) -> Result<Option<Vec<usize>>> {
    if lhs.len() != rank || rhs.len() != rank {
        return Err(Error::InvalidInput(
            "CKYZ series multiplication rank mismatch".into(),
        ));
    }
    let mut degree = Vec::with_capacity(rank);
    for (&lhs_entry, &rhs_entry) in lhs.iter().zip(rhs.iter()) {
        degree.push(lhs_entry.checked_add(rhs_entry).ok_or_else(|| {
            Error::InvalidInput("CKYZ multidegree addition overflowed usize".into())
        })?);
    }
    if ckyz_total_degree(&degree)? > max_total_degree {
        return Ok(None);
    }
    Ok(Some(degree))
}

fn ckyz_series_power_cache(
    series: &BTreeMap<Vec<usize>, Rational>,
    max_exponent: usize,
    rank: usize,
    max_total_degree: usize,
) -> Result<Vec<BTreeMap<Vec<usize>, Rational>>> {
    let mut powers = Vec::with_capacity(max_exponent + 1);
    let mut identity = BTreeMap::new();
    identity.insert(vec![0; rank], Rational::from(1));
    powers.push(identity);
    for exponent in 1..=max_exponent {
        let next = ckyz_series_mul(&powers[exponent - 1], series, rank, max_total_degree)?;
        if next.is_empty() {
            powers.push(next);
            powers.resize_with(max_exponent + 1, BTreeMap::new);
            break;
        }
        powers.push(next);
    }
    Ok(powers)
}

#[cfg(test)]
fn ckyz_series_power_cache_domain(
    series: &BTreeMap<Vec<usize>, Rational>,
    max_exponent: usize,
    domain: &CkyzMonomialDomain,
) -> Result<Vec<BTreeMap<Vec<usize>, Rational>>> {
    let series = CkyzIndexedSeries::from_btree(series, domain, "CKYZ power-cache input")?;
    series
        .power_cache(max_exponent, domain)?
        .into_iter()
        .map(|power| Ok(power.to_btree(domain)))
        .collect()
}

fn ckyz_series_exp(
    series: &BTreeMap<Vec<usize>, Rational>,
    rank: usize,
    max_total_degree: usize,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    validate_ckyz_series_has_zero_constant(series, rank, "CKYZ exponential input")?;
    let max_exponent = ckyz_series_min_total_degree(series, rank, "CKYZ exponential input")?
        .map_or(0, |min_degree| max_total_degree / min_degree);

    let mut out = BTreeMap::new();
    out.insert(vec![0; rank], Rational::from(1));
    let mut power = out.clone();
    let mut factorial = Integer::from(1);
    for exponent in 1..=max_exponent {
        power = ckyz_series_mul(&power, series, rank, max_total_degree)?;
        if power.is_empty() {
            break;
        }
        factorial *= Integer::from(exponent);
        ckyz_series_add_scaled_assign(
            &mut out,
            &power,
            Rational::from(1) / Rational::from(factorial.clone()),
        );
    }
    Ok(out)
}

#[cfg(test)]
fn ckyz_series_exp_domain(
    series: &BTreeMap<Vec<usize>, Rational>,
    domain: &CkyzMonomialDomain,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let series = CkyzIndexedSeries::from_btree(series, domain, "CKYZ exponential input")?;
    Ok(series.exp(domain)?.to_btree(domain))
}

fn ckyz_series_compose(
    series: &BTreeMap<Vec<usize>, Rational>,
    arguments: &[BTreeMap<Vec<usize>, Rational>],
    max_total_degree: usize,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let rank = arguments.len();
    let mut max_exponents = vec![0usize; rank];
    for (degree, coefficient) in series {
        if *coefficient == 0 || ckyz_total_degree(degree)? > max_total_degree {
            continue;
        }
        if degree.len() != rank {
            return Err(Error::InvalidInput(
                "CKYZ series composition rank mismatch".into(),
            ));
        }
        for (coordinate, &exponent) in degree.iter().enumerate() {
            max_exponents[coordinate] = max_exponents[coordinate].max(exponent);
        }
    }
    let power_caches = arguments
        .iter()
        .zip(max_exponents.iter())
        .map(|(argument, &max_exponent)| {
            ckyz_series_power_cache(argument, max_exponent, rank, max_total_degree)
        })
        .collect::<Result<Vec<_>>>()?;

    let mut out = BTreeMap::new();
    for (degree, coefficient) in series {
        if *coefficient == 0 || ckyz_total_degree(degree)? > max_total_degree {
            continue;
        }
        if degree.len() != rank {
            return Err(Error::InvalidInput(
                "CKYZ series composition rank mismatch".into(),
            ));
        }
        let mut monomial = BTreeMap::new();
        monomial.insert(vec![0; rank], Rational::from(1));
        for (coordinate, &exponent) in degree.iter().enumerate() {
            if exponent == 0 {
                continue;
            }
            monomial = ckyz_series_mul(
                &monomial,
                &power_caches[coordinate][exponent],
                rank,
                max_total_degree,
            )?;
            if monomial.is_empty() {
                break;
            }
        }
        ckyz_series_add_scaled_assign(&mut out, &monomial, coefficient.clone());
    }
    Ok(out)
}

fn ckyz_series_compose_domain(
    series: &BTreeMap<Vec<usize>, Rational>,
    arguments: &[BTreeMap<Vec<usize>, Rational>],
    domain: &CkyzMonomialDomain,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let rank = arguments.len();
    if rank != domain.rank {
        return Err(Error::InvalidInput(
            "CKYZ series composition rank mismatch".into(),
        ));
    }
    let series = CkyzIndexedSeries::from_btree(series, domain, "CKYZ series composition")?;
    let arguments = arguments
        .iter()
        .map(|argument| {
            CkyzIndexedSeries::from_btree(argument, domain, "CKYZ series composition argument")
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(series.compose(&arguments, domain)?.to_btree(domain))
}

fn ckyz_nontrivial_covers(degree: &[usize]) -> Vec<usize> {
    let max_cover = degree.iter().copied().max().unwrap_or(0);
    (2..=max_cover)
        .filter(|&cover| degree.iter().all(|entry| entry % cover == 0))
        .collect()
}

fn ckyz_divide_degree(degree: &[usize], divisor: usize) -> Result<Vec<usize>> {
    if divisor == 0 {
        return Err(Error::InvalidInput(
            "CKYZ degree divisor must be nonzero".into(),
        ));
    }
    Ok(degree.iter().map(|entry| entry / divisor).collect())
}

fn ckyz_cover_weight(coefficients: &[i64], degree: &[usize]) -> Result<i64> {
    if coefficients.len() != degree.len() {
        return Err(Error::InvalidInput(
            "CKYZ cover-weight rank does not match degree rank".into(),
        ));
    }
    let mut dot = 0i64;
    for (&coefficient, &degree_entry) in coefficients.iter().zip(degree.iter()) {
        let degree_entry = i64::try_from(degree_entry).map_err(|_| {
            Error::InvalidInput("CKYZ cover-weight degree entry does not fit i64".into())
        })?;
        let contribution = coefficient.checked_mul(degree_entry).ok_or_else(|| {
            Error::InvalidInput("CKYZ cover-weight multiplication overflowed i64".into())
        })?;
        dot = dot.checked_add(contribution).ok_or_else(|| {
            Error::InvalidInput("CKYZ cover-weight addition overflowed i64".into())
        })?;
    }
    dot.checked_neg()
        .ok_or_else(|| Error::InvalidInput("CKYZ cover weight overflowed i64".into()))
}

fn ckyz_grading_vector_from_cover_weights(coefficients: &[i64]) -> Result<Vec<usize>> {
    if coefficients.is_empty() {
        return Err(Error::InvalidInput(
            "CKYZ grading vector requires at least one cover weight".into(),
        ));
    }
    coefficients
        .iter()
        .map(|&coefficient| {
            if coefficient <= 0 {
                return Err(Error::InvalidInput(
                    "CKYZ grading weights must be positive".into(),
                ));
            }
            usize::try_from(coefficient)
                .map_err(|_| Error::InvalidInput("CKYZ grading weight does not fit usize".into()))
        })
        .collect()
}

fn validate_ckyz_target_degrees(target_degrees: &[Vec<usize>], rank: usize) -> Result<()> {
    if target_degrees.is_empty() {
        return Err(Error::InvalidInput(
            "CKYZ targeted GV extraction requires at least one target degree".into(),
        ));
    }
    for degree in target_degrees {
        if degree.len() != rank {
            return Err(Error::InvalidInput(
                "CKYZ target degree rank does not match relation rank".into(),
            ));
        }
        if degree.iter().all(|&entry| entry == 0) {
            return Err(Error::InvalidInput(
                "CKYZ target degree must be nonzero".into(),
            ));
        }
        ckyz_total_degree(degree)?;
    }
    Ok(())
}

fn ckyz_cover_closed_target_degrees(target_degrees: &[Vec<usize>]) -> Result<Vec<Vec<usize>>> {
    let mut out = BTreeSet::new();
    for degree in target_degrees {
        out.insert(degree.clone());
        for cover in ckyz_nontrivial_covers(degree) {
            out.insert(ckyz_divide_degree(degree, cover)?);
        }
    }
    let mut degrees = out.into_iter().collect::<Vec<_>>();
    ckyz_sort_degrees_for_extraction(&mut degrees);
    Ok(degrees)
}

fn ckyz_sort_degrees_for_extraction(degrees: &mut [Vec<usize>]) {
    degrees.sort_by(|lhs, rhs| {
        ckyz_total_degree(lhs)
            .expect("validated degree")
            .cmp(&ckyz_total_degree(rhs).expect("validated degree"))
            .then_with(|| lhs.cmp(rhs))
    });
}

fn ckyz_sort_degrees_for_extraction_with_grading(
    degrees: &mut [Vec<usize>],
    grading_vector: &[usize],
) -> Result<()> {
    if grading_vector.is_empty() {
        return Err(Error::InvalidInput(
            "CKYZ extraction grading requires at least one weight".into(),
        ));
    }
    for degree in degrees.iter() {
        if degree.len() != grading_vector.len() {
            return Err(Error::InvalidInput(
                "CKYZ extraction grading rank does not match degree rank".into(),
            ));
        }
        ckyz_grading_degree(degree, grading_vector)?;
        ckyz_total_degree(degree)?;
    }
    degrees.sort_by(|lhs, rhs| {
        ckyz_grading_degree(lhs, grading_vector)
            .expect("validated degree")
            .cmp(&ckyz_grading_degree(rhs, grading_vector).expect("validated degree"))
            .then_with(|| {
                ckyz_total_degree(lhs)
                    .expect("validated degree")
                    .cmp(&ckyz_total_degree(rhs).expect("validated degree"))
            })
            .then_with(|| lhs.cmp(rhs))
    });
    Ok(())
}

fn ckyz_cygv_previous_qn_level_count(rank: usize) -> usize {
    if rank < 4 {
        2
    } else if rank < 10 {
        5
    } else {
        10
    }
}

#[allow(clippy::cast_precision_loss)]
fn ckyz_qn_history_distance_score(degree: &[usize]) -> f64 {
    degree
        .iter()
        .map(|&entry| {
            if entry == 0 {
                0.0
            } else {
                (entry as f64).log2() + 1.0
            }
        })
        .sum()
}

fn ckyz_qn_history_reuse_profile(
    extraction_degrees: &[Vec<usize>],
    extraction_gradings: &[usize],
    domain: &CkyzMonomialDomain,
) -> Result<CkyzQnHistoryReuseProfile> {
    if extraction_degrees.len() != extraction_gradings.len() {
        return Err(Error::InvalidInput(
            "CKYZ q_N history profile degree/grading length mismatch".into(),
        ));
    }
    let level_count = ckyz_cygv_previous_qn_level_count(domain.rank);
    let mut previous_degrees = VecDeque::from(vec![Vec::<Vec<usize>>::new(); level_count]);
    let mut candidate_hit_count = 0usize;
    let mut candidate_miss_count = 0usize;
    let mut unique_deltas = BTreeSet::<Vec<usize>>::new();

    let mut batch_start = 0usize;
    while batch_start < extraction_degrees.len() {
        let batch_grading = extraction_gradings[batch_start];
        let batch_end = extraction_gradings.partition_point(|&grading| grading <= batch_grading);
        for degree in &extraction_degrees[batch_start..batch_end] {
            let mut closest_delta = degree.clone();
            let mut closest_distance = ckyz_qn_history_distance_score(degree);
            let mut hit_previous = false;
            for previous_batch in &previous_degrees {
                for previous_degree in previous_batch {
                    let Some(delta) = ckyz_subtract_degree_multiple(degree, previous_degree, 1)
                    else {
                        continue;
                    };
                    if !domain.contains(&delta) {
                        continue;
                    }
                    let distance = ckyz_qn_history_distance_score(&delta);
                    if distance < closest_distance {
                        closest_distance = distance;
                        closest_delta = delta;
                        hit_previous = true;
                    }
                }
            }
            if hit_previous {
                candidate_hit_count = candidate_hit_count.checked_add(1).ok_or_else(|| {
                    Error::InvalidInput("CKYZ q_N history hit count overflowed".into())
                })?;
            } else {
                candidate_miss_count = candidate_miss_count.checked_add(1).ok_or_else(|| {
                    Error::InvalidInput("CKYZ q_N history miss count overflowed".into())
                })?;
            }
            if closest_delta.iter().any(|&entry| entry != 0) {
                unique_deltas.insert(closest_delta);
            }
        }
        previous_degrees.pop_front();
        previous_degrees.push_back(extraction_degrees[batch_start..batch_end].to_vec());
        batch_start = batch_end;
    }

    Ok(CkyzQnHistoryReuseProfile {
        level_count,
        candidate_hit_count,
        candidate_miss_count,
        unique_delta_count: unique_deltas.len(),
        unique_delta_degrees: unique_deltas.into_iter().collect(),
    })
}

fn compute_ckyz_log_period_corrections_domain(
    relations: &[Vec<i64>],
    domain: &CkyzMonomialDomain,
) -> Result<Vec<BTreeMap<Vec<usize>, Rational>>> {
    validate_ckyz_relations(relations)?;
    if domain.rank != relations.len() {
        return Err(Error::InvalidInput(
            "CKYZ monomial domain rank does not match relation rank".into(),
        ));
    }

    let mut corrections = vec![BTreeMap::new(); domain.rank];
    for degree in domain.nonzero_degrees() {
        let point_pairings = ckyz_point_pairings(relations, degree)?;
        let values = ckyz_log_period_coefficients_for_degree(relations, &point_pairings)?;
        for (coordinate_index, value) in values.into_iter().enumerate() {
            if value != 0 {
                corrections[coordinate_index].insert(degree.clone(), value);
            }
        }
    }
    Ok(corrections)
}

fn compute_ckyz_inverse_mirror_map_domain(
    log_period_corrections: &[BTreeMap<Vec<usize>, Rational>],
    domain: &CkyzMonomialDomain,
) -> Result<Vec<BTreeMap<Vec<usize>, Rational>>> {
    let rank = log_period_corrections.len();
    if rank == 0 {
        return Err(Error::InvalidInput(
            "CKYZ inverse mirror map requires at least one coordinate".into(),
        ));
    }
    if domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ componentwise cutoff rank does not match relation rank".into(),
        ));
    }
    for correction in log_period_corrections {
        validate_ckyz_series(correction, rank, "CKYZ log-period correction")?;
        validate_ckyz_series_has_zero_constant(correction, rank, "CKYZ log-period correction")?;
    }

    let corrections = log_period_corrections
        .iter()
        .map(|correction| {
            CkyzIndexedSeries::from_btree(correction, domain, "CKYZ log-period correction")
        })
        .collect::<Result<Vec<_>>>()?;

    let mut z_of_q = (0..rank)
        .map(|coordinate| CkyzIndexedSeries::coordinate(rank, coordinate, domain))
        .collect::<Vec<_>>();
    for _ in 0..domain.max_total_degree {
        let mut next = Vec::with_capacity(rank);
        for (coordinate, correction) in corrections.iter().enumerate() {
            let correction_at_z = correction.compose(&z_of_q, domain)?;
            let negative_correction = correction_at_z.scaled(Rational::from(-1));
            let exp_negative_correction = negative_correction.exp(domain)?;
            let q_coordinate = CkyzIndexedSeries::coordinate(rank, coordinate, domain);
            next.push(q_coordinate.mul(&exp_negative_correction, domain)?);
        }
        if next == z_of_q {
            break;
        }
        z_of_q = next;
    }
    Ok(z_of_q
        .into_iter()
        .map(|coordinate| coordinate.to_btree(domain))
        .collect())
}

fn substitute_ckyz_series_in_flat_coordinates_domain(
    series_z: &BTreeMap<Vec<usize>, Rational>,
    z_of_q: &[BTreeMap<Vec<usize>, Rational>],
    domain: &CkyzMonomialDomain,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let rank = z_of_q.len();
    if rank == 0 {
        return Err(Error::InvalidInput(
            "CKYZ flat-coordinate substitution requires at least one coordinate".into(),
        ));
    }
    if domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ componentwise cutoff rank does not match relation rank".into(),
        ));
    }
    validate_ckyz_series(series_z, rank, "CKYZ B-model series")?;
    for argument in z_of_q {
        validate_ckyz_series(argument, rank, "CKYZ inverse mirror-map coordinate")?;
        validate_ckyz_series_has_zero_constant(
            argument,
            rank,
            "CKYZ inverse mirror-map coordinate",
        )?;
    }
    ckyz_series_compose_domain(series_z, z_of_q, domain)
}

fn compute_ckyz_local_instanton_potential_corrections_domain(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    domain: &CkyzMonomialDomain,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    if domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ componentwise cutoff rank does not match relation rank".into(),
        ));
    }
    for term in local_intersection_terms {
        if term.first >= rank || term.second >= rank {
            return Err(Error::InvalidInput(
                "CKYZ local intersection term index is outside the relation rank".into(),
            ));
        }
    }

    let alpha = compute_ckyz_log_period_corrections_domain(relations, domain)?;
    let z_of_q = compute_ckyz_inverse_mirror_map_domain(&alpha, domain)?;
    let mut contracted = BTreeMap::new();

    for term in local_intersection_terms {
        let beta = ckyz_second_log_period_series_for_pair_domain(
            relations,
            term.first,
            term.second,
            domain,
        )?;
        let alpha_product =
            ckyz_series_mul_domain(&alpha[term.first], &alpha[term.second], domain)?;
        let mut f_pair = beta;
        ckyz_series_add_scaled_assign(&mut f_pair, &alpha_product, Rational::from(-1));

        let mut term_coefficient = Rational::from(term.coefficient);
        if term.first == term.second {
            term_coefficient /= Rational::from(2);
        }
        ckyz_series_add_scaled_assign(&mut contracted, &f_pair, term_coefficient);
    }

    substitute_ckyz_series_in_flat_coordinates_domain(&contracted, &z_of_q, domain)
}

fn compute_ckyz_local_instanton_potential_z_domain(
    relations: &[Vec<i64>],
    local_intersection_terms: &[CkyzLocalIntersectionTerm],
    domain: &CkyzMonomialDomain,
) -> Result<(
    Vec<BTreeMap<Vec<usize>, Rational>>,
    BTreeMap<Vec<usize>, Rational>,
)> {
    validate_ckyz_relations(relations)?;
    let rank = relations.len();
    if domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ componentwise cutoff rank does not match relation rank".into(),
        ));
    }
    for term in local_intersection_terms {
        if term.first >= rank || term.second >= rank {
            return Err(Error::InvalidInput(
                "CKYZ local intersection term index is outside the relation rank".into(),
            ));
        }
    }

    let alpha = compute_ckyz_log_period_corrections_domain(relations, domain)?;
    let mut contracted = BTreeMap::new();

    for term in local_intersection_terms {
        let beta = ckyz_second_log_period_series_for_pair_domain(
            relations,
            term.first,
            term.second,
            domain,
        )?;
        let alpha_product =
            ckyz_series_mul_domain(&alpha[term.first], &alpha[term.second], domain)?;
        let mut f_pair = beta;
        ckyz_series_add_scaled_assign(&mut f_pair, &alpha_product, Rational::from(-1));

        let mut term_coefficient = Rational::from(term.coefficient);
        if term.first == term.second {
            term_coefficient /= Rational::from(2);
        }
        ckyz_series_add_scaled_assign(&mut contracted, &f_pair, term_coefficient);
    }

    Ok((alpha, contracted))
}

#[cfg(test)]
fn ckyz_q_degree_series_in_z_domain(
    degree: &[usize],
    alpha: &[BTreeMap<Vec<usize>, Rational>],
    domain: &CkyzMonomialDomain,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let rank = alpha.len();
    if rank == 0 || degree.len() != rank || domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ q-degree series rank mismatch".into(),
        ));
    }
    if !domain.contains(degree) {
        return Ok(BTreeMap::new());
    }

    let mut exponent = BTreeMap::new();
    for (coordinate, &power) in degree.iter().enumerate() {
        if power == 0 {
            continue;
        }
        ckyz_series_add_scaled_assign(
            &mut exponent,
            &alpha[coordinate],
            Rational::from(Integer::from(power)),
        );
    }
    let exp_exponent = ckyz_series_exp_domain(&exponent, domain)?;
    let monomial = BTreeMap::from([(degree.to_vec(), Rational::from(1))]);
    ckyz_series_mul_domain(&monomial, &exp_exponent, domain)
}

#[cfg(test)]
fn ckyz_indexed_alpha_series(
    alpha: &[BTreeMap<Vec<usize>, Rational>],
    domain: &CkyzMonomialDomain,
) -> Result<Vec<CkyzIndexedSeries>> {
    if alpha.len() != domain.rank {
        return Err(Error::InvalidInput(
            "CKYZ indexed alpha rank does not match domain rank".into(),
        ));
    }
    alpha
        .iter()
        .map(|series| {
            validate_ckyz_series(series, domain.rank, "CKYZ log-period correction")?;
            validate_ckyz_series_has_zero_constant(
                series,
                domain.rank,
                "CKYZ log-period correction",
            )?;
            CkyzIndexedSeries::from_btree(series, domain, "CKYZ indexed log-period correction")
        })
        .collect()
}

#[cfg(test)]
fn ckyz_indexed_q_delta_series_in_z_domain(
    degree: &[usize],
    indexed_alpha: &[CkyzIndexedSeries],
    domain: &CkyzMonomialDomain,
) -> Result<CkyzIndexedSeries> {
    let rank = indexed_alpha.len();
    if rank == 0 || degree.len() != rank || domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ indexed q-delta series rank mismatch".into(),
        ));
    }
    if !domain.contains(degree) {
        return Ok(CkyzIndexedSeries { terms: Vec::new() });
    }

    let mut exponent = CkyzIndexedSeries { terms: Vec::new() };
    for (coordinate, &power) in degree.iter().enumerate() {
        if power == 0 {
            continue;
        }
        exponent.add_scaled_assign(
            &indexed_alpha[coordinate],
            Rational::from(Integer::from(power)),
        );
    }
    let exp_exponent = exponent.exp(domain)?;
    let monomial =
        CkyzIndexedSeries::monomial(degree, Rational::from(1), domain, "CKYZ indexed q-delta")?;
    monomial.mul(&exp_exponent, domain)
}

#[cfg(test)]
#[derive(Clone, Debug, PartialEq, Eq)]
struct CkyzIndexedQnBuild {
    series: CkyzIndexedSeries,
    reused_previous: bool,
    delta_degree: Vec<usize>,
}

#[cfg(test)]
fn ckyz_indexed_q_degree_series_with_previous_cache_in_z_domain(
    degree_index: usize,
    previous_qn: &VecDeque<HashMap<usize, CkyzIndexedSeries>>,
    previous_qn_indices: &VecDeque<Vec<usize>>,
    q_delta_cache: &mut HashMap<Vec<usize>, CkyzIndexedSeries>,
    indexed_alpha: &[CkyzIndexedSeries],
    domain: &CkyzMonomialDomain,
) -> Result<CkyzIndexedQnBuild> {
    let degree = domain
        .degrees
        .get(degree_index)
        .ok_or_else(|| Error::InvalidInput("CKYZ q_N degree index is outside the domain".into()))?;
    if degree.len() != domain.rank || indexed_alpha.len() != domain.rank {
        return Err(Error::InvalidInput(
            "CKYZ cached q_N series rank mismatch".into(),
        ));
    }

    let mut closest_previous_index = None;
    let mut closest_delta = degree.clone();
    let mut closest_distance = ckyz_qn_history_distance_score(degree);
    for (prev_indices, prev_qns) in previous_qn_indices.iter().zip(previous_qn.iter()) {
        for &prev_index in prev_indices {
            let previous_degree = domain.degrees.get(prev_index).ok_or_else(|| {
                Error::InvalidInput("CKYZ previous q_N index is outside the domain".into())
            })?;
            let Some(delta) = ckyz_subtract_degree_multiple(degree, previous_degree, 1) else {
                continue;
            };
            if !domain.contains(&delta) {
                continue;
            }
            let distance = ckyz_qn_history_distance_score(&delta);
            if distance < closest_distance {
                if !prev_qns.contains_key(&prev_index) {
                    return Err(Error::InvalidInput(
                        "CKYZ previous q_N index is missing its polynomial".into(),
                    ));
                }
                closest_distance = distance;
                closest_delta = delta;
                closest_previous_index = Some(prev_index);
            }
        }
    }

    if !q_delta_cache.contains_key(&closest_delta) {
        let q_delta =
            ckyz_indexed_q_delta_series_in_z_domain(&closest_delta, indexed_alpha, domain)?;
        q_delta_cache.insert(closest_delta.clone(), q_delta);
    }
    let q_delta = q_delta_cache
        .get(&closest_delta)
        .expect("q-delta series was inserted above");
    let series = if let Some(previous_index) = closest_previous_index {
        let previous_series = previous_qn
            .iter()
            .find_map(|batch| batch.get(&previous_index))
            .ok_or_else(|| Error::InvalidInput("CKYZ previous q_N polynomial vanished".into()))?;
        previous_series.mul(q_delta, domain)?
    } else {
        q_delta.clone()
    };

    Ok(CkyzIndexedQnBuild {
        series,
        reused_previous: closest_previous_index.is_some(),
        delta_degree: closest_delta,
    })
}

#[cfg(test)]
fn ckyz_expalpha_power_caches_domain(
    alpha: &[BTreeMap<Vec<usize>, Rational>],
    degrees: &[Vec<usize>],
    domain: &CkyzMonomialDomain,
) -> Result<Vec<Vec<BTreeMap<Vec<usize>, Rational>>>> {
    let rank = alpha.len();
    if rank == 0 || domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ exp(alpha) cache rank mismatch".into(),
        ));
    }
    for series in alpha {
        validate_ckyz_series(series, rank, "CKYZ log-period correction")?;
        validate_ckyz_series_has_zero_constant(series, rank, "CKYZ log-period correction")?;
    }

    let mut max_powers = vec![0usize; rank];
    for degree in degrees {
        if degree.len() != rank {
            return Err(Error::InvalidInput(
                "CKYZ exp(alpha) cache degree rank mismatch".into(),
            ));
        }
        for (coordinate, &entry) in degree.iter().enumerate() {
            max_powers[coordinate] = max_powers[coordinate].max(entry);
        }
    }

    alpha
        .iter()
        .zip(max_powers.iter())
        .map(|(series, &max_power)| {
            let exp_alpha = ckyz_series_exp_domain(series, domain)?;
            ckyz_series_power_cache_domain(&exp_alpha, max_power, domain)
        })
        .collect()
}

#[cfg(test)]
fn ckyz_q_degree_series_from_expalpha_powers_in_z_domain(
    degree: &[usize],
    expalpha_powers: &[Vec<BTreeMap<Vec<usize>, Rational>>],
    domain: &CkyzMonomialDomain,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    let rank = expalpha_powers.len();
    if rank == 0 || degree.len() != rank || domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ cached q-degree series rank mismatch".into(),
        ));
    }
    if !domain.contains(degree) {
        return Ok(BTreeMap::new());
    }

    let mut out = BTreeMap::from([(degree.to_vec(), Rational::from(1))]);
    for (coordinate, &power) in degree.iter().enumerate() {
        if power == 0 {
            continue;
        }
        let Some(power_series) = expalpha_powers
            .get(coordinate)
            .and_then(|coordinate_powers| coordinate_powers.get(power))
        else {
            return Err(Error::InvalidInput(
                "CKYZ cached q-degree requested an uncached exp(alpha) power".into(),
            ));
        };
        out = ckyz_series_mul_domain(&out, power_series, domain)?;
        if out.is_empty() {
            break;
        }
    }
    Ok(out)
}

fn ckyz_delta_degree(target: &[usize], base: &[usize]) -> Option<Vec<usize>> {
    target
        .iter()
        .zip(base.iter())
        .map(|(&target_entry, &base_entry)| target_entry.checked_sub(base_entry))
        .collect()
}

fn ckyz_q_degree_nonzero_coordinate_key(degree: &[usize]) -> Vec<usize> {
    degree
        .iter()
        .enumerate()
        .filter_map(|(coordinate, &power)| (power != 0).then_some(coordinate))
        .collect()
}

fn ckyz_q_degree_exp_support_for_coordinate_key(
    coordinate_key: &[usize],
    degree: &[usize],
    alpha_supports: &[BTreeSet<usize>],
    domain: &CkyzMonomialDomain,
) -> Result<BTreeSet<usize>> {
    let rank = alpha_supports.len();
    if rank == 0 || degree.len() != rank || domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ q-degree exponent-support rank mismatch".into(),
        ));
    }
    if coordinate_key.iter().any(|&coordinate| coordinate >= rank) {
        return Err(Error::InvalidInput(
            "CKYZ q-degree exponent-support coordinate is outside the relation rank".into(),
        ));
    }

    let mut exponent_support = BTreeSet::new();
    for &coordinate in coordinate_key {
        exponent_support.extend(alpha_supports[coordinate].iter().copied());
    }
    ckyz_support_exp_domain(&exponent_support, domain)
}

#[cfg(test)]
fn ckyz_series_li2_domain(
    series: &BTreeMap<Vec<usize>, Rational>,
    domain: &CkyzMonomialDomain,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    validate_ckyz_series_has_zero_constant(series, domain.rank, "CKYZ Li2 input")?;
    let max_exponent = ckyz_series_min_total_degree(series, domain.rank, "CKYZ Li2 input")?
        .map_or(0, |min_degree| domain.max_total_degree / min_degree);
    if max_exponent == 0 {
        return Ok(BTreeMap::new());
    }

    let mut out = series.clone();
    let mut power = series.clone();
    for exponent in 2..=max_exponent {
        power = ckyz_series_mul_domain(&power, series, domain)?;
        if power.is_empty() {
            break;
        }
        let exponent_squared = exponent
            .checked_mul(exponent)
            .ok_or_else(|| Error::InvalidInput("CKYZ Li2 exponent overflowed usize".into()))?;
        ckyz_series_add_scaled_assign(
            &mut out,
            &power,
            Rational::from_signeds(
                1i64,
                i64::try_from(exponent_squared).map_err(|_| {
                    Error::InvalidInput("CKYZ Li2 exponent does not fit i64".into())
                })?,
            ),
        );
    }
    Ok(out)
}

fn ckyz_componentwise_le(lhs: &[usize], rhs: &[usize]) -> bool {
    lhs.len() == rhs.len()
        && lhs
            .iter()
            .zip(rhs.iter())
            .all(|(&lhs_entry, &rhs_entry)| lhs_entry <= rhs_entry)
}

fn ckyz_subtract_degree_multiple(
    target: &[usize],
    degree: &[usize],
    multiple: usize,
) -> Option<Vec<usize>> {
    target
        .iter()
        .zip(degree.iter())
        .map(|(&target_entry, &degree_entry)| {
            degree_entry
                .checked_mul(multiple)
                .and_then(|scaled| target_entry.checked_sub(scaled))
        })
        .collect()
}

fn ckyz_subtract_degree_multiple_into(
    out: &mut Vec<usize>,
    target: &[usize],
    degree: &[usize],
    multiple: usize,
) -> Option<()> {
    if target.len() != degree.len() {
        return None;
    }
    out.clear();
    out.reserve(target.len());
    for (&target_entry, &degree_entry) in target.iter().zip(degree.iter()) {
        let scaled = degree_entry.checked_mul(multiple)?;
        out.push(target_entry.checked_sub(scaled)?);
    }
    Some(())
}

fn ckyz_q_degree_li2_support_intersects_indices_in_z_domain(
    degree: &[usize],
    alpha_supports: &[BTreeSet<usize>],
    target_indices: &[usize],
    exp_support_cache: &mut HashMap<Vec<usize>, BTreeSet<usize>>,
    domain: &CkyzMonomialDomain,
) -> Result<bool> {
    let coordinate_key = ckyz_q_degree_nonzero_coordinate_key(degree);
    if !exp_support_cache.contains_key(&coordinate_key) {
        let exp_support = ckyz_q_degree_exp_support_for_coordinate_key(
            &coordinate_key,
            degree,
            alpha_supports,
            domain,
        )?;
        exp_support_cache.insert(coordinate_key.clone(), exp_support);
    }
    let exp_support = exp_support_cache
        .get(&coordinate_key)
        .expect("exponential support was inserted above");

    for &target_index in target_indices {
        let target = &domain.degrees[target_index];
        if !ckyz_componentwise_le(degree, target) {
            continue;
        }
        let max_multiple = degree
            .iter()
            .zip(target.iter())
            .filter_map(|(&degree_entry, &target_entry)| {
                (degree_entry != 0).then(|| target_entry / degree_entry)
            })
            .min()
            .expect("candidate degree is nonzero");
        for multiple in 1..=max_multiple {
            let Some(delta) = ckyz_subtract_degree_multiple(target, degree, multiple) else {
                continue;
            };
            let Some(delta_index) = domain.index_of(&delta) else {
                continue;
            };
            if exp_support.contains(&delta_index) {
                return Ok(true);
            }
        }
    }

    Ok(false)
}

#[derive(Clone, Debug)]
struct CkyzScaledAlphaTerm {
    degree: Vec<usize>,
    total_degree: usize,
    coefficients: Vec<Rational>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct CkyzAlphaPredecessor {
    term_index: usize,
    remainder_index: usize,
    term_total_degree: usize,
}

#[derive(Clone, Debug, Default)]
struct CkyzExpCoefficientCache {
    scale_ids: HashMap<Vec<usize>, usize>,
    coefficients_by_scale: Vec<HashMap<usize, Rational>>,
    scaled_alpha_coefficients_by_scale: Vec<HashMap<usize, Rational>>,
    alpha_predecessors_by_delta: HashMap<usize, Vec<CkyzAlphaPredecessor>>,
}

impl CkyzExpCoefficientCache {
    fn scale_id(&mut self, scale_degree: &[usize]) -> usize {
        if let Some(&scale_id) = self.scale_ids.get(scale_degree) {
            return scale_id;
        }
        let scale_id = self.scale_ids.len();
        self.scale_ids.insert(scale_degree.to_vec(), scale_id);
        self.coefficients_by_scale.push(HashMap::new());
        self.scaled_alpha_coefficients_by_scale.push(HashMap::new());
        scale_id
    }

    fn coefficient(&self, scale_id: usize, delta_index: usize) -> Option<&Rational> {
        self.coefficients_by_scale
            .get(scale_id)
            .and_then(|coefficients| coefficients.get(&delta_index))
    }

    fn insert_coefficient(&mut self, scale_id: usize, delta_index: usize, coefficient: Rational) {
        let coefficients = self
            .coefficients_by_scale
            .get_mut(scale_id)
            .expect("scale id was allocated by scale_id");
        coefficients.insert(delta_index, coefficient);
    }

    fn scaled_alpha_coefficient(
        &mut self,
        scale_id: usize,
        scale_degree: &[usize],
        term_index: usize,
        term: &CkyzScaledAlphaTerm,
    ) -> Rational {
        let coefficients = self
            .scaled_alpha_coefficients_by_scale
            .get_mut(scale_id)
            .expect("scale id was allocated by scale_id");
        if let Some(coefficient) = coefficients.get(&term_index) {
            return coefficient.clone();
        }
        let coefficient = ckyz_scaled_alpha_coefficient(scale_degree, term);
        coefficients.insert(term_index, coefficient.clone());
        coefficient
    }

    fn alpha_predecessors(
        &mut self,
        delta_index: usize,
        alpha_terms: &[CkyzScaledAlphaTerm],
        domain: &CkyzMonomialDomain,
    ) -> Result<Vec<CkyzAlphaPredecessor>> {
        if let Some(predecessors) = self.alpha_predecessors_by_delta.get(&delta_index) {
            return Ok(predecessors.clone());
        }

        let delta = domain.degrees.get(delta_index).ok_or_else(|| {
            Error::InvalidInput("CKYZ alpha predecessor delta index is outside the domain".into())
        })?;
        let delta_total_degree = ckyz_total_degree(delta)?;
        let mut predecessors = Vec::new();
        for (term_index, term) in alpha_terms.iter().enumerate() {
            if term.total_degree > delta_total_degree {
                break;
            }
            if term.total_degree == 0 {
                continue;
            }
            if !ckyz_componentwise_le(&term.degree, delta) {
                continue;
            }
            let Some(remainder) = ckyz_delta_degree(delta, &term.degree) else {
                continue;
            };
            let Some(remainder_index) = domain.index_of(&remainder) else {
                continue;
            };
            predecessors.push(CkyzAlphaPredecessor {
                term_index,
                remainder_index,
                term_total_degree: term.total_degree,
            });
        }
        self.alpha_predecessors_by_delta
            .insert(delta_index, predecessors.clone());
        Ok(predecessors)
    }
}

fn ckyz_scaled_alpha_terms(
    alpha: &[BTreeMap<Vec<usize>, Rational>],
    domain: &CkyzMonomialDomain,
) -> Result<Vec<CkyzScaledAlphaTerm>> {
    let rank = alpha.len();
    if rank == 0 || domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ scaled-alpha term rank mismatch".into(),
        ));
    }

    let mut terms_by_degree = BTreeMap::<Vec<usize>, Vec<Rational>>::new();
    for (coordinate, series) in alpha.iter().enumerate() {
        validate_ckyz_series(series, rank, "CKYZ log-period correction")?;
        validate_ckyz_series_has_zero_constant(series, rank, "CKYZ log-period correction")?;
        for (degree, coefficient) in series {
            if *coefficient == 0 {
                continue;
            }
            if !domain.contains(degree) {
                continue;
            }
            let coefficients = terms_by_degree
                .entry(degree.clone())
                .or_insert_with(|| vec![Rational::from(0); rank]);
            coefficients[coordinate] = coefficient.clone();
        }
    }

    let mut terms = terms_by_degree
        .into_iter()
        .map(|(degree, coefficients)| {
            Ok(CkyzScaledAlphaTerm {
                total_degree: ckyz_total_degree(&degree)?,
                degree,
                coefficients,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    terms.sort_by(|lhs, rhs| {
        lhs.total_degree
            .cmp(&rhs.total_degree)
            .then_with(|| lhs.degree.cmp(&rhs.degree))
    });
    Ok(terms)
}

fn ckyz_scaled_alpha_coefficient(scale_degree: &[usize], term: &CkyzScaledAlphaTerm) -> Rational {
    let mut coefficient = Rational::from(0);
    for (&scale, term_coefficient) in scale_degree.iter().zip(term.coefficients.iter()) {
        if scale == 0 || *term_coefficient == 0 {
            continue;
        }
        coefficient += Rational::from(Integer::from(scale)) * term_coefficient.clone();
    }
    coefficient
}

fn ckyz_exp_scaled_alpha_coefficient_by_index_in_z_domain(
    scale_id: usize,
    scale_degree: &[usize],
    delta_index: usize,
    alpha_terms: &[CkyzScaledAlphaTerm],
    domain: &CkyzMonomialDomain,
    cache: &mut CkyzExpCoefficientCache,
) -> Result<Rational> {
    if scale_degree.len() != domain.rank {
        return Err(Error::InvalidInput(
            "CKYZ scaled exponential coefficient rank mismatch".into(),
        ));
    }
    let delta = domain.degrees.get(delta_index).ok_or_else(|| {
        Error::InvalidInput("CKYZ scaled exponential delta index is outside the domain".into())
    })?;
    if delta.iter().all(|&entry| entry == 0) {
        return Ok(Rational::from(1));
    }

    if let Some(coefficient) = cache.coefficient(scale_id, delta_index) {
        return Ok(coefficient.clone());
    }

    let delta_total_degree = ckyz_total_degree(delta)?;
    let mut coefficient = Rational::from(0);
    let predecessors = cache.alpha_predecessors(delta_index, alpha_terms, domain)?;
    for predecessor in predecessors {
        let term = &alpha_terms[predecessor.term_index];
        let scaled_alpha_coefficient =
            cache.scaled_alpha_coefficient(scale_id, scale_degree, predecessor.term_index, term);
        if scaled_alpha_coefficient == 0 {
            continue;
        }
        let remainder_coefficient = ckyz_exp_scaled_alpha_coefficient_by_index_in_z_domain(
            scale_id,
            scale_degree,
            predecessor.remainder_index,
            alpha_terms,
            domain,
            cache,
        )?;
        if remainder_coefficient == 0 {
            continue;
        }
        coefficient += Rational::from(Integer::from(predecessor.term_total_degree))
            * scaled_alpha_coefficient
            * remainder_coefficient;
    }
    coefficient /= Rational::from(Integer::from(delta_total_degree));
    cache.insert_coefficient(scale_id, delta_index, coefficient.clone());
    Ok(coefficient)
}

#[cfg(test)]
fn ckyz_q_degree_li2_coefficient_in_z_domain(
    degree: &[usize],
    target: &[usize],
    alpha_terms: &[CkyzScaledAlphaTerm],
    domain: &CkyzMonomialDomain,
    exp_support: Option<&BTreeSet<usize>>,
    exp_cache: &mut CkyzExpCoefficientCache,
) -> Result<Rational> {
    ckyz_q_degree_li2_coefficient_and_support_in_z_domain(
        degree,
        target,
        alpha_terms,
        domain,
        exp_support,
        exp_cache,
    )
    .map(|(coefficient, _)| coefficient)
}

fn ckyz_q_degree_li2_coefficient_and_support_in_z_domain(
    degree: &[usize],
    target: &[usize],
    alpha_terms: &[CkyzScaledAlphaTerm],
    domain: &CkyzMonomialDomain,
    exp_support: Option<&BTreeSet<usize>>,
    exp_cache: &mut CkyzExpCoefficientCache,
) -> Result<(Rational, bool)> {
    if degree.len() != domain.rank || target.len() != domain.rank {
        return Err(Error::InvalidInput(
            "CKYZ coefficient-level Li2 rank mismatch".into(),
        ));
    }
    if !ckyz_componentwise_le(degree, target) {
        return Ok((Rational::from(0), false));
    }

    let max_multiple = degree
        .iter()
        .zip(target.iter())
        .filter_map(|(&degree_entry, &target_entry)| {
            (degree_entry != 0).then(|| target_entry / degree_entry)
        })
        .min()
        .ok_or_else(|| Error::InvalidInput("CKYZ Li2 degree must be nonzero".into()))?;

    let mut coefficient = Rational::from(0);
    let mut has_supported_delta = false;
    let mut delta = Vec::with_capacity(domain.rank);
    for multiple in 1..=max_multiple {
        if ckyz_subtract_degree_multiple_into(&mut delta, target, degree, multiple).is_none() {
            continue;
        }
        if !domain.contains(&delta) {
            continue;
        }
        let delta_index = domain
            .index_of(&delta)
            .expect("delta was checked to be in domain");
        if exp_support.is_some_and(|support| !support.contains(&delta_index)) {
            continue;
        }
        has_supported_delta = true;
        let scale_degree = degree
            .iter()
            .map(|entry| {
                entry.checked_mul(multiple).ok_or_else(|| {
                    Error::InvalidInput("CKYZ scaled degree overflowed usize".into())
                })
            })
            .collect::<Result<Vec<_>>>()?;
        let scale_id = exp_cache.scale_id(&scale_degree);
        let exp_coefficient = ckyz_exp_scaled_alpha_coefficient_by_index_in_z_domain(
            scale_id,
            &scale_degree,
            delta_index,
            alpha_terms,
            domain,
            exp_cache,
        )?;
        if exp_coefficient == 0 {
            continue;
        }
        let multiple_squared = multiple
            .checked_mul(multiple)
            .ok_or_else(|| Error::InvalidInput("CKYZ Li2 exponent overflowed usize".into()))?;
        coefficient += exp_coefficient / Rational::from(Integer::from(multiple_squared));
    }
    Ok((coefficient, has_supported_delta))
}

fn ckyz_z_residual_dependency_degrees(
    alpha: &[BTreeMap<Vec<usize>, Rational>],
    terminal_degrees: &[Vec<usize>],
    domain: &CkyzMonomialDomain,
) -> Result<Vec<Vec<usize>>> {
    let rank = alpha.len();
    if rank == 0 || domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ z-residual dependency rank mismatch".into(),
        ));
    }
    validate_ckyz_target_degrees(terminal_degrees, rank)?;

    let trace_timing = env::var_os("CYRUS_TRACE_CKYZ_Z_HISTORY").is_some();
    let trace_start = trace_timing.then(std::time::Instant::now);
    let mut needed = BTreeSet::new();
    for degree in terminal_degrees {
        let Some(index) = domain.index_of(degree) else {
            return Err(Error::InvalidInput(format!(
                "CKYZ z-residual terminal degree {degree:?} is outside the monomial domain"
            )));
        };
        needed.insert(index);
    }

    let candidate_indices = domain
        .nonzero_degrees()
        .filter_map(|degree| domain.index_of(degree))
        .collect::<Vec<_>>();
    let mut exp_support_cache = HashMap::<Vec<usize>, BTreeSet<usize>>::new();
    let mut candidate_evaluations = 0usize;
    let alpha_supports = alpha
        .iter()
        .map(|series| ckyz_series_support_indices(series, domain))
        .collect::<Vec<_>>();
    let mut changed = true;
    while changed {
        changed = false;
        let needed_snapshot = needed.iter().copied().collect::<Vec<_>>();
        for &candidate_index in &candidate_indices {
            if needed.contains(&candidate_index) {
                continue;
            }
            let candidate_degree = &domain.degrees[candidate_index];
            if !needed_snapshot.iter().any(|&needed_index| {
                ckyz_componentwise_le(candidate_degree, &domain.degrees[needed_index])
            }) {
                continue;
            }

            candidate_evaluations += 1;
            if ckyz_q_degree_li2_support_intersects_indices_in_z_domain(
                candidate_degree,
                &alpha_supports,
                &needed_snapshot,
                &mut exp_support_cache,
                domain,
            )? {
                needed.insert(candidate_index);
                changed = true;
            }
        }
    }

    let mut out = needed
        .into_iter()
        .map(|index| domain.degrees[index].clone())
        .collect::<Vec<_>>();
    ckyz_sort_degrees_for_extraction(&mut out);
    if let Some(start) = trace_start {
        eprintln!(
            "[CKYZ_Z_HISTORY] domain={} terminals={} selected={} candidate_evaluations={} exp_support_classes={} elapsed={:?}",
            domain.degrees.len(),
            terminal_degrees.len(),
            out.len(),
            candidate_evaluations,
            exp_support_cache.len(),
            start.elapsed(),
        );
    }
    Ok(out)
}

fn extract_ckyz_local_gv_invariants_from_z_potential_for_degrees(
    potential_z_coefficients: &BTreeMap<Vec<usize>, Rational>,
    alpha: &[BTreeMap<Vec<usize>, Rational>],
    cover_weight_coefficients: &[i64],
    target_degrees: &[Vec<usize>],
    domain: &CkyzMonomialDomain,
) -> Result<BTreeMap<Vec<usize>, Integer>> {
    let rank = cover_weight_coefficients.len();
    if rank == 0 || alpha.len() != rank || domain.rank != rank {
        return Err(Error::InvalidInput(
            "CKYZ local GV z-series extraction rank mismatch".into(),
        ));
    }
    validate_ckyz_series(
        potential_z_coefficients,
        rank,
        "CKYZ local instanton potential",
    )?;
    validate_ckyz_target_degrees(target_degrees, rank)?;
    let grading_vector = ckyz_grading_vector_from_cover_weights(cover_weight_coefficients)?;

    let mut extraction_degrees = target_degrees.to_vec();
    ckyz_sort_degrees_for_extraction_with_grading(&mut extraction_degrees, &grading_vector)?;
    extraction_degrees.dedup();
    let extraction_gradings = extraction_degrees
        .iter()
        .map(|degree| ckyz_grading_degree(degree, &grading_vector))
        .collect::<Result<Vec<_>>>()?;
    let extraction_indices = extraction_degrees
        .iter()
        .map(|degree| {
            domain.index_of(degree).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "CKYZ z-series extraction degree {degree:?} is outside the monomial domain"
                ))
            })
        })
        .collect::<Result<Vec<_>>>()?;

    let mut residual_by_index = vec![None::<Rational>; domain.degrees.len()];
    for (&degree_index, degree) in extraction_indices.iter().zip(extraction_degrees.iter()) {
        residual_by_index[degree_index] = Some(
            potential_z_coefficients
                .get(degree)
                .cloned()
                .unwrap_or_else(|| Rational::from(0)),
        );
    }

    #[derive(Clone, Debug)]
    struct CkyzZResidualCandidate {
        position: usize,
        weight: i64,
        gv: Integer,
    }

    let mut invariants = BTreeMap::new();
    let trace_timing = env::var_os("CYRUS_TRACE_CKYZ_Z_HISTORY").is_some();
    let trace_progress = env::var_os("CYRUS_TRACE_CKYZ_Z_EXTRACT_PROGRESS").is_some();
    let extraction_start = trace_timing.then(std::time::Instant::now);
    let mut nonzero_gv_count = 0usize;
    let mut li2_coefficient_evaluations = 0usize;
    let mut li2_coefficient_probes = 0usize;
    let mut li2_supported_probes = 0usize;
    let alpha_terms = ckyz_scaled_alpha_terms(alpha, domain)?;
    let alpha_supports = alpha
        .iter()
        .map(|series| ckyz_series_support_indices(series, domain))
        .collect::<Vec<_>>();
    let mut exp_support_cache = HashMap::<Vec<usize>, BTreeSet<usize>>::new();
    let mut coefficient_cache = CkyzExpCoefficientCache::default();

    let mut batch_start = 0usize;
    while batch_start < extraction_degrees.len() {
        let batch_grading = extraction_gradings[batch_start];
        let batch_end = extraction_gradings.partition_point(|&grading| grading <= batch_grading);
        let mut candidates = Vec::new();
        for position in batch_start..batch_end {
            let degree = &extraction_degrees[position];
            let domain_degree_index = extraction_indices[position];
            let coefficient = residual_by_index[domain_degree_index]
                .as_ref()
                .cloned()
                .unwrap_or_else(|| Rational::from(0));
            if coefficient == 0 {
                continue;
            }
            let weight = ckyz_cover_weight(cover_weight_coefficients, degree)?;
            if weight == 0 {
                return Err(Error::InvalidInput(
                    "CKYZ local GV extraction encountered a nonzero coefficient with zero cover weight"
                        .into(),
                ));
            }
            let gv_rational = coefficient / Rational::from(weight);
            if gv_rational.denominator_ref() != &1u32 {
                return Err(Error::InvalidInput(format!(
                    "CKYZ local GV invariant at degree {degree:?} is not integral: {gv_rational}"
                )));
            }
            let gv = Integer::try_from(gv_rational).map_err(|_| {
                Error::InvalidInput(format!(
                    "CKYZ local GV invariant at degree {degree:?} is not integral"
                ))
            })?;
            if gv == 0 {
                continue;
            }
            invariants.insert(degree.clone(), gv.clone());
            nonzero_gv_count += 1;
            candidates.push(CkyzZResidualCandidate {
                position,
                weight,
                gv,
            });
        }

        let mut computed_qn_count = 0usize;
        for candidate in candidates {
            let subtraction_scale =
                -Rational::from(candidate.weight) * Rational::from(candidate.gv);
            let degree = &extraction_degrees[candidate.position];
            let coordinate_key = ckyz_q_degree_nonzero_coordinate_key(degree);
            if !exp_support_cache.contains_key(&coordinate_key) {
                let exp_support = ckyz_q_degree_exp_support_for_coordinate_key(
                    &coordinate_key,
                    degree,
                    &alpha_supports,
                    domain,
                )?;
                exp_support_cache.insert(coordinate_key.clone(), exp_support);
            }
            let exp_support = exp_support_cache
                .get(&coordinate_key)
                .expect("exponential support was inserted above");
            for (target_position, target) in extraction_degrees.iter().enumerate() {
                if target_position <= candidate.position
                    || extraction_gradings[target_position] <= batch_grading
                {
                    continue;
                };
                li2_coefficient_probes =
                    li2_coefficient_probes.checked_add(1).ok_or_else(|| {
                        Error::InvalidInput("CKYZ Li2 coefficient probe count overflowed".into())
                    })?;
                let (li2_coefficient, has_supported_delta) =
                    ckyz_q_degree_li2_coefficient_and_support_in_z_domain(
                        degree,
                        target,
                        &alpha_terms,
                        domain,
                        Some(exp_support),
                        &mut coefficient_cache,
                    )?;
                if !has_supported_delta {
                    continue;
                }
                li2_supported_probes = li2_supported_probes.checked_add(1).ok_or_else(|| {
                    Error::InvalidInput("CKYZ supported Li2 probe count overflowed".into())
                })?;
                if li2_coefficient == 0 {
                    continue;
                }
                let target_index = extraction_indices[target_position];
                li2_coefficient_evaluations =
                    li2_coefficient_evaluations.checked_add(1).ok_or_else(|| {
                        Error::InvalidInput("CKYZ Li2 coefficient count overflowed".into())
                    })?;
                let entry =
                    residual_by_index[target_index].get_or_insert_with(|| Rational::from(0));
                *entry += subtraction_scale.clone() * li2_coefficient;
            }
            computed_qn_count = computed_qn_count
                .checked_add(1)
                .ok_or_else(|| Error::InvalidInput("CKYZ computed q_N count overflowed".into()))?;
        }
        if trace_progress && computed_qn_count != 0 {
            eprintln!(
                "[CKYZ_Z_EXTRACT_PROGRESS] grading={} batch_qns={} total_qns={} li2_probes={} li2_supported_probes={} li2_updates={} exp_support_classes={} exp_coeff_scales={}",
                batch_grading,
                computed_qn_count,
                nonzero_gv_count,
                li2_coefficient_probes,
                li2_supported_probes,
                li2_coefficient_evaluations,
                exp_support_cache.len(),
                coefficient_cache.scale_ids.len(),
            );
        }
        batch_start = batch_end;
    }
    if let Some(start) = extraction_start {
        eprintln!(
            "[CKYZ_Z_EXTRACT] degrees={} nonzero_gvs={} li2_probes={} li2_supported_probes={} li2_updates={} exp_support_classes={} exp_coeff_scales={} elapsed={:?}",
            extraction_degrees.len(),
            nonzero_gv_count,
            li2_coefficient_probes,
            li2_supported_probes,
            li2_coefficient_evaluations,
            exp_support_cache.len(),
            coefficient_cache.scale_ids.len(),
            start.elapsed(),
        );
    }

    Ok(invariants)
}

fn extract_ckyz_local_gv_invariants_from_potential_for_degrees(
    potential_coefficients: &BTreeMap<Vec<usize>, Rational>,
    cover_weight_coefficients: &[i64],
    target_degrees: &[Vec<usize>],
) -> Result<BTreeMap<Vec<usize>, Integer>> {
    let rank = cover_weight_coefficients.len();
    if rank == 0 {
        return Err(Error::InvalidInput(
            "CKYZ local GV extraction requires at least one coordinate".into(),
        ));
    }
    validate_ckyz_series(
        potential_coefficients,
        rank,
        "CKYZ local instanton potential",
    )?;
    validate_ckyz_target_degrees(target_degrees, rank)?;
    let grading_vector = ckyz_grading_vector_from_cover_weights(cover_weight_coefficients)?;

    let mut extraction_degrees = target_degrees.to_vec();
    ckyz_sort_degrees_for_extraction_with_grading(&mut extraction_degrees, &grading_vector)?;
    extraction_degrees.dedup();

    let mut invariants = BTreeMap::new();
    for degree in &extraction_degrees {
        let mut residual = potential_coefficients
            .get(degree)
            .cloned()
            .unwrap_or_else(|| Rational::from(0));
        for cover in ckyz_nontrivial_covers(degree) {
            let primitive = ckyz_divide_degree(degree, cover)?;
            let Some(gv) = invariants.get(&primitive) else {
                continue;
            };
            let weight = ckyz_cover_weight(cover_weight_coefficients, &primitive)?;
            let cover_squared = cover.checked_mul(cover).ok_or_else(|| {
                Error::InvalidInput("CKYZ cover multiplicity overflowed usize".into())
            })?;
            residual -= Rational::from(weight) * Rational::from(gv)
                / Rational::from(Integer::from(cover_squared));
        }

        if let Some(gv) =
            ckyz_integral_gv_from_residual(residual, cover_weight_coefficients, degree)?
        {
            invariants.insert(degree.clone(), gv);
        }
    }

    Ok(invariants)
}

fn ckyz_integral_gv_from_residual(
    residual: Rational,
    cover_weight_coefficients: &[i64],
    degree: &[usize],
) -> Result<Option<Integer>> {
    if residual == 0 {
        return Ok(None);
    }
    let weight = ckyz_cover_weight(cover_weight_coefficients, degree)?;
    if weight == 0 {
        return Err(Error::InvalidInput(
            "CKYZ local GV extraction encountered a nonzero coefficient with zero cover weight"
                .into(),
        ));
    }
    let gv_rational = residual / Rational::from(weight);
    if gv_rational.denominator_ref() != &1u32 {
        return Err(Error::InvalidInput(format!(
            "CKYZ local GV invariant at degree {degree:?} is not integral: {gv_rational}"
        )));
    }
    let gv = Integer::try_from(gv_rational).map_err(|_| {
        Error::InvalidInput(format!(
            "CKYZ local GV invariant at degree {degree:?} is not integral"
        ))
    })?;
    Ok((gv != 0).then_some(gv))
}

fn integer_row_transform_between_bases(
    target_rows: &[Vec<i64>],
    source_rows: &[Vec<i64>],
) -> Result<Option<Vec<Vec<i64>>>> {
    let mut transform = Vec::with_capacity(target_rows.len());
    for target_row in target_rows {
        let Ok(coordinates) = relation_coordinates_in_row_basis(target_row, source_rows) else {
            return Ok(None);
        };
        transform.push(coordinates);
    }
    let transformed = transform_rows_i64(&transform, source_rows);
    if transformed != target_rows {
        return Ok(None);
    }
    if !is_unimodular_i64_matrix(&transform) {
        return Ok(None);
    }
    Ok(Some(transform))
}

fn transform_rows_i64(transform: &[Vec<i64>], matrix: &[Vec<i64>]) -> Vec<Vec<i64>> {
    transform
        .iter()
        .map(|row| {
            (0..matrix[0].len())
                .map(|col| {
                    row.iter()
                        .zip(matrix.iter())
                        .map(|(coefficient, source_row)| coefficient * source_row[col])
                        .sum()
                })
                .collect()
        })
        .collect()
}

fn is_unimodular_i64_matrix(matrix: &[Vec<i64>]) -> bool {
    if matrix.is_empty() || matrix.iter().any(|row| row.len() != matrix.len()) {
        return false;
    }
    let mut rational_matrix = matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|&value| Rational::from(Integer::from(value)))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let determinant = crate::integer_math::determinant_gaussian(&mut rational_matrix);
    determinant == Rational::from(1) || determinant == Rational::from(-1)
}

fn permutations(size: usize) -> Vec<Vec<usize>> {
    let mut values = (0..size).collect::<Vec<_>>();
    let mut out = Vec::new();
    push_permutations(0, &mut values, &mut out);
    out
}

fn push_permutations(start: usize, values: &mut [usize], out: &mut Vec<Vec<usize>>) {
    if start == values.len() {
        out.push(values.to_vec());
        return;
    }
    for idx in start..values.len() {
        values.swap(start, idx);
        push_permutations(start + 1, values, out);
        values.swap(start, idx);
    }
}

fn permute_matrix_columns(matrix: &[Vec<i64>], permutation: &[usize]) -> Vec<Vec<i64>> {
    matrix
        .iter()
        .map(|row| permutation.iter().map(|&idx| row[idx]).collect())
        .collect()
}

fn permute_vector_values(vector: &[i64], permutation: &[usize]) -> Vec<i64> {
    permutation.iter().map(|&idx| vector[idx]).collect()
}

fn relation_coordinates_in_row_basis(target: &[i64], rows: &[Vec<i64>]) -> Result<Vec<i64>> {
    if rows.is_empty() {
        if target.iter().all(|&value| value == 0) {
            return Ok(Vec::new());
        }
        return Err(Error::InvalidInput(
            "nonzero target relation cannot be expressed in an empty row basis".into(),
        ));
    }

    let dim = rows[0].len();
    if target.len() != dim {
        return Err(Error::InvalidInput(
            "target relation dimension does not match row basis".into(),
        ));
    }
    for row in rows {
        if row.len() != dim {
            return Err(Error::InvalidInput(
                "row basis has inconsistent dimensions".into(),
            ));
        }
    }

    for columns in column_combinations(dim, rows.len()) {
        let matrix = columns
            .iter()
            .map(|&col| {
                rows.iter()
                    .map(|row| Rational::from(Integer::from(row[col])))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let rhs = columns
            .iter()
            .map(|&col| Rational::from(Integer::from(target[col])))
            .collect::<Vec<_>>();
        let Some(coordinates) = solve_linear_system_rational(&matrix, &rhs) else {
            continue;
        };
        if !row_basis_coordinates_match_target(&coordinates, target, rows) {
            continue;
        }
        return coordinates
            .into_iter()
            .map(|coordinate| {
                if coordinate.denominator_ref() != &1u32 {
                    return Err(Error::InvalidInput(
                        "target relation has non-integral row-basis coordinates".into(),
                    ));
                }
                let integer = Integer::try_from(coordinate).map_err(|_| {
                    Error::InvalidInput(
                        "target relation row-basis coordinate is not integral".into(),
                    )
                })?;
                i64::try_from(&integer).map_err(|_| {
                    Error::InvalidInput(
                        "target relation row-basis coordinate does not fit in i64".into(),
                    )
                })
            })
            .collect();
    }

    Err(Error::InvalidInput(
        "target relation is not in the supplied row basis".into(),
    ))
}

fn row_basis_coordinates_match_target(
    coordinates: &[Rational],
    target: &[i64],
    rows: &[Vec<i64>],
) -> bool {
    for col in 0..target.len() {
        let mut reconstructed = Rational::from(0);
        for (coordinate, row) in coordinates.iter().zip(rows.iter()) {
            reconstructed += coordinate * Rational::from(Integer::from(row[col]));
        }
        if reconstructed != Rational::from(Integer::from(target[col])) {
            return false;
        }
    }
    true
}

fn column_combinations(n_columns: usize, size: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut current = Vec::with_capacity(size);
    push_column_combinations(0, n_columns, size, &mut current, &mut out);
    out
}

fn push_column_combinations(
    start: usize,
    n_columns: usize,
    size: usize,
    current: &mut Vec<usize>,
    out: &mut Vec<Vec<usize>>,
) {
    if current.len() == size {
        out.push(current.clone());
        return;
    }
    let remaining = size - current.len();
    if remaining > n_columns.saturating_sub(start) {
        return;
    }
    for col in start..=n_columns - remaining {
        current.push(col);
        push_column_combinations(col + 1, n_columns, size, current, out);
        current.pop();
    }
}

fn canonical_rank_two_signature_entries(
    entries: &[RankTwoLocalSupportSignatureEntry],
) -> Option<Vec<RankTwoLocalSupportSignatureEntry>> {
    if entries.is_empty() {
        return Some(Vec::new());
    }

    let mut candidates = Vec::new();
    push_rank_two_signature_candidate(&mut candidates, entries.to_vec());

    for anchor in entries {
        for first in entries {
            let first_basis = coordinate_2d_difference(&first.coordinates, &anchor.coordinates);
            if first_basis == [0, 0] {
                continue;
            }
            for second in entries {
                let second_basis =
                    coordinate_2d_difference(&second.coordinates, &anchor.coordinates);
                if second_basis == [0, 0] {
                    continue;
                }
                let det = first_basis[0] * second_basis[1] - first_basis[1] * second_basis[0];
                if det == 0 {
                    continue;
                }

                let mut transformed = Vec::with_capacity(entries.len());
                let first_basis = first_basis.to_vec();
                let second_basis = second_basis.to_vec();
                let mut valid = true;
                for entry in entries {
                    let target =
                        coordinate_2d_difference(&entry.coordinates, &anchor.coordinates).to_vec();
                    let Some(coordinates) =
                        solve_in_two_vector_basis(&first_basis, &second_basis, &target)
                    else {
                        valid = false;
                        break;
                    };
                    transformed.push(RankTwoLocalSupportSignatureEntry {
                        coefficient: entry.coefficient,
                        coordinates,
                    });
                }
                if valid {
                    push_rank_two_signature_candidate(&mut candidates, transformed);
                }
            }
        }
    }

    candidates.into_iter().min()
}

fn push_rank_two_signature_candidate(
    candidates: &mut Vec<Vec<RankTwoLocalSupportSignatureEntry>>,
    entries: Vec<RankTwoLocalSupportSignatureEntry>,
) {
    candidates.push(normalize_rank_two_signature_entries(entries.clone(), 1));
    candidates.push(normalize_rank_two_signature_entries(entries, -1));
}

fn coordinate_2d_difference(lhs: &[i64; 2], rhs: &[i64; 2]) -> [i64; 2] {
    [lhs[0] - rhs[0], lhs[1] - rhs[1]]
}

fn normalize_rank_two_signature_entries(
    mut entries: Vec<RankTwoLocalSupportSignatureEntry>,
    sign: i64,
) -> Vec<RankTwoLocalSupportSignatureEntry> {
    let origin = entries
        .iter()
        .map(|entry| entry.coordinates)
        .min()
        .unwrap_or([0, 0]);
    for entry in &mut entries {
        entry.coefficient *= sign;
        entry.coordinates[0] -= origin[0];
        entry.coordinates[1] -= origin[1];
    }
    entries.sort_unstable();
    entries
}

fn local_affine_charge_basis(
    relation_points: &[AffineCircuitRelationPoint],
) -> Result<Vec<Vec<i64>>> {
    let Some(first_point) = relation_points.first() else {
        return Ok(Vec::new());
    };
    let dim = first_point.coordinates.len();
    let support_len = relation_points.len();
    let mut matrix = vec![vec![Integer::from(0); support_len]; dim + 1];
    for (col, point) in relation_points.iter().enumerate() {
        matrix[0][col] = Integer::from(1);
        for (row, &coordinate) in point.coordinates.iter().enumerate() {
            matrix[row + 1][col] = Integer::from(coordinate);
        }
    }

    integer_kernel(&matrix)
        .into_iter()
        .map(|row| {
            let mut converted = row
                .iter()
                .map(|value| {
                    i64::try_from(value).map_err(|_| {
                        Error::InvalidInput(
                            "local affine charge basis entry does not fit in i64".into(),
                        )
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            normalize_relation_orientation(&mut converted);
            Ok(converted)
        })
        .collect()
}

fn normalize_relation_orientation(row: &mut [i64]) {
    let Some(first_nonzero) = row.iter().find(|&&value| value != 0).copied() else {
        return;
    };
    if first_nonzero < 0 {
        for value in row {
            *value = -*value;
        }
    }
}

fn local_affine_coordinates(
    relation_points: &[AffineCircuitRelationPoint],
    affine_rank: usize,
) -> Result<Vec<LocalToricCoordinate>> {
    let Some(base) = relation_points.first() else {
        return Ok(Vec::new());
    };
    if affine_rank == 0 {
        let mut local_coordinates = relation_points
            .iter()
            .map(|point| LocalToricCoordinate {
                point_index: point.point_index,
                coordinates: Vec::new(),
            })
            .collect::<Vec<_>>();
        local_coordinates.sort_by_key(|point| point.point_index);
        return Ok(local_coordinates);
    }

    let difference_rows = relation_points
        .iter()
        .skip(1)
        .map(|point| {
            coordinate_difference(&point.coordinates, &base.coordinates)
                .into_iter()
                .map(Integer::from)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let hnf = hermite_normal_form(&difference_rows);
    let basis = hnf
        .into_iter()
        .filter(|row| row.iter().any(|value| *value != 0))
        .map(|row| {
            row.into_iter()
                .map(|value| {
                    i64::try_from(&value).map_err(|_| {
                        Error::InvalidInput(
                            "local affine coordinate basis entry does not fit in i64".into(),
                        )
                    })
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;
    if basis.len() != affine_rank {
        return Err(Error::InvalidInput(format!(
            "rank-{affine_rank} local support produced {} lattice basis rows",
            basis.len()
        )));
    }

    let mut local_coordinates = Vec::with_capacity(relation_points.len());
    for point in relation_points {
        let target = coordinate_difference(&point.coordinates, &base.coordinates);
        let coordinates = relation_coordinates_in_row_basis(&target, &basis)?;
        if coordinates.len() != affine_rank {
            return Err(Error::InvalidInput(
                "local affine coordinate dimension does not match affine rank".into(),
            ));
        }
        local_coordinates.push(LocalToricCoordinate {
            point_index: point.point_index,
            coordinates,
        });
    }
    local_coordinates.sort_by_key(|point| point.point_index);
    Ok(local_coordinates)
}

fn local_rank_two_coordinates(
    relation_points: &[AffineCircuitRelationPoint],
    affine_rank: usize,
) -> Result<Option<Vec<LocalToricCoordinate2D>>> {
    if affine_rank != 2 {
        return Ok(None);
    }
    let Some(base) = relation_points.first() else {
        return Ok(Some(Vec::new()));
    };

    let difference_rows: Vec<Vec<Integer>> = relation_points
        .iter()
        .skip(1)
        .map(|point| {
            coordinate_difference(&point.coordinates, &base.coordinates)
                .into_iter()
                .map(Integer::from)
                .collect()
        })
        .collect();
    let hnf = hermite_normal_form(&difference_rows);
    let basis_rows = hnf
        .into_iter()
        .filter(|row| row.iter().any(|value| *value != 0))
        .collect::<Vec<_>>();
    if basis_rows.len() != 2 {
        return Err(Error::InvalidInput(format!(
            "rank-two local support produced {} lattice basis rows",
            basis_rows.len()
        )));
    }

    let basis = basis_rows
        .iter()
        .map(|row| {
            row.iter()
                .map(|value| {
                    i64::try_from(value).map_err(|_| {
                        Error::InvalidInput(
                            "rank-two local coordinate basis entry does not fit in i64".into(),
                        )
                    })
                })
                .collect::<Result<Vec<_>>>()
        })
        .collect::<Result<Vec<_>>>()?;

    let mut local_coordinates = Vec::with_capacity(relation_points.len());
    for point in relation_points {
        let target = coordinate_difference(&point.coordinates, &base.coordinates);
        let coordinates =
            solve_in_two_vector_basis(&basis[0], &basis[1], &target).ok_or_else(|| {
                Error::InvalidInput(
                    "rank-two local support point is not integral in reconstructed lattice basis"
                        .into(),
                )
            })?;
        local_coordinates.push(LocalToricCoordinate2D {
            point_index: point.point_index,
            coordinates,
        });
    }
    local_coordinates.sort_by_key(|point| point.point_index);
    Ok(Some(local_coordinates))
}

fn affine_relation_rank(relation_points: &[AffineCircuitRelationPoint]) -> usize {
    let Some(base) = relation_points.first() else {
        return 0;
    };
    let rows: Vec<Vec<i64>> = relation_points
        .iter()
        .skip(1)
        .map(|point| {
            point
                .coordinates
                .iter()
                .zip(base.coordinates.iter())
                .map(|(&coord, &base_coord)| coord - base_coord)
                .collect()
        })
        .collect();
    integer_matrix_rank(&rows)
}

fn classify_local_toric_circuit(
    relation_points: &[AffineCircuitRelationPoint],
    affine_rank: usize,
) -> Option<LocalToricCircuitKind> {
    if affine_rank != 2 {
        return None;
    }

    for vertex_coefficient in [1, -1] {
        let interior_coefficient = -3 * vertex_coefficient;
        let mut vertices: Vec<&AffineCircuitRelationPoint> = relation_points
            .iter()
            .filter(|point| point.coefficient == vertex_coefficient)
            .collect();
        let interior: Vec<&AffineCircuitRelationPoint> = relation_points
            .iter()
            .filter(|point| point.coefficient == interior_coefficient)
            .collect();
        if vertices.len() != 3 || interior.len() != 1 {
            continue;
        }
        if vertices.len() + interior.len() != relation_points.len() {
            continue;
        }

        let interior = interior[0];
        let dim = interior.coordinates.len();
        if vertices
            .iter()
            .any(|vertex| vertex.coordinates.len() != dim)
        {
            continue;
        }
        let is_barycenter = (0..dim).all(|coord_idx| {
            let vertex_sum: i128 = vertices
                .iter()
                .map(|vertex| i128::from(vertex.coordinates[coord_idx]))
                .sum();
            3 * i128::from(interior.coordinates[coord_idx]) == vertex_sum
        });
        if !is_barycenter {
            continue;
        }

        vertices.sort_by_key(|point| point.point_index);
        let local_coordinates = local_p2_triangle_coordinates(interior, &vertices)?;
        return Some(LocalToricCircuitKind::LocalP2Triangle {
            interior_point: interior.point_index,
            vertex_points: vertices
                .into_iter()
                .map(|point| point.point_index)
                .collect(),
            interior_coefficient,
            vertex_coefficient,
            local_coordinates,
        });
    }

    None
}

fn local_p2_triangle_coordinates(
    interior: &AffineCircuitRelationPoint,
    vertices: &[&AffineCircuitRelationPoint],
) -> Option<Vec<LocalToricCoordinate2D>> {
    if vertices.len() != 3 {
        return None;
    }
    let first_basis = coordinate_difference(&vertices[0].coordinates, &interior.coordinates);
    let second_basis = coordinate_difference(&vertices[1].coordinates, &interior.coordinates);
    let mut local_coordinates = vec![LocalToricCoordinate2D {
        point_index: interior.point_index,
        coordinates: [0, 0],
    }];
    for vertex in vertices {
        let target = coordinate_difference(&vertex.coordinates, &interior.coordinates);
        let coordinates = solve_in_two_vector_basis(&first_basis, &second_basis, &target)?;
        local_coordinates.push(LocalToricCoordinate2D {
            point_index: vertex.point_index,
            coordinates,
        });
    }
    local_coordinates.sort_by_key(|point| point.point_index);
    Some(local_coordinates)
}

fn coordinate_difference(lhs: &[i64], rhs: &[i64]) -> Vec<i64> {
    lhs.iter()
        .zip(rhs.iter())
        .map(|(&left, &right)| left - right)
        .collect()
}

fn solve_in_two_vector_basis(
    first_basis: &[i64],
    second_basis: &[i64],
    target: &[i64],
) -> Option<[i64; 2]> {
    for first_coord in 0..target.len() {
        for second_coord in (first_coord + 1)..target.len() {
            let det = i128::from(first_basis[first_coord]) * i128::from(second_basis[second_coord])
                - i128::from(first_basis[second_coord]) * i128::from(second_basis[first_coord]);
            if det == 0 {
                continue;
            }
            let first_num = i128::from(target[first_coord])
                * i128::from(second_basis[second_coord])
                - i128::from(target[second_coord]) * i128::from(second_basis[first_coord]);
            let second_num = i128::from(first_basis[first_coord])
                * i128::from(target[second_coord])
                - i128::from(first_basis[second_coord]) * i128::from(target[first_coord]);
            if first_num % det != 0 || second_num % det != 0 {
                continue;
            }
            let first = i64::try_from(first_num / det).ok()?;
            let second = i64::try_from(second_num / det).ok()?;
            let reconstructs_target = target
                .iter()
                .zip(first_basis.iter().zip(second_basis.iter()))
                .all(
                    |(&target_coord, (&first_basis_coord, &second_basis_coord))| {
                        i128::from(first) * i128::from(first_basis_coord)
                            + i128::from(second) * i128::from(second_basis_coord)
                            == i128::from(target_coord)
                    },
                );
            if reconstructs_target {
                return Some([first, second]);
            }
        }
    }
    None
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

/// Compute the exact rational rank of the span of curve rows.
///
/// This is used for potent-ray convergence audits, where the paper reports the
/// rank of the cone spanned by sampled low-dimensional-face rays.
pub fn curve_row_span_rank(curves: &[Vec<i64>]) -> Result<usize> {
    if curves.is_empty() {
        return Ok(0);
    }
    let dim = curves[0].len();
    for curve in curves {
        if curve.len() != dim {
            return Err(Error::InvalidInput(
                "curve rows have inconsistent dimensions".into(),
            ));
        }
    }
    Ok(integer_matrix_rank(curves))
}

/// Check whether a curve row lies in the rational span of other curve rows.
///
/// This is an exact rank comparison over `Q`, not a floating least-squares
/// test. The zero vector lies in the span of an empty row set; any nonzero
/// vector does not.
pub fn curve_in_rational_row_span(curve: &[i64], rows: &[Vec<i64>]) -> Result<bool> {
    if rows.is_empty() {
        return Ok(curve.iter().all(|&value| value == 0));
    }
    let dim = rows[0].len();
    if curve.len() != dim {
        return Err(Error::InvalidInput(format!(
            "curve dimension {} does not match row span dimension {dim}",
            curve.len()
        )));
    }
    for row in rows {
        if row.len() != dim {
            return Err(Error::InvalidInput(
                "curve rows have inconsistent dimensions".into(),
            ));
        }
    }

    let rank = integer_matrix_rank(rows);
    let mut augmented = rows.to_vec();
    augmented.push(curve.to_vec());
    Ok(integer_matrix_rank(&augmented) == rank)
}

fn checked_i128_dot(lhs: &[i64], rhs: &[i64], context: &str) -> Result<i128> {
    lhs.iter()
        .zip(rhs.iter())
        .try_fold(0i128, |acc, (&left, &right)| {
            let term = i128::from(left)
                .checked_mul(i128::from(right))
                .ok_or_else(|| Error::InvalidInput(format!("{context} dot product overflow")))?;
            acc.checked_add(term)
                .ok_or_else(|| Error::InvalidInput(format!("{context} dot product overflow")))
        })
}

fn i128_vector_to_i64(vector: &[i128], context: &str) -> Result<Vec<i64>> {
    vector
        .iter()
        .map(|&value| {
            i64::try_from(value)
                .map_err(|_| Error::InvalidInput(format!("{context} does not fit in i64")))
        })
        .collect()
}

fn checked_ray_multiple(ray: &[i64], multiple: i128, context: &str) -> Result<Vec<i64>> {
    ray.iter()
        .map(|&entry| {
            let value = i128::from(entry)
                .checked_mul(multiple)
                .ok_or_else(|| Error::InvalidInput(format!("{context} multiple overflow")))?;
            i64::try_from(value)
                .map_err(|_| Error::InvalidInput(format!("{context} multiple does not fit in i64")))
        })
        .collect()
}

fn checked_i64_vector_difference(lhs: &[i64], rhs: &[i64], context: &str) -> Result<Vec<i64>> {
    if lhs.len() != rhs.len() {
        return Err(Error::InvalidInput(format!(
            "{context} vector dimensions do not match"
        )));
    }
    lhs.iter()
        .zip(rhs.iter())
        .map(|(&left, &right)| {
            left.checked_sub(right)
                .ok_or_else(|| Error::InvalidInput(format!("{context} subtraction overflow")))
        })
        .collect()
}

fn validate_nilpotent_slice_against_grading(
    slice: &NilpotentRayDegreeSlice,
    grading_vector: &[i64],
) -> Result<()> {
    if grading_vector.is_empty() {
        return Err(Error::InvalidInput(
            "nilpotent-ray slice distance requires a nonempty grading vector".into(),
        ));
    }
    if slice.primitive_ray.len() != grading_vector.len()
        || slice.slice_origin.len() != grading_vector.len()
    {
        return Err(Error::InvalidInput(
            "nilpotent-ray slice dimension does not match grading dimension".into(),
        ));
    }
    if slice.slice_degree <= 0 {
        return Err(Error::InvalidInput(
            "nilpotent-ray slice degree must be positive".into(),
        ));
    }
    let origin_degree = checked_i128_dot(
        &slice.slice_origin,
        grading_vector,
        "nilpotent-ray slice-origin grading",
    )?;
    if origin_degree != slice.slice_degree {
        return Err(Error::InvalidInput(format!(
            "nilpotent-ray slice origin has degree {origin_degree}, expected {}",
            slice.slice_degree
        )));
    }
    Ok(())
}

fn apply_i64_matrix_to_vector(
    matrix: &[Vec<i64>],
    vector: &[i64],
    context: &str,
) -> Result<Vec<i64>> {
    let mut output = Vec::with_capacity(matrix.len());
    for row in matrix {
        if row.len() != vector.len() {
            return Err(Error::InvalidInput(format!(
                "{context} matrix/vector dimensions do not match"
            )));
        }
        let value = row.iter().zip(vector.iter()).try_fold(
            0i128,
            |acc, (&matrix_entry, &vector_entry)| {
                let term = i128::from(matrix_entry)
                    .checked_mul(i128::from(vector_entry))
                    .ok_or_else(|| Error::InvalidInput(format!("{context} overflow")))?;
                acc.checked_add(term)
                    .ok_or_else(|| Error::InvalidInput(format!("{context} overflow")))
            },
        )?;
        output.push(
            i64::try_from(value).map_err(|_| {
                Error::InvalidInput(format!("{context} result does not fit in i64"))
            })?,
        );
    }
    Ok(output)
}

fn checked_i64_infinity_norm(vector: &[i64], context: &str) -> Result<i64> {
    vector.iter().try_fold(0i64, |current, &entry| {
        let abs = entry
            .checked_abs()
            .ok_or_else(|| Error::InvalidInput(format!("{context} contains i64::MIN")))?;
        Ok(current.max(abs))
    })
}

fn collect_nilpotent_divergence_checks_by_ray<'a>(
    checks: &'a [(Vec<i64>, NilpotentRayDivergenceCheck)],
    context: &str,
) -> Result<BTreeMap<Vec<i64>, &'a NilpotentRayDivergenceCheck>> {
    let mut by_ray = BTreeMap::new();
    for (ray, check) in checks {
        if ray != &check.half_cutoff.slice.primitive_ray
            || ray != &check.full_cutoff.slice.primitive_ray
        {
            return Err(Error::InvalidInput(format!(
                "two-pass nop classification {context} check ray does not match its slices"
            )));
        }
        if by_ray.insert(ray.clone(), check).is_some() {
            return Err(Error::InvalidInput(format!(
                "two-pass nop classification contains duplicate {context} checks"
            )));
        }
    }
    Ok(by_ray)
}

fn finite_gv_table_divergence_check_for_candidate(
    primitive_ray: &[i64],
    grading_vector: &[i64],
    cutoff_degree: i128,
    gv_invariants: &[(Vec<i64>, Integer)],
    comparison_charges: &[Vec<i64>],
    context: &str,
) -> Result<NilpotentRayDivergenceCheck> {
    let half_slice = nilpotent_ray_degree_slice_for_cutoff_fraction(
        primitive_ray,
        grading_vector,
        cutoff_degree,
        1,
        2,
    )?
    .ok_or_else(|| {
        Error::InvalidInput(format!(
            "finite GV nop {context} candidate ray does not reach the half-cutoff slice"
        ))
    })?;
    let full_slice = nilpotent_ray_degree_slice_for_cutoff_fraction(
        primitive_ray,
        grading_vector,
        cutoff_degree,
        1,
        1,
    )?
    .ok_or_else(|| {
        Error::InvalidInput(format!(
            "finite GV nop {context} candidate ray does not reach the full-cutoff slice"
        ))
    })?;

    let mut half_slice_lattice_points = finite_gv_nonzero_degree_slice_points(
        grading_vector,
        half_slice.slice_degree,
        gv_invariants,
    )?;
    half_slice_lattice_points.push(half_slice.slice_origin.clone());
    let mut full_slice_lattice_points = finite_gv_nonzero_degree_slice_points(
        grading_vector,
        full_slice.slice_degree,
        gv_invariants,
    )?;
    full_slice_lattice_points.push(full_slice.slice_origin.clone());

    let half_distance = nilpotent_ray_lll_reduced_slice_distance(
        &half_slice,
        &half_slice_lattice_points,
        comparison_charges,
        grading_vector,
    )?;
    let full_distance = nilpotent_ray_lll_reduced_slice_distance(
        &full_slice,
        &full_slice_lattice_points,
        comparison_charges,
        grading_vector,
    )?;
    nilpotent_ray_divergence_check_from_slice_distances(half_distance, full_distance)
}

fn primitive_i64_gcd(ray: &[i64], context: &str) -> Result<i64> {
    let mut gcd = 0i64;
    for &entry in ray {
        let abs = entry
            .checked_abs()
            .ok_or_else(|| Error::InvalidInput(format!("{context} contains i64::MIN")))?;
        gcd = gcd_i64(gcd, abs);
    }
    Ok(gcd)
}

fn primitive_ray_from_curve_class(class: &[i64], context: &str) -> Result<Vec<i64>> {
    let gcd = primitive_i64_gcd(class, context)?;
    if gcd == 0 {
        return Err(Error::InvalidInput(format!("{context} must not be zero")));
    }
    Ok(class.iter().map(|&entry| entry / gcd).collect())
}

/// Apply the finite-cutoff nilpotent-ray test from the McAllister
/// moduli-space reconstruction appendix.
///
/// For a co-prime curve class `C`, the source criterion adds `C` to the
/// apparently nilpotent set if some multiple `k* C` lies below the grading
/// cutoff, has zero genus-zero GV invariant, and the lower weighted sum
/// `sum_{k=1}^{k*-1} k^2 n^0_{kC}` is positive. Missing classes in the supplied
/// finite GV table are treated as zero at the cutoff already computed by the
/// caller.
pub fn detect_apparent_nilpotent_ray_from_gv_multiples(
    primitive_ray: &[i64],
    grading_vector: &[i64],
    cutoff_degree: i128,
    gv_invariants: &[(Vec<i64>, Integer)],
) -> Result<Option<NilpotentRayCandidate>> {
    if primitive_ray.is_empty() {
        return Err(Error::InvalidInput(
            "nilpotent-ray test requires a nonempty primitive ray".into(),
        ));
    }
    if grading_vector.len() != primitive_ray.len() {
        return Err(Error::InvalidInput(
            "nilpotent-ray grading dimension does not match ray dimension".into(),
        ));
    }
    if cutoff_degree <= 0 {
        return Err(Error::InvalidInput(
            "nilpotent-ray cutoff degree must be positive".into(),
        ));
    }

    let gcd = primitive_i64_gcd(primitive_ray, "nilpotent-ray primitive ray")?;
    if gcd == 0 {
        return Err(Error::InvalidInput(
            "nilpotent-ray primitive ray must not be zero".into(),
        ));
    }
    if gcd != 1 {
        return Err(Error::InvalidInput(format!(
            "nilpotent-ray primitive ray must be co-prime, got gcd {gcd}"
        )));
    }

    let degree = checked_i128_dot(primitive_ray, grading_vector, "nilpotent-ray grading")?;
    if degree <= 0 {
        return Err(Error::InvalidInput(format!(
            "nilpotent-ray primitive ray has non-positive grading degree {degree}"
        )));
    }

    let mut gv_by_class: HashMap<Vec<i64>, Integer> = HashMap::with_capacity(gv_invariants.len());
    for (class, gv) in gv_invariants {
        if class.len() != primitive_ray.len() {
            return Err(Error::InvalidInput(
                "nilpotent-ray GV invariant dimension does not match ray dimension".into(),
            ));
        }
        if let Some(existing) = gv_by_class.get(class) {
            if existing != gv {
                return Err(Error::InvalidInput(
                    "nilpotent-ray GV table contains conflicting duplicate classes".into(),
                ));
            }
        } else {
            gv_by_class.insert(class.clone(), gv.clone());
        }
    }

    let max_multiple = cutoff_degree / degree;
    if max_multiple == 0 {
        return Ok(None);
    }

    let zero = Integer::from(0);
    let mut weighted_lower_gv_sum = Integer::from(0);
    let mut multiple = 1i128;
    while multiple <= max_multiple {
        let class = checked_ray_multiple(primitive_ray, multiple, "nilpotent-ray class")?;
        let gv = gv_by_class
            .get(&class)
            .cloned()
            .unwrap_or_else(|| zero.clone());
        if gv == zero && weighted_lower_gv_sum > zero {
            let first_vanishing_degree = degree.checked_mul(multiple).ok_or_else(|| {
                Error::InvalidInput("nilpotent-ray vanishing degree overflow".into())
            })?;
            return Ok(Some(NilpotentRayCandidate {
                primitive_ray: primitive_ray.to_vec(),
                first_vanishing_multiple: multiple,
                first_vanishing_degree,
                weighted_lower_gv_sum,
            }));
        }

        let weight = Integer::from(multiple) * Integer::from(multiple);
        weighted_lower_gv_sum += weight * gv;
        multiple += 1;
    }

    Ok(None)
}

/// Build the finite-cutoff set of apparently nilpotent primitive rays.
///
/// This is the set `N` construction in the McAllister moduli-space
/// reconstruction appendix. The supplied table represents GV data computed up
/// to `cutoff_degree`; nonzero effective classes seed primitive co-prime rays,
/// and each primitive ray is tested with
/// [`detect_apparent_nilpotent_ray_from_gv_multiples`].
pub fn detect_apparent_nilpotent_rays_from_gv_table(
    grading_vector: &[i64],
    cutoff_degree: i128,
    gv_invariants: &[(Vec<i64>, Integer)],
) -> Result<Vec<NilpotentRayCandidate>> {
    if grading_vector.is_empty() {
        return Err(Error::InvalidInput(
            "nilpotent-ray table scan requires a nonempty grading vector".into(),
        ));
    }
    if cutoff_degree <= 0 {
        return Err(Error::InvalidInput(
            "nilpotent-ray cutoff degree must be positive".into(),
        ));
    }

    let zero = Integer::from(0);
    let mut gv_by_class: HashMap<Vec<i64>, Integer> = HashMap::with_capacity(gv_invariants.len());
    let mut primitive_rays = BTreeSet::new();
    for (class, gv) in gv_invariants {
        if class.len() != grading_vector.len() {
            return Err(Error::InvalidInput(
                "nilpotent-ray GV invariant dimension does not match grading dimension".into(),
            ));
        }
        if let Some(existing) = gv_by_class.get(class) {
            if existing != gv {
                return Err(Error::InvalidInput(
                    "nilpotent-ray GV table contains conflicting duplicate classes".into(),
                ));
            }
            continue;
        }
        gv_by_class.insert(class.clone(), gv.clone());

        if gv == &zero {
            continue;
        }
        let degree = checked_i128_dot(class, grading_vector, "nilpotent-ray table grading")?;
        if degree <= 0 {
            return Err(Error::InvalidInput(format!(
                "nilpotent-ray nonzero GV class has non-positive grading degree {degree}"
            )));
        }
        if degree > cutoff_degree {
            continue;
        }
        let primitive = primitive_ray_from_curve_class(class, "nilpotent-ray nonzero GV class")?;
        primitive_rays.insert(primitive);
    }

    let mut candidates = Vec::new();
    for primitive_ray in primitive_rays {
        if let Some(candidate) = detect_apparent_nilpotent_ray_from_gv_multiples(
            &primitive_ray,
            grading_vector,
            cutoff_degree,
            gv_invariants,
        )? {
            candidates.push(candidate);
        }
    }
    Ok(candidates)
}

/// Partition finite-cutoff GV charges into apparently nilpotent rays and
/// apparently potent charges.
///
/// This prepares the appendix's `C \ N` input for the later nop-divergence
/// test. A nonzero charge is treated as nilpotent if its primitive ray is in
/// the detected nilpotent-ray set; otherwise it remains in `potent_charges`.
pub fn partition_finite_cutoff_gv_charges_by_nilpotence(
    grading_vector: &[i64],
    cutoff_degree: i128,
    gv_invariants: &[(Vec<i64>, Integer)],
) -> Result<FiniteCutoffGvChargePartition> {
    let nilpotent_rays =
        detect_apparent_nilpotent_rays_from_gv_table(grading_vector, cutoff_degree, gv_invariants)?;
    let nilpotent_primitive_rays: BTreeSet<Vec<i64>> = nilpotent_rays
        .iter()
        .map(|candidate| candidate.primitive_ray.clone())
        .collect();

    let zero = Integer::from(0);
    let mut gv_by_class: HashMap<Vec<i64>, Integer> = HashMap::with_capacity(gv_invariants.len());
    let mut potent_charges = Vec::new();
    for (class, gv) in gv_invariants {
        if let Some(existing) = gv_by_class.get(class) {
            if existing != gv {
                return Err(Error::InvalidInput(
                    "nilpotent-ray GV table contains conflicting duplicate classes".into(),
                ));
            }
            continue;
        }
        gv_by_class.insert(class.clone(), gv.clone());

        if gv == &zero {
            continue;
        }
        let degree = checked_i128_dot(class, grading_vector, "nilpotent-ray table grading")?;
        if degree > cutoff_degree {
            continue;
        }
        let primitive = primitive_ray_from_curve_class(class, "nilpotent-ray nonzero GV class")?;
        if !nilpotent_primitive_rays.contains(&primitive) {
            potent_charges.push((class.clone(), gv.clone()));
        }
    }
    potent_charges.sort_by(|lhs, rhs| lhs.0.cmp(&rhs.0));

    Ok(FiniteCutoffGvChargePartition {
        nilpotent_rays,
        potent_charges,
    })
}

/// Extract finite-cutoff nonzero GV charges after removing charges whose
/// primitive ray is in `excluded_primitive_rays`.
///
/// This is the finite-table complement used for the paper's second nop pass:
/// after the provisional nop set `F0` is known, the comparison set becomes
/// `C \ F0` rather than `C \ N`.
pub fn finite_cutoff_gv_charges_excluding_primitive_rays(
    grading_vector: &[i64],
    cutoff_degree: i128,
    gv_invariants: &[(Vec<i64>, Integer)],
    excluded_primitive_rays: &[Vec<i64>],
) -> Result<Vec<(Vec<i64>, Integer)>> {
    if grading_vector.is_empty() {
        return Err(Error::InvalidInput(
            "finite-cutoff GV complement requires a nonempty grading vector".into(),
        ));
    }
    if cutoff_degree <= 0 {
        return Err(Error::InvalidInput(
            "finite-cutoff GV complement cutoff degree must be positive".into(),
        ));
    }

    let mut excluded = BTreeSet::new();
    for ray in excluded_primitive_rays {
        if ray.len() != grading_vector.len() {
            return Err(Error::InvalidInput(
                "finite-cutoff GV complement excluded-ray dimension does not match grading dimension"
                    .into(),
            ));
        }
        let gcd = primitive_i64_gcd(ray, "finite-cutoff GV complement excluded ray")?;
        if gcd == 0 {
            return Err(Error::InvalidInput(
                "finite-cutoff GV complement excluded ray must not be zero".into(),
            ));
        }
        if gcd != 1 {
            return Err(Error::InvalidInput(format!(
                "finite-cutoff GV complement excluded ray must be co-prime, got gcd {gcd}"
            )));
        }
        excluded.insert(ray.clone());
    }

    let zero = Integer::from(0);
    let mut gv_by_class: HashMap<Vec<i64>, Integer> = HashMap::with_capacity(gv_invariants.len());
    let mut retained = Vec::new();
    for (class, gv) in gv_invariants {
        if class.len() != grading_vector.len() {
            return Err(Error::InvalidInput(
                "finite-cutoff GV complement class dimension does not match grading dimension"
                    .into(),
            ));
        }
        if let Some(existing) = gv_by_class.get(class) {
            if existing != gv {
                return Err(Error::InvalidInput(
                    "finite-cutoff GV complement table contains conflicting duplicate classes"
                        .into(),
                ));
            }
            continue;
        }
        gv_by_class.insert(class.clone(), gv.clone());

        if gv == &zero {
            continue;
        }
        let degree =
            checked_i128_dot(class, grading_vector, "finite-cutoff GV complement grading")?;
        if degree <= 0 {
            return Err(Error::InvalidInput(format!(
                "finite-cutoff GV complement nonzero class has non-positive grading degree {degree}"
            )));
        }
        if degree > cutoff_degree {
            continue;
        }
        let primitive = primitive_ray_from_curve_class(class, "finite-cutoff GV complement class")?;
        if !excluded.contains(&primitive) {
            retained.push((class.clone(), gv.clone()));
        }
    }

    retained.sort_by(|lhs, rhs| lhs.0.cmp(&rhs.0));
    Ok(retained)
}

/// Extract nonzero finite-GV charges that lie exactly on one grading-degree
/// slice.
///
/// This is a finite-table utility only. It is useful for diagnostics and for
/// wiring explicit slice-lattice inputs, but it is not by itself a certificate
/// that the complete same-degree slice lattice required by the nop test has
/// been generated.
pub fn finite_gv_nonzero_degree_slice_points(
    grading_vector: &[i64],
    slice_degree: i128,
    gv_invariants: &[(Vec<i64>, Integer)],
) -> Result<Vec<Vec<i64>>> {
    if grading_vector.is_empty() {
        return Err(Error::InvalidInput(
            "finite-GV degree-slice extraction requires a nonempty grading vector".into(),
        ));
    }
    if slice_degree <= 0 {
        return Err(Error::InvalidInput(
            "finite-GV degree-slice extraction requires a positive slice degree".into(),
        ));
    }

    let zero = Integer::from(0);
    let mut gv_by_class: HashMap<Vec<i64>, Integer> = HashMap::with_capacity(gv_invariants.len());
    let mut slice_points = BTreeSet::new();
    for (class, gv) in gv_invariants {
        if class.len() != grading_vector.len() {
            return Err(Error::InvalidInput(
                "finite-GV degree-slice class dimension does not match grading dimension".into(),
            ));
        }
        if let Some(existing) = gv_by_class.get(class) {
            if existing != gv {
                return Err(Error::InvalidInput(
                    "finite-GV degree-slice table contains conflicting duplicate classes".into(),
                ));
            }
            continue;
        }
        gv_by_class.insert(class.clone(), gv.clone());

        if gv == &zero {
            continue;
        }
        let degree = checked_i128_dot(class, grading_vector, "finite-GV degree-slice grading")?;
        if degree <= 0 {
            return Err(Error::InvalidInput(format!(
                "finite-GV degree-slice nonzero class has non-positive grading degree {degree}"
            )));
        }
        if degree == slice_degree {
            slice_points.insert(class.clone());
        }
    }

    Ok(slice_points.into_iter().collect())
}

/// Compute the appendix's degree-slice origin for a candidate nilpotent ray.
///
/// The nop-divergence test evaluates the candidate ray on a half-cutoff slice
/// and again on the full-cutoff slice. This helper computes the largest
/// positive integer `k` with
/// `k * C.v_g <= cutoff_degree * cutoff_numerator / cutoff_denominator`, and
/// returns the corresponding origin `k*C`. It does not compute the potent-ray
/// distance in the slice.
pub fn nilpotent_ray_degree_slice_for_cutoff_fraction(
    primitive_ray: &[i64],
    grading_vector: &[i64],
    cutoff_degree: i128,
    cutoff_numerator: i128,
    cutoff_denominator: i128,
) -> Result<Option<NilpotentRayDegreeSlice>> {
    if primitive_ray.is_empty() {
        return Err(Error::InvalidInput(
            "nilpotent-ray degree slice requires a nonempty primitive ray".into(),
        ));
    }
    if primitive_ray.len() != grading_vector.len() {
        return Err(Error::InvalidInput(
            "nilpotent-ray degree-slice grading dimension does not match ray dimension".into(),
        ));
    }
    if cutoff_degree <= 0 {
        return Err(Error::InvalidInput(
            "nilpotent-ray degree-slice cutoff must be positive".into(),
        ));
    }
    if cutoff_numerator <= 0 || cutoff_denominator <= 0 {
        return Err(Error::InvalidInput(
            "nilpotent-ray degree-slice cutoff fraction must be positive".into(),
        ));
    }

    let gcd = primitive_i64_gcd(primitive_ray, "nilpotent-ray degree-slice primitive ray")?;
    if gcd == 0 {
        return Err(Error::InvalidInput(
            "nilpotent-ray degree-slice primitive ray must not be zero".into(),
        ));
    }
    if gcd != 1 {
        return Err(Error::InvalidInput(format!(
            "nilpotent-ray degree-slice primitive ray must be co-prime, got gcd {gcd}"
        )));
    }

    let primitive_degree = checked_i128_dot(
        primitive_ray,
        grading_vector,
        "nilpotent-ray degree-slice grading",
    )?;
    if primitive_degree <= 0 {
        return Err(Error::InvalidInput(format!(
            "nilpotent-ray degree-slice primitive ray has non-positive grading degree {primitive_degree}"
        )));
    }

    let numerator = cutoff_degree
        .checked_mul(cutoff_numerator)
        .ok_or_else(|| Error::InvalidInput("nilpotent-ray degree-slice cutoff overflow".into()))?;
    let denominator = primitive_degree
        .checked_mul(cutoff_denominator)
        .ok_or_else(|| {
            Error::InvalidInput("nilpotent-ray degree-slice denominator overflow".into())
        })?;
    let slice_multiple = numerator / denominator;
    if slice_multiple == 0 {
        return Ok(None);
    }
    let slice_degree = primitive_degree
        .checked_mul(slice_multiple)
        .ok_or_else(|| Error::InvalidInput("nilpotent-ray degree-slice degree overflow".into()))?;
    let slice_origin =
        checked_ray_multiple(primitive_ray, slice_multiple, "nilpotent-ray degree slice")?;

    Ok(Some(NilpotentRayDegreeSlice {
        primitive_ray: primitive_ray.to_vec(),
        cutoff_numerator,
        cutoff_denominator,
        slice_multiple,
        slice_degree,
        slice_origin,
    }))
}

/// Enumerate comparison-ray integer points on a nilpotent degree slice.
///
/// For each supplied comparison charge, this reduces to its primitive ray and
/// keeps it only if an integer positive multiple lands exactly on
/// `slice.slice_degree`. Duplicate slice points are merged. The returned
/// offsets are the inputs to the later LLL-reduced infinity-norm distance
/// calculation; this helper does not compute that distance.
pub fn nilpotent_ray_slice_comparison_points(
    slice: &NilpotentRayDegreeSlice,
    comparison_charges: &[Vec<i64>],
    grading_vector: &[i64],
) -> Result<Vec<NilpotentRaySliceComparisonPoint>> {
    if slice.primitive_ray.len() != grading_vector.len()
        || slice.slice_origin.len() != grading_vector.len()
    {
        return Err(Error::InvalidInput(
            "nilpotent-ray slice dimension does not match grading dimension".into(),
        ));
    }
    if slice.slice_degree <= 0 {
        return Err(Error::InvalidInput(
            "nilpotent-ray slice degree must be positive".into(),
        ));
    }

    let mut by_point = BTreeMap::new();
    for charge in comparison_charges {
        if charge.len() != grading_vector.len() {
            return Err(Error::InvalidInput(
                "nilpotent-ray comparison charge dimension does not match grading dimension".into(),
            ));
        }
        let primitive = primitive_ray_from_curve_class(charge, "nilpotent-ray comparison charge")?;
        let primitive_degree = checked_i128_dot(
            &primitive,
            grading_vector,
            "nilpotent-ray comparison grading",
        )?;
        if primitive_degree <= 0 {
            return Err(Error::InvalidInput(format!(
                "nilpotent-ray comparison charge has non-positive grading degree {primitive_degree}"
            )));
        }
        if slice.slice_degree % primitive_degree != 0 {
            continue;
        }
        let slice_multiple = slice.slice_degree / primitive_degree;
        let slice_point =
            checked_ray_multiple(&primitive, slice_multiple, "nilpotent-ray comparison point")?;
        let offset_from_origin = checked_i64_vector_difference(
            &slice_point,
            &slice.slice_origin,
            "nilpotent-ray comparison offset",
        )?;
        by_point
            .entry(slice_point.clone())
            .or_insert_with(|| NilpotentRaySliceComparisonPoint {
                primitive_ray: primitive,
                primitive_degree,
                slice_multiple,
                slice_point,
                offset_from_origin,
            });
    }

    Ok(by_point.into_values().collect())
}

/// Compute the LLL-reduced integer distance from a candidate nilpotent ray to
/// comparison rays on one finite-cutoff degree slice.
///
/// The caller must provide the same-degree curve classes that define the
/// affine slice lattice. This keeps the function source-bounded: it implements
/// the paper's LLL/norm calculation once a lattice slice is known, but it does
/// not guess whether that slice came from a compact GV semigroup, a chamber
/// continuation, or another certified source.
pub fn nilpotent_ray_lll_reduced_slice_distance(
    slice: &NilpotentRayDegreeSlice,
    slice_lattice_points: &[Vec<i64>],
    comparison_charges: &[Vec<i64>],
    grading_vector: &[i64],
) -> Result<NilpotentRaySliceDistance> {
    validate_nilpotent_slice_against_grading(slice, grading_vector)?;

    let mut lattice_offsets = BTreeSet::new();
    for point in slice_lattice_points {
        if point.len() != grading_vector.len() {
            return Err(Error::InvalidInput(
                "nilpotent-ray slice lattice point dimension does not match grading dimension"
                    .into(),
            ));
        }
        let degree = checked_i128_dot(
            point,
            grading_vector,
            "nilpotent-ray slice lattice point grading",
        )?;
        if degree != slice.slice_degree {
            return Err(Error::InvalidInput(format!(
                "nilpotent-ray slice lattice point has degree {degree}, expected {}",
                slice.slice_degree
            )));
        }
        let offset = checked_i64_vector_difference(
            point,
            &slice.slice_origin,
            "nilpotent-ray slice lattice offset",
        )?;
        if offset.iter().any(|&value| value != 0) {
            lattice_offsets.insert(offset);
        }
    }

    let lattice_offsets: Vec<Vec<i64>> = lattice_offsets.into_iter().collect();
    if lattice_offsets.is_empty() {
        return Err(Error::InvalidInput(
            "nilpotent-ray slice lattice requires at least one nonzero offset".into(),
        ));
    }

    let lll = lll_reduce(&lattice_offsets, true);
    let lll_transform = lll.transform.ok_or_else(|| {
        Error::InvalidInput("nilpotent-ray LLL reduction did not return a transform".into())
    })?;
    let reduced_lattice_offsets = lll.reduced;

    for (input, reduced) in lattice_offsets.iter().zip(reduced_lattice_offsets.iter()) {
        let transformed =
            apply_i64_matrix_to_vector(&lll_transform, input, "nilpotent-ray LLL transform")?;
        if transformed != *reduced {
            return Err(Error::InvalidInput(
                "nilpotent-ray LLL transform does not reproduce reduced offsets".into(),
            ));
        }
    }

    let comparison_points =
        nilpotent_ray_slice_comparison_points(slice, comparison_charges, grading_vector)?;
    let mut transformed_comparison_offsets = Vec::with_capacity(comparison_points.len());
    let mut minimum_infinity_norm: Option<i64> = None;
    for point in &comparison_points {
        let transformed = apply_i64_matrix_to_vector(
            &lll_transform,
            &point.offset_from_origin,
            "nilpotent-ray comparison LLL transform",
        )?;
        let norm =
            checked_i64_infinity_norm(&transformed, "nilpotent-ray transformed comparison offset")?;
        minimum_infinity_norm =
            Some(minimum_infinity_norm.map_or(norm, |current| current.min(norm)));
        transformed_comparison_offsets.push(transformed);
    }

    Ok(NilpotentRaySliceDistance {
        slice: slice.clone(),
        lattice_offsets,
        lll_transform,
        reduced_lattice_offsets,
        comparison_points,
        transformed_comparison_offsets,
        minimum_infinity_norm,
    })
}

/// Compare the paper's half-cutoff and full-cutoff distances for one candidate
/// nilpotent ray.
pub fn nilpotent_ray_divergence_check_from_slice_distances(
    half_cutoff: NilpotentRaySliceDistance,
    full_cutoff: NilpotentRaySliceDistance,
) -> Result<NilpotentRayDivergenceCheck> {
    if half_cutoff.slice.primitive_ray != full_cutoff.slice.primitive_ray {
        return Err(Error::InvalidInput(
            "nilpotent-ray divergence check requires distances for the same primitive ray".into(),
        ));
    }
    let appears_divergent = match (
        half_cutoff.minimum_infinity_norm,
        full_cutoff.minimum_infinity_norm,
    ) {
        (Some(half_distance), Some(full_distance)) => Some(full_distance > half_distance),
        _ => None,
    };

    Ok(NilpotentRayDivergenceCheck {
        half_cutoff,
        full_cutoff,
        appears_divergent,
    })
}

/// Run the paper's half/full finite-cutoff distance comparison for one
/// candidate ray, using explicit same-degree slice lattices supplied by the
/// caller.
///
/// The returned `None` means the candidate ray does not reach the half or full
/// cutoff slice with a positive integer multiple. This function intentionally
/// does not build the slice lattices: those must come from a certified GV
/// semigroup/chamber source.
pub fn nilpotent_ray_divergence_check_with_explicit_slice_lattices(
    primitive_ray: &[i64],
    grading_vector: &[i64],
    cutoff_degree: i128,
    half_slice_lattice_points: &[Vec<i64>],
    full_slice_lattice_points: &[Vec<i64>],
    comparison_charges: &[Vec<i64>],
) -> Result<Option<NilpotentRayDivergenceCheck>> {
    let Some(half_slice) = nilpotent_ray_degree_slice_for_cutoff_fraction(
        primitive_ray,
        grading_vector,
        cutoff_degree,
        1,
        2,
    )?
    else {
        return Ok(None);
    };
    let Some(full_slice) = nilpotent_ray_degree_slice_for_cutoff_fraction(
        primitive_ray,
        grading_vector,
        cutoff_degree,
        1,
        1,
    )?
    else {
        return Ok(None);
    };

    let half_distance = nilpotent_ray_lll_reduced_slice_distance(
        &half_slice,
        half_slice_lattice_points,
        comparison_charges,
        grading_vector,
    )?;
    let full_distance = nilpotent_ray_lll_reduced_slice_distance(
        &full_slice,
        full_slice_lattice_points,
        comparison_charges,
        grading_vector,
    )?;
    nilpotent_ray_divergence_check_from_slice_distances(half_distance, full_distance).map(Some)
}

/// Apply the paper's two-pass nop classification to already-computed
/// divergence checks.
///
/// The first pass is interpreted as comparison against `C \ N`, producing
/// provisional nop rays `F0`. The second pass is interpreted as comparison
/// against `C \ F0`, and is required exactly for the rays in `F0`.
pub fn classify_nilpotent_rays_from_two_pass_divergence_checks(
    nilpotent_rays: &[NilpotentRayCandidate],
    first_pass_checks: &[(Vec<i64>, NilpotentRayDivergenceCheck)],
    second_pass_checks: &[(Vec<i64>, NilpotentRayDivergenceCheck)],
) -> Result<NilpotentRayTwoPassNopClassification> {
    let mut nilpotent_by_ray = BTreeMap::new();
    for candidate in nilpotent_rays {
        let gcd = primitive_i64_gcd(
            &candidate.primitive_ray,
            "two-pass nop classification nilpotent ray",
        )?;
        if gcd == 0 {
            return Err(Error::InvalidInput(
                "two-pass nop classification nilpotent ray must not be zero".into(),
            ));
        }
        if gcd != 1 {
            return Err(Error::InvalidInput(format!(
                "two-pass nop classification nilpotent ray must be co-prime, got gcd {gcd}"
            )));
        }
        if nilpotent_by_ray
            .insert(candidate.primitive_ray.clone(), candidate)
            .is_some()
        {
            return Err(Error::InvalidInput(
                "two-pass nop classification contains duplicate nilpotent rays".into(),
            ));
        }
    }

    let first_pass_by_ray =
        collect_nilpotent_divergence_checks_by_ray(first_pass_checks, "first-pass")?;
    for ray in first_pass_by_ray.keys() {
        if !nilpotent_by_ray.contains_key(ray) {
            return Err(Error::InvalidInput(
                "two-pass nop classification first-pass check references a non-nilpotent ray"
                    .into(),
            ));
        }
    }

    let mut initial_candidate_nop_rays = Vec::new();
    let mut first_pass_inconclusive_rays = Vec::new();
    for candidate in nilpotent_rays {
        let check = first_pass_by_ray
            .get(&candidate.primitive_ray)
            .ok_or_else(|| {
                Error::InvalidInput(
                    "two-pass nop classification is missing a first-pass check".into(),
                )
            })?;
        match check.appears_divergent {
            Some(true) => initial_candidate_nop_rays.push(candidate.clone()),
            Some(false) => {}
            None => first_pass_inconclusive_rays.push(candidate.clone()),
        }
    }

    let f0_by_ray: BTreeMap<Vec<i64>, &NilpotentRayCandidate> = initial_candidate_nop_rays
        .iter()
        .map(|candidate| (candidate.primitive_ray.clone(), candidate))
        .collect();
    let second_pass_by_ray =
        collect_nilpotent_divergence_checks_by_ray(second_pass_checks, "second-pass")?;
    for ray in second_pass_by_ray.keys() {
        if !f0_by_ray.contains_key(ray) {
            return Err(Error::InvalidInput(
                "two-pass nop classification second-pass check references a non-F0 ray".into(),
            ));
        }
    }

    let mut nop_rays = Vec::new();
    let mut second_pass_inconclusive_rays = Vec::new();
    for candidate in &initial_candidate_nop_rays {
        let check = second_pass_by_ray
            .get(&candidate.primitive_ray)
            .ok_or_else(|| {
                Error::InvalidInput(
                    "two-pass nop classification is missing a second-pass check".into(),
                )
            })?;
        match check.appears_divergent {
            Some(true) => nop_rays.push(candidate.clone()),
            Some(false) => {}
            None => second_pass_inconclusive_rays.push(candidate.clone()),
        }
    }

    Ok(NilpotentRayTwoPassNopClassification {
        initial_candidate_nop_rays,
        nop_rays,
        first_pass_inconclusive_rays,
        second_pass_inconclusive_rays,
    })
}

/// Run the paper's finite-cutoff two-pass nop classification from a supplied
/// finite GV table.
///
/// This function assumes `gv_invariants` is already the caller's computed
/// finite table up to `cutoff_degree`. It does not compute GV invariants or
/// certify that the table came from the correct chamber/semigroup; it only
/// applies the appendix algorithm to that explicit finite source.
pub fn classify_nop_rays_from_finite_gv_table(
    grading_vector: &[i64],
    cutoff_degree: i128,
    gv_invariants: &[(Vec<i64>, Integer)],
) -> Result<FiniteGvTableNopClassification> {
    let partition = partition_finite_cutoff_gv_charges_by_nilpotence(
        grading_vector,
        cutoff_degree,
        gv_invariants,
    )?;
    let first_pass_comparison_charges: Vec<Vec<i64>> = partition
        .potent_charges
        .iter()
        .map(|(charge, _)| charge.clone())
        .collect();

    let mut first_pass_checks = Vec::with_capacity(partition.nilpotent_rays.len());
    for candidate in &partition.nilpotent_rays {
        let check = finite_gv_table_divergence_check_for_candidate(
            &candidate.primitive_ray,
            grading_vector,
            cutoff_degree,
            gv_invariants,
            &first_pass_comparison_charges,
            "first-pass",
        )?;
        first_pass_checks.push((candidate.primitive_ray.clone(), check));
    }

    let f0_candidates: Vec<NilpotentRayCandidate> = partition
        .nilpotent_rays
        .iter()
        .zip(first_pass_checks.iter())
        .filter_map(|(candidate, (_, check))| {
            if check.appears_divergent == Some(true) {
                Some(candidate.clone())
            } else {
                None
            }
        })
        .collect();
    let f0_primitive_rays: Vec<Vec<i64>> = f0_candidates
        .iter()
        .map(|candidate| candidate.primitive_ray.clone())
        .collect();
    let second_pass_comparison_charges = finite_cutoff_gv_charges_excluding_primitive_rays(
        grading_vector,
        cutoff_degree,
        gv_invariants,
        &f0_primitive_rays,
    )?;
    let second_pass_comparison_charge_rows: Vec<Vec<i64>> = second_pass_comparison_charges
        .iter()
        .map(|(charge, _)| charge.clone())
        .collect();

    let mut second_pass_checks = Vec::with_capacity(f0_candidates.len());
    for candidate in &f0_candidates {
        let check = finite_gv_table_divergence_check_for_candidate(
            &candidate.primitive_ray,
            grading_vector,
            cutoff_degree,
            gv_invariants,
            &second_pass_comparison_charge_rows,
            "second-pass",
        )?;
        second_pass_checks.push((candidate.primitive_ray.clone(), check));
    }

    let classification = classify_nilpotent_rays_from_two_pass_divergence_checks(
        &partition.nilpotent_rays,
        &first_pass_checks,
        &second_pass_checks,
    )?;

    Ok(FiniteGvTableNopClassification {
        partition,
        first_pass_checks,
        second_pass_comparison_charges,
        second_pass_checks,
        classification,
    })
}

/// Compute the paper's potent-ray convergence terms.
///
/// For a ray `q` with volume `q.t`, the paper defines
/// `xi_n = N_{nq} exp(-2π n q.t)`. This function returns
/// `log(|xi_n|)` for each supplied GV invariant `N_{nq}`.
pub fn potent_ray_log_xi_terms(
    gv_series: &[Integer],
    curve_volume: F64<Pos>,
) -> Result<Vec<Option<F64<Finite>>>> {
    if gv_series.is_empty() {
        return Err(Error::InvalidInput(
            "potent-ray GV series must not be empty".into(),
        ));
    }

    let mut terms = Vec::with_capacity(gv_series.len());
    for (idx, gv) in gv_series.iter().enumerate() {
        let Some(log_gv) = integer_abs_ln(gv)? else {
            terms.push(None);
            continue;
        };
        let degree = (idx + 1) as f64;
        let log_xi = log_gv - 2.0 * PI * degree * curve_volume.get();
        let finite = F64::<Finite>::new(log_xi)
            .ok_or_else(|| Error::InvalidInput("potent-ray log-xi term is not finite".into()))?;
        terms.push(Some(finite));
    }
    Ok(terms)
}

/// Compute a least-squares slope for finite potent-ray log-xi terms.
pub fn potent_ray_log_xi_slope(log_xi_terms: &[Option<F64<Finite>>]) -> Option<F64<Finite>> {
    let mut n = 0.0;
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_xx = 0.0;
    let mut sum_xy = 0.0;

    for (idx, term) in log_xi_terms.iter().enumerate() {
        let Some(term) = term else {
            continue;
        };
        let x = (idx + 1) as f64;
        let y = term.get();
        n += 1.0;
        sum_x += x;
        sum_y += y;
        sum_xx += x * x;
        sum_xy += x * y;
    }

    if n < 2.0 {
        return None;
    }
    let denom = (n * sum_xx) - (sum_x * sum_x);
    if denom == 0.0 {
        return None;
    }
    F64::<Finite>::new(((n * sum_xy) - (sum_x * sum_y)) / denom)
}

/// Summarize potent-ray convergence from a GV series and curve volume.
pub fn potent_ray_convergence(
    gv_series: &[Integer],
    curve_volume: F64<Pos>,
) -> Result<PotentRayConvergence> {
    let log_xi_terms = potent_ray_log_xi_terms(gv_series, curve_volume)?;
    let log_xi_slope = potent_ray_log_xi_slope(&log_xi_terms);
    Ok(PotentRayConvergence {
        log_xi_terms,
        log_xi_slope,
    })
}

/// Project an ambient curve class to coordinates in a divisor basis.
///
/// The McAllister ancillary curve files store ambient divisor-intersection
/// vectors, while `cygv` expects Kähler-basis curve coordinates. For a vector
/// divisor basis this projection is the CYTools convention used throughout the
/// McAllister pipeline: keep the entries at the selected basis divisor indices.
pub fn project_ambient_curve_to_basis(ambient_class: &[i64], basis: &[usize]) -> Result<Vec<i64>> {
    basis
        .iter()
        .map(|&idx| {
            ambient_class.get(idx).copied().ok_or_else(|| {
                Error::InvalidInput(format!(
                    "basis index {idx} is out of bounds for ambient curve dimension {}",
                    ambient_class.len()
                ))
            })
        })
        .collect()
}

/// Project an ambient curve class to coordinates in a matrix divisor basis.
///
/// This is the single-row analogue of CYTools'
/// `mori_cap_matrix.dot(basis.T)` generic-basis projection. Rows of
/// `basis_matrix` are divisor-basis vectors in ambient coordinates.
pub fn project_ambient_curve_to_basis_matrix(
    ambient_class: &[i64],
    basis_matrix: &[Vec<i64>],
) -> Result<Vec<i64>> {
    validate_basis_matrix(basis_matrix, ambient_class.len())?;
    basis_matrix
        .iter()
        .map(|basis_row| checked_i64_dot(ambient_class, basis_row))
        .collect()
}

fn validate_basis_matrix(basis_matrix: &[Vec<i64>], ambient_dim: usize) -> Result<()> {
    if basis_matrix.is_empty() {
        return Err(Error::InvalidInput("basis matrix is empty".into()));
    }
    if ambient_dim == 0 {
        return Err(Error::InvalidInput(
            "ambient curve dimension is empty".into(),
        ));
    }
    for row in basis_matrix {
        if row.len() != ambient_dim {
            return Err(Error::InvalidInput(format!(
                "basis matrix row dimension {} does not match ambient curve dimension {ambient_dim}",
                row.len()
            )));
        }
    }
    Ok(())
}

fn checked_i64_dot(left: &[i64], right: &[i64]) -> Result<i64> {
    let mut acc = 0i128;
    for (&a, &b) in left.iter().zip(right.iter()) {
        let term = i128::from(a) * i128::from(b);
        acc = acc
            .checked_add(term)
            .ok_or_else(|| Error::InvalidInput("matrix-basis projection overflowed i128".into()))?;
    }
    i64::try_from(acc)
        .map_err(|_| Error::InvalidInput("matrix-basis projection does not fit in i64".into()))
}

/// Check an integer normal as an exact supporting-face certificate.
///
/// This function does not try to find a normal. It verifies a proposed normal
/// using integer arithmetic only:
///
/// - every supplied face generator must pair to zero;
/// - every Mori generator must pair non-negatively;
/// - at least one Mori generator must pair positively.
///
/// A failed certificate returns `Ok(None)`. Malformed dimensions or arithmetic
/// overflow return an error.
pub fn check_supporting_mori_face_normal(
    normal: &[i64],
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
) -> Result<Option<SupportingMoriFaceCertificate>> {
    if normal.is_empty() {
        return Err(Error::InvalidInput(
            "supporting Mori face normal is empty".into(),
        ));
    }
    if mori_generators.is_empty() {
        return Err(Error::InvalidInput(
            "supporting Mori face check requires Mori generators".into(),
        ));
    }
    if normal.iter().all(|&value| value == 0) {
        return Ok(None);
    }

    for generator in face_generators {
        validate_curve_dimension("face generator", generator, normal.len())?;
        if exact_i64_dot_checked(normal, generator)? != 0 {
            return Ok(None);
        }
    }

    let mut zero_generator_count = 0usize;
    let mut positive_generator_count = 0usize;
    for generator in mori_generators {
        validate_curve_dimension("Mori generator", generator, normal.len())?;
        let dot = exact_i64_dot_checked(normal, generator)?;
        if dot < 0 {
            return Ok(None);
        }
        if dot == 0 {
            zero_generator_count += 1;
        } else {
            positive_generator_count += 1;
        }
    }

    if positive_generator_count == 0 {
        return Ok(None);
    }

    Ok(Some(SupportingMoriFaceCertificate {
        normal: normal.to_vec(),
        zero_generator_count,
        positive_generator_count,
    }))
}

/// Find an exact supporting normal from the integer kernel of face generators.
///
/// This is deliberately conservative: it only succeeds for codimension-one
/// faces where the supplied face generators have a one-dimensional integer
/// kernel. Higher-codimension faces need an additional exact selection step
/// inside that kernel, so this returns `Ok(None)` for them instead of guessing.
///
/// The found normal is then verified by [`check_supporting_mori_face_normal`],
/// so malformed dimensions and non-supporting candidate faces fail loudly or
/// return `Ok(None)` using the same integer arithmetic as a caller-supplied
/// certificate.
pub fn certify_supporting_mori_face_by_exact_kernel(
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
) -> Result<Option<SupportingMoriFaceCertificate>> {
    let Some(first_mori_generator) = mori_generators.first() else {
        return Err(Error::InvalidInput(
            "supporting Mori face check requires Mori generators".into(),
        ));
    };
    let dim = first_mori_generator.len();
    if dim == 0 {
        return Err(Error::InvalidInput(
            "supporting Mori face generator dimension is zero".into(),
        ));
    }
    for generator in mori_generators {
        validate_curve_dimension("Mori generator", generator, dim)?;
    }
    if face_generators.is_empty() {
        return Ok(None);
    }

    let mut face_matrix = Vec::with_capacity(face_generators.len());
    for generator in face_generators {
        validate_curve_dimension("face generator", generator, dim)?;
        face_matrix.push(
            generator
                .iter()
                .map(|&entry| Integer::from(entry))
                .collect::<Vec<_>>(),
        );
    }

    let kernel = integer_kernel(&face_matrix);
    if kernel.len() != 1 {
        return Ok(None);
    }
    let normal =
        integer_vector_to_i64(&kernel[0], "supporting Mori face exact kernel normal entry")?;
    if let Some(certificate) =
        check_supporting_mori_face_normal(&normal, face_generators, mori_generators)?
    {
        return Ok(Some(certificate));
    }

    let opposite_normal = checked_neg_i64_vector(&normal, "supporting Mori face normal")?;
    check_supporting_mori_face_normal(&opposite_normal, face_generators, mori_generators)
}

/// Search for an exact supporting Mori face certificate using LP candidates.
///
/// This handles higher-codimension faces where the integer kernel has dimension
/// greater than one. The LP pass chooses a candidate normal inside that kernel,
/// then Cyrus rounds/scales candidate normals and verifies them exactly with
/// [`check_supporting_mori_face_normal`]. A successful return is therefore an
/// exact certificate; `Ok(None)` only means this bounded LP search did not find
/// one.
pub fn certify_supporting_mori_face_by_lp_search(
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<Option<SupportingMoriFaceCertificate>> {
    Ok(
        diagnose_supporting_mori_face_by_lp_search(face_generators, mori_generators, options)?
            .certificate,
    )
}

/// Run the LP-assisted supporting-face search and report each failed phase.
///
/// This is a diagnostic companion to [`certify_supporting_mori_face_by_lp_search`].
/// It exists because a missing certificate can mean several different things:
/// no exact codimension-one kernel, no feasible aggregate LP normal, a real
/// normal that cannot be rounded to an exact integer certificate, or the same
/// outcomes across bounded anchor attempts.
pub fn diagnose_supporting_mori_face_by_lp_search(
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<SupportingMoriFaceLpSearchDiagnostic> {
    validate_supporting_face_lp_inputs(face_generators, mori_generators, options)?;
    let dim = mori_generators[0].len();
    let face_rank = if face_generators.is_empty() {
        0
    } else {
        curve_row_span_rank(face_generators)?
    };
    let mut diagnostic = SupportingMoriFaceLpSearchDiagnostic {
        face_rank,
        dim,
        status: String::new(),
        exact_kernel_status: String::new(),
        full_status: "not_attempted".to_string(),
        aggregate_status: "not_attempted".to_string(),
        anchor_attempt_count: 0,
        anchor_lp_solution_count: 0,
        anchor_status_counts: BTreeMap::new(),
        certificate: None,
    };

    if face_generators.is_empty() {
        diagnostic.status = "empty_support".to_string();
        diagnostic.exact_kernel_status = "not_attempted_empty_support".to_string();
        return Ok(diagnostic);
    }

    match certify_supporting_mori_face_by_exact_kernel(face_generators, mori_generators)? {
        Some(certificate) => {
            diagnostic.status = "certified_exact_kernel".to_string();
            diagnostic.exact_kernel_status = "certified".to_string();
            diagnostic.certificate = Some(certificate);
            return Ok(diagnostic);
        }
        None => {
            diagnostic.exact_kernel_status = "no_certificate".to_string();
        }
    }

    match solve_supporting_face_normal_full_lp(face_generators, mori_generators, options)? {
        SupportingFaceLpNormalSearchOutcome::Found(lp_normal) => {
            if let Some(certificate) = integer_supporting_face_certificate_from_lp(
                &lp_normal,
                face_generators,
                mori_generators,
                options,
            )? {
                diagnostic.status = "certified_full_lp".to_string();
                diagnostic.full_status = "certified".to_string();
                diagnostic.certificate = Some(certificate);
                return Ok(diagnostic);
            }
            diagnostic.full_status = "lp_solution_rounding_no_certificate".to_string();
        }
        other => {
            diagnostic.full_status = other.status();
        }
    }

    match solve_supporting_face_normal_aggregate_lp_detailed(
        face_generators,
        mori_generators,
        options,
    )? {
        SupportingFaceLpNormalSearchOutcome::Found(lp_normal) => {
            if let Some(certificate) = integer_supporting_face_certificate_from_lp(
                &lp_normal,
                face_generators,
                mori_generators,
                options,
            )? {
                diagnostic.status = "certified_aggregate_lp".to_string();
                diagnostic.aggregate_status = "certified".to_string();
                diagnostic.certificate = Some(certificate);
                return Ok(diagnostic);
            }
            diagnostic.aggregate_status = "lp_solution_rounding_no_certificate".to_string();
        }
        other => {
            diagnostic.aggregate_status = other.status();
        }
    }

    let anchors = supporting_face_anchor_candidates(face_generators, mori_generators, options)?;
    for anchor in anchors {
        diagnostic.anchor_attempt_count += 1;
        match solve_supporting_face_normal_lp_detailed(
            face_generators,
            mori_generators,
            &anchor,
            options,
        )? {
            SupportingFaceLpNormalSearchOutcome::Found(lp_normal) => {
                diagnostic.anchor_lp_solution_count += 1;
                if let Some(certificate) = integer_supporting_face_certificate_from_lp(
                    &lp_normal,
                    face_generators,
                    mori_generators,
                    options,
                )? {
                    *diagnostic
                        .anchor_status_counts
                        .entry("certified".to_string())
                        .or_insert(0) += 1;
                    diagnostic.status = "certified_anchor_lp".to_string();
                    diagnostic.certificate = Some(certificate);
                    return Ok(diagnostic);
                }
                *diagnostic
                    .anchor_status_counts
                    .entry("lp_solution_rounding_no_certificate".to_string())
                    .or_insert(0) += 1;
            }
            other => {
                *diagnostic
                    .anchor_status_counts
                    .entry(other.status())
                    .or_insert(0) += 1;
            }
        }
    }

    diagnostic.status = "lp_no_certificate".to_string();
    Ok(diagnostic)
}

fn validate_supporting_face_lp_inputs(
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<()> {
    let Some(first_mori_generator) = mori_generators.first() else {
        return Err(Error::InvalidInput(
            "supporting Mori face LP search requires Mori generators".into(),
        ));
    };
    let dim = first_mori_generator.len();
    if dim == 0 {
        return Err(Error::InvalidInput(
            "supporting Mori face LP search dimension is zero".into(),
        ));
    }
    for generator in face_generators {
        validate_curve_dimension("face generator", generator, dim)?;
    }
    for generator in mori_generators {
        validate_curve_dimension("Mori generator", generator, dim)?;
    }
    if options.cutting_rounds == 0 {
        return Err(Error::InvalidInput(
            "supporting Mori face LP search cutting rounds must be positive".into(),
        ));
    }
    if options.scale_limit <= 0 {
        return Err(Error::InvalidInput(
            "supporting Mori face LP search scale limit must be positive".into(),
        ));
    }
    if !options.variable_bound.is_finite() || options.variable_bound <= 0.0 {
        return Err(Error::InvalidInput(
            "supporting Mori face LP search variable bound must be positive and finite".into(),
        ));
    }
    Ok(())
}

fn supporting_face_anchor_candidates(
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<Vec<Vec<i64>>> {
    if face_generators.is_empty() {
        return Ok(Vec::new());
    }
    let dim = face_generators[0].len();
    let face_rank = curve_row_span_rank(face_generators)?;
    if face_rank >= dim {
        return Ok(Vec::new());
    }

    let mut anchors = Vec::new();
    let mut seen = HashSet::new();
    for ray in mori_generators {
        if face_generators.iter().any(|generator| generator == ray) || !seen.insert(ray.clone()) {
            continue;
        }
        let mut rows = face_generators.to_vec();
        rows.push(ray.clone());
        if curve_row_span_rank(&rows)? > face_rank {
            anchors.push(ray.clone());
            if anchors.len() >= options.anchor_attempts {
                break;
            }
        }
    }
    Ok(anchors)
}

enum SupportingFaceLpNormalSearchOutcome {
    Found(Vec<f64>),
    LpNoSolution,
    RepeatedViolation,
    CuttingRoundLimit,
    SolverOther(String),
}

fn lp_status_error_fragment(error: &str) -> String {
    let mut out = String::new();
    for ch in error.chars() {
        let mapped = if ch.is_ascii_alphanumeric() {
            ch.to_ascii_lowercase()
        } else {
            '_'
        };
        if mapped == '_' && out.ends_with('_') {
            continue;
        }
        out.push(mapped);
    }
    out.trim_matches('_').to_string()
}

impl SupportingFaceLpNormalSearchOutcome {
    fn status(&self) -> String {
        match self {
            Self::Found(_) => "lp_solution".to_string(),
            Self::LpNoSolution => "lp_no_solution".to_string(),
            Self::RepeatedViolation => "lp_repeated_violation".to_string(),
            Self::CuttingRoundLimit => "lp_cutting_round_limit".to_string(),
            Self::SolverOther(message) => {
                format!("lp_solver_other_{}", lp_status_error_fragment(message))
            }
        }
    }
}

fn solve_supporting_face_normal_full_lp(
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<SupportingFaceLpNormalSearchOutcome> {
    let aggregate = aggregate_mori_ray_coefficients(mori_generators)?;
    let normal_vars = supporting_face_normal_vars(aggregate.len(), options)?;
    let vars = normal_vars.variables;
    let mut objective = Expression::from(0.0);
    objective.add_mul(0.0, normal_vars.normal[0]);
    let mut model = vars.minimise(objective).using(default_solver);

    model = add_supporting_face_zero_constraints(model, &normal_vars.normal, face_generators);

    let mut aggregate_expr = Expression::from(0.0);
    for (var, &coefficient) in normal_vars.normal.iter().zip(&aggregate) {
        if coefficient != 0 {
            aggregate_expr.add_mul(coefficient as f64, *var);
        }
    }
    model = model.with(aggregate_expr.geq(1.0));

    let all_ray_indices = (0..mori_generators.len()).collect::<Vec<_>>();
    model = add_supporting_face_enforced_ray_constraints(
        model,
        &normal_vars.normal,
        mori_generators,
        &all_ray_indices,
    )?;

    solve_supporting_face_normal_lp_model(model, &normal_vars.normal)
}

fn solve_supporting_face_normal_aggregate_lp_detailed(
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<SupportingFaceLpNormalSearchOutcome> {
    let aggregate = aggregate_mori_ray_coefficients(mori_generators)?;
    let mut enforced_ray_indices = Vec::new();
    let mut enforced_ray_set = HashSet::new();
    for _ in 0..options.cutting_rounds {
        let normal = match solve_supporting_face_normal_aggregate_lp_with_enforced_rays(
            face_generators,
            mori_generators,
            &aggregate,
            &enforced_ray_indices,
            options,
        )? {
            SupportingFaceLpNormalSearchOutcome::Found(normal) => normal,
            other => return Ok(other),
        };
        let Some(violating_idx) = most_negative_lp_normal_violation(&normal, mori_generators)
        else {
            return Ok(SupportingFaceLpNormalSearchOutcome::Found(normal));
        };
        if !enforced_ray_set.insert(violating_idx) {
            return Ok(SupportingFaceLpNormalSearchOutcome::RepeatedViolation);
        }
        enforced_ray_indices.push(violating_idx);
    }
    Ok(SupportingFaceLpNormalSearchOutcome::CuttingRoundLimit)
}

fn aggregate_mori_ray_coefficients(mori_generators: &[Vec<i64>]) -> Result<Vec<i128>> {
    let Some(first_ray) = mori_generators.first() else {
        return Err(Error::InvalidInput(
            "supporting Mori face aggregate normal LP requires Mori generators".into(),
        ));
    };
    let dim = first_ray.len();
    let mut aggregate = vec![0_i128; dim];
    for ray in mori_generators {
        validate_curve_dimension("Mori generator", ray, dim)?;
        for (slot, &coefficient) in aggregate.iter_mut().zip(ray) {
            *slot = (*slot)
                .checked_add(i128::from(coefficient))
                .ok_or_else(|| {
                    Error::InvalidInput(
                        "supporting Mori face aggregate normal LP coefficient overflowed".into(),
                    )
                })?;
        }
    }
    Ok(aggregate)
}

fn solve_supporting_face_normal_aggregate_lp_with_enforced_rays(
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
    aggregate: &[i128],
    enforced_ray_indices: &[usize],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<SupportingFaceLpNormalSearchOutcome> {
    let normal_vars = supporting_face_normal_vars(aggregate.len(), options)?;
    let vars = normal_vars.variables;
    let mut objective = Expression::from(0.0);
    objective.add_mul(0.0, normal_vars.normal[0]);
    let mut model = vars.minimise(objective).using(default_solver);

    model = add_supporting_face_zero_constraints(model, &normal_vars.normal, face_generators);

    let mut aggregate_expr = Expression::from(0.0);
    for (var, &coefficient) in normal_vars.normal.iter().zip(aggregate) {
        if coefficient != 0 {
            aggregate_expr.add_mul(coefficient as f64, *var);
        }
    }
    model = model.with(aggregate_expr.geq(1.0));

    model = add_supporting_face_enforced_ray_constraints(
        model,
        &normal_vars.normal,
        mori_generators,
        enforced_ray_indices,
    )?;
    solve_supporting_face_normal_lp_model(model, &normal_vars.normal)
}

fn solve_supporting_face_normal_lp_detailed(
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
    anchor: &[i64],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<SupportingFaceLpNormalSearchOutcome> {
    let mut enforced_ray_indices = Vec::new();
    let mut enforced_ray_set = HashSet::new();
    for _ in 0..options.cutting_rounds {
        let normal = match solve_supporting_face_normal_lp_with_enforced_rays(
            face_generators,
            mori_generators,
            anchor,
            &enforced_ray_indices,
            options,
        )? {
            SupportingFaceLpNormalSearchOutcome::Found(normal) => normal,
            other => return Ok(other),
        };
        let Some(violating_idx) = most_negative_lp_normal_violation(&normal, mori_generators)
        else {
            return Ok(SupportingFaceLpNormalSearchOutcome::Found(normal));
        };
        if !enforced_ray_set.insert(violating_idx) {
            return Ok(SupportingFaceLpNormalSearchOutcome::RepeatedViolation);
        }
        enforced_ray_indices.push(violating_idx);
    }
    Ok(SupportingFaceLpNormalSearchOutcome::CuttingRoundLimit)
}

fn solve_supporting_face_normal_lp_with_enforced_rays(
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
    anchor: &[i64],
    enforced_ray_indices: &[usize],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<SupportingFaceLpNormalSearchOutcome> {
    let normal_vars = supporting_face_normal_vars(anchor.len(), options)?;
    let vars = normal_vars.variables;
    let mut objective = Expression::from(0.0);
    objective.add_mul(0.0, normal_vars.normal[0]);
    let mut model = vars.minimise(objective).using(default_solver);

    model = add_supporting_face_zero_constraints(model, &normal_vars.normal, face_generators);

    let mut anchor_expr = Expression::from(0.0);
    for (var, &coefficient) in normal_vars.normal.iter().zip(anchor) {
        if coefficient != 0 {
            anchor_expr.add_mul(coefficient as f64, *var);
        }
    }
    model = model.with(anchor_expr.eq(1.0));

    model = add_supporting_face_enforced_ray_constraints(
        model,
        &normal_vars.normal,
        mori_generators,
        enforced_ray_indices,
    )?;
    solve_supporting_face_normal_lp_model(model, &normal_vars.normal)
}

struct SupportingFaceNormalVars {
    variables: ProblemVariables,
    normal: Vec<Variable>,
}

fn supporting_face_normal_vars(
    dim: usize,
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<SupportingFaceNormalVars> {
    if dim == 0 {
        return Err(Error::InvalidInput(
            "supporting Mori face LP dimension is zero".into(),
        ));
    }
    let mut variables = ProblemVariables::new();
    let normal = (0..dim)
        .map(|_| {
            variables.add(
                variable()
                    .min(-options.variable_bound)
                    .max(options.variable_bound),
            )
        })
        .collect::<Vec<_>>();
    Ok(SupportingFaceNormalVars { variables, normal })
}

fn add_supporting_face_zero_constraints<M: SolverModel>(
    mut model: M,
    normal_vars: &[Variable],
    face_generators: &[Vec<i64>],
) -> M {
    for generator in face_generators {
        let mut expr = Expression::from(0.0);
        for (var, &coefficient) in normal_vars.iter().zip(generator) {
            if coefficient != 0 {
                expr.add_mul(coefficient as f64, *var);
            }
        }
        model = model.with(expr.eq(0.0));
    }
    model
}

fn add_supporting_face_enforced_ray_constraints<M: SolverModel>(
    mut model: M,
    normal_vars: &[Variable],
    mori_generators: &[Vec<i64>],
    enforced_ray_indices: &[usize],
) -> Result<M> {
    for &ray_idx in enforced_ray_indices {
        let ray = mori_generators.get(ray_idx).ok_or_else(|| {
            Error::InvalidInput(format!(
                "supporting Mori face enforced ray index {ray_idx} is out of bounds"
            ))
        })?;
        let mut expr = Expression::from(0.0);
        for (var, &coefficient) in normal_vars.iter().zip(ray) {
            if coefficient != 0 {
                expr.add_mul(coefficient as f64, *var);
            }
        }
        model = model.with(expr.geq(0.0));
    }
    Ok(model)
}

fn solve_supporting_face_normal_lp_model<M: SolverModel<Error = ResolutionError>>(
    model: M,
    normal_vars: &[Variable],
) -> Result<SupportingFaceLpNormalSearchOutcome> {
    let solution = match model.solve() {
        Ok(solution) => solution,
        Err(ResolutionError::Infeasible | ResolutionError::Unbounded) => {
            return Ok(SupportingFaceLpNormalSearchOutcome::LpNoSolution);
        }
        Err(ResolutionError::Other("NoSolutionFound")) => {
            return Ok(SupportingFaceLpNormalSearchOutcome::LpNoSolution);
        }
        Err(ResolutionError::Other(message)) => {
            return Ok(SupportingFaceLpNormalSearchOutcome::SolverOther(
                message.to_string(),
            ));
        }
        Err(ResolutionError::Str(message)) => {
            return Ok(SupportingFaceLpNormalSearchOutcome::SolverOther(message));
        }
    };
    let normal = normal_vars
        .iter()
        .map(|var| solution.value(*var))
        .collect::<Vec<_>>();
    if normal.iter().all(|value| value.is_finite()) {
        Ok(SupportingFaceLpNormalSearchOutcome::Found(normal))
    } else {
        Err(Error::InvalidInput(
            "supporting Mori face normal LP returned a non-finite value".into(),
        ))
    }
}

fn most_negative_lp_normal_violation(
    lp_normal: &[f64],
    mori_generators: &[Vec<i64>],
) -> Option<usize> {
    let mut worst_idx = None;
    let mut worst_dot = -1.0e-7;
    for (idx, ray) in mori_generators.iter().enumerate() {
        let dot = lp_normal
            .iter()
            .zip(ray)
            .map(|(&normal_coeff, &ray_coeff)| normal_coeff * ray_coeff as f64)
            .sum::<f64>();
        if dot < worst_dot {
            worst_dot = dot;
            worst_idx = Some(idx);
        }
    }
    worst_idx
}

fn integer_supporting_face_certificate_from_lp(
    lp_normal: &[f64],
    face_generators: &[Vec<i64>],
    mori_generators: &[Vec<i64>],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<Option<SupportingMoriFaceCertificate>> {
    let mut seen_normals = HashSet::new();
    for scale in 1..=options.scale_limit {
        let Some(normal) = rounded_reduced_i64_normal(lp_normal, scale)? else {
            continue;
        };
        if !seen_normals.insert(normal.clone()) {
            continue;
        }
        if let Some(certificate) =
            check_supporting_mori_face_normal(&normal, face_generators, mori_generators)?
        {
            return Ok(Some(certificate));
        }
    }
    Ok(None)
}

fn rounded_reduced_i64_normal(lp_normal: &[f64], scale: i64) -> Result<Option<Vec<i64>>> {
    let mut normal = Vec::with_capacity(lp_normal.len());
    for &value in lp_normal {
        let scaled = value * scale as f64;
        if !scaled.is_finite() || scaled < i64::MIN as f64 || scaled > i64::MAX as f64 {
            return Err(Error::InvalidInput(
                "supporting Mori face LP normal does not fit in i64 after scaling".into(),
            ));
        }
        normal.push(scaled.round() as i64);
    }
    reduce_i64_vector_preserve_sign(&normal)
}

fn reduce_i64_vector_preserve_sign(values: &[i64]) -> Result<Option<Vec<i64>>> {
    let mut gcd = 0i64;
    for &value in values {
        if value == i64::MIN {
            return Err(Error::InvalidInput(
                "supporting Mori face normal coefficient is i64::MIN".into(),
            ));
        }
        gcd = gcd_i64(gcd, value.abs());
    }
    if gcd == 0 {
        return Ok(None);
    }
    values
        .iter()
        .map(|&value| {
            value.checked_div(gcd).ok_or_else(|| {
                Error::InvalidInput("supporting Mori face normal reduction divided by zero".into())
            })
        })
        .collect::<Result<Vec<_>>>()
        .map(Some)
}

/// Extract the Mori generators cut out by an exact supporting normal.
///
/// The returned generators are the rows of `mori_generators` whose exact
/// integer pairing with `normal` is zero. If the proposed normal is not
/// supporting, this returns `Ok(None)`.
pub fn supporting_mori_face_from_normal(
    normal: &[i64],
    mori_generators: &[Vec<i64>],
) -> Result<Option<SupportingMoriFace>> {
    let Some(certificate) = check_supporting_mori_face_normal(normal, &[], mori_generators)? else {
        return Ok(None);
    };

    let mut generators = Vec::with_capacity(certificate.zero_generator_count);
    for generator in mori_generators {
        let dot = exact_i64_dot_checked(normal, generator)?;
        if dot == 0 {
            generators.push(generator.clone());
        }
    }
    debug_assert_eq!(generators.len(), certificate.zero_generator_count);

    Ok(Some(SupportingMoriFace {
        certificate,
        generators,
    }))
}

/// Extract a supporting Mori face only if it contains `target_curve` rationally.
///
/// This is the exact validation boundary needed before using a low-dimensional
/// face as the HKTY context for a target ray. The normal must be supporting, and
/// `target_curve` must lie in the rational row span of the zero-pairing Mori
/// generators.
pub fn supporting_mori_face_for_curve_from_normal(
    normal: &[i64],
    target_curve: &[i64],
    mori_generators: &[Vec<i64>],
) -> Result<Option<SupportingMoriFace>> {
    let Some(face) = supporting_mori_face_from_normal(normal, mori_generators)? else {
        return Ok(None);
    };
    if curve_in_rational_row_span(target_curve, &face.generators)? {
        Ok(Some(face))
    } else {
        Ok(None)
    }
}

/// Check an integer separator as an exact extremal-ray certificate.
///
/// This verifies a Farkas-style certificate for the finite cone generated by
/// `mori_generators`: after orienting `separator_normal`, the target curve must
/// pair negatively with it, while every Mori generator not on the same positive
/// rational ray pairs non-negatively. If this succeeds, `target_curve` cannot
/// be a non-negative real combination of the other generator rays.
///
/// This is only the cone-theoretic ray certificate. It does not compute a GV
/// invariant and does not identify a shrinking divisor.
pub fn check_extremal_mori_ray_separator(
    separator_normal: &[i64],
    target_curve: &[i64],
    mori_generators: &[Vec<i64>],
) -> Result<Option<ExtremalMoriRayCertificate>> {
    if separator_normal.is_empty() {
        return Err(Error::InvalidInput(
            "extremal Mori ray separator normal is empty".into(),
        ));
    }
    if target_curve.len() != separator_normal.len() {
        return Err(Error::InvalidInput(format!(
            "target curve dimension {} does not match separator normal dimension {}",
            target_curve.len(),
            separator_normal.len()
        )));
    }
    if mori_generators.is_empty() {
        return Err(Error::InvalidInput(
            "extremal Mori ray check requires Mori generators".into(),
        ));
    }

    let target_primitive =
        primitive_ray_from_curve_class(target_curve, "extremal Mori ray target")?;
    if let Some(certificate) = check_extremal_mori_ray_separator_oriented(
        separator_normal,
        target_curve,
        &target_primitive,
        mori_generators,
    )? {
        return Ok(Some(certificate));
    }

    let opposite_normal = checked_neg_i64_vector(separator_normal, "extremal Mori ray separator")?;
    check_extremal_mori_ray_separator_oriented(
        &opposite_normal,
        target_curve,
        &target_primitive,
        mori_generators,
    )
}

/// Find an exact integer separator for a target extremal Mori ray.
///
/// This constructs the exact cone of candidate separators
/// `{n | n.other >= 0, n.target <= 0}` using all Mori generators not on the
/// same positive rational ray as `target_curve`, enumerates its rays via the
/// double-description method, and returns the first ray with
/// `n.target < 0` that passes [`check_extremal_mori_ray_separator`].
///
/// A successful result certifies extremality in the supplied finite
/// generator cone. It does not certify that the generator set is the complete
/// Mori cone, and it does not compute GV invariants or chamber data.
pub fn find_extremal_mori_ray_separator(
    target_curve: &[i64],
    mori_generators: &[Vec<i64>],
) -> Result<Option<ExtremalMoriRayCertificate>> {
    if target_curve.is_empty() {
        return Err(Error::InvalidInput(
            "extremal Mori ray target is empty".into(),
        ));
    }
    if mori_generators.is_empty() {
        return Err(Error::InvalidInput(
            "extremal Mori ray separator search requires Mori generators".into(),
        ));
    }

    let target_primitive =
        primitive_ray_from_curve_class(target_curve, "extremal Mori ray target")?;
    let mut hyperplanes = Vec::with_capacity(mori_generators.len() + 1);
    for generator in mori_generators {
        validate_curve_dimension("Mori generator", generator, target_curve.len())?;
        if same_positive_rational_ray(&target_primitive, generator)? {
            continue;
        }
        hyperplanes.push(
            generator
                .iter()
                .map(|&entry| i128::from(entry))
                .collect::<Vec<_>>(),
        );
    }
    hyperplanes.push(
        target_curve
            .iter()
            .map(|&entry| -i128::from(entry))
            .collect::<Vec<_>>(),
    );

    let mut separator_cone = Cone::from_hyperplanes(hyperplanes);
    let separator_rays = separator_cone.rays().to_vec();
    for ray in separator_rays {
        let separator = i128_vector_to_i64(&ray, "extremal Mori ray separator")?;
        if checked_i128_dot(&separator, target_curve, "extremal Mori ray separator")? >= 0 {
            continue;
        }
        if let Some(certificate) =
            check_extremal_mori_ray_separator(&separator, target_curve, mori_generators)?
        {
            return Ok(Some(certificate));
        }
    }

    Ok(None)
}

/// Find an exact integer separator for a target extremal Mori ray using an LP
/// candidate normal followed by exact verification.
///
/// This is a cheaper companion to [`find_extremal_mori_ray_separator`]. The LP
/// solve only proposes a real normal satisfying `n.target <= -1` and
/// `n.other >= 0`; any returned certificate has still passed
/// [`check_extremal_mori_ray_separator`] with integer arithmetic. A failed LP
/// search is inconclusive and must not be read as proof that no separator
/// exists.
pub fn find_extremal_mori_ray_separator_by_lp_search(
    target_curve: &[i64],
    mori_generators: &[Vec<i64>],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<Option<ExtremalMoriRayCertificate>> {
    Ok(
        diagnose_extremal_mori_ray_separator_by_lp_search(target_curve, mori_generators, options)?
            .certificate,
    )
}

/// Diagnose LP-assisted extremal Mori ray separator search.
///
/// This exposes whether the LP failed to find a real normal, found one that
/// could not be converted to an exact integer certificate within the scale
/// limit, or produced a fully verified separator.
pub fn diagnose_extremal_mori_ray_separator_by_lp_search(
    target_curve: &[i64],
    mori_generators: &[Vec<i64>],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<ExtremalMoriRayLpSearchDiagnostic> {
    validate_extremal_mori_ray_lp_inputs(target_curve, mori_generators, options)?;

    let target_primitive =
        primitive_ray_from_curve_class(target_curve, "extremal Mori ray target")?;
    let normal_vars = supporting_face_normal_vars(target_curve.len(), options)?;
    let vars = normal_vars.variables;
    let mut objective = Expression::from(0.0);
    objective.add_mul(0.0, normal_vars.normal[0]);
    let mut model = vars.minimise(objective).using(default_solver);

    let mut target_expr = Expression::from(0.0);
    for (var, &coefficient) in normal_vars.normal.iter().zip(target_curve) {
        if coefficient != 0 {
            target_expr.add_mul(coefficient as f64, *var);
        }
    }
    model = model.with(target_expr.leq(-1.0));

    for generator in mori_generators {
        if same_positive_rational_ray(&target_primitive, generator)? {
            continue;
        }
        let mut expr = Expression::from(0.0);
        for (var, &coefficient) in normal_vars.normal.iter().zip(generator) {
            if coefficient != 0 {
                expr.add_mul(coefficient as f64, *var);
            }
        }
        model = model.with(expr.geq(0.0));
    }

    match solve_supporting_face_normal_lp_model(model, &normal_vars.normal)? {
        SupportingFaceLpNormalSearchOutcome::Found(normal) => {
            integer_extremal_ray_certificate_from_lp_diagnostic(
                &normal,
                target_curve,
                mori_generators,
                options,
            )
        }
        other => Ok(ExtremalMoriRayLpSearchDiagnostic {
            status: other.status(),
            lp_solution_found: false,
            exact_normal_candidate_count: 0,
            certificate: None,
        }),
    }
}

fn validate_extremal_mori_ray_lp_inputs(
    target_curve: &[i64],
    mori_generators: &[Vec<i64>],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<()> {
    if target_curve.is_empty() {
        return Err(Error::InvalidInput(
            "extremal Mori ray LP target is empty".into(),
        ));
    }
    if mori_generators.is_empty() {
        return Err(Error::InvalidInput(
            "extremal Mori ray LP search requires Mori generators".into(),
        ));
    }
    for generator in mori_generators {
        validate_curve_dimension("Mori generator", generator, target_curve.len())?;
    }
    if options.scale_limit <= 0 {
        return Err(Error::InvalidInput(
            "extremal Mori ray LP search scale limit must be positive".into(),
        ));
    }
    if !options.variable_bound.is_finite() || options.variable_bound <= 0.0 {
        return Err(Error::InvalidInput(
            "extremal Mori ray LP search variable bound must be positive and finite".into(),
        ));
    }
    Ok(())
}

fn integer_extremal_ray_certificate_from_lp_diagnostic(
    lp_normal: &[f64],
    target_curve: &[i64],
    mori_generators: &[Vec<i64>],
    options: &SupportingMoriFaceLpSearchOptions,
) -> Result<ExtremalMoriRayLpSearchDiagnostic> {
    let mut seen_normals = HashSet::new();
    for scale in 1..=options.scale_limit {
        let Some(normal) = rounded_reduced_i64_normal(lp_normal, scale)? else {
            continue;
        };
        if !seen_normals.insert(normal.clone()) {
            continue;
        }
        if let Some(certificate) =
            check_extremal_mori_ray_separator(&normal, target_curve, mori_generators)?
        {
            return Ok(ExtremalMoriRayLpSearchDiagnostic {
                status: "certified_lp_separator".to_string(),
                lp_solution_found: true,
                exact_normal_candidate_count: seen_normals.len(),
                certificate: Some(certificate),
            });
        }
    }
    Ok(ExtremalMoriRayLpSearchDiagnostic {
        status: "lp_solution_no_exact_integer_separator".to_string(),
        lp_solution_found: true,
        exact_normal_candidate_count: seen_normals.len(),
        certificate: None,
    })
}

fn check_extremal_mori_ray_separator_oriented(
    separator_normal: &[i64],
    target_curve: &[i64],
    target_primitive: &[i64],
    mori_generators: &[Vec<i64>],
) -> Result<Option<ExtremalMoriRayCertificate>> {
    if separator_normal.iter().all(|&value| value == 0) {
        return Ok(None);
    }
    if exact_i64_dot_checked(separator_normal, target_curve)? >= 0 {
        return Ok(None);
    }

    let mut same_ray_generator_count = 0usize;
    let mut zero_other_generator_count = 0usize;
    let mut positive_other_generator_count = 0usize;
    for generator in mori_generators {
        validate_curve_dimension("Mori generator", generator, separator_normal.len())?;
        if same_positive_rational_ray(target_primitive, generator)? {
            same_ray_generator_count += 1;
            continue;
        }
        let dot = exact_i64_dot_checked(separator_normal, generator)?;
        if dot < 0 {
            return Ok(None);
        }
        if dot == 0 {
            zero_other_generator_count += 1;
        } else {
            positive_other_generator_count += 1;
        }
    }

    if same_ray_generator_count == 0 {
        return Ok(None);
    }

    Ok(Some(ExtremalMoriRayCertificate {
        separator_normal: separator_normal.to_vec(),
        same_ray_generator_count,
        zero_other_generator_count,
        positive_other_generator_count,
    }))
}

fn same_positive_rational_ray(target_primitive: &[i64], generator: &[i64]) -> Result<bool> {
    if generator.iter().all(|&entry| entry == 0) {
        return Ok(false);
    }
    Ok(primitive_ray_from_curve_class(generator, "Mori generator")? == target_primitive)
}

/// Compute genus-zero GV invariants along a ray inside a supplied generator context.
///
/// The ray and `provided_generators` are expressed in Kähler-basis curve
/// coordinates. The generators are passed to the CYTools-style
/// `mcap_generators` path as the caller-supplied local face/semigroup context,
/// and the values for `q, 2q, ..., max_multiple*q` are extracted from the cygv
/// output. Missing multiples are returned as zero.
///
/// # Errors
/// Returns an error for invalid dimensions, empty generator context,
/// non-positive grading degree, integer overflow, or any cygv failure. In
/// panic-unwind builds, cygv panics are converted into errors rather than being
/// hidden.
pub fn compute_ray_gv_series_with_provided_generators(
    ray: &[i64],
    provided_generators: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    max_multiple: u32,
) -> Result<OneDimensionalRayGvSeries> {
    if cfg!(panic = "abort") {
        return Err(Error::InvalidInput(
            "one-dimensional ray GV series requires a panic=unwind build because upstream cygv can still panic internally".into(),
        ));
    }
    if max_multiple == 0 {
        return Err(Error::InvalidInput(
            "one-dimensional ray GV series requires at least one multiple".into(),
        ));
    }
    if ray.is_empty() {
        return Err(Error::InvalidInput(
            "provided-generator ray GV series ray is empty".into(),
        ));
    }
    if provided_generators.is_empty() {
        return Err(Error::InvalidInput(
            "provided-generator ray GV series requires at least one generator".into(),
        ));
    }
    if grading_vector.len() != ray.len() {
        return Err(Error::InvalidInput(
            "provided-generator ray GV grading dimension does not match ray dimension".into(),
        ));
    }
    for generator in provided_generators {
        if generator.len() != ray.len() {
            return Err(Error::InvalidInput(
                "provided-generator ray GV generator dimension does not match ray dimension".into(),
            ));
        }
    }

    let degree = ray
        .iter()
        .zip(grading_vector.iter())
        .map(|(&coefficient, &weight)| i128::from(coefficient) * i128::from(weight))
        .sum::<i128>();
    if degree <= 0 {
        return Err(Error::InvalidInput(format!(
            "provided-generator ray GV target has non-positive grading degree {degree}"
        )));
    }
    let max_degree = degree
        .checked_mul(i128::from(max_multiple))
        .ok_or_else(|| Error::InvalidInput("one-dimensional ray GV degree overflow".into()))?;
    let max_degree_u32 = u32::try_from(max_degree).map_err(|_| {
        Error::InvalidInput("one-dimensional ray GV max degree does not fit in u32".into())
    })?;

    let previous_panic_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let gvs_result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        compute_gv_invariants_with_provided_generators(
            provided_generators,
            grading_vector,
            q_matrix,
            intnums,
            None,
            Some(max_degree_u32),
        )
    }));
    std::panic::set_hook(previous_panic_hook);

    let gvs = match gvs_result {
        Ok(Ok(gvs)) => gvs,
        Ok(Err(err)) => return Err(err),
        Err(payload) => {
            return Err(Error::InvalidInput(format!(
                "one-dimensional ray GV computation panicked: {}",
                panic_payload_message(payload.as_ref())
            )));
        }
    };

    let mut by_class: HashMap<Vec<i32>, Integer> = gvs.into_iter().collect();
    let mut values = Vec::with_capacity(max_multiple as usize);
    for multiple in 1..=max_multiple {
        let target = ray
            .iter()
            .map(|&value| {
                let scaled = i128::from(value) * i128::from(multiple);
                i32::try_from(scaled).map_err(|_| {
                    Error::InvalidInput(
                        "one-dimensional ray GV target coordinate does not fit in i32".into(),
                    )
                })
            })
            .collect::<Result<Vec<_>>>()?;
        values.push(by_class.remove(&target).unwrap_or_else(|| Integer::from(0)));
    }

    Ok(OneDimensionalRayGvSeries {
        ray: ray.to_vec(),
        degree,
        values,
    })
}

/// Compute genus-zero GV invariants along a one-dimensional ray.
///
/// The ray is expressed in Kähler-basis curve coordinates. The function calls
/// the CYTools-style `mcap_generators` path with this single generator and
/// extracts the values for `q, 2q, ..., max_multiple*q`. Missing multiples are
/// returned as zero.
///
/// # Errors
/// Returns an error for invalid dimensions, non-positive grading degree,
/// integer overflow, or any cygv failure. In panic-unwind builds, cygv panics
/// are converted into errors rather than being hidden.
pub fn compute_one_dimensional_ray_gv_series(
    ray: &[i64],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    max_multiple: u32,
) -> Result<OneDimensionalRayGvSeries> {
    compute_ray_gv_series_with_provided_generators(
        ray,
        &[ray.to_vec()],
        grading_vector,
        q_matrix,
        intnums,
        max_multiple,
    )
}

/// Compute genus-zero GV invariants along an ambient one-dimensional ray.
///
/// This is a convenience boundary for McAllister-style ambient curve files. It
/// projects the ambient class to the supplied divisor basis, then delegates to
/// [`compute_one_dimensional_ray_gv_series`]. It does not change the cygv
/// semigroup semantics: a one-generator call remains a local diagnostic unless
/// the caller has established that this is the intended face context.
pub fn compute_ambient_one_dimensional_ray_gv_series(
    ambient_ray: &[i64],
    basis: &[usize],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    max_multiple: u32,
) -> Result<OneDimensionalRayGvSeries> {
    let ray = project_ambient_curve_to_basis(ambient_ray, basis)?;
    compute_one_dimensional_ray_gv_series(&ray, grading_vector, q_matrix, intnums, max_multiple)
}

/// Compute genus-zero GV invariants for the local `P^2` geometry.
///
/// This is the one-parameter local mirror-symmetry computation for
/// `O(-3) -> P^2`, whose toric charge is `(-3, 1, 1, 1)`. It is not a lookup
/// table: Cyrus expands the Picard-Fuchs mirror map, transforms the local
/// B-model Yukawa coupling to the flat coordinate, and applies the standard
/// multiple-cover inversion
/// `K(Q) = -1/3 + sum_d N_d d^3 Q^d / (1 - Q^d)`.
pub fn compute_local_p2_genus_zero_gv_series(max_degree: usize) -> Result<Vec<Integer>> {
    if max_degree == 0 {
        return Ok(Vec::new());
    }

    let mirror_correction = local_p2_mirror_correction(max_degree);
    let z_of_q = local_p2_inverse_mirror_map(&mirror_correction, max_degree);
    let yukawa_z = local_p2_log_yukawa_in_b_model_coordinate(&mirror_correction, max_degree)?;
    let yukawa_q = rational_series_compose(&yukawa_z, &z_of_q, max_degree);

    extract_genus_zero_gv_from_yukawa(&yukawa_q, max_degree)
}

/// Compute genus-zero GV invariants for a recognized affine toric circuit.
///
/// Currently the only supported local model is the rank-two local `P^2`
/// triangle. Other circuit shapes return `Ok(None)` so callers cannot silently
/// promote an unsupported local model to a zero or fitted GV sequence.
pub fn compute_local_toric_circuit_gv_series(
    kind: &LocalToricCircuitKind,
    max_degree: usize,
) -> Result<Option<Vec<Integer>>> {
    match kind {
        LocalToricCircuitKind::LocalP2Triangle { .. } => {
            compute_local_p2_genus_zero_gv_series(max_degree).map(Some)
        }
    }
}

fn local_p2_mirror_correction(order: usize) -> Vec<Rational> {
    let mut correction = vec![Rational::from(0); order + 1];
    for degree in 1..=order {
        let mut numerator = factorial_integer(3 * degree);
        if degree % 2 == 1 {
            numerator = -numerator;
        }
        let degree_factorial = factorial_integer(degree);
        let mut denominator = degree_factorial.clone() * degree_factorial.clone();
        denominator *= degree_factorial;
        denominator *= Integer::from(degree);
        correction[degree] = Rational::from(numerator) / Rational::from(denominator);
    }
    correction
}

fn local_p2_inverse_mirror_map(mirror_correction: &[Rational], order: usize) -> Vec<Rational> {
    let mut z_of_q = vec![Rational::from(0); order + 1];
    z_of_q[1] = Rational::from(1);

    for _ in 0..order {
        let correction_at_z = rational_series_compose(mirror_correction, &z_of_q, order);
        let exp_correction = rational_series_exp(&correction_at_z, order);
        let inverse_exp_correction = rational_series_inverse(&exp_correction, order)
            .expect("formal exponential has unit constant coefficient");
        let mut next = vec![Rational::from(0); order + 1];
        next[1..=order].clone_from_slice(&inverse_exp_correction[..order]);
        z_of_q = next;
    }

    z_of_q
}

fn local_p2_log_yukawa_in_b_model_coordinate(
    mirror_correction: &[Rational],
    order: usize,
) -> Result<Vec<Rational>> {
    let mut theta_t = vec![Rational::from(0); order + 1];
    theta_t[0] = Rational::from(1);
    for degree in 1..=order {
        theta_t[degree] = Rational::from(degree) * mirror_correction[degree].clone();
    }

    let theta_t_cubed = rational_series_pow(&theta_t, 3, order);
    let mut discriminant = vec![Rational::from(0); order + 1];
    discriminant[0] = Rational::from(1);
    if order >= 1 {
        discriminant[1] = Rational::from(27);
    }
    let denominator = rational_series_mul(&discriminant, &theta_t_cubed, order);
    let mut yukawa = rational_series_inverse(&denominator, order).ok_or_else(|| {
        Error::InvalidInput("local P2 Yukawa denominator has zero constant term".into())
    })?;
    for coefficient in &mut yukawa {
        *coefficient *= Rational::from_signeds(-1, 3);
    }
    Ok(yukawa)
}

fn extract_genus_zero_gv_from_yukawa(
    yukawa: &[Rational],
    max_degree: usize,
) -> Result<Vec<Integer>> {
    if yukawa.len() <= max_degree {
        return Err(Error::InvalidInput(
            "local GV Yukawa series is shorter than requested degree".into(),
        ));
    }

    let mut gvs = Vec::with_capacity(max_degree);
    for degree in 1..=max_degree {
        let mut lower_multiple_cover_sum = Rational::from(0);
        for divisor in 1..degree {
            if degree % divisor != 0 {
                continue;
            }
            let divisor_cubed = divisor * divisor * divisor;
            lower_multiple_cover_sum +=
                Rational::from(&gvs[divisor - 1]) * Rational::from(divisor_cubed);
        }
        let gv_rational =
            (yukawa[degree].clone() - lower_multiple_cover_sum) / Rational::from(degree.pow(3));
        if gv_rational.denominator_ref() != &1u32 {
            return Err(Error::InvalidInput(format!(
                "local P2 GV invariant at degree {degree} is not integral: {gv_rational}"
            )));
        }
        let gv = Integer::try_from(gv_rational).map_err(|_| {
            Error::InvalidInput(format!(
                "local P2 GV invariant at degree {degree} is not integral"
            ))
        })?;
        gvs.push(gv);
    }
    Ok(gvs)
}

fn factorial_integer(n: usize) -> Integer {
    let mut out = Integer::from(1);
    for factor in 2..=n {
        out *= Integer::from(factor);
    }
    out
}

fn rational_series_mul(lhs: &[Rational], rhs: &[Rational], order: usize) -> Vec<Rational> {
    let mut out = vec![Rational::from(0); order + 1];
    for (lhs_degree, lhs_coefficient) in lhs.iter().enumerate().take(order + 1) {
        if *lhs_coefficient == 0 {
            continue;
        }
        for (rhs_degree, rhs_coefficient) in rhs.iter().enumerate().take(order + 1 - lhs_degree) {
            if *rhs_coefficient == 0 {
                continue;
            }
            out[lhs_degree + rhs_degree] += lhs_coefficient * rhs_coefficient;
        }
    }
    out
}

fn rational_series_pow(series: &[Rational], exponent: usize, order: usize) -> Vec<Rational> {
    let mut out = vec![Rational::from(0); order + 1];
    out[0] = Rational::from(1);
    for _ in 0..exponent {
        out = rational_series_mul(&out, series, order);
    }
    out
}

fn rational_series_exp(series: &[Rational], order: usize) -> Vec<Rational> {
    let mut out = vec![Rational::from(0); order + 1];
    out[0] = Rational::from(1);
    for degree in 1..=order {
        let mut coefficient = Rational::from(0);
        for term_degree in 1..=degree {
            if term_degree >= series.len() {
                break;
            }
            coefficient += Rational::from(term_degree)
                * series[term_degree].clone()
                * out[degree - term_degree].clone();
        }
        out[degree] = coefficient / Rational::from(degree);
    }
    out
}

fn rational_series_inverse(series: &[Rational], order: usize) -> Option<Vec<Rational>> {
    let constant = series.first()?;
    if *constant == 0 {
        return None;
    }
    let mut out = vec![Rational::from(0); order + 1];
    out[0] = Rational::from(1) / constant;
    for degree in 1..=order {
        let mut coefficient = Rational::from(0);
        for term_degree in 1..=degree {
            if term_degree >= series.len() {
                break;
            }
            coefficient += series[term_degree].clone() * out[degree - term_degree].clone();
        }
        out[degree] = -coefficient / constant;
    }
    Some(out)
}

fn rational_series_compose(
    series: &[Rational],
    argument: &[Rational],
    order: usize,
) -> Vec<Rational> {
    let mut out = vec![Rational::from(0); order + 1];
    let mut argument_power = vec![Rational::from(0); order + 1];
    argument_power[0] = Rational::from(1);
    for degree in 0..=order {
        if degree > 0 {
            argument_power = rational_series_mul(&argument_power, argument, order);
        }
        let coefficient = series
            .get(degree)
            .cloned()
            .unwrap_or_else(|| Rational::from(0));
        if coefficient == 0 {
            continue;
        }
        for term_degree in 0..=order {
            out[term_degree] += coefficient.clone() * argument_power[term_degree].clone();
        }
    }
    out
}

fn integer_abs_ln(value: &Integer) -> Result<Option<f64>> {
    let magnitude = value.clone().abs();
    if magnitude == 0 {
        return Ok(None);
    }

    let decimal = magnitude.to_string();
    if decimal.len() < 300 {
        let parsed = decimal.parse::<f64>().map_err(|err| {
            Error::InvalidInput(format!("failed to parse GV magnitude as f64: {err}"))
        })?;
        return Ok(Some(parsed.ln()));
    }

    let leading_len = decimal.len().min(16);
    let leading = decimal[..leading_len].parse::<f64>().map_err(|err| {
        Error::InvalidInput(format!(
            "failed to parse leading GV magnitude digits: {err}"
        ))
    })?;
    Ok(Some(
        leading.ln() + (decimal.len() - leading_len) as f64 * LN_10,
    ))
}

fn panic_payload_message(payload: &(dyn std::any::Any + Send)) -> String {
    if let Some(message) = payload.downcast_ref::<&str>() {
        (*message).to_string()
    } else if let Some(message) = payload.downcast_ref::<String>() {
        message.clone()
    } else {
        "non-string panic payload".to_string()
    }
}

fn validate_curve_dimension(label: &str, curve: &[i64], expected_dim: usize) -> Result<()> {
    if curve.len() == expected_dim {
        Ok(())
    } else {
        Err(Error::InvalidInput(format!(
            "{label} dimension {} does not match normal dimension {expected_dim}",
            curve.len()
        )))
    }
}

fn exact_i64_dot_checked(lhs: &[i64], rhs: &[i64]) -> Result<i128> {
    lhs.iter()
        .zip(rhs)
        .try_fold(0_i128, |acc, (&left, &right)| {
            let product = i128::from(left)
                .checked_mul(i128::from(right))
                .ok_or_else(|| {
                    Error::InvalidInput("supporting Mori face dot product overflowed".into())
                })?;
            acc.checked_add(product).ok_or_else(|| {
                Error::InvalidInput("supporting Mori face dot product overflowed".into())
            })
        })
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

/// Find an exact finite-semigroup decomposition of `target` into other selected candidates.
///
/// This solves integer feasibility over the finite set of supplied candidates,
/// excluding rows with the same class as `target`. It is stricter than the
/// pair-only pruning helper: a returned decomposition can contain any
/// non-negative integer multiplicities allowed by the target volume. The LP
/// solver works over floating point constraints, so the candidate solution is
/// always verified exactly in integer arithmetic before being returned.
pub fn find_semigroup_decomposition(
    target: &ToricCurveCandidate,
    candidates: &[ToricCurveCandidate],
) -> Result<Option<Vec<CurveDecompositionTerm>>> {
    let dim = validate_semigroup_inputs(target, candidates)?;
    let mut vars = ProblemVariables::new();
    let solver_vars = semigroup_solver_variables(target, candidates, dim, &mut vars)?;

    if solver_vars.is_empty() {
        return Ok(None);
    }

    solve_semigroup_model(vars, &solver_vars, target, candidates, dim)
}

fn validate_semigroup_inputs(
    target: &ToricCurveCandidate,
    candidates: &[ToricCurveCandidate],
) -> Result<usize> {
    let dim = target.class.len();
    if dim == 0 {
        return Err(Error::InvalidInput("target curve dimension is zero".into()));
    }
    validate_exact_lp_class(&target.class)?;

    for candidate in candidates {
        if candidate.class.len() != dim {
            return Err(Error::InvalidInput(
                "candidate curve dimensions are inconsistent".into(),
            ));
        }
        validate_exact_lp_class(&candidate.class)?;
    }

    Ok(dim)
}

fn semigroup_solver_variables(
    target: &ToricCurveCandidate,
    candidates: &[ToricCurveCandidate],
    dim: usize,
    vars: &mut ProblemVariables,
) -> Result<Vec<(usize, Variable)>> {
    let mut solver_vars = Vec::new();
    let target_volume = target.volume.get();
    for (idx, candidate) in candidates.iter().enumerate() {
        debug_assert_eq!(candidate.class.len(), dim);
        if candidate.class == target.class || candidate.volume >= target.volume {
            continue;
        }
        let candidate_volume = candidate.volume.get();
        let upper_bound = (target_volume / candidate_volume).ceil() + 1.0;
        if upper_bound < 1.0 || !upper_bound.is_finite() {
            return Err(Error::InvalidInput(
                "invalid curve volume bound for semigroup decomposition".into(),
            ));
        }
        let variable = vars.add(variable().integer().min(0.0).max(upper_bound));
        solver_vars.push((idx, variable));
    }
    Ok(solver_vars)
}

fn semigroup_total_expression(solver_vars: &[(usize, Variable)]) -> Expression {
    let mut expr = Expression::from(0.0);
    for &(_, variable) in solver_vars {
        expr.add_mul(1.0, variable);
    }
    expr
}

fn semigroup_coordinate_expression(
    coord: usize,
    solver_vars: &[(usize, Variable)],
    candidates: &[ToricCurveCandidate],
) -> Expression {
    let mut expr = Expression::from(0.0);
    for &(candidate_idx, variable) in solver_vars {
        let coefficient = candidates[candidate_idx].class[coord];
        if coefficient != 0 {
            expr.add_mul(coefficient as f64, variable);
        }
    }
    expr
}

fn solve_semigroup_model(
    vars: ProblemVariables,
    solver_vars: &[(usize, Variable)],
    target: &ToricCurveCandidate,
    candidates: &[ToricCurveCandidate],
    dim: usize,
) -> Result<Option<Vec<CurveDecompositionTerm>>> {
    let mut model = vars
        .minimise(semigroup_total_expression(solver_vars))
        .using(default_solver)
        .with(semigroup_total_expression(solver_vars).geq(2.0));

    for coord in 0..dim {
        let expr = semigroup_coordinate_expression(coord, solver_vars, candidates);
        model = model.with(expr.eq(target.class[coord] as f64));
    }

    let solution = match model.solve() {
        Ok(solution) => solution,
        Err(ResolutionError::Infeasible) => return Ok(None),
        Err(err) => {
            return Err(Error::InvalidInput(format!(
                "curve semigroup decomposition solver failed: {err}"
            )));
        }
    };

    let terms = collect_semigroup_solution_terms(&solution, solver_vars, candidates)?;
    verify_semigroup_decomposition(target, &terms)?;
    if terms.is_empty() {
        Ok(None)
    } else {
        Ok(Some(terms))
    }
}

fn collect_semigroup_solution_terms<S: Solution>(
    solution: &S,
    solver_vars: &[(usize, Variable)],
    candidates: &[ToricCurveCandidate],
) -> Result<Vec<CurveDecompositionTerm>> {
    let mut terms = Vec::new();
    for &(candidate_idx, variable) in solver_vars {
        let value = solution.value(variable);
        if !value.is_finite() {
            return Err(Error::InvalidInput(
                "curve semigroup solver returned a non-finite multiplicity".into(),
            ));
        }
        let rounded = value.round();
        if (value - rounded).abs() > 1.0e-7 {
            return Err(Error::InvalidInput(format!(
                "curve semigroup solver returned non-integral multiplicity {value}"
            )));
        }
        if rounded < 0.0 || rounded > u64::MAX as f64 {
            return Err(Error::InvalidInput(format!(
                "curve semigroup solver returned invalid multiplicity {value}"
            )));
        }
        let multiplicity = rounded as u64;
        if multiplicity > 0 {
            terms.push(CurveDecompositionTerm {
                class: candidates[candidate_idx].class.clone(),
                multiplicity,
            });
        }
    }
    terms.sort_unstable_by(|lhs, rhs| lhs.class.cmp(&rhs.class));
    Ok(terms)
}

/// Remove curve candidates that are finite-semigroup sums of other selected candidates.
///
/// This is the faithful finite-set version of the McAllister "sums of others"
/// pruning rule. It is suitable for selected small-curve sets, not for proving a
/// Hilbert basis of the full Mori cone when the candidate set is incomplete.
pub fn remove_semigroup_decomposable_curve_candidates(
    candidates: &[ToricCurveCandidate],
) -> Result<Vec<ToricCurveCandidate>> {
    let mut out = Vec::with_capacity(candidates.len());
    for candidate in candidates {
        if find_semigroup_decomposition(candidate, candidates)?.is_none() {
            out.push(candidate.clone());
        }
    }
    Ok(out)
}

/// Prune selected toric curve candidates with an explicit pruning strategy.
pub fn prune_decomposable_curve_candidates(
    candidates: &[ToricCurveCandidate],
    strategy: CurvePruningStrategy,
) -> Result<Vec<ToricCurveCandidate>> {
    match strategy {
        CurvePruningStrategy::PairDecomposable => {
            remove_pair_decomposable_curve_candidates(candidates)
        }
        CurvePruningStrategy::FiniteSemigroup => {
            remove_semigroup_decomposable_curve_candidates(candidates)
        }
    }
}

/// Exact bounded semigroup-decomposition diagnostic for toric curve candidates.
///
/// This is intentionally diagnostic-only. It can prove that a curve is a sum
/// of up to four other selected candidates, which is useful for investigating
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
    /// Currently supports `max_terms <= 4`. Terms may repeat, matching
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
        if max_terms > 4 {
            return Err(Error::InvalidInput(
                "bounded curve decomposition currently supports max_terms <= 4".into(),
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

        if max_terms < 4 {
            return Ok(None);
        }

        for (first_pair_sum, first_pairs) in &self.pair_sums {
            let remainder = subtract_curve(&target.class, first_pair_sum);
            let Some(second_pairs) = self.pair_sums.get(&remainder) else {
                continue;
            };
            for &(i, j) in first_pairs {
                let first = &self.candidates[i];
                let second = &self.candidates[j];
                if first.volume >= target.volume || second.volume >= target.volume {
                    continue;
                }
                for &(k, l) in second_pairs {
                    let third = &self.candidates[k];
                    let fourth = &self.candidates[l];
                    if third.volume >= target.volume || fourth.volume >= target.volume {
                        continue;
                    }
                    let mut decomposition = vec![
                        first.class.clone(),
                        second.class.clone(),
                        third.class.clone(),
                        fourth.class.clone(),
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

        let mut simps: Vec<Vec<usize>> = simp_2d.into_iter().collect();
        simps.sort();
        for i in 0..simps.len() {
            for j in i..simps.len() {
                let s1: HashSet<usize> = simps[i].iter().copied().collect();
                let s2: HashSet<usize> = simps[j].iter().copied().collect();
                let mut comm: Vec<usize> = s1.intersection(&s2).copied().collect();
                comm.sort_unstable();
                if comm.len() != 2 {
                    continue;
                }
                let mut diff: Vec<usize> = s1.symmetric_difference(&s2).copied().collect();
                diff.sort_unstable();
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

    for witness in compute_origin_circuit_curve_witnesses(
        &pts_ext,
        &facets,
        &simp_2d_all,
        origin_idx,
        Some(&point_face_dims),
    ) {
        if let Some(gv) = resolved_conifold_origin_circuit_gv(&witness.class, origin_idx) {
            insert_toric_gv(&mut gv_by_class, witness.class, gv)?;
        }
    }

    let mut out: Vec<ToricCurveGvInvariant> = gv_by_class
        .into_iter()
        .map(|(class, gv)| ToricCurveGvInvariant { class, gv })
        .collect();
    out.sort_by(|a, b| a.class.cmp(&b.class));
    Ok(out)
}

/// Compute toric curve GV values with local formula provenance.
pub fn compute_toric_curve_gv_diagnostics(
    tri: &Triangulation,
    points: &[Point],
    polytope: &Polytope,
) -> Result<Vec<ToricCurveGvDiagnostic>> {
    if points.is_empty() {
        return Err(Error::InvalidInput("No points provided".into()));
    }
    if polytope.dim() != 4 {
        return Err(Error::InvalidInput(
            "toric curve GV diagnostics are only implemented for 4D polytopes".into(),
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

    let mut gv_by_class: HashMap<Vec<i64>, ToricCurveGvDiagnostic> = HashMap::new();
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

        let mut simps: Vec<Vec<usize>> = simp_2d.into_iter().collect();
        simps.sort();
        for i in 0..simps.len() {
            for j in i..simps.len() {
                let s1: HashSet<usize> = simps[i].iter().copied().collect();
                let s2: HashSet<usize> = simps[j].iter().copied().collect();
                let mut comm: Vec<usize> = s1.intersection(&s2).copied().collect();
                comm.sort_unstable();
                if comm.len() != 2 {
                    continue;
                }
                let mut diff: Vec<usize> = s1.symmetric_difference(&s2).copied().collect();
                diff.sort_unstable();
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
                let edge_coefficients = (
                    class.get(comm[0]).copied().unwrap_or(0),
                    class.get(comm[1]).copied().unwrap_or(0),
                );
                let source = ToricCurveGvSource::TwoFace {
                    edge: comm.clone(),
                    two_face_points: face.points.clone(),
                    two_face_genus: face.genus,
                    edge_coefficients,
                    edge_face_dimensions: (
                        point_face_dims.get(comm[0]).copied(),
                        point_face_dims.get(comm[1]).copied(),
                    ),
                    edge_one_face_genera: (
                        one_face_genera.get(comm[0]).copied().flatten(),
                        one_face_genera.get(comm[1]).copied().flatten(),
                    ),
                };
                insert_toric_gv_diagnostic(&mut gv_by_class, class, gv, source)?;
            }
        }
    }

    for witness in compute_origin_circuit_curve_witnesses(
        &pts_ext,
        &facets,
        &simp_2d_all,
        origin_idx,
        Some(&point_face_dims),
    ) {
        if let Some(gv) = resolved_conifold_origin_circuit_gv(&witness.class, origin_idx) {
            let source = ToricCurveGvSource::ResolvedConifoldOriginCircuit {
                origin_index: origin_idx,
                witness: witness.clone(),
            };
            insert_toric_gv_diagnostic(&mut gv_by_class, witness.class, gv, source)?;
        }
    }

    let mut out: Vec<ToricCurveGvDiagnostic> = gv_by_class.into_values().collect();
    out.sort_by(|a, b| a.class.cmp(&b.class));
    Ok(out)
}

/// Compute the origin-circuit Mori-cap curve classes separately from the GV
/// formulas that currently support only the resolved-conifold subset.
pub fn compute_origin_circuit_curve_diagnostics(
    tri: &Triangulation,
    points: &[Point],
    polytope: &Polytope,
) -> Result<Vec<OriginCircuitCurveDiagnostic>> {
    if points.is_empty() {
        return Err(Error::InvalidInput("No points provided".into()));
    }
    if polytope.dim() != 4 {
        return Err(Error::InvalidInput(
            "origin-circuit diagnostics are only implemented for 4D polytopes".into(),
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
    let (facets, twofaces) = compute_faces_4d(points, polytope)?;
    let point_face_dims = classify_primal_point_face_dimensions(points, polytope)?;
    let origin_idx = points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .ok_or_else(|| Error::InvalidInput("Origin not found in points".into()))?;

    let mut simp_2d_all: HashSet<Vec<usize>> = HashSet::new();
    for f in &twofaces {
        if f.len() < 4 {
            let mut f_sorted = f.clone();
            f_sorted.sort_unstable();
            simp_2d_all.insert(f_sorted);
            continue;
        }

        let face_pts: HashSet<usize> = f.iter().copied().collect();
        for simplex in tri.simplices() {
            let inter: Vec<usize> = simplex
                .iter()
                .filter(|idx| face_pts.contains(idx))
                .copied()
                .collect();
            if inter.len() == 3 {
                let mut inter_sorted = inter;
                inter_sorted.sort_unstable();
                simp_2d_all.insert(inter_sorted);
            }
        }
    }

    let mut witnesses_by_class: BTreeMap<Vec<i64>, Vec<OriginCircuitCurveWitness>> =
        BTreeMap::new();
    for witness in compute_origin_circuit_curve_witnesses(
        &pts_ext,
        &facets,
        &simp_2d_all,
        origin_idx,
        Some(&point_face_dims),
    ) {
        witnesses_by_class
            .entry(witness.class.clone())
            .or_default()
            .push(witness);
    }

    Ok(witnesses_by_class
        .into_iter()
        .map(|(class, witnesses)| {
            origin_circuit_diagnostic_from_class_and_witnesses(class, origin_idx, witnesses)
        })
        .collect())
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

fn insert_toric_gv_diagnostic(
    gv_by_class: &mut HashMap<Vec<i64>, ToricCurveGvDiagnostic>,
    class: Vec<i64>,
    gv: Integer,
    source: ToricCurveGvSource,
) -> Result<()> {
    match gv_by_class.entry(class.clone()) {
        std::collections::hash_map::Entry::Occupied(mut existing) => {
            if existing.get().gv != gv {
                return Err(Error::InvalidInput(format!(
                    "conflicting toric GV diagnostic values for duplicate curve class: {} vs {gv}",
                    existing.get().gv
                )));
            }
            existing.get_mut().sources.push(source);
        }
        std::collections::hash_map::Entry::Vacant(slot) => {
            slot.insert(ToricCurveGvDiagnostic {
                class,
                gv,
                sources: vec![source],
            });
        }
    }
    Ok(())
}

fn gv_lattice_search_request(
    gen_min_points: usize,
    max_deg: Option<u32>,
    lattice_augmentation: GvLatticeAugmentation,
) -> (Option<usize>, Option<i64>) {
    match lattice_augmentation {
        GvLatticeAugmentation::CytoolsDefault => (Some(gen_min_points), None),
        GvLatticeAugmentation::DegreeBoundedDiagnostic => match max_deg {
            Some(degree) => (None, Some(i64::from(degree))),
            None => (Some(gen_min_points), None),
        },
        GvLatticeAugmentation::None => (None, None),
    }
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

fn validate_exact_lp_class(class: &[i64]) -> Result<()> {
    const MAX_EXACT_F64_INTEGER: i128 = 1_i128 << 53;
    for &value in class {
        if i128::from(value).abs() > MAX_EXACT_F64_INTEGER {
            return Err(Error::InvalidInput(
                "curve class entry is too large for exact LP conversion".into(),
            ));
        }
    }
    Ok(())
}

fn verify_semigroup_decomposition(
    target: &ToricCurveCandidate,
    terms: &[CurveDecompositionTerm],
) -> Result<()> {
    let dim = target.class.len();
    let mut sum = vec![0_i128; dim];
    let mut total_terms = 0_u64;
    for term in terms {
        if term.class.len() != dim {
            return Err(Error::InvalidInput(
                "semigroup decomposition term dimension mismatch".into(),
            ));
        }
        if term.multiplicity == 0 {
            return Err(Error::InvalidInput(
                "semigroup decomposition contains zero multiplicity".into(),
            ));
        }
        total_terms = total_terms.checked_add(term.multiplicity).ok_or_else(|| {
            Error::InvalidInput("semigroup decomposition term count overflow".into())
        })?;
        let multiplicity = i128::from(term.multiplicity);
        for (coord_sum, &coefficient) in sum.iter_mut().zip(term.class.iter()) {
            let contribution = i128::from(coefficient)
                .checked_mul(multiplicity)
                .ok_or_else(|| {
                    Error::InvalidInput("semigroup decomposition coefficient overflow".into())
                })?;
            *coord_sum = coord_sum.checked_add(contribution).ok_or_else(|| {
                Error::InvalidInput("semigroup decomposition coordinate overflow".into())
            })?;
        }
    }
    if total_terms < 2 {
        return Err(Error::InvalidInput(
            "semigroup decomposition must use at least two terms".into(),
        ));
    }
    for (actual, &expected) in sum.iter().zip(target.class.iter()) {
        if *actual != i128::from(expected) {
            return Err(Error::InvalidInput(
                "semigroup decomposition failed exact integer verification".into(),
            ));
        }
    }
    Ok(())
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

fn compute_origin_circuit_curve_witnesses(
    pts_ext: &[Vec<i64>],
    facets: &[Vec<usize>],
    simp_2d_all: &HashSet<Vec<usize>>,
    origin_idx: usize,
    point_face_dims: Option<&[usize]>,
) -> Vec<OriginCircuitCurveWitness> {
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
                let class = normalized_row_from_sparse_relation(pts_ext.len(), full_v.clone());
                out.push(OriginCircuitCurveWitness {
                    class,
                    first_facet_exclusive_point: *p1,
                    second_facet_exclusive_point: *p2,
                    shared_two_simplex: s2d.clone(),
                    first_facet: f1.clone(),
                    second_facet: f2.clone(),
                    relation_points: origin_circuit_relation_points(
                        pts_ext,
                        point_face_dims,
                        &full_v,
                    ),
                    sparse_relation: full_v,
                });
            }
        }
    }
    out
}

fn origin_circuit_relation_points(
    pts_ext: &[Vec<i64>],
    point_face_dims: Option<&[usize]>,
    sparse_relation: &[(usize, i64)],
) -> Vec<OriginCircuitRelationPoint> {
    sparse_relation
        .iter()
        .filter_map(|&(point_index, coefficient)| {
            if coefficient == 0 {
                return None;
            }
            let point = pts_ext.get(point_index)?;
            let coordinates = point
                .get(..point.len().saturating_sub(1))
                .unwrap_or(point)
                .to_vec();
            Some(OriginCircuitRelationPoint {
                point_index,
                coefficient,
                coordinates,
                face_dimension: point_face_dims.and_then(|dims| dims.get(point_index).copied()),
            })
        })
        .collect()
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

    if (non_origin_neg_ones == 1 && pos_ones == 2) || (non_origin_neg_ones == 2 && pos_ones == 3) {
        // An isolated resolved conifold curve has normal bundle
        // O(-1) + O(-1), so section 6 gives GV^0 = 1.
        Some(Integer::from(1))
    } else {
        None
    }
}

fn origin_circuit_diagnostic_from_class_and_witnesses(
    class: Vec<i64>,
    origin_idx: usize,
    mut witnesses: Vec<OriginCircuitCurveWitness>,
) -> OriginCircuitCurveDiagnostic {
    let origin_coefficient = class.get(origin_idx).copied().unwrap_or(0);
    let mut negative_coefficient_counts = BTreeMap::new();
    let mut positive_coefficient_counts = BTreeMap::new();
    for (idx, &coefficient) in class.iter().enumerate() {
        if idx == origin_idx || coefficient == 0 {
            continue;
        }
        if coefficient < 0 {
            *negative_coefficient_counts.entry(coefficient).or_insert(0) += 1;
        } else {
            *positive_coefficient_counts.entry(coefficient).or_insert(0) += 1;
        }
    }
    let is_resolved_conifold_pattern =
        resolved_conifold_origin_circuit_gv(&class, origin_idx).is_some();
    witnesses.sort_by(|left, right| {
        left.shared_two_simplex
            .cmp(&right.shared_two_simplex)
            .then_with(|| {
                left.first_facet_exclusive_point
                    .cmp(&right.first_facet_exclusive_point)
            })
            .then_with(|| {
                left.second_facet_exclusive_point
                    .cmp(&right.second_facet_exclusive_point)
            })
            .then_with(|| left.sparse_relation.cmp(&right.sparse_relation))
    });

    OriginCircuitCurveDiagnostic {
        class,
        origin_coefficient,
        negative_coefficient_counts,
        positive_coefficient_counts,
        is_resolved_conifold_pattern,
        witnesses,
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GvLatticeAugmentation {
    CytoolsDefault,
    DegreeBoundedDiagnostic,
    None,
}

impl GvLatticeAugmentation {
    const fn as_str(self) -> &'static str {
        match self {
            Self::CytoolsDefault => "cytools-default",
            Self::DegreeBoundedDiagnostic => "degree-bounded-diagnostic",
            Self::None => "none",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum GvCachePolicy {
    Enabled,
    Disabled,
}

impl GvCachePolicy {
    fn from_env() -> Self {
        match env::var("CYRUS_GV_CACHE") {
            Ok(value) if value == "0" => Self::Disabled,
            _ => Self::Enabled,
        }
    }
}

fn gv_lattice_augmentation_grading(
    rays: &[Vec<i64>],
    semigroup_grading_vector: &[i64],
    lattice_augmentation: GvLatticeAugmentation,
) -> Result<Vec<i64>> {
    match lattice_augmentation {
        GvLatticeAugmentation::CytoolsDefault => compute_grading_vector(rays).ok_or_else(|| {
            Error::InvalidInput("failed to compute CYTools default lattice grading vector".into())
        }),
        GvLatticeAugmentation::DegreeBoundedDiagnostic | GvLatticeAugmentation::None => {
            Ok(semigroup_grading_vector.to_vec())
        }
    }
}

/// Compute GV invariants using cygv.
///
/// This mirrors the CYTools default wrapper: the supplied Mori-cap rays are
/// augmented with `mori.find_lattice_points(min_points=100*h11)` before calling
/// `cygv`, even when the final `cygv` semigroup truncation uses `max_deg`.
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
    compute_gv_invariants_inner(
        rays,
        grading_vector,
        q_matrix,
        intnums,
        min_points,
        max_deg,
        GvLatticeAugmentation::CytoolsDefault,
        GvCachePolicy::from_env(),
    )
}

/// Compute GV invariants with a degree-bounded lattice augmentation.
///
/// This is a diagnostic shortcut, not the CYTools wrapper contract. When
/// `max_deg` is supplied it enumerates Mori-cone lattice points only up to that
/// degree before calling `cygv`, which can make bounded investigations
/// tractable but may omit generators that CYTools would provide before cygv's
/// own semigroup truncation.
///
/// # Errors
/// Returns an error if the input shapes or numeric ranges are invalid.
pub fn compute_gv_invariants_with_degree_bounded_lattice(
    rays: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<Vec<(Vec<i32>, Integer)>> {
    compute_gv_invariants_inner(
        rays,
        grading_vector,
        q_matrix,
        intnums,
        min_points,
        max_deg,
        GvLatticeAugmentation::DegreeBoundedDiagnostic,
        GvCachePolicy::from_env(),
    )
}

/// Compute GV invariants using exactly the caller-provided semigroup generators.
///
/// This mirrors the CYTools `_compute_gvs_gws(..., mcap_generators=...)` path:
/// the supplied generators are passed to `cygv` without first augmenting them
/// with Mori-cone lattice points. Use this for explicit diagnostics or for
/// callers that have already constructed the desired semigroup elements.
///
/// # Errors
/// Returns an error if the input shapes or numeric ranges are invalid.
pub fn compute_gv_invariants_with_provided_generators(
    generators: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<Vec<(Vec<i32>, Integer)>> {
    compute_gv_invariants_inner(
        generators,
        grading_vector,
        q_matrix,
        intnums,
        min_points,
        max_deg,
        GvLatticeAugmentation::None,
        GvCachePolicy::from_env(),
    )
}

/// Compute GV invariants using exactly the caller-provided semigroup generators
/// and an explicit `cygv` nef partition.
///
/// This is the complete-intersection analogue of
/// [`compute_gv_invariants_with_provided_generators`]. Cyrus still delegates the
/// HKTY computation to upstream `cygv`; this function only prepares the finite
/// semigroup, GLSM matrix, intersection tensor, and nef partition boundary.
///
/// # Errors
/// Returns an error if the input shapes or numeric ranges are invalid, if `cygv`
/// cannot construct the semigroup, or if the nef partition is inconsistent with
/// the supplied charge matrix.
pub fn compute_gv_invariants_with_provided_generators_and_nef_partition(
    generators: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    nef_partition: &[Vec<usize>],
    intnums: &Intersection,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<Vec<(Vec<i32>, Integer)>> {
    let (semigroup, q, intnums_map) = provided_cygv_semigroup_inputs(
        generators,
        grading_vector,
        q_matrix,
        intnums,
        min_points,
        max_deg,
    )?;
    let nefpart = cygv_nef_partition(nef_partition, q.nrows())?;
    compute_cygv_rat_threefold_from_semigroup_with_nefpart(
        semigroup,
        &q,
        &nefpart,
        intnums_map,
        "provided-generator complete-intersection GV",
    )
}

/// Compute GV invariants and compact `q_N` polynomial history from
/// caller-provided semigroup generators.
///
/// This mirrors [`compute_gv_invariants_with_provided_generators`] at the
/// `cygv` boundary, but returns the `q_N` polynomials materialized by cygv's
/// series inversion. It intentionally bypasses the normal GV cache because the
/// trace is diagnostic state, not just the final invariant table.
///
/// # Errors
/// Returns an error if the input shapes or numeric ranges are invalid, if cygv
/// cannot construct the semigroup, or if HKTY finds non-integral output.
pub fn compute_gv_invariants_with_provided_generators_qn_trace(
    generators: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<GvInvariantsWithQnTrace> {
    let (semigroup, q, intnums_map) = provided_cygv_semigroup_inputs(
        generators,
        grading_vector,
        q_matrix,
        intnums,
        min_points,
        max_deg,
    )?;
    compute_cygv_rat_threefold_from_semigroup_with_qn_trace(
        semigroup,
        &q,
        intnums_map,
        "provided-generator GV qN trace",
    )
}

/// Compute raw GW/GV coefficient candidates for caller-provided semigroup
/// generators without enforcing GV integrality.
///
/// This is a diagnostic escape hatch for upstream `cygv` failures: when
/// `FIND_GV=true` rejects a non-integral candidate, the vendored `cygv`
/// diagnostic identifies the first failing candidate. Running the same input
/// with `FIND_GV=false` exposes the surrounding coefficient candidates read
/// from the instanton polynomial.
///
/// # Errors
/// Returns an error if input construction, fundamental-period computation, or
/// instanton-data computation fails.
#[doc(hidden)]
pub fn compute_gw_coefficient_trace_with_provided_generators(
    generators: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<Vec<CygvGvCoefficientTrace>> {
    let (semigroup, q, intnums_map) = provided_cygv_semigroup_inputs(
        generators,
        grading_vector,
        q_matrix,
        intnums,
        min_points,
        max_deg,
    )?;
    compute_cygv_rat_threefold_gw_coefficient_trace_from_semigroup(
        semigroup,
        &q,
        intnums_map,
        "provided-generator GW coefficient trace",
    )
}

/// Source-audit helper for the private `cygv::Semigroup::with_max_degree`
/// seed-reduction step.
///
/// In `cygv` 0.1.2, the supplied degree-trimmed seed elements are first
/// converted to a set, the zero vector is removed, and any seed that can be
/// written as a sum of two seeds is removed before the additive closure is
/// generated. This function exposes that exact pair-sum pruning as a
/// diagnostic primitive without running the expensive closure step. It is not
/// a compact GV entry point; production GV computations should call
/// [`compute_gv_invariants`] or another wrapper that runs upstream `cygv`.
///
/// # Errors
/// Returns an error if no elements are supplied, if the dimension is zero, if
/// row dimensions are inconsistent, or if an intermediate sum overflows `i64`.
#[doc(hidden)]
pub fn cygv_pair_reduced_seed_generators(elements: &[Vec<i64>]) -> Result<Vec<Vec<i64>>> {
    let Some(first) = elements.first() else {
        return Err(Error::InvalidInput(
            "cygv seed reduction requires at least one element".into(),
        ));
    };
    let dimension = first.len();
    if dimension == 0 {
        return Err(Error::InvalidInput(
            "cygv seed reduction dimension is zero".into(),
        ));
    }
    if elements.iter().any(|row| row.len() != dimension) {
        return Err(Error::InvalidInput(
            "cygv seed reduction rows have inconsistent dimensions".into(),
        ));
    }

    let zero = vec![0i64; dimension];
    let mut generator_set = elements.iter().cloned().collect::<HashSet<_>>();
    generator_set.remove(&zero);
    let original_generators = generator_set.iter().cloned().collect::<Vec<_>>();
    let mut to_remove = HashSet::new();
    let mut sum = vec![0i64; dimension];

    for lhs_idx in 0..original_generators.len() {
        for rhs_idx in lhs_idx..original_generators.len() {
            for ((slot, &lhs), &rhs) in sum
                .iter_mut()
                .zip(original_generators[lhs_idx].iter())
                .zip(original_generators[rhs_idx].iter())
            {
                *slot = lhs.checked_add(rhs).ok_or_else(|| {
                    Error::InvalidInput("cygv seed reduction sum overflowed i64".into())
                })?;
            }
            if generator_set.contains(&sum) {
                to_remove.insert(sum.clone());
            }
        }
    }

    for generator in &to_remove {
        generator_set.remove(generator);
    }
    let mut generators = generator_set.into_iter().collect::<Vec<_>>();
    generators.sort();
    Ok(generators)
}

/// Compute GV invariants using an explicitly truncated semigroup.
///
/// Unlike [`compute_gv_invariants_with_provided_generators`], this does not ask
/// `cygv` to close the provided rows under addition. The rows are passed to
/// `Semigroup::from_data` as the exact truncation domain for HKTY, with the
/// identity element inserted if absent. This is the right entry point for
/// diagnostics that construct a causal-diamond or face-local semigroup
/// explicitly.
///
/// The HKTY stages still run inside the upstream `cygv` crate. Cyrus only
/// prepares the explicit semigroup data here.
///
/// # Errors
/// Returns an error if the semigroup, grading, GLSM charge matrix, or
/// intersection numbers are inconsistent, or if HKTY finds non-integral GV
/// output for the supplied truncation.
#[allow(clippy::too_many_lines)]
pub fn compute_gv_invariants_with_explicit_semigroup(
    elements: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
) -> Result<Vec<(Vec<i32>, Integer)>> {
    let (semigroup, q, intnums_map) =
        explicit_cygv_semigroup_inputs(elements, grading_vector, q_matrix, intnums)?;
    compute_cygv_rat_threefold_from_semigroup(semigroup, &q, intnums_map, "explicit GV semigroup")
}

/// Compute GV invariants using an explicitly truncated semigroup and an
/// explicit `cygv` nef partition.
///
/// This is intended for source-derived local complete-intersection diagnostics
/// where the finite HKTY domain is already known. It does not close the supplied
/// elements under addition; it passes the exact semigroup to `cygv`.
///
/// # Errors
/// Returns an error if the semigroup, grading, GLSM charge matrix, nef
/// partition, or intersection numbers are inconsistent.
#[allow(clippy::too_many_lines)]
pub fn compute_gv_invariants_with_explicit_semigroup_and_nef_partition(
    elements: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    nef_partition: &[Vec<usize>],
    intnums: &Intersection,
) -> Result<Vec<(Vec<i32>, Integer)>> {
    let (semigroup, q, intnums_map) =
        explicit_cygv_semigroup_inputs(elements, grading_vector, q_matrix, intnums)?;
    let nefpart = cygv_nef_partition(nef_partition, q.nrows())?;
    compute_cygv_rat_threefold_from_semigroup_with_nefpart(
        semigroup,
        &q,
        &nefpart,
        intnums_map,
        "explicit complete-intersection GV semigroup",
    )
}

/// Compute GV invariants and compact `q_N` polynomial history using an
/// explicitly truncated semigroup.
///
/// This is a diagnostic entry point for understanding cygv's degree-ordered
/// subtraction history. It still runs the upstream cygv HKTY implementation;
/// Cyrus only exports the qN polynomials that cygv materializes internally.
///
/// # Errors
/// Returns an error under the same conditions as
/// [`compute_gv_invariants_with_explicit_semigroup`].
#[allow(clippy::too_many_lines)]
pub fn compute_gv_invariants_with_explicit_semigroup_qn_trace(
    elements: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
) -> Result<GvInvariantsWithQnTrace> {
    let (semigroup, q, intnums_map) =
        explicit_cygv_semigroup_inputs(elements, grading_vector, q_matrix, intnums)?;
    compute_cygv_rat_threefold_from_semigroup_with_qn_trace(
        semigroup,
        &q,
        intnums_map,
        "explicit GV semigroup qN trace",
    )
}

/// Compute raw GW/GV coefficient candidates using an explicitly truncated
/// semigroup without enforcing GV integrality.
///
/// This is the explicit-semigroup analogue of
/// [`compute_gw_coefficient_trace_with_provided_generators`]. It is intended
/// for diagnostics of small source/chamber candidate domains where the normal
/// integral GV path fails with a non-integer candidate.
///
/// # Errors
/// Returns an error if the explicit semigroup, grading, GLSM charge matrix, or
/// intersection numbers are inconsistent, or if cygv fails before coefficient
/// readout.
#[doc(hidden)]
#[allow(clippy::too_many_lines)]
pub fn compute_gw_coefficient_trace_with_explicit_semigroup(
    elements: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
) -> Result<Vec<CygvGvCoefficientTrace>> {
    let (semigroup, q, intnums_map) =
        explicit_cygv_semigroup_inputs(elements, grading_vector, q_matrix, intnums)?;
    compute_cygv_rat_threefold_gw_coefficient_trace_from_semigroup(
        semigroup,
        &q,
        intnums_map,
        "explicit GW semigroup coefficient trace",
    )
}

fn provided_cygv_semigroup_inputs(
    generators: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    min_points: Option<u32>,
    max_deg: Option<u32>,
) -> Result<(
    cygv::Semigroup,
    DMatrix<i32>,
    HashMap<(usize, usize, usize), i32>,
)> {
    if min_points.is_some() == max_deg.is_some() {
        return Err(Error::InvalidInput(
            "Exactly one of min_points or max_deg must be specified".into(),
        ));
    }
    if generators.is_empty() {
        return Err(Error::InvalidInput(
            "GV computation requires at least one generator".into(),
        ));
    }
    let dim = generators[0].len();
    if dim == 0 {
        return Err(Error::InvalidInput(
            "GV generator dimension must be positive".into(),
        ));
    }
    if grading_vector.len() != dim {
        return Err(Error::InvalidInput(
            "grading vector length must match GV generator dimension".into(),
        ));
    }

    let grading_vec_i32: Vec<i32> = grading_vector
        .iter()
        .map(|&v| i32::try_from(v))
        .collect::<std::result::Result<_, _>>()
        .map_err(|_| Error::InvalidInput("grading vector does not fit in i32".into()))?;
    let grading = RowDVector::from_row_slice(&grading_vec_i32);
    let q = cygv_q_matrix(q_matrix, dim)?;
    let intnums_map = cygv_intnums_map(intnums)?;

    let mut unique = HashSet::new();
    let mut filtered_vecs = Vec::new();
    for row in generators {
        if row.len() != dim {
            return Err(Error::InvalidInput(
                "GV generator rows have inconsistent dimensions".into(),
            ));
        }
        if unique.insert(row.clone()) {
            filtered_vecs.push(row.clone());
        }
    }
    if let Some(d) = max_deg {
        let d_i128 = i128::from(d);
        filtered_vecs.retain(|row| {
            row.iter()
                .zip(grading_vector.iter())
                .map(|(&x, &g)| i128::from(x) * i128::from(g))
                .sum::<i128>()
                <= d_i128
        });
        if filtered_vecs.is_empty() {
            return Err(Error::InvalidInput(
                "No generators remain after max_deg filtering".into(),
            ));
        }
    }

    let n_gen = filtered_vecs.len();
    let mut gen_data = Vec::with_capacity(dim * n_gen);
    for col in 0..n_gen {
        for row in 0..dim {
            let val = i32::try_from(filtered_vecs[col][row]).map_err(|_| {
                Error::InvalidInput("Mori cone generator does not fit in i32".into())
            })?;
            gen_data.push(val);
        }
    }
    let generators = DMatrix::from_column_slice(dim, n_gen, &gen_data);
    let semigroup = if let Some(d) = max_deg {
        cygv::Semigroup::with_max_degree(generators, grading, d)
    } else if let Some(n) = min_points {
        cygv::Semigroup::with_min_elements(generators, grading, n as usize)
    } else {
        cygv::Semigroup::from_data(generators, grading)
    }
    .map_err(|e| Error::InvalidInput(format!("cygv semigroup construction failed: {e}")))?;

    Ok((semigroup, q, intnums_map))
}

fn explicit_cygv_semigroup_inputs(
    elements: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
) -> Result<(
    cygv::Semigroup,
    DMatrix<i32>,
    HashMap<(usize, usize, usize), i32>,
)> {
    if elements.is_empty() {
        return Err(Error::InvalidInput(
            "explicit GV semigroup elements are empty".into(),
        ));
    }
    let dim = elements[0].len();
    if dim == 0 {
        return Err(Error::InvalidInput(
            "explicit GV semigroup dimension is zero".into(),
        ));
    }
    if grading_vector.len() != dim {
        return Err(Error::InvalidInput(
            "grading vector length must match explicit semigroup dimension".into(),
        ));
    }

    let grading_vec_i32: Vec<i32> = grading_vector
        .iter()
        .map(|&v| i32::try_from(v))
        .collect::<std::result::Result<_, _>>()
        .map_err(|_| Error::InvalidInput("grading vector does not fit in i32".into()))?;
    let grading = RowDVector::from_row_slice(&grading_vec_i32);
    let q = cygv_q_matrix(q_matrix, dim)?;
    let intnums_map = cygv_intnums_map(intnums)?;

    let mut unique = HashSet::new();
    let zero = vec![0i64; dim];
    unique.insert(zero.clone());
    let mut semigroup_rows = vec![zero];
    for row in elements {
        if row.len() != dim {
            return Err(Error::InvalidInput(
                "explicit GV semigroup rows have inconsistent dimensions".into(),
            ));
        }
        if unique.insert(row.clone()) {
            semigroup_rows.push(row.clone());
        }
    }

    let n_elements = semigroup_rows.len();
    let mut element_data = Vec::with_capacity(dim * n_elements);
    for col in 0..n_elements {
        for row in 0..dim {
            let val = i32::try_from(semigroup_rows[col][row]).map_err(|_| {
                Error::InvalidInput("explicit GV semigroup element does not fit in i32".into())
            })?;
            element_data.push(val);
        }
    }
    let elements = DMatrix::from_column_slice(dim, n_elements, &element_data);
    let semigroup = cygv::Semigroup::from_data(elements, grading)
        .map_err(|e| Error::InvalidInput(format!("explicit GV semigroup is inconsistent: {e}")))?;
    Ok((semigroup, q, intnums_map))
}

fn compute_gv_invariants_inner(
    rays: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    min_points: Option<u32>,
    max_deg: Option<u32>,
    lattice_augmentation: GvLatticeAugmentation,
    cache_policy: GvCachePolicy,
) -> Result<Vec<(Vec<i32>, Integer)>> {
    let t0 = std::time::Instant::now();
    eprintln!(
        "[DEBUG] gv start: rays={}, dim={}, max_deg={:?}, min_points={:?}, lattice_augmentation={}",
        rays.len(),
        rays.first().map_or(0, Vec::len),
        max_deg,
        min_points,
        lattice_augmentation.as_str()
    );
    if min_points.is_some() == max_deg.is_some() {
        return Err(Error::InvalidInput(
            "Exactly one of min_points or max_deg must be specified".into(),
        ));
    }
    if rays.is_empty() {
        return Err(Error::InvalidInput(
            "GV computation requires at least one generator".into(),
        ));
    }

    let dim = rays[0].len();
    if dim == 0 {
        return Err(Error::InvalidInput(
            "GV generator dimension must be positive".into(),
        ));
    }
    if rays.iter().any(|row| row.len() != dim) {
        return Err(Error::InvalidInput(
            "GV generator rows have inconsistent dimensions".into(),
        ));
    }
    if grading_vector.len() != dim {
        return Err(Error::InvalidInput(
            "grading vector length must match GV generator dimension".into(),
        ));
    }

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

    // CYTools default: lattice_pts = mori.find_lattice_points(min_points=100*h11)
    let factor = env::var("CYRUS_LATTICE_MIN_POINTS_FACTOR")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(100);
    let gen_min_points = factor * dim;
    let (lattice_min_points, lattice_max_deg) =
        gv_lattice_search_request(gen_min_points, max_deg, lattice_augmentation);
    eprintln!(
        "[DEBUG] gv lattice request: min_points={:?} max_deg={:?}",
        lattice_min_points, lattice_max_deg
    );

    let lattice_pts = if lattice_augmentation != GvLatticeAugmentation::None {
        let lattice_cache = LatticeCacheControls::from_env(1000, 0);
        let lattice_grading =
            gv_lattice_augmentation_grading(rays, grading_vector, lattice_augmentation)?;
        if lattice_grading.len() != dim {
            return Err(Error::InvalidInput(
                "lattice grading vector length must match GV generator dimension".into(),
            ));
        }
        if lattice_augmentation == GvLatticeAugmentation::CytoolsDefault
            && lattice_grading != grading_vector
        {
            eprintln!(
                "[DEBUG] gv lattice grading differs from cygv semigroup grading: lattice={:?} semigroup={:?}",
                lattice_grading, grading_vector
            );
        }
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        LATTICE_CACHE_VERSION.hash(&mut hasher);
        let mut key_rays: Vec<&[i64]> = rays.iter().map(Vec::as_slice).collect();
        key_rays.sort();
        for row in &key_rays {
            row.hash(&mut hasher);
        }
        lattice_grading.hash(&mut hasher);
        lattice_min_points.hash(&mut hasher);
        lattice_max_deg.hash(&mut hasher);
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
                lattice_min_points,
                lattice_max_deg,
                &lattice_grading,
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
    } else {
        eprintln!("[DEBUG] gv lattice points: skipped (using caller-provided generators only)");
        Vec::new()
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

    let q = cygv_q_matrix(q_matrix, dim)?;
    eprintln!(
        "[DEBUG] gv q shape: {}x{}, elapsed={:.2?}",
        q.nrows(),
        q.ncols(),
        t0.elapsed()
    );

    let intnums_map = cygv_intnums_map(intnums)?;

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
        GV_CACHE_VERSION.hash(&mut hasher);
        lattice_augmentation.as_str().hash(&mut hasher);
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
    let cache_enabled = cache_policy == GvCachePolicy::Enabled;
    eprintln!(
        "[DEBUG] gv cache path: {} ({})",
        gv_cache_path.display(),
        if cache_enabled { "enabled" } else { "disabled" }
    );
    if cache_enabled && gv_cache_path.exists() {
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

    let semigroup = if let Some(d) = max_deg {
        cygv::Semigroup::with_max_degree(generators, grading, d)
    } else if let Some(n) = min_points {
        cygv::Semigroup::with_min_elements(generators, grading, n as usize)
    } else {
        cygv::Semigroup::from_data(generators, grading)
    }
    .map_err(|e| Error::InvalidInput(format!("cygv semigroup construction failed: {e}")))?;

    let out =
        compute_cygv_rat_threefold_from_semigroup(semigroup, &q, intnums_map, "cygv GV wrapper")?;

    if !cache_enabled {
        eprintln!("[DEBUG] gv invariants: cache write skipped");
    } else if let Err(e) = fs::create_dir_all(&cache_dir) {
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

fn compute_cygv_rat_threefold_from_semigroup(
    semigroup: cygv::Semigroup,
    q: &DMatrix<i32>,
    intnums_map: HashMap<(usize, usize, usize), i32>,
    context: &str,
) -> Result<Vec<(Vec<i32>, Integer)>> {
    let nefpart: Vec<DVector<i32>> = Vec::new();
    compute_cygv_rat_threefold_from_semigroup_with_nefpart(
        semigroup,
        q,
        &nefpart,
        intnums_map,
        context,
    )
}

fn compute_cygv_rat_threefold_from_semigroup_with_nefpart(
    semigroup: cygv::Semigroup,
    q: &DMatrix<i32>,
    nefpart: &[DVector<i32>],
    intnums_map: HashMap<(usize, usize, usize), i32>,
    context: &str,
) -> Result<Vec<(Vec<i32>, Integer)>> {
    if cfg!(panic = "abort") {
        return Err(Error::InvalidInput(format!(
            "{context}: cygv HKTY execution requires a panic=unwind build because upstream cygv can still panic internally"
        )));
    }

    let previous_panic_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        compute_cygv_rat_threefold_from_semigroup_unchecked(
            semigroup,
            q,
            nefpart,
            intnums_map,
            context,
        )
    }));
    std::panic::set_hook(previous_panic_hook);

    match result {
        Ok(result) => result,
        Err(payload) => Err(Error::InvalidInput(format!(
            "{context}: cygv HKTY execution panicked: {}",
            panic_payload_message(payload.as_ref())
        ))),
    }
}

fn compute_cygv_rat_threefold_from_semigroup_with_qn_trace(
    semigroup: cygv::Semigroup,
    q: &DMatrix<i32>,
    intnums_map: HashMap<(usize, usize, usize), i32>,
    context: &str,
) -> Result<GvInvariantsWithQnTrace> {
    if cfg!(panic = "abort") {
        return Err(Error::InvalidInput(format!(
            "{context}: cygv HKTY execution requires a panic=unwind build because upstream cygv can still panic internally"
        )));
    }

    let previous_panic_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        compute_cygv_rat_threefold_from_semigroup_with_qn_trace_unchecked(
            semigroup,
            q,
            intnums_map,
            context,
        )
    }));
    std::panic::set_hook(previous_panic_hook);

    match result {
        Ok(result) => result,
        Err(payload) => Err(Error::InvalidInput(format!(
            "{context}: cygv HKTY execution panicked: {}",
            panic_payload_message(payload.as_ref())
        ))),
    }
}

fn compute_cygv_rat_threefold_gw_coefficient_trace_from_semigroup(
    semigroup: cygv::Semigroup,
    q: &DMatrix<i32>,
    intnums_map: HashMap<(usize, usize, usize), i32>,
    context: &str,
) -> Result<Vec<CygvGvCoefficientTrace>> {
    if cfg!(panic = "abort") {
        return Err(Error::InvalidInput(format!(
            "{context}: cygv HKTY execution requires a panic=unwind build because upstream cygv can still panic internally"
        )));
    }

    let previous_panic_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        compute_cygv_rat_threefold_gw_coefficient_trace_from_semigroup_unchecked(
            semigroup,
            q,
            intnums_map,
            context,
        )
    }));
    std::panic::set_hook(previous_panic_hook);

    match result {
        Ok(result) => result,
        Err(payload) => Err(Error::InvalidInput(format!(
            "{context}: cygv HKTY execution panicked: {}",
            panic_payload_message(payload.as_ref())
        ))),
    }
}

fn compute_cygv_rat_threefold_from_semigroup_unchecked(
    semigroup: cygv::Semigroup,
    q: &DMatrix<i32>,
    nefpart: &[DVector<i32>],
    intnums_map: HashMap<(usize, usize, usize), i32>,
    context: &str,
) -> Result<Vec<(Vec<i32>, Integer)>> {
    compute_cygv_rat_threefold_raw_from_semigroup_unchecked(
        semigroup,
        q,
        nefpart,
        intnums_map,
        context,
        false,
    )
    .map(|output| output.invariants)
}

fn compute_cygv_rat_threefold_from_semigroup_with_qn_trace_unchecked(
    semigroup: cygv::Semigroup,
    q: &DMatrix<i32>,
    intnums_map: HashMap<(usize, usize, usize), i32>,
    context: &str,
) -> Result<GvInvariantsWithQnTrace> {
    let nefpart: Vec<DVector<i32>> = Vec::new();
    compute_cygv_rat_threefold_raw_from_semigroup_unchecked(
        semigroup,
        q,
        &nefpart,
        intnums_map,
        context,
        true,
    )
}

fn compute_cygv_rat_threefold_gw_coefficient_trace_from_semigroup_unchecked(
    semigroup: cygv::Semigroup,
    q: &DMatrix<i32>,
    intnums_map: HashMap<(usize, usize, usize), i32>,
    context: &str,
) -> Result<Vec<CygvGvCoefficientTrace>> {
    let zero_cutoff = RugRational::new();
    let poly_props = cygv::PolynomialProperties::new(&semigroup, &zero_cutoff);
    let (intnum_dict, intnum_idxpairs, n_indices) = cygv::misc::process_int_nums(intnums_map, true)
        .map_err(|e| {
            Error::InvalidInput(format!(
                "{context}: cygv intersection preprocessing failed: {e}"
            ))
        })?;

    let n_threads = cygv_thread_count_from_env();
    let pool_size = cygv_pool_size_from_env();
    let main_pool = cygv::NumberPool::new(poly_props.zero_cutoff.clone(), pool_size);
    let thread_pools: Vec<_> = (0..n_threads)
        .map(|_| cygv::NumberPool::new(poly_props.zero_cutoff.clone(), pool_size))
        .collect();
    let mut all_pools = (main_pool, thread_pools);
    let nefpart: Vec<DVector<i32>> = Vec::new();

    let fp = cygv::fundamental_period::compute_omega(
        &poly_props,
        &semigroup,
        q,
        &nefpart,
        &intnum_idxpairs,
        &mut all_pools,
    )
    .map_err(|e| Error::InvalidInput(format!("{context}: cygv fundamental period failed: {e}")))?;

    let inst_data = cygv::instanton::compute_instanton_data(
        fp,
        &poly_props,
        &intnum_idxpairs,
        n_indices,
        &intnum_dict,
        true,
        &mut all_pools,
    )
    .map_err(|e| Error::InvalidInput(format!("{context}: cygv instanton data failed: {e}")))?;

    let (_, _, raw_gv_coefficient_trace) = cygv::series_inversion::invert_series_with_qn_trace::<
        RugRational,
        false,
        true,
    >(inst_data, &poly_props, &mut all_pools)
    .map_err(|e| Error::InvalidInput(format!("{context}: cygv series inversion failed: {e}")))?;

    Ok(convert_cygv_gv_coefficient_trace(raw_gv_coefficient_trace))
}

fn convert_cygv_gv_coefficient_trace(
    raw_gv_coefficient_trace: Vec<cygv::series_inversion::GvCoefficientTrace<RugRational>>,
) -> Vec<CygvGvCoefficientTrace> {
    raw_gv_coefficient_trace
        .into_iter()
        .map(|trace| CygvGvCoefficientTrace {
            element_index: trace.element_index,
            degree: trace.degree,
            element: trace.element,
            insertion_index: trace.insertion_index,
            pivot_component: trace.pivot_component,
            instanton_coefficient: trace
                .instanton_coefficient
                .map(|coefficient| coefficient.to_string()),
            gv_candidate: trace
                .gv_candidate
                .map(|coefficient| coefficient.to_string()),
            rounded_gv_candidate: trace
                .rounded_gv_candidate
                .map(|coefficient| coefficient.to_string()),
            status: trace.status.to_string(),
        })
        .collect()
}

fn compute_cygv_rat_threefold_raw_from_semigroup_unchecked(
    semigroup: cygv::Semigroup,
    q: &DMatrix<i32>,
    nefpart: &[DVector<i32>],
    intnums_map: HashMap<(usize, usize, usize), i32>,
    context: &str,
    collect_qn_trace: bool,
) -> Result<GvInvariantsWithQnTrace> {
    let zero_cutoff = RugRational::new();
    let poly_props = cygv::PolynomialProperties::new(&semigroup, &zero_cutoff);
    let (intnum_dict, intnum_idxpairs, n_indices) = cygv::misc::process_int_nums(intnums_map, true)
        .map_err(|e| {
            Error::InvalidInput(format!(
                "{context}: cygv intersection preprocessing failed: {e}"
            ))
        })?;

    let n_threads = cygv_thread_count_from_env();
    let pool_size = cygv_pool_size_from_env();
    let main_pool = cygv::NumberPool::new(poly_props.zero_cutoff.clone(), pool_size);
    let thread_pools: Vec<_> = (0..n_threads)
        .map(|_| cygv::NumberPool::new(poly_props.zero_cutoff.clone(), pool_size))
        .collect();
    let mut all_pools = (main_pool, thread_pools);
    let fp = cygv::fundamental_period::compute_omega(
        &poly_props,
        &semigroup,
        q,
        nefpart,
        &intnum_idxpairs,
        &mut all_pools,
    )
    .map_err(|e| Error::InvalidInput(format!("{context}: cygv fundamental period failed: {e}")))?;

    let inst_data = cygv::instanton::compute_instanton_data(
        fp,
        &poly_props,
        &intnum_idxpairs,
        n_indices,
        &intnum_dict,
        true,
        &mut all_pools,
    )
    .map_err(|e| Error::InvalidInput(format!("{context}: cygv instanton data failed: {e}")))?;

    let (gv, raw_qn_trace, raw_gv_coefficient_trace) = if collect_qn_trace {
        cygv::series_inversion::invert_series_with_qn_trace::<RugRational, true, true>(
            inst_data,
            &poly_props,
            &mut all_pools,
        )
    } else {
        cygv::series_inversion::invert_series::<RugRational, true, true>(
            inst_data,
            &poly_props,
            &mut all_pools,
        )
        .map(|gv| (gv, Vec::new(), Vec::new()))
    }
    .map_err(|e| Error::InvalidInput(format!("{context}: cygv series inversion failed: {e}")))?;

    let mut gv_sorted: Vec<_> = gv.into_iter().collect();
    gv_sorted
        .sort_unstable_by_key(|((element_idx, insertion_idx), _)| (*element_idx, *insertion_idx));

    let mut out = Vec::with_capacity(gv_sorted.len());
    for ((element_idx, _), gv_value) in gv_sorted {
        let (numer, denom) = gv_value.into_numer_denom();
        if denom != rug::Integer::from(1) {
            return Err(Error::InvalidInput(format!(
                "{context}: cygv produced non-integral invariant at element index {element_idx}: denominator {denom}"
            )));
        }
        let gv_int = numer
            .to_string()
            .parse::<Integer>()
            .map_err(|()| Error::InvalidInput("GV integer conversion failed".into()))?;
        out.push((
            semigroup.elements.column(element_idx).as_slice().to_vec(),
            gv_int,
        ));
    }
    let qn_trace = raw_qn_trace
        .into_iter()
        .map(|poly| CygvQnTracePolynomial {
            element_index: poly.element_index,
            degree: poly.degree,
            element: poly.element,
            terms: poly
                .terms
                .into_iter()
                .map(|term| CygvQnTraceTerm {
                    monomial_index: term.monomial_index,
                    exponent: term.exponent,
                    coefficient: term.coefficient.to_string(),
                })
                .collect(),
            li2_terms: poly
                .li2_terms
                .into_iter()
                .map(|term| CygvQnTraceTerm {
                    monomial_index: term.monomial_index,
                    exponent: term.exponent,
                    coefficient: term.coefficient.to_string(),
                })
                .collect(),
        })
        .collect();
    let gv_coefficient_trace = convert_cygv_gv_coefficient_trace(raw_gv_coefficient_trace);

    Ok(GvInvariantsWithQnTrace {
        invariants: out,
        qn_trace,
        gv_coefficient_trace,
    })
}

fn cygv_thread_count_from_env() -> usize {
    env::var("CYRUS_GV_THREADS")
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .map_or_else(
            || std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get),
            |n| {
                if n == 0 { 1 } else { n as usize }
            },
        )
}

fn cygv_pool_size_from_env() -> usize {
    env::var("CYRUS_GV_POOL_SIZE")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(1000)
}

fn cygv_q_matrix(q_matrix: &[Vec<i64>], dim: usize) -> Result<DMatrix<i32>> {
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
    Ok(q.transpose())
}

fn cygv_nef_partition(
    nef_partition: &[Vec<usize>],
    q_row_count: usize,
) -> Result<Vec<DVector<i32>>> {
    if nef_partition.is_empty() {
        return Err(Error::InvalidInput("nef partition is empty".into()));
    }
    let mut seen = HashSet::new();
    let mut out = Vec::with_capacity(nef_partition.len());
    for part in nef_partition {
        if part.is_empty() {
            return Err(Error::InvalidInput(
                "nef partition contains an empty part".into(),
            ));
        }
        let mut converted = Vec::with_capacity(part.len());
        for &idx in part {
            if idx >= q_row_count {
                return Err(Error::InvalidInput(format!(
                    "nef partition index {idx} is out of range for q row count {q_row_count}"
                )));
            }
            if !seen.insert(idx) {
                return Err(Error::InvalidInput(format!(
                    "nef partition index {idx} appears in more than one part"
                )));
            }
            converted.push(i32::try_from(idx).map_err(|_| {
                Error::InvalidInput("nef partition index does not fit in i32".into())
            })?);
        }
        out.push(DVector::from_column_slice(&converted));
    }
    Ok(out)
}

fn cygv_intnums_map(intnums: &Intersection) -> Result<HashMap<(usize, usize, usize), i32>> {
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
    Ok(intnums_map)
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
    let faces = polytope.faces_4d_for_points(points)?;
    Ok((faces.facets, faces.twofaces))
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
    use std::collections::{BTreeMap, BTreeSet, HashMap, VecDeque};
    use std::f64::consts::PI;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::{
        BoundedCurveDecompositionIndex, CkyzExpCoefficientCache, CkyzIndexedSeries,
        CkyzLocalIntersectionTerm, CkyzLocalSurfaceIdentification, CkyzLocalSurfaceKind,
        CkyzMonomialDomain, CurveDecompositionTerm, CurvePruningStrategy, CygvQnTraceTerm,
        GvCachePolicy, GvLatticeAugmentation, LocalToricCircuitKind, LocalToricCoordinate2D,
        NilpotentRayCandidate, NilpotentRayDegreeSlice, NilpotentRaySliceDistance,
        OriginCircuitCurveWitness, OriginCircuitRelationPoint, SupportingMoriFaceLpSearchOptions,
        ToricCurveCandidate, certify_supporting_mori_face_by_exact_kernel,
        certify_supporting_mori_face_by_lp_search, check_extremal_mori_ray_separator,
        check_supporting_mori_face_normal, ckyz_cover_closed_target_degrees,
        ckyz_cygv_previous_qn_level_count, ckyz_expalpha_power_caches_domain, ckyz_grading_degree,
        ckyz_grading_vector_from_cover_weights, ckyz_indexed_alpha_series,
        ckyz_indexed_q_degree_series_with_previous_cache_in_z_domain,
        ckyz_local_domain_profile_for_degrees,
        ckyz_local_domain_profile_for_degrees_with_causal_domain,
        ckyz_local_surface_causal_domain_spec, ckyz_local_surface_cover_weight_coefficients,
        ckyz_local_surface_domain_profile_for_multiples, ckyz_local_surface_target_degrees,
        ckyz_log_period_support_indices_domain, ckyz_multi_degrees_box,
        ckyz_observed_support_domain_for_degrees, ckyz_predicted_support_domain_for_degrees,
        ckyz_q_degree_li2_coefficient_in_z_domain,
        ckyz_q_degree_series_from_expalpha_powers_in_z_domain, ckyz_q_degree_series_in_z_domain,
        ckyz_scaled_alpha_terms, ckyz_second_log_period_series_for_pair_domain,
        ckyz_second_log_period_support_indices_for_pair_domain, ckyz_series_compose,
        ckyz_series_compose_domain, ckyz_series_exp, ckyz_series_exp_domain,
        ckyz_series_li2_domain, ckyz_series_mul, ckyz_series_mul_domain,
        ckyz_series_support_indices, ckyz_sort_degrees_for_extraction_with_grading,
        ckyz_support_exp_domain, ckyz_support_exp_domain_by_powers, ckyz_total_degree,
        ckyz_z_residual_coefficient_work_profile_for_degrees, ckyz_z_residual_dependency_degrees,
        classify_nilpotent_rays_from_two_pass_divergence_checks,
        classify_nop_rays_from_finite_gv_table, compute_ambient_one_dimensional_ray_gv_series,
        compute_ckyz_flat_prepotential_period_corrections, compute_ckyz_inverse_mirror_map,
        compute_ckyz_local_gv_invariants, compute_ckyz_local_gv_invariants_for_degrees,
        compute_ckyz_local_gv_invariants_for_degrees_with_causal_domain,
        compute_ckyz_local_gv_invariants_for_degrees_with_predicted_support_domain,
        compute_ckyz_local_instanton_potential_corrections,
        compute_ckyz_local_instanton_potential_corrections_domain,
        compute_ckyz_local_instanton_potential_z_domain,
        compute_ckyz_local_prepotential_period_corrections,
        compute_ckyz_local_surface_gv_invariants_for_multiples_with_causal_domain,
        compute_ckyz_log_period_corrections, compute_ckyz_log_period_corrections_domain,
        compute_grading_vector, compute_gv_invariants_inner,
        compute_gv_invariants_with_explicit_semigroup,
        compute_gv_invariants_with_explicit_semigroup_and_nef_partition,
        compute_gv_invariants_with_explicit_semigroup_qn_trace,
        compute_gv_invariants_with_provided_generators,
        compute_gv_invariants_with_provided_generators_and_nef_partition,
        compute_gv_invariants_with_provided_generators_qn_trace,
        compute_gw_coefficient_trace_with_explicit_semigroup,
        compute_gw_coefficient_trace_with_provided_generators,
        compute_local_p2_genus_zero_gv_series, compute_local_toric_circuit_gv_series,
        compute_one_dimensional_ray_gv_series, compute_ray_gv_series_with_provided_generators,
        curve_in_rational_row_span, curve_row_span_rank, curve_volume_in_divisor_basis,
        cygv_pair_reduced_seed_generators, detect_apparent_nilpotent_ray_from_gv_multiples,
        detect_apparent_nilpotent_rays_from_gv_table, diagnose_affine_toric_circuit,
        diagnose_supporting_mori_face_by_lp_search, dump_mori_rays_cdd, exact_i64_dot_checked,
        extract_ckyz_local_gv_invariants_from_potential,
        extract_ckyz_local_gv_invariants_from_potential_for_degrees,
        extract_ckyz_local_gv_invariants_from_z_potential_for_degrees,
        find_extremal_mori_ray_separator, find_extremal_mori_ray_separator_by_lp_search,
        find_pair_decomposition, find_semigroup_decomposition,
        finite_cutoff_gv_charges_excluding_primitive_rays, finite_gv_nonzero_degree_slice_points,
        gv_divisor_basis_data, gv_lattice_augmentation_grading, gv_lattice_search_request,
        intersection_in_divisor_basis, intersection_in_matrix_divisor_basis, load_grading_cache,
        local_p2_inverse_mirror_map, local_p2_mirror_correction,
        map_basis_gv_invariants_to_ambient, nilpotent_ray_degree_slice_for_cutoff_fraction,
        nilpotent_ray_divergence_check_from_slice_distances,
        nilpotent_ray_divergence_check_with_explicit_slice_lattices,
        nilpotent_ray_lll_reduced_slice_distance, nilpotent_ray_slice_comparison_points,
        origin_circuit_diagnostic_from_class_and_witnesses,
        partition_finite_cutoff_gv_charges_by_nilpotence, potent_ray_convergence,
        potent_ray_log_xi_terms, project_ambient_curve_to_basis,
        project_ambient_curve_to_basis_matrix, project_mori_cone_cap_rays_for_divisor_basis,
        project_mori_cone_cap_rays_to_basis, project_mori_cone_cap_rays_to_basis_matrix,
        prune_decomposable_curve_candidates, rank_two_local_charge_model,
        rank_two_local_support_signature, remove_pair_decomposable_curve_candidates,
        remove_semigroup_decomposable_curve_candidates, resolved_conifold_origin_circuit_gv,
        subcutoff_toric_curve_candidates, supporting_mori_face_for_curve_from_normal,
        supporting_mori_face_from_normal, write_grading_cache,
    };
    use crate::Intersection;
    use crate::lattice::Point;
    use crate::{DivisorBasis, f64_finite, f64_pos};
    use malachite::Integer;
    use malachite::Rational;
    use nalgebra::{DMatrix, RowDVector};

    #[test]
    fn cygv_pair_seed_reduction_mirrors_private_pair_pruning() {
        let reduced = cygv_pair_reduced_seed_generators(&[
            vec![0, 0],
            vec![1, 0],
            vec![0, 1],
            vec![1, 1],
            vec![2, 0],
            vec![2, 1],
            vec![1, -1],
        ])
        .unwrap();
        assert_eq!(reduced, vec![vec![0, 1], vec![1, -1]]);
    }

    #[test]
    fn cygv_pair_seed_reduction_rejects_bad_inputs() {
        assert!(cygv_pair_reduced_seed_generators(&[]).is_err());
        assert!(cygv_pair_reduced_seed_generators(&[vec![]]).is_err());
        assert!(cygv_pair_reduced_seed_generators(&[vec![1, 0], vec![1]]).is_err());
        assert!(
            cygv_pair_reduced_seed_generators(&[vec![i64::MAX, 0], vec![1, 0], vec![i64::MIN, 0],])
                .is_err()
        );
    }

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
    fn provided_generator_gv_path_skips_lattice_augmentation() {
        let err = compute_gv_invariants_with_provided_generators(
            &[vec![1]],
            &[1],
            &[],
            &Intersection::new(1),
            None,
            Some(1),
        )
        .unwrap_err();

        assert!(err.to_string().contains("q_matrix is empty"));
    }

    #[test]
    fn provided_generator_gv_path_rejects_empty_generators_without_panicking() {
        let err = compute_gv_invariants_with_provided_generators(
            &[],
            &[1],
            &[vec![1]],
            &Intersection::new(1),
            None,
            Some(1),
        )
        .unwrap_err();

        assert!(err.to_string().contains("at least one generator"));
    }

    #[test]
    fn cygv_wrappers_require_exactly_one_truncation_control() {
        let mut intnums = Intersection::new(1);
        set_intersection_i64(&mut intnums, 0, 0, 0, 5);

        let neither = compute_gv_invariants_with_provided_generators(
            &[vec![1]],
            &[1],
            &[vec![1, 1, 1, 1, 1]],
            &intnums,
            None,
            None,
        )
        .unwrap_err();
        let both = compute_gv_invariants_with_provided_generators(
            &[vec![1]],
            &[1],
            &[vec![1, 1, 1, 1, 1]],
            &intnums,
            Some(10),
            Some(1),
        )
        .unwrap_err();

        assert!(
            neither
                .to_string()
                .contains("Exactly one of min_points or max_deg")
        );
        assert!(
            both.to_string()
                .contains("Exactly one of min_points or max_deg")
        );
    }

    #[test]
    fn provided_generator_gv_path_converts_cygv_semigroup_errors_to_result() {
        let mut intnums = Intersection::new(1);
        set_intersection_i64(&mut intnums, 0, 0, 0, 1);

        let err = compute_gv_invariants_with_provided_generators(
            &[vec![1]],
            &[-1],
            &[vec![1]],
            &intnums,
            None,
            Some(1),
        )
        .unwrap_err();

        assert!(
            err.to_string()
                .contains("cygv semigroup construction failed")
        );
    }

    #[test]
    fn provided_generator_complete_intersection_path_validates_nef_partition() {
        let mut intnums = Intersection::new(1);
        set_intersection_i64(&mut intnums, 0, 0, 0, 1);

        let err = compute_gv_invariants_with_provided_generators_and_nef_partition(
            &[vec![1]],
            &[1],
            &[vec![1, 1, 1, 1, 1, 1]],
            &[vec![0, 1], vec![1, 2]],
            &intnums,
            None,
            Some(1),
        )
        .unwrap_err();

        assert!(
            err.to_string()
                .contains("nef partition index 1 appears in more than one part")
        );
    }

    #[test]
    fn explicit_complete_intersection_path_passes_nef_partition_to_cygv() {
        let mut intnums = Intersection::new(1);
        set_intersection_i64(&mut intnums, 0, 0, 0, 1);

        let err = compute_gv_invariants_with_explicit_semigroup_and_nef_partition(
            &[vec![0], vec![1]],
            &[1],
            &[vec![1, 1, 1, 1, 1, 1]],
            &[vec![0], vec![1], vec![2]],
            &intnums,
        )
        .unwrap_err();

        assert!(
            err.to_string().contains("cygv fundamental period failed"),
            "expected a cygv-stage error after nef partition conversion, got {err}"
        );
    }

    #[test]
    fn explicit_semigroup_gv_path_accepts_exact_truncation_without_degree_bound() {
        let err = compute_gv_invariants_with_explicit_semigroup(
            &[vec![1]],
            &[1],
            &[],
            &Intersection::new(1),
        )
        .unwrap_err();

        assert!(err.to_string().contains("q_matrix is empty"));
    }

    #[test]
    fn ckyz_small_monomial_domain_precomputes_addition_table() {
        let domain = CkyzMonomialDomain::componentwise_box(&[2, 2]).unwrap();
        assert!(
            domain.addition_indices.is_some(),
            "small CKYZ domains should precompute product indices"
        );
        assert!(
            domain.addition_pairs_by_lhs.is_some(),
            "small CKYZ domains should precompute sparse valid product pairs"
        );

        let lhs = BTreeMap::from([(vec![1, 0], Rational::from(2))]);
        let rhs = BTreeMap::from([(vec![0, 2], Rational::from(3))]);
        let product = ckyz_series_mul_domain(&lhs, &rhs, &domain).unwrap();

        assert_eq!(product, BTreeMap::from([(vec![1, 2], Rational::from(6))]));
    }

    #[test]
    fn ckyz_domain_multiplication_is_symmetric_with_indexed_sparse_terms() {
        let domain = CkyzMonomialDomain::componentwise_box(&[2, 1]).unwrap();
        let lhs = BTreeMap::from([
            (vec![1, 0], Rational::from(2)),
            (vec![0, 1], Rational::from(5)),
        ]);
        let rhs = BTreeMap::from([(vec![1, 0], Rational::from(3))]);

        let lhs_rhs = ckyz_series_mul_domain(&lhs, &rhs, &domain).unwrap();
        let rhs_lhs = ckyz_series_mul_domain(&rhs, &lhs, &domain).unwrap();

        assert_eq!(lhs_rhs, rhs_lhs);
        assert_eq!(
            lhs_rhs,
            BTreeMap::from([
                (vec![2, 0], Rational::from(6)),
                (vec![1, 1], Rational::from(15)),
            ])
        );
    }

    #[test]
    fn ckyz_power_cache_fills_empty_higher_powers() {
        let series = BTreeMap::from([(vec![1], Rational::from(2))]);

        let powers = super::ckyz_series_power_cache(&series, 4, 1, 2).unwrap();

        assert_eq!(powers.len(), 5);
        assert_eq!(powers[0], BTreeMap::from([(vec![0], Rational::from(1))]));
        assert_eq!(powers[1], BTreeMap::from([(vec![1], Rational::from(2))]));
        assert_eq!(powers[2], BTreeMap::from([(vec![2], Rational::from(4))]));
        assert!(powers[3].is_empty());
        assert!(powers[4].is_empty());
    }

    #[test]
    fn ckyz_large_monomial_domain_uses_checked_addition_fallback() {
        let domain = CkyzMonomialDomain::componentwise_box(&[5_000]).unwrap();
        assert!(
            domain.addition_indices.is_none(),
            "large CKYZ domains should avoid quadratic addition tables"
        );
        assert!(
            domain.addition_pairs_by_lhs.is_none(),
            "large CKYZ domains should avoid quadratic addition-pair tables"
        );

        let lhs = BTreeMap::from([(vec![125], Rational::from(2))]);
        let rhs = BTreeMap::from([(vec![175], Rational::from(3))]);
        let product = ckyz_series_mul_domain(&lhs, &rhs, &domain).unwrap();

        assert_eq!(product, BTreeMap::from([(vec![300], Rational::from(6))]));
    }

    #[test]
    fn ckyz_polygon5_four_multiple_target_downset_precomputes_addition_table() {
        let target_degrees = [
            vec![4, 3, 2],
            vec![8, 6, 4],
            vec![12, 9, 6],
            vec![16, 12, 8],
        ];
        let extraction_degrees = ckyz_cover_closed_target_degrees(&target_degrees).unwrap();
        let domain = CkyzMonomialDomain::target_downset(&extraction_degrees, 3).unwrap();

        assert_eq!(domain.degrees.len(), 1_989);
        assert!(
            domain.addition_indices.is_some(),
            "polygon-5 N=4 target-downset domains should keep indexed products"
        );
        let pair_count = domain
            .addition_pairs_by_lhs
            .as_ref()
            .expect("indexed products should have sparse pairs")
            .iter()
            .map(Vec::len)
            .sum::<usize>();
        assert!(
            pair_count < domain.degrees.len() * domain.degrees.len(),
            "polygon-5 N=4 target-downset products should be sparse"
        );
    }

    #[test]
    fn ckyz_target_downset_keeps_only_past_monomials_for_requested_degrees() {
        let domain = CkyzMonomialDomain::target_downset(&[vec![3, 0], vec![0, 3]], 2).unwrap();

        assert_eq!(domain.degrees.len(), 7);
        assert!(domain.contains(&[0, 0]));
        assert!(domain.contains(&[3, 0]));
        assert!(domain.contains(&[0, 3]));
        assert!(!domain.contains(&[1, 1]));
        assert_eq!(domain.max_total_degree, 3);
    }

    #[test]
    fn ckyz_explicit_domain_drops_products_outside_monomial_map() {
        let domain =
            CkyzMonomialDomain::from_degrees(2, [vec![0, 0], vec![1, 0], vec![0, 1], vec![2, 0]])
                .unwrap();
        let lhs = BTreeMap::from([
            (vec![1, 0], Rational::from(2)),
            (vec![0, 1], Rational::from(5)),
        ]);
        let rhs = BTreeMap::from([(vec![1, 0], Rational::from(3))]);

        let product = ckyz_series_mul_domain(&lhs, &rhs, &domain).unwrap();

        assert_eq!(product, BTreeMap::from([(vec![2, 0], Rational::from(6))]));
        assert!(
            !product.contains_key(&vec![1, 1]),
            "explicit CKYZ domains must mirror cygv monomial_map semantics by dropping absent sums"
        );
    }

    #[test]
    fn ckyz_exponential_uses_minimum_nonzero_degree_bound() {
        let series = BTreeMap::from([(vec![2, 0], Rational::from(2))]);

        let exponential = super::ckyz_series_exp(&series, 2, 5).unwrap();

        assert_eq!(
            exponential,
            BTreeMap::from([
                (vec![0, 0], Rational::from(1)),
                (vec![2, 0], Rational::from(2)),
                (vec![4, 0], Rational::from(2)),
            ])
        );
    }

    #[test]
    fn affine_toric_circuit_detects_local_p2_triangle() {
        let points = vec![
            Point::new(vec![0, 1, -3, 6]),
            Point::new(vec![-2, -1, -4, 5]),
            Point::new(vec![-1, 0, -3, 5]),
            Point::new(vec![-1, 0, -2, 4]),
        ];

        let diagnostic = diagnose_affine_toric_circuit(&[1, 1, -3, 1], &points)
            .unwrap()
            .expect("local P2 row is an affine circuit");

        assert_eq!(diagnostic.coefficient_sum, 0);
        assert_eq!(diagnostic.coordinate_sum, vec![0, 0, 0, 0]);
        assert_eq!(diagnostic.affine_rank, 2);
        assert_eq!(
            diagnostic.coefficient_counts,
            BTreeMap::from([(-3, 1), (1, 3)])
        );
        assert_eq!(diagnostic.local_charge_basis, vec![vec![1, 1, -3, 1]]);
        assert!(
            diagnostic.local_coordinates_2d.is_some(),
            "rank-two local P2 support should expose local coordinates"
        );
        assert_eq!(
            diagnostic.kind,
            Some(LocalToricCircuitKind::LocalP2Triangle {
                interior_point: 2,
                vertex_points: vec![0, 1, 3],
                interior_coefficient: -3,
                vertex_coefficient: 1,
                local_coordinates: vec![
                    LocalToricCoordinate2D {
                        point_index: 0,
                        coordinates: [1, 0],
                    },
                    LocalToricCoordinate2D {
                        point_index: 1,
                        coordinates: [0, 1],
                    },
                    LocalToricCoordinate2D {
                        point_index: 2,
                        coordinates: [0, 0],
                    },
                    LocalToricCoordinate2D {
                        point_index: 3,
                        coordinates: [-1, -1],
                    },
                ],
            })
        );
    }

    #[test]
    fn local_p2_gv_series_is_computed_from_mirror_map() {
        let gvs = compute_local_p2_genus_zero_gv_series(10).unwrap();

        assert_eq!(
            gvs,
            vec![
                Integer::from(3),
                Integer::from(-6),
                Integer::from(27),
                Integer::from(-192),
                Integer::from(1695),
                Integer::from(-17064),
                Integer::from(188454),
                Integer::from(-2228160),
                Integer::from(27748899),
                Integer::from(-360012150),
            ]
        );
    }

    #[test]
    fn ckyz_log_period_corrections_match_local_p2_mirror_map() {
        let corrections = compute_ckyz_log_period_corrections(&[vec![-3, 1, 1, 1]], 6).unwrap();
        let local_p2 = local_p2_mirror_correction(6);

        for degree in 1..=6 {
            assert_eq!(
                corrections[0].get(&vec![degree]),
                Some(&local_p2[degree]),
                "CKYZ logarithmic period should reproduce local P2 mirror correction at degree {degree}"
            );
        }
    }

    #[test]
    fn ckyz_log_period_corrections_compute_f0_coupled_mirror_map_terms() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
        let corrections = compute_ckyz_log_period_corrections(&relations, 2).unwrap();

        let mut expected = BTreeMap::new();
        expected.insert(vec![1, 0], Rational::from(2));
        expected.insert(vec![0, 1], Rational::from(2));
        expected.insert(vec![2, 0], Rational::from(3));
        expected.insert(vec![1, 1], Rational::from(12));
        expected.insert(vec![0, 2], Rational::from(3));

        assert_eq!(corrections[0], expected);
        assert_eq!(corrections[1], expected);
    }

    #[test]
    fn ckyz_log_period_corrections_compute_f1_source_terms() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-1, 0, 1, -1, 1]];
        let corrections = compute_ckyz_log_period_corrections(&relations, 2).unwrap();

        let mut expected_first = BTreeMap::new();
        expected_first.insert(vec![1, 0], Rational::from(2));
        expected_first.insert(vec![2, 0], Rational::from(3));
        expected_first.insert(vec![1, 1], Rational::from(-4));

        let mut expected_second = BTreeMap::new();
        expected_second.insert(vec![1, 0], Rational::from(1));
        expected_second.insert(vec![2, 0], Rational::from_signeds(3, 2));
        expected_second.insert(vec![1, 1], Rational::from(-2));

        assert_eq!(corrections[0], expected_first);
        assert_eq!(corrections[1], expected_second);
    }

    #[test]
    fn ckyz_local_prepotential_period_computes_p2_double_log_source_term() {
        let corrections = compute_ckyz_local_prepotential_period_corrections(
            &[vec![-3, 1, 1, 1]],
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 0,
                coefficient: 1,
            }],
            1,
        )
        .unwrap();

        assert_eq!(corrections.get(&vec![1]), Some(&Rational::from(-18)));
    }

    #[test]
    fn ckyz_local_prepotential_period_computes_f0_source_terms() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
        let corrections = compute_ckyz_local_prepotential_period_corrections(
            &relations,
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 1,
                coefficient: 1,
            }],
            2,
        )
        .unwrap();

        let mut expected = BTreeMap::new();
        expected.insert(vec![1, 0], Rational::from(4));
        expected.insert(vec![0, 1], Rational::from(4));
        expected.insert(vec![2, 0], Rational::from(13));
        expected.insert(vec![1, 1], Rational::from(40));
        expected.insert(vec![0, 2], Rational::from(13));

        assert_eq!(corrections, expected);
    }

    #[test]
    fn ckyz_inverse_mirror_map_matches_local_p2_inverse() {
        let corrections = compute_ckyz_log_period_corrections(&[vec![-3, 1, 1, 1]], 5).unwrap();
        let z_of_q = compute_ckyz_inverse_mirror_map(&corrections, 5).unwrap();
        let local_p2 = local_p2_inverse_mirror_map(&local_p2_mirror_correction(5), 5);

        for (degree, expected) in local_p2.iter().enumerate().skip(1) {
            assert_eq!(
                z_of_q[0].get(&vec![degree]),
                Some(expected),
                "CKYZ inverse mirror map should match local P2 at degree {degree}"
            );
        }
    }

    #[test]
    fn ckyz_inverse_mirror_map_zero_corrections_is_identity() {
        let corrections = vec![BTreeMap::new(), BTreeMap::new()];
        let z_of_q = compute_ckyz_inverse_mirror_map(&corrections, 5).unwrap();

        assert_eq!(z_of_q[0], super::ckyz_coordinate_series(2, 0, 5));
        assert_eq!(z_of_q[1], super::ckyz_coordinate_series(2, 1, 5));
    }

    #[test]
    fn ckyz_inverse_mirror_map_computes_f0_coupled_terms() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
        let corrections = compute_ckyz_log_period_corrections(&relations, 2).unwrap();
        let z_of_q = compute_ckyz_inverse_mirror_map(&corrections, 2).unwrap();

        let mut expected_first = BTreeMap::new();
        expected_first.insert(vec![1, 0], Rational::from(1));
        expected_first.insert(vec![2, 0], Rational::from(-2));
        expected_first.insert(vec![1, 1], Rational::from(-2));

        let mut expected_second = BTreeMap::new();
        expected_second.insert(vec![0, 1], Rational::from(1));
        expected_second.insert(vec![1, 1], Rational::from(-2));
        expected_second.insert(vec![0, 2], Rational::from(-2));

        assert_eq!(z_of_q[0], expected_first);
        assert_eq!(z_of_q[1], expected_second);
    }

    #[test]
    fn ckyz_flat_prepotential_period_substitutes_f0_mirror_map() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
        let corrections = compute_ckyz_flat_prepotential_period_corrections(
            &relations,
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 1,
                coefficient: 1,
            }],
            2,
        )
        .unwrap();

        let mut expected = BTreeMap::new();
        expected.insert(vec![1, 0], Rational::from(4));
        expected.insert(vec![0, 1], Rational::from(4));
        expected.insert(vec![2, 0], Rational::from(5));
        expected.insert(vec![1, 1], Rational::from(24));
        expected.insert(vec![0, 2], Rational::from(5));

        assert_eq!(corrections, expected);
    }

    #[test]
    fn ckyz_local_instanton_potential_removes_alpha_product_terms() {
        let potential = compute_ckyz_local_instanton_potential_corrections(
            &[vec![-3, 1, 1, 1]],
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 0,
                coefficient: 1,
            }],
            3,
        )
        .unwrap();

        assert_eq!(potential.get(&vec![1]), Some(&Rational::from(-9)));
        assert_eq!(
            potential.get(&vec![2]),
            Some(&Rational::from_signeds(135, 4))
        );
        assert_eq!(potential.get(&vec![3]), Some(&Rational::from(-244)));
    }

    #[test]
    fn ckyz_local_gv_extraction_matches_local_p2_table() {
        let gvs = compute_ckyz_local_gv_invariants(
            &[vec![-3, 1, 1, 1]],
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 0,
                coefficient: 1,
            }],
            &[3],
            6,
        )
        .unwrap();

        let expected = BTreeMap::from([
            (vec![1], Integer::from(3)),
            (vec![2], Integer::from(-6)),
            (vec![3], Integer::from(27)),
            (vec![4], Integer::from(-192)),
            (vec![5], Integer::from(1695)),
            (vec![6], Integer::from(-17064)),
        ]);

        assert_eq!(gvs, expected);
    }

    #[test]
    fn ckyz_targeted_local_gv_extraction_matches_local_p2_table() {
        let target_degrees = (1..=10).map(|degree| vec![degree]).collect::<Vec<_>>();
        let gvs = compute_ckyz_local_gv_invariants_for_degrees(
            &[vec![-3, 1, 1, 1]],
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 0,
                coefficient: 1,
            }],
            &[3],
            &target_degrees,
        )
        .unwrap();

        let expected = BTreeMap::from([
            (vec![1], Integer::from(3)),
            (vec![2], Integer::from(-6)),
            (vec![3], Integer::from(27)),
            (vec![4], Integer::from(-192)),
            (vec![5], Integer::from(1695)),
            (vec![6], Integer::from(-17064)),
            (vec![7], Integer::from(188454)),
            (vec![8], Integer::from(-2228160)),
            (vec![9], Integer::from(27748899)),
            (vec![10], Integer::from(-360012150)),
        ]);

        assert_eq!(gvs, expected);
    }

    #[test]
    fn ckyz_targeted_gv_extraction_is_independent_of_requested_degree_order() {
        let potential = BTreeMap::from([
            (vec![1], Rational::from(-2)),
            (vec![2], Rational::from_signeds(-13, 2)),
        ]);
        let expected = BTreeMap::from([(vec![1], Integer::from(2)), (vec![2], Integer::from(3))]);

        let sorted = extract_ckyz_local_gv_invariants_from_potential_for_degrees(
            &potential,
            &[1],
            &[vec![1], vec![2]],
        )
        .unwrap();
        let shuffled = extract_ckyz_local_gv_invariants_from_potential_for_degrees(
            &potential,
            &[1],
            &[vec![2], vec![1]],
        )
        .unwrap();

        assert_eq!(sorted, expected);
        assert_eq!(shuffled, expected);
    }

    #[test]
    fn ckyz_causal_domain_local_gv_extraction_matches_local_p2_table() {
        let target_degrees = (1..=10).map(|degree| vec![degree]).collect::<Vec<_>>();
        let gvs = compute_ckyz_local_gv_invariants_for_degrees_with_causal_domain(
            &[vec![-3, 1, 1, 1]],
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 0,
                coefficient: 1,
            }],
            &[3],
            &target_degrees,
            &[vec![1]],
            &[1],
        )
        .unwrap();

        let expected = BTreeMap::from([
            (vec![1], Integer::from(3)),
            (vec![2], Integer::from(-6)),
            (vec![3], Integer::from(27)),
            (vec![4], Integer::from(-192)),
            (vec![5], Integer::from(1695)),
            (vec![6], Integer::from(-17064)),
            (vec![7], Integer::from(188454)),
            (vec![8], Integer::from(-2228160)),
            (vec![9], Integer::from(27748899)),
            (vec![10], Integer::from(-360012150)),
        ]);

        assert_eq!(gvs, expected);
    }

    #[test]
    fn ckyz_local_gv_extraction_matches_f0_table_start() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
        let gvs = compute_ckyz_local_gv_invariants(
            &relations,
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 1,
                coefficient: 1,
            }],
            &[2, 2],
            4,
        )
        .unwrap();

        let expected_nonzero = BTreeMap::from([
            (vec![1, 0], Integer::from(-2)),
            (vec![0, 1], Integer::from(-2)),
            (vec![1, 1], Integer::from(-4)),
            (vec![2, 1], Integer::from(-6)),
            (vec![1, 2], Integer::from(-6)),
            (vec![3, 1], Integer::from(-8)),
            (vec![1, 3], Integer::from(-8)),
            (vec![2, 2], Integer::from(-32)),
        ]);

        for (degree, expected) in expected_nonzero {
            assert_eq!(gvs.get(&degree), Some(&expected));
        }
        for zero_degree in [vec![2, 0], vec![0, 2], vec![3, 0], vec![0, 3]] {
            assert!(
                !gvs.contains_key(&zero_degree),
                "zero F0 invariant should be omitted at {zero_degree:?}"
            );
        }
    }

    #[test]
    fn ckyz_targeted_local_gv_extraction_matches_full_f0_ray() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
        let target_degrees = [vec![1, 1], vec![2, 2], vec![3, 3]];
        let targeted = compute_ckyz_local_gv_invariants_for_degrees(
            &relations,
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 1,
                coefficient: 1,
            }],
            &[2, 2],
            &target_degrees,
        )
        .unwrap();
        let full = compute_ckyz_local_gv_invariants(
            &relations,
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 1,
                coefficient: 1,
            }],
            &[2, 2],
            6,
        )
        .unwrap();

        for degree in target_degrees {
            assert_eq!(
                targeted.get(&degree),
                full.get(&degree),
                "targeted F0 extraction should match full extraction at {degree:?}"
            );
        }
    }

    #[test]
    fn ckyz_causal_domain_matches_target_downset_for_full_f0_ray() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
        let target_degrees = [vec![1, 1], vec![2, 2], vec![3, 3]];
        let targeted = compute_ckyz_local_gv_invariants_for_degrees(
            &relations,
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 1,
                coefficient: 1,
            }],
            &[2, 2],
            &target_degrees,
        )
        .unwrap();
        let causal = compute_ckyz_local_gv_invariants_for_degrees_with_causal_domain(
            &relations,
            &[CkyzLocalIntersectionTerm {
                first: 0,
                second: 1,
                coefficient: 1,
            }],
            &[2, 2],
            &target_degrees,
            &[vec![1, 0], vec![0, 1]],
            &[1, 1],
        )
        .unwrap();

        assert_eq!(causal, targeted);
    }

    #[test]
    fn ckyz_observed_support_domain_recomputes_targeted_f0_ray() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
        let local_intersection_terms = [CkyzLocalIntersectionTerm {
            first: 0,
            second: 1,
            coefficient: 1,
        }];
        let cover_weights = [2, 2];
        let target_degrees = [vec![3, 2], vec![6, 4], vec![9, 6]];

        let broad = compute_ckyz_local_gv_invariants_for_degrees(
            &relations,
            &local_intersection_terms,
            &cover_weights,
            &target_degrees,
        )
        .unwrap();
        let extraction_degrees = ckyz_cover_closed_target_degrees(&target_degrees).unwrap();
        let broad_domain = CkyzMonomialDomain::target_downset(&extraction_degrees, 2).unwrap();
        let observed_domain = ckyz_observed_support_domain_for_degrees(
            &relations,
            &local_intersection_terms,
            &target_degrees,
        )
        .unwrap();

        assert!(
            observed_domain.degrees.len() <= broad_domain.degrees.len(),
            "observed support domain must not add terms outside the componentwise target downset"
        );

        let potential = compute_ckyz_local_instanton_potential_corrections_domain(
            &relations,
            &local_intersection_terms,
            &observed_domain,
        )
        .unwrap();
        let mut observed = extract_ckyz_local_gv_invariants_from_potential_for_degrees(
            &potential,
            &cover_weights,
            &extraction_degrees,
        )
        .unwrap();
        observed.retain(|degree, _| target_degrees.iter().any(|target| target == degree));

        assert_eq!(observed, broad);
    }

    #[test]
    fn ckyz_predicted_support_domain_covers_observed_f0_ray_support() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
        let local_intersection_terms = [CkyzLocalIntersectionTerm {
            first: 0,
            second: 1,
            coefficient: 1,
        }];
        let cover_weights = [2, 2];
        let target_degrees = [vec![3, 2], vec![6, 4], vec![9, 6]];
        let extraction_degrees = ckyz_cover_closed_target_degrees(&target_degrees).unwrap();

        let observed_domain = ckyz_observed_support_domain_for_degrees(
            &relations,
            &local_intersection_terms,
            &target_degrees,
        )
        .unwrap();
        let predicted_domain = ckyz_predicted_support_domain_for_degrees(
            &relations,
            &local_intersection_terms,
            &target_degrees,
        )
        .unwrap();
        assert_eq!(observed_domain.degrees.len(), 70);
        assert_eq!(predicted_domain.degrees, observed_domain.degrees);
        assert!(
            observed_domain
                .degrees
                .iter()
                .all(|degree| predicted_domain.contains(degree)),
            "support-predicted CKYZ domain must cover the observed broad-computation support"
        );

        let potential = compute_ckyz_local_instanton_potential_corrections_domain(
            &relations,
            &local_intersection_terms,
            &predicted_domain,
        )
        .unwrap();
        let mut predicted = extract_ckyz_local_gv_invariants_from_potential_for_degrees(
            &potential,
            &cover_weights,
            &extraction_degrees,
        )
        .unwrap();
        predicted.retain(|degree, _| target_degrees.iter().any(|target| target == degree));

        let broad = compute_ckyz_local_gv_invariants_for_degrees(
            &relations,
            &local_intersection_terms,
            &cover_weights,
            &target_degrees,
        )
        .unwrap();
        assert_eq!(predicted, broad);
    }

    #[test]
    fn ckyz_predicted_support_domain_api_matches_target_downset_for_f0_ray() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
        let local_intersection_terms = [CkyzLocalIntersectionTerm {
            first: 0,
            second: 1,
            coefficient: 1,
        }];
        let cover_weights = [2, 2];
        let target_degrees = [vec![3, 2], vec![6, 4], vec![9, 6]];

        let predicted = compute_ckyz_local_gv_invariants_for_degrees_with_predicted_support_domain(
            &relations,
            &local_intersection_terms,
            &cover_weights,
            &target_degrees,
        )
        .unwrap();
        let broad = compute_ckyz_local_gv_invariants_for_degrees(
            &relations,
            &local_intersection_terms,
            &cover_weights,
            &target_degrees,
        )
        .unwrap();

        assert_eq!(predicted, broad);
    }

    #[test]
    fn ckyz_predicted_support_domain_covers_observed_polygon5_ray_support() {
        let relations = ckyz_polygon5_relations();
        let local_intersection_terms = ckyz_polygon5_intersection_terms();
        let cover_weights = [1, 1, 1];
        let target_degrees = [vec![4, 3, 2], vec![8, 6, 4]];
        let extraction_degrees = ckyz_cover_closed_target_degrees(&target_degrees).unwrap();
        let broad_domain = CkyzMonomialDomain::target_downset(&extraction_degrees, 3).unwrap();

        let observed_domain = ckyz_observed_support_domain_for_degrees(
            &relations,
            &local_intersection_terms,
            &target_degrees,
        )
        .unwrap();
        let predicted_domain = ckyz_predicted_support_domain_for_degrees(
            &relations,
            &local_intersection_terms,
            &target_degrees,
        )
        .unwrap();
        assert_eq!(broad_domain.degrees.len(), 315);
        assert_eq!(observed_domain.degrees.len(), 265);
        assert_eq!(predicted_domain.degrees, observed_domain.degrees);
        assert!(
            observed_domain
                .degrees
                .iter()
                .all(|degree| predicted_domain.contains(degree)),
            "support-predicted polygon-5 domain must cover the observed broad-computation support"
        );

        let potential = compute_ckyz_local_instanton_potential_corrections_domain(
            &relations,
            &local_intersection_terms,
            &predicted_domain,
        )
        .unwrap();
        let mut predicted = extract_ckyz_local_gv_invariants_from_potential_for_degrees(
            &potential,
            &cover_weights,
            &extraction_degrees,
        )
        .unwrap();
        predicted.retain(|degree, _| target_degrees.iter().any(|target| target == degree));

        let broad = compute_ckyz_local_gv_invariants_for_degrees(
            &relations,
            &local_intersection_terms,
            &cover_weights,
            &target_degrees,
        )
        .unwrap();
        assert_eq!(predicted, broad);
    }

    #[test]
    fn ckyz_predicted_support_domain_api_matches_target_downset_for_polygon5_ray() {
        let relations = ckyz_polygon5_relations();
        let local_intersection_terms = ckyz_polygon5_intersection_terms();
        let cover_weights = [1, 1, 1];
        let target_degrees = [vec![4, 3, 2], vec![8, 6, 4]];

        let predicted = compute_ckyz_local_gv_invariants_for_degrees_with_predicted_support_domain(
            &relations,
            &local_intersection_terms,
            &cover_weights,
            &target_degrees,
        )
        .unwrap();
        let broad = compute_ckyz_local_gv_invariants_for_degrees(
            &relations,
            &local_intersection_terms,
            &cover_weights,
            &target_degrees,
        )
        .unwrap();

        assert_eq!(predicted, broad);
    }

    #[test]
    fn ckyz_z_series_inversion_matches_flat_coordinate_extraction_for_local_models() {
        let cases = [
            (
                vec![vec![-3, 1, 1, 1]],
                vec![CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 0,
                    coefficient: 1,
                }],
                vec![3],
                (1..=6).map(|degree| vec![degree]).collect::<Vec<_>>(),
            ),
            (
                vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]],
                vec![CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 1,
                    coefficient: 1,
                }],
                vec![2, 2],
                vec![vec![1, 1], vec![2, 2], vec![3, 3]],
            ),
            (
                ckyz_polygon5_relations(),
                ckyz_polygon5_intersection_terms(),
                vec![1, 1, 1],
                vec![vec![4, 3, 2], vec![8, 6, 4]],
            ),
        ];

        for (relations, local_intersection_terms, cover_weights, target_degrees) in cases {
            let expected = compute_ckyz_local_gv_invariants_for_degrees(
                &relations,
                &local_intersection_terms,
                &cover_weights,
                &target_degrees,
            )
            .unwrap();
            let extraction_degrees = ckyz_cover_closed_target_degrees(&target_degrees).unwrap();
            let domain =
                CkyzMonomialDomain::target_downset(&extraction_degrees, relations.len()).unwrap();
            let (alpha, potential_z) = compute_ckyz_local_instanton_potential_z_domain(
                &relations,
                &local_intersection_terms,
                &domain,
            )
            .unwrap();
            let path_history_degrees =
                ckyz_z_residual_dependency_degrees(&alpha, &extraction_degrees, &domain).unwrap();
            let mut actual = extract_ckyz_local_gv_invariants_from_z_potential_for_degrees(
                &potential_z,
                &alpha,
                &cover_weights,
                &path_history_degrees,
                &domain,
            )
            .unwrap();
            actual.retain(|degree, _| target_degrees.iter().any(|target| target == degree));

            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn ckyz_z_residual_dependency_keeps_off_ray_f0_history() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]];
        let local_intersection_terms = [CkyzLocalIntersectionTerm {
            first: 0,
            second: 1,
            coefficient: 1,
        }];
        let target_degrees = [vec![1, 1]];
        let domain = CkyzMonomialDomain::target_downset(&target_degrees, 2).unwrap();
        let (alpha, _) = compute_ckyz_local_instanton_potential_z_domain(
            &relations,
            &local_intersection_terms,
            &domain,
        )
        .unwrap();

        let history = ckyz_z_residual_dependency_degrees(&alpha, &target_degrees, &domain).unwrap();

        assert!(
            history.contains(&vec![1, 0]),
            "F0 [1,1] extraction needs the lower off-ray [1,0] subtraction"
        );
        assert!(
            history.contains(&vec![0, 1]),
            "F0 [1,1] extraction needs the lower off-ray [0,1] subtraction"
        );
        assert!(history.contains(&vec![1, 1]));
    }

    #[test]
    fn ckyz_extraction_order_uses_cover_weight_grading() {
        let mut degrees = vec![vec![2, 0], vec![0, 3], vec![1, 1], vec![0, 1]];
        let grading = ckyz_grading_vector_from_cover_weights(&[2, 1]).unwrap();

        ckyz_sort_degrees_for_extraction_with_grading(&mut degrees, &grading).unwrap();

        assert_eq!(
            degrees,
            vec![vec![0, 1], vec![1, 1], vec![0, 3], vec![2, 0]]
        );
    }

    #[test]
    fn ckyz_extraction_order_rejects_nonpositive_grading_weights() {
        let err = ckyz_grading_vector_from_cover_weights(&[2, 0]).unwrap_err();

        assert!(err.to_string().contains("grading weights must be positive"));
    }

    #[test]
    fn ckyz_z_residual_coefficient_work_profile_counts_polygon5_states() {
        let relations = ckyz_polygon5_relations();
        let local_intersection_terms = ckyz_polygon5_intersection_terms();
        let target_degrees = [vec![4, 3, 2], vec![8, 6, 4]];

        let profile = ckyz_z_residual_coefficient_work_profile_for_degrees(
            &relations,
            &local_intersection_terms,
            &[1, 1, 1],
            &target_degrees,
        )
        .unwrap();

        assert_eq!(profile.rank, 3);
        assert_eq!(profile.domain_degree_count, 265);
        assert_eq!(profile.path_history_degree_count, 79);
        assert_eq!(profile.residual_pair_count, 2_887);
        assert_eq!(profile.same_grading_pair_skip_count, 194);
        assert!(profile.componentwise_pair_count < profile.residual_pair_count);
        assert!(profile.support_pair_count <= profile.componentwise_pair_count);
        assert!(profile.support_li2_delta_term_count <= profile.li2_delta_term_count);
        assert!(profile.unique_delta_count <= profile.domain_degree_count);
        assert!(profile.unique_scale_count <= profile.li2_delta_term_count);
        assert!(profile.unique_exp_state_count <= profile.li2_delta_term_count);
        assert!(profile.support_unique_exp_state_count <= profile.unique_exp_state_count);
        assert_eq!(profile.qn_history_level_count, 2);
        assert_eq!(
            profile.qn_history_candidate_hit_count + profile.qn_history_candidate_miss_count,
            profile.path_history_degree_count
        );
        assert!(profile.qn_history_candidate_hit_count > 0);
        assert!(profile.qn_history_unique_delta_count <= profile.path_history_degree_count);
        assert_eq!(
            profile.q_delta_domain_profiles.len(),
            profile.qn_history_unique_delta_count
        );
        assert!(profile.q_delta_domain_profiles.iter().all(|delta_profile| {
            delta_profile.shiftable_exp_degree_count <= profile.domain_degree_count
                && delta_profile.predecessor_closure_degree_count <= profile.domain_degree_count
                && delta_profile.shiftable_exp_degree_count
                    <= delta_profile.predecessor_closure_degree_count
        }));
        assert_eq!(
            profile
                .support_exp_state_counts_by_scale
                .iter()
                .map(|(_, count)| count)
                .sum::<usize>(),
            profile.support_unique_exp_state_count
        );
    }

    #[test]
    fn ckyz_indexed_domain_series_ops_match_total_degree_reference() {
        let rank = 2;
        let max_total_degree = 4;
        let domain = CkyzMonomialDomain::from_degrees(
            rank,
            ckyz_multi_degrees_box(&vec![max_total_degree; rank])
                .into_iter()
                .filter(|degree| ckyz_total_degree(degree).unwrap() <= max_total_degree),
        )
        .unwrap();
        let lhs = BTreeMap::from([
            (vec![1, 0], Rational::from(2)),
            (vec![0, 1], Rational::from(3)),
            (vec![2, 0], Rational::from(-1)),
        ]);
        let rhs = BTreeMap::from([
            (vec![0, 1], Rational::from(5)),
            (vec![1, 1], Rational::from(7)),
            (vec![0, 2], Rational::from(-2)),
        ]);

        let indexed_product = ckyz_series_mul_domain(&lhs, &rhs, &domain).unwrap();
        let reference_product = ckyz_series_mul(&lhs, &rhs, rank, max_total_degree).unwrap();
        assert_eq!(indexed_product, reference_product);

        let exp_input = BTreeMap::from([
            (vec![1, 0], Rational::from(2)),
            (vec![0, 1], Rational::from(-1)),
            (vec![1, 1], Rational::from(3)),
        ]);
        let indexed_exp = ckyz_series_exp_domain(&exp_input, &domain).unwrap();
        let reference_exp = ckyz_series_exp(&exp_input, rank, max_total_degree).unwrap();
        assert_eq!(indexed_exp, reference_exp);

        let composed_series = BTreeMap::from([
            (vec![1, 0], Rational::from(2)),
            (vec![0, 2], Rational::from(3)),
            (vec![1, 1], Rational::from(-1)),
        ]);
        let arguments = vec![
            BTreeMap::from([
                (vec![1, 0], Rational::from(1)),
                (vec![0, 1], Rational::from(2)),
            ]),
            BTreeMap::from([
                (vec![0, 1], Rational::from(1)),
                (vec![2, 0], Rational::from(-1)),
            ]),
        ];
        let indexed_compose =
            ckyz_series_compose_domain(&composed_series, &arguments, &domain).unwrap();
        let reference_compose =
            ckyz_series_compose(&composed_series, &arguments, max_total_degree).unwrap();
        assert_eq!(indexed_compose, reference_compose);
    }

    #[test]
    fn ckyz_cached_expalpha_q_series_matches_direct_exponentials() {
        let relations = ckyz_polygon5_relations();
        let local_intersection_terms = ckyz_polygon5_intersection_terms();
        let target_degrees = [vec![4, 3, 2], vec![8, 6, 4]];
        let extraction_degrees = ckyz_cover_closed_target_degrees(&target_degrees).unwrap();
        let domain = ckyz_predicted_support_domain_for_degrees(
            &relations,
            &local_intersection_terms,
            &extraction_degrees,
        )
        .unwrap();
        let alpha = compute_ckyz_log_period_corrections_domain(&relations, &domain).unwrap();
        let expalpha_powers =
            ckyz_expalpha_power_caches_domain(&alpha, &extraction_degrees, &domain).unwrap();

        for degree in extraction_degrees {
            let direct = ckyz_q_degree_series_in_z_domain(&degree, &alpha, &domain).unwrap();
            let cached = ckyz_q_degree_series_from_expalpha_powers_in_z_domain(
                &degree,
                &expalpha_powers,
                &domain,
            )
            .unwrap();
            assert_eq!(cached, direct);
        }
    }

    #[test]
    fn ckyz_previous_qn_indexed_series_matches_direct_exponentials() {
        let relations = ckyz_polygon5_relations();
        let local_intersection_terms = ckyz_polygon5_intersection_terms();
        let target_degrees = [vec![4, 3, 2], vec![8, 6, 4]];
        let extraction_degrees = ckyz_cover_closed_target_degrees(&target_degrees).unwrap();
        let domain = ckyz_predicted_support_domain_for_degrees(
            &relations,
            &local_intersection_terms,
            &extraction_degrees,
        )
        .unwrap();
        let alpha = compute_ckyz_log_period_corrections_domain(&relations, &domain).unwrap();
        let mut history_degrees =
            ckyz_z_residual_dependency_degrees(&alpha, &extraction_degrees, &domain).unwrap();
        let grading = ckyz_grading_vector_from_cover_weights(&[1, 1, 1]).unwrap();
        ckyz_sort_degrees_for_extraction_with_grading(&mut history_degrees, &grading).unwrap();
        history_degrees.dedup();
        let history_gradings = history_degrees
            .iter()
            .map(|degree| ckyz_grading_degree(degree, &grading).unwrap())
            .collect::<Vec<_>>();
        let history_indices = history_degrees
            .iter()
            .map(|degree| domain.index_of(degree).unwrap())
            .collect::<Vec<_>>();
        let indexed_alpha = ckyz_indexed_alpha_series(&alpha, &domain).unwrap();

        let previous_level_count = ckyz_cygv_previous_qn_level_count(domain.rank);
        let mut previous_qn = VecDeque::from(vec![
            HashMap::<usize, CkyzIndexedSeries>::new();
            previous_level_count
        ]);
        let mut previous_qn_indices =
            VecDeque::from(vec![Vec::<usize>::new(); previous_level_count]);
        let mut q_delta_cache = HashMap::<Vec<usize>, CkyzIndexedSeries>::new();
        let mut reused_previous_count = 0usize;
        let mut selected_deltas = BTreeSet::new();

        let mut batch_start = 0usize;
        while batch_start < history_indices.len() {
            let batch_grading = history_gradings[batch_start];
            let batch_end = history_gradings.partition_point(|&grading| grading <= batch_grading);
            let mut computed_qn = HashMap::new();
            let mut computed_indices = Vec::new();
            for (&degree_index, degree) in history_indices[batch_start..batch_end]
                .iter()
                .zip(history_degrees[batch_start..batch_end].iter())
            {
                let built = ckyz_indexed_q_degree_series_with_previous_cache_in_z_domain(
                    degree_index,
                    &previous_qn,
                    &previous_qn_indices,
                    &mut q_delta_cache,
                    &indexed_alpha,
                    &domain,
                )
                .unwrap();
                let direct = ckyz_q_degree_series_in_z_domain(degree, &alpha, &domain).unwrap();
                assert_eq!(built.series.to_btree(&domain), direct, "degree={degree:?}");
                let indexed_li2 = built.series.li2(&domain).unwrap().to_btree(&domain);
                let direct_li2 = ckyz_series_li2_domain(&direct, &domain).unwrap();
                assert_eq!(indexed_li2, direct_li2, "Li2 degree={degree:?}");
                if built.reused_previous {
                    reused_previous_count += 1;
                }
                if built.delta_degree.iter().any(|&entry| entry != 0) {
                    selected_deltas.insert(built.delta_degree);
                }
                computed_qn.insert(degree_index, built.series);
                computed_indices.push(degree_index);
            }
            previous_qn.pop_front();
            previous_qn_indices.pop_front();
            previous_qn.push_back(computed_qn);
            previous_qn_indices.push_back(computed_indices);
            batch_start = batch_end;
        }

        assert!(reused_previous_count > 0);
        assert!(selected_deltas.len() < history_degrees.len());
    }

    #[test]
    fn ckyz_coefficient_li2_matches_full_series_on_polygon5() {
        let relations = ckyz_polygon5_relations();
        let local_intersection_terms = ckyz_polygon5_intersection_terms();
        let target_degrees = [vec![4, 3, 2], vec![8, 6, 4]];
        let extraction_degrees = ckyz_cover_closed_target_degrees(&target_degrees).unwrap();
        let domain = ckyz_predicted_support_domain_for_degrees(
            &relations,
            &local_intersection_terms,
            &extraction_degrees,
        )
        .unwrap();
        let alpha = compute_ckyz_log_period_corrections_domain(&relations, &domain).unwrap();
        let alpha_terms = ckyz_scaled_alpha_terms(&alpha, &domain).unwrap();
        let mut coefficient_cache = CkyzExpCoefficientCache::default();

        for degree in &extraction_degrees {
            let q_degree = ckyz_q_degree_series_in_z_domain(degree, &alpha, &domain).unwrap();
            let full_li2 = ckyz_series_li2_domain(&q_degree, &domain).unwrap();
            for target in &extraction_degrees {
                let actual = ckyz_q_degree_li2_coefficient_in_z_domain(
                    degree,
                    target,
                    &alpha_terms,
                    &domain,
                    None,
                    &mut coefficient_cache,
                )
                .unwrap();
                let expected = full_li2
                    .get(target)
                    .cloned()
                    .unwrap_or_else(|| Rational::from(0));
                assert_eq!(actual, expected, "degree={degree:?} target={target:?}");
            }
        }
    }

    #[test]
    fn ckyz_z_series_inversion_matches_predicted_support_domain() {
        let cases = [
            (
                vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]],
                vec![CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 1,
                    coefficient: 1,
                }],
                vec![2, 2],
                vec![vec![3, 2], vec![6, 4], vec![9, 6]],
            ),
            (
                ckyz_polygon5_relations(),
                ckyz_polygon5_intersection_terms(),
                vec![1, 1, 1],
                vec![vec![4, 3, 2], vec![8, 6, 4]],
            ),
        ];

        for (relations, local_intersection_terms, cover_weights, target_degrees) in cases {
            let extraction_degrees = ckyz_cover_closed_target_degrees(&target_degrees).unwrap();
            let domain = ckyz_predicted_support_domain_for_degrees(
                &relations,
                &local_intersection_terms,
                &extraction_degrees,
            )
            .unwrap();
            let potential = compute_ckyz_local_instanton_potential_corrections_domain(
                &relations,
                &local_intersection_terms,
                &domain,
            )
            .unwrap();
            let mut expected = extract_ckyz_local_gv_invariants_from_potential_for_degrees(
                &potential,
                &cover_weights,
                &extraction_degrees,
            )
            .unwrap();
            expected.retain(|degree, _| target_degrees.iter().any(|target| target == degree));

            let (alpha, potential_z) = compute_ckyz_local_instanton_potential_z_domain(
                &relations,
                &local_intersection_terms,
                &domain,
            )
            .unwrap();
            let path_history_degrees =
                ckyz_z_residual_dependency_degrees(&alpha, &extraction_degrees, &domain).unwrap();
            let mut actual = extract_ckyz_local_gv_invariants_from_z_potential_for_degrees(
                &potential_z,
                &alpha,
                &cover_weights,
                &path_history_degrees,
                &domain,
            )
            .unwrap();
            actual.retain(|degree, _| target_degrees.iter().any(|target| target == degree));

            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn ckyz_support_predicates_match_rational_period_supports() {
        let cases = [
            (
                vec![vec![-2, 1, 0, 1, 0], vec![-2, 0, 1, 0, 1]],
                vec![CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 1,
                    coefficient: 1,
                }],
                vec![vec![3, 2], vec![6, 4], vec![9, 6]],
            ),
            (
                ckyz_polygon5_relations(),
                ckyz_polygon5_intersection_terms(),
                vec![vec![4, 3, 2], vec![8, 6, 4]],
            ),
        ];

        for (relations, local_intersection_terms, target_degrees) in cases {
            let extraction_degrees = ckyz_cover_closed_target_degrees(&target_degrees).unwrap();
            let domain =
                CkyzMonomialDomain::target_downset(&extraction_degrees, relations.len()).unwrap();
            let alpha = compute_ckyz_log_period_corrections_domain(&relations, &domain).unwrap();
            let alpha_supports =
                ckyz_log_period_support_indices_domain(&relations, &domain).unwrap();
            assert_eq!(alpha_supports.len(), alpha.len());
            for (support, series) in alpha_supports.iter().zip(alpha.iter()) {
                assert_eq!(*support, ckyz_series_support_indices(series, &domain));
            }

            for term in local_intersection_terms {
                let beta = ckyz_second_log_period_series_for_pair_domain(
                    &relations,
                    term.first,
                    term.second,
                    &domain,
                )
                .unwrap();
                let beta_support = ckyz_second_log_period_support_indices_for_pair_domain(
                    &relations,
                    term.first,
                    term.second,
                    &domain,
                )
                .unwrap();
                assert_eq!(beta_support, ckyz_series_support_indices(&beta, &domain));
            }
        }
    }

    #[test]
    fn ckyz_support_exponential_uses_same_semigroup_as_power_union() {
        let domain = CkyzMonomialDomain::target_downset(&[vec![4, 3, 2]], 3).unwrap();
        let support = [vec![1, 0, 0], vec![0, 1, 1], vec![2, 1, 0]]
            .into_iter()
            .map(|degree| {
                domain
                    .index_of(&degree)
                    .expect("degree should be in domain")
            })
            .collect::<BTreeSet<_>>();

        let by_closure = ckyz_support_exp_domain(&support, &domain).unwrap();
        let by_powers = ckyz_support_exp_domain_by_powers(&support, &domain).unwrap();

        assert_eq!(by_closure, by_powers);
    }

    #[test]
    fn ckyz_large_domain_uses_dense_degree_index_without_addition_table() {
        let domain = CkyzMonomialDomain::componentwise_box(&[49, 49]).unwrap();
        assert!(domain.addition_indices.is_none());
        assert!(domain.dense_degree_indices.is_some());

        let lhs = BTreeMap::from([
            (vec![20, 20], Rational::from(2)),
            (vec![30, 30], Rational::from(5)),
        ]);
        let rhs = BTreeMap::from([(vec![25, 25], Rational::from(3))]);

        let product = ckyz_series_mul_domain(&lhs, &rhs, &domain).unwrap();
        assert_eq!(product, BTreeMap::from([(vec![45, 45], Rational::from(6))]));
    }

    #[test]
    fn ckyz_local_domain_profile_reports_polygon5_support_reduction() {
        let relations = ckyz_polygon5_relations();
        let local_intersection_terms = ckyz_polygon5_intersection_terms();
        let target_degrees = [vec![4, 3, 2], vec![8, 6, 4]];

        let profile = ckyz_local_domain_profile_for_degrees(
            &relations,
            &local_intersection_terms,
            &target_degrees,
        )
        .unwrap();

        assert_eq!(profile.rank, 3);
        assert_eq!(profile.extraction_degree_count, 2);
        assert_eq!(profile.max_target_total_degree, 18);
        assert_eq!(profile.target_downset_degree_count, 315);
        assert_eq!(profile.predicted_support_degree_count, 265);
        assert_eq!(profile.predicted_support_max_total_degree, 18);
        assert!(profile.target_downset_has_addition_table);
        assert!(profile.predicted_support_has_addition_table);
        assert_eq!(profile.causal_semigroup_degree_count, None);
    }

    #[test]
    fn ckyz_local_surface_causal_domain_spec_uses_source_weights() {
        let identification = CkyzLocalSurfaceIdentification {
            kind: CkyzLocalSurfaceKind::HirzebruchF1,
            point_permutation: vec![0, 1, 2, 3, 4],
            row_transform: vec![vec![1, 0], vec![0, 1]],
            source_relations: vec![vec![-2, 1, 0, 1, 0], vec![-1, 0, 1, -1, 1]],
            source_target_direction: vec![5, 4],
            c1_coefficients: vec![3, 2],
            local_intersection_terms: vec![
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 1,
                    coefficient: 1,
                },
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 0,
                    coefficient: 1,
                },
            ],
        };

        let spec = ckyz_local_surface_causal_domain_spec(&identification).unwrap();
        assert_eq!(spec.generators, vec![vec![1, 0], vec![0, 1]]);
        assert_eq!(spec.grading_vector, vec![2, 1]);
        assert_eq!(
            ckyz_local_surface_target_degrees(&identification, 3).unwrap(),
            vec![vec![5, 4], vec![10, 8], vec![15, 12]]
        );
    }

    #[test]
    fn ckyz_local_surface_domain_profile_reports_causal_f1_semigroup() {
        let identification = CkyzLocalSurfaceIdentification {
            kind: CkyzLocalSurfaceKind::HirzebruchF1,
            point_permutation: vec![0, 1, 2, 3, 4],
            row_transform: vec![vec![1, 0], vec![0, 1]],
            source_relations: vec![vec![-2, 1, 0, 1, 0], vec![-1, 0, 1, -1, 1]],
            source_target_direction: vec![2, 1],
            c1_coefficients: vec![3, 2],
            local_intersection_terms: vec![
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 1,
                    coefficient: 1,
                },
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 0,
                    coefficient: 1,
                },
            ],
        };

        let surface_profile =
            ckyz_local_surface_domain_profile_for_multiples(&identification, 3).unwrap();
        let target_degrees = ckyz_local_surface_target_degrees(&identification, 3).unwrap();
        let direct_profile = ckyz_local_domain_profile_for_degrees_with_causal_domain(
            &identification.source_relations,
            &identification.local_intersection_terms,
            &target_degrees,
            &[vec![1, 0], vec![0, 1]],
            &[2, 1],
        )
        .unwrap();

        assert_eq!(surface_profile, direct_profile);
        assert_eq!(surface_profile.rank, 2);
        assert_eq!(surface_profile.extraction_degree_count, 3);
        assert_eq!(surface_profile.target_downset_degree_count, 28);
        assert_eq!(surface_profile.causal_semigroup_degree_count, Some(72));
        assert_eq!(surface_profile.causal_semigroup_max_total_degree, Some(15));
        assert_eq!(
            surface_profile.causal_semigroup_has_addition_table,
            Some(true)
        );
    }

    #[test]
    fn ckyz_local_surface_causal_helper_matches_targeted_f1_ray() {
        let identification = CkyzLocalSurfaceIdentification {
            kind: CkyzLocalSurfaceKind::HirzebruchF1,
            point_permutation: vec![0, 1, 2, 3, 4],
            row_transform: vec![vec![1, 0], vec![0, 1]],
            source_relations: vec![vec![-2, 1, 0, 1, 0], vec![-1, 0, 1, -1, 1]],
            source_target_direction: vec![2, 1],
            c1_coefficients: vec![3, 2],
            local_intersection_terms: vec![
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 1,
                    coefficient: 1,
                },
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 0,
                    coefficient: 1,
                },
            ],
        };
        let target_degrees = ckyz_local_surface_target_degrees(&identification, 3).unwrap();
        let targeted = compute_ckyz_local_gv_invariants_for_degrees(
            &identification.source_relations,
            &identification.local_intersection_terms,
            ckyz_local_surface_cover_weight_coefficients(&identification.kind),
            &target_degrees,
        )
        .unwrap();
        let causal = compute_ckyz_local_surface_gv_invariants_for_multiples_with_causal_domain(
            &identification,
            3,
        )
        .unwrap();

        assert_eq!(causal, targeted);
    }

    #[test]
    fn ckyz_local_gv_extraction_rejects_nonintegral_residuals() {
        let potential = BTreeMap::from([(vec![1], Rational::from_signeds(1, 2))]);
        let err = extract_ckyz_local_gv_invariants_from_potential(&potential, &[1], 1).unwrap_err();

        assert!(err.to_string().contains("not integral"));
    }

    #[test]
    fn ckyz_local_gv_extraction_matches_f1_table_start() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-1, 0, 1, -1, 1]];
        let gvs = compute_ckyz_local_gv_invariants(
            &relations,
            &[
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 1,
                    coefficient: 1,
                },
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 0,
                    coefficient: 1,
                },
            ],
            &[2, 1],
            4,
        )
        .unwrap();

        let expected_nonzero = BTreeMap::from([
            (vec![1, 0], Integer::from(-2)),
            (vec![0, 1], Integer::from(1)),
            (vec![1, 1], Integer::from(3)),
            (vec![2, 1], Integer::from(5)),
            (vec![3, 1], Integer::from(7)),
            (vec![2, 2], Integer::from(-6)),
        ]);

        for (degree, expected) in expected_nonzero {
            assert_eq!(gvs.get(&degree), Some(&expected));
        }
        for zero_degree in [vec![2, 0], vec![0, 2], vec![3, 0], vec![0, 3]] {
            assert!(
                !gvs.contains_key(&zero_degree),
                "zero F1 invariant should be omitted at {zero_degree:?}"
            );
        }
    }

    #[test]
    fn ckyz_local_gv_extraction_rejects_wrong_f1_cover_weights() {
        let relations = vec![vec![-2, 1, 0, 1, 0], vec![-1, 0, 1, -1, 1]];
        let err = compute_ckyz_local_gv_invariants(
            &relations,
            &[
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 1,
                    coefficient: 1,
                },
                CkyzLocalIntersectionTerm {
                    first: 0,
                    second: 0,
                    coefficient: 1,
                },
            ],
            &[3, 2],
            2,
        )
        .unwrap_err();

        assert!(err.to_string().contains("not integral"));
    }

    fn ckyz_polygon5_relations() -> Vec<Vec<i64>> {
        vec![
            vec![-1, 1, -1, 1, 0, 0],
            vec![-1, -1, 1, 0, 0, 1],
            vec![-1, 0, 1, -1, 1, 0],
        ]
    }

    fn ckyz_polygon5_intersection_terms() -> Vec<CkyzLocalIntersectionTerm> {
        vec![
            CkyzLocalIntersectionTerm {
                first: 0,
                second: 0,
                coefficient: 1,
            },
            CkyzLocalIntersectionTerm {
                first: 1,
                second: 0,
                coefficient: 1,
            },
            CkyzLocalIntersectionTerm {
                first: 0,
                second: 2,
                coefficient: 1,
            },
            CkyzLocalIntersectionTerm {
                first: 1,
                second: 2,
                coefficient: 1,
            },
        ]
    }

    #[test]
    fn ckyz_targeted_local_gv_extraction_matches_polygon5_table_start() {
        let target_degrees = [
            vec![0, 0, 1],
            vec![0, 1, 0],
            vec![0, 1, 1],
            vec![1, 0, 0],
            vec![1, 0, 1],
            vec![1, 1, 0],
            vec![1, 1, 1],
            vec![2, 1, 1],
            vec![2, 1, 2],
            vec![2, 2, 1],
            vec![2, 2, 2],
        ];

        let gvs = compute_ckyz_local_gv_invariants_for_degrees(
            &ckyz_polygon5_relations(),
            &ckyz_polygon5_intersection_terms(),
            &[1, 1, 1],
            &target_degrees,
        )
        .unwrap();

        let expected_nonzero = BTreeMap::from([
            (vec![0, 0, 1], Integer::from(1)),
            (vec![0, 1, 0], Integer::from(1)),
            (vec![1, 0, 0], Integer::from(1)),
            (vec![1, 0, 1], Integer::from(-2)),
            (vec![1, 1, 0], Integer::from(-2)),
            (vec![1, 1, 1], Integer::from(3)),
            (vec![2, 1, 1], Integer::from(-4)),
            (vec![2, 1, 2], Integer::from(5)),
            (vec![2, 2, 1], Integer::from(5)),
            (vec![2, 2, 2], Integer::from(-6)),
        ]);

        for (degree, expected) in expected_nonzero {
            assert_eq!(gvs.get(&degree), Some(&expected));
        }
        assert!(
            !gvs.contains_key(&vec![0, 1, 1]),
            "zero polygon-5 invariant should be omitted"
        );
    }

    #[test]
    fn ckyz_local_gv_extraction_rejects_printed_polygon5_c1_weights() {
        let err = compute_ckyz_local_gv_invariants_for_degrees(
            &ckyz_polygon5_relations(),
            &ckyz_polygon5_intersection_terms(),
            &[3, 2, 2],
            &[vec![0, 0, 1]],
        )
        .unwrap_err();

        assert!(err.to_string().contains("not integral"));
    }

    #[test]
    fn local_p2_circuit_dispatches_to_local_mirror_series() {
        let points = vec![
            Point::new(vec![0, 1, -3, 6]),
            Point::new(vec![-2, -1, -4, 5]),
            Point::new(vec![-1, 0, -3, 5]),
            Point::new(vec![-1, 0, -2, 4]),
        ];

        let diagnostic = diagnose_affine_toric_circuit(&[1, 1, -3, 1], &points)
            .unwrap()
            .expect("local P2 row is an affine circuit");
        let kind = diagnostic.kind.expect("local P2 kind should be recognized");
        let gvs = compute_local_toric_circuit_gv_series(&kind, 4)
            .unwrap()
            .expect("local P2 circuit should have a GV series");

        assert_eq!(
            gvs,
            vec![
                Integer::from(3),
                Integer::from(-6),
                Integer::from(27),
                Integer::from(-192),
            ]
        );
    }

    #[test]
    fn rank_two_local_support_signature_ignores_relation_orientation() {
        let points = vec![
            Point::new(vec![0, 1, -3, 6]),
            Point::new(vec![-2, -1, -4, 5]),
            Point::new(vec![-1, 0, -3, 5]),
            Point::new(vec![-1, 0, -2, 4]),
        ];

        let positive = diagnose_affine_toric_circuit(&[1, 1, -3, 1], &points)
            .unwrap()
            .expect("local P2 row is an affine circuit");
        let negative = diagnose_affine_toric_circuit(&[-1, -1, 3, -1], &points)
            .unwrap()
            .expect("opposite local P2 row is an affine circuit");

        assert_eq!(
            rank_two_local_support_signature(&positive),
            rank_two_local_support_signature(&negative)
        );
    }

    #[test]
    fn rank_two_local_charge_model_recovers_affine_kernel() {
        let points = vec![
            Point::new(vec![0, 1, -3, 6]),
            Point::new(vec![-2, -1, -4, 5]),
            Point::new(vec![-1, 0, -3, 5]),
            Point::new(vec![-1, 0, -2, 4]),
        ];

        let diagnostic = diagnose_affine_toric_circuit(&[1, 1, -3, 1], &points)
            .unwrap()
            .expect("local P2 row is an affine circuit");
        let signature = rank_two_local_support_signature(&diagnostic)
            .expect("local P2 row should have a signature");
        let model =
            rank_two_local_charge_model(&signature).expect("local P2 should have a charge model");

        assert_eq!(model.points.len(), 4);
        assert_eq!(model.charge_basis.len(), 1);
        assert_eq!(model.target_relation_in_charge_basis, vec![-1]);
        assert!(
            curve_in_rational_row_span(&model.target_relation, &model.charge_basis)
                .expect("charge span check should be exact")
        );
        let reconstructed_target: Vec<i64> = model.charge_basis[0]
            .iter()
            .map(|&value| model.target_relation_in_charge_basis[0] * value)
            .collect();
        assert_eq!(reconstructed_target, model.target_relation);
        for charge in &model.charge_basis {
            let coefficient_sum: i64 = charge.iter().sum();
            assert_eq!(coefficient_sum, 0);
            let coordinate_sum = charge.iter().zip(model.points.iter()).fold(
                [0i64; 2],
                |mut sum, (&coefficient, point)| {
                    sum[0] += coefficient * point.coordinates[0];
                    sum[1] += coefficient * point.coordinates[1];
                    sum
                },
            );
            assert_eq!(coordinate_sum, [0, 0]);
        }
    }

    #[test]
    fn affine_toric_circuit_records_full_rank_local_coordinates() {
        let points = vec![
            Point::new(vec![0, 0, 0, 0]),
            Point::new(vec![1, 2, 0, 2]),
            Point::new(vec![2, 2, 1, 2]),
            Point::new(vec![2, 3, 1, 3]),
            Point::new(vec![2, 3, 2, 3]),
        ];
        let relation = vec![-1, 2, 1, -3, 1];

        let diagnostic = diagnose_affine_toric_circuit(&relation, &points)
            .unwrap()
            .expect("rank-three support should be an affine circuit");

        assert_eq!(diagnostic.affine_rank, 3);
        assert_eq!(diagnostic.local_coordinates.len(), points.len());
        assert!(diagnostic.local_coordinates_2d.is_none());
        for point in &diagnostic.local_coordinates {
            assert_eq!(point.coordinates.len(), 3);
        }
        for coordinate_idx in 0..3 {
            let weighted_sum: i64 = diagnostic
                .local_coordinates
                .iter()
                .zip(relation.iter())
                .map(|(point, &coefficient)| coefficient * point.coordinates[coordinate_idx])
                .sum();
            assert_eq!(weighted_sum, 0);
        }
    }

    #[test]
    fn affine_toric_circuit_detects_local_p2_triangle_opposite_orientation() {
        let points = vec![
            Point::new(vec![0, 1, -3, 6]),
            Point::new(vec![-2, -1, -4, 5]),
            Point::new(vec![-1, 0, -3, 5]),
            Point::new(vec![-1, 0, -2, 4]),
        ];

        let diagnostic = diagnose_affine_toric_circuit(&[-1, -1, 3, -1], &points)
            .unwrap()
            .expect("opposite orientation is still an affine circuit");

        assert_eq!(
            diagnostic.coefficient_counts,
            BTreeMap::from([(-1, 3), (3, 1)])
        );
        assert_eq!(diagnostic.affine_rank, 2);
        assert_eq!(diagnostic.local_charge_basis, vec![vec![1, 1, -3, 1]]);
        assert!(
            diagnostic.local_coordinates_2d.is_some(),
            "rank-two local P2 support should expose local coordinates"
        );
        assert_eq!(
            diagnostic.kind,
            Some(LocalToricCircuitKind::LocalP2Triangle {
                interior_point: 2,
                vertex_points: vec![0, 1, 3],
                interior_coefficient: 3,
                vertex_coefficient: -1,
                local_coordinates: vec![
                    LocalToricCoordinate2D {
                        point_index: 0,
                        coordinates: [1, 0],
                    },
                    LocalToricCoordinate2D {
                        point_index: 1,
                        coordinates: [0, 1],
                    },
                    LocalToricCoordinate2D {
                        point_index: 2,
                        coordinates: [0, 0],
                    },
                    LocalToricCoordinate2D {
                        point_index: 3,
                        coordinates: [-1, -1],
                    },
                ],
            })
        );
    }

    #[test]
    fn affine_toric_circuit_rejects_non_affine_row() {
        let points = vec![
            Point::new(vec![0, 1, -3, 6]),
            Point::new(vec![-2, -1, -4, 5]),
            Point::new(vec![-1, 0, -3, 5]),
            Point::new(vec![-1, 0, -2, 4]),
        ];

        let diagnostic = diagnose_affine_toric_circuit(&[1, 1, -2, 1], &points).unwrap();

        assert_eq!(diagnostic, None);
    }

    #[test]
    fn affine_toric_circuit_requires_rank_two_for_local_p2() {
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![2, 0]),
            Point::new(vec![3, 0]),
            Point::new(vec![7, 0]),
        ];

        let diagnostic = diagnose_affine_toric_circuit(&[1, 1, -3, 1], &points)
            .unwrap()
            .expect("collinear barycenter row is still an affine circuit");

        assert_eq!(diagnostic.affine_rank, 1);
        assert_eq!(diagnostic.kind, None);
    }

    #[test]
    fn origin_circuit_diagnostic_marks_resolved_conifold_pattern() {
        let diagnostic = origin_circuit_diagnostic_from_class_and_witnesses(
            vec![-1, -1, -1, 1, 1, 1],
            0,
            Vec::new(),
        );

        assert_eq!(diagnostic.origin_coefficient, -1);
        assert_eq!(
            diagnostic.negative_coefficient_counts,
            BTreeMap::from([(-1, 2)])
        );
        assert_eq!(
            diagnostic.positive_coefficient_counts,
            BTreeMap::from([(1, 3)])
        );
        assert!(diagnostic.is_resolved_conifold_pattern);
        assert!(diagnostic.witnesses.is_empty());
    }

    #[test]
    fn origin_circuit_diagnostic_marks_standard_resolved_conifold_charge() {
        let diagnostic =
            origin_circuit_diagnostic_from_class_and_witnesses(vec![-1, -1, 1, 1], 0, Vec::new());

        assert_eq!(diagnostic.origin_coefficient, -1);
        assert_eq!(
            diagnostic.negative_coefficient_counts,
            BTreeMap::from([(-1, 1)])
        );
        assert_eq!(
            diagnostic.positive_coefficient_counts,
            BTreeMap::from([(1, 2)])
        );
        assert!(diagnostic.is_resolved_conifold_pattern);
        assert_eq!(
            resolved_conifold_origin_circuit_gv(&diagnostic.class, 0),
            Some(Integer::from(1))
        );
    }

    #[test]
    fn origin_circuit_diagnostic_keeps_unsupported_patterns_explicit() {
        let diagnostic = origin_circuit_diagnostic_from_class_and_witnesses(
            vec![-2, -1, -2, 1, 2, 2],
            0,
            Vec::new(),
        );

        assert_eq!(diagnostic.origin_coefficient, -2);
        assert_eq!(
            diagnostic.negative_coefficient_counts,
            BTreeMap::from([(-2, 1), (-1, 1)])
        );
        assert_eq!(
            diagnostic.positive_coefficient_counts,
            BTreeMap::from([(1, 1), (2, 2)])
        );
        assert!(!diagnostic.is_resolved_conifold_pattern);
    }

    #[test]
    fn origin_circuit_diagnostic_keeps_sorted_witnesses() {
        let first = OriginCircuitCurveWitness {
            class: vec![-1, -1, -1, 1, 1, 1],
            first_facet_exclusive_point: 5,
            second_facet_exclusive_point: 1,
            shared_two_simplex: vec![2, 3, 4],
            first_facet: vec![2, 3, 4, 5],
            second_facet: vec![1, 2, 3, 4],
            sparse_relation: vec![(0, -1), (1, -1), (5, 1)],
            relation_points: vec![OriginCircuitRelationPoint {
                point_index: 0,
                coefficient: -1,
                coordinates: vec![0, 0, 0, 0],
                face_dimension: Some(4),
            }],
        };
        let second = OriginCircuitCurveWitness {
            class: vec![-1, -1, -1, 1, 1, 1],
            first_facet_exclusive_point: 2,
            second_facet_exclusive_point: 4,
            shared_two_simplex: vec![1, 3, 5],
            first_facet: vec![1, 2, 3, 5],
            second_facet: vec![1, 3, 4, 5],
            sparse_relation: vec![(0, -1), (2, 1), (4, -1)],
            relation_points: vec![OriginCircuitRelationPoint {
                point_index: 2,
                coefficient: 1,
                coordinates: vec![1, 0, 0, 0],
                face_dimension: Some(0),
            }],
        };

        let diagnostic = origin_circuit_diagnostic_from_class_and_witnesses(
            vec![-1, -1, -1, 1, 1, 1],
            0,
            vec![first.clone(), second.clone()],
        );

        assert_eq!(diagnostic.witnesses, vec![second, first]);
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
    fn ambient_curve_projection_keeps_basis_entries() {
        let ambient = vec![7, -3, 0, 11];
        let basis = vec![3, 1];

        let projected = project_ambient_curve_to_basis(&ambient, &basis).unwrap();

        assert_eq!(projected, vec![11, -3]);
    }

    #[test]
    fn ambient_curve_projection_rejects_out_of_bounds_basis() {
        let err = project_ambient_curve_to_basis(&[1, 2], &[0, 3]).unwrap_err();

        assert!(err.to_string().contains("basis index 3 is out of bounds"));
    }

    #[test]
    fn ambient_curve_projection_matrix_matches_cytools_dot_basis_transpose() {
        let ambient = vec![2, -3, 5, 7];
        let basis_matrix = vec![vec![1, 0, 0, 0], vec![0, 2, -1, 1]];

        let projected = project_ambient_curve_to_basis_matrix(&ambient, &basis_matrix).unwrap();

        assert_eq!(projected, vec![2, -4]);
    }

    #[test]
    fn ambient_curve_projection_matrix_rejects_bad_width() {
        let err = project_ambient_curve_to_basis_matrix(&[1, 2], &[vec![1, 0, 0]])
            .expect_err("matrix basis row width must match ambient dimension");

        assert!(
            err.to_string()
                .contains("does not match ambient curve dimension")
        );
    }

    #[test]
    fn supporting_mori_face_normal_certifies_exact_face() {
        let certificate = check_supporting_mori_face_normal(
            &[0, 1],
            &[vec![1, 0]],
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
        )
        .unwrap()
        .unwrap();

        assert_eq!(certificate.normal, vec![0, 1]);
        assert_eq!(certificate.zero_generator_count, 1);
        assert_eq!(certificate.positive_generator_count, 2);
    }

    #[test]
    fn supporting_mori_face_normal_rejects_nonvanishing_face_generator() {
        let certificate = check_supporting_mori_face_normal(
            &[0, 1],
            &[vec![1, 1]],
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
        )
        .unwrap();

        assert!(certificate.is_none());
    }

    #[test]
    fn supporting_mori_face_normal_rejects_negative_mori_pairing() {
        let certificate = check_supporting_mori_face_normal(
            &[0, 1],
            &[vec![1, 0]],
            &[vec![1, 0], vec![0, -1], vec![1, 1]],
        )
        .unwrap();

        assert!(certificate.is_none());
    }

    #[test]
    fn supporting_mori_face_normal_rejects_zero_normal() {
        let certificate =
            check_supporting_mori_face_normal(&[0, 0], &[vec![1, 0]], &[vec![1, 0]]).unwrap();

        assert!(certificate.is_none());
    }

    #[test]
    fn supporting_mori_face_normal_checks_dimensions() {
        let err = check_supporting_mori_face_normal(&[0, 1], &[vec![1, 0, 0]], &[vec![1, 0]])
            .unwrap_err();

        assert!(err.to_string().contains("face generator dimension 3"));
    }

    #[test]
    fn supporting_mori_face_normal_reports_dot_overflow() {
        let err = check_supporting_mori_face_normal(
            &[i64::MAX, i64::MAX, i64::MAX],
            &[],
            &[vec![i64::MAX, i64::MAX, i64::MAX]],
        )
        .unwrap_err();

        assert!(err.to_string().contains("dot product overflowed"));
    }

    #[test]
    fn supporting_mori_face_exact_kernel_certifies_codimension_one_face() {
        let certificate = certify_supporting_mori_face_by_exact_kernel(
            &[vec![1, 0]],
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
        )
        .unwrap()
        .unwrap();

        assert_eq!(certificate.normal, vec![0, 1]);
        assert_eq!(certificate.zero_generator_count, 1);
        assert_eq!(certificate.positive_generator_count, 2);
    }

    #[test]
    fn supporting_mori_face_exact_kernel_orients_the_normal() {
        let certificate = certify_supporting_mori_face_by_exact_kernel(
            &[vec![1, 0]],
            &[vec![1, 0], vec![0, -1], vec![1, -1]],
        )
        .unwrap()
        .unwrap();

        assert_eq!(certificate.normal, vec![0, -1]);
        assert_eq!(certificate.zero_generator_count, 1);
        assert_eq!(certificate.positive_generator_count, 2);
    }

    #[test]
    fn supporting_mori_face_exact_kernel_declines_higher_codimension_faces() {
        let certificate = certify_supporting_mori_face_by_exact_kernel(
            &[vec![1, 0, 0]],
            &[vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]],
        )
        .unwrap();

        assert!(certificate.is_none());
    }

    #[test]
    fn supporting_mori_face_lp_search_certifies_higher_codimension_face() {
        let certificate = certify_supporting_mori_face_by_lp_search(
            &[vec![1, 0, 0]],
            &[vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1], vec![1, 1, 1]],
            &SupportingMoriFaceLpSearchOptions::default(),
        )
        .unwrap()
        .unwrap();

        assert_eq!(
            exact_i64_dot_checked(&certificate.normal, &[1, 0, 0]).unwrap(),
            0
        );
        assert!(certificate.zero_generator_count >= 1);
        assert!(certificate.positive_generator_count >= 1);
        assert_eq!(
            certificate.zero_generator_count + certificate.positive_generator_count,
            4
        );
    }

    #[test]
    fn supporting_mori_face_lp_diagnostic_reports_certified_phase() {
        let diagnostic = diagnose_supporting_mori_face_by_lp_search(
            &[vec![1, 0, 0]],
            &[vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1], vec![1, 1, 1]],
            &SupportingMoriFaceLpSearchOptions::default(),
        )
        .unwrap();

        assert_eq!(diagnostic.face_rank, 1);
        assert_eq!(diagnostic.dim, 3);
        assert_eq!(diagnostic.exact_kernel_status, "no_certificate");
        assert!(diagnostic.certificate.is_some());
        assert!(
            diagnostic.status == "certified_full_lp"
                || diagnostic.status == "certified_aggregate_lp"
                || diagnostic.status == "certified_anchor_lp"
        );
        if diagnostic.status == "certified_full_lp" {
            assert_eq!(diagnostic.full_status, "certified");
            assert_eq!(diagnostic.anchor_attempt_count, 0);
        }
        if diagnostic.status == "certified_aggregate_lp" {
            assert_eq!(diagnostic.aggregate_status, "certified");
            assert_eq!(diagnostic.anchor_attempt_count, 0);
        }
    }

    #[test]
    fn supporting_mori_face_lp_search_rejects_full_dimensional_support() {
        let certificate = certify_supporting_mori_face_by_lp_search(
            &[vec![1, 0], vec![0, 1]],
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
            &SupportingMoriFaceLpSearchOptions::default(),
        )
        .unwrap();

        assert!(certificate.is_none());
    }

    #[test]
    fn supporting_mori_face_lp_diagnostic_reports_full_dimensional_failure() {
        let diagnostic = diagnose_supporting_mori_face_by_lp_search(
            &[vec![1, 0], vec![0, 1]],
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
            &SupportingMoriFaceLpSearchOptions::default(),
        )
        .unwrap();

        assert_eq!(diagnostic.face_rank, 2);
        assert_eq!(diagnostic.dim, 2);
        assert_eq!(diagnostic.status, "lp_no_certificate");
        assert_eq!(diagnostic.exact_kernel_status, "no_certificate");
        assert_eq!(diagnostic.full_status, "lp_no_solution");
        assert_eq!(diagnostic.aggregate_status, "lp_no_solution");
        assert_eq!(diagnostic.anchor_attempt_count, 0);
        assert!(diagnostic.certificate.is_none());
    }

    #[test]
    fn supporting_mori_face_lp_search_declines_empty_support() {
        let certificate = certify_supporting_mori_face_by_lp_search(
            &[],
            &[vec![1, 0], vec![0, 1]],
            &SupportingMoriFaceLpSearchOptions::default(),
        )
        .unwrap();

        assert!(certificate.is_none());
    }

    #[test]
    fn supporting_mori_face_lp_search_checks_options() {
        let options = SupportingMoriFaceLpSearchOptions {
            scale_limit: 0,
            ..SupportingMoriFaceLpSearchOptions::default()
        };
        let err = certify_supporting_mori_face_by_lp_search(
            &[vec![1, 0, 0]],
            &[vec![1, 0, 0], vec![0, 1, 0]],
            &options,
        )
        .unwrap_err();

        assert!(err.to_string().contains("scale limit must be positive"));
    }

    #[test]
    fn supporting_mori_face_exact_kernel_checks_dimensions() {
        let err = certify_supporting_mori_face_by_exact_kernel(
            &[vec![1, 0]],
            &[vec![1, 0], vec![0, 1, 0]],
        )
        .unwrap_err();

        assert!(err.to_string().contains("Mori generator dimension 3"));
    }

    #[test]
    fn supporting_mori_face_from_normal_extracts_zero_generators() {
        let face = supporting_mori_face_from_normal(&[0, 1], &[vec![1, 0], vec![0, 1], vec![1, 1]])
            .unwrap()
            .unwrap();

        assert_eq!(face.generators, vec![vec![1, 0]]);
        assert_eq!(face.certificate.zero_generator_count, 1);
        assert_eq!(face.certificate.positive_generator_count, 2);
    }

    #[test]
    fn supporting_mori_face_from_normal_allows_zero_face() {
        let face = supporting_mori_face_from_normal(&[1, 1], &[vec![1, 0], vec![0, 1]])
            .unwrap()
            .unwrap();

        assert!(face.generators.is_empty());
        assert_eq!(face.certificate.zero_generator_count, 0);
        assert_eq!(face.certificate.positive_generator_count, 2);
    }

    #[test]
    fn supporting_mori_face_from_normal_rejects_non_supporting_normal() {
        let face = supporting_mori_face_from_normal(&[0, 1], &[vec![1, 0], vec![0, -1]]).unwrap();

        assert!(face.is_none());
    }

    #[test]
    fn supporting_mori_face_for_curve_requires_curve_in_face_span() {
        let mori_generators = vec![vec![1, 0], vec![0, 1], vec![1, 1]];

        let face = supporting_mori_face_for_curve_from_normal(&[0, 1], &[2, 0], &mori_generators)
            .unwrap()
            .unwrap();
        assert_eq!(face.generators, vec![vec![1, 0]]);

        let missing =
            supporting_mori_face_for_curve_from_normal(&[0, 1], &[1, 1], &mori_generators).unwrap();
        assert!(missing.is_none());
    }

    #[test]
    fn supporting_mori_face_for_curve_rejects_non_supporting_normal() {
        let face = supporting_mori_face_for_curve_from_normal(
            &[0, 1],
            &[1, 0],
            &[vec![1, 0], vec![0, -1]],
        )
        .unwrap();

        assert!(face.is_none());
    }

    #[test]
    fn extremal_mori_ray_separator_certifies_target_ray_exactly() {
        let certificate = check_extremal_mori_ray_separator(
            &[-1, 1],
            &[1, 0],
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
        )
        .unwrap()
        .unwrap();

        assert_eq!(certificate.separator_normal, vec![-1, 1]);
        assert_eq!(certificate.same_ray_generator_count, 1);
        assert_eq!(certificate.zero_other_generator_count, 1);
        assert_eq!(certificate.positive_other_generator_count, 1);
    }

    #[test]
    fn extremal_mori_ray_separator_orients_normal_and_accepts_multiples() {
        let certificate = check_extremal_mori_ray_separator(
            &[1, -1],
            &[1, 0],
            &[vec![2, 0], vec![0, 1], vec![1, 1]],
        )
        .unwrap()
        .unwrap();

        assert_eq!(certificate.separator_normal, vec![-1, 1]);
        assert_eq!(certificate.same_ray_generator_count, 1);
    }

    #[test]
    fn extremal_mori_ray_separator_rejects_decomposable_target() {
        let certificate = check_extremal_mori_ray_separator(
            &[-1, 0],
            &[1, 1],
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
        )
        .unwrap();

        assert!(certificate.is_none());
    }

    #[test]
    fn extremal_mori_ray_separator_finder_finds_exact_separator() {
        let certificate =
            find_extremal_mori_ray_separator(&[1, 0], &[vec![1, 0], vec![0, 1], vec![1, 1]])
                .unwrap()
                .unwrap();

        assert_eq!(certificate.same_ray_generator_count, 1);
        assert_eq!(certificate.zero_other_generator_count, 1);
        assert_eq!(certificate.positive_other_generator_count, 1);
        assert!(
            certificate.separator_normal[0] < 0,
            "separator must pair negatively with the target curve"
        );
    }

    #[test]
    fn extremal_mori_ray_lp_search_finds_exact_separator() {
        let certificate = find_extremal_mori_ray_separator_by_lp_search(
            &[1, 0],
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
            &SupportingMoriFaceLpSearchOptions::default(),
        )
        .unwrap()
        .unwrap();

        assert_eq!(certificate.same_ray_generator_count, 1);
        assert_eq!(certificate.zero_other_generator_count, 1);
        assert_eq!(certificate.positive_other_generator_count, 1);
        assert!(
            certificate.separator_normal[0] < 0,
            "separator must pair negatively with the target curve"
        );
    }

    #[test]
    fn extremal_mori_ray_separator_finder_rejects_decomposable_target() {
        let certificate =
            find_extremal_mori_ray_separator(&[1, 1], &[vec![1, 0], vec![0, 1], vec![1, 1]])
                .unwrap();

        assert!(certificate.is_none());
    }

    #[test]
    fn extremal_mori_ray_lp_search_rejects_decomposable_target() {
        let certificate = find_extremal_mori_ray_separator_by_lp_search(
            &[1, 1],
            &[vec![1, 0], vec![0, 1], vec![1, 1]],
            &SupportingMoriFaceLpSearchOptions::default(),
        )
        .unwrap();

        assert!(certificate.is_none());
    }

    #[test]
    fn extremal_mori_ray_separator_finder_requires_target_ray_in_generators() {
        let certificate =
            find_extremal_mori_ray_separator(&[1, 0], &[vec![0, 1], vec![1, 1]]).unwrap();

        assert!(certificate.is_none());
    }

    #[test]
    fn extremal_mori_ray_separator_checks_inputs() {
        assert!(
            check_extremal_mori_ray_separator(&[1, 0], &[0, 0], &[vec![1, 0]])
                .unwrap_err()
                .to_string()
                .contains("target must not be zero")
        );

        let err =
            check_extremal_mori_ray_separator(&[1, 0], &[1, 0], &[vec![1, 0, 0]]).unwrap_err();
        assert!(err.to_string().contains("Mori generator dimension 3"));
    }

    #[test]
    fn curve_row_span_rank_is_exact_over_rationals() {
        let rows = vec![vec![2, 0, 0], vec![4, 0, 0], vec![0, 3, 0]];

        assert_eq!(curve_row_span_rank(&rows).unwrap(), 2);
        assert!(curve_row_span_rank(&[vec![1], vec![1, 2]]).is_err());
    }

    #[test]
    fn curve_in_rational_row_span_uses_exact_rank_comparison() {
        let rows = vec![vec![2, 0, 0], vec![0, 3, 0]];

        assert!(curve_in_rational_row_span(&[4, 9, 0], &rows).unwrap());
        assert!(!curve_in_rational_row_span(&[1, 0, 1], &rows).unwrap());
    }

    #[test]
    fn curve_in_rational_row_span_handles_empty_span() {
        assert!(curve_in_rational_row_span(&[0, 0], &[]).unwrap());
        assert!(!curve_in_rational_row_span(&[0, 1], &[]).unwrap());
    }

    #[test]
    fn curve_in_rational_row_span_rejects_dimension_mismatch() {
        let err = curve_in_rational_row_span(&[1, 0, 0], &[vec![1, 0]]).unwrap_err();

        assert!(err.to_string().contains("curve dimension 3"));
    }

    #[test]
    fn nilpotent_ray_candidate_detects_first_zero_after_positive_weighted_history() {
        let candidate = detect_apparent_nilpotent_ray_from_gv_multiples(
            &[1, 0],
            &[2, 5],
            4,
            &[(vec![1, 0], Integer::from(3))],
        )
        .unwrap()
        .unwrap();

        assert_eq!(candidate.primitive_ray, vec![1, 0]);
        assert_eq!(candidate.first_vanishing_multiple, 2);
        assert_eq!(candidate.first_vanishing_degree, 4);
        assert_eq!(candidate.weighted_lower_gv_sum, Integer::from(3));
    }

    #[test]
    fn nilpotent_ray_candidate_uses_explicit_zero_gv_values() {
        let candidate = detect_apparent_nilpotent_ray_from_gv_multiples(
            &[0, 1],
            &[3, 2],
            4,
            &[
                (vec![0, 1], Integer::from(2)),
                (vec![0, 2], Integer::from(0)),
            ],
        )
        .unwrap()
        .unwrap();

        assert_eq!(candidate.first_vanishing_multiple, 2);
        assert_eq!(candidate.weighted_lower_gv_sum, Integer::from(2));
    }

    #[test]
    fn nilpotent_ray_candidate_requires_positive_lower_weighted_sum() {
        let no_candidate = detect_apparent_nilpotent_ray_from_gv_multiples(
            &[1, 0],
            &[2, 5],
            4,
            &[(vec![1, 0], Integer::from(-3))],
        )
        .unwrap();

        assert_eq!(no_candidate, None);
    }

    #[test]
    fn nilpotent_ray_candidate_rejects_nonprimitive_or_bad_degree_inputs() {
        let nonprimitive =
            detect_apparent_nilpotent_ray_from_gv_multiples(&[2, 0], &[1, 1], 4, &[]).unwrap_err();
        assert!(nonprimitive.to_string().contains("co-prime"));

        let bad_degree =
            detect_apparent_nilpotent_ray_from_gv_multiples(&[1, 0], &[0, 1], 4, &[]).unwrap_err();
        assert!(bad_degree.to_string().contains("non-positive grading"));
    }

    #[test]
    fn nilpotent_ray_candidate_rejects_bad_gv_table() {
        let dim = detect_apparent_nilpotent_ray_from_gv_multiples(
            &[1, 0],
            &[1, 1],
            4,
            &[(vec![1], Integer::from(1))],
        )
        .unwrap_err();
        assert!(dim.to_string().contains("dimension"));

        let duplicate = detect_apparent_nilpotent_ray_from_gv_multiples(
            &[1, 0],
            &[1, 1],
            4,
            &[
                (vec![1, 0], Integer::from(1)),
                (vec![1, 0], Integer::from(2)),
            ],
        )
        .unwrap_err();
        assert!(duplicate.to_string().contains("conflicting duplicate"));
    }

    #[test]
    fn nilpotent_ray_table_scan_builds_sorted_candidate_set() {
        let candidates = detect_apparent_nilpotent_rays_from_gv_table(
            &[2, 3],
            6,
            &[
                (vec![0, 1], Integer::from(5)),
                (vec![0, 2], Integer::from(0)),
                (vec![1, 0], Integer::from(3)),
            ],
        )
        .unwrap();

        assert_eq!(candidates.len(), 2);
        assert_eq!(candidates[0].primitive_ray, vec![0, 1]);
        assert_eq!(candidates[0].first_vanishing_multiple, 2);
        assert_eq!(candidates[0].weighted_lower_gv_sum, Integer::from(5));
        assert_eq!(candidates[1].primitive_ray, vec![1, 0]);
        assert_eq!(candidates[1].first_vanishing_multiple, 2);
        assert_eq!(candidates[1].weighted_lower_gv_sum, Integer::from(3));
    }

    #[test]
    fn nilpotent_ray_table_scan_reduces_nonprimitive_seed_classes() {
        let candidates = detect_apparent_nilpotent_rays_from_gv_table(
            &[1, 1],
            3,
            &[(vec![2, 0], Integer::from(7))],
        )
        .unwrap();

        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].primitive_ray, vec![1, 0]);
        assert_eq!(candidates[0].first_vanishing_multiple, 3);
        assert_eq!(candidates[0].weighted_lower_gv_sum, Integer::from(28));
    }

    #[test]
    fn nilpotent_ray_table_scan_ignores_nonzero_classes_above_cutoff() {
        let candidates = detect_apparent_nilpotent_rays_from_gv_table(
            &[1, 1],
            4,
            &[(vec![5, 0], Integer::from(7))],
        )
        .unwrap();

        assert!(candidates.is_empty());
    }

    #[test]
    fn nilpotent_ray_table_scan_rejects_bad_tables() {
        let duplicate = detect_apparent_nilpotent_rays_from_gv_table(
            &[1, 1],
            4,
            &[
                (vec![1, 0], Integer::from(1)),
                (vec![1, 0], Integer::from(2)),
            ],
        )
        .unwrap_err();
        assert!(duplicate.to_string().contains("conflicting duplicate"));

        let bad_degree = detect_apparent_nilpotent_rays_from_gv_table(
            &[0, 1],
            4,
            &[(vec![1, 0], Integer::from(1))],
        )
        .unwrap_err();
        assert!(bad_degree.to_string().contains("non-positive grading"));
    }

    #[test]
    fn finite_cutoff_gv_partition_excludes_charges_on_nilpotent_rays() {
        let partition = partition_finite_cutoff_gv_charges_by_nilpotence(
            &[1, 1],
            6,
            &[
                (vec![1, 0], Integer::from(3)),
                (vec![2, 0], Integer::from(4)),
                (vec![0, 1], Integer::from(5)),
                (vec![0, 2], Integer::from(5)),
                (vec![0, 3], Integer::from(5)),
                (vec![0, 4], Integer::from(5)),
                (vec![0, 5], Integer::from(5)),
                (vec![0, 6], Integer::from(5)),
            ],
        )
        .unwrap();

        assert_eq!(partition.nilpotent_rays.len(), 1);
        assert_eq!(partition.nilpotent_rays[0].primitive_ray, vec![1, 0]);
        assert_eq!(
            partition.potent_charges,
            vec![
                (vec![0, 1], Integer::from(5)),
                (vec![0, 2], Integer::from(5)),
                (vec![0, 3], Integer::from(5)),
                (vec![0, 4], Integer::from(5)),
                (vec![0, 5], Integer::from(5)),
                (vec![0, 6], Integer::from(5)),
            ]
        );
    }

    #[test]
    fn finite_cutoff_gv_partition_sorts_potent_charges_and_ignores_above_cutoff() {
        let partition = partition_finite_cutoff_gv_charges_by_nilpotence(
            &[1, 1],
            1,
            &[
                (vec![1, 0], Integer::from(11)),
                (vec![0, 1], Integer::from(5)),
                (vec![5, 0], Integer::from(7)),
            ],
        )
        .unwrap();

        assert!(partition.nilpotent_rays.is_empty());
        assert_eq!(
            partition.potent_charges,
            vec![
                (vec![0, 1], Integer::from(5)),
                (vec![1, 0], Integer::from(11)),
            ]
        );
    }

    #[test]
    fn finite_cutoff_gv_complement_excludes_primitive_rays() {
        let retained = finite_cutoff_gv_charges_excluding_primitive_rays(
            &[1, 1],
            5,
            &[
                (vec![1, 0], Integer::from(3)),
                (vec![2, 0], Integer::from(7)),
                (vec![0, 1], Integer::from(5)),
                (vec![1, 1], Integer::from(11)),
                (vec![0, 5], Integer::from(13)),
                (vec![6, 0], Integer::from(17)),
                (vec![1, 1], Integer::from(11)),
            ],
            &[vec![1, 0], vec![1, 1]],
        )
        .unwrap();

        assert_eq!(
            retained,
            vec![
                (vec![0, 1], Integer::from(5)),
                (vec![0, 5], Integer::from(13)),
            ]
        );
    }

    #[test]
    fn finite_cutoff_gv_complement_rejects_bad_inputs() {
        let nonprimitive =
            finite_cutoff_gv_charges_excluding_primitive_rays(&[1, 1], 3, &[], &[vec![2, 0]])
                .unwrap_err();
        assert!(nonprimitive.to_string().contains("co-prime"));

        let duplicate = finite_cutoff_gv_charges_excluding_primitive_rays(
            &[1, 1],
            3,
            &[
                (vec![1, 0], Integer::from(3)),
                (vec![1, 0], Integer::from(4)),
            ],
            &[],
        )
        .unwrap_err();
        assert!(duplicate.to_string().contains("conflicting duplicate"));

        let bad_degree = finite_cutoff_gv_charges_excluding_primitive_rays(
            &[0, 1],
            3,
            &[(vec![1, 0], Integer::from(3))],
            &[],
        )
        .unwrap_err();
        assert!(bad_degree.to_string().contains("non-positive grading"));
    }

    #[test]
    fn finite_gv_nonzero_degree_slice_points_extracts_sorted_slice() {
        let points = finite_gv_nonzero_degree_slice_points(
            &[2, 1],
            4,
            &[
                (vec![1, 2], Integer::from(7)),
                (vec![2, 0], Integer::from(5)),
                (vec![0, 4], Integer::from(0)),
                (vec![0, 3], Integer::from(11)),
                (vec![1, 2], Integer::from(7)),
            ],
        )
        .unwrap();

        assert_eq!(points, vec![vec![1, 2], vec![2, 0]]);
    }

    #[test]
    fn finite_gv_nonzero_degree_slice_points_rejects_bad_tables() {
        let duplicate = finite_gv_nonzero_degree_slice_points(
            &[1, 1],
            2,
            &[
                (vec![1, 1], Integer::from(3)),
                (vec![1, 1], Integer::from(4)),
            ],
        )
        .unwrap_err();
        assert!(duplicate.to_string().contains("conflicting duplicate"));

        let bad_degree =
            finite_gv_nonzero_degree_slice_points(&[0, 1], 2, &[(vec![1, 0], Integer::from(3))])
                .unwrap_err();
        assert!(bad_degree.to_string().contains("non-positive grading"));
    }

    fn synthetic_nilpotent_candidate(ray: &[i64]) -> NilpotentRayCandidate {
        NilpotentRayCandidate {
            primitive_ray: ray.to_vec(),
            first_vanishing_multiple: 2,
            first_vanishing_degree: 2,
            weighted_lower_gv_sum: Integer::from(1),
        }
    }

    fn synthetic_divergence_check(
        ray: &[i64],
        appears_divergent: Option<bool>,
    ) -> (Vec<i64>, super::NilpotentRayDivergenceCheck) {
        let half_norm = match appears_divergent {
            Some(false) => Some(2),
            Some(true) | None => Some(1),
        };
        let full_norm = match appears_divergent {
            Some(true) => Some(2),
            Some(false) => Some(1),
            None => None,
        };
        let half_slice = NilpotentRayDegreeSlice {
            primitive_ray: ray.to_vec(),
            cutoff_numerator: 1,
            cutoff_denominator: 2,
            slice_multiple: 1,
            slice_degree: 1,
            slice_origin: ray.to_vec(),
        };
        let full_slice = NilpotentRayDegreeSlice {
            primitive_ray: ray.to_vec(),
            cutoff_numerator: 1,
            cutoff_denominator: 1,
            slice_multiple: 2,
            slice_degree: 2,
            slice_origin: ray.iter().map(|value| value * 2).collect(),
        };
        let half_cutoff = NilpotentRaySliceDistance {
            slice: half_slice,
            lattice_offsets: Vec::new(),
            lll_transform: Vec::new(),
            reduced_lattice_offsets: Vec::new(),
            comparison_points: Vec::new(),
            transformed_comparison_offsets: Vec::new(),
            minimum_infinity_norm: half_norm,
        };
        let full_cutoff = NilpotentRaySliceDistance {
            slice: full_slice,
            lattice_offsets: Vec::new(),
            lll_transform: Vec::new(),
            reduced_lattice_offsets: Vec::new(),
            comparison_points: Vec::new(),
            transformed_comparison_offsets: Vec::new(),
            minimum_infinity_norm: full_norm,
        };
        let check =
            nilpotent_ray_divergence_check_from_slice_distances(half_cutoff, full_cutoff).unwrap();
        (ray.to_vec(), check)
    }

    #[test]
    fn nilpotent_ray_degree_slice_computes_half_and_full_cutoff_origins() {
        let half = nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[2, 5], 5, 1, 2)
            .unwrap()
            .unwrap();

        assert_eq!(half.cutoff_numerator, 1);
        assert_eq!(half.cutoff_denominator, 2);
        assert_eq!(half.slice_multiple, 1);
        assert_eq!(half.slice_degree, 2);
        assert_eq!(half.slice_origin, vec![1, 0]);

        let full = nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[2, 5], 5, 1, 1)
            .unwrap()
            .unwrap();

        assert_eq!(full.slice_multiple, 2);
        assert_eq!(full.slice_degree, 4);
        assert_eq!(full.slice_origin, vec![2, 0]);
    }

    #[test]
    fn nilpotent_ray_degree_slice_returns_none_when_ray_misses_cutoff_fraction() {
        let slice =
            nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[5, 1], 4, 1, 2).unwrap();

        assert_eq!(slice, None);
    }

    #[test]
    fn nilpotent_ray_degree_slice_rejects_invalid_inputs() {
        let nonprimitive =
            nilpotent_ray_degree_slice_for_cutoff_fraction(&[2, 0], &[1, 1], 4, 1, 1).unwrap_err();
        assert!(nonprimitive.to_string().contains("co-prime"));

        let bad_fraction =
            nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[1, 1], 4, 0, 1).unwrap_err();
        assert!(bad_fraction.to_string().contains("fraction"));

        let bad_degree =
            nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[0, 1], 4, 1, 1).unwrap_err();
        assert!(bad_degree.to_string().contains("non-positive grading"));
    }

    #[test]
    fn nilpotent_ray_slice_comparison_points_land_on_slice_and_deduplicate() {
        let slice = nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[1, 1], 4, 1, 1)
            .unwrap()
            .unwrap();
        let points = nilpotent_ray_slice_comparison_points(
            &slice,
            &[vec![0, 1], vec![0, 2], vec![1, 1]],
            &[1, 1],
        )
        .unwrap();

        assert_eq!(points.len(), 2);
        assert_eq!(points[0].primitive_ray, vec![0, 1]);
        assert_eq!(points[0].slice_multiple, 4);
        assert_eq!(points[0].slice_point, vec![0, 4]);
        assert_eq!(points[0].offset_from_origin, vec![-4, 4]);
        assert_eq!(points[1].primitive_ray, vec![1, 1]);
        assert_eq!(points[1].slice_multiple, 2);
        assert_eq!(points[1].slice_point, vec![2, 2]);
        assert_eq!(points[1].offset_from_origin, vec![-2, 2]);
    }

    #[test]
    fn nilpotent_ray_slice_comparison_points_skip_nonintegral_slice_hits() {
        let slice = nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[3, 1], 3, 1, 1)
            .unwrap()
            .unwrap();
        let points = nilpotent_ray_slice_comparison_points(&slice, &[vec![1, 1]], &[3, 1]).unwrap();

        assert!(points.is_empty());
    }

    #[test]
    fn nilpotent_ray_slice_comparison_points_reject_bad_comparison_charges() {
        let slice = nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[1, 1], 4, 1, 1)
            .unwrap()
            .unwrap();

        let dim = nilpotent_ray_slice_comparison_points(&slice, &[vec![1]], &[1, 1]).unwrap_err();
        assert!(dim.to_string().contains("dimension"));

        let bad_degree =
            nilpotent_ray_slice_comparison_points(&slice, &[vec![1, 0]], &[0, 1]).unwrap_err();
        assert!(bad_degree.to_string().contains("non-positive grading"));
    }

    #[test]
    fn nilpotent_ray_lll_slice_distance_transforms_comparison_offsets() {
        let slice = nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0, 0], &[1, 1, 1], 4, 1, 1)
            .unwrap()
            .unwrap();

        let distance = nilpotent_ray_lll_reduced_slice_distance(
            &slice,
            &[
                vec![4, 0, 0],
                vec![3, 1, 0],
                vec![3, 0, 1],
                vec![2, 1, 1],
                vec![0, 4, 0],
            ],
            &[vec![0, 1, 0], vec![0, 0, 1], vec![1, 1, 0]],
            &[1, 1, 1],
        )
        .unwrap();

        assert_eq!(distance.slice.slice_origin, vec![4, 0, 0]);
        assert_eq!(
            distance.lattice_offsets,
            vec![
                vec![-4, 4, 0],
                vec![-2, 1, 1],
                vec![-1, 0, 1],
                vec![-1, 1, 0],
            ]
        );
        assert_eq!(
            distance.reduced_lattice_offsets.len(),
            distance.lattice_offsets.len()
        );
        assert_eq!(distance.comparison_points.len(), 3);
        assert_eq!(
            distance.minimum_infinity_norm,
            distance
                .transformed_comparison_offsets
                .iter()
                .map(|offset| offset.iter().map(|value| value.abs()).max().unwrap())
                .min()
        );
    }

    #[test]
    fn nilpotent_ray_lll_slice_distance_reports_no_comparison_hit() {
        let slice = nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[3, 1], 3, 1, 1)
            .unwrap()
            .unwrap();

        let distance = nilpotent_ray_lll_reduced_slice_distance(
            &slice,
            &[vec![1, 0], vec![0, 3]],
            &[vec![1, 1]],
            &[3, 1],
        )
        .unwrap();

        assert!(distance.comparison_points.is_empty());
        assert_eq!(distance.minimum_infinity_norm, None);
    }

    #[test]
    fn nilpotent_ray_lll_slice_distance_rejects_bad_lattice_points() {
        let slice = nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[1, 1], 4, 1, 1)
            .unwrap()
            .unwrap();

        let only_origin =
            nilpotent_ray_lll_reduced_slice_distance(&slice, &[vec![4, 0]], &[], &[1, 1])
                .unwrap_err();
        assert!(only_origin.to_string().contains("nonzero offset"));

        let wrong_degree =
            nilpotent_ray_lll_reduced_slice_distance(&slice, &[vec![3, 0]], &[], &[1, 1])
                .unwrap_err();
        assert!(wrong_degree.to_string().contains("degree"));
    }

    #[test]
    fn nilpotent_ray_divergence_check_compares_half_and_full_distances() {
        let half_slice = nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[1, 1], 4, 1, 2)
            .unwrap()
            .unwrap();
        let full_slice = nilpotent_ray_degree_slice_for_cutoff_fraction(&[1, 0], &[1, 1], 4, 1, 1)
            .unwrap()
            .unwrap();

        let half = NilpotentRaySliceDistance {
            slice: half_slice,
            lattice_offsets: Vec::new(),
            lll_transform: Vec::new(),
            reduced_lattice_offsets: Vec::new(),
            comparison_points: Vec::new(),
            transformed_comparison_offsets: Vec::new(),
            minimum_infinity_norm: Some(2),
        };
        let full = NilpotentRaySliceDistance {
            slice: full_slice,
            lattice_offsets: Vec::new(),
            lll_transform: Vec::new(),
            reduced_lattice_offsets: Vec::new(),
            comparison_points: Vec::new(),
            transformed_comparison_offsets: Vec::new(),
            minimum_infinity_norm: Some(3),
        };

        let check = nilpotent_ray_divergence_check_from_slice_distances(half, full).unwrap();
        assert_eq!(check.appears_divergent, Some(true));
    }

    #[test]
    fn nilpotent_ray_divergence_check_with_explicit_lattices_builds_slices() {
        let check = nilpotent_ray_divergence_check_with_explicit_slice_lattices(
            &[1, 0],
            &[1, 1],
            4,
            &[vec![2, 0], vec![1, 1], vec![0, 2]],
            &[vec![4, 0], vec![3, 1], vec![0, 4]],
            &[vec![0, 1]],
        )
        .unwrap()
        .unwrap();

        assert_eq!(check.half_cutoff.slice.slice_origin, vec![2, 0]);
        assert_eq!(check.full_cutoff.slice.slice_origin, vec![4, 0]);
        assert!(check.half_cutoff.minimum_infinity_norm.is_some());
        assert!(check.full_cutoff.minimum_infinity_norm.is_some());
    }

    #[test]
    fn nilpotent_ray_divergence_check_with_explicit_lattices_returns_none_when_slice_missing() {
        let check = nilpotent_ray_divergence_check_with_explicit_slice_lattices(
            &[1, 0],
            &[5, 1],
            4,
            &[],
            &[],
            &[],
        )
        .unwrap();

        assert_eq!(check, None);
    }

    #[test]
    fn two_pass_nop_classification_builds_f0_and_final_f() {
        let candidates = vec![
            synthetic_nilpotent_candidate(&[1, 0]),
            synthetic_nilpotent_candidate(&[0, 1]),
            synthetic_nilpotent_candidate(&[1, 1]),
        ];
        let classification = classify_nilpotent_rays_from_two_pass_divergence_checks(
            &candidates,
            &[
                synthetic_divergence_check(&[1, 0], Some(true)),
                synthetic_divergence_check(&[0, 1], Some(false)),
                synthetic_divergence_check(&[1, 1], None),
            ],
            &[synthetic_divergence_check(&[1, 0], Some(true))],
        )
        .unwrap();

        assert_eq!(
            classification
                .initial_candidate_nop_rays
                .iter()
                .map(|candidate| candidate.primitive_ray.clone())
                .collect::<Vec<_>>(),
            vec![vec![1, 0]]
        );
        assert_eq!(
            classification
                .nop_rays
                .iter()
                .map(|candidate| candidate.primitive_ray.clone())
                .collect::<Vec<_>>(),
            vec![vec![1, 0]]
        );
        assert_eq!(
            classification.first_pass_inconclusive_rays[0].primitive_ray,
            vec![1, 1]
        );
        assert!(classification.second_pass_inconclusive_rays.is_empty());
    }

    #[test]
    fn two_pass_nop_classification_tracks_second_pass_inconclusive_rays() {
        let candidates = vec![
            synthetic_nilpotent_candidate(&[1, 0]),
            synthetic_nilpotent_candidate(&[0, 1]),
        ];
        let classification = classify_nilpotent_rays_from_two_pass_divergence_checks(
            &candidates,
            &[
                synthetic_divergence_check(&[1, 0], Some(true)),
                synthetic_divergence_check(&[0, 1], Some(true)),
            ],
            &[
                synthetic_divergence_check(&[1, 0], Some(false)),
                synthetic_divergence_check(&[0, 1], None),
            ],
        )
        .unwrap();

        assert_eq!(classification.initial_candidate_nop_rays.len(), 2);
        assert!(classification.nop_rays.is_empty());
        assert_eq!(
            classification.second_pass_inconclusive_rays[0].primitive_ray,
            vec![0, 1]
        );
    }

    #[test]
    fn two_pass_nop_classification_rejects_mismatched_or_missing_checks() {
        let candidates = vec![synthetic_nilpotent_candidate(&[1, 0])];
        let missing_first =
            classify_nilpotent_rays_from_two_pass_divergence_checks(&candidates, &[], &[])
                .unwrap_err();
        assert!(missing_first.to_string().contains("missing a first-pass"));

        let extra_second = classify_nilpotent_rays_from_two_pass_divergence_checks(
            &candidates,
            &[synthetic_divergence_check(&[1, 0], Some(false))],
            &[synthetic_divergence_check(&[1, 0], Some(true))],
        )
        .unwrap_err();
        assert!(extra_second.to_string().contains("non-F0"));

        let mut mismatched = synthetic_divergence_check(&[1, 0], Some(true));
        mismatched.0 = vec![0, 1];
        let mismatched_err = classify_nilpotent_rays_from_two_pass_divergence_checks(
            &candidates,
            &[mismatched],
            &[],
        )
        .unwrap_err();
        assert!(mismatched_err.to_string().contains("does not match"));
    }

    #[test]
    fn finite_gv_table_nop_classification_runs_two_pass_pipeline() {
        let report = classify_nop_rays_from_finite_gv_table(
            &[1, 1],
            4,
            &[
                (vec![1, 0], Integer::from(1)),
                (vec![0, 1], Integer::from(1)),
                (vec![0, 2], Integer::from(1)),
                (vec![0, 3], Integer::from(1)),
                (vec![0, 4], Integer::from(1)),
                (vec![1, 1], Integer::from(1)),
                (vec![2, 2], Integer::from(1)),
            ],
        )
        .unwrap();

        assert_eq!(report.partition.nilpotent_rays.len(), 1);
        assert_eq!(report.partition.nilpotent_rays[0].primitive_ray, vec![1, 0]);
        assert_eq!(report.first_pass_checks.len(), 1);
        assert!(report.first_pass_checks[0].1.appears_divergent.is_some());
        assert_eq!(
            report.second_pass_checks.len(),
            report.classification.initial_candidate_nop_rays.len()
        );
    }

    #[test]
    fn finite_gv_table_nop_classification_requires_nontrivial_slice_lattice() {
        let err =
            classify_nop_rays_from_finite_gv_table(&[1, 1], 2, &[(vec![1, 0], Integer::from(1))])
                .unwrap_err();

        assert!(err.to_string().contains("nonzero offset"));
    }

    #[test]
    fn potent_ray_log_xi_matches_paper_definition() {
        let gv = vec![Integer::from(3), Integer::from(-6), Integer::from(0)];

        let terms = potent_ray_log_xi_terms(&gv, f64_pos!(1.25)).unwrap();

        let expected_1 = 3.0_f64.ln() - 2.0 * PI * 1.25;
        let expected_2 = 6.0_f64.ln() - 2.0 * PI * 2.0 * 1.25;
        assert!((terms[0].unwrap().get() - expected_1).abs() < 1e-12);
        assert!((terms[1].unwrap().get() - expected_2).abs() < 1e-12);
        assert_eq!(terms[2], None);
    }

    #[test]
    fn potent_ray_convergence_reports_negative_decay_slope() {
        let gv = vec![
            Integer::from(3),
            Integer::from(-6),
            Integer::from(27),
            Integer::from(-192),
        ];

        let report = potent_ray_convergence(&gv, f64_pos!(1.37)).unwrap();

        assert_eq!(report.log_xi_terms.len(), 4);
        assert!(
            report.log_xi_slope.unwrap().get() < 0.0,
            "potent-ray log-xi terms should decay at this volume"
        );
    }

    #[test]
    fn one_dimensional_ray_gv_series_rejects_invalid_inputs() {
        assert!(
            compute_one_dimensional_ray_gv_series(&[1], &[1], &[], &Intersection::new(1), 0,)
                .unwrap_err()
                .to_string()
                .contains("at least one multiple")
        );

        assert!(
            compute_one_dimensional_ray_gv_series(&[1], &[0], &[], &Intersection::new(1), 1,)
                .unwrap_err()
                .to_string()
                .contains("non-positive grading degree")
        );

        assert!(
            compute_one_dimensional_ray_gv_series(&[1], &[1, 2], &[], &Intersection::new(1), 1,)
                .unwrap_err()
                .to_string()
                .contains("grading dimension")
        );
    }

    #[test]
    fn provided_generator_ray_gv_series_rejects_invalid_face_context() {
        assert!(
            compute_ray_gv_series_with_provided_generators(
                &[1],
                &[],
                &[1],
                &[],
                &Intersection::new(1),
                1,
            )
            .unwrap_err()
            .to_string()
            .contains("at least one generator")
        );

        assert!(
            compute_ray_gv_series_with_provided_generators(
                &[1, 0],
                &[vec![1]],
                &[1, 1],
                &[],
                &Intersection::new(2),
                1,
            )
            .unwrap_err()
            .to_string()
            .contains("generator dimension")
        );
    }

    #[test]
    fn provided_generator_ray_gv_series_uses_supplied_context_boundary() {
        let err = compute_ray_gv_series_with_provided_generators(
            &[1],
            &[vec![1], vec![2]],
            &[1],
            &[],
            &Intersection::new(1),
            1,
        )
        .unwrap_err();

        assert!(err.to_string().contains("q_matrix is empty"));
    }

    #[test]
    fn one_dimensional_ray_gv_series_propagates_cygv_input_errors() {
        let err = compute_one_dimensional_ray_gv_series(&[1], &[1], &[], &Intersection::new(1), 1)
            .unwrap_err();

        assert!(err.to_string().contains("q_matrix is empty"));
    }

    #[test]
    fn ambient_one_dimensional_ray_gv_series_projects_before_cygv() {
        let err = compute_ambient_one_dimensional_ray_gv_series(
            &[9, 1, 8],
            &[1],
            &[1],
            &[],
            &Intersection::new(1),
            1,
        )
        .unwrap_err();

        assert!(err.to_string().contains("q_matrix is empty"));
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
    fn bounded_decomposition_index_finds_four_term_sums() {
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
                class: vec![2, 2],
                volume: f64_pos!(1.0),
            },
        ];
        let index = BoundedCurveDecompositionIndex::new(&candidates).unwrap();

        assert_eq!(index.find_decomposition(&candidates[2], 3).unwrap(), None);
        assert_eq!(
            index.find_decomposition(&candidates[2], 4).unwrap(),
            Some(vec![vec![0, 1], vec![0, 1], vec![1, 0], vec![1, 0]])
        );
    }

    #[test]
    fn semigroup_decomposition_finds_three_term_sums() {
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

        let decomposition = find_semigroup_decomposition(&candidates[3], &candidates).unwrap();

        assert_eq!(
            decomposition,
            Some(vec![
                CurveDecompositionTerm {
                    class: vec![0, 0, 1],
                    multiplicity: 1,
                },
                CurveDecompositionTerm {
                    class: vec![0, 1, 0],
                    multiplicity: 1,
                },
                CurveDecompositionTerm {
                    class: vec![1, 0, 0],
                    multiplicity: 1,
                },
            ])
        );
    }

    #[test]
    fn semigroup_decomposition_uses_multiplicities() {
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
                class: vec![2, 2],
                volume: f64_pos!(1.0),
            },
        ];

        let decomposition = find_semigroup_decomposition(&candidates[2], &candidates).unwrap();

        assert_eq!(
            decomposition,
            Some(vec![
                CurveDecompositionTerm {
                    class: vec![0, 1],
                    multiplicity: 2,
                },
                CurveDecompositionTerm {
                    class: vec![1, 0],
                    multiplicity: 2,
                },
            ])
        );
    }

    #[test]
    fn semigroup_pruning_removes_multi_term_sums() {
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

        let filtered = remove_semigroup_decomposable_curve_candidates(&candidates).unwrap();

        assert_eq!(filtered, candidates[..3].to_vec());
    }

    #[test]
    fn explicit_curve_pruning_strategy_selects_rule() {
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

        let pair_pruned = prune_decomposable_curve_candidates(
            &candidates,
            CurvePruningStrategy::PairDecomposable,
        )
        .unwrap();
        let semigroup_pruned =
            prune_decomposable_curve_candidates(&candidates, CurvePruningStrategy::FiniteSemigroup)
                .unwrap();

        assert_eq!(CurvePruningStrategy::PairDecomposable.as_str(), "pair");
        assert_eq!(
            CurvePruningStrategy::FiniteSemigroup.as_str(),
            "finite-semigroup"
        );
        assert_eq!(pair_pruned, candidates);
        assert_eq!(semigroup_pruned, candidates[..3].to_vec());
    }

    #[test]
    fn semigroup_decomposition_does_not_use_target_itself() {
        let candidates = vec![ToricCurveCandidate {
            class: vec![1, 0],
            volume: f64_pos!(0.2),
        }];

        let decomposition = find_semigroup_decomposition(&candidates[0], &candidates).unwrap();

        assert_eq!(decomposition, None);
    }

    #[test]
    fn gv_lattice_search_contracts_are_explicit() {
        assert_eq!(
            gv_lattice_search_request(400, None, GvLatticeAugmentation::CytoolsDefault),
            (Some(400), None)
        );
        assert_eq!(
            gv_lattice_search_request(400, Some(7), GvLatticeAugmentation::CytoolsDefault),
            (Some(400), None)
        );
        assert_eq!(
            gv_lattice_search_request(400, Some(7), GvLatticeAugmentation::DegreeBoundedDiagnostic),
            (None, Some(7))
        );
        assert_eq!(
            gv_lattice_search_request(400, Some(7), GvLatticeAugmentation::None),
            (None, None)
        );
    }

    #[test]
    fn cytools_default_lattice_augmentation_uses_cone_grading() {
        let rays = vec![vec![1, 0], vec![0, 1]];
        let supplied_semigroup_grading = vec![7, 11];

        let cytools_default = gv_lattice_augmentation_grading(
            &rays,
            &supplied_semigroup_grading,
            GvLatticeAugmentation::CytoolsDefault,
        )
        .unwrap();
        let diagnostic = gv_lattice_augmentation_grading(
            &rays,
            &supplied_semigroup_grading,
            GvLatticeAugmentation::DegreeBoundedDiagnostic,
        )
        .unwrap();

        assert_eq!(cytools_default, vec![1, 1]);
        assert_eq!(diagnostic, supplied_semigroup_grading);
    }

    #[test]
    fn project_mori_rays_to_basis_normalizes_and_deduplicates() {
        let ambient = vec![vec![0, 2, 4], vec![0, 1, 2], vec![5, 0, 0]];
        let projected = project_mori_cone_cap_rays_to_basis(&ambient, &[1, 2]).unwrap();

        assert_eq!(projected, vec![vec![1, 2]]);
    }

    #[test]
    fn project_mori_rays_to_basis_rejects_bad_basis_index() {
        let err = project_mori_cone_cap_rays_to_basis(&[vec![1, 2]], &[2])
            .expect_err("basis index must be in range");

        assert!(format!("{err}").contains("out of bounds"));
    }

    #[test]
    fn project_mori_rays_to_basis_matrix_normalizes_and_deduplicates() {
        let ambient = vec![vec![2, -3, 5, 7], vec![4, -6, 10, 14], vec![0, 1, 0, 1]];
        let basis_matrix = vec![vec![1, 0, 0, 0], vec![0, 2, -1, 1]];

        let projected =
            project_mori_cone_cap_rays_to_basis_matrix(&ambient, &basis_matrix).unwrap();

        assert_eq!(projected, vec![vec![0, 1], vec![1, -2]]);
    }

    #[test]
    fn project_mori_rays_to_basis_matrix_rejects_inconsistent_rays() {
        let err =
            project_mori_cone_cap_rays_to_basis_matrix(&[vec![1, 2], vec![1, 2, 3]], &[vec![1, 0]])
                .expect_err("ambient Mori rays must have consistent dimensions");

        assert!(err.to_string().contains("inconsistent dimensions"));
    }

    #[test]
    fn divisor_basis_projection_dispatches_vector_basis() {
        let ambient = vec![vec![0, 2, 4], vec![0, 1, 2], vec![5, 0, 0]];

        let projected =
            project_mori_cone_cap_rays_for_divisor_basis(&ambient, DivisorBasis::Indices(&[1, 2]))
                .unwrap();

        assert_eq!(projected, vec![vec![1, 2]]);
    }

    #[test]
    fn divisor_basis_projection_dispatches_matrix_basis() {
        let ambient = vec![vec![2, -3, 5, 7], vec![4, -6, 10, 14], vec![0, 1, 0, 1]];
        let basis_matrix = vec![
            vec![
                Integer::from(1),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
            ],
            vec![
                Integer::from(0),
                Integer::from(2),
                Integer::from(-1),
                Integer::from(1),
            ],
        ];

        let projected = project_mori_cone_cap_rays_for_divisor_basis(
            &ambient,
            DivisorBasis::Matrix {
                standard_basis: &[0, 1],
                basis_matrix: &basis_matrix,
            },
        )
        .unwrap();

        assert_eq!(projected, vec![vec![0, 1], vec![1, -2]]);
    }

    #[test]
    fn gv_divisor_basis_data_builds_matrix_basis_cygv_inputs() {
        let ambient = vec![vec![0, 0, 1, 0], vec![0, 0, 2, 2]];
        let linrels = vec![
            vec![
                Integer::from(1),
                Integer::from(0),
                Integer::from(-1),
                Integer::from(-1),
            ],
            vec![
                Integer::from(0),
                Integer::from(1),
                Integer::from(-2),
                Integer::from(-3),
            ],
        ];
        let standard_basis = vec![2, 3];
        let basis_matrix = vec![
            vec![
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
                Integer::from(1),
            ],
            vec![
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(1),
            ],
        ];

        let data = gv_divisor_basis_data(
            &ambient,
            &linrels,
            DivisorBasis::Matrix {
                standard_basis: &standard_basis,
                basis_matrix: &basis_matrix,
            },
        )
        .unwrap();

        assert_eq!(data.mori_rays, vec![vec![1, 0], vec![2, 1]]);
        assert_eq!(
            data.curve_basis_matrix,
            vec![
                vec![
                    Integer::from(1),
                    Integer::from(2),
                    Integer::from(1),
                    Integer::from(0),
                ],
                vec![
                    Integer::from(0),
                    Integer::from(1),
                    Integer::from(-1),
                    Integer::from(1),
                ],
            ]
        );
        assert_eq!(data.q_matrix, vec![vec![2, 1, 0], vec![1, -1, 1]]);
    }

    #[test]
    fn matrix_basis_quintic_handoff_runs_actual_cygv_degree_one() {
        let linrels = vec![
            vec![
                Integer::from(1),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
            ],
            vec![
                Integer::from(0),
                Integer::from(1),
                Integer::from(-1),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
            ],
            vec![
                Integer::from(0),
                Integer::from(1),
                Integer::from(0),
                Integer::from(-1),
                Integer::from(0),
                Integer::from(0),
            ],
            vec![
                Integer::from(0),
                Integer::from(1),
                Integer::from(0),
                Integer::from(0),
                Integer::from(-1),
                Integer::from(0),
            ],
            vec![
                Integer::from(0),
                Integer::from(1),
                Integer::from(0),
                Integer::from(0),
                Integer::from(0),
                Integer::from(-1),
            ],
        ];
        let standard_basis = vec![1];
        let basis_matrix = vec![vec![
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
            Integer::from(0),
        ]];
        let ambient_mori_rays = vec![vec![0, 1, 1, 1, 1, 1]];
        let data = gv_divisor_basis_data(
            &ambient_mori_rays,
            &linrels,
            DivisorBasis::Matrix {
                standard_basis: &standard_basis,
                basis_matrix: &basis_matrix,
            },
        )
        .expect("matrix basis should build quintic cygv inputs");
        assert_eq!(data.mori_rays, vec![vec![1]]);
        assert_eq!(data.q_matrix, vec![vec![1, 1, 1, 1, 1]]);

        let mut intnums = Intersection::new(1);
        set_intersection_i64(&mut intnums, 0, 0, 0, 5);
        let gvs = compute_gv_invariants_inner(
            &data.mori_rays,
            &[1],
            &data.q_matrix,
            &intnums,
            None,
            Some(1),
            GvLatticeAugmentation::None,
            GvCachePolicy::Disabled,
        )
        .expect("actual cygv should compute quintic degree-one GV");

        assert!(
            gvs.iter()
                .any(|(charge, value)| charge == &[1] && value == &Integer::from(2875)),
            "degree-one quintic GV 2875 missing from {gvs:?}"
        );
    }

    #[test]
    fn cytools_default_quintic_handoff_runs_actual_cygv_degree_one() {
        let mut intnums = Intersection::new(1);
        set_intersection_i64(&mut intnums, 0, 0, 0, 5);

        let gvs = compute_gv_invariants_inner(
            &[vec![1]],
            &[1],
            &[vec![1, 1, 1, 1, 1]],
            &intnums,
            None,
            Some(1),
            GvLatticeAugmentation::CytoolsDefault,
            GvCachePolicy::Disabled,
        )
        .expect("CYTools-style wrapper should compute quintic degree-one GV");

        assert!(
            gvs.iter()
                .any(|(charge, value)| charge == &[1] && value == &Integer::from(2875)),
            "degree-one quintic GV 2875 missing from {gvs:?}"
        );
    }

    #[test]
    fn explicit_quintic_qn_trace_exports_cygv_materialized_polynomial() {
        let mut intnums = Intersection::new(1);
        set_intersection_i64(&mut intnums, 0, 0, 0, 5);

        let traced = compute_gv_invariants_with_explicit_semigroup_qn_trace(
            &[vec![0], vec![1]],
            &[1],
            &[vec![1, 1, 1, 1, 1]],
            &intnums,
        )
        .expect("explicit quintic semigroup should compute degree-one qN trace");

        assert!(
            traced
                .invariants
                .iter()
                .any(|(charge, value)| charge == &[1] && value == &Integer::from(2875)),
            "degree-one quintic GV 2875 missing from {:?}",
            traced.invariants
        );
        assert_eq!(traced.qn_trace.len(), 1);
        assert_eq!(traced.qn_trace[0].degree, 1);
        assert_eq!(traced.qn_trace[0].element, vec![1]);
        assert_eq!(
            traced.qn_trace[0].terms,
            vec![CygvQnTraceTerm {
                monomial_index: 1,
                exponent: vec![1],
                coefficient: "1".to_string(),
            }]
        );
        assert_eq!(
            traced.qn_trace[0].li2_terms,
            vec![CygvQnTraceTerm {
                monomial_index: 1,
                exponent: vec![1],
                coefficient: "1".to_string(),
            }]
        );
        assert_eq!(traced.gv_coefficient_trace.len(), 1);
        assert_eq!(traced.gv_coefficient_trace[0].element, vec![1]);
        assert_eq!(traced.gv_coefficient_trace[0].insertion_index, 0);
        assert_eq!(traced.gv_coefficient_trace[0].pivot_component, 1);
        assert_eq!(
            traced.gv_coefficient_trace[0]
                .instanton_coefficient
                .as_deref(),
            Some("2875")
        );
        assert_eq!(
            traced.gv_coefficient_trace[0].gv_candidate.as_deref(),
            Some("2875")
        );
        assert_eq!(
            traced.gv_coefficient_trace[0]
                .rounded_gv_candidate
                .as_deref(),
            Some("2875")
        );
        assert_eq!(traced.gv_coefficient_trace[0].status, "integer_nonzero_gv");
    }

    #[test]
    fn provided_generator_quintic_qn_trace_matches_normal_wrapper() {
        let mut intnums = Intersection::new(1);
        set_intersection_i64(&mut intnums, 0, 0, 0, 5);

        let normal = compute_gv_invariants_with_provided_generators(
            &[vec![1]],
            &[1],
            &[vec![1, 1, 1, 1, 1]],
            &intnums,
            None,
            Some(1),
        )
        .expect("provided-generator quintic should compute degree-one GV");
        let traced = compute_gv_invariants_with_provided_generators_qn_trace(
            &[vec![1]],
            &[1],
            &[vec![1, 1, 1, 1, 1]],
            &intnums,
            None,
            Some(1),
        )
        .expect("provided-generator quintic should compute degree-one qN trace");

        assert_eq!(traced.invariants, normal);
        assert!(
            traced
                .invariants
                .iter()
                .any(|(charge, value)| charge == &[1] && value == &Integer::from(2875))
        );
        assert_eq!(traced.qn_trace.len(), 1);
        assert_eq!(traced.qn_trace[0].element, vec![1]);
        assert_eq!(traced.qn_trace[0].terms.len(), 1);
        assert_eq!(traced.qn_trace[0].li2_terms.len(), 1);
        assert_eq!(traced.gv_coefficient_trace.len(), 1);
        assert_eq!(traced.gv_coefficient_trace[0].element, vec![1]);
        assert_eq!(traced.gv_coefficient_trace[0].status, "integer_nonzero_gv");
    }

    #[test]
    fn provided_generator_gw_coefficient_trace_exposes_raw_candidate() {
        let mut intnums = Intersection::new(1);
        set_intersection_i64(&mut intnums, 0, 0, 0, 5);

        let trace = compute_gw_coefficient_trace_with_provided_generators(
            &[vec![1]],
            &[1],
            &[vec![1, 1, 1, 1, 1]],
            &intnums,
            None,
            Some(1),
        )
        .expect("provided-generator quintic should expose raw GW coefficient trace");

        assert_eq!(trace.len(), 1);
        assert_eq!(trace[0].element, vec![1]);
        assert_eq!(trace[0].insertion_index, 0);
        assert_eq!(trace[0].pivot_component, 1);
        assert_eq!(trace[0].instanton_coefficient.as_deref(), Some("2875"));
        assert_eq!(trace[0].gv_candidate.as_deref(), Some("2875"));
        assert_eq!(trace[0].rounded_gv_candidate, None);
        assert_eq!(trace[0].status, "nonzero_gw");
    }

    #[test]
    fn explicit_semigroup_gw_coefficient_trace_exposes_raw_candidate() {
        let mut intnums = Intersection::new(1);
        set_intersection_i64(&mut intnums, 0, 0, 0, 5);

        let trace = compute_gw_coefficient_trace_with_explicit_semigroup(
            &[vec![1]],
            &[1],
            &[vec![1, 1, 1, 1, 1]],
            &intnums,
        )
        .expect("explicit quintic semigroup should expose raw GW coefficient trace");

        assert_eq!(trace.len(), 1);
        assert_eq!(trace[0].element, vec![1]);
        assert_eq!(trace[0].insertion_index, 0);
        assert_eq!(trace[0].pivot_component, 1);
        assert_eq!(trace[0].instanton_coefficient.as_deref(), Some("2875"));
        assert_eq!(trace[0].gv_candidate.as_deref(), Some("2875"));
        assert_eq!(trace[0].rounded_gv_candidate, None);
        assert_eq!(trace[0].status, "nonzero_gw");
    }

    #[test]
    fn cyrus_direct_cygv_chain_matches_upstream_quintic_wrapper() {
        let mut cyrus_intnums = Intersection::new(1);
        set_intersection_i64(&mut cyrus_intnums, 0, 0, 0, 5);

        let cyrus_gvs = compute_gv_invariants_inner(
            &[vec![1]],
            &[1],
            &[vec![1, 1, 1, 1, 1]],
            &cyrus_intnums,
            None,
            Some(1),
            GvLatticeAugmentation::None,
            GvCachePolicy::Disabled,
        )
        .expect("Cyrus direct cygv chain should compute quintic GV");
        let cyrus_map = cyrus_gvs
            .into_iter()
            .map(|(charge, value)| (charge, value.to_string()))
            .collect::<BTreeMap<_, _>>();

        let generators = DMatrix::from_column_slice(1, 1, &[1]);
        let grading = RowDVector::from_row_slice(&[1]);
        let q = DMatrix::from_column_slice(5, 1, &[1, 1, 1, 1, 1]);
        let upstream_intnums = HashMap::from([((0usize, 0usize, 0usize), 5i32)]);
        let upstream_gvs = cygv::compute_gv_rat_threefold(
            generators,
            grading,
            Some(1),
            None,
            q,
            Vec::new(),
            upstream_intnums,
            Some(1),
            1000,
        );
        let upstream_map = upstream_gvs
            .into_iter()
            .map(|(charge, value)| (charge.as_slice().to_vec(), value.to_string()))
            .collect::<BTreeMap<_, _>>();

        assert_eq!(cyrus_map, upstream_map);
        assert_eq!(cyrus_map.get(&vec![1]).map(String::as_str), Some("2875"));
    }

    #[test]
    fn matrix_divisor_basis_intersection_matches_dense_tensor_pullback() {
        let mut kappa = Intersection::new(3);
        set_intersection_i64(&mut kappa, 0, 0, 0, 2);
        set_intersection_i64(&mut kappa, 0, 0, 1, 3);
        set_intersection_i64(&mut kappa, 0, 1, 1, 5);
        set_intersection_i64(&mut kappa, 1, 1, 1, 7);
        set_intersection_i64(&mut kappa, 2, 2, 2, 11);
        let basis_matrix = vec![
            vec![Integer::from(1), Integer::from(1), Integer::from(0)],
            vec![Integer::from(0), Integer::from(0), Integer::from(1)],
        ];

        let transformed = intersection_in_matrix_divisor_basis(&kappa, &basis_matrix).unwrap();

        assert_eq!(transformed.dim(), 2);
        assert_eq!(*transformed.get(0, 0, 0).get(), Rational::from(33));
        assert_eq!(*transformed.get(1, 1, 1).get(), Rational::from(11));
        assert_eq!(*transformed.get(0, 0, 1).get(), Rational::from(0));
    }

    #[test]
    fn divisor_basis_intersection_dispatches_vector_and_matrix_shapes() {
        let mut kappa = Intersection::new(3);
        set_intersection_i64(&mut kappa, 0, 0, 2, 13);
        let vector = intersection_in_divisor_basis(&kappa, DivisorBasis::Indices(&[0, 2]))
            .expect("vector basis should dispatch to index filtering");
        assert_eq!(vector.dim(), 2);
        assert_eq!(*vector.get(0, 0, 1).get(), Rational::from(13));

        let standard_basis = vec![0, 2];
        let basis_matrix = vec![
            vec![Integer::from(1), Integer::from(0), Integer::from(0)],
            vec![Integer::from(0), Integer::from(0), Integer::from(1)],
        ];
        let matrix = intersection_in_divisor_basis(
            &kappa,
            DivisorBasis::Matrix {
                standard_basis: &standard_basis,
                basis_matrix: &basis_matrix,
            },
        )
        .expect("matrix basis should dispatch to tensor pullback");

        assert_eq!(matrix.dim(), 2);
        assert_eq!(*matrix.get(0, 0, 1).get(), Rational::from(13));
    }

    #[test]
    fn divisor_basis_intersection_rejects_malformed_shapes() {
        let kappa = Intersection::new(3);
        assert!(intersection_in_divisor_basis(&kappa, DivisorBasis::Indices(&[0, 0])).is_err());
        let basis_matrix = vec![vec![Integer::from(1), Integer::from(0)]];
        assert!(
            intersection_in_matrix_divisor_basis(&kappa, &basis_matrix)
                .expect_err("matrix width should match intersection dimension")
                .to_string()
                .contains("row width")
        );
    }

    fn set_intersection_i64(kappa: &mut Intersection, i: usize, j: usize, k: usize, value: i64) {
        kappa.set(
            i,
            j,
            k,
            crate::types::rational::Rational::<crate::types::tags::Finite>::from_raw(
                Rational::from(value),
            ),
        );
    }

    #[test]
    fn curve_volume_rejects_basis_dimension_mismatch() {
        let err = curve_volume_in_divisor_basis(&[1, 2], &[0, 1], &[f64_finite!(1.0)])
            .expect_err("basis and Kähler coordinate lengths must match");

        assert!(format!("{err}").contains("basis length"));
    }
}
