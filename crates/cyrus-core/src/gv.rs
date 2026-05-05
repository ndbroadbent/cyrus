//! Gopakumar-Vafa invariant computation utilities.
//!
//! Implements the CYTools mori_cone_cap algorithm and grading vector selection,
//! then wires up the `cygv` crate for actual GV computation.
//!
//! Reference: CYTools `calabiyau.py` (mori_cone_cap) and `cone.py` (grading vector).

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
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

const GRADING_CACHE_VERSION: &str = "grading-vector-cytools-lp-v1";
const LATTICE_CACHE_VERSION: &str = "lattice-points-v2";
const CKYZ_ADDITION_TABLE_MAX_ENTRIES: usize = 1_000_000;

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
    let local_coordinates_2d = local_rank_two_coordinates(&relation_points, affine_rank)?;
    let kind = classify_local_toric_circuit(&relation_points, affine_rank);
    Ok(Some(AffineToricCircuitDiagnostic {
        relation_points,
        affine_rank,
        coefficient_counts,
        coefficient_sum,
        coordinate_sum,
        local_charge_basis,
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

    let mut degrees = ckyz_multi_degrees(rank, max_total_degree);
    degrees.sort_by(|lhs, rhs| {
        ckyz_total_degree(lhs)
            .expect("validated degree")
            .cmp(&ckyz_total_degree(rhs).expect("validated degree"))
            .then_with(|| lhs.cmp(rhs))
    });

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
/// substitution, cygv-style `beta - alpha alpha` conversion, and the
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
/// [`compute_ckyz_local_gv_invariants`], but truncates all formal series by the
/// componentwise maximum of `target_degrees`. It is intended for ray checks
/// such as McAllister potent-ray validation, where a total-degree cutoff would
/// compute many irrelevant monomials.
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
    degree_indices: BTreeMap<Vec<usize>, usize>,
    addition_indices: Option<Vec<Vec<Option<usize>>>>,
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
            .collect::<BTreeMap<_, _>>();
        let addition_entries = degrees.len().saturating_mul(degrees.len());
        let addition_indices = if addition_entries <= CKYZ_ADDITION_TABLE_MAX_ENTRIES {
            let mut addition_indices = vec![vec![None; degrees.len()]; degrees.len()];
            for (lhs_index, lhs_degree) in degrees.iter().enumerate() {
                for (rhs_index, rhs_degree) in degrees.iter().enumerate() {
                    let sum_index =
                        ckyz_sum_degree_index(lhs_degree, rhs_degree, rank, &degree_indices)?;
                    addition_indices[lhs_index][rhs_index] = sum_index;
                }
            }
            Some(addition_indices)
        } else {
            None
        };
        Ok(Self {
            rank,
            degrees,
            degree_indices,
            addition_indices,
            max_total_degree,
        })
    }

    fn contains(&self, degree: &[usize]) -> bool {
        self.degree_indices.contains_key(degree)
    }

    fn nonzero_degrees(&self) -> impl Iterator<Item = &Vec<usize>> {
        self.degrees
            .iter()
            .filter(|degree| degree.iter().any(|&entry| entry != 0))
    }

    fn index_of(&self, degree: &[usize]) -> Option<usize> {
        self.degree_indices.get(degree).copied()
    }

    fn sum_index(&self, lhs_index: usize, rhs_index: usize) -> Result<Option<usize>> {
        if let Some(addition_indices) = &self.addition_indices {
            return Ok(addition_indices
                .get(lhs_index)
                .and_then(|row| row.get(rhs_index))
                .copied()
                .flatten());
        }
        let lhs_degree = self.degrees.get(lhs_index).ok_or_else(|| {
            Error::InvalidInput("CKYZ left monomial index is outside the domain".into())
        })?;
        let rhs_degree = self.degrees.get(rhs_index).ok_or_else(|| {
            Error::InvalidInput("CKYZ right monomial index is outside the domain".into())
        })?;
        ckyz_sum_degree_index(lhs_degree, rhs_degree, self.rank, &self.degree_indices)
    }
}

fn ckyz_sum_degree_index(
    lhs_degree: &[usize],
    rhs_degree: &[usize],
    rank: usize,
    degree_indices: &BTreeMap<Vec<usize>, usize>,
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

fn ckyz_coordinate_series_domain(
    rank: usize,
    coordinate: usize,
    domain: &CkyzMonomialDomain,
) -> BTreeMap<Vec<usize>, Rational> {
    let mut out = BTreeMap::new();
    if coordinate < rank && rank == domain.rank {
        let mut degree = vec![0; rank];
        degree[coordinate] = 1;
        if domain.contains(&degree) {
            out.insert(degree, Rational::from(1));
        }
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
    let lhs_terms = ckyz_indexed_domain_terms(lhs, domain)?;
    let rhs_terms = ckyz_indexed_domain_terms(rhs, domain)?;
    let mut out_by_index = vec![None::<Rational>; domain.degrees.len()];
    for (lhs_index, lhs_coefficient) in lhs_terms.iter().copied() {
        for (rhs_index, rhs_coefficient) in rhs_terms.iter().copied() {
            if let Some(product_index) = domain.sum_index(lhs_index, rhs_index)? {
                let entry = out_by_index[product_index].get_or_insert_with(|| Rational::from(0));
                *entry += lhs_coefficient.clone() * rhs_coefficient.clone();
            }
        }
    }
    let mut out = out_by_index
        .into_iter()
        .enumerate()
        .filter_map(|(index, coefficient)| {
            coefficient.and_then(|coefficient| {
                (coefficient != 0).then(|| (domain.degrees[index].clone(), coefficient))
            })
        })
        .collect::<BTreeMap<_, _>>();
    out.retain(|_, coefficient| *coefficient != 0);
    Ok(out)
}

fn ckyz_indexed_domain_terms<'a>(
    series: &'a BTreeMap<Vec<usize>, Rational>,
    domain: &CkyzMonomialDomain,
) -> Result<Vec<(usize, &'a Rational)>> {
    let mut terms = Vec::with_capacity(series.len());
    for (degree, coefficient) in series {
        if *coefficient == 0 {
            continue;
        }
        if degree.len() != domain.rank {
            return Err(Error::InvalidInput(
                "CKYZ series multiplication rank mismatch".into(),
            ));
        }
        if let Some(index) = domain.index_of(degree) {
            terms.push((index, coefficient));
        }
    }
    Ok(terms)
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
        powers.push(next);
    }
    Ok(powers)
}

fn ckyz_series_power_cache_domain(
    series: &BTreeMap<Vec<usize>, Rational>,
    max_exponent: usize,
    domain: &CkyzMonomialDomain,
) -> Result<Vec<BTreeMap<Vec<usize>, Rational>>> {
    let mut powers = Vec::with_capacity(max_exponent + 1);
    let mut identity = BTreeMap::new();
    identity.insert(vec![0; domain.rank], Rational::from(1));
    powers.push(identity);
    for exponent in 1..=max_exponent {
        let next = ckyz_series_mul_domain(&powers[exponent - 1], series, domain)?;
        powers.push(next);
    }
    Ok(powers)
}

fn ckyz_series_exp(
    series: &BTreeMap<Vec<usize>, Rational>,
    rank: usize,
    max_total_degree: usize,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    validate_ckyz_series_has_zero_constant(series, rank, "CKYZ exponential input")?;

    let mut out = BTreeMap::new();
    out.insert(vec![0; rank], Rational::from(1));
    let mut power = out.clone();
    let mut factorial = Integer::from(1);
    for exponent in 1..=max_total_degree {
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

fn ckyz_series_exp_domain(
    series: &BTreeMap<Vec<usize>, Rational>,
    domain: &CkyzMonomialDomain,
) -> Result<BTreeMap<Vec<usize>, Rational>> {
    validate_ckyz_series_has_zero_constant(series, domain.rank, "CKYZ exponential input")?;

    let mut out = BTreeMap::new();
    out.insert(vec![0; domain.rank], Rational::from(1));
    let mut power = out.clone();
    let mut factorial = Integer::from(1);
    for exponent in 1..=domain.max_total_degree {
        power = ckyz_series_mul_domain(&power, series, domain)?;
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
    let mut max_exponents = vec![0usize; rank];
    for (degree, coefficient) in series {
        if *coefficient == 0 || !domain.contains(degree) {
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
            ckyz_series_power_cache_domain(argument, max_exponent, domain)
        })
        .collect::<Result<Vec<_>>>()?;

    let mut out = BTreeMap::new();
    for (degree, coefficient) in series {
        if *coefficient == 0 || !domain.contains(degree) {
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
            monomial =
                ckyz_series_mul_domain(&monomial, &power_caches[coordinate][exponent], domain)?;
            if monomial.is_empty() {
                break;
            }
        }
        ckyz_series_add_scaled_assign(&mut out, &monomial, coefficient.clone());
    }
    Ok(out)
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
    degrees.sort_by(|lhs, rhs| {
        ckyz_total_degree(lhs)
            .expect("validated degree")
            .cmp(&ckyz_total_degree(rhs).expect("validated degree"))
            .then_with(|| lhs.cmp(rhs))
    });
    Ok(degrees)
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

    let mut z_of_q = (0..rank)
        .map(|coordinate| ckyz_coordinate_series_domain(rank, coordinate, domain))
        .collect::<Vec<_>>();
    for _ in 0..domain.max_total_degree {
        let mut next = Vec::with_capacity(rank);
        for (coordinate, correction) in log_period_corrections.iter().enumerate() {
            let correction_at_z = ckyz_series_compose_domain(correction, &z_of_q, domain)?;
            let negative_correction = ckyz_series_scale(&correction_at_z, Rational::from(-1));
            let exp_negative_correction = ckyz_series_exp_domain(&negative_correction, domain)?;
            let q_coordinate = ckyz_coordinate_series_domain(rank, coordinate, domain);
            next.push(ckyz_series_mul_domain(
                &q_coordinate,
                &exp_negative_correction,
                domain,
            )?);
        }
        z_of_q = next;
    }
    Ok(z_of_q)
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

    let mut invariants = BTreeMap::new();
    for degree in target_degrees {
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
    let denom = n * sum_xx - sum_x * sum_x;
    if denom == 0.0 {
        return None;
    }
    F64::<Finite>::new((n * sum_xy - sum_x * sum_y) / denom)
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
            "one-dimensional ray GV series requires a panic=unwind build until cygv panics are converted to Result".into(),
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

    if non_origin_neg_ones == 2 && pos_ones == 3 {
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
    )
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

    let zero_cutoff = RugRational::new();
    let poly_props = cygv::PolynomialProperties::new(&semigroup, &zero_cutoff);
    let (intnum_dict, intnum_idxpairs, n_indices) = cygv::misc::process_int_nums(intnums_map, true)
        .map_err(|e| Error::InvalidInput(format!("cygv intersection preprocessing failed: {e}")))?;

    let n_threads = env::var("CYRUS_GV_THREADS")
        .ok()
        .and_then(|v| v.parse::<u32>().ok())
        .map_or_else(
            || std::thread::available_parallelism().map_or(1, std::num::NonZeroUsize::get),
            |n| {
                if n == 0 { 1 } else { n as usize }
            },
        );
    let pool_size = env::var("CYRUS_GV_POOL_SIZE")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(1000);
    let main_pool = cygv::NumberPool::new(poly_props.zero_cutoff.clone(), pool_size);
    let thread_pools: Vec<_> = (0..n_threads)
        .map(|_| cygv::NumberPool::new(poly_props.zero_cutoff.clone(), pool_size))
        .collect();
    let mut all_pools = (main_pool, thread_pools);
    let nefpart: Vec<DVector<i32>> = Vec::new();

    let fp = cygv::fundamental_period::compute_omega(
        &poly_props,
        &semigroup,
        &q,
        &nefpart,
        &intnum_idxpairs,
        &mut all_pools,
    )
    .map_err(|e| Error::InvalidInput(format!("cygv fundamental period failed: {e}")))?;

    let inst_data = cygv::instanton::compute_instanton_data(
        fp,
        &poly_props,
        &intnum_idxpairs,
        n_indices,
        &intnum_dict,
        true,
        &mut all_pools,
    )
    .map_err(|e| Error::InvalidInput(format!("cygv instanton data failed: {e}")))?;

    let gv = cygv::series_inversion::invert_series::<RugRational, true, true>(
        inst_data,
        &poly_props,
        &mut all_pools,
    )
    .map_err(|e| Error::InvalidInput(format!("cygv series inversion failed: {e}")))?;

    let mut gv_sorted: Vec<_> = gv.into_iter().collect();
    gv_sorted
        .sort_unstable_by_key(|((element_idx, insertion_idx), _)| (*element_idx, *insertion_idx));
    let mut out = Vec::with_capacity(gv_sorted.len());
    for ((element_idx, _), gv_value) in gv_sorted {
        let (numer, denom) = gv_value.into_numer_denom();
        if denom != rug::Integer::from(1) {
            return Err(Error::InvalidInput(format!(
                "explicit GV semigroup produced non-integral invariant at element index {element_idx}: denominator {denom}"
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
    Ok(out)
}

fn compute_gv_invariants_inner(
    rays: &[Vec<i64>],
    grading_vector: &[i64],
    q_matrix: &[Vec<i64>],
    intnums: &Intersection,
    min_points: Option<u32>,
    max_deg: Option<u32>,
    lattice_augmentation: GvLatticeAugmentation,
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
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        LATTICE_CACHE_VERSION.hash(&mut hasher);
        let mut key_rays: Vec<&[i64]> = rays.iter().map(Vec::as_slice).collect();
        key_rays.sort();
        for row in &key_rays {
            row.hash(&mut hasher);
        }
        grading_vector.hash(&mut hasher);
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
    use std::collections::BTreeMap;
    use std::f64::consts::PI;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::{
        BoundedCurveDecompositionIndex, CkyzLocalIntersectionTerm, CkyzMonomialDomain,
        CurveDecompositionTerm, CurvePruningStrategy, GvLatticeAugmentation, LocalToricCircuitKind,
        LocalToricCoordinate2D, OriginCircuitCurveWitness, OriginCircuitRelationPoint,
        ToricCurveCandidate, check_supporting_mori_face_normal, ckyz_series_mul_domain,
        compute_ambient_one_dimensional_ray_gv_series,
        compute_ckyz_flat_prepotential_period_corrections, compute_ckyz_inverse_mirror_map,
        compute_ckyz_local_gv_invariants, compute_ckyz_local_gv_invariants_for_degrees,
        compute_ckyz_local_instanton_potential_corrections,
        compute_ckyz_local_prepotential_period_corrections, compute_ckyz_log_period_corrections,
        compute_grading_vector, compute_gv_invariants_with_explicit_semigroup,
        compute_gv_invariants_with_provided_generators, compute_local_p2_genus_zero_gv_series,
        compute_local_toric_circuit_gv_series, compute_one_dimensional_ray_gv_series,
        compute_ray_gv_series_with_provided_generators, curve_in_rational_row_span,
        curve_row_span_rank, curve_volume_in_divisor_basis, diagnose_affine_toric_circuit,
        dump_mori_rays_cdd, extract_ckyz_local_gv_invariants_from_potential,
        find_pair_decomposition, find_semigroup_decomposition, gv_lattice_search_request,
        load_grading_cache, local_p2_inverse_mirror_map, local_p2_mirror_correction,
        map_basis_gv_invariants_to_ambient, origin_circuit_diagnostic_from_class_and_witnesses,
        potent_ray_convergence, potent_ray_log_xi_terms, project_ambient_curve_to_basis,
        project_ambient_curve_to_basis_matrix, project_mori_cone_cap_rays_to_basis,
        project_mori_cone_cap_rays_to_basis_matrix, prune_decomposable_curve_candidates,
        rank_two_local_charge_model, rank_two_local_support_signature,
        remove_pair_decomposable_curve_candidates, remove_semigroup_decomposable_curve_candidates,
        subcutoff_toric_curve_candidates, supporting_mori_face_for_curve_from_normal,
        supporting_mori_face_from_normal, write_grading_cache,
    };
    use crate::Intersection;
    use crate::lattice::Point;
    use crate::{f64_finite, f64_pos};
    use malachite::Integer;
    use malachite::Rational;

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

        let lhs = BTreeMap::from([(vec![1, 0], Rational::from(2))]);
        let rhs = BTreeMap::from([(vec![0, 2], Rational::from(3))]);
        let product = ckyz_series_mul_domain(&lhs, &rhs, &domain).unwrap();

        assert_eq!(product, BTreeMap::from([(vec![1, 2], Rational::from(6))]));
    }

    #[test]
    fn ckyz_large_monomial_domain_uses_checked_addition_fallback() {
        let domain = CkyzMonomialDomain::componentwise_box(&[1_000]).unwrap();
        assert!(
            domain.addition_indices.is_none(),
            "large CKYZ domains should avoid quadratic addition tables"
        );

        let lhs = BTreeMap::from([(vec![125], Rational::from(2))]);
        let rhs = BTreeMap::from([(vec![175], Rational::from(3))]);
        let product = ckyz_series_mul_domain(&lhs, &rhs, &domain).unwrap();

        assert_eq!(product, BTreeMap::from([(vec![300], Rational::from(6))]));
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
    fn curve_volume_rejects_basis_dimension_mismatch() {
        let err = curve_volume_in_divisor_basis(&[1, 2], &[0, 1], &[f64_finite!(1.0)])
            .expect_err("basis and Kähler coordinate lengths must match");

        assert!(format!("{err}").contains("basis length"));
    }
}
