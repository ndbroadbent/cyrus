//! KKLT vacuum verification: primal-geometry data model and (incrementally)
//! the full Kähler-moduli stabilization entry point.
//!
//! These types and the verification pipeline are being lifted out of the
//! `mcallister_first_principles` binary so the full KKLT stabilization is a
//! callable library routine - used by both the GA verifier and the CLI -
//! rather than logic trapped in a binary `main`. The binary was always
//! intermediate scaffolding; the library is the product.

use malachite::Integer;

use crate::types::f64::F64;
use crate::types::tags::{Finite, Pos};
use crate::{
    DivisorBasis, Intersection, Point, Polytope, Triangulation, apply_finite_f64_basis_transform,
    compute_curve_basis_matrix, divisor_basis_change_matrix,
    effective_prime_divisors_from_curve_basis, heights_to_kahler, is_unimodular,
    scale_divisor_basis_kklt_branch_initialization_to_target,
};

pub mod cone_walk;

/// Primal toric geometry: the reflexive polytope, its FRST heights, the
/// triangulation point set, and the triangulation itself.
pub struct PrimalGeom {
    /// The (canonicalized) primal reflexive polytope.
    pub polytope: Polytope,
    /// FRST heights for the triangulation points.
    pub heights: Vec<F64<Finite>>,
    /// Triangulation points (non-facet-interior lattice points).
    pub triangulation_points: Vec<Point>,
    /// The fine regular star triangulation.
    pub triangulation: Triangulation,
}

/// Primal intersection data: GLSM charges, linear relations, the divisor
/// basis, and the intersection tensor in both the full ambient and the
/// reduced basis.
pub struct PrimalIntersection {
    /// GLSM charge matrix rows (ambient divisor coordinates).
    pub glsm: Vec<Vec<Integer>>,
    /// Linear relations among the triangulation points.
    pub linrels: Vec<Vec<Integer>>,
    /// Divisor basis column indices.
    pub basis: Vec<usize>,
    /// Full ambient intersection tensor.
    pub kappa_full: Intersection,
    /// Intersection tensor reduced to the divisor basis.
    pub kappa_basis: Intersection,
}

/// An owned divisor basis (index list or matrix), mirroring the borrowed
/// [`DivisorBasis`] so it can be stored across pipeline stages.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum OwnedDivisorBasis {
    /// Ambient divisor column indices.
    Indices(Vec<usize>),
    /// Divisor-basis rows in ambient divisor coordinates.
    Matrix {
        /// Standard GLSM column basis used to reduce the matrix block.
        standard_basis: Vec<usize>,
        /// Divisor-basis rows in ambient divisor coordinates, including origin.
        basis_matrix: Vec<Vec<Integer>>,
    },
}

impl OwnedDivisorBasis {
    /// Borrow as a [`DivisorBasis`] for the library APIs.
    #[must_use]
    pub fn as_divisor_basis(&self) -> DivisorBasis<'_> {
        match self {
            Self::Indices(indices) => DivisorBasis::Indices(indices),
            Self::Matrix {
                standard_basis,
                basis_matrix,
            } => DivisorBasis::Matrix {
                standard_basis,
                basis_matrix,
            },
        }
    }

    /// Number of divisor-basis elements.
    #[must_use]
    pub const fn dimension(&self) -> usize {
        match self {
            Self::Indices(indices) => indices.len(),
            Self::Matrix { basis_matrix, .. } => basis_matrix.len(),
        }
    }

    /// Human-readable description of the basis representation.
    #[must_use]
    pub const fn description(&self) -> &'static str {
        match self {
            Self::Indices(_) => "index divisor basis",
            Self::Matrix { .. } => "matrix divisor basis",
        }
    }
}

/// Integer basis-change matrix carrying ambient divisor coordinates from
/// `from_basis` to `to_basis` (a thin owned-basis wrapper over
/// [`divisor_basis_change_matrix`]).
pub fn basis_change_matrix_between_owned(
    glsm: &[Vec<Integer>],
    from_basis: &OwnedDivisorBasis,
    to_basis: &OwnedDivisorBasis,
) -> Result<Vec<Vec<Integer>>, String> {
    divisor_basis_change_matrix(
        glsm,
        from_basis.as_divisor_basis(),
        to_basis.as_divisor_basis(),
    )
    .map_err(|e| e.to_string())
}

/// Cyrus's own computed primal divisor basis (the index basis stored on the
/// intersection data), as an owned basis.
#[must_use]
pub fn computed_primal_basis(intersection: &PrimalIntersection) -> OwnedDivisorBasis {
    OwnedDivisorBasis::Indices(intersection.basis.clone())
}

/// Transform Kähler-basis values from `source_basis` to `target_basis`.
///
/// Returns the values unchanged when the bases coincide; otherwise computes
/// the (necessarily unimodular) basis change and applies it. `log_transform`
/// emits an `[INFO]` line naming the source/target representations.
pub fn transform_kahler_between_owned_divisor_bases(
    glsm: &[Vec<Integer>],
    target_basis: &OwnedDivisorBasis,
    source_basis: &OwnedDivisorBasis,
    values: &[F64<Finite>],
    label: &str,
    log_transform: bool,
) -> Result<Vec<F64<Finite>>, String> {
    if target_basis == source_basis {
        return Ok(values.to_vec());
    }
    let transform = basis_change_matrix_between_owned(glsm, target_basis, source_basis)
        .map_err(|e| format!("failed to compute {label} basis transform: {e}"))?;
    if !is_unimodular(&transform) {
        return Err(format!("{label} basis transform is not unimodular"));
    }
    if log_transform {
        eprintln!(
            "[INFO] transforming {label} from {} source coordinates to {} target coordinates",
            source_basis.description(),
            target_basis.description()
        );
    }
    apply_finite_f64_basis_transform(&transform, values, label).map_err(|e| e.to_string())
}

/// Transform a Kähler point from the production divisor basis into Cyrus's
/// computed primal basis (the coordinates the curve/correction machinery
/// expects).
pub fn transform_production_primal_kahler_to_computed(
    intersection: &PrimalIntersection,
    production_basis: &OwnedDivisorBasis,
    values: &[F64<Finite>],
    label: &str,
) -> Result<Vec<F64<Finite>>, String> {
    transform_kahler_between_owned_divisor_bases(
        &intersection.glsm,
        &computed_primal_basis(intersection),
        production_basis,
        values,
        label,
        true,
    )
}

/// Transform a Kähler point from Cyrus's computed primal basis into the
/// production divisor basis (the inverse of
/// [`transform_production_primal_kahler_to_computed`]).
pub fn transform_computed_primal_kahler_to_production(
    intersection: &PrimalIntersection,
    production_basis: &OwnedDivisorBasis,
    values: &[F64<Finite>],
    label: &str,
) -> Result<Vec<F64<Finite>>, String> {
    transform_kahler_between_owned_divisor_bases(
        &intersection.glsm,
        production_basis,
        &computed_primal_basis(intersection),
        values,
        label,
        true,
    )
}

/// Project the triangulation's FRST heights to a phase-1 KKLT branch
/// initialization in the production divisor basis.
///
/// Maps the certified interior heights through the effective-cone prime
/// divisors to a raw Kähler point, transforms it into the production basis,
/// and scales it onto the phase-1 tau target. Used to seed the KKLT solve
/// from the same chamber the triangulation already certifies.
pub fn height_projected_branch_initialization(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    production_basis: &OwnedDivisorBasis,
    kklt_basis: &[usize],
    tau_phase1: &[F64<Pos>],
) -> Result<Vec<F64<Finite>>, String> {
    let basis_non_origin = intersection
        .basis
        .iter()
        .map(|&idx| {
            idx.checked_sub(1)
                .ok_or_else(|| "height projection basis unexpectedly contains origin".to_string())
        })
        .collect::<Result<Vec<_>, _>>()?;
    let curve_basis = compute_curve_basis_matrix(&intersection.linrels, &intersection.basis)
        .map_err(|e| format!("failed to compute primal curve basis: {e}"))?;
    let prime_divisors = effective_prime_divisors_from_curve_basis(&curve_basis, &basis_non_origin)
        .ok_or_else(|| "failed to extract effective-cone prime divisor rows".to_string())?;
    let raw = heights_to_kahler(&geom.heights, &basis_non_origin, &prime_divisors)
        .ok_or_else(|| "failed to project triangulation heights to Kähler basis".to_string())?;
    let production_raw = transform_computed_primal_kahler_to_production(
        intersection,
        production_basis,
        &raw,
        "height-projected Kähler",
    )?;
    scale_divisor_basis_kklt_branch_initialization_to_target(
        &intersection.kappa_full,
        production_basis.as_divisor_basis(),
        kklt_basis,
        tau_phase1,
        &production_raw,
    )
    .ok_or_else(|| "failed to scale height-projected Kähler point to phase-1 target".to_string())
}
