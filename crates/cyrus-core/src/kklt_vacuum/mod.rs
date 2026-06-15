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
use crate::types::i64::I64;
use crate::types::range::CheckedRange;
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

/// Solve the full worldsheet-instanton-corrected KKLT system, with the GV
/// correction inside every Newton linearization.
///
/// Returns the converged Kähler point in Cyrus's computed primal basis.
/// (The earlier scheme - freezing the correction at the previous iterate and
/// re-solving the classical system to a fixed point - is not a contraction
/// for every geometry, e.g. 5-81-3213 cycles with steps ~0.1 regardless of
/// update mixing, while the corrected-system Newton converges directly.)
///
/// When the production basis differs from the computed basis, the corrections
/// are evaluated in computed coordinates: the (t-independent, unimodular)
/// production->computed basis change is validated once up front and then
/// applied infallibly inside each correction/Jacobian closure, with the
/// chain-rule columns folded into the Jacobian.
#[allow(clippy::too_many_arguments)] // one-shot stabilization entry point
pub fn solve_gv_corrected_kklt(
    intersection: &PrimalIntersection,
    production_primal_basis: &OwnedDivisorBasis,
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
    chi_divisor: &[I64<Finite>],
    c_tau: crate::types::physics::CTau,
    small_curve_gvs: &[(Vec<i64>, Integer)],
    gamma: &[I64<Finite>],
    gv_iteration_source_t: &[F64<Finite>],
    kklt_steps: usize,
) -> Result<Vec<F64<Finite>>, String> {
    let zero_correction = vec![F64::<Finite>::ZERO; kklt_basis.len()];
    let base_tau_target =
        crate::kklt::compute_gv_corrected_target_tau(c_i, chi_divisor, c_tau, &zero_correction)
            .ok_or_else(|| "base KKLT target construction failed".to_string())?;
    let production_is_computed = matches!(
        production_primal_basis,
        OwnedDivisorBasis::Indices(indices) if *indices == intersection.basis
    );
    // The production->computed Kähler basis change is t-independent: validate
    // it ONCE here (Err on a non-unimodular/failed transform) so the
    // per-iterate closures below can apply it without failing. `None` means
    // production == computed (the identity, skipped).
    let correction_transform: Option<Vec<Vec<Integer>>> = if production_is_computed {
        None
    } else {
        let transform = basis_change_matrix_between_owned(
            &intersection.glsm,
            &computed_primal_basis(intersection),
            production_primal_basis,
        )
        .map_err(|e| {
            format!("failed to compute GV-corrected solve iterate basis transform: {e}")
        })?;
        if !is_unimodular(&transform) {
            return Err("GV-corrected solve iterate basis transform is not unimodular".to_string());
        }
        Some(transform)
    };
    let computed_t_for_corrections = |t_production: &[F64<Finite>]| -> Vec<F64<Finite>> {
        match correction_transform.as_ref() {
            None => t_production.to_vec(),
            Some(transform) => apply_finite_f64_basis_transform(
                transform,
                t_production,
                "GV-corrected solve iterate",
            )
            .expect("pre-validated unimodular transform applies to a consistent-dim iterate"),
        }
    };
    // Columns of d(t_computed)/d(t_production) for the correction
    // Jacobian chain rule; identity (and skipped) in the default case.
    let correction_transform_columns: Option<Vec<Vec<F64<Finite>>>> = if production_is_computed {
        None
    } else {
        let production_dim = gv_iteration_source_t.len();
        Some(
            (0..production_dim)
                .map(|column| {
                    let mut unit = vec![F64::<Finite>::ZERO; production_dim];
                    unit[column] = crate::f64_finite!(1.0);
                    computed_t_for_corrections(&unit)
                })
                .collect(),
        )
    };
    let gv_correction_closure = |t_production: &[F64<Finite>]| {
        let t_computed = computed_t_for_corrections(t_production);
        crate::kklt::compute_gv_target_correction_for_ambient_curves(
            small_curve_gvs,
            &intersection.basis,
            kklt_basis,
            &t_computed,
            Some(gamma),
        )
    };
    let gv_correction_jacobian_closure = |t_production: &[F64<Finite>]| {
        let t_computed = computed_t_for_corrections(t_production);
        let jac_computed = crate::kklt::compute_gv_target_correction_jacobian_for_ambient_curves(
            small_curve_gvs,
            &intersection.basis,
            kklt_basis,
            &t_computed,
            Some(gamma),
        )?;
        match correction_transform_columns.as_ref() {
            None => Some(jac_computed),
            Some(columns) => Some(
                jac_computed
                    .iter()
                    .map(|row| {
                        columns
                            .iter()
                            .map(|column| {
                                row.iter()
                                    .zip(column.iter())
                                    .fold(F64::<Finite>::ZERO, |acc, (a, b)| acc + *a * *b)
                            })
                            .collect()
                    })
                    .collect(),
            ),
        }
    };
    let gv_corrected = crate::kklt::solve_gv_corrected_divisor_basis_path_following(
        &intersection.kappa_full,
        production_primal_basis.as_divisor_basis(),
        kklt_basis,
        &base_tau_target,
        gv_iteration_source_t,
        CheckedRange::new(0, kklt_steps),
        gv_correction_closure,
        gv_correction_jacobian_closure,
    )
    .ok_or_else(|| "GV-corrected KKLT solve failed".to_string())?;
    if !gv_corrected.converged {
        return Err(format!(
            "GV-corrected KKLT solve did not converge: rel_err={}",
            gv_corrected.relative_error.get()
        ));
    }
    eprintln!(
        "[INFO] GV-corrected KKLT solve converged: rel_err={}",
        gv_corrected.relative_error.get()
    );
    // Map the converged point back to computed coordinates with the same
    // validated transform (identity when production == computed).
    match correction_transform.as_ref() {
        None => Ok(gv_corrected.t),
        Some(transform) => apply_finite_f64_basis_transform(
            transform,
            &gv_corrected.t,
            "GV-corrected Kähler point",
        )
        .map_err(|e| e.to_string()),
    }
}
