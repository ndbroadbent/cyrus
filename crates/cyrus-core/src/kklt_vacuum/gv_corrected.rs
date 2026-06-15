//! The worldsheet-instanton-corrected KKLT solve: the GV correction inside
//! every Newton linearization (arXiv:2107.09064 SS6).

use malachite::Integer;

use crate::types::f64::F64;
use crate::types::i64::I64;
use crate::types::range::CheckedRange;
use crate::types::tags::{Finite, Pos};
use crate::{apply_finite_f64_basis_transform, is_unimodular};

use super::basis_transform::{basis_change_matrix_between_owned, computed_primal_basis};
use super::{OwnedDivisorBasis, PrimalIntersection};

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
