//! GV divisor-basis data builders for the KKLT vacuum pipeline: thin
//! owned-basis wrappers over [`gv_divisor_basis_data`].

use malachite::Integer;

use crate::{GvDivisorBasisData, gv_divisor_basis_data};

use super::OwnedDivisorBasis;

/// Build GV divisor-basis data (Mori rays, grading, charge matrix) for an
/// owned divisor basis.
pub fn production_gv_basis_data(
    ambient_mori_rays: &[Vec<i64>],
    linrels: &[Vec<Integer>],
    basis: &OwnedDivisorBasis,
    context: &str,
) -> Result<GvDivisorBasisData, String> {
    gv_divisor_basis_data(ambient_mori_rays, linrels, basis.as_divisor_basis())
        .map_err(|e| format!("failed to build {context} GV divisor-basis data: {e}"))
}

/// Build GV divisor-basis data for an index divisor basis.
pub fn vector_gv_basis_data(
    ambient_mori_rays: &[Vec<i64>],
    linrels: &[Vec<Integer>],
    basis: &[usize],
    context: &str,
) -> Result<GvDivisorBasisData, String> {
    production_gv_basis_data(
        ambient_mori_rays,
        linrels,
        &OwnedDivisorBasis::Indices(basis.to_vec()),
        context,
    )
}
