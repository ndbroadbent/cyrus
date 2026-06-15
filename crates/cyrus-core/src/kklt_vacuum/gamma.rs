//! Orientifold B-field `gamma` derivation from the involution parity.
//!
//! The O7-planes of the simple orientifolds in arXiv:2107.09064 sit on the
//! prime toric divisors whose lattice points share a single nonzero parity
//! class `sigma = p mod 2` (the involution's fixed parity). The declared
//! `c_i = 6` so(8) stacks identify `sigma`; the full O7 set is then EVERY
//! point with that parity, which sets the worldsheet-instanton signs
//! `(-1)^{gamma.q}` for curves meeting them.

use std::collections::HashSet;

use crate::Point;
use crate::types::i64::I64;
use crate::types::tags::{Finite, Pos};

/// Parity class `p mod 2` of the lattice point at `idx`.
fn lattice_point_parity(points: &[Point], idx: usize) -> Vec<i64> {
    points[idx]
        .coords()
        .iter()
        .map(|c| c.rem_euclid(2))
        .collect()
}

/// Derive the involution's fixed parity class `sigma = p mod 2` shared by all
/// declared so(8) (`c_i = 6`) stack divisors.
///
/// Returns `Ok(None)` when there are no so(8) stacks (then the B-field
/// vanishes), `Ok(Some(sigma))` for the shared parity, or `Err` when a KKLT
/// divisor index is out of range, the so(8) stacks disagree on a parity, or
/// the shared parity is trivial.
fn orientifold_involution_parity(
    points: &[Point],
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
) -> Result<Option<Vec<i64>>, String> {
    let ambient_dim = points.len();
    let mut sigma: Option<Vec<i64>> = None;
    for (&divisor_idx, ci) in kklt_basis.iter().zip(c_i.iter()) {
        if divisor_idx >= ambient_dim {
            return Err(format!(
                "KKLT divisor index {divisor_idx} exceeds ambient dimension {ambient_dim}"
            ));
        }
        if ci.get() != 6 {
            continue;
        }
        let parity = lattice_point_parity(points, divisor_idx);
        match sigma.as_ref() {
            None => sigma = Some(parity),
            Some(existing) if *existing == parity => {}
            Some(existing) => {
                return Err(format!(
                    "so(8) stack divisors do not share a single involution parity: {existing:?} vs {parity:?} (divisor {divisor_idx})"
                ));
            }
        }
    }
    if sigma.as_ref().is_some_and(|s| s.iter().all(|&p| p == 0)) {
        return Err(
            "so(8) stack divisors have the trivial parity class; cannot derive the orientifold B-field"
                .to_string(),
        );
    }
    Ok(sigma)
}

/// Cross-check the derived `gamma` against the declared `c_i`: every declared
/// so(8) (`c_i = 6`) divisor must keep its own parity class, and no other
/// declared divisor may carry the O7 involution parity.
fn check_gamma_against_declared_c_i(
    gamma: &[i64],
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
    sigma: &[i64],
) -> Result<(), String> {
    for (&divisor_idx, ci) in kklt_basis.iter().zip(c_i.iter()) {
        if ci.get() == 6 && gamma[divisor_idx] != 1 {
            return Err(format!(
                "declared so(8) divisor {divisor_idx} lost its own parity class; gamma derivation is inconsistent"
            ));
        }
        if ci.get() != 6 && gamma[divisor_idx] == 1 {
            return Err(format!(
                "declared c_i={} divisor {divisor_idx} carries the O7 involution parity {sigma:?}; declared orientifold data is inconsistent",
                ci.get()
            ));
        }
    }
    Ok(())
}

/// Derive the half-integral B-field parity vector `gamma` from the
/// orientifold involution.
///
/// The declared `c_i = 6` so(8) stacks identify the fixed parity `sigma`; the
/// full O7 set is then EVERY point with that parity, including divisors
/// outside the KKLT basis. Returns the all-zero vector when the B-field
/// vanishes (no so(8) stacks).
pub fn compute_b_field_gamma_for_o7_divisors(
    points: &[Point],
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
) -> Result<Vec<I64<Finite>>, String> {
    if kklt_basis.len() != c_i.len() {
        return Err("KKLT basis and c_i length mismatch when computing B-field gamma".to_string());
    }
    let ambient_dim = points.len();
    let Some(sigma) = orientifold_involution_parity(points, kklt_basis, c_i)? else {
        return Ok(vec![I64::<Finite>::new(0); ambient_dim]);
    };

    let kklt_index_set: HashSet<usize> = kklt_basis.iter().copied().collect();
    let mut gamma = vec![0_i64; ambient_dim];
    let mut extra_o7 = Vec::new();
    for (idx, entry) in gamma.iter_mut().enumerate().skip(1) {
        if lattice_point_parity(points, idx) == sigma {
            *entry = 1;
            if !kklt_index_set.contains(&idx) {
                extra_o7.push(idx);
            }
        }
    }

    check_gamma_against_declared_c_i(&gamma, kklt_basis, c_i, &sigma)?;

    eprintln!(
        "[INFO] orientifold involution parity sigma={sigma:?}: O7 divisors={} (beyond KKLT basis: {extra_o7:?})",
        gamma.iter().filter(|&&g| g != 0).count()
    );

    Ok(gamma.into_iter().map(I64::<Finite>::new).collect())
}
