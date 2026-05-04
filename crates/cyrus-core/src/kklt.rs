//! KKLT moduli stabilization for Type IIB string compactifications.
//!
//! The KKLT mechanism stabilizes Kähler moduli by balancing flux and
//! non-perturbative terms in the superpotential.
//!
//! ## Target Divisor Volumes
//!
//! The KKLT stabilization targets specific divisor volumes (eq. 5.7):
//! ```text
//! τ_i = (c_i / 2π) × ln(W₀⁻¹)
//! ```
//!
//! Or equivalently:
//! ```text
//! τ_i = c_i / c_τ
//! ```
//!
//! Where:
//! - `c_i` = dual Coxeter numbers (1 for D3-instanton, 6 for O7-plane with so(8))
//! - `c_τ` = 2π / (`g_s` × ln(W₀⁻¹))
//!
//! Reference: arXiv:2107.09064, Section 5

use std::collections::HashMap;
use std::env;
use std::f64::consts::PI;

use malachite::Integer;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

use crate::f64_pos;
use crate::intersection::Intersection;
use crate::types::f64::F64;
use crate::types::i64::I64;
use crate::types::physics::{CTau, DivisorVolume, RelativeError, StringCoupling};
use crate::types::range::CheckedRange;
use crate::types::tags::{Finite, NonNeg, Pos};

/// 2π as a typed positive constant.
const TWO_PI: F64<Pos> = f64_pos!(2.0 * PI);
const DILOG_TOL: f64 = 1e-16;
const DILOG_MAX_TERMS: usize = 100_000;
const ZETA_3: f64 = 1.202_056_903_159_594;

fn kklt_debug_enabled() -> bool {
    env::var_os("CYRUS_KKLT_DEBUG").is_some()
}

/// Compute target divisor volumes `τ_i` = `c_i` / `c_τ`.
///
/// # Arguments
/// * `c_i` - Dual Coxeter numbers (positive integers: 1 for D3, 6 for O7)
/// * `c_tau` - Parameter relating `g_s` to W₀ (positive)
#[must_use]
pub fn compute_target_tau(c_i: &[I64<Pos>], c_tau: CTau) -> Vec<DivisorVolume> {
    c_i.iter()
        .map(|ci| {
            // Pos / Pos = Pos
            ci.to_f64() / c_tau
        })
        .collect()
}

/// Compute corrected target divisor volumes
/// `τ_i = c_i / c_τ + χ(D_i) / 24`.
///
/// This is the zeroth-order McAllister KKLT target after including the divisor
/// Euler-characteristic term from the corrected divisor volume formula.
///
/// Returns `None` if dimensions do not match or any corrected target is not
/// positive.
#[must_use]
pub fn compute_corrected_target_tau(
    c_i: &[I64<Pos>],
    chi_divisor: &[I64<Finite>],
    c_tau: CTau,
) -> Option<Vec<DivisorVolume>> {
    if c_i.len() != chi_divisor.len() {
        return None;
    }

    let twenty_four = f64_pos!(24.0);

    c_i.iter()
        .zip(chi_divisor.iter())
        .map(|(ci, chi)| {
            let kklt_target = ci.to_f64() / c_tau;
            let euler_shift = chi.to_f64() / twenty_four;
            (kklt_target + euler_shift).try_to_pos()
        })
        .collect()
}

/// Compute corrected target divisor volumes including divisor GV corrections:
/// `τ_i = c_i / c_τ + χ(D_i) / 24 - GV_i(t)`.
///
/// This is the classical target that must be hit by
/// `1/2 κ_ijk t^j t^k` when the corrected divisor volume is
/// `T_i = 1/2 κ_ijk t^j t^k - χ(D_i)/24 + GV_i(t)`.
#[must_use]
pub fn compute_gv_corrected_target_tau(
    c_i: &[I64<Pos>],
    chi_divisor: &[I64<Finite>],
    c_tau: CTau,
    gv_correction: &[F64<Finite>],
) -> Option<Vec<DivisorVolume>> {
    if c_i.len() != chi_divisor.len() || c_i.len() != gv_correction.len() {
        return None;
    }

    let twenty_four = f64_pos!(24.0);

    c_i.iter()
        .zip(chi_divisor.iter())
        .zip(gv_correction.iter())
        .map(|((ci, chi), gv)| {
            let kklt_target = ci.to_f64() / c_tau;
            let euler_shift = chi.to_f64() / twenty_four;
            (kklt_target + euler_shift - *gv).try_to_pos()
        })
        .collect()
}

fn real_dilog_series(x: f64) -> Option<f64> {
    if !x.is_finite() || x.abs() >= 1.0 {
        return None;
    }

    let mut sum = 0.0;
    let mut power = x;
    for n in 1..=DILOG_MAX_TERMS {
        let n_f = n as f64;
        let term = power / (n_f * n_f);
        sum += term;
        if term.abs() < DILOG_TOL {
            return Some(sum);
        }
        power *= x;
    }

    None
}

fn real_dilog_real_axis(x: f64) -> Option<f64> {
    if !x.is_finite() || x > 1.0 {
        return None;
    }
    if (x - 1.0).abs() <= f64::EPSILON {
        return Some(PI * PI / 6.0);
    }
    if x.abs() < 1e-100 {
        return Some(0.0);
    }

    if x > 0.5 {
        let one_minus_x = 1.0 - x;
        let reduced = real_dilog_series(one_minus_x)?;
        Some(PI * PI / 6.0 - x.ln() * one_minus_x.ln() - reduced)
    } else if x < -0.5 {
        let reduced_arg = x / (x - 1.0);
        let reduced = real_dilog_series(reduced_arg)?;
        Some(-reduced - 0.5 * (1.0 - x).ln().powi(2))
    } else {
        real_dilog_series(x)
    }
}

fn real_trilog_series(x: f64) -> Option<f64> {
    if !x.is_finite() || x.abs() > 1.0 {
        return None;
    }

    let mut sum = 0.0;
    let mut power = x;
    for n in 1..=DILOG_MAX_TERMS {
        let n_f = n as f64;
        let term = power / (n_f * n_f * n_f);
        sum += term;
        if term.abs() < DILOG_TOL {
            return Some(sum);
        }
        power *= x;
    }

    None
}

fn real_trilog_real_axis(x: f64) -> Option<f64> {
    if !x.is_finite() || x > 1.0 {
        return None;
    }
    if (x - 1.0).abs() <= f64::EPSILON {
        return Some(ZETA_3);
    }
    if x.abs() < 1e-100 {
        return Some(0.0);
    }
    if x < -1.0 {
        let log_neg_x = (-x).ln();
        let reduced = real_trilog_series(1.0 / x)?;
        Some(reduced - log_neg_x.powi(3) / 6.0 - PI * PI * log_neg_x / 6.0)
    } else {
        real_trilog_series(x)
    }
}

fn gv_polylog_argument(q_dot_t: f64, parity: i128) -> Option<f64> {
    if !q_dot_t.is_finite() || q_dot_t == 0.0 {
        return None;
    }
    let sign = if parity.rem_euclid(2) == 0 { 1.0 } else { -1.0 };
    let arg = sign * (-TWO_PI.get() * q_dot_t).exp();
    if !arg.is_finite() || arg > 1.0 {
        return None;
    }
    Some(arg)
}

fn gv_dilog_from_curve_volume(q_dot_t: f64, parity: i128) -> Option<f64> {
    let arg = gv_polylog_argument(q_dot_t, parity)?;
    if arg.abs() < 1e-100 {
        return Some(0.0);
    }

    real_dilog_real_axis(arg)
}

/// Compute the divisor GV correction
/// `GV_i(t) = 1/(2π)^2 Σ_q q_i N_q Li_2((-1)^(γ·q) exp(-2π q·t))`.
///
/// Returns `None` if dimensions do not match, a GV integer cannot be represented
/// as a finite `f64`, or the signed dilogarithm argument lies on the complex
/// branch cut `arg > 1`.
#[must_use]
pub fn compute_gv_target_correction(
    gv_invariants: &[(Vec<i32>, Integer)],
    t: &[F64<Finite>],
    gamma: Option<&[I64<Finite>]>,
) -> Option<Vec<F64<Finite>>> {
    let dim = t.len();
    if dim == 0 || gamma.is_some_and(|g| g.len() != dim) {
        return None;
    }

    let mut correction = vec![0.0f64; dim];
    for (curve, invariant) in gv_invariants {
        if curve.len() != dim {
            return None;
        }

        let q_dot_t = curve
            .iter()
            .zip(t.iter())
            .map(|(&qi, ti)| f64::from(qi) * ti.get())
            .sum::<f64>();

        let parity = gamma.map_or(0_i128, |g| {
            curve
                .iter()
                .zip(g.iter())
                .map(|(&qi, gi)| i128::from(qi) * i128::from(gi.get()))
                .sum::<i128>()
        });
        let dilog = gv_dilog_from_curve_volume(q_dot_t, parity)?;
        if dilog == 0.0 {
            continue;
        }
        let invariant_f = invariant.to_string().parse::<f64>().ok()?;
        if !invariant_f.is_finite() {
            return None;
        }

        for (entry, &qi) in correction.iter_mut().zip(curve.iter()) {
            *entry += f64::from(qi) * invariant_f * dilog;
        }
    }

    let prefactor = 1.0 / (4.0 * PI * PI);
    correction
        .into_iter()
        .map(|value| F64::<Finite>::new(prefactor * value))
        .collect()
}

/// Compute divisor GV corrections for an explicit ambient divisor list.
///
/// `gv_invariants` are curve classes in Kähler-coordinate basis coordinates.
/// `curve_basis_matrix` maps those basis-coordinate curves back to ambient
/// divisor coordinates, matching CYTools' `curve_basis(..., as_matrix=True)`
/// convention: `q_ambient = q_basis * curve_basis_matrix`.
///
/// This is needed for McAllister's mixed-basis KKLT setup, where Kähler
/// coordinates are in `basis` but target divisor volumes are ordered by
/// `kklt_basis`.
#[must_use]
pub fn compute_gv_target_correction_for_divisors(
    gv_invariants: &[(Vec<i32>, Integer)],
    curve_basis_matrix: &[Vec<Integer>],
    kklt_basis: &[usize],
    t: &[F64<Finite>],
    gamma: Option<&[I64<Finite>]>,
) -> Option<Vec<F64<Finite>>> {
    let dim = t.len();
    if dim == 0
        || curve_basis_matrix.len() != dim
        || kklt_basis.is_empty()
        || gamma.is_some_and(|g| g.len() != dim)
    {
        return None;
    }

    let ambient_dim = curve_basis_matrix.first()?.len();
    if ambient_dim == 0
        || curve_basis_matrix
            .iter()
            .any(|row| row.len() != ambient_dim)
        || kklt_basis.iter().any(|&idx| idx >= ambient_dim)
    {
        return None;
    }

    let selected_curve_basis: Vec<Vec<f64>> = curve_basis_matrix
        .iter()
        .map(|row| {
            kklt_basis
                .iter()
                .map(|&idx| row[idx].to_string().parse::<f64>().ok())
                .collect::<Option<Vec<_>>>()
        })
        .collect::<Option<Vec<_>>>()?;
    if selected_curve_basis
        .iter()
        .flatten()
        .any(|value| !value.is_finite())
    {
        return None;
    }

    let mut correction = vec![0.0f64; kklt_basis.len()];
    for (curve, invariant) in gv_invariants {
        if curve.len() != dim {
            return None;
        }

        let q_dot_t = curve
            .iter()
            .zip(t.iter())
            .map(|(&qi, ti)| f64::from(qi) * ti.get())
            .sum::<f64>();

        let parity = gamma.map_or(0_i128, |g| {
            curve
                .iter()
                .zip(g.iter())
                .map(|(&qi, gi)| i128::from(qi) * i128::from(gi.get()))
                .sum::<i128>()
        });
        let dilog = gv_dilog_from_curve_volume(q_dot_t, parity)?;
        if dilog == 0.0 {
            continue;
        }
        let invariant_f = invariant.to_string().parse::<f64>().ok()?;
        if !invariant_f.is_finite() {
            return None;
        }

        for (div_idx, entry) in correction.iter_mut().enumerate() {
            let q_ambient = curve
                .iter()
                .zip(selected_curve_basis.iter())
                .map(|(&qi, basis_row)| f64::from(qi) * basis_row[div_idx])
                .sum::<f64>();
            if !q_ambient.is_finite() {
                return None;
            }
            *entry += q_ambient * invariant_f * dilog;
        }
    }

    let prefactor = 1.0 / (4.0 * PI * PI);
    correction
        .into_iter()
        .map(|value| F64::<Finite>::new(prefactor * value))
        .collect()
}

/// Compute divisor GV corrections from ambient-coordinate curve classes.
///
/// This is useful for validation data such as McAllister's `small_curves.dat`,
/// where curve classes are already expressed by their intersections with all
/// ambient divisors. The Kähler parameters `t` are still coordinates in
/// `basis`, and the output is ordered by `kklt_basis`.
#[must_use]
pub fn compute_gv_target_correction_for_ambient_curves(
    gv_invariants: &[(Vec<i64>, Integer)],
    basis: &[usize],
    kklt_basis: &[usize],
    t: &[F64<Finite>],
    gamma: Option<&[I64<Finite>]>,
) -> Option<Vec<F64<Finite>>> {
    let dim = t.len();
    if dim == 0
        || basis.len() != dim
        || kklt_basis.is_empty()
        || gamma.is_some_and(|g| g.len() != dim)
    {
        return None;
    }

    let ambient_dim = gv_invariants.first()?.0.len();
    if ambient_dim == 0
        || basis.iter().any(|&idx| idx >= ambient_dim)
        || kklt_basis.iter().any(|&idx| idx >= ambient_dim)
    {
        return None;
    }

    let mut correction = vec![0.0f64; kklt_basis.len()];
    let debug = kklt_debug_enabled();
    for (curve_index, (curve, invariant)) in gv_invariants.iter().enumerate() {
        if curve.len() != ambient_dim {
            if debug {
                eprintln!(
                    "[KKLT] ambient GV curve {curve_index} has dimension {}, expected {ambient_dim}",
                    curve.len()
                );
            }
            return None;
        }

        let q_dot_t = basis
            .iter()
            .zip(t.iter())
            .map(|(&idx, ti)| curve[idx] as f64 * ti.get())
            .sum::<f64>();

        let parity = gamma.map_or(0_i128, |g| {
            basis
                .iter()
                .zip(g.iter())
                .map(|(&idx, gi)| i128::from(curve[idx]) * i128::from(gi.get()))
                .sum::<i128>()
        });
        let Some(dilog) = gv_dilog_from_curve_volume(q_dot_t, parity) else {
            if debug {
                let sign = if parity.rem_euclid(2) == 0 { 1.0 } else { -1.0 };
                let arg = sign * (-TWO_PI.get() * q_dot_t).exp();
                eprintln!(
                    "[KKLT] ambient GV curve {curve_index} has invalid q.t={q_dot_t} parity={parity} arg={arg}; first coefficients={:?}",
                    curve.iter().take(16).collect::<Vec<_>>()
                );
            }
            return None;
        };
        if dilog == 0.0 {
            continue;
        }
        let invariant_f = invariant.to_string().parse::<f64>().ok()?;
        if !invariant_f.is_finite() {
            return None;
        }

        for (entry, &divisor_idx) in correction.iter_mut().zip(kklt_basis.iter()) {
            *entry += curve[divisor_idx] as f64 * invariant_f * dilog;
        }
    }

    let prefactor = 1.0 / (4.0 * PI * PI);
    correction
        .into_iter()
        .map(|value| F64::<Finite>::new(prefactor * value))
        .collect()
}

/// Compute the GV instanton contribution to the corrected string-frame volume:
/// `1/(2(2π)^3) Σ_q N_q (Li_3(arg) + 2π(q·t)Li_2(arg))`.
///
/// Curve classes are ambient divisor intersections, while `t` and `gamma` are
/// expressed in `basis` coordinates.
#[must_use]
pub fn compute_gv_volume_correction_for_ambient_curves(
    gv_invariants: &[(Vec<i64>, Integer)],
    basis: &[usize],
    t: &[F64<Finite>],
    gamma: Option<&[I64<Finite>]>,
) -> Option<F64<Finite>> {
    let dim = t.len();
    if dim == 0 || basis.len() != dim || gamma.is_some_and(|g| g.len() != dim) {
        return None;
    }

    let ambient_dim = gv_invariants.first()?.0.len();
    if ambient_dim == 0 || basis.iter().any(|&idx| idx >= ambient_dim) {
        return None;
    }

    let mut correction = 0.0f64;
    for (curve, invariant) in gv_invariants {
        if curve.len() != ambient_dim {
            return None;
        }

        let q_dot_t = basis
            .iter()
            .zip(t.iter())
            .map(|(&idx, ti)| curve[idx] as f64 * ti.get())
            .sum::<f64>();
        let parity = gamma.map_or(0_i128, |g| {
            basis
                .iter()
                .zip(g.iter())
                .map(|(&idx, gi)| i128::from(curve[idx]) * i128::from(gi.get()))
                .sum::<i128>()
        });
        let arg = gv_polylog_argument(q_dot_t, parity)?;
        if arg.abs() < 1e-100 {
            continue;
        }

        let dilog = real_dilog_real_axis(arg)?;
        let trilog = real_trilog_real_axis(arg)?;
        let invariant_f = invariant.to_string().parse::<f64>().ok()?;
        if !invariant_f.is_finite() {
            return None;
        }
        correction += invariant_f * (trilog + TWO_PI.get() * q_dot_t * dilog);
    }

    F64::<Finite>::new(correction / (2.0 * TWO_PI.get().powi(3)))
}

/// Compute `c_τ` = 2π / (`g_s` × ln(W₀⁻¹)).
///
/// This relates the string coupling to the flux superpotential.
///
/// # Arguments
/// * `g_s` - String coupling (positive)
/// * `w0` - Flux superpotential magnitude |W₀| (positive)
#[must_use]
pub fn compute_c_tau(g_s: StringCoupling, w0: F64<Pos>) -> CTau {
    // ln(1/w0) = -ln(w0) = ln(w0.recip())
    // For KKLT we need W0 < 1, so 1/w0 > 1, so ln(1/w0) > 0
    let ln_w0_inv = w0.recip().ln();

    // 2π / (g_s × ln(1/w0))
    let denominator = g_s * ln_w0_inv;

    // This should be positive for valid KKLT
    let denominator_pos = denominator
        .try_to_pos()
        .expect("c_tau denominator must be positive for valid KKLT");
    TWO_PI / denominator_pos
}

/// Compute divisor volumes `τ_i` = (1/2) `κ_ijk` t^j t^k.
///
/// This is the classical formula for divisor volumes in terms of Kähler moduli.
///
/// # Panics
/// Panics if `t.len() != kappa.dim()`.
#[must_use]
pub fn compute_divisor_volumes(kappa: &Intersection, t: &[F64<Finite>]) -> Vec<F64<Finite>> {
    let dim = kappa.dim();
    assert_eq!(t.len(), dim, "dimension mismatch");

    let mut tau = vec![F64::<Finite>::ZERO; dim];

    // For each divisor i, compute τ_i = (1/2) Σ_{j,k} κ_{ijk} t^j t^k
    for ((i_idx, j_idx, k_idx), val) in kappa.iter() {
        let val_f = val.to_f64();

        // Generate all unique permutations of (i,j,k)
        for (a, b, c) in unique_permutations(*i_idx, *j_idx, *k_idx) {
            // Finite * Finite * Finite = Finite
            tau[a] = tau[a] + val_f * t[b] * t[c];
        }
    }

    // Apply the 1/2 factor
    let half = F64::<Finite>::new(0.5).expect("0.5 is finite");
    tau.into_iter().map(|t| half * t).collect()
}

/// Compute the Jacobian `J_ik` = ∂`τ_i/∂t^k` = `κ_ijk` t^j.
///
/// # Panics
/// Panics if `t.len() != kappa.dim()`.
#[must_use]
pub fn compute_jacobian(kappa: &Intersection, t: &[F64<Finite>]) -> Vec<Vec<F64<Finite>>> {
    let dim = kappa.dim();
    assert_eq!(t.len(), dim, "dimension mismatch");

    let mut j = vec![vec![F64::<Finite>::ZERO; dim]; dim];

    for ((i_idx, j_idx, k_idx), val) in kappa.iter() {
        let val_f = val.to_f64();

        for (a, b, c) in unique_permutations(*i_idx, *j_idx, *k_idx) {
            j[a][c] = j[a][c] + val_f * t[b];
        }
    }

    j
}

fn kklt_to_basis_map(
    basis: &[usize],
    kklt_basis: &[usize],
    kappa_all_dim: usize,
) -> Option<Vec<Option<usize>>> {
    if basis.iter().any(|&idx| idx >= kappa_all_dim)
        || kklt_basis.iter().any(|&idx| idx >= kappa_all_dim)
    {
        return None;
    }

    let pt_to_basis_idx: HashMap<usize, usize> = basis
        .iter()
        .copied()
        .enumerate()
        .map(|(idx, pt)| (pt, idx))
        .collect();

    Some(
        kklt_basis
            .iter()
            .map(|pt| pt_to_basis_idx.get(pt).copied())
            .collect(),
    )
}

/// Compute divisor volumes for KKLT divisors when the KKLT divisor set differs
/// from the Kähler coordinate basis.
///
/// `basis` indexes the divisors used as Kähler coordinates. `kklt_basis`
/// indexes the divisors that appear in the non-perturbative superpotential.
/// If a KKLT divisor is not in `basis`, its volume is computed from the full
/// divisor tensor as `τ_D = 1/2 κ_{D,b_j,b_k} t^j t^k`.
#[must_use]
pub fn compute_kklt_divisor_volumes(
    kappa_basis: &Intersection,
    kappa_all: &Intersection,
    basis: &[usize],
    kklt_basis: &[usize],
    t: &[F64<Finite>],
) -> Option<Vec<F64<Finite>>> {
    let basis_dim = basis.len();
    if kappa_basis.dim() != basis_dim || t.len() != basis_dim {
        return None;
    }
    let kklt_to_basis = kklt_to_basis_map(basis, kklt_basis, kappa_all.dim())?;
    let basis_tau = compute_divisor_volumes(kappa_basis, t);
    let half = F64::<Finite>::new(0.5).expect("0.5 is finite");

    let mut tau = Vec::with_capacity(kklt_basis.len());
    for (&pt, basis_idx) in kklt_basis.iter().zip(kklt_to_basis.iter()) {
        if let Some(idx) = *basis_idx {
            tau.push(basis_tau[idx]);
            continue;
        }

        let mut tau_a = F64::<Finite>::ZERO;
        for (j, &bj) in basis.iter().enumerate() {
            for (k, &bk) in basis.iter().enumerate() {
                let kappa = kappa_all.get(pt, bj, bk).to_f64();
                tau_a = tau_a + kappa * t[j] * t[k];
            }
        }
        tau.push(half * tau_a);
    }

    Some(tau)
}

/// Compute the KKLT mixed-basis Jacobian
/// `J[a][k] = ∂τ_a / ∂t^k`.
#[must_use]
pub fn compute_kklt_jacobian(
    kappa_basis: &Intersection,
    kappa_all: &Intersection,
    basis: &[usize],
    kklt_basis: &[usize],
    t: &[F64<Finite>],
) -> Option<Vec<Vec<F64<Finite>>>> {
    let basis_dim = basis.len();
    if kappa_basis.dim() != basis_dim || t.len() != basis_dim {
        return None;
    }
    let kklt_to_basis = kklt_to_basis_map(basis, kklt_basis, kappa_all.dim())?;
    let basis_jacobian = compute_jacobian(kappa_basis, t);

    let mut jacobian = Vec::with_capacity(kklt_basis.len());
    for (&pt, basis_idx) in kklt_basis.iter().zip(kklt_to_basis.iter()) {
        if let Some(idx) = *basis_idx {
            jacobian.push(basis_jacobian[idx].clone());
            continue;
        }

        let mut row = vec![F64::<Finite>::ZERO; basis.len()];
        for (j, &bj) in basis.iter().enumerate() {
            for (k, &bk) in basis.iter().enumerate() {
                let kappa = kappa_all.get(pt, bj, bk).to_f64();
                row[k] = row[k] + kappa * t[j];
            }
        }
        jacobian.push(row);
    }

    Some(jacobian)
}

/// Compute numerical rank and conditioning diagnostics for a KKLT Jacobian.
///
/// The rank threshold follows the usual SVD scale-dependent convention
/// `eps * max(rows, cols) * sigma_max`. The condition number is reported only
/// when the matrix is full numerical rank; otherwise the branch is explicitly
/// rank deficient and the full-rank condition number is undefined.
#[must_use]
pub fn compute_jacobian_diagnostics(
    jacobian: &[Vec<F64<Finite>>],
) -> Option<KkltJacobianDiagnostics> {
    let rows = jacobian.len();
    if rows == 0 {
        return None;
    }
    let cols = jacobian[0].len();
    if cols == 0 || jacobian.iter().any(|row| row.len() != cols) {
        return None;
    }

    let data: Vec<f64> = jacobian
        .iter()
        .flat_map(|row| row.iter().map(|entry| entry.get()))
        .collect();
    let matrix = nalgebra::DMatrix::<f64>::from_row_slice(rows, cols, &data);
    let svd = matrix.svd(false, false);
    if svd.singular_values.iter().any(|value| !value.is_finite()) {
        return None;
    }

    let max_sv = svd.singular_values.iter().copied().fold(0.0f64, f64::max);
    let max_singular_value = F64::<Pos>::new(max_sv)?;
    let tolerance = f64::EPSILON * usize::max(rows, cols) as f64 * max_sv;
    let nonzero: Vec<f64> = svd
        .singular_values
        .iter()
        .copied()
        .filter(|value| *value > tolerance)
        .collect();
    let rank = nonzero.len();
    let min_sv = nonzero.iter().copied().fold(f64::INFINITY, f64::min);
    let min_nonzero_singular_value = F64::<Pos>::new(min_sv)?;
    let max_rank = usize::min(rows, cols);
    let condition_number = if rank == max_rank {
        F64::<Pos>::new(max_sv / min_sv)
    } else {
        None
    };

    Some(KkltJacobianDiagnostics {
        rank,
        max_rank,
        max_singular_value,
        min_nonzero_singular_value,
        condition_number,
    })
}

/// Generate all unique permutations of three indices.
fn unique_permutations(
    i: usize,
    j: usize,
    k: usize,
) -> impl Iterator<Item = (usize, usize, usize)> {
    let perms = [
        (i, j, k),
        (i, k, j),
        (j, i, k),
        (j, k, i),
        (k, i, j),
        (k, j, i),
    ];

    let mut seen = [(false, (0, 0, 0)); 6];
    let mut count = 0;

    for perm in perms {
        let is_dup = seen[..count].iter().any(|(_, p)| *p == perm);
        if !is_dup {
            seen[count] = (true, perm);
            count += 1;
        }
    }

    seen.into_iter()
        .take(count)
        .filter(|(used, _)| *used)
        .map(|(_, perm)| perm)
}

/// Result of KKLT path-following.
#[derive(Debug, Clone)]
pub struct KkltResult {
    /// Kähler moduli t (can be any finite values during optimization).
    pub t: Vec<F64<Finite>>,
    /// Achieved divisor volumes.
    pub tau: Vec<F64<Finite>>,
    /// Target divisor volumes.
    pub tau_target: Vec<DivisorVolume>,
    /// Whether the solver converged.
    pub converged: bool,
    /// Relative error in tau (non-negative).
    pub relative_error: RelativeError,
}

/// A converged KKLT branch produced from an explicit initial Kähler point.
#[derive(Debug, Clone)]
pub struct KkltBranchSolution {
    /// Index of the initial point in the input candidate list.
    pub init_index: usize,
    /// Path-following solve result for this branch.
    pub result: KkltResult,
    /// Classical Calabi-Yau volume `(1/6) κ_ijk t^i t^j t^k`.
    pub classical_volume: F64<Pos>,
    /// Numerical conditioning diagnostics for the final KKLT Jacobian.
    pub jacobian_diagnostics: KkltJacobianDiagnostics,
}

/// Numerical rank/conditioning diagnostics for a KKLT Jacobian.
#[derive(Debug, Clone)]
pub struct KkltJacobianDiagnostics {
    /// Numerical rank from SVD singular values.
    pub rank: usize,
    /// Maximum possible rank, `min(rows, columns)`.
    pub max_rank: usize,
    /// Largest singular value.
    pub max_singular_value: F64<Pos>,
    /// Smallest singular value counted as numerically nonzero.
    pub min_nonzero_singular_value: F64<Pos>,
    /// Full-rank condition number. `None` means the matrix is rank deficient.
    pub condition_number: Option<F64<Pos>>,
}

/// Summary of an explicit KKLT branch-candidate search.
#[derive(Debug, Clone)]
pub struct KkltBranchSearchResult {
    /// Number of input initial points evaluated.
    pub attempted: usize,
    /// Number of candidates whose path-following solve returned a result.
    pub solved: usize,
    /// Number of solved candidates that failed the convergence threshold.
    pub non_converged: usize,
    /// Number of converged candidates with non-positive classical volume.
    pub non_positive_volume: usize,
    /// Converged positive-volume branch solutions.
    pub positive_volume: Vec<KkltBranchSolution>,
}

/// Solve KKLT using path-following.
///
/// This implements the predictor-corrector algorithm from `McAllister` Section 5.2.
///
/// # Arguments
/// * `kappa` - Intersection numbers
/// * `tau_target` - Target divisor volumes (positive)
/// * `t_init` - Initial Kähler moduli guess
/// * `steps` - Checked range for interpolation steps (e.g., `range!(0..100)`)
///
/// # Returns
/// `Some(result)` if converged, `None` if diverged.
pub fn solve_path_following(
    kappa: &Intersection,
    tau_target: &[DivisorVolume],
    t_init: &[F64<Finite>],
    steps: CheckedRange<usize>,
) -> Option<KkltResult> {
    let dim = kappa.dim();
    assert_eq!(t_init.len(), dim);
    assert_eq!(tau_target.len(), dim);
    if steps.end <= steps.start {
        return None;
    }

    solve_path_following_core(
        tau_target,
        t_init,
        steps,
        |t| Some(compute_divisor_volumes(kappa, t)),
        |t| Some(compute_jacobian(kappa, t)),
    )
}

/// Solve KKLT with a simple two-phase path-following branch heuristic.
///
/// Phase 1 follows from a scaled uniform initial point to the uncorrected
/// integer target `τ_i = c_i`. Phase 2 follows from that branch to the
/// corrected target.
///
/// This is not guaranteed to select the exact McAllister branch. The
/// reproduction notes show that `τ(t)=c_i` has multiple branches, and uniform
/// starts can converge to valid but different Kähler points. Exact paper
/// reproduction still needs a first-principles branch-selection strategy
/// rather than assuming this heuristic is canonical.
///
/// Returns `None` if either phase fails or does not converge.
pub fn solve_two_phase_path_following(
    kappa: &Intersection,
    tau_target: &[DivisorVolume],
    c_i: &[I64<Pos>],
    steps: CheckedRange<usize>,
) -> Option<KkltResult> {
    let dim = kappa.dim();
    if dim == 0 || tau_target.len() != dim || c_i.len() != dim || steps.end <= steps.start {
        return None;
    }

    let tau_phase1: Vec<DivisorVolume> = c_i.iter().map(|ci| ci.to_f64()).collect();

    let half = F64::<Finite>::new(0.5).expect("0.5 is finite");
    let mut t_uniform = vec![half; dim];

    let tau_at_uniform = compute_divisor_volumes(kappa, &t_uniform);
    let n_dim = I64::<Pos>::new(dim as i64)?;

    let mean_phase1 = tau_phase1.iter().copied().reduce(|acc, x| acc + x)? / n_dim.to_f64();
    let mean_uniform = tau_at_uniform
        .iter()
        .copied()
        .fold(F64::<Finite>::ZERO, |acc, x| acc + x)
        / n_dim.to_f64();

    if let Some(mean_uniform_pos) = mean_uniform.try_to_pos() {
        let scale = (mean_phase1 / mean_uniform_pos).sqrt();
        for ti in &mut t_uniform {
            *ti = *ti * scale;
        }
    }

    let phase1 = solve_path_following(kappa, &tau_phase1, &t_uniform, steps)?;
    if !phase1.converged {
        return None;
    }

    let phase2 = solve_path_following(kappa, tau_target, &phase1.t, steps)?;
    if phase2.converged { Some(phase2) } else { None }
}

/// Solve KKLT path-following for the McAllister mixed-basis setup.
///
/// This is the direct Rust counterpart of the Python `SparseKappa` adapter:
/// Kähler coordinates are in `basis`, while target divisor volumes are ordered
/// by `kklt_basis`.
pub fn solve_mixed_basis_path_following(
    kappa_basis: &Intersection,
    kappa_all: &Intersection,
    basis: &[usize],
    kklt_basis: &[usize],
    tau_target: &[DivisorVolume],
    t_init: &[F64<Finite>],
    steps: CheckedRange<usize>,
) -> Option<KkltResult> {
    let basis_dim = basis.len();
    if kappa_basis.dim() != basis_dim
        || t_init.len() != basis_dim
        || tau_target.len() != kklt_basis.len()
        || steps.end <= steps.start
    {
        return None;
    }

    solve_path_following_core(
        tau_target,
        t_init,
        steps,
        |t| compute_kklt_divisor_volumes(kappa_basis, kappa_all, basis, kklt_basis, t),
        |t| compute_kklt_jacobian(kappa_basis, kappa_all, basis, kklt_basis, t),
    )
}

/// Evaluate explicit initial Kähler points as KKLT branch candidates.
///
/// This is the reusable branch-selection hook for GA search: the caller owns
/// the generation/evolution of `t_initializations`, while Cyrus evaluates each
/// candidate by solving the mixed-basis KKLT equation and reporting which
/// branches converge with positive classical volume.
#[must_use]
pub fn solve_mixed_basis_path_following_branch_candidates(
    kappa_basis: &Intersection,
    kappa_all: &Intersection,
    basis: &[usize],
    kklt_basis: &[usize],
    tau_target: &[DivisorVolume],
    t_initializations: &[Vec<F64<Finite>>],
    steps: CheckedRange<usize>,
) -> KkltBranchSearchResult {
    enum BranchEvaluation {
        NoResult,
        NonConverged,
        NonPositive,
        Positive(KkltBranchSolution),
    }

    let evaluations: Vec<BranchEvaluation> = t_initializations
        .par_iter()
        .enumerate()
        .map(|(init_index, t_init)| {
            let Some(result) = solve_mixed_basis_path_following(
                kappa_basis,
                kappa_all,
                basis,
                kklt_basis,
                tau_target,
                t_init,
                steps,
            ) else {
                return BranchEvaluation::NoResult;
            };
            if !result.converged {
                return BranchEvaluation::NonConverged;
            }
            let Some(classical_volume) = classical_volume_for_branch(kappa_basis, &result.t) else {
                return BranchEvaluation::NonPositive;
            };
            let jacobian =
                compute_kklt_jacobian(kappa_basis, kappa_all, basis, kklt_basis, &result.t)
                    .expect("converged mixed-basis KKLT branch has a valid final Jacobian");
            let jacobian_diagnostics = compute_jacobian_diagnostics(&jacobian)
                .expect("converged mixed-basis KKLT branch has finite Jacobian diagnostics");
            BranchEvaluation::Positive(KkltBranchSolution {
                init_index,
                result,
                classical_volume,
                jacobian_diagnostics,
            })
        })
        .collect();

    let mut out = KkltBranchSearchResult {
        attempted: t_initializations.len(),
        solved: 0,
        non_converged: 0,
        non_positive_volume: 0,
        positive_volume: Vec::new(),
    };

    for evaluation in evaluations {
        match evaluation {
            BranchEvaluation::NoResult => {}
            BranchEvaluation::NonConverged => {
                out.solved += 1;
                out.non_converged += 1;
            }
            BranchEvaluation::NonPositive => {
                out.solved += 1;
                out.non_positive_volume += 1;
            }
            BranchEvaluation::Positive(solution) => {
                out.solved += 1;
                out.positive_volume.push(solution);
            }
        }
    }

    out
}

/// Generate deterministic KKLT branch initializations from a random seed.
///
/// The patterns mirror the Python branch-search experiments: scaled uniform
/// starts, positive random starts, sparse large directions, log-uniform starts,
/// and contiguous cluster starts. Each candidate is rescaled, when possible,
/// so the mean absolute initial divisor volume is comparable to the target
/// divisor-volume scale. The caller still evaluates every candidate explicitly;
/// this function does not choose or discard a branch.
#[must_use]
pub fn generate_scaled_kklt_branch_initializations(
    kappa_basis: &Intersection,
    kappa_all: &Intersection,
    basis: &[usize],
    kklt_basis: &[usize],
    tau_target: &[DivisorVolume],
    count: usize,
    seed: u64,
) -> Option<Vec<Vec<F64<Finite>>>> {
    let dim = basis.len();
    if dim == 0
        || kappa_basis.dim() != dim
        || tau_target.len() != kklt_basis.len()
        || kklt_basis.iter().any(|&idx| idx >= kappa_all.dim())
    {
        return None;
    }

    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let mut out = Vec::with_capacity(count);
    let uniform_scales = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0];
    for scale in uniform_scales.into_iter().take(count) {
        let raw = vec![F64::<Finite>::new(scale).expect("finite scale"); dim];
        out.push(scale_branch_initialization_to_target(
            kappa_basis,
            kappa_all,
            basis,
            kklt_basis,
            tau_target,
            &raw,
        )?);
    }

    while out.len() < count {
        let pattern = rng.gen_range(0..5);
        let raw = match pattern {
            0 => {
                let scale = 10.0_f64.powf(rng.gen_range(-0.5..1.5));
                vec![F64::<Finite>::new(scale)?; dim]
            }
            1 => {
                let scale = 10.0_f64.powf(rng.gen_range(-0.5..1.5));
                (0..dim)
                    .map(|_| F64::<Finite>::new(standard_normal(&mut rng).abs() * scale + 0.1))
                    .collect::<Option<Vec<_>>>()?
            }
            2 => {
                let mut t = vec![F64::<Finite>::new(0.1).expect("finite"); dim];
                let n_large = rng.gen_range(1..=usize::max(1, dim / 10));
                for _ in 0..n_large {
                    let idx = rng.gen_range(0..dim);
                    t[idx] = F64::<Finite>::new(rng.gen_range(1.0..25.0))?;
                }
                t
            }
            3 => (0..dim)
                .map(|_| F64::<Finite>::new(10.0_f64.powf(rng.gen_range(-0.5..1.5))))
                .collect::<Option<Vec<_>>>()?,
            _ => {
                let mut t = vec![F64::<Finite>::new(0.5).expect("finite"); dim];
                let cluster_size = usize::max(1, dim / 4);
                let start = if cluster_size >= dim {
                    0
                } else {
                    rng.gen_range(0..=(dim - cluster_size))
                };
                for entry in t.iter_mut().skip(start).take(cluster_size) {
                    *entry = F64::<Finite>::new(rng.gen_range(5.0..20.0))?;
                }
                t
            }
        };
        out.push(scale_branch_initialization_to_target(
            kappa_basis,
            kappa_all,
            basis,
            kklt_basis,
            tau_target,
            &raw,
        )?);
    }

    Some(out)
}

fn scale_branch_initialization_to_target(
    kappa_basis: &Intersection,
    kappa_all: &Intersection,
    basis: &[usize],
    kklt_basis: &[usize],
    tau_target: &[DivisorVolume],
    t: &[F64<Finite>],
) -> Option<Vec<F64<Finite>>> {
    let tau_init = compute_kklt_divisor_volumes(kappa_basis, kappa_all, basis, kklt_basis, t)?;
    let mean_tau_init =
        tau_init.iter().map(|entry| entry.get().abs()).sum::<f64>() / tau_init.len() as f64;
    let mean_tau_target = tau_target
        .iter()
        .map(|entry| entry.get().abs())
        .sum::<f64>()
        / tau_target.len() as f64;
    if !mean_tau_init.is_finite()
        || !mean_tau_target.is_finite()
        || mean_tau_init <= 0.0
        || mean_tau_target <= 0.0
    {
        return Some(t.to_vec());
    }

    let scale = (mean_tau_target / mean_tau_init).sqrt();
    if !scale.is_finite() || !(0.01..100.0).contains(&scale) {
        return Some(t.to_vec());
    }
    let scale = F64::<Finite>::new(scale)?;
    t.iter().map(|entry| Some(*entry * scale)).collect()
}

/// Rescale an explicit mixed-basis KKLT branch initialization to the target
/// divisor-volume scale.
///
/// This is the public counterpart of the scaling used by
/// [`generate_scaled_kklt_branch_initializations`].  It is useful for
/// deterministic initial points, such as Kähler coordinates projected from
/// secondary-fan heights, before they are evaluated alongside random branch
/// candidates.
#[must_use]
pub fn scale_mixed_basis_kklt_branch_initialization_to_target(
    kappa_basis: &Intersection,
    kappa_all: &Intersection,
    basis: &[usize],
    kklt_basis: &[usize],
    tau_target: &[DivisorVolume],
    t: &[F64<Finite>],
) -> Option<Vec<F64<Finite>>> {
    scale_branch_initialization_to_target(kappa_basis, kappa_all, basis, kklt_basis, tau_target, t)
}

fn standard_normal<R: Rng + ?Sized>(rng: &mut R) -> f64 {
    let u1 = rng.gen_range(f64::MIN_POSITIVE..1.0);
    let u2 = rng.gen_range(0.0..1.0);
    (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
}

fn classical_volume_for_branch(kappa_basis: &Intersection, t: &[F64<Finite>]) -> Option<F64<Pos>> {
    Some(kappa_basis.contract_triple_finite(t)? / f64_pos!(6.0))
}

fn solve_path_following_core<TauFn, JacobianFn>(
    tau_target: &[DivisorVolume],
    t_init: &[F64<Finite>],
    steps: CheckedRange<usize>,
    compute_tau: TauFn,
    compute_jacobian_fn: JacobianFn,
) -> Option<KkltResult>
where
    TauFn: Fn(&[F64<Finite>]) -> Option<Vec<F64<Finite>>>,
    JacobianFn: Fn(&[F64<Finite>]) -> Option<Vec<Vec<F64<Finite>>>>,
{
    let mut t = t_init.to_vec();
    let tau_init = compute_tau(&t)?;
    let n_steps = I64::<Pos>::new((steps.end - steps.start) as i64)?;

    for m in steps.iter_non_neg() {
        let step_index = I64::<Pos>::new(m.get() - steps.start as i64 + 1)?;
        let alpha = step_index.to_f64() / n_steps.to_f64();
        let one_minus_alpha = f64_pos!(1.0) - alpha;

        let tau_step: Vec<F64<Finite>> = tau_init
            .iter()
            .zip(tau_target.iter())
            .map(|(ti, tt)| one_minus_alpha * *ti + alpha * *tt)
            .collect();

        let tau_current = compute_tau(&t)?;
        let delta_tau: Vec<F64<Finite>> = tau_step
            .iter()
            .zip(tau_current.iter())
            .map(|(ts, tc)| *ts - *tc)
            .collect();

        let j = compute_jacobian_fn(&t)?;
        let epsilon = solve_least_squares(&j, &delta_tau)?;
        apply_delta(&mut t, &epsilon)?;

        for _ in 0..3 {
            let tau_current = compute_tau(&t)?;
            let residual = relative_l2_error(&tau_current, &tau_step)?;
            if residual < 1e-8 {
                break;
            }
            let delta_tau: Vec<F64<Finite>> = tau_step
                .iter()
                .zip(tau_current.iter())
                .map(|(ts, tc)| *ts - *tc)
                .collect();
            let j = compute_jacobian_fn(&t)?;
            let mut correction = solve_least_squares(&j, &delta_tau)?;
            damp_delta_to_max_norm(&mut correction, 1.0)?;
            apply_delta(&mut t, &correction)?;
        }

        let divergence_threshold = F64::<NonNeg>::new(1e6).expect("threshold is non-negative");
        if t.iter().any(|ti| ti.abs() > divergence_threshold) {
            return None;
        }
    }

    let tau = compute_tau(&t)?;
    let n_targets = I64::<Pos>::new(tau_target.len() as i64)?;
    let error_sq: F64<NonNeg> = tau
        .iter()
        .zip(tau_target.iter())
        .map(|(ta, tt)| (*ta - *tt).square())
        .fold(F64::<NonNeg>::ZERO, |acc, x| acc + x);
    let error = (error_sq / n_targets.to_f64()).sqrt();
    let sum_target: F64<Pos> = tau_target.iter().copied().reduce(|acc, x| acc + x)?;
    let mean_target = sum_target / n_targets.to_f64();
    let relative_error = error / mean_target;
    let convergence_threshold = F64::<NonNeg>::new(0.001).expect("threshold is non-negative");
    let converged = relative_error < convergence_threshold;

    Some(KkltResult {
        t,
        tau,
        tau_target: tau_target.to_vec(),
        converged,
        relative_error,
    })
}

fn apply_delta(t: &mut [F64<Finite>], delta: &[F64<Finite>]) -> Option<()> {
    if t.len() != delta.len() {
        return None;
    }
    for (ti, di) in t.iter_mut().zip(delta.iter()) {
        *ti = *ti + *di;
    }
    Some(())
}

fn damp_delta_to_max_norm(delta: &mut [F64<Finite>], max_norm: f64) -> Option<()> {
    let norm_sq = delta
        .iter()
        .map(|entry| entry.get() * entry.get())
        .sum::<f64>();
    if !norm_sq.is_finite() {
        return None;
    }
    let norm = norm_sq.sqrt();
    if norm <= max_norm {
        return Some(());
    }
    let scale = F64::<Finite>::new(max_norm / norm)?;
    for entry in delta {
        *entry = *entry * scale;
    }
    Some(())
}

fn relative_l2_error(actual: &[F64<Finite>], target: &[F64<Finite>]) -> Option<f64> {
    if actual.len() != target.len() || actual.is_empty() {
        return None;
    }
    let numerator = actual
        .iter()
        .zip(target.iter())
        .map(|(a, b)| {
            let delta = a.get() - b.get();
            delta * delta
        })
        .sum::<f64>()
        .sqrt();
    let denominator = target
        .iter()
        .map(|value| value.get() * value.get())
        .sum::<f64>()
        .sqrt();
    if !numerator.is_finite() || !denominator.is_finite() || denominator == 0.0 {
        return None;
    }
    Some(numerator / denominator)
}

/// Solve KKLT with a two-phase branch heuristic for mixed `basis`/`kklt_basis`
/// divisor sets.
pub fn solve_two_phase_mixed_basis_path_following(
    kappa_basis: &Intersection,
    kappa_all: &Intersection,
    basis: &[usize],
    kklt_basis: &[usize],
    tau_target: &[DivisorVolume],
    c_i: &[I64<Pos>],
    steps: CheckedRange<usize>,
) -> Option<KkltResult> {
    let dim = basis.len();
    if dim == 0
        || kappa_basis.dim() != dim
        || kklt_basis.len() != tau_target.len()
        || kklt_basis.len() != c_i.len()
        || steps.end <= steps.start
    {
        return None;
    }

    let tau_phase1: Vec<DivisorVolume> = c_i.iter().map(|ci| ci.to_f64()).collect();
    let half = F64::<Finite>::new(0.5).expect("0.5 is finite");
    let mut t_uniform = vec![half; dim];
    let tau_at_uniform =
        compute_kklt_divisor_volumes(kappa_basis, kappa_all, basis, kklt_basis, &t_uniform)?;
    let n_kklt = I64::<Pos>::new(kklt_basis.len() as i64)?;

    let mean_phase1 = tau_phase1.iter().copied().reduce(|acc, x| acc + x)? / n_kklt.to_f64();
    let mean_uniform = tau_at_uniform
        .iter()
        .copied()
        .fold(F64::<Finite>::ZERO, |acc, x| acc + x)
        / n_kklt.to_f64();

    if let Some(mean_uniform_pos) = mean_uniform.try_to_pos() {
        let scale = (mean_phase1 / mean_uniform_pos).sqrt();
        for ti in &mut t_uniform {
            *ti = *ti * scale;
        }
    }

    let phase1 = solve_mixed_basis_path_following(
        kappa_basis,
        kappa_all,
        basis,
        kklt_basis,
        &tau_phase1,
        &t_uniform,
        steps,
    )?;
    if !phase1.converged {
        return None;
    }

    let phase2 = solve_mixed_basis_path_following(
        kappa_basis,
        kappa_all,
        basis,
        kklt_basis,
        tau_target,
        &phase1.t,
        steps,
    )?;
    if phase2.converged { Some(phase2) } else { None }
}

/// Solve the mixed-basis KKLT system with a fixed-point update for divisor GV
/// target corrections.
///
/// The first solve uses the zeroth-order target
/// `c_i/c_tau + χ(D_i)/24` to select the branch. Each subsequent iteration
/// recomputes `GV_i(t)`, updates the classical target
/// `c_i/c_tau + χ(D_i)/24 - GV_i(t)`, and follows from the previous solution.
///
/// Returns `None` if branch selection fails, any GV correction is invalid, a
/// corrected target is non-positive, or the fixed-point iteration does not
/// converge within `max_gv_iterations`.
#[allow(clippy::too_many_arguments)]
pub fn solve_two_phase_mixed_basis_path_following_with_gv(
    kappa_basis: &Intersection,
    kappa_all: &Intersection,
    basis: &[usize],
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
    chi_divisor: &[I64<Finite>],
    c_tau: CTau,
    gv_invariants: &[(Vec<i32>, Integer)],
    gamma: Option<&[I64<Finite>]>,
    steps: CheckedRange<usize>,
    max_gv_iterations: usize,
    gv_tolerance: F64<NonNeg>,
) -> Option<KkltResult> {
    let debug = kklt_debug_enabled();
    if basis != kklt_basis {
        if debug {
            eprintln!("[KKLT] basis-coordinate GV solver requires basis == kklt_basis");
        }
        return None;
    }
    if max_gv_iterations == 0 {
        if debug {
            eprintln!("[KKLT] max_gv_iterations is zero");
        }
        return None;
    }

    let Some(tau_target) = compute_corrected_target_tau(c_i, chi_divisor, c_tau) else {
        if debug {
            eprintln!("[KKLT] corrected target tau construction failed");
        }
        return None;
    };
    let Some(mut result) = solve_two_phase_mixed_basis_path_following(
        kappa_basis,
        kappa_all,
        basis,
        kklt_basis,
        &tau_target,
        c_i,
        steps,
    ) else {
        if debug {
            eprintln!("[KKLT] zeroth-order mixed-basis path following failed");
        }
        return None;
    };
    if debug {
        eprintln!(
            "[KKLT] zeroth-order solve converged={} rel_err={}",
            result.converged,
            result.relative_error.get()
        );
    }

    for iter in 0..max_gv_iterations {
        let previous_t = result.t.clone();
        let Some(gv_correction) = compute_gv_target_correction(gv_invariants, &previous_t, gamma)
        else {
            if debug {
                eprintln!("[KKLT] GV correction failed at iteration {iter}");
            }
            return None;
        };
        let Some(tau_target) =
            compute_gv_corrected_target_tau(c_i, chi_divisor, c_tau, &gv_correction)
        else {
            if debug {
                eprintln!("[KKLT] GV-corrected target tau failed at iteration {iter}");
            }
            return None;
        };
        let Some(next) = solve_mixed_basis_path_following(
            kappa_basis,
            kappa_all,
            basis,
            kklt_basis,
            &tau_target,
            &previous_t,
            steps,
        ) else {
            if debug {
                eprintln!("[KKLT] corrected mixed-basis path failed at iteration {iter}");
            }
            return None;
        };
        if !next.converged {
            if debug {
                eprintln!(
                    "[KKLT] corrected path did not converge at iteration {iter}: rel_err={}",
                    next.relative_error.get()
                );
            }
            return None;
        }

        let max_relative_step = next
            .t
            .iter()
            .zip(previous_t.iter())
            .map(|(new, old)| (new.get() - old.get()).abs() / (old.get().abs() + 1e-12))
            .fold(0.0f64, f64::max);
        if debug {
            eprintln!(
                "[KKLT] iteration {iter}: max_relative_step={max_relative_step}, rel_err={}",
                next.relative_error.get()
            );
        }
        result = next;
        if max_relative_step <= gv_tolerance.get() {
            return Some(result);
        }
    }

    if debug {
        eprintln!("[KKLT] GV fixed-point iteration did not converge");
    }
    None
}

/// Solve the mixed-basis KKLT system with GV corrections for explicit KKLT
/// divisors.
///
/// Unlike [`solve_two_phase_mixed_basis_path_following_with_gv`], this handles
/// `basis != kklt_basis` by using `curve_basis_matrix` to evaluate each
/// basis-coordinate GV curve on the requested ambient KKLT divisors.
#[allow(clippy::too_many_arguments)]
pub fn solve_two_phase_mixed_basis_path_following_with_divisor_gv(
    kappa_basis: &Intersection,
    kappa_all: &Intersection,
    basis: &[usize],
    kklt_basis: &[usize],
    c_i: &[I64<Pos>],
    chi_divisor: &[I64<Finite>],
    c_tau: CTau,
    gv_invariants: &[(Vec<i32>, Integer)],
    curve_basis_matrix: &[Vec<Integer>],
    gamma: Option<&[I64<Finite>]>,
    steps: CheckedRange<usize>,
    max_gv_iterations: usize,
    gv_tolerance: F64<NonNeg>,
) -> Option<KkltResult> {
    let debug = kklt_debug_enabled();
    if max_gv_iterations == 0 {
        if debug {
            eprintln!("[KKLT] max_gv_iterations is zero");
        }
        return None;
    }

    let Some(tau_target) = compute_corrected_target_tau(c_i, chi_divisor, c_tau) else {
        if debug {
            eprintln!("[KKLT] corrected target tau construction failed");
        }
        return None;
    };
    let Some(mut result) = solve_two_phase_mixed_basis_path_following(
        kappa_basis,
        kappa_all,
        basis,
        kklt_basis,
        &tau_target,
        c_i,
        steps,
    ) else {
        if debug {
            eprintln!("[KKLT] zeroth-order mixed-basis path following failed");
        }
        return None;
    };
    if debug {
        eprintln!(
            "[KKLT] zeroth-order solve converged={} rel_err={}",
            result.converged,
            result.relative_error.get()
        );
    }

    for iter in 0..max_gv_iterations {
        let previous_t = result.t.clone();
        let Some(gv_correction) = compute_gv_target_correction_for_divisors(
            gv_invariants,
            curve_basis_matrix,
            kklt_basis,
            &previous_t,
            gamma,
        ) else {
            if debug {
                eprintln!("[KKLT] divisor GV correction failed at iteration {iter}");
            }
            return None;
        };
        let Some(tau_target) =
            compute_gv_corrected_target_tau(c_i, chi_divisor, c_tau, &gv_correction)
        else {
            if debug {
                eprintln!("[KKLT] GV-corrected target tau failed at iteration {iter}");
            }
            return None;
        };
        let Some(next) = solve_mixed_basis_path_following(
            kappa_basis,
            kappa_all,
            basis,
            kklt_basis,
            &tau_target,
            &previous_t,
            steps,
        ) else {
            if debug {
                eprintln!("[KKLT] corrected mixed-basis path failed at iteration {iter}");
            }
            return None;
        };
        if !next.converged {
            if debug {
                eprintln!(
                    "[KKLT] corrected path did not converge at iteration {iter}: rel_err={}",
                    next.relative_error.get()
                );
            }
            return None;
        }

        let max_relative_step = next
            .t
            .iter()
            .zip(previous_t.iter())
            .map(|(new, old)| (new.get() - old.get()).abs() / (old.get().abs() + 1e-12))
            .fold(0.0f64, f64::max);
        if debug {
            eprintln!(
                "[KKLT] iteration {iter}: max_relative_step={max_relative_step}, rel_err={}",
                next.relative_error.get()
            );
        }
        result = next;
        if max_relative_step <= gv_tolerance.get() {
            return Some(result);
        }
    }

    if debug {
        eprintln!("[KKLT] GV fixed-point iteration did not converge");
    }
    None
}

/// Solve Ax = b using SVD least squares.
fn solve_least_squares(a: &[Vec<F64<Finite>>], b: &[F64<Finite>]) -> Option<Vec<F64<Finite>>> {
    let m = a.len();
    if m == 0 {
        return None;
    }
    let n = a[0].len();
    if n == 0 || b.len() != m {
        return None;
    }
    if a.iter().any(|row| row.len() != n) {
        return None;
    }

    let data: Vec<f64> = a
        .iter()
        .flat_map(|row| row.iter().map(|entry| entry.get()))
        .collect();
    let matrix = nalgebra::DMatrix::<f64>::from_row_slice(m, n, &data);
    let rhs = nalgebra::DVector::<f64>::from_iterator(m, b.iter().map(|entry| entry.get()));
    let solution = matrix.svd(true, true).solve(&rhs, 1e-8).ok()?;

    if solution.iter().any(|value| !value.is_finite()) {
        return None;
    }

    solution
        .iter()
        .map(|value| F64::<Finite>::new(*value))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::range;
    use crate::types::rational::Rational as TypedRational;
    use malachite::{Integer, Rational};

    fn pos_i64(v: i64) -> I64<Pos> {
        I64::<Pos>::new(v).unwrap()
    }

    fn finite_i64(v: i64) -> I64<Finite> {
        I64::<Finite>::new(v)
    }

    fn finite_f64(v: f64) -> F64<Finite> {
        F64::<Finite>::new(v).unwrap()
    }

    #[test]
    fn test_compute_c_tau() {
        // For McAllister 4-214-647:
        // g_s ≈ 0.00911, W₀ ≈ 2.3e-90, expected c_τ ≈ 3.34
        let g_s = f64_pos!(0.00911134);
        let w0 = f64_pos!(2.3e-90);
        let c_tau = compute_c_tau(g_s, w0);
        assert!((c_tau.get() - 3.34).abs() < 0.01, "c_tau = {}", c_tau.get());
    }

    #[test]
    fn test_compute_target_tau() {
        let c_i = vec![pos_i64(6), pos_i64(1), pos_i64(6), pos_i64(1)];
        let c_tau = f64_pos!(3.0);
        let tau = compute_target_tau(&c_i, c_tau);
        assert_eq!(tau.len(), 4);
        assert!((tau[0].get() - 2.0).abs() < 1e-10);
        assert!((tau[1].get() - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_corrected_target_tau() {
        let c_i = vec![pos_i64(6), pos_i64(1)];
        let chi_divisor = vec![finite_i64(24), finite_i64(24)];
        let c_tau = f64_pos!(3.0);

        let tau = compute_corrected_target_tau(&c_i, &chi_divisor, c_tau).unwrap();

        assert_eq!(tau.len(), 2);
        assert!((tau[0].get() - 3.0).abs() < 1e-10);
        assert!((tau[1].get() - (1.0 / 3.0 + 1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_compute_corrected_target_tau_rejects_invalid_targets() {
        let c_i = vec![pos_i64(1)];
        let chi_divisor = vec![finite_i64(-24)];
        let c_tau = f64_pos!(3.0);

        assert!(compute_corrected_target_tau(&c_i, &chi_divisor, c_tau).is_none());
        assert!(compute_corrected_target_tau(&c_i, &[], c_tau).is_none());
    }

    #[test]
    fn test_compute_gv_corrected_target_tau() {
        let c_i = vec![pos_i64(1)];
        let chi_divisor = vec![finite_i64(24)];
        let c_tau = f64_pos!(2.0);
        let gv_correction = vec![finite_f64(0.25)];

        let tau =
            compute_gv_corrected_target_tau(&c_i, &chi_divisor, c_tau, &gv_correction).unwrap();

        assert_eq!(tau.len(), 1);
        assert!((tau[0].get() - 1.25).abs() < 1e-12);
        assert!(compute_gv_corrected_target_tau(&c_i, &chi_divisor, c_tau, &[]).is_none());
    }

    #[test]
    fn test_compute_gv_corrected_target_tau_rejects_non_positive_target() {
        let c_i = vec![pos_i64(1)];
        let chi_divisor = vec![finite_i64(-24)];
        let c_tau = f64_pos!(1.0);
        let gv_correction = vec![finite_f64(1.0)];

        assert!(
            compute_gv_corrected_target_tau(&c_i, &chi_divisor, c_tau, &gv_correction).is_none()
        );
    }

    #[test]
    fn test_compute_gv_target_correction_real_dilog_signs() {
        let t = vec![finite_f64(2.0_f64.ln() / (2.0 * PI))];
        let invariants = vec![(vec![1], Integer::from(1))];

        let positive = compute_gv_target_correction(&invariants, &t, None).unwrap();
        let li2_half = PI * PI / 12.0 - 0.5 * 2.0_f64.ln().powi(2);
        let expected_positive = li2_half / (4.0 * PI * PI);
        assert!((positive[0].get() - expected_positive).abs() < 1e-14);

        let gamma = vec![finite_i64(1)];
        let negative = compute_gv_target_correction(&invariants, &t, Some(&gamma)).unwrap();
        let expected_negative = -0.448_414_206_923_646_2 / (4.0 * PI * PI);
        assert!((negative[0].get() - expected_negative).abs() < 1e-14);
    }

    #[test]
    fn test_compute_gv_target_correction_accepts_real_negative_branch() {
        let t = vec![finite_f64(-2.0_f64.ln() / (2.0 * PI))];
        let invariants = vec![(vec![1], Integer::from(1))];
        let gamma = vec![finite_i64(1)];

        assert!(compute_gv_target_correction(&invariants, &t, None).is_none());

        let correction = compute_gv_target_correction(&invariants, &t, Some(&gamma)).unwrap();
        let expected = -1.436_746_366_883_680_8 / (4.0 * PI * PI);
        assert!((correction[0].get() - expected).abs() < 1e-14);
    }

    #[test]
    fn test_compute_gv_target_correction_rejects_invalid_curve_volume() {
        let t = vec![finite_f64(1.0)];
        let zero_curve = vec![(vec![0], Integer::from(1))];
        let negative_curve = vec![(vec![-1], Integer::from(1))];

        assert!(compute_gv_target_correction(&zero_curve, &t, None).is_none());
        assert!(compute_gv_target_correction(&negative_curve, &t, None).is_none());
        assert!(
            compute_gv_target_correction(&[(vec![1, 0], Integer::from(1))], &t, None).is_none()
        );
    }

    #[test]
    fn test_compute_gv_target_correction_for_divisors_uses_curve_basis_matrix() {
        let t_entry = 2.0_f64.ln() / (10.0 * PI);
        let t = vec![finite_f64(t_entry), finite_f64(t_entry)];
        let invariants = vec![(vec![2, 3], Integer::from(1))];
        let curve_basis = vec![
            vec![Integer::from(1), Integer::from(0), Integer::from(4)],
            vec![Integer::from(0), Integer::from(1), Integer::from(-1)],
        ];
        let kklt_basis = vec![0, 2];

        let correction = compute_gv_target_correction_for_divisors(
            &invariants,
            &curve_basis,
            &kklt_basis,
            &t,
            None,
        )
        .unwrap();

        let li2_half = PI * PI / 12.0 - 0.5 * 2.0_f64.ln().powi(2);
        let prefactor = 1.0 / (4.0 * PI * PI);
        assert!((correction[0].get() - 2.0 * prefactor * li2_half).abs() < 1e-14);
        assert!((correction[1].get() - 5.0 * prefactor * li2_half).abs() < 1e-14);
    }

    #[test]
    fn test_compute_gv_target_correction_for_ambient_curves_matches_basis_mapping() {
        let t_entry = 2.0_f64.ln() / (10.0 * PI);
        let t = vec![finite_f64(t_entry), finite_f64(t_entry)];
        let basis_invariants = vec![(vec![2, 3], Integer::from(1))];
        let ambient_invariants = vec![(vec![2, 3, 5], Integer::from(1))];
        let curve_basis = vec![
            vec![Integer::from(1), Integer::from(0), Integer::from(4)],
            vec![Integer::from(0), Integer::from(1), Integer::from(-1)],
        ];
        let basis = vec![0, 1];
        let kklt_basis = vec![0, 2];

        let from_basis = compute_gv_target_correction_for_divisors(
            &basis_invariants,
            &curve_basis,
            &kklt_basis,
            &t,
            None,
        )
        .unwrap();
        let from_ambient = compute_gv_target_correction_for_ambient_curves(
            &ambient_invariants,
            &basis,
            &kklt_basis,
            &t,
            None,
        )
        .unwrap();

        assert_eq!(from_ambient.len(), from_basis.len());
        for (ambient, basis) in from_ambient.iter().zip(from_basis.iter()) {
            assert!((ambient.get() - basis.get()).abs() < 1e-14);
        }
    }

    #[test]
    fn test_compute_gv_volume_correction_for_ambient_curves() {
        let t = vec![finite_f64(2.0_f64.ln() / (2.0 * PI))];
        let basis = vec![0];
        let invariants = vec![(vec![1], Integer::from(1))];

        let correction =
            compute_gv_volume_correction_for_ambient_curves(&invariants, &basis, &t, None).unwrap();

        let li2_half = PI * PI / 12.0 - 0.5 * 2.0_f64.ln().powi(2);
        let li3_half = 0.537_213_193_608_040_2;
        let expected = (li3_half + 2.0_f64.ln() * li2_half) / (2.0 * (2.0 * PI).powi(3));
        assert!((correction.get() - expected).abs() < 1e-14);
    }

    #[test]
    fn test_compute_divisor_volumes_simple() {
        // κ_000 = 6, t = [2]
        // τ_0 = (1/2) × 6 × 2 × 2 = 12
        let mut kappa = Intersection::new(1);
        kappa.set(
            0,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );

        let t = vec![finite_f64(2.0)];
        let tau = compute_divisor_volumes(&kappa, &t);
        assert!((tau[0].get() - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_jacobian_simple() {
        // κ_000 = 6, t = [2]
        // ∂τ_0/∂t = 6t = 12
        let mut kappa = Intersection::new(1);
        kappa.set(
            0,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );

        let t = vec![finite_f64(2.0)];
        let j = compute_jacobian(&kappa, &t);
        assert!(
            (j[0][0].get() - 12.0).abs() < 1e-10,
            "J[0][0] = {}",
            j[0][0].get()
        );
    }

    #[test]
    fn test_compute_jacobian_distinct_indices_matches_tensor_contraction() {
        let mut kappa = Intersection::new(3);
        kappa.set(
            0,
            1,
            2,
            TypedRational::<Finite>::from_raw(Rational::from(5)),
        );

        let t = vec![finite_f64(2.0), finite_f64(3.0), finite_f64(7.0)];
        let j = compute_jacobian(&kappa, &t);

        assert!((j[0][1].get() - 35.0).abs() < 1e-10);
        assert!((j[0][2].get() - 15.0).abs() < 1e-10);
        assert!((j[1][0].get() - 35.0).abs() < 1e-10);
        assert!((j[1][2].get() - 10.0).abs() < 1e-10);
        assert!((j[2][0].get() - 15.0).abs() < 1e-10);
        assert!((j[2][1].get() - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_jacobian_diagnostics_detect_full_rank_conditioning() {
        let jacobian = vec![
            vec![finite_f64(2.0), finite_f64(0.0)],
            vec![finite_f64(0.0), finite_f64(0.5)],
        ];

        let diagnostics = compute_jacobian_diagnostics(&jacobian).unwrap();

        assert_eq!(diagnostics.rank, 2);
        assert_eq!(diagnostics.max_rank, 2);
        assert!((diagnostics.max_singular_value.get() - 2.0).abs() < 1e-12);
        assert!((diagnostics.min_nonzero_singular_value.get() - 0.5).abs() < 1e-12);
        assert!((diagnostics.condition_number.unwrap().get() - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_jacobian_diagnostics_reports_rank_deficiency() {
        let jacobian = vec![
            vec![finite_f64(1.0), finite_f64(1.0)],
            vec![finite_f64(2.0), finite_f64(2.0)],
        ];

        let diagnostics = compute_jacobian_diagnostics(&jacobian).unwrap();

        assert_eq!(diagnostics.rank, 1);
        assert_eq!(diagnostics.max_rank, 2);
        assert!(diagnostics.condition_number.is_none());
    }

    #[test]
    fn test_compute_kklt_mixed_basis_tau_and_jacobian() {
        // Kähler coordinates use basis divisors [0, 1].
        // KKLT targets are ordered by [0, 2], so divisor 2 must be evaluated
        // from the full tensor κ_{2,b_j,b_k}.
        let mut kappa_basis = Intersection::new(2);
        kappa_basis.set(
            0,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(2)),
        );

        let mut kappa_all = Intersection::new(3);
        kappa_all.set(
            2,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(4)),
        );
        kappa_all.set(
            2,
            0,
            1,
            TypedRational::<Finite>::from_raw(Rational::from(5)),
        );

        let basis = vec![0, 1];
        let kklt_basis = vec![0, 2];
        let t = vec![finite_f64(2.0), finite_f64(3.0)];

        let tau = compute_kklt_divisor_volumes(&kappa_basis, &kappa_all, &basis, &kklt_basis, &t)
            .unwrap();
        assert!((tau[0].get() - 4.0).abs() < 1e-10);
        assert!((tau[1].get() - 38.0).abs() < 1e-10);

        let jacobian =
            compute_kklt_jacobian(&kappa_basis, &kappa_all, &basis, &kklt_basis, &t).unwrap();
        assert!((jacobian[0][0].get() - 4.0).abs() < 1e-10);
        assert!((jacobian[0][1].get() - 0.0).abs() < 1e-10);
        assert!((jacobian[1][0].get() - 23.0).abs() < 1e-10);
        assert!((jacobian[1][1].get() - 10.0).abs() < 1e-10);
    }

    #[test]
    fn test_least_squares_handles_rank_deficient_system() {
        let a = vec![
            vec![finite_f64(1.0), finite_f64(1.0)],
            vec![finite_f64(2.0), finite_f64(2.0)],
        ];
        let b = vec![finite_f64(2.0), finite_f64(4.0)];

        let x = solve_least_squares(&a, &b).unwrap();
        let residual_0 = x[0].get() + x[1].get() - 2.0;
        let residual_1 = 2.0 * x[0].get() + 2.0 * x[1].get() - 4.0;

        assert!(residual_0.abs() < 1e-10, "residual_0 = {residual_0}");
        assert!(residual_1.abs() < 1e-10, "residual_1 = {residual_1}");
    }

    #[test]
    fn test_solve_two_phase_path_following_one_modulus() {
        // κ_000 = 6 gives τ = 3t².
        // Phase 1 target c_i = 3 lands on t = 1; phase 2 target τ = 12 lands on t = 2.
        let mut kappa = Intersection::new(1);
        kappa.set(
            0,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );

        let c_i = vec![pos_i64(3)];
        let tau_target = vec![f64_pos!(12.0)];

        let result =
            solve_two_phase_path_following(&kappa, &tau_target, &c_i, range!(0, 400)).unwrap();

        assert!(result.converged);
        assert!(
            (result.t[0].get() - 2.0).abs() < 1e-5,
            "t = {}",
            result.t[0].get()
        );
        assert!(
            (result.tau[0].get() - 12.0).abs() < 1e-4,
            "tau = {}",
            result.tau[0].get()
        );
    }

    #[test]
    fn test_solve_two_phase_mixed_basis_path_following_one_modulus() {
        // One Kähler coordinate D0, but the KKLT target is for non-basis D1.
        // κ_100 = 6 still gives τ_D1 = 3t².
        let kappa_basis = Intersection::new(1);
        let mut kappa_all = Intersection::new(2);
        kappa_all.set(
            1,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );

        let basis = vec![0];
        let kklt_basis = vec![1];
        let c_i = vec![pos_i64(3)];
        let tau_target = vec![f64_pos!(12.0)];

        let result = solve_two_phase_mixed_basis_path_following(
            &kappa_basis,
            &kappa_all,
            &basis,
            &kklt_basis,
            &tau_target,
            &c_i,
            range!(0, 400),
        )
        .unwrap();

        assert!(result.converged);
        assert!(
            (result.t[0].get() - 2.0).abs() < 1e-5,
            "t = {}",
            result.t[0].get()
        );
        assert!(
            (result.tau[0].get() - 12.0).abs() < 1e-4,
            "tau = {}",
            result.tau[0].get()
        );
    }

    #[test]
    fn test_mixed_basis_path_following_corrects_single_large_step() {
        // With one interpolation step, the predictor alone lands at t=2.5 for
        // τ=3t² and target τ=12. The Newton corrector should bring it back to
        // the true branch point t=2.
        let kappa_basis = Intersection::new(1);
        let mut kappa_all = Intersection::new(2);
        kappa_all.set(
            1,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );

        let basis = vec![0];
        let kklt_basis = vec![1];
        let tau_target = vec![f64_pos!(12.0)];
        let t_init = vec![finite_f64(1.0)];

        let result = solve_mixed_basis_path_following(
            &kappa_basis,
            &kappa_all,
            &basis,
            &kklt_basis,
            &tau_target,
            &t_init,
            range!(0, 1),
        )
        .unwrap();

        assert!(result.converged);
        assert!(
            (result.t[0].get() - 2.0).abs() < 1e-4,
            "t = {}",
            result.t[0].get()
        );
        assert!(
            (result.tau[0].get() - 12.0).abs() < 1e-4,
            "tau = {}",
            result.tau[0].get()
        );
    }

    #[test]
    fn test_branch_candidates_report_positive_volume_solutions() {
        let mut kappa_basis = Intersection::new(1);
        kappa_basis.set(
            0,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );
        let mut kappa_all = Intersection::new(2);
        kappa_all.set(
            1,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );

        let basis = vec![0];
        let kklt_basis = vec![1];
        let tau_target = vec![f64_pos!(12.0)];
        let candidates = vec![vec![finite_f64(1.0)], vec![finite_f64(-1.0)]];

        let search = solve_mixed_basis_path_following_branch_candidates(
            &kappa_basis,
            &kappa_all,
            &basis,
            &kklt_basis,
            &tau_target,
            &candidates,
            range!(0, 4),
        );

        assert_eq!(search.attempted, 2);
        assert_eq!(search.solved, 2);
        assert_eq!(search.non_converged, 0);
        assert_eq!(search.non_positive_volume, 1);
        assert_eq!(search.positive_volume.len(), 1);
        assert_eq!(search.positive_volume[0].init_index, 0);
        assert!((search.positive_volume[0].result.t[0].get() - 2.0).abs() < 1e-4);
        assert!((search.positive_volume[0].classical_volume.get() - 8.0).abs() < 1e-4);
        assert_eq!(search.positive_volume[0].jacobian_diagnostics.rank, 1);
        assert_eq!(search.positive_volume[0].jacobian_diagnostics.max_rank, 1);
        assert!(
            search.positive_volume[0]
                .jacobian_diagnostics
                .condition_number
                .is_some()
        );
    }

    #[test]
    fn test_branch_initialization_generation_is_seeded_and_scaled() {
        let mut kappa_basis = Intersection::new(1);
        kappa_basis.set(
            0,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );
        let mut kappa_all = Intersection::new(2);
        kappa_all.set(
            1,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );

        let basis = vec![0];
        let kklt_basis = vec![1];
        let tau_target = vec![f64_pos!(12.0)];
        let scaled_explicit = scale_mixed_basis_kklt_branch_initialization_to_target(
            &kappa_basis,
            &kappa_all,
            &basis,
            &kklt_basis,
            &tau_target,
            &[finite_f64(1.0)],
        )
        .unwrap();
        assert!((scaled_explicit[0].get() - 2.0).abs() < 1e-12);

        let first = generate_scaled_kklt_branch_initializations(
            &kappa_basis,
            &kappa_all,
            &basis,
            &kklt_basis,
            &tau_target,
            8,
            123,
        )
        .unwrap();
        let second = generate_scaled_kklt_branch_initializations(
            &kappa_basis,
            &kappa_all,
            &basis,
            &kklt_basis,
            &tau_target,
            8,
            123,
        )
        .unwrap();

        assert_eq!(first, second);
        assert_eq!(first.len(), 8);
        assert!(first.iter().all(|candidate| candidate.len() == 1));
        for candidate in first {
            let tau = compute_kklt_divisor_volumes(
                &kappa_basis,
                &kappa_all,
                &basis,
                &kklt_basis,
                &candidate,
            )
            .unwrap();
            assert!(
                (tau[0].get().abs() - 12.0).abs() < 1e-8,
                "tau = {}",
                tau[0].get()
            );
        }
    }

    #[test]
    fn test_solve_two_phase_mixed_basis_path_following_with_gv_one_modulus() {
        let mut kappa_basis = Intersection::new(1);
        kappa_basis.set(
            0,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );
        let kappa_all = kappa_basis.clone();
        let basis = vec![0];
        let kklt_basis = vec![0];
        let c_i = vec![pos_i64(1)];
        let chi_divisor = vec![finite_i64(0)];
        let c_tau = f64_pos!(1.0);
        let invariants = vec![(vec![1], Integer::from(100))];

        let result = solve_two_phase_mixed_basis_path_following_with_gv(
            &kappa_basis,
            &kappa_all,
            &basis,
            &kklt_basis,
            &c_i,
            &chi_divisor,
            c_tau,
            &invariants,
            None,
            range!(0, 400),
            20,
            F64::<NonNeg>::new(1e-10).unwrap(),
        )
        .unwrap();

        assert!(result.converged);
        let gv_correction = compute_gv_target_correction(&invariants, &result.t, None).unwrap();
        let tau_target =
            compute_gv_corrected_target_tau(&c_i, &chi_divisor, c_tau, &gv_correction).unwrap();
        assert!(
            (result.tau[0].get() - tau_target[0].get()).abs() < 1e-6,
            "tau={}, target={}",
            result.tau[0].get(),
            tau_target[0].get()
        );
    }

    #[test]
    fn test_solve_two_phase_mixed_basis_path_following_with_divisor_gv_one_modulus() {
        let kappa_basis = Intersection::new(1);
        let mut kappa_all = Intersection::new(2);
        kappa_all.set(
            1,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(6)),
        );
        let basis = vec![0];
        let kklt_basis = vec![1];
        let c_i = vec![pos_i64(1)];
        let chi_divisor = vec![finite_i64(0)];
        let c_tau = f64_pos!(1.0);
        let invariants = vec![(vec![1], Integer::from(1))];
        let curve_basis = vec![vec![Integer::from(1), Integer::from(1)]];

        let result = solve_two_phase_mixed_basis_path_following_with_divisor_gv(
            &kappa_basis,
            &kappa_all,
            &basis,
            &kklt_basis,
            &c_i,
            &chi_divisor,
            c_tau,
            &invariants,
            &curve_basis,
            None,
            range!(0, 400),
            20,
            F64::<NonNeg>::new(1e-10).unwrap(),
        )
        .unwrap();

        assert!(result.converged);
        let gv_correction = compute_gv_target_correction_for_divisors(
            &invariants,
            &curve_basis,
            &kklt_basis,
            &result.t,
            None,
        )
        .unwrap();
        let tau_target =
            compute_gv_corrected_target_tau(&c_i, &chi_divisor, c_tau, &gv_correction).unwrap();
        assert!(
            (result.tau[0].get() - tau_target[0].get()).abs() < 1e-6,
            "tau={}, target={}",
            result.tau[0].get(),
            tau_target[0].get()
        );
    }
}
