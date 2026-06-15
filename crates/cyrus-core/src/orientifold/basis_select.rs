//! KKLT divisor-basis selection from first principles (arXiv:2107.09064 SS2).
//!
//! The superpotential is parameterized by a basis of `h^{1,1}` linearly
//! independent RIGID prime toric divisors. The paper's admissibility criterion
//! is that, with the basis divisor volumes set to unity, the remaining four
//! prime toric divisors have strictly positive volume - equivalently the
//! constant covector `tau_* = (1,...,1)` lies in the dual effective cone
//! `E(X)^0`. Because divisor volume is linear in the divisor class and the
//! volume vector obeys the same coordinate relations as the classes, this is a
//! single `dim x dim` solve (no KKLT iteration); see
//! [`complement_remaining_volumes`].

use crate::error::Result;
use crate::{Point, Polytope};

use super::classify_divisor_rigidity;

/// Integer determinant of a small square matrix via Laplace expansion.
/// Exact (no overflow for the small dim-x-dim coordinate matrices here).
fn det_i128(m: &[Vec<i128>]) -> i128 {
    let n = m.len();
    match n {
        0 => 1,
        1 => m[0][0],
        2 => m[0][0] * m[1][1] - m[0][1] * m[1][0],
        _ => {
            let mut det = 0i128;
            for (col, &top) in m[0].iter().enumerate() {
                if top == 0 {
                    continue;
                }
                let minor: Vec<Vec<i128>> = m[1..]
                    .iter()
                    .map(|row| {
                        row.iter()
                            .enumerate()
                            .filter(|&(c, _)| c != col)
                            .map(|(_, &v)| v)
                            .collect()
                    })
                    .collect();
                let sign = if col % 2 == 0 { 1 } else { -1 };
                det += sign * top * det_i128(&minor);
            }
            det
        }
    }
}

/// The volumes of the `complement` (non-basis) prime toric divisors when the
/// basis divisor volumes are all set to unity (`tau_* = (1,...,1)`).
///
/// Divisor volume is linear in the divisor class, so the full volume vector
/// obeys the same `dim` coordinate relations the classes do
/// (`sum_I (v_I)_k tau_I = 0`). With the basis volumes pinned to 1, the
/// `complement` volumes are determined: `tau_C = 1 - V_C^{-1} S`, where `V_C`
/// is the `dim x dim` matrix of the complement points' coordinates and `S` is
/// the sum of all non-origin ray coordinates. Returned exactly as
/// `(det V_C, numerators)` with `tau_C[i] = numerators[i] / det` (Cramer's
/// rule). `det = 0` means the complement is degenerate (the basis is not even
/// linearly independent); `det != 0` means the basis is a genuine `Q`-basis of
/// the divisor classes (it need NOT be unimodular - McAllister's 4-214 basis
/// has `det = 2`).
///
/// This is exactly the paper's (arXiv:2107.09064 SS2) effective-cone test:
/// the basis is admissible iff `det != 0` and every `tau_C[i] > 0` (so
/// `tau_* ∈ E(X)^0`).
#[must_use]
pub fn complement_remaining_volumes(points: &[Point], complement: &[usize]) -> (i128, Vec<i128>) {
    let dim = points.first().map_or(0, Point::dim);
    if complement.len() != dim {
        return (0, Vec::new());
    }
    // S = sum of all non-origin ray coordinates.
    let mut s = vec![0i128; dim];
    for p in points {
        if p.coords().iter().all(|&c| c == 0) {
            continue;
        }
        for (k, &c) in p.coords().iter().enumerate() {
            s[k] += i128::from(c);
        }
    }
    // V_C[k][c] = (v_{complement[c]})_k.
    let vc: Vec<Vec<i128>> = (0..dim)
        .map(|k| {
            complement
                .iter()
                .map(|&idx| i128::from(points[idx].coords()[k]))
                .collect()
        })
        .collect();
    let det = det_i128(&vc);
    if det == 0 {
        return (0, Vec::new());
    }
    // x_i = det(V_C | col i -> S) / det ; tau_C[i] = 1 - x_i = (det - det_i)/det.
    let numerators = (0..dim)
        .map(|i| {
            let mut a = vc.clone();
            for (k, row) in a.iter_mut().enumerate() {
                row[i] = s[k];
            }
            det - det_i128(&a)
        })
        .collect();
    (det, numerators)
}

/// Whether a candidate `complement` yields an admissible KKLT basis.
///
/// The basis (all non-origin prime divisors minus `complement`) must be
/// linearly independent and every remaining divisor volume strictly positive
/// at `tau_* = (1,...,1)`.
#[must_use]
pub fn complement_is_admissible(points: &[Point], complement: &[usize]) -> bool {
    let (det, numerators) = complement_remaining_volumes(points, complement);
    det != 0 && numerators.iter().all(|&n| n != 0 && (n > 0) == (det > 0))
}

/// Advance `idx` to the next k-combination of `0..n` in lexicographic order.
/// Returns false when the last combination has been passed.
fn next_combination(idx: &mut [usize], n: usize) -> bool {
    let k = idx.len();
    if k == 0 {
        return false;
    }
    let mut i = k;
    while i > 0 {
        i -= 1;
        if idx[i] != i + n - k {
            idx[i] += 1;
            for j in (i + 1)..k {
                idx[j] = idx[j - 1] + 1;
            }
            return true;
        }
    }
    false
}

/// Enumerate admissible KKLT divisor bases, per arXiv:2107.09064 SS2, ordered
/// from most to least interior (descending minimum remaining-divisor volume at
/// `tau_* = (1,...,1)`).
///
/// Each basis is a set of `h^{1,1}` linearly independent RIGID prime toric
/// divisors such that, with their volumes set to unity, the remaining prime
/// toric divisors have strictly positive volume (`tau_* ∈ E(X)^0`). The
/// non-rigid prime divisors are necessarily excluded; the remaining
/// `dim - (#non-rigid)` excluded divisors are found by brute force over the
/// rigid divisors. At most `limit` bases are returned (0 = all).
///
/// # Errors
/// Returns an error if rigidity classification fails or if there are more
/// non-rigid prime divisors than the codimension `dim`.
pub fn admissible_kklt_bases(
    polytope: &Polytope,
    points: &[Point],
    limit: usize,
) -> Result<Vec<Vec<usize>>> {
    let dim = points.first().map_or(0, Point::dim);
    let (classes, _components) = classify_divisor_rigidity(polytope, points)?;
    let mut rigid: Vec<usize> = classes
        .iter()
        .filter(|(_, c)| c.is_rigid())
        .map(|(idx, _)| *idx)
        .collect();
    rigid.sort_unstable();
    let non_rigid: Vec<usize> = classes
        .iter()
        .filter(|(_, c)| !c.is_rigid())
        .map(|(idx, _)| *idx)
        .collect();
    if non_rigid.len() > dim {
        return Err(crate::error::Error::InvalidInput(format!(
            "{} non-rigid prime toric divisors exceed codimension {dim}; no all-rigid KKLT basis exists",
            non_rigid.len()
        )));
    }
    let drop_count = dim - non_rigid.len();

    // Score each admissible drop-set by its minimum remaining divisor volume.
    let mut scored: Vec<(f64, Vec<usize>)> = Vec::new();
    let mut score = |dropped: &[usize]| {
        let mut complement = non_rigid.clone();
        complement.extend_from_slice(dropped);
        let (det, numerators) = complement_remaining_volumes(points, &complement);
        if det == 0 || !numerators.iter().all(|&n| n != 0 && (n > 0) == (det > 0)) {
            return;
        }
        #[allow(clippy::cast_precision_loss)]
        let min_vol = numerators
            .iter()
            .map(|&n| n as f64 / det as f64)
            .fold(f64::INFINITY, f64::min);
        scored.push((min_vol, dropped.to_vec()));
    };

    if drop_count == 0 {
        score(&[]);
    } else {
        let mut combo: Vec<usize> = (0..drop_count).collect();
        loop {
            let dropped: Vec<usize> = combo.iter().map(|&i| rigid[i]).collect();
            score(&dropped);
            if !next_combination(&mut combo, rigid.len()) {
                break;
            }
        }
    }
    scored.sort_by(|a, b| b.0.total_cmp(&a.0));
    if limit > 0 {
        scored.truncate(limit);
    }

    let drop_to_basis = |dropped: &[usize]| -> Vec<usize> {
        let drop_set: std::collections::HashSet<usize> = dropped.iter().copied().collect();
        rigid
            .iter()
            .copied()
            .filter(|idx| !drop_set.contains(idx))
            .collect()
    };
    Ok(scored.iter().map(|(_, d)| drop_to_basis(d)).collect())
}

/// Select a single KKLT divisor basis from first principles: the most interior
/// admissible basis (see [`admissible_kklt_bases`]).
///
/// # Errors
/// Returns an error if no admissible basis exists.
pub fn select_kklt_basis(polytope: &Polytope, points: &[Point]) -> Result<Vec<usize>> {
    admissible_kklt_bases(polytope, points, 1)?
        .into_iter()
        .next()
        .ok_or_else(|| {
            crate::error::Error::InvalidInput(
                "no admissible KKLT basis: tau_* leaves the dual effective cone".into(),
            )
        })
}
