//! Bias-free exact-isotropic K enumeration over an integer box.
//!
//! The isotropy condition `K^T adj K == 0` is quadratic in K, so for each
//! choice of the first `dim-1` components it reduces to a quadratic in the
//! last - one solve replaces a full sweep over that component, a
//! `(2*k_box+1)`-fold speedup returning the IDENTICAL set as testing every
//! K. The set-equality property is proven against brute force in the pfv
//! tests; the seed scan uses it.

use crate::pfv::is_isotropic_int;

/// Floor integer square root of a non-negative i128 (0 for n <= 0).
const fn isqrt_i128(n: i128) -> i128 {
    if n <= 0 { 0 } else { n.isqrt() }
}

/// Quadratic coefficients (A, B, C) of the isotropy form in the LAST
/// component, given the first `dim-1` components fixed in `k`:
/// `Q = A*k_last^2 + B*k_last + C`. Returns `None` on i128 overflow (the
/// caller then falls back to exact per-value testing).
fn quad_coeffs_i128(adj: &[Vec<i128>], k: &[i64]) -> Option<(i128, i128, i128)> {
    let last = k.len() - 1;
    let a = adj[last][last];
    let mut b: i128 = 0;
    for (i, &ki) in k.iter().enumerate().take(last) {
        b = b.checked_add(adj[i][last].checked_mul(i128::from(ki))?)?;
    }
    b = b.checked_mul(2)?;
    let mut c: i128 = 0;
    for (i, &ki) in k.iter().enumerate().take(last) {
        for (j, &kj) in k.iter().enumerate().take(last) {
            let term = adj[i][j]
                .checked_mul(i128::from(ki))?
                .checked_mul(i128::from(kj))?;
            c = c.checked_add(term)?;
        }
    }
    Some((a, b, c))
}

/// Integer roots of `A x^2 + B x + C = 0` within `[-k_box, k_box]`.
/// `None` signals i128 overflow in the discriminant (caller falls back).
/// `Some(roots)` may be empty (no integer roots) or, for the identically
/// zero form, every value in the box.
#[allow(clippy::many_single_char_names)] // A, B, C are the quadratic's coeffs
fn solve_quadratic_i128(a: i128, b: i128, c: i128, k_box: i64) -> Option<Vec<i64>> {
    let kb = i128::from(k_box);
    let mut roots = Vec::new();
    if a == 0 {
        if b == 0 {
            if c == 0 {
                // Degenerate: the form vanishes for every last component.
                roots.extend(-k_box..=k_box);
            }
        } else if (-c) % b == 0 {
            let x = (-c) / b;
            if x.abs() <= kb {
                roots.push(x as i64);
            }
        }
        return Some(roots);
    }
    let disc = b
        .checked_mul(b)?
        .checked_sub(4i128.checked_mul(a)?.checked_mul(c)?)?;
    if disc < 0 {
        return Some(roots);
    }
    let s = isqrt_i128(disc);
    if s.checked_mul(s)? != disc {
        return Some(roots); // discriminant is not a perfect square
    }
    let two_a = 2i128.checked_mul(a)?;
    for signed in [s, -s] {
        let num = (-b).checked_add(signed)?;
        if num % two_a == 0 {
            let x = num / two_a;
            if x.abs() <= kb {
                let xi = x as i64;
                if !roots.contains(&xi) {
                    roots.push(xi);
                }
            }
        }
    }
    Some(roots)
}

fn enumerate_rec(
    adj: &[Vec<i128>],
    k_box: i64,
    depth: usize,
    k: &mut [i64],
    out: &mut Vec<Vec<i64>>,
    max_results: usize,
) -> bool {
    if out.len() >= max_results {
        return true;
    }
    let dim = k.len();
    if depth + 1 == dim {
        // Last component: solve the quadratic, falling back to exact
        // per-value testing on i128 overflow.
        let last = dim - 1;
        let solved =
            quad_coeffs_i128(adj, k).and_then(|(a, b, c)| solve_quadratic_i128(a, b, c, k_box));
        match solved {
            Some(roots) => {
                for x in roots {
                    k[last] = x;
                    if k.iter().any(|&v| v != 0) {
                        out.push(k.to_vec());
                        if out.len() >= max_results {
                            return true;
                        }
                    }
                }
            }
            None => {
                for x in -k_box..=k_box {
                    k[last] = x;
                    if k.iter().any(|&v| v != 0) && is_isotropic_int(adj, k) {
                        out.push(k.to_vec());
                        if out.len() >= max_results {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }
    for v in -k_box..=k_box {
        k[depth] = v;
        if enumerate_rec(adj, k_box, depth + 1, k, out, max_results) {
            return true;
        }
    }
    false
}

/// Enumerate the exact-isotropic K's in a box, bias-free.
///
/// Returns every nonzero integer K in `[-k_box, k_box]^dim` satisfying
/// `K^T adj K == 0`, by solving the quadratic in the last component for
/// each choice of the first `dim-1` components.
///
/// This returns the IDENTICAL set to testing every K with
/// [`is_isotropic_int`], but ~`(2*k_box+1)`-fold faster (one quadratic
/// solve replaces a full sweep over the last component). It is therefore a
/// bias-free speedup of the seed scan: same seeds, fewer operations.
/// `max_results` bounds memory for the degenerate (identically-zero) form.
#[must_use]
pub fn enumerate_isotropic_in_box(
    adj: &[Vec<i128>],
    dim: usize,
    k_box: i64,
    max_results: usize,
) -> Vec<Vec<i64>> {
    let mut out = Vec::new();
    if dim == 0 {
        return out;
    }
    let mut k = vec![0i64; dim];
    enumerate_rec(adj, k_box, 0, &mut k, &mut out, max_results);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn isqrt_is_floor_sqrt() {
        for n in [
            0i128, 1, 2, 3, 4, 8, 9, 15, 16, 10_000, 10_001, 999_999, 1_000_000,
        ] {
            let r = isqrt_i128(n);
            assert!(r * r <= n && (r + 1) * (r + 1) > n, "isqrt({n}) = {r}");
        }
    }
}
