//! Bias-free exact-isotropic K enumeration over an integer box.
//!
//! The isotropy condition `K^T adj K == 0` is quadratic in K, so for each
//! choice of the first `dim-1` components it reduces to a quadratic in the
//! last - one solve replaces a full sweep over that component, a
//! `(2*k_box+1)`-fold speedup returning the IDENTICAL set as testing every
//! K. The form's constant and linear terms are maintained incrementally down
//! the recursion (see [`extend_partial`]), so each leaf reads its quadratic
//! coefficients in O(1) instead of recomputing an O(dim^2) sum. The
//! set-equality property is proven against brute force in the pfv tests; the
//! seed scan uses it.

use crate::pfv::is_isotropic_int;

/// Floor integer square root of a non-negative i128 (0 for n <= 0).
const fn isqrt_i128(n: i128) -> i128 {
    if n <= 0 { 0 } else { n.isqrt() }
}

/// Extend the running partial form by fixing `k[depth] = v`, in O(dim).
///
/// The recursion maintains two pieces of state, valid at entry to `depth`
/// (the first `depth` components already fixed):
///   * `c_partial = sum_{a,b < depth} adj[a][b] k_a k_b` - the form's value
///     among the fixed components (the leaf's constant term C), and
///   * `grad[j] = sum_{a < depth} (adj[a][j] + adj[j][a]) k_a` for `j >=
///     depth` - the linear coefficient that component `j` will carry (the
///     leaf's `B = grad[last]`).
///
/// Fixing `k[depth] = v` updates both for the deeper level: `grad[j]` gains
/// `(adj[depth][j] + adj[j][depth]) v` for `j > depth`, and the new constant
/// term is `c_partial + grad[depth] v + adj[depth][depth] v^2`. This is the
/// same sum the old per-leaf double loop computed, just reassociated and
/// shared across the subtree - so the enumerated set is IDENTICAL while the
/// per-leaf cost drops from O(dim^2) to O(1).
///
/// Returns the new constant term, or `None` (leaving `grad` UNCHANGED) on
/// i128 overflow, so the caller falls back to exact per-value testing for
/// that subtree. The `grad` updates are validated before any are committed,
/// so an overflow never leaves the shared array partially mutated.
fn extend_partial(
    adj: &[Vec<i128>],
    depth: usize,
    v: i64,
    c_partial: i128,
    grad: &mut [i128],
) -> Option<i128> {
    let dim = grad.len();
    let vi = i128::from(v);
    let v2 = vi.checked_mul(vi)?;
    let nc = c_partial
        .checked_add(grad[depth].checked_mul(vi)?)?
        .checked_add(adj[depth][depth].checked_mul(v2)?)?;
    for j in (depth + 1)..dim {
        let coeff = adj[depth][j].checked_add(adj[j][depth])?;
        grad[j].checked_add(coeff.checked_mul(vi)?)?;
    }
    for j in (depth + 1)..dim {
        grad[j] += (adj[depth][j] + adj[j][depth]) * vi;
    }
    Some(nc)
}

/// Reverse [`extend_partial`]'s gradient update (subtract the same deltas),
/// restoring `grad` for backtracking to the next value at this depth.
fn revert_partial(adj: &[Vec<i128>], depth: usize, v: i64, grad: &mut [i128]) {
    let dim = grad.len();
    let vi = i128::from(v);
    for j in (depth + 1)..dim {
        grad[j] -= (adj[depth][j] + adj[j][depth]) * vi;
    }
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

/// `c_partial` and `grad` carry the incrementally-maintained form (see
/// [`extend_partial`]); `valid` is false once an accumulation overflowed,
/// in which case the leaf tests exactly instead of trusting the state.
#[allow(clippy::too_many_arguments)] // recursion threads the incremental state
fn enumerate_rec(
    adj: &[Vec<i128>],
    k_box: i64,
    depth: usize,
    k: &mut [i64],
    c_partial: i128,
    valid: bool,
    grad: &mut [i128],
    out: &mut Vec<Vec<i64>>,
    max_results: usize,
) -> bool {
    if out.len() >= max_results {
        return true;
    }
    let dim = k.len();
    let last = dim - 1;
    if depth == last {
        // Leaf: the quadratic in k[last] is A x^2 + B x + C with A, B, C read
        // straight off the maintained state - no O(dim^2) recompute. The solve
        // may still overflow in the discriminant (b^2); on that, or on an
        // earlier accumulation overflow (`!valid`), fall back to exact
        // per-value testing, which recomputes the full form and is correct
        // regardless.
        if valid
            && let Some(roots) = solve_quadratic_i128(adj[last][last], grad[last], c_partial, k_box)
        {
            for x in roots {
                k[last] = x;
                if k.iter().any(|&v| v != 0) {
                    out.push(k.to_vec());
                    if out.len() >= max_results {
                        return true;
                    }
                }
            }
            return false;
        }
        for x in -k_box..=k_box {
            k[last] = x;
            if k.iter().any(|&v| v != 0) && is_isotropic_int(adj, k) {
                out.push(k.to_vec());
                if out.len() >= max_results {
                    return true;
                }
            }
        }
        return false;
    }
    for v in -k_box..=k_box {
        k[depth] = v;
        // Extend the partial form for the deeper level; on overflow drop to
        // the per-value path (`valid = false`) for this subtree.
        let (next_c, next_valid, applied) = if valid {
            match extend_partial(adj, depth, v, c_partial, grad) {
                Some(nc) => (nc, true, true),
                None => (c_partial, false, false),
            }
        } else {
            (c_partial, false, false)
        };
        let stop = enumerate_rec(
            adj,
            k_box,
            depth + 1,
            k,
            next_c,
            next_valid,
            grad,
            out,
            max_results,
        );
        if applied {
            revert_partial(adj, depth, v, grad);
        }
        if stop {
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
    let mut grad = vec![0i128; dim];
    enumerate_rec(
        adj,
        k_box,
        0,
        &mut k,
        0,
        true,
        &mut grad,
        &mut out,
        max_results,
    );
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extend_partial_reproduces_leaf_coeffs() {
        // Fixing k0=1 then k1=-1 on a symmetric 3x3 form must yield the same
        // leaf coefficients the old per-leaf double loop produced: A=6 (=
        // adj[2][2]), B=-4 (= grad[2]), C=1 (= c_partial). revert_partial must
        // restore grad exactly for backtracking.
        let adj = vec![vec![1i128, 2, 3], vec![2, 4, 5], vec![3, 5, 6]];
        let mut grad = vec![0i128; 3];
        let c0 = extend_partial(&adj, 0, 1, 0, &mut grad).expect("no overflow");
        assert_eq!(c0, 1);
        assert_eq!(grad, vec![0, 4, 6]);
        let c1 = extend_partial(&adj, 1, -1, c0, &mut grad).expect("no overflow");
        assert_eq!(c1, 1);
        assert_eq!(grad, vec![0, 4, -4]);
        // Leaf coefficients (A, B, C):
        assert_eq!((adj[2][2], grad[2], c1), (6, -4, 1));
        // Backtrack: revert in reverse order restores the zero gradient.
        revert_partial(&adj, 1, -1, &mut grad);
        assert_eq!(grad, vec![0, 4, 6]);
        revert_partial(&adj, 0, 1, &mut grad);
        assert_eq!(grad, vec![0, 0, 0]);
    }

    #[test]
    fn solve_quadratic_covers_all_branches() {
        // a==0, b==0, c==0: degenerate, every value in the box.
        assert_eq!(
            solve_quadratic_i128(0, 0, 0, 2),
            Some(vec![-2, -1, 0, 1, 2])
        );
        // a==0, b!=0, exact integer root: 3x - 6 = 0 -> 2.
        assert_eq!(solve_quadratic_i128(0, 3, -6, 5), Some(vec![2]));
        // a==0, b!=0, non-integer root: 2x - 3 = 0 -> none.
        assert_eq!(solve_quadratic_i128(0, 2, -3, 5), Some(vec![]));
        // a!=0, perfect-square discriminant: x^2 - 1 = 0 -> +-1.
        let mut r = solve_quadratic_i128(1, 0, -1, 5).expect("some");
        r.sort_unstable();
        assert_eq!(r, vec![-1, 1]);
        // a!=0, negative discriminant: x^2 + 1 = 0 -> none.
        assert_eq!(solve_quadratic_i128(1, 0, 1, 5), Some(vec![]));
        // a!=0, non-perfect-square discriminant: x^2 - 2 = 0 -> none.
        assert_eq!(solve_quadratic_i128(1, 0, -2, 5), Some(vec![]));
        // a!=0, real roots outside the box: x^2 - 100 = 0 -> +-10 excluded.
        assert_eq!(solve_quadratic_i128(1, 0, -100, 5), Some(vec![]));
    }

    #[test]
    fn solve_quadratic_signals_overflow_with_none() {
        // b^2 overflows i128 before the discriminant subtraction, so the solve
        // must return None (the caller then falls back to per-value testing)
        // rather than silently producing a wrong root set.
        let f: i128 = 10_000_000_000_000_000_000; // 1e19; (2f)^2 = 4e38 > i128::MAX
        assert_eq!(solve_quadratic_i128(f, 2 * f, f, 4), None);
    }

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
