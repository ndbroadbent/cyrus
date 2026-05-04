//! Projection between secondary-fan heights and Kähler coordinates.
//!
//! This ports the projection step from CYTools `utils.py`
//! `project_heights_to_kahler` / `heights_to_kahler`.  The caller supplies the
//! divisor basis and the effective-cone rays for the non-basis prime divisors;
//! computing those rays is a separate toric-geometry step.

use std::collections::HashSet;

use crate::types::f64::F64;
use crate::types::tags::Finite;

/// Project an `(h11 + 5)` height vector to the corresponding ambient divisor
/// coordinates with origin coordinate fixed to zero.
///
/// `basis_non_origin` contains indices into `heights[1..]`.  `prime_divisors`
/// must contain one row for every non-origin, non-basis divisor, in ascending
/// non-origin index order.  Each row gives that divisor's coefficients in the
/// supplied basis, matching CYTools' effective-cone ray convention.
#[must_use]
pub fn project_heights_to_kahler(
    heights: &[F64<Finite>],
    basis_non_origin: &[usize],
    prime_divisors: &[Vec<F64<Finite>>],
) -> Option<Vec<F64<Finite>>> {
    let non_origin_count = heights.len().checked_sub(1)?;
    if basis_non_origin.iter().any(|&idx| idx >= non_origin_count) {
        return None;
    }

    let basis_set: HashSet<usize> = basis_non_origin.iter().copied().collect();
    if basis_set.len() != basis_non_origin.len() {
        return None;
    }

    let extra_divisors: Vec<usize> = (0..non_origin_count)
        .filter(|idx| !basis_set.contains(idx))
        .collect();
    if prime_divisors.len() != extra_divisors.len()
        || prime_divisors
            .iter()
            .any(|ray| ray.len() != basis_non_origin.len())
    {
        return None;
    }

    let origin_height = heights[0];
    let mut projected: Vec<F64<Finite>> = heights[1..]
        .iter()
        .map(|height| *height - origin_height)
        .collect();

    for (&prime_idx, prime_ray) in extra_divisors.iter().zip(prime_divisors.iter()) {
        let prime_height = projected[prime_idx];
        for (&basis_idx, &coefficient) in basis_non_origin.iter().zip(prime_ray.iter()) {
            projected[basis_idx] = projected[basis_idx] + prime_height * coefficient;
        }
        projected[prime_idx] = projected[prime_idx] - prime_height;
    }

    let mut with_origin = Vec::with_capacity(heights.len());
    with_origin.push(F64::<Finite>::ZERO);
    with_origin.extend(projected);
    Some(with_origin)
}

/// Project heights directly to the chosen Kähler basis coordinates.
#[must_use]
pub fn heights_to_kahler(
    heights: &[F64<Finite>],
    basis_non_origin: &[usize],
    prime_divisors: &[Vec<F64<Finite>>],
) -> Option<Vec<F64<Finite>>> {
    let projected = project_heights_to_kahler(heights, basis_non_origin, prime_divisors)?;
    basis_non_origin
        .iter()
        .map(|&idx| projected.get(idx + 1).copied())
        .collect()
}

/// Embed Kähler basis coordinates into a height vector with non-basis heights
/// set to zero, matching CYTools `kahler_to_heights`.
#[must_use]
pub fn kahler_to_heights(
    kahler: &[F64<Finite>],
    basis_non_origin: &[usize],
    non_origin_count: usize,
) -> Option<Vec<F64<Finite>>> {
    if kahler.len() != basis_non_origin.len()
        || basis_non_origin.iter().any(|&idx| idx >= non_origin_count)
    {
        return None;
    }
    let basis_set: HashSet<usize> = basis_non_origin.iter().copied().collect();
    if basis_set.len() != basis_non_origin.len() {
        return None;
    }

    let mut heights = vec![F64::<Finite>::ZERO; non_origin_count + 1];
    for (&basis_idx, &value) in basis_non_origin.iter().zip(kahler.iter()) {
        heights[basis_idx + 1] = value;
    }
    Some(heights)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn finite(value: f64) -> F64<Finite> {
        F64::<Finite>::new(value).expect("finite test value")
    }

    #[test]
    fn project_heights_matches_cytools_linear_relation_step() {
        let heights = vec![finite(10.0), finite(12.0), finite(15.0), finite(20.0)];
        let basis = vec![0, 2];
        let prime_divisors = vec![vec![finite(2.0), finite(-3.0)]];

        let projected = project_heights_to_kahler(&heights, &basis, &prime_divisors).unwrap();

        assert!((projected[0].get() - 0.0).abs() < 1e-12);
        assert!((projected[1].get() - 12.0).abs() < 1e-12);
        assert!((projected[2].get() - 0.0).abs() < 1e-12);
        assert!((projected[3].get() + 5.0).abs() < 1e-12);

        let kahler = heights_to_kahler(&heights, &basis, &prime_divisors).unwrap();
        assert_eq!(kahler.len(), 2);
        assert!((kahler[0].get() - 12.0).abs() < 1e-12);
        assert!((kahler[1].get() + 5.0).abs() < 1e-12);
    }

    #[test]
    fn kahler_to_heights_sets_only_basis_entries() {
        let kahler = vec![finite(12.0), finite(-5.0)];
        let heights = kahler_to_heights(&kahler, &[0, 2], 3).unwrap();

        assert_eq!(heights.len(), 4);
        assert!((heights[0].get() - 0.0).abs() < 1e-12);
        assert!((heights[1].get() - 12.0).abs() < 1e-12);
        assert!((heights[2].get() - 0.0).abs() < 1e-12);
        assert!((heights[3].get() + 5.0).abs() < 1e-12);
    }

    #[test]
    fn height_projection_rejects_mismatched_effective_cone_rows() {
        let heights = vec![finite(0.0), finite(1.0), finite(2.0)];
        assert!(project_heights_to_kahler(&heights, &[0], &[]).is_none());
        assert!(
            project_heights_to_kahler(&heights, &[0], &[vec![finite(1.0), finite(2.0)]]).is_none()
        );
        assert!(kahler_to_heights(&[finite(1.0)], &[0, 0], 2).is_none());
    }
}
