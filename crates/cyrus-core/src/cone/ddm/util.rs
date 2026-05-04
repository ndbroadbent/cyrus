use malachite::{Integer, Rational};

use crate::integer_math::lcm_integer;

#[cfg(test)]
use super::types::DdmHyperplane;
use super::types::DdmRay;

pub(super) fn active_set_for_intersection(
    p: &DdmRay,
    n: &DdmRay,
    new_hyperplane_idx: usize,
) -> Vec<usize> {
    let mut active = intersect_active(&p.active, &n.active);
    insert_active(&mut active, new_hyperplane_idx);
    active
}

pub(super) fn active_set_for_existing_ray(
    ray: &[i128],
    processed_hyperplanes: &[Vec<i128>],
) -> Vec<usize> {
    let mut active = Vec::new();
    for (idx, h) in processed_hyperplanes.iter().enumerate() {
        if dot_product(h, ray) == 0 {
            active.push(idx);
        }
    }
    active
}

pub(super) fn rational_vector_to_i128(v: &[Rational]) -> Option<Vec<i128>> {
    let mut lcm = Integer::from(1);
    for value in v {
        let (_, denominator) = value.clone().into_numerator_and_denominator();
        lcm = lcm_integer(&lcm, &Integer::from(denominator));
    }

    let mut out = Vec::with_capacity(v.len());
    for value in v {
        let scaled = value * Rational::from(&lcm);
        let integer = Integer::try_from(scaled).ok()?;
        out.push(i128::try_from(&integer).ok()?);
    }
    Some(normalize_ray_preserving_direction(out))
}

pub(super) fn insert_active(active: &mut Vec<usize>, idx: usize) {
    match active.binary_search(&idx) {
        Ok(_) => {}
        Err(pos) => active.insert(pos, idx),
    }
}

pub(super) fn merge_active(active: &mut Vec<usize>, other: &[usize]) {
    for &idx in other {
        insert_active(active, idx);
    }
}

pub(super) fn intersect_active(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut i = 0usize;
    let mut j = 0usize;
    let mut out = Vec::new();
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    out
}

pub(super) fn intersect_active_if_len_at_least(
    a: &[usize],
    b: &[usize],
    min_len: usize,
) -> Option<Vec<usize>> {
    if a.len().min(b.len()) < min_len {
        return None;
    }
    let mut i = 0usize;
    let mut j = 0usize;
    let mut out = Vec::with_capacity(min_len);
    while i < a.len() && j < b.len() {
        if out.len() + (a.len() - i).min(b.len() - j) < min_len {
            return None;
        }
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    (out.len() >= min_len).then_some(out)
}

pub(super) fn normalize_ray_preserving_direction(mut ray: Vec<i128>) -> Vec<i128> {
    let g = gcd_vec(&ray);
    if g > 0 {
        for x in &mut ray {
            *x /= g;
        }
    }
    ray
}

/// Dot product of two integer vectors.
pub(super) fn dot_product(a: &[i128], b: &[i128]) -> i128 {
    a.iter()
        .zip(b.iter())
        .try_fold(0i128, |acc, (&x, &y)| {
            x.checked_mul(y)
                .and_then(|product| acc.checked_add(product))
        })
        .expect("integer dot product overflowed i128")
}

#[cfg(test)]
pub(super) fn hyperplanes_from_dense(rows: Vec<Vec<i128>>) -> Vec<DdmHyperplane> {
    rows.into_iter().map(DdmHyperplane::new).collect()
}

/// GCD of a vector.
pub(super) fn gcd_vec(v: &[i128]) -> i128 {
    v.iter().fold(0, |acc, &x| gcd(acc, x.abs()))
}

/// GCD of two integers.
fn gcd(a: i128, b: i128) -> i128 {
    if b == 0 { a } else { gcd(b, a % b) }
}
