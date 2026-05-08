//! Complete-intersection intersection-number reduction.
//!
//! This ports the complete-intersection branch of
//! `cytools.calabiyau.CalabiYau.intersection_numbers`: ambient top-form
//! intersections are successively intersected with each nef-partition divisor
//! class until a Calabi-Yau threefold triple-intersection tensor remains.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use malachite::Rational as MalachiteRational;
use malachite::num::basic::traits::Zero;

use crate::Point;
use crate::error::{Error, Result};
use crate::triangulation::Triangulation;
use crate::types::rational::Rational as TypedRational;
use crate::types::tags::Finite;

use super::Intersection;
use super::cytools_algorithm::compute_ambient_intersections_cytools;

/// Compute a CICY threefold triple-intersection tensor from ambient toric data.
///
/// This composes the CYTools generic ambient top-form solver with the CYTools
/// complete-intersection reduction. It still requires a source-derived
/// `nef_parts` input; it does not infer or choose a nef partition.
///
/// # Errors
/// Returns an error if ambient top-form computation or CICY reduction fails.
pub fn compute_complete_intersection_cy3_from_ambient_cytools(
    tri: &Triangulation,
    points: &[Point],
    linear_relations_no_origin: &[Vec<i64>],
    nef_parts: &[Vec<usize>],
) -> Result<Intersection> {
    let ambient = compute_ambient_intersections_cytools(tri, points, linear_relations_no_origin)?;
    let ambient_entries = ambient
        .iter()
        .map(|(key, value)| (key.clone(), value.clone()))
        .collect::<Vec<_>>();
    compute_complete_intersection_cy3_intersection_numbers(
        &ambient_entries,
        nef_parts,
        ambient.divisor_count(),
    )
}

/// Reduce ambient top-form intersections to a CICY threefold intersection tensor.
///
/// `ambient_intersections` must contain sorted ambient intersection tuples of
/// length `3 + nef_parts.len()`. Each nef part is a set of ambient divisor
/// indices whose sum defines one complete-intersection divisor class.
///
/// The reduction is the CYTools algorithm: for each nef part and each current
/// tuple, remove one occurrence of any divisor index belonging to that part,
/// canonicalize the remaining tuple, and add the ambient coefficient once per
/// distinct resulting tuple.
///
/// # Errors
/// Returns an error if tuple lengths, divisor indices, or nef partition data are
/// inconsistent.
pub fn compute_complete_intersection_cy3_intersection_numbers(
    ambient_intersections: &[(Vec<usize>, TypedRational<Finite>)],
    nef_parts: &[Vec<usize>],
    divisor_count: usize,
) -> Result<Intersection> {
    if divisor_count == 0 {
        return Err(Error::InvalidInput("divisor count is zero".into()));
    }
    if nef_parts.is_empty() {
        return Err(Error::InvalidInput(
            "complete-intersection reduction requires at least one nef part".into(),
        ));
    }
    validate_nef_parts(nef_parts, divisor_count)?;

    let ambient_tuple_len = nef_parts
        .len()
        .checked_add(3)
        .ok_or_else(|| Error::InvalidInput("ambient tuple length overflowed".into()))?;
    let mut current =
        ambient_intersection_map(ambient_intersections, ambient_tuple_len, divisor_count)?;
    let nef_part_sets = nef_parts
        .iter()
        .map(|part| part.iter().copied().collect::<HashSet<_>>())
        .collect::<Vec<_>>();
    for (step, part) in nef_part_sets.iter().enumerate() {
        let tuple_len = ambient_tuple_len
            .checked_sub(step)
            .expect("step never exceeds ambient tuple length");
        current = apply_nef_part_reduction(&current, part, tuple_len)?;
    }
    intersection_from_reduced_triples(current, divisor_count)
}

fn ambient_intersection_map(
    ambient_intersections: &[(Vec<usize>, TypedRational<Finite>)],
    ambient_tuple_len: usize,
    divisor_count: usize,
) -> Result<BTreeMap<Vec<usize>, MalachiteRational>> {
    let mut current = BTreeMap::<Vec<usize>, MalachiteRational>::new();
    for (indices, value) in ambient_intersections {
        validate_ambient_tuple(indices, ambient_tuple_len, divisor_count)?;
        if *value.get() == MalachiteRational::ZERO {
            continue;
        }
        let mut key = indices.clone();
        key.sort_unstable();
        let entry = current
            .entry(key)
            .or_insert_with(|| MalachiteRational::ZERO);
        *entry += value.get().clone();
    }
    Ok(current)
}

fn apply_nef_part_reduction(
    current: &BTreeMap<Vec<usize>, MalachiteRational>,
    part: &HashSet<usize>,
    tuple_len: usize,
) -> Result<BTreeMap<Vec<usize>, MalachiteRational>> {
    let mut next = BTreeMap::<Vec<usize>, MalachiteRational>::new();
    for (indices, value) in current {
        if indices.len() != tuple_len {
            return Err(Error::InvalidInput(format!(
                "complete-intersection reduction expected tuple length {tuple_len}, got {}",
                indices.len()
            )));
        }
        for choice in remove_one_nef_part_index_choices(indices, part) {
            let entry = next
                .entry(choice)
                .or_insert_with(|| MalachiteRational::ZERO);
            *entry += value.clone();
        }
    }
    next.retain(|_, value| *value != MalachiteRational::ZERO);
    Ok(next)
}

fn intersection_from_reduced_triples(
    triples: BTreeMap<Vec<usize>, MalachiteRational>,
    divisor_count: usize,
) -> Result<Intersection> {
    let mut out = Intersection::new(divisor_count);
    for (indices, value) in triples {
        if indices.len() != 3 {
            return Err(Error::InvalidInput(format!(
                "complete-intersection reduction should produce triples, got length {}",
                indices.len()
            )));
        }
        out.set(
            indices[0],
            indices[1],
            indices[2],
            TypedRational::<Finite>::from_raw(value),
        );
    }
    Ok(out)
}

fn validate_nef_parts(nef_parts: &[Vec<usize>], divisor_count: usize) -> Result<()> {
    let mut seen = HashSet::new();
    for part in nef_parts {
        if part.is_empty() {
            return Err(Error::InvalidInput(
                "nef partition contains an empty part".into(),
            ));
        }
        for &idx in part {
            if idx >= divisor_count {
                return Err(Error::InvalidInput(format!(
                    "nef partition index {idx} is out of range for divisor count {divisor_count}"
                )));
            }
            if !seen.insert(idx) {
                return Err(Error::InvalidInput(format!(
                    "nef partition index {idx} appears in more than one part"
                )));
            }
        }
    }
    Ok(())
}

fn validate_ambient_tuple(
    indices: &[usize],
    ambient_tuple_len: usize,
    divisor_count: usize,
) -> Result<()> {
    if indices.len() != ambient_tuple_len {
        return Err(Error::InvalidInput(format!(
            "ambient intersection tuple length must be {ambient_tuple_len}, got {}",
            indices.len()
        )));
    }
    for &idx in indices {
        if idx >= divisor_count {
            return Err(Error::InvalidInput(format!(
                "ambient intersection index {idx} is out of range for divisor count {divisor_count}"
            )));
        }
    }
    Ok(())
}

fn remove_one_nef_part_index_choices(
    indices: &[usize],
    part: &HashSet<usize>,
) -> BTreeSet<Vec<usize>> {
    let mut choices = BTreeSet::new();
    for remove_pos in 0..indices.len() {
        if !part.contains(&indices[remove_pos]) {
            continue;
        }
        let mut reduced = indices
            .iter()
            .enumerate()
            .filter_map(|(idx, &value)| (idx != remove_pos).then_some(value))
            .collect::<Vec<_>>();
        reduced.sort_unstable();
        choices.insert(reduced);
    }
    choices
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(value: i64) -> TypedRational<Finite> {
        TypedRational::from_raw(MalachiteRational::from(value))
    }

    #[test]
    fn cicy_reduction_applies_single_nef_part_once_per_distinct_choice() {
        let ambient = vec![
            (vec![0, 0, 0, 3], q(-1)),
            (vec![0, 1, 1, 3], q(2)),
            (vec![0, 1, 3, 3], q(4)),
            (vec![1, 1, 2, 2], q(7)),
        ];

        let reduced =
            compute_complete_intersection_cy3_intersection_numbers(&ambient, &[vec![3]], 4)
                .unwrap();

        assert_eq!(reduced.get(0, 0, 0).get(), &MalachiteRational::from(-1));
        assert_eq!(reduced.get(0, 1, 1).get(), &MalachiteRational::from(2));
        assert_eq!(reduced.get(0, 1, 3).get(), &MalachiteRational::from(4));
        assert_eq!(reduced.get(1, 1, 2).get(), &MalachiteRational::ZERO);
    }

    #[test]
    fn cicy_reduction_threads_multiple_nef_parts_in_order() {
        let ambient = vec![
            (vec![0, 1, 2, 3, 4], q(6)),
            (vec![0, 1, 3, 4, 4], q(2)),
            (vec![0, 2, 2, 3, 3], q(9)),
        ];

        let reduced = compute_complete_intersection_cy3_intersection_numbers(
            &ambient,
            &[vec![3], vec![4]],
            5,
        )
        .unwrap();

        assert_eq!(reduced.get(0, 1, 2).get(), &MalachiteRational::from(6));
        assert_eq!(reduced.get(0, 1, 4).get(), &MalachiteRational::from(2));
        assert_eq!(reduced.get(0, 2, 2).get(), &MalachiteRational::ZERO);
    }

    #[test]
    fn cicy_reduction_rejects_overlapping_nef_parts() {
        let ambient = vec![(vec![0, 1, 2, 3, 4], q(1))];

        let err = compute_complete_intersection_cy3_intersection_numbers(
            &ambient,
            &[vec![3], vec![3]],
            5,
        )
        .unwrap_err();

        assert!(
            err.to_string()
                .contains("nef partition index 3 appears in more than one part")
        );
    }

    #[test]
    fn cicy_from_ambient_cytools_matches_quintic_triples() {
        let points = vec![
            Point::new(vec![0, 0, 0, 0]),
            Point::new(vec![-1, -1, -1, -1]),
            Point::new(vec![1, 0, 0, 0]),
            Point::new(vec![0, 1, 0, 0]),
            Point::new(vec![0, 0, 1, 0]),
            Point::new(vec![0, 0, 0, 1]),
        ];
        let tri = Triangulation::new(vec![
            vec![0, 1, 2, 3, 4],
            vec![0, 1, 2, 3, 5],
            vec![0, 1, 2, 4, 5],
            vec![0, 1, 3, 4, 5],
            vec![0, 2, 3, 4, 5],
        ]);
        let linrels = vec![
            vec![-1, 1, 0, 0, 0],
            vec![-1, 0, 1, 0, 0],
            vec![-1, 0, 0, 1, 0],
            vec![-1, 0, 0, 0, 1],
        ];

        let quintic = compute_complete_intersection_cy3_from_ambient_cytools(
            &tri,
            &points,
            &linrels,
            &[vec![1, 2, 3, 4, 5]],
        )
        .unwrap();

        assert_eq!(quintic.get(1, 1, 1).get(), &MalachiteRational::from(5));
        assert_eq!(quintic.get(1, 2, 3).get(), &MalachiteRational::from(5));
        assert_eq!(quintic.get(2, 4, 5).get(), &MalachiteRational::from(5));
    }
}
