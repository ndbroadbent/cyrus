//! Generic CYTools ambient toric intersection computation.
//!
//! This ports `ToricVariety._construct_intnum_equations` for ambient top-form
//! intersections. The existing `compute_intersection_cytools` path is the
//! optimized four-dimensional hypersurface extraction; this module exposes the
//! upstream ambient dictionary needed by complete-intersection reductions.

use std::collections::{BTreeMap, BTreeSet};

use malachite::Rational;
use malachite::num::arithmetic::traits::Abs;
use malachite::num::conversion::traits::RoundingFrom;
use malachite::rounding_modes::RoundingMode;

use crate::Point;
use crate::error::{Error, Result};
use crate::integer_math::determinant_gaussian;
use crate::intersection::cytools_algorithm::extract::float_to_rational;
use crate::triangulation::Triangulation;
use crate::types::rational::Rational as TypedRational;
use crate::types::tags::Finite;

use super::solver::solve_sparse_system;
use crate::intersection::helpers::combinations;

const ROUND_TO_ZERO_THRESHOLD: f64 = 1e-10;

/// Ambient top-form toric intersection numbers keyed by sorted divisor indices.
#[derive(Debug, Clone)]
pub struct AmbientIntersectionNumbers {
    dimension: usize,
    divisor_count: usize,
    entries: BTreeMap<Vec<usize>, TypedRational<Finite>>,
}

impl AmbientIntersectionNumbers {
    /// Complex dimension of the ambient toric variety.
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    /// Number of divisors represented by the output, including index `0`.
    pub const fn divisor_count(&self) -> usize {
        self.divisor_count
    }

    /// Return a sorted top-form intersection value. Missing entries are zero.
    pub fn get(&self, indices: &[usize]) -> TypedRational<Finite> {
        let mut key = indices.to_vec();
        key.sort_unstable();
        self.entries
            .get(&key)
            .cloned()
            .unwrap_or_else(|| TypedRational::<Finite>::from_raw(Rational::from(0)))
    }

    /// Iterate over nonzero entries in deterministic key order.
    pub fn iter(&self) -> impl Iterator<Item = (&Vec<usize>, &TypedRational<Finite>)> {
        self.entries.iter()
    }
}

/// Compute ambient toric top-form intersections using CYTools' generic system.
///
/// The output includes prime toric divisor intersections and intersections
/// involving index `0`, where `0` denotes the canonical divisor, matching
/// CYTools' default `ToricVariety.intersection_numbers()` convention.
///
/// # Errors
/// Returns an error if the triangulation, point dimensions, or linear relations
/// are inconsistent, or if the linear system cannot be solved.
pub fn compute_ambient_intersections_cytools(
    tri: &Triangulation,
    points: &[Point],
    linear_relations_no_origin: &[Vec<i64>],
) -> Result<AmbientIntersectionNumbers> {
    let dimension = validate_ambient_input(tri, points, linear_relations_no_origin)?;
    let distinct = compute_distinct_top_intersections(tri, points, dimension)?;
    let faces_by_size = collect_faces_by_size(tri, dimension)?;
    let variable_array = build_generic_variable_array(points.len(), dimension, &faces_by_size);
    let variable_dict = variable_array
        .iter()
        .enumerate()
        .map(|(idx, key)| (key.clone(), idx))
        .collect::<BTreeMap<_, _>>();
    let c_dict = build_generic_c_dict(&distinct, dimension);
    let eqn_array = build_generic_equation_array(points.len(), dimension, &faces_by_size);
    let eqn_dict = build_generic_eqn_dict(&variable_array, &variable_dict, dimension);
    let (m_triplets, c_vec) = build_generic_linear_system(
        &eqn_array,
        &eqn_dict,
        &c_dict,
        linear_relations_no_origin,
        dimension,
    )?;

    if !variable_array.is_empty() && m_triplets.is_empty() {
        return Err(Error::SingularMatrix(
            "ambient intersection system has variables but no equations".into(),
        ));
    }

    let solution = solve_sparse_system(&m_triplets, &c_vec, variable_array.len())?;
    let mut entries = merge_prime_intersections(&distinct, &variable_array, &solution);
    add_canonical_intersections(&mut entries, dimension);
    Ok(AmbientIntersectionNumbers {
        dimension,
        divisor_count: points.len(),
        entries: rationalize_entries(entries),
    })
}

fn validate_ambient_input(
    tri: &Triangulation,
    points: &[Point],
    linear_relations_no_origin: &[Vec<i64>],
) -> Result<usize> {
    if points.is_empty() {
        return Err(Error::InvalidInput(
            "ambient intersections need points".into(),
        ));
    }
    if tri.simplices().is_empty() {
        return Err(Error::InvalidInput(
            "ambient intersections need a triangulation".into(),
        ));
    }
    let dimension = points[0].dim();
    if dimension == 0 {
        return Err(Error::InvalidInput(
            "ambient point dimension must be positive".into(),
        ));
    }
    for (idx, point) in points.iter().enumerate() {
        if point.dim() != dimension {
            return Err(Error::InvalidInput(format!(
                "point {idx} has dimension {}, expected {dimension}",
                point.dim()
            )));
        }
    }
    let expected_cols = points.len().saturating_sub(1);
    for (row_idx, row) in linear_relations_no_origin.iter().enumerate() {
        if row.len() != expected_cols {
            return Err(Error::InvalidInput(format!(
                "linear relation row {row_idx} has {} columns, expected {expected_cols}",
                row.len()
            )));
        }
    }
    Ok(dimension)
}

fn compute_distinct_top_intersections(
    tri: &Triangulation,
    points: &[Point],
    dimension: usize,
) -> Result<BTreeMap<Vec<usize>, f64>> {
    let mut distinct = BTreeMap::new();
    for simplex in tri.simplices() {
        validate_star_simplex(simplex, dimension)?;
        let mut key = non_origin_simplex(simplex, dimension)?;
        key.sort_unstable();
        distinct
            .entry(key)
            .or_insert(simplex_inverse_augmented_det(simplex, points)?);
    }
    Ok(distinct)
}

fn validate_star_simplex(simplex: &[usize], dimension: usize) -> Result<()> {
    if simplex.len() != dimension + 1 {
        return Err(Error::InvalidInput(format!(
            "simplex {simplex:?} has length {}, expected {}",
            simplex.len(),
            dimension + 1
        )));
    }
    if !simplex.contains(&0) {
        return Err(Error::InvalidInput(format!(
            "CYTools ambient intersection port expects star simplices containing origin index 0, got {simplex:?}"
        )));
    }
    Ok(())
}

fn non_origin_simplex(simplex: &[usize], dimension: usize) -> Result<Vec<usize>> {
    let non_origin = simplex
        .iter()
        .copied()
        .filter(|&idx| idx != 0)
        .collect::<Vec<_>>();
    if non_origin.len() != dimension {
        return Err(Error::InvalidInput(format!(
            "simplex {simplex:?} has {} non-origin rays, expected {dimension}",
            non_origin.len()
        )));
    }
    Ok(non_origin)
}

fn simplex_inverse_augmented_det(simplex: &[usize], points: &[Point]) -> Result<f64> {
    let mut matrix = simplex
        .iter()
        .map(|&point_idx| {
            let mut row = points[point_idx]
                .coords()
                .iter()
                .map(|&coord| Rational::from(coord))
                .collect::<Vec<_>>();
            row.push(Rational::from(1));
            row
        })
        .collect::<Vec<_>>();
    let det = determinant_gaussian(&mut matrix);
    if det == 0 {
        return Err(Error::InvalidInput(format!(
            "simplex {simplex:?} has zero augmented determinant"
        )));
    }
    let (abs_det_f64, _) = f64::rounding_from(&det.abs(), RoundingMode::Nearest);
    if abs_det_f64 <= 0.0 || !abs_det_f64.is_finite() {
        return Err(Error::InvalidInput(format!(
            "simplex {simplex:?} determinant is not a positive finite value"
        )));
    }
    Ok(1.0 / abs_det_f64)
}

fn collect_faces_by_size(
    tri: &Triangulation,
    dimension: usize,
) -> Result<Vec<BTreeSet<Vec<usize>>>> {
    let mut faces_by_size = vec![BTreeSet::new(); dimension + 1];
    for simplex in tri.simplices() {
        let non_origin = simplex
            .iter()
            .copied()
            .filter(|&idx| idx != 0)
            .collect::<Vec<_>>();
        if non_origin.len() != dimension {
            return Err(Error::InvalidInput(format!(
                "simplex {simplex:?} has {} non-origin rays, expected {dimension}",
                non_origin.len()
            )));
        }
        for size in 2..dimension {
            for mut face in combinations(&non_origin, size) {
                face.sort_unstable();
                faces_by_size[size].insert(face);
            }
        }
    }
    Ok(faces_by_size)
}

fn build_generic_variable_array(
    n_points: usize,
    dimension: usize,
    faces_by_size: &[BTreeSet<Vec<usize>>],
) -> Vec<Vec<usize>> {
    let mut variables = BTreeSet::new();
    for idx in 1..n_points {
        variables.insert(vec![idx; dimension]);
    }
    for support_size in 2..dimension {
        let patterns = repetition_choice_patterns(dimension, support_size);
        for support in &faces_by_size[support_size] {
            for pattern in &patterns {
                variables.insert(pattern.iter().map(|&choice| support[choice]).collect());
            }
        }
    }
    variables.into_iter().collect()
}

fn build_generic_c_dict(
    distinct: &BTreeMap<Vec<usize>, f64>,
    dimension: usize,
) -> BTreeMap<Vec<usize>, Vec<(usize, f64)>> {
    let mut c_dict = BTreeMap::<Vec<usize>, Vec<(usize, f64)>>::new();
    for (indices, &value) in distinct {
        for remove_pos in 0..dimension {
            let mut key = indices
                .iter()
                .enumerate()
                .filter_map(|(idx, &entry)| (idx != remove_pos).then_some(entry))
                .collect::<Vec<_>>();
            key.sort_unstable();
            c_dict
                .entry(key)
                .or_default()
                .push((indices[remove_pos], value));
        }
    }
    c_dict
}

fn build_generic_equation_array(
    n_points: usize,
    dimension: usize,
    faces_by_size: &[BTreeSet<Vec<usize>>],
) -> Vec<Vec<usize>> {
    let mut equations = BTreeSet::new();
    if dimension > 1 {
        equations.extend(faces_by_size[dimension - 1].iter().cloned());
    }
    for idx in 1..n_points {
        equations.insert(vec![idx; dimension - 1]);
    }
    for support_size in 2..dimension.saturating_sub(1) {
        let patterns = repetition_choice_patterns(dimension - 1, support_size);
        for support in &faces_by_size[support_size] {
            for pattern in &patterns {
                equations.insert(pattern.iter().map(|&choice| support[choice]).collect());
            }
        }
    }
    equations.into_iter().collect()
}

fn build_generic_eqn_dict(
    variable_array: &[Vec<usize>],
    variable_dict: &BTreeMap<Vec<usize>, usize>,
    dimension: usize,
) -> BTreeMap<Vec<usize>, Vec<(usize, usize)>> {
    let mut eqn_dict = BTreeMap::<Vec<usize>, Vec<(usize, usize)>>::new();
    for variable in variable_array {
        let var_idx = variable_dict[variable];
        let subtuples = combinations(variable, dimension - 1)
            .into_iter()
            .collect::<BTreeSet<_>>();
        for subtuple in subtuples {
            let omitted = omitted_value(variable, &subtuple)
                .expect("subtuple produced by combinations removes one variable entry");
            eqn_dict
                .entry(subtuple)
                .or_default()
                .push((omitted, var_idx));
        }
    }
    eqn_dict
}

fn omitted_value(variable: &[usize], subtuple: &[usize]) -> Option<usize> {
    for remove_pos in 0..variable.len() {
        let reduced = variable
            .iter()
            .enumerate()
            .filter_map(|(idx, &value)| (idx != remove_pos).then_some(value))
            .collect::<Vec<_>>();
        if reduced == subtuple {
            return Some(variable[remove_pos]);
        }
    }
    None
}

fn build_generic_linear_system(
    eqn_array: &[Vec<usize>],
    eqn_dict: &BTreeMap<Vec<usize>, Vec<(usize, usize)>>,
    c_dict: &BTreeMap<Vec<usize>, Vec<(usize, f64)>>,
    linear_relations_no_origin: &[Vec<i64>],
    dimension: usize,
) -> Result<(Vec<(usize, usize, f64)>, Vec<f64>)> {
    let mut triplets = Vec::new();
    let mut c_vec = vec![0.0; linear_relations_no_origin.len() * eqn_array.len()];
    let mut row_ctr = 0;
    for eqn in eqn_array {
        for lin in linear_relations_no_origin {
            if is_distinct_tuple(eqn, dimension - 1)
                && let Some(c_entries) = c_dict.get(eqn)
            {
                for &(point_idx, value) in c_entries {
                    c_vec[row_ctr] += linear_relation_coeff(lin, point_idx)? * value;
                }
            }
            if let Some(entries) = eqn_dict.get(eqn) {
                for &(point_idx, var_idx) in entries {
                    let coeff = linear_relation_coeff(lin, point_idx)?;
                    if coeff.abs() > 1e-12 {
                        triplets.push((row_ctr, var_idx, coeff));
                    }
                }
            }
            row_ctr += 1;
        }
    }
    Ok((triplets, c_vec))
}

fn is_distinct_tuple(indices: &[usize], expected_len: usize) -> bool {
    indices.len() == expected_len
        && indices.iter().copied().collect::<BTreeSet<_>>().len() == expected_len
}

fn linear_relation_coeff(row: &[i64], point_idx: usize) -> Result<f64> {
    if point_idx == 0 || point_idx > row.len() {
        return Err(Error::InvalidInput(format!(
            "linear relation coefficient requested for point index {point_idx}, but row has {} non-origin columns",
            row.len()
        )));
    }
    Ok(row[point_idx - 1] as f64)
}

fn merge_prime_intersections(
    distinct: &BTreeMap<Vec<usize>, f64>,
    variable_array: &[Vec<usize>],
    solution: &[f64],
) -> BTreeMap<Vec<usize>, f64> {
    let mut entries = distinct.clone();
    for (idx, variable) in variable_array.iter().enumerate() {
        let value = solution[idx];
        if value.abs() > ROUND_TO_ZERO_THRESHOLD {
            entries.insert(variable.clone(), value);
        }
    }
    entries
}

fn add_canonical_intersections(entries: &mut BTreeMap<Vec<usize>, f64>, dimension: usize) {
    let mut previous_level = BTreeMap::<Vec<usize>, f64>::new();
    for (indices, &value) in entries.iter() {
        for choice in remove_one_position_choices(indices) {
            let mut key = vec![0];
            key.extend(choice);
            *previous_level.entry(key).or_insert(0.0) -= value;
        }
    }
    retain_nonzero(&mut previous_level);
    entries.extend(previous_level.clone());

    for zero_count in 2..=dimension {
        let mut next_level = BTreeMap::<Vec<usize>, f64>::new();
        for (indices, &value) in &previous_level {
            let tail = &indices[zero_count - 1..];
            for choice in remove_one_position_choices(tail) {
                let mut key = vec![0; zero_count];
                key.extend(choice);
                *next_level.entry(key).or_insert(0.0) -= value;
            }
        }
        retain_nonzero(&mut next_level);
        entries.extend(next_level.clone());
        previous_level = next_level;
    }
}

fn remove_one_position_choices(indices: &[usize]) -> BTreeSet<Vec<usize>> {
    let mut choices = BTreeSet::new();
    for remove_pos in 0..indices.len() {
        choices.insert(
            indices
                .iter()
                .enumerate()
                .filter_map(|(idx, &value)| (idx != remove_pos).then_some(value))
                .collect(),
        );
    }
    choices
}

fn retain_nonzero(entries: &mut BTreeMap<Vec<usize>, f64>) {
    entries.retain(|_, value| value.abs() > ROUND_TO_ZERO_THRESHOLD);
}

fn rationalize_entries(
    entries: BTreeMap<Vec<usize>, f64>,
) -> BTreeMap<Vec<usize>, TypedRational<Finite>> {
    entries
        .into_iter()
        .filter(|(_, value)| value.abs() > ROUND_TO_ZERO_THRESHOLD)
        .map(|(key, value)| {
            (
                key,
                TypedRational::<Finite>::from_raw(float_to_rational(value)),
            )
        })
        .collect()
}

fn repetition_choice_patterns(output_len: usize, support_size: usize) -> Vec<Vec<usize>> {
    assert!(support_size >= 1);
    assert!(support_size <= output_len);
    let repeated_positions = output_len - support_size;
    combinations(
        &(0..output_len.saturating_sub(1)).collect::<Vec<_>>(),
        repeated_positions,
    )
    .into_iter()
    .map(|repeat_set| {
        let repeat_set = repeat_set.into_iter().collect::<BTreeSet<_>>();
        let mut pattern = vec![0; output_len];
        for idx in 1..output_len {
            pattern[idx] = pattern[idx - 1] + usize::from(!repeat_set.contains(&(idx - 1)));
        }
        pattern
    })
    .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repetition_patterns_match_cytools_dim4_cases() {
        assert_eq!(
            repetition_choice_patterns(4, 3),
            vec![vec![0, 0, 1, 2], vec![0, 1, 1, 2], vec![0, 1, 2, 2]]
        );
        assert_eq!(
            repetition_choice_patterns(4, 2),
            vec![vec![0, 0, 0, 1], vec![0, 0, 1, 1], vec![0, 1, 1, 1]]
        );
    }

    #[test]
    fn p4_ambient_intersections_include_canonical_divisor() {
        let points = p4_points();
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

        let ambient = compute_ambient_intersections_cytools(&tri, &points, &linrels).unwrap();

        assert_eq!(ambient.dimension(), 4);
        assert_eq!(ambient.divisor_count(), 6);
        assert_eq!(ambient.get(&[1, 1, 1, 1]).get(), &Rational::from(1));
        assert_eq!(ambient.get(&[1, 2, 3, 4]).get(), &Rational::from(1));
        assert_eq!(ambient.get(&[0, 2, 3, 4]).get(), &Rational::from(-5));
        assert_eq!(ambient.get(&[0, 0, 2, 3]).get(), &Rational::from(25));
        assert_eq!(ambient.get(&[0, 0, 0, 2]).get(), &Rational::from(-125));
        assert_eq!(ambient.get(&[0, 0, 0, 0]).get(), &Rational::from(625));
    }

    #[test]
    fn p4_anticanonical_cicy_reduction_matches_quintic_triples() {
        let points = p4_points();
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

        let ambient = compute_ambient_intersections_cytools(&tri, &points, &linrels).unwrap();
        let ambient_entries = ambient
            .iter()
            .map(|(key, value)| (key.clone(), value.clone()))
            .collect::<Vec<_>>();
        let quintic = crate::intersection::compute_complete_intersection_cy3_intersection_numbers(
            &ambient_entries,
            &[vec![1, 2, 3, 4, 5]],
            ambient.divisor_count(),
        )
        .unwrap();

        assert_eq!(quintic.get(1, 1, 1).get(), &Rational::from(5));
        assert_eq!(quintic.get(1, 2, 3).get(), &Rational::from(5));
        assert_eq!(quintic.get(2, 4, 5).get(), &Rational::from(5));
    }

    fn p4_points() -> Vec<Point> {
        vec![
            Point::new(vec![0, 0, 0, 0]),
            Point::new(vec![-1, -1, -1, -1]),
            Point::new(vec![1, 0, 0, 0]),
            Point::new(vec![0, 1, 0, 0]),
            Point::new(vec![0, 0, 1, 0]),
            Point::new(vec![0, 0, 0, 1]),
        ]
    }
}
