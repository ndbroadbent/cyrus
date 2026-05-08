//! Secondary-cone hyperplanes for regular triangulations.
//!
//! This ports the native CYTools `Triangulation.secondary_cone` circuit step:
//! every pair of adjacent full-dimensional simplices gives one circuit
//! relation, which is a hyperplane normal of the secondary cone.

use std::collections::{BTreeSet, HashSet};

use malachite::Integer;

use super::Triangulation;
use crate::Point;
use crate::error::{Error, Result};
use crate::integer_math::integer_kernel;
use crate::types::f64::F64;
use crate::types::tags::{Finite, Pos};

/// Compute native secondary-cone hyperplane normals from adjacent simplices.
///
/// The output is in point-index order. For each pair of full-dimensional
/// simplices sharing a codimension-one face, this builds the homogenized
/// circuit matrix `[p_i; 1]` and extracts its primitive integer kernel vector,
/// matching the native CYTools sign convention that the first differing point
/// has positive coefficient.
///
/// # Errors
///
/// Returns an error if points have inconsistent dimensions, simplices are not
/// full-dimensional index sets, or an adjacent circuit does not have a unique
/// integer relation.
pub fn secondary_cone_hyperplanes_native(
    points: &[Point],
    triangulation: &Triangulation,
) -> Result<Vec<Vec<i64>>> {
    let Some(dim) = validate_secondary_cone_input(points, triangulation)? else {
        return Ok(Vec::new());
    };

    let simplices = triangulation
        .simplices()
        .iter()
        .map(|simplex| simplex.iter().copied().collect::<BTreeSet<_>>())
        .collect::<Vec<_>>();
    let mut hyperplanes = BTreeSet::new();
    for (left_idx, left) in simplices.iter().enumerate() {
        for right in simplices.iter().skip(left_idx + 1) {
            let common = left.intersection(right).copied().collect::<Vec<_>>();
            if common.len() != dim {
                continue;
            }
            let diff = left
                .symmetric_difference(right)
                .copied()
                .collect::<Vec<_>>();
            if diff.len() != 2 {
                return Err(Error::InvalidInput(
                    "secondary cone adjacent simplices have invalid symmetric difference"
                        .to_string(),
                ));
            }
            let circuit = secondary_cone_circuit_relation(points, &diff, &common)?;
            hyperplanes.insert(circuit);
        }
    }

    Ok(hyperplanes.into_iter().collect())
}

fn validate_secondary_cone_input(
    points: &[Point],
    triangulation: &Triangulation,
) -> Result<Option<usize>> {
    let Some(first_point) = points.first() else {
        return Ok(None);
    };
    let dim = first_point.dim();
    if points.iter().any(|point| point.dim() != dim) {
        return Err(Error::InvalidInput(
            "secondary cone points have inconsistent dimensions".to_string(),
        ));
    }
    for simplex in triangulation.simplices() {
        if simplex.len() != dim + 1 {
            return Err(Error::InvalidInput(format!(
                "secondary cone simplex has {} vertices, expected {}",
                simplex.len(),
                dim + 1
            )));
        }
        let mut seen = HashSet::new();
        for &point_index in simplex {
            if point_index >= points.len() {
                return Err(Error::InvalidInput(format!(
                    "secondary cone simplex point index {point_index} exceeds point count {}",
                    points.len()
                )));
            }
            if !seen.insert(point_index) {
                return Err(Error::InvalidInput(format!(
                    "secondary cone simplex contains duplicate point index {point_index}"
                )));
            }
        }
    }
    Ok(Some(dim))
}

/// Pair secondary-cone hyperplanes with a height vector.
///
/// CYTools' `check_heights` tests these pairings against zero to decide whether
/// the height vector lies strictly inside the secondary cone chamber. The
/// hyperplanes are expected in the same point-index order as the heights.
///
/// # Errors
///
/// Returns an error if any hyperplane length differs from the height-vector
/// length.
pub fn secondary_cone_height_pairings(
    hyperplanes: &[Vec<i64>],
    heights: &[F64<Finite>],
) -> Result<Vec<F64<Finite>>> {
    let mut pairings = Vec::with_capacity(hyperplanes.len());
    for (hyperplane_idx, hyperplane) in hyperplanes.iter().enumerate() {
        if hyperplane.len() != heights.len() {
            return Err(Error::InvalidInput(format!(
                "secondary cone hyperplane {hyperplane_idx} has length {}, expected {}",
                hyperplane.len(),
                heights.len()
            )));
        }
        let mut pairing = F64::<Finite>::ZERO;
        for (&coefficient, &height) in hyperplane.iter().zip(heights.iter()) {
            let coefficient = F64::<Finite>::new(coefficient as f64)
                .expect("i64 secondary cone coefficient is finite as f64");
            pairing = pairing + coefficient * height;
        }
        pairings.push(pairing);
    }
    Ok(pairings)
}

/// Check whether a height vector is strictly inside a secondary cone.
///
/// This mirrors the native CYTools `check_heights` condition: every
/// hyperplane/height pairing must be larger than a caller-supplied positive
/// tolerance. A height on a wall is rejected.
///
/// # Errors
///
/// Returns an error if the hyperplane and height-vector dimensions are
/// inconsistent.
pub fn secondary_cone_strictly_contains_height_vector(
    hyperplanes: &[Vec<i64>],
    heights: &[F64<Finite>],
    epsilon: F64<Pos>,
) -> Result<bool> {
    Ok(secondary_cone_height_pairings(hyperplanes, heights)?
        .into_iter()
        .all(|pairing| pairing.get() > epsilon.get()))
}

fn secondary_cone_circuit_relation(
    points: &[Point],
    diff_points: &[usize],
    common_points: &[usize],
) -> Result<Vec<i64>> {
    let dim = points
        .first()
        .map(Point::dim)
        .ok_or_else(|| Error::InvalidInput("secondary cone point set is empty".to_string()))?;
    let mut columns = Vec::with_capacity(diff_points.len() + common_points.len());
    columns.extend_from_slice(diff_points);
    columns.extend_from_slice(common_points);
    let mut matrix = vec![vec![Integer::from(0); columns.len()]; dim + 1];
    for (col_idx, &point_idx) in columns.iter().enumerate() {
        for (coord_idx, &coord) in points[point_idx].coords().iter().enumerate() {
            matrix[coord_idx][col_idx] = Integer::from(coord);
        }
        matrix[dim][col_idx] = Integer::from(1);
    }

    let kernel = integer_kernel(&matrix);
    if kernel.len() != 1 {
        return Err(Error::LinearAlgebra(format!(
            "secondary cone circuit kernel dimension {}, expected 1",
            kernel.len()
        )));
    }
    let mut relation = kernel
        .into_iter()
        .next()
        .expect("checked one-dimensional kernel");
    if relation.first().is_some_and(|entry| entry < &0) {
        for entry in &mut relation {
            *entry = -entry.clone();
        }
    }

    let mut full = vec![0i64; points.len()];
    for (&point_idx, coefficient) in columns.iter().zip(relation.iter()) {
        full[point_idx] = i64::try_from(coefficient).map_err(|error| {
            Error::LinearAlgebra(format!(
                "secondary cone circuit coefficient does not fit in i64: {error:?}"
            ))
        })?;
    }
    Ok(full)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn finite(value: f64) -> F64<Finite> {
        F64::<Finite>::new(value).expect("finite test value")
    }

    fn positive(value: f64) -> F64<Pos> {
        F64::<Pos>::new(value).expect("positive test value")
    }

    #[test]
    fn square_diagonal_secondary_cone_matches_circuit_relation() {
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![1, 1]),
            Point::new(vec![0, 1]),
        ];
        let triangulation = Triangulation::new(vec![vec![0, 1, 2], vec![0, 2, 3]]);

        let hyperplanes = secondary_cone_hyperplanes_native(&points, &triangulation).unwrap();

        assert_eq!(hyperplanes, vec![vec![-1, 1, -1, 1]]);
    }

    #[test]
    fn secondary_cone_height_pairings_detect_interior_and_wall() {
        let hyperplanes = vec![vec![-1, 1, -1, 1]];
        let interior_heights = vec![finite(0.0), finite(0.0), finite(0.0), finite(1.0)];
        let wall_heights = vec![finite(0.0), finite(0.0), finite(0.0), finite(0.0)];
        let outside_heights = vec![finite(0.0), finite(0.0), finite(1.0), finite(0.0)];

        let pairings = secondary_cone_height_pairings(&hyperplanes, &interior_heights).unwrap();

        assert_eq!(pairings.len(), 1);
        assert!((pairings[0].get() - 1.0).abs() < 1e-12);
        assert!(
            secondary_cone_strictly_contains_height_vector(
                &hyperplanes,
                &interior_heights,
                positive(1e-9)
            )
            .unwrap()
        );
        assert!(
            !secondary_cone_strictly_contains_height_vector(
                &hyperplanes,
                &wall_heights,
                positive(1e-9)
            )
            .unwrap()
        );
        assert!(
            !secondary_cone_strictly_contains_height_vector(
                &hyperplanes,
                &outside_heights,
                positive(1e-9)
            )
            .unwrap()
        );
    }

    #[test]
    fn secondary_cone_height_pairings_reject_mismatched_lengths() {
        let error = secondary_cone_height_pairings(
            &[vec![-1, 1, -1, 1]],
            &[finite(0.0), finite(1.0), finite(0.0)],
        )
        .unwrap_err();

        assert!(error.to_string().contains("expected 3"));
    }

    #[test]
    fn secondary_cone_rejects_non_full_dimensional_simplices() {
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![0, 1]),
        ];
        let triangulation = Triangulation::new(vec![vec![0, 1]]);

        let error = secondary_cone_hyperplanes_native(&points, &triangulation).unwrap_err();

        assert!(error.to_string().contains("expected 3"));
    }
}
