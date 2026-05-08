//! Secondary-cone hyperplanes for regular triangulations.
//!
//! This ports the native CYTools `Triangulation.secondary_cone` circuit step:
//! every pair of adjacent simplices in the triangulation's affine span gives
//! one circuit relation, which is a hyperplane normal of the secondary cone.

use std::collections::{BTreeSet, HashSet};

use malachite::{Integer, Rational};

use super::Triangulation;
use crate::Point;
use crate::error::{Error, Result};
use crate::integer_math::{integer_kernel, matrix_rank};
use crate::polytope::Polytope;
use crate::types::f64::F64;
use crate::types::tags::{Finite, Pos};

/// Which side of an affine circuit a triangulation selects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CircuitOmissionSide {
    /// The triangulation contains all circuit facets obtained by omitting a
    /// positive-coefficient circuit point.
    PositiveCoefficient,
    /// The triangulation contains all circuit facets obtained by omitting a
    /// negative-coefficient circuit point.
    NegativeCoefficient,
    /// No full circuit-omission side appears in the supplied triangulation.
    None,
    /// The supplied triangulation contains facets from both sides, or only a
    /// partial side.
    MixedOrPartial,
}

impl CircuitOmissionSide {
    /// Stable diagnostic string for JSON reports and tests.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::PositiveCoefficient => "positive_coefficient_omission_side",
            Self::NegativeCoefficient => "negative_coefficient_omission_side",
            Self::None => "no_circuit_omission_side",
            Self::MixedOrPartial => "mixed_or_partial_circuit_omission_side",
        }
    }
}

/// Side classification for one affine circuit against a triangulation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CircuitOmissionSideClassification {
    /// Circuit points with positive relation coefficients.
    pub positive_coefficient_points: Vec<usize>,
    /// Circuit points with negative relation coefficients.
    pub negative_coefficient_points: Vec<usize>,
    /// Positive-coefficient points whose omission facet appears.
    pub positive_coefficient_omission_hits: Vec<usize>,
    /// Negative-coefficient points whose omission facet appears.
    pub negative_coefficient_omission_hits: Vec<usize>,
    /// Classified side selected by the supplied triangulation.
    pub side: CircuitOmissionSide,
    /// The circuit facets on the selected side, if a full side is selected.
    pub selected_side_facets: Vec<Vec<usize>>,
    /// The circuit facets on the opposite side, if a full side is selected.
    pub opposite_side_facets: Vec<Vec<usize>>,
}

/// Compute native secondary-cone hyperplane normals from adjacent simplices.
///
/// The output is in point-index order. For each pair of full-dimensional
/// simplices in the triangulation's own affine span sharing a codimension-one
/// face, this builds the homogenized circuit matrix `[p_i; 1]` and extracts
/// its primitive integer kernel vector, matching the native CYTools sign
/// convention that the first differing point has positive coefficient.
///
/// The point coordinates may live in a higher-dimensional ambient lattice than
/// the simplex dimension. This is needed for CYTools' face-restriction path,
/// where a 2-face triangulation is still expressed using ambient point labels.
///
/// # Errors
///
/// Returns an error if points have inconsistent dimensions, simplices are
/// degenerate in their affine span, or an adjacent circuit does not have a
/// unique integer relation.
pub fn secondary_cone_hyperplanes_native(
    points: &[Point],
    triangulation: &Triangulation,
) -> Result<Vec<Vec<i64>>> {
    let Some(simplex_dim) = validate_secondary_cone_input(points, triangulation)? else {
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
            if common.len() != simplex_dim {
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

/// Compute the native secondary cone of a face skeleton.
///
/// This ports the CYTools `Triangulation.secondary_cone(on_faces_dim=N)` native
/// path: restrict the triangulation to each supplied face, compute that face
/// triangulation's adjacent-simplex circuit inequalities, then deduplicate the
/// ambient point-index hyperplanes. Each face is supplied as ambient point
/// indices; its dimension is inferred from its exact affine rank.
///
/// # Errors
///
/// Returns an error if a face has invalid point indices, the triangulation does
/// not induce any full-dimensional simplex on a supplied positive-dimensional
/// face, or any induced face secondary-cone computation is invalid.
pub fn secondary_cone_hyperplanes_native_on_faces(
    points: &[Point],
    triangulation: &Triangulation,
    faces: &[Vec<usize>],
) -> Result<Vec<Vec<i64>>> {
    validate_point_dimensions(points, "secondary cone face restriction points")?;

    let mut hyperplanes = BTreeSet::new();
    for (face_idx, face) in faces.iter().enumerate() {
        let face_set = validate_face_indices(points, face, face_idx)?;
        let face_dim = affine_rank_for_point_index_set(points, &face_set)?;
        if face_dim == 0 {
            continue;
        }

        let mut restricted_simplices = BTreeSet::new();
        for simplex in triangulation.simplices() {
            let restricted = simplex
                .iter()
                .copied()
                .filter(|point_index| face_set.contains(point_index))
                .collect::<BTreeSet<_>>();
            if restricted.len() == face_dim + 1 {
                restricted_simplices.insert(restricted.into_iter().collect::<Vec<_>>());
            }
        }
        if restricted_simplices.is_empty() {
            return Err(Error::InvalidInput(format!(
                "secondary cone face {face_idx} has no induced {face_dim}-simplices"
            )));
        }

        let restricted = Triangulation::new(restricted_simplices.into_iter().collect());
        for hyperplane in secondary_cone_hyperplanes_native(points, &restricted)? {
            hyperplanes.insert(hyperplane);
        }
    }

    Ok(hyperplanes.into_iter().collect())
}

/// Compute the native secondary cone restricted to 4D polytope two-faces.
///
/// This is the reusable CYTools `on_faces_dim=2` path for four-dimensional
/// reflexive polytopes: first construct CYTools-style two-face point-index
/// sets from the polytope and caller-supplied triangulation points, then
/// compute the induced face-skeleton circuit hyperplanes.
///
/// # Errors
///
/// Returns an error if the polytope face construction fails or if the supplied
/// triangulation does not induce valid full-dimensional simplices on a
/// positive-dimensional two-face.
pub fn secondary_cone_hyperplanes_native_on_polytope_2faces_4d(
    points: &[Point],
    triangulation: &Triangulation,
    polytope: &Polytope,
) -> Result<Vec<Vec<i64>>> {
    let faces = polytope.faces_4d_for_points(points)?;
    secondary_cone_hyperplanes_native_on_faces(points, triangulation, &faces.twofaces)
}

/// Construct the circuit facets for one side of a bistellar flip.
///
/// For a circuit `C = C_+ union C_-`, the two triangulations of the circuit are
/// `{C \ {p}: p in C_+}` and `{C \ {p}: p in C_-}`. This helper returns one
/// of those two facet sets in sorted point-index order.
///
/// # Errors
///
/// Returns an error if the circuit has duplicate points, zero coefficients, or
/// lacks either a positive or negative side.
pub fn circuit_omission_facets(
    circuit: &[(usize, i64)],
    side: CircuitOmissionSide,
) -> Result<Vec<Vec<usize>>> {
    let (all_points, positive_points, negative_points) = validate_sparse_circuit(circuit)?;
    let omitted_points = match side {
        CircuitOmissionSide::PositiveCoefficient => &positive_points,
        CircuitOmissionSide::NegativeCoefficient => &negative_points,
        CircuitOmissionSide::None | CircuitOmissionSide::MixedOrPartial => {
            return Ok(Vec::new());
        }
    };
    let mut facets = omitted_points
        .iter()
        .map(|omitted| {
            all_points
                .iter()
                .copied()
                .filter(|point| point != omitted)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    facets.sort();
    facets.dedup();
    Ok(facets)
}

/// Classify which side of a circuit is selected by a triangulation.
///
/// A side is selected when every omission facet from that side is contained in
/// at least one supplied simplex and no omission facet from the opposite side
/// appears. This is the local GKZ circuit-side test used by bistellar flips and
/// secondary-fan wall crossing.
///
/// # Errors
///
/// Returns an error if the circuit is malformed or if a simplex contains
/// duplicate point indices.
pub fn classify_circuit_omission_side(
    simplices: &[Vec<usize>],
    circuit: &[(usize, i64)],
) -> Result<CircuitOmissionSideClassification> {
    let (all_points, positive_points, negative_points) = validate_sparse_circuit(circuit)?;
    validate_simplices_for_circuit_side(simplices)?;

    let omission_hits = |points: &[usize]| {
        points
            .iter()
            .copied()
            .filter(|omitted| {
                let required = all_points
                    .iter()
                    .copied()
                    .filter(|point| point != omitted)
                    .collect::<BTreeSet<_>>();
                simplices.iter().any(|simplex| {
                    let simplex = simplex.iter().copied().collect::<BTreeSet<_>>();
                    required.iter().all(|point| simplex.contains(point))
                })
            })
            .collect::<Vec<_>>()
    };
    let positive_hits = omission_hits(&positive_points);
    let negative_hits = omission_hits(&negative_points);
    let side = if positive_hits.len() == positive_points.len() && negative_hits.is_empty() {
        CircuitOmissionSide::PositiveCoefficient
    } else if negative_hits.len() == negative_points.len() && positive_hits.is_empty() {
        CircuitOmissionSide::NegativeCoefficient
    } else if positive_hits.is_empty() && negative_hits.is_empty() {
        CircuitOmissionSide::None
    } else {
        CircuitOmissionSide::MixedOrPartial
    };
    let opposite_side = match side {
        CircuitOmissionSide::PositiveCoefficient => CircuitOmissionSide::NegativeCoefficient,
        CircuitOmissionSide::NegativeCoefficient => CircuitOmissionSide::PositiveCoefficient,
        CircuitOmissionSide::None | CircuitOmissionSide::MixedOrPartial => {
            CircuitOmissionSide::None
        }
    };
    Ok(CircuitOmissionSideClassification {
        positive_coefficient_points: positive_points,
        negative_coefficient_points: negative_points,
        positive_coefficient_omission_hits: positive_hits,
        negative_coefficient_omission_hits: negative_hits,
        side,
        selected_side_facets: circuit_omission_facets(circuit, side)?,
        opposite_side_facets: circuit_omission_facets(circuit, opposite_side)?,
    })
}

fn validate_secondary_cone_input(
    points: &[Point],
    triangulation: &Triangulation,
) -> Result<Option<usize>> {
    validate_point_dimensions(points, "secondary cone points")?;
    if points.is_empty() || triangulation.simplices().is_empty() {
        return Ok(None);
    }
    let simplex_len = triangulation.simplices()[0].len();
    if simplex_len == 0 {
        return Err(Error::InvalidInput(
            "secondary cone simplex is empty".to_string(),
        ));
    }
    let simplex_dim = simplex_len - 1;
    for simplex in triangulation.simplices() {
        if simplex.len() != simplex_len {
            return Err(Error::InvalidInput(format!(
                "secondary cone simplex has {} vertices, expected {}",
                simplex.len(),
                simplex_len
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
        let simplex_set = simplex.iter().copied().collect::<BTreeSet<_>>();
        let affine_rank = affine_rank_for_point_index_set(points, &simplex_set)?;
        if affine_rank != simplex_dim {
            return Err(Error::InvalidInput(format!(
                "secondary cone simplex affine rank {affine_rank}, expected {simplex_dim}"
            )));
        }
    }
    Ok(Some(simplex_dim))
}

fn validate_point_dimensions(points: &[Point], context: &str) -> Result<()> {
    let Some(first_point) = points.first() else {
        return Ok(());
    };
    let dim = first_point.dim();
    if points.iter().any(|point| point.dim() != dim) {
        return Err(Error::InvalidInput(format!(
            "{context} have inconsistent dimensions"
        )));
    }
    Ok(())
}

fn validate_face_indices(
    points: &[Point],
    face: &[usize],
    face_idx: usize,
) -> Result<BTreeSet<usize>> {
    if face.is_empty() {
        return Err(Error::InvalidInput(format!(
            "secondary cone face {face_idx} is empty"
        )));
    }
    let mut face_set = BTreeSet::new();
    for &point_index in face {
        if point_index >= points.len() {
            return Err(Error::InvalidInput(format!(
                "secondary cone face {face_idx} point index {point_index} exceeds point count {}",
                points.len()
            )));
        }
        if !face_set.insert(point_index) {
            return Err(Error::InvalidInput(format!(
                "secondary cone face {face_idx} contains duplicate point index {point_index}"
            )));
        }
    }
    Ok(face_set)
}

fn validate_sparse_circuit(
    circuit: &[(usize, i64)],
) -> Result<(Vec<usize>, Vec<usize>, Vec<usize>)> {
    if circuit.is_empty() {
        return Err(Error::InvalidInput("circuit is empty".to_string()));
    }
    let mut seen = BTreeSet::new();
    let mut all_points = Vec::with_capacity(circuit.len());
    let mut positive_points = Vec::new();
    let mut negative_points = Vec::new();
    for &(point, coefficient) in circuit {
        if coefficient == 0 {
            return Err(Error::InvalidInput(format!(
                "circuit point {point} has zero coefficient"
            )));
        }
        if !seen.insert(point) {
            return Err(Error::InvalidInput(format!(
                "circuit contains duplicate point index {point}"
            )));
        }
        all_points.push(point);
        if coefficient > 0 {
            positive_points.push(point);
        } else {
            negative_points.push(point);
        }
    }
    if positive_points.is_empty() || negative_points.is_empty() {
        return Err(Error::InvalidInput(
            "circuit must contain positive and negative coefficients".to_string(),
        ));
    }
    all_points.sort_unstable();
    positive_points.sort_unstable();
    negative_points.sort_unstable();
    Ok((all_points, positive_points, negative_points))
}

fn validate_simplices_for_circuit_side(simplices: &[Vec<usize>]) -> Result<()> {
    for (simplex_idx, simplex) in simplices.iter().enumerate() {
        let mut seen = BTreeSet::new();
        for &point in simplex {
            if !seen.insert(point) {
                return Err(Error::InvalidInput(format!(
                    "circuit side simplex {simplex_idx} contains duplicate point index {point}"
                )));
            }
        }
    }
    Ok(())
}

fn affine_rank_for_point_index_set(points: &[Point], indices: &BTreeSet<usize>) -> Result<usize> {
    let Some(&base_idx) = indices.iter().next() else {
        return Ok(0);
    };
    let base = points
        .get(base_idx)
        .ok_or_else(|| Error::InvalidInput("affine rank base index out of bounds".into()))?;
    let rows = indices
        .iter()
        .skip(1)
        .map(|&point_idx| {
            let point = points.get(point_idx).ok_or_else(|| {
                Error::InvalidInput(format!("affine rank point index {point_idx} out of bounds"))
            })?;
            if point.dim() != base.dim() {
                return Err(Error::InvalidInput(
                    "affine rank points have inconsistent dimensions".into(),
                ));
            }
            Ok(point
                .coords()
                .iter()
                .zip(base.coords().iter())
                .map(|(&coord, &base_coord)| Rational::from(coord - base_coord))
                .collect::<Vec<_>>())
        })
        .collect::<Result<Vec<_>>>()?;
    Ok(matrix_rank(&rows))
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
    fn embedded_square_face_secondary_cone_uses_ambient_point_indices() {
        let points = vec![
            Point::new(vec![0, 0, 0]),
            Point::new(vec![1, 0, 0]),
            Point::new(vec![1, 1, 0]),
            Point::new(vec![0, 1, 0]),
        ];
        let triangulation = Triangulation::new(vec![vec![0, 1, 2], vec![0, 2, 3]]);

        let hyperplanes = secondary_cone_hyperplanes_native(&points, &triangulation).unwrap();

        assert_eq!(hyperplanes, vec![vec![-1, 1, -1, 1]]);
    }

    #[test]
    fn face_skeleton_secondary_cone_restricts_ambient_triangulation() {
        let points = vec![
            Point::new(vec![0, 0, 0]),
            Point::new(vec![1, 0, 0]),
            Point::new(vec![1, 1, 0]),
            Point::new(vec![0, 1, 0]),
            Point::new(vec![0, 0, 1]),
        ];
        let triangulation = Triangulation::new(vec![vec![0, 1, 2, 4], vec![0, 2, 3, 4]]);

        let hyperplanes = secondary_cone_hyperplanes_native_on_faces(
            &points,
            &triangulation,
            &[vec![0, 1, 2, 3]],
        )
        .unwrap();

        assert_eq!(hyperplanes, vec![vec![-1, 1, -1, 1, 0]]);
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
    fn circuit_omission_side_classifies_square_diagonal() {
        let simplices = vec![vec![0, 1, 2], vec![0, 2, 3]];
        let circuit = vec![(0, -1), (1, 1), (2, -1), (3, 1)];

        let classification = classify_circuit_omission_side(&simplices, &circuit).unwrap();

        assert_eq!(
            classification.side,
            CircuitOmissionSide::PositiveCoefficient
        );
        assert_eq!(
            classification.side.as_str(),
            "positive_coefficient_omission_side"
        );
        assert_eq!(classification.positive_coefficient_points, vec![1, 3]);
        assert_eq!(classification.negative_coefficient_points, vec![0, 2]);
        assert_eq!(
            classification.positive_coefficient_omission_hits,
            vec![1, 3]
        );
        assert!(classification.negative_coefficient_omission_hits.is_empty());
        assert_eq!(
            classification.selected_side_facets,
            vec![vec![0, 1, 2], vec![0, 2, 3]]
        );
        assert_eq!(
            classification.opposite_side_facets,
            vec![vec![0, 1, 3], vec![1, 2, 3]]
        );
    }

    #[test]
    fn circuit_omission_side_classifies_flipped_square_diagonal() {
        let simplices = vec![vec![0, 1, 3], vec![1, 2, 3]];
        let circuit = vec![(0, -1), (1, 1), (2, -1), (3, 1)];

        let classification = classify_circuit_omission_side(&simplices, &circuit).unwrap();

        assert_eq!(
            classification.side,
            CircuitOmissionSide::NegativeCoefficient
        );
        assert_eq!(
            classification.negative_coefficient_omission_hits,
            vec![0, 2]
        );
        assert!(classification.positive_coefficient_omission_hits.is_empty());
        assert_eq!(
            circuit_omission_facets(&circuit, CircuitOmissionSide::NegativeCoefficient).unwrap(),
            vec![vec![0, 1, 3], vec![1, 2, 3]]
        );
    }

    #[test]
    fn circuit_omission_side_rejects_malformed_circuit() {
        let duplicate =
            classify_circuit_omission_side(&[vec![0, 1]], &[(0, -1), (0, 1)]).unwrap_err();
        assert!(duplicate.to_string().contains("duplicate point index 0"));

        let one_sided =
            circuit_omission_facets(&[(0, 1), (1, 1)], CircuitOmissionSide::PositiveCoefficient)
                .unwrap_err();
        assert!(one_sided.to_string().contains("positive and negative"));
    }

    #[test]
    fn secondary_cone_rejects_degenerate_simplices() {
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![2, 0]),
        ];
        let triangulation = Triangulation::new(vec![vec![0, 1, 2]]);

        let error = secondary_cone_hyperplanes_native(&points, &triangulation).unwrap_err();

        assert!(error.to_string().contains("affine rank 1, expected 2"));
    }

    #[test]
    fn polytope_2face_secondary_cone_uses_computed_4d_faces() {
        let points = vec![
            Point::new(vec![1, 0, 0, 0]),
            Point::new(vec![0, 1, 0, 0]),
            Point::new(vec![0, 0, 1, 0]),
            Point::new(vec![0, 0, 0, 1]),
            Point::new(vec![-1, -1, -1, -1]),
        ];
        let polytope = Polytope::from_vertices(points.clone()).unwrap();
        let triangulation = Triangulation::new(vec![vec![0, 1, 2, 3, 4]]);

        let hyperplanes = secondary_cone_hyperplanes_native_on_polytope_2faces_4d(
            &points,
            &triangulation,
            &polytope,
        )
        .unwrap();

        assert!(hyperplanes.is_empty());
    }
}
