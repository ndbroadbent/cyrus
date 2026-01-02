//! Triangulation of lattice polytopes.
//!
//! Implements Fine, Regular, Star Triangulations (FRST) needed for
//! constructing smooth Calabi-Yau hypersurfaces.
//!
//! Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 4.

use crate::Point;
use crate::error::{Error, Result};
use crate::integer_math::orientation;
use malachite::Rational;
use malachite::num::conversion::traits::FromSciString;
use std::collections::HashSet;

/// A triangulation of a point set.
///
/// Represented as a collection of simplices, where each simplex is a
/// set of indices into the original point set.
#[derive(Debug, Clone)]
pub struct Triangulation {
    /// Indices of points forming each simplex.
    /// For a d-dimensional triangulation, each simplex has d+1 indices.
    simplices: Vec<Vec<usize>>,
}

impl Triangulation {
    /// Create a new triangulation from a list of simplices.
    pub const fn new(simplices: Vec<Vec<usize>>) -> Self {
        Self { simplices }
    }

    /// Get the simplices.
    pub fn simplices(&self) -> &[Vec<usize>] {
        &self.simplices
    }
}

/// Compute a regular triangulation of a set of points using the lifting map.
///
/// Given points S and heights H, the triangulation is the projection of
/// the lower convex hull of the lifted points (p, h_p).
///
/// # Arguments
/// * `points` - Lattice points to triangulate.
/// * `heights` - Heights used for the lifting map.
///
/// # Errors
/// Returns an error if the number of points and heights do not match.
///
/// Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 4.1
pub fn compute_regular_triangulation(points: &[Point], heights: &[f64]) -> Result<Triangulation> {
    if points.len() != heights.len() {
        return Err(Error::InvalidInput(
            "Number of points and heights must match".into(),
        ));
    }

    if points.is_empty() {
        return Ok(Triangulation::new(Vec::new()));
    }

    // 1. Lift points: (v_i, h_i) in R^{d+1} using Rational for exactness
    let lifted: Vec<Vec<Rational>> = points
        .iter()
        .zip(heights.iter())
        .map(|(p, &h)| {
            let mut v: Vec<Rational> = p.coords().iter().map(|&x| Rational::from(x)).collect();
            v.push(
                Rational::from_sci_string(&format!("{h:e}")).unwrap_or_else(|| Rational::from(0)),
            );
            v
        })
        .collect();

    // 2. Compute Convex Hull in d+1 dimensions
    let facets = convex_hull(&lifted);

    // 3. Extract lower faces
    let mut simplices = Vec::new();
    for facet in facets {
        // We propagate error if is_lower_face fails (e.g. empty points)
        // But here we know points is not empty because we checked above.
        if is_lower_face(&facet, &lifted).unwrap_or(false) {
            simplices.push(facet);
        }
    }

    Ok(Triangulation::new(simplices))
}

/// A facet of the convex hull (indices of vertices).
type Facet = Vec<usize>;

/// Compute the convex hull of a set of points in N dimensions.
fn convex_hull(points: &[Vec<Rational>]) -> Vec<Facet> {
    let n = points.len();
    if n == 0 {
        return Vec::new();
    }
    let dim = points[0].len();

    if n <= dim {
        // Points form a single simplex or less
        return vec![(0..n).collect()];
    }

    // 1. Find initial simplex (first dim points + 1)
    let mut current_facets: Vec<Facet> = Vec::new();

    // Create initial facets from the first dim+1 points
    for i in 0..=dim {
        if i >= n {
            break;
        }
        let face: Vec<usize> = (0..=dim).filter(|&idx| idx != i).collect();
        current_facets.push(face);
    }

    // 2. Incrementally add points
    for i in (dim + 1)..n {
        let p = &points[i];

        let mut visible_facets = Vec::new();
        let mut horizon_ridges = HashSet::new();

        for (f_idx, facet) in current_facets.iter().enumerate() {
            if is_visible(facet, p, points) {
                visible_facets.push(f_idx);
                // Add ridges
                for j in 0..facet.len() {
                    let mut ridge = facet.clone();
                    ridge.remove(j);
                    ridge.sort_unstable();
                    if !horizon_ridges.insert(ridge.clone()) {
                        horizon_ridges.remove(&ridge);
                    }
                }
            }
        }

        if visible_facets.is_empty() {
            continue;
        }

        let mut new_facets = Vec::with_capacity(current_facets.len());
        for (idx, facet) in current_facets.iter().enumerate() {
            if !visible_facets.contains(&idx) {
                new_facets.push(facet.clone());
            }
        }

        for ridge in horizon_ridges {
            let mut new_facet = ridge;
            new_facet.push(i);
            new_facets.push(new_facet);
        }

        current_facets = new_facets;
    }

    current_facets
}

/// Check if a point p is "visible" from the outside of a facet.
fn is_visible(facet_indices: &[usize], p: &[Rational], all_points: &[Vec<Rational>]) -> bool {
    if all_points.is_empty() {
        return false;
    }
    // 1. Compute orientation of (Facet, P)
    let mut mat_p = Vec::new();
    for &idx in facet_indices {
        if idx >= all_points.len() {
            return false; // Invalid index
        }
        mat_p.push(all_points[idx].clone());
    }
    mat_p.push(p.to_vec());
    let vol_p = orientation(&mat_p);

    if vol_p == 0 {
        return false; // Coplanar
    }

    // 2. Compute orientation of (Facet, Center)
    // Use centroid of initial simplex (points 0..dim)
    let dim_space = all_points[0].len();
    let mut center = vec![Rational::from(0); dim_space];
    for i in 0..=dim_space {
        if i < all_points.len() {
            for (j, item) in center.iter_mut().enumerate().take(dim_space) {
                *item += &all_points[i][j];
            }
        }
    }
    // Divide by dim_space + 1
    let count = Rational::from((dim_space + 1) as u64);
    for item in center.iter_mut().take(dim_space) {
        *item /= &count;
    }

    let mut mat_c = Vec::new();
    for &idx in facet_indices {
        mat_c.push(all_points[idx].clone());
    }
    mat_c.push(center);

    let vol_c = orientation(&mat_c);

    // If signs are different, p is on the opposite side of the facet from center => Outside => Visible.
    (vol_p * vol_c) < 0
}

/// Check if a facet is a "lower" face.
fn is_lower_face(facet_indices: &[usize], all_points: &[Vec<Rational>]) -> Result<bool> {
    if all_points.is_empty() {
        return Err(Error::InvalidInput("No points provided".into()));
    }
    let dim = all_points[0].len(); // d+1

    // Find min height in the whole set.
    let min_h = all_points
        .iter()
        .map(|p| {
            p.last()
                .ok_or_else(|| Error::InvalidInput("Point has no height".into()))
        })
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .min()
        .ok_or_else(|| Error::InvalidInput("Empty points list".into()))?;

    // Construct a test point below the facet.
    let mut below_point = all_points[facet_indices[0]].clone();
    below_point[dim - 1] = min_h - Rational::from(1000); // Way below

    Ok(is_visible(facet_indices, &below_point, all_points))
}

#[cfg(test)]
mod tests {
    use super::*;
    use malachite::Rational;

    fn r(n: i64) -> Rational {
        Rational::from(n)
    }

    #[test]
    fn test_regular_triangulation_2d_square() {
        // Square (0,0), (1,0), (1,1), (0,1)
        // Use heights that break coplanarity: h = 2x^2 + 3y^2 + xy
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![1, 1]),
            Point::new(vec![0, 1]),
        ];
        // heights: 0, 2, 6, 3
        let heights = vec![0.0, 2.0, 6.0, 3.0];

        let tri = compute_regular_triangulation(&points, &heights).unwrap();
        // Should produce 2 triangles.
        assert_eq!(tri.simplices().len(), 2);
    }

    #[test]
    fn test_is_visible_2d() {
        let all_points = vec![
            vec![r(0), r(0)],
            vec![r(2), r(0)],
            vec![r(0), r(2)],
            vec![r(3), r(3)],               // Outside
            vec![r(1) / r(2), r(1) / r(2)], // Inside
        ];

        // Test P3 (outside)
        assert!(is_visible(&[1, 2], &all_points[3], &all_points));

        // Test P4 (inside)
        assert!(!is_visible(&[1, 2], &all_points[4], &all_points));
    }

    #[test]
    fn test_is_lower_face_3d() {
        let points = vec![
            vec![r(0), r(0), r(0)], // 0
            vec![r(1), r(0), r(0)], // 1
            vec![r(0), r(1), r(0)], // 2
            vec![r(0), r(0), r(1)], // 3 (Top)
        ];

        // Face 0-1-2 (Bottom)
        assert!(is_lower_face(&[0, 1, 2], &points).unwrap());

        // Face 1-2-3 (Side/Top)
        assert!(!is_lower_face(&[1, 2, 3], &points).unwrap());
    }
}
