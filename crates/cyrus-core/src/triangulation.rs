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

/// A non-empty vector of points. Guarantees at least one element.
struct NonEmptyPoints(Vec<Vec<Rational>>);

impl NonEmptyPoints {
    /// Create from a non-empty vector. Returns None if empty.
    fn new(points: Vec<Vec<Rational>>) -> Option<Self> {
        if points.is_empty() {
            None
        } else {
            Some(Self(points))
        }
    }

    fn as_slice(&self) -> &[Vec<Rational>] {
        &self.0
    }

    const fn len(&self) -> usize {
        self.0.len()
    }

    fn first(&self) -> &Vec<Rational> {
        // Safe: guaranteed non-empty by construction
        &self.0[0]
    }

    fn iter(&self) -> impl Iterator<Item = &Vec<Rational>> {
        self.0.iter()
    }
}

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

    /// Check if this triangulation has the Star property.
    ///
    /// A triangulation is Star if a specified point (typically the origin)
    /// is a vertex of every simplex.
    pub fn is_star(&self, point_idx: usize) -> bool {
        self.simplices.iter().all(|s| s.contains(&point_idx))
    }
}

/// Compute Delaunay heights for a set of points.
///
/// For each point p, the height is h(p) = |p|² = p · p (squared Euclidean norm).
/// This is the standard choice for inducing a Delaunay-like triangulation.
pub fn compute_delaunay_heights(points: &[Point]) -> Vec<f64> {
    points
        .iter()
        .map(|p| {
            p.coords()
                .iter()
                .map(|&x| (x as f64) * (x as f64))
                .sum()
        })
        .collect()
}

/// Compute default FRST (Fine Regular Star Triangulation) heights.
///
/// Starting from Delaunay heights (|p|²), adjusts the origin's height
/// downward until the triangulation has the Star property (origin in every simplex).
///
/// # Arguments
/// * `points` - The lattice points to triangulate
/// * `origin_idx` - Index of the origin point (must be in points)
///
/// # Returns
/// * `Ok((heights, triangulation))` - The adjusted heights and resulting triangulation
/// * `Err` if origin_idx is out of bounds or triangulation fails
///
/// Reference: CYTools algorithm (Demirtas et al. 2022, Sec. 3.2)
pub fn compute_frst_heights(
    points: &[Point],
    origin_idx: usize,
) -> Result<(Vec<f64>, Triangulation)> {
    if origin_idx >= points.len() {
        return Err(Error::InvalidInput(format!(
            "origin_idx {} out of bounds for {} points",
            origin_idx,
            points.len()
        )));
    }

    // Start with Delaunay heights
    let mut heights = compute_delaunay_heights(points);

    // Maximum iterations to prevent infinite loop
    const MAX_ITERATIONS: usize = 100;

    for _ in 0..MAX_ITERATIONS {
        let tri = compute_regular_triangulation(points, &heights)?;

        if tri.is_star(origin_idx) {
            return Ok((heights, tri));
        }

        // Adjust origin height downward
        let min_h = heights
            .iter()
            .copied()
            .fold(f64::INFINITY, |a, b| a.min(b));
        let max_h = heights
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, |a, b| a.max(b));

        // h(origin) -= (max - min + 10)
        heights[origin_idx] -= (max_h - min_h) + 10.0;
    }

    Err(Error::InvalidInput(
        "Failed to achieve Star triangulation after max iterations".into(),
    ))
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

    // 2. Wrap in NonEmptyPoints - if empty, return empty triangulation
    let Some(lifted) = NonEmptyPoints::new(lifted) else {
        return Ok(Triangulation::new(Vec::new()));
    };

    // 3. Compute Convex Hull in d+1 dimensions
    let facets = convex_hull(&lifted);

    // 4. Extract lower faces
    let mut simplices = Vec::new();
    for facet in facets {
        if is_lower_face(&facet, &lifted) {
            simplices.push(facet);
        }
    }

    Ok(Triangulation::new(simplices))
}

/// A facet of the convex hull (indices of vertices).
type Facet = Vec<usize>;

/// Compute the convex hull of a set of points in N dimensions.
fn convex_hull(points: &NonEmptyPoints) -> Vec<Facet> {
    let n = points.len();
    let dim = points.first().len();

    if n <= dim {
        // Points form a single simplex or less
        return vec![(0..n).collect()];
    }

    // 1. Find initial simplex (first dim points + 1)
    let mut current_facets: Vec<Facet> = Vec::new();

    // Create initial facets from the first dim+1 points
    // INVARIANT: n > dim (checked above), so n >= dim+1, meaning i in 0..=dim
    // will always satisfy i < n. No break needed.
    for i in 0..=dim {
        let face: Vec<usize> = (0..=dim).filter(|&idx| idx != i).collect();
        current_facets.push(face);
    }

    // 2. Incrementally add points
    for i in (dim + 1)..n {
        let p = &points.as_slice()[i];

        let mut visible_facets = Vec::new();
        let mut horizon_ridges = HashSet::new();

        for (f_idx, facet) in current_facets.iter().enumerate() {
            if is_visible(facet, p, points.as_slice()) {
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
/// Caller guarantees all_points is non-empty and indices are valid.
fn is_visible(facet_indices: &[usize], p: &[Rational], all_points: &[Vec<Rational>]) -> bool {
    // 1. Compute orientation of (Facet, P)
    let mut mat_p = Vec::new();
    for &idx in facet_indices {
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
fn is_lower_face(facet_indices: &[usize], all_points: &NonEmptyPoints) -> bool {
    let dim = all_points.first().len(); // d+1

    // Find min height in the whole set.
    // Safe: NonEmptyPoints guarantees non-empty, and each point has height (last coord)
    let min_h = all_points
        .iter()
        .filter_map(|p| p.last())
        .min()
        .expect("NonEmptyPoints guarantees non-empty");

    // Construct a test point below the facet.
    let mut below_point = all_points.as_slice()[facet_indices[0]].clone();
    below_point[dim - 1] = min_h.clone() - Rational::from(1000); // Way below

    is_visible(facet_indices, &below_point, all_points.as_slice())
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
        let points = NonEmptyPoints::new(vec![
            vec![r(0), r(0), r(0)], // 0
            vec![r(1), r(0), r(0)], // 1
            vec![r(0), r(1), r(0)], // 2
            vec![r(0), r(0), r(1)], // 3 (Top)
        ])
        .unwrap();

        // Face 0-1-2 (Bottom)
        assert!(is_lower_face(&[0, 1, 2], &points));

        // Face 1-2-3 (Side/Top)
        assert!(!is_lower_face(&[1, 2, 3], &points));
    }

    #[test]
    fn test_non_empty_points_rejects_empty() {
        let empty: Vec<Vec<Rational>> = vec![];
        assert!(NonEmptyPoints::new(empty).is_none());
    }

    #[test]
    fn test_non_empty_points_accepts_non_empty() {
        let points = vec![vec![r(0), r(0)]];
        assert!(NonEmptyPoints::new(points).is_some());
    }

    #[test]
    fn test_compute_delaunay_heights() {
        let points = vec![
            Point::new(vec![0, 0]),    // |p|² = 0
            Point::new(vec![1, 0]),    // |p|² = 1
            Point::new(vec![1, 1]),    // |p|² = 2
            Point::new(vec![3, 4]),    // |p|² = 25
            Point::new(vec![-2, -2]),  // |p|² = 8
        ];

        let heights = compute_delaunay_heights(&points);

        assert_eq!(heights.len(), 5);
        assert!((heights[0] - 0.0).abs() < 1e-10);
        assert!((heights[1] - 1.0).abs() < 1e-10);
        assert!((heights[2] - 2.0).abs() < 1e-10);
        assert!((heights[3] - 25.0).abs() < 1e-10);
        assert!((heights[4] - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_is_star() {
        // Triangulation where point 0 is in all simplices
        let tri_star = Triangulation::new(vec![
            vec![0, 1, 2],
            vec![0, 2, 3],
            vec![0, 3, 1],
        ]);
        assert!(tri_star.is_star(0));
        assert!(!tri_star.is_star(1)); // Point 1 is not in simplex [0, 2, 3]

        // Triangulation where point 0 is missing from one simplex
        let tri_not_star = Triangulation::new(vec![
            vec![0, 1, 2],
            vec![1, 2, 3], // Missing point 0
        ]);
        assert!(!tri_not_star.is_star(0));
    }

    #[test]
    fn test_compute_frst_heights_simple() {
        // Simple 2D case: square with origin at center
        // Origin should already be in all simplices with Delaunay heights
        let points = vec![
            Point::new(vec![0, 0]),   // Origin at index 0
            Point::new(vec![-1, -1]),
            Point::new(vec![1, -1]),
            Point::new(vec![1, 1]),
            Point::new(vec![-1, 1]),
        ];

        let result = compute_frst_heights(&points, 0);
        assert!(result.is_ok());

        let (heights, tri) = result.unwrap();
        assert_eq!(heights.len(), 5);
        assert!(tri.is_star(0));
    }

    #[test]
    fn test_compute_frst_heights_invalid_origin() {
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
        ];

        let result = compute_frst_heights(&points, 10); // Invalid index
        assert!(result.is_err());
    }
}
