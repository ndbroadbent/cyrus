//! Triangulation of lattice polytopes.
//!
//! Implements Fine, Regular, Star Triangulations (FRST) needed for
//! constructing smooth Calabi-Yau hypersurfaces.
//!
//! Uses Integer arithmetic for orientation tests (not Rational) for performance.
//! Heights are scaled by 2^40 and rounded to integers - this preserves the sign
//! of all orientation determinants which is all we need for convex hull.
//!
//! Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 4.

use crate::Point;
use crate::error::{Error, Result};
use malachite::Integer;
use rayon::prelude::*;
use std::collections::HashSet;

/// Scale factor for converting f64 heights to integers.
/// 2^40 ≈ 10^12 gives plenty of precision while staying in reasonable integer range.
const HEIGHT_SCALE: f64 = 1_099_511_627_776.0; // 2^40

/// A non-empty vector of lifted points (Integer coordinates including scaled height).
struct LiftedPoints(Vec<Vec<Integer>>);

impl LiftedPoints {
    /// Create from a non-empty vector. Returns None if empty.
    fn new(points: Vec<Vec<Integer>>) -> Option<Self> {
        if points.is_empty() {
            None
        } else {
            Some(Self(points))
        }
    }

    fn as_slice(&self) -> &[Vec<Integer>] {
        &self.0
    }

    const fn len(&self) -> usize {
        self.0.len()
    }

    fn first(&self) -> &Vec<Integer> {
        // Safe: guaranteed non-empty by construction
        &self.0[0]
    }

    fn iter(&self) -> impl Iterator<Item = &Vec<Integer>> {
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
/// Uses scaled Integer arithmetic for exact orientation tests.
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

    // 1. Lift points: (v_i, h_i) in Z^{d+1} using scaled Integer for exactness
    let lifted: Vec<Vec<Integer>> = points
        .iter()
        .zip(heights.iter())
        .map(|(p, &h)| {
            let mut v: Vec<Integer> = p.coords().iter().map(|&x| Integer::from(x)).collect();
            // Scale height to integer: round(h * 2^40)
            #[allow(clippy::cast_possible_truncation)]
            let h_scaled = (h * HEIGHT_SCALE).round() as i64;
            v.push(Integer::from(h_scaled));
            v
        })
        .collect();

    // 2. Wrap in LiftedPoints - if empty, return empty triangulation
    let Some(lifted) = LiftedPoints::new(lifted) else {
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
/// Uses rayon for parallel visibility checks.
fn convex_hull(points: &LiftedPoints) -> Vec<Facet> {
    let n = points.len();
    let dim = points.first().len();

    if n <= dim {
        // Points form a single simplex or less
        return vec![(0..n).collect()];
    }

    // 1. Find initial simplex (first dim points + 1)
    let mut current_facets: Vec<Facet> = Vec::new();

    // Create initial facets from the first dim+1 points
    for i in 0..=dim {
        let face: Vec<usize> = (0..=dim).filter(|&idx| idx != i).collect();
        current_facets.push(face);
    }

    let all_points = points.as_slice();

    // 2. Incrementally add points
    for i in (dim + 1)..n {
        let p = &all_points[i];

        // Parallel visibility check for all facets
        let visibility: Vec<(usize, bool, Vec<Vec<usize>>)> = current_facets
            .par_iter()
            .enumerate()
            .map(|(f_idx, facet)| {
                let visible = is_visible(facet, p, all_points);
                let ridges = if visible {
                    // Compute ridges for this facet
                    (0..facet.len())
                        .map(|j| {
                            let mut ridge = facet.clone();
                            ridge.remove(j);
                            ridge.sort_unstable();
                            ridge
                        })
                        .collect()
                } else {
                    Vec::new()
                };
                (f_idx, visible, ridges)
            })
            .collect();

        // Collect visible facets and ridges (sequential reduction)
        let mut visible_facets = Vec::new();
        let mut horizon_ridges = HashSet::new();

        for (f_idx, visible, ridges) in visibility {
            if visible {
                visible_facets.push(f_idx);
                for ridge in ridges {
                    if !horizon_ridges.insert(ridge.clone()) {
                        horizon_ridges.remove(&ridge);
                    }
                }
            }
        }

        if visible_facets.is_empty() {
            continue;
        }

        // Build new facet list
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
fn is_visible(facet_indices: &[usize], p: &[Integer], all_points: &[Vec<Integer>]) -> bool {
    // 1. Compute orientation of (Facet, P)
    let mut mat_p = Vec::new();
    for &idx in facet_indices {
        mat_p.push(all_points[idx].clone());
    }
    mat_p.push(p.to_vec());
    let vol_p = orientation_integer(&mat_p);

    if vol_p == Integer::from(0) {
        return false; // Coplanar
    }

    // 2. Compute orientation of (Facet, Center)
    // Use centroid of initial simplex (points 0..dim) scaled to avoid division
    let dim_space = all_points[0].len();
    let n_center = dim_space + 1;

    // Sum instead of average (scaling doesn't affect sign comparison)
    let mut center_sum = vec![Integer::from(0); dim_space];
    for i in 0..n_center {
        if i < all_points.len() {
            for (j, item) in center_sum.iter_mut().enumerate().take(dim_space) {
                *item += &all_points[i][j];
            }
        }
    }

    let mut mat_c = Vec::new();
    // Scale facet points by n_center to match center_sum scale
    for &idx in facet_indices {
        let scaled: Vec<Integer> = all_points[idx]
            .iter()
            .map(|x| x * Integer::from(n_center as i64))
            .collect();
        mat_c.push(scaled);
    }
    mat_c.push(center_sum);

    let vol_c = orientation_integer(&mat_c);

    // If signs are different, p is on opposite side from center => Visible
    (vol_p.clone() * vol_c) < Integer::from(0)
}

/// Check if a facet is a "lower" face (normal points downward in height dimension).
fn is_lower_face(facet_indices: &[usize], all_points: &LiftedPoints) -> bool {
    let dim = all_points.first().len(); // d+1

    // Find min height in the whole set
    let min_h = all_points
        .iter()
        .filter_map(|p| p.last())
        .min()
        .expect("LiftedPoints guarantees non-empty");

    // Construct a test point below the facet (way below min height)
    let mut below_point = all_points.as_slice()[facet_indices[0]].clone();
    below_point[dim - 1] = min_h - Integer::from(1_000_000_000_000_i64);

    is_visible(facet_indices, &below_point, all_points.as_slice())
}

/// Compute orientation (signed volume) of a simplex using Integer determinant.
/// Returns positive, negative, or zero.
fn orientation_integer(points: &[Vec<Integer>]) -> Integer {
    let n = points.len();
    if n == 0 {
        return Integer::from(0);
    }
    let dim = n - 1;
    if dim == 0 {
        return Integer::from(0);
    }

    // Build matrix of differences: M[i][j] = points[i+1][j] - points[0][j]
    let mut mat: Vec<Vec<Integer>> = vec![vec![Integer::from(0); dim]; dim];
    for i in 0..dim {
        for j in 0..dim {
            mat[i][j] = &points[i + 1][j] - &points[0][j];
        }
    }

    determinant_integer(&mut mat)
}

/// Compute determinant of a square Integer matrix using Bareiss algorithm.
/// This is fraction-free and avoids the expression swell of Gaussian elimination.
fn determinant_integer(mat: &mut [Vec<Integer>]) -> Integer {
    let n = mat.len();
    if n == 0 {
        return Integer::from(1);
    }
    if n == 1 {
        return mat[0][0].clone();
    }

    let mut sign = Integer::from(1);
    let mut prev_pivot = Integer::from(1);

    for k in 0..n {
        // Find pivot
        let mut pivot_row = None;
        for i in k..n {
            if mat[i][k] != Integer::from(0) {
                pivot_row = Some(i);
                break;
            }
        }

        let Some(pivot_row) = pivot_row else {
            return Integer::from(0); // Singular matrix
        };

        if pivot_row != k {
            mat.swap(k, pivot_row);
            sign = -sign;
        }

        let pivot = mat[k][k].clone();

        // Bareiss algorithm: M[i][j] = (M[i][j] * pivot - M[i][k] * M[k][j]) / prev_pivot
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = &mat[i][j] * &pivot - &mat[i][k] * &mat[k][j];
                mat[i][j] = numerator / &prev_pivot;
            }
        }

        prev_pivot = pivot;
    }

    sign * &mat[n - 1][n - 1]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn int(n: i64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn test_regular_triangulation_2d_square() {
        // Square (0,0), (1,0), (1,1), (0,1)
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![1, 1]),
            Point::new(vec![0, 1]),
        ];
        // heights: 0, 2, 6, 3
        let heights = vec![0.0, 2.0, 6.0, 3.0];

        let tri = compute_regular_triangulation(&points, &heights).unwrap();
        // Should produce 2 triangles
        assert_eq!(tri.simplices().len(), 2);
    }

    #[test]
    fn test_is_visible_2d() {
        let all_points = vec![
            vec![int(0), int(0)],
            vec![int(2), int(0)],
            vec![int(0), int(2)],
            vec![int(3), int(3)], // Outside
            vec![int(1), int(1)], // Inside (using integer approx)
        ];

        // Test P3 (outside)
        assert!(is_visible(&[1, 2], &all_points[3], &all_points));

        // Test P4 (inside) - should not be visible from edge [1,2]
        assert!(!is_visible(&[1, 2], &all_points[4], &all_points));
    }

    #[test]
    fn test_is_lower_face_3d() {
        let points = LiftedPoints::new(vec![
            vec![int(0), int(0), int(0)], // 0
            vec![int(1), int(0), int(0)], // 1
            vec![int(0), int(1), int(0)], // 2
            vec![int(0), int(0), int(1)], // 3 (Top)
        ])
        .unwrap();

        // Face 0-1-2 (Bottom)
        assert!(is_lower_face(&[0, 1, 2], &points));

        // Face 1-2-3 (Side/Top)
        assert!(!is_lower_face(&[1, 2, 3], &points));
    }

    #[test]
    fn test_lifted_points_rejects_empty() {
        let empty: Vec<Vec<Integer>> = vec![];
        assert!(LiftedPoints::new(empty).is_none());
    }

    #[test]
    fn test_lifted_points_accepts_non_empty() {
        let points = vec![vec![int(0), int(0)]];
        assert!(LiftedPoints::new(points).is_some());
    }

    #[test]
    fn test_compute_delaunay_heights() {
        let points = vec![
            Point::new(vec![0, 0]),   // |p|² = 0
            Point::new(vec![1, 0]),   // |p|² = 1
            Point::new(vec![1, 1]),   // |p|² = 2
            Point::new(vec![3, 4]),   // |p|² = 25
            Point::new(vec![-2, -2]), // |p|² = 8
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
        assert!(!tri_star.is_star(1));

        // Triangulation where point 0 is missing from one simplex
        let tri_not_star = Triangulation::new(vec![
            vec![0, 1, 2],
            vec![1, 2, 3], // Missing point 0
        ]);
        assert!(!tri_not_star.is_star(0));
    }

    #[test]
    fn test_compute_frst_heights_simple() {
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

        let result = compute_frst_heights(&points, 10);
        assert!(result.is_err());
    }

    #[test]
    fn test_determinant_integer_2x2() {
        // det([[1, 2], [3, 4]]) = 1*4 - 2*3 = -2
        let mut mat = vec![
            vec![int(1), int(2)],
            vec![int(3), int(4)],
        ];
        assert_eq!(determinant_integer(&mut mat), int(-2));
    }

    #[test]
    fn test_determinant_integer_3x3() {
        // det([[1, 2, 3], [4, 5, 6], [7, 8, 10]]) = -3
        let mut mat = vec![
            vec![int(1), int(2), int(3)],
            vec![int(4), int(5), int(6)],
            vec![int(7), int(8), int(10)],
        ];
        assert_eq!(determinant_integer(&mut mat), int(-3));
    }

    #[test]
    fn test_determinant_integer_singular() {
        // Singular matrix
        let mut mat = vec![
            vec![int(1), int(2)],
            vec![int(2), int(4)],
        ];
        assert_eq!(determinant_integer(&mut mat), int(0));
    }
}
