//! Computational geometry primitives with exact arithmetic.
//!
//! Implements convex hull and related operations for polytope computations.
//! Uses malachite::Integer for exact arithmetic at decision boundaries.
//!
//! Reference: [[project_docs/clean_room/GEOMETRIC_PRIMITIVES.md]]
//! Reference: [[project_docs/clean_room/CONVEX_HULL_D.md]]

use malachite::Integer;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

/// Result of an orientation test.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    /// Point is above the hyperplane (positive side).
    Positive,
    /// Point is below the hyperplane (negative side).
    Negative,
    /// Point is exactly on the hyperplane.
    Zero,
}

impl Orientation {
    fn from_sign(det: &Integer) -> Self {
        match det.cmp(&Integer::from(0)) {
            Ordering::Greater => Orientation::Positive,
            Ordering::Less => Orientation::Negative,
            Ordering::Equal => Orientation::Zero,
        }
    }
}

/// Compute the orientation of point p relative to hyperplane through v1..v_d.
///
/// Returns the sign of the determinant of the (d+1)×(d+1) matrix:
/// ```text
/// | 1  v1[0]  v1[1]  ...  v1[d-1] |
/// | 1  v2[0]  v2[1]  ...  v2[d-1] |
/// | ...                           |
/// | 1  vd[0]  vd[1]  ...  vd[d-1] |
/// | 1  p[0]   p[1]   ...  p[d-1]  |
/// ```
///
/// Uses exact integer arithmetic (malachite::Integer).
pub fn orientation(hyperplane_points: &[&[i64]], query_point: &[i64]) -> Orientation {
    let d = query_point.len();
    assert_eq!(hyperplane_points.len(), d, "Need exactly d points for hyperplane in d dimensions");
    for p in hyperplane_points {
        assert_eq!(p.len(), d, "All points must have same dimension");
    }

    // Build the (d+1) × (d+1) matrix
    let n = d + 1;
    let mut matrix: Vec<Vec<Integer>> = Vec::with_capacity(n);

    // First d rows from hyperplane points
    for hp in hyperplane_points {
        let mut row = Vec::with_capacity(n);
        row.push(Integer::from(1));
        for &coord in hp.iter() {
            row.push(Integer::from(coord));
        }
        matrix.push(row);
    }

    // Last row from query point
    let mut last_row = Vec::with_capacity(n);
    last_row.push(Integer::from(1));
    for &coord in query_point {
        last_row.push(Integer::from(coord));
    }
    matrix.push(last_row);

    let det = determinant_exact(&matrix);
    Orientation::from_sign(&det)
}

/// Compute exact determinant of an n×n matrix using Bareiss algorithm.
///
/// Bareiss algorithm is fraction-free Gaussian elimination that maintains
/// exact integer arithmetic throughout.
fn determinant_exact(matrix: &[Vec<Integer>]) -> Integer {
    let n = matrix.len();
    if n == 0 {
        return Integer::from(1);
    }

    // Clone matrix for in-place modification
    let mut m: Vec<Vec<Integer>> = matrix.to_vec();

    let mut sign = Integer::from(1);
    let mut prev_pivot = Integer::from(1);

    for k in 0..n {
        // Find pivot
        let mut pivot_row = None;
        for i in k..n {
            if m[i][k] != 0 {
                pivot_row = Some(i);
                break;
            }
        }

        let pivot_row = match pivot_row {
            Some(r) => r,
            None => return Integer::from(0), // Singular matrix
        };

        // Swap rows if needed
        if pivot_row != k {
            m.swap(k, pivot_row);
            sign = -sign;
        }

        let pivot = m[k][k].clone();

        // Bareiss elimination
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                // m[i][j] = (m[k][k] * m[i][j] - m[k][j] * m[i][k]) / prev_pivot
                let numerator = &pivot * &m[i][j] - &m[k][j] * &m[i][k];
                m[i][j] = numerator / &prev_pivot;
            }
            m[i][k] = Integer::from(0);
        }

        prev_pivot = pivot;
    }

    sign * &m[n - 1][n - 1]
}

/// A hyperplane in d dimensions: n · x = c
#[derive(Debug, Clone)]
pub struct Hyperplane {
    /// Normal vector (as integers).
    pub normal: Vec<Integer>,
    /// Constant term.
    pub constant: Integer,
}

impl Hyperplane {
    /// Construct hyperplane through d points in d dimensions.
    ///
    /// Uses cofactor expansion to compute the normal vector.
    pub fn from_points(points: &[&[i64]]) -> Option<Self> {
        let d = points.len();
        if d == 0 {
            return None;
        }
        for p in points {
            if p.len() != d {
                return None;
            }
        }

        // Build (d-1) difference vectors: v_i - v_0
        let v0 = points[0];
        let diffs: Vec<Vec<Integer>> = points[1..]
            .iter()
            .map(|v| {
                v.iter()
                    .zip(v0.iter())
                    .map(|(&a, &b)| Integer::from(a - b))
                    .collect()
            })
            .collect();

        // Compute normal via cofactors
        // n_i = (-1)^i * det(M_i) where M_i is diffs with column i removed
        let mut normal = Vec::with_capacity(d);
        for i in 0..d {
            let minor = extract_minor(&diffs, i);
            let det = determinant_exact(&minor);
            let sign = if i % 2 == 0 {
                Integer::from(1)
            } else {
                Integer::from(-1)
            };
            normal.push(sign * det);
        }

        // Check for degenerate case (all zeros)
        if normal.iter().all(|x| *x == 0) {
            return None;
        }

        // Compute constant: c = n · v0
        let constant: Integer = normal
            .iter()
            .zip(v0.iter())
            .map(|(n, &v)| n * Integer::from(v))
            .sum();

        Some(Hyperplane { normal, constant })
    }

    /// Check which side of the hyperplane a point is on.
    pub fn side(&self, point: &[i64]) -> Orientation {
        let dot: Integer = self
            .normal
            .iter()
            .zip(point.iter())
            .map(|(n, &p)| n * Integer::from(p))
            .sum();

        Orientation::from_sign(&(&dot - &self.constant))
    }
}

/// Extract minor matrix by removing column `col` from a (d-1)×d matrix.
fn extract_minor(matrix: &[Vec<Integer>], col: usize) -> Vec<Vec<Integer>> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .filter(|(j, _)| *j != col)
                .map(|(_, v)| v.clone())
                .collect()
        })
        .collect()
}

/// A facet of a convex hull.
#[derive(Debug, Clone)]
pub struct Facet {
    /// Indices of the d vertices defining this facet.
    pub vertices: Vec<usize>,
    /// Outward normal vector.
    pub normal: Vec<Integer>,
    /// Constant term (normal · x = constant on facet).
    pub constant: Integer,
}

impl Facet {
    /// Get the ridges (d-1 subsets of vertices) of this facet.
    fn ridges(&self) -> Vec<Vec<usize>> {
        let d = self.vertices.len();
        let mut ridges = Vec::with_capacity(d);
        for i in 0..d {
            let mut ridge: Vec<usize> = self.vertices.iter().copied().collect();
            ridge.remove(i);
            ridge.sort_unstable();
            ridges.push(ridge);
        }
        ridges
    }

    /// Check if a point is visible from this facet (strictly above).
    fn is_visible(&self, point: &[i64]) -> bool {
        let dot: Integer = self
            .normal
            .iter()
            .zip(point.iter())
            .map(|(n, &p)| n * Integer::from(p))
            .sum();
        dot > self.constant
    }

    /// Check if a point is on this facet.
    #[allow(dead_code)]
    fn is_on(&self, point: &[i64]) -> bool {
        let dot: Integer = self
            .normal
            .iter()
            .zip(point.iter())
            .map(|(n, &p)| n * Integer::from(p))
            .sum();
        dot == self.constant
    }
}

/// Convex hull of a set of points in d dimensions.
#[derive(Debug, Clone)]
pub struct ConvexHull {
    /// The original points.
    pub points: Vec<Vec<i64>>,
    /// The facets of the hull.
    pub facets: Vec<Facet>,
    /// Indices of points that are vertices of the hull.
    pub vertex_indices: Vec<usize>,
}

impl ConvexHull {
    /// Compute the convex hull of a set of points using the Beneath-Beyond algorithm.
    ///
    /// Returns None if the points are degenerate (< d+1 affinely independent points).
    pub fn compute(points: &[Vec<i64>]) -> Option<Self> {
        if points.is_empty() {
            return None;
        }

        let d = points[0].len();
        let n = points.len();

        if n <= d {
            return None; // Not enough points for a d-dimensional hull
        }

        // Find d+1 affinely independent points for initial simplex
        let initial_simplex = find_initial_simplex(points, d)?;

        // Create initial facets from the simplex
        let mut facets = create_simplex_facets(points, &initial_simplex, d);

        // Track which points are already in the hull
        let mut in_hull: HashSet<usize> = initial_simplex.iter().copied().collect();

        // Incrementally add remaining points
        for (idx, point) in points.iter().enumerate() {
            if in_hull.contains(&idx) {
                continue;
            }

            // Find visible facets
            let visible: Vec<usize> = facets
                .iter()
                .enumerate()
                .filter(|(_, f)| f.is_visible(point))
                .map(|(i, _)| i)
                .collect();

            if visible.is_empty() {
                // Point is inside or on the hull, skip
                continue;
            }

            // Find horizon ridges (shared by exactly one visible facet)
            let mut ridge_count: HashMap<Vec<usize>, (usize, usize)> = HashMap::new();
            for &fi in &visible {
                for ridge in facets[fi].ridges() {
                    ridge_count
                        .entry(ridge)
                        .and_modify(|(count, _)| *count += 1)
                        .or_insert((1, fi));
                }
            }

            let horizon_ridges: Vec<(Vec<usize>, usize)> = ridge_count
                .into_iter()
                .filter(|(_, (count, _))| *count == 1)
                .map(|(ridge, (_, fi))| (ridge, fi))
                .collect();

            // Remove visible facets (in reverse order to maintain indices)
            let mut visible_sorted = visible;
            visible_sorted.sort_unstable();
            for &fi in visible_sorted.iter().rev() {
                facets.swap_remove(fi);
            }

            // Create new facets from horizon ridges to new point
            let (interior_sum, interior_count) = compute_centroid_sum(points, &in_hull);

            for (ridge, _parent_facet) in horizon_ridges {
                let mut vertices = ridge.clone();
                vertices.push(idx);

                // Compute hyperplane for new facet
                let facet_points: Vec<&[i64]> =
                    vertices.iter().map(|&i| points[i].as_slice()).collect();

                if let Some(hp) = Hyperplane::from_points(&facet_points) {
                    // Orient normal outward (away from interior)
                    // Use centroid of existing hull vertices as interior point
                    let (normal, constant) = orient_outward_exact(
                        hp.normal,
                        hp.constant,
                        &interior_sum,
                        &interior_count,
                    );

                    facets.push(Facet {
                        vertices,
                        normal,
                        constant,
                    });
                }
            }

            in_hull.insert(idx);
        }

        // Collect vertex indices
        let mut vertex_indices: HashSet<usize> = HashSet::new();
        for facet in &facets {
            for &v in &facet.vertices {
                vertex_indices.insert(v);
            }
        }
        let mut vertex_indices: Vec<usize> = vertex_indices.into_iter().collect();
        vertex_indices.sort_unstable();

        Some(ConvexHull {
            points: points.to_vec(),
            facets,
            vertex_indices,
        })
    }
}

/// Find d+1 affinely independent points for initial simplex.
fn find_initial_simplex(points: &[Vec<i64>], d: usize) -> Option<Vec<usize>> {
    let n = points.len();
    if n <= d {
        return None;
    }

    // Start with the first point
    let mut simplex = vec![0];

    // Greedily add points that are affinely independent
    for i in 1..n {
        simplex.push(i);

        // Check if still affinely independent using orientation
        if simplex.len() == d + 1 {
            // Check if the simplex has non-zero volume
            let pts: Vec<&[i64]> = simplex[..d].iter().map(|&j| points[j].as_slice()).collect();
            let query = &points[simplex[d]];
            let orient = orientation(&pts, query);
            if orient != Orientation::Zero {
                return Some(simplex);
            }
            simplex.pop();
        } else if simplex.len() > 1 {
            // Check affine independence of current subset
            // For k points, check if they span k-1 dimensional space
            if !are_affinely_independent(points, &simplex) {
                simplex.pop();
            }
        }
    }

    // Couldn't find enough affinely independent points
    if simplex.len() == d + 1 {
        Some(simplex)
    } else {
        None
    }
}

/// Check if a set of points is affinely independent.
fn are_affinely_independent(points: &[Vec<i64>], indices: &[usize]) -> bool {
    if indices.len() <= 1 {
        return true;
    }

    let d = points[0].len();

    // Build matrix of difference vectors
    let v0 = &points[indices[0]];
    let diffs: Vec<Vec<Integer>> = indices[1..]
        .iter()
        .map(|&i| {
            points[i]
                .iter()
                .zip(v0.iter())
                .map(|(&a, &b)| Integer::from(a - b))
                .collect()
        })
        .collect();

    // Check rank by looking for non-zero minor
    // For k-1 vectors, they're independent if some (k-1)×(k-1) minor is non-zero
    let k = diffs.len();
    if k > d {
        return false; // More vectors than dimensions
    }

    // Check if the vectors span the expected dimension
    // Use the first k columns
    let minor: Vec<Vec<Integer>> = diffs
        .iter()
        .map(|row| row.iter().take(k).cloned().collect())
        .collect();

    if minor.is_empty() || minor[0].is_empty() {
        return true;
    }

    determinant_exact(&minor) != 0
}

/// Create initial facets from a simplex.
fn create_simplex_facets(points: &[Vec<i64>], simplex: &[usize], d: usize) -> Vec<Facet> {
    let mut facets = Vec::with_capacity(d + 1);

    // Each facet is the simplex with one vertex removed
    for i in 0..=d {
        let mut vertices: Vec<usize> = simplex.iter().copied().collect();
        let omitted_vertex = vertices.remove(i);

        let facet_points: Vec<&[i64]> = vertices.iter().map(|&j| points[j].as_slice()).collect();

        if let Some(hp) = Hyperplane::from_points(&facet_points) {
            // Orient normal outward (away from omitted vertex)
            // Convert single point to sum format (sum = point, count = 1)
            let omitted_point = &points[omitted_vertex];
            let interior_sum: Vec<Integer> =
                omitted_point.iter().map(|&x| Integer::from(x)).collect();
            let interior_count = Integer::from(1);
            let (normal, constant) =
                orient_outward_exact(hp.normal, hp.constant, &interior_sum, &interior_count);

            facets.push(Facet {
                vertices,
                normal,
                constant,
            });
        }
    }

    facets
}

/// Orient normal to point away from the interior point.
/// Interior is given as (sum, count) where actual point = sum/count.
/// We use exact arithmetic: instead of (sum/n) · normal > constant,
/// we check sum · normal > n * constant.
fn orient_outward_exact(
    normal: Vec<Integer>,
    constant: Integer,
    interior_sum: &[Integer],
    interior_count: &Integer,
) -> (Vec<Integer>, Integer) {
    // Compute sum · normal (this is n * centroid · normal)
    let dot: Integer = normal
        .iter()
        .zip(interior_sum.iter())
        .map(|(n, s)| n * s)
        .sum();

    // Compare: sum · normal vs n * constant
    // (equivalent to centroid · normal vs constant)
    let scaled_constant = interior_count * &constant;

    if dot > scaled_constant {
        // Interior is on positive side, flip normal
        let neg_normal: Vec<Integer> = normal.iter().map(|x| -x).collect();
        (neg_normal, -constant)
    } else {
        (normal, constant)
    }
}

/// Compute centroid sum of points at given indices (as Integer for exact arithmetic).
/// Returns (sum, count) where actual centroid = sum / count.
/// We avoid division to maintain exactness.
fn compute_centroid_sum(points: &[Vec<i64>], indices: &HashSet<usize>) -> (Vec<Integer>, Integer) {
    let d = points[0].len();
    let n = Integer::from(indices.len());

    let mut sum = vec![Integer::from(0); d];
    for &i in indices {
        for (j, &coord) in points[i].iter().enumerate() {
            sum[j] += Integer::from(coord);
        }
    }
    (sum, n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_orientation_2d() {
        // Triangle with vertices (0,0), (1,0), (0,1)
        // Point (0.5, 0.5) should be... let's compute
        // Actually for integer points, test with (1,1) vs the line from (0,0) to (2,0)
        let hp = [&[0i64, 0][..], &[2, 0][..]];

        // Point above line
        assert_eq!(orientation(&hp, &[1, 1]), Orientation::Positive);
        // Point below line
        assert_eq!(orientation(&hp, &[1, -1]), Orientation::Negative);
        // Point on line
        assert_eq!(orientation(&hp, &[1, 0]), Orientation::Zero);
    }

    #[test]
    fn test_determinant_exact() {
        // 2x2 identity
        let m = vec![
            vec![Integer::from(1), Integer::from(0)],
            vec![Integer::from(0), Integer::from(1)],
        ];
        assert_eq!(determinant_exact(&m), Integer::from(1));

        // 2x2 with det = ad - bc = 1*4 - 2*3 = -2
        let m = vec![
            vec![Integer::from(1), Integer::from(2)],
            vec![Integer::from(3), Integer::from(4)],
        ];
        assert_eq!(determinant_exact(&m), Integer::from(-2));
    }

    #[test]
    fn test_hyperplane_from_points_2d() {
        // Line through (0,0) and (1,0)
        let hp = Hyperplane::from_points(&[&[0i64, 0][..], &[1, 0][..]]).unwrap();
        // Normal should be perpendicular to (1,0), i.e., (0, ±1)
        assert_eq!(hp.normal[0], Integer::from(0));
        assert!(hp.normal[1] != Integer::from(0));
    }

    #[test]
    fn test_convex_hull_2d_square() {
        // Square with vertices (-1,-1), (1,-1), (1,1), (-1,1)
        let points = vec![
            vec![-1i64, -1],
            vec![1, -1],
            vec![1, 1],
            vec![-1, 1],
        ];

        let hull = ConvexHull::compute(&points).unwrap();

        // Should have 4 vertices
        assert_eq!(hull.vertex_indices.len(), 4);
        // Should have 4 facets (edges in 2D)
        assert_eq!(hull.facets.len(), 4);
    }

    #[test]
    fn test_convex_hull_2d_with_interior() {
        // Square with an interior point at origin
        let points = vec![
            vec![-1i64, -1],
            vec![1, -1],
            vec![1, 1],
            vec![-1, 1],
            vec![0, 0], // interior point
        ];

        let hull = ConvexHull::compute(&points).unwrap();

        // Should still have only 4 vertices (origin is interior)
        assert_eq!(hull.vertex_indices.len(), 4);
        assert!(!hull.vertex_indices.contains(&4)); // Index 4 (origin) not a vertex
    }

    #[test]
    fn test_convex_hull_3d_tetrahedron() {
        // Regular-ish tetrahedron
        let points = vec![
            vec![0i64, 0, 0],
            vec![1, 0, 0],
            vec![0, 1, 0],
            vec![0, 0, 1],
        ];

        let hull = ConvexHull::compute(&points).unwrap();

        // Should have 4 vertices and 4 facets
        assert_eq!(hull.vertex_indices.len(), 4);
        assert_eq!(hull.facets.len(), 4);
    }
}
