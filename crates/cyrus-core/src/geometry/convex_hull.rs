//! Convex hull computation using the Beneath-Beyond algorithm.
//!
//! Implements exact convex hull computation for d-dimensional point sets.

use malachite::Integer;
use std::collections::{HashMap, HashSet};

use super::facet::Facet;
use super::hyperplane::Hyperplane;
use super::orientation::{Orientation, determinant_exact, orientation};

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
        let (d, _n) = validate_points(points)?;

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
            let visible = visible_facets(&facets, point);

            if visible.is_empty() {
                // Point is inside or on the hull, skip
                continue;
            }

            // Find horizon ridges (shared by exactly one visible facet)
            let horizon_ridges = horizon_ridges(&facets, &visible);

            // Remove visible facets (in reverse order to maintain indices)
            remove_visible_facets(&mut facets, &visible);

            // Create new facets from horizon ridges to new point
            let (interior_sum, interior_count) = compute_centroid_sum(points, &in_hull);
            add_new_facets(
                &mut facets,
                points,
                idx,
                horizon_ridges,
                &interior_sum,
                &interior_count,
            );

            in_hull.insert(idx);
        }

        // Collect vertex indices
        let vertex_indices = collect_vertex_indices(&facets);

        Some(Self {
            points: points.to_vec(),
            facets,
            vertex_indices,
        })
    }
}

fn validate_points(points: &[Vec<i64>]) -> Option<(usize, usize)> {
    if points.is_empty() {
        return None;
    }
    let d = points[0].len();
    let n = points.len();
    if n <= d {
        return None;
    }
    Some((d, n))
}

fn visible_facets(facets: &[Facet], point: &[i64]) -> Vec<usize> {
    facets
        .iter()
        .enumerate()
        .filter(|(_, f)| f.is_visible(point))
        .map(|(i, _)| i)
        .collect()
}

fn horizon_ridges(facets: &[Facet], visible: &[usize]) -> Vec<(Vec<usize>, usize)> {
    let mut ridge_count: HashMap<Vec<usize>, (usize, usize)> = HashMap::new();
    for &fi in visible {
        for ridge in facets[fi].ridges() {
            ridge_count
                .entry(ridge)
                .and_modify(|(count, _)| *count += 1)
                .or_insert((1, fi));
        }
    }
    ridge_count
        .into_iter()
        .filter(|(_, (count, _))| *count == 1)
        .map(|(ridge, (_, fi))| (ridge, fi))
        .collect()
}

fn remove_visible_facets(facets: &mut Vec<Facet>, visible: &[usize]) {
    let mut visible_sorted = visible.to_vec();
    visible_sorted.sort_unstable();
    for &fi in visible_sorted.iter().rev() {
        facets.swap_remove(fi);
    }
}

fn add_new_facets(
    facets: &mut Vec<Facet>,
    points: &[Vec<i64>],
    point_idx: usize,
    horizon_ridges: Vec<(Vec<usize>, usize)>,
    interior_sum: &[Integer],
    interior_count: &Integer,
) {
    for (ridge, _parent_facet) in horizon_ridges {
        let mut vertices = ridge.clone();
        vertices.push(point_idx);

        let facet_points: Vec<&[i64]> = vertices.iter().map(|&i| points[i].as_slice()).collect();
        if let Some(hp) = Hyperplane::from_points(&facet_points) {
            let (normal, constant) =
                orient_outward_exact(hp.normal, hp.constant, interior_sum, interior_count);
            facets.push(Facet {
                vertices,
                normal,
                constant,
            });
        }
    }
}

fn collect_vertex_indices(facets: &[Facet]) -> Vec<usize> {
    let mut vertex_indices: HashSet<usize> = HashSet::new();
    for facet in facets {
        for &v in &facet.vertices {
            vertex_indices.insert(v);
        }
    }
    let mut vertex_indices: Vec<usize> = vertex_indices.into_iter().collect();
    vertex_indices.sort_unstable();
    vertex_indices
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
        let mut vertices: Vec<usize> = simplex.to_vec();
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
    fn test_convex_hull_2d_square() {
        // Square with vertices (-1,-1), (1,-1), (1,1), (-1,1)
        let points = vec![vec![-1i64, -1], vec![1, -1], vec![1, 1], vec![-1, 1]];

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
