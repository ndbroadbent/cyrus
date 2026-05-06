//! Convex hull computation for lifted points.
//!
//! Implements incremental convex hull algorithm with parallel visibility checks.

use super::integer_linalg::orientation_integer;
use malachite::Integer;
use rayon::prelude::*;
use std::collections::HashSet;

/// A non-empty vector of lifted points (Integer coordinates including scaled height).
pub struct LiftedPoints(Vec<Vec<Integer>>);

impl LiftedPoints {
    /// Create from a non-empty vector. Returns None if empty.
    pub fn new(points: Vec<Vec<Integer>>) -> Option<Self> {
        if points.is_empty() {
            None
        } else {
            Some(Self(points))
        }
    }

    pub fn as_slice(&self) -> &[Vec<Integer>] {
        &self.0
    }

    pub const fn len(&self) -> usize {
        self.0.len()
    }

    pub fn first(&self) -> &Vec<Integer> {
        // Safe: guaranteed non-empty by construction
        &self.0[0]
    }

    pub fn iter(&self) -> impl Iterator<Item = &Vec<Integer>> {
        self.0.iter()
    }
}

/// A facet of the convex hull (indices of vertices).
pub type Facet = Vec<usize>;

/// Compute the convex hull of a set of points in N dimensions.
/// Uses rayon for parallel visibility checks.
pub fn convex_hull(points: &LiftedPoints) -> Vec<Facet> {
    let n = points.len();
    let dim = points.first().len();

    if n <= dim {
        // Points form a single simplex or less
        return vec![(0..n).collect()];
    }

    let mut current_facets = initial_facets(dim);

    let all_points = points.as_slice();

    // 2. Incrementally add points
    for i in (dim + 1)..n {
        let p = &all_points[i];

        let visibility = compute_visibility(&current_facets, p, all_points);
        let (visible_facets, horizon_ridges) = collect_horizon(visibility);
        if visible_facets.is_empty() {
            continue;
        }

        current_facets = rebuild_facets(&current_facets, &visible_facets, horizon_ridges, i);
    }

    current_facets
}

fn initial_facets(dim: usize) -> Vec<Facet> {
    let mut facets = Vec::new();
    for i in 0..=dim {
        let face: Vec<usize> = (0..=dim).filter(|&idx| idx != i).collect();
        facets.push(face);
    }
    facets
}

fn compute_visibility(
    facets: &[Facet],
    point: &[Integer],
    all_points: &[Vec<Integer>],
) -> Vec<(usize, bool, Vec<Vec<usize>>)> {
    facets
        .par_iter()
        .enumerate()
        .map(|(f_idx, facet)| {
            let visible = is_visible(facet, point, all_points);
            let ridges = if visible {
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
        .collect()
}

fn collect_horizon(
    visibility: Vec<(usize, bool, Vec<Vec<usize>>)>,
) -> (Vec<usize>, HashSet<Vec<usize>>) {
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
    (visible_facets, horizon_ridges)
}

fn rebuild_facets(
    current_facets: &[Facet],
    visible_facets: &[usize],
    horizon_ridges: HashSet<Vec<usize>>,
    new_point: usize,
) -> Vec<Facet> {
    let mut new_facets = Vec::with_capacity(current_facets.len());
    for (idx, facet) in current_facets.iter().enumerate() {
        if !visible_facets.contains(&idx) {
            new_facets.push(facet.clone());
        }
    }
    for ridge in horizon_ridges {
        let mut new_facet = ridge;
        new_facet.push(new_point);
        new_facets.push(new_facet);
    }
    new_facets
}

/// Check if a point p is "visible" from the outside of a facet.
pub fn is_visible(facet_indices: &[usize], p: &[Integer], all_points: &[Vec<Integer>]) -> bool {
    // 1. Compute orientation of (Facet, P)
    let mut mat_p = Vec::new();
    for &idx in facet_indices {
        mat_p.push(all_points[idx].clone());
    }
    mat_p.push(p.to_vec());
    let vol_p = orientation_integer(&mat_p);

    if vol_p == 0 {
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
    (vol_p * vol_c) < 0
}

/// Check if a facet is a "lower" face (normal points downward in height dimension).
pub fn is_lower_face(facet_indices: &[usize], all_points: &LiftedPoints) -> bool {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn int(n: i64) -> Integer {
        Integer::from(n)
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
}
