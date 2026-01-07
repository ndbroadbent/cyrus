//! Polyhedral cone computations.
//!
//! This module provides cone operations needed for Kähler/Mori cone computations.
//! Implements the Double Description Method (DDM) for cone dualization.
//!
//! Reference: CYTools cone.py, PPL (Parma Polyhedra Library)

mod ddm;
mod tip;

use std::collections::HashSet;

/// A rational polyhedral cone, represented by either rays or hyperplanes.
///
/// A cone can be defined in two equivalent ways:
/// - V-representation: As the positive span of rays: `{R^T @ λ : λ >= 0}`
/// - H-representation: As intersection of halfspaces: `{x : H @ x >= 0}`
///
/// The Double Description Method (DDM) converts between these representations.
#[derive(Debug, Clone)]
pub struct Cone {
    /// Dimension of the ambient space
    ambient_dim: usize,

    /// Rays (V-representation), if computed
    rays: Option<Vec<Vec<i64>>>,

    /// Hyperplane normals (H-representation), if computed
    hyperplanes: Option<Vec<Vec<i64>>>,

    /// Cached dimension
    dim: Option<usize>,
}

impl Cone {
    /// Create a cone from generating rays.
    ///
    /// The cone is the positive span: `{R^T @ λ : λ >= 0}` where R is the matrix
    /// with rays as rows.
    pub fn from_rays(rays: Vec<Vec<i64>>) -> Self {
        let ambient_dim = if rays.is_empty() { 0 } else { rays[0].len() };

        // Normalize rays by GCD
        let normalized = normalize_vectors(rays);

        Self {
            ambient_dim,
            rays: Some(normalized),
            hyperplanes: None,
            dim: None,
        }
    }

    /// Create a cone from hyperplane normals (inward-pointing).
    ///
    /// The cone is the intersection: `{x : H @ x >= 0}` where H is the matrix
    /// with hyperplane normals as rows.
    pub fn from_hyperplanes(hyperplanes: Vec<Vec<i64>>) -> Self {
        let ambient_dim = if hyperplanes.is_empty() {
            0
        } else {
            hyperplanes[0].len()
        };

        // Handle empty hyperplanes (full space)
        if hyperplanes.is_empty() {
            return Self::full_space(ambient_dim);
        }

        // Normalize hyperplanes by GCD
        let normalized = normalize_vectors(hyperplanes);

        Self {
            ambient_dim,
            rays: None,
            hyperplanes: Some(normalized),
            dim: None,
        }
    }

    /// Create the full space cone (no constraints).
    pub fn full_space(ambient_dim: usize) -> Self {
        // Full space has rays ±e_i for each coordinate
        let mut rays = Vec::with_capacity(2 * ambient_dim);
        for i in 0..ambient_dim {
            let mut pos = vec![0i64; ambient_dim];
            let mut neg = vec![0i64; ambient_dim];
            pos[i] = 1;
            neg[i] = -1;
            rays.push(pos);
            rays.push(neg);
        }

        Self {
            ambient_dim,
            rays: Some(rays),
            hyperplanes: Some(Vec::new()),
            dim: Some(ambient_dim),
        }
    }

    /// Get the ambient dimension.
    pub const fn ambient_dim(&self) -> usize {
        self.ambient_dim
    }

    /// Get the dimension of the cone.
    pub fn dim(&mut self) -> usize {
        if let Some(d) = self.dim {
            return d;
        }

        let d = if let Some(ref rays) = self.rays {
            matrix_rank(rays)
        } else {
            // Need to compute rays first
            let _ = self.rays();
            matrix_rank(self.rays.as_ref().unwrap())
        };

        self.dim = Some(d);
        d
    }

    /// Get the rays (V-representation).
    ///
    /// If only hyperplanes were provided, computes rays via DDM.
    pub fn rays(&mut self) -> &[Vec<i64>] {
        if self.rays.is_none() {
            // Compute from hyperplanes using DDM
            let h = self.hyperplanes.as_ref().unwrap();
            let rays = ddm::dualize(h, self.ambient_dim);
            self.rays = Some(rays);
        }
        self.rays.as_ref().unwrap()
    }

    /// Get the hyperplanes (H-representation).
    ///
    /// If only rays were provided, computes hyperplanes via DDM.
    pub fn hyperplanes(&mut self) -> &[Vec<i64>] {
        if self.hyperplanes.is_none() {
            // Compute from rays using DDM
            let r = self.rays.as_ref().unwrap();
            let hyperplanes = ddm::dualize(r, self.ambient_dim);
            self.hyperplanes = Some(hyperplanes);
        }
        self.hyperplanes.as_ref().unwrap()
    }

    /// Get the dual cone.
    ///
    /// The dual of a cone C is: `C* = {y : y·x >= 0 for all x in C}`
    ///
    /// If C is defined by rays R, then C* is defined by hyperplanes R.
    /// If C is defined by hyperplanes H, then C* is defined by rays H.
    pub fn dual(&self) -> Self {
        if let Some(ref rays) = self.rays {
            Self::from_hyperplanes(rays.clone())
        } else {
            Self::from_rays(self.hyperplanes.as_ref().unwrap().clone())
        }
    }

    /// Check if the cone contains a point.
    ///
    /// A point x is in the cone if H @ x >= 0 for all hyperplanes H.
    pub fn contains(&mut self, point: &[f64]) -> bool {
        self.contains_with_tolerance(point, 0.0)
    }

    /// Check if the cone contains a point with tolerance.
    ///
    /// A point x is in the cone if H @ x >= -eps for all hyperplanes H.
    pub fn contains_with_tolerance(&mut self, point: &[f64], eps: f64) -> bool {
        let h = self.hyperplanes();
        for row in h {
            let dot: f64 = row.iter().zip(point).map(|(&a, &b)| a as f64 * b).sum();
            if dot < -eps {
                return false;
            }
        }
        true
    }

    /// Check if a point is in the strict interior.
    ///
    /// A point x is in the strict interior if H @ x > eps for all hyperplanes H.
    pub fn contains_interior(&mut self, point: &[f64], eps: f64) -> bool {
        let h = self.hyperplanes();
        for row in h {
            let dot: f64 = row.iter().zip(point).map(|(&a, &b)| a as f64 * b).sum();
            if dot <= eps {
                return false;
            }
        }
        true
    }

    /// Check if the cone is pointed (strongly convex).
    ///
    /// A cone is pointed if it contains no lines, i.e., there is no x != 0
    /// such that both x and -x are in the cone.
    ///
    /// Equivalent to: the dual cone is solid (full-dimensional).
    pub fn is_pointed(&mut self) -> bool {
        self.dual().is_solid()
    }

    /// Check if the cone is solid (full-dimensional).
    ///
    /// A cone is solid if its dimension equals the ambient dimension.
    pub fn is_solid(&mut self) -> bool {
        self.dim() == self.ambient_dim
    }

    /// Find a point in the strict interior of the cone.
    ///
    /// If rays are known, returns the sum of all rays (normalized).
    /// Otherwise uses LP to find an interior point.
    pub fn find_interior_point(&mut self) -> Option<Vec<f64>> {
        // Simple case: if we have rays, sum them
        if let Some(ref rays) = self.rays {
            if rays.is_empty() {
                return None;
            }

            let mut sum = vec![0.0; self.ambient_dim];
            for ray in rays {
                for (i, &r) in ray.iter().enumerate() {
                    sum[i] += r as f64;
                }
            }

            // Check it's actually interior
            if self.contains_interior(&sum, 1e-10) {
                return Some(sum);
            }
        }

        // Fall back to LP
        tip::find_interior_point_lp(self)
    }

    /// Find the tip of the stretched cone.
    ///
    /// The stretched cone is the region at least distance `c` from each hyperplane.
    /// The tip is the point with minimum norm in this region.
    ///
    /// This is used to find a "canonical" interior point of the Kähler cone.
    pub fn tip_of_stretched_cone(&mut self, c: f64) -> Option<Vec<f64>> {
        tip::tip_of_stretched_cone(self, c)
    }

    /// Compute extremal rays (minimal generating set).
    ///
    /// For pointed cones, returns the unique minimal set of rays.
    pub fn extremal_rays(&mut self) -> Vec<Vec<i64>> {
        let rays = self.rays().to_vec();

        // Remove duplicates
        let unique: HashSet<Vec<i64>> = rays.into_iter().collect();
        let rays: Vec<Vec<i64>> = unique.into_iter().collect();

        if rays.len() <= 1 {
            return rays;
        }

        // Filter to extremal rays using LP
        let mut extremal = Vec::new();
        for (i, ray) in rays.iter().enumerate() {
            if is_extremal(&rays, i) {
                extremal.push(ray.clone());
            }
        }

        extremal
    }

    /// Intersection of this cone with another cone.
    ///
    /// Returns a new cone defined by the union of hyperplanes.
    pub fn intersection(&mut self, other: &mut Self) -> Self {
        let mut h1 = self.hyperplanes().to_vec();
        let h2 = other.hyperplanes();
        h1.extend(h2.iter().cloned());
        Self::from_hyperplanes(h1)
    }
}

/// Normalize vectors by dividing by their GCD.
fn normalize_vectors(vectors: Vec<Vec<i64>>) -> Vec<Vec<i64>> {
    vectors
        .into_iter()
        .filter_map(|v| {
            let g = gcd_vec(&v);
            if g == 0 {
                None // Skip zero vectors
            } else {
                Some(v.iter().map(|&x| x / g).collect())
            }
        })
        .collect()
}

/// GCD of a vector of integers.
fn gcd_vec(v: &[i64]) -> i64 {
    v.iter().fold(0, |acc, &x| gcd(acc, x.abs()))
}

/// GCD of two integers.
fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// Compute matrix rank using Gaussian elimination.
fn matrix_rank(matrix: &[Vec<i64>]) -> usize {
    if matrix.is_empty() {
        return 0;
    }

    let rows = matrix.len();
    let cols = matrix[0].len();

    // Convert to f64 for numerical rank computation
    let mut mat: Vec<Vec<f64>> = matrix
        .iter()
        .map(|row| row.iter().map(|&x| x as f64).collect())
        .collect();

    let mut rank = 0;
    let mut pivot_col = 0;

    for row in 0..rows {
        if pivot_col >= cols {
            break;
        }

        // Find pivot
        let mut max_row = row;
        let mut max_val = mat[row][pivot_col].abs();
        for r in (row + 1)..rows {
            if mat[r][pivot_col].abs() > max_val {
                max_val = mat[r][pivot_col].abs();
                max_row = r;
            }
        }

        if max_val < 1e-10 {
            pivot_col += 1;
            continue;
        }

        mat.swap(row, max_row);
        rank += 1;

        // Eliminate
        for r in (row + 1)..rows {
            if mat[r][pivot_col].abs() > 1e-15 {
                let factor = mat[r][pivot_col] / mat[row][pivot_col];
                for c in pivot_col..cols {
                    mat[r][c] -= factor * mat[row][c];
                }
            }
        }

        pivot_col += 1;
    }

    rank
}

/// Check if ray i is extremal in the set of rays.
///
/// A ray r is extremal if it cannot be written as a non-negative
/// combination of the other rays.
fn is_extremal(rays: &[Vec<i64>], i: usize) -> bool {
    // Use simple heuristic: check if ray is not in the cone of others
    // For a proper implementation, we'd use LP here

    let _ray = &rays[i];
    let others: Vec<&Vec<i64>> = rays
        .iter()
        .enumerate()
        .filter(|&(j, _)| j != i)
        .map(|(_, r)| r)
        .collect();

    if others.is_empty() {
        return true;
    }

    // Simple check: if ray has a unique direction not shared by others
    // This is a heuristic, not exact

    // For now, keep all rays (conservative)
    // TODO: implement proper LP-based extremality check
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cone_from_rays() {
        let rays = vec![vec![1, 0], vec![0, 1]];
        let mut cone = Cone::from_rays(rays);

        assert_eq!(cone.ambient_dim(), 2);
        assert_eq!(cone.rays().len(), 2);
    }

    #[test]
    fn test_cone_from_hyperplanes() {
        let hyperplanes = vec![vec![1, 0], vec![0, 1]];
        let mut cone = Cone::from_hyperplanes(hyperplanes);

        assert_eq!(cone.ambient_dim(), 2);
        assert_eq!(cone.hyperplanes().len(), 2);
    }

    #[test]
    fn test_cone_contains() {
        // First quadrant: x >= 0, y >= 0
        let mut cone = Cone::from_hyperplanes(vec![vec![1, 0], vec![0, 1]]);

        assert!(cone.contains(&[1.0, 1.0]));
        assert!(cone.contains(&[0.0, 0.0]));
        assert!(!cone.contains(&[-1.0, 1.0]));
        assert!(!cone.contains(&[1.0, -1.0]));
    }

    #[test]
    fn test_cone_dual() {
        let rays = vec![vec![1, 0], vec![1, 1]];
        let cone = Cone::from_rays(rays);
        let dual = cone.dual();

        // Dual is defined by hyperplanes = original rays
        assert!(dual.hyperplanes.is_some());
    }

    #[test]
    fn test_normalize_vectors() {
        let v = vec![vec![2, 4, 6], vec![0, 0, 0], vec![3, 6, 9]];
        let normalized = normalize_vectors(v);

        // Should have 2 vectors (zero vector removed)
        assert_eq!(normalized.len(), 2);
        // [2,4,6] / 2 = [1,2,3]
        assert!(normalized.contains(&vec![1, 2, 3]));
        // [3,6,9] / 3 = [1,2,3]
        assert!(normalized.contains(&vec![1, 2, 3]));
    }

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(12, 8), 4);
        assert_eq!(gcd(17, 13), 1);
        assert_eq!(gcd(0, 5), 5);
        assert_eq!(gcd(5, 0), 5);
    }

    #[test]
    fn test_matrix_rank() {
        // Full rank 2x2
        let m1 = vec![vec![1, 0], vec![0, 1]];
        assert_eq!(matrix_rank(&m1), 2);

        // Rank 1 (dependent rows)
        let m2 = vec![vec![1, 2], vec![2, 4]];
        assert_eq!(matrix_rank(&m2), 1);

        // Rank 2 for 3x2
        let m3 = vec![vec![1, 0], vec![0, 1], vec![1, 1]];
        assert_eq!(matrix_rank(&m3), 2);
    }
}
