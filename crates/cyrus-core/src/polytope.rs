//! Convex lattice polytope operations.
//!
//! Implements operations on reflexive polytopes for constructing
//! Calabi-Yau hypersurfaces via toric geometry.
//!
//! Reference: [[project_docs/clean_room/DUAL_POLYTOPE.md]]
//! Reference: [[project_docs/clean_room/REFLEXIVE_POLYTOPE_DUAL.md]]

use crate::error::{Error, Result};
use crate::geometry::ConvexHull;
use crate::lattice::Point;
use malachite::Integer;
use std::collections::HashSet;

/// A convex lattice polytope defined by its vertices.
#[derive(Debug, Clone)]
pub struct Polytope {
    vertices: Vec<Point>,
    dim: usize,
}

impl Polytope {
    /// Create a polytope from vertices.
    ///
    /// # Errors
    /// Returns an error if vertices have inconsistent dimensions or are empty.
    pub fn from_vertices(vertices: Vec<Point>) -> Result<Self> {
        if vertices.is_empty() {
            return Err(Error::InvalidInput("no vertices provided".into()));
        }

        let dim = vertices[0].dim();
        for (i, v) in vertices.iter().enumerate() {
            if v.dim() != dim {
                return Err(Error::InvalidDimension {
                    expected: dim,
                    got: v.dim(),
                });
            }
            // Check for duplicates
            for other in &vertices[i + 1..] {
                if v == other {
                    return Err(Error::InvalidInput("duplicate vertex".into()));
                }
            }
        }

        Ok(Self { vertices, dim })
    }

    /// Get the ambient dimension.
    #[inline]
    pub const fn dim(&self) -> usize {
        self.dim
    }

    /// Get the number of vertices.
    #[inline]
    pub const fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Get the vertices.
    #[inline]
    pub fn vertices(&self) -> &[Point] {
        &self.vertices
    }

    /// Check if the origin is strictly interior to the polytope.
    ///
    /// For reflexivity, the origin must be the unique interior point.
    pub fn origin_is_interior(&self) -> bool {
        // Simple check: for each coordinate, vertices span both sides
        for i in 0..self.dim {
            let has_pos = self.vertices.iter().any(|v| v.coords()[i] > 0);
            let has_neg = self.vertices.iter().any(|v| v.coords()[i] < 0);
            if !has_pos || !has_neg {
                return false;
            }
        }
        true
    }

    /// Compute the dual (polar) polytope and return all its lattice points.
    ///
    /// For a reflexive polytope with the origin as interior point,
    /// the dual polytope Δ° is defined by: m · v >= -1 for all vertices v of Δ.
    ///
    /// This returns ALL lattice points in Δ° (vertices + interior/face points),
    /// not just the vertices.
    ///
    /// # Errors
    /// Returns an error if:
    /// - The origin is not interior (required for reflexivity)
    /// - No valid facets are found
    ///
    /// Reference: [[project_docs/clean_room/DUAL_POLYTOPE.md]]
    pub fn compute_dual(&self) -> Result<Self> {
        if !self.origin_is_interior() {
            return Err(Error::NotReflexive(
                "origin must be interior for duality".into(),
            ));
        }

        // First compute the dual vertices (from primal facets)
        let dual_vertices = self.find_dual_vertices()?;

        // Then enumerate all lattice points in the dual polytope
        let all_dual_points = self.enumerate_dual_lattice_points(&dual_vertices)?;

        Self::from_vertices(all_dual_points)
    }

    /// Enumerate all lattice points in the dual polytope.
    ///
    /// The dual Δ° is defined by: m · v >= -1 for all vertices v of Δ.
    /// We use the dual vertices to bound the search space.
    fn enumerate_dual_lattice_points(&self, dual_vertices: &[Point]) -> Result<Vec<Point>> {
        if dual_vertices.is_empty() {
            return Err(Error::InvalidInput("no dual vertices".into()));
        }

        // Compute bounding box from dual vertices
        let dim = self.dim;
        let mut min_coords = vec![i64::MAX; dim];
        let mut max_coords = vec![i64::MIN; dim];

        for dv in dual_vertices {
            for (i, &c) in dv.coords().iter().enumerate() {
                min_coords[i] = min_coords[i].min(c);
                max_coords[i] = max_coords[i].max(c);
            }
        }

        // Get primal vertices for constraint checking
        let primal_verts: Vec<&[i64]> = self.vertices.iter().map(|v| v.coords()).collect();

        // Enumerate all integer points in bounding box and check constraints
        let mut lattice_points = Vec::new();
        let mut candidate = min_coords.clone();

        loop {
            // Check if candidate satisfies all dual constraints: v · m >= -1
            let is_valid = primal_verts.iter().all(|v| {
                let dot: i64 = v.iter().zip(candidate.iter()).map(|(&a, &b)| a * b).sum();
                dot >= -1
            });

            if is_valid {
                lattice_points.push(Point::new(candidate.clone()));
            }

            // Move to next candidate (increment with carry)
            let mut carry = true;
            for i in 0..dim {
                if carry {
                    candidate[i] += 1;
                    if candidate[i] > max_coords[i] {
                        candidate[i] = min_coords[i];
                    } else {
                        carry = false;
                    }
                }
            }

            if carry {
                // We've wrapped around completely, done
                break;
            }
        }

        if lattice_points.is_empty() {
            return Err(Error::InvalidInput("no lattice points found in dual".into()));
        }

        Ok(lattice_points)
    }

    /// Find dual vertices by computing facet normals using convex hull.
    ///
    /// Uses the Beneath-Beyond algorithm for efficient facet enumeration.
    /// For each facet with equation n·x = c, the dual vertex is -n/c
    /// (normalized so m·x = -1 on the facet, m·x >= -1 elsewhere).
    ///
    /// Reference: [[project_docs/clean_room/REFLEXIVE_POLYTOPE_DUAL.md]]
    fn find_dual_vertices(&self) -> Result<Vec<Point>> {
        // Convert vertices to the format expected by ConvexHull
        let points: Vec<Vec<i64>> = self
            .vertices
            .iter()
            .map(|p| p.coords().to_vec())
            .collect();

        // Compute convex hull to get facets
        let hull = ConvexHull::compute(&points)
            .ok_or_else(|| Error::InvalidInput("failed to compute convex hull".into()))?;

        let mut dual_verts: HashSet<Vec<i64>> = HashSet::new();

        for facet in &hull.facets {
            // Facet equation: n·x = c
            // For reflexive polytope with origin interior, c ≠ 0
            // Dual vertex: m = -n/c (so m·x = -1 on facet)

            let c = &facet.constant;
            if *c == Integer::from(0) {
                // Facet passes through origin - not reflexive
                return Err(Error::NotReflexive(
                    "facet passes through origin".into(),
                ));
            }

            // Compute m = -n/c and check integrality
            let mut dual_vertex = Vec::with_capacity(self.dim);
            let mut is_integral = true;

            for n_i in &facet.normal {
                // m_i = -n_i / c
                let neg_n_i = -n_i;

                // Check if division is exact
                if &neg_n_i % c != Integer::from(0) {
                    is_integral = false;
                    break;
                }

                let m_i = &neg_n_i / c;

                // Convert to i64
                match i64::try_from(&m_i) {
                    Ok(val) => dual_vertex.push(val),
                    Err(_) => {
                        is_integral = false;
                        break;
                    }
                }
            }

            if !is_integral {
                return Err(Error::NotReflexive(
                    "dual vertex is not integral".into(),
                ));
            }

            dual_verts.insert(dual_vertex);
        }

        let vertices: Vec<Point> = dual_verts.into_iter().map(Point::new).collect();

        if vertices.is_empty() {
            return Err(Error::InvalidInput("no facets found".into()));
        }

        Ok(vertices)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square() -> Polytope {
        Polytope::from_vertices(vec![
            Point::new(vec![-1, -1]),
            Point::new(vec![1, -1]),
            Point::new(vec![1, 1]),
            Point::new(vec![-1, 1]),
        ])
        .unwrap()
    }

    #[test]
    fn test_polytope_creation() {
        let p = square();
        assert_eq!(p.dim(), 2);
        assert_eq!(p.num_vertices(), 4);
    }

    #[test]
    fn test_origin_interior() {
        let p = square();
        assert!(p.origin_is_interior());
    }

    #[test]
    fn test_origin_not_interior() {
        // Triangle not containing origin
        let p = Polytope::from_vertices(vec![
            Point::new(vec![1, 0]),
            Point::new(vec![2, 0]),
            Point::new(vec![1, 1]),
        ])
        .unwrap();
        assert!(!p.origin_is_interior());
    }

    #[test]
    fn test_empty_vertices_error() {
        let result = Polytope::from_vertices(vec![]);
        assert!(result.is_err());
    }

    #[test]
    fn test_dimension_mismatch_error() {
        let result =
            Polytope::from_vertices(vec![Point::new(vec![1, 2]), Point::new(vec![1, 2, 3])]);
        assert!(result.is_err());
    }

    #[test]
    fn test_compute_dual_square() {
        // The dual of the unit square [-1,1]^2 with vertices (±1, ±1)
        // is the diamond with lattice points: origin + vertices (±1, 0), (0, ±1)
        let p = square();
        let dual = p.compute_dual().unwrap();

        assert_eq!(dual.dim(), 2);

        let dual_coords: HashSet<Vec<i64>> = dual
            .vertices()
            .iter()
            .map(|v| v.coords().to_vec())
            .collect();

        // Expected: all lattice points in dual (origin + 4 vertices)
        let expected: HashSet<Vec<i64>> = [
            vec![0, 0], // origin (interior)
            vec![1, 0],
            vec![-1, 0],
            vec![0, 1],
            vec![0, -1],
        ]
        .into_iter()
        .collect();

        assert_eq!(dual_coords, expected);
    }

    #[test]
    fn test_compute_dual_not_reflexive() {
        // Triangle not containing origin - should fail
        let p = Polytope::from_vertices(vec![
            Point::new(vec![1, 0]),
            Point::new(vec![2, 0]),
            Point::new(vec![1, 1]),
        ])
        .unwrap();
        assert!(p.compute_dual().is_err());
    }

    #[test]
    fn test_compute_dual_3d_cube() {
        // Unit cube [-1,1]^3 with 8 vertices
        let p = Polytope::from_vertices(vec![
            Point::new(vec![-1, -1, -1]),
            Point::new(vec![1, -1, -1]),
            Point::new(vec![-1, 1, -1]),
            Point::new(vec![1, 1, -1]),
            Point::new(vec![-1, -1, 1]),
            Point::new(vec![1, -1, 1]),
            Point::new(vec![-1, 1, 1]),
            Point::new(vec![1, 1, 1]),
        ])
        .unwrap();

        let dual = p.compute_dual().unwrap();

        // Dual of cube is octahedron with lattice points: origin + 6 vertices
        assert_eq!(dual.dim(), 3);
        assert_eq!(dual.num_vertices(), 7); // 6 vertices + 1 origin

        let dual_coords: HashSet<Vec<i64>> = dual
            .vertices()
            .iter()
            .map(|v| v.coords().to_vec())
            .collect();

        let expected: HashSet<Vec<i64>> = [
            vec![0, 0, 0], // origin (interior)
            vec![1, 0, 0],
            vec![-1, 0, 0],
            vec![0, 1, 0],
            vec![0, -1, 0],
            vec![0, 0, 1],
            vec![0, 0, -1],
        ]
        .into_iter()
        .collect();

        assert_eq!(dual_coords, expected);
    }
}
