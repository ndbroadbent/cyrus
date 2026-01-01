//! Convex lattice polytope operations.
//!
//! Implements operations on reflexive polytopes for constructing
//! Calabi-Yau hypersurfaces via toric geometry.

use crate::error::{Error, Result};
use crate::lattice::Point;

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
}
