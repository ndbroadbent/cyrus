//! Core polytope struct and basic operations.

use crate::error::{Error, Result};
use crate::lattice::Point;

/// A convex lattice polytope defined by its vertices.
#[derive(Debug, Clone)]
pub struct Polytope {
    pub(super) vertices: Vec<Point>,
    pub(super) dim: usize,
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
