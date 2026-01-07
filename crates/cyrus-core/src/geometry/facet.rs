//! Facet representation for convex hulls.
//!
//! A facet is a (d-1)-dimensional face of a d-dimensional polytope.

use malachite::Integer;

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
    pub(crate) fn ridges(&self) -> Vec<Vec<usize>> {
        let d = self.vertices.len();
        let mut ridges = Vec::with_capacity(d);
        for i in 0..d {
            let mut ridge: Vec<usize> = self.vertices.clone();
            ridge.remove(i);
            ridge.sort_unstable();
            ridges.push(ridge);
        }
        ridges
    }

    /// Check if a point is visible from this facet (strictly above).
    pub(crate) fn is_visible(&self, point: &[i64]) -> bool {
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
    pub(crate) fn is_on(&self, point: &[i64]) -> bool {
        let dot: Integer = self
            .normal
            .iter()
            .zip(point.iter())
            .map(|(n, &p)| n * Integer::from(p))
            .sum();
        dot == self.constant
    }
}
