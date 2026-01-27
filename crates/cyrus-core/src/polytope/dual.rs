//! Dual (polar) polytope computation.
//!
//! For a reflexive polytope with the origin as interior point,
//! the dual polytope is defined by: m . v >= -1 for all vertices v.
//!
//! Reference: [[project_docs/clean_room/DUAL_POLYTOPE.md]]
//! Reference: [[project_docs/clean_room/REFLEXIVE_POLYTOPE_DUAL.md]]

use crate::error::{Error, Result};
use crate::geometry::ConvexHull;
use crate::lattice::Point;
use malachite::Integer;
use std::collections::HashSet;

use super::Polytope;

impl Polytope {
    /// Return the dual vertices (facet normals) of the polytope.
    ///
    /// These are the vertices of the dual polytope (not all dual lattice points).
    pub fn dual_vertices(&self) -> Result<Vec<Point>> {
        self.find_dual_vertices()
    }
    /// Compute the dual (polar) polytope and return all its lattice points.
    ///
    /// For a reflexive polytope with the origin as interior point,
    /// the dual polytope delta_dual is defined by: m . v >= -1 for all vertices v of delta.
    ///
    /// This returns ALL lattice points in delta_dual (vertices + interior/face points),
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
    /// The dual delta_dual is defined by: m . v >= -1 for all vertices v of delta.
    /// We use the dual vertices to bound the search space.
    pub(super) fn enumerate_dual_lattice_points(
        &self,
        dual_vertices: &[Point],
    ) -> Result<Vec<Point>> {
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
        let primal_verts: Vec<&[i64]> = self
            .vertices
            .iter()
            .map(crate::lattice::Point::coords)
            .collect();

        // Enumerate all integer points in bounding box and check constraints
        let mut lattice_points = Vec::new();
        let mut candidate = min_coords.clone();

        loop {
            // Check if candidate satisfies all dual constraints: v . m >= -1
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
            return Err(Error::InvalidInput(
                "no lattice points found in dual".into(),
            ));
        }

        Ok(lattice_points)
    }

    /// Find dual vertices by computing facet normals using convex hull.
    ///
    /// Uses the Beneath-Beyond algorithm for efficient facet enumeration.
    /// For each facet with equation n.x = c, the dual vertex is -n/c
    /// (normalized so m.x = -1 on the facet, m.x >= -1 elsewhere).
    ///
    /// Reference: [[project_docs/clean_room/REFLEXIVE_POLYTOPE_DUAL.md]]
    pub(super) fn find_dual_vertices(&self) -> Result<Vec<Point>> {
        // Convert vertices to the format expected by ConvexHull
        let points: Vec<Vec<i64>> = self.vertices.iter().map(|p| p.coords().to_vec()).collect();

        // Compute convex hull to get facets
        let hull = ConvexHull::compute(&points)
            .ok_or_else(|| Error::InvalidInput("failed to compute convex hull".into()))?;

        let mut dual_verts: HashSet<Vec<i64>> = HashSet::new();

        for facet in &hull.facets {
            // Facet equation: n.x = c
            // For reflexive polytope with origin interior, c != 0
            // Dual vertex: m = -n/c (so m.x = -1 on facet)

            let c = &facet.constant;
            if *c == Integer::from(0) {
                // Facet passes through origin - not reflexive
                return Err(Error::NotReflexive("facet passes through origin".into()));
            }

            // Compute m = -n/c and check integrality
            let mut dual_vertex = Vec::with_capacity(self.dim);
            let mut is_integral = true;

            for n_i in &facet.normal {
                // m_i = -n_i / c
                let neg_n_i = -n_i;

                // Check if division is exact
                if &neg_n_i % c != 0 {
                    is_integral = false;
                    break;
                }

                let m_i = &neg_n_i / c;

                // Convert to i64
                if let Ok(val) = i64::try_from(&m_i) {
                    dual_vertex.push(val);
                } else {
                    is_integral = false;
                    break;
                }
            }

            if !is_integral {
                return Err(Error::NotReflexive("dual vertex is not integral".into()));
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
