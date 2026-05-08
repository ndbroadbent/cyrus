//! Facet-related operations for polytopes.
//!
//! Operations for determining which points lie interior to facets,
//! used for triangulation point filtering.

use std::collections::HashSet;

use crate::error::{Error, Result};
use crate::lattice::Point;

use super::Polytope;

/// CYTools-style face index sets for a four-dimensional reflexive polytope.
///
/// The indices refer to the caller-supplied point list, not necessarily to the
/// polytope's stored lattice-point order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PolytopeFaces4d {
    /// Facets as point-index sets.
    pub facets: Vec<Vec<usize>>,
    /// Facet boundary point-index sets, excluding points interior to exactly
    /// one 4D facet. CYTools uses `labels_bdry` for star-constraint circuits.
    pub facet_boundaries: Vec<Vec<usize>>,
    /// Two-dimensional faces as point-index sets.
    pub twofaces: Vec<Vec<usize>>,
}

impl Polytope {
    /// Get points that are NOT interior to facets.
    ///
    /// A point is "interior to a facet" if it lies on exactly ONE facet.
    /// This method returns points that lie on 0 (interior) or 2+ (edges/vertices) facets.
    ///
    /// These are the points used for triangulation - points interior to facets
    /// don't affect the triangulation structure.
    ///
    /// For a reflexive polytope, we use the dual vertices to determine facets.
    /// Each dual vertex m corresponds to a primal facet with equation p . m = -1.
    /// A primal point p "saturates" this facet if p . m = -1.
    ///
    /// Returns points in CYTools order:
    /// 1. Interior points (0 facets) first
    /// 2. Boundary points by decreasing saturation count
    /// 3. Lexicographic within each group
    pub fn points_not_interior_to_facets(&self) -> Result<Vec<Point>> {
        // Get dual vertices - each corresponds to a primal facet
        let dual_vertices = self.find_dual_vertices()?;

        // For each point, count how many facets it saturates
        let mut point_saturations: Vec<(usize, usize, &Point)> = Vec::new();

        for (idx, point) in self.vertices.iter().enumerate() {
            let coords = point.coords();
            let mut saturation_count = 0;

            for dual_v in &dual_vertices {
                // A primal point p saturates facet m if p . m = -1
                let dot: i64 = coords
                    .iter()
                    .zip(dual_v.coords().iter())
                    .map(|(&p, &m)| p * m)
                    .sum();

                if dot == -1 {
                    saturation_count += 1;
                }
            }

            // Keep points with saturation_count != 1
            if saturation_count != 1 {
                point_saturations.push((saturation_count, idx, point));
            }
        }

        // Sort by CYTools ordering:
        // 1. Interior points first (saturation_count == 0)
        // 2. Then by DECREASING saturation count
        // 3. Lexicographic within groups
        point_saturations.sort_by(|a, b| {
            // First compare by saturation count category
            let a_interior = a.0 == 0;
            let b_interior = b.0 == 0;

            if a_interior != b_interior {
                // Interior points come first
                return b_interior.cmp(&a_interior);
            }

            if !a_interior {
                // Both are boundary - sort by DECREASING saturation count
                match b.0.cmp(&a.0) {
                    std::cmp::Ordering::Equal => {}
                    other => return other,
                }
            }

            // Within same saturation count, sort lexicographically by coordinates
            a.2.coords().cmp(b.2.coords())
        });

        Ok(point_saturations
            .into_iter()
            .map(|(_, _, p)| p.clone())
            .collect())
    }

    /// Compute CYTools-style facet and two-face point-index sets in 4D.
    ///
    /// This ports the face construction used by CYTools' 4D Mori-cap and
    /// toric-curve routines. Facets are selected by dual-vertex equations
    /// `p.m = -1`; two-faces are intersections of pairs of facets whose
    /// polytope vertices share at least three points.
    ///
    /// # Errors
    ///
    /// Returns an error if the polytope is not four-dimensional, if the dual
    /// vertices cannot be computed, or if supplied points have the wrong
    /// ambient dimension.
    pub fn faces_4d_for_points(&self, points: &[Point]) -> Result<PolytopeFaces4d> {
        if self.dim() != 4 {
            return Err(Error::InvalidInput(
                "faces_4d only defined for 4D polytopes".into(),
            ));
        }
        if points.iter().any(|point| point.dim() != self.dim()) {
            return Err(Error::InvalidInput(
                "faces_4d point dimensions do not match polytope dimension".into(),
            ));
        }

        let dual_vertices = self.dual_vertices()?;
        if dual_vertices.is_empty() {
            return Err(Error::InvalidInput("no dual vertices found".into()));
        }

        let poly_vertices = self.vertices();
        let point_saturation_counts = points
            .iter()
            .map(|point| {
                dual_vertices
                    .iter()
                    .filter(|dual_vertex| lattice_dot(point.coords(), dual_vertex.coords()) == -1)
                    .count()
            })
            .collect::<Vec<_>>();
        let mut facet_vertex_sets = Vec::with_capacity(dual_vertices.len());
        let mut facets = Vec::with_capacity(dual_vertices.len());
        let mut facet_boundaries = Vec::with_capacity(dual_vertices.len());
        for dual_vertex in &dual_vertices {
            let mut vertex_set = HashSet::new();
            for (idx, vertex) in poly_vertices.iter().enumerate() {
                if lattice_dot(vertex.coords(), dual_vertex.coords()) == -1 {
                    vertex_set.insert(idx);
                }
            }
            facet_vertex_sets.push(vertex_set);

            let mut facet_points = Vec::new();
            for (idx, point) in points.iter().enumerate() {
                if lattice_dot(point.coords(), dual_vertex.coords()) == -1 {
                    facet_points.push(idx);
                }
            }
            facet_points.sort_unstable();
            let mut facet_boundary = facet_points
                .iter()
                .copied()
                .filter(|&point_idx| point_saturation_counts[point_idx] > 1)
                .collect::<Vec<_>>();
            facet_boundary.sort_unstable();
            facets.push(facet_points);
            facet_boundaries.push(facet_boundary);
        }

        let mut twofaces = Vec::new();
        for (left_idx, left_vertex_set) in facet_vertex_sets.iter().enumerate() {
            for (right_vertex_set, right_dual_vertex) in facet_vertex_sets
                .iter()
                .zip(dual_vertices.iter())
                .skip(left_idx + 1)
            {
                let shared_vertex_count = left_vertex_set.intersection(right_vertex_set).count();
                if shared_vertex_count < 3 {
                    continue;
                }

                let mut face_points = Vec::new();
                for (idx, point) in points.iter().enumerate() {
                    if lattice_dot(point.coords(), dual_vertices[left_idx].coords()) == -1
                        && lattice_dot(point.coords(), right_dual_vertex.coords()) == -1
                    {
                        face_points.push(idx);
                    }
                }
                face_points.sort_unstable();
                twofaces.push(face_points);
            }
        }

        Ok(PolytopeFaces4d {
            facets,
            facet_boundaries,
            twofaces,
        })
    }

    /// Debug version that also returns saturation histogram.
    #[cfg(test)]
    pub fn points_not_interior_to_facets_debug(
        &self,
    ) -> Result<(Vec<Point>, std::collections::HashMap<usize, usize>)> {
        use std::collections::HashMap;

        let dual_vertices = self.find_dual_vertices()?;

        let mut histogram: HashMap<usize, usize> = HashMap::new();
        let mut point_saturations: Vec<(usize, usize, &Point)> = Vec::new();

        for (idx, point) in self.vertices.iter().enumerate() {
            let coords = point.coords();
            let mut saturation_count = 0;

            for dual_v in &dual_vertices {
                let dot: i64 = coords
                    .iter()
                    .zip(dual_v.coords().iter())
                    .map(|(&p, &m)| p * m)
                    .sum();

                if dot == -1 {
                    saturation_count += 1;
                }
            }

            *histogram.entry(saturation_count).or_insert(0) += 1;

            if saturation_count != 1 {
                point_saturations.push((saturation_count, idx, point));
            }
        }

        point_saturations.sort_by(|a, b| {
            let a_interior = a.0 == 0;
            let b_interior = b.0 == 0;
            if a_interior != b_interior {
                return b_interior.cmp(&a_interior);
            }
            if !a_interior {
                match b.0.cmp(&a.0) {
                    std::cmp::Ordering::Equal => {}
                    other => return other,
                }
            }
            a.2.coords().cmp(b.2.coords())
        });

        Ok((
            point_saturations
                .into_iter()
                .map(|(_, _, p)| p.clone())
                .collect(),
            histogram,
        ))
    }
}

fn lattice_dot(left: &[i64], right: &[i64]) -> i64 {
    left.iter().zip(right.iter()).map(|(&a, &b)| a * b).sum()
}
