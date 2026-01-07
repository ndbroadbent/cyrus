//! Facet-related operations for polytopes.
//!
//! Operations for determining which points lie interior to facets,
//! used for triangulation point filtering.

use crate::error::Result;
use crate::lattice::Point;

use super::Polytope;

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
