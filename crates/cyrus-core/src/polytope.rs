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

    /// Get points that are NOT interior to facets.
    ///
    /// A point is "interior to a facet" if it lies on exactly ONE facet.
    /// This method returns points that lie on 0 (interior) or 2+ (edges/vertices) facets.
    ///
    /// These are the points used for triangulation - points interior to facets
    /// don't affect the triangulation structure.
    ///
    /// For a reflexive polytope, we use the dual vertices to determine facets.
    /// Each dual vertex m corresponds to a primal facet with equation p · m = -1.
    /// A primal point p "saturates" this facet if p · m = -1.
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
                // A primal point p saturates facet m if p · m = -1
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

    #[test]
    fn test_points_not_interior_to_facets_square() {
        // Square with 4 vertices - all are on 2 facets (edges)
        // No points interior to facets
        let p = square();
        let filtered = p.points_not_interior_to_facets().unwrap();

        // All 4 vertices should be kept (each on 2 edges)
        assert_eq!(filtered.len(), 4);
    }

    #[test]
    fn test_points_not_interior_to_facets_with_edge_point() {
        // Square with an extra point on an edge
        // The edge point saturates exactly 1 facet, so it should be excluded
        let p = Polytope::from_vertices(vec![
            Point::new(vec![-1, -1]),
            Point::new(vec![1, -1]),
            Point::new(vec![1, 1]),
            Point::new(vec![-1, 1]),
            Point::new(vec![0, -1]), // Point on bottom edge - saturates 1 facet
        ])
        .unwrap();

        let filtered = p.points_not_interior_to_facets().unwrap();

        // Should have 4 vertices, edge point excluded
        assert_eq!(filtered.len(), 4);

        // Verify edge point is not in result
        let coords: Vec<Vec<i64>> = filtered.iter().map(|p| p.coords().to_vec()).collect();
        assert!(!coords.contains(&vec![0, -1]));
    }

    #[test]
    fn test_points_not_interior_to_facets_with_interior() {
        // Square with origin inside
        let p = Polytope::from_vertices(vec![
            Point::new(vec![-1, -1]),
            Point::new(vec![1, -1]),
            Point::new(vec![1, 1]),
            Point::new(vec![-1, 1]),
            Point::new(vec![0, 0]), // Origin - interior, saturates 0 facets
        ])
        .unwrap();

        let filtered = p.points_not_interior_to_facets().unwrap();

        // Should have 5 points: 4 vertices + 1 interior
        assert_eq!(filtered.len(), 5);

        // Interior point should be FIRST (CYTools ordering)
        assert_eq!(filtered[0].coords(), &[0, 0]);
    }

    #[test]
    fn test_points_not_interior_to_facets_ordering() {
        // Test that ordering is correct:
        // 1. Interior points first
        // 2. Then by decreasing saturation count
        let p = Polytope::from_vertices(vec![
            Point::new(vec![-1, -1]),
            Point::new(vec![1, -1]),
            Point::new(vec![1, 1]),
            Point::new(vec![-1, 1]),
            Point::new(vec![0, 0]), // Interior (0 saturations)
        ])
        .unwrap();

        let filtered = p.points_not_interior_to_facets().unwrap();

        // First point should be interior (0 saturations)
        assert_eq!(filtered[0].coords(), &[0, 0]);

        // Remaining points should all have saturation count >= 2
        // (they're on 2 edges each for a square)
    }

    #[test]
    fn test_points_not_interior_to_facets_mcallister() {
        // Load McAllister polytope 4-214-647
        // Expected: 294 total points -> 219 after filtering (75 removed)
        let manifest_dir = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let input_path = manifest_dir.join("tests/fixtures/4_214_647/polytope.json");

        #[derive(serde::Deserialize)]
        struct PolytopeInput {
            primal_points: Vec<Vec<i64>>,
        }

        let content = std::fs::read_to_string(&input_path)
            .unwrap_or_else(|e| panic!("Failed to read {}: {e}", input_path.display()));
        let input: PolytopeInput = serde_json::from_str(&content)
            .unwrap_or_else(|e| panic!("Failed to parse {}: {e}", input_path.display()));

        let points: Vec<Point> = input
            .primal_points
            .iter()
            .map(|coords| Point::new(coords.clone()))
            .collect();

        assert_eq!(points.len(), 294, "McAllister polytope should have 294 points");

        let polytope = Polytope::from_vertices(points).unwrap();
        let filtered = polytope.points_not_interior_to_facets().unwrap();

        // McAllister uses 219 points for triangulation (75 are interior to facets)
        assert_eq!(
            filtered.len(),
            219,
            "Expected 219 points not interior to facets, got {}",
            filtered.len()
        );

        // Origin should be first (interior point with 0 saturations)
        assert!(
            filtered[0].coords().iter().all(|&x| x == 0),
            "First point should be the origin"
        );
    }
}
