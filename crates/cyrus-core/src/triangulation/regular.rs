//! Regular triangulation computation using the lifting map.
//!
//! Given points S and heights H, the triangulation is the projection of
//! the lower convex hull of the lifted points (p, h_p).

use super::Triangulation;
use super::convex_hull::{LiftedPoints, convex_hull, is_lower_face};
use crate::Point;
use crate::error::{Error, Result};
use malachite::Integer;

/// Scale factor for converting f64 heights to integers.
/// 2^40 ~ 10^12 gives plenty of precision while staying in reasonable integer range.
const HEIGHT_SCALE: f64 = 1_099_511_627_776.0; // 2^40

/// Compute a regular triangulation of a set of points using the lifting map.
///
/// Given points S and heights H, the triangulation is the projection of
/// the lower convex hull of the lifted points (p, h_p).
///
/// Uses scaled Integer arithmetic for exact orientation tests.
///
/// # Arguments
/// * `points` - Lattice points to triangulate.
/// * `heights` - Heights used for the lifting map.
///
/// # Errors
/// Returns an error if the number of points and heights do not match.
///
/// Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 4.1
pub fn compute_regular_triangulation(points: &[Point], heights: &[f64]) -> Result<Triangulation> {
    if points.len() != heights.len() {
        return Err(Error::InvalidInput(
            "Number of points and heights must match".into(),
        ));
    }

    // 1. Lift points: (v_i, h_i) in Z^{d+1} using scaled Integer for exactness
    let lifted: Vec<Vec<Integer>> = points
        .iter()
        .zip(heights.iter())
        .map(|(p, &h)| {
            let mut v: Vec<Integer> = p.coords().iter().map(|&x| Integer::from(x)).collect();
            // Scale height to integer: round(h * 2^40)
            #[allow(clippy::cast_possible_truncation)]
            let h_scaled = (h * HEIGHT_SCALE).round() as i64;
            v.push(Integer::from(h_scaled));
            v
        })
        .collect();

    // 2. Wrap in LiftedPoints
    let Some(lifted) = LiftedPoints::new(lifted) else {
        return Err(Error::InvalidInput(
            "regular triangulation requires a non-empty point set".into(),
        ));
    };

    // 3. Compute Convex Hull in d+1 dimensions
    let facets = convex_hull(&lifted);

    // 4. Extract lower faces
    let mut simplices = Vec::new();
    for facet in facets {
        if is_lower_face(&facet, &lifted) {
            simplices.push(facet);
        }
    }

    // An empty subdivision is never a valid answer for a genuine input;
    // failing loudly here beats every downstream stage silently agreeing
    // that nothing intersects anything.
    if simplices.is_empty() {
        return Err(Error::InvalidInput(
            "regular triangulation produced no lower faces (degenerate lift)".into(),
        ));
    }

    Ok(Triangulation::new(simplices))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_regular_triangulation_spiky_heights_smoke_repro() {
        // Landscape-smoke regression: CYTools default heights for a KS h11=4
        // polytope; the lower-hull extraction silently returned no simplices.
        let points = vec![
            Point::new(vec![0, 0, 0, 0]),
            Point::new(vec![-5, -1, -3, -2]),
            Point::new(vec![0, 0, 1, 0]),
            Point::new(vec![1, 0, 0, 0]),
            Point::new(vec![1, 1, 1, 2]),
            Point::new(vec![-2, 0, -1, 0]),
            Point::new(vec![0, -1, 0, 0]),
            Point::new(vec![0, 1, 0, 0]),
            Point::new(vec![-1, 0, 0, 0]),
        ];
        let heights = vec![0.0, 39.0, 1.0, 1.0, 7.0, 5.0, 1.0, 1.0, 1.0];
        let tri = compute_regular_triangulation(&points, &heights).unwrap();
        assert_eq!(tri.simplices().len(), 16);
    }

    #[test]
    fn test_regular_triangulation_2d_square() {
        // Square (0,0), (1,0), (1,1), (0,1)
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![1, 1]),
            Point::new(vec![0, 1]),
        ];
        // heights: 0, 2, 6, 3
        let heights = vec![0.0, 2.0, 6.0, 3.0];

        let tri = compute_regular_triangulation(&points, &heights).unwrap();
        // Should produce 2 triangles
        assert_eq!(tri.simplices().len(), 2);
    }
}
