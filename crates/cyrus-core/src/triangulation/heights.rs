//! Height computation for regular triangulations.
//!
//! Heights are used in the lifting map to induce triangulations.

use super::Triangulation;
use super::regular::compute_regular_triangulation;
use crate::Point;
use crate::error::{Error, Result};

/// Compute Delaunay heights for a set of points.
///
/// For each point p, the height is h(p) = |p|^2 = p . p (squared Euclidean norm).
/// This is the standard choice for inducing a Delaunay-like triangulation.
pub fn compute_delaunay_heights(points: &[Point]) -> Vec<f64> {
    points
        .iter()
        .map(|p| p.coords().iter().map(|&x| (x as f64) * (x as f64)).sum())
        .collect()
}

/// Compute default FRST (Fine Regular Star Triangulation) heights.
///
/// Starting from Delaunay heights (|p|^2), adjusts the origin's height
/// downward until the triangulation has the Star property (origin in every simplex).
///
/// # Arguments
/// * `points` - The lattice points to triangulate
/// * `origin_idx` - Index of the origin point (must be in points)
///
/// # Returns
/// * `Ok((heights, triangulation))` - The adjusted heights and resulting triangulation
/// * `Err` if origin_idx is out of bounds or triangulation fails
///
/// Reference: CYTools algorithm (Demirtas et al. 2022, Sec. 3.2)
pub fn compute_frst_heights(
    points: &[Point],
    origin_idx: usize,
) -> Result<(Vec<f64>, Triangulation)> {
    if origin_idx >= points.len() {
        return Err(Error::InvalidInput(format!(
            "origin_idx {} out of bounds for {} points",
            origin_idx,
            points.len()
        )));
    }

    // Start with Delaunay heights
    let mut heights = compute_delaunay_heights(points);

    // Maximum iterations to prevent infinite loop
    const MAX_ITERATIONS: usize = 100;

    for _ in 0..MAX_ITERATIONS {
        let tri = compute_regular_triangulation(points, &heights)?;

        if tri.is_star(origin_idx) {
            return Ok((heights, tri));
        }

        // Adjust origin height downward
        let min_h = heights.iter().copied().fold(f64::INFINITY, f64::min);
        let max_h = heights.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        // h(origin) -= (max - min + 10)
        heights[origin_idx] -= (max_h - min_h) + 10.0;
    }

    Err(Error::InvalidInput(
        "Failed to achieve Star triangulation after max iterations".into(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_compute_delaunay_heights() {
        let points = vec![
            Point::new(vec![0, 0]),   // |p|^2 = 0
            Point::new(vec![1, 0]),   // |p|^2 = 1
            Point::new(vec![1, 1]),   // |p|^2 = 2
            Point::new(vec![3, 4]),   // |p|^2 = 25
            Point::new(vec![-2, -2]), // |p|^2 = 8
        ];

        let heights = compute_delaunay_heights(&points);

        assert_eq!(heights.len(), 5);
        assert!((heights[0] - 0.0).abs() < 1e-10);
        assert!((heights[1] - 1.0).abs() < 1e-10);
        assert!((heights[2] - 2.0).abs() < 1e-10);
        assert!((heights[3] - 25.0).abs() < 1e-10);
        assert!((heights[4] - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_frst_heights_simple() {
        let points = vec![
            Point::new(vec![0, 0]), // Origin at index 0
            Point::new(vec![-1, -1]),
            Point::new(vec![1, -1]),
            Point::new(vec![1, 1]),
            Point::new(vec![-1, 1]),
        ];

        let result = compute_frst_heights(&points, 0);
        assert!(result.is_ok());

        let (heights, tri) = result.unwrap();
        assert_eq!(heights.len(), 5);
        assert!(tri.is_star(0));
    }

    #[test]
    fn test_compute_frst_heights_invalid_origin() {
        let points = vec![Point::new(vec![0, 0]), Point::new(vec![1, 0])];

        let result = compute_frst_heights(&points, 10);
        assert!(result.is_err());
    }
}
