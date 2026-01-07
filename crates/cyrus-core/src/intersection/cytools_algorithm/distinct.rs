//! Distinct intersection numbers from simplices.

use crate::Point;
use crate::integer_math::determinant_gaussian;
use crate::triangulation::Triangulation;
use malachite::Rational;
use std::collections::HashSet;

/// Distinct intersection numbers from simplices.
/// Each entry: (i, j, k, l, 1/|det|) for non-origin indices in sorted order.
#[derive(Debug, Clone)]
pub struct DistinctIntnum {
    /// Non-origin indices, sorted
    pub indices: [usize; 4],
    /// 1/|det(pts_ext)|
    pub inv_det: f64,
}

/// Compute distinct intersection numbers from all simplices.
pub fn compute_distinct_intnums(tri: &Triangulation, points: &[Point]) -> Vec<DistinctIntnum> {
    let mut result = Vec::new();
    let mut seen: HashSet<[usize; 4]> = HashSet::new();

    for simplex in tri.simplices() {
        // Get non-origin indices from simplex (indices != 0)
        let non_origin: Vec<usize> = simplex.iter().copied().filter(|&i| i != 0).collect();

        if non_origin.len() < 4 {
            continue;
        }

        // For 4D, simplex has 5 vertices, so after removing origin we have 4 non-origin vertices
        // Build pts_ext: augment with 1s column
        let mut pts_ext: Vec<Vec<Rational>> = simplex
            .iter()
            .map(|&idx| {
                let mut row: Vec<Rational> = points[idx]
                    .coords()
                    .iter()
                    .map(|&x| Rational::from(x))
                    .collect();
                row.push(Rational::from(1));
                row
            })
            .collect();

        let det = determinant_gaussian(&mut pts_ext);
        let abs_det: f64 = if det < 0 {
            (-det).to_string().parse().unwrap_or(1.0)
        } else {
            det.to_string().parse().unwrap_or(1.0)
        };

        if abs_det < 1e-10 {
            continue;
        }

        let inv_det = 1.0 / abs_det;

        // Sort non-origin indices
        let mut indices: [usize; 4] = [non_origin[0], non_origin[1], non_origin[2], non_origin[3]];
        indices.sort_unstable();

        if seen.insert(indices) {
            result.push(DistinctIntnum { indices, inv_det });
        }
    }

    // Sort by indices for determinism
    result.sort_by(|a, b| a.indices.cmp(&b.indices));
    result
}
