//! Hyperplane representation and operations.
//!
//! A hyperplane in d dimensions is defined by n · x = c where n is the normal
//! vector and c is the constant term.

use malachite::Integer;

use super::orientation::{Orientation, determinant_exact};

/// A hyperplane in d dimensions: n · x = c
#[derive(Debug, Clone)]
pub struct Hyperplane {
    /// Normal vector (as integers).
    pub normal: Vec<Integer>,
    /// Constant term.
    pub constant: Integer,
}

impl Hyperplane {
    /// Construct hyperplane through d points in d dimensions.
    ///
    /// Uses cofactor expansion to compute the normal vector.
    pub fn from_points(points: &[&[i64]]) -> Option<Self> {
        let d = points.len();
        if d == 0 {
            return None;
        }
        for p in points {
            if p.len() != d {
                return None;
            }
        }

        // Build (d-1) difference vectors: v_i - v_0
        let v0 = points[0];
        let diffs: Vec<Vec<Integer>> = points[1..]
            .iter()
            .map(|v| {
                v.iter()
                    .zip(v0.iter())
                    .map(|(&a, &b)| Integer::from(a - b))
                    .collect()
            })
            .collect();

        // Compute normal via cofactors
        // n_i = (-1)^i * det(M_i) where M_i is diffs with column i removed
        let mut normal = Vec::with_capacity(d);
        for i in 0..d {
            let minor = extract_minor(&diffs, i);
            let det = determinant_exact(&minor);
            let sign = if i % 2 == 0 {
                Integer::from(1)
            } else {
                Integer::from(-1)
            };
            normal.push(sign * det);
        }

        // Check for degenerate case (all zeros)
        if normal.iter().all(|x| *x == 0) {
            return None;
        }

        // Compute constant: c = n · v0
        let constant: Integer = normal
            .iter()
            .zip(v0.iter())
            .map(|(n, &v)| n * Integer::from(v))
            .sum();

        Some(Self { normal, constant })
    }

    /// Check which side of the hyperplane a point is on.
    pub fn side(&self, point: &[i64]) -> Orientation {
        let dot: Integer = self
            .normal
            .iter()
            .zip(point.iter())
            .map(|(n, &p)| n * Integer::from(p))
            .sum();

        Orientation::from_sign(&(&dot - &self.constant))
    }
}

/// Extract minor matrix by removing column `col` from a (d-1)×d matrix.
fn extract_minor(matrix: &[Vec<Integer>], col: usize) -> Vec<Vec<Integer>> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .filter(|(j, _)| *j != col)
                .map(|(_, v)| v.clone())
                .collect()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hyperplane_from_points_2d() {
        // Line through (0,0) and (1,0)
        let hp = Hyperplane::from_points(&[&[0i64, 0][..], &[1, 0][..]]).unwrap();
        // Normal should be perpendicular to (1,0), i.e., (0, ±1)
        assert_eq!(hp.normal[0], Integer::from(0));
        assert!(hp.normal[1] != 0);
    }
}
