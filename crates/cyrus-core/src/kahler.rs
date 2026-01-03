//! Kähler cone and Mori cone computations.
//!
//! The Kähler cone is the set of Kähler moduli that yield positive volumes
//! for all curves in the Mori cone.
//!
//! Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 3.

use crate::Point;
use crate::error::{Error, Result};
use crate::f64_pos;
use crate::integer_math::integer_kernel;
use crate::triangulation::Triangulation;
use crate::types::f64::F64;
use crate::types::i64::I64;
use crate::types::tags::Finite;
use malachite::Integer;
use std::collections::HashMap;

/// Mori cone generators for a given triangulation.
#[derive(Debug, Clone)]
pub struct MoriCone {
    /// Generators as integer vectors in the space of divisors.
    generators: Vec<Vec<I64<Finite>>>,
}

impl MoriCone {
    /// Create a new Mori cone from generators.
    pub const fn new(generators: Vec<Vec<I64<Finite>>>) -> Self {
        Self { generators }
    }

    /// Create from raw Integer vectors (for compatibility with integer_kernel).
    ///
    /// # Panics
    /// Panics if any generator value does not fit in an i64.
    pub fn from_integers(generators: Vec<Vec<Integer>>) -> Self {
        let typed_generators = generators
            .into_iter()
            .map(|generator| {
                generator
                    .into_iter()
                    .map(|v| {
                        let val = i64::try_from(&v).expect("generator value must fit in i64");
                        I64::<Finite>::new(val)
                    })
                    .collect()
            })
            .collect();
        Self {
            generators: typed_generators,
        }
    }

    /// Get the generators.
    pub fn generators(&self) -> &[Vec<I64<Finite>>] {
        &self.generators
    }

    /// Check if a Kähler class `p` is in the Kähler cone.
    ///
    /// A class `p` is in the Kähler cone if `p · R > 0` for all Mori generators `R`.
    ///
    /// # Panics
    /// Panics if `p.len()` does not match the generator dimension.
    pub fn contains(&self, p: &[F64<Finite>]) -> bool {
        // Threshold for numerical stability (positive)
        let threshold = f64_pos!(1e-10);

        for generator in &self.generators {
            // Dot product: Finite * Finite = Finite, sum of Finite = Finite
            let dot: F64<Finite> = generator
                .iter()
                .zip(p.iter())
                .map(|(g, pi)| g.to_f64() * *pi)
                .fold(F64::<Finite>::ZERO, |acc, x| acc + x);

            // Must be strictly positive (above threshold)
            if dot.get() <= threshold.get() {
                return false;
            }
        }
        true
    }
}

/// Compute the Mori cone generators for a triangulation.
///
/// For each interior ridge (shared by two simplices), there is a unique
/// relation that generates a ray of the Mori cone.
///
/// # Errors
/// Returns an error if no points are provided.
///
/// # Panics
/// Panics if a simplex has no vertex outside its ridge (invariant violation).
pub fn compute_mori_generators(tri: &Triangulation, points: &[Point]) -> Result<MoriCone> {
    let n_pts = points.len();
    if n_pts == 0 {
        return Err(Error::InvalidInput("No points provided".into()));
    }
    let dim = points[0].dim();

    // 1. Identify ridges and their adjacent simplices
    let ridge_map = build_ridge_map(tri, dim);

    let mut generators = Vec::new();

    // 2. For each interior ridge (shared by exactly 2 simplices)
    for (ridge, adj_simplices) in ridge_map {
        if adj_simplices.len() == 2 {
            let s1 = &tri.simplices()[adj_simplices[0]];
            let s2 = &tri.simplices()[adj_simplices[1]];

            // Opposing vertices u1, u2
            let u1 = *s1
                .iter()
                .find(|&v| !ridge.contains(v))
                .expect("Invariant violated: simplex has no vertex outside ridge");
            let u2 = *s2
                .iter()
                .find(|&v| !ridge.contains(v))
                .expect("Invariant violated: simplex has no vertex outside ridge");

            // Find relation among {u1, u2} U ridge
            let mut subset = vec![u1, u2];
            subset.extend(ridge.iter());

            let mut mat = vec![vec![Integer::from(0); subset.len()]; dim + 1];
            for (j, &idx) in subset.iter().enumerate() {
                let coords = points[idx].coords();
                for (i, &c) in coords.iter().enumerate() {
                    mat[i][j] = Integer::from(c);
                }
                mat[dim][j] = Integer::from(1);
            }

            let kernel = integer_kernel(&mat);
            if let Some(relation) = extract_mori_relation(&kernel) {
                // relation is in order of 'subset': [u1, u2, v1, v2, v3, v4]
                // Expand to full points list
                let mut full_gen = vec![I64::<Finite>::ZERO; n_pts];
                for (i, &idx) in subset.iter().enumerate() {
                    let val = i64::try_from(&relation[i]).expect("relation value must fit in i64");
                    full_gen[idx] = I64::<Finite>::new(val);
                }
                generators.push(full_gen);
            }
        }
    }

    Ok(MoriCone::new(generators))
}

fn build_ridge_map(tri: &Triangulation, dim: usize) -> HashMap<Vec<usize>, Vec<usize>> {
    let mut ridge_map: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();
    for (s_idx, simplex) in tri.simplices().iter().enumerate() {
        if simplex.len() != dim + 1 {
            continue;
        }

        for i in 0..simplex.len() {
            let mut ridge = simplex.clone();
            ridge.remove(i);
            ridge.sort_unstable();
            ridge_map.entry(ridge).or_default().push(s_idx);
        }
    }
    ridge_map
}

/// From the kernel of {u1, u2, ridge}, extract the Mori relation.
fn extract_mori_relation(kernel: &[Vec<Integer>]) -> Option<Vec<Integer>> {
    if !kernel.is_empty() {
        let mut res = kernel[0].clone();
        if res[0] < 0 {
            for val in &mut res {
                *val = -val.clone();
            }
        }
        return Some(res);
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn finite_i64(v: i64) -> I64<Finite> {
        I64::<Finite>::new(v)
    }

    fn finite_f64(v: f64) -> F64<Finite> {
        F64::<Finite>::new(v).unwrap()
    }

    #[test]
    fn test_kahler_cone_contains() {
        // Generator R = [1, -2]
        // p = [2, 1] => p·R = 2 - 2 = 0 (boundary)
        // p = [3, 1] => p·R = 3 - 2 = 1 > 0 (inside)
        let mori = MoriCone::new(vec![vec![finite_i64(1), finite_i64(-2)]]);
        assert!(!mori.contains(&[finite_f64(2.0), finite_f64(1.0)]));
        assert!(mori.contains(&[finite_f64(3.0), finite_f64(1.0)]));
    }

    #[test]
    fn test_mori_cone_from_integers() {
        // Test the from_integers constructor
        let generators = vec![
            vec![Integer::from(1), Integer::from(-1)],
            vec![Integer::from(2), Integer::from(3)],
        ];

        let mori = MoriCone::from_integers(generators);
        assert_eq!(mori.generators().len(), 2);
        assert_eq!(mori.generators()[0][0].get(), 1);
        assert_eq!(mori.generators()[0][1].get(), -1);
        assert_eq!(mori.generators()[1][0].get(), 2);
        assert_eq!(mori.generators()[1][1].get(), 3);
    }

    #[test]
    fn test_extract_mori_relation_empty() {
        // Empty kernel should return None
        let kernel: Vec<Vec<Integer>> = vec![];
        assert!(extract_mori_relation(&kernel).is_none());
    }

    #[test]
    fn test_extract_mori_relation_negative_first() {
        // If first element is negative, all should be negated
        let kernel = vec![vec![Integer::from(-2), Integer::from(3), Integer::from(-1)]];
        let result = extract_mori_relation(&kernel).unwrap();
        assert_eq!(result[0], Integer::from(2));
        assert_eq!(result[1], Integer::from(-3));
        assert_eq!(result[2], Integer::from(1));
    }

    #[test]
    fn test_extract_mori_relation_positive_first() {
        // If first element is positive, should stay the same
        let kernel = vec![vec![Integer::from(2), Integer::from(-3), Integer::from(1)]];
        let result = extract_mori_relation(&kernel).unwrap();
        assert_eq!(result[0], Integer::from(2));
        assert_eq!(result[1], Integer::from(-3));
        assert_eq!(result[2], Integer::from(1));
    }

    #[test]
    fn test_compute_mori_generators_simple() {
        // Two triangles sharing an edge in 2D
        let points = vec![
            Point::new(vec![0, 0]),
            Point::new(vec![1, 0]),
            Point::new(vec![0, 1]),
            Point::new(vec![1, 1]),
        ];

        let tri = Triangulation::new(vec![vec![0, 1, 2], vec![1, 2, 3]]);

        let mori = compute_mori_generators(&tri, &points).unwrap();
        assert!(!mori.generators.is_empty());

        let generator = &mori.generators[0];
        assert_eq!(generator[0].get(), 1);
        assert_eq!(generator[3].get(), 1);
        assert_eq!(generator[1].get(), -1);
        assert_eq!(generator[2].get(), -1);
    }
}
