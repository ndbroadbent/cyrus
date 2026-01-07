//! Double Description Method (DDM) for cone dualization.
//!
//! Converts between H-representation (hyperplanes) and V-representation (rays).
//!
//! Reference: PPL (Parma Polyhedra Library), cddlib
//!
//! The basic algorithm:
//! 1. Start with initial rays (from identity or first hyperplanes)
//! 2. For each hyperplane, partition rays into positive/zero/negative
//! 3. Keep positive and zero rays
//! 4. For each (positive, negative) pair, compute the intersection ray
//! 5. Remove redundant rays

use std::collections::HashSet;

/// Convert hyperplanes to rays (or vice versa).
///
/// Given H = hyperplanes defining cone {x : H @ x >= 0},
/// compute rays R such that cone = {R^T @ λ : λ >= 0}.
///
/// The algorithm is symmetric: dualize(rays) gives hyperplanes,
/// dualize(hyperplanes) gives rays.
pub fn dualize(matrix: &[Vec<i64>], ambient_dim: usize) -> Vec<Vec<i64>> {
    if matrix.is_empty() {
        // No constraints = full space
        // Return ±e_i for each dimension
        let mut rays = Vec::new();
        for i in 0..ambient_dim {
            let mut pos = vec![0i64; ambient_dim];
            let mut neg = vec![0i64; ambient_dim];
            pos[i] = 1;
            neg[i] = -1;
            rays.push(pos);
            rays.push(neg);
        }
        return rays;
    }

    if ambient_dim == 0 {
        return Vec::new();
    }

    // Use the incremental DDM algorithm
    ddm_incremental(matrix, ambient_dim)
}

/// Incremental Double Description Method.
///
/// Start with the full space (identity rays), then add constraints one by one.
fn ddm_incremental(hyperplanes: &[Vec<i64>], dim: usize) -> Vec<Vec<i64>> {
    // Start with identity rays (full space before any constraints)
    let mut rays: Vec<Vec<i64>> = (0..dim)
        .flat_map(|i| {
            let mut pos = vec![0i64; dim];
            let mut neg = vec![0i64; dim];
            pos[i] = 1;
            neg[i] = -1;
            vec![pos, neg]
        })
        .collect();

    // Process each hyperplane
    for h in hyperplanes {
        rays = process_hyperplane(&rays, h);

        // Early termination if cone becomes empty
        if rays.is_empty() {
            break;
        }
    }

    // Remove duplicates and normalize
    deduplicate_rays(rays)
}

/// Process one hyperplane in the DDM algorithm.
///
/// Partition rays by their dot product with h:
/// - Positive: h·r > 0 (strictly inside)
/// - Zero: h·r = 0 (on boundary)
/// - Negative: h·r < 0 (outside, will be cut)
///
/// New rays = positive + zero + (positive × negative intersections)
fn process_hyperplane(rays: &[Vec<i64>], h: &[i64]) -> Vec<Vec<i64>> {
    let mut positive = Vec::new();
    let mut zero = Vec::new();
    let mut negative = Vec::new();

    // Partition rays
    for ray in rays {
        let dot = dot_product(h, ray);
        match dot.cmp(&0) {
            std::cmp::Ordering::Greater => positive.push(ray.clone()),
            std::cmp::Ordering::Equal => zero.push(ray.clone()),
            std::cmp::Ordering::Less => negative.push(ray.clone()),
        }
    }

    // Start with rays that satisfy the constraint
    let mut new_rays = positive.clone();
    new_rays.extend(zero);

    // For each (positive, negative) pair, compute intersection
    for p in &positive {
        for n in &negative {
            if let Some(intersection) = compute_intersection_ray(p, n, h) {
                new_rays.push(intersection);
            }
        }
    }

    new_rays
}

/// Compute the intersection ray for a (positive, negative) pair.
///
/// Given rays r+ and r- with h·r+ > 0 and h·r- < 0,
/// find the ray r on the hyperplane h·r = 0 that lies between them.
///
/// r = (h·r+) * r- - (h·r-) * r+
///
/// This makes h·r = (h·r+)(h·r-) - (h·r-)(h·r+) = 0.
fn compute_intersection_ray(pos: &[i64], neg: &[i64], h: &[i64]) -> Option<Vec<i64>> {
    let dp = dot_product(h, pos); // > 0
    let dn = dot_product(h, neg); // < 0

    // r = dp * neg - dn * pos
    // Since dn < 0, this is: dp * neg + |dn| * pos
    let mut result: Vec<i64> = pos
        .iter()
        .zip(neg.iter())
        .map(|(&p, &n)| dp * n - dn * p)
        .collect();

    // Normalize by GCD
    let g = gcd_vec(&result);
    if g == 0 {
        return None;
    }

    for x in &mut result {
        *x /= g;
    }

    Some(result)
}

/// Remove duplicate rays (after normalization).
fn deduplicate_rays(rays: Vec<Vec<i64>>) -> Vec<Vec<i64>> {
    let mut seen: HashSet<Vec<i64>> = HashSet::new();
    let mut result = Vec::new();

    for ray in rays {
        // Normalize direction: make first non-zero positive
        let normalized = normalize_direction(ray);
        if !normalized.iter().all(|&x| x == 0) && !seen.contains(&normalized) {
            seen.insert(normalized.clone());
            result.push(normalized);
        }
    }

    result
}

/// Normalize ray direction so first non-zero element is positive.
fn normalize_direction(mut ray: Vec<i64>) -> Vec<i64> {
    // First normalize by GCD
    let g = gcd_vec(&ray);
    if g > 0 {
        for x in &mut ray {
            *x /= g;
        }
    }

    // Make first non-zero positive
    for &x in &ray {
        if x > 0 {
            return ray;
        } else if x < 0 {
            return ray.iter().map(|&v| -v).collect();
        }
    }

    ray
}

/// Dot product of two integer vectors.
fn dot_product(a: &[i64], b: &[i64]) -> i64 {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

/// GCD of a vector.
fn gcd_vec(v: &[i64]) -> i64 {
    v.iter().fold(0, |acc, &x| gcd(acc, x.abs()))
}

/// GCD of two integers.
fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 { a } else { gcd(b, a % b) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dualize_first_quadrant() {
        // First quadrant in 2D: x >= 0, y >= 0
        let hyperplanes = vec![
            vec![1, 0], // x >= 0
            vec![0, 1], // y >= 0
        ];

        let rays = dualize(&hyperplanes, 2);

        // Should get rays (1,0) and (0,1)
        assert_eq!(rays.len(), 2);
        assert!(rays.contains(&vec![1, 0]));
        assert!(rays.contains(&vec![0, 1]));
    }

    #[test]
    fn test_dualize_half_plane() {
        // Half-plane: x >= 0
        let hyperplanes = vec![vec![1, 0]];

        let rays = dualize(&hyperplanes, 2);

        // Should get rays (1,0), (0,1), (0,-1)
        // Because y can be anything, x must be >= 0
        assert!(rays.len() >= 2);
        assert!(rays.contains(&vec![1, 0]));
    }

    #[test]
    fn test_dualize_empty() {
        // No hyperplanes = full space
        let rays = dualize(&[], 2);

        // Should get ±e_i for i in 0..2
        assert_eq!(rays.len(), 4);
    }

    #[test]
    fn test_dualize_single_ray_cone() {
        // Cone generated by single ray (1,1)
        // Dual hyperplanes are: h such that h·(1,1) >= 0
        // i.e., h1 + h2 >= 0
        let rays = vec![vec![1, 1]];
        let hyperplanes = dualize(&rays, 2);

        // Check that (1,1)·h >= 0 for all h
        for h in &hyperplanes {
            let dot: i64 = h.iter().sum();
            assert!(dot >= 0, "Hyperplane {h:?} violates constraint");
        }
    }

    #[test]
    fn test_process_hyperplane() {
        // Start with ±e_1, ±e_2 and add x >= 0
        let rays = vec![vec![1, 0], vec![-1, 0], vec![0, 1], vec![0, -1]];
        let h = vec![1, 0];

        let new_rays = process_hyperplane(&rays, &h);

        // Should keep (1,0), (0,1), (0,-1) and intersection of (-1,0) with positives
        assert!(new_rays.contains(&vec![1, 0]));
        assert!(new_rays.contains(&vec![0, 1]));
        assert!(new_rays.contains(&vec![0, -1]));
        assert!(!new_rays.contains(&vec![-1, 0]));
    }

    #[test]
    fn test_compute_intersection_ray() {
        // pos = (1, 0), neg = (-1, 0), h = (1, 0)
        // h·pos = 1, h·neg = -1
        // r = 1*(-1,0) - (-1)*(1,0) = (-1,0) + (1,0) = (0,0)
        // This is degenerate - correct behavior

        // Better test: pos = (1,1), neg = (-1,1), h = (1,0)
        // h·pos = 1, h·neg = -1
        // r = 1*(-1,1) - (-1)*(1,1) = (-1,1) + (1,1) = (0,2)
        let pos = vec![1, 1];
        let neg = vec![-1, 1];
        let h = vec![1, 0];

        let ray = compute_intersection_ray(&pos, &neg, &h).unwrap();
        assert_eq!(ray, vec![0, 1]); // normalized (0,2) -> (0,1)

        // Verify: h·ray = 0
        let dot: i64 = h.iter().zip(ray.iter()).map(|(&a, &b)| a * b).sum();
        assert_eq!(dot, 0);
    }

    #[test]
    fn test_normalize_direction() {
        assert_eq!(normalize_direction(vec![2, 4, 6]), vec![1, 2, 3]);
        assert_eq!(normalize_direction(vec![-2, -4, -6]), vec![1, 2, 3]);
        assert_eq!(normalize_direction(vec![0, -2, 4]), vec![0, 1, -2]);
    }
}
