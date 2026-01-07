//! 3-form computation and origin value derivation.
//!
//! Functions for computing 3-forms from 4-forms and deriving origin-indexed values.

use crate::triangulation::Triangulation;
use crate::types::rational::Rational as TypedRational;
use crate::types::tags::Finite;
use std::collections::{HashMap, HashSet};

use super::Intersection;
use super::helpers::float_to_rational;

/// Compute 3-form values by summing 4-forms: κ_{ijk} = Σ_l κ^V_{ijkl}
///
/// The sum is over ALL non-origin points l, not per-simplex.
/// κ_{ijk} = Σ_{l ∈ rays} κ^V_{ijkl}
pub fn compute_3form_from_4form(
    tri: &Triangulation,
    origin_idx: Option<usize>,
    all_4form: &HashMap<Vec<usize>, f64>,
    n_points: usize,
) -> HashMap<Vec<usize>, f64> {
    let mut kappa_3form: HashMap<Vec<usize>, f64> = HashMap::new();

    // Collect all valid 3-tuples from simplices
    let mut all_3tuples: HashSet<Vec<usize>> = HashSet::new();
    for simplex in tri.simplices() {
        let indices: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| Some(i) != origin_idx)
            .collect();

        // Generate all 3-tuples with replacement
        for &i in &indices {
            for &j in &indices {
                if j < i {
                    continue;
                }
                for &k in &indices {
                    if k < j {
                        continue;
                    }
                    all_3tuples.insert(vec![i, j, k]);
                }
            }
        }
    }

    // Get all non-origin point indices (these are the rays)
    let rays: Vec<usize> = (0..n_points).filter(|&i| Some(i) != origin_idx).collect();

    // For each 3-tuple, sum over ALL rays l
    for tuple3 in all_3tuples {
        let i = tuple3[0];
        let j = tuple3[1];
        let k = tuple3[2];

        let mut sum = 0.0;

        // Sum over ALL rays l (not per-simplex!)
        for &l in &rays {
            let mut key4 = vec![i, j, k, l];
            key4.sort_unstable();

            // If this 4-form exists (either known or solved), add it
            if let Some(&val) = all_4form.get(&key4) {
                sum += val;
            }
            // If not in all_4form, it's 0 (not in any simplex)
        }

        kappa_3form.insert(tuple3, sum);
    }

    kappa_3form
}

/// Derive origin-indexed values using linear equivalence.
///
/// D_0 ~ -Σ_{i>0} D_i implies:
/// - κ_{0jk} = -Σ_{i>0} κ_{ijk}
/// - κ_{00k} = -Σ_{i>0} κ_{i0k}
/// - κ_{000} = -Σ_{i>0} κ_{i00}
pub fn derive_origin_values(
    kappa_3form: &HashMap<Vec<usize>, f64>,
    n_pts: usize,
    origin_idx: Option<usize>,
) -> HashMap<Vec<usize>, f64> {
    let mut full = kappa_3form.clone();

    let origin = match origin_idx {
        Some(o) => o,
        None => return full,
    };

    // Get non-origin indices
    let non_origin: Vec<usize> = (0..n_pts).filter(|&i| i != origin).collect();

    // Step 1: κ_{0jk} = -Σ_{i>0} κ_{ijk} for j,k non-origin
    for &j in &non_origin {
        for &k in &non_origin {
            if k < j {
                continue;
            }
            let sum: f64 = non_origin
                .iter()
                .map(|&i| {
                    let mut key = vec![i, j, k];
                    key.sort_unstable();
                    kappa_3form.get(&key).copied().unwrap_or(0.0)
                })
                .sum();

            let mut key = vec![origin, j, k];
            key.sort_unstable();
            full.insert(key, -sum);
        }
    }

    // Step 2: κ_{00k} = -Σ_{i>0} κ_{i0k} for k non-origin
    for &k in &non_origin {
        let sum: f64 = non_origin
            .iter()
            .map(|&i| {
                let mut key = vec![i, origin, k];
                key.sort_unstable();
                full.get(&key).copied().unwrap_or(0.0)
            })
            .sum();

        let mut key = vec![origin, origin, k];
        key.sort_unstable();
        full.insert(key, -sum);
    }

    // Step 3: κ_{000} = -Σ_{i>0} κ_{i00}
    let sum: f64 = non_origin
        .iter()
        .map(|&i| {
            let mut key = vec![i, origin, origin];
            key.sort_unstable();
            full.get(&key).copied().unwrap_or(0.0)
        })
        .sum();

    full.insert(vec![origin, origin, origin], -sum);

    full
}

/// Convert kappa map to Intersection struct.
pub fn kappa_to_intersection(n_pts: usize, kappa: &HashMap<Vec<usize>, f64>) -> Intersection {
    let mut result = Intersection::new(n_pts);

    for (key, &val) in kappa {
        if val.abs() > 1e-10 && key.len() == 3 {
            let rational = float_to_rational(val);
            result.set(
                key[0],
                key[1],
                key[2],
                TypedRational::<Finite>::from_raw(rational),
            );
        }
    }

    result
}
