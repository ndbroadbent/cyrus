//! Extraction of intersection numbers from solved 4-form values.

use super::distinct::DistinctIntnum;
use crate::intersection::Intersection;
use crate::types::rational::Rational as TypedRational;
use crate::types::tags::Finite;
use malachite::Rational;
use std::collections::{HashMap, HashSet};

/// Extract intersection numbers from solved 4-form values.
///
/// kappa_{ijk} = sum_l kappa^V_{sorted(i,j,k,l)}
pub fn extract_intersection_numbers(
    variable_array: &[[usize; 4]],
    solution: &[f64],
    distintnum_array: &[DistinctIntnum],
    n_points: usize,
) -> Intersection {
    // Build map of all 4-form values (known + solved)
    let mut kappa_4form: HashMap<[usize; 4], f64> = HashMap::new();

    // Known values
    for d in distintnum_array {
        kappa_4form.insert(d.indices, d.inv_det);
    }

    // Solved values
    for (idx, &var) in variable_array.iter().enumerate() {
        kappa_4form.insert(var, solution[idx]);
    }

    // Collect all non-origin 3-tuples that appear in the triangulation
    let mut three_tuples: HashSet<[usize; 3]> = HashSet::new();

    // From 4-forms, extract all 3-tuples
    for &[a, b, c, d] in kappa_4form.keys() {
        // Each 4-form (a,b,c,d) gives us 4 possible 3-tuples
        let tuples = [[b, c, d], [a, c, d], [a, b, d], [a, b, c]];
        for mut t in tuples {
            t.sort_unstable();
            if t[0] > 0 {
                // Non-origin only for now
                three_tuples.insert(t);
            }
        }
    }

    // Compute kappa_{ijk} = sum_l kappa^V_{sorted(i,j,k,l)} for each 3-tuple
    let mut kappa_3form: HashMap<[usize; 3], f64> = HashMap::new();

    for &[i, j, k] in &three_tuples {
        let mut sum = 0.0;
        // Sum over all possible 4th indices l
        for l in 1..n_points {
            let mut key4 = [i, j, k, l];
            key4.sort_unstable();
            if let Some(&val) = kappa_4form.get(&key4) {
                sum += val;
            }
        }
        if sum.abs() > 1e-12 {
            kappa_3form.insert([i, j, k], sum);
        }
    }

    // Now derive origin-indexed values using D_0 ~ -sum_{m>0} D_m
    derive_origin_values(&mut kappa_3form, n_points);

    // Convert to Intersection struct
    let mut result = Intersection::new(n_points);
    for (&[i, j, k], &val) in &kappa_3form {
        if val.abs() > 1e-10 {
            let rational = float_to_rational(val);
            result.set(i, j, k, TypedRational::<Finite>::from_raw(rational));
        }
    }

    result
}

/// Derive origin-indexed values using linear equivalence D_0 ~ -sum_{m>0} D_m.
fn derive_origin_values(kappa: &mut HashMap<[usize; 3], f64>, n_points: usize) {
    let non_origin: Vec<usize> = (1..n_points).collect();

    // kappa_{0jk} = -sum_{m>0} kappa_{mjk}
    for &j in &non_origin {
        for &k in &non_origin {
            if k < j {
                continue;
            }
            let sum: f64 = non_origin
                .iter()
                .map(|&m| {
                    let mut key = [m, j, k];
                    key.sort_unstable();
                    kappa.get(&key).copied().unwrap_or(0.0)
                })
                .sum();

            let mut key = [0, j, k];
            key.sort_unstable();
            kappa.insert(key, -sum);
        }
    }

    // kappa_{00k} = -sum_{m>0} kappa_{m0k} = -sum_{m>0} kappa_{0mk}
    for &k in &non_origin {
        let sum: f64 = non_origin
            .iter()
            .map(|&m| {
                let mut key = [0, m, k];
                key.sort_unstable();
                kappa.get(&key).copied().unwrap_or(0.0)
            })
            .sum();

        kappa.insert([0, 0, k], -sum);
    }

    // kappa_{000} = -sum_{m>0} kappa_{m00} = -sum_{m>0} kappa_{00m}
    let sum: f64 = non_origin
        .iter()
        .map(|&m| kappa.get(&[0, 0, m]).copied().unwrap_or(0.0))
        .sum();
    kappa.insert([0, 0, 0], -sum);
}

/// Convert float to rational.
pub fn float_to_rational(x: f64) -> Rational {
    let rounded = x.round();
    if (x - rounded).abs() < 1e-6 {
        #[allow(clippy::cast_possible_truncation)]
        return Rational::from(rounded as i64);
    }

    for denom in [2, 3, 4, 5, 6, 8, 10, 12] {
        let scaled = x * f64::from(denom);
        let rounded_scaled = scaled.round();
        if (scaled - rounded_scaled).abs() < 1e-6 {
            #[allow(clippy::cast_possible_truncation)]
            return Rational::from(rounded_scaled as i64) / Rational::from(denom);
        }
    }

    // Continued fraction approximation
    continued_fraction_approx(x, 1000)
}

fn continued_fraction_approx(x: f64, max_denom: i64) -> Rational {
    let sign = x.signum();
    let x = x.abs();

    let mut a = x.floor();
    let mut h1 = 1_i64;
    #[allow(clippy::cast_possible_truncation)]
    let mut h0 = a as i64;
    let mut k1 = 0_i64;
    let mut k0 = 1_i64;

    let mut remainder = x - a;

    while remainder > 1e-12 && k0.abs() < max_denom {
        let inv = 1.0 / remainder;
        a = inv.floor();

        #[allow(clippy::cast_possible_truncation)]
        let a_i64 = a as i64;
        let h2 = a_i64 * h0 + h1;
        let k2 = a_i64 * k0 + k1;

        if k2.abs() > max_denom {
            break;
        }

        h1 = h0;
        h0 = h2;
        k1 = k0;
        k0 = k2;

        remainder = inv - a;
    }

    if sign < 0.0 {
        Rational::from(-h0) / Rational::from(k0)
    } else {
        Rational::from(h0) / Rational::from(k0)
    }
}
