//! 4-form computation helpers.
//!
//! Functions for computing and enumerating 4-form intersection values.

use crate::Point;
use crate::integer_math::determinant_gaussian;
use crate::triangulation::Triangulation;
use malachite::Rational;
use malachite::num::conversion::traits::RoundingFrom;
use malachite::rounding_modes::RoundingMode;
use std::collections::{HashMap, HashSet};

use super::helpers::combinations;

/// Compute known 4-form values for DISTINCT non-origin 4-tuples.
///
/// For distinct (i, j, k, l), κ^V_{ijkl} = 1/|det(pts)|
pub fn compute_known_4form_values(
    tri: &Triangulation,
    points: &[Point],
    origin_idx: Option<usize>,
) -> HashMap<Vec<usize>, f64> {
    let mut known_4form: HashMap<Vec<usize>, f64> = HashMap::new();

    for simplex in tri.simplices() {
        let indices: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| Some(i) != origin_idx)
            .collect();

        // Only compute for distinct 4-tuples (combinations, not with-replacement)
        for subset in combinations(&indices, 4) {
            if known_4form.contains_key(&subset) {
                continue;
            }

            let mut mat_pts: Vec<Vec<Rational>> = subset
                .iter()
                .map(|&idx| {
                    points[idx]
                        .coords()
                        .iter()
                        .map(|&x| Rational::from(x))
                        .collect()
                })
                .collect();

            let det = determinant_gaussian(&mut mat_pts);
            if det != 0 {
                let abs_det = if det < 0 { -det.clone() } else { det };
                let (det_f64, _) = f64::rounding_from(&abs_det, RoundingMode::Nearest);
                known_4form.insert(subset, 1.0 / det_f64);
            }
        }
    }

    known_4form
}

/// Enumerate unknown 4-form variables (self-intersection tuples with at least one repeated index).
pub fn enumerate_4form_variables(
    tri: &Triangulation,
    origin_idx: Option<usize>,
    known_4form: &HashMap<Vec<usize>, f64>,
) -> (Vec<Vec<usize>>, HashMap<Vec<usize>, usize>) {
    let mut vars = HashSet::new();

    for simplex in tri.simplices() {
        let indices: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| Some(i) != origin_idx)
            .collect();

        // Generate all 4-tuples with replacement (i <= j <= k <= l)
        for &i in &indices {
            for &j in &indices {
                if j < i {
                    continue;
                }
                for &k in &indices {
                    if k < j {
                        continue;
                    }
                    for &l in &indices {
                        if l < k {
                            continue;
                        }
                        let tuple = vec![i, j, k, l];

                        // Only include if NOT in known_4form (i.e., not distinct)
                        // This means it has at least one repeated index
                        if !known_4form.contains_key(&tuple) {
                            vars.insert(tuple);
                        }
                    }
                }
            }
        }
    }

    let mut var_list: Vec<Vec<usize>> = vars.into_iter().collect();
    var_list.sort_unstable();

    let var_map: HashMap<Vec<usize>, usize> = var_list
        .iter()
        .enumerate()
        .map(|(i, v)| (v.clone(), i))
        .collect();

    (var_list, var_map)
}

/// Collect 3-face probes (triples of non-origin indices) for GLSM constraints.
///
/// Use ALL 3-tuples (with replacement) to get enough constraints.
pub fn collect_3face_probes(tri: &Triangulation, origin_idx: Option<usize>) -> HashSet<Vec<usize>> {
    let mut probes = HashSet::new();

    for simplex in tri.simplices() {
        let indices: Vec<usize> = simplex
            .iter()
            .copied()
            .filter(|&i| Some(i) != origin_idx)
            .collect();

        // Generate ALL 3-tuples with replacement (i <= j <= k)
        for &i in &indices {
            for &j in &indices {
                if j < i {
                    continue;
                }
                for &k in &indices {
                    if k < j {
                        continue;
                    }
                    probes.insert(vec![i, j, k]);
                }
            }
        }
    }

    probes
}
