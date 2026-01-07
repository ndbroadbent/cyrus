//! Equation construction for the CYTools algorithm.

use super::distinct::DistinctIntnum;
use crate::triangulation::Triangulation;
use std::collections::{HashMap, HashSet};

/// Build c_dict: maps 3-tuples to known 4-form contributions.
///
/// For each distinct 4-form (i,j,k,l,inv_det), add contribution to:
/// - c_dict[(j,k,l)] += (i, inv_det)
/// - c_dict[(i,k,l)] += (j, inv_det)
/// - c_dict[(i,j,l)] += (k, inv_det)
/// - c_dict[(i,j,k)] += (l, inv_det)
pub fn build_c_dict(
    tri: &Triangulation,
    distintnum_array: &[DistinctIntnum],
) -> HashMap<[usize; 3], Vec<(usize, f64)>> {
    // First, collect all 3-tuples from simplices
    let mut simp_3: HashSet<[usize; 3]> = HashSet::new();
    for simplex in tri.simplices() {
        let non_origin: Vec<usize> = simplex.iter().copied().filter(|&i| i != 0).collect();
        for i in 0..non_origin.len() {
            for j in (i + 1)..non_origin.len() {
                for k in (j + 1)..non_origin.len() {
                    let mut triple = [non_origin[i], non_origin[j], non_origin[k]];
                    triple.sort_unstable();
                    simp_3.insert(triple);
                }
            }
        }
    }

    // Initialize c_dict
    let mut c_dict: HashMap<[usize; 3], Vec<(usize, f64)>> = HashMap::new();
    for &triple in &simp_3 {
        c_dict.insert(triple, Vec::new());
    }

    // Add contributions from distinct 4-forms
    for d in distintnum_array {
        let [d0, d1, d2, d3] = d.indices;
        let inv_det = d.inv_det;

        // c_dict[(d1,d2,d3)] += (d0, inv_det)
        let mut key1 = [d1, d2, d3];
        key1.sort_unstable();
        if let Some(v) = c_dict.get_mut(&key1) {
            v.push((d0, inv_det));
        }

        // c_dict[(d0,d2,d3)] += (d1, inv_det)
        let mut key2 = [d0, d2, d3];
        key2.sort_unstable();
        if let Some(v) = c_dict.get_mut(&key2) {
            v.push((d1, inv_det));
        }

        // c_dict[(d0,d1,d3)] += (d2, inv_det)
        let mut key3 = [d0, d1, d3];
        key3.sort_unstable();
        if let Some(v) = c_dict.get_mut(&key3) {
            v.push((d2, inv_det));
        }

        // c_dict[(d0,d1,d2)] += (d3, inv_det)
        let mut key4 = [d0, d1, d2];
        key4.sort_unstable();
        if let Some(v) = c_dict.get_mut(&key4) {
            v.push((d3, inv_det));
        }
    }

    c_dict
}

/// Build equation array and eqn_dict.
///
/// eqn_array: list of 3-tuples (i,j,k)
/// eqn_dict: maps each 3-tuple to list of (point_idx, variable_idx) pairs
pub fn build_eqn_structures(
    tri: &Triangulation,
    n_points: usize,
    variable_array: &[[usize; 4]],
    _variable_dict: &HashMap<[usize; 4], usize>,
) -> (Vec<[usize; 3]>, HashMap<[usize; 3], Vec<(usize, usize)>>) {
    // Collect all 2-tuples and 3-tuples
    let mut simp_2: HashSet<[usize; 2]> = HashSet::new();
    let mut simp_3: HashSet<[usize; 3]> = HashSet::new();

    for simplex in tri.simplices() {
        let non_origin: Vec<usize> = simplex.iter().copied().filter(|&i| i != 0).collect();

        for i in 0..non_origin.len() {
            for j in (i + 1)..non_origin.len() {
                let mut pair = [non_origin[i], non_origin[j]];
                pair.sort_unstable();
                simp_2.insert(pair);

                for k in (j + 1)..non_origin.len() {
                    let mut triple = [non_origin[i], non_origin[j], non_origin[k]];
                    triple.sort_unstable();
                    simp_3.insert(triple);
                }
            }
        }
    }

    // eqn_array_1: all 3-tuples
    let mut eqn_set: HashSet<[usize; 3]> = simp_3.clone();

    // eqn_array_2: from 2-tuples -> (s0,s0,s1), (s0,s1,s1)
    for &[s0, s1] in &simp_2 {
        let mut e1 = [s0, s0, s1];
        e1.sort_unstable();
        eqn_set.insert(e1);

        let mut e2 = [s0, s1, s1];
        e2.sort_unstable();
        eqn_set.insert(e2);
    }

    // eqn_array_3: (i,i,i) for each non-origin point
    for i in 1..n_points {
        eqn_set.insert([i, i, i]);
    }

    let mut eqn_array: Vec<[usize; 3]> = eqn_set.into_iter().collect();
    eqn_array.sort_unstable();

    // Build eqn_dict: for each variable (v0,v1,v2,v3), add contributions to equations
    let mut eqn_dict: HashMap<[usize; 3], Vec<(usize, usize)>> = HashMap::new();
    for &eq in &eqn_array {
        eqn_dict.insert(eq, Vec::new());
    }

    for (var_idx, &v) in variable_array.iter().enumerate() {
        let [v0, v1, v2, v3] = v;

        // Determine which equations this variable contributes to
        // Based on CYTools eqn_dict construction logic
        if v0 == v3 {
            // (v0, v1, v2, v0) -> contributes to (v0, v1, v2)
            let mut key = [v0, v1, v2];
            key.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key) {
                vec.push((v3, var_idx));
            }
        } else if v0 == v2 {
            // (v0, v1, v0, v3) -> contributes to (v0, v1, v2) and (v0, v1, v3)
            let mut key1 = [v0, v1, v2];
            key1.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key1) {
                vec.push((v3, var_idx));
            }
            let mut key2 = [v0, v1, v3];
            key2.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key2) {
                vec.push((v2, var_idx));
            }
        } else if v0 == v1 && v2 == v3 {
            // (v0, v0, v2, v2) -> contributes to (v0, v0, v2) and (v0, v2, v2)
            let mut key1 = [v0, v1, v2];
            key1.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key1) {
                vec.push((v3, var_idx));
            }
            let mut key2 = [v0, v2, v3];
            key2.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key2) {
                vec.push((v1, var_idx));
            }
        } else if v0 == v1 {
            // (v0, v0, v2, v3) -> contributes to (v0, v0, v2), (v0, v0, v3), (v0, v2, v3)
            let mut key1 = [v0, v1, v2];
            key1.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key1) {
                vec.push((v3, var_idx));
            }
            let mut key2 = [v0, v1, v3];
            key2.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key2) {
                vec.push((v2, var_idx));
            }
            let mut key3 = [v0, v2, v3];
            key3.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key3) {
                vec.push((v1, var_idx));
            }
        } else if v1 == v3 {
            // (v0, v1, v2, v1) -> contributes to (v0, v1, v2) and (v1, v2, v3)
            let mut key1 = [v0, v1, v2];
            key1.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key1) {
                vec.push((v3, var_idx));
            }
            let mut key2 = [v1, v2, v3];
            key2.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key2) {
                vec.push((v0, var_idx));
            }
        } else if v1 == v2 {
            // (v0, v1, v1, v3) -> contributes to (v0, v1, v1), (v0, v1, v3), (v1, v1, v3)
            let mut key1 = [v0, v1, v2];
            key1.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key1) {
                vec.push((v3, var_idx));
            }
            let mut key2 = [v0, v1, v3];
            key2.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key2) {
                vec.push((v2, var_idx));
            }
            let mut key3 = [v1, v2, v3];
            key3.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key3) {
                vec.push((v0, var_idx));
            }
        } else if v2 == v3 {
            // (v0, v1, v2, v2) -> contributes to (v0, v1, v2), (v0, v2, v2), (v1, v2, v2)
            let mut key1 = [v0, v1, v2];
            key1.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key1) {
                vec.push((v3, var_idx));
            }
            let mut key2 = [v0, v2, v3];
            key2.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key2) {
                vec.push((v1, var_idx));
            }
            let mut key3 = [v1, v2, v3];
            key3.sort_unstable();
            if let Some(vec) = eqn_dict.get_mut(&key3) {
                vec.push((v0, var_idx));
            }
        }
        // Note: if all indices are distinct, variable shouldn't exist (not in variable_array)
    }

    (eqn_array, eqn_dict)
}
