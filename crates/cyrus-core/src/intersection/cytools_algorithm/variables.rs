//! Variable array construction for the CYTools algorithm.

use crate::triangulation::Triangulation;
use std::collections::{HashMap, HashSet};

/// Build the variable array: 4-tuples with at least one repeated index.
///
/// From CYTools:
/// - variable_array_1: from 3-tuples (i,j,k) -> (i,i,j,k), (i,j,j,k), (i,j,k,k)
/// - variable_array_2: from 2-tuples (i,j) -> (i,i,j,j), (i,i,i,j), (i,j,j,j)
/// - variable_array_3: self-intersections (i,i,i,i)
pub fn build_variable_array(
    tri: &Triangulation,
    n_points: usize,
) -> (Vec<[usize; 4]>, HashMap<[usize; 4], usize>) {
    let mut variables: HashSet<[usize; 4]> = HashSet::new();

    // Collect all 2-tuples and 3-tuples from simplices (excluding origin)
    let mut simp_2: HashSet<[usize; 2]> = HashSet::new();
    let mut simp_3: HashSet<[usize; 3]> = HashSet::new();

    for simplex in tri.simplices() {
        let non_origin: Vec<usize> = simplex.iter().copied().filter(|&i| i != 0).collect();

        // 2-combinations
        for i in 0..non_origin.len() {
            for j in (i + 1)..non_origin.len() {
                let mut pair = [non_origin[i], non_origin[j]];
                pair.sort_unstable();
                simp_2.insert(pair);
            }
        }

        // 3-combinations
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

    // variable_array_1: from 3-tuples
    for &[s0, s1, s2] in &simp_3 {
        // (s0, s0, s1, s2)
        let mut v1 = [s0, s0, s1, s2];
        v1.sort_unstable();
        variables.insert(v1);

        // (s0, s1, s1, s2)
        let mut v2 = [s0, s1, s1, s2];
        v2.sort_unstable();
        variables.insert(v2);

        // (s0, s1, s2, s2)
        let mut v3 = [s0, s1, s2, s2];
        v3.sort_unstable();
        variables.insert(v3);
    }

    // variable_array_2: from 2-tuples
    for &[s0, s1] in &simp_2 {
        // (s0, s0, s1, s1)
        let mut v1 = [s0, s0, s1, s1];
        v1.sort_unstable();
        variables.insert(v1);

        // (s0, s0, s0, s1)
        let mut v2 = [s0, s0, s0, s1];
        v2.sort_unstable();
        variables.insert(v2);

        // (s0, s1, s1, s1)
        let mut v3 = [s0, s1, s1, s1];
        v3.sort_unstable();
        variables.insert(v3);
    }

    // variable_array_3: self-intersections (i,i,i,i) for each non-origin point
    for i in 1..n_points {
        variables.insert([i, i, i, i]);
    }

    // Sort and create index map
    let mut var_array: Vec<[usize; 4]> = variables.into_iter().collect();
    var_array.sort_unstable();

    let var_dict: HashMap<[usize; 4], usize> = var_array
        .iter()
        .enumerate()
        .map(|(idx, &v)| (v, idx))
        .collect();

    (var_array, var_dict)
}
