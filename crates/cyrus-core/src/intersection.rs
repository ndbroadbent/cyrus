//! Intersection numbers for Calabi-Yau manifolds.
//!
//! The intersection tensor `κ_ijk` encodes the triple intersection numbers
//! of divisors on a Calabi-Yau threefold. It is fully symmetric:
//! `κ_ijk = κ_jik = κ_kij` etc.
//!
//! ## Performance
//!
//! - Sparse storage via `HashMap` for memory efficiency
//! - Parallel iteration via `rayon` for multi-core contraction
//! - Exact arithmetic via `malachite::Rational`
//!
//! Reference: arXiv:2107.09064, and [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]]

use crate::Point;
use crate::error::{Error, Result};
use crate::integer_math::determinant_gaussian;
use crate::triangulation::Triangulation;
use malachite::Integer;
use malachite::Rational;
use malachite::num::arithmetic::traits::Abs;
use malachite::num::conversion::traits::RoundingFrom;
use malachite::rounding_modes::RoundingMode;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// Sparse representation of intersection numbers.
///
/// Only stores unique entries (i ≤ j ≤ k) since the tensor is symmetric.
#[derive(Debug, Clone, Default)]
pub struct Intersection {
    /// Dimension (number of divisors).
    dim: usize,
    /// Map from (i, j, k) with i ≤ j ≤ k to intersection number.
    entries: HashMap<(usize, usize, usize), Rational>,
}

impl Intersection {
    /// Create a new empty intersection tensor of given dimension.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            entries: HashMap::new(),
        }
    }

    /// Get the dimension.
    pub const fn dim(&self) -> usize {
        self.dim
    }

    /// Set an intersection number `κ_ijk`.
    ///
    /// Automatically handles symmetry - only stores canonical form (i ≤ j ≤ k).
    pub fn set(&mut self, i: usize, j: usize, k: usize, value: Rational) {
        let key = canonical_key(i, j, k);
        if value == 0 {
            self.entries.remove(&key);
        } else {
            self.entries.insert(key, value);
        }
    }

    /// Get an intersection number `κ_ijk`.
    pub fn get(&self, i: usize, j: usize, k: usize) -> Rational {
        let key = canonical_key(i, j, k);
        self.entries
            .get(&key)
            .cloned()
            .unwrap_or_else(|| Rational::from(0))
    }

    /// Number of non-zero entries (in canonical form).
    pub fn num_nonzero(&self) -> usize {
        self.entries.len()
    }

    /// Iterate over non-zero entries as ((i, j, k), value) with i ≤ j ≤ k.
    pub fn iter(&self) -> impl Iterator<Item = ((usize, usize, usize), Rational)> + '_ {
        self.entries.iter().map(|(&k, v)| (k, v.clone()))
    }

    /// Parallel iterator over entries.
    pub fn par_iter(&self) -> impl ParallelIterator<Item = ((usize, usize, usize), Rational)> + '_ {
        // Since entries is a HashMap, we can use par_iter directly if we wrap the values
        self.entries.par_iter().map(|(&k, v)| (k, v.clone()))
    }

    /// Compute `κ_ijk t^i t^j t^k` (triple contraction).
    ///
    /// This is used in volume computations: V = (1/6) `κ_ijk` t^i t^j t^k
    ///
    /// # Panics
    /// Panics if `t.len() != self.dim()`.
    pub fn contract_triple(&self, t: &[f64]) -> f64 {
        assert_eq!(t.len(), self.dim, "dimension mismatch");

        // Convert sparse entries to f64 for physics contraction
        self.entries
            .iter()
            .map(|(&(i, j, k), val)| {
                let mult = symmetry_multiplicity(i, j, k);
                let (val_f, _) = f64::rounding_from(val, RoundingMode::Nearest);
                f64::from(mult) * val_f * t[i] * t[j] * t[k]
            })
            .sum()
    }
}

/// Compute canonical key (i ≤ j ≤ k) for symmetric tensor.
fn canonical_key(i: usize, j: usize, k: usize) -> (usize, usize, usize) {
    let mut v = [i, j, k];
    v.sort_unstable();
    v.into()
}

/// Compute symmetry multiplicity for entry (i, j, k) with i ≤ j ≤ k.
const fn symmetry_multiplicity(i: usize, j: usize, k: usize) -> u8 {
    if i == j && j == k {
        1
    } else if i == j || j == k || i == k {
        3
    } else {
        6
    }
}

/// Compute the intersection numbers for a Calabi-Yau hypersurface.
///
/// # Errors
/// Returns an error if the linear system for intersection numbers is unsolvable or singular.
///
/// Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 1.2
pub fn compute_intersection_numbers(
    tri: &Triangulation,
    points: &[Point],
    glsm: &[Vec<Integer>],
) -> Result<Intersection> {
    let n_pts = points.len();
    let dim_v = points[0].dim();

    let var_list = collect_variables(tri, dim_v);
    let var_map: HashMap<Vec<usize>, usize> = var_list
        .iter()
        .enumerate()
        .map(|(i, v)| (v.clone(), i))
        .collect();
    let n_vars = var_list.len();

    let mut m_mat = Vec::new();
    let mut c_vec = Vec::new();

    build_distinct_equations(tri, points, dim_v, &var_map, n_vars, &mut m_mat, &mut c_vec);
    build_glsm_equations(tri, glsm, dim_v, &var_map, n_vars, &mut m_mat, &mut c_vec);

    let sol = solve_rectangular(&m_mat, &c_vec, n_vars)?;
    Ok(reduce_to_cy(n_pts, &var_list, &sol))
}

fn collect_variables(tri: &Triangulation, dim_v: usize) -> Vec<Vec<usize>> {
    let mut variables = HashSet::new();
    for simplex in tri.simplices() {
        if simplex.len() != dim_v + 1 {
            continue;
        }
        let perms = combinations_with_replacement(simplex, dim_v);
        for mut p in perms {
            p.sort_unstable();
            variables.insert(p);
        }
    }
    let mut var_list: Vec<Vec<usize>> = variables.into_iter().collect();
    var_list.sort_unstable();
    var_list
}

fn build_distinct_equations(
    tri: &Triangulation,
    points: &[Point],
    dim_v: usize,
    var_map: &HashMap<Vec<usize>, usize>,
    n_vars: usize,
    m_mat: &mut Vec<Vec<Rational>>,
    c_vec: &mut Vec<Rational>,
) {
    for simplex in tri.simplices() {
        if simplex.len() != dim_v + 1 {
            continue;
        }
        for i in 0..simplex.len() {
            let mut subset: Vec<usize> = simplex
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(_, &v)| v)
                .collect();

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
            let val = if det == 0 {
                Rational::from(0)
            } else {
                Rational::from(1) / det.abs()
            };

            subset.sort_unstable();
            if let Some(&var_idx) = var_map.get(&subset) {
                let mut row = vec![Rational::from(0); n_vars];
                row[var_idx] = Rational::from(1);
                m_mat.push(row);
                c_vec.push(val);
            }
        }
    }
}

fn build_glsm_equations(
    tri: &Triangulation,
    glsm: &[Vec<Integer>],
    dim_v: usize,
    var_map: &HashMap<Vec<usize>, usize>,
    n_vars: usize,
    m_mat: &mut Vec<Vec<Rational>>,
    c_vec: &mut Vec<Rational>,
) {
    let mut probes = HashSet::new();
    for simplex in tri.simplices() {
        if simplex.len() != dim_v + 1 {
            continue;
        }
        let perms = combinations_with_replacement(simplex, dim_v - 1);
        for mut p in perms {
            p.sort_unstable();
            probes.insert(p);
        }
    }

    for probe in probes {
        for q_row in glsm {
            let mut row = vec![Rational::from(0); n_vars];
            for (m, q_val) in q_row.iter().enumerate() {
                if *q_val == 0 {
                    continue;
                }
                let mut key = probe.clone();
                key.push(m);
                key.sort_unstable();
                if let Some(&var_idx) = var_map.get(&key) {
                    row[var_idx] += Rational::from(q_val.clone());
                }
            }
            if row.iter().any(|x| *x != 0) {
                m_mat.push(row);
                c_vec.push(Rational::from(0));
            }
        }
    }
}

fn reduce_to_cy(n_pts: usize, var_list: &[Vec<usize>], sol: &[Rational]) -> Intersection {
    let mut kappa = Intersection::new(n_pts);
    for (key, val) in var_list.iter().zip(sol.iter()) {
        let v0 = key[0];
        let v1 = key[1];
        let v2 = key[2];
        let v3 = key[3];
        update_kappa(&mut kappa, v0, v1, v2, val);
        update_kappa(&mut kappa, v0, v1, v3, val);
        update_kappa(&mut kappa, v0, v2, v3, val);
        update_kappa(&mut kappa, v1, v2, v3, val);
    }
    kappa
}

fn update_kappa(kappa: &mut Intersection, i: usize, j: usize, k: usize, val: &Rational) {
    let current = kappa.get(i, j, k);
    kappa.set(i, j, k, current + val);
}

fn combinations_with_replacement(pool: &[usize], k: usize) -> Vec<Vec<usize>> {
    if k == 0 {
        return vec![vec![]];
    }
    if pool.is_empty() {
        return vec![];
    }

    let mut res = Vec::new();
    for (i, &item) in pool.iter().enumerate() {
        let sub_combs = combinations_with_replacement(&pool[i..], k - 1);
        for mut c in sub_combs {
            let mut new_c = vec![item];
            new_c.append(&mut c);
            res.push(new_c);
        }
    }
    res
}

/// Solve rectangular system Mx=C using Gaussian elimination.
///
/// # Errors
/// Returns an error if the system is underdetermined or rank deficient.
fn solve_rectangular(m: &[Vec<Rational>], c: &[Rational], n_vars: usize) -> Result<Vec<Rational>> {
    let n_eq = m.len();
    if n_eq < n_vars {
        return Err(Error::SingularMatrix("Underdetermined system".into()));
    }

    let mut mat = vec![vec![Rational::from(0); n_vars + 1]; n_eq];
    for (i, row) in mat.iter_mut().enumerate().take(n_eq) {
        for (j, item) in row.iter_mut().enumerate().take(n_vars) {
            *item = m[i][j].clone();
        }
        row[n_vars] = c[i].clone();
    }

    let mut pivot_row = 0;
    for col in 0..n_vars {
        if pivot_row >= n_eq {
            break;
        }
        if let Some(r) = find_pivot(&mat, pivot_row, col, n_eq) {
            mat.swap(pivot_row, r);
            normalize_row(&mut mat[pivot_row], col, n_vars);
            eliminate_other_rows(&mut mat, pivot_row, col, n_vars);
            pivot_row += 1;
        }
    }

    if pivot_row < n_vars {
        return Err(Error::SingularMatrix("Rank deficient system".into()));
    }

    let mut x = Vec::new();
    for row in mat.iter().take(n_vars) {
        x.push(row[n_vars].clone());
    }
    Ok(x)
}

fn find_pivot(mat: &[Vec<Rational>], pivot_row: usize, col: usize, n_eq: usize) -> Option<usize> {
    mat.iter()
        .enumerate()
        .take(n_eq)
        .skip(pivot_row)
        .find(|(_, row)| row[col] != 0)
        .map(|(r, _)| r)
}

fn normalize_row(row: &mut [Rational], col: usize, n_vars: usize) {
    let div = row[col].clone();
    for item in row.iter_mut().take(n_vars + 1).skip(col) {
        *item /= &div;
    }
}

fn eliminate_other_rows(mat: &mut [Vec<Rational>], pivot_row: usize, col: usize, n_vars: usize) {
    let pivot_row_vals = mat[pivot_row].clone();
    for (r2, row) in mat.iter_mut().enumerate() {
        if r2 != pivot_row {
            let factor = row[col].clone();
            if factor != 0 {
                for (j, item) in row.iter_mut().enumerate().take(n_vars + 1).skip(col) {
                    let sub = &factor * &pivot_row_vals[j];
                    *item -= sub;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use malachite::Rational;

    #[test]
    fn test_set_get() {
        let mut kappa = Intersection::new(3);
        kappa.set(0, 1, 2, Rational::from(5));

        assert_eq!(kappa.get(0, 1, 2), Rational::from(5));
        assert_eq!(kappa.get(2, 1, 0), Rational::from(5));
    }
}
