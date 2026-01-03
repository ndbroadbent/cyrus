//! Intersection number computation algorithm.

use super::{Intersection, canonical_key};
use crate::Point;
use crate::error::{Error, Result};
use crate::integer_math::determinant_gaussian;
use crate::triangulation::Triangulation;
use crate::types::rational::Rational as TypedRational;
use crate::types::tags::Finite;
use malachite::num::arithmetic::traits::Abs;
use malachite::{Integer, Rational};
use std::collections::{HashMap, HashSet};

struct SystemBuilder<'a> {
    tri: &'a Triangulation,
    points: &'a [Point],
    dim_v: usize,
    origin_idx: Option<usize>,
    var_map: &'a HashMap<Vec<usize>, usize>,
    n_vars: usize,
    m_mat: Vec<Vec<Rational>>,
    c_vec: Vec<Rational>,
}

impl<'a> SystemBuilder<'a> {
    const fn new(
        tri: &'a Triangulation,
        points: &'a [Point],
        dim_v: usize,
        origin_idx: Option<usize>,
        var_map: &'a HashMap<Vec<usize>, usize>,
        n_vars: usize,
    ) -> Self {
        Self {
            tri,
            points,
            dim_v,
            origin_idx,
            var_map,
            n_vars,
            m_mat: Vec::new(),
            c_vec: Vec::new(),
        }
    }

    fn build_distinct_equations(&mut self) {
        for simplex in self.tri.simplices() {
            let rays: Vec<usize> = simplex
                .iter()
                .filter(|&&i| Some(i) != self.origin_idx)
                .copied()
                .collect();
            if rays.len() == self.dim_v {
                let mut mat_pts: Vec<Vec<Rational>> = rays
                    .iter()
                    .map(|&idx| {
                        self.points[idx]
                            .coords()
                            .iter()
                            .map(|&x| Rational::from(x))
                            .collect()
                    })
                    .collect();

                let det = determinant_gaussian(&mut mat_pts);
                if det != 0 {
                    let val = Rational::from(1) / det.abs();
                    let mut key = rays;
                    key.sort_unstable();
                    if let Some(&var_idx) = self.var_map.get(&key) {
                        let mut row = vec![Rational::from(0); self.n_vars];
                        row[var_idx] = Rational::from(1);
                        self.m_mat.push(row);
                        self.c_vec.push(val);
                    }
                }
            }
        }
    }

    fn build_glsm_equations(&mut self, glsm: &[Vec<Integer>]) {
        let mut probes = HashSet::new();
        for key in self.var_map.keys() {
            for i in 0..key.len() {
                let mut probe = key.clone();
                probe.remove(i);
                probes.insert(probe);
            }
        }

        for probe in probes {
            for q_row in glsm {
                let mut row = vec![Rational::from(0); self.n_vars];
                for (m, q_val) in q_row.iter().enumerate() {
                    if *q_val == 0 || Some(m) == self.origin_idx {
                        continue;
                    }
                    let mut key = probe.clone();
                    key.push(m);
                    key.sort_unstable();
                    if let Some(&var_idx) = self.var_map.get(&key) {
                        row[var_idx] += Rational::from(q_val.clone());
                    }
                }
                if row.iter().any(|x| *x != 0) {
                    self.m_mat.push(row);
                    self.c_vec.push(Rational::from(0));
                }
            }
        }
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

    let origin_idx = points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0));

    let var_list = collect_variables(tri, dim_v, origin_idx);
    let var_map: HashMap<Vec<usize>, usize> = var_list
        .iter()
        .enumerate()
        .map(|(i, v)| (v.clone(), i))
        .collect();
    let n_vars = var_list.len();

    let mut builder = SystemBuilder::new(tri, points, dim_v, origin_idx, &var_map, n_vars);
    builder.build_distinct_equations();
    builder.build_glsm_equations(glsm);

    let sol = solve_rectangular(&builder.m_mat, &builder.c_vec, n_vars)?;
    Ok(reduce_to_cy(n_pts, &var_map, &sol))
}

fn collect_variables(
    tri: &Triangulation,
    dim_v: usize,
    origin_idx: Option<usize>,
) -> Vec<Vec<usize>> {
    let mut variables = HashSet::new();
    for simplex in tri.simplices() {
        let rays: Vec<usize> = simplex
            .iter()
            .filter(|&&i| Some(i) != origin_idx)
            .copied()
            .collect();
        if rays.len() < dim_v {
            continue;
        }

        let perms = combinations_with_replacement(&rays, dim_v);
        for mut p in perms {
            p.sort_unstable();
            variables.insert(p);
        }
    }
    let mut var_list: Vec<Vec<usize>> = variables.into_iter().collect();
    var_list.sort_unstable();
    var_list
}

fn reduce_to_cy(
    n_pts: usize,
    var_map: &HashMap<Vec<usize>, usize>,
    sol: &[Rational],
) -> Intersection {
    let mut kappa = Intersection::new(n_pts);

    // We want to compute kappa_ijk = sum_l kappa_ijkl
    // We only need to check i,j,k that are subsets of some simplex
    let mut potential_triplets = HashSet::new();
    for key in var_map.keys() {
        // key is [a, b, c, d]
        potential_triplets.insert(canonical_key(key[0], key[1], key[2]));
        potential_triplets.insert(canonical_key(key[0], key[1], key[3]));
        potential_triplets.insert(canonical_key(key[0], key[2], key[3]));
        potential_triplets.insert(canonical_key(key[1], key[2], key[3]));
    }

    for (i, j, k) in potential_triplets {
        let mut sum = Rational::from(0);
        for l in 0..n_pts {
            let mut multiset = vec![i, j, k, l];
            multiset.sort_unstable();
            if let Some(&var_idx) = var_map.get(&multiset) {
                sum += &sol[var_idx];
            }
        }
        if sum != 0 {
            // Wrap in typed Rational<Finite> - intersection numbers can be any sign
            kappa.set(i, j, k, TypedRational::<Finite>::from_raw(sum));
        }
    }
    kappa
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
