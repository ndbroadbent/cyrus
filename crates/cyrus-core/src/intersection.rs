//! Intersection number computation.
//!
//! Implements the clean-room algorithm for computing intersection numbers
//! of divisors in a Calabi-Yau threefold.
//!
//! Algorithm:
//! 1. Identify distinct intersections from the triangulation (facets).
//! 2. Use GLSM linear relations to constrain the remaining intersections.
//! 3. Solve the linear system to find all entries of the intersection tensor.
//!
//! Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 1.2

use crate::Point;
use crate::error::{Error, Result};
use crate::f64_pos;
use crate::integer_math::determinant_gaussian;
use crate::triangulation::Triangulation;
use crate::types::f64::F64;
use crate::types::rational::Rational as TypedRational;
use crate::types::tags::{Finite, IsFinite, Pos};
use malachite::num::arithmetic::traits::Abs;
use malachite::{Integer, Rational};
use std::collections::{HashMap, HashSet};

/// Intersection tensor `κ_ijk`.
///
/// Stores intersection numbers for divisors `D_i, D_j, D_k`.
/// Symmetric in all indices. Values can be positive or negative.
#[derive(Debug, Clone)]
pub struct Intersection {
    /// Dimension of the divisor space (number of divisors).
    dim: usize,
    /// Storage for unique entries (i ≤ j ≤ k).
    /// Key: (i, j, k), Value: Intersection number (any finite value).
    entries: HashMap<(usize, usize, usize), TypedRational<Finite>>,
}

impl Intersection {
    /// Create a new empty intersection tensor.
    pub fn new(dim: usize) -> Self {
        Self {
            dim,
            entries: HashMap::new(),
        }
    }

    /// Get the dimension (number of divisors).
    pub const fn dim(&self) -> usize {
        self.dim
    }

    /// Get the number of non-zero entries.
    pub fn num_nonzero(&self) -> usize {
        self.entries.len()
    }

    /// Get value at indices (i, j, k). Returns zero if not set.
    pub fn get(&self, i: usize, j: usize, k: usize) -> TypedRational<Finite> {
        let key = canonical_key(i, j, k);
        self.entries
            .get(&key)
            .cloned()
            .unwrap_or_else(|| TypedRational::<Finite>::from_raw(Rational::from(0)))
    }

    /// Set value at indices (i, j, k).
    ///
    /// Accepts any rational type that can be widened to Finite (Pos, Neg, etc.).
    pub fn set<T>(&mut self, i: usize, j: usize, k: usize, val: TypedRational<T>)
    where
        TypedRational<T>: IsFinite<Finite = TypedRational<Finite>>,
    {
        let key = canonical_key(i, j, k);
        let val = val.to_finite();
        if *val.get() == 0 {
            self.entries.remove(&key);
        } else {
            self.entries.insert(key, val);
        }
    }

    /// Iterator over all non-zero entries.
    /// Returns ((i, j, k), val) where i ≤ j ≤ k.
    pub fn iter(&self) -> impl Iterator<Item = (&(usize, usize, usize), &TypedRational<Finite>)> {
        self.entries.iter()
    }

    /// Parallel iterator over non-zero entries.
    pub fn par_iter(
        &self,
    ) -> impl rayon::iter::ParallelIterator<Item = (&(usize, usize, usize), &TypedRational<Finite>)>
    where
        HashMap<(usize, usize, usize), TypedRational<Finite>>: rayon::iter::IntoParallelRefIterator<'static>,
    {
        use rayon::prelude::*;
        self.entries.par_iter()
    }

    /// Non-parallel iterator over non-zero entries.
    pub fn iter_entries(
        &self,
    ) -> std::collections::hash_map::Iter<'_, (usize, usize, usize), TypedRational<Finite>> {
        self.entries.iter()
    }

    /// Contract with positive moduli `t` to compute `κ_ijk t^i t^j t^k`.
    ///
    /// Input moduli must be positive (F64<Pos>) - for Kähler moduli.
    /// Returns `None` if the contraction is not positive (invalid physics).
    ///
    /// # Panics
    /// Panics if `t.len() != dim`.
    pub fn contract_triple(&self, t: &[F64<Pos>]) -> Option<F64<Pos>> {
        assert_eq!(t.len(), self.dim, "Vector dimension mismatch");

        // Accumulate as Finite (sum of terms could be any sign)
        // Type algebra: Pos * Finite = Finite, Finite + Finite = Finite
        let mut sum = F64::<Finite>::ZERO;
        for (&(i, j, k), val) in &self.entries {
            let mult = symmetry_multiplicity(i, j, k);
            let kappa = val.to_f64();
            // mult: Pos, kappa: Finite, t[i]: Pos → Pos * Finite * Pos * Pos * Pos = Finite
            let term = mult * kappa * t[i] * t[j] * t[k];
            sum = sum + term;
        }

        sum.try_to_pos()
    }

    /// Contract with finite-valued vector `t` to compute `κ_ijk t^i t^j t^k`.
    ///
    /// For cases where the vector can have any sign (e.g., flat directions).
    /// Returns `None` if the contraction is not positive.
    ///
    /// # Panics
    /// Panics if `t.len() != dim`.
    pub fn contract_triple_finite(&self, t: &[F64<Finite>]) -> Option<F64<Pos>> {
        assert_eq!(t.len(), self.dim, "Vector dimension mismatch");

        // Type algebra: Pos * Finite = Finite, Finite * Finite = Finite
        let mut sum = F64::<Finite>::ZERO;
        for (&(i, j, k), val) in &self.entries {
            let mult = symmetry_multiplicity(i, j, k);
            let kappa = val.to_f64();
            // mult: Pos, kappa: Finite, t[i]: Finite → all Finite
            let term = mult * kappa * t[i] * t[j] * t[k];
            sum = sum + term;
        }

        sum.try_to_pos()
    }

    /// Convert to NonEmptyIntersection if this tensor has entries.
    ///
    /// Returns `None` if the tensor is empty (no non-zero entries).
    #[must_use]
    pub fn into_non_empty(self) -> Option<NonEmptyIntersection> {
        if self.entries.is_empty() {
            None
        } else {
            Some(NonEmptyIntersection(self))
        }
    }
}

/// Non-empty intersection tensor.
///
/// Wrapper around `Intersection` that guarantees at least one non-zero entry.
/// Provides branded dimension for type-safe moduli handling.
#[derive(Debug, Clone)]
pub struct NonEmptyIntersection(Intersection);

impl NonEmptyIntersection {
    /// Get the dimension handle branded to this tensor.
    ///
    /// Use this to create `Moduli` vectors of the correct length.
    #[must_use]
    pub fn dim(&self) -> crate::types::branded::Dim<'_> {
        crate::types::branded::Dim(self.0.dim, std::marker::PhantomData)
    }

    /// Contract with validated moduli to compute `κ_ijk t^i t^j t^k`.
    ///
    /// # Panics
    /// Panics if the contraction is non-positive (moduli outside Kähler cone).
    #[must_use]
    pub fn contract_triple(&self, t: &crate::types::branded::Moduli<'_>) -> F64<Pos> {
        // Type algebra: accumulate as Finite, narrow at boundary
        let mut sum = F64::<Finite>::ZERO;
        for (&(i, j, k), val) in &self.0.entries {
            let mult = symmetry_multiplicity(i, j, k);
            let kappa = val.to_f64();
            let term = mult * kappa * t[i] * t[j] * t[k];
            sum = sum + term;
        }

        sum.try_to_pos()
            .expect("contraction must be positive - moduli outside Kähler cone")
    }

    /// Get the underlying intersection tensor.
    #[must_use]
    pub fn inner(&self) -> &Intersection {
        &self.0
    }
}

/// Compute canonical key (i ≤ j ≤ k) for symmetric tensor.
fn canonical_key(i: usize, j: usize, k: usize) -> (usize, usize, usize) {
    let mut v = [i, j, k];
    v.sort_unstable();
    v.into()
}

/// Compute symmetry multiplicity for entry (i, j, k) with i ≤ j ≤ k.
/// Returns F64<Pos> (1, 3, or 6 are all positive).
fn symmetry_multiplicity(i: usize, j: usize, k: usize) -> F64<Pos> {
    if i == j && j == k {
        f64_pos!(1.0)
    } else if i == j || j == k || i == k {
        f64_pos!(3.0)
    } else {
        f64_pos!(6.0)
    }
}

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
        potential_triplets.insert(canonical_key_2(key[0], key[1], key[2]));
        potential_triplets.insert(canonical_key_2(key[0], key[1], key[3]));
        potential_triplets.insert(canonical_key_2(key[0], key[2], key[3]));
        potential_triplets.insert(canonical_key_2(key[1], key[2], key[3]));
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

fn canonical_key_2(i: usize, j: usize, k: usize) -> (usize, usize, usize) {
    let mut v = [i, j, k];
    v.sort_unstable();
    v.into()
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

    #[test]
    fn test_set_get() {
        let mut kappa = Intersection::new(3);
        let val = TypedRational::<Finite>::from_raw(Rational::from(5));
        kappa.set(0, 1, 2, val.clone());

        assert_eq!(*kappa.get(0, 1, 2).get(), Rational::from(5));
        assert_eq!(*kappa.get(2, 1, 0).get(), Rational::from(5));
    }
}
