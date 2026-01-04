//! Intersection number computation algorithm.
//!
//! Uses sparse linear algebra for efficiency on large triangulations.
//!
//! Algorithm (optimized):
//! 1. Enumerate variables: Only 4-tuples from simplex faces (with replacement)
//! 2. Build sparse system: Distinct equations + GLSM from simplex 3-faces only
//! 3. Solve: Normal equations (M^T M x = M^T C) with Cholesky
//! 4. Reduce: Sum over 4th index to get κ_ijk
//!
//! Reference: [[project_docs/clean_room/INTERSECTION_NUMBERS.md]]

use super::{Intersection, canonical_key};
use crate::Point;
use crate::error::{Error, Result};
use crate::integer_math::determinant_gaussian;
use crate::triangulation::Triangulation;
use crate::types::rational::Rational as TypedRational;
use crate::types::tags::Finite;
use faer::prelude::SpSolver;
use faer::sparse::SparseColMat;
use malachite::num::conversion::traits::RoundingFrom;
use malachite::rounding_modes::RoundingMode;
use malachite::{Integer, Rational};
use std::collections::{HashMap, HashSet};

/// Compute the intersection numbers for a Calabi-Yau hypersurface.
///
/// Uses sparse linear algebra for efficient computation on large triangulations.
///
/// # Errors
/// Returns an error if the linear system for intersection numbers is unsolvable or singular.
///
/// Reference: [[project_docs/clean_room/INTERSECTION_NUMBERS.md]]
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

    // Step 1: Collect simplex rays (excluding origin) and 3-face probes
    let (simplex_rays_list, probes) = collect_simplices_and_probes(tri, dim_v, origin_idx);

    // Step 2: Enumerate variables (4-tuples with replacement from simplex rays)
    let var_list = collect_variables(&simplex_rays_list, dim_v);
    if var_list.is_empty() {
        return Ok(Intersection::new(n_pts));
    }

    let var_map: HashMap<Vec<usize>, usize> = var_list
        .iter()
        .enumerate()
        .map(|(i, v)| (v.clone(), i))
        .collect();
    let n_vars = var_list.len();

    // Step 3: Build sparse system
    let (triplets, c_vec) = build_sparse_system(
        &simplex_rays_list,
        &probes,
        points,
        dim_v,
        origin_idx,
        &var_map,
        glsm,
    );

    let n_equations = c_vec.len();
    if n_equations == 0 {
        return Ok(Intersection::new(n_pts));
    }

    // Step 4: Solve using normal equations with Cholesky
    let sol = solve_normal_equations(&triplets, &c_vec, n_equations, n_vars)?;

    // Step 5: Reduce to CY intersection numbers κ_ijk
    Ok(reduce_to_cy(n_pts, &var_map, &sol))
}

/// Collect simplex rays and unique 3-face probes from the triangulation.
fn collect_simplices_and_probes(
    tri: &Triangulation,
    dim_v: usize,
    origin_idx: Option<usize>,
) -> (Vec<Vec<usize>>, HashSet<Vec<usize>>) {
    let mut simplex_rays_list = Vec::new();
    let mut probes = HashSet::new();

    for simplex in tri.simplices() {
        let mut rays: Vec<usize> = simplex
            .iter()
            .filter(|&&i| Some(i) != origin_idx)
            .copied()
            .collect();

        if rays.len() < dim_v {
            continue;
        }

        rays.sort_unstable();

        // Generate (dim_v - 1)-face probes: all (dim_v-1)-subsets of rays
        for i in 0..rays.len() {
            let mut probe = rays.clone();
            probe.remove(i);
            if probe.len() == dim_v - 1 {
                probes.insert(probe);
            }
        }

        simplex_rays_list.push(rays);
    }

    (simplex_rays_list, probes)
}

/// Collect variables (4-tuples with replacement) from simplex rays.
fn collect_variables(simplex_rays_list: &[Vec<usize>], dim_v: usize) -> Vec<Vec<usize>> {
    let mut variables = HashSet::new();

    for rays in simplex_rays_list {
        // Generate all dim_v-tuples with replacement from the rays of this simplex
        let perms = combinations_with_replacement(rays, dim_v);
        for mut p in perms {
            p.sort_unstable();
            variables.insert(p);
        }
    }

    let mut var_list: Vec<Vec<usize>> = variables.into_iter().collect();
    var_list.sort_unstable();
    var_list
}

/// Build the sparse system: distinct equations + GLSM equations.
/// Returns (triplets, rhs_vector).
#[allow(clippy::too_many_arguments)]
fn build_sparse_system(
    simplex_rays_list: &[Vec<usize>],
    probes: &HashSet<Vec<usize>>,
    points: &[Point],
    dim_v: usize,
    origin_idx: Option<usize>,
    var_map: &HashMap<Vec<usize>, usize>,
    glsm: &[Vec<Integer>],
) -> (Vec<(usize, usize, f64)>, Vec<f64>) {
    let mut triplets = Vec::new();
    let mut c_vec = Vec::new();
    let mut row_idx = 0;

    // Distinct intersection equations: for each simplex with dim_v rays
    for rays in simplex_rays_list {
        if rays.len() == dim_v {
            let mut mat_pts: Vec<Vec<Rational>> = rays
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
                let abs_det = if det < 0 { -det.clone() } else { det.clone() };
                let (det_f64, _) = f64::rounding_from(&abs_det, RoundingMode::Nearest);
                let val = 1.0 / det_f64;

                if let Some(&var_idx) = var_map.get(rays) {
                    triplets.push((row_idx, var_idx, 1.0));
                    c_vec.push(val);
                    row_idx += 1;
                }
            }
        }
    }

    // GLSM equations: only from simplex 3-face probes
    for probe in probes {
        for q_row in glsm {
            let mut row_triplets = Vec::new();

            for (m, q_val) in q_row.iter().enumerate() {
                if *q_val == 0 || Some(m) == origin_idx {
                    continue;
                }

                let mut key = probe.clone();
                key.push(m);
                key.sort_unstable();

                if let Some(&var_idx) = var_map.get(&key) {
                    // Convert Integer to f64
                    let coeff: f64 = q_val
                        .to_string()
                        .parse()
                        .unwrap_or(0.0);
                    if coeff != 0.0 {
                        row_triplets.push((var_idx, coeff));
                    }
                }
            }

            if !row_triplets.is_empty() {
                for (var_idx, coeff) in row_triplets {
                    triplets.push((row_idx, var_idx, coeff));
                }
                c_vec.push(0.0);
                row_idx += 1;
            }
        }
    }

    (triplets, c_vec)
}

/// Solve using normal equations: M^T M x = M^T C with sparse Cholesky.
fn solve_normal_equations(
    triplets: &[(usize, usize, f64)],
    c_vec: &[f64],
    _n_rows: usize,
    n_cols: usize,
) -> Result<Vec<f64>> {
    // Group triplets by row for efficient M^T M computation
    let mut rows_data: HashMap<usize, Vec<(usize, f64)>> = HashMap::new();
    for &(row, col, val) in triplets {
        rows_data.entry(row).or_default().push((col, val));
    }

    // Build M^T M as sparse triplets
    // (M^T M)[i,j] = sum over rows k: M[k,i] * M[k,j]
    let mut mtm_map: HashMap<(usize, usize), f64> = HashMap::new();
    for entries in rows_data.values() {
        for &(col_i, val_i) in entries {
            for &(col_j, val_j) in entries {
                *mtm_map.entry((col_i, col_j)).or_insert(0.0) += val_i * val_j;
            }
        }
    }

    // Add regularization to diagonal
    for i in 0..n_cols {
        *mtm_map.entry((i, i)).or_insert(0.0) += 1e-12;
    }

    // Convert to triplets for sparse matrix construction
    let mtm_triplets: Vec<(usize, usize, f64)> = mtm_map
        .into_iter()
        .map(|((i, j), v)| (i, j, v))
        .collect();

    // Build sparse M^T M matrix
    let mtm_sparse = SparseColMat::<usize, f64>::try_new_from_triplets(n_cols, n_cols, &mtm_triplets)
        .map_err(|e| Error::SingularMatrix(format!("Failed to build sparse M^T M: {e:?}")))?;

    // Build M^T C as dense vector
    let mut mtc = faer::Mat::<f64>::zeros(n_cols, 1);
    for &(row, col, val) in triplets {
        mtc[(col, 0)] += val * c_vec[row];
    }

    // Solve using sparse Cholesky
    let chol = mtm_sparse.sp_cholesky(faer::Side::Lower)
        .map_err(|_| Error::SingularMatrix("Sparse Cholesky factorization failed".into()))?;

    let solution = chol.solve(&mtc);

    // Extract solution
    let sol: Vec<f64> = (0..n_cols).map(|i| solution[(i, 0)]).collect();

    if sol.iter().any(|&x| !x.is_finite()) {
        return Err(Error::SingularMatrix("Solution contains non-finite values".into()));
    }

    Ok(sol)
}

/// Reduce 4-form intersection numbers to CY 3-form.
/// κ_ijk = sum_l κ_ijkl
fn reduce_to_cy(
    n_pts: usize,
    var_map: &HashMap<Vec<usize>, usize>,
    sol: &[f64],
) -> Intersection {
    let mut kappa = Intersection::new(n_pts);

    // Collect potential triplets from the variables
    let mut potential_triplets = HashSet::new();
    for key in var_map.keys() {
        // key is [a, b, c, d] - generate all 3-subsets
        potential_triplets.insert(canonical_key(key[0], key[1], key[2]));
        potential_triplets.insert(canonical_key(key[0], key[1], key[3]));
        potential_triplets.insert(canonical_key(key[0], key[2], key[3]));
        potential_triplets.insert(canonical_key(key[1], key[2], key[3]));
    }

    for (i, j, k) in potential_triplets {
        let mut sum = 0.0_f64;

        for l in 0..n_pts {
            let mut multiset = vec![i, j, k, l];
            multiset.sort_unstable();

            if let Some(&var_idx) = var_map.get(&multiset) {
                sum += sol[var_idx];
            }
        }

        // Round to nearest rational if significant
        if sum.abs() > 1e-10 {
            let rational = float_to_rational(sum);
            kappa.set(i, j, k, TypedRational::<Finite>::from_raw(rational));
        }
    }

    kappa
}

/// Convert f64 to a rational number.
fn float_to_rational(x: f64) -> Rational {
    // First try: is it close to an integer?
    let rounded = x.round();
    if (x - rounded).abs() < 1e-9 {
        #[allow(clippy::cast_possible_truncation)]
        return Rational::from(rounded as i64);
    }

    // Try common denominators
    for denom in [2, 3, 4, 5, 6, 8, 10, 12] {
        let scaled = x * f64::from(denom);
        let rounded_scaled = scaled.round();
        if (scaled - rounded_scaled).abs() < 1e-9 {
            #[allow(clippy::cast_possible_truncation)]
            return Rational::from(rounded_scaled as i64) / Rational::from(denom);
        }
    }

    // Continued fraction approximation
    continued_fraction_approx(x, 1000)
}

/// Approximate a float using continued fractions.
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

/// Generate combinations with replacement of size k from the pool.
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
