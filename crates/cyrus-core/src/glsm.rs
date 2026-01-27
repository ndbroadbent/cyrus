//! Geometric Langlands Super-Multiplet (GLSM) charge matrix computation.
//!
//! Implements the CYTools algorithm for constructing an integral GLSM basis
//! (with fallback to a non-integral basis) using Hermite normal form and LLL
//! ordering heuristics.

use crate::Point;
use crate::error::{Error, Result};
use crate::integer_math::{
    hermite_normal_form, is_row_hnf, smith_normal_form_diag, sublattice_index_snf,
};
use crate::utils::lll_reduce;
use malachite::Integer;
use malachite::num::arithmetic::traits::Abs;
use rand::seq::SliceRandom;
use rand::SeedableRng;

/// Compute the GLSM charge matrix for a set of lattice points.
///
/// This follows the CYTools algorithm for selecting an integral basis of
/// columns (when possible). The origin is always included in the computation;
/// it is removed from the output when `include_origin` is false.
pub fn compute_glsm_charge_matrix(
    points: &[Point],
    include_origin: bool,
) -> Result<Vec<Vec<Integer>>> {
    let (glsm, _linrel, _basis) = compute_glsm_and_linrels(points)?;
    if include_origin {
        Ok(glsm)
    } else {
        Ok(glsm.iter().map(|row| row[1..].to_vec()).collect())
    }
}

/// Compute the GLSM linear relations matrix.
pub fn compute_glsm_linear_relations(
    points: &[Point],
    include_origin: bool,
) -> Result<Vec<Vec<Integer>>> {
    let (_glsm, linrel, _basis) = compute_glsm_and_linrels(points)?;
    if include_origin {
        Ok(linrel)
    } else {
        if linrel.is_empty() {
            return Ok(Vec::new());
        }
        let mut trimmed = Vec::with_capacity(linrel.len().saturating_sub(1));
        for row in linrel.iter().skip(1) {
            trimmed.push(row[1..].to_vec());
        }
        Ok(trimmed)
    }
}

/// Compute GLSM charge matrix, linear relations, and basis indices.
pub fn compute_glsm_and_linrels(
    points: &[Point],
) -> Result<(Vec<Vec<Integer>>, Vec<Vec<Integer>>, Vec<usize>)> {
    let (pts, _origin_idx) = ensure_origin_first(points)?;
    let dim = pts[0].dim();
    let n_pts = pts.len();
    eprintln!("[DEBUG] glsm dim={}, n_pts={}", dim, n_pts);

    // Build linrel = points.T (dim x n)
    let mut linrel: Vec<Vec<Integer>> = vec![vec![Integer::from(0); n_pts]; dim];
    for (j, pt) in pts.iter().enumerate() {
        for (i, &coord) in pt.coords().iter().enumerate() {
            linrel[i][j] = Integer::from(coord);
        }
    }

    // Sublattice index via HNF pivot product of the column lattice.
    let sublat_ind = sublattice_index_snf(&linrel);
    let debug = std::env::var_os("CYRUS_GLSM_DEBUG").is_some();
    let expected_basis = parse_expected_basis();
    let pivot_mode = std::env::var("CYRUS_GLSM_PIVOT_MODE").unwrap_or_else(|_| "first".to_string());
    let scan_expected = expected_basis.is_some();
    if debug {
        let snf_diag = smith_normal_form_diag(&linrel);
        let mut snf_prod = Integer::from(1);
        for d in &snf_diag {
            snf_prod *= d.clone();
        }
        eprintln!(
            "[DEBUG] glsm sublat_ind: hnf={}, snf_prod={}",
            sublat_ind, snf_prod
        );
    }

    // Column norms (L1)
    let mut norms: Vec<(usize, i64)> = Vec::with_capacity(n_pts);
    for col in 0..n_pts {
        let mut acc: i64 = 0;
        for row in &linrel {
            let v = row[col].clone().abs();
            let vi = i64::try_from(&v).map_err(|_| {
                Error::InvalidInput("GLSM norm overflow".into())
            })?;
            acc = acc.saturating_add(vi);
        }
        norms.push((col, acc));
    }
    norms.sort_unstable_by(|a, b| a.1.cmp(&b.1));
    let mut indices: Vec<usize> = norms.iter().map(|&(i, _)| i).collect();

    // Build augmented linrel with ones row: (dim+1) x n
    let mut linrel_aug: Vec<Vec<Integer>> = Vec::with_capacity(dim + 1);
    linrel_aug.push(vec![Integer::from(1); n_pts]);
    linrel_aug.extend(linrel.clone());

    // Ensure first dim+1 indices are sorted
    indices[..linrel_aug.len()].sort_unstable();

    let mut found_good_basis = false;
    let mut first_basis_recorded = false;
    let mut expected_found = false;
    let mut chosen_indices: Vec<usize> = Vec::new();
    let mut chosen_linrel_rand: Vec<Vec<Integer>> = Vec::new();
    let mut chosen_basis_exc: Vec<usize> = Vec::new();
    let mut good_exclusions = 0usize;

    let mut rng = rand::rngs::StdRng::seed_from_u64(1337);
    let n_tries = 14;
    let max_ctr = linrel_aug.len() * n_pts + 1;

    for n_try in 0..n_tries {
        match n_try {
            1 => {
                indices = (0..n_pts).collect();
            }
            2 => {
                // LLL ordering heuristic (approximate CYTools behavior)
                let mut linrel_rows: Vec<Vec<i64>> = Vec::with_capacity(dim);
                for row in &linrel {
                    let mut r = Vec::with_capacity(n_pts);
                    for v in row {
                        r.push(i64::try_from(v).map_err(|_| {
                            Error::InvalidInput("GLSM LLL overflow".into())
                        })?);
                    }
                    linrel_rows.push(r);
                }
                let lll = lll_reduce(&linrel_rows, false);
                let mut lll_norms: Vec<(usize, i64)> = Vec::with_capacity(n_pts);
                for col in 0..n_pts {
                    let mut acc: i64 = 0;
                    for row in &lll.reduced {
                        acc = acc.saturating_add(row[col].abs());
                    }
                    lll_norms.push((col, acc));
                }
                lll_norms.sort_unstable_by(|a, b| a.1.cmp(&b.1));
                indices = lll_norms.iter().map(|&(i, _)| i).collect();
            }
            3 => {
                indices = std::iter::once(0)
                    .chain((1..n_pts).rev())
                    .collect();
            }
            _ if n_try > 3 => {
                indices[1..].shuffle(&mut rng);
            }
            _ => {}
        }

        indices[..linrel_aug.len()].sort_unstable();

        for _ctr in 0..max_ctr {
            // CYTools rotates on every iteration (ctr is incremented before the check).
            let st = good_exclusions.max(1);
            indices[st..].rotate_left(1);
            indices[..linrel_aug.len()].sort_unstable();

            // Permute columns
            let mut linrel_rand = vec![vec![Integer::from(0); n_pts]; linrel_aug.len()];
            for (new_col, &old_col) in indices.iter().enumerate() {
                for r in 0..linrel_aug.len() {
                    linrel_rand[r][new_col] = linrel_aug[r][old_col].clone();
                }
            }

            let linrel_hnf = hermite_normal_form(&linrel_rand);
            if debug && !is_row_hnf(&linrel_hnf) {
                eprintln!("[DEBUG] HNF validation failed for current permutation");
            }

            let mut tmp_sublat = Integer::from(1);
            good_exclusions = 0;
            let mut basis_exc_perm: Vec<usize> = Vec::new();
            found_good_basis = true;

            for row in &linrel_hnf {
                let pivot = if pivot_mode == "last" {
                    row.iter()
                        .rposition(|v| *v != 0)
                        .map(|i| (i, &row[i]))
                } else {
                    row.iter()
                        .enumerate()
                        .find(|(_, v)| **v != 0)
                        .map(|(i, v)| (i, v))
                };
                if let Some((i, val)) = pivot {
                    tmp_sublat *= val.clone().abs();
                    if (&sublat_ind % &tmp_sublat) == 0 {
                        good_exclusions += 1;
                        basis_exc_perm.push(i);
                    } else {
                        found_good_basis = false;
                        break;
                    }
                }
            }

            if !found_good_basis {
                continue;
            }

            let mut linrel_dict = vec![0usize; n_pts];
            for (pos, &idx) in indices.iter().enumerate() {
                linrel_dict[idx] = pos;
            }
            let basis_candidate: Vec<usize> = (0..n_pts)
                .filter(|&i| !basis_exc_perm.contains(&linrel_dict[i]))
                .collect();

            if debug {
                eprintln!(
                    "[DEBUG] glsm basis candidate: len={}, n_try={}",
                    basis_candidate.len(),
                    n_try
                );
                eprintln!("[DEBUG] glsm indices: {:?}", indices);
                eprintln!("[DEBUG] glsm basis_exc (perm): {:?}", basis_exc_perm);
                eprintln!("[DEBUG] glsm basis (orig): {:?}", basis_candidate);
            }

            if !first_basis_recorded {
                chosen_indices = indices.clone();
                chosen_linrel_rand = linrel_hnf.clone();
                chosen_basis_exc = basis_exc_perm.clone();
                first_basis_recorded = true;
            }

            if let Some(expected) = expected_basis.as_ref() {
                if basis_candidate == *expected {
                    expected_found = true;
                    eprintln!(
                        "[DEBUG] glsm expected basis found at n_try={}, ctr={}",
                        n_try, _ctr
                    );
                }
            }

            if !scan_expected {
                break;
            }
        }

        if found_good_basis && !scan_expected {
            break;
        }
    }

    if scan_expected && !expected_found {
        eprintln!("[DEBUG] glsm expected basis not found in scan");
    }

    if !first_basis_recorded {
        return Err(Error::InvalidInput(
            "Integral GLSM basis not found".into(),
        ));
    };

    // Map columns back to original order.
    let mut linrel_dict = vec![0usize; n_pts];
    for (pos, &idx) in chosen_indices.iter().enumerate() {
        linrel_dict[idx] = pos;
    }
    let mut linrel_reordered = vec![vec![Integer::from(0); n_pts]; chosen_linrel_rand.len()];
    for orig_col in 0..n_pts {
        let pos = linrel_dict[orig_col];
        for r in 0..chosen_linrel_rand.len() {
            linrel_reordered[r][orig_col] = chosen_linrel_rand[r][pos].clone();
        }
    }

    let basis_exc: Vec<usize> = chosen_basis_exc
        .iter()
        .map(|&i| chosen_indices[i])
        .collect();

    let basis: Vec<usize> = (0..n_pts).filter(|&i| !basis_exc.contains(&i)).collect();

    let glsm_rows = n_pts - linrel_reordered.len();
    let mut glsm = vec![vec![Integer::from(0); n_pts]; glsm_rows];
    for (r, &col) in basis.iter().enumerate() {
        glsm[r][col] = Integer::from(1);
    }

    for &nb in basis_exc.iter().rev() {
        let mut last: Option<(usize, Integer)> = None;
        for (i, row) in linrel_reordered.iter().enumerate() {
            if row[nb] != 0 {
                last = Some((i, row[nb].clone()));
            }
        }
        let Some((i, ii)) = last else {
            continue;
        };
        if (&sublat_ind % &ii) != 0 {
            return Err(Error::InvalidInput(format!(
                "Problem with linear relations (sublat_ind={}, ii={})",
                sublat_ind, ii
            )));
        }
        for r in 0..glsm_rows {
            let mut acc = Integer::from(0);
            for c in 0..n_pts {
                acc += &glsm[r][c] * &linrel_reordered[i][c];
            }
            let mut val = -acc;
            val /= ii.clone();
            glsm[r][nb] = val;
        }
    }

    // Validate GLSM relations (CYTools consistency checks).
    if !glsm_relations_ok(&glsm, &linrel_reordered, &pts) {
        return Err(Error::InvalidInput("Error finding GLSM basis".into()));
    }

    Ok((glsm, linrel_reordered, basis))
}

fn glsm_relations_ok(
    glsm: &[Vec<Integer>],
    linrel: &[Vec<Integer>],
    points: &[Point],
) -> bool {
    if glsm.is_empty() {
        return false;
    }

    // Check glsm * linrel^T == 0
    let linrel_rows = linrel.len();
    let n_cols = linrel[0].len();
    for row in glsm {
        for lr in 0..linrel_rows {
            let mut acc = Integer::from(0);
            for c in 0..n_cols {
                acc += &row[c] * &linrel[lr][c];
            }
            if acc != 0 {
                return false;
            }
        }
    }

    // Check glsm * points == 0 (points matrix n x dim).
    let dim = points[0].dim();
    for row in glsm {
        for d in 0..dim {
            let mut acc = Integer::from(0);
            for (c, p) in points.iter().enumerate() {
                acc += &row[c] * Integer::from(p.coords()[d]);
            }
            if acc != 0 {
                return false;
            }
        }
    }

    true
}

fn ensure_origin_first(points: &[Point]) -> Result<(Vec<Point>, usize)> {
    if points.is_empty() {
        return Err(Error::InvalidInput("No points provided".into()));
    }
    let dim = points[0].dim();
    let origin_idx = points
        .iter()
        .position(|p| p.coords().iter().all(|&x| x == 0))
        .unwrap_or(usize::MAX);

    if origin_idx == 0 {
        return Ok((points.to_vec(), 0));
    }

    let mut pts = Vec::with_capacity(points.len() + if origin_idx == usize::MAX { 1 } else { 0 });
    pts.push(Point::origin(dim));
    for (idx, p) in points.iter().enumerate() {
        if idx == origin_idx {
            continue;
        }
        pts.push(p.clone());
    }
    Ok((pts, 0))
}

fn parse_expected_basis() -> Option<Vec<usize>> {
    let raw = std::env::var("CYRUS_GLSM_EXPECT_BASIS").ok()?;
    let path = std::path::PathBuf::from(&raw);
    let text = if path.exists() {
        std::fs::read_to_string(&path).ok()?
    } else {
        raw
    };
    let mut out = Vec::new();
    for part in text.split(|c| c == ',' || c == '\n' || c == '\r') {
        let trimmed = part.trim();
        if trimmed.is_empty() {
            continue;
        }
        let val: usize = trimmed.parse().ok()?;
        out.push(val);
    }
    if out.is_empty() {
        None
    } else {
        Some(out)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use malachite::Integer;

    fn i(n: i64) -> Integer {
        Integer::from(n)
    }

    #[test]
    fn test_glsm_p2() {
        // P^2 has vertices (1,0), (0,1), (-1,-1)
        let points = vec![
            Point::new(vec![1, 0]),
            Point::new(vec![0, 1]),
            Point::new(vec![-1, -1]),
        ];

        let q = compute_glsm_charge_matrix(&points, true).unwrap();

        // Expected: 1 relation for P^2
        assert_eq!(q.len(), 1);

        let row = &q[0];
        assert_eq!(row.len(), 4);

        // Verify relation Σ Qi v_i = 0
        let mut x = i(0);
        let mut y = i(0);
        let mut sum_q = i(0);

        // Origin
        sum_q += &row[0];
        // Vertices
        x += &row[1] * i(1);
        y += &row[1] * i(0);
        sum_q += &row[1];

        x += &row[2] * i(0);
        y += &row[2] * i(1);
        sum_q += &row[2];

        x += &row[3] * i(-1);
        y += &row[3] * i(-1);
        sum_q += &row[3];

        assert_eq!(x, i(0));
        assert_eq!(y, i(0));
        assert_eq!(sum_q, i(0));
    }

    #[test]
    fn test_glsm_p1xp1() {
        // P^1 x P^1 has 4 rays: (1,0), (-1,0), (0,1), (0,-1)
        let points = vec![
            Point::new(vec![1, 0]),
            Point::new(vec![-1, 0]),
            Point::new(vec![0, 1]),
            Point::new(vec![0, -1]),
        ];

        let q = compute_glsm_charge_matrix(&points, true).unwrap();

        // Expected: 2 relations for P1xP1
        assert_eq!(q.len(), 2);

        for row in &q {
            assert_eq!(row.len(), 5);
            let mut x = i(0);
            let mut y = i(0);
            let mut sum_q = i(0);

            // Origin
            sum_q += &row[0];
            // Vertices
            x += &row[1] * i(1);
            sum_q += &row[1];
            x += &row[2] * i(-1);
            sum_q += &row[2];
            y += &row[3] * i(1);
            sum_q += &row[3];
            y += &row[4] * i(-1);
            sum_q += &row[4];

            assert_eq!(x, i(0));
            assert_eq!(y, i(0));
            assert_eq!(sum_q, i(0));
        }
    }
}
