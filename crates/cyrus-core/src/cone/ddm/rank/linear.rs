use malachite::{Integer, Rational};

use crate::integer_math::gcd_integer;

use super::super::types::DdmHyperplane;

pub(in crate::cone::ddm) const RANK_MODULUS: i128 = 2_147_483_647;

pub(super) fn active_rank_integer_at_least(
    active: &[usize],
    hyperplanes: &[DdmHyperplane],
    dim: usize,
    target_rank: usize,
) -> bool {
    let mut work: Vec<Vec<Integer>> = active
        .iter()
        .map(|&idx| {
            hyperplanes[idx]
                .dense()
                .iter()
                .map(|&x| Integer::from(x))
                .collect()
        })
        .collect();
    integer_rank_at_least(&mut work, dim, target_rank)
}

pub(super) fn integer_rank_at_least(
    work: &mut [Vec<Integer>],
    dim: usize,
    target_rank: usize,
) -> bool {
    let mut rank = 0usize;

    for col in 0..dim {
        let Some(pivot_row) = find_integer_pivot(work, rank, col) else {
            continue;
        };
        work.swap(rank, pivot_row);
        let pivot = work[rank][col].clone();
        eliminate_integer_column(work, rank, col, dim, &pivot);
        rank += 1;
        if rank >= target_rank {
            return true;
        }
    }

    false
}

fn find_integer_pivot(work: &[Vec<Integer>], start_row: usize, col: usize) -> Option<usize> {
    (start_row..work.len()).find(|&row| work[row][col] != 0)
}

fn eliminate_integer_column(
    work: &mut [Vec<Integer>],
    pivot_row: usize,
    col: usize,
    dim: usize,
    pivot: &Integer,
) {
    let pivot_values = work[pivot_row].clone();
    for row in (pivot_row + 1)..work.len() {
        eliminate_integer_row(work, row, col, dim, pivot, &pivot_values);
    }
}

fn eliminate_integer_row(
    work: &mut [Vec<Integer>],
    row: usize,
    col: usize,
    dim: usize,
    pivot: &Integer,
    pivot_values: &[Integer],
) {
    let row_factor = work[row][col].clone();
    if row_factor == 0 {
        return;
    }
    for c in (col + 1)..dim {
        work[row][c] = &work[row][c] * pivot - &row_factor * &pivot_values[c];
    }
    work[row][col] = Integer::from(0);
    normalize_integer_row(&mut work[row], col + 1);
}

fn normalize_integer_row(row: &mut [Integer], start_col: usize) {
    let content = row_content_gcd(row, start_col);
    if content <= Integer::from(1) {
        return;
    }
    for value in row.iter_mut().skip(start_col) {
        *value /= &content;
    }
}

fn row_content_gcd(row: &[Integer], start_col: usize) -> Integer {
    row.iter()
        .skip(start_col)
        .fold(Integer::from(0), |acc, value| gcd_integer(&acc, value))
}

pub(in crate::cone::ddm) fn active_rank_mod_prime_at_least(
    active: &[usize],
    hyperplanes: &[DdmHyperplane],
    dim: usize,
    target_rank: usize,
) -> bool {
    let mut work: Vec<Vec<i128>> = active
        .iter()
        .map(|&idx| {
            hyperplanes[idx]
                .dense()
                .iter()
                .map(|&x| x.rem_euclid(RANK_MODULUS))
                .collect()
        })
        .collect();

    mod_rank_at_least(&mut work, dim, target_rank)
}

#[cfg(test)]
pub(in crate::cone::ddm) fn quotient_rank_mod_prime_at_least(
    work: &[Vec<Integer>],
    dim: usize,
    target_rank: usize,
) -> bool {
    let mut mod_work: Vec<Vec<i128>> = work
        .iter()
        .map(|row| row.iter().map(integer_mod_prime).collect())
        .collect();
    mod_rank_at_least(&mut mod_work, dim, target_rank)
}

pub(in crate::cone::ddm) fn integer_mod_prime(value: &Integer) -> i128 {
    let modulus = Integer::from(RANK_MODULUS);
    let mut residue = value % &modulus;
    if residue < 0 {
        residue += modulus;
    }
    i128::try_from(&residue).expect("residue is smaller than RANK_MODULUS")
}

pub(super) fn mod_rank_at_least(work: &mut [Vec<i128>], dim: usize, target_rank: usize) -> bool {
    let mut rank = 0usize;
    for col in 0..dim {
        let Some(pivot) = find_mod_pivot(work, rank, col) else {
            continue;
        };
        work.swap(rank, pivot);
        normalize_mod_pivot_row(&mut work[rank], col);
        eliminate_mod_column(work, rank, col, dim);
        rank += 1;
        if rank >= target_rank {
            return true;
        }
    }
    false
}

fn find_mod_pivot(work: &[Vec<i128>], start_row: usize, col: usize) -> Option<usize> {
    (start_row..work.len()).find(|&row| work[row][col] != 0)
}

fn normalize_mod_pivot_row(row: &mut [i128], col: usize) {
    let inv = mod_inverse_prime(row[col]);
    for value in row.iter_mut().skip(col) {
        *value = (*value * inv).rem_euclid(RANK_MODULUS);
    }
}

fn eliminate_mod_column(work: &mut [Vec<i128>], pivot_row: usize, col: usize, dim: usize) {
    for row in (pivot_row + 1)..work.len() {
        let factor = work[row][col];
        if factor == 0 {
            continue;
        }
        for c in col..dim {
            work[row][c] = (work[row][c] - factor * work[pivot_row][c]).rem_euclid(RANK_MODULUS);
        }
    }
}

const fn mod_inverse_prime(value: i128) -> i128 {
    mod_pow(value, RANK_MODULUS - 2)
}

const fn mod_pow(mut base: i128, mut exponent: i128) -> i128 {
    let mut out = 1i128;
    base = base.rem_euclid(RANK_MODULUS);
    while exponent > 0 {
        if exponent % 2 == 1 {
            out = (out * base).rem_euclid(RANK_MODULUS);
        }
        base = (base * base).rem_euclid(RANK_MODULUS);
        exponent /= 2;
    }
    out
}

pub(in crate::cone::ddm) struct RankTracker {
    dim: usize,
    pivot_cols: Vec<usize>,
    basis: Vec<Vec<Rational>>,
}

impl RankTracker {
    pub(in crate::cone::ddm) const fn new(dim: usize) -> Self {
        Self {
            dim,
            pivot_cols: Vec::new(),
            basis: Vec::new(),
        }
    }

    pub(in crate::cone::ddm) const fn rank(&self) -> usize {
        self.basis.len()
    }

    pub(in crate::cone::ddm) fn add_row(&mut self, row: &[i128]) {
        let mut work: Vec<Rational> = row
            .iter()
            .map(|&x| Rational::from(Integer::from(x)))
            .collect();

        for (basis_row, &pivot_col) in self.basis.iter().zip(self.pivot_cols.iter()) {
            if work[pivot_col] == 0 {
                continue;
            }
            let factor = work[pivot_col].clone();
            for col in pivot_col..self.dim {
                let sub = &factor * &basis_row[col];
                work[col] -= sub;
            }
        }

        let Some(pivot_col) = work.iter().position(|value| *value != 0) else {
            return;
        };
        let pivot = work[pivot_col].clone();
        for value in work.iter_mut().skip(pivot_col) {
            *value /= &pivot;
        }

        let insert_at = self
            .pivot_cols
            .binary_search(&pivot_col)
            .expect_err("new rank pivot should not duplicate existing pivot");
        self.pivot_cols.insert(insert_at, pivot_col);
        self.basis.insert(insert_at, work);
    }
}
