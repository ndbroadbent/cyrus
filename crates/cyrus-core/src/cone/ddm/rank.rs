use std::collections::HashMap;

use malachite::{Integer, Rational};

use crate::integer_math::lcm_integer;

mod linear;

#[cfg(test)]
pub(super) use linear::quotient_rank_mod_prime_at_least;
pub(super) use linear::{
    RANK_MODULUS, RankTracker, active_rank_mod_prime_at_least, integer_mod_prime,
};
use linear::{active_rank_integer_at_least, integer_rank_at_least, mod_rank_at_least};

use super::types::{DdmHyperplane, DdmRay, RankStats};
#[cfg(test)]
use super::util::hyperplanes_from_dense;
use super::util::intersect_active_if_len_at_least;

#[cfg(test)]
pub(super) fn are_adjacent(p: &DdmRay, n: &DdmRay, hyperplanes: &[Vec<i128>], dim: usize) -> bool {
    let hyperplanes = hyperplanes_from_dense(hyperplanes.to_vec());
    let target_rank = dim.saturating_sub(2);
    if target_rank == 0 {
        return true;
    }
    let Some(common) = intersect_active_if_len_at_least(&p.active, &n.active, target_rank) else {
        return false;
    };
    let mut stats = RankStats::default();
    let mut rank_context = RankContext::default();
    active_rank_at_least(
        &common,
        &hyperplanes,
        dim,
        target_rank,
        &mut rank_context,
        &mut stats,
    )
}

pub(super) fn are_adjacent_cached(
    p: &DdmRay,
    n: &DdmRay,
    hyperplanes: &[DdmHyperplane],
    dim: usize,
    rank_context: &mut RankContext,
    cache: &mut HashMap<Vec<usize>, bool>,
    rank_stats: &mut RankStats,
) -> bool {
    let target_rank = dim.saturating_sub(2);
    if target_rank == 0 {
        return true;
    }
    let Some(common) = intersect_active_if_len_at_least(&p.active, &n.active, target_rank) else {
        rank_stats.too_short += 1;
        return false;
    };
    if let Some(&value) = cache.get(&common) {
        rank_stats.cache_hits += 1;
        return value;
    }
    let value = active_rank_at_least(
        &common,
        hyperplanes,
        dim,
        target_rank,
        rank_context,
        rank_stats,
    );
    cache.insert(common, value);
    value
}

#[derive(Default)]
pub(super) struct RankContext {
    quotient_ranker: Option<BasisQuotientRanker>,
}

impl RankContext {
    pub(super) fn with_basis_inverse(inverse: &[Vec<Rational>]) -> Option<Self> {
        Some(Self {
            quotient_ranker: Some(BasisQuotientRanker::new(inverse)?),
        })
    }

    const fn has_basis(&self) -> bool {
        self.quotient_ranker.is_some()
    }

    fn quotient_rank_at_least(
        &mut self,
        active: &[usize],
        hyperplanes: &[DdmHyperplane],
        target_rank: usize,
    ) -> Option<bool> {
        self.quotient_ranker
            .as_mut()
            .map(|ranker| ranker.rank_at_least(active, hyperplanes, target_rank))
    }
}

struct BasisQuotientRanker {
    dim: usize,
    scaled_inverse: Vec<Vec<Integer>>,
    scaled_inverse_mod: Vec<Vec<i128>>,
    coordinates: HashMap<usize, Vec<Integer>>,
    coordinates_mod: HashMap<usize, Vec<i128>>,
}

impl BasisQuotientRanker {
    fn new(inverse: &[Vec<Rational>]) -> Option<Self> {
        let scaled_inverse = scaled_inverse_matrix(inverse)?;
        let scaled_inverse_mod = scaled_inverse
            .iter()
            .map(|row| row.iter().map(integer_mod_prime).collect())
            .collect();
        Some(Self {
            dim: inverse.len(),
            scaled_inverse,
            scaled_inverse_mod,
            coordinates: HashMap::new(),
            coordinates_mod: HashMap::new(),
        })
    }

    fn rank_at_least(
        &mut self,
        active: &[usize],
        hyperplanes: &[DdmHyperplane],
        target_rank: usize,
    ) -> bool {
        let (basis_count, missing_cols) = active_basis_quotient(active, self.dim);
        if basis_count >= target_rank {
            return true;
        }

        let residual_target = target_rank - basis_count;
        let nonbasis_count = active.iter().filter(|&&idx| idx >= self.dim).count();
        if nonbasis_count < residual_target {
            return false;
        }
        if residual_target == 1 {
            return self.quotient_rank_one_at_least(active, hyperplanes, &missing_cols);
        }
        if self.quotient_rank_mod_prime_at_least(
            active,
            hyperplanes,
            &missing_cols,
            residual_target,
        ) {
            return true;
        }
        if residual_target == 2 {
            return self.quotient_rank_two_at_least(active, hyperplanes, &missing_cols);
        }
        let mut work = self.integer_quotient_work(active, hyperplanes, &missing_cols);
        integer_rank_at_least(&mut work, missing_cols.len(), residual_target)
    }

    fn quotient_rank_one_at_least(
        &mut self,
        active: &[usize],
        hyperplanes: &[DdmHyperplane],
        missing_cols: &[usize],
    ) -> bool {
        active.iter().any(|&idx| {
            idx >= self.dim && {
                let coords = self.coordinates_for(idx, &hyperplanes[idx]);
                missing_cols.iter().any(|&col| coords[col] != 0)
            }
        })
    }

    fn quotient_rank_two_at_least(
        &mut self,
        active: &[usize],
        hyperplanes: &[DdmHyperplane],
        missing_cols: &[usize],
    ) -> bool {
        let mut first: Option<Vec<Integer>> = None;
        let mut pivot_col = 0usize;

        for &idx in active {
            if idx < self.dim {
                continue;
            }
            let coords = self.coordinates_for(idx, &hyperplanes[idx]);
            if let Some(first_row) = &first {
                if !projected_rows_are_proportional(first_row, pivot_col, coords, missing_cols) {
                    return true;
                }
            } else if let Some((pos, row)) = projected_nonzero_row(coords, missing_cols) {
                pivot_col = pos;
                first = Some(row);
            }
        }
        false
    }

    fn quotient_rank_mod_prime_at_least(
        &mut self,
        active: &[usize],
        hyperplanes: &[DdmHyperplane],
        missing_cols: &[usize],
        target_rank: usize,
    ) -> bool {
        let mut work = Vec::new();
        for &idx in active {
            if idx < self.dim {
                continue;
            }
            let coords = self.coordinates_mod_for(idx, &hyperplanes[idx]);
            work.push(missing_cols.iter().map(|&col| coords[col]).collect());
        }
        mod_rank_at_least(&mut work, missing_cols.len(), target_rank)
    }

    fn integer_quotient_work(
        &mut self,
        active: &[usize],
        hyperplanes: &[DdmHyperplane],
        missing_cols: &[usize],
    ) -> Vec<Vec<Integer>> {
        let mut work = Vec::new();
        for &idx in active {
            if idx < self.dim {
                continue;
            }
            let coords = self.coordinates_for(idx, &hyperplanes[idx]);
            work.push(
                missing_cols
                    .iter()
                    .map(|&col| coords[col].clone())
                    .collect(),
            );
        }
        work
    }

    fn coordinates_for(&mut self, row_idx: usize, row: &DdmHyperplane) -> &[Integer] {
        if !self.coordinates.contains_key(&row_idx) {
            let coordinates = self.compute_coordinates(row);
            self.coordinates.insert(row_idx, coordinates);
        }
        self.coordinates
            .get(&row_idx)
            .expect("coordinate cache was just populated")
    }

    fn coordinates_mod_for(&mut self, row_idx: usize, row: &DdmHyperplane) -> &[i128] {
        if !self.coordinates_mod.contains_key(&row_idx) {
            let coordinates = self.compute_coordinates_mod(row);
            self.coordinates_mod.insert(row_idx, coordinates);
        }
        self.coordinates_mod
            .get(&row_idx)
            .expect("modular coordinate cache was just populated")
    }

    fn compute_coordinates(&self, row: &DdmHyperplane) -> Vec<Integer> {
        let mut out = vec![Integer::from(0); self.dim];
        for &(source_col, coefficient) in &row.sparse {
            let factor = Integer::from(coefficient);
            for (target_col, out_value) in out.iter_mut().enumerate() {
                let term = &factor * &self.scaled_inverse[source_col][target_col];
                *out_value += term;
            }
        }
        out
    }

    fn compute_coordinates_mod(&self, row: &DdmHyperplane) -> Vec<i128> {
        let mut out = vec![0i128; self.dim];
        for &(source_col, coefficient) in &row.sparse {
            let factor = coefficient.rem_euclid(RANK_MODULUS);
            for (target_col, out_value) in out.iter_mut().enumerate() {
                let term = factor * self.scaled_inverse_mod[source_col][target_col];
                *out_value = (*out_value + term).rem_euclid(RANK_MODULUS);
            }
        }
        out
    }
}

fn projected_nonzero_row(
    coords: &[Integer],
    missing_cols: &[usize],
) -> Option<(usize, Vec<Integer>)> {
    let projected = missing_cols
        .iter()
        .map(|&col| coords[col].clone())
        .collect::<Vec<_>>();
    projected
        .iter()
        .position(|value| *value != 0)
        .map(|pos| (pos, projected))
}

fn projected_rows_are_proportional(
    first_row: &[Integer],
    pivot_col: usize,
    coords: &[Integer],
    missing_cols: &[usize],
) -> bool {
    let first_pivot = &first_row[pivot_col];
    let candidate_pivot = &coords[missing_cols[pivot_col]];
    first_row
        .iter()
        .zip(missing_cols.iter())
        .all(|(first_value, &col)| &coords[col] * first_pivot == first_value * candidate_pivot)
}

fn scaled_inverse_matrix(inverse: &[Vec<Rational>]) -> Option<Vec<Vec<Integer>>> {
    let mut scale = Integer::from(1);
    for row in inverse {
        for value in row {
            let (_, denominator) = value.clone().into_numerator_and_denominator();
            scale = lcm_integer(&scale, &Integer::from(denominator));
        }
    }

    let scale_rational = Rational::from(&scale);
    inverse
        .iter()
        .map(|row| {
            row.iter()
                .map(|value| Integer::try_from(value * scale_rational.clone()).ok())
                .collect()
        })
        .collect()
}

fn active_basis_quotient(active: &[usize], dim: usize) -> (usize, Vec<usize>) {
    let mut basis_active = vec![false; dim];
    let mut basis_count = 0usize;
    for &idx in active {
        if idx < dim && !basis_active[idx] {
            basis_active[idx] = true;
            basis_count += 1;
        }
    }

    let missing_cols = basis_active
        .iter()
        .enumerate()
        .filter_map(|(idx, active)| (!active).then_some(idx))
        .collect();
    (basis_count, missing_cols)
}

pub(super) fn active_rank_at_least(
    active: &[usize],
    hyperplanes: &[DdmHyperplane],
    dim: usize,
    target_rank: usize,
    rank_context: &mut RankContext,
    rank_stats: &mut RankStats,
) -> bool {
    rank_stats.checks += 1;
    if rank_context.has_basis() && basis_prefix_rank_at_least(active, dim, target_rank) {
        rank_stats.basis_hits += 1;
        return true;
    }

    if let Some(value) = rank_context.quotient_rank_at_least(active, hyperplanes, target_rank) {
        rank_stats.quotient_checks += 1;
        if value {
            rank_stats.quotient_true += 1;
            return true;
        }
        rank_stats.quotient_false += 1;
        return false;
    }

    if active_rank_mod_prime_at_least(active, hyperplanes, dim, target_rank) {
        rank_stats.modular_hits += 1;
        return true;
    }

    rank_stats.integer_checks += 1;
    let value = active_rank_integer_at_least(active, hyperplanes, dim, target_rank);
    if value {
        rank_stats.integer_true += 1;
        return value;
    }
    rank_stats.integer_false += 1;
    false
}

fn basis_prefix_rank_at_least(active: &[usize], dim: usize, target_rank: usize) -> bool {
    active
        .iter()
        .filter(|&&idx| idx < dim)
        .take(target_rank)
        .count()
        >= target_rank
}
