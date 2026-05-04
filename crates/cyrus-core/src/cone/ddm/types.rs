use std::cmp::Ordering;

use malachite::Integer;

#[derive(Clone, Debug)]
pub(super) struct DdmRay {
    pub(super) coeffs: Vec<i128>,
    pub(super) active: Vec<usize>,
}

#[derive(Clone, Debug)]
pub(super) struct DdmHyperplane {
    pub(super) dense: Vec<i128>,
    pub(super) sparse: Vec<(usize, i128)>,
}

impl DdmHyperplane {
    pub(super) fn new(dense: Vec<i128>) -> Self {
        let sparse = dense
            .iter()
            .enumerate()
            .filter_map(|(idx, &value)| (value != 0).then_some((idx, value)))
            .collect();
        Self { dense, sparse }
    }

    pub(super) fn dense(&self) -> &[i128] {
        &self.dense
    }

    #[cfg(test)]
    pub(super) fn dot(&self, ray: &[i128]) -> i128 {
        self.dot_i128(ray)
            .expect("DDM hyperplane dot product overflowed i128")
    }

    pub(super) fn dot_sign(&self, ray: &[i128]) -> Ordering {
        self.dot_i128(ray).map_or_else(
            || self.dot_integer(ray).cmp(&Integer::from(0)),
            |dot| dot.cmp(&0),
        )
    }

    pub(super) fn dot_integer(&self, ray: &[i128]) -> Integer {
        self.sparse
            .iter()
            .fold(Integer::from(0), |acc, &(idx, value)| {
                acc + Integer::from(value) * Integer::from(ray[idx])
            })
    }

    fn dot_i128(&self, ray: &[i128]) -> Option<i128> {
        self.sparse.iter().try_fold(0i128, |acc, &(idx, value)| {
            value
                .checked_mul(ray[idx])
                .and_then(|product| acc.checked_add(product))
        })
    }
}

#[derive(Default)]
pub(super) struct StepStats {
    pub(super) positive: usize,
    pub(super) zero: usize,
    pub(super) negative: usize,
    pub(super) pairs: usize,
    pub(super) adjacent: usize,
    pub(super) created: usize,
    pub(super) rank_cache_entries: usize,
    pub(super) redundant: bool,
    pub(super) rank: RankStats,
}

#[derive(Default)]
pub(super) struct RankStats {
    pub(super) checks: usize,
    pub(super) too_short: usize,
    pub(super) cache_hits: usize,
    pub(super) basis_hits: usize,
    pub(super) modular_hits: usize,
    pub(super) quotient_checks: usize,
    pub(super) quotient_true: usize,
    pub(super) quotient_false: usize,
    pub(super) integer_checks: usize,
    pub(super) integer_true: usize,
    pub(super) integer_false: usize,
}

impl StepStats {
    pub(super) fn from_partitions(positive: usize, zero: usize, negative: usize) -> Self {
        Self {
            positive,
            zero,
            negative,
            pairs: positive * negative,
            ..Self::default()
        }
    }
}
