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

    pub(super) fn dot(&self, ray: &[i128]) -> i128 {
        self.sparse
            .iter()
            .map(|&(idx, value)| value * ray[idx])
            .sum()
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
