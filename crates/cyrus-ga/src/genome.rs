//! Flux genome and evolutionary operators.
//!
//! A genome is an integer flux pair `(K, M)` in the computed dual divisor
//! basis of a fixed geometry. The operators follow the design validated by
//! the Python prototype (string_theory research notes): two-point and
//! uniform crossover over the concatenated `[K | M]` vector, and adaptive
//! per-lineage mutation whose rate and strength grow with stagnation.

use rand::Rng;
use rand_chacha::ChaCha8Rng;
use serde::{Deserialize, Serialize};

/// An integer flux pair in the dual divisor basis.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Genome {
    /// RR flux vector K (length h21).
    pub k: Vec<i64>,
    /// NSNS flux vector M (length h21).
    pub m: Vec<i64>,
}

impl Genome {
    /// Sample a uniform random genome with entries in `[-flux_range, flux_range]`.
    pub fn random(rng: &mut ChaCha8Rng, h21: usize, flux_range: i64) -> Self {
        let sample = |rng: &mut ChaCha8Rng| -> Vec<i64> {
            (0..h21)
                .map(|_| rng.gen_range(-flux_range..=flux_range))
                .collect()
        };
        Self {
            k: sample(rng),
            m: sample(rng),
        }
    }

    /// Mutate in place: each component independently shifts by up to
    /// `strength` with probability `rate`, clamped to the flux range.
    pub fn mutate(&mut self, rng: &mut ChaCha8Rng, rate: f64, strength: i64, flux_range: i64) {
        let mutate_vec = |rng: &mut ChaCha8Rng, v: &mut Vec<i64>| {
            for entry in v.iter_mut() {
                if rng.gen_bool(rate) {
                    let delta = rng.gen_range(-strength..=strength);
                    *entry = (*entry + delta).clamp(-flux_range, flux_range);
                }
            }
        };
        mutate_vec(rng, &mut self.k);
        mutate_vec(rng, &mut self.m);
    }

    /// Two-point crossover over the concatenated `[K | M]` vector.
    pub fn crossover_two_point(&self, other: &Self, rng: &mut ChaCha8Rng) -> Self {
        let len = self.k.len() + self.m.len();
        let mut a = rng.gen_range(0..len);
        let mut b = rng.gen_range(0..len);
        if a > b {
            std::mem::swap(&mut a, &mut b);
        }
        let concat_self: Vec<i64> = self.k.iter().chain(self.m.iter()).copied().collect();
        let concat_other: Vec<i64> = other.k.iter().chain(other.m.iter()).copied().collect();
        let child: Vec<i64> = (0..len)
            .map(|i| {
                if i >= a && i < b {
                    concat_other[i]
                } else {
                    concat_self[i]
                }
            })
            .collect();
        let (k, m) = child.split_at(self.k.len());
        Self {
            k: k.to_vec(),
            m: m.to_vec(),
        }
    }

    /// Uniform crossover: each component independently from either parent.
    pub fn crossover_uniform(&self, other: &Self, rng: &mut ChaCha8Rng) -> Self {
        let pick = |rng: &mut ChaCha8Rng, a: &[i64], b: &[i64]| -> Vec<i64> {
            a.iter()
                .zip(b.iter())
                .map(|(&x, &y)| if rng.gen_bool(0.5) { x } else { y })
                .collect()
        };
        Self {
            k: pick(rng, &self.k, &other.k),
            m: pick(rng, &self.m, &other.m),
        }
    }
}

/// Adaptive mutation schedule from the validated prototype: rate and
/// strength escalate with the lineage's stagnation (generations since its
/// fitness last improved).
#[must_use]
pub fn adaptive_mutation(rng: &mut ChaCha8Rng, stagnation: u32) -> (f64, i64) {
    let pick = |rng: &mut ChaCha8Rng, options: &[(f64, i64)], weights: &[f64]| {
        let total: f64 = weights.iter().sum();
        let mut roll = rng.gen_range(0.0..total);
        for (option, &weight) in options.iter().zip(weights.iter()) {
            if roll < weight {
                return *option;
            }
            roll -= weight;
        }
        *options.last().expect("options are non-empty")
    };
    match stagnation {
        0..=2 => pick(rng, &[(0.1, 1), (0.2, 2), (0.3, 3)], &[0.6, 0.3, 0.1]),
        3..=6 => pick(rng, &[(0.2, 2), (0.3, 3), (0.5, 5)], &[0.4, 0.4, 0.2]),
        7..=11 => pick(rng, &[(0.3, 3), (0.5, 5), (0.7, 8)], &[0.3, 0.4, 0.3]),
        _ => pick(rng, &[(0.5, 5), (0.7, 8), (1.0, 12)], &[0.3, 0.4, 0.3]),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn random_genome_is_deterministic_for_seed() {
        let mut rng_a = ChaCha8Rng::seed_from_u64(7);
        let mut rng_b = ChaCha8Rng::seed_from_u64(7);
        let a = Genome::random(&mut rng_a, 4, 15);
        let b = Genome::random(&mut rng_b, 4, 15);
        assert_eq!(a, b);
        assert!(
            a.k.iter()
                .chain(a.m.iter())
                .all(|&x| (-15..=15).contains(&x))
        );
    }

    #[test]
    fn mutation_respects_flux_range() {
        let mut rng = ChaCha8Rng::seed_from_u64(11);
        let mut genome = Genome::random(&mut rng, 5, 10);
        for _ in 0..100 {
            genome.mutate(&mut rng, 1.0, 12, 10);
            assert!(
                genome
                    .k
                    .iter()
                    .chain(genome.m.iter())
                    .all(|&x| (-10..=10).contains(&x))
            );
        }
    }

    #[test]
    fn crossover_preserves_lengths_and_mixes_parents() {
        let mut rng = ChaCha8Rng::seed_from_u64(3);
        let a = Genome {
            k: vec![1; 4],
            m: vec![2; 4],
        };
        let b = Genome {
            k: vec![-1; 4],
            m: vec![-2; 4],
        };
        for _ in 0..32 {
            let child = a.crossover_two_point(&b, &mut rng);
            assert_eq!(child.k.len(), 4);
            assert_eq!(child.m.len(), 4);
            for (&x, _) in child.k.iter().zip(a.k.iter()) {
                assert!(x == 1 || x == -1);
            }
            let child_u = a.crossover_uniform(&b, &mut rng);
            assert_eq!(child_u.k.len(), 4);
        }
    }

    #[test]
    fn adaptive_mutation_escalates_with_stagnation() {
        let mut rng = ChaCha8Rng::seed_from_u64(5);
        let calm: Vec<(f64, i64)> = (0..64).map(|_| adaptive_mutation(&mut rng, 0)).collect();
        let desperate: Vec<(f64, i64)> = (0..64).map(|_| adaptive_mutation(&mut rng, 20)).collect();
        let mean_strength =
            |xs: &[(f64, i64)]| xs.iter().map(|&(_, s)| s as f64).sum::<f64>() / xs.len() as f64;
        assert!(mean_strength(&desperate) > mean_strength(&calm));
    }
}
