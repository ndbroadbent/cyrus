//! Population management for the flux GA.
//!
//! Frontier selection, hall of fame, adaptive per-lineage mutation, and
//! asteroid restarts - the design the Python prototype validated. The whole
//! state (including the RNG) serializes to JSON so a run can be stopped
//! with Ctrl-C and resumed exactly.

use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde::{Deserialize, Serialize};

use crate::fitness::FitnessReport;
use crate::genome::{Genome, adaptive_mutation};

/// Non-finite-safe float serialization.
///
/// JSON has no -Infinity: serde_json writes non-finite floats as `null`,
/// which cannot deserialize back into `f64`. Fields that legitimately hold
/// NEG_INFINITY sentinels (never-evaluated lineages/polytopes) round-trip
/// through `Option<f64>` instead.
pub mod serde_finite_or_neg_inf {
    use serde::{Deserialize, Deserializer, Serializer};

    /// Serialize non-finite values as null.
    ///
    /// # Errors
    /// Propagates serializer errors.
    pub fn serialize<S: Serializer>(value: &f64, serializer: S) -> Result<S::Ok, S::Error> {
        if value.is_finite() {
            serializer.serialize_some(value)
        } else {
            serializer.serialize_none()
        }
    }

    /// Deserialize null back to NEG_INFINITY.
    ///
    /// # Errors
    /// Propagates deserializer errors.
    pub fn deserialize<'de, D: Deserializer<'de>>(deserializer: D) -> Result<f64, D::Error> {
        Ok(Option::<f64>::deserialize(deserializer)?.unwrap_or(f64::NEG_INFINITY))
    }
}

/// One member of the population.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Individual {
    /// The flux genome.
    pub genome: Genome,
    /// Last fitness report (None until evaluated).
    pub report: Option<FitnessReport>,
    /// Generations since this lineage last improved.
    pub stagnation: u32,
    /// Best fitness this lineage has seen.
    #[serde(with = "serde_finite_or_neg_inf")]
    pub lineage_best: f64,
}

/// GA hyperparameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GaParams {
    /// Total population per generation.
    pub population_size: usize,
    /// Number of elite/frontier individuals carried over.
    pub frontier_size: usize,
    /// Hall-of-fame capacity (best ever, survives restarts).
    pub hall_of_fame_size: usize,
    /// Tournament size for parent selection.
    pub tournament_size: usize,
    /// Probability of crossover (vs pure mutation).
    pub crossover_rate: f64,
    /// Flux entries are sampled/clamped to `[-flux_range, flux_range]`.
    pub flux_range: i64,
    /// Generations without global improvement before an asteroid restart.
    pub asteroid_threshold: u32,
}

impl Default for GaParams {
    fn default() -> Self {
        Self {
            population_size: 512,
            frontier_size: 48,
            hall_of_fame_size: 20,
            tournament_size: 5,
            crossover_rate: 0.5,
            flux_range: 15,
            asteroid_threshold: 25,
        }
    }
}

/// Full serializable GA state: checkpointing this every generation makes
/// stop/resume exact (the RNG state is part of it).
#[derive(Serialize, Deserialize)]
pub struct GaState {
    /// Hyperparameters the run was started with.
    pub params: GaParams,
    /// Current generation index.
    pub generation: u64,
    /// Total evaluations performed.
    pub evaluations: u64,
    /// Current population.
    pub population: Vec<Individual>,
    /// Best-ever individuals (deduplicated by genome).
    pub hall_of_fame: Vec<Individual>,
    /// Best fitness ever seen.
    #[serde(with = "serde_finite_or_neg_inf")]
    pub best_fitness: f64,
    /// Generations since the global best improved.
    pub global_stagnation: u32,
    /// Seeded RNG (serialized so resume is bit-reproducible).
    pub rng: ChaCha8Rng,
}

impl GaState {
    /// Fresh state with a random initial population.
    #[must_use]
    pub fn new(params: GaParams, h21: usize, seed: u64) -> Self {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let population = (0..params.population_size)
            .map(|_| Individual {
                genome: Genome::random(&mut rng, h21, params.flux_range),
                report: None,
                stagnation: 0,
                lineage_best: f64::NEG_INFINITY,
            })
            .collect();
        Self {
            params,
            generation: 0,
            evaluations: 0,
            population,
            hall_of_fame: Vec::new(),
            best_fitness: f64::NEG_INFINITY,
            global_stagnation: 0,
            rng,
        }
    }

    fn fitness_of(individual: &Individual) -> f64 {
        individual
            .report
            .as_ref()
            .map_or(f64::NEG_INFINITY, |r| r.fitness)
    }

    /// Record evaluated reports back into the population and update the
    /// hall of fame, lineage stagnation, and global stagnation counters.
    pub fn absorb_reports(&mut self, reports: Vec<FitnessReport>) {
        assert_eq!(reports.len(), self.population.len());
        self.evaluations += reports.len() as u64;
        let mut improved = false;
        for (individual, report) in self.population.iter_mut().zip(reports) {
            let fitness = report.fitness;
            if fitness > individual.lineage_best {
                individual.lineage_best = fitness;
                individual.stagnation = 0;
            } else {
                individual.stagnation += 1;
            }
            if fitness > self.best_fitness {
                self.best_fitness = fitness;
                improved = true;
            }
            individual.report = Some(report);
        }
        if improved {
            self.global_stagnation = 0;
        } else {
            self.global_stagnation += 1;
        }
        self.update_hall_of_fame();
    }

    fn update_hall_of_fame(&mut self) {
        let mut pool: Vec<Individual> = self.hall_of_fame.clone();
        pool.extend(self.population.iter().cloned());
        pool.sort_by(|a, b| Self::fitness_of(b).total_cmp(&Self::fitness_of(a)));
        pool.dedup_by(|a, b| a.genome == b.genome);
        pool.truncate(self.params.hall_of_fame_size);
        self.hall_of_fame = pool;
    }

    fn tournament_pick(&mut self) -> Individual {
        let mut best: Option<usize> = None;
        for _ in 0..self.params.tournament_size {
            let idx = self.rng.gen_range(0..self.population.len());
            if best.is_none_or(|b| {
                Self::fitness_of(&self.population[idx]) > Self::fitness_of(&self.population[b])
            }) {
                best = Some(idx);
            }
        }
        self.population[best.expect("tournament ran")].clone()
    }

    /// Replace the worst `genomes.len()` individuals with externally
    /// constructed genomes (e.g. exact PFV isotropic samples), resetting
    /// their lineage so they are evaluated fresh next generation.
    pub fn inject_genomes(&mut self, genomes: Vec<Genome>) {
        if genomes.is_empty() {
            return;
        }
        let mut order: Vec<usize> = (0..self.population.len()).collect();
        order.sort_by(|&a, &b| {
            Self::fitness_of(&self.population[a]).total_cmp(&Self::fitness_of(&self.population[b]))
        });
        for (slot, genome) in order.into_iter().zip(genomes) {
            self.population[slot] = Individual {
                genome,
                report: None,
                stagnation: 0,
                lineage_best: f64::NEG_INFINITY,
            };
        }
    }

    /// Produce the next generation: frontier carry-over, tournament +
    /// crossover + adaptive mutation offspring, or an asteroid restart when
    /// globally stagnant.
    pub fn next_generation(&mut self) {
        self.generation += 1;
        if self.global_stagnation >= self.params.asteroid_threshold && !self.hall_of_fame.is_empty()
        {
            self.asteroid_restart();
            return;
        }

        let mut ranked = self.population.clone();
        ranked.sort_by(|a, b| Self::fitness_of(b).total_cmp(&Self::fitness_of(a)));
        ranked.dedup_by(|a, b| a.genome == b.genome);
        let frontier: Vec<Individual> =
            ranked.into_iter().take(self.params.frontier_size).collect();

        let mut next: Vec<Individual> = frontier;
        while next.len() < self.params.population_size {
            let parent = self.tournament_pick();
            let crossover = self.rng.gen_bool(self.params.crossover_rate);
            let mut genome = if crossover {
                let other = self.tournament_pick();
                if self.rng.gen_bool(0.5) {
                    parent
                        .genome
                        .crossover_two_point(&other.genome, &mut self.rng)
                } else {
                    parent
                        .genome
                        .crossover_uniform(&other.genome, &mut self.rng)
                }
            } else {
                parent.genome.clone()
            };
            let (rate, strength) = adaptive_mutation(&mut self.rng, parent.stagnation);
            genome.mutate(&mut self.rng, rate, strength, self.params.flux_range);
            next.push(Individual {
                genome,
                report: None,
                stagnation: parent.stagnation,
                lineage_best: parent.lineage_best,
            });
        }
        self.population = next;
    }

    /// Asteroid impact: rebuild from heavy mutations of the hall of fame
    /// (70%) plus fresh random genomes (30%).
    fn asteroid_restart(&mut self) {
        let h21 = self.hall_of_fame[0].genome.k.len();
        let n_fresh = self.params.population_size * 3 / 10;
        let mut next = Vec::with_capacity(self.params.population_size);
        for _ in 0..n_fresh {
            next.push(Individual {
                genome: Genome::random(&mut self.rng, h21, self.params.flux_range),
                report: None,
                stagnation: 0,
                lineage_best: f64::NEG_INFINITY,
            });
        }
        while next.len() < self.params.population_size {
            let idx = self.rng.gen_range(0..self.hall_of_fame.len());
            let mut genome = self.hall_of_fame[idx].genome.clone();
            genome.mutate(&mut self.rng, 0.7, 8, self.params.flux_range);
            next.push(Individual {
                genome,
                report: None,
                stagnation: 0,
                lineage_best: f64::NEG_INFINITY,
            });
        }
        self.population = next;
        self.global_stagnation = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_report(fitness: f64) -> FitnessReport {
        FitnessReport {
            fitness,
            tier: "test".into(),
            g_s: None,
            w0: None,
            log10_abs_v0: None,
            q_flux: 0.0,
            cpl_w0: None,
            cpl_wa: None,
            deep_tier: None,
            deep_log10_v0: None,
        }
    }

    #[test]
    fn state_roundtrips_through_json_exactly() {
        let mut state = GaState::new(
            GaParams {
                population_size: 8,
                frontier_size: 2,
                ..GaParams::default()
            },
            4,
            42,
        );
        let reports = (0..8).map(|i| dummy_report(f64::from(i))).collect();
        state.absorb_reports(reports);
        state.next_generation();

        let json = serde_json::to_string(&state).expect("serialize");
        let mut resumed: GaState = serde_json::from_str(&json).expect("deserialize");

        // Continuing both must produce identical populations (RNG included).
        state.next_generation();
        resumed.next_generation();
        let genomes: Vec<&Genome> = state.population.iter().map(|i| &i.genome).collect();
        let genomes_resumed: Vec<&Genome> = resumed.population.iter().map(|i| &i.genome).collect();
        assert_eq!(genomes, genomes_resumed);
    }

    #[test]
    fn unevaluated_state_with_neg_infinity_roundtrips() {
        // JSON null <-> NEG_INFINITY: a fresh (never-absorbed) state and
        // post-asteroid lineages hold -inf sentinels and must survive a
        // checkpoint/resume cycle (regression: resumed landscape runs
        // panicked on "parse summary").
        let state = GaState::new(GaParams::default(), 3, 1);
        let json = serde_json::to_string(&state).expect("serialize");
        let resumed: GaState = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(resumed.best_fitness, f64::NEG_INFINITY);
        assert_eq!(resumed.population[0].lineage_best, f64::NEG_INFINITY);
    }

    #[test]
    fn hall_of_fame_keeps_best_and_dedupes() {
        let mut state = GaState::new(
            GaParams {
                population_size: 6,
                hall_of_fame_size: 3,
                ..GaParams::default()
            },
            3,
            7,
        );
        // Force duplicate genomes with different fitness
        let clone = state.population[0].genome.clone();
        state.population[1].genome = clone;
        let reports = vec![
            dummy_report(10.0),
            dummy_report(10.0),
            dummy_report(5.0),
            dummy_report(1.0),
            dummy_report(-3.0),
            dummy_report(7.0),
        ];
        state.absorb_reports(reports);
        assert_eq!(state.hall_of_fame.len(), 3);
        assert!((state.best_fitness - 10.0).abs() < 1e-12);
        assert!(
            state.hall_of_fame[0].genome != state.hall_of_fame[1].genome,
            "hall of fame must deduplicate genomes"
        );
    }

    #[test]
    fn asteroid_restart_triggers_on_stagnation() {
        let mut state = GaState::new(
            GaParams {
                population_size: 8,
                asteroid_threshold: 2,
                ..GaParams::default()
            },
            4,
            9,
        );
        let reports: Vec<FitnessReport> = (0..8).map(|_| dummy_report(1.0)).collect();
        state.absorb_reports(reports);
        state.next_generation();
        for _ in 0..3 {
            let reports: Vec<FitnessReport> = (0..8).map(|_| dummy_report(0.0)).collect();
            state.absorb_reports(reports);
            state.next_generation();
        }
        // Without the asteroid reset the counter would have reached 3.
        assert!(
            state.global_stagnation < 3,
            "asteroid restart should have reset the stagnation counter, got {}",
            state.global_stagnation
        );
    }
}
