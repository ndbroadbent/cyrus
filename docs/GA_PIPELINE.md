# cyrus-ga: GA landscape search with DESI quintessence targeting

`crates/cyrus-ga` searches integer flux pairs `(K, M)` on a fixed Calabi-Yau
geometry for perturbatively flat vacua (Demirtas-Kim-McAllister-Moritz)
whose energy scale and thawing-quintessence behavior match the DESI
evolving-dark-energy measurements.

## Quick start

```bash
# Start (or resume) a search on the validated 4-214-647 geometry:
cargo run --release -p cyrus-ga -- \
    --data-dir /path/to/paper_data/4-214-647 \
    --run-dir runs/desi-search \
    --population 512 --flux-range 20 \
    --seed-flux-basis 3,4,5,8       # optional: inject the published vacuum

# Stop any time (Ctrl-C / kill / laptop sleep). Resume with the SAME
# command: state.json in the run dir checkpoints the full GA state
# (population, hall of fame, RNG) after every generation, so resumes are
# bit-exact. Outputs:
#   <run-dir>/state.json          full checkpoint (atomic replace)
#   <run-dir>/hall_of_fame.json   best-ever candidates, pretty-printed
#   <run-dir>/improvements.jsonl  append-only log of new bests
#   <run-dir>/config.json         the run's parameters
```

Default geometry cost is a one-time ~seconds (mirror GV invariants cache to
disk); per-candidate evaluation is ~10-100 microseconds, so a laptop does
millions of evaluations per hour on all cores (rayon).

## Architecture

- **Genome** (`genome.rs`): `(K, M)` integer vectors in the geometry's
  computed dual divisor basis. Operators from the validated Python
  prototype: two-point/uniform crossover, adaptive per-lineage mutation
  (rate and strength escalate with stagnation).
- **Geometry** (`geometry.rs`): one-time mirror-side preparation through the
  same cyrus-core entry points the McAllister runner drives (dual polytope,
  default-chamber FRST, dual intersection numbers, mirror GV invariants).
  `flux_pair_from_index_basis` imports fluxes written in other bases (e.g.
  the published CYTools-2021 vectors).
- **Fitness** (`fitness.rs`): `cyrus_core::evaluate_vacuum` provides the
  physics (tadpole -> N-matrix -> orthogonality -> racetrack -> W0 -> V0).
  Failures map to graded tier penalties (orthogonality violations are
  graded by |K.p| so the GA has a slope toward the Diophantine condition;
  exact-W0-cancellation dead ends rank below missing racetracks after they
  were observed colonizing the frontier). Valid vacua score on:
  1. **Height**: `log10|V0|` vs the observed dark-energy scale (-121.54).
  2. **DESI slope**: the candidate scale is dressed as a thawing-axion
     potential `V = |V0| (1 - cos(phi/f))` and integrated through the real
     Friedmann + Klein-Gordon equations (`cyrus_core::cosmology`); the
     resulting w(z) is CPL-fitted and scored against DESI
     `w0 = -0.45 +/- 0.21`, `wa = -1.8 +/- 0.6`. Off-scale candidates
     genuinely freeze at w = -1 (no slope reward), physically coupling the
     two components.
  3. **Weak coupling**: small g_s.
- **Population** (`population.rs`): frontier selection, hall of fame,
  tournament parents, asteroid restarts on global stagnation; fully
  serializable including the ChaCha8 RNG.

## Physics conventions and v1 approximations

- The evaluation chain is McAllister's own scan-stage approximation: the
  volume entering V0 is the mirror-volume proxy at `t = p/g_s`, not the
  full KKLT-solved volume (a few-log effect on |V0|, irrelevant for
  ranking; the validated runner remains the post-hoc refinement for
  selected candidates). The seeded 4-214-647 vacuum reproduces
  g_s = 0.00911134 to <1e-3 through this chain (integration test).
- The toric Mori-cone gate is intentionally OFF: McAllister's own published
  flat direction pairs at -0.2 with a toric cap ray (the curve flops away
  on the CY); the physical gates are e^K0 positivity and racetrack
  viability.
- Quintessence modelling knobs (documented, to be replaced by the
  cyrus-cosmology integration): single axion field, decay constant
  `--decay-constant` (default 0.5 M_pl), initial displacement 2f, and the
  identification of the axion potential height with the candidate's |V0|.

## Roadmap

- Compute the axion decay constant and instanton scales from the Kahler
  metric instead of a knob; multi-field (Kahler-mixing) integration for
  the DESI phantom-crossing best fit.
- Quadric-projection mutation operator (for fixed M, K.N^{-1}K = 0 is a
  quadric: project mutations onto it instead of random-walking).
- Geometry-level outer loop (polytope selection), per the GA_v2 research
  notes: the inner (K, M) landscape is near-random, so the learning
  capacity belongs at the geometry level.
- Full-KKLT refinement stage for hall-of-fame candidates via the
  McAllister runner's branch-search hooks.
