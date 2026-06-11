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

## Landscape mode (multi-polytope)

`--polytope-file <jsonl>` searches a pool of geometries instead of one
(fetcher: `string_theory/landscape_smoke/fetch_ga_polytopes.py`). A UCB
bandit allocates inner-GA rounds: first visits sweep the pool in ascending
h21 (densest isotropic flux lattice first) with every third round reserved
for exploiting the best polytope found so far. Per-polytope GA state
persists under `<run-dir>/polytopes/<name>/`; `summary.json` holds the
bandit state and `improvements.jsonl` the global best trajectory. Geometry
preparation runs under a deadline (`--prep-timeout-secs`) and failures
mark the polytope dead, loudly. Each geometry's tadpole bound is computed
automatically ((h11+h21)/2 + 1, the McAllister involution class) and caps
the fitness gate - a global `--q-max` can only tighten it.

## Constructive PFV sampling

The orthogonality condition K.p = 0 is exactly the statement that K is an
isotropic vector of the rational quadratic form N(M)^{-1} - testable
EXACTLY in rational arithmetic (`pfv.rs`). Blind integer search never hits
this measure-zero set (observed: 557k evaluations, zero valid vacua);
instead, every generation re-seeds the worst slice of the population with
exact-isotropic genomes and evolution explores along the constraint
manifold. Result of the first 30-minute landscape run: 114 valid vacua,
best at log10|V0| = -120.6 with CPL wa inside the DESI band.

## Verifying candidates

`cyrus-ga --emit-verification-dir <dir> --polytope-file <pool>
--polytope-name <name> --k a,b,.. --m a,b,..` constructs the automated
orientifold (involution, O7s, rigidity, c_i, tadpole - see
`cyrus_core::orientifold`, validated against all five published examples),
checks the candidate against the true tadpole bound, writes
runner-compatible declared inputs, and prints the full first-principles
verification command. Two known gates remain between an emitted candidate
and a fully verified vacuum: chamber selection (the default chamber's
small curves may exceed GV coverage; fix by flop-walking the secondary
fan) and purity (Pfaffian constancy; needs the dual-fourfold Hodge
formulas of arXiv:1712.04946). Both fail loudly.

## Roadmap

- **Chamber search** (landed): `--chamber-search-steps <n>` on the
  McAllister runner walks bistellar flips from the failing chamber,
  gating each neighbor on phase-1 KKLT reachability (with the neighbor's
  own intersection numbers) and scoring by two-face GV coverage of the
  selected small curves; the winner's certified heights are emitted to
  `heights_chamber_search.dat` for adoption + re-run. Open: deeper walks
  for curve-dense geometries (h21_4_135 descends 166 -> 138 uncovered in
  14 flips; full coverage needs longer budgets or extended GV formulas).
- **Purity**: port the combinatorial fourfold Hodge formulas
  (arXiv:1712.04946) so h^{2,1}(D-hat) = 0 is checked per divisor; O7
  stacks are already rigorously pure (pointwise fixed). Until then,
  verified candidates are labelled "modulo Pfaffian constancy".
- Verification-friendliness triage: score hall-of-fame candidates by
  uncovered-curve count before attempting the full solve, and feed it
  back into the bandit so the search prefers verifiable geometries.
- Compute the axion decay constant and instanton scales from the Kahler
  metric instead of a knob; multi-field (Kahler-mixing) integration for
  the DESI phantom-crossing best fit.
- Performance: cache N(M)^{-1} per M-lineage in the isotropic sampler;
  exhaustive exact K-enumeration per M at h21 = 4 (a 31^4 box of integer
  tests costs milliseconds).
