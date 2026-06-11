# Cyrus

**High-performance Calabi-Yau computations and landscape search for string
theory research, in Rust.**

Cyrus reimplements the algorithms of [CYTools](https://cy.tools) (polytopes,
triangulations, intersection numbers, cones, GV invariants) together with
the perturbatively-flat-vacuum / racetrack / KKLT stabilization pipeline of
Demirtas-Kim-McAllister-Moritz and McAllister et al. (arXiv:2107.09064),
and builds a genetic-algorithm landscape search on top - targeting vacua
whose dark-energy behavior matches the DESI evolving-dark-energy
("quintessence slope") measurements.

Everything is computed from first principles with a strict no-replay
policy: published data files are either *declared model inputs* (polytope
points, chamber heights, flux integers, orientifold choices) or *validation
checkpoints* that Cyrus must reproduce - never values to feed back into the
computation. See `docs/MCALLISTER_DATA_POLICY.md` for the contract.

## Validation status

### All five published McAllister vacua reproduce end-to-end

From declared inputs alone (points, heights, fluxes, KKLT divisor choices),
Cyrus derives the geometry, mirror Gopakumar-Vafa invariants, racetrack
stabilization, instanton-corrected Kahler moduli, and vacuum energy:

| Example | g_s rel. err | c_tau rel. err | corrected V_string (abs err / checkpoint's own floor) | log10\|V0\| |
|---|---|---|---|---|
| 4-214-647 | 2.3e-8 | 9.1e-7 | 4711.4265 (0.006) | -202.26 |
| 5-81-3213 | 3.9e-7 | 5.9e-7 | 198.1154 (0.001) | — |
| 5-113-4627-main | 1.9e-6 | 8.1e-8 | 944.9236 (0.020 / 0.038) | -143.78 |
| 5-113-4627-alt | 9.7e-7 | 2.4e-7 | 388.3804 (0.071 / 0.325) | -213.48 |
| 7-51-13590 (non-favorable) | 1.1e-6 | 4.8e-4 | 141.3468 (0.073 / 1.089) | -56.49 |

Each V_string agreement sits at the corresponding checkpoint pair's own
internal inconsistency floor (their stored targets vs. their stored Kahler
point), which is the best any re-solve can do. Reproducing these required,
among other things: deriving the orientifold B-field from involution
parity, certifying flopped "corrected chambers" by integer kappa
certificates, Newton on the full instanton-corrected KKLT system, exact
no-decomposition certificates for exotic GV invariants, split-divisor Euler
characteristics and Batyrev h11 for non-favorable polytopes, and a
full-term refinement of the racetrack W0. The notable discrepancies are
the checkpoints' own (e.g. `W_0.dat` files store the two-term racetrack
value while their own `c_tau.dat` is consistent only with the full solve).

### Landscape smoke test: exact CYTools agreement across Kreuzer-Skarke

130+ random KS polytopes at h11 = 2..150, both tools run on identical
chambers: **every polytope matches CYTools exactly** on dual points, FRST
simplices, divisor bases, intersection numbers, and h11 (including all 35+
non-favorable samples) - zero silent divergence, sub-second per polytope
even at h11 = 150. The harness lives in `landscape_smoke` (Rust worker) +
`string_theory/landscape_smoke/run_smoke.py` (CYTools driver). It has
caught and led to fixes for real bugs (degenerate convex-hull seeds,
rank-deficient intersection systems, NaN-unsafe residual gates).

## The GA landscape search (`cyrus-ga`)

A resumable genetic algorithm over integer flux pairs `(K, M)` on a fixed
geometry, scoring candidates on vacuum-energy height, weak coupling, and a
**DESI quintessence slope** component: each candidate's energy scale is
dressed as a thawing-axion potential, integrated through the real
Friedmann + Klein-Gordon equations, CPL-fitted, and scored against the
measured (w0, wa). See `docs/GA_PIPELINE.md` for design, physics
conventions, and roadmap.

```bash
cargo run --release -p cyrus-ga -- \
    --data-dir path/to/paper_data/4-214-647 \
    --run-dir runs/desi-search \
    --population 512 --flux-range 20 \
    --seed-flux-basis 3,4,5,8     # optional: anchor with the published vacuum
```

Runs checkpoint their full state (population, hall of fame, RNG) to the
run directory after every generation: stop with Ctrl-C at any time and
resume bit-exactly by re-running the same command. Evaluation is
~10-100 microseconds per candidate across all cores.

## Running the McAllister validation

```bash
cargo build --release -p cyrus-core --bin mcallister_first_principles

# 4-214-647 (flagship example, full assertions):
CYRUS_FIRST_PRINCIPLES=1 \
CYRUS_MCALLISTER_DATA_DIR=path/to/paper_data/4-214-647 \
./target/release/mcallister_first_principles

# Other examples: pass the example's declared flux coordinate basis
# (recorded in crates/cyrus-core/tests/mcallister_e2e/inputs/flux_basis/),
# e.g. 5-81-3213:
D=path/to/paper_data/5-81-3213
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR=$D \
./target/release/mcallister_first_principles --data-dir $D \
    --skip-mcallister-assertions \
    --dual-basis crates/cyrus-core/tests/mcallister_e2e/inputs/flux_basis/5-81-3213.json

# The 5-113-4627 examples additionally need the bounded-impact deferral for
# two sub-tolerance exotic GV invariants (see docs/MCALLISTER_DATA_POLICY.md):
#   --max-missing-gv-impact 3e-3
```

Test suites:

```bash
cargo test --workspace                          # unit + light e2e (no data needed)
CYRUS_FIRST_PRINCIPLES=1 \
CYRUS_MCALLISTER_RUNNER_HEAVY=1 \
CYRUS_MCALLISTER_DATA_DIR=... \
cargo test --release -p cyrus-core --test mcallister_e2e   # heavy validation
CYRUS_MCALLISTER_DATA_DIR=... \
cargo test --release -p cyrus-ga                # GA + known-vacuum integration
```

## Architecture

Workspace crates:

- **`cyrus-core`** - the physics and geometry library:
  - geometry: `polytope`, `triangulation` (FRSTs, secondary cones, regular
    triangulations with exact integer hulls), `cone`, `lattice`,
    `integer_math`
  - topology/algebra: `intersection` (CYTools algorithm with exact rational
    output), `glsm`, `basis`, `curve_basis`, `divisor` (incl. split-divisor
    chi and Batyrev h11 for non-favorable polytopes), `gv`
    (Gopakumar-Vafa invariants via [cygv], toric-curve formulas, exact
    extremality/no-decomposition certificates)
  - physics: `flat_direction` (Demirtas lemma), `racetrack` (dilogarithm
    superpotential, two-term solve + full-term Newton refinement), `kklt`
    (instanton-corrected Newton path-following, flop continuation), `lvs`,
    `volume` (BBHL), `vacuum`, `quintessence`, `cosmology` (Friedmann +
    Klein-Gordon ODE solver producing w(z))
  - search-facing: `pipeline` (`evaluate_vacuum`: the fast
    Demirtas-McAllister scan chain), `policy` (Strict / ForGA / Abort
    boundary handling)
  - binaries: `mcallister_first_principles` (the validation runner),
    `landscape_smoke`, plus stage-focused debug runners
- **`cyrus-ga`** - the landscape GA: flux genomes and operators, graded
  tiered fitness with the DESI slope component, frontier/hall-of-fame
  population management, persistent resumable runs.

### Type-safe numerics

All physics numbers are phantom-typed (`F64<Pos>`, `F64<NonZero>`,
`I64<Finite>`, ...) with compile-time algebra rules: invalid states
(negative volumes, NaN couplings) are unrepresentable, widening is
automatic through arithmetic, and narrowing is explicit and checked.
Computations that fail must fail loudly - no silent fallbacks, anywhere.
See `CLAUDE.md` for the full philosophy.

## Documentation

- `CLAUDE.md` / `AGENTS.md` - project philosophy, data policy, type system.
- `docs/README.md` - index of working notes.
- `docs/MCALLISTER_DATA_POLICY.md` - declared inputs vs. checkpoints (the
  "no cheating" contract), per-example conventions, missing-GV handling.
- `docs/CORRECTED_CHAMBER_RESOLUTION.md` - the flopped-chamber/B-field
  story behind the corrected-volume reproduction.
- `docs/GA_PIPELINE.md` - GA design, DESI fitness, persistence, roadmap.

## Current frontiers

- **Appropriate-phase GV invariants**: two exotic curves on 5-113-4627
  (GV = 3 each) are provably outside every toric-cap computation in the
  input phase; their effect is below the checkpoints' own precision, but
  the paper's phase-changing method remains to be implemented.
- **Quintessence modelling depth**: compute axion decay constants and
  instanton scales from the Kahler metric (replacing the GA's v1 knobs);
  multi-field Kahler-mixing integration for DESI's phantom-crossing
  best fit.
- **Geometry-level search**: a polytope-selection outer loop above the
  per-geometry flux GA, plus full-KKLT refinement of hall-of-fame
  candidates through the validation runner.

## Key references

- Demirtas, Kim, McAllister, Moritz, Rios-Tascon, "Small Cosmological
  Constants in String Theory" (arXiv:2107.09064)
- Demirtas, Kim, McAllister, Moritz, "Vacua with Small Flux Superpotential"
  (arXiv:1912.10047)
- Demirtas, Rios-Tascon, McAllister, "CYTools" (arXiv:2211.03823)
- DESI Collaboration, DR2 BAO measurements of evolving dark energy
- Cicoli et al., "From Inflation to Quintessence" (arXiv:2407.03405)

[cygv]: https://crates.io/crates/cygv
