# Cyrus

**High-performance Calabi-Yau manifold computations for string theory research.**

Cyrus is a Rust toolkit for evaluating Calabi-Yau compactifications, computing
moduli stabilization, and searching the string theory landscape for
configurations matching observed physics. It reimplements the algorithms of
[CYTools](https://cy.tools) (polytopes, triangulations, intersection numbers,
cones) plus the KKLT/racetrack stabilization pipeline of
Demirtas–Kim–McAllister–Moritz–Rios-Tascon (arXiv:2107.09064), validated by
reproducing the published results from first principles.

## Status

- Single crate: `cyrus-core` (polytopes, FRSTs, intersection numbers,
  GLSM/bases, cones, GV invariants via [cygv], racetrack, KKLT solve,
  volumes, vacuum energy, quintessence/cosmology utilities).
- **Validation**: the no-replay pipeline reproduces McAllister et al.
  example 4-214-647 end-to-end — polytope points + heights + fluxes +
  orientifold data in, `V₀ ≈ -5.5×10⁻²⁰³ M_pl⁴` out, with every published
  intermediate checkpoint matched at its own precision. See
  `docs/CORRECTED_CHAMBER_RESOLUTION.md` and
  `crates/cyrus-core/tests/mcallister_e2e/`.
- The other four published examples are work in progress.

## Running the validation

```bash
cargo build --release -p cyrus-core --bin mcallister_first_principles

CYRUS_FIRST_PRINCIPLES=1 \
CYRUS_MCALLISTER_DATA_DIR=path/to/paper_data/4-214-647 \
./target/release/mcallister_first_principles
```

Test suites:

```bash
cargo test -p cyrus-core                       # unit + light e2e tests
CYRUS_FIRST_PRINCIPLES=1 \
CYRUS_MCALLISTER_RUNNER_HEAVY=1 \
CYRUS_MCALLISTER_DATA_DIR=... \
cargo test --release -p cyrus-core --test mcallister_e2e  # heavy validation
```

## Documentation

- `CLAUDE.md` / `AGENTS.md` — project philosophy, data policy, type system.
- `docs/` — working notes; start with `docs/README.md`.
- `docs/MCALLISTER_DATA_POLICY.md` — what counts as input vs. checkpoint
  (the "no cheating" contract for the reproduction).

## Key References

- Demirtas, Kim, McAllister, Moritz, Rios-Tascon, "Small Cosmological
  Constants in String Theory" (arXiv:2107.09064)
- Demirtas, Rios-Tascon, McAllister, "CYTools" (arXiv:2211.03823)
- Cicoli et al., "From Inflation to Quintessence" (arXiv:2407.03405)

[cygv]: https://crates.io/crates/cygv
