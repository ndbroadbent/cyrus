# GA-Ready McAllister Pipeline Audit

This is a prompt-to-artifact checklist for the current Cyrus implementation of
the McAllister/DKMM compactification pipeline. It is intentionally not a
completion certificate: it records what is implemented, what is verified, and
what still blocks a GA-ready no-replay reproduction.

## Objective Restatement

Cyrus must compute the reusable geometry and physics pipeline from declared
inputs, with McAllister 4-214-647 used only as a validation target. The reusable
path must not load downstream ancillary outputs such as saved Kähler parameters,
saved volumes, precomputed intersection data, or precomputed GV data to make the
run pass. Any remaining mismatch must be explicit and localizable.

## Artifact Checklist

| Requirement | Current evidence | Status |
| --- | --- | --- |
| Separate declared inputs from checkpoints and replay-only outputs | `docs/MCALLISTER_DATA_POLICY.md`; `stage0_artifact_policy_is_explicit_and_complete` in `crates/cyrus-core/tests/mcallister_e2e/stage0_data_integrity.rs` | Implemented policy |
| First-principles mode rejects fixture replay | `crates/cyrus-core/tests/mcallister_e2e.rs` rejects `CYRUS_ALLOW_FIXTURES` with `CYRUS_FIRST_PRINCIPLES`; runner `enforce_modes` rejects `--allow-fixtures` in first-principles mode | Implemented guard |
| Runner can operate with declared inputs only | `stage0_first_principles_runner_accepts_declared_inputs_only_data_dir` copies only declared `.dat` inputs and checks the current no-replay result; this is heavy and opt-in via `CYRUS_MCALLISTER_RUNNER_HEAVY=1` | Implemented, opt-in verified path |
| Polytope, dual polytope, triangulation, basis, intersections, Mori cap | `mcallister_first_principles` computes these and validates against checkpoints; Stage 2/3/5 tests cover component behavior | Largely implemented |
| CYTools/cygv generic GV wrapper contract | `compute_gv_invariants` now uses CYTools-style `min_points=100*h11` lattice augmentation even when `max_deg` is supplied; bounded `max_deg` lattice enumeration is isolated in `compute_gv_invariants_with_degree_bounded_lattice` | Contract clarified and guarded |
| Mirror-side racetrack GV data | `dual_curves*.dat` are classified as low-dimensional validation checkpoints; generic GV path maps basis curves back to ambient classes | Implemented path exists; full large check remains expensive |
| High-dimensional selected small toric curves | Cyrus computes Mori-cap rays, applies the input-chamber volume cutoff, and matches the 4-214-647 `small_curves.dat` checkpoint with pair pruning | Implemented for checkpoint rule |
| Small-curve pruning semantics | `CurvePruningStrategy::{PairDecomposable, FiniteSemigroup}` exposes the checkpoint rule and stricter finite-set semigroup diagnostic separately | Implemented, with documented mismatch |
| Small toric GV values | Toric two-face/origin-circuit formulas match `small_curves_gv.dat` for 4-214-647 | Implemented for covered selected-toric formulas |
| Potent-ray convergence data | `curve_row_span_rank`, `potent_ray_convergence`, `diagnose_affine_toric_circuit`, CKYZ local-surface identification, source-derived local CKYZ GV extraction, and the Stage 5 potent-ray diagnostics compute rank, corrected-Kähler volumes, affine local-circuit structure, `log xi_n` slopes, and first-two GV checks for all 395 rank-two CKYZ rows | Partially implemented; full ten-entry extraction, rank-four contexts, and generated low-dimensional-face ray sampling still missing |
| Flop/corrected-chamber continuation | Negative small-curve volumes and real-axis dilog branch behavior are classified; even-parity branch-cut failures are explicit via `GvDilogFailure` | Diagnosed, not resolved |
| KKLT corrected Kähler solve | Runner reaches a no-replay corrected Kähler vector and corrected volume without loading `corrected_kahler_param.dat` by default | Implemented but not exact |
| Corrected target-volume / GV correction agreement | Diagnostics localize the residual to corrected-chamber GV target corrections, not classical geometry or file semantics | Open blocker |
| Final corrected volume and V0 | Runner computes current no-replay `V_string` and `V0`; downstream comparisons are post-computation only | Computed, residual remains |
| GA suitability on new candidates | Core paths accept upstream geometry/flux/moduli choices and do not require McAllister output files in first-principles mode | Partially ready; GV/corrected-chamber gaps remain |

## Current Blocking Gaps

1. Corrected-chamber GV target correction still does not reproduce the
   checkpoint-implied vector. The leading residual has been ruled out as a
   simple issue with divisor `chi`, gamma indexing, `q.t` branch aggregation,
   local toric formulas, or checkpoint-file semantics.
2. Potent-ray convergence checks now compute rank, volumes, and decay slopes for
   supplied ray/GV samples. Cyrus also reconstructs the first two saved GV
   entries for all 395 rank-two CKYZ potent-ray rows from CKYZ relation rows,
   local period coefficients, mirror-map substitution, and multiple-cover
   extraction, with `potent_rays_gv.dat` used only as the assertion. Cyrus still
   does not generate the full sampled low-dimensional-face ray set, reproduce
   all ten entries efficiently, or handle the rank-four local charge contexts.
3. Pair-pruned selected curves match McAllister's `small_curves.dat`, while a
   stricter finite-semigroup diagnostic removes five additional curves. This is
   exposed as a policy choice, not hidden.
4. The no-replay final volume remains close but not exact; the residual is
   downstream of the corrected-chamber GV/Kähler-coordinate instanton layer.

## Recent Verification

The latest GV contract split was checked with:

```bash
cargo test -p cyrus-core gv::tests:: -- --nocapture
cargo test -p cyrus-core gv::tests::provided_generator_ray_gv_series -- --nocapture
cargo test -p cyrus-core gv::tests::one_dimensional_ray_gv_series -- --nocapture
cargo test -p cyrus-core affine_toric_circuit -- --nocapture
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 cargo test -p cyrus-core --test potent_ray_affine_circuits -- --nocapture
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 cargo test -p cyrus-core --test potent_ray_affine_circuits first_mcallister_local_p2_potent_ray_gvs_are_reconstructed -- --nocapture
cargo test -p cyrus-core local_p2 -- --nocapture
cargo test -p cyrus-core --test mcallister_e2e stage5_gv_computation_roadmap -- --nocapture
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 cargo test -p cyrus-core --test mcallister_e2e stage5_mcallister_potent_ray_checkpoint_quantities_are_computed -- --nocapture
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 cargo test -p cyrus-core --test mcallister_e2e stage5_mcallister_small_toric_curves_match_checkpoint -- --nocapture
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 cargo test -p cyrus-core --test mcallister_e2e stage5_mcallister_small_toric_curve_gvs_match_checkpoint -- --nocapture
cargo check -p cyrus-core --bin mcallister_first_principles
cargo fmt --check
git diff --check -- crates/cyrus-core/src/gv.rs crates/cyrus-core/src/lib.rs docs/CYGV_AUDIT.md crates/cyrus-core/tests/mcallister_e2e/stage5_gv.rs crates/cyrus-core/tests/mcallister_e2e/snapshots/mcallister_e2e__stage5_gv__gv_computation_roadmap.snap
```

## Next Concrete Work

The next productive implementation work is not more final-vector fitting. It is
to make the remaining GV layer more first-principles:

1. Generate potent-ray samples from low-dimensional faces of
   `M_infinity(X)` and compute the ray `N_{nq}` series rather than reading
   `potent_rays*.dat`. Rank-two CKYZ rows now have first-two checks; next make
   the local extractor coefficient-targeted enough for complete rows where
   feasible, then handle the rank-four affine supports.
2. Compare broader corrected-chamber per-curve cygv/general-GV values against
   toric formula values and missing non-toric contributions.
3. Implement explicit chamber/flop continuation rules for the Kähler-coordinate
   instanton sums, including the cases where the original real branch is invalid.
