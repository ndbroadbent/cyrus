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
| CYTools/cygv generic GV wrapper contract | `compute_gv_invariants` now uses CYTools-style `min_points=100*h11` lattice augmentation even when `max_deg` is supplied; bounded `max_deg` lattice enumeration is isolated in `compute_gv_invariants_with_degree_bounded_lattice`; CYTools-style matrix-basis projection helpers now match `mori_cap_matrix.dot(basis.T)`; matrix divisor-basis to dual curve-basis construction is available through `compute_curve_basis_matrix_from_divisor_basis_matrix`; `DivisorBasis` and `compute_curve_basis_matrix_for_divisor_basis` make vector-vs-matrix curve-basis dispatch explicit; `project_mori_cone_cap_rays_for_divisor_basis` and `gv_divisor_basis_data` bundle Mori projection, curve-basis construction, and no-origin q-matrix construction for either basis shape; `mcallister_first_principles`, `mcallister_gv`, and `mcallister_racetrack` vector-basis GV paths now use that bundled handoff instead of independently projecting Mori rays and hand-building q matrices; `curve_basis_matrix_without_origin_i64` and `curve_basis_q_matrix_for_divisor_basis_i64` centralize the `curve_basis(include_origin=False, as_matrix=True)` q-matrix boundary used by cygv call sites; `mcallister_first_principles --dual-basis` now parses matrix-basis JSON but rejects it with an explicit unsupported-vector-path error instead of silently treating it as indices | Core GV basis handoff now supports vector or matrix bases, and the runner/useful binaries use it for vector GV paths; McAllister runner remains vector-basis only and must keep rejecting matrix overrides until flux/intersection transforms are generalized |
| Mirror-side racetrack GV data | `dual_curves*.dat` are classified as low-dimensional validation checkpoints; generic GV path maps basis curves back to ambient classes | Implemented path exists; full large check remains expensive |
| High-dimensional selected small toric curves | Cyrus computes Mori-cap rays, applies the input-chamber volume cutoff, and matches the 4-214-647 `small_curves.dat` checkpoint with pair pruning | Implemented for checkpoint rule |
| Small-curve pruning semantics | `CurvePruningStrategy::{PairDecomposable, FiniteSemigroup}` exposes the checkpoint rule and stricter finite-set semigroup diagnostic separately | Implemented, with documented mismatch |
| Small toric GV values | Toric two-face/origin-circuit formulas match `small_curves_gv.dat` for 4-214-647; the source-derived sparse resolved-conifold origin-circuit pattern `(-1,-1,1,1)` is now covered in corrected-chamber diagnostics | Implemented for covered selected-toric formulas |
| Potent-ray convergence data | `curve_row_span_rank`, `potent_ray_convergence`, `diagnose_affine_toric_circuit`, CKYZ local-surface identification, source-derived local CKYZ GV extraction, and the Stage 5 potent-ray diagnostics compute rank, corrected-Kähler volumes, affine local-circuit structure, `log xi_n` slopes, first-four GV checks for all 395 rank-two CKYZ rows, and all-ten checks for the canonical F0 `[1,1]`/`[1,2]` and F1 `[2,1]`/`[3,1]` families; targeted CKYZ extraction now has the broad past-downset API, a cygv-style generated semigroup API, and a support-predicted API that is used by the McAllister rank-two CKYZ gate; `ckyz_local_surface_causal_domain_spec` derives coordinate generators and finite-limit grading from matched CKYZ source data for guardrail checks | Partially implemented; full ten-entry extraction across all families, rank-four contexts, direct coefficient-history domains for large McAllister supports, and generated low-dimensional-face ray sampling still missing |
| Flop/corrected-chamber continuation | Negative small-curve volumes and real-axis dilog branch behavior are classified; even-parity branch-cut failures are explicit via `GvDilogFailure` | Diagnosed, not resolved |
| KKLT corrected Kähler solve | Runner reaches a no-replay corrected Kähler vector and corrected volume without loading `corrected_kahler_param.dat` by default | Implemented but not exact |
| Corrected target-volume / GV correction agreement | Diagnostics localize the residual to corrected-chamber GV target corrections, not classical geometry or file semantics | Open blocker |
| Final corrected volume and V0 | Runner computes current no-replay `V_string` and `V0`; downstream comparisons are post-computation only | Computed, residual remains |
| GA suitability on new candidates | Core paths accept upstream geometry/flux/moduli choices and do not require McAllister output files in first-principles mode | Partially ready; GV/corrected-chamber gaps and full matrix-basis/general-basis handling remain |

## Current Blocking Gaps

1. Corrected-chamber GV target correction still does not reproduce the
   checkpoint-implied vector. The leading residual has been ruled out as a
   simple issue with divisor `chi`, gamma indexing, `q.t` branch aggregation,
   local toric formulas, or checkpoint-file semantics.
   A source-derived sparse resolved-conifold origin-circuit detector now covers
   one previously missing standard `(-1,-1,1,1)` circuit. The latest short
   corrected-chamber diagnostic reports checkpoint-t `toric_covered=412`,
   `toric_missing=8`, and current solved-t `toric_covered=410`,
   `toric_missing=9`; the remaining current missing targets are still
   nonstandard origin-circuit patterns. The latest affine-support diagnostic
   shows those nine remaining misses have ranks `{3: 4, 4: 5}`, so they are
   not rank-two local toric surface supports and cannot be handled by simply
   reusing the CKYZ rank-two potent-ray extractor. The branch-status diagnostic
   now shows the same nine misses are all real-axis evaluable, with `q.t >= 0.1`
   and parity buckets seven even / two odd, so this specific residual is not a
   near-wall `Li2/Li3` branch-cut failure. Earlier LP-witness face diagnostics
   reduced the target-correction delta but exact supporting-face certification
   remained zero for those LP contexts, so those values are still diagnostic
   evidence, not a reusable GV fallback. A fresh read of the
   McAllister GV section reinforces that boundary: the paper's toric-curve
   method is a selection-and-pruning strategy for important nilpotent curves,
   while `cygv` computes through a finite semigroup and degree-ordered
   subtraction history. The next implementation target must reconstruct the
   missing origin-circuit face/chamber semigroup, not infer a GV value directly
   from the sparse coefficient pattern. A fresh release diagnostic confirms
   the nine solved-t misses are still all higher-rank origin-circuit Mori
   generators (`ambient_rays=561596`, `subcutoff=561`, `pair_pruned=419`,
   `toric_covered=410`, `toric_missing=9`), with target grading degrees
   `10..26`. The CMS-general-divisor candidate checks currently fail exact
   divisor-intersection verification for every missing target, so those
   candidate `1/-2/3` values cannot be promoted.
2. Potent-ray convergence checks now compute rank, volumes, and decay slopes for
   supplied ray/GV samples. Cyrus also reconstructs the first four saved GV
   entries for all 395 rank-two CKYZ potent-ray rows, plus all ten entries for
   the canonical F0 source directions `[1,1]`/`[1,2]` and F1 source directions
   `[2,1]`/`[3,1]`, from CKYZ relation rows, local period coefficients,
   mirror-map substitution, and multiple-cover extraction, with
   `potent_rays_gv.dat` used only as the assertion. Cyrus still
   does not generate the full sampled low-dimensional-face ray set, reproduce
   all ten entries efficiently, or handle the rank-four local charge contexts.
   The source-level reason is now explicit: cygv's series inversion subtracts
   lower-degree `Li2(q_N)` history from a semigroup, so a coefficient-targeted
   local domain is needed before the saved ten-entry rows are a fair comparison.
   The current CKYZ support-predicted domain is a prerequisite, not the full
   history domain for ray-direction rows. The gated test can now be raised with
   `CYRUS_CKYZ_MULTIPLES_TO_CHECK` and narrowed with
   `CYRUS_CKYZ_TARGET_DIRECTION`; `N=4` now passes as the default first-principles
   gate through the support-predicted API in about 132 seconds. The older
   target-downset path completed the same gate in about 110 seconds, so the new
   API is a correctness/structure promotion rather than a performance win.
   A focused polygon-5 `[4,3,2]`, `N=5` diagnostic still exceeded a 300 second
   timeout, and an earlier unfiltered `N=5` run exceeded a 600 second timeout,
   so the next step is direct coefficient/path-history support generation rather
   than pushing the existing all-row gate higher. A domain-only diagnostic for
   the same polygon-5 `[4,3,2]`, `N=5` target also exceeded 300 seconds, which
   places the bottleneck inside support-domain prediction before final GV
   extraction.
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
cargo build --release -p cyrus-core --bin mcallister_first_principles
CARGO_PROFILE_RELEASE_PANIC=unwind cargo build --release -p cyrus-core --bin mcallister_first_principles
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 CYRUS_CORRECTED_CHAMBER_LP_FACE_CERTIFICATE=1 timeout 600 ./target/release/mcallister_first_principles --stop-after volume --allow-invalid-ek0 --diagnose-corrected-chamber-gv --diagnose-corrected-chamber-lp-face-gv
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 timeout 300 ./target/release/mcallister_first_principles --stop-after volume --allow-invalid-ek0 --diagnose-corrected-chamber-gv
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_RUNNER_HEAVY=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 cargo test -p cyrus-core --test mcallister_e2e stage0_first_principles_runner_accepts_declared_inputs_only_data_dir -- --nocapture
cargo check -p cyrus-core --bin mcallister_first_principles
cargo fmt --check
git diff --check -- crates/cyrus-core/src/gv.rs crates/cyrus-core/src/lib.rs docs/CYGV_AUDIT.md crates/cyrus-core/tests/mcallister_e2e/stage5_gv.rs crates/cyrus-core/tests/mcallister_e2e/snapshots/mcallister_e2e__stage5_gv__gv_computation_roadmap.snap
```

## Next Concrete Work

The next productive implementation work is not more final-vector fitting. It is
to make the remaining GV layer more first-principles:

1. Generate potent-ray samples from low-dimensional faces of
   `M_infinity(X)` and compute the ray `N_{nq}` series rather than reading
   `potent_rays*.dat`. Rank-two CKYZ rows now have first-four checks, and the
   F0 `[1,1]`/`[1,2]` and F1 `[2,1]`/`[3,1]` families have explicit all-ten
   first-principles regressions; next make the local extractor
   coefficient-targeted enough for complete rows across the larger directions.
   A first causal-domain extraction API now exists and is checked on local
   P2/F0 examples, and matched CKYZ sources can now derive their coordinate
   generators and finite-limit grading for guardrail computations. That full
   source-weighted semigroup is still too broad for the default all-row gate, so
   the next step is the narrower coefficient/path-history domain needed to raise
   the rank-two potent-ray checks beyond the current first-four gate. After
   that, handle the rank-four affine supports.
2. Close the remaining CYTools basis contract gap for GA use: matrix-basis
   projection, dual curve-basis construction, and no-origin q-matrix
   construction now exist, and the current vector-basis runner GV paths use the
   bundled handoff. Higher-level APIs still need to either accept a generic
   matrix basis end-to-end or reject it loudly until that path is ported.
3. Reconstruct the actual corrected-chamber context for the nine remaining
   origin-circuit misses. This means either an exact supporting face/chamber
   semigroup that can be fed through the cygv/HKTY path, or source-derived
   flop-continuation semantics for the relevant `B_2` branch. The LP-witness
   face diagnostic is useful but currently uncertified, so it must remain a
   diagnostic rather than a fallback. The immediate audit checkpoint is to
   reconstruct the smallest exact semigroup/chamber context for these
   higher-rank origin circuits. The prior checkpoint of distinguishing exact
   integer-semigroup decompositions from rational-cone-only LP witnesses is now
   complete: the nine current solved-t misses split into five integer-semigroup
   LP active-generator decompositions and four rational-cone-only
   decompositions. That split is triage only, because cygv assigns GV values
   from a finite semigroup and degree-ordered subtraction history, not from a
   sparse curve relation in isolation. The immediate concrete artifact should
   be a structured export of the missing origin-circuit semigroup/chamber
   context: local charge basis, target degree, lower-degree generators, exact
   active-generator decompositions, and the CYTools/cygv inputs needed to run a
   finite HKTY check without reading `small_curves_gv.dat`.
4. Implement explicit chamber/flop continuation rules for the Kähler-coordinate
   instanton sums, including the cases where the original real branch is invalid.
