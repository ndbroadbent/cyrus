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
| CYTools/cygv generic GV wrapper contract | `compute_gv_invariants` now uses CYTools-style `min_points=100*h11` lattice augmentation even when `max_deg` is supplied; bounded `max_deg` lattice enumeration is isolated in `compute_gv_invariants_with_degree_bounded_lattice`; compact GV execution still runs the upstream Rust `cygv` HKTY modules, but Cyrus now constructs the `cygv::Semigroup`, maps semigroup/intersection/fundamental-period/instanton/series-inversion failures into `Result` errors, and catches remaining upstream cygv panics in unwind builds instead of relying on `cygv::compute_gv_rat_threefold` unwraps; CYTools-style matrix-basis projection helpers now match `mori_cap_matrix.dot(basis.T)`; matrix divisor-basis to dual curve-basis construction is available through `compute_curve_basis_matrix_from_divisor_basis_matrix`; `DivisorBasis` and `compute_curve_basis_matrix_for_divisor_basis` make vector-vs-matrix curve-basis dispatch explicit; `divisor_basis_glsm_coordinate_matrix` and `divisor_basis_change_matrix` centralize index/matrix basis-coordinate transforms; `project_mori_cone_cap_rays_for_divisor_basis` and `gv_divisor_basis_data` bundle Mori projection, curve-basis construction, and no-origin q-matrix construction for either basis shape; `intersection_in_divisor_basis` now dispatches vector filtering and matrix-basis dense tensor pullback for the in-basis intersection tensor passed to `cygv`; `matrix_basis_quintic_handoff_runs_actual_cygv_degree_one` verifies that matrix-basis-constructed q data reaches the actual Rust `cygv` call and reproduces the quintic degree-one GV `2875`; `mcallister_first_principles`, `mcallister_gv`, and `mcallister_racetrack` vector-basis GV paths now use that bundled handoff instead of independently projecting Mori rays and hand-building q matrices; `curve_basis_matrix_without_origin_i64` and `curve_basis_q_matrix_for_divisor_basis_i64` centralize the `curve_basis(include_origin=False, as_matrix=True)` q-matrix boundary used by cygv call sites; `mcallister_first_principles --dual-basis` now accepts index or matrix source-coordinate bases for K/M flux transforms into the computed vector basis | Core GV input construction now supports vector or matrix bases for Mori rays, q matrices, curve-basis matrices, basis-coordinate transforms, in-basis intersections, and a small actual-cygv matrix-basis compact regression. McAllister runner now accepts matrix source bases for K/M flux-coordinate validation, but its production primal/dual/KKLT basis remains vector-indexed until Kähler and branch-volume transforms are generalized |
| Mirror-side racetrack GV data | `dual_curves*.dat` are classified as low-dimensional validation checkpoints; generic GV path maps basis curves back to ambient classes | Implemented path exists; full large check remains expensive |
| High-dimensional selected small toric curves | Cyrus computes Mori-cap rays, applies the input-chamber volume cutoff, and matches the 4-214-647 `small_curves.dat` checkpoint with pair pruning | Implemented for checkpoint rule |
| Small-curve pruning semantics | `CurvePruningStrategy::{PairDecomposable, FiniteSemigroup}` exposes the checkpoint rule and stricter finite-set semigroup diagnostic separately | Implemented, with documented mismatch |
| Small toric GV values | Toric two-face/origin-circuit formulas match `small_curves_gv.dat` for 4-214-647; the source-derived sparse resolved-conifold origin-circuit pattern `(-1,-1,1,1)` is now covered in corrected-chamber diagnostics | Implemented for covered selected-toric formulas |
| Potent-ray convergence data | `curve_row_span_rank`, `potent_ray_convergence`, `diagnose_affine_toric_circuit`, CKYZ local-surface identification, source-derived local CKYZ GV extraction, and the Stage 5 potent-ray diagnostics compute rank, corrected-Kähler volumes, affine local-circuit structure, `log xi_n` slopes, first-four GV checks for all 395 rank-two CKYZ rows, all-ten checks for the canonical F0 `[1,1]`/`[1,2]` and F1 `[2,1]`/`[3,1]` families, a first-seven debug regression for polygon-5 `[4,3,2]`, and ignored release all-ten regressions for both polygon-5 directions `[4,3,2]` and `[3,2,2]`; targeted CKYZ extraction now has the broad past-downset API, a generated local semigroup-domain API, a support-predicted API used by the McAllister rank-two CKYZ gate, and coefficient-level z-residual `Li2(q_N)` updates for narrowed targets; `ckyz_local_surface_causal_domain_spec` derives coordinate generators and finite-limit grading from matched CKYZ source data for guardrail checks | Partially implemented; full ten-entry extraction across all families in the default debug all-row gate, rank-four contexts, direct coefficient-history domains for large McAllister supports, and generated low-dimensional-face ray sampling still missing |
| Finite-cutoff nilpotent-ray preclassification | `detect_apparent_nilpotent_ray_from_gv_multiples` ports the McAllister appendix criterion for adding a primitive ray to the apparently nilpotent set from finite GV multiples; `detect_apparent_nilpotent_rays_from_gv_table` scans a finite GV table for primitive seed rays and builds the initial `N` set; `partition_finite_cutoff_gv_charges_by_nilpotence` exposes the resulting finite `C \ N` charge partition; `finite_cutoff_gv_charges_excluding_primitive_rays` builds finite-table complements such as `C \ F0` for the second nop pass; `finite_gv_nonzero_degree_slice_points` extracts exact-degree nonzero charges from a finite GV table for diagnostics; `nilpotent_ray_degree_slice_for_cutoff_fraction` computes the half/full cutoff `k*C` slice origins; `nilpotent_ray_slice_comparison_points` enumerates comparison-ray integer points on those slices; `nilpotent_ray_lll_reduced_slice_distance` applies the CYTools-style LLL transform to an explicit same-degree slice lattice and computes the comparison-ray infinity norm; `nilpotent_ray_divergence_check_from_slice_distances` and `nilpotent_ray_divergence_check_with_explicit_slice_lattices` perform the paper's `d' > d` comparison once half/full slice lattices are supplied; `classify_nilpotent_rays_from_two_pass_divergence_checks` maps first-pass and second-pass checks to provisional `F0` and final `F`; `classify_nop_rays_from_finite_gv_table` applies the finite-cutoff algorithm to an already-computed finite GV table; `check_extremal_mori_ray_separator` verifies exact separator certificates for primitive Mori-ray extremality, and `find_extremal_mori_ray_separator` enumerates exact DDM separators for a supplied finite generator cone | Finite-table `N`/`C \ N`, complements/slice extraction, slice construction, LLL/norm comparison, two-pass classification control flow, and exact finite-generator ray-extremality verification/search exist; the reusable blocker is still upstream: producing/certifying the finite GV table and chamber/slice source without replay, plus proving the supplied generator set is the relevant Mori/chamber context and deriving flop data |
| Flop/corrected-chamber continuation | Negative small-curve volumes and real-axis dilog branch behavior are classified; even-parity branch-cut failures are explicit via `GvDilogFailure`; exact flop transforms for `kappa`, `c2`, GV reassignment, the stable-Weyl reflection matrix, the Weyl-vs-flop tensor equality check, the exact divisor-quadratic vanishing check on a curve hyperplane, `check_stable_weyl_candidate_certificate`, and finder-backed `find_stable_weyl_candidate_certificate` are now exposed as algebraic primitives | Diagnosed, not resolved |
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
   gate through the support-predicted API in about 132 seconds. The newer
   coefficient-level z-residual extraction avoids materializing full `Li2(q_N)`
   series for narrowed targets and makes focused polygon-5 `[4,3,2]`, `N=7`
   pass in about 46 seconds in debug; `N=8` and `N=9` also pass as debug
   diagnostics in about 118 seconds and 274 seconds respectively. Focused
   polygon-5 `[4,3,2]`, `N=10` still exceeds a 300 second debug timeout, but it
   passes in release in about 62 seconds and now has an ignored full-row
   regression. The other polygon-5 direction `[3,2,2]` passes `N=10` in release
   in about 11 seconds and has the same ignored full-row coverage. A traced
   debug `[4,3,2]`, `N=10` run shows support prediction is no longer
   the slow step (`broad=26691`, `selected=21721`, about 9 seconds); z-history
   selection produces `5235` residual degrees in about 17 seconds, and the
   coefficient-level residual extraction still reaches only grading 4 before a
   180 second timeout. The remaining task is therefore making the live
   residual/source-history path practical for broad all-row gates and rank-four
   contexts, not a compact-`cygv` reimplementation.
   Compact hypersurface GV work remains an integration with the actual Rust
   `cygv` crate; local CKYZ code is only for noncompact/source-local checks that
   cannot be represented as a valid compact `cygv` input.
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
   `potent_rays*.dat`. Rank-two CKYZ rows now have first-four checks, the F0
   `[1,1]`/`[1,2]` and F1 `[2,1]`/`[3,1]` families have explicit all-ten
   first-principles regressions, and both polygon-5 source directions have
   ignored release-verified all-ten regressions using coefficient-level residual
   updates. Next
   make the local extractor source/history-targeted enough for complete rows
   across the larger directions.
   A first causal-domain extraction API now exists and is checked on local
   P2/F0 examples, and matched CKYZ sources can now derive their coordinate
   generators and finite-limit grading for guardrail computations. That full
   source-weighted semigroup is still too broad for the default all-row gate, so
   the next step is the narrower source/history domain needed to raise the
   rank-two potent-ray checks beyond the current first-four all-row gate and
   make full-row checks practical beyond targeted release diagnostics. After
   that, handle the rank-four affine supports.
2. Close the remaining CYTools basis contract gap for GA use: matrix-basis
   projection, dual curve-basis construction, no-origin q-matrix construction,
   and in-basis intersection tensor pullback now exist, and the current
   vector-basis runner GV paths use the bundled handoff where applicable.
   `mcallister_first_principles --dual-basis` now transforms K/M from either
   index or matrix source coordinates into the computed vector basis.
   Higher-level APIs still need to accept a generic matrix basis end-to-end for
   Kähler coordinates, KKLT/branch-volume logic, and internal production GV
   execution, or reject matrix production bases loudly.
   The compact 4-214-647 dual-polytope handoff has been compared directly
   against CYTools latest at the cygv boundary: grading vector, no-origin
   q-matrix, and in-basis intersection numbers match exactly, and the unique
   augmented generator set is identical (`496` unique rows; CYTools passes
   `505` rows with duplicates). So the standard compact GV wrapper is not
   blocked by a handoff mismatch.
3. Reconstruct the actual corrected-chamber context for the nine remaining
   origin-circuit misses. This means either an exact supporting face/chamber
   semigroup that can be fed through the cygv/HKTY path, or source-derived
   flop/Weyl-continuation semantics for the relevant branch. Cyrus now has the
   exact algebraic flop transforms, the exact Weyl reflection matrix primitive,
   the exact Weyl-transformed intersection tensor comparison against the
   flop-updated tensor, and the exact check that a candidate divisor's
   quadratic volume form vanishes identically on the hyperplane `C.t=0`. It also
   has an exact Farkas-style separator verifier and exact DDM separator finder
   for the cone-theoretic statement that a target curve spans an extremal ray
   of a supplied finite Mori-generator cone, plus
   `check_stable_weyl_candidate_certificate` and
   `find_stable_weyl_candidate_certificate` to require all of these algebraic
   subchecks together. These are necessary pieces for the stable-Weyl test, but
   they are not a chamber certificate by themselves. Cyrus still needs the
   source-level certification step that proves the supplied finite generators
   are the relevant Kähler wall/chamber context and identifies the shrinking
   divisor from geometry rather than from a fitted candidate list. The
   LP-witness face diagnostic is useful but currently uncertified, so it must
   remain a diagnostic rather than a fallback. The immediate audit checkpoint is to
   reconstruct the smallest exact semigroup/chamber context for these
   higher-rank origin circuits. The prior checkpoint of distinguishing exact
   integer-semigroup decompositions from rational-cone-only LP witnesses is now
   complete: the nine current solved-t misses split into five integer-semigroup
   LP active-generator decompositions and four rational-cone-only
   decompositions. That split is triage only, because cygv assigns GV values
   from a finite semigroup and degree-ordered subtraction history, not from a
   sparse curve relation in isolation. A structured export now exists via
   `--dump-corrected-chamber-gv-context`; it records the local charge basis,
   target degree, degree-bounded projected Mori rays, exact active-generator
   decompositions, q-matrix, curve-basis matrix, grading vector, and sparse
   corrected-chamber intersection entries needed for a finite HKTY check
   without reading `small_curves_gv.dat`. `mcallister_gv_context` now consumes
   that artifact and verifies the five integer-semigroup active-generator
   diamonds, but those tiny domains are not valid fallback contexts: three
   return `GV=0` and two produce non-integral cygv series-inversion output.
   A slightly broader active-support generator window has also been checked:
   filtering degree-bounded Mori rows to the union of target and LP-active
   supports gives `4..13` generators per target, but still returns only zeros
   or `NonIntegerGVError` panics in cygv.
   The context report also measures larger support-overlap/closure windows:
   one-hop overlap already gives `48..227` rows, and four closure layers reach
   `366..2814` rows, so naive support closure returns to the intractable broad
   domain. Executing the small overlap-4 windows (`3..19` generators) also
   fails as a source domain: seven `NonIntegerGVError` panics and two zeros.
   `mcallister_gv_context --run-support-overlap-generators 0` now uses the same
   report fields to try all positive degree-bounded generators up to each
   target degree through cygv's provided-generator path; combine it with
   `--support-overlap-max-target-degree` to keep this probe bounded to the
   low-degree targets. This remains diagnostic because the exported rows are
   degree-bounded projected generators, not a certified full chamber semigroup
   or CYTools lattice-augmented Mori-cap context.
   The actual Rust `cygv` crate has now been exercised on that bounded path.
   The release profile now uses `panic=unwind` so candidate-local cygv panic
   guards can return explicit diagnostic errors instead of aborting a GA/search
   process. The optimized release run with
   `--run-support-overlap-generators 0 --support-overlap-max-target-degree 10`
   entered `cygv` with `720` supplied degree-bounded generators and did not
   finish within `1200s`, and the debug/unwind run did not finish within
   `600s`. Passing `--pair-reduce-support-overlap-generators` mirrors cygv's
   own pair-sum seed pruning before the crate call and reduced the same
   degree-10 input to `450` generators, but that optimized release run
   still exceeded `900s`. These probes produced no GV value or panic. This is a
   measured runtime blocker, not a reason to reimplement compact cygv locally.
   A direct cygv semigroup-size diagnostic now confirms why: the degree-10
   missing targets already have `720` seed rows at or below the target degree,
   and an unguarded cygv semigroup closure measurement did not finish within
   `120s`. The missing object is therefore the source-derived generated
   semigroup/chamber history, not a small visible-generator subset.
   A guarded actual-cygv degree ladder now sharpens that scale estimate without
   reimplementing the compact semigroup: in release mode the degree-10 target
   seed set reaches `9, 291, 2412, 42228, 324773` cygv semigroup elements by
   ladder degrees `1..5`, from effective seed counts
   `8, 254, 293, 380, 397`; degree `6` already exceeded a `180s` release cap.
   The bounded path-history probe reinforces this: with
   `--element-limit 100000`, both degree-10 targets are already present in the
   partial closure, but the closure still exceeds `100000` elements before
   completing, so Cyrus cannot yet recover cygv's predecessor-subtraction
   history from this broad domain.
   Cyrus now exposes cygv's private pair-sum seed-reduction stage via
   `cygv_pair_reduced_seed_generators`, and `mcallister_gv_context` reports raw
   seed counts, reduced seed counts, and whether each target survives that
   reduction. On the current corrected-chamber context, the measured
   low/mid-degree targets all survive as cygv-reduced seeds: degree 10
   `720 -> 450`, degree 12 `905 -> 486`, degree 18 `1616 -> 702`, and degree 22
   `2212 -> 949`. So the missing targets are not artifacts of cygv's first
   decomposable-seed pruning pass.
   The same context report now mirrors cygv's `compute_omega` negative-GLSM-
   intersection buckets. None of the nine missing targets are in the
   `>2` ignored bucket: the six-point degree `18/22/24` targets are `neg2`,
   and the five-point degree `10/12/26` targets are `neg1`. For degree 10 the
   seed histogram is `{neg1:228, neg2:492}` and the reduced-seed histogram is
   `{neg1:55, neg2:395}`. Thus the low-degree misses are not simply discarded
   before cygv's fundamental-period coefficient formulas.
   The report now mirrors cygv's first series-inversion coordinate choice too:
   all nine targets have nonzero corrected-chamber intersection support for
   the coordinate cygv would read, with first-coordinate `kappa_{i,a,b}` pair
   counts `20,20,19,19,14,13,13,45,6`. So the target coefficients are not
   ruled out by a zero intersection row at the selected series coordinate.
   The context report also has a bounded path-history probe that mirrors
   cygv's degree-trimmed seed set, pair-sum seed reduction, additive closure,
   and previous-`q_N` `target - previous` monomial-map lookup. On the current
   context, the two degree-10 missing targets are present in the bounded
   closure, but closure already exceeds a 20000-element guard before exact
   predecessor counts can be certified. This reinforces that the missing input
   is the source-defined finite semigroup/path history, not another small
   active-support or overlap window.
   The context report now preserves origin-circuit provenance from the export:
   the nine misses are affine-rank-three/four origin circuits with local charge
   rows `[2,1,2,-1,-2,-2]`, `[1,2,1,-2,-1,-1]`,
   `[1,1,-2,-1,-1,2]`, `[1,-2,-1,-1,3]`,
   `[1,-2,-1,-1,2,1]`, `[1,-1,-1,1,-3,3]`,
   `[1,-1,-1,-2,3]`, `[1,-2,-1,3,-1]`, and
   `[1,-2,3,-1,-1]`. These are not rank-two CKYZ local-surface supports, and
   the CMS-general divisor checks still fail exact divisor-intersection
   verification. A cygv hypersurface shape diagnostic now records
   `q_rows - q_cols - 1`: the five six-point affine-rank-four rows are
   fourfold-shaped (`cy_dim=4`), while the four five-point affine-rank-three
   rows are threefold-shaped (`cy_dim=3`) but still lack a source-derived
   semigroup, grading, intersection tensor, and chamber interpretation. The
   same report now groups local charge-row permutation signatures: all four
   compact-threefold-shaped misses share `[-2,-1,-1,1,3]` and are explicitly
   marked `shape_only_missing_source_derived_cygv_inputs`, with the missing
   local semigroup generators, grading vector, `q` orientation/phase,
   intersection tensor, and target-coordinate map listed in the JSON. The core
   affine-circuit diagnostic and McAllister context export now preserve
   full-rank local affine coordinates for these rank-three/four supports rather
   than only rank-two coordinates. The context tool also reconstructs them from
   relation points in older dumps; on the saved 4-214-647 context the four
   compact-threefold-shaped misses now expose five rank-three local points
   each. This gives the next source-derived construction actual local support
   points to work from. The context report also now emits a local `cygv` input
   skeleton for these supports. For all four compact-threefold-shaped misses,
   the witness relation is integral in the one-row local charge basis with
   coordinate `[-1]`, and the transposed local charge rows are recorded as the
   candidate local `q` matrix shape. The
   direct source read reinforces this: `cygv` obtains GV values from the full
   finite semigroup, HKTY alpha/beta construction, and degree-ordered
   `Li2(q_N)` subtraction history, while the McAllister paper's small-curve
   discussion is a selected-toric-curve control argument rather than a closed
   formula for every origin-circuit coefficient pattern. The next step is
   therefore to reconstruct the broader source-derived finite semigroup/path
   history for these higher-rank origin circuits, or to certify a
   flop/Weyl-continuation chain that supplies the data from another chamber,
   not to promote LP active-generator diamonds or reuse the rank-two
   local-surface machinery.
4. Implement explicit chamber/flop continuation rules for the Kähler-coordinate
   instanton sums, including the cases where the original real branch is invalid.
   The exact classical-data transform is now available as
   `flop_transform_intersection_numbers` and `flop_transform_c2_vector`, using
   `kappa'_{abc}=kappa_{abc}-n_C^0 C_a C_b C_c` and
   `c'_a=c_a+2 n_C^0 C_a`. The exact table transform
   `flop_reassign_gv_invariants` also applies `n'^0_{-C}=n^0_C` and
   `n'^0_C=0` while rejecting conflicting duplicate data. The remaining work is
   the certification layer: Cyrus still needs to identify the shrinking curve
   and certify its `n_C^0` before these transforms can drive a corrected-chamber
   instanton sum.
