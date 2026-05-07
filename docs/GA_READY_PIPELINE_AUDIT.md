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
| CYTools/cygv generic GV wrapper contract | `compute_gv_invariants` now uses CYTools-style `min_points=100*h11` lattice augmentation even when `max_deg` is supplied; bounded `max_deg` lattice enumeration is isolated in `compute_gv_invariants_with_degree_bounded_lattice`; compact GV execution still runs the upstream Rust `cygv` HKTY modules, but Cyrus now constructs the `cygv::Semigroup`, maps semigroup/intersection/fundamental-period/instanton/series-inversion failures into `Result` errors, and catches remaining upstream cygv panics in unwind builds instead of relying on `cygv::compute_gv_rat_threefold` unwraps; `cyrus_direct_cygv_chain_matches_upstream_quintic_wrapper` verifies that this direct module call chain matches `cygv::compute_gv_rat_threefold` on the quintic degree-one `2875` case; CYTools-style matrix-basis projection helpers now match `mori_cap_matrix.dot(basis.T)`; matrix divisor-basis to dual curve-basis construction is available through `compute_curve_basis_matrix_from_divisor_basis_matrix`; `DivisorBasis` and `compute_curve_basis_matrix_for_divisor_basis` make vector-vs-matrix curve-basis dispatch explicit; `divisor_basis_glsm_coordinate_matrix` and `divisor_basis_change_matrix` centralize index/matrix basis-coordinate transforms; `project_mori_cone_cap_rays_for_divisor_basis` and `gv_divisor_basis_data` bundle Mori projection, curve-basis construction, and no-origin q-matrix construction for either basis shape; `intersection_in_divisor_basis` now dispatches vector filtering and matrix-basis dense tensor pullback for the in-basis intersection tensor passed to `cygv`; `matrix_basis_quintic_handoff_runs_actual_cygv_degree_one` verifies that matrix-basis-constructed q data reaches the actual Rust `cygv` call and reproduces the quintic degree-one GV `2875`; `mcallister_first_principles`, `mcallister_gv`, and `mcallister_racetrack` vector-basis GV paths now use that bundled handoff instead of independently projecting Mori rays and hand-building q matrices; `curve_basis_matrix_without_origin_i64` and `curve_basis_q_matrix_for_divisor_basis_i64` centralize the `curve_basis(include_origin=False, as_matrix=True)` q-matrix boundary used by cygv call sites; `mcallister_first_principles --dual-basis` now accepts index or matrix source-coordinate bases for K/M flux transforms, and `--production-dual-basis` carries an index or matrix internal dual divisor basis through flat-direction intersections and compact GV input construction; core KKLT now exposes `compute_kklt_divisor_volumes_in_divisor_basis`, `compute_kklt_jacobian_in_divisor_basis`, `solve_divisor_basis_path_following`, `solve_divisor_basis_path_following_branch_candidates`, `generate_scaled_divisor_basis_branch_initializations`, and `scale_divisor_basis_kklt_branch_initialization_to_target` for matrix Kähler-coordinate bases; `mcallister_first_principles --production-primal-basis` now carries an index or matrix internal primal divisor basis through KKLT initialization scaling, path following, branch search, and branch GV-coverage ranking | Core GV input construction now supports vector or matrix bases for Mori rays, q matrices, curve-basis matrices, basis-coordinate transforms, in-basis intersections, and small actual-cygv compact regressions for matrix-basis handoff and direct-wrapper parity. McAllister runner now accepts matrix source bases for K/M flux-coordinate validation, can execute the dual flat/GV handoff in a matrix production basis, and can solve the primal KKLT branch path in an index or matrix production basis. Corrected-chamber diagnostics still transform production Kähler coordinates back to the computed CYTools index basis before selecting curves/exporting traces, so diagnostic representation cleanup remains |
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
   candidate `1/-2/3` values cannot be promoted. For the lower degree-six
   source-ray queue feeding the degree-10 path history, successful
   CMS-general-divisor checks now carry their divisor-basis and ambient-basis
   solution coefficients into the context report, preserving the chamber/input
   evidence needed by the next local intersection-certificate step. The same
   summaries now include the divisor cubic self-intersection computed from the
   corrected basis κ tensor and serialize it as a structured one-entry local
   tensor candidate. A primitive actual-`cygv` probe now checks those raw
   candidates and records the mismatch against the CMS formula values, so this
   remains candidate evidence until the full local `cygv` intersection tensor
   and chamber certificate are derived. The active-dependency context export
   now also resolves the previously anonymous non-toric active leaves into
   concrete source-ray diagnostics. The context export now includes uncovered
   source-ray toric diagnostics, and the context consumer imports them as
   source-derived known GV values after exact degree and duplicate-conflict
   checks. The latest report has 12 toric-covered active leaves, three
   source-derived leaves, two links to other missing targets, and five matching
   uncovered source rays. The formerly unresolved degree-4 leaf with ambient
   support `[(6,1),(200,1),(210,-2)]` is a two-face toric diagnostic with
   `GV=-2`, so it no longer belongs in the local-phase queue. That makes the
   next corrected-GV step local phase/intersection/chamber construction for
   the remaining higher-rank uncovered leaves, not more leaf discovery. The
   enriched active-leaf report now
   includes the matched source-ray q rows and CMS status counts; two lower
   origin-circuit leaves are primitive `neg2` candidates, while three have
   three negative primitive intersections under the oriented local q row and
   still need a certified alternate phase/chamber model. The unit-tensor phase
   probe matches only the degree-14 lower leaf and mismatches the degree-12 and
   three degree-10 leaves, so the report does not promote any new GV values.
   Non-promotable CMS rational solutions are now retained as diagnostics:
   two integral degree-10 solutions mismatch the inferred normal degree and one
   solution is nonintegral.
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
cargo test -p cyrus-core --bin mcallister_gv_context -- --nocapture
cargo build -p cyrus-core --bin mcallister_first_principles
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 timeout 600 ./target/release/mcallister_first_principles --stop-after volume --allow-invalid-ek0 --diagnose-corrected-chamber-gv --dump-corrected-chamber-gv-context /tmp/cyrus_corrected_chamber_gv_context_schema3.json
cargo build -p cyrus-core --bin mcallister_gv_context
./target/debug/mcallister_gv_context --context /tmp/cyrus_corrected_chamber_gv_context_schema3.json --out /tmp/cyrus_gv_context_schema3_ray_context_report.json
./target/debug/mcallister_gv_context --context /tmp/cyrus_corrected_chamber_gv_context_schema3.json --certify-origin-support-domains --origin-support-certificate-limit 256 --out /tmp/cyrus_gv_context_schema3_origin_support_cert_report.json
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
   index or matrix source coordinates into the selected production dual basis,
   and `--production-dual-basis` supports index or matrix internal dual bases
   for flat-direction intersections and compact GV input construction.
   Core KKLT volume, Jacobian, path-following, branch-candidate, and
   initialization-scaling APIs now accept matrix Kähler-coordinate bases.
   `mcallister_first_principles --production-primal-basis` now wires those APIs
   into the first-principles KKLT solve and branch search. The remaining basis
   work is to make corrected-chamber diagnostic exports matrix-native instead
   of transforming solved production Kähler coordinates back to the computed
   CYTools index basis for existing curve-selection and trace code.
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
   remain a diagnostic rather than a fallback. Cyrus now has a reusable
   LP-assisted supporting-face certificate search for higher-codimension
   supports, but it only promotes normals that pass exact integer verification.
   On the schema-3 corrected-chamber context, the bounded origin-support pass
   found no promoted relation-support certificate. The LP no-solution status is
   now reported as no-certificate rather than a hard diagnostic error, and a
   raised `4096`-generator guard checks the degree-10 facet-union domains
   directly: target `7` has relation/shared/union ranks `1/13/194`, target `8`
   has ranks `1/9/177`, and none is promoted to a supporting face.
   `mcallister_gv_context` now also has an opt-in exact extremal-ray probe for
   target Mori generators. It first verifies any exported exact positive
   decomposition by other degree-bounded generators, and only falls back to the
   DDM separator search when that non-extremality certificate is absent. On the
   fresh schema-3 4-214-647 context, all nine solved-t misses short-circuit:
   five are exact integer-semigroup decompositions and four are exact
   rational-cone decompositions. This is useful cone triage, but it neither
   identifies the HKTY semigroup nor promotes a corrected GV value. The
   immediate audit checkpoint is still to reconstruct the smallest exact
   semigroup/chamber context for these higher-rank origin circuits. The prior
   checkpoint of distinguishing exact
   integer-semigroup decompositions from rational-cone-only LP witnesses is now
   complete: the nine current solved-t misses split into five integer-semigroup
   LP active-generator decompositions and four rational-cone-only
   decompositions. That split is triage only, because cygv assigns GV values
   from a finite semigroup and degree-ordered subtraction history, not from a
   sparse curve relation in isolation. The path-history report now summarizes
   closure and lower-seed diamond statuses at the top level; for the two
   degree-10 targets, the bounded probe reaches the `10000` element cap with
   `7608` previous-window elements, while the small lower-seed diamonds return
   `GV=0` for one target and non-integral cygv output for the other. The
   generation-growth counters show that raising the cap to `100000` still
   truncates in the first closure layer: `720` seeds would generate `130414`
   new elements and `131135` total degree-bounded closure elements before any
   later layer. A new `--closure-generation-limit` option can stop after a full
   generated layer; with generation limit `1` and element cap `150000`, both
   degree-10 targets complete that first layer exactly with `131135` elements,
   `99317` previous-window elements, and eight predecessor differences. The
   path-history report now samples those qN predecessor pairs directly; the
   nearest candidates for both degree-10 targets are degree splits `8+2`
   (distances `4` and `3`). With the corrected-chamber covered-toric context
   included, both degree-10 targets have predecessor-pair coverage counts
   `difference_only_toric_covered=2`, `predecessor_only_toric_covered=2`,
   `neither_toric_covered=4`, and no `both_toric_covered` pairs, so the
   immediate missing history includes lower-degree non-toric terms rather than
   only the degree-10 target formulas. The sampled nearest `8+2` splits now
   also report seed membership: the degree-2 side is a toric-covered
   pair-reduced seed, while the degree-8 side is not a supplied seed, so the
   next source-history task is composite lower-degree semigroup reconstruction.
   The seed-sum sample for that degree-8 side splits it again as a
   toric-covered degree-2 reduced seed plus an uncovered degree-6 seed that does
   not survive pair reduction; tracing that degree-6 seed through the
   pair-reduction relation reaches degree-4 and degree-2 reduced seeds that
   both have toric GV values. The missing compact-GV work is therefore cygv's
   composite degree-ordered subtraction history, not leaf toric formula
   discovery for those nearest paths. Running the pair-expanded reduced-leaf
   diamond through actual `cygv` gives a 12-element explicit semigroup and
   `GV=0` for both degree-10 targets, so this still is not the valid compact
   history domain. The next broader local check,
   `--run-path-support-generators`, builds a supplied generator set from the
   union of the target, predecessor, difference, and sampled decomposition
   supports and runs that through actual `cygv`. On the same degree-10 targets
   this gives support/generator sizes `7/11` for target `7` and `6/7` for
   target `8`, and both return `GV=0`. Therefore the path-support generator
   subset is another diagnostic-only negative result, not the missing compact
   semigroup.
   A structured export now exists via
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
   history from this broad domain. This bounded probe can now be requested
   directly with `mcallister_gv_context --probe-cygv-path-history`, which avoids
   accidentally calling cygv's unbounded semigroup constructor before the
   capped source-history diagnostic.
   Running the visible lower-seed diamonds for the two degree-10 targets does
   not provide a fallback: target `7` returns `GV=0` from a 6-element diamond,
   and target `8` fails cygv series inversion with a non-integral GV from an
   8-element diamond.
   Cyrus now exposes cygv's private pair-sum seed-reduction stage only as the
   hidden source-audit helper `cyrus_core::gv::cygv_pair_reduced_seed_generators`,
   and `mcallister_gv_context` reports raw seed counts, reduced seed counts, and
   whether each target survives that reduction. On the current corrected-chamber
   context, the measured
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
   intersection tensor and chamber interpretation. The one-parameter unit
   semigroup and grading are now derived later in the local skeleton. The
   same report now groups local charge-row permutation signatures: all four
   compact-threefold-shaped misses share `[-2,-1,-1,1,3]` and are explicitly
   marked `shape_only_missing_source_derived_cygv_inputs`, with the remaining
   uncertified local `q` phase, intersection tensor, and chamber interpretation
   listed in the JSON. The core
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
   candidate local `q` matrix shape. The same skeleton now records the two
   overall local charge-basis orientation candidates: sign `-1` turns the
   target coordinate into `[1]` and places the positive unit generator in
   cygv's `neg2` fundamental-period bucket, while sign `+1` leaves the target
   at `[-1]` and puts the positive unit generator in `ignored_gt2`. This
   removes the target-coordinate/orientation ambiguity for the compact-shaped
   local supports. The report now also records the target-coordinate gcd and
   primitive direction; on the current saved context all nine origin-circuit
   misses have primitive local target coordinate, with sign `-1` giving `[1]`
   and sign `+1` giving `[-1]`. This is still not a chamber certificate and
   supplies no GV value by itself. The same orientation report now records a
   target-candidate status from the primitive coordinate and cygv's
   negative-intersection omega buckets: only targets `3`, `6`, `7`, and `8`
   have a primitive positive sign-`-1` candidate in a cygv-supported omega
   bucket, while the other five sign-`-1` candidates are positive but
   `ignored_gt2`, and every sign-`+1` candidate is a negative local coordinate.
   The report now includes these as top-level status counts:
   `target_primitive_positive_supported_by_cygv_omega_bucket=4`,
   `target_positive_but_ignored_by_cygv_omega_bucket=5`, and
   `target_not_in_nonnegative_local_semigroup=9`. The same report now keeps
   actual-call readiness separate from target eligibility; all nine local
   skeletons are still `blocked_missing_source_derived_inputs`. The primitive
   one-parameter grading vector is now source-derived as `[1]` for all nine
   skeletons, and the local `q` orientation is now source-derived as sign `-1`
   for all nine skeletons. The oriented local `q` matrix is now serialized in
   both cygv's divisor-row layout and Cyrus' wrapper layout, with
   `source_derived_oriented_q_matrix_layout=9`. The same report records
   `local_gkz_relation_includes_origin_point_requires_phase_mapping=9`, so the
   remaining local `q` phase/chamber mapping is source data, not a matrix
   transpose issue. It also records CYTools' `mori_cone_cap` origin-circuit
   retention rule: all nine candidates have negative origin coefficient and
   aggregate as `source_cytools_retains_negative_origin_coefficient=9`, so they
   are source-derived Mori-cap origin-circuit rows but not complete local
   `cygv` inputs. The one-parameter unit semigroup generator `[[1]]` is now
   source-derived for these skeletons, so the missing source-input counts are
   now the local `q` phase, local intersection tensor, and local chamber
   certificate; none is a valid actual-`cygv` call yet. The opt-in
   active-support provided-generator
   diagnostic now uses the actual Rust `cygv` crate and reports
   `computed_active_support_generators=6`, `hkty_error=3`; all six successful
   target lookups return `GV=0`, and the three errors are non-integral
   series-inversion failures. This rules out the small active-support windows
   as a valid replacement for the missing source-derived semigroup/chamber
   inputs. The exact supporting-face verifier also rejects those same windows
   as codimension-one chamber faces, with
   `active_support_not_certified_as_codimension_one_face=9`. The per-target
   JSON report now carries this status next to each active-support GV/HKTY
   result, so the six `GV=0` diagnostic lookups are visibly non-promotable. The
   corrected-chamber context export is now schema `3` and includes the full
   source-derived facet pair for each origin-circuit witness plus the
   degree-bounded Mori-ray context needed to relate candidate local semigroups
   back to ambient source curves. Old schema-`1` dumps are still readable and
   report `origin_circuit_missing_full_facet_context=9`; the regenerated
   schema-`3` report verifies `source_derived_full_facet_context=9` and records
   `2963` source-derived ambient/basis ray-context rows. The next step is to
   use those facet sets and ray supports for chamber/semigroup certification.
   The first source-support counts are now in the context report: the exact
   sparse relation support gives one generator per missing target, shared
   facet-pair neighborhoods give `12..367`, and full facet unions give
   `458..2381`, so neither naive support restriction is a valid replacement for
   the missing chamber semigroup.
   The opt-in exact-kernel certificate probe keeps this conservative: all nine
   relation-support domains and the four shared-facet domains below the
   256-generator guard now report their exact row-span rank before any
   codimension-one certificate attempt, while the broader shared/facet-union
   domains are explicitly skipped with their actual row counts. In the
   schema-`3` 4-214-647 report, all nine relation supports are rank `1` in
   dimension `214`; the checked shared-facet domains have ranks `9`, `13`, and
   `26` in dimension `214`, so none of the checked source-support domains is a
   codimension-one chamber face.
   The direct source
   read reinforces this: `cygv`
   obtains GV values from the full
   finite semigroup, HKTY alpha/beta construction, and degree-ordered
   `Li2(q_N)` subtraction history, while the McAllister paper's small-curve
   discussion is a selected-toric-curve control argument rather than a closed
   formula for every origin-circuit coefficient pattern. The next step is
   therefore to reconstruct the broader source-derived finite semigroup/path
   history for these higher-rank origin circuits, or to certify a
   flop/Weyl-continuation chain that supplies the data from another chamber,
   not to promote LP active-generator diamonds or reuse the rank-two
   local-surface machinery. The path-history report now makes this distinction
   explicit by classifying sampled predecessor/difference pairs by whether
   they have source-derived nonzero toric GV values that could enter
   `cygv`'s live `previous_qn` cache. For both degree-10 targets, the eight
   first-layer predecessor differences contain no pair where both sides are
   known nonzero lower-degree history: the split is `2/2/4` across
   known-nonzero/unknown, unknown/known-nonzero, and unknown/unknown. Raw
   predecessor decomposability is therefore not enough to explain the missing
   corrected-chamber GV value. The path-support provided-generator diagnostic
   now also looks up the sampled lower classes inside that same small `cygv`
   domain: it matches the four known nonzero toric degree-two values and gives
   six nonzero plus six zero/absent values across the unknown non-toric lower
   classes, while the target still evaluates to `GV=0`. This confirms that
   small-domain lower-history values are diagnostic only until the source
   compact semigroup/chamber history is certified. The same lookup rows are now
   classified against the schema-`3` source-ray inventory: each degree-10
   target has four known toric source-ray lookups, four non-toric source-ray
   lookups, and eight composite/non-source-ray lookups, with no lower lookup
   matching one of the nine direct target rows. The compact history problem is
   therefore broader than the final target classes themselves. The report now
   promotes the non-toric source-ray subset into a unique queue: each of targets
   `7` and `8` has two unique degree-six source rays, four sampled occurrences,
   and diagnostic small-domain GV counts `{"1":2,"-2":2}`; the source ray
   `[(54,-2),(203,1),(206,1),(209,1)]` appears in both target queues. This gives
   the next compact-GV step a concrete shared lower-history target. A fresh
   release context export now bounds auxiliary source-ray classification below
   the first missing-target degree (`9`), classifies `240` uncovered lower source
   rays, and annotates the target `7`/`8` queue with source-derived origin/CMS
   status. The `GV=1` lower rays are resolved-conifold origin circuits; the
   shared `GV=-2` ray has integral CMS divisor checks matching the inferred
   degree. Those lower source-ray values are now fed back into the qN-history
   classifier only when the CMS checks give a unique source-derived formula
   value. For both degree-10 targets, the sampled predecessor pairs now split
   evenly across source-known/unknown, toric-known/unknown,
   unknown/source-known, and unknown/toric-known histories. The path-support
   lookup matches the four source-known and four toric-known nonzero lower
   values, but the target still evaluates to `GV=0`, so the remaining compact
   history is in the degree-four/degree-eight composite classes and the
   uncertified semigroup/chamber continuation. The seed-sum diagnostic now
   records source-derived statuses too: sampled degree-eight unknowns split as
   toric degree-two plus source-known degree-six, and sampled degree-four
   unknowns split as toric degree-two plus toric degree-two. Those decompositions
   identify the finite-history inputs but do not license a multiplicative GV
   shortcut. The naive larger-domain alternative is currently too expensive:
   the actual Rust `cygv` call with all `720` degree-bounded generators through
   target degree `10` timed out after `1200s` in release mode for target `7`
   after `600s` debug timeouts for targets `7` and `8`. The report now also
   applies `cygv`'s nearest live-`q_N` predecessor rule: target `7` would start
   from toric degree-two `[(44,1),(54,1),(206,-2)]` and target `8` from toric
   degree-two `[(54,1),(203,-2)]`, each leaving an unknown degree-eight
   residual monomial. The residual trace shows both degree-eight residuals then
   start from the shared source-known degree-six ray
   `[(54,-2),(203,1),(206,1),(209,1)]` and leave those toric degree-two
   classes. The next compact-HKTY task is therefore the qN polynomial history
   for that source-derived degree-six ray and its degree-eight composites inside
   a certified semigroup/chamber. The new source-ray readiness fields sharpen
   this to a concrete input list: the shared degree-six ray has local charge
   signature `[-1,-1,-1,1,2]` and matching CMS checks, but still lacks
   `local_q_matrix_phase`, `local_intersection_tensor`, and
   `local_chamber_certificate`. Its one-parameter unit semigroup generator and
   grading vector are now source-derived. The follow-up primitive probes show
   that unit tensor normalization explains the `[-1,-2,1,1,1]` lower ray, and
   that omitting the origin/canonical divisor column explains the
   resolved-conifold-like `[-1,1,-1,1,-1,1]` lower ray. These are still
   phase/chamber diagnostics, not promoted GV values. Across all nine
   unresolved origin circuits, the origin-omitted shape aggregate now reports
   five compact-threefold-shaped no-origin rows and four CY-dimension-2 rows.
   The target-level unit-tensor phase probe matches the four CY-dimension-2
   rows in the origin-included convention and only one compact no-origin row in
   the origin-omitted convention; the other four compact no-origin rows still
   mismatch their formula candidate sets.
   The follow-up opt-in local integer tensor scan confirms this is not just a
   missing scalar normalization: with tensor values scanned over `-8..8`, four
   of the nine targets match their CMS formula candidate at tensor value `1`,
   but five targets have no matching integer tensor in that range. All four
   matches remain marked uncertified because the local q-phase, intersection
   tensor, and chamber certificate are not source-derived yet.
   The first-principles context export
   now serializes the nonzero divisor-basis and ambient-basis coefficients from
   successful CMS-general-divisor solves, preserving the source-derived divisor
   input needed by the next chamber/intersection certificate step.
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
   instanton sum. The schema-`3` context report now exposes that readiness
   directly: all nine remaining targets have CMS-style formula candidates, but
   the 14 exact divisor-intersection checks all fail with
   `cms_general_divisor_no_rational_divisor_solution`, so no source-derived
   shrinking divisor certificate is currently available. The default schema-3
   context report now verifies the target-ray side without enabling the
   expensive separator search: five targets are
   `not_extremal_by_exact_integer_semigroup_decomposition`, four are
   `not_extremal_by_exact_rational_cone_decomposition`, and all nine active
   supports remain `active_support_not_certified_as_codimension_one_face`.
   This rules out treating the missing target rays themselves as certified
   chamber walls; the remaining chamber work has to recover a source-derived
   face/continuation history. The same report now aggregates target statuses
   directly. With `--run-integer-diamonds`, the five integer-semigroup
   decompositions split into three `computed_integer_decomposition_diamond`
   rows with target GV `0` and two `hkty_error` rows with non-integral compact
   HKTY output, while the four rational-cone targets remain skipped. This
   keeps the decomposition-diamond path as a diagnostic, not a promoted GV
   computation. The active decomposition generators now aggregate as 12
   toric-covered leaves, three source-derived GV leaves, two matching other
   missing targets, and five matching uncovered source-ray leaves. This shifts
   the actionable GV work to certifying those lower source-ray leaves and their
   chamber history. The report now
   includes `active_decomposition_unresolved_source_leaf_sample`, which lists
   the exact unresolved dependency leaves and their parent occurrences. In the
   current context this sample contains two links back to other missing targets
   and five not-yet-covered source rays of degrees 10, 10, 10, 12, and 14. The
   previous degree-4 source ray is now classified from an exported two-face
   toric diagnostic with `GV=-2`. The sample also records ambient
   origin-relation patterns, showing that the remaining uncovered source-ray
   dependencies all have origin-circuit charge shapes that need local
   chamber/GV treatment. The first-principles export now preserves every
   origin-circuit facet-pair witness in `origin_circuit_witnesses`, in addition
   to the compatibility `origin_circuit_first_witness` field, so multi-witness
   source rows such as the degree-14 `[-1,-2,1,1,1]` dependency keep the
   alternate local chamber data needed by the next certificate step. The
   context diagnostics now consume all serialized witnesses when forming
   support-domain and facet-context summaries, so multi-witness rows expose
   mixed facet context rather than inheriting the first witness's status. A
   fresh all-witness release export now confirms that the five multi-witness
   missing targets all share their local relation across witnesses; the
   remaining blocker is certifying the facet/chamber phase and local
   intersection tensor for that relation.
   A follow-up per-witness domain diagnostic now checks each serialized
   origin-circuit witness separately. The 14 witness domains have relation
   support rank `1`; shared-facet domains range from `12` to `367` generators
   and ranks `9`, `13`, `26`, or `146`. With the existing LP/exact certificate
   path and a `512` generator guard, no relation or shared-facet witness domain
   certifies a supporting face. The two degree-10 facet-union domains under the
   same guard also fail certification, while larger facet unions remain skipped
   explicitly. This rules out hiding the chamber certificate in one selected
   facet-pair witness. The same per-witness report now classifies each domain's
   generators by source status. The shared-facet domains aggregate to `1727`
   toric-covered generator occurrences and `671` source-derived GV generator
   occurrences, but still contain `228` source rays without a derived GV value,
   `20` occurrences of other missing targets, and `10` uncovered-source-ray
   hits. Facet unions contain a much larger unresolved source-ray bucket, so a
   reusable semigroup handoff still needs lower source-ray closure rather than
   a small all-known witness domain. The report now exposes that queue
   directly: across relation, shared-facet, and facet-union domains there are
   `2017` unique non-known generators over `9611` occurrences, consisting of
   the nine missing targets, `41` uncovered-source-ray matches, and `1967`
   degree-bounded source rays with no toric/source-derived GV value yet.
