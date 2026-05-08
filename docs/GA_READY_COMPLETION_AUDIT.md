# GA-Ready Completion Audit

This is a prompt-to-artifact audit for the active objective:

Build a GA-ready Cyrus implementation of the McAllister/DKMM compactification
pipeline. McAllister 4-214-647 is a validation target, not the implementation
source. Cyrus may use ancillary `.dat` files as declared inputs or checkpoints,
but every reusable geometry and physics quantity needed in a landscape search
must be computed from upstream inputs.

Status: **not complete**. The current codebase has important no-replay guards
and several source-derived stages, but the corrected GV/Kahler instanton layer
and full potent-ray row generation are still open.

Course correction: the compact GV implementation boundary is the existing
`cygv` Rust crate. Cyrus should compute CYTools-equivalent inputs
(`generators`, grading vector, curve-basis matrix, and intersection numbers)
and call `cygv`; it should not reimplement `cygv`'s compact HKTY internals
unless a concrete missing API or correctness bug is isolated. The local CKYZ
rank-two surface code is a diagnostic/validation path for noncompact local
models, not the production compact GV engine.

## Concrete Deliverables

The objective breaks down into these success criteria:

1. Declare which McAllister files are model inputs, validation checkpoints, and
   replay-only downstream outputs.
2. Run a first-principles pipeline without loading replay-only outputs such as
   `kahler_param.dat`, `corrected_kahler_param.dat`, `cy_vol.dat`,
   `corrected_cy_vol.dat`, or precomputed GV/intersection data.
3. Compute geometry from declared inputs: polytopes, duals, triangulations,
   divisor bases, GLSM/curve basis, intersection numbers, Mori/Kahler data.
4. Compute GV data from source algorithms, not from saved `*_gv.dat` files.
5. Compute flux flat direction, racetrack quantities, KKLT target data,
   corrected Kahler moduli, corrected volume, and `V0`.
6. Compare each stage against CYTools and/or McAllister checkpoints with
   standalone tests.
7. Keep remaining discrepancies explicit; do not mask them with overrides,
   fallback values, branch seeding, or downstream replay.
8. Leave a reusable code path suitable for GA search on new polytopes,
   triangulations, fluxes, and moduli choices.

## Prompt-To-Artifact Checklist

| Requirement | Evidence inspected | Current status |
| --- | --- | --- |
| Allowed inputs are separated from checkpoints and replay-only outputs | `docs/MCALLISTER_DATA_POLICY.md`; `ARTIFACT_POLICIES` in `crates/cyrus-core/tests/mcallister_e2e/stage0_data_integrity.rs` | Implemented policy. Declared inputs are `points.dat`, `heights.dat`, `K_vec.dat`, `M_vec.dat`, `kklt_basis.dat`, and `target_volumes.dat`. Replay-only outputs are explicitly classified. |
| First-principles mode rejects JSON fixture fallback | `crates/cyrus-core/tests/mcallister_e2e.rs`; `enforce_modes` in `crates/cyrus-core/src/bin/mcallister_first_principles.rs` | Implemented guard. `CYRUS_ALLOW_FIXTURES` is rejected with `CYRUS_FIRST_PRINCIPLES`; the runner refuses silent fixture fallback. |
| Replay-only corrected Kahler parameters are not loaded in the normal solve | `mcallister_first_principles` volume path around `--allow-downstream-kahler`; warning that `corrected_kahler_param.dat` remains validation-only replay | Implemented guard. Replay is behind an explicit flag, and corrected-chamber diagnostics are rejected when replay is enabled. |
| Replay-only volumes are not production inputs | `compare_against_dat` reads `corrected_cy_vol.dat` only after computing `V_string`; `docs/MCALLISTER_DATA_POLICY.md` | Mostly implemented. Comparison is post-computation, but the runner still tolerates a residual up to `0.1`, so this is a diagnostic gate, not proof of exact reproduction. |
| Geometry, basis, intersections, Mori/Kahler data are computed from upstream inputs | `docs/GA_READY_PIPELINE_AUDIT.md`; `mcallister_first_principles`; stage 2/3/5 tests referenced there | Largely implemented, but not a final completion certificate. Generic matrix-basis GV inputs now include Mori projection, q-matrix construction, curve-basis construction, and in-basis intersection tensor pullback. `mcallister_first_principles --dual-basis` accepts index or matrix source bases for K/M flux-coordinate transforms, and `--production-dual-basis` now carries an index or matrix internal dual divisor basis through flat-direction intersections and compact GV inputs. Core KKLT has generic divisor-basis volume, Jacobian, path-following, branch-candidate, and initialization-scaling APIs for matrix Kähler coordinates. `mcallister_first_principles --production-primal-basis` now carries an index or matrix internal primal divisor basis through KKLT path following, branch initialization, branch search, and branch GV coverage by transforming solved Kähler points back to the computed CYTools index basis for existing chamber/GV diagnostics. Remaining matrix-basis work is diagnostic/export cleanup, not the main corrected-GV blocker. |
| CYTools/cygv compact GV wrapper contract is understood and mirrored | `docs/CYGV_AUDIT.md`; direct source read of `reference/cytools/src/cytools/calabiyau.py::_compute_gvs_gws`; downloaded `cygv-0.1.2` source; `mcallister_gv --min-points 20000`; `stage5_mirror_gv_checkpoint_matches_cygv_min_points` | Implemented for the generic compact handoff boundary. CYTools builds inputs, then delegates to `cygv.compute_gv`; no hidden Python GV algorithm remains. The 4-214-647 mirror-side min-points run computes 10556 ambient invariants and matches all 5177 `dual_curves.dat` / `dual_curves_gv.dat` checkpoint rows without loading those rows as inputs, and the command is now covered by an opt-in heavy e2e regression. |
| High-dimensional small toric curve selection is computed, not read | `docs/GA_READY_PIPELINE_AUDIT.md`; stage 5 small-curve tests | Implemented for the checkpoint rule. Pair pruning matches `small_curves.dat`; finite-semigroup pruning removes five additional input-chamber curves, so both policies remain explicit. A fresh no-replay finite-semigroup run shows this is not the corrected-volume fix: the checkpoint-t corrected chamber prunes `556` subcutoff curves to `400` with zero toric-missing rows, and the solved-t corrected chamber prunes `561` to `399` with zero toric-missing rows, but `V_string=4711.504343103075` still has a corrected-volume residual of about `0.071843868`. |
| High-dimensional small toric GV values are computed, not read | `docs/CYGV_AUDIT.md`; `docs/LOCAL_TORIC_GV_SOURCE_MAP.md`; stage 5 small-toric GV tests | Partially implemented. Covered toric two-face/origin-circuit formulas match the checkpoint, including the source-derived resolved-conifold `(-1,-1,1,1)` pattern. Remaining corrected-chamber misses are unresolved higher-rank origin circuits. |
| Potent-ray GV rows are computed from local source data | `docs/POTENT_RAY_SOURCE_READ.md`; `docs/CKYZ_SERIES_DOMAIN_AUDIT.md`; `potent_ray_affine_circuits` tests | Diagnostic only. Cyrus reconstructs rank/volume/slope diagnostics and first-four GV entries for 395 rank-two CKYZ rows. It also has all-ten checks for selected F0/F1 families and ignored release all-ten regressions for both polygon-5 source directions `[4,3,2]` and `[3,2,2]`. These rank-two CKYZ rows should not drive compact GV implementation work; they are local noncompact validation checks. |
| Local CKYZ diagnostics stay separate from compact `cygv` | `docs/CYGV_AUDIT.md`; `docs/CKYZ_SERIES_DOMAIN_AUDIT.md`; direct source read of `cygv-0.1.2/src/{semigroup,fundamental_period,instanton,series_inversion}.rs` | Implemented as a bounded diagnostic path. Compact GV uses `compute_gv_invariants*` wrappers around `cygv`. The local CKYZ code exists because compact `cygv` rejects naive local surface charge matrices with effective CY dimension below three; it is not a replacement for the crate. |
| Flux flat direction and `eK0` are computed | `docs/GA_READY_PIPELINE_AUDIT.md`; stage 4 tests | Implemented for the audited McAllister basis paths, with remaining general-basis work noted separately. |
| Racetrack `g_s` and `W0` are computed | `mcallister_first_principles`; `mcallister_racetrack`; stage 5 racetrack tests | Implemented path exists. Checkpoints such as `g_s.dat` and `W_0.dat` are comparisons, not declared inputs. |
| Corrected KKLT target and corrected Kahler solve are exact | `mcallister_first_principles` KKLT path; `docs/GA_READY_PIPELINE_AUDIT.md`; `docs/CYGV_AUDIT.md` | Not complete. Cyrus computes a no-replay corrected Kahler vector and volume, but the corrected-chamber GV target correction has an unresolved residual. |
| Final corrected volume and `V0` reproduce McAllister exactly | `compare_against_dat`; `docs/GA_READY_PIPELINE_AUDIT.md` | Not complete. The no-replay volume is close but not exact; the residual is documented as an unresolved instanton/chamber discrepancy. |
| Remaining discrepancies are explicit, not hidden by fallbacks | `docs/MCALLISTER_DATA_POLICY.md`; `docs/GA_READY_PIPELINE_AUDIT.md`; `docs/CYGV_AUDIT.md`; `docs/CKYZ_SERIES_DOMAIN_AUDIT.md` | Implemented as documentation and runner warnings/errors. The code should continue to fail loudly rather than substituting saved outputs. |
| GA-ready path works on new candidates | Crate APIs and runner boundaries; current audit docs | Not complete. The current path is partly reusable, but corrected-chamber GV/Kahler gaps, rank-four potent-ray contexts, full matrix-basis handling, and finite GV table/chamber certification remain blockers. |

## Source Audit Result

The CYTools GV boundary is narrow:

- `CalabiYau._compute_gvs_gws` computes in-basis intersection numbers, the
  no-origin curve-basis matrix, Mori-cap generators, optional lattice-point
  augmentation, and a grading vector.
- It then calls `cygv.compute_gv` or `cygv.compute_gw`.
- The `Invariants` wrapper can re-express output charges in another basis, but
  it is not another GV algorithm.

The `cygv==0.1.2` HKTY path is:

1. Build a finite affine semigroup with `Semigroup::with_max_degree`,
   `with_min_elements`, or `from_data`.
2. Build `PolynomialProperties`, whose monomial map is the finite truncation
   contract. Polynomial products whose exponent sum is absent are dropped.
3. Compute the fundamental period and first/second logarithmic derivatives.
4. Compute `alpha`, `beta`, `beta - alpha alpha`, contracted instanton
   polynomials, and `exp(alpha_i)` / `exp(-alpha_i)`.
5. Run degree-ordered `series_inversion::invert_series`, reading candidate GV
   coefficients, materializing `q_N`, subtracting `Li2(q_N)`, and rolling the
   `previous_qn` cache.

This confirms the production compact direction: keep using the actual `cygv`
crate and focus Cyrus work on producing the correct inputs and basis
transforms. The CKYZ z-domain machinery is useful only for local noncompact
surface diagnostics where compact `cygv` is not the right model.

## Last-12-Hours Assessment

The last 12 hours produced substantial artifacts, but also documented rejected
optimization attempts.

Useful progress:

- exact corrected-chamber diagnostics and stable-Weyl/separator primitives;
- finite GV/nilpotent-ray classification helpers;
- CYTools/cygv source audits for the compact and local GV paths;
- CKYZ support-domain profiling and predicted-domain APIs;
- indexed CKYZ finite-domain series operations;
- z-residual, cover-weight ordered CKYZ extraction;
- cygv-shaped rolling `previous_qn` history for local CKYZ extraction;
- profiling that isolated the current bottleneck to direct `q_delta`
  materialization.

Churn that should not continue:

- reimplementing compact `cygv` HKTY internals in Cyrus;
- recursive per-coefficient Li2 demand;
- hybrid direct/reused Li2 extraction;
- shared direct-coordinate demand graph layered on the current domain.

Each CKYZ attempt passed small regressions but failed or regressed on the
McAllister `[4,3,2]`, N=10 diagnostic. That diagnostic should stay bounded
until the compact `cygv` handoff, corrected-chamber semigroup certification,
and matrix-basis pipeline are no longer open.

## Open Blockers

1. **Corrected-chamber GV target corrections.** Nine current solved-t misses
   remain higher-rank origin-circuit Mori generators; small active-support
   windows and LP witnesses are diagnostic only. A fresh read of the latest
   Python `compute_kklt_iterative.py` confirms Cyrus' target-sign convention,
   but that script maps `c_i` onto `basis.dat` and does not implement the mixed
   `basis.dat`/`kklt_basis.dat` corrected-target vector, so it is not evidence
   that Python had reproduced this layer without replay. The schema-4
   corrected-chamber context export is now accepted by `mcallister_gv_context`;
   bounded target-7 and target-8 path-history reports aggregate the closest
   known `q_N` residual as the same source-derived degree-6 predecessor plus
   toric degree-2 differences. The residual difference side is now a
   first-class report queue as well: target 7 leaves known toric
   `[(44,1),(54,1),(206,-2)]`, target 8 leaves known toric
   `[(54,1),(203,-2)]`, and both point back to the same source-derived degree-6
   predecessor. This makes the remaining degree-8 composite qN-history gap
   explicit instead of hiding it inside occurrence samples. The unfiltered
   bounded report now promotes
   closest-known residual source predecessors into a first-class queue: two
   unique source-derived predecessors appear across three occurrences, target
   7/8 share the degree-6 ray, target 6 has the analogous degree-4 ray, and all
   have shared origin-circuit witness relations plus
   `source_derived_full_facet_context`. They are now blocked specifically on
   extracting/certifying the source-derived `local_intersection_tensor` and
   `local_chamber_certificate`. The source-derived GV-history importer now
   requires that same witness/facet context before accepting a CMS-derived
   scalar GV value, so shape-only CMS matches remain diagnostic rather than
   known qN history.
2. **Exact corrected KKLT volume.** The no-replay path computes a corrected
   volume but still has a documented instanton/chamber residual.
3. **Generic matrix-basis pipeline.** Matrix divisor-basis primitives now cover
   compact `cygv` inputs including intersections, `--dual-basis` matrix source
   coordinates can transform K/M fluxes into the production dual basis, and
   `--production-dual-basis` can carry a matrix internal dual basis through
   flat-direction and compact-GV handoff. Core KKLT matrix-basis volume,
   Jacobian, path-following, branch-candidate, and initialization-scaling APIs
   now exist. The first-principles runner also accepts
   `--production-primal-basis` and uses that index or matrix internal primal
   basis for KKLT path following, branch initialization, branch search, and
   branch GV-coverage ranking. Existing chamber/GV diagnostic internals still
   transform production Kähler coordinates back to the computed CYTools index
   basis before selecting curves or exporting traces, so matrix-basis diagnostic
   cleanup remains, but the reusable KKLT production solve is no longer
   vector-only.
4. **Compact GV semigroup/face certification.** Missing corrected-chamber GV
   classes need a source-derived compact or certified face semigroup that can
   be handed to `cygv`, not an unproven local CKYZ substitute. Cyrus now has an
   LP-assisted supporting-face search whose successful outputs are exact
   integer certificates, but the current schema-3 McAllister origin-support
   pass still finds no promoted certificate for the remaining relation
   supports. The LP diagnostic now treats HiGHS `NoSolutionFound` and other
   finite solver statuses such as `Unknown` as no-certificate outcomes rather
   than hard errors, because only exact integer-verified normals may be
   promoted. Raising the origin-support guard to `4096` for the two degree-10
   targets checks their larger source domains without promotion: target `7`
   has relation/shared/union ranks `1/13/194`, target `8` has ranks
   `1/9/177`, and all checked domains return
   `origin_support_lp_no_certificate_*`.
   The context tool now also has an opt-in exact extremal-ray probe. For the
   current schema-3 context it verifies that all nine remaining targets are
   non-extremal in the finite degree-bounded cone by exact decomposition
   certificates: five integer-semigroup decompositions and four rational-cone
   decompositions. This is not a substitute for the missing compact/chamber
   semigroup, because GV values are not determined by finite-cone extremality
   alone. The path-history probe now aggregates closure/lower-seed status
   counts; the actual-local-cygv lower-diamond report has six raw lower-seed
   diamonds compute `GV=0`, three raw lower-seed domains fail with non-integral
   cygv output, and all nine pair-expanded reduced-seed diamonds compute
   `GV=0`. The generation counters show the first closure layer alone would add
   `130414` elements from `720` seeds, so the first visible decomposition is
   not a valid replacement for the full HKTY history.
   The path-history probe can now stop after a full generated layer; for the
   two degree-10 targets the first layer is complete at `131135` elements and
   `99317` previous-window elements, giving a bounded source-history slice for
   the next qN-history analysis. It also samples the actual predecessor pairs;
   both degree-10 targets have eight first-layer predecessor differences, with
   nearest degree split `8+2`. After exporting corrected-chamber toric-covered
   rows into the context, those predecessor pairs split as
   `difference_only_toric_covered=2`, `predecessor_only_toric_covered=2`, and
   `neither_toric_covered=4`, with no pair where both sides are toric-covered.
   This keeps the next compact-GV task focused on lower-degree non-toric
   history, not just the nine target classes. The nearest `8+2` samples now
   show that the degree-2 side is a toric-covered pair-reduced seed, while the
   degree-8 side is not a supplied seed, so that lower-degree history is
   composite rather than a single uncovered Mori generator. The source-ray
   queue now preserves successful CMS-general-divisor solution coefficients for
   these lower classes, so the next chamber/intersection step can use actual
   divisor-basis and ambient-basis evidence rather than only scalar CMS status
   counts. The context report also contracts those solved divisors with the
   corrected basis intersection tensor and records their cubic
   self-intersections as structured one-entry local-intersection tensor
   candidates, with a status that still requires the local phase and chamber
   certificate before promotion. The report now also runs a primitive actual
   `cygv` probe against those raw cubic candidates; the target-`7` lower-ray
   summaries explicitly mismatch (`-6` versus expected `-2`, and `0` versus
   expected `1`), so the cubic evidence is guarded against promotion as a
   fitted local tensor. A bounded unit-tensor comparison sharpens this: the
   `[-1,-2,1,1,1]` lower ray gives the expected `-2` with unit tensor `1`, but
   the `[-1,1,-1,1,-1,1]` lower ray still gives `0` against expected `1`.
   An explicit origin-omitted phase diagnostic then drops the origin/canonical
   divisor column before the unit-tensor `cygv` call: this matches the
   resolved-conifold-like row (`1` against expected `1`) and rejects the
   `[-1,-2,1,1,1]` row with the compact-HKTY dimension guard. This is evidence
   for the local `q` phase convention, not a chamber certificate.
   The all-target origin-omitted shape aggregate originally split the nine
   unresolved origin circuits into five dimension-three no-origin rows and four
   CY-dimension-2 rows. The phase selector now labels those five no-origin rows
   as non-Calabi-Yau dimension-three diagnostics rather than compact CY
   threefold handoffs, because their charge sums are nonzero. A target-level
   unit-tensor phase probe then shows that the four CY-dimension-2 rows match
   their formula candidate set only in the origin-included convention, while
   only one of the five no-origin diagnostic rows matches after origin
   omission. The other four no-origin rows still mismatch, so a single
   origin-column rule is not enough.
   The first seed-sum
   sample for the degree-8 side splits it as another toric-covered degree-2
   reduced seed plus an uncovered degree-6 seed that does not survive pair
   reduction; that degree-6 seed pair-reduces into degree-4 and degree-2
   reduced seeds that both have toric GV values. The next compact-GV task is
   therefore reproducing the composite semigroup subtraction history over these
   leaves. A pair-expanded reduced-leaf explicit semigroup has only 12 elements
   and returns `GV=0` for both degree-10 targets, so it is not the valid compact
   history domain. A broader path-support provided-generator probe, formed from
   the target plus sampled predecessor/difference/decomposition supports and
   still executed through the actual Rust `cygv` crate, also fails as a
   replacement: target `7` uses support size `7` with `11` generators and
   target `8` uses support size `6` with `7` generators, and both return
   `GV=0`. This keeps the missing object at the certified compact/chamber
   semigroup and degree-ordered subtraction history level, not a small support
   subset. The path-history probe now also separates raw predecessor existence
   from `cygv`-style live qN history: for both degree-10 targets, the eight
   sampled predecessor differences split as two known-nonzero toric plus
   unknown pairs, two unknown plus known-nonzero toric pairs, and four
   unknown/unknown pairs. No sampled pair has both sides certified as known
   nonzero lower-degree history, so the next compact-GV task is still the
   unknown lower-degree non-toric history or a certified chamber-continuation
   source for it. The report now also distinguishes scalar GV evidence from
   compact `q_N` polynomial materialization: known toric/source scalar values
   do not mean the compact mirror-map `q_N` polynomial has been exported from
   cygv's history. Cyrus now has a vendored `cygv` 0.1.2 trace API for explicit
   and provided-generator domains, checked on the quintic `2875` case, so the
   observability blocker is reduced to finding the right corrected-chamber
   compact semigroup/history domain. The same small path-support `cygv` domain
   now reports per-predecessor lookups and qN-trace status: for target `7`, it
   materializes 11 qN polynomials with bounded term samples, has qN for ten
   sampled nonzero lower lookups, and records six zero/absent lower lookups
   with no qN, but the target remains `GV=0` with no target qN polynomial.
   Target `8` repeats the same scalar/qN-status split with 6 materialized qN
   polynomials and no target qN polynomial. These lower-class values are
   therefore diagnostic artifacts of the small domain, not promotable compact
   GV history. A reduced two-target aggregate shows 17 qN-polynomial
   occurrences, 11 unique lower curves, and 6 curves shared by both degree-10
   target-support domains; 4 unique lower curves have domain-dependent qN term
   counts.
   The report now also aggregates the uncovered source-ray subset as a unique
   queue. For each of targets `7` and `8`, the queue has two unique degree-six
   source rays, four sampled occurrences, and diagnostic small-domain GV counts
   `{"1":2,"-2":2}`. One degree-six ray,
   `[(54,-2),(203,1),(206,1),(209,1)]`, appears in both target queues, so the
   next compact-GV task is shared lower-degree history rather than nine
   independent final target rows.
   The first-principles context export now bounds this auxiliary source-ray
   classification below the first missing-target degree (`9` here), keeping the
   heavy dump usable while classifying `240` uncovered lower source rays. The
   target `7`/`8` queue now records source-derived origin/CMS status: the
   `GV=1` lower rays are resolved-conifold origin circuits, and the shared
   `GV=-2` ray has integral CMS divisor checks matching the inferred degree.
   The path-history classifier now consumes those source-derived lower-ray
   values under a conservative uniqueness rule. Rerunning targets `7` and `8`
   changes the sampled predecessor qN-history counts to two
   source-known/unknown pairs, two toric-known/unknown pairs, two
   unknown/source-known pairs, and two unknown/toric-known pairs for each
   target. The small path-support `cygv` lookup matches the four source-known
   nonzero lower values and four toric-known nonzero lower values, but the
   target still returns `GV=0`, so this is not a final GV correction.
   The seed-sum diagnostics now carry the same source-derived status and show
   that the sampled degree-eight unknowns are sums of a toric degree-two seed
   and a source-known degree-six seed, while the sampled degree-four unknowns
   are sums of toric degree-two seeds. These are compact-history inputs, not a
   multiplicative shortcut to composite GV values.
   A direct all-degree-bounded `cygv` call is also not currently practical:
   feeding the `720` generators of degree at most `10` to target `7` timed out
   after `1200s` in release mode after earlier debug runs for targets `7` and
   `8` timed out at `600s`. This leaves certified semigroup/chamber reduction
   as the next required implementation, not brute-force enlargement.
   The same support-overlap provided-generator path can now be run with
   `--trace-support-overlap-qn`; when the `cygv` call completes, the report
   records qN polynomial count, target qN status, target qN term count, and a
   bounded qN term sample. A tiny quintic regression verifies that this traces
   the actual `cygv` history and recovers the degree-one `2875` qN polynomial.
   This improves observability but does not change the broad-domain timeout
   conclusion above. Fresh schema-4 target-specific overlap traces now also
   show that the smaller degree-10 support-overlap windows are not certified
   histories: thresholds `1..4` for target `7` (`64`, `41`, `24`, `10`
   generators) and target `8` (`48`, `40`, `18`, `10` generators) all fail
   inside `cygv` series inversion with non-integer GV output. The path-support
   domains do materialize lower compact `q_N` polynomials, including the shared
   ready source-derived degree-six predecessor, but both target domains still
   return target `GV=0` with no target `q_N` polynomial. The path-support
   report now explicitly lists lower `q_N` polynomials that contain the target
   monomial: for both target `7` and target `8`, these are the toric degree-two
   predecessor, the shared source-derived degree-six row, and the unknown
   degree-eight residual, with target-term coefficients `-2`, `-1`, and `1`.
   The coefficients are now paired with exact `Li2(q_N)` terms traced inside
   `cygv`: Li2 target coefficients `-5/2`, `-1`, and `1`, with pivot
   subtraction coefficients `5`, `0`, `-2` for target `7` and `-10`, `2`, `2`
   for target `8`. The report now also records the target-degree GV candidate
   read directly from cygv's mutable instanton polynomial. In both small
   path-support domains that readout is exactly zero:
   `path_support_target_instanton_coefficient=0`,
   `path_support_target_gv_candidate=0`, and
   `path_support_target_gv_coefficient_status=integer_zero_or_absent_gv`.
   Thus the current path-support domains are internally consistent zero-GV
   histories after lower `Li2(q_N)` subtraction; the remaining blocker is the
   certified corrected-chamber semigroup/history, not qN trace visibility. The
   report now exports all path-support GV coefficient-readout status counts as
   well: target `7` has `79` readouts (`11` nonzero, `25` zero/absent, `43`
   missing), and target `8` has `33` (`6`, `16`, `11`), giving a direct audit
   of the small-domain composite decisions.
   The path-history report now also mirrors `cygv`'s nearest live-`q_N`
   predecessor selection. The closest certified qN predecessor for target `7`
   is the toric degree-two class `[(44,1),(54,1),(206,-2)]` with an unknown
   degree-eight residual; for target `8` it is `[(54,1),(203,-2)]`, again with
   an unknown degree-eight residual. So the next missing object is the
   certified qN history for those degree-eight residuals.
   The residual trace now shows both degree-eight residuals choose the shared
   source-known degree-six ray `[(54,-2),(203,1),(206,1),(209,1)]`, leaving the
   same toric degree-two classes above. This confirms the remaining gap is qN
   polynomial history in the certified compact/chamber semigroup, not scalar GV
   values for the lower leaves.
   The latest residual-subtarget path-support probe makes that statement more
   concrete without closing the blocker. For target `7`, the degree-eight
   residual itself runs through a 7-support, 8-generator diagnostic domain and
   returns `GV=-2` after materializing 10 qN polynomials. For target `8`, the
   analogous residual runs through a 6-support, 6-generator domain and also
   returns `GV=-2` after materializing 6 qN polynomials. The original
   degree-ten target-support probes in the same reports still return `GV=0`
   with no target qN polynomial, so the new result shows that the residual can
   be observed as a lower subtarget but does not certify the McAllister
   corrected-chamber target history. The parent-vs-subtarget comparison now
   makes the domain dependence explicit: target `7`'s residual qN has three
   terms in the degree-ten parent path-support domain but one as a residual
   subtarget; target `8` has two parent-domain terms but one residual-subtarget
   term. Both complete samples are marked `different_qn_term_counts`, so the
   remaining work is constructing the certified chamber history whose qN
   polynomial includes the right ambient terms, not finding a scalar residual
   GV value. The parent-only terms now classify as the missing target itself
   for target `8`, and as the missing target plus one generated degree-ten
   non-source-ray term for target `7`; the cancellation is therefore occurring
   at the degree-ten chamber-history layer. The offsets from the residual to
   those parent-only terms are all known degree-two toric classes with `GV=-2`,
   so the remaining history problem is how the certified chamber couples a
   source-known degree-eight residual to sibling toric offsets, not how to
   assign those lower offsets. The latest regenerated reports also look up each
   parent-only exponent in the parent path-support `cygv` run. The
   generated target-`7` side term is nonzero in that parent domain (`GV=-2`,
   one qN term, `integer_nonzero_gv`), while the missing-target-shaped
   parent-only terms for targets `7` and `8` are both parent-domain `GV=0` with
   no qN polynomial and `integer_zero_or_absent_gv` coefficient readouts. The
   small domains therefore expose qN-history coupling but still do not compute
   the missing target GV. The qN-shape trace further shows the generated
   target-`7` side term is only an `identity_single_term_qn_polynomial`, so it
   is not concealing a lower-history expansion inside its own qN polynomial.
   The context report now aggregates the bounded lower-seed decomposition
   status for these parent-only terms: target `7` has two
   `found_lower_seed_decomposition` terms whose tiny decomposition diamonds
   both compute `GV=0` with no target qN polynomial, while target `8` has one
   `found_lower_seed_decomposition` term whose tiny diamond hits an HKTY error.
   The target `8` bounded-diamond error is now counted at report top level and
   the detailed sample pins it to `[(203,-3),(206,1),(209,1)]` with GV
   candidate `4/3`.
   This keeps the next work focused on the parent compact-history domain, not
   on rediscovering lower scalar leaves.
   The target-monomial source comparison is now sharper: bounded lower-leaf
   diamonds reproduce the local source scalar GV and one-term qN polynomial,
   but the parent domain has extra qN monomials. Those parent-only terms are
   now classified in-context. For targets `7` and `8`, they include the final
   missing target ray itself as a raw/reduced seed, plus non-source composite
   degree-eight/degree-ten terms. A diagnostic semigroup built only from the
   observed parent qN term support is still not a certificate: target `7`
   source rows fail with non-integer HKTY output, while target `8` source rows
   compute the scalar source GV values but still materialize only one qN term
   and mismatch the parent qN term counts. Thus the parent qN term set itself
   is not the missing compact/chamber semigroup; the broader certified
   parent-domain history is still required.
   A seed-expanded variant, which closes those observed parent qN terms through
   bounded lower-seed diamonds before calling actual `cygv`, also remains
   diagnostic-only. Target `7`'s shared degree-six source expands to `17`
   elements and still fails HKTY integrality; its degree-eight composite
   expands to `12` elements and computes `GV=-2` but still has
   `different_qn_term_counts`. Target `8`'s expanded degree-six and
   degree-eight supports expand to `10` and `8` elements and both fail HKTY
   integrality. The missing object is therefore not just a bounded closure of
   observed parent qN terms.
   The generated-term sample from the one computed expanded case is also
   specific: target `7`'s degree-eight composite now materializes the identity
   term and the final target-`7` monomial, but still omits the sibling
   non-source degree-ten cancellation monomial
   `[(44,1),(203,-1),(206,-1),(209,1)]` with coefficient `-1`.
   The parent-only offset classifier now shows these degree-ten monomials are
   source-plus-known-degree-two-toric offsets: target `7`'s final-target
   offset is `[(44,1),(54,1),(206,-2)]` with `GV=-2`, the omitted sibling
   offset is `[(54,1),(203,-2)]` with `GV=-2`, and target `8` uses that same
   `[(54,1),(203,-2)]` offset. The next missing object is the chamber history
   that couples a degree-eight composite to several known degree-two toric
   offsets in one parent qN polynomial.
   A provided-generator diagnostic using the source decomposition plus those
   known offsets now reproduces the parent qN term multiset for target `7`'s
   shared degree-six and degree-eight source rows, and for target `8`'s shared
   degree-six row. Target `8`'s degree-eight row still fails with a non-integer
   `cygv` invariant. These matching rows are explanatory candidates, not
   production inputs, until a chamber/face certificate supplies the generator
   set from geometry.
   The same offset-generator candidates now run through the supporting-face
   verifier and all remain uncertified: target `7` reports LP no-certificate
   statuses at ranks `4` and `3` in dimension `214`, and target `8` reports LP
   no-certificate at rank `3` for both the matching degree-six row and the
   failing degree-eight row.
   The LP diagnostic now splits those no-certificate statuses by phase. Fresh
   target `7` and `8` reports show all four offset-generator candidates have
   `lp_search_status=lp_no_certificate`,
   `exact_kernel_status=no_certificate`, and
   `aggregate_status=lp_cutting_round_limit` at the default `64` cutting
   rounds; all `16` bounded anchor LP attempts per candidate also return the
   cutting limit, with zero anchor LP real-normal solutions. The context
   binary now exposes `--supporting-face-lp-cutting-rounds` and
   `--supporting-face-lp-anchor-attempts` for this diagnostic. Raising target
   `7` and `8` to `256` cutting rounds makes both degree-eight aggregate LPs
   return `lp_no_solution`, while the shared degree-six source domain still
   hits the aggregate cutting limit. Anchor attempts still produce no real
   normal: target `7` splits `13/3` and `11/5` between `lp_no_solution` and
   `lp_cutting_round_limit`; target `8` splits `11/5` for both rows. This
   still does not prove no supporting chamber face exists in a different
   generator domain. A new direct full-constraint LP pass does prove the
   sampled offset-generator domains themselves are not supporting faces of the
   exported degree-bounded Mori cone: target `7` and `8` both report
   `full_status=lp_no_solution` for their offset-generator rows. These qN
   shape matches are therefore explanatory artifacts, not promotable chamber
   semigroups.
   The origin-support certificate path now uses the same status policy:
   rerunning target `7` and `8` on the schema-4 context with the `4096` guard
   reports target `7` relation/shared/union statuses
   `skipped_single_generator_higher_codimension`, `lp_no_certificate_rank_13`,
   and `lp_no_certificate_rank_194`; target `8` reports
   `skipped_single_generator_higher_codimension`, `lp_no_certificate_rank_9`,
   and `lp_no_certificate_rank_177`. The previous target-`7` shared-facet
   optimizer `Unknown` no longer aborts the report, but it also does not
   certify a face.
   The offset-generator report now also records the exact sparse generators,
   generator degree buckets, and a raw `FIND_GV=false` GW coefficient trace for
   the provided-generator domain. Fresh target `7`/`8` runs show that the
   target `7` degree-six and degree-eight candidate domains and target `8`
   degree-six candidate domain have integral source readouts with `GV=-2`, but
   still contain `13`, `8`, and `9` non-integral lower GW candidates,
   respectively. Target `8`'s degree-eight domain remains an integral-run
   `hkty_error`; the raw trace reaches source candidate `2` and exposes `9`
   non-integral lower candidates. The vendored `cygv` error now identifies the
   first failing coefficient as
   `element_nonzero=[(203,-3),(206,1),(209,1)]`, pivot component `-3`,
   instanton coefficient `-4`, GV candidate `4/3`, rounded candidate `1`.
   The classified raw GW sample marks that curve as
   `unknown_not_toric_covered` and `source_ray_matches_missing_target`, i.e.
   the first integrality failure is the target-`8` missing ray inside the
   lower coefficient history of this candidate source domain. Full-trace
   counts classify all `9` fractional candidates as
   `unknown_not_toric_covered`; `8` are `not_source_degree_bounded_ray` and
   `1` is `source_ray_matches_missing_target`.
   These traces explain why the candidate domains cannot be promoted without a
   source-derived supporting-face/chamber certificate.
   The target-coefficient balance reconstructs the small-domain
   pre-subtraction readouts: target `7` has lower-source sum `3`,
   post-subtraction coefficient `0`, pivot component `2`, and
   pre-subtraction candidate `3/2`; target `8` has lower-source sum `-6`,
   pivot component `-3`, and pre-subtraction candidate `2`. These are not
   promoted GV values; they are evidence that the sampled domains cancel to
   zero before target qN materialization. The report also compares these
   pre-subtraction candidates with the CMS-general-divisor formula candidates:
   target `7` is `3/2` versus formula `3`, and target `8` is `2` versus
   formula `3`. So the sampled path-support domains are neither the final
   `cygv` target history nor the local formula source domain.
   This comparison is currently bounded to the
   target `7`/`8` degree-ten pair: the all-target run timed out at `900s`,
   targets `2`-`5` time out under `180s` per-target probes, and targets
   `0`/`1`/`6` reach non-integer HKTY errors in their small path-support
   domains.
   The source-ray summary now shows that shared ray has local charge signature
   `[-1,-1,-1,1,2]` and a matching CMS check. The residual-source predecessor
   queue now shows all three occurrences have shared origin-circuit witness
   relations and `source_derived_full_facet_context`. The source-derived scalar
   GV map enforces this full facet context before a source value enters
   qN-history status, and the actual-local-cygv 4-214-647 report now records
   `20` actual local cygv imports matching CMS formulas, `21` actual local cygv
   imports without CMS formula checks, `32` remaining full-facet CMS formula
   imports, `58` rows with no integral matching CMS formula, and `115` rows
   with no origin-circuit witness. Accepted imports split by local charge
   signature as `-1,-1,-1,1,1,1:28`, `-1,-1,-1,1,2:41`,
   `-1,-1,1,1:3`, and `-1,-1,-1,-1,1,3:1`; the residual predecessors stay in
   the `-1,-1,-1,1,2` family. The residual predecessor sample carries its
   local `cygv` skeleton: both residual source rows use including-origin
   `q=[-1,-2,1,1,1]`, semigroup `[[1]]`, grading `[1]`, and a unit-tensor
   probe that gives `GV=-2` against expected `-2`. Cyrus now certifies this as
   the source-derived local bundle phase `O(-1)+O(-2)->P^2`, with unit
   hyperplane tensor `κ_000=1` and a positive-base chamber certificate; the
   residual queue has no local missing-input counts and all three occurrences
   are `ready_for_actual_cygv_call`. The CMS divisor cubic suggests tensor
   value `3`, but the raw-cubic primitive probe gives `GV=-6`; the certified
   local bundle tensor gives the expected `GV=-2` and is the promoted source
   value. The local unit-phase probe now also exports actual `cygv` qN trace
   metadata for that source-derived unit tensor, with the residual-source
   regression verifying one materialized target qN polynomial for the
   `GV=-2` source row. The context export now also keeps the nonzero rational
   divisor-basis
   and ambient-basis coefficients from successful CMS-general-divisor solves,
   so the next step is consuming these computed source values inside the
   broader finite qN history rather than relying on the CMS scalar formula
   status.
   These lower classes are now classified against the exported source-ray
   context: for both degree-10 targets the sampled lookups split as four known
   toric source rays, four non-toric source rays, and eight composite classes
   that are not source degree-bounded rays. None is one of the nine direct
   missing target rows, so reproducing the target qN history requires broader
   finite-semigroup history, not just assigning values to the nine target
   curves.
5. **Potent-ray local diagnostics.** Rank-two N=10 and rank-four potent-ray
   rows are still incomplete, but they are validation diagnostics rather than
   blockers for the compact GA-ready GV engine.

The alternative flop/stable-Weyl path is also still uncertified. The context
report now aggregates the existing CMS divisor checks: all nine remaining
targets have formula-shaped divisor candidates, but the 14 exact
divisor-intersection checks all report
`cms_general_divisor_no_rational_divisor_solution`. The algebraic flop
transforms are present, but Cyrus still lacks a source-derived shrinking
divisor and certified `n_C^0` for these classes. A regenerated default
schema-3 report now records the cheap target-ray certificate as well:
`not_extremal_by_exact_integer_semigroup_decomposition=5`,
`not_extremal_by_exact_rational_cone_decomposition=4`, and
`active_support_not_certified_as_codimension_one_face=9`. The unresolved GV
classes therefore cannot be handled by directly treating the target rays as
certified chamber walls. The per-target active-support report now serializes
that same face-certificate status next to each active-support GV/HKTY result,
so the six active-support `GV=0` diagnostics remain visibly non-promotable. The
`--run-integer-diamonds` diagnostic is also now aggregated as
`target_status_counts`: three integer-semigroup diamonds compute target GV
`0`, two fail compact HKTY with a non-integral invariant, and the four
rational-cone targets are skipped. That diagnostic does not resolve the
remaining GV rows. The active decomposition generators now classify as 12
toric-covered leaves, three source-derived GV leaves, two other missing
targets, and five concrete uncovered source-ray leaves after the
first-principles export unions active non-covered dependencies into the
source-ray diagnostic set and exports uncovered source-ray toric diagnostics.
The context consumer imports those toric diagnostics as source-derived known GV
values after exact degree and duplicate-conflict checks. The report now carries
`active_decomposition_unresolved_source_leaf_sample`, which gives the exact
curve and ambient support for each unresolved dependency; the current sample
contains two missing-target links plus five matching uncovered source-ray
leaves. The previous degree-4 leaf with ambient support
`[(6,1),(200,1),(210,-2)]` is now classified from a two-face toric diagnostic
with `GV=-2`. The remaining source-ray work is therefore
local phase/intersection/chamber certification for the higher-rank
origin-circuit leaves, not finding which lower leaves are involved. The context
export now keeps the complete `origin_circuit_witnesses` list for sampled
missing/source rows, not just the first witness, so later chamber certification
can compare alternate facet-pair witnesses for multi-witness rows such as the
degree-14 `[-1,-2,1,1,1]` dependency. The report reader now builds
origin-circuit support-domain and facet-context summaries from all serialized
witnesses, including explicit mixed-status summaries for rows whose witnesses
do not agree. A fresh all-witness release export shows all five multi-witness
missing targets share the same local relation across their witnesses, so the
open problem is facet/chamber phase and local intersection tensor
certification, not choosing among competing relation vectors. The
per-witness domain diagnostic now checks the same point more directly: all 14
witness relation domains and all 14 shared-facet domains fail the existing
LP/exact supporting-face certificate, and the two degree-10 facet-union domains
small enough for the `512` generator guard fail as well. Larger facet unions
remain explicitly skipped. The report now has an opt-in rational span-closure
scan for those witness domains, backed by a per-domain exact rational row-echelon
basis rather than repeated full rank recomputation. On the current schema-3
context, all nine targets report relation/shared/union domains as
`span_closed_under_degree_bounded_context` with zero extra same-span generators,
after scanning `720..2963` unique degree-bounded candidates depending on target
degree. This rules out the cheap "small domain forgot forced span rows"
explanation; the blocker remains supporting-face or chamber-semigroup
certification. These witness domains also still contain unresolved
source data: the shared-facet domains aggregate to `228` source rays without a
derived GV value, `20` other missing-target hits, and `10` uncovered-source-ray
hits, despite also containing many toric-covered and source-derived-GV
generators. Across all per-witness relation/shared/facet-union domains, the
new unresolved queue contains `2017` unique non-known generators over `9611`
domain occurrences, dominated by broad facet-union source rays. After exporting
the full degree-bounded toric GV diagnostic context (`1053` rows) from the
first-principles runner, the shared-facet unknown bucket drops from `228` to
`36` occurrences and the all-domain unresolved queue drops to `1631` unique
generators over `7360` occurrences. The shared-facet-only unresolved queue is
now explicit as well: `33` unique non-known generators over `66` occurrences,
made up of the nine missing targets, `4` uncovered-source-ray hits, and `20`
source degree-bounded rays still missing toric/source-derived GV values. A
fresh export now includes complete first-principles stats for those `20` source
rays: all are origin-circuit rows, all are blocked on local q-matrix
phase/intersection tensor/chamber certification, and the CMS checks expose `4`
integral inferred-degree matches with formula value `3`, `2` non-integral
matches, and `22` no-solution checks. The local unit-tensor probe matches only
`5/20` expected formula sets (`GV=-2` cases); formula-`3` rows give unit
candidate `6` and origin-omitted candidate `-9`, so the tensor/chamber model
cannot be guessed from the unit probe. The CMS raw-divisor-cubic primitive
probe also fails for the four integral formula-`3` rows, producing `-36`,
`-42`, `-42`, and `-48`, so the missing tensor is not the raw divisor cubic
either. The
enriched active-leaf report
now also carries the matching source rays' local q rows and CMS check counts:
two origin-circuit lower leaves are primitive `neg2` q-row candidates, while
three have three negative primitive intersections and still need a different
certified phase or chamber model. The unit-tensor phase probe matches only the
degree-14 lower leaf (`GV=-2`); the degree-12 leaf gives `-3` versus expected
`3`, and the three degree-10 leaves give `0` origin-included or `-1`
origin-omitted against expected `{-2,1}`. These are diagnostics, not promoted
GV values. The same active-leaf report now retains non-promotable CMS rational
solutions: two integral solutions mismatch the inferred normal degree, one is
nonintegral, and five leaves have no CMS solution summary. Their candidate
statuses are marked `diagnostic_from_*_not_promoted`.
An opt-in target-level integer tensor scan now checks whether a bounded scalar
normalization of the one-parameter local tensor could explain the nine missing
targets. With scan bound `8`, four targets match their CMS formula candidate
at tensor value `1`, but five do not match any integer tensor in the scan
range. Those four matches are still uncertified, so this narrows the work but
does not close the corrected-chamber GV blocker.
The bounded-decomposition diamond diagnostic now exposes raw GW coefficient
history for explicit semigroups even when the integral compact HKTY run fails.
On
`/tmp/cyrus_gv_context_target8_explicit_diamond_gw_report.json`, target `8`'s
parent-only diamond still has `diamond_status=hkty_error`, but the raw trace
has seven coefficient entries and exactly one nonintegral candidate. That
candidate is classified as
`known_qn_history_status=unknown_not_toric_covered` and
`source_class_status=source_ray_matches_missing_target`; the target coefficient
status is `nonzero_gw` with candidate `4/3`. This is evidence against a hidden
lower-leaf explanation for the target-`8` failure, but it remains a diagnostic
because the semigroup/chamber certificate is still absent.
The same diagnostic on target `7` preserves the earlier conclusion: both
parent-only bounded diamonds compute `GV=0`, and their target-row raw GW
statuses are `missing_instanton_coefficient` and `zero_or_absent_gw`. One
diamond has a fractional lower degree-4 raw candidate, but it is
`not_source_degree_bounded_ray`, not the target row. Target `7` therefore still
needs the correct compact `q_N` history; the bounded diamonds do not provide a
promotable nonzero target coefficient.
Target-level reports now embed the local unit probe and integer tensor scan,
so the local-circuit scalar check can be compared directly with compact
bounded-diamond traces. For target `8`, the local oriented one-parameter model
`[-1,2,-3,1,1]` gives unit-tensor candidate `3`, matching the CMS formula sum
at tensor value `1`, while the compact bounded diamond's target raw GW
candidate is `4/3`. This keeps the q-matrix orientation in the plausible
bucket and moves the unresolved discrepancy to certified chamber/intersection
tensor history.
The latest target-family classifier now reports the four compact one-parameter
rows as
`uncertified_one_parameter_split_bundle_over_weighted_p2:base=1,1,2;bundle=1,3;base_hyperplane_square=1/2`.
Equivalently, their charge data has the shape of an
`O(-1)+O(-3)->P(1,1,2)` local model, not the already-certified
`O(-1)+O(-2)->P2` source family. Certified
`source_derived_local_p2_bundle_family_with_tensor_chamber_certificate` rows are
counted separately, so the local unit-tensor match is explicitly not a
promotable replacement for the missing compact GV target correction. The same
report now gives those four rows the tensor/chamber blocker
`weighted_p2_split_bundle_requires_source_derived_resolution_chamber`, while
the five non-CY no-origin diagnostic rows remain blocked on the no-origin phase
choice.
The latest context report also records a source-witness resolution hint for
this exact weighted family: all four weighted rows have one zero-relation point
in the shared two-simplex of the origin-circuit witness
(`target 3 -> 202`, `target 6 -> 199`, `target 7 -> 55`,
`target 8 -> 55`). The aggregate status is
`source_witness_weighted_p2_split_bundle_has_single_zero_relation_shared_resolution_ray=4`.
This is not yet a chamber certificate, but it identifies the source-derived
extra ray that the resolved local phase/intersection construction must use.
Fresh `mcallister_first_principles` context dumps now serialize
`shared_two_simplex_points`, so that extra ray can be recovered with
coordinates instead of only by point ID. The context consumer now classifies
the resolved-support attempt separately; old dumps are blocked at
`weighted_p2_resolved_shared_support_missing_zero_shared_ray_coordinates`, and
the declared point coordinates show the naive support with the extra ray has
affine rank `4`, not a compact threefold phase. The chamber task is therefore
a projection/chamber-map problem, not an append-the-ray cygv handoff. On the
fresh shared-point context, the projection diagnostic classifies all four
weighted rows as
`weighted_p2_zero_shared_ray_has_primitive_unit_height_off_relation_hyperplane_requires_projection_or_chamber_map`;
the zero-shared ray has signed height `-1` above the relation hyperplane in
each case.
The full-union compatibility report now separates restricted-support evidence
from the actual full chamber:
`/tmp/cyrus_all_full_union_compatibility_report.json` has five
`resolved_full_union_compatibility_skipped:resolved_shared_chamber_global_height_not_weighted_p2_split_bundle`
rows and four
`resolved_global_height_strict_but_full_union_selects_star_extras_requires_chamber_transport`
rows. For targets `7`/`8`, the corrected global height is strict on the
exclusive-pair resolved support, but the full star-union global regular
triangulation selects the serialized star extras. The next object is therefore
chamber transport/full chamber-map history, not promotion of a scalar local GV
value or a direct compact `cygv` domain on the restricted support.
The fresh star export then checks the actual chamber neighbors of those shared
two-simplices: all four weighted rows are
`weighted_p2_zero_shared_star_uses_two_alternate_chamber_points`. Targets `3`
and `6` use alternate extras `[2]`/`[46]`, while targets `7` and `8` use
`[195]`/`[212]`. That sharpens the blocker to constructing the resolved source
projection/chamber map from the actual triangulation star, not from the
exclusive pair that appears in the origin-circuit witness relation.
The next regenerated report includes those star-extra coordinates and computes
the actual six-point star-support charge rows. All four weighted rows have
affine rank `4` with one zero-sum charge row: `[2,-1,-1,0,-2,2]` for targets
`3`/`6` and `[1,0,1,-1,0,-1]` for targets `7`/`8`. This gives the next
first-principles candidate to certify against local chamber/intersection data,
instead of using the old scalar weighted-family formula.
The zero-charge reduction classifier sharpens the split: targets `3`/`6`
become the sign-flipped weighted `O(-2)+O(-2)->P(1,1,2)` family, while targets
`7`/`8` become resolved-conifold rows with spectator zero columns. This is
still not a promoted GV value; it identifies which local phase certificates
must be derived next.
The uncertified reduced-row unit-tensor probe does not close the gap: targets
`3`/`6` compute toy `GV=0`, while targets `7`/`8` are rejected by compact
`cygv` as too low-dimensional. The next step is therefore certificate work, not
using the reduced rows as compact GV inputs.
The target-vs-star support guardrail now shows why: each of the four weighted
actual star supports is missing at least two points from the target
origin-circuit relation. The star rows are neighboring chamber rows, not the
target rows whose GV values are missing.
The union-support comparison narrows the next chamber-history task further:
targets `7`/`8` are integral in the union relation lattice with the star row,
whereas targets `3`/`6` have non-integral target coordinates in that lattice.
That split is now serialized instead of inferred from point IDs. The exact
coordinates in `/tmp/cyrus_gv_context_star_union_rational_report.json` are
`target 3 = [1/2,1/2,-3/2]`, `target 6 = [1/2,-2,3/2]`,
`target 7 = [3,0,-1]`, and `target 8 = [1,0,-1]`; the corresponding actual
star rows are integral (`[0,0,1]`, `[0,2,-1]`, `[0,1,0]`, `[0,1,0]`).
The CYTools no-origin global-basis projection is now explicit: using the
`curve_basis(include_origin=False)` convention, all four weighted target/star
rows and target±star combinations project integrally. The projected target
rows match the missing target basis rows exactly; for example target `7`
projects to `[(44,2),(203,1),(206,-3),(209,1)]` and target `8` to
`[(203,-3),(206,1),(209,1)]`.
The follow-up lookup against the known toric-covered and source-derived GV
maps keeps this result negative: all projected target/star/target±star roles
for targets `3`, `6`, `7`, and `8` are `unknown_not_toric_covered`. Their
degrees are target `26/12/10/10`, star `-36/-24/-4/-4`,
target-minus-star `62/36/14/14`, and target-plus-star `-10/-12/6/6`. The
negative star-side degrees confirm that the star-union object is chamber
history rather than an effective scalar GV row already present in the known
maps.
The sign-reversed lookup in
`/tmp/cyrus_gv_context_star_union_opposite_report.json` rules out the simplest
adjacent-effective-class shortcut for targets `3`/`6`: `-star` and
`-(target+star)` remain `unknown_not_toric_covered`. For targets `7`/`8`,
`-star` is a known source-derived degree-4 class with GV `1`, but the target
and target-plus-star rows remain unknown. This makes the star row useful
adjacent-source evidence, not a replacement for the missing target GV history.
The provenance-enriched report also shows the source relation for that known
target `7`/`8` opposite-star row:
`[(0,-1),(55,-1),(195,1),(212,1)]`, the resolved-conifold origin circuit. The
target `6` opposite-star row is at least a degree-bounded source ray, but has no
known toric/source-derived GV value; target `3`'s opposite star is not in the
degree-bounded source-ray context.
The bounded lower-seed diagnostic in
`/tmp/cyrus_gv_context_star_union_seed_decomposition_report.json` keeps that
resolved-conifold side from being a new blocker: the target `7`/`8`
opposite-star row decomposes into two known degree-2 toric/source rows. The
remaining degree-6 `target+star` row is not a degree-bounded source ray and has
no decomposition into up to four lower-degree source seeds in the current
finite context. So the unresolved compact-history object is the degree-6
`target+star` chamber row, not the known `-star` side.
The coordinate report now serializes that row directly in the union charge
basis: target `7` has `target+star=[3,1,-1]`, target `8` has
`target+star=[1,1,-1]`, and both point-row signatures sort to
`[-3,-1,-1,1,1,1,2]` with zero total charge.
The target-plus-star unit-tensor diagnostic report runs the explicit row
through the one-parameter `cygv` path with tensor value `1`; targets `7` and
`8` both return toy GV `0`, but the status now records
`not_compact_threefold_shape_cy_dim_5`. This is evidence against a naive
normalization shortcut, not a promoted GV source.
A synthetic target-8 `target+star` path-history check makes the scale sharper:
at degree `6`, the broad cygv-style source closure has `463` seeds, `355`
reduced seeds, and `65538` first-layer elements, but the `target+star` row is
not yet in closure and has zero predecessor differences in the previous-degree
window. A second-layer run exceeded a `180s` guard. The remaining object is
therefore not a cheap one-layer closure or scalar normalization issue.
The height-profile report now tags every union-support point by role and
coefficient. Targets `3`/`6` have only a coefficient-zero zero-shared ray off
the relation hyperplane, while targets `7`/`8` have a height `-1` star-side
pair: the zero-shared ray with coefficient `+1` and star extra `212` with
coefficient `-1`. The remaining certificate work is therefore a real
chamber-map/intersection problem, not just a missing projection or lookup.
The off-height projection lookup confirms the pair is not independently
promotable: targets `3`/`6` have zero off-height components, and the target
`7`/`8` off-height star component `[(55,1),(212,-1)]` has no integral global
curve-basis projection. The wall-crossing object must therefore be certified
locally before it can affect compact GV history.
Origin-circuit witness-domain compact `cygv` probes are now opt-in and guarded
by `--origin-witness-cygv-generator-limit`. The traced target `8` report shows
the single-generator relation domain computes `GV=2` with one materialized
target qN term and raw target GW candidate `2`, while the 14-generator
shared-facet domain computes `GV=0` with `24` qN polynomials, no target qN
polynomial, and raw target GW candidate `0`. The traced target `7` report
shows the relation domain fails compact HKTY integrality on the missing target
itself with raw candidate `3/2`, while the 21-generator shared-facet domain
computes `GV=0` with `51` qN polynomials, no target qN polynomial, and raw
target GW candidate `0`. The 458/506-generator facet unions are skipped with
the default diagnostic guard. The new certificate diagnostic fields show the
same LP failure modes explicitly: target `8` shared facets have full/aggregate
`lp_no_solution` and anchors split `14` cutting-round limits to `2`
no-solution outcomes; target `7` shared facets have full LP `Unknown`,
aggregate `lp_no_solution`, and anchors split `9` cutting-round limits to `7`
no-solution outcomes; both facet unions have full/aggregate/anchor
`lp_no_solution`. These are not promoted because the LP/exact witness-domain
face certificates still fail, but they rule out the small relation/shared-facet
domains as the corrected McAllister compact histories. Origin-witness
certificate diagnostics now also honor the `--supporting-face-lp-*` CLI
limits; a target `8` smoke run with `2` anchors and `1` cutting round reports
exactly two anchor attempts in the shared and union domains. Higher-budget
certificate-only target `7`/`8` probes with `64` anchors and `256` cutting
rounds still do not promote a witness domain: target `8` shared/union domains
are all `lp_no_solution`, target `7` union domains are all `lp_no_solution`,
and target `7` shared domains retain full LP `Unknown` but have aggregate and
all `64` anchors at `lp_no_solution`.

Cyrus now has a focused native secondary-cone chamber primitive in
`crates/cyrus-core/src/triangulation/secondary.rs`. It ports the CYTools
adjacent-simplex circuit step for full-dimensional triangulations, emits
integer hyperplane normals in point-index order, and evaluates typed finite
height-vector pairings with a strict positive interior test. The focused test
`cargo test -p cyrus-core triangulation::secondary -- --nocapture` covers a
square diagonal circuit, wall/outside height rejection, dimension mismatch, and
non-full-dimensional simplex rejection. This is not a corrected-chamber GV
certificate yet; it is the reusable chamber-membership primitive needed before
local chamber/intersection data can be promoted.
`mcallister_first_principles --dump-corrected-chamber-gv-context` now exports a
top-level `secondary_cone_height_certificate` for the corrected chamber,
including CYTools' `1e-6` wall tolerance, hyperplane count, pairing min/max,
and strict-interior status. The export test
`corrected_chamber_context_export_includes_uncovered_source_toric_diagnostics`
asserts the field is serialized. `mcallister_gv_context` now deserializes,
validates, and reports that certificate, including rejecting inconsistent
hyperplane/pairing counts. This still certifies only the global exported
chamber; candidate local phases remain unpromoted until they carry their own
source-derived chamber/intersection certificates.
The local `cygv` source-resolution summaries now carry both the global
secondary-cone certificate status and a separate local phase chamber-membership
status. A local source hint can therefore say "global corrected chamber is
strictly inside" without pretending the candidate local phase has been
certified; promotion is still blocked until the local phase supplies its own
phase `q` matrix and chamber certificate. The report also aggregates those
local chamber-membership statuses, so this blocker is visible as a top-level
count rather than only in sampled target rows.
With the validation-only global certificate present, the available 4-214-647
contexts split those local blockers into five no-origin/non-CY phase rows and
four weighted-`P2` rows. If the context does not carry triangulation-star
samples, the weighted rows stop at an explicit missing-star-simplices blocker;
when star samples are present, the same status path follows the actual
star-union relation and records whether the target/star relation is integral
in the union chamber lattice or fails at rational-coordinate integrality.
The report now also aggregates star-union off-height component statuses by
role. This makes the wall-crossing obstruction visible at top level: zero
off-height components are harmless, while nonzero components that are not
integral global curve-basis classes remain local chamber-map data rather than
promotable compact GV history.
The same report now aggregates the star-union global-basis projection and
known-qN lookup statuses. On a current release no-replay 4-214-647 export, the
corrected chamber certificate is strictly inside the secondary cone
(`1258` hyperplanes, minimum pairing `0.0011511487808064658`) and every missing
target has two triangulation-star samples. The local blockers split into five
non-`P2`/origin-omitting phase rows, two weighted-`P2` rows with integral
target/star union coordinates, and two weighted-`P2` rows whose target relation
has non-integral local charge coordinates. All four weighted rows do project
integrally to the global no-origin CYTools basis for target, star, and the
target±star sums, but the positive-direction target/star/sum qN lookups remain
`unknown_not_toric_covered`; two negative star opposites have known
source-derived nonzero GV values. The run still ends at
`V_string=4711.504666533563`, with a corrected-volume residual of about
`0.072167298`, so these counts are blocker evidence, not a reproduced result.
The weighted-row chamber coverage diagnostic is stricter: all four weighted
rows are now blocked because the two star simplices do not cover the full
target/star union support. In each case the uncovered points are precisely the
origin-circuit exclusive target-relation points (`[208,211]`, `[49,52]`,
`[208,214]`, and `[211,214]` for targets `3`, `6`, `7`, and `8`). Thus
integral global-basis projection is not enough to promote the local value; a
source-derived chamber covering those target-relation points, or a different
certified support semigroup, is still required.
A local secondary-cone diagnostic for the exclusive-pair resolved support
(origin + shared two-simplex + the two target-relation exclusive points) also
does not supply that chamber. For all four weighted rows the candidate has one
secondary hyperplane, but the affine-height pairing is exactly `0`, so the
exclusive-pair support lies on the wall rather than strictly inside either the
candidate chamber or its sign-flipped neighbor.
The local phase-membership summary now uses that sharper certificate before
leaving the blocker at the coarser integral star-union relation. Regenerating
target `7` and target `8` reports as
`/tmp/cyrus_target7_refined_local_phase_report.json` and
`/tmp/cyrus_target8_refined_local_phase_report.json` gives the top-level local
phase bucket
`local_phase_chamber_blocked_weighted_p2_resolved_shared_affine_wall_global_height_strict:weighted_p2_resolved_shared_chamber_outside_or_on_wall`
for both targets. The same rows still record
`star_union_global_regular_shared_face_matches_serialized_star_extras`, so the
corrected global height selects the serialized star extras in the full
star-union support, while its restriction to the exclusive-pair resolved
support is strictly inside that local secondary cone. This is chamber-map
evidence, but not a promoted local semigroup or intersection tensor.
The all-target refined report
`/tmp/cyrus_all_phase_certificate_first_report.json` now splits the nine
missing rows as five non-`P2`/origin-omitting phase blockers and four
weighted-`P2` rows blocked by the resolved-shared chamber wall. The separate
star-union relation buckets still preserve the lower-level distinction: two
of those weighted rows have integral target/star union coordinates and two have
nonintegral target coordinates in the star-union charge basis.
The local `cygv` readiness report now has a separate promotion gate that
requires this phase certificate in addition to the local shape/tensor/semigroup
inputs. Regenerating `/tmp/cyrus_all_promotion_gated_report.json` keeps the
shape-level actual-call readiness at `blocked_missing_source_derived_inputs:9`,
but the promotion-readiness split is now five
`blocked_local_phase_chamber_certificate:local_phase_chamber_blocked_local_chamber_certificate_blocked_not_including_origin_p2_bundle_phase`
rows and four
`blocked_local_phase_chamber_certificate:local_phase_chamber_blocked_weighted_p2_resolved_shared_affine_wall_global_height_strict_weighted_p2_resolved_shared_chamber_outside_or_on_wall`
rows. The promotion missing-input aggregate adds
`local_phase_chamber_membership_certificate:9`, so no local support can be
mistaken for a promotable `cygv` domain while the chamber certificate is
absent.
The resolved-support chamber diagnostic now also evaluates the actual
corrected global height on the restricted exclusive-pair support. In
`/tmp/cyrus_all_resolved_global_height_report.json`, the five non-`P2` rows are
`resolved_shared_chamber_global_height_not_weighted_p2_split_bundle`, while all
four weighted rows are
`weighted_p2_resolved_shared_chamber_global_height_strictly_inside_exclusive_pair_secondary_cone`.
The affine projection certificate for those weighted rows remains
`weighted_p2_resolved_shared_chamber_outside_or_on_wall`, so the missing object
is the full chamber map/transport that reconciles the strict restricted support
with the full star-union chamber, not a scalar local GV value.
The full-union compatibility field now records this as a top-level aggregate:
five rows are skipped as non-weighted-P2 split bundles, and four rows are
`resolved_global_height_strict_but_full_union_selects_star_extras_requires_chamber_transport`.
Targets `7` and `8` both pair that strict restricted-support status with
`star_union_global_regular_shared_face_matches_serialized_star_extras`, making
the chamber-transport gap explicit.
A bounded supporting-face rerun for the two integral weighted rows (`target 7`
and `target 8`) does not provide that missing certificate. With
`--certify-origin-witness-domains`, two LP anchors, and the facet-union cap
raised to `1024`, relation supports are still skipped as single-generator
codimension-`213` domains, shared-facet domains return `lp_no_certificate`,
and the full facet-union domains also return `lp_no_certificate` at ranks
`194` and `177` respectively. The cheap face-certification path is therefore
closed for these rows.

A no-replay finite-semigroup pruning run also closes the hypothesis that the
live residual is primarily a pair-pruning artifact. The command

```text
/usr/bin/time -p timeout 360 ./target/release/mcallister_first_principles \
  --data-dir /Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 \
  --small-curve-pruning finite-semigroup \
  --diagnose-corrected-chamber-gv \
  --out /tmp/cyrus_finite_semigroup_pruning_summary.json
```

completed in `real 115.77`. It pruned the input-chamber `419` subcutoff
curves to `339`, the checkpoint-t corrected chamber `556` subcutoff curves to
`400` with `toric_missing=0`, and the solved-t corrected chamber `561`
subcutoff curves to `399` with `toric_missing=0`. The final
`V_string=4711.504343103075` still differs from the corrected-volume
checkpoint by about `0.071843868`. Thus finite selected-set pruning removes
the visible toric-missing rows but does not reproduce the corrected-target GV
vector or the final volume.
The corrected-chamber trace writer now records `small_curve_pruning`, generic
`pruned_*` counts, and each selected toric curve's sparse KKLT target
contribution vector. Regenerating
`/tmp/cyrus_finite_semigroup_pruning_trace_with_contrib.json` with
`CYRUS_CORRECTED_CHAMBER_GV_TRACE_JSON` produced a `2.7M` finite-semigroup
trace over `400` pruned toric-covered curves. This makes the next comparison
curve-local: test or replace the large-contribution toric rows and their
compact qN histories directly, instead of fitting final aggregate weights.
Because the full corrected-chamber GV context export still exceeded a `300s`
debug timeout before writing JSON, the runner also has a narrow
`--dump-corrected-chamber-secondary-certificate` flag. On the validation-only
`--allow-downstream-kahler` checkpoint path for 4-214-647, it writes in about
`10.32s` and reports `1262` hyperplanes with minimum pairing
`0.0002917945770377628`, maximum pairing `87.93218664229505`, and
`strictly_inside=true`. This remains a chamber-certificate validation check,
not a production no-replay completion proof.

## Next Concrete Action

The next implementation should be one of these, in order:

1. Extend the secondary-cone certificate from the global corrected chamber to
   each candidate local phase, so a local intersection tensor or semigroup
   cannot be promoted without an explicit chamber-membership certificate.
2. For corrected-chamber missing GV classes, construct a source-derived compact
   semigroup or certified supporting-face semigroup and hand it to the existing
   `cygv` wrappers. If the semigroup cannot be certified, keep the result
   diagnostic-only.
3. Clean up the remaining matrix-basis diagnostic/export edges in
   `mcallister_first_principles`: chamber diagnostics currently operate after
   transforming production Kähler coordinates back to the computed CYTools
   index basis, which is correct for geometry but not a fully matrix-native
   diagnostic representation.
4. Keep local CKYZ/potent-ray checks bounded as source validation. Do not raise
   the large local N gates or add more CKYZ performance machinery unless a
   compact `cygv` call cannot express a required, certified input.

Do not mark the objective complete until every checklist row above is either
implemented with direct evidence or explicitly removed from the objective.
