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
| CYTools/cygv compact GV wrapper contract is understood and mirrored | `docs/CYGV_AUDIT.md`; direct source read of `reference/cytools/src/cytools/calabiyau.py::_compute_gvs_gws`; downloaded `cygv-0.1.2` source | Implemented for the generic compact handoff boundary. CYTools builds inputs, then delegates to `cygv.compute_gv`; no hidden Python GV algorithm remains. |
| High-dimensional small toric curve selection is computed, not read | `docs/GA_READY_PIPELINE_AUDIT.md`; stage 5 small-curve tests | Implemented for the checkpoint rule. Pair pruning matches `small_curves.dat`; finite-semigroup pruning would remove five additional curves, so both policies remain explicit. |
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
   toric degree-2 differences. The unfiltered bounded report now promotes
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
   supports. The LP diagnostic now treats HiGHS `NoSolutionFound` as
   no-certificate rather than a hard error, and raising the origin-support guard
   to `4096` for the two degree-10 targets checks their larger source domains
   without promotion: target `7` has relation/shared/union ranks `1/13/194`,
   target `8` has ranks `1/9/177`, and all checked domains return
   `origin_support_lp_no_certificate_*`.
   The context tool now also has an opt-in exact extremal-ray probe. For the
   current schema-3 context it verifies that all nine remaining targets are
   non-extremal in the finite degree-bounded cone by exact decomposition
   certificates: five integer-semigroup decompositions and four rational-cone
   decompositions. This is not a substitute for the missing compact/chamber
   semigroup, because GV values are not determined by finite-cone extremality
   alone. The path-history probe now aggregates closure/lower-seed status
   counts; on the two degree-10 targets, bounded closure hits the `10000`
   element cap and the small lower-seed diamonds either give `GV=0` or
   non-integral cygv output. The generation counters show the first closure
   layer alone would add `130414` elements from `720` seeds, so the first
   visible decomposition is not a valid replacement for the full HKTY history.
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
   The all-target origin-omitted shape aggregate splits the nine unresolved
   origin circuits into five compact-threefold-shaped no-origin rows and four
   CY-dimension-2 rows, so the remaining phase work is now split rather than a
   single generic bucket. A target-level unit-tensor phase probe then shows
   that the four CY-dimension-2 rows match their formula candidate set only in
   the origin-included convention, while only one of the five compact
   no-origin rows matches after origin omission. The other four compact
   no-origin rows still mismatch, so a single origin-column rule is not enough.
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
   source for it. The same small path-support `cygv` domain now reports
   per-predecessor lookups: it matches the four known degree-two toric values
   and assigns six nonzero plus six zero/absent values among the unknown
   non-toric lower classes, but the target remains `GV=0`. These lower-class
   values are therefore diagnostic artifacts of the small domain, not
   promotable compact GV history.
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
   The source-ray summary now shows that shared ray has local charge signature
   `[-1,-1,-1,1,2]` and a matching CMS check. The residual-source predecessor
   queue now shows all three occurrences have shared origin-circuit witness
   relations and `source_derived_full_facet_context`, so the actual remaining
   inputs are `local_intersection_tensor` and `local_chamber_certificate`. The
   source-derived scalar GV map now enforces this full facet context before a
   CMS-derived value enters qN-history status, and the guarded 4-214-647 report
   still records the same two source predecessors across three occurrences. It
   also records source-ray import counts: `52` full-facet CMS formula imports,
   `79` rows with no integral matching CMS formula, and `115` rows with no
   origin-circuit witness. Those 52 imports now split by local charge signature
   as `-1,-1,-1,1,1,1:28`, `-1,-1,-1,1,2:20`, `-1,-1,1,1:3`, and
   `-1,-1,-1,-1,1,3:1`; the residual predecessors stay in the
   `-1,-1,-1,1,2` family. The residual predecessor sample now carries its
   local `cygv` skeleton: both residual source rows use including-origin
   `q=[-1,-2,1,1,1]`, semigroup `[[1]]`, grading `[1]`, and a unit-tensor
   probe that gives `GV=-2` against expected `-2`. Cyrus now certifies this as
   the source-derived local bundle phase `O(-1)+O(-2)->P^2`, with unit
   hyperplane tensor `κ_000=1` and a positive-base chamber certificate; the
   residual queue now has no local missing-input counts and all three
   occurrences are `ready_for_actual_cygv_call`. The
   CMS divisor cubic suggests tensor value `3`,
   but the raw-cubic primitive probe gives `GV=-6`; the certified local bundle
   tensor gives the expected `GV=-2`. The context export now also keeps the
   nonzero rational divisor-basis and ambient-basis coefficients from
   successful CMS-general-divisor solves, so the next step is to promote or
   consume the ready actual local-`cygv` source calculation rather than relying
   on the CMS scalar formula status.
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
remain explicitly skipped. These witness domains also still contain unresolved
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

## Next Concrete Action

The next implementation should be one of these, in order:

1. For corrected-chamber missing GV classes, construct a source-derived compact
   semigroup or certified supporting-face semigroup and hand it to the existing
   `cygv` wrappers. If the semigroup cannot be certified, keep the result
   diagnostic-only.
2. Clean up the remaining matrix-basis diagnostic/export edges in
   `mcallister_first_principles`: chamber diagnostics currently operate after
   transforming production Kähler coordinates back to the computed CYTools
   index basis, which is correct for geometry but not a fully matrix-native
   diagnostic representation.
3. Keep local CKYZ/potent-ray checks bounded as source validation. Do not raise
   the large local N gates or add more CKYZ performance machinery unless a
   compact `cygv` call cannot express a required, certified input.

Do not mark the objective complete until every checklist row above is either
implemented with direct evidence or explicitly removed from the objective.
