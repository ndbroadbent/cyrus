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
| CYTools/cygv generic GV wrapper contract | `compute_gv_invariants` now uses CYTools-style `min_points=100*h11` lattice augmentation even when `max_deg` is supplied; bounded `max_deg` lattice enumeration is isolated in `compute_gv_invariants_with_degree_bounded_lattice`; compact GV execution still runs the upstream Rust `cygv` HKTY modules, but Cyrus now constructs the `cygv::Semigroup`, maps semigroup/intersection/fundamental-period/instanton/series-inversion failures into `Result` errors, and catches remaining upstream cygv panics in unwind builds instead of relying on `cygv::compute_gv_rat_threefold` unwraps; `compute_gv_invariants_with_provided_generators_and_nef_partition` and `compute_gv_invariants_with_explicit_semigroup_and_nef_partition` expose the upstream `cygv` nef-partition boundary for source-derived complete-intersection diagnostics; `compute_ambient_intersections_cytools` ports the generic CYTools ambient toric top-form system and returns prime plus canonical-divisor intersections; `compute_complete_intersection_cy3_intersection_numbers` ports the CYTools complete-intersection reduction from ambient top intersections and nef-partition divisor classes to a CY3 triple-intersection tensor; `compute_complete_intersection_cy3_from_ambient_cytools` composes those steps while still requiring a source-derived nef partition; `cyrus_direct_cygv_chain_matches_upstream_quintic_wrapper` verifies that this direct module call chain matches `cygv::compute_gv_rat_threefold` on the quintic degree-one `2875` case; CYTools-style matrix-basis projection helpers now match `mori_cap_matrix.dot(basis.T)`; matrix divisor-basis to dual curve-basis construction is available through `compute_curve_basis_matrix_from_divisor_basis_matrix`; `DivisorBasis` and `compute_curve_basis_matrix_for_divisor_basis` make vector-vs-matrix curve-basis dispatch explicit; `divisor_basis_glsm_coordinate_matrix` and `divisor_basis_change_matrix` centralize index/matrix basis-coordinate transforms; `project_mori_cone_cap_rays_for_divisor_basis` and `gv_divisor_basis_data` bundle Mori projection, curve-basis construction, and no-origin q-matrix construction for either basis shape; `intersection_in_divisor_basis` now dispatches vector filtering and matrix-basis dense tensor pullback for the in-basis intersection tensor passed to `cygv`; `matrix_basis_quintic_handoff_runs_actual_cygv_degree_one` verifies that matrix-basis-constructed q data reaches the actual Rust `cygv` call and reproduces the quintic degree-one GV `2875`; `mcallister_first_principles`, `mcallister_gv`, and `mcallister_racetrack` vector-basis GV paths now use that bundled handoff instead of independently projecting Mori rays and hand-building q matrices; `curve_basis_matrix_without_origin_i64` and `curve_basis_q_matrix_for_divisor_basis_i64` centralize the `curve_basis(include_origin=False, as_matrix=True)` q-matrix boundary used by cygv call sites; `mcallister_first_principles --dual-basis` now accepts index or matrix source-coordinate bases for K/M flux transforms, and `--production-dual-basis` carries an index or matrix internal dual divisor basis through flat-direction intersections and compact GV input construction; core KKLT now exposes `compute_kklt_divisor_volumes_in_divisor_basis`, `compute_kklt_jacobian_in_divisor_basis`, `solve_divisor_basis_path_following`, `solve_divisor_basis_path_following_branch_candidates`, `generate_scaled_divisor_basis_branch_initializations`, and `scale_divisor_basis_kklt_branch_initialization_to_target` for matrix Kähler-coordinate bases; `mcallister_first_principles --production-primal-basis` now carries an index or matrix internal primal divisor basis through KKLT initialization scaling, path following, branch search, and branch GV-coverage ranking | Core GV input construction now supports vector or matrix bases for Mori rays, q matrices, curve-basis matrices, basis-coordinate transforms, in-basis intersections, generic ambient top-form intersections, composed ambient-to-CICY CY3 tensor reduction, and small actual-cygv compact regressions for matrix-basis handoff and direct-wrapper parity. McAllister runner now accepts matrix source bases for K/M flux-coordinate validation, can execute the dual flat/GV handoff in a matrix production basis, and can solve the primal KKLT branch path in an index or matrix production basis. Corrected-chamber diagnostics still transform production Kähler coordinates back to the computed CYTools index basis before selecting curves/exporting traces, so diagnostic representation cleanup remains |
| Mirror-side racetrack GV data | `mcallister_gv` computes compact mirror-side GV invariants with the actual Rust `cygv` min-points handoff, maps basis curves back to ambient classes, and validates every `dual_curves.dat` / `dual_curves_gv.dat` checkpoint row; `stage5_mirror_gv_checkpoint_matches_cygv_min_points` preserves this as an opt-in heavy regression | Implemented for 4-214-647: `--min-points 20000` produces 10556 ambient invariants and matches all 5177 checkpoint rows |
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
   The weighted `P(1,1,2)` target-family report now also records one
   source-witness zero-relation shared ray for each of the four compact
   weighted rows (`202`, `199`, `55`, `55`). That points the next chamber task
   at the resolved local phase containing the extra shared ray, rather than at
   another scalar tensor scan. Fresh context dumps now serialize those shared
   two-simplex point records with coordinates. A direct resolved-support check
   shows that appending the shared ray raises the affine rank to `4` for the
   weighted rows, so the missing input is a projection/chamber map rather than
   a naive six-column compact `cygv` handoff. The fresh projection diagnostic
   now also records the primitive affine hyperplane of the five-point relation:
   targets `3`, `7`, and `8` use `[0,1,0,-1,0]`, target `6` uses
   `[0,0,1,0,0]`, and each zero-shared ray has signed height `-1`. Thus the
   extra ray is unit-height resolution data off the local CY hyperplane; the
   unresolved step is choosing/certifying the source chamber projection, not
   discovering another missing point.
   A follow-up star export shows the actual shared two-simplex chamber star
   uses alternate neighbor pairs in every weighted case:
   targets `3` and `6` use `[2]`/`[46]`, and targets `7` and `8` use
   `[195]`/`[212]`. The origin-circuit exclusive pair is therefore diagnostic
   relation data, not the chamber-star pair that should be fed into a resolved
   source projection.
   The star-support export now carries coordinates for those extras and
   computes actual six-point charge rows from the triangulation star. The four
   weighted rows split into two shared rows, `[2,-1,-1,0,-2,2]` for targets
   `3`/`6` and `[1,0,1,-1,0,-1]` for targets `7`/`8`, each with affine rank
   `4` and zero total charge. These rows are candidate local chamber data, not
   promoted GV inputs until the chamber/intersection certificate is derived.
   A zero-charge reduction pass classifies targets `3`/`6` as a sign-flipped
   weighted `O(-2)+O(-2)->P(1,1,2)` family and targets `7`/`8` as
   resolved-conifold charge rows with spectator zero columns. This narrows the
   next work to certifying the local phase maps for those reduced families.
   A toy unit-tensor compact `cygv` probe on those reduced rows is explicitly
   non-promotable: it gives `0` for targets `3`/`6` and rejects targets `7`/`8`
   as below compact CY dimension three. This keeps the blocker on deriving
   certified local/noncompact phase data.
   The target-vs-star support guardrail confirms the star rows cannot be
   assigned to the missing targets directly: each weighted target relation has
   relation points absent from the actual star support. The chamber map still
   has to explain how these neighboring rows determine the target GV history.
   A union-support diagnostic now separates the four weighted cases: targets
   `7`/`8` are integral in the union relation lattice with their star rows,
   while targets `3`/`6` have non-integral target coordinates. That is the next
   local chamber/history object to certify.
   A fresh read of
   `string_theory/mcallister_2107/latest_cytools/compute_kklt_iterative.py`
   confirms the sign convention used by Cyrus,
   `tau_target = c_i/c_tau + chi/24 - GV(t)`, but also confirms that script is
   not a hidden solution of the 4-214-647 mixed-target problem: its primal path
   maps `c_i` onto `basis.dat` and computes GV in the CYTools divisor basis,
   rather than constructing the `kklt_basis.dat` corrected-target vector that
   currently blocks Cyrus.
   `mcallister_gv_context` now accepts the schema-4 context export written by
   `mcallister_first_principles` and aggregates the closest known `q_N`
   residual statuses and degree splits. Bounded target-7 and target-8 runs on
   the current 4-214-647 export record the live residual split as the same
   source-derived degree-6 predecessor plus toric degree-2 differences,
   sharpening the next source-history task to the certified semigroup/chamber
   history around that shared lower composite term. The unfiltered bounded
   report now promotes closest-known residual source predecessors into their own
   queue: two unique source-derived predecessors occur three times, target 7/8
   share the degree-6 ray, target 6 has the analogous degree-4 ray, and every
   occurrence has shared origin-circuit witness relations plus
   `source_derived_full_facet_context`. The source-derived GV-history importer
   now requires that full witness/facet context before accepting a source value.
   The latest actual-local-cygv report promotes ready source skeletons through
   the `cygv` computation path: 20 imports match CMS formulas, 21 imports have
   no CMS formula check, 32 remain full-facet CMS formula imports, 58 rows have
   no integral matching CMS formula, and 115 rows have no origin-circuit witness.
   Accepted source-derived imports now split by local charge signature as
   `-1,-1,-1,1,1,1` 28 times, `-1,-1,-1,1,2` 41 times,
   `-1,-1,1,1` 3 times, and `-1,-1,-1,-1,1,3` once. The residual predecessors
   remain in the `-1,-1,-1,1,2` family and are now computed from the
   source-derived local bundle phase `O(-1)+O(-2)->P^2`: compact
   including-origin `q=[-1,-2,1,1,1]`, semigroup `[[1]]`, grading `[1]`, unit
   tensor `κ_000=1`, and positive-base chamber. The residual predecessor queue
   has no remaining local missing-input counts and all three occurrences are
   `ready_for_actual_cygv_call`; the remaining blocker is now the higher-degree
   compact qN history that uses those source rows, not the scalar value of this
   residual source family. The path-history report now makes this distinction
   explicit: scalar toric/source GV coverage is reported separately from
   compact `q_N` polynomial materialization, and the degree-10 first-layer
   predecessor pairs still have only scalar coverage or unknown scalar coverage,
   not exported compact `q_N` polynomials. Cyrus now vendors `cygv` 0.1.2 with
   a narrow qN-trace API and has quintic regressions proving the trace path
   preserves the degree-one GV `2875` for explicit and provided-generator
   domains. The target-7 path-support probe now records 11 small-domain qN
   polynomials with bounded term samples and ten sampled nonzero lower-class
   qN lookups; target 8 similarly records 6 small-domain qN polynomials and ten
   sampled nonzero lower-class qN lookups. Both targets still have `GV=0` and no
   target qN polynomial in these small domains. A reduced two-target aggregate
   collapses 17 qN-polynomial occurrences to 11 unique lower curves, with 6
   curves shared by both target-support domains and 4 unique curves showing
   domain-dependent qN term counts. The residual-subtarget probe now also
   takes the closest-known degree-8 residual itself as a small path-support
   target: target 7's residual uses support size 7 with 8 generators and returns
   `GV=-2` after materializing 10 qN polynomials, while target 8's residual
   uses support size 6 with 6 generators and also returns `GV=-2` after
   materializing 6 qN polynomials. The original degree-10 target-support
   domains in those same reports remain `GV=0` with no target qN polynomial, so
   this is a diagnostic for the lower composite history, not a corrected-target
   replacement. The parent-vs-subtarget qN comparison now records the residual
   polynomial domain dependence directly: target 7's residual has 3 qN terms in
   the degree-10 parent path-support domain but 1 as a residual subtarget, and
   target 8 has 2 parent-domain terms but 1 residual-subtarget term. Both
   samples are complete and marked `different_qn_term_counts`. The parent-only
   terms now carry source classification: target 8's extra term is the missing
   target itself, while target 7 has the missing target plus one degree-10
   generated term outside the exported source-ray context. Subtracting the
   degree-8 residual from each parent-only term shows all offsets are known
   degree-2 toric classes with `GV=-2`; target 7 sees both its own and target
   8's residual difference, while target 8 sees target 8's residual difference.
   The context report now aggregates those parent-only offset degrees separately
   from known-qN status, so future reports can distinguish low-degree known
   toric offsets from broader unknown qN history without mining each sample row.
   It also distinguishes disabled qN-domain comparisons from requested
   path-history probes that stopped at a closure limit: a current target-7
   schema-4 smoke report records
   `unavailable_exceeded_element_limit_10000__target_in_closure__no_predecessor_differences`,
   not `not_run`, for the rational-cone residual comparison. That report also
   shows the lower-seed decomposition `A+B+B` is present. The path-history
   report now materializes the bounded subset predecessor pairs from that
   decomposition directly; a target-7 schema-4 smoke report exposes four
   candidate splits with degrees `8+2`, `2+8`, `4+6`, and `6+4`. These rows
   preserve the toric/source-known status and seed-sum evidence for each side,
   and the report now aggregates their candidate-count, degree-split, and
   known-history pair-status buckets at top level. They remain diagnostic
   evidence for the missing compact parent history rather than promotable GV
   inputs. The same report now also promotes unknown sides of those lower-seed
   predecessor candidates into a unique diagnostic queue. Target `7` has two
   unique unknown composites at degrees `4` and `8`; target `8` has three, one
   at degree `4` and two at degree `8`. Every queued row is currently
   `not_source_degree_bounded_ray`, so the next missing object is not another
   uncovered source-ray leaf but the parent chamber/history that couples these
   composite rows to known toric/source sides. The queue now carries the
   first-generation seed-sum evidence for the unknown side: the degree-4
   composites split as known toric plus known toric, while the degree-8
   composites split as known toric plus known source, with the source side
   itself pair-reduced into known toric leaves. It also runs a cheap bounded
   seed-decomposition check for each queued composite without invoking another
   compact `cygv` call. A target-6 schema-4 smoke report shows both degree-8
   unknown composites are found as two-term bounded seed decompositions; one is
   twice the same degree-4 source seed, explaining why the first-generation
   reduced-seed split was absent there. Fresh target-7 and target-8 smokes show
   the same status for every queued lower-seed unknown composite, so the current
   blocker has moved from primitive source discovery to composite-history
   propagation through the parent chamber. The queued rows now also classify
   the bounded-decomposition terms by degree and known qN history: target 6's
   degree-8 composites are built entirely from degree-4 known source-derived
   leaves, while target 7/8 degree-4 composites are built from known toric
   degree-2 leaves and degree-8 composites are built from known source
   degree-6 plus known toric degree-2 leaves. When `--run-lower-seed-diamonds`
   is enabled, the same queue now runs the bounded decomposition diamond qN
   trace for these composites. The current target 6/7/8 smokes show that some
   mixed source+toric degree-8 diamonds materialize nonzero GV (`-2` or `2`),
   while repeated/paired known-leaf composites can collapse to zero/absent qN;
   two repeated cases expose nonintegral GW candidates (`-1/4`) rather than
   promotable integer GV history. The report now aggregates this as
   `integer_nonzero_gv`, `integer_zero_or_absent_gv`, and
   `blocked_noninteger_gw_candidates` candidate statuses. Thus known leaves
   alone are not a sufficient promotion rule; the bounded diamond trace is the
   next gate. The lower-seed predecessor split report now joins each known side
   with the bounded-diamond status of the unknown side, so target-level
   candidate pairs can be read directly. Fresh target-6 and target-8 smokes
   show the expected split between known-source/known-toric partners and
   bounded integer-nonzero, zero/absent, or noninteger-blocked composites.
   The same flag now also runs a two-generator candidate-pair diamond using
   the lower-seed predecessor and difference as the domain generators. Target 6
   reports pair domains that are either zero/absent or blocked by HKTY
   noninteger candidates, while target 8 reports four integer-nonzero pair
   domains with GV `2` and two HKTY-error pair domains. This distinguishes
   "the unknown side has a nonzero bounded diamond" from "the candidate pair
   computes the missing target." The report now aggregates candidate-pair
   diamond GV values, compressed errors, qN-trace statuses, noninteger counts,
   and target GW coefficient statuses; the target-8 smoke records GV `2` for
   the four integer-nonzero pair domains and one noninteger GW candidate for
   each of the two HKTY-error domains. It now also compares the candidate-pair
   diamond GV and target GW candidates against the independent CMS expected
   formula sum using exact rationals. Target `8`'s four `GV=2` pair domains all
   mismatch expected sum `3`, and target `6`'s zero/absent or HKTY-error pair
   domains also mismatch expected sum `3`, so these candidate-pair diamonds are
   explicitly nonpromotable rather than near-misses. The same report now
   records the signed target GW instanton-coefficient gap needed to reach the
   independent formula sum: target `8` uses pivot `203` with component `-3` and
   is short by `-3` or `-5` toward required coefficient `-9`; target `7` uses
   pivot `44` with component `2` and is short by `3` or `16` toward required
   coefficient `6`; target `6` uses pivot `47` with component `1` and is short
   by `6` on the HKTY-error branch toward required coefficient `3`. This keeps
   the next task on parent-domain subtraction history rather than scalar
   candidate-pair values.
   The broader path-support generator probe now serializes the same target
   `GW` formula-instanton balance and aggregate balance statuses even when
   compact HKTY fails. Fresh
   target-specific schema-4 runs show target `8` has path-support coefficient
   `-6`, formula-required coefficient `-9`, and missing coefficient `-3`;
   target `7` has `hkty_error` but still exposes coefficient `3`, required
   `6`, and missing `3`; target `6` has pivot `47` with component `1` but no
   target instanton coefficient, so it is explicitly classified as
   `path_support_gw_formula_instanton_target_coefficient_missing`. This moves
   the signed residual from pair-diamond-only evidence into the parent
   path-support diagnostic boundary.
   The path-history probe now also guards lower-seed pair operations before
   pair reduction or pair-sum enumeration when the target has more than
   `--lower-seed-pair-limit` degree-bounded seeds. The default is `1024`, and
   the flag can be raised deliberately for single-target high-budget probes.
   This keeps oversized targets observable instead of burning the whole timeout
   before JSON is written: a guarded target-`0` smoke and target-`2`
   generation-limited smoke now return
   `skipped_seed_pair_limit_1024` in about four seconds. A target-`8` control
   with `seed_count=720` remains below the default cap and still runs the
   lower-seed decomposition path, so the existing target-`7`/`8` diagnostics
   are not replaced by the guard. A guarded all-target batch with path-support
   generators and lower-seed diamonds now writes all nine per-target reports in
   229.70 seconds: targets `0..5` are explicit seed-pair-limit skips, target
   `6` preserves `hkty_error` with missing target coefficient, target `7`
   preserves `hkty_error` with missing formula coefficient `3`, and target `8`
   preserves the computed path-support mismatch with missing coefficient `-3`.
   The bounded decomposition search now uses sparse first-pair sums internally
   rather than storing dense vectors and every duplicate pair. With an explicit
   high-budget `--lower-seed-pair-limit 2048`, a target-`0` path-history probe
   with `--closure-generation-limit 1` now completes in 83.15 seconds,
   pair-reduces `1616` seeds to `702`, and finds a three-term lower-seed
   decomposition before hitting the bounded closure element limit. The same
   `--closure-generation-limit` option now accepts `0` to stop immediately
   after the initial seed set, so high-budget seed-decomposition probes can be
   separated from closure expansion. A target-`1` zero-generation smoke with
   `--lower-seed-pair-limit 2048` completed in 45.40 seconds, pair-reduced
   `1616` seeds to `702`, found a repeated-term three-seed decomposition, and
   exported four lower-seed predecessor candidates with status
   `stopped_generation_limit_0`. The zero-generation mode now maps the rest
   of the default-skipped target set at decomposition level without expanding
   the broad semigroup closure: target `2` completed in 107.10 seconds
   (`2560 -> 1124`, three terms, four candidates), target `3` in 191.04
   seconds (`2963 -> 1316`, three terms, four candidates), target `4` in
   85.22 seconds (`2212 -> 949`, three terms, four candidates), and target `5`
   in 210.63 seconds (`2560 -> 1124`, four terms, six candidates). All of
   these reports keep the closure status as `stopped_generation_limit_0`, so
   they are decomposition/history triage rather than certified compact HKTY
   histories. The lower-seed decomposition search now also has a
   positive-grading pair-degree filter for corrected-chamber probes. It first
   checks exact two-seed hits directly, then skips pair sums that are too high
   in grading degree to combine with any remaining positive seed. On the
   target-`5` zero-generation smoke this preserved the exact decomposition and
   predecessor sample while reducing runtime from 210.63 seconds to 128.75
   seconds. The path-history report now serializes these lower-seed search
   statistics per target and aggregates each statistic at top level. A
   target-`1` zero-generation smoke preserved the prior decomposition/
   predecessor sample, completed in 25.09 seconds, and reported matching
   per-target and top-level buckets for `direct_pair_seed_scan_count=1616`,
   `pair_sum_degree_bound=17`, `pair_sum_candidate_count=388643`,
   `pair_sum_unique_count=386251`, `pair_sum_duplicate_count=2392`, and
   `pair_sum_overdegree_skipped_count=917893`.
   The parent-only classifications now include parent-path-support runtime
   lookups from the same `cygv` run. On the regenerated target 7 report, the
   generated degree-10 side term is a real parent-domain nonzero object
   (`GV=-2`, one qN term, `integer_nonzero_gv`), but the missing-target-shaped
   parent-only term is `GV=0` with no qN polynomial and an
   `integer_zero_or_absent_gv` readout. Target 8's missing-target-shaped
   parent-only term has the same zero/absent parent-domain status. Thus the
   observed small domains still do not compute the missing target GV; they
   expose how lower residual history and a generated sibling term enter the
   parent residual polynomial. A follow-up qN-shape lookup shows the generated
   target 7 side term is an `identity_single_term_qn_polynomial`, so it has no
   hidden lower-history terms of its own in this parent domain. The
   missing-target-shaped parent-only terms for targets 7 and 8 have no
   parent-domain qN polynomial shape because those same parent domains read
   them as zero/absent. The target-coefficient balance now reconstructs the
   small-domain pre-subtraction target readout: target 7 has lower-source sum
   `3`, post-subtraction coefficient `0`, pivot component `2`, and therefore
   pre-subtraction candidate `3/2`; target 8 has lower-source sum `-6`, pivot
   component `-3`, and pre-subtraction candidate `2`. This pins the small
   domains as internally consistent zero histories, not as the missing
   corrected-chamber target values. The report now compares these reconstructed
   candidates against the CMS-general-divisor formula candidates as well:
   target 7 reports `3/2` versus formula `3`, and target 8 reports `2` versus
   formula `3`, so neither sampled domain is a disguised corrected-chamber
   formula domain.
   The same report is not yet scalable across all nine misses: the all-target
   schema-4 run timed out at 900 seconds, targets 2-5 time out under 180-second
   per-target probes, and targets 0/1/6 hit non-integer HKTY errors in their
   small path-support domains.
   The offset-generator qN-history diagnostic now exports the candidate
   generator degree buckets, sparse generator samples, and raw `FIND_GV=false`
   GW coefficient traces. Regenerated target `7` and `8` reports show that
   the matching target-`7` degree-six/degree-eight and target-`8` degree-six
   domains compute the expected source `GV=-2`, while target `8`'s
   degree-eight domain still fails the integral `cygv` run but reaches raw
   source candidate `2`. All four domains expose non-integral lower
   coefficient history (`13`, `8`, `9`, and `9` non-integer GW candidates,
   respectively), and the vendored `cygv` error identifies the first target-`8`
   degree-eight failure as
   `element_nonzero=[(203,-3),(206,1),(209,1)]` with candidate `4/3`; the
   classified GW sample marks that curve as `unknown_not_toric_covered` and
   `source_ray_matches_missing_target`. The full target-`8` degree-eight
   count map has all `9` fractional candidates in `unknown_not_toric_covered`,
   split as `8` `not_source_degree_bounded_ray` and `1`
   `source_ray_matches_missing_target`. All candidate domains remain without
   supporting-face certificates. The supporting-face diagnostic now records
   why: fresh target `7` and `8` reports classify all four offset-generator
   candidates as `lp_search_status=lp_no_certificate`,
   `exact_kernel_status=no_certificate`, and
   `aggregate_status=lp_cutting_round_limit` at the default `64` cuts; all
   `16` bounded anchor LP attempts per candidate also hit the cutting-round
   limit. The context binary now exposes
   `--supporting-face-lp-cutting-rounds` and
   `--supporting-face-lp-anchor-attempts`; rerunning target `7` and `8` with
   `256` cuts turns both degree-eight aggregate LPs into `lp_no_solution`,
   leaves the shared degree-six aggregate LPs at the cutting limit, and still
   yields zero anchor real-normal solutions. The direct full-constraint LP
   pass is sharper: all four offset-generator rows report
   `full_status=lp_no_solution`, so these sampled qN-shape domains are not
   supporting faces of the exported degree-bounded Mori cone. This keeps the
   blocker at certified chamber-semigroup construction, not scalar source-GV
   recovery.
   Applying the trace to McAllister still requires the certified corrected-chamber
   semigroup/history domain.
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
cargo build --release -p cyrus-core --bin mcallister_gv
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 timeout 600 ./target/release/mcallister_gv --min-points 20000
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_MIRROR_GV_HEAVY=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 cargo test -p cyrus-core --test mcallister_e2e stage5_mirror_gv_checkpoint_matches_cygv_min_points -- --nocapture
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
   subchecks together. The exact helper
   `apply_certified_flop_continuation_sequence` now applies an ordered certified
   path to `kappa`, `c2`, and the GV table together. These are necessary pieces
   for the stable-Weyl test, but
   they are not a chamber certificate by themselves. Cyrus still needs the
   source-level certification step that proves the supplied finite generators
   are the relevant Kähler wall/chamber context and identifies the shrinking
   divisor from geometry rather than from a fitted candidate list. The
   LP-witness face diagnostic is useful but currently uncertified, so it must
   remain a diagnostic rather than a fallback. Cyrus now has a reusable
   LP-assisted supporting-face certificate search for higher-codimension
   supports, but it only promotes normals that pass exact integer verification.
   On the schema-3 corrected-chamber context, the bounded origin-support pass
   found no promoted relation-support certificate. LP no-solution and finite
   solver statuses such as `Unknown` are now reported as no-certificate rather
   than hard diagnostic errors, and a raised `4096`-generator guard checks the
   degree-10 facet-union domains directly: target `7` has
   relation/shared/union ranks `1/13/194`, target `8` has ranks `1/9/177`, and
   none is promoted to a supporting face.
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
   The support-overlap runner now has an opt-in `--trace-support-overlap-qn`
   mode that uses the same provided-generator `cygv` qN trace boundary as the
   smaller path-support probe and exports target qN materialization fields when
   the call finishes. A quintic-sized regression proves the flag reaches actual
   `cygv` qN history and recovers the degree-one `2875` target polynomial; the
   broad McAllister degree-10 support-overlap domains remain too large to use as
   the missing chamber history. Fresh schema-4 target-specific traces also
   rule out the smaller overlap windows for the two degree-10 misses: target
   `7` thresholds `1..4` use `64`, `41`, `24`, and `10` generators, target `8`
   uses `48`, `40`, `18`, and `10`, and all eight traced calls fail in
   `cygv` series inversion with non-integer GV output. By contrast, the
   path-support domains are internally integral and materialize compact lower
   `q_N` polynomials (`11` for target `7`, `6` for target `8`), but both still
   return target `GV=0` with no target `q_N` polynomial. This keeps the blocker
   at certified corrected-chamber history, not qN observability. The path
   support report now also identifies lower `q_N` polynomials containing the
   target monomial directly. For both degree-10 targets there are three such
   lower sources: the toric degree-two predecessor, the shared source-derived
   degree-six row, and the unknown degree-eight residual, with target-term
   coefficients `-2`, `-1`, and `1`. The
   `cygv_closest_known_qn_residual_difference_*` fields now aggregate the
   residual difference side explicitly: target `7` leaves known toric
   `[(44,1),(54,1),(206,-2)]`, while target `8` leaves known toric
   `[(54,1),(203,-2)]`, both paired with the shared source-derived degree-six
   predecessor. Thus the unresolved degree-eight residual is visible as a
   composite qN-history gap rather than an unknown scalar GV leaf. The qN/Li2
   trace now comes from `cygv`'s exact
   `Li2(q_N)` polynomial after finite monomial-map truncation: the exact Li2
   target coefficients are `-5/2`, `-1`, and `1`; target `7`'s pivot
   subtraction coefficients are `5`, `0`, and `-2`, while target `8`'s are
   `-10`, `2`, and `2`. This exposes the series-inversion
   cancellation/history that must be reproduced in a certified chamber. The
   same trace now records the target-degree GV candidate read from cygv's
   mutable instanton polynomial: in both small path-support domains the target
   instanton coefficient and GV candidate are exactly `0`, with
   `integer_zero_or_absent_gv` status. The current small domains are therefore
   internally consistent zero histories, not failed target lookups. The
   path-support report also exposes all GV coefficient readout status counts:
   target `7` has `79` readouts (`11` nonzero, `25` zero/absent, `43`
   missing), while target `8` has `33` (`6`, `16`, `11`). This pins the
   composite lower-degree decisions made in the small domains.
   The refreshed parent-only aggregates now show target `7`'s two extra
   parent-domain terms both have lower-seed decompositions whose bounded
   diamonds compute `GV=0` with no target qN polynomial, while target `8`'s one
   parent-only term has a lower-seed decomposition whose bounded diamond hits
   an HKTY error. This rules out a simple "lower-leaf scalar missing" diagnosis
   for those parent terms. The error bucket is now aggregated as well; the
   detailed target `8` sample is the explicit-semigroup non-integer candidate
   `4/3` for `[(203,-3),(206,1),(209,1)]`.
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
   The degree-ladder report now also annotates seed-count-skipped steps with a
   bounded streaming closure summary. A target-`7` debug smoke with
   `--semigroup-measure-max-seeds 64
   --cygv-degree-ladder-bounded-closure-limit 1024` records degree `2`
   as a completed bounded closure with `291` elements, while degrees `3..6`
   all exceed the `1024` element cap during the first closure generation. This
   gives an inexpensive machine-readable scale witness when the actual cygv
   constructor is deliberately skipped.
   Running the visible lower-seed diamonds does not provide a fallback. On the
   actual-local-cygv lower-diamond report, the raw lower-seed domains split into
   six computed `GV=0` diamonds and three non-integral HKTY errors, while the
   pair-expanded reduced-seed diamonds compute `GV=0` for all nine missing
   targets.
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
   signature `[-1,-1,-1,1,2]`, matching CMS checks, shared origin-circuit
   witness relations, and `source_derived_full_facet_context`. Cyrus now
   recognizes the corresponding source-derived local
   `O(-1)+O(-2)->P^2` phase, supplies tensor `κ_000=1` and a positive-base
   chamber, and imports the value through actual local `cygv`, requiring
   agreement with the CMS formula where one exists. The residual predecessor
   queue also shows target `6` has an analogous degree-four source predecessor
   with the same local signature and ready local-cygv source inputs. The
   actual-local-cygv import-status report has `20` actual local cygv imports
   matching CMS formulas, `21` actual local cygv imports without CMS formula
   checks, `32` remaining full-facet CMS formula imports, `58` rows with no
   integral matching CMS formula, and `115` rows with no origin-circuit witness.
   The same local unit-phase probe now records qN trace metadata for the
   certified unit-tensor source call; the residual-source regression verifies
   one materialized target qN polynomial for the source-local `GV=-2` row.
   This is trace evidence for the source-local model, not a promoted compact
   corrected-chamber qN polynomial.
   The fresh schema-3 local-source qN trace report confirms the distinction:
   all nine target-level unit probes now expose one materialized local qN
   polynomial, but all nine target skeletons are still blocked on
   source-derived tensor/chamber certification.
   The CMS divisor cubic suggested tensor value `3`, but the raw-cubic primitive
   probe gave `GV=-6`; the source-derived unit tensor gives the expected `GV=-2`
   and is now the promoted source value for this lower ray. The follow-up
   primitive probes show that unit tensor normalization explains the
   `[-1,-2,1,1,1]` lower ray, and that omitting the origin/canonical divisor
   column explains the
   resolved-conifold-like `[-1,1,-1,1,-1,1]` lower ray. These are still
   phase/chamber diagnostics, not promoted GV values. Across all nine
   unresolved origin circuits, the origin-omitted shape aggregate reports five
   dimension-three no-origin rows and four CY-dimension-2 rows. The phase
   selector now labels those five no-origin rows as non-Calabi-Yau diagnostics,
   not compact CY threefold handoffs, because their charge sums are nonzero. The
   target-level unit-tensor phase probe matches the four CY-dimension-2 rows in
   the origin-included convention and only one no-origin diagnostic row in the
   origin-omitted convention; the other four no-origin rows still mismatch their
   formula candidate sets.
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
   `n'^0_C=0` while rejecting conflicting duplicate data.
   `apply_certified_flop_continuation_sequence` applies these three exact
   updates in lockstep for an ordered certified path. The remaining work is the
   certification layer: Cyrus still needs to identify the shrinking curve and
   certify its `n_C^0` before these transforms can drive a corrected-chamber
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
   The first-principles context export now includes the full degree-bounded
   toric GV diagnostic context for this missing-GV degree range. On the current
   solved-t release export this adds `1053` source-derived toric diagnostic
   rows and reduces the shared-facet unknown bucket from `228` to `36`
   occurrences. The all-domain unresolved queue drops to `1631` unique
   generators over `7360` occurrences. The report now splits out the narrow
   shared-facet queue: `33` unique non-known generators over `66` occurrences,
   consisting of the nine missing targets, `4` uncovered-source-ray hits, and
   `20` source degree-bounded rays still lacking toric/source-derived GV values.
   The first-principles export now also emits complete stats for those `20`
   source rays. They are all origin-circuit rows, all still lack the certified
   local q-matrix phase/intersection tensor/chamber certificate, and their CMS
   checks expose `4` integral inferred-degree matches with formula value `3`,
   `2` non-integral inferred-degree matches, and `22` no-solution checks. The
   shared-source unit-tensor probe matches only `5/20` expected formula sets
   (`GV=-2` cases); formula-`3` rows give unit candidate `6` and
   origin-omitted candidate `-9`, so tensor/chamber certification remains the
   blocker. The CMS raw-divisor-cubic primitive probe also fails for the four
   integral formula-`3` rows, producing `-36`, `-42`, `-42`, and `-48`, so the
   needed tensor is neither the unit tensor nor the raw divisor cubic.
   The context report now also aggregates raw GW coefficient diagnostics for
   bounded-decomposition diamonds built from explicit semigroups. On the fresh
   target-`8` path-support report, the parent-only bounded diamond still fails
   compact HKTY integrality, but the raw trace localizes the failure to one
   nonintegral target-row candidate: `GV` candidate `4/3`, classified as
   `source_ray_matches_missing_target` and `unknown_not_toric_covered`. This
   keeps the target-`8` diamond in the diagnostic bucket rather than promoting
   it, while narrowing the remaining chamber-history problem.
   The analogous target-`7` run keeps both parent-only bounded diamonds at
   `GV=0`; their target raw GW statuses are `missing_instanton_coefficient` and
   `zero_or_absent_gw`, with only one unrelated fractional lower degree-4
   candidate. This contrast separates target `8`'s fractional target-row
   failure from target `7`'s zero-history failure.
   Target reports now embed the local unit probe and local integer tensor
   scan. For target `8`, that puts the local scalar check next to the compact
   diamond evidence: the oriented local q-matrix `[-1,2,-3,1,1]` with tensor
   value `1` gives candidate `3`, matching the CMS formula sum, while the
   compact bounded diamond reads the target coefficient as `4/3`. The remaining
   problem is therefore not the one-parameter orientation itself, but the
   certified chamber/intersection tensor history that would make such a local
   value promotable in the compact pipeline.
   The report now makes that boundary explicit with
   `local_cygv_one_parameter_family_status_counts`. The regenerated all-target
   report classifies the four compact one-parameter target families as
   `uncertified_one_parameter_split_bundle_over_weighted_p2:base=1,1,2;bundle=1,3;base_hyperplane_square=1/2`,
   while certified
   `source_derived_local_p2_bundle_family_with_tensor_chamber_certificate` rows
   are counted separately. Thus the target `8` unit-tensor/formula-sum match
   remains diagnostic and cannot be confused with the source-derived
   `O(-1)+O(-2)->P^2` import family; it now has the more specific
   `O(-1)+O(-3)->P(1,1,2)`-style weighted-base shape. The target skeletons now
   also report
   `local_intersection_tensor_blocked_weighted_p2_split_bundle_requires_source_derived_resolution_chamber`
   and the analogous chamber-certificate blocker for those four rows, so the
   missing input is labeled as a resolution/chamber certificate rather than a
   generic non-`P^2` mismatch.
   Origin-circuit witness-domain reports now optionally run compact `cygv` on
   relation, shared-facet, and guarded facet-union generator domains. For
   target `8`, relation-only gives `GV=2`, materializes one target qN term,
   and has raw target GW candidate `2`, while the shared facet gives `GV=0`
   after `24` qN polynomials with no target qN and raw target GW candidate
   `0`. For target `7`, relation-only fails integrality on the missing target
   row with raw candidate `3/2`, while the shared facet gives `GV=0` after
   `51` qN polynomials with no target qN and raw target GW candidate `0`.
   Both large facet unions exceed the default `64`-generator compact probe
   guard. The certificate diagnostics now also report the LP phase outcomes:
   target `8` shared facets have full/aggregate `lp_no_solution` and anchors
   split `14` cutting-round limits to `2` no-solution outcomes; target `7`
   shared facets have full LP `Unknown`, aggregate `lp_no_solution`, and
   anchors split `9` cutting-round limits to `7` no-solution outcomes; both
   facet unions have full/aggregate/anchor `lp_no_solution`. These readouts
   are diagnostic because the witness-domain supporting-face certificates
   remain absent, but they show the small witness domains are not the missing
   corrected chamber. A smoke run with
   `--supporting-face-lp-anchor-attempts 2 --supporting-face-lp-cutting-rounds 1`
   now confirms the origin-witness diagnostics honor the CLI LP search limits
   instead of using hardcoded defaults. Raising the certificate-only budget to
   `64` anchors and `256` cutting rounds keeps target `8` shared/union domains
   at full/aggregate/anchor `lp_no_solution`; target `7` union domains match
   that, and target `7` shared domains retain full LP `Unknown` but have
   aggregate and all `64` anchors at `lp_no_solution`. The new
   `--scan-origin-witness-span-closure` diagnostic rules out one cheap
   alternative for the full nine-target queue: after switching to a per-domain
   exact rational row-echelon span basis, all relation/shared/union domains
   report `span_closed_under_degree_bounded_context` with zero extra same-span
   generators across `720..2963` unique degree-bounded candidates per target.
   So the issue is not that the small witness domains forgot forced row-span
   members; it is still the absent supporting-face/chamber-semigroup
   certificate. The star-union diagnostic now also serializes exact rational
   coordinates for the four weighted rows: targets `3`/`6` are half-integral in
   the union relation lattice (`[1/2,1/2,-3/2]` and `[1/2,-2,3/2]`), while
   targets `7`/`8` are integral (`[3,0,-1]` and `[1,0,-1]`). This keeps the
   next chamber-map step focused on a real lattice/projection obstruction
   rather than a missing serialization detail. The context consumer now applies
   CYTools' no-origin `q`-matrix column convention and projects all four
   weighted target/star rows plus target±star combinations to integral global
   curve-basis classes; the target projections match the missing-target rows
   exactly. The remaining reusable fix is therefore not the global column map,
   but the local chamber/intersection certificate needed before those rows can
   determine a GV value.
   The subsequent global-basis lookup confirms this is still a chamber-history
   object, not an already-known scalar GV row: every projected target/star and
   target±star role for targets `3`, `6`, `7`, and `8` is
   `unknown_not_toric_covered`. The target degrees are `26/12/10/10`, the
   star degrees are `-36/-24/-4/-4`, the target-minus-star degrees are
   `62/36/14/14`, and the target-plus-star degrees are `-10/-12/6/6`. The
   negative star-side degrees rule out treating the star relation as an
   effective covered GV contribution.
   The sign-reversed lookup makes that conclusion sharper. For targets `3` and
   `6`, the positive-degree opposites of the negative star and target-plus-star
   rows are still `unknown_not_toric_covered`. For targets `7` and `8`, the
   opposite star row is known source-derived GV `1` at degree `4`, but the
   target and target-plus-star rows still have no known toric/source-derived
   scalar value. So the adjacent source row is evidence for the wall history,
   not the corrected target GV itself.
   The lookup now carries source provenance too. The target `7`/`8`
   opposite-star source is the resolved-conifold origin relation
   `[(0,-1),(55,-1),(195,1),(212,1)]`; target `6`'s opposite star is a
   degree-bounded source ray without known GV, and target `3`'s opposite star is
   not a degree-bounded source ray.
   A bounded lower-seed check now runs only for star-union rows of degree at
   most `6`, keeping the diagnostic cheap. It decomposes the target `7`/`8`
   opposite-star row into two known degree-2 toric/source rows, but the
   degree-6 `target+star` row is `not_found_up_to_4` and is still not a
   degree-bounded source ray. This makes the next compact-history target the
   degree-6 chamber row, not the resolved-conifold side.
   The report also exports the union charge-basis coordinates for the sum row:
   target `7` has `target+star=[3,1,-1]`, target `8` has
   `target+star=[1,1,-1]`; both have sorted point-coefficient signature
   `[-3,-1,-1,1,1,1,2]`.
   A unit-tensor diagnostic on the explicit point-row charges now returns toy
   GV `0` for both targets. The report labels the compressed one-parameter
   row as `not_compact_threefold_shape_cy_dim_5`, so it stays diagnostic and
   does not close the corrected-chamber GV gap.
   The star-union report now has an opt-in target-plus-star path-history probe,
   so this chamber-row check no longer requires a synthetic context. On target
   `8`, `--probe-star-union-path-history --closure-generation-limit 1
   --element-limit 150000 --lower-seed-pair-limit 2048` completed in 16.55s
   and reported the degree-6 `target+star` row with `463` seeds, `355`
   reduced seeds, `65538` first-layer closure elements,
   `target_in_closure=false`, `not_found_up_to_4` lower-seed decomposition,
   and top-level status bucket `stopped_generation_limit_1`. The same probe now
   uses streaming element-limit enforcement, so the formerly timing-out
   generation-2 run completes in 24.50s with status
   `exceeded_element_limit_150000_during_generation_2`, `150000` retained
   closure elements, generation-2 partial count `84462`, and
   `target_in_closure=false`. This pushes the next task back to certified
   chamber semigroup/history construction.
   The shared-face fan over the full star-union support has now been checked
   directly as a candidate chamber certificate: regenerating
   `/tmp/cyrus_gv_context_star_union_shared_face_secondary_report.json`
   completed in 23.40s and reports all four weighted-P2 rows as
   `star_union_shared_face_secondary_on_wall`. Each has four shared-face
   simplices, six native secondary-cone hyperplanes, min/max pairing `0`,
   zero-pairing count `6`, no negative pairings, and neither sign of the
   source-derived height vector strictly inside. Thus the obvious fan combining
   the two target-exclusive simplices with the two star simplices is wall data,
   not the missing corrected-chamber interior certificate. The report now also
   serializes the exact point-index circuit supports for those six walls; for
   targets `7`/`8`, one wall is the resolved-conifold star-side relation
   `[(0,-1),(55,-1),(195,1),(212,1)]`, so the obstruction is visible as
   degenerate wall circuitry rather than just an aggregate count.
   The first-principles corrected-chamber context exporter now also carries the
   full corrected secondary-fan height vector as
   `secondary_cone_heights_for_missing`, and `mcallister_gv_context` validates
   its origin-included length, finite entries, and zero origin normalization
   against the no-origin q-matrix width plus one. That gives the next pass the
   actual chamber height data needed to compare these local wall circuits
   against the corrected global chamber, rather than only against the
   relation-support affine height.
   A fresh first-principles context dump at
   `/tmp/cyrus_corrected_chamber_gv_context_global_heights.json` confirms the
   origin-included height convention (`219` heights vs. `218` no-origin
   `q`-matrix columns). The corresponding
   `/tmp/cyrus_gv_context_global_height_report.json` evaluates those local wall
   circuits on the actual corrected global height vector: all four weighted rows
   cross oriented walls, each with six nonzero global-height pairings split as
   four positive and two negative. The target-`7`/`8` resolved-conifold
   star-side circuit `[(0,-1),(55,-1),(195,1),(212,1)]` has global pairing
   `0.42571113815644424`, so the visible shared-face wall fan is not compatible
   with the corrected global chamber height.
   The follow-up `/tmp/cyrus_gv_context_global_regular_report.json` computes
   the local regular triangulation induced by that same corrected height on the
   star-union support. All four weighted rows match the serialized star extras
   over the shared face (`[2,46]` for targets `3`/`6`, `[195,212]` for targets
   `7`/`8`) and none selects the target-exclusive relation points over that
   face. Thus the serialized star simplices are the actual local chamber chosen
   by the corrected global height, and the missing target relation requires a
   chamber-transport or wall-crossing computation rather than a hidden local
   triangulation swap.
   `/tmp/cyrus_gv_context_target_relation_height_report.json` now also pairs
   each missing target origin-circuit relation with the corrected global height
   vector and compares it to the branch `q.t` diagnostic. All nine missing
   targets are positive and all nine match `q.t` to tolerance. This confirms
   the target relations are genuine positive-volume classes in the corrected
   chamber even when the weighted shared-face regular triangulation selects the
   star extras, so the remaining blocker is the chamber-transport/qN history
   between those data rather than a hidden zero-wall relation.
   The same report now records corrected global-height pairings on each
   star-union global-basis lookup row. For targets `7`/`8`, the target rows are
   degree `10` and positive, the target-minus-star rows are degree `14` and
   positive, the target-plus-star rows are degree `6` and positive, and the
   shared star row is degree `-4` with pairing `-0.42571113815644424`; its
   opposite is degree `4`, positive, and has known source-derived GV `1`. That
   makes the remaining transport object explicit in the same rows as degree
   and known-qN status.
   `/tmp/cyrus_gv_context_transport_report.json` now serializes the algebraic
   transport identity itself. Eight projected rows satisfy
   `target = target_plus_star + opposite(star)`. For targets `7`/`8`, this
   writes the unknown degree-10 target as unknown degree-6 target-plus-star
   plus the known degree-4 positive opposite-star source row with GV `1`, and
   the corrected height pairings add back to the target `q.t`. The remaining
   missing object has therefore narrowed again to the positive degree-6
   target-plus-star qN/chamber history.
   The component extremal-ray probe now runs on these transported rows under
   the same `--certify-target-extremal-rays` guard used for original targets.
   It performs an exact same-positive-ray scan before any expensive separator
   search. The current all-target smoke report,
   `/tmp/cyrus_star_union_extremal_smoke_all.json`, shows all nine original
   targets are non-extremal by exact serialized decompositions. It also shows
   the positive target-plus-star components are not finite Mori generator rays:
   six report `not_certified_no_same_positive_ray_generator`, two are zero
   target-plus-star rows, and one lacks a target projection. For the concrete
   target `7`/`8` blocker, target-plus-star has same-ray count `0`, while the
   known opposite-star component has same-ray count `1` but is skipped by the
   default full-cone separator guard (`2963` generators versus limit `256`).
   This rules out another false path: the target-plus-star row is not an
   extremal finite-cone generator waiting for a scalar GV value.
   The same report now checks whether the target-plus-star point support
   contains a finite degree-bounded Mori generator set that could be certified
   as a supporting face. For targets `7`/`8`, the support generator count is
   `0` up to the degree-6 target-plus-star bound and the face status is
   `origin_support_no_generators`. This also rules out treating the
   seven-point target-plus-star support as an undiscovered finite face
   semigroup in the exported degree-bounded context.
   `/tmp/cyrus_star_union_wall_circuit_smoke_all.json` now links the known
   opposite-star component to the local wall geometry directly. For targets
   `7`/`8`, the opposite-star source-derived GV `1` class exactly matches the
   shared-face wall circuit `[(0,-1),(55,-1),(195,1),(212,1)]`; its local
   secondary pairing is `0` and its corrected global-height pairing is
   `0.42571113815644424`. Thus the known degree-4 component is certified as
   visible wall-circuit data in the star-union support, while the unresolved
   degree-6 component remains a chamber/semigroup transport problem.
   `/tmp/cyrus_star_union_wall_transport_readiness_smoke_all.json` now makes
   that handoff machine-readable at report level. It reports two rows in
   `wall_transport_known_wall_remainder_requires_chamber_semigroup`, six rows
   still blocked by missing known wall-circuit history, and one row blocked
   before the transport identity. The two ready-but-unresolved rows are targets
   `7`/`8`: both have the crossed wall curve
   `[(0,-1),(55,-1),(195,1),(212,1)]`, crossed-wall GV `1`,
   target-plus-star support generator count `0`, support-face status
   `origin_support_no_generators`, and missing inputs
   `target_plus_star_qn_history`,
   `target_plus_star_chamber_semigroup_transport`, and
   `shrinking_divisor_or_flop_certificate`. This is a blocker certificate, not
   a GV computation shortcut.
   `/tmp/cyrus_star_union_wall_branch_smoke_all.json` now additionally records
   the B-field branch status of the same crossed wall. The targets `7`/`8`
   wall curve has `q.t=0.42571113815643002`, parity mod 2 equal to `0`, and
   positive-side dilog status `real_ok`, but the negative-side continuation is
   `crossed_wall_negative_continuation_real_branch_cut`. This rules out a
   naive real-axis flop continuation of the known wall GV into the star-side
   chamber; the missing handoff still has to be a certified chamber/semigroup
   history or a source-level continuation that addresses the branch cut.
   `/tmp/cyrus_target7_stable_weyl_rank_probe_report.json` now threads the crossed
   wall through the stable-Weyl/flop-certificate gate. For target `7`, the
   report bucket is
   `stable_weyl_blocked_cms_general_divisor_no_rational_divisor_solution`:
   the crossed-wall source sample has shape candidate divisor `55`, but its
   CMS divisor check has no rational divisor solution. The exposed linear-system
   diagnostic is exact rank obstruction, not a floating mismatch:
   `row_count=219`, `column_count=214`, `rank=6`, `augmented_rank=7`. The
   opposite-wall extremal check is also still guard-skipped at `2963` generators
   under the smoke limit. Thus the exact stable-Weyl algebra is wired into the
   blocker, but it still has no source-derived shrinking divisor to check or
   promote.
   `/tmp/cyrus_gv_context_target_plus_star_support_report.json` now reconstructs
   the point-support charge bases for that degree-6 component. Targets `7`/`8`
   are integral two-parameter relations on affine-rank-4 seven-point supports:
   target `7` has coordinates `[3,-1]`, target `8` has coordinates `[1,1]`.
   This rules out treating the residual as a one-parameter scalar source
   formula; the missing computation is a certified two-parameter local
   chamber/qN history.
   `/tmp/cyrus_target7_cicy_readiness_report.json` and
   `/tmp/cyrus_target8_cicy_readiness_report.json` now test that two-parameter
   support against the actual `cygv` shape contract. Both target-plus-star
   handoffs are still `blocked_missing_source_derived_inputs`, but the blocker
   is sharper: targets `7`/`8` are not compact CY3 hypersurface inputs
   (`q_rows=7`, `q_cols=2`, hypersurface `cy_dim=4`), yet they are
   codimension-2 complete-intersection CY3 shape candidates (`cy_codim=2`,
   `cy_dim=3`). The readiness status is now
   `target_plus_star_local_cygv_blocked_complete_intersection_nef_partition_not_certified`.
   Target `7` has no positive primitive orientation for `[3,-1]`; target `8`
   has positive `[1,1]` coordinates but still needs a source-derived
   nef partition, complete-intersection intersection tensor, semigroup/grading
   certificate, and chamber/tensor data. The reports enumerate the codimension-2
   local bipartitions: both targets have `63` possible nontrivial splits and
   `6` zero-degree split candidates, so the CICY status is
   `complete_intersection_cy3_shape_ambiguous_zero_degree_nef_partition_candidates_requires_source_rule`.
   The report now also carries target-relation balance data for these
   candidates. For both targets, all six zero-degree candidates have
   `target_relation_balanced_inside_each_part`, so the integral target relation
   does not pick a unique nef split either.
   The reusable CYTools-style CICY tensor-reduction step now exists in core as
   `compute_complete_intersection_cy3_intersection_numbers`, but targets `7`/`8`
   still lack the source-derived nef partition and ambient top intersections
   needed to call it.
   After recomputing the affine charge basis, every single-column omission for
   targets `7`/`8` remains Calabi-Yau charged but has `q_rows=6`, `q_cols=1`,
   and hypersurface `cy_dim=4`, so the next fix is a source-derived CICY
   nef-partition/reduction or chamber-transport history, not a direct raw
   hypersurface `cygv` call or a simple point deletion.
   The full union-support height profile is now serialized as well. For
   targets `3`/`6`, the off-hyperplane zero-shared ray carries zero target and
   star coefficient; for targets `7`/`8`, the height `-1` zero-shared ray and
   height `-1` star extra `212` carry opposite star coefficients. This keeps
   the next source-derived task focused on the actual chamber map and
   intersection transport across that wall.
   A direct off-height global projection lookup now rules out another shortcut:
   the target `7`/`8` star-side off-height pair `[(55,1),(212,-1)]` is not an
   integral global curve-basis class, and targets `3`/`6` have no nonzero
   off-height component. The off-height data is local wall data, not a hidden
   scalar GV row.
