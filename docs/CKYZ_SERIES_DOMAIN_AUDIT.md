# CKYZ Series-Domain Audit

This note records the current source-read result for the McAllister potent-ray
GV bottleneck. It is intentionally narrower than `CYGV_AUDIT.md`: the question
here is what finite coefficient domain is needed to compute a saved ray row
without replaying `potent_rays_gv.dat`.

## Source Facts

The relevant `cygv` path is:

1. `Semigroup::{with_max_degree,with_min_elements}` builds a finite affine
   semigroup by closing supplied generators under addition, then sorts elements
   by grading degree.
2. `PolynomialProperties::new` builds a monomial lookup map from exactly those
   semigroup elements.
3. `Polynomial::mul` only keeps a product if the sum of the two monomial
   exponent vectors is present in that monomial map. Products outside the
   finite semigroup are intentionally dropped.
4. `fundamental_period::compute_omega` computes `c0`, `c1`, and `c2` only for
   semigroup elements with zero, one, or two negative GLSM intersections.
5. `instanton::compute_instanton_data` forms
   `alpha = c0_inv*c1`, `beta = c0_inv*c2`, then
   `f_ab = beta_ab - alpha_a alpha_b` and contracts against the intersection
   form, with the diagonal `1/2` convention.
6. `series_inversion::invert_series` processes semigroup elements by distinct
   grading degree. At each degree it reads the current instanton coefficient,
   rounds an integral GV candidate, computes `qN` and `Li2(qN)` for nonzero
   candidates, and subtracts those dilogarithm contributions from the remaining
   higher-degree instanton polynomials.

The `previous_qn` cache inside `series_inversion.rs` is a performance hint for
choosing a closer starting polynomial when computing `qN`; correctness still
comes from multiplying by the necessary `exp(alpha)` correction. The important
correctness dependency is the degree-ordered subtraction history, not the cache
length.

## Consequence For Cyrus

The current Cyrus CKYZ local path computes a flat-coordinate instanton
potential and then applies the CKYZ finite-limit multiple-cover relation. That
is the right source family for rank-two local surfaces, but the finite monomial
domain still controls whether a target coefficient is exact.

For a target degree `T`, Cyrus must include all source-coordinate degrees that
can affect the coefficient of `q^T` through:

- log-period coefficients used by the inverse mirror map `z(q)`;
- products in `beta - alpha alpha`;
- composition of the contracted B-model series into flat coordinates;
- nontrivial multiple covers of `T`.

The current `target_downset` domain is componentwise safe but often too broad.
The newer generated causal domain matches the cygv "finite semigroup" shape
better, but it is still keyed only by source generators and a grading cutoff.
For large McAllister target directions that can be much larger than the actual
coefficient dependency set.

## Non-Goals

- Do not call compact `cygv` directly on a local surface charge matrix. The
  compact `cygv` fundamental-period code rejects effective CY dimension below
  three, so this is the wrong model for local rank-two surfaces.
- Do not replace the remaining work with a coefficient-pattern table. The local
  support pattern identifies the source model; it does not compute the GV row.
- Do not keep raising `CYRUS_CKYZ_MULTIPLES_TO_CHECK` against the broad domain.
  The unfiltered `N=5` gate exceeded a 600 second timeout, which is a domain
  design issue rather than just a small optimization problem.

## Next Implementation Target

The next useful Rust change should be a coefficient-dependency domain builder
for local CKYZ surfaces:

1. Start from the requested target degrees and their nontrivial covers.
2. Walk the formal operations used by
   `compute_ckyz_local_instanton_potential_corrections_domain` backwards:
   flat-coordinate substitution, inverse mirror-map exponential, and
   `beta - alpha alpha`.
3. Add only degrees whose products/compositions can contribute to the requested
   target coefficients.
4. Assert that the resulting sparse domain produces the same GV values as the
   existing broad domain on the local `P2`, `F0`, `F1`, and polygon-5 tests.
5. Only then raise the McAllister rank-two potent-ray gate beyond the current
   first-four default.

This is separate from the later rank-four task. Rank-four saved potent rays
need their own local face/semigroup context before any GV row comparison is
meaningful.

## Current Guardrail And API

Cyrus now has an internal test-only support tracer,
`ckyz_observed_support_domain_for_degrees`, that records the nonzero supports
of the current broad local CKYZ calculation and rebuilds a finite monomial
domain from those supports. The unit test
`ckyz_observed_support_domain_recomputes_targeted_f0_ray` checks that this
domain recomputes the same targeted F0 skew-ray GV values and does not add
terms outside the componentwise target downset.

Cyrus also has an explicit support-predicted public API,
`compute_ckyz_local_gv_invariants_for_degrees_with_predicted_support_domain`.
It computes exact alpha/beta supports on the broad source-derived downset, then
runs the inverse mirror-map and potential-composition steps at the support level
rather than evaluating the full rational coefficient series before building the
final computation domain. The test
`ckyz_predicted_support_domain_covers_observed_f0_ray_support` verifies that the
support-predicted domain covers the observed support and reproduces the same
targeted F0 skew-ray GV values.
`ckyz_predicted_support_domain_covers_observed_polygon5_ray_support` applies
the same check to the rank-three polygon-5 local model for the first two
multiples of the McAllister-style `[4,3,2]` source direction.
The API-level tests
`ckyz_predicted_support_domain_api_matches_target_downset_for_f0_ray` and
`ckyz_predicted_support_domain_api_matches_target_downset_for_polygon5_ray`
check equality with the older broad target-downset extractor.

The McAllister rank-two CKYZ gate now uses this support-predicted API and still
passes the default all-row `N=4` validation, with `potent_rays_gv.dat` used only
as the assertion. The elapsed time increased to about 132 seconds on the latest
run, so this is a correctness/structure promotion rather than a performance
breakthrough.

This is still not the completed coefficient-domain solution. A fresh source
read of `cygv-0.1.2` shows why: cygv computes `exp(alpha)` once and then uses
degree-ordered `q_N`/`Li2(q_N)` path history during series inversion, while the
current local CKYZ path still materializes the full rational inverse mirror map
`z(q)` before extracting coefficients. The support predictor now avoids a
repeated fixed-point support recomputation when all alpha support closures are
identical and caches duplicate alpha-product support multiplications. On the
focused polygon-5 `[4,3,2]`, `N=10` domain profile this changes the
support-domain run from timing out earlier, then 81.7 seconds, to 26.2 seconds.
Adding a dense bounding-box degree index for explicit domains whose full
addition table is too large brings the same support-domain diagnostic to 12.0
seconds (`z` support about 98 ms, contracted support about 8.7 seconds). The
actual all-ten rational GV reconstruction was still stopped after a few
minutes; stack samples remained in `compute_ckyz_inverse_mirror_map_domain ->
ckyz_series_exp_domain -> ckyz_series_mul_domain`. That confirms the next
builder needs cygv-style coefficient/path-history domains, not just narrower
support sets or faster monomial lookup.

A test-only z-series inversion prototype now mirrors cygv's residual
subtraction more directly: it builds the contracted instanton potential in
`z`, walks the finite downset in degree order, subtracts `Li2(q^d)` expressed in
`z`, and only then reads the target GV values. The local P2, F0, and polygon-5
tests match the existing flat-coordinate extractor when the full lower-degree
path history is included. A failed ray-only attempt was instructive: F0
`[1,1]` is wrong unless `[1,0]` and `[0,1]` are processed first, exactly the
same lower-degree history requirement seen in cygv's `series_inversion`.

After checking out `cygv` tag `v0.1.2` (the version installed in the local
CYTools Python environment), the production switch is source-aligned in one
important respect and still incomplete in another:

- Aligned: cygv keeps the instanton corrections in the original monomial
  variable, then performs degree-ordered residual subtraction with
  `Li2(q_N)`. Cyrus now uses the same residual shape for the
  support-predicted CKYZ API instead of first materializing the full
  flat-coordinate potential.
- Incomplete: cygv's traversal is over a finite semigroup sorted by a grading
  vector. Cyrus currently uses every nonzero element of the predicted monomial
  domain as the path history. This is correct for the tested local models, but
  it is still too broad for McAllister `[4,3,2]` at ten multiples.
- Performance hotspot: the long McAllister all-ten run was stopped after a few
  minutes with stack samples in z-series residual extraction, specifically
  repeated `q^d exp(d.alpha)` and `Li2` products over the 21k-element predicted
  domain. The earlier global inverse-map bottleneck has moved, not disappeared.

The latest source read also confirmed that CYTools' `CalabiYau.compute_gvs`
does not hide another Python-side GV algorithm. It builds CYTools inputs
(`curve_basis`, intersection numbers, Mori cap generators, and a grading
vector) and delegates to `cygv.compute_gv`; the `Invariants` basis option is
only an output charge re-expression. The implementation target remains cygv's
HKTY path, not an extra CYTools wrapper formula.

Cyrus now has a conservative z-residual dependency builder that walks backward
from the requested target degrees and keeps lower-degree subtraction history
when a `Li2(q^d)` contribution can affect a needed coefficient in `z`. The F0
regression `target=[1,1]` explicitly checks that `[1,0]` and `[0,1]` remain in
the path history. The support test
`ckyz_predicted_support_domain_api_matches_target_downset_for_polygon5_ray`
now reports:

```text
[CKYZ_Z_HISTORY] domain=265 terminals=2 selected=79 candidate_evaluations=447 exp_support_classes=6 elapsed=5.749584ms
```

For the McAllister-style `[4,3,2]` target through ten multiples, support
prediction still builds the same 21,721-monomial domain, but z-history pruning
now completes instead of stalling:

```text
[CKYZ_Z_HISTORY] domain=21721 terminals=10 selected=5235 candidate_evaluations=38195 exp_support_classes=6 elapsed=16.800214875s
```

The full N=10 filtered McAllister gate was stopped after about three minutes.
Sampling showed the remaining hot path in
`extract_ckyz_local_gv_invariants_from_z_potential_for_degrees ->
ckyz_q_degree_series_in_z_domain -> ckyz_series_exp_domain ->
ckyz_series_mul_domain`, with most leaf time in Malachite rational arithmetic.
This is now the source-aligned next target: mirror cygv's `InstantonData`
structure more closely by precomputing `exp(alpha_i)`/`exp(-alpha_i)` once and
building each `q_N` via monomial shifts and previous `q_N` cache entries,
instead of recomputing `exp(d.alpha)` from scratch for every history degree.

An attempted full-coordinate `exp(alpha_i)` power cache was correct on the
polygon-5 regression but too broad for McAllister: the N=10 run stalled while
computing `exp(alpha_i)` over the full 21,721-monomial predicted domain before
any residual subtraction began. Cyrus therefore now uses the other cygv idea in
production: when a previous nonzero `q_N` is available, it computes
`q_degree = q_previous * q_delta` for the smallest componentwise nonnegative
delta it can find, and caches those `q_delta` series. This keeps the cache
attached to actual subtraction history rather than eagerly materializing full
coordinate exponentials.

On the small polygon-5 trace this gives:

```text
[CKYZ_Z_EXTRACT] degrees=79 nonzero_gvs=43 previous_qns=43 delta_q_cache=4 elapsed=31.225459ms
```

The McAllister `[4,3,2]`, N=10 gate still did not complete. It reached the
same pruned history (`selected=5235`) and sampling showed the remaining hot
path inside the previous-`q_N` helper while computing a direct `q_delta`
exponential:

```text
ckyz_q_degree_series_with_previous_cache_in_z_domain
  -> ckyz_q_degree_series_in_z_domain
  -> ckyz_series_exp_domain
  -> ckyz_series_mul_domain
```

So the next target is not another whole-domain cache. The coefficient work must
be restricted to the degrees that can actually be read or subtracted in the
selected residual history, so that even the first `q_delta` exponential is not
computed over the full predicted monomial domain.

## May 6 Source Refresh

The current source-read pass rechecked the exact CYTools/cygv boundary before
adding more Rust logic:

- CYTools `CalabiYau.compute_gvs` is only an input builder. It computes
  intersection numbers, the no-origin curve-basis matrix, the Mori cap, optional
  lattice-point augmentation, and a grading vector, then calls
  `cygv.compute_gv`. There is no separate hidden Python GV algorithm to port.
- cygv `run_hkty` constructs the finite semigroup first, then computes
  `omega`, `alpha`, `beta`, the contracted instanton polynomials, and finally
  runs `series_inversion::invert_series`.
- cygv polynomial multiplication is domain-defined: if the sum of two monomial
  exponent vectors is absent from `PolynomialProperties::monomial_map`, the
  product is dropped. This is the real truncation contract.
- cygv `compute_instanton_data` precomputes `exp(alpha_i)` and
  `exp(-alpha_i)` once for each coordinate, but the final GV extraction still
  depends on degree-ordered residual subtraction. The rolling `previous_qn`
  cache is a performance aid for building `q_N`; it is not the correctness
  source.
- cygv's threefold inversion reads a candidate GV from the first nonzero curve
  coordinate, checks integrality, records nonzero invariants, computes
  `Li2(q_N)`, and subtracts `GV * component * Li2(q_N)` from all later
  instanton polynomials.

This confirms that the local CKYZ path should keep working in the z-domain and
update only the residual coefficients that are actually in the selected
history. Materializing whole `q_N` or whole `Li2(q_N)` series over the
21,721-monomial McAllister support-predicted domain is not source-required for
reading 5,235 history coefficients.

The current coefficient-level prototype follows that direction: it computes the
coefficient of `Li2(q^d exp(d.alpha(z)))` at a selected target degree directly,
using the standard exponential coefficient recurrence over the CKYZ finite
domain. The polygon-5 regression checks this coefficient path against the old
full-series `Li2` construction.

The first implementation-shape fix replaced the cloned vector-pair exponential
cache with source-equivalent indexed caches. The exponential coefficient cache
is now keyed by `(scale_id, delta_index)`, and repeated
`scale_degree dot alpha_term` values are cached per scale. The small polygon-5
trace improved monotonically:

```text
vector-pair coefficient cache: [CKYZ_Z_EXTRACT] ... exp_coeff_cache=1382 elapsed=35.5ms
indexed coefficient cache:    [CKYZ_Z_EXTRACT] ... exp_coeff_cache=1382 elapsed=29.4ms
per-scale alpha cache:        [CKYZ_Z_EXTRACT] ... exp_coeff_cache=1382 scaled_alpha_cache=1362 elapsed=26.1ms
alpha predecessor cache:      [CKYZ_Z_EXTRACT] ... exp_coeff_cache=1382 scaled_alpha_cache=1362 predecessor_deltas=213 elapsed=20.7ms
```

This did not make the McAllister `[4,3,2]`, N=10 gate finish. The gate still
reaches the same source-predicted domain and history:

```text
[CKYZ_SUPPORT_DOMAIN] rank=3 broad=26691 selected=21721 ...
[CKYZ_Z_HISTORY] domain=21721 terminals=10 selected=5235 candidate_evaluations=38195 ...
```

The alpha predecessor cache removes repeated source-degree scans inside the
exponential coefficient recurrence. It is a source-preserving memoization of the
same polynomial multiplication condition cygv uses: a term contributes only when
the predecessor and remainder monomials both lie in the finite domain.

Sampling after the cache changes no longer shows cloned vector-pair keys as the
top cost. The remaining hot path is recursive coefficient lookup over too many
selected history targets: per-scale `HashMap<usize>` lookups, cached coefficient
returns and clones, and Malachite rational multiply/GCD. The sampled N=10 run was
stopped after about 2.5 minutes without reaching the final extraction trace. The
next substantive change should therefore reduce the coefficient-history/domain
itself. More cache reshaping may help constants, but it will not change the core
5,235-history-degree workload.

A coefficient-work profiler now measures that workload without computing or
asserting GV values. It uses the same support-predicted domain and selected
z-residual history, then counts the formal
`Li2(q^d exp(d.alpha(z)))` coefficient states before doing any Malachite
rational arithmetic. The first-principles McAllister audit command:

```sh
env CYRUS_FIRST_PRINCIPLES=1 \
    CYRUS_CKYZ_MULTIPLES_TO_CHECK=10 \
    CYRUS_CKYZ_TARGET_DIRECTION=4,3,2 \
    CYRUS_PRINT_CKYZ_DOMAIN_PROFILES=1 \
    CYRUS_PRINT_CKYZ_COEFFICIENT_WORK=1 \
    CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 \
    cargo test -p cyrus-core --test potent_ray_affine_circuits \
    mcallister_rank_two_ckyz_domain_profiles_are_inventoried -- --nocapture
```

originally finished in 112.7 seconds and reported:

```text
[CKYZ_PROFILE] kind=3 direction=[4, 3, 2] multiples=10 target_downset=26691 predicted=21721 causal=Some(129766)
[CKYZ_COEFFICIENT_WORK] kind=3 direction=[4, 3, 2] multiples=10 domain=21721 history=5235 residual_pairs=13699995 componentwise_pairs=7154291 li2_terms=9042668 unique_scales=7065 unique_deltas=16750 unique_exp_states=7166981
```

This explains why a whole-domain exponential cache is the wrong next move:
even after z-history pruning, the current path asks for about 9.0 million
multiple-cover delta terms and 7.16 million distinct exponential coefficient
states. The next source-aligned optimization should build a coefficient-demand
graph for the selected residual history, so Cyrus computes only the exponential
coefficients actually read by the subtraction step rather than recursively
probing the full predicted domain for each `(scale, delta)` state.

The profiler now also counts source-level Li2 support pruning. Re-running the
same McAllister audit finished in 139.2 seconds and reported:

```text
[CKYZ_PROFILE] kind=3 direction=[4, 3, 2] multiples=10 target_downset=26691 predicted=21721 causal=Some(129766)
[CKYZ_COEFFICIENT_WORK] kind=3 direction=[4, 3, 2] multiples=10 domain=21721 history=5235 residual_pairs=13699995 componentwise_pairs=7154291 li2_terms=9042668 support_pairs=2620565 support_li2_terms=3219472 unique_scales=7065 unique_deltas=16750 unique_exp_states=7166981 support_unique_exp_states=2587631
```

That support profile is actionable: it reduces direct coefficient states from
7.16 million to 2.59 million, so Cyrus now applies the same support gate before
rational coefficient extraction. On the filtered McAllister `[4,3,2]`, N=4
validation this changes extraction from about 1.26 seconds in the rejected
batched-demand experiment to:

```text
[CKYZ_Z_EXTRACT] degrees=434 nonzero_gvs=217 li2_coefficients=9443 li2_support_skips=24090 li2_support_classes=6 exp_coeff_cache=23529 scaled_alpha_cache=23529 predecessor_deltas=1061 elapsed=817.675958ms
```

The coefficient routine also receives the support set now, so surviving target
pairs no longer compute unsupported multiple-cover deltas just to discover zero
coefficients recursively. On the small polygon-5 trace this lowers the
exponential cache from 1,382 entries to 1,165 entries while keeping the same
502 nonzero Li2 coefficients:

```text
[CKYZ_Z_EXTRACT] degrees=79 nonzero_gvs=43 li2_coefficients=502 li2_support_skips=815 li2_support_classes=6 exp_coeff_cache=1165 scaled_alpha_cache=1165 predecessor_deltas=174 elapsed=20.084541ms
```

The extraction loop now obtains the support flag and coefficient in the same
multiple-cover scan, rather than scanning once for support and then again for
the coefficient. This keeps the same cache sizes and gives:

```text
polygon-5: [CKYZ_Z_EXTRACT] ... exp_coeff_cache=1165 ... elapsed=19.270875ms
McAllister [4,3,2], N=4: [CKYZ_Z_EXTRACT] ... exp_coeff_cache=23529 ... elapsed=814.916084ms
```

The `[4,3,2]`, N=10 filtered gate still timed out after 240 seconds without an
extraction trace, even after reaching the same 5,235-degree history in 16.7
seconds. Support pruning is therefore a real reduction, but not the final
coefficient-demand solution. The remaining target is to shrink the selected
history or compute the surviving 2.59 million supported states with a
source-aligned dynamic program instead of recursive per-coefficient lookup.

A scale-distribution profile shows that the supported work is not concentrated
in a few huge scale degrees:

```text
[CKYZ_COEFFICIENT_WORK] ... support_unique_exp_states=2587631 support_scales=5235
[CKYZ_COEFFICIENT_SCALE] ... scale=[0, 1, 0] support_unique_exp_states=2519
[CKYZ_COEFFICIENT_SCALE] ... scale=[0, 2, 0] support_unique_exp_states=2457
[CKYZ_COEFFICIENT_SCALE] ... scale=[0, 0, 1] support_unique_exp_states=2419
[CKYZ_COEFFICIENT_SCALE] ... scale=[1, 1, 0] support_unique_exp_states=2419
[CKYZ_COEFFICIENT_SCALE] ... scale=[1, 0, 0] support_unique_exp_states=2399
```

The top scale has only about 2.5k supported deltas, while there are 5,235
support scales total. That makes an eager dense cache per scale unattractive.
The next coefficient evaluator should stay sparse and source-ordered, or the
history builder needs a sharper dependency criterion.

After switching CKYZ extraction and the coefficient profile to the same
cover-weight grading order, the McAllister `[4,3,2]`, N=10 profile finished in
154.6 seconds and reported:

```text
[CKYZ_COEFFICIENT_WORK] kind=3 direction=[4, 3, 2] multiples=10 domain=21721 history=5235 residual_pairs=13475754 same_grading_skips=224241 componentwise_pairs=7154291 li2_terms=9042668 support_pairs=2620565 support_li2_terms=3219472 unique_scales=7065 unique_deltas=16750 unique_exp_states=7166981 support_unique_exp_states=2587631 support_scales=5235
```

The grading-level batch order removes 224,241 same-grading residual pairs from
the profile, but it does not change the 2.59 million supported coefficient
states. This confirms that the next real speedup still needs the indexed
finite-polynomial evaluator described below, not just better pair ordering.

## Deeper cygv Source Pass

The latest no-code source pass re-read the exact `cygv-0.1.2` internals that
matter for the remaining CKYZ bottleneck:

- `run_hkty` constructs one finite `Semigroup`, then creates
  `PolynomialProperties` from that semigroup before any period coefficient is
  computed.
- `Semigroup::with_max_degree` trims supplied elements by grading degree,
  derives a reduced generator set, starts from the trimmed supplied elements,
  and closes under generator addition up to the maximum degree.
- `Semigroup::with_min_elements` derives the generator set, starts from the
  identity only, and increases the grading cutoff until the finite set is large
  enough.
- `Polynomial::mul`, `pow`, `exp_pos_neg`, and `li_2` all use the same
  `PolynomialProperties::monomial_map`; products whose exponent sum is absent
  from that map are discarded. This is the truncation contract.
- `compute_instanton_data` computes `alpha`, `beta`, `beta - alpha alpha`, the
  contracted instanton polynomials, and `exp(alpha_i)` / `exp(-alpha_i)` once.
- `series_inversion::invert_series` batches candidate classes by grading
  degree, reads all candidates at that grading before subtracting their
  `Li2(q_N)` contributions, then rolls the `previous_qn` cache forward.

This clarifies the remaining Cyrus divergence. The current CKYZ z-residual path
is algebraically source-derived, but it is not mechanically shaped like cygv:

- `CkyzMonomialDomain` stores an explicit finite monomial set and drops products
  outside it, which matches cygv's truncation idea.
- The support-predicted domain is not a cygv semigroup closure. It is generated
  from a broad target downset by support propagation through the CKYZ local
  operations. This is why it can be much smaller than the unit-generator causal
  semigroup, but also why it must be validated as a coefficient-demand domain,
  not described as the literal cygv domain.
- Extraction is currently coefficient-recursive: for each supported
  `(scale, delta)` it computes one coefficient of `exp(scale.alpha)` by a
  recurrence. cygv instead stores finite polynomials and builds `q_N` from
  cached `exp(alpha_i)` powers and recent `q_N` history.
- Cyrus now sorts CKYZ extraction by the cover-weight grading vector and reads
  a full grading level before applying its subtractions to later gradings. The
  remaining divergence is representation: residuals and coefficients still need
  to move further toward cygv's indexed finite-polynomial mechanics.

The next implementation should therefore be framed as a cygv-shaped indexed
series evaluator, not as another broad-domain tweak:

1. Give the CKYZ finite domain a source-style grading vector, grouped degree
   levels, and index-addressed coefficient arrays. The z-residual state is now
   domain-indexed; the broader finite-polynomial operations still need the same
   treatment.
2. Compute `alpha`, `beta`, local `beta - alpha alpha`, and the contracted
   instanton polynomial in that indexed domain.
3. Precompute sparse `exp(alpha_i)` / `exp(-alpha_i)` over the retained domain,
   or over a proven coefficient-demand subdomain with the same multiplication
   semantics.
4. Run a local analogue of `series_inversion::invert_series`: read all
   candidates at a grading level, build `q_N` from cached exponentials and
   recent history, subtract `Li2(q_N)` from later residual coefficients, then
   advance to the next level.
5. Keep `potent_rays_gv.dat` as an assertion-only comparison after the local
   support, source relations, cover weights, and target direction have all been
   reconstructed from geometry/source formulas.

Until that structure exists, more small cache changes can improve N=4 timing
but are unlikely to make the McAllister `[4,3,2]`, N=10 row robust. The source
gap is now specific: Cyrus needs cygv's finite-polynomial, grading-level
inversion mechanics in the local CKYZ setting, with the monomial set narrowed by
proved coefficient demand rather than by saved GV data.

## Indexed CKYZ Series Checkpoint

Cyrus now has a `CkyzIndexedSeries` layer for the finite CKYZ domain. It stores
nonzero coefficients as sorted `(monomial_index, coefficient)` pairs and uses
the existing domain addition table / dense degree index for multiplication. The
domain multiplication and exponential routines now convert from the public
`BTreeMap<Vec<usize>, Rational>` representation once, perform the internal
polynomial operation in indexed form, then convert back.

This is not yet the full cygv inversion structure, but it moves the broad period
and exponential operations onto the same representation cygv uses internally:
finite monomial indices plus a shared monomial map. The regression
`ckyz_indexed_domain_series_ops_match_total_degree_reference` checks the indexed
domain product and exponential against the older total-degree reference algebra
on a full total-degree domain.

Validation after this change:

```text
cargo fmt --check
cargo test -p cyrus-core ckyz_indexed_domain_series_ops_match_total_degree_reference -- --nocapture
cargo test -p cyrus-core ckyz_z_series_inversion_matches -- --nocapture
cargo test -p cyrus-core ckyz_coefficient_li2_matches_full_series_on_polygon5 -- --nocapture
cargo test -p cyrus-core ckyz_local_gv -- --nocapture
CYRUS_FIRST_PRINCIPLES=1 CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 CYRUS_CKYZ_TARGET_DIRECTION=4,3,2 CYRUS_CKYZ_MULTIPLES_TO_CHECK=4 CYRUS_TRACE_CKYZ_Z_HISTORY=1 cargo test -p cyrus-core --test potent_ray_affine_circuits mcallister_rank_two_ckyz_potent_ray_gvs_are_reconstructed -- --nocapture
```

The McAllister `[4,3,2]`, N=4 trace still reports the same selected history and
coefficient counts, with z extraction at about 817 ms:

```text
[CKYZ_Z_HISTORY] domain=1641 terminals=4 selected=434 candidate_evaluations=2842 exp_support_classes=6 elapsed=132.780375ms
[CKYZ_Z_EXTRACT] degrees=434 nonzero_gvs=217 li2_coefficients=9443 li2_support_skips=22508 li2_support_classes=6 exp_coeff_cache=23529 scaled_alpha_cache=23529 predecessor_deltas=1061 elapsed=817.10025ms
```
