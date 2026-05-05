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
