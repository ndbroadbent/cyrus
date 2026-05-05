# CYTools / cygv GV Audit

This note records the current read-through of the CYTools `compute_gvs()` path and the
`cygv` HKTY implementation. It is intended to prevent more diagnostic churn around
the McAllister Stage 5 GV residual.

## CYTools Contract

The authoritative CYTools wrapper is
`reference/cytools/src/cytools/calabiyau.py::_compute_gvs_gws`.

For hypersurface threefolds it constructs:

- `kappa = cy.intersection_numbers(in_basis=True, format="coo"/"dok")`
- `q = cy.curve_basis(include_origin=False, as_matrix=True)`
- `mori = cy.mori_cone_cap(in_basis=True)`
- `generators = vstack([mori.rays(), mori.find_lattice_points(min_points=100*h11)])`
- `grading_vector = mori.find_grading_vector()`, unless explicitly supplied

Then it calls:

```python
cygv.compute_gv(
    generators=generators,
    grading_vector=grading_vector,
    q=q,
    intnums=intersection_numbers_dok,
    max_deg=max_deg,
    min_points=min_points,
)
```

If `mcap_generators` is supplied, CYTools passes those rows directly as
`generators`; otherwise it augments Mori-cap rays with lattice points before
calling `cygv`.

One subtlety: CYTools does this augmentation with
`mori.find_lattice_points(min_points=100*h11)` regardless of whether the later
`cygv.compute_gv` call is controlled by `max_deg` or by `min_points`. The
`max_deg` value is a cygv semigroup truncation parameter, not a CYTools
lattice-enumeration cutoff.

The `mori_cone_cap` source is also part of the contract, not a black box. It
constructs ambient curve relations from two families:

- adjacent 2-simplices inside each 2-face of the primal polytope;
- origin circuits built from a 2-simplex plus one point from each of the two
  facets containing that simplex.

Each relation is normalized by the integer nullspace vector, sign-flipped so
that the first coefficient is nonnegative, and then projected either by dropping
the origin entry or by applying the divisor basis. This is the exact source of
the "toric curve" candidates used by the high-dimensional Kähler-side
selection.

CYTools' `Cone.find_lattice_points` is a separate operation. It first checks
that the grading vector has only the origin at nonpositive degree, then uses
CP-SAT to enumerate integer points in degree windows sorted by degree and
lexicographic order. This lattice enumeration is what the generic
`compute_gvs()` wrapper uses to augment Mori-cap rays; it is not the same as
the McAllister paper's selected-small-toric-curve cutoff.

The matrix-orientation convention is easy to get wrong:

- Python `cygv.compute_gv` receives `generators` as rows of curve coordinates and
  `q = cy.curve_basis(include_origin=False, as_matrix=True)` as an `h11 x
  n_divisors` row matrix.
- The Rust/PyO3 extension then stores each Python outer row as a column. Thus
  internally, `generators` becomes `h11 x n_generators` and `q` becomes
  `n_divisors x h11`.
- Cyrus' `cygv_q_matrix` mirrors that boundary by accepting the same `h11 x
  n_divisors` curve-basis matrix and transposing it before calling the Rust
  `cygv` API directly.

So the correct Cyrus input shape for `q_matrix` is "basis curve rows", not
"ambient divisor rows". A mismatch here would still type-check and can produce
plausible-looking but wrong GV data.

The `q` rows themselves come from the CYTools basis machinery:

- `CalabiYau.glsm_charge_matrix(include_origin=True)` uses the origin plus the
  prime toric divisors selected by the triangulation.
- A vector divisor basis is sorted, checked as an integral GLSM column basis,
  and stored with an identity block on the chosen basis columns.
- The dual curve-basis matrix is then built from HNF-normalized GLSM linear
  relations: rows are curve-basis elements, chosen basis columns form the
  identity, and non-basis columns are solved in reverse order using the last
  nonzero pivot in the corresponding linear relation.

Cyrus' `compute_curve_basis_matrix` mirrors this vector-basis path. For `cygv`,
therefore, the object to pass is the CYTools-style curve-basis matrix with
origin removed, and the direct Rust call must transpose it to the crate's
internal `n_divisors x h11` convention.

## cygv Semantics

The `cygv` crate does not interpret `generators` as merely the list of curves to
evaluate. They define the affine semigroup used by the HKTY series computation.

The Rust call chain is:

1. `Semigroup::with_max_degree`, `with_min_elements`, or `from_data`
2. `fundamental_period::compute_omega`
3. `instanton::compute_instanton_data`
4. `series_inversion::invert_series`

Important details:

- `Semigroup::with_min_elements` uses the supplied rows only to derive a
  generator set, starts the semigroup from the identity, and closes under
  addition until at least `min_points` elements exist.
- `Semigroup::with_max_degree` first trims supplied rows by grading degree,
  uses the trimmed rows both as initial elements and to derive generators, then
  closes under addition up to `max_deg`.
- Before closure, `cygv` calls `find_generators`, which removes some supplied
  elements that are sums of two others. The current crate implementation only
  tests sums with `n = 2` because `max_sum_elements = 3` and the Rust range is
  `2..max_sum_elements`.
- `compute_omega` builds `q0` from GLSM charges by summing rows for the
  hypersurface nef partition, then groups semigroup elements by the number of
  negative GLSM intersections. Only the 0-, 1-, and 2-negative-intersection
  groups are used by the current HKTY coefficient code.
- `invert_series` processes classes by grading degree and subtracts previously
  found lower-degree `Li2(qN)` contributions from the instanton corrections.

The derivative layer inside `fundamental_period.rs` is also explicit enough to
port locally when the compact hypersurface numerator is absent:

- `compute_c_0neg` builds the regular fundamental coefficient and its first and
  second `rho`-derivatives from factorial products, harmonic numbers, and
  order-two harmonic numbers;
- `compute_c_1neg` handles the simple-pole residue when exactly one GLSM
  pairing is negative; the first derivative is the residue times the negative
  row entry, and the second derivative multiplies that residue by the two
  regular harmonic sums;
- `compute_c_2neg` handles the double-pole case for second derivatives only,
  using the product of the two negative row entries and the parity sign from
  their pole orders.

Cyrus' CKYZ local-period helpers intentionally stop at this analogous
coefficient layer. They do not perform `compute_instanton_data`'s multiplication
by `c0_inv`, do not build alpha/beta/F polynomials, and do not run
`series_inversion`. A nonzero double-log coefficient in this layer is therefore
not yet a GV invariant.

The CKYZ flat-coordinate helper extends that local path by solving the formal
inverse mirror map `z_i(q) = q_i exp(-alpha_i(z(q)))` and substituting the
double-log/prepotential-period `z`-series into the resulting `q`-series. This
matches the mirror-map composition step but still stops before the
multiple-cover subtraction that cygv performs in `series_inversion`.

The local CKYZ GV helper now ports the relevant intermediate cygv conversion:
it forms `beta - alpha_i alpha_j`, contracts it with the local intersection
expression using the same diagonal `1/2` convention as cygv's symmetric
intersection contraction, and then applies the CKYZ `instbase` cover relation
`A_m = sum_{k d = m} w(d) N_d / k^2`. This is deliberately separate from the
compact cygv semigroup inversion and is currently validated only for the source
normalizations that are pinned by local `P^2`, `F0`, and `F1` tables.

This means local face/ray HKTY checks can be correct for isolated classes while
still not reproducing the global CYTools/cygv output. The global result depends
on the full semigroup truncation and the degree-ordered subtraction history.

It also means `cygv` should not be treated as the implementation of the paper's
"remove curves that are sums of others" step. Its generator reduction is an
internal heuristic for the supplied semigroup, not a documented Hilbert-basis
extraction pass over the selected toric curves.

The direct `cygv` source read also fixes the exact computation boundary. In
`cygv-0.1.2/src/hkty.rs`, `run_hkty` constructs the semigroup first, then builds
the fundamental period, instanton data, and finally inverts the series. In
`src/semigroup.rs`, `with_max_degree` and `with_min_elements` close the supplied
generator set under addition before any HKTY coefficient is computed. In
`src/fundamental_period.rs`, `compute_omega` multiplies the semigroup elements
by the GLSM charge matrix and dispatches only the 0-, 1-, and
2-negative-intersection cases. In `src/series_inversion.rs`, the final GV value
for a class is not just read off independently: lower-degree `Li2(qN)` terms are
subtracted from the remaining instanton polynomials degree by degree.

Therefore, when a saved McAllister ray has a GV sequence in
`potent_rays_gv.dat`, the hidden input is not just the ray vector. It is the
low-dimensional face semigroup in which that ray was evaluated, including the
other elements that affect the degree-ordered subtraction history. Feeding the
single ray and its multiples to `cygv` is a valid diagnostic, but it asks a
different mathematical question.

The corresponding Cyrus boundary is now explicit. A low-dimensional face
context may only replace the full Mori semigroup after an integer supporting
normal certifies that every retained face generator has zero pairing and every
ambient Mori generator has nonnegative pairing. `gv.rs` therefore exposes two
separate operations:

- `check_supporting_mori_face_normal` verifies a caller-supplied integer normal;
- `certify_supporting_mori_face_by_exact_kernel` proves codimension-one faces
  directly from the integer kernel of the supplied face generators, then runs
  the same verifier.

The exact-kernel path is deliberately conservative. If the face has a
higher-dimensional kernel, Cyrus returns no certificate instead of choosing an
LP-rounded normal. That leaves the harder higher-codimension face-selection
problem open, but prevents an approximate diagnostic from becoming a hidden
assumption in the GV pipeline.

The McAllister corrected-chamber LP-witness diagnostic confirms this boundary.
With exact-kernel and aggregate-normal LP certificate attempts enabled, the
current local face diagnostic computes values for several toric-missing curves,
but certifies none of them as supporting Mori faces. Those values remain
unpromoted because the semigroup context is not proven to be a valid face.

A separate source-derived correction is now promoted: the origin-circuit
resolved-conifold detector accepts the sparse standard charge pattern
`(-1,-1,1,1)`. This is not an LP witness. It is the ordinary isolated
resolved-conifold normal bundle `O(-1) + O(-1)` case, exposed in McAllister
diagnostics as an origin coefficient `-1`, one other negative unit coefficient,
and two positive unit coefficients. After this change, the corrected-chamber
diagnostic reports the checkpoint-t target correction as `toric_covered=412`
and `toric_missing=8`, while the current solved-t diagnostic reports
`toric_covered=410` and `toric_missing=9`. The remaining missing curves are
still nonstandard origin-circuit patterns and must stay explicit until their
actual semigroup/local-toric context is reconstructed from source.

A follow-up affine-support diagnostic makes this obstruction more precise. The
nine current solved-t corrected-chamber misses are all Mori generators and all
are origin circuits, but their affine supports are not rank-two local toric
surface diagrams: the rank inventory is `{3: 4, 4: 5}`. Each has a one-row
local affine charge basis because the sampled support has `points = affine_rank
+ 2`, but no rank-two local coordinate model is available. Therefore the
rank-two CKYZ/potent-ray machinery cannot simply be reused for these missing
nilpotent origin circuits. The next source-derived object has to be the
higher-rank local semigroup or flop/chamber context for these origin-circuit
supports.

The missing-target diagnostic now also records the B-field parity and real-axis
polylog branch status for these unresolved classes. For the current solved-t
4-214-647 corrected chamber, all nine missing targets are `real_ok`, all have
`q.t >= 0.1`, and the parity split is seven even classes and two odd classes.
So these specific misses are not caused by the `Li2/Li3` real branch cut. The
open implementation target is their higher-rank origin-circuit GV source or
supporting semigroup context, not a hidden branch-continuation choice for a
near-wall curve.

## McAllister Paper Boundary

The source paper itself gives a useful constraint on what may count as a
first-principles reproduction.  The corrected string-frame volume and Kähler
coordinates are defined by the classical intersection form, BBHL term, divisor
Euler corrections, and genus-zero GV sums.  The iterative KKLT equation then
solves the quadratic classical target with the GV `Li2` correction evaluated at
the previous Kähler point.

The paper's high-dimensional GV strategy is not "assign a local formula to
every small ray." It distinguishes potent rays, nilpotent curves, and toric
curves inherited from the ambient toric variety.  The selected toric curves are
chosen by a volume cutoff and pruned when they are sums of other selected
curves; the paper motivates this because in small-`h11` checks such toric
curves contain the Hilbert-basis-sized contributions.  That is evidence for a
selection algorithm, not a closed GV formula for arbitrary origin-circuit
coefficient patterns.

The paper also constrains corrected-chamber continuation.  If a shrinking curve
has zero `B_2` period, the real `Li2/Li3` expressions hit logarithmic branch
cuts and the simple string-tree formula is not valid through the transition. If
the `B_2` period is `1/2`, the `-e^{-2*pi*t}` polylog branch can be continued
and the result should match the flopped phase with transformed geometric data.
Therefore a corrected-chamber implementation needs explicit branch/flop
semantics; silently evaluating a real continuation or choosing a toric GV value
by coefficient pattern would be another hidden assumption.

## Computational Mirror Symmetry Source Read

The later paper *Computational Mirror Symmetry* (arXiv:2303.00757) is the
closest written source for the `cygv` implementation. It confirms the same
contract visible in the code:

- the master formula determines genus-zero GV invariants order by order from
  the fundamental period, mirror map, and intersection form;
- a finite truncation is consistent only when the selected set of curve classes
  is closed under the relevant "past" additions, called the causality
  constraint in the paper;
- the two practical truncations are the degree method and the past-light-cone
  method;
- at high `h11`, full-dimensional degree truncations are exponentially large,
  so useful computations are often done inside low-dimensional faces of the
  Mori cone;
- a one-dimensional ray computation is meaningful only when the ray is actually
  an effective Mori generator, with integrality/nontriviality used as evidence
  for that status;
- after candidate one-dimensional generators are found, two-dimensional cones
  spanned by pairs can be tested the same way, again using integrality as a
  consistency signal.

This is exactly the conceptual boundary Cyrus needs for potent rays. A saved
ray in `potent_rays.dat` is not enough. A first-principles implementation must
reconstruct the low-dimensional face/ray semigroup context, compute the HKTY
data in that context, and use `potent_rays_gv.dat` only as an assertion.

The paper also has a separate "GV invariants of toric curves" section. Those
formulas explain the already-ported selected-small-curve checks: simple
complete-intersection toric curves in two-faces have genus-zero GV values
computed from their curve moduli spaces, with the local cases controlled by the
normal bundle degrees. This direct toric-curve path is not the potent-ray
infinite-family path.

## Two Distinct GV Regimes

The McAllister pipeline uses GV data in two different regimes:

1. Mirror/racetrack GV data on the low-dimensional mirror side. For 4-214-647
   this is the `dual_curves.dat` / `dual_curves_gv.dat` data, and the Python
   reproduction notes validate it with `cy.compute_gvs(min_points=20000)`.
2. Kähler-coordinate corrections on the high-`h11` primal side. This is the
   `small_curves.dat` / `small_curves_gv.dat` data. The paper explicitly says
   systematic high-degree GV computation at large `h11` is infeasible, and that
   the method starts from toric ambient curves, removes decomposable curves, and
   computes selected toric/nilpotent curves plus low-dimensional face/ray data.

The current Stage 5 residual is in the second regime. Treating it as a plain
full-cone `cy.compute_gvs()` mismatch is therefore misleading.

Older Python reproduction docs and scripts sometimes describe Step 8 as simply
`cy.compute_gvs(max_deg=N)` or `cy.compute_gvs(min_points=N)`. That is accurate
for the low-dimensional mirror/racetrack validation data (`dual_curves*.dat`)
and useful as a CYTools/cygv smoke test, but it is not the production recipe
for the high-dimensional Kähler-side selected curves and potent rays. The
McAllister paper and the computational-mirror paper both separate those cases.

The paper's selected-curve method is more precise than several older Python
notes in `string_theory/mcallister_2107`:

- start from curves inherited from the toric ambient variety;
- cut by volume in the approximate KKLT solution, with homogeneous scaling
  corresponding to `c_tau = 1`;
- remove curves that are sums of other selected curves, because those are not
  Hilbert-basis elements;
- compute GV values of the remainder, and separately probe potent rays in
  low-dimensional Mori faces for convergence control.

This is the method behind `small_curves.dat` and `small_curves_gv.dat`. The
generic `cy.compute_gvs(min_points=...)` scripts in the Python reproduction are
useful for the low-dimensional mirror/racetrack side and for experiments, but
they are not a faithful production recipe for the 214-dimensional Kähler-side
instanton correction.

## Current Cyrus Situation

Cyrus already has wrappers around the same Rust `cygv` crate:

- `compute_gv_invariants`
- `compute_gv_invariants_with_degree_bounded_lattice`
- `compute_gv_invariants_with_provided_generators`
- `compute_gv_invariants_with_explicit_semigroup`
- `potent_ray_convergence` for the paper's low-dimensional-face convergence
  diagnostic once a ray and its `N_{nq}` values are available.
- `compute_one_dimensional_ray_gv_series` for extracting `N_q, N_{2q}, ...`
  along a supplied Kähler-basis ray via the cygv provided-generator path.

Their boundaries are different:

- `compute_gv_invariants` mirrors the CYTools default augmentation: it always
  augments Mori-cap rays with `min_points=100*h11` lattice points before calling
  `cygv`, even when the final semigroup truncation is controlled by `max_deg`.
- `compute_gv_invariants_with_degree_bounded_lattice` is the explicitly named
  diagnostic shortcut: when `max_deg` is supplied it bounds the Mori-cone
  lattice-point enumeration before calling `cygv`. This is useful for controlled
  investigations, but it is not the literal CYTools wrapper contract.
- `compute_gv_invariants_with_provided_generators` mirrors CYTools'
  `mcap_generators=...` entry point: Cyrus supplies rows directly to cygv and
  lets cygv perform its own semigroup closure.
- `compute_gv_invariants_with_explicit_semigroup` bypasses cygv closure with
  `Semigroup::from_data`. This is a diagnostic-only shape for hand-constructed
  local semigroups, not a replacement for CYTools' public `compute_gvs()`.
- `compute_ray_gv_series_with_provided_generators` extracts
  `N_q, N_2q, ...` for a target ray from a caller-supplied local
  face/semigroup generator context. The older one-generator helper is now just
  this function with `[q]` as the supplied context, so the API boundary no
  longer conflates "ray series extraction" with "the one-generator semigroup is
  the correct physics context."
- `diagnose_affine_toric_circuit` recognizes when an ambient curve row is an
  affine toric circuit and currently classifies the local `P^2` triangle
  pattern without assigning any GV values. This is an upstream-geometry
  diagnostic for reconstructing potent-ray face contexts, not a checkpoint
  replay mechanism.
- `rank_two_local_charge_model` reconstructs the affine local charge lattice
  of a normalized rank-two support signature and records the integer
  coordinates of the target potent-ray relation in that lattice. This identifies
  the one-parameter direction to feed into a future source-derived local
  mirror/HKTY calculation without using the saved GV row as input.

But the McAllister corrected-chamber diagnostic path currently relies first on
toric two-face/origin-circuit formula values for selected small curves, then uses
local one-dimensional ray or LP-witness face HKTY diagnostics only for misses.
That is not equivalent to CYTools' full `cy.compute_gvs(min_points=N)` call.

A pass over the current Cyrus call sites separates the GV code into three
different trust levels:

- Production-computed for the current 4-214 checkpoint: the selected small
  toric-curve classes and their `1/-2` GV values, via
  `compute_toric_two_face_curve_gv_invariants` and
  `compute_toric_curve_gv_diagnostics`.
- Validation-only: `potent_rays.dat`, `potent_rays_gv.dat`,
  `corrected_kahler_param.dat`, `corrected_target_volumes.dat`, and
  `corrected_cy_vol.dat`. These files may be loaded by tests or optional
  comparison paths, but they are not production inputs for GA evaluation.
- Diagnostic-only: `compute_one_dimensional_ray_gv_series`,
  `compute_ray_gv_series_with_provided_generators`,
  `compute_gv_invariants_with_explicit_semigroup`, and the
  `mcallister_first_principles` LP-witness/lattice-saturation retries. These
  are useful for isolating possible low-dimensional face contexts, but they are
  not promoted to production unless the supporting face and semigroup context
  are independently reconstructed and certified.

One additional audit caveat: the main runner paths that call `cygv` remove the
origin/canonical column from `compute_curve_basis_matrix` before building the
`q` matrix, matching `CYTools`' `curve_basis(include_origin=False)`. The ignored
`stage5_first_potent_ray_one_generator_gv_diagnostic` currently builds its
diagnostic `q` matrix from the full curve-basis rows. That makes its printed
one-generator sequence even less suitable as physics evidence; it should only
be read as a guard that the saved potent-ray row has not been regenerated.

For potent-ray convergence, Cyrus can now compute the exact span rank, corrected
Kähler volumes, and `log xi_n = log|N_{nq}| - 2π n q.t` decay slopes for a
supplied ray/GV sample. This makes `potent_rays*.dat` checkable rather than
opaque, but it is still not the full paper method: Cyrus does not yet sample
low-dimensional faces of `M_infinity(X)` or compute the `N_{nq}` sequence along
those rays from cygv.

The ray series computation has now been promoted out of the
`mcallister_first_principles` diagnostic into `cyrus-core`. That gives the GA
pipeline a reusable way to ask for `N_{nq}` once it has a candidate ray and the
local/global cygv generator context. It is not yet validated against the saved
`potent_rays_gv.dat` values for 4-214-647, because the upstream face/ray
sampling and exact local input construction still need to be ported.

For the high-`h11` Kähler-correction side, however, exact equivalence to a full
global `compute_gvs()` call is not the McAllister method either. The relevant
comparison is the selected toric-curve/face method described in the paper and
represented by `small_curves.dat` and `small_curves_gv.dat`.

The latest Stage 5 residual is therefore best interpreted as a selected-curve
coverage/chamber/continuation problem, not as a dilog aggregation problem and
not as a problem in the already-checked local toric formulas.

## Ancillary Data Semantics

The paper-data `readme.txt` is explicit about the high-dimensional Kähler-side
files for `4-214-647`:

- `small_curves.dat` contains curves found using toric information below a
  cutoff in the approximate KKLT vacuum.
- `small_curves_cutoff.dat` is `1.0`; the cut is applied before instanton
  corrections.
- `small_curves_gv.dat` contains GV invariants for those selected curves.
- `small_curves_vols.dat` contains their volumes in the KKLT vacuum after the
  instanton-corrected move; negative entries mark curves flopped when moving
  from `kahler_param.dat` to `corrected_kahler_param.dat`.

Direct inspection of the `4-214-647` files gives:

- `small_curves.dat`: `344 x 219`
- `small_curves_gv.dat`: 344 values, with 315 entries equal to `1` and 29
  entries equal to `-2`
- `small_curves_vols.dat`: exactly
  `small_curves[:, basis.dat] * corrected_kahler_param.dat`
- approximate-chamber volumes
  `small_curves[:, basis.dat] * kahler_param.dat` are all positive and below
  the cutoff
- corrected-chamber stored volumes include 10 negative/flopped curves

So the selection chamber and the evaluation chamber are not the same thing. The
paper selected small toric curves in the approximate chamber, then evaluated the
instanton-corrected solution where some of those selected classes have crossed a
flop wall.

## Python Reproduction Boundary

The older Python reproduction contains two different styles of GV use:

- the mirror/racetrack side computes or validates `dual_curves.dat` using
  CYTools `compute_gvs(min_points=20000)`;
- the primal KKLT correction side often accepts a `gv_invariants` dictionary as
  an upstream input, and some latest-CYTools experiments call
  `cy.compute_gvs(min_points=100)` generically.

That generic `compute_gvs` call is not the paper's high-`h11` selected
small-curve method. It can be useful as a small validation experiment, but it is
not an implementation strategy for the 214-dimensional Kähler correction.

Some older notes also describe `small_curves.dat` as if it were obtained by a
plain `cy.compute_gvs(max_deg=N)` call. The paper and ancillary readme are more
specific, and they take precedence: for the large-`h11` Kähler side, the method
starts from toric ambient curves, cuts by approximate-vacuum volume, removes
curves that are sums of others, and then computes the GV values of the retained
small toric/nilpotent curves. The direct `cy.compute_gvs(min_points=20000)`
validation applies cleanly to the low-dimensional mirror/racetrack
`dual_curves` data, not to the high-dimensional `small_curves` data.

Concretely, `string_theory/mcallister_2107/REPRODUCTION_OUTLINE.md` still says
Step 8 is `cy.compute_gvs(max_deg=N)` and points to `small_curves.dat` as its
verification file. Treat that line as a stale Python-era simplification. It is
not the contract Cyrus should port for the high-dimensional Kähler correction.
The better statement is: `dual_curves*.dat` validate generic low-dimensional
`compute_gvs`, while `small_curves*.dat` validate the selected toric-curve
pipeline.

Two source details make this boundary explicit:

- `latest_cytools/compute_gv_invariants.py` loads `dual_curves.dat` and
  `dual_curves_gv.dat`, constructs the dual polytope triangulation, runs
  `cy.compute_gvs(min_points=100)`, and then maps basis coordinates back to
  ambient coordinates with `cy.curve_basis(include_origin=True, as_matrix=True)`.
  That script is a low-dimensional mirror-side check, not the high-`h11`
  selected-small-curve algorithm.
- `mcallister_2107/CLAUDE.md` warns that `min_points`, not `max_deg`, is the
  practical CYTools/cygv control for those generic GV runs and states that
  McAllister's quoted "max degree" was the highest sampled curve, not a direct
  cutoff. This matches the CYTools source: `max_deg` belongs to cygv semigroup
  truncation, while CYTools' wrapper independently augments generators with
  `mori.find_lattice_points(min_points=100*h11)`.

The no-code read-through also confirms that
`string_theory/mcallister_2107/full_pipeline_from_data.py` is still a
checkpoint-driven scaffold, not a GA-ready computation. It loads
`dual_curves.dat` and `dual_curves_gv.dat` as upstream GV data, projects those
saved ambient curves to the current basis, and reads `cy_vol.dat` for the final
volume with an explicit TODO to derive it. That script is useful for checking
the racetrack algebra against McAllister's files, but Cyrus must replace those
loaded objects with:

- dual/mirror-side `compute_gvs` from the CYTools/cygv source contract;
- primal high-`h11` selected toric/nilpotent curve construction;
- potent-ray GV samples recomputed in their low-dimensional face semigroup
  contexts;
- the corrected Kahler volume from the KKLT fixed-point and chamber-continuation
  data, not from `cy_vol.dat`.

## Paper Contract for Large h11 GVs

The McAllister paper separates the large-`h11` Kähler-side GV problem into
potent and nilpotent regimes:

- Potent rays are rays with infinitely many nonzero GV invariants. These govern
  the radius of convergence, so the paper checks them by computing GV invariants
  to high degree in low-dimensional faces of the Mori cone.
- Nilpotent curves have only finitely many nonzero multiples. The important
  small nilpotent classes are found by starting from toric ambient curves,
  cutting by Kähler volume, removing curves that are sums of others because they
  are not Hilbert-basis elements, and computing the GV values of the remainder.
- The paper states that, in examples where systematic low-volume checks are
  possible, the nonzero low-volume curves found this way are the toric curves
  with order-one GV values.

For 4-214-647, the ancillary `potent_rays*.dat` files are therefore not
optional decoration: they are the saved evidence for the potent-ray convergence
check. The `small_curves*.dat` files are the saved evidence for the selected
nilpotent/small-toric contribution to the Kähler-coordinate instanton sums.

## Source-Derived Port Map

The next GV work should be split by the source algorithm, not by the names of
the ancillary files. The paper, ancillary readme, CYTools wrapper, and `cygv`
source imply the following Cyrus contracts.

### Selected Nilpotent / Small-Toric Curves

This is the contract behind `small_curves*.dat`:

1. Build the primal triangulation from `points.dat` and `heights.dat`.
2. Compute ambient Mori-cap curve relations with the CYTools
   `mori_cone_cap()` circuit construction.
3. Select toric curves whose volume is positive and below
   `small_curves_cutoff.dat` at the approximate Kähler point
   `kahler_param.dat`. The paper describes this as homogeneous scaling with
   `c_tau = 1`; the checkpoint value is the resulting validation point for
   4-214-647, not a general production input.
4. Remove selected curves that are sums of other selected curves. The
   4-214-647 checkpoint matches pair-decomposable pruning exactly, while
   Cyrus' stricter finite-semigroup diagnostic removes five additional curves.
   Both policies are now explicit, so a GA production run must choose one
   deliberately rather than inheriting the checkpoint accident silently.
5. Compute GV values for the retained toric/nilpotent curves. For the
   checkpoint-covered classes, Cyrus' local two-face and resolved-conifold
   origin-circuit formulas reproduce all saved `1/-2` values. For retained
   classes not covered by these local formulas in a new chamber or geometry,
   the next production method must be a certified face/local HKTY computation
   or a faithful general CYTools/cygv computation, not a fitted value.

Executable evidence for this contract is
`stage5_mcallister_small_toric_curves_match_checkpoint`,
`stage5_mcallister_small_toric_curve_finite_semigroup_diagnostic`, and
`stage5_mcallister_small_toric_curve_gvs_match_checkpoint`.

### Potent-Ray / `M_infinity(X)` Checks

This is the contract behind `potent_rays*.dat`:

1. Construct or sample low-dimensional faces of `M_infinity(X)` in phases where
   those faces are also faces of the Mori cone.
2. Produce rational rays inside those low-dimensional faces and track their
   ambient curve classes.
3. For each ray, compute `N_q, N_2q, ..., N_10q` by running HKTY in the
   corresponding low-dimensional face/semigroup context. A one-generator
   provided-generator call is only a diagnostic unless it is supplied with the
   same local semigroup context used by the face computation.
4. Validate the sample by recomputing row-span rank, corrected-Kähler volumes,
   and the paper's `log xi_n = log|N_{nq}| - 2 pi n q.t` slopes.

Cyrus currently implements only step 4 plus a reusable
`compute_one_dimensional_ray_gv_series` helper. It does not yet implement steps
1-3. Therefore `potent_rays.dat` and `potent_rays_gv.dat` remain validation
samples, not reusable pipeline inputs.

The first useful production test here is not another final-volume assertion.
It should be a focused test that takes one saved potent ray, reconstructs the
low-dimensional face/semigroup context from upstream geometry, regenerates the
first ten GV values, and compares them to the corresponding row of
`potent_rays_gv.dat`. Only after that passes should Cyrus attempt to generate
the full 411-ray 4-214-647 sample.

A direct one-generator diagnostic has now been run for the first saved
4-214-647 potent ray. The ray is not one of the projected Mori-cap generators,
and `compute_ambient_one_dimensional_ray_gv_series(..., max_multiple=10)`
returns `[4, -11, 60, -478, 4588, -49368, 575896, -7131274, 92429484,
-1241983287]`, while the saved checkpoint row is `[3, -6, 27, -192, 1695,
-17064, 188454, -2228160, 27748899, -360012150]`. Therefore the saved
potent-ray GV sequence is not the one-generator `cygv` semigroup; Cyrus must
reconstruct the local face/semigroup context.

The 4-214-647 checkpoint also gives a more specific hint about that missing
context. The ancillary readme says `potent_rays.dat` uses the same `h11+5`
ambient coordinate convention as `dual_curves.dat`: intersections with the
canonical and prime toric divisors. In those coordinates, the first saved
potent ray has sparse support
`[(43, 1), (155, 1), (168, -3), (169, 1)]`, and its first GV values are
`[3, -6, 27, -192, 1695, ...]`. That is the standard genus-zero local
`P^2` (`O(-3) -> P^2`) sequence, and the sparse relation has the expected
`(1, 1, -3, 1)` local charge pattern up to ordering.

The lattice points themselves confirm that this is an affine circuit, not just
a numerical coincidence. In triangulation-point coordinates:

```text
43  = ( 0,  1, -3, 6)
155 = (-2, -1, -4, 5)
168 = (-1,  0, -3, 5)
169 = (-1,  0, -2, 4)
```

and `p43 + p155 - 3*p168 + p169 = 0`, with coefficient sum zero. The same four
indices do not form adjacent two-face simplices in either Cyrus' default FRST
snapshot or the McAllister `heights.dat` triangulation snapshot: none of
`(43,168,169)`, `(155,168,169)`, `(43,155,168)`, or `(43,155,169)` appears as a
triangulation triple. So the local `P^2` context is not recovered by simply
looking for the current-phase two-face circuit in `heights.dat`; it is part of
the missing low-dimensional face/phase reconstruction described by the paper.
Equivalently, `p168 = (p43 + p155 + p169)/3`, so this is the local picture of an
interior point of a toric triangle rather than an ordinary four-vertex
quadrilateral flip.

This is not a proof that every saved potent ray is local `P^2`, but it is a
strong source-derived clue about the right implementation shape. The 411 saved
potent rays collapse to only 24 unique GV sequences; 337 of the rays have 5
nonzero ambient entries and 56 have 4. A direct coefficient/point audit shows
that all 411 saved rows are affine relations in the triangulation-point
configuration, with 395 rank-2 supports and 16 rank-4 supports. The most common
rank-2 patterns are local toric-surface-like sparse charge vectors, e.g.
`(-3,1,1,1)` for local `P^2`, `(-5,1,1,1,2)`, `(-8,1,1,3,3)`, and
`(-10,2,2,3,3)`. Cyrus now records affine support rank and, for the first local
`P^2` ray, reconstructs rank-two local coordinates
`p43=(1,0)`, `p155=(0,1)`, `p168=(0,0)`, `p169=(-1,-1)` directly from the
ambient lattice points. Cyrus now implements the first local `P^2` instance as
an exact one-parameter local mirror calculation: it expands the Picard-Fuchs
mirror map for `O(-3) -> P^2`, transforms the local B-model Yukawa coupling to
the flat coordinate, and applies the multiple-cover inversion to recover
`N_q, ..., N_10q`. The gated 4-214-647 test reads `potent_rays_gv.dat` only as
the assertion for row 0.

This is still not a proof that all potent-ray rows are solved. The remaining
work is to reconstruct the local toric charge/face contexts for the other
rank-two patterns and to handle the 16 rank-four affine supports without
collapsing them into a local `P^2` special case.

The source-read checkpoint after inspecting CYTools and the installed
`cygv-0.1.2` Rust crate is:

- CYTools does not contain the GV algorithm itself; it normalizes the geometry
  into `mcap_generators`, a grading vector, the CYTools curve-basis matrix, and
  intersection numbers, then delegates to cygv.
- cygv's actual boundary is HKTY on a semigroup. The sparse charge vector of a
  potent ray is only one relation in that local toric support; it is not, by
  itself, the semigroup context used by the degree-ordered subtraction in
  `series_inversion`.
- Cyrus now reconstructs the integer affine charge lattice of each saved potent
  ray support from upstream point coordinates, by taking the integer kernel of
  the matrix with rows `[1; lattice coordinates]`.
- Cyrus also reconstructs integral two-dimensional local coordinates for every
  rank-two potent-ray support, by computing a Hermite-normal-form basis for the
  affine difference lattice. In the 4-214-647 checkpoint this covers all 395
  rank-two supports and deliberately leaves the 16 rank-four supports without
  rank-two coordinates.
- Cyrus now derives a normalized rank-two support signature from those local
  coordinates. The signature removes point-index labels, normalizes relation
  sign, and canonicalizes over integral coordinate systems generated by
  support-point differences. This is for grouping local toric supports before
  deriving local mirror/HKTY input; it is not a GV-value lookup.
- A gated 4-214-647 test verifies that every saved potent-ray sparse relation
  lies in the rational row span of that reconstructed local charge basis. This
  uses `potent_rays.dat` as an assertion only; the charge basis itself is
  computed from `points.dat` and the filtered triangulation-point set.
- The same gated test verifies that the rank-two affine relations remain zero
  after projection to the reconstructed local coordinates.
- The structural potent-ray test verifies that all 56 recognized local `P^2`
  supports collapse to one normalized signature while non-`P^2` rank-two
  supports remain distinguishable.

This is the right intermediate object for the next implementation step. It is
not a coefficient-pattern lookup table, and it is not yet a GV computation for
the non-`P^2` families. The missing step is to convert each reconstructed local
charge basis plus local support into the exact local toric/HKTY input expected
by cygv: local semigroup generators, grading vector, local `q` matrix, and the
intersection/Yukawa data relevant for the local model.

### Corrected-Chamber Continuation

The selected-small-curve checkpoint is an input-chamber construction. The
corrected Kähler point can cross flop walls; for 4-214-647, ten saved
small-curve volumes become negative after moving to
`corrected_kahler_param.dat`.

For a GA-ready corrected-volume solve, Cyrus must choose and test one of these
source-backed paths:

1. Recompute selected toric/nilpotent curves, local GV values, divisor Euler
   characteristics, and the Kähler-coordinate instanton vector in the current
   corrected chamber at each KKLT fixed-point iteration.
2. Or implement the explicit flop/chamber continuation of the original
   input-chamber instanton data, including the accompanying transformations of
   the intersection form and divisor Euler characteristics.

Simply evaluating the saved input-chamber finite curve list at negative
`q.t`, dropping near-flop curves, flipping signs, or fitting per-curve weights
has already been ruled out by diagnostics and is not a production method.

### Related Infinity-Cone / NOP Algorithm

The later McAllister moduli-space reconstruction paper
(`mcallister_moduli_space_reconstruction_2212.10573.tex`) gives useful context
for the corrected-chamber problem. It is not the source of the 4-214-647
ancillary ray sample, but it makes the geometry behind the missing step more
explicit:

- `M_infinity` is the closure of the cone generated by potent rays, and its
  dual is the hyperextended Kähler cone.
- A flop curve is nilpotent and lies strictly outside `M_infinity`; the paper
  calls these nilpotent-outside-potent, or `nop`, curves.
- A nilpotent generator on or inside `M_infinity` is not an ordinary flop curve;
  it corresponds to an unstable Weyl reflection or an infinite-distance/CFT
  boundary depending on the facet physics.
- The finite-cutoff algorithm starts from nonzero GV charges up to a grading
  cutoff, classifies apparently nilpotent primitive rays by checking whether a
  later multiple has zero GV while earlier weighted GV data is nonzero, then
  tests whether those rays diverge from the potent cone using degree slices,
  affine lattices, LLL-reduced bases, and infinity-norm distances at half and
  full cutoff.
- Once nop curves are identified, crossing a chamber wall is not just a sign
  change in a finite instanton list. The flop/stable-Weyl continuation modifies
  the classical geometric data and the GV data of the new chamber.

This reinforces the production boundary above: the corrected-chamber fix should
be a chamber/`M_infinity`/nop continuation implementation or a current-chamber
recomputation, not a local adjustment to saved input-chamber curve weights.

## Flop Continuation

The paper warns that continuing the Kähler-coordinate formula through a flop is
subtle. For curves with half-integral `B_2` field (`gamma dot q` odd), the
relevant real-axis dilogarithm continuation is:

```text
Li2(-exp(-2*pi*t))/(2*pi)^2
  = -Li2(-exp(-2*pi*(-t)))/(2*pi)^2 - 1/2*t^2 - 1/24
```

The local TeX source for arXiv:2107.09064 prints this equation with a
`+ 1/2*t^2` term. Direct evaluation of the real dilogarithm identity gives the
negative sign above:

```text
Li2(-x) + Li2(-1/x) = -pi^2/6 - 1/2 log(x)^2, x > 0
```

This is not a place to blindly copy the printed sign into code. Cyrus'
`real_dilog_real_axis` tests match the standard identity and mpmath/SciPy
reference values in the near-flop regime.

This polynomial correction is tied to the transformations of `chi(D_i)` and
`kappa_ijk` across the flop. Therefore simply evaluating the original chamber's
finite GV list at negative `q.t` is not a first-principles corrected-chamber
implementation unless the associated chamber/topology transformation is handled
explicitly.

Executable checks now pin the boundary:

- the real odd-parity negative branch satisfies the standard flop dilogarithm
  identity used by Cyrus' `real_dilog_real_axis`;
- among the 10 negative `small_curves.dat` classes at
  `corrected_kahler_param.dat`, 8 have odd `gamma dot q` and can be evaluated
  on the real `-exp` branch, while 2 have even parity and lie on the real
  `Li2(exp(...))` branch cut.

So the saved input-chamber small-curve set cannot be repaired by applying the
simple odd-parity flop identity uniformly to all negative entries.

## Verified Checkpoints

The approximate-chamber selected-small-curve path is no longer speculative. With
`CYRUS_FIRST_PRINCIPLES=1`, the focused Stage 5 tests recompute from upstream
geometry:

- 419 raw sub-cutoff toric curve candidates from `points.dat`, `heights.dat`,
  `basis.dat`, `kahler_param.dat`, and `small_curves_cutoff.dat`;
- 344 retained curves after the currently verified pair-decomposable pruning;
- 5 of those 344 pair-pruned checkpoint curves are still finite-semigroup
  decomposable as multi-term sums of other selected raw candidates, so a full
  finite-semigroup pruning of the selected input-chamber set would retain 339
  curves and would no longer match `small_curves.dat`;
- exact equality with `small_curves.dat`, using that file only as an assertion;
- exact equality with all 344 `small_curves_gv.dat` values using Cyrus'
  toric two-face/origin-circuit local formulas, again using the file only as an
  assertion.

The successful command was:

```text
CYRUS_FIRST_PRINCIPLES=1 cargo test -p cyrus-core --test mcallister_e2e \
  stage5_mcallister_small_toric -- --nocapture
```

This narrows the remaining Stage 5 issue: the input-chamber small-curve
selection and its local `1/-2` toric GV values are reproduced. The unresolved
part is how to get a production-correct corrected-chamber instanton correction:
full semigroup/Hilbert reduction, missing non-toric face/ray contributions, or
explicit flop/chamber continuation.

## Next Productive Step

Stop expanding top-contributor local checks. The selected-small-curve
approximate-chamber checkpoint is now reproduced, so the next useful comparisons
are corrected-chamber and generality audits:

1. Decide which selected-curve pruning rule to use for GA production search.
   Cyrus now exposes this as `CurvePruningStrategy` and the
   `mcallister_first_principles --small-curve-pruning <pair|finite-semigroup>`
   flag. The default `pair` reproduces the 4-214 checkpoint; the stricter
   `finite-semigroup` path intentionally removes 5 more input-chamber curves.
2. Apply the selected-small-curve rule to the corrected chamber, or if the
   intended McAllister operation is analytic continuation through flops, port
   that transformation explicitly using the paper's polylog identities.
3. Port the low-dimensional face/ray GV path for potent and nilpotent checks,
   where `cygv` full-semigroup input dumps are actually the right method.
4. Use generic CYTools/cygv `compute_gvs` only for the low-dimensional
   mirror/racetrack side and for isolated face/ray subproblems.

The open question is now sharper: is the corrected-target residual caused by a
missing corrected-chamber/non-toric contribution, an incomplete decomposition
reduction, or by applying input-chamber small curves without the required
flop/chamber continuation?

Older wording of this audit treated the faithful selected-small-curve method as
the next code change:

1. Reconstruct, from code/paper, the exact rule that produced `small_curves.dat`
   and `small_curves_gv.dat`: chamber, volume cutoff, toric-curve enumeration,
   pair-decomposable pruning, and face/ray GV method.
2. Compare Cyrus against that rule in the approximate chamber first, where the
   ancillary files are direct checkpoints.
3. Then apply the same rule to the corrected chamber or, if the intended
   McAllister operation is analytic continuation through flops, port that
   transformation explicitly using the paper's polylog identities.
4. Use `cygv` full-semigroup input dumps only for low-dimensional face/ray
   subproblems where that is the method actually being used.

Items 1 and 2 in that older list are now done for the approximate chamber; do
not reopen that loop unless a new source contradicts the current checkpoint.

Concrete missing implementation pieces after the source read are:

- a policy decision for GA runs: keep checkpoint-faithful pair pruning for
  comparability, or run the implemented finite-semigroup selected-set pruning
  knowing it intentionally differs from `small_curves.dat`;
- the low-dimensional face/ray GV computation path used for large-`h11` potent
  and nilpotent curve checks, distinct from a full 214-dimensional global HKTY
  run;
- an explicit chamber/flop continuation rule for the Kähler-coordinate instanton
  terms, including what to do with even-parity curves that hit the
  `Li2(exp(...))` branch cut.

## May 2026 Source-Read Checkpoint

The latest no-code pass through CYTools and cygv reinforces that the next
implementation should start from source-derived local input construction, not
from another Rust-side GV shortcut.

CYTools' `compute_gvs()` path has only two real responsibilities before it
delegates to cygv:

- build the CYTools curve-basis matrix with `curve_basis(include_origin=False)`;
- build the Mori-cap generator context, either from explicit
  `mcap_generators` or from `mori_cone_cap().rays()` plus
  `mori.find_lattice_points(min_points=100*h11)`.

The important source details are in those inputs. `mori_cone_cap()` constructs
ambient relations from adjacent 2-simplices in primal 2-faces and from
"origin circuits" built from a 2-simplex plus one point from each containing
facet. `Cone.find_lattice_points()` first checks that the proposed grading has
only the origin at nonpositive degree, then CP-SAT enumerates integer points in
degree windows and returns them sorted by degree and lexicographic order. These
choices define the semigroup that HKTY sees; they are not cosmetic
preprocessing.

The cygv crate then treats the supplied rows as a semigroup seed, not as an
independent list of classes. `run_hkty` constructs a `Semigroup`, computes the
fundamental period, computes instanton data, and only then performs the
degree-ordered `Li2(qN)` subtraction in `series_inversion`. The high-level cygv
entry points still `unwrap()` several internal errors, while the lower-level
functions return `Result`; any production Cyrus path that builds custom local
semigroups should keep using the lower-level calls or otherwise convert those
failures explicitly.

The Python reproduction scripts are mixed-trust sources. In particular,
`string_theory/mcallister_2107/full_pipeline_from_data.py` loads
`dual_curves*.dat` and final volume checkpoints, while
`2021_cytools/full_pipeline.py` uses generic `compute_gvs(min_points=20000)`
for the low-dimensional mirror/racetrack side. Those are useful validation
scaffolds, but they are not the high-`h11` selected-small-curve or potent-ray
algorithm. The stale line in `REPRODUCTION_OUTLINE.md` saying high-dimensional
GVs are just `cy.compute_gvs(max_deg=N)` should not guide Cyrus implementation.

The concrete next source-derived object is therefore:

1. inventory the normalized rank-two support signatures that Cyrus now derives;
2. for each signature family, derive the local toric charge/mirror input from
   the reconstructed local coordinates and local charge lattice;
3. compute the local mirror/HKTY series from that input;
4. compare against `potent_rays_gv.dat` only after the local input has been
   derived without using the saved GV values.

This is intentionally different from a coefficient-pattern table. The
coefficient patterns are a diagnostic index into local toric geometry, not the
source of the GV values.

The first item is now pinned down for the `4-214-647` checkpoint: the 395
rank-two saved potent-ray supports collapse to 16 normalized support
signatures. The gated structural test records the coefficient-pattern counts
for that inventory, but still does not use those patterns to assign any GV
value. The next implementation step remains item 2: derive the local
toric/mirror input for a non-`P^2` signature family.

Cyrus now also derives the first part of item 2 for every normalized rank-two
signature: a canonical local charge model. This model consists of the canonical
support points, the target affine relation in that point order, and the integer
kernel of the local `[1; x; y]` matrix. The gated potent-ray tests verify that
the target relation lies in this reconstructed local charge lattice for every
rank-two saved support. This is still pre-GV input construction; the remaining
missing piece is the local mirror/HKTY series attached to those charge models.

The current source-read boundary is therefore deliberately conservative. Cyrus'
local `P^2` routine is a one-parameter local mirror-map calculation for
`O(-3) -> P^2`; it is not a template that can be generalized by changing the
coefficient vector. For a non-`P^2` rank-two support, the local charge model can
have more than one independent charge row, and the saved potent ray is an
integer direction inside that charge lattice. The next real implementation
should build the full local toric/HKTY input for one such support family first:
local semigroup generators, grading vector, local `q` matrix, and the relevant
intersection/Yukawa data. After HKTY has produced local GV data in that context,
Cyrus can restrict to the multiples of the target charge direction and compare
with `potent_rays_gv.dat`.

The gated test `mcallister_rank_two_local_charge_models_are_inventoried` now
pins the exact pre-GV inventory behind this statement: all 395 rank-two saved
supports collapse to 16 canonical charge models, each with a fixed target
direction in its local charge basis. This remains upstream input construction.
It deliberately does not attach GV sequences to the non-`P^2` models.

The companion test `mcallister_rank_two_local_supports_are_reflexive_polygons`
checks the geometry behind that inventory. After shifting each support so the
negative-coefficient compact point is at the origin, every normalized rank-two
family has exactly one interior lattice point and every polygon edge has lattice
distance one. The next local mirror implementation should therefore be organized
around canonical bundles over these reflexive polygon toric surfaces.

This keeps the ancillary GV rows as assertions. It also explains why simply
calling cygv on one ray, or adding a coefficient-pattern-to-GV dispatch table,
would be the wrong abstraction: both would skip the semigroup and lower-degree
subtraction history that cygv's `series_inversion` actually uses.

## May 2026 cygv Source Addendum

A fresh read of the installed `cygv-0.1.2` source keeps the same conclusion but
adds a few implementation details that should constrain the next port:

- `run_hkty` is a pipeline wrapper, not an error-safe public contract. It builds a
  `Semigroup`, calls `compute_omega`, calls `compute_instanton_data`, and calls
  `invert_series`, but it still unwraps each internal `Result`. Cyrus'
  production paths should keep using lower-level calls when custom local inputs
  are being assembled.
- `Semigroup::with_min_elements` does not keep the caller-supplied rows as the
  truncation domain. It reduces the supplied rows to generators, starts from the
  identity, and closes under addition until enough elements exist.
- `Semigroup::with_max_degree` first trims caller-supplied rows by grading degree,
  then uses both the trimmed rows and their reduced generator set before closure.
  That makes the caller's local generator context part of the mathematical input.
- `fundamental_period::compute_omega` infers the CY dimension from the `q` matrix
  shape and nef partition, and rejects `cy_dim < 3`. This is the concrete code
  reason compact cygv cannot be treated as a local toric-surface engine.
- `compute_omega` groups semigroup elements by the number of negative GLSM
  intersections and only computes coefficient formulas for the 0-, 1-, and
  2-negative cases. Classes outside those groups are not independent local
  outputs.
- Threefold `series_inversion` extracts a candidate GV for a class by choosing the
  first nonzero curve-coordinate component, dividing the corresponding instanton
  coefficient by that component, checking integrality, then subtracting
  `GV * component * Li2(q_N)` from every affected instanton polynomial.
- The subtraction uses a rolling cache of previously computed `q_N` polynomials
  whose depth depends on `h11`. This reinforces that the answer for a ray is a
  function of the surrounding semigroup history, not just the sparse ray vector.

For Cyrus, the actionable consequence is unchanged: a local non-`P^2` potent-ray
implementation needs source-derived local semigroup generators, a grading vector,
a local `q`/charge matrix, and the correct local Yukawa/intersection object before
any comparison to `potent_rays_gv.dat` is meaningful.

## May 2026 No-Code Source Read

This pass deliberately read the source again before adding more Rust code. The
important boundary is now sharper:

- CYTools' public `CalabiYau.compute_gvs()` path is only a thin orchestration
  layer. It prepares `intersection_numbers(in_basis=True)`,
  `curve_basis(include_origin=False, as_matrix=True)`, `mori_cone_cap(in_basis=True)`,
  a lattice-point augmentation of the Mori cap, and `mori.find_grading_vector()`;
  then it calls `cygv.compute_gv`.
- CYTools' `mori_cone_cap()` is not a generic cone-enumeration shortcut. It builds
  candidate curve relations from adjacent two-simplices inside primal two-faces
  and from origin circuits formed by a two-simplex plus one point from each
  containing facet. Those relations are normalized integer kernel vectors before
  any divisor-basis projection.
- CYTools' `Cone.find_lattice_points()` is a CP-SAT enumeration by grading degree.
  Before enumeration it checks that the grading vector leaves only the origin at
  nonpositive degree. Results are sorted by degree, then lexicographically inside
  each degree.
- The `cygv` Rust Python boundary stores each outer Python row as a nalgebra
  column. Thus the Python `generators` list is written as rows, while the Rust
  `Semigroup` sees columns. The same row-to-column conversion applies to `q`.
- The packaged `cygv` Python helper forwards `prec` positionally to the Rust
  extension slot where the current signature has `n_threads`; CYTools does not
  pass `prec`, so this is not on the current McAllister path, but precision-mode
  Python cygv output should not be treated as source evidence without checking
  this boundary first.
- Cyrus links the Rust `cygv` crate directly, so for local-HKTY work the relevant
  source is `cygv-0.1.2/src/{hkty,semigroup,fundamental_period,instanton,series_inversion}.rs`,
  not the unavailable default `python3` import in this shell.

The porting implication is that we should not expand Cyrus by guessing a local
surface GV formula from the saved rows. The next non-documentation work should
either port a specific CYTools pre-cygv input-construction step, or construct a
single CKYZ/local-HKTY example all the way from source relations to GV extraction
with `potent_rays_gv.dat` used only as the final assertion.

## May 2026 Targeted CKYZ GV Checkpoint

Cyrus now has a source-derived CKYZ local GV extraction path that mirrors the
cygv instanton-data layer for the CKYZ local surface examples:

- compute CKYZ log-period corrections from the source relation rows;
- solve the inverse mirror map formally;
- form the local `beta - alpha_i alpha_j` instanton potential;
- contract with the CKYZ local intersection expression, including the diagonal
  `1/2` factor;
- apply the `instbase` multiple-cover inversion with finite-limit cover weights.

The targeted CKYZ APIs truncate by requested source degrees instead of a single
total degree. The default McAllister rank-two gate now uses the
support-predicted domain API, which computes alpha/beta supports on the broad
source downset and performs inverse mirror-map/potential domain prediction at
the support level before evaluating the final rational series. That is enough
to validate the first four `potent_rays_gv.dat` entries for all 395 rank-two
CKYZ McAllister potent-ray rows without using those GV rows as inputs. Polygon 5
is included with finite-limit cover weights `[1, 1, 1]`; the printed
`C1 = 3J1 + 2J2 + 2J3` weights are intentionally rejected because they produce
non-integral invariants in this extraction. This is not yet enough for all ten
entries: the current multivariable series still becomes too slow on the largest
source directions, so the remaining work is direct coefficient-support
extraction rather than a new physics shortcut.

The gated McAllister CKYZ test defaults to the first four multiples. Setting
`CYRUS_CKYZ_MULTIPLES_TO_CHECK=N` raises that assertion count explicitly for
diagnostics. Before dense monomial indexing, guarded addition tables, and the
target-downset domain, unfiltered `N=4` and polygon-5 direction `[4,3,2]` were
the first clear performance blockers. With those source-preserving domain
improvements, the full rank-two `N=4` gate passed in about 110 seconds through
the broad target-downset path. After promotion to the support-predicted API,
the same full rank-two `N=4` gate passes in about 132 seconds, so that change is
a correctness/structure boundary rather than a performance breakthrough. A
later flattened addition table, followed by sparse valid-pair multiplication,
also made the narrowed polygon-5 `[4,3,2]` `N=4` gate complete in about 149
seconds. A focused polygon-5 `[4,3,2]`, `N=5` diagnostic through the
support-predicted API exceeded a 300 second timeout, and an earlier unfiltered
`N=5` diagnostic exceeded a 600 second timeout, so higher all-row/all-ten
validation remains a finite-domain construction problem rather than a known
coefficient mismatch. A follow-up diagnostic that only constructed the
support-predicted polygon-5 `[4,3,2]`, `N=5` domain also exceeded 300 seconds,
so the bottleneck is already inside support-domain prediction rather than only
the final rational potential extraction.
The same test can be narrowed with `CYRUS_CKYZ_TARGET_DIRECTION=a,b,...`.
Focused all-ten checks now pass for the F0 directions `[1,1]`/`[1,2]` and F1
directions `[2,1]`/`[3,1]`; polygon-5 direction `[4,3,2]` passes through `N=4`
when narrowed but remains the first slow family. The immediate
coefficient-domain work should therefore focus on the rank-three polygon-5
local model first.

## May 2026 CYTools/cygv Porting Gaps

The deeper source read identifies three concrete gaps that should guide the
next code, and none of them are fixed by broadening the existing Rust shortcut.

First, the CYTools basis contract has two paths. Cyrus mirrors the ordinary
vector divisor-basis path well enough for the 4-214 checkpoint:
`set_divisor_basis([indices])` builds the dual curve-basis matrix by putting an
identity block on the chosen GLSM columns and solving non-basis columns from
HNF-normalized linear relations. But CYTools also has an experimental matrix
basis path. In that path `curve_basis(as_matrix=True)` always returns the
stored row matrix, and `mori_cone_cap(in_basis=True)` projects ambient Mori-cap
rows with `mori_cap_matrix.dot(basis.T)` rather than by selecting basis
columns. Cyrus now has explicit helpers for that projection:
`project_ambient_curve_to_basis_matrix` and
`project_mori_cone_cap_rays_to_basis_matrix`. Cyrus also has a source-derived
matrix divisor-basis to dual curve-basis constructor,
`compute_curve_basis_matrix_from_divisor_basis_matrix`, for the CYTools matrix
branch. The public `DivisorBasis` enum and
`compute_curve_basis_matrix_for_divisor_basis` helper now make this
vector-vs-matrix dispatch explicit for callers.
`project_mori_cone_cap_rays_for_divisor_basis` and `gv_divisor_basis_data`
now bundle the corresponding Mori-ray projection, ambient dual curve basis,
and no-origin q-matrix construction for either basis shape. The
`curve_basis_matrix_without_origin_i64` and
`curve_basis_q_matrix_for_divisor_basis_i64` helpers also centralize the
`curve_basis(include_origin=False, as_matrix=True)` q-matrix boundary for direct
cygv calls, so the origin column is not hand-dropped at each GV call site. The
ordinary vector-basis GV paths in `mcallister_first_principles` now route their
Mori projection, curve-basis matrix, and q-matrix construction through
`gv_divisor_basis_data`, so the runner no longer repeats those pieces by hand
for the default mirror, primal fallback, branch-report degree summary, and
corrected-chamber diagnostics. The one-off `mcallister_gv` and
`mcallister_racetrack` binaries use the same bundled handoff.
`mcallister_first_principles --dual-basis` JSON override now distinguishes index
and matrix basis shapes. Index overrides continue through the current
vector-basis path, while matrix overrides parse and then fail with an explicit
unsupported-vector-path error instead of being ignored or misinterpreted as
indices. The remaining gap is higher-level basis typing: APIs that still only
accept vector basis indices need to either accept a matrix basis end-to-end or
reject it loudly at every GA-facing boundary. Silently treating every basis as a
vector of indices is not GA-ready. The core cygv input handoff now accepts both
basis shapes, but the current McAllister runner still stores its pipeline basis
as `Vec<usize>` and remains intentionally vector-basis only until its flux and
intersection-basis transforms are generalized.

Second, cygv's final GV values are history-dependent. In
`series_inversion.rs`, the threefold path chooses a candidate value from the
first nonzero curve-coordinate component, checks integrality, records it, and
then subtracts `GV * component * Li2(q_N)` from the remaining instanton
polynomials for all affected coordinates. The `q_N` construction reuses a
rolling cache of previous degrees. Therefore a saved ray row plus its first ten
multiples is not enough input unless the surrounding local semigroup has also
been reconstructed. The answer for a ray depends on the semigroup's lower
degree classes and on the order in which they are subtracted.

Third, the current all-ten CKYZ blocker is domain shape, not arithmetic alone.
Dense indexing and guarded addition tables fixed the obvious multiplication
overhead, and the first-four checks now run for all 395 rank-two CKYZ
rows. Cyrus' targeted extractor now uses the union of componentwise past
downsets of the requested degrees instead of the single box of coordinate
maxima, so incomparable target degrees no longer force irrelevant mixed
monomials. The all-ten ray-direction case still needs a sharper
coefficient/path-history domain, because multiples of one positive direction
have a large componentwise past even after this downset refinement. The next
implementation should build that coefficient-targeted local monomial domain for
one canonical rank-two support family, then generalize across the 16 normalized
support signatures. More box optimization would still be asking cygv/CKYZ the
wrong finite-domain question.

This makes the next source-derived implementation boundary precise:

1. wire the CYTools generic matrix-basis construction path through higher-level
   GA-facing APIs, or reject unsupported matrix bases loudly at APIs that still
   accept only index bases;
2. for one non-`P^2` rank-two support signature, construct the local toric/HKTY
   input from local coordinates and charge lattice without using saved GV rows;
3. refine the current target-downset domain into a coefficient/path-history
   domain that contains the target multiples and the lower-degree terms needed
   by series inversion;
4. compare to `potent_rays_gv.dat` only as an assertion after the local input
   and finite domain have been built from geometry.

## May 2026 Polynomial Domain Source Read

The next CKYZ/potent-ray work should be modeled on cygv's finite polynomial
domain, not on another broad componentwise box.

In `cygv-0.1.2`, `PolynomialProperties::new` builds `monomial_map` directly
from the columns of the closed `Semigroup`. All polynomial operations then use
that map as the truncation domain:

- `Polynomial::mul` adds the two semigroup column vectors and silently drops the
  product if the sum is absent from `monomial_map`.
- `recipr`, `pow`, `exp_pos_neg`, and `li_2` all iterate only to
  `semigroup.max_degree / min_degree`, with every intermediate product still
  routed through the same semigroup map.
- `series_inversion` processes distinct grading degrees in order, computes
  `q_N` and `Li2(q_N)` for the nonzero classes at that degree, and subtracts
  their contributions before advancing.

This is stricter than saying "keep every monomial with coordinates below the
target." The domain is an additive semigroup closed enough for the requested
history, with a grading order that determines when lower classes are removed.
The current Cyrus CKYZ `target_downset` is useful because it avoids unrelated
incomparable targets, but it is still a formal componentwise past domain. For
large ray multiples such as polygon-5 directions, that past is enormous even
though the actual cygv-style computation only needs the semigroup elements that
can participate in the target coefficient and its prior subtraction history.

The next implementation should therefore introduce an explicit local
`CkyzCausalDomain` or equivalent structure with:

1. local charge-coordinate elements generated from the reconstructed CKYZ
   relation rows and target ray direction;
2. an integer grading vector and deterministic degree ordering;
3. a monomial lookup/addition map where absent sums are intentionally outside
   the finite domain;
4. cover-closed target multiples plus all lower-degree classes whose
   `Li2(q_N)` subtraction can affect those targets.

That gives a concrete correctness test before performance tuning: for small
targets where the current `target_downset` is affordable, the causal domain
must reproduce the same extracted CKYZ GV values. Only after that equivalence
holds should the domain be used to raise `CYRUS_CKYZ_MULTIPLES_TO_CHECK` toward
the full ten-entry `potent_rays_gv.dat` rows.

The first implementation of this boundary is now present as
`compute_ckyz_local_gv_invariants_for_degrees_with_causal_domain`. It takes
explicit local generators and a grading vector, closes the generated semigroup up
to the largest cover-closed target grading, and then runs the same CKYZ
instanton-potential and multiple-cover extraction on that finite domain. The
small guardrail tests currently prove:

- local `P^2` causal-domain extraction reproduces the known degree-1 through
  degree-10 table;
- the standard-generator F0 ray causal-domain extraction agrees with the
  existing `target_downset` extractor for the first four diagonal targets;
- explicit CKYZ domains drop products whose summed degree is absent from the
  monomial map, matching the cygv polynomial-domain contract.

Cyrus now also exposes
`ckyz_local_surface_causal_domain_spec` and
`compute_ckyz_local_surface_gv_invariants_for_multiples_with_causal_domain` for
matched CKYZ local surfaces. These derive the coordinate-unit generators from
the matched CKYZ source relation basis and use the finite-limit cover weights as
the positive grading vector. This removes a caller-invented generator/grading
choice for guardrail computations and passes the F1 ray equivalence check
against the existing targeted extractor.

That source-derived causal domain is still not narrow enough to become the
default all-row McAllister gate. A full rank-two `N=4` run through the
source-weighted causal helper was stopped after several minutes, while the
targeted first-principles F0 directions `[1,1]`/`[1,2]` and F1 directions
`[2,1]`/`[3,1]` now have explicit all-ten regressions against saved rows and
the existing `N=4` support-predicted all-row gate remains the practical default.
Flattening the target-downset addition table
removes enough lookup overhead for the narrowed polygon-5 `[4,3,2]` `N=4` check
to complete. Iterating only sparse valid product pairs improves that narrowed
run to about 149 seconds, but the sample still spends almost all time in inverse
mirror-map composition and rational polynomial multiplication. The
support-predicted API moves domain construction ahead of full rational
composition, but the focused `[4,3,2]`, `N=5` run still timed out after 300
seconds, and the domain-only diagnostic for the same target also timed out
after 300 seconds. This confirms the next missing object more sharply: Cyrus
needs direct coefficient/path-history support generation, not just the full
semigroup generated by the CKYZ coordinate axes up to a source grading cutoff.

This is still not the full McAllister potent-ray solution. The F0
`[1,1]`/`[1,2]` and F1 `[2,1]`/`[3,1]` all-ten checks prove the CKYZ
source-derived path can reach complete saved rows for multiple non-`P^2`
families, but the missing step is still to derive the local generator set and
grading for each normalized support signature from geometry, then use a
narrower coefficient/history domain to raise the larger rank-two CKYZ families
beyond the first four multiples.

## May 2026 Corrected-Chamber Source Checkpoint

A fresh read of the exact local `cygv-0.1.2` crate source and the CYTools
`CalabiYau._compute_gvs_gws` wrapper confirms that the remaining
corrected-chamber GV misses should not be chased with another per-curve formula.

The CYTools wrapper supplies rows of Mori-cap generators, augments them with
`mori.find_lattice_points(min_points=100*h11)` unless `mcap_generators` is
explicitly supplied, and passes the no-origin curve-basis matrix as `q`.
Inside cygv, those rows are not just a list of requested output classes:
`Semigroup::{with_max_degree,with_min_elements}` derives generators, closes the
semigroup under addition, sorts by grading degree, and then the HKTY pipeline
builds the fundamental period, instanton data, and degree-ordered series
inversion. `series_inversion` subtracts lower-degree `Li2(q_N)` contributions
before later degrees are read. Therefore the GV assigned to a class is a
function of the finite semigroup/chamber context, not of the sparse relation
alone.

This matters for the nine current solved-t corrected-chamber misses. They are
Mori generators and origin circuits, but their affine supports have ranks
`{3: 4, 4: 5}`, not the rank-two CKYZ surface supports where Cyrus already has
source-derived local formulas. They are also real-axis evaluable at the current
point (`real_ok`, `q.t >= 0.1`), so the blocker is not the `Li2/Li3` real branch
case. The next useful diagnostic is exact semigroup context reconstruction:
record whether the LP active-generator witnesses are genuine integer-semigroup
decompositions or only rational-cone decompositions, then use that to identify
the smallest source-derived finite domain that can be passed through the HKTY
path. Any value produced from an uncertified LP face or a sparse coefficient
pattern must remain diagnostic, not a promoted GV value.

The last 12 hours of commits did make progress, but they also show the danger
of continuing with local diagnostics indefinitely. Completed work includes the
CYTools-style GV basis handoff, matrix/vector basis separation, CKYZ local
period and multiple-cover extraction, rank-two potent-ray checks through four
multiples by default, all-ten checks for selected F0/F1 families, exact Mori
face certificate machinery, sparse resolved-conifold origin-circuit coverage,
and corrected-chamber missing-target branch/affine-support diagnostics. The
unsolved core is narrower now, but still unsolved: reconstruct the missing
higher-rank origin-circuit semigroup/chamber source, then run HKTY on that
source-derived context.

The exact LP-witness coefficient check is now recorded as well. For the nine
current solved-t corrected-chamber misses, the floating LP active-generator
witnesses solve exactly over `Q`, but split into five integer-semigroup
decompositions and four rational-cone-only decompositions. Example integer
witnesses have coefficients such as `[1, 1, 1]`, `[2, 1]`, `[1, 2]`, and
`[3, 1]`; example rational-only witnesses include `[1/2, 3/2, 1/2]`,
`[1/2, 3/2]`, and `[1, 1/2, 3/2]`. This is useful source-derived triage, not a
GV assignment: even the integer cases still require the actual finite
semigroup/chamber history before HKTY output can be promoted, and the four
rational-only witnesses cannot be explained as decompositions in the current
integer semigroup generated by those active rows.

The fresh release diagnostic
`/tmp/cyrus_corrected_chamber_source_context.log` confirms that the current
diagnostics are not missing an obvious rank-two or CMS-general-divisor shortcut.
The corrected-chamber solved-t run reports `ambient_rays=561596`,
`subcutoff=561`, `pair_pruned=419`, `toric_covered=410`, and
`toric_missing=9`. Each missing target is an origin-circuit Mori generator; the
required grading degrees range from `10` to `26`, with `720..2963` generators at
or below the target degree. The sampled local charge bases include
`[1,-2,-1,-1,3]`, `[1,-2,3,-1,-1]`, `[1,-1,-1,-2,3]`, and
`[2,1,2,-1,-2,-2]`, so these are higher-rank charge models, not the normalized
rank-two CKYZ surface rows already handled elsewhere.

The same diagnostic runs the current CMS-general-divisor shape checks on the
negative support entries of the nine misses. Those checks produce candidate
formula values such as `1`, `-2`, and `3`, but the exact divisor-intersection
verification fails for every candidate (`has_rational_divisor_solution=false`).
So these candidate values must remain unpromoted. This is the clearest current
boundary: the next implementation should export or reconstruct the finite
origin-circuit semigroup/chamber context and then run the HKTY path on that
source-derived domain. Adding another coefficient-pattern-to-GV dispatch would
be indistinguishable from fitting the saved rows.

The old Python reproduction artifacts explain why this matters for the Rust
goal. `full_pipeline_from_data.py` still loads `dual_curves.dat`,
`dual_curves_gv.dat`, and finally `cy_vol.dat` for `V[0]`. The 2021 CYTools
pipeline improved the geometry side but still used checkpoint data for branch
selection (`t_expected` as the Newton seed) and model/control inputs such as
`target_volumes.dat`, `basis.dat`, and `kklt_basis.dat`. Those uses can be
acceptable as validation scaffolding, but they are not a GA-ready pipeline. In
Cyrus, these files should remain assertions or explicit model choices; they
cannot supply downstream GV values, corrected Kähler parameters, final volumes,
or branch seeds.

The corrected-chamber GV diagnostic now has a source-shaped JSON export via
`--dump-corrected-chamber-gv-context`. The export intentionally writes the
finite degree-bounded semigroup slice needed by the current missing targets,
not the full corrected-chamber Mori-cap generator set. For 4-214-647 solved-t,
the full set has `561596` projected rays, while the missing-target bound
`degree<=26` contains `2963` rays. The export also records the no-origin
`cygv` q-matrix, curve-basis matrix, grading vector, corrected-chamber
intersection entries, and complete diagnostics for all nine missing
origin-circuit targets. This is still diagnostic input for HKTY reconstruction,
not a promoted GV fallback.

The first consumer for that export is `mcallister_gv_context`. In dry-run mode
it validates the context shape and reconstructs the five integer-semigroup
active-generator decompositions as explicit decomposition diamonds with only
`6..8` elements. Running `--run-integer-diamonds` on the 4-214-647 context is
fast, but it confirms the diamonds are not sufficient source contexts: three
integer targets return `GV=0`, while two fail inside cygv series inversion with
non-integral GV output. The four rational-cone-only targets are skipped. This
rules out promoting the LP active-generator diamond as the missing GV fallback;
the HKTY domain still needs the broader source-derived finite semigroup/path
history.

The next wider diagnostic window is also negative. `mcallister_gv_context
--run-active-support-generators` filters the degree-bounded projected Mori rows
to those supported on the union of the missing target support and its LP active
generator support. This gives small `mcap_generators` windows of `4..13` rows
per missing target. On the 4-214-647 context, six windows return `GV=0`; the
other three panic inside the provided-generator cygv path with
`NonIntegerGVError`. This shows that support-local generator closure is still
too small or otherwise not the correct CYTools/cygv source domain. The context
report now records support-overlap and iterative support-closure counts as
well: one-hop overlap windows are already `48..227` rows, and four support
closure layers reach `366..2814` generators with support sizes up to the full
`214` coordinates. A naive support closure therefore explodes back toward the
large corrected-chamber domain.

One-hop support-overlap windows are now executable diagnostics too. With
`--run-support-overlap-generators 4`, the 4-214-647 missing targets use `3..19`
provided generators. This also fails as a source domain: seven targets panic in
cygv with `NonIntegerGVError`, and the remaining two return `GV=0`. So the
correct context is not just "rows sharing several active/target support
coordinates" either.

The compact dual-polytope CYTools-to-cygv handoff is now checked at the source
boundary for 4-214-647. A Cyrus dump from
`CYRUS_GV_DUMP_INPUTS=/tmp/cyrus_gv_handoff_4_214_647.json
--bin mcallister_gv -- --min-points 1` was compared with a direct CYTools latest
dump of `cy.compute_gvs` inputs from `dual_points.dat` and
`dual_simplices.dat`. The grading vector, no-origin `q` matrix, and
in-basis intersection dictionary match exactly. CYTools hands cygv `505`
generator rows, while Cyrus dumps `496`, but CYTools' rows have only `496`
unique vectors and the unique generator sets are identical. This means the
standard compact GV path is not currently blocked by a q-matrix, grading,
intersection, or Mori/lattice generator handoff mismatch. The remaining
corrected-chamber misses are a chamber/local source-domain problem, not an
ordinary compact-dual cygv input mismatch.

`mcallister_gv_context` can now measure cygv's own finite semigroup closure
size without running the full HKTY period/series computation, via
`--measure-cygv-semigroups`. This is deliberately guarded by
`--semigroup-measure-max-seeds` because the first corrected-chamber missing
targets at degree `10` already have `720` seed rows at or below the target
degree. An unguarded degree-10 measurement did not finish within `120s`; the
guarded run with `--semigroup-measure-max-target-degree 10
--semigroup-measure-max-seeds 700` reports both degree-10 targets as
`skipped_seed_limit` with `720` seeds. This rules out treating the remaining
origin-circuit misses as a tiny visible-generator problem: the source-shaped
object is cygv's generated semigroup closure, not only the exported Mori rows
or LP-active rows.

After rebuilding `mcallister_gv_context`, the dry-run report now preserves the
origin-circuit provenance exported by `mcallister_first_principles`. A fresh
first-principles run:

```text
CYRUS_FIRST_PRINCIPLES=1 \
CYRUS_MCALLISTER_DATA_DIR=/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/4-214-647 \
timeout 600 ./target/release/mcallister_first_principles \
  --stop-after volume --allow-invalid-ek0 \
  --diagnose-corrected-chamber-gv \
  --dump-corrected-chamber-gv-context /tmp/cyrus_corrected_chamber_gv_context_20260506_053155.json

./target/debug/mcallister_gv_context \
  --context /tmp/cyrus_corrected_chamber_gv_context_20260506_053155.json \
  --out /tmp/cyrus_corrected_chamber_gv_report_20260506_053155.json
```

confirms that all nine current missing targets are higher-rank origin circuits,
not rank-two CKYZ local-surface supports. Their affine ranks and local charge
rows are:

| target | degree | exact LP kind | affine rank | local charge row | cygv hypersurface `cy_dim` |
|--------|--------|---------------|-------------|------------------|----------------------------|
| 0 | 18 | integer semigroup | 4 | `[2, 1, 2, -1, -2, -2]` | 4 |
| 1 | 18 | integer semigroup | 4 | `[1, 2, 1, -2, -1, -1]` | 4 |
| 2 | 24 | integer semigroup | 4 | `[1, 1, -2, -1, -1, 2]` | 4 |
| 3 | 26 | integer semigroup | 3 | `[1, -2, -1, -1, 3]` | 3 |
| 4 | 22 | rational cone | 4 | `[1, -2, -1, -1, 2, 1]` | 4 |
| 5 | 24 | integer semigroup | 4 | `[1, -1, -1, 1, -3, 3]` | 4 |
| 6 | 12 | rational cone | 3 | `[1, -1, -1, -2, 3]` | 3 |
| 7 | 10 | rational cone | 3 | `[1, -2, -1, 3, -1]` | 3 |
| 8 | 10 | rational cone | 3 | `[1, -2, 3, -1, -1]` | 3 |

The report also preserves the first triangulation witness, branch diagnostic,
and failed CMS-general divisor-intersection checks for each target. The first
target, for example, has relation
`[(0,-2),(55,-1),(56,-2),(202,1),(208,2),(211,2)]`, two witnesses, and both
candidate shrinking divisors fail `has_rational_divisor_solution=false`.
All nine local charge rows have zero row sum, so they pass the Calabi-Yau charge
sum check. However, the five six-point affine-rank-four targets are fourfold-
shaped under cygv's compact hypersurface dimension count
`q_rows - q_cols - 1`, not threefold-shaped. The four five-point
affine-rank-three targets are threefold-shaped by this gate, but this is only a
shape check: Cyrus still needs a source-derived compact/local semigroup, grading,
intersection data, and chamber interpretation before any HKTY output can be
promoted.

This narrows the next source-read task: the existing CKYZ rank-two surface path
is the wrong abstraction for these corrected-chamber misses. The missing
implementation needs a source-derived HKTY/local-chamber construction for
one-parameter, affine-rank-three/four origin circuits, or an explicit
flop/chamber-continuation rule that explains why these source domains should be
continued from another chamber. Until that source is identified, these rows
remain research gaps rather than GV values Cyrus can promote.
