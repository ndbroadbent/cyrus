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
