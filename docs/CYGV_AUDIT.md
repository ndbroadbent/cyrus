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
- `compute_gv_invariants_with_provided_generators`
- `compute_gv_invariants_with_explicit_semigroup`

Their boundaries are different:

- `compute_gv_invariants` tries to mirror the CYTools default augmentation, but
  currently changes the lattice-point request to `max_deg`-bounded enumeration
  when `max_deg` is supplied. That is useful for controlled diagnostics, but it
  is not the literal CYTools wrapper contract described above.
- `compute_gv_invariants_with_provided_generators` mirrors CYTools'
  `mcap_generators=...` entry point: Cyrus supplies rows directly to cygv and
  lets cygv perform its own semigroup closure.
- `compute_gv_invariants_with_explicit_semigroup` bypasses cygv closure with
  `Semigroup::from_data`. This is a diagnostic-only shape for hand-constructed
  local semigroups, not a replacement for CYTools' public `compute_gvs()`.

But the McAllister corrected-chamber diagnostic path currently relies first on
toric two-face/origin-circuit formula values for selected small curves, then uses
local one-dimensional ray or LP-witness face HKTY diagnostics only for misses.
That is not equivalent to CYTools' full `cy.compute_gvs(min_points=N)` call.

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

## Flop Continuation

The paper warns that continuing the Kähler-coordinate formula through a flop is
subtle. For curves with half-integral `B_2` field (`gamma dot q` odd), the
relevant real-axis dilogarithm continuation is:

```text
Li2(-exp(-2*pi*t))/(2*pi)^2
  = -Li2(-exp(-2*pi*(-t)))/(2*pi)^2 - 1/2*t^2 - 1/24
```

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

## Next Productive Step

Stop expanding top-contributor local checks. The next useful comparison is a
faithful selected-small-curve method audit:

1. Reconstruct, from code/paper, the exact rule that produced
   `small_curves.dat` and `small_curves_gv.dat`: chamber, volume cutoff,
   toric-curve enumeration, pair-decomposable pruning, and face/ray GV method.
2. Compare Cyrus against that rule in the approximate chamber first, where the
   ancillary files are direct checkpoints.
3. Then apply the same rule to the corrected chamber or, if the intended
   McAllister operation is analytic continuation through flops, port that
   transformation explicitly using the paper's polylog identities.
4. Use `cygv` full-semigroup input dumps only for low-dimensional face/ray
   subproblems where that is the method actually being used.

The open question is now sharper: is the residual caused by using the wrong
small-curve chamber/continuation convention, or by an incomplete port of the
toric/face GV method that generated the selected `small_curves_gv` values?

The next code change should be one of these two narrow audits, not another broad
GV fallback:

1. The selected-small-curve checkpoint now reproduces the published 344 selected
   classes and their `1/-2` GV values, using `small_curves*.dat` only as
   assertions.
2. Continue the corrected-chamber/flop-continuation audit beyond the 10 saved
   input-chamber negative curves: the simple odd-parity identity is verified,
   but the two even-parity saved curves and the broader corrected-chamber
   target-vector residual still need a first-principles treatment.

Concrete missing implementation pieces after the source read are:

- a faithful "sums of others" reduction, ideally Hilbert-basis/semigroup based
  rather than the current pair-only production pruning plus bounded diagnostic;
- the low-dimensional face/ray GV computation path used for large-`h11` potent
  and nilpotent curve checks, distinct from a full 214-dimensional global HKTY
  run;
- an explicit chamber/flop continuation rule for the Kähler-coordinate instanton
  terms, including what to do with even-parity curves that hit the
  `Li2(exp(...))` branch cut.
