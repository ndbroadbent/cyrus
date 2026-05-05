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

## cygv Semantics

The `cygv` crate does not interpret `generators` as merely the list of curves to
evaluate. They define the affine semigroup used by the HKTY series computation.

The Rust call chain is:

1. `Semigroup::with_max_degree`, `with_min_elements`, or `from_data`
2. `fundamental_period::compute_omega`
3. `instanton::compute_instanton_data`
4. `series_inversion::invert_series`

Important details:

- `Semigroup::with_min_elements` starts from the supplied generators and closes
  the semigroup under addition until at least `min_points` elements exist.
- `compute_omega` builds `q0` from GLSM charges by summing rows for the
  hypersurface nef partition, then groups semigroup elements by the number of
  negative GLSM intersections.
- `invert_series` processes classes by grading degree and subtracts previously
  found lower-degree `Li2(qN)` contributions from the instanton corrections.

This means local face/ray HKTY checks can be correct for isolated classes while
still not reproducing the global CYTools/cygv output. The global result depends
on the full semigroup truncation and the degree-ordered subtraction history.

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

## Current Cyrus Situation

Cyrus already has wrappers around the same Rust `cygv` crate:

- `compute_gv_invariants`
- `compute_gv_invariants_with_provided_generators`
- `compute_gv_invariants_with_explicit_semigroup`

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

## Flop Continuation

The paper warns that continuing the Kähler-coordinate formula through a flop is
subtle. For curves with half-integral `B_2` field (`gamma dot q` odd), the
polylog identity used in the paper is:

```text
Li2(-exp(-2*pi*t))/(2*pi)^2
  = -Li2(-exp(-2*pi*(-t)))/(2*pi)^2 + 1/2*t^2 - 1/24
```

This polynomial correction is tied to the transformations of `chi(D_i)` and
`kappa_ijk` across the flop. Therefore simply evaluating the original chamber's
finite GV list at negative `q.t` is not a first-principles corrected-chamber
implementation unless the associated chamber/topology transformation is handled
explicitly.

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

1. Add a selected-small-curve checkpoint for the approximate chamber that must
   reproduce the published 344 selected classes and their `1/-2` GV values
   without reading `small_curves*.dat` except as assertions.
2. Add an explicit corrected-chamber/flop-continuation diagnostic for the 10
   selected classes whose corrected volumes are negative, verifying the
   required `gamma dot q` parity and the polynomial terms before promoting any
   continuation into the production KKLT correction.
