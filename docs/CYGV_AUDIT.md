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
