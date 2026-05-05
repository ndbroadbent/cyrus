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

## Current Cyrus Situation

Cyrus already has wrappers around the same Rust `cygv` crate:

- `compute_gv_invariants`
- `compute_gv_invariants_with_provided_generators`
- `compute_gv_invariants_with_explicit_semigroup`

But the McAllister corrected-chamber diagnostic path currently relies first on
toric two-face/origin-circuit formula values for selected small curves, then uses
local one-dimensional ray or LP-witness face HKTY diagnostics only for misses.
That is not equivalent to CYTools' full `cy.compute_gvs(min_points=N)` call.

The latest Stage 5 residual is therefore best interpreted as an exact global GV
coverage/convention problem, not as a dilog aggregation problem and not as a
problem in the already-checked local toric formulas.

## Next Productive Step

Stop expanding local per-curve diagnostics. The next useful comparison is a
faithful CYTools-style global GV input/output audit:

1. Dump Cyrus' corrected-chamber global GV inputs:
   Mori-cap basis rays, CYTools-style lattice augmentation rows, grading vector,
   `q` matrix, and intersection-number dictionary.
2. Generate the same dump from CYTools for the same corrected-chamber
   triangulation.
3. Compare those inputs before running `cygv`.
4. Only after the inputs match, compare the `cygv` outputs by basis charge and
   mapped ambient charge.

If the full corrected chamber is too large, the reduced test should still use
the same global semigroup construction rule, not an isolated face-local
semigroup, so that the HKTY subtraction order remains comparable.
