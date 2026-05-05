# Potent-Ray Source Read

This note records the current no-code read-through for the McAllister potent-ray
gap. It is deliberately narrow: the goal is to prevent another round of fitting
`potent_rays_gv.dat` instead of reconstructing the geometry that produced it.

## Sources Read

- `reference/cytools/src/cytools/calabiyau.py`
  - `_compute_gvs_gws`
  - `mori_cone_cap`
- `reference/cytools/src/cytools/cone.py`
  - `find_lattice_points`
- `cygv-0.1.2/src`
  - `hkty.rs`
  - `semigroup.rs`
  - `fundamental_period.rs`
  - `instanton.rs`
  - `series_inversion.rs`
- `research/papers/small_cc_2107.09064_source/AdS4_v3.tex`
  - Sections on Gopakumar-Vafa invariants and radius of convergence
- `research/papers/mcallister_moduli_space_reconstruction_2212.10573.tex`
  - Flop/nop/M_infinity classification and chamber-continuation formulas
- `resources/small_cc_2107.09064_source/anc/paper_data/readme.txt`
  - Semantics of `potent_rays*.dat` and `small_curves*.dat`

## What CYTools/cygv Actually Computes

CYTools' public `cy.compute_gvs()` wrapper is compact-hypersurface HKTY:

1. build `mori = cy.mori_cone_cap(in_basis=True)`;
2. augment Mori-cap rays with `mori.find_lattice_points(min_points=100*h11)`;
3. pass generators, grading vector, curve-basis matrix `q`, and triple
   intersections into `cygv.compute_gv`;
4. let `cygv` close the supplied generators into a semigroup before computing
   the fundamental period, instanton data, and degree-ordered series inversion.

The supplied generator list is not merely a list of curve classes to evaluate.
It defines the semigroup context, and the final GV value for a class depends on
lower-degree subtraction history inside that context.

`cygv` is also not a local toric threefold engine. In `fundamental_period.rs`,
the effective Calabi-Yau dimension is inferred from the GLSM charge matrix and
nef partition, and dimensions below three are rejected. So the local `P^2`
sequence visible in the first potent ray cannot honestly be obtained by
pretending the charge vector `(-3,1,1,1)` is just a compact CYTools/cygv
hypersurface input.

## What The Paper Says For Potent Rays

The McAllister paper does not claim that the 214-dimensional Kähler-side
potent-ray data comes from a plain global `cy.compute_gvs()` call. It says:

- systematic high-degree GV computation at large `h11` is infeasible;
- potent rays control the radius of convergence through infinite GV families;
- high-degree GV data is computed inside low-dimensional Mori faces;
- one should find a phase where a face of `M_infinity(X)` is also a face of the
  Mori cone, then compute GV values in that low-dimensional context.

The ancillary files match this interpretation. `potent_rays.dat` stores saved
ambient curve classes, `potent_rays_gv.dat` stores the first ten GV values on
those saved rays, and `potent_rays_rank.dat` / `potent_rays_vols.dat` store
validation summaries. They are checkpoints, not a reusable production input.

## What The 4-214 Checkpoint Reveals

The first saved potent ray has support
`[(43, 1), (155, 1), (168, -3), (169, 1)]`. The triangulation-point coordinates
satisfy

```text
p43 + p155 - 3*p168 + p169 = 0
```

with coefficient sum zero, equivalently

```text
p168 = (p43 + p155 + p169)/3
```

Its saved GV sequence begins

```text
3, -6, 27, -192, 1695, ...
```

which is the standard local `P^2` sequence. The new Cyrus validation test shows
that all 411 saved potent rays are affine toric circuits. Their support ranks
are not all the same: 395 have affine rank 2, 16 have affine rank 4, and 56 of
the rank-2 circuits have this local `P^2` triangle pattern. For the first ray,
Cyrus now reconstructs deterministic rank-two local coordinates
`p43=(1,0)`, `p155=(0,1)`, `p168=(0,0)`, `p169=(-1,-1)` from the upstream
lattice points.

This is a structural clue, not permission to hardcode GV sequences. The same
four indices are not a current-phase adjacent two-face circuit in the
McAllister `heights.dat` triangulation, so the missing operation is a local
face/phase reconstruction, not a lookup against the current FRST.

## What This Rules Out

- Do not treat `potent_rays_gv.dat` as a runtime source for Cyrus.
- Do not add a hardcoded library of local `P^2` or weighted-projective GV
  sequences to make the checkpoint pass.
- Do not expect a one-generator `cygv` call to reproduce `potent_rays_gv.dat`;
  it asks a different semigroup question.
- Do not keep expanding final-volume fitting diagnostics until the face/ray GV
  source problem is solved.

## Next Implementation Target

The next useful Cyrus work is a focused reconstruction test for one saved
potent ray:

1. start from `points.dat` and `heights.dat`;
2. detect the affine circuit support of the first potent ray without using its
   saved GV values;
3. reconstruct the low-dimensional toric face/surface model or phase context
   associated with that rank-2 circuit, while keeping the rank-4 cases
   explicit as a separate class;
4. compute `N_q, ..., N_10q` from that reconstructed context;
5. compare to row 0 of `potent_rays_gv.dat` only as the assertion.

Only after that single-ray test is understood should Cyrus try to regenerate
the full 411-ray sample or use potent-ray data in a GA fitness pipeline.
