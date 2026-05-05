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
- arXiv:2303.00757, *Computational Mirror Symmetry*
  - consistent truncation / causality constraint
  - degree and past-light-cone methods
  - low-dimensional face/ray computations
  - direct toric-curve GV formulas
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

The local CKYZ path Cyrus is building mirrors only the period-coefficient layer
of that machinery, then performs the local-surface steps that are absent from
compact cygv. The double-log/prepotential-period routine computes the source
CKYZ `rho`-derivatives in B-model `z` coordinates, and the inverse mirror-map
helper composes those coefficients into flat `q` coordinates. The local
instanton-potential helper mirrors cygv's essential `beta - alpha alpha`
conversion and applies the CKYZ `instbase` multiple-cover inversion for
source-derived finite-limit weights.

This is now validated for local `P^2`, `F0`, `F1`, and polygon 5. The gated
McAllister rank-two CKYZ test checks the first three saved GV entries for all
395 rank-two CKYZ potent-ray rows from source relation data, and a focused
first-principles regression checks all ten saved entries for the canonical F1
source direction `[2,1]`. The saved `potent_rays_gv.dat` values remain
validation targets only. The full ten-entry rows across all families still need
a sharper coefficient-targeted extractor; the current cover-closed past-downset
domain is source-derived but still too broad for the largest target directions.

The computational-mirror paper explains the implementation boundary that this
source code implies. Finite HKTY computations need a causally closed truncation:
either a degree cutoff from an interior grading vector, or a past-light-cone
set. At high `h11`, useful computations are often deliberately restricted to
low-dimensional Mori faces. A one-dimensional ray computation is valid only if
the ray is a genuine effective Mori generator; integrality of the HKTY output is
used as evidence for that condition, not as a universal guarantee.

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

This also corrects an older reproduction-script simplification: generic
`cy.compute_gvs(min_points=...)` is appropriate for the low-dimensional mirror
racetrack data and for smoke tests, but the high-dimensional Kähler-side potent
rays require reconstructing the low-dimensional face/ray context used for the
HKTY truncation.

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

## Cross-Example Pattern Audit

A source-reading pass over all five McAllister ancillary examples reinforces the
same conclusion. The saved potent rays are highly structured, but the structure
is a small family of local toric charge patterns, not a license to replay the
saved GV rows.

| Example | Potent rays | Coefficient patterns | GV-sequence starts |
|---------|-------------|----------------------|--------------------|
| `4-214-647` | 411 | 32 | 24 |
| `5-113-4627-main` | 1729 | 136 | 51 |
| `5-113-4627-alternative` | 1729 | 136 | 51 |
| `5-81-3213` | 727 | 80 | 40 |
| `7-51-13590` | 758 | 133 | 76 |

Across the examples, the most common sorted nonzero coefficient patterns are
sparse rank-two local-toric-looking relations:

- `(-3, 1, 1, 1)`, whose GV sequence starts
  `3, -6, 27, -192, 1695`;
- `(-5, 1, 1, 1, 2)`, whose GV sequence starts
  `5, -110, 6885, -672000, 83508575`;
- `(-6, 1, 1, 2, 2)`, whose GV sequence starts
  `-6, -288, -40338, -8757888, -2423174610`;
- `(-7, 1, 1, 2, 3)`, whose GV sequence starts
  `7, -644, 171801, -71340864, 37754196305`;
- `(-8, 1, 1, 3, 3)` and `(-8, 1, 2, 2, 3)`;
- `(-10, 1, 1, 4, 4)` and `(-10, 2, 2, 3, 3)`.

This changes the next implementation target. We should not write a generic
"coefficient pattern -> GV row" table. The next real port is a local toric
surface/face reconstruction that derives the charge context and then computes
the one-parameter or low-dimensional mirror series from that context. The local
`P^2` implementation is one verified instance of that shape, not the general
solution.

## What This Rules Out

- Do not treat `potent_rays_gv.dat` as a runtime source for Cyrus.
- Do not add a hardcoded library of local `P^2` or weighted-projective GV
  sequences to make the checkpoint pass.
- Do not expect a one-generator `cygv` call to reproduce `potent_rays_gv.dat`;
  it asks a different semigroup question.
- Do not overinterpret the existing ignored one-generator test sequence. Its
  current `q`-matrix construction includes the origin/canonical column, whereas
  the CYTools GV path uses `curve_basis(include_origin=False)`. The mismatch
  strengthens the conclusion that the test is diagnostic-only; it does not
  provide a trusted local model.
- Do not keep expanding final-volume fitting diagnostics until the face/ray GV
  source problem is solved.

## Implemented Rank-Two CKYZ Reconstruction

Cyrus now has non-replay potent-ray GV reconstruction for the rank-two CKYZ
surface families in the 4-214-647 checkpoint:

1. detect the affine circuit support of the first potent ray from
   `points.dat`/`potent_rays.dat`, without reading its GV row;
2. reconstruct rank-two local coordinates and the local affine charge lattice;
3. identify local `P^2`, `F0`, `F1`, and CKYZ polygon-5 source relation systems
   by point permutation and unimodular row transformation;
4. compute CKYZ local period corrections, mirror-map substitution,
   `beta - alpha alpha`, and multiple-cover extraction in source coordinates;
5. use `potent_rays_gv.dat` only as the assertion.

The gated test is
`mcallister_rank_two_ckyz_potent_ray_gvs_are_reconstructed` in
`crates/cyrus-core/tests/potent_ray_affine_circuits.rs`.

The remaining potent-ray target is to extend this from first-three rank-two
checks to the complete saved rows where feasible, keep the rank-four cases
explicit, and then generate the full 411-ray sample without treating
`potent_rays*.dat` as production input.

As a narrower checkpoint toward that target, Cyrus now reconstructs integral
two-dimensional local coordinates for every rank-two saved potent-ray support
from the upstream lattice points. The reconstruction uses the affine difference
lattice of the support, not a coefficient-pattern table. Cyrus also derives a
normalized local-support signature by dropping point labels, normalizing the
overall relation sign, and canonicalizing over integral coordinate systems
generated by support-point differences. The gated 4-214-647 tests verify:

- 395 rank-two supports expose local 2D coordinates;
- 16 rank-four supports remain explicit and do not expose fake rank-two
  coordinates;
- each rank-two affine relation still sums to zero in the reconstructed local
  coordinates;
- each saved sparse relation lies in the rational row span of the reconstructed
  local affine charge basis.
- all 56 recognized local `P^2` supports collapse to one normalized support
  signature, while the signature set still distinguishes non-`P^2` rank-two
  families.

The audited rank-two inventory for `4-214-647` now has 16 normalized support
signatures across the 395 rank-two supports. Grouped only by coefficient
pattern, the counts are:

| Coefficients | Count |
|--------------|------:|
| `(-3, 1, 1, 1)` | 56 |
| `(-4, 1, 1, 1, 1)` | 17 |
| `(-5, 1, 1, 1, 2)` | 34 |
| `(-6, 1, 1, 2, 2)` | 31 |
| `(-7, 1, 1, 2, 3)` | 32 |
| `(-7, 1, 1, 1, 2, 2)` | 1 |
| `(-8, 1, 1, 3, 3)` | 34 |
| `(-8, 1, 2, 2, 3)` | 32 |
| `(-9, 1, 1, 3, 4)` | 25 |
| `(-9, 1, 1, 2, 2, 3)` | 1 |
| `(-10, 1, 1, 4, 4)` | 32 |
| `(-10, 2, 2, 3, 3)` | 34 |
| `(-11, 1, 1, 4, 5)` | 19 |
| `(-11, 1, 3, 3, 4)` | 9 |
| `(-12, 1, 1, 5, 5)` | 31 |
| `(-14, 1, 4, 4, 5)` | 7 |

This still stops short of assigning GV sequences to the non-`P^2` patterns. The
next missing source-derived object is the local mirror/HKTY input attached to
each reconstructed support.

As a step toward that input, Cyrus now constructs a canonical local charge
model from each normalized rank-two support signature. The model recomputes the
integer kernel of `[1; x; y]` in the canonical local coordinates and verifies
that the saved potent-ray relation belongs to that local charge lattice. This
removes another point-label and checkpoint-data dependency, but it deliberately
does not decide which local mirror/HKTY series, if any, belongs to a non-`P^2`
family.

The local charge model now also records the integer coordinates of the target
potent-ray relation in the reconstructed local charge basis. This is important
for non-`P^2` supports, where the local charge lattice has rank greater than
one and the saved ray is a specific one-parameter direction inside that
lattice.

Cyrus also has a reusable CKYZ source matcher for those rank-two models. It
does exact column-permutation and unimodular-row-transform matching against the
CKYZ relation systems for local `P^2`, `F0`, `F1`, and polygon 5, then reports
the target ray in CKYZ source coordinates together with the source `c1` pairings
and local intersection expression. This is local mirror input construction, not
GV assignment.

## Deeper Internal Read Boundary

A follow-up source pass confirms the boundary for the next code change:

- CYTools' compact `compute_gvs()` wrapper is not a local toric-surface GV
  engine. It normalizes compact CY data into semigroup generators, a grading
  vector, a curve-basis `q` matrix, and triple intersections, then calls cygv.
- cygv's HKTY result depends on the whole semigroup truncation. The saved
  potent-ray row is not enough to reconstruct the calculation unless Cyrus also
  derives the local face/semigroup context that supplies the lower-degree
  subtraction history.
- CYTools' `mori_cone_cap()` and `Cone.find_lattice_points()` are part of that
  context. They determine which curve relations and lattice points seed the
  semigroup before the fundamental-period and mirror-map computation starts.
- The old Python reproduction wording that treats high-dimensional GV
  computation as one generic `cy.compute_gvs(max_deg=N)` call is too coarse.
  It is acceptable for low-dimensional mirror/racetrack validation, but not for
  the Kähler-side potent-ray production path.

So the next implementation should first produce a source-derived inventory of
rank-two local support families from the normalized signatures, then derive the
local mirror/HKTY input for one non-`P^2` family. It should not assign a GV
sequence merely because a coefficient pattern matches an ancillary row.

## Code-Level Targeted-Extraction Implications

Reading `cygv` more closely gives a sharper implementation constraint for the
rank-two CKYZ performance gap. The current Cyrus targeted CKYZ path is
mathematically source-derived, but it still truncates formal series by the
componentwise box spanned by the requested degrees. That is the source of the
all-ten-row blowup: for a target like `n * (a,b,c)`, the box contains many
monomials that cannot affect the requested coefficients.

`cygv` avoids this shape of waste by doing all polynomial arithmetic over a
finite semigroup object. `Semigroup::with_max_degree` and
`Semigroup::with_min_elements` first close the supplied generators under
addition, then `Polynomial::mul` emits a product only when the summed exponent
is present in `poly_props.monomial_map`. In other words, the truncation is an
explicit finite monomial set with an addition lookup, not an axis-aligned box.

The corresponding CKYZ fix should therefore be a coefficient-targeted monomial
domain:

1. start from the target ray multiples and their cover divisors;
2. close under the lower-degree monomials needed by log-period corrections,
   inverse mirror-map substitution, exponentials, and local
   `beta - alpha alpha` multiplication;
3. perform series addition/multiplication/composition against that domain's
   addition map;
4. extract only the requested GV degrees after the same multiple-cover
   subtraction used by the full CKYZ path.

This is a source-derived analogue of cygv's semigroup arithmetic. It is not a
coefficient-pattern table and it is not a one-dimensional shortcut. The saved
`potent_rays_gv.dat` row should still enter only at the final assertion stage.

The cygv series-inversion code also reinforces that "just compute the ray" is
not a well-defined independent operation. `invert_series` processes classes in
grading order and subtracts the `Li2(q^N)` history of previously found classes.
For the CKYZ local-surface path this history is mirrored by the
multiple-cover/primitive-divisor subtraction; for rank-four local affine
supports, Cyrus still needs the actual local semigroup, `q` matrix, and
intersection/Yukawa data before any ray comparison is meaningful.

A post-cache benchmark confirms that broader componentwise-box expansion is
still not an acceptable substitute for that targeted domain. Changing the
rank-two CKYZ potent-ray test locally to check all ten saved entries and running
it under a 600 second timeout still timed out. A shorter instrumented run showed
that this is not only a polygon-5 problem: local `P2` and `F0` direction `[1,1]`
finished, `F1` direction `[2,1]` finished, but `F0` direction `[1,2]` already
took roughly 90 seconds and the run timed out shortly after starting `F1`
direction `[3,1]`. The committed all-row gate must therefore remain at
first-three checks until the finite monomial-domain extractor exists, with
`[2,1]` kept as a focused all-ten non-`P^2` regression. A later indexed-domain
multiplication refactor removed much of the per-product lookup overhead but
still did not make all-ten validation finish under a 600 second timeout, so the
remaining blocker is the monomial set being too broad, not only the arithmetic
representation. A 240 second instrumented run after indexing progressed through
`P2` `[1]`, `F0` `[1,1]`, `F1` `[2,1]`, `F0` `[1,2]`, and `F1` `[3,1]`, then
timed out after starting `F0` `[1,3]`. That profile is better than the
pre-indexed run but still not acceptable as an all-ten validation path.

Flattening the indexed addition table from nested `Option` rows into a sentinel
array removes enough debug-build lookup overhead for the first hard polygon-5
direction `[4,3,2]` to pass the `N=4` narrowed McAllister check. Iterating only
valid sparse product pairs brings that narrowed run to about 149 seconds. A
stack sample still sits almost entirely in `compute_ckyz_inverse_mirror_map_domain`,
especially composition and rational polynomial multiplication. This is progress
on the source-derived path, but it does not change the conclusion: full
ten-entry and all-row validation still need a coefficient-history domain rather
than a broader downset.

## Current Implementation Boundary

The current Cyrus potent-ray code has reached a useful stopping point for
source reading:

- it can recover the affine support and local charge lattice of every saved
  rank-two potent-ray support in the 4-214-647 checkpoint;
- it can express the saved potent-ray relation as an integer direction in that
  local charge lattice;
- it can compute the first local `P^2` GV row from the one-parameter local
  mirror map, without reading the saved GV row as input.

Executable evidence is now split by layer:

- `mcallister_potent_rays_are_affine_toric_circuits` verifies the affine support
  inventory;
- `mcallister_potent_rays_have_local_affine_charge_contexts` verifies that each
  saved relation lies in the reconstructed local affine charge lattice;
- `mcallister_rank_two_local_supports_are_reflexive_polygons` verifies that the
  16 normalized rank-two families are two-dimensional reflexive polygons, with
  the compact point as the unique interior lattice point;
- `mcallister_rank_two_local_charge_models_are_inventoried` pins the canonical
  charge bases and target charge-lattice directions for all 16 rank-two support
  families;
- `mcallister_rank_two_models_identify_ckyz_sources` verifies the reusable CKYZ
  source matcher for all 16 rank-two support families;
- `mcallister_rank_four_local_affine_supports_are_inventoried` keeps the 16
  affine-rank-four rows explicit instead of projecting them into fake rank-two
  local coordinates;
- `mcallister_rank_two_ckyz_potent_ray_gvs_are_reconstructed` reaches GV values
  for all 395 rank-two CKYZ rows, currently checking the first three saved
  entries from source-derived local mirror data;
- `mcallister_f1_2_1_ckyz_potent_ray_all_ten_gvs_are_reconstructed` checks all
  ten saved entries for the canonical F1 source direction `[2,1]` without using
  the saved GV row as input.

This is not yet a general potent-ray GV engine. The rank-two CKYZ local-surface
path is source-derived, and its current targeted extractor uses the
cover-closed union of componentwise past downsets for the requested degrees.
That is sufficient for first-three validation, small local tables, and the full
F1 `[2,1]` saved row, but it is not the right performance model for all ten
saved entries across larger target directions. The next source-derived object
is a coefficient/path-history local-series extractor, not a one-dimensional
coefficient-pattern rule.

The rank-four checkpoint is deliberately separate from the CKYZ local-surface
path. For 4-214-647, those 16 rows collapse to one seven-point affine-rank-four
configuration with local charge basis
`[(1,3,-1,-1,-1,-1,0), (0,1,-1,-1,0,0,1)]`, plus one six-point face with basis
`(1,3,-1,-1,-1,-1)`. This proves the rows are structured local affine charge
contexts, but not two-dimensional reflexive polygons. Their GV source therefore
has to come from the low-dimensional face/semigroup construction, not from the
rank-two CKYZ surface matching.

The intended next code change should therefore be narrow:

1. choose one non-`P^2` normalized support family;
2. construct the local semigroup generators, grading vector, local `q` matrix,
   and local intersection/Yukawa data from the reconstructed support and charge
   basis;
3. run the lower-level cygv/HKTY path on that local input;
4. only then compare the resulting multiples of the target charge direction
   against the corresponding `potent_rays_gv.dat` row.

Anything that directly maps a coefficient pattern such as `(-5, 1, 1, 1, 2)`
to a saved GV sequence would reintroduce the same cheating this audit is meant
to avoid.

## Immediate Source-Reading Priority

Before writing more potent-ray GV code, the missing object has to be identified
from source, not inferred from the ancillary arrays. The current read-through
puts the next work in this order:

1. Keep CYTools/cygv compact-GV orchestration separate from local toric-surface
   mirror formulas. Compact `cy.compute_gvs()` supplies a semigroup, grading,
   curve-basis `q` matrix, and intersection tensor to cygv; it is not itself the
   CKYZ local-surface formula.
2. For rank-two CKYZ-matched families, keep improving the source-derived local
   mirror series from the CKYZ relation rows, finite-limit cover weights, and
   local intersection expression already returned by
   `identify_ckyz_local_surface`.
3. For rank-four affine supports, derive the actual low-dimensional face
   semigroup context before calling cygv. These rows should not be projected into
   fake rank-two polygons.
4. Use `potent_rays_gv.dat` only after the above inputs are constructed, as an
   assertion that the source-derived calculation produced the saved sequence.

This keeps the project from going in circles: the next code should be driven by
one fully read source path, not by another diagnostic wrapper around the same
saved data.

Item 2 now includes logarithmic periods, second-derivative/local-prepotential
periods, inverse mirror-map substitution, local `beta - alpha alpha`
conversion, and multiple-cover extraction. The remaining rank-two CKYZ work is
coefficient-targeted extraction for larger requested degrees, plus broader
validation beyond the all-row first-three gate and the focused F1 `[2,1]`
all-ten row.

## Source-Derived Implementation Checklist

The next implementation should be allowed to start only after these inputs are
available from geometry or source formulas, not from `potent_rays_gv.dat`:

1. canonical local support coordinates and integral charge basis for the chosen
   normalized rank-two family;
2. source relations, local intersection expression, and finite-limit cover
   weights for the matching CKYZ model;
3. target charge direction derived by reconstructing the saved potent-ray
   relation in the local charge lattice;
4. a finite monomial domain containing the target multiples, their cover
   divisors, and every lower-degree monomial that can change the
   `series_inversion.rs` Li2 subtraction history for those targets;
5. an assertion-only comparison against the corresponding saved GV row.

The fourth item is the current missing object. A componentwise box or ordinary
past downset is easy to construct, but cygv's final GV values are produced after
degree-ordered subtraction of previously discovered classes. The finite domain
therefore has to be closed for the coefficient history that can reach the
target coefficients, not just for monomial multiplication syntax.
