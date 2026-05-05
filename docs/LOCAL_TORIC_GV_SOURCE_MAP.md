# Local Toric GV Source Map

This note records the source-read boundary for the next potent-ray GV
implementation step. It is intentionally pre-code: the goal is to avoid
generalizing the local `P^2` routine by pattern matching on coefficients.

## Why This Exists

Cyrus now reconstructs rank-two potent-ray supports as reflexive polygons and
derives their local affine charge lattices. The remaining missing object is the
local mirror or topological-string computation that assigns GV series to those
local toric Calabi-Yau geometries.

The compact CYTools/cygv path does not provide that object:

- CYTools `compute_gvs()` builds compact hypersurface inputs and delegates to
  cygv.
- cygv computes HKTY over a compact-hypersurface semigroup context.
- cygv's `fundamental_period::compute_omega` infers a hypersurface codimension
  and rejects effective CY dimension below three.

So a noncompact local model such as the canonical bundle over a toric surface
must be sourced separately.

## Primary Algorithm Sources

### CKYZ Local Mirror Symmetry

The closest source for the next Cyrus implementation is Chiang-Klemm-Yau-
Zaslow, "Local Mirror Symmetry: Calculations and Interpretations",
arXiv:hep-th/9903053.

The relevant section is "Local mirror symmetry for the canonical line bundle of
a torically described surface". It gives a direct local-mirror construction
from a two-dimensional reflexive polygon:

1. Start with the base polygon `Delta_B` and points `nu_i`.
2. Lift each point to `bar_nu_i = (1, nu_i)`.
3. Choose a basis of integer relations `l^(a)` among the `bar_nu_i`; these span
   the Mori cone of the toric surface.
4. Define local complex-structure variables
   `z_a = product_i a_i^(l_i^(a))`.
5. Build Picard-Fuchs operators from the relation identities.
6. Use the hypergeometric solution
   `c(n,rho) = 1 / product_i Gamma(sum_a l_i^(a)(n_a+rho_a)+1)`.
7. Obtain the constant, logarithmic, and double-logarithmic periods by
   differentiating in the `rho_a` parameters.
8. Use the local mirror map and the local prepotential derivative to extract
   the instanton numbers.

This matches the shape Cyrus already has for local `P^2`: a relation lattice,
a mirror map, a local Yukawa/prepotential object, and a multiple-cover
inversion. The difference is that non-`P^2` polygons are usually
multi-parameter, and the saved potent ray is a direction inside the local
charge lattice rather than the whole lattice.

CKYZ also lists relation data for several two-dimensional reflexive polygons:

- polygon 1: local `P^2`, charge `(-3, 1, 1, 1)`;
- polygons 2-4: local Hirzebruch surfaces `F0`, `F1`, and `F2`;
- polygons 5, 6, and 11: higher-rank examples with extra Picard-Fuchs
  operators from sums of relation vectors;
- the remaining cases are described as having non-simplicial Mori cones, so
  the implementation must handle chamber/subcone choices rather than assuming a
  single canonical coordinate chart.

This is a better source for Cyrus' next step than a coefficient-pattern table:
the local charge model must first be identified with one of these polygon
relation systems up to point permutation and `GL(r,Z)` basis change.

### Topological Vertex

Aganagic-Klemm-Marino-Vafa, "The Topological Vertex",
arXiv:hep-th/0305132, gives a complete A-model algorithm for arbitrary local
toric Calabi-Yau threefolds:

1. Read the toric diagram as a planar graph whose edges carry primitive
   integral directions.
2. Decompose the graph into trivalent `C^3` vertices.
3. Put a Young-tableau representation on each internal edge.
4. Glue vertices by summing over representations, with edge Kahler weights and
   framing factors.
5. Extract genus-zero GV/BPS integers from the resulting closed-string free
   energy expansion.

This route is more general than CKYZ local mirror formulas, but it is a larger
implementation. It requires a partition/Young-diagram engine, Schur/Hopf-link
building blocks, framing conventions, graph decomposition, and BPS extraction.

For the immediate McAllister potent-ray gap, CKYZ is the smaller source-derived
path because Cyrus already has local polygon supports and charge lattices.
Topological vertex remains the fallback if CKYZ-style local mirror formulas
become too special or ambiguous for some of the 16 rank-two families.

## Consequences For Cyrus

The next Cyrus artifact should be local input construction, not GV assignment.
For one non-`P^2` rank-two family, preferably the common
`(-5, 1, 1, 1, 2)` support:

1. Reconstruct the canonical polygon support and its local charge basis from
   Cyrus' existing rank-two local charge model.
2. Identify the polygon relation system against CKYZ data up to column
   permutation and unimodular row transformation.
3. Record the saved potent-ray target direction in that source-derived local
   Mori coordinate system.
4. Only after this identification, implement the relevant local mirror series
   and compare the multiples of the target direction against
   `potent_rays_gv.dat`.

This keeps the saved GV file as an assertion. It also keeps the failure mode
honest: if a polygon cannot be matched to a source relation system or a chamber
choice is ambiguous, Cyrus should report that as a missing local-mirror input,
not silently fall back to the checkpoint row.

## Immediate Non-Code Checks

Before writing a local GV series for any non-`P^2` family, verify:

- the support polygon is matched by lattice isomorphism, not just by sorted
  coefficient pattern;
- the charge basis is matched by an explicit integer unimodular transform;
- the target relation coordinates are transformed through the same basis
  change;
- the chosen CKYZ chamber/subcone is recorded when the Mori cone is
  non-simplicial;
- the comparison against `potent_rays_gv.dat` is gated as validation-only.

The local `P^2` code is therefore a proof of approach, not a template for
changing one charge vector and expecting correct GV values.
