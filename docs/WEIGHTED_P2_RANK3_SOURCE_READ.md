# Weighted P2 Rank-Three Source Read

This note records the current source boundary for the remaining six-point
target `7`/`8` corrected-chamber blocker. It is deliberately narrow: the goal is
to prevent the named local phase from being promoted into a GV value without a
source-derived threefold model and qN history.

## Current Cyrus Object

The shared target `7`/`8` chamber-generator row has point relation

```text
[(0,-1),(55,-2),(208,1),(211,1),(212,2),(214,-1)]
```

The local one-parameter charge classifier identifies it as

```text
local_toric_weighted_p2_rank_three_split_bundle_charge_family:
base=1,1,2;bundle=1,1,2;base_hyperplane_square=1/2
```

The corrected global height chooses the base phase:

```text
weighted_p2_rank_three_split_bundle_selected_base_phase
base points:   [208,211,212], weights [1,1,2]
bundle points: [0,214,55],    degrees [1,1,2]
```

The visible geometry is therefore the total-space phase of

```text
O(-1) + O(-1) + O(-2) -> P(1,1,2)
```

This is useful chamber data, but it is not a numerical CY3 GV source by itself.
The base has complex dimension `2` and the split bundle has rank `3`, so the
visible total space has complex dimension `5`. A `cygv` threefold input would
need an additional codimension-`2` source model or an equivalent insertion/qN
history.

The current report artifacts are:

```text
/tmp/cyrus_six_point_rank3_source_model_target7_v1.json
/tmp/cyrus_six_point_rank3_source_model_target8_v1.json
/tmp/cyrus_six_point_rank3_source_model_target7_v2.json
/tmp/cyrus_six_point_rank3_source_model_target8_v2.json
/tmp/cyrus_six_point_rank3_partition_samples_target7_v1.json
/tmp/cyrus_six_point_rank3_partition_samples_target8_v1.json
/tmp/cyrus_six_point_rank3_origin_included_partitions_target7_v1.json
/tmp/cyrus_six_point_rank3_origin_included_partitions_target8_v1.json
/tmp/cyrus_six_point_rank3_reflexivity_precondition_target7_v1.json
/tmp/cyrus_six_point_rank3_reflexivity_precondition_target8_v1.json
/tmp/cyrus_six_point_origin_split_roles_target7_v1.json
/tmp/cyrus_six_point_origin_split_roles_target8_v1.json
```

The v2 artifacts record the same source-model obstruction plus the source
invariants of the visible weighted-`P2` bundle:

```text
source_model = weighted_p2_rank_three_source_model_visible_phase_is_local_cy5_total_space_not_promotable_gv_source
total_dim    = 5
required_cygv_codim = 2
numerical_gv = weighted_p2_rank_three_visible_phase_is_not_numerical_cy3_requires_source_codim2_or_insertion_history
base_anticanonical_degree = 4
bundle_degree_sum = 4
total_first_chern_degree = 0
base_hyperplane_square = 1/2
base_tensor_status = weighted_p2_rank_three_base_hyperplane_square_fractional_requires_stack_or_source_tensor_normalization
```

The partition-sample artifacts add a bounded raw sample of the 15 codimension-2
partition candidates. The first sampled partition is

```text
first part  = [55], degree [2]
second part = [208,211,212,214], degree [-3]
target relation part sums = -2 and 3
```

and it fails the CYTools-style nef certificate because the partition union hull
does not equal the ambient support-polytope hull, with ambient vertex `0`
missing. The target `8` report has the same sampled obstruction. This makes the
current source gap concrete: the obvious codimension-2 partition enumeration is
not merely missing a zero-degree candidate; its raw candidates also fail the
support-hull certificate.

The origin-included partition artifacts test the nearest tempting repair:
allow point `0` into a partition part instead of following the CYTools
nef-partition convention. For both target `7` and target `8`, this gives `31`
origin-included bipartitions and `6` zero-degree balanced candidates. All `31`
still fail the CYTools-style certificate, now with
`invalid_input_nef_partition_parts_must_exclude_the_origin_index`. Thus the
missing codimension-2 source is not just the excluded origin; accepting the
zero-degree origin-containing splits would violate the input contract we need
for a faithful `cygv`/CYTools handoff.

The origin-split role artifacts classify those six zero-degree candidates
against the selected weighted-`P2` base/bundle phase. For both target rows the
summary is
`weighted_p2_rank_three_origin_included_zero_degree_splits_mix_base_bundle_and_violate_origin_contract`.
The role-signature counts are:

```text
first=base:2,bundle:1,origin:1;second=base:1,bundle:1  -> 3
first=base:1,origin:1;second=base:2,bundle:2           -> 2
first=base:1,bundle:1,origin:1;second=base:2,bundle:1  -> 1
```

So the zero-degree balanced splits are not clean base-vs-bundle phase
separations. They mix base and bundle roles and still include the origin in a
nef part, which keeps them diagnostic-only.

The current-chamber unresolved-generator report now promotes the important
blocker classifications to top-level unique-generator aggregate counts, not
just sample fields. In addition to the origin-split status, the report counts
the selected weighted-`P2` phase, visible source-model status, fractional
base-tensor status, numerical-GV status, and support-reflexivity precondition.
For this row those aggregates should expose the same obstruction in one place:
selected base phase, visible local CY5 rather than CY3, required codimension
`2`, fractional `H^2=1/2` tensor normalization, and non-reflexive support with
the origin as a hull vertex.

A direct CYTools/PALP probe of the same six support coordinates is also
negative. With experimental CYTools features enabled, `Polytope(pts)` reports
dimension `4`, point count `6`, and `is_reflexive=False`; calling
`nef_partitions(codim=2)` raises `ValueError: The polytope must be reflexive`.
So the visible six-point support is not itself a PALP nef-partition source. A
valid source must enlarge/change the source polytope, provide a stacky/twisted
model, or supply the missing chamber/qN history by another sourced route.
Cyrus now exports the native precondition behind this result as
`support_polytope_not_reflexive_origin_is_hull_vertex`: for both target reports,
the local support hull vertices are `[0,55,208,211,212,214]`, so the lattice
origin is on the hull boundary rather than in the interior.

## Sources Read

### cygv 0.1.2

Relevant file:

```text
cygv-0.1.2/src/fundamental_period.rs
```

`compute_omega` infers

```text
ambient_dim = q.nrows - q.ncols
cy_codim    = 1 if no nef partition is supplied, else nefpart.len()
cy_dim      = ambient_dim - cy_codim
```

and rejects `cy_dim < 3`. Cyrus' compact wrapper also converts intersection
numbers to `i32`, so the current compact `cygv` path cannot accept a fractional
one-parameter tensor such as the naive `P(1,1,2)` hyperplane square `1/2`.

For the six-point rank-one visible charge row, the no-nefpart compact wrapper
sees `ambient_dim = 5` and `cy_dim = 4`, not a CY3. A codimension-2 nef
partition would be needed before this can look like a threefold to `cygv`.
Cyrus' CYTools-style codimension-2 scan currently finds no certified
zero-degree nef partition for this support.

### CYTools CICY validation

Relevant files:

```text
reference/cytools/src/cytools/calabiyau.py
reference/cytools/src/cytools/polytope.py
```

`CalabiYau.__init__` validates a supplied nef partition by requiring the union
of the non-origin parts to reconstruct the ambient polytope, then adding the
origin to each part-polytope and checking that their Minkowski sum is reflexive.
`Polytope.nef_partitions` delegates to PALP, but first requires the polytope
itself to be reflexive. These checks match the Cyrus diagnostics above: the
raw non-origin split misses the origin as an ambient hull vertex, while the
origin-containing split violates the non-origin part convention.

### Coates-Corti-Iritani-Tseng, arXiv:0804.2592

Relevant local source:

```text
/tmp/cyrus_source_read/wallcrossings2.tex
```

The paper is relevant but not a direct source for this row.

It treats canonical weighted-projective examples, including
`K_P(1,1,2)`. In the paper's Example II, the toric action has two charge rows
on five coordinates:

```text
[0,0,1,1,-2]
[1,1,0,-2,0]
```

The chamber giving `K_P(1,1,2)` is therefore a rank-two canonical-bundle model,
not the Cyrus rank-one six-coordinate split-bundle phase
`[-1,-1,-2,1,1,2]`.

The same paper also states the general hypergeometric modification for twisted
Gromov-Witten invariants of total spaces of line-bundle sums over compact
orbifolds. That theorem explains how a vector-bundle source should enter an
I-function. It does not by itself give the missing numerical CY3 GV row for the
Cyrus six-point object, because the visible split-bundle total space is a
five-dimensional local Calabi-Yau phase unless an additional codimension-2
source model or insertion history is specified.

A follow-up local-paper grep for the exact model found no
`O(-1)+O(-1)+O(-2) -> P(1,1,2)` source. The closest adjacent example in the
same paper is the toric flop between `O(-1)+O(-2) -> P^2` and
`O(-1)^3 -> P(1,2)`. That is useful evidence that the hypergeometric
line-bundle-sum machinery is the right mathematical neighborhood, but it still
has different base weights, bundle degrees, and local dimension from the Cyrus
six-coordinate row.

## What This Rules Out

- Do not reuse the canonical `K_P(1,1,2)` formula as if it were the same
  geometry. It has a different charge matrix and chamber.
- Do not call the one-parameter compact `cygv` wrapper with a guessed unit
  tensor. That would compute a different compact/hypersurface problem.
- Do not promote the selected GIT chamber into a scalar GV value. The selected
  chamber identifies a phase, not the threefold source and qN subtraction
  history.
- Do not use the McAllister saved GV row to choose the missing tensor or branch.
  That would be downstream replay.
- Do not replace the six-point source by the nearest degree-2 toric overlap
  row. The target-filtered overlap reports show a unique best toric neighbor on
  the positive height side for targets 7 and 8, but it still misses three of
  the six support points and adds an external point.
- Do not reuse the adjacent `O(-1)^3 -> P(1,2)` flop example as if it were the
  target row. It is a rank-three bundle example, but it is not the same
  weighted base or bundle.

## Next Source Task

The missing object is one of:

1. a source-derived codimension-2 local model for the selected rank-three
   weighted-`P2` phase, including a certified nef/complete-intersection
   description, intersection tensor, and chamber certificate;
2. a compact chamber-semigroup construction that hands the correct finite
   qN-history domain to `cygv`; or
3. a sourced twisted/vector-bundle I-function calculation with the insertions
   needed to extract the same numerical CY3 contribution.

Until one of those exists, this row must remain blocked at
`source_ray_weighted_p2_rank_three_split_bundle_source_import_blocked_missing_source_model_tensor_qn_history`.

Implementation note: Cyrus now has `cygv` qN-trace wrappers that accept an
explicit nef partition for both provided-generator and explicit-semigroup
domains. The blocker here is therefore not the ability to call upstream `cygv`
with codimension-2 qN tracing; it is the absence of a source-certified
codimension-2 model/nef partition/tensor/chamber history for this specific
rank-three weighted-`P2` row.
