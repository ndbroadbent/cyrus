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

The current report artifact is
`/tmp/cyrus_six_point_rank3_source_model_target7_v1.json`. It records:

```text
source_model = weighted_p2_rank_three_source_model_visible_phase_is_local_cy5_total_space_not_promotable_gv_source
total_dim    = 5
required_cygv_codim = 2
numerical_gv = weighted_p2_rank_three_visible_phase_is_not_numerical_cy3_requires_source_codim2_or_insertion_history
```

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
