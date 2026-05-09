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
ckyz_source = weighted_p2_rank_three_ckyz_local_surface_source_not_applicable_visible_phase_not_local_surface_cy3
twisted_ifunction_source = weighted_p2_rank_three_twisted_vector_bundle_ifunction_candidate_requires_insertions_and_qn_history
twisted_ifunction_required_insertion_codim = 2
twisted_ifunction_candidate_insertion = base_hyperplane_power_2
twisted_ifunction_readiness = weighted_p2_rank_three_twisted_ifunction_blocked_missing_stack_normalized_codim2_insertion_qn_history
twisted_ifunction_missing = source_derived_codim2_insertion_or_equivalent_observable, stack_normalized_hyperplane_square_tensor, orbifold_sector_pairing_data, equivariant_residue_or_pairing_normalization, twisted_vector_bundle_ifunction_chamber_certificate, twisted_vector_bundle_ifunction_qn_history
twisted_ifunction_degree_profile = d=1/2 half-degree twisted sector, numerator zero order 1; d=1 untwisted sector, numerator zero order 3 > codim-2 insertion
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
the origin as a hull vertex. The same aggregate now separates the source
algorithm route: CKYZ local-surface extraction is not applicable because the
visible phase is not a local CY3 surface geometry, while a CCIT-style
twisted/vector-bundle I-function remains a candidate only after specifying the
insertions and qN history. The report now computes the basic virtual-dimension
handoff for that candidate: the visible CY5 phase needs a complex codimension
`2` insertion to produce a numerical CY3-style invariant, the natural visible
candidate is `H^2` on the weighted base, and the fractional
`H^2 = 1/2` pairing means the insertion cannot be promoted without
stack/source normalization and orbifold-sector pairing data. Cyrus also now
profiles the first CCIT hypergeometric factor ranges for the selected bundle
degrees `[1,1,2]`: half-degree terms land in the twisted sector and require
orbifold-sector pairing, while integer untwisted terms have numerator zero
order `3`, exceeding the codimension-`2` insertion and therefore requiring
equivariant/source normalization rather than a naive non-equivariant `H^2`
coefficient readout.

The same report now also aggregates the Chen-Ruan source-map blockers for the
CCIT-style route. The top-level counts expose two separate pairings: the
adjacent canonical `K_P(1,1,2)` reference reads the `p^2` one-point correlator
through the dual basis class `2*lambda*fun_0-8*p`, while the actual rank-three
split bundle uses the inverse-Euler twisted pairing and reads the untwisted
`p^2` insertion through
`2*(lambda-p)^2*(lambda-2*p)*fun_0`. The split-bundle ordinary
non-equivariant readout still vanishes before mirror-map extraction, and any
nonzero split-bundle contribution still requires a sourced big-J/pairing
reconstruction before numerical promotion. These aggregates are report-level
bookkeeping, not a promotion rule.

The source-basis facts behind that report are now owned by
`cyrus_core::local_orbifold::weighted_p2_rank_three_split_bundle_chen_ruan_source_basis_readout`
rather than being hard-coded only in `mcallister_gv_context`. The core readout
pins the `K_P(1,1,2)` Chen-Ruan basis, the adjacent canonical class
`2*lambda*fun_0-8*p`, and the split-bundle twisted dual class
`2*(lambda-p)^2*(lambda-2*p)*fun_0 =
2*lambda^3*fun_0-8*lambda^2*p+10*lambda*p^2`, with normalization divisor
`2*lambda^3`. It is still source metadata, not a numerical GV promotion rule.
The context report now also exports
`local_cygv_source_resolution_star_union_current_chamber_unresolved_generator_weighted_p2_twisted_vector_bundle_ifunction_chen_ruan_stack_pairing_normalization_status_counts`
so this blocker is visible at the aggregate level.
The dimension, first-Chern, required-insertion, and missing-input readiness
calculation is likewise owned by
`cyrus_core::local_orbifold::weighted_p2_rank_three_split_bundle_source_readiness`,
so the report no longer has to derive this source boundary locally.

Cyrus now traces that dual-basis readout at the numerator level. Since CCIT's
small `J`-function is written as
`<phi_epsilon/(z(z-psi))> phi^epsilon`, the `p^2` insertion coefficient is
the coefficient of `phi^2`, not the ordinary-basis `p^2` component. In the
`K_P(1,1,2)` Chen-Ruan basis,
`phi^2 = 2*lambda*fun_0 - 8*p`. For the selected split bundle the CCIT local
inverse-Euler pairing replaces that adjacent canonical divisor by
`2*lambda^3`, so the numerator-only diagnostic reads the candidate `p^2`
insertion from the ordinary `fun_0` coefficient divided by `2*lambda^3`. For
the first integer split-bundle sectors this gives finite non-equivariant
lambda-polynomials: `d=1` gives `[-1/2,1/2]` and `d=2` gives
`[-3,23/2,-17,12,-4,1/2]`.

The diagnostic also divides by the weighted-base denominator constants from
the same integer-sector hypergeometric factor. For `d=1` and `d=2`, those
constants are `2` and `96`, giving pre-mirror readouts
`[0,0,-1/4,1/4]` and
`[0,0,-1/32,23/192,-17/96,1/8,-1/24,1/192]`.

Cyrus now also expands the full weighted-base denominator through `p^2`
instead of only using its constant term. For the first two integer sectors,
the denominator polynomials are

```text
d=1: D(p) = 2 + 10p + 18p^2
d=2: D(p) = 96 + 688p + 2072p^2
```

The truncated quotient changes the ordinary-basis `p` and `p^2` coefficients,
but the dual-basis `p^2` readout still comes from the `fun_0` coefficient
divided by the split-bundle `2*lambda^3` normalization, so it remains the same
finite readout above. This rules out the raw `p^2` coefficient, denominator
`p`-terms, and the adjacent canonical `2*lambda` normalization as the missing
extraction, but still does not promote a GV value: the nonzero terms are not at
the primary extraction layer.

The same full-denominator diagnostic now keeps the `z` order instead of
collapsing to a pre-mirror lambda polynomial. For the first two integer sectors
the small-J primary `p^2` readout at `z^-2` vanishes:

```text
d=1: primary z^-2 p^2 readout = 0; first nonzero at z^-3 is -1/4
d=2: primary z^-2 p^2 readout = 0; first nonzero at z^-3 is -1/32
```

This is the remaining obstruction after the split-bundle pairing normalization:
the first descendant terms have a non-equivariant limit, but the naive primary
one-point `p^2` insertion is zero in these sectors. A numerical contribution
therefore needs a source-derived twisted big-J/pairing reconstruction or
equivalent qN history; the small-J primary readout still does not provide a
promotion rule.

This ordinary-sector `z`-profile calculation now lives in
`cyrus-core::local_orbifold` as reusable exact rational code rather than only
inside the McAllister context binary. The core API exposes both a sparse
dual-basis lambda-polynomial map and a materialized readout profile that records
the missing `z^-2` primary as an explicit zero.

The wider rank-three split-bundle degree profile now lives there as well:
`cyrus_core::local_orbifold::weighted_p2_rank_three_twisted_ifunction_degree_profiles`
owns the exact numerator, denominator, full hypergeometric, primary/descendant
`z`-readout, and scalar mirror-map blocker statuses consumed by
`mcallister_gv_context`. The binary keeps only legacy test-only helper wrappers
for low-level regression coverage.

The source-certificate boundary now lives in core as well:
`cyrus_core::local_orbifold::weighted_p2_rank_three_twisted_ifunction_source_certificate_requirements`
combines the source-readiness facts, Chen-Ruan source basis, and checked
degree profiles into a single promotion gate. For the selected
`O(-1)+O(-1)+O(-2) -> P(1,1,2)` phase it records that CCIT supplies a
Lagrangian-cone/I-function handoff, not a numerical CY3 invariant: the checked
integer sectors have zero primary `z^-2` `p^2` readout, the first nonzero
signals are descendant-layer data with a finite non-equivariant limit, and
promotion remains blocked on a source chamber certificate, qN history,
codimension-two CY3 projection/local model, and big-J/pairing reconstruction for
descendant readout. The McAllister context report now serializes that
certificate per unresolved generator and
aggregates the certificate promotion, primary-readout, pairing/residue,
twisted-dual-pairing, and required-input statuses.

A follow-up read of the source paper cited by CCIT makes the promotion gate
more precise. `I(q,z) in L` is only the cone-membership step. To extract
twisted invariants one must first identify a twisted `J`-function by putting
the `I`-function in `F(q) z 1 + G(q) + O(z^-1)` form and inverting the mirror
map; the local examples then extract non-descendant data from the `z^-1` term
using the twisted dual pairing. The core certificate now keeps those as
separate blockers: split-bundle `J` normalization, mirror-map inversion,
`z^-1` extraction, and stack-normalized twisted dual pairing. This prevents a
descendant/equivariant coefficient from being mistaken for a numerical CY3 GV
row.

The split-bundle `J`-normalization gate is now computed rather than assumed.
For each checked degree Cyrus bounds the leading hypergeometric inverse
`z`-power by
`sum(base denominator factors) - (sum(bundle numerator factors) - zero factors)`.
For the selected `[1,1,2]` bundle this minimum is `3`, hence after the outer
leading `z` every positive-degree term starts at inverse power at least `2`.
The checked positive-degree terms therefore cannot modify the CCIT `F z + G`
normalization or mirror map. The remaining obstruction is sharper: the
positive-degree `z^-1` primary layer is absent, while the first nonzero
positive-degree terms are descendant/equivariant and still require the
twisted pairing/big-`J`/qN history before any numerical contribution can be
used.

The z-order convention is checked against the adjacent canonical
`K_P(1,1,2)` first sector by using bundle degree `[4]` in the same
hypergeometric readout. There the primary `z^-2` coefficient is nonzero and
equals `11*lambda/4`, matching the first canonical benchmark before mirror-map
subtleties enter at higher degree. The split-bundle zero is therefore a real
feature of the visible `O(-1)+O(-1)+O(-2)` numerator, not an artifact of reading
the wrong z coefficient.

Cyrus now also has a stronger adjacent-source regression for the canonical
model. The test `canonical_weighted_p2_mirror_map_matches_ccit_b_table`
implements the CCIT mirror map `q = x exp(4 h(x))` and recovers the published
`K_P(1,1,2)` table
`b_d = 11/4, 525/16, 6152/9, 1146765/64, 53305261/100, 51550873/3` for
`d <= 6` from the ordinary-`p` primary coefficient. The companion test
`canonical_weighted_p2_mirror_map_matches_ccit_c_table` now also checks the
twisted half-sector primary readout and reproduces
`c_d = -2, -52/9, -2002/25, -83004/49, -3554552/81, -154984300/121,
-6835086702/169` through half-degree `13/2`. The companion
`split_bundle_kp112_mirror_map_diagnostic_has_zero_primary_p_signal` applies
the same adjacent diagnostic to bundle degrees `[1,1,2]` and still finds zero
primary ordinary-`p` signal through the checked integer sectors; the companion
`split_bundle_kp112_mirror_map_diagnostic_has_zero_half_sector_primary_signal`
does the same for the half-sector primary through half-degree `7/2` and also
finds `[0, 0, 0, 0]`. Thus the canonical mirror-map/pairing benchmark is now
pinned in Rust for both untwisted and twisted primary sectors, and it
reinforces rather than removes the split-bundle primary obstruction. The same
computed diagnostics are now serialized in the target-filtered reports. For
both targets `7` and `8`, the adjacent-canonical status counts are
`kp112_canonical_mirror_map_reproduces_ccit_b_table_to_degree_6` and
`kp112_canonical_mirror_map_reproduces_ccit_c_table_to_half_degree_13_over_2`,
while the split-bundle status counts are
`rank_three_split_kp112_mirror_map_primary_p_signal_zero_to_degree_4` and
`rank_three_split_kp112_mirror_map_half_sector_primary_signal_zero_to_half_degree_7_over_2`;
the sample source map carries the canonical `b_d`/`c_d` values and split
primary coefficients `[0, 0, 0, 0]` in both checked sectors. The split
half-sector is not empty: the first nonzero dual-basis readout starts one
descendant level later, at `z^-3`, with coefficients
`[2, -322/27, -11744/375, -22221448/25725]` through half-degree `7/2`.
That keeps the open problem sharply on twisted big-`J`/pairing
reconstruction, rather than on a missing primary mirror-map coefficient.

The source-map report now also serializes the full checked half-sector
descendant profile, not only the first nonzero term. Each readout records the
inverse `z` power, the raw lambda polynomial, the dual-basis-normalized
polynomial when division by `2*lambda` is valid, and a normalization status.
This preserves the `z^-2` zero primary rows and the later descendant rows in
the report, so downstream source-history work can inspect the exact structured
evidence without re-running or reinterpreting the lossy first-nonzero summary.
The report also aggregates the profile statuses and nonzero inverse-`z` powers:
for the checked split-bundle source map this gives five normalized-zero
readouts and seven normalized-nonzero readouts, with four nonzero readouts at
`z^-3` and three additional checked descendant readouts at `z^-4`.

This does not mean the `z^-3` descendant term is useless. CCIT explicitly
states that one-point descendants can reconstruct genus-zero Gromov-Witten
data in the semipositive/generated-by-degree-two setting. But that is a
reconstruction theorem for the relevant Givental cone/source theory, not a
rule that lets Cyrus read a descendant coefficient as a numerical CY3 GV row.
For the present six-point row the missing data is still the certified source
model plus reconstruction/pairing/qN history that says how this descendant
enters the compact corrected-chamber invariant.

A fresh read of the CCIT source makes that boundary more concrete. Proposition
`stuff` says the Givental cone is reconstructed from the small `J`-function
only when the source is semipositive and Chen-Ruan cohomology is generated by
degree-two classes; the broader statement needs the big `J`-function restricted
to a specified generating degree-two subspace. The line-bundle-sum theorem then
computes total-space invariants as inverse-Euler-twisted invariants of the
base, using a hypergeometric modification of the base small `J`-function. That
is a source construction, not a coefficient-promotion rule. The adjacent
`K_P(1,1,2)` example is exactly of this form: it uses the canonical bundle
degree `[4]`, mirror map `q = x exp(4 h(x))`, and reports
`<p^2>_{0,1,1} = 11 lambda / 4`. Our split bundle has degrees `[1,1,2]` and a
zero primary readout, so the canonical benchmark verifies the z convention but
does not supply the missing reconstruction for the split-bundle row.

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
cygv-0.1.2/src/hkty.rs
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
The regression `cygv_explicit_semigroup_rejects_fractional_intersection_tensor`
pins this boundary: fractional `H^2` is rejected before HKTY handoff as a
non-integral intersection tensor.

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
rank-three weighted-`P2` row. The local context diagnostics now surface this
as a direct qN-trace readiness blocker on both the local CICY candidate and
the unresolved current-chamber generator summary: the current six-point row
has `complete_intersection_qn_trace_blocked_origin_exclusion_total_degree_nonzero` and
ready count `0`, so no complete-intersection `cygv` call is legitimate yet.
The newly exported degree totals make the obstruction algebraic rather than
just enumerative: the CYTools-valid origin-excluded support has total degree
`[-1]`, the origin contributes `[1]`, and the origin-included total is `[0]`.
Thus zero-degree splits appear only after violating the CYTools rule that nef
partition parts must exclude the origin.

Fresh target-filtered context reports confirm that this is not a unit-test-only
artifact. Running

```bash
./target/debug/mcallister_gv_context \
  --context /tmp/cyrus_corrected_chamber_gv_context_induced_face_smoke.json \
  --target-index 7 \
  --skip-star-union-secondary-certificates \
  --out /tmp/cyrus_gv_context_target7_z_readout_report.json

./target/debug/mcallister_gv_context \
  --context /tmp/cyrus_corrected_chamber_gv_context_induced_face_smoke.json \
  --target-index 8 \
  --skip-star-union-secondary-certificates \
  --out /tmp/cyrus_gv_context_target8_z_readout_report.json
```

gives the same z-order profile on both rows: two half-degree sectors are marked
as having no untwisted full-hypergeometric z-readout, and the two integer
sectors have zero primary `z^-2` dual-basis `p^2` readout. After applying the
split-bundle inverse-Euler pairing normalization, the first nonzero terms are
descendants at `z^-3`, with lambda-polynomial samples `[-1/4]` and `[-1/32]`.
The serialized report now records their lambda order as `0` and labels them
`weighted_p2_rank_three_split_first_nonzero_descendant_has_nonequivariant_limit_requires_reconstruction`,
so the report distinguishes "a finite descendant term exists" from "a primary
non-equivariant GV contribution has been extracted." This is also a top-level
aggregate now: both target reports have first-nonzero lambda-order counts
`{"0": 2, "missing": 2}` and non-equivariant-limit status counts with two
twisted-sector non-applicable rows plus two finite descendant rows requiring
reconstruction.
The lightweight command keeps the semigroup generator/decomposition readout but
deliberately leaves the chamber secondary certificates skipped. Deeper
complete-intersection source certificates, support-overlap diagnostics, and
lower-seed diamond qN traces remain opt-in via their existing report flags,
because they are certification paths, not prerequisites for this z-readout
check.

Cyrus also now records the scalar mirror-map obstruction explicitly. For both
target-filtered reports, the half-degree sectors are marked as twisted-sector
non-applicable, while the first two integer sectors report
`weighted_p2_rank_three_split_scalar_mirror_map_primary_zero_preserved_by_scalar_reparametrization_to_this_degree`.
This only rules out scalar degree reparametrization as the missing promotion
step; it does not rule out a genuinely sourced big-`J`, pairing, residue, or
qN-history reconstruction.

The path-support report now also aggregates the candidate qN-source semigroup
probes that explain, but do not certify, the parent-domain qN coupling. The
top-level counts include raw qN-term semigroup status, seed-expanded semigroup
status, offset-generator status, offset-generator supporting-face certificate
status, LP full/aggregate statuses, and raw `FIND_GV=false` nonintegral GW
candidate-count buckets. These fields make it mechanically visible when an
offset-generator domain matches a qN shape but still lacks a supporting-face
certificate, so it remains diagnostic rather than promotable compact chamber
history.
