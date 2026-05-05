# McAllister Data Policy

This note records how Cyrus may use the McAllister ancillary files while building
the GA-ready pipeline. The rule is simple: validation data can check a computed
quantity, but it must not become the way Cyrus computes that quantity.

## Declared Model Inputs

These files describe the validation target and may be loaded by the
first-principles runner:

- `points.dat`: primal polytope points for the validation target.
- `heights.dat`: selected regular triangulation heights.
- `K_vec.dat`, `M_vec.dat`: integer flux vectors.
- `kklt_basis.dat`: declared non-perturbative divisor set.
- `target_volumes.dat`: declared `c_i` values from the orientifold/KKLT model
  choice.

For a generic GA search, these inputs are replaced by the candidate polytope,
triangulation, flux vector, divisor selection, and orientifold data. They are
not downstream physics outputs.

## Validation Checkpoints

These files may be read by tests or comparison code after Cyrus computes the
corresponding object:

- `dual_points.dat`, `dual_simplices.dat`: dual polytope and dual FRST
  checkpoints.
- `basis.dat`: primal divisor basis checkpoint.
- `dual_curves.dat`, `dual_curves_gv.dat`: low-dimensional mirror/racetrack
  curve and GV checkpoints.
- `small_curves.dat`, `small_curves_cutoff.dat`, `small_curves_gv.dat`:
  selected high-dimensional primal toric-curve checkpoint. These validate the
  selected toric-curve method, not a full 214-dimensional generic
  `compute_gvs()` call.
- `g_s.dat`, `W_0.dat`: racetrack outputs.
- `c_tau.dat`: derived KKLT scalar checkpoint computed from the racetrack
  outputs, not a declared model input.
- `corrected_heights.dat`: corrected-chamber FRST checkpoint for diagnostics.
  It may be used to compare the chamber implied by the ancillary data after
  Cyrus has computed its own chamber data, but not to select a production
  chamber.
- `corrected_target_volumes.dat`: corrected-chamber classical KKLT divisor
  volumes. This is a checkpoint used to localize the remaining instanton-layer
  discrepancy, not an input to the production solve.
- `potent_rays.dat`, `potent_rays_gv.dat`, `potent_rays_rank.dat`,
  `potent_rays_vols.dat`: potent-ray validation samples and derived
  convergence quantities. Cyrus may use them as assertions after reconstructing
  the ray contexts, but not as the source of a reusable search algorithm.

The reusable runner must still compute these quantities from upstream inputs
when they are needed for a new landscape candidate.

## Replay-Only Artifacts

These are downstream outputs and must not be used to make the first-principles
run pass:

- `kahler_param.dat`
- `corrected_kahler_param.dat`
- `cy_vol.dat`
- `corrected_cy_vol.dat`
- `small_curves_vols.dat`

`corrected_kahler_param.dat` is only allowed as an explicit replay path behind
`--allow-downstream-kahler`, or as a comparison/diagnostic checkpoint after the
first-principles Kähler vector has already been computed. It must not seed the
normal KKLT branch search or fixed-point solve.

`corrected_cy_vol.dat` may be read by `compare_against_dat` after Cyrus has
computed `V_string`. The current comparison is intentionally not exact: the
runner logs the residual as an unresolved instanton/chamber discrepancy.

## Current Runner Boundary

`mcallister_first_principles` is expected to run with only the declared inputs
above. Stage 0 has a heavy opt-in test that copies only those files into a
temporary directory and verifies that the runner reaches the current no-replay
result.

The runner also has validation diagnostics that read downstream checkpoints
when they are present:

- corrected Kähler vector comparison;
- corrected target-volume comparison;
- checkpoint-`t` corrected-chamber GV target diagnostics.

Those diagnostics print `[COMPARE]` output and are not promoted into the
production `t` or `V_string` result. This distinction is important: diagnostic
reads are allowed only because they explain discrepancies; they are not allowed
to repair them.

## Legacy Python Boundary

The older Python reproduction scripts are useful context, but they are not the
contract for the GA-ready engine.

- `string_theory/mcallister_2107/full_pipeline_from_data.py` computes the
  dual-side flat direction, `e^{K0}`, and a racetrack from `dual_curves*.dat`,
  but then reads `cy_vol.dat` for `V[0]` with a TODO to derive it from Kähler
  stabilization. Despite its success message, that script is a validation
  checkpoint, not a no-replay pipeline.
- `string_theory_project/research/mcallister_reproduction/REPRODUCTION_OUTLINE.md`
  still describes high-dimensional `small_curves*.dat` as if it came from a
  generic `cy.compute_gvs(max_deg=N)` call. The paper and ancillary readme are
  more specific: the Kähler-side data comes from selected toric curves below a
  volume cutoff, followed by removal of curves that are sums of others and GV
  evaluation of the retained curves.

## Open Gaps

This policy is not a completion certificate. The remaining known gaps are:

- the high-dimensional selected-curve checkpoint is pair-pruned, not
  finite-semigroup-pruned: a Cyrus finite integer-feasibility diagnostic finds
  5 additional multi-term decompositions among the 344 retained input-chamber
  curves, so full finite-semigroup pruning would retain 339 curves instead of
  matching `small_curves.dat`. Cyrus exposes both rules explicitly via
  `CurvePruningStrategy` and
  `mcallister_first_principles --small-curve-pruning <pair|finite-semigroup>`;
- corrected-chamber/flop continuation for the Kähler-coordinate instanton terms
  remains unresolved, especially for even-parity branch-cut cases;
- corrected-chamber GV coverage still has a checkpoint-implied target residual;
- potent-ray checkpoints are now validated for rank, volumes, convergence
  slopes, and first-four source-derived CKYZ GV entries for the 395 rank-two
  CKYZ rows, but Cyrus does not yet generate the low-dimensional-face ray
  sample, reproduce all ten saved entries efficiently, or handle the rank-four
  local charge contexts;
- the final McAllister corrected volume is close but not reproduced exactly by
  the no-replay path.

Until those gaps are resolved, McAllister 4-214-647 is a validation target with
documented discrepancies, not a completed reproduction.
