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
- `corrected_heights.dat`: corrected-chamber FRST checkpoint for diagnostics.
  It may be used to compare the chamber implied by the ancillary data after
  Cyrus has computed its own chamber data, but not to select a production
  chamber.
- `corrected_target_volumes.dat`: corrected-chamber classical KKLT divisor
  volumes. This is a checkpoint used to localize the remaining instanton-layer
  discrepancy, not an input to the production solve.

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

## Open Gaps

This policy is not a completion certificate. The remaining known gaps are:

- the high-dimensional selected-curve pruning is still pair-sum production
  pruning plus bounded diagnostics, not a full faithful Hilbert-basis or
  arbitrary "sums of others" implementation;
- corrected-chamber/flop continuation for the Kähler-coordinate instanton terms
  remains unresolved, especially for even-parity branch-cut cases;
- corrected-chamber GV coverage still has a checkpoint-implied target residual;
- the final McAllister corrected volume is close but not reproduced exactly by
  the no-replay path.

Until those gaps are resolved, McAllister 4-214-647 is a validation target with
documented discrepancies, not a completed reproduction.
