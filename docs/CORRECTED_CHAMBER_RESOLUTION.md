# Corrected-Chamber Resolution (4-214-647)

Status: the no-replay corrected `V_string` residual (~0.0718 against
`corrected_cy_vol.dat`) is **resolved**. The cause was a single bug in the
B-field parity vector `gamma`, not missing GV invariants. This note records
the data-level derivation, the certificates, and what changed in Cyrus.

## Summary of findings

1. **Chamber relation (integer-exact certificate).** McAllister's corrected
   chamber (`corrected_heights.dat` FRST, 1003 simplices) is the input chamber
   (`heights.dat` FRST, 1011 simplices) after **10 conifold flops**. Each
   flopped class is a single isolated GV=1 curve, and all ten are rows of
   `small_curves.dat` (rows 1, 3, 9, 16, 17, 22, 24, 25, 32, 33). The
   prime-divisor intersection tensors satisfy
   `kappa_input - kappa_corrected = sum_C (q_C)^{x3}` exactly on all 198
   differing entries (verified in CYTools-latest; zero mismatches; no
   origin-row entries change). Script:
   `string_theory/mcallister_2107/latest_cytools/identify_corrected_chamber_flop.py`
   and `decompose_kappa_flop_diff.py`.

2. **The corrected targets and corrected volume are the input-chamber
   expressions analytically continued.** With the correct `gamma` (below),
   all ten crossed curves are **odd-parity**, so
   `Li2/Li3((-1)^{gamma.q} e^{-2pi q.t})` are real-analytic through their
   walls and the input-chamber formulas evaluate the corrected chamber
   exactly (the polylog inversion identities generate precisely the
   `kappa -= q^3` and `chi(D_i) += q_i` flop shifts; see paper
   arXiv:2107.09064 around eq. (4.11)). Numerically, at
   `t = corrected_kahler_param.dat`:
   - `corrected_target_volumes.dat` is reproduced by
     `c_i/c_tau + chi(D_i)/24 - (1/4pi^2) sum_q N_q q_i Li2(...)` with
     **integer** Euler numbers `m_i` that equal the corrected-chamber Braun
     `chi(D_i) = 12 chi(O_D) - D^3` for **all 214 divisors identically**
     (residual fractional parts < 0.012 = McAllister's own Newton tolerance;
     their classical tau at the checkpoint matches targets to ~5.5e-4).
   - `corrected_cy_vol.dat = 4711.432499235554` is reproduced to
     **1.1e-8** by the input-side evaluation
     `kappa_in t^3/6 - zeta(3)*420/(4(2pi)^3)
      + sum_q N_q [Li3 + 2pi(q.t) Li2]((-1)^{gamma.q} e^{-2pi q.t}) / (2(2pi)^3)`
     over exactly the 344 `small_curves.dat` classes (ten of them at
     negative volume, continued), and equivalently to 4.8e-8 by the
     flopped-side evaluation (corrected kappa', ten classes negated).
   Scripts: `explain_corrected_targets.py`,
   `solve_missing_gv_by_integrality.py`, `final_corrected_chamber_solve.py`.

3. **The GV curve set needs no corrected-chamber re-enumeration.** The
   shipped 344-class input-chamber selection, flop-continued, is the complete
   set McAllister used. The famous "9 missing corrected-chamber GV targets"
   contribute only ~2e-4 to any tau target and ~6e-4 to V (volumes 0.71-1.0
   at the corrected point); they are artifacts of re-enumerating toric curves
   in a (slightly different) chamber and are **not needed** for the
   reproduction. The weighted-P2/orbifold-GW program for computing them is
   unnecessary for this purpose.

4. **The bug: `gamma` (O7 B-field parity) was restricted to the KKLT basis.**
   The orientifold's O7-planes sit on **51** prime toric divisors: the 49
   declared `c_i = 6` so(8) stacks in `kklt_basis.dat` **plus points 2 and
   46**, which are not in the KKLT basis (the four non-KKLT prime divisors
   are {1, 2, 46, 130}). All 51 O7 points are exactly the lattice points with
   parity class `p mod 2 == (1,0,0,0)` — the involution's fixed parity — and
   no other point has it. So `gamma` is **derivable first-principles**: take
   the unique common parity `sigma` of the declared so(8) divisors, then mark
   every lattice point with parity `sigma`. With the old (49-entry) gamma,
   two of the ten crossed curves had even parity: their instanton signs were
   wrong, the KKLT fixed point was displaced (landing in a *different*
   adjacent chamber: 6 near-wall classes flip between Cyrus's solved chamber
   and McAllister's), and the corrected volume came out 0.0718 too high.

5. **Cross-validation.** On the 324 classes where the flop-mapped input GV
   table overlaps Cyrus's corrected-chamber toric formulas, the GV values
   agree on every class. The two methods certify each other; the flop map
   also fixes the GV values of any class either method misses.

## What changed in Cyrus

- `mcallister_first_principles.rs::compute_b_field_gamma_for_o7_divisors`
  now derives the involution parity `sigma` from the declared so(8)
  divisors, marks **all** lattice points with that parity (including
  non-KKLT divisors), and errors loudly if the declared `c_i` conflict with
  the parity classes. For 4-214-647 it reports
  `O7 divisors=51 (beyond KKLT basis: [2, 46])`.
- No other production change was required: the existing odd-parity polylog
  continuation (`gv_dilog_from_curve_volume_checked`,
  `real_trilog_real_axis`), the GV-corrected KKLT fixed-point iteration, and
  the input-chamber volume correction were already correct.

## Follow-ups

- The corrected-chamber re-enumeration diagnostics (toric_covered /
  toric_missing, weighted-P2 / Chen-Ruan / split-bundle probes in
  `mcallister_gv_context`) are now known to be unnecessary for the
  reproduction; they remain diagnostics and candidates for removal.
- The pinned no-replay results in
  `stage0_first_principles_runner_accepts_declared_inputs_only_data_dir`
  must track the corrected values.
- `compare_against_dat` corrected-volume tolerance is tightened accordingly.
