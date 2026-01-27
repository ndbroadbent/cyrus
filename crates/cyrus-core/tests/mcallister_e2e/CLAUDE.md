
### McAllister Validation Purpose

The McAllister paper (arXiv:2107.09064) provides ground truth for ONE specific example (polytope 4-214-647). We use this to validate that our pipeline computes correct results.

But the validation must test the SAME code path we'll use in production:
- Start with polytope points (the only true primitive)
- Compute triangulation
- Compute intersection numbers
- Compute flat direction
- Solve racetrack
- Compute W₀, V_string, V₀

If any stage is "load from file" in tests but "compute" in production, the test is worthless.

### McAllister Intermediate Data Files

**CRITICAL**: We have ALL of McAllister's intermediate data in `.dat` files, not just final results.

Location: `/Users/ndbroadbent/code/string_theory/resources/small_cc_2107.09064_source/anc/paper_data/`

Available polytopes:
- `4-214-647/` - Our primary validation target
- `5-113-4627-main/`
- `5-113-4627-alternative/`
- `5-81-3213/`
- `7-51-13590/`

Each directory contains intermediate data for every pipeline stage:
```
points.dat           - Polytope points (input)
heights.dat          - Triangulation heights
dual_simplices.dat   - Triangulation simplices
basis.dat            - Divisor basis indices (CRITICAL for matching)
K_vec.dat            - K flux vector
M_vec.dat            - M flux vector
kahler_param.dat     - Kähler parameters
g_s.dat              - String coupling
W_0.dat              - Superpotential
cy_vol.dat           - CY volume
small_curves.dat     - Curve classes for GV
small_curves_gv.dat  - GV invariant values
potent_rays.dat      - Rays for potential
... and more
```

**Use these for stage-by-stage verification**: Compute each stage, compare against the corresponding `.dat` file. This lets us pinpoint exactly where our computation diverges from theirs.

**Fixtures are opt-in only**: JSON fixtures are deprecated and only used when explicitly enabled
via `CYRUS_ALLOW_FIXTURES=1`. First-principles runs must use `.dat`.

**Mutual exclusion**: `CYRUS_ALLOW_FIXTURES` must NOT be used with `CYRUS_FIRST_PRINCIPLES`.

**CRITICAL: basis.dat vs kklt_basis.dat are DIFFERENT:**
- `basis.dat`: Snapshot of CYTools 2021's divisor basis (214 indices for h11=214)
- `kklt_basis.dat`: Indices of divisors that **contribute to the superpotential**
- `target_volumes.dat`: The **c_i values** (dual Coxeter numbers: 6 for O7, 1 for D3)

**CRITICAL: corrected_* vs uncorrected files:**
- `kahler_param.dat`: Approximate KKLT (ignoring instanton corrections)
- `corrected_kahler_param.dat`: Actual KKLT vacuum (WITH instanton corrections)
- Same for `heights.dat` vs `corrected_heights.dat`, `cy_vol.dat` vs `corrected_cy_vol.dat`

### CYTools Version Differences (CRITICAL!)

**CYTools 2021 vs latest produces DIFFERENT divisor bases for the same triangulation.**

| Version | Basis | K (flux) | M (flux) |
|---------|-------|----------|----------|
| CYTools 2021 (McAllister) | [3,4,5,8] | [-3, -5, 8, 6] | [10, 11, -11, -5] |
| CYTools 2025 (Latest) | [5,6,7,8] | [8, 5, -8, 6] | [-10, -1, 11, -5] |

**Transformation rules:**
- K transforms as **covariant**: `K_new = T⁻¹ @ K_old`
- M transforms as **contravariant**: `M_new = T.T @ M_old`

See `string_theory/mcallister_2107/LATEST_CYTOOLS_CONVERSION_RESULT.md` for details.

**Physics values (invariant under transformation):**
- e^{K₀} = 0.234393
- g_s = 0.00911134
- W₀ = 2.30 × 10⁻⁹⁰
- **V_string = 4711.83** (our validation target)
- V₀ = -5.5 × 10⁻²⁰³ Mpl⁴

### Critical Implementation Details

#### Divisor Basis (no overrides)

McAllister’s basis indices are stored in `basis.dat`. Our code must compute the divisor
basis deterministically and **match `basis.dat` exactly**. No overrides. If it diverges,
the algorithm is wrong and must be fixed.

#### Curve Discovery (GV Invariants)

Curves for GV invariants are FOUND via enumeration, not known ahead of time. McAllister found a specific set of curves. Our approach:

1. **Compute curves from first principles** - run enumeration from our code.
2. **Compare to McAllister’s data** - require that our results match their `.dat` files for this seed.
3. **If mismatch** - treat as a bug and fix the algorithm (no overrides).

This approach:
- Proves the enumeration is identical for the validation target
- Avoids any hidden substitution of precomputed curves

See research docs:
- `string_theory/mcallister_2107/CURVE_DISCREPANCY.md` - superset proof
- `string_theory/mcallister_2107/BASIS_INDEX_MISMATCH_V1.md` - basis issues

## Testing Philosophy

**Be thorough and honest. Never hide failures. Never cheat.**

### Discrepancies Are GOLD

**When you find a discrepancy between our computation and McAllister's expected values, THIS IS THE ENTIRE POINT OF WHAT WE ARE DOING.**

A discrepancy like "κ_abc p^a p^b p^c comes out negative when it should be positive" is not a problem to work around - it's a BUG to investigate. You do NOT:
- Mark the test as `#[ignore]`
- Add a comment saying "needs investigation" and move on
- Assume it's a "convention difference" without proving it
- Give up and say "we'd need their precomputed values"

You DO:
- Stop everything else
- Investigate until you understand EXACTLY why the discrepancy exists
- Write comprehensive research notes documenting your investigation
- Fix the bug in our code (if it's our bug)
- Document the convention difference (if it's a convention difference) with PROOF
- Keep the test failing until the issue is resolved

**Discrepancies reveal bugs. Bugs must be fixed. This is non-negotiable.**

If you find yourself saying "this doesn't work because X" - that's when the REAL work begins. X is the bug. Find it. Fix it. Understand it completely.

### Never Hide Test Failures

- **NEVER** use `#[ignore]` to hide failing tests
- **NEVER** comment out or delete tests that fail
- **NEVER** use `.unwrap_or()` or similar to make tests pass artificially
- Failing tests are valuable - they show what's broken
- A failing test is infinitely better than no test

### Never Cheat in Tests

- **NEVER** load precomputed values when production code computes them
- **NEVER** skip computation steps that production needs
- If a computation isn't implemented yet, the test should FAIL, not load cached data
- Tests must exercise the exact code path production will use

### Snapshot Testing for Pipeline Stages

Use insta snapshots to verify each pipeline stage independently.
Each stage computes from the previous stage's output (no loading shortcuts):

```
Stage 1:  Polytope primal points → Dual polytope lattice points
Stage 2:  (Points, Heights) → Triangulation
Stage 3:  Triangulation → Intersection numbers κ
Stage 4:  (κ, K, M) → Flat direction p, e^K₀
Stage 5:  (Triangulation, curves) → GV invariants
Stage 6:  (GV, M, p, target_volumes) → Racetrack terms
Stage 7:  Racetrack terms → (g_s, τ)
Stage 8:  (g_s, τ, terms) → W₀
Stage 9:  (τ, κ) → Kähler moduli t
Stage 10: (κ, t, h11, h21) → V_string
Stage 11: (e^K₀, g_s, V_string, W₀) → V₀
```

Benefits:
- Verified output of one stage feeds into the next
- When something breaks, you know exactly which stage failed
- Snapshots serve as documentation of expected intermediate values
- Tests the ACTUAL computation code, not file loading

### E2E Test Directory Structure

```
crates/cyrus-core/tests/mcallister_e2e/
├── inputs/              # SEARCH PARAMETERS - what we vary in GA
│   ├── polytope.json    # Which CY manifold (enumerated from Kreuzer-Skarke)
│   ├── heights.json     # Which triangulation (point in secondary fan)
│   └── flux.json        # Which flux configuration (K, M integers)
│
├── overrides/           # LEGACY OVERRIDES (deprecated)
│   │                    # First-principles runs should use .dat inputs instead.
│   ├── dual_basis.json     # Legacy dual basis indices
│   ├── curves.json         # Legacy curves for GV invariants
│   └── target_volumes.json # Legacy brane setup (c_i values)
│
├── snapshots/           # OUR COMPUTED VALUES (insta snapshots)
│   └── ...              # Each stage writes its output here
│
└── assertions/          # MCALLISTER'S GROUND TRUTH
    └── ...              # Compare our snapshots against these
```

**The Distinction:**

1. **inputs/** - Search parameters. These define the vacuum we're evaluating.
   - In the GA, we SEARCH over these (enumerate polytopes, mutate heights/flux)
   - We CANNOT derive these from other data - they are free parameters
   - McAllister found these specific values through their search

2. **overrides/** - Deprecated legacy overrides.
   - First-principles runs should NOT use these.
   - Use `.dat` files in McAllister data dir instead.
   - Overrides remain only for legacy debugging.

3. **snapshots/** - What WE compute from inputs (no overrides).
   - Each pipeline stage writes its output here
   - Insta manages these automatically

4. **assertions/** - McAllister's ground truth.
   - All intermediate values from their `.dat` files
   - Compare our snapshots against these to validate correctness

**Why snapshots/ duplicates assertions/:**

Both directories contain the same intermediate values (triangulation, intersection numbers, etc.). The difference:
- `assertions/` = what McAllister computed (ground truth)
- `snapshots/` = what WE computed (our output)

By keeping them separate, we GUARANTEE that:
1. Our pipeline actually computes everything from inputs/
2. We never accidentally load precomputed values as shortcuts
3. Tests compare our computed snapshots/ against their assertions/

If we only had assertions/, it would be tempting to read from there during computation. The separate snapshots/ directory makes cheating impossible.

**Overrides vs Assertions:**

Overrides are deprecated. If you enable fixtures, you MUST still match the
corresponding assertions exactly. Prefer `.dat` inputs in all first-principles runs.

If we use McAllister's intermediate values as input, we must reproduce their subsequent outputs exactly. This validates our math is correct.

**Without overrides (required):**

There is no override path. If our computed basis or curves differ from the `.dat` files,
we treat it as a correctness bug and fix it.

### Two Branches: With and Without Overrides

**Branch A: No Overrides (our own search)**
```
polytope → search for our own (heights, flux) → compute everything → V₀
```
- Same polytope as McAllister
- We independently search for heights/flux that produce small V₀
- Uses our own computed divisor basis, curves, etc.
- Final V₀ won't be identical, but should be similarly small
- Validates our full pipeline + search works end-to-end
- Loose assertions (same order of magnitude)
- Already demonstrated in `string_theory/` (Python prototype)

**Branch B: With Overrides (reproduce McAllister exactly)**
```
McAllister's (polytope, heights, flux) → ... → V₀
```
- Uses McAllister's exact inputs (their found heights/flux)
- Overrides divisor basis, curves with their values
- Must EXACTLY match assertions/ at every stage
- Validates our math is correct
- Strict assertions

**Why both branches?**

1. **Branch A** proves our pipeline can independently find small-V₀ vacua
2. **Branch B** proves our computations are mathematically correct by reproducing published results exactly

If Branch B fails, our math is wrong.
If Branch A fails to find small V₀, our search algorithm is broken.

### Test Structure

```rust
// GOOD - computes everything, might fail
#[test]
fn test_w0_against_mcallister() {
    let polytope = load_polytope_points("4-214-647");
    let triangulation = compute_triangulation(&polytope, &heights);
    let kappa = compute_intersection_numbers(&triangulation);
    let p = compute_flat_direction(&kappa, &k, &m);
    // ... compute everything ...
    let w0 = compute_w0(...);

    assert_eq!(w0, expected, "W₀ mismatch");
}

// BAD - loads precomputed intersection numbers
#[test]
fn test_w0_against_mcallister() {
    let kappa = load_from_file("intersection.json");  // CHEATING!
    // ...
}
```
