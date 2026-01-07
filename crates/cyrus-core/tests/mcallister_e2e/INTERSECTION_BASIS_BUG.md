# Intersection Number Basis Transformation Bug

## Discovery Date: 2026-01-05

## Summary

Our `intersection_in_basis` function is **fundamentally wrong**. It extracts a subtensor by selecting entries where all indices are in the basis, but the correct approach is to **transform** the intersection numbers using the GLSM linear relations.

## The Bug

Current implementation (`basis.rs:262-279`):
```rust
pub fn intersection_in_basis(kappa: &Intersection, basis: &[usize]) -> Intersection {
    // This is WRONG - it just extracts entries
    for (&(i, j, k), val) in kappa.iter() {
        let bi = basis.iter().position(|&b| b == i);
        let bj = basis.iter().position(|&b| b == j);
        let bk = basis.iter().position(|&b| b == k);
        if let (Some(bi), Some(bj), Some(bk)) = (bi, bj, bk) {
            kappa_basis.set(bi, bj, bk, val.clone());
        }
    }
}
```

## Why This Is Wrong

For the dual polytope (h21=4), using McAllister's basis [3,4,5,8]:
- The full tensor has 30 non-zero κ_ijk entries (indices 1-8)
- Extracting subtensor for [3,4,5,8] gives only **2 non-zero entries**:
  - κ[3,4,5] = 5/6
  - κ[4,5,8] = 5/6
- This makes the N matrix **SINGULAR** - cannot compute flat direction!

## The Correct Approach

The GLSM defines linear relations among divisors:
```
Σ_i Q^a_i D_i ~ 0  for each GLSM row a
```

This means each divisor D_i can be expressed as a linear combination of basis divisors:
```
D_i = Σ_a T_{ia} D_a  where {D_a} is the basis
```

The intersection numbers must be **transformed**:
```
κ_abc = Σ_{ijk} T_{ai} T_{bj} T_{ck} κ_{ijk}
```

This is what CYTools does when you call `intersection_numbers(in_basis=True)`.

## Experimental Verification

Python script tested 70 different 4-element bases from indices [1-8]:
- **None** produced e^K₀ = 0.234393 (McAllister's value)
- Best result: [1,3,4,5] → e^K₀ = 5.54 (24x off!)
- Basis [3,4,5,8] → N matrix is singular (no solution)

This proves that simple subtensor extraction cannot reproduce McAllister's results.

## Required Fix

1. Compute the transformation matrix T from the GLSM kernel
2. Implement proper tensor transformation: κ_abc = Σ T_ai T_bj T_ck κ_ijk
3. Update `intersection_in_basis` to use this transformation

## Impact

This affects ALL intersection number computations that use `intersection_in_basis`:
- Volume calculations
- Flat direction computation (e^K₀)
- Kähler moduli stabilization

## References

- CYTools paper (arXiv:2211.03823) - describes the basis transformation
- McAllister paper (arXiv:2107.09064) - our validation target

## Investigation Progress (2026-01-05)

### Fixed: GLSM Double-Origin Bug

We discovered the GLSM was being computed with a doubled origin:
- When `include_origin=true`, an origin column is prepended
- But McAllister's points already include origin at index 0
- This caused degenerate GLSM relations and garbage intersection numbers

**Fix applied**: For McAllister's data, compute GLSM for non-origin points only (indices 1-8)
with `include_origin=false`, then use `compute_intersection_numbers_with_offset(offset=1)`.

### Still Wrong: Basis Transformation

Even after fixing the GLSM, the intersection numbers are wrong:
- Our solve gives values like κ[0,1,2] = 130/123 ≈ 1.06
- Expected: κ[0,1,2] = κ_points[3,4,5] = 5/6 ≈ 0.833
- The linear solve is populating self-intersections with wrong values

**Root Cause**: The direct formula only gives non-zero values for triplets that appear in
simplex faces. For basis [3,4,5,8], only 2 triplets (κ[3,4,5] and κ[4,5,8]) have direct
values. Other entries (self-intersections) must come from the GLSM transformation, not
from a least-squares solve.

## Next Steps

### Immediate Actions

1. **Request clean room documentation** for how CYTools computes `intersection_numbers(in_basis=True)` - specifically the GLSM-based transformation

2. **Implement the transformation** once we have the algorithm spec

### Algorithm Sketch (Needs Verification)

The transformation should work as follows:

Given:
- Full intersection tensor κ_{ijk} (point-indexed)
- GLSM charge matrix Q (shape: h+1 × n_pts)
- Basis indices B = {b_1, ..., b_h}

Compute transformation matrix T (shape: n_pts × h):
1. Partition columns of Q into basis (Q_B) and non-basis (Q_N)
2. Compute T_non_basis = -Q_N^T @ (Q_B^T)^{-1} (or similar)
3. T_basis = identity

Then:
κ_{abc} = Σ_{ijk} T_{a,i} T_{b,j} T_{c,k} κ_{ijk}

This is a rank-3 tensor contraction, which is O(h³ × n_pts³) but can be optimized.

### Test Plan

After implementing:
1. Verify e^K₀ = 0.234393 for McAllister's basis [3,4,5,8]
2. Verify that intersection numbers are **integers** (or simple fractions) after transformation
3. Verify N matrix is non-singular and p is reasonable
