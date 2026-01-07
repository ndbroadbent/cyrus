# Intersection Number Algorithm - Help Needed

## Problem Statement

I cannot reproduce CYTools' intersection numbers. Unit tests confirm all inputs are correct, but the GLSM equations don't yield the expected values.

**Target**: κ_{345} = 1
**Computed**: κ_{345} ≈ -0.33 (Method A) or inconsistent system (Method B)

## Verified Components (All Correct)

### 1. GLSM Matrix ✅
```
Row 0: [-6,  2,  3, -1,  1,  1,  0,  0,  0]
Row 1: [-6,  2,  3,  1, -1,  0,  1,  0,  0]
Row 2: [-12, 4,  6,  0,  1,  0,  0,  1,  0]
Row 3: [-6,  2,  3,  0,  0,  0,  0,  0,  1]
```
Column 0 = origin, Columns 1-8 = points 1-8. Matches CYTools exactly.

### 2. Simplex Determinants ✅
```
Simplex 5:  [0,1,3,4,5] → rays [1,3,4,5], |det| = 3
Simplex 10: [0,2,3,4,5] → rays [2,3,4,5], |det| = 2
```

### 3. Known 4-Form Values ✅
```
κ^V_{1,3,4,5} = 1/|det| = 1/3
κ^V_{2,3,4,5} = 1/|det| = 1/2
```

### 4. Required Self-Intersection Sum
```
κ_{345} = κ^V_{1345} + κ^V_{2345} + κ^V_{3345} + κ^V_{3445} + κ^V_{3455}
       = 1/3 + 1/2 + S
       = 5/6 + S

For κ_{345} = 1:  S = 1/6 ≈ 0.167
```

## The Problem: GLSM Equations

### Method A: Effective Coefficients (Q_m - Q_0)

For probe (3,4,5), each GLSM row gives an equation:

| Row | Equation | If S = x+y+z |
|-----|----------|--------------|
| 0 | 5x + 7y + 7z = -7.17 | - |
| 1 | 7x + 5y + 6z = -7.17 | - |
| 2 | 12x + 13y + 12z = -14.33 | - |
| 3 | 6x + 6y + 6z = -7.17 | **S = -1.19** ❌ |

Row 3 directly says `6S = -7.17`, so **S = -1.19**, not 1/6.

### Method B: Raw Q_m (Exclude Origin Column)

| Row | Equation | Problem |
|-----|----------|---------|
| 0 | -x + y + z = -2.17 | OK |
| 1 | x - y = -2.17 | OK |
| 2 | y = -4.33 | OK |
| 3 | **0 = -2.17** | **INCONSISTENT!** |

Row 3 has no variables (Q_3=Q_4=Q_5=0) but RHS ≠ 0.

Solving Rows 0-2:
- Row 2: y = -4.33
- Row 1: x = y - 2.17 = -6.5
- Row 0: z = -2.17 + x - y = -4.34
- **S = x + y + z = -15.17** ❌

## Questions

1. **What is the "homogeneity row"?**

   The documentation says to remove "the GLSM charge matrix row for 'homogeneity' (all 1s)". But none of our GLSM rows are all 1s. All rows sum to 0 (CY condition), but that's different.

2. **How does CYTools handle the inconsistency?**

   Method B's Row 3 is `0 = -2.17`, which is impossible. How does CYTools get κ_{345} = 1 from this?

3. **Is there a different equation formulation?**

   Maybe the GLSM equations I'm constructing are fundamentally wrong. What's the correct relationship between GLSM charges and intersection numbers?

4. **Are there additional constraints?**

   Perhaps there are constraints beyond GLSM (e.g., from topology or positivity) that I'm missing?

## Unit Test Code

The analysis above comes from `tests/intersection_unit_tests.rs`:
- `test_glsm_with_origin_matches_cytools` - verifies GLSM ✅
- `test_simplex_determinants` - verifies |det| values ✅
- `test_known_4form_values` - verifies 1/|det| ✅
- `test_all_glsm_rows_for_probe_345` - shows the equation problem ❌

## Request

I need clarification on the exact algorithm CYTools uses. Specifically:
1. How are equations constructed from GLSM?
2. What rows/columns are removed and why?
3. How is the inconsistency in Method B resolved?
