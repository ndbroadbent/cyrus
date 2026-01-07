# Intersection Number Computation Bug Investigation

## Goal
Compute intersection numbers κ_{ijk} for the dual CY manifold (h21=4) that match CYTools' ground truth.

**Target values (from CYTools):**
- κ[3,4,5] = 1
- κ[1,2,3] = 6
- κ[1,1,1] = 56
- κ[3,3,3] = 1

**Our values (wrong):**
- κ[3,4,5] = -122/365 ≈ -0.334
- κ[1,2,3] = 344/315 ≈ 1.09
- κ[1,1,1] = 3421/355 ≈ 9.64
- κ[3,3,3] = -187/401 ≈ -0.467

## What We Know Works

### 1. Points match CYTools exactly
```
Triangulation points (0-8):
  0 : [0, 0, 0, 0]      <- origin
  1 : [-1, 2, -1, -1]
  2 : [1, -1, 0, 0]
  3 : [-1, -1, 1, 1]
  4 : [-1, -1, 1, 2]
  5 : [-1, -1, 2, 1]
  6 : [-1, -1, 2, 3]
  7 : [-1, -1, 3, 2]
  8 : [-1, -1, 2, 2]
```

### 2. Triangulation matches CYTools exactly
15 simplices, each containing origin (index 0):
```
[0,1,2,3,4], [0,1,2,3,5], [0,1,2,4,6], [0,1,2,5,7], [0,1,2,6,7],
[0,1,3,4,5], [0,1,4,5,8], [0,1,4,6,8], [0,1,5,7,8], [0,1,6,7,8],
[0,2,3,4,5], [0,2,4,5,8], [0,2,4,6,8], [0,2,5,7,8], [0,2,6,7,8]
```

### 3. GLSM matches CYTools exactly
```
Our GLSM (4 rows × 9 columns):
  0: [-6,  2,  3, -1,  1,  1,  0,  0,  0]
  1: [-6,  2,  3,  1, -1,  0,  1,  0,  0]
  2: [-12, 4,  6,  0,  1,  0,  0,  1,  0]
  3: [-6,  2,  3,  0,  0,  0,  0,  0,  1]

CYTools GLSM:
  0: [-6, 2, 3, -1, 1, 1, 0, 0, 0]   <- identical
```

### 4. Simplex determinants computed correctly
```
Simplex 5:  rays [1,3,4,5], |det| = 3  → κ^V_{1345} = 1/3
Simplex 10: rays [2,3,4,5], |det| = 2  → κ^V_{2345} = 1/2
```
This is a **singular** (non-smooth) toric variety because |det| > 1 for some simplices.

## The Algorithm We're Implementing

Based on clean room documentation and user's analysis of CYTools:

1. **Known 4-form values**: For distinct indices (i,j,k,l) in a simplex:
   ```
   κ^V_{ijkl} = 1/|det(v_i, v_j, v_k, v_l)|
   ```

2. **Unknown variables**: Self-intersection 4-tuples with repeated indices:
   - {a,a,b,c}, {a,a,b,b}, {a,a,a,b}, {a,a,a,a}

3. **GLSM constraints**: For each GLSM row and each 3-face probe (j,k,l):
   ```
   Σ_i Q_i κ^V_{ijkl} = 0
   ```

4. **Origin handling**: Using linear equivalence D_0 ~ -Σ_{m>0} D_m:
   ```
   Σ_{i>0} (Q_i - Q_0) κ^V_{ijkl} = 0
   ```

5. **Solve**: Normal equations (least squares) for overdetermined system

6. **Sum to get 3-forms**:
   ```
   κ_{ijk} = Σ_l κ^V_{ijkl}
   ```

## The Fundamental Problem

For κ_{345}, we have:
```
κ_{345} = κ^V_{1345} + κ^V_{2345} + κ^V_{3345} + κ^V_{3445} + κ^V_{3455}
        = 1/3 + 1/2 + (self-intersection sum S)
        = 5/6 + S
```

For κ_{345} = 1, we need **S = 1/6 ≈ 0.167**.

But the GLSM equations say otherwise. For probe (3,4,5), GLSM row 3:
```
GLSM row 3: [-6, 2, 3, 0, 0, 0, 0, 0, 1]

Effective coefficients (Q_i - Q_0):
  Q_1 - Q_0 = 2 - (-6) = 8
  Q_2 - Q_0 = 3 - (-6) = 9
  Q_3 - Q_0 = 0 - (-6) = 6
  Q_4 - Q_0 = 0 - (-6) = 6
  Q_5 - Q_0 = 0 - (-6) = 6
  Q_8 - Q_0 = 1 - (-6) = 7  (but κ^V_{3458} = 0, not in any simplex)

Equation:
  8*(1/3) + 9*(1/2) + 6*(κ^V_{3345} + κ^V_{3445} + κ^V_{3455}) = 0
  8/3 + 9/2 + 6*S = 0
  43/6 + 6*S = 0
  S = -43/36 ≈ -1.19
```

**This directly contradicts the expected value S = 1/6!**

The GLSM equation says S ≈ -1.19, but physics requires S ≈ 0.17. This is a ~7× discrepancy in magnitude and opposite sign.

## What I've Tried

### 1. Fixed GLSM column/point mapping
- Initially had 10 GLSM columns for 9 points (double-counting origin)
- Fixed by filtering origin from points before GLSM computation

### 2. Verified effective coefficient formula
- Using (Q_i - Q_0) to handle origin linear equivalence
- Mathematically derived this from D_0 ~ -Σ_{m>0} D_m

### 3. Verified known 4-form values
- 1/|det| for distinct tuples
- Determinants match expected values

### 4. Verified variable enumeration
- Only self-intersection tuples (repeated indices) are variables
- Distinct tuples are known values moved to RHS

### 5. Verified 3-form summation
- Sum over ALL rays l, not per-simplex

## What I Don't Understand

### 1. Why do GLSM equations contradict expected values?

The GLSM row 3 equation **directly** says S = -43/36, but CYTools produces κ_{345} = 1 which requires S = 1/6. This is not a numerical precision issue - it's a fundamental contradiction.

**Question**: Is there something fundamentally different about how CYTools constructs these equations?

### 2. Is the "effective coefficient" approach correct?

I derived (Q_i - Q_0) from the origin linear equivalence. But maybe CYTools handles this differently?

**Question**: Does CYTools use a different formulation that doesn't involve effective coefficients?

### 3. Are there additional constraints I'm missing?

Maybe the GLSM equations alone don't determine the intersection numbers. Are there:
- Additional linear relations?
- Positivity constraints?
- Topological constraints?

**Question**: What other constraints does CYTools use beyond GLSM?

### 4. Is the 3-form summation correct?

I compute κ_{ijk} = Σ_l κ^V_{ijkl}. But maybe:
- The sum includes origin?
- The sum has weights?
- There's a different formula entirely?

**Question**: Is the relationship between 3-forms and 4-forms exactly κ_{ijk} = Σ_l κ^V_{ijkl}?

### 5. Why does CYTools get correct results?

CYTools uses "Normal Equations" (M^T M x = M^T C) for least squares. But if the equations are inconsistent (as they appear to be), least squares would give a minimum-residual solution, not the correct answer.

**Question**: How does CYTools get exact integer results (κ_{345} = 1) from an apparently inconsistent system?

## Specific Debug Output

For probe (3,4,5), GLSM row 0:
```
[DEBUG] Q_0 = -6
[DEBUG]   m=1: (2 - -6) * κ^V[1,3,4,5]=0.3333 -> RHS -= 2.6667
[DEBUG]   m=2: (3 - -6) * κ^V[2,3,4,5]=0.5000 -> RHS -= 4.5000
[DEBUG]   m=3: (-1 - -6) * κ^V[3,3,4,5][var#120]
[DEBUG]   m=4: (1 - -6) * κ^V[3,4,4,5][var#123]
[DEBUG]   m=5: (1 - -6) * κ^V[3,4,5,5][var#124]
[DEBUG]   m=6,7,8: not in simplex (coefficient 6 but value 0)
[DEBUG]   => RHS = -7.166667
[DEBUG]   => 3 variable terms
```

System statistics:
- 336 equations, 167 variables
- Residual: 4.56, RHS norm: 133
- Relative residual: 3.4%

## Files

- Implementation: `crates/cyrus-core/src/intersection/compute.rs`
- Test: `crates/cyrus-core/tests/mcallister_e2e/stage3_intersection.rs`
- Ground truth: `tests/mcallister_e2e/assertions/dual_intersection_ground_truth.json`
- Clean room docs: `project_docs/clean_room/INTERSECTION_NUMBERS.md`

## Request

I need help understanding:

1. **The exact algorithm CYTools uses** - not just the high-level description, but the precise mathematical formulation
2. **Why the GLSM equations appear inconsistent** with the expected results
3. **What I'm fundamentally misunderstanding** about intersection number computation on singular toric varieties
