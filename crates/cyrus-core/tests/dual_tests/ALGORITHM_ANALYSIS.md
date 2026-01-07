# CYTools vs Our Intersection Number Algorithm

## Summary

Our dual tests confirm:
- **All basic computations match** (determinants, 4-form values, GLSM, linear relations)
- **Intersection number algorithm is fundamentally wrong** (0/83 matches)

## CYTools Algorithm (from `_construct_intnum_equations_4d`)

### 1. Variables: 4-form tuples with repeated indices

```python
# From 3-tuples (s0, s1, s2):
variable_array_1 = [(s0,s0,s1,s2), (s0,s1,s1,s2), (s0,s1,s2,s2)]

# From 2-tuples (s0, s1):
variable_array_2 = [(s0,s0,s1,s1), (s0,s0,s0,s1), (s0,s1,s1,s1)]

# Self-intersections:
variable_array_3 = [(i,i,i,i) for i in range(1, n_points)]
```

### 2. Equations: 3-form tuples

```python
# From 3-tuples:
eqn_array_1 = [(s0,s1,s2)]

# From 2-tuples:
eqn_array_2 = [(s0,s0,s1), (s0,s1,s1)]

# Self-intersections:
eqn_array_3 = [(i,i,i) for i in range(1, n_points)]
```

### 3. The constraint system

For each linear_relations row `lin` and each 3-tuple equation `(i,j,k)`:

```
Σ_m lin[m-1] * κ^V_{m,i,j,k} = 0
```

Where:
- `lin` is from `linear_relations(include_origin=False)` - shape (dim, n_points-1)
- `lin[m-1]` accesses coefficient for point m (1-indexed point → 0-indexed array)
- Known 4-form values contribute to RHS
- Unknown 4-form values contribute to matrix M

### 4. Key insight: linear_relations indexing

CYTools uses `include_origin=False`:
```python
linear_relations = self.glsm_linear_relations(include_origin=False)
# Shape: (dim, n_points-1)
# Column 0 → point 1
# Column 1 → point 2
# ...
```

So `lin[cc[0] - 1]` means: for point index `cc[0]` (1-indexed), get coefficient at `cc[0] - 1` (0-indexed).

## Our Algorithm Issues

1. **Wrong variable enumeration**: We enumerate ALL 4-tuples with replacement, not just those with repeated indices from the triangulation structure.

2. **Wrong equation construction**: We don't use the proper eqn_dict structure that maps variables to their contributing equations.

3. **Wrong linear_relations indexing**: We use `include_origin=True` and have offset confusion.

## What We Need to Fix

1. Port the exact CYTools `_construct_intnum_equations_4d` algorithm:
   - Build `variable_array` from 3-tuples and 2-tuples of simplices
   - Build `eqn_array` from 3-tuples and 2-tuples
   - Build `eqn_dict` mapping variables to equations
   - Build `c_dict` for known 4-form contributions
   - Use `linear_relations(include_origin=False)` with correct indexing

2. Fix divisor basis:
   - h11 = GLSM kernel dimension - 1 (not GLSM rows)
   - Need to handle the scaling relation

## Test Results

```
Simplex determinants: ✓ 15/15 match
Distinct 4-form values: ✓ 15/15 match
GLSM shape: ✓ 4x9 matches
GLSM kernel property: ✓ L @ Q^T = 0 verified
Linear relations structure: ✓ Identity block verified
Intersection numbers: ✗ 0/83 match (42% residual)
```
