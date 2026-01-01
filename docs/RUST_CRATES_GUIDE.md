# Rust Crates Guide for Cyrus

This document catalogs useful Rust crates for the cyrus project, organized by functionality.

## Design Principle: Pure Rust Preferred

We strongly prefer **pure Rust** crates over C/C++ bindings for:

1. **Type & memory safety** - No FFI boundary bugs or segfaults
2. **Formal verification** - Aeneas can only translate pure Rust to Lean. C bindings are a verification dead-end.
3. **Portability** - Single `cargo build` works everywhere
4. **Robustness** - Rust's guarantees extend through the entire codebase

Only use C bindings (qhull, GMP/rug) as fallbacks when pure Rust alternatives have issues.

## Currently Used (in Cargo.toml)

| Crate | Purpose | Pure Rust? | Notes |
|-------|---------|:----------:|-------|
| `nalgebra` | Dense linear algebra | ✅ | Standard choice, good for small fixed-size matrices |
| `ndarray` | N-dimensional arrays | ✅ | NumPy-like, good for dynamic sizes |
| `num-rational` | Rational numbers | ✅ | Exact arithmetic for fractions |
| `num-integer` | Integer utilities | ✅ | GCD, LCM, etc. |
| `faer` | Fast linear algebra | ✅ | Excellent performance, verifiable |
| `malachite` | Arbitrary precision | ✅ | LGPL, derived from GMP algorithms |
| `delaunay` | d-dim triangulations | ✅ | CGAL-inspired, replaces CYTools CGAL binary |
| `lll-rs` | LLL lattice reduction | ✅ | For basis reduction |
| `good_lp` + `microlp` | Linear programming | ✅ | For Kähler cone tests |
| `ndarray_einsum_beta` | Tensor contraction | ✅ | Einstein summation |
| `quadrature` | Numerical integration | ✅ | For period computations |
| `ode_solvers` | ODE integration | ✅ | For differential equations |
| `rayon` | Parallelism | ✅ | Data parallelism made easy |

## Recommended Additions

### Convex Hull / Computational Geometry

**For height-induced triangulations (Phase 6):**

| Crate | Use Case | Notes |
|-------|----------|-------|
| **[qhull](https://docs.rs/qhull)** | Convex hull in arbitrary dimensions | Binds to qhull C library. Best for n-dim convex hulls. Computes Delaunay triangulations, Voronoi diagrams, halfspace intersections. |
| **[delaunay](https://docs.rs/delaunay)** | d-dimensional Delaunay | Pure Rust, inspired by CGAL. Handles arbitrary dimensions. |
| **[convexhull3d](https://docs.rs/convexhull3d)** | 3D convex hull | If only 3D needed |

**Recommendation:** Use `delaunay` - it's a **pure Rust** d-dimensional Delaunay library inspired by CGAL. This is ideal because:
- No C/C++ dependencies (unlike qhull which wraps a C library)
- Handles arbitrary dimensions (essential for 4D polytopes)
- Designed as a lightweight CGAL replacement for Rust

For height-induced triangulations, the algorithm is:
1. Lift points (v_i) to (v_i, h_i) in R^{d+1}
2. Compute Delaunay triangulation
3. Extract lower faces (normal has negative last component)
4. Project back to R^d

This replaces CYTools' CGAL binary (`vendor/cytools_latest/external/cgal/cgal-triangulate.cpp`) which uses `CGAL::Regular_triangulation` with `Dynamic_dimension_tag`.

```toml
delaunay = "0.4"
```

**Fallback:** If `delaunay` has issues with very high dimensions or numerical stability, `qhull = "0.3"` wraps the battle-tested qhull C library.

### Arbitrary Precision / Rational Arithmetic

**Already have `rug`, but alternatives:**

| Crate | Use Case | Notes |
|-------|----------|-------|
| **[rug](https://docs.rs/rug)** ✓ | Arbitrary precision | GMP/MPFR/MPC bindings. Best performance. |
| **[malachite](https://malachite.rs/)** | Pure Rust bigints | LGPL licensed. High performance, derived from GMP/FLINT algorithms. Pure Rust. |
| `num-bigint` | Basic bigints | Pure Rust, simpler API, slower |

**Recommendation:** Keep `rug` for performance-critical code. Consider `malachite` if pure-Rust is needed.

### Lattice Basis Reduction (LLL)

**For Smith/Hermite normal forms and lattice computations:**

| Crate | Use Case | Notes |
|-------|----------|-------|
| **[lll-rs](https://github.com/rust-crypto-labs/lll-rs)** | LLL algorithm | Implements LLL and L² algorithms. Experimental but functional. |
| **[lllreduce](https://lib.rs/crates/lllreduce)** | LLL reduction | Simpler API, drafty |
| **[feanor-math](https://docs.rs/feanor-math)** | Number theory | Includes LLL, FFT, polynomial factoring |

**Recommendation:** Start with `lll-rs` for LLL reduction. For Smith/Hermite normal forms, may need to implement ourselves or port from GAP/Julia's MatInt.jl.

```toml
lll-rs = "0.4"
```

### Linear Programming / Integer Programming

**For Kähler cone point-in-cone tests (Phase 7):**

| Crate | Use Case | Notes |
|-------|----------|-------|
| **[good_lp](https://github.com/rust-or/good_lp)** | LP/MILP modeling | Clean API, multiple solver backends |
| **[microlp](https://lib.rs/crates/microlp)** | Pure Rust LP/MILP | Branch & bound for integers. No external deps. |
| **[lpsolve](https://docs.rs/lpsolve)** | lp_solve bindings | C library bindings |

**Solver backends for good_lp:**
- `highs` - Fast, MIT licensed, C++ (recommended)
- `clarabel` - Pure Rust, no integers
- `scip` - Very fast, academic license
- `microlp` - Pure Rust with integers, slower

**Recommendation:** Use `good_lp` with `highs` backend for production, `microlp` for pure-Rust fallback.

```toml
good_lp = { version = "1.14", features = ["highs"] }
```

### Tensor Contraction

**For intersection number computations:**

| Crate | Use Case | Notes |
|-------|----------|-------|
| **[ndarray_einsum_beta](https://docs.rs/ndarray_einsum_beta)** | Einstein summation | Port of numpy einsum. General tensor contraction. |

**Recommendation:** Use for κ_ijk contractions.

```toml
ndarray_einsum_beta = "0.7"
```

### Sparse Matrices

**For large intersection tensors and Mori cone computations:**

| Crate | Use Case | Notes |
|-------|----------|-------|
| **[faer-sparse](https://lib.rs/crates/faer-sparse)** | Sparse linear algebra | Part of faer ecosystem. Cholesky, LU, QR. |
| **[sprs](https://docs.rs/sprs)** | Sparse matrices | CSR/CSC formats, triplet construction |
| **[nalgebra-sparse](https://crates.io/crates/nalgebra-sparse)** | Sparse matrices | Based on nalgebra |
| **[russell_sparse](https://lib.rs/crates/russell_sparse)** | MUMPS/UMFPACK bindings | For very large systems |

**Recommendation:** Already have `faer`, use `faer-sparse` for consistency.

### Numerical Integration

**For period computations (Phase 8):**

| Crate | Use Case | Notes |
|-------|----------|-------|
| **[quadrature](https://docs.rs/quadrature)** | Adaptive integration | Double exponential method. Fast, no allocations. |
| **[quad_gk](https://lib.rs/crates/quad_gk)** | Gauss-Kronrod | High precision (1e-14 tolerance) |
| **[gauss_quad](https://docs.rs/gauss-quad)** | Gaussian quadrature | Legendre, Laguerre, Hermite, Jacobi |
| **[Integrate](https://mtantaoui.github.io/Integrate/)** | General integration | Romberg, Simpson, Newton, Gauss methods |

**Recommendation:** Use `quadrature` for general adaptive integration, `gauss_quad` for specific Gaussian quadrature rules.

```toml
quadrature = "0.1"
gauss_quad = "0.4"
```

### Symbolic Computation

**For polynomial manipulation and potential CAS features:**

| Crate | Use Case | Notes |
|-------|----------|-------|
| **[symbolica](https://symbolica.io/)** | Computer algebra | World-class polynomial algebra. Source-available, free for hobbyists. |
| **[symbolic_polynomials](https://docs.rs/symbolic_polynomials)** | Integer polynomials | Simple polynomial manipulation |
| **[rusymbols](https://github.com/simensgreen/rusymbols)** | Symbolic math | Aims to be like SymPy |

**Recommendation:** `symbolica` if needed, but may be overkill. Consider `symbolic_polynomials` for simpler needs.

## Missing Functionality (Need to Implement)

### Smith/Hermite Normal Form
No good Rust crate exists. Options:
1. **Implement ourselves** - Port from GAP or Julia's [MatInt.jl](https://github.com/jmichel7/MatInt.jl)
2. **Use LLL** - Hermite NF can be computed via LLL algorithm
3. **C bindings** - Wrap FLINT or other C library

### Vertex/Facet Enumeration
No pure-Rust crate for polyhedral vertex/facet enumeration. Options:
1. **Wrap lrs or cdd** - C library bindings
2. **Wrap qhull** - Already have bindings
3. **Implement** - Avis-Fukuda algorithm is simple

### Picard-Fuchs / Period Computation
No existing crate. This is specialized mathematical physics. Options:
1. **Implement** - Port algorithms from literature
2. **Use precomputed** - Load from fixtures initially
3. **Numerical** - Use `quadrature` with cycle integrals

## Recommended Cargo.toml Additions

```toml
[workspace.dependencies]
# Triangulations (pure Rust, CGAL-inspired, d-dimensional)
delaunay = "0.4"

# Lattice reduction
lll-rs = "0.4"

# Linear programming
good_lp = { version = "1.14", features = ["highs"] }

# Tensor contraction
ndarray_einsum_beta = "0.7"

# Numerical integration
quadrature = "0.1"
gauss_quad = "0.4"

# Optional: qhull (C library bindings, fallback if delaunay has issues)
# qhull = "0.3"
```

## Sources

- [qhull documentation](https://docs.rs/qhull)
- [delaunay crate](https://lib.rs/crates/delaunay) - d-dimensional Delaunay inspired by CGAL
- [RGeometry](https://rgeometry.org/) - 2D computational geometry
- [convexhull3d](https://docs.rs/convexhull3d) - 3D quickhull
- [malachite](https://malachite.rs/) - Pure Rust arbitrary precision
- [lll-rs](https://github.com/rust-crypto-labs/lll-rs) - LLL algorithm implementation
- [MatInt.jl](https://github.com/jmichel7/MatInt.jl) - Julia reference for Smith/Hermite NF
- [good_lp](https://github.com/rust-or/good_lp) - LP/MILP modeling
- [microlp](https://lib.rs/crates/microlp) - Pure Rust LP with branch & bound
- [ndarray_einsum_beta](https://docs.rs/ndarray_einsum_beta) - Einstein summation
- [faer-sparse](https://lib.rs/crates/faer-sparse) - Sparse linear algebra
- [sprs](https://docs.rs/sprs) - Sparse matrices
- [quadrature](https://docs.rs/quadrature) - Double exponential integration
- [quad_gk](https://lib.rs/crates/quad_gk) - Gauss-Kronrod quadrature
- [symbolica](https://symbolica.io/) - Computer algebra system
