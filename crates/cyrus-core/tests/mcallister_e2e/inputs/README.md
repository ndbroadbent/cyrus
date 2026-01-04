# McAllister E2E Test Inputs

Source: arXiv:2107.09064
Polytope: 4-214-647 (h11=214, h21=4)

## Primitives

- **polytope.json** - 294 primal lattice points (we compute dual)
- **flux.json** - K, M vectors (4-dim for h21=4)
- **target_volumes.json** - c_i values (6=O7, 1=D3) - physical setup choice

## Injections (McAllister-specific)

These reproduce their exact results. Our GA will compute these differently.

- **heights.json** - Triangulation heights (219 values)
- **divisor_basis.json** - Divisor basis indices (CYTools 2021)
- **kklt_basis.json** - Indices of divisors contributing to superpotential
- **curves.json** - 344 curve classes + GV invariants
