# McAllister E2E Test Inputs

Source: arXiv:2107.09064
Polytope: 4-214-647 (h11=214, h21=4)

## JSON Fixtures (Deprecated)

These files are legacy JSON fixtures and are **only used when**
`CYRUS_ALLOW_FIXTURES=1` is set. First-principles runs should use
the `.dat` files from the McAllister data directory.

**Note:** `CYRUS_ALLOW_FIXTURES` must not be used with `CYRUS_FIRST_PRINCIPLES`.

- **polytope.json** - 294 primal lattice points (we compute dual)
- **heights.json** - Triangulation heights (219 values)
- **flux.json** - K, M vectors (4-dim for h21=4)
