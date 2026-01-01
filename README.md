# Cyrus

**High-performance Calabi-Yau manifold computations for string theory research.**

Cyrus is a Rust toolkit for evaluating Calabi-Yau compactifications, computing moduli stabilization, and searching the string theory landscape for configurations matching observed physics.

## Crates

| Crate | Description |
|-------|-------------|
| `cyrus-core` | Lattice polytopes, triangulations, intersection numbers, Kähler geometry |
| `cyrus-moduli` | KKLT and LVS moduli stabilization, racetrack mechanism |
| `cyrus-cosmology` | Friedmann equations, quintessence evolution, CPL fitting |
| `cyrus-cache` | Intelligent tiered caching with self-optimization |
| `cyrus-ga` | Genetic algorithm for landscape search |
| `cyrus-verify` | Formally verified core algorithms (via Aeneas) |

## Goals

1. **Reproduce published results** - Validate against McAllister et al. (arXiv:2107.09064) and Cicoli et al. (arXiv:2407.03405)
2. **Performance** - 10-100× faster than Python implementations
3. **Correctness** - 100% test coverage, formal verification for critical paths
4. **Usability** - Clean API, thorough documentation

## Installation

```bash
cargo add cyrus-core cyrus-moduli cyrus-cosmology
```

## Quick Example

```rust
use cyrus_core::{Polytope, Triangulation};
use cyrus_moduli::{KkltSolver, Racetrack};

// Load a polytope from vertices
let polytope = Polytope::from_vertices(&vertices)?;

// Get a fine regular star triangulation
let triangulation = polytope.triangulate()?;

// Compute intersection numbers
let kappa = triangulation.intersection_numbers();

// Solve KKLT stabilization
let solver = KkltSolver::new(&kappa, &c_values);
let solution = solver.solve()?;

println!("V_string = {}", solution.volume);
println!("V₀ = {}", solution.vacuum_energy);
```

## Validation

```bash
# Validate against McAllister's KKLT examples
cargo run --bin cyrus-validate -- mcallister

# Validate against Cicoli's quintessence example
cargo run --bin cyrus-validate -- cicoli
```

## License

Licensed under either of:
- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
- MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## Contributing

Contributions are welcome! Please read the [architecture document](docs/ARCHITECTURE.md) first.

## Citation

If you use Cyrus in your research, please cite:

```bibtex
@software{cyrus,
  author = {Broadbent, Nathan},
  title = {Cyrus: High-Performance Calabi-Yau Computations},
  url = {https://github.com/ndbroadbent/cyrus},
  year = {2026}
}
```

## Acknowledgments

This project builds on decades of work by the string theory community. Key references:
- McAllister et al., "Minimal Flux Compactifications" (arXiv:2107.09064)
- Cicoli et al., "From Inflation to Quintessence" (arXiv:2407.03405)
- Demirtas et al., "Small Cosmological Constants" (arXiv:1912.10047)
