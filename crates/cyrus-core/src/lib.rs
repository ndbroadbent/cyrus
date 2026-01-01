//! Core mathematical primitives for Calabi-Yau manifold computations.
//!
//! This crate provides the foundational algorithms for:
//! - Lattice polytope operations
//! - Intersection number computation (`κ_ijk`)
//! - Kähler cone and volume calculations
//!
//! ## Volume Computation
//!
//! The string frame volume is computed as:
//! ```text
//! V_string = (1/6) κ_ijk t^i t^j t^k - BBHL
//! ```
//!
//! where BBHL = ζ(3) χ(X) / (4(2π)³) is the α' correction.

pub mod error;
pub mod intersection;
pub mod lattice;
pub mod polytope;
pub mod volume;

pub use error::{Error, Result};
pub use intersection::Intersection;
pub use lattice::Point;
pub use polytope::Polytope;
pub use volume::{VolumeResult, bbhl_correction, compute_volume, volume_classical, volume_string};
