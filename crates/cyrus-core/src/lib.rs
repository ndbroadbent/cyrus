//! Core mathematical primitives for Calabi-Yau manifold computations.
//!
//! This crate provides the foundational algorithms for:
//! - Lattice polytope operations
//! - Triangulations (FRST)
//! - Intersection number computation (`κ_ijk`)
//! - Kähler cone and volume calculations

pub mod error;
pub mod lattice;
pub mod polytope;

pub use error::{Error, Result};
pub use lattice::Point;
pub use polytope::Polytope;
