//! Convex lattice polytope operations.
//!
//! Implements operations on reflexive polytopes for constructing
//! Calabi-Yau hypersurfaces via toric geometry.
//!
//! Reference: [[project_docs/clean_room/DUAL_POLYTOPE.md]]
//! Reference: [[project_docs/clean_room/REFLEXIVE_POLYTOPE_DUAL.md]]

mod core;
mod dual;
mod facets;

#[cfg(test)]
mod tests;

pub use self::core::Polytope;
