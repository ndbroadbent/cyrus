//! # cyrus-core
//!
//! Core mathematical primitives for Calabi-Yau manifold computations.
//!
//! This crate provides:
//! - Lattice polytope operations
//! - Triangulations (FRST)
//! - Intersection number computation (κ_ijk)
//! - Kähler cone and volume calculations
//! - Gopakumar-Vafa invariants
//!
//! ## Implementation Philosophy
//!
//! All algorithms are implemented from first principles based on the mathematical
//! literature, not ported from existing codebases. Key references:
//!
//! - Cox, Little, Schenck: "Toric Varieties" (textbook)
//! - Kreuzer, Skarke: "Complete classification of reflexive polyhedra" (hep-th/0002240)
//! - Demirtas et al.: "Computational Mirror Symmetry" (arXiv:2303.00757)
//!
//! ## Example
//!
//! ```rust,ignore
//! use cyrus_core::{Polytope, Triangulation};
//!
//! let polytope = Polytope::from_vertices(&vertices)?;
//! let triangulation = polytope.triangulate()?;
//! let kappa = triangulation.intersection_numbers();
//! ```

#![warn(clippy::pedantic)]
#![warn(missing_docs)]
#![deny(unsafe_code)]

pub mod error;
pub mod lattice;
pub mod polytope;
// pub mod triangulation;  // TODO
// pub mod intersection;   // TODO
// pub mod kahler;         // TODO
// pub mod gv;             // TODO

pub use error::{Error, Result};
pub use lattice::LatticePoint;
pub use polytope::Polytope;
