//! Computational geometry primitives with exact arithmetic.
//!
//! Implements convex hull and related operations for polytope computations.
//! Uses malachite::Integer for exact arithmetic at decision boundaries.
//!
//! Reference: [[project_docs/clean_room/GEOMETRIC_PRIMITIVES.md]]
//! Reference: [[project_docs/clean_room/CONVEX_HULL_D.md]]

mod convex_hull;
mod facet;
mod hyperplane;
mod orientation;

pub use convex_hull::ConvexHull;
pub use facet::Facet;
pub use hyperplane::Hyperplane;
pub use orientation::{Orientation, orientation};
