//! KKLT vacuum verification: primal-geometry data model and (incrementally)
//! the full Kähler-moduli stabilization entry point.
//!
//! These types and the verification pipeline are being lifted out of the
//! `mcallister_first_principles` binary so the full KKLT stabilization is a
//! callable library routine - used by both the GA verifier and the CLI -
//! rather than logic trapped in a binary `main`. The binary was always
//! intermediate scaffolding; the library is the product.

use malachite::Integer;

use crate::types::f64::F64;
use crate::types::tags::Finite;
use crate::{DivisorBasis, Intersection, Point, Polytope, Triangulation};

/// Primal toric geometry: the reflexive polytope, its FRST heights, the
/// triangulation point set, and the triangulation itself.
pub struct PrimalGeom {
    /// The (canonicalized) primal reflexive polytope.
    pub polytope: Polytope,
    /// FRST heights for the triangulation points.
    pub heights: Vec<F64<Finite>>,
    /// Triangulation points (non-facet-interior lattice points).
    pub triangulation_points: Vec<Point>,
    /// The fine regular star triangulation.
    pub triangulation: Triangulation,
}

/// Primal intersection data: GLSM charges, linear relations, the divisor
/// basis, and the intersection tensor in both the full ambient and the
/// reduced basis.
pub struct PrimalIntersection {
    /// GLSM charge matrix rows (ambient divisor coordinates).
    pub glsm: Vec<Vec<Integer>>,
    /// Linear relations among the triangulation points.
    pub linrels: Vec<Vec<Integer>>,
    /// Divisor basis column indices.
    pub basis: Vec<usize>,
    /// Full ambient intersection tensor.
    pub kappa_full: Intersection,
    /// Intersection tensor reduced to the divisor basis.
    pub kappa_basis: Intersection,
}

/// An owned divisor basis (index list or matrix), mirroring the borrowed
/// [`DivisorBasis`] so it can be stored across pipeline stages.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum OwnedDivisorBasis {
    /// Ambient divisor column indices.
    Indices(Vec<usize>),
    /// Divisor-basis rows in ambient divisor coordinates.
    Matrix {
        /// Standard GLSM column basis used to reduce the matrix block.
        standard_basis: Vec<usize>,
        /// Divisor-basis rows in ambient divisor coordinates, including origin.
        basis_matrix: Vec<Vec<Integer>>,
    },
}

impl OwnedDivisorBasis {
    /// Borrow as a [`DivisorBasis`] for the library APIs.
    #[must_use]
    pub fn as_divisor_basis(&self) -> DivisorBasis<'_> {
        match self {
            Self::Indices(indices) => DivisorBasis::Indices(indices),
            Self::Matrix {
                standard_basis,
                basis_matrix,
            } => DivisorBasis::Matrix {
                standard_basis,
                basis_matrix,
            },
        }
    }

    /// Number of divisor-basis elements.
    #[must_use]
    pub const fn dimension(&self) -> usize {
        match self {
            Self::Indices(indices) => indices.len(),
            Self::Matrix { basis_matrix, .. } => basis_matrix.len(),
        }
    }

    /// Human-readable description of the basis representation.
    #[must_use]
    pub const fn description(&self) -> &'static str {
        match self {
            Self::Indices(_) => "index divisor basis",
            Self::Matrix { .. } => "matrix divisor basis",
        }
    }
}
