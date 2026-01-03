//! Physics-specific type aliases for Calabi-Yau computations.
//!
//! These types encode the physical constraints of string theory:
//! - Hodge numbers have specific bounds
//! - Coupling constants have sign constraints
//! - Volumes must be positive
//!
//! Using these aliases makes function signatures self-documenting
//! and ensures physics constraints are enforced at compile time.

use super::f64::F64;
use super::i32::I32;
use super::tags::{Finite, GTEOne, Neg, NonNeg, Pos};

// ============================================================================
// Hodge Numbers
// ============================================================================

/// Hodge number h¹¹ - counts Kähler moduli.
///
/// For Calabi-Yau threefolds: h¹¹ ≥ 1 (at least one Kähler modulus).
pub type H11 = I32<GTEOne>;

/// Hodge number h²¹ - counts complex structure moduli.
///
/// For Calabi-Yau threefolds: h²¹ ≥ 0 (can be zero for rigid CY).
pub type H21 = I32<NonNeg>;

// ============================================================================
// Volumes and Moduli
// ============================================================================

/// Calabi-Yau volume in string frame (positive for valid physics).
pub type Volume = F64<Pos>;

/// Classical volume before α' corrections (always positive).
pub type VolumeClassical = F64<Pos>;

/// BBHL α' correction (can be positive, negative, or zero).
pub type BBHLCorrection = F64<Finite>;

// ============================================================================
// Coupling Constants
// ============================================================================

/// String coupling g_s (positive).
pub type StringCoupling = F64<Pos>;

/// Superpotential W₀ (can be any sign, but |W₀| > 0 for stabilization).
pub type Superpotential = F64<Finite>;

// ============================================================================
// Vacuum Energy
// ============================================================================

/// Scalar potential V₀ at the minimum (negative for AdS).
pub type ScalarPotential = F64<Neg>;

/// e^K₀ factor (positive).
pub type EK0 = F64<Pos>;
