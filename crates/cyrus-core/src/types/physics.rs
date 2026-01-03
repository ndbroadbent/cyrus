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

/// Divisor volume τ_i = (1/2) κ_ijk t^j t^k (positive for valid geometry).
pub type DivisorVolume = F64<Pos>;

/// Small cycle modulus τ_s in LVS (positive).
pub type SmallCycleModulus = F64<Pos>;

// ============================================================================
// Coupling Constants
// ============================================================================

/// String coupling g_s (positive, typically 0 < g_s < 1).
pub type StringCoupling = F64<Pos>;

/// Superpotential W₀ (can be any sign, but |W₀| > 0 for stabilization).
pub type Superpotential = F64<Finite>;

/// Imaginary part of the axio-dilaton Im(τ) (positive).
pub type ImTau = F64<Pos>;

// ============================================================================
// Vacuum Energy
// ============================================================================

/// Scalar potential V₀ at the minimum (negative for AdS).
pub type ScalarPotential = F64<Neg>;

/// e^K₀ factor (positive).
pub type EK0 = F64<Pos>;

/// LVS potential value (can be any sign during optimization).
pub type LvsPotential = F64<Finite>;

/// α' correction parameter ξ (positive for LVS, requires χ < 0).
pub type XiCorrection = F64<Pos>;

// ============================================================================
// Racetrack Parameters
// ============================================================================

/// Racetrack coefficient (M·q) × N_q (can be any sign).
pub type RacetrackCoefficient = F64<Finite>;

/// Racetrack exponent q·p (can be any sign, but sorted positive for dominant terms).
pub type RacetrackExponent = F64<Finite>;

/// Numerical error delta (F-term residual, non-negative).
pub type ResidualError = F64<NonNeg>;

// ============================================================================
// KKLT Parameters
// ============================================================================

/// c_τ parameter relating g_s to W₀: c_τ = 2π / (g_s × ln(W₀⁻¹)).
pub type CTau = F64<Pos>;

/// Relative error in numerical solution (non-negative).
pub type RelativeError = F64<NonNeg>;

// ============================================================================
// Cosmology Parameters
// ============================================================================

/// Matter density parameter Ω_m (0 < Ω < 1 for physical cosmology).
pub type OmegaMatter = F64<Pos>;

/// Dark energy density parameter Ω_DE (0 < Ω < 1 for physical cosmology).
pub type OmegaDarkEnergy = F64<Pos>;

/// Hubble parameter H₀ (positive).
pub type HubbleParameter = F64<Pos>;

/// Redshift z (positive for lookback, z = 0 is today).
pub type Redshift = F64<Pos>;

/// Equation of state parameter w = P/ρ (typically -1 ≤ w ≤ 1 for dark energy).
pub type EquationOfState = F64<Finite>;

/// Scalar field value φ (can be any sign).
pub type ScalarField = F64<Finite>;

/// Scalar field velocity dφ/dN (can be any sign).
pub type ScalarFieldVelocity = F64<Finite>;

// ============================================================================
// GV Invariants
// ============================================================================

/// Gopakumar-Vafa invariant value N_q (can be any sign in general).
pub type GvValue = F64<Finite>;
