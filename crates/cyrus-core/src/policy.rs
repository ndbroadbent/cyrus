//! Policy types for handling physics boundaries.
//!
//! The type system catches bugs (negative κ, dimension mismatches, NaN).
//! Policy types handle physics boundaries where "invalid but useful" values
//! occur during GA exploration.
//!
//! # The Pattern
//!
//! ```rust,ignore
//! // Strict: returns None if non-positive
//! let v: Option<F64<Pos>> = volume_string::<Strict>(&kappa, &t, h11, h21);
//!
//! // GA: returns Finite, can be negative (useful for fitness)
//! let v: F64<Finite> = volume_string::<ForGA>(&kappa, &t, h11, h21);
//!
//! // Abort: panics if non-positive (for debugging)
//! let v: F64<Pos> = volume_string::<Abort>(&kappa, &t, h11, h21);
//! ```
//!
//! # When to Use
//!
//! - `Strict`: Final results, validated pipelines, after GA finds valid candidate
//! - `ForGA`: Fitness computation, parameter space exploration
//! - `Abort`: Debugging, batch processing where failure is unexpected

use crate::types::f64::F64;
use crate::types::tags::{Finite, Neg, Pos};

// ============================================================================
// Volume Policy (positive output expected)
// ============================================================================

/// Policy for computations that should produce positive values.
///
/// Examples: string frame volume, e^{K₀}
pub trait VolumePolicy {
    /// Output type after applying the policy.
    type Output;

    /// Transform a finite value according to the policy.
    fn apply(v: F64<Finite>) -> Self::Output;
}

/// Strict physics: fail if non-positive.
///
/// Use for validated pipelines and final results.
pub enum Strict {}

impl VolumePolicy for Strict {
    type Output = Option<F64<Pos>>;

    #[inline]
    fn apply(v: F64<Finite>) -> Self::Output {
        F64::<Pos>::new(v.get())
    }
}

/// GA fitness: keep the signed value (negatives are useful for optimization).
///
/// Use for fitness computation and parameter space exploration.
pub enum ForGA {}

impl VolumePolicy for ForGA {
    type Output = F64<Finite>;

    #[inline]
    fn apply(v: F64<Finite>) -> Self::Output {
        v
    }
}

/// Debug/batch mode: panic on invalid.
///
/// Use when invalid values indicate a bug in the calling code.
pub enum Abort {}

impl VolumePolicy for Abort {
    type Output = F64<Pos>;

    #[inline]
    fn apply(v: F64<Finite>) -> Self::Output {
        F64::<Pos>::new(v.get()).expect("expected positive value")
    }
}

// ============================================================================
// Vacuum Policy (negative output expected for AdS)
// ============================================================================

/// Policy for computations that should produce negative values (AdS vacua).
///
/// Examples: vacuum energy V₀
pub trait VacuumPolicy {
    /// Output type after applying the policy.
    type Output;

    /// Transform a finite value according to the policy.
    fn apply(v: F64<Finite>) -> Self::Output;
}

impl VacuumPolicy for Strict {
    type Output = Option<F64<Neg>>;

    #[inline]
    fn apply(v: F64<Finite>) -> Self::Output {
        F64::<Neg>::new(v.get())
    }
}

impl VacuumPolicy for ForGA {
    type Output = F64<Finite>;

    #[inline]
    fn apply(v: F64<Finite>) -> Self::Output {
        v
    }
}

impl VacuumPolicy for Abort {
    type Output = F64<Neg>;

    #[inline]
    fn apply(v: F64<Finite>) -> Self::Output {
        F64::<Neg>::new(v.get()).expect("expected negative value (AdS)")
    }
}

// ============================================================================
// Generic Policy (for any finite → constrained transformation)
// ============================================================================

/// Generic policy that can target any tag.
///
/// Use when you need custom constraint checking.
pub trait Policy<Target> {
    /// Output type after applying the policy.
    type Output;

    /// Transform a finite value according to the policy.
    fn apply(v: F64<Finite>) -> Self::Output;
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_volume_strict_positive() {
        let v = F64::<Finite>::new(1.5).unwrap();
        let result = <Strict as VolumePolicy>::apply(v);
        assert!(result.is_some());
        assert!((result.unwrap().get() - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_volume_strict_negative() {
        let v = F64::<Finite>::new(-0.5).unwrap();
        let result = <Strict as VolumePolicy>::apply(v);
        assert!(result.is_none());
    }

    #[test]
    fn test_volume_ga_keeps_negative() {
        let v = F64::<Finite>::new(-0.5).unwrap();
        let result = <ForGA as VolumePolicy>::apply(v);
        assert!((result.get() - (-0.5)).abs() < 1e-10);
    }

    #[test]
    fn test_vacuum_strict_negative() {
        let v = F64::<Finite>::new(-1.5).unwrap();
        let result = <Strict as VacuumPolicy>::apply(v);
        assert!(result.is_some());
    }

    #[test]
    fn test_vacuum_strict_positive_fails() {
        let v = F64::<Finite>::new(0.5).unwrap();
        let result = <Strict as VacuumPolicy>::apply(v);
        assert!(result.is_none());
    }

    #[test]
    #[should_panic(expected = "expected positive value")]
    fn test_volume_abort_panics_on_negative() {
        let v = F64::<Finite>::new(-0.5).unwrap();
        let _ = <Abort as VolumePolicy>::apply(v);
    }
}
