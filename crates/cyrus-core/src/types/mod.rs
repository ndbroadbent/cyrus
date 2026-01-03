//! Type-safe wrappers that encode invariants at compile time.
//!
//! These types make invalid states unrepresentable, eliminating
//! scattered runtime checks and hard-to-test error handling.
//!
//! ## Numeric Wrappers
//!
//! - [`tags`]: Tag types for invariant tracking (`Pos`, `Neg`, `NonZero`, `Finite`, etc.)
//! - [`arithmetic`]: Type-level arithmetic rules (`Pos + Pos = Pos`, etc.)
//! - [`f64`]: Type-safe f64 wrapper (`F64<Pos>`, `F64<Neg>`, `F64<Finite>`)
//! - [`i32`]: Type-safe i32 wrapper (`I32<Pos>`, `I32<Neg>`)
//! - [`i64`]: Type-safe i64 wrapper (`I64<Pos>`, `I64<Neg>`)
//! - [`rational`]: Type-safe malachite Rational wrapper (`Rational<Pos>`, etc.)
//!
//! ## Other Types
//!
//! - [`dimension`]: Dimension-safe vectors and matrices with HRTB branding
//! - [`range`]: Compile-time validated range types (`CheckedRange`)
//! - [`branded`]: Dimension-safe moduli vectors
//! - [`interval`]: Bounded interval types (`UnitInterval`)
//! - [`vec`]: Non-empty collection types (`NonEmptyVec<T>`)

// Numeric type system
pub mod algebra;
pub mod arithmetic;
pub mod f64;
pub mod f64_macros;
pub mod i32;
pub mod i64;
pub mod ops_impl;
pub mod physics;
pub mod rational;
pub mod tags;

// Other type modules
pub mod dimension;
pub mod range;

// Utility modules
pub mod branded;
pub mod interval;
pub mod vec;

// Re-export tags
pub use tags::{
    Finite, IsFinite, IsMinusOne, IsNegative, IsNonNeg, IsNonPos, IsNonZero, IsOne, IsPositive,
    IsZero, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Zero,
};

// Re-export arithmetic rules
pub use arithmetic::{AbsResult, AddResult, DivResult, MulResult, NegResult, SubResult};

// Re-export numeric wrappers
pub use f64::F64;
pub use i32::I32;
pub use i64::I64;
pub use rational::Rational;

// Re-export from dimension module
pub use dimension::{BMat, BSquareMat, BVec, Brand, ColDim, ColIdx, Dim, Index, RowDim, RowIdx};

// Re-export from range module
pub use range::CheckedRange;

// Re-export from utility modules
pub use branded::{Dim as LegacyDim, Moduli};
pub use interval::UnitInterval;
pub use vec::NonEmptyVec;

// Re-export physics types
pub use physics::{H11, H21, Volume, VolumeClassical, BBHLCorrection, StringCoupling, Superpotential, ScalarPotential, EK0};
