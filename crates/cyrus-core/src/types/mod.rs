//! Type-safe wrappers that encode invariants at compile time.
//!
//! These types make invalid states unrepresentable, eliminating
//! scattered runtime checks and hard-to-test error handling.
//!
//! ## Numeric Wrappers
//!
//! - [`tags`]: Tag types for invariant tracking (`Pos`, `Neg`, `NonZero`, `Finite`, etc.)
//! - [`algebra`]: Type-level arithmetic rules (`Pos + Pos = Pos`, `Pos * Finite = Finite`, etc.)
//! - [`f64`]: Type-safe f64 wrapper (`F64<Pos>`, `F64<Neg>`, `F64<Finite>`)
//! - [`i32`]: Type-safe i32 wrapper (`I32<Pos>`, `I32<Neg>`)
//! - [`i64`]: Type-safe i64 wrapper (`I64<Pos>`, `I64<Neg>`)
//! - [`rational`]: Type-safe malachite Rational wrapper (`Rational<Pos>`, etc.)
//!
//! ## Other Types
//!
//! - [`branded`]: Dimension-safe moduli vectors
//! - [`physics`]: Physics-specific type aliases (H11, H21, Volume, etc.)

// Numeric type system
pub mod algebra;
pub mod f64;
pub mod f64_macros;
pub mod i32;
pub mod i64;
pub mod ops_impl;
pub mod physics;
pub mod range;
pub mod rational;
pub mod tags;

// Utility modules
pub mod branded;

// Re-export tags
pub use tags::{
    Finite, IsFinite, IsMinusOne, IsNegative, IsNonNeg, IsNonPos, IsNonZero, IsOne, IsPositive,
    IsZero, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Zero,
};

// Re-export numeric wrappers
pub use f64::F64;
pub use i32::I32;
pub use i64::I64;
pub use rational::Rational;

// Re-export from utility modules
pub use branded::{Dim, Moduli};

// Re-export physics types
pub use physics::{
    BBHLCorrection,
    // KKLT
    CTau,
    DivisorVolume,
    EK0,
    EquationOfState,
    // GV invariants
    GvValue,
    // Hodge numbers
    H11,
    H21,
    HubbleParameter,
    ImTau,
    LvsPotential,
    OmegaDarkEnergy,
    // Cosmology
    OmegaMatter,
    // Racetrack
    RacetrackCoefficient,
    RacetrackExponent,
    Redshift,
    RelativeError,
    ResidualError,
    ScalarField,
    ScalarFieldVelocity,
    // Vacuum energy
    ScalarPotential,
    SmallCycleModulus,
    // Coupling constants
    StringCoupling,
    Superpotential,
    // Volumes and moduli
    Volume,
    VolumeClassical,
    XiCorrection,
};
