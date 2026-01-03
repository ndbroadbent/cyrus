//! Type-level arithmetic algebra.
//!
//! Defines the rules for how tagged numeric types combine under arithmetic operations.
//! These rules are defined ONCE and apply to all wrapper types (F64, I32, I64, Rational).
//!
//! # Example
//!
//! ```text
//! Pos + Pos => Pos      (positive + positive = positive)
//! Pos + Neg => Finite   (could be anything)
//! Pos * Neg => Neg      (positive * negative = negative)
//! Pos / Zero => !       (division by zero doesn't compile - no impl exists)
//! ```

mod add;
mod div;
mod mul;
mod sub;
mod unary;

// ============================================================================
// Output Traits - Define what operations produce
// ============================================================================

/// What does `Self + Rhs` produce?
pub trait AddOutput<Rhs> {
    /// The resulting tag type.
    type Output;
}

/// What does `Self - Rhs` produce?
pub trait SubOutput<Rhs> {
    /// The resulting tag type.
    type Output;
}

/// What does `Self * Rhs` produce?
pub trait MulOutput<Rhs> {
    /// The resulting tag type.
    type Output;
}

/// What does `Self / Rhs` produce? (No impl for Zero divisor = compile error)
pub trait DivOutput<Rhs> {
    /// The resulting tag type.
    type Output;
}

/// What does `-Self` produce?
pub trait NegOutput {
    /// The resulting tag type.
    type Output;
}

/// What does `Self.abs()` produce?
pub trait AbsOutput {
    /// The resulting tag type.
    type Output;
}

// ============================================================================
// Rule Declaration Macros
// ============================================================================

macro_rules! add_rules {
    ($($lhs:ty, $rhs:ty => $out:ty;)*) => {
        $(impl AddOutput<$rhs> for $lhs { type Output = $out; })*
    };
}

macro_rules! sub_rules {
    ($($lhs:ty, $rhs:ty => $out:ty;)*) => {
        $(impl SubOutput<$rhs> for $lhs { type Output = $out; })*
    };
}

macro_rules! mul_rules {
    ($($lhs:ty, $rhs:ty => $out:ty;)*) => {
        $(impl MulOutput<$rhs> for $lhs { type Output = $out; })*
    };
}

macro_rules! div_rules {
    ($($lhs:ty, $rhs:ty => $out:ty;)*) => {
        $(impl DivOutput<$rhs> for $lhs { type Output = $out; })*
    };
}

macro_rules! neg_rules {
    ($($input:ty => $out:ty;)*) => {
        $(impl NegOutput for $input { type Output = $out; })*
    };
}

macro_rules! abs_rules {
    ($($input:ty => $out:ty;)*) => {
        $(impl AbsOutput for $input { type Output = $out; })*
    };
}

pub(crate) use {abs_rules, add_rules, div_rules, mul_rules, neg_rules, sub_rules};

#[cfg(test)]
#[allow(clippy::used_underscore_items)]
mod tests {
    use super::*;
    use crate::types::tags::{GTEOne, Neg, Pos};

    // Compile-time check that the rules make sense
    fn _check_add<L: AddOutput<R>, R>() {}
    #[allow(dead_code)]
    fn _check_sub<L: SubOutput<R>, R>() {}
    #[allow(dead_code)]
    fn _check_mul<L: MulOutput<R>, R>() {}
    #[allow(dead_code)]
    fn _check_div<L: DivOutput<R>, R>() {}

    #[test]
    fn test_add_rules_exist() {
        _check_add::<Pos, Pos>();
        _check_add::<Pos, Neg>();
        _check_add::<GTEOne, GTEOne>();
    }

    #[test]
    fn test_div_by_zero_not_impl() {
        // This would fail to compile:
        // _check_div::<Pos, Zero>();
    }
}
