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

use super::tags::*;

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

// ============================================================================
// Addition Rules
// ============================================================================

add_rules! {
    // Pos + X
    Pos, Pos => Pos;
    Pos, Neg => Finite;
    Pos, Zero => Pos;
    Pos, One => GTEOne;
    Pos, Two => GTEOne;
    Pos, MinusOne => Finite;
    Pos, NonZero => Finite;
    Pos, NonNeg => Pos;
    Pos, NonPos => Finite;
    Pos, Finite => Finite;
    Pos, GTEOne => GTEOne;

    // Neg + X
    Neg, Pos => Finite;
    Neg, Neg => Neg;
    Neg, Zero => Neg;
    Neg, One => Finite;
    Neg, Two => Finite;
    Neg, MinusOne => Neg;
    Neg, NonZero => Finite;
    Neg, NonNeg => Finite;
    Neg, NonPos => Neg;
    Neg, Finite => Finite;
    Neg, GTEOne => Finite;

    // Zero + X
    Zero, Pos => Pos;
    Zero, Neg => Neg;
    Zero, Zero => Zero;
    Zero, One => One;
    Zero, Two => Two;
    Zero, MinusOne => MinusOne;
    Zero, NonZero => NonZero;
    Zero, NonNeg => NonNeg;
    Zero, NonPos => NonPos;
    Zero, Finite => Finite;
    Zero, GTEOne => GTEOne;

    // One + X
    One, Pos => GTEOne;
    One, Neg => Finite;
    One, Zero => One;
    One, One => Two;
    One, Two => GTEOne;
    One, MinusOne => Zero;
    One, NonZero => Finite;
    One, NonNeg => GTEOne;
    One, NonPos => Finite;
    One, Finite => Finite;
    One, GTEOne => GTEOne;

    // Two + X
    Two, Pos => GTEOne;
    Two, Neg => Finite;
    Two, Zero => Two;
    Two, One => GTEOne;
    Two, Two => GTEOne;
    Two, MinusOne => One;
    Two, NonZero => Finite;
    Two, NonNeg => GTEOne;
    Two, NonPos => Finite;
    Two, Finite => Finite;
    Two, GTEOne => GTEOne;

    // MinusOne + X
    MinusOne, Pos => Finite;
    MinusOne, Neg => Neg;
    MinusOne, Zero => MinusOne;
    MinusOne, One => Zero;
    MinusOne, Two => One;
    MinusOne, MinusOne => Neg;
    MinusOne, NonZero => Finite;
    MinusOne, NonNeg => Finite;
    MinusOne, NonPos => Neg;
    MinusOne, Finite => Finite;
    MinusOne, GTEOne => Finite;

    // NonZero + X
    NonZero, Pos => Finite;
    NonZero, Neg => Finite;
    NonZero, Zero => NonZero;
    NonZero, One => Finite;
    NonZero, Two => Finite;
    NonZero, MinusOne => Finite;
    NonZero, NonZero => Finite;
    NonZero, NonNeg => Finite;
    NonZero, NonPos => Finite;
    NonZero, Finite => Finite;
    NonZero, GTEOne => Finite;

    // NonNeg + X
    NonNeg, Pos => Pos;
    NonNeg, Neg => Finite;
    NonNeg, Zero => NonNeg;
    NonNeg, One => GTEOne;
    NonNeg, Two => GTEOne;
    NonNeg, MinusOne => Finite;
    NonNeg, NonZero => Finite;
    NonNeg, NonNeg => NonNeg;
    NonNeg, NonPos => Finite;
    NonNeg, Finite => Finite;
    NonNeg, GTEOne => GTEOne;

    // NonPos + X
    NonPos, Pos => Finite;
    NonPos, Neg => Neg;
    NonPos, Zero => NonPos;
    NonPos, One => Finite;
    NonPos, Two => Finite;
    NonPos, MinusOne => Neg;
    NonPos, NonZero => Finite;
    NonPos, NonNeg => Finite;
    NonPos, NonPos => NonPos;
    NonPos, Finite => Finite;
    NonPos, GTEOne => Finite;

    // Finite + X
    Finite, Pos => Finite;
    Finite, Neg => Finite;
    Finite, Zero => Finite;
    Finite, One => Finite;
    Finite, Two => Finite;
    Finite, MinusOne => Finite;
    Finite, NonZero => Finite;
    Finite, NonNeg => Finite;
    Finite, NonPos => Finite;
    Finite, Finite => Finite;
    Finite, GTEOne => Finite;

    // GTEOne + X
    GTEOne, Pos => GTEOne;
    GTEOne, Neg => Finite;
    GTEOne, Zero => GTEOne;
    GTEOne, One => GTEOne;
    GTEOne, Two => GTEOne;
    GTEOne, MinusOne => Finite;
    GTEOne, NonZero => Finite;
    GTEOne, NonNeg => GTEOne;
    GTEOne, NonPos => Finite;
    GTEOne, Finite => Finite;
    GTEOne, GTEOne => GTEOne;
}

// ============================================================================
// Subtraction Rules
// ============================================================================

sub_rules! {
    // Pos - X
    Pos, Pos => Finite;
    Pos, Neg => Pos;
    Pos, Zero => Pos;
    Pos, One => Finite;
    Pos, Two => Finite;
    Pos, MinusOne => GTEOne;
    Pos, NonZero => Finite;
    Pos, NonNeg => Finite;
    Pos, NonPos => Pos;
    Pos, Finite => Finite;
    Pos, GTEOne => Finite;

    // Neg - X
    Neg, Pos => Neg;
    Neg, Neg => Finite;
    Neg, Zero => Neg;
    Neg, One => Neg;
    Neg, Two => Neg;
    Neg, MinusOne => Finite;
    Neg, NonZero => Finite;
    Neg, NonNeg => Neg;
    Neg, NonPos => Finite;
    Neg, Finite => Finite;
    Neg, GTEOne => Neg;

    // Zero - X
    Zero, Pos => Neg;
    Zero, Neg => Pos;
    Zero, Zero => Zero;
    Zero, One => MinusOne;
    Zero, Two => Neg;
    Zero, MinusOne => One;
    Zero, NonZero => NonZero;
    Zero, NonNeg => NonPos;
    Zero, NonPos => NonNeg;
    Zero, Finite => Finite;
    Zero, GTEOne => Neg;

    // One - X
    One, Pos => Finite;
    One, Neg => GTEOne;
    One, Zero => One;
    One, One => Zero;
    One, Two => MinusOne;
    One, MinusOne => Two;
    One, NonZero => Finite;
    One, NonNeg => Finite;
    One, NonPos => GTEOne;
    One, Finite => Finite;
    One, GTEOne => Finite;

    // Two - X
    Two, Pos => Finite;
    Two, Neg => GTEOne;
    Two, Zero => Two;
    Two, One => One;
    Two, Two => Zero;
    Two, MinusOne => GTEOne;
    Two, NonZero => Finite;
    Two, NonNeg => Finite;
    Two, NonPos => GTEOne;
    Two, Finite => Finite;
    Two, GTEOne => Finite;

    // MinusOne - X
    MinusOne, Pos => Neg;
    MinusOne, Neg => Finite;
    MinusOne, Zero => MinusOne;
    MinusOne, One => Neg;
    MinusOne, Two => Neg;
    MinusOne, MinusOne => Zero;
    MinusOne, NonZero => Finite;
    MinusOne, NonNeg => Neg;
    MinusOne, NonPos => Finite;
    MinusOne, Finite => Finite;
    MinusOne, GTEOne => Neg;

    // NonZero - X
    NonZero, Pos => Finite;
    NonZero, Neg => Finite;
    NonZero, Zero => NonZero;
    NonZero, One => Finite;
    NonZero, Two => Finite;
    NonZero, MinusOne => Finite;
    NonZero, NonZero => Finite;
    NonZero, NonNeg => Finite;
    NonZero, NonPos => Finite;
    NonZero, Finite => Finite;
    NonZero, GTEOne => Finite;

    // NonNeg - X
    NonNeg, Pos => Finite;
    NonNeg, Neg => Pos;
    NonNeg, Zero => NonNeg;
    NonNeg, One => Finite;
    NonNeg, Two => Finite;
    NonNeg, MinusOne => GTEOne;
    NonNeg, NonZero => Finite;
    NonNeg, NonNeg => Finite;
    NonNeg, NonPos => NonNeg;
    NonNeg, Finite => Finite;
    NonNeg, GTEOne => Finite;

    // NonPos - X
    NonPos, Pos => Neg;
    NonPos, Neg => Finite;
    NonPos, Zero => NonPos;
    NonPos, One => Neg;
    NonPos, Two => Neg;
    NonPos, MinusOne => Finite;
    NonPos, NonZero => Finite;
    NonPos, NonNeg => NonPos;
    NonPos, NonPos => Finite;
    NonPos, Finite => Finite;
    NonPos, GTEOne => Neg;

    // Finite - X
    Finite, Pos => Finite;
    Finite, Neg => Finite;
    Finite, Zero => Finite;
    Finite, One => Finite;
    Finite, Two => Finite;
    Finite, MinusOne => Finite;
    Finite, NonZero => Finite;
    Finite, NonNeg => Finite;
    Finite, NonPos => Finite;
    Finite, Finite => Finite;
    Finite, GTEOne => Finite;

    // GTEOne - X
    GTEOne, Pos => Finite;
    GTEOne, Neg => GTEOne;
    GTEOne, Zero => GTEOne;
    GTEOne, One => Finite;
    GTEOne, Two => Finite;
    GTEOne, MinusOne => GTEOne;
    GTEOne, NonZero => Finite;
    GTEOne, NonNeg => Finite;
    GTEOne, NonPos => GTEOne;
    GTEOne, Finite => Finite;
    GTEOne, GTEOne => Finite;
}

// ============================================================================
// Multiplication Rules
// ============================================================================

mul_rules! {
    // Pos * X
    Pos, Pos => Pos;
    Pos, Neg => Neg;
    Pos, Zero => Zero;
    Pos, One => Pos;
    Pos, Two => Pos;
    Pos, MinusOne => Neg;
    Pos, NonZero => NonZero;
    Pos, NonNeg => NonNeg;
    Pos, NonPos => NonPos;
    Pos, Finite => Finite;
    Pos, GTEOne => Pos;

    // Neg * X
    Neg, Pos => Neg;
    Neg, Neg => Pos;
    Neg, Zero => Zero;
    Neg, One => Neg;
    Neg, Two => Neg;
    Neg, MinusOne => Pos;
    Neg, NonZero => NonZero;
    Neg, NonNeg => NonPos;
    Neg, NonPos => NonNeg;
    Neg, Finite => Finite;
    Neg, GTEOne => Neg;

    // Zero * X
    Zero, Pos => Zero;
    Zero, Neg => Zero;
    Zero, Zero => Zero;
    Zero, One => Zero;
    Zero, Two => Zero;
    Zero, MinusOne => Zero;
    Zero, NonZero => Zero;
    Zero, NonNeg => Zero;
    Zero, NonPos => Zero;
    Zero, Finite => Zero;
    Zero, GTEOne => Zero;

    // One * X (identity)
    One, Pos => Pos;
    One, Neg => Neg;
    One, Zero => Zero;
    One, One => One;
    One, Two => Two;
    One, MinusOne => MinusOne;
    One, NonZero => NonZero;
    One, NonNeg => NonNeg;
    One, NonPos => NonPos;
    One, Finite => Finite;
    One, GTEOne => GTEOne;

    // Two * X
    Two, Pos => Pos;
    Two, Neg => Neg;
    Two, Zero => Zero;
    Two, One => Two;
    Two, Two => Pos;
    Two, MinusOne => Neg;
    Two, NonZero => NonZero;
    Two, NonNeg => NonNeg;
    Two, NonPos => NonPos;
    Two, Finite => Finite;
    Two, GTEOne => GTEOne;

    // MinusOne * X
    MinusOne, Pos => Neg;
    MinusOne, Neg => Pos;
    MinusOne, Zero => Zero;
    MinusOne, One => MinusOne;
    MinusOne, Two => Neg;
    MinusOne, MinusOne => One;
    MinusOne, NonZero => NonZero;
    MinusOne, NonNeg => NonPos;
    MinusOne, NonPos => NonNeg;
    MinusOne, Finite => Finite;
    MinusOne, GTEOne => Neg;

    // NonZero * X
    NonZero, Pos => NonZero;
    NonZero, Neg => NonZero;
    NonZero, Zero => Zero;
    NonZero, One => NonZero;
    NonZero, Two => NonZero;
    NonZero, MinusOne => NonZero;
    NonZero, NonZero => NonZero;
    NonZero, NonNeg => Finite;
    NonZero, NonPos => Finite;
    NonZero, Finite => Finite;
    NonZero, GTEOne => NonZero;

    // NonNeg * X
    NonNeg, Pos => NonNeg;
    NonNeg, Neg => NonPos;
    NonNeg, Zero => Zero;
    NonNeg, One => NonNeg;
    NonNeg, Two => NonNeg;
    NonNeg, MinusOne => NonPos;
    NonNeg, NonZero => Finite;
    NonNeg, NonNeg => NonNeg;
    NonNeg, NonPos => NonPos;
    NonNeg, Finite => Finite;
    NonNeg, GTEOne => NonNeg;

    // NonPos * X
    NonPos, Pos => NonPos;
    NonPos, Neg => NonNeg;
    NonPos, Zero => Zero;
    NonPos, One => NonPos;
    NonPos, Two => NonPos;
    NonPos, MinusOne => NonNeg;
    NonPos, NonZero => Finite;
    NonPos, NonNeg => NonPos;
    NonPos, NonPos => NonNeg;
    NonPos, Finite => Finite;
    NonPos, GTEOne => NonPos;

    // Finite * X
    Finite, Pos => Finite;
    Finite, Neg => Finite;
    Finite, Zero => Zero;
    Finite, One => Finite;
    Finite, Two => Finite;
    Finite, MinusOne => Finite;
    Finite, NonZero => Finite;
    Finite, NonNeg => Finite;
    Finite, NonPos => Finite;
    Finite, Finite => Finite;
    Finite, GTEOne => Finite;

    // GTEOne * X
    GTEOne, Pos => Pos;
    GTEOne, Neg => Neg;
    GTEOne, Zero => Zero;
    GTEOne, One => GTEOne;
    GTEOne, Two => GTEOne;
    GTEOne, MinusOne => Neg;
    GTEOne, NonZero => NonZero;
    GTEOne, NonNeg => NonNeg;
    GTEOne, NonPos => NonPos;
    GTEOne, Finite => Finite;
    GTEOne, GTEOne => GTEOne;
}

// ============================================================================
// Division Rules (No division by Zero - that impl doesn't exist!)
// ============================================================================

div_rules! {
    // Pos / X (no Zero!)
    Pos, Pos => Pos;
    Pos, Neg => Neg;
    Pos, One => Pos;
    Pos, Two => Pos;
    Pos, MinusOne => Neg;
    Pos, NonZero => NonZero;
    Pos, GTEOne => Pos;

    // Neg / X
    Neg, Pos => Neg;
    Neg, Neg => Pos;
    Neg, One => Neg;
    Neg, Two => Neg;
    Neg, MinusOne => Pos;
    Neg, NonZero => NonZero;
    Neg, GTEOne => Neg;

    // Zero / X
    Zero, Pos => Zero;
    Zero, Neg => Zero;
    Zero, One => Zero;
    Zero, Two => Zero;
    Zero, MinusOne => Zero;
    Zero, NonZero => Zero;
    Zero, GTEOne => Zero;

    // One / X
    One, Pos => Pos;
    One, Neg => Neg;
    One, One => One;
    One, Two => Pos;
    One, MinusOne => MinusOne;
    One, NonZero => NonZero;
    One, GTEOne => Pos;

    // Two / X
    Two, Pos => Pos;
    Two, Neg => Neg;
    Two, One => Two;
    Two, Two => One;
    Two, MinusOne => Neg;
    Two, NonZero => NonZero;
    Two, GTEOne => Pos;

    // MinusOne / X
    MinusOne, Pos => Neg;
    MinusOne, Neg => Pos;
    MinusOne, One => MinusOne;
    MinusOne, Two => Neg;
    MinusOne, MinusOne => One;
    MinusOne, NonZero => NonZero;
    MinusOne, GTEOne => Neg;

    // NonZero / X
    NonZero, Pos => NonZero;
    NonZero, Neg => NonZero;
    NonZero, One => NonZero;
    NonZero, Two => NonZero;
    NonZero, MinusOne => NonZero;
    NonZero, NonZero => NonZero;
    NonZero, GTEOne => NonZero;

    // NonNeg / X
    NonNeg, Pos => NonNeg;
    NonNeg, Neg => NonPos;
    NonNeg, One => NonNeg;
    NonNeg, Two => NonNeg;
    NonNeg, MinusOne => NonPos;
    NonNeg, NonZero => Finite;
    NonNeg, GTEOne => NonNeg;

    // NonPos / X
    NonPos, Pos => NonPos;
    NonPos, Neg => NonNeg;
    NonPos, One => NonPos;
    NonPos, Two => NonPos;
    NonPos, MinusOne => NonNeg;
    NonPos, NonZero => Finite;
    NonPos, GTEOne => NonPos;

    // Finite / X
    Finite, Pos => Finite;
    Finite, Neg => Finite;
    Finite, One => Finite;
    Finite, Two => Finite;
    Finite, MinusOne => Finite;
    Finite, NonZero => Finite;
    Finite, GTEOne => Finite;

    // GTEOne / X
    GTEOne, Pos => Pos;
    GTEOne, Neg => Neg;
    GTEOne, One => GTEOne;
    GTEOne, Two => Pos;
    GTEOne, MinusOne => Neg;
    GTEOne, NonZero => NonZero;
    GTEOne, GTEOne => Pos;
}

// ============================================================================
// Negation Rules
// ============================================================================

neg_rules! {
    Pos => Neg;
    Neg => Pos;
    Zero => Zero;
    One => MinusOne;
    Two => Neg;
    MinusOne => One;
    NonZero => NonZero;
    NonNeg => NonPos;
    NonPos => NonNeg;
    Finite => Finite;
    GTEOne => Neg;
}

// ============================================================================
// Absolute Value Rules
// ============================================================================

abs_rules! {
    Pos => Pos;
    Neg => Pos;
    Zero => Zero;
    One => One;
    Two => Two;
    MinusOne => One;
    NonZero => Pos;
    NonNeg => NonNeg;
    NonPos => NonNeg;
    Finite => NonNeg;
    GTEOne => GTEOne;
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // Compile-time check that the rules make sense
    fn _check_add<L: AddOutput<R>, R>() {}
    fn _check_sub<L: SubOutput<R>, R>() {}
    fn _check_mul<L: MulOutput<R>, R>() {}
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
