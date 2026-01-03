//! Type-level arithmetic rules for tags.
//!
//! These traits define what tag results from combining two tagged values.
//! The rules are mathematical facts independent of the underlying numeric type.
//!
//! # Example
//! ```text
//! Pos + Pos = Pos       (positive + positive = positive)
//! Pos + Neg = Finite    (could be anything)
//! NonZero * NonZero = NonZero
//! Zero + X = X          (additive identity)
//! NonZero.abs() = Pos   (absolute value of non-zero is positive)
//! ```

use super::tags::{Finite, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Zero};

// ============================================================================
// Type-level arithmetic result traits
// ============================================================================

/// Result tag when adding two tagged values.
pub trait AddResult<Rhs> {
    /// The resulting tag type.
    type Output;
}

/// Result tag when subtracting two tagged values.
pub trait SubResult<Rhs> {
    /// The resulting tag type.
    type Output;
}

/// Result tag when multiplying two tagged values.
pub trait MulResult<Rhs> {
    /// The resulting tag type.
    type Output;
}

/// Result tag when dividing by a tagged value.
/// Note: No impl for Zero divisor - division by zero won't compile.
pub trait DivResult<Rhs> {
    /// The resulting tag type.
    type Output;
}

/// Result tag when negating a tagged value.
pub trait NegResult {
    /// The resulting tag type.
    type Output;
}

/// Result tag for absolute value.
pub trait AbsResult {
    /// The resulting tag type.
    type Output;
}

// ============================================================================
// Addition Rules
// ============================================================================

// Pos + Pos = Pos
impl AddResult<Pos> for Pos { type Output = Pos; }

// Neg + Neg = Neg
impl AddResult<Neg> for Neg { type Output = Neg; }

// Pos + Neg = Finite (could cancel)
impl AddResult<Neg> for Pos { type Output = Finite; }
impl AddResult<Pos> for Neg { type Output = Finite; }

// Zero + X = X (identity)
impl AddResult<Pos> for Zero { type Output = Pos; }
impl AddResult<Neg> for Zero { type Output = Neg; }
impl AddResult<Zero> for Zero { type Output = Zero; }
impl AddResult<Finite> for Zero { type Output = Finite; }
impl AddResult<NonZero> for Zero { type Output = NonZero; }
impl AddResult<NonNeg> for Zero { type Output = NonNeg; }
impl AddResult<NonPos> for Zero { type Output = NonPos; }
impl AddResult<One> for Zero { type Output = One; }
impl AddResult<MinusOne> for Zero { type Output = MinusOne; }

// X + Zero = X (identity, symmetric)
impl AddResult<Zero> for Pos { type Output = Pos; }
impl AddResult<Zero> for Neg { type Output = Neg; }
impl AddResult<Zero> for Finite { type Output = Finite; }
impl AddResult<Zero> for NonZero { type Output = NonZero; }
impl AddResult<Zero> for NonNeg { type Output = NonNeg; }
impl AddResult<Zero> for NonPos { type Output = NonPos; }
impl AddResult<Zero> for One { type Output = One; }
impl AddResult<Zero> for MinusOne { type Output = MinusOne; }

// NonNeg + NonNeg = NonNeg
impl AddResult<NonNeg> for NonNeg { type Output = NonNeg; }

// NonPos + NonPos = NonPos
impl AddResult<NonPos> for NonPos { type Output = NonPos; }

// Finite + Finite = Finite
impl AddResult<Finite> for Finite { type Output = Finite; }

// NonZero + anything (except Zero) = Finite (could hit zero)
impl AddResult<Pos> for NonZero { type Output = Finite; }
impl AddResult<Neg> for NonZero { type Output = Finite; }
impl AddResult<NonZero> for NonZero { type Output = Finite; }
impl AddResult<Finite> for NonZero { type Output = Finite; }
impl AddResult<NonZero> for Pos { type Output = Finite; }
impl AddResult<NonZero> for Neg { type Output = Finite; }
impl AddResult<NonZero> for Finite { type Output = Finite; }

// ============================================================================
// Subtraction Rules
// ============================================================================

// Pos - Neg = Pos (subtracting negative makes more positive)
impl SubResult<Neg> for Pos { type Output = Pos; }

// Neg - Pos = Neg (subtracting positive makes more negative)
impl SubResult<Pos> for Neg { type Output = Neg; }

// Pos - Pos = Finite (could be anything)
impl SubResult<Pos> for Pos { type Output = Finite; }

// Neg - Neg = Finite (could be anything)
impl SubResult<Neg> for Neg { type Output = Finite; }

// X - Zero = X
impl SubResult<Zero> for Pos { type Output = Pos; }
impl SubResult<Zero> for Neg { type Output = Neg; }
impl SubResult<Zero> for Zero { type Output = Zero; }
impl SubResult<Zero> for Finite { type Output = Finite; }
impl SubResult<Zero> for NonZero { type Output = NonZero; }
impl SubResult<Zero> for NonNeg { type Output = NonNeg; }
impl SubResult<Zero> for NonPos { type Output = NonPos; }
impl SubResult<Zero> for One { type Output = One; }
impl SubResult<Zero> for MinusOne { type Output = MinusOne; }

// Zero - X = negation of X
impl SubResult<Pos> for Zero { type Output = Neg; }
impl SubResult<Neg> for Zero { type Output = Pos; }
impl SubResult<One> for Zero { type Output = MinusOne; }
impl SubResult<MinusOne> for Zero { type Output = One; }
impl SubResult<Finite> for Zero { type Output = Finite; }
impl SubResult<NonZero> for Zero { type Output = NonZero; }
impl SubResult<NonNeg> for Zero { type Output = NonPos; }
impl SubResult<NonPos> for Zero { type Output = NonNeg; }

// Finite - Finite = Finite
impl SubResult<Finite> for Finite { type Output = Finite; }

// NonZero - X = Finite (could hit zero)
impl SubResult<Pos> for NonZero { type Output = Finite; }
impl SubResult<Neg> for NonZero { type Output = Finite; }
impl SubResult<NonZero> for NonZero { type Output = Finite; }
impl SubResult<Finite> for NonZero { type Output = Finite; }
impl SubResult<NonZero> for Pos { type Output = Finite; }
impl SubResult<NonZero> for Neg { type Output = Finite; }
impl SubResult<NonZero> for Finite { type Output = Finite; }

// ============================================================================
// Multiplication Rules
// ============================================================================

// Sign rules: pos*pos=pos, neg*neg=pos, pos*neg=neg
impl MulResult<Pos> for Pos { type Output = Pos; }
impl MulResult<Neg> for Neg { type Output = Pos; }
impl MulResult<Neg> for Pos { type Output = Neg; }
impl MulResult<Pos> for Neg { type Output = Neg; }

// Zero * anything = Zero
impl MulResult<Pos> for Zero { type Output = Zero; }
impl MulResult<Neg> for Zero { type Output = Zero; }
impl MulResult<Zero> for Zero { type Output = Zero; }
impl MulResult<Finite> for Zero { type Output = Zero; }
impl MulResult<NonZero> for Zero { type Output = Zero; }
impl MulResult<NonNeg> for Zero { type Output = Zero; }
impl MulResult<NonPos> for Zero { type Output = Zero; }
impl MulResult<One> for Zero { type Output = Zero; }
impl MulResult<MinusOne> for Zero { type Output = Zero; }
impl MulResult<Zero> for Pos { type Output = Zero; }
impl MulResult<Zero> for Neg { type Output = Zero; }
impl MulResult<Zero> for Finite { type Output = Zero; }
impl MulResult<Zero> for NonZero { type Output = Zero; }
impl MulResult<Zero> for NonNeg { type Output = Zero; }
impl MulResult<Zero> for NonPos { type Output = Zero; }
impl MulResult<Zero> for One { type Output = Zero; }
impl MulResult<Zero> for MinusOne { type Output = Zero; }

// One * X = X (multiplicative identity)
impl MulResult<Pos> for One { type Output = Pos; }
impl MulResult<Neg> for One { type Output = Neg; }
impl MulResult<Finite> for One { type Output = Finite; }
impl MulResult<NonZero> for One { type Output = NonZero; }
impl MulResult<NonNeg> for One { type Output = NonNeg; }
impl MulResult<NonPos> for One { type Output = NonPos; }
impl MulResult<One> for One { type Output = One; }
impl MulResult<MinusOne> for One { type Output = MinusOne; }
impl MulResult<One> for Pos { type Output = Pos; }
impl MulResult<One> for Neg { type Output = Neg; }
impl MulResult<One> for Finite { type Output = Finite; }
impl MulResult<One> for NonZero { type Output = NonZero; }
impl MulResult<One> for NonNeg { type Output = NonNeg; }
impl MulResult<One> for NonPos { type Output = NonPos; }
impl MulResult<One> for MinusOne { type Output = MinusOne; }

// MinusOne * X = flip sign
impl MulResult<Pos> for MinusOne { type Output = Neg; }
impl MulResult<Neg> for MinusOne { type Output = Pos; }
impl MulResult<MinusOne> for MinusOne { type Output = One; }
impl MulResult<MinusOne> for Pos { type Output = Neg; }
impl MulResult<MinusOne> for Neg { type Output = Pos; }
impl MulResult<Finite> for MinusOne { type Output = Finite; }
impl MulResult<MinusOne> for Finite { type Output = Finite; }
impl MulResult<NonZero> for MinusOne { type Output = NonZero; }
impl MulResult<MinusOne> for NonZero { type Output = NonZero; }

// NonZero * NonZero = NonZero
impl MulResult<NonZero> for NonZero { type Output = NonZero; }
impl MulResult<Pos> for NonZero { type Output = NonZero; }
impl MulResult<Neg> for NonZero { type Output = NonZero; }
impl MulResult<NonZero> for Pos { type Output = NonZero; }
impl MulResult<NonZero> for Neg { type Output = NonZero; }

// NonNeg * NonNeg = NonNeg
impl MulResult<NonNeg> for NonNeg { type Output = NonNeg; }

// NonPos * NonPos = NonNeg
impl MulResult<NonPos> for NonPos { type Output = NonNeg; }

// NonNeg * NonPos = NonPos
impl MulResult<NonPos> for NonNeg { type Output = NonPos; }
impl MulResult<NonNeg> for NonPos { type Output = NonPos; }

// Finite * Finite = Finite
impl MulResult<Finite> for Finite { type Output = Finite; }

// ============================================================================
// Division Rules (no Zero divisor - won't compile)
// ============================================================================

// Pos / Pos = Pos, Neg / Neg = Pos
impl DivResult<Pos> for Pos { type Output = Pos; }
impl DivResult<Neg> for Neg { type Output = Pos; }

// Pos / Neg = Neg, Neg / Pos = Neg
impl DivResult<Neg> for Pos { type Output = Neg; }
impl DivResult<Pos> for Neg { type Output = Neg; }

// Zero / NonZero = Zero
impl DivResult<Pos> for Zero { type Output = Zero; }
impl DivResult<Neg> for Zero { type Output = Zero; }
impl DivResult<NonZero> for Zero { type Output = Zero; }
impl DivResult<One> for Zero { type Output = Zero; }
impl DivResult<MinusOne> for Zero { type Output = Zero; }

// X / One = X
impl DivResult<One> for Pos { type Output = Pos; }
impl DivResult<One> for Neg { type Output = Neg; }
impl DivResult<One> for Finite { type Output = Finite; }
impl DivResult<One> for NonZero { type Output = NonZero; }
impl DivResult<One> for NonNeg { type Output = NonNeg; }
impl DivResult<One> for NonPos { type Output = NonPos; }
impl DivResult<One> for One { type Output = One; }
impl DivResult<One> for MinusOne { type Output = MinusOne; }

// X / MinusOne = flip sign
impl DivResult<MinusOne> for Pos { type Output = Neg; }
impl DivResult<MinusOne> for Neg { type Output = Pos; }
impl DivResult<MinusOne> for One { type Output = MinusOne; }
impl DivResult<MinusOne> for MinusOne { type Output = One; }
impl DivResult<MinusOne> for Finite { type Output = Finite; }
impl DivResult<MinusOne> for NonZero { type Output = NonZero; }

// NonZero / NonZero = NonZero
impl DivResult<NonZero> for NonZero { type Output = NonZero; }
impl DivResult<Pos> for NonZero { type Output = NonZero; }
impl DivResult<Neg> for NonZero { type Output = NonZero; }
impl DivResult<NonZero> for Pos { type Output = NonZero; }
impl DivResult<NonZero> for Neg { type Output = NonZero; }

// Finite / NonZero = Finite
impl DivResult<Pos> for Finite { type Output = Finite; }
impl DivResult<Neg> for Finite { type Output = Finite; }
impl DivResult<NonZero> for Finite { type Output = Finite; }

// NonNeg / Pos = NonNeg, NonNeg / Neg = NonPos
impl DivResult<Pos> for NonNeg { type Output = NonNeg; }
impl DivResult<Neg> for NonNeg { type Output = NonPos; }

// NonPos / Pos = NonPos, NonPos / Neg = NonNeg
impl DivResult<Pos> for NonPos { type Output = NonPos; }
impl DivResult<Neg> for NonPos { type Output = NonNeg; }

// ============================================================================
// Negation Rules
// ============================================================================

impl NegResult for Pos { type Output = Neg; }
impl NegResult for Neg { type Output = Pos; }
impl NegResult for Zero { type Output = Zero; }
impl NegResult for One { type Output = MinusOne; }
impl NegResult for MinusOne { type Output = One; }
impl NegResult for Finite { type Output = Finite; }
impl NegResult for NonZero { type Output = NonZero; }
impl NegResult for NonNeg { type Output = NonPos; }
impl NegResult for NonPos { type Output = NonNeg; }

// ============================================================================
// Absolute Value Rules
// ============================================================================

impl AbsResult for Pos { type Output = Pos; }
impl AbsResult for Neg { type Output = Pos; }
impl AbsResult for Zero { type Output = Zero; }
impl AbsResult for One { type Output = One; }
impl AbsResult for MinusOne { type Output = One; }
impl AbsResult for Finite { type Output = NonNeg; }
impl AbsResult for NonZero { type Output = Pos; }  // KEY: NonZero.abs() = Pos!
impl AbsResult for NonNeg { type Output = NonNeg; }
impl AbsResult for NonPos { type Output = NonNeg; }

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::any::TypeId;

    fn assert_same<A: 'static, B: 'static>() {
        assert_eq!(TypeId::of::<A>(), TypeId::of::<B>());
    }

    #[test]
    fn test_add_rules() {
        assert_same::<<Pos as AddResult<Pos>>::Output, Pos>();
        assert_same::<<Neg as AddResult<Neg>>::Output, Neg>();
        assert_same::<<Pos as AddResult<Neg>>::Output, Finite>();
        assert_same::<<Zero as AddResult<Pos>>::Output, Pos>();
    }

    #[test]
    fn test_mul_rules() {
        assert_same::<<Pos as MulResult<Pos>>::Output, Pos>();
        assert_same::<<Neg as MulResult<Neg>>::Output, Pos>();
        assert_same::<<Pos as MulResult<Neg>>::Output, Neg>();
        assert_same::<<Zero as MulResult<Pos>>::Output, Zero>();
        assert_same::<<NonZero as MulResult<NonZero>>::Output, NonZero>();
    }

    #[test]
    fn test_div_rules() {
        assert_same::<<Pos as DivResult<Pos>>::Output, Pos>();
        assert_same::<<Neg as DivResult<Neg>>::Output, Pos>();
        assert_same::<<NonZero as DivResult<NonZero>>::Output, NonZero>();
    }

    #[test]
    fn test_abs_rules() {
        assert_same::<<Pos as AbsResult>::Output, Pos>();
        assert_same::<<Neg as AbsResult>::Output, Pos>();
        assert_same::<<NonZero as AbsResult>::Output, Pos>();  // The key insight!
        assert_same::<<Finite as AbsResult>::Output, NonNeg>();
    }

    #[test]
    fn test_neg_rules() {
        assert_same::<<Pos as NegResult>::Output, Neg>();
        assert_same::<<Neg as NegResult>::Output, Pos>();
        assert_same::<<One as NegResult>::Output, MinusOne>();
    }
}
