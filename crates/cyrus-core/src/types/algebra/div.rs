//! Division rules for type-level arithmetic.
//!
//! Note: Division by Zero is intentionally not implemented - it won't compile!

use super::super::tags::{
    Finite, GTEOne, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Two, Zero,
};
use super::{DivOutput, div_rules};

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

    // Zero / X (Zero divided by non-zero = Zero)
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

    // Finite / X (must be NonZero to divide!)
    Finite, Pos => Finite;
    Finite, Neg => Finite;
    Finite, One => Finite;
    Finite, Two => Finite;
    Finite, MinusOne => Finite;
    Finite, NonZero => Finite;
    Finite, GTEOne => Finite;
    // Note: Finite / Finite is intentionally NOT implemented because
    // Finite might be zero. Division by Finite requires explicit NonZero conversion.

    // GTEOne / X
    GTEOne, Pos => Pos;
    GTEOne, Neg => Neg;
    GTEOne, One => GTEOne;
    GTEOne, Two => Pos;
    GTEOne, MinusOne => Neg;
    GTEOne, NonZero => NonZero;
    GTEOne, GTEOne => Pos;
}
