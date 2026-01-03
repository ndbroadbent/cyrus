//! Multiplication rules for type-level arithmetic.

use super::super::tags::{
    Finite, GTEOne, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Two, Zero,
};
use super::{MulOutput, mul_rules};

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
