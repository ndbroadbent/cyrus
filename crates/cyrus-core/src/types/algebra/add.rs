//! Addition rules for type-level arithmetic.

use super::super::tags::{
    Finite, GTEOne, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Two, Zero,
};
use super::{AddOutput, add_rules};

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
