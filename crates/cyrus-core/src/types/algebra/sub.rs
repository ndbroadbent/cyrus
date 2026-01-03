//! Subtraction rules for type-level arithmetic.

use super::super::tags::{
    Finite, GTEOne, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Two, Zero,
};
use super::{SubOutput, sub_rules};

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
