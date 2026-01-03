//! Unary operation rules (negation, absolute value).

use super::super::tags::{
    Finite, GTEOne, MinusOne, Neg, NonNeg, NonPos, NonZero, One, Pos, Two, Zero,
};
use super::{AbsOutput, NegOutput, abs_rules, neg_rules};

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
