//! Test: Negation flips the sign tag.

use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::{Pos, Neg, One, MinusOne};

fn main() {
    let p = F64::<Pos>::new(5.0).unwrap();
    let n = F64::<Neg>::new(-3.0).unwrap();

    // -Pos = Neg
    let _a: F64<Neg> = -p;

    // -Neg = Pos
    let _b: F64<Pos> = -n;

    // -One = MinusOne
    let _c: F64<MinusOne> = -F64::<One>::ONE;

    // -MinusOne = One
    let _d: F64<One> = -F64::<MinusOne>::MINUS_ONE;
}
