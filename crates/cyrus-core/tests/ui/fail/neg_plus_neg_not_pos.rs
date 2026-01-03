//! Test: Neg + Neg = Neg, not Pos.
//!
//! Adding two negatives gives a more negative number.

use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::{Neg, Pos};

fn main() {
    let a = F64::<Neg>::new(-2.0).unwrap();
    let b = F64::<Neg>::new(-3.0).unwrap();

    // ERROR: Neg + Neg = Neg, not Pos
    let _result: F64<Pos> = a + b;
}
