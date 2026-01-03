//! Test: Pos + Neg cannot be assigned to Pos (could be zero or negative).
//!
//! This should fail because Pos + Neg = Finite, not Pos.

use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::{Pos, Neg};

fn main() {
    let a = F64::<Pos>::new(5.0).unwrap();
    let b = F64::<Neg>::new(-10.0).unwrap();

    // ERROR: Pos + Neg = Finite, not Pos
    let _result: F64<Pos> = a + b;
}
