//! Test: Finite.abs() is NonNeg, not Pos.
//!
//! A finite value could be zero, so abs() could be zero.

use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::{Finite, Pos};

fn main() {
    let f = F64::<Finite>::new(0.0).unwrap();

    // ERROR: Finite.abs() = NonNeg, not Pos (could be zero!)
    let _result: F64<Pos> = f.abs();
}
