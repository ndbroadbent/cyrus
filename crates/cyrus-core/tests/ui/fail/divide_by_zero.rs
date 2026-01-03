//! Test: Division by F64<Zero> does not compile.
//!
//! This should fail because there is no `Div<F64<Zero>>` implementation
//! for any scalar type - division by zero is prevented at the type level.

use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::{Pos, Zero};

fn main() {
    let p = F64::<Pos>::new(5.0).unwrap();
    let z = F64::<Zero>::ZERO;

    // ERROR: no implementation of `Div<F64<Zero>>` for `F64<Pos>`
    let _result = p / z;
}
