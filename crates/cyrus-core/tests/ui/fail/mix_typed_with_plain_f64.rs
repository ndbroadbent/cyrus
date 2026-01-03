//! Test: Mixing F64<Tag> with plain f64 does not compile.
//!
//! This should fail because there are no operator implementations
//! between typed F64 values and plain f64 - all values must be typed.

use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::Pos;

fn main() {
    let typed = F64::<Pos>::new(5.0).unwrap();
    let plain: f64 = 3.0;

    // ERROR: no implementation of `Add<f64>` for `F64<Pos>`
    let _result = typed + plain;
}
