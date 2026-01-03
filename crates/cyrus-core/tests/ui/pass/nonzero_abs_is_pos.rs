//! Test: NonZero.abs() returns Pos, not just NonNeg.
//!
//! This is the key insight: |x| where x ≠ 0 is strictly positive.

use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::{NonZero, Pos};

fn main() {
    let nz = F64::<NonZero>::new(-7.0).unwrap();

    // The key rule: NonZero.abs() = Pos (not NonNeg!)
    let abs_val: F64<Pos> = nz.abs();

    assert!(abs_val.get() > 0.0);
}
