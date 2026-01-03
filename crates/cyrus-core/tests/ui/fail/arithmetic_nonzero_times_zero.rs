//! Test: NonZero * Zero = Zero, not NonZero.
//!
//! This should fail because multiplying by zero gives zero.

use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::{NonZero, Zero};

fn main() {
    let nz = F64::<NonZero>::new(5.0).unwrap();
    let z = F64::<Zero>::ZERO;

    // ERROR: NonZero * Zero = Zero, cannot assign to NonZero
    let _result: F64<NonZero> = nz * z;
}
