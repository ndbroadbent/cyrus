//! Test: f64_neg! rejects positive values at compile time.
//!
//! This should fail to compile because 1.0 is not negative.

use cyrus_core::f64_neg;

fn main() {
    // ERROR: f64_neg! requires negative value
    let _x = f64_neg!(1.0);
}
