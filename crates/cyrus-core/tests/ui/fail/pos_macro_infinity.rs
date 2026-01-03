//! Test: f64_pos! rejects infinity at compile time.
//!
//! This should fail to compile because infinity is not finite.

use cyrus_core::f64_pos;

fn main() {
    // ERROR: f64_pos! requires finite positive value
    let _x = f64_pos!(f64::INFINITY);
}
