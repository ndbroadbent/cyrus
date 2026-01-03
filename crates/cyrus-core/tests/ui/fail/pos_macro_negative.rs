//! Test: f64_pos! rejects negative values at compile time.
//!
//! This should fail to compile because -1.0 is not positive.

use cyrus_core::f64_pos;

fn main() {
    // ERROR: f64_pos! requires positive value
    let _x = f64_pos!(-1.0);
}
