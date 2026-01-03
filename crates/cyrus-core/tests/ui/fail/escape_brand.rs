//! Test: Brands cannot escape their closure.
//!
//! This should fail to compile because the HRTB constraint
//! `for<'id> FnOnce(Brand<'id>) -> R` prevents the brand from escaping.

use cyrus_core::types::dimension::{Dim, Brand};

fn main() {
    let dim = Dim::new(3);

    // ERROR: brand cannot escape the closure due to HRTB
    let escaped: Brand<'_> = dim.with_brand(|brand| brand);

    // If this compiled, we could use the escaped brand unsafely
    let _ = escaped;
}
