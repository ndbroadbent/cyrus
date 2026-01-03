//! Test: BVec created with wrong brand cannot access with brand's indices.
//!
//! This demonstrates that even if dimensions match, the brand system
//! prevents mixing vectors from different contexts.

use cyrus_core::types::dimension::{Dim, BVec};

fn main() {
    let dim = Dim::new(3);

    // Create a vector in one brand context
    let result = dim.with_brand(|brand1| {
        let v = BVec::new(&dim, &brand1, vec![1, 2, 3]);

        // Try to use indices from a different brand context
        dim.with_brand(|brand2| {
            // These indices are from brand2
            for i in dim.indices(&brand2) {
                // ERROR: v is BVec<'brand1, _> but i is Index<'brand2>
                let _ = v.get(i);
            }
        });

        v.into_inner()
    });

    let _ = result;
}
