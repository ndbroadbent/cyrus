//! Test: Different brands cannot be mixed.
//!
//! This should fail to compile because `BVec<'id1, _>` and `Index<'id2>`
//! have different brand lifetimes.

use cyrus_core::types::dimension::{Dim, BVec};

fn main() {
    let d1 = Dim::new(3);
    let d2 = Dim::new(3);

    d1.with_brand(|b1| {
        d2.with_brand(|b2| {
            let _a = BVec::new(&d1, &b1, vec![1, 2, 3]);
            let _b = BVec::new(&d2, &b2, vec![4, 5, 6]);

            // ERROR: trying to use index from d1's brand on d2's vector
            for i in d1.indices(&b1) {
                // This should not compile - i is Index<'id1> but _b is BVec<'id2, _>
                let _ = _b.get(i);
            }
        })
    });
}
