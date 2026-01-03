//! Test: Indices cannot be forged from raw usize.
//!
//! This should fail to compile because `Index` has private fields
//! and can only be created through `Dim::indices` or `Dim::index`.

use cyrus_core::types::dimension::{Dim, BVec, Index};
use std::marker::PhantomData;

fn main() {
    let dim = Dim::new(3);

    dim.with_brand(|brand| {
        let v = BVec::new(&dim, &brand, vec![1, 2, 3]);

        // ERROR: can't create Index from raw usize - field is private
        let fake_idx = Index {
            i: 100, // Out of bounds!
            _brand: PhantomData,
        };

        // If this compiled, it would be UB
        let _ = v.get(fake_idx);
    });
}
