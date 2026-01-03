//! Test: Zero is the additive identity.

use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::{Zero, Pos, Neg};

fn main() {
    let z = F64::<Zero>::ZERO;
    let p = F64::<Pos>::new(5.0).unwrap();
    let n = F64::<Neg>::new(-3.0).unwrap();

    // Zero + Pos = Pos
    let _a: F64<Pos> = z + p;

    // Pos + Zero = Pos
    let _b: F64<Pos> = p + z;

    // Zero + Neg = Neg
    let _c: F64<Neg> = z + n;

    // Zero - Pos = Neg (0 - positive = negative)
    let _d: F64<Neg> = z - p;

    // Zero - Neg = Pos (0 - negative = positive)
    let _e: F64<Pos> = z - n;
}
