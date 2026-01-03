//! Test: Basic arithmetic type rules compile correctly.

use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::{Pos, Neg, Finite};

fn main() {
    // Pos + Pos = Pos
    let a = F64::<Pos>::new(2.0).unwrap();
    let b = F64::<Pos>::new(3.0).unwrap();
    let _c: F64<Pos> = a + b;

    // Neg * Neg = Pos
    let x = F64::<Neg>::new(-2.0).unwrap();
    let y = F64::<Neg>::new(-3.0).unwrap();
    let _z: F64<Pos> = x * y;

    // Pos / Pos = Pos
    let p = F64::<Pos>::new(6.0).unwrap();
    let q = F64::<Pos>::new(2.0).unwrap();
    let _r: F64<Pos> = p / q;

    // Pos + Neg = Finite
    let _f: F64<Finite> = a + x;

    // Pos * Neg = Neg
    let _n: F64<Neg> = a * x;
}
