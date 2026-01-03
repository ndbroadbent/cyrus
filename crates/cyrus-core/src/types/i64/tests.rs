#![allow(clippy::float_cmp)]
use super::*;
use std::collections::hash_map::DefaultHasher;

#[test]
fn test_i64_size() {
    assert_eq!(std::mem::size_of::<I64>(), std::mem::size_of::<i64>());
    assert_eq!(std::mem::size_of::<I64<Pos>>(), std::mem::size_of::<i64>());
}

#[test]
fn test_pos_new() {
    assert!(I64::<Pos>::new(1).is_some());
    assert!(I64::<Pos>::new(0).is_none());
    assert!(I64::<Pos>::new(-1).is_none());
}

#[test]
fn test_neg_new() {
    assert!(I64::<Neg>::new(-1).is_some());
    assert!(I64::<Neg>::new(0).is_none());
    assert!(I64::<Neg>::new(1).is_none());
}

#[test]
fn test_nonzero_new() {
    assert!(I64::<NonZero>::new(1).is_some());
    assert!(I64::<NonZero>::new(-1).is_some());
    assert!(I64::<NonZero>::new(0).is_none());
}

#[test]
fn test_nonneg_new() {
    assert!(I64::<NonNeg>::new(0).is_some());
    assert!(I64::<NonNeg>::new(1).is_some());
    assert!(I64::<NonNeg>::new(-1).is_none());
}

#[test]
fn test_nonpos_new() {
    assert!(I64::<NonPos>::new(0).is_some());
    assert!(I64::<NonPos>::new(-1).is_some());
    assert!(I64::<NonPos>::new(1).is_none());
}

#[test]
fn test_zero_new() {
    assert!(I64::<Zero>::new(0).is_some());
    assert!(I64::<Zero>::new(1).is_none());
    assert!(I64::<Zero>::new(-1).is_none());
}

#[test]
fn test_one_new() {
    assert!(I64::<One>::new(1).is_some());
    assert!(I64::<One>::new(0).is_none());
    assert!(I64::<One>::new(2).is_none());
}

#[test]
fn test_minusone_new() {
    assert!(I64::<MinusOne>::new(-1).is_some());
    assert!(I64::<MinusOne>::new(0).is_none());
    assert!(I64::<MinusOne>::new(-2).is_none());
}

#[test]
fn test_constants() {
    assert_eq!(I64::<Zero>::ZERO.get(), 0);
    assert_eq!(I64::<One>::ONE.get(), 1);
    assert_eq!(I64::<MinusOne>::MINUS_ONE.get(), -1);
    assert_eq!(I64::<Finite>::ZERO.get(), 0);
}

#[test]
fn test_finite_new() {
    let x = I64::<Finite>::new(42);
    assert_eq!(x.get(), 42);
    let y = I64::<Finite>::new(-100);
    assert_eq!(y.get(), -100);
}

#[test]
#[allow(clippy::float_cmp)] // Integer-to-f64 conversions are exact for small values
fn test_to_f64_conversions() {
    let finite = I64::<Finite>::new(42);
    assert_eq!(finite.to_f64().get(), 42.0);

    let pos = I64::<Pos>::new(10).unwrap();
    assert_eq!(pos.to_f64().get(), 10.0);

    let neg = I64::<Neg>::new(-5).unwrap();
    assert_eq!(neg.to_f64().get(), -5.0);

    let nonzero = I64::<NonZero>::new(7).unwrap();
    assert_eq!(nonzero.to_f64().get(), 7.0);

    let nonneg = I64::<NonNeg>::new(3).unwrap();
    assert_eq!(nonneg.to_f64().get(), 3.0);

    let nonpos = I64::<NonPos>::new(-8).unwrap();
    assert_eq!(nonpos.to_f64().get(), -8.0);

    let zero = I64::<Zero>::new(0).unwrap();
    assert_eq!(zero.to_f64().get(), 0.0);

    let one = I64::<One>::new(1).unwrap();
    assert_eq!(one.to_f64().get(), 1.0);

    let minusone = I64::<MinusOne>::new(-1).unwrap();
    assert_eq!(minusone.to_f64().get(), -1.0);
}

#[test]
fn test_debug_display() {
    let x = I64::<Pos>::new(42).unwrap();
    assert!(format!("{x:?}").contains("Pos"));
    assert!(format!("{x:?}").contains("42"));
    assert_eq!(format!("{x}"), "42");

    let y = I64::<Finite>::new(-5);
    assert!(format!("{y:?}").contains("Finite"));
}

#[test]
fn test_eq_hash() {
    let a = I64::<Pos>::new(5).unwrap();
    let b = I64::<Pos>::new(5).unwrap();
    let c = I64::<Pos>::new(6).unwrap();
    assert_eq!(a, b);
    assert_ne!(a, c);

    let mut h1 = DefaultHasher::new();
    let mut h2 = DefaultHasher::new();
    a.hash(&mut h1);
    b.hash(&mut h2);
    assert_eq!(h1.finish(), h2.finish());
}

#[test]
fn test_default_zero() {
    let z: I64<Zero> = I64::default();
    assert_eq!(z.get(), 0);
}

#[test]
fn test_macros() {
    let x = i64_pos!(42);
    assert_eq!(x.get(), 42);

    let y = i64_neg!(-5);
    assert_eq!(y.get(), -5);

    let z = i64_nonneg!(0);
    assert_eq!(z.get(), 0);
    let z2 = i64_nonneg!(10);
    assert_eq!(z2.get(), 10);

    let w = i64_nonpos!(0);
    assert_eq!(w.get(), 0);
    let w2 = i64_nonpos!(-10);
    assert_eq!(w2.get(), -10);
}

#[test]
fn test_ord() {
    let a = I64::<Pos>::new(1).unwrap();
    let b = I64::<Pos>::new(2).unwrap();
    assert!(a < b);
    assert!(b > a);
    assert!(a <= b);
    assert!(b >= a);
    assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Less));
}

#[test]
fn test_from_raw() {
    let x = I64::<Pos>::from_raw(100);
    assert_eq!(x.get(), 100);
}
