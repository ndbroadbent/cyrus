#![allow(clippy::float_cmp, clippy::approx_constant)]
use super::*;
use std::collections::hash_map::DefaultHasher;

#[test]
fn test_f64_size() {
    assert_eq!(std::mem::size_of::<F64>(), std::mem::size_of::<f64>());
    assert_eq!(std::mem::size_of::<F64<Pos>>(), std::mem::size_of::<f64>());
}

#[test]
fn test_finite_new() {
    assert!(F64::<Finite>::new(0.0).is_some());
    assert!(F64::<Finite>::new(1.0).is_some());
    assert!(F64::<Finite>::new(-1.0).is_some());
    assert!(F64::<Finite>::new(f64::NAN).is_none());
    assert!(F64::<Finite>::new(f64::INFINITY).is_none());
    assert!(F64::<Finite>::new(f64::NEG_INFINITY).is_none());
}

#[test]
fn test_pos_new() {
    assert!(F64::<Pos>::new(1.0).is_some());
    assert!(F64::<Pos>::new(0.0).is_none());
    assert!(F64::<Pos>::new(-1.0).is_none());
    assert!(F64::<Pos>::new(f64::NAN).is_none());
    assert!(F64::<Pos>::new(f64::INFINITY).is_none());
}

#[test]
fn test_neg_new() {
    assert!(F64::<Neg>::new(-1.0).is_some());
    assert!(F64::<Neg>::new(0.0).is_none());
    assert!(F64::<Neg>::new(1.0).is_none());
    assert!(F64::<Neg>::new(f64::NAN).is_none());
    assert!(F64::<Neg>::new(f64::NEG_INFINITY).is_none());
}

#[test]
fn test_nonzero_new() {
    assert!(F64::<NonZero>::new(1.0).is_some());
    assert!(F64::<NonZero>::new(-1.0).is_some());
    assert!(F64::<NonZero>::new(0.0).is_none());
    assert!(F64::<NonZero>::new(f64::NAN).is_none());
}

#[test]
fn test_nonneg_new() {
    assert!(F64::<NonNeg>::new(0.0).is_some());
    assert!(F64::<NonNeg>::new(1.0).is_some());
    assert!(F64::<NonNeg>::new(-1.0).is_none());
    assert!(F64::<NonNeg>::new(f64::NAN).is_none());
}

#[test]
fn test_nonpos_new() {
    assert!(F64::<NonPos>::new(0.0).is_some());
    assert!(F64::<NonPos>::new(-1.0).is_some());
    assert!(F64::<NonPos>::new(1.0).is_none());
    assert!(F64::<NonPos>::new(f64::NAN).is_none());
}

#[test]
fn test_zero_new() {
    assert!(F64::<Zero>::new(0.0).is_some());
    assert!(F64::<Zero>::new(1.0).is_none());
    assert!(F64::<Zero>::new(-0.0).is_some()); // -0.0 == 0.0 in f64
}

#[test]
fn test_one_new() {
    assert!(F64::<One>::new(1.0).is_some());
    assert!(F64::<One>::new(0.0).is_none());
    assert!(F64::<One>::new(2.0).is_none());
}

#[test]
fn test_minusone_new() {
    assert!(F64::<MinusOne>::new(-1.0).is_some());
    assert!(F64::<MinusOne>::new(0.0).is_none());
    assert!(F64::<MinusOne>::new(-2.0).is_none());
}

#[test]
fn test_two_new() {
    assert!(F64::<Two>::new(2.0).is_some());
    assert!(F64::<Two>::new(1.0).is_none());
    assert!(F64::<Two>::new(0.0).is_none());
}

#[test]
fn test_gteone_new() {
    assert!(F64::<GTEOne>::new(1.0).is_some());
    assert!(F64::<GTEOne>::new(2.0).is_some());
    assert!(F64::<GTEOne>::new(100.0).is_some());
    assert!(F64::<GTEOne>::new(0.5).is_none());
    assert!(F64::<GTEOne>::new(0.0).is_none());
    assert!(F64::<GTEOne>::new(f64::NAN).is_none());
    assert!(F64::<GTEOne>::new(f64::INFINITY).is_none());
}

#[test]
#[allow(clippy::float_cmp)] // Exact constant comparisons
fn test_constants() {
    assert_eq!(F64::<Zero>::ZERO.get(), 0.0);
    assert_eq!(F64::<One>::ONE.get(), 1.0);
    assert_eq!(F64::<MinusOne>::MINUS_ONE.get(), -1.0);
    assert_eq!(F64::<Two>::TWO.get(), 2.0);
    assert_eq!(F64::<Finite>::ZERO.get(), 0.0);
    assert_eq!(F64::<NonNeg>::ZERO.get(), 0.0);
}

#[test]
#[allow(clippy::float_cmp)] // Exact constant comparisons
fn test_abs() {
    let f = F64::<Finite>::new(-5.0).unwrap();
    assert_eq!(f.abs().get(), 5.0);

    let p = F64::<Pos>::new(3.0).unwrap();
    assert_eq!(p.abs().get(), 3.0);

    let n = F64::<Neg>::new(-7.0).unwrap();
    assert_eq!(n.abs().get(), 7.0);

    let nz = F64::<NonZero>::new(-4.0).unwrap();
    assert_eq!(nz.abs().get(), 4.0);
}

#[test]
fn test_square() {
    let f = F64::<Finite>::new(-3.0).unwrap();
    assert!((f.square().get() - 9.0).abs() < 1e-10);
}

#[test]
fn test_sqrt() {
    let p = F64::<Pos>::new(4.0).unwrap();
    assert!((p.sqrt().get() - 2.0).abs() < 1e-10);

    let nn = F64::<NonNeg>::new(9.0).unwrap();
    assert!((nn.sqrt().get() - 3.0).abs() < 1e-10);
}

#[test]
fn test_ln() {
    let p = F64::<Pos>::new(std::f64::consts::E).unwrap();
    assert!((p.ln().get() - 1.0).abs() < 1e-10);
}

#[test]
fn test_recip() {
    let p = F64::<Pos>::new(4.0).unwrap();
    assert!((p.recip().get() - 0.25).abs() < 1e-10);
}

#[test]
fn test_try_narrowing() {
    let f = F64::<Finite>::new(5.0).unwrap();
    assert!(f.try_to_pos().is_some());
    assert!(f.try_to_neg().is_none());
    assert!(f.try_to_non_zero().is_some());

    let n = F64::<Finite>::new(-3.0).unwrap();
    assert!(n.try_to_pos().is_none());
    assert!(n.try_to_neg().is_some());
    assert!(n.try_to_non_zero().is_some());

    let z = F64::<Finite>::new(0.0).unwrap();
    assert!(z.try_to_pos().is_none());
    assert!(z.try_to_neg().is_none());
    assert!(z.try_to_non_zero().is_none());
}

#[test]
fn test_display() {
    let x = F64::<Pos>::new(3.14).unwrap();
    assert_eq!(format!("{x}"), "3.14");
}

#[test]
fn test_debug() {
    let x = F64::<Pos>::new(3.14).unwrap();
    let s = format!("{x:?}");
    assert!(s.contains("Pos"));
    assert!(s.contains("3.14"));

    let y = F64::<Finite>::new(-5.0).unwrap();
    assert!(format!("{y:?}").contains("Finite"));
}

#[test]
fn test_eq() {
    let a = F64::<Pos>::new(3.0).unwrap();
    let b = F64::<Pos>::new(3.0).unwrap();
    let c = F64::<Pos>::new(4.0).unwrap();
    assert_eq!(a, b);
    assert_ne!(a, c);
}

#[test]
fn test_hash() {
    use std::collections::HashSet;
    let mut set = HashSet::new();
    set.insert(F64::<Pos>::new(1.0).unwrap());
    set.insert(F64::<Pos>::new(2.0).unwrap());
    assert_eq!(set.len(), 2);
}

#[test]
fn test_hash_negative_zero() {
    // -0.0 and 0.0 should hash the same
    let zero = F64::<Finite>::new(0.0).unwrap();
    let neg_zero = F64::<Finite>::new(-0.0).unwrap();

    let mut h1 = DefaultHasher::new();
    let mut h2 = DefaultHasher::new();
    zero.hash(&mut h1);
    neg_zero.hash(&mut h2);
    assert_eq!(h1.finish(), h2.finish());
}

#[test]
fn test_partial_ord() {
    let a = F64::<Pos>::new(1.0).unwrap();
    let b = F64::<Pos>::new(2.0).unwrap();
    assert!(a < b);
    assert!(b > a);
    assert!(a <= b);
    assert!(b >= a);
    assert_eq!(a.partial_cmp(&b), Some(std::cmp::Ordering::Less));
}

#[test]
#[allow(clippy::float_cmp)] // Exact zero comparison
fn test_default_zero() {
    let z: F64<Zero> = F64::default();
    assert_eq!(z.get(), 0.0);
}

#[test]
#[allow(clippy::float_cmp)] // Exact constant comparison
fn test_from_raw() {
    let x = F64::<Pos>::from_raw(100.0);
    assert_eq!(x.get(), 100.0);
}

#[test]
fn test_widening_via_algebra() {
    // Widening is automatic through the type algebra.
    // Pos + Zero = Pos (identity for Pos)
    // One + Zero = One (identity for One)
    // The algebra handles cross-type operations automatically.
    let p = F64::<Pos>::new(1.0).unwrap();
    let z = F64::<Zero>::ZERO;

    // Pos + Zero = Pos (algebra preserves strongest type where possible)
    let _: F64<Pos> = p + z;

    // Pos * One = Pos
    let o = F64::<One>::ONE;
    let _: F64<Pos> = p * o;

    // When we need Finite, the algebra produces it from mixed-sign ops
    let n = F64::<Neg>::new(-1.0).unwrap();
    let _: F64<Finite> = p + n; // Pos + Neg = Finite
}
