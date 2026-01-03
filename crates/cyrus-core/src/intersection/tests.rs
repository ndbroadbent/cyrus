use super::*;
use malachite::Rational;

#[test]
fn test_set_get() {
    let mut kappa = Intersection::new(3);
    let val = TypedRational::<Finite>::from_raw(Rational::from(5));
    kappa.set(0, 1, 2, val);

    assert_eq!(*kappa.get(0, 1, 2).get(), Rational::from(5));
    assert_eq!(*kappa.get(2, 1, 0).get(), Rational::from(5));
}
