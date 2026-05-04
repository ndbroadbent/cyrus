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

#[test]
fn test_iter_is_key_sorted() {
    let mut kappa = Intersection::new(4);
    for (i, j, k, value) in [(2, 3, 1, 7), (0, 0, 0, 1), (1, 1, 0, 3), (3, 3, 3, 9)] {
        kappa.set(
            i,
            j,
            k,
            TypedRational::<Finite>::from_raw(Rational::from(value)),
        );
    }

    let keys: Vec<_> = kappa.iter().map(|(key, _)| *key).collect();
    assert_eq!(keys, vec![(0, 0, 0), (0, 1, 1), (1, 2, 3), (3, 3, 3)]);
}
