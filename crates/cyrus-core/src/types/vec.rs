//! Non-empty vector type.
//!
//! `NonEmptyVec<T>` guarantees at least one element, eliminating
//! the need for scattered empty checks throughout the codebase.

/// A vector guaranteed to be non-empty.
#[derive(Debug, Clone, PartialEq)]
pub struct NonEmptyVec<T>(Vec<T>);

impl<T> NonEmptyVec<T> {
    /// Create a new non-empty vector. Returns None if vec is empty.
    pub fn new(vec: Vec<T>) -> Option<Self> {
        if vec.is_empty() {
            None
        } else {
            Some(Self(vec))
        }
    }

    /// Get a reference to the inner slice.
    pub fn as_slice(&self) -> &[T] {
        &self.0
    }

    /// Get the first element. Always succeeds since vec is non-empty.
    pub fn first(&self) -> &T {
        &self.0[0]
    }

    /// Get the length. Always >= 1.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Iterate over elements.
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.0.iter()
    }

    /// Consume and return the inner vec.
    pub fn into_inner(self) -> Vec<T> {
        self.0
    }
}

impl<T> std::ops::Index<usize> for NonEmptyVec<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_non_empty_vec_rejects_empty() {
        let empty: Vec<i32> = vec![];
        assert!(NonEmptyVec::new(empty).is_none());
    }

    #[test]
    fn test_non_empty_vec_accepts_non_empty() {
        let v = NonEmptyVec::new(vec![1, 2, 3]).unwrap();
        assert_eq!(v.len(), 3);
        assert_eq!(v.first(), &1);
    }

    #[test]
    fn test_non_empty_vec_index() {
        let v = NonEmptyVec::new(vec![10, 20, 30]).unwrap();
        assert_eq!(v[0], 10);
        assert_eq!(v[2], 30);
    }

    #[test]
    fn test_non_empty_vec_into_inner() {
        let v = NonEmptyVec::new(vec![1, 2, 3]).unwrap();
        let inner = v.into_inner();
        assert_eq!(inner, vec![1, 2, 3]);
    }

    #[test]
    fn test_non_empty_vec_as_slice() {
        let v = NonEmptyVec::new(vec![1, 2, 3]).unwrap();
        assert_eq!(v.as_slice(), &[1, 2, 3]);
    }

    #[test]
    fn test_non_empty_vec_iter() {
        let v = NonEmptyVec::new(vec![1, 2, 3]).unwrap();
        let sum: i32 = v.iter().sum();
        assert_eq!(sum, 6);
    }
}
