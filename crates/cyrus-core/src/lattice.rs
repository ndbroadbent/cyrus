//! Lattice point operations for toric geometry.
//!
//! Provides efficient lattice point representation and enumeration
//! needed for polytope analysis.

use std::ops::{Add, Neg, Sub};

/// A point in an integer lattice Z^n.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Point {
    coords: Vec<i64>,
}

impl Point {
    /// Create a new lattice point from coordinates.
    #[inline]
    pub const fn new(coords: Vec<i64>) -> Self {
        Self { coords }
    }

    /// Create the origin in dimension n.
    #[inline]
    pub fn origin(dim: usize) -> Self {
        Self::new(vec![0; dim])
    }

    /// Get the dimension of this point.
    #[inline]
    pub const fn dim(&self) -> usize {
        self.coords.len()
    }

    /// Get the coordinates as a slice.
    #[inline]
    pub fn coords(&self) -> &[i64] {
        &self.coords
    }

    /// Check if this is the origin.
    #[inline]
    pub fn is_origin(&self) -> bool {
        self.coords.iter().all(|&x| x == 0)
    }

    /// Compute the dot product with another point.
    ///
    /// # Panics
    /// Panics if the dimensions don't match.
    pub fn dot(&self, other: &Self) -> i64 {
        assert_eq!(self.dim(), other.dim(), "dimension mismatch");
        self.coords
            .iter()
            .zip(other.coords.iter())
            .map(|(&a, &b)| a * b)
            .sum()
    }
}

impl From<Vec<i64>> for Point {
    fn from(coords: Vec<i64>) -> Self {
        Self::new(coords)
    }
}

impl From<&[i64]> for Point {
    fn from(coords: &[i64]) -> Self {
        Self::new(coords.to_vec())
    }
}

impl Add for &Point {
    type Output = Point;

    fn add(self, other: Self) -> Point {
        assert_eq!(self.dim(), other.dim(), "dimension mismatch");
        Point::new(
            self.coords
                .iter()
                .zip(other.coords.iter())
                .map(|(&a, &b)| a + b)
                .collect(),
        )
    }
}

impl Sub for &Point {
    type Output = Point;

    fn sub(self, other: Self) -> Point {
        assert_eq!(self.dim(), other.dim(), "dimension mismatch");
        Point::new(
            self.coords
                .iter()
                .zip(other.coords.iter())
                .map(|(&a, &b)| a - b)
                .collect(),
        )
    }
}

impl Neg for &Point {
    type Output = Point;

    fn neg(self) -> Point {
        Point::new(self.coords.iter().map(|&x| -x).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_creation() {
        let p = Point::new(vec![1, 2, 3]);
        assert_eq!(p.dim(), 3);
        assert_eq!(p.coords(), &[1, 2, 3]);
    }

    #[test]
    fn test_origin() {
        let o = Point::origin(4);
        assert!(o.is_origin());
        assert_eq!(o.dim(), 4);
    }

    #[test]
    fn test_dot_product() {
        let p1 = Point::new(vec![1, 2, 3]);
        let p2 = Point::new(vec![4, 5, 6]);
        assert_eq!(p1.dot(&p2), 4 + 10 + 18); // 1*4 + 2*5 + 3*6
    }

    #[test]
    fn test_add() {
        let p1 = Point::new(vec![1, 2]);
        let p2 = Point::new(vec![3, 4]);
        let sum = &p1 + &p2;
        assert_eq!(sum.coords(), &[4, 6]);
    }

    #[test]
    fn test_sub() {
        let p1 = Point::new(vec![5, 7]);
        let p2 = Point::new(vec![2, 3]);
        let diff = &p1 - &p2;
        assert_eq!(diff.coords(), &[3, 4]);
    }

    #[test]
    fn test_neg() {
        let p = Point::new(vec![1, -2, 3]);
        let neg_p = -&p;
        assert_eq!(neg_p.coords(), &[-1, 2, -3]);
    }

    #[test]
    fn test_from() {
        let p1: Point = vec![1, 2].into();
        assert_eq!(p1.coords(), &[1, 2]);

        let slice: &[i64] = &[3, 4];
        let p2: Point = slice.into();
        assert_eq!(p2.coords(), &[3, 4]);
    }
}
