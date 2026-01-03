//! Compile-time validated range types.
//!
//! This module provides [`CheckedRange`] for ranges that are validated at
//! compile time using const assertions.
//!
//! # Example
//!
//! ```
//! use cyrus_core::checked_range;
//! use cyrus_core::types::range::CheckedRange;
//!
//! // Valid range - compiles
//! const VALID: CheckedRange<usize> = checked_range!(3..10);
//!
//! // Invalid range - would not compile:
//! // const INVALID: CheckedRange<usize> = checked_range!(10..3);
//! ```

use std::ops::Range;

/// A range that has been validated at construction time.
///
/// Unlike `std::ops::Range`, a `CheckedRange` guarantees that `start < end`.
/// Use the [`checked_range!`] macro for compile-time validation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct CheckedRange<T> {
    start: T,
    end: T,
}

impl<T> CheckedRange<T> {
    /// Create a new checked range.
    ///
    /// This is a const fn for use in const contexts, but doesn't validate.
    /// Use [`checked_range!`] macro for compile-time validation.
    #[must_use]
    pub const fn new_unchecked(start: T, end: T) -> Self {
        Self { start, end }
    }

    /// Get the start of the range.
    #[must_use]
    pub const fn start(&self) -> &T {
        &self.start
    }

    /// Get the end of the range (exclusive).
    #[must_use]
    pub const fn end(&self) -> &T {
        &self.end
    }
}

impl<T: Copy> CheckedRange<T> {
    /// Get the start value.
    #[must_use]
    pub const fn start_val(&self) -> T {
        self.start
    }

    /// Get the end value.
    #[must_use]
    pub const fn end_val(&self) -> T {
        self.end
    }
}

impl<T: PartialOrd> CheckedRange<T> {
    /// Try to create a checked range. Returns `None` if `start >= end`.
    #[must_use]
    pub fn try_new(start: T, end: T) -> Option<Self> {
        if start < end {
            Some(Self { start, end })
        } else {
            None
        }
    }

    /// Create a checked range. Panics if `start >= end`.
    ///
    /// # Panics
    ///
    /// Panics if `start >= end`.
    #[must_use]
    pub fn new(start: T, end: T) -> Self {
        assert!(start < end, "CheckedRange requires start < end");
        Self { start, end }
    }

    /// Check if a value is within the range.
    #[must_use]
    pub fn contains(&self, value: &T) -> bool {
        *value >= self.start && *value < self.end
    }
}

impl CheckedRange<usize> {
    /// Get the number of elements in the range.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.end - self.start
    }

    /// Check if the range is empty.
    ///
    /// Note: A `CheckedRange` is never empty by construction (start < end),
    /// but this method is provided for API completeness.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        // A checked range is never empty, but we implement this for clippy
        self.start >= self.end
    }

    /// Iterate over the range.
    pub fn iter(&self) -> impl Iterator<Item = usize> {
        self.start..self.end
    }
}

impl CheckedRange<i64> {
    /// Get the number of elements in the range.
    #[must_use]
    pub const fn len(&self) -> usize {
        (self.end - self.start) as usize
    }

    /// Check if the range is empty.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.start >= self.end
    }

    /// Iterate over the range.
    pub fn iter(&self) -> impl Iterator<Item = i64> {
        self.start..self.end
    }
}

impl CheckedRange<f64> {
    /// Get the width of the range.
    #[must_use]
    pub fn width(&self) -> f64 {
        self.end - self.start
    }

    /// Get the midpoint of the range.
    #[must_use]
    pub fn midpoint(&self) -> f64 {
        (self.start + self.end) / 2.0
    }

    /// Linearly interpolate within the range.
    ///
    /// `t = 0.0` gives `start`, `t = 1.0` gives `end`.
    #[must_use]
    pub fn lerp(&self, t: f64) -> f64 {
        self.start + t * (self.end - self.start)
    }

    /// Clamp a value to the range.
    #[must_use]
    pub fn clamp(&self, value: f64) -> f64 {
        value.clamp(self.start, self.end)
    }
}

// Convert to std Range
impl<T: Clone> From<CheckedRange<T>> for Range<T> {
    fn from(r: CheckedRange<T>) -> Self {
        r.start.clone()..r.end.clone()
    }
}

impl<T: Copy> From<&CheckedRange<T>> for Range<T> {
    fn from(r: &CheckedRange<T>) -> Self {
        r.start..r.end
    }
}

// IntoIterator for usize ranges
impl IntoIterator for CheckedRange<usize> {
    type Item = usize;
    type IntoIter = std::ops::Range<usize>;

    fn into_iter(self) -> Self::IntoIter {
        self.start..self.end
    }
}

impl IntoIterator for &CheckedRange<usize> {
    type Item = usize;
    type IntoIter = std::ops::Range<usize>;

    fn into_iter(self) -> Self::IntoIter {
        self.start..self.end
    }
}

// IntoIterator for i64 ranges
impl IntoIterator for CheckedRange<i64> {
    type Item = i64;
    type IntoIter = std::ops::Range<i64>;

    fn into_iter(self) -> Self::IntoIter {
        self.start..self.end
    }
}

/// Create a [`CheckedRange<usize>`] with compile-time validation.
///
/// # Example
///
/// ```
/// use cyrus_core::checked_range;
/// use cyrus_core::types::range::CheckedRange;
///
/// const R: CheckedRange<usize> = checked_range!(0, 10);
/// assert_eq!(R.len(), 10);
///
/// // This would fail at compile time:
/// // const BAD: CheckedRange<usize> = checked_range!(10, 5);
/// ```
#[macro_export]
macro_rules! checked_range {
    ($start:expr, $end:expr) => {{
        const S: usize = $start;
        const E: usize = $end;
        const _: () = assert!(S < E, "checked_range! requires start < end");
        $crate::types::range::CheckedRange::new_unchecked(S, E)
    }};
}

/// Create a [`CheckedRange<i64>`] with compile-time validation.
#[macro_export]
macro_rules! checked_range_i64 {
    ($start:expr, $end:expr) => {{
        const S: i64 = $start;
        const E: i64 = $end;
        const _: () = assert!(S < E, "checked_range_i64! requires start < end");
        $crate::types::range::CheckedRange::new_unchecked(S, E)
    }};
}

/// Create a [`CheckedRange<f64>`] with compile-time validation.
#[macro_export]
macro_rules! checked_range_f64 {
    ($start:expr, $end:expr) => {{
        const S: f64 = $start;
        const E: f64 = $end;
        // Can't use < for f64 in const context directly, use this trick
        const _: () = assert!(
            S < E && S == S && E == E, // Also checks for NaN
            "checked_range_f64! requires start < end and no NaN"
        );
        $crate::types::range::CheckedRange::new_unchecked(S, E)
    }};
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_checked_range_new() {
        let r = CheckedRange::new(1usize, 5usize);
        assert_eq!(*r.start(), 1);
        assert_eq!(*r.end(), 5);
    }

    #[test]
    #[should_panic(expected = "start < end")]
    fn test_checked_range_new_invalid() {
        let _r = CheckedRange::new(5usize, 1usize);
    }

    #[test]
    fn test_checked_range_try_new() {
        assert!(CheckedRange::try_new(1usize, 5usize).is_some());
        assert!(CheckedRange::try_new(5usize, 1usize).is_none());
        assert!(CheckedRange::try_new(3usize, 3usize).is_none());
    }

    #[test]
    fn test_checked_range_contains() {
        let r = CheckedRange::new(2usize, 5usize);
        assert!(!r.contains(&1));
        assert!(r.contains(&2));
        assert!(r.contains(&3));
        assert!(r.contains(&4));
        assert!(!r.contains(&5)); // exclusive end
        assert!(!r.contains(&6));
    }

    #[test]
    fn test_checked_range_len() {
        let r = CheckedRange::new(3usize, 10usize);
        assert_eq!(r.len(), 7);
        assert!(!r.is_empty());
    }

    #[test]
    fn test_checked_range_iter() {
        let r = CheckedRange::new(1usize, 4usize);
        let v: Vec<usize> = r.iter().collect();
        assert_eq!(v, vec![1, 2, 3]);
    }

    #[test]
    fn test_checked_range_into_iter() {
        let r = CheckedRange::new(1usize, 4usize);
        let v: Vec<usize> = r.into_iter().collect();
        assert_eq!(v, vec![1, 2, 3]);
    }

    #[test]
    fn test_checked_range_ref_into_iter() {
        let r = CheckedRange::new(1usize, 4usize);
        let v: Vec<usize> = (&r).into_iter().collect();
        assert_eq!(v, vec![1, 2, 3]);
    }

    #[test]
    fn test_checked_range_i64() {
        let r = CheckedRange::new(-3i64, 2i64);
        assert_eq!(r.len(), 5);
        let v: Vec<i64> = r.iter().collect();
        assert_eq!(v, vec![-3, -2, -1, 0, 1]);
    }

    #[test]
    fn test_checked_range_i64_into_iter() {
        let r = CheckedRange::new(-1i64, 2i64);
        let v: Vec<i64> = r.into_iter().collect();
        assert_eq!(v, vec![-1, 0, 1]);
    }

    #[test]
    fn test_checked_range_f64() {
        let r = CheckedRange::new(0.0f64, 1.0f64);
        assert!((r.width() - 1.0).abs() < 1e-10);
        assert!((r.midpoint() - 0.5).abs() < 1e-10);
        assert!((r.lerp(0.0) - 0.0).abs() < 1e-10);
        assert!((r.lerp(1.0) - 1.0).abs() < 1e-10);
        assert!((r.lerp(0.5) - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_checked_range_f64_clamp() {
        let r = CheckedRange::new(0.0f64, 1.0f64);
        assert!((r.clamp(-1.0) - 0.0).abs() < 1e-10);
        assert!((r.clamp(0.5) - 0.5).abs() < 1e-10);
        assert!((r.clamp(2.0) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_checked_range_from() {
        let r = CheckedRange::new(1usize, 5usize);
        let std_range: Range<usize> = r.into();
        assert_eq!(std_range, 1..5);
    }

    #[test]
    fn test_checked_range_ref_from() {
        let r = CheckedRange::new(1usize, 5usize);
        let std_range: Range<usize> = (&r).into();
        assert_eq!(std_range, 1..5);
    }

    #[test]
    fn test_checked_range_start_end_val() {
        let r = CheckedRange::new(2usize, 8usize);
        assert_eq!(r.start_val(), 2);
        assert_eq!(r.end_val(), 8);
    }

    #[test]
    fn test_checked_range_macro() {
        const R: CheckedRange<usize> = checked_range!(3, 10);
        assert_eq!(R.start_val(), 3);
        assert_eq!(R.end_val(), 10);
        assert_eq!(R.len(), 7);
    }

    #[test]
    fn test_checked_range_i64_macro() {
        const R: CheckedRange<i64> = checked_range_i64!(-5, 5);
        assert_eq!(R.start_val(), -5);
        assert_eq!(R.end_val(), 5);
        assert_eq!(R.len(), 10);
    }

    #[test]
    fn test_checked_range_f64_macro() {
        const R: CheckedRange<f64> = checked_range_f64!(0.0, 1.0);
        assert!((R.width() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_checked_range_hash() {
        use std::collections::HashSet;
        let r1 = CheckedRange::new(1usize, 5usize);
        let r2 = CheckedRange::new(1usize, 5usize);
        let r3 = CheckedRange::new(2usize, 5usize);

        let mut set = HashSet::new();
        set.insert(r1);
        set.insert(r2);
        set.insert(r3);
        assert_eq!(set.len(), 2);
    }

    #[test]
    fn test_checked_range_eq() {
        let r1 = CheckedRange::new(1usize, 5usize);
        let r2 = CheckedRange::new(1usize, 5usize);
        let r3 = CheckedRange::new(2usize, 5usize);

        assert_eq!(r1, r2);
        assert_ne!(r1, r3);
    }
}
