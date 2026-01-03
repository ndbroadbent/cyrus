//! Type-safe range for compile-time verified iteration bounds.
//!
//! `CheckedRange` ensures range validity at compile time and provides
//! typed iteration over `I64` values.

use super::i64::I64;
use super::tags::{NonNeg, Pos};

/// A compile-time checked range with start < end guarantee.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CheckedRange<T> {
    /// The start of the range (inclusive).
    pub start: T,
    /// The end of the range (exclusive).
    pub end: T,
}

impl<T> CheckedRange<T> {
    /// Create a new checked range.
    #[must_use]
    pub const fn new(start: T, end: T) -> Self {
        Self { start, end }
    }
}

/// Create a `CheckedRange<usize>` with compile-time verification that start < end.
#[macro_export]
macro_rules! range {
    ($start:expr, $end:expr) => {{
        const S: usize = $start;
        const E: usize = $end;
        const _: () = {
            if !(S < E) {
                panic!("range! requires start < end");
            }
        };
        $crate::types::range::CheckedRange::new(S, E)
    }};
}

// ============================================================================
// Iterator for CheckedRange<usize> yielding I64<NonNeg>
// ============================================================================

impl CheckedRange<usize> {
    /// Iterate yielding `I64<NonNeg>` values (0, 1, 2, ...).
    #[must_use]
    pub fn iter_non_neg(self) -> CheckedRangeIterNonNeg {
        CheckedRangeIterNonNeg {
            current: self.start,
            end: self.end,
        }
    }

    /// Iterate yielding `I64<Pos>` values (1, 2, 3, ...).
    /// Skips start if start == 0.
    #[must_use]
    pub fn iter_pos(self) -> CheckedRangeIterPos {
        let current = if self.start == 0 { 1 } else { self.start };
        CheckedRangeIterPos {
            current,
            end: self.end,
        }
    }
}

/// Iterator yielding `I64<NonNeg>`.
#[derive(Clone, Debug)]
pub struct CheckedRangeIterNonNeg {
    current: usize,
    end: usize,
}

impl Iterator for CheckedRangeIterNonNeg {
    type Item = I64<NonNeg>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.end {
            let val = I64::<NonNeg>::from_raw(self.current as i64);
            self.current += 1;
            Some(val)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.end.saturating_sub(self.current);
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for CheckedRangeIterNonNeg {}

/// Iterator yielding `I64<Pos>`.
#[derive(Clone, Debug)]
pub struct CheckedRangeIterPos {
    current: usize,
    end: usize,
}

impl Iterator for CheckedRangeIterPos {
    type Item = I64<Pos>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.end && self.current > 0 {
            let val = I64::<Pos>::from_raw(self.current as i64);
            self.current += 1;
            Some(val)
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let remaining = self.end.saturating_sub(self.current);
        (remaining, Some(remaining))
    }
}

impl ExactSizeIterator for CheckedRangeIterPos {}

// ============================================================================
// IntoIterator for CheckedRange - default yields I64<NonNeg>
// ============================================================================

impl IntoIterator for CheckedRange<usize> {
    type Item = I64<NonNeg>;
    type IntoIter = CheckedRangeIterNonNeg;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_non_neg()
    }
}

#[cfg(test)]
mod tests {
    use crate::range;

    #[test]
    fn test_range_macro() {
        let r = range!(0, 10);
        assert_eq!(r.start, 0);
        assert_eq!(r.end, 10);
    }

    #[test]
    fn test_iter_non_neg() {
        let r = range!(0, 3);
        let vals: Vec<_> = r.iter_non_neg().map(super::super::i64::I64::get).collect();
        assert_eq!(vals, vec![0, 1, 2]);
    }

    #[test]
    fn test_iter_pos() {
        let r = range!(1, 4);
        let vals: Vec<_> = r.iter_pos().map(super::super::i64::I64::get).collect();
        assert_eq!(vals, vec![1, 2, 3]);
    }

    #[test]
    fn test_iter_pos_from_zero() {
        let r = range!(0, 3);
        let vals: Vec<_> = r.iter_pos().map(super::super::i64::I64::get).collect();
        assert_eq!(vals, vec![1, 2]); // skips 0
    }

    #[test]
    fn test_into_iter() {
        let r = range!(0, 3);
        let vals: Vec<_> = r.into_iter().map(super::super::i64::I64::get).collect();
        assert_eq!(vals, vec![0, 1, 2]);
    }
}
