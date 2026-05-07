//! A module for errors that can happen in series inversions.

use core::fmt;

/// An error structure for semigroup errors.
#[derive(Debug, Clone, PartialEq)]
pub enum SeriesInversionError {
    NonIntegerGVError(String),
}

/// Implement Display trait for SeriesInversionError.
impl fmt::Display for SeriesInversionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SeriesInversionError::NonIntegerGVError(details) => {
                write!(f, "A non-integer GV invariant was found: {details}")
            }
        }
    }
}
