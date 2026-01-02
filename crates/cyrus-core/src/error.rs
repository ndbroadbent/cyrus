//! Error types for cyrus-core.

use thiserror::Error;

/// Result type alias for cyrus-core operations.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors that can occur in cyrus-core operations.
#[derive(Error, Debug)]
pub enum Error {
    /// The polytope is not reflexive.
    #[error("polytope is not reflexive: {0}")]
    NotReflexive(String),

    /// Invalid dimension for the operation.
    #[error("invalid dimension: expected {expected}, got {got}")]
    InvalidDimension {
        /// Expected dimension.
        expected: usize,
        /// Actual dimension.
        got: usize,
    },

    /// Linear algebra operation failed.
    #[error("linear algebra error: {0}")]
    LinearAlgebra(String),

    /// Singular matrix encountered in solver.
    #[error("singular matrix: {0}")]
    SingularMatrix(String),

    /// Invalid input data.
    #[error("invalid input: {0}")]
    InvalidInput(String),
}
