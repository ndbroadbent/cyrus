//! Utility functions for CYTools port.
//!
//! This module contains common utility functions used throughout cyrus-core,
//! ported from CYTools `utils.py`.

mod affine;
mod gcd;
mod linear_solve;
mod lll;
mod rational;
mod tensor;

pub use affine::find_new_affinely_independent_points;
pub use gcd::{gcd_float, gcd_float_default, gcd_int, gcd_list, gcd_list_default, gcd_list_int};
pub use linear_solve::{LinearBackend, solve_linear_system, solve_linear_system_sparse};
pub use lll::{LllResult, lll_reduce};
pub use rational::{
    array_from_rational, array_to_rational, float_to_rational, gcd_integers, integral_nullspace,
    rational_to_float,
};
pub use tensor::{
    SparseTensor, filter_tensor_indices, symmetric_dense_to_sparse, symmetric_dense_to_sparse_3d,
    symmetric_sparse_to_dense, symmetric_sparse_to_dense_3d,
};
