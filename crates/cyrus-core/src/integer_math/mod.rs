//! Integer linear algebra primitives.
//!
//! Provides algorithms for working with integer matrices, such as
//! Hermite Normal Form (HNF) and integer kernel computation.
//!
//! These are essential for finding the integral basis of the GLSM charge matrix.
//!
//! ## Modules
//!
//! - [`determinant`] - Determinant and orientation computations
//! - [`kernel`] - Integer kernel (nullspace) computation
//! - [`solve`] - Linear system solving for rational matrices
//! - [`matrix_utils`] - Matrix utilities (rank, inversion, GCD/LCM)
//! - [`linear_relations`] - GLSM and linear relations from polytope points

mod determinant;
mod hnf;
mod kernel;
mod linear_relations;
mod matrix_utils;
mod snf;
mod solve;

#[cfg(test)]
mod tests;

// Re-export all public functions
pub use determinant::{determinant_gaussian, orientation};
pub use hnf::{
    hermite_normal_form, hermite_normal_form_columns, is_row_hnf, pivot_product, sublattice_index,
};
pub use kernel::integer_kernel;
pub use linear_relations::{compute_glsm_and_linear_relations, compute_linear_relations_no_origin};
pub use matrix_utils::{
    gcd_integer, invert_matrix, lcm_integer, matrix_rank, rational_matrix_to_integer,
};
pub use snf::{smith_normal_form_diag, sublattice_index_snf};
pub use solve::solve_linear_system_rational;
