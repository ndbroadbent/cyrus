//! Core mathematical primitives for Calabi-Yau manifold computations.
//!
//! This crate provides the foundational algorithms for:
//! - Lattice polytope operations
//! - Intersection number computation (`κ_ijk`)
//! - Kähler cone and volume calculations
//! - Flat direction and `e^{K₀}` computation
//!
//! ## Volume Computation
//!
//! The string frame volume is computed as:
//! ```text
//! V_string = (1/6) κ_ijk t^i t^j t^k - BBHL
//! ```
//!
//! where BBHL = ζ(3) χ(X) / (4(2π)³) is the α' correction.
//!
//! ## Flat Direction
//!
//! Given flux vectors K and M, the flat direction p and Kähler potential factor are:
//! ```text
//! N_ab = κ_abc M^c
//! p = N⁻¹ K
//! e^{K₀} = (4/3 × κ_abc p^a p^b p^c)⁻¹
//! ```

pub mod cosmology;
pub mod divisor;
pub mod error;
pub mod flat_direction;
pub mod glsm;
pub mod integer_math;
pub mod intersection;
pub mod kahler;
pub mod lvs;
pub use kahler::{MoriCone, compute_mori_generators};
pub mod lattice;
pub mod pipeline;
pub mod polytope;
pub use pipeline::{EvaluationRequest, EvaluationResult, evaluate_vacuum};
pub mod racetrack;
pub mod triangulation;
pub mod vacuum;
pub mod volume;

pub use divisor::{compute_divisor_jacobian, compute_divisor_volumes};
pub use error::{Error, Result};
pub use flat_direction::{
    FlatDirectionResult, compute_ek0, compute_flat_direction, compute_flat_direction_full,
    compute_n_matrix, solve_linear_system,
};
pub use glsm::compute_glsm_charge_matrix;
pub use intersection::{Intersection, compute_intersection_numbers};
// pub mod kklt;

/*
pub use kklt::{
    KkltResult, compute_c_tau as kklt_compute_c_tau, compute_divisor_volumes as kklt_compute_tau,
    compute_jacobian as kklt_compute_jacobian, compute_target_tau, solve_path_following,
};
*/
pub use lattice::Point;
pub use polytope::Polytope;
pub use racetrack::{
    GvInvariant, RacetrackResult, RacetrackTerm, build_racetrack_terms, compute_w0, solve_racetrack,
};
pub use triangulation::{Triangulation, compute_regular_triangulation};
pub use vacuum::{VacuumResult, compute_v0, compute_vacuum};
pub use volume::{VolumeResult, bbhl_correction, compute_volume, volume_classical, volume_string};
