// Clippy allows for intentional patterns in physics code
#![allow(clippy::cast_precision_loss)] // i64 to f64 is intentional for floating point math
#![allow(clippy::cast_possible_wrap)] // usize to i64 is safe for our data sizes
#![allow(clippy::cast_possible_truncation)] // Safe truncations are checked at runtime
#![allow(clippy::cast_sign_loss)] // Safe sign conversions are checked at runtime
#![allow(clippy::missing_panics_doc)] // Not all panic conditions need documentation
#![allow(clippy::missing_errors_doc)] // Not all error conditions need documentation
#![allow(clippy::module_name_repetitions)] // Types like `ConeError` in `cone` module are fine
#![allow(clippy::needless_range_loop)] // Matrix operations are clearer with index loops
#![allow(clippy::branches_sharing_code)] // Sometimes clearer to repeat code in branches
#![allow(clippy::must_use_candidate)] // Not all methods need #[must_use]
#![allow(clippy::if_then_some_else_none)] // Pattern is fine
#![allow(clippy::items_after_statements)] // Test helper structs in functions are fine
#![allow(clippy::type_complexity)] // Complex types are acceptable in physics code
#![allow(clippy::single_match_else)] // Match with single arm is sometimes clearer
#![allow(clippy::match_same_arms)] // Sometimes clearer to be explicit about same behavior
#![allow(clippy::option_if_let_else)] // Pattern is often clearer than map_or_else
#![allow(clippy::needless_collect)] // Sometimes collect is clearer
#![allow(clippy::only_used_in_recursion)] // False positives in some algorithms
#![allow(clippy::redundant_else)] // Explicit else for clarity
#![allow(clippy::manual_let_else)] // let...else not always clearer
#![allow(clippy::cmp_owned)] // Owned comparison is sometimes clearer
#![allow(clippy::return_self_not_must_use)] // Methods returning Self don't always need must_use
#![allow(clippy::no_effect_underscore_binding)] // Debug variables are fine

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

pub mod basis;
pub mod cone;
pub mod config;
pub mod cosmology;
pub mod curve_basis;
pub mod divisor;
pub mod error;
pub mod flat_direction;
pub mod geometry;
pub mod glsm;
pub mod gv;
pub mod integer_math;
pub mod intersection;
pub mod kahler;
pub mod lvs;
pub use cone::Cone;
pub use kahler::{MoriCone, compute_mori_generators};
pub mod lattice;
pub mod pipeline;
pub mod policy;
pub mod polytope;
pub mod quintessence;
pub use pipeline::{EvaluationRequest, EvaluationResult, evaluate_vacuum};
pub use policy::{Abort, ForGA, Strict, VacuumPolicy, VolumePolicy};
pub mod racetrack;
pub mod triangulation;
pub mod types;
pub mod utils;
pub mod vacuum;
pub mod volume;

pub use types::{F64, Finite, H11, H21, I32, I64, Neg, Pos};

pub use basis::{basis_change_matrix, compute_divisor_basis, intersection_in_basis, is_unimodular};
pub use curve_basis::compute_curve_basis_matrix;
pub use divisor::{compute_divisor_jacobian, compute_divisor_volumes, compute_kklt_divisor_chi};
pub use error::{Error, Result};
pub use flat_direction::{
    FlatDirectionResult, compute_ek0, compute_flat_direction, compute_flat_direction_full,
    compute_n_matrix, solve_linear_system,
};
pub use glsm::compute_glsm_and_linrels;
pub use glsm::compute_glsm_charge_matrix;
pub use glsm::compute_glsm_linear_relations;
pub use gv::{
    ToricCurveCandidate, ToricCurveGvInvariant, compute_grading_vector, compute_gv_invariants,
    compute_mori_cone_cap_rays, compute_toric_two_face_curve_gv_invariants,
    curve_volume_in_divisor_basis, find_pair_decomposition,
    remove_pair_decomposable_curve_candidates, subcutoff_toric_curve_candidates,
};
pub use integer_math::{compute_glsm_and_linear_relations, compute_linear_relations_no_origin};
pub use intersection::compute_intersection_cytools;
pub use intersection::{
    Intersection, compute_intersection_numbers, compute_intersection_numbers_with_linear_relations,
    compute_intersection_numbers_with_offset,
};
pub mod kklt;

pub use kklt::{
    KkltResult, compute_c_tau as kklt_compute_c_tau, compute_corrected_target_tau,
    compute_divisor_volumes as kklt_compute_tau, compute_gv_target_correction_for_ambient_curves,
    compute_jacobian as kklt_compute_jacobian, compute_kklt_divisor_volumes, compute_kklt_jacobian,
    compute_target_tau, solve_mixed_basis_path_following, solve_path_following,
    solve_two_phase_mixed_basis_path_following, solve_two_phase_path_following,
};
pub use lattice::Point;
pub use polytope::Polytope;
pub use racetrack::{
    GvInvariant, RacetrackResult, RacetrackTerm, build_racetrack_terms, compute_w0,
    compute_w0_from_terms, solve_racetrack,
};
pub use triangulation::{
    Triangulation, compute_delaunay_heights, compute_frst_heights, compute_regular_triangulation,
};
pub use vacuum::{VacuumResult, compute_v0, compute_vacuum};
pub use volume::{VolumeResult, bbhl_correction, compute_volume, volume_classical, volume_string};
