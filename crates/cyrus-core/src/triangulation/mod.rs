//! Triangulation of lattice polytopes.
//!
//! Implements Fine, Regular, Star Triangulations (FRST) needed for
//! constructing smooth Calabi-Yau hypersurfaces.
//!
//! Uses Integer arithmetic for orientation tests (not Rational) for performance.
//! Heights are scaled by 2^40 and rounded to integers - this preserves the sign
//! of all orientation determinants which is all we need for convex hull.
//!
//! Reference: [[project_docs/CYTOOLS_ALGORITHMS_CLEAN_ROOM.md]] Section 4.

mod convex_hull;
mod heights;
mod integer_linalg;
mod regular;
mod secondary;

pub use heights::{compute_delaunay_heights, compute_frst_heights};
pub use regular::compute_regular_triangulation;
pub use secondary::{
    CircuitFlip, CircuitOmissionSide, CircuitOmissionSideClassification, circuit_omission_facets,
    circuit_triangulation_choices, classify_circuit_omission_side, complete_circuit_flip_links,
    expanded_secondary_chamber_choice_count, expanded_secondary_chamber_choice_digits,
    expanded_secondary_chamber_count_from_face_inequality_choices,
    expanded_secondary_chamber_hyperplanes_from_choice,
    expanded_secondary_chamber_hyperplanes_from_choice_index,
    expanded_secondary_chamber_hyperplanes_from_face_inequality_choice_index,
    expanded_secondary_chamber_hyperplanes_from_face_triangulation_choice_index,
    expanded_secondary_chamber_hyperplanes_from_face_triangulation_choice_index_with_star_4d,
    expanded_secondary_cone_hyperplanes_from_face_triangulations,
    expanded_secondary_cone_hyperplanes_from_face_triangulations_with_star_4d,
    expanded_secondary_face_inequality_choices_from_circuit_faces,
    expanded_secondary_face_inequality_choices_from_triangulations,
    expanded_secondary_face_inequality_choices_from_triangulations_with_star_4d,
    expanded_secondary_fan_hyperplanes_on_faces,
    expanded_secondary_fan_hyperplanes_on_polytope_2faces_4d,
    expanded_secondary_group_boring_chamber_choices,
    expanded_secondary_regular_triangulation_from_face_triangulations,
    expanded_secondary_star_hyperplanes_for_face_triangulation_4d, flip_circuit_in_triangulation,
    flip_circuit_link_in_triangulation, secondary_cone_height_pairings,
    secondary_cone_hyperplanes_native, secondary_cone_hyperplanes_native_on_faces,
    secondary_cone_hyperplanes_native_on_polytope_2faces_4d,
    secondary_cone_strictly_contains_height_vector, triangulation_heights_from_secondary_cone,
    triangulation_is_regular_from_secondary_cone,
};

/// A triangulation of a point set.
///
/// Represented as a collection of simplices, where each simplex is a
/// set of indices into the original point set.
#[derive(Debug, Clone)]
pub struct Triangulation {
    /// Indices of points forming each simplex.
    /// For a d-dimensional triangulation, each simplex has d+1 indices.
    simplices: Vec<Vec<usize>>,
}

impl Triangulation {
    /// Create a new triangulation from a list of simplices.
    pub const fn new(simplices: Vec<Vec<usize>>) -> Self {
        Self { simplices }
    }

    /// Get the simplices.
    pub fn simplices(&self) -> &[Vec<usize>] {
        &self.simplices
    }

    /// Check if this triangulation has the Star property.
    ///
    /// A triangulation is Star if a specified point (typically the origin)
    /// is a vertex of every simplex.
    pub fn is_star(&self, point_idx: usize) -> bool {
        self.simplices.iter().all(|s| s.contains(&point_idx))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_star() {
        // Triangulation where point 0 is in all simplices
        let tri_star = Triangulation::new(vec![vec![0, 1, 2], vec![0, 2, 3], vec![0, 3, 1]]);
        assert!(tri_star.is_star(0));
        assert!(!tri_star.is_star(1));

        // Triangulation where point 0 is missing from one simplex
        let tri_not_star = Triangulation::new(vec![
            vec![0, 1, 2],
            vec![1, 2, 3], // Missing point 0
        ]);
        assert!(!tri_not_star.is_star(0));
    }
}
