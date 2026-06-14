//! Automated orientifold construction for the McAllister (arXiv:2107.09064)
//! model class.
//!
//! The paper's orientifolds are single-coordinate sign involutions
//! `x_i -> -x_i` with `h11_- = h21_+ = 0`. In torus-conjugacy terms an
//! involution is labelled by a parity class `sigma in (Z_2)^4 \ {0}`, and:
//!
//! - the O7-planes wrap exactly the prime toric divisors whose lattice
//!   point satisfies `p = sigma (mod 2)` componentwise (this rule is
//!   verified against the declared c_i data of all five published
//!   examples by the validation runner);
//! - each O7 carries four D7s (an so(8) stack), giving the gaugino
//!   condensate dual Coxeter number `c_i = 6`; every other rigid prime
//!   toric divisor hosts a Euclidean D3-brane with `c_i = 1`;
//! - the Freed-Witten B-field is `gamma = 1/2 [D_O7]` (parity vector);
//! - the D3 tadpole is `chi_f / 4 = (h11 + h21)/2 + 1` (paper eq. around
//!   line 1053), reproducing the published values 110/44/60/60/30.
//!
//! Rigidity of a prime toric divisor is combinatorial (Danilov-
//! Khovanskii, as used in the paper's own classification):
//! 2-face-interior points give rigid toric surfaces; 1-face-interior
//! points are P1 fibrations over genus-g curves with `g = l*(dual
//! 2-face)`, rigid iff `g = 0`; vertices have `h^{2,0} = l*(dual facet)`,
//! rigid iff that vanishes.
//!
//! KNOWN LIMITATION (documented, conservative): *pure* rigidity
//! (`h^{2,1}(D-hat) = 0`, constant Pfaffians) additionally needs
//! `h^{1,1}_-(D) = 0`, which the paper computes via the elliptic-fourfold
//! construction. This module does not verify purity; `purity_checked` is
//! false and downstream consumers must treat Pfaffians as O(1) unknowns.

use crate::Point;
use crate::divisor::{dual_face_interior_lattice_points, face_pairing_context, saturated_facets};
use crate::error::Result;
use crate::intersection::Intersection;
use crate::polytope::Polytope;

/// Combinatorial rigidity classification of one prime toric divisor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DivisorClass {
    /// Point strictly interior to a 2-face: rigid toric surface.
    RigidTwoFaceInterior,
    /// Point interior to a 1-face with genus-0 dual: rigid P1 fibration.
    RigidOneFaceInterior,
    /// Vertex with no dual-facet interior points: rigid.
    RigidVertex,
    /// Not rigid; carries the obstructing count.
    NonRigid {
        /// Interior lattice points of the dual face (genus / h^{2,0}).
        dual_interior: usize,
    },
}

impl DivisorClass {
    /// Whether this divisor supports a nonperturbative term.
    #[must_use]
    pub const fn is_rigid(&self) -> bool {
        !matches!(self, Self::NonRigid { .. })
    }
}

/// Combinatorial rigidity classification of every non-origin prime toric
/// divisor, with the irreducible-component count for non-favorable splits.
///
/// This depends only on the polytope's face combinatorics - NOT on the
/// intersection numbers, the triangulation/chamber, or the flux - so it is a
/// cheap, purely per-polytope quantity. A McAllister-class orientifold
/// requires a nonperturbative term on every divisor-basis element, hence
/// every basis divisor must be rigid; a single non-rigid basis divisor rules
/// out the whole polytope for every involution and every chamber. See
/// [`non_rigid_basis_divisors`] for the resulting pre-filter.
///
/// Returns `(classes, components)` where `classes[i] = (point_index, class)`
/// for each non-origin point and `components[idx]` is the irreducible-piece
/// count (>1 for non-favorable splits).
///
/// # Errors
/// Propagates failures from the polytope face-pairing combinatorics.
pub fn classify_divisor_rigidity(
    polytope: &Polytope,
    points: &[Point],
) -> Result<(Vec<(usize, DivisorClass)>, Vec<usize>)> {
    let (hull, dual_vertices, dual_points) = face_pairing_context(polytope, points)?;
    let mut classes: Vec<(usize, DivisorClass)> = Vec::new();
    let mut components = vec![1usize; points.len()];
    for (idx, point) in points.iter().enumerate() {
        if point.coords().iter().all(|&c| c == 0) {
            continue;
        }
        let sat = saturated_facets(point, &dual_vertices);
        let dual_interior = dual_face_interior_lattice_points(
            points,
            &hull.vertex_indices,
            &dual_points,
            &dual_vertices,
            &sat,
        )?;
        let class = match sat.len() {
            2 => {
                components[idx] = 1 + dual_interior;
                DivisorClass::RigidTwoFaceInterior
            }
            3 if dual_interior == 0 => DivisorClass::RigidOneFaceInterior,
            _ if sat.len() >= 4 && dual_interior == 0 => DivisorClass::RigidVertex,
            _ => DivisorClass::NonRigid { dual_interior },
        };
        classes.push((idx, class));
    }
    Ok((classes, components))
}

/// The `basis` divisors that are NOT rigid.
///
/// Empty means the polytope could host a McAllister-class orientifold (the
/// full O7-disjointness check still runs later); non-empty means it cannot,
/// for any involution or chamber - so the geometry can be pruned before any
/// flux search. Cheap: only the combinatorial [`classify_divisor_rigidity`].
///
/// # Errors
/// Propagates failures from [`classify_divisor_rigidity`].
pub fn non_rigid_basis_divisors(
    polytope: &Polytope,
    points: &[Point],
    basis: &[usize],
) -> Result<Vec<usize>> {
    let (classes, _) = classify_divisor_rigidity(polytope, points)?;
    Ok(basis
        .iter()
        .copied()
        .filter(|&idx| {
            classes
                .iter()
                .find(|(i, _)| *i == idx)
                .is_none_or(|(_, c)| !c.is_rigid())
        })
        .collect())
}

/// One candidate orientifold model on a fixed geometry.
#[derive(Debug, Clone)]
pub struct OrientifoldModel {
    /// Involution parity class in (Z_2)^4.
    pub sigma: Vec<i64>,
    /// Triangulation point indices of O7-plane divisors.
    pub o7_points: Vec<usize>,
    /// Rigidity classification per non-origin triangulation point.
    pub divisor_classes: Vec<(usize, DivisorClass)>,
    /// KKLT divisor set (the computed divisor basis indices).
    pub kklt_basis: Vec<usize>,
    /// Dual Coxeter / instanton coefficients per kklt_basis entry
    /// (6 for so(8) stacks on O7s, 1 for Euclidean D3-branes).
    pub c_i: Vec<i64>,
    /// Freed-Witten B-field parity vector over all triangulation points.
    pub gamma: Vec<i64>,
    /// Irreducible-component count per triangulation point (non-favorable
    /// split divisors have more than one; each component hosts its own
    /// Euclidean D3-brane).
    pub components: Vec<usize>,
    /// D3 tadpole bound chi_f/4 = (h11 + h21)/2 + 1.
    pub q_d3: f64,
    /// All O7 pairs have vanishing pairwise intersection with every divisor.
    pub o7_pairwise_disjoint: bool,
    /// Every kklt_basis divisor is rigid and the O7s are disjoint.
    pub viable: bool,
    /// Purity (h^{2,1}(D-hat) = 0) is NOT verified by this module.
    pub purity_checked: bool,
}

fn point_parity(point: &Point) -> Vec<i64> {
    point.coords().iter().map(|c| c.rem_euclid(2)).collect()
}

/// Enumerate the 15 single-coordinate-class involutions of a geometry and
/// classify each into an [`OrientifoldModel`].
///
/// `points` are the triangulation points (CYTools order, origin included),
/// `kappa_full` the full intersection tensor over those points,
/// `divisor_basis` the computed GLSM divisor basis, and `h11`/`h21` the
/// TRUE Hodge numbers (Batyrev-corrected for non-favorable polytopes).
///
/// # Errors
/// Returns an error if the combinatorial face pairing fails.
pub fn enumerate_involutions(
    polytope: &Polytope,
    points: &[Point],
    kappa_full: &Intersection,
    divisor_basis: &[usize],
    h11: usize,
    h21: usize,
) -> Result<Vec<OrientifoldModel>> {
    // Classify every non-origin point once (combinatorial; no intersection
    // numbers). `components[idx]` counts the irreducible pieces the divisor
    // splits into on X (non-favorable 2-face-interior points split into
    // 1 + l*(dual face) rigid components; see divisor.rs).
    let (classes, components) = classify_divisor_rigidity(polytope, points)?;

    let q_d3 = (h11 + h21) as f64 / 2.0 + 1.0;
    let dim = points.first().map_or(0, Point::dim);
    let mut models = Vec::new();
    for code in 1..(1u32 << dim) {
        let sigma: Vec<i64> = (0..dim).map(|bit| i64::from((code >> bit) & 1)).collect();
        let o7_points: Vec<usize> = classes
            .iter()
            .map(|(idx, _)| *idx)
            .filter(|&idx| point_parity(&points[idx]) == sigma)
            .collect();

        let gamma: Vec<i64> = points
            .iter()
            .map(|p| i64::from(!p.coords().iter().all(|&c| c == 0) && point_parity(p) == sigma))
            .collect();

        // c_i: so(8) gaugino condensation (dual Coxeter 6) on O7 stacks;
        // otherwise one Euclidean D3-brane PER irreducible component, so a
        // split divisor class gets c = number of components (verified
        // against 7-51-13590's declared c_i = 2 split divisors).
        let c_i: Vec<i64> = divisor_basis
            .iter()
            .map(|&idx| {
                if gamma[idx] == 1 {
                    6
                } else {
                    i64::try_from(components[idx]).expect("component count fits i64")
                }
            })
            .collect();

        // O7s must be mutually non-intersecting on X: D_i . D_j . D_k = 0
        // for every O7 pair (i, j) and every k.
        let mut disjoint = true;
        'pairs: for (n, &i) in o7_points.iter().enumerate() {
            for &j in &o7_points[n + 1..] {
                for k in 0..points.len() {
                    if kappa_full.get(i, j, k).to_string() != "0" {
                        disjoint = false;
                        break 'pairs;
                    }
                }
            }
        }

        let basis_all_rigid = divisor_basis.iter().all(|&idx| {
            classes
                .iter()
                .find(|(point_idx, _)| *point_idx == idx)
                .is_some_and(|(_, class)| class.is_rigid())
        });

        models.push(OrientifoldModel {
            sigma,
            o7_points,
            divisor_classes: classes.clone(),
            kklt_basis: divisor_basis.to_vec(),
            c_i,
            gamma,
            components: components.clone(),
            q_d3,
            o7_pairwise_disjoint: disjoint,
            viable: basis_all_rigid && disjoint,
            purity_checked: false,
        });
    }
    Ok(models)
}

/// Pick the best viable model: most O7-planes (strongest racetrack
/// hierarchy from so(8) stacks), then most rigid divisors.
#[must_use]
pub fn select_orientifold(models: &[OrientifoldModel]) -> Option<&OrientifoldModel> {
    models
        .iter()
        .filter(|m| m.viable && !m.o7_points.is_empty())
        .max_by_key(|m| {
            (
                m.o7_points.len(),
                m.divisor_classes
                    .iter()
                    .filter(|(_, c)| c.is_rigid())
                    .count(),
            )
        })
}
