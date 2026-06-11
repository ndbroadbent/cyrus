//! Divisor volume computations.
//!
//! The volume of a divisor `D_i` in a Calabi-Yau threefold is:
//! ```text
//! τ_i = (1/2) κ_ijk t^j t^k
//! ```
//!
//! The Jacobian `∂τ_i/∂t^m` is:
//! ```text
//! J_{im} = κ_ijm t^j
//! ```
//!
//! Reference: arXiv:2107.09064, Section 5

use malachite::num::conversion::traits::RoundingFrom;
use malachite::rounding_modes::RoundingMode;
use std::collections::HashSet;

use crate::error::{Error, Result};
use crate::f64_pos;
use crate::geometry::ConvexHull;
use crate::intersection::Intersection;
use crate::lattice::Point;
use crate::polytope::Polytope;
use crate::types::f64::F64;
use crate::types::i64::I64;
use crate::types::tags::Finite;

/// Compute divisor volumes `τ_i` = (1/2) `κ_ijk` t^j t^k.
///
/// # Panics
/// Panics if `t.len() != kappa.dim()`.
#[must_use]
pub fn compute_divisor_volumes(kappa: &Intersection, t: &[F64<Finite>]) -> Vec<F64<Finite>> {
    let dim = kappa.dim();
    assert_eq!(t.len(), dim, "dimension mismatch");

    let half = f64_pos!(0.5);
    let two = f64_pos!(2.0);

    // Initialize tau vector
    let mut tau = vec![F64::<Finite>::ZERO; dim];

    // Build the contractions manually
    for ((a, b, c), val) in kappa.iter() {
        let (val_raw, _) = f64::rounding_from(val.get(), RoundingMode::Nearest);
        let val_f = F64::<Finite>::new(val_raw).expect("intersection value is finite");

        // Contribution to τ_a from κ_{abc} t^b t^c
        // κ_{abc} contributes to τ_a, τ_b, τ_c based on symmetry

        if a == b && b == c {
            // All same: κ_{aaa} t^a t^a - multiplicity 1
            tau[*a] = tau[*a] + half * val_f * t[*a] * t[*a];
        } else if a == b {
            // (a,a,c) with a < c - multiplicity 3
            tau[*a] = tau[*a] + half * val_f * two * t[*a] * t[*c];
            tau[*c] = tau[*c] + half * val_f * t[*a] * t[*a];
        } else if b == c {
            // (a,b,b) with a < b - multiplicity 3
            tau[*a] = tau[*a] + half * val_f * t[*b] * t[*b];
            tau[*b] = tau[*b] + half * val_f * two * t[*a] * t[*b];
        } else {
            // All different (a,b,c) with a < b < c - multiplicity 6
            tau[*a] = tau[*a] + half * val_f * two * t[*b] * t[*c];
            tau[*b] = tau[*b] + half * val_f * two * t[*a] * t[*c];
            tau[*c] = tau[*c] + half * val_f * two * t[*a] * t[*b];
        }
    }

    tau
}

/// Compute Jacobian of divisor volumes: `J_{im}` = `∂τ_i/∂t^m` = `κ_ijm` t^j.
///
/// Returns a `dim × dim` matrix where `J[i][m] = ∂τ_i/∂t^m`.
///
/// # Panics
/// Panics if `t.len() != kappa.dim()`.
#[must_use]
pub fn compute_divisor_jacobian(kappa: &Intersection, t: &[F64<Finite>]) -> Vec<Vec<F64<Finite>>> {
    let dim = kappa.dim();
    assert_eq!(t.len(), dim, "dimension mismatch");

    let mut jac = vec![vec![F64::<Finite>::ZERO; dim]; dim];

    // J[i][m] = κ_{ijm} t^j (summing over j)
    // From canonical (a,b,c) with a ≤ b ≤ c, we need to fill J entries

    for ((a, b, c), val) in kappa.iter() {
        let (val_raw, _) = f64::rounding_from(val.get(), RoundingMode::Nearest);
        let val_f = F64::<Finite>::new(val_raw).expect("intersection value is finite");

        if a == b && b == c {
            // (a,a,a): J[a][a] += κ_{aaa} t^a
            jac[*a][*a] = jac[*a][*a] + val_f * t[*a];
        } else if a == b {
            // (a,a,c) with a < c
            // J[a][a] += κ_{aac} t^c
            // J[a][c] += κ_{aac} t^a (from j=a)
            // J[c][a] += κ_{caa} t^a = κ_{aac} t^a (symmetric, from j=a)
            jac[*a][*a] = jac[*a][*a] + val_f * t[*c];
            jac[*a][*c] = jac[*a][*c] + val_f * t[*a];
            jac[*c][*a] = jac[*c][*a] + val_f * t[*a];
        } else if b == c {
            // (a,b,b) with a < b
            // J[a][b] += κ_{abb} t^b (from j=b)
            // J[b][a] += κ_{bba} t^b = κ_{abb} t^b (symmetric)
            // J[b][b] += κ_{bba} t^a = κ_{abb} t^a
            jac[*a][*b] = jac[*a][*b] + val_f * t[*b];
            jac[*b][*a] = jac[*b][*a] + val_f * t[*b];
            jac[*b][*b] = jac[*b][*b] + val_f * t[*a];
        } else {
            // (a,b,c) all different
            // J[a][b] += κ_{abc} t^c (from κ_{ajb} with j=c)
            // J[a][c] += κ_{abc} t^b (from κ_{ajc} with j=b)
            // J[b][a] += κ_{bac} t^c = κ_{abc} t^c
            // J[b][c] += κ_{bac} t^a = κ_{abc} t^a
            // J[c][a] += κ_{cab} t^b = κ_{abc} t^b
            // J[c][b] += κ_{cab} t^a = κ_{abc} t^a
            jac[*a][*b] = jac[*a][*b] + val_f * t[*c];
            jac[*a][*c] = jac[*a][*c] + val_f * t[*b];
            jac[*b][*a] = jac[*b][*a] + val_f * t[*c];
            jac[*b][*c] = jac[*b][*c] + val_f * t[*a];
            jac[*c][*a] = jac[*c][*a] + val_f * t[*b];
            jac[*c][*b] = jac[*c][*b] + val_f * t[*a];
        }
    }

    jac
}

fn dot(a: &[i64], b: &[i64]) -> i64 {
    a.iter().zip(b.iter()).map(|(&x, &y)| x * y).sum()
}

fn integral_rational_to_i64(
    value: &crate::types::rational::Rational<Finite>,
    label: &str,
) -> Result<i64> {
    if value.get().denominator_ref() != &1u32 {
        return Err(Error::InvalidInput(format!("{label} is not integral")));
    }
    let signed = malachite::Integer::try_from(value.get().clone())
        .map_err(|_| Error::InvalidInput(format!("{label} is not integral")))?;
    i64::try_from(&signed).map_err(|_| Error::InvalidInput(format!("{label} does not fit in i64")))
}

fn saturated_facets(point: &Point, dual_vertices: &[Point]) -> Vec<usize> {
    dual_vertices
        .iter()
        .enumerate()
        .filter_map(|(idx, dual_vertex)| {
            (dot(point.coords(), dual_vertex.coords()) == -1).then_some(idx)
        })
        .collect()
}

fn dual_face_interior_lattice_points(
    points: &[Point],
    hull_vertices: &[usize],
    dual_points: &[Point],
    dual_vertices: &[Point],
    point_saturated_facets: &[usize],
) -> Result<usize> {
    let face_vertices: HashSet<usize> = hull_vertices
        .iter()
        .copied()
        .filter(|&point_idx| {
            point_saturated_facets.iter().all(|&facet_idx| {
                dot(
                    points[point_idx].coords(),
                    dual_vertices[facet_idx].coords(),
                ) == -1
            })
        })
        .collect();

    if face_vertices.is_empty() {
        return Err(Error::InvalidInput(
            "could not determine primal face vertices for divisor point".into(),
        ));
    }

    let mut count = 0usize;
    'dual_points: for dual_point in dual_points {
        for &point_idx in &face_vertices {
            if dot(dual_point.coords(), points[point_idx].coords()) != -1 {
                continue 'dual_points;
            }
        }
        for &point_idx in hull_vertices {
            if face_vertices.contains(&point_idx) {
                continue;
            }
            if dot(dual_point.coords(), points[point_idx].coords()) == -1 {
                continue 'dual_points;
            }
        }
        count += 1;
    }

    Ok(count)
}

fn face_pairing_context(
    polytope: &Polytope,
    points: &[Point],
) -> Result<(ConvexHull, Vec<Point>, Vec<Point>)> {
    let point_coords: Vec<Vec<i64>> = points.iter().map(|point| point.coords().to_vec()).collect();
    let hull = ConvexHull::compute(&point_coords)
        .ok_or_else(|| Error::InvalidInput("failed to compute primal convex hull".into()))?;
    let dual_vertices = polytope.dual_vertices()?;
    let dual_polytope = polytope.compute_dual()?;
    let dual_points = dual_polytope.vertices().to_vec();
    Ok((hull, dual_vertices, dual_points))
}

/// Count the extra divisor classes of a non-favorable polytope.
///
/// Batyrev's correction term beyond the toric GLSM rank:
/// `h11 = n_toric_classes + sum over 2-faces of l*(face) * l*(dual face)`.
/// Each point strictly interior to a 2-face whose dual 1-face has interior
/// lattice points contributes that interior count; favorable polytopes get 0.
pub fn batyrev_h11_extra_classes(polytope: &Polytope, points: &[Point]) -> Result<usize> {
    if polytope.dim() != 4 {
        return Err(Error::InvalidInput(
            "Batyrev h11 correction is currently implemented for 4D polytopes".into(),
        ));
    }
    let (hull, dual_vertices, dual_points) = face_pairing_context(polytope, points)?;

    let mut extra = 0usize;
    for point in points {
        if point.coords().iter().all(|&coord| coord == 0) {
            continue;
        }
        let sat = saturated_facets(point, &dual_vertices);
        if sat.len() != 2 {
            continue;
        }
        extra += dual_face_interior_lattice_points(
            points,
            &hull.vertex_indices,
            &dual_points,
            &dual_vertices,
            &sat,
        )?;
    }
    Ok(extra)
}

/// Compute topological divisor Euler characteristics for KKLT divisors.
///
/// Braun et al. combinatorial formula: `chi(D) = 12 chi(O_D) - D^3`.
/// `points` must be the CYTools-ordered divisor points used by the
/// intersection tensor, and `kklt_basis` indexes into that point list.
pub fn compute_kklt_divisor_chi(
    polytope: &Polytope,
    points: &[Point],
    kappa_all: &Intersection,
    kklt_basis: &[usize],
) -> Result<Vec<I64<Finite>>> {
    if polytope.dim() != 4 {
        return Err(Error::InvalidInput(
            "divisor chi is currently implemented for 4D CY hypersurface polytopes".into(),
        ));
    }
    if points.len() != kappa_all.dim() {
        return Err(Error::InvalidInput(format!(
            "point/intersection dimension mismatch: points={}, kappa_dim={}",
            points.len(),
            kappa_all.dim()
        )));
    }

    let (hull, dual_vertices, dual_points) = face_pairing_context(polytope, points)?;

    let mut out = Vec::with_capacity(kklt_basis.len());
    for &point_idx in kklt_basis {
        let Some(point) = points.get(point_idx) else {
            return Err(Error::InvalidInput(format!(
                "KKLT basis index {point_idx} out of bounds for {} points",
                points.len()
            )));
        };
        if point.coords().iter().all(|&coord| coord == 0) {
            return Err(Error::InvalidInput(
                "origin does not define a prime toric divisor".into(),
            ));
        }

        let sat = saturated_facets(point, &dual_vertices);
        let chi_o = match sat.len() {
            0 => {
                return Err(Error::InvalidInput(
                    "interior point does not define a prime toric divisor".into(),
                ));
            }
            1 => 1,
            2 => {
                // Non-favorable case: the divisor of a point strictly interior
                // to a 2-face splits into 1 + l*(dual face) disjoint rigid
                // components on the CY; chi(O) adds 1 per component.
                let split = dual_face_interior_lattice_points(
                    points,
                    &hull.vertex_indices,
                    &dual_points,
                    &dual_vertices,
                    &sat,
                )?;
                1_i64 + i64::try_from(split).expect("dual-face interior count fits i64")
            }
            3 => {
                let g = dual_face_interior_lattice_points(
                    points,
                    &hull.vertex_indices,
                    &dual_points,
                    &dual_vertices,
                    &sat,
                )?;
                1_i64 - i64::try_from(g).expect("dual-face interior count fits i64")
            }
            _ => {
                let g = dual_face_interior_lattice_points(
                    points,
                    &hull.vertex_indices,
                    &dual_points,
                    &dual_vertices,
                    &sat,
                )?;
                1_i64 + i64::try_from(g).expect("dual-face interior count fits i64")
            }
        };

        let d_cubed =
            integral_rational_to_i64(&kappa_all.get(point_idx, point_idx, point_idx), "D^3")?;
        out.push(I64::<Finite>::new(12 * chi_o - d_cubed));
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::f64_pos;
    use crate::types::range::CheckedRange;
    use crate::types::rational::Rational as TypedRational;
    use crate::types::tags::Pos;
    use malachite::Rational;

    fn finite(v: f64) -> F64<Finite> {
        F64::<Finite>::new(v).unwrap()
    }

    fn hypercube_4d_vertices() -> Vec<Point> {
        let mut points = Vec::new();
        for x0 in [-1, 1] {
            for x1 in [-1, 1] {
                for x2 in [-1, 1] {
                    for x3 in [-1, 1] {
                        points.push(Point::new(vec![x0, x1, x2, x3]));
                    }
                }
            }
        }
        points
    }

    #[test]
    fn test_divisor_volume_simple() {
        // κ_000 = 6, t = [2]
        // τ_0 = (1/2) × 6 × 2² = 12
        let mut kappa = Intersection::new(1);
        kappa.set(
            0,
            0,
            0,
            TypedRational::<Pos>::new(Rational::from(6)).unwrap(),
        );

        let t = vec![finite(2.0)];
        let tau = compute_divisor_volumes(&kappa, &t);

        assert_eq!(tau.len(), 1);
        assert!((tau[0].get() - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_divisor_volume_mixed() {
        // 2D case: κ_001 = 3, t = [1, 2]
        // τ_0 = (1/2) × 3 × 2 × 1 × 2 = 6 (κ_{001}t^0t^1 + κ_{010}t^1t^0)
        // τ_1 = (1/2) × 3 × 1² = 1.5 (κ_{100}t^0t^0)
        let mut kappa = Intersection::new(2);
        kappa.set(
            0,
            0,
            1,
            TypedRational::<Pos>::new(Rational::from(3)).unwrap(),
        );

        let t = vec![finite(1.0), finite(2.0)];
        let tau = compute_divisor_volumes(&kappa, &t);

        assert_eq!(tau.len(), 2);
        assert!((tau[0].get() - 6.0).abs() < 1e-10);
        assert!((tau[1].get() - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_jacobian_simple() {
        // κ_000 = 6, t = [2]
        // τ_0 = (1/2) × 6 × t² = 3t²
        // J_{00} = dτ_0/dt^0 = 6t = 12
        let mut kappa = Intersection::new(1);
        kappa.set(
            0,
            0,
            0,
            TypedRational::<Pos>::new(Rational::from(6)).unwrap(),
        );

        let t = vec![finite(2.0)];
        let jac = compute_divisor_jacobian(&kappa, &t);

        assert_eq!(jac.len(), 1);
        assert!((jac[0][0].get() - 12.0).abs() < 1e-10);
    }

    #[test]
    fn test_compute_kklt_divisor_chi_vertex_with_no_dual_interior_points() {
        let points = hypercube_4d_vertices();
        let polytope = Polytope::from_vertices(points.clone()).unwrap();
        let mut kappa_all = Intersection::new(points.len());
        kappa_all.set(
            0,
            0,
            0,
            TypedRational::<Finite>::from_raw(Rational::from(5)),
        );

        let chi = compute_kklt_divisor_chi(&polytope, &points, &kappa_all, &[0]).unwrap();

        assert_eq!(chi[0].get(), 7);
    }

    #[test]
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)] // Test code with small indices
    fn test_jacobian_consistency() {
        // Verify J numerically by finite differences
        let mut kappa = Intersection::new(2);
        kappa.set(
            0,
            0,
            1,
            TypedRational::<Pos>::new(Rational::from(3)).unwrap(),
        );
        kappa.set(
            1,
            1,
            1,
            TypedRational::<Pos>::new(Rational::from(2)).unwrap(),
        );

        let t = vec![finite(1.5), finite(2.5)];
        let jac = compute_divisor_jacobian(&kappa, &t);

        // Finite difference check
        let eps = f64_pos!(1e-6);
        let two = f64_pos!(2.0);
        let indices = CheckedRange::new(0, 2);
        for m in indices.iter_non_neg() {
            let m_idx = m.get() as usize;
            let mut t_plus = t.clone();
            let mut t_minus = t.clone();
            // Finite + Pos = Finite (algebra handles this)
            t_plus[m_idx] = t_plus[m_idx] + eps;
            // Finite - Pos = Finite (algebra handles this)
            t_minus[m_idx] = t_minus[m_idx] - eps;

            let tau_plus = compute_divisor_volumes(&kappa, &t_plus);
            let tau_minus = compute_divisor_volumes(&kappa, &t_minus);

            for i in indices.iter_non_neg() {
                let i_idx = i.get() as usize;
                let j_numerical = (tau_plus[i_idx] - tau_minus[i_idx]) / (two * eps);
                // Compare using raw values since abs() returns NonNeg
                assert!(
                    (jac[i_idx][m_idx] - j_numerical).abs().get() < 1e-5,
                    "J[{}][{}] mismatch: {} vs {}",
                    i_idx,
                    m_idx,
                    jac[i_idx][m_idx].get(),
                    j_numerical.get()
                );
            }
        }
    }
}
