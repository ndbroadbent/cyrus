//! Coverage for selected primal small curves whose GV invariants the toric
//! formulas miss.
//!
//! Exact no-decomposition certificates (minimal grading degree, extremal wall
//! rays) computed via the one-dimensional HKTY series, plus the explicitly
//! opted-in bounded-impact deferral.

use std::collections::{HashMap, HashSet};

use crate::compute_grading_vector;
use crate::types::f64::F64;
use crate::types::tags::Finite;

use super::gv_basis::vector_gv_basis_data;
use super::{PrimalGeom, PrimalIntersection};

/// Exact GV for a missing class that is a multiple of an extremal wall Mori
/// generator, via the one-dimensional HKTY series.
///
/// Compute GV invariants for missing classes whose grading degree is below
/// twice the minimal positive generator degree. Such classes admit no
/// decomposition into positive-degree effective classes at all (any
/// decomposition needs two parts of at least the minimal degree), so the
/// HKTY inversion for them involves no lower instanton channels and the
/// one-dimensional series on the class itself is exact regardless of the
/// decomposition domain — the cap-vs-effective-cone caveat that blocks
/// higher-degree exotic origin circuits cannot apply. (Verified against the
/// published GV checkpoint and a CYTools face computation on 5-113-4627's
/// degree-2 origin circuit: GV = -2.)
/// If the missing class is a positive integer multiple of a wall Mori
/// generator that an exact integer separator certifies as EXTREMAL among the
/// wall generators, its minimal face in the toric Mori cone is that single
/// ray: every decomposition is a lower multiple of the same primitive, so
/// the one-dimensional HKTY series is the exact computation. (Verified on
/// 7-51-13590, where all three remaining origin circuits are extremal rays
/// with GV = 1, matching the published checkpoint and CYTools.)
fn try_extremal_wall_ray_gv(
    wall_generators: &[Vec<i64>],
    intersection: &PrimalIntersection,
    q_matrix: &[Vec<i64>],
    grading: &[i64],
    class: &[i64],
) -> Result<Option<malachite::Integer>, String> {
    let Some((generator_index, multiple)) = parallel_ambient_ray_multiple(class, wall_generators)
    else {
        return Ok(None);
    };
    let generators_basis: Vec<Vec<i64>> = wall_generators
        .iter()
        .map(|generator| {
            intersection
                .basis
                .iter()
                .map(|&idx| generator.get(idx).copied().unwrap_or(0))
                .collect()
        })
        .collect();
    let target_basis = &generators_basis[generator_index];
    let lp_options = crate::SupportingMoriFaceLpSearchOptions::default();
    let certificate = crate::gv::diagnose_extremal_mori_ray_separator_by_lp_search(
        target_basis,
        &generators_basis,
        &lp_options,
    )
    .map_err(|e| format!("extremal wall-ray separator search failed: {e}"))?;
    if certificate.certificate.is_none() {
        eprintln!(
            "[INFO] extremal wall-ray GV: class is a wall-generator multiple but no exact extremality certificate found (status={}); leaving for other methods",
            certificate.status
        );
        return Ok(None);
    }
    let series = crate::compute_ambient_one_dimensional_ray_gv_series(
        &wall_generators[generator_index],
        &intersection.basis,
        grading,
        q_matrix,
        &intersection.kappa_basis,
        multiple,
    )
    .map_err(|e| format!("one-dimensional wall-ray GV series failed: {e}"))?;
    series
        .values
        .get(multiple as usize - 1)
        .cloned()
        .map(Some)
        .ok_or_else(|| {
            format!(
                "one-dimensional wall-ray GV series returned {} values; multiple {multiple} requested",
                series.values.len()
            )
        })
}

/// Find a generator that the ambient class is a positive integer multiple of
/// (componentwise `class == k * generator`).
fn parallel_ambient_ray_multiple(class: &[i64], rays: &[Vec<i64>]) -> Option<(usize, u32)> {
    'next: for (index, ray) in rays.iter().enumerate() {
        if ray.len() != class.len() {
            continue;
        }
        let mut multiple: Option<i64> = None;
        for (&c, &r) in class.iter().zip(ray.iter()) {
            if c == 0 && r == 0 {
                continue;
            }
            if c == 0 || r == 0 || c % r != 0 {
                continue 'next;
            }
            let ratio = c / r;
            if ratio <= 0 {
                continue 'next;
            }
            match multiple {
                None => multiple = Some(ratio),
                Some(existing) if existing == ratio => {}
                Some(_) => continue 'next,
            }
        }
        if let Some(k) = multiple
            && let Ok(k) = u32::try_from(k)
        {
            return Some((index, k));
        }
    }
    None
}

/// Split off missing-GV classes whose bounded impact is below the requested
/// `--max-missing-gv-impact` threshold.
///
/// A class whose maximal possible contribution (assuming
/// `|GV| <= missing_gv_abs_bound`) stays under the threshold is excluded from
/// the targets and V_string with a loud warning, and re-verified at the
/// solved Kahler point by `verify_deferred_missing_gv_bounds`.
pub fn defer_bounded_impact_missing_gvs(
    missing_gv_classes: &mut Vec<Vec<i64>>,
    intersection: &PrimalIntersection,
    kklt_basis: &[usize],
    selection_t: &[F64<Finite>],
    max_missing_gv_impact: f64,
    missing_gv_abs_bound: u32,
) -> Vec<Vec<i64>> {
    let mut deferred = Vec::new();
    if missing_gv_classes.is_empty() || max_missing_gv_impact <= 0.0 {
        return deferred;
    }
    let mut still_missing = Vec::new();
    for class in missing_gv_classes.drain(..) {
        let impact =
            missing_gv_unit_impact_bound(&class, &intersection.basis, kklt_basis, selection_t)
                .map(|unit| unit * f64::from(missing_gv_abs_bound));
        match impact {
            Some(bound) if bound <= max_missing_gv_impact => {
                eprintln!(
                    "[WARN] deferring missing GV curve with bounded impact {bound:.3e} <= --max-missing-gv-impact {max_missing_gv_impact:.3e} (assuming |GV| <= {missing_gv_abs_bound}); its contribution is EXCLUDED from targets and V_string; class={class:?}"
                );
                deferred.push(class);
            }
            _ => still_missing.push(class),
        }
    }
    *missing_gv_classes = still_missing;
    deferred
}

/// Re-verify the deferral assumption at the solved Kahler point: curve
/// volumes shrink during the corrected solve, so a bound granted at
/// selection time must be re-established where it actually matters.
pub fn verify_deferred_missing_gv_bounds(
    deferred_missing_classes: &[Vec<i64>],
    intersection: &PrimalIntersection,
    kklt_basis: &[usize],
    solved_t: &[F64<Finite>],
    max_missing_gv_impact: f64,
    missing_gv_abs_bound: u32,
) -> Result<(), String> {
    for class in deferred_missing_classes {
        let realized =
            missing_gv_unit_impact_bound(class, &intersection.basis, kklt_basis, solved_t)
                .map(|unit| unit * f64::from(missing_gv_abs_bound));
        match realized {
            Some(bound) if bound <= max_missing_gv_impact => {}
            other => {
                return Err(format!(
                    "deferred missing-GV curve exceeds the impact bound at the solved Kahler point (realized={other:?}, limit={max_missing_gv_impact:.3e}); class={class:?}"
                ));
            }
        }
    }
    Ok(())
}

/// Conservative per-unit-GV bound on a missing curve's contribution to any
/// KKLT tau target or to V_string at the Kahler point `t` (computed basis).
/// Uses the positive-parity dilogarithm (which dominates the odd-parity one)
/// and Li3(x) <= Li2(x) on (0,1). Returns None when the curve volume is not
/// positive (no bound possible without continuation data).
fn missing_gv_unit_impact_bound(
    class: &[i64],
    basis: &[usize],
    kklt_basis: &[usize],
    t: &[F64<Finite>],
) -> Option<f64> {
    let volume: f64 = basis
        .iter()
        .zip(t.iter())
        .map(|(&idx, ti)| class.get(idx).copied().unwrap_or(0) as f64 * ti.get())
        .sum();
    if volume.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
        return None;
    }
    let dilog = crate::kklt::gv_dilog_from_curve_volume_checked(volume, 0).ok()?;
    let max_kklt_coeff = kklt_basis
        .iter()
        .map(|&idx| class.get(idx).copied().unwrap_or(0).unsigned_abs())
        .max()? as f64;
    let two_pi = std::f64::consts::TAU;
    let tau_bound = dilog * max_kklt_coeff / (2.0 * two_pi.powi(2));
    let volume_bound = (1.0 + two_pi * volume) * dilog / (2.0 * two_pi.powi(3));
    Some(tau_bound.max(volume_bound))
}

/// Compute exact GV invariants for missing ambient classes that provably do
/// not decompose.
///
/// These are the classes below twice the minimal generator grading degree, or
/// positive multiples of an extremal wall Mori generator. For each such class
/// the one-dimensional HKTY series is the exact computation; classes that may
/// decompose are left for other methods.
/// Basis-fixed precompute for the minimal-degree GV pre-pass.
///
/// These are the parts of the minimal-degree computation that depend only on
/// the geometry and the production divisor basis — NOT the Kähler point or the
/// set of missing classes. Hoisting them out lets an admissible-basis scan pay
/// the 561k-ray projection, grading, and Mori-generator costs once instead of
/// on every candidate probe.
pub struct MinimalDegreeGvContext {
    /// GV grading vector in the production divisor basis.
    grading: Vec<i64>,
    /// Curve-basis (charge) matrix for the one-dimensional HKTY series.
    q_matrix: Vec<Vec<i64>>,
    /// Minimal positive grading degree among the Mori cap generators.
    min_generator_degree: i128,
    /// Mori cap generators at exactly `min_generator_degree`, in basis coords.
    minimal_degree_generators: HashSet<Vec<i64>>,
    /// Wall-based (ridge-adjacent) Mori generators in ambient coordinates.
    wall_generators: Vec<Vec<i64>>,
}

/// Build the basis-fixed [`MinimalDegreeGvContext`] once for a geometry and its
/// production divisor basis.
///
/// This carries the dominant per-probe cost (the 561k-ray projection to the
/// divisor basis, the grading vector, and the wall Mori generators), so an
/// admissible-basis scan pays it a single time.
///
/// # Errors
/// Fails if the GV basis data, the grading vector, or the wall Mori generators
/// cannot be built, or if there are no positive-degree Mori generators.
pub fn build_minimal_degree_gv_context(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    ambient_rays: &[Vec<i64>],
) -> Result<MinimalDegreeGvContext, String> {
    let gv_basis_data = vector_gv_basis_data(
        ambient_rays,
        &intersection.linrels,
        &intersection.basis,
        "primal",
    )?;
    let rays = &gv_basis_data.mori_rays;
    let grading = compute_grading_vector(rays)
        .ok_or_else(|| "failed to compute primal GV grading vector".to_string())?;
    let degree_of = |basis_class: &[i64]| -> i128 {
        basis_class
            .iter()
            .zip(grading.iter())
            .map(|(&c, &g)| i128::from(c) * i128::from(g))
            .sum()
    };
    let min_generator_degree = rays
        .iter()
        .map(|ray| degree_of(ray))
        .filter(|&deg| deg > 0)
        .min()
        .ok_or_else(|| "no positive-degree Mori generators".to_string())?;
    let minimal_degree_generators: HashSet<Vec<i64>> = rays
        .iter()
        .filter(|ray| degree_of(ray) == min_generator_degree)
        .cloned()
        .collect();

    // Wall-based Mori generators (ridge-adjacent simplex relations only).
    // Unlike the all-pairs cap list above, this is the honest minimal
    // generating description of the toric Mori cone, so exact extremality
    // certificates are meaningful against it.
    let wall_generators: Vec<Vec<i64>> =
        crate::compute_mori_generators(&geom.triangulation, &geom.triangulation_points)
            .map_err(|e| format!("failed to compute wall Mori generators: {e}"))?
            .generators()
            .iter()
            .map(|generator| generator.iter().map(|value| value.get()).collect())
            .collect();
    Ok(MinimalDegreeGvContext {
        grading,
        q_matrix: gv_basis_data.q_matrix,
        min_generator_degree,
        minimal_degree_generators,
        wall_generators,
    })
}

/// Compute minimal-degree GV invariants for `required_ambient_classes` using a
/// precomputed [`MinimalDegreeGvContext`].
///
/// This is the per-Kähler-point work: the missing-class loop only, with all
/// basis-fixed data supplied by `ctx`.
///
/// # Errors
/// Fails if a one-dimensional HKTY series for a non-decomposable missing class
/// cannot be computed.
pub fn compute_minimal_degree_gv_with_context(
    ctx: &MinimalDegreeGvContext,
    intersection: &PrimalIntersection,
    required_ambient_classes: &[Vec<i64>],
) -> Result<HashMap<Vec<i64>, malachite::Integer>, String> {
    let grading = &ctx.grading;
    let min_generator_degree = ctx.min_generator_degree;
    let degree_of = |basis_class: &[i64]| -> i128 {
        basis_class
            .iter()
            .zip(grading.iter())
            .map(|(&c, &g)| i128::from(c) * i128::from(g))
            .sum()
    };
    let mut out = HashMap::new();
    for class in required_ambient_classes {
        let basis_class: Vec<i64> = intersection
            .basis
            .iter()
            .map(|&idx| class.get(idx).copied().unwrap_or(0))
            .collect();
        let class_degree = degree_of(&basis_class);
        // A class can decompose only as a sum of two semigroup elements, and
        // any semigroup element of the minimal generator degree must be a
        // generator itself. So for class degree < 2*min nothing decomposes,
        // and for degree == 2*min the only candidates are exact generator
        // pairs, checked here over integer vectors.
        let decomposable = if class_degree <= 0 {
            true
        } else if class_degree < 2 * min_generator_degree {
            false
        } else if class_degree == 2 * min_generator_degree {
            ctx.minimal_degree_generators.iter().any(|generator| {
                let remainder: Vec<i64> = basis_class
                    .iter()
                    .zip(generator.iter())
                    .map(|(&c, &g)| c - g)
                    .collect();
                ctx.minimal_degree_generators.contains(&remainder)
            })
        } else {
            true
        };
        if decomposable {
            if let Some(value) = try_extremal_wall_ray_gv(
                &ctx.wall_generators,
                intersection,
                &ctx.q_matrix,
                grading,
                class,
            )? {
                eprintln!(
                    "[INFO] extremal wall-ray GV: missing class at grading degree {class_degree} is a multiple of an extremal wall Mori generator (exact separator certificate); one-dimensional HKTY series gives GV={value}"
                );
                out.insert(class.clone(), value);
                continue;
            }
            eprintln!(
                "[INFO] minimal-degree GV: missing class at grading degree {class_degree} (min generator degree {min_generator_degree}) may decompose; leaving for other methods"
            );
            continue;
        }
        let gcd = basis_class
            .iter()
            .fold(0u64, |acc, &c| gcd_u64(acc, c.unsigned_abs()));
        if gcd != 1 {
            // k-fold multiples decompose into k copies of the primitive and
            // always fail the degree gate; reaching here would be a bug.
            continue;
        }
        let series = crate::compute_ambient_one_dimensional_ray_gv_series(
            class,
            &intersection.basis,
            grading,
            &ctx.q_matrix,
            &intersection.kappa_basis,
            1,
        )
        .map_err(|e| format!("one-dimensional minimal-degree GV series failed: {e}"))?;
        let Some(value) = series.values.first() else {
            return Err("one-dimensional minimal-degree GV series returned no values".to_string());
        };
        eprintln!(
            "[INFO] minimal-degree GV: missing class at grading degree {class_degree} (< 2 x min generator degree {min_generator_degree}) admits no decompositions; one-dimensional HKTY series gives GV={value}"
        );
        out.insert(class.clone(), value.clone());
    }
    Ok(out)
}

/// Compute minimal-degree GV invariants, building the basis-fixed context
/// inline.
///
/// Equivalent to [`build_minimal_degree_gv_context`] followed by
/// [`compute_minimal_degree_gv_with_context`]; prefer the split form when
/// covering many Kähler points over a fixed production basis (the context build
/// dominates the cost).
///
/// # Errors
/// Fails if the context build or the per-class computation fails.
pub fn compute_minimal_degree_gv_by_ambient_class(
    geom: &PrimalGeom,
    intersection: &PrimalIntersection,
    required_ambient_classes: &[Vec<i64>],
    ambient_rays: &[Vec<i64>],
) -> Result<HashMap<Vec<i64>, malachite::Integer>, String> {
    let ctx = build_minimal_degree_gv_context(geom, intersection, ambient_rays)?;
    compute_minimal_degree_gv_with_context(&ctx, intersection, required_ambient_classes)
}

fn gcd_u64(a: u64, b: u64) -> u64 {
    if b == 0 { a } else { gcd_u64(b, a % b) }
}
