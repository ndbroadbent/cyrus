//! Fitness evaluation: tiered validity gates plus DESI-targeted scoring.
//!
//! The hard physics comes from `cyrus_core::evaluate_vacuum` (the validated
//! Demirtas-McAllister scan chain). Failures map to graded tier penalties so
//! the GA can climb toward validity; successes are scored on three
//! components:
//!
//! 1. **Height**: `log10|V0|` against a target (default the observed
//!    dark-energy scale, `log10(2.888e-122) = -121.54`).
//! 2. **Slope (DESI)**: the candidate's energy scale is dressed as a
//!    thawing-axion quintessence potential `V = |V0| (1 - cos(phi/f))`,
//!    integrated through the real Friedmann + Klein-Gordon equations
//!    (`cyrus_core::cosmology`), CPL-fitted, and scored against the DESI
//!    (w0, wa) measurement. Only candidates near the right height CAN
//!    produce the measured slope - the components are physically coupled.
//! 3. **Weak coupling**: small g_s (perturbative control).
//!
//! v1 modelling approximations (documented, revisit with cyrus-cosmology
//! integration): single-field axion, decay constant `f` is a configured
//! knob (computing it from the Kahler metric is follow-up work), and the
//! axion potential height is identified with the candidate's |V0| scale.

use cyrus_core::cosmology::{CosmologyParams, Potential, solve_cosmology};
use cyrus_core::pipeline::{EvaluationRequest, evaluate_vacuum};
use serde::{Deserialize, Serialize};

use crate::genome::Genome;
use crate::geometry::GaGeometry;

/// H0 in Planck units (M_pl), used to express V0 in critical-density units.
pub const H0_PLANCK: f64 = 1.18e-61;

/// Configuration for fitness scoring.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FitnessConfig {
    /// D3 tadpole bound Q_max.
    pub q_max: f64,
    /// Target `log10|V0|` (observed dark energy: -121.54).
    pub target_log10_v0: f64,
    /// DESI CPL equation-of-state today.
    pub desi_w0: f64,
    /// DESI uncertainty on w0.
    pub desi_w0_sigma: f64,
    /// DESI CPL evolution parameter.
    pub desi_wa: f64,
    /// DESI uncertainty on wa.
    pub desi_wa_sigma: f64,
    /// Axion decay constant in M_pl (v1 modelling knob).
    pub decay_constant: f64,
    /// Initial axion displacement in units of f (thawing start).
    pub initial_displacement: f64,
    /// Component weights: height, slope, weak coupling.
    pub weight_height: f64,
    /// Weight of the DESI slope component.
    pub weight_slope: f64,
    /// Weight of the weak-coupling component.
    pub weight_gs: f64,
}

impl Default for FitnessConfig {
    fn default() -> Self {
        Self {
            q_max: 138.0,
            target_log10_v0: -121.54,
            desi_w0: -0.45,
            desi_w0_sigma: 0.21,
            desi_wa: -1.8,
            desi_wa_sigma: 0.6,
            decay_constant: 0.5,
            initial_displacement: 2.0,
            weight_height: 1.0,
            weight_slope: 1.0,
            weight_gs: 0.2,
        }
    }
}

/// Everything the GA records about one evaluated genome.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FitnessReport {
    /// Scalar fitness (higher is better).
    pub fitness: f64,
    /// Pipeline tier reached ("valid" or the failure stage).
    pub tier: String,
    /// String coupling, when the racetrack solved.
    pub g_s: Option<f64>,
    /// Flux superpotential |W0|.
    pub w0: Option<f64>,
    /// `log10|V0|` of the AdS vacuum.
    pub log10_abs_v0: Option<f64>,
    /// Flux tadpole.
    pub q_flux: f64,
    /// CPL fit of the integrated w(z), when computed.
    pub cpl_w0: Option<f64>,
    /// CPL evolution parameter of the integrated w(z).
    pub cpl_wa: Option<f64>,
}

struct AxionPotential {
    /// Height in critical-density units (3 H0^2 M_pl^2 = 3 in solver units).
    height: f64,
    decay_constant: f64,
}

impl Potential for AxionPotential {
    fn value(&self, phi: f64) -> f64 {
        self.height * (1.0 - (phi / self.decay_constant).cos())
    }
    fn deriv(&self, phi: f64) -> f64 {
        self.height / self.decay_constant * (phi / self.decay_constant).sin()
    }
}

/// Least-squares CPL fit `w(a) = w0 + wa (1 - a)` over `z <= 3`.
#[must_use]
pub fn fit_cpl(redshifts: &[f64], w_values: &[f64]) -> Option<(f64, f64)> {
    let points: Vec<(f64, f64)> = redshifts
        .iter()
        .zip(w_values.iter())
        .filter(|(z, w)| **z <= 3.0 && w.is_finite())
        .map(|(z, w)| (1.0 - 1.0 / (1.0 + z), *w))
        .collect();
    if points.len() < 4 {
        return None;
    }
    let n = points.len() as f64;
    let sum_x: f64 = points.iter().map(|(x, _)| x).sum();
    let sum_y: f64 = points.iter().map(|(_, y)| y).sum();
    let sum_xx: f64 = points.iter().map(|(x, _)| x * x).sum();
    let sum_xy: f64 = points.iter().map(|(x, y)| x * y).sum();
    let denom = n.mul_add(sum_xx, -(sum_x * sum_x));
    if denom.abs() < 1e-12 {
        return None;
    }
    let wa = n.mul_add(sum_xy, -(sum_x * sum_y)) / denom;
    let w0 = (sum_y - wa * sum_x) / n;
    Some((w0, wa))
}

/// Integrate the thawing-axion cosmology for a candidate's energy scale and
/// return the CPL fit of w(z).
#[must_use]
pub fn quintessence_cpl(abs_v0: f64, cfg: &FitnessConfig) -> Option<(f64, f64)> {
    // Express the axion height in solver units where H0 = 1, M_pl = 1:
    // critical density = 3, observed dark energy ~ 0.7 * 3.
    let height = abs_v0 / (H0_PLANCK * H0_PLANCK);
    if !height.is_finite() || height <= 0.0 || height > 1e6 {
        return None;
    }
    let params = CosmologyParams {
        omega_m0: 0.3,
        omega_de0: 0.7,
        h0: 1.0,
    };
    let potential = AxionPotential {
        height,
        decay_constant: cfg.decay_constant,
    };
    let phi_i = cfg.initial_displacement * cfg.decay_constant;
    let result = solve_cosmology(&params, &potential, phi_i, 0.0, 100.0).ok()?;
    fit_cpl(&result.redshifts, &result.w_z)
}

fn gaussian_score(value: f64, target: f64, sigma: f64) -> f64 {
    let pull = (value - target) / sigma;
    (-0.5 * pull * pull).exp()
}

/// Evaluate a genome: tiered penalties for invalid candidates, weighted
/// component score for valid ones.
#[must_use]
pub fn evaluate_fitness(geom: &GaGeometry, cfg: &FitnessConfig, genome: &Genome) -> FitnessReport {
    let req = EvaluationRequest {
        kappa: &geom.kappa_basis,
        // The toric cone is too strict for real CY flat directions; the
        // physical gates (e^K0, racetrack) do the filtering. See pipeline.rs.
        mori: None,
        gv: &geom.gv,
        h11: geom.mirror_h11,
        h21: geom.mirror_h21,
        // The binding tadpole is the geometry's own orientifold bound
        // (a global config can only tighten it further). Discovered the
        // hard way: a fixed 138 admitted candidates at 4x the true bound.
        q_max: cfg.q_max.min(geom.q_d3),
    };
    let result = match evaluate_vacuum(&req, &genome.k, &genome.m) {
        Ok(result) => result,
        Err(e) => {
            return FitnessReport {
                fitness: -3000.0,
                tier: format!("evaluation error: {e}"),
                g_s: None,
                w0: None,
                log10_abs_v0: None,
                q_flux: 0.0,
                cpl_w0: None,
                cpl_wa: None,
            };
        }
    };

    let mut report = FitnessReport {
        fitness: 0.0,
        tier: "valid".to_string(),
        g_s: None,
        w0: None,
        log10_abs_v0: None,
        q_flux: result.q_flux,
        cpl_w0: None,
        cpl_wa: None,
    };

    if !result.success {
        let reason = result.reason.unwrap_or_else(|| "unknown".to_string());
        // Graded tiers (validated prototype design): later stages score
        // higher so the GA can climb toward validity.
        report.fitness = if reason.starts_with("Tadpole") {
            // Overshoot above the bound or below zero, graded either way.
            let overshoot = (result.q_flux - cfg.q_max).max(-result.q_flux).max(0.0);
            -2000.0 - overshoot
        } else if reason.starts_with("N matrix") {
            -1800.0
        } else if reason.starts_with("Flat direction") {
            -1500.0
        } else if reason.starts_with("Orthogonality") {
            // Graded: |K.p| -> 0 climbs from -1200 toward the racetrack
            // tier, giving the GA a slope toward the Diophantine condition.
            let violation = result.k_dot_p.map_or(f64::INFINITY, f64::abs);
            400.0f64.mul_add((-violation / 2.0).exp(), -1200.0)
        } else if reason.starts_with("No stable racetrack") {
            -800.0
        } else if reason.contains("cancelled exactly") {
            // Symmetric fluxes with W0 identically zero are dead ends, not
            // progress; rank them below a missing racetrack so they cannot
            // colonize the frontier (observed reward-hacking attractor).
            -1000.0
        } else {
            -600.0
        };
        report.tier = reason;
        return report;
    }

    let racetrack = result.racetrack.expect("success implies racetrack");
    let vacuum = result.vacuum.expect("success implies vacuum");
    let g_s = racetrack.g_s.get();
    let abs_v0 = vacuum.v0.get().abs();
    let log10_abs_v0 = abs_v0.log10();
    report.g_s = Some(g_s);
    report.w0 = Some(vacuum.w0.get());
    report.log10_abs_v0 = Some(log10_abs_v0);

    // Height: each decade away from the target costs e^(-1/25).
    let height_score = 100.0 * (-(log10_abs_v0 - cfg.target_log10_v0).abs() / 25.0).exp();

    // Slope: only meaningful near the right height; the integration itself
    // enforces that (a wildly off-scale potential freezes at w = -1 or
    // never dominates).
    let cpl = quintessence_cpl(abs_v0, cfg);
    let slope_score = cpl.map_or(0.0, |(w0_fit, wa_fit)| {
        report.cpl_w0 = Some(w0_fit);
        report.cpl_wa = Some(wa_fit);
        100.0
            * gaussian_score(w0_fit, cfg.desi_w0, cfg.desi_w0_sigma)
            * gaussian_score(wa_fit, cfg.desi_wa, cfg.desi_wa_sigma)
    });

    // Weak coupling: full marks toward g_s -> 0, zero by g_s = 1.
    let gs_score = 100.0 * (1.0 - g_s).clamp(0.0, 1.0);

    report.fitness = cfg.weight_slope.mul_add(
        slope_score,
        cfg.weight_height
            .mul_add(height_score, cfg.weight_gs * gs_score),
    );
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cpl_fit_recovers_linear_w_of_a() {
        let true_w0: f64 = -0.45;
        let true_wa: f64 = -1.8;
        let redshifts: Vec<f64> = (0..40).map(|i| f64::from(i) * 0.075).collect();
        let w: Vec<f64> = redshifts
            .iter()
            .map(|&z| {
                let a = 1.0 / (1.0 + z);
                true_wa.mul_add(1.0 - a, true_w0)
            })
            .collect();
        let (w0, wa) = fit_cpl(&redshifts, &w).expect("fit");
        assert!((w0 - true_w0).abs() < 1e-9, "w0 = {w0}");
        assert!((wa - true_wa).abs() < 1e-9, "wa = {wa}");
    }

    #[test]
    fn quintessence_cpl_at_observed_scale_thaws() {
        let cfg = FitnessConfig::default();
        // |V0| at the observed dark-energy scale with the default knobs.
        let (w0, wa) = quintessence_cpl(2.888e-122, &cfg).expect("integration");
        // A thawing axion: w0 above -1, drifting; and far-off scales freeze.
        assert!(w0 > -1.0 && w0 < 0.0, "w0 = {w0}");
        assert!(wa.abs() < 5.0, "wa = {wa}");
        let frozen = quintessence_cpl(1e-200, &cfg);
        if let Some((w0_frozen, _)) = frozen {
            assert!(
                (w0_frozen - -1.0).abs() < 0.05 || w0_frozen.is_nan(),
                "deep-AdS-scale candidate should freeze at w = -1, got {w0_frozen}"
            );
        }
    }

    #[test]
    fn gaussian_score_peaks_at_target() {
        assert!((gaussian_score(-0.45, -0.45, 0.21) - 1.0).abs() < 1e-12);
        assert!(gaussian_score(-1.0, -0.45, 0.21) < 0.05);
    }
}
