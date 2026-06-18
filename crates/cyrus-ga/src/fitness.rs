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
use cyrus_core::pipeline::{EvaluationRequest, EvaluationResult, evaluate_vacuum};
use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::Pos;
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
    /// Run the in-process deep KKLT verification (the full primal-side
    /// arXiv:2107.09064 SS6 stabilization, `GaGeometry::deep_verify`) on a
    /// VALID candidate whose proxy CPL fit `(w0, wa)` lands within this many
    /// DESI sigma on BOTH parameters (see [`within_desi_sigma`]). `None`
    /// disables it. The gate is PHYSICAL, not a fitness number: the primal
    /// solve is expensive, so it runs only on candidates whose dark-energy
    /// equation of state is actually DESI-relevant (the per-polytope coverage
    /// cache then bounds the cost to one admissible-basis scan per polytope). A
    /// passing deep-verify lifts the candidate into a strictly higher fitness
    /// band (see [`evaluate_fitness`]).
    #[serde(default)]
    pub deep_verify_desi_sigma: Option<f64>,
}

/// Fitness band floor for a candidate whose primal KKLT vacuum was solved
/// (deep-verify produced a real `V0`) but whose worldsheet-instanton expansion
/// left curves uncovered/deferred. Strictly above any non-deep-verified valid.
const STABILIZED_FLOOR: f64 = 1000.0;
/// Fitness band floor for a fully GV-controlled deep-verified vacuum (no
/// deferred curves): the instanton expansion is under control. Strictly above
/// the stabilized band.
const CONTROLLED_FLOOR: f64 = 2000.0;

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
            deep_verify_desi_sigma: None,
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
    /// Deep-verify outcome when it was triggered: `"controlled"`,
    /// `"kklt-stabilized"`, or a `"deep-verify failed: ..."` reason. `None`
    /// when the candidate never crossed the deep-verify threshold.
    #[serde(default)]
    pub deep_tier: Option<String>,
    /// Real `log10|V0|` from the primal KKLT solve, when deep-verify succeeded
    /// (replaces the mirror proxy for the verified candidate).
    #[serde(default)]
    pub deep_log10_v0: Option<f64>,
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

/// Strictly MONOTONIC failure ladder: each pipeline gate a candidate
/// clears moves it into a higher, non-overlapping band, and within each
/// band a continuous gradient rewards getting closer to the next gate. The
/// pipeline order is tadpole -> N invertible -> orthogonality (K.p = 0) ->
/// mirror cone interior (q.p > 0) -> racetrack -> valid, so the bands must
/// ascend in that exact order.
///
/// The previous ladder was NON-monotonic: a candidate that achieved exact
/// orthogonality but landed on a cone wall scored -1500, BELOW a merely
/// near-orthogonal one at -800, so the GA fled the exact-isotropic region
/// where every real vacuum lives. Now clearing orthogonality always
/// dominates failing it.
fn invalid_fitness(reason: &str, result: &EvaluationResult, cfg: &FitnessConfig) -> f64 {
    // Band floors, ascending with pipeline progress. Width 400 each, so a
    // fully-graded lower band never reaches the next band's floor.
    const TADPOLE: f64 = -3000.0;
    const W0_DEAD_END: f64 = -2800.0; // PFV but W0 identically 0: reward-hack
    const N_SINGULAR: f64 = -2600.0;
    const ORTHOGONALITY: f64 = -2400.0; // + up to 380 as |K.p| -> 0
    const CONE_WALL: f64 = -1900.0; // passed orthogonality; + up to 380
    const NO_RACETRACK: f64 = -1400.0; // genuine PFV; + up to 380
    const RACETRACK_UNCTRL: f64 = -800.0; // racetrack found, drift too big
    const OTHER: f64 = -2000.0;

    if reason.starts_with("Tadpole") {
        let overshoot = (result.q_flux - cfg.q_max).max(-result.q_flux).max(0.0);
        TADPOLE - overshoot
    } else if reason.contains("cancelled exactly") {
        W0_DEAD_END
    } else if reason.starts_with("N matrix") {
        N_SINGULAR
    } else if reason.starts_with("Orthogonality") {
        // |K.p| -> 0 climbs the band toward the cone gate.
        let violation = result.k_dot_p.map_or(f64::INFINITY, f64::abs);
        380.0f64.mul_add((-violation / 2.0).exp(), ORTHOGONALITY)
    } else if reason.starts_with("Flat direction") {
        // EXACT orthogonality achieved; p on/outside a cone wall. The
        // closer the worst curve's q.p is to zero from below, the closer
        // p is to the cone interior - climb toward it.
        let margin = result.cone_margin.unwrap_or(f64::NEG_INFINITY);
        // margin in (-inf, eps); map to [0, 1) increasing as margin -> 0.
        let closeness = (margin.min(0.0) / 1.0).exp(); // exp of a non-positive number
        380.0f64.mul_add(closeness, CONE_WALL)
    } else if reason.starts_with("No stable racetrack") {
        // A genuine perturbatively flat vacuum (all three PFV conditions
        // hold) that lacks an instanton hierarchy for a racetrack. The
        // closest-to-valid failure: top of the invalid ladder, just below
        // a found-but-uncontrolled racetrack.
        NO_RACETRACK + 380.0
    } else if reason.starts_with("Racetrack refinement") {
        // Two-term racetrack solved, full-series refinement out of control.
        RACETRACK_UNCTRL
    } else {
        OTHER
    }
}

/// Weighted DESI-targeted component score for an energy scale `abs_v0` and
/// string coupling `g_s`: a height term (closeness of `log10|V0|` to the
/// observed dark-energy scale) gating a DESI-slope term and a weak-coupling
/// term. Returns the combined score (in roughly `[0, 100*(weights)]`) and the
/// CPL `(w0, wa)` fit when the quintessence integration ran. Shared by the
/// valid (proxy `V0`) and the deep-verified (real `V0`) scoring so both bands
/// reward the same physics.
fn weighted_component_score(
    abs_v0: f64,
    g_s: f64,
    cfg: &FitnessConfig,
) -> (f64, Option<(f64, f64)>) {
    let height_score = 100.0 * (-(abs_v0.log10() - cfg.target_log10_v0).abs() / 25.0).exp();
    let cpl = quintessence_cpl(abs_v0, cfg);
    let slope_norm = cpl.map_or(0.0, |(w0_fit, wa_fit)| {
        gaussian_score(w0_fit, cfg.desi_w0, cfg.desi_w0_sigma)
            * gaussian_score(wa_fit, cfg.desi_wa, cfg.desi_wa_sigma)
    });
    let gs_norm = (1.0 - g_s).clamp(0.0, 1.0);
    let combined = height_score
        * cfg.weight_slope.mul_add(
            slope_norm,
            cfg.weight_gs.mul_add(gs_norm, cfg.weight_height),
        );
    (combined, cpl)
}

/// Apply the deep-verify band to a VALID candidate that cleared the threshold:
/// run the full primal KKLT stabilization and, on success, lift the candidate
/// into the stabilized or controlled band scored on the REAL `V0`. A failure
/// leaves it in the valid band (it still passed the scan) with the reason
/// recorded - the ladder never demotes a valid candidate for being tested.
fn apply_deep_verify(
    geom: &GaGeometry,
    cfg: &FitnessConfig,
    ek0: F64<Pos>,
    g_s: F64<Pos>,
    w0: F64<Pos>,
    report: &mut FitnessReport,
) {
    match geom.deep_verify(ek0, g_s, w0) {
        Ok(verdict) => {
            let real_abs_v0 = verdict.v0.get().abs();
            report.deep_log10_v0 = Some(real_abs_v0.log10());
            let (combined, cpl) = weighted_component_score(real_abs_v0, g_s.get(), cfg);
            if let Some((w0_fit, wa_fit)) = cpl {
                report.cpl_w0 = Some(w0_fit);
                report.cpl_wa = Some(wa_fit);
            }
            let controlled = verdict.gv_controlled && verdict.deferred_missing_gv_count == 0;
            let (floor, tier) = if controlled {
                (CONTROLLED_FLOOR, "controlled")
            } else {
                (STABILIZED_FLOOR, "kklt-stabilized")
            };
            report.deep_tier = Some(tier.to_string());
            report.tier = tier.to_string();
            report.fitness = floor + combined;
        }
        Err(e) => {
            report.deep_tier = Some(format!("deep-verify failed: {e}"));
        }
    }
}

/// Whether a candidate's CPL fit `(w0, wa)` lands within `n_sigma` of the DESI
/// w0-wa constraints on BOTH parameters - the physical gate for deep
/// verification (only a DESI-relevant dark-energy equation of state is worth the
/// primal KKLT solve). A `None` fit (the quintessence integration did not run)
/// never qualifies.
fn within_desi_sigma(
    cpl_w0: Option<f64>,
    cpl_wa: Option<f64>,
    n_sigma: f64,
    cfg: &FitnessConfig,
) -> bool {
    let (Some(w0), Some(wa)) = (cpl_w0, cpl_wa) else {
        return false;
    };
    (w0 - cfg.desi_w0).abs() <= n_sigma * cfg.desi_w0_sigma
        && (wa - cfg.desi_wa).abs() <= n_sigma * cfg.desi_wa_sigma
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
                deep_tier: None,
                deep_log10_v0: None,
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
        deep_tier: None,
        deep_log10_v0: None,
    };

    if !result.success {
        let reason = result
            .reason
            .clone()
            .unwrap_or_else(|| "unknown".to_string());
        report.fitness = invalid_fitness(&reason, &result, cfg);
        report.tier = reason;
        return report;
    }

    let racetrack = result.racetrack.expect("success implies racetrack");
    let vacuum = result.vacuum.expect("success implies vacuum");
    let g_s = racetrack.g_s.get();
    let abs_v0 = vacuum.v0.get().abs();
    report.g_s = Some(g_s);
    report.w0 = Some(vacuum.w0.get());
    report.log10_abs_v0 = Some(abs_v0.log10());

    // Valid band: DESI-targeted weighted component score on the PROXY V0.
    // Height GATES the score (DESI is a target SCALE first, log|V0| at the
    // observed dark-energy density, and a slope second; weak coupling is a
    // control tiebreaker). Multiplicative gating means slope and g_s only
    // reward candidates already near the right scale - the old additive form
    // let a weakly-coupled vacuum 95 orders off-scale top the leaderboard.
    let (combined, cpl) = weighted_component_score(abs_v0, g_s, cfg);
    if let Some((w0_fit, wa_fit)) = cpl {
        report.cpl_w0 = Some(w0_fit);
        report.cpl_wa = Some(wa_fit);
    }
    report.fitness = combined;

    // Deep-verify band: the toil this automates. The gate is PHYSICAL - a valid
    // candidate is worth the expensive primal KKLT solve only if its proxy
    // dark-energy equation of state (w0, wa) already sits within the configured
    // DESI sigma; a pass lifts it into a strictly higher band scored on the real
    // V0. The per-polytope coverage cache keeps the cost bounded.
    if let Some(n_sigma) = cfg.deep_verify_desi_sigma
        && within_desi_sigma(report.cpl_w0, report.cpl_wa, n_sigma, cfg)
    {
        apply_deep_verify(geom, cfg, vacuum.ek0, vacuum.g_s, vacuum.w0, &mut report);
    }
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
