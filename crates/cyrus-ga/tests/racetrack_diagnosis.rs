//! Diagnose a full-term racetrack refinement failure for a GA candidate:
//! does the Newton iteration diverge (numerics), or converge to a
//! stationary point that drifted away from the two-term seed (the
//! truncated racetrack is untrustworthy for this candidate)?
//!
//! Gated on CYRUS_RACETRACK_DIAG = "<pool.jsonl>|<name>|<k csv>|<m csv>".

use cyrus_core::f64_pos;
use cyrus_core::racetrack::{RacetrackTerm, build_racetrack_terms, solve_racetrack};
use cyrus_core::types::i64::I64;
use cyrus_core::types::tags::Finite;
use cyrus_ga::geometry::{DEFAULT_GV_MIN_POINTS, GaGeometry};
use cyrus_ga::multi::load_pool;

const TWO_PI: f64 = std::f64::consts::TAU;

/// Plain reimplementation of the refinement Newton loop with telemetry:
/// returns (iterations, final re, final im, converged, drift_fraction).
fn newton_with_telemetry(
    seed_re: f64,
    seed_im: f64,
    terms: &[RacetrackTerm],
) -> (usize, f64, f64, bool, f64) {
    let (mut re, mut im) = (seed_re, seed_im);
    let mut converged = false;
    let mut iterations = 0;
    for iter in 0..64 {
        iterations = iter + 1;
        let (mut dre, mut dim) = (0.0f64, 0.0f64);
        let (mut d2re, mut d2im) = (0.0f64, 0.0f64);
        for term in terms {
            let a = term.exponent.get();
            let c = term.coefficient.get();
            let magnitude = (-TWO_PI * a * im).exp();
            if magnitude == 0.0 {
                continue;
            }
            let phase = TWO_PI * a * re;
            let (z_re, z_im) = (magnitude * phase.cos(), magnitude * phase.sin());
            let (w_re, w_im) = (1.0 - z_re, -z_im);
            let w_norm2 = w_re.mul_add(w_re, w_im * w_im);
            if w_norm2 <= 0.0 {
                println!("  iter {iter}: w_norm2 <= 0 (z hit 1)");
                return (iterations, re, im, false, (im - seed_im).abs() / seed_im);
            }
            let z_norm2 = z_re.mul_add(z_re, z_im * z_im);
            let ln_re = 0.5 * z_re.mul_add(-2.0, z_norm2).ln_1p();
            let ln_im = w_im.atan2(w_re);
            let f = TWO_PI * a * c;
            dre += f * ln_im;
            dim -= f * ln_re;
            let (q_re, q_im) = (
                (z_re * w_re + z_im * w_im) / w_norm2,
                (z_im * w_re - z_re * w_im) / w_norm2,
            );
            let f2 = -(TWO_PI * a) * (TWO_PI * a) * c;
            d2re += f2 * q_re;
            d2im += f2 * q_im;
        }
        let d2_norm2 = d2re.mul_add(d2re, d2im * d2im);
        if d2_norm2 == 0.0 {
            println!("  iter {iter}: second derivative vanished");
            return (iterations, re, im, false, (im - seed_im).abs() / seed_im);
        }
        let step_re = dre.mul_add(d2re, dim * d2im) / d2_norm2;
        let step_im = dim.mul_add(d2re, -(dre * d2im)) / d2_norm2;
        re -= step_re;
        im -= step_im;
        if !(re.is_finite() && im.is_finite()) {
            println!("  iter {iter}: step went non-finite");
            return (iterations, re, im, false, (im - seed_im).abs() / seed_im);
        }
        if iter < 8 || iter % 8 == 0 {
            println!(
                "  iter {iter}: re={re:.6e} im={im:.6e} |step|={:.3e}",
                step_re.abs() + step_im.abs()
            );
        }
        if step_re.abs() + step_im.abs() < 1e-14 * im.abs() {
            converged = true;
            break;
        }
    }
    (
        iterations,
        re,
        im,
        converged,
        (im - seed_im).abs() / seed_im,
    )
}

#[test]
fn diagnose_racetrack_refinement() {
    let Ok(spec) = std::env::var("CYRUS_RACETRACK_DIAG") else {
        eprintln!("CYRUS_RACETRACK_DIAG not set; skipping");
        return;
    };
    let parts: Vec<&str> = spec.split('|').collect();
    let (pool_path, name) = (parts[0], parts[1]);
    let parse = |s: &str| -> Vec<i64> { s.split(',').map(|x| x.parse().unwrap()).collect() };
    let (k, m) = (parse(parts[2]), parse(parts[3]));

    let pool = load_pool(std::path::Path::new(pool_path)).expect("pool");
    let record = pool.iter().find(|p| p.name == name).expect("polytope");
    let geom =
        GaGeometry::prepare_from_points(&record.points, DEFAULT_GV_MIN_POINTS).expect("geometry");

    // Flat direction p = N^{-1} K, as the pipeline computes it.
    let m_typed: Vec<I64<Finite>> = m.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let k_typed: Vec<I64<Finite>> = k.iter().map(|&x| I64::<Finite>::new(x)).collect();
    let n_mat = cyrus_core::flat_direction::compute_n_matrix(&geom.kappa_basis, &m_typed);
    let p = cyrus_core::flat_direction::solve_linear_system_faer(&n_mat, &k_typed)
        .expect("N invertible");

    let terms = build_racetrack_terms(&geom.gv, &m_typed, &p, f64_pos!(1.0));
    println!("racetrack terms: {}", terms.len());
    let two_term = solve_racetrack(&terms).expect("two-term racetrack converges (GA saw it)");
    println!(
        "two-term seed: g_s={:.6} Im(tau)={:.4} Re(tau)={:.4}",
        two_term.g_s.get(),
        two_term.im_tau.get(),
        two_term.re_tau.get()
    );

    let (iters, re, im, converged, drift) =
        newton_with_telemetry(two_term.re_tau.get(), two_term.im_tau.get(), &terms);
    println!("\nVERDICT:");
    println!("  newton converged: {converged} ({iters} iterations)");
    println!(
        "  final: Re(tau)={re:.6} Im(tau)={im:.6} -> g_s={:.6}",
        1.0 / im
    );
    println!(
        "  Im(tau) drift from two-term seed: {:.2}% (runner rejects > 5%)",
        drift * 100.0
    );
    if converged && drift <= 0.05 {
        println!("  => refinement SHOULD have passed; suspect a numerics bug in the runner path");
    } else if converged {
        println!(
            "  => CONVERGED but the full series moved the vacuum {:.1}%: the GA's two-term values are untrustworthy here",
            drift * 100.0
        );
    } else {
        println!("  => Newton DIVERGED: no controlled stationary point near the two-term seed");
    }

    // Re-evaluate the candidate through the (now refinement-aware) fast
    // chain and print the corrected report.
    let cfg = cyrus_ga::fitness::FitnessConfig::default();
    let genome = cyrus_ga::genome::Genome { k, m };
    let report = cyrus_ga::fitness::evaluate_fitness(&geom, &cfg, &genome);
    println!(
        "\nREFINED EVALUATION: tier={} fitness={:.2} g_s={:?} logV0={:?} cpl=({:?},{:?})",
        report.tier, report.fitness, report.g_s, report.log10_abs_v0, report.cpl_w0, report.cpl_wa
    );
}
