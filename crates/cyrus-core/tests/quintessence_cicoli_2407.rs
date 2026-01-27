//! Reproduce the Cicoli et al. 2407.03405 numerical example scales.

use cyrus_core::{f64_pos, i64_pos};
use cyrus_core::quintessence::{
    axion_masses, derive_cicoli_2407, late_time_gradient_2d, late_time_potential_2d,
    AxionPotential1D, Cicoli2407Params,
};
use cyrus_core::cosmology::{solve_cosmology, CosmologyParams};
use cyrus_core::types::f64::F64;
use cyrus_core::types::tags::Finite;

fn assert_rel_close(actual: f64, expected: f64, rel_tol: f64, label: &str) {
    let denom = expected.abs().max(1.0);
    let rel_err = (actual - expected).abs() / denom;
    assert!(
        rel_err <= rel_tol,
        "{label} rel_err={rel_err} (actual={actual}, expected={expected})"
    );
}

fn assert_in_range(actual: f64, low: f64, high: f64, label: &str) {
    assert!(
        actual >= low && actual <= high,
        "{label} out of range: {actual} not in [{low}, {high}]"
    );
}

#[test]
fn cicoli_2407_example_matches_paper_scales() {
    // Parameters from arXiv:2407.03405, Sec. 4.5 (rounded values).
    let params = Cicoli2407Params {
        k: f64_pos!(1.0),
        hat_k: f64_pos!(1.0),
        g_s: f64_pos!(0.1),
        xi: f64_pos!(0.5),
        w0: f64_pos!(1.0),
        a1_prefactor: f64_pos!(1.0),
        a2_prefactor: f64_pos!(1.0),
        n1: i64_pos!(1),
        n2: i64_pos!(27),
    };

    let derived = derive_cicoli_2407(&params, f64_pos!(0.085), f64_pos!(0.0038));

    // Paper values (rounded): tau1 ~ 1.324, tau2 ~ 1122, V ~ 900.
    assert_rel_close(derived.tau_1.get(), 1.324, 0.02, "tau_1");
    assert_rel_close(derived.tau_2.get(), 1122.0, 0.02, "tau_2");
    assert_rel_close(derived.volume.get(), 900.0, 0.03, "volume");

    // Decay constants should match the target values by construction.
    assert_rel_close(derived.f1.get(), 0.085, 1e-12, "f1");
    assert_rel_close(derived.f2.get(), 0.0038, 1e-12, "f2");

    // Hierarchy and order-of-magnitude checks for the late-time potential.
    assert!(derived.lambda_2_4 > derived.lambda_1_4);
    assert_in_range(derived.lambda_1_4.get(), 1.0e-123, 1.0e-120, "lambda_1_4");
    assert_in_range(derived.lambda_2_4.get(), 1.0e-119, 1.0e-116, "lambda_2_4");

    let (m1, m2) = axion_masses(derived.lambda_1_4, derived.lambda_2_4, derived.f1, derived.f2);
    assert_in_range(m1.get(), 1.0e-60, 2.0e-59, "m1");
    assert_in_range(m2.get(), 1.0e-57, 3.0e-56, "m2");
}

#[test]
fn axion_potential_1d_sanity_checks() {
    let pot = AxionPotential1D {
        lambda_4: f64_pos!(1.0e-120),
        f: f64_pos!(0.1),
    };

    let phi0 = F64::<Finite>::new(0.0).unwrap();
    let v0 = pot.value_typed(phi0).get();
    let dv0 = pot.deriv_typed(phi0).get();
    assert!(v0.abs() < 1e-14, "V(0) should be ~0");
    assert!(dv0.abs() < 1e-14, "V'(0) should be ~0");

    let phi_pi = F64::<Finite>::new(0.1 * std::f64::consts::PI).unwrap();
    let v_pi = pot.value_typed(phi_pi).get();
    assert_rel_close(v_pi, 2.0e-120, 1e-12, "V(pi f)");
}

#[test]
fn axion_potential_2d_minimum_checks() {
    let lambda_1_4 = f64_pos!(1.0e-120);
    let lambda_2_4 = f64_pos!(2.0e-118);
    let f1 = f64_pos!(0.1);
    let f2 = f64_pos!(0.01);

    let phi_0 = F64::<Finite>::new(0.0).unwrap();
    let v0 = late_time_potential_2d(phi_0, phi_0, f1, f2, lambda_1_4, lambda_2_4).get();
    let (dphi1, dphi2) = late_time_gradient_2d(phi_0, phi_0, f1, f2, lambda_1_4, lambda_2_4);
    assert!(v0.abs() < 1e-14, "V(0,0) should be ~0");
    assert!(dphi1.get().abs() < 1e-14, "dV/dphi1 at minimum should be ~0");
    assert!(dphi2.get().abs() < 1e-14, "dV/dphi2 at minimum should be ~0");

    let phi2_pi = F64::<Finite>::new(f2.get() * std::f64::consts::PI).unwrap();
    let v_pi = late_time_potential_2d(phi_0, phi2_pi, f1, f2, lambda_1_4, lambda_2_4).get();
    assert_rel_close(v_pi, 2.0 * lambda_2_4.get(), 1e-12, "V(0, pi f2)");
}

#[test]
fn axion_potential_cosmology_w_is_minus_one() {
    let params = CosmologyParams {
        omega_m0: 0.3,
        omega_de0: 0.7,
        h0: 1.0,
    };
    let f = f64_pos!(0.1);
    let pot = AxionPotential1D {
        lambda_4: f64_pos!(1.0e-2),
        f,
    };
    let phi_i = f.get() * std::f64::consts::PI;
    let result = solve_cosmology(&params, &pot, phi_i, 0.0, 10.0)
        .expect("cosmology solver should succeed");

    let w_start = result.w_z.first().copied().unwrap_or(0.0);
    let w_end = result.w_z.last().copied().unwrap_or(0.0);
    assert_rel_close(w_start, -1.0, 1e-6, "w(z_start)");
    assert_rel_close(w_end, -1.0, 1e-6, "w(z=0)");
}
