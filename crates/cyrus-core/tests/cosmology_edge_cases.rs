#![allow(missing_docs)]
use cyrus_core::cosmology::{CosmologyParams, Potential, solve_cosmology};

struct LcdmPotential;
impl Potential for LcdmPotential {
    fn value(&self, _phi: f64) -> f64 {
        1.0
    }
    fn deriv(&self, _phi: f64) -> f64 {
        0.0
    }
}

struct ZeroPotential;
impl Potential for ZeroPotential {
    fn value(&self, _phi: f64) -> f64 {
        0.0
    }
    fn deriv(&self, _phi: f64) -> f64 {
        0.0
    }
}

struct LinearPotential;
impl Potential for LinearPotential {
    fn value(&self, phi: f64) -> f64 {
        phi
    }
    fn deriv(&self, _phi: f64) -> f64 {
        1.0
    }
}

#[test]
fn test_cosmology_zero_potential() {
    // Zero potential with matter only
    let params = CosmologyParams {
        omega_m0: 1.0,
        omega_de0: 0.0,
        h0: 1.0,
    };
    let pot = ZeroPotential;

    let res = solve_cosmology(&params, &pot, 0.0, 0.0, 1.0).unwrap();
    assert!(!res.redshifts.is_empty());
}

#[test]
fn test_cosmology_low_z_start() {
    // Very low starting redshift
    let params = CosmologyParams {
        omega_m0: 0.3,
        omega_de0: 0.7,
        h0: 1.0,
    };
    let pot = LcdmPotential;

    let res = solve_cosmology(&params, &pot, 0.0, 0.0, 0.1).unwrap();
    assert!(!res.redshifts.is_empty());
}

#[test]
fn test_cosmology_with_velocity() {
    // Start with non-zero field velocity
    let params = CosmologyParams {
        omega_m0: 0.3,
        omega_de0: 0.7,
        h0: 1.0,
    };
    let pot = LcdmPotential;

    // Small velocity
    let res = solve_cosmology(&params, &pot, 1.0, 0.1, 1.0).unwrap();
    assert!(!res.redshifts.is_empty());
}

#[test]
fn test_cosmology_linear_potential() {
    // Linear (quintessence) potential
    let params = CosmologyParams {
        omega_m0: 0.3,
        omega_de0: 0.7,
        h0: 1.0,
    };
    let pot = LinearPotential;

    let res = solve_cosmology(&params, &pot, 1.0, 0.0, 1.0).unwrap();
    assert!(!res.redshifts.is_empty());
    // w should deviate from -1 for rolling field
}

#[test]
fn test_cosmology_high_velocity() {
    // High velocity that triggers early return (phi'^2 >= 6)
    // When 3 - 0.5 * phi'^2 becomes small/negative
    let params = CosmologyParams {
        omega_m0: 0.3,
        omega_de0: 0.7,
        h0: 1.0,
    };
    let pot = LcdmPotential;

    // Very high initial velocity - triggers the h_sq <= 0 branch
    let res = solve_cosmology(&params, &pot, 0.0, 3.0, 1.0).unwrap();
    // Should still return valid result (ODE handles the early return gracefully)
    assert!(!res.redshifts.is_empty());
}

struct NegativePotential;
impl Potential for NegativePotential {
    fn value(&self, _phi: f64) -> f64 {
        -10.0 // Very negative potential
    }
    fn deriv(&self, _phi: f64) -> f64 {
        0.0
    }
}

#[test]
fn test_cosmology_negative_potential() {
    // Negative potential leading to potentially negative H^2
    let params = CosmologyParams {
        omega_m0: 0.1,
        omega_de0: 0.0,
        h0: 1.0,
    };
    let pot = NegativePotential;

    // Should handle gracefully
    let res = solve_cosmology(&params, &pot, 0.0, 0.0, 1.0);
    // May succeed or fail depending on how the ODE handles it
    if let Ok(r) = res {
        assert!(!r.redshifts.is_empty());
    }
}
