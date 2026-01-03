//! Quintessence cosmology solver.
//!
//! Solves the background evolution of a scalar field rolling down a potential
//! in an expanding universe.
//!
//! Equations of motion (Friedmann + Klein-Gordon):
//! 1. H^2 = (rho_m + rho_phi) / 3
//! 2. d^2phi/dt^2 + 3H dphi/dt + V'(phi) = 0
//! 3. rho_m = rho_m0 * a^-3
//! 4. rho_phi = 0.5 * dot^2 + V(phi)
//!
//! We solve for w(z) = P_phi / rho_phi.

use crate::error::{Error, Result};
use ode_solvers::rk4::Rk4;
use ode_solvers::{System, Vector2};

/// Cosmological parameters (Planck units M_p = 1)
#[derive(Debug, Clone)]
pub struct CosmologyParams {
    /// Matter density today
    pub omega_m0: f64,
    /// Dark energy density today (target)
    pub omega_de0: f64,
    /// Hubble parameter today (in eV or similar units, but usually we scale H0=1)
    pub h0: f64,
}

/// Quintessence potential trait
pub trait Potential: Send + Sync {
    /// Value V(phi)
    fn value(&self, phi: f64) -> f64;
    /// Derivative V'(phi)
    fn deriv(&self, phi: f64) -> f64;
}

/// State of the universe: [ln(a), phi, dphi/dln(a)]
/// Using e-folds N = ln(a) as time variable is numerically more stable.
///
/// Variables: y[0] = phi, y[1] = dphi/dN.
///
/// H(N)^2 = (rho_m0 * e^(-3N) + rho_phi) / 3
/// rho_phi = 0.5 * H^2 * (dphi/dN)^2 + V(phi)
///
/// Solving for H:
/// H^2 (1 - (phi')^2/6) = V/3 + rho_m/3
/// H^2 = (V + rho_m) / (3 - 0.5 * (phi')^2)
///
/// KG Equation in N:
/// phi'' + (3 + H'/H) phi' + V'/H^2 = 0
///
/// Friedmann 2:
/// H'/H = -1.5 * (1 + w_eff)
/// w_eff = P_tot / rho_tot = (P_phi) / (rho_m + rho_phi)
/// P_phi = 0.5 * H^2 * (phi')^2 - V
///
type State = Vector2<f64>;

struct QuintessenceSystem<'a> {
    params: &'a CosmologyParams,
    potential: &'a dyn Potential,
}

impl System<f64, State> for QuintessenceSystem<'_> {
    fn system(&self, n: f64, y: &State, dy: &mut State) {
        let phi = y[0];
        let phi_prime = y[1]; // dphi/dN

        let v = self.potential.value(phi);
        let dv = self.potential.deriv(phi);

        // Matter density rho_m = 3 * H0^2 * Om0 * exp(-3N)
        let rho_m = 3.0 * self.params.h0.powi(2) * self.params.omega_m0 * (-3.0 * n).exp();

        // Hubble parameter H^2
        let h_sq = (rho_m + v) / 0.5f64.mul_add(-phi_prime.powi(2), 3.0);

        if h_sq <= 0.0 || (0.5f64.mul_add(-phi_prime.powi(2), 3.0)).abs() < 1e-6 {
            dy[0] = 0.0;
            dy[1] = 0.0;
            return;
        }

        // Kinetic energy = 0.5 * H^2 * phi'^2
        let kinetic = 0.5 * h_sq * phi_prime.powi(2);
        let p_tot = kinetic - v;
        let rho_tot = 3.0 * h_sq;
        let w_eff = p_tot / rho_tot;

        dy[0] = phi_prime;
        // phi'' = - (3 + H'/H) phi' - V'/H^2
        // 3 + H'/H = 1.5 (1 - w_eff) ??
        // H'/H = -1.5 (1 + w_eff)
        // So 3 + H'/H = 3 - 1.5 - 1.5 w_eff = 1.5 (1 - w_eff)
        // Wait, check sign.
        // H^2 ~ a^(-3(1+w)). 2 H H' = -3(1+w) a^{-...} a' = -3(1+w) H^2 * H ?
        // d(H^2)/dN = 2 H dH/dN = -3(1+w) H^2.
        // dH/dN = -1.5 (1+w) H. Correct.
        let friction = 1.5 * (1.0 - w_eff);
        dy[1] = (-friction).mul_add(phi_prime, -(dv / h_sq));
    }
}

/// Result of cosmological evolution
#[derive(Debug)]
pub struct CosmologyResult {
    /// Redshifts z
    pub redshifts: Vec<f64>,
    /// Equation of state parameter w(z)
    pub w_z: Vec<f64>,
    /// Dark energy density parameter Omega_phi(z)
    pub omega_phi_z: Vec<f64>,
}

/// Validate the starting redshift for cosmology integration.
///
/// # Arguments
/// * `z_start` - Starting redshift
///
/// # Errors
/// Returns an error if:
/// - `z_start <= 0`: Zero or negative redshift is unphysical and causes ODE solver issues
///   - `z_start = 0` gives `n_start = n_end = 0` (zero interval, causes ODE solver panic)
///   - `z_start < 0` is unphysical (negative redshift has no cosmological meaning)
///   - `z_start = -1` gives `ln(0) = -inf` (causes ODE solver panic)
pub fn validate_z_start(z_start: f64) -> Result<()> {
    if z_start <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "z_start must be positive, got {z_start}"
        )));
    }
    Ok(())
}

/// Solve the cosmology from z_initial to z_final (usually z=0).
///
/// # Arguments
/// * `params` - Cosmological parameters
/// * `potential` - Scalar potential V(phi)
/// * `phi_i` - Initial field value
/// * `phi_prime_i` - Initial field velocity dphi/dN
/// * `z_start` - Starting redshift (must be > 0, e.g. 1000)
///
/// # Errors
/// Returns an error if z_start <= 0 or if ODE integration fails.
pub fn solve_cosmology(
    params: &CosmologyParams,
    potential: &dyn Potential,
    phi_i: f64,
    phi_prime_i: f64,
    z_start: f64,
) -> Result<CosmologyResult> {
    validate_z_start(z_start)?;

    let n_start = -z_start.ln_1p();
    let n_end = 0.0; // z=0
    let dx = (n_end - n_start) / 100.0; // Ensure step is positive and reasonable

    let system = QuintessenceSystem { params, potential };
    let state = Vector2::new(phi_i, phi_prime_i);

    let mut stepper = Rk4::new(system, n_start, state, n_end, dx);

    // Integrate
    let res = stepper.integrate();

    match res {
        Ok(_) => {
            let mut redshifts = Vec::new();
            let mut w_z = Vec::new();
            let mut omega_phi_z = Vec::new();

            for (n, y) in stepper.x_out().iter().zip(stepper.y_out().iter()) {
                let z = (-n).exp_m1().abs(); // Ensure positive z due to numerical noise around 0
                let phi = y[0];
                let phi_prime = y[1];

                let v = potential.value(phi);
                let rho_m = 3.0 * params.h0.powi(2) * params.omega_m0 * (-3.0 * n).exp();
                let h_sq = (rho_m + v) / 0.5f64.mul_add(-phi_prime.powi(2), 3.0);

                let kinetic = 0.5 * h_sq * phi_prime.powi(2);
                let rho_phi = kinetic + v;
                let p_phi = kinetic - v;

                let w = if rho_phi.abs() > 1e-10 {
                    p_phi / rho_phi
                } else {
                    -1.0
                };
                let omega = rho_phi / (3.0 * h_sq);

                redshifts.push(z);
                w_z.push(w);
                omega_phi_z.push(omega);
            }

            Ok(CosmologyResult {
                redshifts,
                w_z,
                omega_phi_z,
            })
        }
        Err(_) => Err(Error::LinearAlgebra("ODE integration failed".into())), // Reuse error variant
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    struct LcdmPotential;
    impl Potential for LcdmPotential {
        fn value(&self, _phi: f64) -> f64 {
            1.0
        } // Constant Lambda
        fn deriv(&self, _phi: f64) -> f64 {
            0.0
        }
    }

    #[test]
    fn test_lcdm_evolution() {
        let params = CosmologyParams {
            omega_m0: 0.3,
            omega_de0: 0.7,
            h0: 1.0,
        };
        let pot = LcdmPotential;

        let res = solve_cosmology(&params, &pot, 0.0, 0.0, 1.0).unwrap();

        assert!(!res.redshifts.is_empty());
        // w(z) should be close to -1 for constant potential
        for w in res.w_z {
            assert!((w - (-1.0)).abs() < 0.01);
        }
    }

    #[test]
    fn test_validate_z_start_positive() {
        // Valid positive redshifts
        assert!(validate_z_start(1.0).is_ok());
        assert!(validate_z_start(0.1).is_ok());
        assert!(validate_z_start(1000.0).is_ok());
        assert!(validate_z_start(1e-9).is_ok()); // Very small but positive
    }

    #[test]
    fn test_validate_z_start_zero() {
        // z_start = 0 gives n_start = n_end = 0 (zero integration interval)
        let res = validate_z_start(0.0);
        assert!(res.is_err());
        assert!(
            res.unwrap_err()
                .to_string()
                .contains("z_start must be positive")
        );
    }

    #[test]
    fn test_validate_z_start_negative() {
        // Negative redshift is unphysical
        let res = validate_z_start(-0.5);
        assert!(res.is_err());
        assert!(
            res.unwrap_err()
                .to_string()
                .contains("z_start must be positive")
        );
    }

    #[test]
    fn test_validate_z_start_minus_one() {
        // z_start = -1 would give ln(0) = -inf
        let res = validate_z_start(-1.0);
        assert!(res.is_err());
        assert!(
            res.unwrap_err()
                .to_string()
                .contains("z_start must be positive")
        );
    }

    #[test]
    fn test_validate_z_start_large_negative() {
        // Very negative values
        let res = validate_z_start(-100.0);
        assert!(res.is_err());
    }
}
