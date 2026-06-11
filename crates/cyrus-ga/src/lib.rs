//! Genetic-algorithm landscape search over flux vacua.
//!
//! Searches integer flux pairs `(K, M)` on a fixed Calabi-Yau geometry for
//! perturbatively flat vacua whose racetrack-stabilized energy scale and
//! thawing-quintessence behavior match the DESI evolving-dark-energy
//! measurements. Physics evaluation is `cyrus_core::evaluate_vacuum` (the
//! validated Demirtas-McAllister scan chain); cosmological scoring
//! integrates the real Friedmann + Klein-Gordon equations via
//! `cyrus_core::cosmology` and CPL-fits w(z) against DESI (w0, wa).
//!
//! Runs persist their full state (population, hall of fame, RNG) to a run
//! directory every generation, so they can be stopped and resumed exactly.

pub mod fitness;
pub mod genome;
pub mod geometry;
pub mod population;
