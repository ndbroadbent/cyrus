//! Configuration module for Cyrus.
//!
//! This is a direct port of CYTools `config.py`.
//!
//! ## Features
//!
//! - Thread configuration for parallel computations
//! - Experimental features flag
//! - Optimizer backend configuration
//!
//! ## Example
//!
//! ```rust
//! use cyrus_core::config;
//!
//! // Set number of threads (None = use all available)
//! config::set_n_threads(Some(4));
//!
//! // Enable experimental features
//! config::enable_experimental_features();
//! ```

use std::sync::RwLock;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

// =============================================================================
// Thread Configuration
// =============================================================================

/// The number of CPU threads to use in parallel computations.
/// When set to 0, uses all available threads (via rayon's default).
static N_THREADS: AtomicUsize = AtomicUsize::new(0);

/// Get the number of threads to use for parallel computations.
///
/// Returns `None` if all available threads should be used,
/// or `Some(n)` for a specific number of threads.
///
/// # Example
///
/// ```rust
/// use cyrus_core::config;
///
/// let threads = config::n_threads();
/// match threads {
///     None => println!("Using all available threads"),
///     Some(n) => println!("Using {} threads", n),
/// }
/// ```
#[must_use]
pub fn n_threads() -> Option<usize> {
    match N_THREADS.load(Ordering::Relaxed) {
        0 => None,
        n => Some(n),
    }
}

/// Set the number of threads to use for parallel computations.
///
/// Pass `None` to use all available threads, or `Some(n)` for a specific number.
///
/// # Example
///
/// ```rust
/// use cyrus_core::config;
///
/// // Use 4 threads
/// config::set_n_threads(Some(4));
///
/// // Use all available threads
/// config::set_n_threads(None);
/// ```
pub fn set_n_threads(threads: Option<usize>) {
    N_THREADS.store(threads.unwrap_or(0), Ordering::Relaxed);
}

/// Get the effective number of threads to use.
///
/// If `n_threads()` returns `None`, this returns the number of available CPUs.
/// Otherwise returns the configured thread count.
#[must_use]
pub fn effective_n_threads() -> usize {
    n_threads().unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(std::num::NonZeroUsize::get)
            .unwrap_or(1)
    })
}

// =============================================================================
// Experimental Features
// =============================================================================

/// Whether experimental features are enabled.
static EXP_FEATURES_ENABLED: AtomicBool = AtomicBool::new(false);

/// Check if experimental features are enabled.
///
/// # Example
///
/// ```rust
/// use cyrus_core::config;
///
/// if config::experimental_features_enabled() {
///     // Use experimental feature
/// }
/// ```
#[must_use]
pub fn experimental_features_enabled() -> bool {
    EXP_FEATURES_ENABLED.load(Ordering::Relaxed)
}

/// Enable experimental features.
///
/// # Warning
///
/// Experimental features may be broken, not fully tested, or undergo
/// significant changes in future versions.
///
/// # Example
///
/// ```rust
/// use cyrus_core::config;
///
/// config::enable_experimental_features();
/// assert!(config::experimental_features_enabled());
/// ```
pub fn enable_experimental_features() {
    EXP_FEATURES_ENABLED.store(true, Ordering::Relaxed);
    eprintln!(
        "\n**************************************************************\n\
         Warning: You have enabled experimental features of Cyrus.\n\
         Some of these features may be broken or not fully tested,\n\
         and they may undergo significant changes in future versions.\n\
         **************************************************************\n"
    );
}

/// Disable experimental features.
///
/// # Example
///
/// ```rust
/// use cyrus_core::config;
///
/// config::disable_experimental_features();
/// assert!(!config::experimental_features_enabled());
/// ```
pub fn disable_experimental_features() {
    EXP_FEATURES_ENABLED.store(false, Ordering::Relaxed);
}

// =============================================================================
// Optimizer Backend
// =============================================================================

/// Available optimizer backends for cone computations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum OptimizerBackend {
    /// OR-Tools (Google's Operations Research tools) - default
    #[default]
    OrTools,
    /// GLPK (GNU Linear Programming Kit)
    Glpk,
    /// HiGHS optimizer
    HiGhs,
    /// OSQP (Operator Splitting Quadratic Program)
    Osqp,
}

/// The current optimizer backend.
static OPTIMIZER_BACKEND: RwLock<OptimizerBackend> = RwLock::new(OptimizerBackend::OrTools);

/// Get the current optimizer backend.
///
/// # Example
///
/// ```rust
/// use cyrus_core::config::{self, OptimizerBackend};
///
/// let backend = config::optimizer_backend();
/// ```
#[must_use]
pub fn optimizer_backend() -> OptimizerBackend {
    *OPTIMIZER_BACKEND.read().expect("RwLock poisoned")
}

/// Set the optimizer backend.
///
/// # Example
///
/// ```rust
/// use cyrus_core::config::{self, OptimizerBackend};
///
/// config::set_optimizer_backend(OptimizerBackend::HiGhs);
/// ```
pub fn set_optimizer_backend(backend: OptimizerBackend) {
    *OPTIMIZER_BACKEND.write().expect("RwLock poisoned") = backend;
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_n_threads_default() {
        // Reset to default
        set_n_threads(None);
        assert!(n_threads().is_none());
    }

    #[test]
    fn test_n_threads_set() {
        set_n_threads(Some(4));
        assert_eq!(n_threads(), Some(4));
        // Reset
        set_n_threads(None);
    }

    #[test]
    fn test_effective_n_threads() {
        set_n_threads(None);
        let effective = effective_n_threads();
        assert!(effective > 0);

        set_n_threads(Some(2));
        assert_eq!(effective_n_threads(), 2);

        // Reset
        set_n_threads(None);
    }

    #[test]
    fn test_experimental_features_default() {
        // Note: This test may fail if run after enable_experimental_features
        // in another test. Tests should be independent.
        disable_experimental_features();
        assert!(!experimental_features_enabled());
    }

    #[test]
    fn test_experimental_features_enable_disable() {
        disable_experimental_features();
        assert!(!experimental_features_enabled());

        // We avoid calling enable_experimental_features() in tests
        // because it prints a warning to stderr
        EXP_FEATURES_ENABLED.store(true, Ordering::Relaxed);
        assert!(experimental_features_enabled());

        disable_experimental_features();
        assert!(!experimental_features_enabled());
    }

    #[test]
    fn test_optimizer_backend_default() {
        assert_eq!(optimizer_backend(), OptimizerBackend::OrTools);
    }

    #[test]
    fn test_optimizer_backend_set() {
        set_optimizer_backend(OptimizerBackend::HiGhs);
        assert_eq!(optimizer_backend(), OptimizerBackend::HiGhs);

        // Reset to default
        set_optimizer_backend(OptimizerBackend::OrTools);
    }
}
