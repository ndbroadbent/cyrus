//! Interval types for bounded floating-point values.
//!
//! `UnitInterval` guarantees a value in the range (0, 1].

/// A floating-point value guaranteed to be in the interval (0, 1].
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct UnitInterval(f64);

impl UnitInterval {
    /// Create a new unit interval value. Returns None if val <= 0 or val > 1.
    pub fn new(val: f64) -> Option<Self> {
        if val > 0.0 && val <= 1.0 {
            Some(Self(val))
        } else {
            None
        }
    }

    /// Get the inner value.
    pub fn get(self) -> f64 {
        self.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unit_interval_rejects_zero() {
        assert!(UnitInterval::new(0.0).is_none());
    }

    #[test]
    fn test_unit_interval_rejects_negative() {
        assert!(UnitInterval::new(-0.5).is_none());
    }

    #[test]
    fn test_unit_interval_rejects_greater_than_one() {
        assert!(UnitInterval::new(1.1).is_none());
    }

    #[test]
    fn test_unit_interval_accepts_one() {
        let u = UnitInterval::new(1.0).unwrap();
        assert_eq!(u.get(), 1.0);
    }

    #[test]
    fn test_unit_interval_accepts_half() {
        let u = UnitInterval::new(0.5).unwrap();
        assert_eq!(u.get(), 0.5);
    }

    #[test]
    fn test_unit_interval_accepts_small_positive() {
        let u = UnitInterval::new(0.001).unwrap();
        assert_eq!(u.get(), 0.001);
    }
}
