//! Generic operator implementations using the algebra rules.
//!
//! This macro generates Add, Sub, Mul, Div, Neg implementations for any
//! wrapper type that holds a numeric inner value.

/// Generate all operator implementations for a wrapper type.
///
/// Usage: `impl_ops!(F64, f64);` or `impl_ops!(I32, i32);`
#[macro_export]
macro_rules! impl_ops {
    ($Wrapper:ident, $Inner:ty) => {
        // ====================================================================
        // Sum (for iterators) - only for self-additive types (T + T = T)
        // ====================================================================

        impl<T> ::std::iter::Sum for $Wrapper<T>
        where
            T: $crate::types::algebra::AddOutput<T, Output = T>,
            $Inner: Default,
        {
            fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
                // Note: empty sum returns zero, which may violate Pos constraint
                // This is consistent with std::iter::Sum behavior
                $Wrapper(iter.map(|x| x.0).sum(), ::std::marker::PhantomData)
            }
        }
        // ====================================================================
        // Add
        // ====================================================================

        impl<L, R> ::std::ops::Add<$Wrapper<R>> for $Wrapper<L>
        where
            L: $crate::types::algebra::AddOutput<R>,
        {
            type Output = $Wrapper<L::Output>;

            #[inline]
            fn add(self, rhs: $Wrapper<R>) -> Self::Output {
                $Wrapper(self.0 + rhs.0, ::std::marker::PhantomData)
            }
        }

        // ====================================================================
        // Sub
        // ====================================================================

        impl<L, R> ::std::ops::Sub<$Wrapper<R>> for $Wrapper<L>
        where
            L: $crate::types::algebra::SubOutput<R>,
        {
            type Output = $Wrapper<L::Output>;

            #[inline]
            fn sub(self, rhs: $Wrapper<R>) -> Self::Output {
                $Wrapper(self.0 - rhs.0, ::std::marker::PhantomData)
            }
        }

        // ====================================================================
        // Mul
        // ====================================================================

        impl<L, R> ::std::ops::Mul<$Wrapper<R>> for $Wrapper<L>
        where
            L: $crate::types::algebra::MulOutput<R>,
        {
            type Output = $Wrapper<L::Output>;

            #[inline]
            fn mul(self, rhs: $Wrapper<R>) -> Self::Output {
                $Wrapper(self.0 * rhs.0, ::std::marker::PhantomData)
            }
        }

        // ====================================================================
        // Div
        // ====================================================================

        impl<L, R> ::std::ops::Div<$Wrapper<R>> for $Wrapper<L>
        where
            L: $crate::types::algebra::DivOutput<R>,
        {
            type Output = $Wrapper<L::Output>;

            #[inline]
            fn div(self, rhs: $Wrapper<R>) -> Self::Output {
                $Wrapper(self.0 / rhs.0, ::std::marker::PhantomData)
            }
        }

        // ====================================================================
        // Neg
        // ====================================================================

        impl<T> ::std::ops::Neg for $Wrapper<T>
        where
            T: $crate::types::algebra::NegOutput,
        {
            type Output = $Wrapper<T::Output>;

            #[inline]
            fn neg(self) -> Self::Output {
                $Wrapper(-self.0, ::std::marker::PhantomData)
            }
        }
    };
}

pub use impl_ops;
