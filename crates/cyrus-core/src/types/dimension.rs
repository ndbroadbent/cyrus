//! Dimension-safe vectors and matrices with HRTB branding.
//!
//! This module provides compile-time dimension safety using branded types.
//! The key insight is that phantom lifetimes create unforgeable "brands" that
//! tie vectors and indices to their originating dimension.
//!
//! # Core Concepts
//!
//! - **[`Dim`]**: A runtime dimension value that can create branded contexts
//! - **[`Brand`]**: An unforgeable phantom type tied to a specific dimension context
//! - **[`BVec`]**: A vector branded with a dimension - same brand = same length
//! - **[`Index`]**: A branded index that is guaranteed valid for its dimension
//! - **[`BMat`]**: A matrix with separate row and column brands
//!
//! # Example
//!
//! ```
//! use cyrus_core::types::dimension::{Dim, BVec};
//!
//! let dim = Dim::new(3);
//!
//! dim.with_brand(|brand| {
//!     let a = BVec::new(&dim, &brand, vec![1.0, 2.0, 3.0]);
//!     let b = BVec::new(&dim, &brand, vec![4.0, 5.0, 6.0]);
//!
//!     // Safe iteration - no bounds checks needed in release
//!     let dot: f64 = dim.indices(&brand)
//!         .map(|i| a.get(i) * b.get(i))
//!         .sum();
//!
//!     assert!((dot - 32.0).abs() < 1e-10);
//! });
//! ```
//!
//! # Safety Guarantees
//!
//! - `BVec<'id, T>` and `BVec<'id, U>` with the same brand have the same length
//! - `Index<'id>` is always valid for `BVec<'id, T>` (cannot be out of bounds)
//! - Brands cannot escape their `with_brand` closure (HRTB constraint)
//! - Indices cannot be forged from raw `usize` (private field)

use std::marker::PhantomData;
use std::ops::{Add, Mul};

// ============================================================================
// Brand (Unforgeable Phantom Type)
// ============================================================================

/// An unforgeable brand tied to a specific dimension context.
///
/// `Brand<'id>` uses an invariant phantom lifetime to prevent:
/// - Forging brands from outside `with_brand`
/// - Mixing brands from different `with_brand` calls
/// - Escaping brands from their closure
///
/// The `*mut &'id ()` trick makes the lifetime invariant (not covariant),
/// which is critical for soundness.
#[derive(Debug)]
pub struct Brand<'id>(PhantomData<*mut &'id ()>);

// Brand is NOT Clone or Copy - this prevents duplication

// ============================================================================
// Dim (Dimension Handle)
// ============================================================================

/// A runtime dimension value that can create branded contexts.
///
/// `Dim` wraps a `usize` representing the dimension (number of elements).
/// Use [`Dim::with_brand`] to enter a branded context where vectors and
/// indices are statically guaranteed to match this dimension.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Dim {
    n: usize,
}

impl Dim {
    /// Create a new dimension.
    #[must_use]
    pub const fn new(n: usize) -> Self {
        Self { n }
    }

    /// Get the dimension value.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.n
    }

    /// Enter a branded context.
    ///
    /// The closure receives a [`Brand`] that cannot escape. This brand
    /// ties all [`BVec`]s and [`Index`]es created within to this dimension.
    ///
    /// # Example
    ///
    /// ```
    /// use cyrus_core::types::dimension::{Dim, BVec};
    ///
    /// let dim = Dim::new(3);
    /// let sum = dim.with_brand(|brand| {
    ///     let v = BVec::new(&dim, &brand, vec![1, 2, 3]);
    ///     dim.indices(&brand).map(|i| *v.get(i)).sum::<i32>()
    /// });
    /// assert_eq!(sum, 6);
    /// ```
    pub fn with_brand<R, F>(&self, f: F) -> R
    where
        F: for<'id> FnOnce(Brand<'id>) -> R,
    {
        // Create a fresh brand that is tied to this call
        f(Brand(PhantomData))
    }

    /// Create an iterator over valid indices for this dimension.
    ///
    /// Each [`Index`] yielded is guaranteed to be in bounds for any
    /// [`BVec`] created with the same brand.
    pub fn indices<'id>(&self, _brand: &Brand<'id>) -> impl Iterator<Item = Index<'id>> {
        (0..self.n).map(|i| Index {
            i,
            _brand: PhantomData,
        })
    }

    /// Try to create a valid index from a raw `usize`.
    ///
    /// Returns `None` if `i >= self.n`.
    #[must_use]
    pub fn index<'id>(&self, _brand: &Brand<'id>, i: usize) -> Option<Index<'id>> {
        (i < self.n).then(|| Index {
            i,
            _brand: PhantomData,
        })
    }

    /// Create a zero-initialized vector with this dimension.
    #[must_use]
    pub fn zeros<'id, T: Default + Clone>(&self, _brand: &Brand<'id>) -> BVec<'id, T> {
        BVec {
            data: vec![T::default(); self.n],
            _brand: PhantomData,
        }
    }

    /// Create a vector filled with a value.
    #[must_use]
    pub fn fill<'id, T: Clone>(&self, _brand: &Brand<'id>, value: T) -> BVec<'id, T> {
        BVec {
            data: vec![value; self.n],
            _brand: PhantomData,
        }
    }
}

// ============================================================================
// Index (Branded Index)
// ============================================================================

/// A branded index that is guaranteed valid for vectors with the same brand.
///
/// `Index<'id>` can only be created through [`Dim::indices`] or [`Dim::index`],
/// ensuring it's always in bounds for `BVec<'id, T>`.
#[derive(Debug, Clone, Copy)]
pub struct Index<'id> {
    i: usize,
    _brand: PhantomData<Brand<'id>>,
}

impl<'id> Index<'id> {
    /// Get the raw index value.
    ///
    /// This is useful for interfacing with external code, but loses
    /// the safety guarantee.
    #[must_use]
    pub const fn raw(&self) -> usize {
        self.i
    }
}

// PartialEq compares the underlying index
impl<'id> PartialEq for Index<'id> {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i
    }
}

impl<'id> Eq for Index<'id> {}

impl<'id> std::hash::Hash for Index<'id> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.i.hash(state);
    }
}

// ============================================================================
// BVec (Branded Vector)
// ============================================================================

/// A vector branded with a dimension lifetime.
///
/// All `BVec<'id, T>` with the same `'id` have the same length, and
/// `Index<'id>` is always valid for indexing into them.
#[derive(Debug, Clone)]
pub struct BVec<'id, T> {
    data: Vec<T>,
    _brand: PhantomData<Brand<'id>>,
}

impl<'id, T> BVec<'id, T> {
    /// Create a new branded vector.
    ///
    /// # Panics
    ///
    /// Panics if `data.len() != dim.n()`.
    #[must_use]
    pub fn new(dim: &Dim, _brand: &Brand<'id>, data: Vec<T>) -> Self {
        assert_eq!(
            data.len(),
            dim.n(),
            "BVec length {} doesn't match dimension {}",
            data.len(),
            dim.n()
        );
        Self {
            data,
            _brand: PhantomData,
        }
    }

    /// Try to create a branded vector, returning `None` if length mismatches.
    #[must_use]
    pub fn try_new(dim: &Dim, _brand: &Brand<'id>, data: Vec<T>) -> Option<Self> {
        (data.len() == dim.n()).then(|| Self {
            data,
            _brand: PhantomData,
        })
    }

    /// Get the length of the vector.
    #[must_use]
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Check if the vector is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get a reference to an element by branded index.
    ///
    /// The brand system guarantees the index is valid.
    #[must_use]
    pub fn get(&self, idx: Index<'id>) -> &T {
        debug_assert!(idx.i < self.data.len(), "brand system bug: index out of bounds");
        &self.data[idx.i]
    }

    /// Get a mutable reference to an element by branded index.
    #[must_use]
    pub fn get_mut(&mut self, idx: Index<'id>) -> &mut T {
        debug_assert!(idx.i < self.data.len(), "brand system bug: index out of bounds");
        &mut self.data[idx.i]
    }

    /// Set an element by branded index.
    pub fn set(&mut self, idx: Index<'id>, value: T) {
        *self.get_mut(idx) = value;
    }

    /// Get the underlying data as a slice.
    #[must_use]
    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    /// Get the underlying data as a mutable slice.
    #[must_use]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }

    /// Consume and return the underlying data.
    #[must_use]
    pub fn into_inner(self) -> Vec<T> {
        self.data
    }

    /// Iterate over references to elements.
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    /// Iterate over mutable references to elements.
    pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.data.iter_mut()
    }

    /// Map each element to a new value.
    #[must_use]
    pub fn map<U, F: FnMut(&T) -> U>(&self, f: F) -> BVec<'id, U> {
        BVec {
            data: self.data.iter().map(f).collect(),
            _brand: PhantomData,
        }
    }

    /// Zip with another vector of the same brand.
    pub fn zip<'a, U>(&'a self, other: &'a BVec<'id, U>) -> impl Iterator<Item = (&'a T, &'a U)> {
        // Same brand guarantees same length
        self.data.iter().zip(other.data.iter())
    }
}

impl<'id, T: Clone> BVec<'id, T> {
    /// Map each element in place.
    pub fn map_in_place<F: FnMut(&T) -> T>(&mut self, mut f: F) {
        for x in &mut self.data {
            *x = f(x);
        }
    }
}

impl<'id, T: Copy + Add<Output = T>> BVec<'id, T> {
    /// Element-wise addition of two vectors with the same brand.
    #[must_use]
    pub fn add(&self, other: &BVec<'id, T>) -> BVec<'id, T> {
        // Same brand guarantees same length - no check needed
        BVec {
            data: self
                .data
                .iter()
                .zip(other.data.iter())
                .map(|(&a, &b)| a + b)
                .collect(),
            _brand: PhantomData,
        }
    }
}

impl<'id, T: Copy + Mul<Output = T>> BVec<'id, T> {
    /// Element-wise multiplication of two vectors with the same brand.
    #[must_use]
    pub fn mul(&self, other: &BVec<'id, T>) -> BVec<'id, T> {
        BVec {
            data: self
                .data
                .iter()
                .zip(other.data.iter())
                .map(|(&a, &b)| a * b)
                .collect(),
            _brand: PhantomData,
        }
    }
}

impl<'id, T: Copy + Mul<Output = T> + Add<Output = T> + Default> BVec<'id, T> {
    /// Dot product of two vectors with the same brand.
    #[must_use]
    pub fn dot(&self, other: &BVec<'id, T>) -> T {
        self.data
            .iter()
            .zip(other.data.iter())
            .map(|(&a, &b)| a * b)
            .fold(T::default(), |acc, x| acc + x)
    }
}

// Standard indexing with branded index
impl<'id, T> std::ops::Index<Index<'id>> for BVec<'id, T> {
    type Output = T;

    fn index(&self, idx: Index<'id>) -> &Self::Output {
        self.get(idx)
    }
}

impl<'id, T> std::ops::IndexMut<Index<'id>> for BVec<'id, T> {
    fn index_mut(&mut self, idx: Index<'id>) -> &mut Self::Output {
        self.get_mut(idx)
    }
}

// ============================================================================
// Row and Column Indices for Matrices
// ============================================================================

/// A branded row index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct RowIdx<'row> {
    i: usize,
    _brand: PhantomData<Brand<'row>>,
}

impl<'row> RowIdx<'row> {
    /// Get the raw index value.
    #[must_use]
    pub const fn raw(&self) -> usize {
        self.i
    }
}

/// A branded column index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ColIdx<'col> {
    j: usize,
    _brand: PhantomData<Brand<'col>>,
}

impl<'col> ColIdx<'col> {
    /// Get the raw index value.
    #[must_use]
    pub const fn raw(&self) -> usize {
        self.j
    }
}

// ============================================================================
// BMat (Branded Matrix)
// ============================================================================

/// A matrix with separate row and column brands.
///
/// `BMat<'row, 'col, T>` guarantees that:
/// - Row indices from `rows.indices(row_brand)` are valid
/// - Column indices from `cols.indices(col_brand)` are valid
/// - Multiplying with `BVec<'col, T>` produces `BVec<'row, T>`
#[derive(Debug, Clone)]
pub struct BMat<'row, 'col, T> {
    data: Vec<T>,       // row-major: data[row * ncols + col]
    nrows: usize,
    ncols: usize,
    _row_brand: PhantomData<Brand<'row>>,
    _col_brand: PhantomData<Brand<'col>>,
}

impl<'row, 'col, T> BMat<'row, 'col, T> {
    /// Create a new branded matrix.
    ///
    /// # Panics
    ///
    /// Panics if `data.len() != rows.n() * cols.n()`.
    #[must_use]
    pub fn new(
        rows: &Dim,
        cols: &Dim,
        _row_brand: &Brand<'row>,
        _col_brand: &Brand<'col>,
        data: Vec<T>,
    ) -> Self {
        let expected = rows.n() * cols.n();
        assert_eq!(
            data.len(),
            expected,
            "BMat data length {} doesn't match {}x{}={}",
            data.len(),
            rows.n(),
            cols.n(),
            expected
        );
        Self {
            data,
            nrows: rows.n(),
            ncols: cols.n(),
            _row_brand: PhantomData,
            _col_brand: PhantomData,
        }
    }

    /// Create a zero-initialized matrix.
    #[must_use]
    pub fn zeros(
        rows: &Dim,
        cols: &Dim,
        _row_brand: &Brand<'row>,
        _col_brand: &Brand<'col>,
    ) -> Self
    where
        T: Default + Clone,
    {
        Self {
            data: vec![T::default(); rows.n() * cols.n()],
            nrows: rows.n(),
            ncols: cols.n(),
            _row_brand: PhantomData,
            _col_brand: PhantomData,
        }
    }

    /// Get the number of rows.
    #[must_use]
    pub const fn nrows(&self) -> usize {
        self.nrows
    }

    /// Get the number of columns.
    #[must_use]
    pub const fn ncols(&self) -> usize {
        self.ncols
    }

    /// Get a reference to an element.
    #[must_use]
    pub fn get(&self, row: RowIdx<'row>, col: ColIdx<'col>) -> &T {
        debug_assert!(row.i < self.nrows && col.j < self.ncols, "brand system bug: index out of bounds");
        &self.data[row.i * self.ncols + col.j]
    }

    /// Get a mutable reference to an element.
    #[must_use]
    pub fn get_mut(&mut self, row: RowIdx<'row>, col: ColIdx<'col>) -> &mut T {
        debug_assert!(row.i < self.nrows && col.j < self.ncols, "brand system bug: index out of bounds");
        &mut self.data[row.i * self.ncols + col.j]
    }

    /// Set an element.
    pub fn set(&mut self, row: RowIdx<'row>, col: ColIdx<'col>, value: T) {
        *self.get_mut(row, col) = value;
    }

    /// Get a row as a vector.
    #[must_use]
    pub fn row(&self, row: RowIdx<'row>) -> BVec<'col, T>
    where
        T: Clone,
    {
        let start = row.i * self.ncols;
        let end = start + self.ncols;
        BVec {
            data: self.data[start..end].to_vec(),
            _brand: PhantomData,
        }
    }

    /// Get the underlying data as a slice.
    #[must_use]
    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    /// Consume and return the underlying data.
    #[must_use]
    pub fn into_inner(self) -> Vec<T> {
        self.data
    }
}

/// Helper to create row indices from a dimension.
pub struct RowDim {
    dim: Dim,
}

impl RowDim {
    /// Create a new row dimension.
    #[must_use]
    pub const fn new(n: usize) -> Self {
        Self { dim: Dim::new(n) }
    }

    /// Get the dimension value.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.dim.n()
    }

    /// Enter a branded row context.
    pub fn with_brand<R, F>(&self, f: F) -> R
    where
        F: for<'row> FnOnce(Brand<'row>) -> R,
    {
        self.dim.with_brand(f)
    }

    /// Create an iterator over valid row indices.
    pub fn indices<'row>(&self, _brand: &Brand<'row>) -> impl Iterator<Item = RowIdx<'row>> {
        (0..self.dim.n()).map(|i| RowIdx {
            i,
            _brand: PhantomData,
        })
    }

    /// Get the inner Dim reference.
    #[must_use]
    pub const fn dim(&self) -> &Dim {
        &self.dim
    }
}

/// Helper to create column indices from a dimension.
pub struct ColDim {
    dim: Dim,
}

impl ColDim {
    /// Create a new column dimension.
    #[must_use]
    pub const fn new(n: usize) -> Self {
        Self { dim: Dim::new(n) }
    }

    /// Get the dimension value.
    #[must_use]
    pub const fn n(&self) -> usize {
        self.dim.n()
    }

    /// Enter a branded column context.
    pub fn with_brand<R, F>(&self, f: F) -> R
    where
        F: for<'col> FnOnce(Brand<'col>) -> R,
    {
        self.dim.with_brand(f)
    }

    /// Create an iterator over valid column indices.
    pub fn indices<'col>(&self, _brand: &Brand<'col>) -> impl Iterator<Item = ColIdx<'col>> {
        (0..self.dim.n()).map(|j| ColIdx {
            j,
            _brand: PhantomData,
        })
    }

    /// Get the inner Dim reference.
    #[must_use]
    pub const fn dim(&self) -> &Dim {
        &self.dim
    }
}

/// Type alias for a square matrix (same row and column brand).
pub type BSquareMat<'id, T> = BMat<'id, 'id, T>;

// ============================================================================
// Matrix-Vector Operations
// ============================================================================

impl<'row, 'col, T> BMat<'row, 'col, T>
where
    T: Copy + Mul<Output = T> + Add<Output = T> + Default,
{
    /// Multiply matrix by column vector.
    ///
    /// The type system guarantees `vec.len() == self.ncols` via the brand.
    #[must_use]
    pub fn mul_vec(&self, vec: &BVec<'col, T>) -> BVec<'row, T> {
        let result: Vec<T> = (0..self.nrows)
            .map(|i| {
                let row_start = i * self.ncols;
                self.data[row_start..row_start + self.ncols]
                    .iter()
                    .zip(vec.as_slice().iter())
                    .map(|(&m, &v)| m * v)
                    .fold(T::default(), |acc, x| acc + x)
            })
            .collect();
        BVec {
            data: result,
            _brand: PhantomData,
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // --- Dim tests ---

    #[test]
    fn test_dim_new() {
        let d = Dim::new(5);
        assert_eq!(d.n(), 5);
    }

    #[test]
    fn test_dim_with_brand() {
        let d = Dim::new(3);
        let result = d.with_brand(|_brand| 42);
        assert_eq!(result, 42);
    }

    #[test]
    fn test_dim_indices() {
        let d = Dim::new(4);
        d.with_brand(|brand| {
            let indices: Vec<usize> = d.indices(&brand).map(|i| i.raw()).collect();
            assert_eq!(indices, vec![0, 1, 2, 3]);
        });
    }

    #[test]
    fn test_dim_index() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            assert!(d.index(&brand, 0).is_some());
            assert!(d.index(&brand, 2).is_some());
            assert!(d.index(&brand, 3).is_none());
            assert!(d.index(&brand, 100).is_none());
        });
    }

    #[test]
    fn test_dim_zeros() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let v: BVec<'_, i32> = d.zeros(&brand);
            assert_eq!(v.as_slice(), &[0, 0, 0]);
        });
    }

    #[test]
    fn test_dim_fill() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let v = d.fill(&brand, 7);
            assert_eq!(v.as_slice(), &[7, 7, 7]);
        });
    }

    // --- Index tests ---

    #[test]
    fn test_index_raw() {
        let d = Dim::new(5);
        d.with_brand(|brand| {
            let idx = d.index(&brand, 3).unwrap();
            assert_eq!(idx.raw(), 3);
        });
    }

    #[test]
    fn test_index_eq() {
        let d = Dim::new(5);
        d.with_brand(|brand| {
            let i1 = d.index(&brand, 2).unwrap();
            let i2 = d.index(&brand, 2).unwrap();
            let i3 = d.index(&brand, 3).unwrap();
            assert_eq!(i1, i2);
            assert_ne!(i1, i3);
        });
    }

    #[test]
    fn test_index_hash() {
        use std::collections::HashSet;
        let d = Dim::new(5);
        d.with_brand(|brand| {
            let mut set = HashSet::new();
            set.insert(d.index(&brand, 1).unwrap());
            set.insert(d.index(&brand, 2).unwrap());
            set.insert(d.index(&brand, 1).unwrap());
            assert_eq!(set.len(), 2);
        });
    }

    // --- BVec tests ---

    #[test]
    fn test_bvec_new() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let v = BVec::new(&d, &brand, vec![1, 2, 3]);
            assert_eq!(v.len(), 3);
            assert!(!v.is_empty());
        });
    }

    #[test]
    #[should_panic(expected = "doesn't match dimension")]
    fn test_bvec_new_wrong_len() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let _v = BVec::new(&d, &brand, vec![1, 2]); // Wrong length!
        });
    }

    #[test]
    fn test_bvec_try_new() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            assert!(BVec::try_new(&d, &brand, vec![1, 2, 3]).is_some());
            assert!(BVec::try_new(&d, &brand, vec![1, 2]).is_none());
        });
    }

    #[test]
    fn test_bvec_get_set() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let mut v = BVec::new(&d, &brand, vec![1, 2, 3]);
            let idx = d.index(&brand, 1).unwrap();

            assert_eq!(*v.get(idx), 2);
            v.set(idx, 20);
            assert_eq!(*v.get(idx), 20);
        });
    }

    #[test]
    fn test_bvec_index_ops() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let mut v = BVec::new(&d, &brand, vec![1, 2, 3]);
            let idx = d.index(&brand, 0).unwrap();

            assert_eq!(v[idx], 1);
            v[idx] = 10;
            assert_eq!(v[idx], 10);
        });
    }

    #[test]
    fn test_bvec_slices() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let mut v = BVec::new(&d, &brand, vec![1, 2, 3]);

            assert_eq!(v.as_slice(), &[1, 2, 3]);

            v.as_mut_slice()[1] = 20;
            assert_eq!(v.as_slice(), &[1, 20, 3]);
        });
    }

    #[test]
    fn test_bvec_into_inner() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let v = BVec::new(&d, &brand, vec![1, 2, 3]);
            assert_eq!(v.into_inner(), vec![1, 2, 3]);
        });
    }

    #[test]
    fn test_bvec_iter() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let v = BVec::new(&d, &brand, vec![1, 2, 3]);
            let sum: i32 = v.iter().sum();
            assert_eq!(sum, 6);
        });
    }

    #[test]
    fn test_bvec_iter_mut() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let mut v = BVec::new(&d, &brand, vec![1, 2, 3]);
            for x in v.iter_mut() {
                *x *= 2;
            }
            assert_eq!(v.as_slice(), &[2, 4, 6]);
        });
    }

    #[test]
    fn test_bvec_map() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let v = BVec::new(&d, &brand, vec![1, 2, 3]);
            let v2 = v.map(|x| x * 2);
            assert_eq!(v2.as_slice(), &[2, 4, 6]);
        });
    }

    #[test]
    fn test_bvec_map_in_place() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let mut v = BVec::new(&d, &brand, vec![1, 2, 3]);
            v.map_in_place(|x| x * 2);
            assert_eq!(v.as_slice(), &[2, 4, 6]);
        });
    }

    #[test]
    fn test_bvec_zip() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let a = BVec::new(&d, &brand, vec![1, 2, 3]);
            let b = BVec::new(&d, &brand, vec![4, 5, 6]);
            let pairs: Vec<_> = a.zip(&b).map(|(&x, &y)| (x, y)).collect();
            assert_eq!(pairs, vec![(1, 4), (2, 5), (3, 6)]);
        });
    }

    #[test]
    fn test_bvec_add() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let a = BVec::new(&d, &brand, vec![1, 2, 3]);
            let b = BVec::new(&d, &brand, vec![4, 5, 6]);
            let c = a.add(&b);
            assert_eq!(c.as_slice(), &[5, 7, 9]);
        });
    }

    #[test]
    fn test_bvec_mul() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let a = BVec::new(&d, &brand, vec![1, 2, 3]);
            let b = BVec::new(&d, &brand, vec![4, 5, 6]);
            let c = a.mul(&b);
            assert_eq!(c.as_slice(), &[4, 10, 18]);
        });
    }

    #[test]
    fn test_bvec_dot() {
        let d = Dim::new(3);
        d.with_brand(|brand| {
            let a = BVec::new(&d, &brand, vec![1.0_f64, 2.0, 3.0]);
            let b = BVec::new(&d, &brand, vec![4.0_f64, 5.0, 6.0]);
            let dot: f64 = a.dot(&b);
            assert!((dot - 32.0).abs() < 1e-10);
        });
    }

    #[test]
    fn test_bvec_empty() {
        let d = Dim::new(0);
        d.with_brand(|brand| {
            let v: BVec<'_, i32> = BVec::new(&d, &brand, vec![]);
            assert!(v.is_empty());
            assert_eq!(v.len(), 0);
        });
    }

    // --- BMat tests ---

    #[test]
    fn test_bmat_new() {
        let rows = Dim::new(2);
        let cols = Dim::new(3);

        rows.with_brand(|rb| {
            cols.with_brand(|cb| {
                let m = BMat::new(&rows, &cols, &rb, &cb, vec![1, 2, 3, 4, 5, 6]);
                assert_eq!(m.nrows(), 2);
                assert_eq!(m.ncols(), 3);
            });
        });
    }

    #[test]
    #[should_panic(expected = "doesn't match")]
    fn test_bmat_new_wrong_size() {
        let rows = Dim::new(2);
        let cols = Dim::new(3);

        rows.with_brand(|rb| {
            cols.with_brand(|cb| {
                let _m = BMat::new(&rows, &cols, &rb, &cb, vec![1, 2, 3, 4, 5]); // Wrong!
            });
        });
    }

    #[test]
    fn test_bmat_zeros() {
        let rows = Dim::new(2);
        let cols = Dim::new(3);

        rows.with_brand(|rb| {
            cols.with_brand(|cb| {
                let m: BMat<'_, '_, i32> = BMat::zeros(&rows, &cols, &rb, &cb);
                assert_eq!(m.as_slice(), &[0, 0, 0, 0, 0, 0]);
            });
        });
    }

    #[test]
    fn test_bmat_get_set() {
        let rows = RowDim::new(2);
        let cols = ColDim::new(3);

        rows.with_brand(|rb| {
            cols.with_brand(|cb| {
                let mut m = BMat::new(rows.dim(), cols.dim(), &rb, &cb, vec![1, 2, 3, 4, 5, 6]);

                let r0: RowIdx<'_> = rows.indices(&rb).next().unwrap();
                let c1: ColIdx<'_> = cols.indices(&cb).nth(1).unwrap();

                assert_eq!(*m.get(r0, c1), 2);
                m.set(r0, c1, 20);
                assert_eq!(*m.get(r0, c1), 20);
            });
        });
    }

    #[test]
    fn test_bmat_row() {
        let rows = RowDim::new(2);
        let cols = ColDim::new(3);

        rows.with_brand(|rb| {
            cols.with_brand(|cb| {
                let m = BMat::new(rows.dim(), cols.dim(), &rb, &cb, vec![1, 2, 3, 4, 5, 6]);
                let r1 = rows.indices(&rb).nth(1).unwrap();
                let row = m.row(r1);
                assert_eq!(row.as_slice(), &[4, 5, 6]);
            });
        });
    }

    #[test]
    fn test_bmat_mul_vec() {
        let rows = Dim::new(2);
        let cols = Dim::new(3);

        rows.with_brand(|rb| {
            cols.with_brand(|cb| {
                // Matrix: [[1,2,3],[4,5,6]]
                let m = BMat::new(&rows, &cols, &rb, &cb, vec![1, 2, 3, 4, 5, 6]);
                let v = BVec::new(&cols, &cb, vec![1, 0, 1]);

                let result = m.mul_vec(&v);
                // [1*1+2*0+3*1, 4*1+5*0+6*1] = [4, 10]
                assert_eq!(result.as_slice(), &[4, 10]);
            });
        });
    }

    #[test]
    fn test_bmat_into_inner() {
        let rows = Dim::new(2);
        let cols = Dim::new(2);

        rows.with_brand(|rb| {
            cols.with_brand(|cb| {
                let m = BMat::new(&rows, &cols, &rb, &cb, vec![1, 2, 3, 4]);
                assert_eq!(m.into_inner(), vec![1, 2, 3, 4]);
            });
        });
    }

    // --- RowDim/ColDim tests ---

    #[test]
    fn test_rowdim_coldim() {
        let rows = RowDim::new(2);
        let cols = ColDim::new(3);

        assert_eq!(rows.n(), 2);
        assert_eq!(cols.n(), 3);

        rows.with_brand(|rb| {
            let row_indices: Vec<usize> = rows.indices(&rb).map(|r| r.raw()).collect();
            assert_eq!(row_indices, vec![0, 1]);
        });

        cols.with_brand(|cb| {
            let col_indices: Vec<usize> = cols.indices(&cb).map(|c| c.raw()).collect();
            assert_eq!(col_indices, vec![0, 1, 2]);
        });
    }

    // --- Square matrix tests ---

    #[test]
    fn test_square_mat() {
        let d = Dim::new(2);

        d.with_brand(|brand| {
            // For square matrices, use the same brand for rows and cols
            // This requires some restructuring since we can't use the same brand twice
            // In practice, you'd use a single with_brand for square matrices
            let m: BSquareMat<'_, i32> = BMat {
                data: vec![1, 2, 3, 4],
                nrows: 2,
                ncols: 2,
                _row_brand: PhantomData,
                _col_brand: PhantomData,
            };
            assert_eq!(m.nrows(), 2);
            assert_eq!(m.ncols(), 2);

            // Can use same index for row and col
            let v = BVec::new(&d, &brand, vec![1, 2]);
            let result = m.mul_vec(&v);
            // [1*1+2*2, 3*1+4*2] = [5, 11]
            assert_eq!(result.as_slice(), &[5, 11]);
        });
    }

}
