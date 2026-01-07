//! Tensor operations for symmetric sparse tensors.
//!
//! Direct port of CYTools `utils.py` tensor functions.

use std::collections::BTreeMap;

/// A sparse symmetric tensor stored as a map from sorted index tuples to values.
///
/// For a rank-3 tensor, indices (i,j,k) are stored with i <= j <= k.
pub type SparseTensor<T> = BTreeMap<Vec<usize>, T>;

/// Filter a tensor to only include specified indices, reindexing to 0..len(indices).
///
/// This function selects a specific subset of indices from a tensor and reindexes
/// them so that indices are in the range `0..indices.len()` with the ordering
/// specified by the input indices. This is used to convert intersection numbers
/// to a given basis.
///
/// # Arguments
///
/// * `tensor` - Sparse symmetric tensor as `{(i,j,...): value}`
/// * `indices` - List of indices to keep (defines the new ordering)
///
/// # Returns
///
/// A new sparse tensor with only the specified indices, reindexed.
///
/// # Example
///
/// ```rust
/// use std::collections::BTreeMap;
/// use cyrus_core::utils::filter_tensor_indices;
///
/// let mut tensor: BTreeMap<Vec<usize>, i64> = BTreeMap::new();
/// tensor.insert(vec![0, 1], 0);
/// tensor.insert(vec![1, 1], 1);
/// tensor.insert(vec![1, 2], 2);
/// tensor.insert(vec![1, 3], 3);
/// tensor.insert(vec![2, 3], 4);
///
/// let filtered = filter_tensor_indices(&tensor, &[1, 3]);
/// // Original indices [1,3] map to new indices [0,1]
/// // (1,1) -> (0,0) with value 1
/// // (1,3) -> (0,1) with value 3
///
/// assert_eq!(filtered.get(&vec![0, 0]), Some(&1));
/// assert_eq!(filtered.get(&vec![0, 1]), Some(&3));
/// assert_eq!(filtered.len(), 2);
/// ```
pub fn filter_tensor_indices<T: Clone>(
    tensor: &SparseTensor<T>,
    indices: &[usize],
) -> SparseTensor<T> {
    // Build reindex map: old_index -> new_index
    let reindex: BTreeMap<usize, usize> = indices
        .iter()
        .enumerate()
        .map(|(new_idx, &old_idx)| (old_idx, new_idx))
        .collect();

    let mut result = SparseTensor::new();

    for (key, val) in tensor {
        // Check if all indices are in our filter set
        if key.iter().all(|idx| reindex.contains_key(idx)) {
            // Reindex and sort
            let mut new_key: Vec<usize> = key.iter().map(|idx| reindex[idx]).collect();
            new_key.sort_unstable();
            result.insert(new_key, val.clone());
        }
    }

    result
}

/// Convert a sparse symmetric tensor to a dense array.
///
/// The sparse tensor has indices stored in sorted order (i <= j <= ... <= k).
/// The dense output fills in all permutations of each index tuple.
///
/// # Arguments
///
/// * `tensor` - Sparse symmetric tensor
/// * `dim` - Dimension of each axis (if None, inferred from tensor)
///
/// # Returns
///
/// A dense tensor as a nested Vec structure.
///
/// # Example
///
/// ```rust
/// use std::collections::BTreeMap;
/// use cyrus_core::utils::symmetric_sparse_to_dense;
///
/// let mut tensor: BTreeMap<Vec<usize>, i64> = BTreeMap::new();
/// tensor.insert(vec![0, 0], 1);
/// tensor.insert(vec![0, 1], 2);
/// tensor.insert(vec![1, 1], 3);
///
/// let dense = symmetric_sparse_to_dense(&tensor, Some(2));
/// assert_eq!(dense, vec![vec![1, 2], vec![2, 3]]);
/// ```
pub fn symmetric_sparse_to_dense<T: Clone + Default>(
    tensor: &SparseTensor<T>,
    dim: Option<usize>,
) -> Vec<Vec<T>> {
    if tensor.is_empty() {
        return Vec::new();
    }

    // Get rank from first key
    let rank = tensor.keys().next().map_or(0, Vec::len);
    assert!(
        (rank == 2),
        "symmetric_sparse_to_dense currently only supports rank-2 tensors"
    );

    // Determine dimension
    let dim = dim.unwrap_or_else(|| {
        tensor
            .keys()
            .flat_map(|k| k.iter())
            .max()
            .map_or(0, |&m| m + 1)
    });

    // Create dense output
    let mut out = vec![vec![T::default(); dim]; dim];

    // Fill in values (symmetric)
    for (key, val) in tensor {
        let i = key[0];
        let j = key[1];
        out[i][j] = val.clone();
        out[j][i] = val.clone();
    }

    out
}

/// Convert a dense symmetric tensor to sparse format.
///
/// Only stores upper triangular entries (i <= j <= ... <= k) with non-zero values.
///
/// # Arguments
///
/// * `tensor` - Dense symmetric tensor (2D)
///
/// # Returns
///
/// Sparse tensor with sorted index tuples.
///
/// # Example
///
/// ```rust
/// use cyrus_core::utils::symmetric_dense_to_sparse;
///
/// let dense = vec![vec![1, 2], vec![2, 3]];
/// let sparse = symmetric_dense_to_sparse(&dense);
///
/// assert_eq!(sparse.get(&vec![0, 0]), Some(&1));
/// assert_eq!(sparse.get(&vec![0, 1]), Some(&2));
/// assert_eq!(sparse.get(&vec![1, 1]), Some(&3));
/// assert_eq!(sparse.len(), 3);
/// ```
pub fn symmetric_dense_to_sparse(tensor: &[Vec<i64>]) -> SparseTensor<i64> {
    let mut out = SparseTensor::new();

    let dim = tensor.len();
    for i in 0..dim {
        for j in i..dim {
            let val = tensor[i][j];
            if val != 0 {
                out.insert(vec![i, j], val);
            }
        }
    }

    out
}

/// Convert a rank-3 sparse symmetric tensor to dense format.
///
/// # Arguments
///
/// * `tensor` - Sparse rank-3 symmetric tensor
/// * `dim` - Dimension of each axis
///
/// # Returns
///
/// Dense 3D tensor as nested Vecs.
pub fn symmetric_sparse_to_dense_3d<T: Clone + Default>(
    tensor: &SparseTensor<T>,
    dim: usize,
) -> Vec<Vec<Vec<T>>> {
    let mut out = vec![vec![vec![T::default(); dim]; dim]; dim];

    for (key, val) in tensor {
        if key.len() != 3 {
            continue;
        }
        let (i, j, k) = (key[0], key[1], key[2]);

        // Fill all 6 permutations for symmetric tensor
        out[i][j][k] = val.clone();
        out[i][k][j] = val.clone();
        out[j][i][k] = val.clone();
        out[j][k][i] = val.clone();
        out[k][i][j] = val.clone();
        out[k][j][i] = val.clone();
    }

    out
}

/// Convert a dense rank-3 symmetric tensor to sparse format.
///
/// Only stores entries with i <= j <= k.
pub fn symmetric_dense_to_sparse_3d(tensor: &[Vec<Vec<i64>>]) -> SparseTensor<i64> {
    let mut out = SparseTensor::new();

    let dim = tensor.len();
    for i in 0..dim {
        for j in i..dim {
            for k in j..dim {
                let val = tensor[i][j][k];
                if val != 0 {
                    out.insert(vec![i, j, k], val);
                }
            }
        }
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_filter_tensor_indices_basic() {
        let mut tensor: SparseTensor<i64> = BTreeMap::new();
        tensor.insert(vec![0, 1], 0);
        tensor.insert(vec![1, 1], 1);
        tensor.insert(vec![1, 2], 2);
        tensor.insert(vec![1, 3], 3);
        tensor.insert(vec![2, 3], 4);

        let filtered = filter_tensor_indices(&tensor, &[1, 3]);

        assert_eq!(filtered.len(), 2);
        assert_eq!(filtered.get(&vec![0, 0]), Some(&1)); // (1,1) -> (0,0)
        assert_eq!(filtered.get(&vec![0, 1]), Some(&3)); // (1,3) -> (0,1)
    }

    #[test]
    fn test_filter_tensor_indices_empty() {
        let tensor: SparseTensor<i64> = BTreeMap::new();
        let filtered = filter_tensor_indices(&tensor, &[1, 3]);
        assert!(filtered.is_empty());
    }

    #[test]
    fn test_filter_tensor_indices_no_match() {
        let mut tensor: SparseTensor<i64> = BTreeMap::new();
        tensor.insert(vec![0, 1], 5);

        let filtered = filter_tensor_indices(&tensor, &[2, 3]);
        assert!(filtered.is_empty());
    }

    #[test]
    fn test_symmetric_sparse_to_dense_2d() {
        let mut tensor: SparseTensor<i64> = BTreeMap::new();
        tensor.insert(vec![0, 0], 1);
        tensor.insert(vec![0, 1], 2);
        tensor.insert(vec![1, 1], 3);

        let dense = symmetric_sparse_to_dense(&tensor, Some(2));

        assert_eq!(dense, vec![vec![1, 2], vec![2, 3]]);
    }

    #[test]
    fn test_symmetric_dense_to_sparse_2d() {
        let dense = vec![vec![1, 2], vec![2, 3]];
        let sparse = symmetric_dense_to_sparse(&dense);

        assert_eq!(sparse.len(), 3);
        assert_eq!(sparse.get(&vec![0, 0]), Some(&1));
        assert_eq!(sparse.get(&vec![0, 1]), Some(&2));
        assert_eq!(sparse.get(&vec![1, 1]), Some(&3));
    }

    #[test]
    fn test_symmetric_roundtrip_2d() {
        let mut original: SparseTensor<i64> = BTreeMap::new();
        original.insert(vec![0, 0], 5);
        original.insert(vec![0, 1], 3);
        original.insert(vec![1, 1], 7);

        let dense = symmetric_sparse_to_dense(&original, Some(2));
        let sparse = symmetric_dense_to_sparse(&dense);

        assert_eq!(original, sparse);
    }

    #[test]
    fn test_symmetric_sparse_to_dense_3d() {
        let mut tensor: SparseTensor<i64> = BTreeMap::new();
        tensor.insert(vec![0, 0, 0], 1);
        tensor.insert(vec![0, 0, 1], 2);
        tensor.insert(vec![0, 1, 1], 3);
        tensor.insert(vec![1, 1, 1], 4);

        let dense = symmetric_sparse_to_dense_3d(&tensor, 2);

        // Check diagonal
        assert_eq!(dense[0][0][0], 1);
        assert_eq!(dense[1][1][1], 4);

        // Check symmetry
        assert_eq!(dense[0][0][1], 2);
        assert_eq!(dense[0][1][0], 2);
        assert_eq!(dense[1][0][0], 2);
    }

    #[test]
    fn test_symmetric_roundtrip_3d() {
        let mut original: SparseTensor<i64> = BTreeMap::new();
        original.insert(vec![0, 0, 0], 1);
        original.insert(vec![0, 1, 1], 2);
        original.insert(vec![1, 1, 1], 3);

        let dense = symmetric_sparse_to_dense_3d(&original, 2);
        let sparse = symmetric_dense_to_sparse_3d(&dense);

        assert_eq!(original, sparse);
    }
}
