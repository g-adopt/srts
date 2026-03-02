"""Coefficient format conversion between Fortran flat arrays and pyshtools cilm.

Fortran flat: (lmax+1)^2 values ordered [c_00, c_10, c_11c, c_11s, c_20, c_21c, c_21s, ...]
pyshtools cilm: [2, lmax+1, lmax+1] with cilm[0,l,m] = cosine, cilm[1,l,m] = sine.

The .raw file includes a 0.01 scaling and normylm factor (sqrt(2) for m!=0).
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np


@lru_cache(maxsize=8)
def _index_tables(lmax: int) -> dict:
    """Precompute all index mapping arrays using vectorized numpy operations.

    Returns dict with:
        flat_cos_idx: flat indices of cosine terms, shape (n_lm,) where n_lm = (lmax+1)(lmax+2)/2
        flat_sin_idx: flat indices of sine terms (m>0 only), shape (n_mgt0,)
        l_vals:       l value for each (l,m) pair, shape (n_lm,)
        m_vals:       m value for each (l,m) pair, shape (n_lm,)
        mgt0:         mask for m>0 among cosine entries, shape (n_lm,)
        wnorm:        normylm flat array, shape ((lmax+1)^2,)
    """
    # Total number of (l, m) pairs with m >= 0
    n_lm = (lmax + 1) * (lmax + 2) // 2

    # Build l and m for every (l,m) pair: l repeats as [0, 1,1, 2,2,2, ...]
    l_vals = np.repeat(np.arange(lmax + 1), np.arange(1, lmax + 2))
    m_vals = np.concatenate([np.arange(l + 1) for l in range(lmax + 1)])

    # Flat index for the cosine term of each (l,m):
    # For degree l, the block starts at l^2. Within the block, m=0 is at offset 0,
    # m>0 cos is at offset 2*m-1, m>0 sin is at offset 2*m.
    block_start = l_vals ** 2
    cos_offset = np.where(m_vals == 0, 0, 2 * m_vals - 1)
    flat_cos_idx = block_start + cos_offset

    # Flat index for sine terms (m > 0 only)
    mgt0 = m_vals > 0
    sin_offset = 2 * m_vals[mgt0]
    flat_sin_idx = block_start[mgt0] + sin_offset

    # normylm: 1.0 for m=0, sqrt(2) for m!=0 (both cos and sin slots)
    n_flat = (lmax + 1) ** 2
    wnorm = np.ones(n_flat)
    wnorm[flat_cos_idx[mgt0]] = np.sqrt(2.0)
    wnorm[flat_sin_idx] = np.sqrt(2.0)

    return {
        "flat_cos_idx": flat_cos_idx,
        "flat_sin_idx": flat_sin_idx,
        "l_vals": l_vals,
        "m_vals": m_vals,
        "mgt0": mgt0,
        "wnorm": wnorm,
    }


def _normylm_flat(lmax: int) -> np.ndarray:
    """Normylm weights as a flat (lmax+1)^2 array: 1.0 for m=0, sqrt(2) for m!=0."""
    return _index_tables(lmax)["wnorm"]


def fortran_flat_to_shcoeffs(flat: np.ndarray, lmax: int) -> np.ndarray:
    """Convert Fortran flat (lmax+1)^2 array to cilm[2, lmax+1, lmax+1]."""
    t = _index_tables(lmax)
    cilm = np.zeros((2, lmax + 1, lmax + 1))
    cilm[0, t["l_vals"], t["m_vals"]] = flat[t["flat_cos_idx"]]
    cilm[1, t["l_vals"][t["mgt0"]], t["m_vals"][t["mgt0"]]] = flat[t["flat_sin_idx"]]
    return cilm


def shcoeffs_to_fortran_flat(cilm: np.ndarray) -> np.ndarray:
    """Convert pyshtools cilm[2, lmax+1, lmax+1] to Fortran flat array."""
    lmax = cilm.shape[1] - 1
    t = _index_tables(lmax)
    flat = np.zeros((lmax + 1) ** 2)
    flat[t["flat_cos_idx"]] = cilm[0, t["l_vals"], t["m_vals"]]
    flat[t["flat_sin_idx"]] = cilm[1, t["l_vals"][t["mgt0"]], t["m_vals"][t["mgt0"]]]
    return flat


def fortran_flat_raw_to_shcoeffs(flat: np.ndarray, lmax: int) -> np.ndarray:
    """Convert .raw coefficients (with 0.01 and normylm) to pyshtools cilm."""
    cilm = fortran_flat_to_shcoeffs(flat, lmax)
    cilm *= 100.0
    cilm[:, :, 1:] /= np.sqrt(2.0)
    return cilm


def shcoeffs_to_fortran_flat_raw(cilm: np.ndarray) -> np.ndarray:
    """Convert pyshtools cilm to .raw format (with normylm and 0.01)."""
    cilm_scaled = cilm.copy()
    cilm_scaled[:, :, 1:] *= np.sqrt(2.0)
    flat = shcoeffs_to_fortran_flat(cilm_scaled)
    flat *= 0.01
    return flat


def batch_fortran_flat_raw_to_shcoeffs(flat_batch: np.ndarray, lmax: int) -> np.ndarray:
    """Convert a batch of .raw flat arrays to a stack of cilm arrays.

    Vectorized equivalent of calling fortran_flat_raw_to_shcoeffs row-by-row.

    Args:
        flat_batch: shape (nlayers, natd) in .raw convention (0.01 scale + normylm).
        lmax: Maximum SH degree.

    Returns:
        shape (nlayers, 2, lmax+1, lmax+1).
    """
    t = _index_tables(lmax)
    nlayers = flat_batch.shape[0]
    cilm = np.zeros((nlayers, 2, lmax + 1, lmax + 1), dtype=np.float64)
    cilm[:, 0, t["l_vals"], t["m_vals"]] = flat_batch[:, t["flat_cos_idx"]]
    cilm[:, 1, t["l_vals"][t["mgt0"]], t["m_vals"][t["mgt0"]]] = flat_batch[:, t["flat_sin_idx"]]
    cilm *= 100.0
    cilm[:, :, :, 1:] /= np.sqrt(2.0)
    return cilm


def cilm_stack_to_internal(cilm_stack: np.ndarray) -> np.ndarray:
    """Convert a stack of cilm arrays to internal flat format.

    Args:
        cilm_stack: shape (ndepths, 2, lmax+1, lmax+1).

    Returns:
        shape (ndepths, natd) in Fortran flat raw convention.
    """
    ndepths = cilm_stack.shape[0]
    natd = (cilm_stack.shape[2]) ** 2
    result = np.empty((ndepths, natd), dtype=np.float64)
    for i in range(ndepths):
        result[i] = shcoeffs_to_fortran_flat_raw(cilm_stack[i])
    return result


def internal_to_cilm_stack(flat_stack: np.ndarray, lmax: int) -> np.ndarray:
    """Convert internal flat format to a stack of cilm arrays.

    Args:
        flat_stack: shape (ndepths, natd) in Fortran flat raw convention.
        lmax: Maximum spherical harmonic degree.

    Returns:
        shape (ndepths, 2, lmax+1, lmax+1).
    """
    ndepths = flat_stack.shape[0]
    result = np.empty((ndepths, 2, lmax + 1, lmax + 1), dtype=np.float64)
    for i in range(ndepths):
        result[i] = fortran_flat_raw_to_shcoeffs(flat_stack[i], lmax)
    return result
