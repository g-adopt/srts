"""Grid to spherical harmonic expansion via regularized least-squares.

Reimplements mkexpmatxy.f + invexpandxy.f using pyshtools PlmBar and
vectorized numpy operations. No Python loops over grid points.

The design matrix A has rows = ylm(point) * wnorm. We accumulate A^T A
and A^T d without materializing A by grouping points that share a latitude
(PlmBar depends only on colatitude) and vectorizing the longitude
dependence (cos/sin of m*phi).
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pyshtools
import scipy.linalg

from srts.coeffs import (
    _index_tables,
    batch_fortran_flat_raw_to_shcoeffs,
    fortran_flat_raw_to_shcoeffs,
)


@lru_cache(maxsize=8)
def _flat_index_structure(lmax: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precompute the PlmBar index, m-value, and is_sin flag for every flat SH slot.

    Returns:
        plm_idx: PlmBar storage index for each flat slot, shape ((lmax+1)^2,)
        m_val:   azimuthal order m for each flat slot, shape ((lmax+1)^2,)
        is_sin:  True for sine slots, shape ((lmax+1)^2,)
    """
    leny = (lmax + 1) ** 2
    plm_idx = np.empty(leny, dtype=np.intp)
    m_val = np.empty(leny, dtype=np.intp)
    is_sin = np.zeros(leny, dtype=bool)

    t = _index_tables(lmax)
    l_vals = t["l_vals"]
    m_vals = t["m_vals"]
    mgt0 = t["mgt0"]
    cos_idx = t["flat_cos_idx"]
    sin_idx = t["flat_sin_idx"]

    # Cosine slots
    plm_idx[cos_idx] = l_vals * (l_vals + 1) // 2 + m_vals
    m_val[cos_idx] = m_vals

    # Sine slots (m > 0)
    l_mgt0 = l_vals[mgt0]
    m_mgt0 = m_vals[mgt0]
    plm_idx[sin_idx] = l_mgt0 * (l_mgt0 + 1) // 2 + m_mgt0
    m_val[sin_idx] = m_mgt0
    is_sin[sin_idx] = True

    return plm_idx, m_val, is_sin


def _accumulate_ata_atd(
    lon: np.ndarray,
    lat: np.ndarray,
    values: np.ndarray,
    lmax: int,
    wnorm: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Accumulate A^T A and A^T d, fully vectorized within each latitude group.

    Groups points by unique latitude so PlmBar is called once per latitude.
    Within each group, the longitude dependence is handled via broadcasting
    with no Python loops over points.
    """
    leny = (lmax + 1) ** 2
    ata = np.zeros((leny, leny), dtype=np.float64)
    atd = np.zeros(leny, dtype=np.float64)

    plm_idx, m_val, is_sin = _flat_index_structure(lmax)

    unique_lats, inverse = np.unique(lat, return_inverse=True)

    for ilat in range(len(unique_lats)):
        mask = inverse == ilat
        lons = lon[mask]
        vals = values[mask]
        n = len(lons)

        cos_theta = np.cos(np.radians(90.0 - unique_lats[ilat]))
        plm = pyshtools.legendre.PlmBar(lmax, cos_theta, csphase=-1) / np.sqrt(4.0 * np.pi)

        phi = np.radians(lons)  # (n,)

        # m_phi: (leny,) outer (n,) → (leny, n)
        m_phi = np.outer(m_val, phi)
        trig = np.where(is_sin[:, np.newaxis], np.sin(m_phi), np.cos(m_phi))

        # Match Fortran LEGNDR convention: divide m>0 by sqrt(2)
        plm_at_slots = plm[plm_idx]
        plm_at_slots = np.where(m_val > 0, plm_at_slots / np.sqrt(2.0), plm_at_slots)
        scale = plm_at_slots * wnorm  # (leny,)
        At = scale[:, np.newaxis] * trig  # (leny, n)

        # ATA += At @ At^T, ATd += At @ vals
        ata += At @ At.T
        atd += At @ vals

    return ata, atd


def _solve_damped(
    ata: np.ndarray,
    atd: np.ndarray,
    damp: float,
) -> np.ndarray:
    """Eigendecompose ATA and solve with damping, matching invexpandxy.f."""
    eigenvalues, eigenvectors = scipy.linalg.eigh(ata)

    threshold = 1e-7 * eigenvalues[-1]
    active = eigenvalues > threshold

    V = eigenvectors[:, active]
    lam = eigenvalues[active]

    projections = V.T @ atd
    x = V @ (projections / (lam + damp))
    return x


def expand_to_sh(
    lon: np.ndarray,
    lat: np.ndarray,
    values: np.ndarray,
    lmax: int,
    damp: float = 1.0,
) -> np.ndarray:
    """Regularized least-squares expansion of grid data into SH coefficients.

    Returns flat coefficient array in .raw convention (with normylm and 0.01).
    """
    wnorm = _index_tables(lmax)["wnorm"]
    ata, atd = _accumulate_ata_atd(lon, lat, values, lmax, wnorm)
    x = _solve_damped(ata, atd, damp)
    x *= wnorm * 0.01
    return x


def precompute_expansion(
    lon: np.ndarray,
    lat: np.ndarray,
    lmax: int,
    damp: float = 1.0,
) -> dict:
    """Precompute the expansion operator for a fixed grid.

    Since all layers share the same grid, ATA and its eigendecomposition are
    identical. This returns everything needed to expand arbitrary data vectors
    on that grid without recomputing the matrix.

    Returns dict with keys: 'V', 'lam', 'wnorm', 'build_atd' (callable).
    """
    wnorm = _index_tables(lmax)["wnorm"]

    plm_idx, m_val, is_sin = _flat_index_structure(lmax)
    unique_lats, inverse = np.unique(lat, return_inverse=True)

    # Precompute per-latitude quantities
    lat_data = []
    for ilat in range(len(unique_lats)):
        mask = inverse == ilat
        cos_theta = np.cos(np.radians(90.0 - unique_lats[ilat]))
        plm = pyshtools.legendre.PlmBar(lmax, cos_theta, csphase=-1) / np.sqrt(4.0 * np.pi)
        phi = np.radians(lon[mask])
        m_phi = np.outer(m_val, phi)
        trig = np.where(is_sin[:, np.newaxis], np.sin(m_phi), np.cos(m_phi))
        plm_at_slots = plm[plm_idx]
        plm_at_slots = np.where(m_val > 0, plm_at_slots / np.sqrt(2.0), plm_at_slots)
        scale = plm_at_slots * wnorm
        At = scale[:, np.newaxis] * trig  # (leny, n_at_lat)
        lat_data.append((mask, At))

    # Build ATA once
    leny = (lmax + 1) ** 2
    ata = np.zeros((leny, leny), dtype=np.float64)
    for mask, At in lat_data:
        ata += At @ At.T

    eigenvalues, eigenvectors = scipy.linalg.eigh(ata)
    threshold = 1e-7 * eigenvalues[-1]
    active = eigenvalues > threshold
    V = eigenvectors[:, active]
    lam = eigenvalues[active]

    def build_atd(values: np.ndarray) -> np.ndarray:
        atd = np.zeros(leny, dtype=np.float64)
        for mask, At in lat_data:
            atd += At @ values[mask]
        return atd

    def build_atd_batch(values_batch: np.ndarray) -> np.ndarray:
        """Accumulate A^T D for all layers simultaneously via GEMM.

        Args:
            values_batch: shape (nlayers, npoints).

        Returns:
            shape (leny, nlayers).
        """
        nlayers = values_batch.shape[0]
        ATD = np.zeros((leny, nlayers), dtype=np.float64)
        for mask, At in lat_data:
            ATD += At @ values_batch[:, mask].T
        return ATD

    return {
        "V": V,
        "lam": lam,
        "wnorm": wnorm,
        "build_atd": build_atd,
        "build_atd_batch": build_atd_batch,
        "damp": damp,
    }


def expand_with_precomputed(precomp: dict, values: np.ndarray) -> np.ndarray:
    """Expand data values using a precomputed expansion operator.

    Returns flat coefficients in .raw convention (with normylm and 0.01).
    """
    atd = precomp["build_atd"](values)
    V, lam = precomp["V"], precomp["lam"]
    projections = V.T @ atd
    x = V @ (projections / (lam + precomp["damp"]))
    x *= precomp["wnorm"] * 0.01
    return x


class SphericalHarmonicExpansion:
    """Regularized least-squares SH expansion with cached grid operator.

    Wraps precompute_expansion + expand_with_precomputed, returning
    pyshtools cilm arrays instead of internal flat format.
    """

    def __init__(self, lon: np.ndarray, lat: np.ndarray, lmax: int, damp: float = 1.0):
        self._lmax = lmax
        self._precomp = precompute_expansion(lon, lat, lmax, damp)

    @property
    def lmax(self) -> int:
        return self._lmax

    def expand(self, values: np.ndarray) -> np.ndarray:
        """Expand grid values into spherical harmonic coefficients.

        Args:
            values: Grid values, shape (npoints,).

        Returns:
            cilm array, shape (2, lmax+1, lmax+1).
        """
        raw = expand_with_precomputed(self._precomp, values)
        return fortran_flat_raw_to_shcoeffs(raw, self._lmax)

    def expand_batch(self, values_batch: np.ndarray) -> np.ndarray:
        """Expand multiple layers of grid values at once.

        Fully batched: replaces the per-layer GEMV loop with a single GEMM
        per latitude group, then solves all layers simultaneously.

        Args:
            values_batch: Grid values, shape (nlayers, npoints).

        Returns:
            cilm arrays, shape (nlayers, 2, lmax+1, lmax+1).
        """
        precomp = self._precomp
        V, lam, damp, wnorm = (
            precomp["V"], precomp["lam"], precomp["damp"], precomp["wnorm"]
        )

        # (leny, nlayers): one GEMM per latitude group
        ATD = precomp["build_atd_batch"](values_batch)

        # Batch solve — identical algebra to expand_with_precomputed, applied column-wise
        projections = V.T @ ATD                                   # (n_active, nlayers)
        X = V @ (projections / (lam + damp)[:, np.newaxis])      # (leny, nlayers)
        X *= (wnorm * 0.01)[:, np.newaxis]

        # X.T is (nlayers, leny) — convert all layers to cilm in one vectorized pass
        return batch_fortran_flat_raw_to_shcoeffs(X.T, self._lmax)
