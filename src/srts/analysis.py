"""Power spectra, correlations, grid expansion, depth profile analysis.

Replaces depmap.f, raw2xyz-eqdist.f, pwrspecsph.f, correlatorsph.f, sph2v_input.f.
All operations are vectorized — no Python loops over points, coefficients, or degrees.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np
import pyshtools

from srts.coeffs import fortran_flat_to_shcoeffs, fortran_flat_raw_to_shcoeffs, _index_tables
from srts.splines import depth_to_xd, get_spline_basis


# ---------------------------------------------------------------------------
# Precomputed degree-block indexing for flat coefficient arrays
# ---------------------------------------------------------------------------

@lru_cache(maxsize=8)
def _degree_block_indices(lmax: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return arrays that allow reduceat operations over degree blocks.

    Returns:
        block_starts: start index of each degree block, shape (lmax+1,)
        block_sizes:  2*l+1 for each l, shape (lmax+1,)
        degree_label: degree l for each flat coefficient, shape ((lmax+1)^2,)
    """
    l_arr = np.arange(lmax + 1)
    block_sizes = 2 * l_arr + 1
    block_starts = np.cumsum(block_sizes) - block_sizes  # [0, 1, 4, 9, ...]
    degree_label = np.repeat(l_arr, block_sizes)
    return block_starts, block_sizes, degree_label


# ---------------------------------------------------------------------------
# Depth evaluation
# ---------------------------------------------------------------------------

def evaluate_at_depth(sph_coeffs: np.ndarray, lmax: int, depth_km: float) -> np.ndarray:
    """Evaluate .sph model at a given depth → flat (lmax+1)^2 raw coefficients.

    Equivalent to depmap.f.
    """
    xd = depth_to_xd(depth_km)
    weights = get_spline_basis().evaluate_all(np.atleast_1d(xd))[:, 0]  # (21,)
    return weights @ sph_coeffs  # (natd,)


def evaluate_at_depths(sph_coeffs: np.ndarray, lmax: int, depths_km: np.ndarray) -> np.ndarray:
    """Evaluate .sph model at multiple depths at once.

    Args:
        sph_coeffs: shape (21, natd).
        depths_km: shape (ndepths,).

    Returns:
        shape (ndepths, natd).
    """
    xd = depth_to_xd(depths_km)
    weights = get_spline_basis().evaluate_all(xd)  # (21, ndepths)
    return weights.T @ sph_coeffs  # (ndepths, natd)


# ---------------------------------------------------------------------------
# Power spectrum
# ---------------------------------------------------------------------------

def power_spectrum(raw_coeffs: np.ndarray, lmax: int) -> tuple[np.ndarray, float]:
    """Per-degree and total power from flat SH coefficients.

    Matches pwrspecsph.f:
      P(l) = sqrt( sum_m(c_lm^2) / (2l+1) )
      P_total = sqrt( (1/sqrt(4pi)) * sum_{l>=1} P(l)^2 )
    """
    block_starts, block_sizes, _ = _degree_block_indices(lmax)

    sq = raw_coeffs ** 2
    per_degree_sum = np.add.reduceat(sq, block_starts)  # (lmax+1,)
    per_degree = np.sqrt(per_degree_sum / block_sizes)

    total = np.sqrt(np.sum(per_degree[1:] ** 2) / np.sqrt(4.0 * np.pi))
    return per_degree, total


def power_spectrum_batch(raw_batch: np.ndarray, lmax: int) -> tuple[np.ndarray, np.ndarray]:
    """Power spectrum for a batch of coefficient sets at once.

    Args:
        raw_batch: shape (nbatch, (lmax+1)^2).

    Returns:
        per_degree: shape (nbatch, lmax+1).
        total: shape (nbatch,).
    """
    block_starts, block_sizes, _ = _degree_block_indices(lmax)

    sq = raw_batch ** 2  # (nbatch, natd)
    # reduceat along axis=1
    per_degree_sum = np.add.reduceat(sq, block_starts, axis=1)  # (nbatch, lmax+1)
    per_degree = np.sqrt(per_degree_sum / block_sizes[np.newaxis, :])

    total = np.sqrt(np.sum(per_degree[:, 1:] ** 2, axis=1) / np.sqrt(4.0 * np.pi))
    return per_degree, total


# ---------------------------------------------------------------------------
# Correlation
# ---------------------------------------------------------------------------

def correlation(
    coeffs1: np.ndarray,
    coeffs2: np.ndarray,
    lmax: int,
) -> tuple[np.ndarray, float]:
    """Per-degree and total correlation between two flat coefficient sets.

    Matches correlatorsph.f.
    """
    block_starts, _, _ = _degree_block_indices(lmax)

    s1 = np.add.reduceat(coeffs1 ** 2, block_starts)
    s2 = np.add.reduceat(coeffs2 ** 2, block_starts)
    cross = np.add.reduceat(coeffs1 * coeffs2, block_starts)

    denom = np.sqrt(s1) * np.sqrt(s2)
    per_degree = np.where(denom > 0, cross / denom, 0.0)

    # Total (l >= 1)
    x1_total = np.sum(s1[1:])
    x2_total = np.sum(s2[1:])
    cross_total = np.sum(cross[1:])
    denom_total = np.sqrt(x1_total) * np.sqrt(x2_total)
    total = cross_total / denom_total if denom_total > 0 else 0.0

    return per_degree, total


def correlation_batch(
    batch1: np.ndarray,
    batch2: np.ndarray,
    lmax: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Correlation for a batch of coefficient pairs.

    Args:
        batch1, batch2: shape (nbatch, (lmax+1)^2).

    Returns:
        per_degree: shape (nbatch, lmax+1).
        total: shape (nbatch,).
    """
    block_starts, _, _ = _degree_block_indices(lmax)

    s1 = np.add.reduceat(batch1 ** 2, block_starts, axis=1)
    s2 = np.add.reduceat(batch2 ** 2, block_starts, axis=1)
    cross = np.add.reduceat(batch1 * batch2, block_starts, axis=1)

    denom = np.sqrt(s1) * np.sqrt(s2)
    per_degree = np.where(denom > 0, cross / denom, 0.0)

    x1_total = np.sum(s1[:, 1:], axis=1)
    x2_total = np.sum(s2[:, 1:], axis=1)
    cross_total = np.sum(cross[:, 1:], axis=1)
    denom_total = np.sqrt(x1_total) * np.sqrt(x2_total)
    total = np.where(denom_total > 0, cross_total / denom_total, 0.0)

    return per_degree, total


# ---------------------------------------------------------------------------
# Grid expansion
# ---------------------------------------------------------------------------

def expand_to_grid(
    raw_coeffs: np.ndarray,
    lmax: int,
    spacing: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Expand SH coefficients to equidistant lat/lon grid using pyshtools.

    Returns (lon_1d, lat_1d, values_2d).
    """
    cilm = fortran_flat_raw_to_shcoeffs(raw_coeffs, lmax)
    coeffs = pyshtools.SHCoeffs.from_array(cilm, normalization='ortho', csphase=-1)
    grid = coeffs.expand(grid='DH2', extend=True)
    return grid.lons(), grid.lats(), grid.data


def evaluate_at_points(
    sph_coeffs: np.ndarray,
    lmax: int,
    depth_km: float,
    lat: np.ndarray,
    lon: np.ndarray,
) -> np.ndarray:
    """Evaluate .sph model at arbitrary (lat, lon, depth) points.

    Equivalent to sph2v_input.f.
    """
    raw = evaluate_at_depth(sph_coeffs, lmax, depth_km)
    cilm = fortran_flat_to_shcoeffs(raw, lmax)
    return _evaluate_sh_at_points(cilm, lmax, lat, lon) * 100.0


def _evaluate_sh_at_points(
    cilm: np.ndarray,
    lmax: int,
    lat: np.ndarray,
    lon: np.ndarray,
) -> np.ndarray:
    """Evaluate SH expansion at arbitrary (lat, lon) points.

    Vectorized: groups by unique latitude, broadcasts longitude dependence.
    """
    npts = len(lat)
    result = np.zeros(npts)

    t = _index_tables(lmax)
    l_vals = t["l_vals"]
    m_vals = t["m_vals"]
    mgt0 = t["mgt0"]

    flat_cos = cilm[0, l_vals, m_vals]               # (n_lm,)
    flat_sin = cilm[1, l_vals[mgt0], m_vals[mgt0]]   # (n_mgt0,)
    plm_idx_cos = l_vals * (l_vals + 1) // 2 + m_vals
    plm_idx_sin = plm_idx_cos[mgt0]
    m_vals_mgt0 = m_vals[mgt0]

    unique_lats, inverse = np.unique(lat, return_inverse=True)

    for ilat in range(len(unique_lats)):
        mask = inverse == ilat
        cos_theta = np.cos(np.radians(90.0 - unique_lats[ilat]))
        plm = pyshtools.legendre.PlmBar(lmax, cos_theta, csphase=-1) / np.sqrt(4.0 * np.pi)

        phi = np.radians(lon[mask])  # (n_at_lat,)

        # Match Fortran LEGNDR convention: divide m>0 by sqrt(2)
        # LEGNDR = PlmBar(csphase=-1) / sqrt(4pi) / sqrt(2) for m>0
        plm_c = plm[plm_idx_cos]  # (n_lm,)
        plm_c = np.where(m_vals > 0, plm_c / np.sqrt(2.0), plm_c)
        m_phi_c = np.outer(m_vals, phi)  # (n_lm, n_at_lat)
        contrib_cos = (plm_c * flat_cos)[:, np.newaxis] * np.cos(m_phi_c)

        # Sine terms: sum over (l,m>0) of plm * flat_sin * sin(m*phi)
        plm_s = plm[plm_idx_sin] / np.sqrt(2.0)
        m_phi_s = np.outer(m_vals_mgt0, phi)
        contrib_sin = (plm_s * flat_sin)[:, np.newaxis] * np.sin(m_phi_s)

        result[mask] = contrib_cos.sum(axis=0) + contrib_sin.sum(axis=0)

    return result


# ---------------------------------------------------------------------------
# Full depth profile analysis
# ---------------------------------------------------------------------------

def depth_profile_analysis(
    repar_sph: np.ndarray,
    filt_sph: np.ndarray,
    reference_sph: np.ndarray,
    lmax: int,
    depth_range: tuple[float, float] = (25.0, 2875.0),
    depth_step: float = 25.0,
) -> dict:
    """Run the full Step 6 analysis at every depth, fully vectorized.

    Evaluates all depths at once via matrix multiply, then uses batch
    power spectrum and correlation functions.
    """
    depths = np.arange(depth_range[0], depth_range[1] + depth_step / 2, depth_step)

    # Evaluate all three models at all depths at once: (ndepths, natd) each
    raw_repar = evaluate_at_depths(repar_sph, lmax, depths)
    raw_filt = evaluate_at_depths(filt_sph, lmax, depths)
    raw_ref = evaluate_at_depths(reference_sph, lmax, depths)

    # Batch power spectra
    pd_repar, pt_repar = power_spectrum_batch(raw_repar, lmax)
    pd_filt, pt_filt = power_spectrum_batch(raw_filt, lmax)
    pd_ref, pt_ref = power_spectrum_batch(raw_ref, lmax)

    # Batch correlations
    cd_repar, ct_repar = correlation_batch(raw_repar, raw_ref, lmax)
    cd_filt, ct_filt = correlation_batch(raw_filt, raw_ref, lmax)

    return {
        "depths": depths,
        "power_repar": pt_repar,
        "power_filt": pt_filt,
        "power_ref": pt_ref,
        "power_deg_repar": pd_repar,
        "power_deg_filt": pd_filt,
        "power_deg_ref": pd_ref,
        "corr_repar_ref": ct_repar,
        "corr_filt_ref": ct_filt,
        "corr_deg_repar": cd_repar,
        "corr_deg_filt": cd_filt,
    }
