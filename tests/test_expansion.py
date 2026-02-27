"""Tests for spherical harmonic expansion (expansion.py)."""

import numpy as np
import pyshtools
import pytest

from srts.coeffs import fortran_flat_to_shcoeffs, shcoeffs_to_fortran_flat
from srts.expansion import expand_to_sh, precompute_expansion, expand_with_precomputed


def _make_grid(spacing=2.0):
    """Create a regular lon/lat grid."""
    lons_1d = np.arange(0, 360, spacing)
    lats_1d = np.arange(-90, 90 + spacing / 2, spacing)
    lon_grid, lat_grid = np.meshgrid(lons_1d, lats_1d)
    return lon_grid.ravel(), lat_grid.ravel()


class TestExpandKnownSignal:

    def test_y00_recovery(self):
        """A constant field should map to c_00 only."""
        lon, lat = _make_grid(5.0)
        lmax = 4
        values = np.ones_like(lon) * 3.0
        raw = expand_to_sh(lon, lat, values, lmax, damp=1.0)
        cilm = fortran_flat_to_shcoeffs(raw / 0.01, lmax)
        # After removing normylm: cilm should have c_00 dominant
        # The exact value depends on the grid normalization, but all other
        # coefficients should be very small compared to c_00.
        assert abs(cilm[0, 0, 0]) > 0
        ratio = np.max(np.abs(cilm[0, 1:, :])) / abs(cilm[0, 0, 0])
        assert ratio < 5e-3

    def test_y20_recovery(self):
        """Expand Y_20 (zonal degree 2) and check the dominant coefficient."""
        lon, lat = _make_grid(2.0)
        lmax = 10
        colat = np.radians(90.0 - lat)
        plm = np.array([pyshtools.legendre.PlmBar(2, np.cos(c))[3] for c in colat])
        raw = expand_to_sh(lon, lat, plm, lmax, damp=1.0)
        cilm = fortran_flat_to_shcoeffs(raw / 0.01, lmax)
        # c_20 should be much larger than all others (except maybe c_00 due to mean)
        c20 = abs(cilm[0, 2, 0])
        others = np.abs(cilm.copy())
        others[0, 2, 0] = 0
        others[0, 0, 0] = 0  # ignore mean
        assert c20 > 10 * np.max(others)


class TestPrecomputed:

    def test_precomputed_matches_direct(self):
        """precompute_expansion + expand_with_precomputed should give same results as expand_to_sh."""
        lon, lat = _make_grid(5.0)
        lmax = 6
        rng = np.random.default_rng(42)
        values = rng.standard_normal(len(lon))

        direct = expand_to_sh(lon, lat, values, lmax, damp=1.0)
        precomp = precompute_expansion(lon, lat, lmax, damp=1.0)
        precomp_result = expand_with_precomputed(precomp, values)
        np.testing.assert_allclose(precomp_result, direct, atol=1e-14)

    def test_precomputed_different_data(self):
        """Precomputed operator should work with different data on the same grid."""
        lon, lat = _make_grid(5.0)
        lmax = 6
        precomp = precompute_expansion(lon, lat, lmax, damp=1.0)

        rng = np.random.default_rng(10)
        v1 = rng.standard_normal(len(lon))
        v2 = rng.standard_normal(len(lon))

        r1 = expand_with_precomputed(precomp, v1)
        r2 = expand_with_precomputed(precomp, v2)
        # They should differ
        assert np.max(np.abs(r1 - r2)) > 1e-6
