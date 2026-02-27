"""Tests for coefficient format conversion (coeffs.py)."""

import numpy as np
import pytest

from srts.coeffs import (
    _index_tables,
    _normylm_flat,
    fortran_flat_to_shcoeffs,
    shcoeffs_to_fortran_flat,
    fortran_flat_raw_to_shcoeffs,
    shcoeffs_to_fortran_flat_raw,
)


class TestIndexTables:

    def test_sizes(self):
        t = _index_tables(4)
        n_lm = 5 * 6 // 2  # 15
        n_flat = 25
        assert len(t["l_vals"]) == n_lm
        assert len(t["m_vals"]) == n_lm
        assert len(t["flat_cos_idx"]) == n_lm
        assert len(t["wnorm"]) == n_flat

    def test_l_m_values(self):
        t = _index_tables(3)
        expected_l = [0, 1, 1, 2, 2, 2, 3, 3, 3, 3]
        expected_m = [0, 0, 1, 0, 1, 2, 0, 1, 2, 3]
        np.testing.assert_array_equal(t["l_vals"], expected_l)
        np.testing.assert_array_equal(t["m_vals"], expected_m)

    def test_flat_indices_degree2(self):
        """For lmax=2, verify the flat array layout:
        [c_00, c_10, c_11c, c_11s, c_20, c_21c, c_21s, c_22c, c_22s]
        """
        t = _index_tables(2)
        # Cosine indices for (l,m) = (0,0),(1,0),(1,1),(2,0),(2,1),(2,2)
        np.testing.assert_array_equal(t["flat_cos_idx"], [0, 1, 2, 4, 5, 7])
        # Sine indices for m>0: (1,1),(2,1),(2,2)
        np.testing.assert_array_equal(t["flat_sin_idx"], [3, 6, 8])

    def test_normylm_values(self):
        t = _index_tables(2)
        wnorm = t["wnorm"]
        sqrt2 = np.sqrt(2.0)
        # m=0 slots: indices 0, 1, 4
        assert wnorm[0] == 1.0
        assert wnorm[1] == 1.0
        assert wnorm[4] == 1.0
        # m!=0 slots: indices 2,3,5,6,7,8
        for i in [2, 3, 5, 6, 7, 8]:
            assert wnorm[i] == pytest.approx(sqrt2)


class TestConversions:

    @pytest.mark.parametrize("lmax", [2, 4, 10, 20])
    def test_round_trip(self, lmax):
        """Convert Fortran flat → cilm → Fortran flat should be identity."""
        n = (lmax + 1) ** 2
        flat = np.random.default_rng(42).standard_normal(n)
        cilm = fortran_flat_to_shcoeffs(flat, lmax)
        recovered = shcoeffs_to_fortran_flat(cilm)
        np.testing.assert_allclose(recovered, flat, atol=1e-15)

    def test_cilm_shape(self):
        flat = np.zeros(25)
        cilm = fortran_flat_to_shcoeffs(flat, 4)
        assert cilm.shape == (2, 5, 5)

    def test_specific_values(self):
        """Set specific coefficients and verify they end up in the right cilm slots."""
        lmax = 2
        flat = np.zeros(9)
        flat[0] = 1.0   # c_00
        flat[2] = 2.0   # c_11_cos
        flat[3] = 3.0   # c_11_sin
        flat[7] = 4.0   # c_22_cos
        flat[8] = 5.0   # c_22_sin
        cilm = fortran_flat_to_shcoeffs(flat, lmax)
        assert cilm[0, 0, 0] == 1.0  # c_00
        assert cilm[0, 1, 1] == 2.0  # c_11_cos
        assert cilm[1, 1, 1] == 3.0  # c_11_sin
        assert cilm[0, 2, 2] == 4.0  # c_22_cos
        assert cilm[1, 2, 2] == 5.0  # c_22_sin


class TestRawConversions:

    @pytest.mark.parametrize("lmax", [2, 10, 20])
    def test_round_trip_raw(self, lmax):
        """Convert raw flat → cilm → raw flat should be identity."""
        n = (lmax + 1) ** 2
        flat = np.random.default_rng(123).standard_normal(n) * 0.01
        cilm = fortran_flat_raw_to_shcoeffs(flat, lmax)
        recovered = shcoeffs_to_fortran_flat_raw(cilm)
        np.testing.assert_allclose(recovered, flat, atol=1e-15)

    def test_raw_scaling(self):
        """Raw format has 0.01 scaling and normylm sqrt(2) for m!=0."""
        lmax = 1
        flat = np.array([0.01, 0.02, 0.03 * np.sqrt(2), 0.04 * np.sqrt(2)])
        cilm = fortran_flat_raw_to_shcoeffs(flat, lmax)
        assert cilm[0, 0, 0] == pytest.approx(1.0)
        assert cilm[0, 1, 0] == pytest.approx(2.0)
        assert cilm[0, 1, 1] == pytest.approx(3.0)
        assert cilm[1, 1, 1] == pytest.approx(4.0)
