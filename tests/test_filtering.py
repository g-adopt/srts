"""Tests for resolution matrix filtering (filtering.py)."""

import numpy as np
import pytest

from srts.filtering import apply_resolution_matrix, extract_sp_from_spt, DEFAULT_EPS
from srts.model_data import load_model_data


class TestDefaultEps:

    def test_values_exist(self):
        assert 12 in DEFAULT_EPS
        assert 20 in DEFAULT_EPS
        assert 40 in DEFAULT_EPS

    def test_s40_eps(self):
        assert DEFAULT_EPS[40] == pytest.approx(20e-4)


class TestApplyResolutionMatrix:

    @pytest.fixture
    def s40_data(self):
        return load_model_data(40)

    def test_output_shape(self, s40_data):
        natd = s40_data.natd
        ndp = s40_data.ndep
        model = np.random.default_rng(42).standard_normal((ndp, natd)) * 1e-4
        result = apply_resolution_matrix(model, s40_data, eps=20e-4)
        assert result.shape == (ndp, natd)

    def test_zero_input_gives_zero_output(self, s40_data):
        model = np.zeros((s40_data.ndep, s40_data.natd))
        result = apply_resolution_matrix(model, s40_data, eps=20e-4)
        np.testing.assert_allclose(result, 0.0, atol=1e-15)

    def test_linearity(self, s40_data):
        """R(alpha*m) = alpha * R(m) when icrust=0 or for non-first-layer terms."""
        rng = np.random.default_rng(42)
        model = rng.standard_normal((s40_data.ndep, s40_data.natd)) * 1e-4
        r1 = apply_resolution_matrix(model, s40_data, eps=20e-4)
        r2 = apply_resolution_matrix(2.0 * model, s40_data, eps=20e-4)
        # If icrust==1, the first depth layer gets ×1000, but that's also linear
        np.testing.assert_allclose(r2, 2.0 * r1, rtol=1e-12)

    def test_eigenvalue_cutoff(self, s40_data):
        """Verify eigenvalue cutoff computation is sensible."""
        eps = 20e-4
        eta = s40_data.eigenvalues[0] * eps
        cutoff = eta / 5000.0
        n_active = np.sum(s40_data.eigenvalues > cutoff)
        assert n_active > 0
        # For S40RTS, all eigenvectors may be active (cutoff is very low)
        assert n_active <= len(s40_data.eigenvalues)


class TestExtractSp:

    def test_extract_shape(self):
        md = load_model_data(40)
        filtered = np.random.default_rng(42).standard_normal((md.ndep, md.natd))
        coeffs, ndmn, ndmx = extract_sp_from_spt(filtered, md)
        ndp = ndmx - ndmn + 1
        assert coeffs.shape == (ndp, md.natd)

    def test_extract_values(self):
        md = load_model_data(40)
        ndp = int(md.ipardps[1, 0]) - int(md.ipardps[0, 0]) + 1
        filtered = np.arange(md.ndep * md.natd, dtype=np.float64).reshape(md.ndep, md.natd)
        coeffs, ndmn, ndmx = extract_sp_from_spt(filtered, md)
        np.testing.assert_array_equal(coeffs, filtered[:ndp])
