"""Tests for analysis tools (analysis.py)."""

import numpy as np
import pytest

from srts.analysis import (
    correlation,
    correlation_batch,
    evaluate_at_depth,
    evaluate_at_depths,
    power_spectrum,
    power_spectrum_batch,
)


class TestPowerSpectrum:

    def test_zero_input(self):
        lmax = 4
        coeffs = np.zeros((lmax + 1) ** 2)
        pd, total = power_spectrum(coeffs, lmax)
        assert total == 0.0
        np.testing.assert_array_equal(pd, 0.0)

    def test_single_degree(self):
        """Put all power in degree 2 and verify."""
        lmax = 4
        natd = (lmax + 1) ** 2
        coeffs = np.zeros(natd)
        # degree 2 block: indices 4,5,6,7,8 (5 terms)
        coeffs[4:9] = 1.0
        pd, total = power_spectrum(coeffs, lmax)
        assert pd[2] == pytest.approx(1.0)  # sqrt(5/5) = 1.0
        for l in range(lmax + 1):
            if l != 2:
                assert pd[l] == pytest.approx(0.0, abs=1e-15)

    def test_positive(self):
        lmax = 10
        rng = np.random.default_rng(42)
        coeffs = rng.standard_normal((lmax + 1) ** 2)
        pd, total = power_spectrum(coeffs, lmax)
        assert total > 0
        assert np.all(pd >= 0)

    def test_batch_matches_single(self):
        lmax = 6
        rng = np.random.default_rng(42)
        batch = rng.standard_normal((5, (lmax + 1) ** 2))
        pd_batch, total_batch = power_spectrum_batch(batch, lmax)
        for i in range(5):
            pd_single, total_single = power_spectrum(batch[i], lmax)
            np.testing.assert_allclose(pd_batch[i], pd_single, atol=1e-14)
            assert total_batch[i] == pytest.approx(total_single, abs=1e-14)


class TestCorrelation:

    def test_self_correlation(self):
        lmax = 6
        rng = np.random.default_rng(42)
        coeffs = rng.standard_normal((lmax + 1) ** 2)
        pd, total = correlation(coeffs, coeffs, lmax)
        assert total == pytest.approx(1.0, abs=1e-14)
        np.testing.assert_allclose(pd[1:], 1.0, atol=1e-14)

    def test_anticorrelation(self):
        lmax = 6
        rng = np.random.default_rng(42)
        coeffs = rng.standard_normal((lmax + 1) ** 2)
        pd, total = correlation(coeffs, -coeffs, lmax)
        assert total == pytest.approx(-1.0, abs=1e-14)

    def test_uncorrelated(self):
        """Two large random vectors should have low total correlation."""
        lmax = 40
        rng = np.random.default_rng(42)
        c1 = rng.standard_normal((lmax + 1) ** 2)
        c2 = rng.standard_normal((lmax + 1) ** 2)
        _, total = correlation(c1, c2, lmax)
        assert abs(total) < 0.3

    def test_batch_matches_single(self):
        lmax = 6
        rng = np.random.default_rng(42)
        b1 = rng.standard_normal((5, (lmax + 1) ** 2))
        b2 = rng.standard_normal((5, (lmax + 1) ** 2))
        pd_batch, total_batch = correlation_batch(b1, b2, lmax)
        for i in range(5):
            pd_s, total_s = correlation(b1[i], b2[i], lmax)
            np.testing.assert_allclose(pd_batch[i], pd_s, atol=1e-14)
            assert total_batch[i] == pytest.approx(total_s, abs=1e-14)


class TestDepthEvaluation:

    def test_single_depth(self):
        """evaluate_at_depth returns correct shape."""
        lmax = 4
        natd = (lmax + 1) ** 2
        sph = np.random.default_rng(42).standard_normal((21, natd))
        result = evaluate_at_depth(sph, lmax, 1000.0)
        assert result.shape == (natd,)

    def test_multiple_depths(self):
        lmax = 4
        natd = (lmax + 1) ** 2
        sph = np.random.default_rng(42).standard_normal((21, natd))
        depths = np.array([500.0, 1000.0, 1500.0, 2000.0])
        result = evaluate_at_depths(sph, lmax, depths)
        assert result.shape == (4, natd)

    def test_consistency(self):
        """evaluate_at_depths should match evaluate_at_depth for each depth."""
        lmax = 4
        natd = (lmax + 1) ** 2
        sph = np.random.default_rng(42).standard_normal((21, natd))
        depths = np.array([500.0, 1000.0, 2000.0])
        batch = evaluate_at_depths(sph, lmax, depths)
        for i, d in enumerate(depths):
            single = evaluate_at_depth(sph, lmax, d)
            np.testing.assert_allclose(batch[i], single, atol=1e-14)
