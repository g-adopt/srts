"""Tests for the 21-knot cardinal cubic spline basis (splines.py)."""

import numpy as np
import pytest

from srts.splines import (
    NKNOTS,
    RCMB,
    REARTH,
    RMOHO,
    SPLINE_KNOTS,
    SplineBasis,
    depth_to_xd,
    get_spline_basis,
    xd_to_depth,
)


class TestDepthConversion:

    def test_cmb(self):
        depth_cmb = REARTH - RCMB
        assert depth_to_xd(depth_cmb) == pytest.approx(-1.0)

    def test_moho(self):
        depth_moho = REARTH - RMOHO
        assert depth_to_xd(depth_moho) == pytest.approx(1.0)

    def test_round_trip_scalar(self):
        for depth in [100.0, 500.0, 1500.0, 2500.0]:
            xd = depth_to_xd(depth)
            recovered = xd_to_depth(xd)
            assert recovered == pytest.approx(depth, abs=1e-10)

    def test_round_trip_array(self):
        depths = np.array([100.0, 500.0, 1500.0, 2500.0])
        xd = depth_to_xd(depths)
        recovered = xd_to_depth(xd)
        np.testing.assert_allclose(recovered, depths, atol=1e-10)

    def test_monotonicity(self):
        """Deeper = more negative xd."""
        depths = np.array([100.0, 500.0, 1500.0, 2500.0])
        xd = depth_to_xd(depths)
        assert np.all(np.diff(xd) < 0)


class TestSplineBasis:

    def test_identity_at_knots(self):
        """Basis i should equal 1 at knot (20-i) and 0 at all other knots.

        This anti-diagonal pattern matches the Fortran splh(ind) convention:
        splh(0) peaks at the shallowest knot (xd=+1), splh(20) at the deepest (xd=-1).
        """
        basis = SplineBasis()
        vals = basis.evaluate_all(SPLINE_KNOTS)  # (21, 21)
        expected = np.eye(NKNOTS)[::-1, :]
        np.testing.assert_allclose(vals, expected, atol=1e-12)

    def test_partition_of_unity(self):
        """Spline basis functions should sum to ~1 at interior points (cardinal splines)."""
        basis = SplineBasis()
        xd = np.linspace(-0.99, 0.99, 100)
        vals = basis.evaluate_all(xd)  # (21, 100)
        sums = vals.sum(axis=0)
        np.testing.assert_allclose(sums, 1.0, atol=1e-10)

    def test_outside_range(self):
        """Points outside [-1, 1] should return 0."""
        basis = SplineBasis()
        xd = np.array([-1.5, 1.5, -2.0, 3.0])
        vals = basis.evaluate_all(xd)
        np.testing.assert_array_equal(vals, 0.0)

    def test_shape(self):
        basis = SplineBasis()
        xd = np.linspace(-1, 1, 50)
        vals = basis.evaluate_all(xd)
        assert vals.shape == (21, 50)

    def test_singleton_returns_correct(self):
        """get_spline_basis should return the cached singleton."""
        b1 = get_spline_basis()
        b2 = get_spline_basis()
        assert b1 is b2

    def test_evaluate_single_basis(self):
        basis = SplineBasis()
        xd = np.array([SPLINE_KNOTS[10]])
        val = basis.evaluate(10, xd)
        assert val[0] == pytest.approx(1.0, abs=1e-12)
        val_other = basis.evaluate(5, xd)
        assert val_other[0] == pytest.approx(0.0, abs=1e-12)

    def test_smoothness(self):
        """Check that basis functions produce smooth (non-oscillatory) values at fine resolution."""
        basis = SplineBasis()
        xd = np.linspace(-1, 1, 1000)
        vals = basis.evaluate_all(xd)  # (21, 1000)
        for i in range(NKNOTS):
            v = vals[i]
            assert v.min() >= -0.25, f"Basis {i} has unexpected negative value {v.min()}"
            assert v.max() <= 1.05, f"Basis {i} has unexpected large value {v.max()}"
