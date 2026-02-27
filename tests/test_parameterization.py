"""Tests for depth spline parameterization (parameterization.py)."""

import numpy as np
import pytest

from srts.parameterization import _spline_projection_operator, project_layer_to_sph
from srts.splines import get_spline_basis


class TestSplineProjectionOperator:

    def test_operator_shape(self):
        basis = get_spline_basis()
        V, lam, z = _spline_projection_operator(basis)
        assert z.shape[0] == 2892
        assert z.shape[1] == 21
        assert len(lam) == V.shape[1]
        assert V.shape[0] == 21

    def test_eigenvalues_positive(self):
        basis = get_spline_basis()
        V, lam, z = _spline_projection_operator(basis)
        assert np.all(lam > 0)


class TestProjectLayer:

    def test_output_shape(self):
        basis = get_spline_basis()
        V, lam, z = _spline_projection_operator(basis)
        lmax = 4
        natd = (lmax + 1) ** 2
        raw = np.ones(natd)
        result = project_layer_to_sph(raw, lmax, 100.0, 200.0, V, lam, z)
        assert result.shape == (21, natd)

    def test_deep_layer_activates_deep_splines(self):
        """A layer at CMB depth should primarily activate deep spline basis functions.

        Basis ordering: basis 0 = shallowest (Moho), basis 20 = deepest (CMB).
        So a deep layer should activate basis indices near 20.
        """
        basis = get_spline_basis()
        V, lam, z = _spline_projection_operator(basis)
        lmax = 2
        natd = (lmax + 1) ** 2
        raw = np.ones(natd)
        result = project_layer_to_sph(raw, lmax, 2800.0, 2891.0, V, lam, z)
        weights = result[:, 0]
        assert np.sum(np.abs(weights[16:])) > np.sum(np.abs(weights[:5]))

    def test_shallow_layer_activates_shallow_splines(self):
        """A layer near the Moho should primarily activate shallow spline basis functions.

        Basis ordering: basis 0 = shallowest (Moho), basis 20 = deepest (CMB).
        So a shallow layer should activate basis indices near 0.
        """
        basis = get_spline_basis()
        V, lam, z = _spline_projection_operator(basis)
        lmax = 2
        natd = (lmax + 1) ** 2
        raw = np.ones(natd)
        result = project_layer_to_sph(raw, lmax, 25.0, 100.0, V, lam, z)
        weights = result[:, 0]
        assert np.sum(np.abs(weights[:5])) > np.sum(np.abs(weights[16:]))

    def test_coefficient_scaling(self):
        """Result should scale linearly with raw_coeffs."""
        basis = get_spline_basis()
        V, lam, z = _spline_projection_operator(basis)
        lmax = 2
        natd = (lmax + 1) ** 2
        raw = np.random.default_rng(42).standard_normal(natd)
        r1 = project_layer_to_sph(raw, lmax, 500.0, 600.0, V, lam, z)
        r2 = project_layer_to_sph(2.0 * raw, lmax, 500.0, 600.0, V, lam, z)
        np.testing.assert_allclose(r2, 2.0 * r1, atol=1e-14)
