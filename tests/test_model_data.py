"""Tests for model data loading (model_data.py)."""

import numpy as np
import pytest

from srts.model_data import ModelData, load_model_data


class TestLoadModelData:

    @pytest.mark.parametrize("degree", [12, 20, 40])
    def test_load_succeeds(self, degree):
        md = load_model_data(degree)
        assert md.lmax == degree
        assert md.natd == (degree + 1) ** 2

    def test_invalid_degree(self):
        with pytest.raises(ValueError, match="degree must be"):
            load_model_data(30)

    def test_s40rts_metadata(self):
        md = load_model_data(40)
        assert md.lmax == 40
        assert md.ndep == 21
        assert md.natd == 1681
        assert md.lenatd == 21 * 1681

    def test_eigenvalues_positive_decreasing(self):
        md = load_model_data(40)
        assert np.all(md.eigenvalues > 0)
        assert md.eigenvalues[0] > md.eigenvalues[-1]
        # They should be sorted in decreasing order (as read from Fortran)
        assert np.all(np.diff(md.eigenvalues) <= 0)

    def test_eigenvectors_shape(self):
        md = load_model_data(40)
        neig = len(md.eigenvalues)
        assert md.eigenvectors.shape == (neig, md.lenatd)

    def test_weights_shape(self):
        md = load_model_data(40)
        assert md.weights.shape[1] == md.natd

    def test_reference_model(self):
        md = load_model_data(40)
        assert md.reference_lmax == 40
        assert md.reference_ndmn == 4
        assert md.reference_ndmx == 24
        ndp = md.reference_ndmx - md.reference_ndmn + 1
        assert md.reference_coefficients.shape == (ndp, md.natd)

    def test_iparsw_first_active(self):
        """First parameter (SP) should always be active."""
        for degree in [12, 20, 40]:
            md = load_model_data(degree)
            assert md.iparsw[0] == 1
