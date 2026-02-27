"""Tests for the class-based tomographic filtering API."""

import numpy as np
import pytest

from srts.coeffs import (
    cilm_stack_to_internal,
    fortran_flat_raw_to_shcoeffs,
    internal_to_cilm_stack,
    shcoeffs_to_fortran_flat_raw,
)
from srts.expansion import (
    SphericalHarmonicExpansion,
    expand_with_precomputed,
    precompute_expansion,
)
from srts.filtering import apply_resolution_matrix, extract_sp_from_spt
from srts.model_data import load_model_data
from srts.parameterization import DepthParameterization, _spline_projection_operator, project_layer_to_sph
from srts.splines import get_spline_basis
from srts.tomographic_filter import S12RTS, S20RTS, S40RTS, TomographicFilter


# ---------------------------------------------------------------------------
# Round-trip: cilm <-> internal format
# ---------------------------------------------------------------------------

def _random_valid_cilm(rng, ndepths, lmax):
    """Generate random cilm arrays with valid structure (m <= l only)."""
    cilm = np.zeros((ndepths, 2, lmax + 1, lmax + 1))
    for l in range(lmax + 1):
        for m in range(l + 1):
            cilm[:, 0, l, m] = rng.standard_normal(ndepths)
            if m > 0:
                cilm[:, 1, l, m] = rng.standard_normal(ndepths)
    return cilm


class TestCilmRoundTrip:
    """Verify cilm_stack_to_internal and internal_to_cilm_stack are inverses."""

    @pytest.mark.parametrize("lmax", [4, 12, 20])
    def test_round_trip_random(self, lmax):
        rng = np.random.default_rng(42)
        ndepths = 5
        cilm = _random_valid_cilm(rng, ndepths, lmax)
        flat = cilm_stack_to_internal(cilm)
        recovered = internal_to_cilm_stack(flat, lmax)
        np.testing.assert_allclose(recovered, cilm, atol=1e-14)

    @pytest.mark.parametrize("lmax", [4, 12])
    def test_single_coefficients(self, lmax):
        """A single cilm round-trips through the batch functions."""
        cilm = np.zeros((1, 2, lmax + 1, lmax + 1))
        cilm[0, 0, 2, 1] = 3.5
        cilm[0, 1, 3, 2] = -1.2
        flat = cilm_stack_to_internal(cilm)
        recovered = internal_to_cilm_stack(flat, lmax)
        np.testing.assert_allclose(recovered, cilm, atol=1e-14)


# ---------------------------------------------------------------------------
# SphericalHarmonicExpansion class vs functional API
# ---------------------------------------------------------------------------

class TestSphericalHarmonicExpansion:

    @pytest.fixture
    def grid_data(self):
        """Simple lon/lat grid with synthetic data."""
        lons = np.linspace(0, 360, 73, endpoint=False)
        lats = np.linspace(-90, 90, 37)
        lon_grid, lat_grid = np.meshgrid(lons, lats)
        lon = lon_grid.ravel()
        lat = lat_grid.ravel()
        # Simple harmonic pattern
        values = np.cos(np.radians(lat)) * np.sin(2 * np.radians(lon))
        return lon, lat, values

    def test_expand_matches_functional(self, grid_data):
        """Class expand() output, converted to raw, matches expand_to_sh."""
        lon, lat, values = grid_data
        lmax = 8

        # Functional
        precomp = precompute_expansion(lon, lat, lmax)
        raw_func = expand_with_precomputed(precomp, values)

        # Class-based
        expander = SphericalHarmonicExpansion(lon, lat, lmax)
        cilm = expander.expand(values)

        # Convert cilm back to raw to compare
        raw_class = shcoeffs_to_fortran_flat_raw(cilm)
        np.testing.assert_allclose(raw_class, raw_func, atol=1e-14)

    def test_expand_batch(self, grid_data):
        """expand_batch produces same results as repeated expand calls."""
        lon, lat, values = grid_data
        lmax = 8

        # Create multiple "layers" with different patterns
        rng = np.random.default_rng(123)
        nlayers = 4
        batch = np.stack([values * (i + 1) + rng.standard_normal(len(values)) * 0.01
                          for i in range(nlayers)])

        expander = SphericalHarmonicExpansion(lon, lat, lmax)
        batch_result = expander.expand_batch(batch)

        assert batch_result.shape == (nlayers, 2, lmax + 1, lmax + 1)
        for i in range(nlayers):
            single = expander.expand(batch[i])
            np.testing.assert_allclose(batch_result[i], single, atol=1e-14)

    def test_lmax_property(self, grid_data):
        lon, lat, _ = grid_data
        expander = SphericalHarmonicExpansion(lon, lat, lmax=12)
        assert expander.lmax == 12


# ---------------------------------------------------------------------------
# DepthParameterization class vs functional API
# ---------------------------------------------------------------------------

class TestDepthParameterization:

    def test_project_layer_matches_functional(self):
        """project_layer output, converted to internal, matches project_layer_to_sph."""
        lmax = 8
        rng = np.random.default_rng(99)
        cilm = rng.standard_normal((2, lmax + 1, lmax + 1)) * 0.001
        dep_top, dep_bottom = 500.0, 700.0

        # Functional path
        raw = shcoeffs_to_fortran_flat_raw(cilm)
        spline_basis = get_spline_basis()
        spline_V, spline_lam, spline_z = _spline_projection_operator(spline_basis)
        sph_func = project_layer_to_sph(raw, lmax, dep_top, dep_bottom, spline_V, spline_lam, spline_z)

        # Class-based path
        projector = DepthParameterization()
        sph_class = projector.project_layer(cilm, dep_top, dep_bottom)

        # Convert class result to internal for comparison
        sph_class_flat = cilm_stack_to_internal(sph_class)
        np.testing.assert_allclose(sph_class_flat, sph_func, atol=1e-14)

    def test_reparameterize_matches_functional(self):
        """reparameterize sums match manual project_layer calls."""
        lmax = 4
        rng = np.random.default_rng(77)
        nlayers = 3
        depth_boundaries = np.array([200.0, 500.0, 800.0, 1100.0])
        layer_cilms = [rng.standard_normal((2, lmax + 1, lmax + 1)) * 0.001
                       for _ in range(nlayers)]

        projector = DepthParameterization()
        result = projector.reparameterize(layer_cilms, depth_boundaries)

        # Manual sum
        manual = np.zeros_like(result)
        for i in range(nlayers):
            manual += projector.project_layer(
                layer_cilms[i], depth_boundaries[i], depth_boundaries[i + 1]
            )

        np.testing.assert_allclose(result, manual, atol=1e-14)

    def test_evaluate_at_depth(self):
        """evaluate_at_depth returns the correct cilm at a given depth."""
        lmax = 4
        rng = np.random.default_rng(55)
        sph_cilm = rng.standard_normal((21, 2, lmax + 1, lmax + 1)) * 0.001

        depth = 1000.0
        result = DepthParameterization.evaluate_at_depth(sph_cilm, depth)

        # Compare with functional path
        from srts.analysis import evaluate_at_depth as eval_depth_func
        flat = cilm_stack_to_internal(sph_cilm)
        raw_func = eval_depth_func(flat, lmax, depth)
        cilm_func = fortran_flat_raw_to_shcoeffs(raw_func, lmax)

        np.testing.assert_allclose(result, cilm_func, atol=1e-13)

    def test_evaluate_at_depths(self):
        """evaluate_at_depths matches repeated evaluate_at_depth calls."""
        lmax = 4
        rng = np.random.default_rng(33)
        sph_cilm = rng.standard_normal((21, 2, lmax + 1, lmax + 1)) * 0.001
        depths = np.array([500.0, 1000.0, 1500.0, 2000.0])

        result = DepthParameterization.evaluate_at_depths(sph_cilm, depths)
        assert result.shape == (4, 2, lmax + 1, lmax + 1)

        for i, d in enumerate(depths):
            single = DepthParameterization.evaluate_at_depth(sph_cilm, d)
            np.testing.assert_allclose(result[i], single, atol=1e-13)


# ---------------------------------------------------------------------------
# TomographicFilter class
# ---------------------------------------------------------------------------

class TestTomographicFilter:

    def test_construction_and_properties(self):
        filt = TomographicFilter(40)
        assert filt.degree == 40
        assert filt.lmax == 40
        assert filt.eps == 20e-4

    def test_custom_eps(self):
        filt = TomographicFilter(20, eps=0.01)
        assert filt.eps == 0.01

    def test_invalid_degree(self):
        with pytest.raises(ValueError, match="degree must be 12, 20, or 40"):
            TomographicFilter(30)

    def test_factory_functions(self):
        s40 = S40RTS()
        assert s40.degree == 40
        assert s40.eps == 20e-4

        s20 = S20RTS()
        assert s20.degree == 20
        assert s20.eps == 35e-4

        s12 = S12RTS()
        assert s12.degree == 12
        assert s12.eps == 40e-4

    def test_reference_model_shape(self):
        filt = S40RTS()
        ref = filt.reference_model
        assert ref.ndim == 4
        assert ref.shape[1] == 2
        assert ref.shape[2] == 41
        assert ref.shape[3] == 41

    def test_filter_matches_functional(self):
        """Class filter() produces identical results to the functional API."""
        md = load_model_data(40)

        # Use the reference model coefficients as input (internal format)
        # Pad to 21 depth levels if needed
        repar_internal = np.zeros((21, md.natd), dtype=np.float64)
        ndp = md.reference_coefficients.shape[0]
        repar_internal[:ndp] = md.reference_coefficients

        # Functional path
        filtered_func = apply_resolution_matrix(repar_internal, md, eps=20e-4)
        filt_func, _, _ = extract_sp_from_spt(filtered_func, md)

        # Class path
        repar_cilm = internal_to_cilm_stack(repar_internal, md.lmax)
        filt = S40RTS()
        filtered_class = filt.filter(repar_cilm)

        # Convert class result back to internal for comparison
        filtered_class_flat = cilm_stack_to_internal(filtered_class)

        np.testing.assert_allclose(filtered_class_flat, filt_func, atol=1e-12)


# ---------------------------------------------------------------------------
# TomographicFilter against Fortran reference
# ---------------------------------------------------------------------------

class TestTomographicFilterFortranReference:
    """Test class-based API against Fortran reference outputs."""

    def test_filter_fortran_repar(self, fortran_reference_dir):
        from srts.io import read_sph

        repar_path = fortran_reference_dir / "inpm.S40.examplefile.dvs.repar.sph"
        filt_path = fortran_reference_dir / "oupm.S40.examplefile.dvs.filt.sph"
        if not repar_path.exists() or not filt_path.exists():
            pytest.skip("Fortran reference data not available")

        repar = read_sph(repar_path)
        filt_ref = read_sph(filt_path)

        # Convert Fortran reparameterized to cilm
        repar_cilm = internal_to_cilm_stack(repar["coefficients"], repar["lmax"])

        # Pad to 21 depths (Fortran repar has ndmn=4, ndmx=24 → 21 depths)
        assert repar_cilm.shape[0] == 21

        # Filter using class API
        s40 = S40RTS()
        filtered_cilm = s40.filter(repar_cilm)
        filtered_flat = cilm_stack_to_internal(filtered_cilm)

        filt_ref_coeffs = filt_ref["coefficients"]
        max_abs_diff = np.max(np.abs(filtered_flat - filt_ref_coeffs))
        max_val = np.max(np.abs(filt_ref_coeffs))
        assert max_abs_diff < 0.001 * max_val, (
            f"Max abs diff={max_abs_diff:.6e}, max value={max_val:.6e}"
        )


# ---------------------------------------------------------------------------
# Imports from package root
# ---------------------------------------------------------------------------

class TestPackageExports:

    def test_imports(self):
        from srts import (
            DepthParameterization,
            S12RTS,
            S20RTS,
            S40RTS,
            SphericalHarmonicExpansion,
            TomographicFilter,
            tomographic_filter,
        )
        assert callable(tomographic_filter)
        assert callable(S40RTS)
        assert callable(S20RTS)
        assert callable(S12RTS)
