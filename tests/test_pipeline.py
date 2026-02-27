"""Integration tests: validate against Fortran reference outputs."""

import pathlib

import numpy as np
import pytest

from srts.io import read_sph
from srts.filtering import apply_resolution_matrix, extract_sp_from_spt
from srts.model_data import load_model_data
from srts.analysis import (
    evaluate_at_depth,
    power_spectrum,
    correlation,
)


FORTRAN_REPAR = "inpm.S40.examplefile.dvs.repar.sph"
FORTRAN_FILT = "oupm.S40.examplefile.dvs.filt.sph"


def _load_reference(ref_dir, filename):
    filepath = ref_dir / filename
    if not filepath.exists():
        pytest.skip(f"Reference file {filename} not available")
    return read_sph(filepath)


class TestFilteringAgainstFortran:
    """Apply the resolution matrix to the Fortran reparameterized model and
    compare the result against the Fortran filtered output."""

    @pytest.fixture
    def reference_data(self, fortran_reference_dir):
        repar = _load_reference(fortran_reference_dir, FORTRAN_REPAR)
        filt = _load_reference(fortran_reference_dir, FORTRAN_FILT)
        md = load_model_data(40)
        return repar, filt, md

    def test_filtering_matches_fortran(self, reference_data):
        repar, filt_ref, md = reference_data
        repar_coeffs = repar["coefficients"]
        filt_ref_coeffs = filt_ref["coefficients"]

        filtered = apply_resolution_matrix(repar_coeffs, md, eps=20e-4)
        filt_coeffs, ndmn, ndmx = extract_sp_from_spt(filtered, md)

        assert ndmn == filt_ref["ndmn"]
        assert ndmx == filt_ref["ndmx"]
        assert filt_coeffs.shape == filt_ref_coeffs.shape

        # The Fortran .sph file has e12.4 precision (~4 decimal digits), so
        # comparing filtered output (computed from those truncated inputs) to
        # the stored Fortran filtered output should match within that precision.
        max_abs_diff = np.max(np.abs(filt_coeffs - filt_ref_coeffs))
        max_val = np.max(np.abs(filt_ref_coeffs))
        assert max_abs_diff < 0.01 * max_val, (
            f"Max abs diff={max_abs_diff:.6e}, max value={max_val:.6e}"
        )

    def test_power_spectrum_at_depth(self, reference_data):
        """Power spectra of filtered models should be close at a representative depth."""
        repar, filt_ref, md = reference_data
        repar_coeffs = repar["coefficients"]
        filt_ref_coeffs = filt_ref["coefficients"]

        filtered = apply_resolution_matrix(repar_coeffs, md, eps=20e-4)
        filt_coeffs, _, _ = extract_sp_from_spt(filtered, md)

        depth = 1000.0
        raw_python = evaluate_at_depth(filt_coeffs, 40, depth)
        raw_fortran = evaluate_at_depth(filt_ref_coeffs, 40, depth)

        pd_python, total_python = power_spectrum(raw_python, 40)
        pd_fortran, total_fortran = power_spectrum(raw_fortran, 40)

        assert total_python == pytest.approx(total_fortran, rel=0.05)

    def test_correlation_with_reference(self, reference_data):
        """Filtered model should have positive correlation with the real S40RTS at 1000 km."""
        _, filt_ref, md = reference_data
        filt_ref_coeffs = filt_ref["coefficients"]

        depth = 1000.0
        raw_filt = evaluate_at_depth(filt_ref_coeffs, 40, depth)
        raw_ref = evaluate_at_depth(md.reference_coefficients, 40, depth)

        _, total = correlation(raw_filt, raw_ref, 40)
        assert total > 0, f"Expected positive correlation, got {total}"


class TestSphRoundTrip:
    """Verify that reading Fortran .sph, writing it back, and re-reading gives the same data."""

    def test_repar_round_trip(self, fortran_reference_dir, tmp_path):
        from srts.io import write_sph

        ref = _load_reference(fortran_reference_dir, FORTRAN_REPAR)
        out = tmp_path / "round_trip.sph"
        write_sph(out, ref["lmax"], ref["ndmn"], ref["ndmx"], ref["coefficients"])
        reread = read_sph(out)
        np.testing.assert_allclose(
            reread["coefficients"], ref["coefficients"], atol=5e-4
        )
