"""Step-by-step validation against Fortran pipeline reference outputs.

Compares every stage of the Python srts pipeline against the outputs
produced by the original Fortran dofilt_ES_new on the example model
with S40RTS (degree=40).

Pipeline stages tested:
  Steps 1+2: Reparameterization (layer .dat -> .sph)
  Step 3:    Filtering (apply resolution matrix -> filtered .sph)
  Step 6a:   Depth evaluation (.sph -> .raw at specific depths)
  Step 6c:   Power spectrum (.raw -> per-degree and total power)
  Step 6d:   Correlation (.raw x .raw -> cross-correlation)
  Step 5:    Point evaluation (.sph -> values at arbitrary grid points)
  End-to-end: Full pipeline from input layers through analysis
"""

import pathlib

import numpy as np
import pytest

from srts.io import read_raw, read_sph
from srts.model_data import load_model_data
from srts.filtering import apply_resolution_matrix, extract_sp_from_spt
from srts.analysis import (
    evaluate_at_depth,
    evaluate_at_points,
    power_spectrum,
    correlation,
)

DEGREE = 40
FIRST_LAYER = 1
LAST_LAYER = 64
NAME = "examplefile.dvs"


def _skip_if_missing(path):
    if not path.exists():
        pytest.skip(f"Reference data not found: {path}")


# ---------------------------------------------------------------------------
# Helpers for reading Fortran output files
# ---------------------------------------------------------------------------


def _read_power_total(path):
    """Read Fortran pwr.dat: 'depth power*100' per line."""
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def _read_power_per_degree(path):
    """Read Fortran pwr.deg.dat: 'depth l power*100' per line."""
    data = np.loadtxt(path)
    depths = np.unique(data[:, 0])
    ndeg = int(data[:, 1].max()) + 1
    ndepths = len(depths)
    return depths, data[:, 2].reshape(ndepths, ndeg)


def _read_corr_total(path):
    """Read Fortran corr.dat: 'depth correlation' per line."""
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1]


def _read_corr_per_degree(path):
    """Read Fortran corr.deg.dat: 'depth l correlation' per line."""
    data = np.loadtxt(path)
    depths = np.unique(data[:, 0])
    ndeg = int(data[:, 1].max()) + 1
    ndepths = len(depths)
    return depths, data[:, 2].reshape(ndepths, ndeg)


def _read_output_points(path):
    """Read Fortran outputfile: 'lat lon value*100' per line."""
    data = np.loadtxt(path)
    return data[:, 0], data[:, 1], data[:, 2]


# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def fortran_repar(fortran_reference_dir):
    sph_path = fortran_reference_dir / "inpm.S40.examplefile.dvs.repar.sph"
    _skip_if_missing(sph_path)
    return read_sph(sph_path)


@pytest.fixture(scope="module")
def fortran_filt(fortran_reference_dir):
    sph_path = fortran_reference_dir / "oupm.S40.examplefile.dvs.filt.sph"
    _skip_if_missing(sph_path)
    return read_sph(sph_path)


@pytest.fixture(scope="module")
def model_data():
    return load_model_data(DEGREE)


# ---------------------------------------------------------------------------
# Steps 1+2: Reparameterization
# ---------------------------------------------------------------------------


class TestReparameterization:
    """Run reparameterize() on the example model and compare to Fortran .sph.

    This is the most critical test: it validates the entire SH expansion +
    spline projection chain from raw layer .dat files to the reparameterized
    .sph model.
    """

    @pytest.fixture(scope="class")
    def python_repar(self, fortran_reference_dir):
        example_dir = fortran_reference_dir / "examplemodel"
        _skip_if_missing(example_dir)
        from srts.parameterization import reparameterize

        return reparameterize(
            str(example_dir), NAME, FIRST_LAYER, LAST_LAYER, DEGREE
        )

    def test_shape_matches(self, python_repar, fortran_repar):
        assert python_repar.shape == fortran_repar["coefficients"].shape

    def test_coefficients_match(self, python_repar, fortran_repar):
        fort = fortran_repar["coefficients"]
        max_val = np.max(np.abs(fort))
        max_diff = np.max(np.abs(python_repar - fort))
        rel_diff = max_diff / max_val
        assert rel_diff < 0.01, (
            f"Reparameterization relative diff = {rel_diff:.4e} "
            f"(max_diff={max_diff:.4e}, max_val={max_val:.4e})"
        )

    def test_per_depth_correlation(self, python_repar, fortran_repar):
        """Each of the 21 depth levels should be highly correlated."""
        fort = fortran_repar["coefficients"]
        for i in range(21):
            py_row = python_repar[i]
            ft_row = fort[i]
            if np.all(np.abs(ft_row) < 1e-20):
                continue
            r = np.corrcoef(py_row, ft_row)[0, 1]
            assert r > 0.99, f"Depth level {i}: correlation = {r:.6f}"


# ---------------------------------------------------------------------------
# Step 3: Resolution matrix filtering
# ---------------------------------------------------------------------------


class TestFiltering:
    """Apply resolution matrix to Fortran reparameterized model,
    compare against Fortran filtered output."""

    def test_filtering_coefficients(self, fortran_repar, fortran_filt, model_data):
        filtered = apply_resolution_matrix(
            fortran_repar["coefficients"], model_data, eps=20e-4
        )
        filt_coeffs, ndmn, ndmx = extract_sp_from_spt(filtered, model_data)

        assert ndmn == fortran_filt["ndmn"]
        assert ndmx == fortran_filt["ndmx"]

        fort = fortran_filt["coefficients"]
        max_val = np.max(np.abs(fort))
        max_diff = np.max(np.abs(filt_coeffs - fort))
        rel_diff = max_diff / max_val
        assert rel_diff < 0.01, (
            f"Filtering relative diff = {rel_diff:.4e}"
        )

    def test_filtering_per_depth(self, fortran_repar, fortran_filt, model_data):
        """Each depth level of the filtered output should correlate highly."""
        filtered = apply_resolution_matrix(
            fortran_repar["coefficients"], model_data, eps=20e-4
        )
        filt_coeffs, _, _ = extract_sp_from_spt(filtered, model_data)
        fort = fortran_filt["coefficients"]

        for i in range(fort.shape[0]):
            if np.all(np.abs(fort[i]) < 1e-20):
                continue
            r = np.corrcoef(filt_coeffs[i], fort[i])[0, 1]
            assert r > 0.999, f"Filtered depth level {i}: correlation = {r:.6f}"


# ---------------------------------------------------------------------------
# Step 6a: Depth evaluation (.sph -> .raw)
# ---------------------------------------------------------------------------


class TestDepthEvaluation:
    """Evaluate Fortran .sph at specific depths and compare to Fortran .raw files.

    Both Python and Fortran read the same .sph (e12.4 precision) and evaluate
    via spline basis at the same depth, so agreement should be very tight.
    """

    @pytest.mark.parametrize("depth", [100, 500, 675, 1000, 1500, 2000, 2500])
    def test_repar_raw(self, fortran_repar, fortran_reference_dir, depth):
        raw_path = (
            fortran_reference_dir / "rawfiles"
            / f"inpm.S40.{NAME}.repar.{depth:04d}.raw"
        )
        _skip_if_missing(raw_path)

        _, fortran_raw = read_raw(raw_path)
        python_raw = evaluate_at_depth(
            fortran_repar["coefficients"], DEGREE, float(depth)
        )

        assert python_raw.shape == fortran_raw.shape
        max_val = np.max(np.abs(fortran_raw))
        if max_val > 1e-20:
            max_diff = np.max(np.abs(python_raw - fortran_raw))
            # Both read the same .sph file, so differences are purely numerical
            assert max_diff / max_val < 1e-4, (
                f"Depth {depth} km: rel diff = {max_diff / max_val:.4e}"
            )

    @pytest.mark.parametrize("depth", [100, 500, 675, 1000, 1500, 2000, 2500])
    def test_filt_raw(self, fortran_filt, fortran_reference_dir, depth):
        raw_path = (
            fortran_reference_dir / "rawfiles"
            / f"oupm.S40.{NAME}.filt.{depth:04d}.raw"
        )
        _skip_if_missing(raw_path)

        _, fortran_raw = read_raw(raw_path)
        python_raw = evaluate_at_depth(
            fortran_filt["coefficients"], DEGREE, float(depth)
        )

        max_val = np.max(np.abs(fortran_raw))
        if max_val > 1e-20:
            max_diff = np.max(np.abs(python_raw - fortran_raw))
            assert max_diff / max_val < 1e-4, (
                f"Depth {depth} km: rel diff = {max_diff / max_val:.4e}"
            )


# ---------------------------------------------------------------------------
# Step 6c: Power spectrum
# ---------------------------------------------------------------------------


class TestPowerSpectrum:
    """Compute power spectrum from Fortran .sph and compare to Fortran power files.

    Fortran stores power * 100 (the awk multiplication in dofilt_ES_new).
    """

    def test_total_power_repar(self, fortran_repar, fortran_reference_dir):
        pwr_path = (
            fortran_reference_dir / "pwrfiles"
            / f"inpm.S40.{NAME}.repar.pwr.dat"
        )
        _skip_if_missing(pwr_path)

        fortran_depths, fortran_power = _read_power_total(pwr_path)

        # Test every 10th depth for speed
        for i in range(0, len(fortran_depths), 10):
            depth = fortran_depths[i]
            raw = evaluate_at_depth(
                fortran_repar["coefficients"], DEGREE, depth
            )
            _, total = power_spectrum(raw, DEGREE)
            python_power = total * 100  # match Fortran *100 convention
            assert python_power == pytest.approx(fortran_power[i], rel=0.01), (
                f"Depth {depth} km: python={python_power:.6f}, "
                f"fortran={fortran_power[i]:.6f}"
            )

    def test_total_power_filt(self, fortran_filt, fortran_reference_dir):
        pwr_path = (
            fortran_reference_dir / "pwrfiles"
            / f"oupm.S40.{NAME}.filt.pwr.dat"
        )
        _skip_if_missing(pwr_path)

        fortran_depths, fortran_power = _read_power_total(pwr_path)

        for i in range(0, len(fortran_depths), 10):
            depth = fortran_depths[i]
            raw = evaluate_at_depth(
                fortran_filt["coefficients"], DEGREE, depth
            )
            _, total = power_spectrum(raw, DEGREE)
            python_power = total * 100
            assert python_power == pytest.approx(fortran_power[i], rel=0.01), (
                f"Depth {depth} km: python={python_power:.6f}, "
                f"fortran={fortran_power[i]:.6f}"
            )

    def test_per_degree_power_repar(self, fortran_repar, fortran_reference_dir):
        deg_path = (
            fortran_reference_dir / "pwrfiles"
            / f"inpm.S40.{NAME}.repar.pwr.deg.dat"
        )
        _skip_if_missing(deg_path)

        fortran_depths, fortran_per_deg = _read_power_per_degree(deg_path)

        # Test at a few representative depths
        for depth_target in [500.0, 1000.0, 2000.0]:
            idx = np.argmin(np.abs(fortran_depths - depth_target))
            depth = fortran_depths[idx]

            raw = evaluate_at_depth(
                fortran_repar["coefficients"], DEGREE, depth
            )
            per_deg, _ = power_spectrum(raw, DEGREE)
            python_per_deg = per_deg * 100

            np.testing.assert_allclose(
                python_per_deg[: fortran_per_deg.shape[1]],
                fortran_per_deg[idx],
                rtol=0.01,
                err_msg=f"Per-degree power mismatch at {depth} km",
            )


# ---------------------------------------------------------------------------
# Step 6d: Correlation with reference S40RTS
# ---------------------------------------------------------------------------


class TestCorrelation:
    """Cross-correlate Fortran .sph models with S40RTS and compare to
    Fortran correlation files."""

    def test_total_correlation_repar(self, fortran_repar, model_data, fortran_reference_dir):
        corr_path = (
            fortran_reference_dir / "comparefiles"
            / f"corr.S40RTS..inpm.S40.{NAME}.repar.corr.dat"
        )
        _skip_if_missing(corr_path)

        fortran_depths, fortran_corr = _read_corr_total(corr_path)
        ref_sph = model_data.reference_coefficients

        for i in range(0, len(fortran_depths), 10):
            depth = fortran_depths[i]
            raw_repar = evaluate_at_depth(
                fortran_repar["coefficients"], DEGREE, depth
            )
            raw_ref = evaluate_at_depth(ref_sph, DEGREE, depth)
            _, total = correlation(raw_repar, raw_ref, DEGREE)

            assert total == pytest.approx(fortran_corr[i], abs=0.01), (
                f"Depth {depth} km: python={total:.6f}, "
                f"fortran={fortran_corr[i]:.6f}"
            )

    def test_total_correlation_filt(self, fortran_filt, model_data, fortran_reference_dir):
        corr_path = (
            fortran_reference_dir / "comparefiles"
            / f"corr.S40RTS..oupm.S40.{NAME}.filt.corr.dat"
        )
        _skip_if_missing(corr_path)

        fortran_depths, fortran_corr = _read_corr_total(corr_path)
        ref_sph = model_data.reference_coefficients

        for i in range(0, len(fortran_depths), 10):
            depth = fortran_depths[i]
            raw_filt = evaluate_at_depth(
                fortran_filt["coefficients"], DEGREE, depth
            )
            raw_ref = evaluate_at_depth(ref_sph, DEGREE, depth)
            _, total = correlation(raw_filt, raw_ref, DEGREE)

            assert total == pytest.approx(fortran_corr[i], abs=0.01), (
                f"Depth {depth} km: python={total:.6f}, "
                f"fortran={fortran_corr[i]:.6f}"
            )

    def test_per_degree_correlation_repar(self, fortran_repar, model_data, fortran_reference_dir):
        deg_path = (
            fortran_reference_dir / "comparefiles"
            / f"corr.S40RTS..inpm.S40.{NAME}.repar.corr.deg.dat"
        )
        _skip_if_missing(deg_path)

        fortran_depths, fortran_per_deg = _read_corr_per_degree(deg_path)
        ref_sph = model_data.reference_coefficients

        for depth_target in [500.0, 1000.0, 2000.0]:
            idx = np.argmin(np.abs(fortran_depths - depth_target))
            depth = fortran_depths[idx]

            raw_repar = evaluate_at_depth(
                fortran_repar["coefficients"], DEGREE, depth
            )
            raw_ref = evaluate_at_depth(ref_sph, DEGREE, depth)
            per_deg, _ = correlation(raw_repar, raw_ref, DEGREE)

            # Skip l=0 (often degenerate)
            np.testing.assert_allclose(
                per_deg[1: fortran_per_deg.shape[1]],
                fortran_per_deg[idx, 1:],
                atol=0.01,
                err_msg=f"Per-degree correlation mismatch at {depth} km",
            )


# ---------------------------------------------------------------------------
# Step 5: Point evaluation at original grid points
# ---------------------------------------------------------------------------


class TestPointEvaluation:
    """Evaluate .sph at original grid points and compare to Fortran outputfiles.

    Tests a representative subsample of points for speed.
    """

    @pytest.mark.parametrize("depth", [677, 1038, 1806])
    def test_repar_point_values(self, fortran_repar, fortran_reference_dir, depth):
        out_path = (
            fortran_reference_dir / "outputfiles"
            / f"inpm.S40.{NAME}.repar.{depth:04d}.dat"
        )
        _skip_if_missing(out_path)

        fortran_lat, fortran_lon, fortran_vals = _read_output_points(out_path)

        # Subsample every 5000th point for speed
        step = 5000
        lat = fortran_lat[::step]
        lon = fortran_lon[::step]
        expected = fortran_vals[::step]

        python_vals = evaluate_at_points(
            fortran_repar["coefficients"], DEGREE, float(depth), lat, lon
        )

        max_val = np.max(np.abs(expected))
        if max_val > 1e-10:
            max_diff = np.max(np.abs(python_vals - expected))
            rel_diff = max_diff / max_val
            assert rel_diff < 0.05, (
                f"Depth {depth} km: rel diff = {rel_diff:.4e} "
                f"(max_diff={max_diff:.4e}, max_val={max_val:.4e})"
            )

    @pytest.mark.parametrize("depth", [677, 1038])
    def test_filt_point_values(self, fortran_filt, fortran_reference_dir, depth):
        out_path = (
            fortran_reference_dir / "outputfiles"
            / f"oupm.S40.{NAME}.filt.{depth:04d}.dat"
        )
        _skip_if_missing(out_path)

        fortran_lat, fortran_lon, fortran_vals = _read_output_points(out_path)

        step = 5000
        lat = fortran_lat[::step]
        lon = fortran_lon[::step]
        expected = fortran_vals[::step]

        python_vals = evaluate_at_points(
            fortran_filt["coefficients"], DEGREE, float(depth), lat, lon
        )

        max_val = np.max(np.abs(expected))
        if max_val > 1e-10:
            max_diff = np.max(np.abs(python_vals - expected))
            rel_diff = max_diff / max_val
            assert rel_diff < 0.05, (
                f"Depth {depth} km: rel diff = {rel_diff:.4e}"
            )


# ---------------------------------------------------------------------------
# End-to-end pipeline
# ---------------------------------------------------------------------------


class TestEndToEnd:
    """Run the full pipeline from input layers through analysis and compare
    every output against Fortran reference data."""

    @pytest.fixture(scope="class")
    def result(self, fortran_reference_dir):
        example_dir = fortran_reference_dir / "examplemodel"
        _skip_if_missing(example_dir)
        from srts import tomographic_filter

        return tomographic_filter(
            str(example_dir), NAME, FIRST_LAYER, LAST_LAYER, DEGREE,
            run_analysis=True,
        )

    def test_repar_matches_fortran(self, result, fortran_repar):
        fort = fortran_repar["coefficients"]
        py = result["repar_coeffs"]
        max_val = np.max(np.abs(fort))
        max_diff = np.max(np.abs(py - fort))
        assert max_diff / max_val < 0.01, (
            f"End-to-end repar: relative diff = {max_diff / max_val:.4e}"
        )

    def test_filt_matches_fortran(self, result, fortran_filt):
        fort = fortran_filt["coefficients"]
        py = result["filt_coeffs"]
        max_val = np.max(np.abs(fort))
        max_diff = np.max(np.abs(py - fort))
        # Slightly looser: errors compound from reparameterization + filtering
        assert max_diff / max_val < 0.02, (
            f"End-to-end filt: relative diff = {max_diff / max_val:.4e}"
        )

    def test_analysis_power_repar(self, result, fortran_reference_dir):
        pwr_path = (
            fortran_reference_dir / "pwrfiles"
            / f"inpm.S40.{NAME}.repar.pwr.dat"
        )
        _skip_if_missing(pwr_path)

        _, fortran_power = _read_power_total(pwr_path)
        python_power = result["analysis"]["power_repar"] * 100

        np.testing.assert_allclose(
            python_power, fortran_power, rtol=0.05,
            err_msg="End-to-end reparameterized power mismatch",
        )

    def test_analysis_power_filt(self, result, fortran_reference_dir):
        pwr_path = (
            fortran_reference_dir / "pwrfiles"
            / f"oupm.S40.{NAME}.filt.pwr.dat"
        )
        _skip_if_missing(pwr_path)

        _, fortran_power = _read_power_total(pwr_path)
        python_power = result["analysis"]["power_filt"] * 100

        np.testing.assert_allclose(
            python_power, fortran_power, rtol=0.05,
            err_msg="End-to-end filtered power mismatch",
        )

    def test_analysis_correlation_repar(self, result, fortran_reference_dir):
        corr_path = (
            fortran_reference_dir / "comparefiles"
            / f"corr.S40RTS..inpm.S40.{NAME}.repar.corr.dat"
        )
        _skip_if_missing(corr_path)

        _, fortran_corr = _read_corr_total(corr_path)
        python_corr = result["analysis"]["corr_repar_ref"]

        np.testing.assert_allclose(
            python_corr, fortran_corr, atol=0.05,
            err_msg="End-to-end reparameterized correlation mismatch",
        )

    def test_analysis_correlation_filt(self, result, fortran_reference_dir):
        corr_path = (
            fortran_reference_dir / "comparefiles"
            / f"corr.S40RTS..oupm.S40.{NAME}.filt.corr.dat"
        )
        _skip_if_missing(corr_path)

        _, fortran_corr = _read_corr_total(corr_path)
        python_corr = result["analysis"]["corr_filt_ref"]

        np.testing.assert_allclose(
            python_corr, fortran_corr, atol=0.05,
            err_msg="End-to-end filtered correlation mismatch",
        )
