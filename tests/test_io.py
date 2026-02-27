"""Tests for I/O module (io.py)."""

import pathlib
import tempfile

import numpy as np
import pytest

from srts.io import read_raw, write_raw, read_sph, write_sph, _parse_sph_header


class TestRaw:

    def test_write_read_round_trip(self, tmp_path):
        lmax = 4
        n = (lmax + 1) ** 2
        coeffs = np.random.default_rng(42).standard_normal(n)
        filepath = tmp_path / "test.raw"
        write_raw(filepath, lmax, coeffs)
        lmax_read, coeffs_read = read_raw(filepath)
        assert lmax_read == lmax
        np.testing.assert_allclose(coeffs_read, coeffs, atol=1e-7)

    def test_format_5_per_line(self, tmp_path):
        """Verify the Fortran format: 5 values per line, e16.8."""
        lmax = 2
        n = 9
        coeffs = np.arange(1, 10, dtype=np.float64) * 0.001
        filepath = tmp_path / "test.raw"
        write_raw(filepath, lmax, coeffs)
        with open(filepath) as f:
            header = f.readline()
            assert header.strip() == "2"
            line2 = f.readline()
            vals = line2.split()
            assert len(vals) == 5
            line3 = f.readline()
            vals = line3.split()
            assert len(vals) == 4


class TestSphHeader:

    def test_parse_standard_header(self):
        header = "          40 11111111111111111111111111111111111111111  24 000111111111111111111111"
        result = _parse_sph_header(header)
        assert result["lmax"] == 40
        assert result["ndep"] == 24
        assert result["ndmn"] == 4
        assert result["ndmx"] == 24
        assert result["ndp"] == 21

    def test_parse_header_partial_mask(self):
        header = "10 11111111111  24 000000000011111111111111"
        result = _parse_sph_header(header)
        assert result["lmax"] == 10
        assert result["ndmn"] == 11
        assert result["ndmx"] == 24


class TestSph:

    def test_write_read_round_trip(self, tmp_path):
        lmax = 4
        ndmn, ndmx = 4, 24
        ndp = ndmx - ndmn + 1
        natd = (lmax + 1) ** 2
        coeffs = np.random.default_rng(99).standard_normal((ndp, natd)) * 0.001
        filepath = tmp_path / "test.sph"
        write_sph(filepath, lmax, ndmn, ndmx, coeffs)
        result = read_sph(filepath)
        assert result["lmax"] == lmax
        assert result["ndmn"] == ndmn
        assert result["ndmx"] == ndmx
        assert result["ndp"] == ndp
        np.testing.assert_allclose(result["coefficients"], coeffs, atol=1e-3)

    def test_read_fortran_reference(self, fortran_reference_dir):
        """Read the Fortran reference .sph file and verify structure."""
        filepath = fortran_reference_dir / "inpm.S40.examplefile.dvs.repar.sph"
        if not filepath.exists():
            pytest.skip("Reference data not available")
        result = read_sph(filepath)
        assert result["lmax"] == 40
        assert result["ndmn"] == 4
        assert result["ndmx"] == 24
        assert result["ndp"] == 21
        assert result["coefficients"].shape == (21, 1681)

    def test_read_fortran_filt_reference(self, fortran_reference_dir):
        """Read the Fortran filtered reference .sph file."""
        filepath = fortran_reference_dir / "oupm.S40.examplefile.dvs.filt.sph"
        if not filepath.exists():
            pytest.skip("Reference data not available")
        result = read_sph(filepath)
        assert result["lmax"] == 40
        assert result["ndmn"] == 4
        assert result["ndmx"] == 24
        assert result["coefficients"].shape == (21, 1681)
