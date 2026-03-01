# Testing the srts package

The srts test suite has two layers: unit tests that verify each module in isolation, and Fortran validation tests that compare every pipeline stage against the reference Fortran outputs produced by `dofilt_ES_new` on the example model with S40RTS (degree 40).

## Running the tests

```bash
# All unit tests (~14 seconds)
python -m pytest tests/ --ignore=tests/test_fortran_validation.py -v

# Fortran validation, fast subset (~5 seconds)
python -m pytest tests/test_fortran_validation.py -v -k "not Reparameterization and not EndToEnd"

# Full validation including reparameterization and end-to-end (~7 minutes)
python -m pytest tests/test_fortran_validation.py -v

# Everything
python -m pytest tests/ -v
```

## Test data and CI

The test suite depends on two categories of external data that are too large for git:

**HDF5 model data** (S12RTS.h5, S20RTS.h5, S40RTS.h5 — totaling ~4.8 GB) lives in `src/srts/data/` and is required by the unit tests that exercise model loading, filtering, and the pipeline. These files are downloaded automatically from S3 on first use via lazy download in `model_data.py`, so no manual setup is needed for local development.

**Fortran reference data** is organized into two tiers on S3 at `s3://gadopt/srts/`:

- `tier_1.tar.gz` (~29 MB): Contains the reparameterized and filtered .sph files, .raw depth slices at 7 depths, power spectrum files, correlation files, and point evaluation outputs. This covers validation steps 3-6 and the pipeline integration tests.

- `tier_2.tar.gz` (~317 MB): Contains the 64 input layer .dat files, depth_layers.dat, and points.dat for the example model. This is only needed for the reparameterization tests (TestReparameterization) and the end-to-end pipeline test (TestEndToEnd).

The test fixtures resolve the Fortran reference data location in this order:
1. `<repo_root>/test-data/fortran-reference/` — used by CI (downloaded and extracted there)
2. `SRTS_FORTRAN_GEODYN` environment variable — for custom locations
3. Hardcoded fallback to the local Fortran tree — for development on the original machine

Tests that cannot find reference data are automatically skipped rather than failing.

## CI workflow

GitHub Actions runs on every push and nightly:

- **Push**: Downloads tier_1 + all HDF5 files, runs all unit tests and the fast subset of Fortran validation (everything except TestReparameterization and TestEndToEnd).
- **Nightly** (cron): Downloads tier_1 + tier_2 + all HDF5 files, runs the full test suite including reparameterization and end-to-end tests.

Both configurations test against Python 3.12 and 3.13.

## Running tests without S3 data

If you have the Fortran reference data at a custom location, point the `SRTS_FORTRAN_GEODYN` environment variable at it:

```bash
export SRTS_FORTRAN_GEODYN=/path/to/geodyn
python -m pytest tests/ -v
```

Without any reference data, unit tests that need HDF5 files will trigger lazy download from S3, and Fortran validation tests will be skipped.

## Uploading test data to S3

The `tools/upload_test_data.sh` script packages the Fortran reference data into the two tarballs and uploads everything (tarballs + HDF5 files) to the DigitalOcean Spaces bucket. It requires `s3cmd` configured with `~/.s3cfg-gadopt`.


## Unit tests (76 tests across 9 files)

These test individual modules without requiring external Fortran reference data.

**test_splines.py** (11 tests) validates the 21-knot cardinal cubic spline depth basis. Checks depth-to-xd conversion (CMB maps to -1.0, Moho to +1.0), round-trip preservation, monotonicity, the anti-diagonal identity at knots, partition of unity, and smoothness bounds.

**test_coeffs.py** (10 tests) covers coefficient format conversions between the Fortran flat layout and pyshtools cilm arrays. Verifies index table dimensions, (l,m) value extraction, flat index layout for degree 2, normylm values (1.0 for m=0, sqrt(2) for m>0), round-trip conversions for both standard and raw formats at lmax in {2, 4, 10, 20}, and the 0.01 scaling plus sqrt(2) normylm factors in raw format.

**test_analysis.py** (11 tests) exercises the spectral analysis functions. Power spectrum tests check zero input, single-degree isolation, non-negativity, and batch-vs-single consistency. Correlation tests verify self-correlation (1.0), anticorrelation (-1.0), low correlation for random vectors, and batch matching. Depth evaluation tests confirm correct output shapes and that batch evaluation matches repeated single-depth calls.

**test_filtering.py** (8 tests) tests the resolution matrix and parameter extraction. Validates that default damping factors exist for degrees 12, 20, and 40 (S40 uses 20e-4), that output shapes are correct, that zero input produces zero output, that the filtering is linear, and that the eigenvalue cutoff is sensible. Also checks that extracting the S-wave parameter from the multi-parameter output gives the right shape and depth range.

**test_io.py** (8 tests) verifies I/O for the .raw and .sph Fortran formats. Tests write/read round-trips for both formats, the Fortran 5-values-per-line layout, .sph header parsing (standard and partial depth masks), and reading actual Fortran-generated reference files when available.

**test_model_data.py** (7 tests) covers loading of tomographic model metadata from the HDF5 data files. Parametrized across degrees {12, 20, 40}, checks that loading succeeds with correct lmax and natd, that degree 30 raises an error, and that S40RTS has the expected metadata (ndep=21, natd=1681). Also validates that eigenvalues are positive and decreasing, eigenvector shapes are correct, weight matrix dimensions match, and the reference S40RTS model loads properly.

**test_parameterization.py** (7 tests) tests the depth spline projection operator. Verifies operator shape (V is 21 x n_basis), positive eigenvalues, correct output shape from projecting a layer, that deep layers (2800-2891 km) activate deep spline basis functions (indices 16-20), that shallow layers (25-100 km) activate shallow ones (indices 0-5), and that the projection scales linearly with input.

**test_expansion.py** (4 tests) validates spherical harmonic expansion from grid data to coefficients. A constant field should produce only c_00 (all other terms below 0.5% of c_00), and Y_20 should produce a dominant c_20 term. The precomputed expansion operator is verified to match direct expansion exactly and to work for different data vectors on the same grid.

**test_pipeline.py** (4 tests) provides integration-level checks against Fortran reference .sph files. Applies the resolution matrix to the Fortran reparameterized model and verifies the filtered output matches within 1% of the maximum value (accounting for e12.4 format precision). Checks that the power spectrum at 1000 km matches within 5%, that the filtered model shows positive correlation with S40RTS, and that a .sph write/read round-trip preserves data within 5e-4.


## Fortran validation tests (36 tests in test_fortran_validation.py)

These compare every stage of the Python pipeline against Fortran reference outputs
from running `dofilt_ES_new` on the example model. The reference data includes
.sph files, .raw depth slices, power spectra, correlation files, and point
evaluations at the original grid coordinates.

**TestReparameterization** (3 tests) runs the full Python reparameterization (reading 64 layer .dat files, expanding each to spherical harmonics, projecting onto the spline basis, and summing) and compares the resulting .sph coefficients against the Fortran output. Checks that shapes match, overall coefficients agree within 1% relative error, and each of the 21 depth levels correlates above 0.99 with the Fortran reference. This is the slowest test group as it processes 64 layers at lmax=40.

**TestFiltering** (2 tests) applies the resolution matrix to the Fortran-generated reparameterized .sph and compares the filtered output against the Fortran filtered .sph. Validates that ndmn, ndmx, and overall coefficient differences stay within 1% of the maximum value, and that each depth level correlates above 0.999 with the Fortran reference.

**TestDepthEvaluation** (14 tests) evaluates both the reparameterized and filtered .sph files at 7 depths (100, 500, 675, 1000, 1500, 2000, 2500 km) and compares the resulting flat coefficient arrays against Fortran .raw files. Since both Python and Fortran read the same .sph file (e12.4 precision), differences are purely numerical and must stay below 0.01% relative error.

**TestPowerSpectrum** (3 tests) computes total and per-degree power from the .sph files at multiple depths and compares against the Fortran power files. Total power is checked at every 10th depth with 1% relative tolerance. Per-degree spectra are compared at 500, 1000, and 2000 km depth.

**TestCorrelation** (3 tests) cross-correlates the reparameterized and filtered models against the real S40RTS tomography and compares both total and per-degree correlation against Fortran reference files. Total correlation is checked at every 10th depth with absolute tolerance of 0.01.

**TestPointEvaluation** (5 tests) evaluates the .sph files at specific (lat, lon, depth) points drawn from the original geodynamic model grid and compares against Fortran sph2v_input output. Tests use every 5000th point for speed, covering 3 depths for the reparameterized model (677, 1038, 1806 km) and 2 for the filtered model (677, 1038 km). Relative differences must stay below 5%.

**TestEndToEnd** (6 tests) runs the complete Python pipeline from raw input layers through filtering and analysis, then compares every output against Fortran. Checks reparameterized coefficients (1% tolerance), filtered coefficients (2% tolerance to account for compounding errors), total power at all 115 depths (5% tolerance), and total correlation at all 115 depths (0.05 absolute tolerance).


## Key normalization issue and fix

The Fortran LEGNDR subroutine uses fully orthonormalized associated Legendre functions with the Condon-Shortley phase (-1)^m. The pyshtools PlmBar function returns geodesy-normalized functions without that phase and with an extra sqrt(2) factor for m>0. The relationship is:

    LEGNDR(l, m=0) = PlmBar(l, 0) / sqrt(4*pi)
    LEGNDR(l, m>0) = (-1)^m * PlmBar(l, m) / sqrt(2) / sqrt(4*pi)

This required using `csphase=-1` in all PlmBar calls and dividing m>0 values by sqrt(2) in both the expansion code (expansion.py) and the point evaluation code (analysis.py). Without this fix, point evaluations produced correct values at the poles (where only m=0 terms contribute) but wildly wrong values at all other latitudes.
