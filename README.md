# srts — Tomographic Filtering for Geodynamic Models

`srts` is a Python package for applying tomographic resolution filters to geodynamic models, following the methodology of Ritsema et al. (2007). It wraps the S12RTS, S20RTS, and S40RTS resolution operators, letting you reparameterize a geodynamic model onto the seismic tomography basis and then apply the tomographic resolution matrix to produce a filtered model — i.e., what a seismologist would recover from your model given real data coverage.

The package replaces a legacy Fortran/C pipeline with pure Python (NumPy, SciPy, pyshtools, h5py), runs entirely in memory, and is more precise than the original because it avoids the intermediate ASCII file truncation that accumulated through the Fortran pipeline.

## Installation

```bash
pip install -e ".[test]"
```

Dependencies are NumPy, SciPy, pyshtools, and h5py. The HDF5 model files (~4.8 GB for all three models) are downloaded automatically on first use from a public storage bucket, so you need internet access the first time you run anything that touches model data.

## How it works

The filtering pipeline has three logical steps.

First, your model (expressed as depth slices on a lon/lat grid) is expanded into spherical harmonics up to the target degree using regularized least-squares. Second, the per-layer SH coefficients are projected onto the 21-knot cubic spline depth basis used by the SxRTS models — this is the reparameterization. Third, the SxRTS resolution matrix is applied to that reparameterized model, producing a filtered version that reflects the limited resolution of the actual tomographic inversion.

The package offers two interfaces for this: a file-based pipeline function that reads layer `.dat` files from disk, and a class-based API for workflows where the data already lives in memory.

## File-based pipeline

If your model is stored as depth-slice ASCII files (the format described below), `tomographic_filter()` runs the entire pipeline in one call:

```python
from srts import tomographic_filter

result = tomographic_filter(
    "geodyn/examplemodel",   # directory with layer .dat files and depth_layers.dat
    "examplefile.dvs",       # model name (file prefix before .layer.NNN.dat)
    first_layer=1,
    last_layer=64,
    degree=40,               # 12, 20, or 40
    output_dir=".",          # optional: writes .sph files to disk
    run_analysis=True,
)

repar   = result["repar_coeffs"]   # reparameterized model, shape (21, natd)
filtered = result["filt_coeffs"]   # filtered model, shape (ndp, natd)
analysis = result["analysis"]      # power spectra and correlations vs reference
```

When `output_dir` is set, the function writes two files in the Fortran-compatible `.sph` format:

- `inpm.S{degree}.{name}.repar.sph` — the reparameterized model
- `oupm.S{degree}.{name}.filt.sph` — the tomographically filtered model

These `.sph` files can be read by the S20RTS/S40RTS plotting tools available from Jeroen Ritsema's and Paula Koelemeijer's websites.

### Input file format

The layer `.dat` files live in a subdirectory of `geodyn/`:

```
geodyn/
└── examplemodel/
    ├── depth_layers.dat
    ├── examplefile.dvs.layer.001.dat
    ├── examplefile.dvs.layer.002.dat
    └── ...
```

Each layer file contains three columns: longitude (−180 to 180), latitude (−90 to 90), and dVs in percent. The `depth_layers.dat` file lists the depth boundaries in km, one per line, with one more entry than the number of layers — layer `N` represents the depth interval between line `N` and line `N+1`. Depths should not exceed 2890 km (the CMB).

An example dataset is included in `geodyn/examplemodel/`.

### Analysis output

When `run_analysis=True`, the `"analysis"` key in the result dict contains power spectra and correlations evaluated at 115 depths from 25 to 2875 km (25 km steps) for the reparameterized model, the filtered model, and the reference tomographic model:

```python
a = result["analysis"]

a["depths"]           # (115,) depth array in km
a["power_repar"]      # (115,) total RMS power of reparameterized model
a["power_filt"]       # (115,) total RMS power of filtered model
a["power_ref"]        # (115,) total RMS power of reference (SxRTS)
a["corr_repar_ref"]   # (115,) correlation with reference, reparameterized model
a["corr_filt_ref"]    # (115,) correlation with reference, filtered model
a["power_deg_repar"]  # (115, lmax+1) per-degree power, reparameterized
a["power_deg_filt"]   # (115, lmax+1) per-degree power, filtered
a["corr_deg_repar"]   # (115, lmax+1) per-degree correlation, reparameterized
a["corr_deg_filt"]    # (115, lmax+1) per-degree correlation, filtered
```

## Class-based API

For workflows where data lives in memory (e.g., output from a geodynamic simulation), three classes decompose the pipeline into independent steps. All public methods use pyshtools `cilm[2, lmax+1, lmax+1]` arrays as their interface, so the output is directly compatible with pyshtools for visualization.

### SphericalHarmonicExpansion

Expands grid data into spherical harmonic coefficients using regularized least-squares. The normal equations are precomputed on construction and reused across layers, so expanding many layers on the same grid is efficient.

```python
from srts import SphericalHarmonicExpansion

expander = SphericalHarmonicExpansion(lon, lat, lmax=40)

# Expand a single layer: (npoints,) → cilm (2, 41, 41)
cilm = expander.expand(values)

# Expand multiple layers at once: (nlayers, npoints) → (nlayers, 2, 41, 41)
cilm_batch = expander.expand_batch(values_2d)
```

`lon` and `lat` are 1D arrays of coordinates in degrees. For `lmax=40` on a typical geodynamic model grid, construction takes a few seconds; for `lmax=12` it is nearly instant.

### DepthParameterization

Projects per-layer SH coefficients onto the 21-knot cubic spline depth basis and evaluates spline-basis models at arbitrary depths.

```python
from srts import DepthParameterization

projector = DepthParameterization()

# Reparameterize: list of (2, lmax+1, lmax+1) cilm arrays + depth boundaries (km)
# Returns (21, 2, lmax+1, lmax+1)
model = projector.reparameterize(list(cilm_batch), depth_boundaries)

# Evaluate at a single depth (static method, no instance needed)
cilm_1000 = DepthParameterization.evaluate_at_depth(model, 1000.0)

# Evaluate at multiple depths at once: returns (ndepths, 2, lmax+1, lmax+1)
depths = np.arange(100, 2800, 50)
cilm_stack = DepthParameterization.evaluate_at_depths(model, depths)
```

`depth_boundaries` is a numpy array of shape `(nlayers+1,)` giving the top and bottom of each layer in km.

### TomographicFilter

Loads the eigenvectors and eigenvalues of the SxRTS inversion and applies the resolution matrix. Factory functions `S40RTS()`, `S20RTS()`, and `S12RTS()` provide convenient construction with the default damping parameters from the original inversions.

```python
from srts import S40RTS, S20RTS, S12RTS

s40 = S40RTS()        # eps = 20e-4  (default for S40RTS)
s20 = S20RTS()        # eps = 35e-4  (default for S20RTS)
s12 = S12RTS()        # eps = 40e-4  (default for S12RTS)

# Custom damping
s40_soft = S40RTS(eps=0.005)

# Apply the resolution matrix: (21, 2, 41, 41) → (ndp, 2, 41, 41)
filtered = s40.filter(model)

# Access the published reference model as cilm arrays: (ndp, 2, 41, 41)
reference = s40.reference_model
```

The `filter()` method takes a model in the 21-spline basis (the output of `DepthParameterization.reparameterize()`) and returns the filtered model at the depth levels covered by the resolution operator.

### Full example

```python
import numpy as np
from srts import S40RTS, SphericalHarmonicExpansion, DepthParameterization

# lon, lat: (npoints,) coordinate arrays in degrees
# values:   list of (npoints,) arrays, one per depth layer
# depth_boundaries: (nlayers+1,) depth boundaries in km

# Expand each layer to spherical harmonics
expander = SphericalHarmonicExpansion(lon, lat, lmax=40)
layer_cilms = expander.expand_batch(np.stack(values))  # (nlayers, 2, 41, 41)

# Project onto the 21-spline depth basis
projector = DepthParameterization()
model = projector.reparameterize(list(layer_cilms), depth_boundaries)

# Apply the S40RTS resolution matrix
s40 = S40RTS()
filtered = s40.filter(model)

# Evaluate at depths of interest
filtered_at_depths = DepthParameterization.evaluate_at_depths(filtered, np.arange(100, 2800, 50))
# → (ndepths, 2, 41, 41) — cilm arrays ready for pyshtools

# Visualize a slice with pyshtools
import pyshtools
coeffs = pyshtools.SHCoeffs.from_array(filtered_at_depths[10], normalization='ortho', csphase=-1)
coeffs.expand(grid='DH2').plot()
```

## Numerical precision

The Python implementation is strictly more precise than the original Fortran pipeline. The Fortran code communicates between pipeline stages through ASCII files in `e12.4` format, giving roughly 4 significant digits per value. Over 60+ depth layers, this truncation accumulates. The Python package works in float64 throughout with no intermediate file I/O, eliminating this source of error entirely.

When both implementations process the same `.sph` input (isolating precision effects from algorithm differences), they agree to 6–7 significant digits on filtering, depth evaluation, power spectra, and cross-correlation. The reparameterization step shows a ~0.7% global discrepancy that is entirely attributable to accumulated ASCII truncation in the Fortran path, not to any algorithmic difference.

## Testing

The test suite covers both unit tests for each module and Fortran validation tests that compare the Python output against reference files produced by the original `dofilt_ES_new` pipeline. To run the unit tests:

```bash
python -m pytest tests/ --ignore=tests/test_fortran_validation.py -v
```

See [TESTING.md](TESTING.md) for the full description of the test data, validation strategy, and CI configuration. See [DETAILS.md](DETAILS.md) for a thorough walkthrough of the class-based API.

## References

- Ritsema, J., van Heijst, H.J. & Woodhouse, J.H. (1999). Complex shear wave velocity structure imaged beneath Africa and Iceland. *Science*, 286, 1925–1928. (S20RTS)
- Ritsema, J., McNamara, A.K. & Bull, A.L. (2007). Tomographic filtering of geodynamic models: Implications for model interpretation and large-scale mantle structure. *Journal of Geophysical Research*, 112, B01303.
- Ritsema, J., Deuss, A., van Heijst, H.J. & Woodhouse, J.H. (2011). S40RTS: a degree-40 shear-velocity model for the mantle from new Rayleigh wave dispersion, teleseismic traveltime and normal-mode splitting function measurements. *Geophysical Journal International*, 184, 1223–1236. (S40RTS)
- Koelemeijer, P., Ritsema, J., Deuss, A. & van Heijst, H.J. (2016). SP12RTS: a degree-12 model of shear- and compressional-wave velocity for Earth's mantle. *Geophysical Journal International*, 204, 1024–1039. (S12RTS)
