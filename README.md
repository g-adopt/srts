# srts — Tomographic Filtering for Geodynamic Models

`srts` is a Python package for applying tomographic resolution filters to geodynamic models, following the methodology of Ritsema et al. (2007). It wraps the S12RTS, S20RTS, and S40RTS resolution operators, letting you reparameterize a geodynamic model onto the seismic tomography basis and apply the resolution matrix to produce a filtered model — i.e., what a seismologist would recover from your model given real data coverage.

## Installation

```bash
pip install srts
```

Dependencies are NumPy, SciPy, pyshtools, and h5py. The HDF5 model files (~4.8 GB for all three models) are downloaded automatically on first use from a public storage bucket, so you need internet access the first time you run anything that touches model data.

## Class-based API

For workflows where model data lives in memory — for example, output from a geodynamic simulation — three classes decompose the pipeline into independent steps. All public methods use pyshtools `cilm[2, lmax+1, lmax+1]` arrays as their interface, so output is directly compatible with pyshtools for analysis and visualization.

### Full example

```python
import numpy as np
from srts import S40RTS, SphericalHarmonicExpansion, DepthParameterization

# lon, lat: (npoints,) coordinate arrays in degrees
# values:   list of (npoints,) arrays, one per depth layer
# depth_boundaries: (nlayers+1,) depth boundaries in km

# Step 1 — expand each layer to spherical harmonics
expander = SphericalHarmonicExpansion(lon, lat, lmax=40)
layer_cilms = expander.expand_batch(np.stack(values))   # (nlayers, 2, 41, 41)

# Step 2 — project onto the 21-knot cubic spline depth basis
projector = DepthParameterization()
model = projector.reparameterize(list(layer_cilms), depth_boundaries)  # (21, 2, 41, 41)

# Step 3 — apply the S40RTS resolution matrix
s40 = S40RTS()
filtered = s40.filter(model)   # (ndp, 2, 41, 41)

# Evaluate at depths of interest
filtered_at_depths = DepthParameterization.evaluate_at_depths(filtered, np.arange(100, 2800, 50))
# → (ndepths, 2, 41, 41) — cilm arrays ready for pyshtools

# Visualize a slice with pyshtools
import pyshtools
coeffs = pyshtools.SHCoeffs.from_array(filtered_at_depths[10], normalization='ortho', csphase=-1)
coeffs.expand(grid='DH2').plot()
```

### SphericalHarmonicExpansion

Expands grid data into spherical harmonic coefficients using regularized least-squares. The normal equations are precomputed on construction and reused across layers, so expanding many layers on the same grid is efficient.

```python
from srts import SphericalHarmonicExpansion

expander = SphericalHarmonicExpansion(lon, lat, lmax=40)

# Single layer: (npoints,) → cilm (2, 41, 41)
cilm = expander.expand(values)

# All layers at once: (nlayers, npoints) → (nlayers, 2, 41, 41)
cilm_batch = expander.expand_batch(values_2d)
```

`lon` and `lat` are 1D arrays of coordinates in degrees. The setup cost scales with the number of unique latitudes in the grid. On a regular lon/lat mesh the operator is built from one matrix product per latitude band, which is far more efficient than the same number of irregularly scattered points.

### DepthParameterization

Projects per-layer SH coefficients onto the 21-knot cubic spline depth basis and evaluates spline-basis models at arbitrary depths.

```python
from srts import DepthParameterization

projector = DepthParameterization()

# Reparameterize: list of cilm arrays + depth boundaries in km → (21, 2, 41, 41)
model = projector.reparameterize(list(cilm_batch), depth_boundaries)

# Evaluate at a single depth
cilm_1000 = DepthParameterization.evaluate_at_depth(model, 1000.0)

# Evaluate at multiple depths: → (ndepths, 2, 41, 41)
cilm_stack = DepthParameterization.evaluate_at_depths(model, np.arange(100, 2800, 50))
```

`depth_boundaries` is a 1D array of shape `(nlayers+1,)` giving the top and bottom of each layer in km.

### TomographicFilter

Loads the eigenvectors and eigenvalues of the SxRTS inversion and applies the resolution matrix. Factory functions `S40RTS()`, `S20RTS()`, and `S12RTS()` provide construction with the default damping parameters from the original inversions.

```python
from srts import S40RTS, S20RTS, S12RTS

s40 = S40RTS()        # eps = 20e-4
s20 = S20RTS()        # eps = 35e-4
s12 = S12RTS()        # eps = 40e-4

# Custom damping
s40_soft = S40RTS(eps=0.005)

# Apply: (21, 2, 41, 41) → (ndp, 2, 41, 41)
filtered = s40.filter(model)

# Access the published reference model
reference = s40.reference_model   # (ndp, 2, 41, 41)
```

## File-based pipeline

For workflows that match the original Fortran `dofilt_ES_new` pipeline — where the model is stored as depth-slice ASCII files on disk — `tomographic_filter()` runs the entire pipeline in one call and optionally writes `.sph` output files in the same format as the original Fortran code.

```python
from srts import tomographic_filter

result = tomographic_filter(
    "geodyn/examplemodel",   # directory with layer .dat files and depth_layers.dat
    "examplefile.dvs",       # file prefix before .layer.NNN.dat
    first_layer=1,
    last_layer=64,
    degree=40,               # 12, 20, or 40
    output_dir=".",          # optional: writes .sph files to disk
    run_analysis=True,
)

repar    = result["repar_coeffs"]   # reparameterized model, shape (21, natd)
filtered = result["filt_coeffs"]    # filtered model, shape (ndp, natd)
analysis = result["analysis"]       # power spectra and correlations vs reference
```

When `output_dir` is set, two files are written in the Fortran-compatible `.sph` format:

- `inpm.S{degree}.{name}.repar.sph` — the reparameterized model
- `oupm.S{degree}.{name}.filt.sph` — the tomographically filtered model

These `.sph` files are readable by the S20RTS/S40RTS plotting tools available from Jeroen Ritsema's and Paula Koelemeijer's websites.

### Input file format

```
geodyn/
└── examplemodel/
    ├── depth_layers.dat
    ├── examplefile.dvs.layer.001.dat
    ├── examplefile.dvs.layer.002.dat
    └── ...
```

Each layer file contains three columns: longitude (−180 to 180), latitude (−90 to 90), and dVs in percent. The `depth_layers.dat` file lists depth boundaries in km, one per line, with one more entry than the number of layers. Depths should not exceed 2890 km.

### Analysis output

When `run_analysis=True`, the `"analysis"` key contains power spectra and correlations at 115 depths from 25 to 2875 km:

```python
a = result["analysis"]

a["depths"]           # (115,) depth array in km
a["power_repar"]      # (115,) RMS power of reparameterized model
a["power_filt"]       # (115,) RMS power of filtered model
a["power_ref"]        # (115,) RMS power of reference (SxRTS)
a["corr_repar_ref"]   # (115,) correlation with reference, reparameterized
a["corr_filt_ref"]    # (115,) correlation with reference, filtered
a["power_deg_repar"]  # (115, lmax+1) per-degree power, reparameterized
a["power_deg_filt"]   # (115, lmax+1) per-degree power, filtered
a["corr_deg_repar"]   # (115, lmax+1) per-degree correlation, reparameterized
a["corr_deg_filt"]    # (115, lmax+1) per-degree correlation, filtered
```

## Accuracy

The Python implementation is more accurate than the original Fortran pipeline. The Fortran code communicated between pipeline stages through ASCII files in `e12.4` format, giving roughly 4 significant digits per value — truncation that accumulates across 60+ depth layers. This package works in float64 throughout with no intermediate file I/O. When both implementations process the same `.sph` input, they agree to 6–7 significant digits on filtering, depth evaluation, power spectra, and cross-correlation.

## References

- Ritsema, J., van Heijst, H.J. & Woodhouse, J.H. (1999). Complex shear wave velocity structure imaged beneath Africa and Iceland. *Science*, 286, 1925–1928. (S20RTS)
- Ritsema, J., McNamara, A.K. & Bull, A.L. (2007). Tomographic filtering of geodynamic models: Implications for model interpretation and large-scale mantle structure. *Journal of Geophysical Research*, 112, B01303.
- Ritsema, J., Deuss, A., van Heijst, H.J. & Woodhouse, J.H. (2011). S40RTS: a degree-40 shear-velocity model for the mantle from new Rayleigh wave dispersion, teleseismic traveltime and normal-mode splitting function measurements. *Geophysical Journal International*, 184, 1223–1236. (S40RTS)
- Koelemeijer, P., Ritsema, J., Deuss, A. & van Heijst, H.J. (2016). SP12RTS: a degree-12 model of shear- and compressional-wave velocity for Earth's mantle. *Geophysical Journal International*, 204, 1024–1039. (S12RTS)
