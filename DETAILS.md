# Class-Based API — Detailed Guide

The `srts` package exposes three classes that decompose the tomographic filtering pipeline into independent, composable steps. Each class wraps the same validated numerical routines used by the file-based `tomographic_filter()` function, but operates entirely in memory using pyshtools `cilm[2, lmax+1, lmax+1]` arrays as the public interface. Nothing touches disk unless you explicitly write files yourself.

This guide walks through each class, explains the data flow, and shows how to wire them together for common use cases.

## Coefficient format: cilm arrays

Every public method accepts and returns pyshtools-style `cilm` arrays of shape `(2, lmax+1, lmax+1)`. The first axis separates cosine (`cilm[0, l, m]`) and sine (`cilm[1, l, m]`) coefficients. Only entries with `m <= l` carry meaningful values. This is the standard pyshtools convention and can be passed directly to `pyshtools.SHCoeffs.from_array()` for visualization, grid expansion, etc.

Internally, the package uses a Fortran-derived flat format for the actual numerical work. The `cilm_stack_to_internal()` and `internal_to_cilm_stack()` functions in `srts.coeffs` handle the conversion, but you should never need to call them directly — the classes handle it.

## SphericalHarmonicExpansion

This class expands irregularly-sampled grid data (longitude, latitude, values) into spherical harmonic coefficients using a regularized least-squares approach. It precomputes the normal equations and eigendecomposition on construction, so expanding multiple data vectors on the same grid is fast.

### Construction

```python
from srts import SphericalHarmonicExpansion

expander = SphericalHarmonicExpansion(lon, lat, lmax=40, damp=1.0)
```

`lon` and `lat` are 1D numpy arrays of coordinates in degrees (longitude in [0, 360] or [-180, 180], latitude in [-90, 90]). `lmax` is the maximum spherical harmonic degree. `damp` controls the Tikhonov regularization (the default of 1.0 matches the Fortran pipeline).

The constructor calls `precompute_expansion()` internally, which involves an eigendecomposition of the `(lmax+1)^2 × (lmax+1)^2` normal matrix. For `lmax=40` this takes a few seconds; for `lmax=12` it's nearly instant. Once built, the object is reusable for any number of data vectors on that grid.

### Expanding a single layer

```python
cilm = expander.expand(values)  # values: (npoints,) → cilm: (2, 41, 41)
```

### Expanding multiple layers at once

```python
cilm_batch = expander.expand_batch(values_2d)  # (nlayers, npoints) → (nlayers, 2, 41, 41)
```

This simply loops `expand()` internally, but it keeps the code concise.

### Typical context

You'd use this when your data comes from a simulation or some other source as gridded values at discrete depth layers. The output cilm arrays are what `DepthParameterization` expects as input.

## DepthParameterization

This class handles the projection of per-layer spherical harmonic coefficients onto the 21-knot cubic spline depth basis used by SxRTS models. It also provides static methods to evaluate a spline-basis model at arbitrary depths.

### Construction

```python
from srts import DepthParameterization

projector = DepthParameterization()
```

The constructor precomputes the spline basis functions and their projection operator (a small 21×21 eigendecomposition). This is cheap and only done once.

### Projecting a single layer

```python
sph = projector.project_layer(cilm, depth_top=500.0, depth_bottom=700.0)
# cilm: (2, lmax+1, lmax+1) → sph: (21, 2, lmax+1, lmax+1)
```

This takes one layer's SH coefficients and returns its contribution to the 21-spline basis, weighted by how much of the depth interval falls within each spline's support. The output has 21 depth levels, each carrying the spline-weighted coefficients.

### Reparameterizing a full model

```python
model = projector.reparameterize(layer_cilms, depth_boundaries)
```

`layer_cilms` is a list of `nlayers` cilm arrays. `depth_boundaries` is a numpy array of shape `(nlayers + 1,)` giving the depth boundaries in km — `depth_boundaries[i]` and `depth_boundaries[i+1]` define the top and bottom of layer `i`. The method sums all per-layer projections into a single `(21, 2, lmax+1, lmax+1)` array, which is the reparameterized model ready for filtering.

### Evaluating at depth

These are static methods, so you don't need an instance:

```python
# Single depth
cilm_1000 = DepthParameterization.evaluate_at_depth(model, 1000.0)
# → (2, lmax+1, lmax+1)

# Multiple depths at once
depths = np.array([500.0, 1000.0, 1500.0, 2000.0])
cilm_stack = DepthParameterization.evaluate_at_depths(model, depths)
# → (4, 2, lmax+1, lmax+1)
```

These work on any spline-basis model — the reparameterized model before filtering, the filtered model after filtering, or even the published SxRTS reference. The input just needs to have shape `(ndp, 2, lmax+1, lmax+1)` where `ndp` is the number of depth spline levels (21 for a full model, fewer for a filtered output).

## TomographicFilter

This is the resolution matrix operator. It loads the eigenvectors, eigenvalues, and data weights for a given SxRTS model and applies the resolution matrix `R` to produce a filtered model.

### Construction

```python
from srts import TomographicFilter, S40RTS, S20RTS, S12RTS

# Using the class directly
filt = TomographicFilter(degree=40, eps=20e-4)

# Or using factory functions (which use default eps values)
s40 = S40RTS()        # eps = 20e-4
s20 = S20RTS()        # eps = 35e-4
s12 = S12RTS()        # eps = 40e-4

# Custom damping
s40_custom = S40RTS(eps=0.005)
```

Construction loads the HDF5 model data file (downloaded from S3 on first use, cached thereafter). This includes the eigenvectors of the tomographic inversion and the smoothness weights.

### Properties

```python
s40.degree           # 40
s40.lmax             # 40
s40.eps              # 0.002
s40.reference_model  # (ndp, 2, 41, 41) — the published S40RTS model as cilm arrays
```

The `reference_model` property gives you the actual tomographic model in cilm format, which is useful for computing correlations or power spectra against your filtered output.

### Filtering

```python
filtered = s40.filter(model)
# model: (21, 2, 41, 41) → filtered: (ndp, 2, 41, 41)
```

The input must be in the 21-spline basis (i.e., the output of `DepthParameterization.reparameterize()`). The output has `ndp` depth levels, which is determined by the resolution operator's active depth range (typically 21 for S40RTS).

## Putting it all together

### Full pipeline from grid data

This is the most common use case: you have a geodynamic model evaluated on a lon/lat grid at multiple depth layers, and you want to see what a seismic tomographer would recover.

```python
import numpy as np
from srts import S40RTS, SphericalHarmonicExpansion, DepthParameterization

# Your data: lon/lat coordinates, values per layer, depth boundaries
# lon, lat: (npoints,) in degrees
# values: list of (npoints,) arrays, one per layer
# depth_boundaries: (nlayers+1,) in km

# Step 1: Expand each layer to spherical harmonics
expander = SphericalHarmonicExpansion(lon, lat, lmax=40)
layer_cilms = [expander.expand(v) for v in values]

# Step 2: Project onto the 21-spline depth basis
projector = DepthParameterization()
model = projector.reparameterize(layer_cilms, depth_boundaries)

# Step 3: Apply the resolution matrix
s40 = S40RTS()
filtered = s40.filter(model)

# Step 4: Evaluate at depths of interest
depths = np.arange(100, 2800, 50)
filtered_at_depths = DepthParameterization.evaluate_at_depths(filtered, depths)
# → (ndepths, 2, 41, 41) — ready for pyshtools visualization
```

### Comparing against the published tomographic model

```python
import pyshtools

s40 = S40RTS()
ref = s40.reference_model  # published S40RTS as cilm

# Evaluate the reference at 1000 km
ref_1000 = DepthParameterization.evaluate_at_depth(ref, 1000.0)

# Expand to a grid using pyshtools
coeffs = pyshtools.SHCoeffs.from_array(ref_1000, normalization='ortho', csphase=-1)
grid = coeffs.expand(grid='DH2')
grid.plot()
```

### Filtering an already-parameterized model

If you already have a model in the 21-spline SH basis (perhaps from a previous run or from manipulating coefficients directly), you can skip straight to filtering:

```python
s20 = S20RTS()
filtered = s20.filter(my_sph_model)  # my_sph_model: (21, 2, 21, 21)
```

### Using multiple resolution operators

The three filters are independent objects, so you can apply all of them to the same input model to see how data coverage affects what's recovered:

```python
from srts import S12RTS, S20RTS, S40RTS, DepthParameterization

# Assuming `model_12`, `model_20`, `model_40` are reparameterized
# at the appropriate lmax (12, 20, 40 respectively)

filtered_12 = S12RTS().filter(model_12)
filtered_20 = S20RTS().filter(model_20)
filtered_40 = S40RTS().filter(model_40)

# Compare at 1000 km
for label, filt in [("S12", filtered_12), ("S20", filtered_20), ("S40", filtered_40)]:
    c = DepthParameterization.evaluate_at_depth(filt, 1000.0)
    print(f"{label}: max |coeff| = {np.max(np.abs(c)):.4f}")
```

### Using expand_batch for efficiency

When you have many layers, `expand_batch` is cleaner than a list comprehension:

```python
values_2d = np.stack(values)  # (nlayers, npoints)
layer_cilms = expander.expand_batch(values_2d)  # (nlayers, 2, lmax+1, lmax+1)
model = projector.reparameterize(list(layer_cilms), depth_boundaries)
```

## Relationship to the functional API

The class-based API is a thin wrapper around the same functions used by `tomographic_filter()`. No algorithms were changed, no numerical code was duplicated. Specifically:

- `SphericalHarmonicExpansion.expand()` calls `precompute_expansion()` and `expand_with_precomputed()` from `srts.expansion`
- `DepthParameterization.project_layer()` calls `project_layer_to_sph()` from `srts.parameterization`
- `TomographicFilter.filter()` calls `apply_resolution_matrix()` and `extract_sp_from_spt()` from `srts.filtering`

The only added layer is the cilm ↔ internal format conversion at the boundary of each method. This conversion is exact (no floating point loss) and verified by round-trip tests.

The file-based `tomographic_filter()` function remains available and unchanged. Both APIs coexist in the same package and can be mixed freely.
