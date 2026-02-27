# Resolution Matrix - S[40, 20, 12]RTS

## Introduction
Scripts to reparametrise and tomographically filter a geodynamic model using the S20RTS, S40RTS or S12RTS resolution operators, see *Ritsema et al., 1999*, *Ritsema et al., 2011* and *Koelemeijer et al., 2016* for details regarding the models respectively, and *Ritsema et al., 2007* for details regarding the resolution operator.  These resolution filters allow you to filter your model
up to degree 12 (S12RTS), 20 (S20RTS) or degree 40 (S40RTS).

`dofilt_ES_new` is a bash scripts that filters only the S-wave structure. This script will give a reparametrised model file in the `.sph` format,
given by the file starting with `inpm`.  The parametrisation consists of **21 clamped cubic splines in depth and spherical harmonics up to the specified degree laterally**.
The script will also give a tomographically filtered model (again a".sph" file), starting with `oupm`, where the resolution operator has been applied.

These ".sph" files can be processed and plotted using the S20RTS/S40RTS plotting scripts available on the website of Jeroen Ritsema:
http://www.earth.lsa.umich.edu/~jritsema/S20RTS_plotting.tar.gz

Newer versions for both GMT4 and GMT5 are available on the website of Paula Koelemeijer:
https://www.earth.ox.ac.uk/~univ4152/files/SP12RTS_plotting.tar.gz
https://www.earth.ox.ac.uk/~univ4152/files/SP12RTS_plotting_GMT5.tar.gz

New in this script (compared to February 2016) is that the ".sph" files are automatically processed and expanded to give generate slices through the reparameterised ("inpm") and filtered ("oupm") models.  In addition, the correlation between these and the tomographic model
is calculated automatically, as well as their RMS power spectra.

## Installation

### Prerequisites

- **gfortran** (GNU Fortran compiler)
- **gcc** (GNU C compiler)
- **make** (build system)
- **Unix-like environment** (Linux, macOS, or Windows with WSL)

### Quick Start

```bash
# Clone or download the repository
# Navigate to the TOMOFILT directory
cd /path/to/tomofilt_new_ES

# Build everything (libraries and executables)
make all

# Set environment variable (add to your shell profile)
export TOMOFILT=/path/to/tomofilt_new_ES
```

### Detailed Installation Steps

#### 1. Set Environment Variable

Define the `TOMOFILT` environment variable pointing to the installation directory:

**For Bash** (add to `~/.bashrc` or `~/.bash_profile`):
```bash
export TOMOFILT=/path/to/tomofilt_new_ES
```

**For C-shell** (add to `~/.cshrc`):
```csh
setenv TOMOFILT /path/to/tomofilt_new_ES
```

> **Note**: Use `echo $SHELL` to determine which shell you're running.

#### 2. Build Libraries and Executables

The modernized build system allows you to build everything with a single command:

```bash
make all
```

This will:
- Build all four libraries: `libgeo.a`, `libHJVH.a`, `libS20.a`, `libSphFunc.a`
- Compile all 13 executables and place them in the `bin/` directory

**Alternative build options:**
```bash
# Build only libraries
make libraries

# Build only executables (requires libraries)
make executables

# Build individual libraries
make libgeo
make libHJVH
make libS20
make libSphFunc

# Check build status
make status

# Clean build artifacts
make clean
```

#### 3. Verify Installation

Check that all components were built successfully:

```bash
# Check build status
make status

# Verify executables are in bin/
ls bin/

# Test environment variable
echo $TOMOFILT
```

You should see 13 executables in the `bin/` directory and all libraries showing as "built" in the status report.

### Troubleshooting

If you encounter build issues:

```bash
# Check environment and dependencies
make doctor

# Clean and rebuild
make rebuild

# Get help
make help
```

## Input

The filtering script firstly reparametrises a given input model to the
same parametrisation as S20RTS/S40RTS/S12RTS.
The input model must therefore be evaluated in slices at specified depths
in the mantle, present in the directory "geodyn".

### Directory Structure

The input model should be organized in the `geodyn/` directory as follows:

```
geodyn/
└── [model_name]/
    ├── depth_layers.dat
    ├── [name].layer.001.dat
    ├── [name].layer.002.dat
    └── ...
```

**Example**: `geodyn/examplemodel/`

### File Format

Model slices should follow the naming convention:
```
[name].layer.[num].dat
```

Where:
- `[name]` = model run name (e.g., `examplefile.dvs`)
- `[num]` = slice number in 3-digit format (e.g., `001`, `010`, `100`)

**Example**: `geodyn/examplemodel/examplefile.dvs.layer.002.dat`

### Data Format

Each slice file should contain three columns:
```
longitude(-180,180)  latitude(-90,90)  dvs(%)
```

### Depth Configuration

The depth intervals must be specified in `geodyn/[model]/depth_layers.dat`:

- Slice `[num]` represents the interval between depth line `N` and `N+1`
- Example: slice `001` represents the model between depths on lines 1 and 2
- Typically evaluated at the midpoint between the two depths
- Number of depths should be 1 more than the number of model slices
- Maximum depth should not exceed 2890 km (CMB)

### Example Dataset

Sample model slices are provided in `geodyn/examplemodel/`:
- Files: `examplefile.dvs.layer.001.dat` to `examplefile.dvs.layer.126.dat`
- Corresponding depth file: `depth_layers.dat`

## Usage

### S-wave Structure Filtering

Filter S-wave velocity structure using S12RTS (degree 12), S20RTS (degree 20), or S40RTS (degree 40) resolution operators.

#### Command Syntax

```bash
./dofilt_ES_new [model] [name] [firstlay] [lastlay] [degree]
```

#### Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `model` | Subdirectory name in `geodyn/` | `examplemodel` |
| `name` | File prefix (before `.layer.num.dat`) | `examplefile.dvs` |
| `firstlay` | First layer below crust | `1` |
| `lastlay` | Last layer above CMB | `64` |
| `degree` | Resolution degree (12, 20, or 40) | `40` |

#### Example

```bash
./dofilt_ES_new examplemodel examplefile.dvs 1 64 40
```

#### Output Files

The script generates two `.sph` files:

1. **Reparameterized model**: `inpm.S[degree].[name].repar.sph`
2. **Filtered model**: `oupm.S[degree].[name].filt.sph`

**Example output**:
- `inpm.S40.examplefile.dvs.repar.sph`
- `oupm.S40.examplefile.dvs.filt.sph`

### Optional Outputs

Optional model outputs are automatically calculated when the parameters `do_expand` and `do_compare` in the `dofilt_ES_new` script are set to `1`.

#### Model Expansion (`do_expand = 1`)

When enabled, the reparameterized and filtered model files are evaluated at the given geographic input points for each depth in `depth_layers.dat`.

> **Note**: Output values don't exactly correspond to input values, as they represent depth intervals rather than specific depths.

**Output location**: `geodyn/outputfiles/`

**File naming**:
- `inpm.S[degree].[name].repar.[depth].dat`
- `oupm.S[degree].[name].filt.[depth].dat`

**Example**:
- `inpm.S40.examplefile.dvs.repar.2840.dat`
- `oupm.S40.examplefile.dvs.filt.2840.dat`

#### Model Comparison (`do_compare = 1`)

When enabled, performs comparisons between the reparameterized/filtered geodynamic files and the tomographic model (S12RTS/S20RTS/S40RTS).

##### Spherical Harmonic Expansions

Calculated for every 25 km depth and stored in `geodyn/rawfiles/`:

- `S[degree]RTS.[depth].raw`
- `inpm.S[degree].[name].repar.[depth].raw`
- `oupm.S[degree].[name].filt.[depth].raw`

##### Power Spectra

Generated from spherical harmonic files and stored in `geodyn/pwrfiles/`:

| File Type | Description | Example |
|-----------|-------------|---------|
| `*.pwr.dat` | Overall RMS power | `S40RTS.pwr.dat` |
| `*.pwr.deg.dat` | Power per degree | `S40RTS.pwr.deg.dat` |

**Complete file list**:
- `S[degree]RTS.pwr.dat` / `S[degree]RTS.pwr.deg.dat`
- `inpm.S[degree].[name].repar.pwr.dat` / `inpm.S[degree].[name].repar.pwr.deg.dat`
- `oupm.S[degree].[name].filt.pwr.dat` / `oupm.S[degree].[name].filt.pwr.deg.dat`

##### Correlation Analysis

Correlation calculated between tomographic model and reparameterized/filtered `.sph` files.

**Output location**: `geodyn/comparefiles/`

| File Type | Description | Example |
|-----------|-------------|---------|
| `*.corr.dat` | Overall correlation | `corr.S40RTS..inpm.S40.examplefile.dvs.repar.corr.dat` |
| `*.corr.deg.dat` | Correlation per degree | `corr.S40RTS..inpm.S40.examplefile.dvs.repar.corr.deg.dat` |

##### Lateral Expansion

Files are expanded laterally at all depths using equidistant sampling (1°, 2°, or 5° intervals).

**Output location**: `geodyn/slices/`

**File naming**:
- `S[degree]RTS.[depth].raw.dat`
- `inpm.S[degree].[name].repar.[depth].raw.dat`
- `oupm.S[degree].[name].filt.[depth].raw.dat`

**Example**:
- `S40RTS.0300.raw.dat`
- `inpm.S40.examplefile.repar.0300.raw.dat`
- `oupm.S40.examplefile.filt.0300.raw.dat`

## Python Package (srts)

A pure Python reimplementation of the entire Fortran pipeline is available as the `srts` package. It provides the same tomographic filtering functionality with no compiled dependencies beyond NumPy, SciPy, pyshtools, and h5py.

### Installation

```bash
pip install -e ".[test]"
```

### Usage

```python
from srts import tomographic_filter

result = tomographic_filter(
    "geodyn/examplemodel",    # model directory
    "examplefile.dvs",        # model name
    first_layer=1,
    last_layer=64,
    degree=40,                # 12, 20, or 40
    run_analysis=True,
)

repar = result["repar_coeffs"]   # reparameterized model (21 x natd)
filt = result["filt_coeffs"]     # filtered model (ndp x natd)
analysis = result["analysis"]    # power spectra, correlations vs reference
```

Individual pipeline steps are also available:

```python
from srts.parameterization import reparameterize
from srts.model_data import load_model_data
from srts.filtering import apply_resolution_matrix, extract_sp_from_spt
from srts.analysis import evaluate_at_depth, power_spectrum, correlation
```

### Numerical Precision

The Python implementation is more precise than the Fortran pipeline. The original code pipes intermediate results through ASCII files in Fortran `e12.4` format, which provides roughly 4 significant digits per value. Over the course of 60+ layer iterations, these rounding errors accumulate. The Python package uses float64 throughout with no intermediate truncation.

When both implementations start from the same `.sph` file (isolating the effect of intermediate precision), agreement is 6 to 7 significant digits for filtering, depth evaluation, power spectra, and cross-correlation. The reparameterization step, which accumulates all 60+ layers of ASCII truncation in the Fortran path, shows a roughly 0.7% global discrepancy that is entirely attributable to the Fortran's precision loss.

### Validation

The test suite uses a two-tier strategy. 76 unit tests verify each module in isolation (splines, spherical harmonic expansion, I/O round-trips, etc.). 36+ Fortran validation tests compare every pipeline stage against reference outputs produced by the original Fortran `dofilt_ES_new` on an example geodynamic model with S40RTS.

The validation tests are organized by pipeline stage. For isolated steps where both codes read the same `.sph` input file, tolerances are tight (relative differences below 1e-5). For the full end-to-end pipeline, tolerances are slightly looser to accommodate the genuine precision gap from accumulated ASCII truncation in the Fortran path. Additionally, a multi-degree test class validates the filtering pipeline for all three resolutions (S12RTS, S20RTS, S40RTS) using self-consistency checks on the published reference models.

To run the full test suite:

```bash
python -m pytest tests/ -v
```

## Support

If you encounter any problems or have questions, please contact:

**Paula Koelemeijer**
pjkoelemeijer@cantab.net

**Jeroen Ritsema**
jritsema@umich.edu

## References

- *Ritsema et al., 1999* - S20RTS model details
- *Ritsema et al., 2007* - Resolution operator details
- *Ritsema et al., 2011* - S40RTS model details
- *Koelemeijer et al., 2016* - S12RTS model details
