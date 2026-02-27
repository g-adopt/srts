"""End-to-end tomographic filtering pipeline (replaces dofilt_ES_new)."""

from __future__ import annotations

import numpy as np

from srts.analysis import depth_profile_analysis
from srts.filtering import DEFAULT_EPS, apply_resolution_matrix, extract_sp_from_spt
from srts.io import read_sph, write_sph
from srts.model_data import ModelData, load_model_data
from srts.parameterization import reparameterize


def tomographic_filter(
    model_dir: str,
    name: str,
    first_layer: int,
    last_layer: int,
    degree: int,
    eps: float | None = None,
    output_dir: str | None = None,
    run_analysis: bool = True,
) -> dict:
    """Run the complete tomographic filtering pipeline.

    Equivalent to dofilt_ES_new.

    Args:
        model_dir: Path to the geodyn model directory containing layer .dat files
            and depth_layers.dat.
        name: Model name (e.g. "examplefile.dvs").
        first_layer: First layer index (1-based).
        last_layer: Last layer index (exclusive).
        degree: Maximum SH degree (12, 20, or 40).
        eps: Damping parameter. If None, uses the default for the degree.
        output_dir: If set, writes .sph files in Fortran format.
        run_analysis: If True, runs the full depth profile analysis.

    Returns:
        dict with keys:
            'repar_coeffs': reparameterized model, shape (21, natd)
            'filt_coeffs':  filtered model, shape (ndp, natd)
            'filt_ndmn':    first depth index of filtered model
            'filt_ndmx':    last depth index of filtered model
            'model_data':   ModelData used for filtering
            'analysis':     depth profile analysis dict (if run_analysis=True)
    """
    if eps is None:
        eps = DEFAULT_EPS.get(degree)
        if eps is None:
            raise ValueError(f"No default eps for degree={degree}; provide explicitly")

    # Step 1+2: Reparameterize
    repar = reparameterize(model_dir, name, first_layer, last_layer, degree)

    # Load resolution operator data
    model_data = load_model_data(degree)

    # Step 3: Apply resolution matrix
    filtered = apply_resolution_matrix(repar, model_data, eps)
    filt_coeffs, ndmn, ndmx = extract_sp_from_spt(filtered, model_data)

    # Write output files if requested
    if output_dir is not None:
        import pathlib
        out = pathlib.Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)

        write_sph(
            out / f"inpm.S{degree}.{name}.repar.sph",
            degree, 4, 24, repar,
        )
        write_sph(
            out / f"oupm.S{degree}.{name}.filt.sph",
            degree, ndmn, ndmx, filt_coeffs,
        )

    result = {
        "repar_coeffs": repar,
        "filt_coeffs": filt_coeffs,
        "filt_ndmn": ndmn,
        "filt_ndmx": ndmx,
        "model_data": model_data,
    }

    # Step 6: Analysis
    if run_analysis:
        ref_sph_data = read_sph_reference(model_data)
        result["analysis"] = depth_profile_analysis(
            repar, filt_coeffs, ref_sph_data, degree,
        )

    return result


def read_sph_reference(model_data: ModelData) -> np.ndarray:
    """Extract the reference model coefficients in the shape expected by analysis.

    The reference model from HDF5 is stored as (ndp, natd). We need it padded
    to match the .sph 21-depth convention (indices 4..24) if necessary.
    """
    return model_data.reference_coefficients
