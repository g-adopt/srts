"""SH + depth splines → .sph model (replaces sphexp + sphadd).

Projects per-layer SH coefficients onto the 21-knot depth spline basis,
then sums all layers. All operations are vectorized.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg

from srts.expansion import expand_with_precomputed, precompute_expansion
from srts.io import read_layer_data
from srts.splines import SplineBasis, get_spline_basis


def _spline_projection_operator(spline_basis: SplineBasis) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precompute the spline ATA eigendecomposition (shared across all layers).

    Matches sphexp.f: evaluates 21 basis functions on a 2892-point uniform xd grid,
    builds (21x21) ATA matrix, eigendecomposes.

    Returns:
        V: active eigenvectors, shape (21, n_active)
        lam: active eigenvalues, shape (n_active,)
        z: spline basis on grid, shape (2892, 21)
    """
    maxp = 2891
    m = maxp // 2  # 1445
    mp = maxp + 1  # 2892

    xd_grid = np.arange(-m, m + 1, dtype=np.float64) / m  # 2891 points
    z = np.zeros((mp, 21))
    z[:2 * m + 1, :] = spline_basis.evaluate_all(xd_grid).T

    ata = z.T @ z
    eigenvalues, eigenvectors = scipy.linalg.eigh(ata)

    threshold = 1e-7 * eigenvalues[-1]
    active = eigenvalues > threshold
    return eigenvectors[:, active], eigenvalues[active], z


def project_layer_to_sph(
    raw_coeffs: np.ndarray,
    lmax: int,
    dep1: float,
    dep2: float,
    spline_V: np.ndarray,
    spline_lam: np.ndarray,
    spline_z: np.ndarray,
) -> np.ndarray:
    """Project one layer's .raw coefficients onto the depth spline basis.

    Args:
        raw_coeffs: Flat .raw coefficients, shape ((lmax+1)^2,).
        lmax: Maximum SH degree.
        dep1: Shallow depth boundary (km).
        dep2: Deep depth boundary (km).
        spline_V, spline_lam, spline_z: From _spline_projection_operator.

    Returns:
        shape (21, natd) — spline weights times raw coefficients per depth level.
    """
    mp = spline_z.shape[0]

    # Depth indicator vector
    id1 = round(2891 - dep2)
    id2 = round(2891 - dep1)
    d = np.zeros(mp)
    d[id1:id2 + 1] = 1.0

    # ATd and solve (damp=0 in sphexp.f)
    atd = spline_z.T @ d
    projections = spline_V.T @ atd
    x = spline_V @ (projections / spline_lam)  # (21,) spline weights

    # Outer product: each spline weight times the raw coefficients
    return x[:, np.newaxis] * raw_coeffs[np.newaxis, :]


def reparameterize(
    model_dir: str,
    name: str,
    first_layer: int,
    last_layer: int,
    degree: int,
    damp: float = 1.0,
) -> np.ndarray:
    """Full reparameterization: read layers, expand to SH, project to splines, sum.

    Replaces Steps 1+2 of dofilt_ES_new.

    The design matrix eigendecomposition and spline projection operator are
    computed once and reused across all layers.

    Returns:
        Coefficient array of shape (21, natd) — the reparameterized model.
    """
    data = read_layer_data(model_dir, name, first_layer, last_layer)
    depths = data["depths"]

    # Precompute SH expansion operator (shared grid across layers)
    precomp = precompute_expansion(data["lon"], data["lat"], degree, damp)

    # Precompute spline projection operator
    spline_basis = get_spline_basis()
    spline_V, spline_lam, spline_z = _spline_projection_operator(spline_basis)

    natd = (degree + 1) ** 2
    result = np.zeros((21, natd), dtype=np.float64)

    for iz_idx in range(data["nlayers"]):
        iz = first_layer + iz_idx
        dep1 = depths[iz - 1]
        dep2 = depths[iz]

        raw_coeffs = expand_with_precomputed(precomp, data["values"][iz_idx])
        sph_coeffs = project_layer_to_sph(
            raw_coeffs, degree, dep1, dep2, spline_V, spline_lam, spline_z
        )
        result += sph_coeffs

    return result
