"""Tomographic resolution matrix application (replaces mk3d_res.f + spt2sph.f).

Implements m' = R m where R is the tomographic resolution operator:
  eta = eigenvalues[0] * eps
  cutoff = eta / 5000
  For each eigenvector v_i with eigenvalue lambda_i > cutoff:
    projection = dot(weights_inv * x_in, v_i)
    x_out += [lambda_i / (lambda_i + eta)] * projection * v_i
  x_out *= weights
"""

from __future__ import annotations

import numpy as np

from srts.model_data import ModelData

DEFAULT_EPS = {12: 40e-4, 20: 35e-4, 40: 20e-4}


def apply_resolution_matrix(
    model_coeffs: np.ndarray,
    model_data: ModelData,
    eps: float | None = None,
) -> np.ndarray:
    """Apply the tomographic resolution matrix to a reparameterized model.

    Fully vectorized: the eigenvector loop is replaced by matrix operations.

    Args:
        model_coeffs: Reparameterized model, shape (ndp, natd), Fortran convention.
        model_data: ModelData with eigenvectors, eigenvalues, weights.
        eps: Damping parameter. If None, uses DEFAULT_EPS[model_data.lmax].

    Returns:
        Filtered model coefficients, same shape as input.
    """
    if eps is None:
        eps = DEFAULT_EPS[model_data.lmax]

    ndp = model_coeffs.shape[0]
    natd = model_coeffs.shape[1]
    lenatd = ndp * natd

    # Flatten input
    x_in = model_coeffs.ravel()  # (lenatd,)

    # Weights
    if model_data.ismth == 1:
        twts = model_data.weights.ravel()  # (lenatd,)
        twtsinv = 1.0 / twts
    else:
        twts = np.ones(lenatd)
        twtsinv = np.ones(lenatd)

    # Damping
    eigenvalues = model_data.eigenvalues
    eigenvectors = model_data.eigenvectors  # (neig, lenatd)
    eta = eigenvalues[0] * eps
    cutoff = eta / 5000.0

    # Select active eigenvectors
    active = eigenvalues > cutoff
    V = eigenvectors[active]       # (n_active, lenatd)
    lam = eigenvalues[active]      # (n_active,)

    # Vectorized filtering
    #   projections = V @ (twtsinv * x_in)   → (n_active,)
    #   filters = lam / (lam + eta)          → (n_active,)
    #   x_out = V^T @ (filters * projections) → (lenatd,)
    #   x_out *= twts
    projections = V @ (twtsinv * x_in)
    filtered_proj = (lam / (lam + eta)) * projections
    x_out = V.T @ filtered_proj

    if model_data.ismth == 1:
        x_out *= twts

    # Crustal reweighting
    if model_data.icrust == 1:
        x_out[:natd] *= 1000.0

    return x_out.reshape(ndp, natd)


def extract_sp_from_spt(
    filtered: np.ndarray,
    model_data: ModelData,
) -> tuple[np.ndarray, int, int]:
    """Extract the S-wave (SP) parameter from multi-parameter filtered output.

    Equivalent to spt2sph extracting the .SP.sph component.

    Returns:
        (coefficients, ndmn, ndmx) where coefficients has shape (ndp, natd).
    """
    # iparsw[0] == 1 means SP is active; ipardps[0] gives its depth range
    ipardps = model_data.ipardps  # shape (2, mp1)
    ndmn = int(ipardps[0, 0])
    ndmx = int(ipardps[1, 0])
    ndp = ndmx - ndmn + 1

    # SP is always the first parameter block in the .spt ordering
    natd = (model_data.lmax + 1) ** 2
    return filtered[:ndp], ndmn, ndmx
