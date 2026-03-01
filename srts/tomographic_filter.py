"""Class-based tomographic filter wrapping the resolution matrix operator."""

from __future__ import annotations

import numpy as np

from srts.coeffs import cilm_stack_to_internal, internal_to_cilm_stack
from srts.filtering import DEFAULT_EPS, apply_resolution_matrix, extract_sp_from_spt
from srts.model_data import ModelData, load_model_data


class TomographicFilter:
    """Resolution matrix operator for SxRTS tomographic filtering.

    Accepts models in the 21-spline basis as cilm arrays, applies the
    tomographic resolution matrix R, and returns filtered models.
    """

    def __init__(self, degree: int, eps: float | None = None):
        if degree not in (12, 20, 40):
            raise ValueError(f"degree must be 12, 20, or 40, got {degree}")
        self._degree = degree
        self._model_data = load_model_data(degree)
        self._eps = eps if eps is not None else DEFAULT_EPS[degree]

    @property
    def degree(self) -> int:
        return self._degree

    @property
    def eps(self) -> float:
        return self._eps

    @property
    def lmax(self) -> int:
        return self._model_data.lmax

    @property
    def reference_model(self) -> np.ndarray:
        """The published SxRTS model as cilm arrays.

        Returns:
            shape (ndp, 2, lmax+1, lmax+1).
        """
        return internal_to_cilm_stack(
            self._model_data.reference_coefficients, self._model_data.lmax
        )

    def filter(self, coeffs: np.ndarray) -> np.ndarray:
        """Apply the tomographic resolution matrix to a reparameterized model.

        Args:
            coeffs: Model in spline basis, shape (21, 2, lmax+1, lmax+1).

        Returns:
            Filtered model, shape (ndp, 2, lmax+1, lmax+1).
        """
        flat = cilm_stack_to_internal(coeffs)
        filtered = apply_resolution_matrix(flat, self._model_data, self._eps)
        filt_coeffs, _, _ = extract_sp_from_spt(filtered, self._model_data)
        return internal_to_cilm_stack(filt_coeffs, self._model_data.lmax)


def S40RTS(eps: float | None = None) -> TomographicFilter:
    """Create a TomographicFilter for S40RTS."""
    return TomographicFilter(40, eps)


def S20RTS(eps: float | None = None) -> TomographicFilter:
    """Create a TomographicFilter for S20RTS."""
    return TomographicFilter(20, eps)


def S12RTS(eps: float | None = None) -> TomographicFilter:
    """Create a TomographicFilter for S12RTS."""
    return TomographicFilter(12, eps)
