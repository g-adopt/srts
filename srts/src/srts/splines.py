"""21-knot cardinal cubic spline depth basis for SxRTS models.

Reimplements splhsetup/splh/rspln/rsple from the Fortran code, fully vectorized.

splh(ind, xd) uses reversed indexing: ind=0 → CMB (xd=-1), ind=20 → Moho (xd=+1).
"""

from __future__ import annotations

import numpy as np

SPLINE_KNOTS = np.array([
    -1.00000, -0.78631, -0.59207, -0.41550, -0.25499,
    -0.10909,  0.02353,  0.14409,  0.25367,  0.35329,
     0.44384,  0.52615,  0.60097,  0.66899,  0.73081,
     0.78701,  0.83810,  0.88454,  0.92675,  0.96512,
     1.00000,
])

NKNOTS = 21
RCMB = 3480.0
RMOHO = 6346.619
RMOHO_SPHEXP = 6346.0
REARTH = 6371.0


def depth_to_xd(depth_km: float | np.ndarray, rmoho: float = RMOHO) -> float | np.ndarray:
    """Convert depth (km) to normalized coordinate xd in [-1, 1]."""
    r = REARTH - np.asarray(depth_km, dtype=np.float64)
    return -1.0 + 2.0 * (r - RCMB) / (rmoho - RCMB)


def xd_to_depth(xd: float | np.ndarray, rmoho: float = RMOHO) -> float | np.ndarray:
    """Convert normalized coordinate xd to depth (km)."""
    r = RCMB + (np.asarray(xd, dtype=np.float64) + 1.0) / 2.0 * (rmoho - RCMB)
    return REARTH - r


def _rspln(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Compute cubic spline coefficients. Direct translation of Fortran RSPLN.

    Uses 3-point end derivatives (NOT natural splines).
    Returns q of shape (3, n).
    """
    n = len(x)
    q = np.zeros((3, n))
    f = np.zeros((3, n))

    if n < 3:
        if n == 2:
            q[0, 0] = (y[1] - y[0]) / (x[1] - x[0])
        return q

    j1 = 0
    j2 = n - 3

    a0 = 0.0
    h = x[j1 + 1] - x[j1]
    h2 = x[j1 + 2] - x[j1]
    y0 = h * h2 * (h2 - h)
    b0 = (y[j1] * (h * h - h2 * h2) + y[j1 + 1] * h2 * h2 - y[j1 + 2] * h * h) / y0
    b1 = b0

    for i in range(j1, j2 + 1):
        h = x[i + 1] - x[i]
        y0 = y[i + 1] - y[i]
        h2 = h * h
        ha = h - a0
        h2a = h - 2.0 * a0
        h3a = 2.0 * h - 3.0 * a0
        h2b = h2 * b0
        q[0, i] = h2 / ha
        q[1, i] = -ha / (h2a * h2)
        q[2, i] = -h * h2a / h3a
        f[0, i] = (y0 - h * b0) / (h * ha)
        f[1, i] = (h2b - y0 * (2.0 * h - a0)) / (h * h2 * h2a)
        f[2, i] = -(h2b - 3.0 * y0 * ha) / (h * h3a)
        a0 = q[2, i]
        b0 = f[2, i]

    i = j2 + 1
    h = x[i + 1] - x[i]
    y0 = y[i + 1] - y[i]
    h2 = h * h
    ha = h - a0
    h2a = h * ha
    h2b = h2 * b0 - y0 * (2.0 * h - a0)
    q[0, i] = h2 / ha
    f[0, i] = (y0 - h * b0) / h2a
    ha = x[j2] - x[i + 1]
    y0 = -h * ha * (ha + h)
    ha_sq = ha * ha
    y0 = (y[i + 1] * (h2 - ha_sq) + y[i] * ha_sq - y[j2] * h2) / y0
    q[2, i] = (y0 * h2a + h2b) / (h * h2 * (h - 2.0 * a0))
    q[1, i] = f[0, i] - q[0, i] * q[2, i]

    for j in range(j1, j2 + 1):
        k = i - 1
        q[0, i] = f[2, k] - q[2, k] * q[1, i]
        q[2, k] = f[1, k] - q[1, k] * q[0, i]
        q[1, k] = f[0, k] - q[0, k] * q[2, k]
        i = k

    q[0, i] = b1
    q[:, n - 1] = 0.0
    return q


class SplineBasis:
    """21-knot cardinal cubic spline basis functions.

    Precomputes coefficients so that evaluate_all is a pure vectorized operation.
    """

    def __init__(self) -> None:
        self._knots = SPLINE_KNOTS.copy()
        self._n = NKNOTS
        # Precompute spline coefficients for each cardinal basis function
        qq0 = np.eye(NKNOTS)
        qq = np.zeros((3, NKNOTS, NKNOTS))
        for i in range(NKNOTS):
            qq[:, :, i] = _rspln(self._knots, qq0[:, i])

        # Reorder so that index 0 = splh(0) = deepest (CMB), index 20 = shallowest (Moho)
        # splh(ind) uses column NKNOTS-1-ind, so we reverse the column order
        self._y = qq0[:, ::-1].copy()       # (NKNOTS, NKNOTS) — y values per basis
        self._q = qq[:, :, ::-1].copy()      # (3, NKNOTS, NKNOTS) — spline coeffs per basis

    def evaluate_all(self, xd: np.ndarray) -> np.ndarray:
        """Evaluate all 21 basis functions at an array of xd values.

        Fully vectorized: no Python loops over points or basis functions.

        Args:
            xd: shape (npts,)

        Returns:
            shape (21, npts). Row i = splh(i, xd).
        """
        xd = np.atleast_1d(np.asarray(xd, dtype=np.float64))
        npts = len(xd)
        result = np.zeros((self._n, npts))

        valid = (xd >= -1.0) & (xd <= 1.0)
        if not np.any(valid):
            return result

        xv = xd[valid]  # (nv,)

        # Find interval index for each point: searchsorted on increasing knots
        # idx[j] = interval containing xv[j], clipped to [0, NKNOTS-2]
        idx = np.searchsorted(self._knots, xv, side='right') - 1
        np.clip(idx, 0, self._n - 2, out=idx)

        h = xv - self._knots[idx]  # (nv,)

        # Evaluate all 21 basis functions at once using broadcasting
        # self._y: (NKNOTS, 21) — y[knot, basis]
        # self._q: (3, NKNOTS, 21) — q[coeff, knot, basis]
        # We need y[idx, :] and q[:, idx, :] for each point
        y_at_idx = self._y[idx, :]            # (nv, 21)
        q0 = self._q[0, idx, :]               # (nv, 21)
        q1 = self._q[1, idx, :]               # (nv, 21)
        q2 = self._q[2, idx, :]               # (nv, 21)

        # Horner's method: y + h*(q0 + h*(q1 + h*q2))
        h = h[:, np.newaxis]  # (nv, 1) for broadcasting
        vals = y_at_idx + h * (q0 + h * (q1 + h * q2))  # (nv, 21)

        result[:, valid] = vals.T  # (21, nv)
        return result

    def evaluate(self, ind: int, xd: np.ndarray) -> np.ndarray:
        """Evaluate a single basis function ind at array of xd values."""
        return self.evaluate_all(xd)[ind]


_spline_basis: SplineBasis | None = None


def get_spline_basis() -> SplineBasis:
    """Get the cached spline basis singleton."""
    global _spline_basis
    if _spline_basis is None:
        _spline_basis = SplineBasis()
    return _spline_basis


def evaluate_spline_basis(xd: np.ndarray) -> np.ndarray:
    """Evaluate all 21 basis functions at xd. Returns shape (21, npts)."""
    return get_spline_basis().evaluate_all(xd)
