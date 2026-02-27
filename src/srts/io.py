"""I/O for the legacy Fortran file formats (.raw, .sph, .spt, layer .dat)."""

from __future__ import annotations

import pathlib
import re

import numpy as np


def read_raw(filepath: str | pathlib.Path) -> tuple[int, np.ndarray]:
    """Read a .raw file. Returns (lmax, flat coefficient array).

    Format: line 1 = lmax (and optionally other header fields),
    remaining lines = (lmax+1)^2 coefficients, 5 per line in e16.8 format.
    """
    filepath = pathlib.Path(filepath)
    with open(filepath) as f:
        header = f.readline().split()
        lmax = int(header[0])
        n = (lmax + 1) ** 2
        values = []
        for line in f:
            values.extend(float(x) for x in line.split())
    coeffs = np.array(values[:n], dtype=np.float64)
    return lmax, coeffs


def write_raw(filepath: str | pathlib.Path, lmax: int, coeffs: np.ndarray) -> None:
    """Write a .raw file in Fortran format."""
    filepath = pathlib.Path(filepath)
    n = (lmax + 1) ** 2
    assert coeffs.shape == (n,)
    with open(filepath, "w") as f:
        f.write(f"{lmax:3d}\n")
        for i in range(0, n, 5):
            chunk = coeffs[i : i + 5]
            f.write("".join(f"{v:16.8e}" for v in chunk) + "\n")


def _parse_sph_header(line: str) -> dict:
    """Parse a .sph header line following the rsphhead.f logic.

    The header has 4 whitespace-separated groups:
      1. lmax (integer)
      2. binary mask of which degrees are present (up to 100 digits of 0/1)
      3. ndep (integer) — total number of depth slots
      4. binary mask of which depth slots are present (up to 24 digits of 0/1)

    The Fortran wsphhead writes this as:
      format(10x, i4, x, 100i1)  for lmax + degree mask
      format(i4, x, 24i1)        for ndep + depth mask
    concatenated into a single line.

    rsphhead.f uses strgrep to split on whitespace groups.
    """
    # Use the same logic as rsphhead: split into space-separated groups
    groups = line.split()
    if len(groups) < 4:
        raise ValueError(f"Expected 4 groups in sph header, got {len(groups)}: {line!r}")

    lmax = int(groups[0])
    lask_str = groups[1]
    ndep = int(groups[2])
    mask_str = groups[3]

    # Parse degree mask (lask): string of 0/1 digits
    lask = [int(c) for c in lask_str]

    # Parse depth mask: string of 0/1 digits
    mask = [int(c) for c in mask_str]

    # Find actual lmax from lask
    actual_lmax = 0
    for i, v in enumerate(lask):
        if v == 1:
            actual_lmax = i

    # Find ndmn and ndmx from mask (1-indexed in Fortran)
    ndmn = None
    ndmx = None
    for i, v in enumerate(mask):
        if v in (1, 2):
            if ndmn is None:
                ndmn = i + 1  # 1-indexed
            ndmx = i + 1

    if ndmn is None:
        raise ValueError(f"No active depth levels in mask: {mask_str}")

    ndp = ndmx - ndmn + 1

    return {
        "lmax": actual_lmax,
        "ndep": ndep,
        "ndmn": ndmn,
        "ndmx": ndmx,
        "ndp": ndp,
        "lask": lask,
        "mask": mask,
    }


def read_sph(filepath: str | pathlib.Path) -> dict:
    """Read a .sph file. Returns dict with header info and coefficients.

    Keys: 'lmax', 'ndmn', 'ndmx', 'ndp', 'ndep', 'lask', 'mask',
          'coefficients' (2D array: ndp x natd where natd = (lmax+1)^2)
    """
    filepath = pathlib.Path(filepath)
    with open(filepath) as f:
        header_line = f.readline()
        header = _parse_sph_header(header_line)

        lmax = header["lmax"]
        ndp = header["ndp"]
        natd = (lmax + 1) ** 2

        # The Fortran format 11e12.4 writes up to 11 values per line,
        # wrapping to the next line for degrees with >11 terms.
        # Read all remaining data as a flat stream of floats.
        all_values = []
        for line in f:
            all_values.extend(float(x) for x in line.split())

        total_expected = ndp * natd
        coefficients = np.array(all_values[:total_expected], dtype=np.float64).reshape(ndp, natd)

    result = dict(header)
    result["coefficients"] = coefficients
    return result


def _parse_fortran_floats(line: str, expected: int) -> np.ndarray:
    """Parse Fortran-formatted floats from a line.

    Handles both space-separated and fixed-width (11e12.4) formats.
    """
    # Try space-separated first
    parts = line.split()
    if len(parts) == expected:
        return np.array([float(x) for x in parts])

    # Try fixed-width e12.4 format (12 chars per value)
    vals = []
    pos = 0
    width = 12
    while pos + width <= len(line) and len(vals) < expected:
        token = line[pos : pos + width].strip()
        if token:
            vals.append(float(token))
        pos += width

    if len(vals) == expected:
        return np.array(vals)

    # Fallback: just parse all floats we can find
    vals = [float(x) for x in line.split()]
    return np.array(vals[:expected])


def write_sph(
    filepath: str | pathlib.Path,
    lmax: int,
    ndmn: int,
    ndmx: int,
    coeffs: np.ndarray,
) -> None:
    """Write a .sph file in Fortran format.

    Args:
        filepath: Output path.
        lmax: Maximum spherical harmonic degree.
        ndmn: First depth index (1-indexed, typically 4).
        ndmx: Last depth index (1-indexed, typically 24).
        coeffs: Coefficient array, shape (ndp, natd) where ndp = ndmx - ndmn + 1.
    """
    filepath = pathlib.Path(filepath)
    ndp = ndmx - ndmn + 1
    natd = (lmax + 1) ** 2
    assert coeffs.shape == (ndp, natd), f"Expected ({ndp}, {natd}), got {coeffs.shape}"

    header = _make_sph_header(lmax, ndmn, ndmx)

    with open(filepath, "w") as f:
        f.write(f" {header}\n")
        for idep in range(ndp):
            idx = 0
            for l in range(lmax + 1):
                nterms = 2 * l + 1
                vals = coeffs[idep, idx : idx + nterms]
                f.write(_format_e12_4(vals) + "\n")
                idx += nterms


def _make_sph_header(lmax: int, ndmn: int, ndmx: int) -> str:
    """Construct the .sph header string matching wsphhead.f output."""
    # Degree mask: 1 for l=0..lmax, 0 beyond
    lask = "1" * (lmax + 1)

    # Depth mask: 24 slots, 1 for ndmn..ndmx
    ndep_total = 24
    mask = ["0"] * ndep_total
    for i in range(ndmn - 1, ndmx):
        mask[i] = "1"
    mask_str = "".join(mask)

    # Format matching wsphhead: format(10x, i4, x, 100i1) || format(i4, x, 24i1)
    str1 = f"{'':>10s}{lmax:4d} {lask}"
    str2 = f"{ndep_total:4d} {mask_str}"
    return str1 + str2


def _format_e12_4(vals: np.ndarray) -> str:
    """Format values in Fortran 11e12.4 format."""
    parts = []
    for v in vals:
        parts.append(f"{v:12.4e}")
    return "".join(parts)


def read_spt(filepath: str | pathlib.Path) -> dict:
    """Read a .spt file (multi-parameter model).

    Returns dict with 'lmax', 'mp1', 'iparsw', 'ipardps', 'coefficients'.
    """
    filepath = pathlib.Path(filepath)
    with open(filepath) as f:
        # Header: format(i2, x, i2, x, 10i1)
        header_line = f.readline()
        lmax = int(header_line[:2])
        mp1 = int(header_line[3:5])
        iparsw = [int(header_line[6 + i]) for i in range(mp1)]

        # Read depth ranges for active parameters
        ipardps = {}
        npartot = 0
        for i in range(mp1):
            if iparsw[i] == 1:
                line = f.readline()
                parts = line.split()
                idp1, idp2 = int(parts[0]), int(parts[1])
                ipardps[i] = (idp1, idp2)
                npartot += idp2 - idp1 + 1

        # Read coefficients
        natd = (lmax + 1) ** 2
        coefficients = np.zeros((npartot, natd), dtype=np.float64)
        for idep in range(npartot):
            idx = 0
            for l in range(lmax + 1):
                nterms = 2 * l + 1
                line = f.readline()
                vals = _parse_fortran_floats(line, nterms)
                coefficients[idep, idx : idx + nterms] = vals
                idx += nterms

    return {
        "lmax": lmax,
        "mp1": mp1,
        "iparsw": iparsw,
        "ipardps": ipardps,
        "coefficients": coefficients,
    }


def read_layer_data(
    model_dir: str | pathlib.Path,
    name: str,
    first_layer: int,
    last_layer: int,
) -> dict:
    """Read geodyn layer files and depth_layers.dat.

    Returns dict with 'depths', 'lon', 'lat', 'values' (list of arrays), 'nlayers'.
    """
    model_dir = pathlib.Path(model_dir)

    # Read depth boundaries
    depth_file = model_dir / "depth_layers.dat"
    depths = []
    with open(depth_file) as f:
        for line in f:
            line = line.strip()
            if line:
                depths.append(float(line))
    depths = np.array(depths)

    # Read layer data files
    lon = None
    lat = None
    layer_values = []

    for iz in range(first_layer, last_layer):
        num = f"{iz:03d}"
        layer_file = model_dir / f"{name}.layer.{num}.dat"
        data = np.loadtxt(layer_file)
        if lon is None:
            lon = data[:, 0]
            lat = data[:, 1]
        layer_values.append(data[:, 2])

    return {
        "depths": depths,
        "lon": lon,
        "lat": lat,
        "values": layer_values,
        "nlayers": last_layer - first_layer,
    }
