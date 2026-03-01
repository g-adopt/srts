"""Load resolution operator data (eigenvectors, weights) for SxRTS models."""

from __future__ import annotations

import pathlib
import sys
from dataclasses import dataclass

import h5py
import numpy as np

S3_BASE_URL = "https://gadopt.syd1.digitaloceanspaces.com/srts/"


@dataclass
class ModelData:
    """Holds the resolution operator data for one SxRTS model."""

    lmax: int
    ndep: int
    eigenvalues: np.ndarray       # (neig,)
    eigenvectors: np.ndarray      # (neig, lenatd)
    weights: np.ndarray           # (ndep, natd) from .smthp_21
    ismth: int
    icrust: int
    mp1: int
    iparsw: np.ndarray
    parwts: np.ndarray
    ipardps: np.ndarray
    reference_coefficients: np.ndarray  # from SxRTS.sph
    reference_lmax: int
    reference_ndmn: int
    reference_ndmx: int

    @property
    def natd(self) -> int:
        return (self.lmax + 1) ** 2

    @property
    def lenatd(self) -> int:
        return self.ndep * self.natd


def _data_dir() -> pathlib.Path:
    """Return the path to the shipped data directory."""
    return pathlib.Path(__file__).parent / "data"


def _download_file(url: str, destination: pathlib.Path) -> None:
    """Download a file from a public S3 URL.

    Tries boto3 with unsigned config first (better for large files),
    falls back to urllib for environments without boto3.
    """
    size_info = ""
    try:
        import urllib.request
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req) as resp:
            content_length = resp.headers.get("Content-Length")
            if content_length:
                size_mb = int(content_length) / (1024 * 1024)
                size_info = f" ({size_mb:.0f} MB)"
    except Exception:
        pass

    filename = destination.name
    print(f"Downloading {filename}{size_info} from S3...", file=sys.stderr)

    tmp = destination.with_suffix(".download")
    try:
        try:
            import boto3
            from botocore import UNSIGNED
            from botocore.config import Config

            s3 = boto3.client(
                "s3",
                endpoint_url="https://syd1.digitaloceanspaces.com",
                config=Config(signature_version=UNSIGNED),
            )
            s3.download_file("gadopt", f"srts/{filename}", str(tmp))
        except ImportError:
            import urllib.request
            urllib.request.urlretrieve(url, str(tmp))

        tmp.rename(destination)
        print(f"Downloaded {filename} to {destination}", file=sys.stderr)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise


def _ensure_model_data(degree: int) -> pathlib.Path:
    """Return path to the HDF5 file, downloading from S3 if absent."""
    h5path = _data_dir() / f"S{degree}RTS.h5"
    if h5path.exists():
        return h5path
    h5path.parent.mkdir(parents=True, exist_ok=True)
    _download_file(S3_BASE_URL + f"S{degree}RTS.h5", h5path)
    return h5path


def load_model_data(degree: int) -> ModelData:
    """Load HDF5 data for S12RTS (degree=12), S20RTS (20), or S40RTS (40)."""
    if degree not in (12, 20, 40):
        raise ValueError(f"degree must be 12, 20, or 40, got {degree}")

    h5path = _ensure_model_data(degree)

    with h5py.File(h5path, "r") as f:
        meta = f["metadata"]
        lmax = int(meta["lmax"][()])
        ndep = int(meta["ndep"][()])
        icrust = int(meta["icrust"][()])
        ismth = int(meta["ismth"][()])
        mp1 = int(meta["mp1"][()])
        iparsw = meta["iparsw"][:]
        parwts = meta["parwts"][:]
        ipardps = meta["ipardps"][:]

        eig = f["eigenvectors"]
        eigenvalues = eig["eigenvalues"][:]
        eigenvectors = eig["vectors"][:]

        wts = f["weights"]
        weights = wts["twts"][:]

        ref = f["reference_model"]
        ref_coeffs = ref["coefficients"][:]
        ref_lmax = int(ref["lmax"][()])
        ref_ndmn = int(ref["ndmn"][()])
        ref_ndmx = int(ref["ndmx"][()])

    return ModelData(
        lmax=lmax,
        ndep=ndep,
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        weights=weights,
        ismth=ismth,
        icrust=icrust,
        mp1=mp1,
        iparsw=iparsw,
        parwts=parwts,
        ipardps=ipardps,
        reference_coefficients=ref_coeffs,
        reference_lmax=ref_lmax,
        reference_ndmn=ref_ndmn,
        reference_ndmx=ref_ndmx,
    )
