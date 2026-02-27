#!/usr/bin/env python3
"""Convert Fortran binary model files (.evc, .smthp_21, .sph) to HDF5.

Reads unformatted sequential Fortran binary files directly in Python using
numpy/struct, avoiding the need for a separate Fortran dump program.

Usage:
    python convert_to_hdf5.py /path/to/utils/S40RTS /path/to/output/S40RTS.h5
"""

from __future__ import annotations

import argparse
import pathlib
import struct
import sys

import h5py
import numpy as np


def _read_fortran_record(f) -> bytes:
    """Read one record from a Fortran sequential unformatted file.

    Fortran writes: [4-byte length] [data] [4-byte length]
    """
    header = f.read(4)
    if len(header) < 4:
        raise EOFError("End of file reached")
    (nbytes,) = struct.unpack("i", header)
    data = f.read(nbytes)
    trailer = f.read(4)
    (nbytes2,) = struct.unpack("i", trailer)
    if nbytes != nbytes2:
        raise ValueError(f"Record length mismatch: {nbytes} vs {nbytes2}")
    return data


def read_evc(filepath: pathlib.Path) -> dict:
    """Read a .evc file (eigenvectors from tomographic inversion).

    Record 1: lmaxh, numatd2, ndep, icrust, idensi1, idum1, ismth  (7 x int32)
    Record 2: mp1, iparsw(1:mp1), parwts(1:mp1), ipardps(1:2,1:mp1)  (mixed)
    Records 3..N: eigv(i), evc(1:lenatd)  (float64)
    """
    with open(filepath, "rb") as f:
        # Record 1: 7 integers
        data = _read_fortran_record(f)
        ints = struct.unpack("7i", data)
        lmaxh, numatd2, ndep, icrust, idensi1, idum1, ismth = ints

        # Record 2: mp1 (int), iparsw (mp1 ints), parwts (mp1 float32), ipardps (2*mp1 ints)
        # Note: parwts is default Fortran REAL (4 bytes), not DOUBLE PRECISION
        data = _read_fortran_record(f)
        offset = 0
        (mp1,) = struct.unpack_from("i", data, offset)
        offset += 4
        iparsw = np.array(struct.unpack_from(f"{mp1}i", data, offset))
        offset += mp1 * 4
        parwts = np.array(struct.unpack_from(f"{mp1}f", data, offset), dtype=np.float64)
        offset += mp1 * 4
        ipardps_flat = np.array(struct.unpack_from(f"{2 * mp1}i", data, offset))
        # Fortran stores ((ipardps(i,j),i=1,2),j=1,mp1), i.e. column-major pairs
        ipardps = ipardps_flat.reshape(mp1, 2).T  # (2, mp1)

        # Records 3..N: eigenvalues and eigenvectors
        natd = numatd2
        lenatd = ndep * natd
        eigenvalues = []
        eigenvectors = []

        while True:
            try:
                data = _read_fortran_record(f)
            except EOFError:
                break
            if len(data) == 0:
                break
            # Each record: 1 double (eigenvalue) + lenatd doubles (eigenvector)
            expected_size = (1 + lenatd) * 8
            if len(data) < expected_size:
                break
            vals = np.frombuffer(data, dtype=np.float64)
            eigenvalues.append(vals[0])
            eigenvectors.append(vals[1 : 1 + lenatd])

    eigenvalues = np.array(eigenvalues)
    eigenvectors = np.array(eigenvectors)

    return {
        "lmaxh": lmaxh,
        "numatd2": numatd2,
        "ndep": ndep,
        "icrust": icrust,
        "ismth": ismth,
        "mp1": mp1,
        "iparsw": iparsw,
        "parwts": parwts,
        "ipardps": ipardps,
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors,
    }


def read_smthp(filepath: pathlib.Path) -> dict:
    """Read a .smthp_21 file (smoothness/data weights).

    Record 1: lmaxw, nsmn, nsmx, ndepw (int32), etaz, etah, etai (float32), iderh, iderv (int32)
    Records 2..ndepw+1: twts(1:natd) per depth level (float32)
    """
    with open(filepath, "rb") as f:
        data = _read_fortran_record(f)
        # 4 ints + 3 floats + 2 ints (all default Fortran types = 4 bytes each)
        offset = 0
        lmaxw, nsmn, nsmx, ndepw = struct.unpack_from("4i", data, offset)
        offset += 16
        etaz, etah, etai = struct.unpack_from("3f", data, offset)
        offset += 12
        iderh, iderv = struct.unpack_from("2i", data, offset)

        natd = (lmaxw + 1) ** 2
        twts = np.zeros((ndepw, natd), dtype=np.float64)
        for i in range(ndepw):
            data = _read_fortran_record(f)
            twts[i] = np.frombuffer(data, dtype=np.float32)[:natd].astype(np.float64)

    return {
        "lmaxw": lmaxw,
        "nsmn": nsmn,
        "nsmx": nsmx,
        "ndepw": ndepw,
        "etaz": etaz,
        "etah": etah,
        "etai": etai,
        "iderh": iderh,
        "iderv": iderv,
        "twts": twts,
    }


def read_reference_sph(filepath: pathlib.Path) -> dict:
    """Read a reference .sph file (ASCII format)."""
    from srts.io import read_sph
    return read_sph(filepath)


def convert_model(model_dir: pathlib.Path, output_path: pathlib.Path) -> None:
    """Convert all binary files for one SxRTS model to HDF5."""
    model_name = model_dir.name  # e.g. "S40RTS"

    evc_path = model_dir / f"{model_name}.evc"
    smthp_path = model_dir / f"{model_name}.smthp_21"
    sph_path = model_dir / f"{model_name}.sph"

    print(f"Reading {evc_path} ...")
    evc = read_evc(evc_path)

    print(f"Reading {smthp_path} ...")
    smthp = read_smthp(smthp_path)

    print(f"Reading {sph_path} ...")
    ref = read_reference_sph(sph_path)

    print(f"Writing {output_path} ...")
    with h5py.File(output_path, "w") as f:
        meta = f.create_group("metadata")
        meta.create_dataset("lmax", data=evc["lmaxh"])
        meta.create_dataset("ndep", data=evc["ndep"])
        meta.create_dataset("icrust", data=evc["icrust"])
        meta.create_dataset("ismth", data=evc["ismth"])
        meta.create_dataset("mp1", data=evc["mp1"])
        meta.create_dataset("iparsw", data=evc["iparsw"])
        meta.create_dataset("parwts", data=evc["parwts"])
        meta.create_dataset("ipardps", data=evc["ipardps"])

        eig = f.create_group("eigenvectors")
        eig.create_dataset("eigenvalues", data=evc["eigenvalues"])
        eig.create_dataset("vectors", data=evc["eigenvectors"])

        wts = f.create_group("weights")
        wts.create_dataset("lmaxw", data=smthp["lmaxw"])
        wts.create_dataset("ndepw", data=smthp["ndepw"])
        wts.create_dataset("twts", data=smthp["twts"])

        refg = f.create_group("reference_model")
        refg.create_dataset("lmax", data=ref["lmax"])
        refg.create_dataset("ndmn", data=ref["ndmn"])
        refg.create_dataset("ndmx", data=ref["ndmx"])
        refg.create_dataset("coefficients", data=ref["coefficients"])

    print(f"Done: {output_path} "
          f"({evc['eigenvalues'].shape[0]} eigenvectors, "
          f"lmax={evc['lmaxh']}, ndep={evc['ndep']})")


def main():
    parser = argparse.ArgumentParser(description="Convert Fortran binary SxRTS files to HDF5")
    parser.add_argument("model_dir", type=pathlib.Path, help="Directory containing .evc, .smthp_21, .sph")
    parser.add_argument("output", type=pathlib.Path, help="Output .h5 file path")
    args = parser.parse_args()

    convert_model(args.model_dir, args.output)


if __name__ == "__main__":
    main()
