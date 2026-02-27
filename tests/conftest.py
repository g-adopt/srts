import os
import pathlib

import pytest

TESTS_DIR = pathlib.Path(__file__).parent
REPO_ROOT = TESTS_DIR.parent

_LOCAL_FALLBACK = pathlib.Path(
    "/Users/sghelichkhani/Workplace/tomographic_filtering/srts/geodyn"
)


@pytest.fixture(scope="session")
def fortran_reference_dir():
    """Locate Fortran reference data, trying multiple sources.

    Resolution order:
      1. <repo_root>/test-data/fortran-reference/  (CI download location)
      2. SRTS_FORTRAN_GEODYN environment variable
      3. Hardcoded local dev path
    """
    ci_path = REPO_ROOT / "test-data" / "fortran-reference"
    if ci_path.is_dir():
        return ci_path

    env = os.environ.get("SRTS_FORTRAN_GEODYN")
    if env:
        p = pathlib.Path(env)
        if p.is_dir():
            return p

    if _LOCAL_FALLBACK.is_dir():
        return _LOCAL_FALLBACK

    pytest.skip("Fortran reference data not available")
