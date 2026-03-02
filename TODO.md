# Performance Optimisation TODO

Profiled with `python -m cProfile -o output.prof examples/benchmark.py` across three grid sizes
(1k, 10k, 65k points, 65 layers, lmax=40). The findings below are ordered by expected impact.

---

## 1. Batch `build_atd` across layers in `expand_batch` — biggest win

`expand_batch` currently calls the `build_atd` closure once per layer inside a Python loop. Each call
itself loops over latitude groups, accumulating `At @ values[mask]` as a GEMV for each group. For 65
layers and 181 unique latitudes (largest grid) this amounts to roughly 11,700 individual matrix-vector
products. The entire batch could instead be done as one GEMM per latitude group:
`At @ values_batch[:, mask].T` produces the full `(leny, nlayers)` ATD matrix in one shot. BLAS DGEMM
has substantially better arithmetic intensity than DGEMV because it reuses rows of `At` across all 65
layer columns simultaneously. The solve step (`V.T @ ATD` and `V @ (ATD / lam)`) would similarly
collapse from 65 GEMVs to two GEMMs. This is the clearest vectorisation opportunity in the codebase and
the one most likely to give a meaningful wall-time reduction.

Relevant code: `expansion.py`, `precompute_expansion` closure `build_atd`, `SphericalHarmonicExpansion.expand_batch`.

---

## 2. Cache-friendly matrix product in `apply_resolution_matrix`

The resolution matrix `V` is `(16000, 35301)` — about 4.5 GB — and all eigenvectors are active for
S40RTS, so there is no truncation. The two products `V @ x` and `V.T @ y` together read the full
matrix twice per call, making the function memory-bandwidth-bound at around 0.43 s/call. The forward
product `V @ x` is fine because V is C-contiguous (row-major) and rows are read sequentially. The
reverse product `V.T @ y` is the problem: with C-contiguous V, accessing columns is cache-unfriendly.
Replacing `V.T @ filtered_proj` with `np.dot(filtered_proj, V)` achieves the identical result but
reads V in row order, which is more cache-coherent and should reduce effective memory latency.

A more substantial gain would require changing the algorithm — truncating the eigenbasis at load time,
or using a randomised SVD — but that is a bigger undertaking.

Relevant code: `filtering.py`, `apply_resolution_matrix`.

---

## 3. Eliminate redundant format conversions in `reparameterize`

`DepthParameterization.reparameterize` calls `project_layer` once per input layer. Each `project_layer`
call converts the input cilm array to flat format, does the spline projection, then calls
`internal_to_cilm_stack` which loops over all 21 spline-depth entries converting each one back to cilm
individually. For 65 layers this produces 65 × 21 = 1,365 redundant round-trips between flat and cilm
representations, accounting for most of the 4,353 calls to `fortran_flat_to_shcoeffs` visible in the
profile. The fix is to restructure `reparameterize` to work entirely in flat format: convert the full
input stack once at entry, accumulate the spline projection in flat format, and convert to cilm a
single time at the end. The total time lost here is only about 40 ms, so this is an architectural
improvement more than a bottleneck fix, but it will compound with the batch expansion change above.

Relevant code: `parameterization.py`, `DepthParameterization.reparameterize`, `project_layer`.

---

## 4. Investigate double `eigh` call per grid in `precompute_expansion`

`scipy.linalg.eigh` appears 6 times in the profile for 3 benchmark grids, suggesting something triggers
it twice per grid rather than once. The eigendecomposition of the 1681×1681 ATA matrix takes roughly
0.37 s each, so an accidental double-call would waste about 1.1 s total. Worth tracing whether the
`SphericalHarmonicExpansion` constructor is somehow reaching `precompute_expansion` via two paths, or
whether there is a second call site being exercised at startup. If it turns out to be a genuine
double-call, fixing it is essentially free performance.

Even a single call per grid cannot be avoided with the current least-squares formulation, but switching
from `scipy.linalg.eigh` to `numpy.linalg.eigh` occasionally has lower overhead for full
decompositions and is worth a quick comparison.

Relevant code: `expansion.py`, `precompute_expansion`.

---

## 5. Minor: replace `np.outer` with broadcasting in the precompute loop

Inside the latitude loop in `precompute_expansion`, `np.outer(m_val, phi)` computes a
`(1681, n_lons_at_lat)` matrix. Replacing it with `m_val[:, np.newaxis] * phi` avoids the function
call dispatch overhead and is occasionally faster for this shape. The total time attributed to
`np.outer` in the profile is 0.08 s across all three grids, so the absolute gain is small, but the
change is a one-liner.

Relevant code: `expansion.py`, `precompute_expansion`, latitude loop.
