# Numba-Accelerated FCI Hamiltonian Assembly

**Date:** 2026-03-12
**Status:** Complete
**Version:** v0.9.38

---

## Summary

Replaced pure-Python loops in `DirectCISolver.assemble_hamiltonian()` with
Numba `@njit`-compiled kernels. Assembly time for LiH nmax=3 (367,290 SDs)
dropped from **975s to 1.3s** — a **738x speedup**.

---

## Before/After Comparison (LiH nmax=3, R=3.015 bohr)

| Phase | Before | After | Speedup |
|:---|---:|---:|---:|
| Construction (lattice + ERI) | 3.3s | 0.4s | - |
| DirectCI precompute | 0.4s | 0.2s | 2x |
| **H assembly** | **975s** | **1.3s** | **738x** |
| Eigensolve (eigsh) | 5.1s | 5.0s | - |
| **Total FCI** | **984s** | **7.8s** | **126x** |

### Scaling across systems

| System | N_SD | Python (s) | Numba (s) | Speedup |
|:---|---:|---:|---:|---:|
| He nmax=3 | 378 | 0.015 | 0.000 | ~15x |
| Li nmax=3 | 3,276 | 0.31 | 0.007 | 44x |
| Be nmax=3 | 20,475 | 3.6 | 0.09 | 40x |
| Li nmax=5 | 215,820 | 116s | ~2s | ~58x |
| Be nmax=4 | 487,635 | 357s | ~8s | ~45x |
| **LiH nmax=3** | **367,290** | **975s** | **1.3s** | **738x** |

The LiH speedup is larger than atomic systems because LiH has sparser ERIs
(0.55% density), so the Numba kernel skips more empty target pairs without
the Python interpreter overhead.

### Impact on PES scans

| Scenario | Per calc | Per CP point (3 calcs) | 10-point PES |
|:---|---:|---:|---:|
| Before | 984s | 49 min | 8.2 hours |
| **After** | **7.8s** | **23s** | **4 min** |

---

## Implementation

### Files

| File | Role |
|:---|:---|
| `geovac/direct_ci_numba.py` | Numba JIT kernels (new) |
| `geovac/direct_ci.py` | Integration: auto-detects Numba, falls back to Python |
| `tests/test_direct_ci_numba.py` | 8 correctness + speedup tests (new) |

### Architecture

1. **SD basis → 2D int32 array**: `sd_basis_to_array()` converts the list of
   sorted tuples to a contiguous `(n_sd, n_el)` NumPy array.

2. **Combinatorial number system**: `_sd_to_index()` converts any sorted
   occupation tuple to its combinatorial (colexicographic) index via
   `Σ C(occ[i], i+1)`. A precomputed `colex_to_lex` permutation array maps
   this to the sequential index used by `itertools.combinations`.

3. **Spatial target flattening**: The `spatial_targets` dict is converted to
   flat arrays (`offsets`, `lengths`, `values`) for O(1) lookup in Numba.
   Forward and reverse keys are pre-merged to eliminate the swap loop.

4. **Monolithic serial kernel**: `_assemble_kernel()` computes diagonal +
   singles + doubles in a single pass over all SDs. Serial `@njit` alone
   gives the full speedup — no `parallel=True` needed for the target.

5. **Overflow protection**: Output arrays are pre-allocated with a generous
   estimate (50 per SD). If the kernel hits the limit, it returns early and
   the caller retries with 2x allocation.

### Key design decisions

- **Serial @njit, not parallel**: Serial JIT gives 100x+ over CPython for
  tight numerical loops. The parallel two-pass approach (count → allocate →
  fill) was implemented but unnecessary for the target speedup. It can be
  enabled later if needed.

- **Combinatorial indexing instead of dict**: Numba can't use Python dicts.
  The combinatorial number system gives O(1) index computation from any
  sorted SD tuple. The colex→lex permutation array handles the ordering
  mismatch with `itertools.combinations`.

- **Pre-merged spatial targets**: The Python code checked both `(a,b)` and
  `(b,a)` keys in a swap loop. This caused double-counting in the Numba
  kernel. Pre-merging in `prepare_spatial_targets()` eliminates the issue.

---

## Correctness

- **Matrix-level**: `max |H_python - H_numba| < 1e-12` for all tested systems
- **Energy-level**: `|E_python - E_numba| < 1e-10 Ha` for all tested systems
- **NNZ match**: Identical sparsity patterns (e.g., 8,461,158 for LiH nmax=3)
- **Test suite**: 8/8 tests pass (`tests/test_direct_ci_numba.py`)
- **Existing tests**: All `tests/test_direct_ci.py` tests pass through the
  Numba path automatically
- **Topological integrity**: 18/18 symbolic proofs pass (unchanged)

---

## Limitations

- First-time JIT compilation takes ~2.5s (cached thereafter)
- Numba `cache=True` requires write access to `__pycache__`
- Falls back gracefully to Python if Numba is not installed
- Memory: pre-allocates output arrays at ~50 entries per SD (can be tuned)

---

## Dependencies

- `numba >= 0.60` (tested with 0.64.0)
- `llvmlite` (installed automatically with numba)
- No changes to `scipy`, `numpy`, or other dependencies
