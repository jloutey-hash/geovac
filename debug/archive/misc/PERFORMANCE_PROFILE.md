# LiH FCI Performance Profile

**Date:** 2026-03-12
**System:** LiH, R=3.015 bohr, nmax=3, 4 electrons
**Status:** Diagnostic complete

---

## Executive Summary

**The Hamiltonian assembly dominates at 99.1% of total runtime.** Within assembly, singles excitations take 73% and doubles take 25%. The eigensolve (eigsh) is essentially free at 0.5%. The graph Laplacian's O(V) advantage is completely lost in the Python-loop FCI assembly.

| Phase | Time (s) | % of Total |
|:---|---:|---:|
| Construction (lattice + ERI) | 3.3 | 0.3% |
| DirectCI precompute | 0.4 | 0.0% |
| **H assembly — diagonal** | **9.5** | **1.0%** |
| **H assembly — singles** | **716.5** | **72.8%** |
| **H assembly — doubles** | **246.8** | **25.1%** |
| Eigensolve (eigsh) | 5.1 | 0.5% |
| **Total** | **984.3** | **100%** |

---

## Scaling: nmax=2 vs nmax=3

| Metric | nmax=2 | nmax=3 | Ratio |
|:---|---:|---:|---:|
| Spatial orbitals | 10 | 28 | 2.8x |
| Spin-orbitals | 20 | 56 | 2.8x |
| Slater determinants | 4,845 | 367,290 | 75.8x |
| CI NNZ | 46,509 | 8,461,158 | 181.9x |
| Assembly time | 0.35s | 975s | 2,765x |
| Total time | 0.48s | 984s | 2,050x |

**Scaling exponent:** Assembly scales as N_SD^2.3 (log(2765)/log(75.8) = 1.83 in wall time, but the inner loops are O(N_SD x n_el x n_sp) for singles). The 2,765x blowup for 75.8x more SDs indicates **super-linear scaling in the Python loop**, consistent with O(N_SD x n_el x n_sp) = O(367k x 4 x 56) = 82M iterations for singles alone.

---

## Sparsity Analysis

### H1 (one-electron Hamiltonian)

| Property | Value |
|:---|---:|
| Shape | 28 x 28 |
| NNZ | 108 / 784 |
| Density | 13.8% |
| Diagonal NNZ | 28 |
| Off-diagonal NNZ | 80 |
| Memory (dense) | 6.1 KB |

**Verdict:** Tiny matrix. Dense storage is correct (6 KB). No optimization needed here.

### ERI (two-electron integrals)

| Property | Value |
|:---|---:|
| Stored entries (dict) | 3,394 |
| Max possible (28^4) | 614,656 |
| Density | 0.55% |
| Dense 4D memory | 4.69 MB |
| Dict memory (est) | 0.31 MB |
| Same-atom A (Li) | 1,492 |
| Same-atom B (H) | 1,492 |
| Cross-atom | 410 |

**Threshold analysis — all ERIs are significant:**

| Threshold | Entries above | % of stored |
|:---|---:|---:|
| > 1e-4 | 3,334 | 98.2% |
| > 1e-6 | 3,394 | 100% |
| > 1e-8 | 3,394 | 100% |

**Verdict:** ERIs are already extremely sparse (0.55% of dense 4D). All stored entries are significant (no negligible tails to prune). The dense 4D array (4.69 MB) is the right choice for fast element access — dict lookups were the original bottleneck (v0.9.6).

### CI Hamiltonian

| Property | Value |
|:---|---:|
| Shape | 367,290 x 367,290 |
| NNZ | 8,461,158 |
| Full size | 134.9 billion |
| Density | 0.0063% |
| Avg NNZ/row | 23.0 |
| CSR memory | 98.2 MB |

**Verdict:** Extremely sparse. The CI matrix is well-suited for sparse storage. eigsh handles it efficiently (231 matvecs x 9.76 ms = 2.3s total matvec time).

---

## Eigensolve Performance

| Metric | Value |
|:---|---:|
| Algorithm | ARPACK (eigsh, SA) |
| Matvec count | 231 |
| Total matvec time | 2.26s |
| Avg matvec time | 9.76 ms |
| Total eigsh time | 5.08s |
| Overhead (ARPACK) | 2.82s |

**Verdict:** eigsh is already fast. The sparse CI matrix (23 nnz/row avg) makes each matvec cheap. This is NOT the bottleneck.

---

## Memory Breakdown (nmax=3)

| Component | Memory |
|:---|---:|
| H1 dense | 0.006 MB |
| ERI 4D array | 4.69 MB |
| CI matrix (CSR) | 98.2 MB |
| SD basis + index | ~50 MB (est) |
| Peak tracked | 568.3 MB |

The ~568 MB peak is dominated by temporary arrays during COO→CSR conversion and the eigenvector storage.

---

## Root Cause Analysis

### Why is assembly so slow?

The `DirectCISolver.assemble_hamiltonian()` method ([direct_ci.py:161](geovac/direct_ci.py#L161)) has three nested Python loops:

**Singles (73% of time):**
```
for I in range(367,290):           # each SD
    for kp in range(4):            # each occupied orbital
        for r in range(56):        # each virtual orbital
            # Slater-Condon matrix element (4 more inner iterations)
```
= 367,290 x 4 x 56 = **82.3 million** outer iterations, each with an inner loop over 4 occupied orbitals for the Coulomb/exchange terms. Total: **~330 million** Python-level operations.

**Doubles (25% of time):**
```
for I in range(367,290):           # each SD
    for kp in range(4):            # occupied pair
        for kq in range(kp+1, 4):  # ...
            for (sp_lo, sp_hi) in targets:  # ERI targets (~5 avg)
```
The `spatial_targets` dict prunes the virtual space effectively (2,878 target pairs), but the SD loop is still 367k iterations in pure Python.

### The fundamental problem

The assembly is **O(N_SD x n_el x n_sp)** in *Python-interpreted loops*. Each iteration does:
- frozenset membership test (`r in occ_I`)
- Dense array indexing (`H1[sp_p, sp_r]`, `eri[sp_p, sp_q, sp_r, sp_q]`)
- List/tuple manipulation for new SD construction
- Dict lookup for SD→index mapping

These are fast individually (~100ns each) but 330M iterations x ~3us/iteration = ~990s.

---

## Recommendations (ranked by expected impact)

### 1. Cythonize the assembly loops (est. 50-100x speedup)

The innermost loops in `assemble_hamiltonian()` are pure Python with NumPy array indexing — a textbook case for Cython or Numba. Converting the singles loop to C-level iteration would reduce 717s to ~7-15s.

**Approach:** `@numba.njit` on the singles and doubles loops, with SD basis as a 2D int array and pre-allocated output arrays. The dense H1/ERI arrays are already NumPy — no conversion needed.

**Estimated speedup:** 50-100x on assembly → total time drops from 984s to ~15-25s.

### 2. Matrix-free (sigma-vector) direct CI (est. 100x+ for PES scans)

Instead of explicitly building H as a sparse matrix, compute H|c> on-the-fly during each eigsh matvec. This is the standard approach in production CI codes (MOLPRO, PySCF).

**Current:** Build H (975s) + 231 matvecs (2.3s) = 977s
**Matrix-free:** 231 sigma-vector builds x ~4s each = ~924s (no gain without Cython)

This approach only wins if combined with Cython/C for the sigma-vector, where each sigma build would take ~0.04s → 231 x 0.04 = 9.2s total. It also eliminates the 98 MB CSR storage.

### 3. Spin-adapted determinants (est. 2-4x reduction in N_SD)

The current basis enumerates all C(56,4) = 367,290 SDs without spin symmetry. For a singlet ground state (S=0, M_S=0), only determinants with 2 alpha + 2 beta electrons contribute: C(28,2)^2 = 142,884 SDs — a **2.6x reduction**.

Further, configuration state functions (CSFs) would reduce this to ~70k, giving ~5x fewer basis functions and ~25x faster assembly.

### 4. Screening via Schwarz inequality (est. 2-5x speedup on doubles)

For doubles, pre-compute Q(ij) = sqrt(|(ij|ij)|) and skip pairs where Q(ij)*Q(rs) < threshold. With 98% of ERIs above 1e-4, this won't help much at nmax=3, but becomes critical at nmax=4+.

### 5. Parallelization (est. N_cores x speedup)

The SD loop is embarrassingly parallel — each row of H is independent. A simple `multiprocessing.Pool` or `concurrent.futures` could give 4-8x on a typical machine, no algorithmic changes needed.

---

## Impact Estimates for a BSSE-Corrected PES Point

Each CP-corrected PES point requires **3 FCI calculations** (molecule + 2 ghost atoms).

| Scenario | Time per calc | Time per PES point | 10-point PES |
|:---|---:|---:|---:|
| **Current** | 984s | 49 min | 8.2 hours |
| + Cython assembly (50x) | ~20s | 1 min | 10 min |
| + Spin adaptation (2.6x) | ~8s | 24s | 4 min |
| + Both | ~4s | 12s | 2 min |

---

## The Core Question Answered

> Does the graph Laplacian's O(V) advantage survive the molecular/FCI pipeline?

**Yes and no.** The graph Laplacian construction and ERI computation are fast (3.3s). The eigsh solve is fast (5s) because the CI Hamiltonian is extremely sparse (23 nnz/row). The O(V) advantage is preserved in the *physics* — but the *assembly* of the CI Hamiltonian from Slater-Condon rules is a standard computational chemistry bottleneck that has nothing to do with the graph Laplacian.

The ERIs are already 99.45% sparse (only 3,394 of 614,656 possible entries). This is a consequence of the graph topology — only graph-adjacent orbitals have non-zero coupling. **The sparsity is there; it's the Python loop overhead that kills performance.**

The fix is implementation, not physics: Cython/Numba for the assembly loop, or a matrix-free sigma-vector approach. Both are standard in production quantum chemistry codes and would bring nmax=3 LiH from 16 minutes to under 30 seconds.
