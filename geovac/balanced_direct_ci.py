"""
Matrix-free (Davidson) direct CI for the balanced-coupled Hamiltonian.
======================================================================

Specialised to the N_alpha = N_beta = 2 sector (4 electrons) -- exactly the
balanced-LiH sector that geovac.coupled_composition.coupled_fci_energy
diagonalizes in the A/B/C and chemistry-error sprints -- but written so the
sigma vector is pure BLAS-3.

Representation
--------------
Compact amplitude  C[Ia, Ib],  Ia = alpha pair (i<j), Ib = beta pair (k<l),
strings enumerated in itertools.combinations order (identical ordering to
coupled_fci_energy).

H = H_a (x) 1  +  1 (x) H_b  +  V_ab  +  E_core
  H_a  = one-body(alpha) + two-body(alpha,alpha)   dense n_a x n_a, built once
  V_ab = sum_{pqrs} (pq|rs) E^a_pq E^b_rs          one (M^2 x M^2) dgemm

For N=2 the same-spin block has the closed form

  H_a[(i<j),(k<l)] = d_jl h[i,k] - d_il h[j,k] - d_jk h[i,l] + d_ik h[j,l]
                     + s * ( (ik|jl) - (il|jk) )

with s = +1 (physically correct, and the convention of the corrected
geovac.coupled_composition._double_excitation_phase) or, in `faithful` mode,
s = -1 whenever {i,j} and {k,l} are disjoint -- the HISTORICAL pre-correction
sign convention of that function, kept only to reproduce bug-era banked
numbers (the balanced max_n series recorded before the sign fix; the spurious
shift is ~+4 mHa over-binding on balanced LiH).  See
tests/test_balanced_direct_ci.py::test_library_same_spin_double_phase_is_correct
for the element-by-element pin against a brute-force Fock-space FCI.

The alpha-beta term uses the redundant antisymmetric tensor

  Chat[i,j,k,l] = +C[(i<j),(k<l)],  antisymmetric in (i,j) and in (k,l)

so that

  (V_ab C)^[i,j,k,l] = A_a A_b [ sum_{q,s} (i q | k s) Chat[q,j,s,l] ]

which is a single (M^2 x M^2) @ (M^2 x M^2) matrix product, i.e. M^6 flops of
BLAS-3 rather than an explicit sparse-matrix assembly.
"""
from __future__ import annotations

import itertools
import time
from typing import Any, Dict, Optional, Tuple

import numpy as np


# --------------------------------------------------------------------------
def pair_lists(M: int) -> Tuple[np.ndarray, np.ndarray]:
    """(i, j) orbital arrays for the C(M,2) pairs in itertools.combinations order."""
    pairs = np.array(list(itertools.combinations(range(M), 2)), dtype=np.int64)
    return pairs[:, 0].copy(), pairs[:, 1].copy()


# --------------------------------------------------------------------------
def build_same_spin_H(
    h1: np.ndarray,
    eri: np.ndarray,
    M: int,
    faithful: bool = False,
    block: int = 1024,
) -> np.ndarray:
    """Dense n_a x n_a alpha-only Hamiltonian (one-body + same-spin two-body)."""
    I, J = pair_lists(M)
    na = I.size
    H = np.empty((na, na), dtype=np.float64)
    erif = eri.reshape(-1)
    M2, M3 = M * M, M * M * M
    for lo in range(0, na, block):
        hi = min(lo + block, na)
        Ia = I[lo:hi][:, None]
        Ja = J[lo:hi][:, None]
        Ib = I[None, :]
        Jb = J[None, :]
        blk = (np.where(Ja == Jb, h1[Ia, Ib], 0.0)
               - np.where(Ia == Jb, h1[Ja, Ib], 0.0)
               - np.where(Ja == Ib, h1[Ia, Jb], 0.0)
               + np.where(Ia == Ib, h1[Ja, Jb], 0.0))
        two = (erif[(Ia * M3 + Ib * M2 + Ja * M + Jb).ravel()]
               - erif[(Ia * M3 + Jb * M2 + Ja * M + Ib).ravel()]).reshape(hi - lo, na)
        if faithful:
            disjoint = (Ia != Ib) & (Ia != Jb) & (Ja != Ib) & (Ja != Jb)
            two = np.where(disjoint, -two, two)
        H[lo:hi] = blk + two
    return H


# --------------------------------------------------------------------------
class DirectCI4e:
    """Matrix-free CI for N_alpha = N_beta = 2 (4 electrons)."""

    def __init__(self, h1: np.ndarray, eri: np.ndarray, e_core: float,
                 faithful: bool = False, verbose: bool = False):
        M = h1.shape[0]
        assert eri.shape == (M, M, M, M)
        self.M = M
        self.e_core = float(e_core)
        self.faithful = bool(faithful)
        self.verbose = verbose
        self.I, self.J = pair_lists(M)
        self.na = self.nb = int(self.I.size)
        self.ndet = self.na * self.nb

        t0 = time.perf_counter()
        self.Ha = build_same_spin_H(h1, eri, M, faithful=faithful)
        t_ha = time.perf_counter() - t0

        t0 = time.perf_counter()
        self.G2 = np.ascontiguousarray(eri.transpose(0, 2, 1, 3)).reshape(M * M, M * M)
        t_g2 = time.perf_counter() - t0

        t0 = time.perf_counter()
        M2 = M * M
        ia, ja = self.I[:, None], self.J[:, None]
        ib, jb = self.I[None, :], self.J[None, :]
        self._idx = [
            ((ia * M + ib) * M2 + (ja * M + jb)).astype(np.int32).ravel(),
            ((ja * M + ib) * M2 + (ia * M + jb)).astype(np.int32).ravel(),
            ((ia * M + jb) * M2 + (ja * M + ib)).astype(np.int32).ravel(),
            ((ja * M + jb) * M2 + (ia * M + ib)).astype(np.int32).ravel(),
        ]
        self._sgn = (1.0, -1.0, -1.0, 1.0)
        t_idx = time.perf_counter() - t0

        self._C2 = np.zeros(M2 * M2, dtype=np.float64)
        self._T2 = np.empty((M2, M2), dtype=np.float64)
        self._buf = np.empty(self.ndet, dtype=np.float64)

        d_a = np.diag(self.Ha)
        ar = np.arange(M)
        Jmat = np.ascontiguousarray(eri[ar, ar][:, ar, ar])          # J[p,q] = (pp|qq)
        W = Jmat[self.I] + Jmat[self.J]                              # (na, M)
        self.diag = (d_a[:, None] + d_a[None, :]
                     + W[:, self.I] + W[:, self.J] + self.e_core)

        if verbose:
            mem = (self.Ha.nbytes + self.G2.nbytes + self._C2.nbytes
                   + self._T2.nbytes + sum(a.nbytes for a in self._idx)
                   + self.diag.nbytes + self._buf.nbytes)
            print(f"  [setup] M={M} n_a={self.na} n_det={self.ndet:,}  "
                  f"Ha {t_ha:.1f}s  G2 {t_g2:.1f}s  idx {t_idx:.1f}s  "
                  f"persistent mem ~{mem/1e9:.2f} GB", flush=True)
        self.n_sigma = 0
        self.t_sigma = 0.0

    # -------------------------------------------------------------- sigma
    def sigma(self, C: np.ndarray) -> np.ndarray:
        """C : (na, nb) -> H C."""
        t0 = time.perf_counter()
        M2 = self.M * self.M
        out = self.Ha @ C
        out += C @ self.Ha
        out += self.e_core * C
        Cf = C.reshape(-1)
        C2 = self._C2
        C2[:] = 0.0
        for idx, s in zip(self._idx, self._sgn):
            np.multiply(Cf, s, out=self._buf)
            C2[idx] = self._buf
        np.dot(self.G2, C2.reshape(M2, M2), out=self._T2)
        T2f = self._T2.reshape(-1)
        outf = out.reshape(-1)
        for idx, s in zip(self._idx, self._sgn):
            np.take(T2f, idx, out=self._buf)
            if s > 0:
                outf += self._buf
            else:
                outf -= self._buf
        self.n_sigma += 1
        self.t_sigma += time.perf_counter() - t0
        return out

    # ------------------------------------------------------------ davidson
    def ground_state(self, tol: float = 1e-6, max_iter: int = 300,
                     max_sub: int = 10, n_guess: int = 4,
                     verbose: bool = True,
                     time_budget_s: Optional[float] = None) -> Dict[str, Any]:
        nd = self.ndet
        dflat = self.diag.reshape(-1)
        V = np.empty((max_sub, nd), dtype=np.float64)
        W = np.empty((max_sub, nd), dtype=np.float64)

        order = np.argsort(dflat)[:n_guess]
        k = 0
        for o in order:
            v = np.zeros(nd)
            v[o] = 1.0
            if k:
                v -= (V[:k] @ v) @ V[:k]
            nrm = float(np.linalg.norm(v))
            if nrm > 1e-8:
                V[k] = v / nrm
                k += 1
        nW = 0
        theta = float(dflat[order[0]])
        rn = np.inf
        it = -1
        t_start = time.perf_counter()
        hist = []
        for it in range(max_iter):
            while nW < k:
                W[nW] = self.sigma(V[nW].reshape(self.na, self.nb)).reshape(-1)
                nW += 1
            S = V[:k] @ W[:k].T
            S = 0.5 * (S + S.T)
            w, y = np.linalg.eigh(S)
            theta = float(w[0])
            yv = y[:, 0]
            x = yv @ V[:k]
            r = yv @ W[:k] - theta * x
            rn = float(np.linalg.norm(r))
            hist.append((it, theta, rn))
            if verbose:
                print(f"    davidson it={it:3d} k={k:2d} E={theta:+.12f} |r|={rn:.3e} "
                      f"[{time.perf_counter()-t_start:.0f}s]", flush=True)
            if rn < tol:
                break
            if time_budget_s is not None and (time.perf_counter() - t_start) > time_budget_s:
                if verbose:
                    print("    davidson: TIME BUDGET exceeded -- returning current Ritz value",
                          flush=True)
                break
            if k == max_sub:
                V[0] = x / np.linalg.norm(x)
                W[0] = self.sigma(V[0].reshape(self.na, self.nb)).reshape(-1)
                k, nW = 1, 1
            den = dflat - theta
            den = np.where(np.abs(den) < 1e-6, 1e-6, den)
            t = r / den
            for _ in range(2):
                t -= (V[:k] @ t) @ V[:k]
            nrm = float(np.linalg.norm(t))
            if nrm < 1e-12:
                break
            V[k] = t / nrm
            k += 1
        return {
            'E': theta, 'residual': rn, 'n_iter': it + 1,
            'n_sigma': self.n_sigma, 't_sigma_total': self.t_sigma,
            't_sigma_avg': self.t_sigma / max(self.n_sigma, 1),
            'n_det': nd, 'history': hist,
            'wall_s': time.perf_counter() - t_start,
            'converged': bool(rn < tol),
        }


# --------------------------------------------------------------------------
def solve_from_ham(ham: Dict[str, Any], n_electrons: int = 4,
                   faithful: bool = False, verbose_setup: bool = False,
                   **kw) -> Dict[str, Any]:
    assert n_electrons == 4, "DirectCI4e is specialised to N_alpha = N_beta = 2"
    ci = DirectCI4e(ham['h1'], ham['eri'], ham['nuclear_repulsion'],
                    faithful=faithful, verbose=verbose_setup)
    return ci.ground_state(**kw)
