"""Atomic transcorrelated (TC / xTC) Coulomb-Sturmian engine (Paper 14, sec:tc).

Tracked port of the sprint drivers ``debug/xtc_poc_li.py`` (the s-only PoC engine)
and ``debug/xtc_quantum_cost_measure.py`` (the non-Hermitian quantum-cost
diagnostics), promoted so the Paper 14 transcorrelated-operator claims are backed by
a regression suite (``tests/test_paper14_tc_nonhermitian.py``) rather than transient
``debug/`` scripts.  The physics is a faithful copy of ``debug/xtc_poc_li.py`` --- no
re-derivation, no fitted parameters --- with type hints, a dataclass system container,
and the operator-cost diagnostics collected into one place.

Transcorrelation
----------------
With a Jastrow factor ``tau = sum_{i<j} u(r_ij)`` and the Slater geminal
``u(r) = -(1/2 gamma) e^{-gamma r}`` (``u'(0) = 1/2`` -- the singlet Kato cusp), the
similarity-transformed Hamiltonian is

    Htilde = e^{-tau} H e^{tau} = H + D + K + L3 ,

where

    D  = -sum_i [ (1/2) lap_i tau + (1/2) (grad_i tau)^2 (j=k part) ]   (2-body, Hermitian),
    K  = -sum_i (grad_i tau) . grad_i                                    (2-body, NON-Hermitian),
    L3 = -(1/2) sum_i sum_{j!=i, k!=i, j!=k} grad_i u(r_ij) . grad_i u(r_ik)  (3-body).

``xTC`` (Christlmaier-Kats-Alavi, JCP 159 014113 (2023)) contracts the genuine 3-body
``L3`` down to an effective 2-body operator ``v2`` (+ 1-body ``v1`` + scalar ``v0``)
against the reference 1-RDM (diagonal occupation of the single-determinant reference).

Structural facts backed here (see the module functions and the test file):

  * the xTC 3-body -> 2-body contraction ``v2`` is Hermitian to machine precision,
    and so is ``D`` (``asym_w``);  the ONLY non-Hermitian term is the generic 2-body
    convective ``K`` (``asym_K``), present in ANY transcorrelated Hamiltonian;
  * the full non-Hermitian effective operator has an entirely real spectrum, a real
    ground state, benign eigenvector conditioning ``kappa_V = O(1)``, and an LCU
    1-norm within ~1.13x of plain;
  * ``K`` roughly doubles the Jordan-Wigner Pauli count (117 -> 249 for the He/Li
    s-only ns=3 systems); symmetrizing ``(Htilde + Htilde^dag)/2`` restores the plain
    count (117) and an LCU 1-norm below plain, but is an accuracy false economy (it
    overshoots below the exact energy -- see :func:`tc_energies`).

Everything is built on ONE radial grid with the Coulomb-Sturmian s-radials
``R_{n0}(r; n, k) = N e^{-k r} L^1_{n-1}(2 k r)`` (shared decay ``k``),
Loewdin-orthonormalized, then FCI via particle-number-projected second-quantized
operator application.  Non-Hermitian pieces are diagonalized with ``scipy.linalg.eig``
and the real ground eigenvalue is selected.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations, permutations
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.linalg import eig, eigh
from scipy.special import eval_genlaguerre

# exact non-relativistic ground-state energies (Ha), for accuracy gating
EXACT_NR: Dict[str, float] = {"He": -2.90372, "Li": -7.47806}

Det = Tuple[int, ...]


# ----------------------------------------------------------------------------
# Coulomb-Sturmian s-radial R_{n0}(r; n, k) and derivative  (shared decay k)
# ----------------------------------------------------------------------------
def R_and_dR(n: int, r: np.ndarray, k: float) -> Tuple[np.ndarray, np.ndarray]:
    """Radial function ``R_{n0}`` and its derivative on the grid ``r`` (shared decay k)."""
    N = 2.0 * k ** 1.5 / n
    x = 2 * k * r
    L1 = eval_genlaguerre(n - 1, 1, x)
    R = N * np.exp(-k * r) * L1
    L2 = eval_genlaguerre(n - 2, 2, x) if n >= 2 else np.zeros_like(r)
    dR = N * np.exp(-k * r) * (-k * L1 - 2 * k * L2)   # dL^1_{n-1}/dx = -L^2_{n-2}
    return R, dR


def make_grid(k: float, Ng: int = 1400, r_max: Optional[float] = None
              ) -> Tuple[np.ndarray, np.ndarray]:
    """Quadratic radial grid + trapezoidal weights adapted to decay ``k``."""
    if r_max is None:
        r_max = 46.0 / k
    t = np.linspace(0.0, 1.0, Ng)
    r = r_max * t ** 2
    r[0] = 1e-9
    wr = np.zeros(Ng)
    wr[1:-1] = (r[2:] - r[:-2]) / 2.0
    wr[0] = (r[1] - r[0]) / 2.0
    wr[-1] = (r[-1] - r[-2]) / 2.0
    return r, wr


# ----------------------------------------------------------------------------
# Geminal
# ----------------------------------------------------------------------------
def up_of(r: np.ndarray, g: float) -> np.ndarray:
    """``u'(r) = (1/2) e^{-g r}`` for the Slater geminal."""
    return 0.5 * np.exp(-g * r)


def w_kernel(r: np.ndarray, g: float) -> np.ndarray:
    """Hermitian TC effective e-e kernel ``1/r + D`` (finite at the origin)."""
    r = np.asarray(r, float)
    small = r < 1e-8
    rr = np.where(small, 1.0, r)
    out = -np.expm1(-g * rr) / rr + (g / 2) * np.exp(-g * rr) - 0.25 * np.exp(-2 * g * rr)
    out = np.where(small, 1.5 * g - 0.25, out)
    return out


# ----------------------------------------------------------------------------
# (Ng x Ng) angular (L=0) / vertex kernels, integrated in u = r_12
# ----------------------------------------------------------------------------
# These were previously Gauss-Legendre in x = cos(theta_12).  That converges only as
# ~1/nx on any Coulomb-containing integrand, because 1/r_12 is singular at x = 1 when
# r_1 ~ r_2: measured 2.7e-3 relative error at the production nx = 160, 1.1e-3 at
# nx = 400, 3.6e-4 at nx = 1200.  Propagated to energies it over-binds by 20-70 uHa
# (He plain 39, He + geminal 23, Li 69) -- inside chemical accuracy, but 5-15% of the
# recorded R12-CI error figures.
#
# Substituting u = r_12 removes it exactly.  Since dx = -u du / (r_1 r_2),
#
#     (1/2) INT_-1^1 f(r_12) dx  =  (1/(2 r_1 r_2)) INT_{|r1-r2|}^{r1+r2} f(u) u du
#
# and the Jacobian's factor of u CANCELS the 1/r_12, leaving a smooth integrand.  The
# projected vertex kernels benefit the same way, using
#     r_1 - r_2 x = (r_1^2 - r_2^2 + u^2) / (2 r_1)
#     r_1 x - r_2 = (r_1^2 - r_2^2 - u^2) / (2 r_2)
# so their 1/r_12 cancels too.  Result: 8.3e-8 and nx-INDEPENDENT (converged by nx=40,
# where the x-route is still at 3.6e-4 with nx=1200) -- more accurate and cheaper.
# Verified by convergence: the old x-quadrature marches monotonically onto these values
# as nx grows (Li crosses through at nx=3600).  See debug/sprint_r12ci_gamma_transfer_memo.md.
def _u_nodes(r: np.ndarray, nx: int):
    """Per-(r_i, r_j) Gauss-Legendre data for the u = r_12 substitution."""
    R1, R2 = np.meshgrid(r, r, indexing="ij")
    lo = np.abs(R1 - R2)
    half = (R1 + R2 - lo) / 2.0
    ts, ws = np.polynomial.legendre.leggauss(nx)
    pref = 1.0 / (2.0 * R1 * R2)
    return R1, R2, lo, half, ts, ws, pref


def build_coul_kernel(r: np.ndarray, nx: int = 200) -> np.ndarray:
    """gamma-independent L=0 Coulomb kernel only (for fast plain-FCI k-scans).

    Exactly 1/max(r_i, r_j); computed via the u substitution so the same code path
    serves the general case.
    """
    R1, R2, lo, half, ts, ws, pref = _u_nodes(r, nx)
    acc = np.zeros_like(R1)
    for t, wt in zip(ts, ws):
        u = lo + half * (t + 1.0)
        acc += wt * half * (1.0 / np.maximum(u, 1e-30)) * u
    return pref * acc


def build_kernels(r: np.ndarray, g: float, nx: int = 200) -> Dict[str, np.ndarray]:
    """Coulomb, finite-TC ``w``, and convective-vertex (kA, kB) L=0 kernels."""
    R1, R2, lo, half, ts, ws, pref = _u_nodes(r, nx)
    Kcoul = np.zeros_like(R1)
    Kw = np.zeros_like(R1)
    KkA = np.zeros_like(R1)
    KkB = np.zeros_like(R1)
    dR2 = R1 * R1 - R2 * R2
    for t, wt in zip(ts, ws):
        u = lo + half * (t + 1.0)
        u = np.maximum(u, 1e-30)
        jac = wt * half
        up = 0.5 * np.exp(-g * u)
        Kcoul += jac * (1.0 / u) * u
        Kw += jac * w_kernel(u, g) * u
        # (r1 - r2 x)/r12 * u  ->  (r1^2 - r2^2 + u^2) / (2 r1)
        KkA += jac * up * (dR2 + u * u) / (2.0 * R1)
        # (r1 x - r2)/r12 * u  ->  (r1^2 - r2^2 - u^2) / (2 r2)
        KkB += jac * up * (dR2 - u * u) / (2.0 * R2)
    return dict(coul=pref * Kcoul, w=pref * Kw, kA=pref * KkA, kB=pref * KkB)


# ----------------------------------------------------------------------------
# One-body (grid), spatial, non-orthogonal Sturmian basis, s-only
# ----------------------------------------------------------------------------
def build_one_body(ns: int, r: np.ndarray, wr: np.ndarray, k: float, Z: float
                   ) -> Tuple[np.ndarray, np.ndarray, Dict[int, Tuple[np.ndarray, np.ndarray]], np.ndarray]:
    """Overlap ``S``, one-body ``h1`` (kinetic + nuclear), radial table, and ``r^2 dr``."""
    Rtab = {n: R_and_dR(n, r, k) for n in range(1, ns + 1)}
    W = r * r * wr
    S = np.zeros((ns, ns))
    h1 = np.zeros((ns, ns))
    for i in range(ns):
        Ri, dRi = Rtab[i + 1]
        for j in range(ns):
            Rj, dRj = Rtab[j + 1]
            S[i, j] = np.sum(Ri * Rj * W)
            T = 0.5 * np.sum(dRi * dRj * W)             # l=0 gradient form
            Vnuc = -Z * np.sum(Ri * Rj * r * wr)         # -Z/r : r^2 dr * (1/r) = r dr
            h1[i, j] = T + Vnuc
    return S, h1, Rtab, W


# ----------------------------------------------------------------------------
# Two-body spatial integrals <ij|op|kl> (electron1: i,k ; electron2: j,l), s-only
# ----------------------------------------------------------------------------
def two_body(ns: int, Rtab: Dict[int, Tuple[np.ndarray, np.ndarray]], W: np.ndarray,
             Kmats: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Coulomb, finite-TC ``w``, and convective ``K`` two-body tensors (physicist order)."""
    RR = {i: Rtab[i + 1][0] for i in range(ns)}
    dRR = {i: Rtab[i + 1][1] for i in range(ns)}
    D = {(i, k): RR[i] * RR[k] * W for i in range(ns) for k in range(ns)}      # density*W
    Dd = {(i, k): RR[i] * dRR[k] * W for i in range(ns) for k in range(ns)}    # R_i dR_k *W

    def eri(K: np.ndarray) -> np.ndarray:
        out = np.zeros((ns, ns, ns, ns))
        for i in range(ns):
            for j in range(ns):
                for kk in range(ns):
                    for l in range(ns):
                        out[i, j, kk, l] = D[(i, kk)] @ K @ D[(j, l)]
        return out

    eri_coul = eri(Kmats["coul"])
    eri_w = eri(Kmats["w"])
    # convective K2 = -<ij| u'(r) rhat.(grad1-grad2) |kl>
    eri_K = np.zeros((ns, ns, ns, ns))
    for i in range(ns):
        for j in range(ns):
            for kk in range(ns):
                for l in range(ns):
                    t1 = Dd[(i, kk)] @ Kmats["kA"] @ D[(j, l)]
                    t2 = D[(i, kk)] @ Kmats["kB"] @ Dd[(j, l)]
                    eri_K[i, j, kk, l] = -t1 + t2
    return eri_coul, eri_w, eri_K


# ----------------------------------------------------------------------------
# Three-body spatial integral V3[i,j,k; l,m,n]  (vertex particle1: i<->l)
# ----------------------------------------------------------------------------
def three_body(ns: int, Rtab: Dict[int, Tuple[np.ndarray, np.ndarray]], W: np.ndarray,
               Kmats: Dict[str, np.ndarray]) -> np.ndarray:
    """Genuine 3-body ``L3`` spatial tensor (shared-vertex vertex-kernel product)."""
    RR = {i: Rtab[i + 1][0] for i in range(ns)}
    kA = Kmats["kA"]
    Dvert = {(i, l): RR[i] * RR[l] * W for i in range(ns) for l in range(ns)}
    G = {(j, m): kA @ (RR[j] * RR[m] * W) for j in range(ns) for m in range(ns)}  # Ng-vector
    V3 = np.zeros((ns, ns, ns, ns, ns, ns))
    for i in range(ns):
        for l in range(ns):
            dv = Dvert[(i, l)]
            for j in range(ns):
                for m in range(ns):
                    Gjm = G[(j, m)]
                    dvG = dv * Gjm
                    for kk in range(ns):
                        for n in range(ns):
                            V3[i, j, kk, l, m, n] = dvG @ G[(kk, n)]
    return V3


# ----------------------------------------------------------------------------
# Loewdin orthonormalization + basis transforms
# ----------------------------------------------------------------------------
def lowdin(S: np.ndarray) -> np.ndarray:
    """Symmetric orthonormalizer ``X = S^{-1/2}``."""
    ev, U = eigh(S)
    ev = np.maximum(ev, 1e-14)
    return U @ np.diag(1.0 / np.sqrt(ev)) @ U.T


def transform_1(M: np.ndarray, X: np.ndarray) -> np.ndarray:
    return X @ M @ X


def transform_2(E: np.ndarray, X: np.ndarray) -> np.ndarray:
    t = np.einsum("pi,ijkl->pjkl", X, E)
    t = np.einsum("qj,pjkl->pqkl", X, t)
    t = np.einsum("rk,pqkl->pqrl", X, t)
    return np.einsum("sl,pqrl->pqrs", X, t)


def transform_3(V: np.ndarray, X: np.ndarray) -> np.ndarray:
    t = np.einsum("ai,ijklmn->ajklmn", X, V)
    t = np.einsum("bj,ajklmn->abklmn", X, t)
    t = np.einsum("ck,abklmn->abclmn", X, t)
    t = np.einsum("dl,abclmn->abcdmn", X, t)
    t = np.einsum("em,abcdmn->abcden", X, t)
    return np.einsum("fn,abcden->abcdef", X, t)


# ----------------------------------------------------------------------------
# Spin-orbital plumbing.  spin-orbital p = 2*spatial + spin
# ----------------------------------------------------------------------------
def spatial(p: int) -> int:
    return p >> 1


def spin(p: int) -> int:
    return p & 1


def apply_ops(det: Det, ops: List[Tuple[str, int]]) -> Optional[Tuple[int, Det]]:
    """Apply creation/annihilation ops (LEFT->RIGHT written; rightmost acts first)."""
    sign = 1
    d = det
    for kind, idx in reversed(ops):
        if kind == "a":
            if idx not in d:
                return None
            sign *= (-1) ** sum(1 for y in d if y < idx)
            d = tuple(y for y in d if y != idx)
        else:
            if idx in d:
                return None
            sign *= (-1) ** sum(1 for y in d if y < idx)
            d = tuple(sorted(d + (idx,)))
    return sign, d


def make_dets(nso: int, n_elec: int) -> Tuple[List[Det], Dict[Det, int]]:
    dets = list(combinations(range(nso), n_elec))
    return dets, {d: i for i, d in enumerate(dets)}


def h_spin(h1_spatial: np.ndarray, nso: int) -> np.ndarray:
    """Lift a spatial 1-body matrix to spin-orbitals (spin-diagonal)."""
    h = np.zeros((nso, nso))
    for p in range(nso):
        for q in range(nso):
            if spin(p) == spin(q):
                h[p, q] = h1_spatial[spatial(p), spatial(q)]
    return h


def asym_from_phys(eri_spatial: np.ndarray, nso: int) -> np.ndarray:
    """Antisymmetrized spin-orbital 2-body ``<pq||rs>`` from a physicist spatial eri."""
    a = np.zeros((nso, nso, nso, nso))
    for p in range(nso):
        for q in range(nso):
            for rr in range(nso):
                for s in range(nso):
                    v = 0.0
                    if spin(p) == spin(rr) and spin(q) == spin(s):
                        v += eri_spatial[spatial(p), spatial(q), spatial(rr), spatial(s)]
                    if spin(p) == spin(s) and spin(q) == spin(rr):
                        v -= eri_spatial[spatial(p), spatial(q), spatial(s), spatial(rr)]
                    a[p, q, rr, s] = v
    return a


# ----------------------------------------------------------------------------
# FCI matrix from 1-body h[p,q] and ANTISYMMETRIZED 2-body asym[p,q,r,s]=<pq||rs>
#   (convention-free operator application; works for non-Hermitian h/asym)
# ----------------------------------------------------------------------------
def build_H(dets: List[Det], didx: Dict[Det, int], h: np.ndarray, asym: np.ndarray,
            nso: int, v0: float = 0.0) -> np.ndarray:
    """FCI matrix ``sum h a+a + (1/4) sum <pq||rs> a+ a+ a a + v0`` (non-Hermitian OK)."""
    n = len(dets)
    H = np.zeros((n, n))
    hnz = [(p, q) for p in range(nso) for q in range(nso) if abs(h[p, q]) > 1e-14]
    for J, dJ in enumerate(dets):
        H[J, J] += v0
        for (p, q) in hnz:
            if q in dJ:
                res = apply_ops(dJ, [("c", p), ("a", q)])
                if res:
                    sgn, I = res
                    H[didx[I], J] += sgn * h[p, q]
        occ = dJ
        for rr in occ:
            for s in occ:
                if s == rr:
                    continue
                for p in range(nso):
                    for q in range(nso):
                        if q == p:
                            continue
                        val = asym[p, q, rr, s]
                        if abs(val) < 1e-14:
                            continue
                        res = apply_ops(dJ, [("c", p), ("c", q), ("a", s), ("a", rr)])
                        if res:
                            sgn, I = res
                            H[didx[I], J] += 0.25 * sgn * val
    return H


def build_H3(dets: List[Det], didx: Dict[Det, int], V3so: np.ndarray, nso: int,
             prefac: float = -0.5) -> np.ndarray:
    """Exact 3-body FCI matrix from ``V3so`` (for validation of the xTC contraction)."""
    n = len(dets)
    H = np.zeros((n, n))
    ns = V3so.shape[0]
    for J, dJ in enumerate(dets):
        occ = list(dJ)
        for (s, t, u) in permutations(occ, 3):
            ss, sspin = spatial(s), spin(s)
            ts, tspin = spatial(t), spin(t)
            us, uspin = spatial(u), spin(u)
            for pi in range(ns):
                p = 2 * pi + sspin
                for qi in range(ns):
                    q = 2 * qi + tspin
                    if q == p:
                        continue
                    for ri in range(ns):
                        rr = 2 * ri + uspin
                        if rr == p or rr == q:
                            continue
                        val = V3so[pi, qi, ri, ss, ts, us]
                        if abs(val) < 1e-14:
                            continue
                        res = apply_ops(dJ, [("c", p), ("c", q), ("c", rr),
                                             ("a", u), ("a", t), ("a", s)])
                        if res:
                            sgn, I = res
                            H[didx[I], J] += prefac * sgn * val
    return H


# ----------------------------------------------------------------------------
# xTC contraction: 3-body -> effective 0/1/2-body via reference 1-RDM (diag occ)
# ----------------------------------------------------------------------------
def xtc_contract(V3so: np.ndarray, nso: int, occ: Tuple[int, ...], prefac: float = -0.5
                 ) -> Tuple[np.ndarray, np.ndarray, float]:
    """Contract the 3-body ``L3`` to effective (v2, v1_bare, v0_bare) via a diagonal 1-RDM."""
    def T(p: int, q: int, rr: int, s: int, t: int, u: int) -> float:
        if spin(p) != spin(s) or spin(q) != spin(t) or spin(rr) != spin(u):
            return 0.0
        return prefac * V3so[spatial(p), spatial(q), spatial(rr),
                             spatial(s), spatial(t), spatial(u)]

    perms = list(permutations(range(3)))

    def sgn(perm: Tuple[int, ...]) -> int:
        s = 1
        for a in range(3):
            for b in range(a + 1, 3):
                if perm[a] > perm[b]:
                    s = -s
        return s

    def W(idx: Tuple[int, int, int, int, int, int]) -> float:  # idx = (p,q,r,s,t,u)
        bra = idx[:3]
        ket = idx[3:]
        tot = 0.0
        for pb in perms:
            sb = sgn(pb)
            b = (bra[pb[0]], bra[pb[1]], bra[pb[2]])
            for pk in perms:
                sk = sgn(pk)
                kk = (ket[pk[0]], ket[pk[1]], ket[pk[2]])
                v = T(b[0], b[1], b[2], kk[0], kk[1], kk[2])
                if v != 0.0:
                    tot += sb * sk * v
        return tot

    occ = list(occ)
    v2 = np.zeros((nso, nso, nso, nso))
    v1 = np.zeros((nso, nso))
    for p in range(nso):
        for q in range(nso):
            for s in range(nso):
                for t in range(nso):
                    acc = 0.0
                    for o in occ:
                        acc += W((p, q, o, s, t, o))
                    v2[p, q, s, t] = acc
    for p in range(nso):
        for s in range(nso):
            acc = 0.0
            for o in occ:
                for op in occ:
                    acc += W((p, o, op, s, o, op))
            v1[p, s] = 0.5 * acc
    v0 = 0.0
    for o in occ:
        for op in occ:
            for opp in occ:
                v0 += W((o, op, opp, o, op, opp))
    v0 /= 6.0
    # ---- convert normal-ordered (v0,v1,v2) to BARE coefficients (Wick, diagonal gamma)
    h1_bare = v1.copy()
    for p in range(nso):
        for s in range(nso):
            acc = 0.0
            for q in occ:
                acc += -v2[p, q, s, q] + v2[p, q, q, s]
            for o in occ:
                acc += v2[o, p, s, o] - v2[o, p, o, s]
            h1_bare[p, s] += 0.25 * acc
    v0_bare = v0
    for o in occ:
        v0_bare -= v1[o, o]
    for p in occ:
        for q in occ:
            v0_bare += 0.25 * (v2[p, q, p, q] - v2[p, q, q, p])
    return v2, h1_bare, v0_bare


# ----------------------------------------------------------------------------
# Ground states
# ----------------------------------------------------------------------------
def ground(H: np.ndarray, hermitian: bool = False) -> Tuple[float, float]:
    """Physical ground energy (and |Im| of the selected eigenvalue for non-Hermitian H)."""
    if hermitian:
        w = eigh(H, eigvals_only=True)
        return float(w[0]), 0.0
    ev = eig(H, right=False)
    real = ev.real
    imag = ev.imag
    mask = np.abs(imag) < 1e-6 * (1 + np.abs(real))
    cand = real[mask] if mask.any() else real
    E = float(np.min(cand))
    idx = np.argmin(np.where(mask, real, np.inf)) if mask.any() else np.argmin(real)
    return E, float(abs(imag[idx]))


def sym_ground(H: np.ndarray) -> float:
    """Hermitian ground energy of the symmetrized operator ``(H + H^dag)/2``."""
    Hs = 0.5 * (H + H.conj().T)
    return float(eigh(Hs, eigvals_only=True)[0])


# ============================================================================
# Operator-cost diagnostics (non-Hermitian friendly)
# ============================================================================
def herm_dev(A: np.ndarray) -> float:
    """Hermiticity deviation of an antisym 2-body tensor ``<pq||rs>``: max|A - A^dag|."""
    return float(np.max(np.abs(A - np.conjugate(np.transpose(A, (2, 3, 0, 1))))))


def herm_dev_matrix(M: np.ndarray) -> float:
    """Hermiticity deviation of a plain matrix: max|M - M^dag|."""
    return float(np.max(np.abs(M - M.conj().T)))


def dagger_asym(A: np.ndarray) -> np.ndarray:
    """Hermitian conjugate of an antisym 2-body tensor ``<pq||rs> -> <rs||pq>^*``."""
    return np.conjugate(np.transpose(A, (2, 3, 0, 1)))


def measure_matrix(H: np.ndarray, m_low: int = 4, deg_tol: float = 2e-3) -> Dict[str, object]:
    """Non-Hermitian FCI-matrix diagnostics: spectrum reality, non-normality,
    eigenvector condition numbers (full + low-lying block + per-eigenvalue
    Bauer-Fike), and the physical (degeneracy-merged) spectral gap."""
    H = np.asarray(H)
    n = H.shape[0]
    m_low = min(m_low, n)
    Hd = H.conj().T
    comm = Hd @ H - H @ Hd
    fH2 = np.linalg.norm(H, "fro") ** 2
    nonnorm = float(np.linalg.norm(comm, "fro") / fH2) if fH2 > 0 else 0.0

    w, VL, VR = _eig_lr(H)
    order = np.argsort(w.real)
    w = w[order]
    VR = VR[:, order]
    VL = VL[:, order]
    im_spread = float(np.max(np.abs(w.imag)))
    all_real = bool(im_spread < 1e-6 * (1 + np.max(np.abs(w.real))))
    E0 = float(w[0].real)
    im0 = float(abs(w[0].imag))

    re = w.real
    gap_raw = float(re[1] - re[0]) if n > 1 else 0.0
    gap_merged = 0.0
    for i in range(1, n):
        if re[i] - re[0] > deg_tol:
            gap_merged = float(re[i] - re[0])
            break
    in_mani = re[re - re[0] <= deg_tol]
    gs_split = float(in_mani.max() - in_mani.min())
    gs_degeneracy = int(in_mani.size)

    kV_full = float(np.linalg.cond(VR))
    kV_low = float(np.linalg.cond(VR[:, :m_low]))
    kappa_eig: List[float] = []
    for i in range(m_low):
        rr = VR[:, i]
        ll = VL[:, i]
        nr = np.linalg.norm(rr)
        nl = np.linalg.norm(ll)
        ov = abs(np.vdot(ll, rr))
        kappa_eig.append(float(nr * nl / ov) if ov > 0 else float("inf"))
    return dict(n=int(n), nonnormality=nonnorm, im_spread=im_spread,
                all_real=all_real, E0=E0, im0=im0,
                gap_raw=gap_raw, gap_merged=gap_merged,
                gs_split=gs_split, gs_degeneracy=gs_degeneracy,
                kV_full=kV_full, kV_low=kV_low,
                kappa_eig_low=kappa_eig, max_kappa_eig_low=float(max(kappa_eig)))


def _eig_lr(H: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    import scipy.linalg as sla
    return sla.eig(H, left=True, right=True)


def lcu_lambda(hso: np.ndarray, asym: np.ndarray, nso: int, tol: float = 1e-9
               ) -> Dict[str, float]:
    """Jordan-Wigner LCU 1-norm of ``sum h a+a + (1/4) sum <pq||rs> a+ a+ a a``.

    Complex coefficients are allowed (non-Hermitian friendly); the identity term is
    dropped.  Returns the Pauli-string count ``n_pauli``, the 1-norm ``lam``, and the
    largest imaginary Pauli coefficient ``max_imag_coeff``.
    """
    import openfermion as of

    ferm = of.FermionOperator()
    for p in range(nso):
        for q in range(nso):
            c = hso[p, q]
            if abs(c) > tol:
                ferm += of.FermionOperator(((p, 1), (q, 0)), complex(c))
    for p in range(nso):
        for q in range(nso):
            for r in range(nso):
                for s in range(nso):
                    c = 0.25 * asym[p, q, r, s]
                    if abs(c) > tol:
                        ferm += of.FermionOperator(((p, 1), (q, 1), (s, 0), (r, 0)),
                                                   complex(c))
    qop = of.jordan_wigner(ferm)
    qop.compress(1e-9)
    terms = {t: c for t, c in qop.terms.items() if t != ()}
    lam = float(sum(abs(c) for c in terms.values()))
    max_imag = float(max((abs(c.imag) for c in terms.values()), default=0.0))
    return dict(n_pauli=int(len(terms)), lam=lam, max_imag_coeff=max_imag)


# ============================================================================
# High-level assembly of an atomic TC/xTC system (s-only)
# ============================================================================
@dataclass
class AtomicTCSystem:
    """Container for the assembled s-only atomic TC/xTC operators (orthonormal MO)."""
    ns: int
    k: float
    gamma: float
    Z: int
    n_elec: int
    nso: int
    ndet: int
    dets: List[Det]
    didx: Dict[Det, int]
    hso: np.ndarray                 # one-body (spin-orbital)
    asym_coul: np.ndarray           # Coulomb <pq||rs> (Hermitian)
    asym_w: np.ndarray              # D term (Hermitian TC finite kernel)
    asym_K: np.ndarray              # convective K (NON-Hermitian)
    ref_occ: Tuple[int, ...]
    V3o: Optional[np.ndarray] = None    # orthonormal 3-body L3 (if built)
    v2: Optional[np.ndarray] = None     # xTC-contracted L3 -> 2-body (Hermitian)
    v1: Optional[np.ndarray] = None     # xTC 1-body correction
    v0: float = 0.0                     # xTC scalar

    # ---- operator variants: return (hso_eff, asym_eff, v0_eff) ----
    def plain(self) -> Tuple[np.ndarray, np.ndarray, float]:
        return self.hso, self.asym_coul, 0.0

    def d_only(self) -> Tuple[np.ndarray, np.ndarray, float]:
        return self.hso, self.asym_w, 0.0

    def tc2(self) -> Tuple[np.ndarray, np.ndarray, float]:
        """2-body TC ``D + K`` (non-Hermitian; no 3-body)."""
        return self.hso, self.asym_w + self.asym_K, 0.0

    def xtc_noK(self) -> Tuple[np.ndarray, np.ndarray, float]:
        """Hermitian ``D + xTC-L3`` (drops the convective K)."""
        assert self.v2 is not None, "xTC-L3 not built (with_L3=False)"
        return self.hso + self.v1, self.asym_w + self.v2, self.v0

    def xtc_full(self) -> Tuple[np.ndarray, np.ndarray, float]:
        """Full non-Hermitian ``D + K + xTC-L3``."""
        assert self.v2 is not None, "xTC-L3 not built (with_L3=False)"
        return self.hso + self.v1, self.asym_w + self.asym_K + self.v2, self.v0


def build_atomic_system(ns: int, k: float, gamma: float, Z: int, n_elec: int,
                        Ng: int = 800, nx: int = 128, with_L3: bool = True
                        ) -> AtomicTCSystem:
    """Assemble the s-only atomic TC/xTC operators in the orthonormal MO basis.

    Mirrors ``debug/xtc_quantum_cost_measure.build_sonly`` exactly.  The aufbau
    reference occupation is chosen from the one-body diagonal energies (lowest 1 or 2
    spatial orbitals; for ``n_elec == 3`` the doublet ``1s^2 2s^1``).
    """
    r, wr = make_grid(k, Ng=Ng)
    S, h1s, Rtab, W = build_one_body(ns, r, wr, k, Z)
    Km = build_kernels(r, gamma, nx=nx)
    eri_coul, eri_w, eri_K = two_body(ns, Rtab, W, Km)
    V3 = three_body(ns, Rtab, W, Km) if with_L3 else None
    X = lowdin(S)
    h1o = transform_1(h1s, X)
    eri_coul_o = transform_2(eri_coul, X)
    eri_w_o = transform_2(eri_w, X)
    eri_K_o = transform_2(eri_K, X)
    V3o = transform_3(V3, X) if with_L3 else None
    nso = 2 * ns
    dets, didx = make_dets(nso, n_elec)
    hso = h_spin(h1o, nso)
    asym_coul = asym_from_phys(eri_coul_o, nso)
    asym_w = asym_from_phys(eri_w_o, nso)
    asym_K = asym_from_phys(eri_K_o, nso)
    diag_e = np.diag(h1o)
    order = np.argsort(diag_e)
    if n_elec == 2:
        o0 = order[0]
        ref_occ: Tuple[int, ...] = (2 * o0, 2 * o0 + 1)
    else:
        o0, o1 = order[0], order[1]
        ref_occ = (2 * o0, 2 * o0 + 1, 2 * o1)
    v2 = v1 = None
    v0 = 0.0
    if with_L3:
        v2, v1, v0 = xtc_contract(V3o, nso, ref_occ)
    return AtomicTCSystem(ns=ns, k=k, gamma=gamma, Z=Z, n_elec=n_elec, nso=nso,
                          ndet=len(dets), dets=dets, didx=didx, hso=hso,
                          asym_coul=asym_coul, asym_w=asym_w, asym_K=asym_K,
                          ref_occ=ref_occ, V3o=V3o, v2=v2, v1=v1, v0=v0)


def build_fci_matrix(sys: AtomicTCSystem, variant: Tuple[np.ndarray, np.ndarray, float]
                     ) -> np.ndarray:
    """FCI matrix for a variant ``(hso_eff, asym_eff, v0_eff)`` of the system."""
    hso, asym, v0 = variant
    return build_H(sys.dets, sys.didx, hso, asym, sys.nso, v0=v0)


def plain_energy(sys: AtomicTCSystem) -> float:
    """Plain (Coulomb) Hermitian FCI ground energy."""
    H = build_fci_matrix(sys, sys.plain())
    E, _ = ground(H, hermitian=True)
    return E


def tc_energies(sys: AtomicTCSystem) -> Tuple[float, float, float]:
    """Full non-Hermitian TC ground ``E_TC`` (+ |Im|) and symmetrized ground ``E_sym``.

    For ``with_L3`` systems the full operator is ``D + K + xTC-L3``; otherwise the
    2-body ``D + K``.  ``E_sym`` is the Hermitian ground of ``(H + H^dag)/2``.
    """
    variant = sys.xtc_full() if sys.v2 is not None else sys.tc2()
    H = build_fci_matrix(sys, variant)
    E_TC, im = ground(H, hermitian=False)
    E_sym = sym_ground(H)
    return E_TC, im, E_sym


def plain_fci(ns: int, k: float, Z: int = 3, n_elec: int = 3, Ng: int = 1000,
              nx: int = 128) -> Tuple[float, Dict[str, float], int, int]:
    """Fast plain (Coulomb-only) Hermitian FCI + effective-operator LCU metrics."""
    r, wr = make_grid(k, Ng=Ng)
    S, h1s, Rtab, W = build_one_body(ns, r, wr, k, Z)
    Kc = build_coul_kernel(r, nx=nx)
    RR = {i: Rtab[i + 1][0] for i in range(ns)}
    D = {(i, kk): RR[i] * RR[kk] * W for i in range(ns) for kk in range(ns)}
    eri = np.zeros((ns, ns, ns, ns))
    for i in range(ns):
        for j in range(ns):
            for kk in range(ns):
                for l in range(ns):
                    eri[i, j, kk, l] = D[(i, kk)] @ Kc @ D[(j, l)]
    X = lowdin(S)
    h1o = transform_1(h1s, X)
    eri_o = transform_2(eri, X)
    nso = 2 * ns
    dets, didx = make_dets(nso, n_elec)
    hso = h_spin(h1o, nso)
    asym = asym_from_phys(eri_o, nso)
    H = build_H(dets, didx, hso, asym, nso)
    E, _ = ground(H, hermitian=True)
    lam = lcu_lambda(hso, asym, nso)
    return E, lam, nso, len(dets)
