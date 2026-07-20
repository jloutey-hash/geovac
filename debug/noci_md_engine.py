"""N4 engine -- s/p Gaussian integrals + general-N non-orthogonal Slater-Condon.

NOCI sandbox (branch sandbox/noci, repo frozen -- no release path).  Step N4 of
debug/noci_sandbox_notes.md needs two things N1/N2 lacked:

  1. p-type integrals (Na 2p core): a McMurchie-Davidson engine for contracted
     Cartesian Gaussians of arbitrary angular momentum (validated s+p here).
  2. N-electron non-orthogonal matrix elements at N=12 (NaH all-electron):
     N2's explicit permutation expansion dies past ~6 electrons; here we use the
     Loewdin cofactor rules (first- and second-order complementary minors).

Plus a bitstring FCI (orthonormal Loewdin basis) for the in-basis reference,
and a least-squares STO->Gaussian shape fitter for the 2p / 3s shapes (1s and
2s shapes reuse the N1/N2 hardcoded STO-6G/STO-3G expansions).

Running this module directly executes the validation suite:
  V-A  s-only agreement vs N1's closed-form s engine (exact same physics)
  V-B  p-primitive integrals vs center-derivatives of N1 s integrals (FD)
  V-C  fitted-shape atom energies (H 1s, H 2p at zeta=1/2, vs exact -0.5/-0.125)
  V-D  Loewdin-cofactor machinery vs N2's permutation machinery, random tensors
  V-E  Loewdin-cofactor machinery vs N2's stored LiH ladder (real system)
  V-F  bitstring FCI vs gensc(s=I) on random tensors + vs N2's stored LiH FCI
"""

from __future__ import annotations

import json
import math
import os
from itertools import combinations
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.special import hyp1f1

HARTREE_TO_EV = 27.211386

# --------------------------------------------------------------- primitives


def _dfact(n: int) -> int:
    """(2n-1)!! with the n=0 convention (-1)!! = 1."""
    out = 1
    k = 2 * n - 1
    while k > 1:
        out *= k
        k -= 2
    return out


def prim_norm(alphas: np.ndarray, lmn: Tuple[int, int, int]) -> np.ndarray:
    i, j, k = lmn
    L = i + j + k
    return ((2.0 * alphas / np.pi) ** 0.75 * (4.0 * alphas) ** (L / 2.0)
            / math.sqrt(_dfact(i) * _dfact(j) * _dfact(k)))


class BasisFn:
    """Contracted Cartesian Gaussian: sum_i d_i N_i x^l y^m z^n exp(-a_i r^2),
    d_i given on NORMALIZED primitives; the contraction is renormalized."""

    def __init__(self, center, lmn: Tuple[int, int, int],
                 alphas: np.ndarray, dcoeffs: np.ndarray):
        self.center = np.asarray(center, dtype=float)
        self.lmn = tuple(int(x) for x in lmn)
        self.alphas = np.asarray(alphas, dtype=float)
        self.coeffs = np.asarray(dcoeffs, dtype=float) * prim_norm(self.alphas, self.lmn)
        self.coeffs = self.coeffs / math.sqrt(overlap_md(self, self))


def _E(i: int, j: int, t: int, Q: float, a, b, cache) -> np.ndarray:
    """Hermite expansion coefficient E_t^{ij} (arrays over primitive pairs).
    Q = A_x - B_x;  X_PA = -(b/p) Q,  X_PB = +(a/p) Q."""
    if t < 0 or t > i + j or i < 0 or j < 0:
        return 0.0
    key = (i, j, t)
    if key in cache:
        return cache[key]
    p = a + b
    if i == 0 and j == 0 and t == 0:
        val = np.exp(-(a * b / p) * Q * Q)
    elif j == 0:
        val = (_E(i - 1, j, t - 1, Q, a, b, cache) / (2.0 * p)
               - (b * Q / p) * _E(i - 1, j, t, Q, a, b, cache)
               + (t + 1) * _E(i - 1, j, t + 1, Q, a, b, cache))
    else:
        val = (_E(i, j - 1, t - 1, Q, a, b, cache) / (2.0 * p)
               + (a * Q / p) * _E(i, j - 1, t, Q, a, b, cache)
               + (t + 1) * _E(i, j - 1, t + 1, Q, a, b, cache))
    cache[key] = val
    return val


def boys(n: int, x: np.ndarray) -> np.ndarray:
    return hyp1f1(n + 0.5, n + 1.5, -x) / (2.0 * n + 1.0)


def _R_hermite(tmax: int, umax: int, vmax: int, alpha, X, Y, Z) -> Dict:
    """Hermite Coulomb integrals R^0_{tuv} as arrays; standard recursion."""
    r2 = X * X + Y * Y + Z * Z
    nmax = tmax + umax + vmax
    Fn = [boys(n, alpha * r2) for n in range(nmax + 1)]
    cache: Dict[Tuple[int, int, int, int], np.ndarray] = {}

    def R(n, t, u, v):
        if t < 0 or u < 0 or v < 0:
            return 0.0
        key = (n, t, u, v)
        if key in cache:
            return cache[key]
        if t == u == v == 0:
            val = (-2.0 * alpha) ** n * Fn[n]
        elif t > 0:
            val = (t - 1) * R(n + 1, t - 2, u, v) + X * R(n + 1, t - 1, u, v)
        elif u > 0:
            val = (u - 1) * R(n + 1, t, u - 2, v) + Y * R(n + 1, t, u - 1, v)
        else:
            val = (v - 1) * R(n + 1, t, u, v - 2) + Z * R(n + 1, t, u, v - 1)
        cache[key] = val
        return val

    return {(t, u, v): R(0, t, u, v)
            for t in range(tmax + 1) for u in range(umax + 1)
            for v in range(vmax + 1)}


class _Pair:
    """Precomputed primitive-pair data for a contracted orbital pair."""

    def __init__(self, a: BasisFn, b: BasisFn):
        A, B = a.center, b.center
        al, bl = np.meshgrid(a.alphas, b.alphas, indexing="ij")
        self.al = al.ravel()
        self.bl = bl.ravel()
        self.cc = np.outer(a.coeffs, b.coeffs).ravel()
        self.p = self.al + self.bl
        self.P = (self.al[:, None] * A + self.bl[:, None] * B) / self.p[:, None]
        self.Q = A - B
        self.lmn_a = a.lmn
        self.lmn_b = b.lmn
        self.caches = [dict(), dict(), dict()]

    def E_dim(self, d: int, i: int, j: int, t: int):
        return _E(i, j, t, self.Q[d], self.al, self.bl, self.caches[d]
                  if (i, j) == (self.lmn_a[d], self.lmn_b[d]) else
                  self._scratch(d, i, j))

    def _scratch(self, d, i, j):
        # separate cache for shifted (kinetic) calls to avoid key collisions
        key = ("scr", d, i, j)
        if not hasattr(self, "_scr"):
            self._scr = {}
        if key not in self._scr:
            self._scr[key] = dict()
        return self._scr[key]

    def hermite_components(self):
        """List of ((t,u,v), E_t*E_u*E_v arrays) for the pair's own lmn."""
        ia, ja = self.lmn_a, self.lmn_b
        comps = []
        for t in range(ia[0] + ja[0] + 1):
            Ex = self.E_dim(0, ia[0], ja[0], t)
            for u in range(ia[1] + ja[1] + 1):
                Ey = self.E_dim(1, ia[1], ja[1], u)
                for v in range(ia[2] + ja[2] + 1):
                    Ez = self.E_dim(2, ia[2], ja[2], v)
                    comps.append(((t, u, v), Ex * Ey * Ez))
        return comps


def _s_prim(pair: _Pair, lmn_b: Tuple[int, int, int]) -> np.ndarray:
    """Primitive overlap arrays with the KET angular momentum overridden."""
    if min(lmn_b) < 0:
        return np.zeros_like(pair.p)
    ia = pair.lmn_a
    val = (np.pi / pair.p) ** 1.5
    for d in range(3):
        val = val * _E(ia[d], lmn_b[d], 0, pair.Q[d], pair.al, pair.bl,
                       pair._scratch(d, ia[d], lmn_b[d]))
    return val


def overlap_md(a: BasisFn, b: BasisFn) -> float:
    pair = _Pair(a, b)
    return float(np.sum(pair.cc * _s_prim(pair, b.lmn)))


def kinetic_md(a: BasisFn, b: BasisFn) -> float:
    pair = _Pair(a, b)
    j1, j2, j3 = b.lmn
    bl = pair.bl
    term = bl * (2.0 * (j1 + j2 + j3) + 3.0) * _s_prim(pair, b.lmn)
    term -= 2.0 * bl ** 2 * (_s_prim(pair, (j1 + 2, j2, j3))
                             + _s_prim(pair, (j1, j2 + 2, j3))
                             + _s_prim(pair, (j1, j2, j3 + 2)))
    term -= 0.5 * (j1 * (j1 - 1) * _s_prim(pair, (j1 - 2, j2, j3))
                   + j2 * (j2 - 1) * _s_prim(pair, (j1, j2 - 2, j3))
                   + j3 * (j3 - 1) * _s_prim(pair, (j1, j2, j3 - 2)))
    return float(np.sum(pair.cc * term))


def nuclear_md(a: BasisFn, b: BasisFn, C, Z: float) -> float:
    pair = _Pair(a, b)
    C = np.asarray(C, dtype=float)
    comps = pair.hermite_components()
    tmax = max(k[0][0] for k in comps)
    umax = max(k[0][1] for k in comps)
    vmax = max(k[0][2] for k in comps)
    PC = pair.P - C
    Rtab = _R_hermite(tmax, umax, vmax, pair.p, PC[:, 0], PC[:, 1], PC[:, 2])
    total = np.zeros_like(pair.p)
    for (t, u, v), Etuv in comps:
        total += Etuv * Rtab[(t, u, v)]
    return float(-Z * np.sum(pair.cc * (2.0 * np.pi / pair.p) * total))


def eri_md(a: BasisFn, b: BasisFn, c: BasisFn, d: BasisFn) -> float:
    """Chemist (ab|cd) = int a(1)b(1) 1/r12 c(2)d(2)."""
    p1 = _Pair(a, b)
    p2 = _Pair(c, d)
    c1 = p1.hermite_components()
    c2 = p2.hermite_components()
    t1 = max(k[0][0] for k in c1) + max(k[0][0] for k in c2)
    u1 = max(k[0][1] for k in c1) + max(k[0][1] for k in c2)
    v1 = max(k[0][2] for k in c1) + max(k[0][2] for k in c2)
    pp = p1.p[:, None]
    qq = p2.p[None, :]
    alpha = pp * qq / (pp + qq)
    PQ = p1.P[:, None, :] - p2.P[None, :, :]
    Rtab = _R_hermite(t1, u1, v1, alpha, PQ[:, :, 0], PQ[:, :, 1], PQ[:, :, 2])
    pref = 2.0 * np.pi ** 2.5 / (pp * qq * np.sqrt(pp + qq))
    total = np.zeros_like(alpha)
    for (t, u, v), E1 in c1:
        for (tt, uu, vv), E2 in c2:
            sign = (-1) ** (tt + uu + vv)
            total += sign * E1[:, None] * E2[None, :] * Rtab[(t + tt, u + uu, v + vv)]
    cc = p1.cc[:, None] * p2.cc[None, :]
    return float(np.sum(cc * pref * total))


def integral_set_md(orbs: List[BasisFn], nuclei: List[Tuple[np.ndarray, float]]):
    """(s, h, g) over a contracted basis; 8-fold ERI symmetry exploited."""
    n = len(orbs)
    s = np.zeros((n, n))
    h = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            s[i, j] = s[j, i] = overlap_md(orbs[i], orbs[j])
            hij = kinetic_md(orbs[i], orbs[j]) + sum(
                nuclear_md(orbs[i], orbs[j], pos, z) for pos, z in nuclei)
            h[i, j] = h[j, i] = hij
    g = np.zeros((n, n, n, n))
    pair_idx = [(i, j) for i in range(n) for j in range(i, n)]
    for pi_, (i, j) in enumerate(pair_idx):
        for k, l in pair_idx[pi_:]:
            val = eri_md(orbs[i], orbs[j], orbs[k], orbs[l])
            for (x, y) in ((i, j), (j, i)):
                for (z, w) in ((k, l), (l, k)):
                    g[x, y, z, w] = val
                    g[z, w, x, y] = val
    return s, h, g


# ------------------------------------------- general-N non-orthogonal dets


def det_pair_gensc(so_i: Sequence[Tuple[int, int]],
                   so_j: Sequence[Tuple[int, int]],
                   s: np.ndarray, h: np.ndarray, g: np.ndarray
                   ) -> Tuple[float, float]:
    """<D_I|D_J> and <D_I|H|D_J> via Loewdin cofactor rules (any overlap
    pattern, any N).  Spin-orbitals are (spatial, spin); chemist g."""
    N = len(so_i)
    M = np.zeros((N, N))
    for a_, (p, sp) in enumerate(so_i):
        for b_, (q, sq) in enumerate(so_j):
            M[a_, b_] = s[p, q] if sp == sq else 0.0
    S = float(np.linalg.det(M))

    # first-order cofactors
    H1 = 0.0
    rows = np.arange(N)
    for a_ in range(N):
        pa, sa = so_i[a_]
        sub_rows = np.delete(rows, a_)
        for b_ in range(N):
            qb, sb = so_j[b_]
            if sa != sb:
                continue
            sub = M[np.ix_(sub_rows, np.delete(rows, b_))]
            D1 = (-1) ** (a_ + b_) * float(np.linalg.det(sub)) if N > 1 else 1.0
            H1 += h[pa, qb] * D1

    # second-order cofactors
    H2 = 0.0
    if N >= 2:
        for a_ in range(N):
            pa, sa = so_i[a_]
            for c_ in range(a_ + 1, N):
                pc, sc = so_i[c_]
                rkeep = np.delete(rows, [a_, c_])
                for b_ in range(N):
                    qb, sb = so_j[b_]
                    for d_ in range(b_ + 1, N):
                        qd, sd = so_j[d_]
                        w = 0.0
                        if sa == sb and sc == sd:
                            w += g[pa, qb, pc, qd]
                        if sa == sd and sc == sb:
                            w -= g[pa, qd, pc, qb]
                        if w == 0.0:
                            continue
                        ckeep = np.delete(rows, [b_, d_])
                        sub = M[np.ix_(rkeep, ckeep)]
                        D2 = ((-1) ** (a_ + b_ + c_ + d_)
                              * (float(np.linalg.det(sub)) if N > 2 else 1.0))
                        H2 += w * D2
    return S, H1 + H2


def noci_ground_gensc(dets, s, h, g, tol: float = 1e-10):
    """Ground energy over a determinant list (generalized eigenproblem)."""
    nd = len(dets)
    smat = np.zeros((nd, nd))
    hmat = np.zeros((nd, nd))
    for i in range(nd):
        for j in range(i, nd):
            sij, hij = det_pair_gensc(dets[i], dets[j], s, h, g)
            smat[i, j] = smat[j, i] = sij
            hmat[i, j] = hmat[j, i] = hij
    w, v = np.linalg.eigh(smat)
    keep = w > tol * w.max()
    x = v[:, keep] / np.sqrt(w[keep])
    e = float(np.linalg.eigvalsh(x.T @ hmat @ x)[0])
    cond = float(w.max() / w.min()) if w.min() > 0 else np.inf
    return e, cond


# --------------------------------------------------- bitstring FCI (ortho)


def lowdin_orbitals(s: np.ndarray):
    w, v = np.linalg.eigh(s)
    return v @ np.diag(w ** -0.5) @ v.T


def transform_integrals(x: np.ndarray, h: np.ndarray, g: np.ndarray):
    ht = x.T @ h @ x
    gt = np.einsum("pi,qj,rk,sl,pqrs->ijkl", x, x, x, x, g, optimize=True)
    return ht, gt


def fci_ground(h: np.ndarray, g: np.ndarray, n_elec: int,
               n_states: int = 1) -> float:
    """FCI over an ORTHONORMAL spatial basis (h, g chemist), spin-orbitals
    so = (spatial, spin) enumerated as 2*p + spin.  Slater-Condon rules,
    dense H (fine for dim <= a few thousand)."""
    n_sp = h.shape[0]
    sos = [(p, sp) for p in range(n_sp) for sp in (0, 1)]
    dets = [tuple(sorted(c)) for c in combinations(range(2 * n_sp), n_elec)]
    idx = {d: i for i, d in enumerate(dets)}
    dim = len(dets)

    def spat(so):
        return sos[so][0]

    def spin(so):
        return sos[so][1]

    def anti(pq, rs):
        """<pq||rs> physicist antisymmetrized over spin-orbitals."""
        p, q = pq
        r, s_ = rs
        val = 0.0
        if spin(p) == spin(r) and spin(q) == spin(s_):
            val += g[spat(p), spat(r), spat(q), spat(s_)]
        if spin(p) == spin(s_) and spin(q) == spin(r):
            val -= g[spat(p), spat(s_), spat(q), spat(r)]
        return val

    H = np.zeros((dim, dim))
    for I, di in enumerate(dets):
        occ = set(di)
        # diagonal
        e = sum(h[spat(p), spat(p)] for p in di if True)
        e = 0.0
        for p in di:
            if spin(p) == spin(p):
                e += h[spat(p), spat(p)]
        for ai in range(len(di)):
            for bi in range(ai + 1, len(di)):
                e += anti((di[ai], di[bi]), (di[ai], di[bi]))
        H[I, I] = e
        # singles and doubles
        for J in range(I + 1, dim):
            dj = dets[J]
            occj = set(dj)
            diff_i = sorted(occ - occj)
            diff_j = sorted(occj - occ)
            nd = len(diff_i)
            if nd > 2:
                continue
            if nd == 1:
                p, q = diff_i[0], diff_j[0]
                if spin(p) != spin(q):
                    continue
                phase = (-1) ** (di.index(p) + dj.index(q))
                val = h[spat(p), spat(q)]
                for r in di:
                    if r == p:
                        continue
                    val += anti((p, r), (q, r))
                H[I, J] = H[J, I] = phase * val
            else:
                p1, p2 = diff_i
                q1, q2 = diff_j
                phase = (-1) ** (di.index(p1) + di.index(p2)
                                 + dj.index(q1) + dj.index(q2))
                H[I, J] = H[J, I] = phase * anti((p1, p2), (q1, q2))
    evals = np.linalg.eigvalsh(H)
    return float(evals[0]) if n_states == 1 else evals[:n_states]


# ------------------------------------------------------- STO shape fitting


def fit_sto_shape(l: int, n_r: int, n_gauss: int = 6,
                  r_max: float = 40.0, n_grid: int = 6000):
    """Least-squares fit of the zeta=1 Slater radial shape r^{n_r-1} e^{-r}
    by n_gauss Gaussians r^l exp(-a_i r^2) (radial part of an l-primitive).
    Returns (alphas, dcoeffs_on_normalized_primitives, fit_overlap)."""
    from scipy.optimize import minimize

    r = np.linspace(1e-6, r_max, n_grid)
    w = r * r
    target = r ** (n_r - 1) * np.exp(-r)
    tnorm = math.sqrt(np.trapezoid(target * target * w, r))
    target = target / tnorm

    def basis_mat(log_a):
        a = np.exp(log_a)
        return (r[:, None] ** l) * np.exp(-a[None, :] * r[:, None] ** 2), a

    def resid(log_a):
        B, a = basis_mat(log_a)
        G = (B * w[:, None]).T @ B * (r[1] - r[0])
        t = (B * w[:, None]).T @ target * (r[1] - r[0])
        try:
            c = np.linalg.solve(G + 1e-14 * np.eye(len(a)), t)
        except np.linalg.LinAlgError:
            return 1e6
        fit = B @ c
        return float(np.trapezoid((fit - target) ** 2 * w, r))

    x0 = np.log(np.geomspace(0.03, 60.0, n_gauss))
    best = minimize(resid, x0, method="Nelder-Mead",
                    options={"maxiter": 4000, "xatol": 1e-8, "fatol": 1e-14})
    B, a = basis_mat(best.x)
    G = (B * w[:, None]).T @ B * (r[1] - r[0])
    t = (B * w[:, None]).T @ target * (r[1] - r[0])
    c = np.linalg.solve(G + 1e-14 * np.eye(len(a)), t)
    fit = B @ c
    fnorm = math.sqrt(np.trapezoid(fit * fit * w, r))
    quality = float(np.trapezoid(fit * target * w, r) / fnorm)

    # convert raw radial coefficients -> coefficients on normalized 3D primitives.
    # 3D primitive (z-type example) has radial part r^l e^{-a r^2} times the
    # angular factor; the radial-only normalization constant of the primitive is
    # Nrad = sqrt( 2 (2a)^{l+3/2} ... ) -- easiest numerically:
    dco = np.empty(len(a))
    for i_ in range(len(a)):
        prim = r ** l * np.exp(-a[i_] * r ** 2)
        pn = math.sqrt(np.trapezoid(prim * prim * w, r))
        dco[i_] = c[i_] * pn
    dco = dco / np.linalg.norm(dco) if False else dco
    return a, dco, quality


def sto_shape_basis(center, kind: str, zeta: float,
                    shapes: Dict[str, Tuple[np.ndarray, np.ndarray]],
                    lmn: Tuple[int, int, int]) -> BasisFn:
    """Build a BasisFn from a fitted zeta=1 shape, scaled to zeta."""
    a, dco = shapes[kind]
    return BasisFn(center, lmn, a * zeta ** 2, dco)


# hardcoded zeta=1 shapes from N1/N2 (s-type; coefficients on normalized prims)
STO6G_1S = (
    np.array([23.31030, 4.235916, 1.185057, 0.4070989, 0.1580884, 0.06510954]),
    np.array([0.00916360, 0.04936150, 0.16853830, 0.37056280, 0.41649150, 0.13033400]),
)
STO3G_2S = (
    np.array([0.9942030, 0.2310310, 0.0751386]),
    np.array([-0.09996723, 0.39951283, 0.70011547]),
)


# ------------------------------------------------------------- validations


def _random_sym_tensors(m: int, rng) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    s = rng.standard_normal((m, m)) * 0.2
    s = 0.5 * (s + s.T) + np.eye(m) * (1.0 + 0.5)
    h = rng.standard_normal((m, m))
    h = 0.5 * (h + h.T)
    g = rng.standard_normal((m, m, m, m))
    g = (g + g.transpose(1, 0, 2, 3) + g.transpose(0, 1, 3, 2)
         + g.transpose(1, 0, 3, 2) + g.transpose(2, 3, 0, 1)
         + g.transpose(3, 2, 0, 1) + g.transpose(2, 3, 1, 0)
         + g.transpose(3, 2, 1, 0)) / 8.0
    return s, h, g


def run_validations() -> None:
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from noci_h2_probe import (Basis1s, STO6G_1S as N1_6G, eri as eri_s,
                               integral_set, kinetic as kin_s,
                               nuclear as nuc_s, overlap as ov_s)
    from noci_lih_probe import STO3G_2S as N2_2S, det_pair_elements

    print("=== N4 engine validation suite ===\n")
    rng = np.random.default_rng(7)

    # V-A: s-only agreement vs N1 closed forms (LiH configuration, R=3.25)
    pos_li = np.array([0.0, 0.0, 0.0])
    pos_h = np.array([0.0, 0.0, 3.25])
    n1_orbs = [Basis1s(pos_li, 2.69, *N1_6G), Basis1s(pos_li, 0.65, *N2_2S),
               Basis1s(pos_h, 1.0, *N1_6G)]
    md_orbs = [BasisFn(pos_li, (0, 0, 0), N1_6G[0] * 2.69 ** 2, N1_6G[1]),
               BasisFn(pos_li, (0, 0, 0), N2_2S[0] * 0.65 ** 2, N2_2S[1]),
               BasisFn(pos_h, (0, 0, 0), N1_6G[0], N1_6G[1])]
    nuclei = [(pos_li, 3.0), (pos_h, 1.0)]
    s1, h1, g1 = integral_set(n1_orbs, nuclei)
    s2, h2, g2 = integral_set_md(md_orbs, nuclei)
    print(f"V-A  s-only vs N1 closed forms: max|dS| = {np.max(np.abs(s1 - s2)):.2e}"
          f"  max|dh| = {np.max(np.abs(h1 - h2)):.2e}"
          f"  max|dg| = {np.max(np.abs(g1 - g2)):.2e}")

    # V-B: p-primitive vs finite-difference center-derivative of s integrals
    al = np.array([0.8])
    one = np.array([1.0])
    A = np.array([0.1, -0.2, 0.3])
    B = np.array([0.5, 0.4, -1.1])
    eps = 1e-5

    def s_prim_at(Az):
        aa = Basis1s(np.array([A[0], A[1], Az]), 1.0, al, one)
        bb = Basis1s(B, 1.0, np.array([1.3]), one)
        return aa, bb

    # d/dA_z of the N1 s-orbital integrals; p_z = deriv / (2 alpha) on the raw
    # (un-renormalized) primitive -- compare RATIOS via normalized functions:
    a_p = BasisFn(A, (0, 0, 1), al, one)
    b_s = BasisFn(B, (0, 0, 0), np.array([1.3]), one)
    ap_norm = a_p.coeffs[0]
    ovs = []
    for f, fmd in ((ov_s, overlap_md), (kin_s, kinetic_md)):
        ap, bp = s_prim_at(A[2] + eps)
        am, bm = s_prim_at(A[2] - eps)
        raw_p = ap.coeffs[0]  # normalized s coeff; deriv of normalized s prim
        der = (f(ap, bp) - f(am, bm)) / (2 * eps)
        # normalized s primitive: N_s * e^{-a r^2}; d/dA_z -> N_s * 2a (z-Az) e
        # = (N_s / N_p) * 2a/(2a)... : p_z normalized = (1/sqrt(2a... use ratio:
        # d/dA_z [N_s g_s] = N_s * 2a * (z-A_z) g_s ; normalized p prim is
        # N_p (z-A_z) g_s  =>  deriv = (N_s * 2a / N_p) * [normalized p integral]
        Ns = ap.coeffs[0]
        Np = ap_norm
        pred = der * Np / (Ns * 2.0 * al[0])
        got = fmd(a_p, b_s)
        ovs.append(abs(pred - got))
    print(f"V-B  p_z primitive vs FD center-derivative: overlap diff = {ovs[0]:.2e}"
          f"  kinetic diff = {ovs[1]:.2e}")
    # nuclear + eri single checks
    Cn = np.array([-0.3, 0.7, 0.2])
    ap, bp = s_prim_at(A[2] + eps)
    am, bm = s_prim_at(A[2] - eps)
    der = (nuc_s(ap, bp, Cn, 1.0) - nuc_s(am, bm, Cn, 1.0)) / (2 * eps)
    Ns = ap.coeffs[0]
    pred = der * ap_norm / (Ns * 2.0 * al[0])
    got = nuclear_md(a_p, b_s, Cn, 1.0)
    print(f"     nuclear diff = {abs(pred - got):.2e}", end="")
    cc = Basis1s(Cn, 1.2, np.array([0.9]), one)
    dd = Basis1s(np.array([0.0, 0.0, 0.0]), 0.7, np.array([0.5]), one)
    der = (eri_s(ap, bp, cc, dd) - eri_s(am, bm, cc, dd)) / (2 * eps)
    pred = der * ap_norm / (Ns * 2.0 * al[0])
    cmd = BasisFn(Cn, (0, 0, 0), np.array([0.9 * 1.2 ** 2]), one)
    dmd = BasisFn(np.array([0.0, 0.0, 0.0]), (0, 0, 0), np.array([0.5 * 0.7 ** 2]), one)
    got = eri_md(a_p, b_s, cmd, dmd)
    print(f"  eri diff = {abs(pred - got):.2e}")

    # V-C: fitted shapes -> atom energies
    shapes = {}
    for kind, (l, n_r) in (("1s", (0, 1)), ("2p", (1, 2)), ("3s", (0, 3))):
        a, dco, q = fit_sto_shape(l, n_r)
        shapes[kind] = (a, dco)
        print(f"V-C  fit {kind}: quality <fit|STO> = {q:.6f}")
    orig = np.array([0.0, 0.0, 0.0])
    h_nuc = [(orig, 1.0)]
    chi = sto_shape_basis(orig, "1s", 1.0, shapes, (0, 0, 0))
    sA, hA, _ = integral_set_md([chi], h_nuc)
    print(f"     E(H 1s, fitted, zeta=1)   = {hA[0, 0] / sA[0, 0]:.6f}  (exact -0.5)")
    chi = sto_shape_basis(orig, "2p", 0.5, shapes, (0, 0, 1))
    sA, hA, _ = integral_set_md([chi], h_nuc)
    print(f"     E(H 2p_z, fitted, zeta=.5) = {hA[0, 0] / sA[0, 0]:.6f}  (exact -0.125)")

    # V-D: gensc vs permutation machinery on random tensors
    worst = 0.0
    for m, nelec, trials in ((3, 2, 6), (4, 3, 6), (4, 4, 4)):
        s, h, g = _random_sym_tensors(m, rng)
        sos = [(p, sp) for p in range(m) for sp in (0, 1)]
        for _ in range(trials):
            di = list(rng.choice(len(sos), size=nelec, replace=False))
            dj = list(rng.choice(len(sos), size=nelec, replace=False))
            so_i = [sos[k] for k in di]
            so_j = [sos[k] for k in dj]
            s_ref, h_ref = det_pair_elements(so_i, so_j, s, h, g)
            s_new, h_new = det_pair_gensc(so_i, so_j, s, h, g)
            worst = max(worst, abs(s_ref - s_new), abs(h_ref - h_new))
    print(f"V-D  gensc vs permutation machinery (random tensors): worst = {worst:.2e}")

    # V-E: gensc on the real LiH ladder vs stored N2 results at R=3.25
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "data", "noci_lih_probe_results.json")) as fh:
        n2 = json.load(fh)
    z1s, z2s, zh = n2["zetas"]["Li1s"], n2["zetas"]["Li2s"], n2["zetas"]["H"]
    n1_orbs = [Basis1s(pos_li, z1s, *N1_6G), Basis1s(pos_li, z2s, *N2_2S),
               Basis1s(pos_h, zh, *N1_6G)]
    s, h, g = integral_set(n1_orbs, nuclei)
    vnn = 3.0 / 3.25
    UP, DN = 0, 1
    core = [(0, UP), (0, DN)]
    det_cov_a = core + [(1, UP), (2, DN)]
    det_cov_b = core + [(2, UP), (1, DN)]
    det_ion_h = core + [(2, UP), (2, DN)]
    row = next(r for r in n2["rows"] if abs(r["R"] - 3.25) < 1e-9)
    e2, _ = noci_ground_gensc([det_cov_a, det_cov_b], s, h, g)
    e3, _ = noci_ground_gensc([det_cov_a, det_cov_b, det_ion_h], s, h, g)
    print(f"V-E  gensc real LiH R=3.25: |dE(cov2)| = {abs(e2 + vnn - row['cov (2 dets)']):.2e}"
          f"  |dE(cov+ionH)| = {abs(e3 + vnn - row['cov+ionH (3 dets)']):.2e}")

    # V-F: bitstring FCI vs stored N2 FCI + vs gensc(s=I) on random tensors
    x = lowdin_orbitals(s)
    ht, gt = transform_integrals(x, h, g)
    e_fci = fci_ground(ht, gt, 4) + vnn
    print(f"V-F  bitstring FCI vs stored N2 FCI at R=3.25: |dE| = "
          f"{abs(e_fci - row['FCI (15 dets)']):.2e}")
    m = 3
    s_, h_, g_ = _random_sym_tensors(m, rng)
    x_ = lowdin_orbitals(s_)
    ht_, gt_ = transform_integrals(x_, h_, g_)
    e_bit = fci_ground(ht_, gt_, 3)
    sos = [(p, sp) for p in range(m) for sp in (0, 1)]
    all_dets = [[sos[k] for k in c] for c in combinations(range(2 * m), 3)]
    e_gen, _ = noci_ground_gensc(all_dets, s_, h_, g_)
    print(f"     bitstring FCI vs gensc complete space (random, 3e): |dE| = "
          f"{abs(e_bit - e_gen):.2e}")
    print("\n[validation suite complete]")


if __name__ == "__main__":
    run_validations()
