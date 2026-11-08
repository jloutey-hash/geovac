"""Single-centre coupled-channel solver for one electron (shared by the MF and CE registrations).

Records: docs/predictions/2026-10-20_molecular_fns.md (MF1-MF4) and
docs/predictions/2026-10-20_core_eft.md (CE1-CE4), both committed before computing.

Nucleus d sits at the origin.  The wave function is psi = sum_l (u_l(r)/r) Y_l0, m = 0.  On the grid
t = ln r (uniform), u = e^{t/2} w, the radial operator of channel l becomes
    -1/2 w'' + ((l + 1/2)^2 / 2) w + e^{2t} V w = E e^{2t} w,
a symmetric generalized problem A w = E B w with B = diag(e^{2t}).  Off-centre point charges enter
through their multipoles, coupled by Gaunt coefficients <Y_l0|P_lam|Y_l'0>.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
from scipy.special import eval_legendre
from sympy.physics.wigner import wigner_3j

_FD6 = np.array([2 / 180, -27 / 180, 270 / 180, -490 / 180, 270 / 180, -27 / 180, 2 / 180])


def gaunt_m0(l: int, lam: int, lp: int) -> float:
    """<Y_l0 | P_lam | Y_l'0>."""
    if (l + lam + lp) % 2 or lam < abs(l - lp) or lam > l + lp:
        return 0.0
    return float(np.sqrt((2 * l + 1) * (2 * lp + 1)) * float(wigner_3j(l, lam, lp, 0, 0, 0)) ** 2)


@dataclass
class Grid:
    t0: float
    t1: float
    h: float
    t: np.ndarray = field(init=False)
    r: np.ndarray = field(init=False)

    anchor: Optional[float] = None                      # a radius that must sit on a grid point

    def __post_init__(self) -> None:
        if self.anchor is not None:
            ta = float(np.log(self.anchor))
            self.t0 = ta - self.h * np.ceil((ta - self.t0) / self.h)
        n = int(round((self.t1 - self.t0) / self.h))
        self.t = self.t0 + self.h * np.arange(1, n)       # interior points; w = 0 at both ends
        self.r = np.exp(self.t)

    @property
    def n(self) -> int:
        return len(self.t)

    def snap(self, r0: float) -> float:
        """Nearest grid radius to r0 (to put a kink on a grid point)."""
        return float(self.r[np.argmin(abs(np.log(self.r) - np.log(r0)))])


def d2_matrix(n: int, h: float) -> sps.csr_matrix:
    """Sixth-order second derivative with zero Dirichlet values outside the interior."""
    diags = [np.full(n - abs(k), _FD6[k + 3]) for k in range(-3, 4)]
    return sps.diags(diags, list(range(-3, 4)), shape=(n, n), format="csr") / h**2


@dataclass
class System:
    """One electron: nucleus Zd at the origin, point charges (Z, R) on the +z axis, an optional extra
    spherical local potential, per-channel extra local potentials, and an optional per-channel
    non-local operator (callable acting on a channel's w)."""
    grid: Grid
    lmax: int
    Zd: float
    others: List[Tuple[float, float]] = field(default_factory=list)   # (charge, distance)
    vsph: Optional[np.ndarray] = None                                  # extra spherical V(r)
    vchan: Dict[int, np.ndarray] = field(default_factory=dict)         # extra V_l(r) in channel l
    parity: Optional[int] = None                                       # None, 0 (even l), 1 (odd l)

    def channels(self) -> List[int]:
        ls = list(range(self.lmax + 1))
        if self.parity is not None:
            ls = [l for l in ls if l % 2 == self.parity]
        return ls

    def coupling(self) -> Dict[Tuple[int, int], np.ndarray]:
        """V_{l l'}(r) of the off-centre charges (local in r), on the grid."""
        r = self.grid.r
        ls = self.channels()
        out: Dict[Tuple[int, int], np.ndarray] = {}
        for (Z, R) in self.others:
            rl, rg = np.minimum(r, R), np.maximum(r, R)
            for i, l in enumerate(ls):
                for lp in ls[i:]:
                    acc = np.zeros_like(r)
                    for lam in range(abs(l - lp), l + lp + 1, 2):
                        g = gaunt_m0(l, lam, lp)
                        if g:
                            acc += -Z * g * rl**lam / rg ** (lam + 1)
                    out[(l, lp)] = out.get((l, lp), 0) + acc
        return out

    def local_diag(self, l: int) -> np.ndarray:
        r = self.grid.r
        v = -self.Zd / r
        if self.vsph is not None:
            v = v + self.vsph
        if l in self.vchan:
            v = v + self.vchan[l]
        return v

    def matrices(self) -> Tuple[sps.csr_matrix, sps.csr_matrix]:
        """A and B, t-major ordering (all channels at one t are contiguous)."""
        g = self.grid
        n, ls = g.n, self.channels()
        nc = len(ls)
        e2t = np.exp(2 * g.t)
        D2 = d2_matrix(n, g.h).tocoo()
        rows, cols, vals = [], [], []
        idx = lambda i, c: i * nc + c                                  # noqa: E731
        for c, l in enumerate(ls):
            rows += list(idx(D2.row, c)); cols += list(idx(D2.col, c)); vals += list(-0.5 * D2.data)
            diag = 0.5 * (l + 0.5) ** 2 + e2t * self.local_diag(l)
            ii = idx(np.arange(n), c)
            rows += list(ii); cols += list(ii); vals += list(diag)
        cpl = self.coupling()
        for (l, lp), v in cpl.items():
            c, cp = ls.index(l), ls.index(lp)
            ii, jj = idx(np.arange(n), c), idx(np.arange(n), cp)
            rows += list(ii); cols += list(jj); vals += list(e2t * v)
            if c != cp:
                rows += list(jj); cols += list(ii); vals += list(e2t * v)
        N = n * nc
        A = sps.csr_matrix((vals, (rows, cols)), shape=(N, N))
        B = sps.diags(np.repeat(e2t, nc), format="csr")
        return A, B

    # ---- helpers on t-major vectors -------------------------------------------------------------
    def split(self, x: np.ndarray) -> Dict[int, np.ndarray]:
        nc = len(self.channels())
        X = x.reshape(self.grid.n, nc)
        return {l: X[:, c] for c, l in enumerate(self.channels())}

    def join(self, d: Dict[int, np.ndarray]) -> np.ndarray:
        nc = len(self.channels())
        X = np.zeros((self.grid.n, nc))
        for c, l in enumerate(self.channels()):
            if l in d:
                X[:, c] = d[l]
        return X.reshape(-1)

    def inner(self, x: np.ndarray, y: np.ndarray) -> float:
        """<x|y> = sum_l int e^{2t} w_l^x w_l^y dt (trapezoid on the uniform t grid)."""
        e2t = np.repeat(np.exp(2 * self.grid.t), len(self.channels()))
        return float(np.sum(e2t * x * y) * self.grid.h)

    def apply_V(self, x: np.ndarray) -> np.ndarray:
        """V x for the full local potential (diagonal parts + multipole coupling), as w-values."""
        d = self.split(x)
        out = {l: self.local_diag(l) * d[l] for l in d}
        for (l, lp), v in self.coupling().items():
            out[l] = out[l] + v * d[lp]
            if l != lp:
                out[lp] = out[lp] + v * d[l]
        return self.join(out)


def ground_state(sysm: System, sigma: float, k: int = 1, extra: Optional[sps.spmatrix] = None
                 ) -> Tuple[np.ndarray, np.ndarray]:
    A, B = sysm.matrices()
    if extra is not None:
        A = A + extra
    vals, vecs = spla.eigsh(A, k=k, M=B, sigma=sigma, which="LM")
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    for j in range(vecs.shape[1]):
        v = vecs[:, j]
        v /= np.sqrt(sysm.inner(v, v))
        d = sysm.split(v)
        if d[sysm.channels()[0]][len(v) // (2 * len(sysm.channels()))] < 0:
            v *= -1
        vecs[:, j] = v
    return vals, vecs


def reduced_solve(sysm: System, E: float, psi: np.ndarray, q: np.ndarray,
                  extra: Optional[sps.spmatrix] = None) -> np.ndarray:
    """eta with (H - E) eta = -Q q and <psi|eta> = 0, by a bordered system (q as w-values)."""
    A, B = sysm.matrices()
    if extra is not None:
        A = A + extra
    M = (A - E * B).tocsr()
    Bpsi = B @ psi
    qQ = q - psi * sysm.inner(psi, q)
    rhs = np.concatenate([-(B @ qQ), [0.0]])
    K = sps.bmat([[M, Bpsi[:, None]], [Bpsi[None, :], None]], format="csc")
    sol = spla.spsolve(K, rhs)
    return sol[:-1]


def value_at_origin(sysm: System, x: np.ndarray, rlo: float = 1e-8, rhi: float = 1e-6) -> float:
    """psi(0) from the s channel, psi = u_0/r Y00 = e^{-t/2} w_0 / sqrt(4 pi), fitted on [rlo, rhi].
    The window stays far from the inner grid edge, where the Dirichlet condition mixes in a small
    irregular solution (relative size e^{t0 - t})."""
    w0 = sysm.split(x)[0]
    r = sysm.grid.r
    m = (r >= rlo) & (r <= rhi)
    vals = np.exp(-0.5 * np.log(r[m])) * w0[m] / np.sqrt(4 * np.pi)
    p = np.polyfit(r[m] / rhi, vals, 3)
    return float(np.polyval(p, 0.0))


def value_on_axis(sysm: System, x: np.ndarray, R: float) -> float:
    """psi(R z-hat) = sum_l u_l(R)/R * sqrt((2l+1)/(4 pi)), with R on (or interpolated to) the grid."""
    from scipy.interpolate import CubicSpline
    d = sysm.split(x)
    tR = np.log(R)
    out = 0.0
    for l, w in d.items():
        wR = CubicSpline(sysm.grid.t, w)(tR)
        uR = np.exp(0.5 * tR) * wR
        out += uR / R * np.sqrt((2 * l + 1) / (4 * np.pi))
    return float(out)

def reduced_solve_multi(sysm: System, E: float, vs: List[np.ndarray], q: np.ndarray) -> np.ndarray:
    """eta with (H - E) eta = -Q q, Q = 1 - sum |v><v| over the (orthonormal) vs, and <v|eta> = 0."""
    A, B = sysm.matrices()
    M = (A - E * B).tocsr()
    qQ = q.copy()
    for v in vs:
        qQ = qQ - v * sysm.inner(v, q)
    Bv = [B @ v for v in vs]
    k = len(vs)
    col = sps.csr_matrix(np.column_stack(Bv))
    K = sps.bmat([[M, col], [col.T, None]], format="csc")
    rhs = np.concatenate([-(B @ qQ), np.zeros(k)])
    return spla.spsolve(K, rhs)[:-k]

