r"""Prolate separation constant Lambda as an exact commuting label (stage A: one electron).

One electron, charges Z_A at z = -R/2 (eta = -1) and Z_B at z = +R/2 (eta = +1):

    H = (F + G) / (xi^2 - eta^2),
    F = (2/R^2) Ft,  Ft = -d_xi (xi^2-1) d_xi + m^2/(xi^2-1) - R Z_+ xi
    G = (2/R^2) Gt,  Gt = -d_eta (1-eta^2) d_eta + m^2/(1-eta^2) + R Z_- eta
    Z_+ = Z_A + Z_B,  Z_- = Z_A - Z_B.

The separation-constant operator (dimensionless, = the eta-equation eigenvalue A,
-> l(l+1) in the united-atom limit):

    Lambda = (xi^2 Gt + eta^2 Ft) / (xi^2 - eta^2).

In a PRODUCT basis f_a(xi) g_b(eta) e^{i m phi} every matrix is a Kronecker
combination of 1D blocks (the (R/2)^3 dphi factors cancel):

    W   = r2 (x) I   - r0 (x) a2          (the xi^2 - eta^2 Jacobian; metric)
    Hm  = (2/R^2) [ Fx (x) I + r0 (x) Gy ]          H c = E W c
    N   = r2 (x) Gy  + Fx (x) a2                    Lambda c = lam W c

with Fx = Kx - R Z_+ r1,  Gy = diag(l(l+1)) + R Z_- a1.  These are exactly
Atkinson's operator determinants Delta_0 = W, Delta_1 = Hm, Delta_2 = N of the
two-parameter eigenproblem, so Gamma_1 = W^-1 Hm and Gamma_2 = W^-1 N COMMUTE
EXACTLY at every finite truncation (Atkinson, Multiparameter Eigenvalue Problems,
1972) and their joint eigenvectors are decomposable x (x) y.

Blocks.  eta: normalized associated Legendre P_l^m, l = m..m+Ne-1 -- the eta
kinetic block is diag(l(l+1)) exactly and eta, eta^2 are the closed-form
(algebraic, square roots of rationals) tridiagonal / pentadiagonal matrices.
xi: f_a = (xi^2-1)^{m/2} e^{-alpha (xi-1)} L_a^{(m)}(2 alpha (xi-1)).  Every xi
integrand is polynomial x e^{-t}, evaluated with a Gauss-Laguerre rule of enough
nodes to be EXACT for that degree (no approximation; cross-checked against the
mpf moment route of geovac.prolate_recondition in debug/prolate_separation_constant_xcheck.py).
No transcendental enters the matrices beyond sqrt of rationals and the exact
Gauss nodes; alpha is a free basis scale.
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "4")

from typing import Dict, List, Tuple

import numpy as np
from scipy.linalg import cholesky, eigh, solve_triangular
from scipy.special import eval_genlaguerre, roots_laguerre


# ---------------------------------------------------------------- eta blocks (closed form)
def eta_c(l: int, m: int) -> float:
    """<l+1, m | eta | l, m> for normalized associated Legendre functions."""
    return float(np.sqrt(((l + 1) ** 2 - m * m) / ((2 * l + 1) * (2 * l + 3))))


def eta_blocks(Ne: int, m: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ls = np.arange(m, m + Ne)
    c = np.array([eta_c(int(l), m) for l in ls])          # c[i] = <l_i+1|eta|l_i>
    a1 = np.diag(c[:-1], 1) + np.diag(c[:-1], -1)
    cm = np.concatenate([[0.0], c[:-1]])                   # <l|eta|l-1>
    d2 = c ** 2 + cm ** 2                                   # full-space (eta^2)_{ll}
    o2 = c[:-1] * c[1:]                                     # (eta^2)_{l+2,l}
    a2 = np.diag(d2) + np.diag(o2[:-1], 2) + np.diag(o2[:-1], -2)
    Ly = np.diag((ls * (ls + 1)).astype(float))
    return a1, a2, Ly, ls


# ---------------------------------------------------------------- xi blocks (exact Gauss)
def xi_blocks(Nx: int, m: int, alpha: float) -> Tuple[np.ndarray, ...]:
    nq = Nx + m + 6
    t, w = roots_laguerre(nq)
    x = t / (2 * alpha)
    xi = x + 1.0
    u = x * (x + 2.0)
    Lv = np.array([eval_genlaguerre(a, m, t) for a in range(Nx)])            # L_a^{(m)}(t)
    dL = np.array([(-eval_genlaguerre(a - 1, m + 1, t) if a > 0 else 0 * t)
                   for a in range(Nx)])                                         # dL/dt
    wt = w / (2 * alpha)
    um = u ** m

    def blk(fac: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        return (A * (fac * wt)) @ B.T

    r0 = blk(um, Lv, Lv)
    r1 = blk(um * xi, Lv, Lv)
    r2 = blk(um * xi ** 2, Lv, Lv)
    if m == 0:
        qt = -alpha * Lv + 2 * alpha * dL                 # f'/e^{-alpha x}
        Kx = blk(u, qt, qt)
    else:
        q = m * xi * Lv - alpha * u * Lv + 2 * alpha * u * dL
        Kx = blk(u ** (m - 1), q, q) + m * m * blk(u ** (m - 1), Lv, Lv)
    # normalize functions to unit xi^2 moment (conditioning only)
    s = 1.0 / np.sqrt(np.diag(r2))
    S = np.outer(s, s)
    return r0 * S, r1 * S, r2 * S, Kx * S


# ---------------------------------------------------------------- operators
class OneE:
    def __init__(self, R: float, ZA: float, ZB: float, m: int, Nx: int, Ne: int,
                 alpha: float, flip_minus: bool = False) -> None:
        self.R, self.ZA, self.ZB, self.m, self.Nx, self.Ne = R, ZA, ZB, m, Nx, Ne
        r0, r1, r2, Kx = xi_blocks(Nx, m, alpha)
        a1, a2, Ly, ls = eta_blocks(Ne, m)
        Zp, Zm = ZA + ZB, ZA - ZB
        I = np.eye(Ne)
        Fx = Kx - R * Zp * r1
        Gy = Ly + R * Zm * a1
        Gy_lam = Ly + R * (-Zm if flip_minus else Zm) * a1      # fire-test variant
        self.W = np.kron(r2, I) - np.kron(r0, a2)
        self.Hm = (2.0 / R ** 2) * (np.kron(Fx, I) + np.kron(r0, Gy))
        self.N = np.kron(r2, Gy_lam) + np.kron(Fx, a2)
        self.ls = ls
        self.blocks = dict(r0=r0, r1=r1, r2=r2, Kx=Kx, a1=a1, a2=a2, Ly=Ly)

    def orthonormal(self) -> Tuple[np.ndarray, np.ndarray]:
        L = cholesky(self.W, lower=True)
        Hh = solve_triangular(L, solve_triangular(L, self.Hm, lower=True).T, lower=True).T
        Nh = solve_triangular(L, solve_triangular(L, self.N, lower=True).T, lower=True).T
        Hh = 0.5 * (Hh + Hh.T)
        Nh = 0.5 * (Nh + Nh.T)
        self._L = L
        return Hh, Nh

    def commutator(self) -> float:
        Hh, Nh = self.orthonormal()
        C = Hh @ Nh - Nh @ Hh
        return float(np.linalg.norm(C, 2) / (np.linalg.norm(Hh, 2) * np.linalg.norm(Nh, 2)))

    def joint(self, k: int = 30) -> Dict[str, np.ndarray]:
        """Lowest-k joint eigenpairs: diagonalize H, read Lambda (and its off-diag leakage)."""
        Hh, Nh = self.orthonormal()
        E, V = eigh(Hh)
        # resolve any degeneracy of H with Lambda (none expected generically)
        lamM = V[:, :k].T @ Nh @ V[:, :k]
        lam = np.diag(lamM).copy()
        off = np.abs(lamM - np.diag(lam))
        # coefficient vectors in the product basis
        C = solve_triangular(self._L.T, V[:, :k], lower=False)
        rank1 = []
        idx = getattr(self, "_idx", None)
        for i in range(k):
            cf = C[:, i]
            if idx is not None:
                cf = np.zeros(self.Nx * self.Ne)
                cf[idx] = C[:, i]
            s = np.linalg.svd(cf.reshape(self.Nx, self.Ne), compute_uv=False)
            rank1.append(s[1] / s[0])
        return dict(E=E[:k], lam=lam, off=off, C=C, rank1=np.array(rank1))


# ---------------------------------------------------------------- nodal labels
def _assoc_legendre_norm(lmax: int, m: int, eta: np.ndarray) -> np.ndarray:
    """Normalized P_l^m(eta), l = m..lmax (stable recurrence), rows by l - m."""
    out = np.zeros((lmax - m + 1, eta.size))
    # P_m^m normalized
    pmm = np.ones_like(eta) * np.sqrt(0.5)
    s = np.sqrt(np.clip(1 - eta ** 2, 0, None))
    for k in range(1, m + 1):
        pmm = -np.sqrt((2 * k + 1) / (2 * k)) * s * pmm
    out[0] = pmm
    if lmax > m:
        out[1] = np.sqrt(2 * m + 3) * eta * pmm
    for l in range(m + 2, lmax + 1):
        a = np.sqrt((4 * l * l - 1) / (l * l - m * m))
        b = np.sqrt(((l - 1) ** 2 - m * m) / (4 * (l - 1) ** 2 - 1))
        out[l - m] = a * (eta * out[l - m - 1] - b * out[l - m - 2])
    return out


def _sign_changes(v: np.ndarray, rel: float = 1e-6) -> int:
    thr = rel * np.abs(v).max()
    s = np.sign(v[np.abs(v) > thr])
    return int(np.sum(s[1:] != s[:-1]))


def nodal_labels(op: OneE, c: np.ndarray, alpha: float) -> Tuple[int, int]:
    M = c.reshape(op.Nx, op.Ne)
    U, s, Vt = np.linalg.svd(M)
    x_c, y_c = U[:, 0], Vt[0]
    eta = np.cos(np.linspace(1e-4, np.pi - 1e-4, 6000))
    P = _assoc_legendre_norm(op.m + op.Ne - 1, op.m, eta)
    y = y_c @ P
    t = np.linspace(0, 1, 6000) ** 2 * (60.0 + 4 * op.Nx)
    # positive prefactors (xi^2-1)^{m/2} e^{-t/2} dropped: they do not change signs
    xf = xfun_true(op, x_c, alpha, t) * np.exp(-t / 2)
    return _sign_changes(xf, 1e-3), _sign_changes(y, 1e-3)


def xfun_true(op: OneE, x_c: np.ndarray, alpha: float, t: np.ndarray) -> np.ndarray:
    nq = op.Nx + op.m + 6
    tq, wq = roots_laguerre(nq)
    xq = tq / (2 * alpha)
    uq = xq * (xq + 2.0)
    Lq = np.array([eval_genlaguerre(a, op.m, tq) for a in range(op.Nx)])
    r2d = (Lq ** 2 * (uq ** op.m * (xq + 1) ** 2 * wq / (2 * alpha))).sum(axis=1)
    s = 1.0 / np.sqrt(r2d)
    Lv = np.array([eval_genlaguerre(a, op.m, t) for a in range(op.Nx)])
    return (x_c * s) @ Lv


# ---------------------------------------------------------------- helpers for the stage-A driver
def eta_power(Ne: int, m: int, k: int) -> np.ndarray:
    """Matrix of eta^k on the first Ne normalized P_l^m (exact: full-space product, truncated)."""
    a1big = eta_blocks(Ne + k + 1, m)[0]
    return np.linalg.matrix_power(a1big, k)[:Ne, :Ne]


def eta_expectation(op: "OneE", c: np.ndarray) -> float:
    """<eta> in the full measure (xi^2 - eta^2) for a coefficient vector c (W-normalized)."""
    b = op.blocks
    a3 = eta_power(op.Ne, op.m, 3)
    M = np.kron(b["r2"], b["a1"]) - np.kron(b["r0"], a3)
    return float(c @ M @ c / (c @ op.W @ c))


def restrict(op: "OneE", parity: int) -> "OneE":
    """Homonuclear only: restrict to eta-parity (-1)^(l-m) = (-1)^parity (exact invariant subspace)."""
    sel = np.array([((l - op.m) % 2 == parity) for l in op.ls])
    idx = np.array([a * op.Ne + b for a in range(op.Nx) for b in range(op.Ne) if sel[b]])
    new = object.__new__(OneE)
    new.__dict__.update(op.__dict__)
    new.W = op.W[np.ix_(idx, idx)]
    new.Hm = op.Hm[np.ix_(idx, idx)]
    new.N = op.N[np.ix_(idx, idx)]
    new._idx = idx
    new._parity = parity
    new._full = op
    return new


def joint_labeled(op: "OneE", k: int, alpha: float) -> List[Dict]:
    J = op.joint(k)
    out = []
    full = getattr(op, "_full", op)
    for i in range(min(k, J["C"].shape[1])):
        c = J["C"][:, i]
        if hasattr(op, "_idx"):
            cf = np.zeros(full.Nx * full.Ne)
            cf[op._idx] = c
        else:
            cf = c
        nx, ny, dx, dy = index_labels(full, float(J["E"][i]), float(J["lam"][i]),
                                      getattr(op, "_parity", -1))
        out.append(dict(E=float(J["E"][i]), lam=float(J["lam"][i]), nxi=nx, neta=ny,
                        eta=eta_expectation(full, cf), rank1=float(J["rank1"][i]),
                        sep_dist=max(dx, dy)))
    return out


def index_labels(op: "OneE", E: float, lam: float, parity: int = -1
                 ) -> Tuple[int, int, float, float]:
    """Exact nodal labels by Klein's oscillation theorem: at the joint eigenvalue (E, lam) the
    separated 1D Sturm-Liouville problems are
        eta:  (Gy + eps a2) y = lam y                (eps = E R^2 / 2)
        xi :  (Fx - eps r2) x = -lam r0 x
    and the number of nodes of each factor equals the index of its eigenvalue.  Returns
    (n_xi, n_eta, dist_xi, dist_eta) with the distances to the matched 1D eigenvalue
    (both ~0 certify that the 2D joint eigenpair IS a separated product)."""
    b = op.blocks
    R = op.R
    Zp, Zm = op.ZA + op.ZB, op.ZA - op.ZB
    eps = E * R * R / 2.0
    Gy = b["Ly"] + R * Zm * b["a1"] + eps * b["a2"]
    Gy = 0.5 * (Gy + Gy.T)
    if parity >= 0:   # homonuclear: eta-parity block (g/u pairs are near-degenerate at large R)
        sel = np.array([((l - op.m) % 2 == parity) for l in op.ls])
        ey = np.linalg.eigvalsh(Gy[np.ix_(sel, sel)])
    else:
        ey = np.linalg.eigvalsh(Gy)
    Fx = b["Kx"] - R * Zp * b["r1"] - eps * b["r2"]
    ex = eigh(0.5 * (Fx + Fx.T), 0.5 * (b["r0"] + b["r0"].T), eigvals_only=True)
    iy = int(np.argmin(np.abs(ey - lam)))
    dy = float(abs(ey[iy] - lam))
    if parity >= 0:
        iy = 2 * iy + parity
    ix = int(np.argmin(np.abs(ex + lam)))
    return ix, iy, float(abs(ex[ix] + lam)), dy
