"""The prolate separation constant Lambda is an exact commuting label of the one-electron
two-centre problem, interpolating united-atom L^2 and separated-atom Runge-Lenz A_z.

Backs the stage-A result of the separation-constant sprint (Paper 11's prolate geometry;
PM to capture).  One electron, charges Z_A at z = -R/2 (eta = -1), Z_B at z = +R/2:

    H      = (F + G)/(xi^2 - eta^2)
    Lambda = (xi^2 Gt + eta^2 Ft)/(xi^2 - eta^2)      (the eta-equation eigenvalue)
    Ft = -d(xi^2-1)d + m^2/(xi^2-1) - R Z_+ xi,   Gt = -d(1-eta^2)d + m^2/(1-eta^2) + R Z_- eta,
    F, G = (2/R^2) Ft, Gt.

In a product basis these are Atkinson's operator determinants (W, Hm, N) of a two-parameter
eigenproblem, so W^-1 Hm and W^-1 N commute EXACTLY at every truncation.

Checks and the wrong answer each rejects (fire-tested with debug/qa/fire_test.py, see the
sprint memo):
  (a) commutation to rounding, heteronuclear (3,1) and homonuclear (1,1)
      -- rejects a Lambda with the wrong sign on the Z_A - Z_B term (commutator ~1e-3..1e-2).
  (b) united-atom limit R -> 0: lambda -> l(l+1), E -> -(Z_A+Z_B)^2/2n^2, with (n, l) read off
      the exact nodal labels (Klein oscillation theorem), never from lambda itself
      -- rejects a mis-normalized Lambda (e.g. 2 x Lambda, which still commutes).
  (c) separated-atom limit R -> infinity: (lambda + p^2)/R - Z_other -> (Z_self/n) k with k in
      the spectrum of the SO(4) Runge-Lenz A_z on shell n, built from its reduced matrix elements
      (tests/test_graph_s3_nonconvergence.so4_generators; the target integers are never typed)
      -- rejects Lambda + (R^3/40) H, which commutes AND passes (b) but is not Lambda.
"""
from __future__ import annotations

import os
import sys

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

from typing import Dict, Tuple

import numpy as np
import pytest
from scipy.linalg import eigh
from scipy.special import eval_genlaguerre, roots_laguerre

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_graph_s3_nonconvergence import so4_generators  # noqa: E402


# ------------------------------------------------------------------ 1D blocks
def _eta_blocks(Ne: int, m: int) -> Dict[str, np.ndarray]:
    """Normalized associated Legendre P_l^m, l = m..m+Ne-1: closed-form eta, eta^2, L^2."""
    ls = np.arange(m, m + Ne + 2)
    c = np.sqrt(((ls + 1) ** 2 - m * m) / ((2 * ls + 1) * (2 * ls + 3)))
    a1big = np.diag(c[:-1], 1) + np.diag(c[:-1], -1)
    a1 = a1big[:Ne, :Ne]
    a2 = (a1big @ a1big)[:Ne, :Ne]
    Ly = np.diag((ls[:Ne] * (ls[:Ne] + 1)).astype(float))
    return dict(a1=a1, a2=a2, Ly=Ly)


def _xi_blocks(Nx: int, m: int, alpha: float) -> Dict[str, np.ndarray]:
    """(xi^2-1)^{m/2} e^{-alpha(xi-1)} L_a^{(m)}(2 alpha (xi-1)); polynomial x e^{-t}
    integrands, so a Gauss-Laguerre rule with Nx+m+6 nodes is exact."""
    t, w = roots_laguerre(Nx + m + 6)
    x = t / (2 * alpha)
    xi, u, wt = x + 1.0, x * (x + 2.0), w / (2 * alpha)
    L = np.array([eval_genlaguerre(a, m, t) for a in range(Nx)])
    dL = np.array([-eval_genlaguerre(a - 1, m + 1, t) if a else 0 * t for a in range(Nx)])

    def blk(f: np.ndarray, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        return (A * (f * wt)) @ B.T

    r0, r1, r2 = (blk(u ** m * xi ** k, L, L) for k in (0, 1, 2))
    if m == 0:
        q = -alpha * L + 2 * alpha * dL
        Kx = blk(u, q, q)
    else:
        q = m * xi * L - alpha * u * L + 2 * alpha * u * dL
        Kx = blk(u ** (m - 1), q, q) + m * m * blk(u ** (m - 1), L, L)
    s = 1.0 / np.sqrt(np.diag(r2))
    S = np.outer(s, s)
    return dict(r0=r0 * S, r1=r1 * S, r2=r2 * S, Kx=Kx * S)


def build(R: float, ZA: float, ZB: float, m: int, Nx: int, Ne: int, alpha: float
          ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, Dict]:
    x, y = _xi_blocks(Nx, m, alpha), _eta_blocks(Ne, m)
    Zp, Zm = ZA + ZB, ZA - ZB
    I = np.eye(Ne)
    Fx = x["Kx"] - R * Zp * x["r1"]
    Gy = y["Ly"] + R * Zm * y["a1"]
    W = np.kron(x["r2"], I) - np.kron(x["r0"], y["a2"])
    Hm = (2.0 / R ** 2) * (np.kron(Fx, I) + np.kron(x["r0"], Gy))
    Gy_L = y["Ly"] + R * Zm * y["a1"]                       # Lambda's own eta operator
    N = np.kron(x["r2"], Gy_L) + np.kron(Fx, y["a2"])
    return W, Hm, N, dict(x=x, y=y, Fx=Fx, Gy=Gy)


def _orthonormal(W: np.ndarray, M: np.ndarray) -> np.ndarray:
    w, U = np.linalg.eigh(W)
    X = U / np.sqrt(w)
    out = X.T @ M @ X
    return 0.5 * (out + out.T)


def joint(R: float, ZA: float, ZB: float, m: int, Nx: int, Ne: int, alpha: float, k: int):
    W, Hm, N, b = build(R, ZA, ZB, m, Nx, Ne, alpha)
    Hh, Nh = _orthonormal(W, Hm), _orthonormal(W, N)
    E, V = np.linalg.eigh(Hh)
    lam = np.einsum("ik,ij,jk->k", V[:, :k], Nh, V[:, :k])
    labels = []
    for Ei, li in zip(E[:k], lam):
        eps = Ei * R * R / 2.0
        ey = np.linalg.eigvalsh(b["Gy"] + eps * b["y"]["a2"])
        Fx = b["Fx"] - eps * b["x"]["r2"]
        ex = eigh(0.5 * (Fx + Fx.T), b["x"]["r0"], eigvals_only=True)
        labels.append((int(np.argmin(np.abs(ex + li))), int(np.argmin(np.abs(ey - li)))))
    return E[:k], lam, labels


# ------------------------------------------------------------------ (a) commutation
@pytest.mark.parametrize("Z", [(3.0, 1.0), (1.0, 1.0)], ids=["LiH-like_3_1", "H2plus_1_1"])
@pytest.mark.parametrize("R", [0.5, 2.0, 8.0])
@pytest.mark.parametrize("m", [0, 1])
def test_lambda_commutes_with_H_exactly(Z: Tuple[float, float], R: float, m: int) -> None:
    W, Hm, N, _ = build(R, Z[0], Z[1], m, 14, 16, max(0.3, 0.5 * R))
    Hh, Nh = _orthonormal(W, Hm), _orthonormal(W, N)
    rel = np.linalg.norm(Hh @ Nh - Nh @ Hh, 2) / (np.linalg.norm(Hh, 2) * np.linalg.norm(Nh, 2))
    assert rel < 1e-11, f"relative ||[H, Lambda]|| = {rel:.2e}"
    # non-trivial: Lambda is not a function of H (it separates states H alone orders)
    lam = np.linalg.eigvalsh(Nh)
    assert np.ptp(lam) > 1.0


# ------------------------------------------------------------------ (b) united atom
def test_united_atom_limit_is_L2() -> None:
    ZA, ZB, ZUA = 3.0, 1.0, 4.0
    R = 0.01
    seen = set()
    for m in (0, 1):
        E, lam, labels = joint(R, ZA, ZB, m, 16, 12, 0.005, 6)
        for Ei, li, (nx, ne) in zip(E, lam, labels):
            l = ne + m
            n = nx + ne + m + 1
            assert abs(li - l * (l + 1)) < 5e-3, (m, nx, ne, li)
            assert abs(Ei + ZUA ** 2 / (2 * n * n)) < 2e-2 * ZUA ** 2 / (2 * n * n), (n, Ei)
            seen.add((n, l, m))
    # the check must have seen s, p and d united-atom states
    assert {l for (_, l, _) in seen} >= {0, 1, 2}


# ------------------------------------------------------------------ (c) separated atoms
def _az_spec(n: int, m: int) -> np.ndarray:
    _, As, states = so4_generators(n)
    idx = [i for i, (nn, _, mm) in enumerate(states) if nn == n and mm == m]
    return np.linalg.eigvalsh(As[2].real[np.ix_(idx, idx)])


def test_separated_atom_limit_is_runge_lenz_Az() -> None:
    ZA, ZB, m = 3.0, 1.0, 0
    Rs = np.array([20.0, 40.0, 80.0])
    runs = [joint(R, ZA, ZB, m, 16, 36, 0.5 * R, 4) for R in Rs]
    hits = 0
    for s_i in range(4):
        lab = runs[0][2][s_i]
        E = np.array([r[0][s_i] for r in runs])
        lam = np.array([r[1][s_i] for r in runs])
        assert all(r[2][s_i] == lab for r in runs), "state order changed across R"
        y = lam - Rs ** 2 * E / 2.0                               # lambda + p^2
        s, _, _ = np.linalg.solve(np.vstack([Rs, np.ones(3), 1 / Rs]).T, y)
        # lowest four m=0 states of (3,1) at large R all sit on Z_A = 3 (n = 1, 2, 2, 3)
        Eat = E[-1] + ZB / Rs[-1]
        n = int(round(ZA / np.sqrt(-2 * Eat)))
        k = (s - ZB) * n / ZA
        spec = _az_spec(n, m)
        assert np.min(np.abs(spec - k)) < 0.02, (s_i, n, k, spec)
        hits += int(n > 1)
    assert hits >= 2          # at least two states with a non-trivial (n > 1) A_z spectrum
