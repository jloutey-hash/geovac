"""Two-centre Coulomb-Sturmian Shibuya-Wulfman (SW) metric on a full sigma basis.

Backing for Paper 60 sec:molecular's sigma-basis paragraph and sec:obstruction's fixed-S posing
(added 2026-10-09).  The s-only sine-basis family of ``geovac/sturmian_sigma_law.py`` is the
special case ``sonly=True``; this module adds l > 0 (n <= N, l < n, m = 0 per centre).

All integrals are position-space prolate-spheroidal quadrature -- an independent route from the
momentum-space sine-basis code.  With a common scale k, centre a(t) and principal number n_t of
basis function t, the Sturmian equation gives

    S_{t't} = (n_t / k) <chi_t'| 1/r_{a(t)} |chi_t>   (SW metric; intra-centre block = identity)
    W_{t't} = (1 / k)   <chi_t'| 1/r_A + 1/r_B |chi_t>  (Wulfman potential)

Two posings of the one-electron H2+-type problem on the same basis:
  * fixed R (Paper 60 eq:sw): W C = k' S C with isoenergetic self-consistency k' = k;
  * fixed S = kR (Shibuya-Wulfman 1965; Avery & Avery 2017): the standard Hermitian
    eigenproblem K C' = k C', K = D^1/2 S D^1/2, D = diag(Z/n), and R = S/k.
"""
from math import factorial
from typing import List, Tuple

import numpy as np
from scipy.optimize import brentq
from scipy.special import eval_legendre, genlaguerre

NX, NE = 160, 160
_TX, _WX = np.polynomial.laguerre.laggauss(NX)
_EE, _WE = np.polynomial.legendre.leggauss(NE)

CANON_TOL = 1e-12  # canonical orthogonalisation threshold, relative to lambda_max(S)


def basis(N: int, sonly: bool = False) -> List[Tuple[int, int]]:
    """(n, l) labels per centre, m = 0: s-only (l = 0) or full sigma (l < n)."""
    if sonly:
        return [(n, 0) for n in range(1, N + 1)]
    return [(n, l) for n in range(1, N + 1) for l in range(n)]


def radial(n: int, l: int, k: float, r: np.ndarray) -> np.ndarray:
    """L2-normalised Coulomb-Sturmian radial function at scale k (hydrogenic, Z = n k)."""
    nrm = np.sqrt((2 * k) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    x = 2 * k * r
    return nrm * x ** l * np.exp(-k * r) * genlaguerre(n - l - 1, 2 * l + 1)(x)


def build(N: int, k: float, R: float, sonly: bool = False,
          Z: float = 1.0) -> Tuple[np.ndarray, np.ndarray, List[Tuple[int, int]]]:
    """SW metric S and Wulfman matrix W (nuclear charge Z on both centres; A at z=-R/2, B at +R/2).

    Ordering: all centre-A functions, then all centre-B functions, each in ``basis(N, sonly)``
    order.  S is symmetrised.
    """
    alpha = k * R
    xi = 1.0 + _TX[:, None] / alpha
    eta = _EE[None, :]
    wgt = (_WX * np.exp(_TX))[:, None] / alpha * _WE[None, :]
    h = R / 2.0
    r_a, r_b = h * (xi + eta), h * (xi - eta)
    c_a = np.clip((xi * eta + 1) / (xi + eta), -1, 1)
    c_b = np.clip((xi * eta - 1) / (xi - eta), -1, 1)
    vol = 2 * np.pi * h ** 3 * (xi ** 2 - eta ** 2) * wgt
    labels = basis(N, sonly)

    def ang(l: int, c: np.ndarray) -> np.ndarray:
        return np.sqrt((2 * l + 1) / (4 * np.pi)) * eval_legendre(l, c)

    funcs = ([(radial(n, l, k, r_a) * ang(l, c_a), 0, n) for (n, l) in labels]
             + [(radial(n, l, k, r_b) * ang(l, c_b), 1, n) for (n, l) in labels])
    m = len(funcs)
    S, W = np.zeros((m, m)), np.zeros((m, m))
    for i, (fi, _, _) in enumerate(funcs):
        for j, (fj, cj, nj) in enumerate(funcs):
            prod = fi * fj * vol
            v_a, v_b = np.sum(prod / r_a), np.sum(prod / r_b)
            S[i, j] = (nj / k) * (v_a if cj == 0 else v_b)
            W[i, j] = Z * (v_a + v_b) / k
    return 0.5 * (S + S.T), W, labels


def top_generalized(W: np.ndarray, S: np.ndarray) -> float:
    """Largest root of W C = k S C by canonical orthogonalisation (S may be cond ~1e15)."""
    ev, V = np.linalg.eigh(S)
    keep = ev > CANON_TOL * ev[-1]
    X = V[:, keep] / np.sqrt(ev[keep])
    return float(np.linalg.eigvalsh(X.T @ (0.5 * (W + W.T)) @ X)[-1])


def galerkin_energy(N: int, R: float) -> float:
    """H2+ electronic energy, fixed-R isoenergetic posing (Paper 60 eq:sw)."""
    def resid(kk: float) -> float:
        S, W, _ = build(N, kk, R)
        return top_generalized(W, S) - kk

    k = brentq(resid, 0.9, 2.5, xtol=1e-11)
    return -0.5 * k * k


def kernel_kmax(N: int, s: float) -> float:
    """Largest eigenvalue k of the fixed-S standard posing K = D^1/2 S D^1/2."""
    S, _, labels = build(N, 1.0, s)
    d = np.sqrt(np.array([1.0 / n for (n, _) in labels] * 2))
    return float(np.linalg.eigvalsh(d[:, None] * S * d[None, :])[-1])


def kernel_energy(N: int, R: float) -> float:
    """H2+ electronic energy, fixed-S posing, with a 1-D root search for S such that S/k = R."""
    s = brentq(lambda ss: ss / kernel_kmax(N, ss) - R, 1e-2, 40.0, xtol=1e-11)
    k = kernel_kmax(N, s)
    return -0.5 * k * k


def gerade_ungerade(S: np.ndarray, labels: List[Tuple[int, int]]) -> Tuple[np.ndarray, np.ndarray]:
    """Gerade and ungerade blocks of S (inversion maps (n,l) on A to (n,l) on B with (-1)^l)."""
    m = len(labels)
    pi = np.diag([(-1.0) ** l for (_, l) in labels])
    U = np.vstack([np.hstack([np.eye(m), np.eye(m)]), np.hstack([pi, -pi])]) / np.sqrt(2)
    T = U.T @ S @ U
    return T[:m, :m], T[m:, m:]


def bessel_deficit(N: int, s: float, sonly: bool = False) -> float:
    """1 - sum_i <chi^A_i, chi^B_1s>_S^2 (A-set is S-orthonormal), at k = 1, R = s."""
    S, _, labels = build(N, 1.0, s, sonly)
    m = len(labels)
    return float(1.0 - np.sum(S[:m, m] ** 2))
