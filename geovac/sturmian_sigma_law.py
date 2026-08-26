"""Principal-angle law for the two-center Shibuya-Wulfman (SW) metric.

One object -- the singular spectrum {sigma_k} of the cross-center overlap block
C -- carries both structural "walls" of the two-center shared-scale problem:

    cond(S)        = (1 + sigma_max) / (1 - sigma_max)        (metric conditioning)
    ||[P_A, P_B]|| = max_k sigma_k * sqrt(1 - sigma_k^2)      (composition wall, v4.73.0)

valid because the SW intra-center block is EXACTLY the identity,
(2/pi) * int_0^pi sin(a chi) sin(b chi) d chi = delta_ab, so
S = [[I, C], [C^T, I]] and spec(S) = {1 +/- sigma_k}.

In the sine basis the cross-block is the finite section of a multiplication
operator, C_ab = <e_a, M_W e_b> with symbol W_s(chi) = j0(s * cot(chi/2)),
s = k*R.  Since sup|W| = W(pi) = 1, sigma_max -> 1 with the band-limited
concentration rate

    1 - sigma_max = c_sym * pi^2 / n^2 + o(n^-2),   c_sym = s^2 / 24  (SW),

i.e. the conditioning exponent is exactly 2 (asymptotically; finite windows fit
lower slopes such as the N^1.85 / N^1.97 reported in Paper 60, which are
pre-asymptotic readings of this one law).  For equivalent centers the gerade
block is I + C, whose conditioning converges to the R- and N-independent
constant 2 / (1 + min_x j0(x)) = 2.555041...

Backing tests: tests/test_paper60_sigma_law.py.  Sprint chronicle: CHANGELOG
v4.103.0; drivers debug/aha_t1_*.py.
"""
from __future__ import annotations

from typing import Tuple

import numpy as np


def sw_cross_block(s: float, nmax: int, M: int = 300001) -> np.ndarray:
    """SW two-center s-s cross block at reduced separation s = k*R.

    C_ab = (2/pi) int_0^pi sin(a chi) sin(b chi) j0(s cot(chi/2)) dchi,
    a, b = 1..nmax, by uniform trapezoid on M points (vectorized).
    """
    chi = np.linspace(1e-8, np.pi, M)
    w = np.full(M, chi[1] - chi[0])
    w[0] *= 0.5
    w[-1] *= 0.5
    x = s / np.tan(chi / 2.0)
    symbol = np.ones_like(x)
    nz = x != 0.0
    symbol[nz] = np.sin(x[nz]) / x[nz]
    a = np.arange(1, nmax + 1)
    smat = np.sin(np.outer(a, chi))
    return (2.0 / np.pi) * (smat * (w * symbol)) @ smat.T


def assemble_two_center(C: np.ndarray) -> np.ndarray:
    """S = [[I, C], [C^T, I]] (intra-center blocks exactly identity)."""
    n = C.shape[0]
    I = np.eye(n)
    return np.block([[I, C], [C.T, I]])


def sigma_spectrum(C: np.ndarray) -> np.ndarray:
    """Singular values of the cross block = cosines of the principal angles."""
    return np.linalg.svd(C, compute_uv=False)


def cond_from_sigma(sig: np.ndarray) -> float:
    """cond(S) = (1 + sigma_max)/(1 - sigma_max)."""
    smax = float(np.max(sig))
    return (1.0 + smax) / (1.0 - smax)


def commutator_from_sigma(sig: np.ndarray) -> float:
    """||[P_A, P_B]||_2 = max_k sigma_k sqrt(1 - sigma_k^2)."""
    return float(np.max(sig * np.sqrt(np.maximum(0.0, 1.0 - sig ** 2))))


def commutator_direct(C: np.ndarray) -> float:
    """||[P_A, P_B]||_2 from explicitly constructed orthogonal projectors.

    Coordinates: Cholesky-factor the ambient Gram S = L L^T = V^T V with
    V = L^T; the basis vectors' ambient coordinates are the columns of V, so
    the two center subspaces are the spans of the first/last n columns of V.
    """
    S = assemble_two_center(C)
    V = np.linalg.cholesky(S).T
    n = C.shape[0]
    QA, _ = np.linalg.qr(V[:, :n])
    QB, _ = np.linalg.qr(V[:, n:])
    PA = QA @ QA.T
    PB = QB @ QB.T
    return float(np.linalg.norm(PA @ PB - PB @ PA, 2))


def predicted_one_minus_sigma_max(s: float, n: int) -> float:
    """Leading asymptote 1 - sigma_max = (s^2/24) * pi^2 / n^2 (per-center n)."""
    return (s ** 2 / 24.0) * np.pi ** 2 / n ** 2


#: single-parameter collapse constant: (1 - sigma_max) * (n/s)^2 -> pi^2/24
COLLAPSE_CONSTANT: float = np.pi ** 2 / 24.0


def gerade_constant() -> float:
    """cond(I + C) -> 2/(1 + min_x j0(x)), R- and N-independent.

    min_x sin(x)/x = -0.2172336... at tan x = x (x ~ 4.49341).
    """
    from scipy.optimize import minimize_scalar

    res = minimize_scalar(lambda x: np.sin(x) / x, bounds=(np.pi, 2 * np.pi),
                          method="bounded")
    return 2.0 / (1.0 + res.fun)


def sigma_law_report(s: float, nmax: int, M: int = 300001) -> Tuple[float, float, float]:
    """(1 - sigma_max, cond(S), ||[P_A,P_B]||) at (s, nmax) from one build."""
    sig = sigma_spectrum(sw_cross_block(s, nmax, M))
    return 1.0 - float(np.max(sig)), cond_from_sigma(sig), commutator_from_sigma(sig)
