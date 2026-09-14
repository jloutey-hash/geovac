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
pre-asymptotic readings of this one law).

PRIOR ART (recorded 2026-09-11): this asymptotic is NOT derived here.  It is the
Kac-Murdock-Szego extreme-eigenvalue law (c_1 = pi^2, 1953; normal form in
Boettcher-Widom arXiv:math/0412269).  What is ours is the IDENTIFICATION -- that
the two-center Shibuya-Wulfman metric in the sine basis IS such a finite
section, Toeplitz minus Hankel with symbol j0(kR cot(chi/2)).  Do not describe
the law as derived here.  For equivalent centers the gerade
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


def theta2_band_matrix(n: int) -> np.ndarray:
    """<theta^2> in the sine band, in closed form:  no quadrature, no Bessel.

    T_ab = (2/pi) int_0^pi sin(a chi) sin(b chi) (pi - chi)^2 dchi, evaluated
    exactly as a Toeplitz-minus-Hankel matrix with c_j = 2/j^2 (j != 0),
    c_0 = pi^2/3, i.e. symbol (pi - chi)^2 = theta^2:

        T_ab = 2/(a-b)^2 - 2/(a+b)^2   (a != b),
        T_aa = pi^2/3 - 1/(2 a^2).

    Its least eigenvalue is the minimal mean-square spread of a band-limited
    function about the Fock pole, and n^2 lam_min -> pi^2 -- the
    Kac-Murdock-Szego constant obtained with no Bessel function present, which
    is what shows that the pi^2 of the conditioning law is a TRUNCATION
    constant rather than Bessel content.  Backing test:
    tests/test_paper60_contraction_window.py.
    """
    a = np.arange(1, n + 1)
    A = a[:, None]
    B = a[None, :]
    d = A - B
    s = A + B
    T = 2.0 / np.where(d == 0, 1, d) ** 2 - 2.0 / s ** 2
    np.fill_diagonal(T, np.pi ** 2 / 3.0 - 1.0 / (2.0 * a ** 2))
    return T


def weighted_blocks(s: float, n: int, weight=None, M: int = 200001
                    ) -> Tuple[np.ndarray, np.ndarray]:
    """Intra/cross blocks under a radial weight W on the Fock sphere.

    A_ab = (2/pi) int sin(a chi) sin(b chi) W(chi) dchi              (intra)
    B_ab = (2/pi) int sin(a chi) sin(b chi) W(chi) j0(s cot(chi/2)) dchi (cross)

    weight=None means W = 1, for which (A, B) reproduces (I, sw_cross_block).
    The generalized symbol is the quotient W*j0/W = j0, so the conditioning law
    is predicted to be independent of any smooth positive W -- measured, with
    the vanishing-W control, in tests/test_paper60_contraction_window.py.
    """
    chi = np.linspace(1e-9, np.pi, M)
    wq = np.full(M, chi[1] - chi[0])
    wq[0] *= 0.5
    wq[-1] *= 0.5
    Wv = np.ones_like(chi) if weight is None else np.asarray(weight(chi), float)
    x = s / np.tan(chi / 2.0)
    sym = np.ones_like(x)
    nz = x != 0.0
    sym[nz] = np.sin(x[nz]) / x[nz]
    a = np.arange(1, n + 1)
    smat = np.sin(np.outer(a, chi))
    A = (2.0 / np.pi) * (smat * (wq * Wv)) @ smat.T
    B = (2.0 / np.pi) * (smat * (wq * Wv * sym)) @ smat.T
    return A, B


def generalized_sigma_max(A: np.ndarray, B: np.ndarray) -> float:
    """Largest generalized singular value: spectral radius of A^-1/2 B A^-1/2."""
    ev, Q = np.linalg.eigh(A)
    ev = np.maximum(ev, 1e-300)
    half = Q @ np.diag(ev ** -0.5) @ Q.T
    return float(np.max(np.abs(np.linalg.eigvalsh(half @ B @ half))))


def sigma_law_report(s: float, nmax: int, M: int = 300001) -> Tuple[float, float, float]:
    """(1 - sigma_max, cond(S), ||[P_A,P_B]||) at (s, nmax) from one build."""
    sig = sigma_spectrum(sw_cross_block(s, nmax, M))
    return 1.0 - float(np.max(sig)), cond_from_sigma(sig), commutator_from_sigma(sig)
