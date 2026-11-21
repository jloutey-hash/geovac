"""Track-1 core: ONE object -- the cross-center overlap singular spectrum {sigma_k}.

Two exact functionals of the same spectrum are tested here:

    cond(S)          = (1 + sigma_max) / (1 - sigma_max)         [metric-conditioning wall]
    ||[P_A, P_B]||   = max_k  sigma_k * sqrt(1 - sigma_k^2)      [composition wall, v4.73.0]

valid whenever the two intra-center blocks are the identity (orthonormal within a
center), so that S = [[I, C], [C^T, I]] and sigma_k = svd(C) = cos(principal angles).

Provides:
  * SW (Shibuya-Wulfman) shared-scale s-block via the exact momentum-space 1D chi
    reduction (identical integrand to debug/sturmian_sw_momentum.py, cross-checked
    here by an independent Gauss-Legendre quadrature in v = cot(chi/2)).
  * Exact two-center hydrogenic (Goscinskian a = Q/n) s-s overlaps via
    Mulliken/Ruedenberg A_m(p)/B_n(q) auxiliary integrals -- gives the CHARGE
    dependence the SW metric structurally does not have.

Diagnostic only.  No paper / CLAUDE / test edits.
"""
from __future__ import annotations

import numpy as np
from math import comb, factorial
from numpy.polynomial.legendre import leggauss

_trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))

# ======================================================================= SW blocks
# S_{n'n}(R) = (2/pi) int_0^pi sin(n' chi) sin(n chi) sinc(kR cot(chi/2)) dchi
# intra (R=0): (2/pi) int_0^pi sin(n' chi) sin(n chi) dchi = delta   [EXACT]
_M = 600001
_chi = np.linspace(1e-8, np.pi, _M)
_cot = 1.0 / np.tan(_chi / 2.0)
_1mcos = 1.0 - np.cos(_chi)
_sincache: dict = {}


def _sin(j: int) -> np.ndarray:
    if j not in _sincache:
        _sincache[j] = np.sin(j * _chi)
    return _sincache[j]


def _sinc(x: np.ndarray) -> np.ndarray:
    out = np.ones_like(x)
    nz = x != 0.0
    out[nz] = np.sin(x[nz]) / x[nz]
    return out


def sw_block_trapz(s: float, nmax: int, kind: str = "S") -> np.ndarray:
    """nmax x nmax two-center s-s block at reduced separation s = k*R (s=0 -> intra)."""
    sfac = np.ones_like(_chi) if s == 0.0 else _sinc(s * _cot)
    wfac = sfac if kind == "S" else _1mcos * sfac
    B = np.zeros((nmax, nmax))
    for a in range(1, nmax + 1):
        for b in range(a, nmax + 1):
            v = (2.0 / np.pi) * _trapz(_sin(a) * _sin(b) * wfac, _chi)
            B[a - 1, b - 1] = B[b - 1, a - 1] = v
    return B


def sw_block_gl(s: float, nmax: int, kind: str = "S", V: float = 400.0,
                per_panel: int = 24, panels_per_period: int = 4) -> np.ndarray:
    """Independent quadrature: v = cot(chi/2) in (0, inf); dchi = 2 dv/(1+v^2).

    integrand = sin(2a arccot v) sin(2b arccot v) * sinc(s v) * 2/(1+v^2)   [kind 'S']
    and an extra factor (1 - cos chi) = 2/(1+v^2)                          [kind 'm'].
    Panels aligned to the sin(s v) period; algebraic v^-5 tail beyond V.
    """
    period = 2 * np.pi / s if s > 0 else V
    h = min(period / panels_per_period, 0.5)
    edges = np.arange(0.0, V + h, h)
    xg, wg = leggauss(per_panel)
    v = np.concatenate([0.5 * (b - a) * xg + 0.5 * (a + b)
                        for a, b in zip(edges[:-1], edges[1:])])
    w = np.concatenate([0.5 * (b - a) * wg for a, b in zip(edges[:-1], edges[1:])])
    th = np.arctan2(1.0, v)                       # arccot(v)
    Smat = np.array([np.sin(2 * a * th) for a in range(1, nmax + 1)])
    sc = np.ones_like(v) if s == 0.0 else np.sin(s * v) / (s * v)
    W = sc * 2.0 / (1.0 + v ** 2)
    if kind != "S":
        W = W * 2.0 / (1.0 + v ** 2)
    return (2.0 / np.pi) * (Smat * (w * W)) @ Smat.T


def assemble_two_center(C: np.ndarray) -> np.ndarray:
    """S = [[I, C], [C^T, I]] with EXACT identity intra-blocks."""
    n = C.shape[0]
    I = np.eye(n)
    return np.block([[I, C], [C.T, I]])


# ================================================== exact hydrogenic two-center s-s
def _lag1_coeffs(n: int, a: float) -> np.ndarray:
    """coefficients c_j of  L^1_{n-1}(2 a r) = sum_j c_j r^j."""
    return np.array([(-1.0) ** i * comb(n, n - 1 - i) * (2.0 * a) ** i / factorial(i)
                     for i in range(n)], dtype=np.longdouble)


def hyd_s_coeffs(n: int, a: float) -> np.ndarray:
    """R_{n0}(r) = sum_j d_j r^j e^{-a r}, normalized (int R^2 r^2 dr = 1), a = Z/n."""
    norm = np.sqrt(np.longdouble((2 * a) ** 3) * factorial(n - 1) / (2 * n * factorial(n)))
    return norm * _lag1_coeffs(n, a)


def _A_aux(mmax: int, p) -> np.ndarray:
    """A_m(p) = int_1^inf xi^m e^{-p xi} dxi, m = 0..mmax  (all-positive sum)."""
    out = np.zeros(mmax + 1, dtype=np.longdouble)
    ep = np.exp(-p)
    for m in range(mmax + 1):
        acc = np.longdouble(0.0)
        term = np.longdouble(1.0)
        for i in range(m + 1):
            acc += term / p ** (i + 1)
            term *= (m - i)
        out[m] = ep * acc
    return out


def _B_aux(nmax_: int, q) -> np.ndarray:
    """B_n(q) = int_{-1}^{1} eta^n e^{-q eta} deta, n = 0..nmax_."""
    out = np.zeros(nmax_ + 1, dtype=np.longdouble)
    if abs(float(q)) < 0.5:
        for n in range(nmax_ + 1):
            acc = np.longdouble(0.0)
            term = np.longdouble(1.0)          # (-q)^s / s!
            for sdx in range(0, 120):
                if (n + sdx) % 2 == 0:
                    acc += term * 2.0 / (n + sdx + 1)
                term *= (-q) / (sdx + 1)
                if abs(float(term)) < 1e-45 and sdx > 4:
                    break
            out[n] = acc
    else:
        eq, emq = np.exp(q), np.exp(-q)
        out[0] = (eq - emq) / q
        for n in range(1, nmax_ + 1):
            out[n] = (((-1.0) ** n * eq - emq) / q) + (n / q) * out[n - 1]
    return out


def two_center_s_overlap(nA: int, aA: float, nB: int, aB: float, R: float) -> float:
    """<chi^A_{nA,0} | chi^B_{nB,0}>, exact Mulliken/Ruedenberg auxiliary form."""
    Rl = np.longdouble(R)
    half = Rl / 2
    p = half * (np.longdouble(aA) + np.longdouble(aB))
    q = half * (np.longdouble(aA) - np.longdouble(aB))
    dA = hyd_s_coeffs(nA, aA)
    dB = hyd_s_coeffs(nB, aB)
    mmax = (nA - 1) + (nB - 1) + 2
    A = _A_aux(mmax, p)
    B = _B_aux(mmax, q)
    tot = np.longdouble(0.0)
    for j, cj in enumerate(dA):
        for kk, ck in enumerate(dB):
            pref = cj * ck * half ** (3 + j + kk)
            sub = np.longdouble(0.0)
            for u in range(j + 2):
                cu = comb(j + 1, u)
                for vv in range(kk + 2):
                    sub += (cu * comb(kk + 1, vv) * (-1.0) ** vv
                            * A[j + kk + 2 - u - vv] * B[u + vv])
            tot += pref * sub
    return float(0.5 * tot)                    # (1/4pi) * 2pi = 1/2


def goscinskian_cross_block(nmax: int, ZA: float, ZB: float, R: float) -> np.ndarray:
    """C[i,j] = <chi^A_{i+1} | chi^B_{j+1}> for Goscinskian a = Z/n on each center."""
    C = np.zeros((nmax, nmax))
    for i in range(1, nmax + 1):
        for j in range(1, nmax + 1):
            C[i - 1, j - 1] = two_center_s_overlap(i, ZA / i, j, ZB / j, R)
    return C


# ============================================================= the two functionals
def sigma_spectrum(C: np.ndarray) -> np.ndarray:
    return np.linalg.svd(C, compute_uv=False)


def cond_from_sigma(sig: np.ndarray) -> float:
    smax = float(sig.max())
    return (1.0 + smax) / (1.0 - smax)


def commutator_from_sigma(sig: np.ndarray) -> float:
    return float(np.max(sig * np.sqrt(np.maximum(0.0, 1.0 - sig ** 2))))


def commutator_direct(C: np.ndarray) -> float:
    """||[P_A,P_B]|| from explicit orthogonal projectors in the joint span.

    Realize the 2n non-orthogonal functions as columns of X with X^T X = S; take
    X = S^{1/2} (any square root works, the projectors are basis-independent).
    """
    n = C.shape[0]
    S = assemble_two_center(C)
    w, U = np.linalg.eigh(S)
    X = U @ np.diag(np.sqrt(np.maximum(w, 0.0))) @ U.T     # X^T X = S
    QA = np.linalg.qr(X[:, :n])[0]
    QB = np.linalg.qr(X[:, n:])[0]
    PA, PB = QA @ QA.T, QB @ QB.T
    return float(np.linalg.norm(PA @ PB - PB @ PA, 2))


def canonical_correlations(SAA: np.ndarray, SAB: np.ndarray, SBB: np.ndarray) -> np.ndarray:
    """sigma_k = svd(SAA^{-1/2} SAB SBB^{-1/2}) -- principal-angle cosines for
    general (non-identity) intra blocks."""
    def invsqrt(M):
        w, U = np.linalg.eigh(M)
        return U @ np.diag(w ** -0.5) @ U.T
    return np.linalg.svd(invsqrt(SAA) @ SAB @ invsqrt(SBB), compute_uv=False)


def loglog_fit(N, y):
    """power-law fit y ~ N^p; returns (p, R2, max |ln-residual|)."""
    N = np.asarray(N, float)
    y = np.asarray(y, float)
    ln = np.log(y)
    pc = np.polyfit(np.log(N), ln, 1)
    res = ln - np.polyval(pc, np.log(N))
    r2 = 1 - np.sum(res ** 2) / np.sum((ln - ln.mean()) ** 2)
    return float(pc[0]), float(r2), float(np.abs(res).max())
