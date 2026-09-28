r"""Paper 19, well-shape subsection, "Refinement: the atom-centered core is a
convenience, not a necessity" (Table tab:pure_product_core).

Claim backed: a PURE prolate product basis -- xi-Laguerre x eta-Legendre, one
shared exponent alpha, no atom-centered function (the Paper 11/12 basis) --
carries the Li-like 1s core on the TWO-BODY side.  For Z = 3 with the H nucleus
absent, the Galerkin 1s is exactly e^{-Z r_A} in the complete-basis limit and its
self-repulsion (cc|cc) is exactly 5Z/8 = 15/8.  At alpha = R/2 the 15x15 basis
puts (cc|cc) within 0.01 mHa of 15/8 (measured +6.2e-3 mHa at R = 3.015), and
its R-spread over R = 2.5 / 3.015 / 3.5 below 1e-3 mHa (measured 8.5e-4).

Wrong answers this rejects (fire-tested with debug/qa/fire_test.py):
  * the 10x10 value (+0.362 mHa): the basis is NOT trivially exact, so a guard
    that passes at 10x10 would be checking nothing;
  * too low a working precision.  At dps = 62 (instead of 40 + 8N = 160) the
    same computation returns +5.75e-3 mHa -- inside the 0.01 mHa gate, finite and
    plausible, but wrong.  This is Paper 19's precision trap in miniature, and
    only the pinned high-precision value catches it.

The two-body evaluation is the Paper 12 Neumann route built from geovac
primitives only (closed-form A moments and B seeds of
geovac.neumann_vee_general_m, the Laguerre/Legendre coefficients and eta moments
of geovac.prolate_recondition).  The l-sum terminates at l = 2N because the pair
densities are polynomial in eta.  The routine is reproduced here from the sprint
driver debug/pure_product_core_eri.py so that this paper-backing test does not
depend on the prunable debug/ tree (C22 check D).

Transcendentals (tagged): the truncated value is a finite Q-combination of the
Paper 12 Neumann seeds e^{-c}, e^{c}E_1(2c), Euler gamma, ln c, ln 2 at the density
rate c = 2 alpha; they cancel to the rational 15/8 only in the complete-basis limit.
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

import mpmath as mp
import numpy as np
import pytest
from scipy.linalg import eigh

from geovac import neumann_vee_general_m as ngm
from geovac import prolate_recondition as pr

N_GATE = 15
R_GATE = "3.015"
# Converged (cc|cc) - 15/8 at N = 15, alpha = R/2, R = 3.015, in mHa.  Measured
# 2026-09-26 at dps 80, 100, 160 (identical to 17 significant digits).
ERR_15_MHA = 6.2219770788e-3
PIN_TOL_MHA = 1e-6


def _dps(N: int) -> int:
    """Working precision: ~8N digits (Paper 19's precision rule)."""
    return 40 + 8 * N


# --------------------------------------------------------------- eta moments
def _Me(q: int) -> mp.mpf:
    return mp.mpf(2) / (q + 1) if q % 2 == 0 else mp.mpf(0)


_YCACHE: Dict[Tuple[int, int, int], mp.mpf] = {}


def _Ymom(l: int, Q: int) -> mp.mpf:
    """int_{-1}^{1} eta^Q P_l(eta) d eta (exact rational)."""
    key = (l, Q, mp.mp.dps)
    if key not in _YCACHE:
        if l > Q or (Q + l) % 2:
            v = mp.mpf(0)
        else:
            v = pr._mom_eta(ngm._shift(list(ngm._leg_coeffs(l)), Q))
        _YCACHE[key] = v
    return _YCACHE[key]


def _conv2(K1: Sequence[Sequence], K2: Sequence[Sequence]) -> List[List[mp.mpf]]:
    """rho[P][Q] = sum K1[j][l] K2[P-j][Q-l] (monomial product)."""
    n1, m1, n2, m2 = len(K1), len(K1[0]), len(K2), len(K2[0])
    rho = [[mp.mpf(0)] * (m1 + m2 - 1) for _ in range(n1 + n2 - 1)]
    for j in range(n1):
        for l in range(m1):
            a = K1[j][l]
            if a == 0:
                continue
            for jj in range(n2):
                rj = rho[j + jj]
                for ll in range(m2):
                    b = K2[jj][ll]
                    if b != 0:
                        rj[l + ll] += a * b
    return rho


def _overlap_from_rho(rho: List[List[mp.mpf]], c: mp.mpf, R: mp.mpf) -> mp.mpf:
    A = ngm._mono_moments(c, len(rho) + 3)
    s = mp.mpf(0)
    for P, row in enumerate(rho):
        for Q, x in enumerate(row):
            if x != 0:
                s += x * (A[P + 2] * _Me(Q) - A[P] * _Me(Q + 2))
    return 2 * mp.pi * (R / 2) ** 3 * s


# --------------------------------------------------------------- one-body core
def _core_orbital(N: int, alpha: mp.mpf, R: mp.mpf, Z: mp.mpf) -> Tuple[List[List[mp.mpf]], float]:
    """Galerkin 1s of h(Z_A = Z, Z_B = 0) in the N x N product basis.

    Returns the monomial coefficient matrix K[j][l] of the unit-normalized
    orbital divided by e^{-alpha xi}, and the orbital energy.
    """
    Cx = [pr.laguerre_coeffs(a, alpha, N) for a in range(N)]
    Ce = [pr.legendre_coeffs(b, N) for b in range(N)]
    A = ngm._mono_moments(2 * alpha, 2 * N + 4)
    Dx = []
    for a in range(N):
        p = Cx[a]
        d = [(k + 1) * p[k + 1] for k in range(N - 1)] + [mp.mpf(0)]
        Dx.append([d[k] - alpha * p[k] for k in range(N)])
    De = []
    for b in range(N):
        p = Ce[b]
        De.append([(k + 1) * p[k + 1] for k in range(N - 1)] + [mp.mpf(0)])

    def xmat(P, Qm, shift: int, minus1: bool = False) -> np.ndarray:
        M = np.zeros((N, N))
        for a in range(N):
            for b in range(a, N):
                s = mp.mpf(0)
                for i in range(N):
                    if P[a][i] == 0:
                        continue
                    for j in range(N):
                        if Qm[b][j] == 0:
                            continue
                        t = A[i + j + shift]
                        if minus1:
                            t -= A[i + j]
                        s += P[a][i] * Qm[b][j] * t
                M[a, b] = M[b, a] = float(s)
        return M

    def emat(P, Qm, shift: int, oneminus: bool = False) -> np.ndarray:
        M = np.zeros((N, N))
        for a in range(N):
            for b in range(a, N):
                s = mp.mpf(0)
                for i in range(N):
                    if P[a][i] == 0:
                        continue
                    for j in range(N):
                        if Qm[b][j] == 0:
                            continue
                        t = _Me(i + j + shift)
                        if oneminus:
                            t = _Me(i + j) - _Me(i + j + 2)
                        s += P[a][i] * Qm[b][j] * t
                M[a, b] = M[b, a] = float(s)
        return M

    Sx0, Sx1, Sx2 = (xmat(Cx, Cx, k) for k in (0, 1, 2))
    Tx = xmat(Dx, Dx, 2, minus1=True)
    Ey0, Ey1, Ey2 = (emat(Ce, Ce, k) for k in (0, 1, 2))
    Ty = emat(De, De, 0, oneminus=True)
    h3 = float((R / 2) ** 3)
    tp = float(2 * mp.pi)
    Rf, Zf = float(R), float(Z)
    S = h3 * tp * (np.kron(Sx2, Ey0) - np.kron(Sx0, Ey2))
    T = 0.5 * (4 / Rf ** 2) * h3 * tp * (np.kron(Tx, Ey0) + np.kron(Sx0, Ty))
    V = -Zf * (Rf / 2) ** 2 * tp * (np.kron(Sx1, Ey0) - np.kron(Sx0, Ey1))
    H = T + V
    dS = np.sqrt(np.diag(S))
    w, U = eigh(H / np.outer(dS, dS), S / np.outer(dS, dS))
    v = (U[:, 0] / dS).reshape(N, N)
    vm = [[mp.mpf(float(v[a, b])) for b in range(N)] for a in range(N)]
    tmp = [[sum((vm[a][b] * Ce[b][l] for b in range(N)), mp.mpf(0)) for l in range(N)]
           for a in range(N)]
    K = [[sum((Cx[a][j] * tmp[a][l] for a in range(N)), mp.mpf(0)) for l in range(N)]
         for j in range(N)]
    nrm = _overlap_from_rho(_conv2(K, K), 2 * alpha, R)
    sc = 1 / mp.sqrt(nrm)
    K = [[sc * x for x in row] for row in K]
    return K, float(w[0])


# --------------------------------------------------------------- Neumann contraction
def _Gvec(l: int, P: int, c: mp.mpf) -> List[mp.mpf]:
    """G_P[n] = sum_{j>=n} W_j j!/(n! c^{j-n+1}),  W = xi^P P_l(xi)."""
    W = ngm._W_poly(l, 0, 0, P)
    deg = len(W) - 1
    G = [mp.mpf(0)] * (deg + 1)
    ic = 1 / c
    for n in range(deg + 1):
        s = mp.mpf(0)
        f = ic
        for j in range(n, deg + 1):
            if j > n:
                f = f * j * ic
            if W[j] != 0:
                s += W[j] * f
        G[n] = s
    return G


def _eri_self(rho: List[List[mp.mpf]], c: mp.mpf, R: mp.mpf) -> mp.mpf:
    """(rho|rho) for one m = 0 pair density rho(xi, eta) e^{-c xi}."""
    l_hi = len(rho[0]) + 1
    p_max = len(rho) + 1
    p_corr_max = 2 * p_max + l_hi
    Am = ngm._mono_moments(c, p_max + l_hi + 2)
    B = ngm._B_table(0, 0, l_hi, p_max, c)
    Bs = ngm._B_table(0, 0, l_hi, p_corr_max, 2 * c)

    def wvec(l: int) -> List[mp.mpf]:
        w = [mp.mpf(0)] * (p_max + 1)
        for P in range(p_max + 1):
            s = mp.mpf(0)
            if 0 <= P - 2 < len(rho):
                for Q, x in enumerate(rho[P - 2]):
                    if x != 0:
                        s += x * _Ymom(l, Q)
            if P < len(rho):
                for Q, x in enumerate(rho[P]):
                    if x != 0:
                        s -= x * _Ymom(l, Q + 2)
            w[P] = s
        return w

    total = mp.mpf(0)
    for l in range(l_hi + 1):
        w = wvec(l)
        if all(x == 0 for x in w):
            continue
        Al = [ngm._A_moment(l, 0, 0, P, Am) for P in range(p_max + 1)]
        wA = sum(w[P] * Al[P] for P in range(p_max + 1))
        wB = sum(w[P] * B[(l, P)] for P in range(p_max + 1))
        bw = [sum(w[P] * Bs[(l, P + n)] for P in range(p_max + 1))
              for n in range(p_max + l + 1)]
        corr = mp.mpf(0)
        for P in range(p_max + 1):
            if w[P] != 0:
                G = _Gvec(l, P, c)
                corr += w[P] * sum(G[n] * bw[n] for n in range(len(G)))
        total += (2 * l + 1) * 2 * (wA * wB - corr)
    return (2 / R) * (R / 2) ** 6 * (2 * mp.pi) ** 2 * total


def cc_error_mha(N: int, R: str, dps: int) -> Tuple[float, float]:
    """((cc|cc) - 15/8 in mHa, 1s energy error in Ha) for Z = 3, alpha = R/2."""
    old = mp.mp.dps
    try:
        mp.mp.dps = dps
        _YCACHE.clear()
        Rm = mp.mpf(R)
        alpha = Rm / 2
        Z = mp.mpf(3)
        K, E = _core_orbital(N, alpha, Rm, Z)
        J = _eri_self(_conv2(K, K), 2 * alpha, Rm)
        return float((J - mp.mpf(15) / 8) * 1000), E + 4.5
    finally:
        mp.mp.dps = old
        _YCACHE.clear()


# --------------------------------------------------------------- tests
@pytest.fixture(scope="module")
def err15() -> Tuple[float, float]:
    return cc_error_mha(N_GATE, R_GATE, _dps(N_GATE))


def test_one_body_core_is_resolved(err15):
    """The Galerkin 1s energy is -Z^2/2 = -4.5 to well below the two-body error."""
    _, e1 = err15
    assert 0.0 <= e1 < 1e-6


def test_cc_within_gate_at_15x15(err15):
    """(cc|cc) - 15/8 < 0.01 mHa at 15x15, alpha = R/2, R = 3.015 (Paper 19)."""
    err, _ = err15
    assert abs(err) < 0.01


def test_cc_matches_converged_value(err15):
    """Pinned to the high-precision value; rejects a too-low working precision."""
    err, _ = err15
    assert abs(err - ERR_15_MHA) < PIN_TOL_MHA


def test_smaller_basis_is_not_converged():
    """The gate discriminates: 10x10 is off by +0.36 mHa (Table tab:pure_product_core)."""
    err, _ = cc_error_mha(10, R_GATE, _dps(10))
    assert 0.3 < err < 0.4


@pytest.mark.slow
def test_r_spread_at_15x15():
    """alpha proportional to R removes the R-dependence: spread < 1e-3 mHa at 15x15."""
    errs = [cc_error_mha(N_GATE, R, _dps(N_GATE))[0] for R in ("2.5", "3.015", "3.5")]
    assert all(abs(e) < 0.01 for e in errs)
    assert max(errs) - min(errs) < 1e-3
