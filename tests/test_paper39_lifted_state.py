"""Paper 39 main theorem (lifted-state route, Lemma lifted-state-T): backing test.

The theorem is Paper 38's thm:main_unconditional lifted factor-wise to
SU(2) x SU(2): d_GH <= gamma^(ab) <= gamma_a + gamma_b with constant 1.  Its
two quantitative ingredients are checked here (closing the coverage gap the
claim matrix declared on 2026-10-02):

  T1. Lemma lifted-state-T(b), the joint first moment.  gamma^(ab) is the first
      moment of the PRODUCT Fejer measure against the product geodesic
      distance sqrt(d(e,g)^2 + d(e,h)^2), computed here by an independent 2-D
      quadrature over the two conjugacy-class angles.  It must satisfy the
      strict sandwich
          max(gamma_a, gamma_b) < gamma^(ab) < gamma_a + gamma_b .
      This rejects (i) the rate written as max(gamma) (withdrawn 2026-10-02:
      the proved bound is never below max(gamma)), and (ii) an L1 / sum
      product distance (which would give equality with the sum).  The
      single-factor moments are cross-checked against geovac's gamma_rate
      (a second route, 1-D sympy-lambdified mpmath quadrature), and the (2,2)
      joint value is pinned by an adaptive 2-D mpmath quadrature, which
      rejects a helper using another product distance (the L-infinity max
      passes the sandwich but not the pin).

  T2. Lemma lifted-state-T(b), the product-state smoothing identity.  With the
      lifted state xi_ab = xi_a (x) xi_b (the Clifford C^2 factor is inert for
      the state and is omitted), translating by (g, h) and reading the
      expectation of a product multiplier M_a (x) M_b gives exactly
          sigma_J sigma_K D^J(g)_{A,B} D^K(h)_{C,D},
      i.e. the product of the two single-factor smoothings (Paper 38's T2
      identity, tests/test_p38_action_seminorm.py, on each factor).  This is
      what makes the induced measure the PRODUCT Fejer measure.  It rejects
      an entangled or mis-normalised joint state.
"""
from __future__ import annotations

from fractions import Fraction

import mpmath
import numpy as np
import pytest
import sympy as sp

from geovac.central_fejer_su2 import (
    _chi,
    central_fejer_kernel_su2,
    gamma_rate,
)


# ---------------------------------------------------------------------------
# T1: the joint first moment
# ---------------------------------------------------------------------------

def _kernel_numpy(n_max: int):
    return sp.lambdify(_chi, central_fejer_kernel_su2(n_max, _chi), modules="numpy")


def _gauss_legendre(a: float, b: float, n: int):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w


def _class_measure(n_max: int, nodes: int = 400):
    """Nodes chi in [0, 2pi] and weights K(chi) sin^2(chi/2)/pi dchi.
    d(e, g) = chi (dual-Coxeter normalisation, as in gamma_rate)."""
    chi, w = _gauss_legendre(0.0, 2 * np.pi, nodes)
    K = _kernel_numpy(n_max)
    with np.errstate(divide="ignore", invalid="ignore"):
        Kv = np.asarray(K(chi), dtype=float)
    return chi, w * Kv * np.sin(chi / 2) ** 2 / np.pi


def _moments(n_a: int, n_b: int):
    ca, wa = _class_measure(n_a)
    cb, wb = _class_measure(n_b)
    mass_a, mass_b = wa.sum(), wb.sum()
    g_a = float((wa * ca).sum())
    g_b = float((wb * cb).sum())
    dist = np.sqrt(ca[:, None] ** 2 + cb[None, :] ** 2)
    g_ab = float((wa[:, None] * wb[None, :] * dist).sum())
    return mass_a, mass_b, g_a, g_b, g_ab


@pytest.mark.parametrize("n_a,n_b", [(2, 2), (2, 3), (3, 4), (4, 4)])
def test_joint_moment_sandwich(n_a: int, n_b: int) -> None:
    mass_a, mass_b, g_a, g_b, g_ab = _moments(n_a, n_b)
    # unit-mass probability measures (Lemma L2-T(a))
    assert abs(mass_a - 1) < 1e-10 and abs(mass_b - 1) < 1e-10
    # single-factor moments agree with geovac's independent route
    assert abs(g_a - float(gamma_rate(n_a, prec=30))) < 1e-8
    assert abs(g_b - float(gamma_rate(n_b, prec=30))) < 1e-8
    # the strict sandwich: the rate is neither max(gamma) nor the sum
    assert g_ab > max(g_a, g_b) + 1e-3, (g_ab, g_a, g_b)
    assert g_ab < g_a + g_b - 1e-3, (g_ab, g_a, g_b)


def test_joint_moment_independent_route() -> None:
    """Pin gamma^(ab) at (2,2) by a second, independent route: adaptive 2-D
    mpmath quadrature of the product-Riemannian (L2) distance against the
    product kernel.  Rejects a helper that silently uses another product
    distance (e.g. the L-infinity max, which the sandwich alone cannot see)."""
    mpmath.mp.dps = 20
    K = sp.lambdify(_chi, central_fejer_kernel_su2(2, _chi), modules="mpmath")

    def w(c):
        return K(c) * mpmath.sin(c / 2) ** 2 / mpmath.pi

    val = mpmath.quad(lambda x, y: w(x) * w(y) * mpmath.sqrt(x ** 2 + y ** 2),
                      [0, 2 * mpmath.pi], [0, 2 * mpmath.pi])
    _, _, g_a, g_b, g_ab = _moments(2, 2)
    assert abs(float(val) - g_ab) < 1e-6, (float(val), g_ab)
    # and the tightness ratio the claim matrix quotes (~25% below the sum)
    assert 0.70 < g_ab / (g_a + g_b) < 0.78


# ---------------------------------------------------------------------------
# T2: the product-state smoothing identity (scalar window j <= 1/2 per factor)
# ---------------------------------------------------------------------------

_JS = [Fraction(0), Fraction(1, 2)]


def _jm(j):
    k = int(round(2 * j))
    return [Fraction(k, 2) - i for i in range(k + 1)]


def _index():
    index, idx = {}, 0
    for j in _JS:
        for a in _jm(j):
            for b in _jm(j):
                index[(j, a, b)] = idx
                idx += 1
    return index, idx


def _su2_generators(j):
    ms = [float(Fraction(int(round(2 * j)), 2)) - i for i in range(int(round(2 * j)) + 1)]
    d = len(ms)
    Jz = np.diag(ms).astype(complex)
    Jp = np.zeros((d, d), dtype=complex)
    for i in range(1, d):
        m = ms[i]
        Jp[i - 1, i] = np.sqrt(j * (j + 1) - m * (m + 1))
    Jy = (Jp - Jp.conj().T) / 2j
    return Jz, Jy


def _wigner_D(j, a, b, c):
    from scipy.linalg import expm
    Jz, Jy = _su2_generators(j)
    return expm(-1j * a * Jz) @ expm(-1j * b * Jy) @ expm(-1j * c * Jz)


def _multiplier(J, A, B, index, dim):
    from sympy import S
    from sympy.physics.quantum.cg import CG as _CG

    def cg(j1, m1, j2, m2, j3, m3):
        return float(_CG(S(j1), S(m1), S(j2), S(m2), S(j3), S(m3)).doit())

    M = np.zeros((dim, dim), dtype=complex)
    for j2 in _JS:
        for j1 in _JS:
            if abs(J - j2) > j1 or j1 > J + j2 or (j1 + j2 - J) % 1 != 0:
                continue
            pref = np.sqrt((int(2 * j2) + 1) / (int(2 * j1) + 1))
            for a2 in _jm(j2):
                a1 = A + a2
                if abs(a1) > j1:
                    continue
                ca = cg(J, A, j2, a2, j1, a1)
                if ca == 0.0:
                    continue
                for b2 in _jm(j2):
                    b1 = B + b2
                    if abs(b1) > j1:
                        continue
                    cb = cg(J, B, j2, b2, j1, b1)
                    if cb == 0.0:
                        continue
                    M[index[(j1, a1, b1)], index[(j2, a2, b2)]] = pref * ca * cb
    return M


def _U(g):
    from scipy.linalg import block_diag
    blocks = []
    for j in _JS:
        Dg = _wigner_D(float(j), *g)
        blocks.append(np.kron(Dg.conj(), np.eye(Dg.shape[0])))
    return block_diag(*blocks)


def _fejer_vector(index, dim):
    xi = np.zeros(dim, dtype=complex)
    for j in _JS:
        for a in _jm(j):
            xi[index[(j, a, a)]] = 1.0
    return xi / np.linalg.norm(xi)


def test_product_state_smoothing_identity() -> None:
    index, dim = _index()
    xi = _fejer_vector(index, dim)
    xi_ab = np.kron(xi, xi)                      # the lifted product state
    sigma = {Fraction(1, 2): np.sqrt(2) / 3, Fraction(1): 2.0 / 9.0}
    rng = np.random.default_rng(39)
    grid = [(tuple(rng.uniform(0, 2 * np.pi, 3)), tuple(rng.uniform(0, 2 * np.pi, 3)))
            for _ in range(6)]
    for J in sigma:
        for K in sigma:
            msJ, msK = _jm(J), _jm(K)
            MJ = _multiplier(J, msJ[0], msJ[-1], index, dim)
            MK = _multiplier(K, msK[0], msK[-1], index, dim)
            M = np.kron(MJ, MK)
            for g, h in grid:
                v = np.kron(_U(g), _U(h)) @ xi_ab
                ups = np.vdot(v, M @ v)
                expected = (sigma[J] * sigma[K]
                            * _wigner_D(float(J), *g)[0, -1]
                            * _wigner_D(float(K), *h)[0, -1])
                assert abs(ups - expected) < 1e-12, (float(J), float(K), g, h)
