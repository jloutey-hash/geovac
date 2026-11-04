"""Paper 59 sec:reduction -- owed guards (group2 /qa FULL 2026-11-03/04, docs/qa/group2.debt.md).

Three claims of Paper 59 sec:reduction, each held to the number the prose prints:

(1) eq:pf  [SYMBOLIC PROOF].  With f = Q^{-1/2}, Q = (x^2-1)(rho x^2 + 1 - rho), and
    L(D) = int_1^inf e^{-Dx} f dx, so that L^{(n)}(D) = (-1)^n int_1^inf x^n e^{-Dx} f dx,
    applying the operator of eq:pf to L under the integral gives int_1^inf e^{-Dx} f P(x,D) dx
    with P built from the five printed coefficients.  The derivation in the prose ("integration
    by parts ... yields an exact fourth-order Picard--Fuchs equation") is the statement that
        e^{-Dx} f(x) P(x, D)  ==  -d/dx [ sqrt(Q) e^{-Dx} ]           (exact, as functions)
    and that the boundary term sqrt(Q) e^{-Dx} vanishes at x = 1 (branch point) and x -> inf.
    The existing `test_paper59_L4_annihilates_the_physical_master_N` is a 1e-14 NUMERICAL
    residual; this is the exact identity.  The operator coefficients are typed here from the
    paper (eq:pf, Paper 59 ~L402), not taken from code.

(2) The kinematic dictionary at FIVE parameter points each [MEASURED, 1e-30 and 1e-31]:
    D1  eq:K0 = (2 pi / sqrt c) G_a(p, b),  a = m/sqrt c, p = D sqrt c;
    D2  eq:besselmoment = 4 pi a1 a2 int_R db G_{a1}(p1,b) G_{a2}(p2,b).
    `test_paper59_euclidean_dictionary.py` checks one point each at 1e-18 / 1e-14; the prose
    number is "five parameter points each".  Slow (~15 s at 32 digits).

(3) The full position-space fibre [MEASURED, 1.8e-10; "not regression-tested" until now]:
        J(s,t) = int_0^inf dk j0(k|W|) P(s,k) P(t,k)
               = (2/pi) int_0^inf db f_s(b) <f_t>_{|W|}(b),
    f_x(b) = sqrt(c) int_1^inf w(u) u K_1(R)/R du (the propagator form), <.>_{|W|} the uniform
    box average over [b-|W|, b+|W|], |W| = s + t.  Computed here by an independent float64
    tensor-Gauss route (numpy/scipy) against an mpmath momentum-space quadrature.  Slow.
"""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp
from mpmath import mp, mpf, sqrt as msqrt, cos as mcos, exp as mexp, besselk, pi as mpi, quad, \
    quadosc, inf, sin as msin
from scipy.special import k1 as K1


# ---------------------------------------------------------------- (1) eq:pf exact identity
# eq:pf, typed from the paper:
#   D rho L'''' + 2 rho L''' + D(1-2rho) L'' + (1-2rho) L' - D(1-rho) L = 0
def _pf_coeffs(D, rho):
    """[c4, c3, c2, c1, c0] multiplying L'''' .. L in eq:pf."""
    return [D * rho, 2 * rho, D * (1 - 2 * rho), (1 - 2 * rho), -D * (1 - rho)]


def test_paper59_eqpf_integrand_is_exact_derivative():
    x, D, rho = sp.symbols('x D rho', positive=True)
    Q = (x**2 - 1) * (rho * x**2 + 1 - rho)
    f = 1 / sp.sqrt(Q)
    c = _pf_coeffs(D, rho)
    # L^{(n)} -> (-1)^n x^n under the integral
    P = sum(c[4 - n] * (-1)**n * x**n for n in range(5))
    integrand = sp.exp(-D * x) * f * P
    exact = -sp.diff(sp.sqrt(Q) * sp.exp(-D * x), x)
    residual = sp.simplify((integrand - exact) * sp.sqrt(Q) * sp.exp(D * x))
    assert residual == 0, residual
    # the boundary term sqrt(Q) e^{-Dx} vanishes at the branch point and at infinity
    bt = sp.sqrt(Q) * sp.exp(-D * x)
    assert sp.simplify(bt.subs(x, 1)) == 0
    assert sp.limit(bt.subs({rho: sp.Rational(1, 3), D: sp.Rational(7, 10)}), x, sp.oo) == 0


def test_paper59_eqpf_characteristic_polynomial():
    """The prose factorisation rho(s^2-1)(s^2+(1-rho)/rho) is the D-leading part of eq:pf."""
    s, D, rho = sp.symbols('s D rho', positive=True)
    c = _pf_coeffs(D, rho)
    lead = sp.expand(sum(sp.Poly(c[4 - n], D).coeff_monomial(D) * s**n for n in range(5)))
    assert sp.expand(lead - rho * (s**2 - 1) * (s**2 + (1 - rho) / rho)) == 0


# ---------------------------------------------------------------- (2) dictionary, five points
def _G2D(tau, x, m):
    return besselk(0, m * msqrt(tau**2 + x**2)) / (2 * mpi)


def _dampquad(f, A):
    return quad(f, [mpf(0)] + [mpf(q) / A for q in (0.5, 1, 2, 4, 8, 16, 32)] + [inf])


@pytest.mark.slow
def test_paper59_eqK0_propagator_five_points():
    mp.dps = 32
    pts = [('0.2', '1', '1', '0.7'), ('0.05', '1', '1', '1.5'), ('0.25', '1', '2', '0.3'),
           ('0.13', '0.8', '0.6', '2.2'), ('0.24', '1', '1', '0.05')]
    worst = mpf(0)
    for cs, ms, Ds, bs in pts:
        c, m, D, b = (mpf(z) for z in (cs, ms, Ds, bs))
        A = D * msqrt(c)
        f = lambda k: mcos(k * b) * mexp(-D * msqrt(c * k * k + m * m)) / msqrt(c * k * k + m * m)
        lhs = quadosc(f, [0, inf], period=2 * mpi / b) if b / A > 3 else _dampquad(f, A)
        rhs = (2 * mpi / msqrt(c)) * _G2D(D * msqrt(c), b, m / msqrt(c))
        worst = max(worst, abs(lhs - rhs) / abs(rhs))
    assert worst < mpf(10)**(-28), worst


@pytest.mark.slow
def test_paper59_two_propagator_correlator_five_points():
    mp.dps = 32
    pts = [('0.2', '0.1', '1', '1'), ('0.25', '0.05', '1', '1'), ('0.15', '0.15', '1', '1'),
           ('0.24', '0.01', '1', '2'), ('0.1', '0.2', '0.5', '1.5')]
    worst = mpf(0)
    for c1s, c2s, D1s, D2s in pts:
        c1, c2, D1, D2 = (mpf(z) for z in (c1s, c2s, D1s, D2s))
        a1, a2 = 1 / msqrt(c1), 1 / msqrt(c2)
        p1, p2 = D1 * msqrt(c1), D2 * msqrt(c2)
        f = lambda k: (mexp(-D1 * msqrt(c1 * k * k + 1)) / msqrt(c1 * k * k + 1)
                       * mexp(-D2 * msqrt(c2 * k * k + 1)) / msqrt(c2 * k * k + 1))
        lhs = _dampquad(f, p1 + p2)
        seg = sorted({mpf(0), min(p1, p2), max(p1, p2), mpf(1), mpf(3), mpf(10)})
        rhs = 4 * mpi * a1 * a2 * 2 * quad(lambda b: _G2D(p1, b, a1) * _G2D(p2, b, a2),
                                           list(seg) + [inf])
        worst = max(worst, abs(lhs - rhs) / abs(rhs))
    assert worst < mpf(10)**(-28), worst


# ---------------------------------------------------------------- (3) position-space fibre
def _w(u):
    """Laplace weight of the four zeta-derivatives: P/c = int_1^inf w(u) e^{-u Delta} du."""
    v = u - 1
    return v**2 / 2 + v**3 / 2 + v**4 / 8


def _P_mp(x, k):
    c = x * (1 - x)
    Dl = msqrt(c * k * k + 1)
    return c * mexp(-Dl) * (Dl**-3 + 3 * Dl**-4 + 3 * Dl**-5)


def _J_momentum(s, t):
    """J(s,t) = int_0^inf j0(k|W|) P(s,k) P(t,k) dk at 25 digits (mpmath)."""
    mp.dps = 25
    s, t = mpf(s), mpf(t)
    b = s + t
    A = msqrt(s * (1 - s)) + msqrt(t * (1 - t))
    f = lambda k: (msin(k * b) / (k * b) if k * b > mpf('1e-30') else mpf(1)) * _P_mp(s, k) * _P_mp(t, k)
    return float(_dampquad(f, A))


def _gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w


def _panels(edges, n):
    xs, ws = zip(*(_gl(n, a, b) for a, b in zip(edges[:-1], edges[1:])))
    return np.concatenate(xs), np.concatenate(ws)


_U, _WU = _panels([1, 2, 4, 8, 14, 25, 40, 60], 48)


def _f_prop(x, b):
    """Propagator (position-space) form f_x(b) = sqrt(c) int_1^inf w(u) u K_1(R)/R du, vectorised in b."""
    c = x * (1 - x)
    b = np.abs(np.asarray(b, dtype=float))[..., None]
    R = np.sqrt(_U**2 + b * b / c)
    return np.sqrt(c) * np.sum(_WU * _w(_U) * _U * K1(R) / R, axis=-1)


def _J_position(s, t):
    bW = s + t
    c1 = s * (1 - s)
    Bmax = 40 * np.sqrt(c1) + 2 * bW
    edges = np.unique(np.concatenate([np.linspace(0, bW, 5), np.linspace(bW, Bmax, 13)]))
    bs, wb = _panels(list(edges), 40)
    vs, wv = _panels(list(np.linspace(-bW, bW, 9)), 40)
    smeared = (_f_prop(t, bs[:, None] + vs[None, :]) * wv[None, :]).sum(axis=1) / (2 * bW)
    return (2 / np.pi) * np.sum(wb * _f_prop(s, bs) * smeared)


@pytest.mark.slow
@pytest.mark.parametrize("s,t", [(0.3, 0.2), (0.4, 0.4), (0.25, 0.6)])
def test_paper59_position_space_fibre(s, t):
    Jm = _J_momentum(s, t)
    Jp = _J_position(s, t)
    rel = abs(Jm - Jp) / abs(Jm)
    # Paper 59 prints 1.8e-10, the floor of its nested mpmath Gauss run; this float64
    # tensor-Gauss route agrees to ~1e-15 (measured 2026-11-04: 1.3e-15, 9.1e-16, 8.9e-16), so the
    # guard is set at 1e-12 -- tighter than the prose, never looser.
    assert rel < 1e-12, (s, t, Jm, Jp, rel)
