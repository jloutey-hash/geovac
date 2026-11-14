"""Paper 58: the exchange class (AB|AB) at sigma != 0 (l > 0), FAST guard.

The only sigma != 0 value check before this file was
`test_exchange_sigma_nonzero_vs_md` (slow, Gaussian-fit referee at ~1e-5
tolerance). This file checks one sigma = 1 quartet,

    (2p_{+1,A} 1s_B | 1s_A 2p_{+1,B}),  Z_A = 3, Z_B = 1, R = 3,

against an independent reference at ~1e-10 (numerical xi path) and ~1e-12
(closed-form xi path), at the SAME Neumann truncation tau_max.

The reference shares nothing with geovac/two_center_eri.py: own hydrogenic
orbitals (scipy genlaguerre / lpmv, densities evaluated pointwise in prolate
coordinates), own Legendre functions (P via numpy Legendre derivatives, Q via
the hypergeometric form), and the textbook Neumann expansion of 1/r12, whose
transcription -- including the (-1)^sigma sign and the [(tau-s)!/(tau+s)!]^2
weight -- is itself checked pointwise against 1/|r1 - r2|
(`test_reference_neumann_kernel_reproduces_coulomb`). The xi double integral
is done by the ordered split t2 = t1 w with a double-exponential map
t = e^s, self-converged to ~1e-14.

Wrong answers this rejects (fire-tested with debug/qa/fire_test.py):
  * (-1)^sigma weight sign flipped
  * the tau = sigma term of the Neumann sum dropped / the last term dropped
  * the P_{k-1} P_{tau-k} sign in Q_tau flipped
  * [(tau-s)!/(tau+s)!]^2 taken to the first power
  * eta exponent sign flipped (QUARTET_ODD only; an exact symmetry of QUARTET)
"""
from __future__ import annotations

from fractions import Fraction
from math import factorial, gamma, pi, sqrt

import numpy as np
import pytest

special = pytest.importorskip("scipy.special")
mp = pytest.importorskip("mpmath")

from numpy.polynomial import legendre as Leg  # noqa: E402

from geovac.two_center_eri import exchange_value  # noqa: E402

R = 3.0
QUARTET = (3, (2, 1, 1), (1, 0, 0), 1, (1, 0, 0), (2, 1, 1))   # ZA, a, b, ZB, c, d
# Second quartet: 2s_B makes the density polynomial ODD in eta. QUARTET's
# polynomials are even in eta, so an eta -> -eta error (exponent sign in the
# eta half) is an exact symmetry there and invisible; this one sees it.
QUARTET_ODD = (3, (2, 1, 1), (2, 0, 0), 1, (1, 0, 0), (2, 1, 1))
TAU_MAX = 4


# ---------------------------------------------------------------- reference
def _radial(Z, n, l, r):
    rho = 2.0 * Z * r / n
    N = sqrt((2.0 * Z / n) ** 3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return N * np.exp(-rho / 2) * rho ** l * special.eval_genlaguerre(n - l - 1, 2 * l + 1, rho)


def _orb(Z, o, r, c):
    """chi with its e^{i m phi} stripped (the phase convention cancels: each
    2p_{+1} enters once conjugated and once not)."""
    n, l, m = o
    N = sqrt((2 * l + 1) / (4 * pi) * factorial(l - abs(m)) / factorial(l + abs(m)))
    return _radial(Z, n, l, r) * N * special.lpmv(abs(m), l, np.clip(c, -1.0, 1.0))


def _density(Z1, o1, Z2, o2, xi, eta):
    rA, rB = R * (xi + eta) / 2, R * (xi - eta) / 2
    return (_orb(Z1, o1, rA, (1 + xi * eta) / (xi + eta))
            * _orb(Z2, o2, rB, (xi * eta - 1) / (xi - eta)))


def _P_out(tau, s, x):
    """(x^2 - 1)^{s/2} d^s P_tau / dx^s on x >= 1."""
    c = np.zeros(tau + 1)
    c[tau] = 1.0
    return (x * x - 1) ** (s / 2) * Leg.legval(x, Leg.legder(c, s))


def _Q(tau, x):
    """Q_tau(x), x > 1, hypergeometric form (no P Q_0 - W cancellation)."""
    return (sqrt(pi) * gamma(tau + 1) / (gamma(tau + 1.5) * (2 * x) ** (tau + 1))
            * special.hyp2f1((tau + 1) / 2, (tau + 2) / 2, tau + 1.5, 1 / (x * x)))


def _Q_out(tau, s, x):
    """(x^2 - 1)^{s/2} d^s Q_tau / dx^s, s = 0, 1 ((x^2-1) Q' = tau (x Q_tau - Q_{tau-1}))."""
    if s == 0:
        return _Q(tau, x)
    assert s == 1
    return tau * (x * _Q(tau, x) - _Q(tau - 1, x)) / np.sqrt(x * x - 1)


def _neumann_weight(tau, s):
    return (-1) ** s * (2 * tau + 1) * (factorial(tau - s) / factorial(tau + s)) ** 2


def reference_exchange(quartet, tau_max, h=0.12, s_lo=-24.0, s_hi=3.2, n_w=48, n_e=48):
    ZA, oa, ob, ZB, oc, od = quartet
    s = oa[2] - ob[2]
    assert s == od[2] - oc[2] == 1
    t1 = np.exp(np.arange(s_lo, s_hi + 1e-12, h))      # xi - 1, double-exponential map
    wn, ww = Leg.leggauss(n_w)
    wn, ww = (wn + 1) / 2, ww / 2
    en, ew = Leg.leggauss(n_e)
    tw = t1[:, None] * wn[None, :]                       # t2 = t1 * w, w in (0, 1)

    def F(t, Z1, o1, Z2, o2):
        xi = 1.0 + t[..., None]
        return _density(Z1, o1, Z2, o2, xi, en) * (xi * xi - en * en)

    F1a, F1b = F(t1, ZA, oa, ZB, ob), F(tw, ZA, oa, ZB, ob)
    F2a, F2b = F(t1, ZA, oc, ZB, od), F(tw, ZA, oc, ZB, od)
    total = 0.0
    for tau in range(s, tau_max + 1):
        pe = ew * special.lpmv(s, tau, en)
        G1a, G1b, G2a, G2b = F1a @ pe, F1b @ pe, F2a @ pe, F2b @ pe
        Pb = _P_out(tau, s, 1 + tw)
        Qa = _Q_out(tau, s, 1 + t1)
        # xi_< = t1 w, xi_> = t1, both orderings; dt1 = t1 ds, dt2 = t1 dw
        integrand = t1 * t1 * Qa * (G1a * ((G2b * Pb) @ ww) + G2a * ((G1b * Pb) @ ww))
        total += _neumann_weight(tau, s) * h * integrand.sum()
    return (R / 2) ** 6 * (2 * pi) ** 2 * (2 / R) * total


_REF: dict = {}


def _ref(quartet):
    if quartet not in _REF:
        _REF[quartet] = reference_exchange(quartet, TAU_MAX)
    return _REF[quartet]


@pytest.fixture(scope="module")
def ref():
    return _ref(QUARTET)


# ---------------------------------------------------------------- the referee's own checks
def test_reference_neumann_kernel_reproduces_coulomb():
    """The Neumann transcription used above (sign, weight, P/Q conventions) sums
    to 1/|r1 - r2| pointwise, and the fast Q^1 matches mpmath's type-3 Q."""
    def cart(xi, eta, phi):
        rho = R / 2 * sqrt((xi * xi - 1) * (1 - eta * eta))
        return np.array([rho * np.cos(phi), rho * np.sin(phi), R / 2 * (1 + xi * eta)])

    for p1, p2 in (((1.3, 0.4, 0.2), (2.7, -0.5, 1.9)), ((1.9, -0.8, 2.5), (1.15, 0.6, 0.4))):
        exact = 1 / np.linalg.norm(cart(*p1) - cart(*p2))
        (x1, e1, f1), (x2, e2, f2) = p1, p2
        xl, xg = min(x1, x2), max(x1, x2)
        tot = 0.0
        for t in range(0, 41):
            for s in range(0, t + 1):
                tot += ((1 if s == 0 else 2) * _neumann_weight(t, s)
                        * float(mp.legenp(t, s, xl, type=3).real)
                        * float(mp.legenq(t, s, xg, type=3).real)
                        * special.lpmv(s, t, e1) * special.lpmv(s, t, e2)
                        * np.cos(s * (f1 - f2)))
        assert abs(2 / R * tot - exact) < 1e-7 * exact, (p1, p2, 2 / R * tot, exact)
    for t in (1, 3, 4):
        for x in (1.001, 1.3, 2.5, 8.0):
            assert abs(_Q_out(t, 1, x) - float(mp.legenq(t, 1, x, type=3).real)) \
                <= 1e-12 * abs(float(mp.legenq(t, 1, x, type=3).real))
            assert abs(_P_out(t, 1, x) - float(mp.legenp(t, 1, x, type=3).real)) \
                <= 1e-12 * abs(float(mp.legenp(t, 1, x, type=3).real))


def test_reference_self_converged(ref):
    coarse = reference_exchange(QUARTET, TAU_MAX, h=0.18, n_w=36, n_e=36)
    assert abs(coarse - ref) <= 1e-12 * abs(ref), (coarse, ref)


# ---------------------------------------------------------------- the guards
@pytest.mark.parametrize("quartet", [QUARTET, QUARTET_ODD], ids=["even-eta", "odd-eta"])
def test_exchange_sigma1_numerical_xi_vs_reference(quartet):
    """Production path (eta closed form, xi nested quadrature): ~1e-10."""
    ref = _ref(quartet)
    ZA, oa, ob, ZB, oc, od = quartet
    got = exchange_value(Fraction(ZA), oa, ob, Fraction(ZB), oc, od, R, tau_max=TAU_MAX)
    assert abs(got - ref) <= 3e-10 * abs(ref), f"{got!r} vs reference {ref!r}"


def test_exchange_sigma1_exact_xi_vs_reference(ref):
    """Closed-form xi path (ordered_xi_general at sigma = 1): ~1e-12."""
    ZA, oa, ob, ZB, oc, od = QUARTET
    got = exchange_value(Fraction(ZA), oa, ob, Fraction(ZB), oc, od, R, tau_max=TAU_MAX,
                         exact_xi=True)
    assert abs(got - ref) <= 1e-12 * abs(ref), f"{got!r} vs reference {ref!r}"
