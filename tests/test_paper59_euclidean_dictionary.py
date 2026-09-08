"""Papers 59 + 61 -- the 2D-Euclidean-QFT kinematic dictionary.

sec:reduction is Paper 59; sec:modular moved to Paper 61 in the 2026-09-06
split, so this docstring named one paper for two owners (corrected 2026-09-07).

Backs the kinematics paragraphs (2026-08-21, debug/sprint_t2_euclidean_kinematics_memo.md):
  (1) eq:K0 is the frequency representation of the 2D Euclidean propagator
      G_a(tau,x) = (1/2pi) K_0(a sqrt(tau^2+x^2)):  LHS = (2 pi/sqrt c) G_a(p,b),
      a = zeta/sqrt(c), p = D sqrt(c).
  (2) eq:period with the convention-free prefactor:  int dk/sqrt((c1 k^2+1)(c2 k^2+1))
      = K(m)/sqrt(c_max), m = 1 - c_min/c_max  (BOTH orderings of c1,c2).
  (3) The fibre is a two-propagator correlator:
      int_0^inf dk e^{-D1 Dl1 - D2 Dl2}/(Dl1 Dl2) = 4 pi a1 a2 int_R db G_{a1}(p1,b) G_{a2}(p2,b).
  (4) Family invariance a_i p_i = zeta D_i; the Stokes location is the squared complexified
      Euclidean interval |z*| = p1^2 + |W|^2, of which -(sqrt(c1) -+ i b)^2 is the D=1 case.
  (5) Oscillatory-corner closed form:  J(s,t) = sig^4 (pi/2)(7/e)^2 alpha(1-alpha)/b_W + O(sig^5)
      at the (1,0) corner  (J = int j0(k(s+t)) P(s,k) P(t,k) dk, P = c e^{-Dl}(Dl^-3+3Dl^-4+3Dl^-5)).
Heavy multi-point validation lives in debug/t2_euclidean_dictionary.py (sections A-F)."""
import sympy as sp
from mpmath import mp, mpf, mpc, sqrt, cos, exp, besselk, ellipk, pi, e, quad, quadosc, sin, inf


def _G2D(tau, x, m):
    return besselk(0, m * sqrt(tau**2 + x**2)) / (2 * pi)


def test_eqK0_is_2d_propagator():
    mp.dps = 25
    c, zeta, D, b = mpf('0.17'), mpf('1.0'), mpf('0.8'), mpf('0.6')
    lhs = quad(lambda k: cos(k * b) * exp(-D * sqrt(c * k**2 + zeta**2))
               / sqrt(c * k**2 + zeta**2), [0, inf])
    a, p = zeta / sqrt(c), D * sqrt(c)
    rhs = (2 * pi / sqrt(c)) * _G2D(p, b, a)
    assert abs(lhs - rhs) / abs(rhs) < mpf(10)**(-18)


def test_eq_period_prefactor_cmax_both_orderings():
    mp.dps = 25
    for c1, c2 in [(mpf('0.21'), mpf('0.04')), (mpf('0.04'), mpf('0.21'))]:
        lhs = quad(lambda k: 1 / sqrt((c1 * k**2 + 1) * (c2 * k**2 + 1)), [0, inf])
        m = 1 - min(c1, c2) / max(c1, c2)
        rhs = ellipk(m) / sqrt(max(c1, c2))
        assert abs(lhs - rhs) / rhs < mpf(10)**(-18)


def test_two_propagator_correlator():
    mp.dps = 20
    c1, c2, D1, D2 = mpf('0.12'), mpf('0.05'), mpf('1.0'), mpf('1.0')
    Dl = lambda c, k: sqrt(c * k**2 + 1)
    lhs = quad(lambda k: exp(-D1 * Dl(c1, k) - D2 * Dl(c2, k)) / (Dl(c1, k) * Dl(c2, k)),
               [0, inf])
    a1, a2 = 1 / sqrt(c1), 1 / sqrt(c2)
    p1, p2 = D1 * sqrt(c1), D2 * sqrt(c2)
    rhs = 4 * pi * a1 * a2 * 2 * quad(
        lambda b: _G2D(p1, b, a1) * _G2D(p2, b, a2), [0, inf])
    assert abs(lhs - rhs) / abs(rhs) < mpf(10)**(-14)


def test_family_invariance_and_zstar_specialisation():
    zeta, c, D, W = sp.symbols('zeta c D W', positive=True)
    a, p = zeta / sp.sqrt(c), D * sp.sqrt(c)
    assert sp.simplify(a * p - zeta * D) == 0          # mass x time, Feynman-parameter-free
    zstar = -(D * sp.sqrt(c) - sp.I * W)**2            # interval form
    assert sp.simplify(zstar.subs(D, 1) - (-(sp.sqrt(c) - sp.I * W)**2)) == 0   # D=1 case
    assert sp.simplify(sp.Abs(zstar.rewrite(sp.re)) - (p**2 + W**2)) == 0 or \
        sp.simplify(sp.expand(zstar * sp.conjugate(zstar)) - (p**2 + W**2)**2) == 0


def _P(x, k):
    c = x * (1 - x)
    D = sqrt(c * k * k + 1)
    return c * exp(-D) * (D**-3 + 3 * D**-4 + 3 * D**-5)


def test_oscillatory_corner_closed_form():
    mp.dps = 15
    alpha, sig = mpf('0.5'), mpf('0.04')
    s, t = 1 - sig**2 * alpha, sig**2 * (1 - alpha)      # the (1,0) corner
    bW = s + t
    J = quadosc(lambda k: (sin(k * bW) / (k * bW)) * _P(s, k) * _P(t, k),
                [0, inf], period=pi / bW)
    pred = (pi / 2) * (7 / e)**2 * alpha * (1 - alpha) / bW
    assert abs(J / sig**4 - pred) / pred < mpf('0.02')   # O(sig^2) approach; driver: 1e-4 at small sig
