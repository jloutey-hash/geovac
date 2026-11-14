"""Paper 58 sec:qfd (i): the shell-route hybrid (AA|AB) evaluator covers EVERY
rate ordering a_c = Z_A/n_c vs a_d = Z_B/n_d.

Claim guarded (paper_58_abelian_residue.tex, "shell-route hybrid evaluator
covers every ordering"):
  * a_d = a_c (rate coincidence): the closed form exists, is finite, and the
    E_1 of the generic case collapses to the pair gamma_E + ln R (equal
    coefficients);
  * a_d > a_c: the closed form carries Ei(|a_d - a_c| R) and is REAL;
  * both match an independent quadrature;
  * the coincidence value is the limit of its neighbours (continuity).

The reference below is self-contained: its own hydrogenic orbitals, its own
one-centre potential (textbook r_<^L / r_>^{L+1} expansion, radial integrals by
sympy), and Gauss-Legendre quadrature in (r_A, r_B) triangle coordinates. It
shares nothing with geovac/two_center_eri.py -- not the Gaunt re-coupling, not
the radial polynomials, not V_L. Self-convergence of the reference: ~1e-14.

Wrong answers this rejects (fire-tested with debug/qa/fire_test.py when written):
  1. the j = 0 Ein term dropped from log_moment_finite      -> eq value moves
  2. E_1 instead of -Ei for a negative rate in e1_moment      -> closed form goes
     complex (Im = O(pi * coeff)); its REAL part is unchanged, so a guard that
     takes re() -- as the debug driver does -- cannot see this one
  3. the b = 0 branch of finite_power_exp removed (the old 0/0 path) -> no value
"""
from __future__ import annotations

from fractions import Fraction
from math import factorial, pi, sqrt

import numpy as np
import pytest

sp = pytest.importorskip("sympy")
special = pytest.importorskip("scipy.special")

from geovac.two_center_eri import R_s, hybrid_closed_form  # noqa: E402

_r = sp.Symbol("r", positive=True)
_rp = sp.Symbol("rp", positive=True)


# ---------------------------------------------------------------- reference
def _R_sym(Z, n, l):
    Z = sp.Rational(Z)
    rho = 2 * Z * _r / n
    N = sp.sqrt((2 * Z / n) ** 3 * sp.factorial(n - l - 1) / (2 * n * sp.factorial(n + l)))
    return N * sp.exp(-rho / 2) * rho ** l * sp.assoc_laguerre(n - l - 1, 2 * l + 1, rho)


def _ynorm(l):
    return sqrt((2 * l + 1) / (4 * pi))


def _potential(Z, oa, ob):
    """V(r, theta) of chi_a chi_b (m = 0, both on A) = sum_L v_L(r) P_L(cos)."""
    (na, la, _), (nb, lb, _) = oa, ob
    prod = (np.polynomial.legendre.Legendre.basis(la)
            * np.polynomial.legendre.Legendre.basis(lb)).coef
    radial = sp.expand(_R_sym(Z, na, la) * _R_sym(Z, nb, lb))
    out = {}
    for L, cL in enumerate(prod):
        if abs(cL) < 1e-15:
            continue
        f = sp.nsimplify(cL) * radial.subs(_r, _rp)
        inner = sp.integrate(f * _rp ** (L + 2), (_rp, 0, _r))
        outer = sp.integrate(f * _rp ** (1 - L), (_rp, _r, sp.oo))
        vL = 4 * sp.pi / (2 * L + 1) * (inner / _r ** (L + 1) + _r ** L * outer)
        out[L] = sp.lambdify(_r, vL, "numpy")
    return out, _ynorm(la) * _ynorm(lb)


def _orbital(Z, o):
    n, l, m = o
    assert m == 0
    f = sp.lambdify(_r, _R_sym(Z, n, l), "numpy")
    return lambda rr, c: f(rr) * _ynorm(l) * special.eval_legendre(l, c)


def hybrid_reference(ZA, oa, ob, oc, ZB, od, R, n=60, tail=40.0):
    """(ab|cd), a, b, c on A (origin), d on B at +R z, real m = 0 orbitals."""
    V, nab = _potential(ZA, oa, ob)
    chic, chid = _orbital(ZA, oc), _orbital(ZB, od)
    x, w = np.polynomial.legendre.leggauss(n)

    def inner(rb):
        lo, hi = np.abs(rb - R), rb + R
        ra = lo[:, None] + (hi - lo)[:, None] * (x[None, :] + 1) / 2
        wa = (hi - lo)[:, None] * w[None, :] / 2
        cA = np.clip((ra * ra + R * R - rb[:, None] ** 2) / (2 * R * ra), -1, 1)
        cB = np.clip((ra * ra - rb[:, None] ** 2 - R * R) / (2 * R * rb[:, None]), -1, 1)
        Vv = nab * sum(v(ra) * special.eval_legendre(L, cA) for L, v in V.items())
        f = ra * rb[:, None] / R * Vv * chic(ra, cA) * chid(rb[:, None], cB)
        return (f * wa).sum(axis=1)

    total = 0.0
    edges = [0.0, R] + list(R + np.array([1.0, 3.0, 7.0, 15.0, tail]))
    for a, b in zip(edges[:-1], edges[1:]):
        rb = a + (b - a) * (x + 1) / 2
        total += ((b - a) / 2 * w * inner(rb)).sum()
    return 2 * pi * total


# ---------------------------------------------------------------- cases
# (2p0 2p0 | 1s_A 1s_B): a_c = Z_A/1, a_d = Z_B/1.
_PAIR = ((2, 1, 0), (2, 1, 0), (1, 0, 0))
EQ = (3, 3)        # a_d = a_c = 3     (rate coincidence)
GT = (1, 2)        # a_d = 2 > a_c = 1 (Ei branch)
R_VALUES = (sp.Rational(5, 2), sp.Rational(7, 5))


@pytest.fixture(scope="module")
def closed():
    cache = {}

    def get(ZA, ZB):
        key = (Fraction(ZA), Fraction(ZB))
        if key not in cache:
            cache[key] = hybrid_closed_form(key[0], *_PAIR, key[1], (1, 0, 0))
        return cache[key]
    return get


def _value(e, Rv):
    v = complex(sp.N(e.subs(R_s, Rv), 30))
    return v


@pytest.mark.parametrize("Zs", [EQ, GT], ids=["a_d=a_c", "a_d>a_c"])
@pytest.mark.parametrize("Rv", R_VALUES, ids=["R=5/2", "R=7/5"])
def test_value_matches_independent_quadrature(closed, Zs, Rv):
    """Closed form is REAL and equals the self-contained quadrature to 1e-10."""
    ZA, ZB = Zs
    v = _value(closed(ZA, ZB), Rv)
    ref = hybrid_reference(ZA, *_PAIR, ZB, (1, 0, 0), float(Rv))
    assert np.isfinite(v.real) and np.isfinite(v.imag), v
    # A complex value whose real part is right is the signature of E_1 taken on
    # a negative argument (E_1(-x) = -Ei(x) - i pi): reject it explicitly.
    assert abs(v.imag) <= 1e-20 + 1e-14 * abs(v.real), f"not real: {v}"
    assert abs(v.real - ref) <= 1e-10 * abs(ref), f"{v.real!r} vs quadrature {ref!r}"


def test_function_content_per_branch(closed):
    """a_d > a_c carries Ei(|a_d - a_c| R) with positive argument and no E_1 of a
    negative argument; the coincidence carries no Ei at all."""
    e_gt, e_eq = closed(*GT), closed(*EQ)
    ei_args = {f.args[0] for f in e_gt.atoms(sp.Ei)}
    assert ei_args == {sp.Integer(GT[1] - GT[0]) * R_s}, ei_args
    for f in e_gt.atoms(sp.expint) | e_eq.atoms(sp.expint):
        arg = f.args[1].subs(R_s, 1)
        assert arg.is_positive, f"E_1 of a non-positive argument: {f}"
    assert not e_eq.atoms(sp.Ei), e_eq.atoms(sp.Ei)


def test_coincidence_gamma_plus_lnR_pair(closed):
    """At a_d = a_c the E_1 collapses to gamma_E + ln R: both present, with
    IDENTICAL coefficients, and neither survives anywhere else in the form."""
    e = sp.expand(sp.expand_log(closed(*EQ), force=True))
    g = e.coeff(sp.EulerGamma)
    L = e.coeff(sp.log(R_s))
    rest = sp.expand(e - g * sp.EulerGamma - L * sp.log(R_s))
    assert not rest.has(sp.EulerGamma) and not rest.has(sp.log(R_s))
    assert g != 0, "no Euler gamma at the coincidence"
    for Rv in (sp.Rational(5, 2), sp.Rational(7, 5), sp.Integer(4)):
        gv, Lv = sp.N(g.subs(R_s, Rv), 30), sp.N(L.subs(R_s, Rv), 30)
        assert abs(gv) > 1e-6, f"gamma coefficient vanishes at R={Rv}"
        assert abs(gv - Lv) <= 1e-25 * abs(gv), f"R={Rv}: gamma {gv} vs ln R {Lv}"


def test_continuity_across_coincidence(closed):
    """v(a_d = a_c (1 +- eps)) -> v(a_d = a_c) linearly: the exact log branch is
    the limit of its neighbours on BOTH sides (the a_d > a_c side is the Ei
    branch, the a_d < a_c side the pre-existing E_1 branch)."""
    Rv = sp.Rational(5, 2)
    v0 = _value(closed(*EQ), Rv).real
    eps = Fraction(1, 10 ** 6)
    d = {}
    for sgn in (+1, -1):
        e = hybrid_closed_form(Fraction(3), *_PAIR, Fraction(3) * (1 + sgn * eps),
                               (1, 0, 0))
        d[sgn] = _value(e, Rv).real - v0
    slope = (d[+1] - d[-1]) / (2 * float(eps))
    assert 0 < abs(d[+1]) < 1e-6 and 0 < abs(d[-1]) < 1e-6, d
    assert d[+1] * d[-1] < 0, f"not a crossing: {d}"
    assert abs(d[+1] + d[-1]) <= 1e-3 * abs(d[+1]), f"asymmetric (kink/offset): {d}"
    assert abs(slope) < 1.0, slope
