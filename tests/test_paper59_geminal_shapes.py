"""Paper 59 sec:density_ring: geminal shapes on the one-function reading.

Record docs/predictions/2026-10-22_geminal_shapes.md (GA1-GA2, committed 5e40da9a before computing); GA2 and the
unregistered slot diagnostic are guarded here.  Written as a separate step; no debug/ import.

  (a) GA2: L_G = t^4 (theta+1)^2 - t^2 theta^2 (4 theta + 1) + 4 theta^2 (theta - 2)^2 annihilates
      H_G = int K0(tx) K0(x) e^{-x^2} dx, and its monodromy companion J = int I0(tx) K0(x) e^{-x^2} dx grows like
      e^{(t-1)^2/4}: no Fuchsian operator admits that, so every annihilator is irregular at infinity.
  (b) The momentum slot is the geminal's: at gamma = 3/2 the two other role assignments give families whose
      traces of Frobenius differ from C2b's at some good prime, while J1 against its partner J2 and Legendre
      against its 2-isogenous form agree at every good prime.

Wrong answers rejected: a wrong coefficient in L_G; H_G itself offered as the growing companion; the piece's own
family offered as an alternative assignment (the slot as a convention).
"""
from __future__ import annotations

import mpmath as mp
import pytest
import sympy as sp

t = sp.Symbol("t", positive=True)
x = sp.Symbol("x")
LG_INNER = 1          # 4 theta + LG_INNER
COMPANION = "I0"      # the growing companion's Bessel factor
ALT_MOMENTUM = ("1", "t")   # which inverse length sits in the momentum slot in the alternatives


def _theta_poly_to_d(ops):
    """sum_k t^(2k) p_k(theta) -> d/dt coefficient list"""
    f = sp.Function("f")(t)
    expr = 0
    for k, p in ops.items():
        pw = [f]
        for _ in range(sp.Poly(p, x).degree()):
            pw.append(sp.expand(t * sp.diff(pw[-1], t)))
        expr += t ** (2 * k) * sum(sp.Poly(p, x).coeff_monomial(x**j) * pw[j] for j in range(len(pw)))
    expr = sp.expand(expr)
    cs = [sp.expand(expr.coeff(sp.Derivative(f, (t, k)) if k else f)) for k in range(5)]
    return cs


def _LG():
    u = x
    return _theta_poly_to_d({2: (u + 1) ** 2, 1: -u**2 * (4 * u + LG_INNER), 0: 4 * u**2 * (u - 2) ** 2})


def _k0_derivs(z):
    k0, k1 = mp.besselk(0, z), mp.besselk(1, z)
    return [k0, -k1, k0 + k1 / z, -(k1 * (1 + 2 / z**2) + k0 / z), k0 * (1 + 3 / z**2) + k1 * (2 / z + 6 / z**3)]


@pytest.mark.slow
def test_gaussian_avatar_is_annihilated():
    mp.mp.dps = 25
    L = _LG()
    for tv in (mp.mpf("0.9"), mp.mpf("2.3")):
        vals = [mp.quad(lambda y, n=n: y**n * _k0_derivs(tv * y)[n] * mp.besselk(0, y) * mp.exp(-y * y), [0, 0.5, 1, 2, 4, 8, mp.inf])
                for n in range(5)]
        terms = [sp.lambdify(t, c, "mpmath")(tv) * v for c, v in zip(L, vals)]
        assert abs(sum(terms)) / max(abs(v) for v in terms) < mp.mpf(10) ** -18


def test_gaussian_companion_has_gaussian_growth():
    mp.mp.dps = 20
    bes = {"I0": lambda z: mp.besseli(0, z), "K0": lambda z: mp.besselk(0, z)}[COMPANION]
    logJ = {}
    for tv in (mp.mpf(10), mp.mpf(20)):
        logJ[tv] = mp.log(mp.quad(lambda y: bes(tv * y) * mp.besselk(0, y) * mp.exp(-y * y), [0, 1, 2, 4, 8, 16, 32, mp.inf]))
    # a Fuchsian solution grows at most like a power: log J(20) - log J(10) <= d log 2 for some modest d
    assert logJ[mp.mpf(20)] - logJ[mp.mpf(10)] > 50
    for tv, lj in logJ.items():
        assert abs(lj - (tv - 1) ** 2 / 4) < 2 * mp.log(tv) + 2


def _ap2(j, p):
    jq = sp.Rational(j)
    if jq.q % p == 0:
        return None
    jp = jq.p * pow(jq.q, -1, p) % p
    if jp == 0 or (jp - 1728) % p == 0:
        return None
    A, B = 3 * jp * (1728 - jp) % p, 2 * jp * (1728 - jp) ** 2 % p
    if (4 * A**3 + 27 * B**2) % p == 0:
        return None
    sq = {}
    for y in range(p):
        sq[y * y % p] = sq.get(y * y % p, 0) + 1
    n = sum(sq.get((X**3 + A * X + B) % p, 0) for X in range(p)) + 1
    return (p + 1 - n) ** 2


def _agree(j1, j2, x0s, primes):
    agree = total = 0
    for x0 in x0s:
        for p in primes:
            a, b = _ap2(j1.subs(x, x0), p), _ap2(j2.subs(x, x0), p)
            if a is None or b is None:
                continue
            total += 1
            agree += a == b
    return agree, total


def _j(z2):
    return sp.cancel(64 * (z2 + 3) ** 3 / (1 - z2) ** 2)


def test_momentum_slot_is_the_geminals():
    g = sp.Rational(3, 2)
    J1 = _j((x + 1 - g**2) ** 2 / (4 * x))
    J2 = (x**2 - sp.Rational(11, 2) * x + sp.Rational(25, 16)) ** 3 * 4096 / (x**2 * (4 * x - 1) * (4 * x - 25))
    alt = {"1": _j((x + g**2 - 1) ** 2 / (4 * x * g**2)), "t": _j((1 + g**2 - x) ** 2 / (4 * g**2)),
           "gamma": J1}
    primes = list(sp.primerange(5, 60))
    x0s = [sp.Rational(k) for k in (2, 3, 5, 7)]
    a, n = _agree(J1, J2, x0s, primes)
    assert n >= 20 and a == n                                  # the two search solutions are isogenous
    l = sp.Symbol("l")
    a, n = _agree((256 * (l**2 - l + 1) ** 3 / (l**2 * (l - 1) ** 2)).subs(l, x),
                  (16 * (l**2 + 14 * l + 1) ** 3 / (l * (1 - l) ** 4)).subs(l, x), x0s, primes)
    assert n >= 20 and a == n                                  # Legendre and its 2-isogenous form
    for slot in ALT_MOMENTUM:
        a, n = _agree(J1, alt[slot], x0s, primes)
        assert n >= 20 and a < n                               # not isogenous: some good prime disagrees
