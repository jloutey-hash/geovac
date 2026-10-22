"""Paper 59 sec:four_density, the real three-dimensional four-density star.

Record docs/predictions/2026-10-21_four_density_star.md (FD1-FD3, committed 9a33fbef before computing).
Written as a separate step from the result; no debug/ import.  The exact pieces found by search (the
order-one piece X1, the order-four piece R4pp and the gauge A from F_4's operator) are loaded from
tests/data/paper59_ladder_2026_10_21.json and re-verified here.  Shared operator helpers come from
tests/test_paper59_generic_slice.py.

At the fingerprint fibre on the slice (t;1,1,1), S(t) = (16/pi^4) t^3 int_0^oo chi(x)^3 K0(t x) dx / x,
chi = int_0^x s K1(s) ds.  Its Mellin operator L has order 19, and W8 = R4pp o X1 o E2 o (theta - 3)
is an exact right factor that annihilates S, with E2 = (theta-3)^2 - t^2 theta (theta-2) equivalent to
the Legendre operator and R4pp equivalent to F_4's operator (eq:four_density_op).
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

import test_paper59_generic_slice as GS

t, u, z, y = GS.t, GS.u, GS.z, GS.y
DATA = json.loads((Path(__file__).parent / "data" / "paper59_ladder_2026_10_21.json").read_text())["star"]
REL_CHI = [sp.Integer(0), z - 3, sp.Integer(4)]       # theta^3 chi = (z-3) theta chi + 4 theta^2 chi
QC = [-36 * (u - 7) * (u - 5) ** 2 * (u - 3) ** 2 * (u - 2) * (u - 1),   # the chain, Paper 59 sec:density_ring
      (u - 5) * (u - 3) ** 2 * (u - 2) * (45 * u**3 + 54 * u**2 - 128 * u + 164),
      -(u - 3) * (u - 2) * (9 * u**5 + 72 * u**4 + 59 * u**3 + 88 * u**2 - 8 * u - 160),
      -20 * u * (u - 2) * (u + 1) * (u + 2) * (u + 4)]
Q4 = [9 * u**2 * (u - 2) ** 2, -2 * u**2 * (5 * u**2 + 1), (u + 1) ** 4]   # F_4's operator in t
TH3 = [sp.Integer(-3), sp.Integer(1)]
E2 = [sp.Integer(9), -6 + 2 * t**2, 1 - t**2]    # (theta-3)^2 - t^2 theta (theta-2)


def _sym_power(rel, N, m):
    """theta^N f = sum_{j<N} rel[j](z) theta^j f (theta z = 2 z): minimal relation for f^m."""
    keys = list(itertools.combinations_with_replacement(range(N), m))
    idx = {k: i for i, k in enumerate(keys)}

    def theta(vec):
        new = [sp.Integer(0)] * len(keys)
        for k, cf in zip(keys, vec):
            if cf == 0:
                continue
            new[idx[k]] += 2 * z * sp.diff(cf, z)
            for l in range(m):
                kk = list(k)
                kk[l] += 1
                if kk[l] < N:
                    new[idx[tuple(sorted(kk))]] += cf
                else:
                    for j in range(N):
                        k2 = kk.copy()
                        k2[l] = j
                        new[idx[tuple(sorted(k2))]] += cf * rel[j]
        return [sp.expand(c) for c in new]

    rows = [[sp.Integer(1)] + [sp.Integer(0)] * (len(keys) - 1)]
    for _ in range(len(keys)):
        rows.append(theta(rows[-1]))
        ns = sp.Matrix(rows).T.nullspace()
        if ns:
            v = ns[0]
            den = sp.lcm([sp.fraction(sp.together(e))[1] for e in v])
            v = [sp.factor(sp.cancel(e * den)) for e in v]
            g = sp.gcd_list([e for e in v if e != 0])
            return [sp.expand(sp.cancel(e / g)) for e in v]
    raise RuntimeError("no relation")


def _proportional(Qa, Qb):
    if len(Qa) != len(Qb):
        return False
    r = sp.cancel(sp.expand(Qa[0]) / sp.expand(Qb[0]))
    return r.is_number and all(sp.expand(a - r * b) == 0 for a, b in zip(Qa, Qb))


def _x1():
    """X1 = c0 + c1 delta (w = t^2, delta = theta/2) in theta-form."""
    w = sp.Symbol("w")
    c0, c1 = [sp.sympify(c, locals={"w": w}).subs(w, t**2) for c in DATA["X1_delta_w"]]
    return [sp.expand(c0), sp.expand(c1 / 2)]


@pytest.fixture(scope="module")
def Lstar():
    return GS._coeff_list(GS._mellin_t_ode(_sym_power(REL_CHI, 3, 3), -3, 3, lambda sg: sg**2))


def test_star_operator_and_controls(Lstar):
    """Rejects a wrong symmetric-power or Mellin construction (m = 2 must give the chain's operator; K0^3
    must give F_4's) and 'the star's Mellin operator has order at most 11' (FD1 as registered: 19)."""
    QCb = GS._mellin_t_ode(_sym_power(REL_CHI, 3, 2), -2, 3, lambda sg: sg**2)
    assert _proportional(QCb, QC)
    Q4b = GS._mellin_t_ode(_sym_power([z, sp.Integer(0)], 2, 3), 1, 0, lambda sg: sg**2)
    assert _proportional(Q4b, Q4)
    assert len(Lstar) - 1 == 19


def test_star_bottom_pieces(Lstar):
    """Rejects 'the one-leaf piece is algebraic, as in the chain' and a wrong bottom: theta - 3 and then
    E2 are exact right factors, and the injective map f -> t^3 (1 + theta) f sends the Legendre period
    2F1(1/2,1/2;1;t^2) into Sol(E2), so E2 is equivalent to F_2's (irreducible) operator."""
    L18, rem = GS._rdiv(Lstar, TH3)
    assert all(r == 0 for r in rem)
    assert all(sp.cancel(r) == 0 for r in GS._rdiv(L18, E2)[1])
    N = 30
    f = sum(sp.binomial(2 * k, k) ** 2 * sp.Rational(1, 16**k) * t ** (2 * k) for k in range(N))
    th = lambda g: sp.expand(t * sp.diff(g, t))
    g = sp.expand(t**3 * (f + th(f)))
    e = sp.expand(E2[0] * g + E2[1] * th(g) + E2[2] * th(th(g)))
    assert min(m[0] for m in sp.Poly(e, t).monoms()) >= 2 * N
    # the map f -> t^3 (1 + theta) f is injective on Sol(Legendre): its kernel 1/t is not a solution
    h = 1 / t
    assert sp.cancel(GS._cth(GS._cth(h)) - t**2 * (GS._cth(GS._cth(h)) + 2 * GS._cth(h) + h)) != 0


@pytest.mark.slow
def test_star_f4_piece_is_exact(Lstar):
    """Rejects 'no piece of the star's operator is equivalent to F_4's' and a wrong stored piece: X1 and
    then R4pp are exact right factors over theta - 3 and E2, and R4pp is exactly the annihilator of the
    image of Sol(F_4's operator) under the stored gauge A."""
    R4 = GS._compose(GS._compose(_x1(), E2), TH3)
    q, rem = GS._rdiv(Lstar, R4)
    assert all(sp.cancel(r) == 0 for r in rem)
    R4pp = GS._coeff_list([sp.sympify(c, locals={"u": u}) for c in DATA["R4pp_Q"]])
    assert all(sp.cancel(r) == 0 for r in GS._rdiv(q, R4pp)[1])
    A = [sp.sympify(c, locals={"t": t}) for c in DATA["A"]]
    image_ann = GS._annihilator_of_image(A, GS._coeff_list(Q4))
    monic = [sp.cancel(c / R4pp[4]) for c in R4pp]
    assert all(sp.cancel(a - b) == 0 for a, b in zip(image_ann, monic))


# ------------------------------------------------------------------ numerics on S
def _chi(x):
    with mp.extradps(20):
        return mp.pi * x / 2 * (mp.besselk(1, x) * mp.struvel(0, x) + mp.besselk(0, x) * mp.struvel(1, x))


def _grid(dps, hinv):
    mp.mp.dps = dps
    h = mp.mpf(1) / hinv
    out = []
    for k in range(int(-5 * hinv), int(5.2 * hinv) + 1):
        sv = k * h
        x = mp.exp(sv - mp.exp(-sv))
        out.append((x, h * x * (1 + mp.exp(-sv)), _chi(x) ** 3 / x))
    return out


def _theta_k0(n):
    """theta^j K0(y) = a_j K0 + b_j theta K0, theta^2 K0 = y^2 K0."""
    al, be = [sp.Integer(1)], [sp.Integer(0)]
    for _ in range(n):
        a, b = al[-1], be[-1]
        al.append(sp.expand(y * sp.diff(a, y) + b * y**2))
        be.append(sp.expand(a + y * sp.diff(b, y)))
    return al, be


def _apply(nodes, Lc, tv):
    """(sum_k Lc[k] theta^k S)(tv) relative to its largest term; S = t^3 I(t), theta^k (t^3 I) =
    t^3 sum_m C(k,m) 3^(k-m) theta^m I (the constant 16/pi^4 drops out of the ratio)."""
    tv = mp.mpf(tv)
    n = len(Lc) - 1
    al, be = _theta_k0(n)
    fa = [sp.lambdify(y, a, "mpmath") for a in al]
    fb = [sp.lambdify(y, b, "mpmath") for b in be]
    I = [mp.mpf(0)] * (n + 1)
    for x, w, gx in nodes:
        yy = tv * x
        k0, k1 = mp.besselk(0, yy), mp.besselk(1, yy)
        thk0 = -yy * k1
        for j in range(n + 1):
            I[j] += w * gx * (fa[j](yy) * k0 + fb[j](yy) * thk0)
    thS = [tv**3 * sum(sp.binomial(k, m) * 3 ** (k - m) * I[m] for m in range(k + 1)) for k in range(n + 1)]
    cs = [mp.mpf(sp.N(sp.sympify(c).subs(t, sp.Rational(str(tv))), mp.mp.dps + 10)) for c in Lc]
    terms = [c * v for c, v in zip(cs, thS)]
    return abs(sum(terms)) / max(abs(x) for x in terms)


@pytest.mark.slow
def test_w8_annihilates_star_and_needs_f4():
    """Rejects 'the star sits below the F_4 piece' (R4 S = 0) and a wrong order-four piece (a 1/200 change
    in R4pp's theta^2 coefficient must spoil the annihilation): W8 kills S to 1e-25 of its largest term
    at 32 digits, at t = 1.3 and 2.5."""
    R4 = GS._compose(GS._compose(_x1(), E2), TH3)
    R4pp = GS._coeff_list([sp.sympify(c, locals={"u": u}) for c in DATA["R4pp_Q"]])
    W8 = GS._compose(R4pp, R4)
    Rp = list(R4pp)
    Rp[2] = R4pp[2] * sp.Rational(201, 200)
    W8p = GS._compose(Rp, R4)
    nodes = _grid(32, 48)
    for tv in ("1.3", "2.5"):
        assert _apply(nodes, W8, tv) < mp.mpf(10) ** -25
        assert _apply(nodes, W8p, tv) > mp.mpf(10) ** -6
        assert _apply(nodes, R4, tv) > mp.mpf(10) ** -3
