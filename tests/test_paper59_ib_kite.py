"""Paper 59 sec:density_ring: the kite operators from integration by parts (record IB, 2026-10-24).

Record docs/predictions/2026-10-24_ibp_kite.md (committed e7679061 before computing).  Written as a separate step from
the result; no debug/ import.  Data tests/data/paper59_ib_2026_10_24.json.

Five guards.
(1) The off-slice kite operator (order 11, z-degree up to 113, integer coefficients) annihilates the stored Taylor
    series of the kite at (1, 2; 3/2, 3) through the coefficient recurrence.  Wrong answer rejected: one coefficient
    changed by 1/200.
(2) Its singular points: the leading coefficient vanishes at the eight crossed cuts and not at -25/4, -25; the local
    exponents at -1/4 contain -1/2, at -1 the exponent 0 is double, at -25/4 they are 0..10.
(3) Analytic solutions counted modulo a prime: 8 of 11 at z = -1 (genuine), 11 of 11 at z = -25/4 (ordinary) and at
    a root of the irreducible remaining factor of the leading coefficient (apparent).  The same code reports fewer
    than 11 at a genuine point, so the test of "apparent" can fail.
(4) The gamma = 3/2 operator from integration by parts equals the series-guessed one modulo P, after making both
    monic.  Wrong answer rejected: one coefficient of the IBP operator changed.
(5) The reduction method on the one-loop bubble in D = 3, by a separate exact implementation: dB/dz reduced to B and
    the two tadpoles satisfies the arctan closed form.  Wrong answer rejected: the same reduction in D = 4.
Not re-run here: the reductions themselves (minutes each) and the gamma = 1/2, 1, 2, 3 operators (record only).
"""
from __future__ import annotations

import json
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp
import sympy as sp
from flint import fmpq, fmpq_poly, fmpz_poly, nmod_mat, nmod_poly

DATA = json.loads((Path(__file__).parent / "data" / "paper59_ib_2026_10_24.json").read_text())
CROSSED = ["-1/4", "-1", "-9/4", "-4", "-9", "-49/4", "-16", "-81/4"]
DIRECT = ["-25/4", "-25"]


def _pf_ops():
    return [[int(v) for v in p] for p in DATA["pf_op"]]


def _series_residual(op, xs, dps=120):
    mp.mp.dps = dps
    x = [mp.mpf(s) for s in xs]
    n, order, worst = len(x), len(op) - 1, mp.mpf(0)
    for N in range(n - order - 1):
        acc = size = mp.mpf(0)
        for k, co in enumerate(op):
            for i, c in enumerate(co):
                m = N - i
                if m < 0 or c == 0 or m + k >= n:
                    continue
                ff = 1
                for j in range(k):
                    ff *= m + k - j
                t = mp.mpf(c.numerator) / c.denominator * ff * x[m + k] if isinstance(c, F) else c * ff * x[m + k]
                acc += t
                size += abs(t)
        if size:
            worst = max(worst, abs(acc) / size)
    return worst


def test_pf_operator_annihilates_its_series():
    op = _pf_ops()
    assert _series_residual(op, DATA["pf_series"]) < mp.mpf(10) ** -90
    bad = [list(p) for p in op]
    bad[5][-1] = F(bad[5][-1]) * (1 + F(1, 200))
    bad = [[F(v) for v in p] for p in bad]
    assert _series_residual(bad, DATA["pf_series"]) > mp.mpf(10) ** -60


def _exponents(op, z0: F):
    order = len(op) - 1
    sh = fmpq_poly([fmpq(z0.numerator, z0.denominator), 1])
    Ps = [fmpq_poly(p)(sh) for p in op]
    v = [min((i for i, c in enumerate(q.coeffs()) if c != 0), default=10**9) for q in Ps]
    delta = min(v[k] - k for k in range(order + 1))
    r = sp.Symbol("r")
    ind = sum(sp.Rational(int(Ps[k].coeffs()[v[k]].p), int(Ps[k].coeffs()[v[k]].q)) * sp.ff(r, k)
              for k in range(order + 1) if v[k] - k == delta)
    return Ps[-1].coeffs()[0] == 0, sp.roots(sp.Poly(sp.expand(ind), r))


def test_pf_singular_points_are_the_crossed_cuts():
    op = _pf_ops()
    lead = fmpz_poly(op[-1])
    for z in CROSSED:
        q = F(z)
        assert fmpq_poly(op[-1])(fmpq(q.numerator, q.denominator)) == 0, z
    for z in DIRECT:
        q = F(z)
        assert fmpq_poly(op[-1])(fmpq(q.numerator, q.denominator)) != 0, z
    _, e = _exponents(op, F(-1, 4))
    assert sp.Rational(-1, 2) in e
    _, e = _exponents(op, F(-1))
    assert e.get(sp.Integer(0)) == 2
    vanish, e = _exponents(op, F(-25, 4))
    assert not vanish and sorted(e) == list(range(11))
    assert lead.degree() == 113


def _analytic_count(op, rho, p, N=40):
    order = len(op) - 1
    Ps = [nmod_poly([c % p for c in q], p)(nmod_poly([rho, 1], p)) for q in op]
    v = [min((i for i, c in enumerate(q.coeffs()) if int(c)), default=10**9) for q in Ps]
    delta = min(v[k] - k for k in range(order + 1))
    rows = []
    for m in range(N + delta + 1):
        row = [0] * (N + 1)
        for k, q in enumerate(Ps):
            for i, c in enumerate([int(c) for c in q.coeffs()]):
                n = m - i + k
                if c and 0 <= n <= N and n >= k:
                    ff = 1
                    for j in range(k):
                        ff = ff * (n - j) % p
                    row[n] = (row[n] + c * ff) % p
        rows.append(row)
    return (N + 1) - nmod_mat(rows, p).rank()


def test_genuine_ordinary_and_apparent_points():
    op = _pf_ops()
    p = 2305843009213693921
    red = lambda q: q.numerator * pow(q.denominator, -1, p) % p  # noqa: E731
    assert _analytic_count(op, red(F(-1)), p) == 8
    assert _analytic_count(op, red(F(-25, 4)), p) == 11
    fac = fmpz_poly(op[-1]).factor()[1]
    big = [f for f, e in fac if f.degree() > 1]
    assert [f.degree() for f in big] == [92]
    f = big[0]
    from sympy import prevprime
    q = p
    for _ in range(80):
        q = prevprime(q)
        roots = [r for r, m in nmod_poly([int(c) % q for c in f.coeffs()], q).factor()[1] if r.degree() == 1]
        if roots:
            rho = (-int(roots[0].coeffs()[0]) * pow(int(roots[0].coeffs()[1]), -1, q)) % q
            assert _analytic_count(op, rho, q) == 11
            break
    else:
        raise AssertionError("no prime with a root found")


def _guessed_to_D(Qs, P):
    u = sp.Symbol("u")
    coef = {}
    for i, s in enumerate(Qs):
        for (k,), c in sp.Poly(sp.sympify(s.replace("^", "**"), locals={"u": u}), u).terms():
            for j in range(k + 1):
                s2 = int(sp.functions.combinatorial.numbers.stirling(k, j, kind=2))
                if s2:
                    coef[(j, i + j)] = coef.get((j, i + j), 0) + int(c) * s2
    r = max(j for j, _ in coef)
    out = []
    for j in range(r + 1):
        deg = max(e for jj, e in coef if jj == j)
        out.append([coef.get((j, e), 0) % P for e in range(deg + 1)])
    return out


def _same_operator(A, B, P):
    """equal up to a left factor in GF(P)(z): a_k / a_r = b_k / b_r for every k, cross-multiplied."""
    if len(A) != len(B):
        return False
    a = [nmod_poly(p, P) for p in A]
    b = [nmod_poly(p, P) for p in B]
    return all(a[k] * b[-1] == b[k] * a[-1] for k in range(len(a)))


def test_ibp_operator_equals_guessed_at_gamma_three_halves():
    P = DATA["P"]
    ibp = [list(p) for p in DATA["g32_ibp_op"]]
    guessed = _guessed_to_D(DATA["g32_guessed_Q"], P)
    assert len(ibp) - 1 == 6
    assert _same_operator(ibp, guessed, P)
    ibp[3][7] = (ibp[3][7] + 1) % P
    assert not _same_operator(ibp, guessed, P)


# ---- (5) a separate exact reduction of the D = 3 bubble ----
def _bubble_relations(a, b, z, m1s, m2s, D):
    """IBP identities d/dl . v for v = l, p on I(a, b) = int (l^2 + m1^2)^-a ((l+p)^2 + m2^2)^-b, as dicts."""
    # l.l = D1 - m1s ; l.p = (D2 - D1 - z + m1s - m2s)/2 ; D1 lowers a, D2 lowers b
    def lp():
        return {(0, -1): F(1, 2), (-1, 0): F(-1, 2), None: (-z + m1s - m2s) / 2}

    def add(out, key, val):
        out[key] = out.get(key, 0) + val

    rels = []
    for v in ("l", "p"):
        out = {}
        if v == "l":
            add(out, (a, b), F(D))
        # from D1^-a: -2a (v.l) I(a+1, b)
        vl = {(-1, 0): F(1), None: -m1s} if v == "l" else lp()
        for k, c in vl.items():
            key = (a + 1, b) if k is None else (a + 1 + k[0], b + k[1])
            add(out, key, -2 * a * c)
        # from D2^-b: -2b (v.(l+p)) I(a, b+1)
        if v == "l":
            vlp = {(-1, 0): F(1), None: -m1s}
            for k, c in lp().items():
                vlp[k] = vlp.get(k, 0) + c
        else:
            vlp = dict(lp())
            vlp[None] = vlp[None] + z
        for k, c in vlp.items():
            key = (a, b + 1) if k is None else (a + k[0], b + 1 + k[1])
            add(out, key, -2 * b * c)
        rels.append({k: c for k, c in out.items() if c and not (k[0] <= 0 and k[1] <= 0)})
    return rels


def _bubble_dB_dz(z, m1, m2, D):
    m1s, m2s = m1 * m1, m2 * m2
    key = lambda n: (int(n[0] > 0) + int(n[1] > 0), max(n[0], 0) + max(n[1], 0) - min(n[0], 0) - min(n[1], 0),
                     -min(n[0], 0) - min(n[1], 0), n)  # noqa: E731
    rows = {}

    def add(eq):
        eq = dict(eq)
        while eq:
            top = max(eq, key=key)
            if top not in rows:
                c = eq.pop(top)
                rows[top] = {n: v / c for n, v in eq.items()}
                return
            f = eq.pop(top)
            for n, v in rows[top].items():
                eq[n] = eq.get(n, 0) - f * v
                if eq[n] == 0:
                    del eq[n]

    for a in range(-1, 4):
        for b in range(-1, 4):
            if a + b >= 1:
                for r in _bubble_relations(a, b, z, m1s, m2s, D):
                    add(r)
    # dB/dz = (1/(2z)) p.d/dp: -(1/z) (p.(l+p)) I(1, 2), p.(l+p) = l.p + z
    lp = {(0, -1): F(1, 2), (-1, 0): F(-1, 2), None: (-z + m1s - m2s) / 2 + z}
    eq = {}
    for k, c in lp.items():
        n = (1, 2) if k is None else (1 + k[0], 2 + k[1])
        eq[n] = eq.get(n, 0) - c / z
    out = {}
    while eq:
        top = max(eq, key=key)
        c = eq.pop(top)
        if top[0] <= 0 and top[1] <= 0 or c == 0:
            continue
        if top in rows:
            for n, v in rows[top].items():
                eq[n] = eq.get(n, 0) - c * v
        else:
            out[top] = out.get(top, 0) + c
    return out


def _bubble_check(D):
    mp.mp.dps = 40
    z, m1, m2 = F(3), F(1), F(3, 2)
    red = _bubble_dB_dz(z, m1, m2, D)
    M = mp.mpf(5) / 2
    B = lambda zz: 2 * mp.pi**2 / mp.sqrt(zz) * mp.atan(mp.sqrt(zz) / M)  # noqa: E731
    vals = {(1, 1): B(3), (1, 0): -2 * mp.pi**2 * 1, (0, 1): -2 * mp.pi**2 * mp.mpf(3) / 2}
    rhs = sum(mp.mpf(c.numerator) / c.denominator * vals[n] for n, c in red.items())
    lhs = mp.diff(B, 3)
    return abs(lhs - rhs) / abs(lhs)


def test_reduction_method_on_the_bubble():
    assert _bubble_check(3) < mp.mpf(10) ** -30
    assert _bubble_check(4) > mp.mpf(10) ** -6
