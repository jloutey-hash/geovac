"""Paper 59 sec:density_ring: the Gaussian geminal exp(-r^2) on both edges of the correlated triangle (GY1-GY3).

Record docs/predictions/2026-10-23_gaussian_geminal_triangle.md (committed b84e1231 before computing).
Written as a separate step from the result; no debug/ import.  Data tests/data/paper59_gy_2026_10_24.json.

Four guards.
(1) The stored order-nine operator annihilates the kite's Taylor series of Y_G (and its rational part B alone) through
    the coefficient recurrence, to far below the precision a wrong operator leaves.  Wrong answers rejected: one
    coefficient changed by 1/200; the operator applied with the wrong monomial offset (t^2 in place of t^3).
(2) It factors exactly as F2 o F1, orders five and four.  Wrong answer rejected: F1 with one coefficient changed.
(3) It is irregular at x = infinity, with formal growth exp(t^2/4) and exp(t^2/8): the Newton polygon's top edge has
    slope -1 and characteristic polynomial C^5 (2C - 1)^2 (4C - 1)^2 up to a constant.  The same code reports F_3's
    operator as Fuchsian there (it would accept no irregular point it was not shown).
(4) x = 1 is an ordinary point (exponents 0..8) and x = 4 a genuine one (exponent -3/2).  The same code reports x = 4
    as non-ordinary, so a test of "ordinary" can fail.
Not re-run here: minimality, the irreducibility of F1 and F2 and the absence of an F_3 piece (bounded searches,
record only), and the quadrature on Y_G.
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

x = sp.Symbol("x")
u = sp.Symbol("u")
DATA = json.loads((Path(__file__).parent / "data" / "paper59_gy_2026_10_24.json").read_text())


def _polys(lst, var=x):
    return [sp.Poly(sp.sympify(s.replace("^", "**"), locals={var.name: var}), var) for s in lst]


def _recurrence_residual(Q, phis, offset=3, nmax=120):
    """max over N of |sum_i Q_i(offset + 2(N - i)) c_(N-i)| / sum |terms|, c_n = binom(-3/2, n) phi_n."""
    mp.mp.dps = 130
    c, b = [], mp.mpf(1)
    for n, s in enumerate(phis[:nmax]):
        c.append(b * mp.mpf(s))
        b = b * (mp.mpf(-3) / 2 - n) / (n + 1)
    worst = mp.mpf(0)
    for N in range(len(Q), nmax):
        terms = [mp.mpf(int(Q[i].eval(offset + 2 * (N - i)))) * c[N - i] for i in range(len(Q))]
        size = sum(abs(v) for v in terms)
        worst = max(worst, abs(sum(terms)) / size)
    return worst


def test_operator_annihilates_the_gaussian_series():
    Q = _polys(DATA["L9G_Q"], u)
    for key in ("phi_n", "phib_n"):
        assert _recurrence_residual(Q, DATA[key]) < mp.mpf(10) ** -90, key
    # one coefficient changed by 1/200
    q = Q[5].as_expr()
    lead = sp.Poly(q, u).coeffs()[0]
    Qp = list(Q)
    Qp[5] = sp.Poly(q + sp.Rational(1, 200) * lead * u ** sp.degree(q, u), u)
    assert _recurrence_residual(Qp, DATA["phi_n"]) > mp.mpf(10) ** -20
    # the wrong monomial offset
    assert _recurrence_residual(Q, DATA["phi_n"], offset=2) > mp.mpf(10) ** -20


def _th(p, m):
    for _ in range(m):
        p = sp.Poly(2 * x * p.diff(x).as_expr(), x)
    return p


def _compose(A, B):
    out = [sp.Poly(0, x) for _ in range(len(A) + len(B) - 1)]
    for i, a in enumerate(A):
        for j, b in enumerate(B):
            for m in range(i + 1):
                out[i - m + j] += a * sp.binomial(i, m) * _th(b, m)
    return out


def _right_divides(L, R):
    rem = [c.as_expr() for c in L]
    Q = {}
    while len(rem) >= len(R):
        k = len(rem) - len(R)
        q = sp.cancel(rem[-1] / R[-1].as_expr())
        Q[k] = q
        tk = _compose([sp.Poly(0, x)] * k + [sp.Poly(1, x)], R)
        rem = [sp.cancel(rem[i] - q * tk[i].as_expr()) for i in range(len(rem))]
        while rem and rem[-1] == 0:
            rem.pop()
    return all(r == 0 for r in rem), Q


def test_two_pieces_compose_to_the_operator():
    L, F1, F2 = (_polys(DATA[k]) for k in ("L9G_theta", "F1_theta", "F2_theta"))
    assert (len(L) - 1, len(F1) - 1, len(F2) - 1) == (9, 4, 5)
    ok, Q = _right_divides(L, F1)
    assert ok
    ratios = {sp.cancel(Q[k] / F2[k].as_expr()) for k in Q if F2[k].as_expr() != 0}
    assert len(ratios) == 1
    F1p = list(F1)
    F1p[0] = F1[0] + sp.Poly(1, x)
    assert not _right_divides(L, F1p)[0]


def _newton_top_edge(theta_ops):
    """points (k, deg c_k); the edge of maximal k + deg (slope -1 in this normalisation), its characteristic poly."""
    C = sp.Symbol("C")
    pts = [(k, c.degree()) for k, c in enumerate(theta_ops) if not c.is_zero]
    top = max(k + d for k, d in pts)
    edge = [k for k, d in pts if k + d == top]
    height = max(d for _, d in pts) - pts[-1][1]
    return height, edge, sp.Poly(sum(sp.LC(theta_ops[k]) * C**k for k in edge), C)


def test_irregular_at_infinity_with_two_gaussian_rates():
    L = _polys(DATA["L9G_theta"])
    height, edge, cp = _newton_top_edge(L)
    assert height == 4 and edge == [5, 6, 7, 8, 9]
    C = sp.Symbol("C")
    assert sp.factor(cp.as_expr() / sp.LC(cp)) == sp.factor(C**5 * (C - sp.Rational(1, 2)) ** 2 * (C - sp.Rational(1, 4)) ** 2)
    # exp(c x) with c = C/2: rates 1/4 and 1/8
    # the same code on F_3's operator 4 theta^3 - x (theta + 1)^3 (Fuchsian): height 0
    th = sp.Symbol("th")
    f3 = sp.Poly(4 * th**3 - x * (th + 1) ** 3, th)
    f3_theta = [sp.Poly(f3.as_expr().coeff(th, k), x) for k in range(4)]
    assert _newton_top_edge(f3_theta)[0] == 0


def _indicial_roots(theta_ops, x0):
    """indicial polynomial at x = x0 != 0 in the D = d/dx form."""
    order = len(theta_ops) - 1
    thk = {0: {0: sp.Integer(1)}}
    for k in range(1, order + 1):
        new = {}
        for j, c in thk[k - 1].items():
            new[j] = new.get(j, 0) + 2 * x * sp.diff(c, x)
            new[j + 1] = new.get(j + 1, 0) + 2 * x * c
        thk[k] = {j: sp.expand(v) for j, v in new.items()}
    Pd = [0] * (order + 1)
    for k, c in enumerate(theta_ops):
        for j, v in thk[k].items():
            Pd[j] = sp.expand(Pd[j] + c.as_expr() * v)
    s, r = sp.symbols("s r")
    P = [sp.Poly(sp.expand(p.subs(x, x0 + s)), s) for p in Pd]
    v = [min(m[0] for m in p.monoms()) for p in P]
    delta = min(v[k] - k for k in range(order + 1))
    poly = sum(P[k].coeff_monomial(s ** v[k]) * sp.ff(r, k) for k in range(order + 1) if v[k] - k == delta)
    return delta == v[-1] - order, sorted(sp.roots(sp.Poly(sp.expand(poly), r)).keys())


def test_x_equals_one_ordinary_and_x_equals_four_genuine():
    L = _polys(DATA["L9G_theta"])
    fuchs1, roots1 = _indicial_roots(L, sp.Integer(1))
    assert fuchs1 and roots1 == list(range(9))
    fuchs4, roots4 = _indicial_roots(L, sp.Integer(4))
    assert fuchs4 and sp.Rational(-3, 2) in roots4
    assert roots4 != list(range(9))
