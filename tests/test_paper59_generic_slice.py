"""Paper 59 sec:three_electron, the second slice of the chain fibre: scales (t;1,2).

Record docs/predictions/2026-10-21_generic_slice.md (GS1-GS5, committed a606d28f before computing).
Written as a separate step from the result, self-contained (no debug/ import); the exact gauge found by
search is loaded from tests/data/paper59_ladder_2026_10_21.json and re-verified here.

G_b(t) = int_0^oo x K1(t x) J(x) J(2x) dx, J = 3 int_0^x K0 - 3 x K0 - x^2 K1 (the Slater leaf weight of
eq:chain_slater_slice).  The Mellin route gives L of order 17; L = M14 o E2 o (theta+2) exactly, E2 the
least common left multiple of the two leaf pieces; M14 has a right factor R of order 4 equivalent to the
operator R_V of V(t) = int K0(t x) K0(x) K0(2x) dx; W7 = R o E2 o (theta+2) annihilates G_b; and no
rational map splits R_V off W7.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import pytest
import sympy as sp

t, u, z, s, y = sp.symbols("t u z s y")
DATA = json.loads((Path(__file__).parent / "data" / "paper59_ladder_2026_10_21.json").read_text())["generic_slice"]
B_SCALE = sp.Integer(2)

# The registered b = 1 slice operator (tests/test_paper59_chain_slater_slice.py, eq:chain_slater_slice)
QG1 = [-4 * (u - 5) * (u - 3) * (u - 1) ** 2 * (u + 1) * (u + 2) * (u + 5),
       (u - 3) * (u - 1) * (u + 1) * (u + 2) * (5 * u**3 + 88 * u**2 + 460 * u + 708),
       -(u - 1) * (u + 2) * (u**5 + 36 * u**4 + 463 * u**3 + 2718 * u**2 + 7396 * u + 7680),
       -6 * (u + 2) * (u + 4) * (u + 7) * (u + 8) * (u + 12)]


# ------------------------------------------------------------------ Mellin route
def _tensor_rel(rel1, rel2):
    """theta^N1 f = sum rel1[j] theta^j f, theta^N2 g = sum rel2[j] theta^j g (theta z = 2 z): the
    minimal relation sum_k v[k](z) theta^k (f g) = 0 in the product basis."""
    N1, N2 = len(rel1), len(rel2)
    keys = [(i, j) for i in range(N1) for j in range(N2)]
    idx = {k: n for n, k in enumerate(keys)}

    def theta(vec):
        new = [sp.Integer(0)] * len(keys)
        for (i, j), cf in zip(keys, vec):
            if cf == 0:
                continue
            new[idx[(i, j)]] += 2 * z * sp.diff(cf, z)
            if i + 1 < N1:
                new[idx[(i + 1, j)]] += cf
            else:
                for a in range(N1):
                    new[idx[(a, j)]] += cf * rel1[a]
            if j + 1 < N2:
                new[idx[(i, j + 1)]] += cf
            else:
                for b in range(N2):
                    new[idx[(i, b)]] += cf * rel2[b]
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


def _mellin_t_ode(v, alpha, beta, r):
    """Phi(s) = A(alpha - s) E(s + beta), A = M[f], E(sig+2) = r(sig) E(sig): returns Q with
    sum_i t^(2i) Q_i(theta_t) F = 0 (the route of Paper 59 sec:four_density)."""
    K = max(sp.Poly(e, z).degree() for e in v if e != 0)
    P = [sp.expand(sum(sp.Poly(v[j], z).coeff_monomial(z**k) * u**j for j in range(len(v))))
         for k in range(K + 1)]
    R = []
    for i in range(K + 1):
        k = K - i
        S = s + 2 * K
        R.append(sp.expand(P[k].subs(u, -(alpha - S + 2 * k)) * sp.prod([r(S + beta - 2 * l) for l in range(1, k + 1)])))
    g = sp.gcd_list([e for e in R if e != 0])
    R = [sp.expand(sp.cancel(e / g)) for e in R]
    return [sp.factor(sp.expand(R[i].subs(s, -u - 2 * i))) for i in range(K + 1)]


def _L_of(b):
    rel = [sp.Integer(0), z - 5, sp.Integer(6)]
    relb = [sp.Integer(0), b**2 * z - 5, sp.Integer(6)]
    return _mellin_t_ode(_tensor_rel(rel, relb), 2, 0, lambda sg: sg**2 - 1)


def _V_of(b):
    return _mellin_t_ode(_tensor_rel([z, sp.Integer(0)], [b**2 * z, sp.Integer(0)]), 1, 0, lambda sg: sg**2)


# ------------------------------------------------------------------ operators, theta-form lists
def _coeff_list(Q):
    n = max(sp.Poly(q, u).degree() for q in Q) + 1
    out = [sp.Integer(0)] * n
    for i, q in enumerate(Q):
        for k, c in enumerate(sp.Poly(sp.expand(q), u).all_coeffs()[::-1]):
            out[k] += c * t ** (2 * i)
    return [sp.expand(c) for c in out]


def _cth(c):
    return sp.cancel(t * sp.diff(c, t))


def _compose(A, B):
    out = [sp.Integer(0)] * (len(A) + len(B) - 1)
    for k, a in enumerate(A):
        if a == 0:
            continue
        for j, b in enumerate(B):
            bi = b
            for i in range(k + 1):
                out[k - i + j] += a * sp.binomial(k, i) * bi
                bi = _cth(bi)
    return [sp.cancel(c) for c in out]


def _rdiv(L, R):
    L = [sp.cancel(c) for c in L]
    r = len(R) - 1
    Qo = [sp.Integer(0)] * (len(L) - r)
    for k in range(len(L) - 1, r - 1, -1):
        c = sp.cancel(L[k] / R[-1])
        Qo[k - r] = c
        if c == 0:
            continue
        sub = _compose([0] * (k - r) + [c], R)
        L = [sp.cancel(L[i] - (sub[i] if i < len(sub) else 0)) for i in range(len(L))]
    return Qo, L[:r]


def _first_order(h):
    n, d = sp.fraction(sp.cancel(_cth(h) / h))
    return [sp.expand(-n), sp.expand(d)]


TH2 = [sp.Integer(2), sp.Integer(1)]
HA, HB = t / (1 - t**2) ** 4, t / (4 - t**2) ** 4


def _E2():
    return [sp.sympify(c, locals={"t": t}) for c in DATA["E2"]]


def _annihilator_of_image(Bg, RV):
    """R with R o Bg = 0 modulo the right ideal of RV (order 4): the annihilator of Bg(Sol RV)."""
    S = [sp.cancel(-c / RV[4]) for c in RV[:4]]

    def th_vec(p):
        out = [sp.cancel(_cth(p[0]) + p[-1] * S[0])]
        for j in range(1, 4):
            out.append(sp.cancel(p[j - 1] + _cth(p[j]) + p[-1] * S[j]))
        return out
    vecs = [list(Bg)]
    for _ in range(4):
        vecs.append(th_vec(vecs[-1]))
    ns = sp.Matrix([[vecs[k][i] for k in range(5)] for i in range(4)]).nullspace()
    assert len(ns) == 1
    return [sp.cancel(x / ns[0][4]) for x in ns[0]]


@pytest.fixture(scope="module")
def L2():
    return _coeff_list(_L_of(B_SCALE))


def test_mellin_operator_and_controls(L2):
    """Rejects a wrong tensor construction (b = 1 must give the registered slice operator exactly) and
    'the generic slice's Mellin operator has order at most 10' (GS1 as registered: it is 17).  Also V's
    operator: order 4 at b = 2, F_3's operator at b = 1."""
    Q1 = _L_of(sp.Integer(1))
    assert len(Q1) == len(QG1) and all(sp.expand(a - b) == 0 for a, b in zip(Q1, QG1))
    assert len(L2) - 1 == 17
    QV = _V_of(B_SCALE)
    assert [sp.expand(q) for q in QV] == [sp.expand(sp.sympify(q, locals={"u": u})) for q in DATA["QV"]]
    assert max(sp.Poly(q, u).degree() for q in QV) == 4
    assert [sp.expand(q) for q in _V_of(sp.Integer(1))] == [4 * u**3, -sp.expand((u + 1) ** 3)]


@pytest.mark.slow
def test_leaf_pieces_enter_as_their_lclm(L2):
    """Rejects 'the two leaf pieces compose one after the other' and a missing leaf piece: theta + 2 is
    a right factor; each leaf piece is a right factor of L/(theta+2); neither order of composition is;
    their LCLM E2 is, and E2 kills both leaf solutions."""
    L16, rem = _rdiv(L2, TH2)
    assert all(r == 0 for r in rem)
    Ra, Rb = _first_order(HA), _first_order(HB)
    for R1 in (Ra, Rb):
        assert all(sp.cancel(r) == 0 for r in _rdiv(L16, R1)[1])
    for first, second in ((Ra, Rb), (Rb, Ra)):
        q, rem = _rdiv(L16, first)
        assert not all(sp.cancel(r) == 0 for r in _rdiv(q, second)[1])
    E2 = _E2()
    for h in (HA, HB):
        assert sp.cancel(E2[0] * h + E2[1] * _cth(h) + E2[2] * _cth(_cth(h))) == 0
    assert all(sp.cancel(r) == 0 for r in _rdiv(L16, E2)[1])


@pytest.mark.slow
def test_k3_piece_is_v_operator(L2):
    """Rejects 'no piece of M14 = L/E3 is equivalent to V's operator' and a wrong gauge: the stored gauge
    B maps Sol(R_V) into Sol(M14), exactly.  Its image's annihilator R has order 4 and is an exact right
    factor of M14."""
    E3 = _compose(_E2(), TH2)
    M14, rem = _rdiv(L2, E3)
    assert all(sp.cancel(r) == 0 for r in rem)
    Bg = [sp.sympify(c, locals={"t": t}) for c in DATA["B"]]
    RV = _coeff_list([sp.sympify(q, locals={"u": u}) for q in DATA["QV"]])
    R = _annihilator_of_image(Bg, RV)
    assert len(R) == 5
    assert all(sp.cancel(r) == 0 for r in _rdiv(M14, R)[1])


# ------------------------------------------------------------------ numerics on G_b
def _struve_L(X):
    with mp.extradps(20):
        return mp.pi / 2 * X * (mp.besselk(0, X) * mp.struvel(-1, X) + mp.besselk(1, X) * mp.struvel(0, X))


def _J(X):
    return 3 * _struve_L(X) - 3 * X * mp.besselk(0, X) - X * X * mp.besselk(1, X)


def _grid(dps, hinv):
    mp.mp.dps = dps
    h = mp.mpf(1) / hinv
    nodes = []
    for k in range(int(-5 * hinv), int(5.2 * hinv) + 1):
        sv = k * h
        x = mp.exp(sv - mp.exp(-sv))
        nodes.append((x, h * x * (1 + mp.exp(-sv)), x * _J(x) * _J(2 * x)))
    return nodes


def _theta_k1(n):
    al, be = [sp.Integer(1)], [sp.Integer(0)]
    for _ in range(n):
        a, b = al[-1], be[-1]
        al.append(sp.expand(y * sp.diff(a, y) + b * (1 + y**2)))
        be.append(sp.expand(a + y * sp.diff(b, y)))
    return al, be


def _apply(nodes, Lc, tv):
    """(sum_k Lc[k] theta^k G_b)(tv) relative to its largest term."""
    tv = mp.mpf(tv)
    al, be = _theta_k1(len(Lc) - 1)
    fa = [sp.lambdify(y, a, "mpmath") for a in al]
    fb = [sp.lambdify(y, b, "mpmath") for b in be]
    I = [mp.mpf(0)] * len(Lc)
    for x, w, fx in nodes:
        yy = tv * x
        k0, k1 = mp.besselk(0, yy), mp.besselk(1, yy)
        th1 = -yy * k0 - k1
        for j in range(len(Lc)):
            I[j] += w * fx * (fa[j](yy) * k1 + fb[j](yy) * th1)
    cs = [mp.mpf(sp.N(sp.sympify(c).subs(t, sp.Rational(str(tv))), mp.mp.dps + 10)) for c in Lc]
    terms = [c * v for c, v in zip(cs, I)]
    return abs(sum(terms)) / max(abs(x) for x in terms)


@pytest.mark.slow
def test_w7_annihilates_and_needs_k3():
    """Rejects 'G_b lies below the K3 piece' (E3 G_b = 0) and a wrong order-four piece (a 1/200 change
    in one coefficient of R must spoil the annihilation): W7 = R o E3 kills G_b to 1e-25 of its largest
    term at 32 digits, at t = 1.3 and 2.5."""
    E3 = _compose(_E2(), TH2)
    Bg = [sp.sympify(c, locals={"t": t}) for c in DATA["B"]]
    RV = _coeff_list([sp.sympify(q, locals={"u": u}) for q in DATA["QV"]])
    R = _annihilator_of_image(Bg, RV)
    W7 = _compose(R, E3)
    Rp = list(R)
    Rp[1] = R[1] * sp.Rational(201, 200)
    W7p = _compose(Rp, E3)
    nodes = _grid(32, 48)
    for tv in ("1.3", "2.5"):
        assert _apply(nodes, W7, tv) < mp.mpf(10) ** -25
        assert _apply(nodes, W7p, tv) > mp.mpf(10) ** -7
        assert _apply(nodes, E3, tv) > mp.mpf(10) ** -2


# ------------------------------------------------------------------ no split, mod p
P = (1 << 61) - 1


def _fr(x):
    x = Fraction(x)
    return x.numerator % P * pow(x.denominator % P, P - 2, P) % P


def _ev(q, x):
    val = Fraction(0)
    for (m,), c in sp.Poly(q, u).terms():
        val += Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)) * Fraction(x) ** m
    return _fr(val)


def _frob(Q, lam, NT):
    c = [1]
    for n in range(1, NT // 2):
        acc = sum(_ev(Q[i], lam + 2 * (n - i)) * c[n - i] for i in range(1, len(Q)) if n - i >= 0) % P
        d = _ev(Q[0], lam + 2 * n)
        c.append((-acc) * pow(d, P - 2, P) % P)
    return {lam + 2 * n: cn for n, cn in enumerate(c)}


def _apply_L(Q, ser):
    out = {}
    for i, q in enumerate(Q):
        cs = [(m[0], _fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)))) for m, c in sp.Poly(q, u).terms()]
        for e, c in ser.items():
            val = sum(cc * pow(e % P, m, P) for m, cc in cs) % P
            if val:
                out[e + 2 * i] = (out.get(e + 2 * i, 0) + val * c) % P
    return {e: c for e, c in out.items() if c}


def _rank(rows):
    rows = [r[:] for r in rows]
    rk = 0
    for col in range(len(rows[0])):
        piv = next((i for i in range(rk, len(rows)) if rows[i][col]), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        iv = pow(rows[rk][col], P - 2, P)
        rows[rk] = [x * iv % P for x in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][col]:
                f = rows[i][col]
                rows[i] = [(x - f * yv) % P for x, yv in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def _to_Q(Lc):
    den = sp.lcm([sp.fraction(sp.cancel(c))[1] for c in Lc])
    pol = [sp.Poly(sp.expand(sp.cancel(c * den)), t) for c in Lc]
    deg = max(p.degree() for p in pol)
    Q = [sp.expand(sum(p.coeff_monomial(t ** (2 * i)) * u**k for k, p in enumerate(pol))) for i in range(deg // 2 + 1)]
    while Q and Q[-1] == 0:
        Q.pop()
    return Q


def _gauge_nullity(Q, f, Bp, a, Dg, NT):
    den = sp.Poly(sp.expand((t**2 - 1) ** Bp * (t**2 - 4) ** Bp * (t**2 - 9) ** Bp), t).all_coeffs()[::-1]
    cs = [_fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q))) for c in den]
    w = [pow(cs[0], P - 2, P)]
    for n in range(1, NT):
        acc = sum(cs[k] * w[n - k] for k in range(1, min(n, len(cs) - 1) + 1)) % P
        w.append((-acc) * w[0] % P)
    top = min(f) + NT - 2 * len(Q) - Dg - a - 8
    cols = []
    for j in range(4):
        g = {}
        for e1, c1 in f.items():
            cj = c1 * pow(e1 % P, j, P) % P
            for e2, c2 in enumerate(w):
                if c2 and e1 + e2 < NT:
                    g[e1 + e2] = (g.get(e1 + e2, 0) + cj * c2) % P
        for d in range(-a, Dg + 1):
            cols.append(_apply_L(Q, {e + d: c for e, c in g.items()}))
    exps = sorted(e for e in set().union(*cols) if e < top)
    return len(cols) - _rank([[col.get(e, 0) for col in cols] for e in exps])


@pytest.mark.slow
def test_k3_piece_does_not_split():
    """Rejects 'G_b is a rational combination of V-type functions plus elementary ones' (a split of the
    K3 piece off W7), and a search that finds nothing whatever it is given: the planted split
    E3 o R_V is found at the same box, which exceeds the leading-term bound (poles of order 4 at t = 0
    and 6 at t = +-1, +-2, +-3, numerator excess 1 at infinity)."""
    E3 = _compose(_E2(), TH2)
    Bg = [sp.sympify(c, locals={"t": t}) for c in DATA["B"]]
    QV = [sp.sympify(q, locals={"u": u}) for q in DATA["QV"]]
    RV = _coeff_list(QV)
    R = _annihilator_of_image(Bg, RV)
    QW = _to_Q(_compose(R, E3))
    QX = _to_Q(_compose(E3, RV))
    Bp, a, Dg = 10, 8, 65
    NT = 2 * (4 * (Dg + a + 1)) + 300
    f = _frob(QV, 2, NT)
    assert _gauge_nullity(QX, f, Bp, a, Dg, NT) == 1
    assert _gauge_nullity(QW, f, Bp, a, Dg, NT) == 0
