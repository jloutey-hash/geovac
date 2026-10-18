"""Paper 59 sec:density_ring -- a ring of three densities does not climb the period ladder.

Record: docs/predictions/2026-10-17_density_loop.md (committed before the computation).
Each test names the wrong answer it rejects.

  eq:star_triangle      the loop integral, by an independent Feynman-parameter route
  eq:ring_one_center    the (1,1,1) value, against eq:ring_reduced and against Remiddi's master
  (1,2,3) relation      the 100-digit integer relation, re-checked at 40 digits
  eq:ring_op            L_G G = -3 t^2, by the momentum form with t-derivatives taken analytically
                        (independent of the Mellin route and of the Struve position form)
  factorization         L_G = M o (theta - r) exactly, r from (1 - t^2/4)^(-3/2), and the map
                        A = t^3 (1 + theta)/(4 - t^2) sends Legendre solutions to solutions of M
  F_2 operator          the two-density moment F_2(t) = int K0(x) K0(tx) dx solves that Legendre
                        equation, so the quotient is the two-density elliptic rung
  chain control         the chain's operator peels exactly to N2 o R3 o R1 o (theta - 2), with
                        eq:chain_k3_map sending F_3 solutions into R3; W = R3 o R1 o (theta - 2)
                        annihilates the chain itself, so N2 is spurious; and within the searched
                        rational maps the K3 piece does not split off the elementary ones
  ring source           L_G 1 = 3 t^2, so L_G (G + 1) = 0; (theta - r)(G + 1) is nonzero at two t,
                        on the momentum form and independently on the Struve form
  chain order four      R1 o S1 = 1 mod R3 and W = Q1 o P exactly; P C = 0 at 30 digits
  chain no-split proof  the exponent bound table, the 27-unknown system inconsistent, and the
                        same code finds a planted split
  correlated triangle  record docs/predictions/2026-10-17_ring_extensions.md (QG): the geminal
                        triangle's annihilator divides exactly as R3 o (theta - r), R3 tied to F_3's
                        operator by an exact map; on Y itself, from an independent kite evaluator,
                        the annihilator vanishes and (theta - r) Y does not: a K3 constituent
"""
from __future__ import annotations

import random

import mpmath as mp
import numpy as np
import pytest
import sympy as sp

t, u, w = sp.symbols("t u w", positive=True)

# eq:ring_op, as printed in the paper: sum_i t^(2i) Q_i(theta)
QG = [-4 * u * (u - 3) ** 2, 5 * u**3 - 7 * u**2 - 5 * u + 3, -u * (u + 2) * (u + 3)]


def _th(c):
    return sp.together(t * sp.diff(c, t))


def _coeff_list(Q):
    n = max(sp.Poly(q, u).degree() for q in Q) + 1
    out = [sp.Integer(0)] * n
    for i, q in enumerate(Q):
        for k, c in enumerate(sp.Poly(sp.expand(q), u).all_coeffs()[::-1]):
            out[k] += c * t ** (2 * i)
    return [sp.expand(c) for c in out]


def _apply(Lc, f):
    """sum_k Lc[k] theta^k f for a function f(t)."""
    acc, g = 0, f
    for c in Lc:
        acc += c * g
        g = _th(g)
    return sp.simplify(acc)


def _compose(A, B):
    out = [0] * (len(A) + len(B) - 1)
    for k, a in enumerate(A):
        for j, b in enumerate(B):
            bi = b
            for i in range(k + 1):
                out[k - i + j] += a * sp.binomial(k, i) * bi
                bi = _th(bi)
    return [sp.together(c) for c in out]


# ---------------------------------------------------------------- eq:star_triangle
def _triangle_feynman(x1, x2, x3):
    d2 = lambda a, b: sum((ai - bi) ** 2 for ai, bi in zip(a, b))
    s12, s13, s23 = d2(x1, x2), d2(x1, x3), d2(x2, x3)

    def f(a1, a2):
        a3 = 1 - a1 - a2
        if a3 <= 0:
            return mp.mpf(0)
        return (a1 * a2 * s12 + a1 * a3 * s13 + a2 * a3 * s23) ** mp.mpf(-1.5)

    val = mp.quad(lambda a1: mp.quad(lambda a2: f(a1, a2), [0, (1 - a1) / 2, 1 - a1]),
                  [0, mp.mpf(1) / 2, 1])
    return mp.pi**2 / 2 * val


@pytest.mark.slow
def test_star_triangle_normalisation():
    """Rejects a wrong constant in eq:star_triangle (pi^3/2, 2 pi^3, pi^2, ...) and a wrong
    dependence on the distances: the Feynman-parameter value must equal pi^3/(|x12||x23||x31|)."""
    mp.mp.dps = 15
    rng = random.Random(59)
    for _ in range(2):
        xs = [[mp.mpf(rng.uniform(-2, 2)) for _ in range(3)] for _ in range(3)]
        d = lambda a, b: mp.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))
        unique = mp.pi**3 / (d(xs[0], xs[1]) * d(xs[1], xs[2]) * d(xs[2], xs[0]))
        assert abs(_triangle_feynman(*xs) / unique - 1) < 1e-8


# ---------------------------------------------------------------- one centre
def _ring_one_center_111():
    """eq:ring_one_center as printed."""
    L2, L3 = mp.log(2), mp.log(3)
    return (mp.mpf(11) / 3 - mp.pi**2 / 8 - mp.mpf(3) / 2 * mp.polylog(2, -mp.mpf(1) / 2)
            - mp.mpf(46) / 9 * L2 - mp.mpf(3) / 4 * L2**2 + mp.mpf(3) / 2 * L2 * L3)


def _t_reduced(za, zb, zc):
    """eq:ring_reduced for one-centre normalised 1s densities, rho~ = 16 z^4/(k^2+4z^2)^2,
    with the w integral done in closed form."""
    za, zb, zc = map(mp.mpf, (za, zb, zc))
    rt = lambda k, z: 16 * z**4 / (k * k + 4 * z * z) ** 2
    g = 2 * zc
    Phi = lambda x: 16 * zc**4 * (x / (2 * g**2 * (x * x + g * g)) + mp.atan(x / g) / (2 * g**3))

    def inner(p):
        return mp.quad(lambda q: rt(q, zb) * (Phi(p + q) - Phi(abs(p - q))),
                       [0, p / 2, p, 2 * p + 1, mp.inf])
    return mp.quad(lambda p: rt(p, za) * inner(p), [0, mp.mpf(1) / 2, 2, 8, mp.inf]) / mp.pi


def _f000(w1, w2, w3):
    """Remiddi's master integral, in the form of Pachucki-Puchalski-Remiddi (2004) eq. 7."""
    def term(a, b, c):
        x = a / (b + c)
        return mp.log(x) * mp.log(1 + x) + mp.polylog(2, -x) + mp.polylog(2, 1 - x)
    return -(term(w3, w1, w2) + term(w2, w3, w1) + term(w1, w2, w3)) / (2 * w1 * w2 * w3)


def _t_remiddi(za, zb, zc):
    ws = tuple(2 * mp.mpf(z) for z in (za, zb, zc))
    d3 = mp.diff(lambda a, b, c: _f000(a, b, c), ws, (1, 1, 1))
    return 64 * (mp.mpf(za) * zb * zc) ** 3 * (-d3)


@pytest.mark.slow
def test_one_center_closed_form_against_both_routes():
    """Rejects a wrong eq:ring_one_center (a dropped or rescaled Li2(-1/2), pi^2 or log term)
    and a wrong eq:ring_reduced: the printed value must match the reduced momentum formula
    and, independently, Remiddi's master differentiated in the exponents."""
    mp.mp.dps = 20
    closed = _ring_one_center_111()
    assert abs(_t_remiddi(1, 1, 1) / closed - 1) < 1e-15
    assert abs(_t_reduced(1, 1, 1) / closed - 1) < 1e-12


def test_one_center_123_relation():
    """Rejects a wrong integer relation at (1,2,3) and the registered falsifier 'the triangle is
    rational': the dilogarithm relation must hold, and T must not be a small rational."""
    mp.mp.dps = 40
    T = _t_remiddi(1, 2, 3)
    L2, L3, L5 = mp.log(2), mp.log(3), mp.log(5)
    rhs = (3344 - 180 * mp.pi**2
           + 1080 * (mp.polylog(2, mp.mpf(1) / 5) + mp.polylog(2, mp.mpf(1) / 6) + mp.polylog(2, mp.mpf(1) / 3))
           - 3360 * L2 - 525 * L5 + 540 * L2**2 + 1080 * L2 * L3 + 1080 * L3**2
           - 2160 * L2 * L5 + 540 * L5**2)
    assert abs(360 * T - rhs) < mp.mpf(10) ** -30
    assert mp.pslq([T, mp.mpf(1)], maxcoeff=10**8, maxsteps=10**5) is None


# ---------------------------------------------------------------- eq:ring_op
@pytest.fixture(scope="module")
def ring_kernel():
    """L_G applied to the t-dependent factor of the momentum form, h(t,w) = w/Delta_c(w)
    with Delta_c = sqrt(w^2/t^2 + 1)."""
    h = t * w / sp.sqrt(w**2 + t**2)
    K = sp.simplify(_apply(_coeff_list(QG), h))
    return sp.lambdify((t, w), K, "mpmath")


def _LG_G(Kf, tv):
    tv = mp.mpf(tv)
    D3 = lambda k: (k * k + 1) ** mp.mpf(-1.5)
    kk = lambda x: Kf(tv, x) if x > 0 else mp.mpf(0)
    inner = lambda p: mp.quad(lambda q: D3(q) * (kk(p + q) - kk(abs(p - q))), [0, p, 2 * p + 1, mp.inf])
    return mp.quad(lambda p: D3(p) * inner(p), [0, 1, 4, mp.inf]) / mp.pi


@pytest.mark.slow
def test_ring_operator_inhomogeneous(ring_kernel):
    """Rejects a wrong eq:ring_op (any changed coefficient moves the result by O(1)) and the
    reading 'L_G G = 0': the triangle must give exactly -3 t^2."""
    mp.mp.dps = 15
    for tv in ("0.5", "2.5"):
        assert abs(_LG_G(ring_kernel, tv) + 3 * mp.mpf(tv) ** 2) < 1e-9


# ---------------------------------------------------------------- factorization
def test_ring_operator_factors_exactly():
    """Rejects 'L_G is irreducible' and a wrong algebraic factor or wrong quotient:
    (1 - t^2/4)^(-3/2) solves L_G; L_G = M o (theta - r) with zero remainder; and the map
    A = t^3 (1 + theta)/(4 - t^2) sends every Legendre solution to a solution of M."""
    LG = _coeff_list(QG)
    y = (1 - t**2 / 4) ** sp.Rational(-3, 2)
    assert _apply(LG, y) == 0
    r = sp.together(_th(y) / y)
    assert sp.simplify(r - 3 * t**2 / (4 - t**2)) == 0
    R = [-r, sp.Integer(1)]
    rem = list(LG)
    M = [0, 0, 0]
    for k in (3, 2, 1):
        c = sp.factor(sp.together(rem[k]))
        M[k - 1] = c
        sub = _compose([0] * (k - 1) + [c], R)
        rem = [sp.together(rem[i] - (sub[i] if i < len(sub) else 0)) for i in range(4)]
    assert all(sp.simplify(c) == 0 for c in rem)
    assert sp.simplify(M[2] - (-(t**2 - 1) * (t**2 - 4))) == 0
    # Legendre theta^2 - t^2 (theta+1)^2: theta^2 f = s0 f + s1 theta f
    s0, s1 = t**2 / (1 - t**2), 2 * t**2 / (1 - t**2)

    def th_pair(p):
        return (sp.together(_th(p[0]) + p[1] * s0), sp.together(p[0] + _th(p[1]) + p[1] * s1))

    a = t**3 / (4 - t**2)
    cur, acc = (a, a), [0, 0]
    for mk in M:
        acc = [acc[0] + mk * cur[0], acc[1] + mk * cur[1]]
        cur = th_pair(cur)
    assert all(sp.simplify(c) == 0 for c in acc)
    # control: the map is not trivially zero, and a wrong map (1 - theta) fails
    cur, acc = (a, -a), [0, 0]
    for mk in M:
        acc = [acc[0] + mk * cur[0], acc[1] + mk * cur[1]]
        cur = th_pair(cur)
    assert any(sp.simplify(c) != 0 for c in acc)


def test_quotient_is_the_two_density_rung():
    """Rejects identifying the quotient with the wrong elliptic family: the two-density moment
    F_2(t) = int_0^inf K0(x) K0(t x) dx must itself solve theta^2 - t^2 (theta+1)^2."""
    mp.mp.dps = 20
    K0, K1 = (lambda y: mp.besselk(0, y)), (lambda y: mp.besselk(1, y))
    mom = lambda g: mp.quad(lambda x: K0(x) * g(x), [0, 1, mp.inf])
    for tv in (mp.mpf("0.4"), mp.mpf("1.7")):
        # theta acts on K0(t x) only: theta K0(y) = -y K1(y), theta^2 K0(y) = y^2 K0(y)
        F = mom(lambda x: K0(tv * x))
        th1 = mom(lambda x: -tv * x * K1(tv * x))
        th2 = mom(lambda x: (tv * x) ** 2 * K0(tv * x))
        res = th2 - tv**2 * (th2 + 2 * th1 + F)
        assert abs(res) < 1e-14 * abs(F)


# ---------------------------------------------------------------- the chain control
# The 3-D chain at the same fibre, C(t) = (8 t^3/pi^3) int chi(x)^2 K0(t x) dx with
# chi(x) = int_0^x s K1(s) ds.  Its Mellin-route operator, sum_i t^(2i) Q_i(theta):
QC = [-36 * (u - 7) * (u - 5) ** 2 * (u - 3) ** 2 * (u - 2) * (u - 1),
      (u - 5) * (u - 3) ** 2 * (u - 2) * (45 * u**3 + 54 * u**2 - 128 * u + 164),
      -(u - 3) * (u - 2) * (9 * u**5 + 72 * u**4 + 59 * u**3 + 88 * u**2 - 8 * u - 160),
      -20 * u * (u - 2) * (u + 1) * (u + 2) * (u + 4)]
# eq:chain_k3_map on (f, theta f, theta^2 f), and the two elementary pieces as printed
K3_MAP = [t**3 * (t**4 - 18 * t**2 - 4), t**3 * (2 * t**4 - 43 * t**2 - 40),
          t * (t**6 - 27 * t**4 - 96 * t**2 + 32)]
K3_DEN = (t**2 - 4) ** 2
R1 = [-(3 + t**2), 1 - t**2]           # (theta - 3) - t^2 (theta + 1), solution t^3/(1-t^2)^2
TH2 = [sp.Integer(-2), sp.Integer(1)]  # theta - 2, solution t^2
L3 = [4 * u**3, -(u + 1) ** 3]         # F_3's operator 4 theta^3 - t^2 (theta + 1)^3


def _cth(c):
    return sp.cancel(t * sp.diff(c, t))


def _ccompose(A, B):
    out = [sp.Integer(0)] * (len(A) + len(B) - 1)
    for k, a in enumerate(A):
        for j, b in enumerate(B):
            bi = b
            for i in range(k + 1):
                out[k - i + j] += a * sp.binomial(k, i) * bi
                bi = _cth(bi)
    return [sp.cancel(c) for c in out]


def _rdiv(L, R):
    """L = Q o R + rem with ord(rem) < ord(R); exact, rational coefficients."""
    L = [sp.cancel(c) for c in L]
    r = len(R) - 1
    Q = [sp.Integer(0)] * (len(L) - r)
    for k in range(len(L) - 1, r - 1, -1):
        c = sp.cancel(L[k] / R[-1])
        Q[k - r] = c
        sub = _ccompose([0] * (k - r) + [c], R)
        L = [sp.cancel(L[i] - (sub[i] if i < len(sub) else 0)) for i in range(len(L))]
    return Q, L[:r]


def _mod_L3(L3c):
    """theta acting on (c0, c1, c2) ~ c0 f + c1 theta f + c2 theta^2 f, f in Sol(L3)."""
    S = [sp.cancel(-c / L3c[3]) for c in L3c[:3]]

    def th_vec(p):
        c0, c1, c2 = p
        return [sp.cancel(_cth(c0) + c2 * S[0]), sp.cancel(c0 + _cth(c1) + c2 * S[1]),
                sp.cancel(c1 + _cth(c2) + c2 * S[2])]
    return th_vec


def _apply_mod(L, p, th_vec):
    acc, cur = [0, 0, 0], p
    for m in L:
        acc = [sp.cancel(x + m * y) for x, y in zip(acc, cur)]
        cur = th_vec(cur)
    return acc


def _chain_pieces(qc=QC, k3_map=K3_MAP, r1=R1, l3=L3):
    """Peel the chain operator: QC = N2 o R3 o R1 o (theta - 2), R3 the annihilator of the
    image of Sol(L3) under the map.  Returns (remainders, M5 o map mod L3, R3, N2)."""
    LC = _coeff_list(qc)
    L6, rem0 = _rdiv(LC, TH2)
    M5, rem1 = _rdiv(L6, r1)
    th_vec = _mod_L3(_coeff_list(l3))
    A = [sp.cancel(a / K3_DEN) for a in k3_map]
    image = _apply_mod(M5, A, th_vec)
    vecs = [A]
    for _ in range(3):
        vecs.append(th_vec(vecs[-1]))
    ns = sp.Matrix([[vecs[k][i] for k in range(4)] for i in range(3)]).nullspace()
    R3 = [sp.cancel(x / ns[0][3]) for x in ns[0]] if len(ns) == 1 else None
    N2, rem3 = _rdiv(M5, R3) if R3 else (None, [1])
    return rem0 + rem1 + rem3, image, R3, N2


@pytest.fixture(scope="module")
def chain_pieces():
    return _chain_pieces()


@pytest.mark.slow
def test_chain_operator_peels_to_k3(chain_pieces):
    """Rejects 'the chain's K3 piece is only mod-p evidence' and any wrong map in
    eq:chain_k3_map, wrong elementary piece R1, or F_3's operator at argument t^2 instead of
    t^2/4: the chain operator divides exactly by theta - 2 and R1, and the map sends every
    solution of F_3's operator into the quotient, in exact rational arithmetic."""
    rems, image, R3, N2 = chain_pieces
    assert all(sp.simplify(c) == 0 for c in rems)
    assert all(sp.simplify(c) == 0 for c in image)
    assert len(R3) == 4 and len(N2) == 3
    # control: one changed coefficient in the map must break the identity
    bad = [K3_MAP[0], K3_MAP[1], t * (t**6 - 27 * t**4 - 95 * t**2 + 32)]
    _, image_bad, _, _ = _chain_pieces(k3_map=bad)
    assert any(sp.simplify(c) != 0 for c in image_bad)


def _to_Q(L):
    """function-left theta coefficients (rational, even in t) -> sum_i t^(2i) Q_i(u),
    after clearing denominators on the left."""
    den = sp.lcm([sp.denom(sp.together(c)) for c in L])
    Q = {}
    for k, c in enumerate(L):
        pc = sp.Poly(sp.expand(sp.cancel(c * den)), t)
        for (e,), v in zip(pc.monoms(), pc.coeffs()):
            assert e % 2 == 0
            Q[e // 2] = Q.get(e // 2, 0) + v * u**k
    return [sp.expand(Q.get(i, 0)) for i in range(max(Q) + 1)]


def _chi(x):
    x = mp.mpf(x)
    if x == 0:
        return mp.mpf(0)
    with mp.workdps(mp.mp.dps + 10):
        return +(mp.pi * x / 2 * (mp.besselk(1, x) * mp.struvel(0, x) + mp.besselk(0, x) * mp.struvel(1, x)))


def _on_chain(Q, tv):
    """(sum_i t^(2i) Q_i(theta)) C at t = tv.  theta acts on t^3 K0(t x) only:
    theta (t^3 g) = t^3 (theta + 3) g, and on g = a K0(y) + b theta K0(y) with
    theta K0(y) = -y K1(y), theta^2 K0(y) = y^2 K0(y), y = t x."""
    y = sp.Symbol("y", positive=True)
    J = max(sp.Poly(q, u).degree() for q in Q)
    al, be = [sp.Integer(1)], [sp.Integer(0)]
    for _ in range(J):
        a, b = al[-1], be[-1]
        al.append(sp.expand(y * sp.diff(a, y) + b * y**2))
        be.append(sp.expand(a + y * sp.diff(b, y)))
    parts = []
    for i, q in enumerate(Q):
        pc = sp.Poly(sp.expand(q.subs(u, u + 3)), u).all_coeffs()[::-1]
        A = sp.expand(sum(c * al[j] for j, c in enumerate(pc)))
        B = sp.expand(sum(c * be[j] for j, c in enumerate(pc)))
        parts.append((i, sp.lambdify(y, A, "mpmath"), sp.lambdify(y, B, "mpmath")))
    tv = mp.mpf(tv)

    def h(x):
        yy = tv * x
        k0, k1 = mp.besselk(0, yy), mp.besselk(1, yy)
        return _chi(x) ** 2 * sum(tv ** (2 * i) * (a(yy) * k0 - b(yy) * yy * k1) for i, a, b in parts)
    return 8 * tv**3 / mp.pi**3 * mp.quad(h, [0, mp.mpf(1) / 4, 1, 4, 16, 64, mp.inf])


@pytest.mark.slow
def test_chain_top_factor_is_spurious(chain_pieces):
    """Rejects 'the chain needs the order-two top factor N2' (or an unidentified piece
    above its K3 family): W = R3 o R1 o (theta - 2) must annihilate C itself, while W with one
    coefficient of R3 changed must not."""
    mp.mp.dps = 25
    assert abs(_chi(1) - mp.quad(lambda s: s * mp.besselk(1, s), [0, 1])) < mp.mpf(10) ** -22
    _, _, R3, _ = chain_pieces
    W = _ccompose(_ccompose(R3, R1), TH2)
    R3bad = [R3[0], R3[1], sp.cancel(R3[2] * sp.Rational(199, 198)), R3[3]]
    Wbad = _ccompose(_ccompose(R3bad, R1), TH2)
    tv = "1.3"
    C = _on_chain([sp.Integer(1)], tv)
    assert abs(_on_chain(_to_Q(W), tv)) < mp.mpf(10) ** -18 * abs(C)
    assert abs(_on_chain(_to_Q(Wbad), tv)) > mp.mpf(10) ** -3


def _f3_series(N):
    from fractions import Fraction
    out, term = {}, Fraction(1)
    for n in range(N // 2):
        out[2 * n] = term
        term = term * Fraction(2 * n + 1, 2) ** 3 / Fraction(n + 1) ** 3 / 4
    return out


def _split_nullity(Q, b, c, e, N=240, lo=8, hi=16):
    """Dimension of rational maps B = sum_j p_j(t) theta^j / den, deg p_j in [-lo, hi],
    den = (t^2-1)^b (t^2-4)^c (4t^2+11)^e, with (sum t^(2i) Q_i(theta)) (B f3) = 0 as a t-series."""
    from fractions import Fraction
    from sympy.polys.matrices import DomainMatrix
    Qp = [sp.Poly(q, u) for q in Q]
    den = sp.Poly(sp.expand((t**2 - 1) ** b * (t**2 - 4) ** c * (4 * t**2 + 11) ** e), t)
    dc = [Fraction(int(x.p), int(x.q)) for x in (sp.Rational(v) for v in den.all_coeffs()[::-1])]
    w = [1 / dc[0]]
    for n in range(1, N):
        w.append(-sum(dc[k] * w[n - k] for k in range(1, min(n, len(dc) - 1) + 1)) / dc[0])
    f = _f3_series(N)
    top = N - 2 * len(Q) - hi - 6
    cols = []
    for j in range(3):
        g = {}
        for e1, c1 in enumerate(w):
            for e2, c2 in f.items():
                if e1 + e2 < N and c1:
                    g[e1 + e2] = g.get(e1 + e2, 0) + c1 * c2 * e2**j
        for d in range(-lo, hi + 1):
            col = {}
            for e0, v in g.items():
                for i, q in enumerate(Qp):
                    val = int(q.eval(e0 + d))
                    if val and e0 + d + 2 * i < top:
                        col[e0 + d + 2 * i] = col.get(e0 + d + 2 * i, 0) + val * v
            cols.append(col)
    exps = sorted(set().union(*cols))
    rows = [[sp.Rational(col.get(x, 0).numerator, col.get(x, 0).denominator) if col.get(x, 0) else sp.Integer(0)
             for col in cols] for x in exps]
    M = DomainMatrix.from_list_sympy(len(rows), len(cols), rows).convert_to(sp.QQ)
    return M.nullspace().shape[0]


@pytest.mark.slow
def test_chain_k3_does_not_split(chain_pieces):
    """Rejects 'the chain is a rational combination of F_3 plus an elementary term' (registered
    outcomes A and B), within the searched maps: no rational B sends Sol(L3) into Sol(W).
    Positive control: the reversed composition (elementary pieces) o R3 has such a map, and the
    same search must find it."""
    _, _, R3, _ = chain_pieces
    QW = _to_Q(_ccompose(_ccompose(R3, R1), TH2))
    for e in (0, 1):
        for b in range(3):
            for c in range(3):
                assert _split_nullity(QW, b, c, e) == 0
    Qs = _to_Q(_ccompose(_ccompose(R1, TH2), R3))
    assert _split_nullity(Qs, 0, 2, 0) >= 1


# ---------------------------------------------------------------- the ring's source and its elliptic piece
def test_ring_source_is_the_operator_on_a_constant():
    """Rejects a wrong reading of the source term in eq:ring_op ('L_G (G - 1) = 0', a constant
    other than 1, or a source the operator cannot produce from a pure power): L_G applied to the
    constant 1 must be exactly 3 t^2, so with L_G G = -3 t^2 (guarded above) L_G (G + 1) = 0."""
    LG = _coeff_list(QG)
    assert sp.expand(_apply(LG, sp.Integer(1)) - 3 * t**2) == 0


R_RING = 3 * t**2 / (4 - t**2)                 # r = theta log (1 - t^2/4)^(-3/2)
# (theta - r)(G + 1) at t = 0.5 and 1.3, the nonzero solution of M that G carries; the momentum
# and Struve routes below agree on these to ~1e-16 at 15 digits
RING_ELLIPTIC = {"0.5": mp.mpf("-0.1605062101395547"), "1.3": mp.mpf("-2.328885856751796")}
_X0 = 120


def _psi(x):
    """psi = L_{-1} - I_1 (eq:ring_struve), with an asymptotic series where the two cancel."""
    x = mp.mpf(x)
    if x < _X0:
        with mp.workdps(mp.mp.dps + int(x / 2.2) + 15):
            return +(mp.struvel(-1, x) - mp.besseli(1, x))
    tot, m = mp.mpf(0), 0
    while True:
        c = (-1) ** m * mp.factorial(2 * m) * mp.binomial(mp.mpf(-1.5), m) * x ** (-2 * m - 2)
        tot += c
        if abs(c) < mp.eps * abs(tot) or m > x / 2:
            return 2 / mp.pi * tot
        m += 1


def _theta_psi(x):
    """x psi'(x), from the Struve and Bessel derivative identities."""
    x = mp.mpf(x)
    if x < _X0:
        with mp.workdps(mp.mp.dps + int(x / 2.2) + 15):
            return +(x * mp.struvel(-2, x) + mp.struvel(-1, x) - x * mp.besseli(0, x) + mp.besseli(1, x))
    tot, m = mp.mpf(0), 0
    while True:
        c = (-1) ** m * mp.factorial(2 * m) * mp.binomial(mp.mpf(-1.5), m) * x ** (-2 * m - 2) * (-2 * m - 2)
        tot += c
        if abs(c) < mp.eps * abs(tot) or m > x / 2:
            return 2 / mp.pi * tot
        m += 1


def _elliptic_piece_struve(tv):
    """(theta - r)(G + 1) on the Struve form: theta [t^2 psi(t x)] = t^2 [(2 + theta) psi](t x)."""
    rv = mp.mpf(sp.Rational(R_RING.subs(t, sp.Rational(tv))))
    tv = mp.mpf(tv)
    pts = [0, mp.mpf(1) / 4, 1, 4, 16, 64, _X0, 4 * _X0, mp.inf]
    G = mp.pi / 2 * tv**2 * mp.quad(lambda x: x * x * _psi(x) ** 2 * _psi(tv * x), pts)
    thG = mp.pi / 2 * tv**2 * mp.quad(lambda x: x * x * _psi(x) ** 2 * (2 * _psi(tv * x) + _theta_psi(tv * x)), pts)
    return thG - rv * G - rv


@pytest.mark.slow
def test_ring_needs_its_elliptic_piece():
    """Rejects 'G is algebraic plus elementary' (the elliptic piece set to zero: (theta - r)(G+1)
    = 0), a wrong algebraic factor r, and a dropped constant ((theta - r) G in place of
    (theta - r)(G + 1)).  (theta - r)(G + 1) must be -0.1605... at t = 0.5 and -2.3288... at
    t = 1.3, on the momentum form with theta taken analytically, and independently on the Struve
    form of eq:ring_struve.  Tolerance 1e-10: the two routes agree to ~1e-16 at 15 digits, while
    each named wrong answer moves the value by more than 0.1."""
    mp.mp.dps = 15
    hk = t * w / sp.sqrt(w**2 + t**2)
    kern = sp.lambdify((t, w), sp.simplify(_th(hk) - R_RING * hk), "mpmath")
    for tv, want in RING_ELLIPTIC.items():
        rv = mp.mpf(sp.Rational(R_RING.subs(t, sp.Rational(tv))))
        momentum = _LG_G(kern, tv) - rv                  # _LG_G integrates any kernel on the momentum form
        struve = _elliptic_piece_struve(tv)
        assert abs(momentum - want) < 1e-10
        assert abs(struve - want) < 1e-10
        assert abs(want) > 0.1


# ---------------------------------------------------------------- the chain: order-four factor and no split
# Work in x = t^2, delta = x d/dx (theta = 2 delta); all chain operators are even in t.
# Operators are lists [c_0, c_1, ...] meaning sum_k c_k(x) delta^k, functions on the left.
x_, lam_ = sp.symbols("x lam")
Rat = sp.Rational

# the partial section S1 (R1 o S1 = 1 mod R3), on (1, delta, delta^2)
S1_DELTA = [(12 * x_**3 + 229 * x_**2 - 31 * x_ + 12) / (3 * x_ * (x_ - 1) ** 2 * (4 * x_ + 11)),
            -2 * (20 * x_**3 + 619 * x_**2 - 143 * x_ + 80) / (15 * x_ * (x_ - 1) ** 2 * (4 * x_ + 11)),
            -4 * (x_ - 4) * (20 * x_**2 - 7 * x_ + 5) / (15 * x_ * (x_ - 1) ** 2 * (4 * x_ + 11))]
# W = Q1 o P, P = (theta - 2) - S1 o R1 o (theta - 2)
Q1_DELTA = [15 * x_ * (x_ - 1) * (x_ + 3) * (4 * x_ + 11) / ((x_ - 4) * (20 * x_**2 - 7 * x_ + 5)),
            30 * x_ * (x_ - 1) ** 2 * (4 * x_ + 11) / ((x_ - 4) * (20 * x_**2 - 7 * x_ + 5))]
# the no-split bound table: exponents of W and of R3 at each point, and the pole-order bound on
# (s0, s1, s2); 'ord' marks R3 ordinary there
NOSPLIT_TABLE = {
    sp.Integer(0): ([Rat(1, 2), 1, Rat(3, 2), Rat(3, 2), Rat(5, 2)], [Rat(1, 2), Rat(3, 2), Rat(5, 2)], [-2, -2, -2]),
    sp.Integer(1): ([-1, 0, 1, 2, 3], "ord", [-1, -1, -1]),
    sp.Integer(4): ([Rat(-3, 2), 0, 1, 2, 3], [Rat(-7, 2), 0, 1], [-2, -1, 0]),
    Rat(-11, 4): ([0, 1, 2, 3, 5], [0, 1, 3], [-3, -2, -1]),
    sp.oo: ([-1, 0, Rat(1, 2), 1, 2], [-1, 0, 1], [-2, -2, -2]),
}


def _dl(e):
    return sp.cancel(x_ * sp.diff(e, x_))


def _dcompose(A, B):
    out = [sp.Integer(0)] * (len(A) + len(B) - 1)
    for k, a in enumerate(A):
        for j, b in enumerate(B):
            bi = b
            for i in range(k + 1):
                out[k - i + j] += a * sp.binomial(k, i) * bi
                bi = _dl(bi)
    return [sp.cancel(c) for c in out]


def _ddiv(L, Rr):
    """L = Q o Rr + rem, ord rem < ord Rr (delta form, exact)."""
    L = [sp.cancel(c) for c in L]
    m = len(Rr) - 1
    Q = [sp.Integer(0)] * (len(L) - m)
    for k in range(len(L) - 1, m - 1, -1):
        c = sp.cancel(L[k] / Rr[-1])
        Q[k - m] = c
        sub = _dcompose([0] * (k - m) + [c], Rr)
        L = [sp.cancel(L[i] - (sub[i] if i < len(sub) else 0)) for i in range(len(L))]
    return Q, [sp.cancel(c) for c in L[:m]]


def _theta_to_delta(L):
    """theta-form coefficients, even in t -> delta form in x."""
    return [sp.cancel(sp.sympify(c).subs(t, sp.sqrt(x_)) * 2**k) for k, c in enumerate(L)]


def _delta_to_theta(L):
    return [sp.cancel(sp.sympify(c).subs(x_, t**2) / 2**k) for k, c in enumerate(L)]


def _chain_delta(R3):
    """R3, R1, theta - 2, E = R1 o (theta - 2), and F_3's operator L3, in delta form."""
    R3d = _theta_to_delta(R3)
    R1d = [-(3 + x_), 2 * (1 - x_)]
    T2d = [sp.Integer(-2), sp.Integer(2)]
    L3d = [sp.cancel(a + b) for a, b in zip([0, 0, 0, 32], _dcompose([-x_], [1, 6, 12, 8]))]
    return R3d, R1d, T2d, _dcompose(R1d, T2d), L3d


def _d_form(L):
    """delta form -> d/dx form."""
    out = [sp.Integer(0)] * len(L)
    dk = [sp.Integer(1)]
    for c in L:
        for j, a in enumerate(dk):
            out[j] += c * a
        nxt = [sp.Integer(0)] * (len(dk) + 1)
        for j, a in enumerate(dk):
            nxt[j] += x_ * sp.diff(a, x_)
            nxt[j + 1] += x_ * a
        dk = [sp.expand(a) for a in nxt]
    return [sp.cancel(c) for c in out]


def _exps_finite(L, q):
    """Exponents of L at a root of the irreducible q in Q[x] (q = x - p for a rational point),
    exact in Q[x]/(q); returns (exponent multiset, regular singular?, L ordinary there?)."""
    a = _d_form(L)
    n = len(a) - 1
    a = [sp.cancel(c / a[-1]) for c in a]
    qp = sp.Poly(q, x_)
    dq = qp.diff(x_).as_expr()
    I, regular, ordinary = sp.prod([lam_ - j for j in range(n)]), True, True
    for k in range(n):
        num, den = sp.fraction(sp.cancel(a[k]))
        m, dd = 0, sp.Poly(den, x_)
        while dd.degree() > 0 and sp.rem(dd, qp).is_zero:
            dd = sp.quo(dd, qp)
            m += 1
        ordinary = ordinary and m == 0
        regular = regular and m <= n - k
        if m == n - k:
            nm, dn = sp.fraction(sp.cancel(a[k] * q**m))
            v = sp.rem(sp.expand(nm * sp.invert(sp.expand(dn), q, x_)), q, x_)
            v = sp.rem(sp.expand(v * sp.invert(sp.expand(dq**m), q, x_)), q, x_)
            assert sp.Poly(v, x_).degree() <= 0
            I += v * sp.prod([lam_ - j for j in range(k)])
    rts = sp.roots(sp.Poly(sp.expand(I), lam_))
    assert sum(rts.values()) == n
    return sorted((z for z, mm in rts.items() for _ in range(mm)), key=sp.re), regular, ordinary


def _exps_inf(L):
    w_ = sp.Symbol("w_")
    tot = sp.cancel(sp.Add(*[sp.cancel(c.subs(x_, 1 / w_)) * (-lam_) ** k for k, c in enumerate(L)]))
    num, den = sp.fraction(sp.factor(tot))
    nd, dd = sp.Poly(sp.expand(num), w_), sp.Poly(sp.expand(den), w_)
    I = sp.cancel(nd.as_expr().coeff(w_, min(e[0] for e in nd.monoms()))
                  / dd.as_expr().coeff(w_, min(e[0] for e in dd.monoms())))
    rts = sp.roots(sp.Poly(I, lam_))
    assert sum(rts.values()) == len(L) - 1 == sp.Poly(I, lam_).degree()
    return sorted((z for z, mm in rts.items() for _ in range(mm)), key=sp.re)


def _sing_factors(*ops):
    fs = {x_}
    for L in ops:
        den = sp.lcm([sp.denom(sp.together(sp.cancel(c / L[-1]))) for c in L])
        for f, _ in sp.factor_list(sp.Poly(den, x_))[1]:
            if f.degree() > 0:
                fs.add(f.monic().as_expr())
    return sorted(fs, key=sp.default_sort_key)


def _split_bounds(Top, Bot):
    """Pole-order bounds on S = s0 + s1 delta + s2 delta^2 with Bot o S = 1 mod Top, from exact
    exponents: ord_p s_k >= ceil(mu_W - rho_Top + k) where Top is singular (p finite, != 0),
    ceil(mu_W) where Top is ordinary, ceil(mu_W - rho_Top) at 0 and oo (mu_W = least exponent of
    W = Top o Bot, rho_Top = largest exponent of Top)."""
    W = _dcompose(Top, Bot)
    n = len(Top) - 1
    table = {}
    for q in _sing_factors(W, Top):
        eW, regW, _ = _exps_finite(W, q)
        eT, regT, ordT = _exps_finite(Top, q)
        assert regW and regT
        mu, rho = min(sp.re(z) for z in eW), max(sp.re(z) for z in eT)
        if q == x_:
            b = [sp.ceiling(mu - rho)] * n
        elif ordT:
            b = [sp.ceiling(mu)] * n
        else:
            b = [sp.ceiling(mu - rho + k) for k in range(n)]
        p = sp.solve(q, x_)[0] if sp.Poly(q, x_).degree() == 1 else q
        table[p] = (eW, "ord" if ordT else eT, b)
    eW, eT = _exps_inf(W), _exps_inf(Top)
    table[sp.oo] = (eW, eT, [sp.ceiling(min(map(sp.re, eW)) - max(map(sp.re, eT)))] * n)
    return table


def _split_search(Top, Bot, table):
    """Exact linear algebra over Q for Bot o S = 1 mod Top within the bound.
    Returns (unknowns, rank A, rank [A|b], S or None)."""
    from sympy.polys.matrices import DomainMatrix
    n = len(Top) - 1
    dens, degs = [], []
    for k in range(n):
        D = sp.Integer(1)
        for p, (_, _, b) in table.items():
            if p is not sp.oo:
                D *= ((x_ - p) if p.is_Rational else p) ** max(0, -int(b[k]))
        dens.append(sp.expand(D))
        degs.append(sp.Poly(dens[-1], x_).degree() - int(table[sp.oo][2][k]))
    cols, labels = [], []
    for k in range(n):
        for m in range(degs[k] + 1):
            S = [sp.Integer(0)] * n
            S[k] = sp.cancel(x_**m / dens[k])
            cols.append(_ddiv(_dcompose(Bot, S), Top)[1])
            labels.append((k, m))
    target = [sp.Integer(1)] + [sp.Integer(0)] * (n - 1)
    den = sp.lcm([sp.denom(sp.together(c)) for col in cols for c in col])
    rows = []
    for i in range(n):
        polys = [sp.Poly(sp.expand(sp.cancel(col[i] * den)), x_) for col in cols]
        tp = sp.Poly(sp.expand(target[i] * den), x_)
        top = max(p.degree() for p in polys + [tp])
        for e in range(top + 1):
            rows.append([p.coeff_monomial(x_**e) for p in polys] + [tp.coeff_monomial(x_**e)])
    Mx = DomainMatrix.from_list_sympy(len(rows), len(cols) + 1, rows).convert_to(sp.QQ)
    rA, rAb = Mx[:, :len(cols)].rank(), Mx.rank()
    S = None
    if rA == rAb:
        vs, params = Mx[:, :len(cols)].to_Matrix().gauss_jordan_solve(Mx[:, len(cols):].to_Matrix())
        vs = vs.subs({p: 0 for p in params})
        S = [sp.Integer(0)] * n
        for c, (k, m) in zip(vs, labels):
            S[k] += c * x_**m / dens[k]
        S = [sp.cancel(s) for s in S]
    return len(cols), rA, rAb, S


def _annihilator(L, G):
    """Monic annihilator of G(Sol L), L of order 3."""
    n = len(L) - 1
    vecs = [_ddiv(G, L)[1]]
    for _ in range(n):
        vecs.append(_ddiv(_dcompose([0, 1], vecs[-1]), L)[1])
    ns = sp.Matrix([[vecs[k][i] for k in range(n + 1)] for i in range(n)]).nullspace()
    assert len(ns) == 1
    return [sp.cancel(c / ns[0][n]) for c in ns[0]]


def _is_one_mod(Bot, S, Top):
    rem = _ddiv(_dcompose(Bot, S), Top)[1]
    return [sp.cancel(c) for c in rem] == [1] + [0] * (len(Top) - 2)


@pytest.mark.slow
def test_chain_partial_section_and_order_four_factor(chain_pieces):
    """Rejects a wrong partial section S1 (any changed coefficient) and a W that is not a left
    multiple of the order-four P: R1 o S1 = 1 mod R3 and W = Q1 o P with zero remainder, both in
    exact rational arithmetic, with R3 the K3 piece peeled from the chain operator above.
    Control: S1 with one coefficient changed must leave a nonzero remainder."""
    _, _, R3, _ = chain_pieces
    R3d, R1d, T2d, Ed, _ = _chain_delta(R3)
    assert _is_one_mod(R1d, S1_DELTA, R3d)
    assert not _is_one_mod(R1d, [S1_DELTA[0] * Rat(101, 100)] + S1_DELTA[1:], R3d)
    P = [sp.cancel(a - b) for a, b in zip(T2d + [0, 0, 0], _dcompose(S1_DELTA, Ed))]
    assert len(P) == 5 and sp.cancel(P[4]) != 0
    Q, rem = _ddiv(_dcompose(R3d, Ed), P)
    assert all(c == 0 for c in rem)
    assert [sp.cancel(a - b) for a, b in zip(Q, Q1_DELTA)] == [0, 0]
    assert sp.cancel(sum(P) * x_) == 0                      # P kills t^2 = x (delta^k x = x)


@pytest.mark.slow
def test_chain_k3_piece_does_not_split_proof(chain_pieces):
    """Rejects 'a split exists' and a shrunken bound: the exponents of W = R3 o E and of R3 at
    x = 0, 1, 4, -11/4, oo must be the recorded table (so a miscomputed exponent, which would move
    a bound, is caught), the bound must give 27 unknowns, and the exact system E o S = 1 mod R3 must
    be inconsistent (rank 27 < 28).  Positive control, the same code: on the planted split
    Top = ann(E(Sol L3)) it must find S, and S must verify exactly."""
    _, _, R3, _ = chain_pieces
    R3d, _, _, Ed, L3d = _chain_delta(R3)
    Top = R3d
    table = _split_bounds(Top, Ed)
    assert set(table) == set(NOSPLIT_TABLE)
    for p, (eW, eT, b) in NOSPLIT_TABLE.items():
        gW, gT, gb = table[p]
        assert list(gW) == [sp.sympify(z) for z in eW]
        assert gT == eT if eT == "ord" else list(gT) == [sp.sympify(z) for z in eT]
        assert list(gb) == b
    nunk, rA, rAb, S = _split_search(Top, Ed, table)
    assert (nunk, rA, rAb, S) == (27, 27, 28, None)
    Topc = _annihilator(L3d, Ed)
    tc = _split_bounds(Topc, Ed)
    nc, rc, rcb, Sc = _split_search(Topc, Ed, tc)
    assert nc == 34 and rc == rcb == 34 and Sc is not None
    assert _is_one_mod(Ed, Sc, Topc)


def _on_chain_delta(Pd, tv):
    """(P C)(tv) for a delta-form P with rational coefficients: clear the left denominator the
    way _to_Q does, evaluate on the chain integral, divide the denominator back out."""
    Pt = _delta_to_theta(Pd)
    den = sp.lcm([sp.denom(sp.together(c)) for c in Pt])
    return _on_chain(_to_Q(Pt), tv) / mp.mpf(sp.Rational(den.subs(t, sp.Rational(tv))))


@pytest.mark.slow
def test_chain_lies_on_the_order_four_factor(chain_pieces):
    """Rejects 'the chain needs the artanh solution of R1 (theta - 2)' (minimal order five):
    P C / h, h = t^3/(1-t^2)^2, must vanish to 1e-25 at 30 digits at two t, while P built from
    S1 with its first coefficient times 1.01 must give more than 1e-3 at both."""
    mp.mp.dps = 30
    _, _, R3, _ = chain_pieces
    R3d, _, T2d, Ed, _ = _chain_delta(R3)

    def P_of(S):
        return [sp.cancel(a - b) for a, b in zip(T2d + [0, 0, 0], _dcompose(S, Ed))]
    P = P_of(S1_DELTA)
    Pbad = P_of([S1_DELTA[0] * Rat(101, 100)] + S1_DELTA[1:])
    for tv in ("0.5", "1.3"):
        h = mp.mpf(tv) ** 3 / (1 - mp.mpf(tv) ** 2) ** 2
        assert abs(_on_chain_delta(P, tv) / h) < mp.mpf(10) ** -25
        assert abs(_on_chain_delta(Pbad, tv) / h) > mp.mpf(10) ** -3


# ---------------------------------------------------------------- the correlated triangle
# Y(t) = <rho_i rho_j rho_k f_ik f_jk / r_ij>, f = exp(-gamma r), gamma = 1, scales (1, 1, t),
# t on the density carrying both geminals.  Record docs/predictions/2026-10-17_ring_extensions.md
# (QG; committed before the computation, which falsified its prediction QG2).
# The minimal annihilator, sum_i t^(2i) Q_i(theta):
QG4 = [-200960 * (u - 5) * (u - 3) ** 2 * (u - 1),
       64 * (u - 3) * (1648 * u**3 + 6992 * u**2 - 9043 * u - 1010),
       -16 * (932 * u**4 + 12704 * u**3 + 54000 * u**2 + 71790 * u - 45),
       12 * u * (20 * u**3 + 376 * u**2 + 1637 * u + 1875),
       9 * u * (u + 2) * (u + 4) * (u + 6)]
# its order-one right factor theta - r, with the algebraic solution
# t^3 (33 t^6 - 612 t^4 + 5232 t^2 - 49088) / (t^2 - 4)^(9/2)
R_QG = 12 * (3 * t**6 - 316 * t**4 + 15824 * t**2 + 49088) / (
    (t**2 - 4) * (33 * t**6 - 612 * t**4 + 5232 * t**2 - 49088))
# the map on (f, theta f, theta^2 f), f a solution of F_3's operator, up to a constant
_QG_APP = (9 * t**4 + 312 * t**2 - 12560) ** 2
QG_MAP = [(135 * t**10 - 9648 * t**8 + 389792 * t**6 - 2708448 * t**4 - 34259456 * t**2 - 5890560) / t**3,
          2 * (135 * t**10 - 8298 * t**8 + 401336 * t**6 - 2935648 * t**4 - 41430400 * t**2 + 45388800) / t**3,
          (t**2 - 4) * (135 * t**10 - 6678 * t**8 + 422744 * t**6 - 1960864 * t**4 - 67384192 * t**2
                        - 11781120) / t**5]


def _adjQ(Q):
    """(t^(2i) Q_i(theta))* = t^(2i) Q_i(-theta - 2i), for the pairing with dt/t."""
    return [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(Q)]


def _qg_peel(qg4=QG4, qg_map=QG_MAP, r_qg=R_QG):
    th_vec = _mod_L3(_coeff_list(L3))
    A = [sp.cancel(a / _QG_APP) for a in qg_map]
    image = _apply_mod(_coeff_list(_adjQ(qg4)), A, th_vec)
    vecs = [A, th_vec(A), th_vec(th_vec(A))]
    det = sp.Matrix(vecs).det()
    vecs.append(th_vec(vecs[-1]))
    ns = sp.Matrix([[vecs[j][i] for j in range(4)] for i in range(3)]).nullspace()
    R = [sp.cancel(x / ns[0][3]) for x in ns[0]]
    Q1, rem_adj = _rdiv(_coeff_list(_adjQ(qg4)), R)
    _, rem_r1 = _rdiv(_coeff_list(qg4), [-r_qg, sp.Integer(1)])
    return image, det, Q1, rem_adj, rem_r1


@pytest.mark.slow
def test_correlated_triangle_operator_peels_to_k3():
    """Rejects 'the correlated triangle's annihilator has no K3 piece' (registered QG2, elliptic or
    lower) and any wrong map or wrong algebraic factor: the order-four annihilator divides exactly
    by theta - r, whose solution is algebraic; its adjoint has a right factor R whose solutions are
    the image of the solutions of F_3's operator under QG_MAP (injective); F_3's operator is
    anti-self-adjoint, so the order-three left factor of the annihilator, the adjoint of R, is
    equivalent to it.  All in exact rational arithmetic."""
    image, det, Q1, rem_adj, rem_r1 = _qg_peel()
    assert all(sp.simplify(c) == 0 for c in image)
    assert sp.simplify(det) != 0
    assert len(Q1) == 2 and all(sp.simplify(c) == 0 for c in rem_adj)
    assert all(sp.simplify(c) == 0 for c in rem_r1)
    assert all(sp.expand(a + b) == 0 for a, b in zip(_adjQ(L3), L3))
    y = t**3 * (33 * t**6 - 612 * t**4 + 5232 * t**2 - 49088) / (t**2 - 4) ** sp.Rational(9, 2)
    assert sp.simplify(t * sp.diff(y, t) - R_QG * y) == 0
    # control: a 0.1% change in one entry of the map must break the identity
    bad = [QG_MAP[0], QG_MAP[1], QG_MAP[2] * sp.Rational(1001, 1000)]
    assert any(sp.simplify(c) != 0 for c in _qg_peel(qg_map=bad)[0])


def _kite(wv, g, n=120):
    """Phi~(w) of the geminal kite: the two densities rho_i, rho_j and the two geminals around the
    massive loop, Fourier transform in the apex momentum w.  Double-bubble form: Feynman
    parameters u, v on the {rho_i, rho_j} and {f, f} bubbles, the Coulomb line in position space,
    Phi~ = 2 g^2 int int sqrt(u(1-u)) v(1-v) K / (M^3 N^5) du dv, with
    M = sqrt(1 + u(1-u) w^2), N = sqrt(g^2 + v(1-v) w^2), s = M + N, k2 = (u-v)^2 w^2, D = s^2 + k2,
    K = 3/D + 6 s^2/D^2 + 2(N^2 + 3MN)(3s^2 - k2)/D^3 + 24 M N^2 s (s^2 - k2)/D^4."""
    import numpy as np
    k = np.arange(1, n + 1)
    a = np.pi * k / (n + 1)
    uu, Wu = (1 - np.cos(a)) / 2, np.pi / (n + 1) * np.sin(a) ** 2 / 4  # Gauss-Chebyshev, 2nd kind
    x, wx = np.polynomial.legendre.leggauss(n)
    vv, Wv = (1 - x) / 2, wx / 2
    z = (np.asarray(wv, float) ** 2)[:, None, None]
    uu, vv = uu[None, :, None], vv[None, None, :]
    M = np.sqrt(1 + uu * (1 - uu) * z)
    V = vv * (1 - vv)
    N = np.sqrt(g * g + V * z)
    s, k2 = M + N, (uu - vv) ** 2 * z
    s2, iD = s * s, 1 / (s * s + k2)
    K = (3 * iD + 6 * s2 * iD**2 + 2 * (N * N + 3 * M * N) * (3 * s2 - k2) * iD**3
         + 24 * M * N * N * s * (s2 - k2) * iD**4)
    return 2 * g * g * np.einsum("aij,i,j->a", K * V / N**5 / M**3, Wu, Wv)


def _theta_Y(g, ts, J=4, nw=400):
    """theta^j Y(t), j = 0..J, Y = (1/(2 pi^2)) int_0^inf w^2 t^3 (w^2 + t^2)^(-3/2) Phi~(w) dw;
    theta acts on the kernel (1 + x)^(-3/2), x = w^2/t^2, as -2x d/dx."""
    import numpy as np
    xs = sp.Symbol("x", positive=True)
    ker = [(1 + xs) ** sp.Rational(-3, 2)]
    for _ in range(J):
        ker.append(sp.simplify(-2 * xs * sp.diff(ker[-1], xs)))
    kf = [sp.lambdify(xs, c, "numpy") for c in ker]
    xg, wg = np.polynomial.legendre.leggauss(nw)
    sg = (xg + 1) / 2
    wv = np.tan(np.pi * sg / 2)
    dw = wg / 2 * np.pi / 2 / np.cos(np.pi * sg / 2) ** 2
    base = wv**2 * _kite(wv, g) * dw / (2 * np.pi**2)
    return {tv: [float(np.sum(base * f(wv**2 / tv**2))) for f in kf] for tv in ts}


def _apply_Q_num(Q, thY, tv):
    terms = [tv ** (2 * i) * sum(float(c) * thY[j] for j, c in enumerate(sp.Poly(q, u).all_coeffs()[::-1]))
             for i, q in enumerate(Q)]
    return abs(sum(terms)) / max(abs(x) for x in terms)


@pytest.mark.slow
def test_correlated_triangle_sits_on_the_k3_rung():
    """Rejects 'the correlated triangle sits on the elliptic rung or lower' (registered QG2) and
    'Y is algebraic': on Y itself, from an independent evaluation of the kite, the order-four
    annihilator vanishes to 1e-10 of its largest term at three t, a 1e-4 change in one coefficient
    does not, and (theta - r) Y is of order Y.  So (theta - r) Y is a nonzero solution of the
    order-three factor that test_correlated_triangle_operator_peels_to_k3 ties to F_3's operator,
    which is irreducible: Y has a K3 constituent.  Normalisation: Phi~(0) = 3 pi/16 exactly
    (position-space closed form), and Y(0.8) against the 30-digit grid of the record."""
    import numpy as np
    g = 1.0
    assert abs(_kite([0.0], g)[0] - 3 * np.pi / 16) < 1e-13
    ts = (0.5, 0.8, 1.3)
    thY = _theta_Y(g, ts)
    assert abs(thY[0.8][0] / 0.003981259454856025324 - 1) < 1e-9
    bad = [QG4[0], QG4[1], QG4[2] * sp.Rational(10001, 10000), QG4[3], QG4[4]]
    r = sp.lambdify(t, R_QG)
    for tv in ts:
        assert _apply_Q_num(QG4, thY[tv], tv) < 1e-10
        assert _apply_Q_num(bad, thY[tv], tv) > 1e-6
        assert abs(thY[tv][1] - r(tv) * thY[tv][0]) > 0.5 * abs(thY[tv][0])


# ---------------------------------------------------------------- no split (SMALL sweep 2026-10-18)
def _even_rdiv_delta(Q, r):
    """L = Top o (theta - r) exactly, both in delta form."""
    Qq, rem = _rdiv(_coeff_list(Q), [-r, sp.Integer(1)])
    assert all(sp.simplify(c) == 0 for c in rem)
    return _theta_to_delta(Qq), _theta_to_delta([-r, sp.Integer(1)])


def _annihilator_any(L, G):
    """Monic annihilator of G(Sol L) for any order of G (pads the remainders)."""
    n = len(L) - 1
    pad = lambda v: list(v) + [sp.Integer(0)] * (n - len(v))
    vecs = [pad(_ddiv(G, L)[1])]
    for _ in range(n):
        vecs.append(pad(_ddiv(_dcompose([0, 1], vecs[-1]), L)[1]))
    ns = sp.Matrix([[vecs[k][i] for k in range(n + 1)] for i in range(n)]).nullspace()
    assert len(ns) == 1
    return [sp.cancel(c / ns[0][n]) for c in ns[0]]


@pytest.mark.slow
def test_ring_and_correlated_triangle_operators_do_not_split():
    """Rejects 'the factorisation is a direct sum', which would let G be algebraic plus a rational image
    of an elliptic period, and would let an order-three operator annihilate the correlated triangle Y.
    L_G = M o (theta - lambda_G) and the order-four annihilator of Y, R3~ o R1~: the exact system
    Bot o S = 1 mod Top, within the pole and degree bound set by the local exponents, must be
    inconsistent (7 unknowns, rank 7 < 8; 48 unknowns, rank 48 < 49).  Since R3~ is irreducible
    (equivalent to F_3's operator), no split means no order-three right factor of R3~ o R1~ holds Y.
    Positive control, the same code: on the planted split Top' = ann(Bot(Sol Top)) the search must
    find S (16 and 87 unknowns) and S must verify exactly."""
    for Q, r, want, ctrl in ((QG, R_RING, (7, 7, 8), 16), (QG4, R_QG, (48, 48, 49), 87)):
        Top, Bot = _even_rdiv_delta(Q, r)
        nunk, rA, rAb, S = _split_search(Top, Bot, _split_bounds(Top, Bot))
        assert (nunk, rA, rAb) == want and S is None
        Topc = _annihilator_any(Top, Bot)
        nc, rc, rcb, Sc = _split_search(Topc, Bot, _split_bounds(Topc, Bot))
        assert nc == rc == rcb == ctrl and Sc is not None
        assert _is_one_mod(Bot, Sc, Topc)


# ---------------------------------------------------------------- the massive loop (SMALL sweep 2026-10-18)
def _loop_closed(x1, x2, m, nu):
    """Solid-angle closed form: I_111 = pi^2 arctan(sqrt(det G)/D)/sqrt(det G) in the Gram data
    G_ij = (m_i^2 + m_j^2 + s_ij)/2, D = m1 m2 m3 + G12 m3 + G13 m2 + G23 m1; the geminal loop
    (propagators 1 and 2 squared) is d/d(m1^2) d/d(m2^2) of it."""
    x1 = [mp.mpf(v) for v in x1]
    x2 = [mp.mpf(v) for v in x2]
    s12 = sum((a - b)**2 for a, b in zip(x1, x2))
    s13 = sum(a * a for a in x1)
    s23 = sum(b * b for b in x2)

    def i111(a, b, c):
        m1, m2, m3 = mp.sqrt(a), mp.sqrt(b), mp.sqrt(c)
        g12, g13, g23 = (a + b + s12) / 2, (a + c + s13) / 2, (b + c + s23) / 2
        det = a * (b * c - g23**2) - g12 * (g12 * c - g23 * g13) + g13 * (g12 * g23 - b * g13)
        D = m1 * m2 * m3 + g12 * m3 + g13 * m2 + g23 * m1
        return mp.pi**2 * mp.atan(mp.sqrt(det) / D) / mp.sqrt(det)

    a, b, c = [mp.mpf(v)**2 for v in m]
    if nu == (1, 1, 1):
        return i111(a, b, c)
    return mp.diff(lambda aa, bb: i111(aa, bb, c), (a, b), (1, 1))


def _loop_direct(x1, x2, m, nu, Nang=120, Nr=160):
    """int d^3 l prod_i (|l - x_i|^2 + m_i^2)^(-nu_i) with x3 = 0: spherical coordinates, Gauss-Legendre
    in cos(theta) and r (three panels plus the tail r = R/s), trapezoid in phi; float64."""
    x1, x2 = np.array(x1), np.array(x2)
    uu, wu = np.polynomial.legendre.leggauss(Nang)
    ph = np.linspace(0, 2 * np.pi, Nang, endpoint=False)
    wph = 2 * np.pi / Nang
    st = np.sqrt(1 - uu**2)
    n = np.stack([np.outer(st, np.cos(ph)), np.outer(st, np.sin(ph)), np.outer(uu, np.ones_like(ph))], -1)
    d1, d2 = n @ x1, n @ x2
    R = max(np.linalg.norm(x1), np.linalg.norm(x2)) + 2
    sr, wr = np.polynomial.legendre.leggauss(Nr)

    def f(ri):
        return ((ri**2 - 2 * ri * d1 + x1 @ x1 + m[0]**2) ** -nu[0]
                * (ri**2 - 2 * ri * d2 + x2 @ x2 + m[1]**2) ** -nu[1] * (ri**2 + m[2]**2) ** -nu[2] * ri**2)

    tot = 0.0
    for lo, hi in [(0, R / 3), (R / 3, 2 * R / 3), (2 * R / 3, R)]:
        for ri, wi in zip((hi - lo) / 2 * sr + (hi + lo) / 2, (hi - lo) / 2 * wr):
            tot += wi * np.sum(wu[:, None] * f(ri)) * wph
    for si, wi in zip((sr + 1) / 2, wr / 2):
        tot += wi * R / si**2 * np.sum(wu[:, None] * f(R / si)) * wph
    return tot


def test_massive_loop_is_a_solid_angle():
    """Rejects a wrong closed form for the correlated triangle's loop: with masses (gamma, gamma, 0),
    gamma = 1, the solid-angle formula and its (2,2,1) mass derivative (the geminal loop) must match
    direct three-dimensional quadrature to 1e-8 at two external triples.  Controls: the massless
    uniqueness value pi^3/(|x12||x13||x23|) applied to the massive loop, and the formula with
    sqrt(det G) replaced by det G, must each miss the quadrature by more than 1e-3."""
    mp.mp.dps = 30
    triples = [([0.3, -1.1, 0.7], [1.2, 0.4, -0.5]), ([-0.9, 0.2, 1.4], [0.5, -1.3, -0.8])]
    m = (1.0, 1.0, 0.0)
    for x1, x2 in triples:
        for nu in ((1, 1, 1), (2, 2, 1)):
            closed = float(_loop_closed(x1, x2, m, nu))
            direct = _loop_direct(x1, x2, m, nu)
            assert abs(direct / closed - 1) < 1e-8
    x1, x2 = triples[0]
    good = _loop_direct(x1, x2, m, (1, 1, 1))
    s12 = sum((a - b)**2 for a, b in zip(x1, x2))
    s13 = sum(a * a for a in x1)
    s23 = sum(b * b for b in x2)
    g12, g13, g23 = (2 + s12) / 2, (1 + s13) / 2, (1 + s23) / 2
    det = g12 * g23 * g13 * 2 - g23**2 - g13**2
    D = g13 + g23
    unique = float(mp.pi**3 / mp.sqrt(s12 * s13 * s23))
    det_for_sqrt = float(mp.pi**2 * mp.atan(det / D) / det)
    for wrong in (unique, det_for_sqrt):
        assert abs(good / wrong - 1) > 1e-3
