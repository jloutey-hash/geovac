"""Paper 59 sec:four_density -- four linked densities give the first rung beyond elliptic.

Guards (each names the wrong answer it rejects):
  * tree identity eq:tree at n = 4, star and chain: rejects a shape-dependent or mis-normalised
    reduction (a 2^(n-1) prefactor, or star != chain);
  * the Mellin-route ODE for F_n(t) = int K0(x)^(n-1) K0(tx) dx: rejects a wrong operator, both by
    symbolic equality with eq:four_density_op and by direct quadrature of F_4;
  * the quadratic-relation counts 0, 1, 0, 1 for n = 2..5, with instrument controls that must see
    1 (product of two elliptic periods) and 3 (elliptic cube): rejects F_4 being elliptic-built;
  * irreducibility (no order <= 3 annihilator of the holomorphic solution, degree <= 12, rank mod a
    prime) with the positive control that order 4 finds the operator: rejects a reducible L;
  * the Yukawa coupling is not constant (16, 1152, 99520): rejects a symmetric cube.
  * irreducibility for n = 2..5 with no degree bound, by local exponents and Fuchs' relation, with
    the reducible triangle operator L_G as positive control: rejects a factor of any degree;
  * the banana relation in y = 1/(16 t^2): rejects a wrong variable, Hadamard factor, propagator
    count, or AESZ operator.
  * the Ising class at t = 1, n <= 8: rejects a wrong threshold set or parity, a mis-read local
    type at a threshold, a wrong operator at n = 6, 7, 8, and a wrong closed form or L-value for
    c_{3,0}, c_{4,0}, c_{6,0}.
Records: docs/predictions/2026-10-16_linking_ceiling.md; docs/predictions/2026-10-17_ising_class.md."""
import itertools
import pytest
from fractions import Fraction as Fr
import numpy as np
import mpmath as mp
import sympy as sp
from sympy.functions.combinatorial.numbers import stirling

z, s, u, t, h, th = sp.symbols("z s u t h th")


# ---------------------------------------------------------------- tree identity
def _tree3d(c, forms, N=201, U=9.0):
    g = np.linspace(-U, U, N)
    dh = g[1] - g[0]
    k = np.sinh(g)
    w = np.cosh(g) * dh
    P, Q, R = np.meshgrid(k, k, k, indexing="ij", sparse=True)
    tot = w[:, None, None] * w[None, :, None] * w[None, None, :]
    for cv, L in zip(c, forms(P, Q, R)):
        tot = tot * (cv * L**2 + 1.0) ** -1.5
    return float(tot.sum())


def _bessel_side(c):
    mp.mp.dps = 20
    a = [1 / mp.sqrt(cv) for cv in c]
    f = lambda y: y**4 * mp.fprod(mp.besselk(1, av * y) for av in a)
    return float(2**4 / (mp.pi * mp.fprod(c)) * mp.quad(f, [0, 1, mp.inf]))


def test_tree_identity_star_and_chain_n4():
    star = lambda p, q, r: (p + q + r, p, q, r)
    chain = lambda p, q, r: (p, q - p, r - q, r)
    for c in [(0.2, 0.1, 0.15, 0.05), (0.25, 0.08, 0.12, 0.2)]:
        ref = _bessel_side(c)
        for forms in (star, chain):
            assert abs(_tree3d(c, forms) / ref - 1) < 1e-5


# ---------------------------------------------------------------- Mellin-route ODE
def _sym_power_ode(m):
    def theta(vec):
        new = [sp.Integer(0)] * (m + 1)
        for i, cf in enumerate(vec):
            new[i] += z * 2 * sp.diff(cf, z)
            if i < m:
                new[i + 1] += cf * (m - i)
            if i > 0:
                new[i - 1] += cf * i * z
        return [sp.expand(c) for c in new]
    rows = [[sp.Integer(1)] + [sp.Integer(0)] * m]
    for _ in range(m + 1):
        rows.append(theta(rows[-1]))
    ns = sp.Matrix(rows).T.nullspace()
    assert len(ns) == 1
    v = ns[0]
    den = sp.lcm([sp.fraction(sp.together(e))[1] for e in v])
    v = [sp.factor(sp.cancel(e * den)) for e in v]
    g = sp.gcd_list([e for e in v if e != 0])
    return [sp.expand(sp.cancel(e / g)) for e in v]


def _t_ode(n):
    """sum_i t^{2i} Q_i(theta_t) F_n = 0, from the Mellin recurrence of K0^(n-1)."""
    v = _sym_power_ode(n - 1)
    K = max(sp.Poly(e, z).degree() for e in v if e != 0)
    P = [sp.expand(sum(sp.Poly(v[j], z).coeff_monomial(z**k) * u**j for j in range(len(v))))
         for k in range(K + 1)]
    Rt = [sp.expand(P[k].subs(u, s - 1 - 2 * k) * sp.prod([(s - 2 * i)**2 for i in range(1, k + 1)]))
          for k in range(K + 1)]
    R = [sp.expand(Rt[K - i].subs(s, s + 2 * K)) for i in range(K + 1)]
    return [sp.expand(R[i].subs(s, -u - 2 * i)) for i in range(K + 1)]


def test_four_density_operator_matches_eq_four_density_op():
    Q = _t_ode(4)
    Lw = [q.subs(u, 2 * th) for q in Q]                       # theta_t = 2 theta_w, w = t^2
    Lv = {2 - i: sp.expand(c.subs(th, -th)) for i, c in enumerate(Lw)}   # v = 1/w
    Lx = {k: sp.expand(c.subs(th, th + sp.Rational(1, 2)) * 4**k) for k, c in Lv.items()}
    lead = Lx[0].as_poly(th).LC()
    target = {0: th**4, 1: -(2 * th + 1)**2 * (10 * th**2 + 10 * th + 3),
              2: 9 * (2 * th + 1)**2 * (2 * th + 3)**2}
    for k in target:
        assert sp.expand(Lx[k] / lead - target[k]) == 0


def _annihilation_residual(n, T, power=None):
    """|L_n F| / (sum of the sizes of its terms) at t = T, with L_n = _t_ode(n) and
    F = int K0(x)^power K0(t x) dx by quadrature (power = n - 1 gives F_n).  theta^j acts on
    K0(t x) through K0 and x K1."""
    Q = _t_ode(n)
    power = n - 1 if power is None else power
    mp.mp.dps = 20
    J = max(sp.Poly(q, u).degree() for q in Q)
    al, be = [sp.Integer(1)], [sp.Integer(0)]                # theta^j K0 = al K0 + be (u K1)
    for _ in range(J):
        a_, b_ = al[-1], be[-1]
        al.append(sp.expand(u * sp.diff(a_, u) - u**2 * b_))
        be.append(sp.expand(u * sp.diff(b_, u) - a_))
    alf = [sp.lambdify(u, a_, "mpmath") for a_ in al]
    bef = [sp.lambdify(u, b_, "mpmath") for b_ in be]
    T = mp.mpf(T)
    thF = [mp.quad(lambda x, j=j: mp.besselk(0, x)**power * (alf[j](T * x) * mp.besselk(0, T * x)
                   + bef[j](T * x) * T * x * mp.besselk(1, T * x)), [0, 1, 4, mp.inf])
           for j in range(J + 1)]
    tot = mp.mpf(0)
    scale = mp.mpf(0)
    for i, q in enumerate(Q):
        for (deg,), cf in sp.Poly(q, u).terms():
            term = T**(2 * i) * mp.mpf(int(cf)) * thF[deg]
            tot += term
            scale += abs(term)
    return abs(tot) / scale


def test_four_density_operator_annihilates_F4_by_quadrature():
    assert _annihilation_residual(4, "0.7") < 1e-15


# ---------------------------------------------------------------- quadratic relations
def _taylor_basis(Q, t0, N):
    r = max(sp.Poly(q, u).degree() for q in Q)
    a = [sp.Integer(0)] * (r + 1)
    for i, q in enumerate(Q):
        for (j,), cf in sp.Poly(q, u).terms():
            for k in range(j + 1):
                a[k] += cf * t**(2 * i) * stirling(j, k, kind=2) * t**k
    A = [[Fr(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in
          reversed(sp.Poly(sp.expand(x.subs(t, t0 + h)), h).all_coeffs())] for x in a]
    ff = lambda n, k: 1 if k == 0 else n * ff(n - 1, k - 1)
    out = []
    for e in range(r):
        c = [Fr(0)] * (N + r)
        for k in range(r):
            c[k] = Fr(1 if k == e else 0) / ff(k, k)
        for M in range(N):
            tot = Fr(0)
            for k in range(r + 1):
                for l, alk in enumerate(A[k]):
                    idx = M - l + k
                    if alk == 0 or idx < 0 or (k == r and l == 0):
                        continue
                    tot += alk * c[idx] * ff(idx, k)
            c[M + r] = -tot / (A[r][0] * ff(M + r, r))
        out.append(c[:N])
    return out


def _mul(a, b):
    return [sum(a[i] * b[k - i] for i in range(k + 1)) for k in range(len(a))]


def _rank(rows):
    rows = [r[:] for r in rows]
    rk = 0
    for col in range(len(rows[0])):
        piv = next((i for i in range(rk, len(rows)) if rows[i][col] != 0), None)
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][col] != 0:
                f = rows[i][col] / rows[rk][col]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def _relations(basis):
    r = len(basis)
    prods = [_mul(basis[i], basis[j]) for i in range(r) for j in range(i, r)]
    return len(prods) - _rank(prods)


def test_quadric_counter_controls_see_one_and_three():
    N, t0 = 50, sp.Rational(1, 3)
    Q2 = _t_ode(2)
    u1, u2 = _taylor_basis(Q2, t0, N)                             # elliptic family in t
    cube = [_mul(_mul(u1, u1), u1), _mul(_mul(u1, u1), u2), _mul(_mul(u1, u2), u2), _mul(_mul(u2, u2), u2)]
    assert _relations(cube) == 3
    Q2b = [q * 4**i for i, q in enumerate(Q2)]                    # the same family in 2t: a different equation
    v1, v2 = _taylor_basis(Q2b, t0, N)
    tensor = [_mul(a, b) for a in (u1, u2) for b in (v1, v2)]
    assert _relations(tensor) == 1
    assert _relations([u1, u2]) == 0


def test_quadric_counts_0_1_0_1():
    counts = [_relations(_taylor_basis(_t_ode(n), sp.Rational(1, 2), 50)) for n in (2, 3, 4, 5)]
    assert counts == [0, 1, 0, 1]


# ---------------------------------------------------------------- irreducibility + Yukawa
P0 = lambda s_: s_**4
P1 = lambda s_: 4 * (2 * s_ + 1)**2 * (10 * s_**2 + 10 * s_ + 3)
P2 = lambda s_: 144 * (2 * s_ + 1)**2 * (2 * s_ + 3)**2


def _B(NT):
    B = [1]
    for n in range(1, NT):
        r = P1(n - 1) * B[n - 1] - (P2(n - 2) * B[n - 2] if n >= 2 else 0)
        assert r % P0(n) == 0
        B.append(r // P0(n))
    return B


def test_holomorphic_solution_is_the_hadamard_product():
    B = _B(25)
    cf = lambda n: sp.binomial(2 * n, n)**2 * sum(sp.binomial(n, k)**2 * sp.binomial(2 * k, k)
                                                  for k in range(n + 1))
    assert all(B[n] == cf(n) for n in range(25))


def test_four_density_operator_is_irreducible():
    PR = (1 << 61) - 1
    NT = 160
    Bp = [b % PR for b in _B(NT)]

    def rank_mod(rows):
        rows = [x[:] for x in rows]
        rk = 0
        for c in range(len(rows[0])):
            piv = next((i for i in range(rk, len(rows)) if rows[i][c]), None)
            if piv is None:
                continue
            rows[rk], rows[piv] = rows[piv], rows[rk]
            iv = pow(rows[rk][c], PR - 2, PR)
            rows[rk] = [x * iv % PR for x in rows[rk]]
            for i in range(len(rows)):
                if i != rk and rows[i][c]:
                    f = rows[i][c]
                    rows[i] = [(x - f * y) % PR for x, y in zip(rows[i], rows[rk])]
            rk += 1
        return rk

    def annihilators(order, D):
        rows = [[pow(m - d, j, PR) * Bp[m - d] % PR if m - d >= 0 else 0
                 for j in range(order + 1) for d in range(D + 1)] for m in range(NT)]
        return (order + 1) * (D + 1) - rank_mod(rows)

    assert [annihilators(o, 12) for o in (1, 2, 3)] == [0, 0, 0]
    assert annihilators(4, 2) == 1


def test_yukawa_coupling_is_not_constant():
    N, E = 8, 4

    def emul(p, q):
        return [sum(p[i] * q[k - i] for i in range(k + 1)) for k in range(E)]

    def einv(p):
        r_ = [Fr(1) / p[0]] + [Fr(0)] * (E - 1)
        for k in range(1, E):
            r_[k] = -sum(p[i] * r_[k - i] for i in range(1, k + 1)) / p[0]
        return r_

    def ep(f, n):
        x = sp.Symbol("x")
        Pl = sp.Poly(sp.expand(f(x + n)), x)
        return [Fr(int(Pl.coeff_monomial(x**k))) for k in range(E)]

    A = [[Fr(1)] + [Fr(0)] * (E - 1)]
    for n in range(1, N):
        rhs = emul(ep(P1, n - 1), A[n - 1])
        if n >= 2:
            rhs = [p - q for p, q in zip(rhs, emul(ep(P2, n - 2), A[n - 2]))]
        A.append(emul(rhs, einv(ep(P0, n))))
    f = [[A[n][k] for n in range(N)] for k in range(E)]

    def mul(p, q):
        return [sum(p[i] * q[k - i] for i in range(k + 1)) for k in range(N)]

    def inv(p):
        r_ = [Fr(1) / p[0]] + [Fr(0)] * (N - 1)
        for k in range(1, N):
            r_[k] = -sum(p[i] * r_[k - i] for i in range(1, k + 1)) / p[0]
        return r_

    def theta(p):
        return [k * p[k] for k in range(N)]

    def exp0(p):
        r_ = [Fr(1)] + [Fr(0)] * (N - 1)
        for k in range(1, N):
            r_[k] = sum(j * p[j] * r_[k - j] for j in range(1, k + 1)) / k
        return r_

    def compose(p, q):
        r_ = [Fr(0)] * N
        pw = [Fr(1)] + [Fr(0)] * (N - 1)
        for k in range(N):
            if p[k]:
                r_ = [a + p[k] * b for a, b in zip(r_, pw)]
            pw = mul(pw, q)
        return r_

    i0 = inv(f[0])
    g1 = mul(f[1], i0)
    g2 = mul(f[2], i0)
    H = [a - b / 2 for a, b in zip(g2, mul(g1, g1))]
    dtd = inv([Fr(1) + c if k == 0 else c for k, c in enumerate(theta(g1))])
    K = mul(dtd, theta(mul(dtd, theta(H))))
    K[0] += 1
    qz = [Fr(0)] + exp0(g1)[:N - 1]
    zq = [Fr(0), Fr(1)] + [Fr(0)] * (N - 2)
    for _ in range(N):
        c = compose(qz, zq)
        zq = [a - (b - (1 if k == 1 else 0)) for k, (a, b) in enumerate(zip(zq, c))]
    Kq = compose(K, zq)
    assert Kq[:4] == [1, 16, 1152, 99520]



# ---------------------------------------------------------------- irreducibility without a degree bound
# The F_n operators as stored with the Mellin-route record (sum_i t^(2i) Q_i(theta)); the guard
# also re-derives them with _t_ode, so a drift in either is caught.
F_OPS = {2: [-u**2, (u + 1)**2],
         3: [4 * u**3, -(u + 1)**3],
         4: [9 * u**2 * (u - 2)**2, -2 * u**2 * (5 * u**2 + 1), (u + 1)**4],
         5: [-64 * u**2 * (u - 2)**2 * (u - 1), 4 * u**3 * (5 * u**2 + 3), -(u + 1)**5]}
# the ring's triangle operator L_G (eq:ring_op), reducible by construction: the positive control
LG_OP = [-4 * u * (u - 3)**2, 5 * u**3 - 7 * u**2 - 5 * u + 3, -u * (u + 2) * (u + 3)]
lam = sp.Symbol("lam")
R_ = sp.Rational


def _local_exponents(Q):
    """Local exponents of sum_i t^(2i) Q_i(theta) at every singular point of P^1.
    At t = 0 they are the roots of Q_0, at oo the roots of Q_K(-lam) (t = 1/s, theta -> -theta);
    at a finite p != 0 they come from the d/dt form.  Returns ({point: sorted exponents},
    fuchsian_everywhere)."""
    K = len(Q) - 1
    n = max(sp.Poly(q, u).degree() for q in Q)
    ind = {sp.Integer(0): sp.Poly(Q[0].subs(u, lam), lam), sp.oo: sp.Poly(Q[K].subs(u, -lam), lam)}
    a = [sp.Integer(0)] * (n + 1)                              # theta^k = sum_j S(k,j) t^j D^j
    for i, q in enumerate(Q):
        for (k,), c in sp.Poly(sp.expand(q), u).terms():
            for j in range(k + 1):
                a[j] += c * t**(2 * i) * stirling(k, j, kind=2) * t**j
    a = [sp.expand(c) for c in a]
    e = sp.Symbol("e")
    for p in sp.roots(sp.Poly(sp.cancel(a[n] / t**n), t)):
        expr = sp.expand(sum(c.subs(t, p + e) * sp.ff(lam, j) * e**(n - j) for j, c in enumerate(a)))
        pe = sp.Poly(expr, e)
        low = min(m[0] for m in pe.monoms())
        ind[p] = sp.Poly(pe.as_expr().coeff(e, low), lam)
    fuchsian = all(P.degree() == n for P in ind.values())
    exps = {}
    for p, P in ind.items():
        rts = sp.roots(P)
        assert sum(rts.values()) == P.degree()
        exps[p] = sorted((r for r, m in rts.items() for _ in range(m)), key=lambda z: sp.re(z))
    return exps, fuchsian


def _fuchs_factor_orders(exps, n):
    """Orders k (1 <= k < n) for which a right factor is NOT excluded.  A right factor of order k
    takes, at each singular point p, k of the operator's exponents (a sub-multiset), and each
    apparent singularity it adds raises its Fuchs sum by at least 1, so it needs
        sum_p (smallest sum of k exponents at p) <= (|S| - 2) k (k - 1) / 2.
    No bound on degrees or on the number of apparent singularities enters."""
    S = len(exps)
    possible = []
    for k in range(1, n):
        lhs = sum(min(sum(sp.re(x) for x in c) for c in itertools.combinations(ex, k))
                  for ex in exps.values())
        if lhs <= R_((S - 2) * k * (k - 1), 2):
            possible.append(k)
    return possible


# local exponents in t at 0, the finite singular points and oo: the record's Fuchs table
F_EXPONENTS = {
    2: {0: [0, 0], 1: [0, 0], -1: [0, 0], sp.oo: [1, 1]},
    3: {0: [0, 0, 0], 2: [0, R_(1, 2), 1], -2: [0, R_(1, 2), 1], sp.oo: [1, 1, 1]},
    4: {0: [0, 0, 2, 2], 1: [0, 1, 1, 2], -1: [0, 1, 1, 2], 3: [0, 1, 1, 2], -3: [0, 1, 1, 2],
        sp.oo: [1, 1, 1, 1]},
    5: {0: [0, 0, 1, 2, 2], 2: [0, 1, R_(3, 2), 2, 3], -2: [0, 1, R_(3, 2), 2, 3],
        4: [0, 1, R_(3, 2), 2, 3], -4: [0, 1, R_(3, 2), 2, 3], sp.oo: [1, 1, 1, 1, 1]},
}
LG_EXPONENTS = {0: [0, 3, 3], 1: [0, 0, 1], -1: [0, 0, 1], 2: [R_(-3, 2), 0, 1], -2: [R_(-3, 2), 0, 1],
                sp.oo: [0, 2, 3]}


def _check_exponent_table(Q, table):
    n = max(sp.Poly(q, u).degree() for q in Q)
    exps, fuchsian = _local_exponents(Q)
    assert fuchsian
    got = {sp.sympify(p): list(ex) for p, ex in exps.items()}
    want = {sp.sympify(p): [sp.sympify(x) for x in ex] for p, ex in table.items()}
    assert got == want
    total = sum(sum(ex) for ex in exps.values())
    assert total == R_((len(exps) - 2) * n * (n - 1), 2)        # Fuchs' relation for L itself
    return exps, n


def test_fn_operators_irreducible_by_fuchs():
    """Rejects a reducible F_n operator (n = 2..5) with no bound on the degree of a factor (the
    degree-12 rank test above cannot exclude a factor of higher degree; this argument can), and
    rejects a wrong operator or a mis-read exponent: the Mellin-route operator must equal the
    stored one, its local exponents must be the recorded table, Fuchs' relation must close, and
    the exponent bound must exclude every order k < n.  Paired with the positive control below."""
    for n in (2, 3, 4, 5):
        Q = _t_ode(n)
        assert len(Q) == len(F_OPS[n]) and all(sp.expand(a - b) == 0 for a, b in zip(Q, F_OPS[n]))
        exps, order = _check_exponent_table(F_OPS[n], F_EXPONENTS[n])
        assert order == n
        assert _fuchs_factor_orders(exps, n) == []


def test_fuchs_argument_sees_the_reducible_triangle():
    """Positive control: rejects an exponent argument that calls everything irreducible (a max
    in place of the min, a bound with the wrong sign or offset).  On the triangle operator L_G,
    which factors exactly as M o (theta - r), the same function must leave a right factor of
    order 1 possible."""
    exps, n = _check_exponent_table(LG_OP, LG_EXPONENTS)
    assert n == 3
    assert 1 in _fuchs_factor_orders(exps, n)


# ---------------------------------------------------------------- the banana relation
# Printed operators of g(y), F_n = t^-1 g(y), y = 1/(16 t^2), as sum_j y^j P_j(theta_y).
# n = 4 is AESZ #58; n = 5 is the fifth-order operator whose pullback is AESZ #39 (Almkvist-Zudilin).
BANANA_PRINTED = {
    2: {0: th**2, 1: -4 * (2 * th + 1)**2},
    3: {0: th**3, 1: -8 * (2 * th + 1)**3},
    4: {0: th**4, 1: -4 * (2 * th + 1)**2 * (10 * th**2 + 10 * th + 3),
        2: 144 * (2 * th + 1)**2 * (2 * th + 3)**2},
    5: {0: th**5, 1: -8 * (2 * th + 1)**3 * (5 * th**2 + 5 * th + 2),
        2: 2**10 * (th + 1) * (2 * th + 1)**2 * (2 * th + 3)**2},
}


def _g_operator(Q):
    """sum_i t^(2i) Q_i(theta_t) acting on F = t^-1 g(y), y = 1/(16 t^2), as an operator on g:
    theta_t = -2 theta_y, theta_t t^-1 = t^-1 (theta_t - 1), t^(2i) = (16 y)^-i; multiplied on the
    left by (16 y)^K and normalised so the y^0 part is monic.  Returns {j: P_j(theta)}."""
    K = len(Q) - 1
    Lg = {K - i: sp.expand(16**(K - i) * q.subs(u, -2 * th - 1)) for i, q in enumerate(Q)}
    lead = sp.Poly(Lg[0], th).LC()
    return {j: sp.expand(c / lead) for j, c in Lg.items()}


def _walk_moment(N, m):
    """W_N(2m): the sum over a_1 + ... + a_N = m of the squared multinomial (m; a_1, ..., a_N)."""
    row = [Fr(1, int(sp.factorial(k))**2) for k in range(m + 1)]
    acc = [Fr(1)] + [Fr(0)] * m
    for _ in range(N):
        acc = [sum(acc[i] * row[k - i] for i in range(k + 1)) for k in range(m + 1)]
    return acc[m] * int(sp.factorial(m))**2


def test_tree_ladder_is_legendre_hadamard_banana():
    """Rejects a wrong dictionary entry: a wrong variable (y = 1/(4 t^2), or the t^-1 dropped),
    a wrong Hadamard factor (binom(2m,m) to the first power, or none), a banana with the wrong
    number of propagators (W_n in place of W_{n-1}), or a mis-copied AESZ operator.  In exact
    integers: the operator of g, transformed from the Mellin-route operator of F_n, must equal
    the printed operator and must annihilate sum binom(2m,m)^2 W_{n-1}(2m) y^m through m = 30,
    its unique holomorphic solution at the MUM point y = 0.  The operator is irreducible (the Fuchs
    guard above, in t; a factor in y would pull back to one in t), so it is that series' minimal
    operator, i.e. the Hadamard product.  Not guarded here: that sum_k W_N(2k) s^-k is the banana
    period (cited, Vanhove 2014)."""
    Mx = 30
    legendre = BANANA_PRINTED[2]          # theta^2 - 4y(2 theta+1)^2 is the operator of binom(2m,m)^2
    for m in range(1, Mx):
        assert (legendre[0].subs(th, m) * sp.binomial(2 * m, m)**2
                + legendre[1].subs(th, m - 1) * sp.binomial(2 * m - 2, m - 1)**2) == 0
    for n in (2, 3, 4, 5):
        Lg = _g_operator(_t_ode(n))
        assert set(Lg) == set(BANANA_PRINTED[n])
        assert all(sp.expand(Lg[j] - BANANA_PRINTED[n][j]) == 0 for j in Lg)
        c = [sp.binomial(2 * m, m)**2 * _walk_moment(n - 1, m) for m in range(Mx)]
        for m in range(Mx):
            assert sum(Lg[j].subs(th, m - j) * c[m - j] for j in Lg if m - j >= 0) == 0


# ---------------------------------------------------------------- the Ising class, n <= 8
# At t = 1 the ladder is the Ising class of Bessel moments, F_n(1) = c_{n,0} = int K0(x)^n dx.
# Record: docs/predictions/2026-10-17_ising_class.md (committed before the computation).
def test_ising_thresholds_and_their_exponents_n_le_8():
    """Rejects a wrong threshold set for n = 2..8 (a shifted family such as (n - 2j)^2, a missing
    or an extra finite singular point), a flipped parity (t = 1 singular for odd n), and a
    mis-read local type.  At every threshold t = n - 1 - 2j > 0 the exponents must be
    0, 1, ..., n - 2 together with (n - 2)/2, so one integer exponent is doubled (a logarithm)
    for even n and one is a half-integer (a square-root branch) for odd n.  The operators for
    n = 6, 7, 8, which no stored table pins, are checked on F_n by the next test."""
    for n in range(2, 9):
        exps, fuchsian = _local_exponents(_t_ode(n))
        assert fuchsian
        pos = [n - 1 - 2 * j for j in range(n) if n - 1 - 2 * j > 0]
        finite = sorted((p for p in exps if p not in (0, sp.oo)), key=float)
        assert finite == sorted(pos + [-p for p in pos])
        assert (1 in exps) == (n % 2 == 0)
        for p in pos:
            ex = exps[p]
            assert ex == sorted(list(range(n - 1)) + [R_(n - 2, 2)])
            if n % 2 == 0:
                assert all(x.is_integer for x in ex) and len(set(ex)) == n - 1
            else:
                assert len(set(ex)) == n and [x for x in ex if not x.is_integer] == [R_(n - 2, 2)]


@pytest.mark.slow
def test_ising_operators_annihilate_F6_F7_F8_by_quadrature():
    """Rejects a wrong Mellin-route operator at n = 6, 7, 8: applied to F_n by quadrature at
    t = 0.7 it must cancel to the working precision.  Control: the same operator applied to the
    neighbouring integral (one more K0 factor) must not, so the cancellation is not generic."""
    for n in (6, 7, 8):
        assert _annihilation_residual(n, "0.7") < 1e-15
        assert _annihilation_residual(n, "0.7", power=n) > 1e-6


# LMFDB 24.6.a.a, the weight-6 level-24 newform attached to the sixth moment (Lim-Tu-Yu):
# a_p for every prime p <= 97 as printed there; 2 and 3 are the bad primes.
A24_PRIMES = {2: 0, 3: -9, 5: -34, 7: -240, 11: -124, 13: 46, 17: 1954, 19: -1924, 23: 2840,
              29: -8922, 31: -4648, 37: -4362, 41: -2886, 43: 11332, 47: 7008, 53: -22594,
              59: -28, 61: -6386, 67: -39076, 71: -54872, 73: 21034, 79: 26632, 83: 56188,
              89: 64410, 97: -116158}


def _hecke_coeffs(ap, k, M, bad):
    a = [0] * (M + 1)
    a[1] = 1
    for m in range(2, M + 1):
        val = 1
        for p, r in sp.factorint(m).items():
            if p in bad:
                val *= ap[p]**r
            else:
                x0, x1 = 1, ap[p]
                for _ in range(r - 1):
                    x0, x1 = x1, ap[p] * x1 - p**(k - 1) * x0
                val *= x1
        a[m] = val
    return a


def _eta2_4_eta4_4(M):
    """a_1..a_M of eta(2 tau)^4 eta(4 tau)^4 = q prod (1 - q^2m)^4 (1 - q^4m)^4."""
    c = [0] * (M + 1)
    c[0] = 1
    for m in range(1, M + 1):
        for step in (2 * m, 4 * m):
            if step > M:
                continue
            for _ in range(4):
                for i in range(M, step - 1, -1):
                    c[i] -= c[i - step]
    return [0] + c[:M]


def _completed_L(a, N, k, s, eps, cut=1):
    """Lambda(s) = (sqrt N/2 pi)^s Gamma(s) L(s) from the two-cut form, independent of the
    cut > 0 exactly when N, eps and the coefficients are those of a newform."""
    b = 2 * mp.pi / mp.sqrt(N)
    return mp.fsum(a[m] * (mp.gammainc(s, b * m * cut) / (b * m)**s
                           + eps * mp.gammainc(k - s, b * m / cut) / (b * m)**(k - s))
                   for m in range(1, len(a)) if a[m])


def _L_value(a, N, k, s, eps):
    return _completed_L(a, N, k, s, eps) / ((mp.sqrt(N) / (2 * mp.pi))**s * mp.gamma(s))


def _cut_dependence(a, N, k, eps):
    s_ = mp.mpf(k) / 2 + mp.mpf("0.37")
    l1 = _completed_L(a, N, k, s_, eps)
    return abs(_completed_L(a, N, k, s_, eps, cut=mp.mpf("1.25")) / l1 - 1)


def test_ising_moments_are_the_registered_values():
    """Rejects a wrong closed form or L-value for c_{n,0} = int K0^n:
    c_{3,0} = 3 Gamma(1/3)^6 / (32 pi 2^(2/3));  c_{4,0} = (pi^4/4) 4F3(1/2,1/2,1/2,1/2; 1,1,1; 1)
    = 4 pi^2 L(f, 2), f = eta(2t)^4 eta(4t)^4 of level 8;  c_{6,0} = 108 pi^2 L(g, 4), g = 24.6.a.a.
    The L-value instrument is tested before it is used: the two-cut form must not depend on the
    cut for the right level and sign (eps = +1 at level 8, -1 at level 24) and must for a wrong
    sign, a wrong level, or the a_3 = +9 twin of g.  Control: level 16 in place of 8 misses c_{4,0}."""
    with mp.workdps(32):
        moment = {n: mp.quad(lambda x, n=n: mp.besselk(0, x)**n, [0, 0.5, 1, 2, 4, 8, 16, 40])
                  for n in (3, 4, 6)}
        f = _eta2_4_eta4_4(40)
        g = _hecke_coeffs(A24_PRIMES, 6, 97, bad={2, 3})
        g_twin = _hecke_coeffs({**A24_PRIMES, 3: 9}, 6, 97, bad={2, 3})
        assert _cut_dependence(f, 8, 4, +1) < 1e-25
        assert _cut_dependence(g, 24, 6, -1) < 1e-25
        assert _cut_dependence(f, 8, 4, -1) > 1e-3
        assert _cut_dependence(g, 24, 6, +1) > 1e-3
        assert min(_cut_dependence(f, 16, 4, e) for e in (1, -1)) > 1e-3
        assert min(_cut_dependence(g, 12, 6, e) for e in (1, -1)) > 1e-3
        assert min(_cut_dependence(g_twin, 24, 6, e) for e in (1, -1)) > 1e-3

        c30 = 3 * mp.gamma(mp.mpf(1) / 3)**6 / (32 * mp.pi * mp.cbrt(4))
        assert abs(moment[3] / c30 - 1) < 1e-25
        assert abs(moment[4] / (mp.pi**4 / 4 * mp.hyper([0.5] * 4, [1, 1, 1], 1)) - 1) < 1e-25
        assert abs(moment[4] / (4 * mp.pi**2 * _L_value(f, 8, 4, 2, +1)) - 1) < 1e-25
        assert abs(moment[4] / (4 * mp.pi**2 * _L_value(f, 16, 4, 2, +1)) - 1) > 0.1
        assert abs(moment[6] / (108 * mp.pi**2 * _L_value(g, 24, 6, 4, -1)) - 1) < 1e-25
