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
Record: docs/predictions/2026-10-16_linking_ceiling.md."""
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


def test_four_density_operator_annihilates_F4_by_quadrature():
    Q = _t_ode(4)
    mp.mp.dps = 20
    J = max(sp.Poly(q, u).degree() for q in Q)
    al, be = [sp.Integer(1)], [sp.Integer(0)]                # theta^j K0 = al K0 + be (u K1)
    for _ in range(J):
        a_, b_ = al[-1], be[-1]
        al.append(sp.expand(u * sp.diff(a_, u) - u**2 * b_))
        be.append(sp.expand(u * sp.diff(b_, u) - a_))
    alf = [sp.lambdify(u, a_, "mpmath") for a_ in al]
    bef = [sp.lambdify(u, b_, "mpmath") for b_ in be]
    T = mp.mpf("0.7")
    thF = [mp.quad(lambda x, j=j: mp.besselk(0, x)**3 * (alf[j](T * x) * mp.besselk(0, T * x)
                   + bef[j](T * x) * T * x * mp.besselk(1, T * x)), [0, 1, 4, mp.inf])
           for j in range(J + 1)]
    tot = mp.mpf(0)
    scale = mp.mpf(0)
    for i, q in enumerate(Q):
        for (deg,), cf in sp.Poly(q, u).terms():
            term = T**(2 * i) * mp.mpf(int(cf)) * thF[deg]
            tot += term
            scale += abs(term)
    assert abs(tot) / scale < 1e-15


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
