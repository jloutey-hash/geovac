"""Lit-scan companion: Euler-characteristic drop test for the kite (Helmer-Papathanasiou-Tellander criterion,
arXiv:2402.14787: on a Landau component the signed Euler characteristic of the very affine variety
{x in (C*)^5 : G(x) != 0} drops; |chi| = number of master integrals, Bitoun-Bogner-Klausen-Panzer arXiv:1712.09215).
|chi| = number of critical points of log-likelihood  sum_i nu_i log x_i - s log G  for generic nu, s
(Huh; Lee-Pomeransky / Mizera critical-point count).  Count them by a Groebner basis at fixed z, in exact arithmetic
modulo a prime (sympy, modulus=p) to keep coefficient growth down.
usage: python lit_scan_euler_kite.py a b g1 g2 m5sq z1 [z2 ...]
"""
import sys
from fractions import Fraction as Fr

import sympy as sp

a, b, g1, g2, m5sq = [sp.Rational(str(Fr(t))) for t in sys.argv[1:6]]
x = sp.symbols("x1:6")
U = sp.expand((x[0] + x[1]) * (x[2] + x[3]) + x[4] * (x[0] + x[1] + x[2] + x[3]))
k1, k2, q1, q2, w, z = sp.symbols("k1 k2 q1 q2 w z")
mom = {1: (k1, k2), 2: (w - k1, -k2), 3: (k1 - q1, k2 - q2), 4: (w - k1 + q1, -k2 + q2), 5: (q1, q2)}
m2 = {1: a**2, 2: b**2, 3: g1**2, 4: g2**2, 5: m5sq}
Q = sp.expand(sum(x[i - 1] * (mom[i][0]**2 + mom[i][1]**2 + m2[i]) for i in mom))
X = [k1, q1]
H = sp.Matrix(2, 2, lambda i, j: sp.diff(Q, X[i], X[j]) / 2)
gv = sp.Matrix([sp.diff(Q, xx) / 2 for xx in X]).subs({k1: 0, q1: 0})
Q0 = Q.subs({k1: 0, q1: 0, k2: 0, q2: 0})
Fpol = sp.expand(sp.expand(U * Q0 - (gv.T * H.adjugate() * gv)[0]).subs(w**2, z))
assert w not in Fpol.free_symbols
G = sp.expand(U + Fpol)
nu = [sp.Rational(3, 7), sp.Rational(5, 11), sp.Rational(2, 13), sp.Rational(7, 17), sp.Rational(4, 19)]
s = sp.Rational(9, 23)
t = sp.Symbol("t")
p = 32003


def count(zv):
    g = sp.expand(G.subs(z, zv))
    eqs = [sp.expand(nu[i] * g - s * x[i] * sp.diff(g, x[i])) for i in range(5)]
    eqs.append(t * x[0] * x[1] * x[2] * x[3] * x[4] * g - 1)
    # clear denominators, reduce mod p
    eqs = [sp.Poly(sp.together(e).as_numer_denom()[0], t, *x) for e in eqs]
    gb = sp.groebner([e.as_expr() for e in eqs], t, *x, order="grevlex", modulus=p)
    # count standard monomials (zero-dimensional)
    lead = [sp.Poly(e, t, *x, modulus=p).monoms(order="grevlex")[0] for e in gb.exprs]
    nv = 6
    for i in range(nv):
        if not any(all(m[j] == 0 for j in range(nv) if j != i) for m in lead):
            return "not zero-dimensional"
    bound = [min(m[i] for m in lead if all(m[j] == 0 for j in range(nv) if j != i)) for i in range(nv)]
    import itertools
    cnt = 0
    for mono in itertools.product(*[range(bb) for bb in bound]):
        if not any(all(mono[j] >= m[j] for j in range(nv)) for m in lead):
            cnt += 1
    return cnt


for zs in sys.argv[6:]:
    print(zs, "critical points:", count(sp.Rational(zs)), flush=True)
