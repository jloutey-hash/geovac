"""Lit-scan companion: the singular locus of {G = U + F = 0} in the torus for the kite (= the 'full-polytope face'
of the principal-A-determinant scan, equivalently U = F = 0 with grad U parallel to grad F: the classical
second-type / mixed condition).  Print its Groebner basis (dehomogenized) at a random z and at z0, to see what the
positive-dimensional locus is and whether it changes at z0.
usage: python lit_scan_sing_kite.py a b g1 g2 m5sq z..."""
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
G = sp.expand(U + Fpol)
t = sp.Symbol("t")
for zs in sys.argv[6:]:
    g = sp.expand(G.subs(z, sp.Rational(zs)))
    eqs = [g] + [sp.expand(v * sp.diff(g, v)) for v in x]
    gb = sp.groebner(eqs + [t * sp.prod(x) - 1], t, *x, order="lex")
    print("z =", zs)
    for e in gb.exprs:
        if not e.has(t):
            print("   ", sp.factor(e))
    print(flush=True)
