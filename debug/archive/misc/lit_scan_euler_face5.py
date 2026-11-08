"""Lit-scan companion: Euler count of the x5 -> infinity face of the kite's Lee-Pomeransky polynomial
(massless Coulomb line 5: G = x5*h + (terms without x5); h = dG/dx5 is the polynomial of the graph with line 5
deleted, i.e. the soft region q -> 0).  Does h's count drop at z0?
usage: python lit_scan_euler_face5.py zgen ztest"""
import sys, itertools
import sympy as sp
a, b, g1, g2 = sp.Integer(1), sp.Rational(3, 2), sp.Rational(7, 3), sp.Rational(8, 3)
x = sp.symbols("x1:6")
U = sp.expand((x[0] + x[1]) * (x[2] + x[3]) + x[4] * (x[0] + x[1] + x[2] + x[3]))
k1, k2, q1, q2, w, z = sp.symbols("k1 k2 q1 q2 w z")
mom = {1: (k1, k2), 2: (w - k1, -k2), 3: (k1 - q1, k2 - q2), 4: (w - k1 + q1, -k2 + q2), 5: (q1, q2)}
m2 = {1: a**2, 2: b**2, 3: g1**2, 4: g2**2, 5: 0}
Q = sp.expand(sum(x[i - 1] * (mom[i][0]**2 + mom[i][1]**2 + m2[i]) for i in mom))
X = [k1, q1]
H = sp.Matrix(2, 2, lambda i, j: sp.diff(Q, X[i], X[j]) / 2)
gv = sp.Matrix([sp.diff(Q, xx) / 2 for xx in X]).subs({k1: 0, q1: 0})
Q0 = Q.subs({k1: 0, q1: 0, k2: 0, q2: 0})
Fpol = sp.expand(sp.expand(U * Q0 - (gv.T * H.adjugate() * gv)[0]).subs(w**2, z))
G = sp.expand(U + Fpol)
assert sp.degree(G, x[4]) == 1
h = sp.expand(sp.diff(G, x[4]))
print("h =", sp.factor(h))
nu = [sp.Rational(3, 7), sp.Rational(5, 11), sp.Rational(2, 13), sp.Rational(7, 17)]
s = sp.Rational(9, 23)
t = sp.Symbol("t")
vs = list(x[:4])
for zs in sys.argv[1:]:
    g = sp.expand(h.subs(z, sp.Rational(zs)))
    eqs = [sp.expand(nu[i] * g - s * vs[i] * sp.diff(g, vs[i])) for i in range(4)] + [t * sp.prod(vs) * g - 1]
    eqs = [sp.together(e).as_numer_denom()[0] for e in eqs]
    gens = [t] + vs
    gb = sp.groebner(eqs, *gens, order="grevlex", modulus=32003)
    if gb.exprs == [1]:
        print(zs, 0); continue
    lead = [sp.Poly(e, *gens, modulus=32003).monoms(order="grevlex")[0] for e in gb.exprs]
    nv = len(gens); bound = []
    for i in range(nv):
        c = [m[i] for m in lead if all(m[j] == 0 for j in range(nv) if j != i)]
        bound.append(min(c) if c else None)
    if None in bound:
        print(zs, "not zero-dimensional"); continue
    print(zs, sum(1 for mono in itertools.product(*[range(bb) for bb in bound])
                  if not any(all(mono[j] >= m[j] for j in range(nv)) for m in lead)), flush=True)
