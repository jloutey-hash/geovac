"""Lit-scan companion: the singular locus C_z of {G=0} in the torus (a curve, for every z: the degenerate
full-polytope face) with z kept symbolic.  Where does the closure of C_z meet a coordinate hyperplane x_i = 0, and
where is C_z itself singular?  Each gives a polynomial in z; test whether -111/1120 is a root.
(Whitney-stratification-style probe, after Helmer-Papathanasiou-Tellander arXiv:2402.14787.)"""
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
t = sp.Symbol("t")
z0 = sp.Rational(-111, 1120)
eqs = [G] + [sp.expand(v * sp.diff(G, v)) for v in x]
gb = sp.groebner(eqs + [t * sp.prod(x) - 1], t, *x, z, order="lex")
closure = [e for e in gb.exprs if not e.has(t)]
print("closure generators:", len(closure), flush=True)
for e in closure:
    print("   ", sp.factor(e))


def zpolys(extra, label):
    g2b = sp.groebner(closure + extra, *x, z, order="lex")
    zp = [e for e in g2b.exprs if e.free_symbols <= {z}]
    out = "; ".join(str(sp.factor(e)) for e in zp) if zp else ("unit" if g2b.exprs == [1] else "none (z free)")
    hit = any(e.free_symbols and sp.Poly(e, z).eval(z0) == 0 for e in zp)
    print(label, ":", out, "| z0 root:", hit, flush=True)


for i in range(5):
    zpolys([x[i]], f"closure meets x{i+1}=0")
J = sp.Matrix([[sp.diff(e, v) for v in x] for e in closure])
# singular points of the curve: Jacobian rank < 4
minors = []
for rows in itertools.combinations(range(J.rows), 4):
    for cols in itertools.combinations(range(5), 4):
        minors.append(sp.expand(J.extract(list(rows), list(cols)).det()))
        if len(minors) > 60:
            break
    if len(minors) > 60:
        break
gbs = sp.groebner(closure + minors + [t * sp.prod(x) - 1], t, *x, z, order="lex")
zp = [e for e in gbs.exprs if e.free_symbols <= {z}]
print("singular points of C_z in the torus:", "; ".join(str(sp.factor(e)) for e in zp) if zp else gbs.exprs[:2],
      "| z0 root:", any(e.free_symbols and sp.Poly(e, z).eval(z0) == 0 for e in zp), flush=True)
