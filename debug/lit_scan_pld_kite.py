"""Lit-scan companion (topology_rung_ninth_point_memo): principal-Landau-determinant style face scan of the kite.

Lee-Pomeransky polynomial G = U + F of the kite (lines 1 (k,a), 2 (w-k,b), 3 (k-q,g1), 4 (w-k+q,g2), 5 (q,m5)),
Euclidean, z = w^2.  Following Fevola-Mizera-Telen (arXiv:2311.16219), the principal A-determinant of G is the
product over faces Gamma of Newt(G) of the discriminants Delta_Gamma (G_Gamma and all x_i dG_Gamma/dx_i vanish at a
point of the torus).  For each face we test, at fixed z, whether that system has a solution in the torus
(Groebner basis with a saturation variable; unit ideal = no solution).
Run at z0 = -111/1120 and at controls (known genuine points, and a random z).
usage: python lit_scan_pld_kite.py a b g1 g2 [m5sq]
"""
import sys
import itertools
from fractions import Fraction as Fr

import numpy as np
import sympy as sp
from scipy.spatial import ConvexHull

a, b, g1, g2 = [sp.Rational(str(Fr(x))) for x in sys.argv[1:5]]
m5sq = sp.Rational(sys.argv[5]) if len(sys.argv) > 5 else sp.Integer(0)
x = sp.symbols("x1:6")
z = sp.Symbol("z")
U = (x[0] + x[1]) * (x[2] + x[3]) + x[4] * (x[0] + x[1] + x[2] + x[3])
U = sp.expand(U)
# F from the momentum-space quadratic form (same construction as debug/pko_landau_leading.py)
k1, k2, q1, q2, w = sp.symbols("k1 k2 q1 q2 w")
mom = {1: (k1, k2), 2: (w - k1, -k2), 3: (k1 - q1, k2 - q2), 4: (w - k1 + q1, -k2 + q2), 5: (q1, q2)}
m2 = {1: a**2, 2: b**2, 3: g1**2, 4: g2**2, 5: m5sq}
Q = sp.expand(sum(x[i - 1] * (mom[i][0]**2 + mom[i][1]**2 + m2[i]) for i in mom))
X = [k1, q1]
H = sp.Matrix(2, 2, lambda i, j: sp.diff(Q, X[i], X[j]) / 2)
assert sp.expand(H.det() - U) == 0
gv = sp.Matrix([sp.diff(Q, xx) / 2 for xx in X]).subs({k1: 0, q1: 0})
Q0 = Q.subs({k1: 0, q1: 0, k2: 0, q2: 0})
Fpol = sp.expand(sp.expand(U * Q0 - (gv.T * H.adjugate() * gv)[0]).subs(w**2, z))
assert w not in Fpol.free_symbols
G = sp.expand(U + Fpol)
P = sp.Poly(G, *x)
terms = P.terms()
expo = np.array([t[0] for t in terms], dtype=float)
print("U =", U)
print("F =", sp.factor_terms(Fpol))
print("monomials of G:", len(terms))

# face lattice of Newt(G): facets from qhull, then all intersections of facet vertex sets
hull = ConvexHull(expo)
eqs = np.unique(np.round(hull.equations, 9), axis=0)
facets = []
for e in eqs:
    vals = expo @ e[:-1] + e[-1]
    facets.append(frozenset(np.where(np.abs(vals) < 1e-7)[0]))
facets = list(set(facets))
faces = set(facets)
frontier = set(facets)
while frontier:
    new = set()
    for f1 in frontier:
        for f2 in facets:
            f = f1 & f2
            if len(f) >= 1 and f not in faces:
                new.add(f)
    faces |= new
    frontier = new
faces.add(frozenset(range(len(terms))))
print("faces (incl. whole polytope):", len(faces), " facets:", len(facets))

tt = sp.Symbol("tt")


def torus_solution(face, zval, timeout_vars=None):
    gpoly = sum(terms[i][1] * sp.prod([x[j]**terms[i][0][j] for j in range(5)]) for i in face)
    gpoly = sp.expand(gpoly.subs(z, zval))
    if gpoly == 0:
        return None
    vars_in = [v for v in x if gpoly.has(v)]
    eq = [gpoly] + [sp.expand(v * sp.diff(gpoly, v)) for v in vars_in]
    eq = [e for e in eq if e != 0]
    # torus: x_i != 0; dehomogenize nothing (G is not homogeneous); saturate by product
    sat = tt * sp.prod(vars_in) - 1
    gb = sp.groebner(eq + [sat], tt, *vars_in, order="grevlex")
    return not (len(gb.exprs) == 1 and gb.exprs[0] == 1)


zvals = {"z0=-111/1120": sp.Rational(-111, 1120)}
for extra in sys.argv[6:]:
    zvals[extra] = sp.Rational(extra)
zvals["random -37/101"] = sp.Rational(-37, 101)
zface = [f for f in faces if any(sp.Poly(terms[i][1], z).degree() > 0 for i in f)]
print("faces whose restriction depends on z:", len(zface))
results = {}
for name, zv in zvals.items():
    results[name] = {f: torus_solution(f, zv) for f in faces}
rnd = results["random -37/101"]
print("faces with a torus solution at the random z (degenerate: discriminant identically zero):",
      sum(1 for f in faces if rnd[f]))
for name in zvals:
    if name.startswith("random"):
        continue
    hits = [f for f in sorted(faces, key=len) if results[name][f] and not rnd[f]]
    print()
    print(f"{name}: faces with a torus solution here but not at the random z: {len(hits)}")
    for f in hits:
        mons = [sp.prod([x[j]**terms[i][0][j] for j in range(5)]) for i in sorted(f)]
        print("  ", len(f), mons, flush=True)
