"""PK record, stage-1 diagnosis (after scoring PK1): is the extra genuine point of the P2 operator the kite's
leading Landau singularity (all five lines on shell)?
Symanzik form for the kite, lines 1 (k, a), 2 (w-k, b), 3 (k-q, g1), 4 (w-k+q, g2), 5 (q, 0), Euclidean, z = w^2:
F(alpha) = U * C - B^T adj(M) B, homogeneous of degree 3.  Landau: dF/dalpha_i = 0 for all i (then F = 0 by Euler).
Solved by elimination with sympy (Groebner basis, alpha_5 = 1) for the z-values.
usage: python pko_landau_leading.py a b g1 g2"""
import sys
from fractions import Fraction as F

import sympy as sp

a, b, g1, g2 = [sp.Rational(str(F(x))) for x in sys.argv[1:5]]
al = sp.symbols("a1:6")
z = sp.Symbol("z")
k1, k2, q1, q2, w = sp.symbols("k1 k2 q1 q2 w")
# momenta in a plane, w along the first axis; scalar products only through z = w^2
mom = {1: (k1, k2), 2: (w - k1, -k2), 3: (k1 - q1, k2 - q2), 4: (w - k1 + q1, -k2 + q2), 5: (q1, q2)}
m5sq = sp.Rational(sys.argv[5]) if len(sys.argv) > 5 else 0   # optional regulator mass^2 for the Coulomb line
m2 = {1: a**2, 2: b**2, 3: g1**2, 4: g2**2, 5: m5sq}
Q = sum(al[i - 1] * (mom[i][0]**2 + mom[i][1]**2 + m2[i]) for i in mom)
Q = sp.expand(Q)
# minimise over loop momenta: F = U * Q(min); U = det of the Hessian / 4 per component
X = [k1, q1]
H = sp.Matrix(2, 2, lambda i, j: sp.diff(Q, X[i], X[j]) / 2)
U = sp.expand(H.det())
g = sp.Matrix([sp.diff(Q, x) / 2 for x in X]).subs({k1: 0, q1: 0})
Q0 = Q.subs({k1: 0, q1: 0, k2: 0, q2: 0})
Fpol = sp.expand(U * Q0 - (g.T * H.adjugate() * g)[0])
Fpol = sp.expand(Fpol.subs(w**2, z))
assert w not in Fpol.free_symbols
from itertools import combinations
tt = sp.Symbol("tt")
target = sp.Rational(-111, 1120)
print("U =", U)
for ncon in ((0,) if m5sq else (0, 1, 2)):
    for con in combinations(range(5), ncon):
        keep = [j for j in range(5) if j not in con]
        Fr = Fpol.subs({al[j]: 0 for j in con})
        Ur = U.subs({al[j]: 0 for j in con})
        if Ur == 0 or z not in Fr.free_symbols:
            continue
        norm = keep[-1]
        eqs = [sp.diff(Fr, al[j]).subs(al[norm], 1) for j in keep]
        prod = Ur.subs(al[norm], 1)
        for j in keep[:-1]:
            prod = prod * al[j]
        free = [al[j] for j in keep[:-1]]
        G = sp.groebner(eqs + [tt * prod - 1], tt, *free, z, order="lex")
        zp = [q for q in G.exprs if q.free_symbols <= {z}]
        txt = "; ".join(str(sp.factor(q)) for q in zp) if zp else ("unit ideal" if G.exprs == [1] else "no eliminant")
        hit = any(sp.Poly(q, z).eval(target) == 0 for q in zp if q.free_symbols)
        if m5sq:
            for q in zp:
                if q.free_symbols:
                    print("  real roots:", [sp.N(r, 12) for r in sp.Poly(q, z).real_roots()])
        print("contracted", [c + 1 for c in con], ":", txt, "| -111/1120 a root:", hit, flush=True)
