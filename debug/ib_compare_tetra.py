"""IB3: compare the tetrahedron operator (d/dx form on Y/t^3, mod P) with a stored Y operator
sum_i x^i Q_i(theta_t), theta_t = t d/dt = 2 x d/dx, acting on Y.  On G = Y/t^3 the stored operator is
sum_i x^i Q_i(2 theta_x + 3).  usage: python ib_compare_tetra.py <ib op json> <stored Q json>"""
import json, sys
import sympy as sp
from flint import nmod_poly
from ib_reduce import P
from ib_kite_de import RF
from ib_compare import monic, right_divide

u = sp.Symbol("u")


def stored_Y_to_D(path):
    d = json.load(open(path))
    Q = [sp.Poly(sp.expand(sp.sympify(s.replace("^", "**"), locals={"u": u}).subs(u, 2 * u + 3)), u) for s in d["Q"]]
    coef = {}
    for i, q in enumerate(Q):
        for (k,), c in q.terms():
            for j in range(k + 1):
                s2 = sp.functions.combinatorial.numbers.stirling(k, j, kind=2)
                if s2:
                    coef[(j, i + j)] = coef.get((j, i + j), 0) + int(c) * int(s2)
    r = max(j for j, _ in coef)
    return [RF(nmod_poly([coef.get((j, e), 0) % P for e in range(max(e for jj, e in coef if jj == j) + 1)], P))
            for j in range(r + 1)]


ib = json.load(open(sys.argv[1]))
L = [RF(nmod_poly(p, P)) for p in ib["polys"]]
S = stored_Y_to_D(sys.argv[2])
print("IBP order", len(L) - 1, " stored order", len(S) - 1)
if len(L) == len(S):
    print("equal after making monic:", all((x - y).is_zero() for x, y in zip(monic(L), monic(S))))
rem = right_divide(L, S)
print("stored right-divides IBP operator:", len(rem) == 0)
