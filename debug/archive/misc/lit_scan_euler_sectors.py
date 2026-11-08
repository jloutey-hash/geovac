"""Lit-scan companion: Euler-characteristic count (critical points of the Lee-Pomeransky log-likelihood) for
every sector of the kite (lines contracted: x_j = 0), at z0 = -111/1120 and at a random z.  A sector whose count
drops at z0 carries the singularity (Helmer-Papathanasiou-Tellander arXiv:2402.14787 criterion, per sector).
usage: python lit_scan_euler_sectors.py a b g1 g2 m5sq zgen ztest"""
import sys, itertools
from fractions import Fraction as Fr
import sympy as sp

a, b, g1, g2, m5sq = [sp.Rational(str(Fr(t))) for t in sys.argv[1:6]]
zgen, ztest = sp.Rational(sys.argv[6]), sp.Rational(sys.argv[7])
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
p = int(sys.argv[8]) if len(sys.argv) > 8 else 32003


def count(g, vs):
    eqs = [sp.expand(nu[x.index(v)] * g - s * v * sp.diff(g, v)) for v in vs]
    eqs.append(t * sp.prod(vs) * g - 1)
    eqs = [sp.together(e).as_numer_denom()[0] for e in eqs]
    gens = [t] + list(vs)
    gb = sp.groebner(eqs, *gens, order="grevlex", modulus=p)
    if gb.exprs == [1]:
        return 0
    lead = [sp.Poly(e, *gens, modulus=p).monoms(order="grevlex")[0] for e in gb.exprs]
    nv = len(gens)
    bound = []
    for i in range(nv):
        c = [m[i] for m in lead if all(m[j] == 0 for j in range(nv) if j != i)]
        if not c:
            return "inf"
        bound.append(min(c))
    return sum(1 for mono in itertools.product(*[range(bb) for bb in bound])
               if not any(all(mono[j] >= m[j] for j in range(nv)) for m in lead))


for ncon in range(0, 4):
    for con in itertools.combinations(range(5), ncon):
        vs = [x[j] for j in range(5) if j not in con]
        g = sp.expand(G.subs({x[j]: 0 for j in con}))
        if sp.expand(U.subs({x[j]: 0 for j in con})) == 0:
            continue
        c0 = count(sp.expand(g.subs(z, zgen)), vs)
        c1 = count(sp.expand(g.subs(z, ztest)), vs)
        flag = "  <-- DROP" if c0 != c1 else ""
        print("contracted", [j + 1 for j in con], "generic", c0, "test", c1, flag, flush=True)
