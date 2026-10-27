"""Phase-1 cost probe for docs/predictions/2026-10-27_kite_ninth_point_origin.md.
Times, at RANDOM generic (z, eps) only, (1) the Groebner basis of the Lee-Pomeransky critical-point ideal modulo p,
(2) the standard-monomial count, and (3) the traces of multiplication by t, x1..x5 (normal forms of v*b).
Nothing here evaluates anything at z9, at a threshold, or on any predicted curve.
usage: python nq_probe_cost.py p z_residue eps_residue
"""
import itertools
import sys
import time
from fractions import Fraction as Fr

import sympy as sp

p = int(sys.argv[1])
zv = int(sys.argv[2]) % p
ev = int(sys.argv[3]) % p
a, b, g1, g2 = Fr(1), Fr(3, 2), Fr(7, 3), Fr(8, 3)


def res(q):
    q = Fr(q)
    return q.numerator * pow(q.denominator, -1, p) % p


x = sp.symbols("x1:6")
t = sp.Symbol("t")
z, eps = sp.symbols("z eps")
U = sp.expand((x[0] + x[1]) * (x[2] + x[3]) + x[4] * (x[0] + x[1] + x[2] + x[3]))
k1, k2, q1, q2, w = sp.symbols("k1 k2 q1 q2 w")
mom = {1: (k1, k2), 2: (w - k1, -k2), 3: (k1 - q1, k2 - q2), 4: (w - k1 + q1, -k2 + q2), 5: (q1, q2)}
m2 = {1: sp.Rational(a.numerator, a.denominator)**2, 2: sp.Rational(b.numerator, b.denominator)**2,
      3: sp.Rational(g1.numerator, g1.denominator)**2, 4: sp.Rational(g2.numerator, g2.denominator)**2, 5: eps}
Q = sp.expand(sum(x[i - 1] * (mom[i][0]**2 + mom[i][1]**2 + m2[i]) for i in mom))
X = [k1, q1]
H = sp.Matrix(2, 2, lambda i, j: sp.diff(Q, X[i], X[j]) / 2)
gv = sp.Matrix([sp.diff(Q, xx) / 2 for xx in X]).subs({k1: 0, q1: 0})
Q0 = Q.subs({k1: 0, q1: 0, k2: 0, q2: 0})
Fpol = sp.expand(sp.expand(U * Q0 - (gv.T * H.adjugate() * gv)[0]).subs(w**2, z))
G = sp.expand(U + Fpol)
nu = [Fr(3, 7), Fr(5, 11), Fr(2, 13), Fr(7, 17), Fr(4, 19)]
s = Fr(9, 23)

# reduce G's coefficients mod p with z, eps substituted as residues
Gp = sp.Poly(G, *x, z, eps)
terms = {}
for mono, c in Gp.terms():
    c = Fr(int(sp.numer(c)), int(sp.denom(c)))
    key = mono[:5]
    terms[key] = (terms.get(key, 0) + res(c) * pow(zv, mono[5], p) * pow(ev, mono[6], p)) % p
g = sum(sp.Integer(c) * sp.prod([x[j]**m[j] for j in range(5)]) for m, c in terms.items() if c)
t0 = time.process_time()
eqs = [sp.expand(res(nu[i]) * g - res(s) * x[i] * sp.diff(g, x[i])) for i in range(5)]
eqs.append(t * x[0] * x[1] * x[2] * x[3] * x[4] * g - 1)
gb = sp.groebner(eqs, t, *x, order="grevlex", modulus=p)
t1 = time.process_time()
gens = (t,) + tuple(x)
lead = [sp.Poly(e, *gens, modulus=p).monoms(order="grevlex")[0] for e in gb.exprs]
nv = 6
bound = [min(m[i] for m in lead if all(m[j] == 0 for j in range(nv) if j != i)) for i in range(nv)]
std = [mono for mono in itertools.product(*[range(bb) for bb in bound])
       if not any(all(mono[j] >= m[j] for j in range(nv)) for m in lead)]
t2 = time.process_time()
std_set = {m: i for i, m in enumerate(std)}
traces = []
for v in gens:
    tr = 0
    for mono in std:
        bexpr = sp.prod([gens[j]**mono[j] for j in range(nv)])
        _, r = gb.reduce(v * bexpr)
        rp = sp.Poly(r, *gens, modulus=p)
        tr = (tr + int(rp.coeff_monomial(bexpr))) % p
    traces.append(tr)
t3 = time.process_time()
print(f"p={p} z={zv} eps={ev} count={len(std)} gb_cpu={t1 - t0:.1f}s count_cpu={t2 - t1:.1f}s "
      f"traces_cpu={t3 - t2:.1f}s traces={traces}", flush=True)
