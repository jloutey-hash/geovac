import json
from pathlib import Path
from flint import nmod_mat, nmod_poly, nmod
import ib_kite_de as K
from ib_reduce import P
import ib_tetra_pf_degrees as G
m = G.ms[0]
ys = [r[2][m] for r in G.res]
f = G.interp(ys); N, D = G.ratrec(f)
print("ratrec degs", N.degree(), D.degree(), "gcd deg", N.gcd(D).degree())
print("check all pts", all(N(x) == D(x) * nmod(y, P) for x, y in zip(G.xs, ys)))
print("sample x types", type(G.res[0][0]), G.res[0][0] < P, min(r[0] for r in G.res))
dn, dd = N.degree(), D.degree(); nunk = dn + dd + 2
xs = [r[0] for r in G.res][:-10]; ys2 = ys[:-10]
rows = [[pow(x, k, P) for k in range(dn + 1)] + [(-y * pow(x, k, P)) % P for k in range(dd + 1)] for x, y in zip(xs[:nunk-1], ys2[:nunk-1])]
rows.append([0] * (dn + 1) + [K._NORM[k] for k in range(dd + 1)])
Mx = nmod_mat(rows, P)
print("rank", Mx.rank(), "of", nunk)
try:
    sol = Mx.solve(nmod_mat([[0]] * (nunk - 1) + [[1]], P)); print("solved")
except Exception as e:
    print("solve error", type(e).__name__, e)
Mp = nmod_mat(rows[:-1], P)
print("point-rows rank", Mp.rank())
print("norm.D", sum(K._NORM[k] * int(D.coeffs()[k]) for k in range(dd + 1)) % P)
print("distinct xs", len(set(xs)), len(xs))
print("P value", P, "flint modulus ok", nmod(P - 1, P) + 1 == 0)
print("D coeffs", [int(v) for v in D.coeffs()])
print("NORM[:18]", K._NORM[:18])
for m2 in G.ms[:6]:
    y2 = [r[2][m2] for r in G.res]; N2, D2 = G.ratrec(G.interp(y2))
    print(m2, "norm.D", sum(K._NORM[k] * int(c) for k, c in enumerate(D2.coeffs())) % P, "D lead", int(D2.coeffs()[-1]), "D0", int(D2.coeffs()[0]))
