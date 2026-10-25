"""Preview (mod one prime, not a score): singular points of the off-slice tetrahedron operator for Y/t^3 in x = t^2.
Rational roots of the leading coefficient by search, F_P factorisation of the rest, local exponents and the count
of formal power-series solutions at each root (r of r = apparent)."""
import json
from fractions import Fraction as F
from pathlib import Path
from flint import nmod_poly, nmod_mat
import os
d = json.loads((Path(__file__).parent / os.environ.get("IB_OPFILE", "data/ib_tetra_pf_op.json")).read_text())
P = d["P"]; ops = [nmod_poly(p, P) for p in d["polys"]]; order = len(ops) - 1
fr = lambda q: q.numerator * pow(q.denominator, -1, P) % P


def local(rho, N=70):
    sh = nmod_poly([rho, 1], P)
    Ps = [q(sh) for q in ops]
    v = [min((i for i, c in enumerate(q.coeffs()) if int(c)), default=10**9) for q in Ps]
    delta = min(v[k] - k for k in range(order + 1))
    ind = nmod_poly([0], P)
    for k in range(order + 1):
        if v[k] - k == delta:
            ff = nmod_poly([1], P)
            for j in range(k):
                ff = ff * nmod_poly([(-j) % P, 1], P)
            ind = ind + ff * int(Ps[k].coeffs()[v[k]])
    roots, rem = [], ind
    for den in (1, 2, 4, 8):
        for j in range(-80 * den, 80 * den + 1):
            r = F(j, den)
            if r.denominator != den:
                continue
            lin = nmod_poly([(-fr(r)) % P, 1], P)
            while rem.degree() > 0 and rem(fr(r)) == 0:
                rem = rem // lin; roots.append(r)
    fuchs = delta == v[-1] - order
    M = N + delta
    rows = []
    for m in range(M + 1):
        row = [0] * (N + 1)
        for k, q in enumerate(Ps):
            for i, c in enumerate(q.coeffs()):
                c = int(c); n = m - i + k
                if c and 0 <= n <= N and n >= k:
                    ff = 1
                    for j in range(k):
                        ff = ff * (n - j) % P
                    row[n] = (row[n] + c * ff) % P
        rows.append(row)
    dim = (N + 1) - nmod_mat(rows, P).rank()
    return fuchs, sorted(roots), rem.degree(), dim


lead = ops[-1]
print("order", order, "leading degree", lead.degree())
found, rem = [], lead
for den in (1, 2, 3, 4, 8, 9, 16, 25, 36, 64):
    for j in range(-120 * den, 120 * den + 1):
        x = F(j, den)
        if x.denominator != den:
            continue
        lin = nmod_poly([(-fr(x)) % P, 1], P); mult = 0
        while rem.degree() > 0 and rem(fr(x)) == 0:
            rem = rem // lin; mult += 1
        if mult:
            found.append((x, mult))
print("rational roots (x, multiplicity):", [(str(x), m) for x, m in sorted(found)])
fac = rem.factor()[1]
print("rest: degree", rem.degree(), "F_P factor degrees", sorted((f.degree(), e) for f, e in fac))
for x, m in sorted(found) + [(F(25, 4), 0), (F(25), 0), (F(7, 3), 0)]:
    fu, ex, unf, dim = local(fr(x))
    print(f"x = {x} (mult {m}): {'Fuchsian' if fu else 'IRREGULAR'}, exponents {[str(e) for e in ex]}"
          + (f" +{unf} unfound" if unf else "") + f"; analytic solutions {dim} of {order}")
for f, e in fac:
    if f.degree() == 1:
        rho = (-int(f.coeffs()[0]) * pow(int(f.coeffs()[1]), -1, P)) % P
        fu, ex, unf, dim = local(rho)
        print(f"F_P root of the rest: analytic solutions {dim} of {order} -> {'apparent' if dim == order else 'GENUINE'}")
