"""Preview (mod one prime): count formal Laurent-series solutions x^n0 C[[s]] of the off-slice tetrahedron operator
at a root rho of a factor f of the leading coefficient, working in the field F_P[y]/(f) (rho = y).  At an irreducible
factor of degree > 1 this is the apparent test without needing a prime where f has a root."""
import json, sys
from fractions import Fraction as F
from pathlib import Path
from flint import nmod_poly, fq_default_ctx, fq_default_poly_ctx, fmpz_mod_poly_ctx
import os
d = json.loads((Path(__file__).parent / os.environ.get("IB_OPFILE", "data/ib_tetra_pf_op.json")).read_text())
P = d["P"]; ops = [nmod_poly(p, P) for p in d["polys"]]; order = len(ops) - 1
fr = lambda q: q.numerator * pow(q.denominator, -1, P) % P


def count(f, n0=0, N=60):
    K = fq_default_ctx(modulus=fmpz_mod_poly_ctx(P)([int(c) for c in f.coeffs()]))
    R = fq_default_poly_ctx(K)
    y = K.gen()
    shift = R([y, K(1)])
    Ps = []
    for q in ops:
        res = R([K(0)])
        for c in reversed([int(v) for v in q.coeffs()]):
            res = res * shift + R([K(c)])
        Ps.append([res[i] for i in range(res.degree() + 1)])
    v = [min((i for i, c in enumerate(cs) if not c.is_zero()), default=10**9) for cs in Ps]
    delta = min(v[k] - k for k in range(order + 1))
    cols = list(range(n0, N + 1))
    rows = []
    for m in range(n0 + delta, N + delta + 1):
        row = [K(0)] * len(cols)
        for k, cs in enumerate(Ps):
            for i, c in enumerate(cs):
                n = m - i + k
                if c.is_zero() or not (n0 <= n <= N):
                    continue
                ff = 1
                for j in range(k):
                    ff = ff * (n - j) % P
                row[n - n0] = row[n - n0] + c * K(ff)
        rows.append(row)
    # rank by Gaussian elimination
    rank, ncol = 0, len(cols)
    for col in range(ncol):
        piv = next((r for r in range(rank, len(rows)) if not rows[r][col].is_zero()), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = rows[rank][col] ** -1
        rows[rank] = [e * inv for e in rows[rank]]
        for r in range(len(rows)):
            if r != rank and not rows[r][col].is_zero():
                t = rows[r][col]
                rows[r] = [a - t * b for a, b in zip(rows[r], rows[rank])]
        rank += 1
    return ncol - rank, delta


lead = ops[-1]
rats = [F(0), F(1, 4), F(1), F(9, 4), F(4), F(9), F(49, 4), F(16), F(81, 4)]
rem = lead
for x in rats:
    lin = nmod_poly([(-fr(x)) % P, 1], P)
    while rem(fr(x)) == 0:
        rem = rem // lin
if "--rational" in sys.argv:
    for x in rats:
        lin = nmod_poly([(-fr(x)) % P, 1], P)
        for n0 in (0, -1, -2, -3, -4):
            dim, _ = count(lin, n0)
            print(f"x = {x}: solutions in x^({n0}) C[[s]]: {dim} of {order}", flush=True)
else:
    for f, e in rem.factor()[1]:
        fm = f * pow(int(f.coeffs()[-1]), -1, P)
        dim, delta = count(fm)
        print(f"rest factor of degree {f.degree()} (mult {e}) over F_P: analytic solutions {dim} of {order} "
              f"-> {'apparent' if dim == order else 'GENUINE'}", flush=True)
