"""Apparent-singularity test for an integer d/dz-form operator: factor its leading polynomial over Z, and at a
root of each factor (found modulo a prime where the factor has one) count the formal power-series solutions:
r of them means every solution is analytic there (apparent point).  usage: python ib_apparent.py <op json> [skip roots]"""
import json, sys
from fractions import Fraction as F
from flint import fmpz_poly, nmod_poly, nmod_mat
from sympy import prevprime

d = json.load(open(sys.argv[1]))
ops_z = [fmpz_poly([int(v) for v in p]) for p in d["polys"]]
order = len(ops_z) - 1
skip = [F(s) for s in sys.argv[2:]]
fac = ops_z[-1].factor()
primes = []
p = (1 << 61) - 1
for _ in range(60):
    p = prevprime(p); primes.append(p)


def series_dim(rho, p, N=40):
    ops = [nmod_poly([int(c) % p for c in q.coeffs()], p) for q in ops_z]
    sh = nmod_poly([rho, 1], p)
    Ps = [q(sh) for q in ops]
    v = [min((i for i, c in enumerate(q.coeffs()) if int(c)), default=10**9) for q in Ps]
    delta = min(v[k] - k for k in range(order + 1))
    M = N + delta
    rows = []
    for m in range(M + 1):
        row = [0] * (N + 1)
        for k, q in enumerate(Ps):
            co = [int(c) for c in q.coeffs()]
            for i, c in enumerate(co):
                n = m - i + k
                if c and 0 <= n <= N and n - k >= 0:
                    ff = 1
                    for j in range(k):
                        ff = ff * (n - j) % p
                    row[n] = (row[n] + c * ff) % p
        rows.append(row)
    A = nmod_mat(rows, p)
    return (N + 1) - A.rank()


for f, e in fac[1]:
    deg = f.degree()
    if deg == 1:
        root = F(-int(f.coeffs()[0]), int(f.coeffs()[1]))
        if root in skip:
            continue
    for p in primes:
        fp = nmod_poly([int(c) % p for c in f.coeffs()], p)
        roots = [r for r, m in fp.factor()[1] if r.degree() == 1]
        if roots:
            r1 = roots[0]
            rho = (-int(r1.coeffs()[0]) * pow(int(r1.coeffs()[1]), -1, p)) % p
            dim = series_dim(rho, p)
            print(f"factor degree {deg} (multiplicity {e}): analytic solutions {dim} of {order}"
                  f" -> {'apparent' if dim == order else 'GENUINE'} (checked mod a prime with a root)")
            break
    else:
        print(f"factor degree {deg}: no root modulo the primes tried")
