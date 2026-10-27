"""PK record: Taylor coefficients of an operator's analytic solution at z = 0 from its recurrence.
Given phi_0..phi_{f-1} (the free indices, the nonnegative integer exponents at 0), generate phi_f..phi_N and
compare with a direct series.  Optionally perturb one coefficient of P_k by 201/200 (negative control).
usage: python pko_recurrence.py <op json> <series json> <N> [perturb k]
Environment: PKO_DPS (working digits, default 120); PKO_INPUT (series json supplying the free coefficients,
default the comparison series); PKO_SHOW=1 prints per-coefficient errors.
Operator form: sum_k p_k(z) (d/dz)^k phi = 0, p_k given by integer coefficient lists ("polys")."""
import json
import os
import sys
from fractions import Fraction as F

import mpmath as mp

sys.set_int_max_str_digits(10**6)
mp.mp.dps = int(os.environ.get("PKO_DPS", "120"))
op = json.load(open(sys.argv[1]))
ser = json.load(open(sys.argv[2]))
N = int(sys.argv[3])
polys = [[int(v) for v in p] for p in op["polys"]]
if len(sys.argv) > 4 and sys.argv[4] != "none":
    k = int(sys.argv[4])
    i = min(j for j, c in enumerate(polys[k]) if c != 0)   # lowest nonzero coefficient: enters the recurrence from n = 0
    polys[k] = [F(v) for v in polys[k]]
    polys[k][i] = polys[k][i] * F(201, 200)
    print("perturbed P_%d coefficient of z^%d by 201/200" % (k, i))
r = len(polys) - 1
v = [min(j for j, c in enumerate(p) if c != 0) for p in polys]
delta = min(v[k] - k for k in range(r + 1))


def ff(n, k):
    out = 1
    for j in range(k):
        out *= (n - j)
    return out


def ind(n):
    return sum(F(polys[k][v[k]]) * ff(n, k) for k in range(r + 1) if v[k] - k == delta)


free = [n for n in range(0, N + 1) if ind(n) == 0]
direct = [mp.mpf(s) for s in ser["phi_n"]]
inp = [mp.mpf(s) for s in json.load(open(os.environ["PKO_INPUT"]))["phi_n"]] if os.environ.get("PKO_INPUT") else direct
phi = {}
for n in range(N + 1):
    if n in free:
        phi[n] = inp[n]
        continue
    m = n + delta
    s = mp.mpf(0)
    for k in range(r + 1):
        for i, c in enumerate(polys[k]):
            if c == 0:
                continue
            nn = m - i + k
            if nn < 0 or nn >= n:
                continue
            s += (mp.mpf(F(c).numerator) / F(c).denominator) * ff(nn, k) * phi[nn]
    q = ind(n)
    phi[n] = -s * q.denominator / q.numerator
gen = [n for n in range(N + 1) if n not in free and n < len(direct)]
worst = max(abs(phi[n] / direct[n] - 1) for n in gen)
if os.environ.get("PKO_SHOW"):
    print([(n, mp.nstr(abs(phi[n] / direct[n] - 1), 2)) for n in gen[:12]])
print("free indices", free, "; generated", gen[0], "to", gen[-1], "; worst relative error vs direct", mp.nstr(worst, 3))
