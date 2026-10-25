"""Round trip: stored L9 (Q-form on Y) -> integer d/dx form on Y/t^3 -> d_to_Q -> compare with L9 up to scale."""
import json
import sympy as sp
from ib_tetra_lift import d_to_Q
u = sp.Symbol("u")
d = json.load(open("data/ring_ext_QG_5c_ylog_g1.5_r9_d25.json"))
Q0 = [sp.expand(sp.sympify(s.replace("^", "**"), locals={"u": u})) for s in d["Q"]]
coef = {}
for i, q in enumerate(Q0):
    for (k,), c in sp.Poly(sp.expand(q.subs(u, 2 * u + 3)), u).terms():
        for j in range(k + 1):
            s2 = sp.functions.combinatorial.numbers.stirling(k, j, kind=2)
            if s2:
                coef[(j, i + j)] = coef.get((j, i + j), 0) + int(c) * int(s2)
r = max(j for j, _ in coef)
ip = [[coef.get((j, e), 0) for e in range(max(e for jj, e in coef if jj == j) + 1)] for j in range(r + 1)]
Q1, lo = d_to_Q(ip)
print("orders", r, "x-degrees", len(Q0) - 1, len(Q1) - 1, "shift", lo)
# compare up to a common rational scale
ratio = None; ok = len(Q0) == len(Q1)
for a, b in zip(Q0, Q1):
    a, b = sp.Poly(a, u), sp.Poly(b, u)
    if a.is_zero != b.is_zero:
        ok = False; break
    if a.is_zero:
        continue
    for (k,), c in a.terms():
        rr = sp.Rational(b.coeff_monomial(u ** k), c)
        ratio = rr if ratio is None else ratio
        ok &= rr == ratio
    ok &= len(a.terms()) == len(b.terms())
print("round trip equal up to scale", ratio, ":", ok)
# negative control: perturb one coefficient
ip[2][3] += 1
Q2, _ = d_to_Q(ip)
print("perturbed differs:", any(sp.expand(a * ratio - b) != 0 for a, b in zip(Q0, Q2)))
