"""Atlas gap memo, door A: exact first-order Stark matrices in shells n = 2, 3, 4 (z in the (l, m) basis,
exact radial integrals with sympy) and their eigenvalues vs (3/2) n (n1 - n2); also the
cotangent-sum identity sum_{k=1}^{p-1} cot^2(pi k/p) = (p-1)(p-2)/3 behind the lens-space signature eta.
Output: debug/data/atlas_gap_linear_stark.json
"""
import json
import sympy as sp
from sympy.physics.hydrogen import R_nl
from sympy.physics.wigner import gaunt

r = sp.symbols('r', positive=True)

def z_elem(n, l1, m, l2):
    # <n l1 m| r cos(theta) |n l2 m>; cos(theta) = sqrt(4 pi/3) Y_10
    rad = sp.integrate(R_nl(n, l1, r, 1) * R_nl(n, l2, r, 1) * r ** 3, (r, 0, sp.oo))
    # <l1 m|Y_10|l2 m> = (-1)^m * gaunt(l1, 1, l2, -m, 0, m)
    ang = (-1) ** m * gaunt(l1, 1, l2, -m, 0, m)
    return sp.simplify(rad * sp.sqrt(4 * sp.pi / 3) * ang)

out = {}
for n in (2, 3, 4):
    res = {}
    for m in range(0, n):
        ls = list(range(abs(m), n))
        M = sp.Matrix(len(ls), len(ls), lambda i, j: z_elem(n, ls[i], m, ls[j]))
        ev = sorted([sp.nsimplify(sp.simplify(e)) for e in M.eigenvals(multiple=True)])
        # prediction: (3/2) n k, k = n1 - n2, n1 + n2 = n - |m| - 1
        pred = sorted([sp.Rational(3, 2) * n * (n1 - (n - abs(m) - 1 - n1)) for n1 in range(0, n - abs(m))])
        res[str(m)] = {"eigs": [str(e) for e in ev], "pred": [str(x) for x in pred], "match": bool(ev == pred)}
        print(n, m, ev, ev == pred)
    out[str(n)] = res
cot = {}
for p in range(2, 12):
    s = sum(sp.cot(sp.pi * k / p) ** 2 for k in range(1, p))
    val = sp.nsimplify(sp.N(s, 50), rational=True)
    cot[str(p)] = [str(val), str(sp.Rational((p - 1) * (p - 2), 3)), bool(abs(sp.N(s, 50) - sp.Rational((p - 1) * (p - 2), 3)) < 1e-40)]
print(cot)
out["cot2_sum"] = cot
json.dump(out, open("debug/data/atlas_gap_linear_stark.json", "w"), indent=1)
