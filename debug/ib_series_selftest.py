"""Self-test of ib_crt.series_residual: the stored gamma = 3/2 kite operator on the stored series, and the
same operator with one coefficient changed by 1/200 (must fail)."""
import json, sys
from fractions import Fraction
from ib_compare import stored_to_D
import sympy as sp
from ib_crt import series_residual
u = sp.Symbol("u")
d = json.load(open("data/mg_5_op_kite_g1.5_stored.json"))
Q = [sp.Poly(sp.sympify(s.replace("^", "**"), locals={"u": u}), u) for s in d["Q"]]
coef = {}
for i, q in enumerate(Q):
    for (k,), c in q.terms():
        for j in range(k + 1):
            s2 = sp.functions.combinatorial.numbers.stirling(k, j, kind=2)
            if s2:
                coef[(j, i + j)] = coef.get((j, i + j), 0) + int(c) * int(s2)
r = max(j for j, _ in coef)
op = [[Fraction(coef.get((j, e), 0)) for e in range(max(e for jj, e in coef if jj == j) + 1)] for j in range(r + 1)]
ser = json.load(open(sys.argv[1]))["phi_n"][:300]
print("stored op on stored series:", series_residual(op, ser, dps=200))
op[3][5] = op[3][5] * (1 + Fraction(1, 200)) if op[3][5] else Fraction(1)
print("perturbed op:", series_residual(op, ser, dps=200))
