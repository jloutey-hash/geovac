"""Loop test PT2, step D7: exact peeling of the chain control's operator,
L_C = M5 . R1 . (theta - 2),  R1 = (theta - 3) - t^2 (theta + 1)  (solution t^3/(1-t^2)^2),
then M5's structure: f3 / f2 as sub or quotient (rational gauge, as D5), log-free solutions of
M5 and M5* and their minimal orders, quadric count.  Operators are sum_i t^{2i} Q_i(theta)."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C
import density_loop_D_gauge as D
import density_loop_D4_factors as F4
import density_loop_D5_rational as F5
import density_loop_D6_peel as F6
u = sp.Symbol("u")
OUT = Path(__file__).parent / "data" / "density_loop_D7.json"

def rdiv(M, R):
    """M = N . R exactly; returns N (list of theta-polys) or raises."""
    k = len(R) - 1
    n = len(M) - 1 - k
    N = []
    for i in range(n + 1):
        acc = M[i]
        for j in range(1, min(i, k) + 1):
            acc -= N[i - j].subs(u, u + 2 * j) * R[j]
        q, r = sp.div(sp.expand(acc), sp.expand(R[0]), u)
        assert r == 0, ("not divisible at", i, sp.factor(r))
        N.append(sp.expand(q))
    for i in range(n + 1, len(M)):
        acc = M[i] - sum(N[i - j].subs(u, u + 2 * j) * R[j] for j in range(0, k + 1) if 0 <= i - j <= n)
        assert sp.expand(acc) == 0, ("remainder at", i)
    return N

B = json.loads((Path(__file__).parent / "data" / "density_loop_B.json").read_text())
QC = [sp.sympify(q) for q in B["QC"]]
M6 = rdiv(QC, [u - 2])
M5 = rdiv(M6, [u - 3, -(u + 1)])
print("M5:", flush=True)
for i, q in enumerate(M5):
    print(f"  t^{2*i}:", sp.factor(q), flush=True)
lead = sum(sp.Poly(q, u).coeff_monomial(u**5) * sp.Symbol("t") ** (2 * i) for i, q in enumerate(M5))
print("  leading theta^5 coefficient:", sp.factor(lead), flush=True)
f2 = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=F4.NT)
f3 = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=F4.NT, scale=Fraction(1, 4))
rec = {"M5": [str(sp.factor(q)) for q in M5]}
for ref, f, r in [("f3", f3, 3), ("f2", f2, 2)]:
    for side, QQ in [("R", M5), ("L", F4.adjoint(M5))]:
        hit = 0
        for b in range(4):
            for c in range(4):
                hit = F5.test(QQ, f, r, 10, b, c)
                if hit: break
            if hit: break
        rec[f"M5|{ref}|{side}"] = hit
        print(f"M5 vs {ref} ({side}): {'FOUND' if hit else 'none'}", flush=True)
F6.run("M5", M5)
F6.run("M5* (adjoint)", F4.adjoint(M5))
# quadric count at t0 = 1/2
sols, rr = C.taylor_basis(M5, sp.Rational(1, 2), 80)
qc = C.quadric_count(sols, 80)
print("M5 order", rr, " quadric count at t0 = 1/2:", qc, flush=True)
rec["quadric_count"] = qc
OUT.write_text(json.dumps(rec, indent=1))
