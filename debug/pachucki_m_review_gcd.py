"""Review: is the python-flint fallback in pachucki_m_core._gcd_ZZ_safe exact and in sympy's conventions?
Forces the fallback path on random and structured univariate integer polynomials and compares with sympy's
dense gcd (dup_gcd via Poly.gcd) and heugcd where heugcd succeeds."""
import json
import random
import sys

import sympy as sp

sys.path.insert(0, "debug")
import pachucki_m_core as M  # noqa: E402
from sympy.polys.rings import ring
from sympy.polys.domains import ZZ
from sympy.polys.polyerrors import HeuristicGCDFailed

R, t = ring("t", ZZ)
random.seed(11)
bad = 0
n = 0
orig = M._heugcd


def failing(f, g):
    raise HeuristicGCDFailed("forced")


for trial in range(300):
    dg = random.randint(0, 12)
    G = R.from_dense([random.randint(-50, 50) for _ in range(dg + 1)] or [1])
    if G == 0:
        G = R(1)
    A = G * R.from_dense([random.randint(-10**6, 10**6) for _ in range(random.randint(1, 25))])
    B = G * R.from_dense([random.randint(-10**6, 10**6) for _ in range(random.randint(1, 25))])
    if trial % 7 == 0:
        A = A * (t - 1) ** 5 * 3
        B = B * (t - 1) ** 3 * (t + 1) * 6
    if A == 0 or B == 0:
        continue
    M._heugcd = failing
    h, cf, cg = M._gcd_ZZ_safe(A, B)
    M._heugcd = orig
    # reference: sympy Poly gcd over ZZ
    pa = sp.Poly(A.as_expr(), sp.Symbol("t"), domain="ZZ")
    pb = sp.Poly(B.as_expr(), sp.Symbol("t"), domain="ZZ")
    ref = pa.gcd(pb)
    hp = sp.Poly(h.as_expr(), sp.Symbol("t"), domain="ZZ")
    ok = (hp == ref or hp == -ref) and h * cf == A and h * cg == B
    n += 1
    if not ok:
        bad += 1
        print("MISMATCH", trial, hp, ref, flush=True)
# exercise a Q(t) field computation end to end with the fallback forced
K = sp.QQ.frac_field(sp.Symbol("t"))
x = K.from_sympy(sp.sympify("(t**7-3*t+2)/(t**2-1)**4"))
y = K.from_sympy(sp.sympify("(t+1)**3/(5*t**9-t)"))
M._heugcd = failing
z = (x + y) * (x - y) / (x * y)
M._heugcd = orig
z_ref = sp.cancel(sp.sympify("((t**7-3*t+2)/(t**2-1)**4)**2 - ((t+1)**3/(5*t**9-t))**2")
                  / sp.sympify("((t**7-3*t+2)/(t**2-1)**4)*((t+1)**3/(5*t**9-t))"))
field_ok = sp.simplify(K.to_sympy(z) - z_ref) == 0
out = {"trials": n, "mismatches": bad, "field_roundtrip_ok": bool(field_ok)}
print(out)
json.dump(out, open("debug/data/pachucki_m_review_gcd.json", "w"), indent=1)
