"""exch_reverse_review_det.py -- REVIEW of debug/exch_reverse_memo.md Q4 at R values the memo did not use.
det of the 9 x 9 coordinate matrix of (K, ..., K^(8)) on the nine functions, as an integer polynomial in t at
rational R (helper calculus, tests/_exact_order).  Residual = det with t, t-1, t+1 stripped.  Report real roots in
(0, 1) of each residual and of their gcd, plus an explicit value of the residual at t = 1/2, 1/3, 9/10.
Run: python debug/exch_reverse_review_det.py CASE R1 R2 ... -> debug/data/exch_reverse_review_det_<CASE>.json"""
import json
import sys
import time
from fractions import Fraction

import sympy as sp

sys.path.insert(0, "tests")
import _pachucki_extension as X  # noqa: E402
import _exact_order as EO  # noqa: E402
import _exchange_m_nonzero as M  # noqa: E402

M.install_gcd_fallback()
case = sys.argv[1]
Rs = [Fraction(r) for r in sys.argv[2:]]
t = sp.Symbol("t")
T0 = time.time()
if case in ("2psig", "3dsig"):
    cft = json.load(open("tests/data/paper58_pachucki_ext_2026_11_18.json"))["sym"][case]
else:
    cft = json.load(open("tests/data/paper58_exchange_m_2026_11_18.json"))["sym"][case]
ctx = X.Ctx("Qt")
K = EO.stored_K(ctx, cft)
rows = EO.derivative_rows(ctx, K, 8)
per_row, other = EO.denominator_audit(rows)
assert not other, other
cleared = EO.clear_rows(ctx, rows, per_row)
out = {"case": case, "rows": []}
res_all = []
for R0 in Rs:
    N, scale = EO.determinant_polynomial(cleared, R0)
    mult, res = EO.strip_linear(N)
    P = sp.Poly(list(reversed(res)), t, domain="ZZ")
    roots01 = [str(r) for r in sp.Poly(P, t).intervals(eps=sp.Rational(1, 10**6)) if r[0][0] >= 0 and r[0][1] <= 1]
    vals = {s: str(sp.N(P.eval(sp.Rational(s)), 6)) for s in ("1/2", "1/3", "9/10", "99/100")}
    row = {"R": str(R0), "degree": len(N) - 1, "mult": mult, "residual_degree": P.degree(),
           "n_roots_(0,1)": int(P.count_roots(0, 1)), "root_intervals_(0,1)": roots01, "residual_at": vals,
           "time_s": round(time.time() - T0, 1)}
    out["rows"].append(row)
    res_all.append(res)
    print(row, flush=True)
pt = {}
for tv in (Fraction(1, 2), Fraction(1, 3), Fraction(9, 10), Fraction(1, 10)):
    for R0 in Rs:
        d = EO.fraction_det(EO.rows_at(ctx, rows, tv, R0))
        pt["%s@R=%s" % (tv, R0)] = "nonzero" if d != 0 else "ZERO"
out["direct_exact_det"] = pt
print(pt, flush=True)
g = EO.common_factor(res_all)
out["gcd_degree"] = g.degree()
out["gcd_roots_(0,1)"] = int(g.count_roots(0, 1)) if g.degree() > 0 else 0
print("gcd", g.degree(), out["gcd_roots_(0,1)"], flush=True)
json.dump(out, open("debug/data/exch_reverse_review_det_%s.json" % case, "w"), indent=1)
