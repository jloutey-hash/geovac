"""exch_reverse_det.py -- leg (D) on t in (0, 1): det of the 9 x 9 coordinate matrix of (K, ..., K^(8)) on the nine
functions is not identically zero in R for any real t0 in (0, 1).

Input: the frozen Q(t) closed form (a = 1, b = t) of CASE (tests/data).  The coordinates and the connection matrix are
rational in t with denominators t, t - 1, t + 1 only, so they specialise at every t0 in (0, 1), where (by the analytic
continuation / fresh-solve argument of the memo) they are the coordinates of K on the real-branch nine functions.

Algorithm A (tests/_exact_order.py): rows by the helper calculus (Calc.diff + in_nine) over Q(t); rows cleared by
t^i (t-1)^j (t+1)^k; det at rational R as an integer polynomial in t by fraction-free Bareiss at integer t + exact
Newton interpolation (checked on extra nodes).  R in {1, 2, 7/5, 1/3, 5}.
Algorithm B (independent of A apart from the input): rows by the HAND connection matrix (tests/_exact_order.hand_table,
typed from the definitions) on Laurent coefficients over Q(t); entries specialised at R; rows cleared by the lcm of
their t-denominators; det by sympy DomainMatrix over ZZ[t].  R in {3/2, 11/7, 9}.
For each R: strip t, t-1, t+1; count real roots of the residual in (0, 1) and in (1, oo); gcd of the residuals
(within A, within B, and across all eight R); real roots of every gcd in (0, 1).
Run: python debug/exch_reverse_det.py CASE  -> debug/data/exch_reverse_det_<CASE>.json"""
import json
import sys
import time
from fractions import Fraction

import sympy as sp
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, "tests")
sys.path.insert(0, "debug")
import _pachucki_extension as X  # noqa: E402
import _exact_order as EO  # noqa: E402
import _exchange_m_nonzero as M  # noqa: E402
import exch_reverse_core as C  # noqa: E402

M.install_gcd_fallback()
case = sys.argv[1]
T0 = time.time()
t = sp.Symbol("t")


def log(*a):
    print(*a, "[%.1fs]" % (time.time() - T0), flush=True)


def roots_in(P, lo, hi):
    """distinct real roots of sympy Poly P in the open interval (lo, hi) (hi may be oo)."""
    if P.degree() <= 0:
        return 0
    n = int(P.count_roots(lo, hi))
    for e in (lo, hi):
        if e != sp.oo and P.eval(e) == 0:
            n -= 1
    return n


def poly_of(asc):
    return sp.Poly(list(reversed(asc)), t, domain="ZZ")


ctx = X.Ctx("Qt")
cft = C.stored(case)["cft"]
K = EO.stored_K(ctx, cft)
out = {"case": case, "A": [], "B": []}

# ---------------- algorithm A
rows = EO.derivative_rows(ctx, K, 8)
per_row, other = EO.denominator_audit(rows)
assert not other, other
out["A_row_clearing"] = per_row
cleared = EO.clear_rows(ctx, rows, per_row)
log("A rows built; clearing", per_row)
resA = []
for R0 in (Fraction(1), Fraction(2), Fraction(7, 5), Fraction(1, 3), Fraction(5)):
    N, scale = EO.determinant_polynomial(cleared, R0)
    mult, res = EO.strip_linear(N)
    P = poly_of(res)
    row = {"R": str(R0), "degree": len(N) - 1, "mult_t_t-1_t+1": mult, "residual_degree": P.degree(),
           "residual_roots_(0,1)": roots_in(P, 0, 1), "residual_roots_(1,oo)": roots_in(P, 1, sp.oo)}
    out["A"].append(row)
    resA.append(res)
    log("A", row)

# ---------------- algorithm B (hand connection matrix)
FF = ctx.K
G = EO.hand_table(ctx)          # row i: {j: {p: c}}: g_i' = sum_j sum_p c R^p g_j
co0 = X.in_nine(ctx, K)
assert co0 is not None


def dvec(v):
    """d/dR of sum_i v_i(R) g_i, coordinates as Laurent dicts."""
    w = [dict() for _ in range(9)]
    for i, d in enumerate(v):
        for p, c in d.items():
            if p != 0:
                w[i][p - 1] = w[i].get(p - 1, FF.zero) + p * c
            for j, dj in G[i].items():
                for q, cq in dj.items():
                    w[j][p + q] = w[j].get(p + q, FF.zero) + c * cq
    return [{p: c for p, c in d.items() if c != 0} for d in w]


vs = [[dict(d) for d in co0]]
for _ in range(8):
    vs.append(dvec(vs[-1]))
# cross-check: B rows equal A rows identically
same = all({p: c for p, c in vb[i].items() if c != 0} == {p: c for p, c in ra[i].items() if c != 0}
           for vb, ra in zip(vs, rows) for i in range(9))
out["hand_rows_equal_helper_rows"] = same
log("B rows built; equal to helper rows:", same)
resB = []
for R0 in (Fraction(3, 2), Fraction(11, 7), Fraction(9)):
    Rq = FF(R0.numerator) / FF(R0.denominator)
    Mrows = []
    for v in vs:
        ent = []
        for d in v:
            acc = FF.zero
            for p, c in d.items():
                acc = acc + c * Rq ** p
            ent.append(acc)
        dens = [sp.Poly(FF.to_sympy(e.denom), t) if e != 0 else sp.Poly(1, t) for e in ent]
        L = dens[0]
        for dd in dens[1:]:
            L = sp.lcm(L, dd)
        prow = []
        for e in ent:
            if e == 0:
                prow.append(sp.Poly(0, t, domain="QQ"))
            else:
                y = sp.cancel(FF.to_sympy(e) * L.as_expr())
                prow.append(sp.Poly(y, t, domain="QQ"))
        cden = 1
        for pp in prow:
            for cc in pp.all_coeffs():
                cden = sp.ilcm(cden, sp.Rational(cc).q)
        Mrows.append([sp.Poly(pp * cden, t, domain="ZZ") for pp in prow])
    DM = DomainMatrix([[sp.ZZ[t].from_sympy(p.as_expr()) for p in r] for r in Mrows], (9, 9), sp.ZZ[t])
    det = DM.det()
    Pd = sp.Poly(sp.ZZ[t].to_sympy(det), t, domain="ZZ")
    asc = [int(c) for c in reversed(Pd.all_coeffs())]
    mult, res = EO.strip_linear(asc)
    P = poly_of(res)
    row = {"R": str(R0), "degree": Pd.degree(), "mult_t_t-1_t+1": mult, "residual_degree": P.degree(),
           "residual_roots_(0,1)": roots_in(P, 0, 1), "residual_roots_(1,oo)": roots_in(P, 1, sp.oo)}
    out["B"].append(row)
    resB.append(res)
    log("B", row)

# ---------------- gcds
for nm, rs in (("A", resA), ("B", resB), ("all", resA + resB)):
    g = EO.common_factor(rs)
    out["gcd_" + nm] = {"degree": g.degree(), "roots_(0,1)": roots_in(g, 0, 1), "roots_(1,oo)": roots_in(g, 1, sp.oo),
                        "poly": str(g.as_expr()) if g.degree() < 40 else "deg %d" % g.degree()}
    log("gcd", nm, out["gcd_" + nm])
# a pair-wise check too (each pair of R already removes every residual root)
pairs = []
allr = resA + resB
for i in range(len(allr)):
    for j in range(i + 1, len(allr)):
        pairs.append(EO.common_factor([allr[i], allr[j]]).degree())
out["pairwise_gcd_degrees_max"] = max(pairs)
out["time_s"] = round(time.time() - T0, 1)
log("pairwise gcd degree max:", max(pairs))
json.dump(out, open("debug/data/exch_reverse_det_%s.json" % case, "w"), indent=1)
