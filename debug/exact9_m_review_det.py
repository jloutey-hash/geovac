"""exact9_m_review_det.py -- algorithm 2 for leg (D), m != 0 cases (2ppi, 3ddel); adapted from
debug/exact9_review_det.py.  Independent of debug/exact9_m_rank.py:
  * derivative rows from the helper calculus over Q(t) (tests/_pachucki_extension Calc.diff + in_nine), not from
    the hand connection matrix;
  * determinant numerators by evaluation at integer t modulo the prime 2^61 - 1 and Newton interpolation (pure
    Python), not sympy DomainMatrix over Z[t];
  * R values 3/2, 11/7, 9 (none of algorithm 1's five).
After stripping t, t - 1, t + 1, the pairwise gcd of the residual factors mod p is computed.  Degree 0 for a pair
means no common root at all of the two residuals (hence no t0 > 1 with det(t0, R) == 0 identically).

Run: python debug/exact9_m_review_det.py CASE   -> debug/data/exact9_m_review_det_<CASE>.json
"""
import json
import sys
import time
from fractions import Fraction

import sympy as sp

sys.path.insert(0, "tests")
sys.path.insert(0, "debug")
import _pachucki_extension as PE  # noqa: E402
import pachucki_m_core  # noqa: E402,F401  (exact flint fallback for sympy's heugcd)

case = sys.argv[1]
T0 = time.time()
P = (1 << 61) - 1
out = {"case": case, "prime": P}


def log(*a):
    print(*a, "[%.1fs]" % (time.time() - T0), flush=True)


ctx = PE.Ctx("Qt")
tS, KK = ctx.t, ctx.K
calc = PE.Calc(ctx)
nine = PE.nine_functions(ctx)

coords = json.load(open("debug/data/pachucki_m_case_%s_sym_res.json" % case))["nine_coords"]
K: dict = {}
for (g, _), cd in zip(nine, coords):
    for p, s in cd.items():
        K = PE.fadd(K, PE.fscale(g, KK.from_sympy(sp.sympify(s, locals={"t": tS})), int(p)))
ders = calc.derivatives(K, 8)
rows = []
for d in ders:
    co = PE.in_nine(ctx, d)
    assert co is not None, "derivative leaves the nine-function span"
    rows.append(co)
# K re-expanded on the nine equals the stored coordinates
back = PE.in_nine(ctx, K)
ok_back = all({int(p): KK.from_sympy(sp.sympify(s, locals={"t": tS})) for p, s in coords[j].items()
               if sp.sympify(s) != 0} == {p: v for p, v in back[j].items() if v != 0} for j in range(9))
out["stored_coords_roundtrip"] = bool(ok_back)
log("derivative rows via helper calculus built; round trip:", ok_back)


def compile_elt(x):
    def pl(poly):
        d = {}
        for (e,), cc in poly.terms():
            cc = Fraction(int(cc.numerator), int(cc.denominator))
            d[e] = cc.numerator % P * pow(cc.denominator, -1, P) % P
        return [d.get(e, 0) for e in range(max(d) + 1)] if d else [0]
    return pl(x.numer), pl(x.denom)


def horner(cs, x):
    acc = 0
    for cc in reversed(cs):
        acc = (acc * x + cc) % P
    return acc


other_den = set()
row_exp = []
for rw in rows:
    ex = {"t": 0, "t - 1": 0, "t + 1": 0}
    for cdict in rw:
        for x in cdict.values():
            den = KK.to_sympy(KK(x.denom))
            for f, m in sp.factor_list(den)[1]:
                s = str(f)
                if s in ex:
                    ex[s] = max(ex[s], int(m))
                else:
                    other_den.add(s)
    row_exp.append(ex)
out["denominator_factors_other_than_t_tm1_tp1"] = sorted(other_den)
out["row_clearing_exponents"] = row_exp
EX = {k: sum(r[k] for r in row_exp) for k in ("t", "t - 1", "t + 1")}
out["clearing_exponents_total"] = EX
log("denominators: other factors", sorted(other_den), "clearing", EX)
comp = [[{p: compile_elt(x) for p, x in cdict.items()} for cdict in rw] for rw in rows]
degbound = 0
for rw, ex in zip(rows, row_exp):
    mx = 0
    for cdict in rw:
        for x in cdict.values():
            dn = sp.Poly(KK.to_sympy(KK(x.numer)), tS).degree()
            dd = sp.Poly(KK.to_sympy(KK(x.denom)), tS).degree()
            mx = max(mx, dn - dd + ex["t"] + ex["t - 1"] + ex["t + 1"])
    degbound += max(mx, 0)
NPTS = degbound + 40
out["interpolation_points"] = NPTS
log("degree bound", degbound, "points", NPTS)


def det_mod(M):
    M = [r[:] for r in M]
    n = len(M)
    d = 1
    for i in range(n):
        pr = next((k for k in range(i, n) if M[k][i]), None)
        if pr is None:
            return 0
        if pr != i:
            M[i], M[pr] = M[pr], M[i]
            d = -d
        d = d * M[i][i] % P
        inv = pow(M[i][i], -1, P)
        for k in range(i + 1, n):
            if M[k][i]:
                f = M[k][i] * inv % P
                M[k] = [(x - f * y) % P for x, y in zip(M[k], M[i])]
    return d % P


def cleared_det_at(tk, R0):
    rm = R0.numerator % P * pow(R0.denominator, -1, P) % P
    M = []
    for rw in comp:
        r_ = []
        for cdict in rw:
            acc = 0
            for p, (nu, de) in cdict.items():
                acc = (acc + horner(nu, tk) * pow(horner(de, tk), -1, P) % P * pow(rm, p, P)) % P
            r_.append(acc)
        M.append(r_)
    d = det_mod(M)
    return d * pow(tk, EX["t"], P) % P * pow(tk - 1, EX["t - 1"], P) % P * pow(tk + 1, EX["t + 1"], P) % P


def newton(xs, ys):
    n = len(xs)
    dd = ys[:]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) * pow(xs[i] - xs[i - j], -1, P) % P
    poly = [0]
    for i in range(n - 1, -1, -1):
        new = [0] * (len(poly) + 1)
        for k, cc in enumerate(poly):
            new[k + 1] = (new[k + 1] + cc) % P
            new[k] = (new[k] - cc * xs[i]) % P
        new[0] = (new[0] + dd[i]) % P
        poly = new
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def pdivmod(A, B):
    A = A[:]
    q = [0] * max(1, len(A) - len(B) + 1)
    inv = pow(B[-1], -1, P)
    while len(A) >= len(B) and any(A):
        f = A[-1] * inv % P
        sh = len(A) - len(B)
        q[sh] = f
        for i, cc in enumerate(B):
            A[i + sh] = (A[i + sh] - f * cc) % P
        while len(A) > 1 and A[-1] == 0:
            A.pop()
        if A == [0]:
            break
    return q, A


def pgcd(A, B):
    while B != [0]:
        A, B = B, pdivmod(A, B)[1]
    inv = pow(A[-1], -1, P)
    return [x * inv % P for x in A]


def strip(A, lin):
    m = 0
    while True:
        q, r_ = pdivmod(A, lin)
        if r_ != [0]:
            return A, m
        A = q
        while len(A) > 1 and A[-1] == 0:
            A.pop()
        m += 1


R0S = [Fraction(3, 2), Fraction(11, 7), Fraction(9)]
resid, info = {}, {}
for R0 in R0S:
    xs = list(range(2, 2 + NPTS))
    ys = [cleared_det_at(x, R0) for x in xs]
    poly = newton(xs, ys)
    chk = all(horner(poly, x) == cleared_det_at(x, R0) for x in range(2 + NPTS, 2 + NPTS + 25))
    A, m0 = strip(poly, [0, 1])
    A, m1 = strip(A, [P - 1, 1])
    A, mm1 = strip(A, [1, 1])
    resid[R0] = A
    info[str(R0)] = {"degree_cleared": len(poly) - 1, "zero": poly == [0], "extra_points_agree": bool(chk),
                     "degree_below_point_count": len(poly) - 1 < NPTS - 30,
                     "mult_t": m0, "mult_t-1": m1, "mult_t+1": mm1, "residual_degree": len(A) - 1}
    log("R0 =", R0, info[str(R0)])
out["det_numerators_mod_p"] = info
gg = {}
for i in range(len(R0S)):
    for j in range(i + 1, len(R0S)):
        g = pgcd(resid[R0S[i]], resid[R0S[j]])
        gg["%s|%s" % (R0S[i], R0S[j])] = len(g) - 1
out["gcd_degree_of_residual_factors_mod_p"] = gg
log("gcd degrees of residual factors:", gg)


def det_frac(M):
    M = [r[:] for r in M]
    n = len(M)
    d = Fraction(1)
    for i in range(n):
        pr = next((k for k in range(i, n) if M[k][i] != 0), None)
        if pr is None:
            return Fraction(0)
        if pr != i:
            M[i], M[pr] = M[pr], M[i]
            d = -d
        d *= M[i][i]
        for k in range(i + 1, n):
            if M[k][i] != 0:
                f = M[k][i] / M[i][i]
                M[k] = [x - f * y for x, y in zip(M[k], M[i])]
    return d


def stored_det(tv, R0):
    M = [[sum((PE.qt_to_fraction(ctx, x, tv) * R0 ** p for p, x in cdict.items()), Fraction(0)) for cdict in rw]
         for rw in rows]
    return det_frac(M)


# exact det at one (t, R) also seen by algorithm 1 (R = 7/5, t = 3): for the cross-algorithm sanity print
out["exact_det_t3_R7_5_float"] = float(stored_det(Fraction(3), Fraction(7, 5)))
scan = {}
grid = [Fraction(1) + Fraction(k, 100) for k in range(1, 100)] + [Fraction(2) + Fraction(k, 10) for k in range(0, 81)] \
    + [Fraction(k) for k in (11, 13, 17, 25, 40, 77, 150, 400, 1000, 10 ** 4, 10 ** 6)]
for R0 in R0S[:2]:
    signs = set()
    for tv in grid:
        d = stored_det(tv, R0)
        signs.add((d > 0) - (d < 0))
    scan[str(R0)] = {"points": len(grid), "signs_seen": sorted(signs)}
    log("sign scan R0 =", R0, scan[str(R0)])
out["sign_scan_t_in_(1,1e6]"] = scan
out["verdict_no_common_root"] = bool(all(v == 0 for v in gg.values()) and not other_den
                                     and all(v["extra_points_agree"] and not v["zero"] for v in info.values()))
out["time_s"] = round(time.time() - T0, 1)
json.dump(out, open("debug/data/exact9_m_review_det_%s.json" % case, "w"), indent=1)
log("VERDICT (algorithm 2) residual factors pairwise coprime mod p:", out["verdict_no_common_root"])
