"""exact9_review_det.py -- adversarial re-computation of leg (D) of debug/exact9_independence_memo.md.

Independent of debug/exact9_rank.py in three ways:
  * the derivative rows come from the helper's own calculus over Q(t) (Calc.diff + in_nine), not from a
    hand-written connection matrix; the hand matrix of the memo's table is separately compared IDENTICALLY in t;
  * the determinant numerators are obtained by evaluation at integer t modulo a prime and Newton interpolation
    (pure Python), not by sympy DomainMatrix over Z[t];
  * R values 3/2, 11/7, 9 (none of the memo's five).
Also: fresh solves over Q at t not equal to 2 against the stored coefficients, and exact determinants from those
fresh solves compared with the stored-coefficient route.

Run: python debug/exact9_review_det.py CASE   (1s 2psig 3dsig ppsig pppi0)
"""
import json
import sys
import time
from fractions import Fraction

import sympy as sp

sys.path.insert(0, "tests")
import _pachucki_extension as PE  # noqa: E402

case = sys.argv[1]
T0 = time.time()
P = (1 << 61) - 1
out = {"case": case, "prime": P}


def log(*a):
    print(*a, "[%.1fs]" % (time.time() - T0), flush=True)


ctx = PE.Ctx("Qt")
tS = ctx.t
KK = ctx.K
calc = PE.Calc(ctx)
nine = PE.nine_functions(ctx)

# ---------------------------------------------------------------- 1. memo's derivative table, identically in t
tF = ctx.b
one = ctx.one
c = one / (2 * (one - tF * tF))
# table of the memo (a = 1, b = t): row i -> {j: {power of R: coeff}}
TAB = [
    {5: {-1: c}, 6: {-1: -c}, 7: {-1: -c}, 8: {-1: c}},
    {2: {-1: one}},
    {2: {0: -(one + tF)}},
    {3: {0: -2 * one}},
    {4: {0: -2 * tF}},
    {5: {0: 2 * one}, 4: {-1: one}},
    {6: {0: 2 * tF}, 3: {-1: one}},
    {7: {0: -2 * one}, 4: {-1: one}, 2: {-1: -2 * one}},
    {8: {0: -2 * tF}, 3: {-1: one}, 2: {-1: -2 * one}},
]
ok_tab = True
for i, (g, _) in enumerate(nine):
    co = PE.in_nine(ctx, calc.diff(g))
    if co is None:
        ok_tab = False
        continue
    for j in range(9):
        mine = {p: v for p, v in TAB[i].get(j, {}).items() if v != 0}
        theirs = {p: v for p, v in co[j].items() if v != 0}
        if mine != theirs:
            ok_tab = False
            print("TABLE MISMATCH", i, j, mine, theirs)
out["memo_derivative_table_equals_helper_calculus_identically_in_t"] = ok_tab
log("memo table = helper calculus over Q(t):", ok_tab)

# ---------------------------------------------------------------- 2. K from the stored coefficients, as a Func
src = ("debug/data/pachucki_ext_symbolic_1s.json" if case == "1s" else "debug/data/pachucki_ext_sym_%s.json" % case)
cft = json.load(open(src))["closed_form_coefficients_t"]
K: dict = {}
for (g, _), nm in zip(nine, PE.NAMES9):
    for p, s in cft.get(nm, {}).items():
        K = PE.fadd(K, PE.fscale(g, KK.from_sympy(sp.sympify(s, locals={"t": tS})), int(p)))
ders = calc.derivatives(K, 8)
rows = []
for d in ders:
    co = PE.in_nine(ctx, d)
    assert co is not None
    rows.append(co)
log("derivative rows via helper calculus built")

# ---------------------------------------------------------------- 3. fresh solves over Q at other t
fresh = {}
TV = [Fraction(7, 3), Fraction(3), Fraction(17, 10)] if case in ("1s", "2psig") else [Fraction(7, 3)]
fresh_rows = {}
for tv in TV:
    cq = PE.Ctx("Q", 1, tv)
    calq, fS, Fs, P2, nfree = PE.solve_all(cq, PE.case_downset(case))
    Kq = PE.assemble_K(cq, case, fS)
    co = PE.in_nine(cq, Kq)
    stored = PE.coeff_strings_at(cft, tv)
    same = co is not None and all({p: v for p, v in co[j].items() if v != 0} == stored[j] for j in range(9))
    fresh[str(tv)] = bool(same)
    dq = calq.derivatives(Kq, 8)
    fresh_rows[tv] = [PE.in_nine(cq, d) for d in dq]
    log("fresh Q solve at t =", tv, "equals stored:", same)
out["stored_equals_fresh_Q_solve_at_t"] = fresh


# ---------------------------------------------------------------- 4. compile the rows for modular evaluation
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
EX = {k: sum(r[k] for r in row_exp) for k in ("t", "t - 1", "t + 1")}
out["clearing_exponents_total"] = EX
log("denominators: other factors", sorted(other_den), "clearing", EX)
comp = [[{p: compile_elt(x) for p, x in cdict.items()} for cdict in rw] for rw in rows]
degbound = 0
for rw in rows:
    degbound += max((sp.Poly(KK.to_sympy(KK(x.numer)), tS).degree() for cdict in rw for x in cdict.values()),
                    default=0) + 0
NPTS = degbound + EX["t"] + EX["t - 1"] + EX["t + 1"] + 40
out["interpolation_points"] = NPTS


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
    """coefficients (ascending) of the interpolating polynomial mod P"""
    n = len(xs)
    dd = ys[:]
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) * pow(xs[i] - xs[i - j], -1, P) % P
    poly = [0]
    for i in range(n - 1, -1, -1):
        new = [0] * (len(poly) + 1)                    # poly * (x - xs[i]) + dd[i]
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
resid = {}
info = {}
for R0 in R0S:
    xs = list(range(2, 2 + NPTS))
    ys = [cleared_det_at(x, R0) for x in xs]
    poly = newton(xs, ys)
    chk = all(horner(poly, x) == cleared_det_at(x, R0) for x in range(2 + NPTS, 2 + NPTS + 25))
    A, m0 = strip(poly, [0, 1])
    A, m1 = strip(A, [P - 1, 1])
    A, mm1 = strip(A, [1, 1])
    resid[R0] = A
    info[str(R0)] = {"degree_cleared": len(poly) - 1, "extra_points_agree": bool(chk), "mult_t": m0,
                     "mult_t-1": m1, "mult_t+1": mm1, "residual_degree": len(A) - 1}
    log("R0 =", R0, info[str(R0)])
out["det_numerators_mod_p"] = info
gg = {}
for i in range(len(R0S)):
    for j in range(i + 1, len(R0S)):
        g = pgcd(resid[R0S[i]], resid[R0S[j]])
        gg["%s|%s" % (R0S[i], R0S[j])] = len(g) - 1
out["gcd_degree_of_residual_factors_mod_p"] = gg
log("gcd degrees of residual factors:", gg)


# ---------------------------------------------------------------- 5. exact rational determinants
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


cmp_ = {}
for tv, fr in fresh_rows.items():
    for R0 in R0S[:2]:
        Mf = [[sum((Fraction(x) * R0 ** p for p, x in cdict.items()), Fraction(0)) for cdict in rw] for rw in fr]
        d1, d2 = det_frac(Mf), stored_det(tv, R0)
        cmp_["t=%s,R=%s" % (tv, R0)] = {"equal": bool(d1 == d2), "nonzero": bool(d1 != 0), "float": float(d1)}
out["det_fresh_solve_vs_stored_route"] = cmp_
log("fresh-vs-stored determinants:", {k: (v["equal"], v["nonzero"]) for k, v in cmp_.items()})

# sign scan of det(t, R0) on t > 1 (exact rationals)
scan = {}
grid = [Fraction(1) + Fraction(k, 200) for k in range(1, 200)] + [Fraction(2) + Fraction(k, 20) for k in range(0, 161)] \
    + [Fraction(k) for k in (11, 13, 17, 25, 40, 77, 150, 400, 1000, 10 ** 4, 10 ** 6)]
for R0 in R0S[:2]:
    signs = set()
    mn = None
    for tv in grid:
        d = stored_det(tv, R0)
        signs.add((d > 0) - (d < 0))
    scan[str(R0)] = {"points": len(grid), "signs_seen": sorted(signs)}
    log("sign scan R0 =", R0, scan[str(R0)])
out["sign_scan_t_in_(1,1e6]"] = scan
out["time_s"] = round(time.time() - T0, 1)
json.dump(out, open("debug/data/exact9_review_det_%s.json" % case, "w"), indent=1)
log("done")
