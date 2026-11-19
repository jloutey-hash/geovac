"""exact9_rank.py -- for which real t = b/a > 1 does the coordinate matrix of (K, K', ..., K^(8)) on the nine
functions have rank 9 over C(R)?

Input: the stored closed-form coefficients of K on the nine functions (a = 1, b = t), strings in t per power of R
(debug/data/pachucki_ext_symbolic_1s.json, debug/data/pachucki_ext_sym_<case>.json).
The connection matrix GAM (d/dR of the nine functions in the nine functions) is written here by hand from the
definitions, not taken from the helper; it is compared with the helper's Calc.diff at t = 2 and t = 3, and the
stored coefficients are compared with a fresh solve over Q at (1, 2) by the helper.

Bad t0: det M(t0, R) identically 0 in R.  If t0 > 1 is bad then det M(t0, Rj) = 0 for every rational Rj (the
entries have denominators only R^k t^i (t-1)^j (t+1)^k', checked), so t0 is a common root of the numerators
N_j(t) = numer det M(t, Rj).  We compute g = gcd_j N_j and count its real roots in (1, oo) (Sturm).

Run: python debug/exact9_rank.py CASE      (CASE in 1s 2psig 3dsig ppsig pppi0)
Writes debug/data/exact9_rank_<CASE>.json
"""
import json
import sys
import time
from fractions import Fraction

import sympy as sp
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, "tests")
import _pachucki_extension as PE  # noqa: E402

case = sys.argv[1]
T0 = time.time()
t, R = sp.symbols("t R")
FF = sp.QQ.frac_field(t, R)
tF, RF = FF.from_sympy(t), FF.from_sympy(R)
one, zero = FF.one, FF.zero
NAMES = PE.NAMES9

src = ("debug/data/pachucki_ext_symbolic_1s.json" if case == "1s" else "debug/data/pachucki_ext_sym_%s.json" % case)
cft = json.load(open(src))["closed_form_coefficients_t"]

# ---- v0: coordinates of K on the nine functions, elements of Q(t, R)
v0 = []
for nm in NAMES:
    acc = zero
    for p, s in cft.get(nm, {}).items():
        acc = acc + FF.from_sympy(sp.sympify(s, locals={"t": t})) * RF ** int(p)
    v0.append(acc)

# ---- connection matrix, by hand (a = 1, b = t); row i = coordinates of g_i'
# g1 f0; g2 Ei(-(1+t)R); g3 e^{-(1+t)R}; g4 e^{-2R}; g5 e^{-2tR}; g6 e^{2R}Ei(-2(1+t)R); g7 e^{2tR}Ei(-2(1+t)R);
# g8 psi_a = e^{-2R}[Ei(2(1-t)R) - 2Ei((1-t)R) - LC]; g9 psi_b = e^{-2tR}[Ei(2(t-1)R) - 2Ei((t-1)R) - LC]
c = one / (2 * (one - tF * tF))          # label 60 prefactor 1/(2(a^2-b^2))
G = [[zero] * 9 for _ in range(9)]
# f0' = -FM1/R, FM1 = c (g7 - g6 - g9 + g8)
G[0][6], G[0][5], G[0][8], G[0][7] = -c / RF, c / RF, c / RF, -c / RF
G[1][2] = one / RF                                   # Ei(-(1+t)R)' = e^{-(1+t)R}/R
G[2][2] = -(one + tF)
G[3][3] = -2 * one
G[4][4] = -2 * tF
G[5][5], G[5][4] = 2 * one, one / RF                 # e^{2R} e^{-2(1+t)R}/R = e^{-2tR}/R
G[6][6], G[6][3] = 2 * tF, one / RF
G[7][7], G[7][4], G[7][2] = -2 * one, one / RF, -2 / RF      # e^{-2R}[e^{2(1-t)R} - 2 e^{(1-t)R}]/R
G[8][8], G[8][3], G[8][2] = -2 * tF, one / RF, -2 / RF


def dvec(v):
    out = [x.diff(RF) for x in v]
    for i in range(9):
        if v[i] != 0:
            for j in range(9):
                if G[i][j] != 0:
                    out[j] = out[j] + v[i] * G[i][j]
    return out


rows = [v0]
for _ in range(9):
    rows.append(dvec(rows[-1]))
print(case, "derivative rows built %.1fs" % (time.time() - T0), flush=True)


def ev(x, tv=None, Rv=None):
    """evaluate an element of Q(t,R) at t = tv and/or R = Rv -> sympy expr"""
    e = FF.to_sympy(x)
    sub = {}
    if tv is not None:
        sub[t] = sp.Rational(tv.numerator, tv.denominator)
    if Rv is not None:
        sub[R] = sp.Rational(Rv.numerator, Rv.denominator)
    return e.subs(sub)


out = {"case": case}

# ---- check 1: GAM against the helper's calculus at t = 2, 3, 5/2 (helper Func algebra over Q)
ok_gam = True
for tv in (Fraction(2), Fraction(3), Fraction(5, 2)):
    ctx = PE.Ctx("Q", 1, tv)
    calc = PE.Calc(ctx)
    nine = PE.nine_functions(ctx)
    for i, (g, _) in enumerate(nine):
        co = PE.in_nine(ctx, calc.diff(g))
        assert co is not None, ("derivative leaves the span", i)
        for j in range(9):
            mine = sp.Poly(sp.together(ev(G[i][j], tv=tv) * R ** 3), R)   # Laurent -> poly
            theirs = sum(sp.Rational(cc.numerator, cc.denominator) * R ** (p + 3) for p, cc in co[j].items())
            if sp.expand(mine.as_expr() - theirs) != 0:
                ok_gam = False
                print("GAM mismatch", tv, i, j)
out["connection_matrix_equals_helper_at_t_2_3_5/2"] = ok_gam
print("connection matrix = helper calculus:", ok_gam, flush=True)

# ---- check 2: stored coefficients = fresh helper solve over Q at (1, 2)
tv = Fraction(2)
ctx = PE.Ctx("Q", 1, tv)
calc, fS, Fs, P2, nfree = PE.solve_all(ctx, PE.case_downset(case))
Kq = PE.assemble_K(ctx, case, fS)
co = PE.in_nine(ctx, Kq)
ok_v0 = co is not None
if ok_v0:
    for j in range(9):
        theirs = sum(sp.Rational(cc.numerator, cc.denominator) * R ** p for p, cc in co[j].items())
        if sp.simplify(ev(v0[j], tv=tv) - theirs) != 0:
            ok_v0 = False
out["stored_coefficients_equal_fresh_Q_solve_at_t_2"] = ok_v0
print("stored closed form = fresh solve at (1,2):", ok_v0, "%.1fs" % (time.time() - T0), flush=True)

# ---- denominators of all entries of rows 0..8
den_f = {}
for rw in rows[:9]:
    for x in rw:
        for f, m in sp.factor_list(FF.to_sympy(x).as_numer_denom()[1])[1]:
            den_f[str(f)] = max(den_f.get(str(f), 0), int(m))
out["entry_denominator_factors"] = den_f
print("entry denominators:", den_f, flush=True)

# ---- determinant at R = Rj as a rational function of t
ZT = sp.ZZ[t]
Rjs = [Fraction(1), Fraction(2), Fraction(7, 5), Fraction(1, 3), Fraction(5)]
nums = []
dets = {}
for Rj in Rjs:
    M = []
    for rw in rows[:9]:
        es = [sp.together(ev(x, Rv=Rj)) for x in rw]
        den = sp.lcm([sp.denom(e) for e in es])
        M.append([sp.Poly(sp.cancel(e * den), t, domain="ZZ") for e in es])
    dm = DomainMatrix([[ZT.from_sympy(p.as_expr()) for p in rw] for rw in M], (9, 9), ZT)
    d = sp.Poly(ZT.to_sympy(dm.det()), t)
    nums.append(d)
    fl = sp.factor_list(d.as_expr())
    small = [(str(f), int(m)) for f, m in fl[1] if sp.Poly(f, t).degree() <= 2]
    big = [(sp.Poly(f, t).degree(), int(m)) for f, m in fl[1] if sp.Poly(f, t).degree() > 2]
    dets[str(Rj)] = {"degree": d.degree(), "zero": bool(d.is_zero), "small_factors": small,
                     "other_factor_degrees": big,
                     "real_roots_in_(1,oo)": int(sp.Poly(d, t).count_roots(1, sp.oo)) - int(d.eval(1) == 0)}
    print("R =", Rj, "det numerator degree", d.degree(), "small factors", small, "others (deg, mult)", big,
          "%.1fs" % (time.time() - T0), flush=True)
g = nums[0]
for d in nums[1:]:
    g = sp.gcd(g, d)
g = sp.Poly(g, t)
out["det_at_Rj"] = dets
out["gcd_over_Rj"] = str(sp.factor(g.as_expr()))
nroots = int(g.count_roots(1, sp.oo)) - int(g.eval(1) == 0) if g.degree() > 0 else 0
out["gcd_real_roots_in_(1,oo)"] = nroots
print("gcd over R_j:", sp.factor(g.as_expr()), "; real roots in (1, oo):", nroots, flush=True)

# ---- explicit determinant values at the t with coincidences (all at R = 7/5), over Q
spec = {}
for tv in (Fraction(2), Fraction(3), Fraction(3, 2), Fraction(5), Fraction(4, 3)):
    Mq = sp.Matrix(9, 9, lambda i, j: ev(rows[i][j], tv=tv, Rv=Fraction(7, 5)))
    spec[str(tv)] = bool(Mq.det() != 0)
out["det_nonzero_at_R_7/5_for_t"] = spec
print("det != 0 at R = 7/5:", spec, flush=True)
out["verdict_rank9_for_every_real_t_gt_1"] = bool(nroots == 0 and ok_gam and ok_v0
                                                   and set(den_f) <= {"t", "t - 1", "t + 1", "R"})
out["time_s"] = round(time.time() - T0, 1)
json.dump(out, open("debug/data/exact9_rank_%s.json" % case, "w"), indent=1)
print("VERDICT rank 9 for every real t > 1:", out["verdict_rank9_for_every_real_t_gt_1"], flush=True)
