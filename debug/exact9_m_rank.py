"""exact9_m_rank.py -- leg (D) of the exact-order-9 argument for the two m != 0 exchange integrals
2p_{+-1,A} 1s_B (CASE 2ppi) and 3d_{+-2,A} 1s_B (CASE 3ddel).

Algorithm 1 (adapted from debug/exact9_rank.py): coordinates of K, K', ..., K^(8) on the nine functions over
Q(t, R) (a = 1, b = t) with the hand-typed connection matrix; determinant at R = Rj as a polynomial in Z[t]
(sympy DomainMatrix over ZZ[t]); gcd over five Rj (python-flint, exact); real roots of the gcd in (1, oo).

Input: debug/data/pachucki_m_case_<CASE>_sym_res.json (nine_coords, strings in t, a = 1, b = t).
Checks: (i) hand connection matrix = tests/_exact_order.hand_table = helper calculus, identically over Q(t);
(ii) stored Q(t) closed form at t = 2 = the Q run debug/data/pachucki_m_case_<CASE>_1_2.json;
(iii) (4/5, 3/2) Q run = scaled Q(t) form at t = 15/8 (K(R;a,b) = a K(aR;1,b/a); f0 scale exponent solved for);
(iv) entry denominators of rows 0..8 in t.

Run: python debug/exact9_m_rank.py CASE      -> debug/data/exact9_m_rank_<CASE>.json
"""
import json
import sys
import time
from fractions import Fraction

import flint
import sympy as sp
from sympy.polys.matrices import DomainMatrix

sys.path.insert(0, "tests")
sys.path.insert(0, "debug")
import _pachucki_extension as PE  # noqa: E402
import _exact_order as EO  # noqa: E402
import pachucki_m_core  # noqa: E402,F401  (installs the exact flint fallback for sympy's heugcd)

case = sys.argv[1]
T0 = time.time()


def log(*a):
    print(*a, "[%.1fs]" % (time.time() - T0), flush=True)


t, R = sp.symbols("t R")
FF = sp.QQ.frac_field(t, R)
tF, RF = FF.from_sympy(t), FF.from_sympy(R)
one, zero = FF.one, FF.zero

src = "debug/data/pachucki_m_case_%s_sym_res.json" % case
coords = json.load(open(src))["nine_coords"]
assert len(coords) == 9
out = {"case": case, "source": src}

v0 = []
for cd in coords:
    acc = zero
    for p, s in cd.items():
        acc = acc + FF.from_sympy(sp.sympify(s, locals={"t": t})) * RF ** int(p)
    v0.append(acc)

# ---- connection matrix, typed by hand (a = 1, b = t), as in exact9_rank.py
c = one / (2 * (one - tF * tF))
G = [[zero] * 9 for _ in range(9)]
G[0][6], G[0][5], G[0][8], G[0][7] = -c / RF, c / RF, c / RF, -c / RF
G[1][2] = one / RF
G[2][2] = -(one + tF)
G[3][3] = -2 * one
G[4][4] = -2 * tF
G[5][5], G[5][4] = 2 * one, one / RF
G[6][6], G[6][3] = 2 * tF, one / RF
G[7][7], G[7][4], G[7][2] = -2 * one, one / RF, -2 / RF
G[8][8], G[8][3], G[8][2] = -2 * tF, one / RF, -2 / RF

# check (i): G = _exact_order.hand_table = helper calculus, identically over Q(t)
ctxT = PE.Ctx("Qt")
ht = EO.hand_table(ctxT)
hp = EO.helper_table(ctxT)
ok_i = True
for i in range(9):
    hrow = {j: {p: v for p, v in d.items() if v != 0} for j, d in ht[i].items()}
    hrow = {j: d for j, d in hrow.items() if d}
    if hp[i] is None or hrow != {j: {p: v for p, v in d.items() if v != 0} for j, d in hp[i].items()}:
        ok_i = False
    for j in range(9):
        mine = FF.to_sympy(G[i][j])
        theirs = sum(ctxT.K.to_sympy(v) * R ** p for p, v in ht[i].get(j, {}).items())
        if sp.simplify(mine - theirs) != 0:
            ok_i = False
            log("G mismatch", i, j)
out["connection_matrix_equals_hand_table_and_helper_identically_in_t"] = ok_i
log("connection matrix = _exact_order.hand_table = helper calculus over Q(t):", ok_i)


def dvec(v):
    o = [x.diff(RF) for x in v]
    for i in range(9):
        if v[i] != 0:
            for j in range(9):
                if G[i][j] != 0:
                    o[j] = o[j] + v[i] * G[i][j]
    return o


rows = [v0]
for _ in range(8):
    rows.append(dvec(rows[-1]))
log("derivative rows built")


def ev(x, tv=None, Rv=None):
    e = FF.to_sympy(x)
    sub = {}
    if tv is not None:
        sub[t] = sp.Rational(tv.numerator, tv.denominator)
    if Rv is not None:
        sub[R] = sp.Rational(Rv.numerator, Rv.denominator)
    return e.subs(sub)


# ---- check (ii): t = 2 vs the Q run at (1, 2)
q12 = json.load(open("debug/data/pachucki_m_case_%s_1_2.json" % case))["nine_coords"]
ok_ii = True
for j in range(9):
    theirs = sum(sp.Rational(s) * R ** int(p) for p, s in q12[j].items())
    if sp.expand(ev(v0[j], tv=Fraction(2)) - theirs) != 0:
        ok_ii = False
        log("t=2 mismatch in coordinate", j)
out["stored_Qt_at_t_2_equals_Q_run_1_2"] = ok_ii
log("Q(t) closed form at t = 2 = Q run (1,2):", ok_ii)

# ---- check (iii): (4/5, 3/2) Q run vs a^{1-s_i} c_i(aR; t = 15/8), s_i = 0 for i >= 1, s_0 solved
aq, tq = sp.Rational(4, 5), sp.Rational(15, 8)
q45 = json.load(open("debug/data/pachucki_m_case_%s_4-5_3-2.json" % case))["nine_coords"]
ok_iii = True
s0 = None
for j in range(9):
    theirs = sp.expand(sum(sp.Rational(s) * R ** int(p) for p, s in q45[j].items()))
    mine = sp.expand(ev(v0[j], tv=Fraction(15, 8)).subs(R, aq * R) * aq)
    if j == 0:
        ratio = sp.simplify(theirs / mine)
        k = sp.log(ratio, aq)
        s0 = sp.nsimplify(-k)   # theirs = a^{1 - s0} c0(aR) -> ratio = a^{-s0}
        ok_iii = ok_iii and ratio.is_Rational and aq ** (-s0) == ratio
    elif sp.expand(theirs - mine) != 0:
        ok_iii = False
        log("(4/5,3/2) mismatch in coordinate", j)
out["scaled_Qt_at_t_15_8_equals_Q_run_4-5_3-2"] = bool(ok_iii)
out["f0_scaling_exponent_s0 (f0(R;a,b) = a^s0 f0(aR;1,b/a))"] = str(s0)
log("(4/5,3/2) Q run = a K(aR;1,15/8):", ok_iii, " f0 scale exponent s0 =", s0)

# ---- check (iv): entry denominators of rows 0..8
den_f = {}
for rw in rows:
    for x in rw:
        for f, m in sp.factor_list(FF.to_sympy(x).as_numer_denom()[1])[1]:
            den_f[str(f)] = max(den_f.get(str(f), 0), int(m))
out["entry_denominator_factors_max_mult"] = den_f
other = sorted(set(den_f) - {"t", "t - 1", "t + 1", "R"})
out["entry_denominator_other_factors"] = other
log("entry denominators:", den_f, " other:", other)


def to_flint(p):
    return flint.fmpz_poly([int(cc) for cc in reversed(p.all_coeffs())])


def real_roots_gt1(fp):
    """(number of real roots, number certainly > 1, number undecided) -- arb root isolation (real roots exact-real)"""
    if fp.degree() <= 0:
        return 0, 0, 0
    rr = [x for x, _ in fp.complex_roots() if x.imag == 0]
    gt = sum(1 for x in rr if x.real > 1)
    und = sum(1 for x in rr if not (x.real > 1) and not (x.real < 1))
    return len(rr), gt, und


def strip(fp, lin):
    m = 0
    while True:
        q, r = divmod(fp, lin)
        if r != 0:
            return fp, m
        fp, m = q, m + 1


LT, LM, LP = flint.fmpz_poly([0, 1]), flint.fmpz_poly([-1, 1]), flint.fmpz_poly([1, 1])
ZT = sp.ZZ[t]
Rjs = [Fraction(1), Fraction(2), Fraction(7, 5), Fraction(1, 3), Fraction(5)]
nums, dets = [], {}
for Rj in Rjs:
    M = []
    for rw in rows:
        es = [sp.together(ev(x, Rv=Rj)) for x in rw]
        den = sp.lcm([sp.denom(e) for e in es])
        M.append([sp.Poly(sp.cancel(e * den), t, domain="ZZ") for e in es])
    dm = DomainMatrix([[ZT.from_sympy(p.as_expr()) for p in rw] for rw in M], (9, 9), ZT)
    d = sp.Poly(ZT.to_sympy(dm.det()), t)
    fp = to_flint(d)
    nums.append(fp)
    A, m0 = strip(fp, LT)
    A, m1 = strip(A, LM)
    A, mp1 = strip(A, LP)
    nrr, n_gt1, n_und = real_roots_gt1(A)
    dets[str(Rj)] = {"degree": fp.degree(), "zero": fp == 0, "mult_t": m0, "mult_t-1": m1, "mult_t+1": mp1,
                     "residual_degree": A.degree(), "residual_real_roots": nrr,
                     "residual_real_roots_in_(1,oo)": n_gt1, "residual_roots_undecided_vs_1": n_und}
    log("R =", Rj, dets[str(Rj)])
g = nums[0]
for fp in nums[1:]:
    g = g.gcd(fp)
G0, g0 = strip(g, LT)
G0, g1 = strip(G0, LM)
G0, gp1 = strip(G0, LP)
_, nroots, gund = real_roots_gt1(G0)
nroots += gund
out["det_at_Rj"] = dets
out["gcd_over_Rj"] = str(sp.factor(sp.Poly(list(reversed([int(x) for x in g.coeffs()])), t).as_expr()))
out["gcd_mult_t_tm1_tp1_residual_degree"] = [g0, g1, gp1, G0.degree()]
out["gcd_real_roots_in_(1,oo)"] = nroots
log("gcd over Rj:", out["gcd_over_Rj"], " real roots in (1,oo):", nroots)

# spot values (exact) at R = 7/5 -- illustration only, not the argument
spec = {}
for tv in (Fraction(2), Fraction(3), Fraction(3, 2), Fraction(5), Fraction(4, 3)):
    Mq = sp.Matrix(9, 9, lambda i, j: ev(rows[i][j], tv=tv, Rv=Fraction(7, 5)))
    spec[str(tv)] = bool(Mq.det() != 0)
out["det_nonzero_at_R_7/5_for_t"] = spec
log("det != 0 at R = 7/5:", spec)
out["verdict_rank9_for_every_real_t_gt_1"] = bool(nroots == 0 and ok_i and ok_ii and not other
                                                   and not any(v["zero"] for v in dets.values()))
out["time_s"] = round(time.time() - T0, 1)
json.dump(out, open("debug/data/exact9_m_rank_%s.json" % case, "w"), indent=1)
log("VERDICT rank 9 for every real t > 1:", out["verdict_rank9_for_every_real_t_gt_1"])
