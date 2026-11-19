"""Review of debug/pachucki_m_nonzero_memo.md: plant wrong coefficients in (**) and see which control catches them.
Also records the number of free unknowns left by every solve_D2 (uniqueness of the decaying solution) and evaluates
the closed form at a new pair.

(**)  q f_{m+1}'' = c f_{m+1} + m^2 (d^4 - Sig d^2 + e) f_{m-1} + m (m-1)^2 (m-2) f_{m-3}'' + S_m

Run: python debug/pachucki_m_review_plants.py A B R PLANT[,PLANT...] [--cases 2ppi,3ddel]
Writes debug/data/pachucki_m_review_plants_<A>_<B>.json  (one record per plant)."""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction

import mpmath as mp
import sympy as sp

sys.path.insert(0, "debug")
sys.path.insert(0, "tests")
import pachucki_ext_core as E  # noqa: E402
import pachucki_m_core as M  # noqa: E402
import _pachucki_exchange as PX  # noqa: E402

FREE = []


def solve_fn_plant(calc, RJ, D, m, flow, flow3, plant):
    ctx = calc.ctx
    S = RJ.source(m)
    out = {}
    q0 = RJ.q[(0, 0, 0, 0)]
    m2 = ctx.c(m) if plant == "m2_to_m" else ctx.c(m * m)
    k3 = {"drop_f_m3": 0, "f_m3_coef_6": m * (m - 1) * (m - 2)}.get(plant, m * (m - 1) ** 2 * (m - 2))
    ssgn = -1 if plant == "S_sign" else 1
    csgn = -1 if plant == "c_sign" else 1
    esgn = 0 if plant == "drop_e" else 1
    sigsgn = -1 if plant == "Sig_sign" else 1
    d4 = 0 if plant == "drop_d4" else 1
    for al in D:
        rhs = E.fscale(dict(S[al]), ctx.c(ssgn))
        for U in M.subs_of(al):
            T = tuple(s - u for s, u in zip(al, U))
            bn = ctx.c(M.binom(al, U))
            if sum(U) and RJ.q.get(U, 0) != 0:
                rhs = E.fadd(rhs, E.fscale(calc.D2(out[T]), bn * RJ.q[U]), -1)
            if plant == "c_highU_dropped" and sum(U) >= 4:
                pass
            elif sum(U) and RJ.c.get(U, 0) != 0:
                rhs = E.fadd(rhs, E.fscale(out[T], csgn * bn * RJ.c[U]))
            lo = flow[T]
            if sum(U) == 0 and d4:
                rhs = E.fadd(rhs, E.fscale(calc.D2(calc.D2(lo)), m2))
            if RJ.Sig.get(U, 0) != 0:
                rhs = E.fadd(rhs, E.fscale(calc.D2(lo), -sigsgn * m2 * bn * RJ.Sig[U]))
            if esgn and RJ.e.get(U, 0) != 0:
                rhs = E.fadd(rhs, E.fscale(lo, m2 * bn * RJ.e[U]))
            if sum(U) == 0 and m >= 3 and k3:
                rhs = E.fadd(rhs, E.fscale(calc.D2(flow3[T]), ctx.c(k3)))
        out[al] = calc.solve_D2(E.fscale(rhs, 1 / q0))
        FREE.append(calc.last_free)
        assert E.fadd(E.fscale(calc.D2(out[al]), q0), rhs, -1) == {}, al
    return out


def controls(ctx, plant):
    """stage-1 controls with the plant in place: label 61 (n=1, n=3) at the exchange point, control (a), (a')."""
    a, b = ctx.a, ctx.b
    C_ = ctx.c
    D = E.downset([(2, 2, 1, 1)])
    calc0, fS, Fs, P2, nfree = E.solve_all(ctx, D, check=True)
    calc = M.MCalc(ctx)
    RJ = M.RJets(ctx, D, calc)
    Z = (0, 0, 0, 0)
    u, w = a, b
    Wm2 = {("E", -2 * u, 0): {-1: -(-1 / (w - u) - 1 / (u + w))},
           ("E", -2 * w, 0): {-1: -(-1 / (u - w) - 1 / (u + w))},
           ("E", -(u + w), 0): {0: -C_(2), -1: -2 / (u + w)},
           ("Ei", ctx.zero, -(u + w)): {0: C_(2)}}
    f2 = solve_fn_plant(calc, RJ, D, 1, fS, None, plant)
    lab61 = E.fscale(E.fadd(E.fadd(E.fscale(fS[Z], 2 * (u * u + w * w)), calc.D2(fS[Z]), -1), Wm2), 1 / (u * u - w * w) ** 2)
    r = {"label61_n1": E.fadd(lab61, f2[Z], -1) == {}}
    lhs = E.fscale(fS[(1, 1, 1, 1)], C_(16))
    rhs = {}
    for al, co, sh in (((1, 1, 1, 1), a * a + b * b, 0), ((0, 1, 1, 1), 2 * a, 0), ((1, 0, 1, 1), 2 * b, 0),
                       ((2, 0, 1, 1), a * b, 0), ((0, 2, 1, 1), a * b, 0), ((0, 0, 1, 1), -a * b, 2)):
        rhs = E.fadd(rhs, E.fscale(f2[al], 8 * co, sh))
    r["control_a"] = E.fadd(lhs, rhs, -1) == {}
    f4 = solve_fn_plant(calc, RJ, D, 3, f2, fS, plant)
    # label 61 at n = 3 with HAND coefficients 18, -9, -12, source from label 63 is NOT available for n=3 by hand;
    # use the unplanted S_3 (so an S-sign plant is visible here only through f_2/f_4 consistency).
    S3 = RJ.source(3)
    lhs61 = E.fscale(calc.D2(f4[Z]), (u * u - w * w) ** 2)
    rhs61 = E.fadd(E.fadd(E.fadd(E.fscale(calc.D2(f2[Z]), 18 * (u * u + w * w)), E.fscale(calc.D2(calc.D2(f2[Z])), C_(9)), -1),
                          E.fscale(calc.D2(fS[Z]), C_(12)), -1), S3[Z], -1)
    r["label61_n3"] = E.fadd(lhs61, rhs61, -1) == {}
    lhs = E.fscale(f2[(1, 1, 1, 1)], C_(12))
    rhs = {}
    for al, co, sh in (((1, 1, 1, 1), a * a + b * b, 0), ((0, 1, 1, 1), 2 * a, 0), ((1, 0, 1, 1), 2 * b, 0),
                       ((2, 0, 1, 1), a * b, 0), ((0, 2, 1, 1), a * b, 0), ((0, 0, 1, 1), -a * b, 2)):
        rhs = E.fadd(rhs, E.fscale(f4[al], co, sh))
    r["control_a_prime_f4"] = E.fadd(lhs, rhs, -1) == {}
    return r


def case_value(ctx, case, Rnums, plant):
    A_, B_ = sp.Rational(ctx.a.numerator, ctx.a.denominator), sp.Rational(ctx.b.numerator, ctx.b.denominator)
    s1A, s1B, s2A, s2B, v, Rs = sp.symbols("s1A s1B s2A s2B v R", positive=True)
    z1 = (s1A + Rs**2 - s1B) / (2 * Rs)
    z2 = (s2A + Rs**2 - s2B) / (2 * Rs)
    X = (s1A + s2A - v) / 2 - z1 * z2
    SPEC = {"2ppi": (X, -A_**5 * B_**3 / 2),
            "3ddel": (2 * X**2 - (s1A - z1**2) * (s2A - z2**2), A_**7 * B_**3 / 12)}
    Pol, N2 = SPEC[case]
    poly = sp.Poly(sp.expand(Pol), s1A, s1B, s2A, s2B, v)
    terms = {}
    for (i1, j1, i2, j2, k), co in poly.terms():
        terms.setdefault(int(k), []).append(((2 * i1 + 1, 2 * j1 + 1, 2 * i2 + 1, 2 * j2 + 1), co))
    kmax = max(terms)
    Dk, gens = {}, []
    for k in range(kmax, -1, -1):
        gens = gens + [al for al, _ in terms.get(k, [])]
        Dk[k] = E.downset(gens)
    calc = M.MCalc(ctx)
    _, f0, _, _, _ = E.solve_all(ctx, Dk[0], check=True)
    fl = {0: f0}
    for k in range(1, kmax + 1):
        RJ = M.RJets(ctx, Dk[k], calc)
        fl[k] = solve_fn_plant(calc, RJ, Dk[k], 2 * k - 1, fl[k - 1], fl[k - 2] if k >= 2 else None, plant)
    pref = Fraction(str(16 * N2))
    K = {}
    for k, lst in terms.items():
        for al, cR in lst:
            for term in sp.Add.make_args(sp.expand(cR)):
                cc, e = term.as_coeff_exponent(Rs)
                cc = sp.Rational(cc)
                K = E.fadd(K, E.fscale(fl[k][al], pref * Fraction(int(cc.p), int(cc.q)), int(e) - 1))
    alg = PX.Algebra(ctx.a, ctx.b)
    vals = {}
    for Rv in Rnums:
        x = mp.mpf(Rv)
        vals[Rv] = mp.nstr(PX.evaluate(alg, K, x, PX.f0_numeric(alg, x)), 50)
    return vals, M.in_nine(ctx, K) is not None


if __name__ == "__main__":
    mp.mp.dps = 60
    a, b = Fraction(sys.argv[1]), Fraction(sys.argv[2])
    Rnums = sys.argv[3].split(",")
    plants = sys.argv[4].split(",")
    cases = sys.argv[sys.argv.index("--cases") + 1].split(",") if "--cases" in sys.argv else []
    ctx = E.Ctx("Q", a, b)
    tg = f"{a}_{b}".replace("/", "-")
    fn = f"debug/data/pachucki_m_review_plants_{tg}.json"
    try:
        res = json.load(open(fn))
    except Exception:
        res = {}
    for pl in plants:
        t0 = time.time()
        FREE.clear()
        rec = {"plant": pl}
        if "--nocontrols" not in sys.argv:
            rec.update(controls(ctx, pl))
        for cs in cases:
            vals, in9 = case_value(ctx, cs, Rnums, pl)
            rec[cs] = vals
            rec[cs + "_in_nine"] = in9
        rec["max_free_unknowns_solve_D2"] = max(FREE) if FREE else None
        rec["n_solves"] = len(FREE)
        rec["time_s"] = round(time.time() - t0, 1)
        print(rec, flush=True)
        res[pl] = rec
        json.dump(res, open(fn, "w"), indent=1)
