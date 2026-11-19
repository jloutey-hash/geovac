"""exch_reverse_review_closed.py -- REVIEW (cold) of debug/exch_reverse_memo.md.

Own evaluator of the FROZEN Q(t) closed form (tests/data, a = 1, b = t), written from the memo's definitions,
not from exch_reverse_core.  Coefficients: exact numerator/denominator integer polynomials in t, evaluated in
mpmath (so irrational or complex t is allowed).  K(R; a, b) = a K(aR; 1, b/a).

Two branches:
  real : LC = ln((1+t)/|t-1|), Ei = principal value (mpmath ei on the reals)
  cont : the holomorphic continuation FROM t > 1, evaluated at complex t (upper half-plane):
         LC = Log(1+t) - Log(t-1); Ei(w) for keys whose argument is > 0 at t > 1: gamma + Log(w) - Ein(-w);
         for keys whose argument is < 0 at t > 1: gamma + Log(-w) - Ein(-w) = -E1(-w).  Ein entire.
f0(R) = int_R^oo FM1(s)/s ds, FM1 = c (psi_a - psi_b + e^{2tR}Ei(-2(1+t)R) - e^{2R}Ei(-2(1+t)R)),
c = 1/(2(1 - t^2)).
"""
import json
import sys

import mpmath as mp
import sympy as sp

NAMES9 = ["f0", "-E1((a+b)R)", "e^{-(a+b)R}", "e^{-2aR}", "e^{-2bR}", "-e^{2aR}E1(2(a+b)R)", "-e^{2bR}E1(2(a+b)R)",
          "psi_a", "psi_b"]
TS = sp.Symbol("t")


def load(case):
    if case in ("2psig", "3dsig"):
        return json.load(open("tests/data/paper58_pachucki_ext_2026_11_18.json"))["sym"][case]
    return json.load(open("tests/data/paper58_exchange_m_2026_11_18.json"))["sym"][case]


def coeff_polys(case):
    out = []
    cft = load(case)
    for nm in NAMES9:
        d = {}
        for p, s in cft.get(nm, {}).items():
            e = sp.together(sp.sympify(s, locals={"t": TS}))
            n, dd = sp.fraction(e)
            d[int(p)] = ([int(c) for c in sp.Poly(n, TS).all_coeffs()],
                         [int(c) for c in sp.Poly(dd, TS).all_coeffs()])
        out.append(d)
    return out


def Ei_pos(w):   # continuation of Ei from positive w: gamma + Log w - Ein(-w), Ein(z) = E1(z) + Log z + gamma
    return -mp.e1(-w) + mp.log(w) - mp.log(-w)


def Ei_neg(w):   # continuation of Ei from negative w: gamma + Log(-w) - Ein(-w) = -E1(-w)
    return -mp.e1(-w)


def funcs(t, R, branch):
    """g1..g8 (indices 1..8 of NAMES9) and FM1 at (a, b) = (1, t)."""
    if branch == "real":
        LC = mp.log((1 + t) / abs(t - 1))
        E = mp.ei
        Ep = Em = E
    else:
        LC = mp.log(1 + t) - mp.log(t - 1)
        Ep, Em = Ei_pos, Ei_neg
    g = [None] * 9
    g[1] = Em(-(1 + t) * R)
    g[2] = mp.exp(-(1 + t) * R)
    g[3] = mp.exp(-2 * R)
    g[4] = mp.exp(-2 * t * R)
    E2 = Em(-2 * (1 + t) * R)
    g[5] = mp.exp(2 * R) * E2
    g[6] = mp.exp(2 * t * R) * E2
    # at t > 1: 1 - t < 0 (psi_a arguments negative), t - 1 > 0 (psi_b arguments positive)
    g[7] = mp.exp(-2 * R) * (Em(2 * (1 - t) * R) - 2 * Em((1 - t) * R) - LC)
    g[8] = mp.exp(-2 * t * R) * (Ep(2 * (t - 1) * R) - 2 * Ep((t - 1) * R) - LC)
    c = 1 / (2 * (1 - t * t))
    fm1 = c * (g[7] - g[8] + g[6] - g[5])
    return g, fm1


def f0(t, R, branch):
    f = lambda s: funcs(t, s, branch)[1] / s
    return mp.quad(f, [R, R + 1, R + 2, R + 4, R + 10, R + 30, mp.inf])


def K1(case_polys, t, R, branch):
    g, _ = funcs(t, R, branch)
    g[0] = f0(t, R, branch)
    s = 0
    for i, d in enumerate(case_polys):
        for p, (n, dd) in d.items():
            s += mp.polyval(n, t) / mp.polyval(dd, t) * R ** p * g[i]
    return s


def K(case_polys, a, b, R, branch="real"):
    return a * K1(case_polys, b / a, a * R, branch)


if __name__ == "__main__":
    mp.mp.dps = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    q = lambda s: mp.mpf(sp.Rational(s).p) / sp.Rational(s).q
    res = {}
    P = {c: coeff_polys(c) for c in ("2psig", "3dsig", "2ppi", "3ddel")}
    # 1. reproduce the memo's (2,1) values and a known a < b value
    for c, a, b, R in (("2psig", "2", "1", "7/5"), ("3dsig", "2", "1", "7/5"), ("2psig", "1", "2", "7/5"),
                       ("2ppi", "3/2", "4/5", "23/10")):
        v = K(P[c], q(a), q(b), q(R))
        res["repro_%s_%s_%s_%s" % (c, a, b, R)] = mp.nstr(v, 55)
        print(c, a, b, R, mp.nstr(v, 50), flush=True)
    # 2. new points for the hp comparison
    for c, a, b, R in (("3dsig", "5/2", "1", "7/5"), ("2ppi", "5/2", "1", "7/5"), ("2ppi", "1", "9/10", "2"),
                       ("3ddel", "1", "9/10", "7/5"), ("2psig", "5/2", "1", "2"), ("2ppi", "1", "99/100", "3/2"),
                       ("3dsig", "1", "9/10", "2")):
        v = K(P[c], q(a), q(b), q(R))
        res["new_%s_%s_%s_%s" % (c, a, b, R)] = mp.nstr(v, 55)
        print(c, a, b, R, mp.nstr(v, 50), flush=True)
    # 3. irrational ratio, a = sqrt 3, b = 1 (closed form only, real branch) and the continuation check
    for c in ("2ppi", "3dsig"):
        v = K(P[c], mp.sqrt(3), mp.mpf(1), mp.mpf(7) / 5)
        res["sqrt3_%s" % c] = mp.nstr(v, 55)
    # 4. analytic continuation from t > 1 through the upper half-plane: F(s + i eps) -> real branch at s < 1
    cont = {}
    for c in ("2psig", "3dsig", "2ppi", "3ddel"):
        for s in ("1/2", "1/3", "9/10"):
            sv, R = q(s), mp.mpf(7) / 5
            vr = K1(P[c], sv, R, "real")
            row = {"real": mp.nstr(vr, 40)}
            for e in (mp.mpf(10) ** -12, mp.mpf(10) ** -24):
                vc = K1(P[c], sv + 1j * e, R, "cont")
                row["eps%s" % mp.nstr(e, 2)] = [mp.nstr(abs(mp.re(vc) - vr), 3), mp.nstr(abs(mp.im(vc)), 3)]
            # same continuation through the LOWER half-plane must agree as well (branch-free)
            vl = K1(P[c], sv - 1j * mp.mpf(10) ** -24, R, "cont")
            row["lower_eps1e-24"] = [mp.nstr(abs(mp.re(vl) - vr), 3), mp.nstr(abs(mp.im(vl)), 3)]
            # sanity: at t > 1 the cont branch equals the real branch on the real axis
            vr2 = K1(P[c], q("2"), R, "real")
            vc2 = K1(P[c], q("2") + 1j * mp.mpf(10) ** -24, R, "cont")
            row["t2_check"] = mp.nstr(abs(vc2 - vr2), 3)
            cont["%s_%s" % (c, s)] = row
            print(c, s, row, flush=True)
    res["continuation"] = cont
    json.dump(res, open("debug/data/exch_reverse_review_closed.json", "w"), indent=1)
