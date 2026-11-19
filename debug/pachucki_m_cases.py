"""m != 0 two-centre exchange integrals through r12-weighted masters (pachucki_m_core), exact.

Per electron: orbital product = N_A N_B (x + i y)^s e^{-a r_A - b r_B}  (s = 1: 2p_{+-1,A} 1s_B;  s = 2: 3d_{+-2,A} 1s_B,
and also 2p_{+1,A} conj(2p_{-1,B}), the sigma = 2 part of real 2p_x 2p_x exchange -- the SAME density).
Angular reduction (X = x1 x2 + y1 y2, rho^2 = x^2 + y^2 = r_A^2 - z^2, z = (r_A^2 + R^2 - r_B^2)/(2R)):
    Re (x1+iy1)(x2-iy2)       = X = (r1A^2 + r2A^2 - r12^2)/2 - z1 z2
    Re (x1+iy1)^2 (x2-iy2)^2  = 2 X^2 - rho1^2 rho2^2
Monomial r12^{2k} r1A^n1 r1B^n2 r2A^n3 r2B^n4 / r12 (all n even)  ->  (4 pi)^2 / R * d^alpha f_{2k},  alpha = n + 1.
Cases (value conventions = geovac.qfd_core.exchange_hp with the orbital tuples of the committed tests):
  2ppi   (2p+1_A 1s_B | 2p-1_A 1s_B)  = - (a^5 b^3/pi^2) (1/2) int int X psi psi / r12   [= -(2p_x 1s | 2p_x 1s)]
  3ddel  (3d+2_A 1s_B | 3d-2_A 1s_B)  =   (a^7 b^3/(12 pi^2)) int int (2X^2 - rho1^2 rho2^2) psi psi / r12
  1s     control through the Laplacian identity is in pachucki_m_stage1.py; here '1s' is the plain f_{1111}.
Run: python debug/pachucki_m_cases.py CASE (A B | sym) [--op ORDER DEGREE] [--num]
"""
from __future__ import annotations

import json
import pickle
import sys
import time
from fractions import Fraction

import sympy as sp

sys.path.insert(0, "debug")
sys.path.insert(0, "tests")
import pachucki_ext_core as E  # noqa: E402
import pachucki_m_core as M  # noqa: E402
import _pachucki_exchange as PX  # noqa: E402

case = sys.argv[1]
SYM = sys.argv[2] == "sym"
t0 = time.time()
if SYM:
    ctx = E.Ctx("Qt")
    tg = f"{case}_sym"
    A_, B_ = sp.Integer(1), sp.Symbol("t")
else:
    a, b = Fraction(sys.argv[2]), Fraction(sys.argv[3])
    ctx = E.Ctx("Q", a, b)
    tg = f"{case}_{a}_{b}".replace("/", "-")
    A_, B_ = sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator)

s1A, s1B, s2A, s2B, v, Rs = sp.symbols("s1A s1B s2A s2B v R", positive=True)
z1 = (s1A + Rs**2 - s1B) / (2 * Rs)
z2 = (s2A + Rs**2 - s2B) / (2 * Rs)
X = (s1A + s2A - v) / 2 - z1 * z2
SPEC = {
    "1s": (sp.Integer(1), A_**3 * B_**3),
    "2ppi": (X, -A_**5 * B_**3 / 2),
    "3ddel": (2 * X**2 - (s1A - z1**2) * (s2A - z2**2), A_**7 * B_**3 / 12),
}
Pol, N2pi2 = SPEC[case]
poly = sp.Poly(sp.expand(Pol), s1A, s1B, s2A, s2B, v)
terms = {}                       # k -> list of (alpha, Laurent coefficient in R)
for (i1, j1, i2, j2, k), co in poly.terms():
    terms.setdefault(int(k), []).append(((2 * i1 + 1, 2 * j1 + 1, 2 * i2 + 1, 2 * j2 + 1), co))
kmax = max(terms)
# multi-index sets: f_{2k} needs f_{2k-2} (and f_{2k-4}) on the same set
Dk = {}
gens = []
for k in range(kmax, -1, -1):
    gens = gens + [al for al, _ in terms.get(k, [])]
    Dk[k] = E.downset(gens)
print(case, tg, {k: len(d) for k, d in Dk.items()}, flush=True)

calc = M.MCalc(ctx)
if "--from-symres" in sys.argv:
    # jets solved (with residual checks) by the resumable driver pachucki_m_symres.py
    import importlib.util
    raw = pickle.load(open(f"debug/data/pachucki_m_symres_{case}.pkl", "rb"))

    def dec(Fn):
        o = {}
        for key, d in Fn.items():
            kk = tuple(x if (i == 0 or (key[0] == "E" and i == 2)) else ctx.K.from_sympy(sp.sympify(x))
                       for i, x in enumerate(key))
            o[kk] = {p: ctx.K.from_sympy(sp.sympify(c)) for p, c in d.items()}
        return o
    fl = {k: {al: dec(Fn) for al, Fn in d.items()} for k, d in raw.items() if k <= kmax}
    assert all(set(fl[k]) == set(Dk[k]) for k in Dk), "symres state incomplete"
    nfree = None
    tg = tg + "_res"
    print("jets loaded from the resumable run", {k: len(d) for k, d in fl.items()}, flush=True)
else:
    calc0, f0, Fs, P2, nfree = E.solve_all(ctx, Dk[0], check=True)
    print("f_0 jets: %d multi-indices, residuals identically 0, free unknowns %d  (%.0fs)" % (len(Dk[0]), nfree, time.time() - t0),
          flush=True)
    fl = {0: f0}
    for k in range(1, kmax + 1):
        RJ = M.RJets(ctx, Dk[k], calc)
        fl[k] = M.solve_fn(calc, RJ, Dk[k], 2 * k - 1, fl[k - 1], fl[k - 2] if k >= 2 else None, check=True)
        print("f_%d jets: %d multi-indices, residuals identically 0  (%.0fs)" % (2 * k, len(Dk[k]), time.time() - t0), flush=True)

if SYM:
    K_ = ctx.K
    pref = K_.from_sympy(sp.together(16 * N2pi2))
else:
    pref = Fraction(str(16 * N2pi2))


def cfield(x):
    x = sp.Rational(x)
    return ctx.c(Fraction(int(x.p), int(x.q)))


K = {}
for k, lst in terms.items():
    for al, cR in lst:
        for term in sp.Add.make_args(sp.expand(cR)):
            cc, e = term.as_coeff_exponent(Rs)
            K = E.fadd(K, E.fscale(fl[k][al], pref * cfield(cc), int(e) - 1))
keys = sorted(K, key=str)
cK = M.in_nine(ctx, K)
print("keys of K: %d; K in nine-function span: %s" % (len(keys), cK is not None), flush=True)
out = {"case": case, "tag": tg, "sets": {str(k): len(d) for k, d in Dk.items()}, "residuals_identically_zero": True,
       "free_unknowns_f0": nfree, "keys_of_K": [str(k) for k in keys], "K_in_nine_span": cK is not None,
       "laurent_range_of_K": [min(min(d) for d in K.values()), max(max(d) for d in K.values())]}
if not SYM:
    pickle.dump({"K": K, "nine_coords": cK}, open(f"debug/data/pachucki_m_K_{tg}.pkl", "wb"))

ders = [K]
for _ in range(10):
    ders.append(calc.diff(ders[-1]))
coords = [M.in_nine(ctx, A) for A in ders]
out["derivatives_in_nine_span"] = all(c is not None for c in coords)

if SYM:
    t = ctx.t
    dens = set()
    for c in cK:
        for p, x in c.items():
            d = sp.factor_list(ctx.K.to_sympy(x).as_numer_denom()[1])
            for fct, _ in d[1]:
                dens.add(str(fct))
    out["denominator_factors_of_K_coefficients"] = sorted(dens)
    maxpole = 0
    for c in cK:
        for p, x in c.items():
            d = sp.Poly(ctx.K.to_sympy(x).as_numer_denom()[1], t)
            mlt = 0
            while d.eval(1) == 0:
                d = sp.Poly(sp.cancel(d.as_expr() / (t - 1)), t)
                mlt += 1
            maxpole = max(maxpole, mlt)
    out["max_pole_order_at_t_1"] = maxpole
    out["nine_coords"] = [{str(p): str(ctx.K.to_sympy(x)) for p, x in c.items()} for c in cK]
    print("denominator factors:", sorted(dens), " max pole at t=1:", maxpole, flush=True)
    # rank at (t, R) = (3, 7/5)
    Mx = [[sum((Fraction(str(ctx.K.to_sympy(x).subs(t, 3))) * Fraction(7, 5) ** p for p, x in c.items()), Fraction(0))
           for c in row] for row in coords]
    out["rank_at_t3_R7_5"] = PX.rank_q(Mx)
    print("rank of K..K^(10) coordinates at (t,R)=(3,7/5):", out["rank_at_t3_R7_5"], flush=True)
else:
    Rv = Fraction(7, 5)
    Mx = [[sum((x * Rv ** p for p, x in c.items()), Fraction(0)) for c in row] for row in coords]
    out["rank_K_to_K10"] = PX.rank_q(Mx)
    out["nine_coords"] = [{str(p): str(x) for p, x in c.items()} for c in cK]
    print("rank K..K^(10) at R = 7/5:", out["rank_K_to_K10"], flush=True)
    if "--op" in sys.argv:
        i = sys.argv.index("--op")
        order, degree = int(sys.argv[i + 1]), int(sys.argv[i + 2])
        t1 = time.time()
        ops = PX.operators_annihilating(ders, order, degree)
        low = PX.operators_annihilating(ders, order, degree - 1)
        lowo = PX.operators_annihilating(ders, order - 1, degree + 12)
        print("operators (order<=%d, deg<=%d): %d; (deg<=%d): %d; (order<=%d, deg<=%d): %d   (%.0fs)"
              % (order, degree, len(ops), degree - 1, len(low), order - 1, degree + 12, len(lowo), time.time() - t1), flush=True)
        out["n_operators"] = {"order": order, "degree": degree, "count": len(ops), "count_degree_minus_1": len(low),
                              "count_order_minus_1_degree_plus_12": len(lowo)}
        if len(ops) == 1:
            op = [PX.trim(c) for c in ops[0]]
            out["derived_operator"] = [[str(x) for x in c] for c in op]
            out["derived_coeff_degrees"] = [len(c) - 1 for c in op]
    if "--num" in sys.argv:
        import mpmath as mp
        mp.mp.dps = 60
        alg = PX.Algebra(a, b)
        vals = {}
        for Rv_ in sys.argv[sys.argv.index("--num") + 1].split(","):
            x = mp.mpf(Rv_)
            f0n = PX.f0_numeric(alg, x)
            vals[Rv_] = mp.nstr(PX.evaluate(alg, K, x, f0n), 55)
            print("  closed form K(%s) = %s" % (Rv_, vals[Rv_]), flush=True)
        out["closed_form_values_dps60"] = vals
out["time_s"] = round(time.time() - t0, 1)
json.dump(out, open(f"debug/data/pachucki_m_case_{tg}.json", "w"), indent=1)
print("done %.0fs" % (time.time() - t0), flush=True)
