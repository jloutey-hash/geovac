"""PART B over Q(t): the m = 0 exchange densities of debug/pachucki_ext_lgt0.py with a = 1, b = t symbolic.
Writes the closed form's coefficients on the nine functions (same format as pachucki_ext_symbolic_1s.json) to
debug/data/pachucki_ext_sym_<CASE>.json, after checking every residual is identically 0 over Q(t) and that K lies in
the nine-function span identically in t.  Run: python debug/pachucki_ext_lgt0_sym.py CASE
"""
import json
import sys
import time

import sympy as sp

sys.path.insert(0, "debug")
import pachucki_ext_core as E  # noqa: E402

case = sys.argv[1]
t0 = time.time()
ctx = E.Ctx("Qt")
KK, tS = ctx.K, ctx.t
rA2, rB2, Rs = sp.symbols("rA2 rB2 R", positive=True)
zA = (rA2 + Rs ** 2 - rB2) / (2 * Rs)
zB = zA - Rs
SPEC = {   # (P, N2 * pi^2) with a = 1, b = t
    "1s": (sp.Integer(1), tS ** 3),
    "2psig": (zA, tS ** 3),
    "3dsig": (3 * zA ** 2 - rA2, tS ** 3 / 18),
    "ppsig": (zA * zB, tS ** 5),
    "pppi0": (rA2 - zA ** 2, tS ** 5 / 4),
}
P, N2pi2 = SPEC[case]
P = sp.expand(P)
mono = {}
for term in sp.Add.make_args(P):
    c, rest = term.as_independent(rA2, rB2)
    i = sp.degree(rest, rA2) if rest.has(rA2) else 0
    j = sp.degree(rest, rB2) if rest.has(rB2) else 0
    mono[(int(i), int(j))] = mono.get((int(i), int(j)), 0) + c
terms = [((2 * i1 + 1, 2 * j1 + 1, 2 * i2 + 1, 2 * j2 + 1), sp.expand(c1 * c2))
         for (i1, j1), c1 in mono.items() for (i2, j2), c2 in mono.items()]
D = E.downset([al for al, _ in terms])
print(case, "over Q(t); multi-indices:", len(D), flush=True)
calc, fS, Fs, P2, nfree = E.solve_all(ctx, D, check=True)
print("solved; residuals identically 0; free unknowns", nfree, "%.0fs" % (time.time() - t0), flush=True)
pref = KK.from_sympy(16 * N2pi2)
K = {}
for al, cR in terms:
    for term in sp.Add.make_args(cR):
        cc, k = term.as_coeff_exponent(Rs)
        K = E.fadd(K, E.fscale(fS[al], pref * KK.from_sympy(cc), int(k) - 1))
a, b, one = ctx.a, ctx.b, ctx.one


def psi(x, y):
    return {("Ei", -2 * x, 2 * (x - y)): {0: one}, ("Ei", -2 * x, (x - y)): {0: -2 * one}, ("E", -2 * x, 1): {0: -one}}


nine = [{("f",): {0: one}}, {("Ei", ctx.zero, -(a + b)): {0: one}}, {("E", -(a + b), 0): {0: one}},
        {("E", -2 * a, 0): {0: one}}, {("E", -2 * b, 0): {0: one}},
        {("Ei", 2 * a, -2 * (a + b)): {0: one}}, {("Ei", 2 * b, -2 * (a + b)): {0: one}}, psi(a, b), psi(b, a)]
lead = [next(iter(g)) for g in nine[:7]] + [("Ei", -2 * a, 2 * (a - b)), ("Ei", -2 * b, 2 * (b - a))]
names = ["f0", "-E1((a+b)R)", "e^{-(a+b)R}", "e^{-2aR}", "e^{-2bR}", "-e^{2aR}E1(2(a+b)R)", "-e^{2bR}E1(2(a+b)R)",
         "psi_a", "psi_b"]
rest = {k: dict(v) for k, v in K.items()}
co = []
for g, l in zip(nine, lead):
    c = rest.get(l, {})
    co.append(dict(c))
    for p, x in c.items():
        rest = E.fadd(rest, E.fscale(g, x, p), -1)
in_span = rest == {}
print("K in the nine-function span over Q(t):", in_span)
dens = set()
cft = {}
for nm, c in zip(names, co):
    cft[nm] = {}
    for p, x in c.items():
        e = sp.factor(KK.to_sympy(x))
        cft[nm][p] = str(e)
        for f, m in sp.factor_list(sp.denom(e))[1]:
            dens.add(str(f))
print("denominator factors in t:", dens)
json.dump({"case": case, "free_unknowns": nfree, "residuals_identically_zero": True, "K_in_nine_span_over_Qt": in_span,
           "closed_form_coefficients_t": cft, "denominator_factors_in_t": sorted(dens),
           "time_s": round(time.time() - t0, 1)}, open(f"debug/data/pachucki_ext_sym_{case}.json", "w"), indent=1)
print("done %.0fs" % (time.time() - t0))
