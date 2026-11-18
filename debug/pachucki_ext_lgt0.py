"""PART B: m = 0 exchange densities (polynomials in r_A^2, r_B^2, R) through Pachucki's equation, exact over Q.

Map.  rho(r) = N_A N_B P(r_A^2, r_B^2; R) e^{-a r_A - b r_B} per electron;  r_X^{2i} e^{-u r_X} = (-1)^{2i+1}... i.e.
      r^n e^{-u r} = (-1)^{n+1} d_u^{n+1} [e^{-u r}/r],   so a monomial r1A^n1 r1B^n2 r2A^n3 r2B^n4 maps to
      (-1)^{n1+n2+n3+n4} d^alpha f,  alpha = (n1+1, n2+1, n3+1, n4+1) in the exponent order (u3, u2, w2, w3),
      and K = N2 (4 pi)^2 / R * sum (coefficient) (-1)^{|n|} f_alpha   (f of arXiv:0905.1635 label 01).
Cases (a = rate of the A orbital, b of the B orbital; z_A = (r_A^2 + R^2 - r_B^2)/(2R), z_B = z_A - R):
  2psig  (2p0_A 1s_B | 2p0_A 1s_B)    P = z_A               N2 = a^5 b^3 / pi^2
  3dsig  (3d0_A 1s_B | 3d0_A 1s_B)    P = 3 z_A^2 - r_A^2   N2 = a^7 b^3 / (18 pi^2)
  ppsig  (2p0_A 2p0_B | 2p0_A 2p0_B)  P = z_A z_B           N2 = a^5 b^5 / pi^2
  pppi0  (2p+1_A 2p+1_B | same)       P = r_A^2 - z_A^2     N2 = a^5 b^5 / (4 pi^2)   (density (x^2+y^2), m cancels)
  1s     control                      P = 1                 N2 = a^3 b^3 / pi^2
Run: python debug/pachucki_ext_lgt0.py CASE A B [--op] [--num]
"""
from __future__ import annotations

import json
import os
import sys
import time
from fractions import Fraction

import sympy as sp

sys.path.insert(0, "debug")
sys.path.insert(0, "tests")
import pachucki_ext_core as E  # noqa: E402
import _pachucki_exchange as PX  # noqa: E402

case = sys.argv[1]
a, b = Fraction(sys.argv[2]), Fraction(sys.argv[3])
DO_OP, DO_NUM = "--op" in sys.argv, "--num" in sys.argv
t0 = time.time()
rA2, rB2, Rs = sp.symbols("rA2 rB2 R", positive=True)
zA = (rA2 + Rs ** 2 - rB2) / (2 * Rs)
zB = zA - Rs
A_, B_ = sp.Rational(a.numerator, a.denominator), sp.Rational(b.numerator, b.denominator)
SPEC = {
    "1s": (sp.Integer(1), A_ ** 3 * B_ ** 3),
    "2psig": (zA, A_ ** 5 * B_ ** 3),
    "3dsig": (3 * zA ** 2 - rA2, A_ ** 7 * B_ ** 3 / 18),
    "ppsig": (zA * zB, A_ ** 5 * B_ ** 5),
    "pppi0": (rA2 - zA ** 2, A_ ** 5 * B_ ** 5 / 4),
}
P, N2pi2 = SPEC[case]                      # N2 = N2pi2 / pi^2
P = sp.expand(P)
# per-electron monomials: coefficient (Laurent in R) of rA2^i rB2^j
mono = {}
for term in sp.Add.make_args(P):
    c, rest = term.as_independent(rA2, rB2)
    i = sp.degree(rest, rA2) if rest.has(rA2) else 0
    j = sp.degree(rest, rB2) if rest.has(rB2) else 0
    mono[(int(i), int(j))] = mono.get((int(i), int(j)), 0) + c
terms = []
for (i1, j1), c1 in mono.items():
    for (i2, j2), c2 in mono.items():
        al = (2 * i1 + 1, 2 * j1 + 1, 2 * i2 + 1, 2 * j2 + 1)      # (-1)^{|n|} = +1: all n even
        terms.append((al, sp.expand(c1 * c2)))
D = E.downset([al for al, _ in terms])
print(case, (a, b), "multi-indices:", len(D), flush=True)
ctx = E.Ctx("Q", a, b)
calc, fS, Fs, P2, nfree = E.solve_all(ctx, D, check=True)
tsolve = time.time() - t0
print("all f_alpha solved, residuals identically 0; free unknowns %d; %.1fs" % (nfree, tsolve), flush=True)
K = {}
pref = Fraction(str(16 * N2pi2))           # N2 (4 pi)^2 = 16 N2pi2
for al, cR in terms:
    for term in sp.Add.make_args(cR):
        cc, k = term.as_coeff_exponent(Rs)
        K = E.fadd(K, E.fscale(fS[al], pref * Fraction(str(cc)), int(k) - 1))
keys = sorted(K, key=str)
print("keys of K:", len(keys), keys)

# nine-function span
one = Fraction(1)


def psi(x, y):
    return {("Ei", -2 * x, 2 * (x - y)): {0: one}, ("Ei", -2 * x, (x - y)): {0: -2 * one}, ("E", -2 * x, 1): {0: -one}}


nine = [{("f",): {0: one}}, {("Ei", Fraction(0), -(a + b)): {0: one}}, {("E", -(a + b), 0): {0: one}},
        {("E", -2 * a, 0): {0: one}}, {("E", -2 * b, 0): {0: one}},
        {("Ei", 2 * a, -2 * (a + b)): {0: one}}, {("Ei", 2 * b, -2 * (a + b)): {0: one}}, psi(a, b), psi(b, a)]
lead = [next(iter(g)) if len(g) == 1 else ("Ei", -2 * x, 2 * (x - y)) for g, (x, y) in
        zip(nine, [(0, 0)] * 7 + [(a, b), (b, a)])]


def in_nine(A):
    rest = {k: dict(v) for k, v in A.items()}
    co = []
    for g, l in zip(nine, lead):
        c = rest.get(l, {})
        co.append(dict(c))
        for p, x in c.items():
            rest = E.fadd(rest, E.fscale(g, x, p), -1)
    return co if rest == {} else None


cK = in_nine(K)
ders = [K]
for _ in range(10):
    ders.append(calc.diff(ders[-1]))
coords = [in_nine(A) for A in ders]
Rv = Fraction(7, 5)
rank = None
if all(c is not None for c in coords):
    M = [[sum((x * Rv ** p for p, x in c.items()), Fraction(0)) for c in row] for row in coords]
    rank = PX.rank_q(M)
print("K in nine-function span:", cK is not None, " rank K..K^(10) at R = 7/5:", rank, flush=True)
tg = f"{case}_{a}_{b}".replace("/", "-")
out = {"case": case, "pair": [str(a), str(b)], "multi_indices": len(D), "free_unknowns": nfree,
       "residuals_identically_zero": True, "solve_time_s": round(tsolve, 1),
       "keys_of_K": [str(k) for k in keys], "K_in_nine_span": cK is not None, "rank_K_to_K10": rank,
       "laurent_range_of_K": [min(min(d) for d in K.values()), max(max(d) for d in K.values())]}

COMMITTED = {("2psig", "1", "2"): "ep_heldout_2psig_1_2_EPS.json", ("ppsig", "1", "2"): "eq_heldout_ppsig_1_2_EQS.json",
             ("pppi0", "1", "2"): "eq_heldout_pppi0_1_2_EQS.json", ("1s", "1", "2"): None}
if DO_OP:
    fn = COMMITTED.get((case, str(a), str(b)))
    Rsym = sp.Symbol("R")
    if fn:
        J = json.load(open(os.path.join("debug", "data", fn)))
        cop = [sp.Poly(sp.sympify(s_, locals={"R": Rsym}), Rsym) for s_ in J["operator_Q_crt"]]
        order, degree = J["order"], J["degree"]
    else:
        cop, order = None, 9
        degree = int(os.environ.get("PX_DEG", "16"))
    t1 = time.time()
    ops = PX.operators_annihilating(ders, order, degree)
    print("operators of order <= %d, degree <= %d annihilating the closed form: %d (%.0fs)"
          % (order, degree, len(ops), time.time() - t1), flush=True)
    out["n_operators_order_degree"] = [order, degree, len(ops)]
    ops_lower = PX.operators_annihilating(ders, order, degree - 1)
    out["n_operators_degree_minus_1"] = len(ops_lower)
    print("  ... with degree <= %d: %d" % (degree - 1, len(ops_lower)))
    if len(ops) == 1:
        op = [PX.trim(c) for c in ops[0]]
        out["derived_coeff_degrees"] = [len(c) - 1 for c in op]
        if cop is not None:
            mine = [sp.Poly(sum(sp.Integer(x) * Rsym ** j for j, x in enumerate(c)), Rsym) for c in op]
            # find scalar: ratio of leading coefficients of the top-order coefficient
            ratio = sp.Rational(cop[-1].LC()) / sp.Rational(mine[-1].LC())
            eq = all(sp.expand(cop[i].as_expr() - ratio * mine[i].as_expr()) == 0 for i in range(len(cop)))
            out["equals_committed_operator_over_Q_up_to_scalar"] = bool(eq)
            out["scalar"] = str(ratio)
            print("  derived operator == committed operator_Q_crt (up to scalar %s): %s" % (ratio, eq))
if DO_NUM:
    import mpmath as mp
    mp.mp.dps = 50
    alg = PX.Algebra(a, b) if a < b else None
    vals = {}
    for Rv_ in ("1.4", "3.015"):
        x = mp.mpf(Rv_)
        f0 = PX.f0_numeric(alg, x)
        vals[Rv_] = mp.nstr(PX.evaluate(alg, K, x, f0), 48)
        print("  closed form K(%s) =" % Rv_, vals[Rv_], flush=True)
    out["closed_form_values_dps50"] = vals
out["time_s"] = round(time.time() - t0, 1)
json.dump(out, open(f"debug/data/pachucki_ext_lgt0_{tg}.json", "w"), indent=1)
print("done %.0fs" % (time.time() - t0))
