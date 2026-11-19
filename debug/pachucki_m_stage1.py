"""Stage 1 (over Q, one pair): transcription check of the w1-route against label 61/63, f_2 jets, control (a).

Control (a) (no r12 weight):  nabla_1^2 r12 = 2/r12, so with psi = e^{-a r_A - b r_B}
    int int psi(1) psi(2) / r12 = (1/2) int int [nabla^2 psi](1) psi(2) r12,
    nabla^2 psi = [a^2 + b^2 - 2a/r_A - 2b/r_B + ab (r_A/r_B + r_B/r_A - R^2/(r_A r_B))] psi
=>  16 f_{1111} = 8 [ (a^2+b^2) f2_{1111} + 2a f2_{0111} + 2b f2_{1011} + ab (f2_{2011} + f2_{0211}) - ab R^2 f2_{0011} ]
identically (f = f_0 jets from the existing engine, f2 = jets of f_2 from the new route).
Run: python debug/pachucki_m_stage1.py A B
"""
import json
import sys
import time
from fractions import Fraction

sys.path.insert(0, "debug")
import pachucki_ext_core as E  # noqa: E402
import pachucki_m_core as M  # noqa: E402

t0 = time.time()
if sys.argv[1] == "sym":
    ctx = E.Ctx("Qt")
    tg = "sym"
else:
    ctx = E.Ctx("Q", Fraction(sys.argv[1]), Fraction(sys.argv[2]))
    tg = f"{sys.argv[1]}_{sys.argv[2]}".replace("/", "-")
a, b = ctx.a, ctx.b
C_ = ctx.c
D = E.downset([(2, 2, 1, 1)])
calc0, fS, Fs, P2, nfree = E.solve_all(ctx, D, check=True)
calc = M.MCalc(ctx)
RJ = M.RJets(ctx, D, calc)
Z = (0, 0, 0, 0)
out = {"pair": [str(a), str(b)], "field": ctx.kind}

# --- 0. the log decomposition used by the engine against label 13 verbatim:
#        P_full - sum C ln(t + c + eps w1) must be a polynomial in t of degree <= 2 (third t-derivative zero)
import sympy as sp  # noqa: E402
rem = M.P_full() - sum(C * sp.log(ts_ + c + e * M.w1s) for (C, c, e), ts_ in zip(M.P_logs(), [M.ts] * 4))
d3 = sp.diff(rem, M.ts, 3)
pt = {M.u2s: sp.Rational(17, 10), M.u3s: sp.Rational(9, 10), M.w2s: sp.Rational(11, 10), M.w3s: sp.Rational(23, 10),
      M.w1s: sp.Rational(3, 10), M.ts: sp.Rational(7, 5)}
val = abs(sp.N(d3.subs(pt), 40))
ctl = abs(sp.N(sp.diff(M.P_full(), M.ts, 3).subs(pt), 40))
chk0 = bool(val < 1e-30 and ctl > 1e-6)
print("label 13 verbatim minus the engine's log decomposition: third t-derivative %.1e (P itself: %.3e):" % (val, ctl), chk0,
      flush=True)
out["P_log_decomposition_matches_label13"] = chk0

# --- 1. S_1 at the exchange point against label 63 (W_1^{(-2)} differentiated twice)
S1 = RJ.source(1)
u, w = a, b
Wm2 = {("E", -2 * u, 0): {-1: -(-1 / (w - u) - 1 / (u + w))},
       ("E", -2 * w, 0): {-1: -(-1 / (u - w) - 1 / (u + w))},
       ("E", -(u + w), 0): {0: -C_(2), -1: -2 / (u + w)},
       ("Ei", ctx.zero, -(u + w)): {0: C_(2)}}
chk1 = M.fadd(calc.D2(Wm2), S1[Z]) == {}     # S_1 = -W_1
print("S_1(exchange point) == -(label 63 W_1^{(-2)})'':", chk1, flush=True)
out["S1_equals_label63_second_derivative"] = chk1

# --- 2. f_2 jets
f2 = M.solve_fn(calc, RJ, D, 1, fS, check=True)
print("f_2 jets solved on %d multi-indices (residuals identically 0)  %.1fs" % (len(D), time.time() - t0), flush=True)
# label 61 at n = 1:  (u^2-w^2)^2 f_2 = 2 (u^2+w^2) f_0 - f_0'' + W_1^{(-2)}
lab61 = M.fscale(M.fadd(M.fadd(M.fscale(fS[Z], 2 * (u * u + w * w)), calc.D2(fS[Z]), -1), Wm2), 1 / (u * u - w * w) ** 2)
chk2 = M.fadd(lab61, f2[Z], -1) == {}
print("f_2(exchange point) == label 61 (n = 1):", chk2, flush=True)
out["f2_equals_label61"] = chk2

# --- 3. control (a)
lhs = M.fscale(fS[(1, 1, 1, 1)], C_(16))
rhs = {}
for al, co, sh in (((1, 1, 1, 1), a * a + b * b, 0), ((0, 1, 1, 1), 2 * a, 0), ((1, 0, 1, 1), 2 * b, 0),
                   ((2, 0, 1, 1), a * b, 0), ((0, 2, 1, 1), a * b, 0), ((0, 0, 1, 1), -a * b, 2)):
    rhs = M.fadd(rhs, M.fscale(f2[al], 8 * co, sh))
chk3 = M.fadd(lhs, rhs, -1) == {}
print("control (a): 1s exchange closed form from f_2 jets (Laplacian identity) == existing, exactly:", chk3, flush=True)
out["control_a_exact"] = chk3
# --- 4. the same control one rung up (m = 3 machinery):  nabla^2 r12^3 = 12 r12
#        (a^2+b^2) f4_{1111} + 2a f4_{0111} + 2b f4_{1011} + ab (f4_{2011} + f4_{0211}) - ab R^2 f4_{0011} = 12 f2_{1111}
f4 = M.solve_fn(calc, RJ, D, 3, f2, fS, check=True)
print("f_4 jets solved on %d multi-indices (residuals identically 0)  %.1fs" % (len(D), time.time() - t0), flush=True)
# label 61 at n = 3:  (u^2-w^2)^2 f_4 = 18 (u^2+w^2) f_2 - 9 f_2'' - 12 f_0 + W_3^{(-2)};  checked twice differentiated,
# with W_3^{(-2)}'' = -S_3 (so no hand transcription of d^3/dw1^3 of label 63 is needed for the f-dependent part)
S3 = RJ.source(3)
lhs61 = M.fscale(calc.D2(f4[Z]), (u * u - w * w) ** 2)
rhs61 = M.fadd(M.fadd(M.fadd(M.fscale(calc.D2(f2[Z]), 18 * (u * u + w * w)), M.fscale(calc.D2(calc.D2(f2[Z])), C_(9)), -1),
                      M.fscale(calc.D2(fS[Z]), C_(12)), -1), S3[Z], -1)
chk4b = M.fadd(lhs61, rhs61, -1) == {}
print("f_4(exchange point)'' consistent with label 61 (n = 3) coefficients 18, -9, -12:", chk4b, flush=True)
out["f4_consistent_with_label61_n3"] = chk4b
lhs = M.fscale(f2[(1, 1, 1, 1)], C_(12))
rhs = {}
for al, co, sh in (((1, 1, 1, 1), a * a + b * b, 0), ((0, 1, 1, 1), 2 * a, 0), ((1, 0, 1, 1), 2 * b, 0),
                   ((2, 0, 1, 1), a * b, 0), ((0, 2, 1, 1), a * b, 0), ((0, 0, 1, 1), -a * b, 2)):
    rhs = M.fadd(rhs, M.fscale(f4[al], co, sh))
chk4 = M.fadd(lhs, rhs, -1) == {}
print("control (a'), one rung up: f_4 jets reproduce 12 f2_{1111} (Laplacian identity), exactly:", chk4, flush=True)
out["control_a_f4_exact"] = chk4
out["f2_1111_in_nine_span"] = M.in_nine(ctx, f2[(1, 1, 1, 1)]) is not None
out["f4_1111_in_nine_span"] = M.in_nine(ctx, f4[(1, 1, 1, 1)]) is not None
print("f2_{1111}, f4_{1111} in nine-function span:", out["f2_1111_in_nine_span"], out["f4_1111_in_nine_span"])
out["time_s"] = round(time.time() - t0, 1)
json.dump(out, open(f"debug/data/pachucki_m_stage1_{tg}.json", "w"), indent=1)
print("done %.1fs" % (time.time() - t0))
