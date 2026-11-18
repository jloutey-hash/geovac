"""PART A: the 1s exchange integral over Q(t), a = 1, b = t (scaling: K(R; a, b) = a K(aR; 1, b/a)).

1. solve the 15 derivative equations (*) over Q(t) exactly (jets + direct solver, debug/pachucki_ext_core.py);
   every residual is checked to be identically 0 in Q(t)[R, 1/R] x keys;
2. specialise to t = 2, 3/2 and compare with the reviewed helper tests/_pachucki_exchange.py (Fraction route);
3. ties: K's coordinates satisfy the two psi-ties identically in t; rewrite K in the nine functions;
4. the nine-function span is closed under d/dR (exactly, over Q(t)); rank of K..K^(9) = 9 at a specialisation;
5. restore (a, b): coefficient of R^p g(R; a, b) = a^{p+1} c_p(b/a) (f0: a^{p+3} c_p(b/a)).
Run: python debug/pachucki_ext_symbolic.py"""
import json, sys, time
from fractions import Fraction
import sympy as sp
sys.path.insert(0, "debug"); sys.path.insert(0, "tests")
import pachucki_ext_core as E

t0 = time.time()
ctx = E.Ctx("Qt")
K_ = ctx.K
tS = ctx.t
D = E.downset([(1, 1, 1, 1)])
calc, fS, Fs, P2, nfree = E.solve_all(ctx, D, check=True)
print("15 equations solved over Q(t), residuals identically 0; free ansatz unknowns:", nfree, "%.1fs" % (time.time() - t0))
K = E.fscale(fS[(1, 1, 1, 1)], 16 * ctx.b ** 3, -1)
res = {"free_unknowns": nfree, "residuals_identically_zero": True}

# 2. specialisation vs the reviewed helper
import _pachucki_exchange as PX


def spec(A, tv):
    out = {}
    for k, d in A.items():
        kk = tuple(Fraction(str(K_.to_sympy(x).subs(tS, tv))) if (not isinstance(x, (int, str))) else x for x in k)
        for p, c in d.items():
            v = Fraction(str(K_.to_sympy(c).subs(tS, tv)))
            if v:
                out.setdefault(kk, {})[p] = out.get(kk, {}).get(p, 0) + v
    return E.clean(out)


def norm_keys(A):
    return {tuple(Fraction(x) if isinstance(x, (Fraction, int)) and not isinstance(x, bool) and i > 0 else x
                  for i, x in enumerate(k)): d for k, d in A.items()}


for tv in (sp.Integer(2), sp.Rational(3, 2)):
    alg = PX.Algebra(1, Fraction(str(tv)))
    fS2, *_ = alg.solve_all()
    K2 = alg.K(fS2)
    ok = norm_keys(spec(K, tv)) == norm_keys(K2)
    res["equals_helper_at_t=%s" % tv] = ok
    print("specialisation t =", tv, "equals tests/_pachucki_exchange.py K:", ok)

# 3. ties and the nine functions
a, b = ctx.a, ctx.b
def co(k, p=None):
    return K.get(k, {})
kA2, kA1, kAL = ("Ei", -2 * a, 2 * (a - b)), ("Ei", -2 * a, (a - b)), ("E", -2 * a, 1)
kB2, kB1, kBL = ("Ei", -2 * b, 2 * (b - a)), ("Ei", -2 * b, (b - a)), ("E", -2 * b, 1)
tie = all(E.fadd({"x": co(k1)}, {"x": co(k2)}, 2) == {} for k1, k2 in ((kA1, kA2), (kB1, kB2)))
tieL = all(E.fadd({"x": co(k1)}, {"x": co(k2)}) == {} for k1, k2 in ((kAL, kA2), (kBL, kB2)))
res["tie_Ei1_equals_minus2_Ei2_identically"] = tie
res["tie_LC_equals_minus_Ei2_identically"] = tieL
print("ties hold identically in t:", tie, tieL)
keys = sorted(K, key=str)
print("keys of K:", len(keys))

def psi(x, y):
    return {("Ei", -2 * x, 2 * (x - y)): {0: ctx.one}, ("Ei", -2 * x, (x - y)): {0: -2 * ctx.one},
            ("E", -2 * x, 1): {0: -ctx.one}}

nine = [{("f",): {0: ctx.one}}, {("Ei", ctx.zero, -(a + b)): {0: ctx.one}}, {("E", -(a + b), 0): {0: ctx.one}},
        {("E", -2 * a, 0): {0: ctx.one}}, {("E", -2 * b, 0): {0: ctx.one}},
        {("Ei", 2 * a, -2 * (a + b)): {0: ctx.one}}, {("Ei", 2 * b, -2 * (a + b)): {0: ctx.one}}, psi(a, b), psi(b, a)]
lead = [("f",), ("Ei", ctx.zero, -(a + b)), ("E", -(a + b), 0), ("E", -2 * a, 0), ("E", -2 * b, 0),
        ("Ei", 2 * a, -2 * (a + b)), ("Ei", 2 * b, -2 * (a + b)), kA2, kB2]
names = ["f0", "-E1((a+b)R)", "e^{-(a+b)R}", "e^{-2aR}", "e^{-2bR}", "-e^{2aR}E1(2(a+b)R)", "-e^{2bR}E1(2(a+b)R)",
         "psi_a", "psi_b"]


def in_nine(A):
    """coefficients (Laurent dicts) on the nine functions, or None if A is not in their Q(t)[R,1/R]-span."""
    rest = {k: dict(v) for k, v in A.items()}
    coeffs = []
    for g, l in zip(nine, lead):
        c = rest.get(l, {})
        coeffs.append(dict(c))
        for p, x in c.items():
            rest = E.fadd(rest, E.fscale(g, x, p), -1)
    return coeffs if rest == {} else None


cK = in_nine(K)
res["K_in_nine_function_span_over_Qt"] = cK is not None
print("K in the Q(t)[R,1/R]-span of the nine functions:", cK is not None)
closed = all(in_nine(calc.diff(g)) is not None for g in nine[1:]) and in_nine(calc.diff(nine[0])) is not None
res["nine_span_closed_under_d/dR"] = closed
print("nine-function span closed under d/dR:", closed)

# 4. rank of K, K', ..., K^(9) in the nine coordinates: rank over Q(t)(R) >= rank at a specialisation; <= 9
ders = [K]
for _ in range(9):
    ders.append(calc.diff(ders[-1]))
coords = [in_nine(A) for A in ders]
assert all(c is not None for c in coords)
Rv, tv = sp.Rational(7, 5), sp.Integer(3)
Mx = sp.Matrix([[sum(K_.to_sympy(x).subs(tS, tv) * Rv ** p for p, x in c.items()) for c in row] for row in coords])
rk = Mx.rank()
res["rank_K..K9_at_t=3_R=7/5"] = int(rk)
print("rank of K..K^(9) coordinates at t = 3, R = 7/5:", rk, "(hence 9 over Q(t)(R))")

# 5. closed form in (a, b): coefficient of R^p g_i(R; a, b)
A_, B_, R_ = sp.symbols("a b R", positive=True)
cf = {}
for nm, c in zip(names, cK):
    e = 0
    for p, x in c.items():
        sh = 3 if nm == "f0" else 1
        e += A_ ** (p + sh) * K_.to_sympy(x).subs(tS, B_ / A_) * R_ ** p
    cf[nm] = sp.factor(sp.simplify(e))
    print(" ", nm, ":", cf[nm])
res["closed_form_coefficients_ab"] = {k: str(v) for k, v in cf.items()}
res["closed_form_coefficients_t"] = {nm: {p: str(sp.factor(K_.to_sympy(x))) for p, x in c.items()}
                                     for nm, c in zip(names, cK)}
dens = set()
for nm, c in zip(names, cK):
    for x in c.values():
        for f, m in sp.factor_list(sp.denom(sp.factor(K_.to_sympy(x))))[1]:
            dens.add(str(f))
res["denominator_factors_in_t"] = sorted(dens)
print("denominator factors in t:", dens)
res["definitions"] = {"f0": "f(R) of arXiv:0905.1635 label 01 at u3=w2=a, u2=w3=b; f0(R;a,b) = a^-2 f0(aR;1,b/a)",
                      "psi_a": "e^{-2aR}[Ei(2(a-b)R) - 2 Ei((a-b)R) - LC], LC = ln((a+b)/(b-a))",
                      "psi_b": "e^{-2bR}[Ei(2(b-a)R) - 2 Ei((b-a)R) - LC]",
                      "K": "sum over the nine of coefficient * function; K(R;a,b) = a K(aR;1,b/a)"}
res["time_s"] = round(time.time() - t0, 1)
json.dump(res, open("debug/data/pachucki_ext_symbolic_1s.json", "w"), indent=1)
