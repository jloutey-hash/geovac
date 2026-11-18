"""PART A, a -> b: limit of the derived closed form at a = 1, b = t = 1 + eps.

K(R; t) = sum_i c_i(t, R) g_i(R; t) (debug/data/pachucki_ext_symbolic_1s.json).  Every g_i is expanded in eps:
the elementary / E1 pieces by Taylor's theorem, the two psi brackets by Ei(y) = gamma + ln|y| + sum y^k/(k k!)
(the ln eps cancels inside each bracket), and f0(R; t) = sum_k eps^k I_k(R), I_k = int_R^oo FM1_k(s)/s ds with FM1_k
the eps^k coefficient of label 60 (f0 = int_R^oo F^(-1)(s)/s ds, label 59).  Each I_k is reduced by parts to the class
   {R^j e^{-2R}, R^j e^{-2R} Lg, E1(2R), R^j e^{2R}E1(4R)}   (Lg = gamma + ln R)
plus three integrals outside it, J1 = int_R^oo e^{-2s} Lg(s)/s ds, J2 = int_R^oo e^{2s}E1(4s)/s ds,
J3 = int_R^oo E1(2s)/s ds.  Output: the eps-poles (must vanish), the J-coefficients of the limit, and the limit
compared with Sugiura's closed form of the equal-exponent exchange integral.
Run: python debug/pachucki_ext_limit.py
"""
import json
import time
from functools import lru_cache

import sympy as sp

t0 = time.time()
R, eps = sp.symbols("R epsilon", positive=True)
e2, Lg, X2, X4, J1, J2, J3 = sp.symbols("e2 Lg X2 X4 J1 J2 J3")      # e^{-2R}, gamma+lnR, E1(2R), e^{2R}E1(4R)
import sys
SRC = sys.argv[1] if len(sys.argv) > 1 else "debug/data/pachucki_ext_symbolic_1s.json"
CASE = sys.argv[2] if len(sys.argv) > 2 else "1s"
COMMITTED11 = sys.argv[3] if len(sys.argv) > 3 else "debug/data/ep_heldout_1s_1_1_EPT.json"
d = json.load(open(SRC))
T = sp.Symbol("t")
cf_t = {nm: sum(sp.sympify(v) * R ** int(p) for p, v in c.items()) for nm, c in d["closed_form_coefficients_t"].items()}
MAXPOLE = max(dict(sp.factor_list(sp.denom(sp.factor(sp.sympify(v))))[1]).get(T - 1, 0)
              for c in d["closed_form_coefficients_t"].values() for v in c.values())
NORD = MAXPOLE + 2                                                     # eps orders kept
print("max pole order in t - 1:", MAXPOLE)


def ser(e, n=NORD):
    return sp.series(e, eps, 0, n).removeO()


x = sp.Symbol("x", positive=True)
E1d = [None]                       # E1^{(k)}(x), k >= 1
g = -sp.exp(-x) / x
for k in range(1, NORD + 2):
    E1d.append(g)
    g = sp.diff(g, x)


def E1_shift(x0, h):
    """E1(x0 + h) - E1(x0), Taylor polynomial in h (h = O(eps))."""
    return sum(h ** k / sp.factorial(k) * E1d[k].subs(x, x0) for k in range(1, NORD + 1))


def Ssum(y):
    return sum(y ** k / (k * sp.factorial(k)) for k in range(1, NORD + 1))


ER = sp.exp(-2 * R)
fun = {
    "-E1((a+b)R)": -(X2 + E1_shift(2 * R, eps * R)),
    "e^{-(a+b)R}": ER * sp.exp(-eps * R),
    "e^{-2aR}": ER,
    "e^{-2bR}": ER * sp.exp(-2 * eps * R),
    "-e^{2aR}E1(2(a+b)R)": -(X4 + sp.exp(2 * R) * E1_shift(4 * R, 2 * eps * R)),
    "-e^{2bR}E1(2(a+b)R)": -sp.exp(2 * eps * R) * (X4 + sp.exp(2 * R) * E1_shift(4 * R, 2 * eps * R)),
    # Ei(2(a-b)R) - 2 Ei((a-b)R) - LC with Ei(y) = gamma + ln|y| + S(y), LC = ln((2+eps)/eps): ln eps cancels
    "psi_a": ER * (-Lg - sp.log(1 + eps / 2) + Ssum(-2 * eps * R) - 2 * Ssum(-eps * R)),
    "psi_b": ER * sp.exp(-2 * eps * R) * (-Lg - sp.log(1 + eps / 2) + Ssum(2 * eps * R) - 2 * Ssum(eps * R)),
}
fun = {k: ser(sp.expand(v)) for k, v in fun.items()}
c = 1 / (2 * (1 - (1 + eps) ** 2))           # label 60 prefactor
FM1 = c * (fun["psi_a"] - fun["psi_b"] + fun["-e^{2bR}E1(2(a+b)R)"] - fun["-e^{2aR}E1(2(a+b)R)"])
FM1 = ser(sp.expand(FM1), NORD - 1)
print("FM1 expanded %.0fs" % (time.time() - t0), flush=True)


def to_class(e):
    e = sp.expand(e)
    e = sp.expand(e.subs(sp.exp(-2 * R), e2).subs(sp.exp(2 * R), 1 / e2))
    assert not e.has(sp.exp), e.atoms(sp.exp)
    return e


@lru_cache(maxsize=None)
def I_e(j):
    """int_R^oo s^j e^{-2s} ds."""
    if j >= 0:
        return sp.expand(sum(sp.factorial(j) / sp.factorial(i) * R ** i / sp.Integer(2) ** (j - i + 1)
                             for i in range(j + 1)) * e2)
    if j == -1:
        return X2
    return sp.expand(-R ** (j + 1) * e2 / (j + 1) + sp.Rational(2, j + 1) * I_e(j + 1))


@lru_cache(maxsize=None)
def I_eL(j):
    """int_R^oo s^j e^{-2s} Lg(s) ds."""
    if j >= 0:
        P = sp.expand(I_e(j) / e2)          # d/ds[-P e^{-2s}] = s^j e^{-2s}
        Ps = sp.Poly(P, R)
        rest = sum(cf * I_e(m[0] - 1) for m, cf in zip(Ps.monoms(), Ps.coeffs()))
        return sp.expand(P * e2 * Lg + rest)
    if j == -1:
        return J1
    return sp.expand(-R ** (j + 1) * e2 * Lg / (j + 1) + sp.Rational(2, j + 1) * I_eL(j + 1)
                     - sp.Rational(1, j + 1) * I_e(j))


@lru_cache(maxsize=None)
def I_X4(j):
    """int_R^oo s^j X4(s) ds, X4' = 2 X4 - e^{-2s}/s."""
    if j >= 0:
        return sp.expand(-R ** j * X4 / 2 - (sp.Rational(j, 2) * I_X4(j - 1) if j > 0 else 0)
                         + sp.Rational(1, 2) * I_e(j - 1))
    if j == -1:
        return J2
    return sp.expand(-R ** (j + 1) * X4 / (j + 1) - sp.Rational(2, j + 1) * I_X4(j + 1) + sp.Rational(1, j + 1) * I_e(j))


@lru_cache(maxsize=None)
def I_X2(j):
    """int_R^oo s^j X2(s) ds, X2' = -e^{-2s}/s."""
    if j == -1:
        return J3
    return sp.expand(-R ** (j + 1) * X2 / (j + 1) + I_e(j) / (j + 1))


def integrate_over_s(G):
    """int_R^oo G(s)/s ds, G a combination of R^j {e2, e2 Lg, X4, X2}."""
    out = 0
    for term in sp.Add.make_args(sp.expand(G / R)):
        cR, rest = term.as_independent(e2, Lg, X4, X2)
        cc, j = cR.as_coeff_exponent(R)
        assert not cc.has(R), term
        j = int(j)
        table = {e2: I_e, e2 * Lg: I_eL, X4: I_X4, X2: I_X2}
        if rest not in table:
            raise ValueError(("unexpected", term))
        out += cc * table[rest](j)
    return sp.expand(out)


rep = {e2: sp.exp(-2 * R), Lg: sp.EulerGamma + sp.log(R), X2: sp.expint(1, 2 * R),
       X4: sp.exp(2 * R) * sp.expint(1, 4 * R)}
for j in range(-4, 5):
    for fn, ig in ((I_e, e2), (I_eL, e2 * Lg), (I_X4, X4), (I_X2, X2)):
        v = fn(j)
        if v.has(J1, J2, J3):
            continue
        chk = sp.simplify(sp.diff(v.subs(rep), R) + (R ** j * ig).subs(rep))
        assert chk == 0, (fn.__name__, j, chk)
print("reducers verified: d/dR of each = -integrand, j = -4..4", flush=True)

FM1c = to_class(FM1)
Ik = [integrate_over_s(sp.expand(FM1c.coeff(eps, k))) for k in range(NORD - 1)]
f0exp = sum(eps ** k * Ik[k] for k in range(NORD - 1))
print("I_k reduced %.0fs; J content of I_k:" % (time.time() - t0),
      [[str(sp.factor(Ik[k].coeff(J))) for J in (J1, J2, J3)] for k in range(NORD - 1)], flush=True)
total = 0
for nm, cft in cf_t.items():
    cser = sp.series(cft.subs(T, 1 + eps), eps, 0, 1).removeO()
    gi = f0exp if nm == "f0" else to_class(fun[nm])
    total += sp.expand(cser * gi)
total = sp.expand(total)
poles = {m: sp.simplify(total.coeff(eps, -m)) for m in range(1, MAXPOLE + 2)}
lim = sp.expand(total.coeff(eps, 0))
print("poles vanish:", {m: (v == 0) for m, v in poles.items()})
Jc = {str(J): sp.simplify(lim.coeff(J)) for J in (J1, J2, J3)}
print("J coefficients in the limit:", Jc)
S = 1 + R + R ** 2 / 3
Sp = 1 - R + R ** 2 / 3
sug = sp.Rational(1, 5) * (-e2 * (-sp.Rational(25, 8) + sp.Rational(23, 4) * R + 3 * R ** 2 + R ** 3 / 3)
                           + 6 / R * (e2 * S ** 2 * Lg - Sp ** 2 * X4 + 2 * S * Sp * X2))   # Ei(-x) = -E1(x)
dif = sp.simplify(lim - sp.expand(sug)) if CASE == "1s" else "n/a"
print("limit - Sugiura:", dif)
# ---- minimal operator of the limit at a = b = 1, compared with the committed (record) operator over Q
basis = [e2 * Lg, e2, X2, X4]


def coords(e):
    e = sp.expand(e)
    cL = e.coeff(e2 * Lg)
    rest = sp.expand(e - cL * e2 * Lg)
    out_ = [cL, rest.coeff(e2), rest.coeff(X2), rest.coeff(X4)]
    assert sp.expand(rest - out_[1] * e2 - out_[2] * X2 - out_[3] * X4) == 0
    return out_


def dcoords(c4):
    cL, ce, c2, c4_ = c4
    return [sp.diff(cL, R) - 2 * cL, sp.diff(ce, R) - 2 * ce + cL / R - c2 / R - c4_ / R,
            sp.diff(c2, R), sp.diff(c4_, R) + 2 * c4_]


Jop = json.load(open(COMMITTED11))
ORD, DEG = Jop["order"], Jop["degree"]
ders = [coords(lim)]
for _ in range(ORD):
    ders.append([sp.expand(x) for x in dcoords(ders[-1])])
qs = sp.symbols("q0:%d" % ((ORD + 1) * (DEG + 1)))
eqs = []
for bi in range(4):
    expr = sp.expand(sum(qs[i * (DEG + 1) + j] * R ** j * ders[i][bi] for i in range(ORD + 1) for j in range(DEG + 1)))
    num = sp.numer(sp.together(expr))
    eqs += sp.Poly(sp.expand(num), R).coeffs()
Msys = sp.Matrix([[sp.expand(e).coeff(q) for q in qs] for e in eqs])
ns = Msys.nullspace()
print("operators of order <= %d, degree <= %d annihilating the limit: %d" % (ORD, DEG, len(ns)))
match = None
if len(ns) == 1:
    v = ns[0]
    Rsym = sp.Symbol("R")
    mine = [sp.expand(sum(v[i * (DEG + 1) + j] * R ** j for j in range(DEG + 1))) for i in range(ORD + 1)]
    theirs = [sp.sympify(s_, locals={"R": R}) for s_ in Jop["operator_Q_crt"]]
    ratio = sp.Poly(theirs[-1], R).LC() / sp.Poly(mine[-1], R).LC()
    match = all(sp.expand(theirs[i] - ratio * mine[i]) == 0 for i in range(ORD + 1))
    print("limit's operator == committed (1,1) operator over Q up to scalar:", match)
out = {"case": CASE, "source": SRC, "committed_11": COMMITTED11, "order_degree": [ORD, DEG],
       "n_operators": len(ns), "operator_equals_committed_11": match, "poles_vanish": {str(m): bool(v == 0) for m, v in poles.items()},
       "J_coefficients_of_limit": {k: str(v) for k, v in Jc.items()},
       "J_content_of_I_k": [[str(sp.factor(Ik[k].coeff(J))) for J in (J1, J2, J3)] for k in range(NORD - 1)],
       "limit": str(sp.collect(lim, [e2 * Lg, e2, X2, X4])), "limit_minus_sugiura": str(dif),
       "symbols": "e2 = e^{-2R}, Lg = gamma + ln R, X2 = E1(2R), X4 = e^{2R} E1(4R); a = b = 1",
       "time_s": round(time.time() - t0, 1)}
json.dump(out, open("debug/data/pachucki_ext_limit%s.json" % ("" if CASE == "1s" else "_" + CASE), "w"), indent=1)
print("done %.0fs" % (time.time() - t0))
