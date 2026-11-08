"""Round-1d: assemble the explicit symbolic closed form of S_min from the
round-1b identified t-value reductions, simplify, and verify at 200+ dps.

Reads debug/data/smin_dossier_round1b.json (must show fully_reduced or the
per-t-value coefficient dicts), builds the sympy expression, and prints both
sympy and LaTeX forms.  Appends the result to the JSON.
"""
import json
from fractions import Fraction
from pathlib import Path

import mpmath
import sympy as sp

mpmath.mp.dps = 230

data_path = Path(__file__).parent / "data" / "smin_dossier_round1b.json"
data = json.loads(data_path.read_text())

spi = sp.pi
sln2 = sp.log(2)
sz3, sz5, sz7 = sp.zeta(3), sp.zeta(5), sp.zeta(7)
sLi4 = sp.Symbol("Li4half")   # polylog(4, 1/2)
sLi5 = sp.Symbol("Li5half")
sLi6 = sp.Symbol("Li6half")
sLi7 = sp.Symbol("Li7half")
sA1 = sp.Symbol("zeta_m6_1")  # zeta(-6,1)
sA2 = sp.Symbol("zeta_m5_2")  # zeta(-5,2)

SYM = {
    "zeta3": sz3, "zeta5": sz5, "zeta7": sz7,
    "pi2*ln2": spi ** 2 * sln2, "pi4*ln2": spi ** 4 * sln2,
    "pi6*ln2": spi ** 6 * sln2,
    "pi2*zeta3": spi ** 2 * sz3, "pi2*zeta5": spi ** 2 * sz5,
    "pi4*zeta3": spi ** 4 * sz3,
    "ln2^3": sln2 ** 3, "ln2^5": sln2 ** 5, "ln2^7": sln2 ** 7,
    "zeta3*ln2^2": sz3 * sln2 ** 2, "pi2*ln2^3": spi ** 2 * sln2 ** 3,
    "zeta5*ln2^2": sz5 * sln2 ** 2, "zeta3^2*ln2": sz3 ** 2 * sln2,
    "zeta3*pi2*ln2^2": sz3 * spi ** 2 * sln2 ** 2,
    "zeta3*ln2^4": sz3 * sln2 ** 4, "pi4*ln2^3": spi ** 4 * sln2 ** 3,
    "pi2*ln2^5": spi ** 2 * sln2 ** 5,
    "Li4*ln2": sLi4 * sln2, "Li5": sLi5,
    "Li4*zeta3": sLi4 * sz3, "Li4*pi2*ln2": sLi4 * spi ** 2 * sln2,
    "Li4*ln2^3": sLi4 * sln2 ** 3, "Li5*pi2": sLi5 * spi ** 2,
    "Li5*ln2^2": sLi5 * sln2 ** 2, "Li6*ln2": sLi6 * sln2, "Li7": sLi7,
    "zeta(-6,1)": sA1, "zeta(-5,2)": sA2,
}

NUM = {
    sLi4: mpmath.polylog(4, mpmath.mpf(1) / 2),
    sLi5: mpmath.polylog(5, mpmath.mpf(1) / 2),
    sLi6: mpmath.polylog(6, mpmath.mpf(1) / 2),
    sLi7: mpmath.polylog(7, mpmath.mpf(1) / 2),
}


def slam(s):
    zeven = {2: spi ** 2 / 6, 4: spi ** 4 / 90, 6: spi ** 6 / 945,
             8: spi ** 8 / 9450}
    z = zeven[s] if s % 2 == 0 else {3: sz3, 5: sz5, 7: sz7}[s]
    return (1 - sp.Rational(1, 2 ** s)) * z


def build(coeff_dict):
    return sum(sp.Rational(Fraction(fr)) * SYM[nm]
               for nm, fr in coeff_dict.items())


t21 = build(data["t21"])
t23 = build(data["t23"])
t41 = build(data["t41"])
t43 = build(data.get("t43_products_only") or data["t43_with_depth2"])
tsym = {(2, 1): t21, (2, 3): t23, (4, 1): t41, (4, 3): t43,
        (2, 2): (slam(2) ** 2 - slam(4)) / 2,
        (4, 4): (slam(4) ** 2 - slam(8)) / 2}
pair24 = slam(2) * slam(4) - slam(6)

c_coeffs = {1: 1, 2: -3, 3: -1, 4: 3}
b_coeffs = {2: 1, 4: -1}
d_coeffs = {3: 1, 4: -3, 5: -2, 6: 6, 7: 1, 8: -3}

S = sp.Integer(0)
for c, ec in c_coeffs.items():
    for b, eb in b_coeffs.items():
        bnd = ((slam(b) - 1)
               + sp.Rational(1, 3 ** c) * (slam(b) - 1 - sp.Rational(1, 3 ** b)))
        if (b, c) in ((2, 4), (4, 2)):
            S += 64 * ec * eb * (-bnd)
        else:
            S += 64 * ec * eb * (tsym[(b, c)] - bnd)
S += 64 * 3 * pair24
S += 32 * sum(ds * (slam(s) - 1 - sp.Rational(1, 3 ** s))
              for s, ds in d_coeffs.items())

S = sp.nsimplify(sp.expand(S))
S = sp.collect(sp.expand(S), [sz3, sz5, sz7, sln2, sLi4, sLi5, sLi6, sLi7,
                              sA1, sA2])
print("S_min closed form (collected):\n")
sp.pprint(S, use_unicode=False, wrap_line=False)
print("\nLaTeX:\n", sp.latex(S))

# numeric verification
expr_num = S
for sym, val in NUM.items():
    expr_num = expr_num.subs(sym, sp.Float(mpmath.nstr(val, 225), 225))
val = sp.N(expr_num, 220)
S_num = mpmath.mpf(str(val))
S_anchor = mpmath.nsum(
    lambda k: (2 * mpmath.hurwitz(2, mpmath.mpf(int(k)) + mpmath.mpf(3) / 2)
               - mpmath.mpf(1) / 2
               * mpmath.hurwitz(4, mpmath.mpf(int(k)) + mpmath.mpf(3) / 2)) ** 2,
    [1, mpmath.inf], method='levin')
res = abs(S_num - S_anchor)
print("\nresidual vs anchor at 200+ dps:", mpmath.nstr(res, 5))

data["symbolic_closed_form"] = str(S)
data["symbolic_closed_form_latex"] = sp.latex(S)
data["symbolic_residual"] = mpmath.nstr(res, 5)
data_path.write_text(json.dumps(data, indent=2))
print("Appended to", data_path)
