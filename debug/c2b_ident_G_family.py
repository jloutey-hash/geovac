"""EC3 (docs/predictions/2026-10-22_c2b_identification.md): the elliptic family behind C2b, explicitly.
With P = 16x^2 + 152x + 25, R = (4x - 5)(16x^2 - 616x + 25), P^3 - R^2 = 1728 x (4x - 1)^2 (4x - 25)^2:
  E_e : Y^2 = X^3 - 3 e^2 P X - 2 e^3 R   has j = P^3 / (x (4x-1)^2 (4x-25)^2) = J.
Checks: the identity, the discriminant and fibre types for e = 1 and e = x, the pull-back to t (x = t^2),
and numerically that x^(1/2) P^(-1/4) 2F1(1/12, 5/12; 1; 1728/J) is a constant times a period of E_x."""
import json
from pathlib import Path
import mpmath as mp
import sympy as sp

x, t, X = sp.symbols("x t X")
P = 16 * x**2 + 152 * x + 25
R = (4 * x - 5) * (16 * x**2 - 616 * x + 25)
D = x * (4 * x - 1) ** 2 * (4 * x - 25) ** 2
out = {"identity_P3_minus_R2": str(sp.expand(P**3 - R**2 - 1728 * D) == 0)}
for name, e in (("e=1", sp.Integer(1)), ("e=x", x)):
    A, B = -3 * e**2 * P, -2 * e**3 * R
    disc = sp.factor(-16 * (4 * A**3 + 27 * B**2))
    j = sp.factor(sp.cancel(1728 * 4 * A**3 / (4 * A**3 + 27 * B**2)))
    # orders at infinity on a rational surface: weights 4, 6, 12
    oinf = (4 - sp.degree(A, x), 6 - sp.degree(B, x), 12 - sp.degree(sp.expand(4 * A**3 + 27 * B**2), x))
    out[name] = {"A": str(sp.expand(A)), "B": str(sp.expand(B)), "disc": str(disc), "j": str(j),
                 "ord_inf_A_B_Delta": [int(v) for v in oinf]}
    print(name, "disc =", disc, "| j == J:", sp.simplify(j - P**3 / D) == 0, "| ord_inf (A,B,Delta) =", oinf)
A1, B1 = (-3 * P).subs(x, t**2), (-2 * R).subs(x, t**2)
disc_t = sp.factor(-16 * (4 * A1**3 + 27 * B1**2))
out["over_t"] = {"A": str(sp.expand(A1)), "B": str(sp.expand(B1)), "disc": str(disc_t),
                 "ord_inf_Delta": int(12 - sp.degree(sp.expand(4 * A1**3 + 27 * B1**2), t))}
print("over t (e = 1, x = t^2): disc =", disc_t, "| ord_inf Delta =", out["over_t"]["ord_inf_Delta"])

# numerical period check on E_x near x = 0 (where |1728/J| < 1)
mp.mp.dps = 40
def period(xv):
    e = xv
    a = -3 * e**2 * (16 * xv**2 + 152 * xv + 25)
    b = -2 * e**3 * (4 * xv - 5) * (16 * xv**2 - 616 * xv + 25)
    r = sorted(mp.polyroots([1, 0, a, b]), key=lambda z: -mp.re(z))
    e1, e2, e3 = r
    return mp.pi / mp.agm(mp.sqrt(e1 - e3), mp.sqrt(e1 - e2))   # real period of dX/Y, up to a constant
def source(xv):
    Pv = 16 * xv**2 + 152 * xv + 25
    J = Pv**3 / (xv * (4 * xv - 1) ** 2 * (4 * xv - 25) ** 2)
    return mp.sqrt(xv) * Pv ** mp.mpf(-0.25) * mp.hyp2f1(mp.mpf(1) / 12, mp.mpf(5) / 12, 1, 1728 / J)
ratios = []
def periods(xv):
    a = -3 * xv**2 * (16 * xv**2 + 152 * xv + 25)
    b = -2 * xv**3 * (4 * xv - 5) * (16 * xv**2 - 616 * xv + 25)
    e1, e2, e3 = sorted([mp.re(z) for z in mp.polyroots([1, 0, a, b])], reverse=True)
    w1 = mp.pi / mp.agm(mp.sqrt(e1 - e3), mp.sqrt(e1 - e2))   # real period / 2 of dX/Y
    w2 = mp.pi / mp.agm(mp.sqrt(e1 - e3), mp.sqrt(e2 - e3))   # imaginary period / 2i
    return w1, w2
for xv in (mp.mpf("0.004"), mp.mpf("0.01"), mp.mpf("0.02"), mp.mpf("0.04")):
    w1, w2 = periods(xv)
    s0 = source(xv)
    ratios.append((mp.nstr(xv, 4), mp.nstr(xv * w1 / s0, 30), mp.nstr(xv * w2 / s0, 30)))
    print("x =", mp.nstr(xv, 4), " x w1/source =", mp.nstr(xv * w1 / s0, 30), "  x w2/source =", mp.nstr(xv * w2 / s0, 30))
out["x_period_over_source"] = ratios
json.dump(out, open(Path(__file__).parent / "data" / "c2b_ident_G_family.json", "w"), indent=1)
