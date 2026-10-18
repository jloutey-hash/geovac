"""Loop test PT2, step D8: confirm D7's M5 > f3 identification with longer series, and test M5's
top (left) factor against the elliptic family at argument t^2/4 (singular at t = +-2)."""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_D_gauge as D
import density_loop_D4_factors as F4
F4.NT = 900
import density_loop_D5_rational as F5
F5.NT = 900
u = sp.Symbol("u")
rec = {}
M5 = [sp.sympify(q) for q in json.loads((Path(__file__).parent / "data" / "density_loop_D7.json").read_text())["M5"]]
B = json.loads((Path(__file__).parent / "data" / "density_loop_B.json").read_text())
QG = [sp.sympify(q) for q in B["QG"]]
NT = 900
f3 = D.hyper_t([Fraction(1, 2)] * 3, [1, 1], NT=NT, scale=Fraction(1, 4))
f2 = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT)
f2q = D.hyper_t([Fraction(1, 2)] * 2, [1], NT=NT, scale=Fraction(1, 4))
jobs = [("M5 > f3 (R), NT 900", M5, f3, 3, "R"),
        ("triangle > f2 (L), NT 900", QG, f2, 2, "L"),
        ("M5 top vs f2(t^2/4) (L)", M5, f2q, 2, "L"),
        ("M5 bottom vs f2(t^2/4) (R)", M5, f2q, 2, "R"),
        ("triangle vs f2(t^2/4) (L)", QG, f2q, 2, "L"),
        ("triangle vs f2(t^2/4) (R)", QG, f2q, 2, "R")]
for name, Q, f, r, side in jobs:
    QQ = Q if side == "R" else F4.adjoint(Q)
    hit = 0
    for b in range(4):
        for c in range(4):
            hit = F5.test(QQ, f, r, 12, b, c, a=6)
            if hit:
                break
        if hit:
            break
    rec[name] = hit
    print(f"{name:32s}: {'FOUND (nullity %d)' % hit if hit else 'none'}", flush=True)
(Path(__file__).parent / "data" / "density_loop_D8.json").write_text(json.dumps(rec, indent=1))
