"""CF record (d68e6fe5): check the recalled closed form (GR 6.672.6 at nu = 0) at real c, and its continuation
c = i g for the cosh part:  int_0^oo K0(a x) K0(b x) cos(c x) dx = (pi^2/(4 sqrt(ab))) P_(-1/2)(z),
z = (a^2 + b^2 + c^2)/(2ab).  For the cosh part z = (t^2 + 1 - g^2)/(2t) can lie in (-1, 1): Ferrers
(type 2) is tried there, alongside type 3.
Usage: python debug/c2b_feyn_D_formula.py"""
import json
from pathlib import Path

import mpmath as mp

mp.mp.dps = 15
DATA = Path(__file__).parent / "data"
out = {"real_c": [], "cosh": []}
for a, b, c in [(1, 1, 0.5), (0.7, 1.3, 2.0)]:
    a, b, c = mp.mpf(a), mp.mpf(b), mp.mpf(c)
    lhs = mp.quadosc(lambda y: mp.besselk(0, a * y) * mp.besselk(0, b * y) * mp.cos(c * y), [0, mp.inf], omega=c)
    z = (a**2 + b**2 + c**2) / (2 * a * b)
    rhs = mp.pi**2 / (4 * mp.sqrt(a * b)) * mp.legenp(-0.5, 0, z, type=3)
    out["real_c"].append(float(abs(lhs - rhs) / abs(lhs)))
    print("real c", (float(a), float(b), float(c)), "rel diff", out["real_c"][-1])
g = mp.mpf(3) / 2
for tv in [mp.mpf("0.8"), mp.mpf("1.7"), mp.mpf("3.2")]:
    lhs = mp.quad(lambda y: mp.besselk(0, tv * y) * mp.besselk(0, y) * mp.cosh(g * y), [0, 1, 4, 16, 64, mp.inf])
    z = (tv**2 + 1 - g**2) / (2 * tv)
    row = {"z": float(z)}
    for ty in (2, 3):
        try:
            rhs = mp.pi**2 / (4 * mp.sqrt(tv)) * mp.legenp(-0.5, 0, z, type=ty)
            row[f"type{ty}"] = float(abs(lhs - mp.re(rhs)) / abs(lhs))
        except Exception as e:  # noqa: BLE001
            row[f"type{ty}"] = str(e)
    out["cosh"].append(row)
    print("cosh t =", float(tv), row)
json.dump(out, open(DATA / "c2b_feyn_D_formula.json", "w"), indent=1)
