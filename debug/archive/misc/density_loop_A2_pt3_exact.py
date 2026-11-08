"""Loop test PT3: exact one-centre triangle values from Remiddi's master (validated against the
star-triangle route to 1e-21 by density_loop_A_pt1_pt3.py), plus a 100-digit PSLQ against Q.

T(z) = 64 prod z^3 * (-d_w1)(-d_w2)(-d_w3) f000(w),  w = 2 z,
f000 = -(1/(2 w1 w2 w3)) sum_cyc [ln x ln(1+x) + Li2(-x) + Li2(1-x)],  x = w_i/(w_j+w_k).
"""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

OUT = Path(__file__).parent / "data" / "density_loop_A2.json"

w1, w2, w3 = sp.symbols("w1 w2 w3", positive=True)


def term(a, b, c):
    x = a / (b + c)
    return sp.log(x) * sp.log(1 + x) + sp.polylog(2, -x) + sp.polylog(2, 1 - x)


F = -(term(w3, w1, w2) + term(w2, w3, w1) + term(w1, w2, w3)) / (2 * w1 * w2 * w3)
D3 = -sp.diff(F, w1, w2, w3)          # (-d1)(-d2)(-d3) F


def exact_T(z):
    sub = {w1: 2 * z[0], w2: 2 * z[1], w3: 2 * z[2]}
    val = 64 * sp.Integer(z[0] * z[1] * z[2]) ** 3 * D3.subs(sub)
    return sp.nsimplify(sp.expand_log(sp.simplify(val), force=True))


def main() -> None:
    rec = {}
    mp.mp.dps = 110
    for z in [(1, 1, 1), (1, 2, 3)]:
        ex = exact_T(z)
        num = sp.N(ex, 105)
        tv = mp.mpf(str(num))
        rel = mp.pslq([tv, mp.mpf(1)], maxcoeff=10 ** 12, maxsteps=10 ** 7)
        atoms = sorted({str(a) for a in ex.atoms(sp.log, sp.polylog, sp.atan)})
        rec[str(z)] = {"exact": str(ex), "value_60": mp.nstr(tv, 60),
                       "pslq_T_1": str(rel), "transcendental_atoms": atoms}
        print(z, "exact:", ex, flush=True)
        print("   value:", mp.nstr(tv, 40), " pslq[T,1]:", rel, " atoms:", atoms, flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
