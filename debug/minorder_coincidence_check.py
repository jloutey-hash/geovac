"""The second hybrid rate coincidence a_d = a_p + a_c  ((2p0 2p0|1s 1s_B), Z_A = 3, Z_B = 6).

Checks: (1) atom content and constants; (2) ln R coefficient of f at R = 0 (exact); (3) the closed form against
the module's own quadrature and against the closed form at neighbouring Z_B (continuity in Z_B).
Writes debug/data/minorder_coincidence_check.json.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "debug"))

import minorder_core as MC  # noqa: E402
from geovac.two_center_eri import hybrid_closed_form, hybrid_quadrature  # noqa: E402

F = Fraction
RO, R = MC.RO, MC.R
S1, P0, P1 = (1, 0, 0), (2, 1, 0), (2, 1, 1)

out = {}
for tag, mk, zb0 in (("fam1_ZB6", lambda zb: hybrid_closed_form(F(3), P0, P0, S1, F(zb), S1), F(6)),
                     ("fam2_ZB6", lambda zb: hybrid_closed_form(F(2), P1, P1, P1, F(zb), P1), F(6))):
    f = mk(zb0)
    atoms, reps = RO.decompose(f)
    rec = {"atoms": [str(a) for a in atoms], "constants": [str(k) for k in reps]}
    rec["lnR_coefficient_of_f_at_0"] = str(RO.lnR_coefficient(f))
    vals = {}
    for Rv in (sp.Rational(7, 5), sp.Rational(1, 2)):
        v0 = sp.N(f.subs(R, Rv), 30)
        nb = {}
        for eps in (F(1, 1000), F(-1, 1000), F(1, 100000), F(-1, 100000)):
            nb[str(eps)] = str(sp.N(mk(zb0 + eps).subs(R, Rv), 30))
        vals[str(Rv)] = {"closed_form": str(v0), "neighbours": nb}
        print(tag, "R =", Rv, "f =", v0, flush=True)
        for k, v in nb.items():
            print("    Z_B + (%s): %s" % (k, v), flush=True)
    rec["values"] = vals
    if tag == "fam1_ZB6":
        try:
            q = hybrid_quadrature(F(3), P0, P0, S1, F(6), S1, 1.4)
            rec["quadrature_R_1.4"] = repr(q)
            print(tag, "quadrature at R = 1.4:", q, flush=True)
        except Exception as e:      # noqa: BLE001
            rec["quadrature_error"] = repr(e)
            print("quadrature failed:", repr(e), flush=True)
    print(tag, "atoms", rec["atoms"], "constants", rec["constants"], "lnR coeff:", rec["lnR_coefficient_of_f_at_0"],
          flush=True)
    out[tag] = rec
    (ROOT / "debug" / "data" / "minorder_coincidence_check.json").write_text(json.dumps(out, indent=1))
