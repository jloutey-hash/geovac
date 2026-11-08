"""PF2 discriminator (coordinator addition, 2026-10-23): local data of an operator in x = t^2 (theta-form json, key, or
Q-form json with key 'Q') at the direct-pairing thresholds x = 25/4 (t = 1 + 3/2) and x = 25 (t = 2 + 3), which are
absent under the crossed reading, and at the eight predicted points; genuine = a non-integer or negative exponent or a
forced logarithm; regular = not a root of the leading coefficient.

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).
usage: python pf_discriminator.py <op json> <key|Q> <out tag>"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib3 as Lc
import cf_g2_lib as C

d = json.loads(Path(sys.argv[1]).read_text())
L = G.from_Q([sp.sympify(q, locals={"u": G.u_sym}) for q in d["Q"]]) if sys.argv[2] == "Q" else G.parse_op(d[sys.argv[2]])
pts = {"discriminator": [Fraction(25, 4), Fraction(25)],
       "predicted": [Fraction(v) for v in ("1/4", "1", "9/4", "4", "9", "49/4", "16", "81/4")]}
rr = dict(C.rational_roots(L))
out = {"op": sys.argv[1], "key": sys.argv[2], "order": len(L) - 1, "x_degree": G.op_degree(L),
       "rational_roots_of_lead": {str(k): v for k, v in rr.items()}}
for grp, lst in pts.items():
    out[grp] = {}
    for x0 in lst:
        if rr.get(x0, 0) == 0:
            out[grp][str(x0)] = {"lead_multiplicity": 0, "verdict": "regular (not a root of the leading coefficient)"}
            continue
        exps = Lc.exps_finite_t(L, x0)[0]
        cl = Lc.monodromy_type(L, x0)
        nonint = any(Fraction(e).denominator != 1 or Fraction(e) < 0 for e in exps)
        logs = any(not ok for c in cl for _, _, ok in c["exponents"])
        out[grp][str(x0)] = {"lead_multiplicity": rr[x0], "exps": [str(e) for e in exps], "classes": [c["exponents"] for c in cl],
                             "verdict": "GENUINE" if (nonint or logs) else "apparent-type (non-negative integer, log-free)"}
others = [x0 for x0 in rr if x0 != 0 and x0 not in pts["discriminator"] and x0 not in pts["predicted"]]
out["other_rational_roots"] = {}
for x0 in others:
    exps = Lc.exps_finite_t(L, x0)[0]
    cl = Lc.monodromy_type(L, x0)
    nonint = any(Fraction(e).denominator != 1 or Fraction(e) < 0 for e in exps)
    logs = any(not ok for c in cl for _, _, ok in c["exponents"])
    out["other_rational_roots"][str(x0)] = {"exps": [str(e) for e in exps], "verdict": "GENUINE" if (nonint or logs) else "apparent-type"}
for k in ("discriminator", "predicted", "other_rational_roots"):
    print(k, json.dumps(out[k]), flush=True)
(Path(__file__).parent / "data" / f"pf_discriminator_{sys.argv[3]}.json").write_text(json.dumps(out, indent=1, default=str))
