"""Compare the operator derived from the closed form (pachucki_m_case_<tag>.json) with a committed operator_Q_crt.
Run: python debug/pachucki_m_opcompare.py TAG COMMITTED_JSON     (e.g. 2ppi_1_2 debug/data/ep_heldout_2ppi_1_2_EPS.json)"""
import json
import sys

import sympy as sp

tag, fn = sys.argv[1], sys.argv[2]
J = json.load(open(f"debug/data/pachucki_m_case_{tag}.json"))
C = json.load(open(fn))
R = sp.Symbol("R")
cop = [sp.Poly(sp.sympify(s_, locals={"R": R}), R) for s_ in C["operator_Q_crt"]]
mine = [sp.Poly(sum(sp.Integer(x) * R**j for j, x in enumerate(c)), R) if c else sp.Poly(0, R) for c in J["derived_operator"]]
res = {"committed_file": fn, "committed_order_degree": [C.get("order"), C.get("degree")],
       "derived_order": len(mine) - 1, "derived_degree": max(p.degree() for p in mine)}
if len(cop) == len(mine):
    ratio = sp.Rational(cop[-1].LC()) / sp.Rational(mine[-1].LC())
    res["scalar"] = str(ratio)
    res["equal_over_Q_up_to_scalar"] = bool(all(sp.expand(cop[i].as_expr() - ratio * mine[i].as_expr()) == 0
                                                for i in range(len(cop))))
else:
    res["equal_over_Q_up_to_scalar"] = False
print(res)
J["operator_comparison"] = res
json.dump(J, open(f"debug/data/pachucki_m_case_{tag}.json", "w"), indent=1)
