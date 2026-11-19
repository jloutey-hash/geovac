"""Compare derived operators with the operators frozen in the committed tests (read as text; nothing imported).
Run: python debug/pachucki_m_opcompare_tests.py"""
import ast
import json
import re

import sympy as sp

R = sp.Symbol("R")
res = {}
for tag, fn, name in (("2ppi_1_2", "tests/test_paper58_exchange_l1.py", "OP_P_1_2"),
                      ("3ddel_1_2", "tests/test_paper58_exchange_d_pp.py", "OP_3DDEL_1_2")):
    src = open(fn, encoding="utf-8").read()
    m = re.search(r"^%s = (\[.*?^\])" % name, src, re.S | re.M)
    lst = ast.literal_eval(m.group(1))
    cop = [sp.Poly(sp.sympify(s_, locals={"R": R}), R) for s_ in lst]
    try:
        J = json.load(open(f"debug/data/pachucki_m_case_{tag}.json"))
        mine = [sp.Poly(sum(sp.Integer(x) * R**j for j, x in enumerate(c)), R) for c in J["derived_operator"]]
    except (OSError, KeyError) as e:
        res[tag] = "no derived operator: %r" % e
        continue
    ratio = sp.Rational(cop[-1].LC()) / sp.Rational(mine[-1].LC())
    eq = len(cop) == len(mine) and all(sp.expand(cop[i].as_expr() - ratio * mine[i].as_expr()) == 0 for i in range(len(cop)))
    res[tag] = {"test_file": fn, "name": name, "order": len(cop) - 1, "degree": max(p.degree() for p in cop),
                "scalar": str(ratio), "equal_over_Q_up_to_scalar": bool(eq)}
    print(tag, res[tag])
json.dump(res, open("debug/data/pachucki_m_opcompare_tests.json", "w"), indent=1)
