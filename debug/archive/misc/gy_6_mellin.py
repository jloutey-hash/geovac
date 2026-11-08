"""GY stage 6: the Mellin route (ring_ext_QG_ops.mellin_Y_from_phi, record QG) from the kite's exact operator in z
(gy_5_op_kite.json, Q_theta_z) to an annihilator L_Y^M of Y_G(t) (an order bound), and the exact statement that the
Y_log operator L9G (gy_5_op_ylog.json) right-divides it: L_Y^M = Q o L9G exactly (pseudo-division over Q).  Since the
minimal annihilator of Y_log is a right factor of every annihilator of Y (monodromy around t = 0, QG memo), and L9G
annihilates Y itself on the grid (gy_9_onY), L9G is the minimal annihilator of Y_G over Q(x).
usage: python gy_6_mellin.py <kite op json> <ylog op json> <out tag>"""
from __future__ import annotations
import json, sys
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import ring_ext_QG_ops as O
import gy_newton as NW

DATA = Path(__file__).parent / "data"
u = G.u_sym
kite = json.loads(Path(sys.argv[1]).read_text())
ylog = json.loads(Path(sys.argv[2]).read_text())
Qphi = [sp.sympify(q, locals={"u": O.u}) for q in kite["Q_theta_z"]]
QY = O.mellin_Y_from_phi(Qphi)
LM = G.from_Q([sp.sympify(str(q), locals={"u": u}) for q in QY])
L9 = G.from_Q([sp.sympify(q, locals={"u": u}) for q in ylog["Q"]])
print("kite operator: order", max(sp.Poly(q, O.u).degree() for q in Qphi), "z-degree", len(Qphi) - 1, flush=True)
print("Mellin-route operator L_Y^M: order", len(LM) - 1, "x-degree", G.op_degree(LM), flush=True)
print("Y_log operator: order", len(L9) - 1, "x-degree", G.op_degree(L9), flush=True)
m, Qt, rem = G.pseudo_rdiv(LM, L9)
ok = G.is_zero_op(rem)
print("L9G right-divides L_Y^M exactly:", ok, " quotient order", len(Qt) - 1, flush=True)
nwM = NW.newton_infinity(LM)
nw9 = NW.newton_infinity(L9)
out = {"kite_order": max(sp.Poly(q, O.u).degree() for q in Qphi), "kite_zdegree": len(Qphi) - 1,
       "LM_order": len(LM) - 1, "LM_xdegree": G.op_degree(LM), "L9_order": len(L9) - 1, "L9_xdegree": G.op_degree(L9),
       "L9_right_divides_LM": bool(ok), "quotient_order": len(Qt) - 1,
       "LM_theta_form": [str(c) for c in LM], "LM_Q": [str(q) for q in QY],
       "newton_LM": nwM, "newton_L9": nw9}
(DATA / f"gy_6_mellin_{sys.argv[3]}.json").write_text(json.dumps(out, indent=1, default=str))
print("Newton at infinity: L_Y^M height", nwM["height"], "edges", [(e["slope"], e.get("exponents_c_of_exp_c_x_sigma")) for e in nwM["edges"]], flush=True)
print("Newton at infinity: L9G height", nw9["height"], "edges", [(e["slope"], e.get("exponents_c_of_exp_c_x_sigma")) for e in nw9["edges"]], flush=True)
