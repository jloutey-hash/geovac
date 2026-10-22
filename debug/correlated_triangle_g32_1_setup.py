"""GT stage 1: load L9 (gamma = 3/2) and G4 (gamma = 1); sanity: R2 divides L9 exactly; local
exponents of L9 at its singular points; adjoint convention check against QG's Q-form adjoint."""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from pathlib import Path
import sympy as sp
import flint
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G

u = G.u_sym
L9 = G.load_L9()
print("L9 order", len(L9) - 1, "x-degree", G.op_degree(L9))
rec = json.loads((G.DATA / "ring_ext_QG_6e_order2_g1.5.json").read_text())
R2 = G.from_Q([sp.sympify(q, locals={"u": u}) for q in rec["R2_Q"]])
t0 = time.time()
m, Qt, rem = G.pseudo_rdiv(L9, R2)
print("L9 / R2: pseudo steps", m, "remainder zero:", G.is_zero_op(rem), f"{time.time()-t0:.1f}s")
# adjoint check vs Q-form rule Q_i(u) -> Q_i(-u - 2i)
Q = [sp.sympify(q, locals={"u": u}) for q in json.loads((G.DATA / "ring_ext_QG_5c_ylog_g1.5_r9_d25.json").read_text())["Q"]]
Qa = [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(Q)]
A1 = G.from_Q(Qa)
A2 = G.adjoint(L9)
print("adjoint matches QG rule:", all((a - b).is_zero() for a, b in zip(A1, A2)))
print("L9** == L9:", all((a - b).is_zero() for a, b in zip(G.adjoint(A2), L9)))
print("t=0 exponents (t):", G.exps_of(G.indicial_t0(L9)))
print("inf growth exponents (t):", G.exps_of(G.indicial_inf(L9)))
for x0 in (Fraction(1, 4), Fraction(4), Fraction(25, 4), Fraction(9)):
    print("x0 =", x0, "exponents:", G.exps_of(G.indicial_at(L9, x0)))
G4 = G.load_G4()
print("G4: t0", G.exps_of(G.indicial_t0(G4)), "inf", G.exps_of(G.indicial_inf(G4)), "x=4", G.exps_of(G.indicial_at(G4, 4)))
lead = L9[-1]
print("L9 leading coefficient factors:", [(str(f), e) for f, e in lead.factor()[1]])
