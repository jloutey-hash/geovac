"""MG3: local data of the kite's operator in z (Q form sum_i z^i Q_i(theta_z) from mg_5_guess.py, kind phi):
rational roots of the leading coefficient with multiplicities, the non-linear irreducible factors (apparent
candidates) with their degrees and the mod-p apparent-point check (all exponents non-negative integers, no
logs), the local exponents at z = 0, at every rational root and at infinity, and the monodromy classes.

Two routes.  (A) the g32 library: the operator P_i(theta_z) is rewritten as P_i(u/2) 2^r in u = theta_t
(theta_t = 2 x d/dx with x = z), so hform/indicial_at give the h-exponents at finite z0 directly; at z = 0 and
infinity the u-roots are divided by 2.  (B) sympy: the D-form (ring_ext_QG_ops.theta_to_D) and its indicial
polynomials, as in the QG memo.  Both are written; they must agree.
usage: python mg_kite_local.py <kite op json> <gamma> <out tag>"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib3 as Lc
import cf_g2_lib as C
import ring_ext_QG_ops as O

u, z, a = O.u, O.z, O.a
PF, GAMMA, TAG = sys.argv[1], Fraction(sys.argv[2]), sys.argv[3]
rec = json.loads(Path(PF).read_text())
Q = [sp.sympify(q, locals={"u": u}) for q in rec["Q"]]
r = max(sp.Poly(q, u).degree() for q in Q)
Qu = [sp.expand(q.subs(u, u / 2) * 2**r) for q in Q]
L = G.from_Q(Qu)
out = {"gamma": str(GAMMA), "order": r, "z_degree": G.op_degree(L),
       "predicted_points": sorted({str(-v * v) for v in (Fraction(2), 2 * GAMMA, 1 + GAMMA, abs(1 - GAMMA))} | {"0"})}
lead = C.lead_factors(L)
out["leading_factors"] = [(str(f), e, f.degree()) for f, e in lead]
rr = C.rational_roots(L)
out["rational_roots"] = [(str(v), e) for v, e in rr]
print("order", r, "z-degree", out["z_degree"], flush=True)
print("leading coefficient factors (poly, multiplicity, degree):", [(s[:60], e, dg) for s, e, dg in out["leading_factors"]], flush=True)
# route A
halve = lambda es: [str(Fraction(e) / 2) for e in es]
out["A"] = {"z0_exponents": halve(Lc.exps_t0(L)[0]), "inf_exponents_z": halve(Lc.exps_inf_growth(L)[0]),
            "z0_classes": [c["exponents"] for c in Lc.monodromy_type(L, Fraction(0))],
            "inf_classes": [c["exponents"] for c in Lc.monodromy_type(L, "inf")], "points": {}}
for z0, e in rr:
    if z0 == 0:
        continue
    out["A"]["points"][str(z0)] = {"exps": [str(v) for v in Lc.exps_finite_t(L, z0)[0]],
                                   "classes": [c["exponents"] for c in Lc.monodromy_type(L, z0)]}
out["apparent_check"] = {}
for f, e in lead:
    if f.degree() >= 2:
        ac = C.apparent_check(L, f, nprimes=2)
        out["apparent_check"][str(f)[:50]] = {"degree": f.degree(), "multiplicity": e, "check": ac}
        print("apparent check (degree %d factor):" % f.degree(), ac, flush=True)
print("route A: z=0 exponents", out["A"]["z0_exponents"], "; infinity", out["A"]["inf_exponents_z"], flush=True)
for k, v in out["A"]["points"].items():
    print(f"   z = {k}: exponents {v['exps']}; classes {v['classes']}", flush=True)
# route B
A = O.theta_to_D(Q, var=z, step=1)
out["B"] = {"z0": str(O.indicial(A, 0)), "inf": str(O.indicial_inf(A)), "points": {}}
for z0, e in rr:
    if z0 != 0:
        out["B"]["points"][str(z0)] = str(O.indicial(A, sp.Rational(z0.numerator, z0.denominator)))
print("route B indicial at 0:", out["B"]["z0"], flush=True)
print("route B indicial at infinity (f ~ z^-a):", out["B"]["inf"], flush=True)
for k, v in out["B"]["points"].items():
    print(f"   route B z = {k}: {v}", flush=True)
(C.DATA / f"mg_kite_local_{TAG}.json").write_text(json.dumps(out, indent=1, default=str))
