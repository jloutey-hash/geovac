"""PF2 kite clause: local data of the off-slice kite's operator in z (Q form sum_i z^i Q_i(theta_z) from pf_5_guess.py,
kind phi): rational roots of the leading coefficient with multiplicities, the irreducible factors of degree >= 2
(apparent candidates) with the mod-p apparent-point check (all exponents non-negative integers, no logs), the local
exponents at z = 0, at every rational root and at infinity, and the monodromy classes.  Copy of debug/mg_kite_local.py
with the predicted points of the record (z = minus the eight squares) and a verdict per point, plus the direct-pairing
discriminator z = -25/4, -25 (coordinator, 2026-10-23) reported explicitly.

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).

Route A: the g32 library with P_i(theta_z) rewritten as P_i(u/2) 2^r in u = theta_t (theta_t = 2 x d/dx, x = z), so
hform/indicial_at give the h-exponents at finite z0 directly; at z = 0 and infinity the u-roots are divided by 2.
Route B: sympy D-form indicial polynomials (ring_ext_QG_ops).  Both written; they must agree.
usage: python pf_kite_local.py <kite op json> <out tag>"""
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

u, z = O.u, O.z
PF, TAG = sys.argv[1], sys.argv[2]
PRED = sorted({-Fraction(v) for v in ("1/4", "1", "9/4", "4", "9", "49/4", "16", "81/4")})
# the direct-pairing thresholds (1 + 3/2)^2 and (2 + 3)^2, absent under the crossed reading: the one-point discriminator
DISCRIM = [-Fraction(25, 4), -Fraction(25)]
rec = json.loads(Path(PF).read_text())
Q = [sp.sympify(q, locals={"u": u}) for q in rec["Q"]]
r = max(sp.Poly(q, u).degree() for q in Q)
Qu = [sp.expand(q.subs(u, u / 2) * 2**r) for q in Q]
L = G.from_Q(Qu)
out = {"order": r, "z_degree": G.op_degree(L), "predicted_points": [str(v) for v in PRED] + ["0", "inf"]}
lead = C.lead_factors(L)
out["leading_factors"] = [(str(f), e, f.degree()) for f, e in lead]
rr = C.rational_roots(L)
out["rational_roots"] = [(str(v), e) for v, e in rr]
print("order", r, "z-degree", out["z_degree"], flush=True)
print("leading coefficient factors (poly, multiplicity, degree):", [(s[:60], e, dg) for s, e, dg in out["leading_factors"]], flush=True)
halve = lambda es: [str(Fraction(e) / 2) for e in es]
out["A"] = {"z0_exponents": halve(Lc.exps_t0(L)[0]), "inf_exponents_z": halve(Lc.exps_inf_growth(L)[0]),
            "z0_classes": [c["exponents"] for c in Lc.monodromy_type(L, Fraction(0))],
            "inf_classes": [c["exponents"] for c in Lc.monodromy_type(L, "inf")], "points": {}}
verdict = {}
for z0 in PRED:
    mult = dict((v, e) for v, e in rr).get(z0, 0)
    if mult == 0:
        verdict[str(z0)] = "regular (not a root of the leading coefficient)"
        continue
    exps = Lc.exps_finite_t(L, z0)[0]
    classes = Lc.monodromy_type(L, z0)
    out["A"]["points"][str(z0)] = {"lead_multiplicity": mult, "exps": [str(v) for v in exps], "classes": [c["exponents"] for c in classes]}
    nonint = any(Fraction(e).denominator != 1 or Fraction(e) < 0 for e in exps)
    logs = any(not ok for c in classes for _, _, ok in c["exponents"])
    verdict[str(z0)] = "genuine" if (nonint or logs) else "apparent-type exponents (non-negative integers, log-free)"
for z0, e in rr:
    if z0 not in PRED and z0 != 0:
        exps = Lc.exps_finite_t(L, z0)[0]
        classes = Lc.monodromy_type(L, z0)
        out["A"]["points"][str(z0)] = {"lead_multiplicity": e, "exps": [str(v) for v in exps], "classes": [c["exponents"] for c in classes], "outside_prediction": True}
        nonint = any(Fraction(x).denominator != 1 or Fraction(x) < 0 for x in exps)
        logs = any(not ok for c in classes for _, _, ok in c["exponents"])
        verdict[str(z0)] = ("GENUINE OUTSIDE THE PREDICTED SET" if (nonint or logs) else "apparent (rational root outside the set, integer log-free exponents)")
out["verdict_by_point"] = verdict
disc = {}
for z0 in DISCRIM:
    mult = dict((v, e) for v, e in rr).get(z0, 0)
    if mult == 0:
        disc[str(z0)] = {"lead_multiplicity": 0, "verdict": "regular point (not a root of the leading coefficient)"}
        continue
    exps = Lc.exps_finite_t(L, z0)[0]
    classes = Lc.monodromy_type(L, z0)
    nonint = any(Fraction(e).denominator != 1 or Fraction(e) < 0 for e in exps)
    logs = any(not ok for c in classes for _, _, ok in c["exponents"])
    disc[str(z0)] = {"lead_multiplicity": mult, "exps": [str(e) for e in exps], "classes": [c["exponents"] for c in classes],
                     "verdict": "GENUINE" if (nonint or logs) else "apparent-type (non-negative integer, log-free exponents)"}
out["direct_pairing_discriminator"] = disc
print("direct-pairing discriminator (z = -25/4, -25):", disc, flush=True)
out["apparent_check"] = {}
for f, e in lead:
    if f.degree() >= 2:
        ac = C.apparent_check(L, f, nprimes=2)
        out["apparent_check"][str(f)[:50]] = {"degree": f.degree(), "multiplicity": e, "check": ac}
        print("apparent check (degree %d factor):" % f.degree(), ac, flush=True)
print("route A: z=0 exponents", out["A"]["z0_exponents"], "; infinity", out["A"]["inf_exponents_z"], flush=True)
for k, v in out["A"]["points"].items():
    print(f"   z = {k}: exponents {v['exps']}; classes {v['classes']}", flush=True)
print("verdicts:", verdict, flush=True)
A = O.theta_to_D(Q, var=z, step=1)
out["B"] = {"z0": str(O.indicial(A, 0)), "inf": str(O.indicial_inf(A)), "points": {}}
for z0, e in rr:
    if z0 != 0:
        out["B"]["points"][str(z0)] = str(O.indicial(A, sp.Rational(z0.numerator, z0.denominator)))
print("route B indicial at 0:", out["B"]["z0"], flush=True)
print("route B indicial at infinity (f ~ z^-a):", out["B"]["inf"], flush=True)
for k, v in out["B"]["points"].items():
    print(f"   route B z = {k}: {v}", flush=True)
(C.DATA / f"pf_kite_local_{TAG}.json").write_text(json.dumps(out, indent=1, default=str))
