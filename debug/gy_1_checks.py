"""GY stage 1: route A (momentum triangle, closed-form Gaussian loop, 2-D quadrature) against route B (kite by the
convolution theorem, 1-D in u) at three w, gate 1e-14; the exponential-shape control of route A's skeleton
(ga1_common.phi_routeA 'exp' vs ring_ext_QG_common.phi_gem, lineage code, unmodified); phi_ts (tanh-sinh, real w)
against phi_z (Gauss-Chebyshev) at two w including a large one; E phi_a + phi_b = phi on the grid form.
Output debug/data/gy_1_checks.json"""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import mpmath as mp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm
import ga1_common as GA
import gy_common as GY

OUT = Path(__file__).parent / "data" / "gy_1_checks.json"
rec = {}
worst = 0.0
t0 = time.time()
# (d) route A vs route B, three w, 1e-14
mp.mp.dps = 22
ctx.dps = 40
rows = []
for w in ("0.5", "1.3", "3.0"):
    a = GY.phi_routeA(w, "1")
    b = GY.phi_z(acb(w) * acb(w), arb(1), 400, "full").real
    rel = float(abs(a / Cm.to_mp(b) - 1))
    rows.append({"w": w, "routeA": mp.nstr(a, 20), "routeB": b.mid().str(25, radius=False), "rel": rel})
    print(f"(d) Gaussian kite w = {w}: route A {mp.nstr(a, 18)}  route B {b.mid().str(18, radius=False)}  rel {rel:.2e}  ({time.time()-t0:.0f}s)", flush=True)
    worst = max(worst, rel)
rec["routeA_vs_routeB"] = rows
rec["routeA_vs_routeB_gate_1e-14"] = "PASS" if max(r["rel"] for r in rows) < 1e-14 else "FAIL"
# exp-shape control of the route-A skeleton (lineage code)
mp.mp.dps = 15
ctx.dps = 30
a = GA.phi_routeA("1.3", "1", "exp"); b = Cm.phi_gem(arb("1.3"), arb(1), 80, 80).real
rel = float(abs(a / Cm.to_mp(b) - 1)); rec["exp_kite_routeA_control"] = {"w": "1.3", "gamma": "1", "rel": rel}
print(f"(d) exp kite control (lineage code) w = 1.3: rel {rel:.2e}", flush=True)
# (e) tanh-sinh form vs Chebyshev form, and the split on the tanh-sinh form
ctx.dps = 44
rows = []
for w in ("1.7", "6.0"):
    a = GY.phi_ts(arb(w), arb(1), arb(1) / 32, 42, "full").real
    b = GY.phi_z(acb(w) * acb(w), arb(1), 600, "full").real
    rel = float(abs(a - b) / abs(b))
    E = GY.E_const(arb(1))
    pa = GY.phi_ts(arb(w), arb(1), arb(1) / 32, 42, "a").real
    pb = GY.phi_ts(arb(w), arb(1), arb(1) / 32, 42, "b").real
    rel2 = float(abs(E * pa + pb - a) / abs(a))
    rows.append({"w": w, "phi_ts": a.mid().str(30, radius=False), "phi_z": b.mid().str(30, radius=False), "rel": rel, "split_rel": rel2})
    print(f"(e) phi_ts vs phi_z at w = {w}: rel {rel:.2e}; E phi_a + phi_b = phi: {rel2:.2e}", flush=True)
    worst = max(worst, rel, rel2)
rec["phi_ts_vs_phi_z"] = rows
# large-w sanity of the tanh-sinh form (the grid reaches w ~ 1e7): finite, decays like exp(-w^2/8)
vals = {}
for w in ("10", "30", "1e3", "1e6"):
    v = GY.phi_ts(arb(w), arb(1), arb(1) / 16, 42, "full").real
    vals[w] = v.mid().str(12, radius=False) if v.is_finite() else "NOT FINITE"
rec["large_w_values"] = vals
print("(f) large w:", vals, flush=True)
rec["worst"] = worst
print("WORST", worst, flush=True)
OUT.write_text(json.dumps(rec, indent=1))
