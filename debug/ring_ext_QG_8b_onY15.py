"""QG step 8b (gamma = 3/2): the exact order-9 operator L9 (minimal annihilator of Y_log, 5d) applied to
Y itself at t = 0.8, 1.3, 1.7, 4.0 from the cached 32-digit Phi~ w-grid; plus a control with one
coefficient of L9 perturbed."""
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_8_onY as E

DATA = Path(__file__).parent / "data"
u = sp.Symbol("u")
rec = json.loads((DATA / "ring_ext_QG_3_grid_gem_g1.5_d32_h0.0625.json").read_text())
Q = [sp.sympify(q, locals={"u": u}) for q in json.loads((DATA / "ring_ext_QG_5c_ylog_g1.5_r9_d25.json").read_text())["Q"]]
Qbad = list(Q)
Qbad[3] = sp.expand(Q[3] + sp.Poly(Q[3], u).LC() // 10**6 * u**2)
mp.mp.dps = rec["digits"] + 8
out = {}
for t in ("0.8", "1.3", "1.7", "4.0"):
    th = E.thetaY(rec, t, 9)
    th2 = E.thetaY(rec, t, 9, every=2)
    tot, big = E.apply_Q(Q, th, t)
    tb, bb = E.apply_Q(Qbad, th, t)
    out[t] = {"Y": mp.nstr(th[0], 30), "selfcheck_h_vs_2h": mp.nstr(abs(th2[0] / th[0] - 1), 3),
              "L9Y_rel": mp.nstr(abs(tot) / big, 3), "control_rel": mp.nstr(abs(tb) / bb, 3)}
    print(t, out[t], flush=True)
(DATA / "ring_ext_QG_8b_onY_g1.5.json").write_text(json.dumps(out, indent=1))
