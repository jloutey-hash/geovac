"""GT3: K6 (the right factor of L9 below the order-3 piece, L9 = F o K6) applied to Y itself at
t = 0.8, 1.3, 1.7 from QG's cached 32-digit w-grid (gamma = 3/2), relative to the largest term
(QG's convention: the terms are the x-power groups x^i sum_k [x^i]c_k theta^k Y; the largest single
monomial is reported too).  Reference: L9 itself on the same grid.
Negative control (fixed before running): the scalar coefficient of largest magnitude in K6's theta^0
coefficient polynomial is multiplied by (1 + 1/200)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
import flint
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import ring_ext_QG_8_onY as E

rec = json.loads((G.DATA / "ring_ext_QG_3_grid_gem_g1.5_d32_h0.0625.json").read_text())
d = json.loads((G.DATA / "correlated_triangle_g32_series_final.json").read_text())
K6 = G.parse_op(d["chain_K6"])
L9 = G.load_L9()
cs = [int(v) for v in K6[0].coeffs()]
j = max(range(len(cs)), key=lambda i: abs(cs[i]))
mp.mp.dps = rec["digits"] + 8


def apply(L, th, t, perturb=None):
    t2 = mp.mpf(t) ** 2
    D = G.op_degree(L)
    groups, monos = [], []
    for i in range(D + 1):
        acc = mp.mpf(0)
        for k, c in enumerate(L):
            cc = c.coeffs()
            if i < len(cc) and int(cc[i]):
                v = mp.mpf(int(cc[i]))
                if perturb is not None and (k, i) == perturb:
                    v = v * (1 + mp.mpf(1) / 200)
                term = v * t2 ** i * th[k]
                acc += term
                monos.append(abs(term))
        groups.append(acc)
    tot = mp.fsum(groups)
    return tot, max(abs(g) for g in groups), max(monos)


out = {"perturbed_coefficient": {"theta_power": 0, "x_power": j, "value": str(cs[j]), "factor": "1 + 1/200"}}
for t in ("0.8", "1.3", "1.7"):
    th = E.thetaY(rec, t, 9)
    th2 = E.thetaY(rec, t, 9, every=2)
    k6, bg, bm = apply(K6, th, t)
    k6b, bgb, bmb = apply(K6, th, t, perturb=(0, j))
    l9, lg, lm = apply(L9, th, t)
    # same quantities from the h vs 2h grid (quadrature self-check)
    k62, bg2, _ = apply(K6, th2, t)
    row = {"Y": mp.nstr(th[0], 25), "quad_selfcheck_Y": mp.nstr(abs(th2[0] / th[0] - 1), 3),
           "K6Y": mp.nstr(k6, 15), "K6Y_rel_group": mp.nstr(abs(k6) / bg, 6), "K6Y_rel_monomial": mp.nstr(abs(k6) / bm, 6),
           "K6Y_h_vs_2h": mp.nstr(abs(k62 / k6 - 1), 3),
           "control_K6bad_rel_group": mp.nstr(abs(k6b) / bgb, 6), "control_change_rel": mp.nstr(abs(k6b - k6) / abs(k6), 4),
           "L9Y_rel_group": mp.nstr(abs(l9) / lg, 3)}
    out[t] = row
    print(t, row, flush=True)
(G.DATA / "correlated_triangle_g32_9_onY.json").write_text(json.dumps(out, indent=1))
