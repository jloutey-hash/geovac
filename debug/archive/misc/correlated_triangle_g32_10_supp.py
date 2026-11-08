"""GT supplementary (not registered): (a) L9 with its largest theta^0 scalar coefficient times (1 + 1/200),
on Y at the GT3 points, to show the residual instrument separates an annihilator from a perturbed one;
(b) the split search on the gamma = 1 order-3 quotient (pipeline consistency with the control)."""
import json, sys
from pathlib import Path
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import correlated_triangle_g32_lib4 as S
import ring_ext_QG_8_onY as E
import importlib.util
spec = importlib.util.spec_from_file_location("onY", str(Path(__file__).parent / "correlated_triangle_g32_9_onY.py"))
rec = json.loads((G.DATA / "ring_ext_QG_3_grid_gem_g1.5_d32_h0.0625.json").read_text())
L9 = G.load_L9()
cs = [int(v) for v in L9[0].coeffs()]
j = max(range(len(cs)), key=lambda i: abs(cs[i]))
mp.mp.dps = rec["digits"] + 8
out = {}


def apply(L, th, t, perturb=None):
    t2 = mp.mpf(t) ** 2
    groups = []
    for i in range(G.op_degree(L) + 1):
        acc = mp.mpf(0)
        for k, c in enumerate(L):
            cc = c.coeffs()
            if i < len(cc) and int(cc[i]):
                v = mp.mpf(int(cc[i]))
                if perturb == (k, i):
                    v *= 1 + mp.mpf(1) / 200
                acc += v * t2 ** i * th[k]
        groups.append(acc)
    return abs(mp.fsum(groups)) / max(abs(g) for g in groups)


for t in ("0.8", "1.3", "1.7"):
    th = E.thetaY(rec, t, 9)
    out[t] = {"L9Y_rel": mp.nstr(apply(L9, th, t), 3), "L9_perturbed_rel": mp.nstr(apply(L9, th, t, (0, j)), 3)}
    print(t, out[t], flush=True)
c = json.loads((G.DATA / "correlated_triangle_g32_control_g1.json").read_text())
F1 = G.parse_op(c["F"])
out["gamma1_F_split"] = {"sol": S.split_search(F1, "sol"), "wronskian": S.split_search(F1, "wronskian")}
(G.DATA / "correlated_triangle_g32_10_supp.json").write_text(json.dumps(out, indent=1, default=str))
