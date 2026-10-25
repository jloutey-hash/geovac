"""GA3 (copy of ga1_9_onY.py, output ga3_*): an operator (Q-form json, or a theta-form operator from a pieces json) applied to Y_lin itself at
gamma = 2 from the cached w-grid (ga3_3_grid.py), with the g32 1/200 perturbation control: the scalar
coefficient of largest magnitude in the theta^0 part, times (1 + 1/200).  Copy of debug/cf_g2_9_onY.py.
Residual = |L Y| / (largest t^(2i)-group), QG's convention; the largest single monomial is reported too;
quadrature h vs 2h.
usage: python ga1_9_onY.py <grid json> <op json> <key|Q> <out tag> <t list>"""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
import sympy as sp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import ring_ext_QG_8_onY as E

DATA = Path(__file__).parent / "data"


def load(opfile, key):
    d = json.loads(Path(opfile).read_text())
    if key == "Q":
        return G.from_Q([sp.sympify(q, locals={"u": G.u_sym}) for q in d["Q"]])
    return G.parse_op(d[key])


def apply(L, th, t, perturb=None):
    t2 = mp.mpf(t) ** 2
    groups, monos = [], []
    for i in range(G.op_degree(L) + 1):
        acc = mp.mpf(0)
        for k, c in enumerate(L):
            cc = c.coeffs()
            if i < len(cc) and int(cc[i]):
                v = mp.mpf(int(cc[i]))
                if perturb == (k, i):
                    v *= 1 + mp.mpf(1) / 200
                term = v * t2 ** i * th[k]
                acc += term
                monos.append(abs(term))
        groups.append(acc)
    tot = mp.fsum(groups)
    return tot, max(abs(g) for g in groups), max(monos)


def main():
    rec = json.loads(Path(sys.argv[1]).read_text())
    L = load(sys.argv[2], sys.argv[3])
    tag = sys.argv[4]
    ts = sys.argv[5].split(",")
    r = len(L) - 1
    cs = [int(v) for v in L[0].coeffs()]
    j = max(range(len(cs)), key=lambda i: abs(cs[i]))
    mp.mp.dps = rec["digits"] + 8
    out = {"op": sys.argv[2], "key": sys.argv[3], "order": r, "perturbed": {"theta_power": 0, "x_power": j, "factor": "1 + 1/200"}}
    for t in ts:
        th = E.thetaY(rec, t, r)
        th2 = E.thetaY(rec, t, r, every=2)
        v, bg, bm = apply(L, th, t)
        vb, bgb, _ = apply(L, th, t, perturb=(0, j))
        v2, bg2, _ = apply(L, th2, t)
        row = {"Y": mp.nstr(th[0], 25), "quad_h_vs_2h_Y": mp.nstr(abs(th2[0] / th[0] - 1), 3),
               "LY_rel_group": mp.nstr(abs(v) / bg, 3), "LY_rel_monomial": mp.nstr(abs(v) / bm, 3),
               "LY_rel_group_2h": mp.nstr(abs(v2) / bg2, 3),
               "perturbed_rel_group": mp.nstr(abs(vb) / bgb, 3)}
        out[t] = row
        print(t, row, flush=True)
    (DATA / f"ga3_9_onY_{tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
