"""GY stage 9 (ga3_9_onY.py's functions, imported, on the Gaussian grid): an operator (Q-form json, or a theta-form key of a
pieces json) applied to Y_G itself from the cached w-grid (gy_3_grid.py), with the g32 1/200 perturbation control (the
scalar coefficient of largest magnitude in the theta^0 part, times 1 + 1/200), residual relative to the largest
t^(2i)-group (QG convention) and to the largest single monomial; quadrature h vs 2h.  Then the second normalisation
(ga3_9b): at each t the coefficient whose monomial is LARGEST at that t is perturbed by 1/200 instead.
usage: python gy_9_onY.py <grid json> <op json> <key|Q|theta_form> <out tag> <t list>"""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import ring_ext_QG_8_onY as E
import ga3_9_onY as Y9

DATA = Path(__file__).parent / "data"


def load(opfile, key):
    if key == "theta_form":
        d = json.loads(Path(opfile).read_text())
        return G.parse_op(d["theta_form"])
    return Y9.load(opfile, key)


def main():
    rec = json.loads(Path(sys.argv[1]).read_text())
    L = load(sys.argv[2], sys.argv[3])
    tag = sys.argv[4]
    ts = sys.argv[5].split(",")
    r = len(L) - 1
    cs = [int(v) for v in L[0].coeffs()]
    j = max(range(len(cs)), key=lambda i: abs(cs[i]))
    mp.mp.dps = rec["digits"] + 8
    out = {"op": sys.argv[2], "key": sys.argv[3], "order": r, "x_degree": G.op_degree(L),
           "perturbed_A": {"theta_power": 0, "x_power": j, "factor": "1 + 1/200", "rule": "largest theta^0 coefficient"},
           "perturbed_B_rule": "at each t, the coefficient (theta^k, x^i) whose monomial |c t^(2i) theta^k Y| is largest, times 1 + 1/200"}
    for t in ts:
        th = E.thetaY(rec, t, r)
        th2 = E.thetaY(rec, t, r, every=2)
        v, bg, bm = Y9.apply(L, th, t)
        vb, bgb, _ = Y9.apply(L, th, t, perturb=(0, j))
        v2, bg2, _ = Y9.apply(L, th2, t)
        t2 = mp.mpf(t) ** 2
        best, where = mp.mpf(0), None
        for k, c in enumerate(L):
            for i, cf in enumerate(c.coeffs()):
                if int(cf):
                    term = abs(mp.mpf(int(cf)) * t2 ** i * th[k])
                    if term > best:
                        best, where = term, (k, i)
        vc, bgc, _ = Y9.apply(L, th, t, perturb=where)
        row = {"Y": mp.nstr(th[0], 25), "quad_h_vs_2h_Y": mp.nstr(abs(th2[0] / th[0] - 1), 3),
               "LY_rel_group": mp.nstr(abs(v) / bg, 3), "LY_rel_monomial": mp.nstr(abs(v) / bm, 3),
               "LY_rel_group_2h": mp.nstr(abs(v2) / bg2, 3),
               "perturbed_A_rel_group": mp.nstr(abs(vb) / bgb, 3),
               "perturbed_B_where": where, "perturbed_B_rel_group": mp.nstr(abs(vc) / bgc, 3)}
        out[t] = row
        print(t, row, flush=True)
    (DATA / f"gy_9_onY_{tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
