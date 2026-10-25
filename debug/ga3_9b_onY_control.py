"""GA3 supplementary 1/200 control (docs/predictions/2026-10-23_geminal_shapes_g2.md, decision gate (c)): ga3_9_onY.py
perturbs the largest theta^0 coefficient of the operator, whose monomial is small at larger t; here, at each t, the
coefficient whose monomial |c_{k,i} t^(2i) theta^k Y| is LARGEST at that t is changed by 1/200 instead, so the control
measures the operator's sensitivity where it is actually supported.  Same grid, operator and residual convention.
usage: python ga3_9b_onY_control.py <grid json> <op json> <out tag> <t list>"""
from __future__ import annotations
import json, sys
from pathlib import Path
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import correlated_triangle_g32_lib as G
import ring_ext_QG_8_onY as E
import ga3_9_onY as Y9

DATA = Path(__file__).parent / "data"


def main():
    rec = json.loads(Path(sys.argv[1]).read_text())
    L = Y9.load(sys.argv[2], "Q")
    tag = sys.argv[3]
    ts = sys.argv[4].split(",")
    r = len(L) - 1
    mp.mp.dps = rec["digits"] + 8
    out = {"op": sys.argv[2], "order": r, "rule": "perturb by 1/200 the coefficient (theta^k, x^i) whose monomial is largest at this t"}
    for t in ts:
        th = E.thetaY(rec, t, r)
        t2 = mp.mpf(t) ** 2
        best, where = mp.mpf(0), None
        for k, c in enumerate(L):
            for i, cf in enumerate(c.coeffs()):
                if int(cf):
                    term = abs(mp.mpf(int(cf)) * t2 ** i * th[k])
                    if term > best:
                        best, where = term, (k, i)
        v, bg, bm = Y9.apply(L, th, t)
        vb, bgb, _ = Y9.apply(L, th, t, perturb=where)
        row = {"perturbed_theta_power": where[0], "perturbed_x_power": where[1],
               "LY_rel_group": mp.nstr(abs(v) / bg, 3), "perturbed_rel_group": mp.nstr(abs(vb) / bgb, 3)}
        out[t] = row
        print(t, row, flush=True)
    (DATA / f"ga3_9b_onY_control_{tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
