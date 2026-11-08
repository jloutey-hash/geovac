"""MF step B: control C3 (single-centre H2+ energy) and the registered S_mol(R).

Record: docs/predictions/2026-10-20_molecular_fns.md (committed 01f79fa5 before computing).
S_mol(R) = -( lim_{rn -> 0} [ Q_mol(rn) - Q_H1s(rn) ] ), both on the same grid, same rn.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import molfns_solver as S  # noqa: E402
from molfns_A_controls import mixed_Q  # noqa: E402

OUT = Path(__file__).parent / "data" / "molfns_B.json"
E_PROLATE_R2 = -1.1026342144949


def run(R0: float, lmax: int, h: float, rns, t0: float = -40.0, rmax: float = 60.0):
    g = S.Grid(t0=t0, t1=np.log(rmax), h=h, anchor=R0)
    R = g.snap(R0)
    mol = S.System(grid=g, lmax=lmax, Zd=1.0, others=[(1.0, R)])
    t = time.time()
    vals, vecs = S.ground_state(mol, sigma=-1.0 / R - 0.7 if R < 3 else -0.75, k=1)
    E, psi = float(vals[0]), vecs[:, 0]
    atom = S.System(grid=g, lmax=0, Zd=1.0)
    av, avec = S.ground_state(atom, sigma=-0.55, k=1)
    out = {"R": R, "lmax": lmax, "h": h, "E": E, "rows": {}}
    for rn0 in rns:
        rn = g.snap(rn0)
        qm = mixed_Q(mol, E, psi, rn)
        qa = mixed_Q(atom, float(av[0]), avec[:, 0], rn)
        out["rows"][f"{rn:.3e}"] = {"Qmol": qm, "QH": qa, "S": -(qm - qa)}
    out["sec"] = time.time() - t
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "c3"
    rec = {}
    if mode == "c3":
        for (lmax, h) in [(16, 0.02), (16, 0.01), (16, 0.005), (32, 0.005)]:
            o = run(2.0, lmax, h, [1e-3], t0=-30.0, rmax=40.0)
            print(f"R={o['R']:.6f} lmax={lmax} E={o['E']:.10f} dE={o['E']-E_PROLATE_R2:.2e} "
                  f"S(1e-3)={list(o['rows'].values())[0]['S']:.5f} {o['sec']:.0f}s", flush=True)
            rec[f"c3_l{lmax}"] = o
    if mode == "prod":
        plan = [(2.0, [16, 24, 32, 48, 64]), (1.0, [16, 24, 32, 48]), (4.0, [16, 24, 32, 48]),
                (8.0, [16, 24, 32, 48, 64])]
        for R0, lms in plan:
            for lmax in lms:
                o = run(R0, lmax, 0.005, [1e-3, 3e-4, 1e-4], t0=-30.0, rmax=45.0)
                Ss = [v["S"] for v in o["rows"].values()]
                print(f"R={o['R']:.4f} lmax={lmax} E={o['E']:.10f} S(rn)={['%.5f' % x for x in Ss]} "
                      f"{o['sec']:.0f}s", flush=True)
                rec[f"R{R0}_l{lmax}"] = o
                OUT.with_name("molfns_B_prod.json").write_text(json.dumps(rec, indent=1))
    if mode == "c3h":
        for lmax in (32, 48):
            o = run(2.0, lmax, 0.0025, [1e-3, 3e-4, 1e-4], t0=-30.0, rmax=45.0)
            Ss = [v["S"] for v in o["rows"].values()]
            print(f"R={o['R']:.4f} lmax={lmax} h=0.0025 E={o['E']:.10f} S(rn)={['%.5f' % x for x in Ss]} "
                  f"{o['sec']:.0f}s", flush=True)
            rec[f"R2_l{lmax}_h0.0025"] = o
            OUT.with_name("molfns_B_c3h.json").write_text(json.dumps(rec, indent=1))
    OUT.with_name(f"molfns_B_{mode}.json").write_text(json.dumps(rec, indent=1))
