"""PF1-PF4: residual-versus-degree scan of the kite-operator guess on the N = 300 direct series (460x460 substituted
nodes, Horner kernel verified to 8e-300), to separate genuine annihilators (a sharp drop to the data's precision,
~1e-216 at the tail) from threshold solutions.  ring_ext_QG_guess.try_operator, scale 9 (the Taylor radius).

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).
usage: python pf_5s_scan.py"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import pf_5_guess as PG
import ring_ext_QG_guess as Gs

SER = str(Path(__file__).parent / "data" / "pf_series_direct_explore_N300_d280_460x460_sub.json")
N, DPS = 300, 280


def job(args):
    r, d = args
    c, lam, step = PG.series(SER, "phi", N, DPS)
    out = Gs.try_operator(c, lam, step, r, d, n_to=N, scale=9, tol_digits=150)
    return r, d, (r + 1) * (d + 1), (out[1].mid().str(3, radius=False) if out else None)


if __name__ == "__main__":
    t0 = time.time()
    tasks = [(r, d) for r in (6, 9, 12) for d in range(8, 40, 2) if (r + 1) * (d + 1) <= N - 40]
    with Pool(12) as pool:
        res = pool.map(job, tasks, chunksize=1)
    out = {"series": SER, "N": N, "dps": DPS, "scale": 9, "rows": []}
    for r, d, U, v in res:
        print(f"r={r} d={d} U={U} worst check-row residual {v}", flush=True)
        out["rows"].append({"r": r, "d": d, "U": U, "residual": v})
    out["seconds"] = time.time() - t0
    (Path(__file__).parent / "data" / "pf_5s_scan_N300.json").write_text(json.dumps(out, indent=1))
    print(f"({out['seconds']:.0f}s)", flush=True)
