"""PF1-PF4 stage 3: the off-slice kite Phi~(w) on the QG exp-sinh w-grid at 32 digits (copy of debug/ga3_3_grid.py
with pf_common.phi_ts), writing data/pf_grid_d<digits>_h<h>.json.  Y and theta^j Y then come from
ring_ext_QG_3_values.Y_theta / ring_ext_QG_8_onY.thetaY on this grid (the kernel part is shape-independent).

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).
usage: python pf_3_grid.py <digits> <htau> [nproc]"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
from flint import arb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_3_values as V
import pf_common as P

DATA = Path(__file__).parent / "data"


def phi_adaptive(args):
    w, digits = args
    ctx.dps = digits + 12
    w = arb(w)
    h = arb(1) / 16
    prev = P.phi_ts(w, h, digits + 10).real
    while True:
        h = h / 2
        cur = P.phi_ts(w, h, digits + 10).real
        if abs(cur - prev) < abs(cur) * arb(10) ** (-digits) or h < arb(1) / 1024:
            return cur.mid().str(digits + 8, radius=False), str(h)
        prev = cur


def main():
    digits, htau = int(sys.argv[1]), sys.argv[2]
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 15
    fn = DATA / f"pf_grid_d{digits}_h{htau}.json"
    grid = V.wgrid(mp.mpf(htau), digits)
    keep = [(w, dw) for w, dw in grid if w ** 3 * dw > mp.mpf(10) ** (-digits - 4) and dw / w ** 5 > mp.mpf(10) ** (-digits - 4)]
    tasks = [(mp.nstr(w, digits + 8), digits) for w, dw in keep]
    t0 = time.time()
    with Pool(nproc) as pool:
        res = pool.map(phi_adaptive, tasks, chunksize=1)
    rec = {"kind": "pf off-slice exp kite (a=1, b=2, gamma_1=3/2, gamma_2=3)", "gamma": "3/2,3", "digits": digits, "htau": htau,
           "w": [mp.nstr(w, digits + 8) for w, _ in keep], "dw": [mp.nstr(dw, digits + 8) for _, dw in keep],
           "phi": [r[0] for r in res], "h_used": [r[1] for r in res], "seconds": time.time() - t0}
    fn.write_text(json.dumps(rec, indent=0))
    print("nodes", len(keep), "seconds", round(rec["seconds"], 1), flush=True)
    for t in ("0.8", "1.3", "1.7", "2.5"):
        print("t", t, "Y =", mp.nstr(V.Y_theta(rec, t, 0)[0], digits), flush=True)


if __name__ == "__main__":
    main()
