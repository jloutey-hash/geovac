"""GY stage 3 (copy of ga3_3_grid.py with gy_common.phi_ts): the Gaussian kite Phi~(w) on the QG exp-sinh w-grid at
gamma = 1, 32 digits, writing data/gy_grid_g<g>_d<digits>_h<h>.json.  Y and theta^j Y come from
ring_ext_QG_3_values.Y_theta on this grid (the kernel part is shape-independent).
usage: python gy_3_grid.py <gamma> <digits> <htau>"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
from flint import arb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_3_values as V
import gy_common as GY

DATA = Path(__file__).parent / "data"


def phi_adaptive(args):
    w, g, digits = args
    ctx.dps = digits + 12
    w = arb(w)
    h = arb(1) / 16
    prev = GY.phi_ts(w, g, h, digits + 10, "full").real
    while True:
        h = h / 2
        cur = GY.phi_ts(w, g, h, digits + 10, "full").real
        if abs(cur - prev) < abs(cur) * arb(10) ** (-digits) or h < arb(1) / 1024:
            return cur.mid().str(digits + 8, radius=False), str(h)
        prev = cur


def main():
    g, digits, htau = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    fn = DATA / f"gy_grid_g{g}_d{digits}_h{htau}.json"
    grid = V.wgrid(mp.mpf(htau), digits)
    keep = [(w, dw) for w, dw in grid if w ** 3 * dw > mp.mpf(10) ** (-digits - 4) and dw / w ** 5 > mp.mpf(10) ** (-digits - 4)]
    tasks = [(mp.nstr(w, digits + 8), g, digits) for w, dw in keep]
    t0 = time.time()
    with Pool(15) as pool:
        res = pool.map(phi_adaptive, tasks, chunksize=1)
    rec = {"kind": "gauss", "gamma": g, "digits": digits, "htau": htau,
           "w": [mp.nstr(w, digits + 8) for w, _ in keep], "dw": [mp.nstr(dw, digits + 8) for _, dw in keep],
           "phi": [r[0] for r in res], "h_used": [r[1] for r in res], "seconds": time.time() - t0}
    fn.write_text(json.dumps(rec, indent=0))
    print("nodes", len(keep), "seconds", round(rec["seconds"], 1), flush=True)
    for t in ("0.8", "1.3", "1.7", "2.5"):
        print("t", t, "Y_G =", mp.nstr(V.Y_theta(rec, t, 0)[0], digits), flush=True)


if __name__ == "__main__":
    main()
