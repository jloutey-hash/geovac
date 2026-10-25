"""MG1-MG3 stage 3 (copy of cf_g2_3_grid.py): the kite Phi~(w) on an exp-sinh w-grid at gamma = 2 (copy of QG step 3's build_grid,
writing data/mg_grid_g<g>_d<digits>_h<h>.json; the QG script is not modified).  Y and theta^j Y then come
from ring_ext_QG_3_values.Y_theta on this grid.
usage: python cf_g2_3_grid.py <gamma> <digits> <htau>"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_3_values as V

DATA = Path(__file__).parent / "data"


def main():
    g, digits, htau = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    fn = DATA / f"mg_grid_g{g}_d{digits}_h{htau}.json"
    grid = V.wgrid(mp.mpf(htau), digits)
    keep = [(w, dw) for w, dw in grid if w ** 3 * dw > mp.mpf(10) ** (-digits - 4) and dw / w ** 5 > mp.mpf(10) ** (-digits - 4)]
    tasks = [("gem", mp.nstr(w, digits + 8), g, digits) for w, dw in keep]
    t0 = time.time()
    with Pool(15) as pool:
        res = pool.map(V.phi_adaptive, tasks, chunksize=1)
    rec = {"kind": "gem", "gamma": g, "digits": digits, "htau": htau,
           "w": [mp.nstr(w, digits + 8) for w, _ in keep], "dw": [mp.nstr(dw, digits + 8) for _, dw in keep],
           "phi": [r[0] for r in res], "h_used": [r[1] for r in res], "seconds": time.time() - t0}
    fn.write_text(json.dumps(rec, indent=0))
    print("nodes", len(keep), "seconds", round(rec["seconds"], 1), flush=True)
    for t in ("0.8", "1.3", "1.7"):
        print("t", t, "Y =", mp.nstr(V.Y_theta(rec, t, 0)[0], digits), flush=True)


if __name__ == "__main__":
    main()
