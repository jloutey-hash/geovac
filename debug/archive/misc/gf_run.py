"""GF1-GF3 (record 2026-10-25_gaussian_fit_geminal.md): Y_n(t) = sum_ij c_i c_j Y_GG(alpha_i g1^2, alpha_j g2^2; t) from the
fits in data/gf_fits.json (one fit at gamma = 1, rescaled), against the exact Y(t) of the 60-digit grid; E_n = Y_n/Y - 1.
Each Y_GG at the default node set; a sample of pairs is repeated at the larger node set (the quadrature control).
usage: python gf_run.py <nproc>"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
import mpmath as mp
sys.path.insert(0, str(Path(__file__).parent))
import gf_engine as F

DATA = Path(__file__).parent / "data"
TS = [0.5, 1.0, 2.0, 4.0]
G1, G2 = 1.5, 3.0


def task(args):
    A, B, big = args
    if big:
        return F.Y_GG(A, B, TS, nw=36, ns=36, nd=36).tolist()
    return F.Y_GG(A, B, TS).tolist()


def main():
    nproc = int(sys.argv[1])
    fits = json.loads((DATA / "gf_fits.json").read_text())
    jobs, keys = [], []
    for wname, byn in fits.items():
        for n, f in byn.items():
            al = [float(v) for v in f["alpha"]]
            for i in range(len(al)):
                for j in range(len(al)):
                    keys.append((wname, n, i, j))
                    jobs.append((al[i] * G1 ** 2, al[j] * G2 ** 2, False))
    # quadrature control: every pair of the n = 10 W1 fit at the larger node set
    ctrl = [(k, (a, b, True)) for k, (a, b, _) in zip(keys, jobs) if k[0] == "W1" and k[1] == max(fits["W1"], key=int)]
    t0 = time.time()
    with Pool(nproc) as pool:
        vals = pool.map(task, jobs + [c[1] for c in ctrl], chunksize=1)
    print(f"{len(jobs)} pairs + {len(ctrl)} control pairs ({time.time()-t0:.0f}s)", flush=True)
    base = dict(zip(keys, vals[:len(jobs)]))
    worst = 0.0
    for (k, _), v in zip(ctrl, vals[len(jobs):]):
        worst = max(worst, float(np.max(np.abs(np.array(base[k]) / np.array(v) - 1))))
    print("quadrature control (n = 10, W1, all pairs): worst relative difference", f"{worst:.2e}", flush=True)

    import ring_ext_QG_8_onY as E
    mp.mp.dps = 70
    g60 = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
    Yex = [float(E.thetaY(g60, str(t), 0)[0]) for t in TS]
    out = {"t": TS, "Y_exact": Yex, "quadrature_control_worst": worst, "E": {}}
    for wname, byn in fits.items():
        out["E"][wname] = {}
        for n, f in sorted(byn.items(), key=lambda kv: int(kv[0])):
            c = [float(v) for v in f["c"]]
            Yn = np.zeros(len(TS))
            for i in range(len(c)):
                for j in range(len(c)):
                    Yn += c[i] * c[j] * np.array(base[(wname, n, i, j)])
            En = Yn / np.array(Yex) - 1
            out["E"][wname][n] = En.tolist()
            print(wname, f"n={n:>2}", "E_n(t) at t = 0.5, 1, 2, 4:", " ".join(f"{e:+.3e}" for e in En), flush=True)
    (DATA / "gf_results.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
