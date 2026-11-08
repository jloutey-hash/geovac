"""PG1-PG3 (docs/predictions/2026-10-26_production_geminal_fit.md, committed 0a98d7e3 before computing): Psi4's shipped
six-Gaussian fit of exp(-r) (libmints/integralparameters.cc, FittedSlaterCorrelationFactor; magnitudes of its
coefficients) on the exact off-slice Y, by record GF's engine; record GF's W1 n = 6 row recomputed alongside; one
Psi4 pair at the larger node set (quadrature control).
usage: python pg_run.py <nproc>"""
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
PSI4 = {"alpha": [0.2209, 1.004, 3.622, 12.16, 45.87, 254.4],
        "c": [0.3144, 0.3037, 0.1681, 0.09811, 0.06024, 0.03726]}


def task(args):
    A, B, big = args
    if big:
        return F.Y_GG(A, B, TS, nw=36, ns=36, nd=36).tolist()
    return F.Y_GG(A, B, TS).tolist()


def Yfit(fit, nproc):
    al, c = [float(v) for v in fit["alpha"]], [float(v) for v in fit["c"]]
    jobs = [(al[i] * G1 ** 2, al[j] * G2 ** 2, False) for i in range(len(al)) for j in range(len(al))]
    with Pool(nproc) as pool:
        vals = pool.map(task, jobs, chunksize=1)
    Y = np.zeros(len(TS))
    k = 0
    for i in range(len(al)):
        for j in range(len(al)):
            Y += c[i] * c[j] * np.array(vals[k])
            k += 1
    return Y, jobs, vals


def main():
    nproc = int(sys.argv[1])
    t0 = time.time()
    import ring_ext_QG_8_onY as E
    mp.mp.dps = 70
    g60 = json.loads((DATA / "pf_grid_d60_h0.03125.json").read_text())
    Yex = np.array([float(E.thetaY(g60, str(t), 0)[0]) for t in TS])
    w1 = json.loads((DATA / "gf_fits.json").read_text())["W1"]["6"]
    out = {"t": TS, "Y_exact": Yex.tolist(), "E": {}}
    for name, fit in (("Psi4", PSI4), ("W1_n6_repro", w1)):
        Y, jobs, vals = Yfit(fit, nproc)
        E_ = Y / Yex - 1
        out["E"][name] = E_.tolist()
        print(name, "E(t) at t = 0.5, 1, 2, 4:", " ".join(f"{e:+.4e}" for e in E_), flush=True)
        if name == "Psi4":
            A, B, _ = jobs[0]
            big = np.array(task((A, B, True)))
            ctrl = float(np.max(np.abs(np.array(vals[0]) / big - 1)))
            out["quadrature_control_pair0"] = ctrl
            print("quadrature control, smallest-exponent pair, larger node set: worst relative difference", f"{ctrl:.2e}", flush=True)
    out["seconds"] = time.time() - t0
    (DATA / "pg_results.json").write_text(json.dumps(out, indent=1))
    print(f"done {out['seconds']:.0f}s", flush=True)


if __name__ == "__main__":
    main()
