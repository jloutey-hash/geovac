"""GY stage 4: Taylor coefficients phi_n, phi_a,n, phi_b,n of the Gaussian kite at z = 0 by Cauchy's formula on
|z| = rho (copy of the ga3_4_series.py driver with gy_common.phi_all_z; the integrand is analytic for |z| < 4,
M = sqrt(1 + u(1-u) z) vanishing first at z = -4, and the Gaussian pair carries no other finite point).
Node self-check at z = -rho (nn vs nn + 100 nodes).  Mode 'exp' runs the same driver on the exponential kite of record
QG (ring_ext_QG_common.phi_gem_z) for the positive control against the stored gamma = 3/2 series.
usage: python gy_4_series.py <gamma> <Mp> <dps> <nodes> <rho> [gauss|exp]"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import gy_common as GY
import ring_ext_QG_common as Cm

DATA = Path(__file__).parent / "data"


def job(args):
    g, j, Mp, rho, dps, nn, mode = args
    ctx.dps = dps
    ang = 2 * arb.pi() * j / Mp
    z = acb(arb(rho) * ang.cos(), arb(rho) * ang.sin())
    if 2 * j == Mp:
        z = acb(-arb(rho))          # exact: a ball straddling the negative axis would cross the sqrt cut
    if mode == "exp":
        v = Cm.phi_gem_z(z, arb(g), nn, nn)
        vals = [v]
    else:
        vals = GY.phi_all_z(z, arb(g), nn)
    return j, [(v.real.mid().str(dps, radius=False), v.imag.mid().str(dps, radius=False)) for v in vals]


def main():
    g = sys.argv[1]
    Mp, dps, nn, rho_s = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    mode = sys.argv[6] if len(sys.argv) > 6 else "gauss"
    rho = arb(rho_s)
    ctx.dps = dps
    if mode == "exp":
        a = Cm.phi_gem_z(acb(-rho), arb(g), nn, nn); b = Cm.phi_gem_z(acb(-rho), arb(g), nn + 100, nn + 100)
    else:
        a = GY.phi_z(acb(-rho), arb(g), nn, "full"); b = GY.phi_z(acb(-rho), arb(g), nn + 100, "full")
    selfcheck = float(abs(a - b) / abs(b))
    print("node self-check at z = -rho:", selfcheck, flush=True)
    t0 = time.time()
    with Pool(15) as pool:
        res = pool.map(job, [(g, j, Mp, rho_s, dps, nn, mode) for j in range(Mp // 2 + 1)], chunksize=2)
    print("evaluations", len(res), "seconds", round(time.time() - t0, 1), flush=True)
    mp.mp.dps = dps
    nser = len(res[0][1])
    rec = {"gamma": g, "mode": mode, "Mp": Mp, "rho": rho_s, "dps": dps, "nodes": nn, "node_selfcheck": selfcheck, "seconds": time.time() - t0}
    names = ["phi_n"] if mode == "exp" else ["phi_n", "phia_n", "phib_n"]
    for idx, name in enumerate(names):
        vals = {}
        for j, lst in res:
            re, im = lst[idx]
            vals[j] = mp.mpc(re, im)
            if 0 < j < Mp // 2:
                vals[Mp - j] = mp.conj(vals[j])
        coef = []
        for n in range(Mp // 2):
            s = mp.fsum(vals[j] * mp.expjpi(-2 * mp.mpf(j * n) / Mp) for j in range(Mp))
            coef.append((s / Mp / mp.mpf(rho_s) ** n).real)
        rec[name] = [mp.nstr(c, dps - 5) for c in coef]
        for n in (0, 1, 2, 10, 50, 100, 200, 400, 600, 790):
            if n + 1 < len(coef):
                print(name, n, mp.nstr(coef[n], 25), " ratio", mp.nstr(coef[n] / coef[n + 1], 10), flush=True)
    tag = "exp" if mode == "exp" else "gauss"
    (DATA / f"gy_series_{tag}_g{g}_M{Mp}_d{dps}_rho{rho_s}.json").write_text(json.dumps(rec, indent=0))


if __name__ == "__main__":
    main()
