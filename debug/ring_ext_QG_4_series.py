"""QG step 4: Taylor coefficients phi_n of phi(z) = Phi~(sqrt z) at z = 0 (geminal kite, apex), by
Cauchy's formula on |z| = rho = 2 (nearest singularity |z| = 4: z = -4 and z = -(1+gamma)^2, -4 gamma^2),
Mp points, each phi(z_j) a 2-D Gauss (Chebyshev-2 x Legendre) sum at dps digits; conjugate symmetry
halves the work.  Node self-check: n and n + 40 nodes at the worst point z = -rho.
usage: python ring_ext_QG_4_series.py <gamma> [Mp] [dps] [nodes]"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm

DATA = Path(__file__).parent / "data"


def job(args):
    g, j, Mp, rho, dps, nn = args
    ctx.dps = dps
    ang = 2 * arb.pi() * j / Mp
    z = acb(arb(rho) * ang.cos(), arb(rho) * ang.sin())
    v = Cm.phi_gem_z(z, arb(g), nn, nn)
    return j, v.real.mid().str(dps, radius=False), v.imag.mid().str(dps, radius=False)


def main():
    g = sys.argv[1]
    Mp = int(sys.argv[2]) if len(sys.argv) > 2 else 512
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else 150
    nn = int(sys.argv[4]) if len(sys.argv) > 4 else 220
    rho = arb(sys.argv[5]) if len(sys.argv) > 5 else arb(2)
    rho_s = sys.argv[5] if len(sys.argv) > 5 else "2"
    ctx.dps = dps
    # node self-check at z = -rho (worst case)
    a = Cm.phi_gem_z(acb(-rho), arb(g), nn, nn)
    b = Cm.phi_gem_z(acb(-rho), arb(g), nn + 100, nn + 100)
    selfcheck = float(abs(a - b) / abs(b))
    print("node self-check at z = -rho:", selfcheck, flush=True)
    t0 = time.time()
    with Pool(15) as pool:
        res = pool.map(job, [(g, j, Mp, rho_s, dps, nn) for j in range(Mp // 2 + 1)], chunksize=2)
    print("evaluations", len(res), "seconds", round(time.time() - t0, 1), flush=True)
    mp.mp.dps = dps
    vals = {}
    for j, re, im in res:
        vals[j] = mp.mpc(re, im)
        if 0 < j < Mp // 2:
            vals[Mp - j] = mp.conj(vals[j])
    coef = []
    for n in range(Mp // 2):
        s = mp.fsum(vals[j] * mp.expjpi(-2 * mp.mpf(j * n) / Mp) for j in range(Mp))
        coef.append((s / Mp / mp.mpf(rho_s) ** n).real)
    rec = {"gamma": g, "Mp": Mp, "rho": rho_s, "dps": dps, "nodes": nn, "node_selfcheck": selfcheck,
           "phi_n": [mp.nstr(c, dps - 5) for c in coef]}
    tag = ("" if Mp == 512 else f"_M{Mp}") + ("" if rho_s == "2" else f"_rho{rho_s}")
    (DATA / f"ring_ext_QG_4_series_g{g}{tag}.json").write_text(json.dumps(rec, indent=0))
    for n in (0, 1, 2, 10, 50, 100, 200, 400, 600):
        print(n, mp.nstr(coef[n], 25), " ratio", mp.nstr(coef[n] / coef[n + 1], 10) if n + 1 < len(coef) else "")


if __name__ == "__main__":
    main()
