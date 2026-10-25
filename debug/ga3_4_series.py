"""GA3 stage 4 (copy of ga1_4_series.py, output ga3_*; at gamma = 2 the integrand is analytic for |z| < 4 too, N = sqrt(4 + v(1-v) z) vanishes only at |z| >= 16): Taylor coefficients phi_n of the cusp-correct-geminal kite phi(z) at gamma = 3/2, by Cauchy's
formula on |z| = rho.  Copy of debug/cf_g2_4_series.py with ga1_common.phi_z (shape 'lin' by default) and
output debug/data/ga3_series_<shape>_g<g>_M<Mp>_d<dps>_rho<rho>.json (no existing script is modified).

Radius (as at gamma = 3/2; at gamma = 2 likewise): the integrand is analytic in z for |z| < 4 (M = sqrt(1 + u(1-u) z) vanishes first at
z = -4; N = sqrt(9/4 + v(1-v) z) only at |z| >= 9; D = s^2 + (u-v)^2 z stays away from 0 for |z| < 4), so
phi is analytic there, and the record-QG radius rho = 2 applies to the lin shape too; the point
z = -(1 - gamma)^2 = -1/4 is not on the principal sheet.  Checked by a second run at rho = 1.
usage: python ga1_4_series.py <gamma> <Mp> <dps> <nodes> <rho> [shape]"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ga1_common as GA

DATA = Path(__file__).parent / "data"


def job(args):
    g, j, Mp, rho, dps, nn, shape = args
    ctx.dps = dps
    ang = 2 * arb.pi() * j / Mp
    z = acb(arb(rho) * ang.cos(), arb(rho) * ang.sin())
    v = GA.phi_z(z, arb(g), nn, nn, shape)
    return j, v.real.mid().str(dps, radius=False), v.imag.mid().str(dps, radius=False)


def main():
    g = sys.argv[1]
    Mp, dps, nn, rho_s = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    shape = sys.argv[6] if len(sys.argv) > 6 else "lin"
    rho = arb(rho_s)
    ctx.dps = dps
    a = GA.phi_z(acb(-rho), arb(g), nn, nn, shape)
    b = GA.phi_z(acb(-rho), arb(g), nn + 100, nn + 100, shape)
    selfcheck = float(abs(a - b) / abs(b))
    print("node self-check at z = -rho:", selfcheck, flush=True)
    t0 = time.time()
    with Pool(15) as pool:
        res = pool.map(job, [(g, j, Mp, rho_s, dps, nn, shape) for j in range(Mp // 2 + 1)], chunksize=2)
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
    rec = {"gamma": g, "shape": shape, "Mp": Mp, "rho": rho_s, "dps": dps, "nodes": nn, "node_selfcheck": selfcheck,
           "phi_n": [mp.nstr(c, dps - 5) for c in coef]}
    (DATA / f"ga3_series_{shape}_g{g}_M{Mp}_d{dps}_rho{rho_s}.json").write_text(json.dumps(rec, indent=0))
    for n in (0, 1, 2, 10, 50, 100, 200, 400, 600, 790):
        if n + 1 < len(coef):
            print(n, mp.nstr(coef[n], 25), " ratio", mp.nstr(coef[n] / coef[n + 1], 10), flush=True)


if __name__ == "__main__":
    main()
