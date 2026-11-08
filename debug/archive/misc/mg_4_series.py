"""MG1-MG3 (docs/predictions/2026-10-23_gamma1_mechanism.md), stage 4 (copy of cf_g2_4_series.py; gamma, rho as arguments): Taylor coefficients phi_n of the
geminal kite phi(z) at gamma = 2, by Cauchy's formula on |z| = rho.  Copy of debug/ring_ext_QG_4_series.py
with the output renamed to debug/data/mg_series_*.json (the QG script is not modified).

Radius at gamma = 2.  The integrand of Cm.phi_gem_z is analytic in z for |z| < 4: M = sqrt(1 + u(1-u) z)
vanishes first at z = -4 (u = 1/2); N = sqrt(4 + v(1-v) z) only at |z| >= 16; D = s^2 + (u-v)^2 z
stays away from zero (|s| >= Re M + Re N > |u-v| |z|^(1/2) on |z| < 4).  The point z = -(1-gamma)^2 = -1
is therefore not on the principal sheet (as z = -1/4 was not at gamma = 3/2, where rho = 2 worked).
This is checked numerically by a second run at rho = 1 (coefficients must agree).
usage: python cf_g2_4_series.py <gamma> <Mp> <dps> <nodes> <rho>"""
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
    Mp, dps, nn, rho_s = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    rho = arb(rho_s)
    ctx.dps = dps
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
    (DATA / f"mg_series_g{g}_M{Mp}_d{dps}_rho{rho_s}.json").write_text(json.dumps(rec, indent=0))
    for n in (0, 1, 2, 10, 50, 100, 200, 400, 600, 790):
        if n + 1 < len(coef):
            print(n, mp.nstr(coef[n], 25), " ratio", mp.nstr(coef[n] / coef[n + 1], 10), flush=True)


if __name__ == "__main__":
    main()
