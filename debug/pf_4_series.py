"""PF1-PF4 stage 4: Taylor coefficients phi_n of the off-slice kite phi(z) at z = 0 by Cauchy's formula on |z| = rho
(copy of debug/ga3_4_series.py with pf_common.phi_z; no existing script is modified).

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).

Radius (record: 9).  On the principal sheet M = sqrt(u a^2 + (1-u) b^2 + u(1-u) z) vanishes first at
z = -(a+b)^2 = -9 (u = 2/3); N = sqrt((1-v) gamma_1^2 + v gamma_2^2 + v(1-v) z) at -(gamma_1+gamma_2)^2 = -81/4
(v = 1/3); D = s^2 + (1-u-v)^2 z at the two mixed thresholds -(a+gamma_2)^2 = -16 and -(b+gamma_1)^2 = -49/4.
So phi is analytic for |z| < 9 and rho < 9 is admissible; checked by a second radius (the coefficient ratios and the
agreement of the two runs).  The coefficients are NOT in pi Q off the slice (pf_1_checks (e)), so the exact mod-p
route of ga3_5d does not apply; the rows are floats at dps digits.
usage: python pf_4_series.py <Mp> <dps> <nu> <nv> <rho> [tag] [a b g1 g2]"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import pf_common as P

DATA = Path(__file__).parent / "data"


def job(args):
    j, Mp, rho, dps, nu, nv, cfg = args
    ctx.dps = dps
    ang = 2 * arb.pi() * j / Mp
    z = acb(arb(rho) * ang.cos(), arb(rho) * ang.sin())
    v = P.phi_z(z, nu, nv, *cfg)
    return j, v.real.mid().str(dps, radius=False), v.imag.mid().str(dps, radius=False)


def main():
    Mp, dps, nu, nv, rho_s = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    tag = sys.argv[6] if len(sys.argv) > 6 else "main"
    cfg = tuple(sys.argv[7:11]) if len(sys.argv) > 10 else (P.A_SC, P.B_SC, P.G1, P.G2)
    ctx.dps = dps
    rho = arb(rho_s)
    checks = {}
    for name, z in (("-rho", acb(-rho)), ("rho e^{i pi/3}", acb(rho * arb(1) / 2, rho * (arb(3).sqrt() / 2)))):
        a = P.phi_z(z, nu, nv, *cfg)
        b = P.phi_z(z, nu + 200, nv + 200, *cfg)
        checks[name] = float(abs(a - b) / abs(b))
        print(f"node self-check ({nu}x{nv} vs +200) at z = {name}: {checks[name]:.2e}", flush=True)
    t0 = time.time()
    with Pool(15) as pool:
        res = pool.map(job, [(j, Mp, rho_s, dps, nu, nv, cfg) for j in range(Mp // 2 + 1)], chunksize=2)
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
    rec = {"config": {"a": cfg[0], "b": cfg[1], "gamma_1": cfg[2], "gamma_2": cfg[3]}, "Mp": Mp, "rho": rho_s, "dps": dps,
           "nodes": [nu, nv], "node_selfcheck": checks, "seconds": time.time() - t0,
           "phi_n": [mp.nstr(c, dps - 5) for c in coef]}
    (DATA / f"pf_series_{tag}_M{Mp}_d{dps}_rho{rho_s}.json").write_text(json.dumps(rec, indent=0))
    for n in (0, 1, 2, 10, 50, 100, 200, 400, 800, 1200, 1500):
        if n + 1 < len(coef):
            print(n, mp.nstr(coef[n], 25), " ratio c_n/c_{n+1}", mp.nstr(coef[n] / coef[n + 1], 12), flush=True)


if __name__ == "__main__":
    main()
