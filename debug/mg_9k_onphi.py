"""MG3: the kite's operator (Q form in theta_z) applied to phi(z) itself by quadrature, at several real z0 > 0,
with the 1/200 control.  theta_z^k phi(z0) from the Taylor coefficients of phi at z0, obtained by Cauchy's
formula on a small circle |z - z0| = r0 (M points; every point is an independent Feynman-parameter quadrature,
ring_ext_QG_common.phi_gem_z at nu = nv = nodes, self-checked against nodes + 100 at one point).
Residual = |L phi| / max_i |z0^i Q_i(theta) phi| (the QG group convention); the circle's aliasing error is
estimated by the last Taylor coefficient used.
usage: python mg_9k_onphi.py <gamma> <kite op json> <out tag> <z list> [dps=50] [M=64] [r0=0.15] [nodes=300]"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
import sympy as sp
from flint import arb, acb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm

DATA = Path(__file__).parent / "data"
u = sp.Symbol("u")


def job(args):
    g, z0, r0, j, M, dps, nn = args
    ctx.dps = dps
    ang = 2 * arb.pi() * j / M
    zz = acb(arb(z0) + arb(r0) * ang.cos(), arb(r0) * ang.sin())
    v = Cm.phi_gem_z(zz, arb(g), nn, nn)
    return j, v.real.mid().str(dps, radius=False), v.imag.mid().str(dps, radius=False)


def taylor_at(g, z0, r0, M, dps, nn, pool):
    res = pool.map(job, [(g, z0, r0, j, M, dps, nn) for j in range(M)], chunksize=2)
    mp.mp.dps = dps
    vals = {j: mp.mpc(re, im) for j, re, im in res}
    coef = []
    for n in range(M):
        s = mp.fsum(vals[j] * mp.expjpi(-2 * mp.mpf(j * n) / M) for j in range(M))
        coef.append(s / M / mp.mpf(r0) ** n)
    return coef


def main():
    g, opf, tag = sys.argv[1], sys.argv[2], sys.argv[3]
    zs = sys.argv[4].split(",")
    dps = int(sys.argv[5]) if len(sys.argv) > 5 else 50
    M = int(sys.argv[6]) if len(sys.argv) > 6 else 64
    r0 = sys.argv[7] if len(sys.argv) > 7 else "0.15"
    nn = int(sys.argv[8]) if len(sys.argv) > 8 else 300
    rec = json.loads(Path(opf).read_text())
    Q = [sp.Poly(sp.sympify(q, locals={"u": u}), u).all_coeffs()[::-1] for q in rec["Q"]]
    r = max(len(q) for q in Q) - 1
    from sympy.functions.combinatorial.numbers import stirling
    S = [[int(stirling(k, j, kind=2)) for j in range(k + 1)] for k in range(r + 1)]
    # perturbation: the largest |coefficient| in the theta^0 parts
    big = max(((i, abs(int(q[0]))) for i, q in enumerate(Q) if q and q[0] != 0), key=lambda v: v[1])[0]
    ctx.dps = dps
    out = {"gamma": g, "op": opf, "dps": dps, "M": M, "r0": r0, "nodes": nn,
           "perturbed": {"z_power": big, "theta_power": 0, "factor": "1 + 1/200"}}
    za = acb(zs[0])
    a1 = Cm.phi_gem_z(za, arb(g), nn, nn)
    a2 = Cm.phi_gem_z(za, arb(g), nn + 100, nn + 100)
    out["node_selfcheck"] = float(abs(a1 - a2) / abs(a2))
    print("node self-check (nodes vs nodes+100) at z =", zs[0], ":", out["node_selfcheck"], flush=True)
    t0 = time.time()
    with Pool(8) as pool:
        for zs0 in zs:
            coef = taylor_at(g, zs0, r0, M, dps, nn, pool)
            mp.mp.dps = dps
            z0 = mp.mpf(zs0)
            # theta^k phi = sum_j S(k,j) z0^j phi^(j), phi^(j) = j! a_j
            der = [mp.factorial(j) * coef[j] for j in range(r + 1)]
            th = [mp.fsum(S[k][j] * z0**j * der[j] for j in range(k + 1)) for k in range(r + 1)]
            alias = abs(coef[M - 1] * mp.mpf(r0) ** (M - 1)) / abs(coef[0])

            def apply(perturb):
                groups = []
                for i, q in enumerate(Q):
                    acc = mp.mpf(0)
                    for k, c in enumerate(q):
                        c = mp.mpf(int(c))
                        if perturb and i == big and k == 0:
                            c *= 1 + mp.mpf(1) / 200
                        acc += c * z0**i * th[k]
                    groups.append(acc)
                return abs(mp.fsum(groups)) / max(abs(v) for v in groups)
            row = {"phi": mp.nstr(th[0].real, 25), "LY_rel_group": mp.nstr(apply(False).real, 3),
                   "perturbed_rel_group": mp.nstr(apply(True).real, 3), "circle_alias_est": mp.nstr(alias, 3)}
            out[zs0] = row
            print(zs0, row, f"({time.time()-t0:.0f}s)", flush=True)
    (DATA / f"mg_9k_onphi_{tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
