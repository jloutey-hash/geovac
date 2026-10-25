"""PF1-PF4 stage 4b: Taylor coefficients phi_n of the off-slice kite phi(z) at z = 0 by truncated power-series arithmetic
in z at every (u, v) Feynman node (flint arb_series, cap N), instead of Cauchy's formula on a circle.

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).

Why: off the slice the Cauchy route loses (R/rho)^n digits at coefficient n (R = 9) and its complex-z kernel costs
~1 ms per node pair at 700 digits, while the coefficient integrands [z^n](M^-3 K) are real and analytic in (u, v) on
the square (branch points of M = sqrt(4 - 3u + u(1-u) z) at u = 4/3 and of N at v = -1/3 only), so every coefficient
comes out at the full working precision; the price is that [z^n] has polynomial degree ~2n in u and v, so the node
counts grow with N (measured by the convergence test, mode `test`, against the Cauchy coefficients of pf_4_series at
small n).  Gauss-Chebyshev (2nd kind, the sqrt(u(1-u)) weight) in u, Gauss-Legendre in v, ring_ext_QG_common's tables.
Per pair: s = M + N, D = s^2 + (1-u-v)^2 z, K = 3/D + 6 s^2/D^2 + 2(N^2+3MN)(3s^2-k2)/D^3 + 24 M N^2 s (s^2-k2)/D^4,
phi = 2 gamma_1 gamma_2 a^3 b^3 sum_u sum_v w_u w_v v(1-v) M^-3 N^-5 K.
usage: python pf_4b_series_direct.py run  <N> <dps> <nu> <nv> <tag> [nproc=15] [sub|plain]
       python pf_4b_series_direct.py test <N> <dps> <nu list> <nv list> <reference series json> [nproc] [sub|plain]
The optional last argument 'sub' selects the substituted quadrature (u_nodes, v_tables_sub)."""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import mpmath as mp
from flint import arb, arb_series, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_common as Cm
import pf_common as P

DATA = Path(__file__).parent / "data"
_VT = {}


QUAD = "plain"      # "plain": Gauss-Chebyshev in u, Gauss-Legendre in v;  "sub": M0 = mu, N0 = sigma substitution


def v_tables(N, dps, nv):
    if QUAD == "sub":
        return v_tables_sub(N, dps, nv)
    key = (N, dps, nv)
    if key in _VT:
        return _VT[key]
    ctx.dps = dps
    ctx.cap = N
    a, b, g1, g2 = P.params()
    g12, g22 = g1 * g1, g2 * g2
    tabs = []
    for v, Wv in Cm.gl_nodes01(nv):
        V = v * (1 - v)
        Ns = arb_series([g12 * (1 - v) + g22 * v, V]).sqrt()          # N(z)
        iN = 1 / Ns
        iN2 = iN * iN
        N2 = Ns * Ns
        cv = Wv * V
        tabs.append((v, cv * (iN2 * iN2 * iN), Ns, N2))           # weight * N^-5, N, N^2
    _VT[key] = tabs
    return tabs


def v_tables_sub(N, dps, nv):
    """sigma = N0 in [gamma_1, gamma_2], v = (sigma^2 - gamma_1^2)/(gamma_2^2 - gamma_1^2), dv = 2 sigma dsigma/(gamma_2^2 - gamma_1^2);
    Gauss-Legendre in y, sigma = gamma_1 + (gamma_2 - gamma_1) y.  The integrand is then rational in sigma (poles at sigma = 0
    and sigma = -mu only), instead of carrying the branch point of N0 = sqrt(...) at v = -1/3."""
    key = ("sub", N, dps, nv)
    if key in _VT:
        return _VT[key]
    ctx.dps = dps
    ctx.cap = N
    a, b, g1, g2 = P.params()
    g12, g22 = g1 * g1, g2 * g2
    tabs = []
    for y, Wy in Cm.gl_nodes01(nv):
        sg = g1 + (g2 - g1) * y
        v = (sg * sg - g12) / (g22 - g12)
        V = v * (1 - v)
        Ns = arb_series([sg * sg, V]).sqrt()          # N(z), N(0) = sigma
        iN = 1 / Ns
        iN2 = iN * iN
        cv = Wy * (g2 - g1) * 2 * sg / (g22 - g12) * V
        tabs.append((v, cv * (iN2 * iN2 * iN), Ns, Ns * Ns))
    _VT[key] = tabs
    return tabs


def u_nodes(nu):
    """[(u, weight for int_0^1 sqrt(u(1-u)) F(u) du)]; 'sub': mu = M0 in [a, b], u = (b^2 - mu^2)/(b^2 - a^2),
    sqrt(u(1-u)) du = sqrt((b-mu)(mu-a)) sqrt((b+mu)(mu+a)) 2 mu dmu/(b^2-a^2)^2, Gauss-Chebyshev (2nd kind) in y,
    mu = a + (b-a) y."""
    if QUAD != "sub":
        return Cm.cheb2_nodes(nu)
    a, b, g1, g2 = P.params()
    out = []
    for y, Wy in Cm.cheb2_nodes(nu):
        mu = a + (b - a) * y
        u = (b * b - mu * mu) / (b * b - a * a)
        w = Wy * (b - a) ** 2 * ((b + mu) * (mu + a)).sqrt() * 2 * mu / (b * b - a * a) ** 2
        out.append((u, w))
    return out


def u_chunk(args):
    N, dps, nv, nu, uidx, quad = args
    global QUAD
    QUAD = quad
    ctx.dps = dps
    ctx.cap = N
    a, b, g1, g2 = P.params()
    a2, b2 = a * a, b * b
    tabs = v_tables(N, dps, nv)
    us = u_nodes(nu)
    total = arb_series([0])
    for k in uidx:
        u, Wu = us[k]
        M = arb_series([a2 * u + b2 * (1 - u), u * (1 - u)]).sqrt()
        iM = 1 / M
        cu = Wu * (iM * iM * iM)
        acc = arb_series([0])
        for v, cv, Ns, N2 in tabs:
            d = 1 - u - v
            W = d * d
            s = M + Ns
            s2 = s * s
            k2 = arb_series([0, W])
            iD = 1 / (s2 + k2)
            MN = M * Ns
            A = (N2 + 3 * MN) * (3 * s2 - k2)
            B = (M * N2) * (s * (s2 - k2))
            # K = iD (3 + iD (6 s^2 + iD (2 A + 24 B iD)))   (Horner in 1/D; same K as ring_ext_QG_common)
            K = iD * (3 + iD * (6 * s2 + iD * (2 * A + 24 * (B * iD))))
            acc += cv * K
        total += cu * acc
    return [c.mid().str(dps + 5, radius=False) for c in total.coeffs()], max(float(c.rad()) for c in total.coeffs()) if total.length() else 0.0


def run(N, dps, nu, nv, nproc):
    ctx.dps = dps
    ctx.cap = N
    a, b, g1, g2 = P.params()
    pref = 2 * g1 * g2 * a * a * a * b * b * b
    chunks = [list(range(nu))[i::nproc] for i in range(nproc)]
    t0 = time.time()
    with Pool(nproc) as pool:
        parts = pool.map(u_chunk, [(N, dps, nv, nu, ch, QUAD) for ch in chunks if ch], chunksize=1)
    mp.mp.dps = dps + 5
    coef = [mp.mpf(0)] * N
    worst_rad = 0.0
    for cs, rad in parts:
        worst_rad = max(worst_rad, rad)
        for n, c in enumerate(cs[:N]):
            coef[n] += mp.mpf(c)
    pf = mp.mpf(pref.mid().str(dps + 5, radius=False))
    coef = [c * pf for c in coef]
    return coef, time.time() - t0, worst_rad


def main():
    global QUAD
    if sys.argv[-1] in ("sub", "plain"):
        QUAD = sys.argv.pop()
    mode = sys.argv[1]
    if mode == "run":
        N, dps, nu, nv, tag = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
        nproc = int(sys.argv[7]) if len(sys.argv) > 7 else 15
        coef, sec, rad = run(N, dps, nu, nv, nproc)
        rec = {"config": {"a": P.A_SC, "b": P.B_SC, "gamma_1": P.G1, "gamma_2": P.G2}, "route": "direct power series in z at every (u,v) node (arb_series)", "quadrature": QUAD,
               "N": N, "dps": dps, "nodes": [nu, nv], "seconds": sec, "worst_ball_radius": rad,
               "phi_n": [mp.nstr(c, dps) for c in coef]}
        (DATA / f"pf_series_direct_{tag}_N{N}_d{dps}_{nu}x{nv}{'_sub' if QUAD == 'sub' else ''}.json").write_text(json.dumps(rec, indent=0))
        print(f"N {N}, dps {dps}, nodes {nu}x{nv}: {sec:.0f}s; worst ball radius {rad:.1e}", flush=True)
        for n in (0, 1, 2, 10, 50, 100, 200, 400, 600, 800, 1000):
            if n + 1 < N:
                print(n, mp.nstr(coef[n], 25), " ratio", mp.nstr(coef[n] / coef[n + 1], 12), flush=True)
    else:
        N, dps = int(sys.argv[2]), int(sys.argv[3])
        nus = [int(v) for v in sys.argv[4].split(",")]
        nvs = [int(v) for v in sys.argv[5].split(",")]
        ref = json.loads(Path(sys.argv[6]).read_text())
        nproc = int(sys.argv[7]) if len(sys.argv) > 7 else 15
        mp.mp.dps = dps + 5
        refc = [mp.mpf(s) for s in ref["phi_n"]]
        out = {}
        for nu, nv in zip(nus, nvs):
            coef, sec, rad = run(N, dps, nu, nv, nproc)
            row = {}
            for n in (0, 10, 25, 50, 75, 100, 125, 150, 199):
                if n < min(N, len(refc)):
                    row[n] = mp.nstr(abs(coef[n] / refc[n] - 1), 3)
            out[f"{nu}x{nv}"] = {"seconds": sec, "rel_vs_reference": row, "worst_ball_radius": rad}
            print(f"nodes {nu}x{nv} ({sec:.0f}s, ball radius {rad:.1e}): rel error vs reference by n:", row, flush=True)
        (DATA / f"pf_series_direct_test_N{N}_d{dps}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
