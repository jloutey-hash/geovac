"""PF1-PF4 stage 5: minimal annihilators from the off-slice kite series (copy of debug/mg_5_guess.py, parametrised;
ring_ext_QG_guess imported, not modified).

Object: the exponential-geminal correlated triangle Y(t) at the fingerprint fibre with vertex 1 = density rho_i of
scale a = 1, vertex 2 = density rho_j of scale b = 2, vertex 3 = apex density rho_k of scale t; geminal
exp(-gamma_1 r_13) with gamma_1 = 3/2, geminal exp(-gamma_2 r_23) with gamma_2 = 3, Coulomb 1/r_12 (QG conventions).

kind ylog: Y_log = t^3 sum_n binom(-3/2, n) phi_n t^(2n), ansatz sum_{i<=d} t^(2i) P_i(theta_t)   (lam 3, step 2)
kind phi : phi(z) = sum_n phi_n z^n, ansatz sum_{i<=d} z^i P_i(theta_z)                            (lam 0, step 1)
scale: rows weighted scale^n (ring_ext_QG_guess.build_rows); scale = the Taylor radius 9 keeps c_n scale^n O(1).
Search: for each order r the largest degree d with at least `chk` check rows (default N/3 as in the lineage; a
smaller chk is allowed for the large operators and is stated); a candidate passes if the worst relative residual
on the check rows is below 10^-tol; then the minimal d by bisection.  The rational step recognises the ratios of
the null vector (continued fractions, `digits` agreement, denominators <= 10^mde), clears denominators and
RECHECKS the integer operator exactly on all N coefficients (ring_ext_QG_guess.exact_residuals).
usage:
  python pf_5_guess.py search   <series json> <kind> <N> <dps> <tol> <scale> <r list> <chk> <tag> [nproc=8]
  python pf_5_guess.py rational <series json> <kind> <N> <dps> <r> <d> <digits> <maxden-exp> <scale> <out tag>"""
from __future__ import annotations
import json, sys, time
from fractions import Fraction
from math import gcd, lcm
from multiprocessing import Pool
from pathlib import Path
sys.set_int_max_str_digits(10**7)
import mpmath as mp
import sympy as sp
from flint import arb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_guess as Gs

DATA = Path(__file__).parent / "data"
u, z = sp.symbols("u z")


def series(path, kind, N, dps):
    ctx.dps = dps
    rec = json.loads(Path(path).read_text())
    phi = [arb(s) for s in rec["phi_n"]][:N]
    if kind == "phi":
        return phi, 0, 1
    c, b = [], arb(1)
    for n in range(N):
        c.append(b * phi[n])
        b = b * (arb(-3) / 2 - n) / (n + 1)
    return c, 3, 2


def ok(c, lam, step, r, d, N, tol, scale):
    out = Gs.try_operator(c, lam, step, r, d, n_to=N, scale=scale, tol_digits=tol)
    if out is None:
        return False, None
    p, worst = out
    return (worst < arb(10) ** (-tol)), out


def job(args):
    path, kind, N, dps, tol, scale, r, chk = args
    c, lam, step = series(path, kind, N, dps)
    dmax = (N - chk) // (r + 1) - 1          # at least chk check rows
    t0 = time.time()
    good, out = ok(c, lam, step, r, dmax, N, tol, scale)
    if not good:
        return r, None, f"none up to d={dmax}" + (f" (resid {out[1].mid().str(3, radius=False)})" if out else ""), time.time() - t0
    lo, hi = -1, dmax
    best = out
    while hi - lo > 1:
        mid = (lo + hi) // 2
        g_, o_ = ok(c, lam, step, r, mid, N, tol, scale)
        if g_:
            hi, best = mid, o_
        else:
            lo = mid
    return r, hi, {"resid": best[1].mid().str(3, radius=False)}, time.time() - t0


def search(path, kind, N, dps, tol, scale, rs, chk, tag, nproc=8):
    with Pool(min(nproc, len(rs))) as pool:
        res = pool.map(job, [(path, kind, N, dps, tol, scale, r, chk) for r in rs], chunksize=1)
    out = {"series": str(path), "kind": kind, "N": N, "dps": dps, "tol": tol, "scale": scale, "min_check_rows": chk}
    for r, d, info, sec in res:
        print(f"r={r}: d={d}  {info}  ({sec:.0f}s)", flush=True)
        out[r] = {"d": d, "info": info}
    (DATA / f"pf_5_search_{tag}.json").write_text(json.dumps(out, indent=1))


def rational(path, kind, N, dps, r, d, digits, mde, scale, otag):
    c, lam, step = series(path, kind, N, dps)
    p, worst = Gs.try_operator(c, lam, step, r, d, n_to=N, scale=scale, tol_digits=digits)
    print("float residual on the check rows", worst.mid().str(3, radius=False), flush=True)
    mp.mp.dps = dps
    keys = sorted(p)
    vals = {k: mp.mpf(p[k].mid().str(dps, radius=False)) for k in keys}
    big = max(abs(v) for v in vals.values())
    nz = [k for k in keys if abs(vals[k]) > mp.mpf(10) ** (-digits) * big]
    tol = mp.mpf(10) ** (-digits)
    for ref in sorted(nz, key=lambda k: -abs(vals[k])):
        fr, okk = {}, True
        for k in keys:
            xv = vals[k] / vals[ref]
            q = Fraction(mp.nstr(xv, digits + 5, min_fixed=-10**9, max_fixed=10**9)).limit_denominator(10 ** mde)
            if abs(mp.mpf(q.numerator) / q.denominator - xv) > tol * max(1, abs(xv)):
                okk = False
                break
            fr[k] = q
        if not okk:
            continue
        L = 1
        for q in fr.values():
            L = lcm(L, q.denominator)
        ints = {k: int(q * L) for k, q in fr.items()}
        G = 0
        for v in ints.values():
            G = gcd(G, v)
        P = {k: v // G for k, v in ints.items()}
        ex = Gs.exact_residuals(P, c, lam, step, N)
        print("normalised on", ref, ": integer operator, max |coeff| ~ 10^%d, exact recheck on all %d coefficients %s"
              % (len(str(max(abs(v) for v in P.values()))), N, ex.mid().str(3, radius=False)), flush=True)
        Q = [sp.expand(sum(sp.Integer(P[(i, j)]) * u**j for j in range(r + 1) if (i, j) in P)) for i in range(d + 1)]
        lead = sp.factor(sum(Q[i].coeff(u, r) * z ** i for i in range(d + 1)))
        print("leading theta^r polynomial in z (phi) or t^2 (ylog):", lead, flush=True)
        print("Q_0 =", sp.factor(Q[0]), "\nQ_d =", sp.factor(Q[d]), flush=True)
        note = ("operator sum_i t^(2i) Q_i(theta), u = theta = t d/dt, annihilating Y_log" if kind == "ylog"
                else "operator sum_i z^i Q_i(theta_z), u = theta_z = z d/dz, annihilating phi(z)")
        (DATA / f"pf_5_op_{otag}.json").write_text(json.dumps(
            {"series": str(path), "kind": kind, "N": N, "dps": dps, "r": r, "d": d, "Q": [str(q) for q in Q],
             "float_residual_check_rows": worst.mid().str(3, radius=False), "recognition_digits": digits, "maxden_exp": mde,
             "exact_recheck": ex.mid().str(3, radius=False), "lead": str(lead), "note": note}, indent=1))
        return
    print("no normalisation worked", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "search":
        search(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]),
               [int(v) for v in sys.argv[8].split(",")], int(sys.argv[9]), sys.argv[10],
               int(sys.argv[11]) if len(sys.argv) > 11 else 8)
    else:
        rational(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]),
                 int(sys.argv[8]), int(sys.argv[9]), int(sys.argv[10]), sys.argv[11])
