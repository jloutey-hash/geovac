"""GA1 stage 5: the minimal annihilator of Y_log from the kite series (copy of debug/cf_g2_5_guess.py, itself a
copy of the QG step 5b search and the 5d rational scan; file names ga1_*; no existing script is modified).

Y_log = t^3 sum_n binom(-3/2, n) phi_n t^(2n)  (QG memo; any annihilator of Y annihilates Y_log).
Operator ansatz sum_{i<=d} t^(2i) P_i(theta), deg P_i <= r.

usage:
  python ga1_5_guess.py search   <series json> <N> <dps> <tol> <r list> <out tag>   -> data/ga1_5_search_<tag>.json
  python ga1_5_guess.py rational <series json> <N> <dps> <r> <d> <digits> <maxden-exp> <out tag>
                                                                                   -> data/ga1_5_op_<out tag>.json"""
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


def series(path, N, dps):
    ctx.dps = dps
    rec = json.loads(Path(path).read_text())
    phi = [arb(s) for s in rec["phi_n"]][:N]
    c, b = [], arb(1)
    for n in range(N):
        c.append(b * phi[n])
        b = b * (arb(-3) / 2 - n) / (n + 1)
    return c, 3, 2


def ok(c, lam, step, r, d, N, tol):
    out = Gs.try_operator(c, lam, step, r, d, n_to=N, scale=4, tol_digits=tol)
    if out is None:
        return False, None
    p, worst = out
    return (worst < arb(10) ** (-tol)), out


def job(args):
    path, N, dps, tol, r = args
    c, lam, step = series(path, N, dps)
    dmax = (2 * N) // (3 * (r + 1)) - 1          # at least N/3 check rows
    t0 = time.time()
    good, out = ok(c, lam, step, r, dmax, N, tol)
    if not good:
        return r, None, f"none up to d={dmax}", time.time() - t0
    lo, hi = -1, dmax
    best = out
    while hi - lo > 1:
        mid = (lo + hi) // 2
        g_, o_ = ok(c, lam, step, r, mid, N, tol)
        if g_:
            hi, best = mid, o_
        else:
            lo = mid
    return r, hi, {"resid": best[1].mid().str(3, radius=False)}, time.time() - t0


def search(path, N, dps, tol, rs, tag):
    with Pool(min(15, len(rs))) as pool:
        res = pool.map(job, [(path, N, dps, tol, r) for r in rs], chunksize=1)
    out = {"series": str(path), "N": N, "dps": dps, "tol": tol}
    for r, d, info, sec in res:
        print(f"r={r}: d={d}  {info}  ({sec:.0f}s)", flush=True)
        out[r] = {"d": d, "info": info}
    (DATA / f"ga1_5_search_{tag}.json").write_text(json.dumps(out, indent=1))


def rational(path, N, dps, r, d, digits, mde, otag):
    c, lam, step = series(path, N, dps)
    p, worst = Gs.try_operator(c, lam, step, r, d, n_to=N, scale=4, tol_digits=digits)
    print("residual", worst.mid().str(3, radius=False), flush=True)
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
        print("leading theta^r polynomial in t^2:", lead, flush=True)
        print("Q_0 =", sp.factor(Q[0]), "\nQ_d =", sp.factor(Q[d]), flush=True)
        (DATA / f"ga1_5_op_{otag}.json").write_text(json.dumps(
            {"series": str(path), "N": N, "dps": dps, "r": r, "d": d, "Q": [str(q) for q in Q],
             "exact_recheck": ex.mid().str(3, radius=False), "lead": str(lead),
             "note": "operator sum_i t^(2i) Q_i(theta), u = theta = t d/dt, annihilating Y_log"}, indent=1))
        return
    print("no normalisation worked", flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "search":
        search(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), [int(v) for v in sys.argv[6].split(",")], sys.argv[7])
    else:
        rational(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]),
                 int(sys.argv[8]), sys.argv[9])
