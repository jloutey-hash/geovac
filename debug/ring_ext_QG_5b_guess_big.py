"""QG step 5b: parallel minimal-degree search (binary search in d for each order r) for annihilators
of phi or of Y_log, from a long high-precision series (ring_ext_QG_4_series_g<g>_M<Mp>.json).
usage: python ring_ext_QG_5b_guess_big.py <gamma> <phi|ylog> <series-tag> <N> <dps> <tol> <r list>"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
sys.set_int_max_str_digits(10**7)
import mpmath as mp
from flint import arb, ctx
sys.path.insert(0, str(Path(__file__).parent))
import os
import ring_ext_QG_guess as Gs
SCALE = float(os.environ.get("QG_SCALE", "4"))

DATA = Path(__file__).parent / "data"


def series(g, which, tag, N, dps):
    ctx.dps = dps
    rec = json.loads((DATA / f"ring_ext_QG_4_series_g{g}{tag}.json").read_text())
    phi = [arb(s) for s in rec["phi_n"]][:N]
    if which == "phi":
        return phi, 0, 1
    c, b = [], arb(1)
    for n in range(N):
        c.append(b * phi[n])
        b = b * (arb(-3) / 2 - n) / (n + 1)
    return c, 3, 2


def ok(c, lam, step, r, d, N, tol):
    out = Gs.try_operator(c, lam, step, r, d, n_to=N, scale=SCALE if SCALE != 4 else 4, tol_digits=tol)
    if out is None:
        return False, None
    p, worst = out
    return (worst < arb(10) ** (-tol)), out


def job(args):
    g, which, tag, N, dps, tol, r = args
    c, lam, step = series(g, which, tag, N, dps)
    dmax = (2 * N) // (3 * (r + 1)) - 1          # at least N/3 check rows
    t0 = time.time()
    good, _ = ok(c, lam, step, r, dmax, N, tol)
    if not good:
        return r, None, f"none up to d={dmax}", time.time() - t0
    lo, hi = -1, dmax           # invariant: hi works, lo fails
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(c, lam, step, r, mid, N, tol)[0]:
            hi = mid
        else:
            lo = mid
    _, out = ok(c, lam, step, r, hi, N, tol)
    p, worst = out
    P = Gs.rationalize(p, maxden=10 ** (tol // 2 - 5), digits=tol)
    ex = Gs.exact_residuals(P, c, lam, step, N) if P else None
    return r, hi, {"resid": worst.mid().str(3, radius=False), "rational": P is not None,
                   "exact": ex.mid().str(3, radius=False) if ex is not None else None,
                   "P": {f"{i},{j}": str(v) for (i, j), v in P.items()} if P else None}, time.time() - t0


def main():
    g, which, tag, N, dps, tol = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    rs = [int(x) for x in sys.argv[7].split(",")]
    tag = "" if tag == "-" else tag
    with Pool(min(15, len(rs))) as pool:
        res = pool.map(job, [(g, which, tag, N, dps, tol, r) for r in rs], chunksize=1)
    out = {}
    for r, d, info, sec in res:
        short = info if isinstance(info, str) else {k: v for k, v in info.items() if k != "P"}
        print(f"r={r}: d={d}  {short}  ({sec:.0f}s)", flush=True)
        out[r] = {"d": d, "info": info}
    (DATA / f"ring_ext_QG_5b_{which}_g{g}{tag}.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
