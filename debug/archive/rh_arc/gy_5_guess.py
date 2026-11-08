"""GY stage 5 (float route, copy of ga3_5_guess.py made multi-series): minimal annihilator of Y_G,log from the 600-digit
Cauchy series, with the split series stacked: rows from the phi_a-part (restricted to the coefficients above the Cauchy
floor, the a-part being entire with factorially decaying coefficients) and from the phi_b-part; or from phi alone
('phi' mode: the minimal annihilator over C(x) at the degree bound, which may be smaller than the one over Q(x)).
Operator ansatz sum_{i<=d} t^(2i) P_i(theta), deg P_i <= r; U - 1 rows + a random normalisation, residual on the rest.
usage: python gy_5_guess.py search <series json> <N> <dps> <tol> <r list> <out tag> [phi|ab]"""
from __future__ import annotations
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
sys.set_int_max_str_digits(10**7)
import mpmath as mp
import random
from flint import arb, arb_mat, ctx
sys.path.insert(0, str(Path(__file__).parent))
import ring_ext_QG_guess as Gs

DATA = Path(__file__).parent / "data"


def series(path, N, dps, mode):
    ctx.dps = dps
    rec = json.loads(Path(path).read_text())
    out = []
    keys = ["phi_n"] if mode == "phi" else ["phia_n", "phib_n"]
    for key in keys:
        phi = [arb(s) for s in rec[key]][:N]
        nuse = N
        if key == "phia_n":
            rho = mp.mpf(rec["rho"]); mp.mp.dps = 30
            nuse = 0
            for n in range(N):
                if abs(mp.mpf(rec[key][n])) > mp.mpf(10) ** (-dps + 60) / rho ** n:
                    nuse = n + 1
        c, b = [], arb(1)
        for n in range(nuse):
            c.append(b * phi[n])
            b = b * (arb(-3) / 2 - n) / (n + 1)
        out.append(c)
    return out


def try_multi(cs, lam, step, r, d, scale, tol):
    U = (r + 1) * (d + 1)
    rows = []
    for c in cs:
        rows += Gs.build_rows(c, lam, step, r, d, 0, len(c), scale)
    if len(rows) < U + 10:
        return None
    rng = random.Random(1)
    xi = [arb(rng.randint(1, 10**6)) for _ in range(U)]
    # use the first U-1 rows spread over both series: interleave
    A = arb_mat(rows[: U - 1] + [xi])
    b = arb_mat([[0]] * (U - 1) + [[1]])
    try:
        x = A.solve(b, algorithm="approx")
    except Exception:
        return None
    xs = [x[k, 0].mid() for k in range(U)]
    worst = mp.mpf(0)
    for row in rows[U - 1:]:
        res = sum((a.mid() * v for a, v in zip(row, xs)), arb(0))
        size = sum((abs(a.mid() * v) for a, v in zip(row, xs)), arb(0))
        rel = mp.mpf(1) if size == 0 else abs(mp.mpf((res / size).mid().str(20, radius=False)))
        worst = max(worst, rel)
    return worst < mp.mpf(10) ** (-tol), worst


def job(args):
    path, N, dps, tol, r, mode = args
    cs = series(path, N, dps, mode)
    total = sum(len(c) for c in cs)
    dmax = (2 * total) // (3 * (r + 1)) - 1
    t0 = time.time()
    res = try_multi(cs, 3, 2, r, dmax, 4, tol)
    if res is None or not res[0]:
        return r, None, f"none up to d={dmax}", time.time() - t0
    lo, hi = -1, dmax
    best = res
    while hi - lo > 1:
        mid = (lo + hi) // 2
        rr = try_multi(cs, 3, 2, r, mid, 4, tol)
        if rr is not None and rr[0]:
            hi, best = mid, rr
        else:
            lo = mid
    return r, hi, {"resid": mp.nstr(best[1], 3)}, time.time() - t0


if __name__ == "__main__":
    path, N, dps, tol = sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    rs = [int(v) for v in sys.argv[6].split(",")]
    tag = sys.argv[7]
    mode = sys.argv[8] if len(sys.argv) > 8 else "ab"
    with Pool(min(15, len(rs))) as pool:
        res = pool.map(job, [(path, N, dps, tol, r, mode) for r in rs], chunksize=1)
    out = {"series": path, "N": N, "dps": dps, "tol": tol, "mode": mode}
    for r, d, info, sec in res:
        print(f"r={r}: d={d}  {info}  ({sec:.0f}s)", flush=True)
        out[r] = {"d": d, "info": info}
    (DATA / f"gy_5_search_{tag}.json").write_text(json.dumps(out, indent=1))
