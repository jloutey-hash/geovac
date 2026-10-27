"""PK record, stage 2: the kite's differential system with s = g2^2 free, (a, b, g1) = (1, 2, 3/2).
Run with IB_PRIME set (read at import by ib_reduce).  Pool(6).

  python pko_family.py degree <tag> <nsamp>
      PK5: one fixed random z, nsamp random s; every entry of (c, A) rebuilt as a rational function of s
      (nsamp - 5 points to fit, 5 held out).  Writes debug/data/pko_fam_<tag>.json.
  python pko_family.py grid <tag> <ns> <nz> [s list]
      stage 2b: ns random s (or the given rationals) x nz random z; raw samples to debug/data/pko_fam_<tag>.json.
  python pko_family.py points <tag> <s> <nz>
      direct samples at one rational s and nz fresh random z (PK6 (ii)).
"""
from __future__ import annotations

import json
import random
import sys
import time
from fractions import Fraction as F
from multiprocessing import Pool
from pathlib import Path

from ib_reduce import P, Reducer, shifts_box, dz_terms, fr
from ib_kite import kite_family, TARGET

DATA = Path(__file__).parent / "data"
BASE = (F(1), F(2), F(3, 2))


def sample(args):
    zi, si = args
    fam = kite_family(*BASE, F(3))
    m2 = [m for _, m, _ in fam.props]
    m2[3] = si % P
    fam.setup({("w", "w"): zi}, 3, masses2=m2)
    seeds = [n for n in shifts_box(5, 5, 1, fam.integer) if not fam.is_zero(n)]
    Rd = Reducer(fam)
    Rd.seed(seeds)
    c = Rd.reduce({TARGET: 1})
    masters = set(c)
    rows = {}
    todo = list(masters)
    while todo:
        m = todo.pop()
        r = Rd.reduce(dz_terms(fam, m, "w", zi))
        rows[m] = r
        for n in r:
            if n not in masters:
                masters.add(n)
                todo.append(n)
    ms = sorted(masters, key=Rd.key)
    return zi, si, [str(m) for m in ms], {str(m): c.get(m, 0) for m in ms}, \
        {str(m): {str(n): v for n, v in rows[m].items()} for m in ms}


def run(pairs):
    t0 = time.time()
    with Pool(6, maxtasksperchild=1) as pool:
        res = pool.map(sample, pairs, chunksize=1)
    return res, time.time() - t0


def main():
    mode, tag = sys.argv[1], sys.argv[2]
    out = DATA / f"pko_fam_{tag}.json"
    assert not out.exists(), "refusing to overwrite"
    rng = random.Random(f"{mode}|{tag}|{P}")
    if mode == "degree":
        nsamp = int(sys.argv[3])
        z0 = rng.randrange(10**6, P - 1)
        ss = [rng.randrange(10**6, P - 1) for _ in range(nsamp)]
        res, sec = run([(z0, s) for s in ss])
        from ib_kite_de import fit_rational
        from flint import nmod
        bases = {tuple(r[2]) for r in res}
        assert len(bases) == 1, "master basis differs between samples"
        ms = res[0][2]
        fitp, testp = res[:-5], res[-5:]
        xs = [r[1] for r in fitp]
        degs = []
        fails = 0

        def rebuild(get):
            nonlocal fails
            ys = [get(r) for r in fitp]
            if all(y == 0 for y in ys):
                return None
            f = fit_rational(xs, ys, len(xs) - 3)
            if f is None:
                fails += 1
                return None
            N, D = f
            ok = all((N(r[1]) - nmod(get(r), P) * D(r[1])) == 0 for r in testp)
            if not ok:
                fails += 1
                return None
            return N.degree(), D.degree()

        for m in ms:
            d = rebuild(lambda r, m=m: r[3][m])
            if d:
                degs.append(("c", m, d))
            for n in ms:
                d = rebuild(lambda r, m=m, n=n: r[4][m].get(n, 0))
                if d:
                    degs.append(("A", m + "->" + n, d))
        mx = max(a + b for _, _, (a, b) in degs)
        print(f"{nsamp} samples in {sec:.0f} s wall; masters {len(ms)}; entries {len(degs)}; failed fits {fails}; "
              f"max s-degree num {max(a for _, _, (a, b) in degs)}, den {max(b for _, _, (a, b) in degs)}, total {mx}")
        out.write_text(json.dumps({"P": P, "z0": z0, "s": ss, "seconds": sec, "masters": ms, "failed": fails,
                                   "degrees": degs, "max_total": mx}))
    elif mode in ("grid", "points"):
        if mode == "grid":
            ns, nz = int(sys.argv[3]), int(sys.argv[4])
            ss = [fr(F(x)) for x in sys.argv[5].split(",")] if len(sys.argv) > 5 else \
                [rng.randrange(10**6, P - 1) for _ in range(ns)]
        else:
            ss, nz = [fr(F(sys.argv[3]))], int(sys.argv[4])
        pairs = [(rng.randrange(10**6, P - 1), s) for s in ss for _ in range(nz)]
        res, sec = run(pairs)
        print(f"{len(pairs)} samples in {sec:.0f} s wall")
        out.write_text(json.dumps({"P": P, "seconds": sec, "samples": res}))


if __name__ == "__main__":
    main()
