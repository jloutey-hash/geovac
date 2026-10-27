"""NP record (docs/predictions/2026-10-26_kite_ninth_point.md): shared engine.
One reduction sample of the kite at given squared masses (residues mod P, P from IB_PRIME at import), massless or
massive Coulomb line; the z9 readout from the five top-row carriers; rational reconstruction of residues."""
from __future__ import annotations

import time
from fractions import Fraction as F
from multiprocessing import Pool

from ib_reduce import P, Reducer, shifts_box, dz_terms, fr
from ib_kite import kite_family, TARGET

TOP = "(0, 0, 0, 0, 1)"
CARRIERS = ["(0, 0, 0, 0, -1)", "(0, 1, -1, 0, 0)", "(0, 1, 0, 0, -1)", "(1, 0, 0, -1, 0)", "(1, 0, 0, 0, -1)"]
NPROC = 7


def sample(args):
    """args = (zi, m2 (5 residues), massive, keep) -> dict.  keep: 'all' or 'top' (top row + c only)."""
    zi, m2, massive, keep = args
    t0 = time.process_time()
    fam = kite_family()
    if massive:
        fam.zero_sector = lambda pos: not pos[2] and not pos[3] and not pos[4]
    fam.setup({("w", "w"): zi}, 3, masses2=[int(x) % P for x in m2])
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
    A = {str(m): {str(n): v for n, v in rows[m].items()} for m in ms}
    if keep == "top":
        A = {TOP: A.get(TOP, {})}
    return {"z": zi, "m2": [int(x) for x in m2], "masters": [str(m) for m in ms],
            "c": {str(m): c.get(m, 0) for m in ms}, "A": A, "cpu": time.process_time() - t0}


def run(arglist, nproc=NPROC):
    t0 = time.time()
    with Pool(nproc, maxtasksperchild=1) as pool:
        res = pool.map(sample, arglist, chunksize=1)
    return res, time.time() - t0


def readout(samps):
    """samps: >= 4 samples at one rate point.  Returns (r0, ok, note): r0 the common pole of the five carriers
    under y = C/(z - r0), fitted from samples 0,1 and checked on all others."""
    out = None
    for lab in CARRIERS:
        ys = [s["A"][TOP].get(lab, 0) % P for s in samps]
        zs = [s["z"] for s in samps]
        y1, y2 = ys[0], ys[1]
        if (y1 - y2) % P == 0:
            return None, False, f"carrier {lab}: degenerate"
        r = (y1 * zs[0] - y2 * zs[1]) * pow((y1 - y2) % P, P - 2, P) % P
        C = y1 * (zs[0] - r) % P
        if C == 0:
            return None, False, f"carrier {lab}: zero"
        for y, z in zip(ys, zs):
            if y * (z - r) % P != C:
                return None, False, f"carrier {lab}: check fails"
        if out is None:
            out = r
        elif out != r:
            return None, False, "carriers disagree"
    return out, True, "ok"


def ratrec(u, p=None):
    """rational reconstruction of residue u mod p (|num|, den <= sqrt(p/2))."""
    p = p or P
    import math
    bound = math.isqrt(p // 2)
    r0, r1, s0, s1 = p, u % p, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    x = F(r1, s1)
    return x if (x.numerator % p * pow(x.denominator % p, p - 2, p) - u) % p == 0 else None
