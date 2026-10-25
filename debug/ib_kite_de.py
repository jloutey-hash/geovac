"""IB record: the kite's differential equation from IBP, mod P, and its operator for phi.

usage: python ib_kite_de.py <tag> <a> <b> <g1> <g2> [R S nsamp sym gem_index]
Per sample z: seed the reducer, close the master basis under d/dz, return c (target in masters) and
A (dM/dz in masters).  Then: rebuild every entry as a rational function of z mod P (checked on held-out
points), build the cyclic vector of the target, and write the operator in d/dz form (mod P) to
debug/data/ib_<tag>_op.json."""
from __future__ import annotations

import json
import random
import sys
import time
from fractions import Fraction as F
from multiprocessing import Pool
from pathlib import Path

from flint import nmod_poly, nmod_mat, nmod

from ib_reduce import P, Reducer, shifts_box, dz_terms, fr
from ib_kite import kite_family, TARGET

DATA = Path(__file__).parent / "data"


def sample(args):
    zi, masses, R_, S_, sym, gidx = args
    fam = kite_family(*masses, gem_index=gidx)
    fam.setup({("w", "w"): zi}, 3)
    seeds = [n for n in shifts_box(5, R_, S_, fam.integer) if not fam.is_zero(n)]
    Rd = Reducer(fam)
    Rd.seed(seeds, sym=(lambda n: (n[1], n[0], n[3], n[2], n[4])) if sym else None)
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
    return zi, ms, {str(m): c.get(m, 0) for m in ms}, {str(m): {str(n): v for n, v in rows[m].items()} for m in ms}


# ---------- rational functions mod P ----------
_rng_norm = random.Random(7)
_NORM = [_rng_norm.randrange(1, P) for _ in range(400)]
def fit_rational(xs, ys, max_deg):
    """smallest total degree (dn, dd) rational function N/D (D monic-normalised by D(0) or lead) through
    the points; returns (N, D) as nmod_poly or None."""
    npts = len(xs)
    for tot in range(0, max_deg + 1):
        for dd in range(0, tot + 1):
            dn = tot - dd
            nunk = dn + 1 + dd + 1
            if nunk > npts - 2:
                continue
            # N(x) - y D(x) = 0, normalise D's constant... use one random linear normalisation row
            rows = []
            for x, y in zip(xs[: nunk - 1], ys[: nunk - 1]):
                rows.append([pow(x, k, P) for k in range(dn + 1)] + [(-y * pow(x, k, P)) % P for k in range(dd + 1)])
            rows.append([0] * (dn + 1) + [_NORM[k] for k in range(dd + 1)])  # random normalisation of D
            rhs = [0] * (nunk - 1) + [1]
            M = nmod_mat(rows, P)
            try:
                sol = M.solve(nmod_mat([[v] for v in rhs], P))
            except ZeroDivisionError:
                continue
            co = [int(sol[k, 0]) for k in range(nunk)]
            N = nmod_poly(co[: dn + 1], P)
            D = nmod_poly(co[dn + 1:], P)
            ok = all((N(x) - nmod(y, P) * D(x)) == 0 for x, y in zip(xs, ys))
            if ok:
                return N, D
    return None


class RF:
    """rational function mod P as (num, den) nmod_poly."""
    __slots__ = ("n", "d")

    def __init__(self, n, d=None):
        self.n = n if isinstance(n, nmod_poly) else nmod_poly([n], P)
        self.d = d if d is not None else nmod_poly([1], P)
        self._norm()

    def _norm(self):
        if self.n.is_zero():
            self.d = nmod_poly([1], P)
            return
        g = self.n.gcd(self.d)
        if g.degree() > 0:
            self.n, self.d = self.n // g, self.d // g
        lc = self.d.coeffs()[-1]
        if int(lc) != 1:
            il = nmod(1, P) / lc
            self.n, self.d = self.n * il, self.d * il

    def __add__(s, o):
        return RF(s.n * o.d + o.n * s.d, s.d * o.d)

    def __sub__(s, o):
        return RF(s.n * o.d - o.n * s.d, s.d * o.d)

    def __mul__(s, o):
        return RF(s.n * o.n, s.d * o.d)

    def __truediv__(s, o):
        return RF(s.n * o.d, s.d * o.n)

    def diff(s):
        return RF(s.n.derivative() * s.d - s.n * s.d.derivative(), s.d * s.d)

    def is_zero(s):
        return s.n.is_zero()


def cyclic_operator(c, A):
    """c: list RF (target in masters); A: matrix RF with dM_i/dz = sum_j A[i][j] M_j.
    v0 = c, v_{k+1} = v_k' + v_k A.  Returns coefficients a_0..a_r (RF) with sum a_k D^k phi = 0, a_r = 1."""
    n = len(c)
    vs = [c]
    while True:
        v = vs[-1]
        nv = [v[j].diff() for j in range(n)]
        for i in range(n):
            if v[i].is_zero():
                continue
            for j in range(n):
                if not A[i][j].is_zero():
                    nv[j] = nv[j] + v[i] * A[i][j]
        vs.append(nv)
        # is vs[-1] in span of vs[:-1]?  solve over GF(P)(z) by evaluating at random points then rebuild:
        r = len(vs) - 1
        rel = _relation(vs)
        if rel is not None:
            return rel
        if r > n:
            raise RuntimeError("no relation up to the master count")


def _relation(vs):
    """find a_k (RF) with sum_{k<r} a_k vs[k] = -vs[r], by Gaussian elimination over RF."""
    r = len(vs) - 1
    n = len(vs[0])
    # matrix columns = vs[0..r-1], rhs = -vs[r]; rows = master components
    M = [[vs[k][i] for k in range(r)] + [RF(nmod_poly([0], P)) - vs[r][i]] for i in range(n)]
    piv_cols = []
    row = 0
    for col in range(r):
        p = next((i for i in range(row, n) if not M[i][col].is_zero()), None)
        if p is None:
            return None  # dependent columns: the earlier vectors are already dependent (should not happen)
        M[row], M[p] = M[p], M[row]
        pv = M[row][col]
        M[row] = [x / pv for x in M[row]]
        for i in range(n):
            if i != row and not M[i][col].is_zero():
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[row])]
        piv_cols.append(col)
        row += 1
    # consistency: remaining rows must have zero rhs
    for i in range(row, n):
        if not M[i][r].is_zero():
            return None
    a = [M[i][r] for i in range(r)] + [RF(nmod_poly([1], P))]
    return a


def main():
    tag = sys.argv[1]
    masses = [F(x) for x in sys.argv[2:6]]
    R_ = int(sys.argv[6]) if len(sys.argv) > 6 else 5
    S_ = int(sys.argv[7]) if len(sys.argv) > 7 else 1
    nsamp = int(sys.argv[8]) if len(sys.argv) > 8 else 80
    sym = (sys.argv[9] == "sym") if len(sys.argv) > 9 else (masses[0] == masses[1] and masses[2] == masses[3])
    gidx = F(sys.argv[10]) if len(sys.argv) > 10 else F(1)
    rng = random.Random(12345)
    zs = [rng.randrange(10**6, P - 1) for _ in range(nsamp)]
    t0 = time.time()
    with Pool(6, maxtasksperchild=1) as pool:
        res = pool.map(sample, [(z, masses, R_, S_, sym, gidx) for z in zs], chunksize=1)
    print(f"{nsamp} samples in {time.time()-t0:.0f}s", flush=True)
    bases = {tuple(r[1]) for r in res}
    if len(bases) != 1:
        counts = {}
        for r in res:
            counts[tuple(r[1])] = counts.get(tuple(r[1]), 0) + 1
        best = max(counts, key=counts.get)
        res = [r for r in res if tuple(r[1]) == best]
        print("WARNING: master basis differs at", nsamp - len(res), "points; using", len(res))
    ms = [str(m) for m in res[0][1]]
    print("masters", len(ms), ms, flush=True)
    hold = 10
    fitp, testp = res[:-hold], res[-hold:]
    xs = [r[0] for r in fitp]

    def rebuild(getter):
        ys = [getter(r) for r in fitp]
        if all(y == 0 for y in ys):
            return RF(nmod_poly([0], P)), 0
        f = fit_rational(xs, ys, len(xs) - 3)
        if f is None:
            raise RuntimeError("rational fit failed: too few samples")
        N, D = f
        for r in testp:
            assert (N(r[0]) - nmod(getter(r), P) * D(r[0])) == 0, "held-out point fails"
        return RF(N, D), N.degree() + D.degree()

    c = []
    maxdeg = 0
    for m in ms:
        f, d = rebuild(lambda r, m=m: r[2][m])
        c.append(f)
        maxdeg = max(maxdeg, d)
    A = []
    for m in ms:
        row = []
        for n in ms:
            f, d = rebuild(lambda r, m=m, n=n: r[3][m].get(n, 0))
            row.append(f)
            maxdeg = max(maxdeg, d)
        A.append(row)
    print("rebuilt; max total degree", maxdeg, "held-out points pass", flush=True)
    a = cyclic_operator(c, A)
    order = len(a) - 1
    den = nmod_poly([1], P)
    for x in a:
        den = den * x.d // den.gcd(x.d)
    polys = [[int(v) for v in (x.n * (den // x.d)).coeffs()] for x in a]
    out = {"tag": tag, "masses": [str(m) for m in masses], "R": R_, "S": S_, "nsamp": nsamp, "sym": sym,
           "gem_index": str(gidx), "P": P, "masters": ms, "max_entry_degree": maxdeg, "order": order,
           "form": "sum_k polys[k](z) (d/dz)^k phi = 0, mod P, cleared denominators",
           "polys": polys, "seconds": time.time() - t0,
           "c": [[[int(v) for v in x.n.coeffs()], [int(v) for v in x.d.coeffs()]] for x in c],
           "A": [[[[int(v) for v in x.n.coeffs()], [int(v) for v in x.d.coeffs()]] for x in row] for row in A]}
    (DATA / f"ib_{tag}_op.json").write_text(json.dumps(out))
    print("order", order, "degrees", [len(p) - 1 for p in polys], f"total {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
