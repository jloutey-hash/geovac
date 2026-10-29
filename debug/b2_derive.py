"""Record B2: from the first-order system of debug/b2_engine.py (no data), the scalar recurrence in L of the output
by the cyclic-vector method, at consecutive points lam0 + s of F_p (shared matrices), with rational reconstruction
of the normalised coefficients.
    python b2_derive.py dims <p>             : system dimension, cyclic order
    python b2_derive.py rec <p> <npts>       : recurrence -> data/b2_recurrence_<p>.json
"""
import json
import random
import sys
import time

import flint

import b2_engine as E
from bu_derive import ratrecon, stack

PRIMES = [(1 << 61) - 1, 4611686018427387847]


class B2Sys:
    def __init__(self, p, lam_ref=987654321):
        self.p = p
        ctx = E.Ctx(p, lam_ref)
        out = E.output(ctx)
        self.coords, _ = E.closure(ctx, list(out))
        self.idx = {m: t for t, m in enumerate(self.coords)}
        self.N = len(self.coords)
        self.cache = {}

    def at(self, lam):
        lam = int(lam) % self.p
        if lam in self.cache:
            return self.cache[lam]
        ctx = E.Ctx(self.p, lam)
        n = self.N
        ent = [0] * (n * n)
        for t, m in enumerate(self.coords):
            for m2, v in E.mono_shift(ctx, m).items():
                ent[t * n + self.idx[m2]] = int(v)        # KeyError = structure not closed at this point
        A = flint.nmod_mat(n, n, ent, self.p)
        cv = [0] * n
        for m2, v in E.output(ctx).items():
            cv[self.idx[m2]] = int(v)
        c = flint.nmod_mat(1, n, cv, self.p)
        if len(self.cache) > 4000:
            self.cache.clear()
        self.cache[lam] = (A, c)
        return A, c

    def rows(self, lam, kmax):
        n = self.N
        Pm = flint.nmod_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)], self.p)
        out = []
        for k in range(kmax + 1):
            A, c = self.at(lam + k)
            out.append(c * Pm)
            if k < kmax:
                Pm = A * Pm
        return out


def cyclic_order(S, lam):
    rows = S.rows(lam, S.N + 1)
    lo, hi = 0, S.N + 1                       # rank of rows[:m+1] <= m first happens at m = order
    M = stack(rows, S.p)
    # incremental: find the first m with rank(rows[:m+1]) <= m
    for m in range(1, S.N + 2):
        if stack(rows[:m + 1], S.p).rank() <= m:
            return m
    return None


def normalized(S, lam, m):
    M = stack(S.rows(lam, m), S.p)
    X, nul = M.transpose().nullspace()
    if nul != 1:
        return None
    v = [X[i, 0] for i in range(m + 1)]
    return [x / v[m] for x in v]


if __name__ == '__main__':
    mode, p = sys.argv[1], int(sys.argv[2])
    t0 = time.time()
    S = B2Sys(p)
    print('system dimension', S.N, round(time.time() - t0, 1), 's', flush=True)
    m = cyclic_order(S, 123456789)
    print('cyclic order', m, round(time.time() - t0, 1), 's', flush=True)
    if mode == 'rec':
        npts = int(sys.argv[3])
        lam0 = random.Random(11).randrange(10 ** 6, p // 2)
        pts, vals = [], []
        for s in range(npts):
            v = normalized(S, lam0 + s, m)
            if v is None:
                continue
            pts.append(lam0 + s)
            vals.append(v)
            if s % 50 == 0:
                print('pt', s, round(time.time() - t0), 's', flush=True)
        fr = [ratrecon([v[k] for v in vals], pts, p) for k in range(m + 1)]
        den = flint.nmod_poly([1], p)
        for _, s_ in fr:
            den = den * s_ // den.gcd(s_)
        polys = [r * (den // s_) for r, s_ in fr]
        g = polys[0]
        for q in polys[1:]:
            g = g.gcd(q)
        polys = [q // g for q in polys]
        degs = [q.degree() for q in polys]
        print('order', m, 'degrees', degs, 'max', max(degs), 'points', len(pts), round(time.time() - t0), 's', flush=True)
        json.dump(dict(p=p, order=m, dim=S.N, npts=len(pts), polys=[[int(c) for c in q.coeffs()] for q in polys]),
                  open(f'data/b2_recurrence_{p}_{npts}.json', 'w'))
