"""FP step U: R8 = the minimal operator with R8 P in the t-module of P_h4 (the F_4 piece):
x-side  sum_{i<=I, j<=r} a_ij T^i theta^j g + sum_{i<=I, j<=3} b_ij T^i theta^j h4 = 0 with a != 0.
Smallest r (mod p, the a-part of the nullspace must be nonzero), smallest I, then multi-prime CRT + rational
reconstruction and the exact x-identity.  t-side: R8 = sum_i t^(2i) sum_j a_ij (2 - theta)^j,
C = -sum_i t^(2i) sum_j b_ij (2 - theta)^j, R8 P = C P_h4."""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
import flint

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_path_A_operator as A  # noqa: E402
import four_density_path_C_tmodule as C  # noqa: E402
import four_density_path_D_w12 as D  # noqa: E402
import four_density_path_O_h4 as Oh  # noqa: E402

u = A.u
OUT = HERE / "data" / "four_density_path_U.json"


def build(r, rh, I, p):
    M = A.path_module()
    basis, A0, A1 = C.module_matrices(M, ("c0", "h0"))
    mp_ = C.ModP(A0, A1, p, -I - 2, I + max(r, rh) + 6)
    colsg = C.columns(mp_, r, I)
    h = Oh.vec_of(basis, Oh.H4, mp_)
    th = [h]
    for _ in range(rh):
        th.append(mp_.theta(th[-1]))
    colsh = {}
    for j, E in enumerate(th):
        cur = E
        for i in range(I + 1):
            colsh[(i, j)] = cur
            if i < I:
                cur = mp_.T(cur)
    keys = [("g", i, j) for i in range(I + 1) for j in range(r + 1)] + [("h", i, j) for i in range(I + 1) for j in range(rh + 1)]
    vecs = [colsg[(i, j)] if k == "g" else colsh[(i, j)] for (k, i, j) in keys]
    mats = np.stack([v.reshape(-1) for v in vecs], axis=1)
    mats = mats[np.nonzero(mats.any(axis=1))[0]]
    Mm = flint.nmod_mat(mats.shape[0], mats.shape[1], [int(x) for x in mats.reshape(-1)], p)
    return Mm, keys


def gpart_rank(r, rh, I, p=C.P0):
    Mm, keys = build(r, rh, I, p)
    X, nul = Mm.nullspace()
    gi = [n for n, k in enumerate(keys) if k[0] == "g"]
    if not nul:
        return 0, 0
    rows = [[int(X[i, k]) for k in range(nul)] for i in gi]
    return nul, flint.nmod_mat(len(rows), nul, [x for rr in rows for x in rr], p).rank()


def exact_identity(c, r, rh, I):
    M = A.path_module()
    basis, A0, A1 = C.module_matrices(M, ("c0", "h0"))
    n = len(basis)
    F = lambda x: Fraction(int(sp.Rational(x).p), int(sp.Rational(x).q))
    A0 = [[F(x) for x in row] for row in A0]
    A1 = [[F(x) for x in row] for row in A1]

    def theta(E):
        out = [dict() for _ in range(n)]
        for a in range(n):
            for pw, cf in E[a].items():
                if pw:
                    out[a][pw] = out[a].get(pw, 0) + 2 * pw * cf
                for b in range(n):
                    if A0[b][a]:
                        out[b][pw] = out[b].get(pw, 0) + A0[b][a] * cf
                    if A1[b][a]:
                        out[b][pw + 1] = out[b].get(pw + 1, 0) + A1[b][a] * cf
        return [{k: v for k, v in d.items() if v} for d in out]

    def sub(E, G):
        return [{k: E[a].get(k, 0) - G[a].get(k, 0) for k in set(E[a]) | set(G[a])} for a in range(n)]

    def T(E):
        G = sub(theta(E), E)
        G = sub(theta(G), G)
        return [{k - 1: v for k, v in d.items() if v} for d in G]

    acc = [dict() for _ in range(n)]

    def run(start, kind, rr):
        th = [start]
        for _ in range(rr):
            th.append(theta(th[-1]))
        for j in range(rr + 1):
            cur = th[j]
            for i in range(I + 1):
                cf = c.get((kind, i, j), 0)
                if cf:
                    for a in range(n):
                        for k, v in cur[a].items():
                            acc[a][k] = acc[a].get(k, 0) + cf * v
                if i < I:
                    cur = T(cur)
    g = [dict() for _ in range(n)]
    g[0][0] = Fraction(1)
    run(g, "g", r)
    idx = {m: k for k, m in enumerate(basis)}
    h = [dict() for _ in range(n)]
    for m, cf in Oh.H4.items():
        h[idx[tuple(sorted(m))]][0] = Fraction(cf)
    run(h, "h", rh)
    return all(v == 0 for d in acc for v in d.values())


def main() -> None:
    rec = {}
    t0 = time.time()
    rh = 3
    Imax = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    rmin = None
    for r in range(1, 12):
        nul, gr = gpart_rank(r, rh, Imax)
        print(f"r={r}: nullity {nul}, g-part rank {gr} at t-degree {2*Imax} ({time.time()-t0:.0f}s)", flush=True)
        if gr:
            rmin = r
            break
    lo, hi = 0, Imax
    while lo < hi:
        mid = (lo + hi) // 2
        if gpart_rank(rmin, rh, mid)[1]:
            hi = mid
        else:
            lo = mid + 1
    I = lo
    print(f"minimal: r={rmin}, I={I} (t-degree {2*I}); nullity/g-rank {gpart_rank(rmin, rh, I)}", flush=True)
    rec["r"], rec["I"] = rmin, I
    ps = D.primes_below(400, C.P0 + 1)
    acc, mod, prev = None, 1, None
    for n_, p in enumerate(ps):
        Mm, keys = build(rmin, rh, I, p)
        X, nul = Mm.nullspace()
        assert nul == 1, nul
        v = [int(X[k, 0]) for k in range(len(keys))]
        piv = max(k for k, x in enumerate(v) if x and keys[k][0] == "g")
        iv = pow(v[piv], p - 2, p)
        v = [x * iv % p for x in v]
        if acc is None:
            acc, mod = v, p
        else:
            acc = [(a + mod * ((b - a) * pow(mod, -1, p) % p)) % (mod * p) for a, b in zip(acc, v)]
            mod *= p
        rr = [D.ratrec(a, mod) for a in acc]
        ok = all(x is not None for x in rr)
        if ok and rr == prev:
            print(f"stable after {n_+1} primes ({time.time()-t0:.0f}s)", flush=True)
            break
        prev = rr if ok else None
    c = dict(zip(keys, rr))
    ex = exact_identity(c, rmin, rh, I)
    print("exact x-identity:", ex, f"({time.time()-t0:.0f}s)", flush=True)
    RQ = [sp.Integer(0)] * (I + 1)
    CQ = [sp.Integer(0)] * (I + 1)
    for (k, i, j), v in c.items():
        if v:
            term = sp.Rational(v.numerator, v.denominator) * (2 - u) ** j
            if k == "g":
                RQ[i] += term
            else:
                CQ[i] -= term
    rec["R_Q"] = [str(sp.expand(x)) for x in RQ]
    rec["C_Q"] = [str(sp.expand(x)) for x in CQ]
    rec["exact_x_identity"] = ex
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
