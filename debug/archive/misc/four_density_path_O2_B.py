"""FP step O2: the exact operator B (order <= 11) and polynomial q with B P + q(t^2) P_h4 = 0 (step O found it mod p).
x-side: sum b_ij T^i theta^j g + sum q_i T^i h4 = 0; t-side: sum b_ij t^(2i) (2 - theta_t)^j P + sum q_i t^(2i) P_h4 = 0.
Minimal I by bisection (nullity 1), multi-prime CRT + rational reconstruction, exact x-identity check over Q.
Also A4 = Mellin operator of P_h4 (alpha = 4, beta = 3 on K0^3) and its equality with the x-scan order-4 operator."""
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
import density_loop_B_ode as Bm  # noqa: E402
import four_density_path_A_operator as A  # noqa: E402
import four_density_path_C_tmodule as C  # noqa: E402
import four_density_path_D_w12 as D  # noqa: E402
import four_density_path_O_h4 as Oh  # noqa: E402

u, z = A.u, A.z
t = sp.Symbol("t")
OUT = HERE / "data" / "four_density_path_O2.json"


def system(r, I, p):
    M = A.path_module()
    basis, A0, A1 = C.module_matrices(M, ("c0", "h0"))
    mp_ = C.ModP(A0, A1, p, -I - 2, I + r + 6)
    cols = C.columns(mp_, r, I)
    keys = [("g", i, j) for i in range(I + 1) for j in range(r + 1)]
    h = Oh.vec_of(basis, Oh.H4, mp_)
    hc = [h]
    for _ in range(I):
        hc.append(mp_.T(hc[-1]))
    keys += [("h", i, 0) for i in range(I + 1)]
    mats = np.stack([cols[(i, j)].reshape(-1) for (_, i, j) in keys if _ == "g"] + [x.reshape(-1) for x in hc], axis=1)
    mats = mats[np.nonzero(mats.any(axis=1))[0]]
    return flint.nmod_mat(mats.shape[0], mats.shape[1], [int(x) for x in mats.reshape(-1)], p), keys


def kernel(r, I, p):
    Mm, keys = system(r, I, p)
    X, nul = Mm.nullspace()
    return X, nul, keys


def main() -> None:
    rec = {}
    t0 = time.time()
    r = 11
    lo, hi = 21, 40
    while lo < hi:
        mid = (lo + hi) // 2
        if kernel(r, mid, C.P0)[1]:
            hi = mid
        else:
            lo = mid + 1
    I = lo
    X, nul, keys = kernel(r, I, C.P0)
    print(f"minimal I = {I} (t-degree {2*I}), nullity {nul} ({time.time()-t0:.0f}s)", flush=True)
    assert nul == 1
    ps = D.primes_below(400, C.P0 + 1)
    acc, mod, prev = None, 1, None
    for n_, p in enumerate(ps):
        X, nul, keys = kernel(r, I, p)
        assert nul == 1
        v = [int(X[k, 0]) for k in range(len(keys))]
        piv = max(k for k, x in enumerate(v) if x)
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
        if n_ % 20 == 19:
            print(f"  {n_+1} primes: {sum(x is not None for x in rr)}/{len(rr)} ({time.time()-t0:.0f}s)", flush=True)
    c = dict(zip(keys, rr))
    # exact x-identity: sum b T^i theta^j g - sum q T^i h4 == 0
    cg = {(i, j): v for (k, i, j), v in c.items() if k == "g" and v}
    ok_g = D.exact_check(A.path_module(), ("c0", "h0"), list(cg), cg, r, I) if False else None
    # do the exact check with both parts: build g-part and h-part vectors exactly
    exact = exact_identity(cg, {i: v for (k, i, j), v in c.items() if k == "h" and v}, r, I)
    print("exact x-identity:", exact, f"({time.time()-t0:.0f}s)", flush=True)
    Bq = [sp.Integer(0)] * (I + 1)
    for (i, j), v in cg.items():
        Bq[i] += sp.Rational(v.numerator, v.denominator) * (2 - u) ** j
    qq = sum(sp.Rational(v.numerator, v.denominator) * t ** (2 * i) for (k, i, j), v in c.items() if k == "h" and v)
    rec["B_Q"] = [str(sp.expand(x)) for x in Bq]
    rec["q"] = str(sp.factor(qq))
    rec["I"] = I
    rec["exact_x_identity"] = exact
    print("q(t^2) =", sp.factor(qq), flush=True)
    # A4: Mellin operator of P_h4 = N t^3 int K0(t x) x^6 K0^3 dx
    vK3 = Bm.sym_power([z, 0], 2, 3)
    A4 = Bm.mellin_t_ode(vK3, 4, 3, lambda sg: sg**2)
    rec["A4_Q"] = [str(q) for q in A4]
    print("A4 (Mellin of P_h4):", A4, flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


def exact_identity(cg, ch, r, I):
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

    def addto(acc, E, cf):
        for a in range(n):
            for k, v in E[a].items():
                acc[a][k] = acc[a].get(k, 0) + cf * v

    acc = [dict() for _ in range(n)]
    g = [dict() for _ in range(n)]
    g[0][0] = Fraction(1)
    th = [g]
    for _ in range(r):
        th.append(theta(th[-1]))
    for j in range(r + 1):
        cur = th[j]
        for i in range(I + 1):
            if cg.get((i, j)):
                addto(acc, cur, cg[(i, j)])
            if i < I:
                cur = T(cur)
    idx = {m: k for k, m in enumerate(basis)}
    h = [dict() for _ in range(n)]
    for m, cf in Oh.H4.items():
        h[idx[tuple(sorted(m))]][0] = Fraction(cf)
    cur = h
    for i in range(I + 1):
        if ch.get(i):
            addto(acc, cur, ch[i])
        if i < I:
            cur = T(cur)
    return all(v == 0 for d in acc for v in d.values())


if __name__ == "__main__":
    main()
