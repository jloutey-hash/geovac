"""FP step D: the order-12 operator W12 (t-degree 78) from the x-module relation found in step C, exactly:
nullspace mod several primes (nullity 1 each), CRT + rational reconstruction until stable, then an exact check
of the x-module identity w g = 0 over Q (Fractions).  W12 = sum_i t^(2i) q_i(theta_t), q_i(u) = sum_j c_ij (2-u)^j."""
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

u = A.u
OUT = HERE / "data" / "four_density_path_D.json"


def primes_below(n, start):
    out, q = [], start
    while len(out) < n:
        q = sp.prevprime(q)
        out.append(q)
    return out


def ratrec(a, m):
    a %= m
    if a == 0:
        return Fraction(0)
    r0, r1, s0, s1 = m, a, 0, 1
    bound = int((m // 2) ** 0.5)
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if abs(s1) > bound or s1 == 0:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def kernel_vector(M, start, r, I, p):
    basis, A0, A1 = C.module_matrices(M, start)
    mp_ = C.ModP(A0, A1, p, -I - 2, I + r + 4)
    cols = C.columns(mp_, r, I)
    keys = [(i, j) for i in range(I + 1) for j in range(r + 1)]
    mats = np.stack([cols[k].reshape(-1) for k in keys], axis=1)
    mats = mats[np.nonzero(mats.any(axis=1))[0]]
    Mm = flint.nmod_mat(mats.shape[0], mats.shape[1], [int(x) for x in mats.reshape(-1)], p)
    X, nul = Mm.nullspace()
    assert nul == 1, nul
    v = [int(X[i, 0]) for i in range(len(keys))]
    piv = max(i for i, x in enumerate(v) if x)
    iv = pow(v[piv], p - 2, p)
    return keys, [x * iv % p for x in v], piv


def exact_check(M, start, keys, c, r, I):
    """sum c_ij T^i theta^j g == 0 over Q, exactly (module vectors with Fraction Laurent coefficients)."""
    basis, A0, A1 = C.module_matrices(M, start)
    n = len(basis)
    A0 = [[Fraction(int(sp.Rational(x).p), int(sp.Rational(x).q)) for x in row] for row in A0]
    A1 = [[Fraction(int(sp.Rational(x).p), int(sp.Rational(x).q)) for x in row] for row in A1]

    def theta(E):  # E: list (per basis) of dict power -> Fraction
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

    def sub(E, F):
        return [{k: E[a].get(k, 0) - F[a].get(k, 0) for k in set(E[a]) | set(F[a])} for a in range(n)]

    def T(E):
        F = sub(theta(E), E)
        F = sub(theta(F), F)
        return [{k - 1: v for k, v in d.items() if v} for d in F]

    g = [dict() for _ in range(n)]
    g[0][0] = Fraction(1)
    th = [g]
    for _ in range(r):
        th.append(theta(th[-1]))
    acc = [dict() for _ in range(n)]
    for j in range(r + 1):
        cur = th[j]
        for i in range(I + 1):
            cf = c.get((i, j), 0)
            if cf:
                for a in range(n):
                    for k, v in cur[a].items():
                        acc[a][k] = acc[a].get(k, 0) + cf * v
            if i < I:
                cur = T(cur)
    return all(v == 0 for d in acc for v in d.values())


def to_Q(keys, c, I):
    Q = [sp.Integer(0)] * (I + 1)
    for (i, j) in keys:
        cf = c.get((i, j), 0)
        if cf:
            Q[i] += sp.Rational(cf.numerator, cf.denominator) * (2 - u) ** j
    return [sp.expand(q) for q in Q]


def solve(M, start, r, I, label, nprimes_max=200):
    t0 = time.time()
    ps = primes_below(nprimes_max, C.P0 + 1)
    mod, acc, prev = 1, None, None
    for n_, p in enumerate(ps):
        keys, v, piv = kernel_vector(M, start, r, I, p)
        if acc is None:
            acc = v
            mod = p
            piv0 = piv
        else:
            assert piv == piv0
            acc = [(a + mod * ((b - a) * pow(mod, -1, p) % p)) for a, b in zip(acc, v)]
            mod *= p
            acc = [a % mod for a in acc]
        rr = [ratrec(a, mod) for a in acc]
        ok = all(x is not None for x in rr)
        if ok and prev == rr:
            print(f"[{label}] stable after {n_ + 1} primes ({time.time()-t0:.0f}s)", flush=True)
            break
        prev = rr if ok else None
        if n_ % 10 == 9:
            print(f"[{label}] {n_ + 1} primes, reconstructed {sum(x is not None for x in rr)}/{len(rr)} ({time.time()-t0:.0f}s)", flush=True)
    c = dict(zip(keys, rr))
    ex = exact_check(M, start, keys, c, r, I)
    print(f"[{label}] exact x-module identity: {ex} ({time.time()-t0:.0f}s)", flush=True)
    return to_Q(keys, c, I), ex


def main() -> None:
    rec = json.loads(OUT.read_text()) if OUT.exists() else {}
    if "chain_W4" not in rec:
        chi_pt = {"c0": {"c1": 1}, "c1": {"c2": 1}, "c2": {"c1": z - 3, "c2": 4}}
        Q, ex = solve(A.Module(chi_pt), ("c0", "c0"), 4, 4, "chain")
        rec["chain_W4"] = {"Q": [str(q) for q in Q], "exact": ex}
        print("chain W4:", [sp.factor(q) for q in Q], flush=True)
        OUT.write_text(json.dumps(rec, indent=1))
    if "W12" not in rec:
        Q, ex = solve(A.path_module(), ("c0", "h0"), 12, 39, "path")
        rec["W12"] = {"Q": [str(q) for q in Q], "exact_x_identity": ex}
        OUT.write_text(json.dumps(rec, indent=1))
        print("W12 Q_0:", sp.factor(Q[0]), flush=True)
        print("W12 Q_top:", sp.factor(Q[-1]), flush=True)


z = A.z
if __name__ == "__main__":
    main()
