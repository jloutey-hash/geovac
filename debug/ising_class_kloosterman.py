"""Quadratic-twisted symmetric-power Kloosterman moments tm_k(p) (Lim-Tu-Yu, arXiv:2609.13576, sec 1):
  tm_k(p) = sum_{a in F_p^*} (a/p) Kl2^{Sym^k}(p;a),  Kl2(p;a) = sum_x e(x + a/x),  alpha+beta = -Kl2,
  alpha beta = p, Sym^k = sum alpha^{k-i} beta^i.
Used to (i) check Theorem thm:EC for 8.4.a.a (tm_4(p) = -a_f(p) p, f = eta(2t)^4 eta(4t)^4) and
48.3.e.a (twist of eta(2t)^3 eta(6t)^3 by chi_{-4}), and (ii) obtain a_p of 24.6.a.a from
tm_6(p) = -a_f(p) p - p^3 (p >= 5), checked against the printed a(5) = -34, a(7) = -240."""
import json
import numpy as np
import sympy as sp

PMAX = 620

def tm(p, ks):
    x = np.arange(1, p)
    inv = np.array([pow(int(v), -1, p) for v in x])
    res = {}
    # Kl(a) for all a: sum_x cos(2 pi (x + a x^{-1})/p); real for prime p
    A = np.arange(1, p)[:, None]
    ph = (x[None, :] + (A * inv[None, :]) % p) % p
    Kl = np.cos(2 * np.pi * ph / p).sum(axis=1)
    leg = np.array([sp.legendre_symbol(int(a), p) for a in range(1, p)], dtype=float)
    e = -Kl
    S = [np.ones_like(e), e.copy()]
    for k in range(2, max(ks) + 1):
        S.append(e * S[-1] - p * S[-2])
    for k in ks:
        v = float((leg * S[k]).sum())
        res[k] = int(round(v)); res[f"{k}_err"] = abs(v - round(v))
    return res

def eta_product(spec, M):
    """q-expansion coefficients a_1..a_M of prod eta(d tau)^e (integer weight, integral q-order)."""
    E = [0] * (M + 1); 
    # prod (1-q^n) via Euler pentagonal
    for k in range(-60, 61):
        g = k * (3 * k - 1) // 2
        if 0 <= g <= M:
            E[g] += (-1) ** (k % 2)
    def scale(ser, d):
        out = [0] * (M + 1)
        for i, c in enumerate(ser):
            if i * d <= M:
                out[i * d] += c
        return out
    def mul(a, b):
        out = [0] * (M + 1)
        for i, ai in enumerate(a):
            if ai:
                for j in range(M + 1 - i):
                    out[i + j] += ai * b[j]
        return out
    # inverse of E (partition numbers) by recursion
    P = [1] + [0] * M
    for n in range(1, M + 1):
        P[n] = -sum(E[j] * P[n - j] for j in range(1, n + 1))
    ser = [1] + [0] * M
    order24 = 0
    for d, e in spec:
        base = scale(E if e > 0 else P, d)
        for _ in range(abs(e)):
            ser = mul(ser, base)
        order24 += d * e
    assert order24 % 24 == 0
    sh = order24 // 24
    coeffs = [0] * (M + 1)
    for i, c in enumerate(ser):
        if i + sh <= M:
            coeffs[i + sh] = c
    return coeffs  # coeffs[n] = a_n

if __name__ == "__main__":
    primes = list(sp.primerange(3, PMAX))
    f8 = eta_product([(2, 4), (4, 4)], PMAX)
    f12 = eta_product([(2, 3), (6, 3)], PMAX)
    out = {"tm": {}, "checks": {}}
    a24 = {}
    bad8, bad48, notint, ramanujan = [], [], [], []
    maxerr = 0.0
    for p in primes:
        r = tm(p, [3, 4, 6])
        maxerr = max(maxerr, r["3_err"], r["4_err"], r["6_err"])
        out["tm"][p] = {k: r[k] for k in (3, 4, 6)}
        if r[4] != -f8[p] * p:
            bad8.append(p)
        if p >= 5:
            chi4 = 1 if p % 4 == 1 else -1
            if r[3] != -(chi4 * f12[p]) * p:
                bad48.append(p)
            num = -(r[6] + p**3)
            if num % p:
                notint.append(p)
            ap = num // p
            a24[p] = ap
            if ap * ap > 4 * p**5:
                ramanujan.append(p)
    out["checks"] = {"max_rounding_err": maxerr,
                     "8.4.a.a_mismatch_primes": bad8, "48.3.e.a_mismatch_primes": bad48,
                     "24.6.a.a_nonintegral": notint, "24.6.a.a_ramanujan_violations": ramanujan,
                     "a24(5)": a24.get(5), "a24(7)": a24.get(7), "a24(11)": a24.get(11), "a24(13)": a24.get(13)}
    out["a24_primes"] = a24
    print(out["checks"])
    json.dump(out, open("debug/data/ising_class_kloosterman.json", "w"), indent=1)
