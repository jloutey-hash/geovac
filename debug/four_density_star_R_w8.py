"""FD3 step R (copy of step N with another target): the exact map from F_4's operator into
Y5 = W9/R4 (step Q found F_4 as a sub of Y5 mod p, odd parity), giving an order-8 right factor
W8 = R4'' o R4 of W9 that skips the order-1 piece theta - r.  Below: step N's docstring.
FD1 step N: the exact map from F_4's operator into M14 (step M found it mod p, odd parity, in the
box (t^2-1)^4 (t^2-4)^3 (t^2-9)^4 (5t^2+4)^3 (5t^4+2t^2+2)^3, numerators t^[-10..74]).
(1) shrink the box mod p (greedy, one factor at a time, then the numerator range);
(2) the gauge A = sum_{j<4} a_j(t) theta^j, a_j = p_j(t)/den, mod several primes, CRT + rational
    reconstruction;
(3) exact check, no series: M14 o A = 0 modulo the right ideal of Q4 (F_4's operator), in the
    delta = theta/2, w = t^2 toolkit (A is odd: A = t * A_even, and t^-1 M14 t is M14 with
    theta -> theta + 1);
(4) R4' = annihilator of A(Sol Q4) (order 4, exact), M14 = N10 o R4' exactly.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
import four_density_star_ops as O  # noqa: E402

u = sp.Symbol("u")
t = sp.Symbol("t")
OUT = HERE / "data" / "four_density_star_R.json"
Q4 = [9 * u**2 * (u - 2) ** 2, -2 * u**2 * (5 * u**2 + 1), (u + 1) ** 4]
def _primes(n):
    out, q = [], (1 << 61) - 1
    while len(out) < n:
        if sp.isprime(q):
            out.append(q)
        q = sp.prevprime(q)
    return out


PRIMES = _primes(10)
FACT = {"b": t**2 - 1, "c": t**2 - 4, "e": t**2 - 9, "k1": 5 * t**2 + 4, "k2": 5 * t**4 + 2 * t**2 + 2}


def isprime_check():
    for p in PRIMES:
        assert sp.isprime(p), p


class ModP:
    def __init__(self, p):
        self.p = p

    def fr(self, x):
        x = Fraction(x)
        return x.numerator % self.p * pow(x.denominator % self.p, self.p - 2, self.p) % self.p

    def frob(self, Q, lam, NT):
        p = self.p
        Qp = [sp.Poly(q, u) for q in Q]
        cs = [[(m[0], self.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)))) for m, c in zip(qp.monoms(), qp.coeffs())] for qp in Qp]
        ev = lambda i, x: sum(cc * pow(x % p, m, p) for m, cc in cs[i]) % p
        c = [1]
        for n in range(1, NT // 2):
            s_ = sum(ev(i, lam + 2 * (n - i)) * c[n - i] for i in range(1, len(Q)) if n - i >= 0) % p
            dd = ev(0, lam + 2 * n)
            assert dd != 0
            c.append((-s_) * pow(dd, p - 2, p) % p)
        return {lam + 2 * n: cn for n, cn in enumerate(c) if cn}

    def theta(self, ser, j):
        p = self.p
        return {e: c * pow(e % p, j, p) % p for e, c in ser.items()}

    def inv_series(self, den, N):
        p = self.p
        cs = [self.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q))) for c in sp.Poly(den, t).all_coeffs()[::-1]]
        out = [pow(cs[0], p - 2, p)]
        for n in range(1, N):
            s_ = sum(cs[k] * out[n - k] for k in range(1, min(n, len(cs) - 1) + 1)) % p
            out.append((-s_) * out[0] % p)
        return {n: c for n, c in enumerate(out) if c}

    def mul(self, a, b, top):
        p = self.p
        out = {}
        for e1, c1 in a.items():
            for e2, c2 in b.items():
                if e1 + e2 < top:
                    out[e1 + e2] = (out.get(e1 + e2, 0) + c1 * c2) % p
        return out

    def apply_L(self, Qcs, ser):
        p = self.p
        out = {}
        for i, cs in enumerate(Qcs):
            for e, c in ser.items():
                val = sum(cc * pow(e % p, m, p) for m, cc in cs) % p
                if val:
                    out[e + 2 * i] = (out.get(e + 2 * i, 0) + val * c) % p
        return {e: c for e, c in out.items() if c}

    def prepQ(self, Q):
        return [[(m[0], self.fr(Fraction(int(sp.Rational(c).p), int(sp.Rational(c).q)))) for m, c in zip(sp.Poly(q, u).monoms(), sp.Poly(q, u).coeffs())] for q in Q]

    def system(self, Qcs, f, r, den, lo, hi, NT, top):
        w = self.inv_series(den, NT)
        cols, labels = [], []
        for j in range(r):
            g = self.mul(w, self.theta(f, j), NT + min(f))
            for d in range(lo, hi + 1, 2):
                cols.append(self.apply_L(Qcs, {e + d: v for e, v in g.items()}))
                labels.append((j, d))
        exps = sorted(x for x in set().union(*cols) if x < top)
        rows = [[col.get(x, 0) for col in cols] for x in exps]
        return rows, labels

    def nullspace(self, rows, ncols):
        p = self.p
        M = [r[:] for r in rows]
        piv_cols, rk = [], 0
        for col in range(ncols):
            piv = next((i for i in range(rk, len(M)) if M[i][col]), None)
            if piv is None:
                continue
            M[rk], M[piv] = M[piv], M[rk]
            iv = pow(M[rk][col], p - 2, p)
            M[rk] = [x * iv % p for x in M[rk]]
            for i in range(len(M)):
                if i != rk and M[i][col]:
                    fct = M[i][col]
                    M[i] = [(x - fct * y) % p for x, y in zip(M[i], M[rk])]
            piv_cols.append(col)
            rk += 1
        free = [c for c in range(ncols) if c not in piv_cols]
        vecs = []
        for fcol in free:
            v = [0] * ncols
            v[fcol] = 1
            for i, pc in enumerate(piv_cols):
                v[pc] = (-M[i][fcol]) % p
            vecs.append(v)
        return vecs


def make_den(pw):
    d = sp.Integer(1)
    for k, e in pw.items():
        d *= FACT[k] ** e
    return sp.expand(d)


def nullity(mp_, Qcs, f, pw, lo, hi, NT):
    top = min(f) + NT - 2 * len(Qcs) - hi - 8
    rows, labels = mp_.system(Qcs, f, 4, make_den(pw), lo, hi, NT, top)
    return len(mp_.nullspace(rows, len(labels))), len(labels), len(rows)


def ratrec(a, m):
    """rational reconstruction of a mod m."""
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


def main() -> None:
    isprime_check()
    rec = {}
    t0 = time.time()
    import four_density_star_P_score as PS
    ops = PS.build()
    Y5, rem = O.rdiv(ops["W9"], ops["R4"])
    assert all(c == 0 for c in rem)
    M14Q, _ = O.to_Q(Y5)          # target (variable name kept from step N)
    rec["target"] = "Y5 = W9/R4"
    NT = 1100
    m0 = ModP(PRIMES[0])
    Qcs = m0.prepQ(M14Q)
    f4 = m0.frob(Q4, 2, NT)
    pw = {"b": 4, "c": 3, "e": 4, "k1": 3, "k2": 3}
    lo, hi = -9, 73                      # odd parity
    if "fast" in sys.argv:               # the minimal box found by the greedy shrink of a previous run
        pw, lo, hi = {"b": 3, "c": 0, "e": 4, "k1": 1, "k2": 0}, 3, 23  # step N's box (not used here)
    n0 = nullity(m0, Qcs, f4, pw, lo, hi, NT)
    print("start box:", pw, (lo, hi), n0, f"({time.time()-t0:.0f}s)", flush=True)
    assert n0[0] >= 1
    # greedy shrink of the denominator powers, adjusting the numerator top by the degree removed
    for k in ([] if "fast" in sys.argv else ["k2", "k1", "e", "c", "b"]):
        while pw[k] > 0:
            trial = dict(pw)
            trial[k] -= 1
            dd = sp.degree(FACT[k], t)
            n = nullity(m0, Qcs, f4, trial, lo, hi - dd, NT)
            if n[0] >= 1:
                pw, hi = trial, hi - dd
                print("  shrink", k, "->", pw, (lo, hi), n, f"({time.time()-t0:.0f}s)", flush=True)
            else:
                break
    while "fast" not in sys.argv:
        n = nullity(m0, Qcs, f4, pw, lo, hi - 2, NT)
        if n[0] >= 1:
            hi -= 2
        else:
            break
    while "fast" not in sys.argv:
        n = nullity(m0, Qcs, f4, pw, lo + 2, hi, NT)
        if n[0] >= 1:
            lo += 2
        else:
            break
    # minimality of the box (each factor one power lower, and the numerator range one step narrower)
    shr = {}
    for k in pw:
        if pw[k] > 0:
            trial = dict(pw); trial[k] -= 1
            shr[k] = nullity(m0, Qcs, f4, trial, lo, hi, NT)[0]
    shr["hi-2"] = nullity(m0, Qcs, f4, pw, lo, hi - 2, NT)[0]
    shr["lo+2"] = nullity(m0, Qcs, f4, pw, lo + 2, hi, NT)[0]
    print("one-step-smaller boxes (nullity):", shr, flush=True)
    rec["smaller_boxes_nullity"] = shr
    nfin = nullity(m0, Qcs, f4, pw, lo, hi, NT)
    print("minimal box:", pw, (lo, hi), nfin, f"({time.time()-t0:.0f}s)", flush=True)
    rec["box"] = {"powers": pw, "lo": int(lo), "hi": int(hi), "nullity_unknowns_rows": nfin}
    OUT.write_text(json.dumps(rec, indent=1))
    # (2) multi-prime gauge
    den = make_den(pw)
    vecs, labels = [], None
    for p in PRIMES:
        mp_ = ModP(p)
        Qc = mp_.prepQ(M14Q)
        f = mp_.frob(Q4, 2, NT)
        top = min(f) + NT - 2 * len(Qc) - hi - 8
        rows, labels = mp_.system(Qc, f, 4, den, lo, hi, NT, top)
        ns = mp_.nullspace(rows, len(labels))
        assert len(ns) == 1, len(ns)
        v = ns[0]
        # normalize: last nonzero coordinate = 1
        piv = max(i for i, x in enumerate(v) if x)
        iv = pow(v[piv], p - 2, p)
        vecs.append([x * iv % p for x in v])
        # CRT + reconstruction
        mod = 1
        acc = [0] * len(v)
        for pp, vv in zip(PRIMES[:len(vecs)], vecs):
            for i in range(len(v)):
                # combine acc (mod mod) with vv[i] (mod pp)
                x = acc[i] + mod * ((vv[i] - acc[i]) * pow(mod, -1, pp) % pp)
                acc[i] = x
            mod *= pp
        acc = [x % mod for x in acc]
        rr = [ratrec(x, mod) for x in acc]
        print(f"  prime {len(vecs)}: reconstructed {sum(r is not None for r in rr)}/{len(rr)}  ({time.time()-t0:.0f}s)", flush=True)
        if all(r is not None for r in rr) and len(vecs) >= 2:
            # stability: same as with one fewer prime?
            if "prev" in locals() and prev == rr:
                break
            prev = rr
    a = [sp.Integer(0)] * 4
    for r_, (j, d) in zip(rr, labels):
        a[j] += sp.Rational(r_.numerator, r_.denominator) * t**d
    a = [sp.cancel(x / den) for x in a]
    rec["A"] = [str(x) for x in a]
    OUT.write_text(json.dumps(rec, indent=1))
    # (3) exact identity: A odd -> A = t * Ae with Ae even; M14 (t g) = t * M14shift g, theta -> theta + 1
    Ae = [sp.cancel(x / t) for x in a]
    AeD = O.from_theta_t([sp.Integer(0)] * 0 + [x for x in Ae]) if False else None
    # convert A_e = sum_j ae_j(t) theta^j to delta form: theta^j = 2^j delta^j
    AeD = []
    for j, x in enumerate(Ae):
        AeD += [O.from_theta_t([x])[0] * 2**j]
    M14sh = O.from_Q([sp.expand(q.subs(u, u + 1)) for q in M14Q])
    Q4op = O.from_Q(Q4)
    S4 = [-c / Q4op[4] for c in Q4op[:4]]

    def dvec(vv):
        c0, c1, c2, c3 = vv
        return [O.d(c0) + c3 * S4[0], c0 + O.d(c1) + c3 * S4[1], c1 + O.d(c2) + c3 * S4[2], c2 + O.d(c3) + c3 * S4[3]]
    acc, cur = [O.K(0)] * 4, AeD
    for m in M14sh:
        acc = [x + m * y for x, y in zip(acc, cur)]
        cur = dvec(cur)
    rec["M14_A_mod_Q4_zero"] = all(x == 0 for x in acc)
    print("exact: target o A = 0 mod Q4:", rec["M14_A_mod_Q4_zero"], f"({time.time()-t0:.0f}s)", flush=True)
    # (4) R4' on the shifted (even) side: annihilator of Ae(Sol Q4) for M14sh; then conjugate back
    vecs4 = [AeD]
    for _ in range(4):
        vecs4.append(dvec(vecs4[-1]))
    # nullspace of the 4 x 5 matrix over Q(w) by hand (Gaussian elimination on FracElements)
    Mx = [[vecs4[k][i] for k in range(5)] for i in range(4)]
    piv = []
    rk = 0
    for col in range(5):
        pr = next((i for i in range(rk, 4) if Mx[i][col] != 0), None)
        if pr is None:
            continue
        Mx[rk], Mx[pr] = Mx[pr], Mx[rk]
        iv = 1 / Mx[rk][col]
        Mx[rk] = [x * iv for x in Mx[rk]]
        for i in range(4):
            if i != rk and Mx[i][col] != 0:
                fct = Mx[i][col]
                Mx[i] = [x - fct * y for x, y in zip(Mx[i], Mx[rk])]
        piv.append(col)
        rk += 1
    free = [c for c in range(5) if c not in piv]
    assert free == [4], free
    row = [O.K(0)] * 5
    row[4] = O.K(1)
    for i, pc in enumerate(piv):
        row[pc] = -Mx[i][4]
    R4p_sh = row
    N10sh, rem = O.rdiv(M14sh, R4p_sh)
    rec["M14_div_R4p_zero"] = all(x == 0 for x in rem)
    print("target = N o R4'' exact (shifted frame):", rec["M14_div_R4p_zero"], f"({time.time()-t0:.0f}s)", flush=True)
    # back to the unshifted frame: R4' = t R4p_sh t^-1 = R4p_sh with delta -> delta - 1/2
    def shift(Lop, s_):
        out = [O.K(0)] * len(Lop)
        for k, c in enumerate(Lop):
            for j in range(k + 1):
                out[j] += c * sp.binomial(k, j) * sp.Rational(s_) ** (k - j)
        return out
    R4p = shift(R4p_sh, sp.Rational(-1, 2))
    N10 = shift(N10sh, sp.Rational(-1, 2))
    R4pQ, D1 = O.to_Q(R4p)
    N10Q, D2 = O.to_Q(N10)
    rec["R4p_Q"] = [str(q) for q in R4pQ]
    rec["N10_Q"] = [str(q) for q in N10Q]
    print("R4' Q-form:", len(R4pQ), "t-terms; leading/t^0:", sp.factor(R4pQ[0]), flush=True)
    OUT.write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
