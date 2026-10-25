"""IB record: lift an IBP operator from several primes to Q (CRT + rational reconstruction), check it on one
more prime, and apply it in floating point to a stored Taylor series of the kite.

usage: python ib_crt.py <series json> <op json prime 1> <op json prime 2> ... [--check <op json>]
Each op json holds polys (d/dz form, monic top coefficient polynomial, as written by ib_kite_de.py).
The series check: c_n = phi_n; the coefficient of z^N in sum_k p_k(z) D^k phi, relative to the sum of the
absolute values of its terms, for every N the series allows."""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp

sys.set_int_max_str_digits(10**6)


def crt(residues, moduli):
    x, M = 0, 1
    for r, m in zip(residues, moduli):
        t = ((r - x) * pow(M, -1, m)) % m
        x, M = x + M * t, M * m
    return x % M, M


def ratrec(a, M):
    import math
    bound = math.isqrt(M // 2)
    r0, r1, s0, s1 = M, a % M, 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1, s0, s1 = r1, r0 - q * r1, s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)


def lift(ops):
    order = len(ops[0]["polys"]) - 1
    assert all(len(o["polys"]) - 1 == order for o in ops), "orders differ between primes"
    moduli = [o["P"] for o in ops]
    out = []
    for k in range(order + 1):
        deg = max(len(o["polys"][k]) for o in ops)
        co = []
        for i in range(deg):
            res = [o["polys"][k][i] if i < len(o["polys"][k]) else 0 for o in ops]
            a, M = crt(res, moduli)
            q = ratrec(a, M)
            if q is None:
                return None
            co.append(q)
        out.append(co)
    return out


def check_prime(op_q, o):
    p = o["P"]
    for k, co in enumerate(op_q):
        ref = o["polys"][k]
        for i, q in enumerate(co):
            r = ref[i] if i < len(ref) else 0
            if (q.numerator - r * q.denominator) % p:
                return False
    return True


def series_residual(op_q, phis, dps=300):
    mp.mp.dps = dps
    x = [mp.mpf(s) for s in phis]
    Nmax = len(x)
    order = len(op_q) - 1
    worst = mp.mpf(0)
    for N in range(0, Nmax - order - 1):
        acc = mp.mpf(0)
        size = mp.mpf(0)
        for k, co in enumerate(op_q):
            for i, c in enumerate(co):
                m = N - i
                if m < 0 or c == 0 or m + k >= Nmax:
                    continue
                ff = 1
                for j in range(k):
                    ff *= (m + k - j)
                term = mp.mpf(c.numerator) / c.denominator * ff * x[m + k]
                acc += term
                size += abs(term)
        if size:
            worst = max(worst, abs(acc) / size)
    return worst


def main():
    args = sys.argv[1:]
    check = None
    if "--check" in args:
        i = args.index("--check")
        check = json.loads(Path(args[i + 1]).read_text())
        args = args[:i] + args[i + 2:]
    series = json.loads(Path(args[0]).read_text())
    ops = [json.loads(Path(a).read_text()) for a in args[1:]]
    op_q = lift(ops)
    if op_q is None:
        print("rational reconstruction failed: more primes needed")
        return
    hmax = max(max(abs(q.numerator), q.denominator).bit_length() for co in op_q for q in co)
    print("lifted from", len(ops), "primes; largest numerator/denominator", hmax, "bits")
    if check is not None:
        print("agrees with the extra prime:", check_prime(op_q, check))
    phis = series["phi_n"]
    print("series residual (worst relative, over", len(phis), "coefficients):",
          mp.nstr(series_residual(op_q, phis), 3))
    Path(args[1]).with_name(Path(args[1]).stem.split("_p")[0] + "_Q.json").write_text(
        json.dumps({"form": "sum_k polys[k](z) (d/dz)^k phi = 0 over Q, monic top polynomial",
                    "polys": [[str(q) for q in co] for co in op_q]}))


if __name__ == "__main__":
    main()
