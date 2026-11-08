"""IB4: lift the PF kite's differential system (c, A) from several primes to Q, build the operator for phi
exactly over Q(z) by the cyclic vector, check it against every prime's operator, then apply it in floating
point to the stored 300-coefficient series.

usage: python ib_lift_system.py <series json> <op json with c, A> ...
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import mpmath as mp
from flint import fmpq_poly, fmpq, fmpz

from ib_crt import crt, ratrec, series_residual

sys.set_int_max_str_digits(10**6)


class RQ:
    """rational function over Q, (num, den) fmpq_poly, den monic."""
    __slots__ = ("n", "d")

    def __init__(self, n, d=None):
        self.n = n
        self.d = d if d is not None else fmpq_poly([1])
        self._norm()

    def _norm(self):
        if self.n == 0:
            self.d = fmpq_poly([1])
            return
        g = self.n.gcd(self.d)
        if g.degree() > 0:
            self.n, self.d = self.n // g, self.d // g
        lc = self.d.coeffs()[-1]
        if lc != 1:
            self.n, self.d = self.n / lc, self.d / lc

    def __add__(s, o):
        return RQ(s.n * o.d + o.n * s.d, s.d * o.d)

    def __sub__(s, o):
        return RQ(s.n * o.d - o.n * s.d, s.d * o.d)

    def __mul__(s, o):
        return RQ(s.n * o.n, s.d * o.d)

    def __truediv__(s, o):
        return RQ(s.n * o.d, s.d * o.n)

    def diff(s):
        return RQ(s.n.derivative() * s.d - s.n * s.d.derivative(), s.d * s.d)

    def is_zero(s):
        return s.n == 0


def lift_poly(lists, mods):
    deg = max(len(l) for l in lists)
    out = []
    for i in range(deg):
        a, M = crt([l[i] if i < len(l) else 0 for l in lists], mods)
        q = ratrec(a, M)
        if q is None:
            return None
        out.append(q)
    return out


def to_fq(co):
    return fmpq_poly([fmpq(q.numerator, q.denominator) for q in co]) if co else fmpq_poly([0])


def lift_rf(entries, mods):
    """entries: per prime [num coeffs, den coeffs] (den monic)."""
    n = lift_poly([e[0] for e in entries], mods)
    d = lift_poly([e[1] for e in entries], mods)
    if n is None or d is None:
        return None
    return RQ(to_fq(n), to_fq(d))


def check_rf(x: RQ, e, p):
    """x mod p equals the prime's (num, den)?"""
    def red(poly):
        return [(int(c.p) * pow(int(c.q), -1, p)) % p for c in poly.coeffs()]
    return red(x.n) == [v % p for v in e[0]] and red(x.d) == [v % p for v in e[1]]


def cyclic(c, A):
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
        rel = relation(vs)
        if rel is not None:
            return rel
        if len(vs) > n + 1:
            raise RuntimeError("no relation")


def relation(vs):
    r = len(vs) - 1
    n = len(vs[0])
    zero = RQ(fmpq_poly([0]))
    M = [[vs[k][i] for k in range(r)] + [zero - vs[r][i]] for i in range(n)]
    row = 0
    for col in range(r):
        p = next((i for i in range(row, n) if not M[i][col].is_zero()), None)
        if p is None:
            return None
        M[row], M[p] = M[p], M[row]
        pv = M[row][col]
        M[row] = [x / pv for x in M[row]]
        for i in range(n):
            if i != row and not M[i][col].is_zero():
                f = M[i][col]
                M[i] = [x - f * y for x, y in zip(M[i], M[row])]
        row += 1
    for i in range(row, n):
        if not M[i][r].is_zero():
            return None
    return [M[i][r] for i in range(r)] + [RQ(fmpq_poly([1]))]


def main():
    series = json.loads(Path(sys.argv[1]).read_text())
    ops = [json.loads(Path(a).read_text()) for a in sys.argv[2:]]
    ops = [o for o in ops if "A" in o]
    mods = [o["P"] for o in ops]
    print("primes with the system saved:", len(ops))
    hold = ops[-1]
    use, mods_use = ops[:-1], mods[:-1]
    nm = len(use[0]["c"])
    c = [lift_rf([o["c"][i] for o in use], mods_use) for i in range(nm)]
    A = [[lift_rf([o["A"][i][j] for o in use], mods_use) for j in range(nm)] for i in range(nm)]
    if any(x is None for x in c) or any(x is None for row in A for x in row):
        print("system reconstruction failed: more primes needed")
        return
    ok = all(check_rf(c[i], hold["c"][i], hold["P"]) for i in range(nm)) and \
        all(check_rf(A[i][j], hold["A"][i][j], hold["P"]) for i in range(nm) for j in range(nm))
    hb = max(max(abs(int(q.p)).bit_length(), int(q.q).bit_length()) for row in A + [c] for x in row
             for poly in (x.n, x.d) for q in poly.coeffs())
    print(f"system lifted from {len(use)} primes (largest height {hb} bits); held-out prime agrees: {ok}")
    if not ok:
        return
    a = cyclic(c, A)
    order = len(a) - 1
    den = fmpq_poly([1])
    for x in a:
        den = den * x.d // den.gcd(x.d)
    polys = [x.n * (den // x.d) for x in a]
    # clear to integers
    from math import lcm
    L = 1
    for p in polys:
        for q in p.coeffs():
            L = lcm(L, int(q.q))
    ipolys = [[int((q * L).p) for q in p.coeffs()] for p in polys]
    print("operator over Q: order", order, "degrees", [len(p) - 1 for p in ipolys],
          "largest coefficient", max(abs(v).bit_length() for p in ipolys for v in p), "bits")
    # agrees with each prime's operator (monic-normalised there)?
    agree = 0
    for o in ops:
        p = o["P"]
        lead = ipolys[-1][-1] % p
        il = pow(lead, -1, p)
        mine = [[v * il % p for v in poly] for poly in ipolys]
        theirs = [[v % p for v in poly] for poly in o["polys"]]
        agree += mine == theirs
    print("agrees with the operator at", agree, "of", len(ops), "primes")
    op_q = [[Fraction(v) for v in p] for p in ipolys]
    res = series_residual(op_q, series["phi_n"], dps=300)
    print("series residual (worst relative, 300 coefficients at 280 digits):", mp.nstr(res, 3))
    Path("data/ib_pf_op_Q.json").write_text(json.dumps(
        {"form": "sum_k polys[k](z) (d/dz)^k phi = 0, integer coefficients", "order": order,
         "polys": [[str(v) for v in p] for p in ipolys], "series_residual": mp.nstr(res, 5)}))


if __name__ == "__main__":
    main()
