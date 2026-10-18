"""QG step 6h (gamma = 3/2): right factors of L9 through Frobenius solutions at a finite singular
point t0 (mod p): local variable h = t - t0, theta_h form, exponent e (rational); smallest order k of a
polynomial-coefficient operator (h-degree <= Dg) killing the solution.  k = 1: a hyperexponential
(algebraic / elementary) right factor; k = 2: an order-2 piece, etc.
usage: python ring_ext_QG_6h_localfrob.py <op-json> <key> <adj|-> <t0> <e> [Dg list]"""
from __future__ import annotations
import json, sys
from fractions import Fraction
from pathlib import Path
import sympy as sp
from sympy.functions.combinatorial.numbers import stirling
sys.path.insert(0, str(Path(__file__).parent))
import density_loop_C_factor as C
import density_loop_D_gauge as D

P = C.P
u = sp.Symbol("u")


def main():
    path, key, adj, t0s, es = sys.argv[1:6]
    Dgs = [int(x) for x in (sys.argv[6].split(",") if len(sys.argv) > 6 else ["10", "20", "40"])]
    Q = [sp.sympify(q, locals={"u": u}) for q in json.loads(Path(path).read_text())[key]]
    if adj == "adj":
        Q = [sp.expand(q.subs(u, -u - 2 * i)) for i, q in enumerate(Q)]
    r = max(sp.Poly(q, u).degree() for q in Q)
    d = len(Q) - 1
    # D-form a_k(t) mod p (coefficient lists in t)
    A = [[0] * (2 * d + r + 1) for _ in range(r + 1)]
    for i, q in enumerate(Q):
        for (j,), c in sp.Poly(q, u).terms():
            cm = D.fr(Fraction(int(c)))
            for k in range(j + 1):
                A[k][2 * i + k] = (A[k][2 * i + k] + cm * int(stirling(j, k, kind=2))) % P
    x0 = D.fr(Fraction(t0s))

    def shift(a):
        a = a[:]
        out = []
        while a:
            acc, q = 0, [0] * (len(a) - 1)
            for idx in range(len(a) - 1, -1, -1):
                acc = (acc * x0 + a[idx]) % P
                if idx > 0:
                    q[idx - 1] = acc
            out.append(acc)
            a = q
        return out
    Bh = [shift(ak) for ak in A]          # a_k(t0 + h) coefficients in h
    ords = [next((m for m, c in enumerate(b) if c), 10**9) for b in Bh]
    m = max(k - ords[k] for k in range(r + 1))
    # operator sum_s h^s B_s(theta): B_s(x) = sum_k b_{k, s - m + k} ff_k(x)
    maxs = max(len(b) for b in Bh) + m
    def ffk(k, x):
        v = 1
        for i_ in range(k):
            v = v * (x - i_) % P
        return v

    def Bs(s, x):
        tot = 0
        for k in range(r + 1):
            idx = s - m + k
            if 0 <= idx < len(Bh[k]) and Bh[k][idx]:
                tot = (tot + Bh[k][idx] * ffk(k, x)) % P
        return tot
    e = D.fr(Fraction(es))
    print("indicial check B_0(e) =", Bs(0, e), flush=True)
    NT = 700
    c = [1]
    for n in range(1, NT):
        s_ = sum(Bs(s, (e + n - s) % P) * c[n - s] for s in range(1, min(n, maxs) + 1)) % P
        dd = Bs(0, (e + n) % P)
        if dd == 0:
            print("resonance at n =", n, flush=True)
            return
        c.append((-s_) * pow(dd, P - 2, P) % P)
    for Dg in Dgs:
        found = None
        for k in range(1, r + 1):
            cols = []
            for dd_ in range(Dg + 1):
                for j in range(k + 1):
                    cols.append([(pow((e + n - dd_) % P, j, P) * c[n - dd_] % P) if n - dd_ >= 0 else 0 for n in range(NT)])
            U = len(cols)
            if U + 40 > NT:
                break
            rows = [[col[n] for col in cols] for n in range(NT)]
            if U - C.rank_mod(rows) > 0:
                found = k
                break
        print(f"t0 = {t0s}, exponent {es}: minimal order with h-degree <= {Dg}: {found}", flush=True)
        if found == 1:
            break


if __name__ == "__main__":
    main()
