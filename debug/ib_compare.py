"""IB record: compare an IBP operator (d/dz form mod P) with a stored kite operator (sum_i z^i Q_i(theta_z)).

usage: python ib_compare.py <ib op json> <stored op json>
Exact over GF(P)(z): equality after making both monic in D, else right division of the IBP operator by the
stored one (Ore division in d/dz), reporting the remainder."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp
from flint import nmod_poly, nmod

sys.path.insert(0, str(Path(__file__).parent))
from ib_reduce import P  # noqa: E402
from ib_kite_de import RF  # noqa: E402

u = sp.Symbol("u")


def stored_to_D(path):
    d = json.loads(Path(path).read_text())
    Q = [sp.Poly(sp.sympify(s.replace("^", "**"), locals={"u": u}), u) for s in d["Q"]]
    r = max(q.degree() for q in Q)
    # theta^k = sum_j S2(k, j) z^j D^j
    coef = {}  # (j, power of z) -> integer
    for i, q in enumerate(Q):
        for (k,), c in q.terms():
            for j in range(k + 1):
                s2 = sp.functions.combinatorial.numbers.stirling(k, j, kind=2)
                if s2:
                    coef[(j, i + j)] = coef.get((j, i + j), 0) + int(c) * int(s2)
    polys = []
    for j in range(r + 1):
        deg = max((e for (jj, e) in coef if jj == j), default=0)
        polys.append(nmod_poly([coef.get((j, e), 0) % P for e in range(deg + 1)], P))
    return [RF(p) for p in polys]


def monic(op):
    lead = op[-1]
    return [x / lead for x in op]


def compose_D_left(a, op):
    """a(z) * op as operator (a RF), op list of RF."""
    return [a * x for x in op]


def D_times(op):
    """D o op: D (sum c_k D^k) = sum (c_k' D^k + c_k D^(k+1))."""
    out = [RF(nmod_poly([0], P)) for _ in range(len(op) + 1)]
    for k, c in enumerate(op):
        out[k] = out[k] + c.diff()
        out[k + 1] = out[k + 1] + c
    return out


def right_divide(L, R):
    """L = Q o R + rem, operators as lists of RF (index = power of D)."""
    rem = list(L)
    zero = RF(nmod_poly([0], P))
    while len(rem) >= len(R):
        while rem and rem[-1].is_zero():
            rem.pop()
        if len(rem) < len(R):
            break
        k = len(rem) - len(R)
        q = rem[-1] / R[-1]
        t = list(R)
        for _ in range(k):
            t = D_times(t)
        t = compose_D_left(q, t)
        rem = [rem[i] - (t[i] if i < len(t) else zero) for i in range(len(rem))]
        while rem and rem[-1].is_zero():
            rem.pop()
    return rem


def main():
    ib = json.loads(Path(sys.argv[1]).read_text())
    L = [RF(nmod_poly(p, P)) for p in ib["polys"]]
    S = stored_to_D(sys.argv[2])
    print("IBP order", len(L) - 1, " stored order", len(S) - 1)
    if len(L) == len(S):
        a, b = monic(L), monic(S)
        same = all((x - y).is_zero() for x, y in zip(a, b))
        print("equal after making monic:", same)
    rem = right_divide(L, S)
    print("stored right-divides IBP operator:", len(rem) == 0, " remainder order", len(rem) - 1 if rem else None)


if __name__ == "__main__":
    main()
