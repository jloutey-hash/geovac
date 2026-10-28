"""Guard for Paper 59 sec:density_ring, "the operators are also evaluators" (record
docs/predictions/2026-10-26_operator_evaluator.md, OE1-OE2).

The off-slice kite phi(z) at (1, 2; 3/2, 3) from its exact order-11 operator (tests/data/paper59_ib_2026_10_24.json,
pf_op) and its first six Taylor coefficients at z = 0 as the only input: the series at the singular point z = 0, then
Taylor steps of half the distance to the nearest root of the leading coefficient, apparent roots included, along a
path that passes above the apparent root on the positive axis.  Claims guarded:
  * six coefficients and the operator fix the series (coefficients 6..20 against the independently computed series);
  * the continuation reaches grid nodes past the series' disk and past the apparent root and agrees with the
    60-digit quadrature grid behind Y (two nodes, z = 0.500 and z = 1.0; tests/data/paper59_kite_p2_2026_10_27.json).
The evaluator below is a copy of debug/oe_engine.py (no debug/ import).
"""
import json
import math
from pathlib import Path

from flint import acb, acb_poly, arb, ctx, fmpq_poly, fmpz_poly

HERE = Path(__file__).parent / "data"
IB = json.loads((HERE / "paper59_ib_2026_10_24.json").read_text())
NODES = json.loads((HERE / "paper59_kite_p2_2026_10_27.json").read_text())["grid_nodes"]


def _mid(x):
    return acb(x.real.mid(), x.imag.mid()) if isinstance(x, acb) else x.mid()


def _ff(n, k):
    out = 1
    for i in range(k):
        out *= n - i
    return out


class Operator:
    def __init__(self, polys_int):
        self.order = len(polys_int) - 1
        self.q = [fmpq_poly([int(v) for v in p]) for p in polys_int]

    def series_at_zero(self, free, N):
        den = 1
        for p in self.q:
            den = den * int(p.denom()) // math.gcd(den, int(p.denom()))
        cols = [[int(c * den) for c in p.coeffs()] for p in self.q]
        vals = [next(i for i, c in enumerate(cs) if c != 0) for cs in cols]
        delta = min(vals[k] - k for k in range(self.order + 1))
        emax = max(len(cs) - 1 - k - delta for k, cs in enumerate(cols))
        Q = [[(cs[k + delta + e] if 0 <= k + delta + e < len(cs) else 0) for k, cs in enumerate(cols)]
             for e in range(emax + 1)]
        a, freed = [], []
        for n in range(N):
            s = arb(0)
            for e in range(1, min(len(Q), n + 1)):
                qe = sum(Q[e][k] * _ff(n - e, k) for k in range(self.order + 1))
                if qe:
                    s += a[n - e] * qe
            q0 = sum(Q[0][k] * _ff(n, k) for k in range(self.order + 1))
            if q0 == 0:
                freed.append(n)
                a.append(_mid(arb(free[n])))
            else:
                a.append(_mid(-s / q0))
        return a, freed

    def local(self, c, h):
        sub = acb_poly([c, h])
        r = []
        for k, p in enumerate(self.q):
            sh = acb_poly([acb(arb(cf)) for cf in p.coeffs()])(sub)
            r.append([cf * h ** (self.order - k) for cf in sh.coeffs()])
        lead = r[self.order][0]
        return [[_mid(x / lead) for x in row] for row in r]


def step_series(r, order, init, N):
    a = [_mid(x) for x in init]
    b = [[] for _ in range(order + 1)]
    lead = r[order][0]
    degs = [len(row) for row in r]

    def push(n):
        for k in range(order + 1):
            m = n - k
            if m >= 0 and len(b[k]) == m:
                b[k].append(a[n] * _ff(n, k))
    for n in range(len(a)):
        push(n)
    for m in range(0, N - order):
        s = acb(0)
        for k in range(order + 1):
            for j in range(1 if k == order else 0, min(degs[k], m + 1)):
                s += r[k][j] * b[k][m - j]
        a.append(_mid(-s / (lead * _ff(m + order, order))))
        push(m + order)
    return a


def eval_series(a, v, nder=0):
    out = []
    for d in range(nder + 1):
        s = acb(0)
        for n in range(len(a) - 1, d - 1, -1):
            s = _mid(s * v + a[n] * _ff(n, d))
        out.append(s)
    return out


def _roots(polys):
    out = []
    for f, _ in fmpz_poly([int(v) for v in polys[-1]]).factor()[1]:
        out += [complex(r[0]) for r in f.complex_roots()]
    return out


def walk(op, roots, c, der, T):
    """Taylor steps from c to T, each of half the distance to the nearest root of the leading coefficient."""
    while True:
        cc = complex(float(c.real), float(c.imag))
        Tc = complex(float(T.real), float(T.imag))
        r = min(abs(cc - q) for q in roots)
        dist = abs(Tc - cc)
        h = acb(r / 2)
        a = step_series(op.local(c, h), op.order, [der[k] * h ** k / math.factorial(k) for k in range(op.order)], 260)
        last = dist <= r / 2 * 0.999
        nxt = T if last else c + (T - c) * (r / 2 / dist)
        ders = eval_series(a, (nxt - c) / h, op.order - 1)
        c, der = nxt, [ders[k] / h ** k for k in range(op.order)]
        if last:
            return c, der


def test_six_coefficients_and_the_operator_fix_the_kite():
    ctx.dps = 80
    try:
        op = Operator(IB["pf_op"])
        ref = [arb(v) for v in IB["pf_series"]]
        a0, freed = op.series_at_zero({n: ref[n] for n in range(6)}, 300)
        assert freed == [0, 1, 2, 3, 4, 5]
        assert max(float(abs(a0[n] / ref[n] - 1)) for n in range(6, 21)) < 1e-50
        roots = _roots(IB["pf_op"])
        assert any(abs(q.imag) < 1e-12 and 0 < q.real < 0.25 for q in roots)   # the apparent root to pass above
        c = acb(arb(1) / 16)
        der = eval_series([acb(x) for x in a0], c, op.order - 1)
        c, der = walk(op, roots, c, der, acb(arb("0.17"), arb("0.06")))
        c, der = walk(op, roots, c, der, acb(arb("0.30"), arb("0.06")))
        for node in NODES:
            z = arb(node["w"]) ** 2
            assert float(z) > 0.25
            _, dz = walk(op, roots, c, der, acb(z))
            assert float(abs(dz[0].real / arb(node["phi"]) - 1)) < 1e-60, node["index"]
    finally:
        ctx.dps = 15
