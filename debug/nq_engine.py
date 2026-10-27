"""Engine for docs/predictions/2026-10-27_kite_ninth_point_origin.md (phase 2).
Lee-Pomeransky critical-point ideal of the kite at P2 (or given squared rates), eps = m5^2, modulo a prime p:
  nu_i G - s x_i dG/dx_i (i = 1..5),  t x1...x5 G - 1.
sample(p, zres, eres) -> (count, [tr(t), tr(x1), ..., tr(x5)]): count = number of standard monomials of the grevlex
Groebner basis (critical points with multiplicity), traces of multiplication by each variable on the quotient.
Same nu, s and kite construction as debug/lit_scan_euler_kite.py."""
from __future__ import annotations

import itertools
from fractions import Fraction as Fr

import sympy as sp

NU = [Fr(3, 7), Fr(5, 11), Fr(2, 13), Fr(7, 17), Fr(4, 19)]
S = Fr(9, 23)
P2SQ = (Fr(1), Fr(9, 4), Fr(49, 9), Fr(64, 9))   # (A, B, G1, G2)

_x = sp.symbols("x1:6")
_t = sp.Symbol("t")
_z, _e = sp.symbols("z eps")
_CACHE = {}


def res(q, p: int) -> int:
    q = Fr(q)
    return q.numerator * pow(q.denominator, -1, p) % p


def build_G(sq=P2SQ):
    """G = U + F as {x-monomial: {(zdeg, epsdeg): Fraction}}."""
    key = tuple(sq)
    if key in _CACHE:
        return _CACHE[key]
    x = _x
    U = sp.expand((x[0] + x[1]) * (x[2] + x[3]) + x[4] * (x[0] + x[1] + x[2] + x[3]))
    k1, k2, q1, q2, w = sp.symbols("k1 k2 q1 q2 w")
    mom = {1: (k1, k2), 2: (w - k1, -k2), 3: (k1 - q1, k2 - q2), 4: (w - k1 + q1, -k2 + q2), 5: (q1, q2)}
    m2 = {i + 1: sp.Rational(sq[i].numerator, sq[i].denominator) for i in range(4)}
    m2[5] = _e
    Q = sp.expand(sum(x[i - 1] * (mom[i][0]**2 + mom[i][1]**2 + m2[i]) for i in mom))
    X = [k1, q1]
    H = sp.Matrix(2, 2, lambda i, j: sp.diff(Q, X[i], X[j]) / 2)
    assert sp.expand(H.det() - U) == 0
    gv = sp.Matrix([sp.diff(Q, xx) / 2 for xx in X]).subs({k1: 0, q1: 0})
    Q0 = Q.subs({k1: 0, q1: 0, k2: 0, q2: 0})
    Fpol = sp.expand(sp.expand(U * Q0 - (gv.T * H.adjugate() * gv)[0]).subs(w**2, _z))
    assert w not in Fpol.free_symbols
    G = sp.Poly(sp.expand(U + Fpol), *x, _z, _e)
    out = {}
    for mono, c in G.terms():
        out.setdefault(mono[:5], {})[(mono[5], mono[6])] = Fr(int(sp.numer(c)), int(sp.denom(c)))
    _CACHE[key] = out
    return out


def sample(p: int, zres: int, eres: int, sq=P2SQ, want_traces: bool = True):
    Gd = build_G(sq)
    x, t = _x, _t
    terms = {}
    for m, cd in Gd.items():
        c = sum(res(cf, p) * pow(zres, dz, p) * pow(eres, de, p) for (dz, de), cf in cd.items()) % p
        if c:
            terms[m] = c
    g = sum(sp.Integer(c) * sp.prod([x[j]**m[j] for j in range(5)]) for m, c in terms.items())
    eqs = [sp.expand(res(NU[i], p) * g - res(S, p) * x[i] * sp.diff(g, x[i])) for i in range(5)]
    eqs.append(t * x[0] * x[1] * x[2] * x[3] * x[4] * g - 1)
    gens = (t,) + tuple(x)
    gb = sp.groebner(eqs, *gens, order="grevlex", modulus=p)
    lead = [sp.Poly(e, *gens, modulus=p).monoms(order="grevlex")[0] for e in gb.exprs]
    nv = 6
    for i in range(nv):
        if not any(all(m[j] == 0 for j in range(nv) if j != i) for m in lead):
            return -1, None
    bound = [min(m[i] for m in lead if all(m[j] == 0 for j in range(nv) if j != i)) for i in range(nv)]
    std = [mono for mono in itertools.product(*[range(bb) for bb in bound])
           if not any(all(mono[j] >= m[j] for j in range(nv)) for m in lead)]
    if not want_traces:
        return len(std), None
    traces = []
    for v in gens:
        tr = 0
        for mono in std:
            bexpr = sp.prod([gens[j]**mono[j] for j in range(nv)])
            _, r = gb.reduce(v * bexpr)
            tr = (tr + int(sp.Poly(r, *gens, modulus=p).coeff_monomial(bexpr))) % p
        traces.append(tr)
    return len(std), traces
