"""Minimality audit, step 1: irreducibility of the Paper 59 operators by local exponents + Fuchs.

For a Fuchsian operator L of order n on P^1 with singular set S, a right factor B of order k
(1 <= k < n) has, at every p in S, local exponents forming a sub-multiset of L's exponents at p
(the valuation multiset of a subspace of the local solution space), and at each extra (apparent)
singularity an exponent sum >= k(k-1)/2 + 1.  Fuchs' relation for B,
    sum_{p in P^1} (sigma_p(B) - k(k-1)/2) = -k(k-1),
then forces  sum_{p in S} minsub_k(p) <= (|S| - 2) k (k-1)/2,
where minsub_k(p) is the smallest sum of k of L's exponents at p (real parts).  If the inequality
fails for every k, L has no right factor of any order and is irreducible (no degree bound).

Operators are given as sum_i t^(2i) Q_i(theta), theta = t d/dt (function-left), as in Paper 59.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import sympy as sp

t, lam, u, s = sp.symbols("t lambda u s")
OUT = Path(__file__).parent / "data" / "minimality_audit_fuchs.json"


def theta_to_D(Q):
    """sum_i t^(2i) Q_i(theta) -> list a_j(t) with L = sum_j a_j(t) d^j/dt^j."""
    n = max(sp.Poly(q, u).degree() for q in Q)
    a = [sp.Integer(0)] * (n + 1)
    for i, q in enumerate(Q):
        for (k,), c in sp.Poly(sp.expand(q), u).terms():
            # theta^k = sum_j S(k,j) t^j D^j
            for j in range(k + 1):
                a[j] += c * t ** (2 * i) * sp.functions.combinatorial.numbers.stirling(k, j) * t**j
    return [sp.expand(x) for x in a]


def indicial(a, p):
    """Indicial polynomial of sum a_j D^j at t = p (p finite or sp.oo); exponent = valuation."""
    if p is sp.oo:
        # y = t^(-lam): D^j t^(-lam) = ff(-lam, j) t^(-lam-j)
        expr = sum(aj * sp.ff(-lam, j) * t ** (-j) for j, aj in enumerate(a))
        expr = sp.expand(expr)
        # leading behaviour as t -> oo: highest power of t
        terms = sp.Poly(sp.expand(expr * t ** 50), t)
        top = terms.degree()
        return sp.factor(terms.coeff_monomial(t**top))
    x = sp.Symbol("x")
    expr = sum(aj.subs(t, x + p) * sp.ff(lam, j) * x ** (-j) for j, aj in enumerate(a))
    num = sp.Poly(sp.expand(expr * x ** 50), x)
    low = min(m[0] for m in num.monoms())
    return sp.factor(num.coeff_monomial(x**low))


def roots_multiset(poly):
    r = sp.roots(sp.Poly(poly, lam))
    out = []
    for k, m in r.items():
        out += [k] * m
    return sorted(out, key=lambda z: float(sp.re(z)))


def analyse(name, Q):
    a = theta_to_D(Q)
    n = len(a) - 1
    lead = sp.factor(a[-1])
    fin = [r for r in sp.roots(sp.Poly(sp.cancel(a[-1] / t**n), t)).keys()]
    pts = sorted(set([sp.Integer(0)] + fin), key=lambda z: (float(sp.re(z)), float(sp.im(z)))) + [sp.oo]
    rec = {"order": n, "leading_coeff": str(lead), "points": {}}
    ok = True
    exps = {}
    for p in pts:
        I = indicial(a, p)
        e = roots_multiset(I)
        deg = sp.Poly(I, lam).degree()
        if deg != n:
            ok = False  # not Fuchsian at p, or p not singular in the expected way
        exps[p] = e
        rec["points"][str(p)] = {"indicial": str(I), "exponents": [str(x) for x in e]}
    # drop points that are actually ordinary (exponents 0..n-1, no log check needed for bound)
    S = [p for p in pts if sorted(exps[p]) != list(range(n))]
    fuchs_total = sum(sum(exps[p]) for p in pts)
    rec["fuchs_check"] = {"sum_exponents": str(sp.nsimplify(fuchs_total)),
                          "expected": str((len(pts) - 2) * n * (n - 1) // 2)}
    rec["fuchsian_everywhere"] = ok
    verdict = {}
    for k in range(1, n):
        mins = []
        for p in pts:
            e = [float(sp.re(x)) for x in exps[p]]
            mins.append(min(sum(c) for c in itertools.combinations(e, k)))
        lhs = sum(mins)
        rhs = (len(pts) - 2) * k * (k - 1) / 2
        verdict[k] = {"min_sum": lhs, "bound": rhs, "factor_possible": lhs <= rhs + 1e-12}
    rec["right_factor_orders"] = verdict
    rec["irreducible_by_fuchs"] = ok and not any(v["factor_possible"] for v in verdict.values())
    print(f"--- {name}: order {n}, leading coeff {lead}")
    for p in pts:
        print(f"   t={p}: exponents {exps[p]}")
    print(f"   Fuchs: sum {rec['fuchs_check']}")
    for k, v in verdict.items():
        print(f"   right factor of order {k}: min exponent sum {v['min_sum']} vs bound {v['bound']}"
              f" -> {'POSSIBLE' if v['factor_possible'] else 'excluded'}")
    print(f"   irreducible (unbounded degree): {rec['irreducible_by_fuchs']}")
    return rec


OPS = {
    "F2_legendre": [-u**2, (u + 1) ** 2],
    "F3": [4 * u**3, -(u + 1) ** 3],
    "F4_aesz58_t": [9 * u**2 * (u - 2) ** 2, -2 * u**2 * (5 * u**2 + 1), (u + 1) ** 4],
    "F5": [-64 * u**2 * (u - 2) ** 2 * (u - 1), 4 * u**3 * (5 * u**2 + 3), -(u + 1) ** 5],
    "ring_LG": [-4 * u * (u - 3) ** 2, 5 * u**3 - 7 * u**2 - 5 * u + 3, -u * (u + 2) * (u + 3)],
}

if __name__ == "__main__":
    out = {name: analyse(name, Q) for name, Q in OPS.items()}
    OUT.write_text(json.dumps(out, indent=1))
