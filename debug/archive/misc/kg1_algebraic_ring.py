"""KG-1: Algebraic ring of bare Klein-Gordon spectrum on S^3 x R.

Tests whether omega_n^2 = n(n+2)/R^2 + m^2 lives in Q[sqrt(d_1), sqrt(d_2), ...]
for rational m^2 panel, and exhibits ring breakage for irrational m^2.
R = 1 throughout.
"""
from __future__ import annotations

import json
from pathlib import Path

import sympy as sp

OUT = Path(__file__).parent / "data" / "kg1_algebraic_ring.json"


def squarefree_decompose(q: sp.Rational) -> tuple[sp.Rational, int]:
    """Write q = c^2 * d with d square-free positive integer, c rational.

    Works for positive rationals.  Returns (c, d) with sqrt(q) = c * sqrt(d).
    """
    if q == 0:
        return sp.Integer(0), 1
    q = sp.Rational(q)  # ensure exact rational, no nsimplify
    num, den = q.p, q.q
    num = int(num)
    den = int(den)
    # sqrt(num/den) = sqrt(num*den)/den
    radicand = num * den
    # extract largest square divisor
    sq = 1
    d = radicand
    p = 2
    while p * p <= d:
        while d % (p * p) == 0:
            d //= (p * p)
            sq *= p
        p += 1
    c = sp.Rational(sq, den)
    return c, d


def analyze_rational_m2(m2: sp.Rational, n_max: int = 50):
    """For each n, compute omega_n^2 = n(n+2) + m^2 and decompose."""
    generators: dict[int, list[int]] = {}  # n -> sorted list of d
    seen_d: set[int] = set()
    contains_pi = False  # sanity flag (should stay False)
    rows = []
    for n in range(1, n_max + 1):
        omega2 = sp.Rational(n * (n + 2)) + m2
        c, d = squarefree_decompose(omega2)
        omega = c * sp.sqrt(d)
        # verify
        assert sp.simplify(omega**2 - omega2) == 0, (n, omega2, c, d)
        # check transcendental content
        free = omega.free_symbols
        if sp.pi in omega.atoms() or sp.E in omega.atoms() or free:
            contains_pi = True
        if d != 1:
            seen_d.add(d)
        rows.append({
            "n": n,
            "omega2_sym": str(omega2),
            "omega_sym": str(omega),
            "rational_coeff": str(c),
            "squarefree_d": d,
        })
    return {
        "m2": str(m2),
        "n_max": n_max,
        "squarefree_generators": sorted(seen_d),
        "contains_transcendental": contains_pi,
        "rows": rows,
    }


def analyze_irrational_m2(m2_expr: sp.Expr, m2_label: str, n_max: int = 10):
    """Negative-control: irrational m^2 should break the ring at n=1."""
    rows = []
    first_irrational_n = None
    for n in range(1, n_max + 1):
        omega2 = sp.Rational(n * (n + 2)) + m2_expr
        # check whether omega2 is rational
        is_rational = omega2.is_rational
        # check transcendental content
        contains_pi = sp.pi in omega2.atoms()
        contains_sqrt = any(isinstance(a, sp.Pow) and a.exp == sp.Rational(1, 2)
                            for a in sp.preorder_traversal(omega2))
        if (not is_rational) and first_irrational_n is None:
            first_irrational_n = n
        rows.append({
            "n": n,
            "omega2_sym": str(omega2),
            "is_rational": bool(is_rational) if is_rational is not None else None,
            "contains_pi": bool(contains_pi),
            "contains_sqrt": bool(contains_sqrt),
        })
    return {
        "m2_label": m2_label,
        "m2_expr": str(m2_expr),
        "n_max": n_max,
        "first_irrational_n": first_irrational_n,
        "rows": rows,
    }


def main():
    rational_panel = [sp.Integer(0), sp.Integer(1), sp.Rational(1, 4), sp.Integer(2)]
    irrational_panel = [
        (sp.Rational(1, 1) / sp.pi, "1/pi"),
        (sp.sqrt(2), "sqrt(2)"),
    ]

    rational_results = []
    for m2 in rational_panel:
        rational_results.append(analyze_rational_m2(m2, n_max=50))

    irrational_results = []
    for m2_expr, label in irrational_panel:
        irrational_results.append(analyze_irrational_m2(m2_expr, label, n_max=10))

    # Headline
    union_d = sorted({d for r in rational_results for d in r["squarefree_generators"]})
    headline = {
        "union_squarefree_generators_rational_panel": union_d,
        "any_rational_m2_contains_transcendental": any(
            r["contains_transcendental"] for r in rational_results
        ),
        "irrational_m2_first_break_n": {
            r["m2_label"]: r["first_irrational_n"] for r in irrational_results
        },
    }

    out = {
        "script": "kg1_algebraic_ring.py",
        "headline": headline,
        "rational_panel": rational_results,
        "irrational_panel": irrational_results,
    }
    OUT.write_text(json.dumps(out, indent=2))
    print(json.dumps(headline, indent=2))


if __name__ == "__main__":
    main()
