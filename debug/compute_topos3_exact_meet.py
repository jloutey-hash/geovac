"""Sprint Topos-3 (exact redo) — two-center frame meet via Mulliken/
Ruedenberg auxiliary integrals.  Replaces the parked quadrature driver.

Every overlap <phi_{n1 l1 m}(0, Z1) | phi_{n2 l2 m}(R zhat, Z2)> reduces in
prolate spheroidal coordinates to a POLYNOMIAL in (xi, eta) times
e^{-p xi} e^{-q eta}, with p = (a+b)R/2, q = (a-b)R/2 (a = Z1/n1,
b = Z2/n2).  The classical cancellation is ASSERTED (denominator == 1
after sympy cancel) — a built-in correctness check.  With the exact
auxiliary integrals
    A_i(p) = int_1^inf xi^i e^{-p xi} dxi   = e^{-p} * a_i,  a_i rational,
    B_j(q) = int_{-1}^1 eta^j e^{-q eta} deta
           = e^{q} u_j + e^{-q} v_j  (q != 0;  u, v rational)
           = 2/(j+1) (j even) or 0 (j odd)   (q == 0),
the unnormalized overlap is
    S = e^{-p} ( U e^{q} + V e^{-q} )      (U, V exact rationals),
so S = 0  <=>  U = V = 0  (Lindemann independence of e^{q}, e^{-q} over Q
for rational q != 0)  — support is EXACTLY decidable.

Verdicts computed:
  * meet: m-conservation structural; per-m-block support connectivity
    (exact) => is the two-center meet exactly the m-grading algebra?
    Ladder: identical N -> same-center sum(2l+1) -> two-center #m.
  * join: connected support + Burnside => full M_d per m-block; d >= 3
    blocks => KS-obstructed (Topos-1 witness).
  * cross-check: normalized <1s|1s>(R) must equal (1 + R + R^2/3) e^{-R}
    EXACTLY (rational identity at rational R).

Output: debug/data/sprint_topos3_exact_meet.json
Backing test: tests/test_topos3_exact_meet.py (self-contained).
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import comb, factorial

import sympy as sp

xi, eta = sp.symbols("xi eta", positive=False)


# ------------------------------------------------------------ radial data

def radial_poly_coeffs(Z: Fraction, n: int, l: int):
    """Unnormalized hydrogenic radial:  R(r) = sum_k c_k r^k * e^{-Z r/n},
    c_k exact rationals (rho^l L^{2l+1}_{n-l-1}(rho), rho = 2Z r/n)."""
    lam = sp.Rational(2 * Z.numerator, Z.denominator * n)   # 2Z/n
    lag = [sp.Rational((-1) ** j * comb(n - l - 1 + 2 * l + 1, n - l - 1 - j),
                       factorial(j)) for j in range(n - l)]
    return {l + j: c * lam ** j for j, c in enumerate(lag)}, \
        sp.Rational(Z.numerator, Z.denominator * n)          # decay a = Z/n


def radial_norm_sq(Z: Fraction, n: int, l: int) -> sp.Rational:
    """Exact  int_0^inf R(r)^2 r^2 dr  for the unnormalized radial above."""
    coeffs, a = radial_poly_coeffs(Z, n, l)
    total = sp.Integer(0)
    for k1, c1 in coeffs.items():
        for k2, c2 in coeffs.items():
            k = k1 + k2 + 2
            total += c1 * c2 * sp.factorial(k) / (2 * a) ** (k + 1)
    return sp.nsimplify(total)


def theta_norm_sq(l: int, m: int) -> sp.Rational:
    """int_{-1}^1 P_l^m(x)^2 dx = 2 (l+m)! / ((2l+1)(l-m)!)."""
    return sp.Rational(2 * factorial(l + m), (2 * l + 1) * factorial(l - m))


def G_lm(l: int, m: int, x):
    """Polynomial part of the associated Legendre:  d^m/dx^m P_l(x)."""
    t = sp.Symbol("_t")
    return sp.diff(sp.legendre(l, t), t, m).subs(t, x)


# ------------------------------------------------------ integrand assembly

def overlap_UV(Z1: Fraction, n1: int, l1: int,
               Z2: Fraction, n2: int, l2: int,
               m: int, R: Fraction):
    """Exact unnormalized overlap:  returns (p, q, U, V) with
    S = e^{-p} (U e^{q} + V e^{-q})   [q != 0]   or
    S = e^{-p} * U, V = None          [q == 0]."""
    m = abs(m)
    Rs = sp.Rational(R.numerator, R.denominator)
    half = Rs / 2
    c1, a = radial_poly_coeffs(Z1, n1, l1)
    c2, b = radial_poly_coeffs(Z2, n2, l2)
    ra = half * (xi + eta)
    rb = half * (xi - eta)
    ct_a = (1 + xi * eta) / (xi + eta)
    ct_b = (xi * eta - 1) / (xi - eta)
    Rad1 = sum(c * ra ** k for k, c in c1.items())
    Rad2 = sum(c * rb ** k for k, c in c2.items())
    # theta parts: P_l^m(x) = (1-x^2)^{m/2} G_lm(x)  (sign support-neutral);
    # (1-ct_a^2)(1-ct_b^2) = (xi^2-1)^2 (1-eta^2)^2 / ((xi+eta)(xi-eta))^2
    s2 = (xi ** 2 - 1) * (1 - eta ** 2)
    ang = (s2 ** m / ((xi + eta) * (xi - eta)) ** m
           * G_lm(l1, m, ct_a) * G_lm(l2, m, ct_b))
    vol = half ** 3 * (xi ** 2 - eta ** 2)
    expr = sp.cancel(sp.together(Rad1 * Rad2 * ang * vol))
    num, den = sp.fraction(expr)
    assert den == 1, f"classical cancellation FAILED: den = {den}"
    poly = sp.Poly(sp.expand(num), xi, eta)
    p = (a + b) * Rs / 2
    q = (a - b) * Rs / 2

    # exact auxiliary integrals
    max_i = max(mon[0] for mon in poly.monoms())
    max_j = max(mon[1] for mon in poly.monoms())
    a_coef = {}
    a_coef[0] = 1 / p
    for i in range(1, max_i + 1):
        a_coef[i] = (1 + i * a_coef[i - 1]) / p
    if q != 0:
        u = {0: 1 / q}
        v = {0: -1 / q}
        for j in range(1, max_j + 1):
            u[j] = sp.Rational((-1) ** j, 1) / q + j * u[j - 1] / q
            v[j] = sp.Rational(-1, 1) / q + j * v[j - 1] / q
        U = sp.Integer(0)
        V = sp.Integer(0)
        for (i, j), c in zip(poly.monoms(), poly.coeffs()):
            U += c * a_coef[i] * u[j]
            V += c * a_coef[i] * v[j]
        return sp.nsimplify(p), sp.nsimplify(q), sp.nsimplify(U), sp.nsimplify(V)
    else:
        U = sp.Integer(0)
        for (i, j), c in zip(poly.monoms(), poly.coeffs()):
            if j % 2 == 0:
                U += c * a_coef[i] * sp.Rational(2, j + 1)
        return sp.nsimplify(p), sp.Integer(0), sp.nsimplify(U), None


def is_zero(p, q, U, V) -> bool:
    return (U == 0) if V is None else (U == 0 and V == 0)


# ---------------------------------------------------------------- probe

def m_blocks(n_max):
    out = {}
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                out.setdefault(m, []).append((n, l))
    return out


def components(supp):
    d = len(supp)
    parent = list(range(2 * d))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(d):
        for j in range(d):
            if supp[i][j]:
                ra, rb = find(i), find(d + j)
                if ra != rb:
                    parent[ra] = rb
    return len({find(x) for x in range(2 * d)})


def main():
    results = {}

    # exact cross-check: normalized <1s(Z=1)|1s(Z=1)>(R) = (1+R+R^2/3) e^{-R}
    for R in (Fraction(1), Fraction(2), Fraction(7, 2)):
        p, q, U, V = overlap_UV(Fraction(1), 1, 0, Fraction(1), 1, 0, 0, R)
        norm_sq = (radial_norm_sq(Fraction(1), 1, 0) * theta_norm_sq(0, 0))
        S_norm = sp.nsimplify(U / norm_sq)     # q == 0 here; e^{-p} = e^{-R}
        Rq = sp.Rational(R.numerator, R.denominator)
        expected = sp.nsimplify(1 + Rq + Rq ** 2 / 3)
        results.setdefault("crosscheck_1s1s", {})[f"R={R}"] = {
            "computed_rational": str(S_norm), "expected": str(expected),
            "exact_match": bool(sp.simplify(S_norm - expected) == 0)}

    n_max = 3
    blocks = m_blocks(n_max)
    N = sum(len(v) for v in blocks.values())
    same_center = sum(2 * l + 1 for l in range(n_max))
    for (Z2, R) in ((Fraction(1), Fraction(1)), (Fraction(1), Fraction(2)),
                    (Fraction(3), Fraction(1))):
        key = f"Z1=1,Z2={Z2},R={R},n_max={n_max}"
        cell = {"m_blocks": {}}
        all_conn = True
        for m, states in sorted(blocks.items()):
            if m < 0:
                continue
            d = len(states)
            supp = []
            zeros = []
            for aa in states:
                row = []
                for bb in states:
                    puv = overlap_UV(Fraction(1), aa[0], aa[1],
                                     Z2, bb[0], bb[1], m, R)
                    z = is_zero(*puv)
                    row.append(0 if z else 1)
                    if z:
                        zeros.append((aa, bb))
                supp.append(row)
            conn = components(supp) == 1
            all_conn &= conn
            cell["m_blocks"][f"m={m}"] = {
                "states": states, "dim": d, "connected": conn,
                "full_support": all(all(r) for r in supp),
                "exact_zeros": zeros}
        cell["meet_is_m_grading"] = all_conn
        cell["ladder"] = {"identical": N, "same_center_(l,m)": same_center,
                          "two_center_(m)": len(blocks)}
        cell["join_ks_obstructed_dims"] = sorted(
            {len(v) for v in blocks.values() if len(v) >= 3}, reverse=True)
        results[key] = cell

    with open("debug/data/sprint_topos3_exact_meet.json", "w") as fh:
        json.dump(results, fh, indent=1, default=str)

    print("crosscheck 1s1s:", results["crosscheck_1s1s"])
    for key, cell in results.items():
        if not key.startswith("Z1"):
            continue
        print(f"=== {key} ===  meet=m-grading: {cell['meet_is_m_grading']}  "
              f"ladder {cell['ladder']}  KS dims {cell['join_ks_obstructed_dims']}")
        for mk, b in cell["m_blocks"].items():
            print(f"    {mk}: d={b['dim']} connected={b['connected']} "
                  f"full={b['full_support']} zeros={b['exact_zeros']}")


if __name__ == "__main__":
    main()
