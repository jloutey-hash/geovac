"""Backing test for Paper 55 eq:geovac_grothendieck_class_explicit (M2 motive).

Promoted from debug/sprint_a8_grothendieck_class_memo.md (June 2026) so the
numbered general-n Grothendieck-class equation carries a tracked, regression-
protected backing instead of a transient debug/ citation (Clean-Room policy,
group3 re-cert 2026-08-24).

What is verified (genuinely, not tautologically):
  The paper states TWO equations --- the *defining* class eq:geovac_grothendieck_class
  in terms of the quadric class [Z_{1,2n}], and the *explicit* closed form
  eq:geovac_grothendieck_class_explicit.  This test DERIVES the explicit form
  from the defining equation + the Fathizadeh-Marcolli quadric class
  eq:z_quadric_class and asserts symbolic equality (sympy expand), for general
  n = 1..8 and via the two independent forms of [Z_{1,2n}] (palindrome sum vs
  [P^n]*(1+L^n)).  A wrong coefficient anywhere in the explicit form fails.
  L is the Lefschetz motive (a free symbol); all classes must land in Z[L].
"""
from __future__ import annotations

import sympy as sp
import pytest

L = sp.Symbol("L")


def z_quadric_palindrome(n: int) -> sp.Expr:
    """[Z_{1,2n}] as the length-(2n+1) palindrome with doubled middle L^n.

    eq:z_quadric_class LHS: L^{2n}+...+L^{n+1} + 2 L^n + L^{n-1}+...+L+1.
    """
    return sum(L**k for k in range(0, 2 * n + 1)) + L**n


def z_quadric_projective(n: int) -> sp.Expr:
    """[Z_{1,2n}] = [P^n] * (1 + L^n), the F-M factored form."""
    P_n = sum(L**k for k in range(0, n + 1))  # [P^n] = 1 + L + ... + L^n
    return P_n * (1 + L**n)


def geovac_class_defining(n: int) -> sp.Expr:
    """eq:geovac_grothendieck_class: [V_n] in terms of the quadric class [Z_{1,2n}]."""
    Z = z_quadric_palindrome(n)
    return L ** (2 * n + 3) - 2 * L ** (2 * n + 2) - (L - 2) * (L - 1) * Z - (L - 2)


def geovac_class_explicit(n: int) -> sp.Expr:
    """eq:geovac_grothendieck_class_explicit: the fully-expanded closed form."""
    return (L ** (2 * n + 3) - 3 * L ** (2 * n + 2) + 2 * L ** (2 * n + 1)
            - L ** (n + 2) + 3 * L ** (n + 1) - 2 * L ** n)


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6, 7, 8])
def test_two_quadric_forms_agree(n):
    """[Z_{1,2n}] palindrome == [P^n]*(1+L^n) (eq:z_quadric_class RHS)."""
    assert sp.expand(z_quadric_palindrome(n) - z_quadric_projective(n)) == 0


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5, 6, 7, 8])
def test_explicit_class_derives_from_defining(n):
    """The explicit closed form FOLLOWS from the defining equation + quadric class.

    This is the load-bearing check: it re-derives eq:...explicit from
    eq:geovac_grothendieck_class rather than restating it.
    """
    assert sp.expand(geovac_class_defining(n) - geovac_class_explicit(n)) == 0


def test_n1_collapse():
    """At n=1 the two middle terms collapse (2n+1 = n+2 = 3): L^5-3L^4+L^3+3L^2-2L."""
    expected = L**5 - 3 * L**4 + L**3 + 3 * L**2 - 2 * L
    assert sp.expand(geovac_class_explicit(1) - expected) == 0
    # and it agrees with the defining equation at n=1
    assert sp.expand(geovac_class_defining(1) - expected) == 0


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_class_is_in_Z_of_L(n):
    """[V_n^GeoVac] in Z[L] (integer coefficients) -- the mixed-Tate sub-ring claim."""
    poly = sp.Poly(sp.expand(geovac_class_explicit(n)), L)
    for c in poly.all_coeffs():
        assert c == int(c), f"non-integer coefficient {c} at n={n}"


def test_symbolic_general_n_identity():
    """General-n symbolic identity (n a free positive integer symbol).

    Uses the two closed geometric-series forms so the exponents are symbolic,
    proving the identity for ALL n at once, not only the sampled values.
    """
    n = sp.Symbol("n", positive=True, integer=True)
    # [Z_{1,2n}] via geometric series: (L^{2n+1}-1)/(L-1) + L^n
    Z = (L ** (2 * n + 1) - 1) / (L - 1) + L**n
    defining = (L ** (2 * n + 3) - 2 * L ** (2 * n + 2)
                - (L - 2) * (L - 1) * Z - (L - 2))
    explicit = (L ** (2 * n + 3) - 3 * L ** (2 * n + 2) + 2 * L ** (2 * n + 1)
                - L ** (n + 2) + 3 * L ** (n + 1) - 2 * L ** n)
    assert sp.simplify(sp.expand(defining) - explicit) == 0


def test_nontautology_guard_wrong_coefficient_fails():
    """A perturbed explicit form must NOT match the defining equation.

    Guards against the test degenerating into a self-comparison: flip the
    -3 L^{2n+2} coefficient to -4 and confirm the derivation-consistency
    check would FAIL for it.
    """
    n = 3
    wrong = geovac_class_explicit(n) - L ** (2 * n + 2)  # -3 -> -4
    assert sp.expand(geovac_class_defining(n) - wrong) != 0
