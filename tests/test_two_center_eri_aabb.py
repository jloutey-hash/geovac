"""Regression net for the two-center ERI build, (AA|BB) class.

Pins everything increment 1 and 1b validated, so that increment 1c -- replacing
the 2D quadrature with a closed form -- cannot silently drift. Without this the
closed form would be checked only against whatever I happened to remember.

Backs `geovac/two_center_eri.py` (promoted from an exploratory driver once it was
validated four independent ways).

Coverage, weakest to strongest:
  radial integral unit tests      lower/upper_integral vs numerical quadrature,
                                  including the negative-exponent branch where
                                  E_1 enters
  multipole structure             Phase 0' predictions: L range, parity, M rule
  seed structure                  E_1 absent at L=0, present at L=2 -- pinned
                                  because 1c must not lose the seed accounting
  monopole sum rule               int rho d3r = <chi_a|chi_b> exactly
  pointwise reconstruction        the (L,M) sum reproduces the orbital product
  EXACT J(R)                      the headline gate: the full (AA|BB) chain
                                  against the closed-form VB J integral
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

sp = pytest.importorskip("sympy")
pytest.importorskip("scipy")

from geovac.two_center_eri import (  # noqa: E402
    V_L_radial, lower_integral, multipole_decomposition, radial_norm,
    radial_poly, upper_integral, r_s,
)

Z1 = Fraction(1)


# --------------------------------------------------------------- unit: radials

@pytest.mark.parametrize("q", [0, 1, 3, 5])
def test_lower_integral_matches_quadrature(q):
    from scipy import integrate
    b, r = 1.3, 2.1
    got = float(lower_integral(q, sp.Float(b), sp.Float(r)))
    ref, _ = integrate.quad(lambda s: s ** q * np.exp(-b * s), 0.0, r,
                            epsabs=1e-13, epsrel=1e-13)
    assert abs(got - ref) < 1e-10, f"q={q}: {got} vs {ref}"


@pytest.mark.parametrize("p", [0, 2, 4, -1, -2, -3])
def test_upper_integral_matches_quadrature(p):
    """Includes the p <= -1 branch, which is where E_1 enters."""
    from scipy import integrate
    b, r = 1.3, 2.1
    got = complex(sp.N(upper_integral(p, sp.Float(b), sp.Float(r)))).real
    ref, _ = integrate.quad(lambda s: s ** p * np.exp(-b * s), r, np.inf,
                            epsabs=1e-13, epsrel=1e-13)
    assert abs(got - ref) < 1e-10, f"p={p}: {got} vs {ref}"


def test_upper_integral_seed_is_E1_only_and_only_below_zero():
    """The seed appears exactly once, at p = -1, and never for p >= 0.

    Pinned because increment 1c's closed form must preserve this accounting --
    Phase 0 Q2 concluded the seed set is {e^a E_1(a)} and this is where that
    becomes operational.
    """
    b, r = sp.Symbol("b", positive=True), sp.Symbol("r", positive=True)
    for p in (0, 1, 4):
        assert not upper_integral(p, b, r).atoms(sp.expint), \
            f"p={p} should be E_1-free"
    for p in (-1, -2, -3):
        assert upper_integral(p, b, r).atoms(sp.expint), \
            f"p={p} should carry E_1"
    # and no OTHER transcendental sneaks in
    allowed = {"exp", "expint"}
    for p in (-1, -3, 2):
        names = {type(f).__name__ for f in upper_integral(p, b, r).atoms(sp.Function)}
        assert names <= allowed, f"p={p} introduced {names - allowed}"


def test_V_L_has_no_exp_polar_artifact():
    """sympy's `integrate` emitted exp_polar here; the hand-written form must not.

    exp_polar is branch bookkeeping that would propagate into 1c's closed form
    and make simplification unreliable.
    """
    rad, b = {0: sp.Integer(1)}, sp.Integer(2)
    for L in (0, 1, 2, 3):
        expr = V_L_radial(rad, b, L)
        assert not expr.atoms(sp.exp_polar), f"L={L} still emits exp_polar"


# ----------------------------------------------------- multipole structure

@pytest.mark.parametrize("l1,m1,l2,m2", [
    (0, 0, 0, 0), (1, 1, 1, 1), (1, 1, 0, 0), (2, 2, 1, 1), (2, 1, 2, -1),
])
def test_multipole_structure_matches_phase0p(l1, m1, l2, m2):
    """L in {|l1-l2|..l1+l2}, parity l1+l2+L even, M = m2-m1 only."""
    terms = multipole_decomposition(Z1, max(l1 + 1, 1), l1, m1,
                                    Z1, max(l2 + 1, 1), l2, m2)
    Ls = sorted({L for L, _M, _g, _r, _b in terms})
    Ms = {M for _L, M, _g, _r, _b in terms}
    assert Ls, "no surviving multipole"
    assert max(Ls) <= l1 + l2, f"L exceeded l1+l2: {Ls}"
    assert all((l1 + l2 + L) % 2 == 0 for L in Ls), f"parity violated: {Ls}"
    assert Ms == {m2 - m1}, f"M rule violated: {Ms}"


def test_monopole_sum_rule():
    """int rho d3r = <chi_a|chi_b>: 1 on the diagonal, 0 for orthogonal pairs."""
    for (n1, l1), (n2, l2), want in ((((1, 0)), ((1, 0)), 1),
                                     (((2, 0)), ((2, 0)), 1),
                                     (((1, 0)), ((2, 0)), 0),
                                     (((2, 1)), ((2, 1)), 1)):
        terms = multipole_decomposition(Z1, n1, l1, 0, Z1, n2, l2, 0)
        tot = 0
        for L, M, g, rad, b in terms:
            if (L, M) != (0, 0):
                continue
            tot += g * sp.sqrt(4 * sp.pi) * sum(
                c * sp.factorial(k + 2) / b ** (k + 3) for k, c in rad.items())
        assert abs(float(sp.nsimplify(tot)) - want) < 1e-10


def test_pointwise_reconstruction():
    """The (L,M) sum reproduces conj(chi_a)(r) chi_b(r) -- the strong leg.

    Would catch any Gaunt phase or normalization slip.
    """
    def chi(Z, n, l, m, rv, th, ph):
        c, a = radial_poly(Z, n, l)
        N = radial_norm(Z, n, l)
        rad = sum(float(N * cc) * rv ** k for k, cc in c.items()) \
            * np.exp(-float(a) * rv)
        return rad * complex(sp.Ynm(l, m, th, ph).expand(func=True).evalf())

    worst = 0.0
    for (n1, l1, m1), (n2, l2, m2) in (((2, 1, 1), (2, 1, 1)),
                                       ((2, 1, 1), (1, 0, 0)),
                                       ((2, 1, 0), (2, 1, 0))):
        terms = multipole_decomposition(Z1, n1, l1, m1, Z1, n2, l2, m2)
        for rv, th, ph in ((0.8, 0.7, 0.4), (1.9, 2.1, 5.0)):
            direct = np.conj(chi(Z1, n1, l1, m1, rv, th, ph)) \
                * chi(Z1, n2, l2, m2, rv, th, ph)
            recon = 0j
            for L, M, g, rad, b in terms:
                radv = sum(float(c) * rv ** k for k, c in rad.items()) \
                    * np.exp(-float(b) * rv)
                recon += complex(g) * radv * complex(
                    sp.Ynm(L, M, th, ph).expand(func=True).evalf())
            worst = max(worst, abs(direct - recon))
    assert worst < 1e-9, f"reconstruction worst error {worst:.3e}"


# ------------------------------------------------------- the headline gate

def _J_exact(zeta: float, R: float) -> float:
    rho = zeta * R
    return (1.0 / R) * (1.0 - np.exp(-2 * rho)
                        * (1 + 1.375 * rho + 0.75 * rho ** 2 + rho ** 3 / 6.0))


@pytest.mark.slow
@pytest.mark.parametrize("R", [1.5, 2.5, 4.0])
def test_aabb_reproduces_exact_J_integral(R):
    """(1s_A 1s_A | 1s_B 1s_B) against the closed-form VB J integral.

    THE gate for the (AA|BB) class. Increment 1c's closed form must reproduce
    this; it is an independent analytic result, not a self-consistency check.
    """
    from scipy import integrate

    terms = multipole_decomposition(Z1, 1, 0, 0, Z1, 1, 0, 0)
    total = 0.0
    for LA, MA, gA, radA, bA in terms:
        VA = sp.lambdify(r_s, sp.re(V_L_radial(radA, bA, LA)), "numpy")
        for LB, MB, gB, radB, bB in terms:
            if MA + MB != 0:
                continue
            radB_f = sp.lambdify(
                r_s,
                sum(c * r_s ** k for k, c in radB.items()) * sp.exp(-bB * r_s),
                "numpy")
            nA = float(sp.sqrt((2 * LA + 1) / (4 * sp.pi)))
            nB = float(sp.sqrt((2 * LB + 1) / (4 * sp.pi)))

            def integrand(u, rb):
                rA = np.sqrt(rb * rb + R * R + 2 * rb * R * u)
                if rA < 1e-12:
                    return 0.0
                return radB_f(rb) * rb * rb * nB * VA(rA) * nA

            val, _ = integrate.dblquad(integrand, 0.0, 60.0,
                                       lambda _r: -1.0, lambda _r: 1.0,
                                       epsabs=1e-11, epsrel=1e-11)
            total += float(gA) * float(gB) * 2 * np.pi * val

    ref = _J_exact(1.0, R)
    assert abs(total - ref) < 1e-9, (
        f"R={R}: chain gives {total:.12f}, exact J(R) = {ref:.12f}"
    )
