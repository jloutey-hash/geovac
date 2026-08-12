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
    R_s, V_L_radial, aabb_closed_form, aabb_quadrature, lower_integral,
    multipole_decomposition, radial_norm, radial_poly, radial_product,
    shell_kernel_antiderivatives, shell_kernel_quadrature, t_s,
    upper_integral, r_s, x_s, y_s,
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

    This pins the BEHAVIOUR OF `upper_integral` itself, which is correct and
    worth keeping. It says nothing about whether the (AA|BB) class reaches the
    p < 0 branch -- increment 1c showed it never does; see
    `test_no_E1_on_the_physical_multipole_support` below.
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


def test_no_E1_on_the_physical_multipole_support():
    """E_1 cannot enter the (AA|BB) class -- the corrected seed accounting.

    `V_L_radial` reaches its E_1 branch only when the upper-region exponent
    k + 1 - L goes negative. A real orbital product has k >= l1 + l2 while Gaunt
    caps L at l1 + l2, so k + 1 - L >= 1 always.

    Increment 1 reported the opposite ("Ei appears at L = 2") by evaluating
    V_L_radial(rad_1sx1s, L=2) -- but a 1s x 1s product has no L = 2 multipole,
    so that term does not exist. This test pins the question on the support the
    class actually has, and is the reason 1c's output is elementary.
    """
    worst_margin = 10 ** 9
    for n1 in range(1, 5):
        for l1 in range(n1):
            for n2 in range(1, 5):
                for l2 in range(n2):
                    rad, _b = radial_product(Z1, n1, l1, Z1, n2, l2)
                    k_min = min(rad)
                    for L, _M, _g, r, b in multipole_decomposition(
                            Z1, n1, l1, 0, Z1, n2, l2, 0):
                        worst_margin = min(worst_margin, k_min + 1 - L)
                        assert not V_L_radial(r, b, L).atoms(sp.expint), (
                            f"E_1 in V_L for ({n1},{l1})x({n2},{l2}) at L={L}")
    assert worst_margin >= 1, f"upper-region exponent reached {worst_margin}"


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


# ------------------------------------------- increment 1c: the closed form

Z2 = Fraction(2)


@pytest.mark.parametrize("LA,MA,LB,MB", [
    (0, 0, 0, 0), (1, 0, 1, 0), (2, 0, 0, 0),
    (1, 1, 1, -1), (2, 1, 1, -1), (2, 2, 2, -2), (1, -1, 2, 1),
])
def test_shell_kernel_matches_angular_quadrature(LA, MA, LB, MB):
    """The kernel vs direct angular quadrature -- the convention check.

    Catches Condon-Shortley sign errors, signed-M normalization slips, and a
    wrong pairing of the two sin^|M| factors. Shells chosen NOT to intersect so
    1/d is smooth and the quadrature is an honest reference.
    """
    A_in, A_out, pref = shell_kernel_antiderivatives(LA, MA, LB, MB)
    for xv, yv, Rv in ((0.5, 1.3, 3.0), (3.4, 0.6, 1.5), (0.9, 1.2, 4.0)):
        lo, hi = abs(yv - Rv), yv + Rv
        A = A_out if xv <= lo else A_in
        core = A.subs(t_s, hi) - A.subs(t_s, lo)
        got = float(sp.re((pref / (y_s * R_s) * core).subs(
            {x_s: xv, y_s: yv, R_s: Rv}).evalf()))
        ref = shell_kernel_quadrature(LA, MA, LB, MB, xv, yv, Rv)
        assert abs(got - ref) < 1e-9, f"{got} vs {ref} at x={xv} y={yv} R={Rv}"


def test_closed_form_is_symbolically_the_textbook_J_integral():
    """THE headline gate for 1c: not just numerically equal to J(R) -- equal to
    it as an expression, so the classic result is reproduced by derivation.
    """
    expr = aabb_closed_form(Z1, (1, 0, 0), (1, 0, 0), Z1, (1, 0, 0), (1, 0, 0))
    target = (1 / R_s) * (1 - sp.exp(-2 * R_s)
                          * (1 + sp.Rational(11, 8) * R_s
                             + sp.Rational(3, 4) * R_s ** 2 + R_s ** 3 / 6))
    assert sp.simplify(expr - target) == 0


def test_closed_form_is_elementary():
    """`exp` and nothing else -- no expint, no log, at any l or M.

    This is the structural claim of 1c. It holds because the shell-kernel
    formulation never produces a 1/r_A term (all r_A powers are even on the
    outside branch, odd and >= 1 on the inside branch), and because the
    (y +- R) denominators cancel identically -- both of which are asserted
    inside the builder.
    """
    for args in ((Z1, (1, 0, 0), (1, 0, 0), Z1, (1, 0, 0), (1, 0, 0)),
                 (Z1, (2, 1, 0), (2, 1, 0), Z2, (1, 0, 0), (1, 0, 0)),
                 (Z1, (3, 2, 1), (3, 2, 0), Z2, (2, 1, 0), (2, 1, 1))):
        e = aabb_closed_form(*args)
        assert not e.atoms(sp.expint), "E_1 in the closed form"
        assert not e.atoms(sp.log), "log in the closed form"
        assert {type(f).__name__ for f in e.atoms(sp.Function)} <= {"exp"}


@pytest.mark.parametrize("name,ZA,oa,ob,ZB,oc,od", [
    ("2p0 2p0|1s 1s",    Z1, (2, 1, 0), (2, 1, 0), Z2, (1, 0, 0), (1, 0, 0)),
    ("1s 2p0|2s 2p0",    Z1, (1, 0, 0), (2, 1, 0), Z2, (2, 0, 0), (2, 1, 0)),
    ("2p1 2p0|2p0 2p1",  Z1, (2, 1, 1), (2, 1, 0), Z2, (2, 1, 0), (2, 1, 1)),
    ("3d2 3d0|2p-1 2p1", Z1, (3, 2, 2), (3, 2, 0), Z2, (2, 1, -1), (2, 1, 1)),
])
def test_closed_form_matches_quadrature_for_l_gt_0_and_M_ne_0(
        name, ZA, oa, ob, ZB, oc, od):
    """vs the 1b quadrature route, which shares no 1c code.

    The only reference that can referee M != 0: a Cartesian-Gaussian engine
    cannot express a single complex Y_lm with m != 0.
    """
    R = 2.5
    got = float(sp.re(sp.N(aabb_closed_form(ZA, oa, ob, ZB, oc, od,
                                            sp.nsimplify(R)), 30)))
    ref = aabb_quadrature(ZA, oa, ob, ZB, oc, od, R)
    assert abs(got - ref) < 1e-11, f"{name}: {got} vs {ref}"


@pytest.mark.parametrize("ZA,oa,ob,ZB,oc,od", [
    (Z1, (1, 0, 0), (1, 0, 0), Z2, (2, 1, 0), (2, 1, 0)),
    (Z1, (2, 1, 0), (2, 0, 0), Z2, (2, 1, 0), (1, 0, 0)),
])
def test_closed_form_centre_swap_consistency(ZA, oa, ob, ZB, oc, od):
    """(cd|ab) = (-1)^{la+lb+lc+ld} (ab|cd).

    Regions A/B/C are NOT symmetric under x <-> y, so swapping which centre sits
    at the origin routes the computation through different region logic. Any
    error in the region bookkeeping breaks this.
    """
    R = sp.Rational(5, 2)
    v1 = float(sp.re(sp.N(aabb_closed_form(ZA, oa, ob, ZB, oc, od, R), 30)))
    v2 = float(sp.re(sp.N(aabb_closed_form(ZB, oc, od, ZA, oa, ob, R), 30)))
    phase = (-1) ** (oa[1] + ob[1] + oc[1] + od[1])
    assert abs(v1 - phase * v2) < 1e-12, f"{v1} vs {phase * v2}"


def test_closed_form_large_R_multipole_limit():
    """R * (AA|BB) -> q_A q_B as R -> oo; unit charges here, so 1."""
    e = aabb_closed_form(Z1, (2, 1, 0), (2, 1, 0), Z2, (1, 0, 0), (1, 0, 0))
    assert sp.limit(sp.expand(e) * R_s, R_s, sp.oo) == 1


def test_closed_form_respects_the_M_selection_rule():
    """M_A + M_B != 0 is killed by the phi integral, so the quartet vanishes."""
    e = aabb_closed_form(Z1, (2, 1, 1), (2, 1, -1), Z2, (2, 1, 0), (2, 1, 0),
                         sp.Rational(5, 2))
    assert sp.simplify(e) == 0


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
