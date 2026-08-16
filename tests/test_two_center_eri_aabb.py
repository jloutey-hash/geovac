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


# ------------------------------------- increment 2: the hybrid class (AA|AB)

def test_e1_moment_matches_quadrature_both_branches():
    """int_0^X t^n e^{-ct} E_1(at) dt -- the LOG-producing family.

    Includes c < 0, where the combined rate mu = a + c goes negative. That is
    not exotic: the mirror hybrid class (AB|BB) puts the one-center pair on the
    LIGHT centre, so the naive ln((a+c)/a) form would take a log of a negative
    number. The Ein formulation is branch-safe and this pins it.
    """
    from scipy import integrate
    from scipy.special import exp1
    from geovac.two_center_eri import e1_moment

    for n in (0, 1, 2, 3):
        for c, a, X in ((1.3, 2.0, 3.0), (-1.0, 3.0, 2.0),
                        (-3.0, 1.0, 2.5), (-2.0, 2.0, 3.0)):
            got = float(sp.N(e1_moment(n, sp.nsimplify(c), sp.nsimplify(a),
                                       sp.nsimplify(X)), 30))
            ref, _ = integrate.quad(
                lambda t, n=n, c=c, a=a: t ** n * np.exp(-c * t) * exp1(a * t),
                0, X, limit=500, epsabs=1e-16, epsrel=1e-14)
            assert abs(got - ref) <= 1e-11 * max(abs(ref), 1.0), \
                f"n={n} c={c} a={a} X={X}: {got} vs {ref}"


def test_e1_moment_shifted_matches_quadrature():
    """int_0^X t^n e^{-ct} E_1(a(t+s)) dt with s > 0 -- no logarithm.

    The j = 0 term of the sum is an ordinary term here (in the unshifted case it
    is consumed cancelling a divergence). Dropping it was a real bug; this pins
    it.
    """
    from scipy import integrate
    from scipy.special import exp1
    from geovac.two_center_eri import e1_moment_shifted

    for n in (0, 1, 2, 3):
        for c, a, s, X in ((1.3, 2.0, 6.0, 3.0), (0.7, 3.0, 3.0, 1.5)):
            got = float(sp.N(e1_moment_shifted(
                n, sp.nsimplify(c), sp.nsimplify(a), sp.nsimplify(s),
                sp.nsimplify(X)), 30))
            ref, _ = integrate.quad(
                lambda t, n=n, c=c, a=a, s=s:
                    t ** n * np.exp(-c * t) * exp1(a * (t + s)),
                0, X, limit=400, epsabs=1e-18, epsrel=1e-15)
            assert abs(got - ref) <= 1e-12 * abs(ref), \
                f"n={n} c={c} a={a} s={s} X={X}: {got} vs {ref}"


@pytest.mark.parametrize("name,ZA,oa,ob,oc,ZB,od", [
    ("(1s 1s|1s 1s_B)",   Fraction(3), (1, 0, 0), (1, 0, 0), (1, 0, 0),
     Fraction(1), (1, 0, 0)),
    ("(1s 1s|2p0 2p0_B)", Fraction(3), (1, 0, 0), (1, 0, 0), (2, 1, 0),
     Fraction(1), (2, 1, 0)),
    ("(2s 2s|2p1 2p1_B)", Fraction(3), (2, 0, 0), (2, 0, 0), (2, 1, 1),
     Fraction(2), (2, 1, 1)),
    ("(1s 1s|3d1 2p1_B)", Fraction(3), (1, 0, 0), (1, 0, 0), (3, 2, 1),
     Fraction(2), (2, 1, 1)),
])
def test_hybrid_closed_form_matches_quadrature(name, ZA, oa, ob, oc, ZB, od):
    """s-type one-center pair: exact against the independent quadrature route."""
    from geovac.two_center_eri import hybrid_closed_form, hybrid_quadrature
    R = sp.Rational(5, 2)
    got = float(sp.re(sp.N(hybrid_closed_form(ZA, oa, ob, oc, ZB, od, R), 30)))
    ref = hybrid_quadrature(ZA, oa, ob, oc, ZB, od, 2.5)
    assert abs(got - ref) < 1e-11, f"{name}: {got} vs {ref}"


def test_hybrid_s_type_is_elementary():
    """min r_A power = -2(l_a + l_b), so an s-type pair stays E_1-free.

    Depends ONLY on the one-center pair -- l_c and l_d cancel out -- which is
    why d functions on the far centre do not spoil it.
    """
    from geovac.two_center_eri import hybrid_closed_form
    for oc, od in (((1, 0, 0), (1, 0, 0)), ((2, 1, 0), (2, 1, 0)),
                   ((3, 2, 1), (2, 1, 1))):
        e = hybrid_closed_form(Fraction(3), (1, 0, 0), (1, 0, 0), oc,
                               Fraction(2), od, sp.Rational(5, 2))
        assert not e.atoms(sp.expint), f"E_1 in s-type hybrid {oc}/{od}"
        assert not e.atoms(sp.log), f"log in s-type hybrid {oc}/{od}"


@pytest.mark.parametrize("name,oa,ob,oc,ZB,od", [
    ("(2p0 2p0|1s 1s_B)",   (2, 1, 0), (2, 1, 0), (1, 0, 0), Fraction(1), (1, 0, 0)),
    ("(2p0 2p0|2p0 1s_B)",  (2, 1, 0), (2, 1, 0), (2, 1, 0), Fraction(1), (1, 0, 0)),
    ("(2p1 2p1|1s 1s_B)",   (2, 1, 1), (2, 1, 1), (1, 0, 0), Fraction(1), (1, 0, 0)),
    ("(2p1 2p1|2p1 2p1_B)", (2, 1, 1), (2, 1, 1), (2, 1, 1), Fraction(2), (2, 1, 1)),
    ("(3d0 3d0|1s 1s_B)",   (3, 2, 0), (3, 2, 0), (1, 0, 0), Fraction(1), (1, 0, 0)),
])
def test_hybrid_l_gt_0_via_shells(name, oa, ob, oc, ZB, od):
    """l > 0 on the one-center pair, through the shell reformulation.

    The first row is the quartet that BLOCKED the direct route: V_L is regular
    at r -> 0 but its split into q_L r^{-(L+1)} + e^{-br}(Laurent) is not, and
    the r_A lower limit passes through zero at r_B = R. Carrying the shell radius
    as a parameter makes r_A -> 0 fall in the inside branch, power +L, so nothing
    spurious is generated.
    """
    from geovac.two_center_eri import hybrid_closed_form, hybrid_quadrature
    R = sp.Rational(5, 2)
    got = float(sp.re(sp.N(hybrid_closed_form(Fraction(3), oa, ob, oc, ZB, od,
                                              R), 30)))
    ref = hybrid_quadrature(Fraction(3), oa, ob, oc, ZB, od, 2.5)
    assert abs(got - ref) < 1e-11, f"{name}: {got} vs {ref}"


def test_hybrid_two_routes_agree_on_the_s_type_overlap():
    """The direct (V_L) and shell routes overlap at l_a = l_b = 0; they must agree.

    Independent derivations of the same number -- the direct route never forms a
    shell integral and the shell route never forms V_L.
    """
    from geovac.two_center_eri import _hybrid_direct, hybrid_closed_form_shell
    R = sp.Rational(5, 2)
    for oc, od in (((1, 0, 0), (1, 0, 0)), ((2, 1, 0), (2, 1, 0))):
        a = float(sp.re(sp.N(_hybrid_direct(
            Fraction(3), (1, 0, 0), (1, 0, 0), oc, Fraction(1), od, R), 30)))
        b = float(sp.re(sp.N(hybrid_closed_form_shell(
            Fraction(3), (1, 0, 0), (1, 0, 0), oc, Fraction(1), od, R), 30)))
        assert abs(a - b) < 1e-13, f"{oc}/{od}: direct {a} vs shell {b}"


def test_hybrid_seed_set_is_E1_plus_log_for_l_gt_0():
    """The corrected seed set, realized.

    Phase 0-h originally claimed {E_1} alone, having checked only the r_B + R
    endpoint; the |r_B - R| endpoint passes through zero and contributes a
    logarithm. The built closed form carries exactly exp, E_1 and log -- and the
    s-type case stays elementary, since the seed depends only on the one-center
    pair.
    """
    from geovac.two_center_eri import hybrid_closed_form
    R = sp.Rational(5, 2)
    e = hybrid_closed_form(Fraction(3), (2, 1, 0), (2, 1, 0), (1, 0, 0),
                           Fraction(1), (1, 0, 0), R)
    names = {type(f).__name__ for f in e.atoms(sp.Function)}
    assert "expint" in names and "log" in names, names
    assert names <= {"exp", "expint", "log"}, names


@pytest.mark.parametrize("p", [0, 2, -1, -2, -4])
def test_finite_power_exp_both_signs(p):
    """int_lo^hi r^p e^{-br} dr for either sign of b.

    b < 0 is reached whenever the r_B lower limit |r_A - R| contributes
    e^{+a_d r_A}, which is why upper_integral's E_1 branch is not enough there.
    """
    from scipy import integrate
    from geovac.two_center_eri import finite_power_exp
    for b, lo, hi in ((1.3, 0.4, 2.7), (-1.1, 0.5, 2.0), (-2.5, 0.3, 1.8)):
        got = float(sp.N(finite_power_exp(p, sp.nsimplify(b), sp.nsimplify(lo),
                                          sp.nsimplify(hi)), 30))
        ref, _ = integrate.quad(
            lambda r, p=p, b=b: r ** p * np.exp(-b * r), lo, hi,
            limit=400, epsabs=1e-15, epsrel=1e-14)
        assert abs(got - ref) <= 1e-11 * max(abs(ref), 1.0), \
            f"p={p} b={b}: {got} vs {ref}"


# ---------------------- increment 3a: exchange-class spheroidal expansion

@pytest.mark.parametrize("oa,ob", [
    ((1, 0, 0), (1, 0, 0)), ((2, 1, 0), (1, 0, 0)), ((2, 1, 1), (2, 1, 1)),
    ((2, 1, 1), (2, 1, -1)), ((3, 2, 2), (2, 1, 1)), ((3, 1, -1), (3, 2, 1)),
])
def test_spheroidal_product_is_polynomial_and_exact(oa, ob):
    """conj(chi_a^A) chi_b^B is polynomial in (xi, eta) x separable exponentials.

    This is the foundation of the whole exchange class: it is what makes the
    Neumann route tractable at all. Phase 0-e argued it on paper and only ever
    checked 1s, where it is trivial.

    Also pins the parity fact the construction relies on: |m_a| + |m_b| + |sigma|
    is even with sigma = m_a - m_b, so the leftover rho half-power becomes an
    integer once the kernel's own (1-eta^2)^{|sigma|/2} is folded in. Cases with
    a genuine half-integer half_power (e.g. 3/2) are included on purpose.
    """
    from geovac.two_center_eri import (eta_s, two_center_spheroidal_product,
                                       radial_norm, radial_poly, xi_s)
    Rv = 2.5
    P, half, p, q = two_center_spheroidal_product(
        Fraction(3), oa, Fraction(1), ob, sp.Rational(5, 2))
    assert P.is_polynomial(xi_s, eta_s), f"{oa}x{ob}: not polynomial"

    sigma = oa[2] - ob[2]
    assert (abs(oa[2]) + abs(ob[2]) + abs(sigma)) % 2 == 0, "parity claim broken"
    assert (half + sp.Rational(abs(sigma), 2)).is_integer, "half-power not integral"

    def chi(Z, n, l, m, r, th, ph):
        c, a = radial_poly(Z, n, l)
        N = radial_norm(Z, n, l)
        rad = sum(float(N * cc) * r ** k for k, cc in c.items()) \
            * np.exp(-float(a) * r)
        return rad * complex(sp.Ynm(l, m, th, ph).expand(func=True).evalf())

    Pf = sp.lambdify((xi_s, eta_s), P, "numpy")
    for xv, ev, phv in ((1.7, 0.3, 0.7), (2.4, -0.6, 2.1), (1.15, 0.85, 4.4)):
        rA, rB = Rv * (xv + ev) / 2, Rv * (xv - ev) / 2
        zA = Rv * (1 + xv * ev) / 2
        thA = np.arccos(np.clip(zA / rA, -1, 1))
        thB = np.arccos(np.clip((zA - Rv) / rB, -1, 1))
        direct = np.conj(chi(Fraction(3), *oa, rA, thA, phv)) \
            * chi(Fraction(1), *ob, rB, thB, phv)
        built = (Pf(xv, ev) * ((xv ** 2 - 1) * (1 - ev ** 2)) ** float(half)
                 * np.exp(-float(p) * xv - float(q) * ev)
                 * np.exp(1j * (ob[2] - oa[2]) * phv))
        assert abs(direct - built) < 1e-12, f"{oa}x{ob} at ({xv},{ev}): {direct} vs {built}"


# ------------------ increment 3b: exchange class assembled for general (l,m)

def test_exchange_ml_conservation():
    """sigma = m_a - m_b must equal sigma = m_d - m_c, or the quartet vanishes.

    Two independent phi integrals impose it, one per electron.
    """
    from geovac.two_center_eri import exchange_value
    v = exchange_value(Fraction(3), (2, 1, 1), (1, 0, 0),
                       Fraction(1), (1, 0, 0), (1, 0, 0), 3.0, tau_max=4)
    assert v == 0.0


def test_exchange_reproduces_phase0e_partial_sum():
    """The general (l, m) assembly must reduce to the sigma=0 1s result.

    Pinned against the Phase 0-e partial sum at tau = 6, which was itself
    validated against McMurchie-Davidson. This is the regression that catches a
    general-case rewrite silently breaking the case that already worked.
    """
    from geovac.two_center_eri import exchange_value
    v = exchange_value(Fraction(3), (1, 0, 0), (1, 0, 0),
                       Fraction(1), (1, 0, 0), (1, 0, 0), 3.0, tau_max=6)
    assert abs(v - 0.0063039227) < 1e-9, v


@pytest.mark.slow
def test_exchange_sigma_nonzero_vs_md():
    """sigma != 0 against an independent engine -- the leg that actually bites.

    A sigma=0 check cannot catch an error in any sigma-dependent factor: the
    (-1)^sigma, the [(tau-s)!/(tau+s)!]^2, or the P^mu / Q^mu conventions on
    (1,oo) vs (-1,1). Routed through 2p_{+1} = -(px + i py)/sqrt(2); axial
    symmetry cancels the px/py cross terms and equates the diagonal ones, leaving
    (px_A 1s_B | 1s_A px_B).

    n_gauss >= 10 per the Phase 0-h finding that the 6-Gaussian default is far
    too loose for two-centre overlap densities -- and exchange is the worst case,
    since BOTH densities are two-centre.
    """
    from geovac import noci_engine as E
    from geovac.two_center_eri import exchange_value
    shapes = {}
    for kind, (l, nr) in (("1s", (0, 1)), ("2p", (1, 2))):
        arr, dco, _q = E.fit_sto_shape(l, nr, n_gauss=10)
        shapes[kind] = (arr, dco)
    pa, pb = np.array([0., 0., 0.]), np.array([0., 0., 3.])

    def B(c, kind, zeta, lmn):
        return E.sto_shape_basis(c, kind, zeta, shapes, lmn)

    got = exchange_value(Fraction(3), (2, 1, 1), (1, 0, 0),
                         Fraction(1), (1, 0, 0), (2, 1, 1), 3.0, tau_max=10)
    md = E.eri_md(B(pa, "2p", 1.5, (1, 0, 0)), B(pb, "1s", 1.0, (0, 0, 0)),
                  B(pa, "1s", 3.0, (0, 0, 0)), B(pb, "2p", 0.5, (1, 0, 0)))
    assert abs(got - md) < 1e-5, f"sigma=1: {got} vs md {md}"


# ------------- increment 3c: the ordered xi integral closes, at weight 1

@pytest.mark.parametrize("n", [0, 1, 2, 3])
def test_log_moment_matches_quadrature(n):
    """int_0^oo t^n e^{-ct} ln t dt -- the family that carries Euler's gamma."""
    from scipy import integrate
    from geovac.two_center_eri import log_moment
    for c in (0.7, 1.3, 2.5):
        got = float(sp.N(log_moment(n, sp.nsimplify(c)), 30))
        ref, _ = integrate.quad(
            lambda t, n=n, c=c: t ** n * np.exp(-c * t) * np.log(t),
            0, np.inf, limit=400, epsabs=1e-14, epsrel=1e-13)
        assert abs(got - ref) <= 1e-11 * max(abs(ref), 1.0), f"n={n} c={c}"


@pytest.mark.parametrize("n", [0, 1, 2, 3])
def test_log_shift_moment_matches_quadrature(n):
    """int_0^oo t^n e^{-ct} ln(t+s) dt, s > 0 -- no divergence is ever formed."""
    from scipy import integrate
    from geovac.two_center_eri import log_shift_moment
    for c, s in ((1.3, 2.0), (0.7, 2.0), (2.5, 1.0)):
        got = float(sp.N(log_shift_moment(n, sp.nsimplify(c), sp.nsimplify(s)), 30))
        ref, _ = integrate.quad(
            lambda t, n=n, c=c, s=s: t ** n * np.exp(-c * t) * np.log(t + s),
            0, np.inf, limit=400, epsabs=1e-14, epsrel=1e-13)
        assert abs(got - ref) <= 1e-11 * max(abs(ref), 1.0), f"n={n} c={c} s={s}"


def test_ordered_xi_integral_closes_at_weight_one():
    """THE periods result: the ordered double integral closes, and at weight 1.

    It is an iterated integral over a simplex -- the shape that defines a period
    -- so the natural expectation is that it lands one rung up, at weight 2
    (dilogarithms, zeta(2)). It does not. The closed form carries exp, E_1, log
    and Euler's gamma, and nothing higher.

    Established at sigma = 0. For sigma != 0 the d^sigma Q_tau derivatives put
    poles at xi = +-1; increment 3a's parity fact makes the net exponent there
    exactly 0, but that is argued rather than verified.
    """
    from geovac.two_center_eri import ordered_xi_closed
    p1, p2 = sp.symbols("p1 p2", positive=True)
    e = ordered_xi_closed(p1, p2)
    names = {type(f).__name__ for f in e.atoms(sp.Function)}
    assert names <= {"exp", "expint", "log"}, f"unexpected functions {names}"
    assert not (names & {"polylog", "dilog", "lerchphi", "zeta"}), \
        "a weight-2 object reached the closed form"
    assert e.has(sp.EulerGamma), "gamma should be present -- it is the xi=1 endpoint"


@pytest.mark.parametrize("p", [sp.Rational(3, 2), sp.Integer(2), sp.Integer(3)])
def test_ordered_xi_closed_matches_quadrature(p):
    """The closed form against a direct nested quadrature of the same object."""
    from scipy import integrate
    from scipy.special import eval_legendre
    from geovac.two_center_eri import ordered_xi_closed

    pf = float(p)

    def outer(x1):
        lo, _ = integrate.quad(lambda x2: np.exp(-pf * x2), 1.0, x1,
                               epsabs=1e-13, epsrel=1e-12, limit=200)
        hi, _ = integrate.quad(
            lambda x2: np.exp(-pf * x2) * 0.5 * np.log((x2 + 1) / (x2 - 1)),
            x1, np.inf, epsabs=1e-13, epsrel=1e-12, limit=200)
        q0 = 0.5 * np.log((x1 + 1) / (x1 - 1))
        return np.exp(-pf * x1) * (q0 * lo + eval_legendre(0, x1) * hi)

    ref, _ = integrate.quad(outer, 1.0, np.inf, epsabs=1e-12, epsrel=1e-11,
                            limit=200)
    got = float(sp.re(sp.N(ordered_xi_closed(p, p), 30)))
    assert abs(got - ref) <= 1e-10 * max(abs(ref), 1e-3), f"p={pf}: {got} vs {ref}"


# ---------- increment 3d: sigma != 0 does not break the weight-1 result

@pytest.mark.parametrize("tau,sigma,H", [
    (1, 1, 1), (2, 1, 1), (2, 2, 2), (3, 1, 1), (3, 2, 2), (3, 3, 3), (4, 2, 2),
])
def test_sigma_poles_are_absorbed_by_the_prefactor(tau, sigma, H):
    """d^sigma Q_tau has poles of order up to sigma at xi = +-1; (xi^2-1)^H eats them.

    This is the one place a weight-2 object could have entered the exchange
    class. It cannot, because H - |sigma| = (|m_a|+|m_b|-|m_a-m_b|)/2 >= 0 by the
    triangle inequality, so both pieces come out polynomial.
    """
    from geovac.two_center_eri import Q_tau_sigma_split, xi_s
    a, b = Q_tau_sigma_split(tau, sigma)
    pref = (xi_s ** 2 - 1) ** H
    pl = sp.cancel(sp.expand(pref * a))
    pr = sp.cancel(sp.together(pref * b))
    assert pl.is_polynomial(xi_s), f"tau={tau} s={sigma}: log coefficient not polynomial"
    assert pr.is_polynomial(xi_s), f"tau={tau} s={sigma}: pole survived the prefactor"


def test_triangle_bound_forbids_new_poles():
    """H - |sigma| >= 0 for every (m_a, m_b) -- the reason 3d closes.

    H = (|m_a|+|m_b|)/2 + |sigma|/2 with sigma = m_a - m_b, so the net exponent
    against a pole of order |sigma| is (|m_a|+|m_b|-|m_a-m_b|)/2, non-negative by
    the triangle inequality and zero exactly when m_a, m_b have opposite signs.
    """
    worst = None
    for ma in range(-3, 4):
        for mb in range(-3, 4):
            sig = ma - mb
            H = sp.Rational(abs(ma) + abs(mb), 2) + sp.Rational(abs(sig), 2)
            assert H.is_integer, f"half-power not integral at ({ma},{mb})"
            net = H - abs(sig)
            worst = net if worst is None else min(worst, net)
            assert net >= 0, f"({ma},{mb}) gives net exponent {net} < 0"
    assert worst == 0, f"expected the bound to be attained, got {worst}"


# ------- increment 3e: the general (tau, sigma, H, j) assembly loop

def test_ordered_xi_general_reproduces_the_tau0_special_case():
    """The general loop must reduce to 3c's hand-built tau = 0 form."""
    from geovac.two_center_eri import ordered_xi_closed, ordered_xi_general
    for p in (sp.Rational(3, 2), sp.Integer(2)):
        g = float(sp.re(sp.N(ordered_xi_general(0, 0, 0, 0, 0, 0, p, p), 30)))
        c = float(sp.re(sp.N(ordered_xi_closed(p, p), 30)))
        assert abs(g - c) < 1e-13, f"p={p}: general {g} vs special {c}"


@pytest.mark.parametrize("tau,sg,H1,H2,j1,j2,s1,s2", [
    (2, 0, 1, 1, 1, 0, "2", "3/2"),
    (1, 1, 1, 1, 0, 0, "3/2", "2"),
    (2, 1, 1, 1, 1, 1, "2", "2"),
    (2, 2, 2, 2, 0, 0, "3/2", "5/2"),
])
def test_ordered_xi_general_matches_nested_quadrature(tau, sg, H1, H2, j1, j2,
                                                      s1, s2):
    """General parameters, including p1 != p2, against direct nested quadrature."""
    from scipy import integrate
    from geovac.two_center_eri import (Q_tau_sigma_split, ordered_xi_general,
                                       xi_s)
    P1, P2 = sp.Rational(s1), sp.Rational(s2)
    p1, p2 = float(P1), float(P2)
    Pp = sp.diff(sp.legendre(tau, xi_s), xi_s, sg) if sg else sp.legendre(tau, xi_s)
    a, b = Q_tau_sigma_split(tau, sg)
    Pf = sp.lambdify(xi_s, Pp, "numpy")
    af, bf = sp.lambdify(xi_s, a, "numpy"), sp.lambdify(xi_s, b, "numpy")

    def DQ(u):
        return float(af(u)) * 0.5 * np.log((u + 1) / (u - 1)) + float(bf(u))

    def outer(x1):
        lo, _ = integrate.quad(
            lambda x2: x2 ** j2 * (x2 ** 2 - 1) ** H2 * np.exp(-p2 * x2)
            * float(Pf(x2)), 1.0, x1, epsabs=1e-13, epsrel=1e-12, limit=200)
        hi, _ = integrate.quad(
            lambda x2: x2 ** j2 * (x2 ** 2 - 1) ** H2 * np.exp(-p2 * x2) * DQ(x2),
            x1, np.inf, epsabs=1e-13, epsrel=1e-12, limit=200)
        return (x1 ** j1 * (x1 ** 2 - 1) ** H1 * np.exp(-p1 * x1)
                * (DQ(x1) * lo + float(Pf(x1)) * hi))

    ref, _ = integrate.quad(outer, 1.0, np.inf, epsabs=1e-12, epsrel=1e-11,
                            limit=200)
    got = float(sp.re(sp.N(ordered_xi_general(tau, sg, H1, H2, j1, j2, P1, P2), 30)))
    assert abs(got - ref) <= 1e-10 * max(abs(ref), 1e-3), f"{got} vs {ref}"


def test_ordered_xi_general_stays_weight_one():
    """General parameters must not introduce a weight-2 object either."""
    from geovac.two_center_eri import ordered_xi_general
    p1, p2 = sp.symbols("p1 p2", positive=True)
    e = ordered_xi_general(2, 1, 1, 1, 1, 1, p1, p2)
    names = {type(f).__name__ for f in e.atoms(sp.Function)}
    assert names <= {"exp", "expint", "log"}, f"unexpected {names}"
    assert not (names & {"polylog", "dilog", "zeta"}), "weight-2 object appeared"


def test_aabb_closed_form_is_lindemann_separable():
    """(AA|BB) is zero-DECIDABLE by the same argument that decides S and h.

    Paper 58's S and h rows are decided because each entry is
    e^{-p}(U e^q + V e^{-q}) with U, V exact rationals, so it vanishes iff
    U = V = 0 by Lindemann. This pins the analogous structure for (AA|BB):
    every pi CANCELS -- the harmonic normalisations against the 4pi/(2L+1) of
    the multipole potential -- leaving A_0(R) + sum_j A_j(R) e^{-lambda_j R}
    with A_j rational and lambda_j rational. For algebraic R the exponents are
    distinct algebraic numbers, so {1, e^{-lambda_j R}} is linearly independent
    over the algebraics and the sum vanishes iff every A_j does.

    The load-bearing part is that the pi power is COMMON across terms. If it
    mixed, vanishing would need independence of pi against the exponentials,
    which is open (Schanuel territory) rather than classical.
    """
    from geovac.two_center_eri import aabb_closed_form, R_s
    Z2_, Z3_ = Fraction(2), Fraction(3)
    for args in ((Z1, (1, 0, 0), (1, 0, 0), Z1, (1, 0, 0), (1, 0, 0)),
                 (Z3_, (2, 1, 0), (2, 1, 0), Z1, (1, 0, 0), (1, 0, 0)),
                 (Z3_, (2, 1, 1), (2, 1, 1), Z2_, (2, 1, 0), (2, 1, 0))):
        e = sp.expand(aabb_closed_form(*args))
        powers = set()
        for term in sp.Add.make_args(e):
            pw = 0
            for f in sp.Mul.make_args(term):
                if f == sp.pi:
                    pw += 1
                elif f.is_Pow and f.base == sp.pi:
                    pw += f.exp
            powers.add(sp.nsimplify(pw))
        assert len(powers) == 1, f"pi power mixes across terms: {powers}"
        # every exponential rate rational, so the exponents are algebraic in R
        for f in e.atoms(sp.exp):
            rate = sp.simplify(-sp.diff(f.args[0], R_s))
            assert rate.is_rational, f"non-rational exponential rate {rate}"
        assert not e.atoms(sp.expint), "E_1 present -- not Lindemann-separable"
        assert not e.atoms(sp.log), "log present -- not Lindemann-separable"


def test_decided_census_aabb_block_has_no_accidental_zeros():
    """Step 2: every permitted (AA|BB) entry is DECIDED nonzero.

    Paper 58's `g` row is counted -- it says which entries the rules PERMIT,
    never whether each permitted entry is nonzero, and its Gaussian corroboration
    decides zeros by a 1e-10 float threshold. With (AA|BB) in closed form and
    Lindemann-separable, each entry can be DECIDED: group by exponential rate,
    ask whether every rational coefficient vanishes.

    Result on the census configuration (Z_A=3, Z_B=1, n_max=2, R=3): all 195
    permitted entries are genuinely nonzero -- no accidental zeros, no symmetry
    zeros the counting missed. So the counted density is exactly the true
    density on this block.

    Sampled here (the full 195-entry sweep is ~110 s and lives in
    debug/step2_decided_census.py); the sample is chosen to include l>0 on both
    sides and m != 0, where an accidental cancellation would be most likely.
    """
    from geovac.two_center_eri import aabb_closed_form, R_s
    Z3_ = Fraction(3)
    sample = [((1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0)),
              ((2, 1, 0), (2, 1, 0), (2, 1, 0), (2, 1, 0)),
              ((2, 1, 1), (2, 1, 0), (2, 1, 0), (2, 1, 1)),
              ((1, 0, 0), (2, 1, 1), (2, 1, 1), (2, 0, 0))]
    for a, b, c, d in sample:
        e = aabb_closed_form(Z3_, a, b, Z1, c, d)
        groups: dict = {}
        for term in sp.Add.make_args(sp.expand(e)):
            rate, coeff = sp.Integer(0), sp.Integer(1)
            for f in sp.Mul.make_args(term):
                if isinstance(f, sp.exp):
                    arg = sp.expand(f.args[0])
                    dd = -sp.diff(arg, R_s)
                    rate += dd
                    coeff *= sp.exp(sp.expand(arg + dd * R_s))
                else:
                    coeff *= f
            k = sp.nsimplify(rate)
            groups[k] = groups.get(k, 0) + coeff
        nonzero = any(sp.simplify(cc.subs(R_s, 3)) != 0 for cc in groups.values())
        assert nonzero, f"({a}{b}|{c}{d}) decided ZERO -- an accidental zero"


def test_decided_census_aabb_robust_off_census_config():
    """Bet 2: the (AA|BB) 'no accidental zeros' is a class property, not census-specific.

    The permitted set and Gaunt-zero count are pure angular data (charge/R-independent);
    only accidental *radial* zeros can move across configs. A 5x7 charge x R grid (full
    195 entries each, incl. equal-charge and rate-coincidence, plus an n_max=3
    node-bearing sample) is clean in debug/bet2_aabb_accidental_zero_sweep.py; sampled
    here at an off-census config -- equal charges Z=1, R=5/2, the degeneracy case where
    coincident exponential rates are the likeliest source of a cancellation.
    """
    from geovac.two_center_eri import aabb_closed_form, R_s
    Z1_ = Fraction(1)
    Rv = sp.Rational(5, 2)
    sample = [((1, 0, 0), (1, 0, 0), (1, 0, 0), (1, 0, 0)),
              ((2, 1, 0), (2, 1, 0), (2, 1, 0), (2, 1, 0)),
              ((2, 1, 1), (2, 1, 0), (2, 1, 0), (2, 1, 1)),
              ((1, 0, 0), (2, 1, 1), (2, 1, 1), (2, 0, 0))]
    for a, b, c, d in sample:
        e = aabb_closed_form(Z1_, a, b, Z1_, c, d)   # equal charges -> degenerate rates
        groups: dict = {}
        for term in sp.Add.make_args(sp.expand(e)):
            rate, coeff = sp.Integer(0), sp.Integer(1)
            for f in sp.Mul.make_args(term):
                if isinstance(f, sp.exp):
                    arg = sp.expand(f.args[0])
                    dd = -sp.diff(arg, R_s)
                    rate += dd
                    coeff *= sp.exp(sp.expand(arg + dd * R_s))
                else:
                    coeff *= f
            k = sp.nsimplify(rate)
            groups[k] = groups.get(k, 0) + coeff
        nonzero = any(sp.simplify(cc.subs(R_s, Rv)) != 0 for cc in groups.values())
        assert nonzero, f"({a}{b}|{c}{d}) accidental zero at equal-charge/R=5/2"


def test_exchange_gamma_survives_every_tau():
    """Bet 3: standalone Euler gamma appears at EVERY tau, so it does not cancel in
    assembly -- the exchange class is gamma-blocked for decidability (gamma is not known
    irrational), on top of the E_1 wall it shares with the hybrid class. This is what
    makes the cross-class census DECISION transcendence-hard rather than 'ordinary work':
    weight-one already contains E_1 and gamma, and only (AA|BB) (pure {exp}, pi cancelled)
    is Lindemann-decidable.
    """
    from geovac.two_center_eri import ordered_xi_general
    p1, p2 = sp.symbols("p1 p2", positive=True)
    for tau in (0, 1, 2, 3):
        e = ordered_xi_general(tau, 0, 0, 0, 0, 0, p1, p2)
        assert e.has(sp.EulerGamma), f"gamma absent at tau={tau} -- would break the wall claim"


def test_three_center_1e_closes_weight_one_gamma_free():
    """Bet 1: the 3-centre ONE-electron integral <chi_Y|-Z_X/r_X|chi_Z> closes at
    weight one AND gamma-free for a source pinned off the foci axis (xi_X > 1).

    Structure: with the source pinned, the ordered Neumann split is at a FIXED xi_X,
    so per (tau, density-monomial) the xi integral is
        Q_tau(s0) * int_1^{s0} poly*P_tau e^{-p xi} dxi          (finite, elementary)
      + P_tau(s0) * int_{s0}^inf poly*(-W_tau) e^{-p xi} dxi     (elementary)
      + P_tau(s0) * int_{s0}^inf poly*P_tau*Q0 e^{-p xi} dxi     (the only transcendence)
    composing the engine's validated weight-1 moments. The two-centre exchange carried
    Euler gamma from its xi=1 endpoint; a source at xi_X>1 never reaches it, so gamma
    drops out -> {exp, E_1, ln}, a proper subset of the exchange seed set. So weight-1
    (the arc's central property) survives the third centre. Full 3-D validation vs an
    independent reference: debug/bet1_three_center_1e_{reference,assembly,symbolic}.py.
    """
    from scipy import integrate
    from scipy.special import eval_legendre
    from geovac.two_center_eri import (finite_power_exp, log_shift_moment,
                                       upper_integral)
    xi = sp.Symbol("xi", positive=True)

    def xi_pinned(j1, H1, tau, p, s0):
        poly = sp.expand(xi ** j1 * (xi ** 2 - 1) ** H1)
        Ptau = sp.legendre(tau, xi)
        W = sum(sp.legendre(k - 1, xi) * sp.legendre(tau - k, xi) / sp.Integer(k)
                for k in range(1, tau + 1))
        Q0s = sp.log((s0 + 1) / (s0 - 1)) / 2
        Ps = Ptau.subs(xi, s0)
        Qs = Ps * Q0s - (W.subs(xi, s0) if tau else 0)
        pP = sp.Poly(sp.expand(poly * Ptau), xi)
        inner = Qs * sum(b * finite_power_exp(n, p, sp.Integer(1), s0)
                         for (n,), b in pP.terms())
        negW = sp.Poly(sp.expand(-poly * W), xi) if tau else sp.Poly(0, xi)
        oelem = sum(b * upper_integral(n, p, s0) for (n,), b in negW.terms())
        olog = sp.Integer(0)
        for (n,), b in pP.terms():
            Lp = sp.exp(-p * s0) * sum(sp.binomial(n, m) * s0 ** (n - m)
                                       * log_shift_moment(m, p, s0 + 1) for m in range(n + 1))
            Lm = sp.exp(-p * s0) * sum(sp.binomial(n, m) * s0 ** (n - m)
                                       * log_shift_moment(m, p, s0 - 1) for m in range(n + 1))
            olog += b * sp.Rational(1, 2) * (Lp - Lm)
        return inner + Ps * (oelem + olog)

    # (i) weight verdict, symbolic p and s0
    p_s, s0_s = sp.symbols("p s0", positive=True)
    for (j1, H1, tau) in [(0, 0, 0), (0, 0, 1), (2, 1, 2), (0, 0, 3)]:
        e = xi_pinned(j1, H1, tau, p_s, s0_s)
        funcs = {type(f).__name__ for f in e.atoms(sp.Function)}
        assert funcs <= {"exp", "expint", "log"}, f"unexpected {funcs}"
        assert not e.has(sp.EulerGamma), "gamma should drop out for a source at xi_X>1"
        assert not (funcs & {"polylog", "dilog", "zeta"}), "weight-2 object appeared"

    # (ii) correctness vs a direct quadrature of the SAME pinned integral
    p_v, s0_v = sp.Rational(7, 4), sp.Rational(53, 32)     # p, xi_X (>1)
    pf, s0f = float(p_v), float(s0_v)
    for (j1, H1, tau) in [(0, 0, 1), (2, 1, 2)]:
        got = float(sp.re(sp.N(xi_pinned(j1, H1, tau, p_v, s0_v), 30)))

        def integrand(x, use_Q):
            base = x ** j1 * (x * x - 1) ** H1 * np.exp(-pf * x) * eval_legendre(tau, x)
            if not use_Q:
                return base
            Q0 = 0.5 * np.log((x + 1) / (x - 1))
            Wv = sum(eval_legendre(k - 1, x) * eval_legendre(tau - k, x) / k
                     for k in range(1, tau + 1))
            return (x ** j1 * (x * x - 1) ** H1 * np.exp(-pf * x)
                    * (eval_legendre(tau, x) * Q0 - Wv))
        inner, _ = integrate.quad(lambda x: integrand(x, False), 1.0, s0f,
                                  epsabs=1e-13, epsrel=1e-12, limit=200)
        outer, _ = integrate.quad(lambda x: integrand(x, True), s0f, np.inf,
                                  epsabs=1e-13, epsrel=1e-12, limit=200)
        Q0s = 0.5 * np.log((s0f + 1) / (s0f - 1))
        Ps = eval_legendre(tau, s0f)
        Ws = sum(eval_legendre(k - 1, s0f) * eval_legendre(tau - k, s0f) / k
                 for k in range(1, tau + 1))
        ref = (Ps * Q0s - Ws) * inner + Ps * outer
        assert abs(got - ref) < 1e-9, f"(j1={j1},H1={H1},tau={tau}): {got} vs {ref}"
