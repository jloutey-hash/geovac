"""Paper 59 -- the three-electron chain integral is a K3 whose period is
"elliptic squared" (Bailey), and its one-centre limit is rational.

Backs the pre-registered predictions P3a/P3b/P3c/P3d of
docs/predictions/2026-10-13_three_electron_k3.md (registered 2026-10-13
before any three-electron computation; outcomes recorded in CHANGELOG
v5.20.14).  Self-contained (no debug/ import, per the transient-dir policy);
the formulas are re-implemented from debug/three_electron_k3_check.py (F2,
F3, F4, F5) and debug/three_electron_p3d_check.py (Bailey V3), whose recorded
values are the pins below.

The integral.  The chain I = <rho_a(1) rho_b(2) rho_c(3) / (r12 r13)> with
two-centre 1s Slater densities; in momentum space
    I ~ int d3p d3q rho_b(p) rho_c(q) rho_a(p+q) / (p^2 q^2),
two momentum integrations, three dispersion roots Delta_b = sqrt(cb p^2+1),
Delta_c = sqrt(cc q^2+1), Delta_a = sqrt(ca |p+q|^2+1).

P3a  At the D = 0 fingerprint fibre the cos(p,q) angle integral is
     elementary: int dx (ca u+1)^{-5/2}, u = p^2+q^2+2pqx, equals
     (Delta_-^{-3} - Delta_+^{-3}) / (3 ca p q), Delta_+- = sqrt(ca (p+-q)^2+1).
     The reduced 2-D fibre equals the direct 3-D fibre (recorded 1e-15..4e-14).
P3b  Homogenised, the three quadrics are six lines in P^2; the double cover
     branched along them is a K3 iff no four lines are concurrent and no two
     coincide (recorded: none for any positive c).  p-fixed slices are
     genus-1 curves with p-dependent j (non-isotrivial).
P3c  One-centre limit (rational transforms, no root): the chain integral is
     RATIONAL -- 23/54, 133/144, 289/225 for zeta = (1,1,1), (1,2,3), (2,1,3).
P3d  The holomorphic K3 period
        P = int_{R^2} dp dq / sqrt((cb p^2+1)(cc q^2+1)(ca (p+q)^2+1))
          = 8/(pi sqrt(ca cb cc)) V3(1/sqrt ca, 1/sqrt cb, 1/sqrt cc),
        V3(a,b,c) = int_0^inf K0(at) K0(bt) K0(ct) dt,
     and Bailey (1936), restated as Bailey-Borwein-Broadhurst-Glasser,
     J. Phys. A 41, 205203 (2008), eq. (33), for c >= b >= a:
        (2c/pi) V3 = K(k-) K(k+') + K(k+) K(k-'),
        k+- = (sqrt((c+a)^2-b^2) +- sqrt((c-a)^2-b^2)) / (2c),  k' = sqrt(1-k^2).
     Equal-mass anchor, eq. (47): V3(1,1,1) = (pi/2) K(k3) K(k3'),
     k3 = sin(pi/12), = 3 Gamma(1/3)^6 / (32 pi 2^(2/3)).

NOTE: V3 here is int_0^inf K0 K0 K0 dt with NO extra factor of t -- the
normalisation of the P3d driver and of BBBG eq. (33) -- not int t K0 K0 K0.
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
from scipy.integrate import dblquad, tplquad

# recorded values (debug/data/three_electron_k3_check.log, _p3d_check.json)
F2_C = (0.21, 0.09, 0.16)
F2_RECORDED = 2.516012626371          # reduced == direct, |diff| 8.9e-16
P_DIRECT_RECORDED = 121.63429507482526  # P(0.21, 0.09, 0.16), 2-D quadrature
V3_111 = "6.94882278107962978943643644547"   # quad = Gamma form = K3K3' = eq.(33)


def _rho(k: float, c: float) -> float:
    return (c * k * k + 1.0) ** -2.5


def _reduced_fibre(ca: float, cb: float, cc: float) -> float:
    """2-D reduced fibre: int dp dq rho_b rho_c * [angle-integrated rho_a]."""
    def f(tq: float, tp: float) -> float:
        p = tp / (1 - tp)
        q = tq / (1 - tq)
        if p * q < 1e-14:
            ang = 2.0 * (ca * max(p, q) ** 2 + 1.0) ** -2.5
        else:
            Dm = np.sqrt(ca * (p - q) ** 2 + 1.0)
            Dp = np.sqrt(ca * (p + q) ** 2 + 1.0)
            ang = (Dm ** -3 - Dp ** -3) / (3.0 * ca * p * q)
        return _rho(p, cb) * _rho(q, cc) * ang / (1 - tp) ** 2 / (1 - tq) ** 2
    return dblquad(f, 0, 1, 0, 1, epsabs=1e-13, epsrel=1e-12)[0]


def _direct_fibre(ca: float, cb: float, cc: float) -> float:
    """Direct 3-D fibre over (p, q, x = cos theta); no angular reduction."""
    def f(x: float, tq: float, tp: float) -> float:
        p = tp / (1 - tp)
        q = tq / (1 - tq)
        u = p * p + q * q + 2 * p * q * x
        return (_rho(p, cb) * _rho(q, cc) * (ca * u + 1.0) ** -2.5
                / (1 - tp) ** 2 / (1 - tq) ** 2)
    return tplquad(f, 0, 1, 0, 1, -1, 1, epsabs=1e-12, epsrel=1e-10)[0]


# ------------------------------------------------------------- (b) P3a fibre
def test_reduced_2d_fibre_equals_direct_3d_fibre():
    """P3a, falsifier 2: the angle-reduced 2-D fibre equals the direct
    (p, q, cos theta) 3-D fibre at (ca, cb, cc) = (0.21, 0.09, 0.16), and both
    equal the recorded 2.516012626371.

    Rejects: a wrong angular reduction -- e.g. the sign flipped
    (Delta_+^{-3} - Delta_-^{-3}) gives the negative; dropping the 1/3 gives
    3x; the Yukawa-style 1/Delta instead of 1/(3 Delta^3) is off at 1e-1.
    Tolerance 1e-9 (recorded 8.9e-16; the 3-D quadrature floor is ~1e-11)."""
    red = _reduced_fibre(*F2_C)
    dire = _direct_fibre(*F2_C)
    assert abs(red - dire) < 1e-9, (red, dire, abs(red - dire))
    assert abs(red - F2_RECORDED) < 1e-10, (red, F2_RECORDED)


def test_angle_integral_is_elementary_no_log():
    """P3a, falsifier 1 (symbolic): the cos-angle antiderivative of the
    D = 0 density (c u + 1)^{-5/2}, u = p^2+q^2+2pqx, is
    -1/(3 c p q Delta^3) with Delta = sqrt(c u + 1) -- differentiating back
    in x returns the integrand exactly, and no log appears.  Same for the
    Yukawa kernel e^{-D Delta}/Delta: antiderivative -e^{-D Delta}/(D c p q).

    Rejects: an angle integral that produces a logarithm or a new root (the
    residual would be nonzero or the antiderivative would carry log)."""
    import sympy as sp
    c, D, p, q, x = sp.symbols("c D p q x", positive=True)
    u = p**2 + q**2 + 2 * p * q * x
    Delta = sp.sqrt(c * u + 1)
    anti5 = -1 / (3 * c * p * q * Delta**3)
    assert sp.simplify(sp.diff(anti5, x) - (c * u + 1) ** sp.Rational(-5, 2)) == 0
    anti_y = -sp.exp(-D * Delta) / (D * c * p * q)
    assert sp.simplify(sp.diff(anti_y, x) - sp.exp(-D * Delta) / Delta) == 0
    assert not anti5.has(sp.log) and not anti_y.has(sp.log)


# ---------------------------------------------------------- (a) P3c rational
def _one_centre_momentum(za: int, zb: int, zc: int) -> float:
    """Chain integral with one-centre 1s^2 densities via the momentum route."""
    def rt(k: float, z: float) -> float:          # FT of the density
        return 16.0 * z**4 / (k * k + 4.0 * z * z) ** 2

    def f(tq: float, tp: float) -> float:
        p = tp / (1 - tp)
        q = tq / (1 - tq)
        if p * q < 1e-14:
            ang = 2.0 * rt(max(p, q), za)
        else:
            ang = 8.0 * za**4 / (p * q) * (1.0 / ((p - q) ** 2 + 4 * za * za)
                                           - 1.0 / ((p + q) ** 2 + 4 * za * za))
        return rt(p, zb) * rt(q, zc) * ang / (1 - tp) ** 2 / (1 - tq) ** 2

    val = dblquad(f, 0, 1, 0, 1, epsabs=1e-13, epsrel=1e-12)[0]
    return 2.0 / np.pi**2 * val


def _one_centre_position(za: int, zb: int, zc: int):
    """Exact position-space value, mpmath 30 digits: int 4 pi r^2 rho_a V_b V_c."""
    import mpmath as mp
    mp.mp.dps = 30
    dens = lambda z, r: z**3 / mp.pi * mp.exp(-2 * z * r)
    pot = lambda z, r: 1 / r - mp.exp(-2 * z * r) * (z + 1 / r)
    return mp.quad(lambda r: 4 * mp.pi * r**2 * dens(za, r) * pot(zb, r) * pot(zc, r),
                   [0, 1, 4, mp.inf])


def test_one_centre_limit_is_rational():
    """P3c, falsifier 5 + 'stronger than predicted': with one-centre
    densities the chain integral is RATIONAL, and the momentum route (the
    same reduced 2-D fibre machinery) reproduces it.  Recomputed:
    zeta = (1,1,1) -> 23/54, (1,2,3) -> 133/144, (2,1,3) -> 289/225.

    Rejects: any transcendental content (a log, as Fromm-Hill's one-centre
    r12 r23 r31 integrals carry dilogs; 23/54 + log 2 = 1.119 is excluded by
    the exact-rational check), a wrong rational (the nearest denominator
    <= 1000 fraction must be exactly the registered one), and a momentum
    route that disagrees with position space.  Tolerances: position-space
    vs exact rational 1e-25 (30 digits); momentum vs rational 1e-11."""
    import mpmath as mp
    cases = {(1, 1, 1): Fraction(23, 54),
             (1, 2, 3): Fraction(133, 144),
             (2, 1, 3): Fraction(289, 225)}
    for zs, frac in cases.items():
        pos = _one_centre_position(*zs)
        assert abs(pos - mp.mpf(frac.numerator) / frac.denominator) < mp.mpf(10) ** -25, (zs, pos)
        mom = _one_centre_momentum(*zs)
        assert abs(mom - float(frac)) < 1e-11, (zs, mom, float(frac))
        assert Fraction(mom).limit_denominator(1000) == frac, (zs, mom)


# ------------------------------------------------------------ (c) P3d Bailey
def _V3_quad(a, b, c):
    import mpmath as mp
    g = lambda t: mp.besselk(0, a * t) * mp.besselk(0, b * t) * mp.besselk(0, c * t)
    return mp.quad(g, [0, mp.mpf(1) / 8, 1, 4, mp.inf])


def _V3_bailey(a, b, c):
    """BBBG (2008) eq. (33); valid in both the real-k and the complex-conjugate-k regimes."""
    import mpmath as mp
    a, b, c = sorted([mp.mpf(a), mp.mpf(b), mp.mpf(c)])          # c largest
    s1 = mp.sqrt((c + a) ** 2 - b ** 2)
    s2 = mp.sqrt((c - a) ** 2 - b ** 2 + 0j)
    kp = (s1 + s2) / (2 * c)
    km = (s1 - s2) / (2 * c)
    K = lambda k: mp.ellipk(k * k)                               # mpmath takes m = k^2
    Kc = lambda k: mp.ellipk(1 - k * k)                          # K(k')
    return mp.pi / (2 * c) * (K(km) * Kc(kp) + K(kp) * Kc(km))


def test_bailey_v3_closed_form_equals_three_bessel_quadrature():
    """P3d: BBBG eq. (33) equals the quadrature of V3 = int_0^inf K0 K0 K0 dt
    at two unequal scale triples -- (1, 2, 4), real k+-, and the registered
    (ca,cb,cc) = (0.21,0.09,0.16) masses 1/sqrt c, where k+- are complex
    conjugates -- and at equal mass against the Gamma form and K3 K3'.

    Rejects: a V3 normalisation with an extra factor t (int t K0^3 differs
    at O(1)); a wrong modulus convention (K(k) vs K(m=k^2) differs at 1e-2);
    swapping K(k-)K(k+') for K(k-)K(k-') (differs at O(1)); the 30-digit
    equal-mass anchor 6.9488227810796297894364364455 rejects any formula
    error above 1e-25.  Tolerance 1e-17 at 20 digits (recorded 1e-34 at 34)."""
    import mpmath as mp
    mp.mp.dps = 20
    tol = mp.mpf(10) ** -17
    for masses in ((1, 2, 4), tuple(1 / mp.sqrt(mp.mpf(x)) for x in F2_C)):
        vq = _V3_quad(*masses)
        vb = _V3_bailey(*masses)
        assert abs(mp.im(vb)) < tol, vb
        assert abs(vq - mp.re(vb)) < tol, (masses, vq, vb)
    mp.mp.dps = 32
    v111 = _V3_quad(1, 1, 1)
    gamma_form = 3 * mp.gamma(mp.mpf(1) / 3) ** 6 / (32 * mp.pi * mp.cbrt(4))
    k3 = mp.sin(mp.pi / 12)
    k3form = mp.pi / 2 * mp.ellipk(k3**2) * mp.ellipk(1 - k3**2)
    anchor = mp.mpf(V3_111)
    for v in (v111, gamma_form, k3form, mp.re(_V3_bailey(1, 1, 1))):
        assert abs(v - anchor) < mp.mpf(10) ** -25, (v, anchor)


def test_k3_period_reduces_to_three_scale_bessel_moment():
    """P3d, the Fourier bridge: the holomorphic K3 period
    P = int_{R^2} dp dq / sqrt((cb p^2+1)(cc q^2+1)(ca (p+q)^2+1)), evaluated
    by direct 2-D quadrature, equals 8/(pi sqrt(ca cb cc)) V3(1/sqrt ca, ...)
    at (0.21, 0.09, 0.16); recorded 121.63429507482526 vs 121.63429507464818.

    Rejects: a wrong bridge constant (8/pi -> 4/pi halves it; a missing
    1/sqrt(ca cb cc) is off by 20x) and a wrong mass map (c instead of
    1/sqrt c).  Tolerance 1e-8 (recorded 1.8e-10, set by the 2-D quadrature)."""
    import mpmath as mp
    mp.mp.dps = 20
    ca, cb, cc = F2_C

    def f(tq: float, tp: float) -> float:
        p = np.tan(tp)
        q = np.tan(tq)
        return ((1 + p * p) * (1 + q * q)
                / np.sqrt((cb * p * p + 1) * (cc * q * q + 1) * (ca * (p + q) ** 2 + 1)))

    h = np.pi / 2
    p_direct = dblquad(f, -h, h, -h, h, epsabs=1e-11, epsrel=1e-10)[0]
    m = [1 / mp.sqrt(mp.mpf(x)) for x in F2_C]
    p_bessel = float(8 / (mp.pi * mp.sqrt(mp.mpf(ca) * cb * cc)) * _V3_quad(*m))
    assert abs(p_direct - P_DIRECT_RECORDED) < 1e-8, p_direct
    assert abs(p_direct - p_bessel) < 1e-8, (p_direct, p_bessel, abs(p_direct - p_bessel))


# ------------------------------------------------------- (d) P3b six lines
def test_branch_locus_is_six_lines_no_four_concurrent():
    """P3b, falsifier 3 (exact, symbolic): homogenising the three quadrics
    cb p^2 + w^2, cc q^2 + w^2, ca (p+q)^2 + w^2 gives six lines in P^2,
        sb p +- i w = 0,  sc q +- i w = 0,  sa (p+q) +- i w = 0
    (sa, sb, sc the square roots of the scales).  The double cover of P^2
    branched along six lines is a K3 (at worst rational double points) iff
    the six lines are distinct and no four are concurrent.  Checked for every
    positive (sa, sb, sc): the six lines are pairwise distinct, and for each
    of the 15 quadruples at least one 3x3 minor is +-2i times a monomial in
    (sa, sb, sc), hence nonzero for all positive scales.  (Triple points DO
    occur on the loci sa(sb +- sc) = sb sc -- rational double points, still
    K3 -- the test does not claim otherwise.)

    Rejects: a configuration where four lines meet (the cover is not K3) or
    two coincide -- e.g. if the a-lines were sa p +- i w (no q), the quadruple
    (b+, b-, a+, a-) would be concurrent at (0,1,0) and every one of its
    minors would vanish identically."""
    import itertools
    import sympy as sp
    sa, sb, sc = sp.symbols("s_a s_b s_c", positive=True)
    I = sp.I
    lines = {"b+": (sb, 0, I), "b-": (sb, 0, -I),
             "c+": (0, sc, I), "c-": (0, sc, -I),
             "a+": (sa, sa, I), "a-": (sa, sa, -I)}
    # pairwise distinct: cross product of coefficient vectors nonzero
    for l1, l2 in itertools.combinations(lines, 2):
        cross = sp.Matrix(lines[l1]).cross(sp.Matrix(lines[l2]))
        assert any(sp.simplify(x) != 0 for x in cross), (l1, l2)
    # no four concurrent: some minor is a nonzero monomial for every quadruple
    for quad in itertools.combinations(lines, 4):
        monomial_minor = False
        for trio in itertools.combinations(quad, 3):
            det = sp.factor(sp.Matrix([lines[k] for k in trio]).det())
            if det != 0 and not det.has(sp.Add):          # +-2 i * monomial
                monomial_minor = True
                break
        assert monomial_minor, quad


def test_p_fixed_slice_is_genus_one_and_non_isotrivial():
    """P3b, falsifier 4: the p-fixed coordinate slice y^2 = (cc q^2+1)
    (ca (p+q)^2+1) is a quartic in q (genus 1), and its j-invariant depends
    on p (non-isotrivial: the surface is not a product).  Recorded at
    (ca, cc) = (0.21, 0.16): j = 1.2e7 at p = 0 down to 1.75e3 at p = 5.

    Rejects: an isotrivial family (dj/dp == 0, which would make the cover a
    product or rational rather than a genuine K3) and a degenerate slice
    (discriminant 4I^3 - J^2 = 0, j undefined)."""
    import sympy as sp
    p, q = sp.symbols("p q")
    ca, cc = sp.Rational(21, 100), sp.Rational(16, 100)
    quart = sp.Poly(sp.expand((cc * q**2 + 1) * (ca * (p + q) ** 2 + 1)), q)
    a0, a1, a2, a3, a4 = [quart.coeff_monomial(q**k) for k in (4, 3, 2, 1, 0)]
    I4 = 12 * a0 * a4 - 3 * a1 * a3 + a2**2
    J4 = (72 * a0 * a2 * a4 + 9 * a1 * a2 * a3 - 27 * a0 * a3**2
          - 27 * a4 * a1**2 - 2 * a2**3)
    disc = sp.simplify(4 * I4**3 - J4**2)
    j = 1728 * 4 * I4**3 / disc
    assert quart.degree() == 4
    for pv in (0, 1, 5):
        assert disc.subs(p, pv) != 0, pv
    assert sp.simplify(sp.diff(j, p)) != 0
    j0, j5 = float(j.subs(p, 0)), float(j.subs(p, 5))
    assert abs(j0 - 12003923.52) < 1.0, j0
    assert abs(j5 - 1753.5399) < 1e-3, j5
