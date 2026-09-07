"""Paper 61 (companion of Paper 59) -- the resurgent-skeleton pattern across the corpus's two-electron
integral classes (CHANGELOG v4.103.0; drivers debug/aha_t3_*.py).

Pins the two cheap anchors of the 4/4 pattern (algebraic Stokes data up to a
singularity-type-fixed pi-power; boundary period forced at weight <= 1):

  1. The Paper-18 Level-2 exchange seed e^a E_1(a): Borel transform 1/(1+zeta),
     one simple pole, Stokes constant exactly 2*pi*i (rational multiple 1) --
     realized as the E_1 branch cut, Disc E_1(-x) = 2*pi*i.
  2. The Paper-58 hybrid two-center ERI class {E_1, ln}: the sector Stokes
     constants cancel pairwise (P_(c,lam1) = -P_(c,lam2)), so the physical ERI
     is CUT-FREE in R -- the median resummation IS the closed form.
"""
import mpmath as mp
import sympy as sp


def test_e1_seed_stokes_constant():
    """E_1(-x - i0) - E_1(-x + i0) = 2*pi*i exactly (the rank-1 anchor)."""
    mp.mp.dps = 40
    eps = mp.mpf("1e-30")
    for x in (mp.mpf("1.5"), mp.mpf("2.5")):
        jump = mp.e1(mp.mpc(-x, -eps)) - mp.e1(mp.mpc(-x, eps))
        assert abs(jump - 2 * mp.pi * mp.mpc(0, 1)) < mp.mpf("1e-25")
        # and for the corpus object f = e^a E_1(a) the discontinuity is
        # 2*pi*i * e^a -- the trans-series residue is the rational 1
        f_jump = mp.exp(-x) * jump
        assert abs(f_jump - 2 * mp.pi * mp.mpc(0, 1) * mp.exp(-x)) < mp.mpf("1e-25")


def test_hybrid_eri_cut_free():
    """The hybrid {E_1, ln} closed form is single-valued in R: its E_1 cuts
    cancel pairwise across sectors (|Disc F| / |F| at machine-negligible level).
    Quartet (2p0 2p0 | 1s 1s_B), Z_A = 3, Z_B = 1."""
    from geovac.two_center_eri import R_s, hybrid_closed_form

    mp.mp.dps = 50
    expr = hybrid_closed_form(3, (2, 1, 0), (2, 1, 0), (1, 0, 0), 1, (1, 0, 0))
    f = sp.lambdify(R_s, expr, modules="mpmath")
    eps = mp.mpf("1e-35")
    for x in (mp.mpf("1.3"), mp.mpf("2.7")):
        above = f(mp.mpc(-x, eps))
        below = f(mp.mpc(-x, -eps))
        disc = abs(above - below)
        scale = abs(above)
        assert disc / scale < mp.mpf("1e-25")


def test_exchange_class_normal_form():
    """The exchange kernel (Increment 3c) collapses to the single-sector
    reduced normal form

      F = [e^{-AR}(gamma + ln(kappa R)) + e^{(a-b)R}E1(2aR)
           + e^{(b-a)R}E1(2bR) - e^{AR}E1(2AR)] / (2 a b R^2),

    A = a+b, kappa = 2ab/A -- an EXACT symbolic identity.  kappa is the
    charge-weighted product of the Borel positions (the exponents ARE the
    Stokes charges), and gamma appears only in the boundary bundle
    (gamma + ln kappa R): gamma is coordinate bookkeeping, not a Borel-plane
    transcendental."""
    from geovac.two_center_eri import R_s, ordered_xi_closed

    R = R_s
    for a, b in ((sp.Integer(1), sp.Integer(2)),
                 (sp.Rational(3, 2), sp.Rational(5, 2))):
        A = a + b
        kappa = 2 * a * b / A
        nf = (sp.exp(-A * R) * (sp.EulerGamma + sp.log(kappa * R))
              + sp.exp((a - b) * R) * sp.E1(2 * a * R)
              + sp.exp((b - a) * R) * sp.E1(2 * b * R)
              - sp.exp(A * R) * sp.E1(2 * A * R)) / (2 * a * b * R ** 2)
        diff = sp.simplify(sp.expand(ordered_xi_closed(a * R, b * R) - nf))
        assert diff == 0


def test_exchange_class_is_multivalued_in_R():
    """Unlike the hybrid class (cut-free), the exchange object is genuinely
    multivalued in complex R -- the price of the R-dependent log, costing no
    new transcendental (the Borel transform is real-analytic on (0, inf))."""
    mp.mp.dps = 40
    a, b = mp.mpf(1), mp.mpf("1.5")
    A = a + b
    kap = 2 * a * b / A

    def F(R):
        return (mp.exp(-A * R) * (mp.euler + mp.log(kap * R))
                + mp.exp((a - b) * R) * mp.e1(2 * a * R)
                + mp.exp((b - a) * R) * mp.e1(2 * b * R)
                - mp.exp(A * R) * mp.e1(2 * A * R)) / (2 * a * b * R ** 2)

    eps = mp.mpf("1e-30")
    x = mp.mpf("1.3")
    disc = F(mp.mpc(-x, eps)) - F(mp.mpc(-x, -eps))
    assert abs(disc) > mp.mpf("1e-6") * abs(F(mp.mpc(-x, eps)))
