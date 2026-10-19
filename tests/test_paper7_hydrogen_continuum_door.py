"""Paper 7 sec:continuum -- releasing n into the hydrogen continuum.

Record: docs/predictions/2026-10-18_hydrogen_continuum.md (HC1-HC4, committed before the
computation in fb4fdd85). Written as its own step after the result; it does not import the
computation script. Each test names the wrong answer it rejects.

  HC1  the S^3 zonal function continued by n -> i rho, chi -> i chi is i times the H^3 zonal
       function; the Casimir n^2 - 1 continues to -(1 + rho^2).
  HC2  the s-wave contact density per unit energy is Z/pi on the bound side for every n, and
       (Z/pi)/(1 - e^(-2 pi eta)) on the continuum side.
  HC3  the factor at the level of the radial functions: the bound radial density per unit energy,
       continued to n = i eta, times 1/(1 - e^(-2 pi eta)) is the energy-normalized continuum
       density, at every r and every l; the two product identities that make it l-independent;
       the bound 1s -> np closed form against exact integration; Stobbe's formula against direct
       quadrature; the threshold bracket.
  HC4  the 1s continuum density at eta = 1, with the factor and e^(-4 eta arccot eta) stripped off,
       is the rational 16/3: no Coulomb phase survives in the squared quantity.

Coulomb functions are mpmath's, whose eta is DLMF's: eta_DLMF = -Z/k for attraction.
"""
from __future__ import annotations

import mpmath as mp
import sympy as sp
from sympy.physics.hydrogen import R_nl


def _factor(eta):
    return 1 / (1 - mp.e ** (-2 * mp.pi * eta))


# ---------------------------------------------------------------- HC1


def test_hc1_zonal_continuation_and_casimir():
    """Rejects a continuation onto the sphere's own zonal form sin(rho chi)/sin(chi) (no curvature
    flip) and a continued Casimir other than -(1 + rho^2)."""
    n, chi, rho = sp.symbols("n chi rho", positive=True)
    f = sp.sin(n * chi) / sp.sin(chi)
    assert sp.simplify(sp.diff(sp.sin(chi) ** 2 * sp.diff(f, chi), chi) / sp.sin(chi) ** 2
                       + (n**2 - 1) * f) == 0
    g = sp.sin(rho * chi) / sp.sinh(chi)
    assert sp.simplify(sp.diff(sp.sinh(chi) ** 2 * sp.diff(g, chi), chi) / sp.sinh(chi) ** 2
                       + (rho**2 + 1) * g) == 0
    cont = f.subs(n, sp.I * rho).subs(chi, sp.I * chi)
    assert sp.simplify(cont / g - sp.I) == 0
    assert sp.expand((n**2 - 1).subs(n, sp.I * rho)) == -(1 + rho**2)


# ---------------------------------------------------------------- HC2


def test_hc2_contact_density_bound_side_is_Z_over_pi():
    """Rejects a bound contact density per unit energy that depends on n or differs from Z/pi."""
    Z, r = sp.symbols("Z r", positive=True)
    for n in range(1, 7):
        dens = sp.limit(R_nl(n, 0, r, Z), r, 0) ** 2 / (4 * sp.pi) * n**3 / Z**2
        assert sp.simplify(dens - Z / sp.pi) == 0


def test_hc2_contact_density_continuum_side():
    """Rejects (Z/pi) e^(-2 pi eta)/(1 - e^(-2 pi eta)), the form the sign of eta flips into."""
    mp.mp.dps = 30
    for Z in (1, 2):
        for k in (mp.mpf("0.5"), mp.mpf(2)):
            eta = Z / k
            rho = mp.mpf(10) ** -25
            C0 = mp.coulombf(0, -eta, rho) / rho
            contact = (2 / (mp.pi * k)) * (C0 * k) ** 2 / (4 * mp.pi)
            assert abs(contact / (Z / mp.pi * _factor(eta)) - 1) < mp.mpf(10) ** -20


# ---------------------------------------------------------------- HC3


def test_hc3_product_identities_behind_l_independence():
    """Rejects an l-dependent factor: the bound normalization's prod_{s=-l..l}(n + s) continues to
    i eta (-1)^l prod_{s=1..l}(s^2 + eta^2), and |Gamma(l+1+i eta)|^2 carries the same product."""
    n, eta = sp.symbols("n eta", positive=True)
    for l in range(5):
        P = sp.prod([n + s for s in range(-l, l + 1)])
        target = sp.I * eta * (-1) ** l * sp.prod([s**2 + eta**2 for s in range(1, l + 1)])
        assert sp.expand(P.subs(n, sp.I * eta) - target) == 0
    mp.mp.dps = 30
    for e in (mp.mpf("0.3"), mp.mpf(2)):
        for l in range(5):
            g2 = abs(mp.gamma(l + 1 + 1j * e)) ** 2
            rhs = mp.pi * e / mp.sinh(mp.pi * e) * mp.fprod([s**2 + e**2 for s in range(1, l + 1)])
            assert abs(g2 / rhs - 1) < mp.mpf(10) ** -25


def _bound_density_per_energy(n, l: int, Z: int, r):
    """R_nl(r)^2 dn/dE as an analytic function of n: R = N (2Zr/n)^l e^(-Zr/n) M(l+1-n, 2l+2, 2Zr/n),
    N^2 = (2Z/n)^3 prod_{s=-l..l}(n+s) / (2n ((2l+1)!)^2), dn/dE = n^3/Z^2."""
    P = mp.fprod([n + s for s in range(-l, l + 1)])
    N2 = (2 * Z / n) ** 3 * P / (2 * n * mp.factorial(2 * l + 1) ** 2)
    x = 2 * Z * r / n
    R = x**l * mp.e ** (-x / 2) * mp.hyp1f1(l + 1 - n, 2 * l + 2, x)
    return N2 * R**2 * n**3 / Z**2


def test_hc3_continued_radial_density_times_factor_is_continuum():
    """Rejects a factor that depends on l or on r, and the sign-flipped factor
    e^(-2 pi eta)/(1 - e^(-2 pi eta)) (from taking mpmath's eta as +Z/k)."""
    mp.mp.dps = 30
    for Z in (1, 2):
        for k in (mp.mpf("0.5"), mp.mpf(2)):
            eta = Z / k
            for l in range(4):
                for r in (mp.mpf("0.7"), mp.mpf(3)):
                    b = _bound_density_per_energy(1j * eta, l, Z, r)
                    assert abs(mp.im(b)) < mp.mpf(10) ** -20 * abs(b)
                    cont = (2 / (mp.pi * k)) * (mp.coulombf(l, -eta, k * r) / r) ** 2
                    assert abs(cont / (mp.re(b) * _factor(eta)) - 1) < mp.mpf(10) ** -20


def test_hc3_bound_1s_np_closed_form_exact():
    """Rejects a misprinted bound closed form, for example (n-1)^(2n-3) in place of (n-1)^(2n-4)."""
    r = sp.symbols("r", positive=True)
    for n in range(2, 7):
        I = sp.integrate(R_nl(1, 0, r, 1) * R_nl(n, 1, r, 1) * r**3, (r, 0, sp.oo))
        f = sp.Rational(2, 3) * (sp.Rational(1, 2) - sp.Rational(1, 2 * n * n)) * I**2
        closed = sp.Integer(2) ** 8 * n**5 * sp.Integer(n - 1) ** (2 * n - 4) / (3 * sp.Integer(n + 1) ** (2 * n + 4))
        assert sp.simplify(f - closed) == 0


def _stobbe(eta, Z=1):
    return (2**8 * eta**8 * mp.e ** (-4 * eta * mp.acot(eta))
            / (3 * Z**2 * (1 + eta**2) ** 4) * _factor(eta))


def _dfdE_1s_quad(k, Z: int = 1):
    """1s -> E p by direct quadrature against the energy-normalized regular Coulomb function."""
    eta = Z / k
    R1s = lambda r: 2 * mp.mpf(Z) ** 1.5 * mp.e ** (-Z * r)
    u = lambda r: mp.sqrt(2 / (mp.pi * k)) * mp.coulombf(1, -eta, k * r)
    I = mp.quad(lambda r: R1s(r) * u(r) * r**2, [0, 1, 4, 16, 64])
    return mp.mpf(2) / 3 * (k**2 / 2 + mp.mpf(Z) ** 2 / 2) * I**2


def test_hc3_stobbe_against_quadrature_and_hc4_rational_strip():
    """Rejects Stobbe's formula with e^(-2 eta arccot eta) or without the factor, and any Coulomb
    phase in the squared quantity: at eta = 1 the stripped density is exactly 16/3."""
    mp.mp.dps = 20
    for k in (mp.mpf(1), mp.mpf("0.5")):
        q = _dfdE_1s_quad(k)
        assert abs(q / _stobbe(1 / k) - 1) < mp.mpf(10) ** -14
    q1 = _dfdE_1s_quad(mp.mpf(1))
    stripped = q1 / (_factor(mp.mpf(1)) * mp.e ** (-4 * mp.acot(1)))
    assert abs(stripped - mp.mpf(16) / 3) < mp.mpf(10) ** -14


def test_hc3_threshold_bracket():
    """Rejects a discontinuity at threshold: at n = eta = 1000 the bound and continuum 1s densities
    bracket 2^8 e^(-4)/3 from above and below and agree to 1e-5."""
    mp.mp.dps = 30
    n = mp.mpf(1000)
    bound = 2**8 * n**5 * (n - 1) ** (2 * n - 4) / (3 * (n + 1) ** (2 * n + 4)) * n**3
    cont = _stobbe(mp.mpf(1000))
    thr = 2**8 * mp.e ** -4 / 3
    assert bound > thr > cont
    assert abs(bound - cont) < mp.mpf(10) ** -5
