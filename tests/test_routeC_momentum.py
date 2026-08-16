"""Route C -- momentum-space two-body 3-centre ERI (XY|XZ): regression pins.

Backs the OPEN-frontier sprint result (build plan section 10.5,
debug/routeC_momentum_poc.py, memo debug/sprint_routeC_momentum_memo.md).

Self-contained (no debug/ import, per the transient-dir policy).  Pins the two
load-bearing, VALIDATED facts -- the method is exact and the single-dispersion
factor is a Bessel K0.  The transcendence WEIGHT of the coupled two-scale radial
integral is deliberately NOT pinned: it is the genuinely open question.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from geovac import noci_engine as GE

X = np.array([0.0, 0.0, 0.0])
Y = np.array([0.0, 0.0, 2.0])
Z = np.array([1.5, 0.0, 0.5])
GROUND_TRUTH = 0.2049417218262148   # eri_md on 8-Gaussian 1s fits


def _kgrid(n_rad, n_cos, n_phi, k_scale):
    u, wu = leggauss(n_rad); u = 0.5 * (u + 1.0); wu = 0.5 * wu
    k = k_scale * u / (1.0 - u)
    wk = wu * (k_scale / (1.0 - u) ** 2) * k ** 2
    c, wc = leggauss(n_cos); st = np.sqrt(1.0 - c * c)
    phi = 2.0 * np.pi * np.arange(n_phi) / n_phi; wphi = 2.0 * np.pi / n_phi
    kx = k[:, None, None] * (st[None, :, None] * np.cos(phi)[None, None, :])
    ky = k[:, None, None] * (st[None, :, None] * np.sin(phi)[None, None, :])
    kz = k[:, None, None] * (c[None, :, None] * np.ones_like(phi)[None, None, :])
    kv = np.stack([kx.ravel(), ky.ravel(), kz.ravel()], axis=-1)
    w = (wk[:, None, None] * wc[None, :, None]
         * np.ones_like(phi)[None, None, :] * wphi).ravel()
    return kv, w


def _gauss_density_ft(bra, ket):
    A, B = bra.center, ket.center
    ai, aj = bra.alphas[:, None], ket.alphas[None, :]
    ca, cb = bra.coeffs, ket.coeffs
    p = ai + aj
    K = np.exp(-ai * aj / p * float(np.dot(A - B, A - B)))
    coef = (ca[:, None] * cb[None, :]) * K * (np.pi / p) ** 1.5
    P = (ai[:, :, None] * A + aj[:, :, None] * B) / p[:, :, None]
    inv4p = 1.0 / (4.0 * p)

    def ft(kv):
        k2 = np.einsum("ij,ij->i", kv, kv)
        phase = np.exp(1j * np.einsum("ni,abi->nab", kv, P))
        gauss = np.exp(-np.einsum("n,ab->nab", k2, inv4p))
        return np.einsum("ab,nab,nab->n", coef, phase, gauss)

    return ft


def test_momentum_reproduces_three_center_eri():
    """(1/2pi^2) int d3k/k^2 rho1~ conj(rho2~) == eri_md, to k-grid precision.

    Validates the momentum-space formula, the FT/Coulomb conventions, and the
    k-integrator.  The three centres enter only as phases -- the coordinate wall
    (no spheroidal system for three foci) is dissolved.
    """
    sh_a, sh_d, _ = GE.fit_sto_shape(0, 1, n_gauss=8)
    shapes = {"1s": (sh_a, sh_d)}
    a = GE.sto_shape_basis(X, "1s", 1.0, shapes, (0, 0, 0))
    b = GE.sto_shape_basis(Y, "1s", 1.0, shapes, (0, 0, 0))
    c = GE.sto_shape_basis(X, "1s", 1.0, shapes, (0, 0, 0))
    d = GE.sto_shape_basis(Z, "1s", 1.0, shapes, (0, 0, 0))
    ref = GE.eri_md(a, b, c, d)

    kv, w = _kgrid(90, 32, 32, 4.0)
    k2 = np.einsum("ij,ij->i", kv, kv)
    f1 = _gauss_density_ft(a, b)(kv)
    f2 = _gauss_density_ft(c, d)(kv)
    mom = (1.0 / (2.0 * np.pi ** 2)) * np.sum(w * (f1 * np.conj(f2)) / k2)

    assert abs(mom.imag) < 1e-10
    assert abs(mom.real - ref) < 1e-9, (mom.real, ref)


def test_single_dispersion_factor_is_besselK0():
    """int_0^inf cos(kb) e^{-D sqrt(c k^2+m^2)}/sqrt(c k^2+m^2) dk
         = (1/sqrt c) K0( (m/sqrt c) sqrt(c D^2 + b^2) ).

    A single dispersion factor closes in closed form via the Fock substitution
    k*sqrt(c)=m*sinh(theta) -- a momentum-space Coulomb-Sturmian (Bessel K0),
    weight one.  Two DIFFERENT scales c1=s(1-s), c2=t(1-t) coupled by j0(k|W|) is
    the three-centre obstruction; this pins the one-scale baseline.
    """
    from scipy.integrate import quad
    from scipy.special import k0

    b, D, c, m = 1.3, 2.0, 0.21, 1.0
    val, _ = quad(lambda k: math.cos(k * b) * math.exp(-D * math.sqrt(c * k * k + m * m))
                  / math.sqrt(c * k * k + m * m), 0, np.inf, limit=400)
    closed = (1.0 / math.sqrt(c)) * k0((m / math.sqrt(c)) * math.sqrt(c * D * D + b * b))
    assert abs(val - closed) < 1e-10, (val, closed)


def test_three_center_kernel_is_elliptic():
    """The two-scale radial kernel lives on y^2=(c1 k^2+1)(c2 k^2+1): a genus-1
    elliptic curve for c1 != c2, degenerating to genus-0 (rational) at c1 == c2.

    Decisive witness: the D=0 period is a COMPLETE ELLIPTIC INTEGRAL K when the two
    Fock scales differ, and pi/(2 sqrt c) when they coincide.  This is the exact
    transcendence obstruction of the two-body three-centre ERI -- the third centre
    raises the transcendence from the two-centre engine's genus-0 {E1, ln, gamma}
    to genus-1 (elliptic).  m(s,t)=1-c_min/c_max degenerates only on the measure-zero
    locus c1=c2 (equal Fock scales = a shared hypersphere).
    """
    from scipy.integrate import quad
    from scipy.special import ellipk

    # off-diagonal: c1 != c2 -> complete elliptic integral, nondegenerate modulus
    c1, c2 = 0.11, 0.21
    num, _ = quad(lambda k: 1.0 / math.sqrt((c1 * k * k + 1) * (c2 * k * k + 1)),
                  0, np.inf, limit=400)
    a1, a2 = 1.0 / math.sqrt(c1), 1.0 / math.sqrt(c2)
    m = 1.0 - (a2 / a1) ** 2
    closed = (1.0 / (a1 * math.sqrt(c1 * c2))) * ellipk(m)
    assert abs(num - closed) < 1e-9, (num, closed)
    assert 0.0 < m < 1.0                      # nondegenerate elliptic

    # diagonal: c1 == c2 -> rational (genus 0), elementary period pi/(2 sqrt c)
    c = 0.19
    numd, _ = quad(lambda k: 1.0 / ((c * k * k + 1)), 0, np.inf, limit=400)
    assert abs(numd - math.pi / (2.0 * math.sqrt(c))) < 1e-9


# ---------------------------------------------------------------------------
# ABW-obstruction pins (Paper 59 sec:obstruction, memo "(a) ABW push").
# The one-mass slice L(D,rho) = int_1^inf e^{-Dx} dx / sqrt(Q),
# Q=(x^2-1)(rho x^2+1-rho), rho=c2/c1, is NOT an elliptic dilogarithm of
# elementary/period data: the ABW variation-of-parameters mechanism is
# obstructed because the modulus Picard-Fuchs source is in-module.
# ---------------------------------------------------------------------------

def test_modulus_pf_annihilates_periods():
    """M_rho = rho(1-rho) d^2/drho^2 + (1-2rho) d/drho - 1/4 (Legendre PF) has the
    two periods K(rho), K(1-rho) as homogeneous solutions (Paper 59 eq:modpf).

    These are the ABW homogeneous solutions; they live in the MODULUS variable
    (not the Laplace variable D) -- the reason the naive 4th-order-D-PF attack is
    mis-aimed.
    """
    import mpmath as mp
    mp.mp.dps = 30

    def M_rho(f, rho):
        rho = mp.mpf(rho)
        return (rho * (1 - rho) * mp.diff(f, rho, 2)
                + (1 - 2 * rho) * mp.diff(f, rho, 1) - mp.mpf(1) / 4 * f(rho))

    for rho in ['0.37', '0.5', '0.23']:
        assert abs(M_rho(lambda r: mp.ellipk(r), rho)) < mp.mpf(10) ** -10
        assert abs(M_rho(lambda r: mp.ellipk(1 - r), rho)) < mp.mpf(10) ** -10


def test_modulus_source_is_exact_derivative_g1_zero():
    """M_rho[Q^{-1/2}] = d/dx[ g Q^{-1/2} ],  g = -x(x^2-1)/(4(rho x^2+1-rho)),
    with g(1)=0 (Paper 59 eq:exactform).

    g(1)=0 is the mechanism of the obstruction: modulus-differentiation produces
    NO branch-point boundary term, so the source has no elementary "tadpole"
    subtopology (unlike the sunrise) and stays on the elliptic curve.
    """
    import mpmath as mp
    mp.mp.dps = 30
    rho = mp.mpf('0.37')

    def Q(x):
        return (x * x - 1) * (rho * x * x + 1 - rho)

    def F(x):                      # Q^{-1/2}
        return 1 / mp.sqrt(Q(x))

    def g(x):
        return -x * (x * x - 1) / (4 * (rho * x * x + 1 - rho))

    # g vanishes at the branch point x=1
    assert abs(g(mp.mpf(1))) < mp.mpf(10) ** -25

    # exact-form identity at a sample x: apply M_rho (in rho) to F, compare to
    # d/dx[g F] (in x); both at fixed (x, rho).
    x0 = mp.mpf('1.7')

    def F_of_rho(r):
        return 1 / mp.sqrt((x0 * x0 - 1) * (r * x0 * x0 + 1 - r))

    lhs = (rho * (1 - rho) * mp.diff(F_of_rho, rho, 2)
           + (1 - 2 * rho) * mp.diff(F_of_rho, rho, 1) - mp.mpf(1) / 4 * F_of_rho(rho))
    rhs = mp.diff(lambda x: g(x) * F(x), x0)
    assert abs(lhs - rhs) < mp.mpf(10) ** -18, (lhs, rhs)

    # SYMBOLIC pin (backs the paper's [SYMBOLIC PROOF] tag, not just a spot check):
    # M_rho[Q^{-1/2}] - d/dx[g Q^{-1/2}] simplifies to EXACTLY 0, and g(1)=0 exactly.
    import sympy as sp
    xs, rs = sp.symbols('x rho', positive=True)
    Qs = (xs * xs - 1) * (rs * xs * xs + 1 - rs)
    Fs = Qs ** sp.Rational(-1, 2)
    gs = -xs * (xs * xs - 1) / (4 * (rs * xs * xs + 1 - rs))
    M_rho_F = rs * (1 - rs) * sp.diff(Fs, rs, 2) + (1 - 2 * rs) * sp.diff(Fs, rs) - sp.Rational(1, 4) * Fs
    assert sp.simplify(M_rho_F - sp.diff(gs * Fs, xs)) == 0
    assert sp.simplify(gs.subs(xs, 1)) == 0


def test_single_bessels_and_periods_are_not_D_solutions():
    """The 4th-order D-operator O[y] = Dp y'''' + 2p y''' + D(1-2p) y'' + (1-2p) y'
    - D(1-p) y  (p=rho) does NOT annihilate the single Bessels {K0(D),I0(D),
    J0(wD),Y0(wD)} nor the periods (D-constants).

    This is why the naive ABW attack (VoP on the D-PF, expecting periods/Bessels as
    homogeneous solutions) is mis-aimed -- the curve moves with the modulus, not D.
    Backs the eq:modpf-paragraph [MEASURED] claims (were driver-only).
    """
    import mpmath as mp
    mp.mp.dps = 30
    rho = mp.mpf('0.37')
    w = mp.sqrt((1 - rho) / rho)
    D = mp.mpf('1.7')

    def O(y):
        return (D * rho * mp.diff(y, D, 4) + 2 * rho * mp.diff(y, D, 3)
                + D * (1 - 2 * rho) * mp.diff(y, D, 2) + (1 - 2 * rho) * mp.diff(y, D, 1)
                - D * (1 - rho) * y(D))

    for y in (lambda D: mp.besselk(0, D), lambda D: mp.besseli(0, D),
              lambda D: mp.besselj(0, w * D), lambda D: mp.bessely(0, w * D)):
        assert abs(O(y)) > mp.mpf(10) ** -3          # NOT a D-solution (coupled sectors)

    # a period is D-constant: O[const] = -D(1-rho) exactly, not 0
    assert abs(O(lambda D: mp.mpf(1)) - (-D * (1 - rho))) < mp.mpf(10) ** -20


@pytest.mark.slow
def test_modulus_source_is_in_module_not_period_plus_bessel():
    """The modulus source S = M_rho[L] closes EXACTLY only in L's own Gauss-Manin
    D-module {L,L',L'',L'''} (Paper 59 eq:inmodule); it does NOT admit a finite
    period+Bessel (ABW subtopology) decomposition at equal parameter budget.

    In-module (Hyp A) fits far tighter than period+Bessel (Hyp C) with the SAME
    number of parameters -> no elementary ABW inhomogeneity -> VoP is circular.
    """
    import mpmath as mp
    mp.mp.dps = 30
    rho = mp.mpf('0.37')
    w = mp.sqrt((1 - rho) / rho)

    def _q(gu):
        return mp.quad(gu, [0, 1, 3, 8, 20, mp.inf])

    def data(D):
        D = mp.mpf(D)
        def mom(n):
            return _q(lambda u: (1 + u * u) ** n * 2 * mp.e ** (-D * (1 + u * u))
                      / mp.sqrt((2 + u * u) * (rho * (1 + u * u) ** 2 + 1 - rho)))
        Lk = [mom(0), -mom(1), mom(2), -mom(3)]
        dLdr = _q(lambda u: -u * u * mp.e ** (-D * (1 + u * u)) * mp.sqrt(2 + u * u)
                  / (rho * (1 + u * u) ** 2 + 1 - rho) ** mp.mpf('1.5'))
        d2 = _q(lambda u: (mp.mpf(3) / 2) * u ** 4 * mp.e ** (-D * (1 + u * u))
                * (2 + u * u) ** mp.mpf('1.5') / (rho * (1 + u * u) ** 2 + 1 - rho) ** mp.mpf('2.5'))
        S = rho * (1 - rho) * d2 - (2 * rho - 1) * dLdr - mp.mpf(1) / 4 * Lk[0]
        # BOTH Bessel sectors (the paper's claim): K0,K1 at arg D and J0,Y0 at arg wD
        bess = [mp.besselk(0, D), mp.besselk(1, D), mp.besselj(0, w * D), mp.bessely(0, w * D)]
        return Lk, bess, S

    def fit(kind, deg):
        gens = (lambda Lk, bs: Lk) if kind == 'A' else (lambda Lk, bs: [Lk[0], Lk[1]] + list(bs))
        ncols = len((gens([0] * 4, [0] * 4))) * (deg + 1)
        Dall = [mp.mpf('0.4') + mp.mpf('0.37') * i for i in range(ncols + 4)]
        dat = [data(D) for D in Dall]
        def row(D, Lk, bs):
            return [D ** k * gg for gg in gens(Lk, bs) for k in range(deg + 1)]
        A = mp.matrix([row(Dall[i], dat[i][0], dat[i][1]) for i in range(ncols)])
        b = mp.matrix([dat[i][2] for i in range(ncols)])
        coef = mp.lu_solve(A, b)
        return max(abs(sum(coef[j] * row(Dall[i], dat[i][0], dat[i][1])[j]
                           for j in range(ncols)) - dat[i][2])
                   for i in range(ncols, ncols + 4))

    # The discriminator is SATURATION vs degree: exact finite closure is FLAT,
    # a convergent approximation keeps IMPROVING.
    res_A2, res_A3 = fit('A', 2), fit('A', 3)
    res_C2, res_C3 = fit('C', 2), fit('C', 3)

    assert res_A2 < mp.mpf(10) ** -25                    # in-module: closes to the floor
    assert res_A2 / res_A3 < 10                          # Hyp A SATURATES (flat = exact)
    assert res_C2 / res_C3 > mp.mpf(10) ** 3             # Hyp C IMPROVES (approximation)
    assert res_C2 / res_A2 > mp.mpf(10) ** 8             # equal budget: period+Bessel far worse


# ---------------------------------------------------------------------------
# Cosmic-Galois pins (Paper 59 sec:modular): the elliptic family is the
# Legendre / Gamma(2) universal family, its CM-fiber periods are Gamma-values,
# and the integrated observable is a Gamma(2) multiple modular value.
# ---------------------------------------------------------------------------

def test_cosmic_galois_family_is_gamma2():
    """lambda(tau(rho)) = 1 - rho exactly (Paper 59 eq:lambda_rho): the family is the
    Legendre / Gamma(2) universal family, tau(rho)=i K(rho)/K(1-rho)."""
    import mpmath as mp
    mp.mp.dps = 30

    def lam_theta(tau):
        q = mp.e ** (1j * mp.pi * tau)
        return (mp.jtheta(2, 0, q) / mp.jtheta(3, 0, q)) ** 4

    for rho in ['0.37', '0.6', '0.15']:
        r = mp.mpf(rho)
        tau = 1j * mp.ellipk(r) / mp.ellipk(1 - r)
        lam = mp.re(lam_theta(tau))
        assert abs(lam - (1 - r)) < mp.mpf(10) ** -20, (rho, lam)


def test_cosmic_galois_cm_periods_are_gamma_values():
    """CM-fiber periods are Gamma-values (Chowla-Selberg), at two fundamental
    discriminants (Paper 59 sec:modular): disc -4 (tau=i, rho=1/2, in the physical
    domain) and disc -8 (tau=i sqrt2)."""
    import mpmath as mp
    mp.mp.dps = 40
    # disc -4
    K4 = mp.ellipk(mp.mpf('0.5'))
    cf4 = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(mp.pi))
    assert abs(K4 - cf4) < mp.mpf(10) ** -30
    # disc -8: lambda_2 = 3 - 2 sqrt2
    lam2 = 3 - 2 * mp.sqrt(2)
    K8 = mp.ellipk(lam2)
    cf8 = (mp.sqrt(1 + mp.sqrt(2)) * mp.gamma(mp.mpf(1) / 8) * mp.gamma(mp.mpf(3) / 8)
           / (2 ** mp.mpf('3.25') * mp.sqrt(mp.pi)))
    assert abs(K8 - cf8) < mp.mpf(10) ** -30


def test_cosmic_galois_integrated_value():
    """The integrated collinear T2 = 0.39535576590171392... (Paper 59 sec:modular),
    via the fast evaluator (GL tensor + sin^2 substitution + s<->t symmetry).  A
    float64 low-order reproduction pins it to ~6 digits (the high-precision 17-digit
    value is in debug/routeC_fast_evaluator.py)."""
    from numpy.polynomial.legendre import leggauss
    from scipy.integrate import quad

    def Pf(x, k):
        c = x * (1 - x); D = np.sqrt(c * k * k + 1)
        return c * np.exp(-D) * (1 / D ** 3 + 3 / D ** 4 + 3 / D ** 5)

    def Jf(s, t):
        b = s + t
        f = lambda k: (np.sin(k * b) / (k * b) if k * b > 1e-12 else 1.0) * Pf(s, k) * Pf(t, k)
        v, _ = quad(f, 0, np.inf, limit=400)
        return v

    N = 24
    x, w = leggauss(N)
    phi = (np.pi / 4) * (x + 1); wphi = (np.pi / 4) * w
    s = np.sin(phi) ** 2; jac = np.sin(2 * phi)
    tot = 0.0
    for i in range(N):
        for j in range(N):
            tot += wphi[i] * jac[i] * wphi[j] * jac[j] * Jf(s[i], s[j])
    T2 = (8 / np.pi) * tot
    assert abs(T2 - 0.395355766) < 1e-6, T2
