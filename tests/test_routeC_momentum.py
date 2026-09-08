"""Route C -- momentum-space two-body 3-centre ERI (XY|XZ): regression pins.

Backs the OPEN-frontier sprint result (build plan section 10.5,
debug/routeC_momentum_poc.py, memo debug/sprint_routeC_momentum_memo.md).

Self-contained (no debug/ import, per the transient-dir policy).  Pins the two
load-bearing, VALIDATED facts -- the method is exact and the single-dispersion
factor is a Bessel K0.  The transcendence WEIGHT of the coupled two-scale radial
integral is deliberately NOT pinned: it is the genuinely open question.

NOTE (2026-09-07): this file also backs **Paper 61** (sec:modular,
sec:bessel_algebra, eq:lambda_rho), which was split out of Paper 59 on
2026-09-06.  The split's sweep keyed on the `test_paper59_*` filename and so
missed this file entirely -- eight loci still said "Paper 59" for labels
Paper 61 owns.  If either paper moves a section again, grep for BOTH names
here.
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

    # Witness the abstract / sec:momentum C8 headline: on a converged k-grid the
    # momentum-space evaluation reproduces the McMurchie-Davidson closed form to
    # the quoted 1.9e-14 (grid convergence 1.5e-12 @ (90,32,32) -> 7.8e-15 @
    # (128,48,48) -> 2.1e-15 @ (160,64,64); the coarse assertion above validates
    # the method/conventions, this pins the headline residual). ~1.7 s.
    kv2, w2 = _kgrid(128, 48, 48, 4.0)
    k2b = np.einsum("ij,ij->i", kv2, kv2)
    f1b = _gauss_density_ft(a, b)(kv2)
    f2b = _gauss_density_ft(c, d)(kv2)
    mom2 = (1.0 / (2.0 * np.pi ** 2)) * np.sum(w2 * (f1b * np.conj(f2b)) / k2b)
    assert abs(mom2.real - ref) < 1.9e-14, (mom2.real, ref, abs(mom2.real - ref))


def test_single_dispersion_factor_is_besselK0():
    """int_0^inf cos(kb) e^{-D sqrt(c k^2+m^2)}/sqrt(c k^2+m^2) dk
         = (1/sqrt c) K0( (m/sqrt c) sqrt(c D^2 + b^2) ).

    A single dispersion factor closes in closed form via the Fock substitution
    k*sqrt(c)=m*sinh(theta) -- a momentum-space Coulomb-Sturmian (Bessel K0),
    weight one.  Two DIFFERENT scales c1=s(1-s), c2=t(1-t) coupled by j0(k|W|) is
    the three-centre obstruction; this pins the one-scale baseline.

    High-precision (mpmath): the eq:K0 identity is exact -- residual is 0 to the
    dps=40 floor, so the test now witnesses the quoted ~6e-18, not float64 1e-10.
    """
    import mpmath as mp
    mp.mp.dps = 40
    b, D, c, m = mp.mpf('1.3'), mp.mpf(2), mp.mpf('0.21'), mp.mpf(1)
    val = mp.quad(lambda k: mp.cos(k * b) * mp.e ** (-D * mp.sqrt(c * k * k + m * m))
                  / mp.sqrt(c * k * k + m * m), [0, mp.inf])
    closed = (1 / mp.sqrt(c)) * mp.besselk(0, (m / mp.sqrt(c)) * mp.sqrt(c * D * D + b * b))
    assert abs(val - closed) < mp.mpf(10) ** -35, (val, closed)


def test_three_center_kernel_is_elliptic():
    """The two-scale radial kernel lives on y^2=(c1 k^2+1)(c2 k^2+1): a genus-1
    elliptic curve for c1 != c2, degenerating to genus-0 (rational) at c1 == c2.

    Decisive witness: the D=0 period is a COMPLETE ELLIPTIC INTEGRAL K when the two
    Fock scales differ, and pi/(2 sqrt c) when they coincide.  This is the exact
    transcendence obstruction of the two-body three-centre ERI -- the third centre
    raises the transcendence from the two-centre engine's genus-0 {E1, ln, gamma}
    to genus-1 (elliptic).  m(s,t)=1-c_min/c_max degenerates only on the measure-zero
    locus c1=c2 (equal Fock scales = a shared hypersphere).

    High-precision (mpmath): the D=0 period is the complete elliptic integral K
    to the dps=40 floor (residual 0), so the test now witnesses the quoted
    ~31-digit period, not the float64 1e-9 of the earlier scipy version.
    """
    import mpmath as mp
    mp.mp.dps = 40

    # off-diagonal: c1 != c2 -> complete elliptic integral, nondegenerate modulus
    c1, c2 = mp.mpf('0.11'), mp.mpf('0.21')
    num = mp.quad(lambda k: 1 / mp.sqrt((c1 * k * k + 1) * (c2 * k * k + 1)), [0, mp.inf])
    a1, a2 = 1 / mp.sqrt(c1), 1 / mp.sqrt(c2)
    m = 1 - (a2 / a1) ** 2                     # mpmath ellipk takes parameter m
    closed = (1 / (a1 * mp.sqrt(c1 * c2))) * mp.ellipk(m)
    assert abs(num - closed) < mp.mpf(10) ** -30, (num, closed)
    assert 0 < m < 1                           # nondegenerate elliptic

    # diagonal: c1 == c2 -> rational (genus 0), elementary period pi/(2 sqrt c)
    c = mp.mpf('0.19')
    numd = mp.quad(lambda k: 1 / (c * k * k + 1), [0, mp.inf])
    assert abs(numd - mp.pi / (2 * mp.sqrt(c))) < mp.mpf(10) ** -30


def _kspace_grid(n_rad, n_cos, n_phi, k_scale):
    """3D k-grid: radial u->k=k_scale u/(1-u), Gauss-Legendre in cos(theta),
    uniform phi (returns (kvecs, weights) with the d^3k = k^2 dk dcos dphi measure)."""
    u, wu = leggauss(n_rad); u = 0.5 * (u + 1.0); wu = 0.5 * wu
    k = k_scale * u / (1.0 - u); wk = wu * (k_scale / (1.0 - u) ** 2) * k ** 2
    c, wc = leggauss(n_cos); st = np.sqrt(1.0 - c * c)
    phi = 2.0 * np.pi * np.arange(n_phi) / n_phi; wp = 2.0 * np.pi / n_phi
    kx = k[:, None, None] * (st[None, :, None] * np.cos(phi)[None, None, :])
    ky = k[:, None, None] * (st[None, :, None] * np.sin(phi)[None, None, :])
    kz = k[:, None, None] * (c[None, :, None] * np.ones_like(phi)[None, None, :])
    kv = np.stack([kx.ravel(), ky.ravel(), kz.ravel()], axis=-1)
    w = (wk[:, None, None] * wc[None, :, None]
         * np.ones_like(phi)[None, None, :] * wp).ravel()
    return kv, w


def _slater_density_ft(cenB, D, zA, zB, nt=64, h=1e-3):
    """Exact FT of the two-centre 1s Slater density (X at origin, other centre at
    cenB, D=X-cenB): e^{-zr}=-d/dz(e^{-zr}/r), the Yukawa-product FT is a 1D Feynman
    integral, rho~ = N_A N_B d^2/dzA dzB of it (central finite differences)."""
    NA = math.sqrt(zA ** 3 / np.pi); NB = math.sqrt(zB ** 3 / np.pi)
    Dn = float(np.linalg.norm(D))
    x, wx = leggauss(nt); t = 0.5 * (x + 1.0); wt = 0.5 * wx

    def T(kv, zx, zy):
        kdotB = kv @ cenB; kdotD = kv @ D; k2 = np.einsum("ij,ij->i", kv, kv)
        Delta = np.sqrt(np.outer(k2, t * (1 - t)) + (t * zx ** 2 + (1 - t) * zy ** 2)[None, :])
        phase = np.exp(1j * (kdotB[:, None] + np.outer(kdotD, 1 - t)))
        return 2 * np.pi * ((phase * np.exp(-Delta * Dn) / Delta) * wt[None, :]).sum(axis=1)

    def ft(kv):
        return NA * NB * (T(kv, zA + h, zB + h) - T(kv, zA + h, zB - h)
                          - T(kv, zA - h, zB + h) + T(kv, zA - h, zB - h)) / (4 * h * h)
    return ft


@pytest.mark.slow
def test_momentum_true_slater_density_ft():
    """(XY|XZ) via the TRUE 1s Slater density FT (Feynman reduction, not the
    Gaussian fit) reproduces the ground truth to 4.3e-7 (Paper 59 sec:momentum,
    C8 headline 1).  This validates the actual physical density, not the fitted
    Gaussians of test_momentum_reproduces_three_center_eri."""
    kv, w = _kspace_grid(160, 48, 48, 5.0)
    k2 = np.einsum("ij,ij->i", kv, kv)
    f1 = _slater_density_ft(Y, X - Y, 1.0, 1.0)(kv)   # X at origin
    f2 = _slater_density_ft(Z, X - Z, 1.0, 1.0)(kv)
    mom = (1.0 / (2.0 * np.pi ** 2)) * np.sum(w * (f1 * np.conj(f2)) / k2)
    assert abs(mom.imag) < 1e-9
    assert abs(mom.real - GROUND_TRUTH) < 1e-6, (mom.real, abs(mom.real - GROUND_TRUTH))


def test_momentum_angular_reduced_j0_form():
    """The fully angular-reduced form eq:reduced (angular integral done -> single
    j0(k|W|) kernel, W=sY-tZ; 2D Feynman + 1D radial) reproduces the ground truth
    to 1.8e-6 (Paper 59 sec:momentum, C8 headline 1).  Pins the reduction eq:reduced,
    not just the raw 3D momentum integral."""
    from itertools import product
    D1, D2 = float(np.linalg.norm(Y)), float(np.linalg.norm(Z))
    ns = nt = 32; nk = 200; k_scale = 5.0; h = 2e-3
    xs, ws = leggauss(ns); s = 0.5 * (xs + 1); ws = 0.5 * ws
    xt, wt = leggauss(nt); t = 0.5 * (xt + 1); wt = 0.5 * wt
    xk, wk0 = leggauss(nk); u = 0.5 * (xk + 1); wu = 0.5 * wk0
    k = k_scale * u / (1 - u); wk = wu * (k_scale / (1 - u) ** 2); k2 = k * k
    Sg, Tg = np.meshgrid(s, t, indexing="ij")
    W = Sg[..., None] * Y[None, None, :] - Tg[..., None] * Z[None, None, :]
    Wn = np.linalg.norm(W, axis=-1)
    kb = Wn[:, :, None] * k[None, None, :]
    jb = np.where(kb > 1e-12, np.sin(np.where(kb > 1e-12, kb, 1.0)) / np.where(kb > 1e-12, kb, 1.0), 1.0)

    def Jval(za, zb, zc, zd):
        D1a = np.sqrt(s[:, None] * (1 - s[:, None]) * k2[None, :]
                      + (s[:, None] * za ** 2 + (1 - s[:, None]) * zb ** 2))
        D2a = np.sqrt(t[:, None] * (1 - t[:, None]) * k2[None, :]
                      + (t[:, None] * zc ** 2 + (1 - t[:, None]) * zd ** 2))
        g1 = np.exp(-D1 * D1a) / D1a; g2 = np.exp(-D2 * D2a) / D2a
        Ik = (jb * g1[:, None, :] * g2[None, :, :] * wk[None, None, :]).sum(axis=2)
        return (Ik * ws[:, None] * wt[None, :]).sum()

    tot = 0.0
    for sa, sb, sc, sd in product([1, -1], repeat=4):
        tot += sa * sb * sc * sd * Jval(1 + sa * h, 1 + sb * h, 1 + sc * h, 1 + sd * h)
    mC = (8.0 / np.pi) * tot / (2 * h) ** 4
    assert abs(mC - GROUND_TRUTH) < 5e-6, (mC, abs(mC - GROUND_TRUTH))


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
        # measured worst residual ~1e-31 at dps=30; tolerance tightened from the
        # old 1e-10 to actually witness the quoted (paper eq:modpf) precision.
        assert abs(M_rho(lambda r: mp.ellipk(r), rho)) < mp.mpf(10) ** -25
        assert abs(M_rho(lambda r: mp.ellipk(1 - r), rho)) < mp.mpf(10) ** -25


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
# Cosmic-Galois pins (Paper 61 sec:modular; split out of Paper 59
# 2026-09-06): the elliptic family is the
# Legendre / Gamma(2) universal family, its CM-fiber periods are Gamma-values,
# and the integrated observable is a Gamma(2) multiple modular value.
# ---------------------------------------------------------------------------

def test_cosmic_galois_family_is_gamma2():
    """lambda(tau(rho)) = 1 - rho exactly (Paper 61 eq:lambda_rho): the family is the
    Legendre / Gamma(2) universal family, tau(rho)=i K(rho)/K(1-rho)."""
    import mpmath as mp
    mp.mp.dps = 30

    def lam_theta(tau):
        q = mp.e ** (1j * mp.pi * tau)
        return (mp.jtheta(2, 0, q) / mp.jtheta(3, 0, q)) ** 4

    # dps=50 so the test witnesses the ~1e-41 precision the paper (sec:modular)
    # quotes for lambda(tau(rho))=1-rho (worst residual ~2.7e-51 at dps=50; the
    # earlier dps=30 floored the assertion at 1e-20 while the driver reached 1e-41).
    mp.mp.dps = 50
    for rho in ['0.37', '0.6', '0.15']:
        r = mp.mpf(rho)
        tau = 1j * mp.ellipk(r) / mp.ellipk(1 - r)
        lam = mp.re(lam_theta(tau))
        assert abs(lam - (1 - r)) < mp.mpf(10) ** -41, (rho, lam)


def test_cosmic_galois_cm_periods_are_gamma_values():
    """CM-fiber periods are Gamma-values (Chowla-Selberg), at two fundamental
    discriminants (Paper 61 sec:modular): disc -4 (tau=i, rho=1/2, in the physical
    domain) and disc -8 (tau=i sqrt2)."""
    import mpmath as mp
    # dps=55 so the test witnesses the ~1e-51 precision quoted in the driver
    # (at dps=40 the disc-8 residual floors at ~4.6e-41).
    mp.mp.dps = 55
    # disc -4
    K4 = mp.ellipk(mp.mpf('0.5'))
    cf4 = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(mp.pi))
    assert abs(K4 - cf4) < mp.mpf(10) ** -50
    # disc -8: lambda_2 = 3 - 2 sqrt2
    lam2 = 3 - 2 * mp.sqrt(2)
    K8 = mp.ellipk(lam2)
    cf8 = (mp.sqrt(1 + mp.sqrt(2)) * mp.gamma(mp.mpf(1) / 8) * mp.gamma(mp.mpf(3) / 8)
           / (2 ** mp.mpf('3.25') * mp.sqrt(mp.pi)))
    assert abs(K8 - cf8) < mp.mpf(10) ** -50


def test_cosmic_galois_integrated_value():
    """The integrated collinear T2 = 0.39535576590171392... (Paper 61 sec:modular),
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


# ---------------------------------------------------------------------------
# Differential-Galois reducibility scan of the rank-4 Picard-Fuchs operator
# L4 (Paper 59 eq:pf, in the Laplace variable D; rho = c2/c1):
#   D rho L'''' + 2 rho L''' + D(1-2rho) L'' + (1-2rho) L' - D(1-rho) L = 0.
# Backs the sharpened closed-form open problem (Paper 59 sec:obstruction/sec:scope):
# L4 is self-adjoint (Galois group in Sp4) and admits no low-order factor over
# Q(rho)(D) with polar locus at the true singularities {0,inf}.
# Memo: debug/sprint_L4_reducibility_scan_memo.md.
# ---------------------------------------------------------------------------

def _L4_coeffs(D, rho):
    """L4 = sum p[k] d^k/dD^k (non-monic, original eq:pf form)."""
    return {0: -D * (1 - rho), 1: (1 - 2 * rho), 2: D * (1 - 2 * rho),
            3: 2 * rho, 4: D * rho}


def test_L4_is_self_adjoint():
    """L4* = L4 exactly (formally self-adjoint) => differential Galois group in Sp4,
    and left factors coincide with right factors."""
    import sympy as sp
    D, rho = sp.symbols('D rho')
    y = sp.Function('y')(D)
    p = _L4_coeffs(D, rho)
    Ly = sum(p[k] * sp.diff(y, D, k) for k in p)
    Lstar = sum((-1) ** k * sp.diff(p[k] * y, D, k) for k in p)   # formal adjoint
    assert sp.simplify(sp.expand(Ly - Lstar)) == 0


def test_L4_indicial_exponents():
    """Indicial polynomial of L4 at the regular-singular point D=0 is
    rho*s*(s-2)*(s-1)^2 => exponents {0,1,1,2}.  The DOUBLE exponent at s=1 is the
    local-monodromy origin of the D ln D nonanalyticity of N(D)."""
    import sympy as sp
    from collections import defaultdict
    D, s, rho = sp.symbols('D s rho')
    p = _L4_coeffs(D, rho)
    bucket = defaultdict(lambda: sp.Integer(0))
    for k, c in p.items():
        for (deg,), a in sp.Poly(c, D).terms():          # c = sum a D^deg
            bucket[deg - k] += a * sp.ff(s, k)           # term a*ff(s,k)*D^(s+deg-k)
    ind = sp.factor(sp.simplify(bucket[min(bucket)]))
    assert sp.simplify(ind - rho * s * (s - 2) * (s - 1) ** 2) == 0
    roots = sp.roots(sp.Poly(ind, s))                    # {root: multiplicity}
    assert roots == {sp.Integer(0): 1, sp.Integer(1): 2, sp.Integer(2): 1}


def test_L4_symbol_only_two_bessel_sectors():
    """The symbol s^4 + ((1-2rho)/rho)s^2 - (1-rho)/rho factors over Q(rho) ONLY as
    (s^2-1)(s^2+(1-rho)/rho) -- the K0 and J0 Bessel sectors; no mixed split (its
    constant term is negative).  So the only candidate order-2 factorizations are
    those two, both excluded by test_L4_no_order2_rational_factor."""
    import sympy as sp
    s, rho = sp.symbols('s rho')
    symbol = s ** 4 + ((1 - 2 * rho) / rho) * s ** 2 - (1 - rho) / rho
    assert sp.expand(symbol - (s ** 2 - 1) * (s ** 2 + (1 - rho) / rho)) == 0
    fac = sp.factor(symbol)                              # over Q(rho)
    # irreducible-over-Q(rho) quadratic pieces: (s-1),(s+1),(rho s^2 - rho + 1)
    assert sp.simplify(fac * rho - (s - 1) * (s + 1) * (rho * s ** 2 - rho + 1)) == 0


def test_L4_irregular_at_infinity():
    """L4 has an IRREGULAR singularity at D=inf (Poincare rank 1).  Substituting the
    formal exponential ansatz y=e^{lam D} and collecting in D, the leading (top) edge
    sits at D^1 and its coefficient is the degree-4 polynomial
    rho*(lam^2-1)*(lam^2+(1-rho)/rho), whose four roots lam={+-1, +-i*sqrt((1-rho)/rho)}
    are all NONZERO (the K0/I0 and J0/Y0 Bessel sectors).  A nonzero constant lam is the
    defining signature of an irregular singularity.  Contrast: the sunrise/Legendre
    elliptic-period operator M_rho = rho(1-rho)d^2 + (1-2rho)d - 1/4 is FUCHSIAN -- its
    leading edge forces lam=0.  This is the operator-level form of "N(D) is the irregular
    Laplace dual of the (Fuchsian) unequal-mass sunrise", so the elliptic-polylog
    machinery that closes the sunrise (Bogner-Mueller-Stach-Weinzierl, arXiv:1907.01251)
    cannot reach N(D).  Driver: debug/routeC_L4_irregularity.py."""
    import sympy as sp
    D, rho, lam = sp.symbols('D rho lambda')
    # L4: leading exponential balance at D=inf
    p = _L4_coeffs(D, rho)
    expr = sum(c * lam ** k for k, c in p.items())         # e^{lam D} factored out
    poly = sp.Poly(sp.expand(expr), D)
    top = poly.degree()
    lead = sp.expand(poly.nth(top))
    assert top == 1                                        # full-width slope-1 edge
    assert sp.simplify(lead - rho * (lam ** 2 - 1) * (lam ** 2 + (1 - rho) / rho)) == 0
    roots = sp.roots(sp.Poly(lead, lam))
    assert sum(roots.values()) == 4                        # all four exponents accounted for
    assert all(r != 0 for r in roots)                      # every exponent nonzero => irregular
    # Fuchsian contrast: the Legendre elliptic-period operator forces lam = 0 only
    leg = {2: rho * (1 - rho), 1: (1 - 2 * rho), 0: sp.Rational(-1, 4)}
    lexpr = sum(c * lam ** k for k, c in leg.items())      # ODE variable is rho here
    lpoly = sp.Poly(sp.expand(lexpr), rho)
    llead = sp.expand(lpoly.nth(lpoly.degree()))
    assert sp.roots(sp.Poly(llead, lam)) == {sp.Integer(0): 2}


def test_T2_modular_pullback():
    """Backs Paper 61 sec:modular: the fibre period is the Gamma(2) modular quantity
    K(m) = (pi/2) theta3(0,q)^2 at tau = i K(1-m)/K(m), with lambda(tau) = (theta2/theta3)^4 = m.
    Validated to 40 digits -- the tau<->m pullback the Lambert-series route rests on.
    Driver: (modular foundation for the T2 hand-off)."""
    import mpmath as mp
    mp.mp.dps = 40
    for m in [mp.mpf('0.1'), mp.mpf('0.3'), mp.mpf('0.5'), mp.mpf('0.7')]:
        K = mp.ellipk(m); Kp = mp.ellipk(1-m)
        q = mp.exp(-mp.pi*Kp/K)                        # nome, q = e^{i pi tau}, tau = i Kp/K
        th2 = mp.jtheta(2, 0, q); th3 = mp.jtheta(3, 0, q)
        assert abs((mp.pi/2)*th3**2 - K) < mp.mpf('1e-35')      # K(m) = (pi/2) theta3^2
        assert abs((th2/th3)**4 - m)    < mp.mpf('1e-35')      # lambda(tau) = m


def test_T2_fiber_spectral():
    """Corrects Paper 61 sec:modular: the ~20-digit T2 ceiling was a FIXED fiber-grid artifact,
    NOT an intrinsic k-grid limit.  The fiber J(s,t)=int_0^inf j0(k(s+t)) P(s,k)P(t,k) dk, on the
    decay-scaled map k=L u/(1-u), L=1/(sqrt(c_s)+sqrt(c_t)), converges SPECTRALLY: a fixed
    Gauss-Legendre fiber reaches 30+ digits by Nk~120.  So the fiber is not the precision
    bottleneck (the outer-integral GL rate is).  Driver: debug/routeC_T2_highprec.py."""
    import mpmath as mp
    mp.mp.dps = 40
    def Pw(x, k):
        c = x*(1-x); D = mp.sqrt(c*k*k+1)
        return c*mp.e**(-D)*(1/D**3+3/D**4+3/D**5)
    def gl(N):
        from mpmath import legendre
        r=[]; w=[]
        for kk in range(1, N+1):
            x = mp.cos(mp.pi*(kk-mp.mpf('0.25'))/(N+mp.mpf('0.5')))
            for _ in range(100):
                f=legendre(N,x); fp=N*(x*legendre(N,x)-legendre(N-1,x))/(x*x-1); dx=f/fp; x-=dx
                if abs(dx) < mp.mpf(10)**(-mp.mp.dps-6): break
            r.append(x)
        for x in r:
            fp=N*(x*legendre(N,x)-legendre(N-1,x))/(x*x-1); w.append(2/((1-x*x)*fp*fp))
        return r, w
    def J(s, t, Nk):
        cs=s*(1-s); ct=t*(1-t); L=1/(mp.sqrt(cs)+mp.sqrt(ct)); b=s+t
        xs, ws = gl(Nk); tot=mp.mpf(0)
        for x, w in zip(xs, ws):
            u=(x+1)/2; k=L*u/(1-u); dk=L/(1-u)**2
            j0 = mp.sin(k*b)/(k*b) if k*b > mp.mpf('1e-40') else mp.mpf(1)
            tot += (w/2)*j0*Pw(s,k)*Pw(t,k)*dk
        return tot
    s = t = mp.mpf('0.5')                                   # bulk point (slowest fiber)
    ref = mp.quad(lambda k: (mp.sin(k*(s+t))/(k*(s+t)))*Pw(s,k)*Pw(t,k), [0,1,3,8,20,mp.inf])
    e60  = abs(J(s,t,60)  - ref)
    e120 = abs(J(s,t,120) - ref)
    assert e60  < mp.mpf('1e-15')                           # already deep at moderate order
    assert e120 < mp.mpf('1e-20')                           # geometric gain => SPECTRAL fiber


def test_N_stokes_constants_algebraic():
    """The build result: N(D)'s resurgence closes with ALGEBRAIC Stokes constants.  The Borel
    transform psi(z)=[(z+2)(rho(1+z)^2+1-rho)]^{-1/2} has sqrt singularities with algebraic
    amplitudes a_*(-2)=1 and a_*(-1+-iw)^2 = -1/2 -+ (i/2) sqrt(rho/(1-rho)); the elliptic period
    enters only at the boundary: sqrt(c1) N(0) = L(0,rho) = K(1-rho) exactly.  So N(D) is the
    Borel-Laplace transform of an algebraic differential with algebraic Stokes data and an
    elliptic-period normalization -- the closed form at the irregular level.
    Driver: debug/routeC_stokes_constants.py."""
    import mpmath as mp
    mp.mp.dps = 40
    rho = mp.mpf(1)/5; w = mp.sqrt((1 - rho)/rho)
    psi = lambda z: 1/mp.sqrt((z + 2)*(rho*(1 + z)**2 + (1 - rho)))
    a2 = mp.limit(lambda h: psi(mp.mpf(-2) + h)*mp.sqrt(h), 0)          # real-sector amplitude
    assert abs(a2 - 1) < mp.mpf('1e-30')
    ac = mp.limit(lambda h: psi((-1 + 1j*w) + h)*mp.sqrt(h), 0)        # complex-sector amplitude
    assert abs(ac**2 - (mp.mpf(-1)/2 - (1j/2)*mp.sqrt(rho/(1 - rho)))) < mp.mpf('1e-30')
    L0 = mp.quad(lambda x: 1/mp.sqrt((x**2 - 1)*(rho*x**2 + (1 - rho))), [1, mp.inf])
    assert abs(L0 - mp.ellipk(1 - rho)) < mp.mpf('1e-20')             # D=0 boundary = elliptic period


def test_N_resurgent_series_and_borel():
    """N(D) (eq:laplace) about the Bessel point rho=0 has the explicit ELEMENTARY series
    N=(1/sqrt c1) sum_n (-1)^n ((2n-1)!!)^2/(2^n n!) rho^n K_n(D)/D^n; it reproduces the direct
    integral at optimal truncation (a few digits, then diverges -- asymptotic).  Its Borel
    transform psi(zeta)=[(zeta+2)(rho(1+zeta)^2+1-rho)]^{-1/2} is singular at the 3 non-dominant
    branch points {-2, -1 +- i sqrt((1-rho)/rho)}, confirmed by Borel-Pade.  Function-level
    complement to test_L4_irregular_at_infinity.  Drivers: debug/routeC_Nborel.py,
    debug/routeC_Ndiagonal_expansion.py."""
    import mpmath as mp
    mp.mp.dps = 30
    rho = mp.mpf(1)/5; D = mp.mpf(3)   # D=3: optimal-truncation err ~8e-5 (comfortably a few digits)
    def a(n):
        d = mp.mpf(1)
        for k in range(1, 2*n, 2): d *= k
        return (-1)**n * d**2 / (mp.mpf(2)**n * mp.factorial(n))
    L = mp.quad(lambda th: mp.e**(-D*mp.cosh(th))/mp.sqrt(rho*mp.cosh(th)**2+1-rho), [0, 8])
    S = mp.mpf(0); best = None
    for n in range(0, 12):
        S += a(n) * rho**n * mp.besselk(n, D) / D**n
        e = abs(S - L); best = e if best is None else min(best, e)
    assert best < mp.mpf('1e-3')                 # asymptotic series tracks N to a few digits
    w = mp.sqrt((1 - rho)/rho)                    # = 2 at rho=1/5
    coeffs = mp.taylor(lambda x: 1/mp.sqrt((x+2)*(rho*(1+x)**2+(1-rho))), 0, 26)
    _, q = mp.pade(coeffs, 12, 12)
    poles = sorted(mp.polyroots(q[::-1], maxsteps=200, extraprec=200), key=abs)
    assert min(abs(z + 2) for z in poles) < 0.05                 # leading real pole ~ -2
    assert min(abs(z - (-1 + 1j*w)) for z in poles) < 0.05       # leading complex pole ~ -1+2i


def test_L4_no_order1_hyperexponential_factor():
    """No order-1 right factor: no solution y = e^{lam D} Q(D) (Q polynomial) for
    lam in {+1,-1,+iw,-iw,0}, w=2 at rho=1/5.  Since the D=0 exponents are integers,
    this is the COMPLETE order-1-factor test for polar locus {0,inf}."""
    import sympy as sp
    D = sp.symbols('D')
    rho = sp.Rational(1, 5)
    p = _L4_coeffs(D, rho)

    def has_solution(lam, Nmax=6):
        for N in range(Nmax + 1):
            qs = sp.symbols(f'q0:{N + 1}')
            Q = sum(qs[j] * D ** j for j in range(N + 1))
            yv = sp.exp(lam * D) * Q
            e = sp.expand(sum(p[k] * sp.diff(yv, D, k) for k in p) / sp.exp(lam * D))
            sold = sp.solve(sp.Poly(sp.simplify(e), D).all_coeffs(), qs, dict=True)
            if sold and any(sold[0].get(q, q) != 0 for q in qs):
                return True
        return False

    w = sp.Integer(2)                                    # w^2=(1-rho)/rho=4
    for lam in [sp.Integer(1), sp.Integer(-1), sp.I * w, -sp.I * w, sp.Integer(0)]:
        assert not has_solution(lam), f"unexpected hyperexponential solution at lam={lam}"


def test_L4_no_order2_rational_factor_with_positive_control():
    """No order-2 right factor B2=d^2+a d+b over Q(rho)(D) with polar locus {0,inf}
    (a: simple pole at 0; b: double pole at 0; bounded at inf).  The identical
    pipeline RECOVERS a planted factor (positive control), so the negative is real."""
    import sympy as sp
    D = sp.symbols('D')
    am1, a0, a1, bm2, bm1, b0, b1 = sp.symbols('am1 a0 a1 bm2 bm1 b0 b1')
    unks = [am1, a0, a1, bm2, bm1, b0, b1]

    def remainder(p3, p2, p1, p0, a, b):
        ap, bp = sp.diff(a, D), sp.diff(b, D)
        P = a ** 2 - ap - b
        Q0 = a * b - bp
        R1 = (-a * P + sp.diff(P, D) + Q0) + p3 * (a ** 2 - ap - b) + p2 * (-a) + p1
        R0 = (sp.diff(Q0, D) - b * P) + p3 * (a * b - bp) + p2 * (-b) + p0
        return R1, R0

    def search(p3, p2, p1, p0):
        a = am1 / D + a0 + a1 * D
        b = bm2 / D ** 2 + bm1 / D + b0 + b1 * D
        eqs = []
        for R in remainder(p3, p2, p1, p0, a, b):
            num = sp.numer(sp.cancel(sp.together(R)))
            eqs.extend(sp.Poly(sp.expand(num), D).all_coeffs())
        eqs = [sp.expand(e) for e in eqs if e != 0]
        return sp.solve(eqs, unks, dict=True)

    # POSITIVE CONTROL: L = (d^2+(2/D)d+3) o (d^2+(1/D)d-1) has p=(3/D,2,1/D,-3);
    # the pipeline must recover the planted right factor (am1=1,b0=-1, rest 0).
    ctrl = search(3 / D, sp.Integer(2), 1 / D, sp.Integer(-3))
    assert any(s.get(b0) == -1 and s.get(am1) == 1 for s in ctrl), ctrl

    # L4 at rho=1/5, monic: p3=2/D, p2=(1-2rho)/rho=3, p1=3/D, p0=-(1-rho)/rho=-4.
    assert search(2 / D, sp.Integer(3), 3 / D, sp.Integer(-4)) == []


def test_L4_eigenring_is_trivial():
    """The eigenring E(L4) = {R : ord R<4, L R == 0 mod L} = C (scalars only), so L4 is
    INDECOMPOSABLE.  With self-adjointness (Sp4) and no hyperexponential solution, this
    pins reducibility to a Lagrangian (self-dual order-2) factor.  Solution space of a
    solution y is its jet (y,y',y'',y'''); Y4 = -(p0 Y0+p1 Y1+p2 Y2+p3 Y3) [monic]."""
    import sympy as sp
    D = sp.symbols('D')
    rho = sp.Rational(1, 5)
    p0, p1, p2, p3 = -(1 - rho) / rho, (1 - 2 * rho) / (rho * D), (1 - 2 * rho) / rho, 2 / D

    def deriv(v):
        a0, a1, a2, a3 = v
        return [sp.diff(a0, D) - a3 * p0, a0 + sp.diff(a1, D) - a3 * p1,
                a1 + sp.diff(a2, D) - a3 * p2, a2 + sp.diff(a3, D) - a3 * p3]

    def Lapply(v):
        Dm = [v]
        for _ in range(4):
            Dm.append(deriv(Dm[-1]))
        cf = [p0, p1, p2, p3, sp.Integer(1)]            # D0..D4, monic
        return [sum(cf[m] * Dm[m][i] for m in range(5)) for i in range(4)]

    # self-check: R = identity (r=[1,0,0,0]) => L[y] = 0
    assert all(sp.simplify(e) == 0 for e in Lapply([sp.Integer(1), 0, 0, 0]))

    # ansatz r_i = sum_{j=-2}^{2} c_ij D^j ; count free params in the solution space
    cs = {i: sp.symbols(f'c{i}_a c{i}_b c{i}_c c{i}_d c{i}_e') for i in range(4)}
    r = [sum(cs[i][j] * D ** (j - 2) for j in range(5)) for i in range(4)]
    allc = [c for i in range(4) for c in cs[i]]
    eqs = []
    for Ei in Lapply(r):
        eqs.extend(sp.Poly(sp.expand(sp.numer(sp.cancel(sp.together(Ei)))), D).all_coeffs())
    sol = list(sp.linsolve([sp.expand(e) for e in eqs if e != 0], allc))[0]
    free = set().union(*[sp.sympify(s).free_symbols for s in sol]) & set(allc)
    subs = dict(zip(allc, sol))
    nonzero = [i for i in range(4)
               if sp.simplify(sum(subs[cs[i][j]] * D ** (j - 2) for j in range(5))) != 0]
    assert len(free) == 1 and nonzero == [0], (len(free), nonzero)   # E(L4) = C, only r0


# Monodromy of L4 around D=0, in the exponential-line (Lefschetz-thimble) basis, is the
# integer unipotent matrix below -- computed numerically (below / debug/routeC_L4_reducibility.py),
# IDENTICAL for rho=1/5 and rho=1/3 (a topological invariant of the 4-branch-point config of
# y^2=Q).  It closes the base-field reducibility question: no invariant coordinate subspace +
# the exponential torus => L4 is IRREDUCIBLE over C(D), hence over Q(rho)(D).
_L4_MONODROMY = [[-1, 2, 2, 2], [-2, 3, 2, 2], [-2, 2, 3, 2], [2, -2, -2, -1]]


def test_L4_monodromy_certificate_irreducible():
    """Exact certificate from the computed monodromy M0 (exp-line basis): unipotent with a
    single 2x2 Jordan block (the one D ln D log), and NO invariant coordinate subspace.  With
    the exponential torus (invariant subspaces must be coordinate subspaces) => L4 IRREDUCIBLE
    over C(D) => over Q(rho)(D).  Numerically re-derived in test_L4_monodromy_numerically."""
    import sympy as sp
    from itertools import combinations
    M = sp.Matrix(_L4_MONODROMY)
    I4 = sp.eye(4)
    assert M.det() == 1 and M.trace() == 4
    assert M.eigenvals() == {sp.Integer(1): 4}                       # unipotent
    assert (M - I4) ** 2 == sp.zeros(4) and (M - I4).rank() == 1     # single 2x2 Jordan block
    # torus forces invariant subspaces to be coordinate subspaces; none is invariant:
    def invariant(S):
        other = [i for i in range(4) if i not in S]
        return all(M[c, a] == 0 for a in S for c in other)
    assert not any(invariant(S) for d in (1, 2, 3) for S in combinations(range(4), d))


def _l4_monodromy_numeric(rho_str):
    """Numerically re-derive the L4 monodromy around D=0 in the thimble basis for
    a given rho, returning (err, Mint).  Thimble solution y=int_{x_k}^{...}
    e^{-Dx}dx/sqrt(Q), monodromy by RK4 of the companion system around |D|=|D0|;
    C = thimble jets at D0.  Parametrized by rho: w2=(1-rho)/rho sets the
    imaginary branch points +-i sqrt(w2), c2=(1-2rho)/rho sets the monic
    coefficients (p3=2/D, p2=c2, p1=c2/D, p0=-w2)."""
    import mpmath as mp
    mpf = mp.mpf
    rho = mpf(rho_str)
    w2 = (1 - rho) / rho                                     # (1-rho)/rho
    c2 = (1 - 2 * rho) / rho                                 # (1-2rho)/rho
    wr = mp.sqrt(w2)
    BP = [mpf(1), mpf(-1), 1j * wr, -1j * wr]
    Q = lambda x: (x * x - 1) * (x * x + w2)
    Qp = lambda x: 2 * x * (2 * x * x + (w2 - 1))
    Nf = 1500

    def thimble(xk, D):
        d = mp.conj(D) / abs(D); umax = mp.sqrt(45 / abs(D)) + 3; h = umax / Nf
        ref = mp.arg(Qp(xk) * d); pa = ref; vals = []
        for i in range(Nf + 1):
            u = i * h; x = xk + u * u * d
            if i == 0:
                base = 2 * d * mp.e ** (-D * xk) / (mp.sqrt(abs(Qp(xk) * d)) * mp.e ** (1j * ref / 2))
            else:
                qv = Q(x); a = mp.arg(qv)
                while a - pa > mp.pi: a -= 2 * mp.pi
                while a - pa < -mp.pi: a += 2 * mp.pi
                pa = a
                base = mp.e ** (-D * x) / (mp.sqrt(abs(qv)) * mp.e ** (1j * a / 2)) * (2 * u) * d
            vals.append([((-x) ** n) * base for n in range(4)])
        out = []
        for n in range(4):
            s = vals[0][n] + vals[Nf][n]
            for i in range(1, Nf): s += (4 if i % 2 else 2) * vals[i][n]
            out.append(s * h / 3)
        return out

    # monic L4: y'''' = w2*y - (c2/D)y' - c2*y'' - (2/D)y'''
    Am = lambda D: mp.matrix([[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
                              [w2, -c2 / D, -c2, -2 / D]])
    D0 = 3 * mp.e ** (1j * mpf('0.5')); R = abs(D0); th0 = mp.arg(D0)
    Y = mp.eye(4); Ns = 1500; dt = 2 * mp.pi / Ns
    Fn = lambda t, Y: (1j * R * mp.e ** (1j * (th0 + t))) * (Am(R * mp.e ** (1j * (th0 + t))) * Y)
    t = mpf(0)
    for _ in range(Ns):
        k1 = Fn(t, Y); k2 = Fn(t + dt / 2, Y + (dt / 2) * k1)
        k3 = Fn(t + dt / 2, Y + (dt / 2) * k2); k4 = Fn(t + dt, Y + dt * k3)
        Y = Y + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4); t += dt
    C = mp.matrix(4, 4)
    for k, xk in enumerate(BP):
        j = thimble(xk, D0)
        for n in range(4): C[n, k] = j[n]
    Me = C ** -1 * Y * C
    err = max(abs(Me[i, j] - mp.nint(mp.re(Me[i, j]))) for i in range(4) for j in range(4))
    Mint = [[int(mp.nint(mp.re(Me[i, j]))) for j in range(4)] for i in range(4)]
    return err, Mint


@pytest.mark.slow
def test_L4_monodromy_numerically():
    """Numerically re-derive the L4 monodromy around D=0 in the thimble basis and
    check it rounds to the SAME integer matrix _L4_MONODROMY for BOTH rho=1/5 AND
    rho=1/3 -- the topological-invariant claim (Paper 59 sec:obstruction) that
    the monodromy is independent of the modulus, not just asserted in a comment."""
    import mpmath as mp
    mp.mp.dps = 20
    for rho_str in ('0.2', '0.3333333333333333333'):           # rho = 1/5 and 1/3
        err, Mint = _l4_monodromy_numeric(rho_str)
        assert err < mp.mpf(10) ** -4, (rho_str, err)          # unambiguously integer
        assert Mint == _L4_MONODROMY, (rho_str, Mint)


# ---------------------------------------------------------------------------
# Intersection-form theorem (Paper 61 sec:bessel_algebra) -- closes the
# "one step short of a theorem": the concomitant B = pi x (integer intersection
# form) is now the CANONICAL block symplectic Omega, forced by the monodromy +
# the concomitant's sector-block structure.  Driver: debug/routeC_intersection_form.py.
# ---------------------------------------------------------------------------
_OMEGA = [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]]


def test_intersection_form_is_forced_canonical_symplectic():
    """Symbolic closure of sec:bessel_algebra.  Among antisymmetric forms preserved by the
    exact-integer monodromy M0 (_L4_MONODROMY, thimble basis) -- a 4-parameter family --
    the ones with the concomitant's vanishing cross-sector pairings (real sector {0,1} vs
    imaginary {2,3}) are EXACTLY Z*Omega, Omega the canonical block symplectic form.  Omega
    is unimodular (det=1, nondegenerate) and monodromy-preserved (M0^T Omega M0 = Omega =>
    MONODROMY group in Sp(Omega,Z)=Sp4(Z) -- NOT the differential Galois
    group, which is Zariski-closed and cannot lie in a discrete group; it is in
    Sp4(C) by self-adjointness. Corrected 2026-09-07).  The period-cut {K,I,J} sub-block is
    rank 2; the fourth (Y0) thimble completes it to nondegenerate rank 4."""
    import sympy as sp
    M0 = sp.Matrix(_L4_MONODROMY)
    Om = sp.Matrix(_OMEGA)
    # (a) Omega is a nondegenerate integer symplectic form preserved by M0
    assert Om.T == -Om and Om.det() == 1
    assert M0.T * Om * M0 == Om                       # MONODROMY in Sp(Omega,Z)
    # (b) it is FORCED: solve M0^T J M0 = J for antisymmetric J -> 4-parameter family,
    #     then impose the concomitant's block-diagonality -> exactly Z*Omega
    a, b, c, d, e, f = sp.symbols('a b c d e f')
    J = sp.Matrix([[0, a, b, c], [-a, 0, d, e], [-b, -d, 0, f], [-c, -e, -f, 0]])
    fam = sp.solve([(M0.T * J * M0 - J)[i, j] for i in range(4) for j in range(4)],
                   [a, b, c, d, e, f], dict=True)
    assert len(fam) == 1
    Jp = J.subs(fam[0])
    assert len(Jp.free_symbols) == 4                            # 4-parameter preserved family
    block = sp.solve([Jp[0, 2], Jp[0, 3], Jp[1, 2], Jp[1, 3]], list(Jp.free_symbols), dict=True)
    Jb = sp.simplify(Jp.subs(block[0]))
    scale = list(Jb.free_symbols)
    assert len(scale) == 1                                      # exactly a 1-parameter Z*Omega
    assert sp.simplify(Jb.subs({scale[0]: 1})) == Om
    # (c) the period-cut {K,I,J} block is rank 2 (odd-dim antisym; the Y0 thimble completes it)
    assert sp.Matrix([[0, -1, 0], [1, 0, 2], [0, -2, 0]]).rank() == 2


def test_intersection_form_period_cut_values():
    """Numeric backing for B[K,I]=-pi, B[K,J]=0, B[I,J]=2pi -- the concomitant of the L4
    masters over the period cuts {[1,inf),[-1,1], i[-w,w]}, i.e. the paper's measured
    entries of pi*Omega, reproduced here at rho=1/2.  These are the (integer) intersection
    numbers -1, 0, 2 in the period-cut basis."""
    import mpmath as mp
    mp.mp.dps = 25
    rho = mp.mpf(1) / 2
    D = mp.mpf(1)
    w = mp.sqrt((1 - rho) / rho)

    def conc(y, z):
        a2 = D * rho; a1 = D * (1 - 2 * rho)
        return (z[0] * rho * (y[2] + D * y[3]) - z[1] * (a2 * y[2])
                - y[0] * rho * (z[2] + D * z[3]) + y[1] * (a2 * z[2])
                + a1 * (z[0] * y[1] - y[0] * z[1]))

    K = [mp.quad(lambda x: (-x) ** k * mp.e ** (-D * x) / mp.sqrt((x * x - 1) * (rho * x * x + 1 - rho)),
                 [1, mp.mpf('1.02'), mp.mpf('1.3'), 2, 4, 8, 16, mp.inf]) for k in range(4)]
    I = [mp.quad(lambda x: (-x) ** k * mp.e ** (-D * x) / mp.sqrt((1 - x * x) * (rho * x * x + 1 - rho)),
                 [-1, mp.mpf('-0.5'), 0, mp.mpf('0.5'), 1]) for k in range(4)]

    def Jrow(k):
        return mp.quad(lambda u: (((-1j * u) ** k * mp.e ** (-1j * D * u)
                       / mp.sqrt((u * u + 1) * ((1 - rho) - rho * u * u))).real), [-w, 0, w])
    Jm = [Jrow(k) for k in range(4)]
    bki = conc(K, I) / mp.pi
    bkj = conc(K, Jm) / mp.pi
    bij = conc(I, Jm) / mp.pi
    assert abs(bki - (-1)) < mp.mpf(10) ** -12, bki
    assert abs(bkj) < mp.mpf(10) ** -10, bkj
    assert abs(bij - 2) < mp.mpf(10) ** -12, bij


# ---- self-contained ports of the intersection-form driver (was debug/routeC_intersection_form.py;
# inlined per the transient-dir policy so the cited [SYMBOLIC] backing runs from a clean checkout) ----
def _if_thimble(xk, Dc, rho, Nf=6000):
    """Thimble masters [s(0),s(1),s(2),s(3)](D), s(n)=int_g (-x)^n e^{-D x} dx/sqrt(Q) from xk."""
    import mpmath as mp
    w2 = (1 - rho) / rho
    Qf = lambda x: (x * x - 1) * (x * x + w2)
    Qpf = lambda x: 2 * x * (2 * x * x + (w2 - 1))
    d = mp.conj(Dc) / abs(Dc)
    h = (mp.sqrt(60 / abs(Dc)) + 3) / Nf
    ref = mp.arg(Qpf(xk) * d); pa = ref; vals = []
    for i in range(Nf + 1):
        u = i * h; x = xk + u * u * d
        if i == 0:
            base = 2 * d * mp.e ** (-Dc * xk) / (mp.sqrt(abs(Qpf(xk) * d)) * mp.e ** (1j * ref / 2))
        else:
            qv = Qf(x); a = mp.arg(qv)
            while a - pa > mp.pi: a -= 2 * mp.pi
            while a - pa < -mp.pi: a += 2 * mp.pi
            pa = a
            base = mp.e ** (-Dc * x) / (mp.sqrt(abs(qv)) * mp.e ** (1j * a / 2)) * (2 * u) * d
        vals.append([((-x) ** n) * base for n in range(4)])
    out = []
    for n in range(4):
        ssum = vals[0][n] + vals[Nf][n]
        for i in range(1, Nf):
            ssum += (4 if i % 2 else 2) * vals[i][n]
        out.append(ssum * h / 3)
    return out


def _if_concomitant(y, z, D, rho):
    """Lagrange bilinear concomitant of the self-adjoint L4 (eq:pf); constant on solutions."""
    Ly, Ly1, Ly2, Ly3 = y; Lz, Lz1, Lz2, Lz3 = z
    a2 = D * rho; a1 = D * (1 - 2 * rho)
    return (Lz * rho * (Ly2 + D * Ly3) - Lz1 * (a2 * Ly2)
            - Ly * rho * (Lz2 + D * Lz3) + Ly1 * (a2 * Lz2)
            + a1 * (Lz * Ly1 - Ly * Lz1))


def _if_thimble_form(rho, Nf=6000):
    """Full 4x4 concomitant matrix B/pi in the thimble basis."""
    import mpmath as mp
    Dc = 3 * mp.e ** (1j * mp.mpf('0.5'))
    w = mp.sqrt((1 - rho) / rho)
    S = [_if_thimble(xk, Dc, rho, Nf) for xk in (mp.mpf(1), mp.mpf(-1), 1j * w, -1j * w)]
    B = mp.matrix(4, 4)
    for a in range(4):
        for b in range(4):
            B[a, b] = _if_concomitant(S[a], S[b], Dc, rho) / mp.pi
    return B


def _if_branch_point_proof():
    """SYMBOLIC block-diagonality + plane-equality + the pi from branch-point Laplace asymptotics.
    Returns (ident, Bconst, Bsq): the identity Q'(x)Q'(-x)=-4x^2(2rho x^2-2rho+1)^2, the within-
    sector leading value, and B[s_x,s_-x]^2 (= -pi^2, x- and rho-independent)."""
    import sympy as sp
    D, x, rho = sp.symbols('D x rho'); half = sp.Rational(1, 2)
    A, Bc = sp.symbols('A B')
    y = A * sp.exp(-x * D) * D ** (-half); z = Bc * sp.exp(+x * D) * D ** (-half)
    def dd(f, n):
        for _ in range(n): f = sp.diff(f, D)
        return f
    Y = [dd(y, n) for n in range(4)]; Z = [dd(z, n) for n in range(4)]
    a2 = D * rho; a1 = D * (1 - 2 * rho)
    Bform = (Z[0] * rho * (Y[2] + D * Y[3]) - Z[1] * (a2 * Y[2])
             - Y[0] * rho * (Z[2] + D * Z[3]) + Y[1] * (a2 * Z[2])
             + a1 * (Z[0] * Y[1] - Y[0] * Z[1]))
    Bconst = sp.simplify(sp.limit(sp.simplify(Bform), D, sp.oo))
    Q = (x ** 2 - 1) * (rho * x ** 2 + 1 - rho); Qp = sp.diff(Q, x)
    g = 2 * rho * x ** 2 - 2 * rho + 1
    ident = sp.simplify(sp.expand(Qp * Qp.subs(x, -x)) - sp.expand(-4 * x ** 2 * g ** 2)) == 0
    AB2 = sp.pi ** 2 / sp.expand(Qp * Qp.subs(x, -x))
    Bsq = sp.simplify((2 * x * (-g)) ** 2 * AB2)
    return ident, sp.simplify(Bconst), Bsq


@pytest.mark.slow
def test_intersection_form_thimble_block_structure():
    """The concomitant is block-diagonal in the thimble basis and equals pi*(rho i)*Omega:
    the real (K0/I0) and imaginary (J0/Y0) planes each pair within themselves (cross-sector
    pairings ~1e-20), the two planes are EQUAL, and B[g+1,g-1]/pi = rho i exactly -- so the
    integer form realized is Omega (rho i is the thimble normalization; the physical period-cut
    normalization gives the real -pi,0,2pi of the companion test)."""
    import mpmath as mp
    mp.mp.dps = 25
    for rho in (mp.mpf(1) / 5, mp.mpf(1) / 2):
        B = _if_thimble_form(rho, Nf=6000)
        cross = max(abs(B[i, j]) for (i, j) in [(0, 2), (0, 3), (1, 2), (1, 3)])
        assert cross < mp.mpf(10) ** -18, (rho, cross)          # block-diagonal
        assert abs(B[0, 1] - B[2, 3]) < mp.mpf(10) ** -18       # two planes equal
        assert abs(B[0, 1] / (rho * 1j) - 1) < mp.mpf(10) ** -18  # = rho i (=> integer form Omega)


def test_intersection_form_pi_and_planes_are_symbolic():
    """Fully SYMBOLIC origin of the pi and of plane-equality (upgrades the measured inputs of
    the intersection-form theorem).  From the leading branch-point asymptotics
    s_c(D) ~ e^{-x_c D} sqrt(pi/(Q'(x_c) D)) (the D-constant concomitant equals its D->oo
    limit; subleading O(1/D) corrections die), the within-sector value is
    B[s_x,s_-x] = 2 A B x (-2 rho x^2+2rho-1), A B = pi/sqrt(Q'(x)Q'(-x)); the exact identity
    Q'(x)Q'(-x) = -4 x^2 (2 rho x^2-2 rho+1)^2 cancels the algebraic factor, giving
    B[s_x,s_-x]^2 = -pi^2 independent of x AND rho -- the pi (= Gamma(1/2)^2) and the equality
    of the two symplectic planes, proved rather than measured (self-contained port; see
    _if_branch_point_proof above)."""
    import sympy as sp
    ident, Bconst, Bsq = _if_branch_point_proof()
    assert ident is True                                    # Q'(x)Q'(-x) = -4 x^2 (2rho x^2-2rho+1)^2
    assert sp.simplify(Bsq + sp.pi ** 2) == 0              # B[s_x,s_-x]^2 = -pi^2 (x- and rho-independent)
    x, rho, A, B = sp.symbols('x rho A B')
    assert sp.simplify(Bconst - 2 * A * B * x * (-2 * rho * x ** 2 + 2 * rho - 1)) == 0


# ---------------------------------------------------------------------------
# F12 / explicitly-correlated integrals in the Fock momentum representation
# (Paper 59 sec:f12).  Drivers: debug/fock_f12_momentum_probe.py, fock_f12_genus_probe.py.
# ---------------------------------------------------------------------------
# ---- self-contained ports of the F12 momentum probe (was debug/fock_f12_momentum_probe.py;
# inlined per the transient-dir policy so the cited sec:f12 [MEASURED] backing runs from a clean checkout) ----
def _f12_Dft(k, a):
    """FT of the normalized 1s density |(a^3/pi)^{1/2} e^{-a r}|^2 : 16 a^4/(k^2+4a^2)^2."""
    return 16 * a ** 4 / (k * k + 4 * a * a) ** 2


def _f12_kern_coulomb(k):
    import mpmath as mp
    return 4 * mp.pi / (k * k)


def _f12_kern_geminal(k, g):
    import mpmath as mp
    return 8 * mp.pi * g / (k * k + g * g) ** 2


def _f12_J_momentum(kern, a, b, R):
    """(1/2pi^2) int_0^inf k^2 f~(k) D_a(k) D_b(k) j0(kR) dk."""
    import mpmath as mp
    def integrand(k):
        j0 = mp.sin(k * R) / (k * R) if R > 0 else mp.mpf(1)
        return k * k * kern(k) * _f12_Dft(k, a) * _f12_Dft(k, b) * j0
    return mp.quad(integrand, [0, a, 2 * a, 4 * a, 8 * a, mp.inf]) / (2 * mp.pi ** 2)


def _f12_J_direct_samecenter(fpos, a, nlag=60):
    """<rho rho|f(r12)> for two 1s densities at the SAME center via perimetric Gauss-Laguerre."""
    import numpy as np
    from numpy.polynomial.laguerre import laggauss
    xa, wa = laggauss(nlag)
    uu = xa / a; vv = xa / a; ww = xa / (2 * a)
    A, B, C = np.meshgrid(range(nlag), range(nlag), range(nlag), indexing='ij')
    u = uu[A]; v = vv[B]; w = ww[C]
    r1 = (v + w) / 2; r2 = (u + w) / 2; r12 = (u + v) / 2
    dens = (a ** 3 / np.pi) ** 2
    gg = 8 * np.pi ** 2 * dens * fpos(r1, r2, r12) * r1 * r2 * r12
    weight = (wa[A] * wa[B] * wa[C]) / (a * a * 2 * a) * 0.25
    return float(np.sum(weight * gg))


def test_f12_kernel_swap_and_native_geminal():
    """A correlated two-electron integral <rho1|f(r12)|rho2> is the SAME momentum object as
    the ERI with 4pi/k^2 -> f~(k) (kernel-swap).  The Slater geminal f=e^{-g r12} has the
    EXACT regular FT 8pi g/(k^2+g^2)^2 (native; no Gaussian expansion).  Validated: Coulomb
    same-center reproduces the exact 5a/8; the geminal momentum form == a direct
    position-space perimetric quadrature to <1e-10."""
    import numpy as np
    import mpmath as mp
    mp.mp.dps = 25
    a = 1.0
    jc = _f12_J_momentum(_f12_kern_coulomb, a, a, 0.0)
    assert abs(float(jc) - 5*a/8) < 1e-8                              # Coulomb 5a/8
    for g in (0.5, 1.0, 2.0):
        jm = _f12_J_momentum(lambda k, g=g: _f12_kern_geminal(k, g), a, a, 0.0)
        jd = _f12_J_direct_samecenter(lambda r1, r2, r12, g=g: np.exp(-g*r12), a)
        assert abs(float(jm) - jd) < 1e-10, (g, float(jm), jd)       # kernel-swap == direct


def test_f12_geminal_pair_is_elementary():
    """The 2-center geminal PAIR integral (single-center densities) is genus-0/elementary:
    int_0^inf [k^2/(k^2+g^2)^2] j0(kR) dk = (pi/4g) e^{-g R} in closed form."""
    import sympy as sp
    k, R, g = sp.symbols('k R g', positive=True)
    closed = sp.simplify(sp.integrate(k*sp.sin(k*R)/(k**2 + g**2)**2, (k, 0, sp.oo))/R)
    assert sp.simplify(closed - sp.pi*sp.exp(-g*R)/(4*g)) == 0


def test_f12_multicenter_shares_eri_elliptic_genus():
    """The geminal does NOT lower the genus for multi-center integrals.  The two-scale curve
    period is the complete elliptic integral (genus-1 signature), set by the TWO two-center
    densities, not the kernel; a rational kernel over the genus-1 curve is still an elliptic
    integral.  The geminal moment is elementary on the diagonal c1=c2 and NOT off it."""
    import mpmath as mp
    mp.mp.dps = 25

    def P1(c1, c2):
        return mp.quad(lambda k: 1/mp.sqrt((c1*k*k+1)*(c2*k*k+1)), [0, 1, 4, mp.inf])
    for c1, c2 in [(2.0, 1.0), (5.0, 0.3)]:
        cmax, cmin = max(c1, c2), min(c1, c2)
        assert abs(P1(c1, c2) - (1/mp.sqrt(cmax))*mp.ellipk(1 - cmin/cmax)) < mp.mpf(10)**-14
    g = 1.0

    def Mg(c1, c2):
        return mp.quad(lambda k: k*k*8*mp.pi*g/(k*k+g*g)**2/mp.sqrt((c1*k*k+1)*(c2*k*k+1)),
                       [0, g, 1, 4, mp.inf])

    def Melem(c):
        return mp.quad(lambda k: k*k*8*mp.pi*g/(k*k+g*g)**2/(c*k*k+1), [0, g, 1, 4, mp.inf])
    assert abs(Mg(3, 3) - Melem(3)) < mp.mpf(10)**-12                 # diagonal: elementary
    assert abs(Mg(5, 0.3) - Melem(mp.sqrt(5*0.3)))/abs(Mg(5, 0.3)) > 0.05   # off-diagonal: elliptic


# ---------------------------------------------------------------------------
# Paper 61 sec:modular -- the corrected PSLQ ring (quasiperiod + Eisenstein
# L-value G) and the oscillation-free tail-analytic fibre.  These back the
# ring-correction and tail-fibre paragraphs added to sec:modular.
# ---------------------------------------------------------------------------

def test_catalan_G_native_to_X2_modular_integrals():
    """G = Catalan = beta(2) = L(2,chi_-4) is native to X(2) integrals of the
    first/second-kind elliptic integrals: int_0^1 K(k) dk = 2G, int_0^1 E(k) dk
    = G + 1/2.  This justifies including G in the closed-form ring (sec:modular)."""
    import mpmath as mp
    mp.mp.dps = 40
    G = mp.catalan
    IK = mp.quad(lambda k: mp.ellipk(k * k), [0, 1])   # ellipk takes m=k^2
    IE = mp.quad(lambda k: mp.ellipe(k * k), [0, 1])
    assert abs(IK - 2 * G) < mp.mpf(10) ** -35
    assert abs(IE - (G + mp.mpf(1) / 2)) < mp.mpf(10) ** -35


def test_legendre_quasiperiod_at_disc4():
    """E(1/2) = pi/(4 varpi) + varpi/2 with varpi = K(1/2): the second-kind
    quasiperiod is 1/varpi over Q[pi,varpi], so a period-only ring cannot
    represent a quasiperiod closed form (sec:modular ring correction)."""
    import mpmath as mp
    mp.mp.dps = 40
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(mp.pi))
    assert abs(varpi - mp.ellipk(mp.mpf('0.5'))) < mp.mpf(10) ** -35   # varpi = K(1/2)
    lhs = mp.ellipe(mp.mpf('0.5'))
    rhs = mp.pi / (4 * varpi) + varpi / 2
    assert abs(lhs - rhs) < mp.mpf(10) ** -35


def _hp_leggauss(N):
    """Gauss-Legendre nodes/weights at mpmath precision: numpy double seed + Newton
    refinement via the Legendre recurrence (so the fibre tests reach full precision)."""
    import mpmath as mp
    from numpy.polynomial.legendre import leggauss as _npgl
    x0, _ = _npgl(N)
    xs, ws = [], []
    for xd in x0:
        x = mp.mpf(float(xd))
        for _ in range(6):
            p0, p1 = mp.mpf(1), x
            for k in range(2, N + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dP = N * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dP; x -= dx
            if abs(dx) < mp.mpf(10) ** -(mp.mp.dps + 4):
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, N + 1):
            p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
        dP = N * (x * p1 - p0) / (x * x - 1)
        xs.append(x); ws.append(2 / ((1 - x * x) * dP * dP))
    return xs, ws


def _P(x, k):
    import mpmath as mp
    c = x * (1 - x)
    D = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-D) * (1 / D ** 3 + 3 / D ** 4 + 3 / D ** 5)


def _fibre_direct(s, t, Nk):
    """Direct decay-map fibre J(s,t) = int_0^inf j0(k(s+t)) P(s,k)P(t,k) dk."""
    import mpmath as mp
    cs, ct = s * (1 - s), t * (1 - t)
    L = 1 / (mp.sqrt(cs) + mp.sqrt(ct))
    b = s + t
    xs, ws = _hp_leggauss(Nk)
    tot = mp.mpf(0)
    for x, w in zip(xs, ws):
        u = (x + 1) / 2
        k = L * u / (1 - u); dk = L / (1 - u) ** 2
        j0 = mp.sin(k * b) / (k * b)
        tot += (w / 2) * j0 * _P(s, k) * _P(t, k) * dk
    return tot


def _fibre_tail_analytic(s, t, K, Nq, M):
    """Oscillation-free fibre: bounded int_0^K + analytic incomplete-gamma tail.
    g(k)=P1 P2 e^{Ak} k^6 is cancellation-stable; its 1/k-series (Vandermonde fit)
    weights the tail; one gammainc(-6,zK) + downward recurrence."""
    import mpmath as mp
    A = mp.sqrt(s * (1 - s)) + mp.sqrt(t * (1 - t)); b = s + t; z = A - 1j * b
    # bounded part
    xs, ws = _hp_leggauss(Nq); bnd = mp.mpf(0)
    for x, w in zip(xs, ws):
        k = K * (x + 1) / 2; wk = K * w / 2
        j0 = mp.sin(k * b) / (k * b)
        bnd += wk * j0 * _P(s, k) * _P(t, k)
    # tail coefficients from a Vandermonde fit of g(k)
    def g(k):
        return _P(s, k) * _P(t, k) * mp.e ** (A * k) * k ** 6
    k0 = max(mp.mpf(60), 35 / mp.sqrt(min(s * (1 - s), t * (1 - t))))
    ks = [k0 * (i + 1) for i in range(M + 1)]
    V = mp.matrix(M + 1, M + 1)
    for i, k in enumerate(ks):
        for m in range(M + 1):
            V[i, m] = 1 / k ** m
    d = mp.lu_solve(V, mp.matrix([g(k) for k in ks]))
    # tail = (1/b) Im[ sum_m d_m z^{6+m} Gamma(-6-m, zK) ]
    x = z * K; ex = mp.e ** (-x); Gs = []; Gcur = mp.gammainc(-6, x)
    for m in range(M + 1):
        Gs.append(Gcur); a = -6 - m
        Gcur = (Gcur - x ** (a - 1) * ex) / (a - 1)
    tot = mp.mpc(0)
    for m in range(M + 1):
        tot += d[m] * z ** (6 + m) * Gs[m]
    return bnd + (tot / b).imag


def test_tail_analytic_fibre_matches_direct():
    """The oscillation-free tail-analytic fibre reproduces the direct decay-map
    fibre to high precision, at an interior point and an edge point (sec:modular
    'the fibre precision wall is removable')."""
    import mpmath as mp
    mp.mp.dps = 40
    for s, t in [(mp.mpf('0.3'), mp.mpf('0.17')), (mp.mpf('0.05'), mp.mpf('0.28'))]:
        ref = _fibre_direct(s, t, 500)
        val = _fibre_tail_analytic(s, t, mp.mpf(110), 500, 6)
        assert abs(val - ref) < mp.mpf(10) ** -30


def test_corrected_ring_weight2_negative():
    """Guarded integer-relation check: the ~19-digit collinear value W=T2*pi/8 has
    NO low-height relation in the weight-<=2 corrected ring {pi, varpi, 1/varpi, G}
    (signed powers), and a same-magnitude decoy behaves the same -- the decisive
    weight-<=2 negative that pins the closed form to weight-three-requiring-G."""
    import mpmath as mp
    mp.mp.dps = 19
    varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(mp.pi))
    G = mp.catalan; pi = mp.pi
    # weight<=2 monomials pi^a varpi^b G^d, a>=0, b in Z, d in {0,1}, a+|b|+2d in 1..2
    basis = [mp.mpf(1)]
    for a in range(3):
        for bb in range(-2, 3):
            for dd in (0, 1):
                w = a + abs(bb) + 2 * dd
                if 1 <= w <= 2:
                    basis.append(pi ** a * varpi ** bb * G ** dd)
    V = mp.mpf('0.3953557659017139641')        # canonical collinear value (DoD C8 #7, ~19 dig)
    W = V * pi / 8
    decoy = W * (1 + mp.mpf(10) ** -12) + mp.euler / mp.mpf(10) ** 5
    rel_r = mp.pslq([W] + basis, maxcoeff=10 ** 6, maxsteps=10 ** 5)
    rel_d = mp.pslq([decoy] + basis, maxcoeff=10 ** 6, maxsteps=10 ** 5)
    def height(r):
        return None if r is None else max(abs(x) for x in r)
    hr, hd = height(rel_r), height(rel_d)
    # negative = either no relation, or the real relation is no smaller than the decoy's
    assert hr is None or (hd is not None and hr >= hd / 10)


def test_paper59_diagonal_A():
    """The T2 cusp series' leading Fourier-Whittaker coefficient A = J(1/4,1/4,b=1)/4 is the
    genus-0 diagonal value.  Three pins (sec:modular + sprint_T2_zhou_bd_adaptation_memo.md):
      (1) A = 0.0790540321168... (direct 1D quadrature);
      (2) A is a SINGLE-SCALE Bessel moment -- the load-bearing inner identity
          int cos(kb) e^{-t sqrt(k^2/4+1)} dk = 2 t K1(sqrt(t^2+4 b^2))/sqrt(t^2+4 b^2)
          (= -d/dt of eq:K0), which makes A a sigma-integral of K1 (NOT an elementary number);
      (3) guarded, decoy-controlled: A has NO low-height closure in the weight-1 two-center ring
          {pi, e^-a} at the natural arguments a in {2, 2 sqrt2} -- so even the leading cusp datum is
          an irreducible Bessel moment, not a divisor-sum/elementary quantity (one twist-type beyond
          the Broadhurst-Dorigoni character resummation)."""
    import mpmath as mp
    mp.mp.dps = 34

    def P(k):
        d = mp.sqrt(k * k / 4 + 1)
        return (mp.mpf(1) / 4) * mp.e ** (-d) * (d ** -3 + 3 * d ** -4 + 3 * d ** -5)

    # (1) A value
    def integ(k):
        j0 = mp.sin(k) / k if k > mp.mpf('1e-30') else mp.mpf(1)
        return j0 * P(k) ** 2
    A = mp.quad(integ, [0, 1, 2, 4, 8, 16, mp.inf]) / 4
    assert abs(A - mp.mpf('0.07905403211681687671461246470412103')) < mp.mpf(10) ** -30

    # (2) inner K1-moment identity (the exact analytic step of the Bessel-moment reduction)
    for t, b in [(mp.mpf(2), mp.mpf('0.5')), (mp.mpf(3), mp.mpf('0.7'))]:
        lhs = mp.quadosc(lambda k: mp.cos(k * b) * mp.e ** (-t * mp.sqrt(k * k / 4 + 1)),
                         [0, mp.inf], period=2 * mp.pi / b)
        r = mp.sqrt(t * t + 4 * b * b)
        assert abs(lhs - 2 * t * mp.besselk(1, r) / r) < mp.mpf(10) ** -25

    # (3) guarded PSLQ negative in the weight-1 ring at args {2, 2 sqrt2}
    a2 = 2 * mp.sqrt(2)
    basis = [mp.pi, mp.e ** -2, mp.e ** (-a2)]
    decoy = A * mp.mpf('1.0000000000031415926535') + mp.mpf('1e-9')
    rel = mp.pslq([A] + basis, tol=mp.mpf(10) ** -26, maxcoeff=10 ** 6, maxsteps=10 ** 6)
    reld = mp.pslq([decoy] + basis, tol=mp.mpf(10) ** -26, maxcoeff=10 ** 6, maxsteps=10 ** 6)
    ht = None if rel is None else max(abs(c) for c in rel)
    htd = None if reld is None else max(abs(c) for c in reld)
    # negative: no relation, or A's relation is no smaller than the decoy's (basis-coverage artifact)
    assert ht is None or (htd is not None and ht >= htd / 10)
