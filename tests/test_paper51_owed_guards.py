"""Owed guards for Paper 51 sec. 10 ("Owed tests"), 2026-11-07 narrowing.

Paper: papers/group5_qed_gauge/paper_51_gravity_arc.tex.

Every test names the wrong answer it rejects.  Expected values come from
routes that do not share the code under test: exact continuum Bessel zeros
(scipy root finding on J_nu, independent of the finite-difference
substrate), the Sommerfeld--Cheeger / Kabat / Fursaev--Miele closed form,
sympy Bernoulli polynomials, a Weyl closed form, or the printed number.

Loci backed (paper labels):
  * sec:cone "Continuum": tip(t) = c'_{1/2}(1) = +1/6 to 2e-8 (exact Bessel
    zeros at t/R^2 = 0.01, 0.02, 0.05), the finite-disk tip 0.166520
    (-0.088%) at t/R^2 = 0.1, and the total-trace coefficient
    -(1/12)(1/alpha - alpha), both branches;
  * sec:discrete: continuum a_0 = 1.888 at R = 10, t = 0.1, substrate 5.5% high;
  * sec:cone "Excess angle" (B4): recovery 0.675 ... 0.919, 1 - rec ~ sqrt(a);
  * sec:L6 "(B1, rate.)" and "(B2, target.)" at t = 1, and the B2 t = 10
    panel (a = 0.1 ... 0.00625, Hellmann--Feynman tip errors -5.32 ...
    -0.42% against 1/6, converging from below to 0.166520) with the +0.5%
    bias and spurious sign change of a step-0.1 central difference in alpha;
  * rem:two_term_uniqueness: the S^5 / S^7 identity and m+1 power terms;
  * subsec:entropy_route: the rank-4 trace and the total-trace cone
    coefficient as tested inputs to G_action = G_entropy.

The B2 t = 10 panel is guarded (test_paper51_b2_t10_panel_hellmann_feynman,
test_paper51_b2_t10_central_difference_artifact); its tip errors are
relative to 1/6, as printed.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import sympy as sp
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from scipy.special import jv, jvp

SIXTH = 1.0 / 6.0


# ---------------------------------------------------------------------------
# Continuum: exact Bessel zeros of J_nu (Dirichlet disk/cone of radius 1)
# ---------------------------------------------------------------------------

def _bessel_zeros(nu: float, zmax: float) -> np.ndarray:
    """All zeros of J_nu in (nu, zmax].  For nu >= 1/2 the first zero
    exceeds nu and consecutive zeros are more than pi apart, so a scan at
    step 0.5 brackets each zero exactly once."""
    xs = np.arange(nu, zmax + 1.0, 0.5)
    if len(xs) < 2:
        return np.array([])
    f = jv(nu, xs)
    idx = np.nonzero(np.sign(f[:-1]) * np.sign(f[1:]) < 0)[0]
    return np.array([brentq(lambda x: jv(nu, x), xs[i], xs[i + 1],
                            xtol=1e-15, rtol=1e-15, maxiter=200) for i in idx])


def _dzero_dnu(nu: float, z: np.ndarray) -> np.ndarray:
    """d j_{nu,k} / d nu = -(dJ/dnu)/(dJ/dx) at the zero (implicit function
    theorem); dJ/dnu by a Richardson-extrapolated central difference."""
    def d(h):
        return (jv(nu + h, z) - jv(nu - h, z)) / (2 * h)
    dJ = (4 * d(5e-4) - d(1e-3)) / 3
    return -dJ / jvp(nu, z)


def _cone_K(alpha: float, s: float, with_dalpha: bool = False):
    """Continuum spinor cone heat trace K(alpha, t) at s = t/R^2:
    K = 4 sum_j sum_k exp(-s j_{m,k}^2), m = (j + 1/2)/alpha
    (two components x the two signs of m).  Optionally dK/dalpha."""
    zmax = math.sqrt(40.0 / s)
    K = dK = 0.0
    j = 0
    while True:
        m = (j + 0.5) / alpha
        if m > zmax:
            break
        z = _bessel_zeros(m, zmax)
        if len(z) == 0:
            break
        e = np.exp(-s * z ** 2)
        K += 4 * float(e.sum())
        if with_dalpha:
            # dm/dalpha = -m/alpha
            hp = float(np.sum(-2 * s * z * _dzero_dnu(m, z) * e))
            dK += 4 * hp * (-m / alpha)
        j += 1
    return (K, dK) if with_dalpha else K


def test_paper51_bessel_root_finder_control() -> None:
    """Control for the routes below: the zeros of J_{1/2} are k pi exactly.
    Rejects: a root finder that skips or duplicates zeros."""
    z = _bessel_zeros(0.5, 40.0)
    assert np.allclose(z / np.pi, np.arange(1, len(z) + 1), atol=1e-13)
    assert len(z) == 12


def test_paper51_continuum_tip_one_sixth() -> None:
    """sec:cone (eq:tip_def): the continuum tip dK/dalpha|_1 - K(1,t) equals
    c'_{1/2}(1) = +1/6, from exact Bessel zeros at t/R^2 = 0.01, 0.02, 0.05
    to the paper's 2e-8 (measured 1e-11, 4e-12, 9.9e-9).

    Rejects: +1/12 (a per-component reading of the coefficient), -1/6
    (the scalar sign, c = +(1/12)(1/alpha - alpha)), and a finite-disk
    correction at t/R^2 <= 0.05 above the printed 2e-8."""
    for s in (0.01, 0.02, 0.05):
        K, dK = _cone_K(1.0, s, with_dalpha=True)
        tip = dK - K
        assert abs(tip - SIXTH) < 2e-8, (s, tip)


TIP_CONT_S01 = 0.166520        # printed finite-disk tip at t/R^2 = 0.1


def test_paper51_continuum_tip_finite_disk_s01() -> None:
    """sec:cone "Continuum" and sec:L6 (B2, t = 10 panel): at t/R^2 = 0.1
    the finite-disk continuum tip is 0.166520, -0.088% from 1/6 (exact
    Bessel zeros; measured 0.16651956).

    Rejects: 1/6 as the t/R^2 = 0.1 target (off by 1.5e-4)."""
    K, dK = _cone_K(1.0, 0.1, with_dalpha=True)
    tip = dK - K
    assert abs(tip - TIP_CONT_S01) < 5e-7, tip
    assert abs(100 * (tip / SIXTH - 1) - (-0.088)) < 5e-4, tip
    assert abs(tip - SIXTH) > 1e-4, tip


def test_paper51_total_trace_cone_coefficient() -> None:
    """sec:cone eq:spinor_sc_tip, both branches: K(alpha,t) - alpha K(1,t)
    = -(1/12)(1/alpha - alpha) for the two-component spinor, total trace,
    at alpha in {1/2, 2/3, 3/2, 2, 3} from exact Bessel zeros (agreement to
    ~1e-13); alpha > 1 is the plain continuation (sec:cone "Excess angle").

    Rejects: +1/12 (scalar sign) and -1/24 (per component)."""
    s = 0.01
    K1 = _cone_K(1.0, s)
    for al in (0.5, 2.0 / 3.0, 1.5, 2.0, 3.0):
        d = _cone_K(al, s) - al * K1
        c = -(1.0 / 12.0) * (1.0 / al - al)
        assert abs(d - c) < 1e-9, (al, d, c)


def test_paper51_continuum_a0_and_substrate_excess() -> None:
    """sec:discrete: the estimator a_0 = K(t) 4 pi t / (pi R^2) on the
    continuum Dirichlet spinor disk is 1.888 at R = 10, t = 0.1 (exact
    Bessel zeros; Weyl closed form 2 - 2 sqrt(pi t)/R = 1.8879 agrees to the
    O(t) constant term), and the substrate value 1.992 (N_phi = 192,
    N_rho = 200, a = 0.05) is 5.5% above it.

    Rejects: the bulk value 2.0 as the continuum target, and 1.992 as a
    recovery."""
    from geovac.gravity.warped_dirac import DiscreteDiskDirac
    R, t = 10.0, 0.1
    a0 = _cone_K(1.0, t / R ** 2) * 4 * math.pi * t / (math.pi * R ** 2)
    weyl = 2.0 - 2.0 * math.sqrt(math.pi * t) / R
    assert abs(a0 - 1.888) < 5e-4, a0
    assert abs(a0 - weyl) < 1e-3, (a0, weyl)
    sub = DiscreteDiskDirac(N_rho=200, a=0.05, N_phi=192).heat_trace(t) \
        * 4 * math.pi * t / (math.pi * R ** 2)
    assert abs(sub - 1.992) < 1e-3, sub
    assert abs(100 * (sub / a0 - 1) - 5.5) < 0.05, sub / a0


# ---------------------------------------------------------------------------
# Discrete substrate: lean tridiagonal radial solver (rho_i = i a)
# ---------------------------------------------------------------------------

def _radial(m: float, R: float, N: int, t: float, corrected: bool = False,
            vectors: bool = False):
    """Low eigenpairs (lambda < 45/t) of L(m) = -d^2 + (m^2 - 1/4)/rho^2 on
    rho_i = i a, Dirichlet beyond the last site.  Production grid a = R/N
    (wall at R + a); corrected grid a = R/(N+1) (wall at R)."""
    a = R / (N + 1) if corrected else R / N
    rho = np.arange(1, N + 1) * a
    d = 2 / a ** 2 + (m * m - 0.25) / rho ** 2
    e = -np.ones(N - 1) / a ** 2
    out = eigh_tridiagonal(d, e, eigvals_only=not vectors, select="v",
                           select_range=(-1e9, 45.0 / t), lapack_driver="stemr")
    return (out, rho) if vectors else out


def _K_disc(alpha: float, R: float, N: int, t: float) -> float:
    """Discrete spinor wedge heat trace with exact azimuthal momenta."""
    tot = 0.0
    j = 0
    while True:
        v = float(np.exp(-t * _radial((j + 0.5) / alpha, R, N, t)).sum())
        tot += 4 * v
        if v < 1e-18:
            return tot
        j += 1


def _tip_disc(R: float, N: int, t: float) -> float:
    """Discrete tip dK/dalpha|_1 - K(1,t); the m-derivative by
    Hellmann--Feynman, d lambda/dm = 2 m <rho^-2>."""
    K = dK = 0.0
    j = 0
    while True:
        m = j + 0.5
        (w, V), rho = _radial(m, R, N, t, vectors=True)
        e = np.exp(-t * w)
        if e.sum() < 1e-18:
            return dK - K
        dl = 2 * m * np.sum(V ** 2 / rho[:, None] ** 2, axis=0)
        K += 4 * float(e.sum())
        dK += -4 * m * float(np.sum(-t * dl * e))
        j += 1


def _recovery(alpha: float, R: float, N: int, t: float) -> float:
    return ((_K_disc(alpha, R, N, t) - alpha * _K_disc(1.0, R, N, t))
            / (1 / alpha - alpha)) / (-1.0 / 12.0)


def test_paper51_b4_excess_angle_recovery_sqrt_a() -> None:
    """sec:cone "Excess angle" (B4): at alpha = 2, t = 1, R = 10 the
    recovery rises as N_rho = 200 -> 3200 (0.675, 0.769, 0.837, 0.885,
    0.919), with successive error ratios 1.4085, 1.4143, 1.4164 (about
    sqrt 2, i.e. 1 - recovery ~ sqrt a).

    Rejects: the fixed-a plateau (~0.68) read as the continuum excess-angle
    coefficient, and an O(a) error law (ratio 2)."""
    R, t = 10.0, 1.0
    Ns = (200, 400, 800, 1600, 3200)
    rec = [_recovery(2.0, R, N, t) for N in Ns]
    for got, want in zip(rec, (0.675, 0.769, 0.837, 0.885, 0.919)):
        assert abs(got - want) < 6e-4, rec
    ratios = [(1 - rec[i]) / (1 - rec[i + 1]) for i in range(len(rec) - 1)]
    for got, want in zip(ratios[:3], (1.4085, 1.4143, 1.4164)):
        assert abs(got - want) < 5e-4, ratios
    assert all(abs(r - math.sqrt(2)) < 0.01 for r in ratios), ratios


@pytest.mark.slow
def test_paper51_b4_one_cutoff_past() -> None:
    """B4 one cutoff past the paper's data (rule 4): N_rho = 6400 keeps the
    error ratio at sqrt 2 (measured 1.4165, after 1.4168 at 3200/1600;
    the ratios settle about 0.002 above sqrt 2) and the recovery rising.

    Rejects: a sqrt(a) law that was a fit artifact of the five printed
    points (ratio drifting toward 1 or 2 beyond them)."""
    R, t = 10.0, 1.0
    r16, r32, r64 = (_recovery(2.0, R, N, t) for N in (1600, 3200, 6400))
    assert r16 < r32 < r64 < 1.0
    assert abs((1 - r32) / (1 - r64) - math.sqrt(2)) < 3e-3


def test_paper51_b4_lean_solver_matches_production() -> None:
    """The lean solver above is the production substrate: at N_rho = 200,
    a = 0.05, alpha = 2, t = 1, the recovery from
    geovac.gravity.warped_dirac.DiscreteWedgeDiracSpectral /
    DiscreteDiskDiracSpectral agrees with it to 1e-8.

    Rejects: a production change of the radial operator (e.g. the
    centrifugal (m^2 - 1/4) term) that the lean B4 route would not see."""
    from geovac.gravity.warped_dirac import (
        DiscreteDiskDiracSpectral, DiscreteWedgeDiracSpectral)
    t = 1.0
    Kw = DiscreteWedgeDiracSpectral(N_rho=200, a=0.05, N_phi=280, alpha=2.0).heat_trace(t)
    Kd = DiscreteDiskDiracSpectral(N_rho=200, a=0.05, N_phi=140).heat_trace(t)
    rec_prod = ((Kw - 2 * Kd) / (0.5 - 2.0)) / (-1.0 / 12.0)
    assert abs(rec_prod - _recovery(2.0, 10.0, 200, t)) < 1e-8, rec_prod


def _bessel_lambda0(m: float, R: float) -> float:
    return (_bessel_zeros(m, m + 20.0)[0] / R) ** 2


def test_paper51_b1_eigenvalue_order_and_wall_placement() -> None:
    """sec:L6 (B1, rate): lambda_0(m; a) against the exact (j_{m,1}/R)^2
    converges at order q = 1.00 at m = 1/2, 3/2, 5/2 on the production grid
    (error ratio -> 2 per halving, N_rho = 100 ... 3200, one past the
    paper's 1600), including m = 1/2 where the centrifugal term vanishes;
    on the grid a = R/(N_rho + 1) the ratio is 4 (O(a^2)).

    Rejects: O(a^2) on the production grid (the pre-B1 belief) and the
    centrifugal term as the cause of the O(a) rate."""
    R = 10.0
    Ns = (100, 200, 400, 800, 1600, 3200)
    for m in (0.5, 1.5, 2.5):
        ex = _bessel_lambda0(m, R)
        err = [abs(_radial(m, R, N, 1.0)[0] - ex) for N in Ns]
        q = float(np.polyfit(np.log([R / N for N in Ns]), np.log(err), 1)[0])
        assert abs(q - 1.00) < 0.01, (m, q)
        errc = [abs(_radial(m, R, N, 1.0, corrected=True)[0] - ex) for N in Ns[:4]]
        rc = [errc[i] / errc[i + 1] for i in range(3)]
        assert all(abs(r - 4.0) < 0.05 for r in rc), (m, rc)


def test_paper51_b1_b2_tip_order_and_extrapolation() -> None:
    """sec:L6 (B1) tip order p = 0.97 over N_rho = 50 ... 1600 (R = 10,
    t = 1), and (B2) linear-in-a extrapolation of the tip to 1/6: the
    two-point value (N_rho = 800, 1600) 0.166654, -0.008% (measured
    0.16665377), and the least-squares value (four finest) 0.166582,
    -0.05%.  At t/R^2 = 0.01 the continuum tip is 1/6 to 2e-8 (test above),
    so 1/6 is the right target here.  The m-derivative is Hellmann--Feynman
    (exact for the discrete operator).

    Rejects: the retired two-point 0.166656 (-0.006%, from a finite
    difference with step 1e-4 carrying ~1e-6 eigenvalue noise), an
    extrapolation that assumes the wrong order (order-2 Richardson on this
    O(a) sequence lands farther from 1/6), and p = 2."""
    R, t = 10.0, 1.0
    Ns = (50, 100, 200, 400, 800, 1600)
    a = np.array([R / N for N in Ns])
    tips = np.array([_tip_disc(R, N, t) for N in Ns])
    p = float(np.polyfit(np.log(a), np.log(np.abs(SIXTH - tips)), 1)[0])
    assert abs(p - 0.97) < 0.005, p
    two_pt = 2 * tips[-1] - tips[-2]
    assert abs(two_pt - B2_TWO_POINT) < 5e-7, two_pt
    rel2 = 100 * (two_pt - SIXTH) / SIXTH
    assert abs(rel2 - (-0.008)) < 5e-4, rel2
    A = np.vstack([a[-4:], np.ones(4)]).T
    lsq = float(np.linalg.lstsq(A, tips[-4:], rcond=None)[0][1])
    assert abs(lsq - 0.166582) < 1e-6, lsq
    assert abs(100 * (lsq - SIXTH) / SIXTH - (-0.05)) < 0.005, lsq
    rich2 = (4 * tips[-1] - tips[-2]) / 3   # wrong-order Richardson
    assert abs(rich2 - SIXTH) > 10 * abs(two_pt - SIXTH)


B2_TWO_POINT = 0.166654        # printed (B2) two-point extrapolation, t = 1

# (B2) t = 10 panel, R = 10 (t/R^2 = 0.1): a and the printed tip errors in %
# against 1/6, derivative in alpha exact (Hellmann--Feynman).
T10_A = (0.1, 0.05, 0.025, 0.0125, 0.00625)
T10_ERR_PCT = (-5.32, -2.73, -1.41, -0.75, -0.42)
T10_CD_STEP = 0.1              # the retired central-difference step in alpha


def _t10_panel_tips(deriv: str = "hf") -> np.ndarray:
    """The t = 10 panel's tips.  deriv = "hf": exact alpha-derivative
    (_tip_disc, Hellmann--Feynman); "cd": central difference in alpha with
    step T10_CD_STEP on the exact-momentum wedge trace (the retired route)."""
    R, t = 10.0, 10.0
    out = []
    for a in T10_A:
        N = round(R / a)
        if deriv == "hf":
            out.append(_tip_disc(R, N, t))
        else:
            h = T10_CD_STEP
            out.append((_K_disc(1 + h, R, N, t) - _K_disc(1 - h, R, N, t)) / (2 * h)
                       - _K_disc(1.0, R, N, t))
    return np.array(out)


def test_paper51_b2_t10_panel_hellmann_feynman() -> None:
    """sec:L6 (B2, target), t = 10 panel: with the alpha-derivative exact
    (Hellmann--Feynman), the tip errors against 1/6 at a = 0.1 ... 0.00625
    are -5.32, -2.73, -1.41, -0.75, -0.42% (measured -5.320, -2.728,
    -1.414, -0.753, -0.421), negative and shrinking at every a (from
    below), and the linear-in-a extrapolation of the two finest lands on
    the finite-disk continuum tip 0.166520 (within 1e-5), not on 1/6.

    Rejects: the retired sign change (+0.045% at a = 0.00625, -0.30% at
    0.0125, from a step-0.1 central difference), any finite-step
    alpha-derivative in this panel (biased +0.43 ... +0.47%, far outside
    the 0.005-point tolerance), and 1/6 as the panel's limit."""
    tips = _t10_panel_tips()
    err = 100 * (tips / SIXTH - 1)
    for got, want in zip(err, T10_ERR_PCT):
        assert abs(got - want) < 0.005, err
    assert np.all(err < 0) and np.all(np.diff(err) > 0), err
    ex = 2 * tips[-1] - tips[-2]
    assert abs(ex - TIP_CONT_S01) < 1e-5, ex
    assert abs(ex - SIXTH) > 1e-4, ex
    # approach to the finite-disk value is O(a): error ratio ~2 per halving
    e = TIP_CONT_S01 - tips
    assert all(1.9 < e[i] / e[i + 1] < 2.1 for i in range(len(e) - 1)), e


def test_paper51_b2_t10_central_difference_artifact() -> None:
    """sec:L6 (B2), t = 10 panel: "a central difference in alpha with step
    0.1 biases the tip by +0.5% and would show a spurious sign change".
    The step-0.1 central difference sits +0.43 ... +0.47% above the exact
    (Hellmann--Feynman) tip at every a, and its errors against 1/6 at the
    two finest a are the retired -0.30% and +0.045%.

    Rejects: the central-difference values as the tip (they cross 1/6,
    the exact ones do not), and a bias that is a rounding effect (< 0.4%)."""
    hf = _t10_panel_tips("hf")
    cd = _t10_panel_tips("cd")
    bias = 100 * (cd / hf - 1)
    assert np.all((bias > 0.4) & (bias < 0.5)), bias
    assert round(float(np.mean(bias)), 1) == 0.5
    err_cd = 100 * (cd / SIXTH - 1)
    assert abs(err_cd[-2] - (-0.30)) < 0.005 and abs(err_cd[-1] - 0.045) < 0.0005, err_cd
    assert err_cd[-1] > 0 > err_cd[-2]
    assert np.all(hf < SIXTH)


# ---------------------------------------------------------------------------
# rem:two_term_uniqueness -- odd spheres
# ---------------------------------------------------------------------------

def _ch_degeneracy_poly(d: int, shift):
    """CH Dirac on S^d, d = 2m+1: |lambda_n| = n + shift (shift = d/2),
    total degeneracy g_n = 2^{m+1} C(n + 2m, 2m), as a polynomial in
    x = n + shift."""
    m = (d - 1) // 2
    x = sp.symbols("x")
    n = x - shift
    g = 2 ** (m + 1) * sp.prod([(n + i) for i in range(1, 2 * m + 1)]) / sp.factorial(2 * m)
    return sp.Poly(sp.expand(g), x), x


def _zeta_at_minus_k(d: int, k: int, shift) -> sp.Expr:
    """zeta(s) = sum_n g_n x^{-2s} at s = -k via Hurwitz:
    zeta_H(-p, a) = -B_{p+1}(a)/(p+1)."""
    P, x = _ch_degeneracy_poly(d, shift)
    tot = sp.Integer(0)
    for (power,), c in P.terms():
        p = 2 * k + power
        tot += c * (-sp.bernoulli(p + 1, shift) / (p + 1))
    return sp.nsimplify(tot)


def test_paper51_odd_sphere_zeta_vanishes() -> None:
    """rem:two_term_uniqueness: zeta_unit(-k) = 0 on S^3, S^5, S^7 for
    k = 0..7 (sympy Bernoulli polynomials), and the degeneracy is even in
    the shifted eigenvalue.

    Rejects: "S^3-specific" (the identity failing on S^5 or S^7)."""
    for d in (3, 5, 7):
        P, x = _ch_degeneracy_poly(d, sp.Rational(d, 2))
        assert sp.expand(P.as_expr().subs(x, -x) - P.as_expr()) == 0
        for k in range(8):
            assert _zeta_at_minus_k(d, k, sp.Rational(d, 2)) == 0, (d, k)


def test_paper51_odd_sphere_m_plus_1_power_terms() -> None:
    """rem:two_term_uniqueness: on S^{2m+1} the heat trace has exactly m+1
    power-law terms, sum_j c_j (1/2) Gamma(j + 1/2) t^{-j-1/2}, j = 0..m;
    the direct mode sum matches it to 1e-12 at t = 0.02, 0.05, 0.1
    (Poisson remainder ~ e^{-pi^2/t}).

    Rejects: a further power term (any nonzero t^{1/2} coefficient)."""
    for d in (5, 7):
        m = (d - 1) // 2
        P, x = _ch_degeneracy_poly(d, sp.Rational(d, 2))
        coeffs = {pw[0]: float(c) for pw, c in P.terms()}
        assert sorted(coeffs) == [2 * j for j in range(m + 1)]
        f = sp.lambdify(x, P.as_expr(), "numpy")
        n = np.arange(0, 4000)
        xs = n + d / 2.0
        for t in (0.02, 0.05, 0.1):
            direct = float(np.sum(f(xs) * np.exp(-t * xs ** 2)))
            closed = sum(coeffs[2 * j] * 0.5 * math.gamma(j + 0.5) * t ** (-j - 0.5)
                         for j in range(m + 1))
            assert abs(direct / closed - 1) < 1e-12, (d, t, direct, closed)


# ---------------------------------------------------------------------------
# subsec:entropy_route -- the rank-4 trace and the cone coefficient as inputs
# ---------------------------------------------------------------------------

def _product_trace_ratio(t: float) -> tuple:
    """Explicit operators on a truncated S^3 x S^1 spectrum.
    rank 4: D = D3 (x) sigma_1 + w (x) sigma_2 (D3, w commuting diagonals);
    rank 2: D' = D3 + w (no C^2 factor).  Returns
    (Tr e^{-tD^2} / (K3 K1), Tr e^{-tD'^2} / (K3 K1))."""
    lam = np.array([s * (n + 1.5) for n in range(6) for s in (1, -1)
                    for _ in range((n + 1) * (n + 2))])       # g_n/2 per sign
    w = 2 * np.pi * (np.arange(-6, 6) + 0.5) / 3.0             # anti-periodic, beta = 3
    L, W = np.meshgrid(lam, w, indexing="ij")
    s1 = np.array([[0, 1], [1, 0]], dtype=complex)
    s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
    # block per (lambda, omega): commuting scalars -> 2x2 blocks
    ev4 = []
    for l_, w_ in zip(L.ravel(), W.ravel()):
        B = l_ * s1 + w_ * s2
        ev4.extend(np.linalg.eigvalsh(B @ B))
    ev4 = np.array(ev4)
    ev2 = (L + W).ravel() ** 2
    K3K1 = float(np.sum(np.exp(-t * lam ** 2)) * np.sum(np.exp(-t * w ** 2)))
    return (float(np.sum(np.exp(-t * ev4))) / K3K1,
            float(np.sum(np.exp(-t * ev2))) / K3K1)


def test_paper51_rank4_trace_is_derived() -> None:
    """sec:g2: on the rank-4 bundle D^2 = D3^2 + w^2 and the heat trace is
    exactly 2 K3 K1; the rank-2 sum D3 + w has a cross term and does not
    factorize.

    Rejects: rank 2 (factor 1), the trace that gave the withdrawn factor of
    2 between the G routes."""
    for t in (0.05, 0.2):
        r4, r2 = _product_trace_ratio(t)
        assert abs(r4 - 2.0) < 1e-12, r4
        assert abs(r2 - 1.0) > 1e-2, r2


def test_paper51_G_routes_agree_from_tested_inputs() -> None:
    """subsec:entropy_route: G_action = G_entropy = 3 pi/Lambda^2 with both
    inputs measured, not assumed: the bundle factor from the explicit
    product operator (test above) and the cone tip c'(1) from exact Bessel
    zeros (test_paper51_continuum_tip_one_sixth).  Entropy:
    S = c'(1) * [leading r_h^2/t coefficient of the S^2 Dirac trace] *
    Lambda^2 (Gaussian phi(1) = 1); A/(4G) = S.

    Rejects: rank 2 or a per-component cone coefficient (each a factor 2 in
    G_action/G_entropy), and the scalar sign (G_entropy < 0)."""
    Lam = 8.0
    factor = round(_product_trace_ratio(0.1)[0], 9)
    K, dK = _cone_K(1.0, 0.01, with_dalpha=True)
    cprime = dK - K

    def density(R):
        n = np.arange(0, 4000)
        l2 = ((n + 1.5) / R) ** 2
        return factor * float(np.sum(2 * (n + 1) * (n + 2) * np.exp(-l2 / Lam ** 2))) \
            * Lam / (2 * math.sqrt(math.pi))

    _, B = np.linalg.solve(np.array([[1.0, 1.0], [8.0, 2.0]]),
                           np.array([density(1.0), density(2.0)]))
    G_action = -3 * math.pi / (4 * B)
    r_h = 2.0
    ts = np.array([1e-3, 2e-3])
    k = np.arange(1, 4000)
    Ks2 = [float(np.sum(4 * k * np.exp(-k ** 2 * tt / r_h ** 2))) for tt in ts]
    c_rh2, _ = np.linalg.solve(np.array([[1 / ts[0], 1.0], [1 / ts[1], 1.0]]),
                               np.array(Ks2))
    S_BH = cprime * c_rh2 * Lam ** 2
    G_entropy = 4 * math.pi * r_h ** 2 / (4 * S_BH)
    assert G_entropy > 0
    assert abs(G_action / G_entropy - 1) < 1e-6, (G_action, G_entropy)
    assert abs(G_action * Lam ** 2 - 3 * math.pi) < 1e-5
