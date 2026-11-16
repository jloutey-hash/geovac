"""Guard for record AX / AX2b (docs/predictions/2026-11-15_aa_bc_d_mneq0.md; amendment and Results
docs/predictions/2026-11-16_aa_bc_d_mneq0_2b.md, scored GO): the three-centre integrals

    DX = (1s_A 1s_A | 1s_B 3d_xz,C),      3d_xz    = N_xz x_C z_C e^{-zC rC},             m = 1, in the ABC plane
    DD = (1s_A 1s_A | 1s_B 3d_x2-y2,C),   3d_x2-y2 = N_x2y2 (x_C^2 - y_C^2) e^{-zC rC},   m = 2
    N_xz = sqrt(15/(4 pi)) R_n,  N_x2y2 = sqrt(15/(16 pi)) R_n,  R_n = sqrt((2 zC)^7 / 720),

over Slater functions, B at the origin, C at (0, 0, D1), A at D2 (sin g, 0, cos g), x = cos g:

    DX = sum_{L >= 1} G_L P_L^1(x),   DD = sum_{L >= 2} K_L P_L^2(x),   P_L^m = (1 - x^2)^{m/2} d^m P_L (no
    Condon-Shortley phase), against AB's kernel NA - Y + zA dY (the 1s potential of A).

Paper 59 sec:scope (inventory row 8, part 2b).  Written 2026-11-16 as a separate step after the computation
(CLAUDE.md sec. 9).  The record's engine paths (debug/aapc_exact.py, `dens_dxz`, `dens_dx2y2`) are ported here; the
base layer (Field, Bessel tables, exponential integrals, AB's 1s density, the factored kernels, the three-region
integrator) is imported from tests/test_paper59_aa_bc_p.py, which imports tests/test_paper59_aa_bc.py.  Nothing is
imported from debug/.

Claims guarded, each with the wrong answer it rejects (details in each test's docstring):

(0) the normalizations N_xz, N_x2y2 (unit norm by quadrature) and the geometry: the 3d_yz and 3d_xy integrals vanish
    because their densities are odd under y -> -y while V_A, with A in the xz-plane, is even (checked on the azimuthal
    integral at sample (r, theta), where the xz and x^2-y^2 integrands are not zero).
(1) per-channel values: G_L, K_L of the exact split equal a quadrature that projects the m = 1 or m = 2 density
    multipole directly ((r t - D1)(1 - t^2) P_L' and (1 - t^2)^2 P_L'', no T recurrence, no raising identity),
    against the multipole of the exact 1s potential V_A.  Both sides are given enough digits (the VOIDs of AX and
    AX2 were the check's precision: see the tolerance notes).
(2) the summed series at the record's eight points (both integrals, both cases, two angles each) against the R2
    referee (n = 34) within the record's tolerance 1e-10, and against the record's 2000-digit P1 values (L <= 240)
    within a stated truncation bound.  This is what tests the z_C factor of DX (the decay cannot), the
    normalizations and the addition-theorem factor.
(3) the constant set: mod 2^61 - 1, G_L for 1 <= L <= 20 and K_L for 2 <= L <= 20 carry exactly AB's 21 constants
    at the same case, no Euler gamma; G_0, K_0, K_1 are empty.
(4) decay: beta_eff(40) of G_L in P_L^1 and of K_L in P_L^2 between the midpoints of the in-test toy exponents around
    (x* - x)^{7/2} in the same basis (record P3: L^-5 and L^-6 at L = 220); the 9/2 markers (a cancelled 7/2
    coefficient) and the 5/2 branch fall outside.  One sign for 20 <= L <= 41.

Fire-test plants and outcomes: the comment block at the end of this file.
"""
from __future__ import annotations

import importlib.util
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import pytest

_spec = importlib.util.spec_from_file_location('_ap_guard_ax', Path(__file__).with_name('test_paper59_aa_bc_p.py'))
AP = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(AP)
AB = AP.AB

Field, params, _put = AB.Field, AB.params, AB._put
PRIME = (1 << 61) - 1


# ---------------------------------------------------------------------------------------------------------------
# Record-AX engine (port of debug/aapc_exact.py, paths dens_dxz and dens_dx2y2)
# ---------------------------------------------------------------------------------------------------------------


def _dacc(out, dd, w, shift):
    for (e0, b), poly in dd.items():
        for n, v in poly.items():
            _put(out, e0, b, n + shift, w * v)


def _ecache(inner, D1, zC, it, kt, K):
    cache = {}

    def E(k):
        if k not in cache:
            cache[k] = AB.density(k, inner, D1, zC, it, kt, K)
        return cache[k]
    return E


def _tmul(E, L, K):
    """(T E)_L = [L E_{L-1} + (L+1) E_{L+1}] / (2L+1): multiply by t = cos(theta) in AB's normalization."""
    out = {}
    inv = K(1) / K(2 * L + 1)
    if L >= 1:
        _dacc(out, E(L - 1), K(L) * inv, 0)
    _dacc(out, E(L + 1), K(L + 1) * inv, 0)
    return out


def dens_dxz(L, inner, D1, zC, it, kt, K):
    """m = 1 multipole of x (z - D1) e^{-zC r_C} about B (P_L^1 coefficient / (2L+1)).  g_k = r (T e)_k - D1 e_k are
    the AB-normalized Legendre coefficients of (z - D1) e^{-zC s}; then r (g_{L-1} - g_{L+1}) / (2L+1)."""
    out = {}
    if L == 0:
        return out
    E = _ecache(inner, D1, zC, it, kt, K)
    D1q = K(Fraction(D1))

    def g(k):
        d = {}
        _dacc(d, _tmul(E, k, K), K(1), 1)       # r t  (the z_C factor)
        _dacc(d, E(k), -D1q, 0)                 # - D1
        return d
    inv = K(1) / K(2 * L + 1)
    _dacc(out, g(L - 1), inv, 1)
    _dacc(out, g(L + 1), -inv, 1)               # m = 1 raising step
    return out


def dens_dx2y2(L, inner, D1, zC, it, kt, K):
    """m = 2 multipole of r^2 (1 - t^2) cos(2 phi) e^{-zC s} about B (P_L^2 coefficient / (2L+1)):
    r^2 [e_{L-2}/((2L-1)(2L+1)) - 2 e_L/((2L-1)(2L+3)) + e_{L+2}/((2L+1)(2L+3))]."""
    out = {}
    if L < 2:
        return out
    E = _ecache(inner, D1, zC, it, kt, K)
    _dacc(out, E(L - 2), K(1) / K((2 * L - 1) * (2 * L + 1)), 2)
    _dacc(out, E(L), K(-2) / K((2 * L - 1) * (2 * L + 3)), 2)     # m = 2 middle term
    _dacc(out, E(L + 2), K(1) / K((2 * L + 1) * (2 * L + 3)), 2)
    return out


DENS = {'DX': dens_dxz, 'DD': dens_dx2y2}
M = {'DX': 1, 'DD': 2}


def coefficients(case, Lmax, K, integral):
    """List over L of {key: q}; G_L (DX) or K_L (DD) = 2 pi N_B N_Q sum q c."""
    D1, D2, zB, zC, zA = params(case)
    mu = 2 * zA
    it, kt = AB.bessel_tables(K, Lmax + 3)
    I = AB.Integrals(K)
    co = {'NA': K(1), 'Y': K(-1), 'dY': K(Fraction(zA))}
    out = []
    for L in range(Lmax + 1):
        dens = {inner: DENS[integral](L, inner, D1, zC, it, kt, K) for inner in (True, False)}
        kers = {nr: AP.kernel_factored(L, nr, co, it, kt, K, mu, D2, 0) for nr in (True, False)}
        res = AP.integrate(dens, kers, D1, D2, zB, I, K)
        out.append({k: v for k, v in res.items() if v != 0})
    return out


def _mf(q):
    q = Fraction(q)
    return mp.mpf(q.numerator) / q.denominator


def n_Q(zC, integral):
    """N_xz = sqrt(15/(4 pi)) R_n, N_x2y2 = sqrt(15/(16 pi)) R_n, R_n = sqrt((2 zC)^7 / 720) (current precision)."""
    Rn = mp.sqrt((2 * zC) ** 7 / 720)
    ang = mp.sqrt(15 / (4 * mp.pi)) if integral == 'DX' else mp.sqrt(15 / (16 * mp.pi))   # N_xz / N_x2y2
    return ang * Rn


def norm_const(case, integral):
    """N_B N_Q at the current precision."""
    zB, zC = (_mf(v) for v in params(case)[2:4])
    return mp.sqrt(zB ** 3 / mp.pi) * n_Q(zC, integral)


def numeric(Q, case, integral, dps):
    with mp.workdps(dps):
        tw = 2 * mp.pi * norm_const(case, integral)
        keys = set(k for q in Q for k in q)
        cv = {k: AB.const_value(k) for k in keys}
        return [tw * sum((mp.mpf(int(v.p)) / int(v.q) * cv[k] for k, v in q.items()), mp.mpf(0)) for q in Q]


_EXACT = {}


def exact_Q(case, integral, Lmax):
    key = (case, integral)
    if key not in _EXACT or len(_EXACT[key]) <= Lmax:
        _EXACT[key] = coefficients(case, Lmax, Field(None), integral)
    return _EXACT[key][:Lmax + 1]


def PLm(m, Lmax, x):
    """P_L^m(x) = (1 - x^2)^{m/2} d^m P_L(x), L = 0..Lmax, m = 1, 2, no Condon-Shortley phase (fixed-m recurrence
    (L-m+1) Q_{L+1} = (2L+1) x Q_L - (L+m) Q_{L-1} on Q_L = d^m P_L)."""
    out = [mp.mpf(0)] * (Lmax + 1)
    fac = (1 - x * x) ** (mp.mpf(m) / 2)
    q0, q1 = mp.mpf(0), mp.mpf(1) if m == 1 else mp.mpf(3)
    out[m] = fac * q1
    for L in range(m, Lmax):
        q0, q1 = q1, ((2 * L + 1) * x * q1 - (L + m) * q0) / (L - m + 1)
        out[L + 1] = fac * q1
    return out


# ---------------------------------------------------------------------------------------------------------------
# Independent route: direct projection of the m = 1, 2 density multipoles (port of the record's AX0 (ii) check)
# ---------------------------------------------------------------------------------------------------------------


def dmP(m, L, t):
    """d^m P_L(t), m = 1, 2, by the fixed-m recurrence (not P_L' from the plain recurrence)."""
    if L < m:
        return mp.mpf(0)
    q0, q1 = mp.mpf(0), mp.mpf(1) if m == 1 else mp.mpf(3)
    for n in range(m, L):
        q0, q1 = q1, ((2 * n + 1) * t * q1 - (n + m) * q0) / (n - m + 1)
    return q1


def _ang(func, r, D, nodes):
    """int_{-1}^{1} func(s, t) dt with s = |r - D e_z| as the variable: Gauss-Legendre in s."""
    a, b = abs(r - D), r + D
    h, c = (b - a) / 2, (b + a) / 2
    tot = mp.mpf(0)
    for u, w in nodes:
        s = c + h * u
        t = (r * r + D * D - s * s) / (2 * r * D)
        tot += w * func(s, t) * s
    return tot * h / (r * D)


def direct_projection(case, L, integral, dps, deg=5):
    """G_L / K_L = 4 pi/(2L+1) int r^2 rho^m_L(r) v_L(r) dr with
        DX: rho^1_L = (2L+1)/(2L(L+1)) N e^{-zB r} r int e^{-zC s} (r t - D1) (1-t^2) P_L'(t) dt
        DD: rho^2_L = (2L+1)(L-2)!/(2 (L+2)!) N e^{-zB r} r^2 int e^{-zC s} (1-t^2)^2 P_L''(t) dt
    (the m-th associated-Legendre projection, sqrt(1-t^2)^m P_L^m = (1-t^2)^m d^m P_L), v_L the Legendre multipole
    of the exact 1s potential of A about B; Gauss-Legendre (3 * 2^(deg-1) nodes) in s, tanh-sinh in r."""
    m = M[integral]
    with mp.workdps(dps):
        D1, D2, zB, zC, zA = (_mf(v) for v in params(case))
        N = norm_const(case, integral)
        nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec)
        if m == 1:
            c = (2 * L + 1) / mp.mpf(2 * L * (L + 1))

            def rho(r):
                return c * N * mp.exp(-zB * r) * r * _ang(
                    lambda s, t: mp.exp(-zC * s) * (r * t - D1) * (1 - t * t) * dmP(1, L, t), r, D1, nodes)
        else:
            c = (2 * L + 1) * mp.factorial(L - 2) / (2 * mp.factorial(L + 2))

            def rho(r):
                return c * N * mp.exp(-zB * r) * r * r * _ang(
                    lambda s, t: mp.exp(-zC * s) * (1 - t * t) ** 2 * dmP(2, L, t), r, D1, nodes)

        def vL(r):
            return (2 * L + 1) / mp.mpf(2) * _ang(lambda s, t: AB.potential_1s(s, zA) * AB._legendre(L, t),
                                                  r, D2, nodes)
        return mp.quad(lambda r: 4 * mp.pi / (2 * L + 1) * r * r * rho(r) * vL(r), [0, D1, D2, D2 + 4, mp.inf])


# ---------------------------------------------------------------------------------------------------------------
# Toys: (1 - x^2)^{m/2} (x* - x)^alpha in the P_L^m basis (record AX control T: Gauss-Legendre projection of
# (x* - x)^alpha, then m antiderivatives, since (1 - x^2)^{m/2} h^(m) = sum e_L P_L^m when h = sum e_L P_L)
# ---------------------------------------------------------------------------------------------------------------


def x_star(case):
    D1, D2 = params(case)[:2]
    return (D1 ** 2 + D2 ** 2) / (2 * D1 * D2)


def _anti_list(c):
    g = [mp.mpf(0)] * len(c)
    for L in range(1, len(c) - 1):
        g[L] = c[L - 1] / (2 * L - 1) - c[L + 1] / (2 * L + 3)
    return g


_TOYS = {}


def toy(case, alpha, Lmax, dps=60, deg=7):
    """{1: P_L^1 coefficients of sqrt(1 - x^2)(x* - x)^alpha, 2: P_L^2 coefficients of (1 - x^2)(x* - x)^alpha}."""
    key = (case, alpha, Lmax, dps, deg)
    if key not in _TOYS:
        with mp.workdps(dps):
            xs, al = _mf(x_star(case)), _mf(alpha)
            nodes = mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec)
            n = Lmax + 4
            c = [mp.mpf(0)] * (n + 1)
            for x, w in nodes:
                hv = w * (xs - x) ** al
                p0, p1 = mp.mpf(1), x
                c[0] += hv
                c[1] += hv * x
                for k in range(1, n):
                    p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
                    c[k + 1] += hv * p1
            c = [(2 * k + 1) / mp.mpf(2) * v for k, v in enumerate(c)]
            a1 = _anti_list(c)
            _TOYS[key] = {1: a1, 2: _anti_list(a1)}
    return _TOYS[key]


def rho(case):
    xs = _mf(x_star(case))
    return 1 / (xs + mp.sqrt(xs * xs - 1))


def beta_eff(F, case, L):
    """log|g_{L+1}/g_{L-1}| / log((L+1)/(L-1)), g_L = F_L / rho^L."""
    r = rho(case)
    return float(mp.log(abs((F[L + 1] / r ** (L + 1)) / (F[L - 1] / r ** (L - 1)))) / mp.log(mp.mpf(L + 1) / (L - 1)))


ALPHAS = [Fraction(k, 2) for k in (3, 5, 7, 9, 11)]


def window(case, L, m, alpha=Fraction(7, 2)):
    """(lo, hi): midpoints between the toy beta_eff(L) at alpha and those at alpha + 1 (lo) and alpha - 1 (hi)."""
    t = {a: beta_eff(toy(case, a, FAST_L)[m], case, L) for a in ALPHAS}
    i = ALPHAS.index(alpha)
    return (t[ALPHAS[i]] + t[ALPHAS[i + 1]]) / 2, (t[ALPHAS[i]] + t[ALPHAS[i - 1]]) / 2


# ---------------------------------------------------------------------------------------------------------------
# Shared data: both integrals, both cases, L <= FAST_L over Q (~4 s per engine)
# ---------------------------------------------------------------------------------------------------------------

FAST_L, DPS = 81, 600    # the split cancels ~220 digits at A2, L = 60 (measured: 200 vs 700 digits); dps + 100 checked
CASES = ('A1', 'A2')
_SERIES = {}


def series(case, integral, dps=DPS):
    key = (case, integral, dps)
    if key not in _SERIES:
        _SERIES[key] = numeric(exact_Q(case, integral, FAST_L), case, integral, dps)
    return _SERIES[key]


# ---------------------------------------------------------------------------------------------------------------
# (0) normalizations and the vanishing yz / xy integrals
# ---------------------------------------------------------------------------------------------------------------


def test_normalizations_are_unit_norm():
    """int |N_Q Q(x, y, z) e^{-zC r}|^2 d^3r = 1 for Q = xz and x^2 - y^2 (zC = 1 and 6/5), by quadrature of the
    Cartesian polynomial on the sphere (no 4 pi/15 closed form).  Rejects N_xz and N_x2y2 swapped (norm 4 or 1/4), a
    wrong R_n (e.g. (2 zC)^5 or 720 -> 120), a missing sqrt."""
    with mp.workdps(30):
        for zC in (mp.mpf(1), mp.mpf(6) / 5):
            rad = mp.quad(lambda r: r ** 6 * mp.exp(-2 * zC * r), [0, mp.inf])
            for integral, poly in (('DX', lambda st, ct, cp, sp: (st * cp * ct) ** 2),
                                   ('DD', lambda st, ct, cp, sp: (st * st * (cp * cp - sp * sp)) ** 2)):
                angf = mp.quad(lambda th, ph: poly(mp.sin(th), mp.cos(th), mp.cos(ph), mp.sin(ph)) * mp.sin(th),
                               [0, mp.pi], [0, 2 * mp.pi])
                assert abs(n_Q(zC, integral) ** 2 * angf * rad - 1) < mp.mpf('1e-20'), (integral, zC)


def _geometry(case, x):
    """C and A in Cartesian coordinates about B = 0: C = (0, 0, D1), A = D2 (sin g, 0, cos g), x = cos g."""
    D1, D2 = (_mf(v) for v in params(case)[:2])
    return (mp.mpf(0), mp.mpf(0), D1), (D2 * mp.sqrt(1 - x * x), mp.mpf(0), D2 * x)


def _azimuthal(case, x, Qf, r, th):
    """int_0^{2 pi} e^{-zB r} Q(r - C) e^{-zC |r - C|} V_A(|r - A|) dphi at (r, theta) about B."""
    zB, zC, zA = (_mf(v) for v in params(case)[2:])
    C, A = _geometry(case, x)

    def f(ph):
        p = (r * mp.sin(th) * mp.cos(ph), r * mp.sin(th) * mp.sin(ph), r * mp.cos(th))
        dc = [p[i] - C[i] for i in range(3)]
        sc = mp.sqrt(sum(v * v for v in dc))
        sa = mp.sqrt(sum((p[i] - A[i]) ** 2 for i in range(3)))
        return mp.exp(-zB * r) * Qf(*dc) * mp.exp(-zC * sc) * AB.potential_1s(sa, zA)
    return mp.quad(f, [0, mp.pi / 2, mp.pi, 3 * mp.pi / 2, 2 * mp.pi])


SHAPES = {'xz': lambda X, Y, Z: X * Z, 'x2y2': lambda X, Y, Z: X * X - Y * Y,
          'yz': lambda X, Y, Z: Y * Z, 'xy': lambda X, Y, Z: X * Y}


def test_yz_and_xy_vanish_by_reflection():
    """At A1, x = 1/2 and A2, x = -1/3, three (r, theta) about B each: the azimuthal integral of the 3d_yz and 3d_xy
    integrands is below 1e-25 of the 3d_xz one, while those of 3d_xz and 3d_x2-y2 exceed 1e-5 (smallest measured
    6.9e-5).  The yz and xy
    densities are odd under y -> -y and V_A is even with A in the xz-plane, so the full integrals vanish (record AX:
    R2 exactly 0).  Rejects an A placed out of the xz-plane and a yz or xy shape that is not odd in y.  This guards
    the symmetry argument; the engine has no yz or xy path to test."""
    with mp.workdps(30):
        for case, xq in (('A1', Fraction(1, 2)), ('A2', Fraction(-1, 3))):
            for r, th in ((mp.mpf('0.7'), mp.mpf('0.4')), (mp.mpf('1.9'), mp.mpf('1.3')),
                          (mp.mpf('3.1'), mp.mpf('2.2'))):
                v = {k: _azimuthal(case, _mf(xq), f, r, th) for k, f in SHAPES.items()}
                assert abs(v['xz']) > 1e-5 and abs(v['x2y2']) > 1e-5, (case, r, th, v)
                assert abs(v['yz']) < mp.mpf('1e-25') * abs(v['xz']), (case, r, th, v['yz'])
                assert abs(v['xy']) < mp.mpf('1e-25') * abs(v['xz']), (case, r, th, v['xy'])


# ---------------------------------------------------------------------------------------------------------------
# (1) per-channel values against the direct projection
# ---------------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize('integral,case,L', [('DX', 'A1', 3), ('DX', 'A2', 7), ('DD', 'A1', 5), ('DD', 'A2', 2)])
def test_channel_matches_direct_projection(integral, case, L):
    """G_L (DX) / K_L (DD) of the exact split, evaluated at 300 digits (it cancels fewer than 30 digits at L <= 7),
    equal the direct projection at 30 digits within relative 1e-26 (measured 7e-32 to 2e-31; at 22 digits the same
    projection gives 1e-23 here but 2e-16 at A2 DD L = 16, the precision that VOIDed record AX at L = 12).  The
    projection uses (r t - D1)(1 - t^2) P_L' and (1 - t^2)^2 P_L'' as is: no T recurrence, no raising identity.
    Rejects a dropped z_C factor or D1 term, a sign or index error in either raising step, a Condon-Shortley sign on
    the m = 1 step (each moves a channel by far more than 1e-6)."""
    ex = numeric(exact_Q(case, integral, L), case, integral, 300)[L]
    q = direct_projection(case, L, integral, 30)
    with mp.workdps(60):
        rel = abs(q / ex - 1)
        assert rel < mp.mpf('1e-26'), (integral, case, L, mp.nstr(q, 30), mp.nstr(ex, 30), mp.nstr(rel, 3))


@pytest.mark.slow
@pytest.mark.parametrize('integral', ['DX', 'DD'])
def test_channel_matches_direct_projection_record_points(integral):
    """Record AX2b's points L = 16 and 25, both cases: the split at 300 digits (it cancels ~1e65 at L = 20) against
    the projection at 44 digits, relative 1e-36 (record AX2b-0: 8e-46 to 2.9e-40; its bar control |q44/q56 - 1| <=
    2.9e-40 at L = 25)."""
    for case in CASES:
        ex = numeric(exact_Q(case, integral, 25), case, integral, 300)
        for L in (16, 25):
            q = direct_projection(case, L, integral, 44)
            with mp.workdps(80):
                assert abs(q / ex[L] - 1) < mp.mpf('1e-36'), (integral, case, L, mp.nstr(abs(q / ex[L] - 1), 3))


# ---------------------------------------------------------------------------------------------------------------
# (2) the summed series at the record's eight points
# ---------------------------------------------------------------------------------------------------------------

# record AX, control R2: Gaussian-fit referee, n = 34 (tolerance 1e-10)
R2_34 = {('A1', 'DX', '1/2'): '-0.03301838760929426', ('A1', 'DX', '0'): '-0.04177923697376483',
         ('A2', 'DX', '-1/3'): '-0.050507409529967365', ('A2', 'DX', '4/5'): '-0.018247628218714944',
         ('A1', 'DD', '1/2'): '0.011120171409288062', ('A1', 'DD', '0'): '0.010679184288401823',
         ('A2', 'DD', '-1/3'): '0.012344129315206417', ('A2', 'DD', '4/5'): '0.010238840327312952'}
# record AX2b, P1: exact split at 2000 digits summed to L = 240 (debug/data/aadx_p1p3_2000.json, sha256 c9151de7...)
P1_2000 = {('A1', 'DX', '1/2'): '-0.0330183876091728373514774989021',
           ('A1', 'DX', '0'): '-0.0417792369733423956801105007282',
           ('A2', 'DX', '-1/3'): '-0.0505074095304871411360006408648',
           ('A2', 'DX', '4/5'): '-0.0182476282186817620221659846985',
           ('A1', 'DD', '1/2'): '0.011120171409256727037653466317',
           ('A1', 'DD', '0'): '0.0106791842881800061275047154628',
           ('A2', 'DD', '-1/3'): '0.0123441293154449715358575764668',
           ('A2', 'DD', '4/5'): '0.0102388403273425571976224826105'}
TRUNC = mp.mpf('1e-22')    # sum to L = 81 minus the L = 240 value: measured at most 3.0e-24 (A1 DX, x = 1/2)


@pytest.mark.parametrize('integral', ['DX', 'DD'])
def test_summed_series_matches_referee_eight_points(integral):
    """Sum_{L <= 81} G_L P_L^1(x) (DX) or K_L P_L^2(x) (DD), P_L^m without the Condon-Shortley phase, both cases,
    both record angles: within 1e-10 of the R2 Gaussian referee (record: at most 5.2e-13) and within TRUNC of the
    record's 2000-digit P1 values.  This clause tests DX's z_C factor (record AX: the decay cannot), the
    normalizations N_xz, N_x2y2 and the 2 pi addition-theorem factor.  Rejects a series in P_L instead of P_L^m, a
    Condon-Shortley sign on P_L^1, N_xz and N_x2y2 swapped, a referee value off by 2e-10."""
    m = M[integral]
    for case in CASES:
        F = series(case, integral)
        with mp.workdps(DPS):
            for (c_, i_, xs), ref in R2_34.items():
                if (c_, i_) != (case, integral):
                    continue
                Pm = PLm(m, FAST_L, _mf(Fraction(xs)))
                v = sum(F[L] * Pm[L] for L in range(m, FAST_L + 1))
                assert abs(v - mp.mpf(ref)) < 1e-10, (case, integral, xs, mp.nstr(v, 20))
                assert abs(v - mp.mpf(P1_2000[(case, integral, xs)])) < TRUNC, (case, integral, xs, mp.nstr(v, 30))


# ---------------------------------------------------------------------------------------------------------------
# (3) the constant set equals AB's 21
# ---------------------------------------------------------------------------------------------------------------


@pytest.mark.parametrize('case', CASES)
def test_constant_set_is_ABs_21_no_gamma(case):
    """Mod 2^61 - 1: G_L for 1 <= L <= 20 and K_L for 2 <= L <= 20 carry exactly AB's 21 constants at the same case
    (A1: AB's frozen EXPECTED_A1; A2: AB's own engine), no Euler gamma at any L; G_0, K_0, K_1 are empty; and every
    coefficient for L <= 25 is unchanged when the internal L_max goes 30 -> 40.  Rejects a spurious gamma, a lost or
    extra constant, a set that changes with L, and a top-of-range truncation of the T or raising recurrences."""
    K = Field(PRIME)
    ab = AB.EXPECTED_A1 if case == 'A1' else AP._ab_set(case)
    assert len(ab) == 21
    for integral, L0 in (('DX', 1), ('DD', 2)):
        f = coefficients(case, 30, K, integral)
        assert all(f[L] == {} for L in range(L0)), (integral, case)
        for L in range(L0, 21):
            S = set(f[L])
            assert not [k for k in S if k[1] == 'g'], (integral, case, L)
            assert S == ab, (integral, case, L, sorted(map(repr, S ^ ab)))
        assert f[:26] == coefficients(case, 40, K, integral)[:26], (integral, case)


# ---------------------------------------------------------------------------------------------------------------
# (4) decay: the 7/2 branch in P_L^1 (L^-5) and P_L^2 (L^-6)
# ---------------------------------------------------------------------------------------------------------------

# record AX, control T (700 nodes / 120 digits): beta_eff(40) of the toys
T40 = {('A1', 1, Fraction(7, 2)): -5.937, ('A2', 1, Fraction(7, 2)): -5.766,
       ('A1', 1, Fraction(9, 2)): -7.325, ('A2', 1, Fraction(9, 2)): -7.091,
       ('A1', 2, Fraction(7, 2)): -7.325, ('A2', 2, Fraction(7, 2)): -7.091,
       ('A1', 2, Fraction(9, 2)): -8.778, ('A2', 2, Fraction(9, 2)): -8.474}


def test_toys_reproduce_record_control_T():
    """The in-test toys' beta_eff(40) against the record's control T to its printed 3 decimals.  Rejects a wrong
    antiderivative recurrence, one antiderivative where m = 2 needs two, a wrong x*, a missing (2k+1)/2."""
    for (case, m, a), v in T40.items():
        assert abs(beta_eff(toy(case, a, FAST_L)[m], case, 40) - v) < 6e-4, (case, m, a)


def test_decay_seven_halves_branch_in_PLm():
    """beta_eff(L) at L = 40, 60, 80, both cases: G_L between the midpoints of the P_L^1 toys around the 7/2 branch,
    K_L between those of the P_L^2 toys (record at L = 40: DX -5.959 / -5.804, DD -7.175 / -6.964; P3 at L = 220:
    -5.185 / -5.150, -6.241 / -6.190).  A cancelled 7/2 coefficient (the 9/2 markers) or a 5/2 branch falls outside.
    One sign (positive) for 5 <= L <= 81.  Precision: beta at dps and dps + 100 agree to 1e-12.  Not a test of DX's
    z_C factor (record AX: AQ's p density decays on the same branch); the summed-series test is."""
    for integral in ('DX', 'DD'):
        m = M[integral]
        for case in CASES:
            F, F_hi = series(case, integral), series(case, integral, DPS + 100)
            for L in (40, 60, 80):
                b = beta_eff(F, case, L)
                assert abs(b - beta_eff(F_hi, case, L)) < 1e-12, (integral, case, L)
                lo, hi = window(case, L, m)
                assert lo < b < hi, (integral, case, L, b, (lo, hi))
            assert all(v > 0 for v in F[5:FAST_L + 1]), (integral, case)


# TIMINGS (2026-11-16, this machine): fast tests ~47 s in all (four engines to L = 81 over Q ~3 s each, the series at
# 600 and 700 digits, summed tests ~12 s each; direct projections 1.5-3 s each; constant sets mod p ~2.7 s per case;
# norm quadrature 4 s; toys 2 s).  Slow tests (--slow) ~70 s more: the L = 16, 25 projections at 44 digits.
#
# FIRE-TEST PLANTS (python debug/qa/fire_test.py tests/test_paper59_aa_bc_d_mneq0.py --plant "OLD=>NEW" -k <sel>)
# 2026-11-16, written and fired as a separate step after the guard.  CH = "test_channel_matches_direct_projection and
# not record_points".  27 plants fired; 4 did not, each expected and covered elsewhere (marked *).
#   z_C factor dropped (DX density = AQ's x e^{-zC s}; two plants together)
#        "_dacc(d, _tmul(E, k, K), K(1), 1)       # r t  (the z_C factor)=>_dacc(d, E(k), K(1), 0)  # zC dropped"
#        "_dacc(d, E(k), -D1q, 0)                 # - D1=>pass"
#        -k CH FIRED; -k summed_series FIRED; -k decay_seven DID NOT FIRE*  (the record's disclosed weak point: AQ's p
#        density decays on the same 7/2 branch in P_L^1, so the decay cannot see z_C; the values can)
#   r t term dropped       "# r t  (the z_C factor)" line, K(1) -> K(0)        -k "CH or summed_series"   FIRED
#   D1 term dropped        "_dacc(d, E(k), -D1q, 0)=>_dacc(d, E(k), 0 * D1q, 0)"   -k CH FIRED; -k summed_series FIRED;
#                          -k decay_seven DID NOT FIRE* (same branch)
#   m = 1 raising sign     "_dacc(out, g(L + 1), -inv, 1)=>_dacc(out, g(L + 1), inv, 1)"     -k CH   FIRED
#   m = 1 raising index    "_dacc(out, g(L + 1), -inv, 1)=>_dacc(out, g(L), -inv, 1)"        -k CH   FIRED
#   m = 1 step with a Condon-Shortley phase (both terms' signs flipped)          -k CH FIRED; -k summed_series FIRED
#   T recurrence weight    "_dacc(out, E(L + 1), K(L + 1) * inv, 0)=>... K(L) * inv, 0)"      -k CH   FIRED
#   m = 2 middle sign      "K(-2) / K((2 * L - 1) * (2 * L + 3)), 2)=>K(2) / ..."             -k CH   FIRED
#   m = 2 index L+2 -> L+1 "_dacc(out, E(L + 2), K(1)=>_dacc(out, E(L + 1), K(1)"             -k CH   FIRED
#   m = 2 density from a single m = 1 step (E(L-1), -E(L+1), middle term 0; three plants)
#                                                                                -k CH FIRED; -k summed_series FIRED
#     (applying the m = 1 step of dens_px literally twice is NOT a wrong answer: it reduces exactly to the three-term
#     m = 2 formula, since the 1/(2l+1) of the normalization cancels that of the identity)
#   N_xz / N_x2y2 swapped  "if integral == 'DX' else mp.sqrt(15 / (16 * mp.pi))=>if integral == 'DD' else ..."
#                          -k summed_series FIRED; -k normalizations FIRED; -k CH DID NOT FIRE* (the projection shares
#                          the normalization, as record AX2b discloses for AX2b-0)
#   R_n 720 -> 120         "Rn = mp.sqrt((2 * zC) ** 7 / 720)=>... / 120)"   -k "normalizations or summed_series"   FIRED
#   P_L^m with a Condon-Shortley phase  "fac = (1 - x * x) ** (mp.mpf(m) / 2)=>fac = (-1) ** m * ..."
#                          -k "summed_series and DX" FIRED; -k "summed_series and DD" DID NOT FIRE* ((-1)^2 = 1: not
#                          a wrong answer for m = 2)
#   P_L for P_L^m          "v = sum(F[L] * Pm[L] ...=>v = sum(F[L] * AB._legendre(L, _mf(Fraction(xs))) ..."
#                                                                                -k summed_series   FIRED
#   referee off 2e-10      "('A1', 'DD', '1/2'): '0.011120171409288062'=>... '0.011120171609288062'"   -k summed_series FIRED
#   spurious Euler gamma   "out.append({k: v ...})=>res[(Fraction(0), 'g', None)] = K(1); out.append(...)"
#                                                                                -k constant_set   FIRED
#   m = 2 toy with one antiderivative  "_TOYS[key] = {1: a1, 2: _anti_list(a1)}=>_TOYS[key] = {1: a1, 2: a1}"
#                                                       -k toys_reproduce FIRED; -k decay_seven FIRED
#   window on the 9/2 branch / on the 5/2 branch  "alpha=Fraction(7, 2)):=>alpha=Fraction(9, 2)):" (and 5/2)
#                                                                                -k decay_seven   FIRED, FIRED
#   A out of the xz-plane  "mp.mpf(0), D2 * x)=>D2 * mp.mpf('0.1'), D2 * x)"   -k vanish   FIRED
#   yz shape not odd in y  "'yz': lambda X, Y, Z: Y * Z=>'yz': lambda X, Y, Z: (Y + X) * Z"   -k vanish   FIRED
