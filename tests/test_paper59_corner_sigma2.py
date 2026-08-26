"""Backing test for Paper 59 (sec:modular): the (s,t)->(0,0) corner of the
collinear three-center Bessel-moment observable T2 is a rho^{3/2} (rho=s+t)
non-analyticity, and the Duffy angular split + radial substitution rho=sigma^2
renders the corner integrand ANALYTIC in sigma, restoring SPECTRAL
Gauss-Legendre convergence -- the method that lifted the value past its earlier
~16-digit corner-under-resolution ceiling.

Self-contained (does NOT import the transient debug/ drivers, per the clean-room
policy): reimplements the minimal collinear integrand and evaluates the corner
triangle T1 = int_{s+t<=delta} J two ways -- (A) raw Duffy + GL in rho
[baseline, non-analytic rho^{3/2}, only ALGEBRAIC convergence] and (C) Duffy +
rho=sigma^2 + GL [predicted SPECTRAL]. The test asserts that (C)'s successive
degree-refinement differences shrink SUPER-geometrically (spectral) while (A)'s
do not, and pins the corner value.

    T2 = (8/pi) int_[0,1]^2 J(s,t) ds dt,  J(s,t) = int_0^inf dk j0(k(s+t)) P(s,k)P(t,k)
    P(x,k) = c e^{-Del}(1/Del^3 + 3/Del^4 + 3/Del^5),  c=x(1-x), Del=sqrt(c k^2+1)
"""
import mpmath as mp
import pytest


def _P(x, k):
    c = x * (1 - x)
    Del = mp.sqrt(c * k * k + 1)
    return c * mp.e ** (-Del) * (1 / Del ** 3 + 3 / Del ** 4 + 3 / Del ** 5)


def _k_grid(theta_max, panel_w, n_panel):
    """Branch-point-regularizing k = 2 sinh(theta) paneled Gauss-Legendre grid
    (fixed Im(theta)=pi/2 gap for c<=1/4; compresses the large-k reach needed at
    small c)."""
    rule = mp.calculus.quadrature.GaussLegendre(mp.mp)
    base = rule.calc_nodes(n_panel, mp.mp.prec)
    nodes = []
    lo = mp.mpf(0)
    while lo < theta_max - mp.mpf('1e-30'):
        hi = min(lo + panel_w, theta_max)
        half, mid = (hi - lo) / 2, (hi + lo) / 2
        for x, w in base:
            th = mid + half * x
            nodes.append((2 * mp.sinh(th), w * half * 2 * mp.cosh(th)))
        lo = hi
    return nodes


def _J(s, t, kn):
    b = s + t
    tot = mp.mpf(0)
    for k, w in kn:
        kb = k * b
        j0 = mp.sin(kb) / kb if kb > mp.mpf('1e-40') else mp.mpf(1) - kb * kb / 6
        tot += w * j0 * _P(s, k) * _P(t, k)
    return tot


def _gl(deg):
    return mp.calculus.quadrature.GaussLegendre(mp.mp).calc_nodes(deg, mp.mp.prec)


def _T1_sigma2(delta, deg, kn):
    """(C) corner triangle {s+t<=delta} via Duffy s=rho*a, t=rho*(1-a), then
    rho=sigma^2: integrand 2 sigma^3 J on (sigma,a) in [0,sqrt(delta)]x[0,1]
    (a folded to [0,1/2], x2)."""
    ss, aa = _gl(deg), _gl(deg)
    sig_max = mp.sqrt(delta); half = sig_max / 2
    tot = mp.mpf(0)
    for xs, ws in ss:
        sig = half * (xs + 1)
        rho = sig * sig
        pref = (half * ws) * 2 * sig ** 3
        acc = mp.mpf(0)
        for xa, wa in aa:
            a = (xa + 1) / 4
            acc += (wa / 4) * _J(rho * a, rho * (1 - a), kn)
        tot += pref * 2 * acc
    return tot


def _T1_raw(delta, deg, kn):
    """(A) corner triangle via raw Duffy + GL in rho (no substitution)."""
    rr, aa = _gl(deg), _gl(deg)
    half = delta / 2
    tot = mp.mpf(0)
    for xr, wr in rr:
        rho = half * (xr + 1)
        pref = (half * wr) * rho
        acc = mp.mpf(0)
        for xa, wa in aa:
            a = (xa + 1) / 4
            acc += (wa / 4) * _J(rho * a, rho * (1 - a), kn)
        tot += pref * 2 * acc
    return tot


@pytest.mark.slow
def test_paper59_corner_sigma2_is_spectral():
    """rho=sigma^2 makes the corner integrand analytic => GL SPECTRAL, whereas
    raw-rho GL on the rho^{3/2} integrand is only ALGEBRAIC. Assert the sigma^2
    scheme accelerates super-geometrically and does so FAR faster than raw-rho,
    and that both schemes converge to the SAME corner value."""
    mp.mp.dps = 30
    delta = mp.mpf('0.05')
    kn = _k_grid(mp.mpf(16), mp.mpf(2), 4)   # M=192

    vC = [_T1_sigma2(delta, d, kn) for d in (3, 4, 5)]
    vA = [_T1_raw(delta, d, kn) for d in (3, 4, 5)]
    dC1, dC2 = abs(vC[1] - vC[0]), abs(vC[2] - vC[1])
    dA1, dA2 = abs(vA[1] - vA[0]), abs(vA[2] - vA[1])

    # (C) spectral: the second diff is orders below the first (accelerating).
    # Measured: dC2/dC1 ~ 1.0e-4  vs raw dA2/dA1 ~ 8.7e-3.
    assert dC2 < dC1 * mp.mpf('1e-3'), f"sigma^2 not spectral: {dC1} -> {dC2}"
    # sigma^2 accelerates an order of magnitude faster than the algebraic raw-rho.
    assert (dC2 / dC1) < (dA2 / dA1) * mp.mpf('0.05'), \
        f"sigma^2 not >> raw acceleration: {dC2/dC1} vs {dA2/dA1}"
    # both schemes converge to the SAME corner value (independent cross-check).
    assert abs(vC[2] - vA[2]) < mp.mpf('1e-10')


@pytest.mark.slow
def test_paper59_corner_value():
    """Pin the corner-triangle value T1(delta=0.05) (Paper 59 sec:modular
    intermediate), independent of the debug/ evaluator. The ~13-digit pin is
    stable across k-grid size (M>=192 agree there)."""
    mp.mp.dps = 30
    delta = mp.mpf('0.05')
    kn = _k_grid(mp.mpf(16), mp.mpf(2), 4)    # M=192
    T1 = _T1_sigma2(delta, 5, kn)
    ref = mp.mpf('5.075171796722616e-6')      # certified, ~13-digit pin
    assert abs(T1 - ref) < mp.mpf('1e-18'), mp.nstr(T1, 20)


# ---------------------------------------------------------------------------
# The ASSEMBLED collinear value (corner + bulk), pinned against the paper's
# canonical headline 0.3953557659017139641 (Paper 59 sec:modular).  The
# spectral corner T1 (sigma^2 Duffy) is glued to the smooth bulk on a tiling
# [0,1]^2 = T1{s+t<=d} + T2far{[0,d]^2, s+t>d} + 2*RectA + BigSquare; the
# smooth pieces use cosine-clustered Gauss-Legendre.  This is an INDEPENDENT
# code path from the debug/ hp_final evaluator that produced the headline; the
# two agree to ~12 digits at this cost (and to ~16 digits at deg=5, confirming
# the headline value) -- the assembled witness the earlier ~6-digit float64
# test (test_cosmic_galois_integrated_value) lacked.
# ---------------------------------------------------------------------------

# canonical collinear value to the ~19 digits reachable via the spectral corner
# (debug/routeC_hp_final.py; ceiling ~20, set by the pointwise corner k-grid).
_V_COLLINEAR = mp.mpf('0.3953557659017139641')


def _T2far(delta, deg, kn):
    """Far triangle {s,t in [0,delta], s+t>delta} via reflected Duffy
    s=delta-rho*a, t=delta-rho*(1-a) (smooth: s,t near delta, off the corner)."""
    rr, aa = _gl(deg), _gl(deg); half = delta / 2; tot = mp.mpf(0)
    for xr, wr in rr:
        rho = half * (xr + 1); pref = (half * wr) * rho; acc = mp.mpf(0)
        for xa, wa in aa:
            a = (xa + 1) / 4
            acc += (wa / 4) * _J(delta - rho * a, delta - rho * (1 - a), kn)
        tot += pref * 2 * acc
    return tot


def _clustered(deg, lo, hi):
    """GL in u pushed through x=(1-cos(pi u))/2 (clusters BOTH ends, matching the
    c=x(1-x) edge zeros), affine to [lo,hi]; folds in the dx/du and span Jacobians."""
    span = hi - lo; out = []
    for xu, wu in _gl(deg):
        u = (xu + 1) / 2; wU = wu / 2
        x = (1 - mp.cos(mp.pi * u)) / 2
        dxdu = (mp.pi / 2) * mp.sin(mp.pi * u)
        out.append((lo + span * x, wU * dxdu * span))
    return out


def _rect(deg, slo, shi, tlo, thi, kn):
    ns, nt = _clustered(deg, slo, shi), _clustered(deg, tlo, thi); tot = mp.mpf(0)
    for s, ws in ns:
        acc = mp.mpf(0)
        for t, wt in nt:
            acc += wt * _J(s, t, kn)
        tot += ws * acc
    return tot


@pytest.mark.slow
def test_paper59_assembled_collinear_value():
    """[MEASURED] Assemble the FULL collinear T2 (corner triangle + smooth bulk)
    self-contained and check it reproduces the paper's canonical headline
    0.3953557659017139641 (Paper 59 sec:modular) to ~12 digits.  This is the
    assembled witness beyond the corner sub-triangle: an independent tiling
    (vs debug/hp_final) that agrees with the headline, resolving the earlier
    ~6-digit-only backing.  deg=4 / M=384 gives |T2 - headline| ~ 1.5e-13; the
    same tiling at deg=5 reaches ~1e-16 (a 16-digit cross-validation of the
    headline in the sprint record, too slow to pin as a routine test)."""
    mp.mp.dps = 30
    delta = mp.mpf('0.05')
    kn = _k_grid(mp.mpf(16), mp.mpf(2), 5)     # M=384
    deg = 4
    T1 = _T1_sigma2(delta, deg, kn)
    T2far = _T2far(delta, deg, kn)
    RectA = _rect(deg, mp.mpf(0), delta, delta, mp.mpf(1), kn)
    BigSq = _rect(deg, delta, mp.mpf(1), delta, mp.mpf(1), kn)
    T2 = (8 / mp.pi) * (T1 + T2far + 2 * RectA + BigSq)
    assert abs(T2 - _V_COLLINEAR) < mp.mpf('1e-12'), mp.nstr(T2, 20)


# ---------------------------------------------------------------------------
# The guarded PSLQ-NEGATIVE (Paper 59 sec:modular / sec:bessel_algebra).
#
# What is ROBUSTLY DECIDABLE at the reachable precision.  The collinear value is
# known to only ~15-16 cross-confirmed / ~19 best-estimate digits (extending it
# is the ~40-digit collaboration frontier, not a test-time computation), and
# mpmath.pslq needs >=16 digits to run at all.  So a DECISIVE integer-relation
# search is possible ONLY over a SMALL basis: the disc-4 CM ring
#     {1, pi, varpi=K(1/2)=Gamma(1/4)^2/(4 sqrt pi), pi^2, varpi^2, pi*varpi}
# (weight<=2, n=6 -- decidable at ~16 digits).  We fit the RAW period W = V*pi/8
# (the (8/pi) Coulomb prefactor stripped, so W is the natural period; a Laurent
# relation in V is a polynomial relation in W).  A genuine closure is a
# target-coeff-nonzero, LOW-height (<=40) integer relation.  The NEGATIVE: W's
# relation height is ~140-390 (far above low), a magnitude-matched structureless
# decoy matches (also high-height), and a planted low-height combination IS caught
# (positive control -- the detector is not blind).
#
# The WIDER weight-3 / disc-8-inclusive negative (the ~20-element ring of
# debug/routeC_pslq_v2.py) is a bounded, DRIVER-OBSERVED exclusion that is NOT
# decidable at ~19 digits: that basis is over-determined at reachable precision
# (the spurious low-height hits are a DIFFERENT vector at each precision), so a
# decisive test needs ~40 digits.  It is deliberately NOT asserted here.  Paper 59
# sec:modular/sec:bessel_algebra tier those legs as bounded, "~40 digits for a
# definitive decision" (the collaboration frontier).
#
# HONEST TIER: a bounded search-negative, NOT a transcendence proof.
# ---------------------------------------------------------------------------

_HEIGHT_LOW = 40          # "low-height" threshold; reuses _V_COLLINEAR above


def test_paper59_pslq_negative_disc4_guarded():
    """[MEASURED search-negative, disc-4 weight<=2] The raw period W = V*pi/8 is
    NOT a low-height (height<=40) integer combination of the disc-4 CM ring
    {1, pi, varpi, pi^2, varpi^2, pi*varpi} -- its relation height is ~140-390,
    decisively above low -- and a magnitude-matched structureless decoy matches
    (also high-height).  A planted low-height combination (2 pi + 3 varpi) IS
    caught (positive control -> the negative is not vacuous).  This is the part of
    the paper's PSLQ-negative that is ROBUSTLY DECIDABLE at the reachable ~19
    digits; the wider weight-3 / disc-8 exclusion needs ~40 digits and is a
    driver-observed bound, not asserted here (see the module comment)."""
    with mp.workdps(25):          # constants accurate; W trusted to ~19 digits
        pi = mp.pi
        varpi = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(pi))
        basis = [mp.mpf(1), pi, varpi, pi ** 2, varpi ** 2, pi * varpi]
        W = _V_COLLINEAR * pi / 8
        decoy = mp.log(mp.mpf(11)) / mp.mpf('15.4')     # ~0.1557, structureless
        control = 2 * pi + 3 * varpi                     # planted height-3 low closure

        def is_low_height_closure(x):
            rel = mp.pslq([-x] + basis, tol=mp.mpf(10) ** -16,
                          maxcoeff=10 ** 8, maxsteps=10 ** 6)
            if rel is None or rel[0] == 0:
                return False                             # none / basis-internal
            return max(abs(c) for c in rel) <= _HEIGHT_LOW

        assert is_low_height_closure(control), "detector blind to a planted closure"
        assert not is_low_height_closure(W)              # robust disc-4 negative
        assert not is_low_height_closure(decoy)          # matched control
