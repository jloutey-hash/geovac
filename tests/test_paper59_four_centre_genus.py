"""Paper 59 -- the genus ladder stops at one for the four-centre ERI (XY|ZW).

Backs the pre-registered prediction P1 of
docs/predictions/2026-10-13_four_centre_genus.md (registered 2026-10-13
before any four-centre computation; outcome recorded in CHANGELOG v5.20.14).
Self-contained (no debug/ import, per the transient-dir policy); the formula
is re-implemented from debug/four_centre_genus_check.py, whose recorded
values are the pins below.

Claim.  The general four-centre 1s Slater ERI (all exponents zeta = 1)
reduces to a two-parameter Feynman integral over ONE radial k-integral,

    (AB|CD) = (8/pi) int_0^1 ds int_0^1 dt int_0^inf dk
              j0(k |W(s,t)|) G(k,s;D1) G(k,t;D2),

    G(k,t;D) = c e^{-D Delta} (D^2/Delta^3 + 3D/Delta^4 + 3/Delta^5),
    c = t(1-t),  Delta = sqrt(c k^2 + 1),
    D1 = |A-B|, D2 = |C-D|,  W(s,t) = P1(s) - P2(t),  P(t) = tB + (1-t)A,

whose ONLY square roots are Delta1(s) and Delta2(t): the fibre at fixed
(s,t) is the Paper 59 genus-1 curve y^2 = (c1 k^2+1)(c2 k^2+1).  No genus 2
at four centres.  Geometry enters only through |W| (phase) and D1, D2.

Recorded values (debug/data/four_centre_genus_check.json, 2026-10-13):
    three-centre (XY|XZ), Paper 59's own geometry:  reduced 0.204941652531
                                                    MD n=16  0.204941652260
    four-centre  (XY|ZW), four distinct centres:    reduced 0.218814149639
                                                    MD n=16  0.218814149425
The MD (McMurchie-Davidson on n-Gaussian STO fits) column converges to the
Slater value from below as n grows; the reduced column is the true Slater
integral and is n_st-converged to 1e-12 by n_st = 24.
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad

# Paper 59's three-centre geometry (test_routeC_momentum.py) and the
# four-centre geometry of the registered check.
X = [0.0, 0.0, 0.0]
Y = [0.0, 0.0, 2.0]
Z = [1.5, 0.0, 0.5]
THREE_CENTRE = (X, Y, X, Z)
FOUR_CENTRE = ([0.0, 0.0, 0.0], [0.0, 0.0, 1.4], [1.9, 0.4, -0.2], [2.2, 1.1, 1.0])

THREE_CENTRE_REDUCED = 0.204941652531   # Paper 59's own (XY|XZ) Slater value
FOUR_CENTRE_REDUCED = 0.218814149639    # recorded, reduced two-root formula
FOUR_CENTRE_MD16 = 0.218814149425       # recorded, MD on 16-Gaussian fits


def _G(k: float, t: float, D: float) -> float:
    c = t * (1.0 - t)
    De = np.sqrt(c * k * k + 1.0)
    return c * np.exp(-D * De) * (D * D / De**3 + 3.0 * D / De**4 + 3.0 / De**5)


def _j0(x: float) -> float:
    return 1.0 if x < 1e-12 else np.sin(x) / x


def reduced_eri(A, B, C, Dd, n_st: int = 24) -> float:
    """(AB|CD) via the reduced two-square-root formula; sin^2 maps on s, t."""
    A, B, C, Dd = map(np.asarray, (A, B, C, Dd))
    D1 = float(np.linalg.norm(A - B))
    D2 = float(np.linalg.norm(C - Dd))
    x, wx = leggauss(n_st)
    th = 0.25 * np.pi * (x + 1.0)                 # theta in (0, pi/2)
    u = np.sin(th) ** 2                           # u in (0, 1)
    wu = wx * 0.25 * np.pi * 2.0 * np.sin(th) * np.cos(th)
    total = 0.0
    for si, ws in zip(u, wu):
        P1 = si * B + (1.0 - si) * A
        for ti, wt in zip(u, wu):
            P2 = ti * Dd + (1.0 - ti) * C
            Wn = float(np.linalg.norm(P1 - P2))
            val, _ = quad(lambda k: _j0(k * Wn) * _G(k, si, D1) * _G(k, ti, D2),
                          0.0, np.inf, limit=400, epsabs=1e-13, epsrel=1e-11)
            total += ws * wt * val
    return 8.0 / np.pi * total


def _md_eri(A, B, C, Dd, n_gauss: int) -> float:
    """Independent route: McMurchie-Davidson on an n-Gaussian fit of the 1s STO."""
    from geovac import noci_engine as GE
    sh_a, sh_d, _ = GE.fit_sto_shape(0, 1, n_gauss=n_gauss)
    shapes = {"1s": (sh_a, sh_d)}
    f = [GE.sto_shape_basis(np.asarray(P, float), "1s", 1.0, shapes, (0, 0, 0))
         for P in (A, B, C, Dd)]
    return float(GE.eri_md(*f))


# ----------------------------------------------------------------- (a) values
def test_reduced_formula_reproduces_paper59_three_centre_value():
    """The two-root formula, evaluated at Paper 59's (XY|XZ) geometry, gives
    Paper 59's own Slater value 0.204941652531 (recomputed here, n_st = 24).

    Rejects: a formula that is NOT the (XY|XZ) integral in the shared-centre
    limit (e.g. the kernel with the 3D/Delta^4 term dropped gives 0.2049 off
    at the 1e-4 level; a wrong prefactor 8/pi -> 4/pi halves it).
    Tolerance 2e-11 (recorded n_st = 24 -> 56 drift 1.4e-14)."""
    val = reduced_eri(*THREE_CENTRE, n_st=24)
    assert abs(val - THREE_CENTRE_REDUCED) < 2e-11, (val, THREE_CENTRE_REDUCED)


def test_reduced_formula_reproduces_recorded_four_centre_value():
    """The same two-root formula at four DISTINCT centres reproduces the
    recorded value 0.218814149639 (recomputed, n_st = 24).

    Rejects: any drift of the reduced formula from the registered outcome
    (falsifier 2 of the registration, half of it: the formula itself has not
    moved).  Tolerance 2e-11 (recorded n_st = 24 -> 56 drift 4e-15)."""
    val = reduced_eri(*FOUR_CENTRE, n_st=24)
    assert abs(val - FOUR_CENTRE_REDUCED) < 2e-11, (val, FOUR_CENTRE_REDUCED)


def test_reduced_formula_matches_independent_md_route_fast():
    """Falsifier 2 (fast form): the reduced two-root formula agrees with an
    INDEPENDENT evaluation of the same four-centre integral -- Gaussian-fit
    McMurchie-Davidson -- here at n = 8 Gaussians, which the registered run
    places 6.5e-8 below the Slater value (fit error, converging with n).

    Rejects: a two-root formula that is exact-looking but computes a different
    integral (e.g. the four-centre (XY|ZW) with W set to zero gives 0.2x, off
    by 1e-2; the shared-centre value 0.2049 is off by 1.4e-2).
    Tolerance 3e-7 (recorded |reduced - MD8| = 6.5e-8)."""
    red = reduced_eri(*FOUR_CENTRE, n_st=16)
    md8 = _md_eri(*FOUR_CENTRE, n_gauss=8)
    assert abs(red - md8) < 3e-7, (red, md8, abs(red - md8))


@pytest.mark.slow
def test_reduced_formula_matches_independent_md_route_n16():
    """Falsifier 2 (full form, ~12 s): reduced two-root formula vs MD on
    16-Gaussian fits, the registered comparison.  Recorded: 0.218814149639
    vs 0.218814149425, |diff| 2.1e-10 (the residual Gaussian-fit error).

    Rejects: any discrepancy above the fit floor, i.e. a missing third root
    or a wrong kernel would show at >= 1e-8.  Tolerance 1e-9."""
    red = reduced_eri(*FOUR_CENTRE, n_st=24)
    md16 = _md_eri(*FOUR_CENTRE, n_gauss=16)
    assert abs(md16 - FOUR_CENTRE_MD16) < 1e-11, (md16, FOUR_CENTRE_MD16)
    assert abs(red - md16) < 1e-9, (red, md16, abs(red - md16))


# -------------------------------------------------------------- (b) structure
def test_kernel_is_second_z_derivative_of_yukawa_product():
    """Independent derivation of G.  The two-centre 1s Slater density
    e^{-zA rA} e^{-zB rB} is d^2/dzA dzB of the Yukawa product; the Feynman
    parametrisation of the Yukawa-product transform is e^{-D Delta}/Delta
    with Delta^2 = t(1-t) k^2 + t zB^2 + (1-t) zA^2.  So
        G(k,t;D) = d^2/dzA dzB [ e^{-D Delta}/Delta ] |_{zA=zB=1};
    the prefactor c = t(1-t) in G is produced by the two dDelta/dz factors,
    not by an external measure.  sympy must reduce the difference against the
    explicit kernel to 0.

    Rejects: a kernel with a wrong coefficient (3/Delta^5 -> 2/Delta^5 or
    3D/Delta^4 -> 2D/Delta^4), a wrong power, or a spurious extra factor of
    c (residual c(c-1)E, caught when this test was first written) -- the
    residual is then a nonzero rational function of Delta."""
    import sympy as sp
    k, t, D, zA, zB = sp.symbols("k t D z_A z_B", positive=True)
    c = t * (1 - t)
    Delta = sp.sqrt(c * k**2 + t * zB**2 + (1 - t) * zA**2)
    yuk = sp.exp(-D * Delta) / Delta
    derived = sp.diff(yuk, zA, zB).subs({zA: 1, zB: 1})
    De = sp.sqrt(c * k**2 + 1)
    explicit = c * sp.exp(-D * De) * (D**2 / De**3 + 3 * D / De**4 + 3 / De**5)
    resid = sp.simplify(derived - explicit)
    assert resid == 0, resid


def test_integrand_square_roots_are_exactly_the_two_paper59_quadrics():
    """Falsifiers 1 and 3: the k-integrand at fixed (s,t) contains exactly two
    square-root atoms, with radicands c1 k^2+1 and c2 k^2+1 (c1 = s(1-s),
    c2 = t(1-t)), and NOTHING else -- no third root, no root coupling s and t,
    no root in |W|.  With Delta1, Delta2 as symbols and the exponential /
    spherical-Bessel phase factored off, what remains is a rational function
    of (k, Delta1, Delta2).  That is the genus-1 fibre y^2 = (c1 k^2+1)(c2 k^2+1).

    Rejects: a reduction needing a third independent root (e.g. a factor
    1/sqrt(c1 c2 k^4 + (c1+c2) k^2 + 1) written as one radical, or
    sqrt(s t k^2 + 1)) -- the radicand set then differs from the two quadrics;
    and a non-rational remainder (a log or E1 in Delta) is rejected by the
    rational-function check."""
    import sympy as sp
    k, s, t, D1, D2, W = sp.symbols("k s t D_1 D_2 W", positive=True)
    c1, c2 = s * (1 - s), t * (1 - t)
    De1, De2 = sp.sqrt(c1 * k**2 + 1), sp.sqrt(c2 * k**2 + 1)

    def G(De, c, D):
        return c * sp.exp(-D * De) * (D**2 / De**3 + 3 * D / De**4 + 3 / De**5)

    integrand = sp.sin(k * W) / (k * W) * G(De1, c1, D1) * G(De2, c2, D2)
    radicands = {sp.expand(p.base) for p in integrand.atoms(sp.Pow)
                 if p.exp.is_Rational and p.exp.q == 2}
    assert radicands == {sp.expand(c1 * k**2 + 1), sp.expand(c2 * k**2 + 1)}, radicands

    # strip the two roots to symbols and the phase factors; remainder rational
    d1, d2 = sp.symbols("Delta_1 Delta_2", positive=True)
    core = integrand.subs({De1: d1, De2: d2})
    core = sp.simplify(core / (sp.sin(k * W) / (k * W) * sp.exp(-D1 * d1 - D2 * d2)))
    assert core.free_symbols <= {k, s, t, D1, D2, d1, d2}
    assert not core.has(sp.exp, sp.sin, sp.log, sp.sqrt)
    assert core.is_rational_function(k, d1, d2), core

    # the fingerprint fibre (D1 = D2 = 0, W = 0): 9 c1 c2 / (Delta1^5 Delta2^5)
    fibre = sp.simplify(core.subs({D1: 0, D2: 0}))
    assert sp.simplify(fibre - 9 * c1 * c2 / (d1**5 * d2**5)) == 0, fibre


def test_geometry_enters_only_through_W_D1_D2():
    """P2 of the registration: positions enter the reduced integrand only via
    |W(s,t)|, D1, D2.  Two four-centre geometries related by a rigid motion
    (rotation + translation) have identical reduced values; a geometry with
    the same D1, D2 but different relative placement does NOT (so |W| is
    load-bearing, and the agreement is not a trivial constant).

    Rejects: a formula that depends on absolute coordinates (fails the rigid
    motion check) or ignores W (passes the second geometry wrongly)."""
    A, B, C, Dd = (np.asarray(p, float) for p in FOUR_CENTRE)
    th = 0.7
    R = np.array([[np.cos(th), -np.sin(th), 0.0],
                  [np.sin(th), np.cos(th), 0.0],
                  [0.0, 0.0, 1.0]])
    shift = np.array([0.3, -1.1, 2.0])
    moved = [R @ p + shift for p in (A, B, C, Dd)]
    v0 = reduced_eri(A, B, C, Dd, n_st=16)
    v1 = reduced_eri(*moved, n_st=16)
    assert abs(v0 - v1) < 1e-11, (v0, v1)
    # same D1, D2, different W: translate the (C, D) pair away by 0.5 bohr
    far = (A, B, C + np.array([0.5, 0.0, 0.0]), Dd + np.array([0.5, 0.0, 0.0]))
    v2 = reduced_eri(*far, n_st=16)
    assert abs(v0 - v2) > 1e-3, (v0, v2)
