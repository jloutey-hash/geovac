r"""Backing tests for Paper 12, Sec. "Restoring the Azimuthal Channels" [MEASURED]:
the general-m (mu > 0) Neumann V_ee via the MOMENT-RECURRENCE radial engine
(geovac.neumann_vee_general_m), which extends Paper 12's algebraic sigma-only
A_l/B_l/X_l tables to associated Legendre functions and -- the load-bearing win --
is STABLE at mu = 2 (the delta channel) where the differentiation engine
(geovac.prolate_general_m.vee_matrix, d^m Q_l on a grid) suffers the d^4 Q_l
endpoint cancellation and the H2 energy diverges.

Each test names the wrong answer it rejects (Sec. 9 guard rule).  Fire-tested by
`debug/firetest_p12_general_m.py`.

Scope this engine claims (and this file pins): the engine is now fully
QUADRATURE-FREE -- the low-l B_l seeds are CLOSED FORM (_seed_B_closed), the
only transcendental inputs being E_1(2c), Euler gamma and ln c (isolated in the
log-moment primitive _L_moments); the forward recurrence, not differentiation,
removes the mu = 2 instability.  The quadrature reference _seed_B is kept for the
seed-validation test only.
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac import neumann_vee_general_m as gm
from geovac import prolate_general_m as pg

E_EXACT = -1.174475          # Kolos & Wolniewicz
DE_EXACT = pg.DE_EXACT
R = pg.R_DEFAULT


def _de_pct(e: float) -> float:
    return 100.0 * (-1.0 - e) / DE_EXACT


def _energy(engine, j_max, l_max, mu_max, alpha=1.0, l_neumann=14):
    """Ground-state energy with V_ee from `engine` (gm or pg grid)."""
    basis = pg.generate_basis(j_max, l_max, mu_max, alpha)
    mom = pg.Moments(2.0 * alpha, 6 * max(j_max, l_max) + 6 * (mu_max + 2) + 20)
    s, h1 = pg.one_body(basis, R, 1.0, mom)
    if engine == "recur":
        v = gm.vee_matrix(basis, R, l_neumann)
    else:
        v = pg.vee_matrix(basis, R, pg.XiGrid(alpha), l_neumann)
    h = h1 + v + (1.0 / R) * s
    e, _kept, _tot = pg.solve_generalized(h, s)
    return e, len(basis)


# ======================================================================
# 1. The mu = 0 limit must reduce to the corpus's exact m = 0 machinery
# ======================================================================

def test_mu0_reduces_to_neumann_vee():
    """REJECTS: a general-m engine that silently breaks its own m = 0 base case.

    geovac.neumann_vee is the independently-derived, closed-form m = 0 path.  The
    moment-recurrence engine must agree with it elementwise at mu = 0, or the
    associated-Legendre generalization has changed the physics it extends.
    """
    from geovac.hylleraas import HylleraasBasisFunction
    from geovac.neumann_vee import compute_vee_matrix_neumann

    alpha = 1.0
    mine = pg.generate_basis(2, 2, 0, alpha)
    theirs = [HylleraasBasisFunction(b.j, b.k, b.l, b.m, 0, alpha) for b in mine]

    v_mine = gm.vee_matrix(mine, R, l_neumann=14)
    v_them = compute_vee_matrix_neumann(theirs, R, l_max=20)
    rel = np.abs(v_mine - v_them) / np.maximum(np.abs(v_them), 1e-12)
    assert rel.max() < 1e-7, (
        f"moment-recurrence V_ee disagrees with geovac.neumann_vee at m = 0 by "
        f"{rel.max():.2e} relative; the extension broke its own base case"
    )


# ======================================================================
# 2. The radial X table must match a high-precision reference
# ======================================================================

def test_X_blocks_match_high_precision_reference():
    """REJECTS: a silently-wrong recurrence (forward-B instability, a dropped
    weight, a sign error in the IBP correction) -- the failure mode a
    differentiation route would show as garbage at large l or m.

    (No longer @slow: the X-table's B_l seeds are closed-form as of v5.13.4, so
    this build is ~0.2s where the quadrature-seeded version was tens of seconds.)

    The references are computed independently of the engine's own B recurrence:
    the inner P integral in CLOSED FORM (monomial partial integrals) and the
    outer Q integral by direct mpmath quadrature of d^m Q_l (dps = 35).  They
    span a mu <= 1 block, a mu = 1 block, and three delta-channel (m = 4) blocks
    -- the regime the differentiation engine cannot reach.
    """
    # X[l,m,s](P1,P2) at alpha = 1  (per-electron rate c = 2)
    ref = {
        (2, 2, 2, 0, 0): 0.0100142205252448,
        (3, 1, 1, 0, 0): -0.002172642484118069,
        (4, 4, 4, 0, 0): 87.70947318323068,
        (6, 4, 4, 0, 0): 2144.07536568145,
        (8, 4, 4, 0, 2): 114762.1627845131,
    }
    ms_pairs = [(1, 1), (2, 2), (4, 4)]
    Xtab = gm.build_Xtab(ms_pairs, l_neumann=8, p_max=2, basis_alpha=1.0)

    worst = 0.0
    for (l, m, s, P1, P2), r in ref.items():
        got = Xtab[(l, m, s)][P1, P2]
        rel = abs(got - r) / abs(r)
        worst = max(worst, rel)
        assert rel < 1e-9, (
            f"X[{l},{m},{s}]({P1},{P2}) = {got:.12g} differs from the "
            f"high-precision reference {r:.12g} by {rel:.2e} relative"
        )
    assert worst < 1e-9
    # the ordered integral is symmetric in (P1, P2) by construction
    M = Xtab[(4, 4, 4)]
    assert np.allclose(M, M.T, atol=1e-12), "X_l^{m,s} must be symmetric in P1<->P2"


# ======================================================================
# 3. The delta channel (mu = 2) must be STABLE -- the load-bearing claim
# ======================================================================

@pytest.mark.slow
def test_mu2_delta_channel_is_stable():
    """REJECTS: the d^4 Q_l differentiation blow-up.

    The grid engine (differentiation) gives a NON-variational E ~ -14.5 Ha at this
    basis -- the delta contribution is unusable.  The moment-recurrence engine
    must instead give a variational energy (above the exact -1.174475) that adds
    a small delta gain on top of the |m| <= 1 value.

    The guard discriminates by contrast: it asserts the recurrence energy is
    variational AND that the grid engine at the SAME basis is not, so a test that
    trivially passed both ways is excluded.
    """
    e_pi, _ = _energy("recur", 2, 2, 1)      # |m| <= 1
    e_delta, n = _energy("recur", 2, 2, 2)   # |m| <= 2  (delta included)

    assert e_delta > E_EXACT, (
        f"mu <= 2 energy {e_delta:.6f} is BELOW the exact {E_EXACT:.6f}; the "
        f"delta channel is unstable (the differentiation blow-up was not cured)"
    )
    # the delta channel lowers the energy a little further, and stays sane
    assert e_delta <= e_pi + 1e-9, "adding the delta channel must not raise E"
    assert _de_pct(e_delta) > 98.5, (
        f"mu <= 2 recovers only {_de_pct(e_delta):.2f}% of D_e; the delta "
        f"channel build is not sound"
    )

    # contrast: the differentiation engine is NOT variational here (the failure
    # this engine exists to fix).  If the grid engine ever became stable at this
    # basis, this contrast would no longer discriminate and the test says so.
    e_grid, _ = _energy("grid", 2, 2, 2)
    assert e_grid < E_EXACT - 1.0, (
        f"the differentiation (grid) engine gives E = {e_grid:.4f} at mu = 2, "
        f"which is no longer the documented catastrophic blow-up; this contrast "
        f"can no longer prove the recurrence engine fixed anything"
    )


# ======================================================================
# 4. At mu <= 1 the recurrence engine must reproduce the sound grid engine
# ======================================================================

def test_reproduces_grid_engine_at_mu1_where_grid_is_sound():
    """REJECTS: a recurrence engine that drifts from the validated differentiation
    result in the regime where the grid engine is trustworthy.

    At (2,2), |m| <= 1 the Neumann truncation is capped low enough that the grid
    engine's q_deriv has not yet lost precision, so the two engines must agree to
    the conditioning floor of this basis.
    """
    e_recur, n1 = _energy("recur", 2, 2, 1)
    e_grid, n2 = _energy("grid", 2, 2, 1)
    assert n1 == n2 == 54
    assert abs(e_recur - e_grid) < 1e-5, (
        f"recurrence E = {e_recur:.6f} vs grid E = {e_grid:.6f} differ by "
        f"{1e6 * abs(e_recur - e_grid):.1f} uHa at mu <= 1, where the grid "
        f"engine is sound"
    )
    assert _de_pct(e_recur) > 98.5


# ======================================================================
# 5. The engine must not be secretly differentiating Q at high order
# ======================================================================

def test_engine_carries_intact_weight_not_expanded():
    """REJECTS: expanding (xi^2-1)^s into xi-powers (memo finding 2).

    The bare xi^p d^m Q_l moment diverges for m >= 1; only the intact (xi^2-1)^s
    weight (s >= m/2) regularises the xi = 1 endpoint.  The seed B_l^{m=4,s=4}(0)
    must therefore be finite and positive -- a naive monomial expansion would
    return inf/nan.  (Both the quadrature reference _seed_B and the production
    closed form must give the same finite value.)
    """
    import mpmath as mp
    with mp.workdps(30):
        b = gm._seed_B(4, 4, 4, 0, mp.mpf(2.0))
        assert mp.isfinite(b) and b > 0, (
            f"regularised delta-channel seed B_4^(4,4)(0) = {b}; the intact "
            f"weight is not being carried (finding 2)"
        )


# ======================================================================
# 6. The closed-form B seeds (quadrature-free) must match the quadrature ref
# ======================================================================

def test_L_moment_primitive_closed_form():
    """REJECTS: a wrong closed form for the one transcendental primitive
    L_n(c) = int_1^inf xi^n Q_0(xi) e^{-c xi} dxi.

    A sign error, a dropped E_1(2c) (the ln(xi+1) tail) or a missing Euler-gamma
    (the ln(xi-1) branch) would fail against direct mpmath quadrature.
    """
    import mpmath as mp
    with mp.workdps(30):
        for c in (mp.mpf(2), mp.mpf('3.7')):
            L = gm._L_moments(6, c)
            for n in range(7):
                q = mp.quad(lambda u, n=n, c=c:
                            (1 + u) ** n * mp.mpf('0.5') * mp.log((2 + u) / u)
                            * mp.e ** (-c * (1 + u)),
                            [0, mp.mpf('0.1'), mp.mpf(1), mp.mpf(4), mp.inf])
                rel = abs(L[n] - q) / abs(q)
                assert rel < 1e-25, (
                    f"L_{n}({float(c)}) closed form {L[n]} differs from quadrature "
                    f"{q} by {float(rel):.1e}; the E_1/gamma/ln primitive is wrong"
                )


@pytest.mark.slow
def test_closed_form_B_seeds_match_quadrature():
    """REJECTS: a closed-form B seed that disagrees with the mpmath-quadrature
    reference _seed_B.

    (@slow: exercises the quadrature reference _seed_B ~70 times, ~10s -- the
    closed form it validates is itself sub-millisecond.)

    This closed form REPLACED the quadrature that was 95% of vee_mp.  It rests on
    s >= m fully polynomializing every (xi-1)^{-k} pole term of d^mQ_l; a bug in
    that polynomialization (wrong coeff_k, a (xi+1) vs (xi-1) swap, an off-by-one
    in s-k) would show here.  Covers sigma/pi/delta and mixed (m,s) at the two
    seed orders l = m, m+1.
    """
    import mpmath as mp
    with mp.workdps(30):
        c = mp.mpf(2.0)
        Lmom = gm._L_moments(40, c)
        A = gm._mono_moments(c, 60)
        worst = mp.mpf(0)
        for (m, s) in [(0, 0), (1, 1), (2, 2), (1, 2), (2, 3)]:
            for l in (m, m + 1):
                for p in range(7):
                    bc = gm._seed_B_closed(l, m, s, p, Lmom, A)
                    bq = gm._seed_B(l, m, s, p, c)
                    if abs(bq) > mp.mpf('1e-30'):
                        worst = max(worst, abs(bc - bq) / abs(bq))
        assert worst < 1e-20, (
            f"closed-form B seed disagrees with quadrature by {float(worst):.1e} "
            f"relative; the polynomialization (s>=m) has a bug"
        )


def test_B_table_uses_closed_form_not_quadrature():
    """REJECTS: a silent regression to quadrature seeding.

    The production _B_table must not call the quadrature _seed_B (its cost was the
    V_ee bottleneck).  Patch _seed_B to raise; _B_table must still succeed and its
    l = m, m+1 rows must match the closed-form seeds.
    """
    import mpmath as mp
    with mp.workdps(30):
        c = mp.mpf(2.0)
        orig = gm._seed_B
        gm._seed_B = lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("_B_table must not use quadrature _seed_B"))
        try:
            B = gm._B_table(1, 1, 6, 5, c)
        finally:
            gm._seed_B = orig
        Lmom = gm._L_moments(30, c)
        A = gm._mono_moments(c, 40)
        for p in range(4):
            ref = gm._seed_B_closed(1, 1, 1, p, Lmom, A)
            assert abs(B[(1, p)] - ref) < 1e-25
