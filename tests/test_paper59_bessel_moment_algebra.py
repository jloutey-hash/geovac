"""Paper 59 sec:bessel_algebra -- the Bessel-moment period algebra of the
integrated three-centre observable: master family, Wronskian determinant, the
Broadhurst-Mellit / Fresan-Sabbah-Yu quadratic period relations, and the
critical-L-value negative.

Self-contained (no debug/ import, per the transient-dir policy).  Ports the
verified master-family construction from debug/routeC_bessel_moment_algebra.py.

L4 Picard-Fuchs operator (Paper 59 eq:pf), rho fixed, ' = d/dD:
    D rho L'''' + 2 rho L''' + D(1-2rho) L'' + (1-2rho) L' - D(1-rho) L = 0.
The four masters are thimble Laplace integrals s_c(D) = int_{gamma_c} e^{-Dx}/sqrt(Q),
Q = (x^2-1)(rho x^2 + 1 - rho); their k-th D-derivative brings down (-x)^k.
"""
from __future__ import annotations

import mpmath as mp
import pytest


# --- Wronskian: self-adjoint sub-leading ratio p3/p4 = 2/D  =>  W ~ D^-2 ---

def test_paper59_wronskian_is_D_minus_2():
    """[SYMBOLIC] Abel's identity on the (self-adjoint) L4: the sub-leading ratio
    p3/p4 = (2 rho)/(D rho) = 2/D forces the Wronskian of the master matrix to the
    pure monomial W(D) = W0 * D^-2 (the Broadhurst-Mellit / Zhou determinant), the
    elliptic lift of W[K0,I0] = 1/D.  No PSLQ; precision-independent."""
    import sympy as sp
    D, rho = sp.symbols('D rho', positive=True)
    p4, p3 = D * rho, 2 * rho
    ratio = sp.simplify(p3 / p4)
    assert sp.simplify(ratio - 2 / D) == 0
    W = sp.simplify(sp.exp(-sp.integrate(ratio, D)))     # Abel: W'/W = -(p3/p4)
    assert sp.simplify(W * D ** 2) == 1                  # W = W0 * D^-2


# --- dim S_k(Gamma(2)): first cusp form is weight 6 (the L-value negative) ---

def test_paper59_gamma2_first_cusp_form_is_weight_6():
    """[OBSERVATION] M_*(Gamma(2)) = C[theta2^4, theta4^4] is free on two weight-2
    generators (3 cusps => dim Eis = 3 for k>=4), so dim S_k = (k/2 + 1) - dim Eis
    vanishes for k = 2, 4 and first reaches 1 at k = 6.  The observable's weight is
    <= 3, so it cannot be a critical L-value of a weight-6 cusp form (motivic weight
    5 >> 3); with the Eisenstein-flavoured pairing it is an Eisenstein/CM period,
    not a cusp-form L-value."""
    def dim_S(k):
        dim_M = k // 2 + 1
        dim_Eis = 3 if k >= 4 else 2
        return dim_M - dim_Eis
    assert dim_S(2) == 0
    assert dim_S(4) == 0
    assert dim_S(6) == 1          # first cusp form, weight 6 = theta2^4 theta3^4 theta4^4
    assert dim_S(8) == 2


# --- master family + the {-pi, 0, 2pi} Broadhurst-Mellit / FSY period pairing ---

def _sK(D, rho, k):     # [1,inf): K0-sector, the physical N ~ e^{-D}
    def f(x):
        Q = (x * x - 1) * (rho * x * x + 1 - rho)
        return mp.mpf(0) if Q <= 0 else (-x) ** k * mp.e ** (-D * x) / mp.sqrt(Q)
    return mp.quad(f, [1, mp.mpf('1.05'), mp.mpf('1.3'), 2, 4, 8, 16, mp.inf])


def _sI(D, rho, k):     # (-1,1): I0-sector ~ e^{+D}
    def f(x):
        nQ = (1 - x * x) * (rho * x * x + 1 - rho)
        return mp.mpf(0) if nQ <= 0 else (-x) ** k * mp.e ** (-D * x) / mp.sqrt(nQ)
    return mp.quad(f, [-1, mp.mpf('-0.5'), 0, mp.mpf('0.5'), 1])


def _sJ(D, rho, k):     # imaginary axis x=iu: J0-sector ~ cos(wD)
    w = mp.sqrt((1 - rho) / rho)

    def f(u):
        aQ = (u * u + 1) * ((1 - rho) - rho * u * u)
        if aQ <= 0:
            return mp.mpf(0)
        return ((-1j * u) ** k * mp.e ** (-1j * D * u) / mp.sqrt(aQ)).real
    return mp.quad(f, [-w, 0, w])


@pytest.mark.slow
def test_paper59_L4_annihilates_the_physical_master_N():
    """[MEASURED, eq:pf residual] Direct annihilation witness: plug the PHYSICAL
    master period N(D)=s_K(D) (the cut [1,inf) thimble, ~e^{-D}, which is N(D)
    itself) and its D-derivatives straight into the rank-4 Picard-Fuchs operator
    eq:pf and check L4[N] ~ 0.  The other L4 tests certify the OPERATOR
    (self-adjointness / exponents / factorization battery / integer monodromy)
    and rule out the WRONG candidates (single Bessels / periods are NOT
    D-solutions); this is the missing test that the actual master period IS a
    solution.  The k-th D-derivative of s_K is _sK(.,.,k) (it brings down
    (-x)^k).  Measured residual ~1e-18 at dps=30 (adaptive quad; the ~1e-24
    figure of eq:pf is the symbolic-derivation residual)."""
    with mp.workdps(30):
        worst = mp.mpf(0)
        for rho in ['0.5', '0.25', '0.2']:
            r = mp.mpf(rho)
            for D in [mp.mpf('0.7'), mp.mpf(1), mp.mpf('1.6')]:
                N = [_sK(D, r, k) for k in range(5)]   # N, N', N'', N''', N''''
                res = (D * r * N[4] + 2 * r * N[3] + D * (1 - 2 * r) * N[2]
                       + (1 - 2 * r) * N[1] - D * (1 - r) * N[0])
                assert abs(res) < mp.mpf(10) ** -14, (rho, D, res)
                worst = max(worst, abs(res))
        assert worst < mp.mpf(10) ** -14, worst


def _concomitant(y, z, D, rho):
    # Lagrange bilinear concomitant of the self-adjoint L4 (a2=D*rho, a1=D*(1-2rho))
    Ly, Ly1, Ly2, Ly3 = y
    Lz, Lz1, Lz2, Lz3 = z
    a2, a1 = D * rho, D * (1 - 2 * rho)
    return (Lz * rho * (Ly2 + D * Ly3) - Lz1 * (a2 * Ly2)
            - Ly * rho * (Lz2 + D * Lz3) + Ly1 * (a2 * Lz2)
            + a1 * (Lz * Ly1 - Ly * Lz1))


@pytest.mark.slow
def test_paper59_period_pairing_minus_pi_0_2pi():
    """[MEASURED, 25 digits] The Lagrange bilinear concomitant of the (self-adjoint)
    L4 on the master periods is the Broadhurst-Mellit / FSY quadratic relation, D- and
    rho-independent: B[K,I] = -pi, B[K,J] = 0, B[I,J] = 2pi.  Checked at rho=1/2
    (all three entries) and rho=1/4 (rho-independence of B[K,I]).

    dps=55 so the test witnesses the paper's (sec:bessel_algebra) quoted 25-digit
    precision (worst residual ~2.3e-29 at dps=55; the earlier dps=30 floored the
    assertion at 1e-10 while the relations agree to ~1e-17 there and 25+ digits here)."""
    with mp.workdps(55):
        pi = mp.pi
        D = mp.mpf(1)
        for rho in ['0.5', '0.25']:
            r = mp.mpf(rho)
            YK = [_sK(D, r, k) for k in range(4)]
            YI = [_sI(D, r, k) for k in range(4)]
            BKI = _concomitant(YK, YI, D, r)
            assert abs(BKI + pi) < mp.mpf(10) ** -25, (rho, BKI)   # = -pi, rho-independent
            if rho == '0.5':
                YJ = [_sJ(D, r, k) for k in range(4)]
                BKJ = _concomitant(YK, YJ, D, r)
                BIJ = _concomitant(YI, YJ, D, r)
                assert abs(BKJ) < mp.mpf(10) ** -25, BKJ            # = 0
                assert abs(BIJ - 2 * pi) < mp.mpf(10) ** -25, BIJ   # = 2 pi


@pytest.mark.slow
def test_paper59_all_three_regular_masters_satisfy_pf():
    """[MEASURED] All THREE regular thimble masters -- s_K (the physical N(D)),
    s_I, and s_J -- satisfy the rank-4 Picard-Fuchs operator eq:pf (Paper 59
    sec:bessel_algebra: "the three regular masters satisfy Eq.(pf)").  The
    companion test_paper59_L4_annihilates_the_physical_master_N pins only s_K;
    this witnesses the *plural* claim.  Worst residual across the (rho,D) grid:
    s_K ~1e-18, s_I ~3e-17, s_J ~1e-15 (s_J is the worst regular master, so the
    paper's global bound is ~1e-15; guarded here at <1e-14)."""
    with mp.workdps(30):
        def L4res(master, D, r):
            N = [master(D, r, k) for k in range(5)]
            return abs(D * r * N[4] + 2 * r * N[3] + D * (1 - 2 * r) * N[2]
                       + (1 - 2 * r) * N[1] - D * (1 - r) * N[0])
        for master in (_sK, _sI, _sJ):
            for rho in ['0.5', '0.25', '0.2']:
                r = mp.mpf(rho)
                for D in [mp.mpf('0.7'), mp.mpf(1), mp.mpf('1.6')]:
                    assert L4res(master, D, r) < mp.mpf(10) ** -14, (master.__name__, rho, D)
