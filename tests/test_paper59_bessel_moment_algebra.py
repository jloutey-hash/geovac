"""Paper 61 (companion of Paper 59) sec:bessel_algebra -- the Bessel-moment period algebra of the
integrated three-centre observable: master family, Wronskian determinant, the
Broadhurst-Roberts / Fresan-Sabbah-Yu quadratic period relations, and the
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


# --- master family + the {-pi, 0, 2pi} Broadhurst-Roberts / FSY period pairing ---

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
    L4 on the master periods is the Broadhurst-Roberts / FSY quadratic relation, D- and
    rho-independent: B[K,I] = -pi, B[K,J] = 0, B[I,J] = 2pi.  Checked at rho=1/2
    (all three entries) and rho=1/4 (rho-independence of B[K,I]).

    dps=55 so the test witnesses the companion Paper 61's (sec:bessel_algebra) quoted 25-digit
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
    s_I, and s_J -- satisfy the rank-4 Picard-Fuchs operator eq:pf (Paper 61
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


def test_paper61_wronskian_constant_W0_is_pi_squared_over_rho_squared():
    """[SYMBOLIC] Paper 61 eq:W0 -- the Wronskian CONSTANT is W_0 = pi^2/rho^2.

    Abel's identity fixes only the D-dependence (W(D) = W_0 D^-2, the companion
    test above).  W_0 itself is fixed by the branch-point data of
    Q(x) = (x^2-1)(rho x^2 + 1 - rho): with branch points {+-1, +-i w},
    w^2 = (1-rho)/rho,

        sum_c x_c = 0                      (exponentials cancel)
        prod_c Q'(x_c) = -16 w^2           => prod_c A_c = pi^2 / (4 i w)
        Vandermonde prod_{a<b}(x_b - x_a) = 4 i w / rho^2
        =>  W_0 = pi^2 / rho^2

    up to the thimble ordering/orientation sign -- the branch of
    sqrt(1-rho)/sqrt(rho-1) = -i on 0 < rho < 1.

    This RETIRES the claim that "the transcendence of the individual masters
    cancels in their determinant"
    [retracted 2026-09-07: p61-w0-transcendence-cancels].
    It does not cancel:  pi^2 is SQUARED, and W_0 depends on rho.

    Deliberately checks rho-DEPENDENCE and the EXPONENT, not just the value at
    one rho:  the retired reading is exactly a constant W_0, and the existing
    Wronskian test normalises W_0 to 1 by construction, so a single-point value
    check would be satisfied by the very belief this test exists to exclude."""
    import sympy as sp

    x, rho = sp.symbols('x rho', positive=True)
    w = sp.sqrt((1 - rho) / rho)
    Q = (x ** 2 - 1) * (rho * x ** 2 + 1 - rho)
    Qp = sp.diff(Q, x)
    xs = [sp.Integer(1), sp.Integer(-1), sp.I * w, -sp.I * w]

    # the four branch points really are the roots of Q
    for c in xs:
        assert sp.simplify(Q.subs(x, c)) == 0

    # the three ingredients, each independently
    assert sp.simplify(sum(xs)) == 0
    assert sp.simplify(sp.prod([Qp.subs(x, c) for c in xs]) + 16 * w ** 2) == 0
    vdm = sp.prod([xs[b] - xs[a] for a in range(4) for b in range(a + 1, 4)])
    assert sp.simplify(vdm - 4 * sp.I * w / rho ** 2) == 0

    W0 = sp.simplify(sp.pi ** 2 / sp.sqrt(sp.prod([Qp.subs(x, c) for c in xs])) * vdm)

    sym, vals = {}, {}
    for r in (sp.Rational(3, 10), sp.Rational(1, 2), sp.Rational(71, 100)):
        got = sp.nsimplify(sp.simplify(W0.subs(rho, r)))
        want = sp.pi ** 2 / r ** 2
        assert sp.simplify(got - want) == 0, (r, got, want)
        sym[r], vals[r] = got, sp.N(got, 30)

    # (b) rho-DEPENDENCE: rejects W_0 = 1 and every other constant outright
    r1, r2 = sp.Rational(3, 10), sp.Rational(71, 100)
    ratio = vals[r1] / vals[r2]
    assert abs(ratio - 1) > sp.Rational(1, 2), "W_0 must not be rho-independent"

    # (c) the EXPONENT is -2, not -1 or -3: log-ratio of the two samples
    expo = sp.log(ratio) / sp.log(r2 / r1)
    assert abs(sp.N(expo) - 2) < 1e-25, ("rho-exponent", sp.N(expo))

    # (d) pi^2 is present SYMBOLICALLY, i.e. the transcendence did NOT cancel.
    # Checked on the symbolic value, not the Float: W_0/pi^2 must be the exact
    # rational 1/rho^2, which is what "does not cancel" means here.
    half = sp.Rational(1, 2)
    assert sp.simplify(sym[half] / sp.pi ** 2 - 4) == 0
    assert sp.simplify(sym[half] / sp.pi ** 2).is_rational is True
