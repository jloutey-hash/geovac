"""
Symbolic derivation of Paper 33 Theorem 2 (Dirac 8/8): Furry's theorem
from the Dirac spinor phase constraint.

Paper 33 claims V(a, a, q, m_q) = 0 identically for the diagonal Dirac
vertex element.  The mechanism (Theorem thm:dirac_8_of_8, Steps 1-4):

    psi_a = (f * Omega_kappa, i * g * Omega_{-kappa})^T
    alpha_z = ((0, sigma_z), (sigma_z, 0))
    psi_a^dagger alpha_z psi_a = -2 f g * Im[M],
        M = Omega_kappa^dagger sigma_z Omega_{-kappa}

and Im[M] == 0 pointwise because sigma_z conserves m_s (hence m_l), so
the e^{+-i m_l phi} phases in Y_{l}^{m_l *} Y_{l'}^{m_l} cancel exactly,
leaving M a real function of theta alone.  Since the integrand vanishes
POINTWISE, V(a, a, q, m_q) = 0 for EVERY photon mode (q, m_q).

This file ports debug/archive/qed_arc/paper33_furry_verification.py
(assert-free driver, 2026-04 QED arc) into real pytest asserts.  It is
deliberately INDEPENDENT of geovac.vector_qed.dirac_vertex_coupling:
that production function short-circuits the diagonal (`if a == b:
return 0.0`) as a documented kinematic identity, so testing it against
itself would be tautological.  Here the vertex bilinear is DERIVED from
sympy spin-spherical harmonics and shown to vanish symbolically (== 0,
exact), which is the derivation backing the production short-circuit.

Tier: SYMBOLIC PROOF (six concrete (kappa, m_j) cases, sympy-exact).
"""

import pytest
import sympy as sp
from sympy import I, Matrix, Rational, pi, simplify, sin
from sympy.physics.quantum.cg import CG


theta, phi = sp.symbols('theta phi', real=True)

# Representative diagonal cases (kappa, m_j): the full n_max=2 kappa set
# {-1, +1, -2} across distinct |m_j|, matching the archived driver.
DIAGONAL_CASES = [
    (-1, Rational(1, 2)),    # GS   kappa=-1, l=0, j=1/2, m_j=+1/2
    (-1, Rational(-1, 2)),   # GS   kappa=-1, l=0, j=1/2, m_j=-1/2
    (1, Rational(1, 2)),     # p1/2 kappa=+1, l=1, j=1/2, m_j=+1/2
    (1, Rational(-1, 2)),    # p1/2 kappa=+1, l=1, j=1/2, m_j=-1/2
    (-2, Rational(1, 2)),    # p3/2 kappa=-2, l=1, j=3/2, m_j=+1/2
    (-2, Rational(3, 2)),    # p3/2 kappa=-2, l=1, j=3/2, m_j=+3/2
]


def Y(l, m):
    """Spherical harmonic Y_l^m(theta, phi), expanded to functional form."""
    return sp.functions.special.spherical_harmonics.Ynm(
        l, m, theta, phi).expand(func=True)


def kappa_to_l(kappa):
    """Standard Dirac convention: kappa = -l-1 for j=l+1/2; kappa = +l for j=l-1/2.

    Local re-derivation on purpose -- do NOT import the production helper;
    this test must be independent of geovac.vector_qed / dirac_matrix_elements.
    """
    return -1 - kappa if kappa < 0 else kappa


def spin_spherical(l, j, m_j):
    """Spin-spherical harmonic Omega_{l,j,m_j} as a 2x1 sympy Matrix.

    Omega_{l,j,m_j} = sum_{m_s} <l, m_l; 1/2, m_s | j, m_j> Y_l^{m_l} chi_{m_s}
    with m_l = m_j - m_s.
    """
    half = Rational(1, 2)
    omega = Matrix([0, 0])
    for row, m_s in enumerate([half, -half]):
        m_l = m_j - m_s
        if abs(m_l) > l:
            continue
        cg = CG(l, m_l, half, m_s, j, m_j).doit()
        if cg == 0:
            continue
        omega[row] = omega[row] + cg * Y(l, m_l)
    return omega


def Omega_kappa(kappa, m_j):
    """Camporesi-Higuchi convention: Omega_{kappa, m_j} with j = |kappa| - 1/2."""
    l = kappa_to_l(kappa)
    j = abs(kappa) - Rational(1, 2)
    return spin_spherical(l, j, m_j)


def diagonal_bilinear(kappa, m_j):
    """The diagonal Dirac vertex bilinear V_LS + V_SL (pointwise on S^2).

    For a == b, m_j conservation forces m_q = 0, so only the spatial
    z-component of alpha contributes.  With psi = (Omega_kappa,
    i*Omega_{-kappa})^T (real radials f, g suppressed -- they multiply the
    whole bilinear):

        V_LS = Omega_kappa^dagger (i sigma_z Omega_{-kappa})
        V_SL = (-i Omega_{-kappa}^dagger)(sigma_z Omega_kappa)
    """
    sig_z = Matrix([[1, 0], [0, -1]])
    Om_kap = Omega_kappa(kappa, m_j)
    Om_neg = Omega_kappa(-kappa, m_j)  # same m_j by m-conservation
    LS = Om_kap.H * (I * sig_z * Om_neg)
    SL = (-I * Om_neg.H) * (sig_z * Om_kap)
    return (LS + SL)[0, 0]


def M_bilinear(kappa, m_j):
    """M = Omega_kappa^dagger sigma_z Omega_{-kappa} (Step 3 of the theorem)."""
    sig_z = Matrix([[1, 0], [0, -1]])
    return (Omega_kappa(kappa, m_j).H * sig_z * Omega_kappa(-kappa, m_j))[0, 0]


class TestFurryBilinearVanishesPointwise:
    """Theorem Step 2+3: psi^dagger alpha_z psi == 0 POINTWISE on S^2.

    This is the load-bearing identity: because the integrand vanishes
    before integration, V(a, a, q, m_q) = 0 for every q >= 1 and m_q.
    """

    @pytest.mark.parametrize("kappa,m_j", DIAGONAL_CASES)
    def test_bilinear_identically_zero(self, kappa, m_j):
        bilinear = diagonal_bilinear(kappa, m_j)
        assert simplify(bilinear) == 0, (
            f"V_LS + V_SL != 0 for kappa={kappa}, m_j={m_j}: {bilinear}")


class TestFurryMechanismIsMlConservation:
    """Paper 33 Remark rem:mechanism: the cancellation is Im M == 0 via
    m_l conservation, NOT mere Hermiticity of alpha_z."""

    @pytest.mark.parametrize("kappa,m_j", DIAGONAL_CASES)
    def test_im_M_vanishes(self, kappa, m_j):
        """Step 3: M = Omega_kappa^dag sigma_z Omega_{-kappa} is purely real."""
        M = M_bilinear(kappa, m_j)
        assert simplify(sp.im(M)) == 0, (
            f"Im M != 0 for kappa={kappa}, m_j={m_j}")

    def test_M_itself_not_identically_zero(self):
        """Non-vacuity: M != 0 generically, so the vanishing of the bilinear
        is the i-phase / Im-part cancellation, not a trivial zero."""
        nonzero_found = any(
            simplify(M_bilinear(kappa, m_j)) != 0
            for kappa, m_j in DIAGONAL_CASES
        )
        assert nonzero_found, (
            "M == 0 in ALL cases -- the theorem's mechanism (Im M = 0 with "
            "M != 0) would be vacuous; investigate")


class TestFurryVertexIntegralVanishes:
    """Theorem Step 4: the projected vertex integral
    V(a, a, q, 0) = int psi^dagger alpha_z psi * Y_q^0 dOmega == 0 (sympy-exact).

    Redundant given the pointwise zero, but this is the quantity the
    archived driver computed and the quantity Paper 33's proof closes on.
    """

    @pytest.mark.parametrize("kappa,m_j", DIAGONAL_CASES)
    def test_integral_zero_q1(self, kappa, m_j):
        bilinear = diagonal_bilinear(kappa, m_j)
        integrand = bilinear * Y(1, 0) * sin(theta)
        result = sp.integrate(sp.integrate(integrand, (phi, 0, 2 * pi)),
                              (theta, 0, pi))
        assert simplify(result) == 0

    def test_integral_zero_higher_q(self):
        """One higher-q spot check (q=3, GS case): pointwise zero implies
        zero against ANY photon harmonic."""
        bilinear = diagonal_bilinear(-1, Rational(1, 2))
        integrand = bilinear * Y(3, 0) * sin(theta)
        result = sp.integrate(sp.integrate(integrand, (phi, 0, 2 * pi)),
                              (theta, 0, pi))
        assert simplify(result) == 0


def test_production_shortcircuit_agrees_with_derivation():
    """Consistency (NOT the derivation): the production dirac_vertex_coupling
    diagonal short-circuit returns exactly 0.0, the value derived above."""
    from geovac.vector_qed import build_dirac_electron_states, dirac_vertex_coupling

    states = build_dirac_electron_states(2)
    for a in states:
        for q in (1, 2):
            for m_q in range(-q, q + 1):
                assert dirac_vertex_coupling(a, a, q, m_q) == 0.0
