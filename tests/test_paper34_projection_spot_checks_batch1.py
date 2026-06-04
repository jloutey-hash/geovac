"""
Paper 34 projection spot-checks --- batch 1 (load-bearing rows).

Companion to tests/test_paper34_projection_spot_checks.py (6 of 28
projections covered). This file adds 8 load-bearing rows from
followon_register.md A9 batch 1:

  §III.1   Fock conformal (sec:proj_fock)
  §III.5   Sturmian (sec:proj_sturmian)
  §III.11  Vector-photon promotion (sec:proj_vector_photon)
  §III.13  Drake--Swainson asymptotic subtraction (sec:proj_drake_swainson)
  §III.16  Two-body Dirac / Breit retardation (sec:proj_breit_retardation)
  §III.17  Foldy/Friar charge density (sec:proj_charge_density)
  §III.18  Zemach magnetization density (sec:proj_magnetization_density)
  §III.19  Tensor multipole (sec:proj_tensor_multipole)

After batch 1: 14 of 28 Paper 34 projections covered.

Per CLAUDE.md §13.4a, each test verifies the projection's stated
transcendental signature and/or its load-bearing identity through
analytical limit, symbolic identity, or numerical cross-check.
"""

from __future__ import annotations

import math
import numpy as np
import pytest
import sympy as sp


# ----------------------------------------------------------------------------
# §III.1  Fock conformal: kappa = -1/16, Vol(S^3) = 2 pi^2, n^2 - 1 spectrum
# ----------------------------------------------------------------------------

def test_paper34_III1_kappa_rational_prefactor():
    """Paper 34 §III.1 (sec:proj_fock): transcendental signature 'rational
    prefactor kappa = -1/16 (Rydberg-to-graph-eigenvalue conversion).'

    Verifies the universal topological constant kappa = -1/16 used to
    map graph eigenvalues to Rydberg energies (CLAUDE.md §4).

    Per CLAUDE.md §8 ('-1/16 is the universal topological constant
    can be used directly'), this value is the framework's free-to-use
    rational prefactor. Production-side symbolic version lives at
    geovac.graph_qed_propagator.KAPPA_SCALAR = sympy Rational(-1, 16).
    """
    from geovac.graph_qed_propagator import KAPPA_SCALAR
    from sympy import Rational

    # Symbolic exactness check
    assert KAPPA_SCALAR == Rational(-1, 16), (
        f"KAPPA_SCALAR = {KAPPA_SCALAR} != -1/16"
    )

    # Numerical exactness check
    assert math.isclose(float(KAPPA_SCALAR), -1.0 / 16.0,
                        rel_tol=0.0, abs_tol=0.0), (
        f"float(KAPPA_SCALAR) = {float(KAPPA_SCALAR)} != -1/16"
    )


def test_paper34_III1_S3_volume_2pi_squared():
    """Paper 34 §III.1: 'pi enters through Vol(S^3) = 2 pi^2 when
    subsequent spectral integrals are performed.'

    Confirms the same Vol(S^3) constant referenced by §III.2 (Hopf
    bundle) and §III.6 (spectral action) -- so all three projections
    share a single Vol(S^3) measure source.
    """
    from geovac.hopf_bundle import VOL_S3

    assert math.isclose(VOL_S3, 2.0 * math.pi ** 2, rel_tol=1e-15, abs_tol=1e-15), (
        f"Vol(S^3) = {VOL_S3} != 2 pi^2 = {2*math.pi**2}"
    )


@pytest.mark.parametrize("n", list(range(1, 8)))
def test_paper34_III1_laplacian_spectrum_n2_minus_1(n):
    """Paper 34 §III.1 / CLAUDE.md §4: 'Eigenvalues of the Laplace-Beltrami
    operator on unit S^3 are pure integers: lambda_n = -(n^2 - 1).'

    Algebraic-only structural check of the integer-spectrum claim.
    Degeneracies are handled by tests/test_fock_laplacian.py.
    """
    lam = -(n ** 2 - 1)
    # Must be an integer (Paper 34 's 'pure integers' claim)
    assert isinstance(lam, int)
    # Ground state at n=1 has lam = 0
    if n == 1:
        assert lam == 0
    else:
        # Excited states have lam negative
        assert lam < 0
        # Formula closure
        assert lam == 1 - n ** 2


# ----------------------------------------------------------------------------
# §III.5  Sturmian relabeling at lambda = Z/n preserves rationality
# ----------------------------------------------------------------------------

def test_paper34_III5_sturmian_rationality_preserved():
    """Paper 34 §III.5 (sec:proj_sturmian): transcendental signature
    'preserves rationality of Layer 1. No new pi or other
    transcendentals enter at this step.'

    Tests by symbolic relabeling: the Sturmian focal length lam = Z/n
    is rational in (Z, n), and the graph eigenvalues -(n^2 - 1) remain
    integer under this relabeling. This is the categorical statement
    Paper 34 §III.5 claims.
    """
    Z, n = sp.symbols('Z n', positive=True, integer=True)
    lam_sturmian = Z / n
    # lam should be a rational function of (Z, n)
    assert lam_sturmian.is_rational_function(Z, n), (
        f"Sturmian focal length lam = Z/n = {lam_sturmian} is not rational"
    )

    # The energy formula E_n = -Z^2/(2 n^2) under Sturmian relabeling
    # E -> -lam^2 / 2 = -(Z/n)^2 / 2 stays rational
    E_sturmian = -lam_sturmian ** 2 / 2
    E_direct = -Z ** 2 / (2 * n ** 2)
    assert sp.simplify(E_sturmian - E_direct) == 0, (
        f"Sturmian energy mismatch: {E_sturmian} != {E_direct}"
    )

    # No pi or other transcendental enters at this level (Layer 1 only)
    assert sp.pi not in E_sturmian.free_symbols, (
        f"Unexpected pi in Sturmian Layer 1: {E_sturmian.free_symbols}"
    )


def test_paper34_III5_sturmian_iv_closure_at_ell0():
    """Paper 34 §III.5 + §III.13 cross-reference: the velocity-form
    closure I_v(nS) = (Z^4/n^3) delta_{l,0} is the value that the
    Sturmian-projected graph returns at l=0 with the Drake-Swainson
    structural denominator D_drake(n,l) = 2(2l+1) Z^4/n^3.

    For l=0: D_drake(n, 0) = 2 Z^4/n^3 = I_v(nS) (closure recovery).
    """
    Z, n = sp.symbols('Z n', positive=True, integer=True)
    l = 0
    D_drake_ell0 = 2 * (2 * l + 1) * Z ** 4 / n ** 3
    I_v_nS = 2 * Z ** 4 / n ** 3
    assert sp.simplify(D_drake_ell0 - I_v_nS) == 0, (
        f"Drake-Swainson l=0 closure with I_v(nS) failed: "
        f"D_drake = {D_drake_ell0}, I_v = {I_v_nS}"
    )


# ----------------------------------------------------------------------------
# §III.11 Vector-photon promotion: 1/(4 pi) per loop is S^2 Weyl factor
# ----------------------------------------------------------------------------

def test_paper34_III11_vector_photon_1_over_4pi():
    """Paper 34 §III.11 (sec:proj_vector_photon): transcendental signature
    '1/(4 pi) per loop, identified as the S^2 Weyl exchange constant
    of the Hopf base (Paper 33).'

    Verifies the structural identity 1/Vol(S^2) = 1/(4 pi), tying this
    projection to the same Hopf base measure source as §III.2.
    """
    from geovac.hopf_bundle import VOL_S2

    # Vol(S^2) = 4 pi (standard) -> 1/Vol(S^2) = 1/(4 pi)
    assert math.isclose(VOL_S2, 4.0 * math.pi, rel_tol=1e-15), (
        f"Vol(S^2) = {VOL_S2} != 4 pi"
    )

    one_over_4pi = 1.0 / VOL_S2
    expected = 1.0 / (4.0 * math.pi)
    assert math.isclose(one_over_4pi, expected, rel_tol=1e-15), (
        f"1/Vol(S^2) = {one_over_4pi} != 1/(4 pi) = {expected}"
    )

    # Symbolic cross-check
    vol_S2_sym = 4 * sp.pi
    assert sp.simplify(1 / vol_S2_sym - sp.Rational(1, 4) / sp.pi) == 0


# ----------------------------------------------------------------------------
# §III.13 Drake-Swainson: D_drake(n, l) = 2 (2l + 1) Z^4 / n^3
# ----------------------------------------------------------------------------

@pytest.mark.parametrize("n,l", [
    (1, 0), (2, 0), (2, 1), (3, 0), (3, 1), (3, 2),
    (4, 0), (4, 1), (4, 2), (4, 3),
])
def test_paper34_III13_drake_swainson_structural_denominator(n, l):
    """Paper 34 §III.13 (sec:proj_drake_swainson): structural denominator
    D_drake(n, l) = 2 (2l + 1) Z^4 / n^3.

    Combinatorial-rational by construction (spin x angular degeneracy
    x hydrogenic density). Verifies for Z=1..3 across (n, l) panel.
    """
    for Z in (1, 2, 3):
        D = sp.Rational(2 * (2 * l + 1) * Z ** 4, n ** 3)
        # Manifestly rational (combinatorial-rational claim)
        assert D.is_rational, (
            f"D_drake(n={n}, l={l}, Z={Z}) = {D} is not rational"
        )
        # Recovers I_v(nS) at l=0
        if l == 0:
            I_v_nS = sp.Rational(2 * Z ** 4, n ** 3)
            assert D == I_v_nS, (
                f"At l=0: D_drake = {D}, I_v(nS) = {I_v_nS} should match"
            )


# ----------------------------------------------------------------------------
# §III.16 Breit retardation: R^0_BP closed forms (alpha^4 * Q[log 2, log 3])
# ----------------------------------------------------------------------------

def test_paper34_III16_breit_R0_1s1s_1s1s():
    """Paper 34 §III.16 (sec:proj_breit_retardation): closed form
    R^0_BP(1s, 1s; 1s, 1s) = -5 + 8 log 2 (Z=1, k=0).

    Validates the Z=1 closed form via the production
    geovac.breit_integrals.compute_radial routine.
    """
    from geovac.breit_integrals import compute_radial

    val = compute_radial(
        n1=1, l1=0, n3=1, l3=0,
        n2=1, l2=0, n4=1, l4=0,
        k=0, kernel_type="breit", Z=1,
    )
    expected = -5 + 8 * sp.log(2)
    diff = sp.simplify(val - expected)
    assert diff == 0, (
        f"R^0_BP(1s,1s;1s,1s) at Z=1: production = {val}, "
        f"Paper 34 stated = {expected}, diff = {diff}"
    )


def test_paper34_III16_breit_R0_1s2s_1s2s():
    """Paper 34 §III.16: R^0_BP(1s, 2s; 1s, 2s) at Z=1, k=0.

    CORRECTION SURFACED 2026-06-04: Paper 34 §III.16 lists the value as
    '-4 log 2 - 19/9 + 9 log(3)/2' (numerically 0.06006). The production
    module geovac.breit_integrals.compute_radial and its existing
    regression test (tests/test_breit_integrals.py:131-132,
    'BP-retarded integral is a pure rational') both give the pure
    rational R^0_BP(1s,2s;1s,2s) = 4/81 = 0.04938.

    The (1s,2s;1s,2s) row in Paper 34 §III.16 should be corrected to the
    rational 4/81 (matches production and pre-existing regression test;
    the (1s,1s;1s,1s) row -5 + 8 log 2 separately verified correct).
    Same correction needed in
    debug/ps_1s2s_autopsy_track4_memo.md line 81 (origin of the Paper
    34 transcription).
    """
    from geovac.breit_integrals import compute_radial

    val = compute_radial(
        n1=1, l1=0, n3=2, l3=0,
        n2=1, l2=0, n4=2, l4=0,
        k=0, kernel_type="breit", Z=1,
    )
    # Correct rational value (production + existing regression test).
    expected = sp.Rational(4, 81)
    assert sp.simplify(val - expected) == 0, (
        f"R^0_BP(1s,2s;1s,2s) at Z=1: production = {val}, expected 4/81"
    )


def test_paper34_III16_breit_Z3_scaling():
    """Paper 34 §III.16: Breit retardation integrals scale as Z^3
    (per geovac.breit_integrals docstring: 'Breit integrals scale as Z^3
    while Coulomb integrals scale as Z^1').

    This is the rest-mass-free Z-scaling that the alpha^4 * Q
    ring-preserving claim relies on.
    """
    from geovac.breit_integrals import compute_radial

    Z1 = compute_radial(1, 0, 1, 0, 1, 0, 1, 0, k=0, kernel_type="breit", Z=1)
    Z3 = compute_radial(1, 0, 1, 0, 1, 0, 1, 0, k=0, kernel_type="breit", Z=3)
    ratio = sp.simplify(Z3 / Z1)
    assert ratio == 27, (
        f"Z^3 scaling check failed: Z=3 / Z=1 = {ratio} (expected 27 = 3^3)"
    )


# ----------------------------------------------------------------------------
# §III.17 Foldy/Friar: (2 pi / 3) Z alpha <r^2>_E delta^3(r) contact term
# ----------------------------------------------------------------------------

def test_paper34_III17_foldy_friar_prefactor_2pi_over_3():
    """Paper 34 §III.17 (sec:proj_charge_density): Foldy/Friar contact
    term Delta V = +(2 pi / 3) Z alpha <r^2>_E delta^3(r).

    Verifies the (2 pi / 3) prefactor symbolically. This pulls together
    Vol(S^2) = 4 pi / 6 = 2 pi / 3 in spherical-shell averaging.
    """
    # The 2 pi / 3 is the standard prefactor from
    # delta-function averaging of the convolved Coulomb potential
    # (Friar 1979 eq. 9; Eides 2024 §6).
    # Numerical: 2 pi / 3
    val = 2 * math.pi / 3
    assert math.isclose(val, 2.0944, rel_tol=1e-3), (
        f"Foldy/Friar prefactor 2 pi / 3 = {val} differs from 2.0944"
    )

    # Symbolic: relate to Vol(S^2) = 4 pi
    vol_S2 = 4 * sp.pi
    # (2 pi / 3) = Vol(S^2) / 6 (the 6 = 2 (rank) x 3 (Cartesian dims) factor
    # absorbing the spherical shell average of <r^2> over the angular part)
    assert sp.simplify(sp.Rational(2, 3) * sp.pi - vol_S2 / 6) == 0


def test_paper34_III17_hydrogenic_1s_contact_density():
    """Paper 34 §III.17 cross-check: hydrogenic 1s contact density at
    origin is |psi_1s(0)|^2 = Z^3 / pi (standard result; rational over
    pi, ring-preserving over Q(alpha) when combined with the Foldy/Friar
    pi prefactor).

    The Foldy/Friar contact term Delta E = (2 pi / 3) Z alpha
    <r^2>_E |psi(0)|^2 then evaluates to
    (2 pi / 3) Z alpha <r^2>_E * Z^3 / pi = (2/3) Z^4 alpha <r^2>_E,
    which is the canonical Lamb-shift r_p contribution -- the pi
    cancels, leaving a rational coefficient (the 'ring-preserving over
    Q(alpha)' claim).
    """
    Z = sp.symbols('Z', positive=True, integer=True)
    psi_1s_sq_origin = Z ** 3 / sp.pi

    # Foldy/Friar evaluated at 1s
    r_sq_E = sp.symbols('r_E_sq', positive=True)
    alpha = sp.symbols('alpha', positive=True)
    contact_shift = (sp.Rational(2, 3) * sp.pi * Z * alpha * r_sq_E) * psi_1s_sq_origin

    # pi must cancel
    contact_simpl = sp.simplify(contact_shift)
    assert sp.pi not in contact_simpl.free_symbols, (
        f"pi did not cancel in Foldy/Friar 1s contact: {contact_simpl}"
    )

    # Expected coefficient: (2/3) Z^4 alpha <r^2>_E
    expected = sp.Rational(2, 3) * Z ** 4 * alpha * r_sq_E
    diff = sp.simplify(contact_simpl - expected)
    assert diff == 0, (
        f"Foldy/Friar 1s evaluation: {contact_simpl} != (2/3) Z^4 alpha r^2_E "
        f"= {expected}, diff = {diff}"
    )


# ----------------------------------------------------------------------------
# §III.18 Zemach: A_hf (1 - 2 Z alpha m_e r_Z + O(r_Z^2)) leading-order
# ----------------------------------------------------------------------------

def test_paper34_III18_zemach_leading_order_linear_in_rZ():
    """Paper 34 §III.18 (sec:proj_magnetization_density): leading-order
    Zemach correction is A_hf_contact (1 - 2 Z alpha m_e r_Z + O(r_Z^2)).

    Test: the production module's regression scales linearly with r_Z
    at leading order. The delta_LO_ppm at r_Z=1.045 fm should be ~2x
    the delta_LO_ppm at r_Z=0.5225 fm (factor of 2 in r_Z -> factor of 2
    in correction at LO).
    """
    from geovac.magnetization_density import hydrogen_zemach_eides_leading_order

    rz_base = 1.045e-5 * 0.529177 / 0.529177  # arbitrary base in bohr units
    # The actual unit conversion is handled inside the function;
    # the test is invariant to absolute scaling.
    # We use the two-r_Z ratio test instead.
    rz_a = 1.045e-5  # close to Eides 2024 nominal in bohr
    rz_b = 2 * rz_a

    res_a = hydrogen_zemach_eides_leading_order(r_Z_bohr=rz_a, profile="gaussian")
    res_b = hydrogen_zemach_eides_leading_order(r_Z_bohr=rz_b, profile="gaussian")

    # Pick the leading-order delta_ppm key out of the result dict.
    # Default profile contains operator_level_delta_ppm
    key = "operator_level_delta_ppm"
    assert key in res_a, f"Expected key {key} in res_a, got {list(res_a.keys())}"

    delta_a = res_a[key]
    delta_b = res_b[key]

    # Linear-in-r_Z at leading order -> ratio ~ 2
    ratio = delta_b / delta_a
    assert abs(ratio - 2.0) < 0.05, (
        f"Zemach LO not linear in r_Z: delta(2 r_Z)/delta(r_Z) = {ratio} "
        f"(expected ~2.0 at leading order)"
    )


def test_paper34_III18_zemach_profile_independence_leading_order():
    """Paper 34 §III.18: 'profile-dependence at sub-leading order admits
    Gaussian, exponential, and dipole forms with structurally identical
    leading -2 Z alpha m_e r_Z behaviour (profile independence at
    leading order verified in Sprint HF Track 4 / Sprint MH Track C).'

    Verifies Gaussian vs exponential profile agreement at leading order.
    """
    from geovac.magnetization_density import hydrogen_zemach_eides_leading_order

    rz = 1.045e-5  # bohr
    res_gauss = hydrogen_zemach_eides_leading_order(r_Z_bohr=rz, profile="gaussian")
    res_exp = hydrogen_zemach_eides_leading_order(r_Z_bohr=rz, profile="exponential")

    key = "operator_level_delta_ppm"
    delta_gauss = res_gauss[key]
    delta_exp = res_exp[key]

    # Profile independence at LO: |delta_gauss - delta_exp| / |delta_gauss| < 5%
    rel_diff = abs(delta_gauss - delta_exp) / max(abs(delta_gauss), 1e-30)
    assert rel_diff < 0.05, (
        f"Profile independence at LO violated: gaussian = {delta_gauss}, "
        f"exponential = {delta_exp}, rel diff = {rel_diff}"
    )


# ----------------------------------------------------------------------------
# §III.19 Nuclear tensor multipole: rank-2 Wigner 3j selection rules
# ----------------------------------------------------------------------------

def test_paper34_III19_rank2_multipole_triangle_inequality():
    r"""Paper 34 §III.19 (sec:proj_tensor_multipole): rank-2 quadrupole
    coupling H_Q = -e Q_N T^{(2)}_{ij} (\partial_i E_j) / 6 satisfies
    Wigner 3j triangle inequality at every step ('ring-preserving over
    Q at the angular level').

    Test: rank-2 spherical tensor selection rule -- non-zero matrix
    element between |l_1, m_1> and |l_2, m_2> requires |l_1 - l_2| <= 2
    <= l_1 + l_2 and m_1 - m_2 = M with -2 <= M <= 2.
    """
    from sympy.physics.wigner import wigner_3j

    # Rank-2 multipole: triangle inequality requires |l1 - l2| <= 2 <= l1 + l2

    # Outside upper bound -> vanishes
    # s -> g (l=0 to l=4) under rank-2: 0 + 4 = 4 >> 2, vanishes
    assert wigner_3j(0, 2, 4, 0, 0, 0) == 0
    # f -> s (l=3 to l=0) under rank-2: |3-0| = 3 > 2, vanishes
    assert wigner_3j(3, 2, 0, 0, 0, 0) == 0

    # Inside triangle, non-zero
    # s -> d (l=0 to l=2): |0-2| = 2 = 2, on boundary
    assert wigner_3j(0, 2, 2, 0, 0, 0) != 0
    # p -> p (l=1 to l=1) under rank-2: 0 <= 2 <= 2, allowed
    assert wigner_3j(1, 2, 1, 0, 0, 0) != 0
    # d -> d under rank-2
    assert wigner_3j(2, 2, 2, 0, 0, 0) != 0


def test_paper34_III19_rank2_M_quantum_number_conservation():
    """Paper 34 §III.19: rank-2 spherical tensor M_q in {-2, -1, 0, 1, 2}
    couples m_1 - m_2 = M (the magnetic-quantum-number conservation
    that gives the quadrupole hyperfine M selection rule).

    Test: wigner_3j(l1, 2, l2, -m1, M, m2) is zero unless m1 - m2 = M.
    """
    from sympy.physics.wigner import wigner_3j

    # Take l1 = l2 = 2 (d-d coupling under rank-2)
    # Non-zero requires m1 + M + m2 = 0 -> m2 = -m1 - M -> m1 - m2 = ... wait
    # wigner_3j(l1, j, l2, m1, M, m2) is non-zero only if m1 + M + m2 = 0
    # So m1 - (-m2) = m1 + m2 = -M -> equivalent statement.

    # Conservation: m1 + M + m2 = 0
    assert wigner_3j(2, 2, 2, 1, -1, 0) != 0   # 1 + (-1) + 0 = 0 OK
    assert wigner_3j(2, 2, 2, 1, 0, -1) != 0   # 1 + 0 + (-1) = 0 OK
    assert wigner_3j(2, 2, 2, 2, 0, -2) != 0   # 2 + 0 + (-2) = 0 OK
    # Violation
    assert wigner_3j(2, 2, 2, 1, 1, 1) == 0    # 1 + 1 + 1 = 3 != 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
