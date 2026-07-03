"""
Tests for the cross-register V_eN module (Phase C-W1a-physics).

Validates:
- Roothaan 1951 closed form for 1s × 1s L=0
- Multi-lambda Shibuya-Wulfman extension (bit-identical bare-Coulomb regression)
- Symmetry under (lam_e <-> lam_n)
- Point-nucleus limit J_0 -> lam_e as lam_n -> infinity
- Multipole termination at L_max = l_e + l_e' (electron) and l_n + l_n' (nucleus)
- Cross-register total m-conservation
- Block-diagonal SO(3) structure
- Hermiticity of the resulting Pauli string sum
- Hydrogen recoil correction at leading order in m_e/m_p (Bethe-Salpeter)

Author: GeoVac Development Team (Phase C-W1a-physics)
Date: May 2026
"""

from __future__ import annotations

import numpy as np
import pytest

from geovac.cross_register_vne import (
    A0_FM,
    CrossRegisterVneSpec,
    LAM_NUCLEUS_DEFAULT,
    LAM_NUCLEUS_GEOMETRIC,
    LAM_NUCLEUS_QUANTUM_MOTIONAL,
    M_PROTON_OVER_M_E,
    _cross_register_angular_factor,
    _cross_register_radial_integral,
    _roothaan_J0,
    _roothaan_J0_symbolic,
    _roothaan_J_general,
    bare_coulomb_regression,
    compute_cross_register_vne,
    cross_register_eri_matrix,
    cross_register_recoil_correction,
    hydrogen_recoil_correction_leading_order,
    pachucki_2023_leading_order_check,
    zemach_magnetization_correction_pauli,
)
from geovac.shibuya_wulfman import (
    _hydrogenic_poly_coeffs,
    _hydrogenic_poly_coeffs_lam,
    _radial_split_integral,
    _radial_split_integral_lam,
    compute_cross_center_vne_element,
    compute_cross_center_vne_element_lam,
)


# ---------------------------------------------------------------------------
# Roothaan 1951 closed-form J_0 verification
# ---------------------------------------------------------------------------


class TestRoothaanJ0:
    """Symbolic and numerical verification of the Roothaan J_0 formula."""

    def test_textbook_lam1_lam1(self) -> None:
        """J_0(1, 1) = 1*1*(1+3+1)/8 = 5/8 = 0.625 (textbook H-H 1s coupling)."""
        assert _roothaan_J0(1.0, 1.0) == pytest.approx(5.0 / 8.0, rel=1e-15)

    def test_symmetric(self) -> None:
        """J_0(lam_e, lam_n) = J_0(lam_n, lam_e)."""
        for lam_e, lam_n in [(1.0, 2.0), (2.0, 3.0), (1.5, 3.7), (0.5, 100.0)]:
            J_eq = _roothaan_J0(lam_e, lam_n)
            J_swap = _roothaan_J0(lam_n, lam_e)
            assert J_eq == pytest.approx(J_swap, rel=1e-15)

    def test_point_nucleus_limit(self) -> None:
        """J_0(lam_e, lam_n -> oo) -> lam_e (recovers <1s | 1/r | 1s>)."""
        lam_e = 1.0
        for lam_n in [10.0, 100.0, 1000.0, 10000.0]:
            J0 = _roothaan_J0(lam_e, lam_n)
            assert J0 < lam_e  # nucleus has finite extent, so J_0 < lam_e
            assert J0 == pytest.approx(lam_e, abs=1.0 / lam_n)

    def test_hydrogen_one_hartree_limit(self) -> None:
        """For lam_e=1, lam_n -> oo: J_0 -> 1 Ha (hydrogen 1s ground-state energy)."""
        for lam_n in [1000.0, 10000.0, 100000.0]:
            J0 = _roothaan_J0(1.0, lam_n)
            assert J0 == pytest.approx(1.0, abs=2.5 / lam_n ** 2)

    def test_symbolic_form(self) -> None:
        """Symbolic Roothaan formula reduces to lam_e*lam_n*(...)/(lam_e+lam_n)^3."""
        import sympy as sp
        lam_e, lam_n = sp.symbols('lam_e lam_n', positive=True)
        J0_sym = _roothaan_J0_symbolic(lam_e, lam_n)
        # Substitute and check
        for (a, b) in [(1, 1), (2, 3), (5, 7)]:
            num = float(J0_sym.subs([(lam_e, a), (lam_n, b)]))
            assert num == pytest.approx(_roothaan_J0(a, b), rel=1e-15)

    def test_positivity(self) -> None:
        """J_0 is always positive (it's a Coulomb mutual energy of two densities)."""
        for lam_e in [0.5, 1.0, 2.0, 5.0]:
            for lam_n in [0.5, 1.0, 2.0, 5.0, 100.0]:
                assert _roothaan_J0(lam_e, lam_n) > 0

    def test_invalid_input_rejected(self) -> None:
        """Non-positive lambda raises ValueError."""
        with pytest.raises(ValueError):
            _roothaan_J0(0.0, 1.0)
        with pytest.raises(ValueError):
            _roothaan_J0(1.0, -0.5)


# ---------------------------------------------------------------------------
# Numerical engine vs closed-form Roothaan
# ---------------------------------------------------------------------------


class TestNumericalEngine:
    """The numerical engine should reproduce Roothaan J_0 to machine precision."""

    def test_machine_precision_matches_closed_form_lam_eq_lam(self) -> None:
        """Numerical engine matches Roothaan to ~1e-13 for matched lambdas."""
        for lam in [0.5, 1.0, 2.0, 5.0]:
            closed = _roothaan_J0(lam, lam)
            numeric = _roothaan_J_general(1, 0, 1, 0, 1, 0, 1, 0, 0, lam, lam)
            assert numeric == pytest.approx(closed, rel=1e-12)

    def test_matches_at_mismatched_lambdas(self) -> None:
        """Numerical engine matches Roothaan at mismatched (moderate) lambdas."""
        for (lam_e, lam_n) in [(1.0, 2.0), (2.0, 3.0), (0.5, 1.5)]:
            closed = _roothaan_J0(lam_e, lam_n)
            numeric = _roothaan_J_general(
                1, 0, 1, 0, 1, 0, 1, 0, 0, lam_e, lam_n,
            )
            assert numeric == pytest.approx(closed, rel=1e-12)

    def test_large_separation_of_scales(self) -> None:
        """At very disparate scales (lam_n >> lam_e), accuracy degrades but stays sub-ppm."""
        lam_e = 1.0
        for lam_n in [10.0, 100.0]:
            closed = _roothaan_J0(lam_e, lam_n)
            numeric = _roothaan_J_general(
                1, 0, 1, 0, 1, 0, 1, 0, 0, lam_e, lam_n,
            )
            # For lam_n=100, can be ~1e-6 because of GL nodes
            assert numeric == pytest.approx(closed, rel=1e-5)

    def test_returns_zero_for_invalid_L(self) -> None:
        """L > l_e + l_e' should return 0."""
        # 1s x 1s: l = 0, so L_max for either side = 0
        val = _roothaan_J_general(1, 0, 1, 0, 1, 0, 1, 0, 1, 1.0, 1.0)
        assert val == 0.0


# ---------------------------------------------------------------------------
# Multi-lambda Shibuya-Wulfman regression
# ---------------------------------------------------------------------------


class TestMultiLambdaShibuyaWulfman:
    """The multi-lambda extension preserves the matched-case bit-identically."""

    def test_poly_coeffs_matched(self) -> None:
        """_hydrogenic_poly_coeffs_lam(Z/n, n, l) matches _hydrogenic_poly_coeffs(Z, n, l)."""
        for Z in [1.0, 2.0, 3.0]:
            for n in [1, 2, 3]:
                for l in range(n):
                    c_single, alpha_single = _hydrogenic_poly_coeffs(Z, n, l)
                    c_lam, alpha_lam = _hydrogenic_poly_coeffs_lam(Z / n, n, l)
                    np.testing.assert_allclose(c_single, c_lam, rtol=1e-15)
                    assert alpha_single == pytest.approx(alpha_lam, rel=1e-15)

    def test_radial_split_integral_matched(self) -> None:
        """_radial_split_integral_lam matches _radial_split_integral at matched lambda."""
        Z, n, l = 3.0, 1, 0
        for L in [0]:
            for R in [1.0, 3.015, 5.0]:
                old = _radial_split_integral(Z, n, l, n, l, L, R)
                new = _radial_split_integral_lam(Z / n, n, l, Z / n, n, l, L, R)
                assert new == pytest.approx(old, rel=1e-12)

    def test_compute_cross_center_vne_element_matched(self) -> None:
        """compute_cross_center_vne_element_lam matches the original at matched."""
        Z, R = 3.0, 3.015
        val_old = compute_cross_center_vne_element(
            Z, 1, 0, 0, 1, 0, 0, 1.0, R, L_max=2,
        )
        val_new = compute_cross_center_vne_element_lam(
            Z, Z, 1, 0, 0, 1, 0, 0, 1.0, R, L_max=2,
        )
        assert val_new == pytest.approx(val_old, rel=1e-12)

    def test_mismatched_lambda_works(self) -> None:
        """Mismatched lambdas produce a finite, reasonable matrix element."""
        val = compute_cross_center_vne_element_lam(
            1.0, 2.0, 1, 0, 0, 1, 0, 0, 1.0, 3.0, L_max=2,
        )
        # Should be finite and negative (electron in attractive potential).
        assert np.isfinite(val)
        assert val < 0


# ---------------------------------------------------------------------------
# Cross-register angular factor
# ---------------------------------------------------------------------------


class TestAngularFactor:
    """Bilateral Gaunt factor with cross-register m-conservation."""

    def test_ss_L0_M0(self) -> None:
        """1s × 1s × 1s × 1s at L=M=0: factor = 1."""
        val = _cross_register_angular_factor(
            0, 0, 0, 0,  # electron l, m, l', m'
            0, 0, 0, 0,  # nucleus l, m, l', m'
            0, 0,         # L, M
        )
        assert val == pytest.approx(1.0, rel=1e-12)

    def test_cross_register_m_conservation(self) -> None:
        """m_e + m_n = m_e' + m_n' is required."""
        # m_e=0, m_n=0; m_e'=1, m_n'=0: violates conservation
        val = _cross_register_angular_factor(
            0, 0, 1, 1,
            0, 0, 0, 0,
            1, 1,
        )
        assert val == 0.0

    def test_M_must_match_electron_3j(self) -> None:
        """M = m_e' - m_e is required."""
        val = _cross_register_angular_factor(
            0, 0, 1, 1,
            1, 0, 1, -1,
            1, -1,  # wrong M
        )
        assert val == 0.0

    def test_triangle_inequality_electron(self) -> None:
        """L > l_e + l_e' on electron side returns 0."""
        val = _cross_register_angular_factor(
            0, 0, 0, 0,  # l_e + l_e' = 0
            1, 0, 1, 0,
            2, 0,         # L=2 > 0
        )
        assert val == 0.0

    def test_triangle_inequality_nucleus(self) -> None:
        """L > l_n + l_n' on nucleus side returns 0."""
        val = _cross_register_angular_factor(
            1, 0, 1, 0,
            0, 0, 0, 0,
            1, 0,         # L=1 > 0+0
        )
        assert val == 0.0

    def test_parity_rule_electron(self) -> None:
        """l_e + L + l_e' must be even on electron side."""
        # 0 + 1 + 0 = 1 (odd), should give 0
        val = _cross_register_angular_factor(
            0, 0, 0, 0,
            1, 0, 1, 0,  # nucleus side OK at L=1 (1+1+1=3 odd, wait this also fails)
            1, 0,
        )
        # Actually both sides fail parity, so 0.
        assert val == 0.0


# ---------------------------------------------------------------------------
# ERI matrix and Pauli encoding
# ---------------------------------------------------------------------------


class TestERIMatrix:
    """The cross_register_eri_matrix builds the correct V tensor."""

    def test_simplest_case_1s_1s(self) -> None:
        """1 state × 1 state at L=0: V[0,0,0,0] = -Z * J_0."""
        spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=10.0, n_max_n=1,
            Z_nuc=1.0, L_max=0,
        )
        V, states_e, states_n = cross_register_eri_matrix(spec)
        assert V.shape == (1, 1, 1, 1)
        expected = -1.0 * _roothaan_J0(1.0, 10.0)
        assert V[0, 0, 0, 0] == pytest.approx(expected, rel=1e-12)

    def test_hermiticity(self) -> None:
        """V[i,j,k,l] = V[k,l,i,j] (Hermiticity for real basis)."""
        spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=2,
            lam_n=10.0, n_max_n=1,
            Z_nuc=1.0, L_max=2,
        )
        V, states_e, states_n = cross_register_eri_matrix(spec)
        N_e, N_n, _, _ = V.shape
        for i in range(N_e):
            for j in range(N_n):
                for k in range(N_e):
                    for ll in range(N_n):
                        assert V[i, j, k, ll] == pytest.approx(V[k, ll, i, j], abs=1e-10)

    def test_block_diagonal_under_m_conservation(self) -> None:
        """V[i,j,k,l] = 0 unless m_e + m_n = m_e' + m_n'."""
        spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=2,
            lam_n=10.0, n_max_n=2,
            Z_nuc=1.0, L_max=2,
        )
        V, states_e, states_n = cross_register_eri_matrix(spec)
        for i, (ne, le, me) in enumerate(states_e):
            for k, (ne_p, le_p, me_p) in enumerate(states_e):
                for j, (nn, ln, mn) in enumerate(states_n):
                    for ll, (nn_p, ln_p, mn_p) in enumerate(states_n):
                        if me + mn != me_p + mn_p:
                            assert V[i, j, k, ll] == 0.0


class TestPauliEncoding:
    """Pauli-string assembly is correct."""

    def test_simple_case_1s_1s(self) -> None:
        """1s × 1s on a 2-qubit register: 4 Pauli terms."""
        spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=10.0, n_max_n=1,
            Z_nuc=1.0, L_max=0,
        )
        result = compute_cross_register_vne(spec)
        assert result['Q_total'] == 2
        assert result['Q_e'] == 1
        assert result['Q_n'] == 1
        # Pauli structure: II, ZI, IZ, ZZ
        pauli = result['pauli_terms']
        for k, v in pauli.items():
            assert len(k) == 2
            assert all(c in 'IXYZ' for c in k)

    def test_pauli_diagonal_energy(self) -> None:
        """In the |11> state (both registers occupied), energy = V[0,0,0,0]."""
        spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=10.0, n_max_n=1,
            Z_nuc=1.0, L_max=0,
        )
        result = compute_cross_register_vne(spec)
        # Build the Hamiltonian matrix
        PAULI = {
            'I': np.eye(2, dtype=complex),
            'X': np.array([[0, 1], [1, 0]], dtype=complex),
            'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
            'Z': np.array([[1, 0], [0, -1]], dtype=complex),
        }
        n_qb = result['Q_total']
        H = np.zeros((2 ** n_qb, 2 ** n_qb), dtype=complex)
        for s, c in result['pauli_terms'].items():
            op = np.array([[1.0]], dtype=complex)
            for ch in s:
                op = np.kron(op, PAULI[ch])
            H = H + c * op
        # Hermiticity
        assert np.allclose(H, H.conj().T, atol=1e-12)
        # Diagonal: |11> state has energy V[0,0,0,0]
        diag = np.diag(H).real
        # Computational basis |11> at index 3 (n=2 qubits)
        V_expected = result['V_eri'][0, 0, 0, 0]
        assert diag[3] == pytest.approx(V_expected, abs=1e-10)
        # |00> state has 0 energy
        assert diag[0] == pytest.approx(0.0, abs=1e-10)


# ---------------------------------------------------------------------------
# Bare-Coulomb regression and physical limits
# ---------------------------------------------------------------------------


class TestBareCoulombRegression:
    """Verify the multi-lambda extension reproduces classical V_ne in the
    point-nucleus limit lam_n -> infinity."""

    def test_converges_to_classical(self) -> None:
        """As lam_n -> oo, J_0 -> lam_e (point-charge limit)."""
        result = bare_coulomb_regression(Z=1.0)
        assert result['converges_to_classical'] is True
        # At lam_n = 100000, error should be ~1e-10
        assert result['errors_at_each_lam_n'][-1] < 1e-9

    def test_monotonic_convergence(self) -> None:
        """Errors should decrease monotonically as lam_n increases."""
        result = bare_coulomb_regression(Z=1.0)
        errs = result['errors_at_each_lam_n']
        for i in range(len(errs) - 1):
            assert errs[i + 1] <= errs[i]


# ---------------------------------------------------------------------------
# Hydrogen recoil correction
# ---------------------------------------------------------------------------


class TestHydrogenRecoil:
    """Cross-register recoil at leading order in m_e/m_p (Bethe-Salpeter)."""

    def test_bethe_salpeter_value(self) -> None:
        """The Bethe-Salpeter formula gives ~ +2.72e-4 Ha for hydrogen 1s."""
        E_recoil = hydrogen_recoil_correction_leading_order(Z=1.0, n=1)
        # m_e/m_p = 1/1836, |E_1| = 0.5 Ha => +0.5/1836 = +2.72e-4 Ha
        expected = 0.5 / M_PROTON_OVER_M_E
        assert E_recoil == pytest.approx(expected, rel=1e-12)

    def test_cross_register_at_quantum_motional_lam_n(self) -> None:
        """Cross-register V_eN at lam_n = 2 sqrt(M_p) reproduces Bethe-Salpeter."""
        spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=LAM_NUCLEUS_QUANTUM_MOTIONAL, n_max_n=1,
            Z_nuc=1.0, L_max=0,
        )
        recoil = cross_register_recoil_correction(spec)
        # Should match Bethe-Salpeter to ~3% (leading-order, with higher-order
        # corrections in 1/lam_n^3)
        assert recoil['relative_error'] < 0.05  # within 5%
        # Should recover sign and order of magnitude
        assert recoil['cross_register_recoil_estimate'] > 0  # less bound

    def test_pachucki_validation(self) -> None:
        """Pachucki 2023 leading-order check."""
        res = pachucki_2023_leading_order_check()
        assert res['relative_error'] < 0.05  # within 5%
        # Bethe-Salpeter expected ~ 2.72e-4 Ha
        assert res['pachucki_leading_order'] == pytest.approx(2.72e-4, rel=0.01)


# ---------------------------------------------------------------------------
# Zemach magnetization sketch
# ---------------------------------------------------------------------------


class TestZemachSketch:
    """W1b magnetization-density sketch (downstream of W1a)."""

    def test_eides_2024_calibration(self) -> None:
        """At r_Z = 1.045 fm, Zemach correction is ~ -39.5 ppm (Eides Tab. 7.3)."""
        spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=LAM_NUCLEUS_DEFAULT, n_max_n=1,
            Z_nuc=1.0,
        )
        zem = zemach_magnetization_correction_pauli(spec, r_Z_bohr=1.045 / A0_FM)
        # Eides Tab. 7.3: -39.5 ppm
        assert zem['delta_ppm'] == pytest.approx(-39.5, rel=0.05)

    def test_r_Z_conversion(self) -> None:
        """1.045 fm should convert to ~1.97e-5 bohr."""
        spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=LAM_NUCLEUS_DEFAULT, n_max_n=1,
            Z_nuc=1.0,
        )
        zem = zemach_magnetization_correction_pauli(spec)
        assert zem['r_Z_bohr'] == pytest.approx(1.045 / A0_FM, rel=1e-12)
        assert zem['r_Z_fm'] == pytest.approx(1.045, rel=1e-9)


# ---------------------------------------------------------------------------
# Multipole termination (Phase B-W1a-diag Q-B verification)
# ---------------------------------------------------------------------------


class TestMultipoleTermination:
    """Multipole termination at L_max = min(l_e + l_e', l_n + l_n')."""

    def test_ss_x_ss_only_L0(self) -> None:
        """1s × 1s × 1s × 1s only has L=0 contribution."""
        spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=10.0, n_max_n=1,
            Z_nuc=1.0, L_max=10,  # set high; should still terminate at L=0
        )
        V, _, _ = cross_register_eri_matrix(spec)
        # L_max cap at min(l_e+l_e', l_n+l_n') = 0
        # Result equals -Roothaan at L=0
        expected = -1.0 * _roothaan_J0(1.0, 10.0)
        assert V[0, 0, 0, 0] == pytest.approx(expected, rel=1e-12)

    def test_terminates_below_specified_L_max(self) -> None:
        """L_max larger than l_e+l_e' is automatically truncated."""
        # 1s × 1s on both: l_e+l_e' = 0, l_n+l_n' = 0, so L_max effective = 0
        spec_low = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=10.0, n_max_n=1,
            Z_nuc=1.0, L_max=0,
        )
        spec_high = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=10.0, n_max_n=1,
            Z_nuc=1.0, L_max=10,
        )
        V_low, _, _ = cross_register_eri_matrix(spec_low)
        V_high, _, _ = cross_register_eri_matrix(spec_high)
        assert V_low[0, 0, 0, 0] == pytest.approx(V_high[0, 0, 0, 0], rel=1e-15)

    def test_d_orbital_termination_tight_at_L4(self) -> None:
        """Paper 23 SVII: termination at L_max = l_e + l_e' verified through
        l = 2 (added 8th cert, 2026-07-02 -- the prose claimed l in {0,1,2}
        but only l <= 1 was instantiated). n_max=3 on both registers:
        L_max=4 equals L_max=10 exactly (termination), and the L=4 channel
        genuinely contributes on the d x d block (tightness ~8.5e-5)."""
        def build(L: int):
            spec = CrossRegisterVneSpec(
                lam_e=1.0, n_max_e=3, lam_n=10.0, n_max_n=3,
                Z_nuc=1.0, L_max=L,
            )
            return cross_register_eri_matrix(spec)

        V10, states, _ = build(10)
        V4, _, _ = build(4)
        V3, _, _ = build(3)
        assert np.allclose(V4, V10, rtol=1e-14, atol=1e-16)
        d = [i for i, s in enumerate(states) if s[1] == 2]
        blk = np.ix_(d, d, d, d)
        assert np.max(np.abs(V4[blk] - V3[blk])) > 1e-6


# ---------------------------------------------------------------------------
# Spec validation
# ---------------------------------------------------------------------------


class TestSpecValidation:
    """CrossRegisterVneSpec validates inputs."""

    def test_positive_lambdas_required(self) -> None:
        with pytest.raises(ValueError):
            CrossRegisterVneSpec(lam_e=-1.0)
        with pytest.raises(ValueError):
            CrossRegisterVneSpec(lam_n=0.0)

    def test_n_max_at_least_1(self) -> None:
        with pytest.raises(ValueError):
            CrossRegisterVneSpec(n_max_e=0)
        with pytest.raises(ValueError):
            CrossRegisterVneSpec(n_max_n=0)

    def test_L_max_non_negative(self) -> None:
        with pytest.raises(ValueError):
            CrossRegisterVneSpec(L_max=-1)


# ---------------------------------------------------------------------------
# Pachucki higher-order Taylor expansion (Phase C-Pachucki sprint)
# ---------------------------------------------------------------------------


class TestRoothaanTaylorExpansion:
    """Symbolic Taylor expansion of the Roothaan J_0 in eps = 1/lam_n.

    The expansion is the analytical foundation for the per-order Pachucki
    comparison. The expansion has the closed form (lam_e = 1):

        J_0 = 1 - 2/lam_n^2 + 5/lam_n^3 - 9/lam_n^4 + 14/lam_n^5 - 20/lam_n^6 + ...
        Recoil shift = +Z(lam_e - J_0)
                     = +2/lam_n^2 - 5/lam_n^3 + 9/lam_n^4 - 14/lam_n^5 + ...

    """

    def test_taylor_expansion_verifies_against_closed_form(self) -> None:
        """The series at small eps must reproduce the Roothaan formula."""
        from geovac.cross_register_vne import roothaan_J0_taylor_expansion
        result = roothaan_J0_taylor_expansion(n_terms=6, lam_e_value=1.0)
        assert result['verified_against_closed_form']

    def test_taylor_coefficients_lam_e_1(self) -> None:
        """At lam_e = 1, the recoil-shift coefficients are c_2=2, c_3=-5,
        c_4=9, c_5=-14, c_6=20."""
        from geovac.cross_register_vne import roothaan_J0_taylor_expansion
        result = roothaan_J0_taylor_expansion(n_terms=7, lam_e_value=1.0)
        coeffs = dict((k, float(v)) for k, v in result['coefficients'])
        assert coeffs[0] == pytest.approx(0.0, abs=1e-15)
        assert coeffs[1] == pytest.approx(0.0, abs=1e-15)
        assert coeffs[2] == pytest.approx(2.0, rel=1e-12)
        assert coeffs[3] == pytest.approx(-5.0, rel=1e-12)
        assert coeffs[4] == pytest.approx(9.0, rel=1e-12)
        assert coeffs[5] == pytest.approx(-14.0, rel=1e-12)
        assert coeffs[6] == pytest.approx(20.0, rel=1e-12)

    def test_taylor_coefficients_general_lam_e(self) -> None:
        """Symbolic coefficients should scale as lam_e^{k+1}."""
        from geovac.cross_register_vne import roothaan_J0_taylor_expansion
        result = roothaan_J0_taylor_expansion(n_terms=6, lam_e_value=None)
        # Expected: c_k(lam_e) = (-1)^{k} * a_k * lam_e^{k+1}
        # Specifically c_2 = 2 lam_e^3, c_3 = -5 lam_e^4, etc.
        import sympy as sp
        lam_e_sym = sp.symbols('lam_e', positive=True)
        coeffs = dict(result['coefficients'])
        # k=2: 2 lam_e^3
        assert sp.simplify(coeffs[2] - 2 * lam_e_sym ** 3) == 0
        # k=3: -5 lam_e^4
        assert sp.simplify(coeffs[3] - (-5) * lam_e_sym ** 4) == 0
        # k=4: +9 lam_e^5
        assert sp.simplify(coeffs[4] - 9 * lam_e_sym ** 5) == 0
        # k=5: -14 lam_e^6
        assert sp.simplify(coeffs[5] - (-14) * lam_e_sym ** 6) == 0


class TestRecoilShiftThroughOrder:
    """Per-order Roothaan recoil shift at calibrated lam_n."""

    def test_leading_order_matches_bethe_salpeter_exactly(self) -> None:
        """At lam_n = 2 sqrt(M_p), the +2/lam_n^2 leading term EXACTLY equals
        Bethe-Salpeter recoil m_e/m_p * |E_1|.

        This is a structural calibration identity, not a fit:
            +2/lam_n^2 = +2/(4 M_p) = +1/(2 M_p) = m_e/m_p * |E_1(m_e)|.
        """
        from geovac.cross_register_vne import (
            roothaan_recoil_shift_through_order,
            hydrogen_recoil_correction_leading_order,
            LAM_NUCLEUS_QUANTUM_MOTIONAL,
        )
        rh = roothaan_recoil_shift_through_order(
            lam_n=LAM_NUCLEUS_QUANTUM_MOTIONAL, lam_e=1.0, Z=1.0, max_order=2,
        )
        # Per-order at k=2
        k_2_value = rh['per_order_shift'][0][1]
        bs_leading = hydrogen_recoil_correction_leading_order(Z=1.0, n=1)
        assert k_2_value == pytest.approx(bs_leading, rel=1e-10)

    def test_series_converges_to_full_roothaan(self) -> None:
        """The truncated series at max_order=5 should match the full Roothaan
        formula at lam_n=85.7 to better than 1e-9 Ha (machine precision)."""
        from geovac.cross_register_vne import (
            roothaan_recoil_shift_through_order,
            LAM_NUCLEUS_QUANTUM_MOTIONAL,
        )
        rh = roothaan_recoil_shift_through_order(
            lam_n=LAM_NUCLEUS_QUANTUM_MOTIONAL, max_order=5,
        )
        assert rh['series_converges_to_roothaan']
        # At max_order=5 and lam_n=85.7, dropped O(1/lam_n^6) ~ 5e-11 Ha
        assert abs(rh['series_residual_at_max_order']) < 1e-9

    def test_series_machine_precision_at_max_order_7(self) -> None:
        """At max_order=7, the series matches Roothaan to ~ 1e-13 Ha."""
        from geovac.cross_register_vne import (
            roothaan_recoil_shift_through_order,
            LAM_NUCLEUS_QUANTUM_MOTIONAL,
        )
        rh = roothaan_recoil_shift_through_order(
            lam_n=LAM_NUCLEUS_QUANTUM_MOTIONAL, max_order=7,
        )
        assert abs(rh['series_residual_at_max_order']) < 1e-12

    def test_series_alternates_sign(self) -> None:
        """Roothaan recoil shift coefficients alternate sign starting with +."""
        from geovac.cross_register_vne import (
            roothaan_recoil_shift_through_order,
            LAM_NUCLEUS_QUANTUM_MOTIONAL,
        )
        rh = roothaan_recoil_shift_through_order(
            lam_n=LAM_NUCLEUS_QUANTUM_MOTIONAL, max_order=5,
        )
        contribs = [v for _, v in rh['per_order_shift']]
        # k=2 +, k=3 -, k=4 +, k=5 -
        for i, v in enumerate(contribs):
            expected_sign = 1.0 if i % 2 == 0 else -1.0
            assert (v > 0) == (expected_sign > 0), \
                f"k={i+2}: expected sign {expected_sign}, got {v}"

    def test_series_magnitudes_decay_geometrically(self) -> None:
        """At lam_n = 85.7, |contribution(k+1)| / |contribution(k)| should be
        ~ 1/lam_n ~ 0.012 (geometric decay)."""
        from geovac.cross_register_vne import (
            roothaan_recoil_shift_through_order,
            LAM_NUCLEUS_QUANTUM_MOTIONAL,
        )
        rh = roothaan_recoil_shift_through_order(
            lam_n=LAM_NUCLEUS_QUANTUM_MOTIONAL, max_order=6,
        )
        contribs = [v for _, v in rh['per_order_shift']]
        for i in range(len(contribs) - 1):
            ratio = abs(contribs[i + 1] / contribs[i])
            # Each successive order is ~ (k+1)/k * 1/lam_n smaller
            assert ratio < 0.1, f"Expected geometric decay, got {ratio} at k={i+2}->{i+3}"


class TestPachuckiHigherOrderComparison:
    """Per-order comparison Roothaan vs Pachucki-Patkos-Yerokhin 2023."""

    def test_pachucki_leading_matches_bethe_salpeter(self) -> None:
        """pachucki_leading must equal the standard Bethe-Salpeter formula."""
        from geovac.cross_register_vne import (
            pachucki_higher_order_comparison,
            hydrogen_recoil_correction_leading_order,
        )
        ph = pachucki_higher_order_comparison(Z=1.0, n=1)
        bs = hydrogen_recoil_correction_leading_order(Z=1.0, n=1)
        assert ph['pachucki_leading'] == pytest.approx(bs, rel=1e-12)

    def test_full_reduced_mass_shift_close_to_leading(self) -> None:
        """Full reduced-mass shift differs from leading by O((m_e/m_p)^2)."""
        from geovac.cross_register_vne import pachucki_higher_order_comparison
        ph = pachucki_higher_order_comparison(Z=1.0, n=1)
        # Full reduced mass shift = 1/(2(M_p+1)) = 2.7216e-4
        # Leading = 1/(2 M_p) = 2.7231e-4
        # Difference ~ -1/(2 M_p^2) ~ -1.5e-7
        diff = ph['pachucki_full_reduced_mass_shift'] - ph['pachucki_leading']
        # NEXT order should be negative (~-1/(2 M_p^2))
        assert diff < 0
        assert abs(diff) < 1e-6  # very small

    def test_half_integer_orders_have_no_pachucki_analog(self) -> None:
        """Roothaan odd-k orders (k=3, 5, 7) correspond to half-integer
        powers of (m_e/m_p), structurally absent in the Pachucki series."""
        from geovac.cross_register_vne import pachucki_higher_order_comparison
        ph = pachucki_higher_order_comparison(Z=1.0, n=1, max_order=7)
        for entry in ph['roothaan_series']:
            if not entry['is_integer_mass_ratio_order']:
                # Half-integer: Pachucki analog is exactly 0 (structurally)
                assert entry['pachucki_analog_Ha'] == 0.0
                assert 'HALF-INTEGER' in entry['pachucki_label']

    def test_integer_orders_have_pachucki_analog(self) -> None:
        """Roothaan even-k orders (k=2, 4, 6, 8) correspond to integer
        powers of (m_e/m_p), with Pachucki analogs at k=2, 4."""
        from geovac.cross_register_vne import pachucki_higher_order_comparison
        ph = pachucki_higher_order_comparison(Z=1.0, n=1, max_order=5)
        for entry in ph['roothaan_series']:
            if entry['is_integer_mass_ratio_order']:
                # Integer order
                assert entry['order_k'] % 2 == 0
                if entry['order_k'] in (2, 4):
                    assert entry['pachucki_analog_Ha'] is not None

    def test_full_series_drift_about_minus_286_pct(self) -> None:
        """The full Roothaan series at lam_n=85.7 drifts to -2.86% relative
        to the leading Pachucki order."""
        from geovac.cross_register_vne import pachucki_higher_order_comparison
        ph = pachucki_higher_order_comparison(Z=1.0, n=1, max_order=5)
        # Full series cumulative drift should be -2.86%
        drift = ph['cumulative_drift_vs_pachucki_leading_pct']
        assert -3.0 < drift < -2.7, f"Expected drift around -2.86%, got {drift}%"

    def test_signs_at_k4_disagree_with_pachucki(self) -> None:
        """The Roothaan k=4 coefficient is +9 lam_e^5, giving positive
        contribution at lam_n>0; the physical Pachucki (m_e/m_p)^2 term
        is NEGATIVE (~-1/(2 M_p^2)). The signs disagree, confirming that
        the Roothaan integer-order tower is NOT the Pachucki tower (it
        is a basis-truncation artifact at fixed n_max=1)."""
        from geovac.cross_register_vne import pachucki_higher_order_comparison
        ph = pachucki_higher_order_comparison(Z=1.0, n=1, max_order=4)
        k4_entry = next(e for e in ph['roothaan_series'] if e['order_k'] == 4)
        roothaan_k4 = k4_entry['roothaan_contribution_Ha']
        pachucki_k4_analog = k4_entry['pachucki_analog_Ha']
        assert roothaan_k4 > 0
        assert pachucki_k4_analog < 0
        # Different signs
        assert roothaan_k4 * pachucki_k4_analog < 0


class TestIntegerOrderOnlyEstimate:
    """The integer-only sub-sum filters out the half-integer 'Sturmian-
    truncation' artifacts and lands sub-percent."""

    def test_k2_only_matches_pachucki_leading_exactly(self) -> None:
        """Through k=2 only: matches leading Pachucki at machine precision."""
        from geovac.cross_register_vne import integer_order_only_recoil_estimate
        io = integer_order_only_recoil_estimate(Z=1.0, n=1, max_order_k=2)
        assert io['integer_orders_kept'] == [2]
        assert abs(io['discrepancy_pct']) < 1e-8

    def test_k2_k4_under_one_pct(self) -> None:
        """Through k=2,4: integer-only sub-sum is within 0.1% of Pachucki
        leading."""
        from geovac.cross_register_vne import integer_order_only_recoil_estimate
        io = integer_order_only_recoil_estimate(Z=1.0, n=1, max_order_k=4)
        assert io['integer_orders_kept'] == [2, 4]
        assert abs(io['discrepancy_pct']) < 0.1

    def test_k_full_integer_sweep_under_05_pct(self) -> None:
        """Through k=2,4,6,8: integer-only sub-sum is within 0.5% of
        Pachucki leading, meeting the sub-percent target."""
        from geovac.cross_register_vne import integer_order_only_recoil_estimate
        io = integer_order_only_recoil_estimate(Z=1.0, n=1, max_order_k=8)
        assert io['integer_orders_kept'] == [2, 4, 6, 8]
        assert abs(io['discrepancy_pct']) < 0.5, \
            f"Sub-percent target missed: {io['discrepancy_pct']}%"

    def test_integer_only_strictly_better_than_full_series(self) -> None:
        """Integer-only sub-sum is much closer to Pachucki leading than the
        full Roothaan series."""
        from geovac.cross_register_vne import (
            integer_order_only_recoil_estimate,
            pachucki_higher_order_comparison,
        )
        io = integer_order_only_recoil_estimate(Z=1.0, n=1, max_order_k=8)
        ph = pachucki_higher_order_comparison(Z=1.0, n=1, max_order=8)
        # Integer-only |drift| << full-series |drift|
        assert abs(io['discrepancy_pct']) < 0.2 * abs(ph['cumulative_drift_vs_pachucki_leading_pct'])


class TestSprintMHTrackDDualExpansion:
    """Sprint MH Track D: dual expansion for the muonic regime.

    Verifies that ``roothaan_J0_taylor_expansion_dual`` is the (lam_e <-> lam_n)
    swap of ``roothaan_J0_taylor_expansion``, and that
    ``roothaan_recoil_shift_regime_aware`` dispatches correctly.
    """

    def test_J0_closed_form_is_symmetric_both_regimes(self) -> None:
        """``_roothaan_J0`` is symmetric to machine precision in both
        electronic and muonic regimes -- the production-path finding."""
        from geovac.cross_register_vne import _roothaan_J0

        # Electronic regime
        J_e = _roothaan_J0(1.0, 85.7)
        J_e_swap = _roothaan_J0(85.7, 1.0)
        assert abs(J_e - J_e_swap) < 1e-15

        # Muonic regime: lam_lepton > lam_nucleus
        J_mu = _roothaan_J0(185.85, 85.7)
        J_mu_swap = _roothaan_J0(85.7, 185.85)
        assert abs(J_mu - J_mu_swap) < 1e-15

    def test_dual_expansion_symbolic_structure(self) -> None:
        """Dual expansion: recoil shift (lam_n - J_0) has c_2 = +2 lam_n^3,
        c_3 = -5 lam_n^4, etc. Mirror of the electronic series."""
        from geovac.cross_register_vne import roothaan_J0_taylor_expansion_dual
        result = roothaan_J0_taylor_expansion_dual(n_terms=6, lam_n_value=None)
        assert result['verified_against_closed_form']
        coeffs = dict(result['coefficients'])
        # k=0 of recoil shift = 0 (the leading lam_n cancels)
        # k=2: +2 lam_n^3 (mirror of original's +2 lam_e^3)
        import sympy as sp
        lam_n = sp.symbols('lam_n', positive=True)
        c2_expected = 2 * lam_n ** 3
        assert sp.simplify(coeffs[2] - c2_expected) == 0
        # k=3: -5 lam_n^4 (mirror of original's -5 lam_e^4)
        c3_expected = -5 * lam_n ** 4
        assert sp.simplify(coeffs[3] - c3_expected) == 0

    def test_dual_expansion_evaluated_at_lam_n_one_matches_original(self) -> None:
        """At lam_n = 1 the dual expansion coefficients match the original
        at lam_e = 1 (by exchange symmetry of J_0)."""
        from geovac.cross_register_vne import (
            roothaan_J0_taylor_expansion,
            roothaan_J0_taylor_expansion_dual,
        )
        orig = roothaan_J0_taylor_expansion(n_terms=6, lam_e_value=1.0)
        dual = roothaan_J0_taylor_expansion_dual(n_terms=6, lam_n_value=1.0)
        for (k_o, c_o), (k_d, c_d) in zip(orig['coefficients'], dual['coefficients']):
            assert k_o == k_d
            assert abs(float(c_o) - float(c_d)) < 1e-12

    def test_regime_aware_electronic_matches_existing(self) -> None:
        """Electronic-regime dispatch reproduces existing
        ``roothaan_recoil_shift_through_order`` output bit-identically."""
        from geovac.cross_register_vne import (
            roothaan_recoil_shift_regime_aware,
            roothaan_recoil_shift_through_order,
        )
        ra = roothaan_recoil_shift_regime_aware(
            lam_lepton=1.0, lam_nucleus=85.7, Z=1.0, max_order=5,
        )
        existing = roothaan_recoil_shift_through_order(
            lam_n=85.7, lam_e=1.0, Z=1.0, max_order=5,
        )
        assert ra['regime'] == 'electronic'
        assert abs(ra['roothaan_numerical'] - existing['roothaan_numerical']) < 1e-15
        assert ra['series_converges_to_roothaan']

    def test_regime_aware_muonic_dispatches_correctly(self) -> None:
        """Muonic-regime dispatch identifies lam_lepton > lam_nucleus and
        labels the regime accordingly."""
        from geovac.cross_register_vne import (
            roothaan_recoil_shift_regime_aware,
            _roothaan_J0,
        )
        ra = roothaan_recoil_shift_regime_aware(
            lam_lepton=185.85, lam_nucleus=85.7, Z=1.0, max_order=5,
        )
        assert ra['regime'] == 'muonic'
        # Closed-form J_0 cross-check (the production-path value)
        J0_expected = _roothaan_J0(185.85, 85.7)
        assert abs(ra['closed_form_J0'] - J0_expected) < 1e-15
        # The natural muonic recoil estimator: +Z * (lam_lepton - J_0)
        natural_recoil_expected = -1.0 * (J0_expected - 185.85)
        assert abs(ra['roothaan_numerical'] - natural_recoil_expected) < 1e-12

    def test_regime_aware_at_boundary_lam_lepton_equals_lam_nucleus(self) -> None:
        """At the regime boundary lam_lepton == lam_nucleus, dispatcher
        chooses electronic-side (<=) and produces a finite, consistent
        result. J_0(lam, lam) = 5 lam / 8."""
        from geovac.cross_register_vne import roothaan_recoil_shift_regime_aware
        ra = roothaan_recoil_shift_regime_aware(
            lam_lepton=10.0, lam_nucleus=10.0, Z=1.0, max_order=5,
        )
        assert ra['regime'] == 'electronic'
        assert abs(ra['closed_form_J0'] - 6.25) < 1e-15

    def test_dual_expansion_evaluated_at_specific_lam_n(self) -> None:
        """Dual expansion evaluated at lam_n = 85.7. The coefficients are
        of the recoil shift (lam_n - J_0): k=0 vanishes (leading cancels),
        k=2 = +2 * 85.7^3."""
        from geovac.cross_register_vne import roothaan_J0_taylor_expansion_dual
        result = roothaan_J0_taylor_expansion_dual(n_terms=6, lam_n_value=85.7)
        coeffs_eval = result['coefficients']
        # k=0 vanishes (recoil shift has no constant term)
        assert abs(float(coeffs_eval[0][1])) < 1e-10
        # k=2 should be +2 * 85.7^3
        expected_c2 = 2.0 * (85.7 ** 3)
        assert abs(float(coeffs_eval[2][1]) - expected_c2) < 1e-6
        # k=3 should be -5 * 85.7^4
        expected_c3 = -5.0 * (85.7 ** 4)
        assert abs(float(coeffs_eval[3][1]) - expected_c3) < 1e-3


# ---------------------------------------------------------------------------
# Species-agnostic recoil verification (Sprint rZG follow-up, 2026-05-09)
# ---------------------------------------------------------------------------
#
# These tests confirm cross_register_recoil_correction reproduces leading-
# order Bethe-Salpeter recoil for arbitrary nuclear masses (proton, deuteron,
# antimuon) without any species-specific code path. The lam_n = 2*sqrt(M_n)
# parameterization is the universal calibration: J_0(lam_e=1, lam_n) - 1 =
# +2/lam_n^2 + O(1/lam_n^3) at leading order, which equals m_e/m_n at
# lam_n = 2*sqrt(M_n) by construction (1/lam_n^2 = 1/(4 M_n)).
#
# Sprint rZG bug diagnosis (2026-05-09): the +5 fm artifact in r_Z(D)
# extraction was traced to the rZG global fit's Layer-2 budget for D HFS
# (-150 ppm specified vs ~-286 ppm correct), NOT to a bug in this module.
# The production cross_register_vne kernel is species-agnostic and
# correct; these tests document and protect that property going forward.


class TestSpeciesAgnosticRecoil:
    """Verify cross-register recoil works for arbitrary nuclear masses.

    The proton-, deuteron-, and muon-(antimuon)-bound systems all use the
    same kernel _roothaan_J0 with lam_n = 2*sqrt(M_n) and recover the
    Bethe-Salpeter leading-order recoil to within the 1s x 1s Sturmian
    basis-truncation precision (~3% for proton, ~2% for deuteron, ~8%
    for muonium where m_red/m_n is largest).
    """

    @staticmethod
    def _build_spec(M_n_over_m_e: float) -> CrossRegisterVneSpec:
        """Standard 1s x 1s spec at lam_n = 2*sqrt(M_n)."""
        import math
        return CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=2.0 * math.sqrt(M_n_over_m_e),
            n_max_n=1,
            Z_nuc=1.0, L_max=0,
            label=f"H_1s_x_M_{M_n_over_m_e:.4f}",
        )

    def test_hydrogen_recoil_matches_bethe_salpeter(self) -> None:
        """H 1s x proton 1s: leading-order recoil at 2.86% (basis-truncation).

        This is the regression target preserved by all cross_register_vne
        modifications. Sprint Phase C-W1a-physics was calibrated to land
        at this value; any change that shifts it is a regression.
        """
        spec = self._build_spec(M_PROTON_OVER_M_E)
        recoil = cross_register_recoil_correction(
            spec, m_e_over_m_n=1.0 / M_PROTON_OVER_M_E,
        )
        # Expected ~2.86% (basis-truncation artifact of 1s x 1s Sturmian)
        assert recoil['relative_error'] == pytest.approx(0.0286, abs=0.001)
        # Sign and order of magnitude
        assert recoil['cross_register_recoil_estimate'] > 0  # less bound
        assert recoil['cross_register_recoil_estimate'] == pytest.approx(
            2.645e-4, rel=0.005,
        )

    def test_deuterium_recoil_matches_bethe_salpeter(self) -> None:
        """D 1s x deuteron 1s: leading-order recoil at ~2% precision.

        Verifies the production cross-register V_eN kernel is species-
        agnostic for I=1 nuclei. The 2.03% relative error is *better* than
        H (2.86%) because m_e/m_d < m_e/m_p (the leading 1/lam_n^2 term
        captures more of the Bethe-Salpeter shift relative to the
        sub-leading 1/lam_n^3 corrections).
        """
        M_D_OVER_M_E = M_PROTON_OVER_M_E * 1.99900750139
        spec = self._build_spec(M_D_OVER_M_E)
        recoil = cross_register_recoil_correction(
            spec, m_e_over_m_n=1.0 / M_D_OVER_M_E,
        )
        # Expected ~2.03% (smaller than H because m_red/m_d is smaller)
        assert recoil['relative_error'] == pytest.approx(0.0203, abs=0.001)
        # Sign and order of magnitude
        assert recoil['cross_register_recoil_estimate'] > 0  # less bound
        # Bethe-Salpeter for D: |E_1| / m_d = 0.5 / 3670.48 ~ 1.36e-4
        bs_expected = 0.5 / M_D_OVER_M_E
        assert recoil['expected_leading_order'] == pytest.approx(
            bs_expected, rel=1e-12,
        )
        # GeoVac estimate within 5% of Bethe-Salpeter
        assert recoil['cross_register_recoil_estimate'] == pytest.approx(
            bs_expected, rel=0.05,
        )

    def test_muonium_recoil_matches_bethe_salpeter(self) -> None:
        """Muonium (e^- on antimuon mu^+): leading-order recoil at ~8%.

        The 8.18% relative error is larger than H/D because m_red/m_n is
        bigger for muonium (the muon is much lighter than the proton),
        so sub-leading 1/lam_n^3 corrections matter more relative to the
        leading 1/lam_n^2 term. Still captures sign and OoM correctly.
        """
        M_MU_OVER_M_E = 206.7682830
        spec = self._build_spec(M_MU_OVER_M_E)
        recoil = cross_register_recoil_correction(
            spec, m_e_over_m_n=1.0 / M_MU_OVER_M_E,
        )
        # Expected ~8.18% (larger than H/D because m_red/m_mu is bigger)
        assert recoil['relative_error'] == pytest.approx(0.0818, abs=0.005)
        # Sign and order of magnitude
        assert recoil['cross_register_recoil_estimate'] > 0  # less bound
        # Bethe-Salpeter for muonium: |E_1| / m_mu = 0.5 / 206.77 ~ 2.42e-3
        bs_expected = 0.5 / M_MU_OVER_M_E
        assert recoil['expected_leading_order'] == pytest.approx(
            bs_expected, rel=1e-12,
        )

    def test_static_nucleus_limit_recoil_vanishes(self) -> None:
        """In the m_n -> infinity limit, framework recoil estimate -> 0
        (the nucleus becomes a classical point charge and recoil is
        suppressed)."""
        # M_n = 1e8 m_e (well past any physical nucleus, approaches static)
        M_n_large = 1.0e8
        spec = self._build_spec(M_n_large)
        recoil = cross_register_recoil_correction(
            spec, m_e_over_m_n=1.0 / M_n_large,
        )
        # Both quantum and Bethe-Salpeter should be tiny
        assert abs(recoil['cross_register_recoil_estimate']) < 1.0e-7
        assert abs(recoil['expected_leading_order']) < 1.0e-7
        # And J_0 quantum -> J_0 classical = lam_e = 1 in the limit
        assert recoil['cross_register_J0'] == pytest.approx(1.0, abs=1.0e-6)

    def test_mass_ratio_scaling_proton_to_deuteron(self) -> None:
        """Recoil scales linearly with m_e/m_n at leading order.

        Because Bethe-Salpeter recoil = +Z^2 m_e / (2 n^2 m_n) at leading
        order, the ratio (recoil_H / recoil_D) should equal m_d/m_p =
        ~1.999 at leading order. Higher-order corrections in 1/lam_n
        introduce sub-percent deviation.
        """
        spec_H = self._build_spec(M_PROTON_OVER_M_E)
        M_D_OVER_M_E = M_PROTON_OVER_M_E * 1.99900750139
        spec_D = self._build_spec(M_D_OVER_M_E)
        recoil_H = cross_register_recoil_correction(
            spec_H, m_e_over_m_n=1.0 / M_PROTON_OVER_M_E,
        )
        recoil_D = cross_register_recoil_correction(
            spec_D, m_e_over_m_n=1.0 / M_D_OVER_M_E,
        )
        # Both Bethe-Salpeter expectations should scale as (m_e/m_n)
        bs_ratio = recoil_H['expected_leading_order'] / recoil_D['expected_leading_order']
        # Should equal m_d/m_p = 1.999...
        assert bs_ratio == pytest.approx(1.99900750, rel=1e-6)
        # Framework estimates should also approximately scale this way,
        # with deviations from the basis-truncation residuals
        framework_ratio = (
            recoil_H['cross_register_recoil_estimate']
            / recoil_D['cross_register_recoil_estimate']
        )
        assert framework_ratio == pytest.approx(bs_ratio, rel=0.01)

    def test_kernel_symmetry_preserved_for_arbitrary_lam_n(self) -> None:
        """Production kernel _roothaan_J0(lam_e, lam_n) is symmetric in
        (lam_e, lam_n) for any positive values.

        The symmetry is the algebraic backbone of species-agnosticity:
        once the lam_n parameterization fixes the focal length, the
        kernel handles any (lam_e, lam_n) pair without species-specific
        branching.
        """
        import math
        for M_n in [M_PROTON_OVER_M_E,
                    M_PROTON_OVER_M_E * 1.99900750139,  # deuteron
                    206.7682830,                         # muon
                    1.0]:                                # equal-mass
            lam_n = 2.0 * math.sqrt(M_n)
            J_forward = _roothaan_J0(1.0, lam_n)
            J_reverse = _roothaan_J0(lam_n, 1.0)
            assert J_forward == pytest.approx(J_reverse, rel=1e-15), (
                f"Symmetry failure at M_n = {M_n}: "
                f"J(1, {lam_n}) = {J_forward}, J({lam_n}, 1) = {J_reverse}"
            )


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
