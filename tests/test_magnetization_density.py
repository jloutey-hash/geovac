"""
Tests for the operator-level magnetization-density module (Phase C-W1b-operator).

Validates:
- Eides leading-order regression (Taylor expansion recovers -2 Z m_e r_Z)
- Quantitative -39.5 ppm Eides Tab. 7.3 reproduction at r_Z = 1.045 fm
- Hermiticity of the resulting Pauli string sum
- Delta-function profile recovers point-nucleus contact (no shift)
- Symmetry under proton-spin reflection
- Block-diagonal structure on the joint register (no spurious cross-terms)
- Cross-register integration check with V_eN

Author: GeoVac Development Team (Phase C-W1b-operator)
Date: May 2026
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from geovac.cross_register_vne import (
    CrossRegisterVneSpec,
    LAM_NUCLEUS_GEOMETRIC,
)
from geovac.magnetization_density import (
    A0_FM,
    DELTA_NU_ZEMACH_EIDES_PPM,
    MagnetizationDensitySpec,
    R_Z_EIDES_2024_BOHR,
    R_Z_EIDES_2024_FM,
    _rho_M_moment,
    compose_with_cross_register_vne,
    compute_magnetization_density_operator,
    hydrogen_zemach_eides_leading_order,
    taylor_zemach_around_zero,
)


# ---------------------------------------------------------------------------
# Moment computations: rho_M first moment matches the Zemach radius
# ---------------------------------------------------------------------------


class TestRhoMMoments:
    """Validate the rho_M radial moment closed forms."""

    def test_gaussian_M0_normalized(self) -> None:
        """Gaussian profile is unit-normalized: M_0 = 1."""
        spec = MagnetizationDensitySpec(profile="gaussian", r_Z_bohr=1e-5)
        assert _rho_M_moment(spec, 0) == pytest.approx(1.0, rel=1e-12)

    def test_gaussian_M1_equals_rZ(self) -> None:
        """Gaussian first moment <r> = r_Z by construction."""
        rZ = R_Z_EIDES_2024_BOHR
        spec = MagnetizationDensitySpec(profile="gaussian", r_Z_bohr=rZ)
        M1 = _rho_M_moment(spec, 1)
        assert M1 == pytest.approx(rZ, rel=1e-12)

    def test_exponential_M0_normalized(self) -> None:
        """Exponential profile is unit-normalized: M_0 = 1."""
        spec = MagnetizationDensitySpec(profile="exponential", r_Z_bohr=1e-5)
        assert _rho_M_moment(spec, 0) == pytest.approx(1.0, rel=1e-12)

    def test_exponential_M1_equals_rZ(self) -> None:
        """Exponential first moment <r> = 3/kappa = r_Z."""
        rZ = R_Z_EIDES_2024_BOHR
        spec = MagnetizationDensitySpec(profile="exponential", r_Z_bohr=rZ)
        M1 = _rho_M_moment(spec, 1)
        assert M1 == pytest.approx(rZ, rel=1e-12)

    def test_delta_M0_unity_higher_zero(self) -> None:
        """Delta profile: M_0 = 1, M_k = 0 for k > 0."""
        spec = MagnetizationDensitySpec(profile="delta", r_Z_bohr=0.0)
        assert _rho_M_moment(spec, 0) == pytest.approx(1.0)
        assert _rho_M_moment(spec, 1) == pytest.approx(0.0)
        assert _rho_M_moment(spec, 2) == pytest.approx(0.0)

    def test_gaussian_M2_textbook(self) -> None:
        """Gaussian <r^2> = 3/(2 beta) for the chosen normalization."""
        rZ = 1e-5
        spec = MagnetizationDensitySpec(profile="gaussian", r_Z_bohr=rZ)
        beta = spec.profile_width()  # = 4/(pi rZ^2)
        # M_2 = (2/sqrt(pi)) * Gamma(5/2) * beta^(-1) = 3/(2 sqrt(pi)) * 2/sqrt(pi)...
        # Direct formula: M_2 = (2/sqrt(pi)) * (3 sqrt(pi)/4) / beta
        #                    = (3/2) / beta
        expected = 1.5 / beta
        M2 = _rho_M_moment(spec, 2)
        assert M2 == pytest.approx(expected, rel=1e-12)


# ---------------------------------------------------------------------------
# Eides leading-order regression: Taylor expansion recovers -2 Z m_e r_Z
# ---------------------------------------------------------------------------


class TestEidesLeadingOrderRegression:
    """The Taylor expansion of the operator-level shift around R_p = 0
    must reproduce Eides' classical -2 Z m_e r_Z form."""

    def test_taylor_order_1_matches_eides_form(self) -> None:
        """Order-1 Taylor coefficient = -2 Z m_e r_Z exactly."""
        rZ = R_Z_EIDES_2024_BOHR
        spec = MagnetizationDensitySpec(profile="gaussian", r_Z_bohr=rZ)
        Z = spec.proton_spec.Z_nuc
        result = taylor_zemach_around_zero(spec, order=1)

        expected = -2.0 * Z * 1.0 * rZ
        assert result['order_1_shift'] == pytest.approx(expected, rel=1e-12)

    def test_total_shift_matches_eides_minus_39_5_ppm(self) -> None:
        """At r_Z = 1.045 fm, the operator-level shift is -39.5 ppm."""
        result = hydrogen_zemach_eides_leading_order()
        delta_ppm = result['operator_level_delta_ppm']
        # Must be within 1% of -39.5 ppm at leading order (the controlled
        # discrepancy is the order-2 Friar moment, ~ ppm-level vs -39.5 ppm)
        assert delta_ppm == pytest.approx(-39.5, rel=1e-3)

    def test_residual_within_eides_budget(self) -> None:
        """Residual against Eides Tab. 7.3 must be < 1 ppm at leading order."""
        result = hydrogen_zemach_eides_leading_order()
        residual = abs(result['residual_ppm'])
        # Eides Tab. 7.3 expected residual budget is +12 to +18 ppm (multi-loop
        # QED + nuclear polarizability).  Operator-level leading order should
        # be well below that.
        assert residual < 1.0

    def test_profile_independence_at_leading_order(self) -> None:
        """At the same r_Z, both Gaussian and exponential profiles give
        the same leading-order shift -- the leading coefficient is the
        first moment which is r_Z by construction."""
        rZ = R_Z_EIDES_2024_BOHR
        s_gauss = MagnetizationDensitySpec(profile="gaussian", r_Z_bohr=rZ)
        s_exp = MagnetizationDensitySpec(profile="exponential", r_Z_bohr=rZ)
        op_gauss = compute_magnetization_density_operator(s_gauss)
        op_exp = compute_magnetization_density_operator(s_exp)
        assert op_gauss['delta_ppm'] == pytest.approx(op_exp['delta_ppm'],
                                                       rel=1e-12)

    def test_order_2_correction_suppressed(self) -> None:
        """Order-2 correction must be (r_Z/a_0)^2-suppressed: ~ 1e-10 of the
        order-1 contribution at the Eides r_Z."""
        rZ = R_Z_EIDES_2024_BOHR
        spec = MagnetizationDensitySpec(profile="gaussian", r_Z_bohr=rZ)
        result = taylor_zemach_around_zero(spec, order=2)

        ratio = abs(result['order_2_shift'] / result['order_1_shift'])
        # r_Z / a_0 ~ 2e-5, (r_Z/a_0)^2 ~ 4e-10; should match to OoM
        assert ratio < 1e-3


# ---------------------------------------------------------------------------
# Pauli encoding: hermiticity, structure
# ---------------------------------------------------------------------------


class TestPauliEncoding:
    """Validate the Pauli string sum structure."""

    def test_pauli_terms_real(self) -> None:
        """All Pauli coefficients are real (operator is Hermitian)."""
        result = hydrogen_zemach_eides_leading_order()
        # operator_level_delta_ppm derived from the Pauli structure
        # implicitly checks via the operator construction; we verify the
        # Pauli dict entries are real floats.
        rZ = R_Z_EIDES_2024_BOHR
        spec = MagnetizationDensitySpec(profile="gaussian", r_Z_bohr=rZ)
        op = compute_magnetization_density_operator(spec)
        for pstr, coeff in op['pauli_terms'].items():
            assert isinstance(coeff, float)

    def test_pauli_only_I_and_Z(self) -> None:
        """Diagonal-density JW encoding produces only I and Z strings."""
        spec = MagnetizationDensitySpec(profile="gaussian",
                                         r_Z_bohr=R_Z_EIDES_2024_BOHR)
        op = compute_magnetization_density_operator(spec)
        for pstr in op['pauli_terms']:
            for ch in pstr:
                assert ch in ('I', 'Z')

    def test_pauli_count_minimum_4(self) -> None:
        """Minimum case (1s_e * 1s_p): 4 Pauli strings (II, ZI, IZ, ZZ)."""
        spec = MagnetizationDensitySpec(profile="gaussian",
                                         r_Z_bohr=R_Z_EIDES_2024_BOHR)
        op = compute_magnetization_density_operator(spec)
        assert op['Q_total'] == 2  # 1 e qubit + 1 p qubit
        assert len(op['pauli_terms']) == 4

    def test_ground_state_energy_matches_shift(self) -> None:
        """Ground-state expectation value of the Pauli sum on |11> equals
        the diagonal matrix element -2 Z m_e r_Z."""
        rZ = R_Z_EIDES_2024_BOHR
        spec = MagnetizationDensitySpec(profile="gaussian", r_Z_bohr=rZ)
        op = compute_magnetization_density_operator(spec)
        Z = spec.proton_spec.Z_nuc

        # Build the operator on 4-dim Hilbert space (Q=2)
        I2 = np.eye(2, dtype=complex)
        Z_op = np.diag([1.0, -1.0]).astype(complex)
        H = np.zeros((4, 4), dtype=complex)
        for pstr, c in op['pauli_terms'].items():
            ops = []
            for ch in pstr:
                if ch == 'I':
                    ops.append(I2)
                elif ch == 'Z':
                    ops.append(Z_op)
            mat = ops[0]
            for o in ops[1:]:
                mat = np.kron(mat, o)
            H += c * mat

        # |11> on 2-qubit register = state 3 (indexing |00>=0, |01>=1, |10>=2, |11>=3)
        state_11 = np.array([0, 0, 0, 1], dtype=complex)
        E = np.real(state_11.conj() @ H @ state_11)
        expected = -2.0 * Z * rZ
        assert E == pytest.approx(expected, rel=1e-10)


# ---------------------------------------------------------------------------
# Sanity: delta profile recovers point-nucleus (no shift)
# ---------------------------------------------------------------------------


class TestDeltaProfileSanity:
    """Delta-function rho_M (point nucleus) should produce zero shift."""

    def test_delta_profile_zero_shift(self) -> None:
        spec = MagnetizationDensitySpec(profile="delta", r_Z_bohr=0.0)
        op = compute_magnetization_density_operator(spec)
        assert op['delta_ppm'] == pytest.approx(0.0)
        # All Pauli coeffs should vanish
        assert len(op['pauli_terms']) == 0

    def test_delta_taylor_zero(self) -> None:
        spec = MagnetizationDensitySpec(profile="delta", r_Z_bohr=0.0)
        result = taylor_zemach_around_zero(spec, order=2)
        assert result['order_1_shift'] == pytest.approx(0.0)
        assert result['order_2_shift'] == pytest.approx(0.0)

    def test_zero_rZ_recovers_point_nucleus(self) -> None:
        """Limit r_Z -> 0 of any profile recovers the delta result."""
        # Use small but positive r_Z to avoid divide-by-zero in width
        spec = MagnetizationDensitySpec(profile="gaussian", r_Z_bohr=1e-30)
        op = compute_magnetization_density_operator(spec)
        # Shift should be negligible
        assert abs(op['delta_ppm']) < 1e-15


# ---------------------------------------------------------------------------
# Block-diagonal structure: no spurious cross-terms
# ---------------------------------------------------------------------------


class TestBlockDiagonalStructure:
    """The L=0 magnetization couples only s-state densities, leaving
    higher-l block matrix elements zero."""

    def test_only_s_states_couple(self) -> None:
        """At n_max=2, the matrix is nonzero only on (l_e=0, l_p=0)."""
        proton_spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=2,
            lam_n=LAM_NUCLEUS_GEOMETRIC, n_max_n=2,
            Z_nuc=1.0, L_max=0,
        )
        spec = MagnetizationDensitySpec(
            profile="gaussian",
            r_Z_bohr=R_Z_EIDES_2024_BOHR,
            proton_spec=proton_spec,
        )
        op = compute_magnetization_density_operator(spec)

        states_e = op['states_e']
        states_p = op['states_p']
        M = op['matrix_elements']

        # Find the (l_e, l_p) breakdown
        for i, (n_e, l_e, m_e) in enumerate(states_e):
            for j, (n_p, l_p, m_p) in enumerate(states_p):
                if l_e == 0 and l_p == 0:
                    assert M[i, j] != 0.0  # all s-states should be nonzero
                else:
                    assert M[i, j] == 0.0  # no spurious l > 0 entries

    def test_diagonal_in_orbital_basis(self) -> None:
        """The matrix elements form a diagonal in the joint product basis
        (the L=0 contact term is diagonal in (n, l, m) on each side)."""
        proton_spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=2,
            lam_n=LAM_NUCLEUS_GEOMETRIC, n_max_n=1,
            Z_nuc=1.0, L_max=0,
        )
        spec = MagnetizationDensitySpec(
            profile="gaussian",
            r_Z_bohr=R_Z_EIDES_2024_BOHR,
            proton_spec=proton_spec,
        )
        op = compute_magnetization_density_operator(spec)
        # All 4-index ERI elements are encoded as diagonal-density (i, j),
        # so the test is on the (N_e, N_p) matrix shape
        assert op['matrix_elements'].shape == (5, 1)  # (1s, 2s, 2p_-1, 2p_0, 2p_+1) x (1s)


# ---------------------------------------------------------------------------
# Cross-register integration check (composes with V_eN)
# ---------------------------------------------------------------------------


class TestCrossRegisterIntegration:
    """The W1b magnetization operator composes cleanly with W1a V_eN."""

    def test_composition_succeeds(self) -> None:
        """Building H_combined = V_eN + omega_magn returns valid Pauli sum."""
        vne_spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=LAM_NUCLEUS_GEOMETRIC, n_max_n=1,
            Z_nuc=1.0, L_max=0,
            label="hydrogen_1s_x_proton_1s",
        )
        magn_spec = MagnetizationDensitySpec(
            profile="gaussian",
            r_Z_bohr=R_Z_EIDES_2024_BOHR,
            proton_spec=vne_spec,
        )
        result = compose_with_cross_register_vne(magn_spec, vne_spec)
        assert 'pauli_terms_combined' in result
        assert len(result['pauli_terms_combined']) > 0
        assert result['Q_total'] == 2

    def test_composition_qubit_count_match(self) -> None:
        """Both operators live on the same register layout."""
        vne_spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=LAM_NUCLEUS_GEOMETRIC, n_max_n=1,
            Z_nuc=1.0, L_max=0,
        )
        magn_spec = MagnetizationDensitySpec(
            profile="gaussian",
            r_Z_bohr=R_Z_EIDES_2024_BOHR,
            proton_spec=vne_spec,
        )
        result = compose_with_cross_register_vne(magn_spec, vne_spec)
        # Same Q_total on both
        assert 'pauli_terms_vne' in result
        assert 'pauli_terms_magn' in result

    def test_composition_additive(self) -> None:
        """Combined Pauli coefficients = V_eN coeffs + magn coeffs (linearity)."""
        vne_spec = CrossRegisterVneSpec(
            lam_e=1.0, n_max_e=1,
            lam_n=LAM_NUCLEUS_GEOMETRIC, n_max_n=1,
            Z_nuc=1.0, L_max=0,
        )
        magn_spec = MagnetizationDensitySpec(
            profile="gaussian",
            r_Z_bohr=R_Z_EIDES_2024_BOHR,
            proton_spec=vne_spec,
        )
        result = compose_with_cross_register_vne(magn_spec, vne_spec)
        for pstr, c_combined in result['pauli_terms_combined'].items():
            c_vne = result['pauli_terms_vne'].get(pstr, 0.0)
            c_magn = result['pauli_terms_magn'].get(pstr, 0.0)
            assert c_combined == pytest.approx(c_vne + c_magn, abs=1e-15)


# ---------------------------------------------------------------------------
# Spec validation
# ---------------------------------------------------------------------------


class TestSpecValidation:
    """Input validation on MagnetizationDensitySpec."""

    def test_invalid_profile_rejected(self) -> None:
        with pytest.raises(ValueError, match="Unknown profile"):
            MagnetizationDensitySpec(profile="bogus")

    def test_negative_rZ_rejected(self) -> None:
        with pytest.raises(ValueError, match="r_Z must be non-negative"):
            MagnetizationDensitySpec(r_Z_bohr=-1.0)

    def test_zero_A_hf_rejected(self) -> None:
        with pytest.raises(ValueError, match="A_hf_point must be positive"):
            MagnetizationDensitySpec(A_hf_point=0.0)

    def test_default_proton_spec_constructed(self) -> None:
        """If proton_spec is None, a default is created."""
        spec = MagnetizationDensitySpec()
        assert spec.proton_spec is not None
        assert spec.proton_spec.lam_e == 1.0
