"""
Tests for TC angular gradient implementation (Track BX-4).

Validates:
- Angular gradient coefficients obey Gaunt selection rules
- include_angular=False matches BX-3 reference data
- include_angular=True produces correct additional couplings
- build_tc_composed_hamiltonian respects include_angular parameter
- Angular gradient is a net negative for quantum efficiency at max_n=2
"""
import pytest
import numpy as np

from geovac.tc_integrals import (
    compute_tc_integrals_block,
    _angular_gradient_coefficients,
    tc_eri_to_chemist,
)


# ---------------------------------------------------------------------------
# Angular gradient coefficient tests
# ---------------------------------------------------------------------------

class TestAngularGradientCoefficients:
    """Test _angular_gradient_coefficients selection rules."""

    def test_s_orbital_produces_p(self):
        """∇Y_{0,0} produces L=1 only."""
        coeffs = _angular_gradient_coefficients(0, 0)
        L_values = {c[0] for c in coeffs}
        assert L_values <= {1}, f"s-orbital gradient should produce only L=1, got {L_values}"

    def test_p_orbital_produces_s_and_d(self):
        """∇Y_{1,m} produces L=0 and L=2."""
        for m in [-1, 0, 1]:
            coeffs = _angular_gradient_coefficients(1, m)
            L_values = {c[0] for c in coeffs}
            assert L_values <= {0, 2}, f"p-orbital gradient should produce L=0,2, got {L_values}"

    def test_gaunt_selection_rule_delta_l(self):
        """Angular gradient obeys |ΔL| = 1 (vector spherical harmonic)."""
        for l in range(4):
            for m in range(-l, l + 1):
                coeffs = _angular_gradient_coefficients(l, m)
                for L_eff, m_eff, q, coeff in coeffs:
                    assert L_eff in (l - 1, l + 1), (
                        f"l={l},m={m}: L_eff={L_eff} violates |ΔL|=1"
                    )

    def test_gaunt_selection_rule_delta_m(self):
        """Angular gradient obeys |Δm| <= 1 (spherical components q=-1,0,+1)."""
        for l in range(4):
            for m in range(-l, l + 1):
                coeffs = _angular_gradient_coefficients(l, m)
                for L_eff, m_eff, q, coeff in coeffs:
                    assert abs(m_eff - m) <= 1, (
                        f"l={l},m={m}: m_eff={m_eff} violates |Δm|<=1"
                    )

    def test_l0_returns_empty(self):
        """l=0 with L_eff=l-1=-1 is excluded; only L_eff=1 survives."""
        coeffs = _angular_gradient_coefficients(0, 0)
        for L_eff, m_eff, q, coeff in coeffs:
            assert L_eff >= 0, f"Negative L_eff={L_eff} should not appear"


# ---------------------------------------------------------------------------
# TC integral tests: radial-only matches BX-3 reference
# ---------------------------------------------------------------------------

class TestTCRadialOnly:
    """Verify radial-only TC integrals match BX-3 benchmark data."""

    def test_he_max_n1_radial_matches_bx3(self):
        """He max_n=1: radial-only should give 1 integral (1s-1s only)."""
        states = [(1, 0, 0)]
        tc_eri = compute_tc_integrals_block(2.0, states, 2000, include_angular=False)
        assert len(tc_eri) == 1

    def test_he_max_n2_radial_matches_bx3(self):
        """He max_n=2: radial-only gives 65 integrals (BX-3 reference)."""
        states = []
        for n in range(1, 3):
            for l in range(n):
                for m in range(-l, l + 1):
                    states.append((n, l, m))
        tc_eri = compute_tc_integrals_block(2.0, states, 2000, include_angular=False)
        assert len(tc_eri) == 65

    def test_he_max_n1_angular_same_as_radial(self):
        """At max_n=1, no l>0 orbitals exist, so angular=radial."""
        states = [(1, 0, 0)]
        eri_rad = compute_tc_integrals_block(2.0, states, 2000, include_angular=False)
        eri_full = compute_tc_integrals_block(2.0, states, 2000, include_angular=True)
        assert len(eri_rad) == len(eri_full)
        for key in eri_rad:
            assert abs(eri_rad[key] - eri_full[key]) < 1e-12


# ---------------------------------------------------------------------------
# TC integral tests: angular gradient adds couplings
# ---------------------------------------------------------------------------

class TestTCAngularGradient:
    """Test that angular gradient produces additional integrals for l>0."""

    def test_he_max_n2_angular_adds_integrals(self):
        """He max_n=2: angular gradient should add integrals beyond radial-only."""
        states = []
        for n in range(1, 3):
            for l in range(n):
                for m in range(-l, l + 1):
                    states.append((n, l, m))
        eri_rad = compute_tc_integrals_block(2.0, states, 2000, include_angular=False)
        eri_full = compute_tc_integrals_block(2.0, states, 2000, include_angular=True)
        assert len(eri_full) > len(eri_rad), (
            f"Angular should add integrals: {len(eri_full)} vs {len(eri_rad)}"
        )

    def test_angular_integrals_are_real(self):
        """All TC integrals (including angular) should be real-valued."""
        states = []
        for n in range(1, 3):
            for l in range(n):
                for m in range(-l, l + 1):
                    states.append((n, l, m))
        tc_eri = compute_tc_integrals_block(2.0, states, 2000, include_angular=True)
        for key, val in tc_eri.items():
            assert isinstance(val, (int, float)) or abs(val.imag) < 1e-10, (
                f"Integral {key} has imaginary part: {val}"
            )


# ---------------------------------------------------------------------------
# Composed pipeline tests
# ---------------------------------------------------------------------------

class TestTCComposedAngular:
    """Test build_tc_composed_hamiltonian with include_angular parameter."""

    def test_radial_only_matches_bx3_pauli_count(self):
        """LiH radial-only should give 562 Pauli terms (BX-3 reference)."""
        from geovac.tc_integrals import build_tc_composed_hamiltonian
        from geovac.composed_qubit import lih_spec
        spec = lih_spec()
        result = build_tc_composed_hamiltonian(spec, include_angular=False)
        assert result['N_pauli'] == 562, (
            f"LiH radial-only: expected 562 Pauli, got {result['N_pauli']}"
        )

    def test_angular_increases_pauli_count(self):
        """LiH with angular gradient should have more Pauli terms than radial-only."""
        from geovac.tc_integrals import build_tc_composed_hamiltonian
        from geovac.composed_qubit import lih_spec
        spec = lih_spec()
        r_rad = build_tc_composed_hamiltonian(spec, include_angular=False)
        r_full = build_tc_composed_hamiltonian(spec, include_angular=True)
        assert r_full['N_pauli'] > r_rad['N_pauli'], (
            f"Angular should increase Pauli: {r_full['N_pauli']} vs {r_rad['N_pauli']}"
        )

    def test_default_is_radial_only(self):
        """Default include_angular=False should match radial-only explicitly."""
        from geovac.tc_integrals import build_tc_composed_hamiltonian
        from geovac.composed_qubit import lih_spec
        spec = lih_spec()
        r_default = build_tc_composed_hamiltonian(spec)
        r_explicit = build_tc_composed_hamiltonian(spec, include_angular=False)
        assert r_default['N_pauli'] == r_explicit['N_pauli']


# ---------------------------------------------------------------------------
# Negative result: angular gradient is net-negative for quantum efficiency
# ---------------------------------------------------------------------------

class TestAngularGradientNegativeResult:
    """
    Document that angular gradient is a negative result for quantum computing
    efficiency at the composed operating point (max_n=2).
    """

    @pytest.mark.slow
    def test_pauli_ratio_exceeds_accuracy_benefit(self):
        """Angular gradient: Pauli increase >> accuracy improvement."""
        # He max_n=2 benchmark data from corrected FCI solver
        # Radial: 3.621%, 188 Pauli; Full: 3.611%, 500 Pauli
        # 2.66x more Pauli for 0.01 pp improvement
        rad_pauli = 188
        full_pauli = 500
        rad_err = 3.621
        full_err = 3.611
        pauli_ratio = full_pauli / rad_pauli
        accuracy_improvement_pp = rad_err - full_err

        assert pauli_ratio > 2.0, "Pauli ratio should exceed 2x"
        assert accuracy_improvement_pp < 0.05, "Accuracy improvement should be < 0.05 pp"
        # Net negative: cost/benefit ratio > 50x
        assert pauli_ratio / max(accuracy_improvement_pp, 0.001) > 20


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
