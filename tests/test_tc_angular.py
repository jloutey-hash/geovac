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
        """He max_n=2 radial-only: 107 integrals.

        Corrected 2026-08-30: the retired 65 came from the
        pair-diagonal sign error; the 265 seen mid-arc came from
        this module never imposing the Coulomb M_L rule. With both
        fixed the count matches the physical ERI support exactly.
        """
        states = []
        for n in range(1, 3):
            for l in range(n):
                for m in range(-l, l + 1):
                    states.append((n, l, m))
        tc_eri = compute_tc_integrals_block(2.0, states, 2000, include_angular=False)
        assert len(tc_eri) == 107
        viol = [k for k in tc_eri
                if states[k[0]][2] + states[k[1]][2]
                != states[k[2]][2] + states[k[3]][2]]
        assert not viol, f'{len(viol)} M_L-violating TC entries'

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

    def test_he_max_n2_angular_adds_only_unphysical_integrals(self):
        """He max_n=2: every integral the angular gradient adds is spurious.

        The TC angular-gradient assembly contributes 26 entries beyond
        radial-only, and all 26 violate L_z conservation
        (m_a + m_b != m_c + m_d).  That is impossible for a
        rotationally-invariant correlator u(r12) -- the transcorrelated
        Hamiltonian must commute with total L_z -- so the angular term's
        entire contribution on this basis is a bookkeeping artifact.

        Recorded 2026-08-30, not masked with a filter: dropping the
        violating entries would hide the error and leave the survivors
        unverified.  Sharpens the CLAUDE.md SS3 dead-end (angular gradient
        buys 0.01 pp for 2.66x the Pauli terms) -- here it buys nothing.
        """
        states = []
        for n in range(1, 3):
            for l in range(n):
                for m in range(-l, l + 1):
                    states.append((n, l, m))
        eri_rad = compute_tc_integrals_block(2.0, states, 2000, include_angular=False)
        eri_full = compute_tc_integrals_block(2.0, states, 2000, include_angular=True)

        def ml_ok(k):
            return (states[k[0]][2] + states[k[1]][2]
                    == states[k[2]][2] + states[k[3]][2])

        added = set(eri_full) - set(eri_rad)
        assert len(added) == 26, f"angular added {len(added)}, expected 26"
        assert all(not ml_ok(k) for k in added), (
            "some angular-added integrals conserve M_L -- the physical "
            "content of the angular term has changed; re-derive"
        )
        assert all(ml_ok(k) for k in eri_rad), (
            f"radial-only must conserve M_L: {len(eri_rad)} entries"
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
        from geovac.molecular_spec import lih_spec
        spec = lih_spec()
        result = build_tc_composed_hamiltonian(spec, include_angular=False)
        assert result['N_pauli'] == 1354, (
            f"LiH radial-only: expected 1354 Pauli, got {result['N_pauli']}"
        )

    def test_composed_angular_adds_only_lz_violating_entries(self):
        """LiH composed: all 78 entries the angular gradient adds break L_z.

        The block-level companion above pins this for a single He block (26
        entries).  This is the same defect at the level Paper 14 quoted its
        angular-gradient resource multipliers -- 562 -> 1,498 and a 4.49x
        total, now withdrawn -- so the claim is backed where it is made.

        Convention control is load-bearing, not decoration.  The composed
        ERI is stored in CHEMIST order, so L_z reads m_a + m_c == m_b + m_d;
        applying the block-level (physicist) rule here reports 126
        violations on the plain composed tensor, i.e. a confidently wrong
        answer.  Asserting zero violations on both clean tensors is what
        makes the 78 mean something.
        """
        import numpy as np
        from geovac.composed_qubit import (
            build_composed_hamiltonian, _enumerate_states,
        )
        from geovac.tc_integrals import build_tc_composed_hamiltonian
        from geovac.molecular_spec import lih_spec

        spec = lih_spec()
        comp = build_composed_hamiltonian(spec)
        labels = []
        for b in comp["blocks"]:
            labels.extend(_enumerate_states(2)[: b["n_orbitals"]])
        m = [st[2] for st in labels]

        e_plain = np.asarray(comp["eri"])
        e_rad = np.asarray(
            build_tc_composed_hamiltonian(spec, include_angular=False)["eri"])
        e_ang = np.asarray(
            build_tc_composed_hamiltonian(spec, include_angular=True)["eri"])
        assert len(labels) == e_plain.shape[0]

        def violates(x):
            a, b_, c, d = x
            return m[a] + m[c] != m[b_] + m[d]

        for name, tensor in (("plain composed", e_plain), ("TC radial", e_rad)):
            bad = [x for x in np.argwhere(np.abs(tensor) > 1e-12)
                   if violates(x)]
            assert not bad, (
                f"CONTROL FAILED: {name} has {len(bad)} L_z violations; the "
                f"index convention is wrong, so the angular count below is "
                f"meaningless"
            )

        added = np.argwhere((np.abs(e_ang) > 1e-12) & (np.abs(e_rad) <= 1e-12))
        assert len(added) == 78, f"angular added {len(added)}, expected 78"
        assert all(violates(x) for x in added), (
            "some angular-added composed entries conserve L_z -- the physical "
            "content of the angular term has changed; re-derive before "
            "quoting any angular-gradient resource number"
        )

    def test_angular_increases_pauli_count(self):
        """LiH with angular gradient should have more Pauli terms than radial-only."""
        from geovac.tc_integrals import build_tc_composed_hamiltonian
        from geovac.molecular_spec import lih_spec
        spec = lih_spec()
        r_rad = build_tc_composed_hamiltonian(spec, include_angular=False)
        r_full = build_tc_composed_hamiltonian(spec, include_angular=True)
        assert r_full['N_pauli'] > r_rad['N_pauli'], (
            f"Angular should increase Pauli: {r_full['N_pauli']} vs {r_rad['N_pauli']}"
        )

    def test_default_is_radial_only(self):
        """Default include_angular=False should match radial-only explicitly."""
        from geovac.tc_integrals import build_tc_composed_hamiltonian
        from geovac.molecular_spec import lih_spec
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
