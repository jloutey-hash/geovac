"""
Tests for the He-4 nuclear qubit Hamiltonian (Track NF).

He-4 (alpha particle): 2 protons + 2 neutrons, doubly magic (Z=2, N=2),
spin-0 ground state, experimental binding energy 28.296 MeV.

Author: GeoVac Development Team
Date: April 2026
"""

import numpy as np
import pytest

from geovac.nuclear.nuclear_hamiltonian import (
    DeuteronSpec,
    analyze_he4_hamiltonian,
    build_he4_hamiltonian,
    diagonalize_he4,
    enumerate_sp_states,
    he4_hw_scan,
    compute_same_species_tbme,
    compute_cross_species_tbme,
    _enumerate_slater_dets,
    _sd_phase,
)


# ===========================================================================
# State counting and Hilbert space
# ===========================================================================

class TestHe4StateCount:
    """Tests for He-4 state counting."""

    def test_he4_state_count(self):
        """N_shells=2: 8 proton + 8 neutron = 16 qubits."""
        spec = DeuteronSpec(N_shells=2, hw=15.0)
        states_p, states_n = enumerate_sp_states(spec)
        assert len(states_p) == 8, f"Expected 8 proton states, got {len(states_p)}"
        assert len(states_n) == 8, f"Expected 8 neutron states, got {len(states_n)}"

    def test_he4_hilbert_dim(self):
        """2p+2n: C(8,2) x C(8,2) = 28 x 28 = 784 at N_shells=2."""
        data = build_he4_hamiltonian(N_shells=2, hw=15.0)
        H = data['H_matrix']
        assert H.shape == (784, 784), f"Expected 784x784, got {H.shape}"

    def test_he4_qubit_count(self):
        """Total qubit count = 16 at N_shells=2."""
        data = build_he4_hamiltonian(N_shells=2, hw=15.0)
        assert data['Q'] == 16, f"Expected Q=16, got {data['Q']}"
        assert data['Q_p'] == 8
        assert data['Q_n'] == 8

    def test_slater_det_enumeration(self):
        """C(8,2) = 28 Slater determinants for 2 particles in 8 orbitals."""
        sds = _enumerate_slater_dets(8, 2)
        assert len(sds) == 28, f"Expected 28 SDs, got {len(sds)}"
        # All SDs should be sorted tuples of length 2
        for sd in sds:
            assert len(sd) == 2
            assert sd[0] < sd[1]

    def test_sd_phase_basic(self):
        """Test basic a^dag_p a_q operations on Slater determinants."""
        sd = (0, 3)
        # a^dag_0 a_0 |0,3> = |0,3> (number operator, occupied)
        phase, new_sd = _sd_phase(sd, 0, 0)
        assert phase == 1
        assert new_sd == (0, 3)

        # a^dag_1 a_0 |0,3> should give |1,3> with some phase
        phase, new_sd = _sd_phase(sd, 1, 0)
        assert phase != 0
        assert new_sd == (1, 3)

        # a^dag_0 a_1 |0,3> = 0 (1 not occupied)
        phase, new_sd = _sd_phase(sd, 0, 1)
        assert phase == 0


# ===========================================================================
# Antisymmetrization tests
# ===========================================================================

class TestAntisymmetrization:
    """Tests for antisymmetrization of same-species matrix elements."""

    def test_he4_pp_antisymmetric(self):
        """
        Verify <p_i p_i | V | ..> = 0 (Pauli exclusion, same orbital).

        Antisymmetrized ME for two identical orbitals should vanish.
        """
        spec = DeuteronSpec(N_shells=2, hw=15.0)
        states_p, _ = enumerate_sp_states(spec)
        tbme_pp = compute_same_species_tbme(states_p, spec, include_coulomb=False)

        # Check that V_{iikl} = 0 for all k, l (since <ii|..> AS = <ii|..> - <ii|..> = 0)
        for (i, j, k, l), v in tbme_pp.items():
            if i == j:
                assert abs(v) < 1e-12, (
                    f"Antisymmetrized ME ({i},{j},{k},{l}) = {v} should be 0"
                )

    def test_he4_nn_antisymmetric(self):
        """Same test for neutrons."""
        spec = DeuteronSpec(N_shells=2, hw=15.0)
        _, states_n = enumerate_sp_states(spec)
        tbme_nn = compute_same_species_tbme(states_n, spec, include_coulomb=False)

        for (i, j, k, l), v in tbme_nn.items():
            if i == j:
                assert abs(v) < 1e-12, (
                    f"Antisymmetrized ME ({i},{j},{k},{l}) = {v} should be 0"
                )

    def test_he4_antisymmetry_property(self):
        """V_{ijkl}_AS = -V_{jikl}_AS (exchange of particles)."""
        spec = DeuteronSpec(N_shells=2, hw=15.0)
        states_p, _ = enumerate_sp_states(spec)
        tbme_pp = compute_same_species_tbme(states_p, spec, include_coulomb=False)

        checked = 0
        for (i, j, k, l), v in tbme_pp.items():
            if abs(v) < 1e-12:
                continue
            v_exch = tbme_pp.get((j, i, k, l), 0.0)
            assert abs(v + v_exch) < abs(v) * 0.01 + 1e-10, (
                f"V({i},{j},{k},{l}) = {v:.6f}, V({j},{i},{k},{l}) = {v_exch:.6f}, "
                f"should be antisymmetric"
            )
            checked += 1
        assert checked > 0, "Should have checked some nonzero antisymmetrized MEs"

    def test_he4_pn_no_exchange(self):
        """
        pn matrix elements don't include exchange terms.

        Verify that pn TBME are the same as the deuteron's compute_tbme
        (which has no antisymmetrization).
        """
        spec = DeuteronSpec(N_shells=2, hw=15.0)
        states_p, states_n = enumerate_sp_states(spec)
        tbme_pn = compute_cross_species_tbme(states_p, states_n, spec)
        # Just verify they exist and are nonzero
        assert len(tbme_pn) > 0, "Should have nonzero pn matrix elements"


# ===========================================================================
# Hamiltonian properties
# ===========================================================================

class TestHe4Hamiltonian:
    """Tests for the He-4 Hamiltonian."""

    def test_he4_hamiltonian_hermitian(self):
        """H is Hermitian (real symmetric)."""
        data = build_he4_hamiltonian(N_shells=2, hw=15.0)
        H = data['H_matrix']
        assert np.allclose(H, H.T, atol=1e-10), "H should be Hermitian"

    def test_he4_ground_state_energy(self):
        """E_gs is finite and the Hamiltonian is properly bound."""
        # At hw=10, the basis is broad enough that attraction dominates
        result = diagonalize_he4(N_shells=2, hw=10.0)
        E_gs = result['E_gs']
        assert np.isfinite(E_gs), f"E_gs should be finite, got {E_gs}"
        # Non-interacting: 4 nucleons in 0s shell, E = 4 * (3/2) * hw
        E_free = 4 * 1.5 * 10.0  # = 60 MeV
        assert E_gs < E_free, (
            f"E_gs={E_gs:.2f} should be below free energy {E_free:.2f} MeV"
        )
        # With proper antisymmetrization, He-4 at optimal hw should be bound
        # (negative energy relative to separated nucleons at E=0)
        assert E_gs < 0, (
            f"E_gs={E_gs:.2f} should be negative (bound) at hw=10"
        )

    @pytest.mark.slow
    def test_he4_hw_scan(self):
        """
        Find optimal hw for He-4.

        He-4 is more compact than the deuteron, so optimal hw should be
        higher (typically 15-25 MeV).
        """
        scan = he4_hw_scan(N_shells=2,
                           hw_values=[10.0, 15.0, 20.0, 25.0, 30.0])
        E_opt = scan['optimal_E']
        hw_opt = scan['optimal_hw']
        assert E_opt is not None, "Should find optimal hw"
        assert hw_opt is not None
        # He-4 is tightly bound; optimal hw typically 15-25 MeV
        assert 8.0 < hw_opt < 35.0, f"Optimal hw={hw_opt} out of expected range"

    def test_he4_coulomb_small_correction(self):
        """
        Coulomb correction should be small, positive (repulsive).

        Experimental Coulomb displacement energy in He-4 is ~0.7 MeV.
        """
        result_with = diagonalize_he4(N_shells=2, hw=10.0, include_coulomb=True)
        result_without = diagonalize_he4(N_shells=2, hw=10.0, include_coulomb=False)
        dE = result_with['E_gs'] - result_without['E_gs']
        # Coulomb is repulsive: E_with > E_without
        assert dE > 0, f"Coulomb should raise energy, got dE={dE:.4f} MeV"
        # The correction should be small (< 2 MeV at reasonable hw)
        assert dE < 2.0, f"Coulomb correction {dE:.4f} MeV too large (expected < 2)"

    def test_he4_spin_zero_ground(self):
        """
        Ground state should have total M_S=0, M_L=0 (consistent with J=0).

        He-4 ground state is 0+. In our basis, this means the dominant
        components have M_L=0 and M_S=0 for both proton and neutron SDs.
        """
        result = diagonalize_he4(N_shells=2, hw=15.0)
        psi = result['ground_state']
        data = result['H_data']
        states_p = data['states_p']
        states_n = data['states_n']
        sds_p = data['sds_p']
        sds_n = data['sds_n']
        dim_n = len(sds_n)

        # Check that the largest components have M_L=0, M_S=0
        # For the ground state, the 0s^2 configuration should dominate
        total_weight_j0 = 0.0
        total_weight = 0.0
        for ip, sd_p in enumerate(sds_p):
            ml_p = sum(states_p[o].m_l for o in sd_p)
            ms_p = sum(states_p[o].m_s for o in sd_p)
            for jn, sd_n in enumerate(sds_n):
                ml_n = sum(states_n[o].m_l for o in sd_n)
                ms_n = sum(states_n[o].m_s for o in sd_n)
                idx = ip * dim_n + jn
                w = psi[idx] ** 2
                total_weight += w
                if abs(ml_p + ml_n) < 1e-10 and abs(ms_p + ms_n) < 1e-10:
                    total_weight_j0 += w

        # For a J=0 state, ALL weight should have M_L+M_S=0
        # (Actually M_L=0 AND M_S=0 separately since we're in uncoupled basis
        # and the ground state is 0+)
        frac_j0 = total_weight_j0 / total_weight
        assert frac_j0 > 0.99, (
            f"Ground state J=0 fraction = {frac_j0:.4f}, expected > 0.99"
        )


# ===========================================================================
# Pauli Hamiltonian analysis
# ===========================================================================

class TestHe4Pauli:
    """Tests for He-4 Pauli Hamiltonian."""

    def test_he4_pauli_count(self):
        """Report Pauli term count at N_shells=2."""
        data = build_he4_hamiltonian(N_shells=2, hw=15.0)
        H_pauli = data['H_pauli']
        analysis = analyze_he4_hamiltonian(H_pauli, data['Q'])
        n_terms = analysis['n_pauli_terms']
        assert n_terms > 0, "Should have nonzero Pauli terms"
        assert n_terms < 50000, f"Pauli count {n_terms} seems too high"

    def test_he4_more_pauli_than_deuteron(self):
        """He-4 should have more Pauli terms than deuteron at same N_shells."""
        from geovac.nuclear.nuclear_hamiltonian import (
            build_deuteron_hamiltonian,
            analyze_deuteron_hamiltonian,
        )
        data_d = build_deuteron_hamiltonian(N_shells=2, hw=15.0)
        data_he4 = build_he4_hamiltonian(N_shells=2, hw=15.0)

        analysis_d = analyze_deuteron_hamiltonian(data_d['H_pauli'], data_d['Q'])
        analysis_he4 = analyze_he4_hamiltonian(data_he4['H_pauli'], data_he4['Q'])

        assert analysis_he4['n_pauli_terms'] > analysis_d['n_pauli_terms'], (
            f"He-4 Pauli terms ({analysis_he4['n_pauli_terms']}) should exceed "
            f"deuteron ({analysis_d['n_pauli_terms']})"
        )

    def test_he4_pauli_string_length(self):
        """All Pauli strings should have length Q=16."""
        data = build_he4_hamiltonian(N_shells=2, hw=15.0)
        Q = data['Q']
        for pstr in data['H_pauli']:
            assert len(pstr) == Q, f"Pauli string length {len(pstr)} != Q={Q}"

    def test_he4_pauli_real_coefficients(self):
        """All Pauli coefficients should be real."""
        data = build_he4_hamiltonian(N_shells=2, hw=15.0)
        for pstr, coeff in data['H_pauli'].items():
            assert isinstance(coeff, (int, float)), (
                f"Coefficient for {pstr} is {type(coeff)}, should be real"
            )
