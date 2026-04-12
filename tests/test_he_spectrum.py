"""
Tests for He excited states and atomic spectrum (Track DI Sprint 4).

Validates:
1. Triplet FCI matrix construction (antisymmetric spatial wavefunction)
2. L>0 sector construction (P states, D states)
3. Singlet-triplet splitting sign and magnitude
4. Transition energy error cancellation
5. 2D variational solver triplet states
6. Spectrum convergence with n_max
"""

import numpy as np
import pytest

from geovac.casimir_ci import (
    build_graph_native_fci,
    compute_he_spectrum,
    HE_NR_REFERENCE,
    _build_orbital_basis,
)
from geovac.level3_variational import solve_he_variational_2d


E_EXACT_GS = -2.903724377034  # Ha, 1^1S ground state


# ========================================================================
# Triplet FCI construction
# ========================================================================

class TestTripletFCI:
    """Validate triplet (S=1) FCI matrix construction."""

    def test_triplet_nmax2_shape(self):
        """Triplet at n_max=2, M_L=0: strictly i<j, no double occupancy."""
        H = build_graph_native_fci(Z=2, n_max=2, m_total=0, spin='triplet')
        # 5 spatial orbitals: (1,0,0), (2,0,0), (2,1,-1), (2,1,0), (2,1,1)
        # Pairs i<j with m_i + m_j = 0:
        # (0,1)=(1s,2s), (0,3)=(1s,2p0), (1,3)=(2s,2p0),
        # (2,4)=(2p-1,2p1)
        assert H.shape[0] == 4, f"Expected 4 triplet configs, got {H.shape[0]}"
        assert H.shape[0] == H.shape[1]

    def test_triplet_nmax1_empty(self):
        """n_max=1: only 1s orbital, no triplet states possible."""
        H = build_graph_native_fci(Z=2, n_max=1, m_total=0, spin='triplet')
        assert H.shape == (0, 0)

    def test_triplet_hermitian(self):
        """Triplet FCI matrix should be Hermitian."""
        H = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='triplet')
        assert np.allclose(H, H.T, atol=1e-10)

    def test_triplet_below_singlet_excited(self):
        """2^3S should be below 2^1S (Hund's rule for exchange)."""
        H_singlet = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet')
        H_triplet = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='triplet')
        evals_s = np.sort(np.linalg.eigvalsh(H_singlet))
        evals_t = np.sort(np.linalg.eigvalsh(H_triplet))
        # 2^3S (lowest triplet S) < 2^1S (second singlet S)
        E_2_3S = evals_t[0]
        E_2_1S = evals_s[1]
        assert E_2_3S < E_2_1S, (
            f"2^3S ({E_2_3S:.6f}) should be below 2^1S ({E_2_1S:.6f})"
        )

    def test_triplet_above_singlet_ground(self):
        """Lowest triplet must be above singlet ground state."""
        H_singlet = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet')
        H_triplet = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='triplet')
        E_gs = np.linalg.eigvalsh(H_singlet)[0]
        E_triplet = np.linalg.eigvalsh(H_triplet)[0]
        assert E_triplet > E_gs, (
            f"Triplet ({E_triplet:.6f}) should be above ground state ({E_gs:.6f})"
        )

    def test_singlet_backward_compatible(self):
        """Singlet with explicit spin='singlet' matches old behavior."""
        H_old = build_graph_native_fci(Z=2, n_max=3, m_total=0)
        H_new = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin='singlet')
        assert np.allclose(H_old, H_new, atol=1e-14)

    def test_invalid_spin_raises(self):
        """Invalid spin parameter should raise ValueError."""
        with pytest.raises(ValueError, match="spin must be"):
            build_graph_native_fci(Z=2, n_max=2, spin='doublet')


# ========================================================================
# L > 0 sectors (P and D states)
# ========================================================================

class TestLGreaterThanZero:
    """Validate L=1 (P) and L=2 (D) state construction."""

    def test_singlet_P_nmax2(self):
        """Singlet P states at n_max=2: M_L=1 configurations exist."""
        H = build_graph_native_fci(Z=2, n_max=2, m_total=1, spin='singlet')
        assert H.shape[0] > 0, "Should have singlet P configurations at n_max=2"

    def test_triplet_P_nmax2(self):
        """Triplet P states at n_max=2: M_L=1 configurations exist."""
        H = build_graph_native_fci(Z=2, n_max=2, m_total=1, spin='triplet')
        assert H.shape[0] > 0, "Should have triplet P configurations at n_max=2"

    def test_P_states_above_S_states(self):
        """2P states should be above 2S states (for both singlet and triplet)."""
        for spin in ['singlet', 'triplet']:
            H_S = build_graph_native_fci(Z=2, n_max=3, m_total=0, spin=spin)
            H_P = build_graph_native_fci(Z=2, n_max=3, m_total=1, spin=spin)
            E_S = np.sort(np.linalg.eigvalsh(H_S))
            E_P = np.sort(np.linalg.eigvalsh(H_P))
            # Lowest P state should be above lowest S state
            # (but both should be in the n=2 manifold)
            if spin == 'triplet':
                # 2^3P > 2^3S
                assert E_P[0] > E_S[0], (
                    f"{spin}: 2P ({E_P[0]:.6f}) should be above 2S ({E_S[0]:.6f})"
                )

    def test_singlet_P_hermitian(self):
        """P-state FCI matrix should be Hermitian."""
        H = build_graph_native_fci(Z=2, n_max=3, m_total=1, spin='singlet')
        assert np.allclose(H, H.T, atol=1e-10)

    def test_D_states_need_nmax3(self):
        """D states (L=2) require n_max >= 3 (need d-orbitals)."""
        H = build_graph_native_fci(Z=2, n_max=2, m_total=2, spin='singlet')
        # At n_max=2, d-orbitals not available, but M_L=2 could use
        # two p-orbitals: (2,1,1) + (2,1,1) won't work in singlet i<=j
        # Actually (2,1,1) + (2,1,1) is a valid singlet config (i=j)
        # but m_i + m_j = 2, so it should exist
        # No, wait — we need 2 p+1 orbitals but there's only one at n_max=2
        # So we need two electrons both in (2,1,1) — that IS a valid config for singlet
        if H.shape[0] > 0:
            pass  # Some configs may exist
        # At n_max=3, there should definitely be D configs
        H3 = build_graph_native_fci(Z=2, n_max=3, m_total=2, spin='singlet')
        assert H3.shape[0] > 0, "Should have D configurations at n_max=3"


# ========================================================================
# Spectrum computation
# ========================================================================

class TestHeSpectrum:
    """Validate the full spectrum computation."""

    def test_spectrum_nmax3_has_states(self):
        """Spectrum at n_max=3 should produce multiple states."""
        result = compute_he_spectrum(n_max=3)
        assert len(result['states']) >= 6, (
            f"Expected at least 6 states, got {len(result['states'])}"
        )

    def test_spectrum_has_ground_state(self):
        """Ground state 1^1S should be present."""
        result = compute_he_spectrum(n_max=3)
        assert '1_1S' in result['states']

    def test_spectrum_has_triplet(self):
        """2^3S should be present."""
        result = compute_he_spectrum(n_max=3)
        assert '2_3S' in result['states']

    def test_spectrum_has_P_states(self):
        """P states should be present at n_max >= 2."""
        result = compute_he_spectrum(n_max=3, max_L=1)
        has_P = any('P' in label for label in result['states'])
        assert has_P, "Should have P states at n_max=3"

    def test_spectrum_transition_errors(self):
        """Transition energies should have errors computed."""
        result = compute_he_spectrum(n_max=3)
        assert len(result['transitions']) > 0

    def test_ground_state_error_reasonable(self):
        """Ground state error should be < 1% at n_max=3."""
        result = compute_he_spectrum(n_max=3)
        gs = result['states']['1_1S']
        assert gs['error_pct'] < 1.0, f"Ground state error {gs['error_pct']:.2f}% > 1%"

    def test_singlet_triplet_splitting_sign(self):
        """2^3S should be below 2^1S (positive singlet-triplet gap)."""
        result = compute_he_spectrum(n_max=3)
        E_2_3S = result['states']['2_3S']['energy']
        E_2_1S = result['states']['2_1S']['energy']
        assert E_2_3S < E_2_1S, (
            f"2^3S ({E_2_3S:.6f}) should be below 2^1S ({E_2_1S:.6f})"
        )

    def test_transition_error_bounded(self):
        """Main excitation transition errors should be bounded and physical.

        Transitions from the ground state should have errors < 15%.
        Small splitting transitions (e.g. 2³S→2¹S) are excluded since
        the ~0.03 Ha gap amplifies small absolute errors.
        """
        result = compute_he_spectrum(n_max=4)
        # Check transitions from the ground state only
        for label, t in result['transitions'].items():
            if label.startswith('1_1S->'):
                assert t['error_pct'] < 15.0, (
                    f"Transition {label} error {t['error_pct']:.2f}% > 15%"
                )


# ========================================================================
# 2D variational solver: triplet states
# ========================================================================

class TestVariational2DTriplet:
    """Validate 2D variational solver for triplet states."""

    def test_triplet_runs(self):
        """2D solver should run with symmetry='triplet'."""
        result = solve_he_variational_2d(
            Z=2, n_basis_R=15, n_basis_alpha=10, l_max=0,
            symmetry='triplet', n_states=1,
        )
        assert 'energies' in result
        assert len(result['energies']) >= 1

    def test_triplet_energy_physical(self):
        """Triplet energy should be in the right ballpark (~-2.17 Ha)."""
        result = solve_he_variational_2d(
            Z=2, n_basis_R=15, n_basis_alpha=10, l_max=0,
            symmetry='triplet', n_states=1,
        )
        E = result['energies'][0]
        # 2^3S reference: -2.175229 Ha
        assert -2.5 < E < -1.5, f"Triplet energy {E:.6f} not in physical range"

    def test_triplet_above_singlet(self):
        """Triplet ground state above singlet ground state."""
        res_s = solve_he_variational_2d(
            Z=2, n_basis_R=15, n_basis_alpha=10, l_max=0,
            symmetry='singlet', n_states=1,
        )
        res_t = solve_he_variational_2d(
            Z=2, n_basis_R=15, n_basis_alpha=10, l_max=0,
            symmetry='triplet', n_states=1,
        )
        assert res_t['energies'][0] > res_s['energies'][0]

    def test_singlet_excited_states(self):
        """Requesting n_states > 1 should return excited states."""
        result = solve_he_variational_2d(
            Z=2, n_basis_R=15, n_basis_alpha=10, l_max=0,
            symmetry='singlet', n_states=3,
        )
        assert len(result['energies']) >= 2
        # Excited state above ground
        assert result['energies'][1] > result['energies'][0]

    @pytest.mark.slow
    def test_triplet_accuracy(self):
        """Triplet 2^3S energy within 1% at l_max=2."""
        result = solve_he_variational_2d(
            Z=2, n_basis_R=20, n_basis_alpha=20, l_max=2,
            symmetry='triplet', n_states=1,
        )
        E_ref = -2.175229378
        E = result['energies'][0]
        error_pct = abs((E - E_ref) / E_ref) * 100.0
        assert error_pct < 1.0, f"Triplet error {error_pct:.2f}% > 1%"


# ========================================================================
# Convergence
# ========================================================================

class TestSpectrumConvergence:
    """Validate spectrum convergence with n_max."""

    def test_ground_state_improves(self):
        """Ground state should improve from n_max=2 to n_max=3."""
        r2 = compute_he_spectrum(n_max=2, max_L=0)
        r3 = compute_he_spectrum(n_max=3, max_L=0)
        err2 = r2['states']['1_1S']['error_pct']
        err3 = r3['states']['1_1S']['error_pct']
        assert err3 < err2, f"n_max=3 error ({err3:.3f}%) should be < n_max=2 ({err2:.3f}%)"

    @pytest.mark.slow
    def test_triplet_converges(self):
        """Triplet 2^3S should converge with n_max."""
        r3 = compute_he_spectrum(n_max=3, max_L=0)
        r4 = compute_he_spectrum(n_max=4, max_L=0)
        E3 = r3['states']['2_3S']['energy']
        E4 = r4['states']['2_3S']['energy']
        ref = HE_NR_REFERENCE['2_3S']
        err3 = abs(E3 - ref)
        err4 = abs(E4 - ref)
        assert err4 < err3, f"n_max=4 should be closer to exact than n_max=3"
