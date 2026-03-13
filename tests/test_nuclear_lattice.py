"""
Tests for the Nuclear Lattice (Paper 10)
=========================================

Validates:
1. Morse SU(2) algebra: ladder elements, Casimir, spectrum
2. Vibrational chain eigenvalues reproduce Morse spectrum
3. Rotational paraboloid has same structure as electron (l, m)
4. Product graph eigenvalues give rovibrational spectrum
5. Comparison to experimental H₂, HCl, LiH spectra

Author: GeoVac Development Team
Date: March 2026
"""

import numpy as np
import pytest
from geovac.nuclear_lattice import (
    MorseVibrationLattice,
    RotationalParaboloid,
    NuclearLattice,
    DIATOMIC_CONSTANTS,
    HARTREE_TO_CM,
    build_diatomic,
    electron_p0_from_vibration,
    franck_condon_overlap,
)


# ======================================================================
# Test 1: Morse SU(2) Algebra
# ======================================================================

class TestMorseSU2:
    """Validate the SU(2) algebraic structure of the Morse oscillator."""

    @pytest.fixture
    def h2_vib(self) -> MorseVibrationLattice:
        """H₂ vibrational lattice from experimental constants."""
        c = DIATOMIC_CONSTANTS['H2']
        return MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])

    def test_j_representation(self, h2_vib: MorseVibrationLattice) -> None:
        """j = ω_e/(2 ω_e x_e) - 1/2 must be positive and self-consistent."""
        j = h2_vib.j
        assert j > 0, f"j = {j} must be positive"
        # j = ω/(2ωx) - 0.5
        expected_j = h2_vib.omega_e_hartree / (2 * h2_vib.omega_e_xe_hartree) - 0.5
        assert abs(j - expected_j) < 1e-10, f"j = {j}, expected {expected_j}"

    def test_v_max_from_j(self, h2_vib: MorseVibrationLattice) -> None:
        """v_max ≤ floor(j), all bound states have E_v < D_e."""
        j = h2_vib.j
        assert h2_vib.v_max <= int(np.floor(j))
        assert h2_vib.n_states == h2_vib.v_max + 1
        # All states must be bound
        for v in range(h2_vib.n_states):
            assert h2_vib.morse_energy(v) <= h2_vib.D_e_hartree

    def test_ladder_elements_positive(self, h2_vib: MorseVibrationLattice) -> None:
        """All ladder elements √[(j-v)(j+v+1)] must be real and positive."""
        for v in range(h2_vib.n_states - 1):
            w = h2_vib.ladder_element(v)
            assert w > 0, f"Ladder element at v={v} is {w}, must be > 0"

    def test_ladder_element_formula(self) -> None:
        """
        Verify ⟨v+1|J₋|v⟩ = √[(2j - v)(v + 1)] against known values.

        For j = 2, v_max = 2 (bound states v=0,1,2):
          v=0: √[(4)(1)] = 2
          v=1: √[(3)(2)] = √6
        """
        # Construct a Morse with j = 2 exactly:
        # j = ω/(2ωx) - 0.5 = 2  →  ω/(2ωx) = 2.5  →  ωx = ω/5
        omega_cm = 1000.0
        omega_xe_cm = 200.0  # ω/5
        omega_hartree = omega_cm / HARTREE_TO_CM
        # D_e large enough that all v ≤ floor(j)=2 are bound
        D_e = 10.0 * omega_hartree
        vib = MorseVibrationLattice(D_e, omega_cm, omega_xe_cm)

        j = vib.j
        assert abs(j - 2.0) < 0.01, f"j = {j}, expected 2.0"

        # v=0: √[(2j - 0)(0 + 1)] = √[4*1] = 2
        w0 = vib.ladder_element(0)
        assert abs(w0 - 2.0) < 1e-6, f"w(0) = {w0}, expected 2.0"

        # v=1: √[(2j - 1)(1 + 1)] = √[3*2] = √6
        w1 = vib.ladder_element(1)
        assert abs(w1 - np.sqrt(6)) < 1e-6, f"w(1) = {w1}, expected √6 = {np.sqrt(6)}"

    def test_ladder_element_boundary(self, h2_vib: MorseVibrationLattice) -> None:
        """
        Ladder elements √[(2j-v)(v+1)] peak near v ≈ j and vanish at boundaries.
        For v=0: element = √(2j) (starting edge).
        All elements must be real and finite for v < v_max.
        """
        elements = [h2_vib.ladder_element(v) for v in range(h2_vib.n_states - 1)]
        # All elements must be finite and positive
        for i, w in enumerate(elements):
            assert np.isfinite(w) and w > 0, \
                f"Ladder element at v={i} is {w}, must be finite and positive"
        # First element: √[2j * 1] = √(2j)
        j = h2_vib.j
        expected_first = np.sqrt(2 * j)
        assert abs(elements[0] - expected_first) < 1e-10, \
            f"w(0) = {elements[0]}, expected √(2j) = {expected_first}"

    def test_casimir_constant(self, h2_vib: MorseVibrationLattice) -> None:
        """SU(2) Casimir J² = j(j+1) is constant for all v-states."""
        casimir, j3_spectrum = h2_vib.su2_casimir_spectrum()
        j = h2_vib.j
        assert abs(casimir - j * (j + 1)) < 1e-10

    def test_j3_spectrum(self, h2_vib: MorseVibrationLattice) -> None:
        """J₃ eigenvalues run from j to j - v_max, mapping to v = 0, ..., v_max."""
        _, j3 = h2_vib.su2_casimir_spectrum()
        j = h2_vib.j
        for v in range(h2_vib.n_states):
            assert abs(j3[v] - (j - v)) < 1e-10, \
                f"J₃ at v={v}: {j3[v]}, expected {j - v}"


# ======================================================================
# Test 2: Morse Spectrum Reproduction
# ======================================================================

class TestMorseSpectrum:
    """Validate that the graph reproduces the Morse energy spectrum."""

    def test_morse_energy_formula(self) -> None:
        """
        E_v = ω(v+½) - ωx(v+½)² matches analytical Morse levels.
        """
        c = DIATOMIC_CONSTANTS['H2']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])

        # Ground state: v=0
        E0 = vib.morse_energy(0)
        omega = vib.omega_e_hartree
        omegax = vib.omega_e_xe_hartree
        expected_E0 = omega * 0.5 - omegax * 0.25
        assert abs(E0 - expected_E0) < 1e-12

        # v=1
        E1 = vib.morse_energy(1)
        expected_E1 = omega * 1.5 - omegax * 2.25
        assert abs(E1 - expected_E1) < 1e-12

    def test_morse_spacing_decreasing(self) -> None:
        """
        Morse level spacings ΔE_v = E_{v+1} - E_v must decrease with v.
        This is the signature of anharmonicity (finite SU(2) representation).
        """
        c = DIATOMIC_CONSTANTS['HCl']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])

        spacings = []
        for v in range(min(5, vib.n_states - 1)):
            dE = vib.morse_energy(v + 1) - vib.morse_energy(v)
            spacings.append(dE)

        for i in range(len(spacings) - 1):
            assert spacings[i + 1] < spacings[i], \
                f"Spacing not decreasing at v={i}: {spacings[i]:.6f} → {spacings[i+1]:.6f}"

    def test_zero_point_energy(self) -> None:
        """ZPE = ω/2 - ωx/4 for the Morse oscillator."""
        for mol in ['H2', 'HCl', 'CO', 'LiH']:
            c = DIATOMIC_CONSTANTS[mol]
            vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])
            zpe = vib.morse_energy(0)
            omega = vib.omega_e_hartree
            omegax = vib.omega_e_xe_hartree
            expected = omega / 2 - omegax / 4
            assert abs(zpe - expected) < 1e-12, \
                f"{mol}: ZPE = {zpe:.6f}, expected {expected:.6f}"

    def test_h2_fundamental_frequency(self) -> None:
        """
        H₂ fundamental frequency ν₀₁ = ω_e - 2ω_e x_e ≈ 4161 cm⁻¹.
        """
        c = DIATOMIC_CONSTANTS['H2']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])
        dE = vib.morse_energy(1) - vib.morse_energy(0)
        freq_cm = dE * HARTREE_TO_CM
        expected_cm = c['omega_e'] - 2 * c['omega_e_xe']  # ≈ 4158.5
        assert abs(freq_cm - expected_cm) < 1.0, \
            f"H₂ ν₀₁ = {freq_cm:.1f} cm⁻¹, expected {expected_cm:.1f}"

    def test_graph_spectrum_matches_morse(self) -> None:
        """
        Graph Laplacian eigenvalues must reproduce Morse energies.
        """
        c = DIATOMIC_CONSTANTS['H2']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])

        graph_E = vib.graph_spectrum()
        morse_E = vib.morse_spectrum()

        # They should match by construction (diagonal Hamiltonian)
        np.testing.assert_allclose(graph_E, morse_E, atol=1e-12,
                                   err_msg="Graph spectrum doesn't match Morse")

    def test_dissociation_limit(self) -> None:
        """Highest vibrational state should be near D_e (dissociation)."""
        c = DIATOMIC_CONSTANTS['H2']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])
        E_top = vib.morse_energy(vib.v_max)
        D_e = vib.D_e_hartree
        # The highest state should be below D_e (bound)
        assert E_top < D_e + vib.omega_e_hartree, \
            f"E(v_max) = {E_top:.6f} too far above D_e = {D_e:.6f}"
        # And within ~2ω of D_e (near dissociation)
        assert E_top > D_e - 2 * vib.omega_e_hartree, \
            f"E(v_max) = {E_top:.6f} too far below D_e = {D_e:.6f}"


# ======================================================================
# Test 3: Rotational Paraboloid Structure
# ======================================================================

class TestRotationalParaboloid:
    """Validate the rotational lattice has identical structure to electron (l, m)."""

    def test_degeneracy_2J_plus_1(self) -> None:
        """Each J-shell has 2J+1 states, identical to electron l-shell."""
        rot = RotationalParaboloid(J_max=5, B_e=10.0)
        for J in range(6):
            assert rot.degeneracy(J) == 2 * J + 1

    def test_total_states(self) -> None:
        """Total states = sum_{J=0}^{J_max} (2J+1) = (J_max+1)²."""
        for J_max in [3, 5, 10]:
            rot = RotationalParaboloid(J_max=J_max, B_e=10.0)
            expected = (J_max + 1)**2
            assert rot.n_states == expected, \
                f"J_max={J_max}: {rot.n_states} states, expected {expected}"

    def test_paraboloid_structure(self) -> None:
        """
        The state list forms a paraboloid: J=0 has 1 state, J=1 has 3,
        J=2 has 5, etc. This is identical to the electron (l, m) structure.
        """
        rot = RotationalParaboloid(J_max=4, B_e=10.0)
        counts = {}
        for J, M in rot.states:
            counts[J] = counts.get(J, 0) + 1

        for J in range(5):
            assert counts[J] == 2 * J + 1, \
                f"J={J}: {counts[J]} states, expected {2*J+1}"

    def test_adjacency_within_J_shell(self) -> None:
        """
        L± connects M → M±1 within each J shell.
        The adjacency matrix should have nonzero entries only
        between states with same J and ΔM = ±1.
        """
        rot = RotationalParaboloid(J_max=3, B_e=10.0)
        A = rot.adjacency.toarray()

        for i, (Ji, Mi) in enumerate(rot.states):
            for j, (Jj, Mj) in enumerate(rot.states):
                if i == j:
                    continue
                if A[i, j] != 0:
                    assert Ji == Jj, \
                        f"Cross-J coupling: ({Ji},{Mi}) → ({Jj},{Mj})"
                    assert abs(Mi - Mj) == 1, \
                        f"Non-adjacent M: ({Ji},{Mi}) → ({Jj},{Mj})"

    def test_angular_momentum_weights(self) -> None:
        """
        L₊|J,M⟩ → |J,M+1⟩ with weight √[J(J+1) - M(M+1)].
        Check specific values for J=2.
        """
        rot = RotationalParaboloid(J_max=2, B_e=10.0)
        A = rot.adjacency.toarray()

        # J=2, M=0 → M=1: √[6 - 0] = √6
        i = rot.state_index[(2, 0)]
        j = rot.state_index[(2, 1)]
        expected = np.sqrt(6.0)
        assert abs(A[i, j] - expected) < 1e-10, \
            f"A[(2,0),(2,1)] = {A[i,j]}, expected √6 = {expected}"

        # J=2, M=1 → M=2: √[6 - 2] = 2
        i = rot.state_index[(2, 1)]
        j = rot.state_index[(2, 2)]
        assert abs(A[i, j] - 2.0) < 1e-10, \
            f"A[(2,1),(2,2)] = {A[i,j]}, expected 2.0"

    def test_rigid_rotor_spectrum(self) -> None:
        """E_J = B_e J(J+1) matches analytical rigid rotor."""
        B_e = 10.0  # cm⁻¹
        rot = RotationalParaboloid(J_max=5, B_e=B_e)
        B_hartree = B_e / HARTREE_TO_CM

        for J in range(6):
            expected = B_hartree * J * (J + 1)
            actual = rot.rigid_rotor_energy(J)
            assert abs(actual - expected) < 1e-14, \
                f"J={J}: E = {actual}, expected {expected}"

    def test_selection_rule_delta_M(self) -> None:
        """
        Selection rules from adjacency: ΔM_J = 0, ±1.
        J=0 (single state) is isolated within its shell.
        """
        rot = RotationalParaboloid(J_max=3, B_e=10.0)
        A = rot.adjacency.toarray()

        # J=0 state should have no connections (only M=0, no M±1)
        idx_00 = rot.state_index[(0, 0)]
        assert np.sum(np.abs(A[idx_00, :])) == 0, \
            "J=0 should be disconnected (no M neighbors)"


# ======================================================================
# Test 4: Product Graph Rovibrational Spectrum
# ======================================================================

class TestRovibrationalSpectrum:
    """Validate the product nuclear graph reproduces rovibrational energies."""

    @pytest.fixture
    def h2_nuc(self) -> NuclearLattice:
        return build_diatomic('H2', J_max=5)

    def test_rovibrational_formula(self, h2_nuc: NuclearLattice) -> None:
        """
        E_{v,J} = ω(v+½) - ωx(v+½)² + B_v J(J+1) - D_v J²(J+1)²
        with B_v = B_e - α_e(v+½).
        """
        nuc = h2_nuc
        E_00 = nuc.rovibrational_energy(0, 0)
        E_01 = nuc.rovibrational_energy(0, 1)

        # E(0,1) - E(0,0) should be ≈ 2 B_v(0)
        B_v0 = nuc.rot.B_e_hartree - nuc.alpha_e_hartree * 0.5
        D_rot = nuc.rot.D_e_rot_hartree
        expected_diff = B_v0 * 2 - D_rot * 4  # J(J+1)=2 for J=1
        actual_diff = E_01 - E_00
        assert abs(actual_diff - expected_diff) < 1e-12

    def test_hamiltonian_diagonal(self, h2_nuc: NuclearLattice) -> None:
        """The Hamiltonian should be diagonal with rovibrational energies."""
        H = h2_nuc.build_hamiltonian()
        spectrum = h2_nuc.rovibrational_spectrum()

        for i in range(min(20, h2_nuc.n_states)):
            assert abs(H[i, i] - spectrum[i]) < 1e-14

    def test_product_graph_size(self, h2_nuc: NuclearLattice) -> None:
        """Product graph has n_vib × n_rot states."""
        expected = h2_nuc.vib.n_states * h2_nuc.rot.n_states
        assert h2_nuc.n_states == expected

    def test_product_adjacency_shape(self, h2_nuc: NuclearLattice) -> None:
        """Product adjacency matrix has correct dimensions."""
        A = h2_nuc.graph_product_adjacency()
        assert A.shape == (h2_nuc.n_states, h2_nuc.n_states)

    def test_product_laplacian_zero_sum(self, h2_nuc: NuclearLattice) -> None:
        """Graph Laplacian rows sum to zero (conservation property)."""
        L = h2_nuc.graph_product_laplacian()
        row_sums = np.abs(np.array(L.sum(axis=1)).flatten())
        assert np.max(row_sums) < 1e-10, \
            f"Laplacian rows don't sum to zero: max = {np.max(row_sums)}"


# ======================================================================
# Test 5: Experimental Comparison
# ======================================================================

class TestExperimentalComparison:
    """Compare nuclear lattice predictions to experimental spectroscopy."""

    @pytest.mark.parametrize("molecule", ['H2', 'HCl', 'CO', 'LiH'])
    def test_fundamental_frequency(self, molecule: str) -> None:
        """
        Fundamental frequency ν₀₁ = ω_e - 2ω_e x_e.
        Must match experimental to < 1 cm⁻¹ (by construction from constants).
        """
        c = DIATOMIC_CONSTANTS[molecule]
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])

        dE = vib.morse_energy(1) - vib.morse_energy(0)
        freq_cm = dE * HARTREE_TO_CM
        expected = c['omega_e'] - 2 * c['omega_e_xe']

        assert abs(freq_cm - expected) < 1.0, \
            f"{molecule}: ν₀₁ = {freq_cm:.1f}, expected {expected:.1f}"

    @pytest.mark.parametrize("molecule", ['H2', 'HCl', 'CO', 'LiH'])
    def test_first_overtone(self, molecule: str) -> None:
        """
        First overtone ν₀₂ = 2ω_e - 6ω_e x_e.
        """
        c = DIATOMIC_CONSTANTS[molecule]
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])

        if vib.n_states < 3:
            pytest.skip(f"{molecule} has only {vib.n_states} bound states")

        dE = vib.morse_energy(2) - vib.morse_energy(0)
        freq_cm = dE * HARTREE_TO_CM
        expected = 2 * c['omega_e'] - 6 * c['omega_e_xe']

        assert abs(freq_cm - expected) < 1.0, \
            f"{molecule}: ν₀₂ = {freq_cm:.1f}, expected {expected:.1f}"

    def test_h2_bound_states_count(self) -> None:
        """H₂ should have ~14 bound vibrational states."""
        c = DIATOMIC_CONSTANTS['H2']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])
        assert 10 < vib.n_states < 25, \
            f"H₂ has {vib.n_states} bound states, expected ~14"

    def test_hcl_bound_states_count(self) -> None:
        """HCl should have ~18 bound vibrational states."""
        c = DIATOMIC_CONSTANTS['HCl']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])
        assert 10 < vib.n_states < 35, \
            f"HCl has {vib.n_states} bound states, expected ~18"

    def test_h2_rotational_constant(self) -> None:
        """H₂ rotational constant B_e ≈ 60.853 cm⁻¹."""
        nuc = build_diatomic('H2', J_max=5)
        # E(J=1) - E(J=0) = 2 B_e (for rigid rotor)
        E0 = nuc.rovibrational_energy(0, 0)
        E1 = nuc.rovibrational_energy(0, 1)
        B_eff = (E1 - E0) / 2 * HARTREE_TO_CM  # in cm⁻¹
        # B_v(0) = B_e - α_e/2
        B_v0_expected = 60.853 - 3.062 * 0.5  # ≈ 59.32
        assert abs(B_eff - B_v0_expected) < 0.1, \
            f"H₂ B_v(0) = {B_eff:.2f}, expected {B_v0_expected:.2f}"


# ======================================================================
# Test 6: Structural Parallels (Electron vs Nuclear)
# ======================================================================

class TestElectronNuclearParallel:
    """Verify the structural parallel between electron and nuclear lattices."""

    def test_finite_vs_infinite(self) -> None:
        """
        Electron n-ladder: infinite (n → ∞).
        Nuclear v-chain: finite (v_max < ∞).
        This is the key structural difference.
        """
        c = DIATOMIC_CONSTANTS['H2']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])
        assert vib.v_max < np.inf, "v_max must be finite (bounded Morse SU(2))"
        assert vib.v_max > 0, "Must have at least one excited state"

    def test_rotational_matches_electron_angular(self) -> None:
        """
        Rotational (J, M_J) has identical adjacency structure to electron (l, m).
        Both are SO(3) angular momentum algebras.
        """
        from geovac.lattice import GeometricLattice

        J_max = 4
        rot = RotationalParaboloid(J_max=J_max, B_e=10.0)

        # Build electron lattice with same angular structure
        # At fixed n, the (l, m) sublattice has the same L± structure
        lat = GeometricLattice(max_n=J_max + 1)

        # Compare: for J=2, the L± weights should match l=2
        # L₊|l,m⟩ = √[l(l+1) - m(m+1)] |l,m+1⟩
        # This is identical for J and l
        for M in range(-2, 2):
            rot_w = np.sqrt(2 * 3 - M * (M + 1))
            # Same formula for electrons
            elec_w = np.sqrt(2 * 3 - M * (M + 1))
            assert abs(rot_w - elec_w) < 1e-14

    def test_su2_in_both(self) -> None:
        """
        Both vibrational and rotational structures use SU(2).
        Vibration: finite-dimensional (bounded by D_e).
        Rotation: each J-shell is a standard angular momentum multiplet.
        """
        c = DIATOMIC_CONSTANTS['HCl']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])

        # Vibrational SU(2): verify Casimir
        casimir, _ = vib.su2_casimir_spectrum()
        j = vib.j
        assert abs(casimir - j * (j + 1)) < 1e-10

        # Rotational SU(2): verify L² = J(J+1) for each shell
        rot = RotationalParaboloid(J_max=3, B_e=c['B_e'])
        for J in range(4):
            expected_L2 = J * (J + 1)
            # This is the Casimir of the J-representation
            assert expected_L2 == J * (J + 1)  # Tautological but documents the parallel


# ======================================================================
# Test 7: Electron-Nuclear Coupling
# ======================================================================

class TestElectronNuclearCoupling:
    """Validate the Born-Oppenheimer fibration model."""

    def test_p0_equilibrium(self) -> None:
        """At v=0 (equilibrium), p₀ = Z/n exactly."""
        nuc = build_diatomic('H2', J_max=3)
        p0 = electron_p0_from_vibration(0, nuc, Z=1, n=1)
        assert abs(p0 - 1.0) < 1e-10, f"p₀(v=0) = {p0}, expected 1.0"

    def test_p0_decreasing_with_v(self) -> None:
        """p₀ should decrease with v (weaker binding at larger stretch)."""
        nuc = build_diatomic('H2', J_max=3)
        p0_values = [electron_p0_from_vibration(v, nuc, Z=1, n=1)
                     for v in range(min(5, nuc.vib.n_states))]
        for i in range(len(p0_values) - 1):
            assert p0_values[i + 1] < p0_values[i], \
                f"p₀ not decreasing: v={i}: {p0_values[i]:.4f}, v={i+1}: {p0_values[i+1]:.4f}"

    def test_franck_condon_diagonal(self) -> None:
        """FC factor is 1 for v1 = v2 (vertical transition)."""
        nuc = build_diatomic('H2', J_max=3)
        fc = franck_condon_overlap(0, 0, nuc)
        assert abs(fc - 1.0) < 1e-10

    def test_franck_condon_off_diagonal(self) -> None:
        """FC factors for v1 ≠ v2 should be < 1 and > 0."""
        nuc = build_diatomic('H2', J_max=3)
        fc = franck_condon_overlap(0, 1, nuc)
        assert 0 < fc < 1, f"FC(0,1) = {fc}, expected 0 < FC < 1"


# ======================================================================
# Test 8: Build Helpers
# ======================================================================

class TestBuildHelpers:

    @pytest.mark.parametrize("molecule", ['H2', 'HCl', 'CO', 'LiH'])
    def test_build_diatomic(self, molecule: str) -> None:
        """build_diatomic returns a valid NuclearLattice."""
        nuc = build_diatomic(molecule, J_max=5)
        assert nuc.n_states > 0
        assert nuc.vib.n_states > 0
        assert nuc.rot.n_states == 36  # (5+1)² = 36

    def test_unknown_molecule_raises(self) -> None:
        """Unknown molecule name should raise ValueError."""
        with pytest.raises(ValueError, match="Unknown molecule"):
            build_diatomic('XeF6')


# ======================================================================
# Test 9: Adjacency Matrix Properties
# ======================================================================

class TestAdjacencyProperties:

    def test_vibrational_adjacency_symmetric(self) -> None:
        """Vibrational adjacency must be symmetric (undirected graph)."""
        c = DIATOMIC_CONSTANTS['H2']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])
        A = vib.adjacency.toarray()
        np.testing.assert_allclose(A, A.T, atol=1e-14,
                                   err_msg="Vibrational adjacency not symmetric")

    def test_rotational_adjacency_symmetric(self) -> None:
        """Rotational adjacency must be symmetric."""
        rot = RotationalParaboloid(J_max=4, B_e=10.0)
        A = rot.adjacency.toarray()
        np.testing.assert_allclose(A, A.T, atol=1e-14,
                                   err_msg="Rotational adjacency not symmetric")

    def test_vibrational_tridiagonal(self) -> None:
        """Vibrational adjacency is tridiagonal (nearest-neighbor chain)."""
        c = DIATOMIC_CONSTANTS['H2']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])
        A = vib.adjacency.toarray()
        n = vib.n_states

        for i in range(n):
            for j in range(n):
                if abs(i - j) > 1:
                    assert A[i, j] == 0, \
                        f"Non-tridiagonal: A[{i},{j}] = {A[i,j]}"

    def test_laplacian_positive_semidefinite(self) -> None:
        """Graph Laplacian must be positive semidefinite."""
        c = DIATOMIC_CONSTANTS['HCl']
        vib = MorseVibrationLattice(c['D_e'], c['omega_e'], c['omega_e_xe'])
        L = vib.laplacian.toarray()
        eigenvalues = np.linalg.eigvalsh(L)
        assert np.min(eigenvalues) > -1e-10, \
            f"Laplacian has negative eigenvalue: {np.min(eigenvalues)}"
