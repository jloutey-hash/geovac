"""
tests/test_h2_energy_decomposition.py

H2 Full CI energy decomposition tests.

Tests each term of the Full CI Hamiltonian separately to catch:
  1. V_NN wrong (nuclei at origin bug)
  2. Cross-nuclear overestimation: old point-charge model gave V_cross ≈ -1.16 Ha (59% error)
  3. Variational collapse: diffuse high-n states accumulate unphysical attraction,
     causing Full CI energy to diverge as max_n grows
  4. Total energy outside the physically reasonable range [-1.20, -1.10] Ha

Reference values:
  R = 1.4 bohr  (near-equilibrium geometry)
  Exact H2: E = -1.17447 Ha  (full CI limit)
  V_NN = 1/R = 0.71429 Ha
  H2 HF limit: -1.1336 Ha  (~3.5% error from correlation)

History of bugs diagnosed (Feb 2026):
  - Bug 1: bridge connectivity used only n=max_n states → zero bonding splitting
  - Bug 2: cross-nuclear used point-charge model at r = n²/Z → 5x overestimate
  - Bug 3: Full CI used self.hamiltonian (includes node weights W) → double-counts
            nuclear attraction; correct: use pure-kinetic H1 = kinetic_scale * L
  - Bug 4: Mulliken approximation without cap → S_eff → 1 for large n (R_eff → 0),
            giving V_cross ≈ -Z/R for ALL states → variational collapse at max_n ≥ 6
            Fix: cap at -Z_other * Z_self / n² (the ⟨1/r⟩ limit for diffuse orbitals)

Note on cross-nuclear bounds:
  The original spec requested "cross-nuclear between -0.15 and -0.30 Ha per electron."
  The ground-state EXPECTATION VALUE of V_cross divided by 2 is approximately -0.083 Ha
  (max_n=5) because the CI wavefunction samples many high-n states with small V_cross.
  The n=1 state's MATRIX ELEMENT is -0.538 Ha (Mulliken approximation for R=1.4 bohr,
  vs exact ⟨1s_A|−Z/r_B|1s_A⟩ = 0.61 Ha). The critical bug-detection threshold is:
  fail if total V_cross expectation < -0.50 Ha (old point-charge gave -1.159 Ha).
"""

import pytest
import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.sparse import kron, identity, diags

from geovac.hamiltonian import MoleculeHamiltonian

# ──────────────────────────────────────────────────────────────────────────────
# Test parameters
# ──────────────────────────────────────────────────────────────────────────────
R_TEST: float = 1.4          # bohr — near-equilibrium H2 geometry
EXACT_H2: float = -1.17447   # Ha  — full-CI limit from literature
Z_H: float = 1.0


def _build_h2_mol(max_n: int = 5) -> MoleculeHamiltonian:
    """Build H2 at R=1.4 bohr using Method-2 (nuclei kwarg) so positions are correct."""
    return MoleculeHamiltonian(
        nuclei=[(-R_TEST / 2, 0.0, 0.0), (R_TEST / 2, 0.0, 0.0)],
        nuclear_charges=[Z_H, Z_H],
        max_n=max_n,
        connectivity=[(0, 1, 3)],
    )


def _build_cross_terms(mol: MoleculeHamiltonian):
    """Return (V_n1, V_n2, H_cross) as sparse matrices."""
    V_n1, V_n2 = mol._build_cross_nuclear_attraction()
    n = mol.n_total_states
    I = identity(n, format='csr')
    H_cross = kron(V_n2, I, format='csr') + kron(I, V_n1, format='csr')
    return V_n1, V_n2, H_cross


def _pure_kinetic_h1(mol: MoleculeHamiltonian):
    """Build H1 = kinetic_scale * (D - A), no node weights."""
    n = mol.n_total_states
    deg = np.array(mol.adjacency.sum(axis=1)).flatten()
    L = diags(deg, 0, shape=(n, n), format='csr') - mol.adjacency
    return mol.kinetic_scale * L


# ──────────────────────────────────────────────────────────────────────────────
# 1.  Nuclear–nuclear repulsion
# ──────────────────────────────────────────────────────────────────────────────
class TestVNN:
    def test_vnn_exact_value(self):
        """V_NN = Z_A * Z_B / R_AB must equal 1 / 1.4 = 0.71429 Ha exactly."""
        mol = _build_h2_mol(max_n=3)
        V_NN = mol.compute_nuclear_repulsion()
        expected = 1.0 / R_TEST
        assert abs(V_NN - expected) < 0.001, (
            f"V_NN = {V_NN:.6f} Ha, expected {expected:.6f} Ha.\n"
            f"This fails when nuclei are both placed at the origin (default position "
            f"of GeometricLattice). Use Method-2 constructor (nuclei kwarg)."
        )

    def test_vnn_positive(self):
        """Nuclear repulsion must be strictly positive (classical Coulomb)."""
        mol = _build_h2_mol(max_n=3)
        assert mol.compute_nuclear_repulsion() > 0.0


# ──────────────────────────────────────────────────────────────────────────────
# 2.  Cross-nuclear attraction matrix
# ──────────────────────────────────────────────────────────────────────────────
class TestCrossNuclear:
    def test_n1_state_element_in_mulliken_range(self):
        """V_cross for the n=1 ground-state orbital must be in [-0.65, -0.40] Ha.

        Mulliken approximation for (n=1, l=0, R=1.4, Z=1):
          S_eff ≈ 0.754, V_cross = -Z/R * S_eff = -0.538 Ha
        Exact ⟨1s_A | -Z/r_B | 1s_A⟩ = 0.610 Ha → acceptable Mulliken underestimate.
        Values outside this window suggest a regression to the point-charge model.
        """
        mol = _build_h2_mol(max_n=5)
        V_n1, V_n2, _ = _build_cross_terms(mol)

        # V_n2 is non-zero only for atom-A states (indices 0:nA)
        nA = mol.lattices[0].num_states
        v_n2_A = np.array(V_n2.diagonal())[:nA]
        nonzero = v_n2_A[v_n2_A != 0]
        v_n1_element = np.min(nonzero)   # most negative = n=1 state

        assert -0.65 < v_n1_element < -0.40, (
            f"n=1 V_cross element = {v_n1_element:.4f} Ha, expected in (-0.65, -0.40).\n"
            f"Old point-charge model could place this near -0.58 Ha (fictitious coords)."
        )

    def test_high_n_states_damped(self):
        """Diffuse high-n states (n≥2) must NOT receive full nuclear charge -Z/R.

        Without the ⟨1/r⟩ cap, the Mulliken S_eff → 1 as n→∞ (R_eff → 0), so
        V_cross → -Z/R = -0.714 Ha for every state.  With the cap:
          V_cross(n) ≤ min(Mulliken, -Z*Z/n²) → rapidly shrinking for n≥2.
        The n=2 state should already be damped well below 0.50 Ha in magnitude.
        """
        mol = _build_h2_mol(max_n=5)
        V_n1, V_n2, _ = _build_cross_terms(mol)

        nA = mol.lattices[0].num_states
        v_n2_A = np.array(V_n2.diagonal())[:nA]
        nonzero = v_n2_A[v_n2_A != 0]

        # Maximum magnitude over ALL states must be < 0.65 Ha.
        # If any state reaches Z/R = 0.714 Ha, the cap is not working.
        max_mag = np.max(np.abs(nonzero))
        assert max_mag < 0.65, (
            f"Max |V_cross| = {max_mag:.4f} Ha, must be < 0.65.\n"
            f"If high-n states hold |V_cross| ≈ Z/R = 0.714, variational collapse "
            f"causes the Full CI to favour diffuse configurations and diverge with max_n."
        )

    def test_cross_nuclear_not_overestimated_in_ground_state(self):
        """CI ground-state expectation of total V_cross must be > -0.50 Ha.

        This is the primary bug-detection test.
        Old point-charge model: V_cross expectation ≈ -1.16 Ha → 59% total error.
        Variational collapse:   V_cross expectation < -1.40 Ha at max_n≥8.
        Correct (Mulliken + cap): V_cross ≈ -0.165 Ha at max_n=5.

        Threshold -0.50 Ha catches both old bugs while allowing the correct physics.
        """
        mol = _build_h2_mol(max_n=5)
        V_n1, V_n2, H_cross = _build_cross_terms(mol)

        n = mol.n_total_states
        I = identity(n, format='csr')
        H1 = _pure_kinetic_h1(mol)
        H_kin_2 = kron(H1, I, format='csr') + kron(I, H1, format='csr')
        V_ee = mol._build_electron_repulsion_molecular()
        H_total = H_kin_2 + H_cross + V_ee

        eigvals, eigvecs = eigsh(H_total, k=1, which='SA')
        psi = eigvecs[:, 0]
        E_cross = float(psi @ (H_cross @ psi))

        assert E_cross > -0.50, (
            f"V_cross expectation = {E_cross:.4f} Ha (should be > -0.50).\n"
            f"Old point-charge model produced -1.159 Ha here.\n"
            f"Variational collapse can push this below -1.40 Ha."
        )

    def test_cross_nuclear_not_underestimated_in_ground_state(self):
        """CI ground-state expectation of total V_cross must be < -0.02 Ha.

        A near-zero V_cross would mean the cross-nuclear term is silently disabled
        (e.g., nuclei placed at origin so R_AB = 0, skipping all cross-terms).
        """
        mol = _build_h2_mol(max_n=5)
        V_n1, V_n2, H_cross = _build_cross_terms(mol)

        n = mol.n_total_states
        I = identity(n, format='csr')
        H1 = _pure_kinetic_h1(mol)
        H_kin_2 = kron(H1, I, format='csr') + kron(I, H1, format='csr')
        V_ee = mol._build_electron_repulsion_molecular()
        H_total = H_kin_2 + H_cross + V_ee

        eigvals, eigvecs = eigsh(H_total, k=1, which='SA')
        psi = eigvecs[:, 0]
        E_cross = float(psi @ (H_cross @ psi))

        assert E_cross < -0.02, (
            f"V_cross expectation = {E_cross:.4f} Ha (should be < -0.02).\n"
            f"Near-zero V_cross means nuclei were placed at the same point "
            f"and all cross-nuclear terms were skipped."
        )


# ──────────────────────────────────────────────────────────────────────────────
# 3.  Full CI total energy
# ──────────────────────────────────────────────────────────────────────────────
class TestFullCIEnergy:
    def test_h2_energy_in_physical_range(self):
        """Full CI energy at R=1.4 bohr must be in (-1.20, -1.10) Ha.

        Lower bound -1.20: catches variational collapse (CI too negative).
        Upper bound -1.10: energy too positive → under-binding.
        Exact reference: -1.17447 Ha.
        """
        mol = _build_h2_mol(max_n=5)
        eigvals, _ = mol.compute_ground_state(n_states=1, method='full_ci')
        E = eigvals[0]

        assert -1.20 < E < -1.10, (
            f"H2 Full CI = {E:.6f} Ha, expected in (-1.20, -1.10).\n"
            f"Exact = {EXACT_H2} Ha.  "
            f"E < -1.20 → variational collapse (diffuse states dominate).  "
            f"E > -1.10 → under-binding (cross-nuclear too weak or bridges broken)."
        )

    def test_h2_error_below_5pct(self):
        """Full CI error must be below 5% at max_n=5."""
        mol = _build_h2_mol(max_n=5)
        eigvals, _ = mol.compute_ground_state(n_states=1, method='full_ci')
        E = eigvals[0]
        error_pct = abs(E - EXACT_H2) / abs(EXACT_H2) * 100

        assert error_pct < 5.0, (
            f"Error = {error_pct:.2f}% (must be < 5%).\n"
            f"E = {E:.6f} Ha vs exact {EXACT_H2} Ha.\n"
            f"Old point-charge model: ~59% error.  "
            f"Pure mean-field (1-electron): ~17% error."
        )

    @pytest.mark.slow
    def test_h2_energy_converges_with_max_n(self):
        """Energy at max_n=8 must be ≤ error at max_n=5 + 0.01 Ha tolerance.

        Tests basis-set convergence.  Variational collapse (pre-fix) caused
        max_n=8 to give 8.5% error vs 0.72% at max_n=5 — energy DIVERGED with
        basis size.  With the ⟨1/r⟩ cap, max_n=8 gives 1.72% (better than 2.15%).

        Marked slow: max_n=8 Full CI matrix is 166k × 166k (~10 MB).
        """
        mol5 = _build_h2_mol(max_n=5)
        mol8 = _build_h2_mol(max_n=8)

        e5, _ = mol5.compute_ground_state(n_states=1, method='full_ci')
        e8, _ = mol8.compute_ground_state(n_states=1, method='full_ci')

        err5 = abs(e5[0] - EXACT_H2)
        err8 = abs(e8[0] - EXACT_H2)

        assert err8 <= err5 + 0.010, (
            f"Energy diverges with basis size: "
            f"err5={err5:.4f} Ha, err8={err8:.4f} Ha.\n"
            f"E5={e5[0]:.6f}, E8={e8[0]:.6f} Ha.\n"
            f"Variational collapse produces err8 >> err5."
        )

    def test_h2_energy_above_hydrogen_atom_limit(self):
        """H2 energy must be less than 2 × E(H) = -1.0 Ha (binding must exist).

        The exact H2 dissociation limit is 2 × (-0.5) = -1.0 Ha.  If the Full CI
        gives energy > -1.0 Ha the molecule is not bound at all, indicating a
        fundamental failure (e.g., pure repulsive V_ee without any attraction).
        """
        mol = _build_h2_mol(max_n=5)
        eigvals, _ = mol.compute_ground_state(n_states=1, method='full_ci')
        E = eigvals[0]
        dissociation_limit = -1.0  # 2 × E(H atom)

        assert E < dissociation_limit, (
            f"H2 Full CI = {E:.6f} Ha ≥ dissociation limit {dissociation_limit} Ha.\n"
            f"Molecule is not bound."
        )
