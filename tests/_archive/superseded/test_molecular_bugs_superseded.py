"""
tests/test_molecular_bugs.py

Regression tests for two specific molecular Hamiltonian bugs diagnosed Feb 2026.

Bug 1 — Cross-nuclear overestimate (point-charge model):
    The old V_cross used fictitious point-charge coordinates (r = n²/Z per state),
    overestimating attraction by ~5x.  Fixed by Mulliken minimal-basis approximation
    with a ⟨1/r⟩ cap that prevents variational collapse at large n.
    Tested by: TestCrossNuclear in test_h2_energy_decomposition.py (detailed),
               test_mulliken_n1_element / test_mulliken_high_n_damped below (regression).

Bug 2 — Bridge connectivity selects only n=max_n shell:
    The old _get_boundary_states_prioritized returned only outermost-shell states,
    leaving the n=1 core (highest ground-state amplitude) unconnected and killing
    the bonding/antibonding splitting.
    Fixed by returning ALL states sorted by (n, l, |m|) — n=1, l=0 first.
    Tested by: test_bridge_connects_n1_core / test_bonding_antibonding_splitting below.

Reference geometry: H2 at R = 1.4 bohr, max_n = 5 (fast), nuclei constructor.
"""

import pytest
import numpy as np
from geovac.hamiltonian import MoleculeHamiltonian
from geovac.lattice import GeometricLattice

R_TEST: float = 1.4   # bohr
EXACT_H2: float = -1.17447  # Ha, full-CI limit


def _build_h2(max_n: int = 5) -> MoleculeHamiltonian:
    """H2 at R=1.4 bohr using the nuclei constructor."""
    return MoleculeHamiltonian(
        nuclei=[(-R_TEST / 2, 0.0, 0.0), (R_TEST / 2, 0.0, 0.0)],
        nuclear_charges=[1.0, 1.0],
        max_n=max_n,
        connectivity=[(0, 1, 3)],
    )


# =============================================================================
# Bug 1 regression: Mulliken cross-nuclear approximation
# =============================================================================

class TestMullikenRegression:
    def test_mulliken_n1_element(self):
        """V_cross for the n=1, l=0 orbital must be in [-0.65, -0.40] Ha.

        Mulliken approximation at R=1.4 bohr, Z=1:
            R_eff = 1.4, S_1s = exp(-1.4)*(1 + 1.4 + 1.4²/3) ≈ 0.754
            ang = 1.0  (l=0)
            V_cross = -1/1.4 * 0.754 = -0.538 Ha
        The old point-charge model placed states at r = n²/Z = 1 bohr for n=1
        and produced a different (often similar) value, but for higher n the
        point-charge model diverged wildly.  We test the correct range here.
        """
        mol = _build_h2(max_n=5)
        V_n1, V_n2 = mol._build_cross_nuclear_attraction()

        # V_n2 is the matrix for cross-attraction to nucleus 2,
        # non-zero only in the atom-A block (first nA rows).
        nA = mol.lattices[0].num_states
        v_diag = np.array(V_n2.diagonal())[:nA]
        nonzero = v_diag[v_diag != 0.0]

        # Most-negative value = n=1, l=0 ground-state orbital
        v_n1_elem = np.min(nonzero)
        assert -0.65 < v_n1_elem < -0.40, (
            f"n=1 V_cross = {v_n1_elem:.4f} Ha — should be in (-0.65, -0.40). "
            f"Regression to point-charge model?"
        )

    def test_mulliken_high_n_damped(self):
        """Diffuse high-n states must NOT receive the full nuclear charge -Z/R.

        Without the ⟨1/r⟩ cap, the Mulliken overlap → 1 as n grows (R_eff → 0),
        giving V_cross → -Z/R = -0.714 Ha for every state.  With the cap,
        |V_cross(n)| ≤ Z²/n² → shrinks rapidly for n ≥ 2.
        """
        mol = _build_h2(max_n=5)
        V_n1, V_n2 = mol._build_cross_nuclear_attraction()

        nA = mol.lattices[0].num_states
        v_diag = np.array(V_n2.diagonal())[:nA]
        nonzero = v_diag[v_diag != 0.0]

        max_mag = np.max(np.abs(nonzero))
        assert max_mag < 0.65, (
            f"Max |V_cross| = {max_mag:.4f} Ha — must be < 0.65. "
            f"If any state reaches Z/R=0.714 the ⟨1/r⟩ cap is not working."
        )

    def test_no_variational_collapse_with_basis(self):
        """Full CI energy at max_n=8 must be ≤ max_n=5 energy + 0.01 Ha.

        Variational collapse (pre-fix) caused energy to INCREASE with max_n
        because diffuse states dominated the CI.  With the cap, larger basis
        converges monotonically toward the exact value.
        """
        mol5 = _build_h2(max_n=5)
        mol8 = _build_h2(max_n=8)
        e5, _ = mol5.compute_ground_state(n_states=1, method='full_ci')
        e8, _ = mol8.compute_ground_state(n_states=1, method='full_ci')
        # Lower energy = more negative = more bound
        assert e8[0] <= e5[0] + 0.010, (
            f"Energy diverges with basis size: E5={e5[0]:.4f}, E8={e8[0]:.4f} Ha. "
            f"Variational collapse causes E8 >> E5 (higher = less bound = wrong)."
        )


# =============================================================================
# Bug 2 regression: bridge connectivity ordering
# =============================================================================

class TestBridgeConnectivityRegression:
    def test_n1_state_is_first_bridge_endpoint(self):
        """_get_boundary_states_prioritized must return the n=1, l=0, m=0 state first.

        This is the ground-state orbital — highest amplitude, most bonding-relevant.
        The old implementation returned only n=max_n states, so (1,0,0) was never
        connected across atoms, killing the bonding/antibonding splitting.
        """
        lattice = GeometricLattice(max_n=5)
        priorities = lattice._get_boundary_states_prioritized()

        # First index must correspond to the (n=1, l=0, m=0) state
        first_state = lattice.states[priorities[0]]
        n, l, m = first_state
        assert (n, l, m) == (1, 0, 0), (
            f"First bridge state is {first_state}, expected (1, 0, 0). "
            f"Bridge ordering bug: n=1 ground state not prioritized."
        )

    def test_ordering_is_sorted_by_n_then_l_then_abs_m(self):
        """All returned indices must be sorted by (n, l, |m|) ascending.

        This guarantees that lower-n, lower-l states — which dominate the
        ground-state wavefunction — are always connected first regardless of
        the requested bridge count.
        """
        lattice = GeometricLattice(max_n=4)
        priorities = lattice._get_boundary_states_prioritized()

        sort_keys = []
        for idx in priorities:
            n, l, m = lattice.states[idx]
            sort_keys.append((n, l, abs(m)))

        assert sort_keys == sorted(sort_keys), (
            "Bridge states not sorted by (n, l, |m|). "
            "High-n outermost-shell states may be selected before n=1 core states."
        )

    def test_n1_state_connected_in_molecule(self):
        """In H2, the n=1 state on atom A must be connected to a state on atom B.

        Verifies that the bridge actually includes the n=1 core and is not pruned
        by the adaptive sparsity mask (bridge weight should be large for n=1).
        """
        mol = _build_h2(max_n=5)
        # The adjacency matrix has off-diagonal blocks for bridges.
        # Extract the top-left A→B cross-block.
        nA = mol.lattices[0].num_states
        adj_coo = mol.adjacency.tocoo()

        # Find any bridge entry: row in [0, nA) and col in [nA, 2*nA)
        bridge_mask = (adj_coo.row < nA) & (adj_coo.col >= nA)
        bridge_rows = adj_coo.row[bridge_mask]

        # The n=1, l=0, m=0 state is index 0 within atom A (first state added)
        n1_local_idx = 0  # (n=1, l=0, m=0) is always the first state
        assert n1_local_idx in bridge_rows, (
            f"n=1 ground state (local index {n1_local_idx}) is not among the "
            f"bridge rows {sorted(set(bridge_rows.tolist()))[:10]}... "
            f"Bridge connectivity bug: outermost-shell-only selection."
        )

    def test_bonding_antibonding_splitting_exists(self):
        """Mean-field H2 must show bonding/antibonding splitting > 0.01 Ha.

        With only n=max_n outermost states bridged (old bug), the 1s core
        eigenfunctions were isolated on each atom and the mean-field gave
        no bonding/antibonding splitting at all (splitting ≈ 0).
        With n=1 states bridged, the splitting is significant (> 0.01 Ha).
        """
        mol = _build_h2(max_n=5)
        energies, _ = mol.compute_ground_state(n_states=2, method='mean_field')
        splitting = energies[1] - energies[0]
        assert splitting > 0.01, (
            f"Bonding/antibonding splitting = {splitting:.5f} Ha — too small. "
            f"Bridge connectivity bug causes near-zero splitting."
        )
