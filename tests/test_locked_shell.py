"""
Tests for LockedShellMolecule — closed-shell-aware molecular solver.

Tests verify:
1. Orbital classification (locked vs active)
2. SD count matches expected combinatorics
3. LiH energy within 5% of full FCI
4. Variational property (E_locked_shell >= E_full_fci)
5. Auto-detection of closed shells
6. PES produces bound molecule
7. Solve completes in < 1 second
"""

import time
import warnings

import numpy as np
import pytest

warnings.filterwarnings("ignore", message=".*LatticeIndex V_ee.*")
warnings.filterwarnings("ignore", message=".*MolecularLatticeIndex.*")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def lih_full_fci():
    """Full 4-electron LiH FCI at R=3.015 (reference)."""
    from geovac import MolecularLatticeIndex
    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
        R=3.015, n_electrons=4,
        vee_method='slater_full', fci_method='auto',
    )
    E, psi = mol.compute_ground_state(n_states=1)
    return mol, E[0]


@pytest.fixture(scope="module")
def lih_locked():
    """Locked-shell LiH: lock Li 1s, active_nmax=2."""
    from geovac.locked_shell import LockedShellMolecule
    mol = LockedShellMolecule(
        Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
        R=3.015, n_electrons=4,
        locked_config={0: [(1, 0)]},
        active_nmax=2,
    )
    E, psi = mol.solve()
    return mol, E[0]


# ---------------------------------------------------------------------------
# Tests: Classification
# ---------------------------------------------------------------------------

class TestClassification:
    """Test orbital classification into locked and active."""

    def test_locked_count(self, lih_locked):
        """Li 1s is 1 spatial orbital locked."""
        mol, _ = lih_locked
        assert len(mol._locked_spatial) == 1

    def test_active_spatial_count(self, lih_locked):
        """active_nmax=2 gives 9 active spatial orbitals (Li: 2s,2p=4; H: 1s,2s,2p=5)."""
        mol, _ = lih_locked
        assert len(mol._active_spatial) == 9

    def test_locked_electrons(self, lih_locked):
        """2 electrons locked in Li 1s."""
        mol, _ = lih_locked
        assert mol.n_locked_el == 2
        assert mol.n_active_el == 2

    def test_sd_count(self, lih_locked):
        """C(18, 2) = 153 active SDs."""
        from math import comb
        mol, _ = lih_locked
        expected = comb(mol.n_active_sp, mol.n_active_el)
        assert mol.n_sd == expected == 153


# ---------------------------------------------------------------------------
# Tests: Accuracy
# ---------------------------------------------------------------------------

class TestAccuracy:
    """Test energy accuracy vs full FCI."""

    def test_energy_within_5pct(self, lih_full_fci, lih_locked):
        """Locked-shell LiH within 5% of full FCI."""
        _, E_full = lih_full_fci
        _, E_ls = lih_locked
        error_pct = abs(E_ls - E_full) / abs(E_full) * 100
        print(f"\nLocked-shell accuracy: {error_pct:.2f}%")
        print(f"  Full FCI:     {E_full:.6f} Ha")
        print(f"  Locked-shell: {E_ls:.6f} Ha")
        assert error_pct < 5.0

    def test_variational(self, lih_full_fci, lih_locked):
        """Locked-shell energy >= full FCI (less variational freedom)."""
        _, E_full = lih_full_fci
        _, E_ls = lih_locked
        assert E_ls >= E_full - 0.01


# ---------------------------------------------------------------------------
# Tests: Performance
# ---------------------------------------------------------------------------

class TestPerformance:
    """Test that locked-shell is fast."""

    def test_solve_under_1_second(self):
        """LiH locked-shell solves in < 1 second."""
        from geovac.locked_shell import LockedShellMolecule
        t0 = time.perf_counter()
        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_nmax=2,
        )
        E, _ = mol.solve()
        dt = time.perf_counter() - t0
        print(f"\nSolve time: {dt:.3f}s")
        assert dt < 1.0, f"Solve took {dt:.2f}s, expected < 1.0s"

    def test_sd_count_under_200(self):
        """Active SDs should be small (~153)."""
        from geovac.locked_shell import LockedShellMolecule
        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_nmax=2,
        )
        assert mol.n_sd < 200


# ---------------------------------------------------------------------------
# Tests: Auto-detection
# ---------------------------------------------------------------------------

class TestAutoDetect:
    """Test automatic closed-shell detection."""

    def test_auto_detect_li_core(self):
        """Auto-detect should lock Li 1s for LiH."""
        from geovac.locked_shell import LockedShellMolecule
        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            active_nmax=2,  # auto-detect locked_config
        )
        # Should have locked Li 1s (atom 0, n=1, l=0)
        assert 0 in mol._locked_config
        assert (1, 0) in mol._locked_config[0]
        assert mol.n_locked_el == 2

    def test_auto_detect_h_no_core(self):
        """H atom should have no locked shells."""
        from geovac.locked_shell import LockedShellMolecule
        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            active_nmax=2,
        )
        assert 1 not in mol._locked_config


# ---------------------------------------------------------------------------
# Tests: PES
# ---------------------------------------------------------------------------

class TestPES:
    """Test PES properties."""

    def test_bound_molecule(self):
        """Locked-shell LiH should be bound."""
        from geovac.locked_shell import LockedShellMolecule
        energies = []
        for R in [2.0, 3.015, 10.0]:
            mol = LockedShellMolecule(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=R, n_electrons=4,
                locked_config={0: [(1, 0)]},
                active_nmax=2,
            )
            E, _ = mol.solve()
            energies.append(E[0])

        # Energy at short R should be lower than at infinity
        D_e = energies[-1] - min(energies)
        assert D_e > 0, f"Molecule not bound: D_e = {D_e}"
