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
8. Direct CI integration: matches matrix method exactly
9. Auto dispatch to Direct CI for large active spaces
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


# ---------------------------------------------------------------------------
# Tests: Direct CI integration
# ---------------------------------------------------------------------------

class TestDirectCI:
    """Test Direct CI integration with LockedShellMolecule."""

    def test_direct_matches_matrix(self):
        """Direct CI must match matrix method to < 1e-8 Ha on LiH."""
        from geovac.locked_shell import LockedShellMolecule
        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_nmax=2,
        )
        E_mat, _ = mol.solve(fci_method='matrix')
        E_dir, _ = mol.solve(fci_method='direct')
        diff = abs(E_mat[0] - E_dir[0])
        print(f"\nDirect vs matrix: diff = {diff:.2e} Ha")
        assert diff < 1e-8, f"Direct CI mismatch: {diff:.2e}"

    def test_adapter_remapping(self):
        """Adapter remaps to contiguous indices correctly."""
        from geovac.locked_shell import LockedShellMolecule
        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_nmax=2,
        )
        adapter = mol._build_direct_ci_adapter()

        # n_sp should be 2 * n_active_spatial
        assert adapter.n_sp == 2 * len(mol._active_spatial)
        # n_electrons should be active count
        assert adapter.n_electrons == mol.n_active_el
        # SD count preserved
        assert adapter.n_sd == mol.n_sd
        # All SD indices should be in range [0, n_sp)
        for sd in adapter.sd_basis:
            for s in sd:
                assert 0 <= s < adapter.n_sp

    def test_auto_dispatch(self):
        """fci_method='auto' uses matrix for small, direct for large."""
        from geovac.locked_shell import LockedShellMolecule
        # Small system: 153 SDs → matrix
        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_nmax=2,
        )
        # auto should pick 'matrix' for 153 SDs (< 5000)
        E_auto, _ = mol.solve(fci_method='auto')
        E_mat, _ = mol.solve(fci_method='matrix')
        assert abs(E_auto[0] - E_mat[0]) < 1e-10

    def test_lazy_eri_4d_direct(self):
        """Direct CI path must NOT build dense _eri_4d."""
        from geovac.locked_shell import LockedShellMolecule
        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_nmax=2,
        )
        assert mol._eri_4d is None, "Dense ERI should not be built during init"
        mol.solve(fci_method='direct')
        assert mol._eri_4d is None, "Direct CI should not trigger dense ERI"

    def test_matrix_builds_eri_4d(self):
        """Matrix path should lazily build dense _eri_4d."""
        from geovac.locked_shell import LockedShellMolecule
        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_nmax=2,
        )
        assert mol._eri_4d is None
        mol.solve(fci_method='matrix')
        assert mol._eri_4d is not None, "Matrix method should build dense ERI"


# ---------------------------------------------------------------------------
# Tests: ERI sparsity
# ---------------------------------------------------------------------------

class TestERISparsity:
    """Verify ERI density scales as O(1/M^2), not O(1)."""

    def test_eri_density_below_threshold(self):
        """Molecular ERI density should be well below 5%."""
        from geovac.locked_shell import LockedShellMolecule
        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_nmax=2,
        )
        n_eri = len(mol._parent._eri)
        m4 = mol.n_spatial ** 4
        density = n_eri / m4
        print(f"\nERI density: {n_eri}/{m4} = {density:.4%}")
        assert density < 0.05, f"ERI density {density:.2%} too high"

    def test_adapter_sparse_remap(self):
        """Adapter _eri should have fewer entries than O(M_active^4)."""
        from geovac.locked_shell import LockedShellMolecule
        mol = LockedShellMolecule(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_electrons=4,
            locked_config={0: [(1, 0)]},
            active_nmax=2,
        )
        adapter = mol._build_direct_ci_adapter()
        n_active_spatial = len(mol._active_spatial)
        m4_active = n_active_spatial ** 4
        n_adapter_eri = len(adapter._eri)
        density = n_adapter_eri / m4_active
        print(f"\nAdapter ERI: {n_adapter_eri}/{m4_active} = {density:.4%}")
        # Should be sparse — well below 10%
        assert density < 0.10
        # Should match parent entries in active block
        active_set = set(mol._active_spatial)
        expected = sum(
            1 for (a, b, c, d) in mol._parent._eri
            if a in active_set and b in active_set
            and c in active_set and d in active_set
        )
        assert n_adapter_eri == expected
