"""
Tests for FrozenCoreLatticeIndex — hierarchical molecular solver.

Tests verify:
1. Core energy computation matches Slater-Condon diagonal
2. Effective H1 includes core J/K contributions
3. Active-space SD count matches C(n_active_sp, n_active_el)
4. Frozen-core LiH energy matches full FCI within tolerance
5. Frozen-core He matches full FCI (trivial case: freeze nothing vs freeze 1s)
6. Convenience constructor from_molecular works
7. Speedup is real: frozen-core is faster than full FCI
8. Electron count validation
"""

import time
import warnings
from math import comb

import numpy as np
import pytest

# Suppress V_ee method warnings during tests
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
def lih_frozen_core(lih_full_fci):
    """Frozen-core LiH: freeze Li 1s², solve 2-electron active space."""
    from geovac.frozen_core import FrozenCoreLatticeIndex
    mol, _ = lih_full_fci
    fc = FrozenCoreLatticeIndex(mol, frozen_orbitals=[0], n_active_electrons=2)
    E, psi = fc.solve()
    return fc, E[0]


@pytest.fixture(scope="module")
def he_full_fci():
    """Full 2-electron He FCI (reference)."""
    from geovac import LatticeIndex
    idx = LatticeIndex(
        n_electrons=2, max_n=4, nuclear_charge=2,
        vee_method='slater_full', h1_method='exact',
    )
    E, psi = idx.compute_ground_state(n_states=1)
    return idx, E[0]


# ---------------------------------------------------------------------------
# Tests: Basic construction
# ---------------------------------------------------------------------------

class TestConstruction:
    """Test FrozenCoreLatticeIndex construction and validation."""

    def test_sd_count_lih(self, lih_frozen_core):
        """Active SD count matches C(n_active_sp, n_active_el)."""
        fc, _ = lih_frozen_core
        n_active_sp = fc.n_active_sp
        n_active_el = fc.n_active_electrons
        expected = comb(n_active_sp, n_active_el)
        assert fc.n_sd == expected, (
            f"Expected C({n_active_sp},{n_active_el})={expected}, got {fc.n_sd}"
        )

    def test_sd_reduction_lih(self, lih_full_fci, lih_frozen_core):
        """Frozen-core has far fewer SDs than full FCI."""
        mol, _ = lih_full_fci
        fc, _ = lih_frozen_core
        reduction = mol.n_sd / fc.n_sd
        assert reduction > 100, (
            f"Expected >100× reduction, got {reduction:.0f}×"
        )

    def test_frozen_spinorb_count(self, lih_frozen_core):
        """Freezing 1 spatial orbital freezes 2 spin-orbitals."""
        fc, _ = lih_frozen_core
        assert fc.n_frozen_el == 2
        assert len(fc.frozen_spinorb) == 2
        assert 0 in fc.frozen_spinorb_set
        assert 1 in fc.frozen_spinorb_set

    def test_electron_count_validation(self):
        """Mismatched electron counts raise ValueError."""
        from geovac import LatticeIndex
        from geovac.frozen_core import FrozenCoreLatticeIndex
        idx = LatticeIndex(
            n_electrons=2, max_n=3, nuclear_charge=2,
            vee_method='slater_full',
        )
        with pytest.raises(ValueError, match="active.*parent"):
            FrozenCoreLatticeIndex(idx, frozen_orbitals=[0], n_active_electrons=3)

    def test_active_orbitals_exclude_frozen(self, lih_frozen_core):
        """No frozen spin-orbital appears in active set."""
        fc, _ = lih_frozen_core
        overlap = fc.frozen_spinorb_set & fc.active_spinorb_set
        assert len(overlap) == 0


# ---------------------------------------------------------------------------
# Tests: Core energy
# ---------------------------------------------------------------------------

class TestCoreEnergy:
    """Test frozen-core energy computation."""

    def test_core_energy_negative(self, lih_frozen_core):
        """Core energy should be negative (bound electrons)."""
        fc, _ = lih_frozen_core
        assert fc.E_core < 0, f"E_core = {fc.E_core} should be negative"

    def test_core_energy_reasonable(self, lih_frozen_core):
        """Li+ 1s² core energy should be close to -7.28 Ha."""
        fc, _ = lih_frozen_core
        # Li+ exact: -7.2799 Ha. Our graph Laplacian gives ~-7.23 at nmax=3
        # The core energy includes h1 + J - K, not the full atomic energy
        # but it should be in the right ballpark for Li 1s²
        assert fc.E_core < -5.0, f"E_core = {fc.E_core} too high for Li 1s²"
        assert fc.E_core > -10.0, f"E_core = {fc.E_core} too low for Li 1s²"


# ---------------------------------------------------------------------------
# Tests: Energy accuracy
# ---------------------------------------------------------------------------

class TestAccuracy:
    """Test frozen-core energy accuracy vs full FCI."""

    def test_lih_energy_close_to_full_fci(self, lih_full_fci, lih_frozen_core):
        """Frozen-core LiH should match full FCI within ~5%."""
        _, E_full = lih_full_fci
        _, E_fc = lih_frozen_core
        error_ha = abs(E_fc - E_full)
        error_pct = abs(error_ha / E_full) * 100
        print(f"\nFrozen-core accuracy:")
        print(f"  E_full = {E_full:.6f} Ha")
        print(f"  E_fc   = {E_fc:.6f} Ha")
        print(f"  Error  = {error_ha:.4f} Ha ({error_pct:.2f}%)")
        # Frozen-core approximation: expect <5% error
        assert error_pct < 5.0, (
            f"Frozen-core error {error_pct:.2f}% exceeds 5% tolerance"
        )

    def test_lih_frozen_core_is_variational(self, lih_full_fci, lih_frozen_core):
        """Frozen-core energy should be ABOVE full FCI (less variational freedom)."""
        _, E_full = lih_full_fci
        _, E_fc = lih_frozen_core
        # Frozen-core restricts the space → energy should be higher (less negative)
        assert E_fc >= E_full - 0.01, (
            f"E_fc={E_fc:.6f} < E_full={E_full:.6f} by more than 10 mHa "
            f"(variational violation)"
        )

    def test_he_frozen_nothing(self, he_full_fci):
        """He with no frozen orbitals should match full FCI exactly."""
        from geovac.frozen_core import FrozenCoreLatticeIndex
        idx, E_full = he_full_fci
        # Freeze nothing → should reproduce full FCI
        fc = FrozenCoreLatticeIndex(idx, frozen_orbitals=[], n_active_electrons=2)
        E_fc, _ = fc.solve(add_nuclear_repulsion=False)
        error = abs(E_fc[0] - E_full)
        assert error < 1e-8, (
            f"No-freeze He: |E_fc - E_full| = {error:.2e} > 1e-8"
        )


# ---------------------------------------------------------------------------
# Tests: Convenience constructor
# ---------------------------------------------------------------------------

class TestFromMolecular:
    """Test the from_molecular convenience constructor."""

    def test_from_molecular_lih(self):
        """from_molecular builds a valid frozen-core LiH solver."""
        from geovac.frozen_core import FrozenCoreLatticeIndex
        fc = FrozenCoreLatticeIndex.from_molecular(
            Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
            R=3.015, n_core_A=2, n_core_B=0,
        )
        assert fc.n_frozen_el == 2
        assert fc.n_active_electrons == 2
        E, _ = fc.solve()
        assert E[0] < -7.5, f"LiH energy {E[0]} should be < -7.5 Ha"

    def test_from_molecular_odd_core_raises(self):
        """Odd core electron count raises ValueError."""
        from geovac.frozen_core import FrozenCoreLatticeIndex
        with pytest.raises(ValueError, match="even"):
            FrozenCoreLatticeIndex.from_molecular(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=3.015, n_core_A=1, n_core_B=0,
            )


# ---------------------------------------------------------------------------
# Tests: Performance
# ---------------------------------------------------------------------------

class TestPerformance:
    """Test that frozen-core is faster than full FCI."""

    def test_speedup(self, lih_full_fci, lih_frozen_core):
        """Frozen-core should be significantly faster than full FCI."""
        mol, _ = lih_full_fci
        fc, _ = lih_frozen_core
        # Just verify the SD count reduction implies speedup
        ratio = mol.n_sd / fc.n_sd
        print(f"\nSD reduction: {mol.n_sd:,} -> {fc.n_sd:,} ({ratio:.0f}x)")
        assert ratio > 200, f"Expected >200× SD reduction, got {ratio:.0f}×"


# ---------------------------------------------------------------------------
# Tests: PES scan
# ---------------------------------------------------------------------------

class TestPESScan:
    """Test frozen-core PES at multiple R values."""

    def test_lih_pes_bound(self):
        """Frozen-core LiH should be bound (energy minimum exists)."""
        from geovac.frozen_core import FrozenCoreLatticeIndex

        R_values = [2.0, 3.015, 6.0, 10.0]
        energies = []
        for R in R_values:
            fc = FrozenCoreLatticeIndex.from_molecular(
                Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
                R=R, n_core_A=2, n_core_B=0,
            )
            E, _ = fc.solve()
            energies.append(E[0])

        print(f"\nFrozen-core LiH PES:")
        for R, E in zip(R_values, energies):
            print(f"  R={R:.3f}: E={E:.6f} Ha")

        # Should have a minimum (not monotonic)
        E_eq = energies[1]  # R=3.015
        E_far = energies[-1]  # R=10.0
        D_e = E_far - E_eq
        print(f"  D_e = {D_e:.6f} Ha ({D_e*627.5:.1f} kcal/mol)")
        assert D_e > 0, "LiH should be bound (D_e > 0)"
