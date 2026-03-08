#!/usr/bin/env python3
"""
Quick verification that geovac package is properly installed and functional.
"""

print("=" * 70)
print("  GEOVAC INSTALLATION VERIFICATION")
print("=" * 70)
print()

# Test 1: Basic import
print("[TEST 1] Import geovac package...")
try:
    import geovac
    print(f"  OK Successfully imported geovac v{geovac.__version__}")
except ImportError as e:
    print(f"  FAIL Failed to import: {e}")
    exit(1)

# Test 2: Import main classes
print("\n[TEST 2] Import main classes...")
try:
    from geovac import GeometricLattice, MoleculeHamiltonian, DiracHamiltonian
    from geovac import AtomicSolver, LatticeIndex
    print("  OK GeometricLattice")
    print("  OK MoleculeHamiltonian")
    print("  OK DiracHamiltonian")
    print("  OK AtomicSolver")
    print("  OK LatticeIndex")
except ImportError as e:
    print(f"  FAIL Failed to import classes: {e}")
    exit(1)

# Test 3: AtomicSolver (single-electron)
print("\n[TEST 3] AtomicSolver hydrogen ground state...")
try:
    from geovac import AtomicSolver, UNIVERSAL_KINETIC_SCALE

    solver = AtomicSolver(max_n=10, Z=1, kinetic_scale=UNIVERSAL_KINETIC_SCALE)
    energies, wf = solver.compute_ground_state(n_states=1)
    E = energies[0]
    error_pct = abs(E + 0.5) / 0.5 * 100
    print(f"  OK H ground state: {E:.6f} Ha (error {error_pct:.2f}%)")
except Exception as e:
    print(f"  FAIL: {e}")
    exit(1)

# Test 4: LatticeIndex (2-electron FCI)
print("\n[TEST 4] LatticeIndex helium FCI...")
try:
    from geovac import LatticeIndex
    li = LatticeIndex(n_electrons=2, max_n=3, nuclear_charge=2, vee_method='chordal')
    eigvals, eigvecs = li.compute_ground_state(n_states=1)
    print(f"  OK He FCI ground state: {eigvals[0]:.6f} Ha")
except Exception as e:
    print(f"  FAIL: {e}")
    exit(1)

# Test 5: Lattice structure
print("\n[TEST 5] GeometricLattice...")
try:
    from geovac import GeometricLattice

    lattice = GeometricLattice(max_n=3)
    print(f"  OK Single-particle states: {lattice.num_states}")
    print(f"  OK Lattice connectivity: {lattice.adjacency.nnz} edges")
    print(f"  OK Sparsity: {(1 - lattice.adjacency.nnz / lattice.num_states**2) * 100:.1f}%")
except Exception as e:
    print(f"  FAIL: {e}")
    exit(1)

# Test 6: Check calibrated constant
print("\n[TEST 6] Verify calibrated constant...")
try:
    print(f"  OK UNIVERSAL_KINETIC_SCALE = {geovac.UNIVERSAL_KINETIC_SCALE}")
    print(f"  OK CALIBRATED_KINETIC_SCALE = {geovac.CALIBRATED_KINETIC_SCALE}")
except Exception as e:
    print(f"  FAIL: {e}")
    exit(1)

# Summary
print("\n" + "=" * 70)
print("  OK ALL TESTS PASSED!")
print("=" * 70)
print()
print("GeoVac is properly installed and ready to use!")
print()
print("Quick Start:")
print("  >>> from geovac import LatticeIndex")
print("  >>> li = LatticeIndex(n_electrons=2, max_n=3, nuclear_charge=2)")
print("  >>> eigvals, eigvecs = li.compute_ground_state()")
print()
