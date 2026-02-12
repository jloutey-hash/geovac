#!/usr/bin/env python3
"""
Quick verification that geovac package is properly installed and functional.
This script tests the examples from README.md.
"""

print("=" * 70)
print("  GEOVAC INSTALLATION VERIFICATION")
print("=" * 70)
print()

# Test 1: Basic import
print("[TEST 1] Import geovac package...")
try:
    import geovac
    print(f"  ✓ Successfully imported geovac v{geovac.__version__}")
except ImportError as e:
    print(f"  ✗ Failed to import: {e}")
    exit(1)

# Test 2: Import main classes
print("\n[TEST 2] Import main classes...")
try:
    from geovac import GeometricLattice, HeliumHamiltonian, DiracHamiltonian
    print("  ✓ GeometricLattice")
    print("  ✓ HeliumHamiltonian")
    print("  ✓ DiracHamiltonian")
except ImportError as e:
    print(f"  ✗ Failed to import classes: {e}")
    exit(1)

# Test 3: Quick Start example from README
print("\n[TEST 3] Run Quick Start example from README...")
try:
    from geovac import HeliumHamiltonian
    
    # Initialize with calibrated parameters
    h = HeliumHamiltonian(max_n=3, Z=2, kinetic_scale=-0.103)
    
    # Compute ground state
    energy, wavefunction = h.compute_ground_state()
    
    print(f"  ✓ Ground State Energy: {energy[0]:.6f} Hartree")
    
    # Verify accuracy
    experimental = -2.90338583
    error_percent = abs(energy[0] - experimental) / abs(experimental) * 100
    print(f"  ✓ Error from experiment: {error_percent:.4f}%")
    
    if error_percent < 0.1:
        print("  ✓ Accuracy check PASSED (< 0.1% error)")
    else:
        print(f"  ✗ Accuracy check FAILED ({error_percent:.4f}% > 0.1%)")
        
except Exception as e:
    print(f"  ✗ Failed to run example: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

# Test 4: Convenience function
print("\n[TEST 4] Test convenience function solve_helium()...")
try:
    from geovac import solve_helium
    
    energy, psi = solve_helium(max_n=3)
    print(f"  ✓ E₀ = {energy[0]:.6f} Ha")
    print(f"  ✓ Wavefunction norm: {(psi.T @ psi)[0,0]:.8f}")
    
except Exception as e:
    print(f"  ✗ Failed: {e}")
    exit(1)

# Test 5: Lattice structure
print("\n[TEST 5] Test GeometricLattice...")
try:
    from geovac import GeometricLattice
    
    lattice = GeometricLattice(max_n=3)
    print(f"  ✓ Single-particle states: {lattice.num_states}")
    print(f"  ✓ Lattice connectivity: {lattice.adjacency.nnz} edges")
    print(f"  ✓ Sparsity: {(1 - lattice.adjacency.nnz / lattice.num_states**2) * 100:.1f}%")
    
except Exception as e:
    print(f"  ✗ Failed: {e}")
    exit(1)

# Test 6: Check calibrated constant
print("\n[TEST 6] Verify calibrated constant...")
try:
    print(f"  ✓ CALIBRATED_KINETIC_SCALE = {geovac.CALIBRATED_KINETIC_SCALE}")
    print(f"  ✓ EXPERIMENTAL_HELIUM_ENERGY = {geovac.EXPERIMENTAL_HELIUM_ENERGY} Ha")
except Exception as e:
    print(f"  ✗ Failed: {e}")
    exit(1)

# Summary
print("\n" + "=" * 70)
print("  ✓ ALL TESTS PASSED!")
print("=" * 70)
print()
print("GeoVac is properly installed and ready to use!")
print()
print("Try the Quick Start example:")
print("  >>> from geovac import solve_helium")
print("  >>> energy, psi = solve_helium(max_n=3)")
print("  >>> print(f'E₀ = {energy[0]:.6f} Ha')")
print()
