"""
Hydrogen Atom Convergence Debug
================================

CRITICAL: H atom should be EXACT with proper kinetic scale.
Current 4.9% error suggests fundamental issue.

This script diagnoses:
1. Basis set convergence
2. Graph Laplacian eigenvalues
3. Kinetic scale calibration
4. Comparison to analytic solution
"""

import numpy as np
import sys
import os
import io

# Set UTF-8 encoding for Windows
if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8')

# Add parent to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from geovac import GeometricLattice, HeliumHamiltonian, UNIVERSAL_KINETIC_SCALE
from scipy.sparse.linalg import eigsh
import time


def test_hydrogen_convergence():
    """Test H atom convergence with increasing basis size"""

    print("="*80)
    print("HYDROGEN ATOM CONVERGENCE TEST")
    print("="*80)
    print("\nAnalytic solution: E₀ = -Z²/(2n²) = -0.5 Ha (exact)")
    print("Universal kinetic scale: -1/16 = -0.0625")
    print("\n" + "="*80)

    print(f"\n{'max_n':>6s} {'States':>8s} {'Energy (Ha)':>14s} {'Error (%)':>12s} {'Time (ms)':>12s}")
    print("-"*80)

    results = []

    for max_n in [3, 5, 7, 10, 12, 15, 20, 25, 30]:
        try:
            # Build Hamiltonian for H atom (Z=1)
            t_start = time.time()
            h = HeliumHamiltonian(max_n=max_n, Z=1, kinetic_scale=UNIVERSAL_KINETIC_SCALE)

            # Get single-particle eigenvalue (H atom has 1 electron)
            eigenvalues, wavefunctions = eigsh(h.h1, k=1, which='SA')
            E = eigenvalues[0]
            t_elapsed = (time.time() - t_start) * 1000

            # Error vs exact
            E_exact = -0.5
            error = abs((E - E_exact) / E_exact) * 100

            n_states = h.lattice.num_states

            results.append({
                'max_n': max_n,
                'n_states': n_states,
                'E': E,
                'error': error,
                'time': t_elapsed
            })

            status = "✓" if error < 0.1 else "⚠"
            print(f"{status} {max_n:>4d} {n_states:>8d} {E:>14.8f} {error:>12.6f} {t_elapsed:>12.2f}")

        except Exception as e:
            print(f"✗ {max_n:>4d} FAILED: {str(e)}")

    print("-"*80)

    # Analyze convergence
    if len(results) > 3:
        print("\n" + "="*80)
        print("CONVERGENCE ANALYSIS")
        print("="*80)

        # Check if converging to -0.5
        final_E = results[-1]['E']
        final_error = results[-1]['error']

        print(f"\nFinal result (max_n={results[-1]['max_n']}):")
        print(f"  Energy:  {final_E:.8f} Ha")
        print(f"  Error:   {final_error:.6f}%")
        print(f"  Target:  -0.500000 Ha (exact)")

        if final_error < 0.01:
            print(f"\n  Status: ✓✓✓ EXCELLENT - Converged to machine precision!")
        elif final_error < 0.1:
            print(f"\n  Status: ✓✓ VERY GOOD - Nearly converged")
        elif final_error < 1.0:
            print(f"\n  Status: ✓ GOOD - Converging")
        else:
            print(f"\n  Status: ⚠ PROBLEM - Not converging properly")
            print(f"\n  DIAGNOSIS: Kinetic scale may be incorrect or basis incomplete")

        # Fit convergence
        if len(results) >= 5:
            print(f"\nExtrapolation to infinite basis:")
            max_ns = np.array([r['max_n'] for r in results[-5:]])
            energies = np.array([r['E'] for r in results[-5:]])

            # Fit E(n) = E_∞ + A/n^α
            from scipy.optimize import curve_fit

            def power_law(n, E_inf, A, alpha):
                return E_inf + A / n**alpha

            try:
                popt, pcov = curve_fit(power_law, max_ns, energies,
                                     p0=[-0.5, -1, 2], maxfev=10000)
                E_inf, A, alpha = popt

                print(f"  Fit: E(n) = E_∞ + A/n^α")
                print(f"  E_∞ = {E_inf:.8f} Ha")
                print(f"  α   = {alpha:.3f}")
                print(f"  Error at ∞: {abs((E_inf + 0.5)/0.5)*100:.6f}%")
            except:
                print("  (Fitting failed)")

    return results


def test_laplacian_eigenvalues():
    """Test graph Laplacian eigenvalues directly"""

    print("\n" + "="*80)
    print("GRAPH LAPLACIAN EIGENVALUE TEST")
    print("="*80)
    print("\nFor H atom, eigenvalues should match: E_n = -Z²/(2n²)")
    print("n=1: E₁ = -0.5 Ha")
    print("n=2: E₂ = -0.125 Ha")
    print("n=3: E₃ = -0.0556 Ha")
    print("\n" + "-"*80)

    max_n = 20
    h = HeliumHamiltonian(max_n=max_n, Z=1, kinetic_scale=UNIVERSAL_KINETIC_SCALE)

    # Get lowest 10 eigenvalues
    k_states = min(10, h.lattice.num_states - 1)
    eigenvalues, wavefunctions = eigsh(h.h1, k=k_states, which='SA')

    # Expected (hydrogenic) eigenvalues
    expected = [-0.5 / n**2 for n in range(1, k_states+1)]

    print(f"\n{'State':>6s} {'Computed':>14s} {'Expected':>14s} {'Error (%)':>12s} {'n_eff':>8s}")
    print("-"*80)

    for i, (E_computed, E_expected) in enumerate(zip(eigenvalues, expected)):
        error = abs((E_computed - E_expected) / E_expected) * 100

        # Effective quantum number: E = -0.5/n²
        n_eff = np.sqrt(-0.5 / E_computed) if E_computed < 0 else np.nan

        status = "✓" if error < 1.0 else "⚠"
        print(f"{status} {i+1:>4d} {E_computed:>14.8f} {E_expected:>14.8f} {error:>12.6f} {n_eff:>8.3f}")

    print("-"*80)
    print("\nInterpretation:")
    print("  If errors are small (<1%): Laplacian correctly represents kinetic energy")
    print("  If errors are large: Graph structure or weights incorrect")


def test_kinetic_scale_optimization():
    """Find optimal kinetic scale for H atom"""

    print("\n" + "="*80)
    print("KINETIC SCALE OPTIMIZATION")
    print("="*80)
    print("\nFinding optimal kinetic scale to minimize H atom error...")
    print(f"Current universal scale: {UNIVERSAL_KINETIC_SCALE:.8f}")
    print("\n" + "-"*80)

    max_n = 30
    scales_to_test = np.linspace(-0.08, -0.04, 41)

    best_scale = None
    best_error = float('inf')

    results = []

    for scale in scales_to_test:
        h = HeliumHamiltonian(max_n=max_n, Z=1, kinetic_scale=scale)
        eigenvalues, _ = eigsh(h.h1, k=1, which='SA')
        E = eigenvalues[0]
        error = abs(E + 0.5)

        results.append((scale, E, error))

        if error < best_error:
            best_error = error
            best_scale = scale

    print(f"\nOptimal kinetic scale: {best_scale:.8f}")
    print(f"Ground state energy:   {-0.5 - best_error:.8f} Ha")
    print(f"Error:                 {best_error*2*100:.6f}%")
    print(f"\nComparison:")
    print(f"  Universal (-1/16):  {UNIVERSAL_KINETIC_SCALE:.8f}")
    print(f"  H-optimized:        {best_scale:.8f}")
    print(f"  Difference:         {abs(best_scale - UNIVERSAL_KINETIC_SCALE):.8f}")

    if abs(best_scale - UNIVERSAL_KINETIC_SCALE) / abs(UNIVERSAL_KINETIC_SCALE) < 0.01:
        print(f"\n  Status: ✓ Universal scale is optimal (<1% difference)")
    else:
        percent_diff = abs(best_scale - UNIVERSAL_KINETIC_SCALE) / abs(UNIVERSAL_KINETIC_SCALE) * 100
        print(f"\n  Status: ⚠ Universal scale differs by {percent_diff:.2f}%")
        print(f"  → May need system-specific calibration")

    # Plot convergence near optimum
    print(f"\nError landscape near optimum:")
    print(f"{'Scale':>12s} {'Energy (Ha)':>14s} {'Error (Ha)':>14s}")
    print("-"*80)

    for scale, E, error in results[::5]:  # Every 5th point
        status = "✓" if abs(scale - best_scale) < 0.001 else " "
        print(f"{status} {scale:>10.8f} {E:>14.8f} {error:>14.8e}")


def test_lattice_structure():
    """Examine lattice construction"""

    print("\n" + "="*80)
    print("LATTICE STRUCTURE ANALYSIS")
    print("="*80)

    max_n = 10
    lattice = GeometricLattice(max_n=max_n)

    print(f"\nLattice parameters:")
    print(f"  max_n:       {max_n}")
    print(f"  States:      {lattice.num_states}")
    print(f"  Edges:       {lattice.num_edges}")
    print(f"  Sparsity:    {lattice.sparsity():.4f}")

    # Check state distribution by n
    print(f"\nStates by principal quantum number n:")
    print(f"{'n':>4s} {'States':>8s} {'Cumulative':>12s} {'Expected (n²)':>15s}")
    print("-"*80)

    cumulative = 0
    for n in range(1, max_n + 1):
        # For each n, have states with l=0,1,...,n-1 and m=-l,...,l
        n_states_for_n = n**2  # Sum over l: (2l+1) from l=0 to n-1
        cumulative += n_states_for_n

        print(f"{n:>4d} {n_states_for_n:>8d} {cumulative:>12d} {n**2:>15d}")

    print("-"*80)
    print(f"Total: {cumulative} states (should match {lattice.num_states})")

    if cumulative == lattice.num_states:
        print("✓ State count matches expected n² scaling")
    else:
        print(f"⚠ Mismatch: {lattice.num_states} vs {cumulative} expected")


def run_all_diagnostics():
    """Run complete diagnostic suite"""

    print("\n" + "#"*80)
    print("#" + " "*78 + "#")
    print("#" + "HYDROGEN ATOM DIAGNOSTIC SUITE".center(78) + "#")
    print("#" + " "*78 + "#")
    print("#"*80 + "\n")

    # Test 1: Convergence
    print("\nTEST 1: BASIS SET CONVERGENCE")
    results = test_hydrogen_convergence()

    # Test 2: Eigenvalues
    print("\n\nTEST 2: LAPLACIAN EIGENVALUES")
    test_laplacian_eigenvalues()

    # Test 3: Lattice structure
    print("\n\nTEST 3: LATTICE STRUCTURE")
    test_lattice_structure()

    # Test 4: Kinetic scale
    print("\n\nTEST 4: KINETIC SCALE OPTIMIZATION")
    test_kinetic_scale_optimization()

    # Final summary
    print("\n\n" + "="*80)
    print("DIAGNOSTIC SUMMARY")
    print("="*80)

    if results:
        final = results[-1]
        print(f"\nFinal H atom result (max_n={final['max_n']}):")
        print(f"  Energy:  {final['E']:.8f} Ha")
        print(f"  Target:  -0.500000 Ha")
        print(f"  Error:   {final['error']:.6f}%")

        if final['error'] < 0.1:
            print(f"\n✓✓✓ PASS: Hydrogen atom converges to exact solution!")
            print(f"    → GeoVac framework is fundamentally sound")
        elif final['error'] < 1.0:
            print(f"\n✓ ACCEPTABLE: Near convergence, may need larger basis")
        else:
            print(f"\n⚠ FAIL: Significant error remains")
            print(f"   → Investigate:")
            print(f"     1. Kinetic scale calibration")
            print(f"     2. Graph Laplacian construction")
            print(f"     3. Potential energy formulation")

    print("\n" + "#"*80 + "\n")


if __name__ == "__main__":
    run_all_diagnostics()
