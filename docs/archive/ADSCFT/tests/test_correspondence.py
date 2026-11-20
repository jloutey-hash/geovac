"""
Test AdS/CFT Correspondence

Validates the boundary-bulk translation and holographic duality.

Tests:
------
1. Paraboloid embedding validation
2. Symplectic plaquette calculations
3. Boundary-bulk consistency
4. State space correspondence
5. Coordinate mapping correctness

Author: GeoVac Development Team
Date: February 14, 2026
Status: Experimental
"""

import sys
import os
import numpy as np

# Add parent directories to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

from ADSCFT.boundary_to_bulk import BoundaryBulkTranslator
from ADSCFT.bulk.paraboloid_lattice import ParaboloidLattice
from ADSCFT.bulk.symplectic import (
    compute_plaquette_area,
    compute_shell_capacity,
    compute_triangle_area
)


def test_paraboloid_embedding():
    """Test that paraboloid embedding satisfies geometric properties"""
    print("\n" + "="*70)
    print("TEST 1: Paraboloid Embedding Validation")
    print("="*70)

    lattice = ParaboloidLattice(max_n=10)

    # Validate embedding
    errors = lattice.validate_embedding(verbose=True)

    # Check all errors are negligible
    for key, error in errors.items():
        assert error < 1e-12, f"{key} error too large: {error}"

    print("\n✓ PASS: Paraboloid embedding validated")


def test_symplectic_plaquettes():
    """Test symplectic plaquette area calculations"""
    print("\n" + "="*70)
    print("TEST 2: Symplectic Plaquette Calculations")
    print("="*70)

    lattice = ParaboloidLattice(max_n=10)

    # Test single plaquette area
    area = compute_plaquette_area(lattice, n=2, l=1, m=0)
    print(f"\nPlaquette area at (n=2, l=1, m=0): {area:.6f}")

    assert area is not None, "Plaquette area should not be None"
    assert area > 0, "Plaquette area must be positive"

    # Test shell capacity
    S_1 = compute_shell_capacity(lattice, n=1)
    print(f"Shell capacity S_1: {S_1:.6f}")

    # S_1 should be close to π/7.2 ≈ 0.436 (empirical from old research)
    expected_S_1 = 0.433
    rel_error = abs(S_1 - expected_S_1) / expected_S_1
    print(f"Relative error vs expected: {rel_error*100:.1f}%")

    assert S_1 > 0, "Shell capacity must be positive"
    assert rel_error < 0.1, f"S_1 too far from expected: {rel_error*100:.1f}% error"

    print("\n✓ PASS: Symplectic calculations validated")


def test_boundary_bulk_consistency():
    """Test that boundary and bulk have consistent state spaces"""
    print("\n" + "="*70)
    print("TEST 3: Boundary-Bulk State Space Consistency")
    print("="*70)

    translator = BoundaryBulkTranslator(max_n=10)

    # Validate correspondence
    results = translator.validate_correspondence(verbose=True)

    # All checks must pass
    for key, result in results.items():
        assert result, f"Correspondence check failed: {key}"

    print("\n✓ PASS: Boundary-bulk consistency verified")


def test_coordinate_mapping():
    """Test that quantum states map to correct 3D coordinates"""
    print("\n" + "="*70)
    print("TEST 4: Coordinate Mapping Correctness")
    print("="*70)

    translator = BoundaryBulkTranslator(max_n=5)

    # Test specific mappings
    test_cases = [
        # (n, l, m) → expected properties
        ((1, 0, 0), {'r': 1.0, 'z': -1.0}),  # Ground state
        ((2, 0, 0), {'r': 4.0, 'z': -0.25}),  # First excited s-state
        ((2, 1, 0), {'r': 4.0, 'z': -0.25}),  # First p-state
    ]

    print("\nState Mappings:")
    print("  (n,ℓ,m)  →  (x,y,z)  [r, z]")
    print("  " + "-"*60)

    for (n, l, m), expected in test_cases:
        coord = translator.get_coordinate_for_state(n, l, m)
        r_actual = np.sqrt(coord[0]**2 + coord[1]**2)
        z_actual = coord[2]

        print(f"  ({n},{l},{m:+d})   →  ({coord[0]:7.3f}, {coord[1]:7.3f}, {coord[2]:7.3f})  "
              f"r={r_actual:.3f}, z={z_actual:.3f}")

        # Check radial distance r = n²
        r_expected = expected['r']
        r_error = abs(r_actual - r_expected) / r_expected if r_expected > 0 else 0
        assert r_error < 1e-10, f"Radial distance error: {r_error}"

        # Check energy depth z = -1/n²
        z_expected = expected['z']
        z_error = abs(z_actual - z_expected) / abs(z_expected) if z_expected != 0 else 0
        assert z_error < 1e-10, f"Energy depth error: {z_error}"

    print("\n✓ PASS: Coordinate mapping validated")


def test_capacity_scaling():
    """Test that symplectic capacity scales correctly with n"""
    print("\n" + "="*70)
    print("TEST 5: Symplectic Capacity Scaling")
    print("="*70)

    translator = BoundaryBulkTranslator(max_n=10)

    # Compute capacities for several shells
    capacities = translator.compute_capacity_series(n_max=5)

    print("\nShell Capacities:")
    print("  n  →  S_n        S_n/n⁴")
    print("  " + "-"*40)

    scaling_factors = []
    for n, S_n in enumerate(capacities, start=1):
        scaling = S_n / (n ** 4)
        scaling_factors.append(scaling)
        print(f"  {n}  →  {S_n:8.6f}  {scaling:.6f}")

    # S_n should scale roughly as n⁴ (geometric area ~ radius²)
    # Scaling factor S_n/n⁴ should be roughly constant
    scaling_variance = np.std(scaling_factors) / np.mean(scaling_factors)
    print(f"\nScaling variance: {scaling_variance*100:.1f}%")

    # Allow 50% variance in scaling (geometric effects)
    assert scaling_variance < 0.5, f"Scaling too irregular: {scaling_variance*100:.1f}% variance"

    print("\n✓ PASS: Capacity scaling validated")


def test_triangle_area():
    """Test basic triangle area calculation"""
    print("\n" + "="*70)
    print("TEST 6: Triangle Area Calculation")
    print("="*70)

    # Test with known triangle (right triangle)
    p0 = np.array([0, 0, 0])
    p1 = np.array([3, 0, 0])
    p2 = np.array([0, 4, 0])

    area = compute_triangle_area(p0, p1, p2)
    expected_area = 0.5 * 3 * 4  # = 6.0

    print(f"\nRight triangle (3-4-5):")
    print(f"  Computed area: {area:.6f}")
    print(f"  Expected area: {expected_area:.6f}")

    rel_error = abs(area - expected_area) / expected_area
    assert rel_error < 1e-10, f"Triangle area error: {rel_error}"

    # Test with equilateral triangle (side length = 2)
    s = 2.0
    h = s * np.sqrt(3) / 2  # Height
    p0 = np.array([0, 0, 0])
    p1 = np.array([s, 0, 0])
    p2 = np.array([s/2, h, 0])

    area = compute_triangle_area(p0, p1, p2)
    expected_area = (s**2 * np.sqrt(3)) / 4  # = √3

    print(f"\nEquilateral triangle (side = {s}):")
    print(f"  Computed area: {area:.6f}")
    print(f"  Expected area: {expected_area:.6f}")

    rel_error = abs(area - expected_area) / expected_area
    assert rel_error < 1e-10, f"Triangle area error: {rel_error}"

    print("\n✓ PASS: Triangle area calculations validated")


def run_all_tests():
    """Run all AdS/CFT correspondence tests"""
    print("\n" + "="*70)
    print("AdS/CFT CORRESPONDENCE TEST SUITE")
    print("="*70)

    tests = [
        test_paraboloid_embedding,
        test_symplectic_plaquettes,
        test_boundary_bulk_consistency,
        test_coordinate_mapping,
        test_capacity_scaling,
        test_triangle_area,
    ]

    passed = 0
    failed = 0

    for test_func in tests:
        try:
            test_func()
            passed += 1
        except AssertionError as e:
            print(f"\n✗ FAIL: {e}")
            failed += 1
        except Exception as e:
            print(f"\n✗ ERROR: {e}")
            failed += 1

    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    print(f"Passed: {passed}/{len(tests)}")
    print(f"Failed: {failed}/{len(tests)}")

    if failed == 0:
        print("\n✓✓✓ ALL TESTS PASSED ✓✓✓")
    else:
        print(f"\n✗✗✗ {failed} TESTS FAILED ✗✗✗")

    print("="*70)

    return failed == 0


if __name__ == "__main__":
    import sys
    import io

    # Set UTF-8 encoding for Windows
    if sys.platform == 'win32':
        sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

    # Run all tests
    success = run_all_tests()

    # Exit with appropriate code
    sys.exit(0 if success else 1)
