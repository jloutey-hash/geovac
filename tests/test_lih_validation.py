"""
LiH Validation Test — Dynamic Focal Length p₀(R) Benchmark
============================================================

v0.9.2: Tests the S³ chordal bridge framework with a heteronuclear
diatomic molecule. LiH is the ideal testbed because:

1. The Li-H bond is strongly heteronuclear (Z=3 vs Z=1), creating
   asymmetric conformal factors Ω_Li ≠ Ω_H at the bridge endpoints.
2. The dynamic focal length p₀(R) differs between the two atoms:
   - p₀_Li = sqrt(Z_Li² + V_nn/R) shifts with bond length
   - p₀_H  = sqrt(Z_H²  + V_nn/R) shifts more dramatically
3. Core states (n=1) for Li have p = Z/1 = 3, while H has p = 1,
   producing genuinely different conformal factors.

Experimental target:
    R_eq(LiH) = 1.5957 Å = 3.015 Bohr
    E(LiH)    ≈ -7.987 Ha (estimated Full CI limit)

Date: February 2026
"""

import numpy as np
import sys
import time

sys.path.insert(0, '.')

from geovac import (
    MoleculeHamiltonian,
    GeometricLattice,
    UNIVERSAL_KINETIC_SCALE,
    CALIBRATED_KINETIC_SCALE,
)


# ==============================================================================
# Reference values
# ==============================================================================

LIH_BOND_LENGTH_BOHR = 3.015   # Experimental R_eq (1.5957 Angstrom)
LIH_ENERGY_TARGET = -7.987     # Ha (estimated Full CI, Cade & Huo 1967)


def test_lih_bridge_r_dependence() -> dict:
    """
    Verify that the dynamic focal length p₀(R) produces R-dependent
    bridge weights, demonstrating that the conformal correction is
    active and non-trivial for the heteronuclear LiH system.

    This is the key validation: at each R, the bridge weight should
    differ from the flat-space STO value because p₀ ≠ Z.
    """
    print("\n" + "=" * 70)
    print("LiH BRIDGE R-DEPENDENCE — DYNAMIC p₀(R) VERIFICATION")
    print("=" * 70)

    max_n = 5
    Z_Li = 3
    Z_H = 1

    R_values = [2.0, 3.015, 4.0, 6.0, 10.0, 100.0]

    print(f"\n  {'R (Bohr)':>10}  {'p₀(Li)':>8}  {'p₀(H)':>8}  {'Ω_Li·Ω_H':>10}  {'W_base':>10}")
    print(f"  {'-'*10}  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*10}")

    t0 = time.time()
    results = []
    for R in R_values:
        atom_Li = GeometricLattice(max_n=max_n, nucleus_position=(0.0, 0.0, 0.0), nuclear_charge=Z_Li)
        atom_H = GeometricLattice(max_n=max_n, nucleus_position=(R, 0.0, 0.0), nuclear_charge=Z_H)
        mol = MoleculeHamiltonian(
            lattices=[atom_Li, atom_H],
            connectivity=[(0, 1, 2)],  # minimal bridges for diagnostics
            kinetic_scale=UNIVERSAL_KINETIC_SCALE,
            bridge_decay_rate=1.0,
        )
        bi = mol.bridge_info[0]
        p0_A = bi['p0_A']
        p0_B = bi['p0_B']

        # Compute Ω product for n=1 core states
        p_Li = Z_Li / 1.0
        p_H = Z_H / 1.0
        omega_Li = 2.0 * p0_A / (p_Li**2 + p0_A**2)
        omega_H = 2.0 * p0_B / (p_H**2 + p0_B**2)
        omega_prod = omega_Li * omega_H

        print(f"  {R:10.3f}  {p0_A:8.4f}  {p0_B:8.4f}  {omega_prod:10.6f}  {bi['bridge_weight']:10.6f}")
        results.append({'R': R, 'p0_A': p0_A, 'p0_B': p0_B, 'omega_prod': omega_prod})

    t1 = time.time()

    # Validation checks
    # 1. p₀ should decrease monotonically as R increases (V_nn shrinks)
    p0_A_vals = [r['p0_A'] for r in results]
    monotone_A = all(p0_A_vals[i] >= p0_A_vals[i+1] for i in range(len(p0_A_vals)-1))

    # 2. p₀(Li) should converge to Z_Li=3 at large R
    converge_Li = abs(results[-1]['p0_A'] - Z_Li) < 0.02

    # 3. p₀(H) should converge to Z_H=1 at large R
    converge_H = abs(results[-1]['p0_B'] - Z_H) < 0.02

    # 4. Ω product should be < 1.0 at equilibrium (conformal correction active)
    omega_eq = results[1]['omega_prod']  # R = 3.015
    correction_active = omega_eq < 0.95

    passed = monotone_A and converge_Li and converge_H and correction_active

    print(f"\n  Time: {(t1-t0)*1000:.0f} ms")
    print(f"\n  Validation:")
    print(f"    [{'PASS' if monotone_A else 'FAIL'}] p₀(Li) decreases monotonically with R")
    print(f"    [{'PASS' if converge_Li else 'FAIL'}] p₀(Li) → {Z_Li} as R → ∞  (actual: {results[-1]['p0_A']:.4f})")
    print(f"    [{'PASS' if converge_H else 'FAIL'}] p₀(H) → {Z_H} as R → ∞  (actual: {results[-1]['p0_B']:.4f})")
    print(f"    [{'PASS' if correction_active else 'FAIL'}] Ω_Li·Ω_H < 0.95 at R_eq  (actual: {omega_eq:.4f})")
    print(f"    Status: {'OK PASS' if passed else 'FAIL'}")

    return {
        'test': 'lih_bridge_r_dep',
        'results': results,
        'passed': passed,
    }


def test_lih_conformal_factors() -> dict:
    """
    Verify that the dynamic focal length produces asymmetric conformal
    factors for the heteronuclear Li-H bond, and that these factors
    vary with bond length R as expected from p₀²(R) = Z² + V_nn.
    """
    print("\n" + "=" * 70)
    print("LiH CONFORMAL FACTOR VERIFICATION")
    print("=" * 70)

    Z_Li = 3
    Z_H = 1
    max_n = 5

    R_test = LIH_BOND_LENGTH_BOHR
    V_nn = Z_Li * Z_H / R_test

    # Expected dynamic focal lengths
    p0_Li_expected = np.sqrt(Z_Li**2 + V_nn)
    p0_H_expected = np.sqrt(Z_H**2 + V_nn)

    # Conformal factors for n=1 core states
    p_Li = Z_Li / 1.0  # n=1
    p_H = Z_H / 1.0    # n=1
    omega_Li = 2.0 * p0_Li_expected / (p_Li**2 + p0_Li_expected**2)
    omega_H = 2.0 * p0_H_expected / (p_H**2 + p0_H_expected**2)

    print(f"\n  At R = {R_test:.3f} Bohr:")
    print(f"    V_nn = {V_nn:.4f} Ha")
    print(f"    p₀(Li) = sqrt({Z_Li}² + {V_nn:.4f}) = {p0_Li_expected:.4f}")
    print(f"    p₀(H)  = sqrt({Z_H}² + {V_nn:.4f}) = {p0_H_expected:.4f}")
    print(f"    Ω(Li, n=1) = {omega_Li:.4f}")
    print(f"    Ω(H,  n=1) = {omega_H:.4f}")
    print(f"    Ω_Li · Ω_H = {omega_Li * omega_H:.4f}")

    # Build molecule and check bridge info
    atom_Li = GeometricLattice(max_n=max_n, nucleus_position=(0.0, 0.0, 0.0), nuclear_charge=Z_Li)
    atom_H = GeometricLattice(max_n=max_n, nucleus_position=(R_test, 0.0, 0.0), nuclear_charge=Z_H)
    mol = MoleculeHamiltonian(
        lattices=[atom_Li, atom_H],
        connectivity=[(0, 1, 4 * max_n)],
        kinetic_scale=UNIVERSAL_KINETIC_SCALE,
        bridge_decay_rate=1.0,
    )

    bi = mol.bridge_info[0]
    p0_A_actual = bi.get('p0_A', None)
    p0_B_actual = bi.get('p0_B', None)

    print(f"\n  Computed by MoleculeHamiltonian:")
    print(f"    p₀_A = {p0_A_actual:.4f}  (expect {p0_Li_expected:.4f})")
    print(f"    p₀_B = {p0_B_actual:.4f}  (expect {p0_H_expected:.4f})")

    # Verify
    tol = 1e-6
    p0_A_ok = p0_A_actual is not None and abs(p0_A_actual - p0_Li_expected) < tol
    p0_B_ok = p0_B_actual is not None and abs(p0_B_actual - p0_H_expected) < tol
    asymmetric = p0_A_actual is not None and abs(p0_A_actual - p0_B_actual) > 0.1

    # Verify R-dependence: at R=10 Bohr, p0 should be closer to Z
    R_far = 10.0
    V_nn_far = Z_Li * Z_H / R_far
    p0_Li_far = np.sqrt(Z_Li**2 + V_nn_far)
    p0_H_far = np.sqrt(Z_H**2 + V_nn_far)

    atom_Li_far = GeometricLattice(max_n=max_n, nucleus_position=(0.0, 0.0, 0.0), nuclear_charge=Z_Li)
    atom_H_far = GeometricLattice(max_n=max_n, nucleus_position=(R_far, 0.0, 0.0), nuclear_charge=Z_H)
    mol_far = MoleculeHamiltonian(
        lattices=[atom_Li_far, atom_H_far],
        connectivity=[(0, 1, 4 * max_n)],
        kinetic_scale=UNIVERSAL_KINETIC_SCALE,
        bridge_decay_rate=1.0,
    )
    bi_far = mol_far.bridge_info[0]
    converging = abs(bi_far['p0_A'] - Z_Li) < abs(p0_A_actual - Z_Li)

    print(f"\n  R-dependence check (R={R_far} Bohr):")
    print(f"    p₀(Li) = {bi_far['p0_A']:.4f}  (should approach {Z_Li})")
    print(f"    p₀(H)  = {bi_far['p0_B']:.4f}  (should approach {Z_H})")

    passed = p0_A_ok and p0_B_ok and asymmetric and converging

    print(f"\n  Validation:")
    print(f"    [{'PASS' if p0_A_ok else 'FAIL'}] p₀(Li) matches analytic")
    print(f"    [{'PASS' if p0_B_ok else 'FAIL'}] p₀(H) matches analytic")
    print(f"    [{'PASS' if asymmetric else 'FAIL'}] Conformal factors are asymmetric")
    print(f"    [{'PASS' if converging else 'FAIL'}] p₀(R) → Z as R → ∞")
    print(f"    Status: {'OK PASS' if passed else 'FAIL'}")

    return {
        'test': 'lih_conformal',
        'p0_A': p0_A_actual,
        'p0_B': p0_B_actual,
        'passed': passed,
    }


# ==============================================================================
# MAIN
# ==============================================================================
if __name__ == '__main__':
    t_start = time.time()

    print("=" * 70)
    print("LiH VALIDATION SUITE — Dynamic Focal Length p₀(R)")
    print("=" * 70)

    r1 = test_lih_conformal_factors()
    r2 = test_lih_bridge_r_dependence()

    print(f"\n{'=' * 70}")
    print("LiH SUMMARY")
    print(f"{'=' * 70}")

    n_passed = sum(1 for r in [r1, r2] if r['passed'])
    print(f"\n  Result: {n_passed}/2 tests passed")
    print(f"  Total time: {time.time() - t_start:.1f}s")
    print(f"{'=' * 70}")
