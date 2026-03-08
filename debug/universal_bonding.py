"""
Universal Bonding Law — Dynamic Bridge Universality Test
========================================================

Tests whether the bridge amplitude A is a universal constant
or scales with atomic number Z.

Method:
    W_bridge = A * exp(-lambda * R)
    E_total = E_electronic(R) + V_NN(R)

Systems:
    H2:  Z=1+1, 2 electrons, R_eq = 1.40 Bohr, E_target = -1.174 Ha
    LiH: Z=3+1, 4 electrons, R_eq = 3.015 Bohr

Hypothesis:
    lambda_H >> lambda_Li  (hydrogen is tighter, faster decay)
    A = ???  (universal constant, or Z-dependent?)

Test Matrix:
    Test 1: H2 with A=8.5  (same as LiH calibration)
    Test 2: H2 with A=1.0  (standard overlap)
    Test 3: H2 lambda sweep (find optimal lambda at fixed A)
    Test 4: H2 A sweep     (find optimal A at fixed lambda)

Date: February 15, 2026
"""

import numpy as np
import sys
import time

sys.path.insert(0, '.')

from geovac import MoleculeHamiltonian


# ======================================================================
# Targets
# ======================================================================
H2_R_EQ = 1.40       # Bohr (experimental equilibrium)
H2_E_TARGET = -1.174  # Ha (experimental total energy at R_eq)
H2_BINDING = -0.17    # Ha (binding energy relative to 2H)


def scan_h2(bridge_amplitude: float = 8.5,
            bridge_decay_rate: float = 1.0,
            max_n: int = 5,
            R_min: float = 0.5,
            R_max: float = 4.0,
            R_step: float = 0.1) -> dict:
    """
    Scan H2 bond length with dynamic bridges.

    H2: 2 electrons, E_total = 2*E_0 + V_NN
    """
    distances = np.arange(R_min, R_max + R_step/2, R_step)
    results = []

    for R in distances:
        mol = MoleculeHamiltonian(
            nuclei=[(0.0, 0.0, 0.0), (R, 0.0, 0.0)],
            nuclear_charges=[1, 1],
            max_n=max_n,
            bridge_amplitude=bridge_amplitude,
            bridge_decay_rate=bridge_decay_rate,
        )

        E, _ = mol.compute_ground_state(n_states=2, method='mean_field')

        # 2-electron filling: E_total = 2*E_0 + V_NN
        E_electronic = 2 * E[0]
        V_NN = mol.compute_nuclear_repulsion()
        E_total = E_electronic + V_NN
        bridge_weight = mol.bridge_info[0]['bridge_weight'] if mol.bridge_info else 0.0

        results.append({
            'R': R,
            'E_electronic': E_electronic,
            'V_NN': V_NN,
            'E_total': E_total,
            'bridge_weight': bridge_weight,
        })

    # Find minimum
    E_min = min(r['E_total'] for r in results)
    R_eq = [r for r in results if r['E_total'] == E_min][0]['R']

    # Separated atom limit (largest R)
    E_inf = results[-1]['E_total']
    binding = E_min - E_inf

    return {
        'results': results,
        'R_eq': R_eq,
        'E_min': E_min,
        'E_inf': E_inf,
        'binding': binding,
        'A': bridge_amplitude,
        'lambda': bridge_decay_rate,
    }


def print_scan(scan: dict, label: str = "") -> None:
    """Print a bond length scan table."""
    results = scan['results']
    print(f"\n{'='*70}")
    print(f"H2 BOND SCAN: {label}")
    print(f"{'='*70}")
    print(f"  A = {scan['A']:.2f}, lambda = {scan['lambda']:.2f}")
    print(f"\n  {'R':>6}  {'W_bridge':>10}  {'E_elec':>10}  {'V_NN':>8}"
          f"  {'E_total':>10}  {'Note':>8}")
    print(f"  {'-'*6}  {'-'*10}  {'-'*10}  {'-'*8}  {'-'*10}  {'-'*8}")

    E_min = scan['E_min']
    for r in results:
        note = '<-- min' if r['E_total'] == E_min else ''
        print(f"  {r['R']:6.2f}  {r['bridge_weight']:10.6f}"
              f"  {r['E_electronic']:10.4f}  {r['V_NN']:8.4f}"
              f"  {r['E_total']:10.4f}  {note:>8}")

    print(f"\n  R_eq:     {scan['R_eq']:.2f} Bohr  (target: {H2_R_EQ:.2f})")
    print(f"  E_min:    {scan['E_min']:.6f} Ha  (target: {H2_E_TARGET:.3f})")
    print(f"  Binding:  {scan['binding']:.4f} Ha  (target: {H2_BINDING:.2f})")

    has_min = scan['R_eq'] != results[0]['R'] and scan['R_eq'] != results[-1]['R']
    print(f"  Status:   {'BOUND STATE' if has_min and scan['binding'] < 0 else 'NO BOND'}")


def test_1_universal_A():
    """Test 1: H2 with A=8.5 (LiH calibration). Is A universal?"""
    print("\n" + "#" * 70)
    print("TEST 1: H2 with A=8.5 (LiH Universal Constant?)")
    print("#" * 70)

    scan = scan_h2(bridge_amplitude=8.5, bridge_decay_rate=1.0, max_n=5)
    print_scan(scan, "A=8.5, lambda=1.0")
    return scan


def test_2_standard_A():
    """Test 2: H2 with A=1.0 (Standard Overlap?)."""
    print("\n" + "#" * 70)
    print("TEST 2: H2 with A=1.0 (Standard Overlap?)")
    print("#" * 70)

    scan = scan_h2(bridge_amplitude=1.0, bridge_decay_rate=1.0, max_n=5)
    print_scan(scan, "A=1.0, lambda=1.0")
    return scan


def test_3_lambda_sweep():
    """Test 3: Sweep lambda at fixed A=8.5 to find optimal decay rate."""
    print("\n" + "#" * 70)
    print("TEST 3: Lambda Sweep (A=8.5 fixed)")
    print("#" * 70)

    print(f"\n  {'lambda':>8}  {'R_eq':>8}  {'E_min':>10}  {'Binding':>10}  {'Status':>12}")
    print(f"  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*12}")

    best_R_err = 999
    best_lam = 0

    for lam in [0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]:
        scan = scan_h2(bridge_amplitude=8.5, bridge_decay_rate=lam, max_n=5)
        results = scan['results']
        has_min = scan['R_eq'] != results[0]['R'] and scan['R_eq'] != results[-1]['R']
        status = 'BOUND' if has_min and scan['binding'] < 0 else 'no bond'

        R_err = abs(scan['R_eq'] - H2_R_EQ)
        if has_min and R_err < best_R_err:
            best_R_err = R_err
            best_lam = lam

        print(f"  {lam:8.2f}  {scan['R_eq']:8.2f}  {scan['E_min']:10.4f}"
              f"  {scan['binding']:10.4f}  {status:>12}")

    print(f"\n  Best lambda: {best_lam:.2f} (R_eq closest to {H2_R_EQ})")
    return best_lam


def test_4_A_sweep():
    """Test 4: Sweep A at fixed lambda=1.0 to find optimal amplitude."""
    print("\n" + "#" * 70)
    print("TEST 4: Amplitude Sweep (lambda=1.0 fixed)")
    print("#" * 70)

    print(f"\n  {'A':>8}  {'R_eq':>8}  {'E_min':>10}  {'Binding':>10}  {'Status':>12}")
    print(f"  {'-'*8}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*12}")

    best_R_err = 999
    best_A = 0

    for A in [1.0, 2.0, 4.0, 6.0, 8.0, 8.5, 10.0, 15.0, 20.0, 30.0, 50.0]:
        scan = scan_h2(bridge_amplitude=A, bridge_decay_rate=1.0, max_n=5)
        results = scan['results']
        has_min = scan['R_eq'] != results[0]['R'] and scan['R_eq'] != results[-1]['R']
        status = 'BOUND' if has_min and scan['binding'] < 0 else 'no bond'

        R_err = abs(scan['R_eq'] - H2_R_EQ)
        if has_min and R_err < best_R_err:
            best_R_err = R_err
            best_A = A

        print(f"  {A:8.1f}  {scan['R_eq']:8.2f}  {scan['E_min']:10.4f}"
              f"  {scan['binding']:10.4f}  {status:>12}")

    print(f"\n  Best A: {best_A:.1f} (R_eq closest to {H2_R_EQ})")
    return best_A


def test_5_comparison():
    """Test 5: Side-by-side comparison of H2 and LiH optimal parameters."""
    print("\n" + "#" * 70)
    print("TEST 5: UNIVERSALITY COMPARISON")
    print("#" * 70)

    # H2 optimal (from Test 3/4 — we'll scan a grid)
    print("\n  Searching for H2 optimal (A, lambda)...")
    best_R_err = 999
    best_A = 0
    best_lam = 0
    best_binding = 0

    for A in [5.0, 8.0, 8.5, 10.0, 12.0, 15.0, 20.0]:
        for lam in [0.3, 0.5, 0.7, 1.0, 1.5, 2.0]:
            scan = scan_h2(bridge_amplitude=A, bridge_decay_rate=lam, max_n=5,
                           R_step=0.05)
            results = scan['results']
            has_min = scan['R_eq'] != results[0]['R'] and scan['R_eq'] != results[-1]['R']
            if not has_min:
                continue
            R_err = abs(scan['R_eq'] - H2_R_EQ)
            if R_err < best_R_err:
                best_R_err = R_err
                best_A = A
                best_lam = lam
                best_binding = scan['binding']
                best_E = scan['E_min']
                best_R = scan['R_eq']

    print(f"\n  {'':>20}  {'H2':>16}  {'LiH':>16}")
    print(f"  {'-'*20}  {'-'*16}  {'-'*16}")
    print(f"  {'A (amplitude)':>20}  {best_A:>16.1f}  {8.5:>16.1f}")
    print(f"  {'lambda (decay)':>20}  {best_lam:>16.2f}  {0.20:>16.2f}")
    print(f"  {'R_eq (Bohr)':>20}  {best_R:>16.2f}  {2.75:>16.2f}")
    print(f"  {'R_exp (Bohr)':>20}  {H2_R_EQ:>16.2f}  {3.015:>16.3f}")
    print(f"  {'Binding (Ha)':>20}  {best_binding:>16.4f}  {-0.19:>16.2f}")

    # Check universality
    print(f"\n  Universality Analysis:")
    if best_A == 8.5:
        print(f"  --> A = 8.5 for BOTH systems: A may be UNIVERSAL!")
    else:
        ratio = best_A / 8.5
        print(f"  --> A_H2 / A_LiH = {best_A:.1f} / 8.5 = {ratio:.2f}")
        if abs(ratio - 1.0) < 0.3:
            print(f"      Similar order of magnitude -> A may be approximately universal")
        else:
            print(f"      Different scales -> A likely depends on Z or system size")

    lam_ratio = best_lam / 0.2
    print(f"  --> lambda_H2 / lambda_LiH = {best_lam:.2f} / 0.20 = {lam_ratio:.1f}")
    print(f"      Hydrogen decays {lam_ratio:.0f}x faster (smaller atom, tighter wavefunction)")


if __name__ == '__main__':
    t0 = time.time()

    print("=" * 70)
    print("UNIVERSAL BONDING LAW: Dynamic Bridge Universality Test")
    print("=" * 70)
    print(f"\nQuestion: Is bridge amplitude A a universal constant?")
    print(f"Method:   W = A * exp(-lambda * R)")
    print(f"System:   H2 (Z=1+1, 2e-, R_eq=1.40 Bohr)")

    # Run all tests
    scan_1 = test_1_universal_A()
    scan_2 = test_2_standard_A()
    best_lambda = test_3_lambda_sweep()
    best_A = test_4_A_sweep()
    test_5_comparison()

    t_total = time.time() - t0
    print(f"\n{'='*70}")
    print(f"Total time: {t_total:.1f}s")
    print(f"{'='*70}")
