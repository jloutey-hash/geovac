"""
Diagnostic: Cross-n overlap bridges for LiH.

Tests whether Wolfsberg-Helmholz coupling between Li 2s and H 1s
shifts R_eq toward the experimental value (3.015 bohr).

Date: 2026-03-20
"""
import numpy as np
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from geovac.lattice_index import MolecularLatticeIndex

R_eq_ref = 3.015  # bohr
D_e_ref_eV = 2.515  # eV


def test_coupling_elements():
    """Verify Li 2s <-> H 1s coupling is nonzero."""
    print("=" * 60)
    print("1. Cross-n coupling matrix elements at R=3.015")
    print("=" * 60)

    mol = MolecularLatticeIndex(
        Z_A=3, Z_B=1,
        nmax_A=3, nmax_B=3,
        R=3.015, n_electrons=4,
        vee_method='slater_full',
        use_dmatrix='hybrid',
        use_cross_n_bridges=True,
    )

    nA = mol._n_spatial_A
    H1 = mol._H1_spatial.toarray()

    print(f"\nnA = {nA} (Li spatial orbitals)")
    print(f"States A: {mol._li_A.lattice.states[:5]}")
    print(f"States B: {mol._li_B.lattice.states[:5]}")

    # Li states: (1,0,0)=idx0, (2,0,0)=idx1, (2,1,-1)=idx2, ...
    # H states: (1,0,0)=idx0, (2,0,0)=idx1, ...
    print(f"\nCross-atom H1 coupling elements:")
    print(f"  Li 1s (idx 0) <-> H 1s (idx {nA}): {H1[0, nA]:.6f} Ha  (same-n=1, D-matrix)")
    print(f"  Li 2s (idx 1) <-> H 1s (idx {nA}): {H1[1, nA]:.6f} Ha  (cross-n, SIGMA BOND)")
    print(f"  Li 2s (idx 1) <-> H 2s (idx {nA+1}): {H1[1, nA+1]:.6f} Ha  (same-n=2, D-matrix)")
    print(f"  Li 3s (idx 4) <-> H 1s (idx {nA}): {H1[4, nA]:.6f} Ha  (cross-n)")
    print(f"  Li 1s (idx 0) <-> H 2s (idx {nA+1}): {H1[0, nA+1]:.6f} Ha  (cross-n)")

    # Verify the key coupling is nonzero
    li2s_h1s = H1[1, nA]
    assert abs(li2s_h1s) > 1e-6, f"Li 2s <-> H 1s should be nonzero, got {li2s_h1s}"
    print(f"\n  SUCCESS: Li 2s <-> H 1s coupling = {li2s_h1s:.6f} Ha")


def test_comparison():
    """Compare energies with and without cross-n bridges."""
    print("\n" + "=" * 60)
    print("2. Energy comparison at R=3.015 bohr")
    print("=" * 60)

    R_test = 3.015
    common = dict(
        Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
        R=R_test, n_electrons=4,
        vee_method='slater_full',
        use_dmatrix='hybrid',
    )

    print("\n--- Without cross-n bridges ---")
    mol_no = MolecularLatticeIndex(**common, use_cross_n_bridges=False)
    E_no, _ = mol_no.compute_ground_state(n_states=1)

    print("\n--- With cross-n bridges ---")
    mol_yes = MolecularLatticeIndex(**common, use_cross_n_bridges=True)
    E_yes, _ = mol_yes.compute_ground_state(n_states=1)

    dE = (E_yes[0] - E_no[0]) * 1000  # mHa
    print(f"\nE(no bridges):   {E_no[0]:.6f} Ha")
    print(f"E(with bridges): {E_yes[0]:.6f} Ha")
    print(f"Delta E:         {dE:.2f} mHa")


def test_pes_scan():
    """PES scan to find R_eq with and without cross-n bridges."""
    print("\n" + "=" * 60)
    print("3. PES scan: R_eq and D_e")
    print("=" * 60)

    R_values = [2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0]
    common = dict(
        Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
        n_electrons=4,
        vee_method='slater_full',
        use_dmatrix='hybrid',
    )

    for label, cross_n in [("WITHOUT cross-n", False), ("WITH cross-n", True)]:
        print(f"\n--- {label} bridges ---")
        energies = []
        for R in R_values:
            mol = MolecularLatticeIndex(**common, R=R, use_cross_n_bridges=cross_n)
            E, _ = mol.compute_ground_state(n_states=1)
            energies.append(E[0])
            print(f"  R={R:5.2f}  E={E[0]:.6f} Ha")

        E_inf = energies[-1]
        idx_min = np.argmin(energies)
        R_min = R_values[idx_min]
        E_min = energies[idx_min]
        D_e_ha = E_inf - E_min
        D_e_eV = D_e_ha * 27.2114

        print(f"  R_eq ~ {R_min:.1f} bohr (ref: {R_eq_ref})")
        print(f"  D_e  ~ {D_e_eV:.3f} eV (ref: {D_e_ref_eV})")


def test_k_sensitivity():
    """Test sensitivity to the Wolfsberg-Helmholz K constant."""
    print("\n" + "=" * 60)
    print("4. K sensitivity at R=3.015 bohr")
    print("=" * 60)

    K_values = [1.0, 1.5, 1.75, 2.0, 2.5]
    common = dict(
        Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
        R=3.015, n_electrons=4,
        vee_method='slater_full',
        use_dmatrix='hybrid',
        use_cross_n_bridges=True,
    )

    print(f"\n{'K':>6}  {'E (Ha)':>12}  {'dE from K=0 (mHa)':>18}")
    print("-" * 42)

    # Baseline: no bridges
    mol_base = MolecularLatticeIndex(
        Z_A=3, Z_B=1, nmax_A=3, nmax_B=3,
        R=3.015, n_electrons=4,
        vee_method='slater_full',
        use_dmatrix='hybrid',
        use_cross_n_bridges=False,
    )
    E_base, _ = mol_base.compute_ground_state(n_states=1)
    print(f"{'0.00':>6}  {E_base[0]:>12.6f}  {'0.000':>18}")

    for K in K_values:
        mol = MolecularLatticeIndex(**common, cross_n_K=K)
        E, _ = mol.compute_ground_state(n_states=1)
        dE = (E[0] - E_base[0]) * 1000
        print(f"{K:>6.2f}  {E[0]:>12.6f}  {dE:>18.3f}")


if __name__ == '__main__':
    test_coupling_elements()
    test_comparison()
    test_pes_scan()
    test_k_sensitivity()
