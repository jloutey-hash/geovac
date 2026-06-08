"""Diagnostic 3 — isolate the bug.

Rotate only h1, only eri, or both — see which produces the FCI shift on LiH.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as la

from geovac.molecular_spec import lih_spec
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy


def main():
    R = 3.015
    spec = lih_spec(R=R, max_n=2)
    result = build_balanced_hamiltonian(spec, R=R, cross_block_h1=True)
    h1 = result['h1']
    eri = result['eri']
    nuc_rep = result['nuclear_repulsion']
    M = result['M']

    fake = {'M': M, 'h1': h1, 'eri': eri, 'nuclear_repulsion': nuc_rep}
    fci_base = coupled_fci_energy(fake, n_electrons=4)['E_coupled']
    print(f"Baseline E_FCI = {fci_base:.10f}")

    # Rotation in Li_core block
    theta = 0.1
    R_mat = np.eye(M)
    R_mat[0, 0] = np.cos(theta)
    R_mat[0, 1] = -np.sin(theta)
    R_mat[1, 0] = np.sin(theta)
    R_mat[1, 1] = np.cos(theta)
    h1_r = R_mat @ h1 @ R_mat.T
    eri_r = np.einsum('pa,qb,rc,sd,abcd->pqrs', R_mat, R_mat, R_mat, R_mat, eri,
                      optimize='optimal')

    # Case A: rotate h1 only, leave eri alone
    fake_a = {'M': M, 'h1': h1_r, 'eri': eri, 'nuclear_repulsion': nuc_rep}
    fci_a = coupled_fci_energy(fake_a, n_electrons=4)['E_coupled']
    print(f"[A] h1 rotated, eri unchanged:  {fci_a:.10f}, "
          f"diff = {abs(fci_a - fci_base):.4e}")

    # Case B: leave h1, rotate eri only
    fake_b = {'M': M, 'h1': h1, 'eri': eri_r, 'nuclear_repulsion': nuc_rep}
    fci_b = coupled_fci_energy(fake_b, n_electrons=4)['E_coupled']
    print(f"[B] h1 unchanged, eri rotated:  {fci_b:.10f}, "
          f"diff = {abs(fci_b - fci_base):.4e}")

    # Case C: both
    fake_c = {'M': M, 'h1': h1_r, 'eri': eri_r, 'nuclear_repulsion': nuc_rep}
    fci_c = coupled_fci_energy(fake_c, n_electrons=4)['E_coupled']
    print(f"[C] both rotated (consistent): {fci_c:.10f}, "
          f"diff = {abs(fci_c - fci_base):.4e}")

    # Reality check: is the rotation actually orthogonal?
    rot_resid = float(np.max(np.abs(R_mat @ R_mat.T - np.eye(M))))
    print(f"\nR_mat orthogonality: {rot_resid:.2e}")

    # Tighter eri rotation check via different einsum order
    print("\nAlternative einsum for eri rotation:")
    eri_r2 = R_mat @ eri.reshape(M, -1)
    eri_r2 = eri_r2.reshape(M, M, M, M)
    eri_r2 = eri_r2.transpose(1, 0, 2, 3)
    eri_r2 = (R_mat @ eri_r2.reshape(M, -1)).reshape(M, M, M, M)
    eri_r2 = eri_r2.transpose(1, 0, 2, 3)
    # Now rotate axis 2 and 3
    eri_r2 = eri_r2.transpose(0, 1, 3, 2)
    eri_r2 = (R_mat @ eri_r2.reshape(-1, M).T).T.reshape(M, M, M, M)
    eri_r2 = eri_r2.transpose(0, 1, 3, 2)
    # Compare to einsum version
    diff = float(np.max(np.abs(eri_r - eri_r2)))
    print(f"  einsum vs sequential matmul diff: {diff:.4e}")


if __name__ == '__main__':
    main()
