"""Diagnostic 2: cross-check rotated FCI via JW + eigenvalue.

Test 3 of diag1 said: rotating Li_core orbitals (0, 1) breaks FCI by 1.2e-2.
Hypothesis: coupled_fci_energy has a bug. Compare to JW + sparse eigenvalue.
"""

from __future__ import annotations

import numpy as np
import scipy.linalg as la
from scipy.sparse.linalg import eigsh

from openfermion import jordan_wigner
from openfermion.linalg import get_sparse_operator

from geovac.molecular_spec import lih_spec, MolecularSpec, OrbitalBlock
from geovac.balanced_coupled import build_balanced_hamiltonian
from geovac.coupled_composition import coupled_fci_energy
from geovac.qubit_encoding import build_fermion_op_from_integrals


def jw_groundstate_energy(h1, eri, nuc_rep, M, n_electrons):
    """Build JW qubit op, get ground state eigenvalue."""
    fermion_op = build_fermion_op_from_integrals(h1, eri, nuc_rep)
    qubit_op = jordan_wigner(fermion_op)
    sparse = get_sparse_operator(qubit_op, n_qubits=2 * M)
    # Particle-number subspace projection: just find lowest eigvalue
    # in full Hilbert space (will be in N-electron sector if HF init is closest)
    eigs, _ = eigsh(sparse, k=1, which='SA')
    return float(eigs[0])


def make_h2_two_center_spec(R: float = 1.4, max_n: int = 2) -> MolecularSpec:
    nuclei = [
        {'Z': 1.0, 'position': (0.0, 0.0, -R/2), 'label': 'H_a'},
        {'Z': 1.0, 'position': (0.0, 0.0, +R/2), 'label': 'H_b'},
    ]
    blocks = [
        OrbitalBlock(label='H_a_atomic', block_type='lone_pair',
                     Z_center=1.0, n_electrons=1, max_n=max_n,
                     center_nucleus_idx=0),
        OrbitalBlock(label='H_b_atomic', block_type='lone_pair',
                     Z_center=1.0, n_electrons=1, max_n=max_n,
                     center_nucleus_idx=1),
    ]
    return MolecularSpec(name='H2', blocks=blocks,
                         nuclear_repulsion_constant=1.0/R,
                         description='H2 two-center', nuclei=nuclei, R=R)


def main():
    # Use H₂ for fast test (M=10, 2 electrons)
    R = 1.4
    spec = make_h2_two_center_spec(R=R, max_n=2)
    result = build_balanced_hamiltonian(spec, R=R, nuclei=spec.nuclei,
                                        cross_block_h1=True)
    h1 = result['h1']
    eri = result['eri']
    nuc_rep = result['nuclear_repulsion']
    M = result['M']
    n_el = 2
    print(f"H2: M = {M}, n_el = {n_el}")

    # Baseline
    fake = {'M': M, 'h1': h1, 'eri': eri, 'nuclear_repulsion': nuc_rep}
    fci_cc = coupled_fci_energy(fake, n_electrons=n_el)['E_coupled']
    print(f"\nBaseline (coupled_fci_energy):  {fci_cc:.10f}")
    fci_jw = jw_groundstate_energy(h1, eri, nuc_rep, M, n_el)
    print(f"Baseline (JW + sparse eig):     {fci_jw:.10f}")
    print(f"diff: {abs(fci_cc - fci_jw):.4e}")

    # Apply rotation in Li_core block (orbitals 0, 1)
    theta = 0.1
    R_mat = np.eye(M)
    R_mat[0, 0] = np.cos(theta)
    R_mat[0, 1] = -np.sin(theta)
    R_mat[1, 0] = np.sin(theta)
    R_mat[1, 1] = np.cos(theta)
    h1_r = R_mat @ h1 @ R_mat.T
    eri_r = np.einsum('pa,qb,rc,sd,abcd->pqrs', R_mat, R_mat, R_mat, R_mat, eri,
                      optimize='optimal')

    fake_r = {'M': M, 'h1': h1_r, 'eri': eri_r, 'nuclear_repulsion': nuc_rep}
    fci_r_cc = coupled_fci_energy(fake_r, n_electrons=n_el)['E_coupled']
    print(f"\nRotated (coupled_fci_energy):   {fci_r_cc:.10f}, "
          f"diff vs baseline = {abs(fci_r_cc - fci_cc):.4e}")
    fci_r_jw = jw_groundstate_energy(h1_r, eri_r, nuc_rep, M, n_el)
    print(f"Rotated (JW + sparse eig):      {fci_r_jw:.10f}, "
          f"diff vs baseline = {abs(fci_r_jw - fci_jw):.4e}")
    print(f"diff between methods: {abs(fci_r_cc - fci_r_jw):.4e}")


if __name__ == '__main__':
    main()
