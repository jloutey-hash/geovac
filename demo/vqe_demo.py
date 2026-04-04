"""
VQE Demo — Run a variational quantum eigensolver on a GeoVac Hamiltonian
=========================================================================

Demonstrates the full pipeline:
  1. Build a qubit Hamiltonian using GeoVac's ecosystem_export module
  2. Export to Qiskit SparsePauliOp
  3. Run VQE with EfficientSU2 ansatz on statevector simulator
  4. Compare to exact diagonalization

Requirements:
  pip install qiskit qiskit-algorithms openfermion

Usage:
  python demo/vqe_demo.py

Author: GeoVac Development Team
Date: April 2026
"""

import time
import warnings

import numpy as np
from scipy.optimize import minimize


def main() -> None:
    # ---------------------------------------------------------------
    # 1. Build Hamiltonian via ecosystem export
    # ---------------------------------------------------------------
    from geovac.ecosystem_export import hamiltonian

    print("Building H2 STO-3G Hamiltonian via ecosystem_export...")
    h = hamiltonian('H2')
    print(f"  System:      {h.metadata['system']}")
    print(f"  Basis:       {h.metadata['basis']}")
    print(f"  Qubits:      {h.n_qubits}")
    print(f"  Pauli terms: {h.n_terms}")
    print(f"  1-norm:      {h.one_norm:.4f}")

    # ---------------------------------------------------------------
    # 2. Export to Qiskit
    # ---------------------------------------------------------------
    spo = h.to_qiskit()
    print(f"\n  Exported to Qiskit SparsePauliOp ({spo.num_qubits} qubits)")

    # ---------------------------------------------------------------
    # 3. Exact diagonalization (reference)
    # ---------------------------------------------------------------
    mat = spo.to_matrix()
    if hasattr(mat, 'toarray'):
        mat = mat.toarray()
    exact_energy = float(np.real(np.linalg.eigvalsh(mat)[0]))
    print(f"\n  Exact ground state energy: {exact_energy:.6f} Ha")

    # ---------------------------------------------------------------
    # 4. Run VQE
    # ---------------------------------------------------------------
    from qiskit.circuit.library import efficient_su2
    from qiskit.quantum_info import Statevector

    n_qubits = h.n_qubits
    reps = 3
    ansatz = efficient_su2(n_qubits, reps=reps, entanglement='linear')
    n_params = ansatz.num_parameters

    print(f"\n  Ansatz: EfficientSU2 (reps={reps}, {n_params} parameters)")
    print("  Optimizer: COBYLA (5 random restarts, 1000 iterations each)")

    def evaluate(params: np.ndarray) -> float:
        """Statevector expectation value <psi|H|psi>."""
        bound = ansatz.assign_parameters(params)
        sv = Statevector(bound).data
        return float(np.real(sv.conj() @ mat @ sv))

    best_energy = np.inf
    np.random.seed(42)

    t0 = time.time()
    for restart in range(5):
        x0 = np.random.uniform(-np.pi, np.pi, n_params)
        result = minimize(evaluate, x0, method='COBYLA',
                          options={'maxiter': 1000, 'rhobeg': 0.5})
        if result.fun < best_energy:
            best_energy = result.fun
    wall_time = time.time() - t0

    error_mha = abs(best_energy - exact_energy) * 1000

    print(f"\n  VQE energy:   {best_energy:.6f} Ha")
    print(f"  Exact energy: {exact_energy:.6f} Ha")
    print(f"  Error:        {error_mha:.3f} mHa")
    print(f"  Wall time:    {wall_time:.1f} s")
    print(f"  Converged within 1 mHa (chemical accuracy): "
          f"{'YES' if error_mha < 1.0 else 'NO'}")

    # ---------------------------------------------------------------
    # 5. Also available: export to PennyLane and OpenFermion
    # ---------------------------------------------------------------
    print("\n  Other export formats available:")
    print(f"    h.to_openfermion() -> QubitOperator ({h.n_terms} terms)")
    try:
        pl_h = h.to_pennylane()
        print(f"    h.to_pennylane()   -> qml.Hamiltonian ({len(pl_h.coeffs)} terms)")
    except ImportError:
        print("    h.to_pennylane()   -> (pennylane not installed)")


if __name__ == '__main__':
    main()
