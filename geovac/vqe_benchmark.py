"""
VQE Benchmark — Head-to-head GeoVac vs Gaussian encoding comparison
====================================================================

Provides utilities for running VQE benchmarks comparing GeoVac lattice
encodings against Gaussian-basis (STO-3G) encodings on the same chemical
systems.  Collects metrics: qubit count, Pauli terms, CNOT count, circuit
depth, QWC measurement groups, VQE energy, and convergence iterations.

Requires: qiskit, qiskit-algorithms, openfermion

Author: GeoVac Development Team
Date: March 2026
"""

import json
import time
import warnings
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from openfermion import QubitOperator

from qiskit.quantum_info import SparsePauliOp
from qiskit.circuit import QuantumCircuit
from qiskit.circuit.library import efficient_su2
from qiskit.primitives import StatevectorEstimator
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import COBYLA

from geovac.qubit_encoding import JordanWignerEncoder
from geovac.gaussian_reference import he_sto3g, h2_sto3g, build_qubit_hamiltonian
from geovac.measurement_grouping import count_qwc_groups


# ---------------------------------------------------------------------------
# Conversion: OpenFermion QubitOperator -> Qiskit SparsePauliOp
# ---------------------------------------------------------------------------

def openfermion_to_sparse_pauli_op(
    qubit_op: QubitOperator,
    n_qubits: Optional[int] = None,
) -> SparsePauliOp:
    """
    Convert an OpenFermion QubitOperator to a Qiskit SparsePauliOp.

    Parameters
    ----------
    qubit_op : QubitOperator
        OpenFermion Pauli string Hamiltonian.
    n_qubits : int, optional
        Number of qubits.  If None, inferred from the operator.

    Returns
    -------
    SparsePauliOp
    """
    if n_qubits is None:
        max_q = -1
        for term in qubit_op.terms:
            for q, _ in term:
                if q > max_q:
                    max_q = q
        n_qubits = max_q + 1 if max_q >= 0 else 1

    pauli_labels: List[str] = []
    coeffs: List[complex] = []

    for term, coeff in qubit_op.terms.items():
        label = ['I'] * n_qubits
        for qubit_idx, pauli_char in term:
            label[qubit_idx] = pauli_char
        # Qiskit uses reversed qubit ordering (qubit 0 is rightmost)
        pauli_labels.append(''.join(reversed(label)))
        coeffs.append(complex(coeff))

    return SparsePauliOp(pauli_labels, coeffs).simplify()


def exact_ground_state_energy(spo: SparsePauliOp) -> float:
    """
    Compute exact ground state energy by diagonalizing the SparsePauliOp.

    Uses sparse eigensolver for > 12 qubits to avoid dense 2^Q memory.

    Parameters
    ----------
    spo : SparsePauliOp

    Returns
    -------
    float
        Lowest eigenvalue.
    """
    from scipy.sparse.linalg import eigsh

    n_qubits = spo.num_qubits
    if n_qubits > 12:
        mat = spo.to_matrix(sparse=True)
        eigvals, _ = eigsh(mat.tocsc(), k=1, which='SA')
        return float(np.real(eigvals[0]))

    mat = spo.to_matrix()
    if hasattr(mat, 'toarray'):
        mat = mat.toarray()
    eigvals = np.linalg.eigvalsh(mat)
    return float(np.real(eigvals[0]))


# ---------------------------------------------------------------------------
# Circuit metrics
# ---------------------------------------------------------------------------

def get_circuit_metrics(
    ansatz: QuantumCircuit,
    optimization_level: int = 3,
) -> Dict[str, int]:
    """
    Transpile a circuit and extract gate count metrics.

    Parameters
    ----------
    ansatz : QuantumCircuit
        The variational ansatz circuit.
    optimization_level : int
        Qiskit transpiler optimization level (0-3).

    Returns
    -------
    dict
        Keys: cx_count, depth, n_params.
    """
    decomposed = ansatz.decompose()
    # For large circuits (>14 qubits), level 3 can OOM in the Rust transpiler.
    # Fall back to level 1 for safety.
    eff_level = min(optimization_level, 1) if ansatz.num_qubits > 14 else optimization_level
    pm = generate_preset_pass_manager(optimization_level=eff_level)
    transpiled = pm.run(decomposed)
    ops = transpiled.count_ops()

    return {
        'cx_count': ops.get('cx', 0),
        'depth': transpiled.depth(),
        'n_params': ansatz.num_parameters,
    }


# ---------------------------------------------------------------------------
# VQE runner
# ---------------------------------------------------------------------------

@dataclass
class VQEResult:
    """Results from a single VQE run."""
    system: str
    encoding: str
    n_qubits: int
    n_pauli_terms: int
    cx_count: int
    circuit_depth: int
    n_qwc_groups: int
    vqe_energy: float
    exact_energy: float
    error_ha: float
    error_pct: float
    n_iterations: int
    wall_time_s: float

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def collect_static_metrics(
    spo: SparsePauliOp,
    of_qubit_op: QubitOperator,
    n_qubits: int,
    system_name: str,
    encoding_name: str,
    exact_energy: float,
    reps: int = 1,
) -> VQEResult:
    """
    Collect all metrics except VQE (for systems too large for statevector sim).

    Returns a VQEResult with vqe_energy = NaN.
    """
    ansatz = efficient_su2(n_qubits, reps=reps, entanglement='linear')
    metrics = get_circuit_metrics(ansatz, optimization_level=3)
    n_qwc = count_qwc_groups(of_qubit_op)

    return VQEResult(
        system=system_name,
        encoding=encoding_name,
        n_qubits=n_qubits,
        n_pauli_terms=len(spo),
        cx_count=metrics['cx_count'],
        circuit_depth=metrics['depth'],
        n_qwc_groups=n_qwc,
        vqe_energy=float('nan'),
        exact_energy=exact_energy,
        error_ha=float('nan'),
        error_pct=float('nan'),
        n_iterations=0,
        wall_time_s=0.0,
    )


def run_vqe(
    spo: SparsePauliOp,
    of_qubit_op: QubitOperator,
    n_qubits: int,
    system_name: str,
    encoding_name: str,
    exact_energy: float,
    reps: int = 1,
    maxiter: int = 500,
) -> VQEResult:
    """
    Run VQE on a SparsePauliOp and collect all metrics.

    Parameters
    ----------
    spo : SparsePauliOp
        Qiskit Hamiltonian.
    of_qubit_op : QubitOperator
        OpenFermion operator (for QWC grouping).
    n_qubits : int
        Number of qubits.
    system_name : str
        e.g. "He", "H2".
    encoding_name : str
        e.g. "GeoVac nmax=2", "STO-3G".
    exact_energy : float
        Exact ground state energy for error calculation.
    reps : int
        EfficientSU2 repetitions.
    maxiter : int
        Maximum optimizer iterations.

    Returns
    -------
    VQEResult
    """
    # Build ansatz
    ansatz = efficient_su2(n_qubits, reps=reps, entanglement='linear')

    # Circuit metrics
    metrics = get_circuit_metrics(ansatz, optimization_level=3)

    # QWC groups
    n_qwc = count_qwc_groups(of_qubit_op)

    # Iteration counter via callback
    iteration_count = [0]

    def callback(nfev, parameters, energy, stepsize):
        iteration_count[0] = nfev

    # Run VQE
    estimator = StatevectorEstimator()
    optimizer = COBYLA(maxiter=maxiter)

    vqe = VQE(
        estimator,
        ansatz,
        optimizer,
        callback=callback,
    )

    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        result = vqe.compute_minimum_eigenvalue(spo)
    wall_time = time.time() - t0

    vqe_energy = float(np.real(result.eigenvalue))
    error_ha = abs(vqe_energy - exact_energy)
    error_pct = 100.0 * error_ha / abs(exact_energy) if exact_energy != 0 else 0.0

    return VQEResult(
        system=system_name,
        encoding=encoding_name,
        n_qubits=n_qubits,
        n_pauli_terms=len(spo),
        cx_count=metrics['cx_count'],
        circuit_depth=metrics['depth'],
        n_qwc_groups=n_qwc,
        vqe_energy=vqe_energy,
        exact_energy=exact_energy,
        error_ha=error_ha,
        error_pct=error_pct,
        n_iterations=iteration_count[0],
        wall_time_s=round(wall_time, 2),
    )


# ---------------------------------------------------------------------------
# System builders
# ---------------------------------------------------------------------------

def build_geovac_he(max_n: int = 2) -> Tuple[SparsePauliOp, QubitOperator, int, float]:
    """
    Build GeoVac He qubit Hamiltonian.

    Returns (SparsePauliOp, QubitOperator, n_qubits, exact_energy).
    """
    from geovac.lattice_index import LatticeIndex

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        li = LatticeIndex(
            n_electrons=2, max_n=max_n, nuclear_charge=2,
            vee_method='slater_full', h1_method='hybrid',
        )

    enc = JordanWignerEncoder(li)
    of_op = enc.build_qubit_operator()
    n_qubits = li.n_sp
    spo = openfermion_to_sparse_pauli_op(of_op, n_qubits)

    # For small systems, diagonalize the qubit Hamiltonian directly.
    # For larger systems (>14 qubits), use the FCI solver instead —
    # Qiskit to_matrix() OOMs on 2^Q × 2^Q matrices for Q > ~20.
    if n_qubits <= 14:
        exact_e = exact_ground_state_energy(spo)
    else:
        energies, _ = li.compute_ground_state(n_states=1)
        exact_e = float(energies[0])

    return spo, of_op, n_qubits, exact_e


def build_gaussian_he() -> Tuple[SparsePauliOp, QubitOperator, int, float]:
    """
    Build Gaussian STO-3G He qubit Hamiltonian.

    Returns (SparsePauliOp, QubitOperator, n_qubits, exact_energy).
    """
    sys_data = he_sto3g()
    _, of_op, _ = build_qubit_hamiltonian(sys_data)
    n_qubits = 2  # 1 spatial orbital -> 2 spin-orbitals
    spo = openfermion_to_sparse_pauli_op(of_op, n_qubits)
    exact_e = exact_ground_state_energy(spo)
    return spo, of_op, n_qubits, exact_e


def build_gaussian_h2() -> Tuple[SparsePauliOp, QubitOperator, int, float]:
    """
    Build Gaussian STO-3G H2 qubit Hamiltonian at R=1.4 bohr.

    Returns (SparsePauliOp, QubitOperator, n_qubits, exact_energy).
    """
    sys_data = h2_sto3g()
    _, of_op, _ = build_qubit_hamiltonian(sys_data)
    n_qubits = 4  # 2 spatial orbitals -> 4 spin-orbitals
    spo = openfermion_to_sparse_pauli_op(of_op, n_qubits)
    exact_e = exact_ground_state_energy(spo)
    return spo, of_op, n_qubits, exact_e


# ---------------------------------------------------------------------------
# Comparison table formatting
# ---------------------------------------------------------------------------

def format_comparison_table(results: List[VQEResult]) -> str:
    """
    Format a list of VQEResults into a comparison table.

    Parameters
    ----------
    results : list of VQEResult

    Returns
    -------
    str
        Formatted table string.
    """
    header = (
        f"{'System':<12} {'Encoding':<18} {'Qubits':>6} {'Pauli':>7} "
        f"{'CX':>5} {'Depth':>6} {'QWC':>5} "
        f"{'VQE E (Ha)':>12} {'Exact E':>12} {'Err%':>8} {'Iters':>6}"
    )
    sep = '-' * len(header)

    lines = [
        "",
        "=" * len(header),
        "VQE Head-to-Head Benchmark: GeoVac vs Gaussian",
        "=" * len(header),
        header,
        sep,
    ]

    for r in results:
        import math
        vqe_str = f"{r.vqe_energy:>12.6f}" if not math.isnan(r.vqe_energy) else f"{'N/A':>12}"
        err_str = f"{r.error_pct:>7.3f}%" if not math.isnan(r.error_pct) else f"{'N/A':>8}"
        iter_str = f"{r.n_iterations:>6}" if r.n_iterations > 0 else f"{'N/A':>6}"
        lines.append(
            f"{r.system:<12} {r.encoding:<18} {r.n_qubits:>6} {r.n_pauli_terms:>7} "
            f"{r.cx_count:>5} {r.circuit_depth:>6} {r.n_qwc_groups:>5} "
            f"{vqe_str} {r.exact_energy:>12.6f} {err_str} {iter_str}"
        )

    lines.append(sep)
    lines.append("")

    return "\n".join(lines)


def save_results(
    results: List[VQEResult],
    path: str = "debug/data/vqe_benchmark_results.json",
) -> None:
    """Save VQE benchmark results to JSON."""
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    data = [r.to_dict() for r in results]
    with open(out_path, 'w') as f:
        json.dump(data, f, indent=2)
