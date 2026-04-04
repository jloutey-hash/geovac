#!/usr/bin/env python
"""
IBM Quantum VQE Demo -- GeoVac Hamiltonian on Aer simulator or IBM hardware
============================================================================

Demonstrates the full quantum computing pipeline:
  1. Build H2 qubit Hamiltonian via GeoVac ecosystem_export
  2. Export to Qiskit SparsePauliOp
  3. Run VQE with EfficientSU2 ansatz
  4. Compare to exact diagonalization

Modes:
  --simulator   Local Aer statevector simulation (DEFAULT, no account needed)
  --token TOKEN IBM Quantum hardware via qiskit-ibm-runtime Estimator V2

Optional:
  --shots N         Shot-based simulation (default: statevector / exact)
  --compare-gaussian Side-by-side with Gaussian STO-3G baseline
  --restarts N      Number of COBYLA restarts (default: 5)
  --maxiter N       Max COBYLA iterations per restart (default: 1000)

Usage:
  python demo/ibm_quantum_demo.py                        # statevector sim
  python demo/ibm_quantum_demo.py --shots 10000          # shot-based sim
  python demo/ibm_quantum_demo.py --compare-gaussian     # side-by-side
  python demo/ibm_quantum_demo.py --token YOUR_TOKEN     # IBM hardware

Requirements:
  pip install qiskit qiskit-aer openfermion
  (For hardware: pip install qiskit-ibm-runtime)

Author: GeoVac Development Team
Date: April 2026
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import warnings
from typing import Any, Dict, Optional, Tuple

import numpy as np
from scipy.optimize import minimize


# ---------------------------------------------------------------------------
# Hamiltonian builders
# ---------------------------------------------------------------------------

def build_geovac_h2() -> Tuple[Any, Dict[str, Any]]:
    """Build GeoVac H2 bond-pair Hamiltonian. Returns (SparsePauliOp, info)."""
    from geovac.ecosystem_export import hamiltonian

    h = hamiltonian('H2')
    spo = h.to_qiskit()
    info = {
        'label': 'GeoVac bond-pair',
        'n_qubits': h.n_qubits,
        'n_terms': h.n_terms,
        'one_norm': h.one_norm,
        'metadata': h.metadata,
    }

    # QWC groups
    try:
        from geovac.measurement_grouping import analyze_measurement_cost
        analysis = analyze_measurement_cost(h.to_openfermion())
        info['n_qwc_groups'] = analysis.n_qwc_groups
    except Exception:
        info['n_qwc_groups'] = None

    return spo, info


def build_gaussian_h2() -> Tuple[Any, Dict[str, Any]]:
    """Build Gaussian STO-3G H2 Hamiltonian. Returns (SparsePauliOp, info)."""
    from geovac.gaussian_reference import h2_sto3g, build_qubit_hamiltonian
    from geovac.ecosystem_export import GeoVacHamiltonian

    sys_data = h2_sto3g()
    _, qubit_op, _ = build_qubit_hamiltonian(sys_data)
    h = GeoVacHamiltonian(qubit_op, metadata={
        'system': 'H2', 'basis': 'STO-3G', 'R_bohr': 1.4, 'Q': 4,
    })
    spo = h.to_qiskit()
    info = {
        'label': 'Gaussian STO-3G',
        'n_qubits': h.n_qubits,
        'n_terms': h.n_terms,
        'one_norm': h.one_norm,
        'metadata': h.metadata,
    }

    try:
        from geovac.measurement_grouping import analyze_measurement_cost
        analysis = analyze_measurement_cost(h.to_openfermion())
        info['n_qwc_groups'] = analysis.n_qwc_groups
    except Exception:
        info['n_qwc_groups'] = None

    return spo, info


# ---------------------------------------------------------------------------
# Exact diagonalization
# ---------------------------------------------------------------------------

def exact_ground_state(spo: Any) -> float:
    """Compute exact ground state energy via full diagonalization."""
    mat = spo.to_matrix()
    if hasattr(mat, 'toarray'):
        mat = mat.toarray()
    eigvals = np.linalg.eigvalsh(np.array(mat, dtype=complex))
    return float(np.real(eigvals[0]))


# ---------------------------------------------------------------------------
# VQE: statevector (scipy COBYLA)
# ---------------------------------------------------------------------------

def run_vqe_statevector(
    spo: Any,
    n_restarts: int = 5,
    maxiter: int = 1000,
    seed: int = 42,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Run VQE using statevector simulation + scipy COBYLA.

    This is the fastest local approach: no shots, no sampling noise.
    """
    from qiskit.circuit.library import efficient_su2
    from qiskit.quantum_info import Statevector

    n_qubits = spo.num_qubits
    ansatz = efficient_su2(n_qubits, reps=3, entanglement='linear')
    n_params = ansatz.num_parameters

    # Build dense matrix for fast expectation values
    mat = spo.to_matrix()
    if hasattr(mat, 'toarray'):
        mat = mat.toarray()
    mat = np.array(mat, dtype=complex)

    def evaluate(params: np.ndarray) -> float:
        bound = ansatz.assign_parameters(params)
        sv = Statevector(bound).data
        return float(np.real(sv.conj() @ mat @ sv))

    if verbose:
        print(f"  Ansatz: EfficientSU2 (reps=3, {n_params} parameters)")
        print(f"  Optimizer: COBYLA ({n_restarts} restarts, {maxiter} iter each)")
        print(f"  Method: statevector (exact expectation values)")

    best_energy = np.inf
    best_params = None
    rng = np.random.RandomState(seed)

    t0 = time.time()
    for restart in range(n_restarts):
        x0 = rng.uniform(-np.pi, np.pi, n_params)
        result = minimize(evaluate, x0, method='COBYLA',
                          options={'maxiter': maxiter, 'rhobeg': 0.5})
        if result.fun < best_energy:
            best_energy = result.fun
            best_params = result.x
    wall_time = time.time() - t0

    return {
        'energy': best_energy,
        'params': best_params,
        'wall_time': wall_time,
        'method': 'statevector',
        'n_restarts': n_restarts,
        'n_params': n_params,
    }


# ---------------------------------------------------------------------------
# VQE: shot-based (Aer sampling via StatevectorEstimator with precision)
# ---------------------------------------------------------------------------

def run_vqe_shots(
    spo: Any,
    shots: int = 10000,
    n_restarts: int = 5,
    maxiter: int = 1000,
    seed: int = 42,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Run VQE with shot-based estimation using Qiskit StatevectorEstimator.

    Uses finite precision to simulate shot noise. The precision is set as
    1/sqrt(shots) to approximate the statistical error from N shots.
    """
    from qiskit.circuit.library import efficient_su2
    from qiskit.quantum_info import Statevector

    n_qubits = spo.num_qubits
    ansatz = efficient_su2(n_qubits, reps=3, entanglement='linear')
    n_params = ansatz.num_parameters

    # For shot-based simulation, we add Gaussian noise to the exact
    # expectation value to simulate finite sampling. This is faster than
    # running a full measurement circuit simulation for 10-qubit systems.
    mat = spo.to_matrix()
    if hasattr(mat, 'toarray'):
        mat = mat.toarray()
    mat = np.array(mat, dtype=complex)

    # Shot noise standard deviation: 1-norm / sqrt(shots) is an upper bound
    one_norm = float(np.sum(np.abs(np.array([c for c in spo.coeffs]))))
    noise_std = one_norm / np.sqrt(shots)
    rng_noise = np.random.RandomState(seed + 1000)

    def evaluate(params: np.ndarray) -> float:
        bound = ansatz.assign_parameters(params)
        sv = Statevector(bound).data
        exact = float(np.real(sv.conj() @ mat @ sv))
        noise = rng_noise.normal(0, noise_std)
        return exact + noise

    if verbose:
        print(f"  Ansatz: EfficientSU2 (reps=3, {n_params} parameters)")
        print(f"  Optimizer: COBYLA ({n_restarts} restarts, {maxiter} iter each)")
        print(f"  Method: shot-based ({shots:,} shots, noise_std={noise_std:.4f})")

    best_energy = np.inf
    best_params = None
    rng = np.random.RandomState(seed)

    t0 = time.time()
    for restart in range(n_restarts):
        x0 = rng.uniform(-np.pi, np.pi, n_params)
        result = minimize(evaluate, x0, method='COBYLA',
                          options={'maxiter': maxiter, 'rhobeg': 0.5})
        if result.fun < best_energy:
            best_energy = result.fun
            best_params = result.x
    wall_time = time.time() - t0

    return {
        'energy': best_energy,
        'params': best_params,
        'wall_time': wall_time,
        'method': f'shot-based ({shots:,} shots)',
        'shots': shots,
        'noise_std': noise_std,
        'n_restarts': n_restarts,
        'n_params': n_params,
    }


# ---------------------------------------------------------------------------
# VQE: IBM Quantum hardware
# ---------------------------------------------------------------------------

def run_vqe_hardware(
    spo: Any,
    token: str,
    n_restarts: int = 3,
    maxiter: int = 200,
    verbose: bool = True,
) -> Dict[str, Any]:
    """
    Run VQE on IBM Quantum hardware via qiskit-ibm-runtime Estimator V2.

    WARNING: This submits real jobs to IBM Quantum. Each function evaluation
    is a remote job, so VQE with many iterations can consume significant
    queue time and credits. Use --simulator for development and testing.
    """
    try:
        from qiskit_ibm_runtime import QiskitRuntimeService, EstimatorV2
    except ImportError:
        raise ImportError(
            "qiskit-ibm-runtime is required for hardware mode. "
            "Install with: pip install qiskit-ibm-runtime"
        )

    from qiskit.circuit.library import efficient_su2
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

    n_qubits = spo.num_qubits

    # Connect to IBM Quantum
    if verbose:
        print(f"  Connecting to IBM Quantum...")

    service = QiskitRuntimeService(
        channel='ibm_quantum',
        token=token,
    )

    # Find least-busy backend with enough qubits
    backends = service.backends(
        min_num_qubits=n_qubits,
        simulator=False,
        operational=True,
    )
    if not backends:
        raise RuntimeError(
            f"No available IBM Quantum backends with >= {n_qubits} qubits. "
            f"The GeoVac H2 encoding requires {n_qubits} qubits."
        )

    backend = service.least_busy(
        min_num_qubits=n_qubits,
        simulator=False,
        operational=True,
    )
    if verbose:
        print(f"  Backend: {backend.name} ({backend.num_qubits} qubits)")

    # Build ansatz
    ansatz = efficient_su2(n_qubits, reps=3, entanglement='linear')
    n_params = ansatz.num_parameters

    # Transpile for target backend
    pm = generate_preset_pass_manager(
        optimization_level=2,
        backend=backend,
    )
    isa_circuit = pm.run(ansatz)

    if verbose:
        print(f"  Ansatz: EfficientSU2 (reps=3, {n_params} parameters)")
        print(f"  Optimizer: COBYLA ({n_restarts} restarts, {maxiter} iter each)")
        print(f"  Method: IBM Quantum hardware")
        print(f"  NOTE: Each VQE iteration submits a remote job. This may be slow.")

    # Transform observable for ISA
    isa_observable = spo.apply_layout(isa_circuit.layout)

    # Create Estimator V2
    estimator = EstimatorV2(mode=backend)

    def evaluate(params: np.ndarray) -> float:
        # Bind parameters to the ISA circuit
        bound = isa_circuit.assign_parameters(dict(zip(ansatz.parameters, params)))
        job = estimator.run([(bound, isa_observable)])
        result = job.result()
        return float(result[0].data.evs)

    best_energy = np.inf
    best_params = None
    rng = np.random.RandomState(42)

    t0 = time.time()
    for restart in range(n_restarts):
        if verbose:
            print(f"  Restart {restart + 1}/{n_restarts}...")
        x0 = rng.uniform(-np.pi, np.pi, n_params)
        result = minimize(evaluate, x0, method='COBYLA',
                          options={'maxiter': maxiter, 'rhobeg': 0.5})
        if result.fun < best_energy:
            best_energy = result.fun
            best_params = result.x
    wall_time = time.time() - t0

    return {
        'energy': best_energy,
        'params': best_params,
        'wall_time': wall_time,
        'method': f'IBM Quantum ({backend.name})',
        'backend': backend.name,
        'n_restarts': n_restarts,
        'n_params': n_params,
    }


# ---------------------------------------------------------------------------
# Results formatting
# ---------------------------------------------------------------------------

def print_results(
    label: str,
    info: Dict[str, Any],
    vqe_result: Dict[str, Any],
    exact_energy: float,
) -> str:
    """Print and return formatted VQE results."""
    error_mha = abs(vqe_result['energy'] - exact_energy) * 1000
    chem_acc = error_mha < 1.0

    lines = [
        f"",
        f"  === {label} ===",
        f"  Qubits:       {info['n_qubits']}",
        f"  Pauli terms:  {info['n_terms']}",
        f"  1-norm:       {info['one_norm']:.4f} Ha",
    ]
    if info.get('n_qwc_groups') is not None:
        lines.append(f"  QWC groups:   {info['n_qwc_groups']}")
    lines += [
        f"  VQE method:   {vqe_result['method']}",
        f"  VQE energy:   {vqe_result['energy']:.6f} Ha",
        f"  Exact energy: {exact_energy:.6f} Ha",
        f"  Error:        {error_mha:.3f} mHa",
        f"  Wall time:    {vqe_result['wall_time']:.1f} s",
        f"  Chemical accuracy (<1 mHa): {'YES' if chem_acc else 'NO'}",
    ]
    text = '\n'.join(lines)
    print(text)
    return text


def print_comparison(
    geovac_info: Dict[str, Any],
    geovac_vqe: Dict[str, Any],
    geovac_exact: float,
    gaussian_info: Dict[str, Any],
    gaussian_vqe: Dict[str, Any],
    gaussian_exact: float,
) -> str:
    """Print side-by-side comparison table."""
    g_err = abs(geovac_vqe['energy'] - geovac_exact) * 1000
    s_err = abs(gaussian_vqe['energy'] - gaussian_exact) * 1000

    lines = [
        "",
        "  " + "=" * 60,
        f"  {'Metric':<25} {'GeoVac':>15} {'Gaussian STO-3G':>15}",
        "  " + "-" * 60,
        f"  {'Qubits':<25} {geovac_info['n_qubits']:>15} {gaussian_info['n_qubits']:>15}",
        f"  {'Pauli terms':<25} {geovac_info['n_terms']:>15} {gaussian_info['n_terms']:>15}",
        f"  {'1-norm (Ha)':<25} {geovac_info['one_norm']:>15.4f} {gaussian_info['one_norm']:>15.4f}",
    ]
    if geovac_info.get('n_qwc_groups') and gaussian_info.get('n_qwc_groups'):
        lines.append(
            f"  {'QWC groups':<25} {geovac_info['n_qwc_groups']:>15} {gaussian_info['n_qwc_groups']:>15}"
        )
    lines += [
        f"  {'VQE energy (Ha)':<25} {geovac_vqe['energy']:>15.6f} {gaussian_vqe['energy']:>15.6f}",
        f"  {'Exact energy (Ha)':<25} {geovac_exact:>15.6f} {gaussian_exact:>15.6f}",
        f"  {'VQE error (mHa)':<25} {g_err:>15.3f} {s_err:>15.3f}",
        f"  {'Wall time (s)':<25} {geovac_vqe['wall_time']:>15.1f} {gaussian_vqe['wall_time']:>15.1f}",
        "  " + "=" * 60,
        "",
        f"  Note: GeoVac uses {geovac_info['n_qubits']} qubits (bond-pair encoding) vs",
        f"  Gaussian STO-3G using {gaussian_info['n_qubits']} qubits. The GeoVac basis captures",
        f"  more correlation physics at the cost of more qubits.",
    ]
    text = '\n'.join(lines)
    print(text)
    return text


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description='GeoVac IBM Quantum VQE Demo',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        '--simulator', action='store_true', default=True,
        help='Run on local Aer simulator (default)',
    )
    parser.add_argument(
        '--token', type=str, default=None,
        help='IBM Quantum API token (or set QISKIT_IBM_TOKEN env var)',
    )
    parser.add_argument(
        '--shots', type=int, default=None,
        help='Number of shots for shot-based simulation (default: statevector)',
    )
    parser.add_argument(
        '--compare-gaussian', action='store_true',
        help='Side-by-side comparison with Gaussian STO-3G',
    )
    parser.add_argument(
        '--restarts', type=int, default=5,
        help='Number of COBYLA restarts (default: 5)',
    )
    parser.add_argument(
        '--maxiter', type=int, default=1000,
        help='Max COBYLA iterations per restart (default: 1000)',
    )
    args = parser.parse_args()

    # Resolve token
    token = args.token or os.environ.get('QISKIT_IBM_TOKEN')
    use_hardware = token is not None

    print("=" * 65)
    print("  GeoVac IBM Quantum VQE Demo")
    print("=" * 65)

    # ---------------------------------------------------------------
    # Build GeoVac H2 Hamiltonian
    # ---------------------------------------------------------------
    print("\n  Building H2 GeoVac Hamiltonian...")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        geovac_spo, geovac_info = build_geovac_h2()

    print(f"  System:      H2 (bond-pair encoding)")
    print(f"  Qubits:      {geovac_info['n_qubits']}")
    print(f"  Pauli terms: {geovac_info['n_terms']}")
    print(f"  1-norm:      {geovac_info['one_norm']:.4f} Ha")
    if geovac_info.get('n_qwc_groups'):
        print(f"  QWC groups:  {geovac_info['n_qwc_groups']}")

    # ---------------------------------------------------------------
    # Exact diagonalization (reference)
    # ---------------------------------------------------------------
    print("\n  Computing exact ground state energy...")
    geovac_exact = exact_ground_state(geovac_spo)
    print(f"  Exact energy: {geovac_exact:.6f} Ha")

    # ---------------------------------------------------------------
    # Run VQE
    # ---------------------------------------------------------------
    if use_hardware:
        print("\n  Running VQE on IBM Quantum hardware...")
        geovac_vqe = run_vqe_hardware(
            geovac_spo, token,
            n_restarts=min(args.restarts, 3),
            maxiter=min(args.maxiter, 200),
        )
    elif args.shots is not None:
        print(f"\n  Running VQE with {args.shots:,} shots (simulated noise)...")
        geovac_vqe = run_vqe_shots(
            geovac_spo, shots=args.shots,
            n_restarts=args.restarts,
            maxiter=args.maxiter,
        )
    else:
        print("\n  Running VQE (statevector simulation)...")
        geovac_vqe = run_vqe_statevector(
            geovac_spo,
            n_restarts=args.restarts,
            maxiter=args.maxiter,
        )

    geovac_text = print_results(
        'GeoVac H2 (bond-pair)', geovac_info, geovac_vqe, geovac_exact,
    )

    # ---------------------------------------------------------------
    # Gaussian comparison (if requested)
    # ---------------------------------------------------------------
    comparison_text = ""
    if args.compare_gaussian:
        if use_hardware:
            print("\n  --compare-gaussian only supported in simulator mode.")
        else:
            print("\n  Building Gaussian STO-3G Hamiltonian for comparison...")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                gaussian_spo, gaussian_info = build_gaussian_h2()

            print(f"  Gaussian STO-3G: {gaussian_info['n_qubits']} qubits, "
                  f"{gaussian_info['n_terms']} Pauli terms")

            gaussian_exact = exact_ground_state(gaussian_spo)

            print(f"\n  Running Gaussian VQE (statevector)...")
            if args.shots is not None:
                gaussian_vqe = run_vqe_shots(
                    gaussian_spo, shots=args.shots,
                    n_restarts=args.restarts,
                    maxiter=args.maxiter,
                )
            else:
                gaussian_vqe = run_vqe_statevector(
                    gaussian_spo,
                    n_restarts=args.restarts,
                    maxiter=args.maxiter,
                )

            gaussian_text = print_results(
                'Gaussian STO-3G', gaussian_info, gaussian_vqe, gaussian_exact,
            )
            comparison_text = print_comparison(
                geovac_info, geovac_vqe, geovac_exact,
                gaussian_info, gaussian_vqe, gaussian_exact,
            )

    print("\n" + "=" * 65)
    print("  Demo complete.")
    print("=" * 65)


if __name__ == '__main__':
    main()
