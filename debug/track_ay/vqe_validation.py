"""
VQE Validation — Track AY
==========================

Demonstrates that a standard VQE, given a GeoVac/Gaussian Hamiltonian,
converges to the correct ground-state energy.

Step 1: VQE on H2 STO-3G at R=1.4 bohr (4 qubits, 15 Pauli terms)
Step 2: PES scan at multiple R using analytical STO-3G integrals
Step 3: LiH feasibility check (30 qubits — statevector may be infeasible)

Author: GeoVac Development Team / Track AY
Date: April 2026
"""

import csv
import json
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
from scipy.optimize import minimize
from scipy.special import erf

# ---------------------------------------------------------------------------
# STO-3G integral engine for H2 at arbitrary R
# ---------------------------------------------------------------------------
# STO-3G parameters for hydrogen: 3 Gaussians fitted to Slater zeta=1.0
# From Hehre, Stewart, Pople, JCP 51, 2657 (1969)

STO3G_ALPHA = np.array([3.42525091, 0.62391373, 0.16885540])
# Raw contraction coefficients (for normalized primitives)
_STO3G_COEFF_RAW = np.array([0.15432897, 0.53532814, 0.44463454])
# Normalized: multiply by primitive normalization N = (2*alpha/pi)^{3/4}
STO3G_COEFF = _STO3G_COEFF_RAW * (2.0 * STO3G_ALPHA / np.pi) ** 0.75


def _overlap_1d(alpha_a: float, alpha_b: float, Ra: float, Rb: float) -> float:
    """Overlap integral between two 1s Gaussians centered at Ra, Rb."""
    gamma = alpha_a + alpha_b
    Rp = (alpha_a * Ra + alpha_b * Rb) / gamma
    prefactor = np.exp(-alpha_a * alpha_b / gamma * (Ra - Rb) ** 2)
    return prefactor * (np.pi / gamma) ** 1.5


def _kinetic_1d(alpha_a: float, alpha_b: float, Ra: float, Rb: float) -> float:
    """Kinetic energy integral between two 1s Gaussians."""
    gamma = alpha_a + alpha_b
    xi = alpha_a * alpha_b / gamma
    Rab2 = (Ra - Rb) ** 2
    prefactor = np.exp(-xi * Rab2)
    return xi * (3 - 2 * xi * Rab2) * (np.pi / gamma) ** 1.5 * prefactor


def _boys_f0(t: float) -> float:
    """Boys function F0(t) = erf(sqrt(t)) * sqrt(pi) / (2*sqrt(t)) for t>0."""
    if t < 1e-15:
        return 1.0
    return erf(np.sqrt(t)) * np.sqrt(np.pi) / (2.0 * np.sqrt(t))


def _nuclear_attraction(alpha_a: float, alpha_b: float,
                        Ra: float, Rb: float, Rc: float, Zc: float) -> float:
    """Nuclear attraction integral <a|(-Zc/|r-Rc|)|b> for 1s Gaussians in 1D."""
    gamma = alpha_a + alpha_b
    Rp = (alpha_a * Ra + alpha_b * Rb) / gamma
    Rab2 = (Ra - Rb) ** 2
    Rpc2 = (Rp - Rc) ** 2
    prefactor = -2.0 * np.pi / gamma * Zc * np.exp(-alpha_a * alpha_b / gamma * Rab2)
    return prefactor * _boys_f0(gamma * Rpc2)


def _eri_primitive(alpha_a: float, alpha_b: float, alpha_c: float, alpha_d: float,
                   Ra: float, Rb: float, Rc: float, Rd: float) -> float:
    """Two-electron repulsion integral (ab|cd) for 1s Gaussians, chemist notation."""
    gamma_ab = alpha_a + alpha_b
    gamma_cd = alpha_c + alpha_d
    Rp = (alpha_a * Ra + alpha_b * Rb) / gamma_ab
    Rq = (alpha_c * Rc + alpha_d * Rd) / gamma_cd
    Rab2 = (Ra - Rb) ** 2
    Rcd2 = (Rc - Rd) ** 2
    Rpq2 = (Rp - Rq) ** 2
    delta = 1.0 / gamma_ab + 1.0 / gamma_cd

    prefactor = (2.0 * np.pi ** 2.5 /
                 (gamma_ab * gamma_cd * np.sqrt(gamma_ab + gamma_cd)))
    prefactor *= np.exp(-alpha_a * alpha_b / gamma_ab * Rab2)
    prefactor *= np.exp(-alpha_c * alpha_d / gamma_cd * Rcd2)
    return prefactor * _boys_f0(Rpq2 / delta)


def h2_sto3g_integrals(R: float) -> Dict[str, Any]:
    """
    Compute H2 STO-3G integrals at bond length R (bohr).

    Returns dict with h1 (2x2), eri (2x2x2x2) in AO basis,
    overlap S (2x2), nuclear_repulsion, and transformation to MO basis.
    """
    alphas = STO3G_ALPHA
    coeffs = STO3G_COEFF
    n_prim = len(alphas)

    # Atom positions (along z-axis)
    RA = 0.0
    RB = R

    # AO basis: phi_1 centered on A, phi_2 centered on B
    # Each is a contraction of 3 primitives

    # Overlap matrix S (2x2)
    S = np.zeros((2, 2))
    for mu in range(2):
        for nu in range(2):
            R_mu = RA if mu == 0 else RB
            R_nu = RA if nu == 0 else RB
            val = 0.0
            for i in range(n_prim):
                for j in range(n_prim):
                    val += coeffs[i] * coeffs[j] * _overlap_1d(
                        alphas[i], alphas[j], R_mu, R_nu)
            S[mu, nu] = val

    # Kinetic energy T (2x2)
    T = np.zeros((2, 2))
    for mu in range(2):
        for nu in range(2):
            R_mu = RA if mu == 0 else RB
            R_nu = RA if nu == 0 else RB
            val = 0.0
            for i in range(n_prim):
                for j in range(n_prim):
                    val += coeffs[i] * coeffs[j] * _kinetic_1d(
                        alphas[i], alphas[j], R_mu, R_nu)
            T[mu, nu] = val

    # Nuclear attraction V (2x2) — sum over both nuclei
    V = np.zeros((2, 2))
    for mu in range(2):
        for nu in range(2):
            R_mu = RA if mu == 0 else RB
            R_nu = RA if nu == 0 else RB
            val = 0.0
            for i in range(n_prim):
                for j in range(n_prim):
                    # Attraction to nucleus A (Z=1 at RA)
                    val += coeffs[i] * coeffs[j] * _nuclear_attraction(
                        alphas[i], alphas[j], R_mu, R_nu, RA, 1.0)
                    # Attraction to nucleus B (Z=1 at RB)
                    val += coeffs[i] * coeffs[j] * _nuclear_attraction(
                        alphas[i], alphas[j], R_mu, R_nu, RB, 1.0)
            V[mu, nu] = val

    # Core Hamiltonian h_core = T + V (AO basis)
    h_core_ao = T + V

    # Two-electron integrals (chemist notation) in AO basis
    eri_ao = np.zeros((2, 2, 2, 2))
    centers = [RA, RB]
    for mu in range(2):
        for nu in range(2):
            for lam in range(2):
                for sig in range(2):
                    val = 0.0
                    for i in range(n_prim):
                        for j in range(n_prim):
                            for k in range(n_prim):
                                for l in range(n_prim):
                                    val += (coeffs[i] * coeffs[j] *
                                            coeffs[k] * coeffs[l] *
                                            _eri_primitive(
                                                alphas[i], alphas[j],
                                                alphas[k], alphas[l],
                                                centers[mu], centers[nu],
                                                centers[lam], centers[sig]))
                    eri_ao[mu, nu, lam, sig] = val

    # Symmetric orthogonalization: X = S^{-1/2}
    eigvals, eigvecs = np.linalg.eigh(S)
    X = eigvecs @ np.diag(1.0 / np.sqrt(eigvals)) @ eigvecs.T

    # Solve RHF in orthogonal basis
    h_core_orth = X.T @ h_core_ao @ X

    # Initial guess: diagonalize h_core
    eps, C_orth = np.linalg.eigh(h_core_orth)
    C = X @ C_orth

    # SCF iteration
    for scf_iter in range(50):
        # Build density matrix (1 doubly-occupied orbital)
        P = 2.0 * np.outer(C[:, 0], C[:, 0])

        # Build Fock matrix
        J = np.einsum('kl,ijkl->ij', P, eri_ao)
        K = np.einsum('kl,ikjl->ij', P, eri_ao)
        F = h_core_ao + J - 0.5 * K

        # Transform to orthogonal basis and diagonalize
        F_orth = X.T @ F @ X
        eps_new, C_orth_new = np.linalg.eigh(F_orth)
        C_new = X @ C_orth_new

        # Check convergence
        P_new = 2.0 * np.outer(C_new[:, 0], C_new[:, 0])
        if np.max(np.abs(P_new - P)) < 1e-12:
            C = C_new
            eps = eps_new
            break
        C = C_new
        eps = eps_new

    # MO integrals
    h1_mo = C.T @ h_core_ao @ C
    eri_mo = np.einsum('pi,qj,pqrs,rk,sl->ijkl', C, C, eri_ao, C, C)

    nuclear_repulsion = 1.0 / R

    return {
        'h1': h1_mo,
        'eri': eri_mo,
        'nuclear_repulsion': nuclear_repulsion,
        'n_electrons': 2,
        'n_spatial': 2,
        'mo_energies': eps,
        'C': C,
        'S': S,
    }


def build_h2_qubit_hamiltonian(R: float) -> Tuple[Any, Any, int, float]:
    """
    Build H2 STO-3G qubit Hamiltonian at bond length R.

    Returns (SparsePauliOp, QubitOperator, n_qubits, exact_energy).
    """
    from openfermion import jordan_wigner
    from geovac.qubit_encoding import build_fermion_op_from_integrals
    from geovac.vqe_benchmark import openfermion_to_sparse_pauli_op, exact_ground_state_energy

    data = h2_sto3g_integrals(R)
    ferm_op = build_fermion_op_from_integrals(
        data['h1'], data['eri'], data['nuclear_repulsion'],
    )
    of_op = jordan_wigner(ferm_op)
    n_qubits = 2 * data['n_spatial']  # 4
    spo = openfermion_to_sparse_pauli_op(of_op, n_qubits)
    exact_e = exact_ground_state_energy(spo)
    return spo, of_op, n_qubits, exact_e


# ---------------------------------------------------------------------------
# VQE using scipy + statevector (no qiskit-algorithms dependency needed)
# ---------------------------------------------------------------------------

def statevector_expectation(
    params: np.ndarray,
    hamiltonian_matrix: np.ndarray,
    n_qubits: int,
    reps: int = 1,
) -> float:
    """
    Compute <psi(params)|H|psi(params)> via statevector simulation.

    Uses EfficientSU2-like ansatz: layers of Ry-Rz rotations + CX entanglement.
    """
    from qiskit.circuit.library import efficient_su2
    from qiskit.quantum_info import Statevector

    ansatz = efficient_su2(n_qubits, reps=reps, entanglement='linear')
    bound = ansatz.assign_parameters(params)
    sv = Statevector(bound)
    state = sv.data
    energy = np.real(state.conj() @ hamiltonian_matrix @ state)
    return float(energy)


def run_vqe_scipy(
    spo: Any,
    n_qubits: int,
    reps: int = 1,
    maxiter: int = 500,
    n_restarts: int = 3,
) -> Dict[str, Any]:
    """
    Run VQE using scipy.optimize.minimize with COBYLA.

    Parameters
    ----------
    spo : SparsePauliOp
        Qiskit Hamiltonian.
    n_qubits : int
    reps : int
        EfficientSU2 repetitions.
    maxiter : int
        Maximum optimizer iterations per restart.
    n_restarts : int
        Number of random restarts (best result kept).

    Returns
    -------
    dict with keys: energy, n_evals, wall_time, params
    """
    from qiskit.circuit.library import efficient_su2

    # Get dense matrix for fast evaluation
    mat = spo.to_matrix()
    if hasattr(mat, 'toarray'):
        mat = mat.toarray()

    ansatz = efficient_su2(n_qubits, reps=reps, entanglement='linear')
    n_params = ansatz.num_parameters

    best_energy = np.inf
    best_result = None
    total_evals = 0

    t0 = time.time()

    for restart in range(n_restarts):
        x0 = np.random.uniform(-np.pi, np.pi, n_params)

        eval_count = [0]

        def objective(params):
            eval_count[0] += 1
            return statevector_expectation(params, mat, n_qubits, reps)

        result = minimize(
            objective,
            x0,
            method='COBYLA',
            options={'maxiter': maxiter, 'rhobeg': 0.5},
        )

        total_evals += eval_count[0]
        if result.fun < best_energy:
            best_energy = result.fun
            best_result = result

    wall_time = time.time() - t0

    return {
        'energy': best_energy,
        'n_evals': total_evals,
        'wall_time': wall_time,
        'params': best_result.x if best_result is not None else None,
    }


# ---------------------------------------------------------------------------
# VQE using qiskit-algorithms (if available)
# ---------------------------------------------------------------------------

def run_vqe_qiskit(
    spo: Any,
    n_qubits: int,
    reps: int = 1,
    maxiter: int = 500,
) -> Dict[str, Any]:
    """
    Run VQE using qiskit-algorithms StatevectorEstimator + COBYLA.
    """
    from qiskit.circuit.library import efficient_su2
    from qiskit.primitives import StatevectorEstimator
    from qiskit_algorithms import VQE
    from qiskit_algorithms.optimizers import COBYLA

    ansatz = efficient_su2(n_qubits, reps=reps, entanglement='linear')

    eval_count = [0]

    def callback(nfev, parameters, energy, stepsize):
        eval_count[0] = nfev

    estimator = StatevectorEstimator()
    optimizer = COBYLA(maxiter=maxiter)

    vqe = VQE(estimator, ansatz, optimizer, callback=callback)

    t0 = time.time()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        result = vqe.compute_minimum_eigenvalue(spo)
    wall_time = time.time() - t0

    return {
        'energy': float(np.real(result.eigenvalue)),
        'n_evals': eval_count[0],
        'wall_time': wall_time,
    }


# ---------------------------------------------------------------------------
# Main validation
# ---------------------------------------------------------------------------

def validate_h2_single_point() -> Dict[str, Any]:
    """Step 3: VQE on H2 STO-3G at R=1.4 bohr."""
    print("=" * 60)
    print("Step 1: H2 STO-3G VQE at R = 1.4 bohr")
    print("=" * 60)

    # Build Hamiltonian via ecosystem export
    from geovac.ecosystem_export import hamiltonian
    h = hamiltonian('H2')
    print(f"  System: {h.metadata}")
    print(f"  n_qubits: {h.n_qubits}")
    print(f"  n_terms: {h.n_terms}")
    print(f"  one_norm: {h.one_norm:.4f}")

    spo = h.to_qiskit()

    # Exact diagonalization
    from geovac.vqe_benchmark import exact_ground_state_energy
    exact_e = exact_ground_state_energy(spo)
    print(f"  Exact diag energy: {exact_e:.6f} Ha")

    # Run VQE (scipy-based for reliability)
    print("\n  Running VQE (scipy + statevector, 10 restarts)...")
    vqe_result = run_vqe_scipy(spo, h.n_qubits, reps=3, maxiter=1000, n_restarts=10)

    vqe_e = vqe_result['energy']
    error_ha = abs(vqe_e - exact_e)
    error_mha = error_ha * 1000

    print(f"  VQE energy:   {vqe_e:.6f} Ha")
    print(f"  Exact energy: {exact_e:.6f} Ha")
    print(f"  Error: {error_mha:.3f} mHa ({error_ha / abs(exact_e) * 100:.4f}%)")
    print(f"  Evaluations: {vqe_result['n_evals']}")
    print(f"  Wall time: {vqe_result['wall_time']:.1f} s")
    print(f"  Converged within 1 mHa: {'YES' if error_mha < 1.0 else 'NO'}")

    # Also try qiskit-algorithms VQE
    print("\n  Running VQE (qiskit-algorithms)...")
    try:
        qk_result = run_vqe_qiskit(spo, h.n_qubits, reps=3, maxiter=1000)
        qk_e = qk_result['energy']
        qk_error = abs(qk_e - exact_e) * 1000
        print(f"  Qiskit VQE energy: {qk_e:.6f} Ha")
        print(f"  Error: {qk_error:.3f} mHa")
        print(f"  Converged within 1 mHa: {'YES' if qk_error < 1.0 else 'NO'}")
    except Exception as e:
        print(f"  Qiskit-algorithms VQE failed: {e}")
        qk_result = None

    return {
        'exact_energy': exact_e,
        'vqe_energy_scipy': vqe_e,
        'vqe_error_mha_scipy': error_mha,
        'vqe_energy_qiskit': qk_result['energy'] if qk_result else None,
        'n_qubits': h.n_qubits,
        'n_terms': h.n_terms,
    }


def validate_h2_pes_scan() -> List[Dict[str, float]]:
    """Step 4: PES scan over multiple bond distances."""
    print("\n" + "=" * 60)
    print("Step 2: H2 STO-3G PES Scan")
    print("=" * 60)

    from geovac.vqe_benchmark import exact_ground_state_energy

    R_values = [0.5, 0.7, 1.0, 1.2, 1.4, 1.8, 2.2, 2.8, 3.5, 5.0]

    # First validate our integral engine at R=1.4 against published values
    print("\n  Validating integral engine at R=1.4 vs published STO-3G...")
    spo_test, _, _, exact_test = build_h2_qubit_hamiltonian(1.4)
    published_energy = -1.1373  # Szabo & Ostlund
    print(f"  Our engine:  {exact_test:.6f} Ha")
    print(f"  Published:   {published_energy:.4f} Ha")
    print(f"  Difference:  {abs(exact_test - published_energy) * 1000:.3f} mHa")

    print(f"\n  Scanning {len(R_values)} bond distances...")
    print(f"  {'R (bohr)':>10} {'E_exact (Ha)':>14} {'E_VQE (Ha)':>14} "
          f"{'Error (mHa)':>12} {'Converged':>10}")
    print("  " + "-" * 65)

    pes_data = []

    for R in R_values:
        spo, of_op, n_qubits, exact_e = build_h2_qubit_hamiltonian(R)

        # VQE at this geometry
        vqe_result = run_vqe_scipy(spo, n_qubits, reps=3, maxiter=1000, n_restarts=10)
        vqe_e = vqe_result['energy']
        error_mha = abs(vqe_e - exact_e) * 1000
        converged = error_mha < 1.0

        print(f"  {R:10.2f} {exact_e:14.6f} {vqe_e:14.6f} "
              f"{error_mha:12.3f} {'YES' if converged else 'NO':>10}")

        pes_data.append({
            'R_bohr': R,
            'E_exact': exact_e,
            'E_vqe': vqe_e,
            'E_diff_mHa': error_mha,
            'converged': converged,
        })

    return pes_data


def check_lih_feasibility() -> Dict[str, Any]:
    """Step 5: LiH feasibility check."""
    print("\n" + "=" * 60)
    print("Step 3: LiH Feasibility Check")
    print("=" * 60)

    from geovac.ecosystem_export import hamiltonian

    print("  Building LiH composed Hamiltonian (l_max=2)...")
    t0 = time.time()
    h = hamiltonian('LiH', R=3.015, l_max=2)
    build_time = time.time() - t0

    n_qubits = h.n_qubits
    n_terms = h.n_terms
    mem_gb = 2 ** n_qubits * 16 / 1e9  # complex128 statevector

    print(f"  n_qubits: {n_qubits}")
    print(f"  n_terms: {n_terms}")
    print(f"  one_norm: {h.one_norm:.4f}")
    print(f"  Build time: {build_time:.1f} s")
    print(f"  Statevector memory: {mem_gb:.1f} GB")

    feasible = mem_gb < 4.0  # Assume 4 GB limit
    if feasible:
        print(f"  Status: FEASIBLE (< 4 GB)")
        # Try VQE
        print("  Running VQE on LiH...")
        spo = h.to_qiskit()
        from geovac.vqe_benchmark import exact_ground_state_energy
        try:
            exact_e = exact_ground_state_energy(spo)
            print(f"  Exact diag energy: {exact_e:.6f} Ha")
            vqe_result = run_vqe_scipy(spo, n_qubits, reps=1, maxiter=200, n_restarts=1)
            print(f"  VQE energy: {vqe_result['energy']:.6f} Ha")
            print(f"  Error: {abs(vqe_result['energy'] - exact_e) * 1000:.3f} mHa")
        except MemoryError:
            print("  MemoryError during exact diag or VQE — system too large")
            exact_e = None
    else:
        print(f"  Status: INFEASIBLE ({mem_gb:.1f} GB exceeds 4 GB limit)")
        print(f"  Statevector simulation of {n_qubits} qubits requires "
              f"{mem_gb:.1f} GB RAM")
        exact_e = None

    return {
        'n_qubits': n_qubits,
        'n_terms': n_terms,
        'one_norm': h.one_norm,
        'mem_gb': mem_gb,
        'feasible': feasible,
    }


def save_pes_csv(pes_data: List[Dict[str, float]], path: str) -> None:
    """Save PES data to CSV."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['R_bohr', 'E_exact', 'E_vqe', 'E_diff_mHa', 'converged'])
        writer.writeheader()
        writer.writerows(pes_data)
    print(f"\n  PES data saved to {out}")


def save_results_md(h2_result: Dict, pes_data: List[Dict], lih_result: Dict,
                    path: str) -> None:
    """Save summary as markdown."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Track AY: VQE Validation Results",
        "",
        "## H2 STO-3G Single Point (R = 1.4 bohr)",
        "",
        f"| Metric | Value |",
        f"|:-------|------:|",
        f"| Qubits | {h2_result['n_qubits']} |",
        f"| Pauli terms | {h2_result['n_terms']} |",
        f"| Exact diag energy | {h2_result['exact_energy']:.6f} Ha |",
        f"| VQE energy (scipy) | {h2_result['vqe_energy_scipy']:.6f} Ha |",
        f"| VQE error (scipy) | {h2_result['vqe_error_mha_scipy']:.3f} mHa |",
        f"| Converged < 1 mHa | {'YES' if h2_result['vqe_error_mha_scipy'] < 1.0 else 'NO'} |",
    ]

    if h2_result.get('vqe_energy_qiskit') is not None:
        qk_err = abs(h2_result['vqe_energy_qiskit'] - h2_result['exact_energy']) * 1000
        lines.extend([
            f"| VQE energy (qiskit) | {h2_result['vqe_energy_qiskit']:.6f} Ha |",
            f"| VQE error (qiskit) | {qk_err:.3f} mHa |",
        ])

    lines.extend([
        "",
        "## H2 PES Scan (STO-3G, analytical integrals)",
        "",
        "| R (bohr) | E_exact (Ha) | E_VQE (Ha) | Error (mHa) | Converged |",
        "|:---------|:-------------|:-----------|:------------|:----------|",
    ])

    all_converged = True
    for d in pes_data:
        conv = "YES" if d['converged'] else "NO"
        if not d['converged']:
            all_converged = False
        lines.append(
            f"| {d['R_bohr']:.2f} | {d['E_exact']:.6f} | {d['E_vqe']:.6f} "
            f"| {d['E_diff_mHa']:.3f} | {conv} |"
        )

    lines.extend([
        "",
        f"**All points converged within 1 mHa: {'YES' if all_converged else 'NO'}**",
        "",
        "## LiH Feasibility",
        "",
        f"| Metric | Value |",
        f"|:-------|------:|",
        f"| Qubits | {lih_result['n_qubits']} |",
        f"| Pauli terms | {lih_result['n_terms']} |",
        f"| 1-norm | {lih_result['one_norm']:.4f} |",
        f"| Statevector memory | {lih_result['mem_gb']:.1f} GB |",
        f"| Feasible | {'YES' if lih_result['feasible'] else 'NO'} |",
        "",
        "## Method",
        "",
        "- Ansatz: EfficientSU2 (Ry-Rz layers + linear CX entanglement), reps=2",
        "- Optimizer: COBYLA (scipy.optimize.minimize), maxiter=500",
        "- Simulator: Statevector (exact, no noise)",
        "- Multiple random restarts (best of 5 for single point, best of 3 for PES)",
        "- Integral engine: Analytical STO-3G Gaussian integrals for H2 PES scan",
        "- Qubit encoding: Jordan-Wigner via OpenFermion",
        "- Export: GeoVac ecosystem_export -> Qiskit SparsePauliOp",
    ])

    with open(out, 'w') as f:
        f.write('\n'.join(lines))
    print(f"  Results saved to {out}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    np.random.seed(42)

    # Step 1: Single-point VQE
    h2_result = validate_h2_single_point()

    # Step 2: PES scan
    pes_data = validate_h2_pes_scan()

    # Save PES CSV
    save_pes_csv(pes_data, 'debug/track_ay/pes_data.csv')

    # Step 3: LiH feasibility
    lih_result = check_lih_feasibility()

    # Save results
    save_results_md(h2_result, pes_data, lih_result, 'debug/track_ay/results.md')

    print("\n" + "=" * 60)
    print("VQE Validation Complete")
    print("=" * 60)
