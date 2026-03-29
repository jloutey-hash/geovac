"""
Trotter Error Bounds — Analytical Cost Metrics for Hamiltonian Simulation
=========================================================================

Computes Trotter error bounds and related cost metrics for Pauli
Hamiltonians.  These metrics quantify the cost of simulating the
Hamiltonian on a fault-tolerant quantum computer.

The key quantity is the Pauli 1-norm:

    lambda = sum_i |c_i|

where c_i are the coefficients in the Pauli decomposition H = sum_i c_i P_i.
This determines the number of LCU queries in qubitization and bounds
the number of Trotter steps required for a given error tolerance.

For first-order Trotterization, the error per step is bounded by:

    error <= (t/r)^2 * lambda^2 / 2

where t is the simulation time and r is the number of Trotter steps.
Inverting for r:

    r >= t * lambda / sqrt(2 * epsilon)

References:
  Childs, Su, Tran, Wiebe & Zhu, PRX 11, 011020 (2021), Eq. 48.
  Campbell, PRL 123, 070503 (2019) — randomized compilation.

Provides:
  - pauli_1norm             : sum of |c_i|
  - max_coefficient         : max |c_i|
  - coefficient_statistics  : full coefficient distribution summary
  - first_order_trotter_bound : error bound for given steps
  - trotter_steps_required  : minimum steps for target error
  - TrotterAnalysis         : dataclass with cost metrics
  - analyze_trotter_cost    : full analysis
  - compare_trotter_cost    : formatted comparison table

Author: GeoVac Development Team
Date: March 2026
"""

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from openfermion import QubitOperator


# ---------------------------------------------------------------------------
# Coefficient metrics
# ---------------------------------------------------------------------------

def pauli_1norm(qubit_op: QubitOperator) -> float:
    """
    Pauli 1-norm: sum of absolute values of all coefficients.

    lambda = sum_i |c_i|  where H = sum_i c_i P_i.

    This is the single most important cost metric for quantum simulation.
    It bounds LCU query counts, Trotter step counts, and randomized
    compilation overhead.

    Parameters
    ----------
    qubit_op : QubitOperator
        Pauli string Hamiltonian.

    Returns
    -------
    float
        The 1-norm (always non-negative).

    References
    ----------
    Childs et al., PRX 11, 011020 (2021).
    """
    return sum(abs(c) for c in qubit_op.terms.values())


def max_coefficient(qubit_op: QubitOperator) -> float:
    """
    Maximum absolute coefficient across all Pauli terms.

    Parameters
    ----------
    qubit_op : QubitOperator

    Returns
    -------
    float
        max |c_i|, or 0.0 if the operator has no terms.
    """
    if not qubit_op.terms:
        return 0.0
    return max(abs(c) for c in qubit_op.terms.values())


def coefficient_statistics(qubit_op: QubitOperator) -> Dict[str, float]:
    """
    Summary statistics of the Pauli coefficient distribution.

    Parameters
    ----------
    qubit_op : QubitOperator

    Returns
    -------
    dict
        Keys: one_norm, max_coeff, median_coeff, mean_coeff,
              n_terms, n_qubits.
    """
    coeffs = np.array([abs(c) for c in qubit_op.terms.values()])

    # Determine n_qubits from operator
    max_qubit = -1
    for term in qubit_op.terms:
        for q, _ in term:
            if q > max_qubit:
                max_qubit = q
    n_qubits = max_qubit + 1 if max_qubit >= 0 else 0

    if len(coeffs) == 0:
        return {
            'one_norm': 0.0,
            'max_coeff': 0.0,
            'median_coeff': 0.0,
            'mean_coeff': 0.0,
            'n_terms': 0,
            'n_qubits': 0,
        }

    return {
        'one_norm': float(np.sum(coeffs)),
        'max_coeff': float(np.max(coeffs)),
        'median_coeff': float(np.median(coeffs)),
        'mean_coeff': float(np.mean(coeffs)),
        'n_terms': len(coeffs),
        'n_qubits': n_qubits,
    }


# ---------------------------------------------------------------------------
# Trotter error bounds
# ---------------------------------------------------------------------------

def first_order_trotter_bound(
    qubit_op: QubitOperator,
    time: float,
    n_steps: int,
) -> float:
    """
    Upper bound on first-order Trotter-Suzuki error.

    Uses the standard lambda^2 bound:

        error <= (t/r)^2 * lambda^2 / 2

    where lambda is the Pauli 1-norm, t is the total simulation time,
    and r is the number of Trotter steps.

    This is a loose but universally applicable bound.  Tighter bounds
    require computing pairwise commutator norms, which is O(N^2) in
    the number of Pauli terms.

    Parameters
    ----------
    qubit_op : QubitOperator
        Pauli string Hamiltonian.
    time : float
        Total simulation time (atomic units).
    n_steps : int
        Number of Trotter steps (r >= 1).

    Returns
    -------
    float
        Upper bound on the operator-norm error.

    References
    ----------
    Childs et al., PRX 11, 011020 (2021), Eq. 48.
    """
    if n_steps < 1:
        raise ValueError(f"n_steps must be >= 1, got {n_steps}")

    lam = pauli_1norm(qubit_op)
    dt = time / n_steps
    return dt**2 * lam**2 / 2.0


def trotter_steps_required(
    qubit_op: QubitOperator,
    time: float,
    epsilon: float,
) -> int:
    """
    Minimum first-order Trotter steps for error <= epsilon.

    Inverts the bound:  r >= t * lambda / sqrt(2 * epsilon).

    Parameters
    ----------
    qubit_op : QubitOperator
        Pauli string Hamiltonian.
    time : float
        Total simulation time (atomic units).
    epsilon : float
        Target error tolerance (> 0).

    Returns
    -------
    int
        Minimum number of Trotter steps (ceiling).
    """
    if epsilon <= 0:
        raise ValueError(f"epsilon must be > 0, got {epsilon}")

    lam = pauli_1norm(qubit_op)
    r = time * lam / math.sqrt(2.0 * epsilon)
    return math.ceil(r)


# ---------------------------------------------------------------------------
# Analysis dataclass
# ---------------------------------------------------------------------------

@dataclass
class TrotterAnalysis:
    """Cost metrics for Trotter simulation of a Pauli Hamiltonian."""

    one_norm: float
    max_coefficient: float
    n_terms: int
    n_qubits: int
    trotter_steps_eps3: int   # steps for epsilon=1e-3 at t=1
    trotter_steps_eps6: int   # steps for epsilon=1e-6 at t=1

    def summary(self) -> str:
        """One-line summary string."""
        return (
            f"Q={self.n_qubits}, terms={self.n_terms:,}, "
            f"lambda={self.one_norm:.4f}, "
            f"r(1e-3)={self.trotter_steps_eps3:,}, "
            f"r(1e-6)={self.trotter_steps_eps6:,}"
        )


def analyze_trotter_cost(qubit_op: QubitOperator) -> TrotterAnalysis:
    """
    Full Trotter cost analysis for a Pauli Hamiltonian.

    Parameters
    ----------
    qubit_op : QubitOperator

    Returns
    -------
    TrotterAnalysis
    """
    stats = coefficient_statistics(qubit_op)

    return TrotterAnalysis(
        one_norm=stats['one_norm'],
        max_coefficient=stats['max_coeff'],
        n_terms=stats['n_terms'],
        n_qubits=stats['n_qubits'],
        trotter_steps_eps3=trotter_steps_required(qubit_op, time=1.0, epsilon=1e-3),
        trotter_steps_eps6=trotter_steps_required(qubit_op, time=1.0, epsilon=1e-6),
    )


# ---------------------------------------------------------------------------
# Comparison table
# ---------------------------------------------------------------------------

def compare_trotter_cost(
    op_a: QubitOperator,
    op_b: QubitOperator,
    labels: Tuple[str, str] = ('GeoVac', 'Gaussian'),
) -> str:
    """
    Formatted comparison table of Trotter costs for two Hamiltonians.

    Parameters
    ----------
    op_a, op_b : QubitOperator
        The two Pauli Hamiltonians to compare.
    labels : tuple of str
        Display names for the two systems.

    Returns
    -------
    str
        Multi-line formatted comparison table.
    """
    a = analyze_trotter_cost(op_a)
    b = analyze_trotter_cost(op_b)

    # Ratios (avoid division by zero)
    norm_ratio = a.one_norm / b.one_norm if b.one_norm > 0 else float('inf')
    steps3_ratio = a.trotter_steps_eps3 / max(1, b.trotter_steps_eps3)
    steps6_ratio = a.trotter_steps_eps6 / max(1, b.trotter_steps_eps6)

    header = (
        f"{'Metric':<25} {labels[0]:>14} {labels[1]:>14} {'Ratio':>10}"
    )
    sep = f"{'-'*25} {'-'*14} {'-'*14} {'-'*10}"

    rows = [
        f"{'Qubits':<25} {a.n_qubits:>14} {b.n_qubits:>14} {'':>10}",
        f"{'Pauli terms':<25} {a.n_terms:>14,} {b.n_terms:>14,} {'':>10}",
        f"{'1-norm (lambda)':<25} {a.one_norm:>14.4f} {b.one_norm:>14.4f} {norm_ratio:>10.2f}",
        f"{'Max |c_i|':<25} {a.max_coefficient:>14.4f} {b.max_coefficient:>14.4f} {'':>10}",
        f"{'Trotter r (eps=1e-3)':<25} {a.trotter_steps_eps3:>14,} {b.trotter_steps_eps3:>14,} {steps3_ratio:>10.2f}",
        f"{'Trotter r (eps=1e-6)':<25} {a.trotter_steps_eps6:>14,} {b.trotter_steps_eps6:>14,} {steps6_ratio:>10.2f}",
        f"{'lambda/Q':<25} {a.one_norm / max(1, a.n_qubits):>14.4f} {b.one_norm / max(1, b.n_qubits):>14.4f} {'':>10}",
    ]

    lines = ["", "=" * 68, header, sep] + rows + ["=" * 68, ""]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Commutator-based Trotter bounds (symplectic Pauli algebra)
# ---------------------------------------------------------------------------

# Binary encoding: I->(0,0), X->(1,0), Y->(1,1), Z->(0,1)
_PAULI_TO_BINARY = {
    'X': (1, 0),
    'Y': (1, 1),
    'Z': (0, 1),
}


def _pauli_string_to_binary(
    pauli_term: tuple,
    n_qubits: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert an OpenFermion Pauli term tuple to binary (x, z) vectors.

    Parameters
    ----------
    pauli_term : tuple
        OpenFermion term key, e.g. ((0, 'X'), (2, 'Z')).
    n_qubits : int
        Total number of qubits (determines vector length).

    Returns
    -------
    x_vec, z_vec : np.ndarray of uint8, shape (n_qubits,)
        Binary symplectic representation.
    """
    x_vec = np.zeros(n_qubits, dtype=np.uint8)
    z_vec = np.zeros(n_qubits, dtype=np.uint8)
    for qubit_idx, pauli_char in pauli_term:
        xb, zb = _PAULI_TO_BINARY[pauli_char]
        x_vec[qubit_idx] = xb
        z_vec[qubit_idx] = zb
    return x_vec, z_vec


def _symplectic_inner_product_matrix(
    X: np.ndarray,
    Z: np.ndarray,
) -> np.ndarray:
    """
    Compute the N×N anticommutativity matrix from binary arrays.

    Two Pauli strings P_j, P_k anticommute iff
        s(j,k) = sum_q [x_j[q]*z_k[q] + z_j[q]*x_k[q]] mod 2 = 1

    Parameters
    ----------
    X : np.ndarray, shape (N, Q), uint8
        Binary x-components of N Pauli strings on Q qubits.
    Z : np.ndarray, shape (N, Q), uint8
        Binary z-components.

    Returns
    -------
    S : np.ndarray, shape (N, N), uint8
        S[j,k] = 1 if P_j and P_k anticommute, 0 if they commute.
    """
    N = X.shape[0]

    # For large N, compute in chunks to avoid memory blowup
    # Each chunk of X @ Z.T is N_chunk × N, stored as uint16 to avoid overflow
    CHUNK = 8192
    if N <= CHUNK:
        # Direct computation
        # Use uint16 to avoid uint8 overflow in matmul (max possible value = Q)
        XZt = X.astype(np.uint16) @ Z.astype(np.uint16).T
        ZXt = Z.astype(np.uint16) @ X.astype(np.uint16).T
        S = ((XZt + ZXt) % 2).astype(np.uint8)
    else:
        S = np.zeros((N, N), dtype=np.uint8)
        for i0 in range(0, N, CHUNK):
            i1 = min(i0 + CHUNK, N)
            Xi = X[i0:i1].astype(np.uint16)
            Zi = Z[i0:i1].astype(np.uint16)
            Z_full = Z.astype(np.uint16)
            X_full = X.astype(np.uint16)
            XZt_chunk = Xi @ Z_full.T        # (chunk, N)
            ZXt_chunk = Zi @ X_full.T        # (chunk, N)
            S[i0:i1, :] = ((XZt_chunk + ZXt_chunk) % 2).astype(np.uint8)

    return S


def pauli_commutator_bound(
    qubit_op: QubitOperator,
    time: float = 1.0,
    epsilon: float = 1e-3,
) -> Dict[str, float]:
    """
    Commutator-based first-order Trotter error bound.

    Uses the symplectic inner product to identify anticommuting pairs,
    then computes the tighter bound:

        ε_comm = (t²/r) Σ_{j<k: anticommuting} |c_j|·|c_k|

    compared to the standard 1-norm bound:

        ε_1norm = (t²/2r) λ²

    Everything is computed in Pauli algebra (binary vectors), with no
    Hilbert-space matrices.

    Parameters
    ----------
    qubit_op : QubitOperator
        Pauli string Hamiltonian.
    time : float
        Simulation time (default 1.0 a.u.).
    epsilon : float
        Target error tolerance for step count computation.

    Returns
    -------
    dict with keys:
        comm_bound : float
            Commutator error bound (for given t, r=1).
        onenorm_bound : float
            Standard λ² bound (for given t, r=1).
        tightening_ratio : float
            comm_bound / onenorm_bound (< 1 means tighter).
        n_anticommuting_pairs : int
        n_total_pairs : int
        anticommuting_fraction : float
        comm_trotter_steps : int
            Steps for given epsilon using commutator bound.
        onenorm_trotter_steps : int
            Steps for given epsilon using 1-norm bound.

    References
    ----------
    Childs, Su, Tran, Wiebe & Zhu, PRX 11, 011020 (2021).
    """
    terms = list(qubit_op.terms.items())
    N = len(terms)

    if N == 0:
        return {
            'comm_bound': 0.0,
            'onenorm_bound': 0.0,
            'tightening_ratio': 0.0,
            'n_anticommuting_pairs': 0,
            'n_total_pairs': 0,
            'anticommuting_fraction': 0.0,
            'comm_trotter_steps': 0,
            'onenorm_trotter_steps': 0,
        }

    # Determine n_qubits
    max_qubit = -1
    for term, _ in terms:
        for q, _ in term:
            if q > max_qubit:
                max_qubit = q
    n_qubits = max_qubit + 1 if max_qubit >= 0 else 0

    # Build binary arrays and coefficient vector
    X = np.zeros((N, n_qubits), dtype=np.uint8)
    Z_arr = np.zeros((N, n_qubits), dtype=np.uint8)
    coeffs = np.zeros(N, dtype=np.float64)

    for i, (term, coeff) in enumerate(terms):
        coeffs[i] = abs(coeff)
        for qubit_idx, pauli_char in term:
            xb, zb = _PAULI_TO_BINARY[pauli_char]
            X[i, qubit_idx] = xb
            Z_arr[i, qubit_idx] = zb

    # Compute anticommutativity matrix
    S = _symplectic_inner_product_matrix(X, Z_arr)

    # Commutator bound: Σ_{j<k: S[j,k]=1} |c_j|·|c_k|
    # Use vectorized: sum of upper triangle of (S * outer(coeffs, coeffs))
    coeff_outer = np.outer(coeffs, coeffs)
    anticomm_weighted = S.astype(np.float64) * coeff_outer

    # Sum upper triangle only (j < k)
    comm_sum = np.sum(np.triu(anticomm_weighted, k=1))

    # Count anticommuting pairs (upper triangle of S)
    n_anticomm = int(np.sum(np.triu(S, k=1)))
    n_total = N * (N - 1) // 2

    # 1-norm bound
    lam = float(np.sum(coeffs))
    onenorm_bound = time**2 * lam**2 / 2.0

    # Commutator bound (factor of 1, not 1/2, because ||[cP, cQ]|| = 2|c||c'|
    # and the sum is over j<k, so comm_error = (t²/r) * Σ_{j<k} 2|c_j||c_k| * (1 if anticomm)
    # But we absorbed the factor: comm_bound = t² * comm_sum
    # (with the factor of 2 from ||[P,Q]|| = 2 already in the bound formula)
    comm_bound = time**2 * comm_sum

    # Tightening ratio
    tightening = comm_bound / onenorm_bound if onenorm_bound > 0 else 0.0

    # Trotter steps: r >= sqrt(comm_bound / epsilon) for comm,
    # r >= t*lambda/sqrt(2*eps) for 1-norm
    if comm_bound > 0 and epsilon > 0:
        comm_steps = math.ceil(math.sqrt(comm_bound / epsilon))
    else:
        comm_steps = 0

    onenorm_steps = math.ceil(time * lam / math.sqrt(2.0 * epsilon)) if epsilon > 0 else 0

    return {
        'comm_bound': comm_bound,
        'onenorm_bound': onenorm_bound,
        'tightening_ratio': tightening,
        'n_anticommuting_pairs': n_anticomm,
        'n_total_pairs': n_total,
        'anticommuting_fraction': n_anticomm / n_total if n_total > 0 else 0.0,
        'comm_trotter_steps': comm_steps,
        'onenorm_trotter_steps': onenorm_steps,
    }


def analyze_commutator_cost(
    qubit_op: QubitOperator,
    epsilons: List[float] = None,
    time: float = 1.0,
) -> Dict[str, Any]:
    """
    Full commutator vs 1-norm Trotter cost comparison.

    Parameters
    ----------
    qubit_op : QubitOperator
        Pauli string Hamiltonian.
    epsilons : list of float
        Error tolerances to compute step counts for.
        Default: [1e-3, 1e-6].
    time : float
        Simulation time (default 1.0 a.u.).

    Returns
    -------
    dict with:
        base : dict from pauli_commutator_bound (at first epsilon)
        step_comparison : list of dicts, one per epsilon, with keys
            epsilon, comm_steps, onenorm_steps, step_ratio
    """
    from typing import Any  # noqa: F811

    if epsilons is None:
        epsilons = [1e-3, 1e-6]

    # Compute the base bound once (the heavy computation)
    base = pauli_commutator_bound(qubit_op, time=time, epsilon=epsilons[0])

    # Recompute step counts for each epsilon
    comm_bound_val = base['comm_bound']
    lam = math.sqrt(2.0 * base['onenorm_bound'] / time**2) if time > 0 else 0.0

    step_comparison = []
    for eps in epsilons:
        if comm_bound_val > 0 and eps > 0:
            c_steps = math.ceil(math.sqrt(comm_bound_val / eps))
        else:
            c_steps = 0

        o_steps = math.ceil(time * lam / math.sqrt(2.0 * eps)) if eps > 0 else 0

        step_comparison.append({
            'epsilon': eps,
            'comm_steps': c_steps,
            'onenorm_steps': o_steps,
            'step_ratio': c_steps / o_steps if o_steps > 0 else 0.0,
        })

    return {
        'base': base,
        'step_comparison': step_comparison,
    }
