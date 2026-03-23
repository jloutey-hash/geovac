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
from typing import Dict, List, Tuple

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
