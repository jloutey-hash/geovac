"""
Measurement Grouping — Qubitwise-Commuting (QWC) Pauli Partitioning
====================================================================

Partitions a Pauli string Hamiltonian into qubitwise-commuting (QWC)
groups for efficient measurement on near-term quantum hardware.

Two Pauli strings qubitwise-commute if, on every qubit, either they
share the same Pauli operator or at least one has identity.  The number
of QWC groups determines the minimum number of distinct measurement
circuits needed to estimate all expectation values.

Algorithm: greedy sorted insertion (Verteletskyi et al., JCP 152,
124114, 2020).  Terms are sorted by descending |coefficient| and
inserted into the first compatible existing group.  This is O(N^2)
in the number of Pauli terms but fast enough for Q < 200.

Provides:
  - qwc_compatible : pure function checking QWC compatibility
  - qwc_groups     : greedy QWC partitioning
  - count_qwc_groups : convenience wrapper
  - MeasurementAnalysis : dataclass with group statistics
  - analyze_measurement_cost : full analysis
  - compare_measurement_cost : formatted comparison table

Author: GeoVac Development Team
Date: March 2026
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

from openfermion import QubitOperator


# Type alias: a Pauli term is ((qubit, pauli_char), ...) as stored by OpenFermion
PauliTerm = Tuple[Tuple[int, str], ...]


# ---------------------------------------------------------------------------
# QWC compatibility check — pure function, no OpenFermion dependency
# ---------------------------------------------------------------------------

def qwc_compatible(term_a: PauliTerm, term_b: PauliTerm) -> bool:
    """
    Check if two Pauli terms qubitwise-commute (QWC).

    Two terms QWC if, on every qubit where both act non-trivially,
    they have the same Pauli operator.

    Parameters
    ----------
    term_a, term_b : PauliTerm
        Tuples of (qubit_index, pauli_char) pairs, e.g.
        ((0, 'X'), (2, 'Z')).  Identity qubits are omitted
        (OpenFermion convention).

    Returns
    -------
    bool
        True if the terms qubitwise-commute.

    Examples
    --------
    >>> qwc_compatible(((0, 'X'), (2, 'Z')), ((0, 'X'), (3, 'Y')))
    True   # qubit 0: both X; qubits 2,3: one is I
    >>> qwc_compatible(((0, 'X'), (1, 'Z')), ((0, 'Y'), (1, 'Z')))
    False  # qubit 0: X vs Y conflict
    """
    # Build qubit -> pauli map for term_a
    ops_a: Dict[int, str] = {q: p for q, p in term_a}

    # Check every qubit in term_b against term_a
    for q, p in term_b:
        if q in ops_a and ops_a[q] != p:
            return False

    return True


def _qwc_compatible_with_group(
    term: PauliTerm,
    group_ops: Dict[int, str],
) -> bool:
    """
    Check if a term is QWC-compatible with a group's merged operator map.

    Instead of checking against every term in the group, we maintain a
    merged dict mapping qubit -> required Pauli.  A new term is compatible
    iff it doesn't conflict on any qubit.

    Parameters
    ----------
    term : PauliTerm
        Candidate term.
    group_ops : dict
        Merged {qubit: pauli_char} for all terms already in the group.

    Returns
    -------
    bool
    """
    for q, p in term:
        if q in group_ops and group_ops[q] != p:
            return False
    return True


# ---------------------------------------------------------------------------
# Greedy sorted-insertion QWC grouping
# ---------------------------------------------------------------------------

def qwc_groups(
    qubit_op: QubitOperator,
) -> List[List[Tuple[PauliTerm, complex]]]:
    """
    Partition a QubitOperator into qubitwise-commuting (QWC) groups.

    Uses greedy sorted insertion: terms are sorted by descending
    |coefficient|, then each is placed into the first existing group
    where it is QWC-compatible.  If no group fits, a new group is created.

    Parameters
    ----------
    qubit_op : QubitOperator
        Pauli string Hamiltonian from Jordan-Wigner transformation.

    Returns
    -------
    list of list of (PauliTerm, complex)
        Each inner list is a QWC group.  Each element is
        (pauli_tuple, coefficient).

    References
    ----------
    Verteletskyi, Yen & Izmaylov, JCP 152, 124114 (2020).
    """
    # Extract and sort terms by descending |coefficient|
    terms = list(qubit_op.terms.items())
    terms.sort(key=lambda x: abs(x[1]), reverse=True)

    # Groups: list of (term_list, merged_ops_dict)
    groups: List[Tuple[List[Tuple[PauliTerm, complex]], Dict[int, str]]] = []

    for term, coeff in terms:
        placed = False
        for group_terms, group_ops in groups:
            if _qwc_compatible_with_group(term, group_ops):
                # Add term to this group and update merged ops
                group_terms.append((term, coeff))
                for q, p in term:
                    group_ops[q] = p
                placed = True
                break

        if not placed:
            # Create new group
            new_ops: Dict[int, str] = {q: p for q, p in term}
            groups.append(([(term, coeff)], new_ops))

    return [group_terms for group_terms, _ in groups]


def count_qwc_groups(qubit_op: QubitOperator) -> int:
    """Count the number of QWC measurement groups."""
    return len(qwc_groups(qubit_op))


# ---------------------------------------------------------------------------
# Analysis dataclass
# ---------------------------------------------------------------------------

@dataclass
class MeasurementAnalysis:
    """Statistics on QWC measurement grouping."""

    n_qwc_groups: int
    n_pauli_terms: int
    max_group_size: int
    min_group_size: int
    mean_group_size: float
    n_qubits: int

    def summary(self) -> str:
        """One-line summary string."""
        return (
            f"Q={self.n_qubits}, terms={self.n_pauli_terms:,}, "
            f"QWC_groups={self.n_qwc_groups}, "
            f"mean_size={self.mean_group_size:.1f}"
        )


def analyze_measurement_cost(qubit_op: QubitOperator) -> MeasurementAnalysis:
    """
    Full QWC measurement cost analysis.

    Parameters
    ----------
    qubit_op : QubitOperator
        Pauli string Hamiltonian.

    Returns
    -------
    MeasurementAnalysis
    """
    groups = qwc_groups(qubit_op)
    group_sizes = [len(g) for g in groups]

    # Determine n_qubits from the operator
    max_qubit = -1
    for term in qubit_op.terms:
        for q, _ in term:
            if q > max_qubit:
                max_qubit = q
    n_qubits = max_qubit + 1 if max_qubit >= 0 else 0

    return MeasurementAnalysis(
        n_qwc_groups=len(groups),
        n_pauli_terms=len(qubit_op.terms),
        max_group_size=max(group_sizes) if group_sizes else 0,
        min_group_size=min(group_sizes) if group_sizes else 0,
        mean_group_size=sum(group_sizes) / len(groups) if groups else 0.0,
        n_qubits=n_qubits,
    )


# ---------------------------------------------------------------------------
# Comparison table
# ---------------------------------------------------------------------------

def compare_measurement_cost(
    op_a: QubitOperator,
    op_b: QubitOperator,
    labels: Tuple[str, str] = ('GeoVac', 'Gaussian'),
) -> str:
    """
    Formatted comparison table of measurement costs for two Hamiltonians.

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
    a = analyze_measurement_cost(op_a)
    b = analyze_measurement_cost(op_b)

    # Ratios (avoid division by zero)
    term_ratio = a.n_pauli_terms / max(1, b.n_pauli_terms)
    group_ratio = a.n_qwc_groups / max(1, b.n_qwc_groups)

    header = (
        f"{'Metric':<25} {labels[0]:>14} {labels[1]:>14} {'Ratio':>10}"
    )
    sep = f"{'-'*25} {'-'*14} {'-'*14} {'-'*10}"

    rows = [
        f"{'Qubits':<25} {a.n_qubits:>14} {b.n_qubits:>14} {'':>10}",
        f"{'Pauli terms':<25} {a.n_pauli_terms:>14,} {b.n_pauli_terms:>14,} {term_ratio:>10.2f}",
        f"{'QWC groups':<25} {a.n_qwc_groups:>14,} {b.n_qwc_groups:>14,} {group_ratio:>10.2f}",
        f"{'Max group size':<25} {a.max_group_size:>14} {b.max_group_size:>14} {'':>10}",
        f"{'Mean group size':<25} {a.mean_group_size:>14.1f} {b.mean_group_size:>14.1f} {'':>10}",
        f"{'Groups/term ratio':<25} "
        f"{a.n_qwc_groups / max(1, a.n_pauli_terms):>14.3f} "
        f"{b.n_qwc_groups / max(1, b.n_pauli_terms):>14.3f} {'':>10}",
    ]

    lines = ["", "=" * 68, header, sep] + rows + ["=" * 68, ""]
    return "\n".join(lines)
