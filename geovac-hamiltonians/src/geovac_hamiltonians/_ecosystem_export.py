"""
Ecosystem Export — Convert GeoVac qubit Hamiltonians to OpenFermion, Qiskit, PennyLane
======================================================================================

Provides a ``GeoVacHamiltonian`` wrapper around the native OpenFermion QubitOperator
produced by the GeoVac pipeline, with lazy export methods for Qiskit and PennyLane.

Each export dependency is optional: the package does not hard-depend on all three
frameworks.  A missing dependency raises ``ImportError`` only when the corresponding
``.to_*()`` method is called.

Convenience entry point::

    from geovac.ecosystem_export import hamiltonian
    H = hamiltonian('LiH', R=3.015, l_max=2)
    print(H.n_terms, H.one_norm)
    qiskit_op  = H.to_qiskit()
    pl_op      = H.to_pennylane()
    of_op      = H.to_openfermion()

Author: GeoVac Development Team
Date: April 2026
"""

from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Optional dependency sentinels
# ---------------------------------------------------------------------------

try:
    from openfermion import QubitOperator as _OFQubitOperator
    _HAS_OPENFERMION = True
except ImportError:
    _HAS_OPENFERMION = False

try:
    from qiskit.quantum_info import SparsePauliOp as _SparsePauliOp
    _HAS_QISKIT = True
except ImportError:
    _HAS_QISKIT = False

try:
    import pennylane as _qml
    _HAS_PENNYLANE = True
except ImportError:
    _HAS_PENNYLANE = False


# ---------------------------------------------------------------------------
# GeoVacHamiltonian wrapper
# ---------------------------------------------------------------------------

class GeoVacHamiltonian:
    """
    Thin wrapper around a GeoVac qubit Hamiltonian (OpenFermion QubitOperator)
    that exposes ecosystem export methods and summary properties.

    For composed systems (LiH, BeH2, H2O), the Phillips-Kleinman
    pseudopotential (PK) is separated from the quantum Hamiltonian by
    default. The PK contribution can be computed classically from the
    1-RDM via ``.pk_classical_energy(one_rdm)``. The ``h1_pk`` attribute
    provides the raw PK one-electron integral matrix.

    Parameters
    ----------
    qubit_op : openfermion.QubitOperator
        The Pauli-string Hamiltonian produced by the GeoVac pipeline
        (Jordan-Wigner encoding). For composed systems, this is the
        electronic-only Hamiltonian (PK excluded).
    metadata : dict, optional
        Additional information (system name, R, l_max, etc.).
    h1_pk : ndarray or None, optional
        PK one-electron integral matrix (M x M, spatial orbitals).
        If provided, enables ``pk_classical_energy()`` and
        ``one_norm_full`` properties.
    qubit_op_full : QubitOperator or None, optional
        The full qubit operator including PK (for one_norm_full).
    """

    def __init__(
        self,
        qubit_op: Any,
        metadata: Optional[Dict[str, Any]] = None,
        h1_pk: Optional[np.ndarray] = None,
        qubit_op_full: Optional[Any] = None,
    ) -> None:
        if not _HAS_OPENFERMION:
            raise ImportError(
                "openfermion is required for GeoVacHamiltonian. "
                "Install with: pip install openfermion"
            )
        self._qubit_op = qubit_op
        self._metadata = metadata or {}
        self._h1_pk = h1_pk
        self._qubit_op_full = qubit_op_full

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def n_qubits(self) -> int:
        """Number of qubits in the encoding."""
        max_q = -1
        for term in self._qubit_op.terms:
            for q, _ in term:
                if q > max_q:
                    max_q = q
        return max_q + 1 if max_q >= 0 else 0

    @property
    def n_terms(self) -> int:
        """Number of Pauli terms (including identity)."""
        return len(self._qubit_op.terms)

    @property
    def one_norm(self) -> float:
        """Pauli 1-norm of the electronic Hamiltonian (PK excluded)."""
        return sum(abs(c) for c in self._qubit_op.terms.values())

    @property
    def one_norm_full(self) -> float:
        """Pauli 1-norm including PK (for reference).

        Returns the same as ``one_norm`` if no PK is present.
        """
        if self._qubit_op_full is not None:
            return sum(abs(c) for c in self._qubit_op_full.terms.values())
        return self.one_norm

    @property
    def h1_pk(self) -> Optional[np.ndarray]:
        """PK one-electron integral matrix (M x M), or None if not applicable."""
        return self._h1_pk

    @property
    def metadata(self) -> Dict[str, Any]:
        """Metadata dict (system, R, l_max, etc.)."""
        return dict(self._metadata)

    # ------------------------------------------------------------------
    # PK classical correction
    # ------------------------------------------------------------------

    def pk_classical_energy(self, one_rdm: np.ndarray) -> float:
        """
        Compute the PK energy contribution from the 1-RDM.

        E_PK = Tr(h1_PK . gamma)

        This is algebraically exact: E_total = E_quantum + E_PK.

        Parameters
        ----------
        one_rdm : ndarray of shape (M, M)
            Spatial 1-RDM from VQE measurement.

        Returns
        -------
        float
            E_PK in Hartree.

        Raises
        ------
        ValueError
            If no PK matrix is available (atomic systems).
        """
        if self._h1_pk is None:
            raise ValueError(
                "No PK matrix available. PK partitioning is only "
                "applicable to composed systems (LiH, BeH2, H2O)."
            )
        from geovac_hamiltonians._pk_partitioning import pk_classical_energy
        return pk_classical_energy(self._h1_pk, one_rdm)

    # ------------------------------------------------------------------
    # Export: OpenFermion
    # ------------------------------------------------------------------

    def to_openfermion(self) -> Any:
        """
        Return the native OpenFermion ``QubitOperator``.

        Returns
        -------
        openfermion.QubitOperator

        Raises
        ------
        ImportError
            If openfermion is not installed.
        """
        if not _HAS_OPENFERMION:
            raise ImportError(
                "openfermion is required. Install with: pip install openfermion"
            )
        return self._qubit_op

    # ------------------------------------------------------------------
    # Export: Qiskit
    # ------------------------------------------------------------------

    def to_qiskit(self) -> Any:
        """
        Convert to a Qiskit ``SparsePauliOp``.

        Uses the reversed-qubit convention that Qiskit expects
        (qubit 0 is the rightmost character in the label string).

        Returns
        -------
        qiskit.quantum_info.SparsePauliOp

        Raises
        ------
        ImportError
            If qiskit is not installed.
        """
        if not _HAS_QISKIT:
            raise ImportError(
                "qiskit is required. Install with: pip install qiskit"
            )
        return _openfermion_to_qiskit(self._qubit_op, self.n_qubits)

    # ------------------------------------------------------------------
    # Export: PennyLane
    # ------------------------------------------------------------------

    def to_pennylane(self) -> Any:
        """
        Convert to a PennyLane ``Hamiltonian``.

        Returns
        -------
        pennylane.Hamiltonian
            A ``qml.Hamiltonian(coeffs, ops)`` constructed from the
            Pauli terms.

        Raises
        ------
        ImportError
            If pennylane is not installed.
        """
        if not _HAS_PENNYLANE:
            raise ImportError(
                "pennylane is required. Install with: pip install pennylane"
            )
        return _openfermion_to_pennylane(self._qubit_op)

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        sys_name = self._metadata.get('system', 'unknown')
        return (
            f"GeoVacHamiltonian(system={sys_name!r}, "
            f"n_qubits={self.n_qubits}, n_terms={self.n_terms}, "
            f"one_norm={self.one_norm:.4f})"
        )


# ---------------------------------------------------------------------------
# Internal converters
# ---------------------------------------------------------------------------

def _openfermion_to_qiskit(
    qubit_op: Any,
    n_qubits: int,
) -> Any:
    """Convert OpenFermion QubitOperator -> Qiskit SparsePauliOp."""
    from qiskit.quantum_info import SparsePauliOp  # type: ignore[import]

    pauli_labels: list[str] = []
    coeffs: list[complex] = []

    for term, coeff in qubit_op.terms.items():
        label = ['I'] * n_qubits
        for qubit_idx, pauli_char in term:
            label[qubit_idx] = pauli_char
        # Qiskit uses reversed qubit ordering (qubit 0 = rightmost)
        pauli_labels.append(''.join(reversed(label)))
        coeffs.append(complex(coeff))

    return SparsePauliOp(pauli_labels, coeffs).simplify()


_OF_PAULI_TO_PL = {
    'X': lambda q: _qml.PauliX(q),
    'Y': lambda q: _qml.PauliY(q),
    'Z': lambda q: _qml.PauliZ(q),
}


def _openfermion_to_pennylane(qubit_op: Any) -> Any:
    """Convert OpenFermion QubitOperator -> PennyLane Hamiltonian."""
    import pennylane as qml  # type: ignore[import]

    coeffs: list[float] = []
    ops: list[Any] = []

    for term, coeff in qubit_op.terms.items():
        coeffs.append(float(np.real(coeff)))
        if len(term) == 0:
            # Identity term
            ops.append(qml.Identity(0))
        else:
            pauli_ops = []
            for qubit_idx, pauli_char in term:
                pauli_ops.append(_OF_PAULI_TO_PL[pauli_char](qubit_idx))
            if len(pauli_ops) == 1:
                ops.append(pauli_ops[0])
            else:
                ops.append(qml.prod(*pauli_ops))

    return qml.Hamiltonian(coeffs, ops)


# ---------------------------------------------------------------------------
# Convenience entry point
# ---------------------------------------------------------------------------

_SYSTEM_REGISTRY: Dict[str, str] = {
    'lih': 'LiH',
    'beh2': 'BeH2',
    'h2o': 'H2O',
    'he': 'He',
    'h2': 'H2',
}


def hamiltonian(
    system: str,
    R: Optional[float] = None,
    l_max: int = 2,
    max_n: int = 2,
    verbose: bool = False,
) -> GeoVacHamiltonian:
    """
    Build a GeoVac composed/atomic qubit Hamiltonian and return as a
    ``GeoVacHamiltonian`` object with export methods.

    Parameters
    ----------
    system : str
        Chemical system: 'LiH', 'BeH2', 'H2O', 'He', or 'H2'.
        Case-insensitive.
    R : float, optional
        Internuclear distance in bohr.  Uses experimental equilibrium
        defaults if not specified (LiH: 3.015, BeH2: 2.502, H2O: 1.809).
    l_max : int
        Maximum angular momentum for composed systems (default 2).
        For atomic systems (He), this controls max_n.
    max_n : int
        Maximum principal quantum number for atomic / H2 systems.
    verbose : bool
        Print build progress.

    Returns
    -------
    GeoVacHamiltonian
        Wrapper with ``.n_terms``, ``.one_norm``, ``.to_openfermion()``,
        ``.to_qiskit()``, ``.to_pennylane()``.
    """
    key = system.strip().lower()
    canonical = _SYSTEM_REGISTRY.get(key)
    if canonical is None:
        raise ValueError(
            f"Unknown system {system!r}. "
            f"Supported: {list(_SYSTEM_REGISTRY.values())}"
        )

    if canonical == 'LiH':
        return _build_lih(R=R, l_max=l_max, verbose=verbose)
    elif canonical == 'BeH2':
        return _build_beh2(R=R, l_max=l_max, verbose=verbose)
    elif canonical == 'H2O':
        return _build_h2o(R=R, l_max=l_max, verbose=verbose)
    elif canonical == 'He':
        return _build_he(max_n=max_n, verbose=verbose)
    elif canonical == 'H2':
        return _build_h2(max_n=max_n, R=R, verbose=verbose)
    else:
        raise ValueError(f"System {canonical!r} not yet implemented.")


# ---------------------------------------------------------------------------
# Per-system builders (thin wrappers around existing pipeline)
# ---------------------------------------------------------------------------

def _build_lih(
    R: Optional[float] = None,
    l_max: int = 2,
    verbose: bool = False,
) -> GeoVacHamiltonian:
    """Build LiH composed Hamiltonian (PK separated by default)."""
    from geovac_hamiltonians._composed_qubit import build_composed_lih

    if R is None:
        R = 3.015
    # Electronic-only Hamiltonian (PK computed but not included)
    result = build_composed_lih(
        max_n_core=l_max, max_n_val=l_max, R=R,
        include_pk=True, pk_in_hamiltonian=False, verbose=verbose,
    )
    # Full Hamiltonian (for one_norm_full reference)
    result_full = build_composed_lih(
        max_n_core=l_max, max_n_val=l_max, R=R,
        include_pk=True, pk_in_hamiltonian=True, verbose=False,
    )
    meta = {
        'system': 'LiH',
        'R_bohr': R,
        'l_max': l_max,
        'M': result.get('M'),
        'Q': result.get('Q'),
        'N_pauli': result.get('N_pauli'),
    }
    return GeoVacHamiltonian(
        result['qubit_op'], metadata=meta,
        h1_pk=result['h1_pk'], qubit_op_full=result_full['qubit_op'],
    )


def _build_beh2(
    R: Optional[float] = None,
    l_max: int = 2,
    verbose: bool = False,
) -> GeoVacHamiltonian:
    """Build BeH2 composed Hamiltonian (PK separated by default)."""
    from geovac_hamiltonians._composed_qubit import build_composed_beh2

    if R is None:
        R = 2.502
    result = build_composed_beh2(
        max_n_core=l_max, max_n_val=l_max, R=R,
        include_pk=True, pk_in_hamiltonian=False, verbose=verbose,
    )
    result_full = build_composed_beh2(
        max_n_core=l_max, max_n_val=l_max, R=R,
        include_pk=True, pk_in_hamiltonian=True, verbose=False,
    )
    meta = {
        'system': 'BeH2',
        'R_bohr': R,
        'l_max': l_max,
        'M': result.get('M'),
        'Q': result.get('Q'),
        'N_pauli': result.get('N_pauli'),
    }
    return GeoVacHamiltonian(
        result['qubit_op'], metadata=meta,
        h1_pk=result['h1_pk'], qubit_op_full=result_full['qubit_op'],
    )


def _build_h2o(
    R: Optional[float] = None,
    l_max: int = 2,
    verbose: bool = False,
) -> GeoVacHamiltonian:
    """Build H2O composed Hamiltonian (PK separated by default)."""
    from geovac_hamiltonians._composed_qubit import build_composed_h2o

    if R is None:
        R = 1.809
    result = build_composed_h2o(
        max_n_core=l_max, max_n_val=l_max, R_OH=R,
        include_pk=True, pk_in_hamiltonian=False, verbose=verbose,
    )
    result_full = build_composed_h2o(
        max_n_core=l_max, max_n_val=l_max, R_OH=R,
        include_pk=True, pk_in_hamiltonian=True, verbose=False,
    )
    meta = {
        'system': 'H2O',
        'R_bohr': R,
        'l_max': l_max,
        'M': result.get('M'),
        'Q': result.get('Q'),
        'N_pauli': result.get('N_pauli'),
    }
    return GeoVacHamiltonian(
        result['qubit_op'], metadata=meta,
        h1_pk=result['h1_pk'], qubit_op_full=result_full['qubit_op'],
    )


def _build_he(
    max_n: int = 2,
    verbose: bool = False,
) -> GeoVacHamiltonian:
    """Build He atomic Hamiltonian via JordanWigner encoding."""
    import warnings
    from geovac_hamiltonians._lattice_index import LatticeIndex
    from geovac_hamiltonians._qubit_encoding import JordanWignerEncoder

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        li = LatticeIndex(
            n_electrons=2, max_n=max_n, nuclear_charge=2,
            vee_method='slater_full', h1_method='hybrid',
        )
    enc = JordanWignerEncoder(li)
    qubit_op = enc.build_qubit_operator()

    meta = {
        'system': 'He',
        'max_n': max_n,
        'Q': li.n_sp,
    }
    return GeoVacHamiltonian(qubit_op, metadata=meta)


def _build_h2(
    max_n: int = 2,
    R: Optional[float] = None,
    basis: Optional[str] = None,
    verbose: bool = False,
) -> GeoVacHamiltonian:
    """
    Build H2 qubit Hamiltonian.

    By default, uses the GeoVac bond-pair encoding (single-center
    hydrogenic basis at Z_eff=1). Pass ``basis='sto-3g'`` to fall back
    to the Gaussian STO-3G reference.

    Parameters
    ----------
    max_n : int
        Maximum principal quantum number for bond-pair encoding.
    R : float or None
        Internuclear distance in bohr (default 1.4).
    basis : str or None
        If 'sto-3g', use Gaussian reference instead of bond-pair.
    verbose : bool
        Print progress.
    """
    if basis is not None and basis.lower().replace('-', '') == 'sto3g':
        # Gaussian STO-3G fallback
        from geovac_hamiltonians._gaussian_reference import h2_sto3g, build_qubit_hamiltonian
        sys_data = h2_sto3g()
        _, qubit_op, _ = build_qubit_hamiltonian(sys_data)
        meta = {
            'system': 'H2',
            'basis': 'STO-3G',
            'R_bohr': 1.4,
            'Q': 4,
        }
        return GeoVacHamiltonian(qubit_op, metadata=meta)

    # GeoVac bond-pair encoding
    from geovac_hamiltonians._composed_qubit import build_h2_bond_pair

    if R is None:
        R = 1.4
    result = build_h2_bond_pair(max_n=max_n, R=R, verbose=verbose)
    meta = {
        'system': 'H2',
        'encoding': 'bond-pair',
        'R_bohr': R,
        'max_n': max_n,
        'M': result['M'],
        'Q': result['Q'],
        'N_pauli': result['N_pauli'],
    }
    return GeoVacHamiltonian(result['qubit_op'], metadata=meta)
