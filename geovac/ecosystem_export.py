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
        # Cached derived quantities (computed lazily)
        self._cached_n_qubits: Optional[int] = None
        self._cached_one_norm: Optional[float] = None

    # ------------------------------------------------------------------
    # Properties (cached after first access)
    # ------------------------------------------------------------------

    @property
    def n_qubits(self) -> int:
        """Number of qubits in the encoding."""
        if self._cached_n_qubits is None:
            max_q = -1
            for term in self._qubit_op.terms:
                for q, _ in term:
                    if q > max_q:
                        max_q = q
            self._cached_n_qubits = max_q + 1 if max_q >= 0 else 0
        return self._cached_n_qubits

    @property
    def n_terms(self) -> int:
        """Number of Pauli terms (including identity)."""
        return len(self._qubit_op.terms)

    @property
    def one_norm(self) -> float:
        """Pauli 1-norm of the electronic Hamiltonian (PK excluded)."""
        if self._cached_one_norm is None:
            self._cached_one_norm = sum(
                abs(c) for c in self._qubit_op.terms.values()
            )
        return self._cached_one_norm

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
        from geovac.pk_partitioning import pk_classical_energy
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

# Maps lowercase key -> canonical name for dispatch
_SYSTEM_REGISTRY: Dict[str, str] = {
    # Atomic / diatomic
    'he': 'He', 'h2': 'H2',
    # First-row main-group hydrides
    'lih': 'LiH', 'beh2': 'BeH2', 'ch4': 'CH4',
    'nh3': 'NH3', 'h2o': 'H2O', 'hf': 'HF',
    # Second-row main-group hydrides ([Ne] frozen core)
    'nah': 'NaH', 'mgh2': 'MgH2', 'sih4': 'SiH4',
    'ph3': 'PH3', 'h2s': 'H2S', 'hcl': 'HCl',
    # Third-row s-block hydrides ([Ar] frozen core)
    'kh': 'KH', 'cah2': 'CaH2',
    # Third-row p-block hydrides ([Ar]3d10 frozen core)
    'geh4': 'GeH4', 'ash3': 'AsH3', 'h2se': 'H2Se', 'hbr': 'HBr',
    # Multi-center diatomics
    'lif': 'LiF', 'co': 'CO', 'n2': 'N2', 'f2': 'F2', 'nacl': 'NaCl',
    # Transition metal hydrides (Z=21-30)
    'sch': 'ScH', 'tih': 'TiH', 'vh': 'VH', 'crh': 'CrH', 'mnh': 'MnH',
    'feh': 'FeH', 'coh': 'CoH', 'nih': 'NiH', 'cuh': 'CuH', 'znh': 'ZnH',
}

# Canonical name -> Z for main-group hydrides (used by _build_hydride)
_HYDRIDE_Z: Dict[str, int] = {
    'LiH': 3, 'BeH2': 4, 'CH4': 6, 'NH3': 7, 'H2O': 8, 'HF': 9,
    'NaH': 11, 'MgH2': 12, 'SiH4': 14, 'PH3': 15, 'H2S': 16, 'HCl': 17,
    'KH': 19, 'CaH2': 20,
    'GeH4': 32, 'AsH3': 33, 'H2Se': 34, 'HBr': 35,
}

# Multi-center diatomics: canonical name → spec factory name
_MULTI_CENTER: Dict[str, str] = {
    'LiF': 'lif_spec', 'CO': 'co_spec', 'N2': 'n2_spec',
    'F2': 'f2_spec', 'NaCl': 'nacl_spec',
}

# Canonical name -> Z for transition metal hydrides
_TM_HYDRIDE_Z: Dict[str, int] = {
    'ScH': 21, 'TiH': 22, 'VH': 23, 'CrH': 24, 'MnH': 25,
    'FeH': 26, 'CoH': 27, 'NiH': 28, 'CuH': 29, 'ZnH': 30,
}


def hamiltonian(
    system: str,
    R: Optional[float] = None,
    max_n: int = 2,
    verbose: bool = False,
    core_method: str = 'pk',
) -> GeoVacHamiltonian:
    """
    Build a GeoVac qubit Hamiltonian and return as a
    ``GeoVacHamiltonian`` object with ecosystem export methods.

    Supports 28 molecular systems: 14 main-group hydrides (first/second/
    third row), 10 transition metal hydrides, plus He and H2.

    Parameters
    ----------
    system : str
        Chemical system name (case-insensitive). Examples:
        'LiH', 'BeH2', 'H2O', 'HF', 'NH3', 'CH4',
        'NaH', 'MgH2', 'SiH4', 'PH3', 'H2S', 'HCl',
        'KH', 'CaH2', 'GeH4', 'AsH3', 'H2Se', 'HBr',
        'ScH', ..., 'ZnH', 'He', 'H2'.
    R : float, optional
        Internuclear distance in bohr. Uses experimental equilibrium
        default if not specified.
    max_n : int
        Maximum principal quantum number (default 2).
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
            f"Supported: {sorted(set(_SYSTEM_REGISTRY.values()))}"
        )

    if canonical in _HYDRIDE_Z:
        return _build_hydride(_HYDRIDE_Z[canonical], R=R, max_n=max_n,
                              verbose=verbose, core_method=core_method)
    elif canonical in _MULTI_CENTER:
        return _build_multi_center(canonical, R=R, max_n=max_n,
                                   verbose=verbose)
    elif canonical in _TM_HYDRIDE_Z:
        return _build_tm_hydride(_TM_HYDRIDE_Z[canonical], R=R,
                                 verbose=verbose)
    elif canonical == 'He':
        return _build_he(max_n=max_n, verbose=verbose)
    elif canonical == 'H2':
        return _build_h2(max_n=max_n, R=R, verbose=verbose)
    else:
        raise ValueError(f"System {canonical!r} not yet implemented.")


# ---------------------------------------------------------------------------
# Unified builder for main-group hydrides (spec-driven)
# ---------------------------------------------------------------------------

def _build_hydride(
    Z: int,
    R: Optional[float] = None,
    max_n: int = 2,
    verbose: bool = False,
    core_method: str = 'pk',
) -> GeoVacHamiltonian:
    """Build any main-group hydride via hydride_spec + general builder."""
    from geovac.molecular_spec import hydride_spec
    from geovac.composed_qubit import build_composed_hamiltonian

    spec = hydride_spec(Z, R=R, max_n=max_n, core_method=core_method)

    # For PK: build without PK in Hamiltonian (partitioned classically)
    # For downfolded: include in Hamiltonian (it IS the effective potential)
    include_in_ham = (core_method == 'downfolded')
    result = build_composed_hamiltonian(
        spec, pk_in_hamiltonian=include_in_ham, verbose=verbose,
    )

    h1_pk = result.get('h1_pk')

    meta = {
        'system': spec.name,
        'R_bohr': R,
        'max_n': max_n,
        'M': result['M'],
        'Q': result['Q'],
        'N_pauli': result['N_pauli'],
    }
    return GeoVacHamiltonian(
        result['qubit_op'], metadata=meta,
        h1_pk=h1_pk,
    )


# ---------------------------------------------------------------------------
# TM hydride builder (spec-driven)
# ---------------------------------------------------------------------------

def _build_multi_center(
    canonical: str,
    R: Optional[float] = None,
    max_n: int = 2,
    verbose: bool = False,
) -> GeoVacHamiltonian:
    """Build a multi-center diatomic via spec factory + general builder."""
    from geovac import molecular_spec as ms
    from geovac.composed_qubit import build_composed_hamiltonian

    factory = getattr(ms, _MULTI_CENTER[canonical])
    kwargs = {'max_n': max_n}
    if R is not None:
        kwargs['R'] = R
    spec = factory(**kwargs)

    result = build_composed_hamiltonian(
        spec, pk_in_hamiltonian=False, verbose=verbose,
    )
    meta = {
        'system': spec.name,
        'R_bohr': R,
        'max_n': max_n,
        'M': result['M'],
        'Q': result['Q'],
        'N_pauli': result['N_pauli'],
    }
    return GeoVacHamiltonian(result['qubit_op'], metadata=meta)


def _build_tm_hydride(
    Z: int,
    R: Optional[float] = None,
    verbose: bool = False,
) -> GeoVacHamiltonian:
    """Build a transition metal hydride (Z=21-30) via spec + general builder."""
    from geovac.molecular_spec import transition_metal_hydride_spec
    from geovac.composed_qubit import build_composed_hamiltonian

    spec = transition_metal_hydride_spec(Z, R=R)
    result = build_composed_hamiltonian(
        spec, pk_in_hamiltonian=False, verbose=verbose,
    )
    meta = {
        'system': spec.name,
        'R_bohr': R,
        'M': result['M'],
        'Q': result['Q'],
        'N_pauli': result['N_pauli'],
    }
    return GeoVacHamiltonian(result['qubit_op'], metadata=meta)


# ---------------------------------------------------------------------------
# He atomic builder
# ---------------------------------------------------------------------------

def _build_he(
    max_n: int = 2,
    verbose: bool = False,
) -> GeoVacHamiltonian:
    """Build He atomic Hamiltonian via JordanWigner encoding."""
    import warnings
    from geovac.lattice_index import LatticeIndex
    from geovac.qubit_encoding import JordanWignerEncoder

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


# ---------------------------------------------------------------------------
# H2 builder (bond-pair encoding)
# ---------------------------------------------------------------------------

def _build_h2(
    max_n: int = 2,
    R: Optional[float] = None,
    verbose: bool = False,
) -> GeoVacHamiltonian:
    """Build H2 qubit Hamiltonian using bond-pair encoding (Z_eff=1)."""
    from geovac.molecular_spec import MolecularSpec, OrbitalBlock
    from geovac.composed_qubit import build_composed_hamiltonian

    if R is None:
        R = 1.4
    # H2 bond-pair: single bond block with Z_eff=1, 2 electrons
    spec = MolecularSpec(
        name='H2',
        blocks=[
            OrbitalBlock(
                label='H2_bond',
                block_type='bond_pair',
                Z_center=1.0,
                n_electrons=2,
                max_n=max_n,
            ),
        ],
        nuclear_repulsion_constant=1.0 / R,
        description=f'H2 bond-pair encoding at R={R:.3f} bohr',
    )
    result = build_composed_hamiltonian(spec, verbose=verbose)
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
