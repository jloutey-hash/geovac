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
    print(H.propinquity_bound)  # Paper 38 basis-truncation error estimate
    qiskit_op  = H.to_qiskit()
    pl_op      = H.to_pennylane()
    of_op      = H.to_openfermion()

Propinquity-bound metadata
--------------------------
Every ``GeoVacHamiltonian`` carries a basis-truncation error estimate
``propinquity_bound`` derived from Paper 38's main theorem (Thm.~\\ref{thm:main},
arXiv:2505.xxxxx, Zenodo 2026-05-07).  For the truncated Camporesi--Higuchi
spectral triple at cutoff ``max_n``,

    Lambda(T_{max_n}, T_{S^3}) <= C_3 * gamma_{max_n}

with ``C_3 = 1`` (Lemma L3, sharp at all cutoffs) and ``gamma_{max_n}`` the
central spectral Fejer mass-concentration moment on SU(2) (Lemma L2).  The
asymptotic rate ``gamma_{max_n} ~ (4/pi) log(max_n) / max_n`` carries the
M1 Hopf-base measure signature ``4/pi = Vol(S^2)/pi^2`` (master Mellin engine,
Paper 18 §III.7).  See ``geovac.central_fejer_su2.gamma_n_via_sum_rule`` for
the closed-form sum-rule evaluator.

This is metadata only: the Pauli coefficients of the qubit Hamiltonian are
bit-identical to the un-instrumented build.  Wired 2026-06-07 (v3.85.0).

Author: GeoVac Development Team
Date: April 2026
"""

from __future__ import annotations

import math
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
    h1 : ndarray or None, optional
        One-electron integral matrix (M x M, spatial orbitals, Hartree).
        When provided alongside ``eri`` and ``ecore``, enables the
        ``to_fcidump()`` exporter (Sprint P1, v3.86.0).
    eri : ndarray or None, optional
        Two-electron integral tensor (M, M, M, M) in **chemist notation**
        ``(pq|rs) = eri[p, q, r, s]`` (Hartree).  Matches
        ``build_fermion_op_from_integrals`` and the FCIDUMP standard.
    ecore : float or None, optional
        Frozen-core + nuclear-repulsion constant (Hartree). For composed
        systems this is ``spec.nuclear_repulsion_constant`` which already
        absorbs frozen-core energy and cross-center attraction.
    n_electrons : int or None, optional
        Number of active electrons in the (h1, eri) Hamiltonian
        (frozen-core electrons excluded). Required for the ``NELEC``
        FCIDUMP header field.
    """

    def __init__(
        self,
        qubit_op: Any,
        metadata: Optional[Dict[str, Any]] = None,
        h1_pk: Optional[np.ndarray] = None,
        qubit_op_full: Optional[Any] = None,
        h1: Optional[np.ndarray] = None,
        eri: Optional[np.ndarray] = None,
        ecore: Optional[float] = None,
        n_electrons: Optional[int] = None,
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
        self._h1 = h1
        self._eri = eri
        self._ecore = ecore
        self._n_electrons = n_electrons
        # Cached derived quantities (computed lazily)
        self._cached_n_qubits: Optional[int] = None
        self._cached_one_norm: Optional[float] = None
        self._cached_propinquity_bound: Optional[float] = None

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
        """Metadata dict (system, R, l_max, propinquity bound, etc.).

        Includes ``propinquity_bound`` (Paper 38 Thm.~\\ref{thm:main}
        basis-truncation error estimate) if ``max_n`` is recorded in
        the underlying metadata.
        """
        meta = dict(self._metadata)
        if 'max_n' in meta and 'propinquity_bound' not in meta:
            try:
                meta['propinquity_bound'] = self.propinquity_bound
                meta['propinquity_bound_C3'] = 1.0
                meta['propinquity_bound_asymptotic'] = 4.0 / math.pi
                meta['propinquity_bound_source'] = 'Paper 38 Thm 1; gamma_n via L2 closed-form sum rule'
            except (ValueError, ImportError):
                pass
        return meta

    @property
    def propinquity_bound(self) -> float:
        """Paper 38 basis-truncation error estimate (qualitative-rate).

        Returns the Latrémolière quantum Gromov--Hausdorff propinquity
        upper bound

            Lambda(T_{max_n}, T_{S^3}) <= C_3 * gamma_{max_n}

        with ``C_3 = 1`` (sharp at every cutoff; Lemma L3) and
        ``gamma_{max_n}`` the central spectral Fejer mass-concentration
        moment on SU(2) (Lemma L2), evaluated in closed form via
        ``geovac.central_fejer_su2.gamma_n_via_sum_rule``.

        Requires ``max_n`` to be recorded in the underlying metadata.
        For atomic / molecular ``max_n in {2, 3, 4}`` (production), the
        bound is finite, positive, and strictly monotone-decreasing.

        The bound is *qualitative-rate* — it converges to zero with
        explicit asymptotic constant ``4/pi`` (Paper 38 Thm.~1(ii)) but
        is loose in absolute terms at the production cutoff
        ``max_n = 2`` (gamma_2 ~ 2.07 is the structural-truncation
        bound from the L2 sum rule, not the practically-observed
        accuracy). Use it for cutoff-to-cutoff *improvement claims* and
        for asymptotic statements; do NOT identify it with the
        absolute relative energy error.

        Caches the computed value after first access. The first call
        runs the closed-form sum rule (Theorem 1(i) of the L2
        quantitative-rate memo); O(max_n^2) operations.

        Returns
        -------
        float
            Numerical value of ``gamma_{max_n}``.  ``C_3 = 1`` so this
            *is* the propinquity bound.

        Raises
        ------
        ValueError
            If ``max_n`` is not recorded in metadata, or if ``max_n``
            is below the threshold for the asymptotic statement.
        """
        if self._cached_propinquity_bound is not None:
            return self._cached_propinquity_bound
        max_n = self._metadata.get('max_n')
        if max_n is None:
            raise ValueError(
                "max_n is not recorded in metadata; propinquity_bound "
                "is undefined. Rebuild via geovac.ecosystem_export.hamiltonian "
                "or supply max_n explicitly in metadata."
            )
        if not isinstance(max_n, (int,)) or max_n < 1:
            raise ValueError(
                f"max_n must be a positive integer, got {max_n!r}"
            )
        # Closed-form sum rule from central_fejer_su2 (Paper 38 L2,
        # quantitative-rate memo Thm 1).
        from geovac.central_fejer_su2 import gamma_n_via_sum_rule
        gamma = float(gamma_n_via_sum_rule(max_n, prec=50))
        # C_3 = 1 (Paper 38 Lemma L3, sharp at all cutoffs).
        bound = 1.0 * gamma
        self._cached_propinquity_bound = bound
        return bound

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
    # Pre-JW integrals (h1, eri, ecore) for classical chemistry consumers
    # ------------------------------------------------------------------

    @property
    def h1(self) -> Optional[np.ndarray]:
        """One-electron integral matrix (M, M) in chemist notation, or None."""
        return self._h1

    @property
    def eri(self) -> Optional[np.ndarray]:
        """Two-electron integral tensor (M, M, M, M) in chemist notation (pq|rs), or None."""
        return self._eri

    @property
    def ecore(self) -> Optional[float]:
        """Frozen-core + nuclear-repulsion constant (Hartree), or None."""
        return self._ecore

    @property
    def n_electrons(self) -> Optional[int]:
        """Number of active electrons (frozen-core excluded), or None."""
        return self._n_electrons

    @property
    def n_orbitals(self) -> Optional[int]:
        """Number of spatial orbitals M, derived from ``h1`` shape, or None."""
        if self._h1 is None:
            return None
        return int(self._h1.shape[0])

    # ------------------------------------------------------------------
    # Export: FCIDUMP (classical chemistry interface, Sprint P1)
    # ------------------------------------------------------------------

    def to_fcidump(
        self,
        filename: str,
        *,
        tol: float = 1e-14,
        ms2: int = 0,
        isym: int = 1,
        orbsym: Optional[list] = None,
    ) -> Dict[str, Any]:
        """
        Export the pre-Jordan-Wigner integrals as a Knowles--Handy FCIDUMP file.

        This is the canonical text format for classical quantum-chemistry
        consumers: pyscf, Block2 (DMRG), ipie (AFQMC), Molpro, and others
        all read it directly. The integrals exported here are bit-identical
        to those that ``build_composed_hamiltonian`` / ``build_balanced_hamiltonian``
        feed into ``build_fermion_op_from_integrals`` before Jordan-Wigner
        transformation, so a downstream FCI / CASCI on the FCIDUMP matches
        the GeoVac qubit-Hamiltonian spectrum modulo the qubit-encoding map.

        Conventions
        -----------
        * Chemist notation ``(pq|rs) = eri[p, q, r, s]``.
        * Hartree atomic units throughout.
        * 1-based orbital indices (FCIDUMP standard; converted from
          GeoVac's internal 0-based indexing on write).
        * Eight-fold permutation symmetry assumed for the two-electron
          integrals (real orbitals): only ``(pq|rs)`` with ``p>=q``,
          ``r>=s``, and ``(p,q)>=(r,s)`` lexicographically is written.
        * ``ecore`` (frozen-core + nuclear-repulsion) is written as the
          final all-zeros index line.

        Parameters
        ----------
        filename : str
            Output file path. Standard extension is ``.fcidump``.
        tol : float, optional
            Drop integrals with absolute value below ``tol`` (default
            ``1e-14``). Setting ``tol=0`` writes every integral including
            structural zeros.
        ms2 : int, optional
            ``2 * S_z`` for the reference (default 0, i.e. closed-shell
            singlet).
        isym : int, optional
            ``ISYM`` header field (default 1, totally-symmetric A1).
        orbsym : list of int, optional
            Per-orbital point-group irrep labels. Defaults to all 1s
            (no symmetry, FCIDUMP convention).

        Returns
        -------
        dict
            Metadata about what was written:

            * ``n_orbitals``: M
            * ``n_electrons``: NELEC field
            * ``n_one_body_terms``: number of nonzero h1 entries written
            * ``n_two_body_terms``: number of symmetry-unique eri terms
            * ``ecore``: the constant written
            * ``filename``: path written

        Raises
        ------
        ValueError
            If ``h1``, ``eri``, ``ecore``, or ``n_electrons`` is missing.

        Notes
        -----
        Sprint P1 (v3.86.0, 2026-06-07) unblocks DMRG (Block2), CCSD(T)
        (pyscf), and AFQMC (ipie) consumers at once.  See the hybrid
        pipeline scoping memo for the multi-month roadmap and §1.3 of
        that memo for the convention-pinning checklist.
        """
        if self._h1 is None or self._eri is None:
            raise ValueError(
                "to_fcidump requires the pre-JW (h1, eri) tensors. "
                "Atomic-builder paths (e.g. He via LatticeIndex) currently "
                "do not surface them. Use a composed builder "
                "(hamiltonian('LiH'), hamiltonian('NaH'), etc.) or wire "
                "h1/eri through manually."
            )
        if self._ecore is None:
            raise ValueError(
                "to_fcidump requires ecore (nuclear repulsion + frozen-core "
                "energy). Pass it via the constructor or check builder wiring."
            )
        if self._n_electrons is None:
            raise ValueError(
                "to_fcidump requires n_electrons. Pass it via the "
                "constructor or check builder wiring."
            )

        h1 = np.asarray(self._h1)
        eri = np.asarray(self._eri)
        M = int(h1.shape[0])
        if h1.shape != (M, M):
            raise ValueError(f"h1 must be (M, M); got {h1.shape}")
        if eri.shape != (M, M, M, M):
            raise ValueError(f"eri must be (M, M, M, M); got {eri.shape}")

        if orbsym is None:
            orbsym = [1] * M
        else:
            if len(orbsym) != M:
                raise ValueError(
                    f"orbsym must have length M={M}; got {len(orbsym)}"
                )

        n_one_body = 0
        n_two_body = 0

        with open(filename, 'w', encoding='ascii') as fh:
            # Header: Knowles-Handy &FCI namelist
            fh.write(
                f" &FCI NORB={M:4d},NELEC={self._n_electrons:4d},"
                f"MS2={ms2:4d},\n"
            )
            orbsym_str = ','.join(str(s) for s in orbsym)
            fh.write(f"  ORBSYM={orbsym_str},\n")
            fh.write(f"  ISYM={isym:d},\n")
            fh.write(" &END\n")

            # Two-electron block: 8-fold permutation symmetry for real
            # orbitals. Canonical loop p>=q, r>=s, (p,q)>=(r,s).
            for p in range(M):
                for q in range(p + 1):
                    pq = p * (p + 1) // 2 + q
                    for r in range(M):
                        for s in range(r + 1):
                            rs = r * (r + 1) // 2 + s
                            if rs > pq:
                                continue
                            val = float(eri[p, q, r, s])
                            if abs(val) < tol:
                                continue
                            # FCIDUMP indices are 1-based
                            fh.write(
                                f"{val: .16E}"
                                f"{p+1:5d}{q+1:5d}{r+1:5d}{s+1:5d}\n"
                            )
                            n_two_body += 1

            # One-electron block: Hermitian symmetry h1[p,q] = h1[q,p];
            # write only p>=q.
            for p in range(M):
                for q in range(p + 1):
                    val = float(h1[p, q])
                    if abs(val) < tol:
                        continue
                    fh.write(
                        f"{val: .16E}"
                        f"{p+1:5d}{q+1:5d}{0:5d}{0:5d}\n"
                    )
                    n_one_body += 1

            # Core energy line: p=q=r=s=0
            fh.write(
                f"{float(self._ecore): .16E}"
                f"{0:5d}{0:5d}{0:5d}{0:5d}\n"
            )

        return {
            'n_orbitals': M,
            'n_electrons': int(self._n_electrons),
            'n_one_body_terms': n_one_body,
            'n_two_body_terms': n_two_body,
            'ecore': float(self._ecore),
            'filename': filename,
        }

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
        bound_str = ''
        if 'max_n' in self._metadata:
            try:
                bound_str = f", prop_bound={self.propinquity_bound:.4f}"
            except (ValueError, ImportError):
                pass
        return (
            f"GeoVacHamiltonian(system={sys_name!r}, "
            f"n_qubits={self.n_qubits}, n_terms={self.n_terms}, "
            f"one_norm={self.one_norm:.4f}{bound_str})"
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
# FCIDUMP reader (pure-Python, for round-trip validation without pyscf)
# ---------------------------------------------------------------------------

def read_fcidump(filename: str) -> Dict[str, Any]:
    """
    Parse a Knowles--Handy FCIDUMP file into ``(h1, eri, ecore)`` tensors.

    Used primarily for round-trip validation of ``GeoVacHamiltonian.to_fcidump``
    without requiring pyscf as a hard dependency. Implements only the
    subset of the FCIDUMP grammar that ``to_fcidump`` writes:

    * ``&FCI`` namelist header with ``NORB``, ``NELEC``, ``MS2``,
      ``ORBSYM``, ``ISYM``.
    * Free-form integral records ``val p q r s`` with 1-based indices.
    * Eight-fold permutation symmetry on ``eri`` (real-orbital convention).
    * One-electron entries flagged by ``r = s = 0``.
    * Core energy flagged by ``p = q = r = s = 0``.

    Parameters
    ----------
    filename : str
        Input FCIDUMP file.

    Returns
    -------
    dict
        Keys: ``n_orbitals``, ``n_electrons``, ``ms2``, ``isym``,
        ``orbsym``, ``h1`` (M, M), ``eri`` (M, M, M, M, chemist
        notation, eight-fold symmetrised), ``ecore``.
    """
    import re

    with open(filename, 'r', encoding='ascii') as fh:
        text = fh.read()

    # Header: everything up to &END (case-insensitive)
    header_match = re.search(r'&FCI(.*?)&END', text, re.IGNORECASE | re.DOTALL)
    if header_match is None:
        raise ValueError(f"{filename}: missing &FCI ... &END header")
    header = header_match.group(1)
    body = text[header_match.end():]

    def _grab_scalar(field: str) -> Optional[str]:
        # Scalar fields (NORB, NELEC, MS2, ISYM) end at comma or newline.
        m = re.search(
            rf'{field}\s*=\s*([^,\n]+)', header, re.IGNORECASE,
        )
        return m.group(1).strip() if m else None

    def _grab_list(field: str) -> Optional[str]:
        # List fields (ORBSYM) are comma-separated integers; consume up
        # to the next keyword (NORB/NELEC/MS2/ISYM) or &END.
        m = re.search(
            rf'{field}\s*=\s*([0-9,\s]+?)(?=[A-Za-z]|&|$)',
            header, re.IGNORECASE | re.DOTALL,
        )
        return m.group(1).strip().rstrip(',') if m else None

    norb_str = _grab_scalar('NORB')
    nelec_str = _grab_scalar('NELEC')
    if norb_str is None or nelec_str is None:
        raise ValueError(f"{filename}: NORB / NELEC missing from header")
    M = int(norb_str)
    n_electrons = int(nelec_str)
    ms2 = int(_grab_scalar('MS2') or '0')
    isym = int(_grab_scalar('ISYM') or '1')
    orbsym_str = _grab_list('ORBSYM')
    if orbsym_str:
        orbsym = [int(s) for s in orbsym_str.split(',') if s.strip()]
    else:
        orbsym = [1] * M

    h1 = np.zeros((M, M))
    eri = np.zeros((M, M, M, M))
    ecore = 0.0

    for line in body.splitlines():
        line = line.strip()
        if not line:
            continue
        toks = line.split()
        if len(toks) < 5:
            continue
        val = float(toks[0])
        i = int(toks[1])
        j = int(toks[2])
        k = int(toks[3])
        l = int(toks[4])
        if i == 0 and j == 0 and k == 0 and l == 0:
            ecore = val
            continue
        if k == 0 and l == 0:
            # One-electron: 1-based -> 0-based; expand Hermitian symmetry
            p = i - 1
            q = j - 1
            h1[p, q] = val
            h1[q, p] = val
            continue
        # Two-electron: 1-based -> 0-based; expand 8-fold symmetry on
        # the chemist-notation (pq|rs) tensor.
        p = i - 1
        q = j - 1
        r = k - 1
        s = l - 1
        for (a, b, c, d) in (
            (p, q, r, s), (q, p, r, s),
            (p, q, s, r), (q, p, s, r),
            (r, s, p, q), (s, r, p, q),
            (r, s, q, p), (s, r, q, p),
        ):
            eri[a, b, c, d] = val

    return {
        'n_orbitals': M,
        'n_electrons': n_electrons,
        'ms2': ms2,
        'isym': isym,
        'orbsym': orbsym,
        'h1': h1,
        'eri': eri,
        'ecore': ecore,
    }


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
    # Fifth/sixth-row s-block alkaline-earth monohydrides
    # (Sprint 3 HA-C, v2.12.0, unblocks Sunaga 2025 comparison)
    'srh': 'SrH', 'bah': 'BaH',
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
    tapered: Optional[str] = None,
) -> GeoVacHamiltonian:
    """
    Build a GeoVac qubit Hamiltonian and return as a
    ``GeoVacHamiltonian`` object with ecosystem export methods.

    Supports 37 molecular systems (the ``_SYSTEM_REGISTRY``): 18 main-group
    hydrides (first through fourth row), 5 multi-center molecules (LiF, CO,
    N2, F2, NaCl), 10 transition metal hydrides (ScH-ZnH), 2 alkaline-earth
    monohydrides (SrH, BaH), plus He and H2.

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
    tapered : {None, 'global', 'per_block', 'extended', 'full'}, optional
        If set, apply Z2 tapering on top of the standard alpha/beta
        parity Z-strings (Paper 14 §sec:hopf_tapering, Paper 29 §5.3,
        plus v3.89.0 / v3.92.0 symmetry extensions):

        - ``'global'``: Hopf-U(1) global ``m -> -m`` reflection; ΔQ=3.
        - ``'per_block'``: Hopf per-sub-block; ΔQ=``2 + n_sub_blocks``
          (recommended baseline; range 3 to 12).
        - ``'extended'``: Hopf per-block + Gaunt-parity ℓ-Z2 per-block
          (v3.89.0; strict improvement over ``'per_block'`` on the
          composed builder: ΔQ = ``2 + 2·n_sub_blocks_with_lodd`` plus
          3% Pauli reduction).  Atom-swap and inversion stabilizers
          are NOT applied (opt-in only via
          ``geovac.extended_tapering.extended_tapered_from_spec``
          due to their Pauli-inflation tradeoff).
        - ``'full'``: ``'extended'`` plus per-sub-block
          particle-conservation Z₂s ``(-1)^N_b`` found by the
          symmetry-adapted basis meta-investigation (v3.92.0).
          Saves an additional 1–8 qubits per molecule on top of
          ``'extended'`` (LiH +2, HF +5, BeH₂ +1, H₂O +1, NH₃ +7,
          CH₄ +8 in the verification panel).  Strict improvement
          over ``'extended'`` in qubit count.  Backed by
          ``geovac.symmetry_adapted_basis.extended_plus_hidden_tapered_from_spec``.

        Default ``None`` preserves the historical Pauli/qubit counts
        of Paper 14 Tables I/II.  The tapering is a similarity
        transform onto the fully-symmetric sector.

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
        ham = _build_hydride(_HYDRIDE_Z[canonical], R=R, max_n=max_n,
                             verbose=verbose, core_method=core_method)
    elif canonical in _MULTI_CENTER:
        ham = _build_multi_center(canonical, R=R, max_n=max_n,
                                  verbose=verbose)
    elif canonical in _TM_HYDRIDE_Z:
        ham = _build_tm_hydride(_TM_HYDRIDE_Z[canonical], R=R,
                                verbose=verbose)
    elif canonical in ('SrH', 'BaH'):
        ham = _build_alkaline_earth_monohydride(
            canonical, R=R, max_n=max_n, verbose=verbose,
        )
    elif canonical == 'He':
        ham = _build_he(max_n=max_n, verbose=verbose)
    elif canonical == 'H2':
        ham = _build_h2(max_n=max_n, R=R, verbose=verbose)
    else:
        raise ValueError(f"System {canonical!r} not yet implemented.")

    if tapered is not None:
        ham = _apply_hopf_tapering_to_geovac_hamiltonian(
            ham, canonical, R=R, max_n=max_n, mode=tapered,
            core_method=core_method, verbose=verbose,
        )
    return ham


def _apply_hopf_tapering_to_geovac_hamiltonian(
    ham: GeoVacHamiltonian,
    canonical: str,
    R: Optional[float],
    max_n: int,
    mode: str,
    core_method: str,
    verbose: bool,
) -> GeoVacHamiltonian:
    """Rebuild the molecule and apply the requested tapering, wrapping
    the result in a ``GeoVacHamiltonian``.

    Supported modes:
      - ``'global'``, ``'per_block'``: via
        ``geovac.z2_tapering.hopf_tapered_from_spec``
      - ``'extended'``: via
        ``geovac.extended_tapering.extended_tapered_from_spec``
        (Hopf + ℓ-parity per-block; atom-swap/inversion off here)

    Re-runs the spec build (cheap relative to the JW + rotation cost)
    rather than trying to taper an already-JW-encoded operator that was
    not built in the P-eigenbasis.
    """
    if mode not in ('global', 'per_block', 'extended', 'full'):
        raise ValueError(
            f"tapered must be None, 'global', 'per_block', "
            f"'extended', or 'full'; got {mode!r}"
        )

    spec = _rebuild_spec(canonical, R=R, max_n=max_n, core_method=core_method)
    pk_in_ham = (core_method == 'downfolded')

    if mode in ('global', 'per_block'):
        from geovac.z2_tapering import hopf_tapered_from_spec
        result = hopf_tapered_from_spec(
            spec, mode=mode, pk_in_hamiltonian=pk_in_ham, verbose=verbose,
        )
        meta_update = {
            'tapered_mode': mode,
            'Q_tapered': result['Q_tapered'],
            'delta_Q': result['delta_Q'],
            'n_sub_blocks': result['n_sub_blocks'],
            'n_sub_blocks_with_antisym': result['n_sub_blocks_with_antisym'],
            'dropped_P_indices': result['dropped_P_indices'],
        }
        qop_tapered = result['qubit_op_tapered']
    elif mode == 'extended':
        from geovac.extended_tapering import extended_tapered_from_spec
        result = extended_tapered_from_spec(
            spec,
            use_hopf=True,
            use_ell_parity=True,
            use_atom_swap=False,
            use_inversion=False,
            pk_in_hamiltonian=pk_in_ham,
            builder='composed',
            verbose=verbose,
        )
        meta_update = {
            'tapered_mode': mode,
            'Q_tapered': result['Q_tapered'],
            'delta_Q': result['delta_Q'],
            'n_stabs_kept': result['n_stabs_kept'],
            'kinds_kept': result['kinds_kept'],
            'dropped': result['dropped'],
        }
        qop_tapered = result['qubit_op_tapered']
    else:  # 'full' — extended + hidden particle-conservation Z2s (v3.92.0)
        from geovac.symmetry_adapted_basis import (
            extended_plus_hidden_tapered_from_spec,
        )
        result = extended_plus_hidden_tapered_from_spec(
            spec,
            use_hopf=True,
            use_ell_parity=True,
            use_atom_swap=False,
            use_inversion=False,
            pk_in_hamiltonian=pk_in_ham,
            builder='composed',
            verbose=verbose,
        )
        meta_update = {
            'tapered_mode': mode,
            'Q_tapered': result['Q_tapered_hidden'],
            'Q_tapered_extended_only': result['Q_tapered'],
            'delta_Q': result['Q_naive'] - result['Q_tapered_hidden'],
            'delta_Q_hidden': result['delta_Q_hidden'],
            'n_hidden_z2_kept': len(result.get('hidden_z2_kept', [])),
            'n_stabs_kept': result['n_stabs_kept'],
            'kinds_kept': result['kinds_kept'],
        }
        qop_tapered = result['qubit_op_tapered_plus_hidden']

    meta = dict(ham._metadata)
    meta.update(meta_update)
    return GeoVacHamiltonian(
        qop_tapered, metadata=meta,
        h1_pk=ham._h1_pk,
        # The pre-JW integrals are pre-tapering; they remain valid for
        # classical chemistry consumers regardless of qubit tapering.
        h1=ham._h1,
        eri=ham._eri,
        ecore=ham._ecore,
        n_electrons=ham._n_electrons,
    )


def _rebuild_spec(
    canonical: str, R: Optional[float], max_n: int, core_method: str,
) -> Any:
    """Rebuild the MolecularSpec for a given canonical system name.
    Mirrors the dispatch logic of :func:`hamiltonian` at the spec level.
    """
    from geovac import molecular_spec as ms

    if canonical in _HYDRIDE_Z:
        return ms.hydride_spec(_HYDRIDE_Z[canonical], R=R, max_n=max_n,
                               core_method=core_method)
    if canonical in _MULTI_CENTER:
        factory = getattr(ms, _MULTI_CENTER[canonical])
        kwargs: Dict[str, Any] = {'max_n': max_n}
        if R is not None:
            kwargs['R'] = R
        return factory(**kwargs)
    if canonical in _TM_HYDRIDE_Z:
        return ms.transition_metal_hydride_spec(_TM_HYDRIDE_Z[canonical], R=R)
    if canonical == 'SrH':
        return ms.srh_spec(R=R, max_n=max_n)
    if canonical == 'BaH':
        return ms.bah_spec(R=R, max_n=max_n)
    if canonical == 'H2':
        from geovac.molecular_spec import MolecularSpec, OrbitalBlock
        Rh = R if R is not None else 1.4
        return MolecularSpec(
            name='H2',
            blocks=[OrbitalBlock(
                label='H2_bond', block_type='bond_pair', Z_center=1.0,
                n_electrons=2, max_n=max_n,
            )],
            nuclear_repulsion_constant=1.0 / Rh,
        )
    if canonical == 'He':
        from geovac.molecular_spec import MolecularSpec, OrbitalBlock
        return MolecularSpec(
            name='He',
            blocks=[OrbitalBlock(
                label='He_core', block_type='atomic', Z_center=2.0,
                n_electrons=2, max_n=max_n,
            )],
            nuclear_repulsion_constant=0.0,
        )
    raise ValueError(
        f"_rebuild_spec: no spec factory for {canonical!r}"
    )


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
        h1=result.get('h1'),
        eri=result.get('eri'),
        ecore=result.get('nuclear_repulsion'),
        n_electrons=sum(b.n_electrons for b in spec.blocks),
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
    return GeoVacHamiltonian(
        result['qubit_op'], metadata=meta,
        h1=result.get('h1'),
        eri=result.get('eri'),
        ecore=result.get('nuclear_repulsion'),
        n_electrons=sum(b.n_electrons for b in spec.blocks),
    )


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
    # transition_metal_hydride_spec uses the standard composed cutoff;
    # max_n=2 is the spec-factory default. Record it so the propinquity
    # bound is well-defined.
    meta = {
        'system': spec.name,
        'R_bohr': R,
        'max_n': 2,
        'M': result['M'],
        'Q': result['Q'],
        'N_pauli': result['N_pauli'],
    }
    return GeoVacHamiltonian(
        result['qubit_op'], metadata=meta,
        h1=result.get('h1'),
        eri=result.get('eri'),
        ecore=result.get('nuclear_repulsion'),
        n_electrons=sum(b.n_electrons for b in spec.blocks),
    )


def _build_alkaline_earth_monohydride(
    canonical: str,
    R: Optional[float] = None,
    max_n: int = 2,
    verbose: bool = False,
) -> GeoVacHamiltonian:
    """Build a heavy-atom alkaline-earth monohydride (SrH, BaH).

    Sprint 3 HA-C (v2.12.0): these monohydrides use the dihydride-minus-one-
    bond block template with a [Kr] (SrH) or [Xe] (BaH) frozen core.
    Unblocks the Sunaga 2025 (PRA 111, 022817) matched-Q comparison.
    """
    from geovac.molecular_spec import srh_spec, bah_spec
    from geovac.composed_qubit import build_composed_hamiltonian

    if canonical == 'SrH':
        spec = srh_spec(R=R, max_n=max_n)
    elif canonical == 'BaH':
        spec = bah_spec(R=R, max_n=max_n)
    else:
        raise ValueError(f"Unexpected monohydride {canonical!r}")

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
    return GeoVacHamiltonian(
        result['qubit_op'], metadata=meta,
        h1=result.get('h1'),
        eri=result.get('eri'),
        ecore=result.get('nuclear_repulsion'),
        n_electrons=sum(b.n_electrons for b in spec.blocks),
    )


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
    return GeoVacHamiltonian(
        result['qubit_op'], metadata=meta,
        h1=result.get('h1'),
        eri=result.get('eri'),
        ecore=result.get('nuclear_repulsion'),
        n_electrons=sum(b.n_electrons for b in spec.blocks),
    )
