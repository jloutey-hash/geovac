"""
Qubit Encoding — Jordan-Wigner transformation for GeoVac Hamiltonians
======================================================================

Converts sparse second-quantized Hamiltonians from LatticeIndex or
MolecularLatticeIndex into qubit (Pauli string) representations via
the Jordan-Wigner transformation.

Provides:
  - JordanWignerEncoder: main encoder class
  - PauliAnalysis: dataclass with term counts, sparsity, max weight

The key result (Paper 13, Sec XII.2): GeoVac Hamiltonians scale as
O(Q^3.15) in Pauli terms vs O(Q^4.60) for Gaussian bases, due to
angular momentum selection rules enforcing ERI sparsity ~ 1/M^2.

Requires: openfermion (pip install openfermion)

Author: GeoVac Development Team
Date: March 2026
"""

from dataclasses import dataclass
from typing import Dict, Optional, Tuple, Union

import numpy as np

try:
    from openfermion import FermionOperator, QubitOperator, jordan_wigner, count_qubits
except ImportError:
    raise ImportError(
        "openfermion is required for qubit encoding. "
        "Install with: pip install openfermion"
    )


@dataclass
class PauliAnalysis:
    """Analysis of a qubit Hamiltonian's Pauli string decomposition."""

    n_qubits: int
    n_pauli_terms: int
    max_pauli_weight: int
    sparsity: float  # fraction of possible 4^Q terms that are nonzero
    n_h1_nonzero: int
    n_eri_nonzero: int
    n_eri_possible: int
    eri_density: float

    @property
    def pauli_density(self) -> float:
        """Fraction of possible Pauli strings that appear."""
        return self.sparsity

    def summary(self) -> str:
        """One-line summary string."""
        return (
            f"Q={self.n_qubits}, Pauli={self.n_pauli_terms:,}, "
            f"max_weight={self.max_pauli_weight}, "
            f"ERI_density={self.eri_density:.1%}"
        )


class JordanWignerEncoder:
    """
    Encodes a GeoVac Hamiltonian as a qubit Hamiltonian via Jordan-Wigner.

    Parameters
    ----------
    lattice_index : LatticeIndex or MolecularLatticeIndex
        A constructed GeoVac lattice index with one-body (H1) and
        two-body (ERI) integrals already computed.
    nuclear_repulsion : float
        Nuclear repulsion energy (V_NN). For atoms this is 0.
        For MolecularLatticeIndex, extracted automatically.

    Examples
    --------
    >>> from geovac import LatticeIndex
    >>> idx = LatticeIndex(n_electrons=2, max_n=3, nuclear_charge=2,
    ...                    vee_method='slater_full', h1_method='hybrid')
    >>> enc = JordanWignerEncoder(idx)
    >>> analysis = enc.analyze()
    >>> print(analysis.n_pauli_terms)
    """

    def __init__(
        self,
        lattice_index: object,
        nuclear_repulsion: Optional[float] = None,
    ) -> None:
        self._li = lattice_index
        self._fermion_op: Optional[FermionOperator] = None
        self._qubit_op: Optional[QubitOperator] = None
        self._analysis: Optional[PauliAnalysis] = None

        # Extract nuclear repulsion
        if nuclear_repulsion is not None:
            self._v_nn = nuclear_repulsion
        elif hasattr(lattice_index, 'V_NN'):
            self._v_nn = lattice_index.V_NN
        else:
            self._v_nn = 0.0

        # Extract dimensions
        self._n_spinorb = lattice_index.n_sp
        self._n_spatial = self._n_spinorb // 2

    def build_fermion_operator(self) -> FermionOperator:
        """
        Build second-quantized FermionOperator from lattice integrals.

        Uses the one-body matrix H1 and two-electron integral dictionary
        from the lattice index. GeoVac ERIs are in physicist notation
        <ab|cd> = integral phi_a(1) phi_b(2) (1/r12) phi_c(1) phi_d(2).

        Returns
        -------
        FermionOperator
            The full electronic Hamiltonian in second quantization.
        """
        if self._fermion_op is not None:
            return self._fermion_op

        li = self._li
        n_spatial = self._n_spatial

        # Extract H1 as dense array
        H1_dense = np.asarray(li._H1_spatial.todense())

        # Extract ERI dictionary
        eri_dict = li._eri

        fermion_op = FermionOperator((), self._v_nn)

        # One-body: H1[p,q] -> sum_sigma a+_{2p+sigma} a_{2q+sigma}
        for p in range(n_spatial):
            for q in range(n_spatial):
                h1_val = H1_dense[p, q]
                if abs(h1_val) < 1e-12:
                    continue
                for sigma in range(2):
                    sp_p = 2 * p + sigma
                    sp_q = 2 * q + sigma
                    fermion_op += FermionOperator(
                        ((sp_p, 1), (sp_q, 0)), h1_val
                    )

        # Two-body: GeoVac ERI in physicist notation <ab|cd>
        # OpenFermion: 0.5 * <ab|cd> * sum_{sigma,tau} a+_a,s a+_b,t a_d,t a_c,s
        for (a, b, c, d), val in eri_dict.items():
            if abs(val) < 1e-12:
                continue
            coeff = 0.5 * val
            for sigma in range(2):
                for tau in range(2):
                    sp_a = 2 * a + sigma
                    sp_b = 2 * b + tau
                    sp_c = 2 * c + sigma
                    sp_d = 2 * d + tau
                    if sp_a == sp_b:
                        continue  # Pauli exclusion
                    fermion_op += FermionOperator(
                        ((sp_a, 1), (sp_b, 1), (sp_d, 0), (sp_c, 0)),
                        coeff,
                    )

        self._fermion_op = fermion_op
        return fermion_op

    def build_qubit_operator(self) -> QubitOperator:
        """
        Apply Jordan-Wigner transformation to get qubit Hamiltonian.

        Returns
        -------
        QubitOperator
            Pauli string decomposition of the Hamiltonian.
        """
        if self._qubit_op is not None:
            return self._qubit_op

        fermion_op = self.build_fermion_operator()
        self._qubit_op = jordan_wigner(fermion_op)
        return self._qubit_op

    def count_pauli_terms(self) -> int:
        """Count distinct Pauli string terms in the qubit Hamiltonian."""
        qubit_op = self.build_qubit_operator()
        return len(qubit_op.terms)

    def max_pauli_weight(self) -> int:
        """
        Maximum Pauli weight (number of non-identity operators in any term).

        For Jordan-Wigner, this is bounded by 2*Q for two-body terms
        but typically much smaller due to sparsity.
        """
        qubit_op = self.build_qubit_operator()
        max_w = 0
        for term in qubit_op.terms:
            if len(term) > max_w:
                max_w = len(term)
        return max_w

    def analyze(self) -> PauliAnalysis:
        """
        Full analysis of the Pauli decomposition.

        Returns
        -------
        PauliAnalysis
            Dataclass with term counts, sparsity metrics, and ERI statistics.
        """
        if self._analysis is not None:
            return self._analysis

        qubit_op = self.build_qubit_operator()
        n_qubits = self._n_spinorb
        n_pauli = len(qubit_op.terms)

        # Max Pauli weight
        max_w = 0
        for term in qubit_op.terms:
            if len(term) > max_w:
                max_w = len(term)

        # ERI statistics
        n_eri_nz = len(self._li._eri)
        n_eri_possible = self._n_spatial ** 4

        # H1 nonzero count
        H1_dense = np.asarray(self._li._H1_spatial.todense())
        n_h1_nz = int(np.count_nonzero(np.abs(H1_dense) > 1e-12))

        # Sparsity: fraction of 4^Q possible Pauli strings
        # (for large Q this is astronomically small)
        if n_qubits <= 20:
            sparsity = n_pauli / (4 ** n_qubits)
        else:
            sparsity = 0.0  # too large to represent

        self._analysis = PauliAnalysis(
            n_qubits=n_qubits,
            n_pauli_terms=n_pauli,
            max_pauli_weight=max_w,
            sparsity=sparsity,
            n_h1_nonzero=n_h1_nz,
            n_eri_nonzero=n_eri_nz,
            n_eri_possible=n_eri_possible,
            eri_density=n_eri_nz / max(1, n_eri_possible),
        )
        return self._analysis


def build_fermion_op_from_integrals(
    h1: np.ndarray,
    eri: np.ndarray,
    nuclear_repulsion: float = 0.0,
) -> FermionOperator:
    """
    Build a FermionOperator from dense one- and two-electron integrals.

    Useful for constructing Gaussian-basis reference Hamiltonians.

    Parameters
    ----------
    h1 : np.ndarray, shape (M, M)
        One-electron integrals in spatial orbital basis.
    eri : np.ndarray, shape (M, M, M, M)
        Two-electron integrals in CHEMIST notation: (pq|rs).
    nuclear_repulsion : float
        Nuclear repulsion energy constant.

    Returns
    -------
    FermionOperator
        Second-quantized Hamiltonian.
    """
    n_spatial = h1.shape[0]
    fermion_op = FermionOperator((), nuclear_repulsion)

    # One-body
    for p in range(n_spatial):
        for q in range(n_spatial):
            if abs(h1[p, q]) < 1e-12:
                continue
            for sigma in range(2):
                sp_p = 2 * p + sigma
                sp_q = 2 * q + sigma
                fermion_op += FermionOperator(
                    ((sp_p, 1), (sp_q, 0)), h1[p, q]
                )

    # Two-body: chemist notation (pq|rs)
    for p in range(n_spatial):
        for q in range(n_spatial):
            for r in range(n_spatial):
                for s in range(n_spatial):
                    if abs(eri[p, q, r, s]) < 1e-12:
                        continue
                    coeff = 0.5 * eri[p, q, r, s]
                    for sigma in range(2):
                        for tau in range(2):
                            sp_p = 2 * p + sigma
                            sp_q = 2 * q + sigma
                            sp_r = 2 * r + tau
                            sp_s = 2 * s + tau
                            if sp_p == sp_r:
                                continue
                            fermion_op += FermionOperator(
                                ((sp_p, 1), (sp_r, 1), (sp_s, 0), (sp_q, 0)),
                                coeff,
                            )

    return fermion_op


def fit_pauli_scaling(
    qubits: np.ndarray,
    pauli_counts: np.ndarray,
) -> Tuple[float, float]:
    """
    Fit power law: pauli_count = a * qubits^exponent.

    Parameters
    ----------
    qubits : array of qubit counts
    pauli_counts : array of Pauli term counts

    Returns
    -------
    exponent : float
        Power-law exponent.
    prefactor : float
        Prefactor a.
    """
    log_q = np.log(qubits.astype(float))
    log_p = np.log(pauli_counts.astype(float))
    exponent, log_a = np.polyfit(log_q, log_p, 1)
    return float(exponent), float(np.exp(log_a))
