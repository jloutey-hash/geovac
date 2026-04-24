"""
Native Dirac Lattice on S³
===========================

Spinorial analog of GeometricLattice. Nodes are Dirac states
|n, κ, m_j⟩ on the unit three-sphere S³ = SU(2). Spin is geometric,
not bolted on — it comes from the spinor bundle over S³.

The scalar lattice (Paper 7) uses nodes (n, l, m) — the scalar
harmonics on S³. This lattice uses the FULL harmonic content of
S³ = SU(2), including the spinor sector. The chirality ℤ₂ (the
sign of the Dirac eigenvalue) pairs states geometrically: two
nodes with the same orbital content but opposite chirality are
the spin-up and spin-down in the geometric sense.

Adjacency: E1 dipole selection rules (Rule B from ihara_zeta_dirac):
  - Δn ∈ {-1, 0, +1}
  - Δl = ±1 (parity flip)
  - |Δj| ≤ 1
  - |Δm_j| ≤ 1

Diagonal weights: Camporesi-Higuchi Dirac eigenvalues |λ_n| = n + 3/2
on unit S³, with chirality sign χ = ±1.

The Hamiltonian is H = D + W where D encodes the Dirac adjacency
and W encodes the diagonal spectrum. The graph IS the Dirac equation.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse import lil_matrix, csr_matrix
from typing import List, Dict, Tuple, Optional
from sympy import Rational

from dataclasses import dataclass
from sympy import Rational as SympyRational

from geovac.dirac_matrix_elements import (
    DiracLabel,
    kappa_to_l,
    kappa_to_j,
    iter_dirac_labels,
)


@dataclass(frozen=True)
class S3SpinorLabel:
    """Spinor harmonic label on S³, including boundary states with l = n_fock.

    The atomic DiracLabel enforces l < n_fock (hydrogen convention).
    S³ spinor harmonics include negative-chirality states at l = n_fock
    that are geometrically native to the three-sphere but do not
    correspond to atomic orbitals.

    These boundary states carry κ > 0 (negative chirality) with
    l = n_fock, j = l - 1/2 = n_fock - 1/2.
    """
    n_fock: int
    kappa: int
    two_m_j: int

    @property
    def l(self) -> int:
        return kappa_to_l(self.kappa)

    @property
    def j_times_2(self) -> int:
        return 2 * abs(self.kappa) - 1

    @property
    def j(self) -> SympyRational:
        return SympyRational(self.j_times_2, 2)

    @property
    def m_j(self) -> SympyRational:
        return SympyRational(self.two_m_j, 2)

    @property
    def is_boundary(self) -> bool:
        """True if this is a boundary state (l = n_fock), absent in atomic convention."""
        return self.l == self.n_fock

    def __post_init__(self):
        if self.n_fock < 1:
            raise ValueError(f"n_fock must be >= 1, got {self.n_fock}")
        if self.kappa == 0:
            raise ValueError("κ = 0 is not allowed")
        two_j = 2 * abs(self.kappa) - 1
        if abs(self.two_m_j) > two_j:
            raise ValueError(f"|two_m_j|={abs(self.two_m_j)} exceeds 2j={two_j}")
        if (self.two_m_j - two_j) % 2 != 0:
            raise ValueError(f"two_m_j parity mismatch")
        l = kappa_to_l(self.kappa)
        if l > self.n_fock:
            raise ValueError(f"l={l} exceeds n_fock={self.n_fock} (max for S³ is l <= n_fock)")
        if l == self.n_fock and self.kappa < 0:
            raise ValueError(f"Boundary state l=n_fock must have κ > 0 (negative chirality)")


def iter_s3_spinor_labels(n_max: int):
    """Generate ALL S³ spinor harmonic labels up to n_max.

    Includes both:
    - Atomic states (l < n_fock): from iter_dirac_labels
    - Boundary states (l = n_fock, κ > 0): geometrically native to S³

    Yields S3SpinorLabel in canonical order.
    """
    for n in range(1, n_max + 1):
        for l in range(n + 1):  # l = 0, ..., n (not n-1!)
            if l == n:
                # Boundary: only κ = +l (negative chirality, j = l - 1/2)
                if l == 0:
                    continue  # κ = 0 forbidden
                kappa = l
                two_j = 2 * l - 1
                for two_mj in range(-two_j, two_j + 1, 2):
                    yield S3SpinorLabel(n, kappa, two_mj)
            else:
                # Atomic: both κ values
                # κ = -(l+1): j = l + 1/2 (positive chirality)
                kappa_neg = -(l + 1)
                two_j_neg = 2 * (l + 1) - 1
                for two_mj in range(-two_j_neg, two_j_neg + 1, 2):
                    yield S3SpinorLabel(n, kappa_neg, two_mj)
                # κ = +l: j = l - 1/2 (negative chirality), only if l >= 1
                if l >= 1:
                    kappa_pos = l
                    two_j_pos = 2 * l - 1
                    for two_mj in range(-two_j_pos, two_j_pos + 1, 2):
                        yield S3SpinorLabel(n, kappa_pos, two_mj)


class DiracLattice:
    """
    Native spinorial lattice on S³.

    Nodes are Dirac states |n, κ, m_j⟩. Spin is encoded in κ:
        κ < 0  →  j = l + 1/2  (spin aligned with orbital AM)
        κ > 0  →  j = l - 1/2  (spin anti-aligned)

    Each node carries a chirality χ = ±1, the geometric ℤ₂ that
    pairs spin states. For the Camporesi-Higuchi spectrum on S³,
    positive-chirality states have eigenvalue +(n+3/2) and
    negative-chirality have -(n+3/2).

    Attributes
    ----------
    n_max : int
        Maximum Fock principal quantum number.
    labels : List[DiracLabel]
        All Dirac states in canonical order.
    adjacency : csr_matrix
        Sparse adjacency matrix (E1 dipole selection rules).
    chirality : np.ndarray
        Chirality sign (+1 or -1) for each node.
    dirac_eigenvalues : np.ndarray
        Camporesi-Higuchi eigenvalue λ = χ · (n + 3/2) per node.
    """

    def __init__(self, n_max: int, nuclear_charge: int = 1, mode: str = 'atomic'):
        if n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {n_max}")
        if mode not in ('atomic', 's3'):
            raise ValueError(f"mode must be 'atomic' or 's3', got {mode!r}")

        self.n_max = n_max
        self.nuclear_charge = nuclear_charge
        self.mode = mode

        self.labels: list = []
        self._label_index: Dict = {}
        self.chirality: np.ndarray = None
        self.dirac_eigenvalues: np.ndarray = None
        self.adjacency: csr_matrix = None

        self._generate_states()
        self._assign_chirality()
        self._build_adjacency()
        self._compute_dirac_eigenvalues()

    @property
    def num_states(self) -> int:
        return len(self.labels)

    @property
    def num_edges(self) -> int:
        return self.adjacency.nnz // 2

    @property
    def scalar_degeneracy(self) -> int:
        """Number of scalar (n,l,m) states up to n_max: Σ n² = n_max(n_max+1)(2n_max+1)/6."""
        n = self.n_max
        return n * (n + 1) * (2 * n + 1) // 6

    @property
    def s3_dirac_degeneracy(self) -> int:
        """Full S³ Dirac state count: Σ 2n(n+1) = 2n_max(n_max+1)(n_max+2)/3."""
        n = self.n_max
        return 2 * n * (n + 1) * (n + 2) // 3

    def sparsity(self) -> float:
        n = self.num_states
        if n == 0:
            return 1.0
        return 1.0 - self.adjacency.nnz / (n * n)

    def _generate_states(self) -> None:
        """Generate all Dirac labels up to n_max in canonical order.

        mode='atomic': l < n_fock only (hydrogen convention, 2n² per shell).
        mode='s3': full S³ spinor harmonics including l = n_fock boundary
                   states (2n(n+1) per shell, matching Camporesi-Higuchi).
        """
        if self.mode == 's3':
            self.labels = list(iter_s3_spinor_labels(self.n_max))
        else:
            self.labels = list(iter_dirac_labels(self.n_max))
        self._label_index = {lab: i for i, lab in enumerate(self.labels)}

    def _assign_chirality(self) -> None:
        """Assign chirality χ = ±1 to each node.

        Convention: κ < 0 (j = l + 1/2, "spin-up") gets χ = +1,
                    κ > 0 (j = l - 1/2, "spin-down") gets χ = -1.

        This is the geometric ℤ₂ from the spinor bundle: the two
        chiralities are related by orientation reversal on S³.
        """
        self.chirality = np.array(
            [-1 if lab.kappa > 0 else +1 for lab in self.labels],
            dtype=np.int8,
        )

    def _build_adjacency(self) -> None:
        """Build E1 dipole adjacency (Rule B from ihara_zeta_dirac).

        Selection rules:
          - Δn_fock ∈ {-1, 0, +1}
          - Δl = ±1 (parity flip, required for E1)
          - |Δj| ≤ 1
          - |Δm_j| ≤ 1
        """
        n = self.num_states
        adj = lil_matrix((n, n), dtype=np.float64)

        for i, a in enumerate(self.labels):
            l_a = a.l
            j2_a = a.j_times_2
            for j_idx, b in enumerate(self.labels):
                if j_idx <= i:
                    continue
                dn = abs(b.n_fock - a.n_fock)
                if dn > 1:
                    continue
                l_b = b.l
                if abs(l_b - l_a) != 1:
                    continue
                j2_b = b.j_times_2
                if abs(j2_b - j2_a) > 2:
                    continue
                dm2 = abs(b.two_m_j - a.two_m_j)
                if dm2 > 2:
                    continue
                adj[i, j_idx] = 1.0
                adj[j_idx, i] = 1.0

        self.adjacency = adj.tocsr()

    def _compute_dirac_eigenvalues(self) -> None:
        """Camporesi-Higuchi eigenvalues: λ = χ · (n_fock + 3/2).

        The unsigned magnitude |λ| = n + 3/2 is the Dirac eigenvalue
        on unit S³. The sign comes from chirality.
        """
        self.dirac_eigenvalues = np.array(
            [float(self.chirality[i]) * (lab.n_fock + 1.5)
             for i, lab in enumerate(self.labels)],
            dtype=np.float64,
        )

    def fock_weights(self) -> np.ndarray:
        """Diagonal Fock potential weights: -Z/n² per node.

        The scalar lattice uses these as the full potential. In the
        Dirac lattice, they serve as the non-relativistic limit.
        """
        Z = self.nuclear_charge
        return np.array(
            [-Z / (lab.n_fock ** 2) for lab in self.labels],
            dtype=np.float64,
        )

    def degree_sequence(self) -> np.ndarray:
        """Number of neighbors per node."""
        return np.array(self.adjacency.sum(axis=1)).flatten().astype(int)

    def per_kappa_blocks(self) -> Dict[int, List[int]]:
        """Group node indices by κ value."""
        blocks: Dict[int, List[int]] = {}
        for i, lab in enumerate(self.labels):
            blocks.setdefault(lab.kappa, []).append(i)
        return blocks

    def chirality_partner(self, label):
        """Find the chirality partner of a state.

        For (n, κ, m_j), the partner has the same n and l but
        opposite κ sign (flipped j = l ± 1/2), with the closest
        valid m_j.

        Returns None if no partner exists (e.g., l=0 has only κ=-1).
        """
        l = label.l
        if label.kappa < 0:
            partner_kappa = l if l > 0 else None
        else:
            partner_kappa = -(l + 1)

        if partner_kappa is None or partner_kappa == 0:
            return None

        two_j_partner = 2 * abs(partner_kappa) - 1
        two_mj = label.two_m_j
        if abs(two_mj) > two_j_partner:
            two_mj = two_j_partner if two_mj > 0 else -two_j_partner
        if (two_mj - two_j_partner) % 2 != 0:
            two_mj += 1 if two_mj < two_j_partner else -1

        try:
            if self.mode == 's3':
                partner = S3SpinorLabel(label.n_fock, partner_kappa, two_mj)
            else:
                partner = DiracLabel(label.n_fock, partner_kappa, two_mj)
        except ValueError:
            return None

        return partner if partner in self._label_index else None

    def shell_structure(self) -> Dict[int, Dict]:
        """Summary of each n-shell: state count, κ values, edge count."""
        shells = {}
        for n in range(1, self.n_max + 1):
            indices = [i for i, lab in enumerate(self.labels) if lab.n_fock == n]
            kappas = sorted(set(self.labels[i].kappa for i in indices))
            sub_adj = self.adjacency[np.ix_(indices, indices)]
            shells[n] = {
                "num_states": len(indices),
                "kappa_values": kappas,
                "internal_edges": sub_adj.nnz // 2,
            }
        return shells

    def __repr__(self) -> str:
        return (
            f"DiracLattice(n_max={self.n_max}, "
            f"states={self.num_states}, "
            f"edges={self.num_edges}, "
            f"Z={self.nuclear_charge})"
        )
