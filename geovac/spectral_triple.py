"""
Spectral triple on the Fock-projected S^3 graph.
=================================================

Constructs the data (A, H, D) of a real spectral triple in the sense
of Connes on the finite GeoVac Dirac graph at a given n_max.

Mathematical structure
----------------------
A spectral triple is a triple (A, H, D) where:
  - H is a (finite-dimensional) Hilbert space.
  - A is a *-algebra represented on H (here: commutative, = functions on
    a finite set of nodes or sectors).
  - D is a self-adjoint operator on H (the Dirac operator).

Additional structures:
  - gamma: a Z/2 grading (chirality), diagonal with chi = +/-1.
  - J: a real structure (antilinear isometry), here a real matrix C
    since all entries are real.

KO-dimension and odd spectral triples
--------------------------------------
S^3 is 3-dimensional, so the spectral triple is ODD.  In an odd
spectral triple, the grading gamma exists as additional structure
(from the volume form on S^3) but is NOT required to anticommute
with D.  In odd dimensions, the volume form COMMUTES with D
(see Connes, NCG Ch. VI; Gracia-Bondia, Varilly, Figueroa Ch. 11).

On the FINITE Dirac graph, even the commutation [D, gamma] = 0 fails
because the E1 dipole adjacency connects states of BOTH same and
opposite chirality (Delta l = +/-1 flips parity, but kappa can be
either -(l+1) or +l for the same l, giving both chirality signs).
The diagonal part Lambda COMMUTES with gamma (since Lambda and gamma
are both diagonal with chi*|lam| and chi respectively).

Real structure J
----------------
Two J constructions are provided:

1. **Permutation J** (j_type='permutation'): pure m_j reversal,
   J: |n, kappa, m_j> -> |n, kappa, -m_j>.  This is a permutation
   matrix with J^2 = +I.  It commutes with D (JD = DJ) because the
   E1 selection rules are symmetric under m_j -> -m_j.  This gives
   a real spectral triple with epsilon = +1.

2. **Kramers J** (j_type='kramers'): time-reversal with Kramers phase,
   J: |n, kappa, m_j> -> (-1)^{j-m_j} |n, kappa, -m_j>.  This gives
   J^2 = -I (since (-1)^{2j} = -1 for half-integer j), the quaternionic
   structure expected for S^3 = SU(2).  However, JD != +/-DJ on the
   finite graph because the uniform adjacency weights do not respect
   the m_j-dependent Kramers phases.  This J requires CG-weighted
   adjacency to achieve exact JD = -DJ (flagged as future work).

Default is j_type='permutation' (the one that satisfies all axioms).

Order-one condition
-------------------
The order-one condition [[D, pi(a)], J pi(b) J^{-1}] = 0 FAILS for
the full sector algebra because D connects different sectors.  This is
a structural feature of finite spectral triples where the Dirac operator
has long-range hops between algebra points.  In the Connes-van Suijlekom
spectral truncation framework (CMP 2021, Hekkelman 2022+2024), the
order-one condition is relaxed for truncated spectral triples.

The order-zero condition [pi(a), J pi(b) J^{-1}] = 0 HOLDS for the
permutation J because J preserves sectors and A is commutative.

Implementation
--------------
- H: full Dirac states on S^3 at finite n_max, from DiracLattice in
  's3' mode.  Dimension = sum_{n=1}^{n_max} 2*n*(n+1).
  At n_max=2: dim=16; at n_max=3: dim=40 (= Delta^{-1} from Paper 2).

- A: functions on the scalar Fock graph nodes (n, l, m) at the same
  n_max.  The representation pi: A -> B(H) sends function f to the
  diagonal matrix where each Dirac state |n, kappa, m_j> gets the value
  f at its orbital node (n, l).  When the algebra is defined on (n,l)
  sectors (the default), all m_j states in the same (n,l) sector share
  the same function value.

- D: the Dirac graph operator D = Lambda + kappa * A_graph where
  Lambda carries the chirality-signed Camporesi-Higuchi eigenvalues
  chi * (n_fock + 1/2) and A_graph is the E1 dipole adjacency from
  DiracLattice.  All entries are exact sympy Rationals.  kappa = -1/16
  is the universal topological constant from Paper 0.

- gamma: the chirality grading, diagonal with chi_i in {+1, -1}
  (kappa < 0 -> chi = +1; kappa > 0 -> chi = -1).

All matrices are built in exact sympy arithmetic (Rational, Integer).
No floating-point numbers are used in any construction.

References
----------
- A. Connes, "Noncommutative Geometry" (1994), Ch. VI.
- A. Connes, "Geometry from the spectral point of view",
  Lett. Math. Phys. 34 (1995) 203-238.
- J. W. Barrett, "Matrix geometries and fuzzy spaces as finite
  spectral triples", J. Math. Phys. 56 (2015) 082301.
- Connes, van Suijlekom, "Spectral truncations in noncommutative
  geometry and operator systems", Commun. Math. Phys. 383 (2021).
- GeoVac Paper 25 (Hopf gauge structure on the Fock graph).
- GeoVac CLAUDE.md section 1.7 WH1 (almost-commutative spectral triple).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import sympy as sp
from sympy import Integer, Matrix, Rational, eye as sp_eye, sqrt as sp_sqrt, zeros as sp_zeros
from sympy.matrices.exceptions import MatrixError as SympyMatrixError
from sympy.physics.wigner import wigner_3j

from geovac.dirac_lattice import DiracLattice, S3SpinorLabel
from geovac.dirac_matrix_elements import kappa_to_l, kappa_to_j

__all__ = ["FockSpectralTriple"]


class FockSpectralTriple:
    """Spectral triple (A, H, D) on the Fock-projected S^3 Dirac graph.

    Parameters
    ----------
    n_max : int
        Maximum Fock principal quantum number (>= 1).  The Dirac
        Hilbert space dimension is sum_{n=1}^{n_max} 2*n*(n+1).
        n_max=2 gives dim=16; n_max=3 gives dim=40.
    kappa : sympy.Rational, optional
        The off-diagonal hopping strength in D = Lambda + kappa * A.
        Default: Rational(-1, 16) (the universal topological constant
        from Paper 0).
    j_type : str, optional
        Type of real structure J to build.  'permutation' (default)
        gives J^2 = +I and JD = DJ.  'kramers' gives J^2 = -I but
        JD relation does not hold cleanly on the finite graph with
        uniform adjacency.
    adjacency_weights : str, optional
        Edge weighting scheme for the E1 dipole adjacency in D.
        'uniform' (default): all edges carry weight kappa (Paper 0).
        'cg': edges carry kappa * w_{ij} where w_{ij} is the Wigner 3j
        CG coupling coefficient for the E1 (q=1) transition between
        Dirac states i and j.  CG weights are exact sympy expressions
        (algebraic, pi-free).

    Attributes
    ----------
    dim_H : int
        Dimension of the Hilbert space H.
    n_sectors : int
        Number of (n, l) orbital sectors (= dim of the algebra A).
    sectors : list of (int, int)
        The (n_fock, l) sector labels in canonical order.
    labels : list of S3SpinorLabel
        The Dirac state labels in canonical order.
    """

    def __init__(
        self,
        n_max: int = 3,
        kappa: object = Rational(-1, 16),
        j_type: str = "permutation",
        adjacency_weights: str = "uniform",
    ):
        if n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {n_max}")
        if j_type not in ("permutation", "kramers"):
            raise ValueError(f"j_type must be 'permutation' or 'kramers', got {j_type!r}")
        if adjacency_weights not in ("uniform", "cg"):
            raise ValueError(
                f"adjacency_weights must be 'uniform' or 'cg', got {adjacency_weights!r}"
            )

        self.n_max = n_max
        self._kappa = Rational(kappa) if not isinstance(kappa, sp.Expr) else kappa
        self._j_type = j_type
        self._adjacency_weights = adjacency_weights

        # Build the underlying Dirac lattice in S3 mode
        self._lattice = DiracLattice(n_max=n_max, mode="s3")
        self.labels: List[S3SpinorLabel] = self._lattice.labels
        self.dim_H: int = len(self.labels)

        # Build the sector structure: (n_fock, l) pairs
        self._sector_map: Dict[Tuple[int, int], int] = {}
        self.sectors: List[Tuple[int, int]] = []
        self._state_to_sector: List[int] = []
        self._build_sectors()
        self.n_sectors: int = len(self.sectors)

        # Pre-compute CG weight matrix if requested
        self._cg_weights: Optional[Matrix] = None
        if adjacency_weights == "cg":
            self._cg_weights = self._build_cg_weight_matrix()

        # Cached matrices (built lazily)
        self._D: Optional[Matrix] = None
        self._gamma: Optional[Matrix] = None
        self._J: Optional[Matrix] = None
        self._Lambda: Optional[Matrix] = None

    # ------------------------------------------------------------------
    # Sector structure
    # ------------------------------------------------------------------

    def _build_sectors(self) -> None:
        """Build the (n, l) sector decomposition of the Hilbert space."""
        seen_order: List[Tuple[int, int]] = []
        for lab in self.labels:
            key = (lab.n_fock, lab.l)
            if key not in self._sector_map:
                self._sector_map[key] = len(seen_order)
                seen_order.append(key)
        self.sectors = seen_order

        self._state_to_sector = []
        for lab in self.labels:
            key = (lab.n_fock, lab.l)
            self._state_to_sector.append(self._sector_map[key])

    # ------------------------------------------------------------------
    # CG-weighted adjacency
    # ------------------------------------------------------------------

    def _build_cg_weight_matrix(self) -> Matrix:
        """Build the CG-weighted adjacency matrix using Wigner 3j symbols.

        For each E1 dipole edge (i, j) in the adjacency, the weight is
        the CG coupling coefficient for an E1 (q=1) transition:

            w(a, b) = (-1)^{j_a - m_j_a} * 3j(j_a, 1, j_b; -m_j_a, m_q, m_j_b)

        where m_q = m_j_b - m_j_a and |m_q| <= 1.

        The 3j symbol is computed in exact sympy arithmetic; the result
        is a rational number or zero.  No overall normalization prefactor
        is applied — the raw 3j serves as the dimensionless edge weight.
        The matrix is symmetric because the 3j satisfies:
            3j(j_a,1,j_b; -m_a,q,m_b) = (-1)^{j_a+1+j_b} * 3j(j_b,1,j_a; -m_b,-q,m_a)
        and the phase/sign is absorbed by the (-1)^{j-m} factor.

        Returns
        -------
        sympy.Matrix
            N x N matrix of CG weights (zero where no edge exists).
        """
        N = self.dim_H
        W = sp_zeros(N, N)

        adj = self._lattice.adjacency.tocoo()
        edge_set = set()
        for r, c, v in zip(adj.row, adj.col, adj.data):
            if r < c and v != 0:
                edge_set.add((int(r), int(c)))

        for (i, j) in edge_set:
            lab_a = self.labels[i]
            lab_b = self.labels[j]

            j_a = sp.Rational(2 * abs(lab_a.kappa) - 1, 2)
            j_b = sp.Rational(2 * abs(lab_b.kappa) - 1, 2)
            m_a = sp.Rational(lab_a.two_m_j, 2)
            m_b = sp.Rational(lab_b.two_m_j, 2)

            m_q = m_b - m_a
            if abs(m_q) > 1:
                continue

            # Check j-triangle for q=1
            if 1 < abs(j_a - j_b) or 1 > j_a + j_b:
                continue

            # Wigner 3j: (j_a, 1, j_b; -m_a, m_q, m_b)
            threej_val = wigner_3j(j_a, 1, j_b, -m_a, m_q, m_b)
            if threej_val == 0:
                continue

            # Phase: (-1)^{j_a - m_a}
            phase_exp = j_a - m_a  # always integer for half-integer j,m
            phase = Integer(-1) ** int(phase_exp)

            weight = phase * threej_val

            W[i, j] = weight
            W[j, i] = weight

        return W

    # ------------------------------------------------------------------
    # Hilbert space H
    # ------------------------------------------------------------------

    def hilbert_space_summary(self) -> Dict:
        """Summary of the Hilbert space H."""
        n_pos = sum(1 for c in self._lattice.chirality if c > 0)
        n_neg = sum(1 for c in self._lattice.chirality if c < 0)
        return {
            "dim_H": self.dim_H,
            "n_max": self.n_max,
            "n_sectors": self.n_sectors,
            "sectors": self.sectors,
            "chi_plus": n_pos,
            "chi_minus": n_neg,
        }

    # ------------------------------------------------------------------
    # Algebra A and representation pi
    # ------------------------------------------------------------------

    def algebra_representation(self, f: List) -> Matrix:
        """Build pi(f): the representation of f in A on H.

        Parameters
        ----------
        f : list of sympy expressions
            A function on the sectors, f[k] = value at sector k.
            Length must equal self.n_sectors.

        Returns
        -------
        sympy.Matrix
            A dim_H x dim_H diagonal matrix with f(sector(i)) on the
            diagonal for each state i.
        """
        if len(f) != self.n_sectors:
            raise ValueError(
                f"f must have length {self.n_sectors} (one per sector), "
                f"got {len(f)}"
            )
        N = self.dim_H
        M = sp_zeros(N, N)
        for i in range(N):
            M[i, i] = f[self._state_to_sector[i]]
        return M

    # ------------------------------------------------------------------
    # Dirac operator D
    # ------------------------------------------------------------------

    @property
    def dirac_operator(self) -> Matrix:
        """The Dirac operator D = Lambda + kappa * A_graph.

        Lambda is diagonal with chirality-signed Camporesi-Higuchi
        eigenvalues chi * (n_fock + 1/2) as exact Rationals.

        A_graph is the E1 dipole adjacency from DiracLattice (s3 mode).

        Returns an exact sympy Matrix of dimension dim_H x dim_H.
        """
        if self._D is not None:
            return self._D

        N = self.dim_H
        D = sp_zeros(N, N)

        # Diagonal: Lambda
        for i, lab in enumerate(self.labels):
            chi = Integer(int(self._lattice.chirality[i]))
            lam_abs = Rational(2 * lab.n_fock + 1, 2)
            D[i, i] = chi * lam_abs

        # Off-diagonal: kappa * adjacency (uniform or CG-weighted)
        if self._kappa != Rational(0):
            if self._adjacency_weights == "cg" and self._cg_weights is not None:
                # CG-weighted: D_{ij} = kappa * w_{ij} for each edge
                for i in range(N):
                    for j in range(N):
                        if i != j and self._cg_weights[i, j] != 0:
                            D[i, j] = self._kappa * self._cg_weights[i, j]
            else:
                # Uniform: all edges get weight kappa
                adj = self._lattice.adjacency.tocoo()
                for r, c, v in zip(adj.row, adj.col, adj.data):
                    if r != c and v != 0:
                        D[int(r), int(c)] = self._kappa

        self._D = D
        return D

    @property
    def diagonal_part(self) -> Matrix:
        """The diagonal part Lambda of D (Camporesi-Higuchi eigenvalues).

        Lambda[i, i] = chi_i * (n_fock_i + 1/2).
        """
        if self._Lambda is not None:
            return self._Lambda

        N = self.dim_H
        Lam = sp_zeros(N, N)
        for i, lab in enumerate(self.labels):
            chi = Integer(int(self._lattice.chirality[i]))
            lam_abs = Rational(2 * lab.n_fock + 1, 2)
            Lam[i, i] = chi * lam_abs
        self._Lambda = Lam
        return Lam

    # ------------------------------------------------------------------
    # Chirality grading gamma
    # ------------------------------------------------------------------

    @property
    def grading(self) -> Matrix:
        """The chirality grading gamma (Z/2-grading on H).

        gamma[i, i] = chi_i, where chi = +1 for kappa < 0 (positive
        chirality) and chi = -1 for kappa > 0 (negative chirality).

        Properties:
          - gamma^2 = I
          - gamma is self-adjoint (gamma = gamma^dagger)

        Note: {gamma, D} = 0 does NOT hold because this is an odd-
        dimensional spectral triple.  The diagonal Lambda commutes
        with gamma, and the adjacency connects both same-chirality
        and opposite-chirality states.
        """
        if self._gamma is not None:
            return self._gamma

        N = self.dim_H
        g = sp_zeros(N, N)
        for i in range(N):
            g[i, i] = Integer(int(self._lattice.chirality[i]))
        self._gamma = g
        return g

    # ------------------------------------------------------------------
    # Real structure J
    # ------------------------------------------------------------------

    @property
    def real_structure(self) -> Matrix:
        """The real structure J as a real matrix C.

        J is an antilinear isometry J = C * K where K is complex
        conjugation.  Since all our matrices are real, K acts trivially,
        so J = C.

        For j_type='permutation':
          C maps |n, kappa, m_j> to |n, kappa, -m_j> with no phase.
          C is a permutation matrix with C^2 = +I, C^T = C.
          CDE = DC (J commutes with D).

        For j_type='kramers':
          C maps |n, kappa, m_j> to (-1)^{j-m_j} |n, kappa, -m_j>.
          C^2 = -I (quaternionic, from (-1)^{2j} = -1).
          CD + DC != 0 on the finite graph (phases mismatch adjacency).
        """
        if self._J is not None:
            return self._J

        N = self.dim_H
        C = sp_zeros(N, N)
        label_to_idx = {lab: i for i, lab in enumerate(self.labels)}

        for i, lab in enumerate(self.labels):
            partner_two_mj = -lab.two_m_j
            partner = S3SpinorLabel(lab.n_fock, lab.kappa, partner_two_mj)

            if partner not in label_to_idx:
                raise RuntimeError(
                    f"Time-reversal partner of state {lab} not found."
                )
            j_idx = label_to_idx[partner]

            if self._j_type == "kramers":
                two_j = 2 * abs(lab.kappa) - 1
                j_minus_mj_times_2 = two_j - lab.two_m_j
                j_minus_mj = j_minus_mj_times_2 // 2
                phase = Integer((-1) ** j_minus_mj)
            else:
                phase = Integer(1)

            C[j_idx, i] = phase

        self._J = C
        return C

    # ------------------------------------------------------------------
    # Commutator [D, pi(f)]
    # ------------------------------------------------------------------

    def commutator(self, f: List) -> Matrix:
        """Compute [D, pi(f)] for a function f on the sectors.

        Parameters
        ----------
        f : list of sympy expressions
            Function values on the sectors.

        Returns
        -------
        sympy.Matrix
            The commutator [D, pi(f)] = D*pi(f) - pi(f)*D.
        """
        pi_f = self.algebra_representation(f)
        D = self.dirac_operator
        return D * pi_f - pi_f * D

    # ------------------------------------------------------------------
    # Connes distance
    # ------------------------------------------------------------------

    def connes_distance(self, sector_i: int, sector_j: int) -> sp.Expr:
        """Compute the Connes spectral distance between two sectors.

        d(i, j) = sup { |f(i) - f(j)| : ||[D, pi(f)]||_op <= 1 }

        For the indicator function f = e_i - e_j, this gives a LOWER
        BOUND on the Connes distance:

            d(i, j) >= 2 / ||[D, pi(e_i - e_j)]||_op

        The exact Connes distance for a finite commutative spectral
        triple requires solving a semidefinite program (optimization
        over all Lipschitz-1 functions).  This method returns the lower
        bound from the indicator function, computed via the operator
        norm (largest singular value).

        Parameters
        ----------
        sector_i, sector_j : int
            Indices into self.sectors.

        Returns
        -------
        sympy expression
            A lower bound on the Connes distance.
        """
        if sector_i == sector_j:
            return Integer(0)

        N_sec = self.n_sectors
        f = [Integer(0)] * N_sec
        f[sector_i] = Integer(1)
        f[sector_j] = Integer(-1)

        comm = self.commutator(f)

        # Compute operator norm numerically
        import numpy as np
        rows, cols = comm.shape
        arr = np.zeros((rows, cols), dtype=np.float64)
        for i in range(rows):
            for j in range(cols):
                arr[i, j] = float(comm[i, j])
        svds = np.linalg.svd(arr, compute_uv=False)
        op_norm = float(max(svds))

        if op_norm < 1e-15:
            return sp.oo

        return Rational(2) / sp.nsimplify(op_norm, rational=False)

    # ------------------------------------------------------------------
    # Axiom checks
    # ------------------------------------------------------------------

    def check_selfadjoint(self) -> bool:
        """Check that D = D^T (self-adjointness for real matrices)."""
        D = self.dirac_operator
        return D.equals(D.T)

    def check_grading_square(self) -> bool:
        """Check that gamma^2 = I."""
        g = self.grading
        return (g * g).equals(sp_eye(self.dim_H))

    def check_grading_selfadjoint(self) -> bool:
        """Check that gamma = gamma^T."""
        g = self.grading
        return g.equals(g.T)

    def check_D_gamma_anticommute(self) -> Tuple[bool, str]:
        """Check whether {D, gamma} = 0.

        Returns (result, explanation).  On the finite Dirac graph this
        is expected to FAIL because:
        1. Lambda commutes (not anticommutes) with gamma (both diagonal,
           chi^2 = 1 gives Lambda*gamma = gamma*Lambda).
        2. The adjacency connects both same-chirality and opposite-
           chirality states.
        """
        D = self.dirac_operator
        g = self.grading
        anti = D * g + g * D
        N = self.dim_H
        result = anti.equals(sp_zeros(N, N))
        if result:
            msg = "PASS: {D, gamma} = 0"
        else:
            n_nonzero = sum(1 for i in range(N) for j in range(N) if anti[i, j] != 0)
            msg = (
                f"EXPECTED FAIL: {{D, gamma}} has {n_nonzero} nonzero entries. "
                f"This is correct for an odd-dimensional spectral triple where "
                f"the volume-form chirality is not required to anticommute with D."
            )
        return result, msg

    def check_D_gamma_commute_diagonal(self) -> bool:
        """Check that the diagonal part of D commutes with gamma.

        [Lambda, gamma] = 0 holds because Lambda[i,i] = chi_i * |lam_i|
        and gamma[i,i] = chi_i are both diagonal, so their product
        commutes: chi_i * |lam_i| * chi_i = |lam_i| = chi_i * chi_i * |lam_i|.
        """
        Lam = self.diagonal_part
        g = self.grading
        return (Lam * g - g * Lam).equals(sp_zeros(self.dim_H, self.dim_H))

    def check_J_squared(self) -> Tuple[bool, int]:
        """Check J^2.

        Returns (pass, epsilon) where J^2 = epsilon * I.
        For j_type='permutation': epsilon = +1.
        For j_type='kramers': epsilon = -1.
        """
        C = self.real_structure
        N = self.dim_H
        C2 = C * C
        if C2.equals(sp_eye(N)):
            return True, 1
        elif C2.equals(-sp_eye(N)):
            return True, -1
        else:
            return False, 0

    def check_J_D_relation(self) -> Tuple[bool, str]:
        """Check the relation between J and D.

        For j_type='permutation': checks JD = DJ (commutation).
        For j_type='kramers': checks JD = -DJ (anticommutation).

        Returns (result, explanation).
        """
        C = self.real_structure
        D = self.dirac_operator
        N = self.dim_H

        if self._j_type == "permutation":
            comm = C * D - D * C
            result = comm.equals(sp_zeros(N, N))
            if result:
                msg = "PASS: JD = DJ (permutation J commutes with D)"
            else:
                n_nz = sum(1 for i in range(N) for j in range(N) if comm[i, j] != 0)
                msg = f"FAIL: [J, D] has {n_nz} nonzero entries"
        else:
            anti = C * D + D * C
            result = anti.equals(sp_zeros(N, N))
            if result:
                msg = "PASS: JD = -DJ (Kramers J anticommutes with D)"
            else:
                n_nz = sum(1 for i in range(N) for j in range(N) if anti[i, j] != 0)
                msg = (
                    f"EXPECTED FAIL: {{J, D}} has {n_nz} nonzero entries. "
                    f"Kramers J requires CG-weighted adjacency for exact "
                    f"anticommutation (flagged as future work)."
                )
        return result, msg

    def check_kramers_D_relation(self) -> Dict:
        """Characterize the Kramers JD + DJ residual.

        Computes R = JD + DJ (which should be zero for exact KO-dim 3
        anticommutation JD = -DJ).  Reports the Frobenius norm, number
        of nonzero entries, and whether R is exactly zero.

        This method works for ANY j_type, but is most informative for
        j_type='kramers'.

        Returns
        -------
        dict with keys:
            'exact_zero' : bool
                True if JD + DJ = 0 exactly.
            'n_nonzero' : int
                Number of nonzero entries in JD + DJ.
            'frobenius_norm' : float
                Frobenius norm of JD + DJ (computed numerically).
            'max_entry' : sympy expression
                Maximum absolute entry in JD + DJ.
            'residual_matrix' : sympy.Matrix
                The full JD + DJ matrix.
        """
        C = self.real_structure
        D = self.dirac_operator
        N = self.dim_H

        R = C * D + D * C

        n_nonzero = sum(1 for i in range(N) for j in range(N) if R[i, j] != 0)
        exact_zero = (n_nonzero == 0)

        # Frobenius norm: sqrt(sum |R_ij|^2)
        frob_sq = Integer(0)
        max_entry = Integer(0)
        for i in range(N):
            for j in range(N):
                val = R[i, j]
                frob_sq += val * val
                abs_val = sp.Abs(val)
                if abs_val > max_entry:
                    max_entry = abs_val

        return {
            "exact_zero": exact_zero,
            "n_nonzero": n_nonzero,
            "frobenius_norm": float(sp.sqrt(frob_sq)),
            "max_entry": max_entry,
            "residual_matrix": R,
        }

    def kramers_D_residual_analysis(self) -> Dict:
        """Per-sector and per-component analysis of the Kramers JD residual.

        Decomposes the residual R = JD + DJ into contributions from the
        diagonal part (JLambda + LambdaJ) and the off-diagonal part
        (J*kappa*A + kappa*A*J), and reports sector-resolved diagnostics.

        Returns
        -------
        dict with keys:
            'diagonal_residual' : dict
                {'exact_zero': bool, 'n_nonzero': int, 'frobenius_norm': float}
                for the diagonal part JLambda + LambdaJ.
            'offdiag_residual' : dict
                {'exact_zero': bool, 'n_nonzero': int, 'frobenius_norm': float}
                for the off-diagonal part J*(kappa*A) + (kappa*A)*J.
            'per_sector_max' : list of (sector_label, float)
                Maximum |R_{ij}| where state i belongs to each sector.
            'constraint_test' : dict
                Tests the necessary constraint {J, A} = -(2/kappa)*Lambda*J
                on the adjacency part.  Reports Frobenius residual.
        """
        C = self.real_structure
        D = self.dirac_operator
        Lam = self.diagonal_part
        N = self.dim_H

        # Build adjacency-only part of D
        A_part = D - Lam  # = kappa * A_graph (or kappa * W for CG-weighted)

        # Diagonal residual: J*Lambda + Lambda*J
        R_diag = C * Lam + Lam * C
        n_nz_diag = sum(1 for i in range(N) for j in range(N) if R_diag[i, j] != 0)
        frob_diag = float(sp.sqrt(sum(R_diag[i, j] ** 2 for i in range(N) for j in range(N))))

        # Off-diagonal residual: J*A_part + A_part*J
        R_offdiag = C * A_part + A_part * C
        n_nz_off = sum(1 for i in range(N) for j in range(N) if R_offdiag[i, j] != 0)
        frob_off = float(sp.sqrt(sum(R_offdiag[i, j] ** 2 for i in range(N) for j in range(N))))

        # Full residual
        R_full = R_diag + R_offdiag

        # Per-sector max
        per_sector_max = []
        for s_idx, s_label in enumerate(self.sectors):
            max_val = 0.0
            for i in range(N):
                if self._state_to_sector[i] == s_idx:
                    for j in range(N):
                        val = float(sp.Abs(R_full[i, j]))
                        if val > max_val:
                            max_val = val
            per_sector_max.append((s_label, max_val))

        # Constraint test: {J, A} should equal -(2/kappa)*Lambda*J
        # for exact JD = -DJ.
        # A_part = kappa * A, so {J, A} = {J, A_part}/kappa
        # Target: {J, A_part}/kappa = -(2/kappa)*Lambda*J
        # => {J, A_part} = -2*Lambda*J
        if self._kappa != 0:
            target = Integer(-2) * Lam * C
            constraint_residual = R_offdiag - target  # {J,A_part} - (-2*Lam*J)
            # But R_offdiag = {J, A_part} = J*A_part + A_part*J
            # We want {J, A_part} = -2*Lam*J
            # So residual = {J, A_part} + 2*Lam*J = R_offdiag + 2*Lam*J
            constraint_residual = R_offdiag + Integer(2) * Lam * C
            frob_constraint = float(sp.sqrt(
                sum(constraint_residual[i, j] ** 2 for i in range(N) for j in range(N))
            ))
            constraint_zero = all(
                constraint_residual[i, j] == 0 for i in range(N) for j in range(N)
            )
        else:
            frob_constraint = 0.0
            constraint_zero = True

        return {
            "diagonal_residual": {
                "exact_zero": n_nz_diag == 0,
                "n_nonzero": n_nz_diag,
                "frobenius_norm": frob_diag,
            },
            "offdiag_residual": {
                "exact_zero": n_nz_off == 0,
                "n_nonzero": n_nz_off,
                "frobenius_norm": frob_off,
            },
            "per_sector_max": per_sector_max,
            "constraint_test": {
                "exact_zero": constraint_zero,
                "frobenius_norm": frob_constraint,
            },
        }

    def check_J_gamma_commute(self) -> bool:
        """Check that J*gamma = gamma*J.

        This holds for both J constructions because J maps states
        within the same (n, kappa) sector, and gamma depends only on
        the sign of kappa.
        """
        C = self.real_structure
        g = self.grading
        return (C * g - g * C).equals(sp_zeros(self.dim_H, self.dim_H))

    def check_order_zero(self) -> bool:
        """Check the order-zero condition: [pi(a), J pi(b) J^{-1}] = 0.

        For the permutation J (which preserves sectors) and commutative A,
        this is automatic: J pi(b) J^{-1} = pi(b), so [pi(a), pi(b)] = 0.
        """
        C = self.real_structure
        N = self.dim_H

        if self._j_type == "permutation":
            # J is an involution: J^{-1} = J
            C_inv = C
        else:
            # J^2 = -I => J^{-1} = -J
            C_inv = -C

        for a_idx in range(self.n_sectors):
            f_a = [Integer(0)] * self.n_sectors
            f_a[a_idx] = Integer(1)
            pi_a = self.algebra_representation(f_a)

            for b_idx in range(self.n_sectors):
                f_b = [Integer(0)] * self.n_sectors
                f_b[b_idx] = Integer(1)
                pi_b = self.algebra_representation(f_b)

                conj = C * pi_b * C_inv
                comm = pi_a * conj - conj * pi_a
                if not comm.equals(sp_zeros(N, N)):
                    return False
        return True

    def check_order_one(self) -> Tuple[bool, str]:
        """Check the order-one condition: [[D, pi(a)], J pi(b) J^{-1}] = 0.

        Returns (result, explanation).  This is expected to FAIL
        because D connects different sectors.  The double commutator
        [[D, pi(a)], pi(b)]_{ij} = D_{ij} * (a_{s(j)}-a_{s(i)}) * (b_{s(j)}-b_{s(i)})
        is nonzero whenever D connects states in different sectors
        and a, b separate those sectors.
        """
        D = self.dirac_operator
        C = self.real_structure
        N = self.dim_H

        if self._j_type == "permutation":
            C_inv = C
        else:
            C_inv = -C

        failures = 0
        for a_idx in range(self.n_sectors):
            f_a = [Integer(0)] * self.n_sectors
            f_a[a_idx] = Integer(1)
            pi_a = self.algebra_representation(f_a)
            comm_D_a = D * pi_a - pi_a * D

            for b_idx in range(self.n_sectors):
                f_b = [Integer(0)] * self.n_sectors
                f_b[b_idx] = Integer(1)
                pi_b = self.algebra_representation(f_b)

                conj = C * pi_b * C_inv
                double_comm = comm_D_a * conj - conj * comm_D_a
                if not double_comm.equals(sp_zeros(N, N)):
                    failures += 1

        total = self.n_sectors ** 2
        if failures == 0:
            return True, "PASS: order-one holds for all algebra pairs"
        else:
            return False, (
                f"EXPECTED FAIL: order-one fails for {failures}/{total} "
                f"algebra pairs. D connects different sectors, so the "
                f"double commutator is generically nonzero. This is a "
                f"standard feature of spectral truncations (Connes-vS 2021)."
            )

    # ------------------------------------------------------------------
    # Spectral action
    # ------------------------------------------------------------------

    def spectral_action(self, f_cutoff=None) -> sp.Expr:
        """Compute the spectral action Tr(f(D^2 / Lambda^2)).

        If f_cutoff is None, return Tr(I) = dim(H).
        If f_cutoff is a number, count eigenvalues with |lambda_i| <= cutoff.

        Parameters
        ----------
        f_cutoff : number, optional
            The energy cutoff Lambda.

        Returns
        -------
        sympy expression
        """
        if f_cutoff is None:
            return Integer(self.dim_H)

        cutoff_sq = Rational(f_cutoff) ** 2
        D = self.dirac_operator
        eigenvalues_dict = D.eigenvals()

        count = Integer(0)
        for ev, mult in eigenvalues_dict.items():
            if ev ** 2 <= cutoff_sq:
                count += Integer(mult)
        return count

    def spectral_action_heat(self, t_param: object = Rational(1)) -> sp.Expr:
        """Compute the heat-kernel spectral action Tr(exp(-t * D^2)).

        Parameters
        ----------
        t_param : sympy expression
            The "time" parameter t > 0.

        Returns
        -------
        sympy expression or float
            Exact sympy expression for uniform weights; float for CG weights
            (where exact eigenvalue computation may fail over algebraic fields).
        """
        t = sp.sympify(t_param)
        if self._adjacency_weights == "cg":
            # Numeric path for CG weights (exact eigenvalues unavailable)
            t_float = float(t)
            D_num = np.array(self.dirac_operator.tolist(), dtype=float)
            evals = np.linalg.eigvalsh(D_num)
            return float(np.sum(np.exp(-t_float * evals**2)))
        D = self.dirac_operator
        eigenvalues_dict = D.eigenvals()
        result = Integer(0)
        for ev, mult in eigenvalues_dict.items():
            result += Integer(mult) * sp.exp(-t * ev ** 2)
        return result

    # ------------------------------------------------------------------
    # Spectrum of D
    # ------------------------------------------------------------------

    def dirac_spectrum(self) -> Dict[sp.Expr, int]:
        """Eigenvalues and multiplicities of D.

        Returns exact sympy eigenvalues for uniform weights.
        For CG weights, uses numeric eigenvalues (numpy) directly
        since sympy cannot solve characteristic polynomials over
        the algebraic extension field generated by Wigner 3j symbols.
        """
        if self._adjacency_weights == "cg":
            return self._dirac_spectrum_numeric()
        try:
            return self.dirac_operator.eigenvals()
        except (SympyMatrixError, NotImplementedError):
            return self._dirac_spectrum_numeric()

    def _dirac_spectrum_numeric(self) -> Dict[sp.Expr, int]:
        """Numeric eigenvalues with multiplicity grouping (tolerance-based)."""
        D_num = np.array(self.dirac_operator.tolist(), dtype=float)
        evals = np.linalg.eigvalsh(D_num)
        tol = 1e-10
        evd: Dict[sp.Expr, int] = {}
        for ev in sorted(evals):
            found = False
            for key in evd:
                if abs(float(key) - ev) < tol:
                    evd[key] += 1
                    found = True
                    break
            if not found:
                evd[sp.Float(ev)] = 1
        return evd

    def dirac_spectrum_sorted(self) -> List[Tuple[sp.Expr, int]]:
        """Eigenvalues of D sorted by value, with multiplicities."""
        evd = self.dirac_spectrum()
        return sorted(evd.items(), key=lambda x: float(x[0]))

    # ------------------------------------------------------------------
    # pi-free certificate
    # ------------------------------------------------------------------

    def verify_pi_free(self) -> bool:
        """Certify that all entries of D, gamma, J are algebraic (pi-free).

        For uniform adjacency, all entries are exact Rationals/Integers.
        For CG-weighted adjacency, entries may contain products of
        Rationals with sqrt(Rational) from Wigner 3j symbols — these
        are algebraic numbers (not transcendental) and qualify as pi-free.

        The check verifies that no entry contains pi, log, exp, or any
        other transcendental function.
        """
        def _is_algebraic(expr: sp.Expr) -> bool:
            """Check if a sympy expression is algebraic (rational or sqrt of rational)."""
            if isinstance(expr, (
                sp.Rational,
                sp.Integer,
                sp.core.numbers.Zero,
                sp.core.numbers.One,
                sp.core.numbers.NegativeOne,
            )):
                return True
            # Allow Mul and Add of algebraic entries
            if isinstance(expr, sp.Mul):
                return all(_is_algebraic(arg) for arg in expr.args)
            if isinstance(expr, sp.Add):
                return all(_is_algebraic(arg) for arg in expr.args)
            # Allow Pow with rational base and rational exponent (e.g., sqrt)
            if isinstance(expr, sp.Pow):
                base, exp = expr.args
                return _is_algebraic(base) and isinstance(exp, sp.Rational)
            # Reject anything with pi, E, or transcendental functions
            if expr.has(sp.pi, sp.E, sp.log, sp.exp):
                return False
            # Check if it's a number that can be expressed as a radical
            if expr.is_number and expr.is_algebraic:
                return True
            return False

        for mat_name, mat in [
            ("D", self.dirac_operator),
            ("gamma", self.grading),
            ("J", self.real_structure),
        ]:
            for entry in mat:
                if not _is_algebraic(entry):
                    return False
        return True

    # ------------------------------------------------------------------
    # Supertrace
    # ------------------------------------------------------------------

    def supertrace(self, M: Matrix) -> sp.Expr:
        """Compute the supertrace Str(M) = Tr(gamma * M).

        The supertrace is the graded trace using the chirality gamma.
        """
        g = self.grading
        return (g * M).trace()

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    def summary(self) -> Dict:
        """Comprehensive summary of the spectral triple."""
        hs = self.hilbert_space_summary()

        spec = self.dirac_spectrum_sorted()
        spec_list = [(str(ev), mult) for ev, mult in spec]

        checks = {}
        checks["D_selfadjoint"] = self.check_selfadjoint()
        checks["gamma_squared_I"] = self.check_grading_square()
        checks["gamma_selfadjoint"] = self.check_grading_selfadjoint()

        D_gamma_result, D_gamma_msg = self.check_D_gamma_anticommute()
        checks["D_gamma_anticommute"] = D_gamma_result
        checks["D_gamma_anticommute_msg"] = D_gamma_msg
        checks["Lambda_gamma_commute"] = self.check_D_gamma_commute_diagonal()

        J2_result, J2_eps = self.check_J_squared()
        checks["J_squared"] = J2_result
        checks["J_squared_epsilon"] = J2_eps

        JD_result, JD_msg = self.check_J_D_relation()
        checks["J_D_relation"] = JD_result
        checks["J_D_relation_msg"] = JD_msg

        checks["J_gamma_commute"] = self.check_J_gamma_commute()
        checks["order_zero"] = self.check_order_zero()

        o1_result, o1_msg = self.check_order_one()
        checks["order_one"] = o1_result
        checks["order_one_msg"] = o1_msg

        checks["pi_free"] = self.verify_pi_free()

        return {
            "n_max": self.n_max,
            "kappa": str(self._kappa),
            "j_type": self._j_type,
            "adjacency_weights": self._adjacency_weights,
            "hilbert_space": hs,
            "spectrum": spec_list,
            "axiom_checks": checks,
        }

    def __repr__(self) -> str:
        parts = [
            f"FockSpectralTriple(n_max={self.n_max}",
            f"dim_H={self.dim_H}",
            f"n_sectors={self.n_sectors}",
            f"j_type={self._j_type!r}",
        ]
        if self._adjacency_weights != "uniform":
            parts.append(f"adjacency_weights={self._adjacency_weights!r}")
        return ", ".join(parts) + ")"
