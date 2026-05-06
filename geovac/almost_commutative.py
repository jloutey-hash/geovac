"""Almost-commutative extension of the GeoVac spectral triple — Sprint H1.

Sprint H1 (Higgs-from-inner-fluctuation) — IMPLEMENTATION.
This module replaces the Track 3 stub interface with the actual minimal
electroweak almost-commutative extension of the GeoVac spectral triple.

ARCHITECTURE (PI directive 2026-05-06)
======================================

Per the PI's Sprint H1 dispatch (post Track 2 closure, post §8 framing):

  - Use TRUTHFUL Camporesi-Higuchi Dirac D_GV
    (NOT offdiag CH; offdiag breaks JD = +DJ at residual ~2.0,
     see debug/real_structure_finite_nmax_memo.md §3D).
  - Internal couplings live on the M_n(C) = A_F factor itself:
    D = D_GV (x) 1_F + gamma_GV (x) D_F
    where gamma_GV is the chirality grading and D_F is a 4x4
    Hermitian operator on H_F = C^2_L (+) C^2_R, with the off-diagonal
    C ↔ H block carrying the Yukawa data.
  - J = J_GV (x) J_F with J_F the standard quaternionic charge
    conjugation J_F = i*sigma_2 K on the C^2 component representing H.
  - This is closer to canonical Chamseddine-Connes than the offdiag
    "Candidate A" originally floated in the scoping memo (which is
    ruled out by Track 2's clean negative on offdiag CH J-symmetry).

Mathematical setup
==================

A_F = C (+) H ----------------------------------------------------------

A_F = C ⊕ H is a real *-algebra of real dimension 5 (1 from C, 4 from H);
viewed as a complex *-algebra acting on H_F = C^4 it has C-dimension 3
(realizing q ∈ H as a 2x2 complex matrix via H ⊗_R C ≅ M_2(C)).

We label the elements by (lambda, q) with lambda ∈ C and q ∈ H realized
as a 2x2 complex matrix satisfying q = q1*1 + q2*i*sigma_1 + q3*i*sigma_2
+ q4*i*sigma_3 (i.e. q^* = q1*1 - q2*i*sigma_1 - q3*i*sigma_2 -
q4*i*sigma_3, the quaternionic involution).

H_F = matter ⊕ antimatter doubled --------------------------------------

The standard Connes-Chamseddine construction doubles the Hilbert space
into "matter" and "antimatter" sectors:

    H_F = H_F^matter (+) H_F^antimatter
        = C^4 (+) C^4
        = (nu_L, e_L, nu_R, e_R) (+) (nu_L, e_L, nu_R, e_R)^bar

with dim_C(H_F) = 8.

The algebra A_F = C ⊕ H acts on the MATTER sector only via

    pi_F^matter(lambda, q) = block_diag(q, diag(lambda, conj(lambda)))

extended to the full H_F by acting as 0 on the antimatter sector
(or equivalently, via the OPPOSITE algebra action on antimatter).

The faithful action on matter:
    pi_F(lambda, q) = block_diag(pi_F^matter, 0_4)

The 'opposite-algebra' action that J realizes:
    pi_F^opp(lambda, q) = block_diag(0_4, J_F pi_F^matter(lambda, q) J_F^{-1})

D_F: the finite Dirac --------------------------------------------------

D_F is an 8x8 Hermitian operator on H_F:

    D_F = [ M       0     ]
          [ 0   M^bar     ]

where M is the 4x4 matter-sector Yukawa-Dirac:

    M = [ 0       Y    ]
        [ Y^dag   0    ]   on (L, R)

with Y a 2x2 complex matrix Y: C^2_L -> C^2_R (Yukawa matrix). For
minimal one-generation lepton:

    Y = diag(y_nu, y_e)

The eigenvalues of D_F are +/-y_nu, +/-y_e (each doubled from
matter+antimatter sectors).

J_F: finite real structure --------------------------------------------

J_F swaps matter <-> antimatter with complex conjugation:

    J_F (psi_matter, psi_antimatter) = (conj(psi_antimatter), conj(psi_matter))

Equivalently, J_F = sigma_x (x) K on the matter/antimatter Z_2 grading,
where sigma_x is the swap matrix.

Verification: J_F^2 = sigma_x conj(sigma_x) = sigma_x sigma_x = +I.

So J_F^2 = +I, giving KO-dim of T_F in {0, 4, 6, 7}. The standard
Connes-Marcolli SM convention places T_F at KO-dim 6 (Connes-Marcolli
Table 13.1), with J_F D_F = +D_F J_F.

For our electroweak slice C ⊕ H this convention gives KO-dim 6
matching the Standard Model spectral triple's KO-dim.

Combined KO-dim
---------------

For T_combined = T_GV ⊗ T_F with KO(GV) = 3, KO(F) = 6:
    KO(combined) = 3 + 6 = 9 ≡ 1 (mod 8)
    epsilon_combined = epsilon_GV * epsilon_F = (-1)(+1) = -1
    epsilon'_combined = epsilon'_GV * epsilon'_F = (+1)(+1) = +1

KO-dim 1: (epsilon, epsilon') = (-, +) per the standard sign table.
So J^2 = -I and JD = +DJ.

This is consistent with the multiplication rules for tensor products of
real spectral triples (Connes-Marcolli 2008 Sec. 13.4).

CRUCIAL STRUCTURAL CONSEQUENCE OF MATTER/ANTIMATTER DOUBLING:

Because pi_F acts only on the matter sector (and J_F pi_F J_F^{-1} only
on the antimatter sector), order-zero [a, J b J^{-1}] = 0 holds
AUTOMATICALLY: a is supported on matter, J b J^{-1} is supported on
antimatter, they commute trivially.

This is the well-known reason CC's order-zero condition holds for the
finite SM triple, and it transfers to the GeoVac extension at the
finite-sector level.

Inner fluctuations
==================

For A = A_GV ⊗ A_F, omega = sum_i a_i [D, b_i] decomposes:

[D, b_i] = [D_GV ⊗ 1_F + gamma_GV ⊗ D_F, b_i^GV ⊗ b_i^F]
        = [D_GV, b_i^GV] ⊗ b_i^F
          + gamma_GV b_i^GV ⊗ [D_F, b_i^F]
          - b_i^GV gamma_GV ⊗ b_i^F D_F (subtracted)

Reorganizing using (gamma b)^GV = (-1)^|b| b gamma when gamma is a
Z_2 grading and b is graded:
For UNGRADED b (the case here: A_GV is bosonic / scalar, so
gamma_GV b_i^GV = b_i^GV gamma_GV) we get:

[D, b_i^GV ⊗ b_i^F]
    = [D_GV, b_i^GV] ⊗ b_i^F + b_i^GV gamma_GV ⊗ [D_F, b_i^F]

(because gamma_GV commutes with scalar multiplication).

Summing over i with multiplication on the left by a_i = a_i^GV ⊗ a_i^F:

omega = sum_i a_i^GV [D_GV, b_i^GV] ⊗ a_i^F b_i^F                  (gauge piece)
      + sum_i a_i^GV b_i^GV gamma_GV ⊗ a_i^F [D_F, b_i^F]            (Higgs piece)

The gauge piece is a 1-form in Omega^1(A_GV) tensored with an
algebra-element on A_F. Its trace over A_F gives the U(1) piece
(coefficient of 1_F) and its traceless part gives the SU(2) piece
(coefficient of pure-quaternion in H).

The Higgs piece has the off-diagonal C <-> H component of A_F
appearing (from [D_F, b_i^F] which is generically off-diagonal). This
is the Higgs scalar field Phi.

Higgs sector localization in matter/antimatter doubled H_F
----------------------------------------------------------

With H_F = matter (+) antimatter, A_F acts on matter sector only.
The matter-sector commutator [M_matter, b_F^matter] generates the
matter-sector Higgs scalar Phi_m: C^2_L^matter -> C^2_R^matter.
The antimatter sector picks up Phi_m^bar = conj(Phi_m) via
J_F-conjugation (the J fluctuation J omega J^{-1}).

So the Higgs scalar lives in the LR block of the MATTER-sector
4x4 sub-block of each 8x8 A_F slice. We extract it as
omega[:, :, 2:4, 0:2] within the matter sector (off-diagonal elements
of indices ranging over (0..3) on each side, then within the L<->R
sub-decomposition).

References
==========

Design memo: debug/almost_commutative_scoping_memo.md (2026-05-06)
Sprint H1 memo: debug/h1_ac_extension_memo.md (this sprint)

Track 2 (J at finite n_max): debug/real_structure_finite_nmax_memo.md

A. Connes & A. Chamseddine, "The spectral action principle," Comm. Math.
  Phys. 186 (1997) 731-750.
A. Chamseddine, A. Connes, M. Marcolli, "Gravity and the standard model
  with neutrino mixing," Adv. Theor. Math. Phys. 11 (2007) 991-1089.
A. Connes & M. Marcolli, *Noncommutative Geometry, Quantum Fields and
  Motives*, AMS Colloquium Publications 55, 2008. Especially Ch. 13.
W. D. van Suijlekom, *Noncommutative Geometry and Particle Physics*,
  Springer, 2015. Chapters 8, 9.
M. Marcolli & W. D. van Suijlekom, "Gauge networks in noncommutative
  geometry," J. Geom. Phys. 75 (2014), arXiv:1301.3480.
C. I. Perez-Sanchez, arXiv:2401.03705 (2024), arXiv:2508.17338 (2025).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Literal, Optional, Sequence, Tuple

import numpy as np

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
)
from geovac.real_structure import RealStructure, build_J_full_dirac
from geovac.spinor_operator_system import (
    SpinorLabel,
    build_spinor_multiplier_matrix,
    spinor_basis,
)


# ===========================================================================
# Pauli matrices (as building blocks for H sector)
# ===========================================================================

SIGMA_0 = np.eye(2, dtype=np.complex128)
SIGMA_1 = np.array([[0, 1], [1, 0]], dtype=np.complex128)
SIGMA_2 = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
SIGMA_3 = np.array([[1, 0], [0, -1]], dtype=np.complex128)


# ===========================================================================
# Finite electroweak data (A_F, H_F, D_F, J_F)
# ===========================================================================


def quaternion_to_matrix(q1: complex, q2: complex, q3: complex, q4: complex) -> np.ndarray:
    """Realize a quaternion q = q1 + q2*i + q3*j + q4*k as a 2x2 complex matrix.

    Standard embedding H -> M_2(C):
        1 -> sigma_0
        i -> i sigma_1
        j -> i sigma_2     (Note: depends on convention)
        k -> i sigma_3

    This makes H ⊗_R C ≅ M_2(C); the *-involution on H is realized as
    Hermitian adjoint on M_2(C).

    Parameters
    ----------
    q1, q2, q3, q4 : complex
        Quaternion components. For a real quaternion these are real.

    Returns
    -------
    2x2 complex array.
    """
    return q1 * SIGMA_0 + 1j * q2 * SIGMA_1 + 1j * q3 * SIGMA_2 + 1j * q4 * SIGMA_3


@dataclass(frozen=True)
class ElectroweakFiniteTriple:
    """Finite electroweak spectral triple (A_F, H_F, D_F, J_F).

    Following Connes-Marcolli 2008 Ch. 13, H_F is doubled into matter and
    antimatter sectors:

        H_F = H_F^matter (+) H_F^antimatter
            = C^4 (+) C^4

    Each sector has internal order (nu_L, e_L, nu_R, e_R).
    A_F = C ⊕ H acts on the matter sector only:
        pi_F(lambda, q) = block_diag(M_matter(lambda, q), 0_4)
        with M_matter(lambda, q) = block_diag(q, diag(lambda, conj(lambda))).

    D_F is an 8x8 Hermitian matrix:
        D_F = block_diag(M, M^bar)
    with M the matter-sector 4x4 Yukawa-Dirac (off-diag L<->R block).

    J_F swaps matter <-> antimatter with complex conjugation:
        J_F = sigma_x (x) I_4 followed by complex conjugation.
        J_F^2 = +I (KO-dim 6 of A_F).

    Attributes
    ----------
    yukawa_nu : complex
        The (nu_R, nu_L) Yukawa entry. Default 0.0 (no neutrino mass).
    yukawa_e : complex
        The (e_R, e_L) Yukawa entry. Default 0.0.
    """

    yukawa_nu: complex = 0.0 + 0.0j
    yukawa_e: complex = 0.0 + 0.0j

    @property
    def dim_H_F(self) -> int:
        """Dimension of H_F = 4 (matter) + 4 (antimatter) = 8."""
        return 8

    @property
    def dim_matter(self) -> int:
        """Dimension of one sector (matter or antimatter): 4."""
        return 4

    def matter_action(
        self, lam: complex, q_components: Tuple[complex, complex, complex, complex]
    ) -> np.ndarray:
        """4x4 representation of A_F on the matter sector ONLY.

        Block-diagonal: L-block = q, R-block = diag(lam, conj(lam)).
        """
        q_mat = quaternion_to_matrix(*q_components)
        M = np.zeros((4, 4), dtype=np.complex128)
        M[0:2, 0:2] = q_mat
        M[2, 2] = lam
        M[3, 3] = np.conj(lam)
        return M

    def algebra_action(
        self, lam: complex, q_components: Tuple[complex, complex, complex, complex]
    ) -> np.ndarray:
        """Faithful *-representation pi_F: A_F -> B(H_F) at element (lambda, q).

        pi_F acts on the MATTER sector only (with antimatter sector zero).
        The opposite-algebra action on antimatter is realized via
        J pi_F J^{-1}.

        Returns
        -------
        8x8 complex array (block_diag(M_matter, 0_4)).
        """
        M_matter = self.matter_action(lam, q_components)
        action = np.zeros((8, 8), dtype=np.complex128)
        action[0:4, 0:4] = M_matter
        # antimatter sector: zero (algebra acts only on matter)
        return action

    def matter_dirac(self) -> np.ndarray:
        """4x4 matter-sector Dirac M with off-diagonal L<->R block.

            M = [ 0      Y^dag ]
                [ Y      0     ]
            with Y = diag(y_nu, y_e).
        """
        M = np.zeros((4, 4), dtype=np.complex128)
        Y = np.diag([self.yukawa_nu, self.yukawa_e]).astype(np.complex128)
        M[0:2, 2:4] = Y.conj().T
        M[2:4, 0:2] = Y
        return M

    def dirac_F(self) -> np.ndarray:
        """Build the 8x8 finite Dirac D_F.

            D_F = block_diag(M, M^bar) on (matter, antimatter)

        where M^bar = conj(M) on antimatter sector.

        Hermitian: D_F^dag = D_F.
        """
        M = self.matter_dirac()
        D = np.zeros((8, 8), dtype=np.complex128)
        D[0:4, 0:4] = M
        D[4:8, 4:8] = np.conj(M)  # antimatter Dirac is the conjugate
        # Verify Hermitian
        assert np.allclose(D, D.conj().T, atol=1e-14), "D_F not Hermitian"
        return D

    def real_structure_F(self) -> np.ndarray:
        """Finite real structure J_F on H_F.

        J_F swaps matter <-> antimatter with complex conjugation:

            J_F (psi_matter, psi_antimatter)
                = (conj(psi_antimatter), conj(psi_matter))

        As J_F = U_F K with U_F a unitary acting on the matter/antimatter
        Z_2 grading:

            U_F = sigma_x (x) I_4 = [[0_4, I_4], [I_4, 0_4]]

        Verification: J_F^2 = U_F conj(U_F) = sigma_x sigma_x (x) I_4 = +I_8.
        So J_F^2 = +I.
        """
        I_4 = np.eye(4, dtype=np.complex128)
        U_F = np.zeros((8, 8), dtype=np.complex128)
        U_F[0:4, 4:8] = I_4
        U_F[4:8, 0:4] = I_4
        return U_F


# ===========================================================================
# Combined almost-commutative spectral triple (A_GV ⊗ A_F, H_GV ⊗ H_F, D)
# ===========================================================================


CandidateD_F_Source = Literal[
    "candidate_a_chirality",   # off-diagonal CH chirality bridging (RULED OUT by Track 2)
    "candidate_b_sturmian",    # DUCC / Sturmian core-valence bridging (speculative)
    "imposed",                 # CC default: D_F off-diagonal entries imposed by hand
    "geovac_natural",          # Test: see whether GeoVac data forces Y = 0
]


class AlmostCommutativeTriple:
    """Almost-commutative extension T = T_GV ⊗ T_F at finite n_max.

    H = H_GV ⊗ H_F, A = A_GV ⊗ A_F, D = D_GV ⊗ 1_F + gamma_GV ⊗ D_F.

    Parameters
    ----------
    n_max : int
        Fock cutoff of the GeoVac sector.
    finite : ElectroweakFiniteTriple, optional
        Finite electroweak triple. Default zero Yukawa (no Higgs vev).
    candidate_source : CandidateD_F_Source
        How D_F is selected. See module docstring §3.3.
    """

    def __init__(
        self,
        n_max: int,
        finite: Optional[ElectroweakFiniteTriple] = None,
        candidate_source: CandidateD_F_Source = "imposed",
    ) -> None:
        if n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {n_max}")
        self.n_max = n_max
        self.finite = finite or ElectroweakFiniteTriple()
        self.candidate_source = candidate_source

        # Build the GeoVac sector data.
        self._op_sys_GV = FullDiracTruncatedOperatorSystem(n_max)
        self._basis_GV = self._op_sys_GV.basis
        self._dim_GV = self._op_sys_GV.dim_H

        # Truthful CH Dirac (per Track 2 + PI directive: NOT offdiag).
        self._D_GV = camporesi_higuchi_full_dirac_matrix(self._basis_GV)

        # Chirality grading gamma_GV: diagonal +/-1 according to chirality
        # label of basis. This satisfies gamma_GV^2 = I and {gamma_GV, D_GV} = 0
        # ONLY IF D_GV is off-diagonal in chirality. Truthful CH is DIAGONAL
        # in chirality with eigenvalue chi*(n+1/2), so gamma_GV commutes with
        # D_GV — gamma_GV is the SIGN of D_GV here, not an anticommuting
        # grading.
        self._gamma_GV = np.diag(
            [float(b.chirality) for b in self._basis_GV]
        ).astype(np.complex128)

        # Finite electroweak data.
        self._D_F = self.finite.dirac_F()
        self._J_F_U = self.finite.real_structure_F()

        # GeoVac real structure.
        self._J_GV = build_J_full_dirac(n_max)

        # Combined dimensions.
        self._dim_H = self._dim_GV * self.finite.dim_H_F  # = dim_GV * 4

    # ------------------------------------------------------------------
    # Hilbert space and operators
    # ------------------------------------------------------------------

    @property
    def dim_GV(self) -> int:
        return self._dim_GV

    @property
    def dim_F(self) -> int:
        return self.finite.dim_H_F

    @property
    def dim_H(self) -> int:
        return self._dim_H

    def dirac_combined(self) -> np.ndarray:
        """D = D_GV ⊗ 1_F + gamma_GV ⊗ D_F."""
        I_F = np.eye(self.finite.dim_H_F, dtype=np.complex128)
        D = np.kron(self._D_GV, I_F) + np.kron(self._gamma_GV, self._D_F)
        return D

    def D_GV(self) -> np.ndarray:
        return self._D_GV.copy()

    def D_F(self) -> np.ndarray:
        return self._D_F.copy()

    def gamma_GV(self) -> np.ndarray:
        return self._gamma_GV.copy()

    def real_structure_combined(self) -> np.ndarray:
        """J = J_GV ⊗ J_F as a unitary U so that J(psi) = U @ conj(psi).

        Combined J^2 = (-1)*(-1) = +1 (KO-dim shifts 3+2 = 5).
        """
        return np.kron(self._J_GV.U, self._J_F_U)

    def algebra_element(
        self, M_GV: np.ndarray, lam: complex, q_components: Tuple[complex, complex, complex, complex],
    ) -> np.ndarray:
        """Build pi(a) = pi_GV(M_GV) ⊗ pi_F(lam, q).

        Parameters
        ----------
        M_GV : ndarray of shape (dim_GV, dim_GV)
            A multiplier matrix from the truncated operator system A_GV.
        lam : complex
            C-summand element.
        q_components : (q1, q2, q3, q4)
            Quaternion components.
        """
        a_F = self.finite.algebra_action(lam, q_components)
        return np.kron(M_GV, a_F)

    def gv_multiplier(self, k: int) -> np.ndarray:
        """Get the k-th multiplier matrix from the GV operator system."""
        return self._op_sys_GV.multiplier_matrices[k]

    @property
    def n_gv_multipliers(self) -> int:
        return len(self._op_sys_GV.multiplier_matrices)

    # ------------------------------------------------------------------
    # Inner fluctuation
    # ------------------------------------------------------------------

    def inner_fluctuation_one_form(
        self,
        generators: Sequence[Tuple[np.ndarray, complex, Tuple[complex, complex, complex, complex],
                                  np.ndarray, complex, Tuple[complex, complex, complex, complex]]],
    ) -> np.ndarray:
        """Build omega = sum_i a_i [D, b_i] for a list of (a, b) pairs.

        Each entry is (M_GV_a, lam_a, q_a, M_GV_b, lam_b, q_b) specifying
        a_i = M_GV_a ⊗ pi_F(lam_a, q_a), b_i = M_GV_b ⊗ pi_F(lam_b, q_b).

        Returns
        -------
        omega : ndarray of shape (dim_H, dim_H), generically not Hermitian.
        """
        D = self.dirac_combined()
        omega = np.zeros((self._dim_H, self._dim_H), dtype=np.complex128)
        for (M_a, lam_a, q_a, M_b, lam_b, q_b) in generators:
            a = self.algebra_element(M_a, lam_a, q_a)
            b = self.algebra_element(M_b, lam_b, q_b)
            omega = omega + a @ (D @ b - b @ D)
        return omega

    def fluctuated_dirac(
        self,
        omega: np.ndarray,
        epsilon_prime: int = -1,
    ) -> np.ndarray:
        """D_omega = D + omega + epsilon' * J omega J^{-1}.

        For combined KO-dim 5, epsilon' = -1.
        For an antilinear J = U K, J op J^{-1} = U conj(op) U^T.
        """
        D = self.dirac_combined()
        U = self.real_structure_combined()
        J_omega_Jinv = U @ np.conj(omega) @ U.T
        return D + omega + epsilon_prime * J_omega_Jinv

    # ------------------------------------------------------------------
    # Fluctuation classification: extract gauge / Higgs sectors
    # ------------------------------------------------------------------

    def decompose_fluctuation(
        self,
        omega: np.ndarray,
    ) -> dict:
        """Decompose omega into matter, antimatter, gauge, and Higgs blocks.

        omega has shape (dim_GV * 8, dim_GV * 8). Reshape to
        (dim_GV, 8, dim_GV, 8) -> (dim_GV, dim_GV, 8, 8), and within each
        8x8 fiber block at GV index pair (i, j), decompose:

            8x8 block = [ 4x4 matter |  4x4 mat-antimat ]
                        [ 4x4 antimat-mat | 4x4 antimat ]

        Within the matter (and antimatter) 4x4 sub-block, decompose into
        L (rows/cols 0:2) and R (rows/cols 2:4):

            matter 4x4 = [ 2x2 LL (gauge) | 2x2 LR (Higgs Phi^dag) ]
                         [ 2x2 RL (Higgs Phi) | 2x2 RR (gauge) ]

        Returns:
          - 'higgs_matter_LR': (d, d, 2, 2) — Phi^dag in matter
          - 'higgs_matter_RL': (d, d, 2, 2) — Phi in matter
          - 'higgs_antimatter_LR': (d, d, 2, 2)
          - 'higgs_antimatter_RL': (d, d, 2, 2)
          - 'gauge_matter_L': (d, d, 2, 2) — quaternion content of L-block in matter
          - 'gauge_matter_R': (d, d, 2, 2) — diag(lam, conj lam) content of R-block
          - 'gauge_antimatter_L', 'gauge_antimatter_R': same on antimatter sector
          - 'matter_antimatter_off': (d, d, 4, 4) — generally zero for legitimate
                                    inner fluctuations (matter and antimatter
                                    sectors decouple).
        """
        d = self._dim_GV
        # Reshape omega: (d * 8, d * 8) -> (d, 8, d, 8) -> (d, d, 8, 8)
        om_4d = omega.reshape(d, 8, d, 8).transpose(0, 2, 1, 3)

        # Matter (rows/cols 0:4) and antimatter (rows/cols 4:8) sub-blocks
        matter_block = om_4d[:, :, 0:4, 0:4]              # (d, d, 4, 4)
        antimatter_block = om_4d[:, :, 4:8, 4:8]          # (d, d, 4, 4)
        ma_off = om_4d[:, :, 0:4, 4:8]                    # matter <- antimatter (4x4)
        am_off = om_4d[:, :, 4:8, 0:4]                    # antimatter <- matter (4x4)

        # Within matter_block, decompose by L (0:2) / R (2:4)
        higgs_matter_LR = matter_block[:, :, 0:2, 2:4]    # L-from-R, contains Phi^dag content
        higgs_matter_RL = matter_block[:, :, 2:4, 0:2]    # R-from-L, contains Phi content
        gauge_matter_L = matter_block[:, :, 0:2, 0:2]     # L-block (quaternion sector)
        gauge_matter_R = matter_block[:, :, 2:4, 2:4]     # R-block (lambda sector)

        # Same for antimatter
        higgs_antimatter_LR = antimatter_block[:, :, 0:2, 2:4]
        higgs_antimatter_RL = antimatter_block[:, :, 2:4, 0:2]
        gauge_antimatter_L = antimatter_block[:, :, 0:2, 0:2]
        gauge_antimatter_R = antimatter_block[:, :, 2:4, 2:4]

        return {
            "higgs_matter_LR": higgs_matter_LR,
            "higgs_matter_RL": higgs_matter_RL,
            "higgs_antimatter_LR": higgs_antimatter_LR,
            "higgs_antimatter_RL": higgs_antimatter_RL,
            "gauge_matter_L": gauge_matter_L,
            "gauge_matter_R": gauge_matter_R,
            "gauge_antimatter_L": gauge_antimatter_L,
            "gauge_antimatter_R": gauge_antimatter_R,
            "matter_antimatter_off": ma_off,
            "antimatter_matter_off": am_off,
        }

    def extract_higgs(self, omega: np.ndarray) -> np.ndarray:
        """Extract the matter-sector Higgs scalar field Phi as (dim_GV, dim_GV, 2, 2).

        Phi is the (R, L) block of omega in the matter sector. A non-zero
        Phi signals a non-trivial Higgs sector from the inner-fluctuation
        decomposition.
        """
        return self.decompose_fluctuation(omega)["higgs_matter_RL"]

    def higgs_norm(self, omega: np.ndarray) -> float:
        """Return the Frobenius norm of the matter-sector Higgs scalar block.

        We use the matter-sector L<->R off-diagonal as the canonical Higgs
        signature; non-zero matter higgs is the relevant test for the
        falsifier.
        """
        decomp = self.decompose_fluctuation(omega)
        # Combine matter LR + RL and antimatter LR + RL (Higgs is the full
        # off-L<->R content)
        higgs_total = (
            np.linalg.norm(decomp["higgs_matter_LR"]) ** 2
            + np.linalg.norm(decomp["higgs_matter_RL"]) ** 2
            + np.linalg.norm(decomp["higgs_antimatter_LR"]) ** 2
            + np.linalg.norm(decomp["higgs_antimatter_RL"]) ** 2
        )
        return float(np.sqrt(higgs_total))

    def gauge_norm(self, omega: np.ndarray) -> float:
        """Return the Frobenius norm of the gauge piece (L-block + R-block).

        Combines matter and antimatter gauge sectors.
        """
        decomp = self.decompose_fluctuation(omega)
        return float(
            np.sqrt(
                np.linalg.norm(decomp["gauge_matter_L"]) ** 2
                + np.linalg.norm(decomp["gauge_matter_R"]) ** 2
                + np.linalg.norm(decomp["gauge_antimatter_L"]) ** 2
                + np.linalg.norm(decomp["gauge_antimatter_R"]) ** 2
            )
        )

    def matter_antimatter_off_norm(self, omega: np.ndarray) -> float:
        """Frobenius norm of the matter<->antimatter off-block of omega.

        For inner fluctuations on the standard CC almost-commutative triple,
        this should be ZERO: the algebra acts on matter only, and so does
        [D, .] (D is block-diagonal in matter/antimatter); J flips sectors.
        """
        decomp = self.decompose_fluctuation(omega)
        return float(
            np.sqrt(
                np.linalg.norm(decomp["matter_antimatter_off"]) ** 2
                + np.linalg.norm(decomp["antimatter_matter_off"]) ** 2
            )
        )

    # ------------------------------------------------------------------
    # Falsifier check (memo §5)
    # ------------------------------------------------------------------

    def check_natural_negative(
        self,
        n_random_generators: int = 50,
        seed: int = 42,
    ) -> Tuple[bool, str, dict]:
        """Test whether candidate D_F sources produce only gauge fluctuations.

        Strategy: build a random sample of (a, b) pairs and check whether
        omega = sum a_i [D, b_i] has any non-trivial Higgs-block content.

        Returns (negative_holds, reason, data) where:
          negative_holds = True if Higgs sector identically zero
          reason = explanation string
          data = dict with sampled norms

        Distinguishes three regimes:
          (a) D_F = 0 (zero Yukawa): higgs_norm should be 0 — sanity check.
          (b) D_F off-diagonal nonzero: higgs_norm > 0 generically.

        For the GEOVAC-NATURAL falsifier test:
          We ask whether GeoVac structure FORCES D_F off-diagonal = 0
          (closing G2 negatively) or admits non-trivial D_F.

        Currently neither candidate A nor candidate B selects a D_F.
        With candidate_source = "imposed", D_F is whatever finite says.
        """
        rng = np.random.default_rng(seed)
        n_mult = self.n_gv_multipliers

        higgs_norms = []
        gauge_norms = []
        for _ in range(n_random_generators):
            i_a = rng.integers(0, n_mult)
            i_b = rng.integers(0, n_mult)
            # Random complex lambda and quaternion components (small magnitude)
            lam_a = (rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
            lam_b = (rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
            q_a = tuple((rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
                        for _ in range(4))
            q_b = tuple((rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
                        for _ in range(4))

            generators = [(
                self.gv_multiplier(i_a), lam_a, q_a,
                self.gv_multiplier(i_b), lam_b, q_b,
            )]
            omega = self.inner_fluctuation_one_form(generators)
            higgs_norms.append(self.higgs_norm(omega))
            gauge_norms.append(self.gauge_norm(omega))

        higgs_max = max(higgs_norms)
        gauge_max = max(gauge_norms)
        higgs_mean = float(np.mean(higgs_norms))
        gauge_mean = float(np.mean(gauge_norms))

        # The natural negative holds if Higgs is identically zero.
        epsilon = 1e-12
        higgs_zero = higgs_max < epsilon * max(gauge_max, 1.0)

        if higgs_zero:
            reason = (
                "Higgs sector identically zero across all sampled generators "
                f"(higgs_max = {higgs_max:.3e}, gauge_max = {gauge_max:.3e}); "
                "D_F off-diagonal does not produce a Higgs scalar."
            )
        else:
            ratio = higgs_max / gauge_max if gauge_max > 0 else float('inf')
            reason = (
                f"Higgs sector non-trivial (higgs_max = {higgs_max:.3e}, "
                f"gauge_max = {gauge_max:.3e}, higgs/gauge = {ratio:.3f})"
            )

        data = {
            "higgs_norms": higgs_norms,
            "gauge_norms": gauge_norms,
            "higgs_max": higgs_max,
            "gauge_max": gauge_max,
            "higgs_mean": higgs_mean,
            "gauge_mean": gauge_mean,
            "n_samples": n_random_generators,
        }

        return higgs_zero, reason, data


# ===========================================================================
# Convenience constructors
# ===========================================================================


def minimal_electroweak_triple(
    n_max: int,
    yukawa_e: complex = 0.0,
    yukawa_nu: complex = 0.0,
    candidate_source: CandidateD_F_Source = "imposed",
) -> AlmostCommutativeTriple:
    """Return the minimal electroweak almost-commutative extension at given n_max.

    Default Yukawa = 0 (sanity-check construction; reproduces gauge sector
    only). Setting yukawa_e or yukawa_nu nonzero switches on the Higgs sector.
    """
    finite = ElectroweakFiniteTriple(yukawa_nu=yukawa_nu, yukawa_e=yukawa_e)
    return AlmostCommutativeTriple(
        n_max=n_max,
        finite=finite,
        candidate_source=candidate_source,
    )


def gauge_only_triple(n_max: int) -> AlmostCommutativeTriple:
    """Return the gauge-only sanity triple: D_F = 0, no Higgs."""
    return minimal_electroweak_triple(n_max, yukawa_e=0.0, yukawa_nu=0.0)


# ===========================================================================
# Spectral action (leading-order expansion)
# ===========================================================================


def heat_kernel_trace(D: np.ndarray, t: float) -> float:
    """Tr exp(-t D^2). Leading-order spectral action input.

    For a finite Hermitian operator D, this is sum_i exp(-t * lambda_i^2)
    where lambda_i are the eigenvalues.
    """
    eigs = np.linalg.eigvalsh(D)
    return float(np.sum(np.exp(-t * eigs ** 2)))


def spectral_action_eigenvalues(D: np.ndarray) -> np.ndarray:
    """Diagonalize D and return its eigenvalues."""
    return np.linalg.eigvalsh(D)


def yukawa_dependence_scan(
    n_max: int,
    yukawa_grid: Sequence[float],
) -> dict:
    """Compute the Higgs-mass-related quantity f(D_omega) over a Yukawa grid.

    For each y in yukawa_grid, build T(y) with yukawa_e = y, yukawa_nu = 0,
    compute the eigenvalue spectrum of D_combined, and report
    Tr(D^2) and ||D - D_0||_F as functions of y.
    """
    Tr_D2 = []
    Frob_shift = []
    eigs_list = []
    for y in yukawa_grid:
        T = minimal_electroweak_triple(n_max, yukawa_e=y, yukawa_nu=0.0)
        D = T.dirac_combined()
        Tr_D2.append(float(np.real(np.trace(D @ D))))
        if abs(y) < 1e-15:
            Frob_shift.append(0.0)
        else:
            T0 = gauge_only_triple(n_max)
            D0 = T0.dirac_combined()
            Frob_shift.append(float(np.linalg.norm(D - D0)))
        eigs_list.append(spectral_action_eigenvalues(D))

    return {
        "yukawa_grid": list(yukawa_grid),
        "Tr_D_squared": Tr_D2,
        "Frob_shift": Frob_shift,
        "eigenvalues_per_yukawa": [list(e) for e in eigs_list],
    }
