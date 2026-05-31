"""Full Connes Standard Model on the GeoVac spectral triple — Sprint G4a.

Extends Sprint H1's electroweak slice (A_F = C (+) H, dim H_F = 8) to the
full Connes-Chamseddine-Marcolli finite algebra:

    A_F = C (+) H (+) M_3(C)

with dim H_F = 32 (16 matter + 16 antimatter).

MATTER HILBERT SPACE (16D, one generation)
==========================================

    H_F^matter = H_lepton (+) H_quark

Leptons (4D):  (nu_L, e_L, nu_R, e_R)
Quarks (12D):  (u_L, d_L, u_R, d_R) (x) (r, g, b)

Ordering: leptons first (indices 0-3), then quarks in flavor (x) color
layout (indices 4-15) with color running faster:
    4-6:   (u_L^r, u_L^g, u_L^b)
    7-9:   (d_L^r, d_L^g, d_L^b)
    10-12: (u_R^r, u_R^g, u_R^b)
    13-15: (d_R^r, d_R^g, d_R^b)

Full H_F = matter (+) antimatter = C^32.

ALGEBRA ACTION
==============

pi(lam, q, m) on matter sector:
    Leptons:  block_diag(q, diag(lam, conj(lam)))               [4x4]
    Quarks:   block_diag(q, diag(lam, conj(lam))) (x) m         [12x12]

M_3(C) acts trivially on leptons (color singlet) and as m on quark
colors. This is the standard Connes-Marcolli representation.

GAUGE GROUP
===========

Inner automorphisms of A_F = C (+) H (+) M_3(C) yield:
    U(1)_Y x SU(2)_L x SU(3)_c

The inner fluctuation omega = sum a_i [D, b_i] decomposes into:
    - Gauge: U(1) x SU(2) x SU(3) connection 1-forms
    - Higgs: off-diagonal L<->R blocks (lepton + quark Yukawa)

PREDICTED VERDICT: POSITIVE-THIN
=================================

Construction works, SM gauge group recovered. BUT: Yukawa, hypercharge
assignments, chirality assignment, and generation count are imposed
inputs (Paper 18 calibration tier), not GeoVac-derived. This is the
G2-G3 corollary generalized: the framework is SM-consistent, not
SM-selecting.

References
----------
A. Chamseddine, A. Connes, M. Marcolli, Adv. Theor. Math. Phys. 11
    (2007) 991-1089.
W. D. van Suijlekom, Noncommutative Geometry and Particle Physics,
    Springer, 2015, Ch. 11.
A. Connes & M. Marcolli, NCG, QF and Motives, AMS, 2008, Ch. 13.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np

from geovac.almost_commutative import (
    AlmostCommutativeTriple,
    SIGMA_0,
    SIGMA_1,
    SIGMA_2,
    SIGMA_3,
    quaternion_to_matrix,
)
from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
)
from geovac.real_structure import build_J_full_dirac


# ======================================================================
# Gell-Mann matrices (SU(3) generators)
# ======================================================================

def gell_mann_matrices() -> list[np.ndarray]:
    """The eight Gell-Mann matrices lambda_1 through lambda_8."""
    lam = [np.zeros((3, 3), dtype=np.complex128) for _ in range(8)]

    # lambda_1
    lam[0][0, 1] = 1; lam[0][1, 0] = 1
    # lambda_2
    lam[1][0, 1] = -1j; lam[1][1, 0] = 1j
    # lambda_3
    lam[2][0, 0] = 1; lam[2][1, 1] = -1
    # lambda_4
    lam[3][0, 2] = 1; lam[3][2, 0] = 1
    # lambda_5
    lam[4][0, 2] = -1j; lam[4][2, 0] = 1j
    # lambda_6
    lam[5][1, 2] = 1; lam[5][2, 1] = 1
    # lambda_7
    lam[6][1, 2] = -1j; lam[6][2, 1] = 1j
    # lambda_8
    lam[7][0, 0] = 1 / np.sqrt(3)
    lam[7][1, 1] = 1 / np.sqrt(3)
    lam[7][2, 2] = -2 / np.sqrt(3)

    return lam


GELL_MANN = gell_mann_matrices()


# ======================================================================
# Standard Model finite triple (A_F, H_F, D_F, J_F, gamma_F)
# ======================================================================

@dataclass(frozen=True)
class StandardModelFiniteTriple:
    """Finite Standard Model spectral triple at KO-dim 6.

    A_F = C (+) H (+) M_3(C) acting on H_F = C^32 (matter + antimatter).

    Attributes
    ----------
    yukawa_nu, yukawa_e : complex
        Lepton Yukawa couplings.
    yukawa_u, yukawa_d : complex
        Quark Yukawa couplings (one generation).
    """

    yukawa_nu: complex = 0.0 + 0.0j
    yukawa_e: complex = 0.0 + 0.0j
    yukawa_u: complex = 0.0 + 0.0j
    yukawa_d: complex = 0.0 + 0.0j

    @property
    def dim_matter(self) -> int:
        return 16

    @property
    def dim_H_F(self) -> int:
        return 32

    def _electroweak_block(
        self, lam: complex, q_components: Tuple[complex, complex, complex, complex]
    ) -> np.ndarray:
        """4x4 electroweak action: block_diag(q, diag(lam, conj(lam)))."""
        q_mat = quaternion_to_matrix(*q_components)
        M = np.zeros((4, 4), dtype=np.complex128)
        M[0:2, 0:2] = q_mat
        M[2, 2] = lam
        M[3, 3] = np.conj(lam)
        return M

    def matter_action(
        self,
        lam: complex,
        q_components: Tuple[complex, complex, complex, complex],
        m: np.ndarray,
    ) -> np.ndarray:
        """16x16 representation of A_F on the matter sector.

        Leptons (0:4):   block_diag(q, diag(lam, conj(lam)))
        Quarks (4:16):   block_diag(q, diag(lam, conj(lam))) (x) m
        """
        ew = self._electroweak_block(lam, q_components)
        M = np.zeros((16, 16), dtype=np.complex128)
        M[0:4, 0:4] = ew
        M[4:16, 4:16] = np.kron(ew, m)
        return M

    def algebra_action(
        self,
        lam: complex,
        q_components: Tuple[complex, complex, complex, complex],
        m: np.ndarray,
    ) -> np.ndarray:
        """32x32 representation pi_F: A_F -> B(H_F).

        Acts on matter sector only; antimatter is via J_F pi J_F^{-1}.
        """
        M_matter = self.matter_action(lam, q_components, m)
        action = np.zeros((32, 32), dtype=np.complex128)
        action[0:16, 0:16] = M_matter
        return action

    def _yukawa_block_lepton(self) -> np.ndarray:
        """4x4 lepton-sector Yukawa-Dirac."""
        M = np.zeros((4, 4), dtype=np.complex128)
        Y = np.diag([self.yukawa_nu, self.yukawa_e]).astype(np.complex128)
        M[0:2, 2:4] = Y.conj().T
        M[2:4, 0:2] = Y
        return M

    def _yukawa_block_quark(self) -> np.ndarray:
        """12x12 quark-sector Yukawa-Dirac (diagonal in color)."""
        M_flavor = np.zeros((4, 4), dtype=np.complex128)
        Y = np.diag([self.yukawa_u, self.yukawa_d]).astype(np.complex128)
        M_flavor[0:2, 2:4] = Y.conj().T
        M_flavor[2:4, 0:2] = Y
        return np.kron(M_flavor, np.eye(3, dtype=np.complex128))

    def dirac_F(self) -> np.ndarray:
        """32x32 finite Dirac D_F = block_diag(M_matter, conj(M_matter))."""
        M = np.zeros((16, 16), dtype=np.complex128)
        M[0:4, 0:4] = self._yukawa_block_lepton()
        M[4:16, 4:16] = self._yukawa_block_quark()
        D = np.zeros((32, 32), dtype=np.complex128)
        D[0:16, 0:16] = M
        D[16:32, 16:32] = np.conj(M)
        assert np.allclose(D, D.conj().T, atol=1e-14), "D_F not Hermitian"
        return D

    def real_structure_F(self) -> np.ndarray:
        """32x32 real structure J_F (unitary part, J = U K).

        Swaps matter <-> antimatter: U_F = [[0, I_16], [I_16, 0]].
        J_F^2 = +I (KO-dim 6).
        """
        I_16 = np.eye(16, dtype=np.complex128)
        U = np.zeros((32, 32), dtype=np.complex128)
        U[0:16, 16:32] = I_16
        U[16:32, 0:16] = I_16
        return U

    def chirality_F(self) -> np.ndarray:
        """32x32 chirality gamma_F (KO-dim 6 convention).

        Matter:      L = +1, R = -1 (both leptons and quarks)
        Antimatter:  sign-flipped (to maintain {J_F, gamma_F} = 0)
        """
        g_lepton = np.diag([1.0, 1.0, -1.0, -1.0])
        g_quark = np.kron(
            np.diag([1.0, 1.0, -1.0, -1.0]),
            np.eye(3)
        )
        g_matter = np.zeros((16, 16), dtype=np.float64)
        g_matter[0:4, 0:4] = g_lepton
        g_matter[4:16, 4:16] = g_quark

        G = np.zeros((32, 32), dtype=np.complex128)
        G[0:16, 0:16] = g_matter
        G[16:32, 16:32] = -g_matter
        return G


# ======================================================================
# Full SM almost-commutative triple
# ======================================================================

class StandardModelACTriple:
    """Almost-commutative extension T = T_GV (x) T_F^SM at finite n_max.

    H = H_GV (x) H_F, A = A_GV (x) A_F, D = D_GV (x) 1_F + gamma_GV (x) D_F.

    With A_F = C (+) H (+) M_3(C) and dim H_F = 32.
    """

    def __init__(
        self,
        n_max: int,
        finite: Optional[StandardModelFiniteTriple] = None,
    ) -> None:
        if n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {n_max}")
        self.n_max = n_max
        self.finite = finite or StandardModelFiniteTriple()

        self._op_sys_GV = FullDiracTruncatedOperatorSystem(n_max)
        self._basis_GV = self._op_sys_GV.basis
        self._dim_GV = self._op_sys_GV.dim_H

        self._D_GV = camporesi_higuchi_full_dirac_matrix(self._basis_GV)
        self._gamma_GV = np.diag(
            [float(b.chirality) for b in self._basis_GV]
        ).astype(np.complex128)

        self._D_F = self.finite.dirac_F()
        self._J_F_U = self.finite.real_structure_F()
        self._gamma_F = self.finite.chirality_F()

        self._J_GV = build_J_full_dirac(n_max)

        self._dim_F = self.finite.dim_H_F
        self._dim_H = self._dim_GV * self._dim_F

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def dim_GV(self) -> int:
        return self._dim_GV

    @property
    def dim_F(self) -> int:
        return self._dim_F

    @property
    def dim_H(self) -> int:
        return self._dim_H

    def D_GV(self) -> np.ndarray:
        return self._D_GV.copy()

    def D_F(self) -> np.ndarray:
        return self._D_F.copy()

    def gamma_GV(self) -> np.ndarray:
        return self._gamma_GV.copy()

    def gamma_F(self) -> np.ndarray:
        return self._gamma_F.copy()

    # ------------------------------------------------------------------
    # Combined operators
    # ------------------------------------------------------------------

    def dirac_combined(self) -> np.ndarray:
        """D = D_GV (x) 1_F + gamma_GV (x) D_F."""
        I_F = np.eye(self._dim_F, dtype=np.complex128)
        return np.kron(self._D_GV, I_F) + np.kron(self._gamma_GV, self._D_F)

    def real_structure_combined(self) -> np.ndarray:
        """J = J_GV (x) J_F (unitary part)."""
        return np.kron(self._J_GV.U, self._J_F_U)

    def chirality_combined(self) -> np.ndarray:
        """gamma = gamma_GV (x) gamma_F."""
        return np.kron(self._gamma_GV, self._gamma_F)

    # ------------------------------------------------------------------
    # Algebra elements
    # ------------------------------------------------------------------

    def algebra_element(
        self,
        M_GV: np.ndarray,
        lam: complex,
        q_components: Tuple[complex, complex, complex, complex],
        m: np.ndarray,
    ) -> np.ndarray:
        """Build pi(a) = pi_GV(M_GV) (x) pi_F(lam, q, m)."""
        a_F = self.finite.algebra_action(lam, q_components, m)
        return np.kron(M_GV, a_F)

    def gv_multiplier(self, k: int) -> np.ndarray:
        return self._op_sys_GV.multiplier_matrices[k]

    @property
    def n_gv_multipliers(self) -> int:
        return len(self._op_sys_GV.multiplier_matrices)

    # ------------------------------------------------------------------
    # Inner fluctuation
    # ------------------------------------------------------------------

    def inner_fluctuation_one_form(
        self,
        generators: Sequence[Tuple[
            np.ndarray, complex, Tuple[complex, ...], np.ndarray,
            np.ndarray, complex, Tuple[complex, ...], np.ndarray,
        ]],
    ) -> np.ndarray:
        """Build omega = sum_i a_i [D, b_i].

        Each entry: (M_GV_a, lam_a, q_a, m_a, M_GV_b, lam_b, q_b, m_b).
        """
        D = self.dirac_combined()
        omega = np.zeros((self._dim_H, self._dim_H), dtype=np.complex128)
        for (M_a, lam_a, q_a, m_a, M_b, lam_b, q_b, m_b) in generators:
            a = self.algebra_element(M_a, lam_a, q_a, m_a)
            b = self.algebra_element(M_b, lam_b, q_b, m_b)
            omega += a @ (D @ b - b @ D)
        return omega

    def fluctuated_dirac(
        self,
        omega: np.ndarray,
        epsilon_prime: int = -1,
    ) -> np.ndarray:
        """D_omega = D + omega + epsilon' * J omega J^{-1}."""
        D = self.dirac_combined()
        U = self.real_structure_combined()
        J_omega_Jinv = U @ np.conj(omega) @ U.T
        return D + omega + epsilon_prime * J_omega_Jinv

    # ------------------------------------------------------------------
    # Fluctuation decomposition
    # ------------------------------------------------------------------

    def decompose_fluctuation(self, omega: np.ndarray) -> dict:
        """Decompose omega into SM gauge and Higgs sectors.

        Reshapes omega from (dim_GV * 32, dim_GV * 32) to per-GV-pair
        32x32 fiber blocks, then extracts:
        - Matter/antimatter sectors (16x16 each)
        - Within matter: lepton (4x4) and quark (12x12) blocks
        - Within each: L-L (gauge), R-R (gauge), L-R / R-L (Higgs)
        - Within quark gauge: color structure (SU(3) content)
        """
        d = self._dim_GV
        f = self._dim_F  # 32

        om_4d = omega.reshape(d, f, d, f).transpose(0, 2, 1, 3)

        matter = om_4d[:, :, 0:16, 0:16]
        antimatter = om_4d[:, :, 16:32, 16:32]
        mat_anti_off = om_4d[:, :, 0:16, 16:32]
        anti_mat_off = om_4d[:, :, 16:32, 0:16]

        # Lepton sector (indices 0:4 within matter)
        lepton = matter[:, :, 0:4, 0:4]
        lepton_gauge_L = lepton[:, :, 0:2, 0:2]
        lepton_gauge_R = lepton[:, :, 2:4, 2:4]
        lepton_higgs_LR = lepton[:, :, 0:2, 2:4]
        lepton_higgs_RL = lepton[:, :, 2:4, 0:2]

        # Quark sector (indices 4:16 within matter)
        quark = matter[:, :, 4:16, 4:16]
        # Quark is 12x12 = (4 flavor) x (3 color)
        # flavor blocks: L (0:6) = (u_L, d_L) x 3 colors, R (6:12) = (u_R, d_R) x 3 colors
        quark_gauge_L = quark[:, :, 0:6, 0:6]
        quark_gauge_R = quark[:, :, 6:12, 6:12]
        quark_higgs_LR = quark[:, :, 0:6, 6:12]
        quark_higgs_RL = quark[:, :, 6:12, 0:6]

        # Lepton-quark off-diagonal (should be zero for SM)
        lepton_quark_off = matter[:, :, 0:4, 4:16]
        quark_lepton_off = matter[:, :, 4:16, 0:4]

        return {
            "matter": matter,
            "antimatter": antimatter,
            "mat_anti_off": mat_anti_off,
            "anti_mat_off": anti_mat_off,
            "lepton": lepton,
            "lepton_gauge_L": lepton_gauge_L,
            "lepton_gauge_R": lepton_gauge_R,
            "lepton_higgs_LR": lepton_higgs_LR,
            "lepton_higgs_RL": lepton_higgs_RL,
            "quark": quark,
            "quark_gauge_L": quark_gauge_L,
            "quark_gauge_R": quark_gauge_R,
            "quark_higgs_LR": quark_higgs_LR,
            "quark_higgs_RL": quark_higgs_RL,
            "lepton_quark_off": lepton_quark_off,
            "quark_lepton_off": quark_lepton_off,
        }

    def higgs_norm(self, omega: np.ndarray) -> float:
        """Frobenius norm of all Higgs (L<->R) blocks."""
        dc = self.decompose_fluctuation(omega)
        return float(np.sqrt(
            np.linalg.norm(dc["lepton_higgs_LR"]) ** 2
            + np.linalg.norm(dc["lepton_higgs_RL"]) ** 2
            + np.linalg.norm(dc["quark_higgs_LR"]) ** 2
            + np.linalg.norm(dc["quark_higgs_RL"]) ** 2
        ))

    def gauge_norm(self, omega: np.ndarray) -> float:
        """Frobenius norm of all gauge (L-L + R-R) blocks."""
        dc = self.decompose_fluctuation(omega)
        return float(np.sqrt(
            np.linalg.norm(dc["lepton_gauge_L"]) ** 2
            + np.linalg.norm(dc["lepton_gauge_R"]) ** 2
            + np.linalg.norm(dc["quark_gauge_L"]) ** 2
            + np.linalg.norm(dc["quark_gauge_R"]) ** 2
        ))

    def extract_color_content(self, omega: np.ndarray) -> dict:
        """Extract SU(3) color content from the quark gauge blocks.

        The quark gauge-L block (d, d, 6, 6) has structure
        (2 flavor) x (3 color) on each side. For a pure SU(3) gauge field,
        the flavor part is proportional to I_2 and the color part is in su(3).
        """
        dc = self.decompose_fluctuation(omega)
        q_L = dc["quark_gauge_L"]  # (d, d, 6, 6)
        q_R = dc["quark_gauge_R"]

        d = self._dim_GV
        # Reshape 6x6 -> (2, 3, 2, 3) to separate flavor and color
        q_L_4d = q_L.reshape(d, d, 2, 3, 2, 3)
        q_R_4d = q_R.reshape(d, d, 2, 3, 2, 3)

        # Color content: trace over flavor (average of diagonal flavor blocks)
        color_L = np.einsum('ijafbg,ab->ijfg', q_L_4d, np.eye(2)) / 2.0
        color_R = np.einsum('ijafbg,ab->ijfg', q_R_4d, np.eye(2)) / 2.0

        # Project onto Gell-Mann basis for SU(3) content
        su3_coeffs_L = np.zeros((d, d, 8), dtype=np.complex128)
        su3_coeffs_R = np.zeros((d, d, 8), dtype=np.complex128)
        for k in range(8):
            su3_coeffs_L[:, :, k] = np.einsum(
                'ijfg,gf->ij', color_L, GELL_MANN[k]
            ) / 2.0  # Tr(lambda_k A) / 2
            su3_coeffs_R[:, :, k] = np.einsum(
                'ijfg,gf->ij', color_R, GELL_MANN[k]
            ) / 2.0

        # U(1) part: trace of color matrix
        u1_L = np.einsum('ijff->ij', color_L) / 3.0
        u1_R = np.einsum('ijff->ij', color_R) / 3.0

        return {
            "color_L": color_L,
            "color_R": color_R,
            "su3_coeffs_L": su3_coeffs_L,
            "su3_coeffs_R": su3_coeffs_R,
            "u1_color_L": u1_L,
            "u1_color_R": u1_R,
        }

    # ------------------------------------------------------------------
    # Axiom verification
    # ------------------------------------------------------------------

    def verify_axioms(self) -> dict:
        """Verify all Connes axioms for the combined triple at finite n_max.

        Returns dict with residuals for each axiom check.
        """
        D = self.dirac_combined()
        U_J = self.real_structure_combined()
        gamma = self.chirality_combined()

        dim = self._dim_H
        I = np.eye(dim, dtype=np.complex128)

        # J^2 = epsilon * I (KO-dim 1: epsilon = -1)
        J_sq = U_J @ np.conj(U_J)
        j_sq_residual = float(np.linalg.norm(J_sq - (-1) * I))

        # JD = epsilon' * DJ (KO-dim 1: epsilon' = +1)
        JD = U_J @ np.conj(D) @ np.linalg.inv(U_J)
        jd_residual = float(np.linalg.norm(JD - D))

        # {gamma, D} = 0 (odd spectral triple)
        gamma_d_anticomm = gamma @ D + D @ gamma
        gamma_d_residual = float(np.linalg.norm(gamma_d_anticomm))

        # gamma^2 = I
        gamma_sq_residual = float(np.linalg.norm(gamma @ gamma - I))

        # J gamma = epsilon'' * gamma J (KO-dim 1: epsilon'' = +1)
        J_gamma = U_J @ np.conj(gamma)
        gamma_J = gamma @ U_J
        j_gamma_residual = float(np.linalg.norm(J_gamma - gamma_J))

        # D Hermitian
        d_herm_residual = float(np.linalg.norm(D - D.conj().T))

        # Order-zero: [a, J b J^{-1}] = 0 for random a, b
        rng = np.random.default_rng(123)
        order_zero_max = 0.0
        I_3 = np.eye(3, dtype=np.complex128)
        for _ in range(20):
            lam_a = rng.standard_normal() + 1j * rng.standard_normal()
            q_a = tuple(rng.standard_normal() for _ in range(4))
            m_a = I_3 + 0.1 * (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))
            lam_b = rng.standard_normal() + 1j * rng.standard_normal()
            q_b = tuple(rng.standard_normal() for _ in range(4))
            m_b = I_3 + 0.1 * (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))

            k_a = rng.integers(0, self.n_gv_multipliers)
            k_b = rng.integers(0, self.n_gv_multipliers)
            a = self.algebra_element(self.gv_multiplier(k_a), lam_a, q_a, m_a)
            b = self.algebra_element(self.gv_multiplier(k_b), lam_b, q_b, m_b)
            JbJinv = U_J @ np.conj(b) @ U_J.T
            comm = a @ JbJinv - JbJinv @ a
            order_zero_max = max(order_zero_max, float(np.linalg.norm(comm)))

        # Order-one: [[D, a], J b J^{-1}] = 0 for random a, b
        order_one_max = 0.0
        for _ in range(20):
            lam_a = rng.standard_normal() + 1j * rng.standard_normal()
            q_a = tuple(rng.standard_normal() for _ in range(4))
            m_a = I_3 + 0.1 * (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))
            lam_b = rng.standard_normal() + 1j * rng.standard_normal()
            q_b = tuple(rng.standard_normal() for _ in range(4))
            m_b = I_3 + 0.1 * (rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))

            k_a = rng.integers(0, self.n_gv_multipliers)
            k_b = rng.integers(0, self.n_gv_multipliers)
            a = self.algebra_element(self.gv_multiplier(k_a), lam_a, q_a, m_a)
            b = self.algebra_element(self.gv_multiplier(k_b), lam_b, q_b, m_b)
            Da_comm = D @ a - a @ D
            JbJinv = U_J @ np.conj(b) @ U_J.T
            order_one = Da_comm @ JbJinv - JbJinv @ Da_comm
            order_one_max = max(order_one_max, float(np.linalg.norm(order_one)))

        return {
            "J_squared_residual": j_sq_residual,
            "JD_relation_residual": jd_residual,
            "gamma_D_anticommutator_residual": gamma_d_residual,
            "gamma_squared_residual": gamma_sq_residual,
            "J_gamma_residual": j_gamma_residual,
            "D_hermitian_residual": d_herm_residual,
            "order_zero_max_residual": order_zero_max,
            "order_one_max_residual": order_one_max,
        }

    # ------------------------------------------------------------------
    # Falsifier: does GeoVac structure force Y = 0?
    # ------------------------------------------------------------------

    def check_natural_negative(
        self,
        n_random_generators: int = 50,
        seed: int = 42,
    ) -> Tuple[bool, str, dict]:
        """Test whether inner fluctuations produce only gauge or also Higgs.

        Same logic as H1's falsifier: if Higgs sector is identically zero
        across random (a, b) pairs, the construction has no Higgs — fails.
        """
        rng = np.random.default_rng(seed)
        n_mult = self.n_gv_multipliers
        I_3 = np.eye(3, dtype=np.complex128)

        higgs_norms = []
        gauge_norms = []
        for _ in range(n_random_generators):
            i_a = rng.integers(0, n_mult)
            i_b = rng.integers(0, n_mult)
            lam_a = (rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
            lam_b = (rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
            q_a = tuple(
                (rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
                for _ in range(4)
            )
            q_b = tuple(
                (rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
                for _ in range(4)
            )
            m_a = I_3 + 0.3 * (
                rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
            )
            m_b = I_3 + 0.3 * (
                rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
            )

            generators = [(
                self.gv_multiplier(i_a), lam_a, q_a, m_a,
                self.gv_multiplier(i_b), lam_b, q_b, m_b,
            )]
            omega = self.inner_fluctuation_one_form(generators)
            higgs_norms.append(self.higgs_norm(omega))
            gauge_norms.append(self.gauge_norm(omega))

        higgs_max = max(higgs_norms)
        gauge_max = max(gauge_norms)
        eps = 1e-12
        higgs_zero = higgs_max < eps * max(gauge_max, 1.0)

        if higgs_zero:
            reason = (
                f"Higgs identically zero (higgs_max={higgs_max:.3e}, "
                f"gauge_max={gauge_max:.3e}): D_F off-diagonal blocked."
            )
        else:
            ratio = higgs_max / gauge_max if gauge_max > 0 else float('inf')
            reason = (
                f"Higgs non-trivial (higgs_max={higgs_max:.3e}, "
                f"gauge_max={gauge_max:.3e}, ratio={ratio:.3f})"
            )

        return higgs_zero, reason, {
            "higgs_max": higgs_max,
            "gauge_max": gauge_max,
            "higgs_mean": float(np.mean(higgs_norms)),
            "gauge_mean": float(np.mean(gauge_norms)),
            "n_samples": n_random_generators,
        }

    def gauge_group_census(
        self, n_samples: int = 30, seed: int = 99
    ) -> dict:
        """Probe the inner fluctuation to identify gauge group content.

        For random (a, b) pairs, compute omega and check:
        1. Does the lepton gauge block have SU(2) structure?
        2. Does the quark gauge block have SU(2) x SU(3) structure?
        3. Is there a U(1) component?

        Returns a census of which gauge sectors are populated.
        """
        rng = np.random.default_rng(seed)
        I_3 = np.eye(3, dtype=np.complex128)

        has_su2_lepton = False
        has_su2_quark = False
        has_su3 = False
        has_u1 = False
        lepton_quark_mixing = False

        for _ in range(n_samples):
            i_a = rng.integers(0, self.n_gv_multipliers)
            i_b = rng.integers(0, self.n_gv_multipliers)
            lam_a = (rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
            lam_b = (rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
            q_a = tuple(
                (rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
                for _ in range(4)
            )
            q_b = tuple(
                (rng.standard_normal() + 1j * rng.standard_normal()) * 0.5
                for _ in range(4)
            )
            m_a = I_3 + 0.3 * (
                rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
            )
            m_b = I_3 + 0.3 * (
                rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3))
            )

            generators = [(
                self.gv_multiplier(i_a), lam_a, q_a, m_a,
                self.gv_multiplier(i_b), lam_b, q_b, m_b,
            )]
            omega = self.inner_fluctuation_one_form(generators)
            dc = self.decompose_fluctuation(omega)

            eps = 1e-10
            if np.linalg.norm(dc["lepton_gauge_L"]) > eps:
                has_su2_lepton = True
            if np.linalg.norm(dc["quark_gauge_L"]) > eps:
                has_su2_quark = True
            if np.linalg.norm(dc["lepton_gauge_R"]) > eps:
                has_u1 = True

            color = self.extract_color_content(omega)
            for k in range(8):
                if np.linalg.norm(color["su3_coeffs_L"][:, :, k]) > eps:
                    has_su3 = True
                    break

            if np.linalg.norm(dc["lepton_quark_off"]) > eps:
                lepton_quark_mixing = True

        return {
            "SU2_lepton": has_su2_lepton,
            "SU2_quark": has_su2_quark,
            "SU3": has_su3,
            "U1": has_u1,
            "lepton_quark_mixing": lepton_quark_mixing,
            "gauge_group": _format_gauge_group(
                has_su2_lepton, has_su2_quark, has_su3, has_u1
            ),
        }


def _format_gauge_group(
    su2_l: bool, su2_q: bool, su3: bool, u1: bool
) -> str:
    parts = []
    if u1:
        parts.append("U(1)")
    if su2_l or su2_q:
        parts.append("SU(2)")
    if su3:
        parts.append("SU(3)")
    return " x ".join(parts) if parts else "trivial"


# ======================================================================
# Convenience constructors
# ======================================================================

def standard_model_triple(
    n_max: int,
    yukawa_e: complex = 0.0,
    yukawa_nu: complex = 0.0,
    yukawa_u: complex = 0.0,
    yukawa_d: complex = 0.0,
) -> StandardModelACTriple:
    """Build the full SM almost-commutative extension at given n_max."""
    finite = StandardModelFiniteTriple(
        yukawa_nu=yukawa_nu,
        yukawa_e=yukawa_e,
        yukawa_u=yukawa_u,
        yukawa_d=yukawa_d,
    )
    return StandardModelACTriple(n_max=n_max, finite=finite)


def sm_gauge_only(n_max: int) -> StandardModelACTriple:
    """SM triple with zero Yukawa (gauge sector only, no Higgs)."""
    return standard_model_triple(n_max)


def sm_one_generation(
    n_max: int,
    yukawa_e: float = 0.00029,
    yukawa_u: float = 0.0012,
    yukawa_d: float = 0.0026,
) -> StandardModelACTriple:
    """SM triple with physical-order one-generation Yukawa couplings.

    Default values are rough m/v ratios (mass / Higgs vev) for first
    generation. These are IMPOSED inputs (Paper 18 calibration tier).
    """
    return standard_model_triple(
        n_max,
        yukawa_e=yukawa_e,
        yukawa_nu=0.0,
        yukawa_u=yukawa_u,
        yukawa_d=yukawa_d,
    )
