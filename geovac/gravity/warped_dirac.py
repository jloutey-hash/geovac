"""Warped Dirac operator on the discrete cigar substrate.

Sprint G4-4a first move (2026-05-28): constant-warp Dirac on
D^2 x S^2 at sprint scale. Implements the load-bearing falsifiers
F1 (factorization at constant warp) and F2 (chirality grading) as
operator-level identities. Variable-warp (G4-4b), conical-defect
(G4-4c), and continuum recovery at small t (G4-4a falsifier F3) are
deferred to subsequent weeks of the G4-4a 4-8 week sprint.

Architecture
------------
The continuum Dirac on the cigar at constant warp r(rho) = r_h:

  D_cigar = D_{D^2} (x) I_{S^2} + gamma^5_{D^2} (x) D_{S^2} / r_h

where gamma^5_{D^2} = -i gamma^1 gamma^2 = sigma_3 is the 2D chirality
matrix. The squared operator:

  D_cigar^2 = D_{D^2}^2 (x) I + I (x) D_{S^2}^2 / r_h^2

(cross term cancels because {gamma^5_{D^2}, D_{D^2}} = 0).

On flat polar 2D, D_{D^2}^2 = -nabla^2 (x) I_2 (rank-2 spinor bundle).
The angular sector uses anti-periodic BC (fermionic), so the m_eff
spectrum is half-integer-shifted.

S^2 spinor spectrum (Camporesi-Higuchi 1996, Paper 28 sec 4.14):

  |lambda_n^{S^2}| = (n + 1) / r_h
  g_n^Dirac = 4 (n + 1)

The cigar squared-Dirac spectrum is the outer sum of the disk-Dirac
spectrum and the S^2-Dirac spectrum; the heat trace factorizes
identically (load-bearing falsifier F1).

Conventions
-----------
- Pauli matrices: gamma^1 = sigma_1, gamma^2 = sigma_2,
  gamma^5 = -i gamma^1 gamma^2 = sigma_3.
- 2D Cl(2, 0) algebra: {gamma^a, gamma^b} = 2 delta^{ab} I_2.
- Anti-periodic phi BC: fermionic, m_eff = (k + 1/2) shifted index.
- Hermitian polar Laplacian via u = sqrt(rho) f substitution
  (G4-3a-cleanup convention).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

import numpy as np

# =================================================================
# 2D Cl(2, 0) gamma matrices and chirality
# =================================================================

GAMMA_2D_1 = np.array([[0, 1], [1, 0]], dtype=complex)
"""sigma_1 = gamma^1 in 2D Cl(2, 0)."""

GAMMA_2D_2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
"""sigma_2 = gamma^2 in 2D Cl(2, 0)."""

GAMMA_2D_5 = np.array([[1, 0], [0, -1]], dtype=complex)
"""sigma_3 = -i gamma^1 gamma^2 = chirality in 2D Cl(2, 0)."""

I_2 = np.eye(2, dtype=complex)
"""2x2 identity (spinor identity)."""


def verify_gamma_algebra_2d(tol: float = 1e-14) -> Dict[str, float]:
    """Verify the Cl(2, 0) gamma algebra identities.

    Returns
    -------
    dict
        Maximum absolute residual for each algebraic identity. Each
        value should be ~ 1e-16 (machine precision); failures
        indicate a convention bug in the gamma matrix definitions.
    """
    g1, g2, g5 = GAMMA_2D_1, GAMMA_2D_2, GAMMA_2D_5
    results = {
        "g1_squared_is_I": float(np.max(np.abs(g1 @ g1 - I_2))),
        "g2_squared_is_I": float(np.max(np.abs(g2 @ g2 - I_2))),
        "g5_squared_is_I": float(np.max(np.abs(g5 @ g5 - I_2))),
        "g1_g2_anticom": float(np.max(np.abs(g1 @ g2 + g2 @ g1))),
        "g5_definition": float(np.max(np.abs(g5 - (-1j * g1 @ g2)))),
        "g5_g1_anticom": float(np.max(np.abs(g5 @ g1 + g1 @ g5))),
        "g5_g2_anticom": float(np.max(np.abs(g5 @ g2 + g2 @ g5))),
    }
    return results


# =================================================================
# Discrete disk Dirac (2D polar)
# =================================================================


@dataclass
class DiscreteDiskDirac:
    """Squared 2D disk Dirac on the discrete polar substrate.

    Parameters
    ----------
    N_rho : int
        Number of radial sites; rho_k = k * a for k = 1, ..., N_rho.
        The apex rho = 0 is excluded (singular point on flat polar 2D).
    a : float
        Radial lattice spacing. IR cutoff R = N_rho * a.
    N_phi : int
        Number of azimuthal sites; the phi-circle is fermionic
        anti-periodic with h_phi = 2 pi / N_phi.

    Notes
    -----
    The disk-Dirac squared operator on flat polar 2D is
        D_{D^2}^2 = - nabla^2 (x) I_2
    where the rank-2 spinor bundle comes from the 2D Cl(2, 0) Clifford
    algebra. The scalar Laplacian on the polar disk has the form

        nabla^2 = partial_rho^2 + (1/rho) partial_rho
                                + (1/rho^2) partial_phi^2

    which we discretize via the symmetric Hermitian polar Laplacian
    (G4-3a-cleanup convention with u = sqrt(rho) f).

    The anti-periodic phi BC (fermionic) shifts the discrete azimuthal
    Laplacian eigenvalues to

        lambda_k^{anti-per} = (2 / h_phi)^2 * sin^2(pi (k + 1/2) / N_phi)

    Half-integer continuum limit at small k: m_eff -> k + 1/2.
    """

    N_rho: int
    a: float
    N_phi: int

    def __post_init__(self) -> None:
        if self.N_rho < 1:
            raise ValueError(f"N_rho must be >= 1, got {self.N_rho}")
        if self.a <= 0:
            raise ValueError(f"a must be > 0, got {self.a}")
        if self.N_phi < 2:
            raise ValueError(f"N_phi must be >= 2, got {self.N_phi}")

    @property
    def h_phi(self) -> float:
        return 2 * np.pi / self.N_phi

    @property
    def R(self) -> float:
        """IR cutoff (largest radial site)."""
        return self.N_rho * self.a

    @property
    def hilbert_dim(self) -> int:
        """Spinor Hilbert space dimension: 2 (spin) * N_rho * N_phi."""
        return 2 * self.N_rho * self.N_phi

    def _hermitian_radial_laplacian(self, m_eff: float) -> np.ndarray:
        """Symmetric tridiagonal radial Laplacian on L^2(rho drho).

        Uses u = sqrt(rho) f substitution from G4-3a-cleanup. The
        centrifugal term becomes (m_eff^2 - 1/4) / rho^2 in the
        u-representation (NOT m_eff^2).
        """
        N = self.N_rho
        a = self.a
        H = np.zeros((N, N))
        for i in range(N):
            k = i + 1
            rho = k * a
            H[i, i] = 2.0 / a**2 + (m_eff**2 - 0.25) / rho**2
            if i > 0:
                H[i, i - 1] = -1.0 / a**2
            if i < N - 1:
                H[i, i + 1] = -1.0 / a**2
        return H

    def squared_eigenvalues(self) -> np.ndarray:
        """All eigenvalues of D_{D^2}^2 on the spinor bundle.

        The block structure: D^2 = -nabla^2 (x) I_2, so each scalar
        Laplacian eigenvalue appears with multiplicity 2 in the
        spinor spectrum.
        """
        scalar_eigs: List[float] = []
        for k_idx in range(self.N_phi):
            if k_idx <= self.N_phi // 2:
                k = k_idx
            else:
                k = k_idx - self.N_phi
            # Fermionic anti-periodic: shift k -> k + 1/2
            m_eff_sq = (2.0 / self.h_phi) ** 2 * np.sin(
                np.pi * (k + 0.5) / self.N_phi
            ) ** 2
            m_eff = float(np.sqrt(m_eff_sq))
            H_rad = self._hermitian_radial_laplacian(m_eff)
            evals = np.linalg.eigvalsh(H_rad)
            scalar_eigs.extend(evals.tolist())
        scalar_arr = np.array(scalar_eigs)
        # Rank-2 spinor bundle: each scalar eigenvalue doubled
        return np.array(sorted(np.concatenate([scalar_arr, scalar_arr])))

    def heat_trace(self, t: float) -> float:
        """Tr exp(-t D^2) for the disk-Dirac squared operator."""
        eigs = self.squared_eigenvalues()
        return float(np.sum(np.exp(-eigs * t)))


@dataclass
class DiscreteWedgeDirac:
    """Squared 2D Dirac on a wedge with apex angle 2*pi*alpha (G4-4c).

    Same construction as DiscreteDiskDirac but with azimuthal period
    2*pi*alpha instead of 2*pi. Fermionic anti-periodic BC at the seam.

    The conical-defect parameter alpha:
      - alpha = 1: smooth disk (no defect); reduces bit-exact to
        DiscreteDiskDirac at matching (N_rho, a, N_phi).
      - alpha < 1: deficit angle 2*pi*(1-alpha) > 0 (pointed cone).
      - alpha > 1: excess angle > 0 (saddle-like cone / covering).

    Parameters
    ----------
    N_rho : int
        Number of radial sites; rho_k = k*a.
    a : float
        Radial lattice spacing.
    N_phi : int
        Number of azimuthal sites on the wedge.
    alpha : float
        Apex-angle multiple. Apex angle = 2*pi*alpha. Must be > 0.

    For the "proper wedge lattice" convention (G4-3c-proper / T1, T3):
    fix h_phi = 2*pi/N_0 (a reference smooth-disk count), and choose
    N_phi(alpha) = alpha * N_0. The caller should enforce this if
    direct comparison with constant-h_phi scaling is desired.

    Anti-periodic azimuthal eigenvalues:
        h_phi = 2*pi*alpha / N_phi
        lambda_k^{disc} = (2/h_phi)^2 * sin^2(pi (k+1/2) / N_phi)
    Continuum: lambda_k -> ((k+1/2)/alpha)^2 (half-integer / alpha).
    """

    N_rho: int
    a: float
    N_phi: int
    alpha: float

    def __post_init__(self) -> None:
        if self.N_rho < 1:
            raise ValueError(f"N_rho must be >= 1, got {self.N_rho}")
        if self.a <= 0:
            raise ValueError(f"a must be > 0, got {self.a}")
        if self.N_phi < 2:
            raise ValueError(f"N_phi must be >= 2, got {self.N_phi}")
        if self.alpha <= 0:
            raise ValueError(f"alpha must be > 0, got {self.alpha}")

    @property
    def h_phi(self) -> float:
        return 2 * np.pi * self.alpha / self.N_phi

    @property
    def R(self) -> float:
        return self.N_rho * self.a

    @property
    def hilbert_dim(self) -> int:
        """Spinor Hilbert space dim = 2 * N_rho * N_phi."""
        return 2 * self.N_rho * self.N_phi

    def _hermitian_radial_laplacian(self, m_eff: float) -> np.ndarray:
        """Same convention as DiscreteDiskDirac."""
        N = self.N_rho
        a = self.a
        H = np.zeros((N, N))
        for i in range(N):
            k = i + 1
            rho = k * a
            H[i, i] = 2.0 / a**2 + (m_eff**2 - 0.25) / rho**2
            if i > 0:
                H[i, i - 1] = -1.0 / a**2
            if i < N - 1:
                H[i, i + 1] = -1.0 / a**2
        return H

    def squared_eigenvalues(self) -> np.ndarray:
        """All eigenvalues of D^2 on the wedge spinor bundle.

        Same rank-2 doubling as DiscreteDiskDirac. At alpha = 1, this
        reduces bit-exact to DiscreteDiskDirac.squared_eigenvalues.
        """
        scalar_eigs: List[float] = []
        for k_idx in range(self.N_phi):
            if k_idx <= self.N_phi // 2:
                k = k_idx
            else:
                k = k_idx - self.N_phi
            m_eff_sq = (2.0 / self.h_phi) ** 2 * np.sin(
                np.pi * (k + 0.5) / self.N_phi
            ) ** 2
            m_eff = float(np.sqrt(m_eff_sq))
            H_rad = self._hermitian_radial_laplacian(m_eff)
            evals = np.linalg.eigvalsh(H_rad)
            scalar_eigs.extend(evals.tolist())
        scalar_arr = np.array(scalar_eigs)
        return np.array(sorted(np.concatenate([scalar_arr, scalar_arr])))

    def heat_trace(self, t: float) -> float:
        eigs = self.squared_eigenvalues()
        return float(np.sum(np.exp(-eigs * t)))


# ---------------------------------------------------------------------
# Spectral azimuthal discretization (G4-6d)
# ---------------------------------------------------------------------
#
# Per the v3.20.0 task #28 finding and the G4-6a first-move reframing
# (debug/g4_6a_multi_substrate_uv_first_move_memo.md):
# the FD azimuthal Laplacian underestimates the UV target 1/(24*pi*t)
# by factor ~1/(4*pi^2) at the truncation edge. The spectral (DST/Fourier)
# discretization replaces the FD eigenvalue
#     m_eff^{FD} = (2 / h_phi) * sin(pi (k + 1/2) / N_phi)
# with the exact continuum eigenvalue
#     m_eff^{spec} = (k + 1/2) / alpha.
# This recovers the UV asymptote at sub-percent precision at substrate
# scales where FD recovers <1% (task #28 measured 25.6% vs 0.04%).


@dataclass
class DiscreteDiskDiracSpectral:
    """Disk-Dirac with SPECTRAL (exact) azimuthal eigenvalues (G4-6d).

    Same radial structure as DiscreteDiskDirac (Hermitian polar
    Laplacian with u = sqrt(rho) f substitution); replaces the FD
    azimuthal eigenvalues with the exact continuum m_eff = k + 1/2.

    Parameters
    ----------
    N_rho : int
        Number of radial sites.
    a : float
        Radial lattice spacing.
    N_phi : int
        Azimuthal mode count truncation. Modes are k = -N_phi/2, ...,
        N_phi/2 - 1 with eigenvalues m_eff = k + 1/2.
    """

    N_rho: int
    a: float
    N_phi: int

    def __post_init__(self) -> None:
        if self.N_rho < 1:
            raise ValueError(f"N_rho must be >= 1, got {self.N_rho}")
        if self.a <= 0:
            raise ValueError(f"a must be > 0, got {self.a}")
        if self.N_phi < 2:
            raise ValueError(f"N_phi must be >= 2, got {self.N_phi}")

    @property
    def R(self) -> float:
        return self.N_rho * self.a

    @property
    def hilbert_dim(self) -> int:
        return 2 * self.N_rho * self.N_phi

    def _hermitian_radial_laplacian(self, m_eff: float) -> np.ndarray:
        """Same as DiscreteDiskDirac._hermitian_radial_laplacian."""
        N = self.N_rho
        a = self.a
        H = np.zeros((N, N))
        for i in range(N):
            k = i + 1
            rho = k * a
            H[i, i] = 2.0 / a**2 + (m_eff**2 - 0.25) / rho**2
            if i > 0:
                H[i, i - 1] = -1.0 / a**2
            if i < N - 1:
                H[i, i + 1] = -1.0 / a**2
        return H

    def squared_eigenvalues(self) -> np.ndarray:
        """Eigenvalues with EXACT azimuthal m_eff = k + 1/2."""
        scalar_eigs: List[float] = []
        for k_idx in range(self.N_phi):
            if k_idx <= self.N_phi // 2:
                k = k_idx
            else:
                k = k_idx - self.N_phi
            m_eff = float(k + 0.5)  # SPECTRAL: exact continuum eigenvalue
            H_rad = self._hermitian_radial_laplacian(m_eff)
            evals = np.linalg.eigvalsh(H_rad)
            scalar_eigs.extend(evals.tolist())
        scalar_arr = np.array(scalar_eigs)
        return np.array(sorted(np.concatenate([scalar_arr, scalar_arr])))

    def heat_trace(self, t: float) -> float:
        eigs = self.squared_eigenvalues()
        return float(np.sum(np.exp(-eigs * t)))


@dataclass
class DiscreteWedgeDiracSpectral:
    """Wedge-Dirac with SPECTRAL (exact) azimuthal eigenvalues (G4-6d).

    Same radial structure as DiscreteWedgeDirac; replaces the FD
    azimuthal eigenvalues with the exact continuum m_eff = (k + 1/2)/alpha.

    Parameters
    ----------
    N_rho : int
        Number of radial sites.
    a : float
        Radial lattice spacing.
    N_phi : int
        Azimuthal mode count truncation. Modes are k = -N_phi/2, ...,
        N_phi/2 - 1 with eigenvalues m_eff = (k + 1/2) / alpha.
    alpha : float
        Apex-angle multiple. Apex angle = 2*pi*alpha. Must be > 0.

    F6 bit-exact check: at alpha = 1, this class reduces to
    DiscreteDiskDiracSpectral at matching (N_rho, a, N_phi).
    """

    N_rho: int
    a: float
    N_phi: int
    alpha: float

    def __post_init__(self) -> None:
        if self.N_rho < 1:
            raise ValueError(f"N_rho must be >= 1, got {self.N_rho}")
        if self.a <= 0:
            raise ValueError(f"a must be > 0, got {self.a}")
        if self.N_phi < 2:
            raise ValueError(f"N_phi must be >= 2, got {self.N_phi}")
        if self.alpha <= 0:
            raise ValueError(f"alpha must be > 0, got {self.alpha}")

    @property
    def R(self) -> float:
        return self.N_rho * self.a

    @property
    def hilbert_dim(self) -> int:
        return 2 * self.N_rho * self.N_phi

    def _hermitian_radial_laplacian(self, m_eff: float) -> np.ndarray:
        """Same as DiscreteWedgeDirac._hermitian_radial_laplacian."""
        N = self.N_rho
        a = self.a
        H = np.zeros((N, N))
        for i in range(N):
            k = i + 1
            rho = k * a
            H[i, i] = 2.0 / a**2 + (m_eff**2 - 0.25) / rho**2
            if i > 0:
                H[i, i - 1] = -1.0 / a**2
            if i < N - 1:
                H[i, i + 1] = -1.0 / a**2
        return H

    def squared_eigenvalues(self) -> np.ndarray:
        """Eigenvalues with EXACT azimuthal m_eff = (k + 1/2)/alpha."""
        scalar_eigs: List[float] = []
        for k_idx in range(self.N_phi):
            if k_idx <= self.N_phi // 2:
                k = k_idx
            else:
                k = k_idx - self.N_phi
            m_eff = float((k + 0.5) / self.alpha)  # SPECTRAL exact eigenvalue
            H_rad = self._hermitian_radial_laplacian(m_eff)
            evals = np.linalg.eigvalsh(H_rad)
            scalar_eigs.extend(evals.tolist())
        scalar_arr = np.array(scalar_eigs)
        return np.array(sorted(np.concatenate([scalar_arr, scalar_arr])))

    def heat_trace(self, t: float) -> float:
        eigs = self.squared_eigenvalues()
        return float(np.sum(np.exp(-eigs * t)))


# Companion: scalar disk Laplacian with periodic phi BC for F3 ratio test


@dataclass
class DiscreteDiskScalar:
    """Scalar Laplacian on the discrete polar disk (periodic phi BC).

    Used for the F3-rough check that D_Dirac^2 / D_scalar -> rank
    factor 2 at small t in the continuum limit. This is the bosonic
    sibling of DiscreteDiskDirac.
    """

    N_rho: int
    a: float
    N_phi: int

    @property
    def h_phi(self) -> float:
        return 2 * np.pi / self.N_phi

    def _hermitian_radial_laplacian(self, m_eff: float) -> np.ndarray:
        N = self.N_rho
        a = self.a
        H = np.zeros((N, N))
        for i in range(N):
            k = i + 1
            rho = k * a
            H[i, i] = 2.0 / a**2 + (m_eff**2 - 0.25) / rho**2
            if i > 0:
                H[i, i - 1] = -1.0 / a**2
            if i < N - 1:
                H[i, i + 1] = -1.0 / a**2
        return H

    def eigenvalues(self) -> np.ndarray:
        eigs: List[float] = []
        for k_idx in range(self.N_phi):
            if k_idx <= self.N_phi // 2:
                k = k_idx
            else:
                k = k_idx - self.N_phi
            # Periodic bosonic: integer m
            m_eff_sq = (2.0 / self.h_phi) ** 2 * np.sin(
                np.pi * k / self.N_phi
            ) ** 2
            m_eff = float(np.sqrt(m_eff_sq))
            H_rad = self._hermitian_radial_laplacian(m_eff)
            evals = np.linalg.eigvalsh(H_rad)
            eigs.extend(evals.tolist())
        return np.array(sorted(eigs))

    def heat_trace(self, t: float) -> float:
        eigs = self.eigenvalues()
        return float(np.sum(np.exp(-eigs * t)))


# =================================================================
# S^2 Dirac spectrum (Camporesi-Higuchi 1996)
# =================================================================


@dataclass
class S2DiracSpectrum:
    """Camporesi-Higuchi spinor spectrum on S^2_{r_h}.

    Eigenvalues |lambda_n^{S^2}| = (n + 1) / r_h, n = 0, 1, ...,
    multiplicity g_n^Dirac = 4 (n + 1). Truncated at n <= l_max.

    Squared spectrum: lambda_n^2 = (n + 1)^2 / r_h^2 with multiplicity
    8 (n + 1) (counting +/- eigenvalues, each appearing 4(n+1) times).
    """

    l_max: int
    r_h: float = 1.0

    def __post_init__(self) -> None:
        if self.l_max < 0:
            raise ValueError(f"l_max must be >= 0, got {self.l_max}")
        if self.r_h <= 0:
            raise ValueError(f"r_h must be > 0, got {self.r_h}")

    def squared_eigenvalues(self) -> np.ndarray:
        """All squared eigenvalues of D_{S^2} / r_h.

        Each |lambda_n|^2 = ((n+1)/r_h)^2 appears with multiplicity
        2 * g_n^Dirac = 8(n+1) (signs +/- of lambda each give 4(n+1)).
        """
        eigs: List[float] = []
        for n in range(self.l_max + 1):
            lam_sq = ((n + 1) / self.r_h) ** 2
            mult = 2 * 4 * (n + 1)  # both signs, each multiplicity 4(n+1)
            eigs.extend([lam_sq] * mult)
        return np.array(eigs)

    def n_modes(self) -> int:
        return sum(2 * 4 * (n + 1) for n in range(self.l_max + 1))


# =================================================================
# Constant-warp tensor product
# =================================================================


@dataclass
class WarpedDiracConstant:
    """Squared Dirac on the cigar at constant warp r(rho) = r_h.

    Parameters
    ----------
    disk : DiscreteDiskDirac
        Discrete 2D disk-Dirac substrate.
    sphere : S2DiracSpectrum
        Camporesi-Higuchi S^2 spinor spectrum at radius r_h.

    The squared Dirac on the cigar at constant warp:

        D_cigar^2 = D_{D^2}^2 (x) I_{S^2} + I_{D^2} (x) D_{S^2}^2

    Total Hilbert space dimension:
        dim H_cigar = (2 N_rho N_phi) * (sum_n 8(n+1))

    The factorization at the heat-trace level (load-bearing falsifier
    F1) is the outer-sum identity:

        K_cigar(t) = K_{D^2}(t) * K_{S^2}(t)
                   = Sum_{i, j} exp(-(a_i + b_j) t)
                   = Sum_i exp(-a_i t) * Sum_j exp(-b_j t)
    """

    disk: DiscreteDiskDirac
    sphere: S2DiracSpectrum

    @property
    def hilbert_dim(self) -> int:
        return self.disk.hilbert_dim * self.sphere.n_modes()

    def squared_eigenvalues(self) -> np.ndarray:
        """All eigenvalues of D_cigar^2 via direct outer sum.

        Beware: total mode count can be very large; for sprint-scale
        parameters this returns ~ 10^6 - 10^7 eigenvalues.
        """
        disk_eigs = self.disk.squared_eigenvalues()
        s2_eigs = self.sphere.squared_eigenvalues()
        return (disk_eigs[:, None] + s2_eigs[None, :]).ravel()

    def heat_trace_factorized(self, t: float) -> float:
        """K_cigar(t) via factorization (cheap path)."""
        K_disk = self.disk.heat_trace(t)
        s2_eigs = self.sphere.squared_eigenvalues()
        K_S2 = float(np.sum(np.exp(-s2_eigs * t)))
        return K_disk * K_S2

    def heat_trace_direct(self, t: float) -> float:
        """K_cigar(t) via direct sum over joint eigenvalues.

        Used to verify F1 (factorization) against the factorized path.
        Memory-bounded by hilbert_dim; use heat_trace_factorized for
        large parameters.
        """
        return float(np.sum(np.exp(-self.squared_eigenvalues() * t)))


# =================================================================
# Load-bearing falsifiers
# =================================================================


def verify_F1_factorization(
    disk: DiscreteDiskDirac,
    sphere: S2DiracSpectrum,
    t_values: List[float],
    tol: float = 1e-10,
) -> Dict[str, Dict]:
    """F1 falsifier: K_cigar(t) = K_disk(t) * K_{S^2}(t) at constant warp.

    Verifies the outer-sum identity numerically at each t_value.
    Should pass to ~ 1e-13 (float64 sum precision over ~ 10^4 - 10^5
    eigenvalues) at sprint-scale parameters.

    Returns
    -------
    dict
        Per-t results with K_factorized, K_direct, rel_err, passed.
        Plus aggregate 'all_passed' boolean.
    """
    warp = WarpedDiracConstant(disk, sphere)
    results: Dict[str, Dict] = {}
    all_passed = True
    for t in t_values:
        K_fact = warp.heat_trace_factorized(t)
        K_dir = warp.heat_trace_direct(t)
        denom = abs(K_dir) if K_dir != 0 else 1.0
        rel_err = abs(K_fact - K_dir) / denom
        passed = bool(rel_err < tol)
        all_passed = all_passed and passed
        results[str(t)] = {
            "K_factorized": K_fact,
            "K_direct": K_dir,
            "abs_err": float(abs(K_fact - K_dir)),
            "rel_err": float(rel_err),
            "passed": passed,
        }
    results["all_passed"] = all_passed
    results["tolerance"] = tol
    return results


def verify_F2_chirality(tol: float = 1e-14) -> Dict[str, Dict]:
    """F2 falsifier: {gamma^5, D_disk} = 0 at the gamma-matrix level.

    Since D_disk = gamma^1 partial_rho + gamma^2 (1/rho) nabla_phi and
    gamma^5 anticommutes with both gamma^1 and gamma^2, the chirality
    grading is verified at the algebraic (Pauli matrix) level. This is
    a bit-exact identity for the 2D Cl(2, 0) algebra.

    Returns
    -------
    dict
        algebra : dict of gamma algebra residuals
        passed : boolean (all residuals < tol)
    """
    algebra = verify_gamma_algebra_2d()
    passed = bool(
        all(v < tol for v in algebra.values())
    )
    return {"algebra": algebra, "passed": passed, "tolerance": tol}


@dataclass
class DiscreteDirac2D:
    """Explicit linear Dirac operator on the 2D disk (G4-4a week 2).

    Construction: canonical chirality-graded form. In each anti-periodic
    azimuthal Fourier mode k, the Dirac block is

        D_k = [[0, sqrt(L_k)], [sqrt(L_k), 0]]

    where L_k is the Hermitian radial Laplacian (from G4-3a-cleanup)
    at the half-integer-shifted m_eff for mode k.

    Properties:
      - D is Hermitian (sqrt(L_k) is Hermitian positive)
      - D^2 is block-diagonal with each L_k eigenvalue doubled
      - {gamma^5, D} = 0 where gamma^5 = diag(I, -I) on the spin grading
      - Spectrum: +/- sqrt(mu_n) for each eigenvalue mu_n of L_k

    The two spin components correspond to gamma^5 = +1 (upper) and -1
    (lower) Weyl spinor subspaces. This is operator-level Dirac, not
    just the squared spectrum.
    """

    N_rho: int
    a: float
    N_phi: int

    @property
    def h_phi(self) -> float:
        return 2 * np.pi / self.N_phi

    @property
    def R(self) -> float:
        return self.N_rho * self.a

    @property
    def hilbert_dim(self) -> int:
        return 2 * self.N_rho * self.N_phi

    def _hermitian_radial_laplacian(self, m_eff: float) -> np.ndarray:
        """Same convention as DiscreteDiskDirac."""
        N = self.N_rho
        a = self.a
        H = np.zeros((N, N))
        for i in range(N):
            k = i + 1
            rho = k * a
            H[i, i] = 2.0 / a**2 + (m_eff**2 - 0.25) / rho**2
            if i > 0:
                H[i, i - 1] = -1.0 / a**2
            if i < N - 1:
                H[i, i + 1] = -1.0 / a**2
        return H

    def fourier_modes(self) -> List[float]:
        """List of m_eff values for each anti-periodic Fourier mode."""
        modes = []
        for k_idx in range(self.N_phi):
            if k_idx <= self.N_phi // 2:
                k = k_idx
            else:
                k = k_idx - self.N_phi
            m_eff_sq = (2.0 / self.h_phi) ** 2 * np.sin(
                np.pi * (k + 0.5) / self.N_phi
            ) ** 2
            modes.append(float(np.sqrt(m_eff_sq)))
        return modes

    def L_k(self, k_idx: int) -> np.ndarray:
        """Scalar Laplacian block in Fourier mode k_idx."""
        modes = self.fourier_modes()
        return self._hermitian_radial_laplacian(modes[k_idx])

    def sqrt_L_k(self, k_idx: int) -> np.ndarray:
        """Hermitian positive square root of L_k."""
        L = self.L_k(k_idx)
        # Eigendecomposition
        evals, V = np.linalg.eigh(L)
        # Clip tiny negative eigenvalues (numerical noise)
        evals_clipped = np.maximum(evals, 0.0)
        sqrt_evals = np.sqrt(evals_clipped)
        # sqrt(L) = V diag(sqrt) V^dagger
        return (V * sqrt_evals) @ V.T.conj()

    def dirac_block(self, k_idx: int) -> np.ndarray:
        """Dirac block in Fourier mode k_idx: D_k = [[0, sqrt(L)], [sqrt(L), 0]]."""
        sL = self.sqrt_L_k(k_idx)
        N = self.N_rho
        D = np.zeros((2 * N, 2 * N), dtype=complex)
        D[:N, N:] = sL
        D[N:, :N] = sL
        return D

    def chirality_block(self) -> np.ndarray:
        """gamma^5 in each Fourier-mode block: diag(I_N, -I_N)."""
        N = self.N_rho
        g5 = np.zeros((2 * N, 2 * N), dtype=complex)
        g5[:N, :N] = np.eye(N)
        g5[N:, N:] = -np.eye(N)
        return g5

    def eigenvalues(self) -> np.ndarray:
        """Signed eigenvalues of D across all Fourier modes.

        For each k_idx: eigenvalues are +/- sqrt(mu_n) for each mu_n
        in the spectrum of L_k. Total = 2 * N_rho * N_phi.
        """
        all_evals = []
        for k_idx in range(self.N_phi):
            L = self.L_k(k_idx)
            mu = np.linalg.eigvalsh(L)
            mu_clipped = np.maximum(mu, 0.0)
            sm = np.sqrt(mu_clipped)
            # +/- sqrt(mu) for each scalar eigenvalue
            all_evals.extend(sm.tolist())
            all_evals.extend((-sm).tolist())
        return np.array(sorted(all_evals))

    def verify_hermitian_per_block(self, tol: float = 1e-12) -> bool:
        """Verify each Dirac block is Hermitian."""
        for k_idx in range(self.N_phi):
            D_k = self.dirac_block(k_idx)
            residual = np.max(np.abs(D_k - D_k.T.conj()))
            if residual > tol:
                return False
        return True

    def verify_chirality_anticommute_per_block(self, tol: float = 1e-12) -> Dict:
        """Verify {gamma^5, D_k} = 0 for each Fourier block."""
        g5 = self.chirality_block()
        results = {}
        all_residuals = []
        for k_idx in range(self.N_phi):
            D_k = self.dirac_block(k_idx)
            anticom = g5 @ D_k + D_k @ g5
            residual = float(np.max(np.abs(anticom)))
            all_residuals.append(residual)
            results[k_idx] = residual
        return {
            "per_block": results,
            "max_residual": float(np.max(all_residuals)),
            "passed": bool(np.max(all_residuals) < tol),
        }

    def verify_D_squared_matches_factorized(
        self, factorized_eigs: np.ndarray, tol: float = 1e-10
    ) -> Dict:
        """Verify D^2 spectrum matches the factorized scalar spectrum.

        Each Fourier-mode block contributes mu_n eigenvalues of L_k,
        each appearing TWICE (rank-2 spinor bundle).
        """
        # Compute D^2 eigenvalues from explicit D
        D_eigs = self.eigenvalues()
        D_sq_eigs = np.sort(D_eigs**2)
        # Compare to factorized
        factorized_sorted = np.sort(factorized_eigs)
        max_diff = float(np.max(np.abs(D_sq_eigs - factorized_sorted)))
        return {
            "max_diff": max_diff,
            "passed": bool(max_diff < tol),
        }

    def lowest_positive_eigenvalue(self) -> float:
        """Smallest positive |D| eigenvalue.

        Continuum prediction at the half-integer ground state m_eff = 1/2:
            |lambda_min| -> j_{1/2, 1} / R = pi / R
        """
        evals = self.eigenvalues()
        pos = evals[evals > 0]
        return float(np.min(pos))


@dataclass
class VariableWarpDirac:
    """Squared Dirac on cigar at variable warp r(rho), Level 1 approximation.

    Spectrum-level Level 1 per G4-4b scoping memo: includes the
    leading position-dependent S^2 mass (n+1)^2 / r(rho_k)^2 at each
    radial site, but DEFERS the spin-connection cross term
    r'(rho)/r(rho) * gamma^rho to G4-4b-d (Level 2 operator-level).

    Squared operator (Level 1):

        D_cigar^2 = (D_{D^2})^2 (x) I + I (x) D_{S^2}^2 / r(rho)^2

    At constant warp r(rho) = r_h, this reduces bit-exact to G4-4a's
    WarpedDiracConstant (load-bearing falsifier F6).

    Per (S^2 mode n, azimuthal Fourier mode k_phi), the radial
    Hamiltonian is

        H_{n, k_phi} = L_disk(k_phi) + diag((n+1)^2 / r(rho_k)^2)

    Each eigenvalue carries multiplicity 16(n+1):
        2 (rank-2 spinor on disk) * 8(n+1) (S^2 Dirac: both signs *
        angular degeneracy 4(n+1)).
    """

    disk: DiscreteDiskDirac
    sphere: S2DiracSpectrum
    warp_profile: np.ndarray
    r_h: float

    def __post_init__(self) -> None:
        if len(self.warp_profile) != self.disk.N_rho:
            raise ValueError(
                f"warp_profile length {len(self.warp_profile)} must match "
                f"N_rho = {self.disk.N_rho}"
            )
        if np.any(np.asarray(self.warp_profile) <= 0):
            raise ValueError("warp_profile must be positive at every site")
        if self.r_h <= 0:
            raise ValueError(f"r_h must be > 0, got {self.r_h}")

    @classmethod
    def smooth_tip(
        cls,
        disk: DiscreteDiskDirac,
        sphere: S2DiracSpectrum,
        r_h: float,
    ) -> "VariableWarpDirac":
        """Smooth-tip warp r(rho) = r_h * sqrt(1 + (rho/r_h)^2).

        Asymptotic-Schwarzschild far field (r(rho) -> rho at large rho)
        + smooth tip at rho = 0 (no conical defect, r'(0) = 0).
        """
        rho_k = np.array(
            [(k + 1) * disk.a for k in range(disk.N_rho)]
        )
        r_profile = r_h * np.sqrt(1.0 + (rho_k / r_h) ** 2)
        return cls(
            disk=disk, sphere=sphere, warp_profile=r_profile, r_h=r_h
        )

    @classmethod
    def constant(
        cls,
        disk: DiscreteDiskDirac,
        sphere: S2DiracSpectrum,
        r_h: float,
    ) -> "VariableWarpDirac":
        """Constant warp r(rho) = r_h. For F6 Riemannian-limit testing.

        Reduces bit-exact to WarpedDiracConstant(disk, sphere) when
        sphere.r_h == r_h.
        """
        r_profile = np.full(disk.N_rho, r_h)
        return cls(
            disk=disk, sphere=sphere, warp_profile=r_profile, r_h=r_h
        )

    @property
    def rho_array(self) -> np.ndarray:
        """rho_k = k * a for k = 1, ..., N_rho."""
        return np.array(
            [(k + 1) * self.disk.a for k in range(self.disk.N_rho)]
        )

    def warp_derivative_over_warp(self) -> np.ndarray:
        """r'(rho_k) / r(rho_k) array (the spin connection magnitude).

        Computed via centered finite difference on the actual warp_profile
        so the result is correct for arbitrary warp profiles, including
        constant (gives 0) and smooth-tip.

        For smooth-tip r(rho) = r_h * sqrt(1 + (rho/r_h)^2):
            r(rho)^2 = r_h^2 + rho^2
            d(r^2)/drho = 2 r r' = 2 rho
            => r' = rho / r
            => r' / r = rho / r^2 = rho / (rho^2 + r_h^2)

        Tip behavior: rho -> 0, r'/r -> rho / r_h^2 (linear, regular)
        Asymptotic: rho -> infty, r'/r -> 1 / rho (1/rho falloff)

        For constant warp r(rho) = r_h: r' = 0 everywhere, r'/r = 0.
        """
        # Centered finite difference of warp_profile
        r = self.warp_profile
        a = self.disk.a
        N = self.disk.N_rho
        r_prime = np.zeros(N)
        # Interior: centered FD
        r_prime[1:-1] = (r[2:] - r[:-2]) / (2 * a)
        # Boundaries: one-sided FD
        if N >= 2:
            r_prime[0] = (r[1] - r[0]) / a
            r_prime[-1] = (r[-1] - r[-2]) / a
        return r_prime / r

    def H_block(
        self,
        n: int,
        k_phi_idx: int,
        include_spin_connection: bool = False,
    ) -> np.ndarray:
        """Radial Hamiltonian L_disk + diag((n+1)^2 / r(rho)^2).

        Block-diagonal-in-Fourier substructure: each (n, k_phi) pair
        gives a single N_rho x N_rho Hermitian matrix.

        Parameters
        ----------
        n : int
            S^2 mode index.
        k_phi_idx : int
            Azimuthal Fourier mode index.
        include_spin_connection : bool, default False
            If True, include the (r'/r)^2 scalar correction to D^2.
            This is the leading scalar contribution from the warped-product
            spin connection r'(rho)/r(rho) * gamma^rho (Camporesi 1996).
            The (r'/r)^2 term is n-independent and acts as a universal
            scalar shift on the radial diagonal. Full Level 2 (including
            gamma^rho mixing cross terms) is deferred to a separate
            multi-week sprint; this is "Level 1.5".
        """
        if k_phi_idx <= self.disk.N_phi // 2:
            k = k_phi_idx
        else:
            k = k_phi_idx - self.disk.N_phi
        m_eff_sq = (2.0 / self.disk.h_phi) ** 2 * np.sin(
            np.pi * (k + 0.5) / self.disk.N_phi
        ) ** 2
        m_eff = float(np.sqrt(m_eff_sq))
        L_k = self.disk._hermitian_radial_laplacian(m_eff)
        mass_diag = (n + 1) ** 2 / self.warp_profile ** 2
        if include_spin_connection:
            spin_conn_sq = self.warp_derivative_over_warp() ** 2
            return L_k + np.diag(mass_diag + spin_conn_sq)
        return L_k + np.diag(mass_diag)

    def heat_trace(
        self, t: float, include_spin_connection: bool = False
    ) -> float:
        """Tr exp(-t D_cigar^2) at variable warp.

        Direct summation over (n, k_phi) blocks with multiplicity
        16(n+1) per radial eigenvalue. Pass
        ``include_spin_connection=True`` to upgrade Level 1 to Level 1.5
        (adds the (r'/r)^2 scalar correction to D^2).
        """
        K = 0.0
        for n in range(self.sphere.l_max + 1):
            mult = 16 * (n + 1)
            for k_phi_idx in range(self.disk.N_phi):
                H = self.H_block(
                    n, k_phi_idx,
                    include_spin_connection=include_spin_connection,
                )
                evals = np.linalg.eigvalsh(H)
                K += mult * float(np.sum(np.exp(-evals * t)))
        return K


def verify_F4_tip_regular(
    disk: DiscreteDiskDirac,
    r_h: float,
    factor_tol: float = 0.50,
) -> Dict:
    """F4: warp-derivative term r'(rho)/r(rho) is regular at the apex.

    For smooth-tip warp:
      r'/r = rho / (rho^2 + r_h^2)
      Tip limit: -> rho/r_h^2 (linear, regular)

    Test:
      - all values finite (no NaN/inf)
      - smallest-rho value matches rho_1 / r_h^2 within factor_tol
    """
    sphere = S2DiracSpectrum(l_max=2, r_h=r_h)
    var = VariableWarpDirac.smooth_tip(disk=disk, sphere=sphere, r_h=r_h)
    deriv = var.warp_derivative_over_warp()
    rho_1 = disk.a
    expected_tip = rho_1 / r_h ** 2
    measured_tip = float(deriv[0])
    finite = bool(np.all(np.isfinite(deriv)))
    max_value = float(np.max(deriv))
    rel_err_tip = abs(measured_tip - expected_tip) / expected_tip
    return {
        "deriv_array_finite": finite,
        "deriv_at_rho_1": measured_tip,
        "expected_at_rho_1": expected_tip,
        "rel_err_tip": float(rel_err_tip),
        "max_deriv_value": max_value,
        "passed": bool(finite and rel_err_tip < factor_tol),
    }


def verify_F6_riemannian_limit(
    disk: DiscreteDiskDirac,
    sphere: S2DiracSpectrum,
    r_h: float,
    t_values: List[float],
    tol: float = 1e-10,
) -> Dict:
    """F6: at constant warp r(rho) = r_h, variable = constant bit-exact.

    Load-bearing falsifier: the leading-position-dependent-mass
    Level 1 approximation must reduce to G4-4a's constant-warp
    factorization bit-exactly when r(rho) is held constant.
    """
    if not np.isclose(sphere.r_h, r_h):
        raise ValueError(
            f"sphere.r_h = {sphere.r_h} must match r_h = {r_h}"
        )
    var = VariableWarpDirac.constant(disk=disk, sphere=sphere, r_h=r_h)
    const = WarpedDiracConstant(disk=disk, sphere=sphere)
    results = {}
    all_passed = True
    for t in t_values:
        K_var = var.heat_trace(t)
        K_const = const.heat_trace_factorized(t)
        denom = abs(K_const) if K_const != 0 else 1.0
        rel_err = abs(K_var - K_const) / denom
        passed = bool(rel_err < tol)
        all_passed = all_passed and passed
        results[str(t)] = {
            "K_var": float(K_var),
            "K_const_factorized": float(K_const),
            "rel_err": float(rel_err),
            "passed": passed,
        }
    results["all_passed"] = all_passed
    results["tolerance"] = tol
    return results


def verify_F7_factorization_loss(
    disk: DiscreteDiskDirac,
    sphere: S2DiracSpectrum,
    r_h: float,
    t_values: List[float],
) -> Dict:
    """F7: factorization-loss at variable warp vs constant-warp.

    At variable warp r(rho), K_cigar^var is NOT K_disk * K_S^2_r_h
    in general. The difference Delta_fact(t) = K_var(t) - K_const(t)
    is the structural signature of the variable warp.

    Predicted sign for smooth-tip warp (r(rho) > r_h for rho > 0):
        1/r(rho)^2 < 1/r_h^2 (smaller S^2 mass at non-tip sites)
        -> smaller eigenvalues -> LARGER heat trace contribution
        -> K_var(t) > K_const(t)
        -> Delta_fact > 0

    Reports K_var, K_const, Delta, ratio per t.
    """
    var = VariableWarpDirac.smooth_tip(disk=disk, sphere=sphere, r_h=r_h)
    const = WarpedDiracConstant(disk=disk, sphere=sphere)
    results = {}
    for t in t_values:
        K_var = var.heat_trace(t)
        K_const = const.heat_trace_factorized(t)
        Delta = K_var - K_const
        ratio = K_var / K_const if K_const != 0 else float("inf")
        results[str(t)] = {
            "K_var": float(K_var),
            "K_const_factorized": float(K_const),
            "Delta_fact": float(Delta),
            "ratio": float(ratio),
            "sign_positive": bool(Delta > 0),
        }
    return results


def verify_F3_continuum_recovery_rough(
    disk: DiscreteDiskDirac,
    t_values: List[float],
    tol_rank_2: float = 0.10,
) -> Dict[str, Dict]:
    """F3 rough check: K_Dirac(t) / K_scalar(t) -> rank factor 2 at small t.

    The continuum Weyl law leading order:
        K_scalar(t) ~ A_{D^2} / (4 pi t)         (rank 1, periodic phi)
        K_Dirac(t) ~ 2 A_{D^2} / (4 pi t)        (rank 2, anti-periodic phi)

    Ratio K_Dirac / K_scalar -> 2 in continuum limit as t -> 0 (UV
    regime). This is a first-pass F3 check; the full F3 with
    quantitative continuum recovery requires UV sweeps and is deferred
    to subsequent G4-4a weeks.

    Parameters
    ----------
    disk : DiscreteDiskDirac
    t_values : list of float
    tol_rank_2 : float
        Tolerance for ratio - 2. At sprint-scale N_phi = 24, ratio
        should be within ~ 0.20 of 2.0 at modest t.
    """
    scalar = DiscreteDiskScalar(N_rho=disk.N_rho, a=disk.a, N_phi=disk.N_phi)
    results: Dict[str, Dict] = {}
    for t in t_values:
        K_Dirac = disk.heat_trace(t)
        K_scalar = scalar.heat_trace(t)
        ratio = K_Dirac / K_scalar
        passed = bool(abs(ratio - 2.0) < tol_rank_2 * 2.0)
        results[str(t)] = {
            "K_Dirac": K_Dirac,
            "K_scalar": K_scalar,
            "ratio": ratio,
            "target": 2.0,
            "passed": passed,
        }
    return results
