"""Compact-temporal Lorentzian Dirac on S^3 x S^1_T at signature (3, 1).

Sprint L3b first move (2026-05-17): the compact-temporal variant of the
L2-C Lorentzian Dirac (`geovac.lorentzian_dirac`).  Replaces the
centered finite-difference d/dt on a bounded grid with the Fourier-
diagonal d/dt = i*omega on the compact circle S^1_T.

Architectural motivation
========================

The L2-C construction uses centered FD with Dirichlet zero BCs, which is
anti-Hermitian on the bounded grid but breaks periodicity.  On the
compact S^1_T, the natural choice is the Fourier-diagonal representation:

    (d/dt) e_k(t) = (2 pi i k / T) e_k(t),

so D_t is exactly diagonal in the momentum basis with eigenvalues
2 pi i k / T.  This is the SHARP discretization (no truncation error
within the kept modes) and anti-Hermitian by construction
((2 pi i k / T)^* = -2 pi i k / T).

Construction
============

The Lorentzian Dirac on S^3 x S^1_T:

    D_L = i * [ gamma^0 (x) d/dt + D_GV (x) I_{N_t} ]
        = i * [ J_spatial (x) D_t_fourier + D_GV (x) I_{N_t} ],

where:
  - J_spatial is the spatial lift of gamma^0 (chirality swap),
  - D_t_fourier is the diagonal momentum-space derivative with
    eigenvalues {2 pi i k / T : k in momentum_grid},
  - D_GV is the Camporesi-Higuchi full-Dirac matrix on H_GV.

Krein-self-adjointness D_L^x = J D_L^dagger J = D_L follows from the
SAME derivation as L2-C:
  - gamma^0^2 = +I (after lift to J_spatial^2 = +I).
  - {gamma^0, D_GV}_H_GV = 0 (chirality swap vs chirality-diagonal).
  - (D_t)^dagger = -D_t (anti-Hermitian, since D_t = i * diag(omega_k)).
  - The temporal D_t commutes with J_spatial because they act on
    different tensor factors.

The proof at the matrix level is identical to L2-C, with the only
change being the temporal D_t implementation.

Riemannian-limit check
======================

At N_t = 1 (singleton k = 0 mode), D_t = (0) (the 1x1 zero matrix,
because omega_0 = 0).  The Lorentzian Dirac reduces to

    D_L | _{N_t=1} = i * D_GV (x) I_1 = i * D_GV.

Same as L2-C: spectrum |spec(D_L)| = |spec(D_GV)|, identical
multiplicities, modulo the global i factor.

API
===

  fourier_d_dt_matrix(N_t, T)
      Diagonal momentum-space (N_t, N_t) derivative matrix:
      diag(2 pi i k / T : k in momentum_grid).

  lorentzian_dirac_compact_matrix(krein, dirac_diag=None)
      Build D_L on a CompactTemporalKreinSpace.

  verify_krein_self_adjoint, verify_anti_hermitian, etc.
      Same checks as L2-C.

References
==========

  Same as L2-C; see geovac/lorentzian_dirac.py docstring.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from geovac.full_dirac_operator_system import (
    camporesi_higuchi_full_dirac_matrix,
)
from geovac.krein_space_compact_temporal import (
    CompactTemporalKreinSpace,
    fourier_momentum_grid,
)
from geovac.lorentzian_dirac import (
    krein_adjoint,
    spatial_chirality_grading,
    verify_anti_hermitian,
    verify_krein_self_adjoint,
)


# ---------------------------------------------------------------------------
# Fourier temporal derivative
# ---------------------------------------------------------------------------


def fourier_d_dt_matrix(N_t: int, T: float = 2.0 * np.pi) -> np.ndarray:
    """Build the Fourier-diagonal d/dt matrix on S^1_T.

    The matrix is diagonal in the momentum basis with eigenvalues
    omega_k = 2 pi i k / T, where k ranges over the symmetric Fourier
    momentum grid (see `fourier_momentum_grid`).

    The matrix is ANTI-HERMITIAN by construction:
      (D_t)^dagger_{kl} = (2 pi i l / T)^* delta_{lk} = -2 pi i l / T
                        = -(D_t)_{kl}.

    At N_t = 1 the singleton k=0 gives D_t = (0).  This matches the
    Riemannian-limit recovery exactly.

    Parameters
    ----------
    N_t : int
        Number of Fourier modes (>= 1).
    T : float
        Circumference of S^1_T (> 0).  Default 2 pi.

    Returns
    -------
    D_t : np.ndarray, shape (N_t, N_t), complex128
        Diagonal anti-Hermitian matrix.
    """
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    if T <= 0:
        raise ValueError(f"T must be > 0, got {T}")
    k_grid = fourier_momentum_grid(N_t)
    omegas = 2.0 * np.pi * k_grid.astype(np.float64) / T
    diag = 1j * omegas  # anti-Hermitian, purely imaginary
    return np.diag(diag.astype(np.complex128))


# ---------------------------------------------------------------------------
# Lorentzian Dirac on S^3 x S^1_T
# ---------------------------------------------------------------------------


def lorentzian_dirac_compact_matrix(
    krein: CompactTemporalKreinSpace,
    dirac_diag: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Build D_L on the compact-temporal Krein space.

    D_L = i * [ J_spatial (x) D_t_fourier + D_GV (x) I_{N_t} ],

    where D_t_fourier = diag(2 pi i k / T) is the Fourier-diagonal
    d/dt.

    Parameters
    ----------
    krein : CompactTemporalKreinSpace
    dirac_diag : np.ndarray, optional
        Override the spatial Dirac.  Default: truthful Camporesi-Higuchi.

    Returns
    -------
    D_L : np.ndarray, shape (dim, dim), complex128
        Lorentzian Dirac, Krein-self-adjoint w.r.t. krein.J.
    """
    if dirac_diag is None:
        D_spatial = camporesi_higuchi_full_dirac_matrix(krein.basis_spatial)
    else:
        D_spatial = np.asarray(dirac_diag, dtype=np.complex128)
        if D_spatial.shape != (krein.dim_spatial, krein.dim_spatial):
            raise ValueError(
                f"dirac_diag must have shape "
                f"({krein.dim_spatial}, {krein.dim_spatial}), "
                f"got {D_spatial.shape}"
            )

    D_t = fourier_d_dt_matrix(krein.N_t, krein.T)
    I_t = np.eye(krein.N_t, dtype=np.complex128)

    # Time part: gamma^0 (x) d/dt -> J_spatial (x) D_t on K
    time_part = np.kron(krein.J_spatial, D_t)

    # Space part: D_GV (x) I_{N_t}
    space_part = np.kron(D_spatial, I_t)

    # Full Lorentzian Dirac (i^t = i for t = 1)
    D_L = 1j * (time_part + space_part)
    return D_L


# ---------------------------------------------------------------------------
# Convenience: structural verifications
# ---------------------------------------------------------------------------


def verify_riemannian_limit_compact(
    n_max: int, tol: float = 1e-12
) -> Tuple[bool, dict]:
    """LOAD-BEARING: at N_t = 1, D_L_compact = i * D_GV bit-identically.

    Same load-bearing check as L2-C, but for the compact-temporal
    construction.  Must pass bit-exact (residual = 0.0) for the
    foundation to be sound.

    Parameters
    ----------
    n_max : int
    tol : float
        Frobenius tolerance.

    Returns
    -------
    (ok, details_dict) with keys 'residual_F_norm',
    'spectral_match_residual'.
    """
    krein_compact = CompactTemporalKreinSpace(n_max=n_max, N_t=1)
    D_L = lorentzian_dirac_compact_matrix(krein_compact)
    D_GV = camporesi_higuchi_full_dirac_matrix(krein_compact.basis_spatial)

    # D_L at N_t=1 should equal i * D_GV
    residual = float(np.linalg.norm(D_L - 1j * D_GV))

    # Spectrum check: |spec(D_L)| should match |spec(D_GV)|
    spec_D_L = np.linalg.eigvalsh(D_L.conj().T @ D_L)
    spec_D_GV_sq = np.linalg.eigvalsh(D_GV.conj().T @ D_GV)
    spec_diff = float(np.max(np.abs(np.sort(spec_D_L) - np.sort(spec_D_GV_sq))))

    details = {
        "n_max": n_max,
        "residual_F_norm": residual,
        "spectral_match_residual": spec_diff,
    }
    ok = residual < tol and spec_diff < tol
    return ok, details


def krein_self_adjoint_residual(
    krein: CompactTemporalKreinSpace,
    dirac_diag: Optional[np.ndarray] = None,
) -> float:
    """Return the Krein-self-adjointness residual ||D_L^x - D_L||_F.

    For truthful Camporesi-Higuchi spatial Dirac, this should be
    machine-precision zero by the algebraic identity proved in L2-C.
    """
    D_L = lorentzian_dirac_compact_matrix(krein, dirac_diag=dirac_diag)
    D_L_krein = krein_adjoint(D_L, krein.J)
    return float(np.linalg.norm(D_L_krein - D_L))
