"""Compact-temporal Krein space on S^3 x S^1_T at signature (3, 1).

Sprint L3b first move (2026-05-17): the compact-temporal variant of the
L2-B Krein space (`geovac.krein_space_construction.KreinSpace`).
Replaces the bounded uniform grid on R_t with a Fourier basis on the
periodic circle S^1_T of circumference T.

Architectural motivation (debug/l3_scoping_memo.md path-of-least-resistance)
---------------------------------------------------------------------------

The strong-form L3 deliverable (continuum Lorentzian propinquity on
S^3 x R) is blocked by the joint compact x non-compact Fejer kernel
question (Q1 of the scoping memo, 2-4 months of original NCG-math work).
The path-of-least-resistance is to first compactify the temporal
direction to S^1_T with periodic BCs, getting a compact x compact
product where both factors admit standard Peter-Weyl / Fourier
machinery.  The de-compactification T -> infinity is a SEPARATE limit
(Sprint L3c), not handled here.

WARNING (CTC obstruction, L2-A audit Section 3.8): periodic temporal
boundary conditions on a Lorentzian manifold normally create closed
timelike curves (Geroch's theorem on causality and topology).  For
the propinquity-FOUNDATION purpose of this sprint we accept this:
the compact-temporal manifold S^3 x S^1_T is a Euclidean / Wick-
rotated geometry (NOT a physical Lorentzian spacetime).  The Krein
structure J = gamma^0 is preserved as an algebraic object, but the
spacetime interpretation shifts: this is the "Wick-rotated propinquity"
candidate workaround (i) of the scoping memo Section 2 L5 entry.

Mathematical setup
==================

The compact-temporal Krein space is

    K_{n_max, N_t, T} = H_GV^{n_max} (x) C^{N_t},

where C^{N_t} is the truncation of L^2(S^1_T) to the first N_t Fourier
modes (the Dirichlet-style truncation of the temporal Hilbert space).
The Fourier basis is

    e_k(t) = exp(2 pi i k t / T) / sqrt(T),   k = -N_t//2, ..., N_t//2,

evaluated at the discrete momentum index k.  The truncation keeps
|k| <= K_max = (N_t - 1) // 2 for odd N_t, or |k| <= K_max with a
single Nyquist mode for even N_t.

The fundamental symmetry is the SAME as in L2-B:

    J = gamma^0_spatial (x) I_{N_t},

with gamma^0_spatial the chirality-swap on H_GV inherited from
`geovac.krein_space_construction.spatial_fundamental_symmetry`.  The
temporal slot is the identity I_{N_t} in the Fourier basis (the
chirality swap is purely spatial, independent of momentum basis).

This matches L2-B verbatim at N_t = 1: the singleton Fourier mode k=0
gives I_1, and the Krein space reduces to H_GV bit-identically.

Conventions
===========

- T > 0 is the circumference of S^1_T.  T = 2 pi is canonical (matches
  the modular period 2 pi from L2-E's BW-alpha construction).
- N_t is the number of Fourier modes kept.  For symmetric truncation
  on odd N_t we keep k in {-K_max, ..., 0, ..., +K_max} with
  K_max = (N_t - 1) / 2.  For even N_t we keep an asymmetric range.
- The Krein space carries no continuous-time structure; the temporal
  multipliers act diagonally in momentum.

Riemannian-limit check
======================

At N_t = 1 (singleton k=0 mode), C^1 is trivial and the Krein space
reduces to H_GV bit-identically.  J | _{N_t=1} = J_spatial (bit-exact).

API
===

  fourier_momentum_grid(N_t)
      Symmetric integer momentum grid k in {-K_max, ..., +K_max}.

  CompactTemporalKreinSpace(n_max, N_t=11, T=2pi)
      .dim                  : N_Dirac(n_max) * N_t
      .J                     : fundamental symmetry (dim, dim)
      .momentum_grid         : the Fourier momentum index array
      .k_to_omega(k)         : k -> 2 pi k / T (continuous frequency)
      .riemannian_limit_check()
      .positive_negative_split()

References
==========

  van den Dungen, K. "Krein spectral triples and the fermionic action."
    Math. Phys. Anal. Geom. 19, 4 (2016). arXiv:1505.01939.  Prop 4.1.

  GeoVac internal:
    geovac/krein_space_construction.py  L2-B substrate
    debug/l3_scoping_memo.md            path-of-least-resistance
    papers/standalone/paper_43_lorentzian_extension.tex  L2 closure
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Tuple

import numpy as np

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    full_dirac_basis,
    full_dirac_dim,
)
from geovac.krein_space_construction import spatial_fundamental_symmetry


# ---------------------------------------------------------------------------
# Symmetric Fourier momentum grid on S^1_T
# ---------------------------------------------------------------------------


def fourier_momentum_grid(N_t: int) -> np.ndarray:
    """Symmetric integer momentum grid for the S^1_T Fourier truncation.

    For odd N_t: returns k in {-K_max, ..., 0, ..., +K_max} with
    K_max = (N_t - 1) // 2.  For even N_t: returns an asymmetric range
    {-K_max, ..., 0, ..., K_max - 1} with K_max = N_t // 2 (the Nyquist
    mode lives on the +k side by convention).

    At N_t = 1: returns the singleton [0] (the constant Fourier mode,
    Riemannian-limit reduction).

    Parameters
    ----------
    N_t : int
        Number of Fourier modes kept.  Must be >= 1.

    Returns
    -------
    k_grid : np.ndarray of shape (N_t,), int64
    """
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    if N_t == 1:
        return np.array([0], dtype=np.int64)
    if N_t % 2 == 1:
        K_max = (N_t - 1) // 2
        return np.arange(-K_max, K_max + 1, dtype=np.int64)
    else:
        K_max = N_t // 2
        return np.arange(-K_max, K_max, dtype=np.int64)


# ---------------------------------------------------------------------------
# CompactTemporalKreinSpace
# ---------------------------------------------------------------------------


@dataclass
class CompactTemporalKreinSpace:
    """Compact-temporal Krein space K = H_GV (x) C^{N_t} on S^3 x S^1_T.

    Same algebraic structure as `geovac.krein_space_construction.KreinSpace`,
    but with the temporal slot interpreted as the first N_t Fourier
    modes on the periodic circle S^1_T of circumference T.

    Parameters
    ----------
    n_max : int
        Spatial spinor-bundle truncation.  Must be >= 1.
    N_t : int
        Number of Fourier modes kept on S^1_T.  Default 11.  Must be >= 1.
    T : float
        Circumference of S^1_T.  Default 2 pi (canonical, matches the
        BW-alpha modular period).  Must be > 0.

    Attributes
    ----------
    basis_spatial : list of FullDiracLabel
    dim_spatial : int
    momentum_grid : np.ndarray, shape (N_t,)
        Integer momentum indices.
    dim : int
        Total: dim_spatial * N_t.
    J_spatial : np.ndarray, shape (dim_spatial, dim_spatial)
    J : np.ndarray, shape (dim, dim)
        Full fundamental symmetry, = J_spatial (x) I_{N_t}.
    """

    n_max: int
    N_t: int = 11
    T: float = 2.0 * np.pi

    basis_spatial: list[FullDiracLabel] = field(init=False)
    dim_spatial: int = field(init=False)
    momentum_grid: np.ndarray = field(init=False)
    dim: int = field(init=False)
    J_spatial: np.ndarray = field(init=False)
    J: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        if self.n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {self.n_max}")
        if self.N_t < 1:
            raise ValueError(f"N_t must be >= 1, got {self.N_t}")
        if self.T <= 0:
            raise ValueError(f"T must be > 0, got {self.T}")

        # Spatial: full Camporesi-Higuchi spinor space (same as L2-B)
        self.basis_spatial = full_dirac_basis(self.n_max)
        self.dim_spatial = len(self.basis_spatial)
        assert self.dim_spatial == full_dirac_dim(self.n_max)

        # Temporal: symmetric Fourier momentum grid
        self.momentum_grid = fourier_momentum_grid(self.N_t)

        # Total dimension
        self.dim = self.dim_spatial * self.N_t

        # Spatial fundamental symmetry (lift of gamma^0 in chiral basis)
        self.J_spatial = spatial_fundamental_symmetry(self.basis_spatial)

        # Full J = J_spatial (x) I_{N_t}
        I_t = np.eye(self.N_t, dtype=np.complex128)
        self.J = np.kron(self.J_spatial, I_t)

    # -------- Frequency accessor --------

    def k_to_omega(self, k: int) -> float:
        """Map momentum index k to continuous frequency omega = 2 pi k / T."""
        return 2.0 * np.pi * float(k) / self.T

    def omegas(self) -> np.ndarray:
        """Vector of continuous frequencies omega_k = 2 pi k / T."""
        return 2.0 * np.pi * self.momentum_grid.astype(np.float64) / self.T

    # -------- Krein axiom checks (identical structure to L2-B) --------

    def verify_J_squared_identity(
        self, tol: float = 1e-12
    ) -> Tuple[bool, float]:
        """J^2 = +I."""
        I = np.eye(self.dim, dtype=np.complex128)
        residual = float(np.linalg.norm(self.J @ self.J - I))
        return residual < tol, residual

    def verify_J_unitary(self, tol: float = 1e-12) -> Tuple[bool, float]:
        """J^* J = I."""
        I = np.eye(self.dim, dtype=np.complex128)
        residual = float(np.linalg.norm(self.J.conj().T @ self.J - I))
        return residual < tol, residual

    def verify_J_hermitian(self, tol: float = 1e-12) -> Tuple[bool, float]:
        """J = J^*."""
        residual = float(np.linalg.norm(self.J - self.J.conj().T))
        return residual < tol, residual

    # -------- Krein inner product --------

    def krein_inner_product(
        self, psi: np.ndarray, phi: np.ndarray
    ) -> complex:
        if psi.shape != (self.dim,) or phi.shape != (self.dim,):
            raise ValueError(
                f"shapes must be ({self.dim},); got {psi.shape}, {phi.shape}"
            )
        return complex(np.vdot(psi, self.J @ phi))

    # -------- K^+ / K^- split --------

    def positive_negative_split(
        self, tol: float = 1e-10
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Projectors P_plus = (I + J)/2, P_minus = (I - J)/2."""
        ok, residual = self.verify_J_squared_identity(tol=tol)
        if not ok:
            raise ValueError(
                f"J^2 != I (residual {residual:.3e}); cannot split"
            )
        I = np.eye(self.dim, dtype=np.complex128)
        P_plus = 0.5 * (I + self.J)
        P_minus = 0.5 * (I - self.J)
        return P_plus, P_minus

    def split_dimensions(self) -> Tuple[int, int]:
        """Return (dim K^+, dim K^-) via trace of projectors."""
        P_plus, P_minus = self.positive_negative_split()
        d_plus = int(np.round(np.real(np.trace(P_plus))))
        d_minus = int(np.round(np.real(np.trace(P_minus))))
        return d_plus, d_minus

    # -------- Riemannian-limit check (LOAD-BEARING) --------

    def riemannian_limit_check(
        self, tol: float = 1e-14
    ) -> Tuple[bool, dict]:
        """At N_t = 1, K reduces to H_GV bit-identically.

        Constructs a fresh CompactTemporalKreinSpace at N_t = 1 (same
        n_max) and verifies basis-by-basis match against
        `full_dirac_basis(n_max)` and `spatial_fundamental_symmetry`.

        Returns
        -------
        (ok, details_dict)
        """
        K_rie = CompactTemporalKreinSpace(n_max=self.n_max, N_t=1, T=self.T)

        ref_basis = full_dirac_basis(self.n_max)
        ref_dim = full_dirac_dim(self.n_max)
        ref_J = spatial_fundamental_symmetry(ref_basis)

        dim_match = K_rie.dim == ref_dim
        basis_match = all(
            K_rie.basis_spatial[i] == ref_basis[i] for i in range(ref_dim)
        )
        J_residual = float(np.linalg.norm(K_rie.J - ref_J))

        # The current self at N_t > 1 has J = J_spatial (x) I_{N_t};
        # the (0, 0) block in temporal momentum should equal J_spatial.
        time_slice_res = 0.0
        if self.N_t >= 1:
            slice_k0 = self.J[0 :: self.N_t, 0 :: self.N_t]
            time_slice_res = float(np.linalg.norm(slice_k0 - ref_J))

        details = {
            "n_max": self.n_max,
            "N_t": self.N_t,
            "T": self.T,
            "dim_match": dim_match,
            "basis_match": basis_match,
            "J_match_residual": J_residual,
            "time_slice_check_residual": time_slice_res,
            "ref_dim": ref_dim,
            "compact_rie_dim": K_rie.dim,
        }
        ok = (
            dim_match
            and basis_match
            and J_residual < tol
            and time_slice_res < tol
        )
        return ok, details

    # -------- Compatibility check vs L2-B (matching algebraic structure) --------

    def matches_l2b_at_N_t(self, l2b_krein) -> Tuple[bool, dict]:
        """Verify J matches the L2-B KreinSpace J for same n_max, N_t.

        Both constructions use J = J_spatial (x) I_{N_t} where J_spatial
        is from spatial_fundamental_symmetry.  At equal n_max and N_t,
        the J matrices should be bit-identical regardless of whether the
        temporal slot is interpreted as grid (L2-B) or momentum (L3b).

        Parameters
        ----------
        l2b_krein : geovac.krein_space_construction.KreinSpace
            An L2-B KreinSpace instance at matching n_max, N_t.

        Returns
        -------
        (ok, details)
        """
        if l2b_krein.n_max != self.n_max or l2b_krein.N_t != self.N_t:
            return False, {
                "n_max_match": l2b_krein.n_max == self.n_max,
                "N_t_match": l2b_krein.N_t == self.N_t,
            }
        residual = float(np.linalg.norm(self.J - l2b_krein.J))
        return residual < 1e-14, {"J_residual": residual}

    def __repr__(self) -> str:
        return (
            f"CompactTemporalKreinSpace(n_max={self.n_max}, "
            f"N_t={self.N_t}, T={self.T:.4f}, "
            f"dim_spatial={self.dim_spatial}, dim={self.dim})"
        )
