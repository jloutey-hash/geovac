"""Compact-temporal Lorentzian truncated operator system on S^3 x S^1_T.

Sprint L3b first move (2026-05-17): compact-temporal analog of the L3a-1
operator system (`geovac.operator_system_lorentzian`).  Replaces the
bounded uniform temporal grid with the Fourier-mode truncation on the
periodic circle S^1_T, building the truncated operator system on the
compact-temporal Krein space (`CompactTemporalKreinSpace`).

Architectural motivation
========================

The L3a-1 substrate works on K_{n_max, N_t} = H_GV (x) C^{N_t} with the
temporal slot interpreted as a bounded grid (point-evaluation algebra).
L3b lifts this to the compact-temporal setting where the temporal slot
is the Fourier momentum basis on S^1_T.  Both interpretations carry the
same fundamental symmetry J = J_spatial (x) I_{N_t}; the algebra acting
on the temporal slot is what changes.

In the L3a-1 grid setting, the natural multiplier basis was the
polynomial g_p(t) = t^p evaluated pointwise; the algebra is the
commutative *-algebra of diagonal matrices on C^{N_t}.  In the
compact-temporal Fourier setting, the natural multiplier basis is the
Fourier mode e_q(theta) = exp(i q theta) acting as the shift / cyclic
matrix in momentum space (Toeplitz-like), of cyclic structure on C^{N_t}.
However, for the propinquity-foundation purpose of this sprint we keep
the multiplier algebra COMMUTATIVE: we represent the Fourier multiplier
e_q as the action of multiplication by e^{i q theta} in the POSITION
basis (via inverse Fourier transform), which when truncated to the kept
modes acts as a banded matrix in momentum space.

To preserve the bit-exact Riemannian-limit reduction at N_t = 1 and to
match the load-bearing falsifier pattern of L3a-1, we adopt the
following commutative-temporal-algebra convention:

  - Temporal multipliers are diagonal matrices in the MOMENTUM basis
    with diagonal pattern g_p(omega_k) = (omega_k)^p, p = 0, ..., N_t-1,
    where omega_k = 2 pi k / T are the continuous frequencies.

This is the EXACT analog of the L3a-1 grid construction, with omega_k
replacing the grid points t_k.  The Vandermonde structure on N_t distinct
omega_k values guarantees the polynomial basis spans the full N_t-
dimensional commutative diagonal subalgebra.

At N_t = 1 (singleton omega_0 = 0): only the constant multiplier g_0 = 1
is non-degenerate (g_p(0) = 0 for p > 0, but we still keep the basis
formally; the linear-span dim collapses to 1).  The Riemannian-limit
reduction matches L3a-1 verbatim.

The Fourier-mode-shift alternative (non-commutative temporal algebra) is
a candidate Sprint L3b-2 extension; we do not implement it here.

Note: at the propinquity FOUNDATION level we only need a well-defined
operator-system construction, not the propinquity-bound computation
itself.  The two conventions (grid-polynomial / momentum-polynomial)
give isomorphic constructions for the purposes of L3b foundation.

API
===

  CompactTemporalTruncatedOperatorSystem(n_max, N_t=11, T=2*pi)
      .multiplier_labels : list of (N, L, M, p)
      .multiplier_matrices : list of (dim_K, dim_K) complex arrays
      .dim                : linear-span dim over C
      .envelope_dim       : dim_K**2
      .achievable_envelope_dim : dim_Weyl^2 * N_t
      .contains(M)        : least-squares membership test
      .identity_in_O()    : verifies I in O
      .is_star_closed()
      .compute_propagation_number(envelope='achievable' or 'full')
      .verify_riemannian_limit()
      .compare_to_l3a1(O_grid)  : compare structure to grid-version
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    FullDiracTruncatedOperatorSystem,
    full_dirac_dim,
)
from geovac.krein_space_compact_temporal import (
    CompactTemporalKreinSpace,
    fourier_momentum_grid,
)
from geovac.operator_system import allowed_multiplier_labels
from geovac.spinor_operator_system import build_spinor_multiplier_matrix


# ---------------------------------------------------------------------------
# Temporal multiplier basis on S^1_T (Fourier-momentum-diagonal convention)
# ---------------------------------------------------------------------------


def compact_temporal_multiplier_matrices(
    N_t: int, T: float = 2.0 * np.pi
) -> List[np.ndarray]:
    """Polynomial-in-frequency basis for the compact-temporal multiplier algebra.

    Returns N_t diagonal (N_t, N_t) matrices

        diag(omega_{k_0}^p, omega_{k_1}^p, ..., omega_{k_{N_t-1}}^p),
        p = 0, ..., N_t - 1,

    where {omega_k} = {2 pi k / T : k in momentum_grid}.

    At p = 0 the matrix is the identity I_{N_t}.

    Convention: this is the momentum-space analog of the L3a-1 grid
    polynomial basis.  The Vandermonde matrix on N_t distinct omega_k
    (which holds when N_t > 1) is invertible, so the basis spans the
    full N_t-dim commutative diagonal subalgebra.

    At N_t = 1 the singleton omega_0 = 0 gives g_p(0) = delta_{p,0}; we
    keep all N_t formal polynomial basis elements, but only g_0 = I_1
    is non-degenerate.  For consistency with L3a-1 we return exactly
    N_t matrices.

    Parameters
    ----------
    N_t : int
    T : float

    Returns
    -------
    list of np.ndarray, each shape (N_t, N_t), dtype complex128.
    """
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    if T <= 0:
        raise ValueError(f"T must be > 0, got {T}")
    k_grid = fourier_momentum_grid(N_t)
    omegas = 2.0 * np.pi * k_grid.astype(np.float64) / T
    mats = []
    for p in range(N_t):
        if p == 0:
            diag = np.ones(N_t, dtype=np.complex128)
        else:
            diag = (omegas ** p).astype(np.complex128)
        mats.append(np.diag(diag))
    return mats


# ---------------------------------------------------------------------------
# CompactTemporalTruncatedOperatorSystem
# ---------------------------------------------------------------------------


class CompactTemporalTruncatedOperatorSystem:
    """Truncated operator system on the compact-temporal Krein space.

    Construction
    ------------
    Multipliers are tensor products of:
      - Spatial chirality-doubled scalar 3-Y multipliers M^spat_{N,L,M}
        (Avery-Wen-Avery via the spinor lift, identical to
        FullDiracTruncatedOperatorSystem).
      - Temporal momentum-polynomial multipliers g_p(omega_k) = omega_k^p,
        p = 0..N_t-1.

    The tensor-product construction at the matrix level uses np.kron;
    O = span_C {M^spat (x) g_p}.

    Parameters
    ----------
    n_max : int
    N_t : int
    T : float
    rank_tol : float

    Attributes
    ----------
    krein : CompactTemporalKreinSpace
    basis_spatial : list of FullDiracLabel
    dim_spatial : int
    dim_K : int
    momentum_grid : ndarray
    spat_labels : list of (N, L, M)
    multiplier_labels : list of (N, L, M, p)
    multiplier_matrices : list of complex ndarrays of shape (dim_K, dim_K)
    basis_matrices : list of (label, matrix)
    dim : int
        Linear-span dim of O over C as a subspace of M_{dim_K}(C).
    envelope_dim : int
        dim_K ** 2.
    achievable_envelope_dim : int
        dim_Weyl ** 2 * N_t  (= scalar-multiplier-reachable subspace dim).
    """

    def __init__(
        self,
        n_max: int,
        N_t: int = 11,
        T: float = 2.0 * np.pi,
        *,
        rank_tol: float = 1e-12,
    ) -> None:
        if n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {n_max}")
        if N_t < 1:
            raise ValueError(f"N_t must be >= 1, got {N_t}")
        if T <= 0:
            raise ValueError(f"T must be > 0, got {T}")
        self.n_max = n_max
        self.N_t = N_t
        self.T = T
        self._rank_tol = rank_tol

        # Build compact-temporal Krein space
        self.krein = CompactTemporalKreinSpace(n_max=n_max, N_t=N_t, T=T)
        self.basis_spatial = self.krein.basis_spatial
        self.dim_spatial = self.krein.dim_spatial
        self.dim_K = self.krein.dim
        self.momentum_grid = self.krein.momentum_grid

        # Build spatial multipliers (chirality-doubled scalar 3-Y),
        # same as FullDiracTruncatedOperatorSystem.
        all_spat_candidates = allowed_multiplier_labels(n_max)
        weyl_basis_in_order = [
            b.to_spinor() for b in self.basis_spatial[: self.dim_spatial // 2]
        ]
        self._weyl_basis_in_order = weyl_basis_in_order

        spat_labels: List[Tuple[int, int, int]] = []
        spat_matrices: List[np.ndarray] = []
        for (N, L, M) in all_spat_candidates:
            weyl_mat = build_spinor_multiplier_matrix(
                N, L, M, weyl_basis_in_order
            )
            if np.allclose(weyl_mat, 0, atol=1e-15):
                continue
            full_spat = np.zeros(
                (self.dim_spatial, self.dim_spatial), dtype=np.complex128
            )
            d_w = self.dim_spatial // 2
            full_spat[:d_w, :d_w] = weyl_mat
            full_spat[d_w:, d_w:] = weyl_mat
            spat_labels.append((N, L, M))
            spat_matrices.append(full_spat)
        self.spat_labels = spat_labels
        self._spat_matrices = spat_matrices

        # Temporal multipliers (momentum-polynomial diagonal matrices)
        temp_mats = compact_temporal_multiplier_matrices(N_t, T)
        self._temp_matrices = temp_mats

        # Tensor-product multipliers
        labels: List[Tuple[int, int, int, int]] = []
        matrices: List[np.ndarray] = []
        for (N, L, M), s_mat in zip(spat_labels, spat_matrices):
            for p, t_mat in enumerate(temp_mats):
                labels.append((N, L, M, p))
                matrices.append(np.kron(s_mat, t_mat))
        self.multiplier_labels = labels
        self.multiplier_matrices = matrices
        self.basis_matrices = list(zip(labels, matrices))

        # Cache vec-stack for span / rank computations
        self._vec_cache = self._build_vec_matrix(matrices)
        if len(matrices) == 0:
            self._dim_cache = 0
        else:
            self._dim_cache = int(
                np.linalg.matrix_rank(self._vec_cache, tol=rank_tol)
            )

    @property
    def dim(self) -> int:
        return self._dim_cache

    @property
    def envelope_dim(self) -> int:
        return self.dim_K ** 2

    @property
    def achievable_envelope_dim(self) -> int:
        """Scalar-multiplier achievable subspace dim = dim_Weyl^2 * N_t."""
        dim_Weyl = self.dim_spatial // 2
        return dim_Weyl * dim_Weyl * self.N_t

    # -----------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------

    @staticmethod
    def _build_vec_matrix(matrices: Sequence[np.ndarray]) -> np.ndarray:
        if len(matrices) == 0:
            return np.zeros((0, 0), dtype=np.complex128)
        dim = matrices[0].shape[0]
        K = len(matrices)
        cols = np.zeros((dim * dim, K), dtype=np.complex128)
        for k, M in enumerate(matrices):
            cols[:, k] = M.reshape(-1)
        return cols

    # -----------------------------------------------------------------
    # Membership test
    # -----------------------------------------------------------------

    def contains(
        self, matrix: np.ndarray, *, tol: float = 1e-10
    ) -> Tuple[bool, float]:
        """Least-squares membership test."""
        if matrix.shape != (self.dim_K, self.dim_K):
            raise ValueError(
                f"matrix shape {matrix.shape} != ({self.dim_K}, {self.dim_K})"
            )
        target = matrix.reshape(-1).astype(np.complex128)
        norm_target = np.linalg.norm(target)
        if norm_target < 1e-30:
            return True, 0.0
        V = self._vec_cache
        if V.size == 0:
            return False, float(norm_target / max(norm_target, 1.0))
        c, *_ = np.linalg.lstsq(V, target, rcond=None)
        residual = np.linalg.norm(V @ c - target)
        ratio = float(residual / norm_target)
        return ratio < tol, ratio

    # -----------------------------------------------------------------
    # Self-checks
    # -----------------------------------------------------------------

    def identity_in_O(self, *, tol: float = 1e-10) -> Tuple[bool, float]:
        """Verify I in O.  The constant function on S^3 x S^1 is in O (p=0)."""
        I = np.eye(self.dim_K, dtype=np.complex128)
        return self.contains(I, tol=tol)

    def is_star_closed(
        self, *, tol: float = 1e-10
    ) -> Tuple[bool, List[Tuple[int, float]]]:
        """Check that each generator's conjugate transpose lies in O."""
        failures = []
        for i, M in enumerate(self.multiplier_matrices):
            ok, residual = self.contains(M.conj().T, tol=tol)
            if not ok:
                failures.append((i, residual))
        return (len(failures) == 0, failures)

    # -----------------------------------------------------------------
    # Krein-positive restriction
    # -----------------------------------------------------------------

    def krein_positive_preservers(
        self, *, tol: float = 1e-10
    ) -> Tuple[List[int], List[Tuple[int, int, int, int]]]:
        """Indices and labels of multipliers preserving K^+.

        Same mechanism as L3a-1: scalar multipliers M ⊕ M commute with
        J_spatial = chirality-swap; temporal diagonal multipliers commute
        with I_{N_t}.  So every multiplier preserves K^+.

        Returns the FULL list (every generator passes).  This is the
        structural finding of L3a-1: the Krein-positive restriction is
        trivial in the natural construction.
        """
        J = self.krein.J
        preserving = []
        preserving_labels = []
        for i, (label, M) in enumerate(self.basis_matrices):
            comm = J @ M - M @ J
            residual = float(np.linalg.norm(comm))
            if residual < tol:
                preserving.append(i)
                preserving_labels.append(label)
        return preserving, preserving_labels

    # -----------------------------------------------------------------
    # Propagation number
    # -----------------------------------------------------------------

    def compute_propagation_number(
        self,
        *,
        max_k: int = 4,
        tol: float = 1e-10,
        verbose: bool = False,
        envelope: str = "achievable",
    ) -> Tuple[int, List[int]]:
        """Compute prop(O) = smallest k such that dim(O^k) = target_dim.

        Parameters
        ----------
        max_k : int
            Iteration cap (default 4 -- prop=2 should land at k=2).
        tol : float
        verbose : bool
        envelope : {'achievable', 'full'}
            Target subspace for the propagation question.

        Returns
        -------
        (prop, dim_sequence). prop = -1 if not reached within max_k.
        """
        N = self.dim_K
        if envelope == "achievable":
            target_dim = self.achievable_envelope_dim
        elif envelope == "full":
            target_dim = N * N
        else:
            raise ValueError(
                f"envelope must be 'achievable' or 'full', got {envelope!r}"
            )
        matrices_O = self.multiplier_matrices

        if not matrices_O:
            return -1, []

        dim_sequence: List[int] = []

        # k = 1
        dim_O1 = self._operator_system_dim(matrices_O, tol=tol)
        dim_sequence.append(dim_O1)
        if verbose:
            print(
                f"  k=1: |gens|={len(matrices_O)}, dim(O^1)={dim_O1}, "
                f"target={target_dim}"
            )
        if dim_O1 == target_dim:
            return 1, dim_sequence

        current_generators = list(matrices_O)
        for k in range(2, max_k + 1):
            basis_Ok_minus_1 = self._extract_matrix_basis(
                current_generators, tol=tol
            )
            new_generators = []
            for A in basis_Ok_minus_1:
                for B in matrices_O:
                    new_generators.append(A @ B)
            dim_Ok = self._operator_system_dim(new_generators, tol=tol)
            dim_sequence.append(dim_Ok)
            if verbose:
                print(
                    f"  k={k}: |basis(O^{k-1})|={len(basis_Ok_minus_1)}, "
                    f"|gens|={len(new_generators)}, "
                    f"dim(O^{k})={dim_Ok}, target={target_dim}"
                )
            if dim_Ok == target_dim:
                return k, dim_sequence
            if dim_Ok == dim_sequence[-2]:
                return -1, dim_sequence
            current_generators = new_generators

        return -1, dim_sequence

    @staticmethod
    def _operator_system_dim(
        matrices: Sequence[np.ndarray], *, tol: float = 1e-10
    ) -> int:
        if len(matrices) == 0:
            return 0
        dim = matrices[0].shape[0]
        K = len(matrices)
        cols = np.zeros((dim * dim, K), dtype=np.complex128)
        for k, M in enumerate(matrices):
            cols[:, k] = M.reshape(-1)
        return int(np.linalg.matrix_rank(cols, tol=tol))

    @staticmethod
    def _extract_matrix_basis(
        matrices: Sequence[np.ndarray], *, tol: float = 1e-10
    ) -> List[np.ndarray]:
        if len(matrices) == 0:
            return []
        dim = matrices[0].shape[0]
        K = len(matrices)
        cols = np.zeros((dim * dim, K), dtype=np.complex128)
        for i, M in enumerate(matrices):
            cols[:, i] = M.reshape(-1)
        U, S, _ = np.linalg.svd(cols, full_matrices=False)
        if S.size == 0:
            return []
        rank = int(np.sum(S > tol * max(S.max(), 1.0)))
        basis_vecs = U[:, :rank]
        return [basis_vecs[:, k].reshape((dim, dim)) for k in range(rank)]

    # -----------------------------------------------------------------
    # LOAD-BEARING: Riemannian limit at N_t = 1
    # -----------------------------------------------------------------

    def verify_riemannian_limit(
        self, *, tol_struct: float = 0.0
    ) -> Tuple[bool, dict]:
        """At N_t = 1, the construction must reduce bit-identically to
        FullDiracTruncatedOperatorSystem(n_max).

        Same load-bearing structure as the L3a-1 grid version: at N_t=1,
        the momentum-polynomial basis collapses to the singleton p=0 with
        omega_0 = 0; only g_0 = I_1 is non-degenerate.  The spatial
        construction is identical to FullDirac, and the tensor product
        with I_1 is bit-identity.

        Returns
        -------
        (ok, details)
        """
        if self.N_t != 1:
            sys_1 = CompactTemporalTruncatedOperatorSystem(
                n_max=self.n_max, N_t=1, T=self.T
            )
            return sys_1.verify_riemannian_limit(tol_struct=tol_struct)

        ref = FullDiracTruncatedOperatorSystem(self.n_max)
        ref_labels_to_idx = {
            lbl: i for i, lbl in enumerate(ref.multiplier_labels)
        }

        max_residual = 0.0
        n_compared = 0
        n_unmatched = 0
        for (N, L, M, p), self_mat in zip(
            self.multiplier_labels, self.multiplier_matrices
        ):
            if p != 0:
                continue
            if (N, L, M) not in ref_labels_to_idx:
                n_unmatched += 1
                continue
            ref_idx = ref_labels_to_idx[(N, L, M)]
            ref_mat = ref.multiplier_matrices[ref_idx]
            residual = float(np.linalg.norm(self_mat - ref_mat))
            max_residual = max(max_residual, residual)
            n_compared += 1

        details = {
            "n_max": self.n_max,
            "N_t": self.N_t,
            "dim_K": self.dim_K,
            "full_dirac_dim": full_dirac_dim(self.n_max),
            "n_multipliers_self": len(self.multiplier_labels),
            "n_multipliers_ref": len(ref.multiplier_labels),
            "n_compared": n_compared,
            "n_unmatched": n_unmatched,
            "max_residual": max_residual,
            "dim_self": self.dim,
            "dim_ref": ref.dim,
        }

        dim_ok = self.dim_K == full_dirac_dim(self.n_max)
        n_mult_ok = (
            len(self.multiplier_labels) == len(ref.multiplier_labels)
            and n_unmatched == 0
        )
        residual_ok = max_residual <= max(tol_struct, 1e-14)
        dim_O_ok = self.dim == ref.dim

        ok = dim_ok and n_mult_ok and residual_ok and dim_O_ok
        details["load_bearing_pass"] = ok
        details["dim_match"] = dim_ok
        details["multiplier_count_match"] = n_mult_ok
        details["matrix_residual_pass"] = residual_ok
        details["dim_O_match"] = dim_O_ok
        return ok, details

    # -----------------------------------------------------------------
    # Compatibility check vs L3a-1 grid version
    # -----------------------------------------------------------------

    def compare_to_l3a1_grid(self, op_grid) -> dict:
        """Compare structure to the L3a-1 grid-based operator system.

        The two constructions share spatial multipliers (chirality-doubled
        Avery-Wen-Avery 3-Y) and differ only in the temporal basis (grid
        vs momentum).  They have the same number of generators and the
        same linear-span dim over C of the spatial sub-factor.

        Parameters
        ----------
        op_grid : LorentzianTruncatedOperatorSystem
            From geovac.operator_system_lorentzian at matching (n_max, N_t).

        Returns
        -------
        dict with comparison data.
        """
        return {
            "n_max_match": op_grid.n_max == self.n_max,
            "N_t_match": op_grid.N_t == self.N_t,
            "dim_K_match": op_grid.dim_K == self.dim_K,
            "n_multipliers_match": (
                len(op_grid.multiplier_labels)
                == len(self.multiplier_labels)
            ),
            "spat_labels_match": op_grid.spat_labels == self.spat_labels,
            "dim_O_grid": op_grid.dim,
            "dim_O_compact": self.dim,
        }

    def __repr__(self) -> str:
        return (
            f"CompactTemporalTruncatedOperatorSystem(n_max={self.n_max}, "
            f"N_t={self.N_t}, T={self.T:.4f}, dim_K={self.dim_K}, "
            f"dim={self.dim})"
        )
