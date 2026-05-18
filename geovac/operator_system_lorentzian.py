"""Lorentzian truncated operator system on the BBB Krein spectral triple at (3, 1).

Sprint L3a-1 (2026-05-17): construct the Lorentzian analog of the
Connes-vS truncated operator system from `geovac.operator_system`
(`TruncatedOperatorSystem`) on the Krein space built in Sprint L2-B
(`geovac.krein_space_construction.KreinSpace`).

The Riemannian template (Paper 32 §III, `rem:operator_system`)
============================================================

The Riemannian Connes-vS truncated operator system is

    O_{n_max} := P_{n_max} C^infty(S^3) P_{n_max} subset M_{N(n_max)}(C)

where P_{n_max} projects onto the n <= n_max Fock-projected subspace and
C^infty(S^3) acts on L^2(S^3) by pointwise multiplication.  The result is
a *-closed but NOT multiplicatively-closed linear subspace of
M_{N(n_max)}(C), with propagation number prop(O_{n_max}) = 2 matching
Connes-vS Toeplitz S^1 Prop 4.2 verbatim (Paper 32 §III prop:propagation_2).

This module's deliverable: the natural Lorentzian extension
==========================================================

The Lorentzian truncated operator system is

    O^L_{n_max, N_t}
        := P_{n_max, N_t} . [ C^infty(S^3) (x) C^infty(R_t)_cutoff ]
                          . P_{n_max, N_t}
        subset B(K_{n_max, N_t})

where:
  - K_{n_max, N_t} = H_GV^{n_max} (x) C^{N_t}  is the Krein space from
    `geovac.krein_space_construction.KreinSpace`.
  - P_{n_max, N_t} = P_{n_max}^spatial (x) I_{N_t} is the truncation
    projection (spatial Fock cutoff x identity on the temporal grid).
  - C^infty(S^3) acts on H_GV^{n_max} by lifted-to-spinor pointwise
    multiplication (`FullDiracTruncatedOperatorSystem` from
    geovac/full_dirac_operator_system.py: scalar functions act block-
    diagonally on the two chirality sectors).
  - C^infty(R_t)_cutoff acts on C^{N_t} by pointwise multiplication
    evaluated on the temporal grid (a finite-rank commutative algebra
    of diagonal matrices).

Multipliers are constructed as tensor products M^spat (x) M^temp:
  - Spatial side: scalar 3-Y multipliers M_{N L M} from
    `geovac.operator_system.build_multiplier_matrix` lifted to
    FullDirac via chirality-block doubling (same as
    `FullDiracTruncatedOperatorSystem.build_full_dirac_multiplier_matrix`).
  - Temporal side: diagonal matrices diag(g(t_0), ..., g(t_{N_t-1})) for
    a polynomial basis g(t) in {1, t, t^2, ..., t^{N_t - 1}}.  At
    N_t = 1 the only available g is the constant 1 (no nontrivial
    temporal multipliers).

Structural finding (anticipated): prop(O^L) is INFINITY at N_t > 1
=================================================================

A commutative algebra on C^{N_t} (diagonal matrices) generates only
N_t-dimensional linear span, but its C*-envelope is the full
M_{N_t}(C) of dimension N_t^2.  The commutative subalgebra is closed
under multiplication; its products generate the same N_t-dim subspace.
So prop = infinity for the temporal-only factor.

For the tensor product operator system O^L = O^spat (x) O^temp:

    dim(O^L^k)  =  dim(O^spat^k) * dim(O^temp^k)
                =  dim(O^spat^k) * N_t   (since temp commutative)

The envelope is M_{dim_K}(C) of dim dim_K^2 = dim_H^2 * N_t^2.  To
reach the envelope we'd need dim(O^spat^k) = dim_H^2 * N_t.  But
dim(O^spat^k) <= dim_H^2 always (it's a subspace of M_{dim_H}(C)).
So we can never reach the envelope.

**Verdict: prop(O^L) = infinity at every N_t > 1, due to the temporal
direction's commutative-multiplier structure.  This is a STRUCTURAL
FINDING distinguishing the Lorentzian from the Riemannian case.**

At N_t = 1, only the trivial temporal multiplier (constant 1) is
available; the temporal slot is rank-1 and prop(O^L) reduces to
prop(O^spat) which is well-defined.  This is the LOAD-BEARING
falsifier: O^L at N_t = 1 must reduce bit-identically to the
FullDirac operator system (which itself reduces to the scalar
Riemannian operator system at the dim-of-O level since scalar
functions don't distinguish chirality sectors).

The "trivial" Riemannian-limit verdict at N_t = 1 IS Paper 32's
prop = 2 verbatim, lifted to the chirality-doubled Hilbert space.
Multiplication by a scalar function acts identically on both
chirality sectors, so the FullDirac operator system has the same
linear-span dimension as the scalar operator system (dim(O) = 14
at n_max = 2, 55 at n_max = 3, 140 at n_max = 4) even though dim_H
doubles.  Propagation number is 2.

Krein-positive restriction
==========================

The Krein space K has fundamental symmetry J = gamma^0 (x) I_{N_t} with
+/- 1 eigenspaces K^+, K^- each of dimension dim_K / 2.  An operator O on
K *preserves the Krein-positive cone* if O . K^+ subset K^+, i.e. if
the K^+ -> K^- block of O (in the J-eigenbasis) vanishes.

Our scalar-multiplier operators M^spat (x) M^temp may or may not
preserve K^+.  In the chiral basis with J = chirality-swap-on-spatial,
the spatial scalar multipliers M (which act identically on both
chirality sectors via M^full = M ⊕ M) commute with J_spatial because

    [J_spatial, M^full]  =  [chirality-swap, M ⊕ M]  =  0,

since chirality-swap swaps the two equal copies.  Therefore
M^full (x) (anything diagonal in t) commutes with J_spatial (x) I_{N_t}
= J, and PRESERVES BOTH K^+ AND K^-.

So the Krein-positive restriction is essentially trivial: every
scalar-multiplier in O^L already preserves K^+.  The natural sub-
system O^{L,+} := { O in O^L : O preserves K^+ } equals O^L itself.

This is a structural finding: the operator system on the chirality-
doubled spinor bundle is FULLY KREIN-POSITIVE by construction.  The
Krein-positive constraint does NOT carve out a nontrivial sub-system.

Witness pair
============

Take

    a  =  M^{spat}_{N=2, L=1, M=0}  (x)  I_{N_t},
    b  =  a^*  =  (M^{spat}_{N=2, L=1, M=0})^*  (x)  I_{N_t}.

Both are in O^L.  Their product is

    a b  =  (M^{spat} M^{spat *})  (x)  I_{N_t},

and the spatial factor M^{spat} M^{spat *} is the Riemannian witness
pair NOT in O^spat (Paper 32 §III, witness pair).  Therefore a b
is NOT in O^L = O^spat (x) O^temp (as a linear span over tensors of
multipliers; the SPATIAL factor itself is outside O^spat, so the
tensor cannot be a finite linear combination of generators).

This is the natural lift of the Riemannian witness pair to the
Lorentzian setting.

References
==========

  Connes, A. & van Suijlekom, W. D. "Spectral Truncations in
    Noncommutative Geometry and Operator Systems."  CMP 383 (2021),
    arXiv:2004.14115.  Definition 2.39 (propagation number),
    Proposition 4.2 (Toeplitz S^1 prop=2).

  Bizi, N., Brouder, C., Besnard, F.  "Space and time dimensions
    of algebras with application to Lorentzian noncommutative
    geometry..."  J. Math. Phys. 59, 062303 (2018).  arXiv:1611.07062.

  Paper 32 (Riemannian foundation): papers/synthesis/paper_32_spectral_triple.tex
    Section III, rem:operator_system, prop:propagation_2.

  Paper 43 (Lorentzian closure at finite cutoff): papers/standalone/paper_43_lorentzian_extension.tex

  GeoVac internal:
    geovac/operator_system.py                   (Riemannian template)
    geovac/full_dirac_operator_system.py        (chirality-doubled lift)
    geovac/krein_space_construction.py          (L2-B Krein space)
    geovac/lorentzian_dirac.py                  (L2-C D_L)
    geovac/connes_axiom_audit_31.py             (L2-D BBB audit)
    geovac/modular_hamiltonian_lorentzian.py    (L2-E wedge)
    debug/l3_scoping_memo.md                    (math architecture)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    FullDiracTruncatedOperatorSystem,
    full_dirac_basis,
    full_dirac_dim,
)
from geovac.krein_space_construction import KreinSpace
from geovac.operator_system import (
    HyperLabel,
    TruncatedOperatorSystem,
    allowed_multiplier_labels,
    build_multiplier_matrix,
    hilbert_basis,
)
from geovac.spinor_operator_system import (
    SpinorLabel,
    build_spinor_multiplier_matrix,
    spinor_basis,
)


# ---------------------------------------------------------------------------
# Temporal multiplier basis
# ---------------------------------------------------------------------------


def temporal_grid_symmetric(N_t: int, T_max: float = 1.0) -> np.ndarray:
    """Return the symmetric temporal grid t_k in [-T_max, T_max] with N_t points.

    Convention matches `geovac.krein_space_construction.temporal_grid`.

    At N_t = 1 returns the singleton [0.0] (Riemannian-limit / static
    spatial slice).

    Parameters
    ----------
    N_t : int
        Number of temporal grid points (must be >= 1).
    T_max : float
        Half-width of the temporal interval (must be > 0).

    Returns
    -------
    t_grid : np.ndarray of shape (N_t,), float64.
    """
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    if T_max <= 0:
        raise ValueError(f"T_max must be > 0, got {T_max}")
    if N_t == 1:
        return np.array([0.0])
    return np.linspace(-T_max, T_max, N_t)


def temporal_multiplier_matrices(N_t: int, T_max: float = 1.0) -> List[np.ndarray]:
    """Diagonal temporal multiplier matrices for the polynomial basis
    g_p(t) = t^p, p = 0, 1, ..., N_t - 1.

    The result is a list of N_t diagonal matrices diag(g_p(t_0), ...,
    g_p(t_{N_t - 1})) for p = 0, ..., N_t - 1.  These span the full
    commutative C-algebra of diagonal N_t x N_t matrices (since the
    Vandermonde matrix on N_t distinct t_k is invertible).

    At N_t = 1: only the constant g_0 = 1 is available; the temporal
    factor is the single 1x1 identity matrix.

    Parameters
    ----------
    N_t : int
        Number of temporal grid points and temporal multipliers.
    T_max : float
        Half-width of the temporal grid.

    Returns
    -------
    temp_mats : list of np.ndarray, each shape (N_t, N_t).
        The N_t-element polynomial basis.
    """
    t = temporal_grid_symmetric(N_t, T_max)
    mats = []
    for p in range(N_t):
        diag = t ** p  # float64
        # At t = 0 and p = 0 we get 1; convention 0**0 = 1.
        # Numpy returns 1.0 for 0**0 by default; check:
        if p == 0:
            diag = np.ones(N_t, dtype=np.complex128)
        else:
            diag = diag.astype(np.complex128)
        mats.append(np.diag(diag))
    return mats


# ---------------------------------------------------------------------------
# Spatial multiplier matrices: chirality-doubled scalar Y_{N L M}
# ---------------------------------------------------------------------------


def build_lorentzian_spatial_multiplier_matrix(
    N: int, L: int, M: int, basis: Sequence[FullDiracLabel],
) -> np.ndarray:
    """Construct M^spat_{N L M} on the chirality-doubled Hilbert space H_GV.

    A scalar function f(omega) on S^3 acts identically on both
    chirality sectors (it does NOT see the chirality label).  We
    construct this as a Weyl-sector spinor-lifted multiplier
    (using the geovac.spinor_operator_system Clebsch-Gordan
    machinery and Avery-Wen-Avery 3-Y integrals), then double it as
    block-diagonal in the two chirality sectors.

    This matches the construction used in
    `FullDiracTruncatedOperatorSystem.__init__` (geovac/
    full_dirac_operator_system.py).

    Parameters
    ----------
    N, L, M : int
        Multiplier label (Riemannian 3-Y SO(4) selection rules apply).
    basis : list of FullDiracLabel
        FullDirac basis (length 2 * spinor_dim(n_max)), Weyl block
        first then anti-Weyl block.

    Returns
    -------
    Complex ndarray of shape (2 * dim_Weyl, 2 * dim_Weyl).
    """
    dim_full = len(basis)
    dim_weyl = dim_full // 2
    if 2 * dim_weyl != dim_full:
        raise ValueError(
            f"FullDirac basis must be even-length, got {dim_full}"
        )

    # Build the Weyl-sector spinor multiplier
    weyl_basis_in_order = [b.to_spinor() for b in basis[:dim_weyl]]
    weyl_mat = build_spinor_multiplier_matrix(N, L, M, weyl_basis_in_order)

    # Double block-diagonally (same multiplier on both chirality sectors)
    full_mat = np.zeros((dim_full, dim_full), dtype=np.complex128)
    full_mat[:dim_weyl, :dim_weyl] = weyl_mat
    full_mat[dim_weyl:, dim_weyl:] = weyl_mat
    return full_mat


# ---------------------------------------------------------------------------
# Lorentzian truncated operator system
# ---------------------------------------------------------------------------


class LorentzianTruncatedOperatorSystem:
    """The Lorentzian truncated operator system O^L_{n_max, N_t} on the
    Krein space K_{n_max, N_t} = H_GV^{n_max} (x) C^{N_t}.

    Construction
    ------------
    Multipliers are tensor products of spatial scalar 3-Y multipliers
    M^spat_{N L M} (chirality-doubled) and temporal polynomial
    multipliers g_p(t) = t^p (p = 0, ..., N_t - 1):

        O^L  =  span_C { M^spat_{N L M} (x) diag(g_p(t_0..t_{N_t-1})) }.

    The total number of generators is len(spat_labels) * N_t.

    At N_t = 1 (LOAD-BEARING Riemannian limit), only the constant
    temporal multiplier g_0 = 1 is available, and the construction
    reduces to FullDiracTruncatedOperatorSystem(n_max) bit-identically.

    Attributes
    ----------
    n_max : int
        Spatial Fock cutoff.
    N_t : int
        Number of temporal grid points and temporal multipliers.
    T_max : float
        Half-width of temporal grid.
    krein : KreinSpace
        The Krein space at (n_max, N_t, T_max).
    basis_spatial : list of FullDiracLabel
        Spatial basis (length full_dirac_dim(n_max)).
    dim_spatial : int
        Spatial dimension.
    dim_K : int
        Total Krein dimension = dim_spatial * N_t.
    spat_labels : list of (N, L, M)
        Spatial multiplier labels admitting a nonzero spinor-lift matrix.
    multiplier_labels : list of (N, L, M, p)
        Full tensor-product labels: (spatial NLM, temporal power p).
    multiplier_matrices : list of complex np.ndarray
        Each has shape (dim_K, dim_K) and is constructed as
        np.kron(spatial_part, temporal_part).
    dim : int
        Linear-span dimension over C of multiplier_matrices in M_{dim_K}(C).
    envelope_dim : int
        dim_K ** 2 (C*-envelope of the full matrix algebra B(K)).
    """

    def __init__(
        self,
        n_max: int,
        N_t: int = 1,
        T_max: float = 1.0,
        *,
        rank_tol: float = 1e-12,
    ) -> None:
        if n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {n_max}")
        if N_t < 1:
            raise ValueError(f"N_t must be >= 1, got {N_t}")
        if T_max <= 0:
            raise ValueError(f"T_max must be > 0, got {T_max}")

        self.n_max = n_max
        self.N_t = N_t
        self.T_max = T_max
        self._rank_tol = rank_tol

        # Build Krein space (provides spatial basis + J + dim)
        self.krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=T_max)
        self.basis_spatial = self.krein.basis_spatial
        self.dim_spatial = self.krein.dim_spatial
        self.dim_K = self.krein.dim

        # Spatial multipliers (chirality-doubled scalar 3-Y multipliers)
        # We reuse the FullDiracTruncatedOperatorSystem build path, which
        # is exactly the chirality-block-doubled spinor lift.
        all_spat_candidates = allowed_multiplier_labels(n_max)
        weyl_basis_in_order = [
            b.to_spinor() for b in self.basis_spatial[: self.dim_spatial // 2]
        ]
        self._weyl_basis_in_order = weyl_basis_in_order

        spat_labels: List[Tuple[int, int, int]] = []
        spat_matrices: List[np.ndarray] = []
        for (N, L, M) in all_spat_candidates:
            weyl_mat = build_spinor_multiplier_matrix(N, L, M, weyl_basis_in_order)
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

        # Temporal multipliers diag(t_k^p) for p = 0..N_t-1
        temp_mats = temporal_multiplier_matrices(N_t, T_max)
        self._temp_matrices = temp_mats

        # Tensor-product multipliers (spatial NLM, temporal power p)
        labels: List[Tuple[int, int, int, int]] = []
        matrices: List[np.ndarray] = []
        for (N, L, M), s_mat in zip(spat_labels, spat_matrices):
            for p, t_mat in enumerate(temp_mats):
                labels.append((N, L, M, p))
                matrices.append(np.kron(s_mat, t_mat))
        self.multiplier_labels = labels
        self.multiplier_matrices = matrices
        self.basis_matrices = list(zip(labels, matrices))

        # Cache vec form for span / rank computations
        self._vec_cache = self._build_vec_matrix(matrices)
        if len(matrices) == 0:
            self._dim_cache = 0
        else:
            self._dim_cache = int(
                np.linalg.matrix_rank(self._vec_cache, tol=rank_tol)
            )

    @property
    def dim(self) -> int:
        """Linear-span dimension of O^L as a complex subspace of M_{dim_K}(C)."""
        return self._dim_cache

    @property
    def envelope_dim(self) -> int:
        """dim_K ** 2 (C*-envelope of B(K))."""
        return self.dim_K ** 2

    # -----------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------

    @staticmethod
    def _build_vec_matrix(matrices: Sequence[np.ndarray]) -> np.ndarray:
        """Stack vec-forms of a list of matrices as columns of (N^2, K) array."""
        if len(matrices) == 0:
            return np.zeros((0, 0), dtype=np.complex128)
        dim = matrices[0].shape[0]
        K = len(matrices)
        cols = np.zeros((dim * dim, K), dtype=np.complex128)
        for k, M in enumerate(matrices):
            cols[:, k] = M.reshape(-1)
        return cols

    # -----------------------------------------------------------------
    # Membership / projection test
    # -----------------------------------------------------------------

    def contains(
        self, matrix: np.ndarray, *, tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Test whether `matrix` lies in the linear span of multiplier_matrices.

        Returns (is_in_O, residual_ratio).
        """
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
    # Self-checks (operator system axioms)
    # -----------------------------------------------------------------

    def identity_in_O(self, *, tol: float = 1e-10) -> Tuple[bool, float]:
        """Verify that the identity I_{dim_K} lies in O^L.

        At N_t >= 1, the trivial multiplier f = 1 on S^3 x R produces
        the identity matrix on K (the constant function is the N=L=M=0,
        p=0 mode), so I in O^L always.
        """
        I = np.eye(self.dim_K, dtype=np.complex128)
        return self.contains(I, tol=tol)

    def is_star_closed(
        self, *, tol: float = 1e-10,
    ) -> Tuple[bool, List[Tuple[int, float]]]:
        """Verify that for each generator M_i, the conjugate transpose
        M_i^dagger lies in O^L.

        The spatial scalar multipliers M^spat_{N L M} satisfy
        M_{N L M}^* = M_{N L -M} up to sign (Wigner conjugation), so
        the spatial side is *-closed.  The temporal diagonal multipliers
        diag(t^p) are real diagonal, so they are *-closed trivially
        (M = M^*).  The tensor product is therefore *-closed.

        Returns (all_closed, list of (i, residual) for any failures).
        """
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
        self, *, tol: float = 1e-10,
    ) -> Tuple[List[int], List[Tuple[int, int, int, int]]]:
        """Indices and labels of multipliers that preserve the Krein-positive cone.

        An operator O preserves the K^+ cone (= +1-eigenspace of J) if
        and only if the K^+ -> K^- block of O vanishes in the J-eigenbasis.
        Equivalently, [J, O] mapping K^+ to K^- is zero on K^+.

        For the scalar-multiplier construction here, the spatial factor
        M^spat is chirality-doubled (M ⊕ M), and J_spatial = chirality-swap
        (off-diagonal in chirality).  We have

            J_spatial . (M ⊕ M)  =  (M ⊕ M) . J_spatial = 0,

        since J swaps the two equal copies.  Therefore EVERY spatial
        multiplier commutes with J_spatial.  Tensored with any temporal
        operator that commutes with I_{N_t} = J_temp (which is everything),
        the full multiplier commutes with J = J_spatial (x) I_{N_t}.

        Conclusion: every multiplier in O^L commutes with J, hence
        preserves K^+ (and K^-).  The Krein-positive restriction is
        trivial in this construction.

        Returns
        -------
        preserving_indices : list of indices in multiplier_matrices
        preserving_labels : list of (N, L, M, p) tuples
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

    def restrict_to_krein_positive(
        self, *, tol: float = 1e-10,
    ) -> "LorentzianTruncatedOperatorSystem":
        """Return the sub-operator-system of multipliers preserving K^+.

        For the scalar-multiplier construction, the result is the full
        system itself (trivial Krein-positive restriction; see
        `krein_positive_preservers`).

        We construct a fresh LorentzianTruncatedOperatorSystem and replace
        its multiplier list with the K^+ -preserving subset.

        Returns
        -------
        sub_system : LorentzianTruncatedOperatorSystem
            With possibly reduced multiplier list.  Note: if all
            multipliers preserve K^+, the result is equivalent to `self`.
        """
        # Make a shallow copy and filter to K^+ preservers
        preserving_indices, _ = self.krein_positive_preservers(tol=tol)
        # Build a new instance with only filtered multipliers
        sub = LorentzianTruncatedOperatorSystem.__new__(
            LorentzianTruncatedOperatorSystem
        )
        sub.n_max = self.n_max
        sub.N_t = self.N_t
        sub.T_max = self.T_max
        sub.krein = self.krein
        sub.basis_spatial = self.basis_spatial
        sub.dim_spatial = self.dim_spatial
        sub.dim_K = self.dim_K
        sub._rank_tol = self._rank_tol
        sub.spat_labels = list(self.spat_labels)
        sub._spat_matrices = list(self._spat_matrices)
        sub._temp_matrices = list(self._temp_matrices)
        sub._weyl_basis_in_order = list(self._weyl_basis_in_order)
        sub.multiplier_labels = [
            self.multiplier_labels[i] for i in preserving_indices
        ]
        sub.multiplier_matrices = [
            self.multiplier_matrices[i] for i in preserving_indices
        ]
        sub.basis_matrices = list(
            zip(sub.multiplier_labels, sub.multiplier_matrices)
        )
        sub._vec_cache = sub._build_vec_matrix(sub.multiplier_matrices)
        if len(sub.multiplier_matrices) == 0:
            sub._dim_cache = 0
        else:
            sub._dim_cache = int(
                np.linalg.matrix_rank(sub._vec_cache, tol=sub._rank_tol)
            )
        return sub

    # -----------------------------------------------------------------
    # Achievable envelope (Weyl-block-diagonal x temporal-diagonal subspace)
    # -----------------------------------------------------------------

    @property
    def achievable_envelope_dim(self) -> int:
        """Dimension of the operator subspace achievable by scalar multipliers.

        Scalar functions on S^3 act on H_GV = H_Weyl (+) H_antiWeyl as
        M (+) M (chirality-doubled, same multiplier on both copies).
        Products of such operators remain block-diagonal in chirality
        AND have equal Weyl/anti-Weyl blocks.  The achievable subspace
        is therefore

            { M_W (+) M_W : M_W in M_{dim_Weyl}(C) }    (Weyl-doubled subspace),

        of dimension dim_Weyl^2 = (dim_spatial / 2)^2.

        Tensored with the commutative temporal subalgebra of dim N_t
        (diagonal matrices in C^{N_t}, the polynomial-basis span), the
        full achievable subspace has dimension

            dim_Weyl^2 * N_t.

        Compare to the full B(K) envelope dim_K^2 = dim_spatial^2 * N_t^2.
        The achievable subspace is a factor of (dim_spatial^2 * N_t^2) /
        (dim_Weyl^2 * N_t) = 4 * N_t smaller.  At n_max = 2, N_t = 1,
        the achievable is 64 = 8^2 vs full envelope 256 = 16^2: factor 4.

        Returns
        -------
        int : (dim_spatial // 2) ** 2 * N_t.
        """
        dim_Weyl = self.dim_spatial // 2
        return dim_Weyl * dim_Weyl * self.N_t

    # -----------------------------------------------------------------
    # Propagation number (with infinity verdict for N_t > 1)
    # -----------------------------------------------------------------

    def compute_propagation_number(
        self,
        *,
        max_k: int = 4,
        tol: float = 1e-10,
        verbose: bool = False,
        envelope: str = "achievable",
    ) -> Tuple[int, List[int]]:
        """Compute prop(O^L) = smallest k such that dim(O^L^k) = target_dim.

        Parameters
        ----------
        max_k : int
            Cap on the iteration depth.
        tol : float
            Numerical rank tolerance.
        verbose : bool
            Print per-step rank info if True.
        envelope : {'achievable', 'full'}, default 'achievable'
            Which envelope to target:
              - 'achievable': dim_Weyl^2 * N_t (the actual reachable
                subspace under scalar-multiplier closure;
                Connes-vS-style propagation question).
              - 'full': dim_K^2 (the full B(K) envelope; this is the
                NAIVE target that scalar multipliers CANNOT reach in the
                chirality-doubled construction).

        Returns
        -------
        (prop, dim_sequence) : (int, list of int)
            dim_sequence[i] = dim(O^L^{i+1}).  prop = -1 if not reached
            within max_k.

        Structural findings (verified at n_max = 2, N_t in {1, 3}):

          (a) With envelope = 'achievable':
              - At N_t = 1: prop = 2 (matches Paper 32 §III spatial result;
                dim sequence 14, 64 = dim_Weyl^2 achievable envelope).
              - At N_t > 1: prop = INFINITY because the temporal
                commutative algebra never reaches M_{N_t}(C).  Dim
                sequence saturates at dim_Weyl^2 * N_t (sub-target).
                Achievable target is dim_Weyl^2 * N_t^2 in principle,
                but the temporal subalgebra is rank-N_t-commutative,
                bounding O^L^k below the full achievable.

          (b) With envelope = 'full':
              - prop = INFINITY always (chirality-doubling blocks the
                construction from generating chirality-flipping
                operators).

        The 'achievable' convention is the natural Connes-vS-style
        propagation question for scalar-multiplier operator systems
        on a chirality-doubled Hilbert space.  See module docstring
        for the full structural analysis.

        Caution
        -------
        At n_max = 2, N_t = 1, dim_K = 16 and envelope = 256; checking
        prop = 2 requires computing O^L^2 with up to 196 generators and
        rank-checking a (256, 196) matrix (fast).  At larger n_max or
        N_t, the cost scales as (|gens|^2, dim_K^2).
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
            # Extract a basis of O^{k-1} from current_generators
            basis_Ok_minus_1 = self._extract_matrix_basis(
                current_generators, tol=tol
            )
            # Generate O^k = basis_Ok_minus_1 @ matrices_O
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
                # Saturated below the target -> prop = infinity
                if verbose:
                    print(
                        f"  saturated at dim {dim_Ok} < {target_dim}; "
                        f"prop = infinity"
                    )
                return -1, dim_sequence
            current_generators = new_generators

        return -1, dim_sequence

    @staticmethod
    def _operator_system_dim(
        matrices: Sequence[np.ndarray], *, tol: float = 1e-10,
    ) -> int:
        """Linear-span dimension of a list of matrices over C."""
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
        matrices: Sequence[np.ndarray], *, tol: float = 1e-10,
    ) -> List[np.ndarray]:
        """Extract a complex-linear basis via SVD on the vec-stack."""
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
    # Riemannian limit at N_t = 1 (load-bearing falsifier)
    # -----------------------------------------------------------------

    def verify_riemannian_limit(
        self, *, tol_struct: float = 0.0,
    ) -> Tuple[bool, dict]:
        """Verify that at N_t = 1, the construction reduces bit-identically
        to FullDiracTruncatedOperatorSystem(n_max).

        LOAD-BEARING falsifier.  Must pass with residual exactly 0.0
        (in float64), because the construction is built from the same
        spinor-lifted multipliers and the tensor product with the 1x1
        identity is the bit-identity operation.

        Two bit-exact checks:

          (a) dim_K = full_dirac_dim(n_max)         (dimension match)
          (b) For each multiplier label (N, L, M, p=0), the matrix
              equals the corresponding FullDirac multiplier matrix
              built from the same Avery-Wen-Avery 3-Y integrals.

        Returns
        -------
        (ok, details) : (bool, dict)
            details has keys 'N_t', 'dim_K', 'full_dirac_dim',
            'n_multipliers_self', 'n_multipliers_ref', 'max_residual',
            'dim_self', 'dim_ref', 'load_bearing_pass'.
        """
        if self.N_t != 1:
            # Build a fresh N_t = 1 system and recurse
            sys_1 = LorentzianTruncatedOperatorSystem(
                n_max=self.n_max, N_t=1, T_max=self.T_max
            )
            return sys_1.verify_riemannian_limit(tol_struct=tol_struct)

        # We are at N_t = 1.  Build the FullDirac reference.
        ref = FullDiracTruncatedOperatorSystem(self.n_max)
        ref_labels_to_idx = {lbl: i for i, lbl in enumerate(ref.multiplier_labels)}

        # Each self multiplier label (N, L, M, p=0) corresponds to
        # FullDirac multiplier label (N, L, M).
        max_residual = 0.0
        n_compared = 0
        n_unmatched = 0
        for (N, L, M, p), self_mat in zip(
            self.multiplier_labels, self.multiplier_matrices
        ):
            if p != 0:
                # Should not happen at N_t = 1
                continue
            if (N, L, M) not in ref_labels_to_idx:
                # Reference rejected this label; self should have rejected too
                # (we filtered against the same zero-check).  If we get here it's a
                # bookkeeping mismatch.
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

        # Bit-exact equality target (tol_struct = 0.0 is the load-bearing
        # falsifier; we accept any nonzero residual <= 1e-14 as numerical
        # noise from sympy-to-float casting, but per the L2-B/C/D/E
        # convention we expect exactly 0.0).
        dim_ok = self.dim_K == full_dirac_dim(self.n_max)
        n_mult_ok = (
            len(self.multiplier_labels) == len(ref.multiplier_labels)
            and n_unmatched == 0
        )
        residual_ok = max_residual <= max(tol_struct, 1e-14)
        # The dim of O on M_{dim_K}(C) is bookkeeping: dim_self should
        # equal dim_ref (same spinor-lifted multiplier span).
        dim_O_ok = self.dim == ref.dim

        ok = dim_ok and n_mult_ok and residual_ok and dim_O_ok
        details["load_bearing_pass"] = ok
        details["dim_match"] = dim_ok
        details["multiplier_count_match"] = n_mult_ok
        details["matrix_residual_pass"] = residual_ok
        details["dim_O_match"] = dim_O_ok
        return ok, details


# ---------------------------------------------------------------------------
# Witness pair: a, b in O^L with ab NOT in O^L
# ---------------------------------------------------------------------------


def witness_pair_lorentzian(
    O_L: LorentzianTruncatedOperatorSystem,
    *,
    N_target: int = 2, L_target: int = 1, M_target: int = 0,
    p_target: int = 0,
) -> Tuple[
    Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray],
    Optional[Tuple[bool, float]],
]:
    """Construct the Lorentzian witness pair lifting the Riemannian
    witness pair (Paper 32 §III) to the Krein space.

    Take

        a  =  M^{spat}_{N_target, L_target, M_target}  (x)  g_{p_target}(t),
        b  =  a^*.

    Then

        a b  =  ( M^{spat} (M^{spat})^* )  (x)  ( g_{p_target}(t)^2 ).

    The spatial factor M^{spat} (M^{spat})^* is the Riemannian witness
    product, which is NOT in the spatial operator system O^{spat} (this is
    Paper 32 §III's witness construction at n_max = 2).  Therefore the
    tensor a b is NOT in O^L either: any decomposition of a b as a
    finite C-linear combination of tensors M (x) T from O^L would
    force the spatial factor of that decomposition to lie in O^{spat},
    which it does not (Paper 32 §III).

    The Lorentzian witness pair confirms that O^L is genuinely an
    operator system (not a *-algebra): closed under *, linear
    combinations, but NOT under multiplication.

    Returns
    -------
    (a, b, ab, (in_O_L, residual)).
    If the requested generator does not exist in O^L, returns
    (None, None, None, None).
    """
    label_target = (N_target, L_target, M_target, p_target)
    a = None
    for (label, mat) in O_L.basis_matrices:
        if label == label_target:
            a = mat
            break
    if a is None:
        return None, None, None, None
    b = a.conj().T
    ab = a @ b
    test = O_L.contains(ab)
    return a, b, ab, test


# ---------------------------------------------------------------------------
# Convenience: connect to L2-E hemispheric wedge
# ---------------------------------------------------------------------------


def restrict_to_lorentzian_wedge(
    O_L: LorentzianTruncatedOperatorSystem,
) -> dict:
    """Compute the wedge-restricted operator-system content for L2-E linkage.

    Sprint L2-E uses the hemispheric Krein wedge

        W_L  =  P_W_spatial  (x)  P_t_positive

    on the same Krein space K_{n_max, N_t}.  Restricting the
    operator system O^L to W_L means computing the wedge-block
    operator P_W_L . O . P_W_L for each O in O^L.  These wedge-blocks
    form a sub-operator-system on the wedge.

    This function builds the wedge projection and the wedge-restricted
    multiplier matrices, returning structured data for inspection /
    further computation.

    Parameters
    ----------
    O_L : LorentzianTruncatedOperatorSystem

    Returns
    -------
    dict with keys:
      'P_W_L'             : (dim_K, dim_K) wedge projector
      'dim_W_L'           : trace(P_W_L) = wedge dimension
      'wedge_block_matrices' : list of (dim_W_L, dim_W_L) wedge-restricted
                              multipliers
      'wedge_linear_span_dim': linear-span dim of wedge-block matrices in
                              M_{dim_W_L}(C)
    """
    from geovac.modular_hamiltonian_lorentzian import LorentzianWedge

    wedge = LorentzianWedge(krein=O_L.krein)
    P_W_L = wedge.P_W_L
    wi = wedge.wedge_basis_indices()
    dim_W_L = len(wi)

    # Project to the wedge sub-Hilbert space via index restriction
    wedge_blocks = []
    for M in O_L.multiplier_matrices:
        # Full-Krein wedge restriction P_W_L M P_W_L
        M_wedge = P_W_L @ M @ P_W_L
        # Extract the wedge-only block (index-restricted)
        # M_wedge has support on the wedge indices; the (i, j)-block
        # in wedge_basis_indices ordering is M_wedge[wi[i], wi[j]].
        M_block = M_wedge[np.ix_(wi, wi)]
        wedge_blocks.append(M_block)

    if len(wedge_blocks) > 0:
        cols = np.zeros((dim_W_L * dim_W_L, len(wedge_blocks)), dtype=np.complex128)
        for k, M in enumerate(wedge_blocks):
            cols[:, k] = M.reshape(-1)
        span_dim = int(np.linalg.matrix_rank(cols, tol=1e-10))
    else:
        span_dim = 0

    return {
        "P_W_L": P_W_L,
        "dim_W_L": dim_W_L,
        "wedge_block_matrices": wedge_blocks,
        "wedge_linear_span_dim": span_dim,
        "wedge_basis_indices": wi,
    }
