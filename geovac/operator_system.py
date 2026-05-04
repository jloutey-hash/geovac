"""Truncated operator system O_{n_max} = P_{n_max} C^infty(S^3) P_{n_max}.

This module realizes the Connes & van Suijlekom (CMP 2021, arXiv:2004.14115)
spectral truncation of the unit-S^3 spectral triple at finite Fock cutoff
n_max. It implements:

  - The basis of O_{n_max} as a *-closed but NOT multiplicatively-closed
    linear subspace of M_{N(n_max)}(C), where N(n_max) = sum_{n=1..n_max} n^2.
  - A `contains` method that tests whether a given matrix lies in O via
    least-squares projection in vec-form.
  - The propagation number prop(O_{n_max}) via incremental computation of
    dim(O^k) until it equals N^2 (the C*-envelope dimension).
  - A witness-pair check: there exist a, b in O with a*b NOT in O.

Mathematical setup
==================

Hilbert space
-------------

The truncated Hilbert space is

    H_{n_max} = span{ |n, l, m> : 1 <= n <= n_max, 0 <= l <= n-1, -l <= m <= l }

with dimension N(n_max) = sum_{n=1..n_max} n^2 = 1, 5, 14, 30, 55, ...

The basis vectors are the (real-orthonormalized) hyperspherical 3-spherical
harmonics Y^{(3)}_{n l m}(omega) on S^3, the Fock-projection images of the
(non-relativistic) hydrogen wavefunctions. Here n is the SO(4) principal
quantum number; equivalently, in the SU(2)_L x SU(2)_R = Spin(4) double cover
language, the n-th SO(4) irrep is the symmetric tensor representation
((n-1)/2, (n-1)/2) of dimension n^2.

Operator system
---------------

The Connes-vS truncated operator system is

    O_{n_max} := P_{n_max} C^infty(S^3) P_{n_max} subset M_N(C)

where P_{n_max} is the orthogonal projection onto H_{n_max} viewed inside
L^2(S^3), and C^infty(S^3) acts on L^2(S^3) by pointwise multiplication.

A function f in C^infty(S^3) expands in the hyperspherical-harmonic basis as

    f(omega) = sum_{N, L, M} c_{N L M} * Y^{(3)}_{N L M}(omega).

Each multiplier Y^{(3)}_{N L M} produces a matrix M_{N L M} in M_N(C) with
entries

    <n, l, m | Y^{(3)}_{N L M} | n', l', m'>
       = integral_{S^3} conj(Y^{(3)}_{n l m}) * Y^{(3)}_{N L M} * Y^{(3)}_{n' l' m'} dOmega.

The linear span (over C) of all such matrices, as (N, L, M) range over all
SO(4) irrep labels for which the matrix is nonzero somewhere in the cutoff,
is O_{n_max}.

Selection rules for the 3-Y integral on S^3
-------------------------------------------

The 3-spherical-harmonic integral factorizes through SO(4) Clebsch-Gordan
into an angular S^2 part (a standard SU(2) Gaunt integral on the 2-sphere)
times an SO(4) "principal-QN coupling" part. Concretely:

The hyperspherical harmonic Y^{(3)}_{n l m} factorizes as

    Y^{(3)}_{n l m}(chi, theta, phi) = R_{n l}(chi) * Y_{l m}(theta, phi)

where chi in [0, pi] is the polar S^3 angle, (theta, phi) are the standard
S^2 coordinates, R_{n l}(chi) is a Gegenbauer-polynomial-based radial-on-the-
sphere factor, and Y_{l m} is the ordinary 2-sphere spherical harmonic.

The 3-Y integral therefore splits as

    integral Y^{(3) *}_{n l m} Y^{(3)}_{N L M} Y^{(3)}_{n' l' m'} dOmega_3
       = G_{l L l'}^{m M m'}    (angular S^2 Gaunt)
         * I_{n N n'}^{l L l'}   (SO(4) radial-on-sphere overlap)

where:

  G_{l L l'}^{m M m'} = integral Y_{l m}^* Y_{L M} Y_{l' m'} dOmega_2

is the standard 2-sphere Gaunt integral, and I_{n N n'}^{l L l'} is the
overlap integral over chi of three Gegenbauer-based radial factors.

The 2-sphere Gaunt integral has the closed form (see Wikipedia
"Gaunt coefficient" or Edmonds 1957):

    G_{l L l'}^{m M m'} = (-1)^m * sqrt((2 l + 1)(2 L + 1)(2 l' + 1)/(4 pi))
                          * (l L l' ; 0 0 0) * (l L l' ; -m M m')

where (j1 j2 j3 ; m1 m2 m3) is the Wigner 3j symbol. Selection rules:

    M = m - m'  (conservation of m, here m_l - m_l' on the bra-ket convention)
    |l - l'| <= L <= l + l'
    l + l' + L is even   (parity / 3j(0,0,0) selection)

For the SO(4) radial overlap on S^3, the selection rule from SO(4)
representation theory is (Avery 1989, Eq. 2.5; Wen-Avery JMP 26, 396 (1985)):

    |n - n'| + 1 <= N <= n + n' - 1     (SO(4) triangle on principal QN)
    L <= N - 1                          (L is bounded by N as l is bounded by n)

When the triangle and parity selection rules are simultaneously satisfied, the
integral I_{n N n'}^{l L l'} is a rational number times a possibly-rational
factor depending on conventions. For the WH1-round-2 task we need only:

  (a) WHETHER the matrix M_{N L M} has any nonzero entries in the truncated
      basis (selection rules);
  (b) The off-diagonal STRUCTURE of M_{N L M} (which (n, l, m) -> (n', l', m')
      pairs are coupled), since the propagation number depends only on the
      LINEAR SPAN of the M matrices, not their precise nonzero values.

For the propagation number computation, the precise transcendental values of
the radial overlaps I_{n N n'}^{l L l'} cancel out: the dimension of the
linear span of {M_{N L M}} as a subspace of M_N(C) depends only on the
support pattern, not on the values, modulo accidental linear dependences.

Witness pair construction
=========================

We will take a, b to be the matrices

    a = M_{N=2, L=1, M=0}      (a "raising" multiplier connecting n -> n+1)
    b = a^*                     (the conjugate "lowering" multiplier)

For the witness check, a*b is "raise then lower." On the interior of the
truncation (n, l, m) -> (n+1, l +/- 1, m) -> (n, l, m) is allowed, but at the
top shell n = n_max the raising step sends the state outside the truncation,
where it is killed by P, so the diagonal entry of a*b at the top shell is
"missing" the contribution from the unreachable n_max + 1 shell. This produces
a deficit on the top-shell diagonal that no element of O can supply (since
every M_{N L M} couples specific shell pairs symmetrically).

Concretely we will exhibit a 5x5 example at n_max = 2.

Propagation number computation
==============================

We compute prop(O_{n_max}) operationally via:

    prop = smallest k >= 1 such that dim(O^k) = N^2

where O^k is the linear span (in M_N(C)) of all products of <= k elements of
O. The C*-envelope is M_N(C) here because Toeplitz-style 3-Y multipliers
generate full matrix algebras (mirroring Connes-vS Prop 4.2 for S^1).

Conjecture (round 2): prop(O_{n_max}) = 2 for all n_max >= 2, matching
Toeplitz S^1.

Numerical approach
==================

We exact-arithmetically compute the symbolic 3-Y matrices using
sympy.physics.wigner.wigner_3j. For the SO(4) radial overlap I_{n N n'}^{l L l'}
we use a CONVENTION-INDEPENDENT placeholder: any nonzero rational seed
suffices for span / rank computations, and for the dimension of O^k the
specific value cancels. We use sympy Rationals and cast to float for rank
computation only. Witness pair details:

  - For "ab in O?" testing we need the actual values of M, but only up to
    overall normalization of each row/column (so any nonzero rational
    placeholder again suffices for the qualitative span test).
  - We are not computing the Connes distance here; that requires actual
    radial-overlap values (round 3).

Self-checks built in:

  - 1 in O                         (the identity multiplier, L=N=M=0 mode)
  - O is *-closed                  (each M_{N L M}^* equals M_{N L -M} up to
                                    a sign factor)
  - dim(O^prop) = N^2              (C*-envelope is full matrix algebra)
  - dim(O^k) < N^2 for k < prop    (propagation is correct)

References
==========

A. Connes & W. D. van Suijlekom, "Spectral Truncations in Noncommutative
Geometry and Operator Systems," CMP 383 (2021), arXiv:2004.14115. See
Definition 2.39 (propagation number) and Proposition 4.2 (Toeplitz S^1
prop = 2 independent of n).

J. Avery, "Hyperspherical Harmonics: Applications in Quantum Theory," Kluwer
1989, Eq. 2.5 and 3.5 for the SO(4) selection rules and 3-Y integral
factorization.

Z. Wen & J. Avery, "Some properties of hyperspherical harmonics," J. Math.
Phys. 26, 396-403 (1985), for the closed-form 3-Y integral on S^d.

Round-1 mapping memo: debug/wh1_connes_vs_mapping_memo.md (with PDF-erratum
at top, 2026-05-03).

Round-1 PDF verification: debug/wh1_connes_vs_pdf_verification.md
(2026-05-03).
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
import sympy as sp
from sympy.physics.wigner import wigner_3j


# ---------------------------------------------------------------------------
# Basis indexing for H_{n_max}
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HyperLabel:
    """A basis label (n, l, m) for the hyperspherical harmonic Y^{(3)}_{nlm}.

    n: SO(4) principal quantum number, n >= 1.
    l: angular momentum on S^2, 0 <= l <= n - 1.
    m: magnetic quantum number, -l <= m <= l.
    """

    n: int
    l: int
    m: int

    def __post_init__(self) -> None:  # pragma: no cover - sanity only
        if not (self.n >= 1 and 0 <= self.l <= self.n - 1 and -self.l <= self.m <= self.l):
            raise ValueError(f"invalid HyperLabel {self}")


def hilbert_basis(n_max: int) -> List[HyperLabel]:
    """Return the basis of H_{n_max} in canonical order.

    Order: by n first (1, 2, ..., n_max); within each n by l (0, ..., n-1);
    within each (n, l) by m (-l, ..., +l).

    The dimension is N(n_max) = sum_{n=1..n_max} n^2.
    """
    basis = []
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                basis.append(HyperLabel(n=n, l=l, m=m))
    return basis


def hilbert_dim(n_max: int) -> int:
    """N(n_max) = sum_{n=1..n_max} n^2."""
    return sum(n * n for n in range(1, n_max + 1))


# ---------------------------------------------------------------------------
# 3-Y integral on S^3: angular Gaunt + SO(4) selection rules
# ---------------------------------------------------------------------------


def _s2_gaunt_symbolic(l: int, L: int, lp: int, m: int, M: int, mp: int) -> sp.Expr:
    """Standard 2-sphere Gaunt integral G_{l L l'}^{m M m'}.

    G = (-1)^m * sqrt((2l+1)(2L+1)(2l'+1)/(4 pi))
         * (l L l' ; 0 0 0) * (l L l' ; -m M m')

    Returns a sympy expression. Returns 0 if any selection rule fails.
    """
    if M != m - mp:
        return sp.Integer(0)
    if not (abs(l - lp) <= L <= l + lp):
        return sp.Integer(0)
    if (l + L + lp) % 2 != 0:
        return sp.Integer(0)
    if abs(M) > L or abs(m) > l or abs(mp) > lp:
        return sp.Integer(0)

    threej_zero = wigner_3j(l, L, lp, 0, 0, 0)
    threej_m = wigner_3j(l, L, lp, -m, M, mp)
    if threej_zero == 0 or threej_m == 0:
        return sp.Integer(0)

    sign = sp.Integer(-1) ** m
    pref = sp.sqrt(sp.Rational((2 * l + 1) * (2 * L + 1) * (2 * lp + 1)) / (4 * sp.pi))
    return sign * pref * threej_zero * threej_m


def _so4_radial_overlap_placeholder(
    n: int, N: int, np_: int, l: int, L: int, lp: int
) -> sp.Expr:
    """Convention-independent nonzero placeholder for I_{n N n'}^{l L l'}.

    LEGACY/round-2 placeholder. Kept for the regression test
    `test_propagation_number_robust_to_placeholder` and for fallback when
    a fast (non-Avery) value is wanted. The DEFAULT for matrix construction
    is now the genuine Avery-Wen-Avery Gegenbauer 3-integral (sprint
    WH1-R3.1, see _so4_radial_overlap_avery below).

    The SO(4) "naive" selection rule we used for the placeholder

        |n - n'| + 1 <= N <= n + n' - 1     (SO(4) triangle)
        L <= N - 1                          (subhood)

    is too restrictive on the (n, l) labels: the genuine Gegenbauer
    triple integral (Avery 1989 Eq. 2.5 + Wen-Avery JMP 26, 396 1985) is
    nonzero on more (n, l, N, L, n', l') configurations than this triangle
    permits, because the sin^l(chi) prefactors in the radial functions
    R_{n,l}(chi) interact with the Gegenbauer polynomials to produce
    additional nonzero couplings. The placeholder was intentionally
    over-restrictive: it captured the support pattern that suffices for
    the propagation-number computation (verified to give prop = 2 robust
    to placeholder choice, round 2).

    We choose:

        I_{n, N, n'}^{l, L, l'} =
            1                          if N == 1 (identity-multiplier branch)
            1 + 1/(N+1) * (n + n')     otherwise (generic distinctness)

    when the (over-restrictive) SO(4) selection rules fire, else 0.
    """
    if not (abs(n - np_) + 1 <= N <= n + np_ - 1):
        return sp.Integer(0)
    if not L <= N - 1:
        return sp.Integer(0)
    if N == 1:
        # Trivial irrep: constant multiplier, identity matrix.
        return sp.Integer(1)
    # Generic positive-rational placeholder, distinct on most index-tuples.
    # Symmetric in (n, n') to keep the operator system *-closed.
    return sp.Integer(1) + sp.Rational(n + np_, N + 1)


def _so4_radial_overlap_avery(
    n: int, N: int, np_: int, l: int, L: int, lp: int
) -> sp.Expr:
    """Genuine Avery-Wen-Avery Gegenbauer 3-integral (sprint WH1-R3.1).

    Drop-in replacement for `_so4_radial_overlap_placeholder` returning the
    closed-form physical value

        I_{n N n'}^{l L l'} = int_0^pi R_{n l}(chi) R_{N L}(chi)
                                       R_{n' l'}(chi) sin^2(chi) dchi

    with R_{n l}(chi) = N_{n l} sin^l(chi) C^{l+1}_{n-l-1}(cos chi) and
    Avery 1989 Eq. 1.6.3 normalization. See `geovac.so4_three_y_integral`
    for the closed-form derivation and orthonormality verification.

    Returns sympy expression in exact arithmetic (rational + sqrt(rational)
    times powers of pi).
    """
    # Imported lazily to keep operator_system standalone-loadable.
    from geovac.so4_three_y_integral import gegenbauer_triple_integral
    return gegenbauer_triple_integral(n, l, N, L, np_, lp)


# Active radial-overlap function used by `three_y_integral`. Switch to
# the placeholder by setting this to `_so4_radial_overlap_placeholder`
# (e.g. to reproduce round-2 placeholder-convention numerical values).
_so4_radial_overlap = _so4_radial_overlap_avery


def three_y_integral(
    n: int, l: int, m: int,
    N: int, L: int, M: int,
    np_: int, lp: int, mp: int,
) -> sp.Expr:
    """The 3-Y integral on S^3 using factorization into S^2 Gaunt and SO(4)
    radial overlap.

    integral conj(Y^{(3)}_{n l m}) * Y^{(3)}_{N L M} * Y^{(3)}_{n' l' m'}
        = I_{n N n'}^{l L l'} * G_{l L l'}^{m M m'}.

    Conjugation convention: for real-form spherical harmonics this is just
    the algebraic product. For complex Y_{l m} we have
    conj(Y_{l m}) = (-1)^m Y_{l, -m}. The factor (-1)^m is absorbed into the
    Gaunt coefficient by the standard formula (which has the (-1)^m
    prefactor).
    """
    # Compute angular factor first (typically faster + more often zero by
    # the strict m_a = m_b + m_c rule). Saves cost when angular is zero.
    angular = _s2_gaunt_symbolic(l, L, lp, m, M, mp)
    if angular == 0:
        return sp.Integer(0)
    radial = _so4_radial_overlap(n, N, np_, l, L, lp)
    if radial == 0:
        return sp.Integer(0)
    return radial * angular


# ---------------------------------------------------------------------------
# Allowed multiplier labels (N, L, M)
# ---------------------------------------------------------------------------


def allowed_multiplier_labels(n_max: int) -> List[Tuple[int, int, int]]:
    """List all (N, L, M) triples for which M_{N L M} could be nonzero in the
    truncation H_{n_max}.

    The bounds:
      N: any integer with N >= 1 and N <= 2 * n_max - 1 (since the SO(4)
         triangle |n - n'| + 1 <= N <= n + n' - 1 with n, n' <= n_max gives
         N <= 2 * n_max - 1).
      L: 0 <= L <= N - 1.
      M: -L <= M <= L.

    Note: this list is the maximal set; some labels may produce all-zero
    matrices in M_N(C) due to the joint angular-S^2 + SO(4) selection rules.
    Those are filtered out at matrix-construction time.
    """
    labels = []
    N_max = 2 * n_max - 1
    for N in range(1, N_max + 1):
        for L in range(N):  # 0 <= L <= N - 1
            for M in range(-L, L + 1):
                labels.append((N, L, M))
    return labels


# ---------------------------------------------------------------------------
# Multiplier matrices M_{N L M} as numpy arrays
# ---------------------------------------------------------------------------


def build_multiplier_matrix(
    N: int, L: int, M: int, basis: Sequence[HyperLabel],
) -> np.ndarray:
    """Construct the matrix M_{N L M} in M_dim(C) with entries
    <n, l, m | Y^{(3)}_{N L M} | n', l', m'>.

    Returns a complex numpy array of shape (dim, dim).

    NOTE: We use the convention-independent placeholder for the radial
    overlap (see _so4_radial_overlap_placeholder); the propagation-number
    computation depends only on the support pattern of M and not on the
    placeholder values.
    """
    dim = len(basis)
    mat = np.zeros((dim, dim), dtype=np.complex128)
    for i, bra in enumerate(basis):
        for j, ket in enumerate(basis):
            val = three_y_integral(
                bra.n, bra.l, bra.m,
                N, L, M,
                ket.n, ket.l, ket.m,
            )
            if val == 0:
                continue
            # Cast sympy expression (with sqrt(...) and 1/sqrt(pi)) to float.
            mat[i, j] = complex(sp.N(val, 30))
    return mat


# ---------------------------------------------------------------------------
# Truncated operator system
# ---------------------------------------------------------------------------


class TruncatedOperatorSystem:
    """The Connes-vS truncated operator system O_{n_max} on S^3.

    Attributes
    ----------
    n_max : int
        The Fock cutoff.
    basis : list of HyperLabel
        Ordered orthonormal basis of H_{n_max}, dim = N(n_max).
    dim_H : int
        N(n_max) = sum_{n=1..n_max} n^2.
    multiplier_labels : list of (N, L, M)
        Allowed multiplier labels (those producing a nonzero matrix).
    multiplier_matrices : list of complex numpy arrays
        Corresponding M_{N L M} matrices, each shape (dim_H, dim_H).
    basis_matrices : list of (label, matrix) pairs
        Synonym pairing for inspection.
    dim : int
        dim of O as a complex linear subspace of M_{dim_H}(C).
    """

    def __init__(self, n_max: int, *, rank_tol: float = 1e-12) -> None:
        if n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {n_max}")
        self.n_max = n_max
        self.basis = hilbert_basis(n_max)
        self.dim_H = len(self.basis)
        self._rank_tol = rank_tol

        all_labels = allowed_multiplier_labels(n_max)
        labels = []
        matrices = []
        for (N, L, M) in all_labels:
            mat = build_multiplier_matrix(N, L, M, self.basis)
            if np.allclose(mat, 0, atol=1e-15):
                continue
            labels.append((N, L, M))
            matrices.append(mat)

        self.multiplier_labels = labels
        self.multiplier_matrices = matrices
        self.basis_matrices = list(zip(labels, matrices))

        # Cache the vec form for span / rank computations.
        self._vec_cache = self._build_vec_matrix(matrices)
        self._dim_cache = int(np.linalg.matrix_rank(self._vec_cache, tol=rank_tol))

    @property
    def dim(self) -> int:
        """Dimension of O as a complex linear subspace of M_{dim_H}(C)."""
        return self._dim_cache

    @property
    def envelope_dim(self) -> int:
        """Dimension of the C*-envelope = M_N(C) = N^2."""
        return self.dim_H ** 2

    # -----------------------------------------------------------------
    # Internal helpers
    # -----------------------------------------------------------------

    @staticmethod
    def _build_vec_matrix(matrices: Sequence[np.ndarray]) -> np.ndarray:
        """Stack the vec-forms of a list of matrices as columns of a tall
        matrix. Returns a (N^2, K) complex array, where K = len(matrices).
        """
        if len(matrices) == 0:
            return np.zeros((0, 0), dtype=np.complex128)
        dim = matrices[0].shape[0]
        K = len(matrices)
        cols = np.zeros((dim * dim, K), dtype=np.complex128)
        for k, M in enumerate(matrices):
            cols[:, k] = M.reshape(-1)
        return cols

    @staticmethod
    def _real_rank(complex_vec_matrix: np.ndarray, tol: float) -> int:
        """Treat the complex vec-form as a real matrix (real and imag parts
        as separate columns) and compute its rank as a *complex* span, i.e.
        rank over C of the original.

        np.linalg.matrix_rank works on complex matrices natively, so we just
        delegate.
        """
        if complex_vec_matrix.size == 0:
            return 0
        return int(np.linalg.matrix_rank(complex_vec_matrix, tol=tol))

    # -----------------------------------------------------------------
    # Membership / projection test
    # -----------------------------------------------------------------

    def contains(self, matrix: np.ndarray, *, tol: float = 1e-10) -> Tuple[bool, float]:
        """Test whether `matrix` lies in the linear span of self.multiplier_matrices.

        Returns (is_in_O, residual_ratio), where residual_ratio is
        ||matrix - proj_O(matrix)|| / ||matrix|| in Frobenius norm.

        Method: solve a least-squares problem c = argmin ||V c - vec(matrix)||
        where V is the vec-stack, then compute the residual.
        """
        if matrix.shape != (self.dim_H, self.dim_H):
            raise ValueError(
                f"matrix shape {matrix.shape} != ({self.dim_H}, {self.dim_H})"
            )
        target = matrix.reshape(-1).astype(np.complex128)
        norm_target = np.linalg.norm(target)
        if norm_target < 1e-30:
            return True, 0.0
        V = self._vec_cache
        # least-squares solve V c = target
        c, *_ = np.linalg.lstsq(V, target, rcond=None)
        residual = np.linalg.norm(V @ c - target)
        ratio = float(residual / norm_target)
        return ratio < tol, ratio

    # -----------------------------------------------------------------
    # Self-checks
    # -----------------------------------------------------------------

    def identity_in_O(self, *, tol: float = 1e-10) -> Tuple[bool, float]:
        """Verify that the identity matrix lies in O.

        The identity comes from the multiplier f = constant 1, which is
        proportional to Y^{(3)}_{1, 0, 0} (the unique L=N-1=0 mode at N=1).
        Equivalently, the L=N=M=0 mode is the constant function on S^3 (up to
        normalization), and its matrix in the orthonormal basis is the
        identity (up to overall scale).
        """
        I = np.eye(self.dim_H, dtype=np.complex128)
        return self.contains(I, tol=tol)

    def is_star_closed(self, *, tol: float = 1e-10) -> Tuple[bool, List[Tuple[int, float]]]:
        """Verify that for each generator M_i, the conjugate transpose M_i^*
        lies in O.

        Returns (all_closed, list of (i, residual) for any failures).
        """
        failures = []
        for i, M in enumerate(self.multiplier_matrices):
            ok, residual = self.contains(M.conj().T, tol=tol)
            if not ok:
                failures.append((i, residual))
        return (len(failures) == 0, failures)


# ---------------------------------------------------------------------------
# Operator system products O^k and propagation number
# ---------------------------------------------------------------------------


def operator_system_products(
    matrices: Sequence[np.ndarray], k: int,
) -> List[np.ndarray]:
    """Compute the linear-span generators of O^k = span{ M_{i1} M_{i2} ... M_{ik}
    : 1 <= ij <= K } as a list of matrices.

    The list is the cartesian product of length k of indices into `matrices`.
    Note: this is in general redundant (many products will be linearly
    dependent), but rank computation handles that.

    For k=1, returns matrices itself.
    """
    if k < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    if k == 1:
        return list(matrices)
    products: List[np.ndarray] = []
    # Iteratively build up: start with O^1, multiply by O on the right each step.
    current = list(matrices)
    for _ in range(k - 1):
        next_level = []
        for A in current:
            for B in matrices:
                next_level.append(A @ B)
        current = next_level
    return current


def operator_system_dim(
    matrices: Sequence[np.ndarray], *, tol: float = 1e-10,
) -> int:
    """Dimension of the linear span of `matrices` as a complex subspace of
    M_N(C). Uses np.linalg.matrix_rank on the vec-stack.
    """
    if len(matrices) == 0:
        return 0
    dim = matrices[0].shape[0]
    K = len(matrices)
    cols = np.zeros((dim * dim, K), dtype=np.complex128)
    for k, M in enumerate(matrices):
        cols[:, k] = M.reshape(-1)
    return int(np.linalg.matrix_rank(cols, tol=tol))


def propagation_number(
    O: TruncatedOperatorSystem,
    *,
    max_k: int = 10,
    tol: float = 1e-10,
    verbose: bool = False,
) -> Tuple[int, List[int]]:
    """Compute prop(O) = smallest k such that dim(O^k) = N^2.

    Returns (prop, dim_sequence) where dim_sequence[i] = dim(O^{i+1}).

    Strategy: compute O^1, O^2, ... incrementally, taking advantage of the
    fact that O^{k+1} can be generated by multiplying a basis of O^k by O on
    the right. This avoids the combinatorial blowup of K^k products for
    large k by deduplicating at each step.

    The cap max_k is a safety bound. If prop > max_k, returns (-1, dims) with
    dim_sequence reported up to max_k.
    """
    N = O.dim_H
    target_dim = N * N
    matrices_O = O.multiplier_matrices
    if not matrices_O:
        return -1, []

    dim_sequence: List[int] = []

    # Track a basis of O^k as a list of matrices (not deduplicated to a true
    # linear-span basis; we just track enough to generate O^{k+1}).
    # We use a more efficient strategy: always work with the list of
    # generators of O^k as { O^{k-1} basis } @ { O basis }, then compute
    # rank.

    # First level: O^1
    current_generators = list(matrices_O)
    dim_O1 = operator_system_dim(current_generators, tol=tol)
    dim_sequence.append(dim_O1)
    if verbose:
        print(f"  k=1: |gens| = {len(current_generators)}, dim(O^1) = {dim_O1}, target = {target_dim}")
    if dim_O1 == target_dim:
        return 1, dim_sequence

    # For k >= 2: at each step, generators of O^{k+1} are { g @ M : g in O^k_generators, M in O }.
    # To keep memory bounded we extract a basis of O^k as matrices via SVD
    # then multiply by all matrices_O.

    for k in range(2, max_k + 1):
        # Extract a basis of O^k from current_generators (compress).
        basis_Ok = _extract_matrix_basis(current_generators, tol=tol)
        # Now generate O^{k+1} = basis_Ok @ matrices_O.
        new_generators = []
        for A in basis_Ok:
            for B in matrices_O:
                new_generators.append(A @ B)
        dim_Ok_plus_1 = operator_system_dim(new_generators, tol=tol)
        dim_sequence.append(dim_Ok_plus_1)
        if verbose:
            print(
                f"  k={k}: |basis(O^{k-1})| = {len(basis_Ok)}, |gens| = {len(new_generators)}, "
                f"dim(O^{k}) = {dim_Ok_plus_1}, target = {target_dim}"
            )
        if dim_Ok_plus_1 == target_dim:
            return k, dim_sequence
        if dim_Ok_plus_1 == dim_sequence[-2]:
            # Saturated below the target; prop = infinity (or at least > max_k).
            if verbose:
                print(
                    f"  saturated at dim {dim_Ok_plus_1} < {target_dim}; "
                    f"prop = infinity (or > {max_k})"
                )
            return -1, dim_sequence
        current_generators = new_generators

    return -1, dim_sequence


def _extract_matrix_basis(
    matrices: Sequence[np.ndarray], *, tol: float = 1e-10,
) -> List[np.ndarray]:
    """Extract a complex-linear basis of span{matrices} as a list of
    matrices, via SVD on the vec-stack.
    """
    if len(matrices) == 0:
        return []
    dim = matrices[0].shape[0]
    K = len(matrices)
    cols = np.zeros((dim * dim, K), dtype=np.complex128)
    for i, M in enumerate(matrices):
        cols[:, i] = M.reshape(-1)
    # SVD: cols = U S Vh; rank = number of singular values > tol * max(S).
    U, S, _ = np.linalg.svd(cols, full_matrices=False)
    if S.size == 0:
        return []
    rank = int(np.sum(S > tol * max(S.max(), 1.0)))
    # Take the first `rank` columns of U as an orthonormal basis of the
    # column space, reshape to matrices.
    basis_vecs = U[:, :rank]
    return [basis_vecs[:, k].reshape((dim, dim)) for k in range(rank)]


# ---------------------------------------------------------------------------
# Witness pair construction
# ---------------------------------------------------------------------------


def witness_pair(
    O: TruncatedOperatorSystem,
    *,
    N_target: int = 2, L_target: int = 1, M_target: int = 0,
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], Optional[np.ndarray], Optional[Tuple[bool, float]]]:
    """Construct the witness pair (a, b) = (M_{N L M}, M_{N L M}^*) and test
    whether ab is in O.

    Returns (a, b, ab, (in_O, residual)). If the requested generator does not
    exist, returns (None, None, None, None).
    """
    label_target = (N_target, L_target, M_target)
    a = None
    for label, mat in O.basis_matrices:
        if label == label_target:
            a = mat
            break
    if a is None:
        return None, None, None, None
    b = a.conj().T
    ab = a @ b
    test = O.contains(ab)
    return a, b, ab, test
