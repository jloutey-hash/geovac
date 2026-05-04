"""Spinor-lifted truncated operator system on the Weyl sector of S^3.

Sprint WH1-R3.2: lift the Connes-vS truncated operator system from the
scalar Fock basis (geovac/operator_system.py) to the Weyl-spinor sector
of the Camporesi-Higuchi spinor bundle on S^3.

Mathematical setup
==================

Spinor basis (Bar 1996, Theorem 1; Camporesi-Higuchi 1996)
-----------------------------------------------------------

The Weyl sector of the spinor bundle on S^3 at level n_ch = n_fock - 1
decomposes under the diagonal SU(2) into total-angular-momentum irreps
with j = 1/2, 3/2, ..., n_ch + 1/2. We use the j = l + 1/2 chain (single
channel per l) with l = 0, ..., n_ch and m_j running -j, ..., +j in
half-integer steps. Per-level dimension:

    dim(level n_ch) = sum_{l=0}^{n_ch} (2l + 2) = (n_ch+1)(n_ch+2).

Cumulative: dim(n_max=1) = 2, dim(n_max=2) = 8, dim(n_max=3) = 20,
dim(n_max=4) = 40 (Fock convention, n_max = n_fock_max).

The basis label is

    SpinorLabel(n_fock, l, two_m_j)

with n_fock in {1, ..., n_max}, l in {0, ..., n_fock - 1}, and two_m_j
an odd integer in {-(2l+1), -(2l-1), ..., +(2l+1)} so that
m_j = two_m_j / 2 is in {-l - 1/2, -l + 1/2, ..., l + 1/2}.

Clebsch-Gordan decomposition (j = l + 1/2 chain)
------------------------------------------------

Each spinor harmonic decomposes into the scalar Y_{l,m} ⊗ |sigma> basis:

    |ψ_{l, j=l+1/2, m_j}> = α_+ |Y_{l, m_j - 1/2}> ⊗ |+1/2>
                          + α_- |Y_{l, m_j + 1/2}> ⊗ |-1/2>

with the standard CG coefficients (real, positive)

    α_+ = sqrt((l + m_j + 1/2) / (2l + 1))
    α_- = sqrt((l - m_j + 1/2) / (2l + 1))

Edge cases:
- m_j = +l + 1/2 (max): α_+ = 1, α_- = 0 (only spin-up, m_l = l).
- m_j = -l - 1/2 (min): α_+ = 0, α_- = 1 (only spin-down, m_l = -l).

When the formula would produce α_± involving Y_{l, m_l} with |m_l| > l,
the corresponding alpha vanishes by the formula (sqrt of zero), so no
manual range-check is needed.

Scalar multiplier matrix elements on the spinor bundle
------------------------------------------------------

A scalar function f(omega) on S^3 acts trivially on the spin index. With
f = sum_{NLM} c_{NLM} Y^{(3)}_{NLM}, the matrix element between spinor
states is

    <ψ_a | M_{NLM} | ψ_b>
        = α_+^a α_+^b <Y_{l_a, m_j^a-1/2} | M | Y_{l_b, m_j^b-1/2}>
        + α_-^a α_-^b <Y_{l_a, m_j^a+1/2} | M | Y_{l_b, m_j^b+1/2}>

i.e. a sum of TWO scalar 3-Y integrals weighted by products of CG
coefficients. The scalar 3-Y integrals are exactly the
geovac.operator_system.three_y_integral results computed via the Avery-
Wen-Avery 3-Y integral on S^3 (sprint R3.1).

Native Camporesi-Higuchi Dirac on the Weyl sector
-------------------------------------------------

The Camporesi-Higuchi Dirac operator on the unit S^3 acts diagonally
in the (n_fock, l, m_j) basis with eigenvalue (Camporesi-Higuchi 1996
Eq. 4.1):

    D~ |ψ_{n_fock, l, m_j}> = (n_fock + 1/2) |ψ_{n_fock, l, m_j}>

(equivalently (n_ch + 3/2) in CH convention; same value by n_fock = n_ch + 1
 minus 1 / 2 + 1).

The full Dirac sector would also include the (-1)-chirality block with
the opposite-sign eigenvalue; this module restricts to the +chirality
Weyl sector for Connes-distance computation.

Why the scalar multipliers do NOT all commute with D~
-----------------------------------------------------

The diagonal Dirac D~ in the (n, l, m_j) basis would naively have a
large [D, .]-kernel: any multiplier diagonal in (n, l, m_j) commutes
with it. The scalar M^{scalar}_{N=1, L=0, M=0} multiplier is the
constant function on S^3 (= identity); but for N >= 2 the M_{N, L=0, M=0}
multipliers are AXISYMMETRIC SCALARS that change n while preserving
(l, m). Even though they preserve the (l, m) labels, they have nonzero
shell-coupling entries (n != n'), and these break the commutation with
D~ (which is diagonal in n). The kernel of [D~, .] is therefore just
the identity multiples (analogous to the round-3 'offdiag' Dirac with
alpha = 1, which had ker = C * I), and the SDP gives finite distances
on every pair (modulo the m_j-reflection forced zeros from operator-
system SO(3) symmetry).

API compatibility
=================

SpinorTruncatedOperatorSystem exposes the same attribute surface as
TruncatedOperatorSystem:

    .basis              list of SpinorLabel
    .dim_H              int = sum_{n=1..n_max} n*(n+1)
    .multiplier_matrices  list of complex np.ndarray of shape (dim_H, dim_H)
    .multiplier_labels    list of (N, L, M)
    .dim                int (linear-span dim)
    .envelope_dim       int = dim_H ** 2
    .contains(matrix, tol=...)
    .identity_in_O(tol=...)
    .is_star_closed(tol=...)

so it drops directly into geovac.connes_distance.compute_distance_matrix.

References
==========

C. Bär, "The Dirac operator on space forms of positive curvature",
J. Math. Soc. Japan 48 (1996) 69-83.

R. Camporesi & A. Higuchi, "On the eigenfunctions of the Dirac operator
on spheres and real hyperbolic spaces", J. Geom. Phys. 20 (1996) 1-18.

GeoVac sprint records:
- WH1-R3.1 memo: debug/wh1_r31_avery_wen_avery_memo.md
- Dirac-on-S^3 Tier 1 D1 module: geovac/dirac_s3.py
- Operator system: geovac/operator_system.py
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np
import sympy as sp
from sympy import Rational, simplify, sqrt

from geovac.operator_system import (
    HyperLabel,
    allowed_multiplier_labels,
    three_y_integral,
)


# ---------------------------------------------------------------------------
# Basis indexing
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SpinorLabel:
    """Label for a Weyl-spinor harmonic on unit S^3.

    Fields
    ------
    n_fock : int
        Fock principal quantum number, n_fock >= 1. CH index n_ch = n_fock - 1.
    l : int
        Orbital angular momentum of the embedded scalar harmonic factor;
        0 <= l <= n_fock - 1. Total angular momentum j = l + 1/2.
    two_m_j : int
        2 * m_j, an odd integer in {-(2l+1), -(2l-1), ..., +(2l+1)}.
        m_j = two_m_j / 2 is half-integer in {-l - 1/2, ..., +l + 1/2}.

    Conventions
    -----------
    Single-chirality (Weyl) sector. j = l + 1/2 chain (no j = l - 1/2
    states; the j = l - 1/2 chain is the other chirality and is omitted
    here). This gives (n_ch+1)(n_ch+2) states per level.
    """

    n_fock: int
    l: int
    two_m_j: int

    def __post_init__(self) -> None:  # pragma: no cover - sanity
        if self.n_fock < 1:
            raise ValueError(f"n_fock must be >= 1, got {self.n_fock}")
        if not (0 <= self.l <= self.n_fock - 1):
            raise ValueError(
                f"l must be in [0, {self.n_fock - 1}], got {self.l}"
            )
        # two_m_j must be odd integer, |two_m_j| <= 2l + 1
        if self.two_m_j % 2 == 0:
            raise ValueError(f"two_m_j must be odd, got {self.two_m_j}")
        if abs(self.two_m_j) > 2 * self.l + 1:
            raise ValueError(
                f"|two_m_j| must be <= {2*self.l + 1}, got {self.two_m_j}"
            )

    @property
    def m_j(self) -> Rational:
        """Half-integer m_j = two_m_j / 2."""
        return Rational(self.two_m_j, 2)

    @property
    def j(self) -> Rational:
        """Total angular momentum j = l + 1/2."""
        return Rational(2 * self.l + 1, 2)


def spinor_basis(n_max: int) -> List[SpinorLabel]:
    """All Weyl-spinor labels on S^3 up to n_fock <= n_max.

    Order: by n_fock first, then l, then two_m_j (ascending).

    Total dim: sum_{n=1..n_max} n*(n+1).
    """
    if n_max < 1:
        raise ValueError(f"n_max must be >= 1, got {n_max}")
    basis = []
    for n_fock in range(1, n_max + 1):
        for l in range(n_fock):  # 0 <= l <= n_fock - 1
            for two_m_j in range(-(2 * l + 1), 2 * l + 2, 2):
                basis.append(SpinorLabel(n_fock=n_fock, l=l, two_m_j=two_m_j))
    return basis


def spinor_dim(n_max: int) -> int:
    """sum_{n=1..n_max} n*(n+1) = (n_max)(n_max+1)(n_max+2)/3 (Weyl sector)."""
    return sum(n * (n + 1) for n in range(1, n_max + 1))


# ---------------------------------------------------------------------------
# Clebsch-Gordan coefficients for j = l + 1/2 chain
# ---------------------------------------------------------------------------


def cg_coefficients(label: SpinorLabel) -> Tuple[Rational, Rational, int, int]:
    """Return (alpha_plus, alpha_minus, m_l_plus, m_l_minus) for a spinor label.

    The spinor harmonic decomposes as

        |ψ_{l, j=l+1/2, m_j}> = alpha_plus |Y_{l, m_l_plus}> |+1/2>
                              + alpha_minus |Y_{l, m_l_minus}> |-1/2>

    with m_l_plus = m_j - 1/2 and m_l_minus = m_j + 1/2 (both integers).

    CG coefficients (real, sympy-exact):
        alpha_plus  = sqrt((l + m_j + 1/2) / (2l + 1))
        alpha_minus = sqrt((l - m_j + 1/2) / (2l + 1))

    Returns Rationals (alpha_plus^2, alpha_minus^2 are rational) for the
    *squared* alpha values, plus the integer m_l indices. We return the
    sqrt-rational alpha values directly as sympy.Expr.

    Edge cases handled by the formula:
    - m_j = +l + 1/2 (max): alpha_plus = 1, alpha_minus = 0.
    - m_j = -l - 1/2 (min): alpha_plus = 0, alpha_minus = 1.
    """
    l = label.l
    two_m_j = label.two_m_j
    # alpha_plus^2 = (l + m_j + 1/2) / (2l + 1) = (2l + two_m_j + 1) / (2(2l + 1))
    plus_num = 2 * l + two_m_j + 1
    minus_num = 2 * l - two_m_j + 1
    denom = 2 * (2 * l + 1)
    alpha_plus_sq = Rational(plus_num, denom)
    alpha_minus_sq = Rational(minus_num, denom)
    # Sanity: sum to 1
    assert alpha_plus_sq + alpha_minus_sq == 1, (
        f"alpha_+^2 + alpha_-^2 != 1 for {label}: "
        f"{alpha_plus_sq} + {alpha_minus_sq}"
    )
    # m_l for sigma = +1/2: m_l = m_j - 1/2 = (two_m_j - 1) / 2 must be integer
    if (two_m_j - 1) % 2 != 0:
        raise AssertionError(f"non-integer m_l_plus for {label}")
    m_l_plus = (two_m_j - 1) // 2
    if (two_m_j + 1) % 2 != 0:
        raise AssertionError(f"non-integer m_l_minus for {label}")
    m_l_minus = (two_m_j + 1) // 2
    return (sqrt(alpha_plus_sq), sqrt(alpha_minus_sq), m_l_plus, m_l_minus)


# ---------------------------------------------------------------------------
# Spinor multiplier matrix builder
# ---------------------------------------------------------------------------


def build_spinor_multiplier_matrix(
    N: int, L: int, M: int, basis: Sequence[SpinorLabel],
) -> np.ndarray:
    """Construct M~_{NLM} for a scalar multiplier on the Weyl spinor bundle.

    Matrix element formula:

        <ψ_a | M_{NLM} | ψ_b>
            = alpha_+^a alpha_+^b <Y_{l_a, m_l_plus^a} | M | Y_{l_b, m_l_plus^b}>
            + alpha_-^a alpha_-^b <Y_{l_a, m_l_minus^a} | M | Y_{l_b, m_l_minus^b}>

    where m_l_plus = m_j - 1/2 and m_l_minus = m_j + 1/2, and the scalar
    matrix elements are the geovac.operator_system.three_y_integral values
    (computed via the Avery-Wen-Avery 3-Y integral on S^3).

    Returns a complex numpy array of shape (dim, dim), dtype complex128.
    """
    dim = len(basis)
    mat = np.zeros((dim, dim), dtype=np.complex128)

    # Pre-compute CG data for every basis label.
    cg_data = [cg_coefficients(b) for b in basis]

    for i, bra in enumerate(basis):
        ap_a, am_a, mlp_a, mlm_a = cg_data[i]
        for j, ket in enumerate(basis):
            ap_b, am_b, mlp_b, mlm_b = cg_data[j]
            total = sp.Integer(0)
            # sigma = +1/2 contribution
            if abs(mlp_a) <= bra.l and abs(mlp_b) <= ket.l:
                if ap_a != 0 and ap_b != 0:
                    s_val = three_y_integral(
                        bra.n_fock, bra.l, mlp_a,
                        N, L, M,
                        ket.n_fock, ket.l, mlp_b,
                    )
                    if s_val != 0:
                        total = total + ap_a * ap_b * s_val
            # sigma = -1/2 contribution
            if abs(mlm_a) <= bra.l and abs(mlm_b) <= ket.l:
                if am_a != 0 and am_b != 0:
                    s_val = three_y_integral(
                        bra.n_fock, bra.l, mlm_a,
                        N, L, M,
                        ket.n_fock, ket.l, mlm_b,
                    )
                    if s_val != 0:
                        total = total + am_a * am_b * s_val
            if total == 0:
                continue
            mat[i, j] = complex(sp.N(total, 30))
    return mat


# ---------------------------------------------------------------------------
# Native Camporesi-Higuchi Dirac matrix
# ---------------------------------------------------------------------------


def camporesi_higuchi_dirac_matrix(basis: Sequence[SpinorLabel]) -> np.ndarray:
    """Native Camporesi-Higuchi Dirac D~ on the Weyl sector at unit S^3.

    Diagonal in (n_fock, l, m_j) basis with eigenvalue n_fock + 1/2
    (equivalently n_ch + 3/2 in CH convention). Real, symmetric.

    Returns a complex numpy array of shape (dim, dim).
    """
    dim = len(basis)
    diag = np.array(
        [float(b.n_fock) + 0.5 for b in basis], dtype=np.complex128,
    )
    return np.diag(diag)


# ---------------------------------------------------------------------------
# SpinorTruncatedOperatorSystem
# ---------------------------------------------------------------------------


class SpinorTruncatedOperatorSystem:
    """Spinor-lifted truncated operator system on the Weyl sector of S^3.

    API-compatible with TruncatedOperatorSystem so it drops into
    geovac.connes_distance.compute_connes_distance unchanged.

    Attributes
    ----------
    n_max : int
    basis : list of SpinorLabel, len = dim_H
    dim_H : int
    multiplier_labels : list of (N, L, M)
    multiplier_matrices : list of complex np.ndarray
    basis_matrices : list of (label, matrix) pairs
    dim : int (linear-span dim of O over C in M_{dim_H}(C))
    envelope_dim : int = dim_H ** 2

    Notes
    -----
    The multiplier_labels are inherited from the scalar
    allowed_multiplier_labels (same SO(4) selection rules); the matrices
    are the spinor-lifted versions. Some multipliers that produce
    nonzero matrices on the scalar bundle may produce all-zero matrices
    on the spinor bundle (because the CG-weighted sum can cancel); we
    filter those out at construction time.
    """

    def __init__(self, n_max: int, *, rank_tol: float = 1e-12) -> None:
        if n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {n_max}")
        self.n_max = n_max
        self.basis = spinor_basis(n_max)
        self.dim_H = len(self.basis)
        self._rank_tol = rank_tol

        all_labels = allowed_multiplier_labels(n_max)
        labels = []
        matrices = []
        for (N, L, M) in all_labels:
            mat = build_spinor_multiplier_matrix(N, L, M, self.basis)
            if np.allclose(mat, 0, atol=1e-15):
                continue
            labels.append((N, L, M))
            matrices.append(mat)

        self.multiplier_labels = labels
        self.multiplier_matrices = matrices
        self.basis_matrices = list(zip(labels, matrices))

        # Cache vec form for span / rank computations.
        self._vec_cache = self._build_vec_matrix(matrices)
        self._dim_cache = int(
            np.linalg.matrix_rank(self._vec_cache, tol=rank_tol)
        )

    @property
    def dim(self) -> int:
        return self._dim_cache

    @property
    def envelope_dim(self) -> int:
        return self.dim_H ** 2

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

    def contains(
        self, matrix: np.ndarray, *, tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Test whether `matrix` lies in span(self.multiplier_matrices)."""
        if matrix.shape != (self.dim_H, self.dim_H):
            raise ValueError(
                f"matrix shape {matrix.shape} != ({self.dim_H}, {self.dim_H})"
            )
        target = matrix.reshape(-1).astype(np.complex128)
        norm_target = np.linalg.norm(target)
        if norm_target < 1e-30:
            return True, 0.0
        V = self._vec_cache
        c, *_ = np.linalg.lstsq(V, target, rcond=None)
        residual = np.linalg.norm(V @ c - target)
        ratio = float(residual / norm_target)
        return ratio < tol, ratio

    def identity_in_O(
        self, *, tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Verify that I_{dim_H} lies in O.

        The constant function f = 1 on S^3 corresponds to the M_{1, 0, 0}
        multiplier (the Y^{(3)}_{1, 0, 0} mode). On the spinor bundle this
        acts as a scalar on each spinor component, so the lifted matrix
        is proportional to the identity.
        """
        I = np.eye(self.dim_H, dtype=np.complex128)
        return self.contains(I, tol=tol)

    def is_star_closed(
        self, *, tol: float = 1e-10,
    ) -> Tuple[bool, List[Tuple[int, float]]]:
        """Verify that for each generator M_i, M_i^* is in O."""
        failures = []
        for i, M in enumerate(self.multiplier_matrices):
            ok, residual = self.contains(M.conj().T, tol=tol)
            if not ok:
                failures.append((i, residual))
        return (len(failures) == 0, failures)


# ---------------------------------------------------------------------------
# Convenience: spinor-graph distance and pretty labels
# ---------------------------------------------------------------------------


def spinor_graph_distance(v: SpinorLabel, w: SpinorLabel) -> int:
    """Combinatorial graph distance on the Weyl-spinor labels.

    Use the L^1 metric on (n_fock, l, two_m_j / 2):

        d(v, w) = |Δn_fock| + |Δl| + |Δm_j|
                = |Δn_fock| + |Δl| + |Δtwo_m_j| / 2

    where Δm_j is half-integer (so |Δm_j| = |Δtwo_m_j| / 2 always integer
    for spinor-spinor pairs since two_m_j parities match within the same
    chirality).
    """
    return (
        abs(v.n_fock - w.n_fock)
        + abs(v.l - w.l)
        + abs(v.two_m_j - w.two_m_j) // 2
    )


def spinor_graph_distance_matrix(
    op_sys: "SpinorTruncatedOperatorSystem",
) -> np.ndarray:
    """N x N integer graph-distance matrix on the spinor basis."""
    N = op_sys.dim_H
    dist = np.zeros((N, N), dtype=int)
    for i in range(N):
        for j in range(N):
            dist[i, j] = spinor_graph_distance(op_sys.basis[i], op_sys.basis[j])
    return dist


def spinor_label_strings(op_sys: "SpinorTruncatedOperatorSystem") -> List[str]:
    """Pretty label strings |n, l, m_j> using m_j as a fraction."""
    out = []
    for b in op_sys.basis:
        # Format m_j as an integer or +/- 1/2, +/- 3/2 etc.
        if b.two_m_j % 2 == 0:
            mj_str = f"{b.two_m_j // 2}"
        else:
            sign = "+" if b.two_m_j > 0 else "-"
            mj_str = f"{sign}{abs(b.two_m_j)}/2"
        out.append(f"|{b.n_fock},{b.l},{mj_str}>")
    return out
