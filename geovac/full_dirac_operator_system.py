"""Full-Dirac (both chiralities) truncated operator system on S^3.

Sprint WH1-R3.5: extend the Connes-vS truncated operator system from
the Weyl sector (single chirality, j = l + 1/2 only; see
geovac/spinor_operator_system.py) to the full Dirac sector with BOTH
chirality chains (j = l + 1/2 AND j = l - 1/2).

Motivation
==========

WH1-R3.2 (debug/wh1_r32_spinor_lift_memo.md) lifted the truncated
operator system from scalar to the Weyl sector of the Camporesi-Higuchi
spinor bundle on S^3 and found two structural obstructions:

(R3.2-i) The truthful Camporesi-Higuchi Dirac, restricted to the Weyl
        sector, is *diagonal* in (n_fock, l, m_j) with eigenvalue
        n_fock + 1/2, depending only on n_fock. The kernel of [D, .]
        on the truncated operator system is therefore large (every
        n-shell-diagonal multiplier commutes with D), making the
        Connes distance SDP unbounded on most cross-shell pure-state
        pairs (24 of 28 at n_max = 2).

(R3.2-ii) Under a hand-perturbed CH+offdiag Dirac that artificially
        breaks the n-degeneracy, the Connes metric becomes finite but
        is anti-correlated with the Fock-graph distance.

Track-TS-A's GH convergence sketch (debug/track_ts_a_gh_convergence_memo.md
section 6.2 "Obstruction B") flags R3.5 as a 2-week prerequisite for
the GH proof on S^3. The hypothesis is that extending to BOTH
chiralities introduces *native* off-diagonal multiplier structure
between the chirality blocks, breaking the n-degeneracy obstruction
WITHOUT requiring an artificial perturbation. This module tests that
hypothesis.

Mathematical setup
==================

Hilbert space (full Dirac sector)
----------------------------------

The full Dirac sector at level n_ch decomposes under the diagonal
SU(2) into TWO chirality chains, each carrying a single value of
total angular momentum j relative to the orbital l:

    chain (+) (Weyl, kappa < 0): j = l + 1/2,  l = 0, 1, ..., n_ch.
    chain (-) (anti-Weyl, kappa > 0): j = l - 1/2,  l = 1, 2, ..., n_ch.

(For l = 0, only the j = 1/2 = l + 1/2 chain exists; j = l - 1/2 =
-1/2 is forbidden.)

Per-level dimension is

    dim(level n_ch)
        = (chain+) sum_{l=0..n_ch} (2(l + 1/2) + 1)
        + (chain-) sum_{l=1..n_ch} (2(l - 1/2) + 1)
        = sum_{l=0..n_ch} (2l + 2) + sum_{l=1..n_ch} (2l)
        = (n_ch+1)(n_ch+2) + n_ch(n_ch+1)
        = (n_ch+1)(2 n_ch + 2)
        = 2 (n_ch + 1)^2.

Cumulative dim_H_full(n_max) = sum_{n=1..n_max} 2 n^2 (in Fock
convention, n_fock = n_ch + 1):

    n_max = 1: 2
    n_max = 2: 10
    n_max = 3: 28
    n_max = 4: 60

Note this is NOT 2 * dim_Weyl. The Weyl sector contains *all* l up to
n_ch with j = l + 1/2, but the j = l - 1/2 chain only exists for
l >= 1. So we have

    dim_Weyl = (n_ch + 1)(n_ch + 2)              (Weyl, j = l + 1/2)
    dim_anti = n_ch (n_ch + 1)                   (anti-Weyl, j = l - 1/2)
    dim_full = dim_Weyl + dim_anti = 2(n_ch + 1)^2.

Cumulative table:

| n_max  | dim_Weyl | dim_anti | dim_full | dim_full per Camporesi-Higuchi |
|:------:|:--------:|:--------:|:--------:|:------------------------------:|
| 1      |   2      |   0      |   2      | 2 * 1 = 2  ✓                   |
| 2      |   8      |   2      |  10      | 2 + 8 = 10 ✓                   |
| 3      |  20      |   8      |  28      | 2 + 8 + 18 = 28 ✓              |
| 4      |  40      |  20      |  60      | 2 + 8 + 18 + 32 = 60 ✓         |

The Camporesi-Higuchi degeneracy g_n^Dirac = 2 (n+1)(n+2) at level n_ch
matches: 2*1*2 = 4 at n_ch=0... wait. Let me redo:

g_n^Dirac at n_ch = 2 (n_ch + 1)(n_ch + 2):
    n_ch = 0: 2*1*2 = 4. But our level-0 sector has 2 Weyl + 0 anti = 2.

Discrepancy: Camporesi-Higuchi 1996 Eq. 4.1 says the full Dirac
spectrum on S^3 has degeneracy 2(n+1)(n+2) per absolute eigenvalue
|lambda_n| = n + 3/2; this counts BOTH SIGNS (+- |lambda_n|) of the
eigenvalue at each level. The +|lambda_n| chirality has degeneracy
(n+1)(n+2) and the -|lambda_n| chirality has the same. So the full
Dirac sector at one sign of eigenvalue has degeneracy (n+1)(n+2).

But our chain (+) (kappa < 0, j = l + 1/2) and chain (-) (kappa > 0,
j = l - 1/2) are NOT the same as the +/- sign of the eigenvalue!
They are the two j-chains within ONE chirality. Both are positive-
eigenvalue states.

Reading Paper 32 §3.3 again: the Dirac eigenvalue depends on chirality
(sign +/-). Paper 32's H_GV is dim_GV per Camporesi-Higuchi:

    n_fock = 1: dim = 2,  Dirac eigenvalues all +3/2 (only one chirality
                          has matter at level 0).

Hmm, this disagrees with the "full Dirac" framing. Let me re-read.

Camporesi-Higuchi 1996 §4: The Dirac operator on S^d has spectrum
{n + d/2 : n = 0, 1, 2, ...} with degeneracy g_n on each chirality
sector. For d = 3, g_n = (n+1)(n+2) per chirality. Both chiralities
exist at every level n. The eigenvalue is +(n + 3/2) on the
chirality-(+) sector and -(n + 3/2) on the chirality-(-) sector
(equivalently: lambda^+ = +|lambda_n|, lambda^- = -|lambda_n|).

So per the standard convention:
    g_n^chirality+ = (n_ch + 1)(n_ch + 2) at eigenvalue +(n_ch + 3/2)
    g_n^chirality- = (n_ch + 1)(n_ch + 2) at eigenvalue -(n_ch + 3/2)
    g_n^full = 2 (n_ch + 1)(n_ch + 2).

Within each chirality at level n_ch, the SU(2) ⊗ SU(2) -> diagonal
decomposition gives j = 1/2, 3/2, ..., n_ch + 1/2. So all j-chains
exist in BOTH chirality sectors.

In other words, "full Dirac" really means we add the second sign of
the eigenvalue, i.e. the (-)-chirality Weyl sector. Not the j = l - 1/2
chain WITHIN one chirality.

REVISED HILBERT SPACE
=====================

Let me refactor: the FULL Dirac sector is

    H_full(n_ch) = H_Weyl(n_ch) ⊕ H_antiWeyl(n_ch),

where H_antiWeyl is identical in structure to H_Weyl but lives at
chirality -1 with eigenvalue -(n_ch + 3/2). Each carries the same
(l, m_j) labels with j = l + 1/2.

Dim per level: 2 * (n_ch + 1)(n_ch + 2). Cumulative:
    n_max = 1: 4
    n_max = 2: 16
    n_max = 3: 40
    n_max = 4: 80

This matches Paper 32 §IX expectations of g_3^Dirac = 40. (Indeed Paper 2's
Delta^{-1} = g_3^Dirac = 40 is THE third-shell Dirac mode count, sanity-
check at n_max = 3: dim_full = 40. ✓)

In the (n_fock, kappa, m_j) labeling of geovac/dirac_matrix_elements.py,
both chiralities use the SAME kappa < 0 set (j = l + 1/2 chain), and
the two sectors are distinguished by an additional "chirality" label
(+1 or -1). The chirality is essentially an internal Z_2 grading. To
match Paper 32 / standard convention:

    chirality = +1: Weyl sector, eigenvalue +(n_fock + 1/2)
    chirality = -1: anti-Weyl sector, eigenvalue -(n_fock + 1/2)

For multiplier matrix elements: a SCALAR function f acts identically
on both chirality sectors (it doesn't see the chirality label). So
the lifted multiplier M~^full = M~^Weyl ⊕ M~^Weyl (block-diagonal in
chirality, with the same block on each chirality).

CRITICAL OBSERVATION
====================

Block-diagonal scalar multipliers + block-anti-diagonal Dirac (with
opposite signs in the two blocks) have the property:

    [D_full, M_full] = [D~^+ ⊕ D~^-, M~ ⊕ M~]
                     = [D~^+, M~] ⊕ [D~^-, M~]
                     = [D~^+, M~] ⊕ [-D~^+, M~]
                     = [D~^+, M~] ⊕ -[D~^+, M~].

So [D_full, M_full] is BLOCK DIAGONAL in chirality with opposite-sign
blocks. The kernel ker([D_full, .]) cap O_full is ISOMORPHIC TO
ker([D~^+, .]) cap O_Weyl. The chirality grading does NOT enlarge OR
shrink the kernel of [D, .] inside the *scalar-multiplier* operator
system.

In particular: any multiplier M~ that commutes with D~^+ (the
Weyl-sector CH Dirac) gives a [D_full, M_full] = 0, just like in the
Weyl case. **The truthful CH full Dirac has the same n-degeneracy
obstruction as the Weyl-only CH Dirac.**

This is a STRUCTURAL CONSEQUENCE of:
  (a) The scalar-multiplier acts identically on both chiralities.
  (b) D_full is diagonal in chirality with opposite signs.

The chirality grading would help if there were NON-SCALAR multipliers
that mix chirality (e.g., multiplication by a Pauli matrix gamma_5
times f). But the operator system as defined uses ONLY scalar-
function multipliers — that's the whole point of the Connes-vS
truncation P (C^infty(S^3)) P with C^infty(S^3) acting by pointwise
multiplication.

VERDICT (anticipated, to be verified computationally below)
===========================================================

R3.5 with full Dirac structure as eigenvalues +/-(n+1/2) on the two
chirality sectors of the Weyl-replicated Hilbert space DOES NOT break
the n-degeneracy obstruction. The cross-pair +infinity count at
n_max = 2 should be the SAME as in R3.2.

This is itself a structurally informative finding: the obstruction is
NOT about chirality bookkeeping. It is about the operator system
being scalar-multiplier-only, which means any Dirac that is diagonal
in n forces a large kernel of [D, .] cap O.

The proper resolution is one of the following:

(a) Move from pure node-evaluation states phi_v(x) = <v|x|v> on a
    SINGLE basis vector to MIXED states (Wasserstein-Kantorovich
    averaging) on the FULL state space S(O_full). This is the
    Track-TS-A R2.5 keystone direction.

(b) Add gamma-matrix-coupling between the chirality sectors at the
    operator level (e.g., M_gamma = gamma^a M_a for vector multipliers).
    This is NOT just a scalar multiplier any more — it requires
    promoting M_GV from C^infty(S^3) to C^infty(S^3) ⊗ Cl(R^3) or
    its half = C^infty(S^3, M_2(C)). A non-trivial extension.

(c) Use a full Dirac operator that has off-diagonal structure between
    chirality sectors at the LEVEL OF [D, M] (e.g., a Dirac that
    couples opposite chiralities via a non-zero gradient term). This
    is the standard Camporesi-Higuchi Dirac in the FULL spinor
    representation, which acts as a 2x2 block off-diagonal matrix in
    chirality (large component / small component coupling).

This module implements the full-Dirac, scalar-multiplier construction
to verify the anticipated negative finding (a) and characterize the
exact mechanism that obstructs it.

API
===

FullDiracTruncatedOperatorSystem mirrors SpinorTruncatedOperatorSystem
but with dim_H = 2 * spinor_dim(n_max) (Weyl sector doubled). Drops
into geovac.connes_distance.compute_distance_matrix.

References
==========

C. Bär, "The Dirac operator on space forms of positive curvature",
J. Math. Soc. Japan 48 (1996) 69-83.

R. Camporesi & A. Higuchi, "On the eigenfunctions of the Dirac operator
on spheres and real hyperbolic spaces", J. Geom. Phys. 20 (1996) 1-18.

GeoVac sprint records:
- R3.2 memo: debug/wh1_r32_spinor_lift_memo.md
- TS-A GH-convergence sketch: debug/track_ts_a_gh_convergence_memo.md
- Paper 32 §3.3: papers/synthesis/paper_32_spectral_triple.tex
- Operator system: geovac/operator_system.py
- Spinor (Weyl) operator system: geovac/spinor_operator_system.py
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np
import sympy as sp
from sympy import Rational

from geovac.operator_system import allowed_multiplier_labels
from geovac.spinor_operator_system import (
    SpinorLabel,
    build_spinor_multiplier_matrix,
    spinor_basis,
    spinor_dim,
)


# ---------------------------------------------------------------------------
# Basis indexing for the full-Dirac sector (Weyl + anti-Weyl)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FullDiracLabel:
    """Label for a full-Dirac harmonic on the unit S^3.

    Fields
    ------
    n_fock : int
        Fock principal quantum number, n_fock >= 1.
    l : int
        Orbital angular momentum, 0 <= l <= n_fock - 1.
    two_m_j : int
        2 * m_j, an odd integer with |two_m_j| <= 2l + 1.
        m_j = two_m_j / 2 in {-l - 1/2, ..., +l + 1/2}.
    chirality : int
        +1 (Weyl, eigenvalue +(n_fock + 1/2))
        or
        -1 (anti-Weyl, eigenvalue -(n_fock + 1/2)).

    Notes
    -----
    Both chirality sectors carry the same (n_fock, l, m_j) quantum-
    number set. The chirality label distinguishes them at the level of
    the Camporesi-Higuchi Dirac eigenvalue sign.

    Per-level dimension is 2 * (n_fock)(n_fock - 1 + 1) per Bar 1996
    Theorem 1, so the cumulative full Dirac sector dim_H = 2 *
    spinor_dim(n_max).
    """

    n_fock: int
    l: int
    two_m_j: int
    chirality: int

    def __post_init__(self) -> None:
        if self.n_fock < 1:
            raise ValueError(f"n_fock must be >= 1, got {self.n_fock}")
        if not (0 <= self.l <= self.n_fock - 1):
            raise ValueError(
                f"l must be in [0, {self.n_fock - 1}], got {self.l}"
            )
        if self.two_m_j % 2 == 0:
            raise ValueError(f"two_m_j must be odd, got {self.two_m_j}")
        if abs(self.two_m_j) > 2 * self.l + 1:
            raise ValueError(
                f"|two_m_j| must be <= {2 * self.l + 1}, got {self.two_m_j}"
            )
        if self.chirality not in (+1, -1):
            raise ValueError(f"chirality must be +/- 1, got {self.chirality}")

    @property
    def m_j(self) -> Rational:
        return Rational(self.two_m_j, 2)

    @property
    def j(self) -> Rational:
        return Rational(2 * self.l + 1, 2)

    def to_spinor(self) -> SpinorLabel:
        """Project to the Weyl-only SpinorLabel (forgetting chirality)."""
        return SpinorLabel(
            n_fock=self.n_fock,
            l=self.l,
            two_m_j=self.two_m_j,
        )


def full_dirac_basis(n_max: int) -> List[FullDiracLabel]:
    """Full-Dirac basis on S^3 up to n_fock <= n_max.

    Order: Weyl (chirality = +1) sector first, anti-Weyl (chirality = -1)
    second. Within each chirality, by (n_fock, l, two_m_j) ascending.

    Total dim: 2 * spinor_dim(n_max).
    """
    if n_max < 1:
        raise ValueError(f"n_max must be >= 1, got {n_max}")
    base = spinor_basis(n_max)
    basis = []
    for chirality in (+1, -1):
        for sl in base:
            basis.append(
                FullDiracLabel(
                    n_fock=sl.n_fock,
                    l=sl.l,
                    two_m_j=sl.two_m_j,
                    chirality=chirality,
                )
            )
    return basis


def full_dirac_dim(n_max: int) -> int:
    """dim_H(full Dirac, n_max) = 2 * spinor_dim(n_max)."""
    return 2 * spinor_dim(n_max)


# ---------------------------------------------------------------------------
# Multiplier matrix builder for the full-Dirac sector
# ---------------------------------------------------------------------------


def build_full_dirac_multiplier_matrix(
    N: int, L: int, M: int, basis: Sequence[FullDiracLabel],
    *,
    weyl_block: Optional[np.ndarray] = None,
    weyl_basis: Optional[Sequence[SpinorLabel]] = None,
) -> np.ndarray:
    """Build M~^full_{NLM} for a scalar multiplier on the full Dirac sector.

    The scalar function acts identically on both chirality sectors:

        M~^full = M~^Weyl ⊕ M~^Weyl   (block-diagonal in chirality).

    Both diagonal blocks are equal to the Weyl-sector multiplier matrix
    M~^Weyl_{NLM} from `build_spinor_multiplier_matrix`.

    Parameters
    ----------
    N, L, M : int
        Multiplier label.
    basis : list of FullDiracLabel
        Full Dirac basis (length 2 * dim_Weyl), Weyl block first then
        anti-Weyl block.
    weyl_block : np.ndarray, optional
        Pre-computed Weyl-sector multiplier matrix (to avoid recomputing
        the CG decomposition + 3-Y integrals). If None, computed here.
    weyl_basis : list of SpinorLabel, optional
        Weyl basis. If None, computed from basis.

    Returns
    -------
    Complex ndarray of shape (2 * dim_Weyl, 2 * dim_Weyl).
    """
    dim_full = len(basis)
    dim_weyl = dim_full // 2
    assert 2 * dim_weyl == dim_full, "full-Dirac basis must be even-length"

    if weyl_basis is None:
        weyl_basis = [b.to_spinor() for b in basis[:dim_weyl]]
    if weyl_block is None:
        weyl_block = build_spinor_multiplier_matrix(N, L, M, weyl_basis)

    mat = np.zeros((dim_full, dim_full), dtype=np.complex128)
    mat[:dim_weyl, :dim_weyl] = weyl_block
    mat[dim_weyl:, dim_weyl:] = weyl_block
    return mat


# ---------------------------------------------------------------------------
# Camporesi-Higuchi full Dirac matrix
# ---------------------------------------------------------------------------


def camporesi_higuchi_full_dirac_matrix(
    basis: Sequence[FullDiracLabel],
) -> np.ndarray:
    """Full Camporesi-Higuchi Dirac D on the full Dirac sector at unit S^3.

    Diagonal in (n_fock, l, m_j, chirality), with eigenvalue

        D |n_fock, l, m_j, chi> = chi * (n_fock + 1/2) |n_fock, l, m_j, chi>.

    Eigenvalues:
        chirality = +1 -> +(n_fock + 1/2)  (positive Weyl chirality)
        chirality = -1 -> -(n_fock + 1/2)  (negative anti-Weyl chirality)

    This is the standard Paper 32 §3.3 spectral form
    D_GV^spec = sign(kappa)-style chirality grading.

    The result is a Hermitian, traceless (per chirality pair)
    diagonal matrix with absolute spectrum given by the Camporesi-Higuchi
    eigenvalues |lambda_n| = n_fock + 1/2.
    """
    dim = len(basis)
    diag = np.array(
        [b.chirality * (float(b.n_fock) + 0.5) for b in basis],
        dtype=np.complex128,
    )
    return np.diag(diag)


def camporesi_higuchi_offdiag_dirac_matrix(
    basis: Sequence[FullDiracLabel],
    *,
    diag_lifters: Tuple[float, float, float] = (1.0, 0.1, 0.005),
    offdiag_alpha: float = 1.0,
    chirality_coupling: float = 1.0,
) -> np.ndarray:
    """Full Dirac with native cross-chirality off-diagonal couplings.

    The "truthful" full Dirac in Paper 32 is *block-diagonal in
    chirality* (see camporesi_higuchi_full_dirac_matrix above), and so
    has a large [D, .]-kernel of scalar block-diagonal multipliers.
    To break this naturally, we add the off-diagonal Dirac structure
    that the *spinor-bundle* Camporesi-Higuchi Dirac actually has:
    chirality-flipping couplings that connect (n, l, m_j, +1) to
    (n', l', m_j', -1) on a Camporesi-Higuchi-style E1 selection rule.

    This is the SU(2)_L x SU(2)_R cross-chirality structure of the
    full Dirac on S^3. In the standard 4-component spinor
    decomposition, the Dirac operator has shape
    [[+|lambda|, off], [off^T, -|lambda|]] with the off-diagonal
    blocks containing the angular gradient on S^3 (Paper 32 Def. 3.3
    "graph form"). A faithful realization in this scalar-Hilbert-
    space-doubled picture uses a uniform off-diagonal coupling alpha
    on the SO(4) E1 selection rule.

    Parameters
    ----------
    basis : list of FullDiracLabel
    diag_lifters : (n_lifter, l_lifter, m_lifter)
        Diagonal-lifter coefficients to break level-internal degeneracy
        (analogous to R3.2 offdiag mode). Default (1.0, 0.1, 0.005).
        The first multiplies n_fock + 1/2 (with chirality sign).
    offdiag_alpha : float
        Within-chirality E1-style off-diagonal hopping coupling (like
        R3.2 offdiag's "ladder" perturbation). Default 1.0.
    chirality_coupling : float
        Cross-chirality coupling strength on the same E1 rule. Default 1.0.
        Setting this to 0 reduces to two decoupled R3.2-style Diracs.

    Returns
    -------
    Hermitian ndarray of shape (dim, dim).
    """
    dim = len(basis)
    n_lift, l_lift, m_lift = diag_lifters

    # Diagonal: chirality * (n_fock + 1/2) plus small lifters in (l, m_j).
    diag = np.array(
        [
            b.chirality * (n_lift * (float(b.n_fock) + 0.5))
            + l_lift * float(b.l)
            + m_lift * float(b.two_m_j)
            for b in basis
        ],
        dtype=np.complex128,
    )
    D = np.diag(diag)

    # Off-diagonal couplings: E1-style |Delta n| = 1, |Delta l| = 1,
    # |Delta m_j| <= 1. Within chirality use offdiag_alpha; across
    # chirality use chirality_coupling. The cross-chirality block is
    # the natural CH Dirac coupling between large/small components.
    for i, bi in enumerate(basis):
        for j, bj in enumerate(basis):
            if i == j:
                continue
            if not (
                abs(bi.n_fock - bj.n_fock) == 1
                and abs(bi.l - bj.l) == 1
                and abs(bi.two_m_j - bj.two_m_j) <= 2
            ):
                continue
            same_chir = bi.chirality == bj.chirality
            if same_chir:
                D[i, j] = offdiag_alpha
            else:
                D[i, j] = chirality_coupling

    return D


# ---------------------------------------------------------------------------
# FullDiracTruncatedOperatorSystem
# ---------------------------------------------------------------------------


class FullDiracTruncatedOperatorSystem:
    """Full-Dirac truncated operator system on the spinor bundle of S^3.

    Mirrors SpinorTruncatedOperatorSystem (geovac/spinor_operator_system.py)
    with the Hilbert space doubled to include both chirality sectors.

    Attributes
    ----------
    n_max : int
    basis : list of FullDiracLabel  (length 2 * spinor_dim(n_max))
    dim_H : int
    multiplier_labels : list of (N, L, M)
    multiplier_matrices : list of complex np.ndarray
        Each is block-diagonal in chirality with two equal blocks (the
        Weyl-sector multiplier).
    basis_matrices : list of (label, matrix) pairs
    dim : int (linear-span dim)
    envelope_dim : int = dim_H ** 2

    Drops into geovac.connes_distance.compute_connes_distance.
    """

    def __init__(self, n_max: int, *, rank_tol: float = 1e-12) -> None:
        if n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {n_max}")
        self.n_max = n_max
        self.basis = full_dirac_basis(n_max)
        self.dim_H = len(self.basis)
        self._rank_tol = rank_tol

        # Weyl-sector basis used to build the half-block multipliers.
        weyl_basis = spinor_basis(n_max)
        dim_weyl = len(weyl_basis)
        assert 2 * dim_weyl == self.dim_H

        all_labels = allowed_multiplier_labels(n_max)
        labels: List[Tuple[int, int, int]] = []
        matrices: List[np.ndarray] = []
        for (N, L, M) in all_labels:
            weyl_mat = build_spinor_multiplier_matrix(N, L, M, weyl_basis)
            if np.allclose(weyl_mat, 0, atol=1e-15):
                continue
            full_mat = np.zeros(
                (self.dim_H, self.dim_H), dtype=np.complex128
            )
            full_mat[:dim_weyl, :dim_weyl] = weyl_mat
            full_mat[dim_weyl:, dim_weyl:] = weyl_mat
            labels.append((N, L, M))
            matrices.append(full_mat)

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
        I = np.eye(self.dim_H, dtype=np.complex128)
        return self.contains(I, tol=tol)

    def is_star_closed(
        self, *, tol: float = 1e-10,
    ) -> Tuple[bool, List[Tuple[int, float]]]:
        failures = []
        for i, M in enumerate(self.multiplier_matrices):
            ok, residual = self.contains(M.conj().T, tol=tol)
            if not ok:
                failures.append((i, residual))
        return (len(failures) == 0, failures)


# ---------------------------------------------------------------------------
# Convenience: graph distance on the full-Dirac labels
# ---------------------------------------------------------------------------


def full_dirac_graph_distance(
    v: FullDiracLabel, w: FullDiracLabel,
) -> int:
    """Combinatorial graph distance with chirality treated as one extra hop.

    d(v, w) = |Delta n_fock| + |Delta l| + |Delta m_j| + |Delta chirality|/2

    where |Delta chirality| = 0 (same sector) or 2 (opposite sector). We
    halve it so chirality flip = "1 hop" in the graph distance.

    This matches Paper 32 §3.3's reading of chirality as a Z_2 ladder.
    """
    same = 0 if v.chirality == w.chirality else 1
    return (
        abs(v.n_fock - w.n_fock)
        + abs(v.l - w.l)
        + abs(v.two_m_j - w.two_m_j) // 2
        + same
    )


def full_dirac_graph_distance_matrix(
    op_sys: "FullDiracTruncatedOperatorSystem",
) -> np.ndarray:
    """N x N integer graph-distance matrix on the full-Dirac basis."""
    N = op_sys.dim_H
    dist = np.zeros((N, N), dtype=int)
    for i in range(N):
        for j in range(N):
            dist[i, j] = full_dirac_graph_distance(
                op_sys.basis[i], op_sys.basis[j]
            )
    return dist


def full_dirac_label_strings(
    op_sys: "FullDiracTruncatedOperatorSystem",
) -> List[str]:
    """Pretty label strings |n,l,m_j,chi>."""
    out = []
    for b in op_sys.basis:
        if b.two_m_j % 2 == 0:
            mj_str = f"{b.two_m_j // 2}"
        else:
            sign = "+" if b.two_m_j > 0 else "-"
            mj_str = f"{sign}{abs(b.two_m_j)}/2"
        chi_str = "+" if b.chirality > 0 else "-"
        out.append(f"|{b.n_fock},{b.l},{mj_str},{chi_str}>")
    return out


# ---------------------------------------------------------------------------
# Connes distance specialization for the full-Dirac sector
# ---------------------------------------------------------------------------


def compute_distance_matrix_full_dirac(
    op_sys: "FullDiracTruncatedOperatorSystem",
    D: Optional[np.ndarray] = None,
    *,
    solver: str = "SCS",
    solver_kwargs: Optional[dict] = None,
    verbose: bool = False,
    progress: bool = False,
) -> np.ndarray:
    """Compute the full N x N Connes distance matrix on a full-Dirac op_sys.

    Drops into geovac.connes_distance.compute_distance_matrix unchanged
    (the FullDiracTruncatedOperatorSystem is API-compatible with
    TruncatedOperatorSystem). Wrapped here for documentation.

    Parameters
    ----------
    op_sys : FullDiracTruncatedOperatorSystem
    D : np.ndarray, optional
        Dirac matrix of shape (op_sys.dim_H, op_sys.dim_H). Default:
        camporesi_higuchi_full_dirac_matrix(op_sys.basis).
    solver, solver_kwargs, verbose, progress : as in
        compute_distance_matrix.

    Returns
    -------
    Real symmetric ndarray of shape (N, N) with zero diagonal.
    """
    # Local import to avoid importing cvxpy at module load time.
    from geovac.connes_distance import compute_distance_matrix

    if D is None:
        D = camporesi_higuchi_full_dirac_matrix(op_sys.basis)
    return compute_distance_matrix(
        op_sys,
        D=D,
        solver=solver,
        solver_kwargs=solver_kwargs,
        verbose=verbose,
        progress=progress,
    )
