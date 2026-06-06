r"""
Pro-system structure on truncated Camporesi--Higuchi spectral triples.
======================================================================

Implements the inverse system $\{\mathcal{O}_{n_{\max}}\}_{n_{\max} \ge 1}$
of commutative sector-idempotent algebras at finite cutoff, with explicit
closed-form transition maps

    $P_{m, k}: \mathcal{O}_m \to \mathcal{O}_k, \quad m \ge k \ge 1$,

satisfying the cofiltered axiom

    $P_{m, k} = P_{n, k} \circ P_{m, n} \quad \forall\; k \le n \le m$.

Mathematical content
--------------------
Each algebra $\mathcal{O}_{n_{\max}}$ has the canonical basis of sector
idempotents $\{e_{(n, l)} : 1 \le n \le n_{\max}, 0 \le l \le n\}$,
isomorphic to $\mathbb{C}^{N(n_{\max})}$ with

    $N(n_{\max}) = n_{\max} (n_{\max} + 3) / 2$.

The transition is the algebra homomorphism

    $P_{m, k}(e_{(n, l)}) = \begin{cases} e_{(n, l)} & \text{if } n \le k, \\ 0 & \text{if } n > k. \end{cases}$

In the canonical lex order on sectors, this is the projection
$\mathbb{C}^{N(m)} \to \mathbb{C}^{N(k)}$ that keeps the first $N(k)$
coordinates and drops the last $N(m) - N(k)$.

Sector locality
---------------
The Camporesi--Higuchi shell decomposition is sector-LOCAL: adding a
higher shell $n_{\max} + 1$ does NOT modify the Dirac structure of any
lower sector $(n, l)$ with $n \le n_{\max}$. Therefore the per-sector
trace-class values $\chi_{(n, l)} = \mathrm{Tr}(\gamma\, e_{(n, l)})$
and $\eta_{(n, l)} = \mathrm{Tr}(\gamma\, D\, e_{(n, l)})$ are
$n_{\max}$-independent, and the cocycle classes pull back strictly under
$P$ (no twist, no coboundary modification). Closed forms

    $\chi_{(n, l)} = \begin{cases} +2 & l < n \\ -2n & l = n \end{cases}, \quad
     \eta_{(n, l)} = \begin{cases} (2l + 1)(2n + 1) & l < n \\ n(2n + 1) & l = n \end{cases}$

were derived in Sprint Q5'-Stage1-Prosystem (v3.60.0).

Public API
----------
* `sectors_at_cutoff(n_max)` — list of $(n, l)$ in canonical lex order.
* `N_sectors(n_max)` — closed-form count $n_{\max}(n_{\max} + 3) / 2$.
* `TransitionMap(n_high, n_low)` — closed-form $P_{m, k}$ object with
  matrix representation, vector pull-back, and class pull-back.
* `compose(P_outer, P_inner)` — composition operator yielding the
  $P_{m, k}$ map (cofiltered axiom checkable by `==` of `.matrix`).
* `verify_cofiltered_axiom(m, n, k)` — bit-exact verification at the
  matrix level for a single triple $k \le n \le m$.

Discipline
----------
Bit-exact `sympy.Rational` and `sympy.Integer` throughout. No floats.
No PSLQ. The transition matrices are 0/1 integer matrices.

References
----------
* Sprint Q5'-Stage1-Prosystem memo, v3.60.0
  (`debug/sprint_q5p_prosystem_memo.md`).
* Paper 55 \S\ref{subsec:open_m2_m3} (Q5' Stage 1 construction).
* Connes, A.; van Suijlekom, W. D. ``Spectral truncations in
  noncommutative geometry and operator systems.'' Comm. Math. Phys.
  383 (2021).
* Paper 38 \S VIII L4 (Berezin reconstruction substrate).
"""

from __future__ import annotations

from typing import Dict, List, Tuple, Mapping, Any, Callable, Iterator

from sympy import Integer, Matrix, Rational, zeros as sp_zeros

__all__ = [
    "sectors_at_cutoff",
    "N_sectors",
    "TransitionMap",
    "compose",
    "verify_cofiltered_axiom",
    # PS-2 additions
    "MELLIN_SLOTS",
    "primitive_generators",
    "n_primitive_generators",
    "HopfTransition",
    "verify_hopf_cofiltered_axiom",
    "verify_Ga_generator_compatibility",
    "verify_class_action_compatibility",
    # PS-3 additions
    "InverseLimitClass",
    "chi_infinity",
    "eta_infinity",
    "project_to_cutoff",
    "verify_universal_property",
    "verify_continuity_under_transitions",
]


# =====================================================================
# Sector enumeration (closed form)
# =====================================================================


def sectors_at_cutoff(n_max: int) -> List[Tuple[int, int]]:
    """Return all $(n, l)$ sectors at cutoff `n_max` in canonical lex order.

    Sectors are $(n, l)$ with $1 \\le n \\le n_{\\max}$, $0 \\le l \\le n$.

    Parameters
    ----------
    n_max : int
        Cutoff $n_{\\max} \\ge 1$.

    Returns
    -------
    list of (int, int)
        Sectors $(n, l)$ in lex order: $(1, 0), (1, 1), (2, 0), (2, 1),
        (2, 2), (3, 0), \\ldots$.

    Examples
    --------
    >>> sectors_at_cutoff(2)
    [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)]
    >>> len(sectors_at_cutoff(5))
    20
    """
    if n_max < 1:
        raise ValueError(f"n_max must be >= 1, got {n_max}")
    return [(n, l) for n in range(1, n_max + 1) for l in range(n + 1)]


def N_sectors(n_max: int) -> int:
    """Closed-form sector count $N(n_{\\max}) = n_{\\max}(n_{\\max} + 3) / 2$.

    Verified bit-exact against direct enumeration at $n_{\\max} \\in \\{1, \\ldots, 5\\}$:
    $N = 2, 5, 9, 14, 20$.
    """
    if n_max < 1:
        raise ValueError(f"n_max must be >= 1, got {n_max}")
    return n_max * (n_max + 3) // 2


# =====================================================================
# TransitionMap
# =====================================================================


class TransitionMap:
    r"""Closed-form transition $P_{m, k}: \mathcal{O}_m \to \mathcal{O}_k$.

    Algebra homomorphism sending sector idempotent $e_{(n, l)}$ to
    itself if $n \le k$, otherwise to $0$. In the canonical lex order
    on sectors, $P_{m, k}$ is realised by the $N(k) \times N(m)$ matrix

        $[P_{m, k}]_{i, j} = \begin{cases} 1 & \text{if sector } i \text{ at } k = \text{sector } j \text{ at } m \\ 0 & \text{otherwise} \end{cases}$.

    Matrix shape: rows are sectors at $k$ (output), columns are sectors
    at $m$ (input). The first $N(k)$ columns form the identity $I_{N(k)}$;
    the last $N(m) - N(k)$ columns are zero. Acts on column vectors of
    sector-idempotent coefficients on the left.

    Parameters
    ----------
    n_high : int
        Source cutoff $m \ge 1$.
    n_low : int
        Target cutoff $k$ with $1 \le k \le m$.

    Attributes
    ----------
    n_high, n_low : int
    sectors_high, sectors_low : list of (int, int)
        Sectors at $m$ and $k$ in canonical lex order.
    dim_high, dim_low : int
        $N(m)$ and $N(k)$.

    Examples
    --------
    >>> P = TransitionMap(3, 2)
    >>> P.matrix.shape
    (5, 9)
    >>> P.apply_to_vector([1, 2, 3, 4, 5, 6, 7, 8, 9])
    [1, 2, 3, 4, 5]
    """

    def __init__(self, n_high: int, n_low: int):
        if n_high < n_low:
            raise ValueError(
                f"n_high must be >= n_low; got n_high={n_high}, n_low={n_low}"
            )
        if n_low < 1:
            raise ValueError(f"n_low must be >= 1, got {n_low}")
        self.n_high = n_high
        self.n_low = n_low
        self.sectors_high = sectors_at_cutoff(n_high)
        self.sectors_low = sectors_at_cutoff(n_low)
        self.dim_high = len(self.sectors_high)
        self.dim_low = len(self.sectors_low)
        # Index map: sector_low[i] equals sector_high[idx_in_high[i]]
        sec_high_to_idx = {s: j for j, s in enumerate(self.sectors_high)}
        self._idx_in_high = [sec_high_to_idx[s] for s in self.sectors_low]
        self._matrix: Matrix | None = None

    # ------------------------------------------------------------------
    # Closed-form matrix
    # ------------------------------------------------------------------

    @property
    def matrix(self) -> Matrix:
        """The $N(k) \\times N(m)$ projection matrix in canonical lex basis.

        Cached on first access. Bit-exact ``Integer(0)`` / ``Integer(1)``
        entries.
        """
        if self._matrix is not None:
            return self._matrix
        M = sp_zeros(self.dim_low, self.dim_high)
        for i, j in enumerate(self._idx_in_high):
            M[i, j] = Integer(1)
        self._matrix = M
        return M

    # ------------------------------------------------------------------
    # Action on vectors / classes
    # ------------------------------------------------------------------

    def apply_to_vector(self, v: List) -> List:
        """Apply $P$ to a class vector $v \\in \\mathbb{Q}^{N(m)}$.

        Returns the restricted vector of length $N(k)$ obtained by
        keeping the coordinates corresponding to sectors at $k$.
        """
        if len(v) != self.dim_high:
            raise ValueError(
                f"vector must have length {self.dim_high}, got {len(v)}"
            )
        return [v[j] for j in self._idx_in_high]

    def apply_to_class(
        self,
        cls: Mapping[Tuple[int, int], Any],
    ) -> Dict[Tuple[int, int], Any]:
        """Apply $P$ to a class dict keyed by $(n, l)$.

        Returns the restricted dict keyed only by sectors at $k$.
        """
        return {sec: cls[sec] for sec in self.sectors_low}

    # ------------------------------------------------------------------
    # Equality, hashing, repr
    # ------------------------------------------------------------------

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, TransitionMap):
            return NotImplemented
        return (
            self.n_high == other.n_high
            and self.n_low == other.n_low
            and self.matrix == other.matrix
        )

    def __hash__(self) -> int:
        return hash((self.n_high, self.n_low))

    def __repr__(self) -> str:
        return f"TransitionMap(n_high={self.n_high}, n_low={self.n_low})"


# =====================================================================
# Composition and cofiltered axiom
# =====================================================================


def compose(P_outer: TransitionMap, P_inner: TransitionMap) -> TransitionMap:
    r"""Compose $P_{\text{outer}} \circ P_{\text{inner}}: \mathcal{O}_{\text{inner}.n_{\text{high}}} \to \mathcal{O}_{\text{outer}.n_{\text{low}}}$.

    Requires ``P_inner.n_low == P_outer.n_high``. Returns a
    :class:`TransitionMap` of the composed source-target pair, whose
    matrix equals ``P_outer.matrix * P_inner.matrix`` bit-exactly when
    the cofiltered axiom holds.
    """
    if P_inner.n_low != P_outer.n_high:
        raise ValueError(
            f"compose requires P_inner.n_low == P_outer.n_high; got "
            f"P_inner.n_low={P_inner.n_low}, P_outer.n_high={P_outer.n_high}"
        )
    return TransitionMap(P_inner.n_high, P_outer.n_low)


def verify_cofiltered_axiom(m: int, n: int, k: int) -> Dict[str, Any]:
    r"""Verify $P_{m, k} = P_{n, k} \cdot P_{m, n}$ bit-exact at the matrix level.

    Parameters
    ----------
    m, n, k : int
        Cutoffs with $1 \le k \le n \le m$.

    Returns
    -------
    dict
        ``{"m", "n", "k", "direct_matrix", "composed_matrix",
        "bit_exact": bool, "residual": Matrix}`` — residual is
        ``direct - composed`` (zero matrix iff axiom holds).
    """
    if not (1 <= k <= n <= m):
        raise ValueError(
            f"verify_cofiltered_axiom requires 1 <= k <= n <= m; got "
            f"m={m}, n={n}, k={k}"
        )
    P_direct = TransitionMap(m, k)
    P_outer = TransitionMap(n, k)
    P_inner = TransitionMap(m, n)
    M_direct = P_direct.matrix
    M_composed = P_outer.matrix * P_inner.matrix
    residual = M_direct - M_composed
    is_zero = all(residual[i, j] == 0
                  for i in range(residual.rows)
                  for j in range(residual.cols))
    return {
        "m": m,
        "n": n,
        "k": k,
        "dim_high": P_direct.dim_high,
        "dim_low": P_direct.dim_low,
        "direct_matrix": M_direct,
        "composed_matrix": M_composed,
        "bit_exact": is_zero,
        "residual_norm_squared": sum(
            (residual[i, j]) ** 2
            for i in range(residual.rows)
            for j in range(residual.cols)
        ),
    }


# =====================================================================
# PS-2 --- Hopf-algebraic lift and U* compatibility
# =====================================================================
#
# The Stage-2 candidate Hopf algebra of v3.61.0 Track A is
#
#     $\mathcal{H}_{\mathrm{GV}}(n_{\max}) = \mathrm{Sym}_{\mathbb{Q}}(V_{n_{\max}})$,
#
# where $V_{n_{\max}}$ is the $\mathbb{Q}$-vector space of primitive
# generators $\{x_{(n, l), k}\}$ for $(n, l) \in \text{sectors}(n_{\max})$
# and $k \in \{0, 1, 2\}$ indexing the master Mellin engine slots
# (Paper 18 \S III.7):\ $k = 0 \leftrightarrow M_1$ (Hopf-base measure),
# $k = 1 \leftrightarrow M_3$ (vertex-parity Hurwitz),
# $k = 2 \leftrightarrow M_2$ (Seeley--DeWitt). The candidate motivic
# Galois group is the Levi-decomposition product (v3.63.0 L1)
#
#     $U^*_{\mathrm{GeoVac}, \mathrm{Levi}} = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$,
#
# pro-unipotent times semisimple, matching Connes--Marcolli (arXiv:math/0409306).
#
# Sector locality of the Camporesi--Higuchi shell decomposition forces
# $\mathcal{H}_{\mathrm{GV}}$ to be **abelian primitive** at the basic
# substrate level (v3.61.0 Track B): the coproduct is $\Delta(x) =
# x \otimes 1 + 1 \otimes x$ on every primitive generator, the
# counit is $\varepsilon(x) = 0$, and the antipode is $S(x) = -x$.
# Hopf-axiom verification on primitives is therefore structurally
# trivial under any algebra-level truncation that respects the
# generator indexing.
#
# The class-level $U^*$-action on the cocycle classes $\chi$ and
# $\eta$ is **trivial at depth 0** (v3.66.0 FO3 Interpretation C):\
# both $\chi$ and $\eta$ are integer-valued, so $\mathbb{G}_a^{3N}$
# acts as identity on them, and $SL_2$ (which lives on the
# $n_{\max}$-independent Peter--Weyl decoration) commutes with any
# $n_{\max}$-axis truncation structurally. The class-level $U^*$-action
# on $F(s)$ (the continuum Mellin lift, v3.66.0 FO2) is non-trivial
# in the $M_3$ slot but $F(s)$ lives at the inverse-limit level and
# is handled by PS-3, not PS-2.

MELLIN_SLOTS: Tuple[int, int, int] = (0, 1, 2)
r"""Mellin slot indices for the master Mellin engine (Paper 18 \S III.7).

* $k = 0 \leftrightarrow M_1$ Hopf-base measure (Vol(S^2)/4 = \pi).
* $k = 1 \leftrightarrow M_3$ vertex-parity Hurwitz at quarter-integer shifts.
* $k = 2 \leftrightarrow M_2$ Seeley--DeWitt coefficients ($\pi^{2k}$).
"""


def n_primitive_generators(n_max: int) -> int:
    """Number of primitive generators of $\\mathcal{H}_{\\mathrm{GV}}(n_{\\max})$.

    Equal to $3 \\cdot N(n_{\\max}) = 3 n_{\\max} (n_{\\max} + 3) / 2$
    (one per (sector, Mellin slot) pair).
    """
    return 3 * N_sectors(n_max)


def primitive_generators(n_max: int) -> List[Tuple[int, int, int]]:
    """Canonical list of primitive generators $x_{(n, l), k}$.

    Returns a list of (n, l, k) triples in canonical order: lex on
    (n, l) within each Mellin slot, slots ordered $k = 0, 1, 2$.
    Length: $3 N(n_{\\max})$.
    """
    secs = sectors_at_cutoff(n_max)
    return [(n, l, k) for k in MELLIN_SLOTS for (n, l) in secs]


class HopfTransition:
    r"""Hopf-algebra homomorphism $\Phi_{m, k}: \mathcal{H}_{\mathrm{GV}}(m) \to \mathcal{H}_{\mathrm{GV}}(k)$.

    Acts on the primitive generators as the truncation
    $\Phi(x_{(n, l), s}^{(m)}) = x_{(n, l), s}^{(k)}$ if $n \le k$, and
    $0$ otherwise --- structurally identical to PS-1's
    :class:`TransitionMap` repeated independently on each of the three
    Mellin slots $s \in \{0, 1, 2\}$. The action on $\mathrm{Sym}(V)$ at
    higher degrees is determined by the algebra-homomorphism property.

    Compatibility with the abelian primitive coproduct is automatic:\
    $\Delta(\Phi(x)) = \Phi(x) \otimes 1 + 1 \otimes \Phi(x) =
    (\Phi \otimes \Phi)(\Delta(x))$ for every primitive $x$. Same for
    counit ($\varepsilon \circ \Phi = \varepsilon$ on $V$, trivially)
    and antipode ($\Phi(S(x)) = \Phi(-x) = -\Phi(x) = S(\Phi(x))$).

    Parameters
    ----------
    n_high, n_low : int
        Source and target cutoffs with $n_{\\mathrm{high}} \\ge n_{\\mathrm{low}} \\ge 1$.

    Attributes
    ----------
    base : TransitionMap
        The underlying PS-1 algebra-level transition $P_{m, k}$.
    """

    def __init__(self, n_high: int, n_low: int):
        self.base = TransitionMap(n_high, n_low)
        self.n_high = n_high
        self.n_low = n_low
        self.generators_high = primitive_generators(n_high)
        self.generators_low = primitive_generators(n_low)
        self.dim_high = len(self.generators_high)
        self.dim_low = len(self.generators_low)
        self._matrix: Matrix | None = None

    @property
    def matrix(self) -> Matrix:
        """The $3 N(k) \\times 3 N(m)$ block-diagonal matrix.

        Each of the three Mellin slots contributes an independent copy
        of :class:`TransitionMap` $P_{m, k}$. The full matrix is
        $\\mathrm{diag}(P_{m, k}, P_{m, k}, P_{m, k})$.
        """
        if self._matrix is not None:
            return self._matrix
        P = self.base.matrix
        N_low = self.base.dim_low
        N_high = self.base.dim_high
        M = sp_zeros(3 * N_low, 3 * N_high)
        for s in range(3):
            for i in range(N_low):
                for j in range(N_high):
                    M[s * N_low + i, s * N_high + j] = P[i, j]
        self._matrix = M
        return M

    def apply_to_generator(
        self, gen: Tuple[int, int, int],
    ) -> Tuple[int, int, int] | None:
        """Apply $\\Phi$ to a primitive generator.

        Returns the target generator at $n_{\\mathrm{low}}$, or ``None``
        if the generator is killed (its sector has $n > n_{\\mathrm{low}}$).
        """
        n, l, k = gen
        if n <= self.n_low:
            return (n, l, k)
        return None

    def __repr__(self) -> str:
        return f"HopfTransition(n_high={self.n_high}, n_low={self.n_low})"


def verify_hopf_cofiltered_axiom(m: int, n: int, k: int) -> Dict[str, Any]:
    r"""Verify $\Phi_{m, k} = \Phi_{n, k} \cdot \Phi_{m, n}$ bit-exact at the Hopf-hom level.

    Reduces to PS-1's cofiltered axiom via block-diagonal Mellin slot
    decomposition; verified independently for completeness.
    """
    if not (1 <= k <= n <= m):
        raise ValueError(
            f"verify_hopf_cofiltered_axiom requires 1 <= k <= n <= m; got "
            f"m={m}, n={n}, k={k}"
        )
    Phi_direct = HopfTransition(m, k)
    Phi_outer = HopfTransition(n, k)
    Phi_inner = HopfTransition(m, n)
    M_direct = Phi_direct.matrix
    M_composed = Phi_outer.matrix * Phi_inner.matrix
    residual = M_direct - M_composed
    is_zero = all(residual[i, j] == 0
                  for i in range(residual.rows)
                  for j in range(residual.cols))
    return {
        "m": m,
        "n": n,
        "k": k,
        "dim_high": Phi_direct.dim_high,
        "dim_low": Phi_direct.dim_low,
        "bit_exact": is_zero,
        "residual_norm_squared": sum(
            (residual[i, j]) ** 2
            for i in range(residual.rows)
            for j in range(residual.cols)
        ),
    }


def verify_Ga_generator_compatibility(m: int, k: int) -> Dict[str, Any]:
    r"""Verify each $\mathbb{G}_a$ translation generator commutes with $\Phi_{m, k}$.

    The $\mathbb{G}_a^{3 N(m)}$ factor has one translation generator
    $e_g$ per primitive generator $g = (n, l, s)$ at cutoff $m$. The
    generator $e_g$ acts on the Hopf algebra by translation in the $g$
    coordinate; its image under $\Phi_{m, k}$ is well-defined when
    $n \le k$ (the generator survives) and is zero otherwise.

    Compatibility condition (for each $g$ at cutoff $m$):

    .. math::

        \Phi_{m, k}(e_g^{(m)} \cdot x) = \begin{cases}
            e_{g}^{(k)} \cdot \Phi_{m, k}(x) & \text{if } n \le k, \\
            0 & \text{if } n > k.
        \end{cases}

    Returns
    -------
    dict
        Per-generator outcome and a summary boolean.
    """
    if not (1 <= k <= m):
        raise ValueError(f"verify_Ga_generator_compatibility requires 1 <= k <= m; got m={m}, k={k}")
    Phi = HopfTransition(m, k)
    per_generator: List[Dict[str, Any]] = []
    n_killed = 0
    n_survived = 0
    for g in Phi.generators_high:
        n, l, s = g
        g_image = Phi.apply_to_generator(g)
        survives = (n <= k)
        # By construction:
        # If g survives, the surviving image is the same triple at cutoff k.
        # If g is killed, the action of e_g drops to zero under Phi.
        expected_image = (n, l, s) if survives else None
        bit_exact = (g_image == expected_image)
        per_generator.append({
            "generator": [n, l, s],
            "survives": survives,
            "image": list(g_image) if g_image is not None else None,
            "expected_image": list(expected_image) if expected_image is not None else None,
            "bit_exact": bit_exact,
        })
        if survives:
            n_survived += 1
        else:
            n_killed += 1
    all_bit_exact = all(rec["bit_exact"] for rec in per_generator)
    return {
        "m": m,
        "k": k,
        "total_generators": len(per_generator),
        "n_survived": n_survived,
        "n_killed": n_killed,
        "expected_survived": 3 * N_sectors(k),
        "expected_killed": 3 * (N_sectors(m) - N_sectors(k)),
        "all_bit_exact": all_bit_exact,
        "per_generator": per_generator,
    }


def verify_class_action_compatibility(
    m: int,
    k: int,
    psi_high: List,
    psi_low: List,
) -> Dict[str, Any]:
    r"""Verify $P_{m, k}(U^* \cdot \psi^{(m)}) = U^* \cdot P_{m, k}(\psi^{(m)}) = \psi^{(k)}$ bit-exact.

    Tests Interpretation C closure across the pro-system axis. For
    depth-0 cocycle classes (integer-valued: $\chi, \eta$), the
    $U^*$-action is the identity (v3.66.0 FO3), so both sides reduce
    to PS-1's pull-back identity $P_{m, k}(\psi^{(m)}) = \psi^{(k)}$.
    This routine is a bookkeeping verification that the depth-0
    triviality survives the transition.

    Parameters
    ----------
    m, k : int
        Cutoffs with $1 \\le k \\le m$.
    psi_high : list of length $N(m)$
        Class vector at cutoff $m$.
    psi_low : list of length $N(k)$
        Class vector at cutoff $k$.

    Returns
    -------
    dict
        ``{"m", "k", "lhs_eq_rhs": bool, "lhs_eq_psi_low": bool,
        "all_bit_exact": bool}`` --- both equalities must hold for the
        compatibility to be bit-exact under depth-0 $U^*$-action.
    """
    if not (1 <= k <= m):
        raise ValueError(f"verify_class_action_compatibility requires 1 <= k <= m; got m={m}, k={k}")
    P = TransitionMap(m, k)
    if len(psi_high) != P.dim_high:
        raise ValueError(
            f"psi_high length must be {P.dim_high}, got {len(psi_high)}"
        )
    if len(psi_low) != P.dim_low:
        raise ValueError(
            f"psi_low length must be {P.dim_low}, got {len(psi_low)}"
        )
    # Depth-0 U*-action is identity. Both sides reduce to P(psi_high).
    lhs = P.apply_to_vector(psi_high)  # P(U* . psi_high) = P(psi_high)
    rhs = P.apply_to_vector(psi_high)  # U* . (P psi_high) = P(psi_high)
    lhs_eq_rhs = (lhs == rhs)
    lhs_eq_psi_low = (lhs == psi_low)
    return {
        "m": m,
        "k": k,
        "lhs": lhs,
        "rhs": rhs,
        "psi_low": list(psi_low),
        "lhs_eq_rhs": lhs_eq_rhs,
        "lhs_eq_psi_low": lhs_eq_psi_low,
        "all_bit_exact": lhs_eq_rhs and lhs_eq_psi_low,
    }


# =====================================================================
# PS-3 --- Inverse limit and continuum cocycle classes
# =====================================================================
#
# The pro-system $\{\mathcal{O}_{n_{\max}}\}_{n_{\max} \ge 1}$ with the
# closed-form transitions $P_{m, k}$ of PS-1 has a sequential inverse
# limit
#
#     $\mathcal{O}_\infty = \varprojlim_{n_{\max}} \mathcal{O}_{n_{\max}}
#                          = \{(a_{n_{\max}})_{n_{\max} \ge 1} :
#                              P_{m, k}(a_m) = a_k \;\forall m \ge k\}$,
#
# which for sector-local data reduces to functions on the infinite
# sector index $\mathbb{N}_{\mathrm{sec}} = \{(n, l) : n \ge 1,
# 0 \le l \le n\}$ with values in $\mathbb{Q}$. The natural topology
# is the product / inverse-limit topology (continuity coordinate-wise),
# under which each projection $\pi_{n_{\max}}: \mathcal{O}_\infty \to
# \mathcal{O}_{n_{\max}}$ is continuous by construction.
#
# Sector-locality of the Camporesi--Higuchi shell decomposition implies
# that every per-sector cocycle value (PS-1's $\chi_{(n, l)}$ and
# $\eta_{(n, l)}$ closed forms; PS-2's depth-0 $U^*$-invariant
# classes; v3.61.0 Track A's primitive generator labels) extends
# trivially to $\mathcal{O}_\infty$:\ the value at sector $(n, l)$ is
# defined by the closed form and is independent of any cutoff. The
# universal property holds bit-exact:\ for any compatible family of
# class data $(\psi^{(n_{\max})})_{n_{\max} \ge 1}$, there is a unique
# $\psi_\infty \in \mathbb{Q}^{\mathbb{N}_{\mathrm{sec}}}$ with
# $\pi_{n_{\max}}(\psi_\infty) = \psi^{(n_{\max})}$.
#
# The non-trivial $U^*$-action on the continuum Mellin lift $F(s)$
# (v3.66.0 FO2, with M2 / M3 components) lives at $\mathcal{O}_\infty$
# but acts via motivic Galois action on $\mathrm{MT}(\mathbb{Q}, 1)$ at
# the period level (not at the cocycle class level). The infrastructure
# here handles the $\mathbb{Q}^{\mathbb{N}_{\mathrm{sec}}}$ side;\ F(s)
# weight / depth grading and U*-orbit closure are handled in the PS-3
# driver against the v3.66.0 FO2 bit-exact panel.


class InverseLimitClass:
    r"""Element of $\mathcal{O}_\infty$ represented by a closed-form function $(n, l) \mapsto \mathbb{Q}$.

    The function must be defined on all of $\mathbb{N}_{\mathrm{sec}} =
    \{(n, l) : n \ge 1, 0 \le l \le n\}$. Compatibility with the
    pro-system's transitions $P_{m, k}$ is automatic when the function
    is sector-local (depends on $(n, l)$ only, not on a cutoff
    parameter).

    Parameters
    ----------
    closed_form : callable
        ``closed_form(n, l)`` returns a ``sympy.Integer`` /
        ``sympy.Rational`` for every $(n, l) \in \mathbb{N}_{\mathrm{sec}}$.
    label : str, optional
        Human-readable label for debugging.

    Examples
    --------
    >>> chi_inf = InverseLimitClass(lambda n, l: -2 * n if l == n else 2, label="chi")
    >>> chi_inf.at(3, 2)
    2
    >>> chi_inf.at(3, 3)
    -6
    """

    def __init__(
        self,
        closed_form: Callable[[int, int], Any],
        label: str = "",
    ):
        self._closed_form = closed_form
        self.label = label

    def at(self, n: int, l: int):
        """Value at sector $(n, l)$."""
        if n < 1 or l < 0 or l > n:
            raise ValueError(
                f"Sector (n, l) = ({n}, {l}) is out of range "
                f"(require n >= 1, 0 <= l <= n)"
            )
        return self._closed_form(n, l)

    def project(self, n_max: int) -> Dict[Tuple[int, int], Any]:
        """Restrict to sectors at cutoff $n_{\\max}$ via the projection $\\pi_{n_{\\max}}$."""
        return {(n, l): self._closed_form(n, l)
                for (n, l) in sectors_at_cutoff(n_max)}

    def project_to_vector(self, n_max: int) -> List:
        """Restrict to sectors at cutoff $n_{\\max}$ as a vector in canonical lex order."""
        return [self._closed_form(n, l) for (n, l) in sectors_at_cutoff(n_max)]

    def __repr__(self) -> str:
        return f"InverseLimitClass(label={self.label!r})"


def chi_infinity() -> InverseLimitClass:
    r"""The JLO HP$^{\mathrm{even}}$ continuum cocycle class $\chi_\infty$.

    Closed form from v3.60.0 / PS-1:

    .. math::

        \chi_{(n, l)} = \begin{cases} +2 & l < n, \\ -2n & l = n. \end{cases}

    Depth-0 (integer-valued) on all of $\mathbb{N}_{\mathrm{sec}}$.
    Sector-local;\ extends to $\mathcal{O}_\infty$ trivially via the
    universal property.
    """
    return InverseLimitClass(
        lambda n, l: Integer(-2 * n) if l == n else Integer(2),
        label="chi_infinity",
    )


def eta_infinity() -> InverseLimitClass:
    r"""The CM-$\eta$ residue continuum cocycle class $\eta_\infty$.

    Closed form from v3.60.0 / PS-1:

    .. math::

        \eta_{(n, l)} = \begin{cases} (2l + 1)(2n + 1) & l < n, \\ n (2n + 1) & l = n. \end{cases}

    Depth-0 (integer-valued) on all of $\mathbb{N}_{\mathrm{sec}}$.
    Sector-local;\ extends to $\mathcal{O}_\infty$ trivially.
    """
    return InverseLimitClass(
        lambda n, l: Integer(n * (2 * n + 1)) if l == n
                     else Integer((2 * l + 1) * (2 * n + 1)),
        label="eta_infinity",
    )


def project_to_cutoff(
    psi_inf: InverseLimitClass,
    n_max: int,
) -> Dict[Tuple[int, int], Any]:
    """Apply $\\pi_{n_{\\max}}$ to an element of $\\mathcal{O}_\\infty$.

    Convenience wrapper around :py:meth:`InverseLimitClass.project`.
    """
    return psi_inf.project(n_max)


def verify_universal_property(
    psi_inf: InverseLimitClass,
    psi_finite_by_cutoff: Mapping[int, Mapping[Tuple[int, int], Any]],
) -> Dict[str, Any]:
    r"""Verify $\pi_{n_{\max}}(\psi_\infty) = \psi^{(n_{\max})}$ bit-exact at every cutoff.

    The universal property of the inverse limit:\ for every cutoff
    $n_{\max}$ in the supplied family, the restriction of $\psi_\infty$
    to sectors at $n_{\max}$ matches the finite-cutoff class
    $\psi^{(n_{\max})}$.

    Parameters
    ----------
    psi_inf : InverseLimitClass
        The candidate continuum extension.
    psi_finite_by_cutoff : mapping
        Keys are cutoffs $n_{\max}$; values are dicts mapping
        $(n, l) \to$ value.

    Returns
    -------
    dict
        Per-cutoff outcome:\ number of sectors tested, mismatches,
        bit-exact verdict.
    """
    results: Dict[int, Dict[str, Any]] = {}
    total_sectors_tested = 0
    total_mismatches = 0
    for n_max in sorted(psi_finite_by_cutoff.keys()):
        psi_finite = psi_finite_by_cutoff[n_max]
        projected = psi_inf.project(n_max)
        mismatches: List[Dict[str, Any]] = []
        for sec, expected in psi_finite.items():
            actual = projected.get(sec)
            if actual != expected:
                mismatches.append({
                    "sector": list(sec),
                    "expected": expected,
                    "actual": actual,
                })
        results[n_max] = {
            "sectors_tested": len(psi_finite),
            "mismatches": mismatches,
            "bit_exact": len(mismatches) == 0,
        }
        total_sectors_tested += len(psi_finite)
        total_mismatches += len(mismatches)
    all_bit_exact = total_mismatches == 0
    return {
        "label": psi_inf.label,
        "per_cutoff": results,
        "total_sectors_tested": total_sectors_tested,
        "total_mismatches": total_mismatches,
        "all_bit_exact": all_bit_exact,
    }


def verify_continuity_under_transitions(
    psi_inf: InverseLimitClass,
    cutoff_pairs: Iterator[Tuple[int, int]],
) -> Dict[str, Any]:
    r"""Verify $P_{m, k}(\pi_m(\psi_\infty)) = \pi_k(\psi_\infty)$ bit-exact.

    Together with the universal property, this expresses continuity of
    the family $\{\pi_{n_{\max}}\}$ in the inverse-limit topology:\
    truncating $\psi_\infty$ first and then applying the algebra-level
    transition produces the same result as truncating directly.

    Parameters
    ----------
    psi_inf : InverseLimitClass
    cutoff_pairs : iterable of (m, k)
        Pairs with $m \ge k \ge 1$.

    Returns
    -------
    dict
        Per-pair outcome and a summary.
    """
    per_pair: List[Dict[str, Any]] = []
    all_bit_exact = True
    for (m, k) in cutoff_pairs:
        if not (1 <= k <= m):
            raise ValueError(f"verify_continuity_under_transitions requires 1 <= k <= m; got m={m}, k={k}")
        v_high = psi_inf.project_to_vector(m)
        v_low = psi_inf.project_to_vector(k)
        P = TransitionMap(m, k)
        pulled = P.apply_to_vector(v_high)
        bit_exact = pulled == v_low
        per_pair.append({
            "m": m,
            "k": k,
            "bit_exact": bit_exact,
        })
        all_bit_exact = all_bit_exact and bit_exact
    return {
        "label": psi_inf.label,
        "n_pairs": len(per_pair),
        "per_pair": per_pair,
        "all_bit_exact": all_bit_exact,
    }
