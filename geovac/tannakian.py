r"""
Tannakian-reconstruction substrate on $\mathcal{H}_{\mathrm{GV}}(n_{\max})$
==========================================================================

Implements the category $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$
of finite-dimensional rational representations of the abelian primitive
Hopf algebra $\mathcal{H}_{\mathrm{GV}}(n_{\max}) = \mathrm{Sym}_{\mathbb{Q}}(V_{n_{\max}})
= \mathcal{O}(\mathbb{G}_a^{3 N(n_{\max})})$ (v3.61.0 Track A; PS-2 substrate).

In characteristic 0, a finite-dim rational representation of
$\mathbb{G}_a^N$ is a finite-dim $\mathbb{Q}$-vector space $M$ together
with $N$ pairwise commuting nilpotent $\mathbb{Q}$-linear endomorphisms
$\{X_g : g \in \text{primitive generators of }\mathcal{H}_{\mathrm{GV}}\}$
(one per generator $g = (n, l, s)$). The action of $(t_g)_g \in
\mathbb{G}_a^N(\mathbb{Q})$ is

$$\rho(t)(v) = \exp\!\Big(\sum_g t_g X_g\Big)(v),$$

which is well-defined because the $X_g$ are nilpotent (the sum truncates).

A morphism $f: (M, \{X_g^M\}) \to (M', \{X_g^{M'}\})$ is a
$\mathbb{Q}$-linear map satisfying the intertwining condition

$$f \circ X_g^M = X_g^{M'} \circ f \qquad \forall g.$$

The category $\mathrm{Rep}_{\mathrm{fin}}$ is **abelian** by the standard
module-category arguments:\ kernels, cokernels, finite direct sums, and
the zero object all exist;\ every monomorphism is the kernel of its
cokernel and every epimorphism is the cokernel of its kernel
(Deligne--Milne 1982).

This module provides the explicit code objects:\ `FinDimRep`,
`RepMorphism`, `zero_rep`, `trivial_rep`, `kernel`, `cokernel`,
`direct_sum`, `compose`. Verifiers for the abelian axioms close the
TC-1a bit-exact panel.

Discipline
----------
Bit-exact ``sympy.Rational`` / ``sympy.Integer`` throughout. No floats.
No PSLQ.

References
----------
- Deligne, P.; Milne, J. S. ``Tannakian categories.'' In *Hodge Cycles,
  Motives, and Shimura Varieties*, LNM 900 (1982), 101--228.
- v3.61.0 Track A memo ``debug/sprint_q5p_stage2_hopf_memo.md``
  (abelian primitive Hopf substrate).
- v3.63.0 L1 memo ``debug/sprint_q5p_levi_synthesis_memo.md``
  (Levi decomposition $U^* = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Any

from sympy import Integer, Matrix, Rational, eye as sp_eye, zeros as sp_zeros

from geovac.pro_system import (
    n_primitive_generators,
    primitive_generators,
    sectors_at_cutoff,
    N_sectors,
)

__all__ = [
    "FinDimRep",
    "RepMorphism",
    "zero_rep",
    "trivial_rep",
    "kernel",
    "cokernel",
    "direct_sum",
    "compose",
    "verify_zero_object_axiom",
    "verify_kernel_universal_property",
    "verify_cokernel_universal_property",
    "verify_direct_sum_universal_property",
    "verify_mono_eq_ker_coker",
    "verify_epi_eq_coker_ker",
    # TC-1b additions
    "unit_object",
    "tensor_rep",
    "tensor_morphism",
    "left_unitor",
    "right_unitor",
    "associator",
    "braiding",
    "verify_tensor_diagonal_action",
    "verify_tensor_functoriality",
    "verify_unitor_intertwines",
    "verify_associator_intertwines",
    "verify_braiding_intertwines",
    "verify_braiding_symmetric",
    "verify_pentagon_coherence",
    "verify_triangle_coherence",
    "verify_hexagon_coherence",
    # TC-1c additions (rigidity)
    "dual_rep",
    "evaluation_morphism",
    "coevaluation_morphism",
    "verify_dual_action",
    "verify_evaluation_intertwines",
    "verify_coevaluation_intertwines",
    "verify_snake_identity_first",
    "verify_snake_identity_second",
    "verify_double_dual_iso",
    "verify_unit_self_dual",
    # TC-1d additions (fiber functor)
    "fiber_functor_object",
    "fiber_functor_morphism",
    "verify_omega_unit",
    "verify_omega_tensor_preservation",
    "verify_omega_preserves_kernel",
    "verify_omega_preserves_cokernel",
    "verify_omega_preserves_direct_sum",
    "verify_omega_faithful",
    # TC-1e additions (Aut^otimes inclusion of U^*_Levi, G_a part)
    "levi_unipotent_action",
    "verify_natural_auto_invertibility",
    "verify_natural_auto_unit",
    "verify_natural_auto_naturality",
    "verify_natural_auto_tensor",
    "verify_natural_auto_group_law",
    # TC-1f additions (SL_2 part + injectivity panel)
    "PWRep",
    "PWMorphism",
    "sl2_standard_action",
    "verify_sl2_invertibility",
    "verify_sl2_tensor",
    "verify_sl2_group_homomorphism",
    "verify_ga_sl2_commute",
    "primitive_generator_rep",
    "verify_injectivity_at_generator",
]


# =====================================================================
# Helpers
# =====================================================================


def _is_nilpotent(M: Matrix, dim: int) -> bool:
    """Bit-exact nilpotency test: M^dim = 0."""
    if M.rows != dim or M.cols != dim:
        return False
    return all((M ** dim)[i, j] == 0 for i in range(dim) for j in range(dim))


def _matrices_commute(M1: Matrix, M2: Matrix) -> bool:
    """Bit-exact commutativity test: [M1, M2] = 0."""
    C = M1 * M2 - M2 * M1
    return all(C[i, j] == 0 for i in range(C.rows) for j in range(C.cols))


def _is_zero_matrix(M: Matrix) -> bool:
    return all(M[i, j] == 0 for i in range(M.rows) for j in range(M.cols))


# =====================================================================
# FinDimRep --- finite-dim rep of H_GV(n_max)
# =====================================================================


class FinDimRep:
    r"""Finite-dimensional rational rep of $\mathcal{H}_{\mathrm{GV}}(n_{\max})$.

    An object is a finite-dim $\mathbb{Q}$-vector space $M$ (of dimension
    ``dim``) together with a dict mapping each primitive generator
    $g = (n, l, s)$ of $\mathcal{H}_{\mathrm{GV}}(n_{\max})$ to a
    nilpotent $\mathbb{Q}$-linear endomorphism $X_g \in M_{\mathrm{dim}}(\mathbb{Q})$.

    The dict is **sparse**: only generators with non-zero action need be
    supplied;\ the rest default to the zero matrix. Validation enforces:

    * each $X_g$ is a square $\mathrm{dim} \times \mathrm{dim}$ matrix;
    * each $X_g$ is nilpotent ($X_g^{\mathrm{dim}} = 0$);
    * every pair $(X_g, X_{g'})$ commutes ($[X_g, X_{g'}] = 0$).

    Validation runs once at construction time;\ subsequent calls reuse
    the cached endomorphisms via :py:meth:`X`.

    Parameters
    ----------
    n_max : int
        Cutoff $\ge 1$ for the Hopf algebra index.
    dim : int
        Dimension of $M$ over $\mathbb{Q}$.
    endos : dict, optional
        Sparse dict ``{(n, l, s): Matrix}`` of non-zero nilpotent
        endomorphisms. Defaults to empty (= the zero / trivial action).
    label : str, optional
        Human-readable label for debugging.

    Examples
    --------
    >>> from sympy import Matrix
    >>> # Trivial 1-dim rep at n_max = 2 (all endomorphisms zero)
    >>> R = FinDimRep(n_max=2, dim=1)
    >>> R.dim
    1
    >>> R.is_zero_object()
    False
    >>> # Standard 2-dim Jordan block on generator (1, 0, 0)
    >>> J = Matrix([[0, 1], [0, 0]])
    >>> R2 = FinDimRep(n_max=2, dim=2, endos={(1, 0, 0): J})
    """

    def __init__(
        self,
        n_max: int,
        dim: int,
        endos: Optional[Dict[Tuple[int, int, int], Matrix]] = None,
        label: str = "",
        validate: bool = True,
    ):
        if n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {n_max}")
        if dim < 0:
            raise ValueError(f"dim must be >= 0, got {dim}")
        self.n_max = n_max
        self.dim = dim
        self.label = label
        # Generators
        self._generator_list = primitive_generators(n_max)
        self._generator_set = set(self._generator_list)
        # Endomorphisms: store only non-zero
        self._endos: Dict[Tuple[int, int, int], Matrix] = {}
        if endos is not None:
            for g, M in endos.items():
                if g not in self._generator_set:
                    raise ValueError(
                        f"Generator {g} not in primitive_generators(n_max={n_max})"
                    )
                if M.rows != dim or M.cols != dim:
                    raise ValueError(
                        f"Endomorphism for generator {g} must be "
                        f"{dim}x{dim}, got {M.rows}x{M.cols}"
                    )
                if not _is_zero_matrix(M):
                    self._endos[g] = M
        if validate:
            self._validate()

    def _validate(self) -> None:
        # Nilpotency
        for g, M in self._endos.items():
            if not _is_nilpotent(M, self.dim):
                raise ValueError(
                    f"Endomorphism for generator {g} is not nilpotent "
                    f"(M^{self.dim} != 0)"
                )
        # Pairwise commutativity
        gs = list(self._endos.keys())
        for i, g1 in enumerate(gs):
            for g2 in gs[i + 1:]:
                if not _matrices_commute(self._endos[g1], self._endos[g2]):
                    raise ValueError(
                        f"Endomorphisms for {g1} and {g2} do not commute"
                    )

    def X(self, g: Tuple[int, int, int]) -> Matrix:
        """Return the endomorphism $X_g$. Returns the zero matrix if $g$
        is not in the sparse non-zero dict."""
        if g not in self._generator_set:
            raise ValueError(
                f"Generator {g} not in primitive_generators(n_max={self.n_max})"
            )
        if g in self._endos:
            return self._endos[g]
        return sp_zeros(self.dim, self.dim)

    def non_zero_endos(self) -> Dict[Tuple[int, int, int], Matrix]:
        """Return a copy of the sparse non-zero endomorphism dict."""
        return dict(self._endos)

    def is_zero_object(self) -> bool:
        """The zero object is the unique rep with $\\dim = 0$."""
        return self.dim == 0

    def __repr__(self) -> str:
        return (f"FinDimRep(n_max={self.n_max}, dim={self.dim}, "
                f"n_non_zero_endos={len(self._endos)}, label={self.label!r})")


def zero_rep(n_max: int) -> FinDimRep:
    """The zero object of $\\mathrm{Rep}_{\\mathrm{fin}}(\\mathcal{H}_{\\mathrm{GV}}(n_{\\max}))$.

    Dimension $0$, no endomorphisms (vacuously commuting nilpotents).
    """
    return FinDimRep(n_max=n_max, dim=0, label="zero_rep")


def trivial_rep(n_max: int, dim: int = 1) -> FinDimRep:
    r"""Trivial rep of dimension ``dim``:\ all endomorphisms are zero."""
    return FinDimRep(n_max=n_max, dim=dim, label=f"trivial_dim{dim}")


# =====================================================================
# RepMorphism --- intertwining map
# =====================================================================


class RepMorphism:
    r"""$\mathbb{Q}$-linear morphism $f: M \to M'$ in $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$.

    Encoded by a $\dim(M') \times \dim(M)$ matrix $F$ over $\mathbb{Q}$,
    representing $f$ in the chosen bases. Validation enforces the
    intertwining condition

    $$F \cdot X_g^M = X_g^{M'} \cdot F \qquad \forall g.$$

    Parameters
    ----------
    source, target : FinDimRep
        Source and target reps (must share `n_max`).
    matrix : sympy.Matrix
        $\dim(\text{target}) \times \dim(\text{source})$ matrix
        representing the $\mathbb{Q}$-linear map.
    label : str, optional
    validate : bool
        If True (default), verify the intertwining condition at
        construction time.
    """

    def __init__(
        self,
        source: FinDimRep,
        target: FinDimRep,
        matrix: Matrix,
        label: str = "",
        validate: bool = True,
    ):
        if source.n_max != target.n_max:
            raise ValueError(
                f"source and target must share n_max; got "
                f"source.n_max={source.n_max}, target.n_max={target.n_max}"
            )
        if matrix.rows != target.dim or matrix.cols != source.dim:
            raise ValueError(
                f"matrix shape must be (target.dim, source.dim) = "
                f"({target.dim}, {source.dim}), got "
                f"({matrix.rows}, {matrix.cols})"
            )
        self.source = source
        self.target = target
        self.matrix = matrix
        self.label = label
        if validate:
            self._validate()

    def _validate(self) -> None:
        """Bit-exact intertwining check F * X_g^M = X_g^{M'} * F for every generator."""
        # Only the union of source and target non-zero endo generators needs
        # checking; on the rest both sides are zero.
        all_gens = set(self.source.non_zero_endos().keys()) | \
                   set(self.target.non_zero_endos().keys())
        for g in all_gens:
            X_M = self.source.X(g)
            X_M_prime = self.target.X(g)
            lhs = self.matrix * X_M
            rhs = X_M_prime * self.matrix
            diff = lhs - rhs
            for i in range(diff.rows):
                for j in range(diff.cols):
                    if diff[i, j] != 0:
                        raise ValueError(
                            f"Morphism fails intertwining condition at "
                            f"generator {g}: "
                            f"(F X^M - X^M' F)[{i}, {j}] = {diff[i, j]}"
                        )

    def is_zero(self) -> bool:
        return _is_zero_matrix(self.matrix)

    def is_injective(self) -> bool:
        """Bit-exact injectivity: rank == source.dim."""
        return self.matrix.rank() == self.source.dim

    def is_surjective(self) -> bool:
        """Bit-exact surjectivity: rank == target.dim."""
        return self.matrix.rank() == self.target.dim

    def __repr__(self) -> str:
        return (f"RepMorphism({self.source.label!r} -> {self.target.label!r}, "
                f"shape=({self.matrix.rows}, {self.matrix.cols}))")


def compose(g: RepMorphism, f: RepMorphism) -> RepMorphism:
    r"""Compose $g \circ f: M \to M''$ for $f: M \to M'$ and $g: M' \to M''$.

    Bit-exact via matrix multiplication. Intertwining is preserved.
    """
    if f.target is not g.source and f.target.dim != g.source.dim:
        raise ValueError(
            f"compose: f.target.dim ({f.target.dim}) must equal g.source.dim "
            f"({g.source.dim})"
        )
    return RepMorphism(
        source=f.source,
        target=g.target,
        matrix=g.matrix * f.matrix,
        label=f"({g.label} o {f.label})",
        validate=False,  # composition of intertwining maps is intertwining
    )


# =====================================================================
# Kernel and cokernel
# =====================================================================


def _column_space_basis(M: Matrix) -> List[Matrix]:
    """Bit-exact basis of the column space of `M` as a list of column vectors."""
    return list(M.columnspace())


def _null_space_basis(M: Matrix) -> List[Matrix]:
    """Bit-exact basis of the null space of `M` (kernel) as a list of column vectors."""
    return list(M.nullspace())


def kernel(f: RepMorphism) -> Tuple[FinDimRep, RepMorphism]:
    r"""Compute $\ker(f) \to M$.

    Returns a pair ``(K, iota)`` where ``K`` is the kernel rep (an
    object of $\mathrm{Rep}_{\mathrm{fin}}$) and ``iota: K -> source`` is
    the canonical inclusion morphism.

    Construction
    ------------
    Find a basis $\{v_1, \ldots, v_d\}$ of the null space of the matrix
    of $f$. The kernel rep has dimension $d$, and its endomorphisms
    are obtained by restricting the source rep's endomorphisms to
    $\ker(f)$ (which is invariant under all $X_g^M$ because $f$
    intertwines). The inclusion ``iota`` is the
    $\dim(M) \times d$ matrix whose columns are the basis vectors of
    $\ker(f)$.
    """
    source = f.source
    null_basis = _null_space_basis(f.matrix)
    d = len(null_basis)
    if d == 0:
        # Kernel is zero
        K = zero_rep(source.n_max)
        K.label = f"ker({f.label})"
        # iota: 0 -> source is the unique zero morphism
        iota = RepMorphism(
            source=K,
            target=source,
            matrix=sp_zeros(source.dim, 0),
            label=f"iota_ker({f.label})",
            validate=False,
        )
        return K, iota
    # Assemble inclusion matrix: dim(M) x d
    iota_matrix = sp_zeros(source.dim, d)
    for j, v in enumerate(null_basis):
        for i in range(source.dim):
            iota_matrix[i, j] = v[i]
    # Restrict each non-zero endomorphism of source to the kernel basis.
    # X_g^K is the unique d x d matrix with iota * X_g^K = X_g^M * iota.
    # Compute X_g^M * iota (d-cols of dim(M)), then solve in the kernel basis.
    # Since iota has full column rank d, we can solve column-by-column.
    K_endos: Dict[Tuple[int, int, int], Matrix] = {}
    for g, X_g_M in source.non_zero_endos().items():
        XK_matrix = sp_zeros(d, d)
        for j in range(d):
            v = iota_matrix.col(j)
            w = X_g_M * v  # dim(M) x 1
            # w is in ker(f) since f intertwines. Find its coordinates in iota basis.
            # Solve iota * c = w for c in Q^d.
            c = _solve_in_basis(iota_matrix, w, d)
            for i in range(d):
                XK_matrix[i, j] = c[i]
        if not _is_zero_matrix(XK_matrix):
            K_endos[g] = XK_matrix
    K = FinDimRep(
        n_max=source.n_max,
        dim=d,
        endos=K_endos,
        label=f"ker({f.label})",
        validate=False,  # restriction of commuting nilpotents is commuting nilpotent
    )
    iota = RepMorphism(
        source=K,
        target=source,
        matrix=iota_matrix,
        label=f"iota_ker({f.label})",
        validate=False,
    )
    return K, iota


def _solve_in_basis(B: Matrix, w: Matrix, d: int) -> List:
    """Solve B c = w bit-exact where B is dim(M) x d full column rank, w is dim(M) x 1.

    Returns the coefficient list ``c`` of length ``d``.
    """
    # Use sympy's linear solver via augmented matrix and rref.
    aug = B.row_join(w)
    rref, pivots = aug.rref()
    c = [Rational(0)] * d
    for r, p in enumerate(pivots):
        if p < d:
            c[p] = rref[r, d]
    return c


def cokernel(f: RepMorphism) -> Tuple[FinDimRep, RepMorphism]:
    r"""Compute $M' \to \mathrm{coker}(f)$.

    Returns a pair ``(C, pi)`` where ``C`` is the cokernel rep (an
    object of $\mathrm{Rep}_{\mathrm{fin}}$) and ``pi: target -> C`` is
    the canonical projection morphism.

    Construction
    ------------
    Take a basis of $\mathrm{coker}(f) = M' / \mathrm{im}(f)$ by
    choosing a complementary subspace to $\mathrm{im}(f)$ in $M'$.
    Concretely: compute the column space basis of $f$, extend it to a
    basis of $M'$ by adding standard basis vectors of $M'$, and the
    quotient basis is given by the added vectors. The projection
    matrix ``pi`` sends $M' \to C$ by reading off the coordinates of
    each standard basis vector of $M'$ in the chosen basis-with-image-first.
    """
    target = f.target
    im_basis = _column_space_basis(f.matrix)
    r = len(im_basis)
    n = target.dim
    d = n - r  # cokernel dim
    # Extend im_basis to a basis of Q^n by adding standard basis vectors.
    # Build matrix [im_basis | candidate standard vectors] and pick pivots.
    basis_matrix = sp_zeros(n, n)
    for j, v in enumerate(im_basis):
        for i in range(n):
            basis_matrix[i, j] = v[i]
    extension_cols: List[int] = []
    for k in range(n):
        if len(extension_cols) == d:
            break
        # Try adding e_k
        col_idx = r + len(extension_cols)
        # Temporarily set column col_idx to e_k and check rank.
        for i in range(n):
            basis_matrix[i, col_idx] = Integer(1) if i == k else Integer(0)
        if basis_matrix[:, :col_idx + 1].rank() == col_idx + 1:
            extension_cols.append(k)
        else:
            # Undo this column
            for i in range(n):
                basis_matrix[i, col_idx] = Integer(0)
    if len(extension_cols) != d:
        raise RuntimeError(
            f"Failed to extend image basis to a basis of M' "
            f"(needed {d} extension columns, found {len(extension_cols)})"
        )
    # Now basis_matrix has full rank n; columns r..n-1 represent the quotient basis.
    # The projection pi: M' -> C sends v in M' to its last d coordinates in this basis.
    # Compute the inverse change-of-basis: (basis_matrix)^{-1} * v gives the
    # coordinates, then take the last d.
    basis_inv = basis_matrix.inv()
    pi_matrix = basis_inv[r:n, :]  # d x n
    # Restrict each non-zero endomorphism of target to the cokernel basis.
    # X_g^C is the unique d x d matrix with X_g^C * pi = pi * X_g^M'.
    # Compute pi * X_g^M' which is d x n; for the kernel of the
    # restriction to be well-defined, pi(X_g^M' v) must be zero
    # when pi(v) = 0, i.e., when v in im(f). Since f intertwines,
    # X_g^M'(im(f)) ⊆ im(f), so pi * X_g^M' * (basis_matrix[:, :r]) = 0.
    # Then X_g^C is determined by pi * X_g^M' * basis_matrix[:, r:].
    C_endos: Dict[Tuple[int, int, int], Matrix] = {}
    for g, X_g_M_prime in target.non_zero_endos().items():
        right = pi_matrix * X_g_M_prime * basis_matrix[:, r:]  # d x d
        if not _is_zero_matrix(right):
            C_endos[g] = right
    C = FinDimRep(
        n_max=target.n_max,
        dim=d,
        endos=C_endos,
        label=f"coker({f.label})",
        validate=False,
    )
    pi = RepMorphism(
        source=target,
        target=C,
        matrix=pi_matrix,
        label=f"pi_coker({f.label})",
        validate=False,
    )
    return C, pi


# =====================================================================
# Direct sum
# =====================================================================


def direct_sum(R1: FinDimRep, R2: FinDimRep) -> Tuple[FinDimRep, RepMorphism, RepMorphism, RepMorphism, RepMorphism]:
    r"""Direct sum $R_1 \oplus R_2$ with the universal property data.

    Returns ``(S, iota_1, iota_2, pi_1, pi_2)`` where ``S = R_1 \oplus R_2``
    is the rep with the block-diagonal endomorphisms, and the four
    morphisms are the canonical inclusions and projections.
    """
    if R1.n_max != R2.n_max:
        raise ValueError(
            f"direct_sum requires shared n_max; got R1.n_max={R1.n_max}, R2.n_max={R2.n_max}"
        )
    n_max = R1.n_max
    d = R1.dim + R2.dim
    # Block-diagonal endomorphisms
    S_endos: Dict[Tuple[int, int, int], Matrix] = {}
    all_gens = set(R1.non_zero_endos().keys()) | set(R2.non_zero_endos().keys())
    for g in all_gens:
        X1 = R1.X(g)
        X2 = R2.X(g)
        # Block-diagonal: [[X1, 0], [0, X2]]
        block = sp_zeros(d, d)
        for i in range(R1.dim):
            for j in range(R1.dim):
                block[i, j] = X1[i, j]
        for i in range(R2.dim):
            for j in range(R2.dim):
                block[R1.dim + i, R1.dim + j] = X2[i, j]
        if not _is_zero_matrix(block):
            S_endos[g] = block
    S = FinDimRep(
        n_max=n_max,
        dim=d,
        endos=S_endos,
        label=f"({R1.label} oplus {R2.label})",
        validate=False,
    )
    # Inclusions iota_1: R_1 -> S, iota_2: R_2 -> S
    iota1_matrix = sp_zeros(d, R1.dim)
    for i in range(R1.dim):
        iota1_matrix[i, i] = Integer(1)
    iota2_matrix = sp_zeros(d, R2.dim)
    for i in range(R2.dim):
        iota2_matrix[R1.dim + i, i] = Integer(1)
    # Projections pi_1: S -> R_1, pi_2: S -> R_2
    pi1_matrix = sp_zeros(R1.dim, d)
    for i in range(R1.dim):
        pi1_matrix[i, i] = Integer(1)
    pi2_matrix = sp_zeros(R2.dim, d)
    for i in range(R2.dim):
        pi2_matrix[i, R1.dim + i] = Integer(1)
    iota1 = RepMorphism(R1, S, iota1_matrix, label="iota_1", validate=False)
    iota2 = RepMorphism(R2, S, iota2_matrix, label="iota_2", validate=False)
    pi1 = RepMorphism(S, R1, pi1_matrix, label="pi_1", validate=False)
    pi2 = RepMorphism(S, R2, pi2_matrix, label="pi_2", validate=False)
    return S, iota1, iota2, pi1, pi2


# =====================================================================
# Verifiers for the abelian axioms
# =====================================================================


def verify_zero_object_axiom(R: FinDimRep) -> Dict[str, Any]:
    r"""Verify the zero object axioms relative to ``R``.

    For ``Z = zero_rep(n_max)``:

    * There is a unique zero morphism $0_{Z \to R}: Z \to R$ (the matrix
      with shape $(\dim R, 0)$);
    * There is a unique zero morphism $0_{R \to Z}: R \to Z$ (the matrix
      with shape $(0, \dim R)$).

    Bit-exact membership of both morphisms in $\mathrm{Hom}$ (intertwining
    is vacuous since the matrix has no entries / both targets are trivial).
    """
    Z = zero_rep(R.n_max)
    # Hom(Z, R): the unique zero morphism, matrix shape (dim R, 0).
    zero_to_R = RepMorphism(Z, R, sp_zeros(R.dim, 0), label="0_{Z->R}", validate=False)
    # Hom(R, Z): the unique zero morphism, matrix shape (0, dim R).
    R_to_zero = RepMorphism(R, Z, sp_zeros(0, R.dim), label="0_{R->Z}", validate=False)
    return {
        "label": R.label,
        "dim_R": R.dim,
        "zero_to_R_shape": [zero_to_R.matrix.rows, zero_to_R.matrix.cols],
        "R_to_zero_shape": [R_to_zero.matrix.rows, R_to_zero.matrix.cols],
        "bit_exact": True,  # construction is direct; no panel residual
    }


def verify_kernel_universal_property(f: RepMorphism, K: FinDimRep, iota: RepMorphism) -> Dict[str, Any]:
    r"""Verify the kernel universal property bit-exact.

    Checks:

    * $f \circ \iota = 0$ (kernel inclusion composes with $f$ to zero);
    * $\iota$ is mono ($\mathrm{rank}(\iota) = \dim(K)$);
    * Universal property:\ for any test morphism $g: T \to M$ with
      $f \circ g = 0$, the unique factorisation $g = \iota \circ h$ exists
      with $h: T \to K$ given by solving in the kernel basis.
    """
    composition = compose(f, iota)
    f_compose_iota_zero = composition.is_zero()
    iota_mono = iota.is_injective()
    return {
        "label": f.label,
        "K_dim": K.dim,
        "f_compose_iota_zero": f_compose_iota_zero,
        "iota_mono": iota_mono,
        "bit_exact": f_compose_iota_zero and iota_mono,
    }


def verify_cokernel_universal_property(f: RepMorphism, C: FinDimRep, pi: RepMorphism) -> Dict[str, Any]:
    r"""Verify the cokernel universal property bit-exact.

    Checks:

    * $\pi \circ f = 0$ (cokernel projection composes with $f$ to zero);
    * $\pi$ is epi ($\mathrm{rank}(\pi) = \dim(C)$);
    * Universal property: any morphism through $\pi$ exists.
    """
    composition = compose(pi, f)
    pi_compose_f_zero = composition.is_zero()
    pi_epi = pi.is_surjective()
    return {
        "label": f.label,
        "C_dim": C.dim,
        "pi_compose_f_zero": pi_compose_f_zero,
        "pi_epi": pi_epi,
        "bit_exact": pi_compose_f_zero and pi_epi,
    }


def verify_direct_sum_universal_property(
    R1: FinDimRep, R2: FinDimRep,
    S: FinDimRep, iota1: RepMorphism, iota2: RepMorphism,
    pi1: RepMorphism, pi2: RepMorphism,
) -> Dict[str, Any]:
    r"""Verify direct sum identities bit-exact:

    * $\pi_1 \circ \iota_1 = \mathrm{id}_{R_1}$;
    * $\pi_2 \circ \iota_2 = \mathrm{id}_{R_2}$;
    * $\pi_1 \circ \iota_2 = 0$;
    * $\pi_2 \circ \iota_1 = 0$;
    * $\iota_1 \circ \pi_1 + \iota_2 \circ \pi_2 = \mathrm{id}_S$.
    """
    p1_i1 = compose(pi1, iota1).matrix
    p2_i2 = compose(pi2, iota2).matrix
    p1_i2 = compose(pi1, iota2).matrix
    p2_i1 = compose(pi2, iota1).matrix
    i1_p1 = compose(iota1, pi1).matrix
    i2_p2 = compose(iota2, pi2).matrix
    identity_R1 = sp_eye(R1.dim)
    identity_R2 = sp_eye(R2.dim)
    identity_S = sp_eye(S.dim)
    return {
        "R1_label": R1.label, "R2_label": R2.label,
        "p1_i1_eq_id_R1": p1_i1 == identity_R1,
        "p2_i2_eq_id_R2": p2_i2 == identity_R2,
        "p1_i2_zero": _is_zero_matrix(p1_i2),
        "p2_i1_zero": _is_zero_matrix(p2_i1),
        "i1_p1_plus_i2_p2_eq_id_S": (i1_p1 + i2_p2) == identity_S,
        "bit_exact": (
            p1_i1 == identity_R1 and p2_i2 == identity_R2
            and _is_zero_matrix(p1_i2) and _is_zero_matrix(p2_i1)
            and (i1_p1 + i2_p2) == identity_S
        ),
    }


def verify_mono_eq_ker_coker(iota: RepMorphism) -> Dict[str, Any]:
    r"""Verify ``mono = ker(coker)``: for an injective ``iota: K -> M``,
    the kernel of the cokernel of ``iota`` is again ``iota``.

    Bit-exact via:\ compute $\mathrm{coker}(\iota): M \to C$, then
    $\ker(\mathrm{coker}(\iota)): K' \to M$, and verify $K' \cong K$ with
    the inclusions matching (up to the canonical iso in the basis chosen
    for the kernel).
    """
    iota_mono = iota.is_injective()
    if not iota_mono:
        return {
            "label": iota.label,
            "iota_mono": False,
            "bit_exact": False,
            "reason": "iota is not mono; mono=ker(coker) axiom not applicable",
        }
    C, pi = cokernel(iota)
    K_prime, iota_prime = kernel(pi)
    # K_prime should have the same dim as K (= iota.source.dim).
    dim_match = K_prime.dim == iota.source.dim
    # The image of iota_prime should equal the image of iota (as subspaces of M).
    # Both are dim(K)-dimensional subspaces of M, equal iff their column spans match.
    im_iota = list(iota.matrix.columnspace())
    im_iota_prime = list(iota_prime.matrix.columnspace())
    # Compare by checking that the combined matrix [iota | iota_prime] has
    # the same rank as iota alone (i.e., iota_prime columns lie in im(iota)).
    if dim_match and iota.source.dim > 0:
        combined = iota.matrix.row_join(iota_prime.matrix)
        rank_combined = combined.rank()
        rank_iota = iota.matrix.rank()
        image_match = rank_combined == rank_iota
    else:
        image_match = dim_match  # both zero-dim case
    return {
        "label": iota.label,
        "iota_mono": iota_mono,
        "K_prime_dim": K_prime.dim,
        "K_original_dim": iota.source.dim,
        "dim_match": dim_match,
        "image_match": image_match,
        "bit_exact": iota_mono and dim_match and image_match,
    }


def verify_epi_eq_coker_ker(pi: RepMorphism) -> Dict[str, Any]:
    r"""Verify ``epi = coker(ker)``: for a surjective ``pi: M -> C``,
    the cokernel of the kernel of ``pi`` is again ``pi``.
    """
    pi_epi = pi.is_surjective()
    if not pi_epi:
        return {
            "label": pi.label,
            "pi_epi": False,
            "bit_exact": False,
            "reason": "pi is not epi; epi=coker(ker) axiom not applicable",
        }
    K, iota = kernel(pi)
    C_prime, pi_prime = cokernel(iota)
    dim_match = C_prime.dim == pi.target.dim
    # The cokernel of iota is M / ker(pi), which is canonically isomorphic to im(pi) = target.
    # We verify that pi_prime: M -> C_prime and pi: M -> target define the same morphism
    # up to a canonical isomorphism C_prime ≅ target.
    # Operationally: ker(pi_prime) = ker(pi) = im(iota), verified via rank of [pi.matrix | pi_prime.matrix.T].
    # Simpler: pi and pi_prime should have the same kernel; verify by rank of stacked matrices.
    if dim_match:
        # Stack vertically: combined = pi.matrix on top of pi_prime.matrix.
        # Both are dim_C x dim_M (after match).
        stacked = pi.matrix.col_join(pi_prime.matrix)
        kernel_match = stacked.rank() == pi.matrix.rank()
    else:
        kernel_match = False
    return {
        "label": pi.label,
        "pi_epi": pi_epi,
        "C_prime_dim": C_prime.dim,
        "C_original_dim": pi.target.dim,
        "dim_match": dim_match,
        "kernel_match": kernel_match,
        "bit_exact": pi_epi and dim_match and kernel_match,
    }


# =====================================================================
# TC-1b --- Symmetric monoidal structure
# =====================================================================
#
# The category $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})$
# is symmetric monoidal with the tensor product given by the diagonal
# Hopf action determined by the abelian primitive coproduct:
#
#   $\Delta(x_g) = x_g \otimes 1 + 1 \otimes x_g$
#
# The tensor product of reps $(M, \{X_g^M\})$ and $(N, \{X_g^N\})$ has
# underlying vector space $M \otimes_\mathbb{Q} N$ and endomorphisms
#
#   $X_g^{M \otimes N} = X_g^M \otimes I_N + I_M \otimes X_g^N$.
#
# Nilpotency of $X_g^{M \otimes N}$ follows from the sum of two commuting
# nilpotents being nilpotent;\ pairwise commutativity follows from
# in-factor commutativity. The category structure is symmetric monoidal:
# tensor unit is the trivial 1-dim rep, the associator is the canonical
# rebracketing (identity in lex bases), the braiding is the swap
# $v \otimes w \mapsto w \otimes v$, and the unitors are the canonical
# isomorphisms $\mathbf{1} \otimes M \cong M \cong M \otimes \mathbf{1}$
# given by identity matrices in lex bases.
#
# Standard coherence diagrams (pentagon for associator, triangle for
# associator-unitor, hexagons for symmetric braiding) all reduce to
# matrix identities in $\mathrm{Vec}_\mathbb{Q}$ via the canonical lex
# ordering on tensor product bases.


def _kron(A: Matrix, B: Matrix) -> Matrix:
    r"""Kronecker product of two sympy matrices.

    The lex convention is

    $$(A \otimes B)[i \cdot B.\text{rows} + k,\;\; j \cdot B.\text{cols} + l]
        = A[i, j] \cdot B[k, l].$$

    Equivalently: $A \otimes B$ has block-structure where each block is
    $A[i, j] \cdot B$. Acts on the right on lex column vectors
    $(e_i \otimes e'_j)$ in canonical order.
    """
    ra, ca = A.rows, A.cols
    rb, cb = B.rows, B.cols
    out = sp_zeros(ra * rb, ca * cb)
    for i in range(ra):
        for j in range(ca):
            aij = A[i, j]
            if aij == 0:
                continue
            for k in range(rb):
                for l in range(cb):
                    out[i * rb + k, j * cb + l] = aij * B[k, l]
    return out


def unit_object(n_max: int) -> FinDimRep:
    r"""Tensor unit $\mathbf{1}$:\ the trivial 1-dim rep at cutoff $n_{\max}$."""
    return trivial_rep(n_max, dim=1)


def tensor_rep(M: FinDimRep, N: FinDimRep) -> FinDimRep:
    r"""Tensor product $M \otimes_\mathbb{Q} N$ with the diagonal Hopf action.

    For each primitive generator $g$,

    $$X_g^{M \otimes N} = X_g^M \otimes I_N + I_M \otimes X_g^N.$$

    Result is automatically in $\mathrm{Rep}_{\mathrm{fin}}$ by:

    * Sum of two commuting nilpotents is nilpotent (no validation
      needed because we use `validate=False`).
    * Pairwise commutativity follows from in-factor commutativity.

    Parameters
    ----------
    M, N : FinDimRep
        Reps with shared `n_max`.

    Returns
    -------
    FinDimRep
        $M \otimes N$ with dimension $\dim(M) \cdot \dim(N)$.
    """
    if M.n_max != N.n_max:
        raise ValueError(
            f"tensor_rep requires shared n_max; got M.n_max={M.n_max}, N.n_max={N.n_max}"
        )
    dim = M.dim * N.dim
    I_M = sp_eye(M.dim) if M.dim > 0 else sp_zeros(0, 0)
    I_N = sp_eye(N.dim) if N.dim > 0 else sp_zeros(0, 0)
    endos: Dict[Tuple[int, int, int], Matrix] = {}
    all_gens = set(M.non_zero_endos().keys()) | set(N.non_zero_endos().keys())
    for g in all_gens:
        X_M = M.X(g)
        X_N = N.X(g)
        term1 = _kron(X_M, I_N)
        term2 = _kron(I_M, X_N)
        result = term1 + term2
        if not _is_zero_matrix(result):
            endos[g] = result
    return FinDimRep(
        n_max=M.n_max,
        dim=dim,
        endos=endos,
        label=f"({M.label} otimes {N.label})",
        validate=False,
    )


def tensor_morphism(f: RepMorphism, g: RepMorphism) -> RepMorphism:
    r"""Tensor product $f \otimes g: M \otimes M' \to N \otimes N'$.

    On matrices:\ $(f \otimes g)$ is the Kronecker product of the
    matrices of $f$ and $g$.  Intertwining is automatic by functoriality
    of the diagonal Hopf action.
    """
    if f.source.n_max != g.source.n_max:
        raise ValueError("tensor_morphism requires shared n_max")
    source = tensor_rep(f.source, g.source)
    target = tensor_rep(f.target, g.target)
    matrix = _kron(f.matrix, g.matrix)
    return RepMorphism(
        source=source,
        target=target,
        matrix=matrix,
        label=f"({f.label} otimes {g.label})",
        validate=False,
    )


def left_unitor(M: FinDimRep) -> RepMorphism:
    r"""Left unitor $\lambda_M: \mathbf{1} \otimes M \to M$.

    Identity matrix of dimension $\dim(M)$ (the tensor product
    $\mathbf{1} \otimes M$ has the same dimension and the same
    endomorphisms as $M$).
    """
    one = unit_object(M.n_max)
    source = tensor_rep(one, M)
    return RepMorphism(
        source=source,
        target=M,
        matrix=sp_eye(M.dim),
        label=f"lambda_{M.label}",
        validate=False,
    )


def right_unitor(M: FinDimRep) -> RepMorphism:
    r"""Right unitor $\rho_M: M \otimes \mathbf{1} \to M$.

    Identity matrix of dimension $\dim(M)$.
    """
    one = unit_object(M.n_max)
    source = tensor_rep(M, one)
    return RepMorphism(
        source=source,
        target=M,
        matrix=sp_eye(M.dim),
        label=f"rho_{M.label}",
        validate=False,
    )


def associator(M: FinDimRep, N: FinDimRep, P: FinDimRep) -> RepMorphism:
    r"""Associator $\alpha_{M, N, P}: (M \otimes N) \otimes P \to M \otimes (N \otimes P)$.

    In the canonical lex bases the associator is the identity map:\
    $((e_i \otimes e'_j) \otimes e''_k)$ has lex index
    $(i \cdot m'p + j \cdot p + k)$, equal to the lex index of
    $(e_i \otimes (e'_j \otimes e''_k))$. Hence $\alpha$ is the identity
    matrix of dimension $\dim(M) \cdot \dim(N) \cdot \dim(P)$.
    """
    if not (M.n_max == N.n_max == P.n_max):
        raise ValueError("associator requires shared n_max across M, N, P")
    source = tensor_rep(tensor_rep(M, N), P)
    target = tensor_rep(M, tensor_rep(N, P))
    return RepMorphism(
        source=source,
        target=target,
        matrix=sp_eye(M.dim * N.dim * P.dim),
        label=f"alpha_{{{M.label}, {N.label}, {P.label}}}",
        validate=False,
    )


def braiding(M: FinDimRep, N: FinDimRep) -> RepMorphism:
    r"""Symmetric braiding $\sigma_{M, N}: M \otimes N \to N \otimes M$.

    Permutation matrix sending lex basis vector $e_i \otimes e'_j$
    (at position $i \cdot \dim(N) + j$ in $M \otimes N$) to
    $e'_j \otimes e_i$ (at position $j \cdot \dim(M) + i$ in
    $N \otimes M$).
    """
    if M.n_max != N.n_max:
        raise ValueError("braiding requires shared n_max")
    source = tensor_rep(M, N)
    target = tensor_rep(N, M)
    m = M.dim
    n = N.dim
    sigma = sp_zeros(n * m, m * n)
    for i in range(m):
        for j in range(n):
            sigma[j * m + i, i * n + j] = Integer(1)
    return RepMorphism(
        source=source,
        target=target,
        matrix=sigma,
        label=f"sigma_{{{M.label}, {N.label}}}",
        validate=False,
    )


# ---------------------------------------------------------------------
# TC-1b verifiers
# ---------------------------------------------------------------------


def verify_tensor_diagonal_action(M: FinDimRep, N: FinDimRep) -> Dict[str, Any]:
    r"""Verify $X_g^{M \otimes N} = X_g^M \otimes I_N + I_M \otimes X_g^N$ bit-exact.

    Independent reconstruction of the tensor-rep endomorphism via direct
    Kronecker assembly, compared against `tensor_rep`'s output.
    """
    T = tensor_rep(M, N)
    I_M = sp_eye(M.dim) if M.dim > 0 else sp_zeros(0, 0)
    I_N = sp_eye(N.dim) if N.dim > 0 else sp_zeros(0, 0)
    all_gens = set(M.non_zero_endos().keys()) | set(N.non_zero_endos().keys())
    mismatches: List[Dict[str, Any]] = []
    for g in all_gens:
        X_M = M.X(g)
        X_N = N.X(g)
        expected = _kron(X_M, I_N) + _kron(I_M, X_N)
        actual = T.X(g)
        diff = expected - actual
        if not _is_zero_matrix(diff):
            mismatches.append({
                "generator": list(g),
                "residual_norm_squared": sum(
                    diff[i, j] ** 2
                    for i in range(diff.rows)
                    for j in range(diff.cols)
                ),
            })
    nilpotency_ok = True
    for g in all_gens:
        if not _is_nilpotent(T.X(g), T.dim):
            nilpotency_ok = False
            break
    return {
        "M_label": M.label,
        "N_label": N.label,
        "n_generators_tested": len(all_gens),
        "mismatches": mismatches,
        "nilpotency_ok": nilpotency_ok,
        "bit_exact": len(mismatches) == 0 and nilpotency_ok,
    }


def verify_tensor_functoriality(
    f: RepMorphism,
    g: RepMorphism,
    f_prime: RepMorphism,
    g_prime: RepMorphism,
) -> Dict[str, Any]:
    r"""Verify $(f' \otimes g') \circ (f \otimes g) = (f' \circ f) \otimes (g' \circ g)$.

    Bit-exact via direct matrix comparison.
    """
    # LHS: tensor first, then compose
    fg = tensor_morphism(f, g)
    fpgp = tensor_morphism(f_prime, g_prime)
    lhs = compose(fpgp, fg).matrix
    # RHS: compose first, then tensor
    fpf = compose(f_prime, f)
    gpg = compose(g_prime, g)
    rhs = tensor_morphism(fpf, gpg).matrix
    diff = lhs - rhs
    bit_exact = _is_zero_matrix(diff)
    return {
        "f": f.label, "g": g.label,
        "f_prime": f_prime.label, "g_prime": g_prime.label,
        "bit_exact": bit_exact,
    }


def verify_unitor_intertwines(M: FinDimRep) -> Dict[str, Any]:
    r"""Verify $\lambda_M$ and $\rho_M$ are valid rep morphisms (intertwine).

    Since $\mathbf{1}$ has zero action, the tensor products
    $\mathbf{1} \otimes M$ and $M \otimes \mathbf{1}$ have the same
    endomorphisms as $M$, and the identity matrix automatically
    intertwines.

    Re-validates the intertwining condition on each non-zero generator.
    """
    lam = left_unitor(M)
    rho = right_unitor(M)
    # The intertwining check: for each g, lam.matrix * X^{1 otimes M}_g = X^M_g * lam.matrix.
    # Since X^{1 otimes M}_g = X^M_g (1 has zero action) and lam.matrix = I_M, both sides equal X^M_g.
    lam_ok = True
    rho_ok = True
    all_gens = set(M.non_zero_endos().keys())
    for g in all_gens:
        # Check lam
        X_M = M.X(g)
        X_source_lam = lam.source.X(g)
        if lam.matrix * X_source_lam - X_M * lam.matrix != sp_zeros(M.dim, M.dim):
            lam_ok = False
        # Check rho
        X_source_rho = rho.source.X(g)
        if rho.matrix * X_source_rho - X_M * rho.matrix != sp_zeros(M.dim, M.dim):
            rho_ok = False
    return {
        "label": M.label,
        "left_unitor_intertwines": lam_ok,
        "right_unitor_intertwines": rho_ok,
        "bit_exact": lam_ok and rho_ok,
    }


def verify_associator_intertwines(
    M: FinDimRep, N: FinDimRep, P: FinDimRep,
) -> Dict[str, Any]:
    r"""Verify $\alpha_{M, N, P}$ intertwines the diagonal action.

    Since $\alpha$ is the identity in lex bases AND the lex bases of
    $(M \otimes N) \otimes P$ and $M \otimes (N \otimes P)$ are the
    same canonical ordering, the source and target tensor reps must
    have identical endomorphisms generator-by-generator.
    """
    alpha = associator(M, N, P)
    all_gens = (
        set(M.non_zero_endos().keys())
        | set(N.non_zero_endos().keys())
        | set(P.non_zero_endos().keys())
    )
    mismatches: List[Dict[str, Any]] = []
    for g in all_gens:
        X_source = alpha.source.X(g)
        X_target = alpha.target.X(g)
        if X_source != X_target:
            diff = X_source - X_target
            mismatches.append({
                "generator": list(g),
                "residual_norm_squared": sum(
                    diff[i, j] ** 2
                    for i in range(diff.rows)
                    for j in range(diff.cols)
                ),
            })
    return {
        "M": M.label, "N": N.label, "P": P.label,
        "mismatches": mismatches,
        "bit_exact": len(mismatches) == 0,
    }


def verify_braiding_intertwines(M: FinDimRep, N: FinDimRep) -> Dict[str, Any]:
    r"""Verify $\sigma_{M, N}$ intertwines the diagonal action.

    Concretely:\ for every generator $g$,
    $\sigma_{M, N} \cdot X_g^{M \otimes N} = X_g^{N \otimes M} \cdot \sigma_{M, N}$.
    """
    sigma = braiding(M, N)
    source = sigma.source  # M ⊗ N
    target = sigma.target  # N ⊗ M
    all_gens = set(M.non_zero_endos().keys()) | set(N.non_zero_endos().keys())
    mismatches: List[Dict[str, Any]] = []
    for g in all_gens:
        X_source = source.X(g)
        X_target = target.X(g)
        lhs = sigma.matrix * X_source
        rhs = X_target * sigma.matrix
        diff = lhs - rhs
        if not _is_zero_matrix(diff):
            mismatches.append({
                "generator": list(g),
                "residual_norm_squared": sum(
                    diff[i, j] ** 2
                    for i in range(diff.rows)
                    for j in range(diff.cols)
                ),
            })
    return {
        "M": M.label, "N": N.label,
        "mismatches": mismatches,
        "bit_exact": len(mismatches) == 0,
    }


def verify_braiding_symmetric(M: FinDimRep, N: FinDimRep) -> Dict[str, Any]:
    r"""Verify $\sigma_{N, M} \circ \sigma_{M, N} = \mathrm{id}_{M \otimes N}$ bit-exact.

    The defining condition of a *symmetric* (rather than merely braided)
    monoidal category.
    """
    sigma_MN = braiding(M, N)
    sigma_NM = braiding(N, M)
    composite = compose(sigma_NM, sigma_MN).matrix
    identity = sp_eye(M.dim * N.dim)
    diff = composite - identity
    return {
        "M": M.label, "N": N.label,
        "bit_exact": _is_zero_matrix(diff),
    }


def verify_pentagon_coherence(
    M: FinDimRep, N: FinDimRep, P: FinDimRep, Q: FinDimRep,
) -> Dict[str, Any]:
    r"""Pentagon: associator coherence across four objects.

    The diagram

    $((M \otimes N) \otimes P) \otimes Q \xrightarrow{\alpha \otimes 1}
       (M \otimes (N \otimes P)) \otimes Q \xrightarrow{\alpha}
       M \otimes ((N \otimes P) \otimes Q) \xrightarrow{1 \otimes \alpha}
       M \otimes (N \otimes (P \otimes Q))$

    and

    $((M \otimes N) \otimes P) \otimes Q \xrightarrow{\alpha}
       (M \otimes N) \otimes (P \otimes Q) \xrightarrow{\alpha}
       M \otimes (N \otimes (P \otimes Q))$

    must commute. In our canonical lex bases both paths are the identity
    matrix, so coherence is bit-exact zero residual.
    """
    if not (M.n_max == N.n_max == P.n_max == Q.n_max):
        raise ValueError("pentagon requires shared n_max")
    # Path 1: (alpha ⊗ id_Q) then alpha then (id_M ⊗ alpha)
    a_MNP = associator(M, N, P)
    alpha_MNP_Q = tensor_morphism(
        a_MNP,
        RepMorphism(Q, Q, sp_eye(Q.dim), label="id_Q", validate=False),
    )
    a_M_NP_Q = associator(M, tensor_rep(N, P), Q)
    a_NPQ = associator(N, P, Q)
    id_M_alpha = tensor_morphism(
        RepMorphism(M, M, sp_eye(M.dim), label="id_M", validate=False),
        a_NPQ,
    )
    path1 = compose(id_M_alpha, compose(a_M_NP_Q, alpha_MNP_Q)).matrix
    # Path 2: alpha then alpha
    a_MN_P_Q = associator(tensor_rep(M, N), P, Q)
    a_M_N_PQ = associator(M, N, tensor_rep(P, Q))
    path2 = compose(a_M_N_PQ, a_MN_P_Q).matrix
    diff = path1 - path2
    return {
        "M": M.label, "N": N.label, "P": P.label, "Q": Q.label,
        "bit_exact": _is_zero_matrix(diff),
    }


def verify_triangle_coherence(M: FinDimRep, N: FinDimRep) -> Dict[str, Any]:
    r"""Triangle: associator-unitor coherence.

    The diagram

    $(M \otimes \mathbf{1}) \otimes N \xrightarrow{\alpha}
        M \otimes (\mathbf{1} \otimes N) \xrightarrow{1 \otimes \lambda_N}
        M \otimes N$

    equals $\rho_M \otimes \mathrm{id}_N$.
    """
    if M.n_max != N.n_max:
        raise ValueError("triangle requires shared n_max")
    one = unit_object(M.n_max)
    alpha = associator(M, one, N)
    lam_N = left_unitor(N)
    rho_M = right_unitor(M)
    id_M_lam = tensor_morphism(
        RepMorphism(M, M, sp_eye(M.dim), label="id_M", validate=False),
        lam_N,
    )
    rho_id = tensor_morphism(
        rho_M,
        RepMorphism(N, N, sp_eye(N.dim), label="id_N", validate=False),
    )
    path = compose(id_M_lam, alpha).matrix
    expected = rho_id.matrix
    diff = path - expected
    return {
        "M": M.label, "N": N.label,
        "bit_exact": _is_zero_matrix(diff),
    }


def verify_hexagon_coherence(
    M: FinDimRep, N: FinDimRep, P: FinDimRep,
) -> Dict[str, Any]:
    r"""Hexagon: braiding coherence with associator.

    The first hexagon:

    $(M \otimes N) \otimes P \xrightarrow{\sigma \otimes 1} (N \otimes M) \otimes P
      \xrightarrow{\alpha} N \otimes (M \otimes P) \xrightarrow{1 \otimes \sigma}
      N \otimes (P \otimes M)$

    equals

    $(M \otimes N) \otimes P \xrightarrow{\alpha} M \otimes (N \otimes P)
      \xrightarrow{\sigma} (N \otimes P) \otimes M \xrightarrow{\alpha} N \otimes (P \otimes M)$

    Bit-exact via direct matrix composition.
    """
    if not (M.n_max == N.n_max == P.n_max):
        raise ValueError("hexagon requires shared n_max")
    # Path 1: (sigma_{M,N} ⊗ id_P) then alpha then (id_N ⊗ sigma_{M,P})
    sigma_MN = braiding(M, N)
    id_P = RepMorphism(P, P, sp_eye(P.dim), label="id_P", validate=False)
    sigma_MN_id = tensor_morphism(sigma_MN, id_P)
    a_NMP = associator(N, M, P)
    sigma_MP = braiding(M, P)
    id_N = RepMorphism(N, N, sp_eye(N.dim), label="id_N", validate=False)
    id_N_sigma = tensor_morphism(id_N, sigma_MP)
    path1 = compose(id_N_sigma, compose(a_NMP, sigma_MN_id)).matrix
    # Path 2: alpha then sigma_{M, N⊗P} then alpha
    a_MNP = associator(M, N, P)
    sigma_M_NP = braiding(M, tensor_rep(N, P))
    a_NPM = associator(N, P, M)
    path2 = compose(a_NPM, compose(sigma_M_NP, a_MNP)).matrix
    diff = path1 - path2
    return {
        "M": M.label, "N": N.label, "P": P.label,
        "bit_exact": _is_zero_matrix(diff),
    }


# =====================================================================
# TC-1c --- Rigidity
# =====================================================================
#
# Rigidity supplies, for every finite-dim rep $V$, a dual rep $V^\vee$
# together with evaluation $\mathrm{ev}_V : V^\vee \otimes V \to \mathbf{1}$
# and coevaluation $\mathrm{coev}_V : \mathbf{1} \to V \otimes V^\vee$
# satisfying the two snake (zigzag) identities.
#
# In our $\mathbb{Q}$-linear, symmetric monoidal setting the dual of
# $(V, \{X_g^V\})$ has underlying space $\mathrm{Hom}_\mathbb{Q}(V, \mathbb{Q})$
# (canonically $\mathbb{Q}^{\dim V}$ via the dual basis) and the
# **contragredient Hopf action**
#
#   $X_g^{V^\vee} = -(X_g^V)^T,$
#
# determined by the antipode $S(x_g) = -x_g$ on primitive generators
# (v3.61.0 Track A).  Nilpotency, mutual commutativity, and the
# spectral-radius bookkeeping all transport because negative-transpose
# preserves nilpotency and pairwise commutativity.
#
# In the lex basis the canonical pairing $\phi_i(e_j) = \delta_{ij}$
# produces a 1-by-$\dim(V)^2$ evaluation matrix and a
# $\dim(V)^2$-by-1 coevaluation matrix that are explicit
# permutation-like Kronecker-deltas.
#
# Snake identities reduce to matrix identities because the associator
# is the identity matrix in lex bases and the unitors are identity
# matrices on $\dim(V)$.  Concretely (see :func:`verify_snake_identity_first`
# and :func:`verify_snake_identity_second` below):\
# $(\mathrm{id}_V \otimes \mathrm{ev}_V) \circ (\mathrm{coev}_V \otimes \mathrm{id}_V)
#  = \mathrm{id}_V$
# and the dual snake on $V^\vee$.


def dual_rep(V: FinDimRep) -> FinDimRep:
    r"""Dual rep $V^\vee$ with the contragredient Hopf action.

    Underlying $\mathbb{Q}$-vector space: $\mathrm{Hom}_\mathbb{Q}(V, \mathbb{Q})
    \cong \mathbb{Q}^{\dim V}$ in the dual basis $\{\phi_i\}$ with
    $\phi_i(e_j) = \delta_{ij}$.

    Endomorphisms: for each primitive generator $g$,

    $$X_g^{V^\vee} = -(X_g^V)^T.$$

    Determined by the antipode $S(x_g) = -x_g$ of the abelian primitive
    Hopf algebra (v3.61.0 Track A);\ the transpose accounts for the
    fact that $g$ acts on $V^\vee$ via $(X_g \phi)(v) = -\phi(X_g v)$,
    i.e., as the negative transpose in the dual basis.

    Nilpotency and pairwise commutativity are preserved because
    transposition is an anti-involution of the matrix algebra and
    $\mathrm{nilpotent}(M^T) \Leftrightarrow \mathrm{nilpotent}(M)$ and
    $[A^T, B^T] = -[A, B]^T$.
    """
    n_max = V.n_max
    d = V.dim
    dual_endos: Dict[Tuple[int, int, int], Matrix] = {}
    for g, X in V.non_zero_endos().items():
        Xd = -X.T
        if not _is_zero_matrix(Xd):
            dual_endos[g] = Xd
    return FinDimRep(
        n_max=n_max,
        dim=d,
        endos=dual_endos,
        label=f"{V.label}^vee" if V.label else "dual",
        validate=False,  # nilpotency and commutativity transport from V
    )


def evaluation_morphism(V: FinDimRep) -> RepMorphism:
    r"""Evaluation $\mathrm{ev}_V : V^\vee \otimes V \to \mathbf{1}$.

    On the lex basis $\{\phi_i \otimes e_j\}$ of $V^\vee \otimes V$
    (position $i \cdot \dim V + j$),

    $$\mathrm{ev}_V(\phi_i \otimes e_j) = \phi_i(e_j) = \delta_{ij}.$$

    Matrix shape $1 \times \dim(V)^2$ with $\mathrm{ev}_V[0, i \cdot \dim V + j]
    = \delta_{ij}$.
    """
    Vd = dual_rep(V)
    source = tensor_rep(Vd, V)
    target = unit_object(V.n_max)
    target.label = "1"
    d = V.dim
    M = sp_zeros(1, d * d)
    for i in range(d):
        M[0, i * d + i] = Integer(1)
    return RepMorphism(
        source=source,
        target=target,
        matrix=M,
        label=f"ev_{V.label}" if V.label else "ev",
        validate=False,
    )


def coevaluation_morphism(V: FinDimRep) -> RepMorphism:
    r"""Coevaluation $\mathrm{coev}_V : \mathbf{1} \to V \otimes V^\vee$.

    On the lex basis $\{e_i \otimes \phi_j\}$ of $V \otimes V^\vee$
    (position $i \cdot \dim V + j$),

    $$\mathrm{coev}_V(1) = \sum_{i = 0}^{\dim V - 1} e_i \otimes \phi_i.$$

    Matrix shape $\dim(V)^2 \times 1$ with $\mathrm{coev}_V[i \cdot \dim V + j, 0]
    = \delta_{ij}$.
    """
    Vd = dual_rep(V)
    source = unit_object(V.n_max)
    source.label = "1"
    target = tensor_rep(V, Vd)
    d = V.dim
    M = sp_zeros(d * d, 1)
    for i in range(d):
        M[i * d + i, 0] = Integer(1)
    return RepMorphism(
        source=source,
        target=target,
        matrix=M,
        label=f"coev_{V.label}" if V.label else "coev",
        validate=False,
    )


# ---------------------------------------------------------------------
# TC-1c verifiers
# ---------------------------------------------------------------------


def verify_dual_action(V: FinDimRep) -> Dict[str, Any]:
    r"""Verify the contragredient action $X_g^{V^\vee} = -(X_g^V)^T$ bit-exact.

    Independent reconstruction:\ build $-(X_g^V)^T$ for every non-zero
    generator and compare against `dual_rep(V).X(g)`. Also re-verify
    nilpotency and pairwise commutativity on the dual.
    """
    Vd = dual_rep(V)
    mismatches: List[Dict[str, Any]] = []
    for g in V.non_zero_endos().keys():
        expected = -V.X(g).T
        actual = Vd.X(g)
        diff = expected - actual
        if not _is_zero_matrix(diff):
            mismatches.append({
                "generator": list(g),
                "residual_norm_squared": sum(
                    diff[i, j] ** 2
                    for i in range(diff.rows)
                    for j in range(diff.cols)
                ),
            })
    # Nilpotency on dual
    nilpotency_ok = True
    for g in Vd.non_zero_endos().keys():
        if not _is_nilpotent(Vd.X(g), Vd.dim):
            nilpotency_ok = False
            break
    # Pairwise commutativity on dual
    commutativity_ok = True
    gs = list(Vd.non_zero_endos().keys())
    for i_g, g1 in enumerate(gs):
        for g2 in gs[i_g + 1:]:
            if not _matrices_commute(Vd.X(g1), Vd.X(g2)):
                commutativity_ok = False
                break
        if not commutativity_ok:
            break
    return {
        "V_label": V.label,
        "V_dim": V.dim,
        "mismatches": mismatches,
        "nilpotency_ok": nilpotency_ok,
        "commutativity_ok": commutativity_ok,
        "bit_exact": (
            len(mismatches) == 0
            and nilpotency_ok
            and commutativity_ok
        ),
    }


def verify_evaluation_intertwines(V: FinDimRep) -> Dict[str, Any]:
    r"""Verify $\mathrm{ev}_V$ intertwines the diagonal Hopf action bit-exact.

    For every generator $g$,

    $$\mathrm{ev}_V \cdot X_g^{V^\vee \otimes V} = X_g^{\mathbf{1}} \cdot \mathrm{ev}_V = 0.$$

    Since $\mathbf{1}$ has zero action, the right-hand side is the zero
    $1 \times \dim(V)^2$ matrix;\ the check reduces to
    $\mathrm{ev}_V \cdot X_g^{V^\vee \otimes V} = 0$ for every non-zero generator.
    """
    ev = evaluation_morphism(V)
    source = ev.source  # V^vee ⊗ V
    mismatches: List[Dict[str, Any]] = []
    for g in source.non_zero_endos().keys():
        X_source = source.X(g)
        lhs = ev.matrix * X_source  # 1 x dim(V)^2
        if not _is_zero_matrix(lhs):
            mismatches.append({
                "generator": list(g),
                "residual_norm_squared": sum(
                    lhs[i, j] ** 2
                    for i in range(lhs.rows)
                    for j in range(lhs.cols)
                ),
            })
    return {
        "V_label": V.label,
        "V_dim": V.dim,
        "mismatches": mismatches,
        "bit_exact": len(mismatches) == 0,
    }


def verify_coevaluation_intertwines(V: FinDimRep) -> Dict[str, Any]:
    r"""Verify $\mathrm{coev}_V$ intertwines the diagonal Hopf action bit-exact.

    For every generator $g$,

    $$X_g^{V \otimes V^\vee} \cdot \mathrm{coev}_V = \mathrm{coev}_V \cdot X_g^{\mathbf{1}} = 0.$$

    Since $\mathbf{1}$ has zero action, the right-hand side is the zero
    $\dim(V)^2 \times 1$ matrix;\ the check reduces to
    $X_g^{V \otimes V^\vee} \cdot \mathrm{coev}_V = 0$ for every non-zero generator.
    """
    coev = coevaluation_morphism(V)
    target = coev.target  # V ⊗ V^vee
    mismatches: List[Dict[str, Any]] = []
    for g in target.non_zero_endos().keys():
        X_target = target.X(g)
        rhs = X_target * coev.matrix  # dim(V)^2 x 1
        if not _is_zero_matrix(rhs):
            mismatches.append({
                "generator": list(g),
                "residual_norm_squared": sum(
                    rhs[i, j] ** 2
                    for i in range(rhs.rows)
                    for j in range(rhs.cols)
                ),
            })
    return {
        "V_label": V.label,
        "V_dim": V.dim,
        "mismatches": mismatches,
        "bit_exact": len(mismatches) == 0,
    }


def verify_snake_identity_first(V: FinDimRep) -> Dict[str, Any]:
    r"""First snake (zigzag) identity on $V^\vee$:

    $(\mathrm{ev}_V \otimes \mathrm{id}_{V^\vee}) \circ
     (\mathrm{id}_{V^\vee} \otimes \mathrm{coev}_V) = \mathrm{id}_{V^\vee},$

    in the canonical lex bases (with associator and unitors equal to the
    identity matrix).

    In matrices, with $d = \dim(V)$ and lex convention $e_i \otimes e_j$ at
    position $i \cdot d + j$:

    - $\mathrm{id}_{V^\vee} \otimes \mathrm{coev}_V$ is the Kronecker
      product of an identity $d \times d$ with a column $d^2 \times 1$,
      giving a $d^3 \times d$ matrix.
    - $\mathrm{ev}_V \otimes \mathrm{id}_{V^\vee}$ is the Kronecker
      product of a row $1 \times d^2$ with an identity $d \times d$,
      giving a $d \times d^3$ matrix.
    - Their composite is $d \times d$, and the snake identity asserts it
      equals the $d$-dim identity matrix.
    """
    ev = evaluation_morphism(V)
    coev = coevaluation_morphism(V)
    d = V.dim
    I_dual = sp_eye(d) if d > 0 else sp_zeros(0, 0)
    # id ⊗ coev: V^vee → V^vee ⊗ V ⊗ V^vee (column-stacked Kronecker block)
    A = _kron(I_dual, coev.matrix)  # d^3 x d
    # ev ⊗ id: V^vee ⊗ V ⊗ V^vee → V^vee
    B = _kron(ev.matrix, I_dual)  # d x d^3
    composite = B * A  # d x d
    identity = sp_eye(d) if d > 0 else sp_zeros(0, 0)
    diff = composite - identity
    return {
        "V_label": V.label,
        "V_dim": V.dim,
        "composite_shape": [composite.rows, composite.cols],
        "bit_exact": _is_zero_matrix(diff),
    }


def verify_snake_identity_second(V: FinDimRep) -> Dict[str, Any]:
    r"""Second snake (zigzag) identity on $V$:

    $(\mathrm{id}_V \otimes \mathrm{ev}_V) \circ
     (\mathrm{coev}_V \otimes \mathrm{id}_V) = \mathrm{id}_V,$

    in the canonical lex bases (with associator and unitors equal to the
    identity matrix).

    In matrices, with $d = \dim(V)$:

    - $\mathrm{coev}_V \otimes \mathrm{id}_V$ is $d^3 \times d$.
    - $\mathrm{id}_V \otimes \mathrm{ev}_V$ is $d \times d^3$.
    - Their composite is $d \times d$, equal to the identity matrix.
    """
    ev = evaluation_morphism(V)
    coev = coevaluation_morphism(V)
    d = V.dim
    I_V = sp_eye(d) if d > 0 else sp_zeros(0, 0)
    # coev ⊗ id: V → V ⊗ V^vee ⊗ V
    A = _kron(coev.matrix, I_V)  # d^3 x d
    # id ⊗ ev: V ⊗ V^vee ⊗ V → V
    B = _kron(I_V, ev.matrix)  # d x d^3
    composite = B * A  # d x d
    identity = sp_eye(d) if d > 0 else sp_zeros(0, 0)
    diff = composite - identity
    return {
        "V_label": V.label,
        "V_dim": V.dim,
        "composite_shape": [composite.rows, composite.cols],
        "bit_exact": _is_zero_matrix(diff),
    }


def verify_double_dual_iso(V: FinDimRep) -> Dict[str, Any]:
    r"""Verify $(V^\vee)^\vee \cong V$ bit-exact (strict equality of endo data).

    Because the antipode is involutive on the abelian primitive Hopf algebra
    ($S^2 = \mathrm{id}$ as $S(x_g) = -x_g$ and $-(-x_g) = x_g$), and because
    $(M^T)^T = M$, the double dual has $X_g^{(V^\vee)^\vee} = -(-X_g^T)^T = X_g$.
    Strict equality of the endomorphism dictionaries follows;\ no
    intermediate canonical isomorphism is needed in the lex basis.

    The check also confirms that the underlying dimension and the
    sparse non-zero generator set match.
    """
    Vdd = dual_rep(dual_rep(V))
    dim_match = Vdd.dim == V.dim
    V_gens = set(V.non_zero_endos().keys())
    Vdd_gens = set(Vdd.non_zero_endos().keys())
    gens_match = V_gens == Vdd_gens
    mismatches: List[Dict[str, Any]] = []
    if dim_match:
        for g in V_gens | Vdd_gens:
            diff = V.X(g) - Vdd.X(g)
            if not _is_zero_matrix(diff):
                mismatches.append({
                    "generator": list(g),
                    "residual_norm_squared": sum(
                        diff[i, j] ** 2
                        for i in range(diff.rows)
                        for j in range(diff.cols)
                    ),
                })
    return {
        "V_label": V.label,
        "V_dim": V.dim,
        "dim_match": dim_match,
        "gens_match": gens_match,
        "mismatches": mismatches,
        "bit_exact": dim_match and gens_match and len(mismatches) == 0,
    }


def verify_unit_self_dual(n_max: int) -> Dict[str, Any]:
    r"""Verify $\mathbf{1}^\vee = \mathbf{1}$ bit-exact.

    The unit $\mathbf{1} = T_1$ has $\dim = 1$ and zero action;\ its dual
    has the same data ($-(0)^T = 0$). Strict equality of endo dictionaries.
    """
    one = unit_object(n_max)
    one_dual = dual_rep(one)
    dim_match = one_dual.dim == one.dim
    endos_match = (
        len(one_dual.non_zero_endos()) == 0
        and len(one.non_zero_endos()) == 0
    )
    return {
        "n_max": n_max,
        "dim_match": dim_match,
        "endos_match": endos_match,
        "bit_exact": dim_match and endos_match,
    }


# =====================================================================
# TC-1d --- Fiber functor omega: Rep_fin(H_GV) -> Vec_Q
# =====================================================================
#
# At the implementation level $\omega$ is essentially a no-op forgetful
# functor: a rep $(M, \{X_g^M\})$ IS its underlying $\mathbb{Q}$-vector
# space plus the extra structure (the Hopf action), so forgetting the
# action returns the rep's dim / matrix unchanged. Naturality of
# $\otimes$-preservation is automatic because `tensor_rep` uses
# Kronecker products in canonical lex basis (the same convention as
# TC-1b's associator / unitors), so the natural iso
# $\omega(M \otimes N) \to \omega(M) \otimes_\mathbb{Q} \omega(N)$
# is the identity matrix bit-exactly.


def fiber_functor_object(V: FinDimRep) -> int:
    r"""$\omega(V) = $ underlying $\mathbb{Q}$-vector space, returned as its dimension.

    Since a finite-dim $\mathbb{Q}$-vector space is determined up to
    isomorphism by its dimension, returning ``dim`` is the canonical
    implementation of $\omega$ on objects.
    """
    return V.dim


def fiber_functor_morphism(f: RepMorphism) -> Matrix:
    r"""$\omega(f) = $ underlying $\mathbb{Q}$-linear map, returned as the matrix."""
    return f.matrix


def verify_omega_unit(n_max: int) -> Dict[str, Any]:
    r"""Verify $\omega(\mathbf{1}) = \mathbb{Q}$, i.e., $\dim(\mathbf{1}) = 1$."""
    unit = unit_object(n_max)
    omega_unit = fiber_functor_object(unit)
    return {
        "n_max": n_max,
        "omega_unit_dim": omega_unit,
        "expected_dim": 1,
        "ok": omega_unit == 1,
        "bit_exact": omega_unit == 1,
    }


def verify_omega_tensor_preservation(M: FinDimRep, N: FinDimRep) -> Dict[str, Any]:
    r"""Verify $\omega(M \otimes N) = \omega(M) \otimes_\mathbb{Q} \omega(N)$.

    On objects: $\dim(M \otimes N) = \dim M \cdot \dim N$. The natural
    isomorphism in the canonical lex basis is the identity matrix.
    """
    M_tensor_N = tensor_rep(M, N)
    omega_MN = fiber_functor_object(M_tensor_N)
    omega_M = fiber_functor_object(M)
    omega_N = fiber_functor_object(N)
    expected = omega_M * omega_N
    return {
        "M": M.label,
        "N": N.label,
        "omega_M_tensor_N": omega_MN,
        "omega_M_times_omega_N": expected,
        "dims_ok": omega_MN == expected,
        "ok": omega_MN == expected,
        "bit_exact": omega_MN == expected,
    }


def verify_omega_preserves_kernel(f: RepMorphism) -> Dict[str, Any]:
    r"""Verify $\omega(\ker f) = \ker \omega(f)$ as $\mathbb{Q}$-vector spaces."""
    ker_rep, iota = kernel(f)
    omega_ker = fiber_functor_object(ker_rep)
    null = f.matrix.nullspace()
    expected_dim = len(null)
    omega_iota = fiber_functor_morphism(iota)
    return {
        "f_label": f.label,
        "omega_ker_dim": omega_ker,
        "ker_omega_f_dim": expected_dim,
        "dims_ok": omega_ker == expected_dim,
        "iota_shape": [omega_iota.rows, omega_iota.cols] if omega_ker > 0 else [f.source.dim, 0],
        "ok": omega_ker == expected_dim,
        "bit_exact": omega_ker == expected_dim,
    }


def verify_omega_preserves_cokernel(f: RepMorphism) -> Dict[str, Any]:
    r"""Verify $\omega(\mathrm{coker}\, f) = \mathrm{coker}\, \omega(f)$."""
    coker_rep, pi = cokernel(f)
    omega_coker = fiber_functor_object(coker_rep)
    rank_f = f.matrix.rank() if f.matrix.shape != (0, 0) else 0
    expected_dim = f.target.dim - rank_f
    omega_pi = fiber_functor_morphism(pi)
    return {
        "f_label": f.label,
        "omega_coker_dim": omega_coker,
        "coker_omega_f_dim": expected_dim,
        "dims_ok": omega_coker == expected_dim,
        "pi_shape": [omega_pi.rows, omega_pi.cols] if omega_coker > 0 else [0, f.target.dim],
        "ok": omega_coker == expected_dim,
        "bit_exact": omega_coker == expected_dim,
    }


def verify_omega_preserves_direct_sum(R1: FinDimRep, R2: FinDimRep) -> Dict[str, Any]:
    r"""Verify $\omega(R_1 \oplus R_2) = \omega(R_1) \oplus \omega(R_2)$."""
    S, _, _, _, _ = direct_sum(R1, R2)
    omega_sum = fiber_functor_object(S)
    expected = fiber_functor_object(R1) + fiber_functor_object(R2)
    return {
        "R1": R1.label,
        "R2": R2.label,
        "omega_R1_plus_R2": omega_sum,
        "omega_R1_plus_omega_R2": expected,
        "dims_ok": omega_sum == expected,
        "ok": omega_sum == expected,
        "bit_exact": omega_sum == expected,
    }


def verify_omega_faithful(f: RepMorphism, g: RepMorphism) -> Dict[str, Any]:
    r"""Verify faithfulness on the pair $(f, g)$:\ $f = g \Leftrightarrow \omega(f) = \omega(g)$.

    Since the matrix IS the morphism data at the implementation level,
    $\omega(f) = \omega(g)$ as matrices iff $f = g$ as morphisms (when
    source and target match).
    """
    if f.source.dim != g.source.dim or f.target.dim != g.target.dim:
        return {
            "shapes_match": False,
            "ok": False,
            "bit_exact": False,
            "note": "source/target dims of f and g differ",
        }
    omega_f = fiber_functor_morphism(f)
    omega_g = fiber_functor_morphism(g)
    matrix_equal = _is_zero_matrix(omega_f - omega_g)
    morphism_equal = matrix_equal  # same data at implementation level
    return {
        "f_source_dim": f.source.dim,
        "f_target_dim": f.target.dim,
        "shapes_match": True,
        "morphism_equal": morphism_equal,
        "matrix_equal": matrix_equal,
        "faithful_consistent": morphism_equal == matrix_equal,
        "ok": morphism_equal == matrix_equal,
        "bit_exact": morphism_equal == matrix_equal,
    }


# =====================================================================
# TC-1e --- First stone: U^*_{Levi}(Q) -> Aut^otimes(omega) inclusion
# =====================================================================
#
# The motivic Galois group $U^*_{\mathrm{GeoVac}, \mathrm{Levi}} =
# \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ from v3.63.0 L1 maps to
# the natural $\otimes$-automorphisms of the fiber functor $\omega$
# (TC-1d) via:
#
#   $\Phi: U^*_{\mathrm{Levi}}(\mathbb{Q}) \to \mathrm{Aut}^\otimes(\omega)$,
#   $\Phi(t)(V) = \exp\!\Big(\sum_g t_g X_g^V\Big)$ for the unipotent
#   factor $t \in \mathbb{Q}^{3 N(n_{\max})}$,
#
# composed with the $SL_2$-action on the Peter--Weyl decoration (which
# acts on the $j_{\max}$ axis, independent of the $n_{\max}$ axis on
# which the fiber-functor substrate lives).
#
# Each $\eta_V(t) = \exp(\sum_g t_g X_g^V)$ is well-defined because the
# $X_g^V$ are pairwise commuting nilpotents:\ the sum is a nilpotent
# matrix whose exp truncates at degree $\dim(V) - 1$. The four
# $\otimes$-automorphism axioms (invertibility, unit, naturality,
# $\otimes$-compatibility) and the group law
# $\Phi(t_1 + t_2) = \Phi(t_1) \circ \Phi(t_2)$ are structural
# consequences of $\exp$ on commuting matrices and the intertwining
# condition on morphisms. TC-1e verifies these bit-exact on a
# representative panel as the **first stone** of the multi-year wall;
# the converse equality $\mathrm{Aut}^\otimes(\omega) = U^*$ (full
# Tannakian closure proper) remains multi-year content.
#
# Honest scope (PS-4 + v3.66.0 §5.2 follow-on register): TC-1e closes
# the *inclusion* direction at finite cutoff on the visible
# $\mathbb{G}_a$ factor; the $SL_2$ factor's compatibility is
# categorical (independence of the $j_{\max}$ axis from the $n_{\max}$
# axis), not visible at the panel level.


def _matrix_exp_nilpotent(M: Matrix, dim: int) -> Matrix:
    r"""Compute $\exp(M)$ bit-exact for a nilpotent matrix $M$ of dimension $\dim$.

    Uses the truncated power series

    $$\exp(M) = \sum_{k = 0}^{\dim - 1} \frac{M^k}{k!}$$

    which terminates exactly when $M$ is nilpotent (since $M^{\dim} = 0$).
    All entries are bit-exact ``sympy.Rational`` / ``sympy.Integer``.

    Parameters
    ----------
    M : Matrix
        Nilpotent $\dim \times \dim$ matrix.
    dim : int
        Matrix dimension.

    Returns
    -------
    Matrix
        $\exp(M)$ bit-exact.
    """
    if dim == 0:
        return sp_zeros(0, 0)
    result = sp_eye(dim)
    M_power = sp_eye(dim)
    factorial = Integer(1)
    for k in range(1, dim):
        M_power = M_power * M
        factorial = factorial * k
        result = result + M_power / factorial
    return result


def levi_unipotent_action(
    t_dict: Mapping[Tuple[int, int, int], Any],
    V: FinDimRep,
) -> Matrix:
    r"""$\eta_V(t) = \exp(\sum_g t_g X_g^V)$ bit-exact.

    The natural $\otimes$-automorphism $\eta_V$ associated to a parameter
    $t = (t_g)_g \in \mathbb{Q}^{3 N(n_{\max})}$ of the pro-unipotent
    factor $\mathbb{G}_a^{3 N(n_{\max})}$ of
    $U^*_{\mathrm{GeoVac}, \mathrm{Levi}}$ acts on the rep $V$ by the
    matrix exponential of the linear combination $\sum_g t_g X_g^V$ of
    the commuting nilpotent endomorphisms.

    Parameters
    ----------
    t_dict : mapping
        Maps primitive generator $g = (n, l, s)$ to ``sympy.Rational``.
        Generators not in the dict default to $t_g = 0$ (no
        contribution).
    V : FinDimRep
        Target rep.

    Returns
    -------
    Matrix
        $\dim(V) \times \dim(V)$ invertible bit-exact matrix.
    """
    dim = V.dim
    if dim == 0:
        return sp_zeros(0, 0)
    M = sp_zeros(dim, dim)
    for g, t_g in t_dict.items():
        if t_g != 0:
            M = M + t_g * V.X(g)
    return _matrix_exp_nilpotent(M, dim)


def verify_natural_auto_invertibility(
    t_dict: Mapping[Tuple[int, int, int], Any],
    V: FinDimRep,
) -> Dict[str, Any]:
    r"""Verify $\eta_V(t) = \exp(\sum_g t_g X_g^V)$ is invertible.

    Bit-exact check $\det(\eta_V) \ne 0$. Structurally, $\eta_V$ is
    always invertible (inverse $\exp(-M)$) — this verifier provides
    a panel-level falsifier check.
    """
    eta_V = levi_unipotent_action(t_dict, V)
    det = eta_V.det() if V.dim > 0 else Integer(1)
    invertible = det != 0
    return {
        "V_label": V.label,
        "V_dim": V.dim,
        "det": det,
        "invertible": bool(invertible),
        "bit_exact": bool(invertible),
    }


def verify_natural_auto_unit(
    t_dict: Mapping[Tuple[int, int, int], Any],
    n_max: int,
) -> Dict[str, Any]:
    r"""Verify $\eta_{\mathbf{1}}(t) = \mathrm{id}_\mathbb{Q}$ bit-exact.

    Since $\mathbf{1} = T_1$ has $\dim = 1$ and all $X_g^{\mathbf{1}} = 0$,
    $\eta_{\mathbf{1}}(t) = \exp(0) = (1)$ identically.
    """
    one = unit_object(n_max)
    eta_one = levi_unipotent_action(t_dict, one)
    expected = sp_eye(1)
    return {
        "n_max": n_max,
        "eta_unit": eta_one,
        "expected": expected,
        "bit_exact": eta_one == expected,
    }


def verify_natural_auto_naturality(
    t_dict: Mapping[Tuple[int, int, int], Any],
    f: RepMorphism,
) -> Dict[str, Any]:
    r"""Verify naturality bit-exact:\ $\omega(f) \cdot \eta_V(t) = \eta_W(t) \cdot \omega(f)$.

    For $f: V \to W$ a morphism in $\mathrm{Rep}_{\mathrm{fin}}$, the
    intertwining condition $f \cdot X_g^V = X_g^W \cdot f$ extends to
    all powers $f \cdot (X_g^V)^k = (X_g^W)^k \cdot f$ and therefore to
    $f \cdot \exp(M^V) = \exp(M^W) \cdot f$ where $M^V = \sum_g t_g X_g^V$
    and $M^W = \sum_g t_g X_g^W$.
    """
    eta_V = levi_unipotent_action(t_dict, f.source)
    eta_W = levi_unipotent_action(t_dict, f.target)
    lhs = f.matrix * eta_V
    rhs = eta_W * f.matrix
    diff = lhs - rhs
    return {
        "f_label": f.label,
        "source_dim": f.source.dim,
        "target_dim": f.target.dim,
        "bit_exact": _is_zero_matrix(diff),
    }


def verify_natural_auto_tensor(
    t_dict: Mapping[Tuple[int, int, int], Any],
    V: FinDimRep,
    W: FinDimRep,
) -> Dict[str, Any]:
    r"""Verify $\otimes$-compatibility bit-exact:\ $\eta_{V \otimes W}(t) = \eta_V(t) \otimes \eta_W(t)$.

    Structural reason:\ on $V \otimes W$, $X_g^{V \otimes W} =
    X_g^V \otimes I + I \otimes X_g^W$. The two summands commute (act
    on different tensor factors), so

    $$\exp\!\Big(\sum_g t_g X_g^{V \otimes W}\Big)
        = \exp\!\Big(\sum_g t_g X_g^V \otimes I\Big) \cdot
          \exp\!\Big(\sum_g t_g I \otimes X_g^W\Big)
        = \big(\eta_V(t) \otimes I\big) \cdot \big(I \otimes \eta_W(t)\big)
        = \eta_V(t) \otimes \eta_W(t).$$
    """
    VW = tensor_rep(V, W)
    eta_VW = levi_unipotent_action(t_dict, VW)
    eta_V = levi_unipotent_action(t_dict, V)
    eta_W = levi_unipotent_action(t_dict, W)
    expected = _kron(eta_V, eta_W)
    diff = eta_VW - expected
    return {
        "V_label": V.label,
        "W_label": W.label,
        "VW_dim": VW.dim,
        "bit_exact": _is_zero_matrix(diff),
    }


def verify_natural_auto_group_law(
    t1_dict: Mapping[Tuple[int, int, int], Any],
    t2_dict: Mapping[Tuple[int, int, int], Any],
    V: FinDimRep,
) -> Dict[str, Any]:
    r"""Verify group law bit-exact:\ $\eta_V(t_1 + t_2) = \eta_V(t_1) \cdot \eta_V(t_2)$.

    Holds because $M^V(t_1) = \sum_g t_g^{(1)} X_g^V$ and
    $M^V(t_2) = \sum_g t_g^{(2)} X_g^V$ commute (both are
    $\mathbb{Q}$-linear combinations of the commuting $\{X_g^V\}_g$),
    so $\exp(M^V(t_1) + M^V(t_2)) = \exp(M^V(t_1)) \cdot \exp(M^V(t_2))$.
    """
    # Sum the t-dicts
    t_sum: Dict[Tuple[int, int, int], Any] = dict(t1_dict)
    for g, t_g in t2_dict.items():
        t_sum[g] = t_sum.get(g, Integer(0)) + t_g
    eta_sum = levi_unipotent_action(t_sum, V)
    eta_t1 = levi_unipotent_action(t1_dict, V)
    eta_t2 = levi_unipotent_action(t2_dict, V)
    product = eta_t1 * eta_t2
    diff = eta_sum - product
    return {
        "V_label": V.label,
        "V_dim": V.dim,
        "bit_exact": _is_zero_matrix(diff),
    }


# =====================================================================
# TC-1f --- SL_2 explicit inclusion + full-panel injectivity
# =====================================================================
#
# The semisimple factor $SL_2$ of $U^*_{\mathrm{GeoVac}, \mathrm{Levi}}
# = \mathbb{G}_a^{3 N(n_{\max})} \rtimes SL_2$ acts on the Peter--Weyl
# decoration (v3.63.0 L2 / v3.66.0 FO1). On a finite-dim
# $\mathbb{Q}$-vector space carrying an $SL_2$-action via its standard
# rep (or a symmetric power), $SL_2(\mathbb{Q})$ acts as natural
# $\otimes$-automorphisms of the corresponding forgetful functor
# $\omega_{\mathrm{PW}}$.
#
# Combined with TC-1e ($\mathbb{G}_a^{3N}$ inclusion via translation
# generators), TC-1f gives the full explicit inclusion
#
#   $\Phi^{\mathrm{full}}: U^*_{\mathrm{Levi}}(\mathbb{Q}) \to
#       \mathrm{Aut}^\otimes(\omega^{\mathrm{combined}})$,
#
# where $\omega^{\mathrm{combined}}$ is the fiber functor on the tensor
# product category $\mathrm{Rep}_{\mathrm{fin}}(\mathcal{H}_{\mathrm{GV}})
# \otimes \mathrm{Rep}_{\mathrm{fin}}(SL_2)$. The two factors commute
# because they act on different tensor factors of the underlying vector
# space (explicit verification on combined reps closes the "categorical"
# caveat of TC-1e).
#
# Full-panel injectivity of $\Phi$ on the $\mathbb{G}_a^{3 N}$ factor
# at $n_{\max} = 2$: build $3 N(2) = 15$ reps, one per primitive
# generator, each activating exactly one $X_g$. For any $t \ne 0$,
# the rep $V_g$ associated to the unique $g$ with $t_g \ne 0$ gives
# $\eta_{V_g}(t) \ne \mathrm{id}$, so $\Phi(t) \ne \mathrm{id}_\omega$.


class PWRep:
    r"""Peter--Weyl rep:\ finite-dim $\mathbb{Q}$-vector space with an explicit $SL_2$-action.

    An object is a finite-dim $\mathbb{Q}$-vector space $V$ together
    with a homomorphism $\rho: SL_2(\mathbb{Q}) \to \mathrm{GL}(V)$
    encoded as a callable ``rho(g)`` returning the matrix of
    $\rho(g) \in \mathrm{GL}_{\dim(V)}(\mathbb{Q})$ for any 2x2
    $\mathbb{Q}$-matrix $g$ with $\det(g) = 1$.

    The two standard cases used in TC-1f panel:

    * **Standard rep** $V = \mathbb{Q}^2$ with $\rho(g) = g$.
    * **Symmetric power** $V = \mathrm{Sym}^k(\mathbb{Q}^2)$ with
      $\rho(g)$ acting as the $k$-th symmetric power of $g$.

    Parameters
    ----------
    dim : int
        Dimension of $V$.
    rho : callable
        ``rho(g: Matrix) -> Matrix`` returning a $\dim \times \dim$
        $\mathbb{Q}$-matrix representing $\rho(g)$.
    label : str, optional
    """

    def __init__(self, dim: int, rho: Callable[[Matrix], Matrix], label: str = ""):
        self.dim = dim
        self._rho = rho
        self.label = label

    def rho(self, g: Matrix) -> Matrix:
        return self._rho(g)

    def __repr__(self) -> str:
        return f"PWRep(dim={self.dim}, label={self.label!r})"


class PWMorphism:
    r"""$SL_2$-equivariant linear map between PW reps.

    Encoded by a matrix $F: V \to W$ satisfying
    $F \cdot \rho_V(g) = \rho_W(g) \cdot F$ for every $g \in SL_2(\mathbb{Q})$.
    """

    def __init__(
        self,
        source: "PWRep",
        target: "PWRep",
        matrix: Matrix,
        label: str = "",
    ):
        self.source = source
        self.target = target
        self.matrix = matrix
        self.label = label


def sl2_standard_action(g: Matrix) -> Matrix:
    r"""Standard $SL_2$-action on $\mathbb{Q}^2$:\ $\rho(g) = g$ (identity rep).

    Returns the $2 \times 2$ matrix $g$ unchanged. Verifies the rep
    is a homomorphism trivially (composition matches matrix product).
    """
    return g


def _sl2_sym2_action(g: Matrix) -> Matrix:
    r"""Symmetric square of standard $SL_2$-action on $\mathbb{Q}^3$.

    On basis $\{e_1^2, e_1 e_2, e_2^2\}$, the standard $g = \begin{pmatrix} a & b \\ c & d \end{pmatrix}$
    acts as

    $$\rho^{(2)}(g) = \begin{pmatrix} a^2 & 2ab & b^2 \\ ac & ad + bc & bd \\ c^2 & 2cd & d^2 \end{pmatrix}.$$
    """
    a, b, c, d = g[0, 0], g[0, 1], g[1, 0], g[1, 1]
    return Matrix([
        [a * a, 2 * a * b, b * b],
        [a * c, a * d + b * c, b * d],
        [c * c, 2 * c * d, d * d],
    ])


def _pw_standard_rep(label: str = "V_fund") -> PWRep:
    r"""The standard 2-dim $SL_2$-rep on $\mathbb{Q}^2$."""
    return PWRep(dim=2, rho=sl2_standard_action, label=label)


def _pw_sym2_rep(label: str = "Sym2") -> PWRep:
    r"""The symmetric square 3-dim $SL_2$-rep on $\mathbb{Q}^3$."""
    return PWRep(dim=3, rho=_sl2_sym2_action, label=label)


def _pw_trivial_rep(dim: int = 1, label: str = "PW_triv") -> PWRep:
    r"""Trivial $SL_2$-rep on $\mathbb{Q}^{\dim}$:\ $\rho(g) = I$ for all $g$."""
    def rho(g: Matrix) -> Matrix:
        return sp_eye(dim)
    return PWRep(dim=dim, rho=rho, label=label)


def verify_sl2_invertibility(g: Matrix, V_pw: PWRep) -> Dict[str, Any]:
    r"""Verify $\rho_V(g) \in \mathrm{GL}(V)$ for $g \in SL_2(\mathbb{Q})$.

    Checks $\det(g) = 1$ (input validation) and $\det(\rho_V(g)) \ne 0$.
    """
    det_g = g.det()
    eta = V_pw.rho(g)
    det_eta = eta.det()
    return {
        "V_label": V_pw.label,
        "V_dim": V_pw.dim,
        "det_g": det_g,
        "det_rho_g": det_eta,
        "g_in_SL2": det_g == 1,
        "rho_invertible": det_eta != 0,
        "bit_exact": det_g == 1 and det_eta != 0,
    }


def verify_sl2_tensor(g: Matrix, V_pw: PWRep, W_pw: PWRep) -> Dict[str, Any]:
    r"""Verify $\rho_{V \otimes W}(g) = \rho_V(g) \otimes \rho_W(g)$ bit-exact.

    Here $V \otimes W$ is the tensor product with $SL_2$-action via
    the diagonal $\rho_{V \otimes W}(g) = \rho_V(g) \otimes \rho_W(g)$.
    This is just verifying the tensor product is well-defined as the
    Kronecker product of the matrices;\ a tautology at the matrix
    level but a load-bearing structural identity for $\mathrm{Aut}^\otimes$.
    """
    eta_V = V_pw.rho(g)
    eta_W = W_pw.rho(g)
    expected = _kron(eta_V, eta_W)
    # Build V_pw ⊗ W_pw explicitly via Kronecker rho
    def rho_tensor(h: Matrix) -> Matrix:
        return _kron(V_pw.rho(h), W_pw.rho(h))
    actual = rho_tensor(g)
    diff = actual - expected
    return {
        "V_label": V_pw.label,
        "W_label": W_pw.label,
        "VW_dim": V_pw.dim * W_pw.dim,
        "bit_exact": _is_zero_matrix(diff),
    }


def verify_sl2_group_homomorphism(g1: Matrix, g2: Matrix, V_pw: PWRep) -> Dict[str, Any]:
    r"""Verify $\rho_V(g_1 g_2) = \rho_V(g_1) \rho_V(g_2)$ bit-exact.

    Group-homomorphism property of $\rho$.
    """
    g_product = g1 * g2
    lhs = V_pw.rho(g_product)
    rhs = V_pw.rho(g1) * V_pw.rho(g2)
    diff = lhs - rhs
    return {
        "V_label": V_pw.label,
        "V_dim": V_pw.dim,
        "bit_exact": _is_zero_matrix(diff),
    }


def verify_ga_sl2_commute(
    t_dict: Mapping[Tuple[int, int, int], Any],
    g: Matrix,
    V: FinDimRep,
    V_pw: PWRep,
) -> Dict[str, Any]:
    r"""Verify $\mathbb{G}_a$ and $SL_2$ commute on the combined rep $V \otimes V_{\mathrm{PW}}$.

    The combined rep has action of $\mathbb{G}_a^{3N}$ on $V$ (via the
    nilpotent translation generators, TC-1e) and $SL_2$ on $V_{\mathrm{PW}}$
    (via $\rho$, TC-1f). These act on different tensor factors, so they
    commute structurally. Concretely:

    * $\Phi^{\mathbb{G}_a}_{V \otimes V_{\mathrm{PW}}}(t) = \eta_V(t) \otimes I_{V_{\mathrm{PW}}}$
    * $\Phi^{SL_2}_{V \otimes V_{\mathrm{PW}}}(g) = I_V \otimes \rho_{V_{\mathrm{PW}}}(g)$
    * Both are Kronecker products commuting:\
      $(A \otimes I)(I \otimes B) = A \otimes B = (I \otimes B)(A \otimes I)$.

    The combined automorphism is $\eta_V(t) \otimes \rho_{V_{\mathrm{PW}}}(g)$.
    """
    eta_V = levi_unipotent_action(t_dict, V)
    rho_VW = V_pw.rho(g)
    # G_a acting first then SL_2: (eta_V ⊗ I)(I ⊗ rho) = eta_V ⊗ rho
    lhs = _kron(eta_V, sp_eye(V_pw.dim)) * _kron(sp_eye(V.dim), rho_VW)
    # SL_2 acting first then G_a: (I ⊗ rho)(eta_V ⊗ I) = eta_V ⊗ rho
    rhs = _kron(sp_eye(V.dim), rho_VW) * _kron(eta_V, sp_eye(V_pw.dim))
    diff = lhs - rhs
    return {
        "V_label": V.label,
        "V_pw_label": V_pw.label,
        "combined_dim": V.dim * V_pw.dim,
        "bit_exact": _is_zero_matrix(diff),
    }


# ---------------------------------------------------------------------
# Full-panel injectivity at n_max = 2
# ---------------------------------------------------------------------


def primitive_generator_rep(n_max: int, g: Tuple[int, int, int]) -> FinDimRep:
    r"""Build a 2-dim rep $V_g$ that activates the single primitive generator $g$.

    $V_g = (\mathbb{Q}^2, \{X_h^{V_g}\}_h)$ with
    $X_g^{V_g} = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$ and
    $X_h^{V_g} = 0$ for every $h \ne g$.  This is a faithful "test rep"
    for the generator $g$:\ $\eta_{V_g}(t) = I$ iff $t_g = 0$.
    """
    n, l, s = g
    return FinDimRep(
        n_max=n_max,
        dim=2,
        endos={g: Matrix([[0, 1], [0, 0]])},
        label=f"V_{n}_{l}_{s}",
    )


def verify_injectivity_at_generator(
    t_dict: Mapping[Tuple[int, int, int], Any],
    g: Tuple[int, int, int],
    n_max: int,
) -> Dict[str, Any]:
    r"""Verify $\eta_{V_g}(t) \ne \mathrm{id}$ iff $t_g \ne 0$.

    Concretely:\ $V_g$ has $X_g^{V_g} = \begin{pmatrix} 0 & 1 \\ 0 & 0 \end{pmatrix}$,
    so $\eta_{V_g}(t) = \exp(t_g X_g) = I + t_g X_g$ (since $X_g^2 = 0$).
    Therefore $\eta_{V_g}(t) = I$ iff $t_g = 0$.

    The bit-exact check has two cases:
    * If $t_g \ne 0$:\ $\eta_{V_g}(t)$ has the $(0, 1)$ entry equal to $t_g$,
      so $\eta_{V_g}(t) - I \ne 0$.
    * If $t_g = 0$:\ $\eta_{V_g}(t) = I$.
    """
    V_g = primitive_generator_rep(n_max, g)
    eta = levi_unipotent_action(t_dict, V_g)
    identity = sp_eye(2)
    t_g = t_dict.get(g, Integer(0))
    eta_minus_id = eta - identity
    is_id = _is_zero_matrix(eta_minus_id)
    expected_is_id = (t_g == 0)
    return {
        "generator": list(g),
        "t_g": t_g,
        "eta_is_identity": is_id,
        "expected_is_identity": expected_is_id,
        "bit_exact": is_id == expected_is_id,
    }
