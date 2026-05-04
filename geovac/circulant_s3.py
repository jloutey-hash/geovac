"""Circulant-style S^3 truncation: WH1 round-3 falsification comparator.

This module implements an alternative finite-dimensional truncation of S^3 that
plays the role on S^3 that the *circulant matrix algebra* C(C_m) plays on S^1
in Connes & van Suijlekom (CMP 2021, arXiv:2004.14115; see Outlook §6,
lines 2611-2618 of `debug/connes_vs_2004_14115.txt`):

  > "The examples we have considered show that spectral truncations allows one
  > to work with finite-dimensional operator systems, while keeping in tact the
  > full symmetry of the original space.  For instance, the Toeplitz operator
  > systems possess an S^1 symmetry, and as a consequence have a very rich
  > extremal and pure state space structure.  This is in contrast with the
  > circulant matrices, where the symmetry is reduced to a discrete group.
  > So, even though both spaces converge in Gromov-Hausdorff distance to the
  > circle, for the second one loses a lot of structure in the
  > finite-dimensional reduction."

WH1 round-2 (`debug/wh1_round2_propagation_memo.md`) showed that the Connes-vS
truncated operator system O_{n_max} = P_{n_max} C^infty(S^3) P_{n_max} on the
Fock-projected S^3 has propagation number prop(O_{n_max}) = 2 at every tested
n_max in {2, 3, 4}, *exactly* matching Connes-vS Proposition 4.2 for the
Toeplitz S^1 case verbatim.

WH1 round-3 (this module + the falsification memo) tests whether prop = 2 is a
STRUCTURAL property of the Connes-vS construction (in which case alternative
truncations should give different prop values), or a generic property of any
reasonable finite truncation of S^3 (in which case the round-2 finding is
weak).

The natural comparator suggested by Connes-vS Outlook §6 is the *circulant*
analog: a *commutative* C*-subalgebra of B(H_{n_max}) of the same dimension as
either the GeoVac envelope (N) or the GeoVac operator system (dim(O)).

Construction
============

We build the algebra of *diagonal matrices* in some chosen orthonormal basis
of C^N. Concretely:

    A_circ_N := { D in M_N(C) : D is diagonal in the canonical basis }
              = span{ E_{i,i} : i = 0, ..., N-1 }
              isomorphic to C^N (as a *-algebra)

This is a *commutative C*-subalgebra* of M_N(C), of complex dimension N. It is
the algebra of multiplication operators by functions on a finite N-point
"discrete S^3 sample." (Whether the points actually live on S^3 is irrelevant
for the propagation-number invariant; what matters is the multiplicative
closure of the algebra, which is automatic by diagonality.)

Crucially:

  - A_circ_N is closed under the matrix product (D_1 * D_2 is diagonal).
  - A_circ_N is closed under the adjoint (D^* is diagonal).
  - 1 in A_circ_N.

So A_circ_N is a *commutative finite-dimensional C*-algebra*. By Definition
2.39 of Connes-vS, since A_circ_N = (A_circ_N)^1 is itself a C*-algebra, we
have

    prop(A_circ_N) = 1   (trivially, by the definition).

This contrasts maximally with the GeoVac result prop(O_{n_max}) = 2 from
round 2.

Why this is the right comparator
================================

Connes-vS in Outlook §6 explicitly identify circulant matrices C(C_m) as the
*non-spectral-truncation* alternative to the Toeplitz operator system C(S^1)^{(n)}
on S^1:

  - C(S^1)^{(n)}: operator system, dim n^2 in M_n(C), prop = 2, S^1-symmetric.
  - C(C_m):       C*-algebra,    dim m  in M_m(C), prop = 1, only Z/m symmetric.

Both of these objects sit in M_N(C) for some N, and both are reasonable finite
discretizations of S^1, but they have *categorically different* propagation
invariants because the operator-system construction *preserves more of the
S^1 symmetry* than the circulant construction does.

Our A_circ_N on S^3 is the analogous "lose-the-symmetry" comparator: we keep
only diagonality (= Z/N-style cyclic symmetry, or even less if we don't bother
to put the diagonal in a Z/N-cyclic order) and discard the SO(4) selection
rules of the Connes-vS construction. The result, by C*-algebra triviality, is
prop = 1, consistent with the C(C_m) story.

Therefore: confirming numerically (and conceptually) that prop(A_circ_N) = 1
contrasts with prop(O_{n_max}) = 2, and demonstrates that the round-2 prop = 2
result is NOT a generic feature of "any finite truncation of S^3" but rather a
structurally-specific feature of the Connes-vS spectral truncation that
preserves the SO(4) symmetry.

Implementation note
===================

Most of the work in this module is therefore *bookkeeping*: we build A_circ_N
with the same propagation-number algorithm used in
`geovac/operator_system.py` (vec-stack rank), so that the comparison is
apples-to-apples (same numerical methodology, same tolerances).  We expect and
verify prop(A_circ_N) = 1 at every N >= 1.

If the numerical computation gives prop(A_circ_N) > 1, that means the
construction is buggy (either A_circ_N is not actually multiplicatively closed,
or the rank computation is broken). We provide an explicit
`verify_multiplicative_closure` method to debug this case.

References
==========

A. Connes & W. D. van Suijlekom, "Spectral Truncations in Noncommutative
Geometry and Operator Systems," CMP 383 (2021), arXiv:2004.14115. See
Definition 2.39 (propagation number), Proposition 4.2 (Toeplitz S^1 prop = 2
independent of n), Section 5 (Toeplitz vs circulant duality), and Outlook §6
(the explicit Toeplitz vs circulant comparison).

WH1 round-2 memo: debug/wh1_round2_propagation_memo.md (2026-05-03).

WH1 round-3 (this round) memo: debug/wh1_round3_falsification_memo.md.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np

# Import only for cross-comparison helper; do NOT modify the imported module.
from geovac.operator_system import (
    TruncatedOperatorSystem,
    hilbert_dim,
    operator_system_dim,
    propagation_number,
)


# ---------------------------------------------------------------------------
# CirculantS3Truncation: the commutative C*-subalgebra of diagonal matrices
# ---------------------------------------------------------------------------


@dataclass
class CirculantPropagationResult:
    """Bundle of values returned by CirculantS3Truncation.compute_propagation_number.

    prop : the propagation number (= 1 for any C*-algebra; the algorithm
        should always recover this).
    dim_sequence : list of dim(A^k) values for k = 1, 2, ... up to convergence.
    """

    prop: int
    dim_sequence: List[int]


class CirculantS3Truncation:
    """The commutative C*-subalgebra A_circ_N of diagonal matrices in M_N(C).

    This is the WH1 round-3 falsification comparator for the Connes-vS
    truncated operator system.  See module docstring for context.

    Attributes
    ----------
    n_points : int
        N, the dimension of the underlying "discrete S^3 sample." Equivalently,
        the dimension of A_circ as a complex *-algebra and the size of the
        ambient matrix algebra M_N(C) it sits inside.
    generators : list of complex numpy arrays
        The N elementary diagonal matrices E_{i, i}, each of shape (N, N).
    dim : int
        N, the dimension of A_circ as a complex linear subspace of M_N(C).
    envelope_dim : int
        N^2, the dimension of M_N(C). The C*-envelope of A_circ_N is M_N(C)
        only if we *embed* A_circ_N as a (proper) subspace; the *internal*
        C*-envelope of A_circ_N (as an abstract C*-algebra) is just
        A_circ_N itself.

        For the propagation-number test, what matters is whether
        (A_circ_N)^1 is *itself* a C*-algebra. It is (trivially, by closure
        under matrix multiplication of diagonal matrices), so prop = 1
        regardless of what the ambient C*-envelope is taken to be.

        We report envelope_dim = N (the *intrinsic* C*-envelope dimension of
        the algebra) by default, with an option to also compute "envelope =
        N^2" propagation following the convention used for
        TruncatedOperatorSystem (treat the ambient M_N(C) as the envelope and
        check whether powers of A_circ_N fill it).  In the latter case
        prop(A_circ_N -> M_N(C)) = infinity for N > 1, because A_circ_N is
        commutative and M_N(C) for N > 1 is non-commutative, so no power of
        A_circ_N can ever equal M_N(C).
    """

    def __init__(self, n_points: int, *, rank_tol: float = 1e-12) -> None:
        if n_points < 1:
            raise ValueError(f"n_points must be >= 1, got {n_points}")
        self.n_points = n_points
        self._rank_tol = rank_tol

        # The N elementary diagonal matrices E_{i,i}.
        self.generators: List[np.ndarray] = []
        for i in range(n_points):
            E = np.zeros((n_points, n_points), dtype=np.complex128)
            E[i, i] = 1.0
            self.generators.append(E)

        # Cache vec-stack and dim.
        self._vec_cache = self._build_vec_matrix(self.generators)
        self._dim_cache = int(np.linalg.matrix_rank(self._vec_cache, tol=rank_tol))

    # ---- properties ----

    @property
    def dim(self) -> int:
        """Dimension of A_circ as a complex linear subspace of M_N(C). Should
        equal n_points (the N elementary diagonals are linearly independent)."""
        return self._dim_cache

    @property
    def envelope_dim(self) -> int:
        """Dimension of the C*-envelope of A_circ as an *abstract* C*-algebra.

        Since A_circ_N is itself a finite-dimensional commutative C*-algebra
        isomorphic to C^N, its C*-envelope is itself, dim = N.

        If you want the *ambient*-envelope dimension (treating M_N(C) as the
        target), use ambient_envelope_dim.
        """
        return self.n_points

    @property
    def ambient_envelope_dim(self) -> int:
        """Ambient envelope dimension N^2.

        Used for the strict apples-to-apples comparison with
        TruncatedOperatorSystem.envelope_dim, which treats M_N(C) as the
        envelope.  A_circ_N never reaches dim = N^2 by powers (since it is
        always commutative), so the resulting "ambient propagation" is
        infinity.
        """
        return self.n_points ** 2

    # ---- helpers ----

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

    # ---- self-checks ----

    def identity_in_algebra(self, *, tol: float = 1e-10) -> Tuple[bool, float]:
        """The identity is sum_i E_{i, i}, exactly in A_circ_N for any N."""
        I = np.eye(self.n_points, dtype=np.complex128)
        target = I.reshape(-1).astype(np.complex128)
        norm_target = np.linalg.norm(target)
        if norm_target < 1e-30:
            return True, 0.0
        c, *_ = np.linalg.lstsq(self._vec_cache, target, rcond=None)
        residual = np.linalg.norm(self._vec_cache @ c - target)
        ratio = float(residual / norm_target)
        return ratio < tol, ratio

    def verify_multiplicative_closure(self, *, tol: float = 1e-12) -> Tuple[bool, List[Tuple[int, int]]]:
        """Verify that for every pair of generators (E_i, E_j), the product
        E_i E_j lies in A_circ_N.

        Since E_i E_j = delta_{i, j} E_i for elementary diagonal matrices,
        this is exact at machine precision.  Returns (all_closed, failures).
        """
        failures: List[Tuple[int, int]] = []
        for i, A in enumerate(self.generators):
            for j, B in enumerate(self.generators):
                prod = A @ B
                target = prod.reshape(-1).astype(np.complex128)
                norm_target = np.linalg.norm(target)
                if norm_target < 1e-30:
                    continue  # zero product, trivially in algebra
                c, *_ = np.linalg.lstsq(self._vec_cache, target, rcond=None)
                residual = np.linalg.norm(self._vec_cache @ c - target)
                if residual / norm_target > tol:
                    failures.append((i, j))
        return (len(failures) == 0, failures)

    def verify_star_closure(self, *, tol: float = 1e-12) -> Tuple[bool, List[int]]:
        """Verify that for every generator E_i, E_i^* lies in A_circ_N.

        Since E_i is real-diagonal, E_i^* = E_i exactly.  Returns
        (all_closed, failures).
        """
        failures: List[int] = []
        for i, A in enumerate(self.generators):
            adj = A.conj().T
            target = adj.reshape(-1).astype(np.complex128)
            norm_target = np.linalg.norm(target)
            if norm_target < 1e-30:
                continue
            c, *_ = np.linalg.lstsq(self._vec_cache, target, rcond=None)
            residual = np.linalg.norm(self._vec_cache @ c - target)
            if residual / norm_target > tol:
                failures.append(i)
        return (len(failures) == 0, failures)

    # ---- propagation number ----

    def compute_propagation_number(
        self,
        *,
        max_k: int = 5,
        tol: float = 1e-10,
        verbose: bool = False,
        envelope: str = "intrinsic",
    ) -> CirculantPropagationResult:
        """Compute prop(A_circ_N) using the same vec-stack rank algorithm as
        TruncatedOperatorSystem (apples-to-apples comparison).

        Parameters
        ----------
        envelope : {"intrinsic", "ambient"}
            "intrinsic" : treat A_circ_N (as an abstract C*-algebra) as its
                own envelope.  Then the target dim is N (= self.envelope_dim).
                Since A_circ_N is multiplicatively closed, dim(A^1) = N
                = target, so prop = 1.
            "ambient"   : treat the ambient matrix algebra M_N(C) as the
                envelope (mirroring the convention used for
                TruncatedOperatorSystem).  Then the target dim is N^2
                (= self.ambient_envelope_dim).  Since A_circ_N is commutative
                and never fills M_N(C) for N > 1, this gives prop = infinity
                (returned as -1).

        Returns
        -------
        CirculantPropagationResult
            (prop, dim_sequence) — same format as
            geovac.operator_system.propagation_number().
        """
        if envelope not in {"intrinsic", "ambient"}:
            raise ValueError(f"envelope must be 'intrinsic' or 'ambient', got {envelope!r}")

        target_dim = self.envelope_dim if envelope == "intrinsic" else self.ambient_envelope_dim

        # k = 1
        dim_A1 = operator_system_dim(self.generators, tol=tol)
        dim_seq = [dim_A1]
        if verbose:
            print(
                f"  [circulant N={self.n_points}, env={envelope}] "
                f"k=1: |gens|={len(self.generators)}, dim(A^1)={dim_A1}, target={target_dim}"
            )
        if dim_A1 == target_dim:
            return CirculantPropagationResult(prop=1, dim_sequence=dim_seq)

        # k >= 2: build A^k by multiplying basis of A^{k-1} on the right by A.
        from geovac.operator_system import _extract_matrix_basis  # internal helper, OK to reuse

        current_generators = list(self.generators)
        for k in range(2, max_k + 1):
            basis_Ak_minus_1 = _extract_matrix_basis(current_generators, tol=tol)
            new_generators = [A @ B for A in basis_Ak_minus_1 for B in self.generators]
            dim_Ak = operator_system_dim(new_generators, tol=tol)
            dim_seq.append(dim_Ak)
            if verbose:
                print(
                    f"  [circulant N={self.n_points}, env={envelope}] "
                    f"k={k}: |basis(A^{k-1})|={len(basis_Ak_minus_1)}, "
                    f"|gens|={len(new_generators)}, dim(A^{k})={dim_Ak}, target={target_dim}"
                )
            if dim_Ak == target_dim:
                return CirculantPropagationResult(prop=k, dim_sequence=dim_seq)
            if dim_Ak == dim_seq[-2]:
                # Saturated below the target: prop = infinity (consistent with
                # commutativity blocking ambient envelope when envelope='ambient').
                return CirculantPropagationResult(prop=-1, dim_sequence=dim_seq)
            current_generators = new_generators

        return CirculantPropagationResult(prop=-1, dim_sequence=dim_seq)


# ---------------------------------------------------------------------------
# Convenience: build a circulant comparator dimensioned to a GeoVac case
# ---------------------------------------------------------------------------


def circulant_for_geovac(n_max: int, *, match: str = "envelope") -> CirculantS3Truncation:
    """Build a CirculantS3Truncation with dimension matched to a GeoVac case.

    Parameters
    ----------
    n_max : int
        The GeoVac Fock cutoff to compare against.
    match : {"envelope", "operator_system"}
        "envelope":        Pick N(circulant) = N(GeoVac) = sum_{n=1..n_max} n^2.
            This makes the *ambient* M_N(C) match GeoVac's M_N(C) exactly.
        "operator_system": Pick N(circulant) = dim(O_{n_max} GeoVac).
            This makes the dimensions of A_circ and O_{n_max} match as
            *abstract C*-subalgebras* / operator systems.

    Returns
    -------
    CirculantS3Truncation
    """
    if match == "envelope":
        N = hilbert_dim(n_max)
    elif match == "operator_system":
        # We have to construct GeoVac to read off dim(O); this is a one-shot
        # cost.
        gv = TruncatedOperatorSystem(n_max)
        N = gv.dim
    else:
        raise ValueError(f"match must be 'envelope' or 'operator_system', got {match!r}")
    return CirculantS3Truncation(n_points=N)


# ---------------------------------------------------------------------------
# Cross-comparison helper
# ---------------------------------------------------------------------------


def compare_to_geovac(
    n_max: int,
    *,
    match: str = "envelope",
    verbose: bool = False,
) -> dict:
    """Build both the GeoVac TruncatedOperatorSystem at n_max and the
    matched-dim circulant analog, compute both propagation numbers, and return
    a comparison dict.

    Parameters
    ----------
    n_max : int
        GeoVac Fock cutoff.
    match : str
        See circulant_for_geovac.
    verbose : bool
        If True, print intermediate steps of both prop computations.

    Returns
    -------
    dict with keys
        "n_max", "match",
        "geovac_N", "geovac_dim_O", "geovac_envelope_dim", "geovac_prop",
        "geovac_dim_sequence",
        "circulant_N", "circulant_dim", "circulant_intrinsic_envelope",
        "circulant_ambient_envelope",
        "circulant_intrinsic_prop", "circulant_intrinsic_dim_sequence",
        "circulant_ambient_prop",  "circulant_ambient_dim_sequence",
        "structural_difference"  (boolean: True iff geovac_prop != 1 and
            circulant_intrinsic_prop == 1, i.e. the round-3 falsification is
            CONFIRMED in the apples-to-apples sense.)
    """
    if verbose:
        print(f"\n=== compare_to_geovac(n_max={n_max}, match={match!r}) ===\n")
        print("Building GeoVac TruncatedOperatorSystem...")

    gv = TruncatedOperatorSystem(n_max)
    if verbose:
        print(f"  GeoVac: N={gv.dim_H}, dim(O)={gv.dim}, envelope=N^2={gv.envelope_dim}")
        print("Computing GeoVac prop...")

    gv_prop, gv_dimseq = propagation_number(gv, verbose=verbose)

    if verbose:
        print(f"\nBuilding circulant comparator (match={match!r})...")

    circ = circulant_for_geovac(n_max, match=match)
    if verbose:
        print(
            f"  Circulant: N={circ.n_points}, dim(A_circ)={circ.dim}, "
            f"intrinsic envelope=N={circ.envelope_dim}, ambient envelope=N^2={circ.ambient_envelope_dim}"
        )
        print("Computing circulant prop (intrinsic envelope)...")

    circ_intrinsic = circ.compute_propagation_number(envelope="intrinsic", verbose=verbose)

    if verbose:
        print("Computing circulant prop (ambient envelope)...")

    circ_ambient = circ.compute_propagation_number(envelope="ambient", verbose=verbose)

    structural_difference = (gv_prop != 1) and (circ_intrinsic.prop == 1)

    return {
        "n_max": n_max,
        "match": match,
        "geovac_N": gv.dim_H,
        "geovac_dim_O": gv.dim,
        "geovac_envelope_dim": gv.envelope_dim,
        "geovac_prop": gv_prop,
        "geovac_dim_sequence": gv_dimseq,
        "circulant_N": circ.n_points,
        "circulant_dim": circ.dim,
        "circulant_intrinsic_envelope": circ.envelope_dim,
        "circulant_ambient_envelope": circ.ambient_envelope_dim,
        "circulant_intrinsic_prop": circ_intrinsic.prop,
        "circulant_intrinsic_dim_sequence": circ_intrinsic.dim_sequence,
        "circulant_ambient_prop": circ_ambient.prop,
        "circulant_ambient_dim_sequence": circ_ambient.dim_sequence,
        "structural_difference": structural_difference,
    }


# ---------------------------------------------------------------------------
# CLI / direct invocation: print a comparison table
# ---------------------------------------------------------------------------


def _print_comparison_table(n_max_values: Sequence[int] = (2, 3, 4)) -> None:
    """Print a side-by-side comparison table for the round-3 memo."""
    print("\nWH1 Round 3 — Falsification table")
    print("=" * 100)
    print(
        f"{'n_max':>5} | {'N(GV)':>5} | {'dim(O_GV)':>9} | {'prop(O_GV)':>10} || "
        f"{'N(circ)':>7} | {'dim(A_circ)':>11} | {'prop intrinsic':>14} | {'prop ambient':>12}"
    )
    print("-" * 100)
    for n_max in n_max_values:
        result = compare_to_geovac(n_max, match="envelope")
        prop_amb_str = "inf" if result["circulant_ambient_prop"] == -1 else str(result["circulant_ambient_prop"])
        print(
            f"{result['n_max']:>5} | {result['geovac_N']:>5} | {result['geovac_dim_O']:>9} | "
            f"{result['geovac_prop']:>10} || "
            f"{result['circulant_N']:>7} | {result['circulant_dim']:>11} | "
            f"{result['circulant_intrinsic_prop']:>14} | {prop_amb_str:>12}"
        )
    print("=" * 100)
    print(
        "\nReading: GeoVac O_{n_max} has prop = 2 (matching Toeplitz S^1).  Circulant\n"
        "A_circ_N is a commutative C*-algebra and has intrinsic prop = 1 (since A^1 is\n"
        "already its own envelope).  In the strict ambient-envelope convention, A_circ\n"
        "never fills M_N(C) and so has prop = infinity.  Either reading confirms that\n"
        "prop(GeoVac) = 2 is structurally specific to the Connes-vS construction and\n"
        "NOT a generic feature of any finite truncation of S^3.\n"
    )


if __name__ == "__main__":
    _print_comparison_table()
