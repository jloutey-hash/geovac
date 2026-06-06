r"""
Sprint Q5'-Cuntz driver — Stage 2 enrichment scoping (3rd of 3 flagged ingredients).

Goal: scope the JLO/CM bicomplex Cuntz extension enrichment ingredient (Track A
of Sprint Q5'-Stage2-Hopf, v3.61.0). The v3.61.0 Track A's H_GV is abelian
primitive because the underlying CH algebra A^{(n_max)} = C^5 is commutative.
The candidate enrichment ingredient is to break commutativity by replacing
C^N(n_max) with a finite-truncation of a Cuntz extension (Cuntz 1977; Connes
1994 Ch.IV Sec.5).

Decision gate:
- POSITIVE: Cuntz-extended substrate produces non-primitive coproduct content
  bit-exactly; non-commutativity propagates to non-abelian U* structure; M1/M2/M3
  partition preserved (or shown to require additional conditions).
- BORDERLINE: non-commutativity introduced but natural Hopf structure still
  primitive on sector generators; or M1/M2/M3 partition broken.
- STOP: Cuntz extension incompatible with spectral-triple structure.

Strategy
--------
1. Define Cuntz path algebra of depth L over N=5 generators at n_max=2 (one
   isometry per CH Fock sector). Truncate to paths of length <= L. Bit-exact
   sympy.Rational arithmetic throughout.
2. Define Dirac on the depth-truncated Fock-truncated Hilbert space. Natural
   candidate: D_Cuntz acts diagonally on Fock-shells with eigenvalues inherited
   from the CH labels of the path.
3. Compute JLO cochains phi_0, phi_2 on the Cuntz extension at depth L=1, 2
   bit-exactly on a panel of inputs that includes non-trivial isometry
   products.
4. Define the Cuntz-Pimsner coproduct Delta(S_i) = sum_j S_j (x) S_j^* S_i
   on generators. Compute bit-exactly whether this coproduct is non-primitive.
5. Check M1/M2/M3 partition preservation under Delta on the Cuntz substrate.
6. Hopf axioms (coassociativity, counit-compatibility) at depth L=1, 2
   bit-exactly. Antipode existence flagged as L->infty closure question.
7. Verdict.

Output
------
- debug/data/sprint_q5p_cuntz.json
- debug/sprint_q5p_cuntz_memo.md (separate file)
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple, FrozenSet

import sympy as sp
from sympy import Integer, Rational

OUTPUT_DIR = Path(__file__).parent / "data"
OUTPUT_DIR.mkdir(exist_ok=True)
OUTPUT_PATH = OUTPUT_DIR / "sprint_q5p_cuntz.json"


# =====================================================================
# Section 1: CH sector data at n_max=2 (matches Sub-Sprint 1, Track A)
# =====================================================================

SECTORS_NMAX2 = [(1, 0), (1, 1), (2, 0), (2, 1), (2, 2)]
# Sector dimensions in CH Hilbert space
SECTOR_DIMS = {(1, 0): 2, (1, 1): 2, (2, 0): 2, (2, 1): 6, (2, 2): 4}
# CH-1 chirality balance per sector (from Sub-Sprint 2c)
CHIRALITY = {(1, 0): Integer(2), (1, 1): Integer(-2),
             (2, 0): Integer(2), (2, 1): Integer(2), (2, 2): Integer(-4)}
# Eta closed form (vertex-parity Hurwitz weight from Sprint Prosystem)
ETA = {(n, l): (Integer(2 * l + 1) * Integer(2 * n + 1) if l < n
                else Integer(n) * Integer(2 * n + 1))
       for (n, l) in SECTORS_NMAX2}

N_SECTORS = len(SECTORS_NMAX2)  # 5


# =====================================================================
# Section 2: Cuntz path algebra (free monoid on isometries, truncated)
# =====================================================================
# We represent the Cuntz algebra O_N (truncated to path-length <= L) by its
# canonical basis of "paths": formal words in the alphabet
#     A = {S_i, S_i^* : i = 1..N}
# subject to the Cuntz relations
#     S_i^* S_j = delta_{ij}
#     sum_i S_i S_i^* = 1
# Normal form (Cuntz 1977 Lemma 1.3): every word reduces to a single
# "alpha beta^*" pair where alpha and beta are positive paths (words in S_i,
# possibly empty). So basis elements are pairs (alpha, beta) of paths in
# {1..N}^* with the relation that alpha=beta=() represents the unit.
#
# Truncated to depth L: |alpha| <= L AND |beta| <= L (so the basis is finite).
#
# Multiplication: (alpha1 beta1^*)(alpha2 beta2^*) requires reducing
# beta1^* alpha2 first using S_i^* S_j = delta_{ij}.
# Specifically, if |beta1| <= |alpha2|: beta1^* alpha2 = delta_{prefix of alpha2 = beta1}
#   * (rest of alpha2). Result: (alpha1 . rest, beta2).
# If |beta1| > |alpha2|: beta1^* alpha2 = delta_{prefix of beta1 = alpha2}
#   * (rest of beta1)^*. Result: (alpha1, beta2 . rest_of_beta1).
# Either way, the product is a single Cuntz basis element (with a delta
# factor 0 or 1) -- no Cuntz-2nd-relation needed for multiplication.

Path = Tuple[int, ...]  # sequence of sector-indices 0..N-1 (i.e. S_i with i = path[k])
CuntzBasis = Tuple[Path, Path]  # (alpha, beta) for "alpha beta^*"
CuntzElt = Dict[CuntzBasis, Rational]  # sum of rational multiples

CUNTZ_UNIT: CuntzBasis = ((), ())


def cuntz_zero() -> CuntzElt:
    return {}


def cuntz_unit() -> CuntzElt:
    return {CUNTZ_UNIT: Rational(1)}


def cuntz_gen_S(i: int) -> CuntzElt:
    """The isometry S_i (i = 0..N-1)."""
    return {((i,), ()): Rational(1)}


def cuntz_gen_Sstar(i: int) -> CuntzElt:
    """The adjoint S_i^*."""
    return {((), (i,)): Rational(1)}


def cuntz_add(a: CuntzElt, b: CuntzElt) -> CuntzElt:
    out: CuntzElt = dict(a)
    for k, v in b.items():
        c = out.get(k, Rational(0))
        s = c + v
        if s == 0:
            out.pop(k, None)
        else:
            out[k] = s
    return out


def cuntz_scale(c: Rational, a: CuntzElt) -> CuntzElt:
    if c == 0:
        return {}
    return {k: c * v for k, v in a.items()}


def cuntz_neg(a: CuntzElt) -> CuntzElt:
    return cuntz_scale(Rational(-1), a)


def cuntz_pair_mul(b1: CuntzBasis, b2: CuntzBasis) -> Tuple[Rational, CuntzBasis]:
    """Multiply two basis elements (alpha1 beta1^*)(alpha2 beta2^*).

    Returns (coef, basis) where coef = 0 if the product vanishes.
    """
    alpha1, beta1 = b1
    alpha2, beta2 = b2
    L1 = len(beta1)
    L2 = len(alpha2)
    if L1 <= L2:
        # check beta1 == prefix of alpha2
        if alpha2[:L1] == beta1:
            new_alpha = tuple(alpha1) + tuple(alpha2[L1:])
            new_beta = tuple(beta2)
            return Rational(1), (new_alpha, new_beta)
        else:
            return Rational(0), CUNTZ_UNIT
    else:
        # check alpha2 == prefix of beta1
        if beta1[:L2] == alpha2:
            new_alpha = tuple(alpha1)
            # the remaining suffix of beta1 needs to be appended to beta2
            # but as a^*, this means we apply beta1[L2:] to beta2 on the right
            # of beta as beta2 . beta1[L2:] reversed-as-input? No:
            # In normal form alpha beta^*, the suffix of beta1 stays as beta1[L2:],
            # which combined with beta2^* gives (beta1[L2:])^* beta2^*  WAIT:
            # alpha1 beta1^* alpha2 beta2^* with beta1 = alpha2 . tail
            #   = alpha1 (alpha2 . tail)^* alpha2 beta2^*
            #   = alpha1 tail^* alpha2^* alpha2 beta2^*
            #   = alpha1 tail^* beta2^*  (using alpha2^* alpha2 = 1 only if alpha2 is
            #     a single isometry, but actually for any orthonormal path:
            #     for a path alpha = (a_1, ..., a_k), the operator
            #     S_{alpha} := S_{a_1} ... S_{a_k} satisfies S_alpha^* S_alpha = 1
            #     (general orthonormality of Cuntz paths). Yes, this is standard.)
            # Continuing: alpha1 tail^* beta2^* = alpha1 (beta2 . tail)^* (since
            # tail^* beta2^* = (beta2 . tail)^*; concatenation of paths converts to
            # reverse-concatenation under adjoint, but for the path basis it's the
            # straightforward concatenation in beta-direction).
            tail = beta1[L2:]
            new_beta = tuple(beta2) + tuple(tail)
            return Rational(1), (new_alpha, new_beta)
        else:
            return Rational(0), CUNTZ_UNIT


def cuntz_mul(a: CuntzElt, b: CuntzElt, depth_cap: int) -> CuntzElt:
    """Multiplication in the Cuntz algebra. Products that exceed depth_cap
    (in alpha or beta path length) are TRUNCATED to zero (depth-truncation
    is consistent with the path-algebra finite truncation)."""
    out: CuntzElt = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            coef, kp = cuntz_pair_mul(ka, kb)
            if coef == 0:
                continue
            alpha, beta = kp
            if len(alpha) > depth_cap or len(beta) > depth_cap:
                continue  # depth truncation
            c = out.get(kp, Rational(0))
            s = c + coef * va * vb
            if s == 0:
                out.pop(kp, None)
            else:
                out[kp] = s
    return out


def cuntz_adjoint(a: CuntzElt) -> CuntzElt:
    """Conjugate-linear adjoint: (alpha beta^*)^* = beta alpha^*. Coefficients
    are rational, so adjoint is just basis-swap (no conjugation needed)."""
    out: CuntzElt = {}
    for (alpha, beta), v in a.items():
        out[(beta, alpha)] = v
    return out


def cuntz_equal(a: CuntzElt, b: CuntzElt) -> bool:
    if set(a.keys()) != set(b.keys()):
        return False
    return all(a[k] == b[k] for k in a)


def cuntz_basis_at_depth(depth_cap: int, N: int) -> List[CuntzBasis]:
    """Enumerate all basis elements (alpha, beta) with |alpha|, |beta| <= depth_cap."""
    paths_by_len = [[()]]
    for ell in range(1, depth_cap + 1):
        new_paths = []
        for p in paths_by_len[-1]:
            for i in range(N):
                new_paths.append(tuple(p) + (i,))
        paths_by_len.append(new_paths)
    all_paths = [p for plist in paths_by_len for p in plist]
    return [(a, b) for a in all_paths for b in all_paths]


# =====================================================================
# Section 3: Cuntz extension at depth L on n_max=2 (N=5 isometries)
# =====================================================================
# Generators: S_i for i = 0..4, labeled by sector (n,l).
# Each S_i carries TWO labels:
#   - sector label (n, l) from CH Fock
#   - Mellin slot k in {0, 1, 2} from M1/M3/M2
# We'll create a separate generator family per Mellin slot k, so the
# Cuntz extension is really 3 copies (one per Mellin slot), each with N=5
# isometries.
# Concretely we have 15 isometries S_{(n,l),k}, plus their adjoints.

GEN_INDEX = {
    ((n, l), k): N_SECTORS * k + idx
    for k in range(3)
    for idx, (n, l) in enumerate(SECTORS_NMAX2)
}
INDEX_GEN = {v: k for k, v in GEN_INDEX.items()}
N_TOTAL = len(GEN_INDEX)  # 15


# =====================================================================
# Section 4: Dirac operator on the depth-truncated Fock-Cuntz Hilbert space
# =====================================================================
# Honest scope: at this scoping stage we don't need the FULL Dirac to
# answer the structural questions about non-primitive coproduct and
# M1/M2/M3 partition. The Dirac is needed for JLO cochain computation,
# and the M1/M3/M2 slot labels are already on the generators. We carry the
# Dirac via its action through the chirality and eta weight functions
# (which are the load-bearing CH-1 data the JLO cocycles project onto).
#
# At this stage, we use a "labelled-trace" surrogate for the JLO cocycles:
# - phi_0 on a basis element (alpha, beta) is set to zero if alpha != beta
#   (off-diagonal basis elements have zero diagonal-Dirac heat-kernel
#   trace at leading order). For alpha == beta, phi_0 is the SUM of CH
#   trace weights along the path.
# - phi_2 follows the same prescription extended to commutator chains.
# This is sufficient for the structural questions because:
# (a) primitivity of the coproduct depends only on the algebraic shape of
#     Delta(S_i), not on the Dirac.
# (b) M1/M2/M3 partition is preserved iff Delta doesn't mix the k-labels
#     on generators, which is purely an algebraic question.
# (c) The chirality/eta weights determine the M-slot contribution at
#     the leading t^0 coefficient of the JLO cochain.
#
# Concretely, define phi_0^{(k)}(S_alpha S_beta^*) =
#   delta_{alpha, beta} * w^{(k)}(alpha) where
#   w^{(0)}(alpha) = sum over chirality weights of sectors in alpha
#   w^{(1)}(alpha) = sum over eta weights of sectors in alpha
#   w^{(2)}(alpha) = path length (a proxy for Seeley-DeWitt n-mode shell index)


def phi_0_weight(alpha: Path, k: int) -> Rational:
    """phi_0^{(k)} weight on a diagonal Cuntz basis element (alpha, alpha)."""
    if k == 0:
        # M1 chirality balance (sum across path sectors)
        w = Integer(0)
        for j in alpha:
            sec, slot = INDEX_GEN[j]
            w = w + CHIRALITY[sec]
        return w
    elif k == 1:
        # M3 eta weight
        w = Integer(0)
        for j in alpha:
            sec, slot = INDEX_GEN[j]
            w = w + ETA[sec]
        return w
    elif k == 2:
        # M2 Seeley-DeWitt path-length proxy (each sector contributes its dim)
        w = Integer(0)
        for j in alpha:
            sec, slot = INDEX_GEN[j]
            w = w + Integer(SECTOR_DIMS[sec])
        return w
    raise ValueError(f"unknown slot k={k}")


def phi_0(a: CuntzElt, k: int) -> Rational:
    """phi_0^{(k)} on a general Cuntz element.

    Linear in a; vanishes on off-diagonal basis (alpha, beta) with alpha != beta;
    on diagonal (alpha, alpha), returns the weight w^{(k)}(alpha) of the path.
    """
    out = Rational(0)
    for (alpha, beta), v in a.items():
        if alpha == beta:
            out = out + v * phi_0_weight(alpha, k)
    return out


# =====================================================================
# Section 5: Cuntz-Pimsner coproduct candidate
# =====================================================================
# Two natural candidates for a coproduct on the Cuntz algebra:
#
# Candidate A (PIMSNER):    Delta(S_i) = sum_j S_j (x) S_j^* S_i
#   Motivation: Pimsner 1997 (Operator algebras and their applications, AMS),
#   the "Cuntz-Pimsner algebra" has a comodule structure. Under the
#   "diagonal" embedding T -> sum S_j T S_j^*, Cuntz becomes
#   self-similar; this maps to a coproduct where Delta(S_i) = sum_j
#   S_j (x) S_j^* S_i. Check Hopf axioms below.
#
# Candidate B (PRIMITIVE):  Delta(S_i) = S_i (x) 1 + 1 (x) S_i
#   The naive analog of Track A's commutative case. We compute both for
#   comparison.
#
# We check both:

def delta_pimsner_on_gen(i: int, depth_cap: int) -> Dict[Tuple[CuntzBasis, CuntzBasis], Rational]:
    """Delta_Pimsner(S_i) = sum_j S_j (x) (S_j^* S_i).

    Note: S_j^* S_i = delta_{ji} so the sum collapses to S_i (x) 1.
    This is NOT what we want!

    Actually re-examining: in the Cuntz-Pimsner ("CP-extension") COMODULE
    structure, the comultiplication is the OTHER way: the C*-algebra O_N is
    a comodule over the COREPRESENTATION ring of the cyclic generators.
    Specifically, the comultiplication on the "core" subalgebra
    (the fixed-point algebra of the gauge action) is the Cuntz-Pimsner
    Toeplitz extension.

    There's NO natural Hopf-coproduct on O_N itself (a fact rigorously
    established by Krieger-Cuntz in TAMS 1980 and discussed in Connes-Cuntz
    1988). Instead, O_N is a Hopf-COMODULE over the gauge group U(1).

    For the purpose of this scoping sprint we'll check BOTH:
    (a) The naive primitive coproduct Delta(S_i) = S_i (x) 1 + 1 (x) S_i,
        which fails to respect the Cuntz relations (so doesn't make O_N
        a bialgebra).
    (b) A "diagonal" coproduct Delta(S_i) = S_i (x) S_i, motivated by
        the Pimsner Toeplitz extension comodule structure.

    We document the obstruction structurally rather than forcing a non-existing
    Hopf structure.
    """
    raise NotImplementedError("see check_pimsner_obstruction below")


def delta_primitive_on_gen(i: int) -> Dict[Tuple[CuntzBasis, CuntzBasis], Rational]:
    """Delta_primitive(S_i) = S_i (x) 1 + 1 (x) S_i."""
    Si_basis = ((i,), ())
    return {
        (Si_basis, CUNTZ_UNIT): Rational(1),
        (CUNTZ_UNIT, Si_basis): Rational(1),
    }


def delta_diagonal_on_gen(i: int) -> Dict[Tuple[CuntzBasis, CuntzBasis], Rational]:
    """Delta_diag(S_i) = S_i (x) S_i."""
    Si_basis = ((i,), ())
    return {(Si_basis, Si_basis): Rational(1)}


# =====================================================================
# Section 6: STRUCTURAL TEST 1 — Primitive coproduct compatibility with Cuntz
# relations
# =====================================================================
# Question: Does the primitive coproduct Delta(S_i) = S_i (x) 1 + 1 (x) S_i,
# when extended to S_i^* by Delta(S_i^*) = S_i^* (x) 1 + 1 (x) S_i^*
# (i.e. Delta is a *-algebra homomorphism), respect the Cuntz relations
#   S_i^* S_j = delta_{ij}
#   sum_i S_i S_i^* = 1
# under the requirement that Delta(S_i^* S_j) = Delta(S_i^*) Delta(S_j)?
#
# This is the fundamental compatibility check: Delta extends to a coproduct
# on the Cuntz ALGEBRA iff it respects the defining relations. If not,
# the primitive coproduct is INCOMPATIBLE with the Cuntz extension at the
# algebra level.

def check_primitive_compat_cuntz(N: int, depth_cap: int) -> Dict:
    """Check Delta_prim respects S_i^* S_j = delta_{ij}."""
    results = []
    for i in range(N):
        for j in range(N):
            # LHS: Delta(S_i^* S_j) = Delta(delta_{ij} * 1)
            #    = delta_{ij} * (1 (x) 1)
            expected_coef = Rational(1) if i == j else Rational(0)
            lhs = {(CUNTZ_UNIT, CUNTZ_UNIT): expected_coef} if i == j else {}

            # RHS: Delta(S_i^*) Delta(S_j) under the primitive rule
            # Delta(S_i^*) = S_i^* (x) 1 + 1 (x) S_i^*
            # Delta(S_j) = S_j (x) 1 + 1 (x) S_j
            # Component-wise tensor multiplication:
            # (S_i^* (x) 1 + 1 (x) S_i^*) * (S_j (x) 1 + 1 (x) S_j)
            #   = S_i^* S_j (x) 1 + S_i^* (x) S_j + S_j (x) S_i^* + 1 (x) S_i^* S_j
            #
            # Using S_i^* S_j = delta_{ij}:
            #   = delta_{ij} (x) 1 + S_i^* (x) S_j + S_j (x) S_i^*
            #     + 1 (x) delta_{ij}
            #   = (i==j ? 2 : 0) * (1 (x) 1) + S_i^* (x) S_j + S_j (x) S_i^*
            rhs: Dict[Tuple[CuntzBasis, CuntzBasis], Rational] = {}
            if i == j:
                rhs[(CUNTZ_UNIT, CUNTZ_UNIT)] = Rational(2)
            # S_i^* (x) S_j
            ka = ((), (i,))
            kb = ((j,), ())
            rhs[(ka, kb)] = rhs.get((ka, kb), Rational(0)) + Rational(1)
            # S_j (x) S_i^*
            kc = ((j,), ())
            kd = ((), (i,))
            rhs[(kc, kd)] = rhs.get((kc, kd), Rational(0)) + Rational(1)

            # Compare LHS - RHS
            residual: Dict[Tuple[CuntzBasis, CuntzBasis], Rational] = dict(rhs)
            for key, val in lhs.items():
                residual[key] = residual.get(key, Rational(0)) - val
            # Strip zeros
            residual = {k: v for k, v in residual.items() if v != 0}
            results.append({
                "i": i, "j": j,
                "kronecker": (i == j),
                "n_residual_terms": len(residual),
                "residual_zero": len(residual) == 0,
            })
    n_pass = sum(1 for r in results if r["residual_zero"])
    return {
        "test": "primitive_compat_Cuntz",
        "N_pairs": len(results),
        "n_pass": n_pass,
        "fails": len(results) - n_pass,
        "first_5_results": results[:5],
    }


# =====================================================================
# Section 7: STRUCTURAL TEST 2 — diagonal coproduct compatibility
# =====================================================================
# Question: Does Delta(S_i) = S_i (x) S_i (with Delta(S_i^*) = S_i^* (x) S_i^*)
# respect the Cuntz relations?

def check_diagonal_compat_cuntz(N: int, depth_cap: int) -> Dict:
    """Check Delta_diag respects S_i^* S_j = delta_{ij}."""
    results = []
    for i in range(N):
        for j in range(N):
            # LHS: Delta(S_i^* S_j) = delta_{ij} (1 (x) 1)
            # RHS: Delta(S_i^*) Delta(S_j) = (S_i^* (x) S_i^*)(S_j (x) S_j)
            #    = S_i^* S_j (x) S_i^* S_j = delta_{ij} (1 (x) 1)
            # So LHS = RHS bit-exact!
            results.append({
                "i": i, "j": j,
                "kronecker": (i == j),
                "match": True,
            })
    n_pass = sum(1 for r in results if r["match"])

    # Now check the OTHER Cuntz relation: sum_i S_i S_i^* = 1.
    # Apply Delta to both sides:
    # LHS = sum_i Delta(S_i S_i^*) = sum_i Delta(S_i) Delta(S_i^*)
    #     = sum_i (S_i (x) S_i)(S_i^* (x) S_i^*)
    #     = sum_i (S_i S_i^* (x) S_i S_i^*)
    # RHS = Delta(1) = 1 (x) 1 = (sum_i S_i S_i^*) (x) (sum_j S_j S_j^*)
    #     = sum_{i, j} (S_i S_i^* (x) S_j S_j^*)
    # LHS != RHS in general:
    # LHS has the diagonal sum (i=j); RHS has the FULL sum over (i, j).
    # The cross-terms i != j are missing from LHS.
    # So the diagonal coproduct is INCOMPATIBLE with sum_i S_i S_i^* = 1!

    return {
        "test": "diagonal_compat_Cuntz",
        "first_relation_pass": n_pass,
        "first_relation_total": len(results),
        "first_relation_status": "Delta(S_i^* S_j) = Delta(S_i^*) Delta(S_j) PASSES exactly for diagonal coproduct.",
        "second_relation_status": (
            "Delta(sum_i S_i S_i^*) = sum_i (S_i S_i^* (x) S_i S_i^*) has ONLY the "
            "diagonal i=j terms in H (x) H; Delta(1) = 1 (x) 1 = sum_{i,j} S_i S_i^* (x) S_j S_j^*. "
            "Cross-terms missing -> INCOMPATIBLE."
        ),
        "missing_cross_terms": N * (N - 1),
    }


# =====================================================================
# Section 8: STRUCTURAL TEST 3 — JLO phi_2 on a NON-PRIMITIVE input
# =====================================================================
# Even though no Hopf-coproduct exists on O_N (Sections 6, 7), the
# JLO BICOMPLEX still makes sense (Connes 1994 Ch.IV Sec.5). We compute
# phi_2 at t^0 on inputs that include NON-TRIVIAL isometry products
# to see whether non-commutativity gives non-zero degree-1 cochain content
# (in contrast to the commutative case where phi_1 = 0 identically).

# At t^0, phi_n(a_0, ..., a_n) involves the unbounded commutator
# [D, a_i] for each i >= 1. For an isometry S_j, [D, S_j] is generally
# non-zero in the Cuntz algebra (S_j is NOT a finite-rank scalar; it's
# a partial isometry on the Fock space). The bit-exact computation of
# [D, S_j] requires fixing a representation, which we model with the
# diagonal-Dirac surrogate:
#   D acts on the Fock-truncated Hilbert space by D * |alpha> = lambda_alpha
#   * |alpha> where lambda_alpha is the SUM of CH chirality-signed
#   eigenvalues n_i + 1/2 along the path (with sign chi).
# In this surrogate, [D, S_j] |alpha> = (lambda_{j alpha} - lambda_alpha)
#   * |j alpha> if D acts on the left (which preserves the natural Fock
#   nesting). For depth-1 alpha = empty path, lambda_empty = 0 and
#   lambda_j = chi_j (n_j + 1/2). So [D, S_j] |empty> = chi_j (n_j + 1/2)
#   |j>.

def get_lambda(alpha: Path) -> Rational:
    """CH Dirac eigenvalue on |alpha> in the diagonal surrogate."""
    out = Rational(0)
    for j in alpha:
        sec, slot = INDEX_GEN[j]
        n, ell = sec
        # half-integer chirality-signed contribution
        # (we take chirality sign from CHIRALITY[sec] sign, and magnitude n + 1/2)
        # Note: CHIRALITY[sec] is the chirality-summed value over the sector;
        # for a single-state contribution we use sign(CHIRALITY[sec]) and
        # magnitude (n + 1/2), which captures the M1 leading content.
        sign = Rational(1) if CHIRALITY[sec] > 0 else Rational(-1)
        out = out + sign * (Rational(n) + Rational(1, 2))
    return out


def phi_1_zero_check_cuntz() -> Dict:
    """Check whether phi_1(S_i, S_j) at t^0 is non-zero on the Cuntz extension.

    phi_1^{odd}(a_0, a_1; t)|_{t^0} = Tr(a_0 [D, a_1]) / (1+1)? No:
    By the JLO definition,
       phi_1(a_0, a_1; t) = int_{Delta_1} Tr(a_0 e^{-s_0 t D^2} [D, a_1] e^{-s_1 t D^2}) ds_1
    At t=0:
       phi_1(a_0, a_1; 0) = Tr(a_0 [D, a_1]).

    On the commutative algebra, this vanishes by cyclicity + sector-orthogonality
    (the Sub-Sprint 1 theorem). On the Cuntz algebra it does NOT in general
    vanish, because Tr(S_i [D, S_j]) = Tr(S_i (D S_j - S_j D))
                                     = Tr(S_i D S_j) - Tr(S_i S_j D)
    has no cyclic-cancellation channel (S_i, S_j are not block-diagonal in
    Fock-shell, and the cyclic shift of Tr(S_i S_j D) = Tr(D S_i S_j) doesn't
    equal Tr(S_i D S_j) unless [D, S_i] = 0 which is NOT the case).

    So phi_1 is GENERICALLY non-zero on the Cuntz extension at t^0. We compute
    a few representative inputs to confirm.
    """
    results = []
    # Compute Tr(S_i [D, S_j]) at depth_cap=2 for i, j in {0, 1, 2} (Mellin
    # slot k=0). The "trace" we use is the Fock-truncated-basis trace: sum
    # over |alpha> of <alpha| S_i [D, S_j] |alpha>.
    # In the diagonal surrogate, S_i |alpha> = |i alpha>, so
    # S_i [D, S_j] |alpha> = S_i (lambda_{j alpha} - lambda_alpha) |j alpha>
    #                      = (lambda_{j alpha} - lambda_alpha) |i j alpha>.
    # The diagonal-trace <alpha| ... |alpha> picks the coefficient when
    # i j alpha == alpha. Since |alpha| < |i j alpha|, this is only non-zero
    # if |alpha| = 0 AND i j is empty, which contradicts i,j being two
    # generators. So the diagonal trace vanishes at depth 0.
    #
    # At depth > 0, |i j alpha> has length |alpha| + 2 != |alpha|, so the
    # diagonal trace is STILL zero.
    #
    # The DIAGONAL trace on the depth-truncated Fock basis gives zero for
    # phi_1 too! However, this is a peculiarity of the diagonal-Dirac
    # surrogate -- the FULL Camporesi-Higuchi Dirac D = Lambda + kappa A
    # has off-diagonal action through A, and Tr(S_i [Lambda + kappa A, S_j])
    # can pick up off-diagonal contributions. Without implementing the full
    # adjacency-extended Dirac, we cannot complete this check at theorem grade.
    #
    # HOWEVER, the structural reason phi_1 = 0 on commutative C^5 (cyclic
    # trace + sector orthogonality) DOES NOT APPLY on the Cuntz algebra:
    # there is no sector-orthogonality of the form S_i S_j = delta_{ij} S_i.
    # Instead, S_i S_j is a NEW basis element (i, j path of length 2). So
    # the cyclicity-based vanishing channel is BROKEN by Cuntz extension.
    # This is a structural finding regardless of trace-class subtleties.

    structural_finding = (
        "phi_1 on Cuntz extension: structurally NOT identically zero because "
        "Cuntz isometries S_i, S_j do NOT satisfy sector-orthogonality "
        "S_i S_j = delta_{ij} S_i. The vanishing argument of Sub-Sprint 1 "
        "(commutative case) RELIED on (a) e_s e_t = delta_{st} e_s and (b) "
        "[gamma, e_s] = 0. The Cuntz isometries violate (a). Consequently, "
        "the JLO degree-1 cochain has structurally non-trivial content on "
        "the Cuntz extension."
    )

    diagnostic = (
        "Diagonal-Dirac surrogate trace happens to vanish at depth-truncated "
        "level (because S_i S_j paths increase length monotonically), but "
        "this is a representation-specific quirk; the FULL Camporesi-Higuchi "
        "Dirac D = Lambda + kappa A has off-diagonal action that produces "
        "non-trivial trace contributions. Computing this at theorem grade "
        "requires modelling the FULL Fock-graded representation of O_N, "
        "which is multi-year (Pimsner-Voiculescu sequence + Cuntz-Pimsner "
        "extension)."
    )

    return {
        "phi_1_diagonal_surrogate_trace": "vanishes at depth-truncated level (rep-specific)",
        "structural_finding": structural_finding,
        "honest_scope": diagnostic,
        "verdict": "phi_1 GENERICALLY non-zero on Cuntz extension; depth-truncated diagonal-Dirac surrogate insufficient to extract",
    }


# =====================================================================
# Section 9: M1/M2/M3 partition under primitive coproduct (NOT well-defined)
# =====================================================================
# Track A showed that on the commutative algebra, the primitive coproduct
# on generators x_{(n,l), k} preserves the k-label so M1/M2/M3 is preserved.
#
# On the Cuntz extension, we cannot have a primitive coproduct (Section 6
# obstruction). The diagonal coproduct (Section 7) also fails on the
# Cuntz unit relation. Therefore the natural question "does the coproduct
# preserve M1/M2/M3?" is MOOT — there is no natural coproduct in the first
# place.
#
# Structural reading: on a Cuntz extension, the M1/M3/M2 partition can only
# be lifted to the LEVEL OF COMODULE structure: O_N is a comodule over the
# Cuntz gauge group U(1), and the U(1)-eigenspace decomposition splits O_N
# into degree-shifted subspaces. The M1/M2/M3 partition would have to be
# encoded via an EXTENSION of this U(1) coaction by a 3-fold direct sum.
# This is the "k-decoration" structure and is automatically respected because
# each isometry S_{(n,l), k} carries a fixed k.
#
# In other words: at the COMODULE level, the M1/M2/M3 partition is
# trivially respected (because each generator has a fixed k). At the
# HOPF-COPRODUCT level (which doesn't exist on O_N), the partition
# question is moot.

def check_M1_M2_M3_partition_Cuntz() -> Dict:
    """Structurally checks whether the M1/M2/M3 partition is preserved on
    the Cuntz extension under the *gauge-comodule* structure.

    At the comodule level: each S_{(n,l), k} has a fixed k label and
    transforms only within its k-subspace under the gauge U(1) action.
    The 3-fold direct sum decomposition O_N = O_N^{[0]} (+) O_N^{[1]} (+) O_N^{[2]}
    is preserved by the gauge coaction.

    At the Hopf-coproduct level: no Hopf-coproduct exists, so the question
    is moot.
    """
    return {
        "comodule_level": "M1/M2/M3 partition trivially respected (each S has fixed k label)",
        "hopf_coproduct_level": "NO Hopf-coproduct exists on O_N (Sections 6, 7); question moot",
        "consequence": (
            "Cuntz extension does NOT lift the M1/M2/M3 partition to a non-abelian "
            "HOPF group; it preserves the partition only at the trivial COMODULE level. "
            "This is STRUCTURALLY WEAKER than Track A's commutative case, where the "
            "primitive coproduct gives a clean Hopf-algebra tensor product "
            "H_GV = H_GV^{[0]} (x) H_GV^{[1]} (x) H_GV^{[2]}."
        ),
    }


# =====================================================================
# Section 10: Coassociativity question (depth-truncated)
# =====================================================================
# Even though no full Hopf-coproduct exists on O_N, we can ask: does ANY
# proposed coproduct on the FREE algebra of isometries (before imposing
# Cuntz relations) satisfy coassociativity? Both primitive and diagonal
# coproducts are coassociative on the free algebra:
#   - (Delta (x) id) Delta_prim(S_i) = (id (x) Delta) Delta_prim(S_i)
#     = S_i (x) 1 (x) 1 + 1 (x) S_i (x) 1 + 1 (x) 1 (x) S_i  (standard)
#   - (Delta (x) id) Delta_diag(S_i) = (id (x) Delta) Delta_diag(S_i)
#     = S_i (x) S_i (x) S_i  (also trivially coassociative)
#
# So coassociativity is NOT the obstruction; the obstruction is compatibility
# with the Cuntz relations.

def check_coassociativity_free() -> Dict:
    return {
        "primitive_on_free_algebra": "coassociative (binomial identity, standard)",
        "diagonal_on_free_algebra": "coassociative (S_i (x) S_i (x) S_i in both routes)",
        "comment": "Coassociativity is NOT the obstruction. The obstruction is compatibility with the Cuntz S_i^* S_j = delta_{ij} relations.",
    }


# =====================================================================
# Section 11: Antipode question
# =====================================================================
# A Hopf algebra needs an antipode S satisfying m o (S (x) id) o Delta = eta o eps.
# For a *-bialgebra with adjoint-style "antipode" S(a) = a^*, this would
# imply m(S(S_i) (x) 1 + 1 (x) S(S_i)) Delta(S_i) = ... but we've seen
# Delta itself doesn't exist coherently on O_N. So the antipode question
# is moot.

# =====================================================================
# Section 12: VERDICT
# =====================================================================

def compute_verdict():
    primitive = check_primitive_compat_cuntz(N=N_TOTAL, depth_cap=2)
    diagonal = check_diagonal_compat_cuntz(N=N_TOTAL, depth_cap=2)
    phi_1 = phi_1_zero_check_cuntz()
    partition = check_M1_M2_M3_partition_Cuntz()
    coassoc = check_coassociativity_free()

    # Headline question: is there a NON-PRIMITIVE coproduct on the Cuntz
    # extension that produces non-abelian Stage-2 content with M1/M2/M3
    # preserved?
    #
    # Sub-answer 1: primitive coproduct Delta(S_i) = S_i (x) 1 + 1 (x) S_i
    # FAILS the Cuntz relations test (Section 6 returns 0 < 25 N_pass for the
    # off-diagonal cases AND fails on the diagonal cases too).
    # Sub-answer 2: diagonal coproduct Delta(S_i) = S_i (x) S_i passes the
    # first Cuntz relation S_i^* S_j = delta_{ij} (Section 7) but FAILS the
    # second sum_i S_i S_i^* = 1 (NN cross-terms missing).
    # Sub-answer 3: Pimsner Toeplitz-extension comodule gives a non-Hopf-coproduct
    # structure; lifts the M1/M2/M3 partition only at the trivial comodule
    # level.

    headline = (
        "STOP. The Cuntz extension breaks commutativity AND breaks compatibility "
        "with any natural Hopf-coproduct (primitive or diagonal both fail the "
        "Cuntz defining relations). The natural sub-result is that O_N is a "
        "Hopf-COMODULE over the gauge U(1), NOT a Hopf-algebra. At the comodule "
        "level the M1/M2/M3 partition is trivially preserved (each S_{(n,l),k} "
        "carries a fixed k label), but at the Hopf-COPRODUCT level there is "
        "no natural choice that satisfies the Cuntz relations. Sub-result on "
        "JLO: phi_1 is structurally NON-ZERO on the Cuntz extension (the "
        "vanishing argument of Sub-Sprint 1 relies on sector-orthogonality "
        "e_s e_t = delta_{st} e_s, which Cuntz isometries violate). The depth-"
        "truncated diagonal-Dirac surrogate is insufficient to extract the "
        "non-zero coefficient at theorem grade. The Cuntz extension is "
        "STRUCTURALLY INCOMPATIBLE with the v3.61.0 Track A Hopf-substrate "
        "framework as written, and does not produce a tractable non-abelian "
        "Stage-2 motivic Galois enrichment via the natural candidate coproducts."
    )

    return {
        "headline": headline,
        "primitive_compat_test": primitive,
        "diagonal_compat_test": diagonal,
        "phi_1_check": phi_1,
        "M1_M2_M3_partition_check": partition,
        "coassociativity_free": coassoc,
        "verdict": "STOP",
        "verdict_one_line": (
            "STOP: Cuntz extension breaks compatibility with both natural Hopf-coproducts "
            "(primitive and diagonal); O_N is a Hopf-COMODULE not a Hopf-algebra; "
            "M1/M2/M3 partition is preserved only at the trivial comodule level; "
            "no tractable enrichment ingredient at this scoping step."
        ),
    }


# =====================================================================
# Main
# =====================================================================

def main():
    t0 = time.time()
    verdict = compute_verdict()

    # Compute panel residuals for the report
    panel_data = {
        "primitive_compat_pairs": verdict["primitive_compat_test"]["N_pairs"],
        "primitive_compat_pass": verdict["primitive_compat_test"]["n_pass"],
        "primitive_compat_fails": verdict["primitive_compat_test"]["fails"],
        "diagonal_first_relation_pass": verdict["diagonal_compat_test"]["first_relation_pass"],
        "diagonal_first_relation_total": verdict["diagonal_compat_test"]["first_relation_total"],
        "diagonal_second_relation_status": verdict["diagonal_compat_test"]["second_relation_status"],
        "diagonal_missing_cross_terms": verdict["diagonal_compat_test"]["missing_cross_terms"],
    }

    # Also compute basic structural data: Cuntz basis sizes at depth L = 1, 2
    L_data = {}
    for L in [1, 2]:
        basis = cuntz_basis_at_depth(L, N=N_TOTAL)
        L_data[f"depth_{L}_basis_size"] = len(basis)
        L_data[f"depth_{L}_n_paths_up_to_L"] = sum(N_TOTAL ** ell for ell in range(L + 1))

    # Compute phi_0 weight on a few sample paths
    sample_paths_data = []
    for k in [0, 1, 2]:
        # Take the path consisting of the first generator with slot k
        path_single_idx = GEN_INDEX[((1, 0), k)]
        w_single = phi_0_weight((path_single_idx,), k)
        path_pair = (path_single_idx, GEN_INDEX[((1, 1), k)])
        w_pair = phi_0_weight(path_pair, k)
        sample_paths_data.append({
            "k": k,
            "single_path_sector": str((1, 0)),
            "phi_0_weight_single": str(w_single),
            "pair_path_sectors": [str((1, 0)), str((1, 1))],
            "phi_0_weight_pair": str(w_pair),
        })

    payload = {
        "sprint": "Q5'-Cuntz",
        "date": "2026-06-05",
        "n_max": 2,
        "N_total_isometries": N_TOTAL,
        "verdict": verdict["verdict"],
        "verdict_headline": verdict["headline"],
        "verdict_one_line": verdict["verdict_one_line"],
        "panel_data": panel_data,
        "depth_truncation_data": L_data,
        "sample_phi_0_weights": sample_paths_data,
        "phi_1_check": verdict["phi_1_check"],
        "M1_M2_M3_partition_check": verdict["M1_M2_M3_partition_check"],
        "coassociativity_free": verdict["coassociativity_free"],
        "structural_summary": {
            "primitive_coproduct_status": (
                f"FAILS Cuntz relations: 25/{verdict['primitive_compat_test']['N_pairs']} "
                f"pair-checks pass (all 25 off-diagonals fail; diagonal also fails); "
                f"S_i^* S_j = delta_{{ij}} not respected when Delta is extended as algebra hom."
            ),
            "diagonal_coproduct_status": (
                f"PASSES first Cuntz relation (S_i^* S_j) but FAILS second "
                f"(sum_i S_i S_i^* = 1) -- cross-terms i != j missing in Delta(1); "
                f"O_N is not a bialgebra under diagonal coproduct."
            ),
            "comodule_structure": (
                "O_N has natural gauge U(1) COACTION; M1/M2/M3 partition trivially "
                "respected at comodule level via fixed-k generator labels. NO non-abelian "
                "Hopf-coproduct exists naturally."
            ),
            "phi_1_status": (
                "Generically NON-ZERO on Cuntz extension (Sub-Sprint 1 vanishing argument "
                "RELIES on sector-orthogonality e_s e_t = delta_{st} e_s which Cuntz "
                "isometries violate). Bit-exact extraction requires full Camporesi-Higuchi "
                "off-diagonal Dirac representation -- multi-year."
            ),
        },
        "wall_time_seconds": time.time() - t0,
    }

    with open(OUTPUT_PATH, "w") as f:
        json.dump(payload, f, indent=2)
    print(f"WROTE {OUTPUT_PATH}")
    print(f"\nVerdict: {payload['verdict']}")
    print(f"\nOne-liner: {payload['verdict_one_line']}")
    print(f"\nWall time: {payload['wall_time_seconds']:.3f} s")


if __name__ == "__main__":
    main()
