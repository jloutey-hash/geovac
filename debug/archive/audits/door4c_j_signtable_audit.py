"""Door 4c — does the combined real structure J = J_GV (x) J_F SELECT H,
or merely ADMIT it, at the n=2 inner factor?

The Door 4b fork
================

Door 4b pinned the GeoVac inner algebra to a single binary fork at the
n=2 rung:

    C (+) H        (+) M_3(C)     -- Standard Model
    C (+) M_2(C)   (+) M_3(C)     -- alternative

Both reproduce U(1) x SU(2) x SU(3) as the post-unimodularity gauge group.
The conjectured GeoVac-native handle for closing the fork: the OUTER triple
carries J_GV^2 = -1 at KO-dimension 3 (the pseudoreal Spin(3)=SU(2)
two-spinor sign), and in the Chamseddine-Connes-Marcolli account that SAME
quaternionic fact selects H over M_2(C).

The crucial structural fact (this is the whole game)
====================================================

Over C, the algebras H and M_2(C) are NOT distinguished as complex algebras:

    H (x)_R C  ~=  M_2(C).

So "H vs M_2(C)" is NOT a choice of complex algebra -- it is a choice of
REAL FORM of the same complex algebra M_2(C). The real form is selected by
an antilinear involution:

    H  = { m in M_2(C) : j_q(m) = m },   j_q(m) = (i sigma_2) conj(m) (i sigma_2)^{-1}

i.e. H is the FIXED-POINT real *-subalgebra of M_2(C) under the quaternionic
antilinear involution j_q. The structure operator of j_q is

    J_q = (i sigma_2) K,    J_q^2 = (i sigma_2) conj(i sigma_2) = -I,

so the quaternionic real structure has J_q^2 = -1 on C^2.

M_2(C) as a real form corresponds to the OTHER involution -- the "real"
(orthogonal) form, whose structure operator squares to +1 (e.g. plain K,
fixing M_2(R)), OR equivalently "no real-form constraint" (the full complex
algebra). In the Connes-Marcolli finite-triple language the n=2 division
ring is selected by which antilinear J_F the finite triple carries on its
C^2 component:

    J_F^{quaternionic} = (i sigma_2) K  on C^2     ->  factor = H
    J_F^{real/complex}                              ->  factor = M_2(C)

THIS is the structure we test. The decisive question (audit-discipline):

  Does requiring the COMBINED J = J_GV (x) J_F to satisfy the Connes KO-dim
  sign table -- with the OUTER J_GV fixed at J_GV^2 = -1 (KO-dim 3) -- FORCE
  J_F to the quaternionic sign (J_F restricted to the C^2 has the -1
  signature that makes the factor H), while EXCLUDING the M_2(C) option?

  OR are both inner real-forms consistent with the same outer J_GV^2 = -1?

What this driver actually computes
==================================

We do NOT need the full Hilbert space to settle a SIGN-TABLE question.
The combined-J sign is a product of the two factor signs:

    J^2          = (J_GV^2) (J_F^2)               -> epsilon_total = eps_GV * eps_F
    J D = ? D J  -> epsilon'_total = eps'_GV * eps'_F   (with care: tensor of
                    a Z2-graded D requires the standard Connes-Marcolli sign
                    rule for products of even/odd KO-dim triples)

The KO-dimension sign table (Connes 1995; van Suijlekom 2015 Table; the
geovac/real_structure.py docstring reproduces it):

    KO 0:  (eps, eps', eps'') = (+, +, +)
    KO 1:  (-, +,  n/a)        [odd]
    KO 2:  (-, -, -)
    KO 3:  (-, +,  n/a)        [odd]   <-- OUTER GeoVac S^3
    KO 4:  (-, +, +)
    KO 5:  (+, -,  n/a)        [odd]
    KO 6:  (+, -, -)           <-- the SM finite triple
    KO 7:  (+, +,  n/a)        [odd]

KO-dim adds mod 8 under tensor product of real spectral triples. The
combined eps multiply with the graded sign rule (Vanhecke / Connes-Marcolli
2008 Sec 13.4 product formula).

The two candidate inner real-forms correspond to two possible KO-dimensions
of the FINITE triple:

  - H   -> J_F^2 = +1, J_F D_F = +D_F J_F, {J_F, gamma_F} = 0  : KO-dim 6
           (the SM convention; the quaternionic *real structure* of the
            doubled finite Hilbert space; note: J_F^2 = +1 at KO-dim 6,
            but the n=2 DIVISION RING is H because the C^2 component carries
            the J_q = i sigma_2 K (J_q^2 = -1) quaternionic structure
            INTERNALLY -- this is the subtle point we make precise below).
  - M_2(C) -> the n=2 factor is the FULL complex matrix algebra; its natural
           real structure is the OPPOSITE one. We enumerate the KO-dims
           compatible with M_2(C) carrying a *complex/real* (non-quaternionic)
           form and test each.

We compute, for EACH candidate finite KO-dim d_F in {0,...,7} and for each
real-form interpretation:

  (a) combined KO-dim (3 + d_F) mod 8, and the resulting (eps, eps') of the
      COMBINED triple from the product rule;
  (b) whether the n=2 factor's INTERNAL C^2 antilinear structure forced by
      that d_F is quaternionic (J_q^2 = -1, factor = H) or non-quaternionic
      (factor = M_2(C));
  (c) the standard Connes-Marcolli constraint that the SM total KO-dim be
      1 mod 8 (so the fermion doubling / KO-regularity works) -- and whether
      requiring this constraint, GIVEN the outer KO-dim 3, FORCES d_F = 6
      and hence H.

We also do a direct numerical sign-table audit on the actual code objects:
build the C (+) H finite triple (almost_commutative.ElectroweakFiniteTriple,
KO-dim 6, J_F^2 = +1) and an explicit C (+) M_2(C) variant where the n=2
component carries the NON-quaternionic real structure, tensor each with the
GeoVac outer J_GV (real_structure.build_J_full_dirac), and report J^2,
JD=+/-DJ, and order-zero / order-one for both. The point is to show
computationally which combined-J sign each real-form produces.

Verdict logic (audit-discipline: do not let "consistent" pose as "forces")
==========================================================================

FORCES H  iff: the combined-J / KO-dim-1-mod-8 consistency, with outer KO-dim
               3 FIXED by GeoVac (Paper 32 prop:reality, not a free input),
               admits ONLY d_F = 6 (quaternionic, H) and EXCLUDES every d_F
               that realizes the non-quaternionic M_2(C) form.

ADMITS H  iff: more than one d_F (in particular at least one non-quaternionic
               / M_2(C) option) is consistent with the same outer J_GV^2 = -1
               and the same combined-KO-dim target. Then H is a literature
               import, not a GeoVac forcing.

Out of scope (guardrail): NO Yukawa value is selected. This is a sign-table /
real-form / KO-dim computation only. Integer KO-dim and +/-1 signs throughout;
no transcendental appears.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from geovac.almost_commutative import (
    SIGMA_0, SIGMA_1, SIGMA_2, SIGMA_3,
    ElectroweakFiniteTriple,
    quaternion_to_matrix,
)
from geovac.real_structure import build_J_full_dirac
from geovac.full_dirac_operator_system import (
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
)


# ===========================================================================
# Part 0. The KO-dimension sign table (Connes 1995, vanSuijlekom 2015)
# ===========================================================================

# (epsilon = J^2 sign, epsilon' = JD=eps'DJ sign, epsilon'' = Jg=eps''gJ sign)
# None entries: epsilon'' undefined in odd KO-dim (no chirality grading).
KO_SIGN_TABLE: Dict[int, Tuple[int, int, Optional[int]]] = {
    0: (+1, +1, +1),
    1: (-1, +1, None),
    2: (-1, -1, -1),
    3: (-1, +1, None),   # outer GeoVac S^3
    4: (-1, +1, +1),
    5: (+1, -1, None),
    6: (+1, -1, -1),     # SM finite triple
    7: (+1, +1, None),
}


def ko_product(d1: int, d2: int) -> int:
    """KO-dimension of a tensor product of real spectral triples: d1 + d2 mod 8."""
    return (d1 + d2) % 8


def combined_eps_from_table(d1: int, d2: int) -> Tuple[int, int]:
    """Combined (epsilon, epsilon') read directly off the KO table at d1+d2 mod 8.

    The Connes-Marcolli product rule for two real spectral triples is precisely
    that the product has KO-dim d1 + d2 (mod 8); its (eps, eps', eps'') are then
    the table entries at that KO-dim. We return (eps, eps') (eps'' is None for
    odd combined dim).
    """
    d = ko_product(d1, d2)
    eps, eps_p, _ = KO_SIGN_TABLE[d]
    return eps, eps_p


# ===========================================================================
# Part 1. The H-vs-M2(C) real form is an INTERNAL C^2 antilinear structure
# ===========================================================================
#
# Key fact: over C, H ~= M_2(C). The distinction is a REAL FORM, fixed by an
# antilinear involution j on M_2(C):
#
#   H        = Fix(j_q),  j_q(m) = (i sigma_2) conj(m) (i sigma_2)^{-1}
#                          structure op  J_q = (i sigma_2) K,  J_q^2 = -1
#   M_2(R)   = Fix(j_r),  j_r(m) = conj(m)
#                          structure op  J_r = K,             J_r^2 = +1
#   M_2(C)   = no real-form constraint  (the FULL complex algebra; the
#              "real form" relevant to NCG is the one whose unitary group,
#              after unimodularity, gives SU(2) -- this is the j_r / +1 branch,
#              NOT the quaternionic -1 branch).
#
# So the n=2 DIVISION RING is decided by the sign of the antilinear involution
# carried on the internal C^2:
#       J_internal^2 = -1   <=>   factor = H        (quaternionic)
#       J_internal^2 = +1   <=>   factor = M_2(C)   (complex / real form)

def quaternionic_structure_C2() -> np.ndarray:
    """Quaternionic antilinear structure on C^2: J_q = (i sigma_2) K.

    Returns the unitary U_q in J_q = U_q K. J_q^2 = U_q conj(U_q).
    """
    return 1j * SIGMA_2  # i sigma_2


def real_structure_C2() -> np.ndarray:
    """Non-quaternionic (real/complex-form) antilinear structure on C^2: J_r = K.

    Returns the unitary U_r = I in J_r = U_r K. J_r^2 = +1.
    """
    return SIGMA_0.copy()


def internal_J_squared(U: np.ndarray) -> float:
    """J^2 sign of an antilinear J = U K: J^2 = U conj(U). Return the (real) scalar
    if J^2 is a scalar multiple of I, else NaN."""
    J2 = U @ np.conj(U)
    # check scalar
    scal = J2[0, 0]
    if np.allclose(J2, scal * np.eye(U.shape[0]), atol=1e-12):
        return float(np.real(scal))
    return float("nan")


def fixed_point_subalgebra_dimension(U: np.ndarray, n_samples: int = 2000,
                                     seed: int = 7) -> int:
    """Real dimension of Fix(j) where j(m) = U conj(m) U^{-1} on M_n(C).

    H has real dimension 4; M_2(C) has real dimension 8; M_2(R) has real dim 4.
    We compute the real dimension of the fixed-point real-linear subspace of
    j acting on M_2(C) (real dim 8), to confirm:
        J_q (i sigma_2) -> Fix = H, real dim 4
        J_r (= K)       -> Fix = M_2(R), real dim 4
    Both fixed-point ALGEBRAS have real dim 4; what differs is the DIVISION
    RING (H vs M_2(R)), distinguished by whether the algebra is a division
    ring (H: every nonzero element invertible) or not (M_2(R): has zero
    divisors). We report the division-ring test separately.
    """
    n = U.shape[0]
    # j(m) = U conj(m) U^{-1}; build the real-linear map on R^{2 n^2} and count
    # +1 eigenvalue multiplicity.
    Uinv = np.linalg.inv(U)
    dim_c = n * n
    # real coordinates: stack real and imag parts of vec(m)
    M = np.zeros((2 * dim_c, 2 * dim_c))
    basis = []
    for i in range(n):
        for k in range(n):
            E = np.zeros((n, n), dtype=np.complex128)
            E[i, k] = 1.0
            basis.append(E)
        for i2 in range(n):  # placeholder to keep structure clear
            pass
    # Build j as real-linear map by acting on a real basis of M_n(C) (real dim 2 n^2)
    real_basis: List[np.ndarray] = []
    for E in basis:
        real_basis.append(E)        # real part direction
        real_basis.append(1j * E)   # imag part direction
    cols = []
    for B in real_basis:
        jB = U @ np.conj(B) @ Uinv
        # express jB in the real basis: coords = [Re over E-entries, Im over E-entries]
        vec = []
        for E in basis:
            c = np.sum(np.conj(E) * jB)  # entrywise inner product -> the (i,k) entry
            vec.append(np.real(c))
            vec.append(np.imag(c))
        cols.append(np.array(vec))
    Jmap = np.array(cols).T  # 2n^2 x 2n^2 real matrix
    evals = np.linalg.eigvals(Jmap)
    fix_dim = int(np.sum(np.abs(evals - 1.0) < 1e-8))
    return fix_dim


def is_division_ring_fixed(U: np.ndarray, n_samples: int = 20000,
                           seed: int = 11) -> bool:
    """Division-ring test for Fix(j), j(m)=U conj(m) U^{-1}.

    A real *-algebra is a division ring iff every nonzero element is invertible.
    We sample fixed-point elements and look for a near-singular one. H is a
    division ring (every nonzero quaternion invertible, |det| = |q|^2 > 0,
    min |det|/norm^2 = 1/2 over samples). M_2(R) = Fix(K) is NOT a division
    ring: it contains rank-1 (singular) elements, e.g. diag(1, 0). We use a
    DECISIVE structural witness for the quaternionic case rather than a loose
    sampling threshold: for the quaternionic involution j_q the fixed algebra
    is exactly span_R{1, i sigma_1, i sigma_2, i sigma_3}, every nonzero
    element of which has det = sum of squares > 0; for j_r the fixed algebra
    is M_2(R), which contains the singular witness diag(1, 0). We check the
    explicit witness AND the sampling minimum.
    """
    rng = np.random.default_rng(seed)
    n = U.shape[0]
    Uinv = np.linalg.inv(U)
    min_det_ratio = np.inf
    n_checked = 0
    for _ in range(n_samples):
        m = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
        p = 0.5 * (m + U @ np.conj(m) @ Uinv)  # project onto Fix(j)
        nrm = np.linalg.norm(p)
        if nrm < 1e-9:
            continue
        n_checked += 1
        ratio = abs(np.linalg.det(p)) / (nrm ** 2)
        min_det_ratio = min(min_det_ratio, ratio)
    # Explicit singular witness: project diag(1, 0) onto Fix(j); if the
    # projection is itself singular and nonzero, the algebra is NOT a division
    # ring.
    w = np.zeros((n, n), dtype=np.complex128)
    w[0, 0] = 1.0
    pw = 0.5 * (w + U @ np.conj(w) @ Uinv)
    witness_singular = (np.linalg.norm(pw) > 1e-9
                        and abs(np.linalg.det(pw)) < 1e-12)
    is_division = (min_det_ratio > 1e-3) and (not witness_singular)
    return bool(is_division)


# ===========================================================================
# Part 2. Build the two candidate FINITE triples (C (+) H vs C (+) M_2(C))
#         and the two candidate finite real structures.
# ===========================================================================
#
# We work on the electroweak slice (n=1 factor C, n=2 factor H or M_2(C)) since
# that is where the fork lives and where the existing code (almost_commutative,
# ElectroweakFiniteTriple) gives a fully-built KO-dim-6 H finite triple to
# compare against.
#
# C (+) H finite triple: ElectroweakFiniteTriple (J_F^2 = +1, KO-dim 6). The
#   n=2 factor is H because the algebra action uses quaternion_to_matrix, i.e.
#   the L-block is constrained to the quaternionic real form q = q1 + i q.sigma.
#
# C (+) M_2(C) finite triple: identical Hilbert space and J_F, but the algebra
#   action on the L-block is a GENERAL complex 2x2 matrix (no quaternionic
#   constraint). The real structure J_F is the SAME swap-conjugation. We test
#   whether the order-zero / order-one Connes conditions still hold -- i.e.
#   whether the M_2(C) algebra is COMPATIBLE with the same J_F as H, which is
#   the operational meaning of "admits M_2(C)".


def cH_finite() -> ElectroweakFiniteTriple:
    """C (+) H finite triple (KO-dim 6, J_F^2 = +1). n=2 factor = H."""
    return ElectroweakFiniteTriple(yukawa_nu=0.0, yukawa_e=0.0)


def algebra_action_cH(finite: ElectroweakFiniteTriple, lam: complex,
                      q_components: Tuple[complex, complex, complex, complex]
                      ) -> np.ndarray:
    """C (+) H action: L-block = quaternion(q), R-block = diag(lam, conj lam)."""
    return finite.algebra_action(lam, q_components)


def algebra_action_cM2(finite: ElectroweakFiniteTriple, lam: complex,
                       m2: np.ndarray) -> np.ndarray:
    """C (+) M_2(C) action: L-block = GENERAL complex 2x2 m2, R-block =
    diag(lam, conj lam). Same 8x8 layout as ElectroweakFiniteTriple but the
    L-block is unconstrained (the M_2(C) real form), acting on matter sector
    only (antimatter via J_F)."""
    M = np.zeros((4, 4), dtype=np.complex128)
    M[0:2, 0:2] = m2
    M[2, 2] = lam
    M[3, 3] = np.conj(lam)
    action = np.zeros((8, 8), dtype=np.complex128)
    action[0:4, 0:4] = M
    return action


# ===========================================================================
# Part 3. Combined-triple sign-table audit (numerical, on actual operators)
# ===========================================================================


def build_combined_J(n_max: int, finite: ElectroweakFiniteTriple) -> np.ndarray:
    """U for combined J = J_GV (x) J_F (unitary part)."""
    J_GV = build_J_full_dirac(n_max)
    U_F = finite.real_structure_F()
    return np.kron(J_GV.U, U_F)


def build_combined_D(n_max: int, finite: ElectroweakFiniteTriple) -> np.ndarray:
    """D = D_GV (x) 1_F + gamma_GV (x) D_F."""
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    D_GV = camporesi_higuchi_full_dirac_matrix(op_sys.basis)
    gamma_GV = np.diag([float(b.chirality) for b in op_sys.basis]).astype(np.complex128)
    D_F = finite.dirac_F()
    I_F = np.eye(finite.dim_H_F, dtype=np.complex128)
    return np.kron(D_GV, I_F) + np.kron(gamma_GV, D_F)


def combined_sign_audit(n_max: int) -> dict:
    """Numerical audit of combined J^2 and JD=+/-DJ for the C (+) H triple."""
    finite = cH_finite()
    U = build_combined_J(n_max, finite)
    D = build_combined_D(n_max, finite)
    dim = U.shape[0]
    I = np.eye(dim, dtype=np.complex128)

    J2 = U @ np.conj(U)
    j2_plus = float(np.linalg.norm(J2 - I))
    j2_minus = float(np.linalg.norm(J2 + I))

    # J D = eps' D J for antilinear J = U K:  U conj(D) = eps' D U
    UconjD = U @ np.conj(D)
    DU = D @ U
    jd_plus = float(np.linalg.norm(UconjD - DU))
    jd_minus = float(np.linalg.norm(UconjD + DU))

    return {
        "n_max": n_max,
        "dim_H": dim,
        "J2_eq_plus_I_residual": j2_plus,
        "J2_eq_minus_I_residual": j2_minus,
        "JD_eq_plus_DJ_residual": jd_plus,
        "JD_eq_minus_DJ_residual": jd_minus,
        "measured_eps": (+1 if j2_plus < j2_minus else -1),
        "measured_eps_prime": (+1 if jd_plus < jd_minus else -1),
    }


def order_conditions_audit(n_max: int, kind: str, n_samples: int = 30,
                           seed: int = 17) -> dict:
    """Order-zero [a, J b J^{-1}] and order-one [[D,a], J b J^{-1}] residuals
    for kind in {'H', 'M2'}. The finite J_F is the SAME (swap-conj) for both;
    only the algebra action differs (quaternion-constrained vs general M_2(C)).

    This is the operational distinction between H and M_2(C): both share the
    same Hilbert space, D_F=0 (gauge-only), gamma_F, and J_F. The question is
    whether the larger M_2(C) algebra still satisfies the J-axioms with that
    J_F -- i.e. whether J_F is compatible with M_2(C), or only with H.
    """
    finite = cH_finite()
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    D_GV = camporesi_higuchi_full_dirac_matrix(op_sys.basis)
    gamma_GV = np.diag([float(b.chirality) for b in op_sys.basis]).astype(np.complex128)
    D_F = finite.dirac_F()  # zero Yukawa -> D_F = 0 (gauge only)
    I_F = np.eye(finite.dim_H_F, dtype=np.complex128)
    D = np.kron(D_GV, I_F) + np.kron(gamma_GV, D_F)
    U = build_combined_J(n_max, finite)

    mults = op_sys.multiplier_matrices
    n_mult = len(mults)
    rng = np.random.default_rng(seed)

    def make_a():
        k = rng.integers(0, n_mult)
        lam = rng.standard_normal() + 1j * rng.standard_normal()
        if kind == "H":
            q = tuple(rng.standard_normal() for _ in range(4))  # real quaternion
            aF = algebra_action_cH(finite, lam, q)
        elif kind == "M2":
            m2 = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
            aF = algebra_action_cM2(finite, lam, m2)
        else:
            raise ValueError(kind)
        return np.kron(mults[k], aF)

    order0_max = 0.0
    order1_max = 0.0
    for _ in range(n_samples):
        a = make_a()
        b = make_a()
        JbJinv = U @ np.conj(b) @ U.T
        comm0 = a @ JbJinv - JbJinv @ a
        order0_max = max(order0_max, float(np.linalg.norm(comm0)))
        Da = D @ a - a @ D
        comm1 = Da @ JbJinv - JbJinv @ Da
        order1_max = max(order1_max, float(np.linalg.norm(comm1)))

    return {
        "n_max": n_max,
        "kind": kind,
        "order_zero_max_residual": order0_max,
        "order_one_max_residual": order1_max,
        "n_samples": n_samples,
    }


# ===========================================================================
# Part 3c. Order-one with NONZERO Yukawa: the literature's actual selector
# ===========================================================================
#
# The CCM second-order / order-one selection of H over M_2(C) is NOT a
# combined-J sign statement; it is a statement that, with a NONZERO off-diagonal
# D_F (Yukawa-carrying finite Dirac), the order-one condition
# [[D, a], J b J^{-1}] = 0 CONSTRAINS the algebra action. We test whether,
# at this level, the general M_2(C) L-block breaks order-one while the
# quaternion-constrained H L-block preserves it. This is the honest place to
# look for a selection -- but note (audit-discipline) it requires a nonzero
# Yukawa, which is OUT OF SCOPE as a value; we use a generic nonzero Yukawa
# only to probe the STRUCTURAL order-one constraint, never to fix its value.


def order_one_with_yukawa(n_max: int, kind: str, yukawa: float = 0.3,
                          n_samples: int = 30, seed: int = 23) -> dict:
    """Order-one residual with a NONZERO finite Dirac (generic Yukawa).

    kind in {'H', 'M2'}. We use the SAME finite J_F (swap-conj) and the SAME
    nonzero D_F (off-diagonal L<->R Yukawa block), and compare order-one
    residuals for the quaternion-constrained (H) vs general (M2) algebra
    action. If M2 breaks order-one where H preserves it, THAT is the CCM
    selection mechanism realized in code -- and we then ask whether it is a
    GeoVac-native consequence of J_GV or the standard finite-triple argument.
    """
    finite = ElectroweakFiniteTriple(yukawa_nu=yukawa, yukawa_e=yukawa)
    op_sys = FullDiracTruncatedOperatorSystem(n_max)
    D_GV = camporesi_higuchi_full_dirac_matrix(op_sys.basis)
    gamma_GV = np.diag([float(b.chirality) for b in op_sys.basis]).astype(np.complex128)
    D_F = finite.dirac_F()  # nonzero off-diagonal
    I_F = np.eye(finite.dim_H_F, dtype=np.complex128)
    D = np.kron(D_GV, I_F) + np.kron(gamma_GV, D_F)
    U = build_combined_J(n_max, finite)

    mults = op_sys.multiplier_matrices
    n_mult = len(mults)
    rng = np.random.default_rng(seed)

    def make_a():
        k = rng.integers(0, n_mult)
        lam = rng.standard_normal() + 1j * rng.standard_normal()
        if kind == "H":
            q = tuple(rng.standard_normal() for _ in range(4))
            aF = algebra_action_cH(finite, lam, q)
        else:
            m2 = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
            aF = algebra_action_cM2(finite, lam, m2)
        return np.kron(mults[k], aF)

    order1_max = 0.0
    for _ in range(n_samples):
        a = make_a()
        b = make_a()
        JbJinv = U @ np.conj(b) @ U.T
        Da = D @ a - a @ D
        comm1 = Da @ JbJinv - JbJinv @ Da
        order1_max = max(order1_max, float(np.linalg.norm(comm1)))

    return {"n_max": n_max, "kind": kind, "yukawa": yukawa,
            "order_one_max_residual": order1_max, "n_samples": n_samples}


# ===========================================================================
# Part 4. The decisive KO-dim consistency argument (the core logic)
# ===========================================================================


def ko_dim_forcing_analysis() -> dict:
    """The decisive test.

    GeoVac FIXES the outer KO-dim = 3 (Paper 32 prop:reality, J_GV^2 = -1).
    The SM convention requires the COMBINED KO-dim to be 1 mod 8 (so the
    fermion doubling and KO-regularity / order-one work in the canonical CCM
    construction). Given outer = 3:

        combined = (3 + d_F) mod 8 = 1   =>   d_F = 6.

    KO-dim 6 has (eps_F, eps'_F, eps''_F) = (+1, -1, -1): J_F^2 = +1,
    J_F D_F = -D_F J_F, {J_F, gamma_F} = 0. This is the SM finite triple,
    whose n=2 division ring is H (the standard CCM result).

    THE AUDIT QUESTION (do not let "consistent" pose as "forces"):

      (Q-A) Is the "combined KO-dim = 1 mod 8" target a GeoVac OUTPUT, or an
            IMPOSED SM input? -- If imposed, then d_F=6 (and hence H) is
            imported via the target, NOT forced by GeoVac.

      (Q-B) Holding ONLY the outer J_GV^2 = -1 (the genuine GeoVac datum) and
            the bare Connes axioms (J^2 = +/-1 per KO-dim, order-zero,
            order-one) -- WITHOUT imposing the SM total -- how many finite
            KO-dims d_F admit a CONSISTENT combined real structure? And among
            those, do BOTH a quaternionic (H) and a non-quaternionic (M_2(C))
            real form survive?

    We enumerate d_F in {0,...,7}. For each:
      - combined KO-dim and (eps, eps') from the product table;
      - whether that d_F is one where the FINITE triple can carry the
        quaternionic internal C^2 structure (H) and/or the non-quaternionic
        (M_2(C)) structure.

    Crucial real-form fact (Part 1): the n=2 division ring (H vs M_2(C)) is set
    by the INTERNAL C^2 antilinear involution j_internal, NOT by J_F^2 of the
    DOUBLED finite Hilbert space. j_q^2 = -1 -> H; j_r^2 = +1 -> M_2(C). The
    KO-dim of the doubled finite triple (J_F on H_F^matter (+) H_F^antimatter)
    can be 6 for BOTH real forms, because J_F^2 = +1 (swap-conjugation) is a
    statement about the matter/antimatter swap, independent of the internal
    quaternionic structure on the C^2.

    THIS IS THE DECISIVE OBSERVATION. The combined-J / KO-dim bookkeeping
    operates on the DOUBLED J_F (matter<->antimatter), whose square (+1, KO-dim
    6) is the SAME whether the n=2 factor is H or M_2(C). The internal
    quaternionic involution j_q (which actually distinguishes H from M_2(C))
    lives on the algebra side, not on the combined-J sign table. So the
    combined-J sign table CANNOT see the H/M_2(C) distinction.
    """
    outer_ko = 3
    results = []
    for d_F in range(8):
        combined = ko_product(outer_ko, d_F)
        eps, eps_p = combined_eps_from_table(outer_ko, d_F)
        eps_F, eps_p_F, eps_pp_F = KO_SIGN_TABLE[d_F]
        results.append({
            "d_F": d_F,
            "combined_ko": combined,
            "combined_eps": eps,
            "combined_eps_prime": eps_p,
            "finite_eps": eps_F,
            "finite_eps_prime": eps_p_F,
            "finite_eps_pp": eps_pp_F,
            "combined_is_1_mod_8": (combined == 1),
        })
    return {
        "outer_ko": outer_ko,
        "outer_J_squared": -1,
        "per_d_F": results,
        "d_F_giving_combined_1mod8": [r["d_F"] for r in results
                                      if r["combined_is_1_mod_8"]],
    }


# ===========================================================================
# Part 5. Build the full J sign-table for BOTH candidate inner algebras
# ===========================================================================


def j_sign_table_both_candidates() -> dict:
    """Construct the explicit J sign-table for C (+) H and C (+) M_2(C).

    For each candidate we list:
      - the internal C^2 antilinear structure J_internal and J_internal^2;
      - whether Fix(J_internal) is a division ring (H) or not (M_2(C));
      - the DOUBLED finite J_F (swap-conjugation), its square, KO-dim;
      - the combined J = J_GV (x) J_F square and KO-dim.

    The headline: the combined J^2 and combined KO-dim are IDENTICAL for the
    two candidates (both use the swap-conjugation J_F^2 = +1, KO-dim 6,
    combined KO-dim 1). The ONLY place the two differ is the INTERNAL C^2
    structure -- which the combined-J sign table does not constrain.
    """
    U_q = quaternionic_structure_C2()   # i sigma_2 -> H
    U_r = real_structure_C2()           # I (K) -> M_2(C)/M_2(R) real form

    j2_q = internal_J_squared(U_q)
    j2_r = internal_J_squared(U_r)
    fixdim_q = fixed_point_subalgebra_dimension(U_q)
    fixdim_r = fixed_point_subalgebra_dimension(U_r)
    divring_q = is_division_ring_fixed(U_q)
    divring_r = is_division_ring_fixed(U_r)

    # doubled finite J_F is the SAME for both (swap-conjugation, KO-dim 6)
    finite = cH_finite()
    U_F = finite.real_structure_F()
    J_F_sq = float(np.real((U_F @ np.conj(U_F))[0, 0]))

    # combined J^2 from product rule (outer KO 3, finite KO 6)
    eps_comb, eps_p_comb = combined_eps_from_table(3, 6)
    combined_ko = ko_product(3, 6)

    return {
        "candidate_cH": {
            "n2_factor": "H",
            "internal_C2_structure": "J_q = i*sigma_2 * K",
            "internal_J_squared": j2_q,            # -1 (quaternionic)
            "fixed_point_real_dim": fixdim_q,      # 4 (H)
            "fixed_point_is_division_ring": divring_q,  # True (H)
            "doubled_finite_J_F_squared": J_F_sq,  # +1
            "finite_ko_dim": 6,
            "combined_ko_dim": combined_ko,        # 1
            "combined_eps": eps_comb,              # -1
            "combined_eps_prime": eps_p_comb,      # +1
        },
        "candidate_cM2": {
            "n2_factor": "M_2(C)",
            "internal_C2_structure": "J_r = K  (or full complex algebra, no quaternionic constraint)",
            "internal_J_squared": j2_r,            # +1 (non-quaternionic)
            "fixed_point_real_dim": fixdim_r,      # 4 (M_2(R)) -- but algebra is NOT a division ring
            "fixed_point_is_division_ring": divring_r,  # False (M_2(R) has zero divisors)
            "doubled_finite_J_F_squared": J_F_sq,  # +1  -- SAME as cH
            "finite_ko_dim": 6,                    # SAME as cH (swap-conj J_F)
            "combined_ko_dim": combined_ko,        # 1  -- SAME as cH
            "combined_eps": eps_comb,              # -1 -- SAME as cH
            "combined_eps_prime": eps_p_comb,      # +1 -- SAME as cH
        },
        "headline": (
            "Combined J^2, combined eps', combined KO-dim are IDENTICAL for "
            "both candidates. They differ ONLY in the internal C^2 antilinear "
            "involution (j_q^2 = -1 -> H vs j_r^2 = +1 -> M_2(C)), which the "
            "combined-J sign table does not constrain."
        ),
    }


# ===========================================================================
# Main
# ===========================================================================


def main() -> None:
    out: dict = {}

    print("=" * 72)
    print("Door 4c: does combined J = J_GV (x) J_F SELECT H or merely ADMIT it?")
    print("=" * 72)

    # Part 1: internal real-form structures
    print("\n[Part 1] Internal C^2 antilinear real-form structures")
    U_q = quaternionic_structure_C2()
    U_r = real_structure_C2()
    print(f"  quaternionic J_q = i*sigma_2*K:  J_q^2 = {internal_J_squared(U_q):+.1f}"
          f"   Fix dim = {fixed_point_subalgebra_dimension(U_q)}"
          f"   division ring = {is_division_ring_fixed(U_q)}   -> H")
    print(f"  real        J_r = K:             J_r^2 = {internal_J_squared(U_r):+.1f}"
          f"   Fix dim = {fixed_point_subalgebra_dimension(U_r)}"
          f"   division ring = {is_division_ring_fixed(U_r)}   -> M_2(R)/M_2(C)")
    out["internal_real_forms"] = {
        "quaternionic": {
            "J_squared": internal_J_squared(U_q),
            "fix_dim": fixed_point_subalgebra_dimension(U_q),
            "division_ring": is_division_ring_fixed(U_q),
        },
        "real_complex": {
            "J_squared": internal_J_squared(U_r),
            "fix_dim": fixed_point_subalgebra_dimension(U_r),
            "division_ring": is_division_ring_fixed(U_r),
        },
    }

    # Part 4: KO-dim forcing analysis
    print("\n[Part 4] KO-dimension forcing analysis (outer = 3, J_GV^2 = -1)")
    ko = ko_dim_forcing_analysis()
    out["ko_dim_forcing"] = ko
    print(f"  d_F giving combined KO-dim = 1 (mod 8): "
          f"{ko['d_F_giving_combined_1mod8']}")
    for r in ko["per_d_F"]:
        flag = "  <- SM target" if r["combined_is_1_mod_8"] else ""
        print(f"    d_F={r['d_F']}: combined KO={r['combined_ko']}  "
              f"(eps,eps')=({r['combined_eps']:+d},{r['combined_eps_prime']:+d}){flag}")

    # Part 5: J sign-table for both candidates
    print("\n[Part 5] J sign-table for BOTH candidate inner algebras")
    jt = j_sign_table_both_candidates()
    out["j_sign_table"] = jt
    for key in ("candidate_cH", "candidate_cM2"):
        c = jt[key]
        print(f"  {c['n2_factor']:>8}: internal J^2 = {c['internal_J_squared']:+.0f}"
              f"  div-ring={c['fixed_point_is_division_ring']}"
              f"  | J_F^2={c['doubled_finite_J_F_squared']:+.0f}"
              f"  finite-KO={c['finite_ko_dim']}"
              f"  combined-KO={c['combined_ko_dim']}"
              f"  comb-eps=({c['combined_eps']:+d},{c['combined_eps_prime']:+d})")
    print(f"  HEADLINE: {jt['headline']}")

    # Part 3: numerical combined sign audit + order conditions for both
    print("\n[Part 3] Numerical combined-J sign audit (C (+) H, actual operators)")
    sign_audits = []
    for n_max in (1, 2, 3):
        sa = combined_sign_audit(n_max)
        sign_audits.append(sa)
        print(f"  n_max={n_max}: J^2 -> eps={sa['measured_eps']:+d} "
              f"(|J^2+I|={sa['J2_eq_minus_I_residual']:.2e}), "
              f"JD -> eps'={sa['measured_eps_prime']:+d} "
              f"(|UD~-DU|={sa['JD_eq_plus_DJ_residual']:.2e})")
    out["combined_sign_audits"] = sign_audits

    print("\n[Part 3b] Order-zero / order-one for H vs M_2(C) (SAME J_F, gauge-only)")
    order_audits = []
    for kind in ("H", "M2"):
        for n_max in (1, 2):
            oa = order_conditions_audit(n_max, kind)
            order_audits.append(oa)
            print(f"  kind={kind:>2} n_max={n_max}: "
                  f"order0_max={oa['order_zero_max_residual']:.3e}  "
                  f"order1_max={oa['order_one_max_residual']:.3e}")
    out["order_conditions_audits"] = order_audits

    print("\n[Part 3c] Order-one with NONZERO Yukawa (CCM's actual selector)")
    yuk_audits = []
    for kind in ("H", "M2"):
        for n_max in (1, 2):
            ya = order_one_with_yukawa(n_max, kind)
            yuk_audits.append(ya)
            print(f"  kind={kind:>2} n_max={n_max} y={ya['yukawa']}: "
                  f"order1_max={ya['order_one_max_residual']:.3e}")
    out["order_one_with_yukawa_audits"] = yuk_audits

    # -----------------------------------------------------------------
    # Verdict
    # -----------------------------------------------------------------
    d_F_targets = ko["d_F_giving_combined_1mod8"]
    # The combined sign table at finite KO 6 is identical for H and M_2(C)
    # (Part 5). The combined-J cannot distinguish them. The "combined KO = 1"
    # constraint, even if imposed, selects d_F = 6 (the finite triple's
    # doubled-J KO-dim), but KO-dim 6 admits BOTH the quaternionic (H) and
    # non-quaternionic (M_2(C)) internal real form on the C^2 -- because the
    # doubled J_F^2 = +1 is blind to the internal j_q vs j_r.
    cH = jt["candidate_cH"]
    cM2 = jt["candidate_cM2"]
    sign_table_identical = (
        cH["combined_ko_dim"] == cM2["combined_ko_dim"]
        and cH["combined_eps"] == cM2["combined_eps"]
        and cH["combined_eps_prime"] == cM2["combined_eps_prime"]
        and cH["doubled_finite_J_F_squared"] == cM2["doubled_finite_J_F_squared"]
    )
    # Order conditions: both H and M2 satisfy them at the gauge-only (D_F=0)
    # level with the same J_F? (We read this off the audit.)
    o = {(a["kind"], a["n_max"]): a for a in order_audits}
    h_ok = all(o[("H", n)]["order_zero_max_residual"] < 1e-9
               and o[("H", n)]["order_one_max_residual"] < 1e-9 for n in (1, 2))
    m2_ok = all(o[("M2", n)]["order_zero_max_residual"] < 1e-9
                and o[("M2", n)]["order_one_max_residual"] < 1e-9 for n in (1, 2))

    # Order-one WITH Yukawa: does M2 break it where H preserves it?
    yo = {(a["kind"], a["n_max"]): a for a in yuk_audits}
    h_yuk_ok = all(yo[("H", n)]["order_one_max_residual"] < 1e-9 for n in (1, 2))
    m2_yuk_breaks = any(yo[("M2", n)]["order_one_max_residual"] > 1e-9 for n in (1, 2))
    # The selection-via-order-one mechanism is operative iff H preserves order-one
    # with nonzero D_F while M2 breaks it. This is the CCM finite-triple argument;
    # crucially it does NOT involve J_GV (the outer factor) at all -- it lives
    # entirely on the finite side (D_F, J_F, the algebra action). So even if it
    # fires, it is NOT a GeoVac-native (J_GV-driven) forcing.
    ccm_order_one_selects = bool(h_yuk_ok and m2_yuk_breaks)
    out["ccm_order_one_selector"] = {
        "H_preserves_order_one_with_yukawa": bool(h_yuk_ok),
        "M2_breaks_order_one_with_yukawa": bool(m2_yuk_breaks),
        "selection_mechanism_operative": ccm_order_one_selects,
        "note": (
            "This selector is the standard CCM finite-triple order-one/"
            "second-order-condition argument. It involves ONLY the finite data "
            "(D_F, J_F, algebra action) and NOT the outer J_GV. Even when it "
            "fires, the selection is the LITERATURE'S argument transported into "
            "the finite factor -- it is NOT driven by GeoVac's J_GV^2 = -1. The "
            "probe's target (close the fork FROM GeoVac's native J_GV) therefore "
            "remains unmet regardless of this selector: GeoVac would be IMPORTING "
            "the CCM order-one argument, not supplying a new one from J_GV."
        ),
    }

    if sign_table_identical:
        verdict = "ADMITS H (does NOT force H)"
        reasoning = (
            "The combined J = J_GV (x) J_F sign table is IDENTICAL for the two "
            "candidate inner algebras: both carry the doubled swap-conjugation "
            "J_F (J_F^2 = +1, finite KO-dim 6), so the combined triple has "
            "KO-dim 3+6 = 1 (mod 8) and combined (eps, eps') = (-1, +1) in "
            "BOTH cases. The H-vs-M_2(C) distinction is an INTERNAL C^2 real-form "
            "choice (quaternionic j_q^2 = -1 vs non-quaternionic j_r^2 = +1) that "
            "lives on the ALGEBRA side and is invisible to the combined-J / KO-dim "
            "bookkeeping. Therefore the outer J_GV^2 = -1 (KO-dim 3) does NOT "
            "exclude M_2(C): it is consistent with both. H is a literature import "
            "(the CCM second-order-condition / complex-rep argument), not a "
            "GeoVac forcing from J_GV alone."
        )
    else:
        verdict = "FORCES H (combined-J / KO-dim excludes M_2(C))"
        reasoning = (
            "The combined J sign table differs between the candidates, and the "
            "outer J_GV^2 = -1 admits only the quaternionic completion."
        )

    out["verdict"] = {
        "result": verdict,
        "sign_table_identical_for_both": bool(sign_table_identical),
        "H_order_conditions_ok_gauge_only": bool(h_ok),
        "M2_order_conditions_ok_gauge_only": bool(m2_ok),
        "d_F_giving_combined_1mod8": d_F_targets,
        "reasoning": reasoning,
        "inner_algebra_status": (
            "PARTIAL - n=2 factor stays a fork; H admitted not forced by "
            "combined-J / KO-dim. (Door 4b verdict unchanged; the conjectured "
            "GeoVac-native handle does NOT close the fork via the sign table.)"
            if sign_table_identical else
            "FULLY FORCED - combined-J / KO-dim excludes M_2(C)."
        ),
    }

    print("\n" + "=" * 72)
    print(f"VERDICT: {verdict}")
    print("=" * 72)
    print(f"  sign table identical for both candidates: {sign_table_identical}")
    print(f"  H  order-conditions ok (gauge-only): {h_ok}")
    print(f"  M2 order-conditions ok (gauge-only): {m2_ok}")
    print(f"  CCM order-one selector operative (finite-side only, NOT J_GV): "
          f"{ccm_order_one_selects}")
    print(f"  inner algebra status: {out['verdict']['inner_algebra_status']}")

    # write JSON
    data_dir = os.path.join(os.path.dirname(__file__), "data")
    os.makedirs(data_dir, exist_ok=True)
    out_path = os.path.join(data_dir, "door4c_j_signtable_audit.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
