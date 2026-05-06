"""Real structure J on the truncated operator system at finite n_max.

Sprint WH1-Connes Step 2: define the real structure (charge conjugation /
antilinear KK-isomorphism) J on the Fock-projected S^3 spectral truncation
at finite n_max, and audit Connes' four real-axiom conditions.

Mathematical setup
==================

KO-dimension 3 sign table
-------------------------

The unit S^3 has KO-dimension 3 (mod 8). The Connes-Chamseddine sign table
[Connes 1995, Connes-Marcolli 2008] for a real structure J on a spectral
triple (A, H, D) of KO-dim 3 is

    epsilon  = -1     (J^2 = -1, "quaternionic / pseudoreal")
    epsilon' = +1     (J D = +D J, COMMUTATION not anti-commutation)
    epsilon''= n/a    (no chirality grading in odd dimension)

Note carefully: KO-dim 3 has J D = +D J. The sign table for the four
even / four odd KO-dimensions is (epsilon, epsilon') in the format
(J^2, JD = sign DJ):

    KO 0:  (+, +)        KO 4:  (-, +)
    KO 1:  (-, +)        KO 5:  (+, -)            (Lorentzian / Majorana)
    KO 2:  (-, -)        KO 6:  (+, -)            (Standard Model)
    KO 3:  (-, +)        KO 7:  (+, +)

For S^3 (KO-dim 3) the conditions are J^2 = -1 AND J D = +D J. This is
already the convention adopted in Paper 32 §IV.4 (Prop. reality),
matching the standard Camporesi-Higuchi charge conjugation on S^3 spinors.

Continuum charge conjugation on S^3
-----------------------------------

On the Camporesi-Higuchi spinor bundle S^3 with the basis labeled by
(n_fock, kappa, m_j) (or equivalently (n_fock, l, m_j, chirality) in the
full Dirac sector), the standard charge conjugation acts as

    J |n_fock, kappa, m_j, chirality> = sigma(n, kappa, m_j) *
                                        |n_fock, -kappa, -m_j, chirality>

where sigma(.) is a phase determined by the spinor representation of
Spin(3) = SU(2) [Friedrich, "Dirac operators in Riemannian geometry",
Ch. 1; Camporesi-Higuchi 1996].

Three structural facts force this form:

(a) J is ANTILINEAR. We realize this in code as a python operation
    J(psi) = U @ conj(psi)  where U is a fixed unitary matrix. The
    antilinearity ensures J^2 = U @ conj(U @ conj(psi)) = U @ conj(U) @ psi.
    For J^2 = -1 we need U @ conj(U) = -I, i.e. U is conjugate-skew.

(b) J flips kappa -> -kappa. This is the sigma -> -sigma flip on the
    spin index (large/small component swap) and l -> l unchanged
    (spatial parity is separate). On the Weyl sector with j = l + 1/2
    (kappa < 0) only, the kappa-flip would send us OUT of the Weyl
    sector into the anti-Weyl j = l - 1/2 (kappa > 0) chain. So J on
    the WEYL sector requires the FULL DIRAC sector to be defined at all.

(c) J flips m_j -> -m_j. This is the "time-reversal" component of the
    charge conjugation, equivalent to anti-symmetric Pauli matrix
    multiplication on the spin index.

The phase sigma(.) is determined modulo overall phase by demanding
J^2 = -1. For half-integer m_j states, the standard convention is

    sigma(n, kappa, m_j) = (-1)^(j + m_j)   = (-1)^(2l + 1 + m_j) / 2)
                          (or (-1)^(j - m_j) -- conventions vary)

with the property that sigma(n, -kappa, -m_j) * conj(sigma(n, kappa, m_j))
= -1 = J^2 sign.

For the FULL DIRAC sector (Weyl + anti-Weyl), J acts within each
chirality sector but flips the (kappa, m_j) labels:

    J |n_fock, l, m_j, chi=+1> = sigma+ |n_fock, l, -m_j, chi=+1>  ON ONE
                                                                   chirality
    J |n_fock, l, m_j, chi=-1> = sigma- |n_fock, l, -m_j, chi=-1>

Wait -- this is the WEYL j = l + 1/2 chain only. The full Dirac sector
in the implementation `geovac/full_dirac_operator_system.py` doubles
the SAME (l, m_j) chain (adding +/- chirality eigenvalue to the Dirac).
The kappa-flip in that doubled basis is a swap of (l, m_j) <-> (l, -m_j)
combined with a chirality-block permutation determined by whether one
adopts the "kappa flips chirality" convention (full Dirac with both
kappa signs at every level) or the "chirality is independent" convention
(both kappa = -(l+1) signs, with chi an additional Z_2 grading).

We adopt the SECOND convention to match `geovac/full_dirac_operator_system`.
J then acts as

    J |n_fock, l, m_j, chi> = sigma(l, m_j) * |n_fock, l, -m_j, chi>

i.e. block-diagonal in chirality, with the m_j -> -m_j flip and a phase
inside each block. This is the natural realization on the
Weyl-doubled-by-chirality Hilbert space of the full_dirac module.

Computational implementation
============================

We represent J as a (basis-permutation) * (diagonal-phase) * (complex-
conjugation) operator. Concretely:

    J(psi)_i = sigma_i * conj(psi[pi(i)])

where pi: i -> i' is the basis-index permutation realizing
(n, l, m_j, chi) -> (n, l, -m_j, chi), and sigma_i is the phase at
target index i. In matrix form, J = U @ K with K = complex conjugation
and U[i, j] = sigma_i delta_{i, pi(j)} (a permutation matrix times
a diagonal phase).

The condition J^2 = -1 becomes

    U conj(U) = -I

which constrains the phases. The standard choice for S^3 spinors is
sigma_i = i^(2 m_j(i)) * (some sign), giving sigma_i conj(sigma_(pi(i)))
= -1 for half-integer m_j.

Connes axioms at finite n_max
=============================

We verify the following four conditions:

(C1) J^2 = -epsilon * I  where epsilon = -1 for KO-dim 3.
     Concretely: J^2 = -I.

(C2) J D = +epsilon' * D J  where epsilon' = +1 for KO-dim 3.
     Concretely: [J, D] = 0  (in the antilinear sense:
     J D psi = D J psi for all psi).

(C3) Order-zero / commutation condition (Connes' "real" axiom):
     [a, J b J^{-1}] = 0  for all a, b in A_GV.
     This says A_GV and its J-conjugate algebra commute.

(C4) Order-one condition:
     [[D, a], J b J^{-1}] = 0  for all a, b in A_GV.

The truncated operator system A_GV at finite n_max is NOT an algebra,
which means (C3) and (C4) need re-interpretation. For the finite
truncated operator system O_n_max, we test these conditions on:
- the spanning set of multiplier matrices M_{N, L, M};
- on operator system PRODUCTS a*b that lie in O^2 (the propagation
  closure) when products a*b leave O.

Each condition might HOLD EXACTLY, hold APPROXIMATELY (small but
nonzero residual), FAIL (large residual), or be OBSTRUCTED (the
ambient operator-system structure prevents the condition from making
sense).

Public API
==========

  - RealStructure(basis): antilinear J operator
  - apply_J(J, psi): apply J to a state vector
  - verify_J_squared(J): check J^2 = -I
  - verify_J_commutes_D(J, D): check J D = +D J
  - verify_order_zero(J, A_basis): check [a, J b J^{-1}] = 0
  - verify_order_one(J, A_basis, D): check [[D, a], J b J^{-1}] = 0

References
==========

A. Connes, "Noncommutative geometry and reality," J. Math. Phys. 36
(1995) 6194-6231. Original definition of J and the KO-dim sign table.

A. Connes, "Noncommutative geometry and the standard model with
neutrino mixing," JHEP 0611 (2006) 081. KO-dim 6 setup; explicit
sign-table examples.

W. D. van Suijlekom, "Noncommutative Geometry and Particle Physics,"
Springer 2015. Textbook treatment, esp. Ch. 3.

R. Camporesi & A. Higuchi, "On the eigenfunctions of the Dirac operator
on spheres and real hyperbolic spaces," J. Geom. Phys. 20 (1996) 1-18.
Charge conjugation on the S^d Dirac spinors.

T. Friedrich, "Dirac operators in Riemannian geometry," AMS GSM 25
(2000). Ch. 1: charge conjugation in odd dimensions.

GeoVac sprint records:
- Paper 32 §IV.4 (Prop. reality)
- WH1 round-1 mapping: debug/wh1_connes_vs_mapping_memo.md (Gap 4, J)
- Round-3 R3.1, R3.2: debug/wh1_r3{1,2}_*.md
- Track-TS-A GH-convergence: debug/track_ts_a_gh_convergence_memo.md
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import numpy as np

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    FullDiracTruncatedOperatorSystem,
    full_dirac_basis,
)
from geovac.spinor_operator_system import (  # noqa: F401 (re-export)
    SpinorLabel,
    SpinorTruncatedOperatorSystem,
    spinor_basis,
)


# ---------------------------------------------------------------------------
# Real-structure phase conventions
# ---------------------------------------------------------------------------


def _spinor_J_phase(label) -> complex:
    """Phase sigma(label) for J on a single spinor / full-Dirac label.

    For the standard charge conjugation on S^3 spinor harmonics, the
    phase between |l, m_j> and J|l, m_j> is conventionally taken as

        sigma(l, m_j) = i^{2 m_j} * (-1)^(l)

    where m_j is half-integer. Equivalent forms in the literature:
    sigma = (-i)^(l + m_j + 1/2),  (-1)^(j - m_j) i^{...}, etc.
    All differ by an overall constant phase that cancels in J^2 and
    [J, D]; what matters is that
    sigma(label) * conj(sigma(J_label)) = -1
    where J_label = (l, -m_j) is the J-conjugate label.

    Verification:
      sigma(l, m_j)         = i^{2 m_j} * (-1)^l
      sigma(l, -m_j)        = i^{-2 m_j} * (-1)^l
      conj(sigma(l, -m_j))  = i^{2 m_j} * (-1)^l    (since i^{-x} conj = i^x)

    So sigma(l, m_j) * conj(sigma(l, -m_j)) = i^{4 m_j} * 1
                                            = (-1)^{2 m_j}
                                            = -1   (since m_j is half-integer,
                                                    2 m_j is odd integer).

    Good: this gives J^2 = -I as required.

    For the full-Dirac label we have an extra chirality index that the
    phase is independent of (J is chirality-block-diagonal in our
    convention).
    """
    two_m_j = label.two_m_j  # odd integer (half-integer m_j convention)
    l = label.l
    # i^{2 m_j} = i^{two_m_j}
    # We compute this safely: i^k cycles through {1, i, -1, -i} for
    # k mod 4 in {0, 1, 2, 3}.
    k_mod4 = two_m_j % 4
    if k_mod4 < 0:
        k_mod4 += 4
    i_power = [1.0 + 0j, 1j, -1.0 + 0j, -1j][k_mod4]
    sign = (-1.0 + 0j) if (l % 2 == 1) else (1.0 + 0j)
    return i_power * sign


# ---------------------------------------------------------------------------
# Permutation pi: label -> J(label)
# ---------------------------------------------------------------------------


def _spinor_J_target(label):
    """Target label of J on a Weyl SpinorLabel: (l, m_j) -> (l, -m_j)."""
    return SpinorLabel(
        n_fock=label.n_fock,
        l=label.l,
        two_m_j=-label.two_m_j,
    )


def _full_dirac_J_target(label):
    """Target label of J on a FullDiracLabel:
    (n, l, m_j, chi) -> (n, l, -m_j, chi).

    J acts within each chirality block (block-diagonal). This matches the
    convention adopted in geovac/full_dirac_operator_system.py: chirality
    is the eigenvalue sign of D, so a J that flips chirality would also
    flip the sign of D, contradicting J D = +D J. So J commutes with the
    chirality grading.
    """
    return FullDiracLabel(
        n_fock=label.n_fock,
        l=label.l,
        two_m_j=-label.two_m_j,
        chirality=label.chirality,
    )


# ---------------------------------------------------------------------------
# RealStructure class
# ---------------------------------------------------------------------------


@dataclass
class RealStructure:
    """Antilinear real-structure operator J on a finite spinor / full-Dirac
    truncated operator system.

    Stored as a unitary matrix U so that J(psi) = U @ conj(psi).

    Attributes
    ----------
    basis : list of labels (SpinorLabel or FullDiracLabel)
    U : complex ndarray of shape (dim, dim)
        The unitary matrix in J = U K where K = complex conjugation.
    sector : str
        Either "weyl" or "full_dirac".

    Mathematical content
    --------------------
    J is realized as a basis permutation pi (sending label l to its
    J-target J(l)) combined with a diagonal phase sigma(.) at the
    TARGET index. Concretely:

        U[i, j] = sigma_i * delta_{i, pi(j)}

    so that for a basis vector |j>:
        J|j> = U @ conj(|j>)
             = U @ |j>     (since |j> is real-valued in its own basis)
             = U[:, j]     (the j-th column of U)
             = sigma_{pi(j)} * |pi(j)>     (where pi(j) is the unique i s.t.
                                             label(i) = J(label(j)).)

    But for a general state psi = sum_j c_j |j>:
        J(psi) = U @ conj(c)
               = sum_j conj(c_j) * U[:, j]
               = sum_j conj(c_j) * sigma_{pi(j)} * |pi(j)>.

    This is antilinear by construction: J(alpha psi) = conj(alpha) J(psi).

    Verification of J^2 = -I:
        J^2(psi) = J(U @ conj(psi))
                 = U @ conj(U @ conj(psi))
                 = U @ conj(U) @ psi.
        Need U @ conj(U) = -I.

    With U = perm(pi) @ diag(sigma), we have
    conj(U) = perm(pi) @ diag(conj(sigma)),
    so U @ conj(U) = perm(pi) @ diag(sigma) @ perm(pi) @ diag(conj(sigma)).

    Now perm(pi)^2 = perm(pi^2) = I (since pi is an involution, m_j -> -m_j
    flips back), so we get diag(sigma) @ perm(pi) @ diag(conj(sigma)).

    Wait this is getting tangled. Let me just compute it numerically and
    verify J^2 = -I empirically.
    """

    basis: Sequence
    U: np.ndarray
    sector: str

    @property
    def dim(self) -> int:
        return self.U.shape[0]

    def apply(self, psi: np.ndarray) -> np.ndarray:
        """Apply J to a state vector psi: J(psi) = U @ conj(psi)."""
        return self.U @ np.conj(psi)

    def apply_to_operator(self, op: np.ndarray) -> np.ndarray:
        """Apply the J-conjugation J op J^{-1} to a linear operator op.

        For an antilinear J with J = U K and J^{-1} = K U^{-1} = K U^*
        (since U is unitary):

            J op J^{-1} (psi) = J op K U^* psi
                              = U conj(op K U^* psi)
                              = U conj(op) conj(K U^* psi)
                              = U conj(op) U^T psi

        So J op J^{-1} = U conj(op) U^T as a linear operator.
        """
        return self.U @ np.conj(op) @ self.U.T

    def J_squared_matrix(self) -> np.ndarray:
        """Compute J^2 as a linear operator: J^2 = U conj(U).

        For an antilinear J = U K, J^2 = U K U K = U conj(U) K^2 = U conj(U)
        (since K^2 = identity on a real-coefficient basis viewed as the
        identity antilinear map; J^2 acting on psi gives U conj(U) psi).
        """
        return self.U @ np.conj(self.U)


# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------


def build_J_weyl(n_max: int) -> RealStructure:
    """Build the J operator on the Weyl-sector spinor basis at cutoff n_max.

    Returns RealStructure with self.U a (dim, dim) unitary matrix realizing
    J as antilinear operator J(psi) = U @ conj(psi) with action

        J |n_fock, l, m_j> = sigma(l, m_j) |n_fock, l, -m_j>.

    On the Weyl sector ALONE, the J is a well-defined antilinear
    bijection because the Weyl-spinor basis {|n_fock, l, m_j> :
    j = l + 1/2 only} is closed under m_j -> -m_j.
    """
    basis = spinor_basis(n_max)
    label_to_index = {b: i for i, b in enumerate(basis)}
    dim = len(basis)
    U = np.zeros((dim, dim), dtype=np.complex128)
    for j, label in enumerate(basis):
        target = _spinor_J_target(label)
        i = label_to_index[target]
        sigma = _spinor_J_phase(target)
        U[i, j] = sigma
    return RealStructure(basis=basis, U=U, sector="weyl")


def build_J_full_dirac(n_max: int) -> RealStructure:
    """Build the J operator on the full-Dirac sector basis at cutoff n_max.

    Returns RealStructure with self.U a (dim, dim) unitary matrix realizing
    J on the doubled (Weyl + anti-Weyl chirality) Hilbert space:

        J |n_fock, l, m_j, chi> = sigma(l, m_j) |n_fock, l, -m_j, chi>.

    Block-diagonal in chirality: J does not flip chirality. This matches
    Paper 32's convention that J D = +D J (since D|.., chi> =
    chi*(n+1/2)|.., chi>, J flipping chi would give J D = -D J, wrong sign
    for KO-dim 3).
    """
    basis = full_dirac_basis(n_max)
    label_to_index = {b: i for i, b in enumerate(basis)}
    dim = len(basis)
    U = np.zeros((dim, dim), dtype=np.complex128)
    for j, label in enumerate(basis):
        target = _full_dirac_J_target(label)
        i = label_to_index[target]
        sigma = _spinor_J_phase(target)  # phase formula independent of chi
        U[i, j] = sigma
    return RealStructure(basis=basis, U=U, sector="full_dirac")


# ---------------------------------------------------------------------------
# Verification of Connes axioms
# ---------------------------------------------------------------------------


def verify_J_unitary(J: RealStructure, *, tol: float = 1e-12) -> Tuple[bool, float]:
    """Check that U is unitary: U U^dagger = I.

    Returns (is_unitary, max |U U^* - I|).
    """
    UUH = J.U @ J.U.conj().T
    I = np.eye(J.dim, dtype=np.complex128)
    err = float(np.max(np.abs(UUH - I)))
    return err < tol, err


def verify_J_squared(J: RealStructure, *, expected_sign: int = -1,
                     tol: float = 1e-12) -> Tuple[bool, float]:
    """Check J^2 = expected_sign * I.

    For KO-dim 3 (S^3), expected_sign = -1.

    Returns (is_correct, max |J^2 - expected_sign * I|).
    """
    J2 = J.J_squared_matrix()
    target = expected_sign * np.eye(J.dim, dtype=np.complex128)
    err = float(np.max(np.abs(J2 - target)))
    return err < tol, err


def verify_J_D_relation(J: RealStructure, D: np.ndarray, *,
                        expected_sign: int = +1,
                        tol: float = 1e-10) -> Tuple[bool, float]:
    """Check J D = expected_sign * D J as antilinear operators.

    For KO-dim 3, expected_sign = +1 (J commutes with D in the antilinear
    sense).

    The condition (J D - expected_sign D J)(psi) = 0 for all psi unpacks
    to: J D psi = expected_sign * D J psi for all psi.

    With J(psi) = U conj(psi), D linear:
        J D psi = U conj(D psi) = U conj(D) conj(psi)
        D J psi = D U conj(psi)

    Equality (up to sign) for all psi requires:
        U conj(D) = expected_sign * D U.

    Returns (is_correct, max |U conj(D) - expected_sign * D U|).
    """
    LHS = J.U @ np.conj(D)
    RHS = expected_sign * D @ J.U
    err = float(np.max(np.abs(LHS - RHS)))
    return err < tol, err


def verify_order_zero(J: RealStructure,
                      A_basis: Sequence[np.ndarray],
                      *,
                      tol: float = 1e-10) -> Tuple[bool, List[Tuple[int, int, float]]]:
    """Order-zero condition: [a, J b J^{-1}] = 0 for all a, b in A.

    For an antilinear J realized as J = U K, the J-conjugation of a
    linear operator b is

        J b J^{-1} = U conj(b) U^T.

    The order-zero condition is

        a (U conj(b) U^T) - (U conj(b) U^T) a = 0.

    Returns (all_pass, list of (i, j, residual) for any failing pair).
    """
    failures = []
    for i, a in enumerate(A_basis):
        for j, b in enumerate(A_basis):
            JbJinv = J.apply_to_operator(b)
            comm = a @ JbJinv - JbJinv @ a
            res = float(np.max(np.abs(comm)))
            if res > tol:
                failures.append((i, j, res))
    return (len(failures) == 0, failures)


def verify_order_one(J: RealStructure,
                     A_basis: Sequence[np.ndarray],
                     D: np.ndarray,
                     *,
                     tol: float = 1e-10) -> Tuple[bool, List[Tuple[int, int, float]]]:
    """Order-one condition: [[D, a], J b J^{-1}] = 0 for all a, b in A.

    Returns (all_pass, list of (i, j, residual) for any failing pair).
    """
    failures = []
    # Pre-compute commutators [D, a]
    Da_commutators = [D @ a - a @ D for a in A_basis]
    for i, Da in enumerate(Da_commutators):
        for j, b in enumerate(A_basis):
            JbJinv = J.apply_to_operator(b)
            comm = Da @ JbJinv - JbJinv @ Da
            res = float(np.max(np.abs(comm)))
            if res > tol:
                failures.append((i, j, res))
    return (len(failures) == 0, failures)


def verify_J_preserves_O(J: RealStructure,
                         op_sys,
                         *,
                         tol: float = 1e-9) -> Tuple[bool, List[Tuple[int, float]]]:
    """Check that J O J^{-1} = O (or equivalently each J a J^{-1} is in O
    when a is in O).

    This is an ancillary property: even if (C3) and (C4) hold or fail,
    the question of whether J b J^{-1} stays IN the operator system O
    is independent and structurally important.

    Returns (all_in_O, list of (i, residual) for failures).
    """
    failures = []
    for i, a in enumerate(op_sys.multiplier_matrices):
        JaJinv = J.apply_to_operator(a)
        in_O, residual = op_sys.contains(JaJinv, tol=tol)
        if not in_O:
            failures.append((i, residual))
    return (len(failures) == 0, failures)


# ---------------------------------------------------------------------------
# Convenience: high-level audit driver
# ---------------------------------------------------------------------------


def audit_J(n_max: int, sector: str = "full_dirac",
            D_mode: str = "truthful",
            *, verbose: bool = False) -> dict:
    """Run the full Connes-axiom audit for J at cutoff n_max.

    Parameters
    ----------
    n_max : int
        Fock cutoff.
    sector : str
        "weyl" or "full_dirac".
    D_mode : str
        For full_dirac sector, the Dirac mode:
        - "truthful": camporesi_higuchi_full_dirac_matrix (block-diagonal
          in chirality with eigenvalues +/- (n + 1/2)).
        - "offdiag": camporesi_higuchi_offdiag_dirac_matrix (with
          chirality-mixing off-diagonal couplings).
        For weyl sector: D = camporesi_higuchi_dirac_matrix (diagonal
        n + 1/2; D_mode is ignored).

    Returns
    -------
    dict with keys:
        "n_max", "sector", "D_mode",
        "dim_H", "dim_O" (operator system dimension),
        "U_unitary", "J_squared", "J_D_relation",
        "J_preserves_O",
        "order_zero", "order_one",
        each value (passes, max_residual_or_count).
    """
    if sector == "weyl":
        op_sys = SpinorTruncatedOperatorSystem(n_max)
        from geovac.spinor_operator_system import camporesi_higuchi_dirac_matrix
        D = camporesi_higuchi_dirac_matrix(op_sys.basis)
        J = build_J_weyl(n_max)
    elif sector == "full_dirac":
        op_sys = FullDiracTruncatedOperatorSystem(n_max)
        from geovac.full_dirac_operator_system import (
            camporesi_higuchi_full_dirac_matrix,
            camporesi_higuchi_offdiag_dirac_matrix,
        )
        if D_mode == "truthful":
            D = camporesi_higuchi_full_dirac_matrix(op_sys.basis)
        elif D_mode == "offdiag":
            D = camporesi_higuchi_offdiag_dirac_matrix(op_sys.basis)
        else:
            raise ValueError(f"unknown D_mode {D_mode!r}")
        J = build_J_full_dirac(n_max)
    else:
        raise ValueError(f"unknown sector {sector!r}")

    if verbose:
        print(f"audit_J(n_max={n_max}, sector={sector}, D_mode={D_mode})")
        print(f"  dim_H = {op_sys.dim_H}, dim_O = {op_sys.dim}, "
              f"|generators| = {len(op_sys.multiplier_matrices)}")

    U_ok, U_res = verify_J_unitary(J)
    J2_ok, J2_res = verify_J_squared(J)
    JD_ok, JD_res = verify_J_D_relation(J, D)
    Jpreserve_ok, Jpreserve_failures = verify_J_preserves_O(J, op_sys)
    order0_ok, order0_failures = verify_order_zero(J, op_sys.multiplier_matrices)
    order1_ok, order1_failures = verify_order_one(J, op_sys.multiplier_matrices, D)

    result = {
        "n_max": n_max,
        "sector": sector,
        "D_mode": D_mode,
        "dim_H": op_sys.dim_H,
        "dim_O": op_sys.dim,
        "n_generators": len(op_sys.multiplier_matrices),
        "U_unitary": (U_ok, U_res),
        "J_squared_minus_I": (J2_ok, J2_res),
        "J_commutes_D": (JD_ok, JD_res),
        "J_preserves_O": (Jpreserve_ok, len(Jpreserve_failures),
                          (max(r for _, r in Jpreserve_failures)
                           if Jpreserve_failures else 0.0)),
        "order_zero": (order0_ok, len(order0_failures),
                       (max(r for _, _, r in order0_failures)
                        if order0_failures else 0.0)),
        "order_one": (order1_ok, len(order1_failures),
                      (max(r for _, _, r in order1_failures)
                       if order1_failures else 0.0)),
    }

    if verbose:
        print(f"  U unitary: {U_ok} (max err {U_res:.3e})")
        print(f"  J^2 = -I: {J2_ok} (max err {J2_res:.3e})")
        print(f"  J D = +D J: {JD_ok} (max err {JD_res:.3e})")
        print(f"  J preserves O: {Jpreserve_ok} "
              f"(failures {len(Jpreserve_failures)}, "
              f"max res {result['J_preserves_O'][2]:.3e})")
        print(f"  order zero [a, J b J^-1] = 0: {order0_ok} "
              f"(failures {len(order0_failures)}, "
              f"max res {result['order_zero'][2]:.3e})")
        print(f"  order one [[D, a], J b J^-1] = 0: {order1_ok} "
              f"(failures {len(order1_failures)}, "
              f"max res {result['order_one'][2]:.3e})")

    return result
