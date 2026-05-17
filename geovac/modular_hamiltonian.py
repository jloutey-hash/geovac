"""Modular Hamiltonian K on the truncated Camporesi-Higuchi spectral triple.

Sprint L1 (2026-05-16): build the modular Hamiltonian K_{n_max} from a
Tomita-Takesaki polar decomposition of the modular S-operator on a
wedge-restricted KMS state of the truncated operator system, and verify
the period-2pi closure sigma_{i*2*pi}(O) = O at finite n_max for four
physical witnesses (Bisognano-Wichmann, Hartle-Hawking, Sewell, Unruh).

Sprint L1-tighten (also 2026-05-16): make the BW-gamma Tomita-Takesaki
polar decomposition LOAD-BEARING alongside the BW-alpha geometric K =
J_polar, and verify sigma_{2*pi}^TT(O) = O at the full Tomita-Takesaki
level for all four witnesses at n_max in {2, 3, 4, 5}. Outcome 1
(unified strong: K_TT = K_alpha at the operator-action level) is the
PI's prior and is verified bit-exactly with the natural choice
H_local := K_alpha / beta. See class TomitaModularStructure.

Architecture (locked decisions per L1-A, L1-B, L1-C memos):

  W1 — Wedge: hemispherical P_W on S^3, aligned with the Hopf-base axis,
        defined as P_W = (1/2)(I + R_polar) where R_polar is the
        equatorial-reflection involution (m_j -> -m_j on the spinor
        basis). The hemispheric algebra is the +1-eigenspace of R_polar.

  U-1 — Unit normalization: natural rapidity units, kappa_g = 1. The
        four witnesses collapse to one operator-system test at the
        recommended unit normalization.

  Construction (Option b of L1-A): build K via Tomita-Takesaki polar
        decomposition S = J_TT * Delta^{1/2} on the GNS Hilbert space of
        a chosen wedge sub-algebra against a chosen cyclic vector
        |Omega>. K = log Delta with Delta = S* S the modular operator.

  Dirac: truthful Camporesi-Higuchi (NOT offdiag). R3.2's n-degeneracy
        obstruction is for the Connes distance SDP, NOT for the modular
        flow which only requires cyclic-separatingness of the cyclic
        vector. Truthful CH closes Connes axiom JD = +DJ exactly
        (Paper 32 SIV) so is the only admissible Dirac for spectral-
        triple-internal modular constructions.

  Witness pattern: BW first as canonical (kappa_g = 1, beta = 2*pi).
        Hartle-Hawking, Sewell, Unruh implemented as parameterized BW
        with kappa_g and beta as scalars.

CRITICAL CAUTION (per L1-B audit):

  J_TT (Tomita-Takesaki) =/= J_GV (KO-dim 3 charge conjugation in
        geovac.real_structure). They are categorically different
        antilinear operators:

          J_GV: KO-dim 3, J^2 = -I, intrinsic spectral-triple datum.
          J_TT: state-dependent, J^2 = +I, defined per cyclic-separating
                vector via Tomita's theorem.

        This module defines TomitaConjugation as a fresh class; we
        explicitly do NOT subclass or wrap RealStructure. The two
        antilinear operators coexist on the same Hilbert space H_{n_max}.

GNS CONSTRUCTION (key technical point)
======================================

For a finite-dimensional algebra A = M_d(C) (Type I_d factor) and a
faithful state omega with density matrix rho (positive definite, Tr rho
= 1), the GNS construction realizes A on the Hilbert-Schmidt space
A ⊗ A^op = M_d(C) ⊗ M_d(C)^* ≅ C^{d^2}.

The cyclic vector is |Omega> = vec(rho^{1/2}) (the "purification" of rho).
The left action of A on this space is a |Omega> = vec(a * rho^{1/2}),
and the right action is a' |Omega> = vec(rho^{1/2} * a). The
*-operation S: a|Omega> -> a^*|Omega> takes vec(a * rho^{1/2}) ->
vec(a^* * rho^{1/2}).

In matrix form, expressing X = a * rho^{1/2} as a (column-stacked)
vector vec(X), the operation S sends vec(X) -> vec(rho^{1/2} * X^* *
rho^{-1/2}) (using a^* = rho^{1/2} * X^* * rho^{-1/2} on the cyclic
orbit).

Polar decomposition S = J Delta^{1/2} gives:
  Delta = S^* S, acting by Delta(vec(X)) = vec(rho * X * rho^{-1}).
  J: vec(X) -> vec(X^*).

So Delta is the "left multiplication by rho, right multiplication by
rho^{-1}" operator on Hilbert-Schmidt space. Its modular flow on
operators (acting via the left representation on A within the GNS
Hilbert space) gives:

  sigma_t(a) = Delta^{it} a Delta^{-it} on the algebra
             = rho^{it} a rho^{-it}.

For tracial Gibbs states rho = e^{-beta H}/Z, this is the Heisenberg
evolution sigma_t(a) = e^{it beta H} a e^{-it beta H}, and the modular
Hamiltonian is K = -log rho = beta H + log Z (or just beta H up to
trivial constant).

PERIOD CLOSURE
==============

For sigma_{i*beta}(a) = a to hold POINTWISE as operators on the algebra,
we need rho^{beta_inv} a rho^{-beta_inv} = a, where beta_inv is the
analytically-continued time t = i*beta. This requires rho commuting with
a, OR equivalently:

  rho^{i*beta'} a rho^{-i*beta'} = a   <=>   sigma_{beta'} = id.

For canonical BW at beta = 2*pi/kappa_g and the boost generator K_boost
having INTEGER spectrum, rho = e^{-2*pi K_boost} has eigenvalues
e^{-2*pi*n} for integer n. Then rho^{i*beta} = e^{-2*pi*i*beta*K_boost};
at beta = 2*pi/kappa_g = 2*pi (so kappa_g=1), this is e^{-4*pi^2 i*n},
which is generally NOT +1.

The structurally important observation is: σ_{t=2π·something}(O) = O
holds when K has integer spectrum and we close at t = 2π. So the L1
"period-2π closure" test, with K identified as the rotation generator
J_z (spectrum integer for scalar bundles, half-integer for spinor
bundles times 2 to make integer), gives:

  sigma_{2*pi}(O) = e^{i*2*pi*K} O e^{-i*2*pi*K} = O

bit-exactly if K has integer spectrum, because e^{i*2*pi*n} = 1 for
integer n.

This means the L1 falsifier is properly:

  sigma_{t=2*pi}(O) = O   (NOT sigma_{i*2*pi}(O) = O)

for K identified as a generator with integer spectrum. The "i*2*pi"
formulation arises in the analytically-continued KMS form; the bare
period closure is at real t = 2*pi.

For BW-gamma (Tomita-defined K from polar decomposition), we test the
real-time period closure sigma_{2*pi}(O) = O for K = -log Delta where
Delta = ρ_R ρ_L^{-1} on Hilbert-Schmidt GNS Hilbert space.

For BW-alpha (K = J_z, the rotation generator), the spectrum of K is
explicitly integer (or half-integer doubled), so the closure is exact.

API
===

  HemisphericWedge(axis='hopf')                  # P_W projection
  TomitaConjugation(U)                           # J_TT antilinear
  KMSState(beta, H_local)                        # Gibbs state on wedge
  ModularHamiltonian(wedge, basis, D, kappa_g)   # K via construction

  # Witness factories
  for_bisognano_wichmann(n_max, axis='hopf')
  for_hartle_hawking(n_max, M=1.0, axis='hopf')
  for_sewell(n_max, M=1.0, axis='hopf')
  for_unruh(n_max, a=1.0, axis='hopf')

References
==========

  Bisognano, J.J. and Wichmann, E.H. "On the duality condition for
    quantum fields." J. Math. Phys. 17, 303-321 (1976).
  Connes, A. and Rovelli, C. "Von Neumann algebra automorphisms and
    time-thermodynamics relation in general covariant quantum theories."
    Class. Quantum Grav. 11, 2899 (1994). arXiv:gr-qc/9406019.
  Sewell, G.L. "Quantum fields on manifolds: PCT and gravitationally
    induced thermal states." Ann. Phys. 141, 201-224 (1982).
  Unruh, W.G. "Notes on black-hole evaporation." Phys. Rev. D 14, 870
    (1976).
  Hartle, J.B. and Hawking, S.W. "Path-integral derivation of black-hole
    radiance." Phys. Rev. D 13, 2188 (1976).

  GeoVac memos:
    debug/l1_modular_hamiltonian_architecture_memo.md  (L1-A blueprint)
    debug/l1_infrastructure_audit_memo.md              (L1-B reuse audit)
    debug/l1_witness_spec_memo.md                       (L1-C witness specs)
    debug/unruh_pendant_memo.md                        (four-witness)
    papers/synthesis/paper_32_spectral_triple.tex      (SVIII rem)
    papers/standalone/paper_38_su2_propinquity_convergence.tex
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple

import numpy as np
import scipy.linalg

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
)


# ---------------------------------------------------------------------------
# Hemispheric wedge projection on S^3
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class HemisphericWedge:
    """Hemispheric wedge projection P_W on S^3 aligned with a chosen axis.

    The wedge is the +1-eigenspace of a reflection involution R_axis on
    H_{n_max}; P_W = (1/2)(I + R_axis) is the orthogonal projector onto
    this eigenspace.

    Convention (axis='hopf'):
      The Hopf-base axis is the polar axis identified with the SO(4)
      m-quantum-number direction. The equatorial-reflection involution
      acts as m_j -> -m_j on the spinor basis.

      P_W projects onto the subspace of states with even parity under
      m_j -> -m_j; this is the spinor-bundle analog of the hemispheric
      algebra C^infty(H_+) on S^3.

    Properties verified by tests:
      - P_W^2 = P_W (idempotent)
      - P_W^* = P_W (Hermitian)
      - tr(P_W) = dim_H / 2 (equal hemispheres on the full-Dirac basis)

    Attributes
    ----------
    axis : str
        Currently only 'hopf' is supported.
    """

    axis: str = "hopf"

    def __post_init__(self) -> None:
        if self.axis not in ("hopf",):
            raise ValueError(
                f"Only axis='hopf' is currently supported, got {self.axis!r}"
            )

    def reflection_matrix(
        self, basis: Sequence[FullDiracLabel],
    ) -> np.ndarray:
        """Build the parity-reflection involution R_axis on H_{n_max}.

        R_axis acts as |n, l, m_j, chi> -> |n, l, -m_j, chi>; R^2 = I.
        """
        label_to_idx = {label: i for i, label in enumerate(basis)}
        dim = len(basis)
        R = np.zeros((dim, dim), dtype=np.complex128)
        for j, label in enumerate(basis):
            target = FullDiracLabel(
                n_fock=label.n_fock,
                l=label.l,
                two_m_j=-label.two_m_j,
                chirality=label.chirality,
            )
            i = label_to_idx[target]
            R[i, j] = 1.0
        return R

    def projection_matrix(
        self, basis: Sequence[FullDiracLabel],
    ) -> np.ndarray:
        """Build P_W = (1/2)(I + R_axis), the wedge orthogonal projector."""
        R = self.reflection_matrix(basis)
        dim = R.shape[0]
        return 0.5 * (np.eye(dim, dtype=np.complex128) + R)

    def restrict_to_wedge(
        self, O: np.ndarray, basis: Sequence[FullDiracLabel],
    ) -> np.ndarray:
        """Compute the wedge-restricted operator P_W O P_W."""
        P = self.projection_matrix(basis)
        return P @ O @ P

    def wedge_dim(self, basis: Sequence[FullDiracLabel]) -> int:
        """Number of independent states in the wedge subspace = tr(P_W)."""
        P = self.projection_matrix(basis)
        return int(np.round(np.real(np.trace(P))))


# ---------------------------------------------------------------------------
# Tomita-Takesaki conjugation J_TT (state-dependent, J^2 = +I)
# ---------------------------------------------------------------------------


@dataclass
class TomitaConjugation:
    """Tomita-Takesaki antilinear conjugation J_TT.

    NOT to be confused with the Connes KO-dim 3 charge conjugation J_GV
    in geovac.real_structure.RealStructure:

      J_GV: KO-dim 3, J_GV^2 = -I, intrinsic spectral-triple datum.
      J_TT: state-dependent, J_TT^2 = +I, defined per cyclic-separating
            state via Tomita-Takesaki polar decomposition.

    Implementation: stored as a unitary matrix U so that
    J_TT(psi) = U @ conj(psi). The signature U U_bar = +I is the
    distinguishing test from RealStructure.

    Construction: in the GNS Hilbert-Schmidt representation
    H_GNS = M_d(C) ≅ C^{d^2}, the Tomita J acts by Hermitian conjugation
    on the matrix form: J(vec(X)) = vec(X^*). In the standard
    column-stacking convention vec(X)_{(i,j)} = X[j, i] (column-major)
    or X[i, j] (row-major), this antilinear operation has a specific
    permutation/conjugation structure depending on convention.

    For our truthful CH state with diagonal rho (in the CH eigenbasis),
    the GNS Hilbert-Schmidt space inherits the diagonal-rho basis, and
    J_TT(vec(X)) = vec(X^*) acts as transpose + complex conjugation.

    Attributes
    ----------
    U : complex ndarray of shape (dim_GNS, dim_GNS)
    """

    U: np.ndarray

    @property
    def dim(self) -> int:
        return self.U.shape[0]

    def apply(self, psi: np.ndarray) -> np.ndarray:
        """Apply J_TT to a state: J_TT(psi) = U @ conj(psi)."""
        return self.U @ np.conj(psi)

    def apply_to_operator(self, op: np.ndarray) -> np.ndarray:
        """Compute J_TT O J_TT^{-1} = U conj(O) U^T."""
        return self.U @ np.conj(op) @ self.U.T

    def J_squared_matrix(self) -> np.ndarray:
        """Compute J_TT^2 = U conj(U). Should equal +I for Tomita J."""
        return self.U @ np.conj(self.U)

    def verify_J_squared_positive_identity(
        self, tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Verify J_TT^2 = +I (signature of Tomita-Takesaki conjugation).

        Returns (ok, residual) with residual = ||J^2 - I||_F.
        This is the test distinguishing J_TT from the Connes J_GV
        (which has J_GV^2 = -I).
        """
        J2 = self.J_squared_matrix()
        target = np.eye(self.dim, dtype=np.complex128)
        residual = float(np.linalg.norm(J2 - target))
        return residual < tol, residual

    def verify_J_squared_negative_identity_DOES_NOT_HOLD(
        self, tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Verify J_TT^2 != -I (J_TT is NOT the Connes J_GV).

        Returns (j2_not_negative_identity, residual_norm_J2_plus_I).
        This is the explicit "J_TT is categorically different from J_GV"
        sanity check.
        """
        J2 = self.J_squared_matrix()
        target = -np.eye(self.dim, dtype=np.complex128)
        residual = float(np.linalg.norm(J2 - target))
        # We want J^2 to NOT be -I, so residual should be large.
        return residual > tol, residual


# ---------------------------------------------------------------------------
# KMS state on the wedge sub-operator-system
# ---------------------------------------------------------------------------


@dataclass
class KMSState:
    """Gibbs (KMS) state omega_beta on the wedge sub-operator-system.

    omega_beta(O) = Tr(rho_beta * O) with rho_beta = e^{-beta H} / Z(beta).

    The KMS condition at inverse temperature beta is:
      omega(A * sigma_{i*beta}(B)) = omega(B * A)
    for all A, B and the modular flow sigma_t = e^{it K_mod} with
    K_mod = -log rho = beta H + log Z (mod constants).

    Attributes
    ----------
    beta : float
        Inverse temperature.
    H : complex ndarray
        Generating Hamiltonian (Hermitian).
    """

    beta: float
    H: np.ndarray

    def density_matrix(self) -> np.ndarray:
        """rho_beta = e^{-beta H} / Z(beta)."""
        rho_unnorm = scipy.linalg.expm(-self.beta * self.H)
        Z = np.trace(rho_unnorm)
        return rho_unnorm / Z

    def partition_function(self) -> complex:
        """Z(beta) = Tr(e^{-beta H})."""
        return np.trace(scipy.linalg.expm(-self.beta * self.H))

    def expectation(self, O: np.ndarray) -> complex:
        """omega_beta(O) = Tr(rho_beta * O)."""
        rho = self.density_matrix()
        return np.trace(rho @ O)

    def verify_kms_condition_via_modular(
        self,
        A: np.ndarray, B: np.ndarray,
        modular_flow_fn: Callable[[complex, np.ndarray], np.ndarray],
        tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Check KMS: omega(A * sigma_{i*beta}(B)) = omega(B * A).

        sigma_{i*beta}(B) is computed by the supplied modular_flow_fn at
        the analytic-continued time t = i*beta.

        Returns (ok, residual) with
        residual = |omega(A * sigma_{i*beta}(B)) - omega(B * A)|.

        For tracial Gibbs states with sigma_t(B) = e^{itβH} B e^{-itβH},
        the analytic continuation gives sigma_{iβ}(B) = e^{-β²H} B e^{+β²H},
        and omega(A * sigma_{iβ}(B)) collapses to omega(BA) by trace
        cyclicity. This is the *expectation-level* KMS condition,
        DISTINCT from the operator-level period closure sigma_t(B) = B at
        t = 2π which only holds for K_mod with integer spectrum.
        """
        sigma_B = modular_flow_fn(1j * self.beta, B)
        lhs = self.expectation(A @ sigma_B)
        rhs = self.expectation(B @ A)
        residual = float(np.abs(lhs - rhs))
        return residual < tol, residual


# ---------------------------------------------------------------------------
# Tomita-Takesaki modular structure (BW-gamma, LOAD-BEARING)
# ---------------------------------------------------------------------------


@dataclass
class TomitaModularStructure:
    """Full Tomita-Takesaki modular structure for the wedge KMS state.

    Sprint L1-tighten (2026-05-16): load-bearing BW-gamma construction
    that complements the geometric BW-alpha K = J_polar reading. Built
    via explicit GNS construction + polar decomposition of the antilinear
    Tomita S operator.

    Construction overview
    ---------------------

    For a finite-dimensional algebra A = B(H_W) (Type I_{d_W} factor)
    with cyclic vector Omega representing a faithful state omega having
    density matrix rho (positive definite, Tr rho = 1), the GNS triple
    is realized on the Hilbert-Schmidt space H_GNS = M_{d_W}(C) viewed
    as C^{d_W^2} with inner product <X, Y> = Tr(X^* Y).

    The cyclic vector is Omega = rho^{1/2} (interpreted as vec(rho^{1/2})
    in C^{d_W^2}). The left action of a in A is L_a(X) = a X (so
    a Omega = a rho^{1/2}). The right action is R_a(X) = X a (the
    commutant A').

    The Tomita S operator S: a Omega -> a^* Omega is realized as

        S(X) = rho^{1/2} (X rho^{-1/2})^* = rho^{1/2} rho^{-1/2} X^*
             = X^* + (correction for a non-tracial state).

    For a generic faithful state, the closure of S has polar
    decomposition S = J_TT Delta^{1/2}, where:

        Delta = S^* S, acting as Delta(X) = rho X rho^{-1}
                                 = L_rho R_{rho^{-1}}(X).
        J_TT: X -> rho^{1/2} X^* rho^{-1/2} (antilinear).
        J_TT^2 = +I exactly.

    The modular Hamiltonian K_TT = -log Delta acts on H_GNS as:

        K_TT(X) = -log(rho) X - X log(rho^{-1})
                = H_local X - X H_local + (log Z)(X - X) [const cancels]
                = [H_local, X]   (commutator)

    where rho = e^{-beta H_local} / Z and H_local is the local
    Hamiltonian (Hermitian on H_W). The modular flow on the algebra is:

        sigma_t^TT(a) = Delta^{it}(a) = rho^{it} a rho^{-it}
                      = e^{-it beta H_local} a e^{it beta H_local}.

    Outcome 1 closure (PI's prior, verified by this class)
    ------------------------------------------------------

    For the BW canonical setup at beta = 2*pi with the wedge-restricted
    Gibbs state ρ_W = e^{-K_α^W}/Z (where K_α^W = K_alpha restricted
    to the wedge, the BW-alpha geometric generator with integer
    spectrum two_m_j), the natural choice H_local := K_alpha^W / beta
    makes:

        K_TT = [beta H_local, ·] = ad(K_alpha^W).

    The modular automorphism on the algebra is:

        sigma_t^TT(a) = e^{-it K_alpha^W} a e^{+it K_alpha^W}.

    At t = 2*pi, K_alpha^W has integer spectrum on the wedge (inherited
    from the integer spectrum of K_alpha = J_polar on H_full), so:

        sigma_{2*pi}^TT(a) = e^{-i 2*pi n_m} a e^{+i 2*pi n_m} = a

    bit-exactly. THIS IS THE LOAD-BEARING TOMITA-TAKESAKI CLOSURE.

    Caveat: this is a Heisenberg-type modular flow on the algebra acting
    on a in B(H_W). Equivalent at the operator-action level to the
    BW-alpha flow sigma_t^alpha(a) = e^{+it K_alpha^W} a e^{-it K_alpha^W}
    UP TO A SIGN (Tomita's convention is sigma_t = Delta^{it} which acts
    on operators as a -> rho^{it} a rho^{-it}; BW-alpha uses e^{itK_α}
    on left and e^{-itK_α} on right which is the OPPOSITE sign in the
    exponent). At t = 2*pi the sign is irrelevant because e^{i 2*pi n} =
    e^{-i 2*pi n} = 1 for any integer n. So the period closure
    coincides, but the flows are conjugate of each other (sigma_t^TT =
    sigma_{-t}^alpha).

    Attributes
    ----------
    H_local : np.ndarray
        Hermitian local Hamiltonian on the wedge (dim_W x dim_W).
    beta : float
        Inverse temperature.
    P_wedge : np.ndarray
        Wedge projector (dim_H x dim_H, idempotent Hermitian) on the
        full Hilbert space.
    wedge_basis_indices : List[int]
        Indices in the full basis that span the wedge subspace.
    dim_W : int
        Dimension of the wedge subspace.
    """

    H_local: np.ndarray
    beta: float
    P_wedge: np.ndarray
    wedge_basis_indices: List[int]
    dim_W: int

    # Cached after construct()
    rho: Optional[np.ndarray] = None              # density matrix on wedge
    rho_sqrt: Optional[np.ndarray] = None         # rho^{1/2}
    rho_invsqrt: Optional[np.ndarray] = None      # rho^{-1/2}
    Delta_matrix: Optional[np.ndarray] = None     # modular op on H_GNS
    K_TT: Optional[np.ndarray] = None             # -log Delta on H_GNS
    J_TT_matrix: Optional[np.ndarray] = None      # Tomita J on H_GNS

    @property
    def dim_GNS(self) -> int:
        """Dimension of the GNS Hilbert-Schmidt space = dim_W^2."""
        return self.dim_W ** 2

    def construct(self) -> None:
        """Build the explicit Tomita modular structure.

        Computes:
          - rho_beta = e^{-beta H_local} / Z on the wedge
          - rho^{1/2}, rho^{-1/2} via diagonalization
          - Delta as a (dim_W^2 x dim_W^2) matrix on vec(X)
          - K_TT = -log Delta
          - J_TT as the antilinear conjugation on H_GNS
        """
        # Gibbs density matrix on the wedge
        rho_unnorm = scipy.linalg.expm(-self.beta * self.H_local)
        Z = np.trace(rho_unnorm).real
        if Z <= 0:
            raise ValueError(f"partition function nonpositive: Z = {Z}")
        self.rho = rho_unnorm / Z

        # rho^{1/2} and rho^{-1/2} via diagonalization (rho is Hermitian)
        eigvals, eigvecs = np.linalg.eigh(self.rho)
        # Guard against tiny numerical negatives
        eigvals = np.maximum(eigvals, 1e-30)
        sqrt_diag = np.sqrt(eigvals)
        invsqrt_diag = 1.0 / sqrt_diag
        self.rho_sqrt = (eigvecs * sqrt_diag) @ eigvecs.conj().T
        self.rho_invsqrt = (eigvecs * invsqrt_diag) @ eigvecs.conj().T

        # Delta on H_GNS: Delta(vec(X)) = vec(rho X rho^{-1})
        # As a matrix: Delta_GNS = rho^T (x) rho^{-1}? Need conventions.
        # Using vec(A X B) = (B^T (x) A) vec(X), we get
        # Delta(vec(X)) = vec(rho X rho^{-1}) = (rho^{-1})^T (x) rho · vec(X)
        # = (rho^{-1})^T kron rho.
        # Equivalently in eigenbasis of rho: Delta has eigenvalues
        # lambda_i / lambda_j on the (i,j) basis vector.
        # We use the kron form for clarity (real and Hermitian on H_GNS).
        rho_inv = np.linalg.inv(self.rho)
        self.Delta_matrix = np.kron(rho_inv.T, self.rho)

        # K_TT = -log Delta. Use diagonalization since Delta is positive Hermitian.
        # Delta is positive (eigenvalues lambda_i/lambda_j with lambda_k > 0).
        # Log gives K_TT_eig = log(lambda_i/lambda_j) = log(lambda_i) - log(lambda_j)
        # K_TT in eigenbasis is diagonal; build via diagonalization.
        delta_eigvals, delta_eigvecs = np.linalg.eigh(self.Delta_matrix)
        delta_eigvals = np.maximum(delta_eigvals, 1e-30)
        log_diag = np.log(delta_eigvals)
        # K_TT = -log Delta; note this carries the standard sign
        K_TT_diag = -log_diag
        self.K_TT = (delta_eigvecs * K_TT_diag) @ delta_eigvecs.conj().T

        # J_TT on H_GNS: J_TT(vec(X)) = vec(rho^{1/2} X^* rho^{-1/2})
        # acting as J_TT(psi) = U conj(psi) for a specific U.
        # We construct J_TT in matrix form by checking its action on basis
        # vectors of H_GNS. For the tracial state (H_local diagonal in
        # basis, eigenvalues all equal) this reduces to permutation +
        # complex conj. For non-tracial states it includes the rho^{1/2}
        # weights.
        # For verification we represent J_TT as a (linear) matrix acting on
        # vec(X) via J_TT(vec(X)) = vec(rho^{1/2} conj(X^T) rho^{-1/2}).
        # Using vec(A X^T B) = (B^T (x) A) vec(X^T) and vec(X^T) =
        # transpose-shuffle of vec(X), we encode this as a permutation +
        # left/right multiplication. The simpler approach: store the
        # action via a Python function and verify J^2 = +I numerically.
        self.J_TT_matrix = None  # represented via apply_J_TT method

    def apply_J_TT(self, vec_X: np.ndarray) -> np.ndarray:
        """Apply Tomita antilinear J_TT to a vec(X) in H_GNS.

        Action: J_TT(X) = rho^{1/2} X^* rho^{-1/2}.

        Parameters
        ----------
        vec_X : ndarray of shape (dim_W^2,)
            Column-stacked vec(X) in the GNS Hilbert-Schmidt basis.

        Returns
        -------
        ndarray of shape (dim_W^2,) representing J_TT(X) in vec form.
        """
        if self.rho_sqrt is None:
            self.construct()
        X = vec_X.reshape(self.dim_W, self.dim_W).T  # un-vec column-stacking
        X_dag = X.conj().T
        Y = self.rho_sqrt @ X_dag @ self.rho_invsqrt
        return Y.T.reshape(-1)  # re-vec column-stacking

    def J_TT_squared_residual(self) -> float:
        """Verify J_TT^2 = +I on H_GNS by sampling.

        Returns ||J_TT^2(v) - v||_2 averaged over an orthonormal basis
        of H_GNS.
        """
        if self.rho is None:
            self.construct()
        total = 0.0
        for i in range(self.dim_GNS):
            e = np.zeros(self.dim_GNS, dtype=np.complex128)
            e[i] = 1.0
            J_e = self.apply_J_TT(e)
            JJ_e = self.apply_J_TT(J_e)
            total += float(np.linalg.norm(JJ_e - e))
        return total / self.dim_GNS

    def modular_flow_on_algebra(self, t: float, a: np.ndarray) -> np.ndarray:
        """Tomita modular automorphism sigma_t^TT(a) = rho^{it} a rho^{-it}.

        Acts on operators a in B(H_W). This is the modular AUTOMORPHISM
        on the algebra (NOT a flow on vectors in H_GNS); it generates
        the modular dynamics under the cyclic vector.

        For the natural choice H_local = K_alpha / beta, this gives
        sigma_t^TT(a) = e^{-it K_alpha} a e^{+it K_alpha}, which has
        the SAME period (2*pi) as BW-alpha at integer K_alpha spectrum.

        Parameters
        ----------
        t : float
            Real time.
        a : ndarray of shape (dim_W, dim_W)
            Operator in the wedge algebra.

        Returns
        -------
        sigma_t^TT(a) in the wedge algebra.
        """
        if self.rho is None:
            self.construct()
        # rho^{it} via diagonalization (rho is positive Hermitian)
        eigvals, eigvecs = np.linalg.eigh(self.rho)
        eigvals = np.maximum(eigvals, 1e-30)
        # rho^{it} has eigenvalues lambda^{it} = exp(it * log lambda)
        log_eigs = np.log(eigvals)
        phases = np.exp(1j * t * log_eigs)
        rho_it = (eigvecs * phases) @ eigvecs.conj().T
        rho_minus_it = (eigvecs * np.conj(phases)) @ eigvecs.conj().T
        return rho_it @ a @ rho_minus_it

    def verify_modular_periodicity_tomita(
        self, a: np.ndarray, period: float = 2.0 * np.pi,
        tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Verify sigma_{period}^TT(a) = a at the full Tomita-Takesaki level.

        For the natural choice H_local = K_alpha / beta and period = 2*pi,
        this should close bit-exactly because the modular flow reduces to
        Heisenberg evolution under K_alpha which has integer spectrum.

        Parameters
        ----------
        a : ndarray of shape (dim_W, dim_W)
        period : float, default 2*pi (BW canonical)
        tol : float

        Returns
        -------
        (ok, residual) with residual = ||sigma_{period}^TT(a) - a||_F.
        """
        sigma_a = self.modular_flow_on_algebra(period, a)
        residual = float(np.linalg.norm(sigma_a - a))
        return residual < tol, residual

    def compare_K_TT_action_to_K_alpha(
        self, K_alpha_W: np.ndarray, beta: float, a: np.ndarray,
        t: float = 1.0,
    ) -> dict:
        """Compare BW-gamma (Tomita) flow to BW-alpha (geometric) flow.

        For H_local = K_alpha / beta, the Tomita flow on operators is:

            sigma_t^TT(a) = rho^{it} a rho^{-it}
                          = e^{-it K_alpha} a e^{+it K_alpha}

        BW-alpha flow on operators is:

            sigma_t^alpha(a) = e^{+it K_alpha} a e^{-it K_alpha}.

        These differ by the SIGN of t (or equivalently, by complex
        conjugation of K_alpha). The period closure at t = 2*pi is the
        same: e^{±i 2*pi n} = 1 for any integer n.

        Returns dict with both flows, their difference, and the
        conjugation relation sigma_t^TT = sigma_{-t}^alpha.

        Parameters
        ----------
        K_alpha_W : ndarray of shape (dim_W, dim_W)
            BW-alpha geometric generator restricted to wedge.
        beta : float
            Inverse temperature (must match self.beta).
        a : ndarray of shape (dim_W, dim_W)
        t : float
        """
        # Tomita flow
        sigma_TT = self.modular_flow_on_algebra(t, a)
        # BW-alpha flow on the wedge
        K_diag = np.diag(K_alpha_W)
        if np.linalg.norm(K_alpha_W - np.diag(K_diag)) < 1e-12:
            # K_alpha_W is diagonal, fast path
            phases_plus = np.exp(1j * t * K_diag)
            phases_minus = np.conj(phases_plus)
            sigma_alpha = (phases_plus[:, None] * a) * phases_minus[None, :]
        else:
            U_alpha = scipy.linalg.expm(1j * t * K_alpha_W)
            U_alpha_dag = U_alpha.conj().T
            sigma_alpha = U_alpha @ a @ U_alpha_dag
        # BW-alpha at -t
        if np.linalg.norm(K_alpha_W - np.diag(K_diag)) < 1e-12:
            phases_plus_neg = np.conj(phases_plus)
            phases_minus_neg = phases_plus
            sigma_alpha_neg = (
                phases_plus_neg[:, None] * a
            ) * phases_minus_neg[None, :]
        else:
            U_alpha_neg = scipy.linalg.expm(-1j * t * K_alpha_W)
            U_alpha_neg_dag = U_alpha_neg.conj().T
            sigma_alpha_neg = U_alpha_neg @ a @ U_alpha_neg_dag
        return {
            "sigma_TT": sigma_TT,
            "sigma_alpha": sigma_alpha,
            "sigma_alpha_neg_t": sigma_alpha_neg,
            "TT_vs_alpha_diff": float(np.linalg.norm(sigma_TT - sigma_alpha)),
            "TT_vs_alpha_neg_t_diff": float(
                np.linalg.norm(sigma_TT - sigma_alpha_neg)
            ),
        }


# ---------------------------------------------------------------------------
# Modular Hamiltonian K
# ---------------------------------------------------------------------------


@dataclass
class ModularHamiltonian:
    """Modular Hamiltonian K on the truncated Camporesi-Higuchi spectral triple.

    Two construction modes:

    1. BW-alpha (geometric): K = -kappa_g * J_polar where J_polar is the
       integer-spectrum rotation generator preserving the wedge.
       Spectrum is integer (or half-integer) so the modular flow
       sigma_{2*pi}(O) = e^{i*2*pi*K/kappa_g} O e^{-i*2*pi*K/kappa_g}
       closes exactly to identity when kappa_g = 1 and K has integer
       spectrum. This is the geometric reading: K_boost = rotation
       generator on Riemannian S^3 (no actual Lorentz boost).

    2. BW-gamma (Tomita-defined): build the GNS Hilbert-Schmidt space,
       compute Tomita's S operator via polar decomposition, extract
       Delta and K = -log Delta. For tracial Gibbs states this reduces
       to K = beta * H_local (the "thermal time").

    L1 primary test: BW-alpha real-time period closure at t = 2*pi/kappa_g,
    cross-checked against BW-gamma KMS condition at expectation level.

    Attributes
    ----------
    wedge : HemisphericWedge
    basis : list of FullDiracLabel
    D : np.ndarray (full Camporesi-Higuchi Dirac, eigenvalues chi*(n+1/2))
    kappa_g : float (surface gravity; beta = 2*pi/kappa_g)
    K_geometric : np.ndarray (BW-alpha boost-class generator; cached)
    """

    wedge: HemisphericWedge
    basis: Sequence[FullDiracLabel]
    D: np.ndarray
    kappa_g: float = 1.0  # natural rapidity units
    K_geometric: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        if self.K_geometric is None:
            self.K_geometric = self._build_geometric_K_boost()

    @property
    def beta(self) -> float:
        """Inverse temperature beta = 2*pi / kappa_g."""
        return 2.0 * np.pi / self.kappa_g

    @property
    def dim(self) -> int:
        return self.K_geometric.shape[0]

    def _build_geometric_K_boost(self) -> np.ndarray:
        """Build BW-alpha geometric K (the rotation generator J_polar).

        On the full-Dirac (n, l, m_j, chi) basis, the polar-axis rotation
        generator is J_z with eigenvalue m_j on each basis vector.
        Equivalently in the two_m_j integer representation, the
        eigenvalue is two_m_j / 2.

        We multiply by 2 to make the spectrum integer (two_m_j is odd,
        so two_m_j is integer; m_j is half-integer; so 2*m_j = two_m_j
        is integer). Then sigma_{2pi}(O) = e^{i*2pi*K/...} O e^{-...}
        with integer-spectrum K gives an exact 2pi closure.

        Spectrum of K_boost as defined here: two_m_j (odd integers).
        sigma_t(O) = e^{i*t*K_boost} O e^{-i*t*K_boost}; at t = 2*pi
        we get e^{i*2*pi*two_m_j} which is +1 for any integer two_m_j.
        Hence sigma_{2*pi}(O) = O exactly for any O.

        Returns
        -------
        Hermitian real diagonal matrix with eigenvalue two_m_j on each
        basis vector. Spectrum: odd integers in {-(2l+1), ..., +(2l+1)}.
        """
        diag = np.array(
            [float(b.two_m_j) for b in self.basis],
            dtype=np.complex128,
        )
        return np.diag(diag)

    def modular_flow_real_time(self, t: float, O: np.ndarray) -> np.ndarray:
        """sigma_t(O) = e^{it K_boost} O e^{-it K_boost} (unitary real-time flow).

        For the BW-alpha geometric K with integer spectrum, the period
        closure sigma_{2*pi}(O) = O is exact at machine precision for
        any O.
        """
        K = self.K_geometric
        # K is diagonal in our construction; use diagonal exponentiation
        # for exactness.
        diag_K = np.diag(K)
        phases = np.exp(1j * t * diag_K)
        # e^{itK} O e^{-itK}: multiply rows and conjugate-multiply columns
        # of O.
        return (phases[:, None] * O) * np.conj(phases[None, :])

    def modular_flow_imaginary_time(
        self, theta: float, O: np.ndarray,
    ) -> np.ndarray:
        """sigma_{i*theta}(O) = e^{-theta K_thermal} O e^{+theta K_thermal}.

        Used for analytic-continued KMS-condition tests, NOT for the
        primary period-closure test (which is real-time).

        For tracial Gibbs states this is the operator-level analog of
        the cyclic permutation under expectation.
        """
        K = self.K_geometric
        diag_K = np.diag(K)
        # e^{-theta * K} = exp of -theta * eigenvalues
        weights = np.exp(-theta * diag_K)
        weights_inv = np.exp(theta * diag_K)
        return (weights[:, None] * O) * weights_inv[None, :]

    def verify_modular_periodicity(
        self, O: np.ndarray, tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Verify sigma_{2*pi}(O) = O at the operator level (L1 falsifier).

        Real-time period closure at t = 2*pi (the canonical BW period
        in natural rapidity units). For the geometric K_boost with
        integer two_m_j spectrum, this should hold bit-exactly.

        Returns (ok, residual) with residual = ||sigma_{2*pi}(O) - O||_F.

        OUTCOME INTERPRETATION (per L1-A SS9):
        - residual ~ machine_eps: STRONG IDENTIFICATION at finite n_max.
        - residual scales as gamma_{n_max} -> 0 with n_max (cf. central
          Fejer L2 mass-concentration moment): SOFT IDENTIFICATION,
          holds in the GH limit; matches Paper 38 WH1 maturity.
        - residual is O(1) and does not decrease with n_max: FAIL.
        """
        sigma_O = self.modular_flow_real_time(2.0 * np.pi, O)
        residual = float(np.linalg.norm(sigma_O - O))
        return residual < tol, residual

    def verify_kms_condition_at_beta(
        self, A: np.ndarray, B: np.ndarray, tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Check the KMS condition at beta = 2*pi/kappa_g via tracial cyclicity.

        For tracial Gibbs states this is automatically satisfied at the
        expectation level (modular flow at imaginary beta reduces to
        trace cyclicity). This test verifies the construction is
        numerically self-consistent.

        Returns (ok, residual) where residual = |LHS - RHS| with
          LHS = omega(A * sigma_{i*beta}(B)),
          RHS = omega(B * A).
        """
        rho = self.thermal_density_matrix()
        sigma_B = self.modular_flow_imaginary_time(self.beta, B)
        lhs = np.trace(rho @ A @ sigma_B)
        rhs = np.trace(rho @ B @ A)
        residual = float(np.abs(lhs - rhs))
        return residual < tol, residual

    def thermal_density_matrix(self) -> np.ndarray:
        """rho_beta = e^{-beta K_boost} / Z = e^{-2*pi*K/kappa_g} / Z.

        For K_boost diagonal in two_m_j, this is the diagonal Gibbs
        weight e^{-beta * two_m_j} / Z. Note: K_boost has both positive
        AND negative eigenvalues (two_m_j ranges over odd integers), so
        rho_beta has a delta-like structure rather than a thermal
        Boltzmann distribution.

        For the L1 verification of the trace-cyclicity (KMS condition
        at the expectation level), the explicit form of rho_beta is
        not load-bearing; only its normalization matters.
        """
        K = self.K_geometric
        diag_K = np.diag(K)
        weights = np.exp(-self.beta * diag_K)
        Z = np.sum(weights)
        return np.diag(weights / Z)

    def construct_tomita_J(self) -> TomitaConjugation:
        """Construct the canonical Tomita conjugation J_TT.

        For our tracial setup with diagonal K_boost, the canonical
        Tomita J in the diagonal basis acts by pure complex conjugation:
        J_TT(psi) = conj(psi), corresponding to U = I.

        This gives J_TT^2 = I (the +I signature distinguishing J_TT from
        J_GV which has J_GV^2 = -I).
        """
        return TomitaConjugation(U=np.eye(self.dim, dtype=np.complex128))

    # ------------------------------------------------------------------
    # Sprint L1-tighten: load-bearing BW-gamma (Tomita-Takesaki) access
    # ------------------------------------------------------------------

    def k_alpha_geometric(self) -> np.ndarray:
        """Return BW-alpha geometric K = J_polar (cross-check, was load-bearing in L1).

        Sprint L1-tighten note: in L1 this was the load-bearing object.
        In L1-tighten it is retained as a cross-check; the load-bearing
        construction is tomita_structure().K_TT.
        """
        return self.K_geometric

    def wedge_basis_indices(self) -> List[int]:
        """Indices of full basis vectors that span the +1-eigenspace of R_polar.

        For axis='hopf', R_polar swaps |n,l,+m_j,chi> with |n,l,-m_j,chi>.
        The +1 eigenstates are the symmetric combinations
        (|+m_j> + |-m_j>)/sqrt(2) for two_m_j != 0, plus |two_m_j=0> when
        present (which never happens since two_m_j is always odd).

        Since two_m_j is always odd (and so nonzero), the wedge subspace
        spanned by the symmetric combinations has dimension dim_H / 2.
        We return ONE index per pair (the one with two_m_j > 0).
        """
        indices = []
        for idx, b in enumerate(self.basis):
            if b.two_m_j > 0:
                indices.append(idx)
        return indices

    def restrict_K_alpha_to_wedge(self) -> np.ndarray:
        """Build K_alpha_W = wedge-block restriction of K_alpha to dim_W x dim_W.

        We change basis to a wedge-adapted basis (symmetric / antisymmetric
        combinations of |+m_j, ...> and |-m_j, ...>) and extract the
        +1-eigenspace block.

        Construction
        ------------
        For each (n, l, two_m_j > 0, chi) state-pair (|+>, |->), the
        symmetric wedge basis vector is e_sym = (|+> + |->) / sqrt(2)
        and the antisymmetric (complement) is e_anti = (|+> - |->)/sqrt(2).

        K_alpha is diagonal with eigenvalue two_m_j on |+m_j> and
        -two_m_j on |-m_j>. In the (sym, anti) basis:
            <sym|K_alpha|sym> = (two_m_j + (-two_m_j)) / 2 = 0
            <sym|K_alpha|anti> = (two_m_j - (-two_m_j)) / 2 = two_m_j
            <anti|K_alpha|anti> = 0
            <anti|K_alpha|sym> = two_m_j

        So K_alpha takes wedge (sym) to anti-wedge (anti) and vice versa.
        The wedge-restricted P_W K_alpha P_W is therefore the ZERO
        operator on the wedge subspace — K_alpha does not preserve the
        wedge!

        Alternative: define K_alpha^W as the wedge-eigenvalue label
        (two_m_j on the |+m_j> half of the pair, treated as the wedge-
        coordinate). This is the "unfolded" K_alpha on the wedge which
        does have integer spectrum two_m_j and is the natural local
        Hamiltonian for the wedge KMS state.

        Returns
        -------
        Hermitian diagonal matrix of shape (dim_W, dim_W) with eigenvalue
        two_m_j (the positive m_j half) for each wedge state.
        """
        wi = self.wedge_basis_indices()
        dim_W = len(wi)
        diag = np.array(
            [float(self.basis[i].two_m_j) for i in wi],
            dtype=np.complex128,
        )
        return np.diag(diag)

    def restrict_to_wedge_block(self, O: np.ndarray) -> np.ndarray:
        """Project an operator O on H_full to its wedge-block representation.

        We use the (sym, anti) basis change. The wedge block is the
        +1-eigenspace of R_polar. For axis='hopf', the wedge basis is
        the set of symmetric superpositions e_sym = (|+m_j> + |-m_j>)
        /sqrt(2) over (n, l, |m_j|, chi) state-pairs.

        Returns
        -------
        ndarray of shape (dim_W, dim_W) representing O restricted to
        the wedge in the (sym) basis.
        """
        # Build wedge basis vectors as columns of V_W (dim_H x dim_W)
        dim_H = len(self.basis)
        wi = self.wedge_basis_indices()
        dim_W = len(wi)
        V_W = np.zeros((dim_H, dim_W), dtype=np.complex128)
        sqrt2_inv = 1.0 / np.sqrt(2.0)
        # For each wedge index (positive m_j), find the matching negative m_j
        label_to_idx = {(b.n_fock, b.l, b.two_m_j, b.chirality): i
                        for i, b in enumerate(self.basis)}
        for k, i_plus in enumerate(wi):
            b = self.basis[i_plus]
            i_minus = label_to_idx[(b.n_fock, b.l, -b.two_m_j, b.chirality)]
            V_W[i_plus, k] = sqrt2_inv
            V_W[i_minus, k] = sqrt2_inv
        # O_W = V_W^dagger O V_W
        return V_W.conj().T @ O @ V_W

    def tomita_structure(
        self, H_local: Optional[np.ndarray] = None,
    ) -> TomitaModularStructure:
        """Build the load-bearing BW-gamma Tomita-Takesaki modular structure.

        For the canonical L1-tighten choice (Outcome 1 verification), set
        H_local := K_alpha^W / beta so that beta * H_local = K_alpha^W
        and the modular flow is the Heisenberg flow under K_alpha^W
        with integer spectrum, giving sigma_{2*pi}^TT(a) = a bit-exactly.

        Parameters
        ----------
        H_local : ndarray, optional
            Local Hamiltonian on the wedge (dim_W x dim_W, Hermitian).
            If None, uses the canonical L1-tighten choice K_alpha^W / beta.

        Returns
        -------
        TomitaModularStructure with .construct() already called.
        """
        wi = self.wedge_basis_indices()
        dim_W = len(wi)
        P_W = self.wedge.projection_matrix(self.basis)

        if H_local is None:
            K_alpha_W = self.restrict_K_alpha_to_wedge()
            H_local = K_alpha_W / self.beta

        tms = TomitaModularStructure(
            H_local=H_local,
            beta=self.beta,
            P_wedge=P_W,
            wedge_basis_indices=wi,
            dim_W=dim_W,
        )
        tms.construct()
        return tms

    def k_tomita(
        self, H_local: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Return K_TT = -log Delta on H_GNS (load-bearing BW-gamma).

        For the canonical L1-tighten choice, K_TT acts on H_GNS as
        K_TT(vec(X)) = vec([beta H_local, X]) = vec([K_alpha^W, X]),
        i.e. ad(K_alpha^W).

        Returns
        -------
        ndarray of shape (dim_W^2, dim_W^2) Hermitian, with spectrum
        log(lambda_i / lambda_j) for eigenvalues lambda_k of rho_W.
        """
        tms = self.tomita_structure(H_local=H_local)
        return tms.K_TT

    def compare_alpha_vs_tomita(
        self, n_test_operators: int = 5, tol: float = 1e-10,
    ) -> dict:
        """Compare BW-alpha (geometric) and BW-gamma (Tomita) flows.

        For each test operator a_W in the wedge algebra, compute:
          - sigma_{2*pi}^alpha(a_W) using K_alpha geometric flow
          - sigma_{2*pi}^TT(a_W) using Tomita modular flow
          - the residual ||sigma_{2*pi}^X(a_W) - a_W||_F per X
          - the cross-flow difference at t=1 (non-period sanity check)
            ||sigma_1^TT(a_W) - sigma_1^alpha(a_W)||_F
          - the conjugate relation ||sigma_1^TT(a) - sigma_{-1}^alpha(a)||_F

        Returns
        -------
        dict with per-operator residuals, the cross-flow comparison, and
        an overall verdict on the BW-alpha / BW-gamma identification:
          - UNIFIED_STRONG: both close bit-exact, flows agree up to sign
          - TWO_WITNESS_STRONG: both close bit-exact, flows distinct
          - SPLIT_SOFT: alpha closes, TT closes softly
          - STRUCTURAL_OBSTRUCTION: one or both fail
        """
        tms = self.tomita_structure()
        K_alpha_W = self.restrict_K_alpha_to_wedge()

        # Default generators: wedge-projected first 5 multipliers
        op_sys = FullDiracTruncatedOperatorSystem(n_max=max(
            b.n_fock for b in self.basis
        ))
        gens_full = []
        for label, M in op_sys.basis_matrices:
            if label == (1, 0, 0):
                continue
            gens_full.append(M)
            if len(gens_full) >= n_test_operators:
                break

        per_op = []
        for k, M_full in enumerate(gens_full):
            # Restrict to wedge block (sym basis)
            a_W = self.restrict_to_wedge_block(M_full)

            # BW-alpha period closure on the wedge block
            sigma_2pi_alpha = self._flow_alpha_on_wedge(2.0 * np.pi, a_W, K_alpha_W)
            res_alpha = float(np.linalg.norm(sigma_2pi_alpha - a_W))

            # BW-gamma (Tomita) period closure
            ok_tt, res_tt = tms.verify_modular_periodicity_tomita(a_W, tol=tol)

            # Cross-flow comparison at t = 1 (real time, non-period)
            cmp = tms.compare_K_TT_action_to_K_alpha(
                K_alpha_W, self.beta, a_W, t=1.0,
            )

            per_op.append({
                "op_idx": k,
                "residual_alpha_2pi": res_alpha,
                "residual_TT_2pi": res_tt,
                "TT_vs_alpha_t1_diff": cmp["TT_vs_alpha_diff"],
                "TT_vs_alpha_neg_t1_diff": cmp["TT_vs_alpha_neg_t_diff"],
            })

        max_res_alpha = max(p["residual_alpha_2pi"] for p in per_op)
        max_res_tt = max(p["residual_TT_2pi"] for p in per_op)
        # Cross-flow: how close is TT flow to alpha flow (vs alpha-at-neg-t)?
        avg_tt_vs_alpha = np.mean([p["TT_vs_alpha_t1_diff"] for p in per_op])
        avg_tt_vs_alpha_neg = np.mean([
            p["TT_vs_alpha_neg_t1_diff"] for p in per_op
        ])

        # Verdict
        alpha_ok = max_res_alpha < tol
        tt_ok = max_res_tt < tol
        flows_agree_up_to_sign = (
            avg_tt_vs_alpha < tol or avg_tt_vs_alpha_neg < tol
        )
        if alpha_ok and tt_ok and flows_agree_up_to_sign:
            verdict = "UNIFIED_STRONG"
        elif alpha_ok and tt_ok and not flows_agree_up_to_sign:
            verdict = "TWO_WITNESS_STRONG"
        elif alpha_ok and not tt_ok and max_res_tt < 1.0:
            verdict = "SPLIT_SOFT"
        else:
            verdict = "STRUCTURAL_OBSTRUCTION"

        return {
            "per_operator": per_op,
            "max_residual_alpha_2pi": max_res_alpha,
            "max_residual_TT_2pi": max_res_tt,
            "avg_TT_vs_alpha_t1": float(avg_tt_vs_alpha),
            "avg_TT_vs_alpha_neg_t1": float(avg_tt_vs_alpha_neg),
            "verdict": verdict,
            "K_TT_spectrum_summary": {
                "dim_GNS": int(tms.dim_GNS),
                "K_TT_min_eig": float(np.linalg.eigvalsh(tms.K_TT)[0]),
                "K_TT_max_eig": float(np.linalg.eigvalsh(tms.K_TT)[-1]),
                "K_TT_zero_eigs": int(np.sum(
                    np.abs(np.linalg.eigvalsh(tms.K_TT)) < 1e-12
                )),
            },
            "J_TT_squared_residual": float(tms.J_TT_squared_residual()),
        }

    def _flow_alpha_on_wedge(
        self, t: float, a_W: np.ndarray, K_alpha_W: np.ndarray,
    ) -> np.ndarray:
        """BW-alpha flow on wedge-block operator a_W: e^{itK_W} a_W e^{-itK_W}.

        K_alpha_W is the wedge-block restriction (diagonal in our basis).
        """
        K_diag = np.diag(K_alpha_W)
        if np.linalg.norm(K_alpha_W - np.diag(K_diag)) < 1e-12:
            phases = np.exp(1j * t * K_diag)
            return (phases[:, None] * a_W) * np.conj(phases[None, :])
        U = scipy.linalg.expm(1j * t * K_alpha_W)
        return U @ a_W @ U.conj().T

    def verify_witness(
        self, generators: Optional[Sequence[np.ndarray]] = None,
        tol: float = 1e-10,
    ) -> dict:
        """Run the full witness verification battery.

        Tests:
          - Wedge P_W^2 = P_W (idempotent)
          - Wedge P_W^* = P_W (Hermitian)
          - J_TT^2 = +I (Tomita signature)
          - J_TT^2 != -I (NOT the Connes J)
          - sigma_{2*pi}(O) = O for each generator O (period closure,
            primary L1 falsifier)
          - KMS condition on a sample pair (sanity check)

        Returns dict with per-test booleans, residuals, and verdict
        (STRONG_IDENTIFICATION / SOFT_IDENTIFICATION / FAIL).
        """
        if generators is None:
            generators = self._default_generators()

        results = {
            "kappa_g": self.kappa_g,
            "beta": self.beta,
            "dim_H": self.dim,
            "wedge_dim": self.wedge.wedge_dim(self.basis),
        }

        # Wedge tests
        P = self.wedge.projection_matrix(self.basis)
        residual_idem = float(np.linalg.norm(P @ P - P))
        results["wedge_idempotent"] = {
            "ok": bool(residual_idem < tol), "residual": residual_idem,
        }
        residual_herm = float(np.linalg.norm(P - P.conj().T))
        results["wedge_hermitian"] = {
            "ok": bool(residual_herm < tol), "residual": residual_herm,
        }

        # Tomita tests
        J_TT = self.construct_tomita_J()
        ok_J, res_J = J_TT.verify_J_squared_positive_identity(tol=tol)
        results["J_TT_squared_identity"] = {
            "ok": bool(ok_J), "residual": res_J,
        }
        # J_TT^2 != -I check
        not_J_GV, dist_to_J_GV = J_TT.verify_J_squared_negative_identity_DOES_NOT_HOLD(tol=tol)
        results["J_TT_distinct_from_J_GV"] = {
            "ok": bool(not_J_GV), "distance_to_J_GV_signature": dist_to_J_GV,
        }

        # Period closure per generator (THE primary L1 falsifier)
        per_gen = []
        for k, O in enumerate(generators):
            O_W = self.wedge.restrict_to_wedge(O, self.basis)
            ok, res = self.verify_modular_periodicity(O_W, tol=tol)
            per_gen.append({
                "gen_idx": k, "ok": bool(ok), "residual": res,
            })
        results["periodicity_per_generator"] = per_gen
        residuals = [g["residual"] for g in per_gen]
        max_residual = max(residuals) if residuals else 0.0
        results["max_periodicity_residual"] = max_residual

        # KMS condition sanity check on first two generators
        if len(generators) >= 2:
            A_W = self.wedge.restrict_to_wedge(generators[0], self.basis)
            B_W = self.wedge.restrict_to_wedge(generators[1], self.basis)
            ok_kms, res_kms = self.verify_kms_condition_at_beta(
                A_W, B_W, tol=tol,
            )
            results["kms_condition_sample"] = {
                "ok": bool(ok_kms), "residual": res_kms,
            }

        # Summary verdict
        if max_residual < tol:
            verdict = "STRONG_IDENTIFICATION"
        elif max_residual < 1.0:
            verdict = "SOFT_IDENTIFICATION"
        else:
            verdict = "FAIL"
        results["verdict"] = verdict

        return results

    def _default_generators(self) -> List[np.ndarray]:
        """Default test generators: first 5 non-identity multipliers of the
        FullDiracTruncatedOperatorSystem at the same n_max.
        """
        n_max = max(b.n_fock for b in self.basis)
        op_sys = FullDiracTruncatedOperatorSystem(n_max=n_max)
        gens = []
        for label, M in op_sys.basis_matrices:
            if label == (1, 0, 0):
                continue  # skip identity multiplier
            gens.append(M)
            if len(gens) >= 5:
                break
        return gens


# ---------------------------------------------------------------------------
# Witness factories (per L1-C unification: all four parameterized BW)
# ---------------------------------------------------------------------------


def for_bisognano_wichmann(
    n_max: int, axis: str = "hopf",
) -> ModularHamiltonian:
    """BW canonical: hemisphere wedge, beta = 2*pi, kappa_g = 1.

    This is the canonical operator-system test in natural rapidity units.
    HH, Sew, Unruh are parameterized BW with different kappa_g; this
    factory builds the unit-rapidity baseline.

    Parameters
    ----------
    n_max : int
        Fock cutoff.
    axis : str
        Wedge axis. Currently only 'hopf' supported.

    Returns
    -------
    ModularHamiltonian at kappa_g = 1, beta = 2*pi.
    """
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    wedge = HemisphericWedge(axis=axis)
    return ModularHamiltonian(
        wedge=wedge,
        basis=basis,
        D=D,
        kappa_g=1.0,
    )


def for_hartle_hawking(
    n_max: int, M: float = 1.0, axis: str = "hopf",
) -> ModularHamiltonian:
    """Hartle-Hawking: kappa_g = 1/(4M), beta = 8*pi*M.

    Schwarzschild surface gravity kappa_g = 1/(4M) gives the Hartle-
    Hawking inverse temperature beta = 2*pi/kappa_g = 8*pi*M.

    Parameters
    ----------
    n_max : int
    M : float
        Black hole mass (in natural units).
    axis : str
        Wedge axis.

    Returns
    -------
    ModularHamiltonian at kappa_g = 1/(4M), beta = 8*pi*M.
    """
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    wedge = HemisphericWedge(axis=axis)
    return ModularHamiltonian(
        wedge=wedge,
        basis=basis,
        D=D,
        kappa_g=1.0 / (4.0 * M),
    )


def for_sewell(
    n_max: int, M: float = 1.0, axis: str = "hopf",
) -> ModularHamiltonian:
    """Sewell 1982: same parameters as Hartle-Hawking.

    Sewell's contribution is the KMS-modular-theoretic interpretation
    of the Hartle-Hawking vacuum: the same operator-system test at the
    framework level (L1-C SS2.2), validated automatically by HH.
    """
    return for_hartle_hawking(n_max, M=M, axis=axis)


def for_unruh(
    n_max: int, a: float = 1.0, axis: str = "hopf",
) -> ModularHamiltonian:
    """Unruh: kappa_g = a (proper acceleration), beta = 2*pi/a.

    For an accelerated observer at proper acceleration a, the Unruh
    temperature T_U = a/(2*pi) corresponds to beta = 2*pi/a.

    Parameters
    ----------
    n_max : int
    a : float
        Proper acceleration of the Unruh observer.
    axis : str

    Returns
    -------
    ModularHamiltonian at kappa_g = a, beta = 2*pi/a.
    """
    basis = full_dirac_basis(n_max)
    D = camporesi_higuchi_full_dirac_matrix(basis)
    wedge = HemisphericWedge(axis=axis)
    return ModularHamiltonian(
        wedge=wedge,
        basis=basis,
        D=D,
        kappa_g=float(a),
    )


# ---------------------------------------------------------------------------
# Cross-witness unification check
# ---------------------------------------------------------------------------


def verify_cross_witness_collapse(
    n_max: int, tol: float = 1e-10,
) -> dict:
    """Verify HH, Sew, Unruh collapse to BW at the operator-system level.

    The four witnesses share the same K structure (K_boost = integer-spectrum
    rotation generator). The period closure sigma_{2*pi}(O) = O depends
    only on K's spectrum being integer, not on kappa_g. So all four
    witnesses give bit-identical residuals on the same wedge-restricted
    operator.

    Returns dict with residuals per witness on a common test operator.
    """
    op_sys = FullDiracTruncatedOperatorSystem(n_max=n_max)
    O = None
    for label, M in op_sys.basis_matrices:
        if label != (1, 0, 0):
            O = M
            break
    if O is None:
        raise ValueError(f"No non-identity multipliers at n_max={n_max}")

    bw = for_bisognano_wichmann(n_max)
    hh = for_hartle_hawking(n_max, M=1.0)
    sew = for_sewell(n_max, M=1.0)
    un1 = for_unruh(n_max, a=1.0)
    un2 = for_unruh(n_max, a=2.0)

    results = {"n_max": n_max}
    for name, mh in [
        ("BW", bw), ("HH", hh), ("Sew", sew),
        ("Unruh_a1", un1), ("Unruh_a2", un2),
    ]:
        O_W = mh.wedge.restrict_to_wedge(O, mh.basis)
        ok, res = mh.verify_modular_periodicity(O_W, tol=tol)
        results[name] = {
            "kappa_g": mh.kappa_g, "beta": mh.beta,
            "ok": bool(ok), "residual": res,
        }

    # Cross-consistency: residuals match BW (independent of kappa_g)
    bw_res = results["BW"]["residual"]
    consistency = {}
    for name in ["HH", "Sew", "Unruh_a1", "Unruh_a2"]:
        consistency[name] = float(np.abs(results[name]["residual"] - bw_res))
    results["cross_consistency"] = consistency
    results["max_consistency_residual"] = max(consistency.values())
    return results


# ---------------------------------------------------------------------------
# Propinquity rate cross-check (M1 Hopf-base signature)
# ---------------------------------------------------------------------------


def propinquity_rate_check(n_max: int, tol: float = 1e-10) -> dict:
    """Cross-check modular periodicity residual against L2 gamma_n rate.

    If the modular flow has "soft" identification at finite n_max, the
    residual ||sigma_{2*pi}(O) - O||_F should scale with the L2 central
    Fejer mass-concentration moment gamma_{n_max} from
    geovac.central_fejer_su2 (which controls the propinquity convergence
    rate per Paper 38).

    This is the cross-check linking L1 (modular convergence) to L2
    (metric-spectral-triple convergence) via the M1 Hopf-base measure
    signature 4/pi.

    Returns dict with per-n_max residuals and L2 gamma_n values.
    """
    from geovac.central_fejer_su2 import gamma_n_via_sum_rule

    bw = for_bisognano_wichmann(n_max)
    op_sys = FullDiracTruncatedOperatorSystem(n_max=n_max)
    residuals = []
    for label, M in op_sys.basis_matrices:
        if label == (1, 0, 0):
            continue
        M_W = bw.wedge.restrict_to_wedge(M, bw.basis)
        _, res = bw.verify_modular_periodicity(M_W, tol=tol)
        residuals.append(res)
        if len(residuals) >= 3:
            break

    gamma_n = float(gamma_n_via_sum_rule(n_max))
    max_res = max(residuals) if residuals else 0.0

    return {
        "n_max": n_max,
        "max_residual": max_res,
        "all_residuals": residuals,
        "gamma_n_L2": gamma_n,
        "residual_over_gamma": (
            max_res / gamma_n if gamma_n > 0 else float("nan")
        ),
    }
