"""Krein-level modular Hamiltonian on the hemispheric wedge of S^3 x R.

Sprint L2-E (2026-05-16): re-derive the Paper 42 unified-strong four-witness
Wick-rotation theorem at the Krein level on a hemispheric wedge of
S^3 x R at signature (3, 1).  This is the Lorentzian-side companion of
geovac.modular_hamiltonian (which closes the Riemannian-side BW-alpha
and BW-gamma constructions of Paper 42 at finite n_max).

Architecture (locked decisions from Sprint L2-A scoping memo Section 4.7,
5.4 and Sprint L2-F falsifier catalogue Section 3 L2E-FALS-1/2/3):

  Wedge W_L on the Krein space K = H_GV (x) C^{N_t}:

      W_L  =  P_W_spatial  (x)  P_t_positive

  with P_W_spatial = (1/2)(I + R_polar) the Paper 42 hemispheric wedge
  on H_GV (m_j -> -m_j reflection), and P_t_positive the projector onto
  positive-t grid points.  At N_t = 1 (Riemannian limit), the temporal
  grid is the singleton t = 0.0 and P_t_positive is taken to be the
  full identity on C^1, recovering Paper 42's spatial wedge exactly.

  At general N_t > 1, the symmetric grid t = linspace(-T_max, T_max, N_t)
  has half the grid points in t > 0; for N_t odd the central point t = 0
  is included in the wedge (with weight 1) and the structural results
  are independent of this convention at the period-closure level.

  Boost generator K_L_alpha:

      K_L_alpha  =  K_alpha^W (x) I_{N_t}

  where K_alpha^W is the Paper 42 BW-alpha geometric generator
  (diag(two_m_j) on the wedge), and the temporal slot acts as identity.
  This is the natural Killing-vector-preserving choice: the spatial
  rotation around the Hopf-base axis preserves the wedge and is
  insensitive to the temporal slot (the Lorentzian boost acts on the
  spacetime wedge via the spatial rotation generator at the Riemannian
  level after Wick rotation).  Integer spectrum is preserved.

  Wedge KMS state at beta = 2*pi (BW canonical):

      rho_W^L  =  e^{- K_L_alpha^W} / Z

  via the BW choice H_local := K_L_alpha^W / beta (Paper 42 Section 4.2
  scope finding O3 lifted to (3, 1)).  This is beta-independent at the
  algebra-action level (same property as Riemannian side), so the
  six-witness collapse is automatic.

  Krein structural conventions per Sprint L2-D verdict:

    R1 (truthful CH): primary path, preserves Riemannian-limit recovery
        (load-bearing per L2-C).  Used for BW-alpha (which only requires
        JD = +DJ via Connes axiom (iv); satisfies bit-exact at all
        tested cutoffs per L2-D).

    R2 (offdiag CH): alternative path, satisfies the BBB universal
        {chi, D} = 0 axiom.  Not used at the wedge-construction level
        because the wedge structure is driven by the polar reflection,
        not by D, and the four-witness theorem closes on the integer
        K_alpha spectrum regardless of D choice.

The construction is the LOAD-BEARING falsifier of Sprint L2-E (L2E-FALS-1
and L2E-FALS-2 per debug/sprint_l2_falsifiers.md):

    sigma_{2*pi}^{L, alpha}(O)  =  O          (BW-alpha period closure)
    sigma_{2*pi}^{L, TT}(O)     =  O          (BW-gamma Tomita closure)
    sigma_t^{L, TT}(O)          =  sigma_{-t}^{L, alpha}(O)   (conjugacy)

and the LOAD-BEARING falsifier #4 of the sprint:

    K_L_alpha | _{N_t = 1}  =  K_alpha       (Riemannian-limit recovery)
    K_L_TT    | _{N_t = 1}  =  K_TT          (Tomita Riemannian limit)

bit-identically at every n_max in {1, 2, 3} x N_t in {1, 11, 21}.

CRITICAL CAUTION (per Sprint L2-D structural finding §5):

  In GeoVac at (3, 1), the BBB universal axiom {chi, D} = 0 FAILS on
  truthful Camporesi-Higuchi D_GV.  This is a load-bearing scope
  finding, NOT a basis bug.  For the L2-E modular Hamiltonian
  construction, this finding is bypassed because:

    (a) BW-alpha uses K_alpha^W = diag(two_m_j), which is a kinematic
        wedge generator (not the Dirac).  Period closure is
        D-independent at the operator-action level: it depends only on
        the integer spectrum of K_alpha^W.

    (b) BW-gamma uses Delta = e^{-K_alpha^W} (left) tensored with
        e^{+K_alpha^W} (right) on the GNS Hilbert-Schmidt space.  This
        is also D-independent because rho_W = e^{-K_alpha^W}/Z under
        the BW choice H_local = K_alpha^W / beta.

  Therefore both Lorentzian constructions inherit the truthful-CH
  Riemannian-limit recovery (load-bearing per L2-C) without needing
  the BBB universal chi-D anticommutation.  Documented in L2-D §5.4
  resolution R1 (accept structural finding, proceed).

  H_LOCAL VERDICT TO TEST (Paper 42 §7.2 open question O3 at (3, 1)):

  Compute K_L_alpha^W / beta at beta = 2*pi (using the Lorentzian
  wedge construction) and compare to the wedge-restricted truthful
  D_L_W (the Lorentzian Dirac from L2-C, restricted to W_L).  Three
  possible verdicts:

    (i)  K_L_alpha^W / beta = D_L_W at (3, 1): H_local choice IS the
         Lorentzian Dirac at signature (3, 1) but not at (3, 0).
         Refines Paper 42 §7.2 cleanly.
    (ii) K_L_alpha^W / beta != D_L_W at (3, 1), same as Riemannian:
         O3 holds universally; the spectral-action-vs-modular-
         Hamiltonian generator distinction is signature-independent.
    (iii) Some intermediate: document the structural relationship.

  This is the most structurally important question this sprint can
  answer.

API
===

  LorentzianWedge(krein)                  # P_W_L spatial x positive-t
  build_K_L_alpha(krein, wedge)           # geometric boost generator
  LorentzianTomitaStructure(krein, wedge) # BW-gamma on Krein-GNS
  LorentzianModularHamiltonian(krein)     # full driver class

  # Six-witness factories (parameterized BW per Paper 42 §3.3 unification)
  for_bisognano_wichmann_lorentzian(n_max, N_t, T_max)
  for_hartle_hawking_lorentzian(n_max, N_t, T_max, M=1.0)
  for_sewell_lorentzian(n_max, N_t, T_max, M=1.0)
  for_unruh_lorentzian(n_max, N_t, T_max, a=1.0)

  # Six-witness collapse verification
  verify_cross_witness_collapse_lorentzian(n_max, N_t, T_max)

References
==========

  van den Dungen, K.  "Krein spectral triples and the fermionic action."
    Math. Phys. Anal. Geom. 19, 4 (2016).  arXiv:1505.01939.

  Bizi, N., Brouder, C., Besnard, F.  "Space and time dimensions of
    algebras..."  J. Math. Phys. 59, 062303 (2018).  arXiv:1611.07062.

  Paper 42 (Riemannian-side closure):
    papers/standalone/paper_42_modular_hamiltonian_four_witness.tex

  Sprint L2 architecture:
    debug/sprint_l2a_scoping_memo.md            Section 4.7, 5.4
    debug/sprint_l2_falsifiers.md              Section 3 L2E-*
    debug/l2_b_krein_construction_memo.md       Krein conventions
    debug/l2_c_lorentzian_dirac_memo.md         D_L construction
    debug/l2_d_connes_axiom_audit_31_memo.md    R1+R2 resolution

  GeoVac code:
    geovac/krein_space_construction.py    (L2-B Krein space)
    geovac/lorentzian_dirac.py            (L2-C Lorentzian Dirac)
    geovac/connes_axiom_audit_31.py       (L2-D BBB audit, J_L)
    geovac/modular_hamiltonian.py         (Riemannian Paper 42 template)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import numpy as np
import scipy.linalg

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    FullDiracTruncatedOperatorSystem,
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
)
from geovac.krein_space_construction import KreinSpace
from geovac.lorentzian_dirac import lorentzian_dirac_matrix
from geovac.modular_hamiltonian import (
    HemisphericWedge,
    ModularHamiltonian,
)


# ---------------------------------------------------------------------------
# Temporal half-line projector P_t_positive on L^2(R_t)_cutoff
# ---------------------------------------------------------------------------


def temporal_positive_half_projector(N_t: int) -> np.ndarray:
    """Projector onto positive-t grid points on a symmetric uniform grid.

    Convention on the symmetric grid t_k = linspace(-T_max, +T_max, N_t):

      - At N_t = 1, the singleton grid is t = 0.0 and the projector is
        the identity on C^1 (Riemannian-limit reduction).
      - At N_t > 1, the grid is symmetric around 0.  For N_t odd the
        center index is included with weight 1; for N_t even no center
        exists.  The convention is t >= 0 inclusive: indices k with
        t_k >= 0 get weight 1; indices with t_k < 0 get weight 0.

    Choice rationale (lock-in): t >= 0 inclusive is the natural
    Rindler-wedge analog (the wedge boundary at t = 0 is included in
    the wedge closure).  At the period-closure level the answer is
    insensitive to whether t = 0 is in the wedge or not, because the
    boost generator K_alpha^W is diagonal in the temporal slot via the
    identity factor I_{N_t}.

    Parameters
    ----------
    N_t : int
        Number of temporal grid points. Must be >= 1.

    Returns
    -------
    P_t : np.ndarray of shape (N_t, N_t), complex128
        Diagonal projector with entries in {0, 1}.

    Raises
    ------
    ValueError if N_t < 1.
    """
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    if N_t == 1:
        # Riemannian limit: singleton t = 0 grid, include in wedge.
        return np.eye(1, dtype=np.complex128)
    # Symmetric grid t = linspace(-T_max, T_max, N_t); P_t selects t >= 0.
    t = np.linspace(-1.0, 1.0, N_t)  # T_max-scale irrelevant for sign
    diag = np.array([1.0 if t_k >= 0 else 0.0 for t_k in t], dtype=np.complex128)
    return np.diag(diag)


# ---------------------------------------------------------------------------
# Lorentzian wedge W_L on K = H_GV (x) C^{N_t}
# ---------------------------------------------------------------------------


@dataclass
class LorentzianWedge:
    """Hemispheric Lorentzian wedge W_L on the Krein space.

    Construction:

      P_W_L  =  P_W_spatial  (x)  P_t_positive

    where:
      - P_W_spatial = (1/2)(I + R_polar) is the Paper 42 hemispheric
        wedge on H_GV (m_j -> -m_j reflection on the spinor basis).
      - P_t_positive selects positive-t grid points (Rindler-wedge analog
        in the temporal direction); at N_t = 1 reduces to identity.

    At N_t = 1 (Riemannian limit) the construction reduces to Paper 42's
    P_W = (1/2)(I + R_polar) on H_GV bit-identically.

    Properties:
      - P_W_L is idempotent and Hermitian (both factors are).
      - dim(W_L) = dim_W_spatial * count(t_k >= 0).

    Attributes
    ----------
    krein : KreinSpace
        The Krein space from Sprint L2-B at the desired (n_max, N_t, T_max).
    P_W_spatial : np.ndarray of shape (dim_spatial, dim_spatial), complex128
        Spatial wedge projector from Paper 42.
    P_t_positive : np.ndarray of shape (N_t, N_t), complex128
        Temporal half-line projector.
    P_W_L : np.ndarray of shape (dim_K, dim_K), complex128
        Full Krein wedge projector P_W_spatial (x) P_t_positive.
    """

    krein: KreinSpace
    P_W_spatial: np.ndarray = field(init=False)
    P_t_positive: np.ndarray = field(init=False)
    P_W_L: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        spatial_wedge = HemisphericWedge(axis="hopf")
        self.P_W_spatial = spatial_wedge.projection_matrix(
            self.krein.basis_spatial
        )
        self.P_t_positive = temporal_positive_half_projector(self.krein.N_t)
        self.P_W_L = np.kron(self.P_W_spatial, self.P_t_positive)

    @property
    def dim_W_spatial(self) -> int:
        """Spatial wedge dimension = dim_H / 2."""
        return self.krein.dim_spatial // 2

    @property
    def dim_W_L(self) -> int:
        """Lorentzian wedge dimension (trace of P_W_L)."""
        return int(np.round(np.real(np.trace(self.P_W_L))))

    @property
    def N_t_positive(self) -> int:
        """Number of positive-t grid points (= dim_W_L / dim_W_spatial)."""
        return int(np.round(np.real(np.trace(self.P_t_positive))))

    def verify_idempotent(self, tol: float = 1e-12) -> Tuple[bool, float]:
        """P_W_L^2 = P_W_L."""
        residual = float(np.linalg.norm(self.P_W_L @ self.P_W_L - self.P_W_L))
        return residual < tol, residual

    def verify_hermitian(self, tol: float = 1e-12) -> Tuple[bool, float]:
        """P_W_L^* = P_W_L."""
        residual = float(np.linalg.norm(self.P_W_L - self.P_W_L.conj().T))
        return residual < tol, residual

    def restrict_to_wedge(self, O: np.ndarray) -> np.ndarray:
        """Compute the wedge-restricted operator P_W_L O P_W_L on K."""
        return self.P_W_L @ O @ self.P_W_L

    def wedge_basis_indices(self) -> List[int]:
        """Indices in the full Krein basis spanning the +1-eigenspace of R_polar
        in the spatial slot AND positive-t in the temporal slot.

        The full Krein basis is (spatial_index, t_index) with stride N_t
        (kron(spatial, temporal) convention).  We return indices of
        states with two_m_j > 0 AND t_k >= 0.
        """
        N_t = self.krein.N_t
        # Spatial wedge: same logic as ModularHamiltonian.wedge_basis_indices
        spatial_wedge_indices = []
        for k, b in enumerate(self.krein.basis_spatial):
            if b.two_m_j > 0:
                spatial_wedge_indices.append(k)
        # Temporal wedge: t_k >= 0 indices
        if N_t == 1:
            temporal_wedge_indices = [0]
        else:
            t = np.linspace(-1.0, 1.0, N_t)
            temporal_wedge_indices = [k for k, tk in enumerate(t) if tk >= 0]
        # Combine
        out = []
        for s_idx in spatial_wedge_indices:
            for t_idx in temporal_wedge_indices:
                # kron(spatial, temporal) ordering: full_idx = s*N_t + t
                out.append(s_idx * N_t + t_idx)
        return out


# ---------------------------------------------------------------------------
# Geometric Lorentzian boost generator K_L_alpha on the wedge
# ---------------------------------------------------------------------------


def build_K_L_alpha_full(krein: KreinSpace) -> np.ndarray:
    """Build the geometric Lorentzian boost generator on the full Krein space.

    Construction (per Sprint L2-A scoping memo Section 4.7):

        K_L_alpha  =  K_alpha_spatial  (x)  I_{N_t}

    where K_alpha_spatial = diag(two_m_j) on H_GV is the Paper 42 BW-alpha
    geometric generator (rotation around the Hopf-base axis, integer
    spectrum two_m_j on the spinor basis).  The temporal slot acts as
    identity because the boost (Wick-rotated rotation) is spatial.

    Properties:
      - K_L_alpha is Hermitian (real diagonal).
      - Spec(K_L_alpha) is contained in odd integers (inherited from
        spatial K_alpha) times N_t multiplicity.

    At N_t = 1 (Riemannian limit), this reduces to K_alpha bit-identically.

    Parameters
    ----------
    krein : KreinSpace
        Krein space from L2-B.

    Returns
    -------
    K_L_alpha : np.ndarray of shape (dim_K, dim_K), complex128 Hermitian diagonal.
    """
    # Spatial K_alpha = diag(two_m_j) on the H_GV basis
    diag_spatial = np.array(
        [float(b.two_m_j) for b in krein.basis_spatial],
        dtype=np.complex128,
    )
    K_alpha_spatial = np.diag(diag_spatial)
    I_t = np.eye(krein.N_t, dtype=np.complex128)
    return np.kron(K_alpha_spatial, I_t)


def restrict_K_L_alpha_to_wedge(
    krein: KreinSpace, wedge: LorentzianWedge,
) -> np.ndarray:
    """Build K_L_alpha^W on the Lorentzian wedge sub-Hilbert space.

    Like Paper 42 §5.1, we extract the "unfolded" K_alpha on the wedge:
    for each (n, l, |m_j|, chi, t_k) state with t_k >= 0, the eigenvalue
    is +two_m_j (the positive half of the two_m_j-pair) and the diagonal
    matrix on this wedge is built explicitly.

    The result is a Hermitian diagonal matrix of shape
    (dim_W_L, dim_W_L) with integer eigenvalues.

    Parameters
    ----------
    krein : KreinSpace
    wedge : LorentzianWedge

    Returns
    -------
    K_L_alpha_W : np.ndarray of shape (dim_W_L, dim_W_L), Hermitian diagonal.
    """
    wedge_indices = wedge.wedge_basis_indices()
    N_t = krein.N_t
    dim_spatial = krein.dim_spatial
    diag = []
    for idx in wedge_indices:
        # idx = spatial_idx * N_t + t_idx (kron convention)
        spatial_idx = idx // N_t
        # two_m_j is the spatial wedge eigenvalue (positive)
        label = krein.basis_spatial[spatial_idx]
        diag.append(float(label.two_m_j))
    return np.diag(np.array(diag, dtype=np.complex128))


def restrict_operator_to_wedge_block(
    O_full: np.ndarray, krein: KreinSpace, wedge: LorentzianWedge,
) -> np.ndarray:
    """Project a full-Krein operator O to its wedge-block representation.

    Uses the wedge-adapted basis change: for each (n, l, two_m_j > 0,
    chi, t_k >= 0), the symmetric wedge basis vector is
    e_sym = (|+m_j, t_k> + |-m_j, t_k>) / sqrt(2).  Then
    O_W = V_W^dagger O V_W in this basis.

    Parameters
    ----------
    O_full : np.ndarray of shape (dim_K, dim_K)
    krein : KreinSpace
    wedge : LorentzianWedge

    Returns
    -------
    O_W : np.ndarray of shape (dim_W_L, dim_W_L)
    """
    N_t = krein.N_t
    dim_K = krein.dim
    wi = wedge.wedge_basis_indices()
    dim_W_L = len(wi)
    V_W = np.zeros((dim_K, dim_W_L), dtype=np.complex128)
    sqrt2_inv = 1.0 / np.sqrt(2.0)
    # Map (spatial_basis_label) -> full spatial index for finding the partner
    label_to_idx = {
        (b.n_fock, b.l, b.two_m_j, b.chirality): i
        for i, b in enumerate(krein.basis_spatial)
    }
    for k_W, full_idx in enumerate(wi):
        # Decompose full_idx = spatial_idx * N_t + t_idx
        spatial_idx_plus = full_idx // N_t
        t_idx = full_idx % N_t
        b = krein.basis_spatial[spatial_idx_plus]
        # Find the partner with two_m_j -> -two_m_j (same n, l, chi, t_idx)
        spatial_idx_minus = label_to_idx[
            (b.n_fock, b.l, -b.two_m_j, b.chirality)
        ]
        full_idx_minus = spatial_idx_minus * N_t + t_idx
        V_W[full_idx, k_W] = sqrt2_inv
        V_W[full_idx_minus, k_W] = sqrt2_inv
    return V_W.conj().T @ O_full @ V_W


# ---------------------------------------------------------------------------
# Lorentzian Tomita-Takesaki modular structure on Krein-GNS space
# ---------------------------------------------------------------------------


@dataclass
class LorentzianTomitaStructure:
    """Tomita-Takesaki modular structure on the Krein wedge KMS state.

    Construction (Lorentzian-side mirror of Paper 42 §6):

    For the wedge sub-algebra A_W = B(K_W) on the Lorentzian wedge
    K_W = P_W_L K and the wedge KMS state

        rho_W^L  =  e^{- K_L_alpha^W} / Z

    under the BW choice H_local := K_L_alpha^W / beta (Paper 42 §4.2
    lifted to (3, 1)), the GNS Hilbert-Schmidt space is

        K_GNS  =  M_{dim_W_L}(C)  =~  C^{dim_W_L^2}

    with cyclic vector Omega = vec(rho^{1/2}).  The Tomita S operator
    is S: a Omega -> a^* Omega, antilinear on K_GNS.  Polar decomposition
    S = J_L_TT * Delta_L^{1/2} gives:

        Delta_L          (vec(X))   =   vec(rho X rho^{-1})
        J_L_TT           (vec(X))   =   vec(rho^{1/2} X^* rho^{-1/2})    (antilinear)
        K_L_TT           :=  - log Delta_L
        sigma_t^{L,TT}   (a)        =   rho^{it} a rho^{-it}.

    The modular flow on the algebra at the period is

        sigma_{2pi}^{L,TT}(a)
            =  e^{-i 2pi K_L_alpha^W} a e^{+i 2pi K_L_alpha^W}
            =  a

    bit-exactly because K_L_alpha^W has integer spectrum on the wedge.

    Lorentzian convention:  the Krein-positive cone reading.  Per
    van den Dungen 2016 §2, the Krein-positive completion of the GNS
    space is the natural setting for Tomita-Takesaki on a Krein
    spectral triple.  Here the wedge KMS state rho_W^L is defined as
    a Hilbert-space trace state (= Tr_{K_W}), so the natural GNS
    triple is the Hilbert-space GNS of the algebra acting on K_W; the
    "Krein-positive" framing matters only for the J_L_TT signature
    distinction (= +I always on positive-trace GNS).

    Attributes
    ----------
    H_local : np.ndarray
        Hermitian local Hamiltonian on the wedge (dim_W_L x dim_W_L).
    beta : float
        Inverse temperature (BW canonical: 2*pi).
    K_L_alpha_W : np.ndarray
        BW-alpha geometric generator on the wedge (diagonal, integer
        spectrum), satisfying K_L_alpha_W = beta * H_local.
    dim_W_L : int
        Dimension of the wedge subspace.

    Cached after construct():
        rho, rho_sqrt, rho_invsqrt, Delta_matrix, K_TT
    """

    H_local: np.ndarray
    beta: float
    K_L_alpha_W: np.ndarray
    dim_W_L: int

    rho: Optional[np.ndarray] = None
    rho_sqrt: Optional[np.ndarray] = None
    rho_invsqrt: Optional[np.ndarray] = None
    Delta_matrix: Optional[np.ndarray] = None
    K_TT: Optional[np.ndarray] = None

    @property
    def dim_GNS(self) -> int:
        """Dimension of the GNS Hilbert-Schmidt space = dim_W_L^2."""
        return self.dim_W_L ** 2

    def construct(self) -> None:
        """Build the explicit Tomita modular structure.

        Computes:
          - rho_W^L = e^{-beta H_local} / Z on the wedge
          - rho^{1/2}, rho^{-1/2} via diagonalization
          - Delta_L as a (dim_W_L^2 x dim_W_L^2) matrix on vec(X)
          - K_L_TT = -log Delta_L
        """
        rho_unnorm = scipy.linalg.expm(-self.beta * self.H_local)
        Z = np.trace(rho_unnorm).real
        if Z <= 0:
            raise ValueError(f"partition function nonpositive: Z = {Z}")
        self.rho = rho_unnorm / Z

        # rho^{1/2} and rho^{-1/2} via diagonalization
        eigvals, eigvecs = np.linalg.eigh(self.rho)
        eigvals = np.maximum(eigvals, 1e-30)
        sqrt_diag = np.sqrt(eigvals)
        invsqrt_diag = 1.0 / sqrt_diag
        self.rho_sqrt = (eigvecs * sqrt_diag) @ eigvecs.conj().T
        self.rho_invsqrt = (eigvecs * invsqrt_diag) @ eigvecs.conj().T

        # Delta_L on K_GNS: vec(rho X rho^{-1}) = ((rho^{-1})^T (x) rho) vec(X)
        rho_inv = np.linalg.inv(self.rho)
        self.Delta_matrix = np.kron(rho_inv.T, self.rho)

        # K_L_TT = -log Delta_L
        delta_eigvals, delta_eigvecs = np.linalg.eigh(self.Delta_matrix)
        delta_eigvals = np.maximum(delta_eigvals, 1e-30)
        log_diag = np.log(delta_eigvals)
        K_TT_diag = -log_diag
        self.K_TT = (delta_eigvecs * K_TT_diag) @ delta_eigvecs.conj().T

    def apply_J_L_TT(self, vec_X: np.ndarray) -> np.ndarray:
        """Apply the Lorentzian-side Tomita antilinear J_L_TT to vec(X).

        Action: J_L_TT(X) = rho^{1/2} X^* rho^{-1/2}.

        The Tomita J on the Krein-positive GNS space has the same
        formal expression as the Riemannian-side Tomita J (Paper 42
        Prop 6.2, J_TT^2 = +I), because the Krein-positive completion
        coincides with the standard Hilbert-space GNS on the wedge-
        trace state rho_W^L.

        Returns vec(rho^{1/2} X^* rho^{-1/2}) in column-stacking form.
        """
        if self.rho_sqrt is None:
            self.construct()
        X = vec_X.reshape(self.dim_W_L, self.dim_W_L).T
        X_dag = X.conj().T
        Y = self.rho_sqrt @ X_dag @ self.rho_invsqrt
        return Y.T.reshape(-1)

    def J_L_TT_squared_residual(self) -> float:
        """Verify J_L_TT^2 = +I on the Krein-GNS space by sampling.

        Returns average ||J_L_TT^2(v) - v||_2 over the standard
        orthonormal basis of K_GNS.  Should be ~ machine precision.
        """
        if self.rho is None:
            self.construct()
        total = 0.0
        for i in range(self.dim_GNS):
            e = np.zeros(self.dim_GNS, dtype=np.complex128)
            e[i] = 1.0
            J_e = self.apply_J_L_TT(e)
            JJ_e = self.apply_J_L_TT(J_e)
            total += float(np.linalg.norm(JJ_e - e))
        return total / self.dim_GNS

    def modular_flow_on_algebra(self, t: float, a: np.ndarray) -> np.ndarray:
        """Lorentzian Tomita modular automorphism on the algebra.

        sigma_t^{L,TT}(a)  =  rho^{it} a rho^{-it}

        Acts on operators a in B(K_W).  At t = 2*pi this closes
        bit-exactly because beta H_local = K_L_alpha^W has integer
        spectrum.

        Parameters
        ----------
        t : float
        a : ndarray of shape (dim_W_L, dim_W_L)

        Returns
        -------
        sigma_t^{L,TT}(a)
        """
        if self.rho is None:
            self.construct()
        eigvals, eigvecs = np.linalg.eigh(self.rho)
        eigvals = np.maximum(eigvals, 1e-30)
        log_eigs = np.log(eigvals)
        phases = np.exp(1j * t * log_eigs)
        rho_it = (eigvecs * phases) @ eigvecs.conj().T
        rho_minus_it = (eigvecs * np.conj(phases)) @ eigvecs.conj().T
        return rho_it @ a @ rho_minus_it

    def verify_modular_periodicity_tomita_lorentzian(
        self,
        a: np.ndarray,
        period: float = 2.0 * np.pi,
        tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Verify sigma_{period}^{L,TT}(a) = a on the algebra (L2E-FALS-2).

        For period = 2*pi and the BW choice H_local = K_L_alpha^W / beta,
        the closure is bit-exact because beta H_local = K_L_alpha^W has
        integer spectrum.

        Returns (ok, residual = ||sigma_{period}^{L,TT}(a) - a||_F).
        """
        sigma_a = self.modular_flow_on_algebra(period, a)
        residual = float(np.linalg.norm(sigma_a - a))
        return residual < tol, residual


# ---------------------------------------------------------------------------
# Full Lorentzian Modular Hamiltonian driver class
# ---------------------------------------------------------------------------


@dataclass
class LorentzianModularHamiltonian:
    """Modular Hamiltonian on the Lorentzian wedge of K = H_GV (x) C^{N_t}.

    Two construction modes (mirroring Paper 42 §5-§6):

      BW-alpha (geometric):  K_L_alpha = K_alpha_spatial (x) I_{N_t}
          with integer spectrum two_m_j on the wedge.  Period closure
          sigma_{2pi}^{L,alpha}(O) = O bit-exact.

      BW-gamma (Tomita-Takesaki):  K_L_TT = -log Delta_L on the
          Krein-GNS Hilbert-Schmidt space.  Period closure
          sigma_{2pi}^{L,TT}(O) = O bit-exact.

    Plus flow conjugacy at general t:
          sigma_t^{L,TT}(O) = sigma_{-t}^{L,alpha}(O) (Paper 42 §7.1).

    Plus Riemannian-limit recovery (LOAD-BEARING falsifier):
          K_L_alpha | _{N_t = 1} = K_alpha (Paper 42 §5)
          K_L_TT    | _{N_t = 1} = K_TT    (Paper 42 §6)

    Plus H_local verdict (Paper 42 §7.2 / O3 at (3, 1)):
          compare H_local := K_L_alpha^W / beta to the wedge-restricted
          truthful D_L_W, report bit-exact whether they agree.

    Attributes
    ----------
    krein : KreinSpace
    wedge : LorentzianWedge (built in __post_init__)
    kappa_g : float (default 1.0, BW canonical)
    K_L_alpha : ndarray of shape (dim_K, dim_K) (BW-alpha full Krein-level)
    """

    krein: KreinSpace
    kappa_g: float = 1.0

    # Built in __post_init__
    wedge: LorentzianWedge = field(init=False)
    K_L_alpha: np.ndarray = field(init=False)
    K_L_alpha_W: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.wedge = LorentzianWedge(krein=self.krein)
        self.K_L_alpha = build_K_L_alpha_full(self.krein)
        self.K_L_alpha_W = restrict_K_L_alpha_to_wedge(
            self.krein, self.wedge
        )

    @property
    def beta(self) -> float:
        """Inverse temperature beta = 2*pi / kappa_g."""
        return 2.0 * np.pi / self.kappa_g

    @property
    def dim_K(self) -> int:
        return self.krein.dim

    @property
    def dim_W_L(self) -> int:
        return self.wedge.dim_W_L

    # ------------------------------------------------------------------
    # BW-alpha geometric flow on the algebra
    # ------------------------------------------------------------------

    def modular_flow_alpha_on_wedge(
        self, t: float, a_W: np.ndarray,
    ) -> np.ndarray:
        """sigma_t^{L,alpha}(a_W) = e^{it K_L_alpha^W} a_W e^{-it K_L_alpha^W}.

        At t = 2*pi the closure is exact because K_L_alpha^W has
        integer spectrum.

        Parameters
        ----------
        t : float
        a_W : ndarray of shape (dim_W_L, dim_W_L)
        """
        K_diag = np.diag(self.K_L_alpha_W)
        # K_L_alpha_W is diagonal by construction
        phases = np.exp(1j * t * K_diag)
        return (phases[:, None] * a_W) * np.conj(phases[None, :])

    def verify_modular_periodicity_alpha_lorentzian(
        self, a_W: np.ndarray, tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Verify sigma_{2pi}^{L,alpha}(a_W) = a_W (L2E-FALS-1, LOAD-BEARING).

        For BW-alpha geometric K with integer spectrum two_m_j on the
        wedge, this should hold bit-exactly because e^{i 2pi n} = 1 for
        integer n.

        Returns (ok, residual).
        """
        sigma_a = self.modular_flow_alpha_on_wedge(2.0 * np.pi, a_W)
        residual = float(np.linalg.norm(sigma_a - a_W))
        return residual < tol, residual

    # ------------------------------------------------------------------
    # BW-gamma (Tomita) modular structure on Krein-GNS
    # ------------------------------------------------------------------

    def tomita_structure(
        self, H_local: Optional[np.ndarray] = None,
    ) -> LorentzianTomitaStructure:
        """Build the Tomita-Takesaki modular structure on K_GNS.

        For the BW choice H_local := K_L_alpha^W / beta, the modular
        flow on the algebra is the Heisenberg flow under K_L_alpha^W
        with integer spectrum, giving sigma_{2pi}^{L,TT}(a) = a
        bit-exactly.

        Parameters
        ----------
        H_local : ndarray of shape (dim_W_L, dim_W_L), optional
            Local Hamiltonian.  If None, uses the BW choice
            K_L_alpha^W / beta.
        """
        if H_local is None:
            H_local = self.K_L_alpha_W / self.beta
        ts = LorentzianTomitaStructure(
            H_local=H_local,
            beta=self.beta,
            K_L_alpha_W=self.K_L_alpha_W,
            dim_W_L=self.dim_W_L,
        )
        ts.construct()
        return ts

    def k_tomita(
        self, H_local: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """Return K_L_TT = -log Delta_L on K_GNS (load-bearing BW-gamma)."""
        return self.tomita_structure(H_local=H_local).K_TT

    # ------------------------------------------------------------------
    # Flow conjugacy
    # ------------------------------------------------------------------

    def compare_alpha_vs_tomita_lorentzian(
        self, a_W: np.ndarray, t: float = 1.0,
    ) -> dict:
        """Compare BW-alpha and BW-gamma flows at general t.

        For H_local = K_L_alpha^W / beta:
          sigma_t^{L,TT}(a)    =  e^{-it K_L_alpha^W} a e^{+it K_L_alpha^W}
          sigma_t^{L,alpha}(a) =  e^{+it K_L_alpha^W} a e^{-it K_L_alpha^W}

        These differ by the sign of t (Tomita convention: rho^{it} on
        the left).  At t = 2*pi the sign is irrelevant
        (e^{i 2pi n} = e^{-i 2pi n} = 1 for integer n).

        Returns dict with both flows, their difference, and the
        conjugation relation sigma_t^{L,TT} = sigma_{-t}^{L,alpha}.
        """
        ts = self.tomita_structure()
        sigma_TT = ts.modular_flow_on_algebra(t, a_W)
        sigma_alpha = self.modular_flow_alpha_on_wedge(t, a_W)
        sigma_alpha_neg = self.modular_flow_alpha_on_wedge(-t, a_W)
        return {
            "sigma_TT": sigma_TT,
            "sigma_alpha": sigma_alpha,
            "sigma_alpha_neg_t": sigma_alpha_neg,
            "TT_vs_alpha_diff": float(np.linalg.norm(sigma_TT - sigma_alpha)),
            "TT_vs_alpha_neg_t_diff": float(
                np.linalg.norm(sigma_TT - sigma_alpha_neg)
            ),
        }

    def verify_flow_conjugacy_lorentzian(
        self, a_W: np.ndarray, t: float = 1.0, tol: float = 1e-10,
    ) -> Tuple[bool, float]:
        """Verify sigma_t^{L,TT}(O) = sigma_{-t}^{L,alpha}(O) (L2E-FALS-3).

        Returns (ok, residual = ||sigma_t^{L,TT} - sigma_{-t}^{L,alpha}||_F).
        """
        cmp = self.compare_alpha_vs_tomita_lorentzian(a_W, t=t)
        residual = cmp["TT_vs_alpha_neg_t_diff"]
        return residual < tol, residual

    # ------------------------------------------------------------------
    # Riemannian-limit recovery (LOAD-BEARING falsifier)
    # ------------------------------------------------------------------

    def riemannian_limit_recovery(
        self, tol: float = 1e-14,
    ) -> Tuple[bool, dict]:
        """At N_t = 1, K_L_alpha and K_L_TT reduce bit-identically to K_alpha
        and K_TT of geovac.modular_hamiltonian.

        If self.N_t > 1, this method builds a fresh
        LorentzianModularHamiltonian at N_t = 1 (same n_max) and verifies
        that its K_L_alpha and K_L_TT match the Riemannian ModularHamiltonian
        constructions of Paper 42.  This is a LOAD-BEARING falsifier:
        failure indicates structural mismatch between the Lorentzian and
        Riemannian wedge constructions.

        Returns
        -------
        (ok, details) where details has 'k_alpha_residual',
        'k_tomita_residual', 'n_max', 'N_t_test'.
        """
        n_max = self.krein.n_max
        # Build the N_t = 1 reduction
        krein_rie = KreinSpace(n_max=n_max, N_t=1, T_max=1.0)
        lmh_rie = LorentzianModularHamiltonian(krein=krein_rie, kappa_g=1.0)

        # Riemannian reference from geovac.modular_hamiltonian
        from geovac.modular_hamiltonian import for_bisognano_wichmann
        mh_rie_ref = for_bisognano_wichmann(n_max=n_max, axis="hopf")
        K_alpha_W_rie = mh_rie_ref.restrict_K_alpha_to_wedge()

        # Compare K_L_alpha_W (Lorentzian at N_t = 1) vs K_alpha^W (Riemannian)
        k_alpha_residual = float(
            np.linalg.norm(lmh_rie.K_L_alpha_W - K_alpha_W_rie)
        )

        # Compare K_L_TT (Lorentzian at N_t = 1) vs K_TT (Riemannian)
        ts_lor = lmh_rie.tomita_structure()
        ts_rie = mh_rie_ref.tomita_structure()
        k_tomita_residual = float(np.linalg.norm(ts_lor.K_TT - ts_rie.K_TT))

        details = {
            "n_max": n_max,
            "N_t_test": 1,
            "k_alpha_W_residual": k_alpha_residual,
            "k_tomita_residual": k_tomita_residual,
            "dim_W_lorentzian_N_t_1": lmh_rie.dim_W_L,
            "dim_W_riemannian": int(np.round(np.real(
                np.trace(mh_rie_ref.wedge.projection_matrix(mh_rie_ref.basis))
            ))) // 1,  # ad-hoc
        }
        ok = (k_alpha_residual < tol) and (k_tomita_residual < tol)
        return ok, details

    # ------------------------------------------------------------------
    # H_local verdict (Paper 42 §7.2 / O3 at (3, 1))
    # ------------------------------------------------------------------

    def h_local_verdict_at_3_1(self, tol: float = 1e-12) -> dict:
        """Compare H_local := K_L_alpha^W / beta to wedge-restricted D_L_W.

        This is the Sprint L2-E headline structural test: the Lorentzian-side
        analog of Paper 42 §7.2 open question O3.

        Riemannian-side finding (Paper 42 §7.2):
          H_local := K_alpha^W / beta is NOT the wedge-restricted truthful
          Camporesi-Higuchi Dirac D_W. The framework's intrinsic Dirac is
          NOT the right local Hamiltonian for the BW vacuum at beta = 2*pi.

        Lorentzian-side question:
          Does the same finding hold at signature (3, 1)?  Specifically:
          is K_L_alpha^W / beta equal to D_L_W (the L2-C Lorentzian Dirac
          restricted to the wedge), or does the same H_local != D_W
          finding extend?

        Verdicts:
          (i)  K_L_alpha^W / beta = D_L_W at (3, 1):
               H_local choice IS the Lorentzian Dirac at signature (3, 1)
               but NOT at (3, 0).  Refines Paper 42 §7.2 cleanly.
          (ii) K_L_alpha^W / beta != D_L_W at (3, 1), same as Riemannian:
               O3 holds universally; the spectral-action-vs-modular-
               Hamiltonian generator distinction is signature-independent.
          (iii) Some intermediate: document structural relationship.

        Comparison metrics (multiple, because D_L_W has i-factor):
          - residual_full     = ||H_local - D_L_W||_F
                               (raw Frobenius; affected by D_L = i*Hermitian)
          - residual_vs_iDGV  = ||H_local - i*D_GV_W||_F
                               (compares to i * D_GV_W which is the natural
                                "Lorentzian Dirac restricted to wedge"
                                content at signature (3, 1))
          - residual_RIE_check = ||H_local_RIE - D_GV_W||_F
                               (Riemannian-side baseline at the same n_max)

        At N_t = 1 the Lorentzian Dirac reduces to D_L = i * D_GV exactly
        (because the temporal derivative vanishes on a singleton grid), so
        D_L_W = i * D_GV_W.  At N_t > 1 the temporal derivative contributes
        and D_L_W gains additional matrix structure off-diagonal in the
        temporal slot.

        Returns dict with verdict ('i', 'ii', or 'iii'), the comparison
        data, and an explicit structural reading.
        """
        # H_local := K_L_alpha^W / beta on the wedge
        H_local_wedge = self.K_L_alpha_W / self.beta

        # D_L on the full Krein space, then restrict to the wedge block
        D_L_full = lorentzian_dirac_matrix(self.krein)
        D_L_W = restrict_operator_to_wedge_block(
            D_L_full, self.krein, self.wedge,
        )

        residual_full = float(np.linalg.norm(H_local_wedge - D_L_W))

        # Normalize by sqrt(dim) to enable cross-N_t comparison
        residual_full_normalized = residual_full / np.sqrt(self.dim_W_L)

        # Riemannian-side baseline: H_local_RIE - D_GV_W (Paper 42 §7.2)
        from geovac.modular_hamiltonian import for_bisognano_wichmann
        mh_rie = for_bisognano_wichmann(n_max=self.krein.n_max, axis="hopf")
        K_alpha_W_rie = mh_rie.restrict_K_alpha_to_wedge()
        beta_rie = mh_rie.beta
        H_local_rie = K_alpha_W_rie / beta_rie
        D_GV_W_rie = mh_rie.restrict_to_wedge_block(mh_rie.D)
        residual_RIE_check = float(np.linalg.norm(H_local_rie - D_GV_W_rie))
        dim_W_rie = D_GV_W_rie.shape[0]
        residual_RIE_normalized = residual_RIE_check / np.sqrt(dim_W_rie)

        # Compare to i * D_GV_W on the wedge (Lorentzian content at N_t = 1)
        # at general N_t this needs to be tensored with the temporal slot.
        if H_local_wedge.shape == D_GV_W_rie.shape:
            residual_vs_iDGV = float(
                np.linalg.norm(H_local_wedge - 1j * D_GV_W_rie)
            )
        else:
            # N_t > 1: compare H_local to i * D_GV_W (x) (positive-t identity).
            # The Lorentzian wedge has stride N_t_positive on temporal slot.
            N_t_pos = self.wedge.N_t_positive
            iDGV_extended = 1j * np.kron(
                D_GV_W_rie, np.eye(N_t_pos, dtype=np.complex128)
            )
            residual_vs_iDGV = float(
                np.linalg.norm(H_local_wedge - iDGV_extended)
            )

        # Spectrum sample of H_local
        H_local_eigs = np.sort(np.real(np.linalg.eigvalsh(H_local_wedge)))

        # Verdict logic
        # Pattern of findings (verified across panel n_max in {1,2,3} x N_t in {1,11,21}):
        #   - At N_t = 1: residual_full_normalized == residual_RIE_normalized
        #     bit-exact. The Lorentzian wedge restriction at the Riemannian
        #     limit IS structurally the Riemannian wedge restriction.
        #     Paper 42 §7.2 O3 holds at signature (3, 1) as a signature-
        #     independent baseline.
        #   - At N_t > 1: residual_full grows with N_t, because D_L has the
        #     temporal-derivative piece i * gamma^0 (x) d/dt which H_local
        #     does not capture (H_local is purely spatial via K_L_alpha^W).
        if residual_full < tol:
            verdict = "i"
            verdict_text = (
                "(i) K_L_alpha^W / beta = D_L_W at (3,1) bit-exact: "
                "H_local choice IS the Lorentzian Dirac at signature (3,1) "
                "but NOT at (3,0). Refines Paper 42 §7.2 to a "
                "signature-DEPENDENT structural finding."
            )
        else:
            # Significant disagreement; classify by Riemannian-baseline matching
            relative_mismatch = (
                abs(residual_full_normalized - residual_RIE_normalized)
                / max(residual_full_normalized, 1e-30)
            )
            if relative_mismatch < 1e-10:
                # Bit-exact match to Riemannian-side residual
                verdict = "ii"
                verdict_text = (
                    "(ii) K_L_alpha^W / beta != D_L_W at (3,1), with the "
                    "EXACTLY SAME structural residual as the Riemannian "
                    "side (signature-INDEPENDENT). Paper 42 §7.2 O3 holds "
                    "UNIVERSALLY at the Riemannian-limit reduction (N_t = 1). "
                    "The framework's intrinsic Dirac is NOT the right local "
                    "Hamiltonian for the BW vacuum at beta = 2*pi at EITHER "
                    "signature."
                )
            else:
                # Lorentzian residual differs from Riemannian baseline
                # (typical at N_t > 1 due to temporal-derivative contribution
                # in D_L = i (gamma^0 (x) d/dt + D_GV (x) I))
                verdict = "iii"
                verdict_text = (
                    f"(iii) Intermediate: residual_full = {residual_full:.4e} "
                    f"(normalized {residual_full_normalized:.4e}) vs "
                    f"Riemannian baseline {residual_RIE_check:.4e} "
                    f"(normalized {residual_RIE_normalized:.4e}). The "
                    f"Lorentzian-side residual differs from the Riemannian "
                    f"baseline structurally: D_L has the temporal-derivative "
                    f"piece i * gamma^0 (x) d/dt which H_local := "
                    f"K_L_alpha^W / beta does not capture. Paper 42 §7.2 "
                    f"O3 PRESERVED at the Riemannian limit (N_t = 1); "
                    f"REFINED at N_t > 1 by the temporal-derivative "
                    f"contribution."
                )

        return {
            "n_max": self.krein.n_max,
            "N_t": self.krein.N_t,
            "beta": self.beta,
            "dim_W_L": self.dim_W_L,
            "residual_full": residual_full,
            "residual_full_normalized": residual_full_normalized,
            "residual_RIE_baseline": residual_RIE_check,
            "residual_RIE_normalized": residual_RIE_normalized,
            "residual_vs_iDGV": residual_vs_iDGV,
            "H_local_spectrum_sample": [float(x) for x in H_local_eigs[:5]],
            "verdict": verdict,
            "verdict_text": verdict_text,
            "lorentzian_eq_riemannian_baseline": bool(
                abs(residual_full_normalized - residual_RIE_normalized) < 1e-10
            ),
        }

    # ------------------------------------------------------------------
    # Witness verification battery (Paper 42 §5-§8 lifted)
    # ------------------------------------------------------------------

    def verify_witness_lorentzian(
        self,
        generators: Optional[Sequence[np.ndarray]] = None,
        tol: float = 1e-10,
    ) -> dict:
        """Run the full Sprint L2-E witness verification battery.

        Tests:
          - Wedge P_W_L^2 = P_W_L (idempotent)
          - Wedge P_W_L^* = P_W_L (Hermitian)
          - L2E-FALS-1: sigma_{2pi}^{L,alpha}(O) = O bit-exact (LOAD-BEARING)
          - L2E-FALS-2: sigma_{2pi}^{L,TT}(O) = O bit-exact
          - L2E-FALS-3: flow conjugacy sigma_t^{L,TT} = sigma_{-t}^{L,alpha}

        Returns dict with per-test booleans, residuals, and verdict.
        """
        if generators is None:
            generators = self._default_generators()

        results = {
            "kappa_g": self.kappa_g,
            "beta": self.beta,
            "n_max": self.krein.n_max,
            "N_t": self.krein.N_t,
            "dim_K": self.dim_K,
            "dim_W_L": self.dim_W_L,
        }

        # Wedge tests
        ok_idem, res_idem = self.wedge.verify_idempotent(tol=tol)
        results["wedge_idempotent"] = {
            "ok": bool(ok_idem), "residual": res_idem,
        }
        ok_herm, res_herm = self.wedge.verify_hermitian(tol=tol)
        results["wedge_hermitian"] = {
            "ok": bool(ok_herm), "residual": res_herm,
        }

        # Period closures per generator
        ts = self.tomita_structure()
        per_gen = []
        for k, O in enumerate(generators):
            # Restrict to wedge block
            O_W = restrict_operator_to_wedge_block(
                O, self.krein, self.wedge,
            )
            # L2E-FALS-1: BW-alpha period closure
            ok_alpha, res_alpha = (
                self.verify_modular_periodicity_alpha_lorentzian(O_W, tol=tol)
            )
            # L2E-FALS-2: BW-gamma Tomita period closure
            ok_tt, res_tt = ts.verify_modular_periodicity_tomita_lorentzian(
                O_W, tol=tol,
            )
            # L2E-FALS-3: flow conjugacy at t = 1
            ok_conj, res_conj = self.verify_flow_conjugacy_lorentzian(
                O_W, t=1.0, tol=tol,
            )
            per_gen.append({
                "gen_idx": k,
                "L2E_FALS_1_alpha": {"ok": bool(ok_alpha), "residual": res_alpha},
                "L2E_FALS_2_TT": {"ok": bool(ok_tt), "residual": res_tt},
                "L2E_FALS_3_conjugacy_t1": {"ok": bool(ok_conj), "residual": res_conj},
            })
        results["per_generator"] = per_gen

        max_alpha = max(
            (g["L2E_FALS_1_alpha"]["residual"] for g in per_gen),
            default=0.0,
        )
        max_tt = max(
            (g["L2E_FALS_2_TT"]["residual"] for g in per_gen),
            default=0.0,
        )
        max_conj = max(
            (g["L2E_FALS_3_conjugacy_t1"]["residual"] for g in per_gen),
            default=0.0,
        )

        results["max_alpha_residual"] = max_alpha
        results["max_TT_residual"] = max_tt
        results["max_conjugacy_residual"] = max_conj

        # Summary verdict
        all_pass = (max_alpha < tol) and (max_tt < tol) and (max_conj < tol)
        results["verdict"] = (
            "STRONG_IDENTIFICATION_LORENTZIAN" if all_pass
            else ("SOFT_IDENTIFICATION" if max(max_alpha, max_tt, max_conj) < 1.0
                  else "FAIL")
        )

        return results

    def _default_generators(
        self, n_test: int = 5,
    ) -> List[np.ndarray]:
        """Default test generators: 5 non-identity multipliers extended to K.

        Builds spatial multipliers from FullDiracTruncatedOperatorSystem
        and extends to the full Krein space by tensoring with I_{N_t}.

        At n_max = 1 the truncated operator system has only the identity;
        we fall back to a structurally-faithful diagonal weight to give
        a non-trivial test object. The period closure is trivial for any
        diagonal operator commuting with K_L_alpha^W, but the test
        machinery still exercises the wedge restriction + flow.
        """
        op_sys = FullDiracTruncatedOperatorSystem(n_max=self.krein.n_max)
        I_t = np.eye(self.krein.N_t, dtype=np.complex128)
        gens = []
        for label, M_spatial in op_sys.basis_matrices:
            if label == (1, 0, 0):
                continue
            M_full = np.kron(M_spatial, I_t)
            gens.append(M_full)
            if len(gens) >= n_test:
                break
        if not gens:
            # Fallback at n_max = 1
            basis = full_dirac_basis(self.krein.n_max)
            diag_weights = np.array(
                [float(b.n_fock) + 0.7 * float(b.l) + 0.3 * float(b.two_m_j)
                 for b in basis],
                dtype=np.complex128,
            )
            M_spatial = np.diag(diag_weights)
            gens.append(np.kron(M_spatial, I_t))
        return gens


# ---------------------------------------------------------------------------
# Witness factories (parameterized BW per Paper 42 §3.3)
# ---------------------------------------------------------------------------


def for_bisognano_wichmann_lorentzian(
    n_max: int, N_t: int = 11, T_max: float = 1.0,
) -> LorentzianModularHamiltonian:
    """BW canonical: hemisphere wedge of S^3 x R, beta = 2*pi, kappa_g = 1.

    Returns LorentzianModularHamiltonian at kappa_g = 1, beta = 2*pi.
    """
    krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=T_max)
    return LorentzianModularHamiltonian(krein=krein, kappa_g=1.0)


def for_hartle_hawking_lorentzian(
    n_max: int, N_t: int = 11, T_max: float = 1.0, M: float = 1.0,
) -> LorentzianModularHamiltonian:
    """Hartle-Hawking: kappa_g = 1/(4M), beta = 8*pi*M."""
    krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=T_max)
    return LorentzianModularHamiltonian(krein=krein, kappa_g=1.0 / (4.0 * M))


def for_sewell_lorentzian(
    n_max: int, N_t: int = 11, T_max: float = 1.0, M: float = 1.0,
) -> LorentzianModularHamiltonian:
    """Sewell 1982: same parameters as Hartle-Hawking."""
    return for_hartle_hawking_lorentzian(n_max, N_t=N_t, T_max=T_max, M=M)


def for_unruh_lorentzian(
    n_max: int, N_t: int = 11, T_max: float = 1.0, a: float = 1.0,
) -> LorentzianModularHamiltonian:
    """Unruh: kappa_g = a, beta = 2*pi/a."""
    krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=T_max)
    return LorentzianModularHamiltonian(krein=krein, kappa_g=float(a))


# ---------------------------------------------------------------------------
# Six-witness collapse at the Krein level (Paper 42 §8 lifted)
# ---------------------------------------------------------------------------


def verify_cross_witness_collapse_lorentzian(
    n_max: int, N_t: int = 11, T_max: float = 1.0, tol: float = 1e-10,
) -> dict:
    """Verify HH, Sew, Unruh collapse to BW at the Lorentzian operator-system level.

    The six physical witnesses (BW canonical at kappa_g = 1; HH at M = 1,
    kappa_g = 1/4; HH at M = 2, kappa_g = 1/8; Sew at M = 1, kappa_g = 1/4;
    Unruh at a = 1, kappa_g = 1; Unruh at a = 2, kappa_g = 2) all give
    bit-identical period closures because K_L_alpha^W has integer
    spectrum independent of kappa_g.

    Cross-witness consistency residual is the difference between the
    period residual of each witness and the BW canonical period
    residual at the same (n_max, N_t).

    Returns dict with per-witness residuals and overall collapse verdict.
    """
    # Pick a sample test operator from the spatial multiplier basis.
    # At n_max = 1 the truncated operator system has only the identity;
    # in that case, build a structurally-faithful test multiplier as the
    # K_alpha_W generator itself (which preserves the wedge by construction).
    op_sys = FullDiracTruncatedOperatorSystem(n_max=n_max)
    O_spatial = None
    for label, M in op_sys.basis_matrices:
        if label != (1, 0, 0):
            O_spatial = M
            break
    if O_spatial is None:
        # n_max = 1 fallback: use the identity multiplier (only available choice)
        # plus a diagonal "test operator" that doesn't break the wedge
        # structure. The period closure is trivial for any operator
        # commuting with K_L_alpha^W; we use a diagonal weight to give a
        # non-trivial sanity check.
        from geovac.full_dirac_operator_system import full_dirac_basis
        basis = full_dirac_basis(n_max)
        diag_weights = np.array(
            [float(b.n_fock) + 0.7 * float(b.l) + 0.3 * float(b.two_m_j)
             for b in basis],
            dtype=np.complex128,
        )
        O_spatial = np.diag(diag_weights)
    I_t = np.eye(N_t, dtype=np.complex128)
    O_full = np.kron(O_spatial, I_t)

    bw = for_bisognano_wichmann_lorentzian(n_max, N_t=N_t, T_max=T_max)
    hh1 = for_hartle_hawking_lorentzian(n_max, N_t=N_t, T_max=T_max, M=1.0)
    hh2 = for_hartle_hawking_lorentzian(n_max, N_t=N_t, T_max=T_max, M=2.0)
    sew = for_sewell_lorentzian(n_max, N_t=N_t, T_max=T_max, M=1.0)
    un1 = for_unruh_lorentzian(n_max, N_t=N_t, T_max=T_max, a=1.0)
    un2 = for_unruh_lorentzian(n_max, N_t=N_t, T_max=T_max, a=2.0)

    results = {"n_max": n_max, "N_t": N_t, "T_max": T_max}
    bw_residuals = {}
    for name, mh in [
        ("BW", bw), ("HH_M1", hh1), ("HH_M2", hh2), ("Sew_M1", sew),
        ("Unruh_a1", un1), ("Unruh_a2", un2),
    ]:
        O_W = restrict_operator_to_wedge_block(O_full, mh.krein, mh.wedge)
        # BW-alpha closure
        ok_alpha, res_alpha = (
            mh.verify_modular_periodicity_alpha_lorentzian(O_W, tol=tol)
        )
        # BW-gamma closure
        ts = mh.tomita_structure()
        ok_tt, res_tt = ts.verify_modular_periodicity_tomita_lorentzian(
            O_W, tol=tol,
        )
        results[name] = {
            "kappa_g": mh.kappa_g, "beta": mh.beta,
            "alpha_ok": bool(ok_alpha), "alpha_residual": res_alpha,
            "TT_ok": bool(ok_tt), "TT_residual": res_tt,
        }
        bw_residuals[name] = (res_alpha, res_tt)

    # Cross-consistency: residuals match BW (independent of kappa_g)
    bw_alpha = results["BW"]["alpha_residual"]
    bw_tt = results["BW"]["TT_residual"]
    consistency = {}
    for name in ["HH_M1", "HH_M2", "Sew_M1", "Unruh_a1", "Unruh_a2"]:
        consistency[name] = {
            "alpha_diff": float(np.abs(results[name]["alpha_residual"] - bw_alpha)),
            "TT_diff": float(np.abs(results[name]["TT_residual"] - bw_tt)),
        }
    results["cross_consistency"] = consistency
    max_alpha_cons = max(c["alpha_diff"] for c in consistency.values())
    max_tt_cons = max(c["TT_diff"] for c in consistency.values())
    results["max_alpha_cross_consistency"] = max_alpha_cons
    results["max_TT_cross_consistency"] = max_tt_cons

    # Cross-witness collapse verdict
    collapse_ok = (max_alpha_cons < tol) and (max_tt_cons < tol)
    results["six_witness_collapse_ok"] = bool(collapse_ok)

    return results
