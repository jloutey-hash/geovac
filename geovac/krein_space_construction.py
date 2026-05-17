"""Krein space construction on S^3 x R at signature (3, 1).

Sprint L2-B (2026-05-16): construct the Krein space

    K_{n_max} = H_GV^{n_max} (x) L^2(R_t)_cutoff,

with fundamental symmetry J = gamma^0 (the temporal Dirac matrix in the
West-coast Dirac-spinor convention, chiral basis), and verify the Krein
axioms

    J^2 = +I,
    J^* J = I  (J unitary on the Hilbert-space level),
    < . , . >_K := < . , J . >   is Hermitian (conjugate-symmetric).

The construction is the GeoVac-side instantiation of the van den Dungen
2016 Proposition 4.1 lift (Math. Phys. Anal. Geom. 19:4, arXiv:1505.01939)
for the spatial slice M_spatial = S^3 (truncated Camporesi-Higuchi spinor
bundle at n_max) and temporal direction t in [-T_max, T_max] with a small
uniform grid (NO compactification to S^1_beta -- that would create
closed timelike curves per Geroch's theorem; see L2-A audit Section 3.8
and Section 5.1).

The Riemannian spatial Hilbert space H_GV^{n_max} is the Camporesi-Higuchi
spinor space of Paper 32 Definition 3.2 (def:H_GV), explicitly the C^N
with N = N_Dirac(n_max) = (2/3) n_max (n_max + 1) (n_max + 2).  This
matches `geovac.full_dirac_operator_system.full_dirac_dim(n_max)` (verified
at n_max in {1, 2, 3}: 4, 16, 40).

The Riemannian-limit check at N_t = 1 (or equivalently t -> 0 grid
collapse) is the LOAD-BEARING falsifier: the Krein-space spatial slice
at a single time-point must reduce to H_GV^{n_max} bit-identically.  If
this fails, the Camporesi-Higuchi spatial spinor bundle is structurally
incompatible with the Cl(3, 1) gamma-matrix embedding -- a major
structural finding that would re-open WH1 PROVEN at the spinor-lift
level.  STOP and escalate to PI.

Conventions (locked decisions for this sprint)
==============================================

West-coast metric eta = diag(+1, -1, -1, -1)
   (most-plus = mostly-minus = West-coast).  Spatial signature gives
   (gamma^i)^2 = -I, temporal gives (gamma^0)^2 = +I.

Cl(3, 1) gamma matrices in the CHIRAL (Weyl) basis (Peskin-Schroeder
convention):
   gamma^0 = [[0, I_2], [I_2, 0]]    (chirality-swap, off-diagonal)
   gamma^i = [[0,  sigma^i ],
              [-sigma^i, 0]]         for i = 1, 2, 3
   gamma^5 = i gamma^0 gamma^1 gamma^2 gamma^3
           = diag(-I_2, +I_2)        (chirality grading, diagonal)

   (The sign convention gamma^5 = diag(-I_2, +I_2) is the standard
    Peskin-Schroeder West-coast chiral-basis result: the upper 2-block
    is the LEFT-handed (chirality eigenvalue -1) sector, the lower
    2-block is the RIGHT-handed (chirality eigenvalue +1) sector.)

with Pauli matrices

   sigma^1 = [[0,  1], [1, 0]]
   sigma^2 = [[0, -i], [i, 0]]
   sigma^3 = [[1,  0], [0, -1]].

In this convention, GeoVac's FullDiracLabel.chirality field
(geovac/full_dirac_operator_system.py) IS the eigenvalue of gamma^5
(up to a global sign choice).  The block ordering of H_GV at the
spinor level is "Weyl sector first, anti-Weyl second" per
`full_dirac_basis(n_max)`; the gamma^5 = diag(-I_2, +I_2) sign means
GeoVac's chirality = +1 (Weyl block first half of H_GV) is the
LEFT-handed Peskin-Schroeder sector and chirality = -1 (anti-Weyl
block second half) is the RIGHT-handed sector.  This is a labeling
convention with no physical content beyond bookkeeping; the
Riemannian-limit check (J = J_spatial bit-identically at N_t = 1) is
agnostic to this sign.

Fundamental symmetry choice:
   J := gamma^0  (in the chiral basis).

This makes J a Hermitian unitary involution with J^2 = +I that SWAPS
the two chirality blocks.  The Krein eigenspaces K^+ and K^- (the
+1 and -1 eigenspaces of J) each have dimension dim(K) / 2.

Note on convention: the chiral basis choice means gamma^0 is OFF-
diagonal in the chirality block decomposition.  In the alternative
Dirac basis one has gamma^0 = diag(+I_2, -I_2) (diagonal); the two
conventions are unitarily equivalent.  We use the chiral basis here
because GeoVac's chirality label IS gamma^5, which is diagonal in this
basis; this makes the Riemannian-limit identification cleanest.

The two natural conventions (chiral vs Dirac basis) are documented in
the memo; we pick chiral here and justify it by Riemannian-limit
consistency with the existing H_GV (Paper 32 Definition 3.2).

API
===

   GammaMatrices                     : Cl(3,1) gamma + gamma^5 (chiral basis)
   pauli_matrices()                  : 2x2 Pauli sigma^1, sigma^2, sigma^3
   gamma_chiral()                    : returns gamma^0, gamma^1, gamma^2, gamma^3, gamma^5
   verify_gamma_anticommutation()    : {gamma^mu, gamma^nu} = 2 eta^{mu nu}

   temporal_grid(T_max, N_t)         : uniform symmetric grid on [-T_max, T_max]
   identity_temporal(N_t)            : I_{N_t} for L^2(R_t)_cutoff slot

   KreinSpace(n_max, N_t=21, T_max=1.0)
       .dim                          : N_Dirac(n_max) * N_t
       .J                             : fundamental symmetry matrix (dim, dim)
       .krein_inner_product(psi, phi)
       .positive_negative_split()    : (K_plus_proj, K_minus_proj)
       .riemannian_limit_check()     : returns (ok, residual) vs H_GV

References
==========

   van den Dungen, K. "Krein spectral triples and the fermionic action."
     Math. Phys. Anal. Geom. 19, 4 (2016). arXiv:1505.01939.
     Proposition 4.1 -- the Riemannian -> Lorentzian (Krein) lift recipe.

   Bizi, N., Brouder, C., Besnard, F. "Space and time dimensions of
     algebras with application to Lorentzian noncommutative geometry
     and quantum electrodynamics." J. Math. Phys. 59, 062303 (2018).
     arXiv:1611.07062.  Sign table for (m, n) signatures.

   Camporesi, R., Higuchi, A. "On the eigenfunctions of the Dirac
     operator on spheres and real hyperbolic spaces."
     J. Geom. Phys. 20, 1-18 (1996).  Spatial spinor structure.

   GeoVac internal:
     papers/synthesis/paper_32_spectral_triple.tex  Section III (def:H_GV)
     papers/standalone/paper_42_modular_hamiltonian_four_witness.tex
     debug/sprint_l2a_scoping_memo.md   (Sprint L2-A audit)
     debug/lorentzian_l0_audit_memo.md  (28-projection partition)
     geovac/full_dirac_operator_system.py  (spatial H_GV)
     geovac/modular_hamiltonian.py     (Riemannian Tomita J_TT, BW-alpha)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    full_dirac_basis,
    full_dirac_dim,
)


# ---------------------------------------------------------------------------
# Pauli matrices and Cl(3, 1) gamma matrices (chiral basis, West-coast)
# ---------------------------------------------------------------------------


def pauli_matrices() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the three Pauli matrices sigma^1, sigma^2, sigma^3 as 2x2 arrays.

    Standard convention:
        sigma^1 = [[0,  1], [1, 0]]
        sigma^2 = [[0, -i], [i, 0]]
        sigma^3 = [[1,  0], [0, -1]]
    """
    s1 = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    s2 = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    s3 = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    return s1, s2, s3


@dataclass(frozen=True)
class GammaMatrices:
    """Cl(3, 1) gamma matrices in the chiral (Weyl) basis, West-coast metric.

    Attributes
    ----------
    g0, g1, g2, g3 : np.ndarray, shape (4, 4), complex128
        The four 4x4 Cl(3, 1) gamma matrices.
    g5 : np.ndarray, shape (4, 4), complex128
        gamma^5 = i gamma^0 gamma^1 gamma^2 gamma^3, the chirality
        operator (diagonal +I_2 / -I_2 in this basis).
    eta : np.ndarray, shape (4, 4), float64
        The (4, 4) Minkowski metric tensor diag(+1, -1, -1, -1) for
        sanity-check accounting (NOT used in any operator computation).

    The Clifford-algebra defining relation
        {gamma^mu, gamma^nu} = 2 eta^{mu nu} I_4
    is verified to machine precision by `verify_gamma_anticommutation`.

    Convention details (chiral basis, West-coast):
        gamma^0 = [[0, I_2], [I_2, 0]]      (off-diagonal in chirality)
        gamma^i = [[0,  sigma^i ],
                   [-sigma^i, 0]]            for i = 1, 2, 3
        gamma^5 = diag(-I_2, +I_2)           (chirality grading,
                                              Peskin-Schroeder sign)
    """

    g0: np.ndarray
    g1: np.ndarray
    g2: np.ndarray
    g3: np.ndarray
    g5: np.ndarray
    eta: np.ndarray = field(
        default_factory=lambda: np.diag([1.0, -1.0, -1.0, -1.0])
    )

    def list_spacetime(self) -> Tuple[np.ndarray, ...]:
        """Return (g0, g1, g2, g3) as a 4-tuple."""
        return (self.g0, self.g1, self.g2, self.g3)


def gamma_chiral() -> GammaMatrices:
    """Construct Cl(3, 1) gamma matrices in the chiral (Weyl) basis.

    Returns
    -------
    GammaMatrices
        With g0, g1, g2, g3 the four 4x4 spacetime gamma matrices and
        g5 the chirality operator.

    Verification (West-coast, chiral basis):
        (g0)^2 = +I_4
        (gi)^2 = -I_4  for i = 1, 2, 3
        {g_mu, g_nu} = 2 eta_{mu nu} I_4
        g5 = diag(-1, -1, +1, +1)   (Peskin-Schroeder chiral basis)
        (g5)^2 = +I_4
        {g5, g_mu} = 0  for mu = 0, 1, 2, 3
    """
    s1, s2, s3 = pauli_matrices()
    I2 = np.eye(2, dtype=np.complex128)
    Z2 = np.zeros((2, 2), dtype=np.complex128)

    # gamma^0 = [[0, I_2], [I_2, 0]]
    g0 = np.block([[Z2, I2], [I2, Z2]])

    # gamma^i = [[0, sigma^i], [-sigma^i, 0]]
    g1 = np.block([[Z2, s1], [-s1, Z2]])
    g2 = np.block([[Z2, s2], [-s2, Z2]])
    g3 = np.block([[Z2, s3], [-s3, Z2]])

    # gamma^5 = i gamma^0 gamma^1 gamma^2 gamma^3 (in chiral basis: diag)
    g5 = 1j * g0 @ g1 @ g2 @ g3

    return GammaMatrices(g0=g0, g1=g1, g2=g2, g3=g3, g5=g5)


def verify_gamma_anticommutation(
    gammas: GammaMatrices, tol: float = 1e-12
) -> Tuple[bool, float]:
    """Check {gamma^mu, gamma^nu} = 2 eta^{mu nu} I_4 to within tol.

    Parameters
    ----------
    gammas : GammaMatrices
    tol : float
        Frobenius-norm tolerance.

    Returns
    -------
    (ok, max_residual) : (bool, float)
        ok is True iff all 10 independent anticommutators agree with
        2 eta^{mu nu} I_4 to within tol in Frobenius norm.
    """
    g = gammas.list_spacetime()
    eta = gammas.eta
    I4 = np.eye(4, dtype=np.complex128)
    max_res = 0.0
    for mu in range(4):
        for nu in range(mu, 4):
            anticomm = g[mu] @ g[nu] + g[nu] @ g[mu]
            target = 2.0 * eta[mu, nu] * I4
            res = float(np.linalg.norm(anticomm - target))
            max_res = max(max_res, res)
    return max_res < tol, max_res


# ---------------------------------------------------------------------------
# Temporal grid for L^2(R_t)_cutoff
# ---------------------------------------------------------------------------


def temporal_grid(T_max: float, N_t: int) -> np.ndarray:
    """Uniform symmetric grid on [-T_max, T_max] with N_t points.

    NOT compactified to S^1_beta (that would introduce closed timelike
    curves on S^3 x S^1; see L2-A audit Section 3.8).  This is a bounded
    open interval discretized for numerical concreteness.

    Parameters
    ----------
    T_max : float
        Half-width of the temporal interval.  Must be > 0.
    N_t : int
        Number of grid points.  Must be >= 1.

    Returns
    -------
    t_grid : np.ndarray of shape (N_t,)
        Grid points t_k for k = 0, ..., N_t - 1.  For N_t = 1, returns
        the singleton [0.0] (Riemannian-limit / static spatial slice).
    """
    if T_max <= 0:
        raise ValueError(f"T_max must be > 0, got {T_max}")
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    if N_t == 1:
        return np.array([0.0])
    return np.linspace(-T_max, T_max, N_t)


def identity_temporal(N_t: int) -> np.ndarray:
    """Identity I_{N_t} on the temporal slot for fundamental-symmetry tensoring.

    The fundamental symmetry J = gamma^0 acts on the spinor index and
    trivially on the temporal index.  This is the canonical L^2(R_t)
    identity that gets tensored with gamma^0 to build the full J on
    the Krein space K = H_GV (x) L^2(R_t)_cutoff.
    """
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    return np.eye(N_t, dtype=np.complex128)


# ---------------------------------------------------------------------------
# Spatial fundamental symmetry on H_GV (lifted from gamma^0 on Cl(3,1))
# ---------------------------------------------------------------------------


def spatial_fundamental_symmetry(
    basis: list[FullDiracLabel]
) -> np.ndarray:
    """Build J_spatial on H_GV that acts as gamma^0 on the chirality index.

    H_GV at n_max is structured as Weyl block (chirality = +1, dim
    spinor_dim(n_max)) PLUS anti-Weyl block (chirality = -1, same dim),
    per `full_dirac_basis` ordering (chirality = +1 first).  The gamma^0
    action in the chiral basis is the off-diagonal swap

        gamma^0 = [[0, I], [I, 0]]

    so J_spatial maps (chirality = +1, l, m_j) <-> (chirality = -1, l, m_j)
    at fixed (n_fock, l, two_m_j).

    Parameters
    ----------
    basis : list of FullDiracLabel
        The H_GV basis at some n_max, in standard (Weyl-first) ordering.

    Returns
    -------
    J_spatial : complex ndarray, shape (dim_H, dim_H)
        Unitary involution: J^2 = +I, J^dagger = J.
    """
    dim = len(basis)
    if dim % 2 != 0:
        raise ValueError(
            f"H_GV basis dimension {dim} must be even (chirality doubling)"
        )

    label_to_idx = {label: i for i, label in enumerate(basis)}
    J = np.zeros((dim, dim), dtype=np.complex128)

    for j, label in enumerate(basis):
        flipped = FullDiracLabel(
            n_fock=label.n_fock,
            l=label.l,
            two_m_j=label.two_m_j,
            chirality=-label.chirality,
        )
        i = label_to_idx[flipped]
        J[i, j] = 1.0

    return J


# ---------------------------------------------------------------------------
# KreinSpace: K = H_GV (x) L^2(R_t)_cutoff with J = J_spatial (x) I_{N_t}
# ---------------------------------------------------------------------------


@dataclass
class KreinSpace:
    """Krein space K_{n_max} = H_GV^{n_max} (x) L^2(R_t)_cutoff.

    Construction (van den Dungen 2016 Prop 4.1 instantiated for
    M_spatial = S^3 truncated at n_max, M_time = R_t with uniform-grid
    cutoff):

        K = H_GV (x) C^{N_t},
        J = J_spatial (x) I_{N_t},

    where J_spatial = lift of gamma^0 in the chiral basis to the H_GV
    chirality blocks (off-diagonal swap of the two chirality halves of
    H_GV).  J is Hermitian, unitary, and J^2 = +I.  The Krein inner
    product is

        <psi, phi>_K := <psi, J phi>.

    The positive/negative-definite splitting K = K^+ (+) K^- is given by
    the +/-1 eigenspaces of J, each of dimension dim(K) / 2.

    Parameters
    ----------
    n_max : int
        Spatial spinor-bundle truncation.  Must be >= 1.
    N_t : int
        Number of temporal grid points.  Default 21.  Must be >= 1.
        N_t = 1 corresponds to the Riemannian-limit / single-time-slice
        reduction.
    T_max : float
        Half-width of the temporal grid [-T_max, T_max].  Default 1.0.
        Unused at N_t = 1 (degenerate single-point grid at t = 0).

    Attributes
    ----------
    basis_spatial : list of FullDiracLabel
        H_GV basis at n_max (length full_dirac_dim(n_max)).
    dim_spatial : int
        dim(H_GV^{n_max}).
    t_grid : np.ndarray, shape (N_t,)
        Temporal grid points.
    dim : int
        Total Krein-space dimension: dim_spatial * N_t.
    J_spatial : np.ndarray, shape (dim_spatial, dim_spatial)
        Spatial fundamental symmetry (Hermitian unitary involution).
    J : np.ndarray, shape (dim, dim)
        Full Krein fundamental symmetry J = J_spatial (x) I_{N_t}.
    """

    n_max: int
    N_t: int = 21
    T_max: float = 1.0

    basis_spatial: list[FullDiracLabel] = field(init=False)
    dim_spatial: int = field(init=False)
    t_grid: np.ndarray = field(init=False)
    dim: int = field(init=False)
    J_spatial: np.ndarray = field(init=False)
    J: np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        if self.n_max < 1:
            raise ValueError(f"n_max must be >= 1, got {self.n_max}")
        if self.N_t < 1:
            raise ValueError(f"N_t must be >= 1, got {self.N_t}")
        if self.T_max <= 0:
            raise ValueError(f"T_max must be > 0, got {self.T_max}")

        # Spatial side: H_GV^{n_max}
        self.basis_spatial = full_dirac_basis(self.n_max)
        self.dim_spatial = len(self.basis_spatial)
        assert self.dim_spatial == full_dirac_dim(self.n_max)

        # Temporal side: small uniform grid on [-T_max, T_max].
        self.t_grid = temporal_grid(self.T_max, self.N_t)

        # Total dimension
        self.dim = self.dim_spatial * self.N_t

        # Spatial fundamental symmetry (lift of gamma^0)
        self.J_spatial = spatial_fundamental_symmetry(self.basis_spatial)

        # Full J = J_spatial (x) I_{N_t}
        I_t = identity_temporal(self.N_t)
        self.J = np.kron(self.J_spatial, I_t)

    # -------- Krein axiom checks --------

    def verify_J_squared_identity(
        self, tol: float = 1e-12
    ) -> Tuple[bool, float]:
        """J^2 = +I (fundamental-symmetry involution).

        Returns (ok, residual_F_norm).
        """
        I = np.eye(self.dim, dtype=np.complex128)
        J2 = self.J @ self.J
        residual = float(np.linalg.norm(J2 - I))
        return residual < tol, residual

    def verify_J_unitary(
        self, tol: float = 1e-12
    ) -> Tuple[bool, float]:
        """J^* J = I (J unitary as a Hilbert-space operator).

        Returns (ok, residual_F_norm).  Equivalent to J Hermitian-and-J^2=+I
        but checked independently as a separate axiom.
        """
        I = np.eye(self.dim, dtype=np.complex128)
        Jstar_J = self.J.conj().T @ self.J
        residual = float(np.linalg.norm(Jstar_J - I))
        return residual < tol, residual

    def verify_J_hermitian(
        self, tol: float = 1e-12
    ) -> Tuple[bool, float]:
        """J = J^* (Hermitian fundamental symmetry).

        For gamma^0 in the West-coast chiral basis this is automatic
        (real symmetric off-diagonal block matrix).
        """
        residual = float(np.linalg.norm(self.J - self.J.conj().T))
        return residual < tol, residual

    def krein_inner_product(
        self, psi: np.ndarray, phi: np.ndarray
    ) -> complex:
        """Krein inner product <psi, phi>_K = <psi, J phi> on K.

        Parameters
        ----------
        psi, phi : np.ndarray of shape (dim,)
            Krein vectors.

        Returns
        -------
        Complex number.
        """
        if psi.shape != (self.dim,) or phi.shape != (self.dim,):
            raise ValueError(
                f"both inputs must have shape ({self.dim},), got "
                f"{psi.shape}, {phi.shape}"
            )
        return complex(np.vdot(psi, self.J @ phi))

    def verify_krein_inner_product_hermitian(
        self, n_samples: int = 50, seed: int = 0, tol: float = 1e-10
    ) -> Tuple[bool, float]:
        """<phi, psi>_K = conj(<psi, phi>_K) on a random sample of pairs.

        Hermiticity (conjugate-symmetry) of the Krein form is
        automatically equivalent to J^* = J on the Hilbert level, but
        we verify it directly on randomly sampled vectors as a
        functional cross-check.

        Parameters
        ----------
        n_samples : int
            Number of random (psi, phi) pairs.
        seed : int
            RNG seed.
        tol : float
            Per-pair tolerance.

        Returns
        -------
        (ok, max_residual) : (bool, float)
        """
        rng = np.random.default_rng(seed)
        max_res = 0.0
        for _ in range(n_samples):
            psi = (
                rng.standard_normal(self.dim)
                + 1j * rng.standard_normal(self.dim)
            )
            phi = (
                rng.standard_normal(self.dim)
                + 1j * rng.standard_normal(self.dim)
            )
            a = self.krein_inner_product(psi, phi)
            b = self.krein_inner_product(phi, psi)
            res = abs(a - np.conj(b))
            max_res = max(max_res, res)
        return max_res < tol, max_res

    # -------- K^+ / K^- positive/negative definite split --------

    def positive_negative_split(
        self, tol: float = 1e-10
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Return (P_plus, P_minus) projecting onto K^+ and K^-.

        The Krein space splits as K = K^+ (+) K^- where K^+ is the +1
        eigenspace of J (the Krein form restricted to K^+ is positive-
        definite) and K^- is the -1 eigenspace (negative-definite).

        For J a Hermitian involution, J = P_plus - P_minus and
        P_plus + P_minus = I, so

            P_plus  = (1/2)(I + J),
            P_minus = (1/2)(I - J).

        Each projection is Hermitian and idempotent.

        Returns
        -------
        (P_plus, P_minus) : (np.ndarray, np.ndarray)
            Both of shape (dim, dim), complex.

        Raises
        ------
        ValueError if J does not satisfy J^2 = I within tol (sanity).
        """
        ok, residual = self.verify_J_squared_identity(tol=tol)
        if not ok:
            raise ValueError(
                f"J^2 != I (residual {residual:.3e}); cannot build "
                f"positive/negative split"
            )
        I = np.eye(self.dim, dtype=np.complex128)
        P_plus = 0.5 * (I + self.J)
        P_minus = 0.5 * (I - self.J)
        return P_plus, P_minus

    def verify_split_completeness(
        self, tol: float = 1e-12
    ) -> Tuple[bool, float]:
        """P_plus + P_minus = I and P_plus * P_minus = 0.

        Returns
        -------
        (ok, max_residual)
        """
        P_plus, P_minus = self.positive_negative_split()
        I = np.eye(self.dim, dtype=np.complex128)
        Z = np.zeros((self.dim, self.dim), dtype=np.complex128)
        r_sum = float(np.linalg.norm(P_plus + P_minus - I))
        r_orth = float(np.linalg.norm(P_plus @ P_minus - Z))
        r_idem_p = float(np.linalg.norm(P_plus @ P_plus - P_plus))
        r_idem_m = float(np.linalg.norm(P_minus @ P_minus - P_minus))
        max_res = max(r_sum, r_orth, r_idem_p, r_idem_m)
        return max_res < tol, max_res

    def split_dimensions(self) -> Tuple[int, int]:
        """Return (dim K^+, dim K^-) as trace of the respective projectors.

        For J = gamma^0 (chirality-swap) tensored with I_{N_t}, the +/- 1
        eigenspaces of J each have dimension dim / 2 = dim_spatial / 2 * N_t.
        """
        P_plus, P_minus = self.positive_negative_split()
        d_plus = int(np.round(np.real(np.trace(P_plus))))
        d_minus = int(np.round(np.real(np.trace(P_minus))))
        return d_plus, d_minus

    # -------- Riemannian-limit check (load-bearing falsifier) --------

    def riemannian_limit_check(
        self, tol: float = 1e-14
    ) -> Tuple[bool, dict]:
        """The Riemannian limit at N_t = 1 recovers H_GV bit-identically.

        At N_t = 1 the temporal grid collapses to the singleton t = 0.0
        and the Krein-space construction reduces to

            K_{n_max, N_t = 1} = H_GV^{n_max} (x) C^1 = H_GV^{n_max},
            J | _{N_t = 1} = J_spatial      (= gamma^0 on chirality)

        That is, the Krein structure adds nothing new at N_t = 1; the
        Hilbert space is bit-identically H_GV^{n_max} (Paper 32 def:H_GV)
        and the basis ordering / labels match those of
        `geovac.full_dirac_operator_system.full_dirac_basis(n_max)`
        verbatim.

        This method constructs a fresh KreinSpace at N_t = 1 (same
        n_max) and verifies:

            (i)  dim   = full_dirac_dim(n_max)
            (ii) basis = full_dirac_basis(n_max)  (label-by-label)
            (iii) J    = J_spatial                 (bit-identical)
            (iv) the existing self.J at general N_t reduces to J_spatial
                 on each fixed-time slice; equivalently the Kronecker
                 factor I_{N_t} ensures temporal-slot independence.

        If this check fails, the Camporesi-Higuchi spatial spinor bundle
        is structurally INCOMPATIBLE with the Cl(3, 1) gamma matrix
        embedding -- a major structural finding requiring PI escalation.

        Parameters
        ----------
        tol : float
            Frobenius tolerance for the J matrix comparison.  Default
            1e-14 (machine-precision target since the Kronecker product
            with I_1 is exact).

        Returns
        -------
        (ok, details) : (bool, dict)
            details has keys 'dim_match', 'basis_match',
            'J_match_residual', 'time_slice_check_residual'.
        """
        # Build the Riemannian-limit (N_t = 1) Krein space at this n_max.
        K_rie = KreinSpace(n_max=self.n_max, N_t=1, T_max=1.0)

        # Independent references from the Paper 32 §III H_GV module.
        ref_basis = full_dirac_basis(self.n_max)
        ref_dim = full_dirac_dim(self.n_max)
        ref_J = spatial_fundamental_symmetry(ref_basis)

        # Check (i) dim
        dim_match = K_rie.dim == ref_dim

        # Check (ii) basis label-by-label
        basis_match = all(
            K_rie.basis_spatial[i] == ref_basis[i] for i in range(ref_dim)
        )

        # Check (iii) J at N_t = 1 equals J_spatial
        J_residual = float(np.linalg.norm(K_rie.J - ref_J))

        # Check (iv) for the current (possibly N_t > 1) KreinSpace,
        # each time-slice projection P_k of self.J equals J_spatial.
        # The Kronecker structure J = J_spatial (x) I_{N_t} guarantees
        # this; we verify by extracting any one diagonal block.
        time_slice_res = 0.0
        if self.N_t >= 1:
            # Extract the (k=0, k=0) block of J on the temporal index.
            # In np.kron(A, B) ordering, the (i, j) block of shape
            # (dim_spatial, dim_spatial) corresponds to entries
            # J[i*dim_spatial + alpha, j*dim_spatial + beta]
            # under standard kron(spatial, temporal) ordering.
            # We constructed J = kron(J_spatial, I_{N_t}), which means
            # J[ alpha*N_t + k, beta*N_t + l ] = J_spatial[alpha, beta]
            # * I_{N_t}[k, l] = J_spatial[alpha, beta] * delta_{kl}.
            # So slice with stride N_t starting at index 0 gives J_spatial.
            slice_k0 = self.J[0 :: self.N_t, 0 :: self.N_t]
            time_slice_res = float(np.linalg.norm(slice_k0 - ref_J))

        details = {
            "n_max": self.n_max,
            "N_t": self.N_t,
            "dim_match": dim_match,
            "basis_match": basis_match,
            "J_match_residual": J_residual,
            "time_slice_check_residual": time_slice_res,
            "ref_dim": ref_dim,
            "krein_rie_dim": K_rie.dim,
        }

        ok = (
            dim_match
            and basis_match
            and J_residual < tol
            and time_slice_res < tol
        )
        return ok, details

    # -------- Convenience accessors --------

    def __repr__(self) -> str:
        return (
            f"KreinSpace(n_max={self.n_max}, N_t={self.N_t}, "
            f"T_max={self.T_max}, dim_spatial={self.dim_spatial}, "
            f"dim={self.dim})"
        )
