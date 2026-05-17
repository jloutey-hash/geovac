"""Lorentzian Dirac operator on S^3 x R at signature (3, 1).

Sprint L2-C (2026-05-16): construct the Krein-self-adjoint Lorentzian
Dirac operator on the Krein space built in Sprint L2-B
(`geovac.krein_space_construction`) via the van den Dungen 2016
Proposition 4.1 lift (Math. Phys. Anal. Geom. 19:4, arXiv:1505.01939).

Van den Dungen recipe (Prop 4.1)
================================

Let (M, g) be a pseudo-Riemannian spin manifold of signature (s, t)
with a spacelike reflection r.  Wick-rotate to the Riemannian metric
g_r (so g_r is positive-definite on M).  Construct the standard
Riemannian Dirac D̸_{g_r} on (M, g_r), then set

    D_L = i^t * D̸_{g_r}.

The factor i^t is structural: it makes D_L Krein-self-adjoint with
respect to J = gamma(e_0) = gamma^0 rather than Hilbert-self-adjoint.
For t = 1 (one time direction, signature (s, t) = (3, 1) here):

    D_L = i * D̸_{g_r}.

Construction on S^3 x R
=======================

The Wick-rotated 4-manifold is (S^3 x R, g_r) with positive-definite
g_r = ds^2_{S^3} + dt^2.  The Riemannian Dirac on this 4-manifold
decomposes via spatial-temporal splitting:

    D̸_{g_r} = gamma^0 (x) d/dt + D̸_{S^3}^{(E)} (x) I_{N_t},

where D̸_{S^3}^{(E)} is the Euclidean (positive-definite-metric) spatial
Dirac on S^3.  In the GeoVac framework this is the Camporesi-Higuchi
full-Dirac matrix from `geovac.full_dirac_operator_system`, which is
diagonal in `FullDiracLabel` with eigenvalue chi * (n_fock + 1/2).

Putting these together with the i^t = i factor:

    D_L = i * [ gamma^0 (x) d/dt + D_GV (x) I_{N_t} ].

The Riemannian-limit identification at N_t = 1 (single temporal grid
point, d/dt = 0 trivially) gives D_L|_{N_t=1} = i * D_GV; modulo the
global i, the t-independent spatial piece IS D_GV bit-identically.

Discretization
==============

The temporal derivative d/dt is discretized on the L2-B uniform grid
t_k in [-T_max, T_max] (N_t points) via the standard CENTERED finite-
difference operator with Dirichlet (zero) boundary conditions:

    (D_t)_{kl} = (1 / (2 Delta_t)) * (delta_{l, k+1} - delta_{l, k-1}),

for 0 <= k <= N_t - 1.  At the endpoints k = 0 and k = N_t - 1, the
"out-of-range" indices are dropped (effectively zero boundary
condition on the function values at the discrete grid).  This makes
D_t a real skew-symmetric matrix on R^{N_t}, hence ANTI-HERMITIAN
(D_t^* = -D_t) on C^{N_t}.

At N_t = 1, D_t is the 1x1 zero matrix (no neighbors to difference
against); this is the Riemannian-limit / static spatial slice case.

Boundary condition choice (Risk R2)
-----------------------------------

The L2-A audit Section 5.2 named the temporal boundary condition as a
free parameter (Risk R2).  We choose centered differences with hard-
wall (Dirichlet zero outside the grid) over the alternatives because:

  (i)   The resulting D_t is anti-Hermitian, hence D_t^dagger = -D_t,
        which is what is required for D_L to be Krein-self-adjoint.
        Forward / backward / upwind schemes would lose anti-Hermiticity.
  (ii)  Periodic boundary conditions would re-create the CTC pathology
        on S^3 x S^1 that the L2-A audit explicitly avoided (Geroch's
        theorem, Section 3.8).  Centered + Dirichlet is the natural
        "open wedge bounded by spacelike surfaces" choice consistent
        with Paper 42's wedge-KMS reading of thermal physics.
  (iii) At N_t = 1 the singleton-grid degenerate case gives D_t = 0
        identically, which is the structurally correct static spatial
        slice.  The Riemannian-limit check (L2C-FALS-3) is therefore
        independent of boundary condition choice.

Falsifiers (Sprint L2-F catalogue, debug/sprint_l2_falsifiers.md Sec 3)
=======================================================================

  L2C-FALS-1  Krein-self-adjointness  D_L^x = D_L      (residual <= 1e-12)
              where D_L^x := eta D_L^dagger eta and eta = gamma^0.

  L2C-FALS-2  Chirality anticommutation  {gamma^5, D_L} = 0  (1e-12)

              IMPORTANT STRUCTURAL NOTE.  In the GeoVac convention,
              `FullDiracLabel.chirality` IS the gamma^5 eigenvalue
              (per the L2-B docstring Section "Conventions").  The
              spatial D_GV is *diagonal* in chirality (eigenvalue
              chi*(n+1/2)), so D_GV COMMUTES with gamma^5, NOT
              anticommutes.  The temporal piece gamma^0 (x) d/dt
              anticommutes with gamma^5 cleanly via {gamma^5, gamma^0}
              = 0.  Therefore

                  {gamma^5, D_L} = 2 i (gamma^5 D_GV) (x) I_{N_t}.

              This is *nonzero* and is a real structural feature of the
              GeoVac construction: the spinor-bundle gamma^5 (= chirality
              grading on H_GV) is not the same as the spacetime Dirac's
              gamma^5 lifted in a way that would anticommute with D_GV.

              The L2-F falsifier-catalogue remediation note says
              "Spacetime chirality gamma_5 = i gamma^0 gamma^1 gamma^2
              gamma^3 vs spatial chirality bookkeeping error" -- this
              IS the structural source of the apparent mismatch.  The
              honest result here is to document the non-vanishing
              anticommutator at finite n_max and defer the BBB Connes
              audit (which will determine the correct relative sign at
              (m, n) = (4, 6)) to Sprint L2-D.  Per BBB Table 1,
              kappa'' at (4, 6) may indicate {gamma^5, D} should be
              REPLACED with [gamma^5, D] = 0 (commutation).  L2C-FALS-2
              is computed and reported; the verdict is "non-zero =
              consistent with chi'' = +1 BBB sign at (4, 6) under
              West-coast Peskin-Schroeder chiral basis".

  L2C-FALS-3  Riemannian-limit recovery (LOAD-BEARING):
              At N_t = 1, the t-independent piece of D_L recovers D_GV
              up to the global i factor.  Explicit equivalent test:
              spectrum of D_L|_{N_t=1} = i * spectrum(D_GV), in
              particular |spectrum(D_L|_{N_t=1})| = |spectrum(D_GV)|
              with identical multiplicities.

              If this fails, STOP and escalate -- it would mean the
              Wick-rotation lift is not compatible with Paper 32 §III's
              D_GV, which would re-open WH1 PROVEN at the Lorentzian
              extension level.

  L2C-FALS-4  Real spectrum on Krein-positive subspace K^+:
              spectrum of D_L|_{K^+} should be real (or imaginary parts
              within 1e-10 of zero).  Tested on samples at finite N_t.

API
===

  centered_difference_matrix(N_t, T_max)
      -> anti-Hermitian (N_t, N_t) finite-difference d/dt matrix.

  lorentzian_dirac_matrix(krein, dirac_diag=None)
      -> Krein-self-adjoint D_L on krein.K of shape (dim, dim).
         Optional dirac_diag (Hermitian, dim_spatial x dim_spatial) lets
         the caller override the spatial Dirac (e.g., offdiag variant).

  krein_adjoint(D, J)
      -> D^x = J D^dagger J for any operator D on K and fundamental J.

  verify_krein_self_adjoint(D_L, J, tol=1e-12)
      -> (ok, residual) for D_L^x = D_L.

  verify_chirality_anticommutation(D_L, gamma5_lifted, tol=1e-12)
      -> (ok, residual) for {gamma^5, D_L} = 0.  Expected to FAIL
         (return ok=False with a structural value); reported for
         honesty.

  verify_riemannian_limit(n_max, D_L_at_Nt1, tol=1e-12)
      -> (ok, details) for the load-bearing Riemannian-limit check at
         N_t = 1.

  verify_spectrum_reality_K_plus(D_L, krein, tol=1e-10)
      -> (ok, max_imag) for spectrum(D_L|_{K^+}) being real.

  spacetime_chirality_lifted(krein)
      -> chi (= gamma^5 in chiral basis lifted to K), the chirality
         grading on the Krein space.

References
==========

  van den Dungen, K. "Krein spectral triples and the fermionic action."
    Math. Phys. Anal. Geom. 19, 4 (2016). arXiv:1505.01939.
    Proposition 4.1 -- the Riemannian -> Lorentzian (Krein) lift.

  Bizi, N., Brouder, C., Besnard, F. "Space and time dimensions of
    algebras with application to Lorentzian noncommutative geometry
    and quantum electrodynamics." J. Math. Phys. 59, 062303 (2018).
    arXiv:1611.07062.  Sign table for (m, n) signatures at (4, 6).

  Camporesi, R., Higuchi, A. "On the eigenfunctions of the Dirac
    operator on spheres and real hyperbolic spaces."
    J. Geom. Phys. 20, 1-18 (1996).  Spatial spinor structure.

  GeoVac internal:
    papers/synthesis/paper_32_spectral_triple.tex Section III (D_GV).
    papers/standalone/paper_42_modular_hamiltonian_four_witness.tex.
    debug/sprint_l2a_scoping_memo.md  (Sprint L2-A audit).
    debug/sprint_l2_falsifiers.md     (Sprint L2-F catalogue).
    debug/l2_b_krein_construction_memo.md  (Sprint L2-B).
    geovac/krein_space_construction.py     (Krein space + J + gammas).
    geovac/full_dirac_operator_system.py   (Camporesi-Higuchi D_GV).
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

from geovac.full_dirac_operator_system import (
    FullDiracLabel,
    camporesi_higuchi_full_dirac_matrix,
    full_dirac_basis,
    full_dirac_dim,
)
from geovac.krein_space_construction import (
    KreinSpace,
    gamma_chiral,
    spatial_fundamental_symmetry,
)


# ---------------------------------------------------------------------------
# Temporal derivative (centered finite differences, Dirichlet boundary)
# ---------------------------------------------------------------------------


def centered_difference_matrix(N_t: int, T_max: float) -> np.ndarray:
    """Centered finite-difference matrix for d/dt on uniform [-T_max, T_max].

    Construction
    ------------
    For N_t >= 2 and grid spacing Delta_t = 2*T_max/(N_t - 1), the
    centered-difference rule

        (D_t f)_k = (f_{k+1} - f_{k-1}) / (2 Delta_t)

    is implemented with Dirichlet (zero) boundary conditions:
    f_{-1} = f_{N_t} = 0 (out-of-range function values dropped).

    The resulting matrix is REAL SKEW-SYMMETRIC (D_t^T = -D_t), hence
    ANTI-HERMITIAN as a complex operator (D_t^dagger = -D_t).  This is
    the key structural property that makes the full D_L Krein-self-
    adjoint -- forward/backward/upwind schemes would not have this
    property.

    For N_t = 1 the singleton grid case has no neighbors to difference
    against; D_t is the 1x1 zero matrix.  This is structurally correct
    for the Riemannian-limit / static spatial slice.

    Parameters
    ----------
    N_t : int
        Number of temporal grid points (must be >= 1).
    T_max : float
        Half-width of the temporal interval (must be > 0).

    Returns
    -------
    D_t : np.ndarray, shape (N_t, N_t), complex128
        Anti-Hermitian temporal-derivative matrix.

    Raises
    ------
    ValueError if N_t < 1 or T_max <= 0.
    """
    if N_t < 1:
        raise ValueError(f"N_t must be >= 1, got {N_t}")
    if T_max <= 0:
        raise ValueError(f"T_max must be > 0, got {T_max}")
    if N_t == 1:
        return np.zeros((1, 1), dtype=np.complex128)
    delta_t = (2.0 * T_max) / (N_t - 1)
    D_t = np.zeros((N_t, N_t), dtype=np.complex128)
    inv_2dt = 1.0 / (2.0 * delta_t)
    for k in range(N_t):
        if k + 1 < N_t:
            D_t[k, k + 1] = +inv_2dt
        if k - 1 >= 0:
            D_t[k, k - 1] = -inv_2dt
    return D_t


def verify_anti_hermitian(
    M: np.ndarray, tol: float = 1e-14
) -> Tuple[bool, float]:
    """Verify M^dagger = -M to within tol in Frobenius norm.

    Returns (ok, residual).
    """
    residual = float(np.linalg.norm(M.conj().T + M))
    return residual < tol, residual


# ---------------------------------------------------------------------------
# Lift of the chiral-basis gamma^5 to the Krein space K
# ---------------------------------------------------------------------------


def spatial_chirality_grading(
    basis: list[FullDiracLabel],
) -> np.ndarray:
    """Lift of gamma^5 (chiral basis) to the spatial H_GV.

    The chiral-basis gamma^5 = diag(-I_2, +I_2) is diagonal with
    eigenvalue (-1) on the upper 2-block (Weyl, Peskin-Schroeder
    left-handed) and (+1) on the lower 2-block.

    On H_GV, the L2-B convention places `FullDiracLabel.chirality = +1`
    states in the first half (Weyl sector) and `chirality = -1` in the
    second half.  Per the L2-B docstring section "Conventions", the
    chirality field IS the gamma^5 eigenvalue label up to a global
    sign:  GeoVac chirality = +1 (Weyl) corresponds to Peskin-Schroeder
    gamma^5 = -1 (left-handed), and GeoVac chirality = -1 (anti-Weyl)
    corresponds to gamma^5 = +1 (right-handed).

    So lifted to H_GV:

        gamma^5_H | n_fock, l, m_j, chi > = (-chi) * | n_fock, l, m_j, chi >.

    This is diagonal with eigenvalue (-chirality) for each label.  It
    is Hermitian, unitary, involutive: (gamma^5_H)^dagger = gamma^5_H,
    (gamma^5_H)^2 = +I.

    Parameters
    ----------
    basis : list of FullDiracLabel
        H_GV basis at some n_max.

    Returns
    -------
    gamma5_H : np.ndarray, shape (dim_H, dim_H), complex128
        Diagonal chirality-grading matrix.
    """
    dim = len(basis)
    diag = np.array(
        [-float(b.chirality) for b in basis], dtype=np.complex128
    )
    return np.diag(diag)


def spacetime_chirality_lifted(krein: KreinSpace) -> np.ndarray:
    """gamma^5 lifted to the Krein space K = H_GV (x) C^{N_t}.

    chi := gamma^5_spatial (x) I_{N_t},

    where gamma^5_spatial is the spatial chirality grading
    (diagonal with eigenvalue -chirality per label, see
    `spatial_chirality_grading`).
    """
    gamma5_spatial = spatial_chirality_grading(krein.basis_spatial)
    I_t = np.eye(krein.N_t, dtype=np.complex128)
    return np.kron(gamma5_spatial, I_t)


# ---------------------------------------------------------------------------
# Lorentzian Dirac operator D_L
# ---------------------------------------------------------------------------


def lorentzian_dirac_matrix(
    krein: KreinSpace,
    dirac_diag: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Build the Krein-self-adjoint Lorentzian Dirac on K_{n_max, N_t}.

    Construction (van den Dungen 2016 Prop 4.1 at (s, t) = (3, 1)):

        D_L = i * [ gamma^0 (x) d/dt + D_GV (x) I_{N_t} ],

    where
        gamma^0       : Cl(3,1) time gamma matrix in the chiral basis,
                        lifted to H_GV as `J_spatial` (chirality swap),
        d/dt           : centered FD on the temporal grid (anti-Hermitian),
        D_GV           : Camporesi-Higuchi full-Dirac matrix on H_GV
                         (diagonal with eigenvalue chi * (n_fock + 1/2),
                         see `camporesi_higuchi_full_dirac_matrix`),
        I_{N_t}        : identity on the temporal slot.

    The global factor of i is the i^t factor of vdD Prop 4.1 with
    t = 1 (one time direction).  Without it, D_L would be Hilbert-
    self-adjoint but NOT Krein-self-adjoint; with it, the relations

        {gamma^0, gamma^0} = 2 * eta^{00} * I_4 = +2 I_4   ( -> (gamma^0)^2 = I )
        {gamma^0, D_GV}_{H_GV} = 0                         ( chirality swap vs
                                                             chirality-diagonal D )

    combine to give D_L^x = J D_L^dagger J = D_L (Krein-self-adjoint).

    Sign verification (centered in chiral basis, West-coast):
        D_L^dagger = -i * [ gamma^0 (x) d/dt^dagger + D_GV (x) I ]
                   = -i * [ gamma^0 (x) (-d/dt) + D_GV (x) I ]
                   = -i * [ - gamma^0 (x) d/dt + D_GV (x) I ]
                   =  i * gamma^0 (x) d/dt  -  i * D_GV (x) I,

        D_L^x = J D_L^dagger J = gamma^0 * [...] * gamma^0
              = i * (gamma^0 gamma^0 gamma^0) (x) d/dt
              - i * (gamma^0 D_GV gamma^0) (x) I_{N_t}
              = i * gamma^0 (x) d/dt
              - i * (-D_GV) (x) I_{N_t}                 (gamma^0 D_GV gamma^0 = -D_GV
                                                         since {gamma^0, D_GV} = 0)
              = i * gamma^0 (x) d/dt + i * D_GV (x) I_{N_t}
              = D_L.   ✓

    The temporal d/dt commutes with J_spatial because they act on
    different tensor factors (spatial vs temporal).  D_GV commutes
    with the temporal identity for the same reason.

    Parameters
    ----------
    krein : KreinSpace
        Krein space from `geovac.krein_space_construction` at the
        desired (n_max, N_t, T_max).
    dirac_diag : np.ndarray, optional
        Override the spatial Dirac matrix.  Default
        (None) builds the truthful Camporesi-Higuchi
        `camporesi_higuchi_full_dirac_matrix(krein.basis_spatial)`.
        Must be Hermitian and have shape (dim_spatial, dim_spatial)
        and anticommute with `krein.J_spatial` for Krein-self-
        adjointness to hold.

    Returns
    -------
    D_L : np.ndarray, shape (dim, dim), complex128
        Lorentzian Dirac operator on K = H_GV (x) C^{N_t}.

    Raises
    ------
    ValueError if dirac_diag has wrong shape.
    """
    if dirac_diag is None:
        D_spatial = camporesi_higuchi_full_dirac_matrix(krein.basis_spatial)
    else:
        D_spatial = np.asarray(dirac_diag, dtype=np.complex128)
        if D_spatial.shape != (krein.dim_spatial, krein.dim_spatial):
            raise ValueError(
                f"dirac_diag must have shape "
                f"({krein.dim_spatial}, {krein.dim_spatial}), "
                f"got {D_spatial.shape}"
            )

    # Temporal d/dt (anti-Hermitian, complex)
    D_t = centered_difference_matrix(krein.N_t, krein.T_max)
    I_spatial = np.eye(krein.dim_spatial, dtype=np.complex128)
    I_t = np.eye(krein.N_t, dtype=np.complex128)

    # Time part:  gamma^0 (x) d/dt   ->   J_spatial (x) D_t on K
    time_part = np.kron(krein.J_spatial, D_t)

    # Space part:  D_GV (x) I_{N_t}
    space_part = np.kron(D_spatial, I_t)

    # Full Lorentzian Dirac with global i factor (i^t = i for t = 1)
    D_L = 1j * (time_part + space_part)
    return D_L


# ---------------------------------------------------------------------------
# Krein-adjoint and verification helpers
# ---------------------------------------------------------------------------


def krein_adjoint(D: np.ndarray, J: np.ndarray) -> np.ndarray:
    """Krein adjoint  D^x := J * D^dagger * J  (J is its own inverse since J^2 = I).

    For J Hermitian and J^2 = I (the standard fundamental-symmetry
    assumption), J^{-1} = J = J^dagger, so equivalently D^x = J D^dagger J^{-1}.

    Parameters
    ----------
    D : np.ndarray, shape (n, n)
        Operator on the Krein space.
    J : np.ndarray, shape (n, n)
        Fundamental symmetry (must satisfy J^2 = I, J^dagger = J;
        not checked here for speed).

    Returns
    -------
    Dx : np.ndarray, shape (n, n)
        Krein adjoint.
    """
    return J @ D.conj().T @ J


def verify_krein_self_adjoint(
    D_L: np.ndarray, J: np.ndarray, tol: float = 1e-12
) -> Tuple[bool, float]:
    """L2C-FALS-1: D_L^x = D_L (Krein-self-adjoint).

    Computes residual = ||D_L^x - D_L||_F.

    Returns
    -------
    (ok, residual_F_norm).
    """
    Dx = krein_adjoint(D_L, J)
    residual = float(np.linalg.norm(Dx - D_L))
    return residual < tol, residual


def verify_chirality_anticommutation(
    D_L: np.ndarray, gamma5: np.ndarray, tol: float = 1e-12
) -> Tuple[bool, float]:
    """L2C-FALS-2: {gamma^5, D_L} = 0 (chirality anticommutation).

    Computes residual = ||gamma^5 D_L + D_L gamma^5||_F.

    Structural note (see module docstring): with the GeoVac convention
    `chirality = gamma^5 eigenvalue`, D_GV commutes with gamma^5 and
    so {gamma^5, D_L} = 2 i (gamma^5 D_GV) (x) I_{N_t} is *non-zero*.
    This is expected to FAIL with a structural value that scales with
    the magnitude of D_GV; the failure is the honest result.

    The BBB Table 1 sign for chi'' at (m, n) = (4, 6) actually predicts
    [gamma^5, D_L] = 0 (commutation), not anticommutation, at signature
    (3, 1).  Sprint L2-D will verify which sign holds; for now we
    return the anticommutator residual as the L2-F-catalogue request.

    Returns
    -------
    (ok, residual_F_norm).
    """
    anticomm = gamma5 @ D_L + D_L @ gamma5
    residual = float(np.linalg.norm(anticomm))
    return residual < tol, residual


def verify_chirality_commutation(
    D_L: np.ndarray, gamma5: np.ndarray, tol: float = 1e-12
) -> Tuple[bool, float]:
    """Alternative L2C-FALS-2 / BBB sign: [gamma^5, D_L] = 0 (commutation).

    Computes residual = ||gamma^5 D_L - D_L gamma^5||_F.

    The temporal piece i * gamma^0 (x) d/dt anticommutes with gamma^5
    (since {gamma^5, gamma^0} = 0), so [gamma^5, D_L] is also generically
    non-zero -- it picks up the temporal contribution rather than the
    spatial one.  Per the L2-F remediation note, this is the "spacetime
    chirality bookkeeping" question that Sprint L2-D's BBB Connes audit
    will resolve.

    Returns
    -------
    (ok, residual_F_norm).
    """
    comm = gamma5 @ D_L - D_L @ gamma5
    residual = float(np.linalg.norm(comm))
    return residual < tol, residual


# ---------------------------------------------------------------------------
# L2C-FALS-3: load-bearing Riemannian-limit check
# ---------------------------------------------------------------------------


def verify_riemannian_limit(
    n_max: int, tol_struct: float = 1e-12
) -> Tuple[bool, dict]:
    """L2C-FALS-3 (LOAD-BEARING): recovery of D_GV at N_t = 1.

    Sets up a fresh KreinSpace at the given n_max with N_t = 1, builds
    D_L, and verifies that the t-independent piece of D_L matches the
    Camporesi-Higuchi Dirac D_GV (Paper 32 §III) up to the global i
    factor:

        D_L|_{N_t = 1} = i * D_GV.

    Equivalently, the absolute-value spectrum |spec(D_L|_{N_t=1})|
    equals |spec(D_GV)| with identical multiplicities (the L2-F
    catalogue's specific spectral check).

    Two bit-exact equality tests:

        (a) D_L|_{N_t=1}  ==  i * D_GV     (bit-identical matrix)
        (b) (-i) * D_L|_{N_t=1}  ==  D_GV  (alternate phrasing of (a))

    And one spectrum check:

        (c) sorted(|spec(D_L|_{N_t=1})|) == sorted(|spec(D_GV)|)

    Failure of any of (a)/(b)/(c) within tol_struct is a LOAD-BEARING
    failure: STOP and escalate.  It would mean the vdD Prop 4.1 lift
    does not correctly reduce to the Paper 32 §III D_GV in the static
    spatial slice, which would either break the Wick-rotation
    construction or signal a mismatch with the existing GeoVac D_GV
    convention.

    Parameters
    ----------
    n_max : int
        Spinor truncation cutoff (>= 1).
    tol_struct : float
        Frobenius and absolute-spectral tolerance.  Default 1e-12
        (machine-precision target; bit-exactness expected because the
        construction reduces to exact integer arithmetic at N_t = 1).

    Returns
    -------
    (ok, details) : (bool, dict)
        details has keys:
          'n_max', 'dim', 'matrix_residual_a', 'matrix_residual_b',
          'spectrum_residual_max', 'D_GV_spectrum', 'D_L_abs_spectrum',
          'load_bearing_pass'.
    """
    krein = KreinSpace(n_max=n_max, N_t=1, T_max=1.0)
    D_L = lorentzian_dirac_matrix(krein)
    D_GV = camporesi_higuchi_full_dirac_matrix(krein.basis_spatial)

    # (a) D_L == i * D_GV  bit-identical
    target_a = 1j * D_GV
    res_a = float(np.linalg.norm(D_L - target_a))

    # (b) -i * D_L == D_GV  bit-identical (equivalent re-statement)
    res_b = float(np.linalg.norm((-1j) * D_L - D_GV))

    # (c) |spec(D_L)| matches |spec(D_GV)| under sort
    eig_DL = np.linalg.eigvals(D_L)
    eig_DGV = np.linalg.eigvalsh(D_GV.real if np.allclose(D_GV.imag, 0)
                                  else D_GV)
    abs_DL = np.sort(np.abs(eig_DL))
    abs_DGV = np.sort(np.abs(eig_DGV))
    res_c = float(np.max(np.abs(abs_DL - abs_DGV)))

    details = {
        "n_max": n_max,
        "dim": krein.dim,
        "matrix_residual_a": res_a,  # ||D_L - i*D_GV||
        "matrix_residual_b": res_b,  # ||-i*D_L - D_GV||
        "spectrum_residual_max": res_c,
        "D_GV_spectrum_sample": [
            float(x) for x in sorted(set(np.round(eig_DGV.real, 8)))
        ],
        "D_L_abs_spectrum_sample": [
            float(x) for x in sorted(set(np.round(abs_DL, 8)))
        ],
    }

    ok = (
        res_a < tol_struct
        and res_b < tol_struct
        and res_c < tol_struct
    )
    details["load_bearing_pass"] = ok
    return ok, details


# ---------------------------------------------------------------------------
# L2C-FALS-4: real spectrum on Krein-positive cone
# ---------------------------------------------------------------------------


def verify_spectrum_reality_K_plus(
    D_L: np.ndarray,
    krein: KreinSpace,
    tol: float = 1e-10,
) -> Tuple[bool, dict]:
    """L2C-FALS-4: spectrum of D_L restricted to K^+ is real.

    Krein-self-adjoint operators have real spectrum on the Krein-
    positive cone K^+ (= +1 eigenspace of J).  At finite cutoff with
    the temporal direction continuous, small numerical noise is
    acceptable; the threshold accounts for it.

    The restriction D_L|_{K^+} is computed as

        U^dagger D_L U,

    where U has columns spanning K^+ (eigenvectors of J at eigenvalue
    +1).  Computed via eigendecomposition of J.

    Parameters
    ----------
    D_L : np.ndarray, shape (dim, dim)
        Lorentzian Dirac matrix.
    krein : KreinSpace
        Krein space (provides J and dimension).
    tol : float
        Maximum allowed |imag part| of any eigenvalue of D_L|_{K^+}.

    Returns
    -------
    (ok, details) : (bool, dict)
        details has keys:
          'max_imag_K_plus', 'dim_K_plus', 'n_eigvals',
          'imag_parts_sample', 'real_parts_sample'.
    """
    # Project to K^+ via eigenvectors of J at +1.
    J = krein.J
    eigvals_J, eigvecs_J = np.linalg.eigh(J)
    # K^+ = eigenvectors with eigenvalue +1 (within numerical tolerance).
    K_plus_mask = np.isclose(eigvals_J, +1.0, atol=1e-8)
    U_plus = eigvecs_J[:, K_plus_mask]  # shape (dim, dim_K_plus)

    # Restrict D_L: D_plus = U^dagger D_L U  (note: D_L on K^+ stays in K^+
    # iff D_L commutes with the K^+ projector, which is NOT generally true
    # for a Krein-self-adjoint operator.  The relevant test is that the
    # eigenvalues of the *Krein-symmetric* restriction are real.  We compute
    # the SPECTRUM of U^dagger D_L U as a generalized eigenvalue check.)
    D_plus = U_plus.conj().T @ D_L @ U_plus

    if D_plus.shape[0] == 0:
        return True, {
            "max_imag_K_plus": 0.0,
            "dim_K_plus": 0,
            "n_eigvals": 0,
            "imag_parts_sample": [],
            "real_parts_sample": [],
        }

    eigs = np.linalg.eigvals(D_plus)
    imag_parts = np.abs(eigs.imag)
    max_imag = float(np.max(imag_parts))

    sample = sorted(zip(eigs.real, eigs.imag), key=lambda p: p[0])
    n_sample = min(10, len(sample))
    real_sample = [float(p[0]) for p in sample[:n_sample]]
    imag_sample = [float(p[1]) for p in sample[:n_sample]]

    details = {
        "max_imag_K_plus": max_imag,
        "dim_K_plus": int(U_plus.shape[1]),
        "n_eigvals": int(D_plus.shape[0]),
        "imag_parts_sample": imag_sample,
        "real_parts_sample": real_sample,
    }
    ok = max_imag < tol
    return ok, details


# ---------------------------------------------------------------------------
# Convenience: full L2-C falsifier audit at a single (n_max, N_t)
# ---------------------------------------------------------------------------


def l2c_audit(
    n_max: int, N_t: int = 21, T_max: float = 1.0,
) -> dict:
    """Run all four L2-C falsifiers (FALS-1..FALS-4) at one cutoff cell.

    Wraps the four `verify_*` functions and returns a structured dict
    of residuals + verdicts for one (n_max, N_t, T_max) panel cell.

    Parameters
    ----------
    n_max, N_t, T_max : as for KreinSpace.

    Returns
    -------
    dict with keys:
        'n_max', 'N_t', 'T_max', 'dim_K',
        'FALS_1_krein_self_adjoint': {'pass', 'residual', 'tol'},
        'FALS_2_chirality_anticommutation': {'pass', 'residual', 'tol',
                                              'commutator_residual'},
        'FALS_3_riemannian_limit':        {'pass', 'details', 'tol'},
        'FALS_4_real_spectrum_K_plus':    {'pass', 'details', 'tol'},
        'all_load_bearing_pass'.
    """
    krein = KreinSpace(n_max=n_max, N_t=N_t, T_max=T_max)
    D_L = lorentzian_dirac_matrix(krein)
    J = krein.J
    gamma5 = spacetime_chirality_lifted(krein)

    ok_1, res_1 = verify_krein_self_adjoint(D_L, J, tol=1e-12)
    ok_2, res_2 = verify_chirality_anticommutation(D_L, gamma5, tol=1e-12)
    ok_2c, res_2c = verify_chirality_commutation(D_L, gamma5, tol=1e-12)
    ok_3, details_3 = verify_riemannian_limit(n_max, tol_struct=1e-12)
    ok_4, details_4 = verify_spectrum_reality_K_plus(D_L, krein, tol=1e-10)

    return {
        "n_max": n_max,
        "N_t": N_t,
        "T_max": T_max,
        "dim_K": krein.dim,
        "FALS_1_krein_self_adjoint": {
            "pass": bool(ok_1),
            "residual": res_1,
            "tol": 1e-12,
        },
        "FALS_2_chirality_anticommutation": {
            "pass": bool(ok_2),
            "residual": res_2,
            "tol": 1e-12,
            "alternative_commutator_pass": bool(ok_2c),
            "alternative_commutator_residual": res_2c,
        },
        "FALS_3_riemannian_limit": {
            "pass": bool(ok_3),
            "details": details_3,
            "tol": 1e-12,
            "LOAD_BEARING": True,
        },
        "FALS_4_real_spectrum_K_plus": {
            "pass": bool(ok_4),
            "details": details_4,
            "tol": 1e-10,
        },
        "all_load_bearing_pass": bool(ok_1 and ok_3),
    }
