"""Sprint F3 Step 1 — Cross-block h1 algebraic diagnostic.

Compute the missing cross-block h1 matrix elements at NaH R = R_eq = 3.566 bohr
using the F2 3D-quadrature infrastructure (closed-form theta integral for s-s).

The cross-block h1 element for orbital A on center A and orbital B on center B is:
    h1_cross[A,B] = <psi_A | T | psi_B>  + sum_C <psi_A | -Z_C/|r - R_C| | psi_B>

For axially symmetric NaH (Na at z=0, H at z=R), and two s-orbitals on different
centers, the matrix element factorizes via prolate spheroidal coordinates or
equivalently 3D Cartesian quadrature in axial symmetry.

We compute:
  (1) <Na 3s | V_ne(H) | H 1s>           [the off-diagonal cross-center attraction]
  (2) <Na 3s | V_ne(Na) | H 1s>          [the off-diagonal own-center attraction]
  (3) <Na 3s | -1/2 nabla^2 | H 1s>      [the kinetic cross-block element]
  (4) Total cross-block h1[Na 3s, H 1s] = (1) + (2) + (3)
  (5) Bonding/antibonding eigenvalues from a 2x2 diagonalization

Both Na 3s hydrogenic (Z_orb=1, the framework baseline) and Na 3s multi-zeta
variants are tested.

The 3D integral for two s-orbitals on different centers separated by R along z:
    psi_A(r_A) * V(r) * psi_B(r_B), where r_A = |r - 0|, r_B = |r - R z_hat|

Axial symmetry: use cylindrical coordinates (rho, z, phi) with phi-integration
giving 2*pi. The s-orbital phi/theta dependence is just 1/sqrt(4*pi). So:
    psi_A psi_B = (1/4pi) R_A(r_A) R_B(r_B)
    matrix element = (1/4pi) * 2*pi * integral_0^inf rho drho * integral_-inf^inf dz *
                     R_A(sqrt(rho^2 + z^2)) * R_B(sqrt(rho^2 + (z-R)^2)) * O(r)
                   = (1/2) * integral_0^inf rho drho * integral_-inf^inf dz * R_A R_B O(r)

For operators:
  V_ne(C) with C at z=Z_C_z: O(r) = -Z_C_charge / sqrt(rho^2 + (z - Z_C_z)^2)
  Kinetic: more involved (Laplacian), but we can use prolate spheroidal closed form
           OR compute via [T = E - V] applied to one side then integrate.

For the kinetic element on two H-like s-orbitals at different centers, the
classical molecular-orbital result gives a small finite value at moderate R.
At R = 3.566 bohr with Na 3s (Z_orb=1) and H 1s (Z_orb=1), both orbitals are
~hydrogenic-like-extended and overlap is non-negligible.

Practical kinetic-element computation:
  T |psi_B>:  if psi_B is hydrogenic at center B with Z_B, then
              -1/2 nabla^2 psi_B = (E_B + Z_B/r_B) psi_B = (-Z_B^2/(2 n_B^2) + Z_B/r_B) psi_B
  So <psi_A | T | psi_B> = E_B <psi_A | psi_B> + Z_B <psi_A | 1/r_B | psi_B>

The overlap <psi_A | psi_B> and <psi_A | 1/r_B | psi_B> are both pure radial-axial
quadrature integrals reducible to (rho, z) integration.

Sign convention check: We use H = T + V with T positive kinetic and V = -Z/r the
attractive electron-nuclear potential. h1[A, B] = T[A,B] + sum_C V_ne(C)[A, B] +
PK[A,B] (PK not relevant for cross-block at this level).

OUTPUT: Step 1 verdict per decision gate:
  - bonding favored by > 50 mHa: PROCEED to Step 2
  - bonding favored 0-50 mHa: PROCEED with caveat
  - bonding NOT favored: STOP, wall has another mechanism
"""

from __future__ import annotations

import json
import math
import time
from typing import Dict, Tuple

import numpy as np
from scipy.special import roots_legendre

from scipy.special import genlaguerre, factorial
from geovac.multi_zeta_orbitals import (
    MultiZetaOrbital,
    get_physical_valence_orbitals,
)


# ---------------------------------------------------------------------------
# Analytical hydrogenic radial wavefunctions (grid-independent normalization)
# ---------------------------------------------------------------------------

def hydrogenic_R_nl(Z: float, n: int, l: int, r: np.ndarray) -> np.ndarray:
    """Normalized hydrogenic radial wavefunction R_{nl}(r), analytical form.

    R_{nl}(r) = N * (2 Z r / n)^l * exp(-Z r / n) * L^{2l+1}_{n-l-1}(2 Z r / n)
    where N = sqrt[(2 Z / n)^3 * (n - l - 1)! / (2 n (n + l)!)]

    Normalized such that integral_0^inf R^2 r^2 dr = 1.
    """
    rho = 2.0 * Z * r / n
    L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    # Analytical normalization constant
    norm = math.sqrt(
        (2.0 * Z / n) ** 3
        * float(factorial(n - l - 1))
        / (2.0 * n * float(factorial(n + l)))
    )
    return norm * rho ** l * np.exp(-rho / 2.0) * L_poly


# ---------------------------------------------------------------------------
# Orbital evaluators (s-orbitals only for this diagnostic)
# ---------------------------------------------------------------------------

def na_3s_multizeta() -> MultiZetaOrbital:
    """Z=11 Na 3s multi-zeta (physical-fit) orbital."""
    orbs = get_physical_valence_orbitals(11)
    for o in orbs:
        if o.n_orbital == 3 and o.l_orbital == 0:
            return o
    raise RuntimeError("Na 3s not in physical-fit registry")


def na_3s_hydrogenic(r: np.ndarray, Z_orb: float = 1.0) -> np.ndarray:
    """Hydrogenic Na 3s at Z_orb=1 (framework baseline)."""
    return hydrogenic_R_nl(Z_orb, 3, 0, r)


def h_1s(r: np.ndarray, Z_orb: float = 1.0) -> np.ndarray:
    """Hydrogenic H 1s."""
    return hydrogenic_R_nl(Z_orb, 1, 0, r)


# ---------------------------------------------------------------------------
# 3D quadrature for cross-block s-s matrix elements via 2D (rho, z) reduction.
# ---------------------------------------------------------------------------

def _quadrature_grid(n_rho: int = 60, n_z: int = 80,
                     rho_max: float = 15.0, z_max: float = 15.0):
    """Build (rho, z, phi) -> (rho, z) 2D grid.

    Gauss-Legendre in (rho, z) on [0, rho_max] x [-z_max, z_max].
    phi-integration done analytically (2*pi factor) since all integrands
    are axially symmetric for s-orbitals.

    Returns:
        rho_pts (n_rho,), rho_w (n_rho,),
        z_pts (n_z,), z_w (n_z,)
    """
    u, wu = roots_legendre(n_rho)
    # rho in [0, rho_max]: rho = rho_max * (u + 1) / 2
    rho_pts = rho_max * (u + 1.0) / 2.0
    rho_w = wu * rho_max / 2.0

    v, wv = roots_legendre(n_z)
    # z in [-z_max, z_max]: z = z_max * v
    z_pts = z_max * v
    z_w = wv * z_max

    return rho_pts, rho_w, z_pts, z_w


def _quadrature_grid_split(
    n_rho_inner: int, n_rho_outer: int,
    n_z_inner: int, n_z_outer: int,
    rho_inner: float, rho_max: float,
    z_inner: float, z_max: float,
    z_center_inner: float = 0.0,
):
    """Composite quadrature: dense Gauss-Legendre near origin, sparser outside.

    Returns (rho_pts, rho_w, z_pts, z_w) where the grid covers
    rho in [0, rho_max] (denser in [0, rho_inner]) and
    z in [z_center_inner - z_max, z_center_inner + z_max] (denser in
    [z_center_inner - z_inner, z_center_inner + z_inner]).
    """
    # rho in [0, rho_inner]
    u1, w1 = roots_legendre(n_rho_inner)
    rho_pts_inner = rho_inner * (u1 + 1.0) / 2.0
    rho_w_inner = w1 * rho_inner / 2.0
    # rho in [rho_inner, rho_max]
    u2, w2 = roots_legendre(n_rho_outer)
    rho_pts_outer = rho_inner + (rho_max - rho_inner) * (u2 + 1.0) / 2.0
    rho_w_outer = w2 * (rho_max - rho_inner) / 2.0

    rho_pts = np.concatenate([rho_pts_inner, rho_pts_outer])
    rho_w = np.concatenate([rho_w_inner, rho_w_outer])

    # z grid: combine an inner dense region [z_low, z_high] = [zc - z_inner, zc + z_inner]
    # with outer regions [-z_max + zc, zc - z_inner] and [zc + z_inner, zc + z_max]
    zc = z_center_inner
    # inner region symmetric around z_center
    v_in, wz_in = roots_legendre(n_z_inner)
    z_pts_inner = zc + z_inner * v_in
    z_w_inner = wz_in * z_inner

    # outer left: [zc - z_max, zc - z_inner]
    v_ol, wz_ol = roots_legendre(n_z_outer)
    z_pts_outer_left = (zc - z_max) + (z_max - z_inner) * (v_ol + 1.0) / 2.0
    z_w_outer_left = wz_ol * (z_max - z_inner) / 2.0

    # outer right: [zc + z_inner, zc + z_max]
    v_or, wz_or = roots_legendre(n_z_outer)
    z_pts_outer_right = (zc + z_inner) + (z_max - z_inner) * (v_or + 1.0) / 2.0
    z_w_outer_right = wz_or * (z_max - z_inner) / 2.0

    z_pts = np.concatenate([z_pts_outer_left, z_pts_inner, z_pts_outer_right])
    z_w = np.concatenate([z_w_outer_left, z_w_inner, z_w_outer_right])

    return rho_pts, rho_w, z_pts, z_w


def _evaluate_radial_at_distance(orbital_callable,
                                  rho_pts: np.ndarray, z_pts: np.ndarray,
                                  z_center: float) -> np.ndarray:
    """Evaluate R(r) on (rho, z) grid where r = sqrt(rho^2 + (z - z_center)^2)."""
    # rho_pts shape (n_rho,), z_pts shape (n_z,)
    rho_grid, z_grid = np.meshgrid(rho_pts, z_pts, indexing='ij')
    r_grid = np.sqrt(rho_grid ** 2 + (z_grid - z_center) ** 2)
    # Flatten for orbital evaluator (which expects 1D), then reshape
    R_flat = orbital_callable(r_grid.flatten())
    return R_flat.reshape(r_grid.shape)


def cross_block_overlap_3d(
    orbital_A,  # callable
    z_A: float,
    orbital_B,  # callable
    z_B: float,
    n_rho: int = 60,
    n_z: int = 80,
    rho_max: float = 15.0,
    z_max: float = 15.0,
) -> float:
    """Compute <psi_A | psi_B> via 2D axial quadrature for s-orbitals.

    Volume element: dV = rho drho dz dphi.
    phi-integration: (1/sqrt(4*pi))^2 * 2*pi = 1/2 prefactor.
    Net: (1/2) * integral rho drho dz * R_A(r_A) * R_B(r_B)
    """
    rho_pts, rho_w, z_pts, z_w = _quadrature_grid(n_rho, n_z, rho_max, z_max)
    R_A = _evaluate_radial_at_distance(orbital_A, rho_pts, z_pts, z_A)
    R_B = _evaluate_radial_at_distance(orbital_B, rho_pts, z_pts, z_B)
    integrand = R_A * R_B
    # Tensor product weights: w_rho[i] * rho_pts[i] * w_z[j]
    weight_2d = np.outer(rho_w * rho_pts, z_w)
    return 0.5 * float(np.sum(weight_2d * integrand))


def cross_block_vne_3d(
    orbital_A, z_A: float,
    orbital_B, z_B: float,
    z_C: float,  # nucleus C location (on z-axis by axial symmetry)
    Z_C: float,  # charge of nucleus C
    n_rho: int = 60,
    n_z: int = 80,
    rho_max: float = 15.0,
    z_max: float = 15.0,
    eps: float = 1e-12,
) -> float:
    """Compute <psi_A | -Z_C/|r - R_C z_hat| | psi_B> via 2D axial quadrature."""
    rho_pts, rho_w, z_pts, z_w = _quadrature_grid(n_rho, n_z, rho_max, z_max)
    R_A = _evaluate_radial_at_distance(orbital_A, rho_pts, z_pts, z_A)
    R_B = _evaluate_radial_at_distance(orbital_B, rho_pts, z_pts, z_B)
    # r_C distance
    rho_grid, z_grid = np.meshgrid(rho_pts, z_pts, indexing='ij')
    r_C = np.sqrt(rho_grid ** 2 + (z_grid - z_C) ** 2 + eps ** 2)
    operator = -Z_C / r_C
    integrand = R_A * R_B * operator
    weight_2d = np.outer(rho_w * rho_pts, z_w)
    return 0.5 * float(np.sum(weight_2d * integrand))


def cross_block_inv_r_3d(
    orbital_A, z_A: float,
    orbital_B, z_B: float,
    z_C: float,
    n_rho: int = 60,
    n_z: int = 80,
    rho_max: float = 15.0,
    z_max: float = 15.0,
    eps: float = 1e-12,
) -> float:
    """Compute <psi_A | 1/|r - R_C z_hat| | psi_B> via 2D axial quadrature."""
    rho_pts, rho_w, z_pts, z_w = _quadrature_grid(n_rho, n_z, rho_max, z_max)
    R_A = _evaluate_radial_at_distance(orbital_A, rho_pts, z_pts, z_A)
    R_B = _evaluate_radial_at_distance(orbital_B, rho_pts, z_pts, z_B)
    rho_grid, z_grid = np.meshgrid(rho_pts, z_pts, indexing='ij')
    r_C = np.sqrt(rho_grid ** 2 + (z_grid - z_C) ** 2 + eps ** 2)
    integrand = R_A * R_B * (1.0 / r_C)
    weight_2d = np.outer(rho_w * rho_pts, z_w)
    return 0.5 * float(np.sum(weight_2d * integrand))


def cross_block_kinetic_via_eigenvalue(
    orbital_A_callable, z_A: float,
    orbital_B_callable, z_B: float,
    Z_B_hydrogenic: float,
    n_B: int,
    n_rho: int = 60,
    n_z: int = 80,
    rho_max: float = 15.0,
    z_max: float = 15.0,
) -> float:
    """Compute <psi_A | T | psi_B> using T|psi_B> = (E_B + Z_B/r_B) |psi_B>.

    For psi_B a HYDROGENIC orbital on center B with charge Z_B:
        H_B |psi_B> = E_B |psi_B>  where  H_B = -1/2 nabla^2 - Z_B/r_B
        T |psi_B> = (E_B + Z_B/r_B) |psi_B>
        E_B = -Z_B^2 / (2 n_B^2)

    Therefore:
        <psi_A | T | psi_B> = E_B * <psi_A | psi_B> + Z_B * <psi_A | 1/r_B | psi_B>
    """
    E_B = -Z_B_hydrogenic ** 2 / (2.0 * n_B ** 2)
    S_AB = cross_block_overlap_3d(orbital_A_callable, z_A, orbital_B_callable, z_B,
                                   n_rho, n_z, rho_max, z_max)
    inv_rB = cross_block_inv_r_3d(orbital_A_callable, z_A, orbital_B_callable, z_B,
                                    z_B, n_rho, n_z, rho_max, z_max)
    return E_B * S_AB + Z_B_hydrogenic * inv_rB


# ---------------------------------------------------------------------------
# 1D sanity check: same-center overlap and inv_r for grid resolution validation.
# ---------------------------------------------------------------------------

def sanity_check_h1s_normalization(n_rho: int = 60, n_z: int = 80,
                                    rho_max: float = 15.0, z_max: float = 15.0):
    """Check that <H 1s | H 1s> = 1 (same-center) via 2D quadrature."""
    S = cross_block_overlap_3d(h_1s, 0.0, h_1s, 0.0, n_rho, n_z, rho_max, z_max)
    return S


# ---------------------------------------------------------------------------
# Diagnostic driver
# ---------------------------------------------------------------------------

def step1_diagnostic(
    R_AB: float = 3.566,    # NaH R_eq
    n_rho: int = 80,
    n_z: int = 100,
    rho_max: float = 20.0,
    z_max: float = 20.0,
    use_multi_zeta_for_na3s: bool = False,
) -> Dict:
    """Run the Step 1 cross-block h1 diagnostic at NaH R = R_eq.

    Coordinate frame: Na at z = 0, H at z = R_AB. Orbitals are s-only for
    this diagnostic.

    Steps:
      1. Compute <Na 3s | V_ne(H) | H 1s>
      2. Compute <Na 3s | V_ne(Na) | H 1s>
      3. Compute <Na 3s | T | H 1s>
      4. Sum: h1_cross[Na 3s, H 1s] = (1) + (2) + (3)
      5. Get diagonal h1 elements (on-site eigenvalues)
      6. Diagonalize 2x2 h1 -> bonding / antibonding eigenvalues
      7. Verdict per decision gate
    """
    t0 = time.time()

    # Orbital callables
    if use_multi_zeta_for_na3s:
        na_3s_mz = na_3s_multizeta()
        def na_3s_call(r): return na_3s_mz.evaluate(r)
        na_3s_n = 3  # for kinetic via eigenvalue we still need to use a
                     # reasonable proxy; mz isn't hydrogenic so kinetic
                     # is harder. We use the hydrogenic kinetic as an
                     # approximation for the mz case in this diagnostic.
    else:
        def na_3s_call(r): return na_3s_hydrogenic(r, Z_orb=1.0)
        na_3s_n = 3

    def h_1s_call(r): return h_1s(r, Z_orb=1.0)
    h_1s_n = 1

    # Sanity check: H 1s normalization
    S_HH = cross_block_overlap_3d(h_1s_call, R_AB, h_1s_call, R_AB,
                                    n_rho, n_z, rho_max, z_max)
    S_NaNa = cross_block_overlap_3d(na_3s_call, 0.0, na_3s_call, 0.0,
                                      n_rho, n_z, rho_max, z_max)

    # 1. Cross-block overlap <Na 3s | H 1s>
    S_AB = cross_block_overlap_3d(na_3s_call, 0.0, h_1s_call, R_AB,
                                    n_rho, n_z, rho_max, z_max)

    # 2. <Na 3s | V_ne(H) | H 1s> -- nucleus H is at z = R_AB, Z = 1
    Vne_H_AB = cross_block_vne_3d(na_3s_call, 0.0, h_1s_call, R_AB,
                                    z_C=R_AB, Z_C=1.0,
                                    n_rho=n_rho, n_z=n_z, rho_max=rho_max, z_max=z_max)

    # 3. <Na 3s | V_ne(Na) | H 1s> -- nucleus Na is at z = 0, Z = 11
    Vne_Na_AB = cross_block_vne_3d(na_3s_call, 0.0, h_1s_call, R_AB,
                                     z_C=0.0, Z_C=11.0,
                                     n_rho=n_rho, n_z=n_z, rho_max=rho_max, z_max=z_max)

    # 4. Kinetic <Na 3s | T | H 1s>
    # For both orbitals hydrogenic at Z_orb=1:
    # T|H 1s> = (E_H + Z_H/r_H) |H 1s>, where E_H = -1/2, Z_H = 1
    T_AB = cross_block_kinetic_via_eigenvalue(
        na_3s_call, 0.0, h_1s_call, R_AB,
        Z_B_hydrogenic=1.0, n_B=h_1s_n,
        n_rho=n_rho, n_z=n_z, rho_max=rho_max, z_max=z_max,
    )

    # 5. Diagonal elements via the SAME 3D quadrature (for internal
    # consistency of the 2x2 diagonalization).
    # h1[Na 3s, Na 3s] = E_Na 3s on-site (hydrogenic -1/18 for Z_orb=1, n=3)
    #                    + integral <Na 3s | V_ne(H) | Na 3s>
    # h1[H 1s, H 1s]   = E_H 1s on-site (-1/2)
    #                    + integral <H 1s | V_ne(Na) | H 1s>
    Vne_H_AA = cross_block_vne_3d(na_3s_call, 0.0, na_3s_call, 0.0,
                                     z_C=R_AB, Z_C=1.0,
                                     n_rho=n_rho, n_z=n_z, rho_max=rho_max, z_max=z_max)
    Vne_Na_BB = cross_block_vne_3d(h_1s_call, R_AB, h_1s_call, R_AB,
                                      z_C=0.0, Z_C=11.0,
                                      n_rho=n_rho, n_z=n_z, rho_max=rho_max, z_max=z_max)
    # On-site eigenvalues
    if use_multi_zeta_for_na3s:
        E_Na_onsite = -0.170  # FrozenCore-screened Na 3s eigenvalue (F1 memo)
    else:
        E_Na_onsite = -1.0 / (2.0 * na_3s_n ** 2)  # = -1/18 for Z_orb=1, n=3
    E_H_onsite = -0.5  # hydrogenic 1s

    h11 = E_Na_onsite + Vne_H_AA  # h1[Na 3s, Na 3s]
    h22 = E_H_onsite + Vne_Na_BB  # h1[H 1s, H 1s]
    h12 = T_AB + Vne_H_AB + Vne_Na_AB  # h1[Na 3s, H 1s] (cross-block)

    # Diagonalize 2x2 h matrix (treating basis as orthogonal; S_AB
    # in principle requires generalized eigenvalue, but we report both)
    H_matrix = np.array([[h11, h12], [h12, h22]])
    S_matrix = np.array([[1.0, S_AB], [S_AB, 1.0]])

    # Orthogonal-basis diagonalization
    eigs_ortho, vecs_ortho = np.linalg.eigh(H_matrix)
    E_bonding_ortho = eigs_ortho[0]
    E_antibonding_ortho = eigs_ortho[1]
    splitting_ortho = E_bonding_ortho - E_antibonding_ortho  # negative if bonding favored

    # Generalized eigenvalue (with overlap)
    from scipy.linalg import eigh as scipy_eigh
    eigs_gen, vecs_gen = scipy_eigh(H_matrix, S_matrix)
    E_bonding_gen = eigs_gen[0]
    E_antibonding_gen = eigs_gen[1]
    splitting_gen = E_bonding_gen - E_antibonding_gen

    # Determine which eigenvector is "bonding" (same sign on both centers)
    # vs antibonding (opposite signs)
    def classify(v):
        if v[0] * v[1] > 0:
            return "bonding"
        else:
            return "antibonding"

    label_low_ortho = classify(vecs_ortho[:, 0])
    label_low_gen = classify(vecs_gen[:, 0])

    # Verdict per gate
    if splitting_ortho < -0.05:
        verdict = "PROCEED_TO_STEP_2_STRONG_BONDING"
    elif abs(splitting_ortho) < 0.05:
        verdict = "PROCEED_TO_STEP_2_WEAK_SPLITTING"
    else:
        # bonding NOT favored (splitting > 0.05, antibonding lower)
        verdict = "STOP_BONDING_NOT_FAVORED"

    result = {
        'R_AB_bohr': R_AB,
        'na_3s_basis': 'multi_zeta' if use_multi_zeta_for_na3s else 'hydrogenic',
        'quadrature': {
            'n_rho': n_rho, 'n_z': n_z,
            'rho_max': rho_max, 'z_max': z_max,
        },
        'sanity_check': {
            'S_NaNa_should_be_1': S_NaNa,
            'S_HH_should_be_1': S_HH,
        },
        'matrix_elements': {
            'overlap_S_AB': S_AB,
            'Vne_H_NaH_offdiag_Ha': Vne_H_AB,
            'Vne_Na_NaH_offdiag_Ha': Vne_Na_AB,
            'kinetic_T_NaH_offdiag_Ha': T_AB,
            'cross_block_h1_NaH_total_Ha': h12,
            'Vne_H_Na3s_diag_Ha': Vne_H_AA,
            'Vne_Na_H1s_diag_Ha': Vne_Na_BB,
        },
        'on_site_eigenvalues': {
            'E_Na_3s_onsite_Ha': E_Na_onsite,
            'E_H_1s_onsite_Ha': E_H_onsite,
        },
        'h1_2x2_matrix': {
            'h11_Na3s_diag_Ha': h11,
            'h22_H1s_diag_Ha': h22,
            'h12_offdiag_Ha': h12,
        },
        'eigenvalues_orthogonal_basis': {
            'E_lower_Ha': E_bonding_ortho,
            'E_upper_Ha': E_antibonding_ortho,
            'splitting_Ha': splitting_ortho,
            'lower_label': label_low_ortho,
        },
        'eigenvalues_generalized_with_overlap': {
            'E_lower_Ha': E_bonding_gen,
            'E_upper_Ha': E_antibonding_gen,
            'splitting_Ha': splitting_gen,
            'lower_label': label_low_gen,
        },
        'verdict': verdict,
        'wall_time_s': time.time() - t0,
    }
    return result


def main():
    print("=" * 78)
    print("Sprint F3 Step 1 — Cross-block h1 algebraic diagnostic")
    print("=" * 78)
    print()
    print("System: NaH at R = R_eq = 3.566 bohr")
    print("Coordinate: Na at z = 0, H at z = R_AB")
    print("Basis: hydrogenic (Z_orb=1) Na 3s and H 1s (framework baseline)")
    print()

    # Sanity check first
    S_HH = sanity_check_h1s_normalization(n_rho=80, n_z=100,
                                            rho_max=20.0, z_max=20.0)
    print(f"Sanity check: <H 1s | H 1s> = {S_HH:.6f} (should be ~1.0)")
    print()

    # Main diagnostic, hydrogenic basis
    print("---- Main diagnostic: hydrogenic Na 3s ----")
    result_hyd = step1_diagnostic(
        R_AB=3.566, n_rho=80, n_z=100, rho_max=20.0, z_max=20.0,
        use_multi_zeta_for_na3s=False,
    )

    print(f"  Cross-block overlap <Na 3s | H 1s>: {result_hyd['matrix_elements']['overlap_S_AB']:+.6f}")
    print(f"  <Na 3s | V_ne(H) | H 1s>:           {result_hyd['matrix_elements']['Vne_H_NaH_offdiag_Ha']:+.6f} Ha")
    print(f"  <Na 3s | V_ne(Na) | H 1s>:          {result_hyd['matrix_elements']['Vne_Na_NaH_offdiag_Ha']:+.6f} Ha")
    print(f"  <Na 3s | T | H 1s>:                 {result_hyd['matrix_elements']['kinetic_T_NaH_offdiag_Ha']:+.6f} Ha")
    print(f"  h1[Na 3s, H 1s] TOTAL:              {result_hyd['matrix_elements']['cross_block_h1_NaH_total_Ha']:+.6f} Ha")
    print()
    print(f"  h11 (Na 3s diag): {result_hyd['h1_2x2_matrix']['h11_Na3s_diag_Ha']:+.6f} Ha")
    print(f"  h22 (H 1s diag):  {result_hyd['h1_2x2_matrix']['h22_H1s_diag_Ha']:+.6f} Ha")
    print(f"  h12 (cross):      {result_hyd['h1_2x2_matrix']['h12_offdiag_Ha']:+.6f} Ha")
    print()
    print("  ORTHOGONAL-basis 2x2 diagonalization:")
    print(f"    E_lower: {result_hyd['eigenvalues_orthogonal_basis']['E_lower_Ha']:+.6f} Ha ({result_hyd['eigenvalues_orthogonal_basis']['lower_label']})")
    print(f"    E_upper: {result_hyd['eigenvalues_orthogonal_basis']['E_upper_Ha']:+.6f} Ha")
    print(f"    Splitting (lower - upper): {result_hyd['eigenvalues_orthogonal_basis']['splitting_Ha']:+.6f} Ha")
    print()
    print("  GENERALIZED-basis (with overlap S):")
    print(f"    E_lower: {result_hyd['eigenvalues_generalized_with_overlap']['E_lower_Ha']:+.6f} Ha ({result_hyd['eigenvalues_generalized_with_overlap']['lower_label']})")
    print(f"    E_upper: {result_hyd['eigenvalues_generalized_with_overlap']['E_upper_Ha']:+.6f} Ha")
    print(f"    Splitting: {result_hyd['eigenvalues_generalized_with_overlap']['splitting_Ha']:+.6f} Ha")
    print()
    print(f"  VERDICT (hydrogenic): {result_hyd['verdict']}")
    print()

    # Multi-zeta diagnostic
    print("---- Main diagnostic: multi-zeta Na 3s ----")
    result_mz = step1_diagnostic(
        R_AB=3.566, n_rho=80, n_z=100, rho_max=20.0, z_max=20.0,
        use_multi_zeta_for_na3s=True,
    )

    print(f"  Cross-block overlap <Na 3s_mz | H 1s>: {result_mz['matrix_elements']['overlap_S_AB']:+.6f}")
    print(f"  <Na 3s_mz | V_ne(H) | H 1s>:           {result_mz['matrix_elements']['Vne_H_NaH_offdiag_Ha']:+.6f} Ha")
    print(f"  <Na 3s_mz | V_ne(Na) | H 1s>:          {result_mz['matrix_elements']['Vne_Na_NaH_offdiag_Ha']:+.6f} Ha")
    print(f"  <Na 3s_mz | T | H 1s>:                 {result_mz['matrix_elements']['kinetic_T_NaH_offdiag_Ha']:+.6f} Ha")
    print(f"  h1[Na 3s_mz, H 1s] TOTAL:              {result_mz['matrix_elements']['cross_block_h1_NaH_total_Ha']:+.6f} Ha")
    print()
    print(f"  h11 (Na 3s_mz diag): {result_mz['h1_2x2_matrix']['h11_Na3s_diag_Ha']:+.6f} Ha")
    print(f"  h22 (H 1s diag):     {result_mz['h1_2x2_matrix']['h22_H1s_diag_Ha']:+.6f} Ha")
    print(f"  h12 (cross):         {result_mz['h1_2x2_matrix']['h12_offdiag_Ha']:+.6f} Ha")
    print()
    print("  ORTHOGONAL-basis 2x2 diagonalization:")
    print(f"    E_lower: {result_mz['eigenvalues_orthogonal_basis']['E_lower_Ha']:+.6f} Ha ({result_mz['eigenvalues_orthogonal_basis']['lower_label']})")
    print(f"    E_upper: {result_mz['eigenvalues_orthogonal_basis']['E_upper_Ha']:+.6f} Ha")
    print(f"    Splitting: {result_mz['eigenvalues_orthogonal_basis']['splitting_Ha']:+.6f} Ha")
    print()
    print(f"  VERDICT (multi-zeta): {result_mz['verdict']}")
    print()

    out_path = "debug/data/sprint_f3_step1_diagnostic.json"
    out = {
        'sprint': 'F3 Step 1 — cross-block h1 algebraic diagnostic',
        'date': '2026-05-23',
        'hydrogenic': result_hyd,
        'multi_zeta': result_mz,
        'sanity_S_HH': S_HH,
    }
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
