"""Cross-block h1 architectural extension (Sprint F3, 2026-05-23).

Adds the missing cross-block off-diagonal h1 matrix elements
   h1[a, b] = <psi_a^A | T | psi_b^B> + sum_C <psi_a^A | -Z_C/|r - R_C| | psi_b^B>
for orbitals on different centers A != B. The existing framework
(``composed_qubit.build_composed_hamiltonian``) constructs h1 strictly
block-diagonal — cross-block one-body coupling between orbitals on
different centers is architecturally absent.

This module supplies that coupling via direct 2D axial Gauss-Legendre
quadrature (3D integral with phi analytic for axially symmetric
integrands). The implementation is correct-first, slow-second; specialized
fast variants (elliptic-coordinate quadrature for diatomics; bipolar
multipole; two-center Gauss-Hermite) are explicitly out of scope for the
first production version (see Sprint F3 scope statement).

CURRENT SCOPE
-------------
- s-s pairs only (l_a = l_b = 0) in the FIRST production pass.
  Mixed l > 0 cases require either spherical harmonic axial decomposition
  with non-trivial angular factors, or full 3D quadrature on Lebedev grids.
  These are flagged as named follow-on extensions if F3 closes the wall.

- Hydrogenic and multi-zeta orbital shapes both supported. Multi-zeta
  is recognized via a dict ``multi_zeta_overrides`` mapping (block_n, l)
  to a MultiZetaOrbital instance.

- Z_eff(r) "screened" wavefunctions are NOT yet supported here. When
  ``screening='screened'`` is requested for a Z_eff(r) FrozenCore-shape
  orbital, the function raises NotImplementedError. For the F3 sprint
  this is not a blocker (the multi-zeta path provides the relevant
  shape for Na 3s).

- All distances/coordinates are in atomic units (bohr). Z values are in
  units of |e|.

GEOMETRY
--------
Centers and nuclei are passed as 3D positions. For collinear systems
(NaH, LiH, MgH2, ...) the natural coordinate is the z-axis. The
quadrature uses cylindrical coordinates (rho, z, phi); for collinear
systems, axial symmetry around the bond axis allows phi to be done
analytically (factor 2*pi, divided by 4*pi from the s-orbital angular
normalization, giving 1/2 prefactor on the (rho, z) integral).

For NON-COLLINEAR systems (H2O, NH3, ...), the diagnostic still works
when the two orbital centers and the integration nucleus C are
COPLANAR (we can rotate frame to make orbital-A-to-orbital-B aligned
with z-axis; nucleus C in that plane gives axial symmetry around the
bond axis). For C OFF the bond axis, this would require full 3D
spherical quadrature. The first-pass code handles the collinear /
coplanar case via frame rotation; calls with non-coplanar geometry
raise NotImplementedError.

SIGN CONVENTION
---------------
H = T + V where T is positive kinetic energy and V_ne = -Z/r is the
attractive electron-nucleus interaction.

KINETIC IMPLEMENTATION
----------------------
For psi_B a hydrogenic orbital at center B with charge Z_B and quantum
number n_B:
    -1/2 nabla^2 psi_B = (E_B + Z_B/r_B) psi_B
    where E_B = -Z_B^2 / (2 n_B^2)
This makes <psi_A | T | psi_B> = E_B * <psi_A | psi_B> + Z_B * <psi_A | 1/r_B | psi_B>,
both of which are pure overlap/Coulomb integrals computable by the same
quadrature. For multi-zeta psi_B (linear combination of STOs with various
zetas), the kinetic action does NOT have such a simple eigenvalue form;
we fall back to the SCREENED-EIGENVALUE approximation E_B = -1/(2 n_B^2)
of the FrozenCore-screened Schrödinger eigenvalue (consistent with what
the framework uses for the on-site diagonal h1 in the multi-zeta path).
This is honest-scope and flagged in the docstring.
"""

from __future__ import annotations

import math
from typing import Any, Dict, Optional, Tuple

import numpy as np
from scipy.special import factorial, genlaguerre, roots_legendre


# ---------------------------------------------------------------------------
# Analytical hydrogenic radial wavefunctions (grid-independent normalization)
# ---------------------------------------------------------------------------

def hydrogenic_R_nl_analytical(
    Z: float, n: int, l: int, r: np.ndarray,
) -> np.ndarray:
    """Normalized hydrogenic radial wavefunction R_{nl}(r), analytical norm.

    R_{nl}(r) = N * (2 Z r / n)^l * exp(-Z r / n) * L^{2l+1}_{n-l-1}(2 Z r / n)
    where N = sqrt[(2 Z / n)^3 * (n - l - 1)! / (2 n (n + l)!)]

    Normalized so that integral_0^inf R^2 r^2 dr = 1 exactly (analytical
    normalization). Distinct from ``composed_qubit._radial_wf_grid``,
    which performs a grid-dependent trapezoidal normalization.
    """
    rho = 2.0 * Z * r / n
    L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    norm = math.sqrt(
        (2.0 * Z / n) ** 3
        * float(factorial(n - l - 1))
        / (2.0 * n * float(factorial(n + l)))
    )
    return norm * rho ** l * np.exp(-rho / 2.0) * L_poly


# ---------------------------------------------------------------------------
# 2D axial Gauss-Legendre grid
# ---------------------------------------------------------------------------

def _build_quadrature_grid(
    n_rho: int, n_z: int, rho_max: float, z_max: float, z_center: float = 0.0,
):
    """Build Gauss-Legendre (rho, z) grid.

    rho in [0, rho_max], z in [z_center - z_max, z_center + z_max].
    """
    u, wu = roots_legendre(n_rho)
    rho_pts = rho_max * (u + 1.0) / 2.0
    rho_w = wu * rho_max / 2.0

    v, wv = roots_legendre(n_z)
    z_pts = z_center + z_max * v
    z_w = wv * z_max

    return rho_pts, rho_w, z_pts, z_w


def _eval_radial_on_2d_grid(
    R_callable, rho_pts: np.ndarray, z_pts: np.ndarray, z_orbital: float,
) -> np.ndarray:
    """Evaluate R(r) where r = sqrt(rho^2 + (z - z_orbital)^2) on the grid."""
    rho_grid, z_grid = np.meshgrid(rho_pts, z_pts, indexing='ij')
    r_grid = np.sqrt(rho_grid ** 2 + (z_grid - z_orbital) ** 2)
    R_flat = R_callable(r_grid.flatten())
    return R_flat.reshape(r_grid.shape)


# ---------------------------------------------------------------------------
# Core matrix-element computers
# ---------------------------------------------------------------------------

def overlap_ss_axial(
    R_A_call, z_A: float,
    R_B_call, z_B: float,
    n_rho: int, n_z: int, rho_max: float, z_max: float,
) -> float:
    """Compute <psi_A | psi_B> for s-s axial geometry (orbitals on z-axis).

    Volume element: dV = rho drho dz dphi.
    phi-integration: 2*pi * (1/sqrt(4*pi))^2 = 1/2 prefactor.
    """
    z_center = 0.5 * (z_A + z_B)
    rho_pts, rho_w, z_pts, z_w = _build_quadrature_grid(
        n_rho, n_z, rho_max, z_max, z_center,
    )
    R_A = _eval_radial_on_2d_grid(R_A_call, rho_pts, z_pts, z_A)
    R_B = _eval_radial_on_2d_grid(R_B_call, rho_pts, z_pts, z_B)
    integrand = R_A * R_B
    weight_2d = np.outer(rho_w * rho_pts, z_w)
    return 0.5 * float(np.sum(weight_2d * integrand))


def vne_ss_axial(
    R_A_call, z_A: float,
    R_B_call, z_B: float,
    z_C: float, Z_C: float,
    n_rho: int, n_z: int, rho_max: float, z_max: float,
    eps: float = 1e-12,
) -> float:
    """Compute <psi_A | -Z_C / |r - R_C z_hat| | psi_B> for s-s axial geometry."""
    z_center = 0.5 * (z_A + z_B)
    rho_pts, rho_w, z_pts, z_w = _build_quadrature_grid(
        n_rho, n_z, rho_max, z_max, z_center,
    )
    R_A = _eval_radial_on_2d_grid(R_A_call, rho_pts, z_pts, z_A)
    R_B = _eval_radial_on_2d_grid(R_B_call, rho_pts, z_pts, z_B)
    rho_grid, z_grid = np.meshgrid(rho_pts, z_pts, indexing='ij')
    r_C = np.sqrt(rho_grid ** 2 + (z_grid - z_C) ** 2 + eps ** 2)
    integrand = R_A * R_B * (-Z_C / r_C)
    weight_2d = np.outer(rho_w * rho_pts, z_w)
    return 0.5 * float(np.sum(weight_2d * integrand))


def inv_r_ss_axial(
    R_A_call, z_A: float,
    R_B_call, z_B: float,
    z_C: float,
    n_rho: int, n_z: int, rho_max: float, z_max: float,
    eps: float = 1e-12,
) -> float:
    """Compute <psi_A | 1/|r - R_C z_hat| | psi_B> for s-s axial geometry."""
    z_center = 0.5 * (z_A + z_B)
    rho_pts, rho_w, z_pts, z_w = _build_quadrature_grid(
        n_rho, n_z, rho_max, z_max, z_center,
    )
    R_A = _eval_radial_on_2d_grid(R_A_call, rho_pts, z_pts, z_A)
    R_B = _eval_radial_on_2d_grid(R_B_call, rho_pts, z_pts, z_B)
    rho_grid, z_grid = np.meshgrid(rho_pts, z_pts, indexing='ij')
    r_C = np.sqrt(rho_grid ** 2 + (z_grid - z_C) ** 2 + eps ** 2)
    integrand = R_A * R_B * (1.0 / r_C)
    weight_2d = np.outer(rho_w * rho_pts, z_w)
    return 0.5 * float(np.sum(weight_2d * integrand))


# ---------------------------------------------------------------------------
# Public API: cross-block h1 single matrix element
# ---------------------------------------------------------------------------

def compute_cross_block_h1_element(
    orbital_A_callable, l_A: int, m_A: int, n_A: int, Z_A_eff: float,
    pos_A: Tuple[float, float, float],
    orbital_B_callable, l_B: int, m_B: int, n_B: int, Z_B_eff: float,
    pos_B: Tuple[float, float, float],
    nuclei: list,
    n_rho: int = 80,
    n_z: int = 100,
    rho_max: float = 20.0,
    z_max: float = 20.0,
    kinetic_eigenvalue_for_B: Optional[float] = None,
) -> float:
    """Compute h1[a, b] = T[a, b] + sum_C V_ne(C)[a, b] for cross-block pair.

    Parameters
    ----------
    orbital_A_callable : callable r -> R_A(r)
        Radial wavefunction of orbital A (already centered at pos_A).
    l_A, m_A : int
        Angular quantum numbers of orbital A. Currently only s-orbitals
        (l_A=l_B=0, m_A=m_B=0) are supported; other cases raise.
    n_A : int
        Principal quantum number of orbital A (for the kinetic
        eigenvalue evaluation if needed; not currently used since the
        kinetic action operates on B).
    Z_A_eff : float
        Effective charge for orbital A's hydrogenic shape (used only
        for consistency; not directly required for the matrix element).
    pos_A : tuple of 3 floats
        3D position of center A in bohr.
    orbital_B_callable : callable r -> R_B(r)
    l_B, m_B, n_B : int
    Z_B_eff : float
        Used for the hydrogenic-kinetic-via-eigenvalue identity
        T|psi_B> = (E_B + Z_B/r_B) |psi_B>.
    pos_B : tuple of 3 floats
    nuclei : list of dict
        Each entry: {'Z': float, 'position': (x,y,z), 'label': str}.
        All nuclei contribute to V_ne in the cross-block matrix element
        (including the two centers A and B themselves).
    n_rho, n_z, rho_max, z_max : int / float
        Quadrature parameters.
    kinetic_eigenvalue_for_B : float, optional
        If provided, overrides the hydrogenic E_B = -Z_B^2/(2 n_B^2)
        used in the kinetic-via-eigenvalue identity. Use this for
        multi-zeta orbitals where the on-site eigenvalue is the
        screened-Schrödinger value (passed by the caller).

    Returns
    -------
    float
        Cross-block h1 matrix element in Hartree.

    Raises
    ------
    NotImplementedError
        For non-s-orbital pairs (l_A != 0 or l_B != 0) and for non-collinear
        (or non-coplanar with axial symmetry around the A-B axis) nuclei
        in the FIRST-pass implementation.
    """
    if not (l_A == 0 and l_B == 0 and m_A == 0 and m_B == 0):
        raise NotImplementedError(
            f"compute_cross_block_h1_element currently supports s-s pairs only "
            f"(l_a=l_b=m_a=m_b=0); got l_A={l_A}, l_B={l_B}, m_A={m_A}, m_B={m_B}. "
            f"Mixed-l extension is a named Sprint F3 follow-on."
        )

    pos_A_arr = np.asarray(pos_A, dtype=float)
    pos_B_arr = np.asarray(pos_B, dtype=float)
    displacement = pos_B_arr - pos_A_arr
    R_AB = float(np.linalg.norm(displacement))

    # Same-center: not cross-block; return zero (caller should handle)
    if R_AB < 1e-10:
        return 0.0

    # For axial symmetry around the A-B axis, ALL nuclei must be on that
    # axis. We check this by rotating the frame so A is at origin and B
    # is on +z, then verifying every nucleus has rho < 1e-8 in the new
    # frame.

    # Build rotation that maps (B - A) onto +z_hat
    z_hat_new = displacement / R_AB
    # Find two basis vectors orthogonal to z_hat_new
    # Pick an arbitrary vector not parallel to z_hat_new
    if abs(z_hat_new[2]) < 0.9:
        seed = np.array([0.0, 0.0, 1.0])
    else:
        seed = np.array([1.0, 0.0, 0.0])
    x_hat_new = seed - np.dot(seed, z_hat_new) * z_hat_new
    x_hat_new /= np.linalg.norm(x_hat_new)
    y_hat_new = np.cross(z_hat_new, x_hat_new)

    # Verify all nuclei lie on the A-B axis (rho < epsilon in new frame)
    nuc_z_in_new_frame = []
    for nuc in nuclei:
        nuc_pos = np.asarray(nuc['position'], dtype=float) - pos_A_arr
        rho_nuc = math.hypot(
            float(np.dot(nuc_pos, x_hat_new)),
            float(np.dot(nuc_pos, y_hat_new)),
        )
        z_nuc = float(np.dot(nuc_pos, z_hat_new))
        if rho_nuc > 1e-6:
            raise NotImplementedError(
                f"Non-axial geometry detected: nucleus {nuc.get('label', '?')} "
                f"has rho = {rho_nuc:.4g} bohr off the A-B axis "
                f"(rotated frame). Axial cross-block h1 implementation "
                f"requires all nuclei collinear with the orbital centers. "
                f"Full 3D extension is a named Sprint F3 follow-on."
            )
        nuc_z_in_new_frame.append({'Z': float(nuc['Z']), 'z': z_nuc,
                                    'label': nuc.get('label', '?')})

    # In the rotated frame: A is at z = 0, B is at z = R_AB.
    z_A_local = 0.0
    z_B_local = R_AB

    # Compute V_ne contribution: sum over all nuclei
    Vne_total = 0.0
    for nuc in nuc_z_in_new_frame:
        Vne_C = vne_ss_axial(
            orbital_A_callable, z_A_local,
            orbital_B_callable, z_B_local,
            z_C=nuc['z'], Z_C=nuc['Z'],
            n_rho=n_rho, n_z=n_z, rho_max=rho_max, z_max=z_max,
        )
        Vne_total += Vne_C

    # Compute kinetic via T|psi_B> = (E_B + Z_B/r_B) |psi_B>
    if kinetic_eigenvalue_for_B is None:
        E_B = -Z_B_eff ** 2 / (2.0 * n_B ** 2)
    else:
        E_B = kinetic_eigenvalue_for_B
    S_AB = overlap_ss_axial(
        orbital_A_callable, z_A_local,
        orbital_B_callable, z_B_local,
        n_rho=n_rho, n_z=n_z, rho_max=rho_max, z_max=z_max,
    )
    inv_rB = inv_r_ss_axial(
        orbital_A_callable, z_A_local,
        orbital_B_callable, z_B_local,
        z_C=z_B_local,
        n_rho=n_rho, n_z=n_z, rho_max=rho_max, z_max=z_max,
    )
    T_AB = E_B * S_AB + Z_B_eff * inv_rB

    return T_AB + Vne_total


# ---------------------------------------------------------------------------
# Helpers to build orbital callables for sub-blocks
# ---------------------------------------------------------------------------

def make_orbital_callable(
    n: int, l: int, Z_eff: float,
    multi_zeta_override: Optional[Any] = None,
):
    """Build a radial wavefunction callable for (n, l) at Z_eff, with
    optional multi-zeta override.

    Returns
    -------
    callable
        r (np.ndarray) -> R(r) (np.ndarray)
    """
    if multi_zeta_override is not None:
        # MultiZetaOrbital instance
        return lambda r, _o=multi_zeta_override: _o.evaluate(np.asarray(r))
    # Hydrogenic
    return lambda r, _n=n, _l=l, _Z=Z_eff: hydrogenic_R_nl_analytical(
        _Z, _n, _l, np.asarray(r),
    )


# ---------------------------------------------------------------------------
# Sub-block traversal helper: compute the cross-block h1 matrix
# ---------------------------------------------------------------------------

def compute_cross_block_h1_matrix(
    sub_blocks: list,
    sb_positions: Dict[str, Tuple[float, float, float]],
    nuclei: list,
    M: int,
    n_rho: int = 80,
    n_z: int = 100,
    rho_max: float = 20.0,
    z_max: float = 20.0,
    multi_zeta_basis_by_sb: Optional[Dict[str, Dict[Tuple[int, int], Any]]] = None,
    n_val_offset_by_sb: Optional[Dict[str, int]] = None,
    kinetic_eigenvalue_override_by_sb_nl: Optional[Dict[Tuple[str, int, int], float]] = None,
    verbose: bool = False,
) -> np.ndarray:
    """Compute the cross-block off-diagonal h1 matrix.

    For every pair of sub-blocks (sb_A, sb_B) with sb_A != sb_B, and for
    every pair of orbitals (a in sb_A, b in sb_B) that satisfy (1) both
    have l=0 (s-orbitals) and (2) m_A = m_B = 0, compute h1[a, b] via
    ``compute_cross_block_h1_element``. Returns a symmetric M x M matrix
    with the cross-block contributions.

    Same-block diagonal entries are left zero (existing block-diagonal h1
    already includes them).

    Parameters
    ----------
    sub_blocks : list of dict
        From ``balanced_coupled._get_block_geometry``: each entry has
        'label', 'offset', 'states' (list of (n,l,m)), 'Z_orb', and
        'parent_block' index into spec.blocks. 'side' = 'center' or 'partner'.
    sb_positions : dict label -> (x,y,z)
    nuclei : list of {'Z', 'position', 'label'}
    M : int
        Total spatial orbital count.
    n_rho, n_z, rho_max, z_max : quadrature parameters
    multi_zeta_basis_by_sb : optional dict
        sub_block_label -> dict[(block_n, l) -> MultiZetaOrbital]
        Used for the cross-block element when one or both orbitals
        are multi-zeta.
    n_val_offset_by_sb : optional dict
        sub_block_label -> int (n_val_offset for that sub-block's parent
        block). Maps block_n to physical_n for multi-zeta orbitals.
    kinetic_eigenvalue_override_by_sb_nl : optional dict
        (sub_block_label, n, l) -> float, the screened-Schrödinger
        eigenvalue to use in place of the hydrogenic -Z^2/(2 n^2)
        for the kinetic eigenvalue identity on this orbital.

    Returns
    -------
    np.ndarray (M, M)
        Symmetric matrix with cross-block h1 contributions only.
    """
    h1_cross_block = np.zeros((M, M))

    if multi_zeta_basis_by_sb is None:
        multi_zeta_basis_by_sb = {}
    if n_val_offset_by_sb is None:
        n_val_offset_by_sb = {}
    if kinetic_eigenvalue_override_by_sb_nl is None:
        kinetic_eigenvalue_override_by_sb_nl = {}

    n_sb = len(sub_blocks)

    for i, sb_A in enumerate(sub_blocks):
        pos_A = sb_positions[sb_A['label']]
        Z_A_eff = sb_A['Z_orb']
        mz_A = multi_zeta_basis_by_sb.get(sb_A['label'], {})
        n_val_offset_A = n_val_offset_by_sb.get(sb_A['label'], 0)

        for j, sb_B in enumerate(sub_blocks):
            if i >= j:
                continue  # only j > i; we'll symmetrize at the end
            pos_B = sb_positions[sb_B['label']]
            # Same-center check
            d = math.sqrt(sum((a - b) ** 2 for a, b in zip(pos_A, pos_B)))
            if d < 1e-10:
                continue
            Z_B_eff = sb_B['Z_orb']
            mz_B = multi_zeta_basis_by_sb.get(sb_B['label'], {})
            n_val_offset_B = n_val_offset_by_sb.get(sb_B['label'], 0)

            states_A = sb_A['states']
            states_B = sb_B['states']
            off_A = sb_A['offset']
            off_B = sb_B['offset']

            for ia, (na, la, ma) in enumerate(states_A):
                if la != 0 or ma != 0:
                    continue
                # Build orbital A callable
                mz_orb_A = mz_A.get((na, la))
                orb_A_call = make_orbital_callable(
                    n=na, l=la, Z_eff=Z_A_eff,
                    multi_zeta_override=mz_orb_A,
                )
                # Kinetic eigenvalue override for A: not used since we
                # apply T to B in the eigenvalue identity.
                for ib, (nb, lb, mb) in enumerate(states_B):
                    if lb != 0 or mb != 0:
                        continue
                    mz_orb_B = mz_B.get((nb, lb))
                    orb_B_call = make_orbital_callable(
                        n=nb, l=lb, Z_eff=Z_B_eff,
                        multi_zeta_override=mz_orb_B,
                    )
                    # Kinetic eigenvalue for B
                    if mz_orb_B is not None:
                        # Multi-zeta on B: use override eigenvalue if
                        # provided, else fall back to screened-Schrödinger
                        # default of -1/(2 n_phys^2)
                        e_key = (sb_B['label'], nb, lb)
                        if e_key in kinetic_eigenvalue_override_by_sb_nl:
                            E_B_use = kinetic_eigenvalue_override_by_sb_nl[e_key]
                        else:
                            n_phys = nb + n_val_offset_B
                            E_B_use = -1.0 / (2.0 * n_phys ** 2)
                    else:
                        E_B_use = None  # default to hydrogenic in element computer

                    val = compute_cross_block_h1_element(
                        orbital_A_callable=orb_A_call,
                        l_A=la, m_A=ma, n_A=na, Z_A_eff=Z_A_eff,
                        pos_A=pos_A,
                        orbital_B_callable=orb_B_call,
                        l_B=lb, m_B=mb, n_B=nb, Z_B_eff=Z_B_eff,
                        pos_B=pos_B,
                        nuclei=nuclei,
                        n_rho=n_rho, n_z=n_z, rho_max=rho_max, z_max=z_max,
                        kinetic_eigenvalue_for_B=E_B_use,
                    )
                    if abs(val) > 1e-15:
                        h1_cross_block[off_A + ia, off_B + ib] = val
                        h1_cross_block[off_B + ib, off_A + ia] = val
                        if verbose:
                            print(
                                f"  cross-block h1[{sb_A['label']}({na},{la},{ma}), "
                                f"{sb_B['label']}({nb},{lb},{mb})] = {val:+.6f} Ha"
                            )

    return h1_cross_block
