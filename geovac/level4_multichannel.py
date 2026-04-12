"""
Level 4 multichannel solver for H2.

Extends Phase 1 (l1=l2=0 only) to arbitrary l_max by coupling (l1, l2)
partial-wave channels through the anisotropic nuclear field and the
electron-electron Gaunt integrals.

The angular eigenvalue problem at fixed (rho, R_e):

    [-1/2 d^2/dalpha^2 + V_diag(alpha) + W_coupling(alpha)] u = mu u

where channels are labeled by (l1, l2) pairs with l1+l2 even (Sigma_g).

References:
  - Level 4 design: papers/core/level4_geometry_design.md
  - Lin, Phys. Rep. 257, 1 (1995)
  - hyperspherical_angular.py for Gaunt integral machinery
"""

import numpy as np
from scipy.linalg import eigh
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal
from typing import Tuple, List, Union
from math import sqrt

from geovac.hyperspherical_angular import gaunt_integral


def _channel_list(
    l_max: int,
    homonuclear: bool = True,
) -> List[Tuple[int, int]]:
    """
    Build list of (l1, l2) channels for sigma orbitals (m=0).

    Constraints:
    - m1 = m2 = 0 (sigma orbitals)
    - l1 + l2 even (gerade symmetry) when homonuclear=True
    - All l1, l2 combinations when homonuclear=False
    - Both orderings (l1, l2) and (l2, l1) included when l1 != l2,
      because they have DIFFERENT centrifugal potentials and the
      nuclear field couples (0,0) -> (0,2) and (0,0) -> (2,0)
      through different electrons.

    Returns list of (l1, l2) tuples sorted by l1+l2 then l1.
    """
    channels = []
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            if homonuclear and (l1 + l2) % 2 != 0:
                continue
            channels.append((l1, l2))
    channels.sort(key=lambda x: (x[0] + x[1], x[0]))
    return channels


def _channel_list_extended(
    l_max: int, m_max: int,
    l_max_per_m: Union[dict, None] = None,
    homonuclear: bool = True,
) -> List[Tuple[int, int, int, int]]:
    """
    Build list of (l1, m1, l2, m2) channels with M = 0.

    Constraints:
    - m1 + m2 = 0 (total M = 0)
    - l1 + l2 even (gerade symmetry) when homonuclear=True
    - All l1+l2 parities when homonuclear=False
    - |m1| <= l1, |m2| = |m1| <= l2
    - |m1| <= m_max

    When m_max=0, returns [(l1, 0, l2, 0)] matching _channel_list order.

    Parameters
    ----------
    l_max : int
        Default maximum angular momentum per electron.
    m_max : int
        Maximum |m| per electron.
    l_max_per_m : dict or None
        Per-|m| angular momentum limit. E.g., {0: 4, 1: 3} uses l_max=4
        for sigma (m=0) and l_max=3 for pi (|m|=1). If None, uses l_max
        for all m values. This allows keeping the total channel count in a
        regime where the single-channel adiabatic approximation is valid.
    homonuclear : bool
        If True, enforce l1+l2 even (gerade symmetry). Default True.

    Returns list sorted by (l1+l2, l1, |m1|, m1).
    """
    channels: List[Tuple[int, int, int, int]] = []
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            if homonuclear and (l1 + l2) % 2 != 0:
                continue
            m_limit = min(l1, l2, m_max)
            for m1 in range(-m_limit, m_limit + 1):
                m2 = -m1
                am = abs(m1)
                # Apply per-m angular momentum limit
                if l_max_per_m is not None and am in l_max_per_m:
                    lm = l_max_per_m[am]
                    if l1 > lm or l2 > lm:
                        continue
                channels.append((l1, m1, l2, m2))
    channels.sort(key=lambda x: (x[0] + x[2], x[0], abs(x[1]), x[1]))
    return channels


def compute_nuclear_coupling(
    l1p: int, l2p: int, l1: int, l2: int,
    m1: int, m2: int,
    alpha: np.ndarray, rho: float, Z: float = 1.0,
    Z_A: float = None, Z_B: float = None,
    rho_A: float = None, rho_B: float = None,
) -> np.ndarray:
    """
    Nuclear attraction coupling via algebraic multipole expansion.

    For nucleus A at distance R_A from origin (+z) and nucleus B at R_B (-z),
    the potential on electron i is:

        V_i = -Z_A sum_k f_k(s_i, rho_A) P_k
              -Z_B sum_k (-1)^k f_k(s_i, rho_B) P_k

    where f_k(s, rho) = (min(s,rho)/max(s,rho))^k / max(s,rho) and
    rho_X = R_X / R_e is the nuclear distance in hyperradial units.

    When rho_A = rho_B = rho (symmetric/midpoint origin), the two terms
    combine into coefficient -(Z_A + (-1)^k Z_B), recovering the original
    formula. When rho_A != rho_B (shifted origin), the radial factors
    differ and two separate sums are required.

    The 3j symbol (l' k l; 0 0 0) requires l'+k+l even, so for a given
    (l', l) pair, only one parity of k contributes.

    Parameters
    ----------
    l1p, l2p : int
        Bra channel angular momenta.
    l1, l2 : int
        Ket channel angular momenta.
    m1, m2 : int
        Shared magnetic quantum numbers (bra = ket, diagonal in m).
    alpha : ndarray
        Alpha grid points.
    rho : float
        R / (2 R_e). Used when rho_A, rho_B not specified.
    Z : float
        Nuclear charge per nucleus (used when Z_A, Z_B not specified).
    Z_A, Z_B : float or None
        Per-nucleus charges. If None, uses Z for both.
    rho_A, rho_B : float or None
        Nuclear distances from origin in R_e units. If None, uses rho
        for both (symmetric midpoint origin).

    Returns
    -------
    V : ndarray of shape (n_alpha,)
        Nuclear coupling at each alpha (charge function units, / R_e).
    """
    from geovac.angular_integrals import wigner3j

    # Resolve charges: backward compatible
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z
    # Resolve nuclear positions: backward compatible
    if rho_A is None:
        rho_A = rho
    if rho_B is None:
        rho_B = rho

    n_alpha = len(alpha)
    cos_a = np.cos(alpha)
    sin_a = np.sin(alpha)
    am1 = abs(m1)
    am2 = abs(m2)

    V = np.zeros(n_alpha)

    # Electron 1: couples (l1p, l1) with |m1|, acts when l2p == l2
    if l2p == l2:
        norm1 = sqrt((2 * l1p + 1) * (2 * l1 + 1))
        phase1 = (-1) ** am1
        s1 = cos_a

        # 3j requires l1p + k + l1 even → k has parity of (l1p + l1)
        k_min = abs(l1p - l1)
        k_max = l1p + l1
        required_parity = (l1p + l1) % 2
        if k_min % 2 != required_parity:
            k_min += 1

        # Precompute radial factors for both nuclei
        min_s1_A = np.minimum(s1, rho_A)
        max_s1_A = np.maximum(s1, rho_A)
        min_s1_B = np.minimum(s1, rho_B)
        max_s1_B = np.maximum(s1, rho_B)

        V1 = np.zeros(n_alpha)
        for k in range(k_min, k_max + 1, 2):
            g0 = wigner3j(l1p, k, l1, 0, 0, 0)
            if abs(g0) < 1e-15:
                continue
            gm = wigner3j(l1p, k, l1, -am1, 0, am1)
            coeff = g0 * gm
            if abs(coeff) < 1e-30:
                continue
            # Nucleus A at +R_A: -Z_A * f_k(s, rho_A)
            c_k_A = -Z_A * (min_s1_A / max_s1_A) ** k / max_s1_A
            # Nucleus B at -R_B: -Z_B * (-1)^k * f_k(s, rho_B)
            c_k_B = -Z_B * (-1)**k * (min_s1_B / max_s1_B) ** k / max_s1_B
            V1 += coeff * (c_k_A + c_k_B)

        V += phase1 * norm1 * V1

    # Electron 2: couples (l2p, l2) with |m2|, acts when l1p == l1
    if l1p == l1:
        norm2 = sqrt((2 * l2p + 1) * (2 * l2 + 1))
        phase2 = (-1) ** am2
        s2 = sin_a

        # 3j requires l2p + k + l2 even → k has parity of (l2p + l2)
        k_min = abs(l2p - l2)
        k_max = l2p + l2
        required_parity = (l2p + l2) % 2
        if k_min % 2 != required_parity:
            k_min += 1

        # Precompute radial factors for both nuclei
        min_s2_A = np.minimum(s2, rho_A)
        max_s2_A = np.maximum(s2, rho_A)
        min_s2_B = np.minimum(s2, rho_B)
        max_s2_B = np.maximum(s2, rho_B)

        V2 = np.zeros(n_alpha)
        for k in range(k_min, k_max + 1, 2):
            g0 = wigner3j(l2p, k, l2, 0, 0, 0)
            if abs(g0) < 1e-15:
                continue
            gm = wigner3j(l2p, k, l2, -am2, 0, am2)
            coeff = g0 * gm
            if abs(coeff) < 1e-30:
                continue
            # Nucleus A at +R_A: -Z_A * f_k(s, rho_A)
            c_k_A = -Z_A * (min_s2_A / max_s2_A) ** k / max_s2_A
            # Nucleus B at -R_B: -Z_B * (-1)^k * f_k(s, rho_B)
            c_k_B = -Z_B * (-1)**k * (min_s2_B / max_s2_B) ** k / max_s2_B
            V2 += coeff * (c_k_A + c_k_B)

        V += phase2 * norm2 * V2

    return V


def _ee_coupling(
    l1p: int, l2p: int, l1: int, l2: int,
    alpha: np.ndarray, l_max: int,
) -> np.ndarray:
    """
    Compute e-e repulsion coupling between channels (l1',l2') and (l1,l2).

    Uses the multipole expansion of 1/r12 in hyperspherical coordinates:
        1/r12 = (1/R_e) sum_k f_k(alpha) * [angular part]

    where f_k(alpha) = (min/max)^k / max with min=min(cos a, sin a),
    max=max(cos a, sin a).

    The angular coupling involves a double Gaunt integral:
        sum_k G(l1', k, l1) * G(l2', k, l2) * f_k / max

    For the He solver (l1=l2=l, l1'=l2'=l'), the Gaunt structure
    simplifies. Here we need the general two-index form.

    In H2, channels have independent (l1, l2), so the coupling is:
        W = sum_k sqrt((2l1'+1)(2l2'+1)(2l1+1)(2l2+1))/4
            * G(l1',k,l1) * G(l2',k,l2) * f_k(alpha) / max_sc

    Parameters
    ----------
    l1p, l2p, l1, l2 : int
        Channel angular momenta.
    alpha : ndarray
        Alpha grid.
    l_max : int
        Maximum l for Gaunt integral range.

    Returns
    -------
    W : ndarray of shape (n_alpha,)
        E-e coupling at each alpha point (charge function units / R_e).
    """
    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)
    min_sc = np.minimum(sin_a, cos_a)
    max_sc = np.maximum(sin_a, cos_a)

    W = np.zeros(len(alpha))

    k_min = max(abs(l1p - l1), abs(l2p - l2))
    k_max = min(l1p + l1, l2p + l2)

    for k in range(k_min, k_max + 1):
        g1 = gaunt_integral(l1p, k, l1)
        g2 = gaunt_integral(l2p, k, l2)
        if abs(g1) < 1e-15 or abs(g2) < 1e-15:
            continue
        f_k = (min_sc / max_sc) ** k
        W += g1 * g2 * f_k / max_sc

    norm = sqrt((2*l1p+1) * (2*l2p+1) * (2*l1+1) * (2*l2+1)) / 4.0
    W *= norm

    return W


def _ee_coupling_general(
    l1p: int, m1p: int, l2p: int, m2p: int,
    l1: int, m1: int, l2: int, m2: int,
    alpha: np.ndarray,
) -> np.ndarray:
    """
    Generalized e-e coupling for arbitrary m using Wigner 3j symbols.

    From the multipole expansion of 1/r12 in spherical harmonics:

        W = norm * sum_k (-1)^q *
            (l1' k l1; 0 0 0)(l1' k l1; -m1' -q m1) *
            (l2' k l2; 0 0 0)(l2' k l2; -m2' q m2) * f_k(a)/max

    where q = m1 - m1' is the azimuthal transfer quantum number, and
    norm = sqrt((2l1'+1)(2l1+1)(2l2'+1)(2l2+1)).

    Selection rules:
    - m1 + m2 = m1' + m2' (M conservation, automatic for M=0)
    - l1' + k + l1 even (from (000) 3j symbol)
    - |l1'-l1| <= k <= l1'+l1, |l2'-l2| <= k <= l2'+l2

    For m1=m1'=m2=m2'=0, reduces to _ee_coupling (verified algebraically:
    gaunt_integral = 2*(3j;000)^2, and the /4 norm cancels with the factor 4).

    Parameters
    ----------
    l1p, m1p, l2p, m2p : int
        Bra channel quantum numbers.
    l1, m1, l2, m2 : int
        Ket channel quantum numbers.
    alpha : ndarray
        Alpha grid.

    Returns
    -------
    W : ndarray of shape (n_alpha,)
        E-e coupling (charge function units / R_e).
    """
    from geovac.angular_integrals import wigner3j

    q = m1 - m1p
    # M conservation check: q must also equal m2p - m2
    if q != m2p - m2:
        return np.zeros(len(alpha))

    sin_a = np.sin(alpha)
    cos_a = np.cos(alpha)
    min_sc = np.minimum(sin_a, cos_a)
    max_sc = np.maximum(sin_a, cos_a)

    W = np.zeros(len(alpha))

    k_min = max(abs(l1p - l1), abs(l2p - l2))
    k_max = min(l1p + l1, l2p + l2)

    for k in range(k_min, k_max + 1):
        # Parity check via (l k l'; 0 0 0) — zero unless l+k+l' even
        g1_0 = wigner3j(l1p, k, l1, 0, 0, 0)
        if abs(g1_0) < 1e-15:
            continue
        g2_0 = wigner3j(l2p, k, l2, 0, 0, 0)
        if abs(g2_0) < 1e-15:
            continue

        g1_m = wigner3j(l1p, k, l1, -m1p, -q, m1)
        g2_m = wigner3j(l2p, k, l2, -m2p, q, m2)

        coeff = g1_0 * g1_m * g2_0 * g2_m
        if abs(coeff) < 1e-30:
            continue

        f_k = (min_sc / max_sc) ** k
        W += coeff * f_k / max_sc

    phase = (-1) ** q
    norm = sqrt((2*l1p+1) * (2*l1+1) * (2*l2p+1) * (2*l2+1))
    W *= phase * norm

    return W


def build_angular_hamiltonian(
    alpha_grid: np.ndarray,
    rho: float,
    R_e: float,
    l_max: int = 2,
    Z: float = 1.0,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
) -> np.ndarray:
    """
    Build the full coupled-channel angular Hamiltonian.

    Matrix structure: block (n_ch x n_ch) at each alpha grid point,
    with tridiagonal kinetic coupling within each channel along alpha.

    Total dimension: N = n_ch * n_alpha.
    Indexing: global = ch * n_alpha + i_alpha.

    Parameters
    ----------
    alpha_grid : ndarray of shape (n_alpha,)
        Interior alpha grid points (Dirichlet BCs at boundaries).
    rho : float
        R / (2 R_e).
    R_e : float
        Electronic hyperradius.
    l_max : int
        Maximum angular momentum.
    Z : float
        Nuclear charge per nucleus (homonuclear default).
    m_max : int
        Maximum |m| per electron. 0 = sigma-only (Phase 2),
        1 = sigma + pi (Phase 3).
    Z_A, Z_B : float or None
        Per-nucleus charges for heteronuclear systems. If None, uses Z.
    z0 : float
        Origin shift along internuclear axis. Positive = toward nucleus A.
        Nuclear positions: A at R/2 - z0, B at R/2 + z0 from origin.

    Returns
    -------
    H : ndarray of shape (N, N)
        Full angular Hamiltonian matrix.
    """
    # Resolve charges
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z
    homonuclear = (Z_A == Z_B)

    # Compute per-nucleus rho values for shifted origin
    # rho = R/(2*R_e), so R_A = R/2 - z0 → rho_A = rho - z0/R_e
    rho_A = rho - z0 / R_e if z0 != 0.0 else rho
    rho_B = rho + z0 / R_e if z0 != 0.0 else rho

    # --- Determine channel list ---
    if m_max == 0:
        channels_2 = _channel_list(l_max, homonuclear=homonuclear)
        n_ch = len(channels_2)
        # Convert to 4-tuple form for unified indexing below
        channels_4 = [(l1, 0, l2, 0) for l1, l2 in channels_2]
    else:
        channels_4 = _channel_list_extended(
            l_max, m_max, l_max_per_m, homonuclear=homonuclear,
        )
        n_ch = len(channels_4)

    n_alpha = len(alpha_grid)
    N = n_ch * n_alpha

    h = alpha_grid[1] - alpha_grid[0]  # uniform grid spacing

    cos_a = np.cos(alpha_grid)
    sin_a = np.sin(alpha_grid)

    H = np.zeros((N, N))

    kinetic_diag = 1.0 / h**2
    kinetic_off = -0.5 / h**2

    def idx(ch: int, i: int) -> int:
        return ch * n_alpha + i

    # --- Diagonal blocks: kinetic + centrifugal + Liouville for each channel ---
    # Centrifugal l(l+1)/(2 cos^2 a) depends on l, not m (m is absorbed
    # into the l(l+1) eigenvalue of the angular Laplacian on each sphere).
    for ic, (l1, m1, l2, m2) in enumerate(channels_4):
        V_cent = (0.5 * l1 * (l1 + 1) / cos_a**2
                  + 0.5 * l2 * (l2 + 1) / sin_a**2)
        V_liouville = -2.0
        V_diag = V_cent + V_liouville

        for i in range(n_alpha):
            ii = idx(ic, i)
            H[ii, ii] = kinetic_diag + V_diag[i]

        # Tridiagonal kinetic within channel
        for i in range(n_alpha - 1):
            ii = idx(ic, i)
            jj = idx(ic, i + 1)
            H[ii, jj] = kinetic_off
            H[jj, ii] = kinetic_off

    # --- Nuclear coupling (diagonal in m1, m2 due to axial symmetry) ---
    for ic, (l1p, m1p, l2p, m2p) in enumerate(channels_4):
        for jc, (l1, m1, l2, m2) in enumerate(channels_4):
            if jc < ic:
                continue  # fill symmetric

            # Nuclear potential is phi-independent: requires m1' = m1, m2' = m2
            if m1p != m1 or m2p != m2:
                continue

            V_nuc = compute_nuclear_coupling(
                l1p, l2p, l1, l2, m1, m2, alpha_grid, rho,
                Z=Z, Z_A=Z_A, Z_B=Z_B,
                rho_A=rho_A, rho_B=rho_B,
            )

            for i in range(n_alpha):
                ii = idx(ic, i)
                jj = idx(jc, i)
                val = R_e * V_nuc[i]
                H[ii, jj] += val
                if ic != jc:
                    H[jj, ii] += val

    # --- E-e coupling ---
    # For m_max=0: original Gaunt-based coupling (exact backward compatibility).
    # For m_max>0: generalized Wigner 3j coupling for all channel pairs.
    for ic, (l1p, m1p, l2p, m2p) in enumerate(channels_4):
        for jc, (l1, m1, l2, m2) in enumerate(channels_4):
            if jc < ic:
                continue

            if m_max == 0:
                W_ee = _ee_coupling(l1p, l2p, l1, l2, alpha_grid, l_max)
            else:
                W_ee = _ee_coupling_general(
                    l1p, m1p, l2p, m2p, l1, m1, l2, m2, alpha_grid
                )

            if np.max(np.abs(W_ee)) < 1e-15:
                continue

            for i in range(n_alpha):
                ii = idx(ic, i)
                jj = idx(jc, i)
                val = R_e * W_ee[i]
                H[ii, jj] += val
                if ic != jc:
                    H[jj, ii] += val

    return H


def solve_angular_multichannel(
    rho: float,
    R_e: float,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 200,
    n_eig: int = 1,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, list]:
    """
    Solve the coupled-channel angular eigenvalue problem.

    Parameters
    ----------
    rho : float
        R / (2 R_e).
    R_e : float
        Electronic hyperradius.
    l_max : int
        Maximum angular momentum per electron.
    Z : float
        Nuclear charge (homonuclear default).
    n_alpha : int
        FD grid points in alpha.
    n_eig : int
        Number of eigenvalues to return.
    m_max : int
        Maximum |m| per electron. 0 = sigma-only, 1 = sigma+pi.
    Z_A, Z_B : float or None
        Per-nucleus charges for heteronuclear systems.
    z0 : float
        Origin shift along internuclear axis.

    Returns
    -------
    mu : ndarray of shape (n_eig,)
        Lowest angular eigenvalues.
    vecs : ndarray of shape (n_eig, N)
        Eigenvectors.
    alpha_grid : ndarray
        Alpha grid.
    channels : list
        Channel labels: (l1, l2) if m_max=0, (l1, m1, l2, m2) if m_max>0.
    """
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z
    homonuclear = (Z_A == Z_B)

    h = (np.pi / 2) / (n_alpha + 1)
    alpha = (np.arange(n_alpha) + 1) * h

    H = build_angular_hamiltonian(alpha, rho, R_e, l_max, Z, m_max,
                                   l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0)

    evals, evecs = eigh(H)

    if m_max == 0:
        channels = _channel_list(l_max, homonuclear=homonuclear)
    else:
        channels = _channel_list_extended(
            l_max, m_max, l_max_per_m, homonuclear=homonuclear,
        )

    return evals[:n_eig], evecs[:, :n_eig].T, alpha, channels


def compute_adiabatic_curve_mc(
    R: float,
    R_e_grid: np.ndarray,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 200,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
) -> np.ndarray:
    """
    Compute adiabatic potential U(R_e) with multichannel angular solver.

    U(R_e; R) = [mu(R_e; R) + 15/8] / R_e^2

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        Electronic hyperradius grid.
    l_max : int
        Maximum angular momentum.
    Z : float
        Nuclear charge (homonuclear default).
    n_alpha : int
        FD grid points.
    m_max : int
        Maximum |m| per electron.
    Z_A, Z_B : float or None
        Per-nucleus charges for heteronuclear systems.
    z0 : float
        Origin shift along internuclear axis.

    Returns
    -------
    U : ndarray
        Effective potential (Ha).
    """
    n_Re = len(R_e_grid)
    mu_vals = np.zeros(n_Re)

    for i, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        evals, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max, Z, n_alpha, m_max=m_max,
            l_max_per_m=l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0,
        )
        mu_vals[i] = evals[0]

    U = (mu_vals + 15.0 / 8.0) / R_e_grid**2
    return U


def compute_adiabatic_curve_dboc(
    R: float,
    R_e_grid: np.ndarray,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 200,
    m_max: int = 0,
    n_states: int = 5,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute adiabatic potential U(R_e) with diagonal Born-Oppenheimer
    correction (DBOC).

    The DBOC adds a positive correction to the adiabatic potential that
    accounts for the non-adiabatic coupling between angular channels:

        U_corr(R_e) = U_0(R_e) + (1/2) sum_j |P_0j(R_e)|^2

    where P_0j = <phi_0 | d phi_j / dR_e> is the first derivative coupling.
    This correction is always positive, preventing variational violations
    that occur when the single-channel adiabatic approximation overcounts
    the angular channel mixing for large channel sets (m_max > 0).

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        Electronic hyperradius grid.
    l_max : int
        Maximum angular momentum.
    Z : float
        Nuclear charge.
    n_alpha : int
        FD grid points.
    m_max : int
        Maximum |m| per electron.
    n_states : int
        Number of adiabatic states for DBOC computation.

    Returns
    -------
    U_corrected : ndarray
        DBOC-corrected effective potential (Ha).
    U_bare : ndarray
        Uncorrected adiabatic potential (Ha).
    """
    n_Re = len(R_e_grid)
    mu_vals = np.zeros(n_Re)
    dboc_vals = np.zeros(n_Re)
    vecs_list: List[np.ndarray] = []

    for i, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        evals, evecs, _, _ = solve_angular_multichannel(
            rho, R_e, l_max, Z, n_alpha,
            n_eig=n_states, m_max=m_max,
            l_max_per_m=l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0,
        )
        mu_vals[i] = evals[0]

        # Phase-align eigenvectors with previous R_e
        if i > 0:
            for j in range(min(n_states, len(evals))):
                if np.dot(vecs_list[-1][j], evecs[j]) < 0:
                    evecs[j] *= -1

        vecs_list.append(evecs.copy())

    # Compute DBOC from central finite differences
    for i in range(1, n_Re - 1):
        h_left = R_e_grid[i] - R_e_grid[i - 1]
        h_right = R_e_grid[i + 1] - R_e_grid[i]
        h_avg = (h_left + h_right) / 2.0

        dboc = 0.0
        n_actual = min(n_states, vecs_list[i].shape[0])
        for j in range(1, n_actual):
            # P_0j = <phi_0(R_e) | d phi_j / dR_e>
            dphi_j = (vecs_list[i + 1][j] - vecs_list[i - 1][j]) / (2 * h_avg)
            P_0j = np.dot(vecs_list[i][0], dphi_j)
            dboc += P_0j**2

        dboc_vals[i] = dboc / 2.0

    # Boundary: extrapolate from interior
    dboc_vals[0] = dboc_vals[1]
    dboc_vals[-1] = dboc_vals[-2]

    U_bare = (mu_vals + 15.0 / 8.0) / R_e_grid**2
    U_corrected = U_bare + dboc_vals

    return U_corrected, U_bare


def solve_coupled_channel_radial(
    R: float,
    R_e_grid: np.ndarray,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 200,
    m_max: int = 0,
    n_coupled: int = 3,
    l_max_per_m: Union[dict, None] = None,
    n_Re_radial: int = 300,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
) -> Tuple[float, np.ndarray]:
    """
    Diabatic coupled-channel radial solver.

    Uses the angular eigenvectors at a reference R_e (the minimum of the
    adiabatic potential) as a fixed diabatic basis. The angular Hamiltonian
    is projected onto this basis at each R_e to give a smooth n_coupled ×
    n_coupled potential matrix V_ij(R_e). The coupled radial Schrödinger
    equation is then solved as a block-matrix eigenvalue problem.

    This properly handles non-adiabatic coupling between angular channels,
    avoiding the variational violations of the single-channel adiabatic
    approximation and the convergence issues of DBOC.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        R_e grid for angular sweeps (used to find reference point).
    l_max : int
        Maximum angular momentum per electron.
    Z : float
        Nuclear charge.
    n_alpha : int
        FD grid points for alpha.
    m_max : int
        Maximum |m| per electron.
    n_coupled : int
        Number of coupled radial channels.
    l_max_per_m : dict or None
        Per-|m| angular momentum limit.
    n_Re_radial : int
        Grid points for radial equation.
    R_e_min, R_e_max : float
        Radial boundaries.

    Returns
    -------
    E_elec : float
        Electronic energy eigenvalue.
    F : ndarray
        Radial wavefunction (n_coupled * n_Re_radial,).
    """
    # Step 1: Quick adiabatic sweep to find reference R_e
    n_Re_ang = len(R_e_grid)
    mu_vals = np.zeros(n_Re_ang)
    for i, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        evals, _, _, _ = solve_angular_multichannel(
            rho, R_e, l_max, Z, n_alpha, m_max=m_max,
            l_max_per_m=l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0,
        )
        mu_vals[i] = evals[0]

    U_adia = (mu_vals + 15.0 / 8.0) / R_e_grid**2
    i_ref = np.argmin(U_adia)
    R_e_ref = R_e_grid[i_ref]

    # Step 2: Get reference eigenvectors at R_e_ref
    rho_ref = R / (2.0 * R_e_ref)
    _, ref_vecs, _, _ = solve_angular_multichannel(
        rho_ref, R_e_ref, l_max, Z, n_alpha,
        n_eig=n_coupled, m_max=m_max, l_max_per_m=l_max_per_m,
        Z_A=Z_A, Z_B=Z_B, z0=z0,
    )
    # ref_vecs: shape (n_coupled, N_angular) — rows are eigenvectors

    # Step 3: Build diabatic potential matrix on the radial grid
    h_Re = (R_e_max - R_e_min) / (n_Re_radial + 1)
    R_e_radial = R_e_min + (np.arange(n_Re_radial) + 1) * h_Re

    # Precompute alpha grid (same as in solve_angular_multichannel)
    h_alpha = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h_alpha

    V_diabatic = np.zeros((n_Re_radial, n_coupled, n_coupled))

    for k, R_e in enumerate(R_e_radial):
        rho = R / (2.0 * R_e)
        H_ang = build_angular_hamiltonian(
            alpha_grid, rho, R_e, l_max, Z, m_max, l_max_per_m,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
        )
        # Project: V_ij = phi_i^T H_ang phi_j
        # H_ang is (N_ang x N_ang), ref_vecs is (n_coupled x N_ang)
        Hv = H_ang @ ref_vecs.T  # (N_ang x n_coupled)
        V_diabatic[k] = ref_vecs @ Hv  # (n_coupled x n_coupled)

        # Add 15/8 to the diagonal (hypercentrifugal term)
        for i in range(n_coupled):
            V_diabatic[k, i, i] = (V_diabatic[k, i, i] + 15.0 / 8.0) / R_e**2

        # Off-diagonal: also divide by R_e^2
        for i in range(n_coupled):
            for j in range(n_coupled):
                if i != j:
                    V_diabatic[k, i, j] /= R_e**2

    # Step 4: Build coupled-channel radial Hamiltonian
    # Dimension: (n_coupled * n_Re_radial) x (n_coupled * n_Re_radial)
    # Block structure: kinetic (tridiagonal within each channel) + V_diabatic
    N_total = n_coupled * n_Re_radial
    H_rad = np.zeros((N_total, N_total))

    kinetic_diag = 1.0 / h_Re**2
    kinetic_off = -0.5 / h_Re**2

    for ic in range(n_coupled):
        offset = ic * n_Re_radial
        # Kinetic: tridiagonal within channel ic
        for k in range(n_Re_radial):
            H_rad[offset + k, offset + k] = kinetic_diag
        for k in range(n_Re_radial - 1):
            H_rad[offset + k, offset + k + 1] = kinetic_off
            H_rad[offset + k + 1, offset + k] = kinetic_off

    # Potential coupling: at each R_e point, couple channels
    for k in range(n_Re_radial):
        for ic in range(n_coupled):
            for jc in range(n_coupled):
                H_rad[ic * n_Re_radial + k,
                      jc * n_Re_radial + k] += V_diabatic[k, ic, jc]

    # Step 5: Solve for lowest eigenvalue
    evals, evecs = eigh(H_rad, subset_by_index=[0, 0])
    E_elec = evals[0]
    F = evecs[:, 0]

    return E_elec, F


def solve_direct_2d(
    R: float,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 100,
    n_Re: int = 200,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    m_max: int = 0,
    l_max_per_m: Union[dict, None] = None,
    verbose: bool = True,
    Z_A: float = None,
    Z_B: float = None,
    z0: float = 0.0,
) -> float:
    """
    Direct 2D solver: solve the full (alpha, R_e) problem without
    the adiabatic approximation.

    Builds the complete Hamiltonian matrix in the product space of
    alpha grid × R_e grid × channel space, then finds the lowest
    eigenvalue using sparse iterative methods.

    H = T_Re ⊗ I_ang + [H_ang(R_e) + 15/8 I] / R_e²

    This avoids the adiabatic approximation entirely and gives a
    properly variational result.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    l_max : int
        Maximum angular momentum.
    Z : float
        Nuclear charge.
    n_alpha : int
        FD grid points for alpha.
    n_Re : int
        Grid points for R_e.
    R_e_min, R_e_max : float
        R_e boundaries.
    m_max : int
        Maximum |m| per electron.
    l_max_per_m : dict or None
        Per-|m| angular momentum limit.
    verbose : bool
        Print diagnostics.

    Returns
    -------
    E_elec : float
        Lowest electronic energy eigenvalue.
    """
    from scipy.sparse.linalg import eigsh

    # Resolve charges
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z
    homonuclear = (Z_A == Z_B)

    # Determine channels
    if m_max == 0:
        channels_2 = _channel_list(l_max, homonuclear=homonuclear)
        channels_4 = [(l1, 0, l2, 0) for l1, l2 in channels_2]
    else:
        channels_4 = _channel_list_extended(
            l_max, m_max, l_max_per_m, homonuclear=homonuclear,
        )
    n_ch = len(channels_4)

    # Grids
    h_alpha = (np.pi / 2) / (n_alpha + 1)
    alpha_grid = (np.arange(n_alpha) + 1) * h_alpha

    h_Re = (R_e_max - R_e_min) / (n_Re + 1)
    R_e_grid = R_e_min + (np.arange(n_Re) + 1) * h_Re

    N_ang = n_ch * n_alpha
    N_total = n_Re * N_ang

    if verbose:
        print(f"  Direct 2D solver: {N_total} x {N_total} sparse matrix")
        print(f"  ({n_Re} R_e points x {n_ch} channels x {n_alpha} alpha points)")

    # Build sparse Hamiltonian in COO format
    # Index: r * N_ang + (ch * n_alpha + ia)
    kinetic_diag_Re = 1.0 / h_Re**2
    kinetic_off_Re = -0.5 / h_Re**2

    # Angular blocks: [H_ang(R_e) + 15/8 I] / R_e^2
    # Collect all COO entries for efficiency
    rows_list: List[np.ndarray] = []
    cols_list: List[np.ndarray] = []
    vals_list: List[np.ndarray] = []

    # Add radial kinetic entries
    for r in range(n_Re):
        ang_idx = np.arange(N_ang)
        diag_idx = r * N_ang + ang_idx
        rows_list.append(diag_idx)
        cols_list.append(diag_idx)
        vals_list.append(np.full(N_ang, kinetic_diag_Re))

        if r < n_Re - 1:
            idx1 = r * N_ang + ang_idx
            idx2 = (r + 1) * N_ang + ang_idx
            rows_list.extend([idx1, idx2])
            cols_list.extend([idx2, idx1])
            vals_list.extend([np.full(N_ang, kinetic_off_Re)] * 2)

    for r, R_e in enumerate(R_e_grid):
        rho = R / (2.0 * R_e)
        H_ang = build_angular_hamiltonian(
            alpha_grid, rho, R_e, l_max, Z, m_max, l_max_per_m,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
        )
        # Add 15/8 to diagonal and scale
        np.fill_diagonal(H_ang, H_ang.diagonal() + 15.0 / 8.0)
        H_ang *= (1.0 / R_e**2)

        # Extract nonzero entries
        nz_i, nz_j = np.nonzero(np.abs(H_ang) > 1e-15)
        if len(nz_i) > 0:
            offset = r * N_ang
            rows_list.append(offset + nz_i)
            cols_list.append(offset + nz_j)
            vals_list.append(H_ang[nz_i, nz_j])

    from scipy.sparse import coo_matrix
    all_rows = np.concatenate(rows_list)
    all_cols = np.concatenate(cols_list)
    all_vals = np.concatenate(vals_list)
    H_csr = coo_matrix((all_vals, (all_rows, all_cols)),
                        shape=(N_total, N_total)).tocsr()
    if verbose:
        print(f"  Nonzeros: {H_csr.nnz}")

    # Find lowest eigenvalue using shift-invert mode.
    # Direct 'SA' (smallest algebraic) can find spurious boundary modes
    # at large matrix sizes. Shift-invert solves (H - sigma*I)^{-1}
    # and finds eigenvalues nearest to sigma.
    # Note: this function returns E_elec (no V_NN), so sigma must target
    # the electronic energy: E_elec ~ E_total - V_NN.
    # Estimate E_atoms for sigma: same logic as main solver
    if Z_A == Z_B:
        E_atoms_est = -Z_A**2
    else:
        Z_max = max(Z_A, Z_B)
        if Z_max == 2:
            E_atoms_est = -2.9037
        else:
            E_atoms_est = -Z_max**2 + 5 * Z_max / 8
    V_NN = Z_A * Z_B / R
    sigma = E_atoms_est - V_NN - 0.15
    try:
        evals, _ = eigsh(H_csr, k=3, sigma=sigma, which='LM')
        E_elec = np.min(evals)
    except Exception as e:
        if verbose:
            print(f"  Shift-invert failed ({e}), falling back to SA")
        evals, _ = eigsh(H_csr, k=6, which='SA')
        E_elec = evals[0]

    if verbose:
        print(f"  Eigenvalues found: {np.sort(evals)}")
        print(f"  sigma = {sigma:.4f}")

    return E_elec


def solve_level4_h2_multichannel(
    R: float,
    l_max: int = 2,
    Z: float = 1.0,
    n_alpha: int = 200,
    n_Re: int = 400,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    verbose: bool = True,
    m_max: int = 0,
    n_coupled: int = 1,
    l_max_per_m: Union[dict, None] = None,
    Z_A: float = None,
    Z_B: float = None,
    E_exact: float = None,
    D_e_exact: float = None,
    origin: str = 'midpoint',
) -> dict:
    """
    Full Level 4 multichannel solver for two-electron diatomics.

    Supports both homonuclear (H2, He2^2+) and heteronuclear (HeH+)
    systems via Z_A, Z_B parameters.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    l_max : int
        Maximum angular momentum per electron.
    Z : float
        Nuclear charge per nucleus (homonuclear default).
    n_alpha : int
        FD grid points for alpha.
    n_Re : int
        Grid points for hyperradial equation.
    R_e_min, R_e_max : float
        Hyperradial boundaries.
    verbose : bool
        Print diagnostics.
    m_max : int
        Maximum |m| per electron. 0 = sigma-only, 1 = sigma+pi.
    n_coupled : int
        Radial solver mode. 1 = single-channel adiabatic (default).
        -1 = direct 2D solver (properly variational, no adiabatic
        approximation, but slower). >1 = DBOC-corrected adiabatic
        (experimental, numerically unstable for large channel counts).
    Z_A, Z_B : float or None
        Per-nucleus charges. If None, uses Z for both (homonuclear).
    E_exact, D_e_exact : float or None
        Known exact values for comparison. If None, uses H2 defaults
        for Z_A=Z_B=1, otherwise omits comparison.
    origin : str or float
        Origin for hyperspherical coordinates.
        'midpoint': geometric center (default, z0=0).
        'charge_center': charge-weighted center, z0 = R(Z_A-Z_B)/(2(Z_A+Z_B)).
        float: explicit z0 value.

    Returns
    -------
    result : dict
        Keys: E_elec, E_total, D_e, D_e_pct, R, l_max, m_max, channels, etc.
    """
    import time
    t0 = time.time()

    # Resolve charges
    if Z_A is None:
        Z_A = Z
    if Z_B is None:
        Z_B = Z
    homonuclear = (Z_A == Z_B)

    # Resolve origin shift
    if isinstance(origin, (int, float)):
        z0 = float(origin)
    elif origin == 'charge_center':
        z0 = R * (Z_A - Z_B) / (2.0 * (Z_A + Z_B))
    else:  # 'midpoint'
        z0 = 0.0

    if m_max == 0:
        channels = _channel_list(l_max, homonuclear=homonuclear)
    else:
        channels = _channel_list_extended(
            l_max, m_max, l_max_per_m, homonuclear=homonuclear,
        )
    n_ch = len(channels)

    # Non-uniform R_e grid for adiabatic curve
    R_e_angular = np.concatenate([
        np.linspace(R_e_min, 1.0, 40),
        np.linspace(1.0, 3.0, 40),
        np.linspace(3.0, 6.0, 30),
        np.linspace(6.0, R_e_max, 20),
    ])
    R_e_angular = np.unique(R_e_angular)

    # System label
    if homonuclear:
        sys_label = f"Z={Z_A}" if Z_A != 1 else "H2"
    else:
        sys_label = f"Z_A={Z_A}, Z_B={Z_B}"

    if verbose:
        print(f"Level 4 multichannel solver for {sys_label} at R = {R:.4f} bohr")
        print(f"  l_max={l_max}, m_max={m_max}, n_ch={n_ch}")
        if not homonuclear:
            print(f"  Heteronuclear: Z_A={Z_A}, Z_B={Z_B}")
        if z0 != 0.0:
            print(f"  Origin shift: z0={z0:.4f} (R_A={R/2 - z0:.4f}, R_B={R/2 + z0:.4f})")
        if n_ch <= 15:
            print(f"  channels={channels}")
        print(f"  n_alpha={n_alpha}, n_Re={n_Re}")
        print(f"  Angular matrix size: {n_ch * n_alpha} x {n_ch * n_alpha}")
        if n_coupled > 1:
            print(f"  DBOC radial solver: {n_coupled} states")
        elif n_coupled == -1:
            print(f"  Direct 2D solver (no adiabatic approximation)")

    if n_coupled == -1:
        # --- Direct 2D solver (properly variational) ---
        E_elec = solve_direct_2d(
            R, l_max, Z, n_alpha, n_Re, R_e_min, R_e_max,
            m_max, l_max_per_m, verbose=verbose,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
        )
        F = np.zeros(n_Re)  # no radial wavefunction in 2D mode

        t2 = time.time()
        if verbose:
            print(f"  2D solve: {t2 - t0:.2f}s")

    elif n_coupled > 1:
        # --- DBOC-corrected adiabatic solver (experimental) ---
        if verbose:
            print(f"  Computing DBOC-corrected adiabatic curve "
                  f"({n_coupled} states) on "
                  f"{len(R_e_angular)} R_e points...")

        U_corrected, U_bare = compute_adiabatic_curve_dboc(
            R, R_e_angular, l_max, Z, n_alpha, m_max,
            n_states=n_coupled, l_max_per_m=l_max_per_m,
            Z_A=Z_A, Z_B=Z_B, z0=z0,
        )

        t1 = time.time()
        if verbose:
            i_min = np.argmin(U_corrected)
            print(f"  DBOC sweep: {t1 - t0:.2f}s")

        U_spline = CubicSpline(R_e_angular, U_corrected,
                               extrapolate=True)
        h_Re = (R_e_max - R_e_min) / (n_Re + 1)
        R_e_radial = R_e_min + (np.arange(n_Re) + 1) * h_Re
        V_radial = U_spline(R_e_radial)

        diag = np.ones(n_Re) / h_Re**2 + V_radial
        off_diag = -0.5 * np.ones(n_Re - 1) / h_Re**2

        evals, evecs = eigh_tridiagonal(
            diag, off_diag,
            select='i', select_range=(0, 0),
        )

        E_elec = evals[0]
        F = evecs[:, 0]

        t2 = time.time()
        if verbose:
            print(f"  Radial solve: {t2 - t1:.2f}s")

    else:
        # --- Single-channel adiabatic solver (original) ---
        if verbose:
            print(f"  Computing adiabatic curve on "
                  f"{len(R_e_angular)} R_e points...")

        U_angular = compute_adiabatic_curve_mc(
            R, R_e_angular, l_max, Z, n_alpha, m_max=m_max,
            l_max_per_m=l_max_per_m, Z_A=Z_A, Z_B=Z_B, z0=z0,
        )

        t1 = time.time()
        if verbose:
            i_min = np.argmin(U_angular)
            print(f"  Angular sweep: {t1 - t0:.2f}s")
            print(f"  U_min = {U_angular[i_min]:.6f} Ha at "
                  f"R_e = {R_e_angular[i_min]:.3f}")

        # Interpolate and solve radial equation
        U_spline = CubicSpline(R_e_angular, U_angular, extrapolate=True)

        h_Re = (R_e_max - R_e_min) / (n_Re + 1)
        R_e_radial = R_e_min + (np.arange(n_Re) + 1) * h_Re
        V_radial = U_spline(R_e_radial)

        diag = np.ones(n_Re) / h_Re**2 + V_radial
        off_diag = -0.5 * np.ones(n_Re - 1) / h_Re**2

        evals, evecs = eigh_tridiagonal(
            diag, off_diag,
            select='i', select_range=(0, 0),
        )

        E_elec = evals[0]
        F = evecs[:, 0]

        t2 = time.time()
        if verbose:
            print(f"  Radial solve: {t2 - t1:.2f}s")

    norm = np.sqrt(np.sum(F**2))
    if norm > 0:
        F /= norm

    V_NN = Z_A * Z_B / R
    E_total = E_elec + V_NN

    # Dissociation limit: where do the 2 electrons go?
    if homonuclear:
        # Each atom gets one electron: 2 × E(Z, 1e) = 2 × (-Z²/2) = -Z²
        E_atoms = -Z_A**2
    else:
        # Both electrons go to the larger nucleus (better binding)
        Z_max = max(Z_A, Z_B)
        if Z_max == 1:
            # H⁻ is barely bound; use H + H limit
            E_atoms = -1.0
        elif Z_max == 2:
            # He atom (exact nonrelativistic)
            E_atoms = -2.9037
        else:
            # Variational 2-electron estimate: E ≈ -Z² + 5Z/8
            E_atoms = -Z_max**2 + 5 * Z_max / 8

    D_e = E_atoms - E_total

    # Determine exact reference values
    if E_exact is None and D_e_exact is None:
        if homonuclear and Z_A == 1.0:
            # H2 defaults
            E_exact = -1.17447
            D_e_exact = 0.17447
        elif not homonuclear and Z_A == 2.0 and Z_B == 1.0:
            # HeH+ reference (Kolos & Peek 1976)
            E_exact = -2.9787
            D_e_exact = E_atoms - E_exact  # 0.0750
    elif D_e_exact is None and E_exact is not None:
        D_e_exact = E_atoms - E_exact

    D_e_pct = (D_e / D_e_exact * 100) if D_e_exact else None

    if verbose:
        print(f"\n  === Results (l_max={l_max}, {n_ch} channels) ===")
        if E_exact is not None:
            print(f"  E_total    = {E_total:.6f} Ha  (exact: {E_exact:.6f})")
            print(f"  D_e        = {D_e:.6f} Ha  (exact: {D_e_exact:.6f})")
            print(f"  D_e / exact = {D_e_pct:.1f}%")
            if homonuclear and Z_A == 1.0:
                if D_e_pct > 92.4:
                    print(f"  ** IMPROVES on Paper 12 Neumann V_ee (92.4%) **")
                else:
                    print(f"  Paper 12 Neumann V_ee: 92.4%")
        else:
            print(f"  E_total    = {E_total:.6f} Ha")
            print(f"  D_e        = {D_e:.6f} Ha")
            print(f"  E_atoms    = {E_atoms:.6f} Ha")
        print(f"  Total time: {t2 - t0:.2f}s")

    result = {
        'E_elec': E_elec,
        'E_total': E_total,
        'D_e': D_e,
        'D_e_pct': D_e_pct,
        'R': R,
        'Z': Z,
        'Z_A': Z_A,
        'Z_B': Z_B,
        'homonuclear': homonuclear,
        'l_max': l_max,
        'm_max': m_max,
        'n_ch': n_ch,
        'n_coupled': n_coupled,
        'channels': channels,
        'R_e_grid_angular': R_e_angular,
        'wavefunction': F,
        'E_atoms': E_atoms,
        'V_NN': V_NN,
        'z0': z0,
        'origin': origin,
    }

    if n_coupled == 1:
        result['U_adiabatic'] = U_angular
    h_Re = (R_e_max - R_e_min) / (n_Re + 1)
    result['R_e_grid_radial'] = R_e_min + (np.arange(n_Re) + 1) * h_Re

    return result
