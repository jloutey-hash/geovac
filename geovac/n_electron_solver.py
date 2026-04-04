"""
4-electron mol-frame hyperspherical solver for LiH.

Track AJ: First PK-free, composition-free GeoVac molecular PES
for a core-valence system.

Coordinates (Jacobi tree for 4 electrons):
    R_e = sqrt(r1^2 + r2^2 + r3^2 + r4^2)   [hyperradius]
    alpha1 = arctan(r2/r1)                      [pair 1-2]
    alpha2 = arctan(rho_34/rho_12)              [inter-pair]
    alpha3 = arctan(r4/r3)                      [pair 3-4]
    + 8 direction angles (theta_i, phi_i)

where rho_12 = sqrt(r1^2 + r2^2), rho_34 = sqrt(r3^2 + r4^2).

The electron radial magnitudes in terms of hyperangles:
    s1 = r1/R_e = cos(alpha2) * cos(alpha1)
    s2 = r2/R_e = cos(alpha2) * sin(alpha1)
    s3 = r3/R_e = sin(alpha2) * cos(alpha3)
    s4 = r4/R_e = sin(alpha2) * sin(alpha3)

At l_max=0 (all s-wave), the angular directions are integrated out
analytically, leaving a 3D eigenvalue problem in (alpha1, alpha2, alpha3).

The angular eigenvalue problem at fixed R_e:
    [T_alpha + R_e * C(Omega; rho)] Phi = mu Phi

where T_alpha is the hyperangular kinetic energy (including centrifugal),
C is the charge function (V_ee + V_nuc) in 1/R_e units, and rho = R/(2*R_e).

The effective potential is then:
    U(R_e) = [mu(R_e) + centrifugal_constant] / R_e^2

with centrifugal_constant = (3N-1)(3N-3)/8 = 11*9/8 = 99/8 for N=4.

Liouville substitution for the 3-hyperangle kinetic energy:
For each alpha_i with volume element weight w_i = sin^{p_i}(a)*cos^{q_i}(a),
the Liouville centrifugal potential (from absorbing the weight into the
wavefunction) is:
    V_cent_i = p_i(p_i-2)/(8*sin^2) + q_i(q_i-2)/(8*cos^2) - (p_i+q_i)^2/8

Pair hyperangles alpha1, alpha3 (p=q=2): V_cent = -2
Inter-pair hyperangle alpha2 (p=q=5): V_cent = 15/(8*sin^2) + 15/(8*cos^2) - 25/2

References:
  - Track AG scoping: geovac/n_electron_scope.py
  - Level 4 template: geovac/level4_multichannel.py
  - Lin, Phys. Rep. 257, 1 (1995)

Algebraic vs quadrature assessment:
  - Kinetic energy (T_alpha): ALGEBRAIC (FD stencil, exact)
  - Liouville centrifugal: ALGEBRAIC (exact formula from volume element exponents)
  - V_ee (6 pairs): QUADRATURE on 3D grid (1/max(s_i, s_j) at l=0)
    Could be algebraic via Gaunt integrals, but at l_max=0 only k=0 contributes,
    so it reduces to 1/max -- still piecewise on the 3D grid.
  - V_nuc (8 terms): QUADRATURE on 3D grid (-Z/max(s_i, rho) at l=0)
    Same as V_ee: piecewise on grid due to max() function.
  - Radial solve: ALGEBRAIC (FD tridiagonal) or spectral Laguerre
"""

import numpy as np
from scipy.linalg import eigh, eigh_tridiagonal
from scipy.interpolate import CubicSpline
from typing import Tuple, Dict, List, Optional
import time


# ==========================================================================
# COORDINATES: 4-electron hyperspherical tree
# ==========================================================================

def electron_radial_magnitudes(
    alpha1: np.ndarray, alpha2: np.ndarray, alpha3: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute s_i = r_i / R_e for each electron.

    Parameters
    ----------
    alpha1, alpha2, alpha3 : ndarray
        Hyperangles, each in [0, pi/2].

    Returns
    -------
    s1, s2, s3, s4 : ndarray
        Fractional radial distances (r_i / R_e).
        Satisfy s1^2 + s2^2 + s3^2 + s4^2 = 1.
    """
    c1, s1_a = np.cos(alpha1), np.sin(alpha1)
    c2, s2_a = np.cos(alpha2), np.sin(alpha2)
    c3, s3_a = np.cos(alpha3), np.sin(alpha3)

    s1 = c2 * c1  # r1/R_e
    s2 = c2 * s1_a  # r2/R_e
    s3 = s2_a * c3  # r3/R_e
    s4 = s2_a * s3_a  # r4/R_e

    return s1, s2, s3, s4


# ==========================================================================
# MULTICHANNEL: l_max >= 1 angular channel machinery
# ==========================================================================

def _enumerate_channels_4e(l_max: int) -> List[Tuple[int, int, int, int]]:
    """Enumerate sigma (M=0) channels for 4 electrons.

    Returns list of (l1, l2, l3, l4) with each l_i in [0, l_max].
    At l_max=0: 1 channel. At l_max=1: 16 channels.
    """
    channels = []
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            for l3 in range(l_max + 1):
                for l4 in range(l_max + 1):
                    channels.append((l1, l2, l3, l4))
    return channels


def _s4_permutation_matrices_channels(
    channels: List[Tuple[int, int, int, int]],
) -> List[np.ndarray]:
    """Build S_4 permutation matrices acting on the channel basis.

    S_4 acts by permuting the l-labels: sigma(l1,l2,l3,l4) = (l_{sigma(1)},...).

    Returns list of 24 permutation matrices (n_ch x n_ch).
    """
    from itertools import permutations

    n_ch = len(channels)
    ch_index = {ch: i for i, ch in enumerate(channels)}
    perms_s4 = list(permutations(range(4)))  # 24 permutations of (0,1,2,3)

    matrices = []
    for perm in perms_s4:
        P = np.zeros((n_ch, n_ch))
        for i, ch in enumerate(channels):
            permuted = (ch[perm[0]], ch[perm[1]], ch[perm[2]], ch[perm[3]])
            j = ch_index[permuted]
            P[i, j] = 1.0
        matrices.append(P)
    return matrices


def _s4_characters_22() -> List[int]:
    """Characters of the [2,2] irrep of S_4.

    S_4 has 5 conjugacy classes:
    e: 1 element, chi=2
    (12): 6 transpositions, chi=0
    (123): 8 3-cycles, chi=-1
    (1234): 6 4-cycles, chi=0
    (12)(34): 3 double transpositions, chi=2

    Returns chi(sigma) for each of the 24 permutations, in the
    order matching itertools.permutations.
    """
    from itertools import permutations

    perms = list(permutations(range(4)))
    chis = []
    for perm in perms:
        # Determine cycle type
        cycle_type = _cycle_type(perm)
        if cycle_type == (1, 1, 1, 1):
            chis.append(2)      # identity
        elif cycle_type == (2, 1, 1):
            chis.append(0)      # transposition
        elif cycle_type == (3, 1):
            chis.append(-1)     # 3-cycle
        elif cycle_type == (4,):
            chis.append(0)      # 4-cycle
        elif cycle_type == (2, 2):
            chis.append(2)      # double transposition
        else:
            chis.append(0)
    return chis


def _cycle_type(perm: Tuple[int, ...]) -> Tuple[int, ...]:
    """Compute cycle type of a permutation."""
    n = len(perm)
    visited = [False] * n
    cycles = []
    for i in range(n):
        if visited[i]:
            continue
        cycle_len = 0
        j = i
        while not visited[j]:
            visited[j] = True
            j = perm[j]
            cycle_len += 1
        cycles.append(cycle_len)
    return tuple(sorted(cycles, reverse=True))


def _build_s4_22_projector(
    channels: List[Tuple[int, int, int, int]],
) -> np.ndarray:
    """Build the S_4 [2,2] projection operator on the channel space.

    P_{[2,2]} = (dim/|G|) * sum_{g in S_4} chi_{[2,2]}(g) * D(g)
              = (2/24) * sum_{g} chi(g) * D(g)

    Returns the projector matrix (n_ch x n_ch).
    """
    perm_matrices = _s4_permutation_matrices_channels(channels)
    chis = _s4_characters_22()

    n_ch = len(channels)
    P = np.zeros((n_ch, n_ch))
    dim_irrep = 2
    order_G = 24

    for D, chi in zip(perm_matrices, chis):
        P += chi * D

    P *= dim_irrep / order_G
    return P


def _s4_channel_basis(
    channels: List[Tuple[int, int, int, int]],
) -> np.ndarray:
    """Compute the S_4 [2,2] basis vectors in channel space.

    Diagonalizes the 16x16 channel projector (trivial) and returns
    the eigenvectors with eigenvalue 1.

    Returns
    -------
    ch_basis : ndarray of shape (n_ch, rank)
        Column vectors spanning the [2,2] subspace in channel space.
    """
    projector = _build_s4_22_projector(channels)
    evals, evecs = np.linalg.eigh(projector)
    mask = evals > 0.5
    return evecs[:, mask]


def _s4_reduce_hamiltonian(
    H_full: np.ndarray,
    projector: np.ndarray,
    n_grid_total: int,
    n_ch: int,
) -> Tuple[np.ndarray, np.ndarray]:
    """Apply S_4 [2,2] projection to the multichannel Hamiltonian.

    Optimized: diagonalizes the small channel projector (n_ch x n_ch)
    to find the [2,2] basis, then constructs the reduced Hamiltonian
    by contracting channel indices. Avoids building and diagonalizing
    the full (n_ch * n_grid_total)^2 projector.

    Returns
    -------
    H_reduced : ndarray
        Hamiltonian in the S_4 [2,2] subspace.
    ch_basis : ndarray
        Channel-space basis vectors.
    """
    # Get channel-space basis (n_ch x rank) -- trivial 16x16 eigensolve
    evals_p, evecs_p = np.linalg.eigh(projector)
    mask = evals_p > 0.5
    ch_basis = evecs_p[:, mask]  # (n_ch, rank)
    rank = ch_basis.shape[1]

    if rank == 0:
        return np.zeros((0, 0)), ch_basis

    # Build the full basis: Q = ch_basis x I_grid
    # Q has shape (n_ch * n_grid, rank * n_grid)
    # Instead of building Q explicitly, compute H_reduced = Q^T H Q
    # by block operations.
    #
    # H_reduced[p*ng + i, q*ng + j] = sum_{a,b} ch_basis[a,p] * H[a*ng+i, b*ng+j] * ch_basis[b,q]
    #
    # H has block structure: H[a*ng+i, b*ng+j] is nonzero only when i==j for
    # the potential part (diagonal in grid), and has kinetic coupling for i~j
    # only within the same channel block (a==b).

    ng = n_grid_total
    dim_red = rank * ng
    H_reduced = np.zeros((dim_red, dim_red))

    for p in range(rank):
        for q in range(p, rank):
            # Compute the (p,q) block of H_reduced: a ng x ng matrix
            block = np.zeros((ng, ng))
            for a in range(n_ch):
                for b in range(n_ch):
                    w = ch_basis[a, p] * ch_basis[b, q]
                    if abs(w) < 1e-15:
                        continue
                    # Extract the (a,b) block of H_full
                    H_ab = H_full[a * ng:(a + 1) * ng, b * ng:(b + 1) * ng]
                    block += w * H_ab

            # Place in reduced Hamiltonian
            H_reduced[p * ng:(p + 1) * ng, q * ng:(q + 1) * ng] = block
            if p != q:
                H_reduced[q * ng:(q + 1) * ng, p * ng:(p + 1) * ng] = block.T

    return H_reduced, ch_basis


def _centrifugal_4e_multichannel(
    l_vals: Tuple[int, int, int, int],
    alpha1: np.ndarray, alpha2: np.ndarray, alpha3: np.ndarray,
) -> np.ndarray:
    """Compute centrifugal potential for a specific channel (l1,l2,l3,l4).

    Each electron i with angular momentum l_i contributes a centrifugal
    barrier l_i(l_i+1)/(2*s_i^2) = l_i(l_i+1)/(2*R_e^2*s_i^2).
    In the charge function formulation (multiplied by R_e), this gives
    angular centrifugal contributions.

    In the Liouville-transformed hyperangle coordinates, the individual
    electron centrifugal terms become:
    - Electron 1 (s1 = cos(a2)*cos(a1)):
      l1(l1+1) / (2*cos^2(a2)*cos^2(a1))
    - Electron 2 (s2 = cos(a2)*sin(a1)):
      l2(l2+1) / (2*cos^2(a2)*sin^2(a1))
    - Electron 3 (s3 = sin(a2)*cos(a3)):
      l3(l3+1) / (2*sin^2(a2)*cos^2(a3))
    - Electron 4 (s4 = sin(a2)*sin(a3)):
      l4(l4+1) / (2*sin^2(a2)*sin^2(a3))

    These are the additional centrifugal terms beyond the Liouville centrifugal
    (which is already included in the l_max=0 solver).

    Parameters
    ----------
    l_vals : tuple of 4 ints
        (l1, l2, l3, l4) angular momentum labels.
    alpha1, alpha2, alpha3 : ndarray
        Flat grid arrays.

    Returns
    -------
    V_cent_channel : ndarray
        Additional centrifugal potential on the grid.
    """
    l1, l2, l3, l4 = l_vals
    V = np.zeros_like(alpha1)

    ca1 = np.cos(alpha1)
    sa1 = np.sin(alpha1)
    ca2 = np.cos(alpha2)
    sa2 = np.sin(alpha2)
    ca3 = np.cos(alpha3)
    sa3 = np.sin(alpha3)

    # Each term: l(l+1) / (2 * s_i^2) where s_i is the fractional radius
    if l1 > 0:
        V += l1 * (l1 + 1) / (2.0 * ca2**2 * ca1**2)
    if l2 > 0:
        V += l2 * (l2 + 1) / (2.0 * ca2**2 * sa1**2)
    if l3 > 0:
        V += l3 * (l3 + 1) / (2.0 * sa2**2 * ca3**2)
    if l4 > 0:
        V += l4 * (l4 + 1) / (2.0 * sa2**2 * sa3**2)

    return V


def _vee_coupling_4e(
    ch_a: Tuple[int, int, int, int],
    ch_b: Tuple[int, int, int, int],
    s_vals: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    gaunt: np.ndarray,
) -> np.ndarray:
    """Compute V_ee coupling between two channels via Gaunt integrals.

    For each pair (i,j), the multipole expansion of 1/r_{ij} gives:
    C_ee^{ij} = sum_k (min(s_i,s_j)/max(s_i,s_j))^k / max(s_i,s_j)
                * G[l_i, k, l_i'] * G[l_j, k, l_j']
                * sqrt((2l_i+1)(2l_i'+1)) * sqrt((2l_j+1)(2l_j'+1)) / 4

    where the Gaunt factors arise from integrating P_{l_i}(cos theta_i)
    P_k(cos theta_{ij}) P_{l_i'}(cos theta_i) over the electron directions.

    The other two electrons (not in the pair) must have the same l: l_m = l_m'
    for m not in {i,j}.

    Parameters
    ----------
    ch_a, ch_b : tuple
        Channel labels (l1, l2, l3, l4).
    s_vals : tuple of 4 ndarrays
        (s1, s2, s3, s4) fractional radii on the grid.
    gaunt : ndarray
        Precomputed Gaunt integrals G[l, k, l'].

    Returns
    -------
    C_ee : ndarray
        V_ee coupling charge function on the grid (units 1/R_e).
    """
    from math import sqrt

    s1, s2, s3, s4 = s_vals
    s_arr = [s1, s2, s3, s4]
    N = len(s1)
    C = np.zeros(N)

    # 6 electron pairs
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

    for (i, j) in pairs:
        # Check: the electrons NOT in this pair must have same l
        others = [m for m in range(4) if m != i and m != j]
        if ch_a[others[0]] != ch_b[others[0]] or ch_a[others[1]] != ch_b[others[1]]:
            continue  # no coupling from this pair

        li_a, lj_a = ch_a[i], ch_a[j]
        li_b, lj_b = ch_b[i], ch_b[j]

        si = s_arr[i]
        sj = s_arr[j]
        min_sij = np.minimum(si, sj)
        max_sij = np.maximum(si, sj)

        # Sum over multipole order k
        k_min = max(abs(li_a - li_b), abs(lj_a - lj_b))
        k_max = min(li_a + li_b, lj_a + lj_b)
        # Also limited by Gaunt array size
        k_max = min(k_max, gaunt.shape[1] - 1)

        for k in range(k_min, k_max + 1):
            g_i = gaunt[li_a, k, li_b]
            g_j = gaunt[lj_a, k, lj_b]
            if abs(g_i) < 1e-15 or abs(g_j) < 1e-15:
                continue

            ratio_k = (min_sij / max_sij) ** k / max_sij
            norm = sqrt((2 * li_a + 1) * (2 * li_b + 1)) * \
                   sqrt((2 * lj_a + 1) * (2 * lj_b + 1)) / 4.0
            C += g_i * g_j * norm * ratio_k

    return C


def _vnuc_coupling_4e(
    ch_a: Tuple[int, int, int, int],
    ch_b: Tuple[int, int, int, int],
    s_vals: Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    rho_X: float, Z_X: float,
    gaunt: np.ndarray,
) -> np.ndarray:
    """Compute nuclear attraction coupling between two channels.

    For nucleus X at distance rho_X from origin (in R_e units),
    the attraction for electron i is:
    -Z_X * sum_k (min(s_i, rho_X)/max(s_i, rho_X))^k / max(s_i, rho_X)
           * G[l_i, k, l_i'] * sqrt((2l_i+1)(2l_i'+1)) / 2

    The other 3 electrons must have the same l: l_m = l_m' for m != i.

    Parameters
    ----------
    ch_a, ch_b : tuple
        Channel labels.
    s_vals : tuple of 4 ndarrays
        Fractional radii.
    rho_X : float
        Nuclear distance from origin in R_e units.
    Z_X : float
        Nuclear charge.
    gaunt : ndarray
        Precomputed Gaunt integrals.

    Returns
    -------
    C_nuc : ndarray
        Nuclear coupling charge function (units 1/R_e).
    """
    from math import sqrt

    s_arr = [s_vals[0], s_vals[1], s_vals[2], s_vals[3]]
    N = len(s_arr[0])
    C = np.zeros(N)

    for i in range(4):
        # Check: all electrons except i must have same l
        others = [m for m in range(4) if m != i]
        if any(ch_a[m] != ch_b[m] for m in others):
            continue

        li_a = ch_a[i]
        li_b = ch_b[i]
        si = s_arr[i]

        min_sr = np.minimum(si, rho_X)
        max_sr = np.maximum(si, rho_X)

        k_min = abs(li_a - li_b)
        k_max = min(li_a + li_b, gaunt.shape[1] - 1)

        for k in range(k_min, k_max + 1):
            g_val = gaunt[li_a, k, li_b]
            if abs(g_val) < 1e-15:
                continue

            ratio_k = (min_sr / max_sr) ** k / max_sr
            norm = sqrt((2 * li_a + 1) * (2 * li_b + 1)) / 2.0
            C += -Z_X * g_val * norm * ratio_k

    return C


def build_angular_hamiltonian_4e_parity_multichannel(
    rho_A: float, rho_B: float,
    R_e: float,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 10,
    l_max: int = 1,
) -> Tuple[np.ndarray, int, int]:
    """Build multichannel angular Hamiltonian with parity (pair antisymmetry).

    Uses half-grids for alpha1 and alpha3 (Dirichlet at pi/4) to enforce
    odd parity under (12) and (34) pair exchanges. This is the simplest
    antisymmetry consistent with the singlet ground state.

    All (l_max+1)^4 channels are included. The parity constraint is in
    the alpha space, not the channel space, so the (0,0,0,0) channel
    is present with its antisymmetric alpha wavefunction.

    Parameters
    ----------
    rho_A, rho_B : float
        Nuclear distances from origin in R_e units.
    R_e : float
        Hyperradius.
    Z_A, Z_B : float
        Nuclear charges.
    n_grid : int
        FD grid points per hyperangle.
    l_max : int
        Maximum angular momentum per electron.

    Returns
    -------
    H : ndarray
        Hamiltonian matrix (n_ch * n_grid^3, n_ch * n_grid^3).
    dim : int
        Matrix dimension.
    n_ch : int
        Number of channels.
    """
    from geovac.hyperspherical_angular import _precompute_gaunt

    channels = _enumerate_channels_4e(l_max)
    n_ch = len(channels)

    # Half-grid for alpha1, alpha3 (pair angles): (0, pi/4)
    # Full grid for alpha2 (inter-pair): (0, pi/2)
    h1 = (np.pi / 4.0) / (n_grid + 1)
    h2 = (np.pi / 2.0) / (n_grid + 1)

    a1 = (np.arange(n_grid) + 1) * h1
    a2 = (np.arange(n_grid) + 1) * h2
    a3 = (np.arange(n_grid) + 1) * h1

    n_grid_total = n_grid ** 3

    A1, A2, A3 = np.meshgrid(a1, a2, a3, indexing='ij')
    a1_flat = A1.ravel()
    a2_flat = A2.ravel()
    a3_flat = A3.ravel()

    s1, s2, s3, s4 = electron_radial_magnitudes(a1_flat, a2_flat, a3_flat)
    s_vals = (s1, s2, s3, s4)

    gaunt = _precompute_gaunt(l_max)

    N_full = n_ch * n_grid_total
    H = np.zeros((N_full, N_full))

    # 1D kinetic matrices
    T1d_pair = np.zeros((n_grid, n_grid))
    for i in range(n_grid):
        T1d_pair[i, i] = 1.0 / h1**2
        if i + 1 < n_grid:
            T1d_pair[i, i + 1] = -0.5 / h1**2
            T1d_pair[i + 1, i] = -0.5 / h1**2

    T1d_inter = np.zeros((n_grid, n_grid))
    for i in range(n_grid):
        T1d_inter[i, i] = 1.0 / h2**2
        if i + 1 < n_grid:
            T1d_inter[i, i + 1] = -0.5 / h2**2
            T1d_inter[i + 1, i] = -0.5 / h2**2

    I1 = np.eye(n_grid)
    T_3d = (np.kron(T1d_pair, np.kron(I1, I1)) +
            np.kron(I1, np.kron(T1d_inter, I1)) +
            np.kron(I1, np.kron(I1, T1d_pair)))

    # Liouville centrifugal
    V_liouville = np.zeros(n_grid_total)
    V_liouville += -2.0  # alpha1
    V_liouville += (15.0 / (8.0 * np.sin(a2_flat)**2) +
                    15.0 / (8.0 * np.cos(a2_flat)**2) - 12.5)
    V_liouville += -2.0  # alpha3

    # Fill diagonal blocks
    for a_idx, ch_a in enumerate(channels):
        block_start = a_idx * n_grid_total
        sl = slice(block_start, block_start + n_grid_total)

        H[sl, sl] = T_3d.copy()
        H[block_start:block_start + n_grid_total,
          block_start:block_start + n_grid_total] += np.diag(V_liouville)

        # Channel centrifugal
        V_ch = _centrifugal_4e_multichannel(ch_a, a1_flat, a2_flat, a3_flat)
        for k in range(n_grid_total):
            H[block_start + k, block_start + k] += V_ch[k]

    # Fill coupling blocks
    for a_idx, ch_a in enumerate(channels):
        for b_idx, ch_b in enumerate(channels):
            if b_idx < a_idx:
                continue

            C_ee = _vee_coupling_4e(ch_a, ch_b, s_vals, gaunt)
            C_nuc_A = _vnuc_coupling_4e(ch_a, ch_b, s_vals, rho_A, Z_A, gaunt)
            C_nuc_B = _vnuc_coupling_4e(ch_a, ch_b, s_vals, rho_B, Z_B, gaunt)
            C_total = R_e * (C_ee + C_nuc_A + C_nuc_B)

            if np.max(np.abs(C_total)) < 1e-15:
                continue

            a_start = a_idx * n_grid_total
            b_start = b_idx * n_grid_total
            for k in range(n_grid_total):
                H[a_start + k, b_start + k] += C_total[k]
                if a_idx != b_idx:
                    H[b_start + k, a_start + k] += C_total[k]

    return H, N_full, n_ch


def build_angular_hamiltonian_4e_multichannel(
    rho_A: float, rho_B: float,
    R_e: float,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 10,
    l_max: int = 1,
    s4_projection: bool = True,
) -> Tuple[np.ndarray, int, int]:
    """Build multichannel angular Hamiltonian for 4 electrons at l_max >= 1.

    The Hamiltonian has block structure:
    H_{ab}(i,j) = T_kinetic * delta_{ab} * delta_{ij}
                + V_cent_ab * delta_{ij}
                + R_e * (C_ee_{ab} + C_nuc_{ab}) * delta_{ij}

    where a, b are channel indices and i, j are 3D grid indices.
    The kinetic energy is the same for all channels (3D FD Laplacian).
    The centrifugal and potential terms are diagonal in the grid but may
    couple different channels.

    Parameters
    ----------
    rho_A, rho_B : float
        Nuclear distances from origin in R_e units.
    R_e : float
        Hyperradius.
    Z_A, Z_B : float
        Nuclear charges.
    n_grid : int
        FD grid points per hyperangle.
    l_max : int
        Maximum angular momentum per electron.
    s4_projection : bool
        If True, project onto S_4 [2,2] symmetry (singlet spatial).

    Returns
    -------
    H : ndarray
        Hamiltonian matrix. If s4_projection, in the [2,2] subspace.
    dim : int
        Matrix dimension.
    n_ch : int
        Number of channels (before projection).
    """
    from geovac.hyperspherical_angular import _precompute_gaunt

    channels = _enumerate_channels_4e(l_max)
    n_ch = len(channels)
    n_grid_total = n_grid ** 3

    # Build the S_4 projector FIRST to determine the subspace dimension
    if s4_projection:
        projector = _build_s4_22_projector(channels)
        # Find rank of projector
        evals_p = np.linalg.eigvalsh(projector)
        rank = int(np.sum(evals_p > 0.5))
        if rank == 0:
            return np.zeros((0, 0)), 0, n_ch

    # 3D FD grid (same as l_max=0 solver)
    h = np.pi / 2.0 / (n_grid + 1)
    a1 = (np.arange(n_grid) + 1) * h
    a2 = (np.arange(n_grid) + 1) * h
    a3 = (np.arange(n_grid) + 1) * h

    A1, A2, A3 = np.meshgrid(a1, a2, a3, indexing='ij')
    a1_flat = A1.ravel()
    a2_flat = A2.ravel()
    a3_flat = A3.ravel()

    s1, s2, s3, s4 = electron_radial_magnitudes(a1_flat, a2_flat, a3_flat)
    s_vals = (s1, s2, s3, s4)

    # Precompute Gaunt integrals
    gaunt = _precompute_gaunt(l_max)

    # Full Hamiltonian: n_ch * n_grid_total
    N_full = n_ch * n_grid_total
    H = np.zeros((N_full, N_full))

    # --- 1D FD kinetic matrix ---
    T1d = np.zeros((n_grid, n_grid))
    for i in range(n_grid):
        T1d[i, i] = 1.0 / h**2
        if i + 1 < n_grid:
            T1d[i, i + 1] = -0.5 / h**2
            T1d[i + 1, i] = -0.5 / h**2

    I1 = np.eye(n_grid)
    T_3d = (np.kron(T1d, np.kron(I1, I1)) +
            np.kron(I1, np.kron(T1d, I1)) +
            np.kron(I1, np.kron(I1, T1d)))

    # Liouville centrifugal (from volume element, same for all channels)
    V_liouville = np.zeros(n_grid_total)
    V_liouville += -2.0  # alpha1 (p=2, q=2)
    V_liouville += (15.0 / (8.0 * np.sin(a2_flat)**2) +
                    15.0 / (8.0 * np.cos(a2_flat)**2) - 12.5)  # alpha2
    V_liouville += -2.0  # alpha3 (p=2, q=2)

    # --- Fill diagonal blocks (kinetic + Liouville centrifugal + channel centrifugal) ---
    for a_idx, ch_a in enumerate(channels):
        block_start = a_idx * n_grid_total
        block_end = block_start + n_grid_total
        sl = slice(block_start, block_end)

        # Kinetic + Liouville centrifugal (same for all channels)
        H[sl, sl] = T_3d.copy()
        H[block_start:block_end, block_start:block_end] += np.diag(V_liouville)

        # Channel-specific centrifugal from individual electron angular momenta
        V_ch = _centrifugal_4e_multichannel(ch_a, a1_flat, a2_flat, a3_flat)
        for k in range(n_grid_total):
            H[block_start + k, block_start + k] += V_ch[k]

    # --- Fill coupling blocks (V_ee + V_nuc) ---
    for a_idx, ch_a in enumerate(channels):
        for b_idx, ch_b in enumerate(channels):
            if b_idx < a_idx:
                continue  # fill upper triangle, then symmetrize

            # V_ee coupling
            C_ee = _vee_coupling_4e(ch_a, ch_b, s_vals, gaunt)

            # V_nuc coupling (both nuclei)
            C_nuc_A = _vnuc_coupling_4e(ch_a, ch_b, s_vals, rho_A, Z_A, gaunt)
            C_nuc_B = _vnuc_coupling_4e(ch_a, ch_b, s_vals, rho_B, Z_B, gaunt)

            C_total = R_e * (C_ee + C_nuc_A + C_nuc_B)

            if np.max(np.abs(C_total)) < 1e-15:
                continue

            # Add to Hamiltonian block (diagonal in grid index)
            a_start = a_idx * n_grid_total
            b_start = b_idx * n_grid_total
            for k in range(n_grid_total):
                H[a_start + k, b_start + k] += C_total[k]
                if a_idx != b_idx:
                    H[b_start + k, a_start + k] += C_total[k]

    # --- S_4 [2,2] projection ---
    if s4_projection:
        H_red, basis = _s4_reduce_hamiltonian(H, projector, n_grid_total, n_ch)
        return H_red, H_red.shape[0], n_ch
    else:
        return H, N_full, n_ch


def solve_angular_4e_multichannel(
    rho_A: float, rho_B: float,
    R_e: float,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 10,
    l_max: int = 1,
    s4_projection: bool = True,
    symmetry: str = 'parity',
    n_states: int = 1,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Solve the multichannel angular eigenvalue problem for 4 electrons.

    Parameters
    ----------
    rho_A, rho_B : float
        Nuclear distances from origin in R_e units.
    R_e : float
        Hyperradius.
    Z_A, Z_B : float
        Nuclear charges.
    n_grid : int
        FD grid points per hyperangle.
    l_max : int
        Maximum angular momentum per electron.
    s4_projection : bool
        If True and symmetry='s4', project onto [2,2] singlet.
    symmetry : str
        'parity': half-grid pair antisymmetry (includes all channels).
        's4': channel-space [2,2] projection (excludes all-same-l channels).
        'none': no symmetry enforcement.
    n_states : int
        Number of eigenvalues to return.

    Returns
    -------
    evals : ndarray of shape (n_states,)
    evecs : ndarray of shape (n_states, dim) or None if dim=0
    """
    if symmetry == 'parity':
        H, dim, n_ch = build_angular_hamiltonian_4e_parity_multichannel(
            rho_A, rho_B, R_e, Z_A, Z_B, n_grid, l_max,
        )
    elif symmetry == 's4':
        H, dim, n_ch = build_angular_hamiltonian_4e_multichannel(
            rho_A, rho_B, R_e, Z_A, Z_B, n_grid, l_max, s4_projection=True,
        )
    else:
        H, dim, n_ch = build_angular_hamiltonian_4e_multichannel(
            rho_A, rho_B, R_e, Z_A, Z_B, n_grid, l_max, s4_projection=False,
        )

    if dim == 0:
        return np.array([0.0]), None

    # Use sparse eigensolve for large matrices, with dense fallback
    if dim > 2000 and n_states <= dim // 4:
        from scipy.sparse.linalg import eigsh as _eigsh
        from scipy.sparse import csr_matrix
        try:
            H_sp = csr_matrix(H)
            evals, evecs = _eigsh(H_sp, k=min(n_states, dim - 2), which='SA')
            idx = np.argsort(evals)
            evals = evals[idx[:n_states]]
            evecs = evecs[:, idx[:n_states]].T
        except Exception:
            # ARPACK convergence failure -- fall back to dense solver
            evals_all, evecs_all = eigh(H)
            evals = evals_all[:n_states]
            evecs = evecs_all[:, :n_states].T
    else:
        evals_all, evecs_all = eigh(H)
        evals = evals_all[:n_states]
        evecs = evecs_all[:, :n_states].T

    return evals, evecs


def compute_adiabatic_curve_4e_multichannel(
    R: float,
    R_e_grid: np.ndarray,
    Z_A: float = 3.0, Z_B: float = 1.0,
    z0: float = 0.0,
    n_grid: int = 8,
    l_max: int = 1,
    s4_projection: bool = True,
    symmetry: str = 'parity',
    pauli_correction: float = 0.0,
    verbose: bool = True,
) -> np.ndarray:
    """Compute adiabatic effective potential at l_max >= 1.

    U(R_e) = [mu(R_e) + 99/8 + pauli_correction] / R_e^2

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        Hyperradius grid.
    Z_A, Z_B : float
        Nuclear charges.
    z0 : float
        Origin shift.
    n_grid : int
        FD grid per hyperangle.
    l_max : int
        Maximum angular momentum.
    s4_projection : bool
        Project onto [2,2] symmetry (only for symmetry='s4').
    symmetry : str
        'parity', 's4', or 'none'.
    pauli_correction : float
        Additional centrifugal correction for Pauli exclusion.
    verbose : bool
        Print progress.

    Returns
    -------
    U : ndarray
        Effective potential (Ha).
    """
    n_Re = len(R_e_grid)
    mu_vals = np.zeros(n_Re)

    for i, R_e in enumerate(R_e_grid):
        R_A = R / 2.0 - z0
        R_B = R / 2.0 + z0
        rho_A = R_A / R_e
        rho_B = R_B / R_e

        evals, _ = solve_angular_4e_multichannel(
            rho_A, rho_B, R_e, Z_A, Z_B, n_grid, l_max,
            s4_projection=s4_projection, symmetry=symmetry,
        )
        mu_vals[i] = evals[0]

        if verbose and (i % 5 == 0 or i == n_Re - 1):
            print(f"  R_e={R_e:.3f}: mu={mu_vals[i]:.4f}")

    U = (mu_vals + CENTRIFUGAL_4E + pauli_correction) / R_e_grid**2
    return U


def solve_4e_lih_multichannel(
    R: float,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 8,
    l_max: int = 1,
    n_Re: int = 300,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    origin: str = 'midpoint',
    s4_projection: bool = True,
    symmetry: str = 'parity',
    pauli_correction: float = 0.0,
    verbose: bool = True,
) -> Dict:
    """Full 4-electron solver for LiH at l_max >= 1.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    n_grid : int
        FD grid points per hyperangle.
    l_max : int
        Maximum angular momentum per electron.
    n_Re : int
        Radial grid points.
    R_e_min, R_e_max : float
        Hyperradial range.
    origin : str
        'midpoint' or 'charge_center'.
    s4_projection : bool
        Project onto [2,2] symmetry (only for symmetry='s4').
    symmetry : str
        'parity', 's4', or 'none'.
    pauli_correction : float
        Additional centrifugal correction.
    verbose : bool
        Print diagnostics.

    Returns
    -------
    result : dict
    """
    t0 = time.time()

    if origin == 'charge_center':
        z0 = R * (Z_A - Z_B) / (2.0 * (Z_A + Z_B))
    else:
        z0 = 0.0

    channels = _enumerate_channels_4e(l_max)
    n_ch = len(channels)

    if verbose:
        print(f"4-electron LiH solver at R = {R:.4f} bohr")
        print(f"  Z_A={Z_A}, Z_B={Z_B}, l_max={l_max}")
        print(f"  n_grid={n_grid}, grid_total={n_grid**3}")
        print(f"  channels={n_ch}, symmetry={symmetry}")
        if z0 != 0.0:
            print(f"  Origin shift: z0={z0:.4f}")

    # R_e grid for adiabatic curve
    R_e_angular = np.concatenate([
        np.linspace(R_e_min, 1.5, 15),
        np.linspace(1.5, 4.0, 20),
        np.linspace(4.0, 8.0, 10),
        np.linspace(8.0, R_e_max, 5),
    ])
    R_e_angular = np.unique(R_e_angular)

    if verbose:
        print(f"  Computing adiabatic curve on {len(R_e_angular)} R_e points...")

    U_angular = compute_adiabatic_curve_4e_multichannel(
        R, R_e_angular, Z_A, Z_B, z0, n_grid, l_max,
        s4_projection=s4_projection, symmetry=symmetry,
        pauli_correction=pauli_correction, verbose=verbose,
    )

    t1 = time.time()
    if verbose:
        i_min = np.argmin(U_angular)
        print(f"  Angular sweep: {t1 - t0:.2f}s")
        print(f"  U_min = {U_angular[i_min]:.6f} Ha at R_e = {R_e_angular[i_min]:.3f}")

    # Interpolate and solve radial
    U_spline = CubicSpline(R_e_angular, U_angular, extrapolate=True)

    h_Re = (R_e_max - R_e_min) / (n_Re + 1)
    R_e_radial = R_e_min + (np.arange(n_Re) + 1) * h_Re
    V_radial = U_spline(R_e_radial)

    diag = np.ones(n_Re) / h_Re**2 + V_radial
    off_diag = -0.5 * np.ones(n_Re - 1) / h_Re**2

    evals_r, evecs_r = eigh_tridiagonal(
        diag, off_diag, select='i', select_range=(0, 0),
    )
    E_elec = evals_r[0]

    t2 = time.time()
    if verbose:
        print(f"  Radial solve: {t2 - t1:.2f}s")

    V_NN = Z_A * Z_B / R
    E_total = E_elec + V_NN

    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    D_e = E_atoms - E_total
    E_exact = -8.0705
    D_e_exact = E_atoms - E_exact

    if verbose:
        print(f"\n  === Results (l_max={l_max}, 4-electron) ===")
        print(f"  E_elec  = {E_elec:.6f} Ha")
        print(f"  V_NN    = {V_NN:.6f} Ha")
        print(f"  E_total = {E_total:.6f} Ha  (exact: {E_exact:.6f})")
        print(f"  E_atoms = {E_atoms:.6f} Ha")
        print(f"  D_e     = {D_e:.6f} Ha  (exact: {D_e_exact:.6f})")
        if D_e_exact > 0:
            print(f"  D_e %   = {D_e / D_e_exact * 100:.1f}%")
        print(f"  Total time: {t2 - t0:.2f}s")

    return {
        'E_elec': E_elec,
        'E_total': E_total,
        'D_e': D_e,
        'D_e_pct': (D_e / D_e_exact * 100) if D_e_exact > 0 else None,
        'R': R,
        'Z_A': Z_A,
        'Z_B': Z_B,
        'l_max': l_max,
        'n_grid': n_grid,
        'n_channels': n_ch,
        'angular_dim_full': n_ch * n_grid**3,
        'V_NN': V_NN,
        'E_atoms': E_atoms,
        'E_exact': E_exact,
        'D_e_exact': D_e_exact,
        'z0': z0,
        'origin': origin,
        'U_adiabatic': U_angular,
        'R_e_grid_angular': R_e_angular,
        'time_angular': t1 - t0,
        'time_radial': t2 - t1,
        'time_total': t2 - t0,
    }


def scan_pes_4e_lih_multichannel(
    R_values: np.ndarray = None,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 8,
    l_max: int = 1,
    origin: str = 'midpoint',
    s4_projection: bool = True,
    symmetry: str = 'parity',
    pauli_correction: float = 0.0,
    verbose: bool = True,
    output_dir: str = None,
) -> Dict:
    """Scan LiH PES at multiple R values with multichannel solver.

    Parameters
    ----------
    R_values : ndarray or None
        Internuclear distances.
    n_grid : int
        FD grid per hyperangle.
    l_max : int
        Maximum angular momentum.
    origin : str
        Origin choice.
    s4_projection : bool
        Apply [2,2] projection (only for symmetry='s4').
    symmetry : str
        'parity', 's4', or 'none'.
    pauli_correction : float
        Additional centrifugal correction.
    verbose : bool
        Print progress.
    output_dir : str or None
        Save results directory.

    Returns
    -------
    result : dict
    """
    if R_values is None:
        R_values = np.array([1.0, 1.5, 2.0, 2.5, 3.0, 3.2, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0])

    n_R = len(R_values)
    E_total = np.zeros(n_R)
    E_elec = np.zeros(n_R)
    D_e = np.zeros(n_R)
    times = np.zeros(n_R)

    t_start = time.time()
    dim_reported = False

    for j, R in enumerate(R_values):
        if verbose:
            print(f"\n{'='*60}")
            print(f"R = {R:.4f} bohr ({j+1}/{n_R})")
            print(f"{'='*60}")

        result = solve_4e_lih_multichannel(
            R, Z_A, Z_B, n_grid=n_grid, l_max=l_max,
            origin=origin, s4_projection=s4_projection,
            symmetry=symmetry, pauli_correction=pauli_correction,
            verbose=verbose,
        )
        E_total[j] = result['E_total']
        E_elec[j] = result['E_elec']
        D_e[j] = result['D_e']
        times[j] = result['time_total']

        if not dim_reported and verbose:
            print(f"  Angular dim (full): {result['angular_dim_full']}")
            dim_reported = True

    t_total = time.time() - t_start

    # Find minimum
    i_min = np.argmin(E_total)
    R_eq = R_values[i_min]
    E_min = E_total[i_min]
    D_e_min = D_e[i_min]

    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    E_exact = -8.0705
    D_e_exact = E_atoms - E_exact
    R_eq_expt = 3.015

    has_minimum = (i_min > 0 and i_min < n_R - 1)

    if verbose:
        print(f"\n{'='*60}")
        print(f"PES SCAN SUMMARY (l_max={l_max})")
        print(f"{'='*60}")
        print(f"  {'R':>6s}  {'E_total':>12s}  {'D_e':>10s}  {'time':>8s}")
        print(f"  {'bohr':>6s}  {'Ha':>12s}  {'Ha':>10s}  {'s':>8s}")
        print(f"  {'-'*6}  {'-'*12}  {'-'*10}  {'-'*8}")
        for j in range(n_R):
            marker = " <-- min" if j == i_min else ""
            print(f"  {R_values[j]:6.3f}  {E_total[j]:12.6f}  {D_e[j]:10.6f}  {times[j]:8.1f}{marker}")
        print()
        if has_minimum:
            R_eq_err = abs(R_eq - R_eq_expt) / R_eq_expt * 100
            print(f"  EQUILIBRIUM FOUND: R_eq = {R_eq:.3f} bohr")
            print(f"    (expt: {R_eq_expt:.3f}, error: {R_eq_err:.1f}%)")
            print(f"  E_min  = {E_min:.6f} Ha  (exact: {E_exact:.6f})")
            print(f"  D_e    = {D_e_min:.6f} Ha  (exact: {D_e_exact:.6f})")
            if D_e_exact > 0:
                print(f"  D_e %  = {D_e_min / D_e_exact * 100:.1f}%")
        else:
            print(f"  NO EQUILIBRIUM FOUND (minimum at boundary R={R_eq:.3f})")
        print(f"  Total time: {t_total:.1f}s")

    scan_result = {
        'R': R_values,
        'E_total': E_total,
        'E_elec': E_elec,
        'D_e': D_e,
        'R_eq': R_eq,
        'E_min': E_min,
        'D_e_min': D_e_min,
        'has_minimum': has_minimum,
        'R_eq_expt': R_eq_expt,
        'E_exact': E_exact,
        'D_e_exact': D_e_exact,
        'time_per_point': times,
        'time_total': t_total,
        'n_grid': n_grid,
        'l_max': l_max,
        'n_channels': (l_max + 1)**4,
        'origin': origin,
        'symmetry': symmetry,
    }

    if output_dir is not None:
        import os
        os.makedirs(output_dir, exist_ok=True)
        np.savez(os.path.join(output_dir, f'pes_4e_lih_lmax{l_max}.npz'), **scan_result)

    return scan_result


# ==========================================================================
# POTENTIALS: V_ee and V_nuc at l_max=0 (charge function form, units 1/R_e)
# ==========================================================================

def vee_charge_function_lmax0(
    s1: np.ndarray, s2: np.ndarray,
    s3: np.ndarray, s4: np.ndarray,
) -> np.ndarray:
    """Total e-e charge function for 4 electrons at l_max=0.

    At l_max=0, the angle-averaged 1/r_{ij} for s-wave electrons is:
        <1/r_{ij}> = 1/(R_e * max(s_i, s_j))

    The charge function (R_e times the potential) is:
        C_ee = sum_{i<j} 1/max(s_i, s_j)

    6 pairs: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4).

    Returns C_ee (dimensionless charge function, will be multiplied by R_e
    in the angular Hamiltonian).
    """
    C = np.zeros_like(s1)
    for si, sj in [(s1, s2), (s1, s3), (s1, s4),
                    (s2, s3), (s2, s4), (s3, s4)]:
        C += 1.0 / np.maximum(si, sj)
    return C


def vnuc_charge_function_lmax0(
    s1: np.ndarray, s2: np.ndarray,
    s3: np.ndarray, s4: np.ndarray,
    rho_A: float, rho_B: float,
    Z_A: float = 3.0, Z_B: float = 1.0,
) -> np.ndarray:
    """Total nuclear charge function for 4 electrons at l_max=0.

    At l_max=0, the angle-averaged nuclear attraction is:
        <-Z_X/|r_i - R_X|> = -Z_X / (R_e * max(s_i, rho_X))

    The charge function (R_e times the potential) is:
        C_nuc = sum_i [-Z_A/max(s_i, rho_A) - Z_B/max(s_i, rho_B)]

    4 electrons x 2 nuclei = 8 terms.

    Parameters
    ----------
    s1, s2, s3, s4 : ndarray
        Fractional radial distances.
    rho_A, rho_B : float
        Nuclear distances from origin in R_e units.
    Z_A, Z_B : float
        Nuclear charges.

    Returns C_nuc (dimensionless charge function).
    """
    C = np.zeros_like(s1)
    for si in [s1, s2, s3, s4]:
        C -= Z_A / np.maximum(si, rho_A)
        C -= Z_B / np.maximum(si, rho_B)
    return C


# ==========================================================================
# ANGULAR HAMILTONIAN: FD discretization in 3 hyperangles
# ==========================================================================

def build_angular_hamiltonian_4e(
    rho_A: float, rho_B: float,
    R_e: float,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 20,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build the 3D angular Hamiltonian for 4 electrons at l_max=0.

    The angular eigenvalue equation at fixed R_e:
        [T_alpha + R_e * C(Omega; rho)] Phi = mu Phi

    where C = C_ee + C_nuc is the charge function (dimensionless).

    T_alpha uses the Liouville substitution to absorb the volume element
    into the wavefunction, giving a standard eigenvalue problem with
    Dirichlet BCs at alpha_i = 0 and pi/2.

    Centrifugal potentials from Liouville (formula: p(p-2)/(8*sin^2) +
    q(q-2)/(8*cos^2) - (p+q)^2/8):
    - alpha1 (p=2, q=2): V = 0 + 0 - 2 = -2
    - alpha2 (p=5, q=5): V = 15/(8*sin^2) + 15/(8*cos^2) - 25/2
    - alpha3 (p=2, q=2): V = -2

    Parameters
    ----------
    rho_A, rho_B : float
        Nuclear distances from origin in R_e units.
    R_e : float
        Electronic hyperradius (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    n_grid : int
        Number of interior FD grid points per hyperangle.

    Returns
    -------
    H : ndarray of shape (N, N) where N = n_grid^3
        Angular Hamiltonian matrix.
    a1_flat, a2_flat, a3_flat : ndarray
        Flattened grid coordinates.
    """
    # Interior grid for each alpha in (0, pi/2) with Dirichlet BCs
    h = np.pi / 2.0 / (n_grid + 1)
    a1 = (np.arange(n_grid) + 1) * h
    a2 = (np.arange(n_grid) + 1) * h
    a3 = (np.arange(n_grid) + 1) * h

    N = n_grid ** 3

    # 3D grid via tensor product (ij indexing for consistent flattening)
    A1, A2, A3 = np.meshgrid(a1, a2, a3, indexing='ij')
    a1_flat = A1.ravel()
    a2_flat = A2.ravel()
    a3_flat = A3.ravel()

    # Electron radial magnitudes
    s1, s2, s3, s4 = electron_radial_magnitudes(a1_flat, a2_flat, a3_flat)

    # Build Hamiltonian
    H = np.zeros((N, N))

    # --- Kinetic energy: -1/2 d^2/da_i^2 for each dimension ---
    # 1D FD kinetic matrix (Dirichlet BCs)
    T1d = np.zeros((n_grid, n_grid))
    for i in range(n_grid):
        T1d[i, i] = 1.0 / h**2
        if i + 1 < n_grid:
            T1d[i, i + 1] = -0.5 / h**2
            T1d[i + 1, i] = -0.5 / h**2

    I1 = np.eye(n_grid)

    # T_total = T_a1 x I x I + I x T_a2 x I + I x I x T_a3
    H += np.kron(T1d, np.kron(I1, I1))
    H += np.kron(I1, np.kron(T1d, I1))
    H += np.kron(I1, np.kron(I1, T1d))

    # --- Centrifugal potential from Liouville substitution ---
    V_cent = np.zeros(N)

    # alpha1: p=2, q=2 -> V = -2
    V_cent += -2.0

    # alpha2: p=5, q=5 -> V = 15/(8*sin^2(a2)) + 15/(8*cos^2(a2)) - 25/2
    V_cent += 15.0 / (8.0 * np.sin(a2_flat)**2) + \
              15.0 / (8.0 * np.cos(a2_flat)**2) - 12.5

    # alpha3: p=2, q=2 -> V = -2
    V_cent += -2.0

    # --- Charge function: V_ee + V_nuc, multiplied by R_e ---
    C_ee = vee_charge_function_lmax0(s1, s2, s3, s4)
    C_nuc = vnuc_charge_function_lmax0(s1, s2, s3, s4, rho_A, rho_B, Z_A, Z_B)
    C_total = R_e * (C_ee + C_nuc)

    # Add potential to diagonal
    V_total = V_cent + C_total
    np.fill_diagonal(H, H.diagonal() + V_total)

    return H, a1_flat, a2_flat, a3_flat


def solve_angular_4e(
    rho_A: float, rho_B: float,
    R_e: float,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 20,
    n_states: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve the angular eigenvalue problem for 4 electrons at l_max=0.

    Parameters
    ----------
    rho_A, rho_B : float
        Nuclear distances from origin in R_e units.
    R_e : float
        Electronic hyperradius (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    n_grid : int
        FD grid points per hyperangle.
    n_states : int
        Number of eigenvalues to return.

    Returns
    -------
    evals : ndarray of shape (n_states,)
        Lowest eigenvalues mu.
    evecs : ndarray of shape (n_states, N)
        Eigenvectors.
    """
    H, _, _, _ = build_angular_hamiltonian_4e(
        rho_A, rho_B, R_e, Z_A, Z_B, n_grid,
    )

    N = H.shape[0]
    if n_states < N // 4 and N > 100:
        from scipy.sparse.linalg import eigsh
        from scipy.sparse import csr_matrix
        H_sp = csr_matrix(H)
        evals, evecs = eigsh(H_sp, k=min(n_states, N - 2), which='SA')
        idx = np.argsort(evals)
        evals = evals[idx]
        evecs = evecs[:, idx].T
    else:
        evals_all, evecs_all = eigh(H)
        evals = evals_all[:n_states]
        evecs = evecs_all[:, :n_states].T

    return evals, evecs


# ==========================================================================
# ADIABATIC CURVE AND PES
# ==========================================================================

# Centrifugal constant for N=4: (3N-1)(3N-3)/8 = 11*9/8 = 99/8
CENTRIFUGAL_4E = 99.0 / 8.0


def solve_angular_4e_parity(
    rho_A: float, rho_B: float,
    R_e: float,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 20,
    n_states: int = 1,
    parity: str = 'odd_12_34',
) -> Tuple[np.ndarray, np.ndarray]:
    """Solve angular eigenvalue problem with parity constraints.

    Enforces specific exchange symmetry by restricting the basis to
    functions with definite parity under hyperangle reflections:
    - alpha1 -> pi/2 - alpha1 (swaps electrons 1 and 2)
    - alpha3 -> pi/2 - alpha3 (swaps electrons 3 and 4)

    For 4 electrons at l_max=0:
    - 'odd_12_34': odd under both (12) and (34) exchanges.
      This creates a Pauli node at alpha1 = pi/4 and alpha3 = pi/4,
      preventing same-pair electron collapse.
      Corresponds to antisymmetric within each pair.

    Implementation: use a half-grid for alpha1 and alpha3 (0 to pi/4)
    with Dirichlet BCs at alpha = pi/4 (zero node for odd functions).

    Parameters
    ----------
    rho_A, rho_B : float
        Nuclear distances from origin in R_e units.
    R_e : float
        Electronic hyperradius.
    Z_A, Z_B : float
        Nuclear charges.
    n_grid : int
        FD grid points per hyperangle.
    n_states : int
        Number of eigenvalues.
    parity : str
        'odd_12_34': odd under (12) and (34).

    Returns
    -------
    evals : ndarray of shape (n_states,)
    evecs : ndarray of shape (n_states, N)
    """
    # Half-grid for alpha1 and alpha3: (0, pi/4) with Dirichlet at both ends
    # Full grid for alpha2: (0, pi/2) with Dirichlet at both ends
    h1 = (np.pi / 4.0) / (n_grid + 1)  # half-range for pair angles
    h2 = (np.pi / 2.0) / (n_grid + 1)  # full range for inter-pair angle

    a1 = (np.arange(n_grid) + 1) * h1   # in (0, pi/4)
    a2 = (np.arange(n_grid) + 1) * h2   # in (0, pi/2)
    a3 = (np.arange(n_grid) + 1) * h1   # in (0, pi/4)

    N = n_grid ** 3

    A1, A2, A3 = np.meshgrid(a1, a2, a3, indexing='ij')
    a1_flat = A1.ravel()
    a2_flat = A2.ravel()
    a3_flat = A3.ravel()

    s1, s2, s3, s4 = electron_radial_magnitudes(a1_flat, a2_flat, a3_flat)

    H = np.zeros((N, N))

    # 1D kinetic matrices with different grid spacings
    T1d_pair = np.zeros((n_grid, n_grid))
    for i in range(n_grid):
        T1d_pair[i, i] = 1.0 / h1**2
        if i + 1 < n_grid:
            T1d_pair[i, i + 1] = -0.5 / h1**2
            T1d_pair[i + 1, i] = -0.5 / h1**2

    T1d_inter = np.zeros((n_grid, n_grid))
    for i in range(n_grid):
        T1d_inter[i, i] = 1.0 / h2**2
        if i + 1 < n_grid:
            T1d_inter[i, i + 1] = -0.5 / h2**2
            T1d_inter[i + 1, i] = -0.5 / h2**2

    I1 = np.eye(n_grid)

    H += np.kron(T1d_pair, np.kron(I1, I1))   # T_a1
    H += np.kron(I1, np.kron(T1d_inter, I1))  # T_a2
    H += np.kron(I1, np.kron(I1, T1d_pair))   # T_a3

    # Centrifugal
    V_cent = np.zeros(N)
    V_cent += -2.0  # alpha1
    V_cent += 15.0 / (8.0 * np.sin(a2_flat)**2) + \
              15.0 / (8.0 * np.cos(a2_flat)**2) - 12.5  # alpha2
    V_cent += -2.0  # alpha3

    # Charge function
    C_ee = vee_charge_function_lmax0(s1, s2, s3, s4)
    C_nuc = vnuc_charge_function_lmax0(s1, s2, s3, s4, rho_A, rho_B, Z_A, Z_B)
    C_total = R_e * (C_ee + C_nuc)

    V_total = V_cent + C_total
    np.fill_diagonal(H, H.diagonal() + V_total)

    # Solve
    if n_states < N // 4 and N > 100:
        from scipy.sparse.linalg import eigsh
        from scipy.sparse import csr_matrix
        H_sp = csr_matrix(H)
        evals, evecs = eigsh(H_sp, k=min(n_states, N - 2), which='SA')
        idx = np.argsort(evals)
        evals = evals[idx]
        evecs = evecs[:, idx].T
    else:
        evals_all, evecs_all = eigh(H)
        evals = evals_all[:n_states]
        evecs = evecs_all[:, :n_states].T

    return evals, evecs


def compute_adiabatic_curve_4e(
    R: float,
    R_e_grid: np.ndarray,
    Z_A: float = 3.0, Z_B: float = 1.0,
    z0: float = 0.0,
    n_grid: int = 15,
    pauli_projection: str = 'node',
    verbose: bool = True,
) -> np.ndarray:
    """Compute adiabatic effective potential U(R_e) for 4 electrons.

    U(R_e) = [mu(R_e) + 99/8] / R_e^2

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        Electronic hyperradius grid.
    Z_A, Z_B : float
        Nuclear charges.
    z0 : float
        Origin shift along internuclear axis.
    n_grid : int
        FD grid points per hyperangle.
    pauli_projection : str
        Antisymmetry enforcement method:
        'none': no antisymmetry (unphysical, totally symmetric state).
        'parity': enforce odd parity under alpha1 -> pi/2 - alpha1
            and alpha3 -> pi/2 - alpha3 (excludes same-pair overlap).
        'node': use minimum nu=2 Casimir correction (Pauli centrifugal).
    verbose : bool
        Print progress.

    Returns
    -------
    U : ndarray
        Effective potential (Ha) at each R_e.
    """
    n_Re = len(R_e_grid)
    mu_vals = np.zeros(n_Re)

    for i, R_e in enumerate(R_e_grid):
        # Nuclear positions in R_e units
        R_A = R / 2.0 - z0  # distance from origin to nucleus A
        R_B = R / 2.0 + z0  # distance from origin to nucleus B
        rho_A = R_A / R_e
        rho_B = R_B / R_e

        if pauli_projection == 'parity':
            evals, _ = solve_angular_4e_parity(
                rho_A, rho_B, R_e, Z_A, Z_B, n_grid,
            )
        else:
            evals, _ = solve_angular_4e(
                rho_A, rho_B, R_e, Z_A, Z_B, n_grid,
            )
        mu_vals[i] = evals[0]

        if verbose and (i % 10 == 0 or i == n_Re - 1):
            print(f"  R_e={R_e:.3f}: mu={mu_vals[i]:.4f}")

    # Pauli centrifugal correction: the physical singlet [2,2] state has
    # nu_min = 2, so mu_free(2) = 12 is the minimum free eigenvalue.
    # Without antisymmetry, the solver finds the nu=0 state (mu_free=0).
    # Adding mu_pauli = 12 to the centrifugal corrects for this.
    if pauli_projection == 'node':
        mu_pauli = 12.0  # nu=2 Casimir - nu=0 Casimir = 12 - 0
        U = (mu_vals + CENTRIFUGAL_4E + mu_pauli) / R_e_grid**2
    else:
        U = (mu_vals + CENTRIFUGAL_4E) / R_e_grid**2
    return U


def solve_4e_lih(
    R: float,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 15,
    n_Re: int = 400,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    origin: str = 'midpoint',
    radial_method: str = 'fd',
    n_basis_radial: int = 25,
    alpha_radial: float = 0.8,
    pauli_projection: str = 'node',
    verbose: bool = True,
) -> Dict:
    """Full 4-electron solver for LiH at l_max=0.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges (Li=3, H=1).
    n_grid : int
        FD grid points per hyperangle (matrix dim = n_grid^3).
    n_Re : int
        Radial grid points (FD method).
    R_e_min, R_e_max : float
        Hyperradial boundaries.
    origin : str
        'midpoint' or 'charge_center'.
    radial_method : str
        'fd' or 'spectral'.
    n_basis_radial : int
        Laguerre basis functions (spectral method).
    alpha_radial : float
        Laguerre decay parameter (spectral method).
    pauli_projection : str
        Antisymmetry enforcement: 'none', 'parity', or 'node'.
    verbose : bool
        Print diagnostics.

    Returns
    -------
    result : dict
    """
    t0 = time.time()

    # Origin shift
    if origin == 'charge_center':
        z0 = R * (Z_A - Z_B) / (2.0 * (Z_A + Z_B))
    else:
        z0 = 0.0

    angular_dim = n_grid ** 3
    if verbose:
        print(f"4-electron LiH solver at R = {R:.4f} bohr")
        print(f"  Z_A={Z_A}, Z_B={Z_B}, l_max=0")
        print(f"  n_grid={n_grid}, angular_dim={angular_dim}")
        print(f"  pauli_projection={pauli_projection}")
        if z0 != 0.0:
            print(f"  Origin shift: z0={z0:.4f}")

    # R_e grid for adiabatic curve (non-uniform, denser near minimum)
    R_e_angular = np.concatenate([
        np.linspace(R_e_min, 1.5, 20),
        np.linspace(1.5, 4.0, 25),
        np.linspace(4.0, 8.0, 15),
        np.linspace(8.0, R_e_max, 10),
    ])
    R_e_angular = np.unique(R_e_angular)

    if verbose:
        print(f"  Computing adiabatic curve on {len(R_e_angular)} R_e points...")

    # Adiabatic angular sweep
    U_angular = compute_adiabatic_curve_4e(
        R, R_e_angular, Z_A, Z_B, z0, n_grid,
        pauli_projection=pauli_projection, verbose=verbose,
    )

    t1 = time.time()
    if verbose:
        i_min = np.argmin(U_angular)
        print(f"  Angular sweep: {t1 - t0:.2f}s")
        print(f"  U_min = {U_angular[i_min]:.6f} Ha at R_e = {R_e_angular[i_min]:.3f}")

    # Interpolate
    U_spline = CubicSpline(R_e_angular, U_angular, extrapolate=True)

    # Solve radial equation
    if radial_method == 'spectral':
        from geovac.hyperspherical_radial import (
            solve_radial_spectral as _solve_radial_spectral,
        )

        def V_safe(R_e_arr):
            R_arr = np.atleast_1d(R_e_arr)
            V = np.where(R_arr < R_e_max,
                         U_spline(np.clip(R_arr, R_e_min, R_e_max)),
                         0.0)
            return V

        evals, F_arr, R_eval = _solve_radial_spectral(
            V_safe, n_basis=n_basis_radial, alpha=alpha_radial,
            R_min=R_e_min, n_states=1,
        )
        E_elec = evals[0]
        F = F_arr[0]
    else:
        h_Re = (R_e_max - R_e_min) / (n_Re + 1)
        R_e_radial = R_e_min + (np.arange(n_Re) + 1) * h_Re
        V_radial = U_spline(R_e_radial)

        diag = np.ones(n_Re) / h_Re**2 + V_radial
        off_diag = -0.5 * np.ones(n_Re - 1) / h_Re**2

        evals_r, evecs_r = eigh_tridiagonal(
            diag, off_diag,
            select='i', select_range=(0, 0),
        )
        E_elec = evals_r[0]
        F = evecs_r[:, 0]

    t2 = time.time()
    if verbose:
        print(f"  Radial solve: {t2 - t1:.2f}s")

    # Nuclear repulsion
    V_NN = Z_A * Z_B / R
    E_total = E_elec + V_NN

    # Dissociation limit: Li(3e) + H(1e)
    # Li ground state: E_Li = -7.4781 Ha (exact nonrelativistic)
    # H ground state: E_H = -0.5 Ha
    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    D_e = E_atoms - E_total

    # LiH exact: -8.0705 Ha (nonrelativistic)
    E_exact = -8.0705
    D_e_exact = E_atoms - E_exact  # 0.0924 Ha

    if verbose:
        print(f"\n  === Results (l_max=0, 4-electron) ===")
        print(f"  E_elec  = {E_elec:.6f} Ha")
        print(f"  V_NN    = {V_NN:.6f} Ha")
        print(f"  E_total = {E_total:.6f} Ha  (exact: {E_exact:.6f})")
        print(f"  E_atoms = {E_atoms:.6f} Ha")
        print(f"  D_e     = {D_e:.6f} Ha  (exact: {D_e_exact:.6f})")
        if D_e_exact > 0:
            print(f"  D_e %   = {D_e / D_e_exact * 100:.1f}%")
        print(f"  Total time: {t2 - t0:.2f}s")

    result = {
        'E_elec': E_elec,
        'E_total': E_total,
        'D_e': D_e,
        'D_e_pct': (D_e / D_e_exact * 100) if D_e_exact > 0 else None,
        'R': R,
        'Z_A': Z_A,
        'Z_B': Z_B,
        'l_max': 0,
        'n_grid': n_grid,
        'angular_dim': angular_dim,
        'V_NN': V_NN,
        'E_atoms': E_atoms,
        'E_exact': E_exact,
        'D_e_exact': D_e_exact,
        'z0': z0,
        'origin': origin,
        'U_adiabatic': U_angular,
        'R_e_grid_angular': R_e_angular,
        'time_angular': t1 - t0,
        'time_radial': t2 - t1,
        'time_total': t2 - t0,
    }

    return result


def scan_pes_4e_lih(
    R_values: np.ndarray = None,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 15,
    origin: str = 'midpoint',
    pauli_projection: str = 'node',
    verbose: bool = True,
    output_dir: str = None,
) -> Dict:
    """Scan LiH PES at multiple R values.

    Parameters
    ----------
    R_values : ndarray or None
        Internuclear distances. If None, uses default 10-point scan.
    Z_A, Z_B : float
        Nuclear charges.
    n_grid : int
        FD grid points per hyperangle.
    origin : str
        'midpoint' or 'charge_center'.
    verbose : bool
        Print diagnostics.
    output_dir : str or None
        Directory to save results.

    Returns
    -------
    result : dict with keys 'R', 'E_total', 'D_e', etc.
    """
    if R_values is None:
        R_values = np.array([2.0, 2.5, 3.0, 3.2, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0])

    n_R = len(R_values)
    E_total = np.zeros(n_R)
    E_elec = np.zeros(n_R)
    D_e = np.zeros(n_R)
    times = np.zeros(n_R)

    t_start = time.time()

    for j, R in enumerate(R_values):
        if verbose:
            print(f"\n{'='*60}")
            print(f"R = {R:.4f} bohr ({j+1}/{n_R})")
            print(f"{'='*60}")

        result = solve_4e_lih(
            R, Z_A, Z_B, n_grid=n_grid,
            origin=origin, pauli_projection=pauli_projection,
            verbose=verbose,
        )
        E_total[j] = result['E_total']
        E_elec[j] = result['E_elec']
        D_e[j] = result['D_e']
        times[j] = result['time_total']

    t_total = time.time() - t_start

    # Find minimum
    i_min = np.argmin(E_total)
    R_eq = R_values[i_min]
    E_min = E_total[i_min]
    D_e_min = D_e[i_min]

    # Reference values
    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    E_exact = -8.0705
    D_e_exact = E_atoms - E_exact
    R_eq_expt = 3.015

    has_minimum = True
    if i_min == 0 or i_min == n_R - 1:
        has_minimum = False  # Minimum at boundary -- not a true minimum

    if verbose:
        print(f"\n{'='*60}")
        print(f"PES SCAN SUMMARY")
        print(f"{'='*60}")
        print(f"  {'R':>6s}  {'E_total':>12s}  {'D_e':>10s}  {'time':>6s}")
        print(f"  {'bohr':>6s}  {'Ha':>12s}  {'Ha':>10s}  {'s':>6s}")
        print(f"  {'-'*6}  {'-'*12}  {'-'*10}  {'-'*6}")
        for j in range(n_R):
            marker = " <-- min" if j == i_min else ""
            print(f"  {R_values[j]:6.3f}  {E_total[j]:12.6f}  {D_e[j]:10.6f}  {times[j]:6.1f}{marker}")
        print()
        if has_minimum:
            R_eq_err = abs(R_eq - R_eq_expt) / R_eq_expt * 100
            print(f"  EQUILIBRIUM FOUND: R_eq = {R_eq:.3f} bohr")
            print(f"    (expt: {R_eq_expt:.3f}, error: {R_eq_err:.1f}%)")
            print(f"  E_min  = {E_min:.6f} Ha  (exact: {E_exact:.6f})")
            print(f"  D_e    = {D_e_min:.6f} Ha  (exact: {D_e_exact:.6f})")
            if D_e_exact > 0:
                print(f"  D_e %  = {D_e_min / D_e_exact * 100:.1f}%")
        else:
            print(f"  NO EQUILIBRIUM FOUND (minimum at boundary R={R_eq:.3f})")
            print(f"  PES is {'monotonically decreasing' if i_min == 0 else 'monotonically increasing'}")
        print(f"  Total time: {t_total:.1f}s")

    scan_result = {
        'R': R_values,
        'E_total': E_total,
        'E_elec': E_elec,
        'D_e': D_e,
        'R_eq': R_eq,
        'E_min': E_min,
        'D_e_min': D_e_min,
        'has_minimum': has_minimum,
        'R_eq_expt': R_eq_expt,
        'E_exact': E_exact,
        'D_e_exact': D_e_exact,
        'time_per_point': times,
        'time_total': t_total,
        'n_grid': n_grid,
        'angular_dim': n_grid**3,
        'l_max': 0,
        'origin': origin,
    }

    if output_dir is not None:
        import os
        os.makedirs(output_dir, exist_ok=True)
        np.savez(os.path.join(output_dir, 'pes_4e_lih_lmax0.npz'), **scan_result)

    return scan_result


# ==========================================================================
# COUPLED-CHANNEL: Non-adiabatic radial solver for 4-electron system
# ==========================================================================
# Track AO: Extends the adiabatic solver with P-matrix coupling between
# multiple adiabatic channels, following the Level 3 pattern
# (algebraic_coupled_channel.py).
#
# Key difference from Level 3: The angular Hamiltonian H_ang(R_e) depends
# nonlinearly on R_e (through rho_X = R_X/R_e), so there is no R-independent
# dH/dR. The P-matrix must be computed by finite differences on the
# eigenvectors:
#   P_nu_mu(R_e) ≈ <Phi_nu(R_e)|Phi_mu(R_e + dR_e)> / dR_e
# ==========================================================================


def _enforce_sign_consistency_4e(
    vecs_prev: np.ndarray,
    vecs_curr: np.ndarray,
) -> np.ndarray:
    """Enforce eigenvector sign consistency between adjacent R_e points.

    Ensures <Phi_mu(R_e_i)|Phi_mu(R_e_{i-1})> > 0 for each channel mu.

    Parameters
    ----------
    vecs_prev : ndarray of shape (n_ch, dim)
        Sign-fixed eigenvectors at previous R_e point.
    vecs_curr : ndarray of shape (n_ch, dim)
        Eigenvectors at current R_e point.

    Returns
    -------
    vecs_fixed : ndarray of shape (n_ch, dim)
    """
    vecs_fixed = vecs_curr.copy()
    for mu in range(vecs_curr.shape[0]):
        if np.dot(vecs_prev[mu], vecs_curr[mu]) < 0:
            vecs_fixed[mu] *= -1
    return vecs_fixed


def compute_multichannel_angular_sweep_4e(
    R: float,
    R_e_grid: np.ndarray,
    Z_A: float = 3.0, Z_B: float = 1.0,
    z0: float = 0.0,
    n_grid: int = 6,
    l_max: int = 2,
    n_channels: int = 3,
    symmetry: str = 's4',
    verbose: bool = True,
) -> dict:
    """Compute multiple adiabatic curves + eigenvectors at each R_e point.

    Returns eigenvalues, eigenvectors, and grid for P-matrix computation.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    R_e_grid : ndarray
        Hyperradius grid.
    Z_A, Z_B : float
        Nuclear charges.
    z0 : float
        Origin shift.
    n_grid : int
        FD grid per hyperangle.
    l_max : int
        Maximum angular momentum per electron.
    n_channels : int
        Number of adiabatic channels to track.
    symmetry : str
        'parity' or 's4'.
    verbose : bool
        Print progress.

    Returns
    -------
    result : dict with keys:
        'mu': ndarray (n_channels, n_Re) -- adiabatic eigenvalues
        'evecs': list of ndarray (n_channels, dim) at each R_e
        'R_e_grid': ndarray
    """
    n_Re = len(R_e_grid)
    mu_all = np.zeros((n_channels, n_Re))
    evecs_all = []
    prev_vecs = None

    for i, R_e in enumerate(R_e_grid):
        R_A = R / 2.0 - z0
        R_B = R / 2.0 + z0
        rho_A = R_A / R_e
        rho_B = R_B / R_e

        evals, evecs = solve_angular_4e_multichannel(
            rho_A, rho_B, R_e, Z_A, Z_B, n_grid, l_max,
            s4_projection=True, symmetry=symmetry,
            n_states=n_channels,
        )

        # evecs shape: (n_channels, dim) -- rows are eigenvectors
        if evecs is None or len(evals) < n_channels:
            # Pad with zeros if fewer states available
            dim = evecs.shape[1] if evecs is not None else 1
            ev_padded = np.zeros((n_channels, dim))
            if evecs is not None:
                n_avail = min(len(evals), n_channels)
                ev_padded[:n_avail] = evecs[:n_avail]
                evals_padded = np.zeros(n_channels)
                evals_padded[:n_avail] = evals[:n_avail]
                evals_padded[n_avail:] = evals[-1] + 100.0  # push unused high
                evals = evals_padded
            evecs = ev_padded

        # Sign consistency
        if prev_vecs is not None:
            evecs = _enforce_sign_consistency_4e(prev_vecs, evecs)
        prev_vecs = evecs.copy()

        mu_all[:, i] = evals[:n_channels]
        evecs_all.append(evecs.copy())

        if verbose and (i % 5 == 0 or i == n_Re - 1):
            print(f"  R_e={R_e:.3f}: mu_0={mu_all[0, i]:.4f}, "
                  f"mu_1={mu_all[1, i]:.4f}" if n_channels > 1 else
                  f"  R_e={R_e:.3f}: mu_0={mu_all[0, i]:.4f}")

    return {
        'mu': mu_all,
        'evecs': evecs_all,
        'R_e_grid': R_e_grid,
    }


def compute_p_matrix_fd_4e(
    angular_data: dict,
    n_channels: int = 3,
) -> dict:
    """Compute P-matrix and DBOC via finite differences on eigenvectors.

    P_nu_mu(R_e_i) = <Phi_nu(R_e_i)|[Phi_mu(R_e_{i+1}) - Phi_mu(R_e_{i-1})]>
                     / (R_e_{i+1} - R_e_{i-1})

    using central differences. At boundaries, uses forward/backward differences.

    Also computes the DBOC (diagonal Born-Oppenheimer correction):
        Q_nu_nu = sum_mu |P_nu_mu|^2

    Parameters
    ----------
    angular_data : dict
        Output from compute_multichannel_angular_sweep_4e.
    n_channels : int
        Number of channels.

    Returns
    -------
    result : dict with keys:
        'P': ndarray (n_channels, n_channels, n_Re) -- P-matrix
        'DBOC': ndarray (n_channels, n_Re) -- diagonal Q
        'P_frob_norm': ndarray (n_Re,) -- Frobenius norm of P at each R_e
        'P_max_offdiag': ndarray (n_Re,) -- max |P_nu_mu| for nu != mu
    """
    R_e_grid = angular_data['R_e_grid']
    evecs_all = angular_data['evecs']
    n_Re = len(R_e_grid)
    n_ch = n_channels

    P = np.zeros((n_ch, n_ch, n_Re))

    for i in range(n_Re):
        if i == 0:
            # Forward difference
            dR = R_e_grid[1] - R_e_grid[0]
            for nu in range(n_ch):
                for mu in range(n_ch):
                    P[nu, mu, i] = np.dot(
                        evecs_all[i][nu], evecs_all[i + 1][mu]
                    ) / dR
                    if nu == mu:
                        # P_mu_mu should be ~0 (norm preservation)
                        # But FD gives <Phi|dPhi/dR> = d<Phi|Phi>/dR / 2 = 0
                        # Actually <Phi_mu|Phi_mu(R+dR)> ~ 1, so subtract 1
                        P[nu, mu, i] = (np.dot(
                            evecs_all[i][nu], evecs_all[i + 1][mu]
                        ) - 1.0) / dR
        elif i == n_Re - 1:
            # Backward difference
            dR = R_e_grid[i] - R_e_grid[i - 1]
            for nu in range(n_ch):
                for mu in range(n_ch):
                    if nu == mu:
                        P[nu, mu, i] = (1.0 - np.dot(
                            evecs_all[i - 1][nu], evecs_all[i][mu]
                        )) / dR
                    else:
                        P[nu, mu, i] = np.dot(
                            evecs_all[i][nu],
                            evecs_all[i][mu] - evecs_all[i - 1][mu]
                        ) / dR
        else:
            # Central difference
            dR = R_e_grid[i + 1] - R_e_grid[i - 1]
            for nu in range(n_ch):
                for mu in range(n_ch):
                    if nu == mu:
                        # Diagonal: should be ~0 by norm conservation
                        P[nu, mu, i] = 0.0
                    else:
                        P[nu, mu, i] = np.dot(
                            evecs_all[i][nu],
                            evecs_all[i + 1][mu] - evecs_all[i - 1][mu]
                        ) / dR

    # DBOC: Q_nu_nu = sum_mu |P_nu_mu|^2
    DBOC = np.zeros((n_ch, n_Re))
    for nu in range(n_ch):
        for mu in range(n_ch):
            if mu != nu:
                DBOC[nu] += P[nu, mu] ** 2

    # Diagnostic quantities
    P_frob_norm = np.zeros(n_Re)
    P_max_offdiag = np.zeros(n_Re)
    for i in range(n_Re):
        P_i = P[:, :, i].copy()
        np.fill_diagonal(P_i, 0.0)
        P_frob_norm[i] = np.linalg.norm(P_i, 'fro')
        if n_ch > 1:
            P_max_offdiag[i] = np.max(np.abs(P_i))

    return {
        'P': P,
        'DBOC': DBOC,
        'P_frob_norm': P_frob_norm,
        'P_max_offdiag': P_max_offdiag,
    }


def solve_4e_lih_coupled(
    R: float,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 6,
    l_max: int = 2,
    n_channels: int = 3,
    n_Re_angular: int = 50,
    n_Re_radial: int = 300,
    R_e_min: float = 0.3,
    R_e_max: float = 15.0,
    origin: str = 'midpoint',
    symmetry: str = 's4',
    q_mode: str = 'diagonal',
    verbose: bool = True,
) -> Dict:
    """Coupled-channel 4-electron solver for LiH.

    Extends the adiabatic solver with P-matrix non-adiabatic coupling
    between multiple adiabatic channels.

    Parameters
    ----------
    R : float
        Internuclear distance (bohr).
    Z_A, Z_B : float
        Nuclear charges.
    n_grid : int
        FD grid points per hyperangle.
    l_max : int
        Maximum angular momentum per electron.
    n_channels : int
        Number of coupled adiabatic channels.
    n_Re_angular : int
        R_e points for angular eigenvalue sweep.
    n_Re_radial : int
        R_e grid points for radial FD solver.
    R_e_min, R_e_max : float
        Hyperradial range.
    origin : str
        'midpoint' or 'charge_center'.
    symmetry : str
        'parity' or 's4'.
    q_mode : str
        Treatment of second-derivative coupling:
        'none': P coupling only (no Q).
        'diagonal': DBOC on V_eff diagonal.
        'full': Full Q from truncated P.
    verbose : bool
        Print diagnostics.

    Returns
    -------
    result : dict
    """
    t0 = time.time()

    if origin == 'charge_center':
        z0 = R * (Z_A - Z_B) / (2.0 * (Z_A + Z_B))
    else:
        z0 = 0.0

    channels = _enumerate_channels_4e(l_max)
    n_ch_total = len(channels)

    if verbose:
        print(f"Coupled-channel 4e LiH solver at R = {R:.4f} bohr")
        print(f"  Z_A={Z_A}, Z_B={Z_B}, l_max={l_max}")
        print(f"  n_grid={n_grid}, total channels={n_ch_total}")
        print(f"  Coupled channels: {n_channels}, q_mode={q_mode}")
        if z0 != 0.0:
            print(f"  Origin shift: z0={z0:.4f}")

    # R_e grid for angular sweep (non-uniform, denser near minimum)
    R_e_angular = np.concatenate([
        np.linspace(R_e_min, 1.0, n_Re_angular // 5),
        np.linspace(1.0, 2.5, n_Re_angular // 3),
        np.linspace(2.5, 5.0, n_Re_angular // 4),
        np.linspace(5.0, R_e_max, n_Re_angular // 5 + 1),
    ])
    R_e_angular = np.unique(R_e_angular)

    if verbose:
        print(f"  Angular sweep: {len(R_e_angular)} R_e points, "
              f"{n_channels} channels...")

    # Step 1: Multi-channel angular sweep
    angular_data = compute_multichannel_angular_sweep_4e(
        R, R_e_angular, Z_A, Z_B, z0, n_grid, l_max,
        n_channels=n_channels, symmetry=symmetry, verbose=verbose,
    )

    t1 = time.time()
    if verbose:
        print(f"  Angular sweep: {t1 - t0:.2f}s")

    # Step 2: Compute P-matrix from finite differences
    p_data = compute_p_matrix_fd_4e(angular_data, n_channels)
    P = p_data['P']
    DBOC = p_data['DBOC']

    if verbose:
        print(f"\n  P-matrix diagnostics:")
        print(f"    ||P||_F max: {np.max(p_data['P_frob_norm']):.6f}")
        print(f"    |P_offdiag| max: {np.max(p_data['P_max_offdiag']):.6f}")
        print(f"    DBOC_0 max: {np.max(DBOC[0]):.6f}")
        if n_channels > 1:
            print(f"    DBOC_1 max: {np.max(DBOC[1]):.6f}")

    # Step 3: Build effective potentials and splines
    mu_curves = angular_data['mu']
    R_e_grid = angular_data['R_e_grid']

    V_eff_all = np.zeros((n_channels, len(R_e_grid)))
    for ch in range(n_channels):
        V_eff_all[ch] = (mu_curves[ch] + CENTRIFUGAL_4E) / R_e_grid**2

    # Add DBOC to diagonal if requested
    if q_mode == 'diagonal':
        for ch in range(n_channels):
            V_eff_all[ch] += DBOC[ch]

    # Create splines
    V_eff_splines = []
    for ch in range(n_channels):
        V_eff_splines.append(
            CubicSpline(R_e_grid, V_eff_all[ch], extrapolate=True)
        )

    P_splines = [[None] * n_channels for _ in range(n_channels)]
    for nu in range(n_channels):
        for mu in range(n_channels):
            P_splines[nu][mu] = CubicSpline(
                R_e_grid, P[nu, mu], extrapolate=True
            )

    # Q splines for 'full' mode
    Q_splines_obj = None
    include_Q = (q_mode == 'full')
    if include_Q:
        Q_splines_obj = [[None] * n_channels for _ in range(n_channels)]
        for i_re in range(len(R_e_grid)):
            pass  # Q computed inline below
        # Q_mu_nu = sum_kappa P_mu_kappa * P_kappa_nu
        Q_data = np.zeros((n_channels, n_channels, len(R_e_grid)))
        for i_re in range(len(R_e_grid)):
            P_i = P[:, :, i_re]
            Q_data[:, :, i_re] = P_i @ P_i.T
        for nu in range(n_channels):
            for mu in range(n_channels):
                Q_splines_obj[nu][mu] = CubicSpline(
                    R_e_grid, Q_data[nu, mu], extrapolate=True
                )

    # Step 4: Solve coupled radial equation using FD
    from scipy.sparse import lil_matrix, csr_matrix
    from scipy.sparse.linalg import eigsh as sparse_eigsh

    n_ch = n_channels
    h_Re = (R_e_max - R_e_min) / (n_Re_radial + 1)
    R_e_radial = R_e_min + (np.arange(n_Re_radial) + 1) * h_Re

    # Evaluate on radial grid
    V_grid = np.zeros((n_ch, n_Re_radial))
    for ch in range(n_ch):
        V_grid[ch] = V_eff_splines[ch](R_e_radial)

    P_grid = np.zeros((n_ch, n_ch, n_Re_radial))
    for nu in range(n_ch):
        for mu in range(n_ch):
            P_grid[nu, mu] = P_splines[nu][mu](R_e_radial)

    Q_grid = np.zeros((n_ch, n_ch, n_Re_radial))
    if include_Q and Q_splines_obj is not None:
        for nu in range(n_ch):
            for mu in range(n_ch):
                Q_grid[nu, mu] = Q_splines_obj[nu][mu](R_e_radial)

    # Build block-tridiagonal sparse matrix
    dim_rad = n_Re_radial * n_ch
    H_rad = lil_matrix((dim_rad, dim_rad))

    kinetic_diag = 1.0 / h_Re**2
    kinetic_off = -0.5 / h_Re**2

    for i in range(n_Re_radial):
        for mu_idx in range(n_ch):
            row = i * n_ch + mu_idx

            # Diagonal: kinetic + V_eff (+ Q diagonal if full)
            H_rad[row, row] = kinetic_diag + V_grid[mu_idx, i]
            if include_Q:
                H_rad[row, row] -= 0.5 * Q_grid[mu_idx, mu_idx, i]

            # Off-diagonal channels at same grid point (Q coupling)
            if include_Q:
                for nu_idx in range(n_ch):
                    if nu_idx != mu_idx:
                        col = i * n_ch + nu_idx
                        H_rad[row, col] -= 0.5 * Q_grid[mu_idx, nu_idx, i]

        # Off-diagonal in R_e: coupling to i+1
        if i < n_Re_radial - 1:
            for mu_idx in range(n_ch):
                row = i * n_ch + mu_idx
                for nu_idx in range(n_ch):
                    col = (i + 1) * n_ch + nu_idx
                    val = 0.0
                    if mu_idx == nu_idx:
                        val += kinetic_off
                    val -= P_grid[mu_idx, nu_idx, i] / (2.0 * h_Re)
                    if abs(val) > 1e-16:
                        H_rad[row, col] = val

        # Off-diagonal in R_e: coupling to i-1
        if i > 0:
            for mu_idx in range(n_ch):
                row = i * n_ch + mu_idx
                for nu_idx in range(n_ch):
                    col = (i - 1) * n_ch + nu_idx
                    val = 0.0
                    if mu_idx == nu_idx:
                        val += kinetic_off
                    val += P_grid[mu_idx, nu_idx, i] / (2.0 * h_Re)
                    if abs(val) > 1e-16:
                        H_rad[row, col] = val

    H_csr = csr_matrix(H_rad)

    # Shift-invert eigensolve
    sigma_target = np.min(V_grid[0]) - 1.0
    n_states = min(5, dim_rad - 2)
    evals_r, evecs_r = sparse_eigsh(
        H_csr, k=n_states, sigma=sigma_target, which='LM',
    )
    idx = np.argsort(evals_r)
    evals_r = evals_r[idx]
    evecs_r = evecs_r[:, idx]
    E_elec = evals_r[0]

    # Channel weights for ground state
    F_gs = evecs_r[:, 0].reshape(n_Re_radial, n_ch)
    weights = np.array([
        h_Re * np.sum(F_gs[:, ch]**2) for ch in range(n_ch)
    ])
    if weights.sum() > 0:
        weights /= weights.sum()

    t2 = time.time()

    # Also compute single-channel (adiabatic) for comparison
    V_single = (mu_curves[0] + CENTRIFUGAL_4E) / R_e_grid**2
    V_single_spline = CubicSpline(R_e_grid, V_single, extrapolate=True)
    V_single_grid = V_single_spline(R_e_radial)
    diag_single = np.ones(n_Re_radial) / h_Re**2 + V_single_grid
    off_single = -0.5 * np.ones(n_Re_radial - 1) / h_Re**2
    evals_single, _ = eigh_tridiagonal(
        diag_single, off_single, select='i', select_range=(0, 0),
    )
    E_single = evals_single[0]

    V_NN = Z_A * Z_B / R
    E_total = E_elec + V_NN
    E_total_single = E_single + V_NN

    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    D_e = E_atoms - E_total
    D_e_single = E_atoms - E_total_single
    E_exact = -8.0705
    D_e_exact = E_atoms - E_exact

    if verbose:
        print(f"\n  === Coupled-Channel Results (l_max={l_max}, "
              f"{n_channels} channels, q_mode={q_mode}) ===")
        print(f"  E_elec (coupled)    = {E_elec:.6f} Ha")
        print(f"  E_elec (adiabatic)  = {E_single:.6f} Ha")
        print(f"  E_total (coupled)   = {E_total:.6f} Ha")
        print(f"  E_total (adiabatic) = {E_total_single:.6f} Ha")
        print(f"  Delta E             = {E_total - E_total_single:+.6f} Ha")
        direction = "UP (correct)" if E_total > E_total_single else "DOWN"
        print(f"  Direction:           {direction}")
        print(f"  D_e (coupled)       = {D_e:.6f} Ha  (exact: {D_e_exact:.6f})")
        print(f"  D_e (adiabatic)     = {D_e_single:.6f} Ha")
        if D_e_exact > 0:
            print(f"  D_e% (coupled)      = {D_e / D_e_exact * 100:.1f}%")
            print(f"  D_e% (adiabatic)    = {D_e_single / D_e_exact * 100:.1f}%")
        print(f"  Channel weights:     {weights}")
        print(f"  Total time:          {t2 - t0:.2f}s")

    return {
        'E_elec': E_elec,
        'E_total': E_total,
        'E_elec_single': E_single,
        'E_total_single': E_total_single,
        'D_e': D_e,
        'D_e_single': D_e_single,
        'D_e_pct': (D_e / D_e_exact * 100) if D_e_exact > 0 else None,
        'D_e_pct_single': (D_e_single / D_e_exact * 100) if D_e_exact > 0 else None,
        'R': R,
        'Z_A': Z_A,
        'Z_B': Z_B,
        'l_max': l_max,
        'n_grid': n_grid,
        'n_channels': n_channels,
        'q_mode': q_mode,
        'V_NN': V_NN,
        'E_atoms': E_atoms,
        'E_exact': E_exact,
        'D_e_exact': D_e_exact,
        'z0': z0,
        'origin': origin,
        'channel_weights': weights,
        'P_data': p_data,
        'mu_curves': mu_curves,
        'R_e_grid_angular': R_e_grid,
        'DBOC': DBOC,
        'time_angular': t1 - t0,
        'time_total': t2 - t0,
    }


def scan_pes_4e_lih_coupled(
    R_values: np.ndarray = None,
    Z_A: float = 3.0, Z_B: float = 1.0,
    n_grid: int = 6,
    l_max: int = 2,
    n_channels: int = 3,
    origin: str = 'midpoint',
    symmetry: str = 's4',
    q_mode: str = 'diagonal',
    verbose: bool = True,
    output_dir: str = None,
) -> Dict:
    """Scan LiH PES with coupled-channel solver.

    Parameters
    ----------
    R_values : ndarray or None
        Internuclear distances.
    n_grid : int
        FD grid per hyperangle.
    l_max : int
        Maximum angular momentum.
    n_channels : int
        Number of coupled channels.
    origin : str
        Origin choice.
    symmetry : str
        'parity' or 's4'.
    q_mode : str
        'none', 'diagonal', or 'full'.
    verbose : bool
        Print progress.
    output_dir : str or None
        Save results directory.

    Returns
    -------
    result : dict
    """
    if R_values is None:
        R_values = np.array([0.7, 1.0, 1.5, 2.0, 3.0])

    n_R = len(R_values)
    E_total = np.zeros(n_R)
    E_total_single = np.zeros(n_R)
    D_e = np.zeros(n_R)
    D_e_single = np.zeros(n_R)
    times = np.zeros(n_R)

    t_start = time.time()

    for j, R_val in enumerate(R_values):
        if verbose:
            print(f"\n{'='*60}")
            print(f"R = {R_val:.4f} bohr ({j+1}/{n_R})")
            print(f"{'='*60}")

        result = solve_4e_lih_coupled(
            R_val, Z_A, Z_B, n_grid=n_grid, l_max=l_max,
            n_channels=n_channels, origin=origin,
            symmetry=symmetry, q_mode=q_mode, verbose=verbose,
        )
        E_total[j] = result['E_total']
        E_total_single[j] = result['E_total_single']
        D_e[j] = result['D_e']
        D_e_single[j] = result['D_e_single']
        times[j] = result['time_total']

    t_total = time.time() - t_start

    # Find minima
    i_min = np.argmin(E_total)
    i_min_s = np.argmin(E_total_single)

    E_Li_exact = -7.4781
    E_H_exact = -0.5
    E_atoms = E_Li_exact + E_H_exact
    E_exact = -8.0705
    D_e_exact = E_atoms - E_exact
    R_eq_expt = 3.015

    if verbose:
        print(f"\n{'='*60}")
        print(f"COUPLED-CHANNEL PES SUMMARY (l_max={l_max}, "
              f"{n_channels} ch, q_mode={q_mode})")
        print(f"{'='*60}")
        print(f"  {'R':>6s}  {'E_coupled':>12s}  {'E_adiab':>12s}  "
              f"{'Delta':>10s}  {'time':>8s}")
        print(f"  {'bohr':>6s}  {'Ha':>12s}  {'Ha':>12s}  "
              f"{'Ha':>10s}  {'s':>8s}")
        print(f"  {'-'*6}  {'-'*12}  {'-'*12}  {'-'*10}  {'-'*8}")
        for j in range(n_R):
            delta = E_total[j] - E_total_single[j]
            marker = " <-- min" if j == i_min else ""
            print(f"  {R_values[j]:6.3f}  {E_total[j]:12.6f}  "
                  f"{E_total_single[j]:12.6f}  {delta:+10.6f}  "
                  f"{times[j]:8.1f}{marker}")
        print()
        print(f"  Coupled:   R_eq ~ {R_values[i_min]:.3f}, "
              f"D_e = {D_e[i_min]:.6f} Ha "
              f"({D_e[i_min]/D_e_exact*100:.1f}%)")
        print(f"  Adiabatic: R_eq ~ {R_values[i_min_s]:.3f}, "
              f"D_e = {D_e_single[i_min_s]:.6f} Ha "
              f"({D_e_single[i_min_s]/D_e_exact*100:.1f}%)")
        print(f"  Exact:     D_e = {D_e_exact:.6f} Ha")
        print(f"  Total time: {t_total:.1f}s")

    scan_result = {
        'R': R_values,
        'E_total': E_total,
        'E_total_single': E_total_single,
        'D_e': D_e,
        'D_e_single': D_e_single,
        'n_channels': n_channels,
        'l_max': l_max,
        'q_mode': q_mode,
        'time_per_point': times,
        'time_total': t_total,
    }

    if output_dir is not None:
        import os
        os.makedirs(output_dir, exist_ok=True)
        np.savez(
            os.path.join(output_dir, f'pes_4e_coupled_lmax{l_max}_nch{n_channels}.npz'),
            **scan_result,
        )

    return scan_result
