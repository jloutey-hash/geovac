"""
Spectral compression analysis for the N-electron angular eigenproblem.

Track AK: Determines the spectral compression ratio for 4-electron
mol-frame hyperspherical coordinates, following the Track K pattern
that achieved 20x compression at Level 4 (2 electrons).

THE TRACK K PATTERN (Level 4, 2 electrons):
  1. Free spectrum: SO(6) Casimir eigenvalues nu(nu+4)/2 (exact, integer)
  2. Spectral basis: Jacobi polynomials P_k^{(a,b)} in hyperangle alpha
  3. V_ee coupling: precomputed once via Gaunt-type integrals (rho-independent)
  4. V_nuc coupling: scales with rho (computed at each rho-point)
  5. Result: 1000 FD points -> 50 spectral basis functions -> 20x compression

THE 4-ELECTRON EXTENSION:
  1. Free spectrum: SO(12) Casimir eigenvalues nu(nu+10)/2
  2. Spectral basis: tensor product of Jacobi polynomials in 3 hyperangles
  3. V_ee coupling: 6 pairs (intra-pair + cross-pair)
  4. V_nuc coupling: 8 terms (4 electrons x 2 nuclei)

This module does NOT implement a solver. It analyzes spectral basis
dimensions and compression ratios for feasibility assessment.
"""

import numpy as np
from math import comb, factorial
from typing import Dict, List, Tuple, Optional
from itertools import product as itertools_product

from geovac.n_electron_scope import (
    casimir_eigenvalue,
    _so_d_degeneracy,
    four_electron_channel_count_atomic,
    four_electron_channel_count_molecular,
    level4_channel_count,
    _count_m_zero_states,
    _count_L0_couplings,
)


# ==========================================================================
# (a) LEVEL 4 REFERENCE: verify Track K pattern
# ==========================================================================

def level4_spectral_dimensions(l_max: int, n_basis: int = 10,
                                homonuclear: bool = True) -> Dict:
    """Compute spectral dimensions for Level 4 (2-electron) angular problem.

    This reproduces the Track K result for comparison with the 4-electron case.

    Parameters
    ----------
    l_max : int
        Maximum angular momentum per electron.
    n_basis : int
        Number of Jacobi polynomial basis functions per channel.
    homonuclear : bool
        Whether the molecule is homonuclear (gerade constraint l1+l2 even).

    Returns
    -------
    dict with dimension analysis.
    """
    # Channel list: (l1, l2) with l1+l2 even for homonuclear, sigma only
    channels = []
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            if homonuclear and (l1 + l2) % 2 != 0:
                continue
            channels.append((l1, l2))

    n_ch = len(channels)

    # FD dimension: n_ch * n_alpha_grid (typically 200)
    n_alpha_grid = 200
    fd_dim = n_ch * n_alpha_grid

    # Spectral dimension: n_ch * n_basis
    spectral_dim = n_ch * n_basis

    # Free spectrum: SO(6) Casimir
    free_spectrum = []
    for l1, l2 in channels:
        for k in range(n_basis):
            nu = l1 + l2 + 2 * k
            mu = nu * (nu + 4) / 2.0
            free_spectrum.append({
                'channel': (l1, l2), 'k': k, 'nu': nu, 'mu_free': mu
            })

    return {
        'n_electrons': 2,
        'angular_dim': 5,  # S^5
        'n_hyperangles': 1,  # alpha
        'n_direction_angles': 4,  # theta1, phi1, theta2, phi2
        'n_channels': n_ch,
        'channels': channels,
        'fd_dim': fd_dim,
        'spectral_dim': spectral_dim,
        'compression': fd_dim / spectral_dim,
        'n_basis_per_channel': n_basis,
    }


# ==========================================================================
# (b) 4-ELECTRON SPECTRAL DIMENSIONS
# ==========================================================================

def four_electron_spectral_dimensions(
    l_max: int,
    n_basis_per_angle: int = 5,
    homonuclear: bool = False,
    sigma_only: bool = False,
) -> Dict:
    """Compute spectral dimensions for 4-electron angular problem.

    The 4-electron angular space S^11 has 3 hyperangles (alpha1, alpha2, alpha3)
    plus 8 direction angles (theta_i, phi_i for i=1..4).

    The channel expansion integrates out the 8 direction angles into angular
    momentum labels (l_i, m_i), leaving 3 hyperangles as continuous variables.

    The spectral basis is a tensor product of Jacobi polynomials in each
    hyperangle, with basis size n_basis_per_angle^3 per channel.

    Parameters
    ----------
    l_max : int
        Maximum angular momentum per electron.
    n_basis_per_angle : int
        Number of Jacobi polynomial basis functions per hyperangle per channel.
    homonuclear : bool
        If True, imposes sum(l_i) even constraint.
    sigma_only : bool
        If True, restricts to all m_i = 0.
    """
    # Channel count
    if sigma_only:
        result = four_electron_channel_count_molecular(
            l_max, sigma_only=True, homonuclear=homonuclear
        )
        n_ch = result['sigma_channels']
    else:
        result = four_electron_channel_count_molecular(
            l_max, sigma_only=False, homonuclear=homonuclear
        )
        n_ch = result['full_channels']

    # S4 singlet [2,2] reduction: ~1/6 of channels survive
    # More precise: use Burnside's lemma dim^2/|G| = 4/24 = 1/6
    n_ch_antisym = max(1, n_ch // 6)

    # For a more precise count at small l_max, we enumerate
    n_ch_antisym_exact = _count_s4_singlet_channels(l_max, sigma_only, homonuclear)

    # Spectral basis per channel: tensor product of 3 Jacobi bases
    n_alpha_basis = n_basis_per_angle ** 3

    # Total spectral dimension
    spectral_dim_raw = n_ch * n_alpha_basis
    spectral_dim_antisym_approx = n_ch_antisym * n_alpha_basis
    spectral_dim_antisym_exact = n_ch_antisym_exact * n_alpha_basis

    # FD reference: discretize each of 3 hyperangles with n_grid points
    # and n_ch channels -> n_ch * n_grid^3
    n_grid = 50  # conservative FD grid per angle
    fd_dim = n_ch * n_grid ** 3
    fd_dim_antisym = n_ch_antisym_exact * n_grid ** 3

    return {
        'n_electrons': 4,
        'angular_dim': 11,  # S^11
        'n_hyperangles': 3,  # alpha1, alpha2, alpha3
        'n_direction_angles': 8,  # 4 * (theta, phi)
        'l_max': l_max,
        'n_channels_raw': n_ch,
        'n_channels_s4_approx': n_ch_antisym,
        'n_channels_s4_exact': n_ch_antisym_exact,
        'n_basis_per_angle': n_basis_per_angle,
        'n_alpha_basis': n_alpha_basis,
        'spectral_dim_raw': spectral_dim_raw,
        'spectral_dim_s4_approx': spectral_dim_antisym_approx,
        'spectral_dim_s4_exact': spectral_dim_antisym_exact,
        'fd_grid_per_angle': n_grid,
        'fd_dim_raw': fd_dim,
        'fd_dim_s4': fd_dim_antisym,
        'compression_fd_to_spectral': (
            fd_dim_antisym / spectral_dim_antisym_exact
            if spectral_dim_antisym_exact > 0 else float('inf')
        ),
        'homonuclear': homonuclear,
        'sigma_only': sigma_only,
    }


def _count_s4_singlet_channels(
    l_max: int,
    sigma_only: bool = False,
    homonuclear: bool = False,
) -> int:
    """Count S4 singlet [2,2] channels for 4-electron problem.

    For the singlet state, the spatial wavefunction has [2,2] symmetry
    under S4 (the permutation group of 4 electrons).

    The [2,2] irrep has dimension 2. The number of [2,2] channels
    in the tensor product of 4 angular momentum representations is:

    n_{[2,2]} = (1/|S4|) sum_{g in S4} chi_{[2,2]}(g) * Tr(D(g))

    where D(g) is the representation matrix of g on the channel space.

    For practical computation, we use the character table of S4 and
    count multiplicities.

    At l_max=0 (all l_i=0, m_i=0): Only 1 channel, and S4 acts only
    on the hyperangles. The [2,2] content depends on the hyperangle
    symmetry adaptation, which affects the alpha basis, not the channel
    count. So channel count is 1.

    For general l_max, we use the 1/6 approximation refined by
    explicit enumeration of symmetry-distinct l-configurations.
    """
    if sigma_only:
        # All m_i = 0. Channels are (l1, l2, l3, l4).
        # S4 permutes the l-labels. [2,2] content:
        # Orbits under S4 weighted by [2,2] character.
        return _count_s4_22_sigma_channels(l_max, homonuclear)
    else:
        # Full (l_i, m_i) channels with M=0.
        # Use approximate 1/6 factor.
        result = four_electron_channel_count_molecular(
            l_max, sigma_only=False, homonuclear=homonuclear
        )
        n_full = result['full_channels']
        # Better estimate: exact for l_max <= 2
        return _count_s4_22_full_channels(l_max, homonuclear)


def _count_s4_22_sigma_channels(l_max: int, homonuclear: bool) -> int:
    """Count [2,2]-symmetric sigma channels.

    For sigma only (all m_i=0), channels are tuples (l1,l2,l3,l4).
    S4 acts by permuting the l-labels.

    The [2,2] representation of S4 has character:
      e: 2,  (12): 0,  (12)(34): 2,  (123): -1,  (1234): 0

    The multiplicity of [2,2] in the permutation representation on
    tuples (l1,l2,l3,l4) is computed via:
      n = (1/24) sum_{g in S4} chi_{[2,2]}(g) * (number of tuples fixed by g)

    Conjugacy classes of S4:
      {e}: 1 element, chi=2, fixed = all tuples = (l_max+1)^4
      {(12)}: 6 elements, chi=0, fixed = tuples with l_{sigma(i)}=l_i
      {(12)(34)}: 3 elements, chi=2, fixed = tuples with l1=l2 AND l3=l4
      {(123)}: 8 elements, chi=-1, fixed = tuples with l1=l2=l3
      {(1234)}: 6 elements, chi=0

    So n = (1/24) [2*(l_max+1)^4 + 0 + 2*3*n_{l1=l2,l3=l4} + (-1)*8*n_{l1=l2=l3} + 0]
    """
    L = l_max + 1

    # (l_max+1)^4 = total tuples
    # BUT with homonuclear constraint sum(l_i) even if applicable
    if homonuclear:
        n_total = _count_even_sum_tuples(l_max, 4)
        n_12_34 = _count_paired_even(l_max)
        n_123 = _count_triple_even(l_max)
    else:
        n_total = L ** 4
        # Fixed by (12)(34): l1=l2 AND l3=l4 -> L^2 choices
        n_12_34 = L ** 2
        # Fixed by (123): l1=l2=l3 -> L^2 choices (l1=l2=l3 free, l4 free)
        n_123 = L ** 2

    n_22 = (2 * n_total + 2 * 3 * n_12_34 - 1 * 8 * n_123) // 24

    # Ensure non-negative (rounding from integer arithmetic)
    return max(1, n_22)


def _count_even_sum_tuples(l_max: int, n: int) -> int:
    """Count n-tuples of values 0..l_max with even sum."""
    L = l_max + 1
    # Number of even-sum tuples = (L^n + (L_even - L_odd)^n) / 2
    # where L_even = number of even values in 0..l_max
    # and L_odd = number of odd values
    L_even = (l_max // 2) + 1  # 0, 2, 4, ...
    L_odd = L - L_even  # 1, 3, 5, ...
    diff = L_even - L_odd  # always 0 or 1
    return (L ** n + diff ** n) // 2


def _count_paired_even(l_max: int) -> int:
    """Count tuples (l,l,l',l') with l+l+l'+l' even, i.e. always even."""
    # l1=l2=l, l3=l4=l'. Sum = 2l+2l' always even.
    return (l_max + 1) ** 2


def _count_triple_even(l_max: int) -> int:
    """Count tuples (l,l,l,l4) with 3l+l4 even."""
    L = l_max + 1
    count = 0
    for l in range(L):
        if (3 * l) % 2 == 0:
            # l4 must be even
            count += (l_max // 2) + 1
        else:
            # l4 must be odd
            count += L - ((l_max // 2) + 1)
    return count


def _count_s4_22_full_channels(l_max: int, homonuclear: bool) -> int:
    """Count [2,2]-symmetric full (l_i, m_i) channels with M=0.

    Uses the same Burnside approach but on (l_i, m_i) tuples.
    S4 permutes both l-labels AND m-labels simultaneously.

    For the full M=0 space, we use the Burnside formula:
    n = (1/24) sum_{g in S4} chi_{[2,2]}(g) * |Fix(g)|

    where Fix(g) = number of M=0 tuples fixed by permutation g.
    """
    # Conjugacy classes of S4 with element counts and [2,2] characters
    # {e}: 1 elem, chi=2
    # {(ij)}: 6 elem, chi=0
    # {(ij)(kl)}: 3 elem, chi=2
    # {(ijk)}: 8 elem, chi=-1
    # {(ijkl)}: 6 elem, chi=0

    # Fix(e) = total M=0 count
    fix_e = 0
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            for l3 in range(l_max + 1):
                for l4 in range(l_max + 1):
                    if homonuclear and (l1 + l2 + l3 + l4) % 2 != 0:
                        continue
                    fix_e += _count_m_zero_states(l1, l2, l3, l4)

    # Fix((12)(34)): tuples with (l1,m1)=(l2,m2) AND (l3,m3)=(l4,m4)
    # AND m1+m2+m3+m4 = 2*m1+2*m3 = 0 so m3 = -m1
    fix_1234_pair = 0
    for l1 in range(l_max + 1):
        for l3 in range(l_max + 1):
            if homonuclear and (2 * l1 + 2 * l3) % 2 != 0:
                continue
            # m1 in [-l1,l1], m3=-m1 requires |m1| <= l3
            m_range = min(l1, l3)
            fix_1234_pair += 2 * m_range + 1

    # Fix((123)): tuples with (l1,m1)=(l2,m2)=(l3,m3), l4 free
    # M=0: 3*m1+m4=0 so m4=-3*m1
    fix_123 = 0
    for l1 in range(l_max + 1):
        for l4 in range(l_max + 1):
            if homonuclear and (3 * l1 + l4) % 2 != 0:
                continue
            for m1 in range(-l1, l1 + 1):
                m4 = -3 * m1
                if abs(m4) <= l4:
                    fix_123 += 1

    # Burnside formula
    n_22 = (2 * fix_e + 0 + 2 * 3 * fix_1234_pair - 1 * 8 * fix_123 + 0) // 24
    return max(1, n_22)


# ==========================================================================
# (c) SO(12) FREE SPECTRUM WITH CHANNEL DECOMPOSITION
# ==========================================================================

def so12_free_spectrum_by_channel(
    l_max: int,
    n_basis_per_angle: int = 5,
    sigma_only: bool = True,
) -> Dict:
    """Build the free SO(12) spectrum decomposed by channel.

    At Level 4, each channel (l1, l2) contributes Jacobi polynomial
    modes k=0,1,2,...,n_basis-1 with free eigenvalue:
      mu_free = nu*(nu+4)/2 where nu = l1 + l2 + 2k

    At 4 electrons, each channel (l1, l2, l3, l4) (sigma only)
    contributes tensor-product modes (k1, k2, k3) where:
      - k1 indexes Jacobi basis for alpha1 (pair 1-2 hyperangle)
      - k2 indexes Jacobi basis for alpha2 (inter-pair hyperangle)
      - k3 indexes Jacobi basis for alpha3 (pair 3-4 hyperangle)

    The free eigenvalue generalizes the SO(6) pattern. For the Jacobi
    tree, the exact free eigenvalue per mode depends on the tree structure.

    For the democratic/balanced tree:
      Pair (1,2): internal quantum number n12 = l1 + l2 + 2*k1
      Pair (3,4): internal quantum number n34 = l3 + l4 + 2*k3
      Inter-pair: grand quantum number nu = n12 + n34 + 2*k2
                  (equivalently, via SO(6) x SO(6) -> SO(12) coupling)
      mu_free = nu*(nu+10)/2  [SO(12) Casimir]

    BUT this is the TOTAL Casimir. The individual mode contributes
    a specific nu = (l1+l2+2k1) + (l3+l4+2k3) + 2*k2.
    """
    if not sigma_only:
        raise NotImplementedError("Full (l,m) channels not yet implemented")

    channels = []
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            for l3 in range(l_max + 1):
                for l4 in range(l_max + 1):
                    channels.append((l1, l2, l3, l4))

    modes = []
    nb = n_basis_per_angle
    for ch in channels:
        l1, l2, l3, l4 = ch
        for k1 in range(nb):
            for k2 in range(nb):
                for k3 in range(nb):
                    # Pair internal quantum numbers
                    n12 = l1 + l2 + 2 * k1
                    n34 = l3 + l4 + 2 * k3
                    # Grand quantum number
                    nu = n12 + n34 + 2 * k2
                    mu_free = nu * (nu + 10) / 2.0
                    modes.append({
                        'channel': ch,
                        'k': (k1, k2, k3),
                        'n12': n12,
                        'n34': n34,
                        'nu': nu,
                        'mu_free': mu_free,
                    })

    # Sort by energy
    modes.sort(key=lambda x: x['mu_free'])

    # Count unique nu values and their multiplicities
    nu_values = {}
    for m in modes:
        nu = m['nu']
        if nu not in nu_values:
            nu_values[nu] = 0
        nu_values[nu] += 1

    return {
        'l_max': l_max,
        'n_basis_per_angle': nb,
        'n_channels': len(channels),
        'n_modes_total': len(modes),
        'n_modes_per_channel': nb ** 3,
        'modes_sorted': modes,
        'nu_multiplicities': dict(sorted(nu_values.items())),
    }


# ==========================================================================
# (d) SPECTRAL BASIS STRUCTURE: Jacobi tree for 3 hyperangles
# ==========================================================================

def jacobi_tree_basis_structure(
    l_max: int,
    n_basis_per_angle: int = 5,
) -> Dict:
    """Describe the tensor-product Jacobi polynomial basis for 3 hyperangles.

    The democratic Jacobi tree for 4 electrons:
      alpha1 = arctan(r2/r1)        -> Jacobi P_k1^{(l2+1/2, l1+1/2)}(cos 2*alpha1)
      alpha3 = arctan(r4/r3)        -> Jacobi P_k3^{(l4+1/2, l3+1/2)}(cos 2*alpha3)
      alpha2 = arctan(rho34/rho12)  -> Jacobi P_k2^{(a2, b2)}(cos 2*alpha2)

    The alpha1 and alpha3 bases are IDENTICAL to the Level 4 single-alpha
    basis (each pair of electrons sees the same hyperangle structure).

    The alpha2 basis couples the two pairs. Its Jacobi parameters are:
      a2 = n34 + 2 = (l3+l4+2k3) + 2
      b2 = n12 + 2 = (l1+l2+2k1) + 2
    These depend on the pair quantum numbers, creating coupling between
    the three hyperangle bases.

    IMPORTANT: The tensor product P_k1 * P_k2 * P_k3 is NOT a simple
    direct product because the alpha2 Jacobi parameters depend on k1
    and k3. This means the alpha2 basis must be rebuilt for each
    (channel, k1, k3) combination.

    For practical computation, we can:
    (a) Fix (k1, k3) and solve the alpha2 problem -> nested 1D solves
    (b) Use a fixed n_basis_per_angle for alpha2 regardless of (k1, k3)
    (c) Truncate the tensor product by total nu cutoff
    """
    nb = n_basis_per_angle

    # For l_max=0, all l_i=0, so Jacobi parameters are:
    # alpha1: P_k1^{(1/2, 1/2)} = Chebyshev-like
    # alpha3: P_k3^{(1/2, 1/2)} = same
    # alpha2: P_k2^{(2k3+2, 2k1+2)} <- depends on k1, k3!

    # This means the basis is NOT a tensor product in the standard sense.
    # However, we can still count the total dimension as n_ch * nb^3
    # because we use nb basis functions for EACH hyperangle regardless
    # of the Jacobi parameters.

    info = {
        'alpha1_basis': 'Jacobi P_k1^{(l2+1/2, l1+1/2)}(cos 2*alpha1)',
        'alpha3_basis': 'Jacobi P_k3^{(l4+1/2, l3+1/2)}(cos 2*alpha3)',
        'alpha2_basis': 'Jacobi P_k2^{(n34+2, n12+2)}(cos 2*alpha2)',
        'alpha2_depends_on_k1_k3': True,
        'alpha2_note': (
            'The alpha2 Jacobi parameters depend on the pair quantum '
            'numbers n12=l1+l2+2k1 and n34=l3+l4+2k3. This creates '
            'non-trivial coupling between the three hyperangle bases. '
            'In practice, we use a fixed n_basis for alpha2 at each '
            '(channel, k1, k3) point, giving nb^3 modes per channel.'
        ),
        'tensor_product_dim': nb ** 3,
        'basis_per_channel': nb ** 3,
    }

    return info


# ==========================================================================
# (e) COMPRESSION RATIO ANALYSIS
# ==========================================================================

def compression_analysis(
    l_max_values: List[int] = None,
    n_basis_per_angle: int = 5,
    n_fd_per_angle: int = 50,
) -> List[Dict]:
    """Compute compression ratios at multiple l_max values.

    Compares:
    - FD dimension: n_ch_antisym * n_fd^3  (3D finite-difference grid)
    - Spectral dimension: n_ch_antisym * n_basis^3  (tensor product Jacobi)
    - Level 4 reference: n_ch_l4 * n_basis_l4

    Parameters
    ----------
    l_max_values : list of int
        Angular momentum cutoffs to analyze.
    n_basis_per_angle : int
        Spectral basis size per hyperangle.
    n_fd_per_angle : int
        FD grid points per hyperangle.

    Returns
    -------
    List of dicts, one per l_max.
    """
    if l_max_values is None:
        l_max_values = [0, 1, 2, 3]

    results = []
    for l_max in l_max_values:
        # 4-electron dimensions
        dims = four_electron_spectral_dimensions(
            l_max,
            n_basis_per_angle=n_basis_per_angle,
            homonuclear=False,  # LiH
            sigma_only=True,
        )

        # Level 4 reference (2-electron, homonuclear H2)
        l4 = level4_spectral_dimensions(l_max, n_basis=10, homonuclear=True)

        # Full (non-sigma) 4-electron
        dims_full = four_electron_spectral_dimensions(
            l_max,
            n_basis_per_angle=n_basis_per_angle,
            homonuclear=False,
            sigma_only=False,
        )

        results.append({
            'l_max': l_max,
            # Level 4 reference
            'l4_channels': l4['n_channels'],
            'l4_spectral_dim': l4['spectral_dim'],
            'l4_fd_dim': l4['fd_dim'],
            'l4_compression': l4['compression'],
            # 4-electron sigma
            '4e_channels_raw': dims['n_channels_raw'],
            '4e_channels_s4': dims['n_channels_s4_exact'],
            '4e_spectral_dim': dims['spectral_dim_s4_exact'],
            '4e_fd_dim': dims['fd_dim_s4'],
            '4e_compression': dims['compression_fd_to_spectral'],
            # 4-electron full
            '4e_full_channels_raw': dims_full['n_channels_raw'],
            '4e_full_channels_s4': dims_full['n_channels_s4_exact'],
            '4e_full_spectral_dim': dims_full['spectral_dim_s4_exact'],
            # Ratio to Level 4
            'ratio_to_l4': (
                dims['spectral_dim_s4_exact'] / l4['spectral_dim']
                if l4['spectral_dim'] > 0 else float('inf')
            ),
        })

    return results


# ==========================================================================
# (f) COUPLING MATRIX CLASSIFICATION
# ==========================================================================

def coupling_matrix_classification() -> List[Dict]:
    """Classify all coupling matrix elements for 4-electron angular problem.

    For each potential term, determine:
    - Count: how many terms of this type
    - rho-dependent: whether it depends on the parameter rho = R/(2*R_e)
    - Algebraic: whether it can be precomputed via Gaunt integrals
    - Structure: what mathematical structure the coupling has
    """
    terms = []

    # 1. Free Hamiltonian (Casimir)
    terms.append({
        'term': 'H_free (SO(12) Casimir)',
        'count': 1,
        'rho_dependent': False,
        'algebraic': True,
        'structure': (
            'Diagonal. mu_free(nu) = nu*(nu+10)/2. '
            'Exact from quantum numbers, no computation needed.'
        ),
    })

    # 2. Intra-pair V_ee: pairs (1,2) and (3,4)
    terms.append({
        'term': 'V_ee intra-pair',
        'count': 2,
        'rho_dependent': False,
        'algebraic': True,
        'structure': (
            'Pairs (1,2) and (3,4). Each has IDENTICAL structure to '
            'Level 4 V_ee: 1/r_{12} expanded in Legendre multipoles, '
            'coupling (l_i, l_j) -> (l_i\', l_j\') via Gaunt integrals. '
            'Acts on alpha1 or alpha3 respectively, diagonal in alpha2. '
            'Precomputed once via Gaunt integrals (rho-independent).'
        ),
    })

    # 3. Cross-pair V_ee: pairs (1,3), (1,4), (2,3), (2,4)
    terms.append({
        'term': 'V_ee cross-pair',
        'count': 4,
        'rho_dependent': False,
        'algebraic': 'PARTIAL',
        'structure': (
            'Pairs (1,3), (1,4), (2,3), (2,4). Each 1/r_{ij} couples '
            'electrons from DIFFERENT pairs in the Jacobi tree. '
            'The multipole expansion still uses Legendre polynomials, '
            'and Gaunt integrals still apply for the angular part. '
            'BUT the radial part (min/max of hyperradial distances) '
            'involves ALL THREE hyperangles simultaneously. '
            'The coupling operator for cross-pair (i,j) where i is in '
            'pair (1,2) and j is in pair (3,4) has the form:\n'
            '  1/r_{ij} = 1/(R_e * d_{ij}(alpha1, alpha2, alpha3, angles))\n'
            'where d_{ij} depends on all 3 hyperangles. '
            'This is NOT separable into alpha1 * alpha2 * alpha3 factors. '
            'The Gaunt integrals handle the direction-angle part, but '
            'the hyperangle part requires 3D quadrature. '
            'Precomputed once (rho-independent) but requires O(n_quad^3) '
            'quadrature points instead of O(n_quad) for intra-pair.'
        ),
    })

    # 4. Nuclear attraction: 4 electrons x 2 nuclei = 8 terms
    terms.append({
        'term': 'V_nuc (nuclear attraction)',
        'count': 8,
        'rho_dependent': True,
        'algebraic': True,
        'structure': (
            '4 electrons x 2 nuclei. Each -Z_X/|r_i - R_X| has the '
            'split-region Legendre expansion (same structure as Level 4). '
            'Nuclear coupling is diagonal in the angular momenta of '
            'ALL electrons EXCEPT electron i (selection rule). '
            'The rho-dependence enters through the split-region boundary '
            'at s_i = rho. Acts on ONE hyperangle (the one containing '
            'electron i) and is diagonal in the other two. '
            'Recomputed at each rho-point. Gaunt integrals for angular part.'
        ),
    })

    # 5. Centrifugal barrier
    terms.append({
        'term': 'Centrifugal barrier',
        'count': 1,
        'rho_dependent': False,
        'algebraic': True,
        'structure': (
            'Absorbed into the Jacobi basis envelope functions: '
            '(cos alpha_i)^{a_i+1} (sin alpha_i)^{b_i+1} for each '
            'hyperangle. The remaining centrifugal potential is '
            '(3N-1)(3N-3)/(8*R_e^2) = 99/8*R_e^2 (constant). '
            'Fully algebraic, diagonal.'
        ),
    })

    return terms


# ==========================================================================
# (g) CROSS-PAIR V_EE ANALYSIS
# ==========================================================================

def cross_pair_vee_analysis() -> Dict:
    """Detailed analysis of cross-pair electron-electron coupling.

    The cross-pair V_ee is the structurally new element at 4 electrons.
    This function analyzes whether it can be computed algebraically.
    """
    return {
        'problem': (
            'For electrons i (in pair 1-2) and j (in pair 3-4), the '
            'distance r_{ij} involves coordinates from BOTH pairs:\n'
            '  r_i = R_e * sin(alpha2) * rho_12_hat * direction_i\n'
            '  r_j = R_e * cos(alpha2) * rho_34_hat * direction_j\n'
            'where rho_12_hat involves alpha1, rho_34_hat involves alpha3.\n\n'
            'Explicitly:\n'
            '  r_{13}^2 = R_e^2 [sin^2(alpha2)*cos^2(alpha1) + '
            'cos^2(alpha2)*cos^2(alpha3) - 2*sin(alpha2)*cos(alpha1)*'
            'cos(alpha2)*cos(alpha3)*cos(theta_{13})]\n\n'
            'This depends on ALL THREE hyperangles AND the inter-electron '
            'angle theta_{13}.'
        ),
        'multipole_expansion': (
            'The standard multipole expansion of 1/r_{13} is:\n'
            '  1/r_{13} = sum_k (r_</r_>)^k P_k(cos theta_{13}) / r_>\n\n'
            'where r_< = min(r_1, r_3), r_> = max(r_1, r_3).\n\n'
            'In hyperspherical coords:\n'
            '  r_1 = R_e * sin(alpha2) * cos(alpha1)\n'
            '  r_3 = R_e * cos(alpha2) * cos(alpha3)\n\n'
            'The ratio r_1/r_3 = sin(alpha2)*cos(alpha1) / '
            '(cos(alpha2)*cos(alpha3)) depends on ALL THREE hyperangles.'
        ),
        'gaunt_structure': (
            'YES, the angular (direction) part still uses Gaunt integrals:\n'
            '  Integral of Y_{l1,m1}^* Y_{k,q} Y_{l3,m3} d(omega) = '
            'Gaunt(l1,k,l3; m1,q,m3)\n\n'
            'This gives the usual selection rules: |l1-l3| <= k <= l1+l3, '
            'm1-m3 = q.\n\n'
            'The radial/hyperangular part is the 3D integral:\n'
            '  I = int d(alpha1) d(alpha2) d(alpha3) '
            'phi_i(alpha1,alpha2,alpha3) * (r_</r_>)^k / r_> * '
            'phi_j(alpha1,alpha2,alpha3) * volume_element\n\n'
            'This is a genuine 3D integral over the three hyperangles, '
            'NOT separable into three 1D integrals.'
        ),
        'separability': (
            'PARTIAL separation is possible:\n'
            '1. The split-region boundary r_1 = r_3 defines a 2D surface '
            'in (alpha1, alpha2, alpha3) space.\n'
            '2. For fixed alpha2, the ratio r_1/r_3 = '
            'tan(alpha2)*cos(alpha1)/cos(alpha3) separates into '
            'alpha1 and alpha3 dependence.\n'
            '3. The split region: r_1 < r_3 when '
            'tan(alpha2) < cos(alpha3)/cos(alpha1).\n\n'
            'This means the integral can be written as a 1D integral '
            'in alpha2, with the integrand being a PRODUCT of two 1D '
            'integrals in alpha1 and alpha3 (for each region). '
            'The split boundary in (alpha1, alpha3) depends on alpha2, '
            'so the inner integrals change at each alpha2 quadrature point.\n\n'
            'Conclusion: REDUCIBLE from 3D to nested 1D integrals, '
            'but not fully algebraic. Cost: O(n_quad * n_quad * n_quad) '
            'where n_quad is 1D quadrature size. Precomputed once.'
        ),
        'algebraic_verdict': (
            'Cross-pair V_ee has Gaunt integral structure for the ANGULAR '
            'part (directions), but the HYPERANGULAR part requires 3D '
            'numerical integration (reducible to nested 1D). This is '
            'analogous to the Level 4 split-region Legendre expansion '
            'but in 3D instead of 1D. It is rho-INDEPENDENT and can be '
            'precomputed once, but the precomputation cost scales as '
            'O(n_ch^2 * n_basis^3 * n_quad^3) instead of '
            'O(n_ch^2 * n_basis * n_quad) at Level 4.'
        ),
    }


# ==========================================================================
# (h) TRACTABILITY VERDICT
# ==========================================================================

def tractability_verdict(
    l_max: int = 2,
    n_basis_per_angle: int = 5,
) -> Dict:
    """Assess tractability of spectral 4-electron angular solver.

    Computes the angular matrix dimension and estimates wall time
    for eigensolves at representative rho-points.
    """
    dims = four_electron_spectral_dimensions(
        l_max,
        n_basis_per_angle=n_basis_per_angle,
        homonuclear=False,
        sigma_only=True,
    )

    spectral_dim = dims['spectral_dim_s4_exact']

    # Wall time estimate for dense eigensolve O(dim^3)
    # Assume ~1 GFLOP/s effective for LAPACK dsyev
    gflops = 1e9
    cost_per_rho = spectral_dim ** 3 / gflops  # seconds
    n_rho = 130  # standard grid
    total_time = cost_per_rho * n_rho

    # Cross-pair V_ee precomputation cost
    n_ch = dims['n_channels_s4_exact']
    nb3 = n_basis_per_angle ** 3
    n_quad = 50  # per hyperangle
    vee_precomp = n_ch ** 2 * nb3 ** 2 * n_quad ** 3 / gflops

    # Level 4 reference
    l4 = level4_spectral_dimensions(l_max, n_basis=10, homonuclear=True)
    l4_time = l4['spectral_dim'] ** 3 / gflops * n_rho

    tractable = total_time < 3600  # less than 1 hour
    practical = total_time < 60  # less than 1 minute

    return {
        'l_max': l_max,
        'spectral_dim': spectral_dim,
        'cost_per_rho_s': cost_per_rho,
        'total_time_s': total_time,
        'total_time_hours': total_time / 3600,
        'vee_precomp_s': vee_precomp,
        'l4_spectral_dim': l4['spectral_dim'],
        'l4_total_time_s': l4_time,
        'ratio_to_l4_dim': spectral_dim / l4['spectral_dim'] if l4['spectral_dim'] > 0 else float('inf'),
        'ratio_to_l4_time': total_time / l4_time if l4_time > 0 else float('inf'),
        'tractable_1hr': tractable,
        'practical_1min': practical,
        'verdict': (
            'PRACTICAL' if practical else
            'TRACTABLE' if tractable else
            'INTRACTABLE'
        ),
    }


# ==========================================================================
# (i) FULL ANALYSIS RUNNER
# ==========================================================================

def run_full_analysis(
    n_basis_per_angle: int = 5,
    output_file: Optional[str] = None,
) -> str:
    """Run the complete Track AK spectral compression analysis.

    Returns a formatted report string.
    """
    lines = []
    lines.append("=" * 72)
    lines.append("TRACK AK: 4-Electron Spectral Compression Analysis")
    lines.append("=" * 72)

    # 1. Compression table
    lines.append("\n1. COMPRESSION TABLE")
    lines.append("-" * 72)
    results = compression_analysis(
        l_max_values=[0, 1, 2, 3],
        n_basis_per_angle=n_basis_per_angle,
    )

    header = (
        f"{'l_max':>5} | {'Full dim':>8} | {'S4-reduced':>10} | "
        f"{'Spectral':>10} | {'Ratio':>8} | {'L4 ref':>8} | {'Verdict':>12}"
    )
    lines.append(header)
    lines.append("-" * len(header))

    for r in results:
        v = tractability_verdict(r['l_max'], n_basis_per_angle)
        lines.append(
            f"{r['l_max']:>5} | "
            f"{r['4e_channels_raw']:>8} | "
            f"{r['4e_channels_s4']:>10} | "
            f"{r['4e_spectral_dim']:>10} | "
            f"{r['4e_compression']:>8.1f} | "
            f"{r['l4_spectral_dim']:>8} | "
            f"{v['verdict']:>12}"
        )

    # 2. Coupling matrix classification
    lines.append("\n\n2. COUPLING MATRIX CLASSIFICATION")
    lines.append("-" * 72)
    terms = coupling_matrix_classification()
    lines.append(
        f"{'Term':>30} | {'Count':>5} | {'rho-dep?':>8} | "
        f"{'Algebraic?':>10}"
    )
    lines.append("-" * 72)
    for t in terms:
        lines.append(
            f"{t['term']:>30} | {t['count']:>5} | "
            f"{'Yes' if t['rho_dependent'] else 'No':>8} | "
            f"{t['algebraic']!s:>10}"
        )

    # 3. Cross-pair V_ee
    lines.append("\n\n3. CROSS-PAIR V_EE ANALYSIS")
    lines.append("-" * 72)
    xp = cross_pair_vee_analysis()
    lines.append(f"Algebraic verdict: {xp['algebraic_verdict']}")

    # 4. Tractability at l_max=2
    lines.append("\n\n4. TRACTABILITY AT l_max=2")
    lines.append("-" * 72)
    v2 = tractability_verdict(2, n_basis_per_angle)
    lines.append(f"  Spectral dim: {v2['spectral_dim']}")
    lines.append(f"  Cost per rho: {v2['cost_per_rho_s']:.2f} s")
    lines.append(f"  Total (130 rho): {v2['total_time_s']:.1f} s ({v2['total_time_hours']:.2f} hrs)")
    lines.append(f"  V_ee precomp: {v2['vee_precomp_s']:.1f} s")
    lines.append(f"  Level 4 ref dim: {v2['l4_spectral_dim']}")
    lines.append(f"  Level 4 ref time: {v2['l4_total_time_s']:.4f} s")
    lines.append(f"  Ratio to L4 (dim): {v2['ratio_to_l4_dim']:.1f}x")
    lines.append(f"  Ratio to L4 (time): {v2['ratio_to_l4_time']:.1f}x")
    lines.append(f"  VERDICT: {v2['verdict']}")

    # 5. Recommended n_basis
    lines.append("\n\n5. RECOMMENDED SPECTRAL BASIS")
    lines.append("-" * 72)
    for nb in [3, 4, 5, 6, 7]:
        v = tractability_verdict(2, nb)
        lines.append(
            f"  n_basis_per_angle={nb}: dim={v['spectral_dim']}, "
            f"time={v['total_time_s']:.1f}s, {v['verdict']}"
        )

    # 6. SO(12) spectrum at l_max=0
    lines.append("\n\n6. SO(12) FREE SPECTRUM (l_max=0)")
    lines.append("-" * 72)
    spec = so12_free_spectrum_by_channel(0, n_basis_per_angle=3)
    nu_mult = spec['nu_multiplicities']
    for nu, mult in list(nu_mult.items())[:10]:
        mu = nu * (nu + 10) / 2.0
        so12_degen = _so_d_degeneracy(12, nu)
        lines.append(
            f"  nu={nu}: mu_free={mu:.1f}, "
            f"SO(12) degeneracy={so12_degen}, "
            f"spectral multiplicity={mult}"
        )

    report = "\n".join(lines)

    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)

    return report
