"""
Scoping module for N-electron mol-frame hyperspherical coordinates.

Track AG: Feasibility analysis for full 4-electron LiH solver
without PK, Z_eff, or composed geometry.

This module computes dimension counts, symmetry reductions,
and cost estimates for N-electron hyperspherical solvers.
It does NOT implement a solver -- it is a scoping tool.
"""

import numpy as np
from math import factorial, comb
from typing import Tuple, List, Dict, Optional
from itertools import product as itertools_product


# ==========================================================================
# (a) COORDINATE ANALYSIS
# ==========================================================================

def n_electron_dimensions(N: int) -> Dict[str, int]:
    """Compute coordinate dimensions for N electrons + fixed nuclei.

    Parameters
    ----------
    N : int
        Number of electrons.

    Returns
    -------
    dict with keys:
        config_dim : total configuration space dimension (3N)
        hyperradius : 1 (always)
        angular_dim : 3N - 1 (angular manifold on S^{3N-1})
        sphere : "S^{3N-1}"
        isometry_group : "SO(3N)"
        casimir_formula : string
    """
    d = 3 * N
    return {
        'N': N,
        'config_dim': d,
        'hyperradius': 1,
        'angular_dim': d - 1,
        'sphere': f'S^{d - 1}',
        'isometry_group': f'SO({d})',
        'casimir_formula': f'nu*(nu+{d - 2})/2',
        'jacobian_centrifugal': f'{(d-1)*(d-3)}/8 = {(d-1)*(d-3)/8}',
    }


def casimir_eigenvalue(N: int, nu: int) -> float:
    """SO(3N) Casimir eigenvalue for grand angular momentum nu.

    mu_free(nu) = nu * (nu + 3N - 2) / 2

    Parameters
    ----------
    N : int
        Number of electrons.
    nu : int
        Grand angular momentum quantum number.
    """
    d = 3 * N
    return nu * (nu + d - 2) / 2.0


def effective_angular_momentum(N: int, nu: int) -> float:
    """Effective angular momentum l_eff from Casimir eigenvalue.

    After extracting R^{-(3N-1)/2}, the centrifugal term is
    l_eff(l_eff+1)/(2R^2) where l_eff = nu + (3N-4)/2.

    Parameters
    ----------
    N : int
        Number of electrons.
    nu : int
        Grand angular momentum quantum number.
    """
    return nu + (3 * N - 4) / 2.0


# ==========================================================================
# (b) ANGULAR HILBERT SPACE: Channel counting
# ==========================================================================

def level4_channel_count(l_max: int, sigma_only: bool = True,
                         homonuclear: bool = True) -> int:
    """Count channels for 2-electron Level 4 solver.

    Channels: (l1, m1, l2, m2) with m1 + m2 = 0, l1+l2 even (homonuclear).
    """
    count = 0
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            if homonuclear and (l1 + l2) % 2 != 0:
                continue
            if sigma_only:
                count += 1  # m1 = m2 = 0 only
            else:
                m_limit = min(l1, l2)
                count += 2 * m_limit + 1  # m1 = -m_limit..m_limit
    return count


def four_electron_channel_count_atomic(l_max: int) -> Dict:
    """Count angular channels for 4-electron atomic (1-center) problem.

    For 4 electrons at a single center, each electron has angular momentum
    (l_i, m_i). The total L = 0, M = 0 state requires:
      sum of all m_i = 0
      Total L from coupling: L = 0

    The angular coordinates are:
      - 3 hyperangles (alpha1, alpha2, alpha3) parameterizing the
        four radial magnitudes via the hyperspherical tree
      - 4 x 2 = 8 angles (theta_i, phi_i) for electron directions
      Total: 11 angular coordinates on S^11

    Channels are labeled by (l1, l2, l3, l4) with coupled total L=0.
    For each set of l-values, the number of M_L=0 couplings to L=0
    is determined by the Clebsch-Gordan series.

    This function counts the UNCOUPLED channels (l1, l2, l3, l4, m1, m2, m3, m4)
    with m1+m2+m3+m4 = 0 as an upper bound.
    """
    # Uncoupled count: all (l1,l2,l3,l4) with each l_i <= l_max,
    # all (m1,m2,m3,m4) with m_i in [-l_i, l_i], sum = 0
    uncoupled = 0
    coupled_l0 = 0  # L=0 coupled channels

    l_sets = []
    for l1 in range(l_max + 1):
        for l2 in range(l_max + 1):
            for l3 in range(l_max + 1):
                for l4 in range(l_max + 1):
                    l_sets.append((l1, l2, l3, l4))

    # For each l-set, count M=0 uncoupled states
    for (l1, l2, l3, l4) in l_sets:
        m_count = _count_m_zero_states(l1, l2, l3, l4)
        uncoupled += m_count

        # Count L=0 coupled channels from this l-set
        # Couple l1+l2 -> L12, l3+l4 -> L34, then L12+L34 -> L=0
        # Requires L12 = L34, so L=0 count = sum over L12 of
        # min(l1,l2,L12) * min(l3,l4,L34) with L34=L12
        n_l0 = _count_L0_couplings(l1, l2, l3, l4)
        coupled_l0 += n_l0

    return {
        'l_max': l_max,
        'n_l_sets': len(l_sets),
        'uncoupled_M0': uncoupled,
        'coupled_L0': coupled_l0,
    }


def _count_m_zero_states(l1: int, l2: int, l3: int, l4: int) -> int:
    """Count states with m1+m2+m3+m4 = 0."""
    count = 0
    for m1 in range(-l1, l1 + 1):
        for m2 in range(-l2, l2 + 1):
            for m3 in range(-l3, l3 + 1):
                m4 = -(m1 + m2 + m3)
                if abs(m4) <= l4:
                    count += 1
    return count


def _count_L0_couplings(l1: int, l2: int, l3: int, l4: int) -> int:
    """Count number of L=0 coupled angular channels.

    Couple (l1,l2) -> L12, (l3,l4) -> L34, require L12 = L34 for L=0.
    Each L12 value from |l1-l2| to l1+l2 can couple with L34=L12
    from |l3-l4| to l3+l4 if L12 is in range for both.
    """
    L12_min = abs(l1 - l2)
    L12_max = l1 + l2
    L34_min = abs(l3 - l4)
    L34_max = l3 + l4
    # L=0 requires L12 = L34
    L_min = max(L12_min, L34_min)
    L_max = min(L12_max, L34_max)
    if L_min > L_max:
        return 0
    return L_max - L_min + 1


def four_electron_channel_count_molecular(
    l_max: int,
    sigma_only: bool = False,
    homonuclear: bool = False,  # LiH is heteronuclear
) -> Dict:
    """Count angular channels for 4-electron mol-frame hyperspherical.

    For a molecule with axial symmetry (C_inf_v for LiH), each electron
    has (l_i, m_i) with the constraint M_total = sum(m_i) = 0.

    For Sigma state: M_total = 0.

    Gerade constraint for homonuclear: sum(l_i) even.
    No such constraint for heteronuclear (LiH).

    Channels: (l1, m1, l2, m2, l3, m3, l4, m4) with sum(m_i) = 0.

    With sigma_only (all m_i = 0), this reduces to (l1, l2, l3, l4).
    """
    if sigma_only:
        count = 0
        for l1 in range(l_max + 1):
            for l2 in range(l_max + 1):
                for l3 in range(l_max + 1):
                    for l4 in range(l_max + 1):
                        if homonuclear and sum([l1, l2, l3, l4]) % 2 != 0:
                            continue
                        count += 1
        return {
            'l_max': l_max,
            'sigma_channels': count,
            'homonuclear': homonuclear,
        }
    else:
        # Full m-coupled count
        count = 0
        for l1 in range(l_max + 1):
            for l2 in range(l_max + 1):
                for l3 in range(l_max + 1):
                    for l4 in range(l_max + 1):
                        if homonuclear and sum([l1, l2, l3, l4]) % 2 != 0:
                            continue
                        m_count = _count_m_zero_states(l1, l2, l3, l4)
                        count += m_count
        return {
            'l_max': l_max,
            'full_channels': count,
            'homonuclear': homonuclear,
        }


# ==========================================================================
# (c) SYMMETRY REDUCTION
# ==========================================================================

def s4_irrep_dimensions() -> Dict[str, int]:
    """Irreducible representations of S_4 (permutation group of 4 objects).

    Young diagrams:
      [4]     -> fully symmetric       -> dim 1
      [3,1]   -> standard              -> dim 3
      [2,2]   -> two-row               -> dim 2
      [2,1,1] -> standard conjugate    -> dim 3
      [1,1,1,1] -> fully antisymmetric -> dim 1
    """
    return {
        '[4]': 1,       # fully symmetric
        '[3,1]': 3,     # standard
        '[2,2]': 2,     # two-row
        '[2,1,1]': 3,   # standard conjugate
        '[1,1,1,1]': 1, # fully antisymmetric (sign representation)
    }


def spin_spatial_pairing() -> Dict[str, Dict]:
    """Map total spin S to required spatial S_4 irrep.

    For 4 electrons with total spin S, the spin state transforms
    under the conjugate Young diagram. The spatial part must
    transform under the original diagram to make the total
    wavefunction antisymmetric.

    Spin Young diagrams (columns = spin-down assignments):
      S=2: [4]         -> spatial must be [1,1,1,1] (antisym)  -- quintet
      S=1: [3,1]       -> spatial must be [2,1,1]              -- triplet
      S=0: [2,2]       -> spatial must be [2,2]                -- singlet (ground state of LiH)
      S=0: [2,2] is the ONLY S=0 irrep for 4 electrons

    LiH ground state: 1Sigma+ (S=0), so spatial irrep is [2,2].
    """
    return {
        'S=0 (singlet)': {
            'spin_irrep': '[2,2]',
            'spatial_irrep': '[2,2]',
            'spatial_dim': 2,
            'description': 'Spatial wavefunction has [2,2] symmetry under S4',
        },
        'S=1 (triplet)': {
            'spin_irrep': '[3,1]',
            'spatial_irrep': '[2,1,1]',
            'spatial_dim': 3,
        },
        'S=2 (quintet)': {
            'spin_irrep': '[4]',
            'spatial_irrep': '[1,1,1,1]',
            'spatial_dim': 1,
        },
    }


def so3_reduction_factor(l_max: int, N: int = 4) -> float:
    """Estimate SO(3) spatial rotation symmetry reduction.

    SO(3) allows restricting to L=0 total angular momentum.
    For N electrons each with l up to l_max, the total number
    of uncoupled (l,m) states is sum_{l=0}^{l_max} (2l+1) = (l_max+1)^2
    per electron.

    Total uncoupled space: ((l_max+1)^2)^N
    L=0, M=0 subspace: roughly ((l_max+1)^2)^N / (2*0+1) per L value,
    but more precisely it's the number of L=0 irreps in the tensor product.

    For a rough estimate: L=0 subspace is ~1/(l_max+1)^2 of the total
    (removing 2 angular degrees of freedom from the 3N-1 = 11 dimensional space).
    """
    total = ((l_max + 1) ** 2) ** N
    # The L=0 restriction removes 3 degrees of freedom (Euler angles)
    # from the 11-dimensional angular space, but the channel count
    # reduction is more subtle.
    return total


def s4_reduction_factor_singlet() -> str:
    """Describe S_4 [2,2] symmetry reduction for singlet.

    The [2,2] irrep has dimension 2, out of 4! = 24 total S_4 elements.
    The projection onto [2,2] retains dim([2,2])^2 / |S_4| = 4/24 = 1/6
    of the full (unsymmetrized) channel space.

    More precisely: the S_4 decomposition of the angular tensor product
    gives a multiplicity of L=0 [2,2] channels that is roughly 1/6 of
    the total uncoupled L=0 M=0 count. The factor of 2 (irrep dimension)
    means we need 2 copies of each spatial function.
    """
    return (
        "S_4 [2,2] projection: spatial irrep dimension 2, |S_4| = 24. "
        "Burnside fraction: dim^2/|G| = 4/24 = 1/6. "
        "Approximately 1/6 of uncoupled L=0 channels survive antisymmetrization."
    )


# ==========================================================================
# (d) SEPARATION STRUCTURE
# ==========================================================================

def separation_analysis_4e() -> Dict[str, str]:
    """Analyze parameter dependence for 4-electron mol-frame problem.

    At Level 4 (2e), the angular problem depends on (R, R_e) only through
    rho = R/(2*R_e), collapsing a 2D sweep to 1D.

    For 4 electrons, we have R (internuclear) and the hyperradius
    R_e = sqrt(r1^2 + r2^2 + r3^2 + r4^2).

    The potential terms:
    - V_ee: 6 pairs, depends only on angular coords (rho-independent for
      atom-centered origin, but rho-dependent through electron-nuclear
      distances in the molecular case)
    - V_nuc: 4x2 = 8 terms, each -Z_A/r_{iA} - Z_B/r_{iB}

    The nuclear distance for electron i:
      r_{iA} = |r_i - R_A| where R_A = (0,0,R/2)
      In hyperspherical: r_i = R_e * s_i(alpha1,alpha2,alpha3) * direction

    The ratio R/(2*R_e) = rho still appears, but now s_i depends on
    THREE hyperangles (not one). The angular problem still depends on
    (R, R_e) only through rho = R/(2*R_e), so the 2D -> 1D collapse
    DOES generalize.
    """
    return {
        'rho_collapse': True,
        'reason': (
            'All nuclear terms are -Z_X / |r_i - R_X|. In hyperspherical '
            'coords, r_i = R_e * g_i(Omega) where g_i depends only on '
            'angular coords. So r_{iA} = R_e * h_iA(Omega, R/R_e). '
            'Since R enters only as R/R_e = 2*rho, the angular eigenvalue '
            'problem depends on (R, R_e) only through rho, just as at Level 4.'
        ),
        'rho_grid_points': '~100-200 (same as Level 4)',
        'V_ee_rho_independent': True,
        'V_nuc_rho_dependent': True,
    }


# ==========================================================================
# (e) SPECTRAL BASIS: SO(12) Casimir
# ==========================================================================

def so12_casimir_spectrum(nu_max: int = 10) -> List[Dict]:
    """Compute SO(12) Casimir eigenvalues for 4-electron problem.

    mu_free(nu) = nu * (nu + 10) / 2

    The ground state for S=0 (singlet, [2,2] spatial irrep) has
    nu_min = N - lambda_1 = 4 - 2 = 2 (from Paper 16).

    So the lowest allowed free eigenvalue is:
    mu_free(2) = 2 * 12 / 2 = 12
    """
    result = []
    for nu in range(nu_max + 1):
        mu = nu * (nu + 10) / 2.0
        l_eff = nu + 4.0  # (3N-4)/2 = 4 for N=4
        result.append({
            'nu': nu,
            'mu_free': mu,
            'l_eff': l_eff,
            'degeneracy': _so_d_degeneracy(12, nu),
        })
    return result


def _so_d_degeneracy(d: int, nu: int) -> int:
    """Degeneracy of SO(d) representation with angular momentum nu.

    For SO(d), the number of linearly independent hyperspherical
    harmonics of degree nu is:

    g(d, nu) = (2*nu + d - 2) / (d - 2) * C(nu + d - 3, nu)

    Reference: Avery, "Hyperspherical Harmonics", eq. (1.3.3).
    """
    if nu == 0:
        return 1
    if d == 2:
        return 2  # SO(2) has degeneracy 2 for all nu >= 1
    # g(d, nu) = (2*nu + d - 2) * C(nu + d - 3, nu) / (d - 2)
    numerator = (2 * nu + d - 2) * comb(nu + d - 3, nu)
    assert numerator % (d - 2) == 0
    return numerator // (d - 2)


# ==========================================================================
# (f) COUPLING MATRICES
# ==========================================================================

def coupling_analysis_4e() -> Dict[str, str]:
    """Analyze V_ee and V_nuc coupling structure for 4 electrons."""
    return {
        'V_ee_pairs': 6,
        'V_ee_description': (
            'C(4,2) = 6 pairs: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4). '
            'Each pair 1/r_{ij} expands in multipoles: '
            '1/r_{ij} = (1/R_e) sum_k f_k(alpha_ij) P_k(cos theta_{ij}). '
            'In the hyperspherical tree, each pair involves different '
            'combinations of the 3 hyperangles. The multipole expansion '
            'and Gaunt integral structure still applies but involves '
            'coupling among all 4 electron angular momenta, not just 2.'
        ),
        'V_ee_gaunt': (
            'Yes, Gaunt integrals still apply for each pair. But the '
            'selection rules involve FOUR angular momenta instead of two. '
            'Each pair (i,j) couples (l_i, l_j) -> (l_i\', l_j\') while '
            'leaving the other two l-values unchanged. This creates a '
            'MUCH denser coupling structure than Level 4.'
        ),
        'V_nuc_terms': 8,
        'V_nuc_description': (
            '4 electrons x 2 nuclei = 8 nuclear attraction terms. '
            'Each -Z_X/|r_i - R_X| has the split-region Legendre form: '
            'sum_k (min(s_i,rho_X)/max(s_i,rho_X))^k P_k(cos theta_iX) / max. '
            'This is the SAME structure as Level 4 but with 4 electrons '
            'instead of 2. The nuclear coupling is diagonal in the '
            'angular momenta of ALL electrons except electron i.'
        ),
        'V_nuc_split_region': True,
    }


# ==========================================================================
# (g) FEASIBILITY ESTIMATES
# ==========================================================================

def feasibility_estimate(l_max: int, n_basis_alpha: int = 10,
                         n_rho: int = 130) -> Dict:
    """Estimate computational cost for 4-electron mol-frame solver.

    Parameters
    ----------
    l_max : int
        Maximum angular momentum per electron.
    n_basis_alpha : int
        Number of basis functions per hyperangle per channel.
    n_rho : int
        Number of rho-grid points.
    """
    # Channel count for LiH (heteronuclear, sigma only first)
    sigma = four_electron_channel_count_molecular(l_max, sigma_only=True,
                                                  homonuclear=False)
    full = four_electron_channel_count_molecular(l_max, sigma_only=False,
                                                 homonuclear=False)

    n_sigma_ch = sigma['sigma_channels']
    n_full_ch = full['full_channels']

    # For the atomic case, get L=0 coupled count
    atomic = four_electron_channel_count_atomic(l_max)
    n_coupled_L0 = atomic['coupled_L0']

    # S_4 antisymmetry reduction: roughly 1/6 for [2,2]
    n_antisym_sigma = max(1, n_sigma_ch // 6)
    n_antisym_full = max(1, n_full_ch // 6)

    # Angular matrix dimension:
    # For 4 electrons, the alpha-space is 3-dimensional (alpha1, alpha2, alpha3).
    # Spectral basis: n_basis^3 per channel (tensor product of 3 Jacobi bases).
    # Total angular matrix dim = n_channels * n_basis^3
    n_alpha_basis = n_basis_alpha ** 3  # 3 hyperangles
    angular_dim_sigma = n_antisym_sigma * n_alpha_basis
    angular_dim_full = n_antisym_full * n_alpha_basis

    # Cost per rho-point eigensolve: O(dim^3) for dense, O(dim^2) for sparse
    # (lowest eigenvalue only)
    cost_dense_sigma = angular_dim_sigma ** 3
    cost_dense_full = angular_dim_full ** 3

    # Wall time estimate (assuming ~1 GFLOP/s effective for eigensolves)
    gflops = 1e9
    time_per_rho_sigma = cost_dense_sigma / gflops  # seconds
    time_per_rho_full = cost_dense_full / gflops
    total_sigma = time_per_rho_sigma * n_rho
    total_full = time_per_rho_full * n_rho

    # Level 4 comparison
    l4_channels = level4_channel_count(l_max, sigma_only=True, homonuclear=False)
    l4_dim = l4_channels * n_basis_alpha  # 1 alpha dimension
    l4_cost = l4_dim ** 3 / gflops * n_rho

    return {
        'l_max': l_max,
        'n_basis_alpha': n_basis_alpha,
        'n_rho': n_rho,
        # Channel counts
        'sigma_channels_raw': n_sigma_ch,
        'full_channels_raw': n_full_ch,
        'coupled_L0_atomic': n_coupled_L0,
        'sigma_channels_antisym': n_antisym_sigma,
        'full_channels_antisym': n_antisym_full,
        # Angular dimensions
        'n_alpha_basis_per_channel': n_alpha_basis,
        'angular_dim_sigma': angular_dim_sigma,
        'angular_dim_full': angular_dim_full,
        # Cost
        'time_per_rho_sigma_s': time_per_rho_sigma,
        'time_per_rho_full_s': time_per_rho_full,
        'total_time_sigma_s': total_sigma,
        'total_time_full_s': total_full,
        # Comparison
        'level4_channels': l4_channels,
        'level4_dim': l4_dim,
        'level4_total_s': l4_cost,
    }


# ==========================================================================
# (h) ANTISYMMETRY ANALYSIS
# ==========================================================================

def antisymmetry_analysis() -> Dict[str, str]:
    """Analyze how 4-electron antisymmetry manifests in hyperspherical coords.

    At Level 4 (2 electrons), the gerade constraint l1+l2=even enforces
    exchange symmetry alpha -> pi/2 - alpha (swapping electrons 1 and 2).

    For 4 electrons, the permutation group is S_4 (24 elements).
    The exchange symmetry acts on the THREE hyperangles that parameterize
    the 4 radial magnitudes, plus permutations of the angular momentum labels.
    """
    return {
        'permutation_group': 'S_4 (order 24)',
        'hyperangle_action': (
            'S_4 acts on the Jacobi tree of hyperangles. For the standard '
            'tree r1,r2,r3,r4 -> (R_e, alpha1, alpha2, alpha3):\n'
            '  alpha1 = arctan(r2/r1) -- swaps electrons 1,2\n'
            '  alpha2 = arctan(r34/r12) where r12=sqrt(r1^2+r2^2), '
            'r34=sqrt(r3^2+r4^2) -- swaps pairs\n'
            '  alpha3 = arctan(r4/r3) -- swaps electrons 3,4\n'
            'Transpositions (12) and (34) act as alpha_i -> pi/2 - alpha_i '
            'on individual hyperangles (analogous to the 2e gerade constraint). '
            'The transposition (13) mixes ALL three hyperangles plus the '
            'angular momenta -- this is the structurally new element.'
        ),
        'gerade_analog': (
            'For 2 electrons: l1+l2 even (gerade) from (12) exchange.\n'
            'For 4 electrons, [2,2] symmetry of S_4 imposes:\n'
            '  - l1+l2 even under (12) exchange\n'
            '  - l3+l4 even under (34) exchange\n'
            '  - A MIXED constraint under (13), (14), (23), (24) that '
            'entangles the two pairs. This mixed constraint CANNOT be '
            'decomposed into independent constraints on each pair.\n'
            'The [2,2] projector acts on the full 4-electron channel '
            'space and eliminates roughly 5/6 of channels.'
        ),
        'implementation': (
            'The S_4 projection must be implemented as an explicit '
            'projection operator P_{[2,2]} = (dim/|G|) sum_{g in S_4} '
            'chi_{[2,2]}(g) * D(g), where D(g) is the representation '
            'matrix of g acting on the channel basis. This is a one-time '
            'setup cost (build the projector), after which the projected '
            'basis has ~1/6 the dimension.'
        ),
    }


# ==========================================================================
# (i) PAULI EXCLUSION / CORE EXCLUSION
# ==========================================================================

def pauli_exclusion_analysis() -> Dict[str, str]:
    """Analyze whether antisymmetry automatically excludes 1s^2 core.

    Key question: In a 4-electron hyperspherical calculation with full
    antisymmetry, would the 1s^2 core electrons automatically separate
    from the valence electrons?
    """
    return {
        'answer': 'PARTIAL — antisymmetry provides the CONSTRAINT but not the SEPARATION',
        'explanation': (
            'Antisymmetry (the [2,2] projector) ensures that no two electrons '
            'can occupy the same single-particle state. This means the 1s^2 '
            'pair is correctly excluded from the valence space at the level '
            'of the wavefunction.\n\n'
            'HOWEVER, the 4-electron hyperspherical solver does NOT separate '
            'core and valence electrons into independent subspaces. All 4 '
            'electrons live in the SAME 12-dimensional configuration space. '
            'The core electrons contribute to the full coupled-channel problem.\n\n'
            'This is both the ADVANTAGE and the COST of the full 4-electron '
            'approach:\n'
            '  ADVANTAGE: No PK pseudopotential needed. No Z_eff approximation. '
            'Core-valence correlation is exact. Core polarization is included.\n'
            '  COST: The 1s^2 core electrons require l_max adequate to describe '
            'BOTH core (tight, l~0) and valence (diffuse, needs higher l) '
            'simultaneously. This is the same challenge that makes Gaussian '
            'basis sets need core functions (cc-pCVXZ vs cc-pVXZ).\n\n'
            'At Level 5 (composed geometry), the PK pseudopotential APPROXIMATES '
            'the effect of antisymmetry between core and valence by enforcing '
            'orthogonality. The full 4-electron calculation replaces this '
            'approximation with exact antisymmetry, at the cost of dramatically '
            'larger angular Hilbert space.'
        ),
        'minimum_nu': (
            'From Paper 16: nu_min = N - lambda_1 = 4 - 2 = 2 for singlet [2,2]. '
            'The mu_free = 2*(2+10)/2 = 12 provides a "Pauli centrifugal barrier" '
            'that effectively keeps the 4 electrons from collapsing to the same '
            'spatial region. This is the group-theoretic encoding of the exclusion '
            'principle.'
        ),
    }


# ==========================================================================
# HYPERSPHERICAL TREE: 4-electron angular coordinate decomposition
# ==========================================================================

def hyperspherical_tree_4e() -> Dict[str, str]:
    """Describe the Jacobi tree for 4-electron hyperspherical coordinates.

    The standard Jacobi tree for 4 particles groups them hierarchically:

    Level 0: r1, r2, r3, r4 (individual electron positions)
    Level 1: pair (1,2) -> rho_12 = sqrt(r1^2 + r2^2), alpha_12 = arctan(r2/r1)
             pair (3,4) -> rho_34 = sqrt(r3^2 + r4^2), alpha_34 = arctan(r4/r3)
    Level 2: (1,2)+(3,4) -> R_e = sqrt(rho_12^2 + rho_34^2),
             alpha_1234 = arctan(rho_34/rho_12)
    """
    return {
        'tree_structure': (
            'Standard democratic tree for 4 electrons:\n'
            '  R_e = sqrt(r1^2 + r2^2 + r3^2 + r4^2)  [hyperradius]\n'
            '  alpha1 = arctan(r2/r1)                    [pair 1-2 correlation]\n'
            '  alpha2 = arctan(sqrt(r3^2+r4^2)/sqrt(r1^2+r2^2))  [inter-pair]\n'
            '  alpha3 = arctan(r4/r3)                    [pair 3-4 correlation]\n'
            'Plus 8 direction angles: theta_i, phi_i for i=1..4\n'
            'Total: 1 (R_e) + 3 (alphas) + 8 (directions) = 12 = 3*4 coordinates'
        ),
        'angular_decomposition': (
            'The 11-dimensional angular space S^11 decomposes as:\n'
            '  3 hyperangles alpha_1, alpha_2, alpha_3 in [0, pi/2]\n'
            '  8 direction angles (theta_i, phi_i) for 4 electrons\n'
            'The channel expansion integrates out the 8 direction angles '
            'into angular momentum labels (l_i, m_i), leaving the 3 '
            'hyperangles as the continuous variables.'
        ),
        'alpha_range': '[0, pi/2] for each of the 3 hyperangles',
        'volume_element': (
            'dV = R_e^11 * sin^2(alpha1)*cos^2(alpha1) * '
            'sin^2(alpha2)*cos^2(alpha2) * sin^2(alpha3)*cos^2(alpha3) * '
            'prod_i sin(theta_i) * dR_e * d(alphas) * d(angles)'
        ),
        'cusp_locations': (
            'e-e coalescence r_{ij} = 0:\n'
            '  Pair (1,2): alpha1 = pi/4, theta_{12} = 0\n'
            '  Pair (3,4): alpha3 = pi/4, theta_{34} = 0\n'
            '  Cross pairs (1,3), (1,4), (2,3), (2,4): involve ALL 3 '
            'hyperangles simultaneously plus inter-electron angles. '
            'These cross-pair cusps are the structurally new challenge: '
            'they live in a 2D submanifold of the 11D angular space '
            'that cuts across the Jacobi tree.'
        ),
    }
