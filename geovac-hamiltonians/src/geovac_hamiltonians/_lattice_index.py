"""
LatticeIndex — Pre-calculated topology for N-electron graph Hamiltonians
=========================================================================

Architecture: the many-body state space is treated as a relational database.

  - sp_states: single-particle states (n, l, m, sigma) — primary keys
  - conformal_weights[i]: Omega_i = 2*p0/(p_i^2 + p0^2), pre-computed once
  - h1_diag[sp], h1_offdiag[sp]: one-electron integrals indexed by spatial state
  - vee_matrix[sp_i, sp_j]: two-electron repulsion (Mulliken approximation)
  - sd_basis: enumerated Slater determinants (sorted tuples of occupied sp indices)

Hamiltonian assembly is a sparse JOIN query:
  For each occupied orbital p in each SD I:
    For each graph-adjacent r of spatial(p):
      If r unoccupied: H[I, J] += phase * h1(spatial_p, spatial_r)
  plus diagonal V_ee correction.

Elements below sparsity_threshold (default 1e-8) are dropped before the
matrix constructor is called — this is the topological sparsity mask.

KNOWN SYSTEMATIC ERRORS (logged at construction, not papered over):
  1. V_ee graph approximation (default: 'chordal'):
     Uses S³ chordal distance V_ee = κ_ee / d²_chord, where d_chord is the
     Euclidean distance in R⁴ between Fock-projected S³ coordinates. This is
     architecturally consistent with the graph Laplacian kinetic term (Paper 7,
     Eq. 14) but does not reproduce exact Hartree two-electron integrals.
     Fallback 'mulliken': uses mean orbital radius r_eff=n²/Z with ~12%
     systematic overestimation; factor ~10x error for same-orbital pairs.
  2. Topological bond contraction: the graph Laplacian uses uniform edge
     weights (W=1), not the exact radial matrix elements. This under-represents
     transitions between shells of very different n.

Reference implementation: MoleculeHamiltonian._solve_full_ci() for 2-electron
systems (geovac/hamiltonian.py). This module extends that architecture to
arbitrary N.

Author: GeoVac Development Team
Date: February 2026
"""

import os
import warnings
import time
from itertools import combinations
from math import factorial
from typing import Dict, List, Optional, Tuple

import numpy as np
import scipy.sparse as sp
from scipy.sparse import csr_matrix, diags, identity
from scipy.sparse.linalg import eigsh
from scipy.integrate import quad

from geovac_hamiltonians._lattice import GeometricLattice

KINETIC_SCALE: float = -1.0 / 16.0   # Universal topological constant


# ---------------------------------------------------------------------------
# S³ density-overlap V_ee (Paper 7, Section VI)
# ---------------------------------------------------------------------------

def _phi_s_orbital(n: int, t: float) -> float:
    """
    Fock-projected density for s-orbital (l=0) with principal quantum number n.

    Phi_n(t) = rho_tilde_{ns,ns}(2Z*t) evaluated with Z=1 (Z scaling handled
    separately by the master formula prefactor 4Z/pi).

    Parameters
    ----------
    n : int
        Principal quantum number (n >= 1)
    t : float
        Dimensionless momentum transfer t = q/(2Z)

    Returns
    -------
    float
        Projected density value at t
    """
    if n == 1:
        # Phi_1s(t) = 1/(1+t^2)^2
        return 1.0 / (1.0 + t * t) ** 2
    elif n == 2:
        # Phi_2s(t) = (1-12t^2+32t^4)/(4t^2+1)^4
        t2 = t * t
        return (1.0 - 12.0 * t2 + 32.0 * t2 * t2) / (4.0 * t2 + 1.0) ** 4
    else:
        # General n: compute from hydrogen density FT with p0=2 (Z=1 units)
        # rho_ns(q) for hydrogen, then substitute q=2t
        # Use numerical evaluation via position-space radial wavefunctions
        return _phi_s_orbital_general(n, t)


def _phi_s_orbital_general(n: int, t: float) -> float:
    """
    General s-orbital projected density via numerical FT of |R_{n0}(r)|^2.

    Phi_n(t) = rho_tilde(q=2t) with Z=1, where
    rho_tilde(q) = integral_0^inf |R_{n0}(r)|^2 * sin(qr)/(qr) * r^2 dr

    Delegates to _form_factor_nl which has the correct hydrogenic
    normalization for all (n, l).

    For n=1,2 the closed-form expressions in _phi_s_orbital are preferred.
    """
    # q = 2*t in Z=1 units
    return _form_factor_nl(n, 0, 1.0, 2.0 * t)


def _form_factor_nl(n: int, l: int, Z: float, q: float) -> float:
    """
    Spherically averaged charge form factor for hydrogenic (n,l) orbital.

    rho_tilde(q) = integral_0^inf |R_{nl,Z}(r)|^2 * sin(qr)/(qr) * r^2 dr

    For s-orbitals (l=0) with n<=2, uses closed-form expressions.
    Otherwise, uses numerical integration of hydrogenic radial wavefunctions.

    The form factor satisfies rho_tilde(0) = 1 (normalization) and
    rho_tilde(q) -> 0 as q -> inf.

    Parameters
    ----------
    n : int
        Principal quantum number
    l : int
        Angular momentum quantum number
    Z : float
        Nuclear charge
    q : float
        Momentum transfer in atomic units (1/bohr)

    Returns
    -------
    float
        Form factor value (dimensionless)
    """
    # Closed-form for n=1,2 s-orbitals (verified against known Slater F0)
    if l == 0 and n <= 2:
        return _phi_s_orbital(n, q / (2.0 * Z))

    # General (n, l, Z) via numerical integration of |R_{nl,Z}(r)|²
    from scipy.special import assoc_laguerre
    from math import factorial

    # Hydrogenic radial wavefunction R_{nl,Z}(r):
    # R_{nl}(r) = N * (2Zr/n)^l * L_{n-l-1}^{2l+1}(2Zr/n) * exp(-Zr/n)
    # N = (2Z/n)^{3/2} * sqrt((n-l-1)! / (2n * ((n+l)!)^3)) * (n+l)!
    nr = n - l - 1  # radial quantum number
    N_sq = (2.0 * Z / n) ** 3 * factorial(nr) / (2.0 * n * (factorial(n + l)) ** 3)
    N_sq *= factorial(n + l) ** 2  # from the associated Laguerre normalization
    # Full: N² = (2Z/n)³ × (n-l-1)! / (2n × (n+l)!) — standard Griffiths convention

    def integrand(r: float) -> float:
        if r < 1e-30:
            return 0.0
        x = 2.0 * Z * r / n
        lag = assoc_laguerre(x, nr, 2 * l + 1)
        # |R_{nl}|² = N² × x^{2l} × L² × exp(-x)
        R_sq = N_sq * x ** (2 * l) * lag * lag * np.exp(-x)
        rho_r = R_sq * r * r  # |R_{nl}|² r²
        if q < 1e-14:
            return rho_r  # sin(qr)/(qr) -> 1
        qr = q * r
        return rho_r * np.sin(qr) / qr

    val, _ = quad(integrand, 0, np.inf, limit=300)
    return val


def _form_factor_nl_grid(
    n: int, l: int, Z: float, q_grid: np.ndarray,
    n_r: int = 500,
) -> np.ndarray:
    """
    Vectorized form factor evaluation on a q-grid via direct radial integration.

    Computes rho_tilde(q) = integral_0^inf |R_{nl,Z}(r)|^2 * sin(qr)/(qr) * r^2 dr
    for all q values simultaneously, using trapezoidal integration on a dense
    exponentially-scaled r-grid.

    For s-orbitals (l=0) with n<=2, uses closed-form expressions (vectorized).

    Parameters
    ----------
    n : int
        Principal quantum number
    l : int
        Angular momentum quantum number
    Z : float
        Nuclear charge
    q_grid : np.ndarray
        1D array of momentum transfer values (1/bohr)
    n_r : int
        Number of radial grid points (default 500)

    Returns
    -------
    np.ndarray
        Form factor values at each q point (same shape as q_grid)
    """
    # Closed-form for n=1,2 s-orbitals (vectorized)
    if l == 0 and n <= 2:
        t = q_grid / (2.0 * Z)
        if n == 1:
            return 1.0 / (1.0 + t * t) ** 2
        else:  # n == 2
            t2 = t * t
            return (1.0 - 12.0 * t2 + 32.0 * t2 * t2) / (4.0 * t2 + 1.0) ** 4

    # General (n, l, Z): dense r-grid with exponential scaling
    from scipy.special import assoc_laguerre

    nr = n - l - 1  # radial quantum number
    N_sq = (2.0 * Z / n) ** 3 * factorial(nr) / (2.0 * n * (factorial(n + l)) ** 3)
    N_sq *= factorial(n + l) ** 2

    # Radial grid: r from 0 to r_max where |R_nl|^2 is negligible
    # The wavefunction decays as exp(-Z*r/n), so exp(-2*Z*r_max/n) ~ 1e-20
    r_max = n * 50.0 / Z  # ~23 decay lengths
    r_grid = np.linspace(1e-12, r_max, n_r)

    # Evaluate |R_{nl}(r)|^2 * r^2 on the grid
    rho_r = 2.0 * Z * r_grid / n  # dimensionless variable
    lag_vals = np.array([assoc_laguerre(xi, nr, 2 * l + 1) for xi in rho_r])
    R_sq = N_sq * rho_r ** (2 * l) * lag_vals ** 2 * np.exp(-rho_r)
    radial_density = R_sq * r_grid ** 2  # |R_nl|^2 * r^2, shape (n_r,)

    # sinc(qr) for all (q, r) pairs: shape (n_q, n_r)
    qr = np.outer(q_grid, r_grid)
    sinc_qr = np.ones_like(qr)
    mask = qr > 1e-14
    sinc_qr[mask] = np.sin(qr[mask]) / qr[mask]

    # Trapezoidal integration: rho(q) = integral |R_nl|^2 * r^2 * sinc(qr) dr
    dr = r_grid[1] - r_grid[0]
    result = (sinc_qr * radial_density[np.newaxis, :]).sum(axis=1) * dr

    return result


def compute_cross_atom_J_batch(
    pairs: List[Tuple[int, int, int, int]],
    R: float,
    ZA: float,
    ZB: float,
    n_q: int = 2000,
    q_max_factor: float = 8.0,
) -> Dict[Tuple[int, int, int, int], float]:
    """
    Batch-compute cross-atom direct Coulomb integrals for multiple orbital pairs.

    Uses a shared q-grid with tabulated form factors to avoid nested numerical
    integration. ~100x faster than individual compute_cross_atom_J calls for
    l>0 orbitals.

    J_AB(R) = (2/pi) integral_0^inf rho_A(q) rho_B(q) sin(qR)/(qR) dq

    Parameters
    ----------
    pairs : list of (na, la, nb, lb) tuples
        Orbital quantum number pairs to compute
    R : float
        Internuclear distance in bohr
    ZA, ZB : float
        Nuclear charges of centers A and B
    n_q : int
        Number of q-grid points for trapezoidal integration
    q_max_factor : float
        q_max = q_max_factor * max(ZA, ZB) — controls grid extent

    Returns
    -------
    dict
        Maps (na, la, nb, lb) -> J value in Hartree
    """
    if not pairs:
        return {}

    # Check disk cache first — return cached values, compute only missing
    cache_dir = os.path.join(os.path.dirname(__file__), 'cache')
    os.makedirs(cache_dir, exist_ok=True)

    results: Dict[Tuple[int, int, int, int], float] = {}
    uncached_pairs: List[Tuple[int, int, int, int]] = []

    for key in pairs:
        na, la, nb, lb = key
        cache_file = os.path.join(
            cache_dir,
            f"cross_atom_J_{na}{la}_{nb}{lb}_R{R:.4f}_Z{ZA:.0f}_{ZB:.0f}.npy",
        )
        if os.path.exists(cache_file):
            results[key] = float(np.load(cache_file))
        else:
            uncached_pairs.append(key)

    if not uncached_pairs:
        return results

    # Build q-grid: denser near origin where form factors peak
    Z_max = max(ZA, ZB)
    q_max = q_max_factor * Z_max
    # Use sqrt-spaced grid for better resolution near q=0
    q_lin = np.linspace(0, np.sqrt(q_max), n_q) ** 2

    # Collect unique (n, l, Z) combinations needed
    ff_keys_A: set = set()
    ff_keys_B: set = set()
    for na, la, nb, lb in uncached_pairs:
        ff_keys_A.add((na, la))
        ff_keys_B.add((nb, lb))

    # Tabulate form factors on q-grid (one vectorized call per unique orbital)
    ff_table_A: Dict[Tuple[int, int], np.ndarray] = {}
    for n, l in ff_keys_A:
        ff_table_A[(n, l)] = _form_factor_nl_grid(n, l, ZA, q_lin)

    ff_table_B: Dict[Tuple[int, int], np.ndarray] = {}
    for n, l in ff_keys_B:
        ff_table_B[(n, l)] = _form_factor_nl_grid(n, l, ZB, q_lin)

    # Precompute sin(qR)/(qR) on grid
    qR = q_lin * R
    sinc_qR = np.ones_like(qR)
    mask = qR > 1e-14
    sinc_qR[mask] = np.sin(qR[mask]) / qR[mask]

    # Compute J for each pair via trapezoidal integration on the q-grid
    for key in uncached_pairs:
        na, la, nb, lb = key
        integrand = ff_table_A[(na, la)] * ff_table_B[(nb, lb)] * sinc_qR
        j_ab = (2.0 / np.pi) * np.trapezoid(integrand, q_lin)
        results[key] = j_ab

        # Save to disk cache
        cache_file = os.path.join(
            cache_dir,
            f"cross_atom_J_{na}{la}_{nb}{lb}_R{R:.4f}_Z{ZA:.0f}_{ZB:.0f}.npy",
        )
        np.save(cache_file, j_ab)

    return results


def compute_vee_s3_overlap(
    n_a: int,
    l_a: int,
    n_b: int,
    l_b: int,
    Z: float,
    p0: float | None = None,
) -> float:
    """
    Compute the direct Coulomb (Slater F0) integral via S3 density-overlap.

    V_ee is a NODE property on S3: each orbital contributes a density profile
    Phi(t) and F0 = (4Z/pi) * integral_0^inf Phi_a(t) * Phi_b(t) dt.

    This is NOT an edge property -- do not use chordal distance d2_chord for
    l>0 orbitals (overestimates 1s-2s by ~30x).

    Currently implemented for s-orbital pairs (l_a = l_b = 0).
    For l>0, the full 3D angular convolution with spherical harmonic multipoles
    is required (deferred to future work -- see Paper 7 Section VI.E).

    Parameters
    ----------
    n_a, l_a : int
        Quantum numbers for orbital a
    n_b, l_b : int
        Quantum numbers for orbital b
    Z : float
        Nuclear charge
    p0 : float or None
        Fock projection scale (default: 2Z for s-orbitals)

    Returns
    -------
    float
        Slater F0 direct Coulomb integral in Hartree atomic units

    Raises
    ------
    NotImplementedError
        If l_a > 0 or l_b > 0 (angular extension pending, see Paper 7 Sec VI.E)
    """
    if l_a > 0 or l_b > 0:
        raise NotImplementedError(
            f"S3 density-overlap V_ee not implemented for l>0 orbitals "
            f"(got l_a={l_a}, l_b={l_b}). The full 3D angular convolution "
            f"with spherical harmonic multipoles is required. "
            f"See Paper 7 Section VI.E. Do NOT fall back to chordal d2 -- "
            f"it overestimates inter-shell integrals by ~30x."
        )

    def integrand(t: float) -> float:
        return _phi_s_orbital(n_a, t) * _phi_s_orbital(n_b, t)

    integral, _ = quad(integrand, 0, np.inf, limit=200)
    return (4.0 * Z / np.pi) * integral


def _wavefunction_form_factor(n: int, l: int, Z: float, q: float) -> float:
    """
    Spherically averaged wavefunction form factor for hydrogenic (n,l) orbital.

    g(q) = integral_0^inf R_{nl,Z}(r) * sin(qr)/(qr) * r^2 dr

    This is the Fourier-Bessel transform of the RADIAL WAVEFUNCTION R_{nl}(r),
    NOT the density |R_{nl}|^2. Used for computing overlap integrals between
    orbitals on different centers.

    Unlike the density form factor (rho_tilde(0) = 1 always), g(0) = integral R r^2 dr
    and is NOT normalized to 1 (e.g., g_{1s}(0) = 4Z^{5/2}/(Z^2)^2 = 4/Z^{3/2}).

    Parameters
    ----------
    n : int
        Principal quantum number
    l : int
        Angular momentum quantum number (must be 0 for now)
    Z : float
        Nuclear charge
    q : float
        Momentum transfer in atomic units (1/bohr)

    Returns
    -------
    float
        Wavefunction form factor value
    """
    # Closed-form expressions using Laplace transforms of r^n e^{-ar} sin(qr)
    # General formula: g(q) = int_0^inf R_{nl}(r) j_0(qr) r^2 dr
    #   = (1/q) int_0^inf R_{nl}(r) r sin(qr) dr  for l=0

    if l == 0 and n == 1:
        # R_{10}(r) = 2Z^{3/2} e^{-Zr}
        # g(q) = 2Z^{3/2} * 2Z/(Z^2+q^2)^2 = 4Z^{5/2}/(Z^2+q^2)^2
        # using int_0^inf r e^{-ar} sin(qr) dr = 2aq/(a^2+q^2)^2
        # and g = (2Z^{3/2}/q) * 2Zq/(Z^2+q^2)^2
        return 4.0 * Z ** 2.5 / (Z * Z + q * q) ** 2

    if l == 0 and n == 2:
        # R_{20}(r) = Z^{3/2}/(2*sqrt(2)) * (2-Zr) * e^{-Zr/2}
        # g(q) = Z^{3/2}/(2*sqrt(2)) * [2*I1 - Z*I2]
        # where a = Z/2,
        #   I1 = int r^2 e^{-ar} sin(qr)/(qr) dr = 2a/(a^2+q^2)^2
        #   I2 = int r^3 e^{-ar} sin(qr)/(qr) dr = 2(3a^2-q^2)/(a^2+q^2)^3
        from math import sqrt
        a = Z / 2.0
        a2q2 = a * a + q * q
        I1 = 2.0 * a / a2q2 ** 2
        I2 = 2.0 * (3.0 * a * a - q * q) / a2q2 ** 3
        return Z ** 1.5 / (2.0 * sqrt(2.0)) * (2.0 * I1 - Z * I2)

    if l == 0 and n == 3:
        # g(q) = 4*a^{5/2} * (3a^4 - 10a^2*q^2 + 3q^4) / (a^2+q^2)^4
        # where a = Z/3. Derived from Laplace transforms:
        # R_{30}(r) = N*(L_2^1(2Zr/3))*e^{-Zr/3}, expanded into r^k e^{-ar} terms,
        # then int r^k e^{-ar} sin(qr)/(qr) dr evaluated analytically.
        # Verified against numerical integration to machine precision.
        a = Z / 3.0
        a2 = a * a
        q2 = q * q
        D = a2 + q2
        num = 3.0 * a2 * a2 - 10.0 * a2 * q2 + 3.0 * q2 * q2
        return 4.0 * a ** 2.5 * num / D ** 4

    # General case: numerical integration
    from scipy.special import assoc_laguerre
    from math import factorial, sqrt

    nr = n - l - 1
    N_coeff = (2.0 * Z / n) ** 1.5 * sqrt(factorial(nr) / (2.0 * n * factorial(n + l)))

    def integrand(r: float) -> float:
        if r < 1e-30:
            return 0.0
        x = 2.0 * Z * r / n
        lag = assoc_laguerre(x, nr, 2 * l + 1)
        R_nl = N_coeff * x ** l * lag * np.exp(-x / 2.0)
        if q < 1e-14:
            return R_nl * r * r
        qr = q * r
        return R_nl * np.sin(qr) / qr * r * r

    val, _ = quad(integrand, 0, np.inf, limit=300)
    return val


def compute_overlap_element(
    na: int, la: int, nb: int, lb: int,
    ZA: float, ZB: float, R: float,
) -> float:
    """
    Overlap integral between hydrogenic orbitals on different centers.

    S_AB = (2/pi) integral_0^inf g_A(q) g_B(q) sin(qR)/(qR) q^2 dq

    where g(q) = integral_0^inf R_{nl,Z}(r) sin(qr)/(qr) r^2 dr
    is the wavefunction form factor (NOT the density form factor).

    For R=0, same-center same-orbital: S = 1 (orthonormal atomic basis).
    For R>0, cross-center: S measures basis overlap between atom-centered
    hydrogenic orbitals, driving BSSE.

    Currently restricted to l=0 (s-orbital) pairs. The l>0 generalization
    requires angular-dependent form factors.

    Parameters
    ----------
    na, la : int
        Quantum numbers for orbital on center A
    nb, lb : int
        Quantum numbers for orbital on center B
    ZA, ZB : float
        Nuclear charges
    R : float
        Internuclear distance in bohr

    Returns
    -------
    float
        Overlap integral (dimensionless)

    Raises
    ------
    NotImplementedError
        If la > 0 or lb > 0
    """
    if la > 0 or lb > 0:
        raise NotImplementedError(
            f"Overlap form factor not implemented for l>0 orbitals "
            f"(got la={la}, lb={lb}). Requires angular-dependent form factors."
        )

    def integrand(q: float) -> float:
        gA = _wavefunction_form_factor(na, la, ZA, q)
        gB = _wavefunction_form_factor(nb, lb, ZB, q)
        qR = q * R
        if qR < 1e-10:
            sinc = 1.0
        else:
            sinc = np.sin(qR) / qR
        return gA * gB * q * q * sinc

    result, _ = quad(integrand, 0, np.inf, limit=300)
    return (2.0 / np.pi) * result


def compute_exact_cross_nuclear(
    n_B: int,
    l_B: int,
    m_B: int,
    Z_B: float,
    Z_A: float,
    R: float,
    n_radial: int = 100,
    n_angular: int = 50,
) -> float:
    """
    Exact two-center nuclear attraction integral via 2D quadrature.

    Computes V = -Z_A * integral |chi^B_{nlm}(r)|^2 / |r - R*zhat| d^3r
    where chi^B is a hydrogen-like orbital centered at B with nuclear charge Z_B,
    and nucleus A with charge Z_A sits at distance R along the z-axis.

    Uses Gauss-Laguerre (radial) x Gauss-Legendre (angular) quadrature.
    For l=0 (s-orbitals), this is mathematically equivalent to the shell-theorem
    formula in _fourier_cross_attraction.  For l>0, the angular dependence of
    |Y_l^m|^2 produces corrections beyond the spherical average.

    Parameters
    ----------
    n_B, l_B, m_B : int
        Quantum numbers of the orbital on atom B.
    Z_B : float
        Nuclear charge of atom B (determines orbital shape).
    Z_A : float
        Nuclear charge of atom A (the attracting nucleus).
    R : float
        Internuclear distance in Bohr.
    n_radial : int
        Number of Gauss-Laguerre radial quadrature points (default 100).
    n_angular : int
        Number of Gauss-Legendre angular quadrature points (default 50).

    Returns
    -------
    float
        Cross-nuclear attraction energy in Hartree (negative).
    """
    if R < 1e-10 or Z_A == 0:
        return 0.0

    from scipy.special import lpmv

    # --- Radial quadrature: Gauss-Laguerre with hydrogenic scaling ---
    # The dominant decay of |R_nl(r)|^2 ~ exp(-2*Z_B*r/n_B).
    # Substitution x = 2*Z_B*r/n_B maps this to exp(-x).
    x_nodes, x_weights = np.polynomial.laguerre.laggauss(n_radial)
    scale = float(n_B) / (2.0 * float(Z_B))  # r = scale * x
    r_nodes = scale * x_nodes  # physical radial coordinates

    # --- Angular quadrature: Gauss-Legendre over cos(theta) ---
    cos_nodes, cos_weights = np.polynomial.legendre.leggauss(n_angular)

    # --- Radial wavefunction |R_{nl}(r)|^2 at quadrature nodes ---
    from scipy.special import assoc_laguerre

    Z = float(Z_B)
    n = n_B
    l = l_B
    m_abs = abs(m_B)
    nr = n - l - 1  # radial quantum number

    # Normalization: N^2 = (2Z/n)^3 * (n-l-1)! / (2n * (n+l)!)
    N_sq = (2.0 * Z / n) ** 3 * factorial(nr) / (2.0 * n * factorial(n + l))

    # Compute |R_nl(r)|^2 at each radial node
    rho_r = 2.0 * Z * r_nodes / n  # dimensionless radial variable
    lag_vals = np.array([assoc_laguerre(xi, nr, 2 * l + 1) for xi in rho_r])
    R_sq = N_sq * rho_r ** (2 * l) * lag_vals ** 2 * np.exp(-rho_r)
    # Note: exp(-rho) = exp(-x) is absorbed by Laguerre weights,
    # so the effective integrand multiplies by exp(+x).
    R_sq_eff = N_sq * rho_r ** (2 * l) * lag_vals ** 2  # without exp(-x)

    # --- Angular density factor ---
    # |Y_l^m(theta,phi)|^2 integrated over phi gives:
    #   angular_density(cos_theta) = (2l+1)/2 * (l-|m|)!/(l+|m|)! * [P_l^|m|(cos_theta)]^2
    # Normalized so that integral_{-1}^{1} angular_density(u) du = 1.
    if l == 0:
        ang_factor = 0.5 * np.ones(n_angular)
    else:
        c_lm = (2 * l + 1) / 2.0 * factorial(l - m_abs) / factorial(l + m_abs)
        plm_vals = np.array([lpmv(m_abs, l, u) for u in cos_nodes])
        ang_factor = c_lm * plm_vals ** 2

    # --- 2D quadrature: sum over (r, cos_theta) ---
    # V = -Z_A * integral |psi|^2 / |r - R*zhat| d^3r
    #   = -Z_A * integral_0^inf integral_{-1}^{1}
    #       |R_nl(r)|^2 * angular_density(u) / sqrt(r^2 + R^2 - 2*r*R*u) * r^2 dr du
    #
    # With substitution r = scale * x, dr = scale * dx, and Laguerre weight e^{-x}:
    #   integral ... dr = scale * sum_i w_i * f(scale*x_i) * e^{x_i} * e^{-x_i}
    #   The e^{-x} in |R_nl|^2 is cancelled by the e^{+x} from the weight change.

    total = 0.0
    R2 = R * R
    for i in range(n_radial):
        r_i = r_nodes[i]
        r2_i = r_i * r_i
        R_sq_i = R_sq_eff[i]  # |R_nl|^2 without exp(-x), Laguerre weight handles it
        r2_factor = r2_i * R_sq_i * scale  # r^2 * |R_nl|^2_eff * dr/dx

        # Angular sum
        ang_sum = 0.0
        for j in range(n_angular):
            u = cos_nodes[j]
            dist2 = r2_i + R2 - 2.0 * r_i * R * u
            if dist2 < 1e-30:
                dist2 = 1e-30
            dist = np.sqrt(dist2)
            ang_sum += cos_weights[j] * ang_factor[j] / dist

        total += x_weights[i] * r2_factor * ang_sum

    return -float(Z_A) * total


_f0_memo: Dict[Tuple[int, int, int, int, float], float] = {}


def _compute_single_f0(
    n1: int, l1: int, n2: int, l2: int, Z: float
) -> float:
    """
    Compute a single Slater F0 integral for hydrogenic orbitals at charge Z.

    F0(n1 l1, n2 l2; Z) = integral |R_{n1l1,Z}(r1)|^2 |R_{n2l2,Z}(r2)|^2
                            * (1/r_>) * r1^2 r2^2 dr1 dr2

    Used for recomputing same-atom V_ee with scaled orbital exponents.
    Results are memoized by (n1, l1, n2, l2, round(Z, 8)).
    """
    key = (n1, l1, n2, l2, round(Z, 8))
    if key in _f0_memo:
        return _f0_memo[key]
    val = _compute_single_f0_impl(n1, l1, n2, l2, Z)
    _f0_memo[key] = val
    # Also cache reverse
    _f0_memo[(n2, l2, n1, l1, round(Z, 8))] = val
    return val


def _compute_single_f0_impl(
    n1: int, l1: int, n2: int, l2: int, Z: float
) -> float:
    """Unmemorized implementation of _compute_single_f0."""
    from scipy.special import genlaguerre

    r_max = 80.0 / Z

    def _radial_wf_unnorm(r: np.ndarray, n: int, l: int) -> np.ndarray:
        rho = 2.0 * Z * r / n
        L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
        return rho ** l * np.exp(-rho / 2.0) * L_poly

    # Normalization constants
    norm1_sq, _ = quad(
        lambda r: _radial_wf_unnorm(np.array(r), n1, l1) ** 2 * r ** 2,
        0, r_max, limit=200
    )
    norm1 = 1.0 / np.sqrt(norm1_sq)
    norm2_sq, _ = quad(
        lambda r: _radial_wf_unnorm(np.array(r), n2, l2) ** 2 * r ** 2,
        0, r_max, limit=200
    )
    norm2 = 1.0 / np.sqrt(norm2_sq)

    def R_wf(r: float, n: int, l: int, norm: float) -> float:
        return float(norm * _radial_wf_unnorm(np.array(r), n, l))

    def _inner(r1: float) -> float:
        if r1 < 1e-30:
            return 0.0
        p1, _ = quad(lambda r2: R_wf(r2, n2, l2, norm2) ** 2 * r2 ** 2,
                     0, r1, limit=100)
        p2, _ = quad(lambda r2: R_wf(r2, n2, l2, norm2) ** 2 * r2,
                     r1, r_max, limit=100)
        return p1 / r1 + p2

    val, _ = quad(
        lambda r1: R_wf(r1, n1, l1, norm1) ** 2 * _inner(r1) * r1 ** 2,
        0, r_max, limit=200
    )
    return val


def compute_cross_atom_J(
    na: int,
    la: int,
    nb: int,
    lb: int,
    R: float,
    ZA: float,
    ZB: float,
) -> float:
    """
    Cross-atom direct Coulomb integral via Fourier convolution.

    J_AB(R) = (2/pi) integral_0^inf rho_A(q) rho_B(q) sin(qR)/(qR) dq

    where rho_A(q) = _form_factor_nl(na, la, ZA, q) is the spherically
    averaged charge form factor of a hydrogenic orbital on center A.

    This gives the classical Coulomb repulsion between two charge
    distributions centered on nuclei separated by distance R.
    For R -> inf, J -> 1/R (point charge limit).

    For l>0 orbitals, the density |Y_l^m|^2 is anisotropic, but the
    form factor uses the spherical average (monopole L=0 term of the
    multipole expansion). This captures 80-90% of the integral at
    typical bond lengths (R ~ 3 bohr) and is exact for s-orbitals.

    Parameters
    ----------
    na, la : int
        Quantum numbers for orbital on center A
    nb, lb : int
        Quantum numbers for orbital on center B
    R : float
        Internuclear distance in bohr
    ZA, ZB : float
        Nuclear charges of centers A and B

    Returns
    -------
    float
        Direct Coulomb integral J_AB in Hartree
    """

    # Check disk cache
    cache_dir = os.path.join(os.path.dirname(__file__), 'cache')
    cache_file = os.path.join(
        cache_dir,
        f"cross_atom_J_{na}{la}_{nb}{lb}_R{R:.4f}_Z{ZA:.0f}_{ZB:.0f}.npy",
    )
    if os.path.exists(cache_file):
        return float(np.load(cache_file))

    def integrand(q: float) -> float:
        rho_a = _form_factor_nl(na, la, ZA, q)
        rho_b = _form_factor_nl(nb, lb, ZB, q)
        qR = q * R
        if qR < 1e-10:
            sinc = 1.0
        else:
            sinc = np.sin(qR) / qR
        return rho_a * rho_b * sinc

    result, _ = quad(integrand, 0, np.inf, limit=300)
    j_ab = (2.0 / np.pi) * result

    # Save to disk cache
    os.makedirs(cache_dir, exist_ok=True)
    np.save(cache_file, j_ab)

    return j_ab


def compute_cross_atom_K(
    na: int,
    la: int,
    nb: int,
    lb: int,
    R: float,
    ZA: float,
    ZB: float,
    eri_A: Dict[Tuple[int, int, int, int], float],
    eri_B: Dict[Tuple[int, int, int, int], float],
    states_A: List[Tuple[int, int, int]],
    states_B: List[Tuple[int, int, int]],
) -> float:
    """
    Cross-atom exchange integral via Mulliken approximation.

    K(aA, bB; R) = S(aA, bB; R)^2 * [F0(a,a; ZA) + F0(b,b; ZB)] / 2

    where S is the overlap integral between orbitals on different centers,
    and F0(a,a) is the same-atom self-Coulomb integral (Slater F^0).

    The Mulliken approximation is used because the exact cross-atom exchange
    integral is a two-center integral that does not factorize into
    single-center form factors (unlike the direct Coulomb J_AB).

    Properties:
    - K -> 0 as R -> inf (S -> 0, correct dissociation)
    - K > 0 (exchange integrals are positive)
    - K < J for the same orbital pair (exchange < Coulomb)
    - Standard in quantum chemistry (Mulliken, 1949)

    Restricted to l=0 (s-orbital) pairs, consistent with the cross-atom
    J restriction.

    Parameters
    ----------
    na, la : int
        Quantum numbers for orbital on center A
    nb, lb : int
        Quantum numbers for orbital on center B
    R : float
        Internuclear distance in bohr
    ZA, ZB : float
        Nuclear charges of centers A and B
    eri_A, eri_B : dict
        Same-atom ERI tables from LatticeIndex._eri for atoms A and B
    states_A, states_B : list of (n, l, m) tuples
        State lists for atoms A and B

    Returns
    -------
    float
        Exchange integral K_AB in Hartree

    Raises
    ------
    NotImplementedError
        If la > 0 or lb > 0
    """
    if la > 0 or lb > 0:
        raise NotImplementedError(
            f"Cross-atom K restricted to s-orbitals (got la={la}, lb={lb})."
        )

    # Overlap integral S(aA, bB)
    S_ab = compute_overlap_element(na, la, nb, lb, ZA, ZB, R)

    # Find spatial indices for the (n, l=0, m=0) states
    idx_a = None
    for i, (n, l, m) in enumerate(states_A):
        if n == na and l == la and m == 0:
            idx_a = i
            break
    idx_b = None
    for i, (n, l, m) in enumerate(states_B):
        if n == nb and l == lb and m == 0:
            idx_b = i
            break

    if idx_a is None or idx_b is None:
        return 0.0

    # Same-atom self-Coulomb F0(a,a) and F0(b,b)
    f0_aa = eri_A.get((idx_a, idx_a, idx_a, idx_a), 0.0)
    f0_bb = eri_B.get((idx_b, idx_b, idx_b, idx_b), 0.0)

    return S_ab * S_ab * (f0_aa + f0_bb) / 2.0


class LatticeIndex:
    """
    Pre-calculated topology index for N-electron graph Hamiltonians.

    For Lithium (n_electrons=3, nuclear_charge=3), the state space is the
    antisymmetrized tensor product of single-particle lattices truncated at
    max_n. Nodes are Slater determinant configurations. Edges are non-zero
    off-diagonal Hamiltonian elements derived from the graph Laplacian.

    Parameters
    ----------
    n_electrons : int
        Number of electrons (3 for Li, 4 for LiH, etc.)
    max_n : int
        Maximum principal quantum number (truncation)
    nuclear_charge : int
        Nuclear charge Z (default 3 for Li)
    sparsity_threshold : float
        Matrix elements below this are dropped during assembly (default 1e-8)

    Attributes
    ----------
    n_sp : int
        Number of spin-orbitals (2 * number of spatial states)
    n_sd : int
        Number of Slater determinants = C(n_sp, n_electrons)
    sp_states : list of (n, l, m, sigma)
        Single-particle basis indexed by integer sp_idx
    conformal_weights : np.ndarray, shape (n_spatial,)
        Conformal factor Omega_i = 2*p0 / (p_i^2 + p0^2) per spatial state
    """

    def __init__(
        self,
        n_electrons: int,
        max_n: int,
        nuclear_charge: int = 3,
        sparsity_threshold: float = 1e-8,
        vee_method: str = 'chordal',
        h1_method: Optional[str] = None,
        fci_method: str = 'auto',
        kinetic_scale: float = KINETIC_SCALE,
        adjacency: Optional[csr_matrix] = None,
        node_weights: Optional[np.ndarray] = None,
        enumerate_sds: bool = True,
    ) -> None:
        if vee_method == 'chordal':
            warnings.warn(
                "LatticeIndex V_ee uses S³ chordal distance: "
                "V_ee = κ_ee / d²_S3,  κ_ee = 5Z/2 Ha,  d²_S3 = 2(1 − ξ_i·ξ_j). "
                "Same-orbital: antipodal d²=4, V_ee = 5Z/8 Ha exact. "
                "Off-diagonal d²(1s,2s)=0.4 → V_ee too large for multi-shell systems.",
                UserWarning,
                stacklevel=2,
            )
        elif vee_method == 'slater':
            warnings.warn(
                "LatticeIndex V_ee uses exact Slater F⁰ integrals computed "
                "numerically from hydrogen-like radial wavefunctions R_{nl}(r). "
                "Exact for all orbital pairs. Construction slower (numerical integration).",
                UserWarning,
                stacklevel=2,
            )
        elif vee_method == 'slater_full':
            warnings.warn(
                "LatticeIndex V_ee uses full Slater integrals (F^k direct + G^k "
                "exchange) with proper angular coupling via Wigner 3j symbols. "
                "Off-diagonal two-electron matrix elements enabled for FCI.",
                UserWarning,
                stacklevel=2,
            )
        elif vee_method == 's3_overlap':
            warnings.warn(
                "LatticeIndex V_ee uses S3 density-overlap formula: "
                "F0(a,b) = (4Z/pi) int Phi_a(t) Phi_b(t) dt (Paper 7 Sec VI). "
                "V_ee is a NODE property (density overlap), NOT an edge property. "
                "s-orbital pairs only; l>0 raises NotImplementedError.",
                UserWarning,
                stacklevel=2,
            )
        else:
            warnings.warn(
                "LatticeIndex V_ee uses Mulliken approximation: ~12% systematic "
                "overestimation of electron repulsion for well-separated orbitals; "
                "larger errors for same-orbital pairs (factor ~10x, V_ee=Z/n²). "
                "Li ground state energy will NOT be quantitatively accurate.",
                UserWarning,
                stacklevel=2,
            )

        self.n_electrons = n_electrons
        self.max_n = max_n
        self.Z = nuclear_charge
        self.threshold = sparsity_threshold
        self.vee_method = vee_method
        self.kinetic_scale = kinetic_scale
        self._adjacency_override = adjacency
        self._node_weights_override = node_weights
        # Default h1_method: 'exact' for slater_full/s3_overlap, 'hybrid' for slater, 'graph' otherwise
        if h1_method:
            self.h1_method = h1_method
        elif vee_method in ('slater_full', 's3_overlap'):
            self.h1_method = 'exact'
        elif vee_method == 'slater':
            self.h1_method = 'hybrid'
        else:
            self.h1_method = 'graph'

        self.fci_method = fci_method

        # --- Build topology in dependency order ---
        self._build_sp_lattice()
        self._compute_conformal_weights()
        self._build_sp_hamiltonian()
        self._build_vee_index()
        if enumerate_sds:
            self._enumerate_sd_basis()
        else:
            self.sd_basis = []
            self.n_sd = 0
            self._sd_index = {}

    # ------------------------------------------------------------------
    # Private construction methods
    # ------------------------------------------------------------------

    def _build_sp_lattice(self) -> None:
        """Build the single-particle geometric lattice and spin-orbital index."""
        self.lattice = GeometricLattice(max_n=self.max_n, nuclear_charge=self.Z)
        n_spatial = self.lattice.num_states

        # sp_states[2*i] = (n,l,m, 0=up)   sp_states[2*i+1] = (n,l,m, 1=down)
        self.sp_states: List[Tuple[int, int, int, int]] = []
        for n, l, m in self.lattice.states:
            self.sp_states.append((n, l, m, 0))
            self.sp_states.append((n, l, m, 1))
        self.n_sp = len(self.sp_states)

        self._sp_index: Dict[Tuple[int, int, int, int], int] = {
            s: i for i, s in enumerate(self.sp_states)
        }

    def _compute_conformal_weights(self) -> None:
        """
        Pre-compute conformal factor Omega_i for each spatial state.

        Omega_i = 2*p0 / (p_i^2 + p0^2)
        where p_i = Z / n_i   (state momentum on S^3)
              p0  = Z          (focal length = nuclear charge for atomic case)

        These are stored in LatticeIndex for future use as conformal edge
        modulation W_ij = W_flat * Omega_i * Omega_j. For the current
        single-atom assembly, they inform the physical interpretation but
        the graph edges use uniform weights (consistent with AtomicSolver).
        """
        p0 = float(self.Z)
        n_spatial = self.lattice.num_states
        self.conformal_weights = np.zeros(n_spatial)

        for idx, (n, l, m) in enumerate(self.lattice.states):
            p_i = float(self.Z) / float(n)
            self.conformal_weights[idx] = 2.0 * p0 / (p_i ** 2 + p0 ** 2)

    def _fock_s3_coords(self) -> np.ndarray:
        """
        Map each spatial state (n,l,m) to its S³ coordinate via Fock
        stereographic projection (Paper 7, Eq. 14).

        With normalised focal length p₀ = 1 and momentum magnitude p_i = 1/n_i
        (Z cancels identically), the S³ embedding is:

            ξ₁ = fac × sin(θ) × cos(φ)
            ξ₂ = fac × sin(θ) × sin(φ)
            ξ₃ = fac × cos(θ)
            ξ₄ = (1 − n²) / (n²+1)

        where  fac = 2n / (n²+1)  and the momentum direction (θ, φ) is:

            θ  encodes the z-projection:
                l=0  →  θ=0
                l>0  →  cos(θ) = m / √(l(l+1))

            φ  encodes the orbital angular momentum magnitude:
                l=0  →  φ=0  (irrelevant: sin(θ)=0)
                l>0  →  φ = 2π × l / n

        The azimuthal angle φ breaks the degeneracy between states with the
        same (n, m=0) but different l, which previously all mapped to
        θ=π/2, φ=0 → identical S³ coordinates.  Distributing l-subshells
        around the azimuthal circle ensures each state gets a unique S³ point.

        All ξ satisfy |ξ|=1.

        Returns
        -------
        coords : ndarray, shape (n_spatial, 4)
        """
        n_spatial = self.lattice.num_states
        coords = np.zeros((n_spatial, 4))
        for idx, (n, l, m) in enumerate(self.lattice.states):
            if l > 0:
                cos_th = np.clip(float(m) / np.sqrt(float(l * (l + 1))), -1.0, 1.0)
                sin_th = np.sqrt(max(0.0, 1.0 - cos_th * cos_th))
                phi = 2.0 * np.pi * l / n
            else:
                cos_th = 1.0
                sin_th = 0.0
                phi = 0.0
            fac = 2.0 * n / (n * n + 1.0)
            coords[idx, 0] = fac * sin_th * np.cos(phi)        # ξ₁
            coords[idx, 1] = fac * sin_th * np.sin(phi)        # ξ₂
            coords[idx, 2] = fac * cos_th                      # ξ₃
            coords[idx, 3] = (1.0 - n * n) / (1.0 + n * n)   # ξ₄
        return coords

    def _build_sp_hamiltonian(self) -> None:
        """
        Build single-particle Hamiltonian in spatial basis.

        Two modes:

        h1_method == 'graph' (default for vee_method='chordal'/'mulliken'):
            H1 = KINETIC_SCALE * (D - A) + diag(node_weights)
               = KINETIC_SCALE * L + V_en
            Off-diagonal elements from graph Laplacian adjacency.

        h1_method == 'exact' (default for vee_method='slater'):
            H1 = diag(-Z²/(2n²)) — exact hydrogen eigenbasis.
            States (n,l,m) ARE the eigenstates, so h1 is purely diagonal.
            No off-diagonal hopping: all excitations come from V_ee.

        Diagonal elements: h_diag[sp] = H1[sp, sp]
        Off-diagonal elements: h_offdiag[sp] = [(r_sp, H1[sp,r_sp]), ...]
        """
        n_spatial = self.lattice.num_states

        if self.h1_method == 'exact':
            # Exact hydrogen eigenbasis: h(n,l,m) = -Z²/(2n²), no off-diagonal
            Z = float(self.Z)
            self._h1_diag = np.array([
                -Z * Z / (2.0 * n * n)
                for n, l, m in self.lattice.states
            ])
            self._H1_spatial = diags(self._h1_diag).tocsr()
            self._h1_offdiag: Dict[int, List[Tuple[int, float]]] = {
                i: [] for i in range(n_spatial)
            }
        elif self.h1_method == 'hybrid':
            # Exact diagonal + graph Laplacian off-diagonal for CI mixing
            Z = float(self.Z)
            self._h1_diag = np.array([
                -Z * Z / (2.0 * n * n)
                for n, l, m in self.lattice.states
            ])

            # Off-diagonal from graph Laplacian adjacency (use override if provided)
            A = self._adjacency_override if self._adjacency_override is not None else self.lattice.adjacency
            H1_offdiag = self.kinetic_scale * (-A)
            self._H1_spatial = (diags(self._h1_diag) + H1_offdiag).tocsr()

            H1_coo = H1_offdiag.tocoo()
            self._h1_offdiag = {i: [] for i in range(n_spatial)}
            for r, c, v in zip(H1_coo.row, H1_coo.col, H1_coo.data):
                if r != c and abs(v) >= self.threshold:
                    self._h1_offdiag[r].append((c, float(v)))
        else:
            # Graph Laplacian mode (use overrides if provided)
            A = self._adjacency_override if self._adjacency_override is not None else self.lattice.adjacency
            deg = np.array(A.sum(axis=1)).flatten()
            L = diags(deg) - A
            nw = self._node_weights_override if self._node_weights_override is not None else self.lattice.node_weights
            V = diags(nw)
            H1 = self.kinetic_scale * L + V
            self._H1_spatial = H1.tocsr()

            self._h1_diag = np.array(H1.diagonal())

            H1_coo = H1.tocoo()
            self._h1_offdiag = {i: [] for i in range(n_spatial)}
            for r, c, v in zip(H1_coo.row, H1_coo.col, H1_coo.data):
                if r != c and abs(v) >= self.threshold:
                    self._h1_offdiag[r].append((c, float(v)))

    def _momentum_vectors(self) -> np.ndarray:
        """
        Assign a representative 3D momentum vector to each spatial state.

        For state (n, l, m):
            |p| = Z/n  (on-shell momentum for energy E_n = -Z²/(2n²))
            direction = (sin θ, 0, cos θ) with:
                l=0: θ=0 (z-axis)
                l>0: cos θ = m / √(l(l+1))

        Returns
        -------
        p_vecs : ndarray, shape (n_spatial, 3)
        """
        n_spatial = self.lattice.num_states
        p_vecs = np.zeros((n_spatial, 3))
        Z = float(self.Z)
        for idx, (n, l, m) in enumerate(self.lattice.states):
            p_mag = Z / float(n)
            if l > 0:
                cos_th = np.clip(float(m) / np.sqrt(float(l * (l + 1))), -1.0, 1.0)
                sin_th = np.sqrt(max(0.0, 1.0 - cos_th * cos_th))
            else:
                cos_th = 1.0
                sin_th = 0.0
            p_vecs[idx, 0] = p_mag * sin_th
            p_vecs[idx, 1] = 0.0
            p_vecs[idx, 2] = p_mag * cos_th
        return p_vecs

    def _conformal_factors_from_momenta(self, p_vecs: np.ndarray) -> np.ndarray:
        """
        Compute conformal factor Ω_k = 2p₀_k / (|p_k|² + p₀_k²) for each state.

        p₀_k = Z/n_k (energy-shell focal length for the kth state).
        |p_k| = Z/n_k (on-shell), so |p_k|² = p₀_k² = Z²/n_k².
        Therefore Ω_k = 2p₀ / (2p₀²) = 1/p₀ = n_k/Z.

        Returns
        -------
        omega : ndarray, shape (n_spatial,)
        """
        n_spatial = self.lattice.num_states
        omega = np.zeros(n_spatial)
        Z = float(self.Z)
        for idx, (n, l, m) in enumerate(self.lattice.states):
            p0 = Z / float(n)
            p_sq = float(np.dot(p_vecs[idx], p_vecs[idx]))
            omega[idx] = 2.0 * p0 / (p_sq + p0 * p0)
        return omega

    @staticmethod
    def _slater_f0_cache_path(Z: float, max_n: int) -> str:
        """Path to disk-cached Slater F0 integrals."""
        import os
        cache_dir = os.path.join(os.path.dirname(__file__), 'cache')
        os.makedirs(cache_dir, exist_ok=True)
        return os.path.join(cache_dir, f'slater_f0_Z{int(Z)}_nmax{max_n}.npz')

    def _compute_slater_f0(self) -> np.ndarray:
        """
        Compute exact Slater F0 direct Coulomb integrals for all spatial pairs.

        F0(n1 l1, n2 l2) = integral |R_{n1l1}(r1)|^2 |R_{n2l2}(r2)|^2
                            * (1/r_>) * r1^2 r2^2 dr1 dr2

        F0 depends only on (n, l), not on m, so we cache by (n, l) pairs.
        Results are cached to disk at geovac/cache/slater_f0_Z{Z}_nmax{N}.npz.

        Returns
        -------
        vee : ndarray, shape (n_spatial, n_spatial)
        """
        import os

        Z = float(self.Z)
        n_sp = self.lattice.num_states
        unique_nl = sorted(set((n, l) for n, l, m in self.lattice.states))

        # --- Try disk cache ---
        cache_path = self._slater_f0_cache_path(Z, self.max_n)
        f0_cache: Dict[Tuple[int, int, int, int], float] = {}

        if os.path.exists(cache_path):
            data = np.load(cache_path)
            keys = data['keys']     # shape (N, 4): n1, l1, n2, l2
            vals = data['values']   # shape (N,)
            for row, val in zip(keys, vals):
                f0_cache[(int(row[0]), int(row[1]),
                          int(row[2]), int(row[3]))] = float(val)

        # Check if all needed pairs are cached
        needed = set()
        for n1, l1 in unique_nl:
            for n2, l2 in unique_nl:
                key = (n1, l1, n2, l2)
                key_rev = (n2, l2, n1, l1)
                if key not in f0_cache and key_rev not in f0_cache:
                    needed.add(key if key <= key_rev else key_rev)

        if needed:
            # Compute missing integrals
            self._compute_f0_integrals(Z, needed, f0_cache, unique_nl)

            # Save to disk
            all_keys = []
            all_vals = []
            for k, v in f0_cache.items():
                all_keys.append(list(k))
                all_vals.append(v)
            np.savez(cache_path,
                     keys=np.array(all_keys, dtype=int),
                     values=np.array(all_vals))

        # Propagate symmetry
        for n1, l1 in unique_nl:
            for n2, l2 in unique_nl:
                key = (n1, l1, n2, l2)
                key_rev = (n2, l2, n1, l1)
                if key not in f0_cache and key_rev in f0_cache:
                    f0_cache[key] = f0_cache[key_rev]

        # --- Fill the V_ee matrix ---
        vee = np.zeros((n_sp, n_sp))
        for i in range(n_sp):
            ni, li, _ = self.lattice.states[i]
            for j in range(i, n_sp):
                nj, lj, _ = self.lattice.states[j]
                v = f0_cache[(ni, li, nj, lj)]
                vee[i, j] = v
                vee[j, i] = v
        return vee

    def _compute_s3_overlap(self) -> np.ndarray:
        """
        Compute V_ee via S3 density-overlap formula (Paper 7, Section VI).

        F0(a,b) = (4Z/pi) * integral_0^inf Phi_a(t) * Phi_b(t) dt

        where t = q/(2Z) is dimensionless momentum transfer and Phi_a(t) is the
        Fock-projected orbital density on S3.

        V_ee is a NODE property (density overlap on S3), NOT an edge property
        (pairwise chordal distance). The chordal ansatz kappa/d2 overestimates
        inter-shell integrals (e.g. 1s-2s) by ~30x.

        Currently supports s-orbital pairs only (l=0). For l>0 pairs, falls
        back to the position-space Slater F0 integral (_compute_f0_integrals).

        Returns
        -------
        vee : ndarray, shape (n_spatial, n_spatial)
        """
        Z = float(self.Z)
        n_sp = self.lattice.num_states
        vee = np.zeros((n_sp, n_sp))
        unique_nl = sorted(set((n, l) for n, l, m in self.lattice.states))

        # Pre-compute F0 for each unique (n1,l1)-(n2,l2) pair
        f0_cache: Dict[Tuple[int, int, int, int], float] = {}
        for n1, l1 in unique_nl:
            for n2, l2 in unique_nl:
                key = (n1, l1, n2, l2)
                key_rev = (n2, l2, n1, l1)
                if key in f0_cache or key_rev in f0_cache:
                    continue
                if l1 == 0 and l2 == 0:
                    # Use S3 density-overlap formula
                    val = compute_vee_s3_overlap(n1, l1, n2, l2, Z)
                else:
                    # l>0: fall back to position-space Slater F0
                    computed = self._compute_f0_integrals(Z, {key})
                    val = computed[key]
                f0_cache[key] = val
                f0_cache[key_rev] = val

        # Fill matrix
        for i in range(n_sp):
            ni, li, _ = self.lattice.states[i]
            for j in range(i, n_sp):
                nj, lj, _ = self.lattice.states[j]
                v = f0_cache[(ni, li, nj, lj)]
                vee[i, j] = v
                vee[j, i] = v
        return vee

    @staticmethod
    def _compute_f0_integrals(
        Z: float,
        needed: set,
        f0_cache: Dict[Tuple[int, int, int, int], float],
        unique_nl: list,
    ) -> None:
        """Compute Slater F0 integrals via nested numerical integration."""
        from scipy.special import genlaguerre
        from scipy.integrate import quad

        r_max = 80.0 / Z

        def _radial_wf_unnorm(r: np.ndarray, n: int, l: int) -> np.ndarray:
            rho = 2.0 * Z * r / n
            L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
            return rho ** l * np.exp(-rho / 2.0) * L_poly

        # Normalization constants
        norm_c: Dict[Tuple[int, int], float] = {}
        for n, l in unique_nl:
            norm_sq, _ = quad(
                lambda r: _radial_wf_unnorm(np.array(r), n, l) ** 2 * r ** 2,
                0, r_max, limit=200
            )
            norm_c[(n, l)] = 1.0 / np.sqrt(norm_sq)

        def R(r: float, n: int, l: int) -> float:
            return float(norm_c[(n, l)] * _radial_wf_unnorm(np.array(r), n, l))

        for n1, l1, n2, l2 in needed:
            def _inner(r1: float, _n2: int = n2, _l2: int = l2) -> float:
                if r1 < 1e-30:
                    return 0.0
                p1, _ = quad(lambda r2: R(r2, _n2, _l2) ** 2 * r2 ** 2,
                             0, r1, limit=100)
                p2, _ = quad(lambda r2: R(r2, _n2, _l2) ** 2 * r2,
                             r1, r_max, limit=100)
                return p1 / r1 + p2

            val, _ = quad(
                lambda r1, _n1=n1, _l1=l1: R(r1, _n1, _l1) ** 2 * _inner(r1) * r1 ** 2,
                0, r_max, limit=200
            )
            f0_cache[(n1, l1, n2, l2)] = val
            f0_cache[(n2, l2, n1, l1)] = val  # symmetry

    # ------------------------------------------------------------------
    # R^k Slater integrals and full two-electron machinery
    # ------------------------------------------------------------------

    @staticmethod
    def _rk_cache_path(Z: float, max_n: int) -> str:
        """Path to disk-cached R^k Slater integrals."""
        cache_dir = os.path.join(os.path.dirname(__file__), 'cache')
        os.makedirs(cache_dir, exist_ok=True)
        return os.path.join(cache_dir, f'slater_rk_Z{int(Z)}_nmax{max_n}.npz')

    @staticmethod
    def _wigner3j(j1: int, j2: int, j3: int,
                  m1: int, m2: int, m3: int) -> float:
        """
        Wigner 3j symbol for integer arguments.

        Uses the Racah formula. All arguments are integers (not half-integers).
        Returns 0 if triangle or m-selection rules are violated.
        """
        # Selection rules
        if m1 + m2 + m3 != 0:
            return 0.0
        if abs(j1 - j2) > j3 or j3 > j1 + j2:
            return 0.0
        if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
            return 0.0

        # Racah formula
        def _tri(a: int, b: int, c: int) -> float:
            return (factorial(a + b - c) * factorial(a - b + c)
                    * factorial(-a + b + c)
                    / factorial(a + b + c + 1))

        pre = ((-1) ** (j1 - j2 - m3)
               * np.sqrt(_tri(j1, j2, j3)
                         * factorial(j1 + m1) * factorial(j1 - m1)
                         * factorial(j2 + m2) * factorial(j2 - m2)
                         * factorial(j3 + m3) * factorial(j3 - m3)))

        t_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
        t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)

        s = 0.0
        for t in range(t_min, t_max + 1):
            s += ((-1) ** t
                  / (factorial(t)
                     * factorial(j3 - j2 + m1 + t)
                     * factorial(j3 - j1 - m2 + t)
                     * factorial(j1 + j2 - j3 - t)
                     * factorial(j1 - m1 - t)
                     * factorial(j2 + m2 - t)))
        return pre * s

    def _compute_rk_integrals(self) -> Dict[Tuple[int, ...], float]:
        """
        Compute R^k(n1l1, n2l2, n3l3, n4l4) Slater radial integrals.

        R^k(abcd) = ∫∫ R_a(r1) R_c(r1) (r_<^k / r_>^{k+1})
                       R_b(r2) R_d(r2) r1^2 r2^2 dr1 dr2

        These are the radial parts of the general two-electron integrals.
        F^k(ab) = R^k(ab,ab) (direct), G^k(ab) = R^k(ab,ba) (exchange).

        Cached to disk at geovac/cache/slater_rk_Z{Z}_nmax{N}.npz.

        Returns
        -------
        rk_cache : dict mapping (n1,l1,n2,l2,n3,l3,n4,l4,k) -> float
        """
        from scipy.special import genlaguerre
        from scipy.integrate import quad

        Z = float(self.Z)
        unique_nl = sorted(set((n, l) for n, l, m in self.lattice.states))
        r_max = 80.0 / Z

        # --- Try disk cache ---
        cache_path = self._rk_cache_path(Z, self.max_n)
        rk_cache: Dict[Tuple[int, ...], float] = {}

        if os.path.exists(cache_path):
            data = np.load(cache_path)
            keys = data['keys']    # shape (N, 9): n1,l1,n2,l2,n3,l3,n4,l4,k
            vals = data['values']
            for row, val in zip(keys, vals):
                rk_cache[tuple(int(x) for x in row)] = float(val)

        # Determine ALL needed R^k(n1l1,n2l2,n3l3,n4l4,k) tuples.
        # For full FCI Slater-Condon rules we need general R^k(a,b,c,d)
        # for ALL (nl) quadruples, not just F^k and G^k patterns.
        needed: set = set()
        for n1, l1 in unique_nl:
            for n2, l2 in unique_nl:
                for n3, l3 in unique_nl:
                    for n4, l4 in unique_nl:
                        k_max = min(l1 + l3, l2 + l4)
                        for k in range(0, k_max + 1):
                            # Parity selection: (l1+l3+k) even AND (l2+l4+k) even
                            if (l1 + l3 + k) % 2 != 0:
                                continue
                            if (l2 + l4 + k) % 2 != 0:
                                continue
                            key = (n1, l1, n2, l2, n3, l3, n4, l4, k)
                            if key not in rk_cache:
                                needed.add(key)

        if not needed:
            return rk_cache

        print(f"[LatticeIndex] Computing {len(needed)} R^k integrals "
              f"(grid method)...")
        t0 = time.perf_counter()

        # --- Grid-based integration (vectorized, ~100x faster than quad) ---
        N_GRID = 2000
        # Log-linear grid: denser near origin where wavefunctions peak
        r_grid = np.linspace(0, r_max, N_GRID + 1)[1:]  # exclude r=0
        dr = r_grid[1] - r_grid[0]

        # Pre-compute normalized radial wavefunctions on grid
        def _radial_wf_unnorm(r: np.ndarray, n: int, l: int) -> np.ndarray:
            rho = 2.0 * Z * r / n
            L_poly = genlaguerre(n - l - 1, 2 * l + 1)(rho)
            return rho ** l * np.exp(-rho / 2.0) * L_poly

        # Store R_{nl}(r) * r on grid for each (n,l)
        R_grid: Dict[Tuple[int, int], np.ndarray] = {}
        for n, l in unique_nl:
            wf = _radial_wf_unnorm(r_grid, n, l)
            norm_sq = np.trapezoid(wf ** 2 * r_grid ** 2, r_grid)
            R_grid[(n, l)] = wf / np.sqrt(norm_sq)

        # Pre-compute Y^k(r) = ∫_0^r f(r') (r'/r)^k r'^2 dr'
        #                     + ∫_r^∞ f(r') (r/r')^k r' dr'
        # where f(r') = R_b(r') * R_d(r')
        # This is the Yk potential for the (b,d) pair at multipole k.

        # Group needed integrals by (n2,l2,n4,l4,k) to reuse Yk
        from collections import defaultdict
        yk_groups: Dict[Tuple[int, ...], List[Tuple[int, ...]]] = defaultdict(list)
        for key in needed:
            n1, l1, n2, l2, n3, l3, n4, l4, k = key
            yk_key = (n2, l2, n4, l4, k)
            yk_groups[yk_key].append(key)

        # Compute Yk potentials
        yk_cache: Dict[Tuple[int, ...], np.ndarray] = {}
        for yk_key in yk_groups:
            n2, l2, n4, l4, k = yk_key
            f_r = R_grid[(n2, l2)] * R_grid[(n4, l4)] * r_grid ** 2
            yk = np.zeros(N_GRID)
            # Cumulative integral for inner part: ∫_0^r f(r') (r'/r)^k dr'
            for i in range(N_GRID):
                r1 = r_grid[i]
                # Inner: sum f(r') * (r'/r1)^k for r' <= r1
                inner = np.sum(f_r[:i + 1] * (r_grid[:i + 1] / r1) ** k) * dr
                # Outer: sum f(r') * (r1/r')^k * (1/r') for r' > r1
                if i + 1 < N_GRID:
                    outer = np.sum(
                        f_r[i + 1:]
                        * (r1 / r_grid[i + 1:]) ** k
                        / r_grid[i + 1:]
                    ) * dr
                else:
                    outer = 0.0
                yk[i] = inner / r1 + outer
            yk_cache[yk_key] = yk

        # Now compute R^k integrals using pre-computed Yk
        for yk_key, keys in yk_groups.items():
            n2, l2, n4, l4, k = yk_key
            yk = yk_cache[yk_key]
            for key in keys:
                n1, l1, _, _, n3, l3, _, _, _ = key
                integrand = (R_grid[(n1, l1)] * R_grid[(n3, l3)]
                             * yk * r_grid ** 2)
                val = np.trapezoid(integrand, r_grid)
                rk_cache[key] = val

        elapsed = time.perf_counter() - t0
        print(f"[LatticeIndex] R^k integrals computed in {elapsed:.1f}s")

        # Save to disk
        all_keys = [list(k) for k in rk_cache]
        all_vals = [rk_cache[k] for k in rk_cache]
        np.savez(cache_path,
                 keys=np.array(all_keys, dtype=int),
                 values=np.array(all_vals))

        return rk_cache

    def _build_two_electron_integrals(self) -> None:
        """
        Build two-electron integral table <ab|g|cd> using Slater integrals.

        <ab|g|cd> = Σ_k c^k(la,ma,lc,mc) c^k(lb,mb,ld,md) R^k(abcd)

        Optimization: pre-compute c^k coefficients, then iterate only over
        pairs with non-zero angular coupling. This avoids the naive O(n^4)
        loop over all spatial quadruples.

        Stores results in self._eri[a,b,c,d] where a,b,c,d are spatial indices.
        """
        rk_cache = self._compute_rk_integrals()
        n_sp = self.lattice.num_states
        states = self.lattice.states

        self._eri: Dict[Tuple[int, int, int, int], float] = {}

        # Pre-compute all c^k coefficients for pairs of spatial states
        # c^k(a,c) depends only on angular quantum numbers (l,m)
        ck_table: Dict[Tuple[int, int, int], float] = {}
        for a in range(n_sp):
            la, ma = states[a][1], states[a][2]
            for c in range(n_sp):
                lc, mc = states[c][1], states[c][2]
                k_max = la + lc
                for k in range(0, k_max + 1):
                    if (la + lc + k) % 2 != 0:
                        continue
                    val = self._ck_coefficient(la, ma, lc, mc, k)
                    if abs(val) > 1e-15:
                        ck_table[(a, c, k)] = val

        print(f"[LatticeIndex] c^k table: {len(ck_table)} non-zero entries")

        # Group c^k entries by (a,c) for efficient inner loop
        ac_k_map: Dict[Tuple[int, int], List[Tuple[int, float]]] = {}
        for (a, c, k), val in ck_table.items():
            key = (a, c)
            if key not in ac_k_map:
                ac_k_map[key] = []
            ac_k_map[key].append((k, val))

        count = 0
        for (a, c), ck_ac_list in ac_k_map.items():
            na, la, ma = states[a]
            nc, lc, mc = states[c]
            for (b, d), ck_bd_list in ac_k_map.items():
                nb, lb, mb = states[b]
                nd, ld, md = states[d]

                # m-selection rule: ma + mb = mc + md
                if ma + mb != mc + md:
                    continue

                val = 0.0
                for k_ac, c_ac in ck_ac_list:
                    for k_bd, c_bd in ck_bd_list:
                        if k_ac != k_bd:
                            continue
                        k = k_ac
                        rk_key = (na, la, nb, lb, nc, lc, nd, ld, k)
                        rk_val = rk_cache.get(rk_key)
                        if rk_val is None:
                            continue
                        val += c_ac * c_bd * rk_val

                if abs(val) > 1e-15:
                    self._eri[(a, b, c, d)] = val
                    count += 1

        print(f"[LatticeIndex] ERI table: {count} non-zero entries")

    def _ck_coefficient(self, la: int, ma: int, lc: int, mc: int, k: int) -> float:
        """
        Gaunt angular coupling coefficient.

        c^k(l,m,l',m') = (-1)^m √((2l+1)(2l'+1)) * (l k l'; 0 0 0) * (l k l'; -m q m')

        where q = mc - ma (m-transfer).
        """
        q = mc - ma
        pre = ((-1) ** ma
               * np.sqrt((2 * la + 1) * (2 * lc + 1)))
        w1 = self._wigner3j(la, k, lc, 0, 0, 0)
        if abs(w1) < 1e-15:
            return 0.0
        w2 = self._wigner3j(la, k, lc, -ma, q, mc)
        return pre * w1 * w2

    def _build_vee_index(self) -> None:
        """
        Pre-compute V_ee(i,j) for all spatial state pairs.

        Dispatches to the method selected by self.vee_method:

        'chordal' (default):
            V_ee(i,j) = κ_ee / d²_S3,  κ_ee = 5Z/2 Ha
            where d²_S3 = 2(1 − ξ_i · ξ_j) is the S³ chordal distance squared.
            Same-orbital pairs (i==j): antipodal override d² = 4,
                V_ee = κ_ee/4 = 5Z/8 Ha = F⁰(1s,1s) exactly.
            Good for He (0.97% error), but inter-shell d²(1s,2s)=0.4 is too
            small → catastrophic V_ee for multi-shell systems like Li.

        'slater' (exact):
            V_ee(i,j) = F⁰(n_i l_i, n_j l_j) — exact Slater direct Coulomb
            integrals computed numerically from hydrogen-like radial wavefunctions.
            F⁰ = ∫∫ |R_{nl}(r₁)|² |R_{n'l'}(r₂)|² (1/r_>) r₁² r₂² dr₁ dr₂.

        's3_overlap' (Paper 7, Sec VI):
            F⁰(a,b) = (4Z/π) ∫₀^∞ Φ_a(t)·Φ_b(t) dt, where Φ_a is the
            Fock-projected orbital density on S³. NODE property, not edge.
            s-orbital pairs via closed-form densities; l>0 falls back to
            position-space Slater F⁰.

        'mulliken' (legacy):
            V_ee(i,j) = 1 / r_ij,  r_ij = |r_i − r_j| with r_i = n_i²/Z.
            Same-orbital: r_eff = n²/Z  →  V_ee = Z/n².
            Known ~12% overestimation for well-separated pairs; ~10x for same n=1.
        """
        n_sp = self.lattice.num_states
        self._vee_matrix: np.ndarray = np.zeros((n_sp, n_sp))

        if self.vee_method == 'chordal':
            KAP_EE: float = 5.0 * float(self.Z) / 2.0
            ANTIPODAL_D_SQ: float = 4.0

            s3 = self._fock_s3_coords()           # (n_sp, 4)
            gram = s3 @ s3.T                       # dot products
            d_sq = 2.0 * (1.0 - gram)              # chordal distance²

            np.fill_diagonal(d_sq, ANTIPODAL_D_SQ)
            np.clip(d_sq, 1e-14, None, out=d_sq)
            self._vee_matrix = KAP_EE / d_sq

        elif self.vee_method == 'slater':
            self._vee_matrix = self._compute_slater_f0()

        elif self.vee_method == 'slater_full':
            # Full Slater integrals: F0 for diagonal V_ee + ERI table for
            # off-diagonal two-electron matrix elements
            self._vee_matrix = self._compute_slater_f0()
            self._build_two_electron_integrals()

        elif self.vee_method == 's3_overlap':
            self._vee_matrix = self._compute_s3_overlap()

        else:  # 'mulliken'
            coords = self._compute_sp_coordinates()
            for i in range(n_sp):
                for j in range(i, n_sp):
                    if i == j:
                        n_i = self.lattice.states[i][0]
                        v = float(self.Z) / float(n_i ** 2)
                    else:
                        rij = float(np.linalg.norm(coords[i] - coords[j]))
                        if rij < 1e-10:
                            n_i = self.lattice.states[i][0]
                            v = float(self.Z) / float(n_i ** 2)
                        else:
                            v = 1.0 / rij
                    self._vee_matrix[i, j] = v
                    self._vee_matrix[j, i] = v

    def _compute_sp_coordinates(self) -> np.ndarray:
        """
        Compute 3D spatial coordinates for each spatial state.
        r = n^2/Z (effective Bohr radius), orbital aligned on z-axis.
        """
        n_sp = self.lattice.num_states
        coords = np.zeros((n_sp, 3))
        for idx, (n, l, m) in enumerate(self.lattice.states):
            r = n ** 2 / float(self.Z)
            if l > 0:
                cos_th = np.clip(m / np.sqrt(l * (l + 1)), -1.0, 1.0)
                theta = np.arccos(cos_th)
            else:
                theta = 0.0
            coords[idx] = [r * np.sin(theta), 0.0, r * np.cos(theta)]
        return coords

    def _enumerate_sd_basis(self) -> None:
        """
        Enumerate all valid Slater determinants.

        Each SD is a sorted tuple of n_electrons distinct spin-orbital indices.
        Pauli exclusion is automatic: combinations() never repeats indices.
        Total: C(n_sp, n_electrons) determinants.
        """
        t0 = time.perf_counter()
        self.sd_basis: List[Tuple[int, ...]] = list(
            combinations(range(self.n_sp), self.n_electrons)
        )
        self.n_sd: int = len(self.sd_basis)
        self._sd_index: Dict[Tuple[int, ...], int] = {
            sd: i for i, sd in enumerate(self.sd_basis)
        }
        elapsed = time.perf_counter() - t0
        n_spatial = self.lattice.num_states
        print(
            f"[LatticeIndex] max_n={self.max_n}: {n_spatial} spatial states, "
            f"{self.n_sp} spin-orbitals, {self.n_sd:,} Slater determinants "
            f"(enumerated in {elapsed:.3f}s)"
        )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def _get_eri(self, a: int, b: int, c: int, d: int) -> float:
        """Look up two-electron integral <ab|cd> from ERI table."""
        return self._eri.get((a, b, c, d), 0.0)

    def _slater_condon_diagonal(self, sd: Tuple[int, ...]) -> float:
        """
        Slater-Condon diagonal: <I|H|I> with full two-electron integrals.

        = Σ_p h1(p,p) + Σ_{p<q} [<pq|pq> - δ(σ_p,σ_q)<pq|qp>]
        """
        h_diag = 0.0
        for p in sd:
            h_diag += self._h1_diag[p >> 1]

        n = len(sd)
        for i in range(n):
            pi = sd[i]
            sp_i = pi >> 1
            sig_i = pi & 1
            for j in range(i + 1, n):
                pj = sd[j]
                sp_j = pj >> 1
                sig_j = pj & 1
                # Coulomb: <ij|ij>
                h_diag += self._get_eri(sp_i, sp_j, sp_i, sp_j)
                # Exchange: -δ(σi,σj) <ij|ji>
                if sig_i == sig_j:
                    h_diag -= self._get_eri(sp_i, sp_j, sp_j, sp_i)
        return h_diag

    def _slater_condon_single(self, sd: Tuple[int, ...],
                              kp: int, r: int) -> float:
        """
        Slater-Condon single excitation: sd[kp] → r.

        = h1(p,r) + Σ_{q in occ, q≠p} [<pq|rq> - δ(σ_p,σ_q)<pq|qr>]

        Returns the UNSIGNED matrix element (caller applies phase).
        """
        p = sd[kp]
        sp_p = p >> 1
        sp_r = r >> 1
        sig_p = p & 1
        sig_r = r & 1

        # One-body part
        val = self._H1_spatial[sp_p, sp_r]

        # Two-body part
        for q in sd:
            if q == p:
                continue
            sp_q = q >> 1
            sig_q = q & 1
            # Coulomb: <pq|rq>
            val += self._get_eri(sp_p, sp_q, sp_r, sp_q)
            # Exchange: -δ(σ_p,σ_q) <pq|qr>
            if sig_p == sig_q:
                val -= self._get_eri(sp_p, sp_q, sp_q, sp_r)

        return val

    def assemble_hamiltonian(self) -> csr_matrix:
        """
        Assemble the N-electron Hamiltonian via Slater-Condon rules.

        For vee_method='slater_full', uses full two-electron integrals:
          Diagonal:  Σ_p h1(p,p) + Σ_{p<q} [<pq|pq> - δ(σ)·<pq|qp>]
          Single exc (p→r): h1(p,r) + Σ_q [<pq|rq> - δ(σ)·<pq|qr>]
          Double exc (pq→rs): <pq|rs> - δ(σ)·<pq|sr>

        For other vee_methods, uses diagonal-only V_ee (original behavior).

        Returns
        -------
        H : csr_matrix, shape (n_sd, n_sd)
        """
        if self.vee_method == 'slater_full':
            return self._assemble_hamiltonian_full()
        return self._assemble_hamiltonian_diag_vee()

    def _assemble_hamiltonian_diag_vee(self) -> csr_matrix:
        """Original assembly: diagonal V_ee only, off-diagonal from h1."""
        t0 = time.perf_counter()

        diag_rows: List[int] = []
        diag_vals: List[float] = []
        off_rows: List[int] = []
        off_cols: List[int] = []
        off_vals: List[float] = []

        for I, sd in enumerate(self.sd_basis):
            occupied_set = set(sd)

            # ---- Diagonal: one-body sum + V_ee ----
            h_diag = 0.0
            for p in sd:
                sp = p >> 1
                h_diag += self._h1_diag[sp]
            h_diag += self._compute_sd_vee(sd)

            if abs(h_diag) >= self.threshold:
                diag_rows.append(I)
                diag_vals.append(h_diag)

            # ---- Off-diagonal: single excitations (h1 only) ----
            for kp, p in enumerate(sd):
                sigma_p = p & 1
                sp_p = p >> 1

                for r_sp, h_val in self._h1_offdiag[sp_p]:
                    if abs(h_val) < self.threshold:
                        continue
                    r = (r_sp << 1) | sigma_p
                    if r in occupied_set:
                        continue

                    new_sd = tuple(sorted(
                        sd[:kp] + (r,) + sd[kp + 1:]
                    ))
                    if new_sd not in self._sd_index:
                        continue
                    J = self._sd_index[new_sd]
                    if J <= I:
                        continue

                    phase = self._compute_phase(sd, kp, r)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * h_val)

        H_diag = csr_matrix(
            (diag_vals, (diag_rows, diag_rows)),
            shape=(self.n_sd, self.n_sd),
        )
        H_upper = csr_matrix(
            (off_vals, (off_rows, off_cols)),
            shape=(self.n_sd, self.n_sd),
        )
        H = H_upper + H_upper.T + H_diag

        elapsed = time.perf_counter() - t0
        print(
            f"[LatticeIndex] Hamiltonian assembled: shape={H.shape}, "
            f"nnz={H.nnz:,}, time={elapsed:.3f}s"
        )
        return H.tocsr()

    def _assemble_hamiltonian_full(self) -> csr_matrix:
        """
        Full Slater-Condon assembly with off-diagonal two-electron integrals.

        Handles diagonal, single excitations, and double excitations.
        Optimized: iterates over SD pairs (I,J) with J>I and classifies
        excitation order by set difference, avoiding the O(n_sp^2) inner
        loop for double excitations.
        """
        t0 = time.perf_counter()

        diag_rows: List[int] = []
        diag_vals: List[float] = []
        off_rows: List[int] = []
        off_cols: List[int] = []
        off_vals: List[float] = []

        n_el = self.n_electrons

        # --- Phase 1: Diagonal ---
        for I, sd_I in enumerate(self.sd_basis):
            h_diag = self._slater_condon_diagonal(sd_I)
            if abs(h_diag) >= self.threshold:
                diag_rows.append(I)
                diag_vals.append(h_diag)

        t_diag = time.perf_counter() - t0

        # --- Phase 2: Off-diagonal (single + double excitations) ---
        # Pre-convert SDs to frozensets for fast set operations
        sd_sets = [frozenset(sd) for sd in self.sd_basis]

        for I in range(self.n_sd):
            sd_I = self.sd_basis[I]
            set_I = sd_sets[I]

            for J in range(I + 1, self.n_sd):
                set_J = sd_sets[J]

                # Determine excitation order
                diff_I = set_I - set_J   # orbitals in I but not J (removed)
                diff_J = set_J - set_I   # orbitals in J but not I (added)
                n_diff = len(diff_I)

                if n_diff == 0:
                    continue  # identical SDs
                elif n_diff > 2:
                    continue  # triple+ excitation: zero by Slater-Condon

                elif n_diff == 1:
                    # Single excitation: p → r
                    p = next(iter(diff_I))
                    r = next(iter(diff_J))

                    # Spin conservation
                    if (p & 1) != (r & 1):
                        continue

                    # Find kp (position of p in sd_I)
                    kp = sd_I.index(p)

                    me = self._slater_condon_single(sd_I, kp, r)
                    if abs(me) < self.threshold:
                        continue

                    phase = self._compute_phase(sd_I, kp, r)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * me)

                elif n_diff == 2:
                    # Double excitation: p,q → r,s
                    removed = sorted(diff_I)
                    added = sorted(diff_J)
                    p, q = removed[0], removed[1]
                    r, s = added[0], added[1]

                    sig_p, sig_q = p & 1, q & 1
                    sig_r, sig_s = r & 1, s & 1
                    sp_p, sp_q = p >> 1, q >> 1
                    sp_r, sp_s = r >> 1, s >> 1

                    # Spin conservation
                    if sig_p + sig_q != sig_r + sig_s:
                        continue

                    # <pq|rs> with proper spin matching
                    me = 0.0
                    if sig_p == sig_r and sig_q == sig_s:
                        me += self._get_eri(sp_p, sp_q, sp_r, sp_s)
                    if sig_p == sig_s and sig_q == sig_r:
                        me -= self._get_eri(sp_p, sp_q, sp_s, sp_r)

                    if abs(me) < self.threshold:
                        continue

                    kp = sd_I.index(p)
                    kq = sd_I.index(q)
                    phase = self._compute_double_phase(sd_I, kp, kq, r, s)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * me)

        H_diag = csr_matrix(
            (diag_vals, (diag_rows, diag_rows)),
            shape=(self.n_sd, self.n_sd),
        )
        H_upper = csr_matrix(
            (off_vals, (off_rows, off_cols)),
            shape=(self.n_sd, self.n_sd),
        )
        H = H_upper + H_upper.T + H_diag

        elapsed = time.perf_counter() - t0
        print(
            f"[LatticeIndex] Full Hamiltonian assembled: shape={H.shape}, "
            f"nnz={H.nnz:,}, time={elapsed:.3f}s "
            f"(diag={t_diag:.1f}s)"
        )
        return H.tocsr()

    def compute_ground_state(
        self, n_states: int = 1
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the N-electron ground state via sparse diagonalisation.

        Uses the method selected by ``fci_method``:
        - ``'auto'``: direct CI for N_SD >= 5000, matrix otherwise
        - ``'direct'``: excitation-driven assembly (O(N_SD x N_connected))
        - ``'matrix'``: explicit pairwise assembly (O(N^2_SD), reference)

        Returns
        -------
        eigvals : np.ndarray, shape (n_states,)
        eigvecs : np.ndarray, shape (n_sd, n_states)
        """
        method = self.fci_method
        if method == 'auto':
            method = 'direct' if self.n_sd >= 5000 else 'matrix'

        if method == 'direct' and self.vee_method == 'slater_full':
            from .direct_ci import DirectCISolver
            solver = DirectCISolver(self)
            return solver.solve(n_states=n_states)

        # Fall back to matrix method (original implementation)
        H = self.assemble_hamiltonian()
        k = min(n_states, self.n_sd - 2)
        # Use deterministic initial vector for reproducible ARPACK convergence
        rng = np.random.RandomState(42)
        v0 = rng.randn(H.shape[0])
        eigvals, eigvecs = eigsh(H, k=k, which="SA", v0=v0)
        order = np.argsort(eigvals)
        return eigvals[order], eigvecs[:, order]

    # ------------------------------------------------------------------
    # Helper methods
    # ------------------------------------------------------------------

    def _compute_sd_vee(self, sd: Tuple[int, ...]) -> float:
        """Sum V_ee(sp_i, sp_j) for all pairs i < j in the Slater determinant."""
        total = 0.0
        n = len(sd)
        for i in range(n):
            for j in range(i + 1, n):
                total += self._vee_matrix[sd[i] >> 1, sd[j] >> 1]
        return total

    @staticmethod
    def _compute_phase(sd: Tuple[int, ...], kp: int, r: int) -> float:
        """
        Fermionic sign for the single excitation sd[kp] → r.

        Phase = (-1)^{kp + kr}
          kp = index of removed orbital in original SD (0-based)
          kr = #{a in SD : a < r, a ≠ sd[kp]}
        """
        p = sd[kp]
        kr = sum(1 for a in sd if a < r and a != p)
        return (-1.0) ** (kp + kr)

    @staticmethod
    def _compute_double_phase(sd: Tuple[int, ...],
                              kp: int, kq: int,
                              r: int, s: int) -> float:
        """
        Fermionic sign for double excitation sd[kp],sd[kq] → r,s.

        Computed by applying two sequential single excitations and
        combining phases.
        """
        p = sd[kp]
        q = sd[kq]

        # First excitation: remove p, insert r
        # Count swaps to move p out
        n_swap_p = kp
        # Count position for r among remaining (after removing p)
        remaining_1 = [a for a in sd if a != p]
        kr = sum(1 for a in remaining_1 if a < r)
        phase1 = (-1) ** (n_swap_p + kr)

        # Second excitation: remove q, insert s from the intermediate SD
        intermediate = sorted(remaining_1[:kr] + [r] + remaining_1[kr:])
        # Find position of q in intermediate
        kq_new = intermediate.index(q)
        remaining_2 = [a for a in intermediate if a != q]
        ks = sum(1 for a in remaining_2 if a < s)
        phase2 = (-1) ** (kq_new + ks)

        return float(phase1 * phase2)

    def assembly_stats(self) -> dict:
        """
        Return statistics about the topology without assembling the full H.
        Useful for scaling benchmarks.
        """
        nnz_offdiag = 0
        for sp, neighbors in self._h1_offdiag.items():
            nnz_offdiag += len(neighbors)

        return {
            "max_n": self.max_n,
            "n_spatial": self.lattice.num_states,
            "n_sp": self.n_sp,
            "n_sd": self.n_sd,
            "h1_nnz": int(self._H1_spatial.nnz),
            "h1_offdiag_edges": nnz_offdiag,
            "vee_method": self.vee_method,
            "h1_method": self.h1_method,
        }


# ===========================================================================
# MolecularLatticeIndex — Heteronuclear diatomic FCI
# ===========================================================================


class MolecularLatticeIndex:
    """
    FCI solver for diatomic molecules using two atomic hydrogenic lattices.

    Constructs the molecular Hamiltonian as:
        H_mol = H1_mol + V_ee(same-atom + cross-atom s-orbital) + V_NN

    The one-electron Hamiltonian combines:
    - Exact atomic eigenvalues -Z²/(2n²) on each atom (diagonal)
    - Cross-nuclear attraction via electrostatic potential (s-orbitals only)
    - Inter-atomic kinetic hopping via graph bridges (off-diagonal)

    V_ee combines exact same-atom Slater integrals with cross-atom direct
    Coulomb J_AB via Fourier convolution of S³ momentum-space densities.
    Only s-orbital cross-atom pairs are included (l>0 requires angular-
    dependent form factors for variational rigour).

    Parameters
    ----------
    Z_A, Z_B : int
        Nuclear charges of atoms A and B
    nmax_A, nmax_B : int
        Basis size for each atom (can differ)
    R : float
        Internuclear distance in Bohr
    n_electrons : int
        Total electron count (= Z_A + Z_B for neutral molecule)
    n_bridges : int
        Number of inter-atomic bridge connections (default 20)
    vee_method : str
        Two-electron integral method ('slater_full' recommended)
    fci_method : str
        FCI assembly method ('auto', 'direct', or 'matrix')
    kinetic_scale : float
        Universal kinetic scale factor (default -1/16)
    cross_atom_vee : bool or str
        Cross-atom V_ee mode: True = all (n,l) pairs (default),
        's_only' = s-orbital pairs only (v0.9.11), False = disabled.
    cross_nuclear_method : str
        Cross-nuclear attraction method: 'exact' = full 2D quadrature
        for all (n,l,m) orbitals (default), 'fourier' = shell-theorem
        formula for s-orbitals only (v0.9.35 baseline).
    zeta_A : float
        Orbital exponent scale for atom A (default 1.0). Each orbital on
        atom A uses effective charge Z_eff = zeta_A * Z_A for its radial
        shape. All integrals (H1 diagonal, Slater F0, cross-nuclear,
        cross-atom J/K) are recomputed with scaled exponents.
    zeta_B : float
        Orbital exponent scale for atom B (default 1.0). Same as zeta_A
        but for atom B orbitals.
    overlap_edges : str, optional
        Add R-dependent cross-atom edges to the adjacency matrix
        proportional to orbital overlap. Makes the graph topology
        (and hence kinetic energy via the degree matrix) R-dependent.
        Modes: 's2' = S², 'fock' = S²×(Z/n)⁴, 'abs' = |S|,
        'kappa' = S²×|κ|. Default: None (disabled).
    overlap_edge_scale : float
        Multiplicative scale for overlap edge weights (default 1.0).
    """

    # Reuse LatticeIndex static methods for fermionic phase computation
    _compute_phase = staticmethod(LatticeIndex._compute_phase)
    _compute_double_phase = staticmethod(LatticeIndex._compute_double_phase)

    def __init__(
        self,
        Z_A: int,
        Z_B: int,
        nmax_A: int,
        nmax_B: int,
        R: float,
        n_electrons: int,
        n_bridges: int = 20,
        vee_method: str = 'slater_full',
        fci_method: str = 'auto',
        kinetic_scale: float = KINETIC_SCALE,
        sparsity_threshold: float = 1e-8,
        orthogonalize: bool = False,
        use_dmatrix: object = False,
        use_sturmian: object = False,
        sturmian_p0: Optional[float] = None,
        cross_atom_vee: object = True,
        cross_nuclear_method: str = 'exact',
        zeta_A: float = 1.0,
        zeta_B: float = 1.0,
        t_corr_lambda: float = 0.0,
        t_corr_fock_weighted: bool = False,
        overlap_edges: Optional[str] = None,
        overlap_edge_scale: float = 1.0,
        cross_nuclear_attenuation: Optional[str] = None,
        cross_nuclear_alpha: float = 1.0,
        cross_nuclear_beta: Optional[float] = None,
        use_cross_n_bridges: bool = False,
        cross_n_K: float = 1.75,
        enumerate_sds: bool = True,
    ) -> None:
        self.Z_A = Z_A
        self.Z_B = Z_B
        self.nmax_A = nmax_A
        self.nmax_B = nmax_B
        self.R = R
        self.n_electrons = n_electrons
        self.n_bridges = n_bridges
        self.vee_method = vee_method
        self.fci_method = fci_method
        self.kinetic_scale = kinetic_scale
        self.threshold = sparsity_threshold
        self.orthogonalize = orthogonalize
        self.use_dmatrix = use_dmatrix
        self.use_sturmian = use_sturmian
        self.cross_atom_vee = cross_atom_vee
        self.cross_nuclear_method = cross_nuclear_method
        self.zeta_A = zeta_A
        self.zeta_B = zeta_B
        self.t_corr_lambda = t_corr_lambda
        self.t_corr_fock_weighted = t_corr_fock_weighted
        self.overlap_edges = overlap_edges
        self.overlap_edge_scale = overlap_edge_scale
        self.cross_nuclear_attenuation = cross_nuclear_attenuation
        self.cross_nuclear_alpha = cross_nuclear_alpha
        self.cross_nuclear_beta = cross_nuclear_beta
        self.use_cross_n_bridges = use_cross_n_bridges
        self.cross_n_K = cross_n_K
        # Effective orbital charges: zeta_scale * Z determines orbital shape
        self._Z_orb_A = zeta_A * float(Z_A) if Z_A > 0 else 0.0
        self._Z_orb_B = zeta_B * float(Z_B) if Z_B > 0 else 0.0

        # Sturmian p0: shared momentum scale (Paper 9)
        if use_sturmian == 'atomic':
            # Atom-dependent p0 (v0.9.22): each atom uses its own
            # self-consistent p0 from isolated-atom FCI.
            self._sturmian_p0_A = compute_atomic_p0(Z_A, nmax_A)
            self._sturmian_p0_B = compute_atomic_p0(Z_B, nmax_B)
            self._sturmian_p0_AB = np.sqrt(
                self._sturmian_p0_A * self._sturmian_p0_B
            )
            # Store a nominal p0 for compatibility
            self._sturmian_p0 = self._sturmian_p0_AB
            self.use_dmatrix = 'hybrid'
        elif use_sturmian == 'molecular':
            # MO-projected Sturmian (v0.9.30): prolate spheroidal MO betas
            # projected onto atom-centered basis via dominant-overlap.
            if sturmian_p0 is not None:
                self._sturmian_p0 = sturmian_p0
            else:
                self._sturmian_p0 = np.sqrt(
                    (float(Z_A)**2 + float(Z_B)**2) / 2.0
                )
            self.use_dmatrix = 'hybrid'

        elif use_sturmian:
            if sturmian_p0 is not None:
                self._sturmian_p0 = sturmian_p0
            else:
                # Default: geometric mean of atomic scales (Paper 9, Sec. VI.B)
                self._sturmian_p0 = np.sqrt(
                    (float(Z_A)**2 + float(Z_B)**2) / 2.0
                )
            # Force hybrid-like D-matrix path for Sturmian
            self.use_dmatrix = 'hybrid'

        # Ghost atom support: Z=0 means basis functions are present
        # but carry no nuclear attraction (Boys-Bernardi counterpoise).
        self._ghost_A = (Z_A == 0)
        self._ghost_B = (Z_B == 0)

        # --- Build atomic lattice indices for ERI computation ---
        # For ghost atoms (Z=0), we still need the lattice topology (states,
        # adjacency) but set Z=1 internally so GeometricLattice doesn't fail.
        # The h1_diag will be overwritten to 0 in _build_molecular_h1.
        Z_A_eff = Z_A if Z_A > 0 else 1
        Z_B_eff = Z_B if Z_B > 0 else 1
        label = "ghost" if (self._ghost_A or self._ghost_B) else "LiH"
        zeta_str = ""
        if abs(zeta_A - 1.0) > 1e-12 or abs(zeta_B - 1.0) > 1e-12:
            zeta_str = f", zeta_A={zeta_A:.4f}, zeta_B={zeta_B:.4f}"
        print(f"[MolecularLatticeIndex] {label}: Z_A={Z_A}, Z_B={Z_B}, "
              f"nmax_A={nmax_A}, nmax_B={nmax_B}, R={R:.2f} bohr, "
              f"Ne={n_electrons}{zeta_str}")

        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self._li_A = LatticeIndex(
                n_electrons=1, max_n=nmax_A, nuclear_charge=Z_A_eff,
                vee_method=vee_method, h1_method='exact',
                fci_method='matrix', kinetic_scale=kinetic_scale,
            )
            self._li_B = LatticeIndex(
                n_electrons=1, max_n=nmax_B, nuclear_charge=Z_B_eff,
                vee_method=vee_method, h1_method='exact',
                fci_method='matrix', kinetic_scale=kinetic_scale,
            )

        self._n_spatial_A = self._li_A.lattice.num_states
        self._n_spatial_B = self._li_B.lattice.num_states
        self._n_spatial = self._n_spatial_A + self._n_spatial_B
        self.n_sp = 2 * self._n_spatial

        # --- Build molecular data structures ---
        self._build_combined_states()
        self._build_combined_adjacency()
        self._build_molecular_h1()
        self._build_molecular_vee()

        # --- Optional Lowdin orthogonalization ---
        if self.orthogonalize:
            S = self._compute_overlap_matrix()
            self._lowdin_orthogonalize(S)

        if enumerate_sds:
            self._enumerate_sd_basis()
        else:
            # Skip SD enumeration (used by FrozenCoreLatticeIndex)
            self.sd_basis = []
            self.n_sd = 0
            self._sd_index = {}

        # Nuclear repulsion (exact; zero for ghost atoms)
        self.V_NN = float(Z_A * Z_B) / R if (Z_A > 0 and Z_B > 0) else 0.0
        print(f"[MolecularLatticeIndex] V_NN = {self.V_NN:.6f} Ha")

    # ------------------------------------------------------------------
    # Construction methods
    # ------------------------------------------------------------------

    def _build_combined_states(self) -> None:
        """Build combined spin-orbital state list from both atoms."""
        self.sp_states: List[Tuple[int, int, int, int]] = []
        self._spatial_atom: List[int] = []  # which atom each spatial state belongs to

        for n, l, m in self._li_A.lattice.states:
            self.sp_states.append((n, l, m, 0))
            self.sp_states.append((n, l, m, 1))
            self._spatial_atom.append(0)

        for n, l, m in self._li_B.lattice.states:
            self.sp_states.append((n, l, m, 0))
            self.sp_states.append((n, l, m, 1))
            self._spatial_atom.append(1)

    def _build_combined_adjacency(self) -> None:
        """
        Build combined adjacency matrix: block-diagonal intra-atom edges
        plus inter-atomic bridges with STO overlap and conformal weighting.
        """
        from scipy.sparse import block_diag as sp_block_diag, coo_matrix

        A_A = self._li_A.lattice.adjacency
        A_B = self._li_B.lattice.adjacency
        nA = self._n_spatial_A
        N = self._n_spatial

        # Block-diagonal intra-atom adjacency
        block = sp_block_diag([A_A, A_B], format='coo')
        rows = list(block.row)
        cols = list(block.col)
        data = list(block.data.astype(np.float64))

        # Inter-atomic bridges (replicates MoleculeHamiltonian bridge logic)
        lat_A = self._li_A.lattice
        lat_B = self._li_B.lattice

        boundary_A = lat_A._get_boundary_states_prioritized()
        boundary_B = lat_B._get_boundary_states_prioritized()
        n_actual = 0  # reset; only set if bridges are actually created

        # Skip bridges if either atom is a ghost (Z=0): no physical hopping
        n_candidate = min(len(boundary_A), len(boundary_B), self.n_bridges)
        if n_candidate > 0 and self.R > 1e-10 and not self._ghost_A and not self._ghost_B:
            R_AB = self.R
            S_R = (1.0 + R_AB + R_AB**2 / 3.0) * np.exp(-R_AB)

            V_nn = float(self.Z_A * self.Z_B) / R_AB
            p0_A = np.sqrt(float(self.Z_A)**2 + V_nn)
            p0_B = np.sqrt(float(self.Z_B)**2 + V_nn)

            local_A = np.array(boundary_A[:n_candidate], dtype=np.intp)
            local_B = np.array(boundary_B[:n_candidate], dtype=np.intp)

            n_vals_A = np.array([lat_A.states[k][0] for k in local_A],
                                dtype=np.float64)
            n_vals_B = np.array([lat_B.states[k][0] for k in local_B],
                                dtype=np.float64)

            p_A = float(self.Z_A) / n_vals_A
            p_B = float(self.Z_B) / n_vals_B
            omega_A = 2.0 * p0_A / (p_A**2 + p0_A**2)
            omega_B = 2.0 * p0_B / (p_B**2 + p0_B**2)

            weights = S_R * omega_A * omega_B

            mask = np.abs(weights) >= 1e-8
            local_A = local_A[mask]
            local_B_offset = local_B[mask] + nA
            weights = weights[mask]
            n_actual = int(mask.sum())

            # Symmetric bridges
            rows.extend(local_A.tolist())
            rows.extend(local_B_offset.tolist())
            cols.extend(local_B_offset.tolist())
            cols.extend(local_A.tolist())
            data.extend(weights.tolist())
            data.extend(weights.tolist())

        # --- Overlap-induced edges ---
        # Add cross-atom edges proportional to orbital overlap.
        # This makes the degree matrix (kinetic energy) R-dependent:
        # as atoms approach, overlap grows → degree increases → kinetic
        # energy rises → creates repulsive wall.
        n_overlap_edges = 0
        if self.overlap_edges and not self._ghost_A and not self._ghost_B:
            n_overlap_edges = self._add_overlap_edges(
                rows, cols, data, nA, N)

        self._adjacency_combined = coo_matrix(
            (data, (rows, cols)), shape=(N, N)
        ).tocsr()
        self._n_bridges_actual = n_actual
        msg = (f"[MolecularLatticeIndex] Bridges: {n_actual} "
               f"(requested {self.n_bridges})")
        if n_overlap_edges > 0:
            msg += f", overlap edges: {n_overlap_edges} (mode={self.overlap_edges})"
        print(msg)

    def _add_overlap_edges(
        self,
        rows: list, cols: list, data: list,
        nA: int, N: int,
    ) -> int:
        """
        Add overlap-induced cross-atom edges to the adjacency matrix.

        These edges make the graph topology R-dependent: as atoms
        approach, overlap grows, degree increases, kinetic energy rises.

        Edge formulas (selected by self.overlap_edges):
          's2'   : A[a,b] = scale * S²(a,b)
          'fock' : A[a,b] = scale * S²(a,b) * (Z_A/n_a)² * (Z_B/n_b)²
          'abs'  : A[a,b] = scale * |S(a,b)|
          'kappa': A[a,b] = scale * S²(a,b) * |KINETIC_SCALE|

        Only s-orbital (l=0) pairs are computed (overlap for l>0 not
        yet implemented).

        Parameters
        ----------
        rows, cols, data : list
            Adjacency COO arrays, extended in-place.
        nA : int
            Number of spatial orbitals on atom A.
        N : int
            Total spatial orbitals.

        Returns
        -------
        int
            Number of overlap edges added.
        """
        states_A = self._li_A.lattice.states
        states_B = self._li_B.lattice.states
        Z_orb_A = self._Z_orb_A
        Z_orb_B = self._Z_orb_B
        R = self.R
        scale = self.overlap_edge_scale
        mode = self.overlap_edges
        Z_A = float(self.Z_A)
        Z_B = float(self.Z_B)

        # Cache overlap by (n_a, n_b)
        s_cache: Dict[Tuple[int, int], float] = {}

        def get_overlap(n_a: int, n_b: int) -> float:
            key = (n_a, n_b)
            if key not in s_cache:
                s_cache[key] = compute_overlap_element(
                    n_a, 0, n_b, 0, Z_orb_A, Z_orb_B, R
                )
            return s_cache[key]

        n_edges = 0
        for i_a, (n_a, l_a, m_a) in enumerate(states_A):
            if l_a > 0:
                continue  # only s-orbital overlaps
            for i_b, (n_b, l_b, m_b) in enumerate(states_B):
                if l_b > 0:
                    continue
                s_ab = get_overlap(n_a, n_b)
                if abs(s_ab) < 1e-12:
                    continue

                # Compute edge weight based on mode
                if mode == 's2':
                    w = scale * s_ab * s_ab
                elif mode == 'fock':
                    w = scale * s_ab * s_ab * (Z_A / n_a)**2 * (Z_B / n_b)**2
                elif mode == 'abs':
                    w = scale * abs(s_ab)
                elif mode == 'kappa':
                    w = scale * s_ab * s_ab * abs(self.kinetic_scale)
                else:
                    raise ValueError(
                        f"Unknown overlap_edges mode '{mode}'. "
                        f"Valid: 's2', 'fock', 'abs', 'kappa'"
                    )

                if abs(w) < 1e-12:
                    continue

                j_b = nA + i_b
                # Symmetric edges
                rows.extend([i_a, j_b])
                cols.extend([j_b, i_a])
                data.extend([w, w])
                n_edges += 1

        return n_edges

    def _build_molecular_h1(self) -> None:
        """
        Build one-electron Hamiltonian for the molecule.

        Diagonal: exact atomic eigenvalues -Z²/(2n²) + cross-nuclear attraction
        (electrostatic potential, s-orbitals only for variational balance).
        Off-diagonal: kinetic hopping from combined adjacency (with bridges),
        or D-matrix cross-atom coupling when use_dmatrix=True.
        """
        n = self._n_spatial
        nA = self._n_spatial_A

        # --- Atom-dependent Sturmian path (v0.9.22): per-atom p0 ---
        if self.use_sturmian == 'atomic' and not self._ghost_A and not self._ghost_B:
            self._build_atomic_sturmian_h1()
            return

        # --- MO-projected Sturmian path (v0.9.30): prolate spheroidal betas ---
        if self.use_sturmian == 'molecular' and not self._ghost_A and not self._ghost_B:
            self._build_molecular_sturmian_h1(self._sturmian_p0)
            return

        # --- Sturmian path (Paper 9): Sturmian diagonal + exact D-matrix ---
        if self.use_sturmian and not self._ghost_A and not self._ghost_B:
            self._build_sturmian_h1(self._sturmian_p0)
            return

        # Diagonal: atomic eigenvalues with orbital exponent scaling.
        # For orbital exponent zeta_scale * Z / n (Z_eff = zeta * Z):
        #   E_n = Z_eff^2/(2n^2) - Z*Z_eff/n^2 = Z_eff*(Z_eff - 2Z)/(2n^2)
        # At zeta=1: E_n = Z*(Z-2Z)/(2n^2) = -Z^2/(2n^2) (standard).
        h1_diag = np.zeros(n)
        if not self._ghost_A:
            Z_eff_A = self._Z_orb_A
            Z_A_f = float(self.Z_A)
            for i, (ni, li, mi) in enumerate(self._li_A.lattice.states):
                h1_diag[i] = Z_eff_A * (Z_eff_A - 2.0 * Z_A_f) / (2.0 * ni**2)
        if not self._ghost_B:
            Z_eff_B = self._Z_orb_B
            Z_B_f = float(self.Z_B)
            for j, (nj, lj, mj) in enumerate(self._li_B.lattice.states):
                h1_diag[nA + j] = Z_eff_B * (Z_eff_B - 2.0 * Z_B_f) / (2.0 * nj**2)

        if self.use_dmatrix == 'hybrid' and not self._ghost_A and not self._ghost_B:
            # --- Hybrid path (v0.9.18): cross-nuclear diagonal + SW off-diagonal ---
            self._apply_cross_nuclear_diagonal(h1_diag)
            self._h1_diag = h1_diag

            # Off-diagonal: intra-atom graph Laplacian only (no inter-atomic bridges)
            from scipy.sparse import block_diag as sp_block_diag
            A_intra = sp_block_diag(
                [self._li_A.lattice.adjacency, self._li_B.lattice.adjacency],
                format='csr'
            )
            H1_offdiag_intra = self.kinetic_scale * (-A_intra)

            # Cross-atom off-diagonal: SW D-matrix elements (same-n)
            H1_dmatrix = self._cross_atom_h1_dmatrix(self.R)

            # Cross-n bridges: Wolfsberg-Helmholz overlap coupling (cross-n)
            if self.use_cross_n_bridges:
                H1_cross_n = self._cross_n_overlap_bridges(self.R)
            else:
                H1_cross_n = csr_matrix((n, n))

            self._H1_spatial = (
                diags(h1_diag) + H1_offdiag_intra + H1_dmatrix + H1_cross_n
            ).tocsr()

            # Build offdiag dict for matrix assembly path
            H1_offdiag_full = (
                H1_offdiag_intra + H1_dmatrix + H1_cross_n
            ).tocoo()
            self._h1_offdiag: Dict[int, List[Tuple[int, float]]] = {
                i: [] for i in range(n)
            }
            for r, c, v in zip(H1_offdiag_full.row, H1_offdiag_full.col,
                               H1_offdiag_full.data):
                if r != c and abs(v) >= self.threshold:
                    self._h1_offdiag[r].append((c, float(v)))

            nnz_cross = H1_dmatrix.nnz
            nnz_cross_n = H1_cross_n.nnz
            print(f"[MolecularLatticeIndex] H1 (hybrid): {n} spatial states, "
                  f"intra nnz={H1_offdiag_intra.nnz}, "
                  f"cross nnz={nnz_cross}, cross-n nnz={nnz_cross_n}")

        elif self.use_dmatrix and self.use_dmatrix != 'hybrid' and not self._ghost_A and not self._ghost_B:
            # --- D-matrix path (Paper 8): cross-atom H1 via SO(4) rotation ---
            # No Fourier cross-nuclear attraction on diagonal.
            # All cross-atom coupling enters through D-matrix off-diagonal.
            self._h1_diag = h1_diag

            # Intra-atom off-diagonal: kinetic hopping within each atom only
            # (block-diagonal part of adjacency, no inter-atomic bridges)
            from scipy.sparse import block_diag as sp_block_diag
            A_intra = sp_block_diag(
                [self._li_A.lattice.adjacency, self._li_B.lattice.adjacency],
                format='csr'
            )
            H1_offdiag_intra = self.kinetic_scale * (-A_intra)

            # Cross-atom off-diagonal: D-matrix elements (same-n)
            H1_dmatrix = self._cross_atom_h1_dmatrix(self.R)

            # Cross-n bridges: Wolfsberg-Helmholz overlap coupling (cross-n)
            if self.use_cross_n_bridges:
                H1_cross_n = self._cross_n_overlap_bridges(self.R)
            else:
                H1_cross_n = csr_matrix((n, n))

            self._H1_spatial = (
                diags(h1_diag) + H1_offdiag_intra + H1_dmatrix + H1_cross_n
            ).tocsr()

            # Build offdiag dict for matrix assembly path
            H1_offdiag_full = (
                H1_offdiag_intra + H1_dmatrix + H1_cross_n
            ).tocoo()
            self._h1_offdiag: Dict[int, List[Tuple[int, float]]] = {
                i: [] for i in range(n)
            }
            for r, c, v in zip(H1_offdiag_full.row, H1_offdiag_full.col,
                               H1_offdiag_full.data):
                if r != c and abs(v) >= self.threshold:
                    self._h1_offdiag[r].append((c, float(v)))

            nnz_cross = H1_dmatrix.nnz
            nnz_cross_n = H1_cross_n.nnz
            print(f"[MolecularLatticeIndex] H1 (D-matrix): {n} spatial states, "
                  f"intra nnz={H1_offdiag_intra.nnz}, "
                  f"cross nnz={nnz_cross}, cross-n nnz={nnz_cross_n}")
        else:
            # --- Standard path: cross-nuclear + bridge hopping ---
            self._apply_cross_nuclear_diagonal(h1_diag)
            self._apply_overlap_kinetic_correction(h1_diag)
            self._h1_diag = h1_diag

            # Off-diagonal: kinetic hopping from combined graph Laplacian
            H1_offdiag = self.kinetic_scale * (-self._adjacency_combined)
            self._H1_spatial = (diags(h1_diag) + H1_offdiag).tocsr()

            # Build offdiag dict for matrix assembly path
            H1_coo = H1_offdiag.tocoo()
            self._h1_offdiag: Dict[int, List[Tuple[int, float]]] = {
                i: [] for i in range(n)
            }
            for r, c, v in zip(H1_coo.row, H1_coo.col, H1_coo.data):
                if r != c and abs(v) >= self.threshold:
                    self._h1_offdiag[r].append((c, float(v)))

            # Report H1 statistics
            n_cross_A = sum(1 for i in range(nA) if h1_diag[i] !=
                            -float(self.Z_A)**2 / (2.0 * self._li_A.lattice.states[i][0]**2))
            print(f"[MolecularLatticeIndex] H1: {n} spatial states, "
                  f"off-diag nnz={H1_offdiag.nnz}")

    def _compute_overlap_attenuation(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute per-orbital overlap-based attenuation factors for cross-nuclear.

        For orbital a on atom A: S_a = sum_b S(a,b)^2 over all s-orbitals b on B.
        Returns (f_A, f_B) arrays of attenuation factors, one per spatial orbital.

        If cross_nuclear_beta is set, alpha becomes n-dependent:
            alpha_eff(a) = alpha_0 * (n_a / Z_a)^beta
        Core orbitals (low n/Z) are attenuated less; valence (high n/Z) more.
        """
        nA = self._n_spatial_A
        nB = self._n_spatial_B
        states_A = self._li_A.lattice.states
        states_B = self._li_B.lattice.states
        Z_orb_A = self._Z_orb_A
        Z_orb_B = self._Z_orb_B
        R = self.R
        alpha0 = self.cross_nuclear_alpha
        beta = self.cross_nuclear_beta
        mode = self.cross_nuclear_attenuation

        # Compute per-orbital effective alpha
        if beta is not None:
            # n-dependent: alpha_eff = alpha0 * (n/Z)^beta
            alpha_A = np.array([
                alpha0 * (float(n) / float(self.Z_A)) ** beta
                for (n, l, m) in states_A
            ])
            alpha_B = np.array([
                alpha0 * (float(n) / float(self.Z_B)) ** beta
                for (n, l, m) in states_B
            ])
        else:
            alpha_A = np.full(nA, alpha0)
            alpha_B = np.full(nB, alpha0)

        # Compute total squared overlap for each orbital
        S2_A = np.zeros(nA)
        S2_B = np.zeros(nB)

        # Only s-orbital (l=0) overlaps are implemented; l>0 contribute
        # negligibly to cross-nuclear attraction for ground states.
        s_cache: Dict[Tuple[int, int], float] = {}
        for i_a, (n_a, l_a, m_a) in enumerate(states_A):
            if l_a != 0:
                continue
            for i_b, (n_b, l_b, m_b) in enumerate(states_B):
                if l_b != 0:
                    continue
                key = (n_a, n_b)
                if key not in s_cache:
                    s_cache[key] = compute_overlap_element(
                        n_a, 0, n_b, 0, Z_orb_A, Z_orb_B, R
                    )
                s2 = s_cache[key] ** 2
                S2_A[i_a] += s2
                S2_B[i_b] += s2

        # Apply attenuation function with per-orbital alpha
        def attenuate(S: np.ndarray, alpha: np.ndarray) -> np.ndarray:
            if mode == 'linear':
                return np.maximum(1.0 - alpha * S, 0.0)
            elif mode == 'exp':
                return np.exp(-alpha * S)
            elif mode == 'pade':
                return 1.0 / (1.0 + alpha * S)
            elif mode == 'quadratic':
                return np.maximum(1.0 - S, 0.0) ** 2
            else:
                return np.ones_like(S)

        f_A = attenuate(S2_A, alpha_A)
        f_B = attenuate(S2_B, alpha_B)
        return f_A, f_B

    def _apply_cross_nuclear_diagonal(self, h1_diag: np.ndarray) -> None:
        """
        Apply cross-nuclear attraction to the h1 diagonal in-place.

        When cross_nuclear_method='exact', uses full 2D quadrature via
        compute_exact_cross_nuclear() for ALL (n,l,m) orbitals.
        When cross_nuclear_method='fourier', uses the shell-theorem formula
        for s-orbitals only (v0.9.35 baseline).

        If cross_nuclear_attenuation is set, each orbital's cross-nuclear
        attraction is multiplied by f(S_a) where S_a is the total squared
        overlap with the partner atom's orbitals.

        Orbital shapes use Z_orb (= zeta * Z_nuclear) to account for
        exponent relaxation, while the attracting charge remains Z_other.

        Skips ghost atoms (Z=0).
        """
        if self._ghost_A or self._ghost_B:
            return

        nA = self._n_spatial_A
        # Orbital shape charges (zeta-scaled); attracting charges are bare Z
        Z_orb_A = self._Z_orb_A
        Z_orb_B = self._Z_orb_B

        # Compute attenuation factors if requested
        if self.cross_nuclear_attenuation:
            f_A, f_B = self._compute_overlap_attenuation()
        else:
            f_A = np.ones(nA)
            f_B = np.ones(self._n_spatial_B)

        if self.cross_nuclear_method == 'exact':
            # Exact 2D quadrature for all orbitals
            for i, (ni, li, mi) in enumerate(self._li_A.lattice.states):
                h1_diag[i] += f_A[i] * compute_exact_cross_nuclear(
                    ni, li, mi, Z_orb_A, float(self.Z_B), self.R)
            for j, (nj, lj, mj) in enumerate(self._li_B.lattice.states):
                h1_diag[nA + j] += f_B[j] * compute_exact_cross_nuclear(
                    nj, lj, mj, Z_orb_B, float(self.Z_A), self.R)
        else:
            # Fourier/shell-theorem: s-orbitals only (l=0)
            for i, (ni, li, mi) in enumerate(self._li_A.lattice.states):
                if li == 0:
                    h1_diag[i] += f_A[i] * self._fourier_cross_attraction(
                        ni, li, Z_orb_A, float(self.Z_B), self.R)
            for j, (nj, lj, mj) in enumerate(self._li_B.lattice.states):
                if lj == 0:
                    h1_diag[nA + j] += f_B[j] * self._fourier_cross_attraction(
                        nj, lj, Z_orb_B, float(self.Z_A), self.R)

    def _apply_overlap_kinetic_correction(self, h1_diag: np.ndarray) -> None:
        """
        Add overlap-dependent kinetic correction to the H1 diagonal.

        For each orbital a on atom A, adds:
            h1_diag[a] += t_corr_lambda * sum_b S(a,b)^2 * w(a,b)

        where S(a,b) is the STO overlap between orbital a and orbital b
        on the other atom, and w(a,b) is:

        - If t_corr_fock_weighted: w = (Z_A/n_a)^2 * (Z_B/n_b)^2
          This weights by the Fock energy-shell p0^2 = Z^2/n^2, suppressing
          diffuse orbitals and concentrating repulsion on compact cores.

        - Otherwise: w = 1 (uniform weighting).

        Only s-orbital (l=0) cross-atom overlaps are computed.

        Parameters
        ----------
        h1_diag : np.ndarray
            Diagonal of H1 matrix, modified in-place.
        """
        if abs(self.t_corr_lambda) < 1e-15:
            return
        if self._ghost_A or self._ghost_B:
            return

        nA = self._n_spatial_A
        states_A = self._li_A.lattice.states
        states_B = self._li_B.lattice.states
        Z_orb_A = self._Z_orb_A
        Z_orb_B = self._Z_orb_B
        R = self.R
        fock = self.t_corr_fock_weighted
        Z_A = float(self.Z_A)
        Z_B = float(self.Z_B)

        # Cache overlap by (n_a, n_b) since only l=0 pairs are computed
        s_cache: Dict[Tuple[int, int], float] = {}

        def get_overlap(n_a: int, n_b: int) -> float:
            key = (n_a, n_b)
            if key not in s_cache:
                s_cache[key] = compute_overlap_element(
                    n_a, 0, n_b, 0, Z_orb_A, Z_orb_B, R
                )
            return s_cache[key]

        # Atom A orbitals: sum over B partners
        for i_a, (n_a, l_a, m_a) in enumerate(states_A):
            if l_a > 0:
                continue
            sum_s2w = 0.0
            for (n_b, l_b, m_b) in states_B:
                if l_b > 0:
                    continue
                s_ab = get_overlap(n_a, n_b)
                w = (Z_A / n_a)**2 * (Z_B / n_b)**2 if fock else 1.0
                sum_s2w += s_ab * s_ab * w
            h1_diag[i_a] += self.t_corr_lambda * sum_s2w

        # Atom B orbitals: sum over A partners
        for i_b, (n_b, l_b, m_b) in enumerate(states_B):
            if l_b > 0:
                continue
            sum_s2w = 0.0
            for (n_a, l_a, m_a) in states_A:
                if l_a > 0:
                    continue
                s_ab = get_overlap(n_a, n_b)
                w = (Z_A / n_a)**2 * (Z_B / n_b)**2 if fock else 1.0
                sum_s2w += s_ab * s_ab * w
            h1_diag[nA + i_b] += self.t_corr_lambda * sum_s2w

    def _build_atomic_sturmian_cross_nuclear(
        self, R: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Per-atom cross-nuclear attraction with atom-dependent p0 (v0.9.22).

        Each atom's orbitals feel the other nucleus through D-matrix elements
        evaluated at that atom's own p0 and corresponding bond angle gamma:

            A-A block: -(Z_B / p0_A) * D^(n)(gamma_A)
            B-B block: -(Z_A / p0_B) * D^(n)(gamma_B)

        where gamma_alpha = arccos((p0_alpha^2 - p_R^2)/(p0_alpha^2 + p_R^2))
        uses the atom-specific p0.  This breaks the single-S3 geometry but
        restores binding for heteronuclear molecules.

        The H cross-nuclear (1s feeling Li nucleus) can be very strong:
        -(Z_A/p0_B) = -3 Ha when Z_A=3, p0_B=1.  If this exceeds -Z_A/R
        (the point-charge limit), it is capped at -Z_A/R.

        Parameters
        ----------
        R : float
            Internuclear distance in Bohr.

        Returns
        -------
        V_cross_A : np.ndarray
            Shape (n_spatial_A, n_spatial_A). Nucleus B acting on A orbitals.
        V_cross_B : np.ndarray
            Shape (n_spatial_B, n_spatial_B). Nucleus A acting on B orbitals.
        """
        from .wigner_so4 import bond_angle, d_matrix_block

        nA = self._n_spatial_A
        nB = self._n_spatial_B
        Z_A = float(self.Z_A)
        Z_B = float(self.Z_B)
        p0_A = self._sturmian_p0_A
        p0_B = self._sturmian_p0_B

        gamma_A = bond_angle(R, p0_A)
        gamma_B = bond_angle(R, p0_B)
        self._gamma_A = gamma_A
        self._gamma_B = gamma_B

        print(f"[Atomic Sturmian] p0_A={p0_A:.4f}, p0_B={p0_B:.4f}, "
              f"gamma_A={gamma_A:.4f} rad ({np.degrees(gamma_A):.1f} deg), "
              f"gamma_B={gamma_B:.4f} rad ({np.degrees(gamma_B):.1f} deg)")

        # Cap: H cross-nuclear -(Z_A/p0_B) must not exceed -Z_A/R
        cap_B = Z_A / R  # maximum magnitude for B-block diagonal

        # V_cross_A: nucleus B on A orbitals = -(Z_B/p0_A) * D_sym(gamma_A)
        V_cross_A = np.zeros((nA, nA))
        offset = 0
        for n_shell in range(1, self.nmax_A + 1):
            n_sq = n_shell * n_shell
            D_raw = d_matrix_block(n_shell, gamma_A)
            D_sym = (D_raw + D_raw.T) / 2.0
            scale = -(Z_B / p0_A)
            V_cross_A[offset:offset + n_sq, offset:offset + n_sq] = (
                scale * D_sym
            )
            offset += n_sq

        # V_cross_B: nucleus A on B orbitals = -(Z_A/p0_B) * D_sym(gamma_B)
        # with capping: diagonal elements capped at -Z_A/R
        V_cross_B = np.zeros((nB, nB))
        V_cross_B_uncapped = np.zeros((nB, nB))
        offset = 0
        n_capped = 0
        for n_shell in range(1, self.nmax_B + 1):
            n_sq = n_shell * n_shell
            D_raw = d_matrix_block(n_shell, gamma_B)
            D_sym = (D_raw + D_raw.T) / 2.0
            scale = -(Z_A / p0_B)
            block_uncapped = scale * D_sym
            V_cross_B_uncapped[offset:offset + n_sq, offset:offset + n_sq] = (
                block_uncapped
            )
            # Cap diagonal elements
            block_capped = block_uncapped.copy()
            for k in range(n_sq):
                if block_capped[k, k] < -cap_B:
                    block_capped[k, k] = -cap_B
                    n_capped += 1
            V_cross_B[offset:offset + n_sq, offset:offset + n_sq] = (
                block_capped
            )
            offset += n_sq

        self._V_cross_B_uncapped = V_cross_B_uncapped
        if n_capped > 0:
            print(f"[Atomic Sturmian] WARNING: {n_capped} B-block diagonal "
                  f"elements capped at -Z_A/R = -{cap_B:.4f} Ha "
                  f"(uncapped: -{Z_A/p0_B:.4f} Ha)")

        return V_cross_A, V_cross_B

    def _build_atomic_sturmian_h1(self) -> None:
        """
        Build one-electron Hamiltonian with atom-dependent Sturmian p0 (v0.9.22).

        Each atom uses its own self-consistent momentum scale p0_alpha
        from the isolated-atom FCI.  The Sturmian diagonal for atom alpha's
        orbital (n,l,m) is:

            epsilon_n(p0_alpha) = p0_alpha^2/2 - Z_alpha * p0_alpha / n

        Cross-nuclear uses per-atom D-matrix with atom-specific gamma.
        Off-diagonal A-B hopping uses geometric-mean p0_AB = sqrt(p0_A * p0_B).

        This breaks the single-S3 geometry but restores binding for
        heteronuclear molecules with large Z_A/Z_B ratios.
        """
        n = self._n_spatial
        nA = self._n_spatial_A
        p0_A = self._sturmian_p0_A
        p0_B = self._sturmian_p0_B
        p0_AB = self._sturmian_p0_AB

        print(f"\n[Atomic Sturmian H1] p0_A={p0_A:.4f} (Z_A={self.Z_A}), "
              f"p0_B={p0_B:.4f} (Z_B={self.Z_B}), "
              f"p0_AB={p0_AB:.4f}")

        # --- Backward compatibility assertion ---
        for Z_val, p0_val in [(self.Z_A, p0_A), (self.Z_B, p0_B)]:
            if Z_val == 0:
                continue
            Z = float(Z_val)
            # At p0 = Z/n, the Sturmian diagonal equals -Z^2/(2n^2)
            for n_shell in range(1, max(self.nmax_A, self.nmax_B) + 1):
                p0_test = Z / n_shell
                eps_s = p0_test**2 / 2.0 - Z * p0_test / n_shell
                eps_a = -Z**2 / (2.0 * n_shell**2)
                assert abs(eps_s - eps_a) < 1e-10, (
                    f"Backward compat FAILED: Z={Z_val}, n={n_shell}"
                )

        # --- Per-atom Sturmian diagonal ---
        h1_diag = np.zeros(n)
        for i, (ni, li, mi) in enumerate(self._li_A.lattice.states):
            h1_diag[i] = p0_A**2 / 2.0 - float(self.Z_A) * p0_A / ni
        for j, (nj, lj, mj) in enumerate(self._li_B.lattice.states):
            h1_diag[nA + j] = p0_B**2 / 2.0 - float(self.Z_B) * p0_B / nj

        # --- Per-atom cross-nuclear via D-matrix ---
        V_cross_A, V_cross_B = self._build_atomic_sturmian_cross_nuclear(
            self.R
        )

        # Store diagonal for diagnostics
        self._h1_diag = h1_diag.copy()

        # --- Off-diagonal: intra-atom graph Laplacian (no bridges) ---
        from scipy.sparse import block_diag as sp_block_diag
        A_intra = sp_block_diag(
            [self._li_A.lattice.adjacency, self._li_B.lattice.adjacency],
            format='csr'
        )
        H1_offdiag_intra = self.kinetic_scale * (-A_intra)

        # --- Cross-atom off-diagonal: SW D-matrix with p0_AB ---
        H1_dmatrix = self._cross_atom_h1_dmatrix(self.R, p0=p0_AB)

        # --- Assemble ---
        from scipy.sparse import block_diag as sp_block_diag2
        V_cross_sparse = sp_block_diag2(
            [csr_matrix(V_cross_A), csr_matrix(V_cross_B)],
            format='csr'
        )

        self._H1_spatial = (
            diags(h1_diag) + V_cross_sparse + H1_offdiag_intra + H1_dmatrix
        ).tocsr()

        # Update h1_diag to include cross-nuclear diagonal for FCI
        for i in range(nA):
            self._h1_diag[i] += V_cross_A[i, i]
        for j in range(self._n_spatial_B):
            self._h1_diag[nA + j] += V_cross_B[j, j]

        # Build offdiag dict for matrix/direct assembly
        H1_offdiag_full = (
            V_cross_sparse + H1_offdiag_intra + H1_dmatrix
        ).tocoo()
        self._h1_offdiag: Dict[int, List[Tuple[int, float]]] = {
            i: [] for i in range(n)
        }
        for r, c, v in zip(H1_offdiag_full.row, H1_offdiag_full.col,
                           H1_offdiag_full.data):
            if r != c and abs(v) >= self.threshold:
                self._h1_offdiag[r].append((c, float(v)))

        # Diagnostics
        print(f"[Atomic Sturmian H1] Li 1s diag: {h1_diag[0]:.6f} Ha "
              f"(cross-nuc: {V_cross_A[0,0]:.6f} Ha)")
        print(f"[Atomic Sturmian H1] H 1s diag: {h1_diag[nA]:.6f} Ha "
              f"(cross-nuc: {V_cross_B[0,0]:.6f} Ha)")
        print(f"[Atomic Sturmian H1] {n} spatial, "
              f"cross_AA={V_cross_sparse.nnz}, "
              f"cross_AB={H1_dmatrix.nnz}, intra={H1_offdiag_intra.nnz}")

    def _build_sturmian_cross_nuclear(
        self, p0: float, R: float
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Sturmian cross-nuclear attraction via D-matrix (Paper 9, Eq. 23).

        In the Sturmian basis, the cross-nuclear matrix element for
        orbitals on center A in the field of nucleus B is:

            ⟨S^A_{n'l'm'} | -Z_B/r_B | S^A_{nlm}⟩
                = -(Z_B/p0) · D^(n)_{(l'm'),(lm)}(γ)

        Block-diagonal in n (same-shell only).  The D-matrix blocks are
        symmetrized (D+D^T)/2 for Hermiticity.

        This replaces the frozen Fourier cross-nuclear of v0.9.20.  The
        result is fully p0-dependent through both the 1/p0 prefactor and
        γ(p0, R), resolving the Fourier–Sturmian inconsistency.

        Parameters
        ----------
        p0 : float
            Shared Sturmian momentum scale.
        R : float
            Internuclear distance in Bohr.

        Returns
        -------
        V_cross_A : np.ndarray
            Shape (n_spatial_A, n_spatial_A).  Nucleus B acting on A orbitals.
        V_cross_B : np.ndarray
            Shape (n_spatial_B, n_spatial_B).  Nucleus A acting on B orbitals.
        """
        from .wigner_so4 import bond_angle, d_matrix_block

        nA = self._n_spatial_A
        nB = self._n_spatial_B
        Z_A = float(self.Z_A)
        Z_B = float(self.Z_B)

        gamma = bond_angle(R, p0)

        # Pre-compute symmetrized D-matrix blocks per shell
        n_max = max(self.nmax_A, self.nmax_B)
        D_blocks: Dict[int, np.ndarray] = {}
        for n_shell in range(1, n_max + 1):
            D_raw = d_matrix_block(n_shell, gamma)
            D_blocks[n_shell] = (D_raw + D_raw.T) / 2.0

        # V_cross_A: nucleus B on A orbitals = -(Z_B/p0) * D_sym
        V_cross_A = np.zeros((nA, nA))
        offset = 0
        for n_shell in range(1, self.nmax_A + 1):
            n_sq = n_shell * n_shell
            scale = -(Z_B / p0)
            V_cross_A[offset:offset + n_sq, offset:offset + n_sq] = (
                scale * D_blocks[n_shell]
            )
            offset += n_sq

        # V_cross_B: nucleus A on B orbitals = -(Z_A/p0) * D_sym
        V_cross_B = np.zeros((nB, nB))
        offset = 0
        for n_shell in range(1, self.nmax_B + 1):
            n_sq = n_shell * n_shell
            scale = -(Z_A / p0)
            V_cross_B[offset:offset + n_sq, offset:offset + n_sq] = (
                scale * D_blocks[n_shell]
            )
            offset += n_sq

        return V_cross_A, V_cross_B

    def _build_molecular_sturmian_h1(self, p0: float) -> None:
        """Build H1 with MO-projected Sturmian diagonal (v0.9.30).

        Uses prolate spheroidal molecular orbital betas projected onto
        atom-centered hydrogenic orbitals via dominant-overlap assignment.

        Diagonal: eps^alpha_{nlm} = p0^2/2 - beta^alpha_{nlm} * Z_alpha * p0
        Off-diagonal: intra-atom graph Laplacian + SW D-matrix cross-atom.

        Parameters
        ----------
        p0 : float
            Shared Sturmian momentum scale.
        """
        from .molecular_sturmian import (
            compute_molecular_sturmian_betas,
            project_mo_betas_to_atom_centers,
        )

        n = self._n_spatial
        nA = self._n_spatial_A
        self._sturmian_p0 = p0

        Z_A = float(self.Z_A)
        Z_B = float(self.Z_B)
        nmax = max(self.nmax_A, self.nmax_B)

        # --- Compute MO betas and project ---
        mo_results = compute_molecular_sturmian_betas(
            Z_A=Z_A, Z_B=Z_B, R=self.R, p0=p0, nmax=nmax
        )
        beta_A, beta_B = project_mo_betas_to_atom_centers(
            mo_results, Z_A=Z_A, Z_B=Z_B, R=self.R, p0=p0, nmax=nmax
        )

        # --- MO-projected Sturmian diagonal ---
        h1_diag = np.zeros(n)
        for i, (ni, li, mi) in enumerate(self._li_A.lattice.states):
            key = (ni, li, mi)
            beta = beta_A.get(key, 1.0 / ni)
            h1_diag[i] = p0**2 / 2.0 - beta * Z_A * p0

        for j, (nj, lj, mj) in enumerate(self._li_B.lattice.states):
            key = (nj, lj, mj)
            beta = beta_B.get(key, 1.0 / nj)
            h1_diag[nA + j] = p0**2 / 2.0 - beta * Z_B * p0

        self._h1_diag = h1_diag.copy()

        # --- Off-diagonal: intra-atom graph Laplacian ---
        from scipy.sparse import block_diag as sp_block_diag
        A_intra = sp_block_diag(
            [self._li_A.lattice.adjacency, self._li_B.lattice.adjacency],
            format='csr'
        )
        H1_offdiag_intra = self.kinetic_scale * (-A_intra)

        # --- Cross-atom off-diagonal: SW D-matrix ---
        H1_dmatrix = self._cross_atom_h1_dmatrix(self.R, p0=p0)

        # --- Assemble ---
        self._H1_spatial = (
            diags(h1_diag) + H1_offdiag_intra + H1_dmatrix
        ).tocsr()

        # Build offdiag dict for matrix/direct assembly
        H1_offdiag_full = (H1_offdiag_intra + H1_dmatrix).tocoo()
        self._h1_offdiag: Dict[int, List[Tuple[int, float]]] = {
            i: [] for i in range(n)
        }
        for r, c, v in zip(H1_offdiag_full.row, H1_offdiag_full.col,
                           H1_offdiag_full.data):
            if r != c and abs(v) >= self.threshold:
                self._h1_offdiag[r].append((c, float(v)))

        nnz_cross_AB = H1_dmatrix.nnz
        print(f"[MolecularLatticeIndex] H1 (MO-projected Sturmian): "
              f"{n} spatial, p0={p0:.4f}, "
              f"cross_AB={nnz_cross_AB}, intra={H1_offdiag_intra.nnz}")

    def _build_sturmian_h1(self, p0: float) -> None:
        """
        Build one-electron Hamiltonian with Sturmian diagonal (Paper 9).

        The Sturmian diagonal replaces the atomic eigenvalue -Z²/(2n²) with:

            ε_n(p0) = p0²/2 - Z·p0/n

        where p0 is the shared molecular momentum scale.  This is the
        expectation value ⟨S_nlm| -½∇² - Z/r |S_nlm⟩ evaluated via the
        Sturmian eigenvalue equation and virial theorem (Paper 9, Sec. VII).

        Note: Paper 9 Eq. (22) states ε_n = -p0²/2 - Zp0/n², which has a
        sign/exponent typo.  The correct formula ε_n = +p0²/2 - Zp0/n is
        derived from the Sturmian equation (Sec. IV) and passes the backward
        compatibility check ε_n(Z/n) = -Z²/(2n²) exactly.

        Cross-nuclear attraction uses the p0-dependent D-matrix formula
        (Paper 9, Eq. 23) in the A-A and B-B diagonal blocks.  The A-B
        off-diagonal block uses SW D-matrix hopping with sturmian=True.

        Parameters
        ----------
        p0 : float
            Shared Sturmian momentum scale.
        """
        n = self._n_spatial
        nA = self._n_spatial_A
        self._sturmian_p0 = p0

        # --- Backward compatibility assertion (Paper 9, Sec. VII) ---
        # At p0 = Z/n, the Sturmian diagonal must equal the atomic
        # eigenvalue -Z²/(2n²).
        for Z_val in [self.Z_A, self.Z_B]:
            if Z_val == 0:
                continue
            Z = float(Z_val)
            for n_shell in range(1, max(self.nmax_A, self.nmax_B) + 1):
                p0_atomic = Z / n_shell
                eps_sturmian = p0_atomic**2 / 2.0 - Z * p0_atomic / n_shell
                eps_atomic = -Z**2 / (2.0 * n_shell**2)
                assert abs(eps_sturmian - eps_atomic) < 1e-10, (
                    f"Backward compatibility FAILED: Z={Z_val}, n={n_shell}, "
                    f"ε_Sturmian(Z/n)={eps_sturmian:.12f} vs "
                    f"ε_atomic={eps_atomic:.12f}"
                )

        # --- Sturmian diagonal: ε_n(p0) = p0²/2 - Z·p0/n ---
        h1_diag = np.zeros(n)
        for i, (ni, li, mi) in enumerate(self._li_A.lattice.states):
            h1_diag[i] = p0**2 / 2.0 - float(self.Z_A) * p0 / ni
        for j, (nj, lj, mj) in enumerate(self._li_B.lattice.states):
            h1_diag[nA + j] = p0**2 / 2.0 - float(self.Z_B) * p0 / nj

        # --- Cross-nuclear via D-matrix (Paper 9 Eq. 23, v0.9.21) ---
        # Replaces frozen Fourier cross-nuclear of v0.9.20.
        V_cross_A, V_cross_B = self._build_sturmian_cross_nuclear(p0, self.R)

        # Store diagonal for diagnostics before adding cross-nuclear
        self._h1_diag = h1_diag.copy()

        # --- Off-diagonal: intra-atom graph Laplacian (no bridges) ---
        from scipy.sparse import block_diag as sp_block_diag
        A_intra = sp_block_diag(
            [self._li_A.lattice.adjacency, self._li_B.lattice.adjacency],
            format='csr'
        )
        H1_offdiag_intra = self.kinetic_scale * (-A_intra)

        # --- Cross-atom off-diagonal: SW D-matrix (sturmian=True) ---
        H1_dmatrix = self._cross_atom_h1_dmatrix(self.R, p0=p0)

        # --- Assemble: diagonal + cross-nuclear (A-A, B-B) + intra + AB ---
        # Build cross-nuclear as sparse block-diagonal matrix
        from scipy.sparse import block_diag as sp_block_diag2
        V_cross_sparse = sp_block_diag2(
            [csr_matrix(V_cross_A), csr_matrix(V_cross_B)],
            format='csr'
        )

        self._H1_spatial = (
            diags(h1_diag) + V_cross_sparse + H1_offdiag_intra + H1_dmatrix
        ).tocsr()

        # Update h1_diag to include cross-nuclear diagonal for FCI
        for i in range(nA):
            self._h1_diag[i] += V_cross_A[i, i]
        for j in range(self._n_spatial_B):
            self._h1_diag[nA + j] += V_cross_B[j, j]

        # Build offdiag dict for matrix/direct assembly
        H1_offdiag_full = (
            V_cross_sparse + H1_offdiag_intra + H1_dmatrix
        ).tocoo()
        self._h1_offdiag: Dict[int, List[Tuple[int, float]]] = {
            i: [] for i in range(n)
        }
        for r, c, v in zip(H1_offdiag_full.row, H1_offdiag_full.col,
                           H1_offdiag_full.data):
            if r != c and abs(v) >= self.threshold:
                self._h1_offdiag[r].append((c, float(v)))

        nnz_cross_AB = H1_dmatrix.nnz
        nnz_cross_AA = V_cross_sparse.nnz
        print(f"[MolecularLatticeIndex] H1 (Sturmian v2): {n} spatial, "
              f"p0={p0:.4f}, cross_AA={nnz_cross_AA}, "
              f"cross_AB={nnz_cross_AB}, intra={H1_offdiag_intra.nnz}")

    def solve_sturmian_p0(
        self,
        R: float,
        p0_init: Optional[float] = None,
        tol: float = 1e-6,
        max_iter: int = 50,
        damping: float = 1.0,
    ) -> Tuple[float, float, int, bool]:
        """
        Self-consistent determination of p0 via p0 = sqrt(-2 * E_mol(p0)).

        Iteration (Paper 9, Sec. VI.B):
          1. Build H(p0), run FCI, get E_mol
          2. p0_target = sqrt(-2 * E_mol)
          3. p0_new = (1-damping)*p0 + damping*p0_target
          4. Repeat until |p0_new - p0| < tol

        Parameters
        ----------
        R : float
            Internuclear distance in Bohr.
        p0_init : float or None
            Initial p0.  If None, uses sqrt((Z_A² + Z_B²) / 2).
        tol : float
            Convergence tolerance on |p0_new - p0|.
        max_iter : int
            Maximum iterations.
        damping : float
            Mixing parameter in (0, 1].  1.0 = no damping (bare iteration),
            0.3 = 30% update per step.  Default 1.0.

        Returns
        -------
        p0_converged : float
            Final p0 value.
        E_mol_converged : float
            Final molecular energy (electronic + V_NN).
        n_iterations : int
            Number of iterations performed.
        converged : bool
            Whether the iteration converged to within tol.
        """
        if not self.use_sturmian:
            raise ValueError("solve_sturmian_p0 requires use_sturmian=True")
        if damping <= 0 or damping > 1.0:
            raise ValueError(f"damping must be in (0, 1], got {damping}")

        # Initialize p0 (Paper 9, Sec. VI.B, step 1)
        if p0_init is not None:
            p0 = p0_init
        else:
            p0 = np.sqrt((float(self.Z_A)**2 + float(self.Z_B)**2) / 2.0)

        damp_str = f", damping={damping:.2f}" if damping < 1.0 else ""
        print(f"\n{'='*60}")
        print(f"Sturmian self-consistency loop: R={R:.3f} bohr{damp_str}")
        print(f"p0_init = {p0:.6f}")
        print(f"{'='*60}")

        E_mol = 0.0
        converged = False

        # Dispatch to correct H1 builder
        def _rebuild_h1(p0_val: float) -> None:
            if self.use_sturmian == 'molecular':
                self._build_molecular_sturmian_h1(p0_val)
            else:
                self._build_sturmian_h1(p0_val)

        for iteration in range(1, max_iter + 1):
            # Rebuild H1 with current p0
            self.R = R
            _rebuild_h1(p0)

            # Run FCI
            eigvals, eigvecs = self.compute_ground_state(n_states=1)
            E_mol = eigvals[0]

            # Update p0 = sqrt(-2 * E_mol) with optional damping
            if E_mol >= 0:
                print(f"  iter {iteration}: p0={p0:.6f}, E_mol={E_mol:.6f} "
                      f"(POSITIVE — unbound, stopping)")
                break

            p0_target = np.sqrt(-2.0 * E_mol)
            p0_new = (1.0 - damping) * p0 + damping * p0_target
            delta = abs(p0_new - p0)

            print(f"  iter {iteration}: p0={p0:.6f}, E_mol={E_mol:.6f} Ha, "
                  f"p0_target={p0_target:.6f}, p0_new={p0_new:.6f}, "
                  f"delta={delta:.2e}")

            if delta < tol:
                p0 = p0_new
                converged = True
                print(f"  CONVERGED at iteration {iteration}: "
                      f"p0={p0:.6f}, E_mol={E_mol:.6f} Ha")
                break

            # Divergence check: if delta is growing after 5 iterations
            if iteration >= 5 and delta > 1.0:
                print(f"  DIVERGING at iteration {iteration}: "
                      f"delta={delta:.4f} > 1.0, stopping early")
                p0 = p0_new
                break

            p0 = p0_new

        # Rebuild with final p0
        self._sturmian_p0 = p0
        _rebuild_h1(p0)

        n_iterations = min(iteration, max_iter)
        return p0, E_mol, n_iterations, converged

    @staticmethod
    def _mulliken_cross_attraction(
        n: int, l: int, Z_self: int, Z_other: int, R: float
    ) -> float:
        """
        Cross-nuclear attraction via Mulliken approximation with variational cap.

        Replicates MoleculeHamiltonian._build_cross_nuclear_attraction logic:
            V_cross = max(-Z_other/R * S_eff, -Z_other * Z_self / n²)
        where the cap prevents variational collapse for diffuse orbitals.
        """
        if R < 1e-10:
            return 0.0
        R_eff = R * float(Z_self) / (n**2)
        S_1s = np.exp(-R_eff) * (1.0 + R_eff + R_eff**2 / 3.0)
        ang = 1.0 / (2 * l + 1)
        S_eff = min(S_1s * ang, 1.0)
        v_mulliken = (-float(Z_other) / R) * S_eff
        v_limit = -float(Z_other) * float(Z_self) / (n**2)
        return max(v_mulliken, v_limit)

    @staticmethod
    def _fourier_cross_attraction(
        n: int, l: int, Z_self: int, Z_other: int, R: float
    ) -> float:
        """
        Cross-nuclear attraction via direct electrostatic potential.

        For the spherically averaged orbital density, the potential of
        nucleus B at distance R from the electron's own nucleus is:

            V = -Z_other × [(1/R) ∫₀^R |R_{nl}(r)|² r² dr
                           + ∫_R^∞ |R_{nl}(r)|² r dr]

        Exact for the spherically averaged density. Correct limits:
        - V → -Z_other/R as R → ∞ (since ∫|R|²r²dr = 1)
        - Smooth, finite at all R > 0

        Uses direct r-space integration (fast, single quadrature).

        Parameters
        ----------
        n, l : int
            Quantum numbers of the orbital on atom A
        Z_self : int
            Nuclear charge of atom A (where the electron sits)
        Z_other : int
            Nuclear charge of atom B (the attracting nucleus)
        R : float
            Internuclear distance in Bohr

        Returns
        -------
        float
            Cross-nuclear attraction energy in Hartree (negative)
        """
        if R < 1e-10:
            return 0.0

        from scipy.integrate import quad
        from scipy.special import assoc_laguerre
        from math import factorial

        Z = float(Z_self)

        # Build |R_{nl,Z}(r)|² using closed-form for common cases
        if n == 1 and l == 0:
            def R_sq(r: float) -> float:
                return 4.0 * Z**3 * np.exp(-2.0 * Z * r)
        elif n == 2 and l == 0:
            def R_sq(r: float) -> float:
                x = Z * r
                return (Z**3 / 8.0) * (2.0 - x)**2 * np.exp(-x)
        elif n == 2 and l == 1:
            def R_sq(r: float) -> float:
                x = Z * r
                return (Z**3 / 24.0) * x**2 * np.exp(-x)
        else:
            nr = n - l - 1
            N_sq = (2.0 * Z / n)**3 * factorial(nr) / (2.0 * n * factorial(n + l))

            def R_sq(r: float) -> float:
                x = 2.0 * Z * r / n
                lag = assoc_laguerre(x, nr, 2 * l + 1)
                return N_sq * x**(2 * l) * lag * lag * np.exp(-x)

        # Electrostatic potential: Φ(R) = (1/R)∫₀^R ρr²dr + ∫_R^∞ ρr dr
        inner, _ = quad(lambda r: R_sq(r) * r * r, 0, R, limit=100)
        outer, _ = quad(lambda r: R_sq(r) * r, R, np.inf, limit=100)

        return -float(Z_other) * (inner / R + outer)

    def _cross_atom_h1_dmatrix(
        self, R: float, p0: Optional[float] = None
    ) -> csr_matrix:
        """
        Cross-atom one-electron coupling via Shibuya-Wulfman nuclear
        attraction integrals on S3.

        Uses the SW formula:
            V^SW = -(Z_B/p0) * D^(n)_{(l'm'),(lm)}(gamma) * f(n, gamma)

        where D is the SO(4) Wigner D-matrix (symmetrized: (D+D^T)/2)
        and f(n, gamma) = sin(gamma) is the form factor encoding
        R-dependence for fixed p0.

        The coupling for each (i_A, j_B) pair is placed symmetrically
        in both off-diagonal blocks to ensure Hermiticity. The effective
        nuclear charge uses the arithmetic mean (Z_A + Z_B)/2 for
        heteronuclear molecules.

        Parameters
        ----------
        R : float
            Internuclear distance in Bohr.
        p0 : float or None
            Momentum scale for bond angle. If None, uses
            sqrt(Z_A^2 + Z_B^2).

        Returns
        -------
        csr_matrix
            Shape (n_spatial, n_spatial). Cross-atom H1 block.
        """
        from .wigner_so4 import bond_angle, d_matrix_block
        from .shibuya_wulfman import sw_form_factor

        nA = self._n_spatial_A
        n = self._n_spatial
        Z_A = float(self.Z_A)
        Z_B = float(self.Z_B)

        # Determine p0 from atomic charges if not provided
        if p0 is None:
            p0 = np.sqrt(Z_A**2 + Z_B**2)
        self._dmatrix_p0 = p0

        gamma = bond_angle(R, p0)
        self._dmatrix_gamma = gamma

        # Effective nuclear charge: arithmetic mean for heteronuclear
        Z_eff = (Z_A + Z_B) / 2.0

        print(f"[SW] p0={p0:.4f}, gamma={gamma:.4f} rad "
              f"({np.degrees(gamma):.1f} deg), R={R:.3f} bohr, "
              f"Z_eff={Z_eff:.2f}")

        # --- OLD geometric-mean coupling (v0.9.16, FAILED) ---
        # kappa = abs(self.kinetic_scale)  # 1/16
        # scale = kappa * sqrt(Z_A * Z_B) / n_shell^2
        # h_val = D_elem * scale
        # Result: coupling ~0.108 Ha, 5x too weak, molecule unbound
        # ---

        states_A = self._li_A.lattice.states  # list of (n, l, m)
        states_B = self._li_B.lattice.states

        # Pre-compute symmetrized D-matrix blocks and form factors per shell
        n_max = max(s[0] for s in states_A)
        D_blocks: Dict[int, np.ndarray] = {}
        f_factors: Dict[int, float] = {}
        for n_shell in range(1, n_max + 1):
            D_raw = d_matrix_block(n_shell, gamma)
            D_blocks[n_shell] = (D_raw + D_raw.T) / 2.0  # Hermitian
            f_factors[n_shell] = sw_form_factor(n_shell, gamma,
                                               sturmian=self.use_sturmian)

        # Build per-shell state index maps for efficient lookup
        # shell_indices_A[n] = list of (local_idx_in_shell, global_spatial_idx)
        shell_indices_A: Dict[int, List[Tuple[int, int]]] = {}
        shell_indices_B: Dict[int, List[Tuple[int, int]]] = {}
        for i_a, (n_a, l_a, m_a) in enumerate(states_A):
            shell_indices_A.setdefault(n_a, []).append((len(shell_indices_A.get(n_a, [])), i_a))
        for i_b, (n_b, l_b, m_b) in enumerate(states_B):
            shell_indices_B.setdefault(n_b, []).append((len(shell_indices_B.get(n_b, [])), i_b))

        rows: List[int] = []
        cols: List[int] = []
        data: List[float] = []

        for n_shell in range(1, n_max + 1):
            if n_shell not in shell_indices_A or n_shell not in shell_indices_B:
                continue

            D_sym = D_blocks[n_shell]
            f = f_factors[n_shell]
            scale = -(Z_eff / p0) * f

            idxs_A = shell_indices_A[n_shell]
            idxs_B = shell_indices_B[n_shell]

            for local_a, global_a in idxs_A:
                for local_b, global_b in idxs_B:
                    h_val = scale * D_sym[local_a, local_b]

                    if abs(h_val) >= self.threshold:
                        j_b = global_b + nA
                        # A->B coupling
                        rows.append(global_a)
                        cols.append(j_b)
                        data.append(h_val)
                        # B->A coupling (Hermitian)
                        rows.append(j_b)
                        cols.append(global_a)
                        data.append(h_val)

        H_cross = csr_matrix((data, (rows, cols)), shape=(n, n))
        max_val = max(abs(v) for v in data) if data else 0.0
        print(f"[SW] Cross-atom H1: {H_cross.nnz} nnz, "
              f"|max|={max_val:.6f} Ha, "
              f"f(1,gamma)={f_factors.get(1, 0):.4f}")
        return H_cross

    def _cross_n_overlap_bridges(self, R: float) -> csr_matrix:
        """
        Cross-n inter-atomic coupling via SW-consistent overlap formula.

        For orbital pairs (i_A, j_B) where n_A != n_B, compute:
            H1[i_A, j_B] = -K * (Z_eff/p0) * sin(gamma) * S(i_A, j_B)

        This uses the same scale as the Shibuya-Wulfman D-matrix formula
        for same-n pairs, but replaces the D-matrix element (which is zero
        for cross-n) with the STO overlap integral. The sin(gamma) form
        factor ensures proper R->inf decay.

        K (cross_n_K) scales the coupling strength relative to same-n
        D-matrix coupling (default 1.0 = same scale as D-matrix).

        Only s-orbital (l=0) pairs are computed (dominant for sigma bonds;
        compute_overlap_element currently supports l=0 only).

        Parameters
        ----------
        R : float
            Internuclear distance in bohr.

        Returns
        -------
        csr_matrix
            Shape (n_spatial, n_spatial). Cross-n bridge coupling matrix.
        """
        from scipy.sparse import lil_matrix
        from .wigner_so4 import bond_angle
        from .shibuya_wulfman import sw_form_factor

        nA = self._n_spatial_A
        n_total = self._n_spatial
        states_A = self._li_A.lattice.states
        states_B = self._li_B.lattice.states
        Z_orb_A = self._Z_orb_A
        Z_orb_B = self._Z_orb_B
        K = self.cross_n_K
        Z_A_f = float(self.Z_A)
        Z_B_f = float(self.Z_B)

        # Use same p0 and Z_eff as D-matrix path
        p0 = getattr(self, '_dmatrix_p0', np.sqrt(Z_A_f**2 + Z_B_f**2))
        gamma = bond_angle(R, p0)
        Z_eff = (Z_A_f + Z_B_f) / 2.0
        f_gamma = sw_form_factor(1, gamma)  # sin(gamma)

        # SW-consistent scale: -(Z_eff/p0) * sin(gamma), same as D-matrix
        sw_scale = -(Z_eff / p0) * f_gamma

        H_bridge = lil_matrix((n_total, n_total))
        n_cross_n = 0
        max_abs = 0.0

        for i_a, (n_a, l_a, m_a) in enumerate(states_A):
            if l_a != 0:
                continue

            for i_b, (n_b, l_b, m_b) in enumerate(states_B):
                if l_b != 0:
                    continue
                # Skip same-n pairs (handled by D-matrix)
                if n_a == n_b:
                    continue

                S_ab = compute_overlap_element(
                    n_a, l_a, n_b, l_b,
                    Z_orb_A, Z_orb_B, R
                )
                if abs(S_ab) < 1e-12:
                    continue

                H_ij = K * sw_scale * S_ab

                if abs(H_ij) < self.threshold:
                    continue

                j_b = nA + i_b
                H_bridge[i_a, j_b] = H_ij
                H_bridge[j_b, i_a] = H_ij
                n_cross_n += 1
                max_abs = max(max_abs, abs(H_ij))

        result = H_bridge.tocsr()
        print(f"[Cross-n bridges] {n_cross_n} s-orbital pairs, "
              f"K={K}, scale={sw_scale:.6f}, |max|={max_abs:.6f} Ha")
        return result

    def _recompute_slater_f0_scaled(
        self, li: 'LatticeIndex', Z_orb: float
    ) -> Tuple[np.ndarray, Dict[Tuple[int, int, int, int], float]]:
        """
        Recompute Slater F0 integrals with scaled orbital exponent Z_orb.

        When zeta != 1.0, the same-atom V_ee must use orbitals with
        effective charge Z_orb instead of Z_nuclear.

        Returns (vee_matrix, eri_dict) for the atom's spatial states.
        """
        states = li.lattice.states
        n_sp = li.lattice.num_states
        unique_nl = sorted(set((n, l) for n, l, m in states))

        # Compute F0 integrals with Z_orb
        f0_cache: Dict[Tuple[int, int, int, int], float] = {}
        for n1, l1 in unique_nl:
            for n2, l2 in unique_nl:
                key = (n1, l1, n2, l2)
                key_rev = (n2, l2, n1, l1)
                if key in f0_cache or key_rev in f0_cache:
                    continue
                val = _compute_single_f0(n1, l1, n2, l2, Z_orb)
                f0_cache[key] = val
                f0_cache[key_rev] = val

        # Fill vee matrix and ERI dict
        vee = np.zeros((n_sp, n_sp))
        eri: Dict[Tuple[int, int, int, int], float] = {}
        for i in range(n_sp):
            ni, li_val, _ = states[i]
            for j in range(i, n_sp):
                nj, lj, _ = states[j]
                v = f0_cache[(ni, li_val, nj, lj)]
                vee[i, j] = v
                vee[j, i] = v

        # Build ERI table matching LatticeIndex format
        if hasattr(li, '_eri'):
            for key_orig in li._eri:
                a, b, c, d = key_orig
                if a < n_sp and b < n_sp and c < n_sp and d < n_sp:
                    na, la, _ = states[a]
                    nb, lb, _ = states[b]
                    nc, lc, _ = states[c]
                    nd, ld, _ = states[d]
                    # Recompute this ERI entry with scaled Z
                    # For F0 (direct): (ab|ab) = F0(a,b)
                    if a == c and b == d:
                        eri[key_orig] = f0_cache[(na, la, nb, lb)]
                    elif a == d and b == c:
                        # Exchange: (ab|ba) — use same F0 for Mulliken approx
                        eri[key_orig] = f0_cache[(na, la, nb, lb)]
                    else:
                        # Keep original for non-diagonal ERIs
                        eri[key_orig] = li._eri[key_orig]

        return vee, eri

    def _build_molecular_vee(self) -> None:
        """
        Build combined V_ee matrix and ERI table.

        Same-atom Slater integrals are exact (cached). When zeta scaling is
        active, same-atom integrals are recomputed with Z_orb = zeta * Z.
        Cross-atom two-electron integrals use Fourier convolution of
        momentum-space densities with sin(qR)/(qR) structure factor.

        Controlled by self.cross_atom_vee:
        - True: all (n,l) orbital pairs (monopole J + s-s Mulliken K)
        - 's_only': s-orbital pairs only (v0.9.11 behavior)
        - False: no cross-atom V_ee (diagnostic mode)
        """
        n = self._n_spatial
        nA = self._n_spatial_A

        self._vee_matrix = np.zeros((n, n))

        # Same-atom V_ee: use original if zeta=1.0, else recompute
        self._eri: Dict[Tuple[int, int, int, int], float] = {}

        if not self._ghost_A:
            if abs(self.zeta_A - 1.0) < 1e-12:
                self._vee_matrix[:nA, :nA] = self._li_A._vee_matrix
                if hasattr(self._li_A, '_eri'):
                    for key, val in self._li_A._eri.items():
                        self._eri[key] = val
            else:
                vee_A, eri_A = self._recompute_slater_f0_scaled(
                    self._li_A, self._Z_orb_A)
                self._vee_matrix[:nA, :nA] = vee_A
                for key, val in eri_A.items():
                    self._eri[key] = val

        if not self._ghost_B:
            if abs(self.zeta_B - 1.0) < 1e-12:
                self._vee_matrix[nA:, nA:] = self._li_B._vee_matrix
                if hasattr(self._li_B, '_eri'):
                    for (a, b, c, d), val in self._li_B._eri.items():
                        self._eri[(a + nA, b + nA, c + nA, d + nA)] = val
            else:
                vee_B, eri_B = self._recompute_slater_f0_scaled(
                    self._li_B, self._Z_orb_B)
                self._vee_matrix[nA:, nA:] = vee_B
                for (a, b, c, d), val in eri_B.items():
                    self._eri[(a + nA, b + nA, c + nA, d + nA)] = val

        # Cross-atom V_ee (skip for ghost atoms or if disabled)
        n_cross = 0
        if not self._ghost_A and not self._ghost_B and self.cross_atom_vee:
            n_cross = self._build_cross_atom_vee()

        n_eri = len(self._eri)
        m4 = n ** 4
        density = n_eri / m4 if m4 > 0 else 0.0
        print(f"[MolecularLatticeIndex] V_ee: {n_eri} ERI entries "
              f"({n_cross} cross-atom J+K), "
              f"density {n_eri}/{m4} = {density:.2%}")

    def _build_cross_atom_vee(self) -> int:
        """
        Compute cross-atom two-electron integrals.

        J_AB uses Fourier convolution of spherically averaged density
        form factors (monopole term, exact for s-orbitals, ~90% for l>0):
            J_AB(R) = (2/pi) int_0^inf rho_A(q) rho_B(q) sin(qR)/(qR) dq

        K_AB uses the Mulliken approximation (s-orbital pairs only):
            K(aA, bB; R) = S(aA, bB)^2 * [F0(a,a; ZA) + F0(b,b; ZB)] / 2

        Controlled by self.cross_atom_vee:
        - True: J for all (n,l) pairs + K for s-s pairs
        - 's_only': J and K for s-s pairs only (v0.9.11 behavior)

        Uses batch computation (compute_cross_atom_J_batch) to evaluate all
        unique J integrals on a shared q-grid, avoiding nested numerical
        integration. ~100x faster than per-integral quad for l>0 orbitals.

        Returns
        -------
        int
            Number of cross-atom ERI entries added (J + K).
        """
        nA = self._n_spatial_A
        states_A = self._li_A.lattice.states
        states_B = self._li_B.lattice.states
        R = self.R
        s_only = (self.cross_atom_vee == 's_only')
        count_J = 0
        count_K = 0

        # Collect unique (na, la, nb, lb) pairs needed
        unique_pairs: set = set()
        for n_a, l_a, m_a in states_A:
            if s_only and l_a > 0:
                continue
            for n_b, l_b, m_b in states_B:
                if s_only and l_b > 0:
                    continue
                unique_pairs.add((n_a, l_a, n_b, l_b))

        # Batch-compute all J integrals
        j_cache = compute_cross_atom_J_batch(
            list(unique_pairs), R,
            self._Z_orb_A, self._Z_orb_B
        )

        # Cache S by (n_a, n_b) — overlap only for s-orbitals
        s_cache: Dict[Tuple[int, int], float] = {}

        for i_a, (n_a, l_a, m_a) in enumerate(states_A):
            if s_only and l_a > 0:
                continue
            for i_b, (n_b, l_b, m_b) in enumerate(states_B):
                if s_only and l_b > 0:
                    continue
                j_b = i_b + nA

                # --- Direct Coulomb J (monopole term) ---
                j_ab = j_cache[(n_a, l_a, n_b, l_b)]

                if abs(j_ab) > 1e-15:
                    self._vee_matrix[i_a, j_b] = j_ab
                    self._vee_matrix[j_b, i_a] = j_ab
                    self._eri[(i_a, j_b, i_a, j_b)] = j_ab
                    self._eri[(j_b, i_a, j_b, i_a)] = j_ab
                    count_J += 2

                # --- Exchange K via Mulliken (s-orbital pairs only) ---
                if l_a == 0 and l_b == 0 and abs(j_ab) > 1e-15:
                    s_key = (n_a, n_b)
                    if s_key not in s_cache:
                        s_cache[s_key] = compute_overlap_element(
                            n_a, 0, n_b, 0,
                            self._Z_orb_A, self._Z_orb_B, R
                        )
                    s_ab = s_cache[s_key]

                    # F0 self-repulsion: use zeta-scaled values from ERI table
                    f0_aa = self._eri.get((i_a, i_a, i_a, i_a), 0.0)
                    f0_bb = self._eri.get(
                        (i_b + nA, i_b + nA, i_b + nA, i_b + nA), 0.0)
                    k_ab = s_ab * s_ab * (f0_aa + f0_bb) / 2.0

                    if abs(k_ab) > 1e-15:
                        self._eri[(i_a, j_b, j_b, i_a)] = k_ab
                        self._eri[(j_b, i_a, i_a, j_b)] = k_ab
                        count_K += 2

        print(f"[MolecularLatticeIndex] Cross-atom V_ee: "
              f"{count_J} Coulomb J, {count_K} exchange K entries "
              f"(mode={'s_only' if s_only else 'all_l'})")
        return count_J + count_K

    def _compute_overlap_matrix(self) -> np.ndarray:
        """
        Compute the overlap matrix S between all spatial orbitals.

        Same-atom blocks are identity (orthonormal hydrogenic basis).
        Cross-atom blocks use Fourier-Bessel overlap:
            S_{ij} = (2/pi) int_0^inf g_A(q) g_B(q) sin(qR)/(qR) q^2 dq

        Only s-orbital (l=0) cross-atom pairs are computed; l>0 cross-atom
        overlaps are set to zero (consistent with l=0 restriction on
        cross-atom V_ee and cross-nuclear attraction).

        Returns
        -------
        S : np.ndarray, shape (n_spatial, n_spatial)
            Overlap matrix (symmetric positive semi-definite)
        """
        n = self._n_spatial
        nA = self._n_spatial_A
        S = np.eye(n)

        # Skip cross-atom overlaps if either atom is a ghost
        if self._ghost_A or self._ghost_B:
            return S

        states_A = self._li_A.lattice.states
        states_B = self._li_B.lattice.states

        # Cache by (na, nb) since l=0, m doesn't matter for s-orbitals
        s_cache: Dict[Tuple[int, int], float] = {}

        for i_a, (n_a, l_a, m_a) in enumerate(states_A):
            if l_a > 0:
                continue
            for i_b, (n_b, l_b, m_b) in enumerate(states_B):
                if l_b > 0:
                    continue
                j_b = i_b + nA

                cache_key = (n_a, n_b)
                if cache_key not in s_cache:
                    s_cache[cache_key] = compute_overlap_element(
                        n_a, 0, n_b, 0,
                        self._Z_orb_A, self._Z_orb_B, self.R
                    )
                s_ab = s_cache[cache_key]

                if abs(s_ab) < 1e-15:
                    continue

                S[i_a, j_b] = s_ab
                S[j_b, i_a] = s_ab

        return S

    def _lowdin_orthogonalize(
        self, S: np.ndarray, eigenvalue_threshold: float = 1e-6
    ) -> None:
        """
        Apply Lowdin symmetric orthogonalization S^{-1/2} to all integrals.

        Transforms the one-electron Hamiltonian and two-electron integrals
        from the non-orthogonal atomic orbital basis to the orthogonalized
        Lowdin basis, removing inter-center linear dependencies that cause
        BSSE.

        The transformation is:
            H1' = X^T H1 X
            (pq|rs)' = sum_{abcd} X_{ap} X_{bq} (ab|cd) X_{cr} X_{ds}

        where X = S^{-1/2} computed via eigendecomposition:
            S = U Lambda U^T  =>  X = U Lambda^{-1/2} U^T

        Eigenvalues below eigenvalue_threshold are dropped (near-linear
        dependencies), reducing the effective basis size.

        Parameters
        ----------
        S : np.ndarray, shape (n_spatial, n_spatial)
            Overlap matrix
        eigenvalue_threshold : float
            Eigenvalues below this are dropped (default 1e-6)
        """
        n = self._n_spatial

        # Eigendecompose overlap matrix
        eigvals, eigvecs = np.linalg.eigh(S)

        # Drop near-zero eigenvalues (linear dependencies)
        mask = eigvals > eigenvalue_threshold
        n_kept = int(mask.sum())
        n_dropped = n - n_kept

        if n_dropped > 0:
            print(f"[Lowdin] Dropping {n_dropped} near-dependent basis functions "
                  f"(threshold={eigenvalue_threshold:.1e})")

        # Build X = S^{-1/2} using only significant eigenvalues
        # For full basis: X = U @ diag(1/sqrt(lambda)) @ U^T
        # With truncation: X is n×n_kept (rectangular if any dropped)
        lam_inv_sqrt = 1.0 / np.sqrt(eigvals[mask])
        U_kept = eigvecs[:, mask]

        # Full S^{-1/2} matrix (square, even with truncation the result
        # projects back into the original n-dimensional space)
        X = U_kept @ np.diag(lam_inv_sqrt) @ U_kept.T

        # --- Transform H1 ---
        H1_dense = self._H1_spatial.toarray()
        H1_orth = X.T @ H1_dense @ X
        self._h1_diag = np.diag(H1_orth).copy()
        self._H1_spatial = csr_matrix(H1_orth)

        # Rebuild off-diagonal dict
        self._h1_offdiag = {i: [] for i in range(n)}
        for r in range(n):
            for c in range(n):
                if r != c and abs(H1_orth[r, c]) >= self.threshold:
                    self._h1_offdiag[r].append((c, float(H1_orth[r, c])))

        # --- Transform V_ee (4-index) ---
        # Build dense 4D ERI array from sparse dict
        eri_4d = np.zeros((n, n, n, n))
        for (a, b, c, d), val in self._eri.items():
            eri_4d[a, b, c, d] = val

        # 4-index transformation: (pq|rs) = sum_{abcd} X_{ap} X_{bq} (ab|cd) X_{cr} X_{ds}
        # Done as 4 sequential contractions for efficiency
        tmp = np.einsum('ap,abcd->pbcd', X, eri_4d)
        tmp = np.einsum('bq,pbcd->pqcd', X, tmp)
        tmp = np.einsum('cr,pqcd->pqrd', X, tmp)
        eri_orth = np.einsum('ds,pqrd->pqrs', X, tmp)

        # Rebuild ERI dict and vee_matrix from transformed integrals
        self._eri = {}
        self._vee_matrix = np.zeros((n, n))
        for p in range(n):
            for q in range(n):
                val_J = eri_orth[p, q, p, q]
                if abs(val_J) > 1e-15:
                    self._vee_matrix[p, q] = val_J
                    self._eri[(p, q, p, q)] = val_J
                # Also store exchange-type integrals
                val_K = eri_orth[p, q, q, p]
                if abs(val_K) > 1e-15 and (p, q, q, p) != (p, q, p, q):
                    self._eri[(p, q, q, p)] = val_K

        # Store all non-negligible ERIs for full Slater-Condon
        for a in range(n):
            for b in range(n):
                for c in range(n):
                    for d in range(n):
                        val = eri_orth[a, b, c, d]
                        if abs(val) > 1e-15:
                            self._eri[(a, b, c, d)] = val

        self._overlap_matrix = S
        self._lowdin_X = X
        n_eri = len(self._eri)
        print(f"[Lowdin] Orthogonalized: {n_kept}/{n} basis functions kept, "
              f"{n_eri} ERI entries")

    def _enumerate_sd_basis(self) -> None:
        """Enumerate all N-electron Slater determinants over combined basis."""
        t0 = time.perf_counter()
        self.sd_basis: List[Tuple[int, ...]] = list(
            combinations(range(self.n_sp), self.n_electrons)
        )
        self.n_sd: int = len(self.sd_basis)
        self._sd_index: Dict[Tuple[int, ...], int] = {
            sd: i for i, sd in enumerate(self.sd_basis)
        }
        elapsed = time.perf_counter() - t0
        print(
            f"[MolecularLatticeIndex] {self._n_spatial} spatial, "
            f"{self.n_sp} spin-orbitals, {self.n_sd:,} SDs "
            f"(enumerated in {elapsed:.3f}s)"
        )

    # ------------------------------------------------------------------
    # ERI access (duck-type compatible with LatticeIndex)
    # ------------------------------------------------------------------

    def _get_eri(self, a: int, b: int, c: int, d: int) -> float:
        """Look up two-electron integral <ab|cd> from ERI table."""
        return self._eri.get((a, b, c, d), 0.0)

    # ------------------------------------------------------------------
    # Hamiltonian assembly (matrix method, for small systems)
    # ------------------------------------------------------------------

    def assemble_hamiltonian(self) -> csr_matrix:
        """
        Assemble the molecular FCI Hamiltonian via Slater-Condon rules.

        Uses full two-electron integrals (same-atom approximation).
        """
        if self.vee_method == 'slater_full':
            return self._assemble_hamiltonian_full()
        return self._assemble_hamiltonian_diag_vee()

    def _assemble_hamiltonian_full(self) -> csr_matrix:
        """Full Slater-Condon assembly with off-diagonal two-electron integrals."""
        t0 = time.perf_counter()

        diag_rows: List[int] = []
        diag_vals: List[float] = []
        off_rows: List[int] = []
        off_cols: List[int] = []
        off_vals: List[float] = []

        # Diagonal
        for I, sd_I in enumerate(self.sd_basis):
            h_diag = 0.0
            for p in sd_I:
                h_diag += self._h1_diag[p >> 1]
            n = len(sd_I)
            for i in range(n):
                pi = sd_I[i]
                sp_i = pi >> 1
                sig_i = pi & 1
                for j in range(i + 1, n):
                    pj = sd_I[j]
                    sp_j = pj >> 1
                    h_diag += self._get_eri(sp_i, sp_j, sp_i, sp_j)
                    if (pj & 1) == sig_i:
                        h_diag -= self._get_eri(sp_i, sp_j, sp_j, sp_i)
            if abs(h_diag) >= self.threshold:
                diag_rows.append(I)
                diag_vals.append(h_diag)

        t_diag = time.perf_counter() - t0

        # Off-diagonal (single + double excitations)
        sd_sets = [frozenset(sd) for sd in self.sd_basis]
        for I in range(self.n_sd):
            sd_I = self.sd_basis[I]
            set_I = sd_sets[I]
            for J in range(I + 1, self.n_sd):
                set_J = sd_sets[J]
                diff_I = set_I - set_J
                diff_J = set_J - set_I
                n_diff = len(diff_I)
                if n_diff == 0 or n_diff > 2:
                    continue
                if n_diff == 1:
                    p = next(iter(diff_I))
                    r = next(iter(diff_J))
                    if (p & 1) != (r & 1):
                        continue
                    kp = sd_I.index(p)
                    sp_p = p >> 1
                    sp_r = r >> 1
                    me = self._H1_spatial[sp_p, sp_r]
                    for q in sd_I:
                        if q == p:
                            continue
                        sp_q = q >> 1
                        me += self._get_eri(sp_p, sp_q, sp_r, sp_q)
                        if (p & 1) == (q & 1):
                            me -= self._get_eri(sp_p, sp_q, sp_q, sp_r)
                    if abs(me) < self.threshold:
                        continue
                    phase = self._compute_phase(sd_I, kp, r)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * me)
                else:
                    removed = sorted(diff_I)
                    added = sorted(diff_J)
                    p, q = removed
                    r, s = added
                    if (p & 1) + (q & 1) != (r & 1) + (s & 1):
                        continue
                    sp_p, sp_q = p >> 1, q >> 1
                    sp_r, sp_s = r >> 1, s >> 1
                    me = 0.0
                    if (p & 1) == (r & 1) and (q & 1) == (s & 1):
                        me += self._get_eri(sp_p, sp_q, sp_r, sp_s)
                    if (p & 1) == (s & 1) and (q & 1) == (r & 1):
                        me -= self._get_eri(sp_p, sp_q, sp_s, sp_r)
                    if abs(me) < self.threshold:
                        continue
                    kp = sd_I.index(p)
                    kq = sd_I.index(q)
                    phase = self._compute_double_phase(sd_I, kp, kq, r, s)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * me)

        H_diag = csr_matrix(
            (diag_vals, (diag_rows, diag_rows)), shape=(self.n_sd, self.n_sd))
        H_upper = csr_matrix(
            (off_vals, (off_rows, off_cols)), shape=(self.n_sd, self.n_sd))
        H = H_upper + H_upper.T + H_diag

        elapsed = time.perf_counter() - t0
        print(f"[MolecularLatticeIndex] H assembled: shape={H.shape}, "
              f"nnz={H.nnz:,}, time={elapsed:.3f}s")
        return H.tocsr()

    def _assemble_hamiltonian_diag_vee(self) -> csr_matrix:
        """Diagonal-only V_ee assembly (for non-slater_full methods)."""
        t0 = time.perf_counter()
        diag_rows: List[int] = []
        diag_vals: List[float] = []
        off_rows: List[int] = []
        off_cols: List[int] = []
        off_vals: List[float] = []

        for I, sd in enumerate(self.sd_basis):
            occupied_set = set(sd)
            h_diag = 0.0
            for p in sd:
                h_diag += self._h1_diag[p >> 1]
            # Diagonal V_ee
            n = len(sd)
            for i in range(n):
                for j in range(i + 1, n):
                    h_diag += self._vee_matrix[sd[i] >> 1, sd[j] >> 1]
            if abs(h_diag) >= self.threshold:
                diag_rows.append(I)
                diag_vals.append(h_diag)
            for kp, p in enumerate(sd):
                sigma_p = p & 1
                sp_p = p >> 1
                for r_sp, h_val in self._h1_offdiag[sp_p]:
                    if abs(h_val) < self.threshold:
                        continue
                    r = (r_sp << 1) | sigma_p
                    if r in occupied_set:
                        continue
                    new_sd = tuple(sorted(sd[:kp] + (r,) + sd[kp + 1:]))
                    J = self._sd_index.get(new_sd)
                    if J is None or J <= I:
                        continue
                    phase = self._compute_phase(sd, kp, r)
                    off_rows.append(I)
                    off_cols.append(J)
                    off_vals.append(phase * h_val)

        H_diag = csr_matrix(
            (diag_vals, (diag_rows, diag_rows)), shape=(self.n_sd, self.n_sd))
        H_upper = csr_matrix(
            (off_vals, (off_rows, off_cols)), shape=(self.n_sd, self.n_sd))
        H = H_upper + H_upper.T + H_diag
        elapsed = time.perf_counter() - t0
        print(f"[MolecularLatticeIndex] H assembled: shape={H.shape}, "
              f"nnz={H.nnz:,}, time={elapsed:.3f}s")
        return H.tocsr()

    # ------------------------------------------------------------------
    # Energy Decomposition (diagnostic instrumentation, v0.9.24)
    # ------------------------------------------------------------------

    def _build_1rdm_diagonal(self, civec: np.ndarray) -> np.ndarray:
        """
        Build diagonal of the one-particle reduced density matrix.

        gamma_pp = sum_I |c_I|^2 * n_p(I)

        where n_p(I) is the occupation of spatial orbital p in SD I
        (0, 1, or 2 counting both spins).

        Parameters
        ----------
        civec : np.ndarray, shape (n_sd,)
            FCI ground-state coefficient vector.

        Returns
        -------
        gamma_diag : np.ndarray, shape (n_spatial,)
            Diagonal of spin-summed 1-RDM.
        """
        n_sp = self._n_spatial
        gamma_diag = np.zeros(n_sp)
        c_sq = civec * civec
        for I, sd in enumerate(self.sd_basis):
            for occ in sd:
                gamma_diag[occ >> 1] += c_sq[I]
        return gamma_diag

    def _compute_h1_expectation(self, civec: np.ndarray) -> float:
        """
        Compute <psi|H1|psi> via sparse single-excitation walk.

        Avoids building the full 1-RDM. Instead, iterates over SDs and
        their single excitations through nonzero H1 off-diagonal elements.

        Parameters
        ----------
        civec : np.ndarray, shape (n_sd,)
            FCI coefficient vector.

        Returns
        -------
        float
            One-electron energy expectation value <H1>.
        """
        # Diagonal H1 contribution
        h1_diag_val = 0.0
        c_sq = civec * civec
        for I, sd in enumerate(self.sd_basis):
            for occ in sd:
                h1_diag_val += c_sq[I] * self._h1_diag[occ >> 1]

        # Off-diagonal H1 contribution (single excitations only)
        h1_offdiag_val = 0.0
        sd_index = self._sd_index
        for I, sd in enumerate(self.sd_basis):
            c_I = civec[I]
            if abs(c_I) < 1e-15:
                continue
            occ_set = frozenset(sd)
            for kp, p in enumerate(sd):
                sigma_p = p & 1
                sp_p = p >> 1
                for r_sp, h_val in self._h1_offdiag[sp_p]:
                    r = (r_sp << 1) | sigma_p
                    if r in occ_set:
                        continue
                    new_sd = tuple(sorted(sd[:kp] + (r,) + sd[kp + 1:]))
                    J = sd_index.get(new_sd)
                    if J is None:
                        continue
                    c_J = civec[J]
                    if abs(c_J) < 1e-15:
                        continue
                    phase = self._compute_phase(sd, kp, r)
                    h1_offdiag_val += c_I * c_J * phase * h_val

        return h1_diag_val + h1_offdiag_val

    def _compute_bridge_expectation(self, civec: np.ndarray) -> Tuple[float, int, float]:
        """
        Compute <psi|V_bridge|psi> using bridge-only adjacency.

        The bridge contribution is the off-diagonal inter-atomic part of
        the combined adjacency matrix, scaled by kinetic_scale.

        Returns
        -------
        v_bridge : float
            Bridge hopping expectation value.
        n_bridge_active : int
            Number of active bridge connections.
        max_bridge_elem : float
            Largest bridge matrix element magnitude.
        """
        from scipy.sparse import block_diag as sp_block_diag

        A_intra = sp_block_diag(
            [self._li_A.lattice.adjacency, self._li_B.lattice.adjacency],
            format='csr'
        )
        A_full = self._adjacency_combined
        # Bridge adjacency = full - intra (only inter-atomic connections)
        A_bridge = (A_full - A_intra).tocoo()

        n_bridge_active = A_bridge.nnz // 2
        if n_bridge_active == 0:
            return 0.0, 0, 0.0

        # Build bridge off-diagonal dict (spatial index → list of (target, weight))
        n = self._n_spatial
        bridge_offdiag: Dict[int, List[Tuple[int, float]]] = {i: [] for i in range(n)}
        max_elem = 0.0
        for r, c, v in zip(A_bridge.row, A_bridge.col, A_bridge.data):
            if r != c and abs(v) >= 1e-15:
                h_val = self.kinetic_scale * (-v)
                bridge_offdiag[r].append((c, h_val))
                max_elem = max(max_elem, abs(h_val))

        # Single-excitation walk for bridge-only operator
        v_bridge = 0.0
        sd_index = self._sd_index
        for I, sd in enumerate(self.sd_basis):
            c_I = civec[I]
            if abs(c_I) < 1e-15:
                continue
            occ_set = frozenset(sd)
            for kp, p in enumerate(sd):
                sigma_p = p & 1
                sp_p = p >> 1
                for r_sp, h_val in bridge_offdiag[sp_p]:
                    r = (r_sp << 1) | sigma_p
                    if r in occ_set:
                        continue
                    new_sd = tuple(sorted(sd[:kp] + (r,) + sd[kp + 1:]))
                    J = sd_index.get(new_sd)
                    if J is None:
                        continue
                    c_J = civec[J]
                    if abs(c_J) < 1e-15:
                        continue
                    phase = self._compute_phase(sd, kp, r)
                    v_bridge += c_I * c_J * phase * h_val

        return v_bridge, n_bridge_active, max_elem

    def decompose_energy(
        self, civec: np.ndarray, E_total: float
    ) -> Dict[str, float]:
        """
        Decompose the FCI ground-state energy into physical components.

        Uses the diagonal 1-RDM for diagonal operators and sparse
        single-excitation walks for off-diagonal operators (T, bridges).
        V_ee is obtained as a residual: V_ee = E_total - V_NN - <H1>.

        Components
        ----------
        T : kinetic energy (graph Laplacian, intra-atom only)
        V_nA : nuclear attraction of A electrons to nucleus A
        V_nB : nuclear attraction of B electrons to nucleus B
        V_cross_A : cross-nuclear attraction of A electrons to nucleus B
        V_cross_B : cross-nuclear attraction of B electrons to nucleus A
        V_bridge : bridge hopping contribution (off-diagonal A-B block)
        V_ee : electron-electron repulsion (residual)
        V_NN : nuclear-nuclear repulsion

        Parameters
        ----------
        civec : np.ndarray, shape (n_sd,)
            FCI ground-state coefficient vector.
        E_total : float
            Total FCI energy (electronic + V_NN) from compute_ground_state.

        Returns
        -------
        components : dict
            Keys: 'T', 'V_nA', 'V_nB', 'V_cross_A', 'V_cross_B',
                  'V_bridge', 'V_ee', 'V_NN', 'E_total', 'E_check'
        """
        n = self._n_spatial
        nA = self._n_spatial_A

        # Diagonal 1-RDM: gamma_pp = sum_I |c_I|^2 * n_p(I)
        gamma_diag = self._build_1rdm_diagonal(civec)

        # --- Nuclear attraction V_nA: -Z_A^2/(2n^2) for A-block ---
        VnA_diag = np.zeros(n)
        if not self._ghost_A:
            for i, (ni, li, mi) in enumerate(self._li_A.lattice.states):
                VnA_diag[i] = -float(self.Z_A)**2 / (2.0 * ni**2)
        VnA_val = np.dot(VnA_diag, gamma_diag)

        # --- Nuclear attraction V_nB: -Z_B^2/(2n^2) for B-block ---
        VnB_diag = np.zeros(n)
        if not self._ghost_B:
            for j, (nj, lj, mj) in enumerate(self._li_B.lattice.states):
                VnB_diag[nA + j] = -float(self.Z_B)**2 / (2.0 * nj**2)
        VnB_val = np.dot(VnB_diag, gamma_diag)

        # --- Cross-nuclear V_cross_A and V_cross_B ---
        # Uses same method as _build_molecular_h1 (exact or fourier).
        Vcross_diag = np.zeros(n)
        self._apply_cross_nuclear_diagonal(Vcross_diag)
        Vcross_A_diag = np.zeros(n)
        Vcross_A_diag[:nA] = Vcross_diag[:nA]
        Vcross_B_diag = np.zeros(n)
        Vcross_B_diag[nA:] = Vcross_diag[nA:]
        Vcross_A_val = np.dot(Vcross_A_diag, gamma_diag)
        Vcross_B_val = np.dot(Vcross_B_diag, gamma_diag)

        # --- H1 total via sparse single-excitation walk ---
        H1_total = self._compute_h1_expectation(civec)

        # --- Bridge hopping (separate excitation walk) ---
        V_bridge_val, n_bridge_active, max_bridge_elem = \
            self._compute_bridge_expectation(civec)

        # --- T = <H1> - diagonal_components - V_bridge ---
        diag_sum = VnA_val + VnB_val + Vcross_A_val + Vcross_B_val
        T_val = H1_total - diag_sum - V_bridge_val

        # --- V_ee as residual: E_total = <H1> + <V_ee> + V_NN ---
        E_elec = E_total - self.V_NN
        V_ee_val = E_elec - H1_total

        return {
            'T': float(T_val),
            'V_nA': float(VnA_val),
            'V_nB': float(VnB_val),
            'V_cross_A': float(Vcross_A_val),
            'V_cross_B': float(Vcross_B_val),
            'V_bridge': float(V_bridge_val),
            'V_ee': float(V_ee_val),
            'V_NN': float(self.V_NN),
            'E_total': float(E_total),
            'E_check': float(H1_total + V_ee_val + self.V_NN),
            'H1_total': float(H1_total),
            'H1_diag_sum': float(diag_sum),
            'n_bridge_active': n_bridge_active,
            'max_bridge_elem': float(max_bridge_elem),
        }

    # ------------------------------------------------------------------
    # Solver
    # ------------------------------------------------------------------

    def compute_ground_state(
        self, n_states: int = 1
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute the molecular ground state via FCI.

        Adds nuclear repulsion V_NN = Z_A * Z_B / R to eigenvalues.

        Returns
        -------
        eigvals : np.ndarray, shape (n_states,)
            Total energies (electronic + nuclear repulsion)
        eigvecs : np.ndarray, shape (n_sd, n_states)
        """
        method = self.fci_method
        if method == 'auto':
            method = 'direct' if self.n_sd >= 5000 else 'matrix'

        if method == 'direct' and self.vee_method == 'slater_full':
            from .direct_ci import DirectCISolver
            solver = DirectCISolver(self)
            eigvals, eigvecs = solver.solve(n_states=n_states)
        else:
            H = self.assemble_hamiltonian()
            k = min(n_states, self.n_sd - 2)
            rng = np.random.RandomState(42)
            v0 = rng.randn(H.shape[0])
            eigvals, eigvecs = eigsh(H, k=k, which="SA", v0=v0)
            order = np.argsort(eigvals)
            eigvals, eigvecs = eigvals[order], eigvecs[:, order]

        # Add nuclear repulsion
        eigvals = eigvals + self.V_NN
        print(f"[MolecularLatticeIndex] E_elec = {eigvals[0] - self.V_NN:.6f} Ha, "
              f"V_NN = {self.V_NN:.6f} Ha, E_total = {eigvals[0]:.6f} Ha")
        return eigvals, eigvecs


# ---------------------------------------------------------------------------
# Boys-Bernardi Counterpoise Correction for BSSE
# ---------------------------------------------------------------------------

def compute_bsse_correction(
    Z_A: int,
    Z_B: int,
    nmax_A: int,
    nmax_B: int,
    R: float,
    n_electrons_A: int,
    n_electrons_B: int,
    vee_method: str = 'slater_full',
    fci_method: str = 'auto',
    orthogonalize: bool = False,
    zeta_A: float = 1.0,
    zeta_B: float = 1.0,
) -> dict:
    """
    Boys-Bernardi counterpoise correction for Basis Set Superposition Error.

    Computes four energies:
      E_A_own   = atom A with its own basis (nmax_A)
      E_B_own   = atom B with its own basis (nmax_B)
      E_A_ghost = atom A in full molecular basis (nmax_A + ghost nmax_B)
      E_B_ghost = atom B in full molecular basis (ghost nmax_A + nmax_B)

    BSSE = (E_A_ghost - E_A_own) + (E_B_ghost - E_B_own)
         = artificial lowering due to basis borrowing (negative number)

    When zeta_A/zeta_B != 1.0, ghost atom fragment calculations use the
    same zeta values as the supermolecule to ensure consistent basis space.

    Parameters
    ----------
    Z_A, Z_B : int
        Nuclear charges
    nmax_A, nmax_B : int
        Basis sizes for each atom
    R : float
        Internuclear distance in Bohr (defines ghost orbital placement)
    n_electrons_A, n_electrons_B : int
        Electron counts for each atom
    vee_method : str
        Two-electron integral method (default 'slater_full')
    fci_method : str
        FCI assembly method (default 'auto')
    orthogonalize : bool
        If True, use Lowdin orthogonalization (default False)
    zeta_A : float
        Orbital exponent scale for atom A (default 1.0)
    zeta_B : float
        Orbital exponent scale for atom B (default 1.0)

    Returns
    -------
    dict
        Keys: E_A_own, E_B_own, E_A_ghost, E_B_ghost, BSSE, BSSE_A, BSSE_B
    """
    import warnings

    # Own-basis energies (unscaled — isolated atoms at zeta=1.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        li_A = LatticeIndex(
            n_electrons=n_electrons_A, max_n=nmax_A, nuclear_charge=Z_A,
            vee_method=vee_method, h1_method='exact', fci_method=fci_method,
        )
        E_A_own = li_A.compute_ground_state(n_states=1)[0][0]

        li_B = LatticeIndex(
            n_electrons=n_electrons_B, max_n=nmax_B, nuclear_charge=Z_B,
            vee_method=vee_method, h1_method='exact', fci_method=fci_method,
        )
        E_B_own = li_B.compute_ground_state(n_states=1)[0][0]

    # Ghost-basis energies: atom A with ghost B orbitals
    # Use same zeta_A, zeta_B as supermolecule for consistent basis space
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol_A_ghost = MolecularLatticeIndex(
            Z_A=Z_A, Z_B=0, nmax_A=nmax_A, nmax_B=nmax_B,
            R=R, n_electrons=n_electrons_A,
            vee_method=vee_method, fci_method=fci_method,
            orthogonalize=orthogonalize,
            zeta_A=zeta_A, zeta_B=zeta_B,
        )
        E_A_ghost = mol_A_ghost.compute_ground_state(n_states=1)[0][0]

    # Ghost-basis energies: atom B with ghost A orbitals
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol_B_ghost = MolecularLatticeIndex(
            Z_A=0, Z_B=Z_B, nmax_A=nmax_A, nmax_B=nmax_B,
            R=R, n_electrons=n_electrons_B,
            vee_method=vee_method, fci_method=fci_method,
            orthogonalize=orthogonalize,
            zeta_A=zeta_A, zeta_B=zeta_B,
        )
        E_B_ghost = mol_B_ghost.compute_ground_state(n_states=1)[0][0]

    BSSE_A = E_A_ghost - E_A_own
    BSSE_B = E_B_ghost - E_B_own
    BSSE = BSSE_A + BSSE_B

    return {
        'E_A_own': E_A_own,
        'E_B_own': E_B_own,
        'E_A_ghost': E_A_ghost,
        'E_B_ghost': E_B_ghost,
        'BSSE': BSSE,
        'BSSE_A': BSSE_A,
        'BSSE_B': BSSE_B,
    }


# Cache for compute_atomic_p0 results
_atomic_p0_cache: Dict[Tuple[int, int], float] = {}


def compute_atomic_p0(Z: int, nmax: int) -> float:
    """
    Compute the self-consistent Sturmian p0 for an isolated atom.

    For a single atom with nuclear charge Z and basis truncation nmax,
    the self-consistent p0 satisfies p0 = sqrt(-2 * E_atom), where
    E_atom is the ground-state eigenvalue of the atom at that p0.

    For hydrogen (Z=1), p0 = 1.0 exactly at any nmax >= 1, since
    E_H = -0.5 Ha and sqrt(1.0) = 1.0.

    For Z=0 (ghost atoms), returns 0.0.

    Results are cached by (Z, nmax) so each computation is done once.

    Parameters
    ----------
    Z : int
        Nuclear charge.
    nmax : int
        Maximum principal quantum number.

    Returns
    -------
    float
        Self-consistent p0 for the isolated atom.
    """
    if Z == 0:
        return 0.0

    key = (Z, nmax)
    if key in _atomic_p0_cache:
        return _atomic_p0_cache[key]

    # Run single-atom FCI to get ground-state energy
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        li = LatticeIndex(
            n_electrons=Z,
            max_n=nmax,
            nuclear_charge=Z,
            vee_method='slater_full',
            fci_method='auto',
        )
        eigvals, _ = li.compute_ground_state(n_states=1)
        E_atom = eigvals[0]

    if E_atom >= 0:
        # Atom is unbound at this nmax — fall back to Z/1
        p0 = float(Z)
        print(f"[compute_atomic_p0] Z={Z}, nmax={nmax}: E_atom={E_atom:.6f} "
              f"(positive, using p0={p0:.4f})")
    else:
        p0 = np.sqrt(-2.0 * E_atom)
        print(f"[compute_atomic_p0] Z={Z}, nmax={nmax}: E_atom={E_atom:.6f} "
              f"Ha, p0={p0:.6f}")

    _atomic_p0_cache[key] = p0
    return p0
