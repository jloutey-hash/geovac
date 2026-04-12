"""
Sparsity Characterization Across Radial Potentials (Track NC)
=============================================================

Measures how ERI sparsity depends on the choice of radial potential.

The key factorization:
    <n1l1m1 n2l2m2|V|n3l3m3 n4l4m4> = sum_k R^k(radial) * C^k(angular)

Angular selection rules (Gaunt/Wigner 3j) depend ONLY on (l, m, k) and are
potential-INDEPENDENT.  Radial integrals R^k depend on the wavefunctions,
which change with the potential.

Hypothesis: angular selection rules produce most of the sparsity, so the
O(Q^2.5) Pauli scaling is largely potential-independent.

All computations in atomic units (hbar = m_e = e = 1) unless noted.
Nuclear potentials (Woods-Saxon) are converted to atomic units.

Author: GeoVac Development Team
Date: April 2026
"""

from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
from scipy.linalg import eigh_tridiagonal

from geovac.angular_integrals import wigner3j


# ---------------------------------------------------------------------------
# Angular coefficient (Gaunt c^k) — reuses the same formula as composed_qubit
# ---------------------------------------------------------------------------

def ck_coefficient(la: int, ma: int, lc: int, mc: int, k: int) -> float:
    """
    c^k(l,m; l',m') = (-1)^m sqrt((2l+1)(2l'+1)) * (l k l'; 0 0 0) * (l k l'; -m q m')
    where q = mc - ma.

    Potential-INDEPENDENT angular coupling coefficient.
    """
    q = mc - ma
    pre = ((-1) ** ma) * np.sqrt((2 * la + 1) * (2 * lc + 1))
    w1 = wigner3j(la, k, lc, 0, 0, 0)
    if abs(w1) < 1e-15:
        return 0.0
    w2 = wigner3j(la, k, lc, -ma, q, mc)
    return pre * w1 * w2


# ---------------------------------------------------------------------------
# Radial Schrodinger solver — Numerov method
# ---------------------------------------------------------------------------

def solve_radial_schrodinger(
    V_func: Callable[[np.ndarray], np.ndarray],
    n_r: int,
    l: int,
    r_max: float = 50.0,
    n_grid: int = 2000,
) -> Tuple[float, np.ndarray, np.ndarray]:
    """
    Solve the 1D radial Schrodinger equation for the n_r-th bound state
    with angular momentum l.

    Equation: [-1/2 d^2u/dr^2 + V_eff(r) u] = E u
    where u(r) = r R(r) and V_eff(r) = V(r) + l(l+1)/(2r^2).

    Uses a grid eigensolver (finite-difference Hamiltonian).

    Parameters
    ----------
    V_func : callable
        V(r) — the radial potential as a function of r array.
    n_r : int
        Radial node count (0 = ground, 1 = first excited, ...).
    l : int
        Angular momentum quantum number.
    r_max : float
        Outer boundary of the radial grid.
    n_grid : int
        Number of interior grid points.

    Returns
    -------
    (E, R_grid, r_grid) : (float, ndarray, ndarray)
        Energy, radial wavefunction R(r), and grid points.
    """
    # Build uniform grid (interior points only; u(0)=u(r_max)=0 boundary)
    dr = r_max / (n_grid + 1)
    r = np.linspace(dr, r_max - dr, n_grid)

    # Effective potential
    V_r = V_func(r)
    V_eff = V_r + 0.5 * l * (l + 1) / r**2

    # Finite-difference kinetic: -1/2 d^2/dr^2 on uniform grid
    # T_{ii} = 1/dr^2, T_{i,i+1} = T_{i+1,i} = -1/(2 dr^2)
    diag = V_eff + 1.0 / dr**2
    off_diag = np.full(n_grid - 1, -0.5 / dr**2)

    # Solve for the lowest few eigenvalues
    n_eig = min(n_r + 5, n_grid)
    # Use full tridiagonal solver for small grids, subset for large
    if n_grid <= 500:
        evals, evecs = eigh_tridiagonal(diag, off_diag)
    else:
        # Use sparse solver for large grids
        from scipy.sparse import diags
        from scipy.sparse.linalg import eigsh
        H_sparse = diags(
            [off_diag, diag, off_diag], [-1, 0, 1],
            shape=(n_grid, n_grid), format='csr',
        )
        n_eig = min(n_r + 10, n_grid - 2)
        evals, evecs = eigsh(H_sparse, k=n_eig, which='SA')
        idx = np.argsort(evals)
        evals = evals[idx]
        evecs = evecs[:, idx]

    # Pick the n_r-th state.
    # For confining potentials (harmonic, square well, Woods-Saxon), all
    # eigenvalues are "bound" even if E > 0.  For Coulomb/Yukawa, bound
    # states have E < 0 and continuum states E > 0.
    # Strategy: use E < 0 filter for potentials that support continuum,
    # but fall back to all eigenvalues if no negative eigenvalues exist
    # (i.e., confining potentials).
    bound_mask = evals < 0
    if np.any(bound_mask):
        bound_evals = evals[bound_mask]
        bound_evecs = evecs[:, bound_mask]
    else:
        # All eigenvalues positive — confining potential, all states are bound
        bound_evals = evals
        bound_evecs = evecs

    if n_r >= len(bound_evals):
        raise ValueError(
            f"Only {len(bound_evals)} bound/discrete states found for l={l}, "
            f"requested n_r={n_r}"
        )

    E = bound_evals[n_r]
    u = bound_evecs[:, n_r]

    # Convert u(r) = r*R(r) to R(r)
    R_grid = u / r

    # Normalize: integral |R(r)|^2 r^2 dr = 1
    norm_sq = np.trapezoid(R_grid**2 * r**2, r)
    if norm_sq > 1e-30:
        R_grid /= np.sqrt(norm_sq)

    # Fix phase convention: first maximum should be positive
    abs_R = np.abs(R_grid)
    first_max_idx = np.argmax(abs_R[:len(abs_R)//2 + 1])
    if first_max_idx < len(R_grid) and R_grid[first_max_idx] < 0:
        R_grid = -R_grid

    return E, R_grid, r


# ---------------------------------------------------------------------------
# Potential definitions
# ---------------------------------------------------------------------------

def _coulomb_potential(Z: float) -> Callable[[np.ndarray], np.ndarray]:
    """V(r) = -Z/r"""
    def V(r: np.ndarray) -> np.ndarray:
        return -Z / r
    return V


def _harmonic_potential(omega: float) -> Callable[[np.ndarray], np.ndarray]:
    """V(r) = 0.5 * omega^2 * r^2  (mass=1 in atomic units)"""
    def V(r: np.ndarray) -> np.ndarray:
        return 0.5 * omega**2 * r**2
    return V


def _woods_saxon_potential(
    V0: float, R0: float, a: float,
) -> Callable[[np.ndarray], np.ndarray]:
    """V(r) = -V0 / (1 + exp((r - R0) / a))"""
    def V(r: np.ndarray) -> np.ndarray:
        x = (r - R0) / a
        # Clip to avoid overflow
        x = np.clip(x, -500, 500)
        return -V0 / (1.0 + np.exp(x))
    return V


def _square_well_potential(
    V0: float, R: float,
) -> Callable[[np.ndarray], np.ndarray]:
    """V(r) = -V0 for r < R, 0 for r >= R"""
    def V(r: np.ndarray) -> np.ndarray:
        return np.where(r < R, -V0, 0.0)
    return V


def _yukawa_potential(
    V0: float, mu: float,
) -> Callable[[np.ndarray], np.ndarray]:
    """V(r) = -V0 * exp(-mu*r) / r"""
    def V(r: np.ndarray) -> np.ndarray:
        return -V0 * np.exp(-mu * r) / r
    return V


def make_potential(
    name: str, params: Dict,
) -> Tuple[Callable[[np.ndarray], np.ndarray], float, int]:
    """
    Create a potential function, suitable r_max, and n_grid.

    Returns (V_func, r_max, n_grid).
    """
    if name == 'coulomb':
        Z = params.get('Z', 1.0)
        return _coulomb_potential(Z), 60.0 / Z, 2000
    elif name == 'harmonic':
        omega = params.get('omega', 1.0)
        if 'hw' in params:
            omega = params['hw']  # hbar*omega with hbar=1
        # r_max ~ 3 * classical turning radius at moderate energy
        r_max = 3.0 * np.sqrt(2.0 * 20.0) / omega
        return _harmonic_potential(omega), min(r_max, 30.0), 2000
    elif name == 'woods_saxon':
        V0 = params['V0']
        R0 = params['R0']
        a = params['a']
        r_max = R0 + 10.0 * a
        return _woods_saxon_potential(V0, R0, a), r_max, 2000
    elif name == 'square_well':
        V0 = params['V0']
        R = params['R']
        return _square_well_potential(V0, R), R * 4.0, 2000
    elif name == 'yukawa':
        V0 = params.get('V0', 1.0)
        mu = params.get('mu', 0.5)
        r_max = min(60.0, 10.0 / mu)
        return _yukawa_potential(V0, mu), r_max, 2000
    else:
        raise ValueError(f"Unknown potential: {name}")


# ---------------------------------------------------------------------------
# Radial wavefunctions for a given potential
# ---------------------------------------------------------------------------

def radial_wavefunctions_for_potential(
    potential_name: str,
    params: Dict,
    n_max: int,
    l_max: int,
) -> Dict[Tuple[int, int], Tuple[float, np.ndarray, np.ndarray]]:
    """
    Compute all radial wavefunctions R_{n_r, l}(r) for a given potential
    up to n_max radial nodes and l_max angular momentum.

    Parameters
    ----------
    potential_name : str
        One of 'coulomb', 'harmonic', 'woods_saxon', 'square_well', 'yukawa'.
    params : dict
        Potential-specific parameters.
    n_max : int
        Maximum number of radial nodes (inclusive).
    l_max : int
        Maximum angular momentum (inclusive).

    Returns
    -------
    dict mapping (n_r, l) -> (E, R_grid, r_grid)
    """
    V_func, r_max, n_grid = make_potential(potential_name, params)

    wavefunctions: Dict[Tuple[int, int], Tuple[float, np.ndarray, np.ndarray]] = {}

    for l_val in range(l_max + 1):
        for n_r in range(n_max + 1):
            try:
                E, R_grid, r_grid = solve_radial_schrodinger(
                    V_func, n_r, l_val, r_max=r_max, n_grid=n_grid,
                )
                wavefunctions[(n_r, l_val)] = (E, R_grid, r_grid)
            except ValueError:
                # Not enough bound states for this (n_r, l) — skip
                pass

    return wavefunctions


# ---------------------------------------------------------------------------
# Slater R^k integrals — cumulative charge method
# ---------------------------------------------------------------------------

def compute_slater_rk(
    R_a: np.ndarray,
    R_b: np.ndarray,
    R_c: np.ndarray,
    R_d: np.ndarray,
    k: int,
    r_grid: np.ndarray,
) -> float:
    """
    Compute Slater R^k(ab; cd) integral via the cumulative-charge (Y^k) method.

    R^k(ab; cd) = int int R_a(r1) R_c(r1) [r_<^k / r_>^{k+1}]
                       R_b(r2) R_d(r2) r1^2 r2^2 dr1 dr2

    Parameters
    ----------
    R_a, R_b, R_c, R_d : ndarray
        Radial wavefunctions on the same grid.
    k : int
        Multipole order.
    r_grid : ndarray
        Radial grid points.

    Returns
    -------
    float
        Value of the Slater integral.
    """
    dr = r_grid[1] - r_grid[0]
    n = len(r_grid)

    # Build the Y^k potential from (b, d) pair using the same loop approach
    # as composed_qubit._compute_rk_integrals_block for accuracy.
    f_bd = R_b * R_d * r_grid**2
    yk = np.zeros(n)
    for i in range(n):
        r1 = r_grid[i]
        inner = np.sum(f_bd[:i + 1] * (r_grid[:i + 1] / r1) ** k) * dr
        if i + 1 < n:
            outer = np.sum(
                f_bd[i + 1:]
                * (r1 / r_grid[i + 1:]) ** k
                / r_grid[i + 1:]
            ) * dr
        else:
            outer = 0.0
        yk[i] = inner / r1 + outer

    # Now integrate R_a(r1) R_c(r1) Y^k(r1) r1^2 dr1
    integrand = R_a * R_c * yk * r_grid**2
    return float(np.trapezoid(integrand, r_grid))


# ---------------------------------------------------------------------------
# State enumeration
# ---------------------------------------------------------------------------

def enumerate_states(
    wavefunctions: Dict[Tuple[int, int], Tuple[float, np.ndarray, np.ndarray]],
) -> List[Tuple[int, int, int]]:
    """
    Enumerate spatial orbital states (n_r, l, m) from available wavefunctions.

    For each (n_r, l) with a computed wavefunction, generates states for
    m = -l, ..., +l.

    Returns sorted list of (n_r, l, m) tuples.
    """
    states = []
    for (n_r, l_val) in sorted(wavefunctions.keys()):
        for m_val in range(-l_val, l_val + 1):
            states.append((n_r, l_val, m_val))
    return states


# ---------------------------------------------------------------------------
# Full ERI tensor with angular + radial factorization
# ---------------------------------------------------------------------------

def compute_eri_tensor(
    wavefunctions: Dict[Tuple[int, int], Tuple[float, np.ndarray, np.ndarray]],
    threshold: float = 1e-10,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """
    Build the full ERI tensor using angular (Gaunt) + radial (Slater) factorization.

    ERI[a,b,c,d] = sum_k c^k(la,ma; lc,mc) * c^k(lb,mb; ld,md) * R^k(a,b; c,d)

    In PHYSICIST notation: <ab|cd>.

    Parameters
    ----------
    wavefunctions : dict
        Mapping (n_r, l) -> (E, R_grid, r_grid).
    threshold : float
        Zero threshold for ERI values.

    Returns
    -------
    (eri_tensor, stats) : (ndarray, dict)
        eri_tensor[a,b,c,d] in physicist notation.
        stats contains sparsity counts.
    """
    states = enumerate_states(wavefunctions)
    n_sp = len(states)

    # Ensure all wavefunctions share the same grid (interpolate if needed)
    # For simplicity, assume they do (same potential solver parameters)
    # Pick the r_grid from the first wavefunction
    first_key = next(iter(wavefunctions))
    _, _, r_grid_ref = wavefunctions[first_key]

    # Pre-compute c^k table
    ck_table: Dict[Tuple[int, int, int], float] = {}
    for a_idx in range(n_sp):
        la, ma = states[a_idx][1], states[a_idx][2]
        for c_idx in range(n_sp):
            lc, mc = states[c_idx][1], states[c_idx][2]
            k_max = la + lc
            for k in range(0, k_max + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                if abs(la - lc) > k:
                    continue
                val = ck_coefficient(la, ma, lc, mc, k)
                if abs(val) > 1e-15:
                    ck_table[(a_idx, c_idx, k)] = val

    # Group by (a, c) pair
    from collections import defaultdict
    ac_k_map: Dict[Tuple[int, int], List[Tuple[int, float]]] = defaultdict(list)
    for (a_idx, c_idx, k), val in ck_table.items():
        ac_k_map[(a_idx, c_idx)].append((k, val))

    # Pre-compute R^k integrals for all needed (n_r_a, l_a, n_r_b, l_b, ..., k)
    # Group needed R^k by (n_r_b, l_b, n_r_d, l_d, k) for Y^k caching
    rk_needed = set()
    for (a_idx, c_idx), ck_ac_list in ac_k_map.items():
        n_r_a, l_a, _ = states[a_idx]
        n_r_c, l_c, _ = states[c_idx]
        for (b_idx, d_idx), ck_bd_list in ac_k_map.items():
            n_r_b, l_b, _ = states[b_idx]
            n_r_d, l_d, _ = states[d_idx]
            m_a, m_b = states[a_idx][2], states[b_idx][2]
            m_c, m_d = states[c_idx][2], states[d_idx][2]
            if m_a + m_b != m_c + m_d:
                continue
            for k_ac, _ in ck_ac_list:
                for k_bd, _ in ck_bd_list:
                    if k_ac == k_bd:
                        rk_needed.add((n_r_a, l_a, n_r_b, l_b,
                                       n_r_c, l_c, n_r_d, l_d, k_ac))

    # Compute R^k integrals
    rk_cache: Dict[Tuple[int, ...], float] = {}
    # Group by (n_r_b, l_b, n_r_d, l_d, k) for Y^k caching
    yk_groups: Dict[Tuple[int, ...], List[Tuple[int, ...]]] = defaultdict(list)
    for key in rk_needed:
        n_r_a, l_a, n_r_b, l_b, n_r_c, l_c, n_r_d, l_d, k = key
        yk_key = (n_r_b, l_b, n_r_d, l_d, k)
        yk_groups[yk_key].append(key)

    for yk_key, keys in yk_groups.items():
        n_r_b, l_b, n_r_d, l_d, k = yk_key
        _, R_b, r_b = wavefunctions[(n_r_b, l_b)]
        _, R_d, r_d = wavefunctions[(n_r_d, l_d)]

        # Build Y^k potential (loop approach for accuracy)
        dr = r_b[1] - r_b[0]
        n_pts = len(r_b)
        f_bd = R_b * R_d * r_b**2
        yk = np.zeros(n_pts)
        for i in range(n_pts):
            r1 = r_b[i]
            inner = np.sum(f_bd[:i + 1] * (r_b[:i + 1] / r1) ** k) * dr
            if i + 1 < n_pts:
                outer = np.sum(
                    f_bd[i + 1:]
                    * (r1 / r_b[i + 1:]) ** k
                    / r_b[i + 1:]
                ) * dr
            else:
                outer = 0.0
            yk[i] = inner / r1 + outer

        for key in keys:
            n_r_a, l_a, _, _, n_r_c, l_c, _, _, _ = key
            _, R_a, r_a = wavefunctions[(n_r_a, l_a)]
            _, R_c, r_c = wavefunctions[(n_r_c, l_c)]
            integrand = R_a * R_c * yk * r_a**2
            rk_cache[key] = float(np.trapezoid(integrand, r_a))

    # Build ERI tensor
    eri = np.zeros((n_sp, n_sp, n_sp, n_sp))

    # Counting stats
    total_elements = 0
    angular_zero = 0
    radial_zero = 0
    combined_nonzero = 0

    for a_idx in range(n_sp):
        for b_idx in range(n_sp):
            for c_idx in range(n_sp):
                for d_idx in range(n_sp):
                    total_elements += 1

                    m_a = states[a_idx][2]
                    m_b = states[b_idx][2]
                    m_c = states[c_idx][2]
                    m_d = states[d_idx][2]

                    # m-conservation check
                    if m_a + m_b != m_c + m_d:
                        angular_zero += 1
                        continue

                    # Get angular (a,c) and (b,d) couplings
                    ck_ac_list = ac_k_map.get((a_idx, c_idx), [])
                    ck_bd_list = ac_k_map.get((b_idx, d_idx), [])

                    if not ck_ac_list or not ck_bd_list:
                        angular_zero += 1
                        continue

                    # Check if any k values match
                    k_ac_set = {kk for kk, _ in ck_ac_list}
                    k_bd_set = {kk for kk, _ in ck_bd_list}
                    if not k_ac_set.intersection(k_bd_set):
                        angular_zero += 1
                        continue

                    # Compute ERI
                    n_r_a, l_a, _ = states[a_idx]
                    n_r_b, l_b, _ = states[b_idx]
                    n_r_c, l_c, _ = states[c_idx]
                    n_r_d, l_d, _ = states[d_idx]

                    val = 0.0
                    for k_ac, c_ac in ck_ac_list:
                        for k_bd, c_bd in ck_bd_list:
                            if k_ac != k_bd:
                                continue
                            rk_key = (n_r_a, l_a, n_r_b, l_b,
                                      n_r_c, l_c, n_r_d, l_d, k_ac)
                            rk_val = rk_cache.get(rk_key, 0.0)
                            val += c_ac * c_bd * rk_val

                    if abs(val) < threshold:
                        radial_zero += 1
                    else:
                        eri[a_idx, b_idx, c_idx, d_idx] = val
                        combined_nonzero += 1

    stats = {
        'total_elements': total_elements,
        'angular_zero': angular_zero,
        'radial_zero': radial_zero,
        'combined_nonzero': combined_nonzero,
        'n_spatial_orbitals': n_sp,
        'n_qubits': 2 * n_sp,
        'angular_sparsity_pct': 100.0 * angular_zero / total_elements,
        'radial_sparsity_pct': 100.0 * radial_zero / total_elements,
        'total_sparsity_pct': 100.0 * (1.0 - combined_nonzero / total_elements),
        'eri_density_pct': 100.0 * combined_nonzero / total_elements,
    }

    return eri, stats


# ---------------------------------------------------------------------------
# Angular-only sparsity (potential-independent)
# ---------------------------------------------------------------------------

def angular_zero_count(n_max: int, l_max: int) -> Tuple[int, int]:
    """
    Count angular zeros purely from Gaunt selection rules.

    Enumerates all (n_r, l, m) states with n_r = 0..n_max, l = 0..l_max,
    m = -l..+l and counts how many ERI elements are zero by angular selection
    rules alone (triangle, parity, m-conservation).

    Returns (angular_zeros, total_elements).
    """
    states = []
    for n_r in range(n_max + 1):
        for l_val in range(l_max + 1):
            for m_val in range(-l_val, l_val + 1):
                states.append((n_r, l_val, m_val))

    n_sp = len(states)
    total = n_sp**4
    angular_zeros = 0

    # Pre-compute which (a, c) pairs have nonzero c^k for any k
    has_angular = {}
    for a_idx in range(n_sp):
        la, ma = states[a_idx][1], states[a_idx][2]
        for c_idx in range(n_sp):
            lc, mc = states[c_idx][1], states[c_idx][2]
            k_vals = []
            k_max = la + lc
            for k in range(0, k_max + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                if abs(la - lc) > k:
                    continue
                val = ck_coefficient(la, ma, lc, mc, k)
                if abs(val) > 1e-15:
                    k_vals.append(k)
            has_angular[(a_idx, c_idx)] = set(k_vals)

    for a_idx in range(n_sp):
        for b_idx in range(n_sp):
            for c_idx in range(n_sp):
                for d_idx in range(n_sp):
                    m_a = states[a_idx][2]
                    m_b = states[b_idx][2]
                    m_c = states[c_idx][2]
                    m_d = states[d_idx][2]

                    if m_a + m_b != m_c + m_d:
                        angular_zeros += 1
                        continue

                    k_ac = has_angular.get((a_idx, c_idx), set())
                    k_bd = has_angular.get((b_idx, d_idx), set())

                    if not k_ac or not k_bd or not k_ac.intersection(k_bd):
                        angular_zeros += 1

    return angular_zeros, total


# ---------------------------------------------------------------------------
# Sparsity table — main analysis
# ---------------------------------------------------------------------------

# Default potential parameters
DEFAULT_POTENTIALS = {
    'coulomb': {'Z': 1.0},
    'harmonic': {'omega': 1.0},
    'woods_saxon': {'V0': 50.0, 'R0': 3.0, 'a': 0.65},
    'square_well': {'V0': 50.0, 'R': 3.0},
    'yukawa': {'V0': 1.0, 'mu': 0.5},
}


def sparsity_table(
    n_max: int = 3,
    l_max: int = 2,
    potentials: Optional[Dict[str, Dict]] = None,
    threshold: float = 1e-10,
) -> Dict[str, Dict[str, float]]:
    """
    For each potential, compute radial wavefunctions, build the ERI tensor,
    and report sparsity statistics.

    Parameters
    ----------
    n_max : int
        Maximum radial node count.
    l_max : int
        Maximum angular momentum.
    potentials : dict, optional
        Mapping potential_name -> params. Defaults to DEFAULT_POTENTIALS.
    threshold : float
        Zero threshold for ERI values.

    Returns
    -------
    dict mapping potential_name -> stats dict
    """
    if potentials is None:
        potentials = DEFAULT_POTENTIALS

    results = {}

    for pot_name, params in potentials.items():
        try:
            wfs = radial_wavefunctions_for_potential(pot_name, params, n_max, l_max)
        except Exception as e:
            results[pot_name] = {'error': str(e)}
            continue

        if len(wfs) == 0:
            results[pot_name] = {'error': 'No bound states found'}
            continue

        _, stats = compute_eri_tensor(wfs, threshold=threshold)
        results[pot_name] = stats

    return results


def print_sparsity_table(results: Dict[str, Dict[str, float]]) -> str:
    """Format sparsity results as a readable table."""
    lines = []
    header = (
        f"{'Potential':<15} {'N_orb':>6} {'Q':>4} {'Total':>10} "
        f"{'Ang Zero%':>10} {'Rad Zero%':>10} {'Total Spar%':>12} {'ERI Dens%':>10}"
    )
    lines.append(header)
    lines.append('-' * len(header))

    for name, stats in results.items():
        if 'error' in stats:
            lines.append(f"{name:<15} ERROR: {stats['error']}")
            continue
        lines.append(
            f"{name:<15} {stats['n_spatial_orbitals']:>6} "
            f"{stats['n_qubits']:>4} "
            f"{stats['total_elements']:>10} "
            f"{stats['angular_sparsity_pct']:>10.2f} "
            f"{stats['radial_sparsity_pct']:>10.2f} "
            f"{stats['total_sparsity_pct']:>12.2f} "
            f"{stats['eri_density_pct']:>10.2f}"
        )

    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# Pauli scaling estimate
# ---------------------------------------------------------------------------

def pauli_scaling_estimate(
    potential_name: str,
    params: Dict,
    n_max_values: Optional[List[int]] = None,
    l_max: int = 2,
    threshold: float = 1e-10,
) -> Dict[str, List[float]]:
    """
    For each n_max, build second-quantized Hamiltonian (using OpenFermion)
    and extract Q, N_Pauli.

    Parameters
    ----------
    potential_name : str
        Potential name.
    params : dict
        Potential parameters.
    n_max_values : list of int, optional
        Values of n_max to scan. Defaults to [1, 2, 3].
    l_max : int
        Maximum angular momentum.
    threshold : float
        Zero threshold.

    Returns
    -------
    dict with keys 'n_max', 'Q', 'N_Pauli', 'scaling_exponent'.
    """
    if n_max_values is None:
        n_max_values = [1, 2, 3]

    from openfermion import jordan_wigner, count_qubits
    from geovac.qubit_encoding import build_fermion_op_from_integrals

    Q_list = []
    N_Pauli_list = []
    n_max_used = []

    for n_max_val in n_max_values:
        try:
            wfs = radial_wavefunctions_for_potential(
                potential_name, params, n_max_val, l_max,
            )
        except Exception:
            continue

        if len(wfs) < 2:
            continue

        states = enumerate_states(wfs)
        n_sp = len(states)

        # Build ERI tensor (physicist notation)
        eri_phys, _ = compute_eri_tensor(wfs, threshold=threshold)

        # Convert to chemist notation: (pq|rs) = <pr|qs>
        eri_chem = np.zeros((n_sp, n_sp, n_sp, n_sp))
        for a in range(n_sp):
            for b in range(n_sp):
                for c in range(n_sp):
                    for d in range(n_sp):
                        eri_chem[a, c, b, d] = eri_phys[a, b, c, d]

        # Build one-body integrals (diagonal energies)
        h1 = np.zeros((n_sp, n_sp))
        for i, (n_r, l_val, m_val) in enumerate(states):
            E_i, _, _ = wfs[(n_r, l_val)]
            h1[i, i] = E_i

        # Build fermion operator and JW transform
        fermion_op = build_fermion_op_from_integrals(h1, eri_chem, 0.0)
        qubit_op = jordan_wigner(fermion_op)

        Q = 2 * n_sp
        n_pauli = len(qubit_op.terms) - (1 if () in qubit_op.terms else 0)

        Q_list.append(Q)
        N_Pauli_list.append(n_pauli)
        n_max_used.append(n_max_val)

    # Fit scaling exponent
    scaling_exponent = None
    if len(Q_list) >= 2:
        log_Q = np.log(np.array(Q_list, dtype=float))
        log_N = np.log(np.array(N_Pauli_list, dtype=float))
        if log_Q[-1] != log_Q[0]:
            coeffs = np.polyfit(log_Q, log_N, 1)
            scaling_exponent = coeffs[0]

    return {
        'n_max': n_max_used,
        'Q': Q_list,
        'N_Pauli': N_Pauli_list,
        'scaling_exponent': scaling_exponent,
    }
