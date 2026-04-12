"""
Harmonic oscillator shell model on S³-style graph.

Implements the 3D isotropic harmonic oscillator shell structure, verifies
harmonic magic numbers, and compares angular/radial sparsity to the Coulomb
case.  All angular integrals reuse the GeoVac Wigner 3j infrastructure.

References:
  - Mayer & Jensen, "Elementary Theory of Nuclear Shell Structure" (1955)
  - GeoVac Paper 0 (packing construction), Paper 7 (S³ proof)
"""

from __future__ import annotations

from collections import defaultdict
from math import factorial, sqrt as msqrt
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import sparse
from scipy.special import eval_genlaguerre, gamma as gamma_fn

from geovac.angular_integrals import wigner3j


# ── helpers ─────────────────────────────────────────────────────────────────

def _ho_shell_degeneracy(N: int) -> int:
    """Degeneracy of HO shell N (including spin): (N+1)(N+2)."""
    return (N + 1) * (N + 2)


def _l_values_in_shell(N: int) -> List[int]:
    """Allowed l values in HO shell N: l = N, N-2, N-4, ... >= 0."""
    return list(range(N, -1, -2))


# ── 1. Exact HO energies ───────────────────────────────────────────────────

def harmonic_oscillator_energies(
    n_max: int,
    l_max: Optional[int] = None,
    hw: float = 1.0,
) -> Dict[Tuple[int, int, int, int], float]:
    """
    Return exact 3D isotropic HO single-particle energies.

    Parameters
    ----------
    n_max : int
        Number of HO shells (N = 0 .. n_max-1).
    l_max : int or None
        If given, restrict orbital angular momentum.  Default: no restriction.
    hw : float
        Oscillator energy hbar*omega (Ha).

    Returns
    -------
    dict  (n_r, l, m_l, sigma) -> E   where E = hw*(2*n_r + l + 3/2).
    """
    energies: Dict[Tuple[int, int, int, int], float] = {}
    for N in range(n_max):
        E_N = hw * (N + 1.5)
        for l in _l_values_in_shell(N):
            if l_max is not None and l > l_max:
                continue
            n_r = (N - l) // 2
            for m_l in range(-l, l + 1):
                for sigma in (0, 1):
                    energies[(n_r, l, m_l, sigma)] = E_N
    return energies


# ── 2. Nuclear graph with HO node weights ──────────────────────────────────

def build_nuclear_graph_harmonic(
    n_max: int,
    hw: float = 1.0,
) -> Dict:
    """
    Construct an S³-style graph with harmonic oscillator node weights.

    Parameters
    ----------
    n_max : int
        Number of HO shells (N = 0 .. n_max-1).
    hw : float
        Oscillator energy (Ha).

    Returns
    -------
    dict with keys:
        'energies'        – 1-D array of single-particle energies
        'states'          – list of (N, n_r, l, m_l, sigma) tuples
        'n_states'        – total number of single-particle states
        'shell_closures'  – list of cumulative particle counts at each N
        'hamiltonian'     – scipy.sparse.csr_matrix (diagonal node weights
                            + off-diagonal transitions matching GeometricLattice
                            pattern: Delta-N = +-1 OR Delta-l = +-1)
    """
    states: List[Tuple[int, int, int, int, int]] = []
    energies_list: List[float] = []

    for N in range(n_max):
        E_N = hw * (N + 1.5)
        for l in _l_values_in_shell(N):
            n_r = (N - l) // 2
            for m_l in range(-l, l + 1):
                for sigma in (0, 1):
                    states.append((N, n_r, l, m_l, sigma))
                    energies_list.append(E_N)

    n_states = len(states)
    energies_arr = np.array(energies_list)

    # Build adjacency: connect states sharing spin with Delta-N = +-1
    # (radial transition) OR Delta-l = +-1 with Delta-m_l in {-1,0,1}
    # (angular transition), same shell allowed.
    rows, cols, vals = [], [], []
    for i in range(n_states):
        Ni, nri, li, mi, si = states[i]
        for j in range(i + 1, n_states):
            Nj, nrj, lj, mj, sj = states[j]
            if si != sj:
                continue  # same spin block
            connected = False
            # Delta-N = +-1 transitions (same l, m)
            if abs(Ni - Nj) == 1 and li == lj and mi == mj:
                connected = True
            # Delta-l = +-1 transitions (same N, Delta-m in {-1,0,1})
            if Ni == Nj and abs(li - lj) == 2 and abs(mi - mj) <= 1:
                # l changes by 2 within a shell (since l parity is fixed within
                # a shell, Delta-l=2 is the minimal intra-shell angular change)
                connected = True
            if connected:
                rows.extend([i, j])
                cols.extend([j, i])
                vals.extend([1.0, 1.0])

    # Degree matrix – adjacency
    A = sparse.csr_matrix((vals, (rows, cols)), shape=(n_states, n_states))
    degree = np.array(A.sum(axis=1)).ravel()
    D = sparse.diags(degree)
    L = D - A  # graph Laplacian

    # Hamiltonian = graph Laplacian scaled to HO + diagonal node weights
    # Use the GeoVac pattern: H = kinetic_scale * L + diag(energies)
    # For the HO the kinetic scale is not kappa=-1/16 (that is Coulomb-specific).
    # We use a nominal kinetic_scale that gives small off-diagonal coupling;
    # the physics is in the node weights (energies).
    kinetic_scale = hw * 0.1  # small coupling relative to shell spacing
    H = kinetic_scale * L + sparse.diags(energies_arr)
    H = sparse.csr_matrix(H)

    # Shell closures
    shell_closures: List[int] = []
    cumulative = 0
    for N in range(n_max):
        cumulative += _ho_shell_degeneracy(N)
        shell_closures.append(cumulative)

    return {
        "energies": energies_arr,
        "states": states,
        "n_states": n_states,
        "shell_closures": shell_closures,
        "hamiltonian": H,
    }


# ── 3. Verify magic numbers ────────────────────────────────────────────────

def verify_magic_numbers(n_max: int = 6) -> Dict:
    """
    Compute HO shell closures and compare to known magic numbers.

    Returns
    -------
    dict with keys 'expected', 'computed', 'match', 'details'.
    """
    expected = [2, 8, 20, 40, 70, 112][:n_max]

    details: List[Dict] = []
    cumulative = 0
    computed: List[int] = []
    for N in range(n_max):
        deg = _ho_shell_degeneracy(N)
        cumulative += deg
        computed.append(cumulative)
        details.append({
            "N": N,
            "l_values": _l_values_in_shell(N),
            "degeneracy": deg,
            "cumulative": cumulative,
        })

    return {
        "expected": expected,
        "computed": computed,
        "match": computed == expected,
        "details": details,
    }


# ── 4. Harmonic radial wavefunctions ───────────────────────────────────────

def harmonic_radial_wavefunctions(
    n_r: int,
    l: int,
    r_grid: np.ndarray,
    b: float = 1.0,
) -> np.ndarray:
    """
    Compute the HO radial wavefunction R_{n_r,l}(r).

    R_{n_r,l}(r) = N_nl * (r/b)^l * L_{n_r}^{l+1/2}((r/b)^2)
                   * exp(-r^2 / (2 b^2))

    Normalization: integral_0^infty |R|^2 r^2 dr = 1.

    Parameters
    ----------
    n_r : int   Radial quantum number (>= 0).
    l   : int   Orbital angular momentum (>= 0).
    r_grid : ndarray   Radial grid points.
    b   : float  Oscillator length sqrt(hbar/(m*omega)).

    Returns
    -------
    ndarray  R_{n_r,l}(r) on the grid.
    """
    x = (r_grid / b) ** 2  # dimensionless variable
    # Associated Laguerre polynomial L_{n_r}^{l+1/2}(x)
    L_poly = eval_genlaguerre(n_r, l + 0.5, x)

    # Unnormalized wavefunction
    wf = (r_grid / b) ** l * L_poly * np.exp(-x / 2.0)

    # Normalization constant from analytical formula:
    # N^2 = 2 * n_r! / (b^3 * Gamma(n_r + l + 3/2))
    norm_sq_analytic = 2.0 * factorial(n_r) / (b ** 3 * gamma_fn(n_r + l + 1.5))
    N_const = np.sqrt(norm_sq_analytic)

    return N_const * wf


# ── 5. Harmonic Slater integrals ───────────────────────────────────────────

def harmonic_slater_integrals(
    n_max: int,
    b: float = 1.0,
    n_grid: int = 2000,
) -> Dict[Tuple[int, ...], float]:
    """
    Compute R^k Slater radial integrals using HO radial wavefunctions.

    Uses the cumulative-charge (Y^k) method from composed_qubit.py.

    Parameters
    ----------
    n_max : int
        Number of HO shells (N = 0 .. n_max-1).
    b : float
        Oscillator length.
    n_grid : int
        Number of radial grid points.

    Returns
    -------
    dict mapping (n_r1,l1, n_r2,l2, n_r3,l3, n_r4,l4, k) -> float
    """
    # Enumerate unique (n_r, l) pairs
    unique_nl: List[Tuple[int, int]] = []
    for N in range(n_max):
        for l in _l_values_in_shell(N):
            n_r = (N - l) // 2
            if (n_r, l) not in unique_nl:
                unique_nl.append((n_r, l))
    unique_nl.sort()

    # Radial grid
    r_max = 6.0 * b * np.sqrt(n_max + 2)  # extend well beyond classical turning point
    r_grid = np.linspace(0, r_max, n_grid + 1)[1:]  # exclude r=0
    dr = r_grid[1] - r_grid[0]

    # Pre-compute radial wavefunctions
    R_on_grid: Dict[Tuple[int, int], np.ndarray] = {}
    for n_r, l in unique_nl:
        R_on_grid[(n_r, l)] = harmonic_radial_wavefunctions(n_r, l, r_grid, b)

    # Determine needed integrals (angular selection rules)
    needed: List[Tuple[int, ...]] = []
    for n1, l1 in unique_nl:
        for n2, l2 in unique_nl:
            for n3, l3 in unique_nl:
                for n4, l4 in unique_nl:
                    k_max = min(l1 + l3, l2 + l4)
                    for k in range(0, k_max + 1):
                        if (l1 + l3 + k) % 2 != 0:
                            continue
                        if (l2 + l4 + k) % 2 != 0:
                            continue
                        needed.append((n1, l1, n2, l2, n3, l3, n4, l4, k))

    if not needed:
        return {}

    # Group by (n2,l2,n4,l4,k) for Y^k potential reuse
    yk_groups: Dict[Tuple[int, ...], List[Tuple[int, ...]]] = defaultdict(list)
    for key in needed:
        n1, l1, n2, l2, n3, l3, n4, l4, k = key
        yk_key = (n2, l2, n4, l4, k)
        yk_groups[yk_key].append(key)

    # Compute Y^k potentials (cumulative-charge method)
    yk_cache: Dict[Tuple[int, ...], np.ndarray] = {}
    for yk_key in yk_groups:
        n2, l2, n4, l4, k = yk_key
        f_r = R_on_grid[(n2, l2)] * R_on_grid[(n4, l4)] * r_grid ** 2
        yk = np.zeros(n_grid)
        for i in range(n_grid):
            r1 = r_grid[i]
            inner = np.sum(f_r[: i + 1] * (r_grid[: i + 1] / r1) ** k) * dr
            if i + 1 < n_grid:
                outer = (
                    np.sum(f_r[i + 1 :] * (r1 / r_grid[i + 1 :]) ** k / r_grid[i + 1 :])
                    * dr
                )
            else:
                outer = 0.0
            yk[i] = inner / r1 + outer
        yk_cache[yk_key] = yk

    # Compute R^k integrals
    rk_cache: Dict[Tuple[int, ...], float] = {}
    for yk_key, keys in yk_groups.items():
        yk = yk_cache[yk_key]
        for key in keys:
            n1, l1, _, _, n3, l3, _, _, _ = key
            integrand = R_on_grid[(n1, l1)] * R_on_grid[(n3, l3)] * yk * r_grid ** 2
            val = np.trapezoid(integrand, r_grid)
            rk_cache[key] = val

    return rk_cache


# ── 6. Sparsity analysis ───────────────────────────────────────────────────

def _ck_coefficient(la: int, ma: int, lc: int, mc: int, k: int) -> float:
    """
    Gaunt angular coupling: c^k(l,m; l',m').

    c^k = (-1)^m sqrt((2l+1)(2l'+1)) * (l k l'; 0 0 0) * (l k l'; -m q m')
    """
    q = mc - ma
    pre = (-1) ** ma * msqrt((2 * la + 1) * (2 * lc + 1))
    w1 = wigner3j(la, k, lc, 0, 0, 0)
    if abs(w1) < 1e-15:
        return 0.0
    w2 = wigner3j(la, k, lc, -ma, q, mc)
    return pre * w1 * w2


def _build_eri_block_ho(
    states_nlm: List[Tuple[int, int, int]],
    rk_cache: Dict[Tuple[int, ...], float],
    threshold: float = 1e-15,
) -> Dict[Tuple[int, int, int, int], float]:
    """
    Build ERI dict <ab|cd> for HO states using Gaunt selection rules.

    Parameters
    ----------
    states_nlm : list of (n_r, l, m_l)
        Spatial orbital quantum numbers.
    rk_cache : dict
        R^k integrals from harmonic_slater_integrals.
    threshold : float
        Threshold for non-zero ERI.

    Returns
    -------
    dict  (a, b, c, d) -> float
    """
    n_sp = len(states_nlm)
    eri: Dict[Tuple[int, int, int, int], float] = {}

    # Pre-compute c^k table
    ck_table: Dict[Tuple[int, int, int], float] = {}
    for a in range(n_sp):
        la, ma = states_nlm[a][1], states_nlm[a][2]
        for c in range(n_sp):
            lc, mc = states_nlm[c][1], states_nlm[c][2]
            k_max = la + lc
            for k in range(0, k_max + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                val = _ck_coefficient(la, ma, lc, mc, k)
                if abs(val) > 1e-15:
                    ck_table[(a, c, k)] = val

    # Group by (a, c) pair
    ac_k_map: Dict[Tuple[int, int], List[Tuple[int, float]]] = defaultdict(list)
    for (a, c, k), val in ck_table.items():
        ac_k_map[(a, c)].append((k, val))

    for (a, c), ck_ac_list in ac_k_map.items():
        na, la, ma = states_nlm[a]
        nc, lc, mc = states_nlm[c]
        for (b, d), ck_bd_list in ac_k_map.items():
            nb, lb, mb = states_nlm[b]
            nd, ld, md = states_nlm[d]
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
            if abs(val) > threshold:
                eri[(a, b, c, d)] = val

    return eri


def _coulomb_rk_integrals(
    states_nlm: List[Tuple[int, int, int]],
    Z: float = 1.0,
    n_grid: int = 2000,
) -> Dict[Tuple[int, ...], float]:
    """
    Compute Coulomb R^k Slater integrals for comparison.

    Uses hydrogenic radial wavefunctions with the same cumulative-charge method.
    """
    from scipy.special import genlaguerre as _genlaguerre

    unique_nl = sorted(set((n_r, l) for n_r, l, m in states_nlm))

    # Map HO (n_r, l) to hydrogenic n = n_r + l + 1
    def _hydro_n(n_r: int, l: int) -> int:
        return n_r + l + 1

    r_max = 80.0 / max(Z, 0.5)
    r_grid = np.linspace(0, r_max, n_grid + 1)[1:]
    dr = r_grid[1] - r_grid[0]

    # Hydrogenic radial wavefunctions
    R_on_grid: Dict[Tuple[int, int], np.ndarray] = {}
    for n_r, l in unique_nl:
        n = _hydro_n(n_r, l)
        rho = 2.0 * Z * r_grid / n
        L_poly = _genlaguerre(n - l - 1, 2 * l + 1)(rho)
        wf = rho ** l * np.exp(-rho / 2.0) * L_poly
        norm_sq = np.trapezoid(wf ** 2 * r_grid ** 2, r_grid)
        if norm_sq < 1e-30:
            R_on_grid[(n_r, l)] = np.zeros_like(r_grid)
        else:
            R_on_grid[(n_r, l)] = wf / np.sqrt(norm_sq)

    # Determine needed integrals
    needed: List[Tuple[int, ...]] = []
    for n1, l1 in unique_nl:
        for n2, l2 in unique_nl:
            for n3, l3 in unique_nl:
                for n4, l4 in unique_nl:
                    k_max = min(l1 + l3, l2 + l4)
                    for k in range(0, k_max + 1):
                        if (l1 + l3 + k) % 2 != 0:
                            continue
                        if (l2 + l4 + k) % 2 != 0:
                            continue
                        needed.append((n1, l1, n2, l2, n3, l3, n4, l4, k))

    if not needed:
        return {}

    yk_groups: Dict[Tuple[int, ...], List[Tuple[int, ...]]] = defaultdict(list)
    for key in needed:
        n1, l1, n2, l2, n3, l3, n4, l4, k = key
        yk_key = (n2, l2, n4, l4, k)
        yk_groups[yk_key].append(key)

    yk_cache: Dict[Tuple[int, ...], np.ndarray] = {}
    for yk_key in yk_groups:
        n2, l2, n4, l4, k = yk_key
        f_r = R_on_grid[(n2, l2)] * R_on_grid[(n4, l4)] * r_grid ** 2
        yk = np.zeros(n_grid)
        for i in range(n_grid):
            r1 = r_grid[i]
            inner = np.sum(f_r[: i + 1] * (r_grid[: i + 1] / r1) ** k) * dr
            if i + 1 < n_grid:
                outer = (
                    np.sum(f_r[i + 1 :] * (r1 / r_grid[i + 1 :]) ** k / r_grid[i + 1 :])
                    * dr
                )
            else:
                outer = 0.0
            yk[i] = inner / r1 + outer
        yk_cache[yk_key] = yk

    rk_cache: Dict[Tuple[int, ...], float] = {}
    for yk_key, keys in yk_groups.items():
        yk = yk_cache[yk_key]
        for key in keys:
            n1, l1, _, _, n3, l3, _, _, _ = key
            integrand = R_on_grid[(n1, l1)] * R_on_grid[(n3, l3)] * yk * r_grid ** 2
            val = np.trapezoid(integrand, r_grid)
            rk_cache[key] = val

    return rk_cache


def sparsity_analysis(
    n_max: int,
    b: float = 1.0,
    n_grid: int = 2000,
    threshold: float = 1e-10,
) -> Dict:
    """
    Analyse ERI sparsity for the harmonic oscillator vs Coulomb potential.

    Parameters
    ----------
    n_max : int
        Number of HO shells.
    b : float
        Oscillator length for HO wavefunctions.
    n_grid : int
        Grid points for radial integration.
    threshold : float
        Threshold below which an ERI is counted as zero.

    Returns
    -------
    dict with sparsity metrics for HO and Coulomb.
    """
    # Build spatial state list (n_r, l, m_l)
    states_nlm: List[Tuple[int, int, int]] = []
    for N in range(n_max):
        for l in _l_values_in_shell(N):
            n_r = (N - l) // 2
            for m_l in range(-l, l + 1):
                states_nlm.append((n_r, l, m_l))

    n_sp = len(states_nlm)
    total_possible = n_sp ** 4

    # ── Angular sparsity (identical for both potentials) ──────────────
    # Count pairs (a, c) with at least one non-zero c^k
    angular_nonzero_pairs = 0
    angular_total_pairs = n_sp * n_sp
    ck_table_ho: Dict[Tuple[int, int, int], float] = {}
    for a in range(n_sp):
        la, ma = states_nlm[a][1], states_nlm[a][2]
        for c in range(n_sp):
            lc, mc = states_nlm[c][1], states_nlm[c][2]
            found = False
            k_max = la + lc
            for k in range(0, k_max + 1):
                if (la + lc + k) % 2 != 0:
                    continue
                val = _ck_coefficient(la, ma, lc, mc, k)
                if abs(val) > 1e-15:
                    ck_table_ho[(a, c, k)] = val
                    found = True
            if found:
                angular_nonzero_pairs += 1

    angular_zero_fraction = 1.0 - angular_nonzero_pairs / angular_total_pairs

    # ── HO radial + full ERI ──────────────────────────────────────────
    rk_ho = harmonic_slater_integrals(n_max, b=b, n_grid=n_grid)
    eri_ho = _build_eri_block_ho(states_nlm, rk_ho, threshold=threshold)
    ho_nonzero = len(eri_ho)

    # ── Coulomb radial + full ERI ─────────────────────────────────────
    rk_coulomb = _coulomb_rk_integrals(states_nlm, Z=1.0, n_grid=n_grid)
    eri_coulomb = _build_eri_block_ho(states_nlm, rk_coulomb, threshold=threshold)
    coulomb_nonzero = len(eri_coulomb)

    # ── Count angular-only zeros (same for both) ──────────────────────
    # An ERI <ab|cd> is angular-zero if there is no k with both
    # c^k(a,c) != 0 and c^k(b,d) != 0
    ac_k_map: Dict[Tuple[int, int], List[int]] = defaultdict(list)
    for (a, c, k) in ck_table_ho:
        ac_k_map[(a, c)].append(k)

    angular_allowed = 0
    for (a, c), ks_ac in ac_k_map.items():
        la, ma = states_nlm[a][1], states_nlm[a][2]
        lc, mc = states_nlm[c][1], states_nlm[c][2]
        for (b, d), ks_bd in ac_k_map.items():
            lb, mb = states_nlm[b][1], states_nlm[b][2]
            ld, md = states_nlm[d][1], states_nlm[d][2]
            if ma + mb != mc + md:
                continue
            # Check if any common k exists
            if set(ks_ac) & set(ks_bd):
                angular_allowed += 1

    angular_zeros = total_possible - angular_allowed

    # Radial zeros: angular-allowed but integral < threshold
    ho_radial_zeros = angular_allowed - ho_nonzero
    coulomb_radial_zeros = angular_allowed - coulomb_nonzero

    return {
        "n_max": n_max,
        "n_spatial_orbitals": n_sp,
        "total_possible_eris": total_possible,
        "angular_zeros": angular_zeros,
        "angular_zero_fraction": angular_zeros / max(total_possible, 1),
        "angular_allowed": angular_allowed,
        "ho": {
            "nonzero_eris": ho_nonzero,
            "radial_zeros": ho_radial_zeros,
            "total_sparsity": 1.0 - ho_nonzero / max(total_possible, 1),
            "eri_density": ho_nonzero / max(total_possible, 1),
        },
        "coulomb": {
            "nonzero_eris": coulomb_nonzero,
            "radial_zeros": coulomb_radial_zeros,
            "total_sparsity": 1.0 - coulomb_nonzero / max(total_possible, 1),
            "eri_density": coulomb_nonzero / max(total_possible, 1),
        },
    }
