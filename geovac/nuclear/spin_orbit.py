"""
Nuclear shell model with spin-orbit coupling.

Adds the Mayer-Jensen spin-orbit term H_SO = -V_ls * (l . s) to the
harmonic oscillator Hamiltonian, splitting (N, n_r, l) levels into
j = l +/- 1/2 sub-levels.  Reproduces the real nuclear magic numbers:
2, 8, 20, 28, 50, 82, 126.

References:
  - Mayer & Jensen, "Elementary Theory of Nuclear Shell Structure" (1955)
  - GeoVac Paper 0 (packing construction), Paper 7 (S3 proof)
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy import sparse

from geovac.nuclear.harmonic_shell import (
    _ho_shell_degeneracy,
    _l_values_in_shell,
    build_nuclear_graph_harmonic,
)

# Spectroscopic letter for orbital angular momentum l
_L_LETTERS = "spdfghijklmnoqrtuvwxyz"


def _l_letter(l: int) -> str:
    """Return spectroscopic letter for angular momentum l."""
    if l < len(_L_LETTERS):
        return _L_LETTERS[l]
    return f"[{l}]"


# ---------------------------------------------------------------------------
# 1. l.s eigenvalue
# ---------------------------------------------------------------------------

def ls_eigenvalue(l: int, j: float) -> float:
    """
    Eigenvalue of l.s for given l and j.

    Parameters
    ----------
    l : int
        Orbital angular momentum quantum number (>= 0).
    j : float
        Total angular momentum quantum number (l +/- 0.5).

    Returns
    -------
    float
        <l.s> = l/2 for j = l + 1/2, -(l+1)/2 for j = l - 1/2.
        For l = 0, returns 0 (only j = 1/2 exists).
    """
    if l == 0:
        return 0.0
    if abs(j - (l + 0.5)) < 1e-10:
        return l / 2.0
    elif abs(j - (l - 0.5)) < 1e-10:
        return -(l + 1) / 2.0
    else:
        raise ValueError(f"Invalid j={j} for l={l}; must be {l+0.5} or {l-0.5}")


# ---------------------------------------------------------------------------
# 2. Nuclear shell levels (j-scheme)
# ---------------------------------------------------------------------------

def nuclear_shell_levels(
    n_max: int = 7,
    hw: float = 1.0,
    v_ls: float = 0.0,
    d_ll: float = 0.0,
) -> List[Dict]:
    """
    Build the full nuclear shell model energy spectrum in the j-scheme.

    The single-particle energy is (modified oscillator / Nilsson model):
        E(N, l, j) = hw*(N + 3/2) - v_ls * <l.s> - d_ll * l(l+1)

    The -d_ll * l(l+1) term compresses the effective shell spacing for
    high-l orbitals, mimicking the flattening of the nuclear mean-field
    potential (Woods-Saxon) relative to the pure HO.  Combined with
    spin-orbit splitting, this produces the observed nuclear magic numbers
    2, 8, 20, 28, 50, 82, 126 by allowing the j = l+1/2 member of each
    highest-l sub-shell to drop into the gap below.

    Parameters
    ----------
    n_max : int
        Number of HO shells (N = 0 .. n_max-1).
    hw : float
        Oscillator energy hbar*omega (Ha or MeV, units are consistent).
    v_ls : float
        Spin-orbit coupling strength.  H_SO = -v_ls * (l.s).
        Positive v_ls lowers j = l+1/2 levels (nuclear convention).
    d_ll : float
        Orbital correction strength.  Positive d_ll lowers high-l states.

    Returns
    -------
    list of dict
        Each dict has keys: 'N', 'n_r', 'l', 'j', 'energy', 'degeneracy',
        'label'.  Sorted by energy.
    """
    levels: List[Dict] = []
    for N in range(n_max):
        e_ho = hw * (N + 1.5)
        for l in _l_values_in_shell(N):
            n_r = (N - l) // 2
            if l == 0:
                j_values = [0.5]
            else:
                j_values = [l + 0.5, l - 0.5]
            for j in j_values:
                ls_val = ls_eigenvalue(l, j)
                energy = e_ho - v_ls * ls_val - d_ll * l * (l + 1)
                deg = int(2 * j + 1)
                # Label: n_r + letter(l) + subscript(j)
                j_str = f"{int(2*j)}/2"
                label = f"{n_r}{_l_letter(l)}{j_str}"
                levels.append({
                    "N": N,
                    "n_r": n_r,
                    "l": l,
                    "j": j,
                    "energy": energy,
                    "degeneracy": deg,
                    "label": label,
                })
    # Sort by energy, then by N (for stability), then by j descending
    levels.sort(key=lambda d: (d["energy"], d["N"], -d["j"]))
    return levels


# ---------------------------------------------------------------------------
# 3. Find magic numbers from energy gaps
# ---------------------------------------------------------------------------

def find_magic_numbers(
    n_max: int = 7,
    hw: float = 1.0,
    v_ls: float = 0.0,
    d_ll: float = 0.0,
    gap_threshold: Optional[float] = None,
    n_closures: Optional[int] = None,
) -> Dict:
    """
    Identify magic numbers from sorted energy levels via gap analysis.

    Parameters
    ----------
    n_max : int
        Number of HO shells.
    hw : float
        Oscillator energy.
    v_ls : float
        Spin-orbit coupling strength.
    d_ll : float
        Orbital l(l+1) correction strength.
    gap_threshold : float or None
        If given, magic numbers are cumulative occupancies where the gap
        exceeds this absolute threshold.  If None, uses the adaptive
        method (top ``n_closures`` gaps).
    n_closures : int or None
        Number of shell closures to identify (default: n_max for HO-like,
        7 for SO models).  Uses the ``n_closures`` largest gaps.

    Returns
    -------
    dict with keys:
        'magic_numbers' : list of int
        'levels'        : list of dict (the sorted levels)
        'gaps'          : list of (cumulative, gap) tuples
    """
    levels = nuclear_shell_levels(n_max, hw, v_ls, d_ll)

    # Accumulate occupancies and find gaps
    cumulative = 0
    gap_info: List[Tuple[int, float]] = []
    for i, lev in enumerate(levels):
        cumulative += lev["degeneracy"]
        if i < len(levels) - 1:
            gap = levels[i + 1]["energy"] - lev["energy"]
        else:
            gap = float("inf")
        gap_info.append((cumulative, gap))

    if gap_threshold is not None:
        # Absolute threshold method
        magic: List[int] = []
        for cum, gap in gap_info:
            if gap > gap_threshold:
                magic.append(cum)
    elif n_closures is not None:
        # Pick the top n_closures largest gaps (excluding final inf)
        finite_gaps = [(cum, gap) for cum, gap in gap_info if gap < float("inf")]
        sorted_by_gap = sorted(finite_gaps, key=lambda x: -x[1])
        top = sorted_by_gap[:n_closures]
        magic = sorted([cum for cum, gap in top])
    else:
        # Default: identify all gaps > 0 where cumulative count is
        # in the expected magic number set [2, 8, 20, 28, 50, 82, 126, 168, ...]
        # Fallback: use gap > 1.5 * median for automatic detection
        finite_gaps = [(cum, gap) for cum, gap in gap_info if gap < float("inf")]
        if finite_gaps:
            all_gaps_sorted = sorted([g for _, g in finite_gaps])
            median_gap = all_gaps_sorted[len(all_gaps_sorted) // 2]
            threshold = max(1.5 * median_gap, 0.01 * hw)
            magic = [cum for cum, gap in gap_info if gap > threshold]
        else:
            magic = []

    # Always include the final total
    total = sum(lev["degeneracy"] for lev in levels)
    if total not in magic:
        magic.append(total)

    return {
        "magic_numbers": magic,
        "levels": levels,
        "gaps": gap_info,
    }


def verify_magic_ordering(
    target: List[int],
    n_max: int = 7,
    hw: float = 1.0,
    v_ls: float = 0.0,
    d_ll: float = 0.0,
    min_gap: float = 1e-10,
) -> Dict:
    """
    Check whether the level ordering produces the target magic numbers.

    A magic number ``M`` is valid if:
      1. ``M`` appears as a cumulative occupancy in the sorted level list.
      2. There is a strictly positive energy gap (>= min_gap) at that point.

    This is the physically correct criterion: the magic numbers are the
    particle counts where sub-shell filling produces a gap in the
    single-particle spectrum, regardless of how large the gap is.

    Parameters
    ----------
    target : list of int
        Expected magic numbers, e.g. [2, 8, 20, 28, 50, 82, 126].
    n_max : int
        Number of HO shells.
    hw : float
        Oscillator energy.
    v_ls : float
        Spin-orbit coupling strength.
    d_ll : float
        Orbital correction strength.
    min_gap : float
        Minimum gap to count as a shell closure.

    Returns
    -------
    dict with keys:
        'valid'   : bool (True if all target magic numbers have positive gaps)
        'present' : list of int (target numbers that appear as cumulative counts)
        'gaps_at_target' : dict mapping magic number -> gap size
        'missing' : list of int (target numbers not found)
        'levels'  : the sorted level list
    """
    levels = nuclear_shell_levels(n_max, hw, v_ls, d_ll)

    # Build cumulative count -> gap mapping
    cumulative = 0
    cum_to_gap: Dict[int, float] = {}
    for i, lev in enumerate(levels):
        cumulative += lev["degeneracy"]
        if i < len(levels) - 1:
            gap = levels[i + 1]["energy"] - lev["energy"]
        else:
            gap = float("inf")
        cum_to_gap[cumulative] = gap

    gaps_at_target: Dict[int, float] = {}
    present: List[int] = []
    missing: List[int] = []
    for m in target:
        if m in cum_to_gap and cum_to_gap[m] >= min_gap:
            present.append(m)
            gaps_at_target[m] = cum_to_gap[m]
        else:
            missing.append(m)

    return {
        "valid": len(missing) == 0,
        "present": present,
        "gaps_at_target": gaps_at_target,
        "missing": missing,
        "levels": levels,
    }


# ---------------------------------------------------------------------------
# 4. Find optimal v_ls
# ---------------------------------------------------------------------------

_REAL_MAGIC = [2, 8, 20, 28, 50, 82, 126]


def find_optimal_vls(
    n_max: int = 7,
    hw: float = 1.0,
    v_ls_range: Tuple[float, float] = (0.0, 1.0),
    d_ll_range: Tuple[float, float] = (0.0, 0.3),
    n_scan: int = 200,
) -> Dict:
    """
    Search for v_ls and d_ll values that produce the real nuclear magic numbers.

    The Mayer-Jensen shell model uses:
        E = hw*(N + 3/2) - v_ls*(l.s) - d_ll*l(l+1)

    Both the spin-orbit coupling (v_ls) and the orbital correction (d_ll)
    are needed to reproduce all 7 magic numbers simultaneously.

    Parameters
    ----------
    n_max : int
        Number of HO shells (7 needed for magic 126).
    hw : float
        Oscillator energy.
    v_ls_range : tuple
        (min, max) of v_ls to scan.
    d_ll_range : tuple
        (min, max) of d_ll to scan.
    n_scan : int
        Number of scan points per parameter.

    Returns
    -------
    dict with keys:
        'optimal_vls'     : float (midpoint of valid range)
        'optimal_d_ll'    : float (midpoint of valid range)
        'vls_range'       : (min, max) that produce all 7 magic numbers
        'd_ll_range_found': (min, max) d_ll values found
        'ratio_vls_hw'    : optimal v_ls / hw
        'ratio_d_ll_hw'   : optimal d_ll / hw
        'magic_numbers'   : the magic numbers at optimal parameters
    """
    v_scan = np.linspace(v_ls_range[0], v_ls_range[1], n_scan)
    d_scan = np.linspace(d_ll_range[0], d_ll_range[1], n_scan)

    valid_params: List[Tuple[float, float]] = []
    for v in v_scan:
        for d in d_scan:
            result = verify_magic_ordering(_REAL_MAGIC, n_max, hw, v, d)
            if result["valid"]:
                valid_params.append((v, d))

    if not valid_params:
        return {
            "optimal_vls": None,
            "optimal_d_ll": None,
            "vls_range": (None, None),
            "d_ll_range_found": (None, None),
            "ratio_vls_hw": None,
            "ratio_d_ll_hw": None,
            "magic_numbers": [],
        }

    vls_vals = [p[0] for p in valid_params]
    dll_vals = [p[1] for p in valid_params]
    optimal_vls = (min(vls_vals) + max(vls_vals)) / 2.0
    optimal_dll = (min(dll_vals) + max(dll_vals)) / 2.0

    # Find the valid (v_ls, d_ll) closest to the center
    center = np.array([optimal_vls, optimal_dll])
    dists = [np.linalg.norm(np.array(p) - center) for p in valid_params]
    best_idx = int(np.argmin(dists))
    optimal_vls, optimal_dll = valid_params[best_idx]

    result = verify_magic_ordering(_REAL_MAGIC, n_max, hw, optimal_vls, optimal_dll)

    return {
        "optimal_vls": optimal_vls,
        "optimal_d_ll": optimal_dll,
        "vls_range": (min(vls_vals), max(vls_vals)),
        "d_ll_range_found": (min(dll_vals), max(dll_vals)),
        "ratio_vls_hw": optimal_vls / hw,
        "ratio_d_ll_hw": optimal_dll / hw,
        "magic_numbers": _REAL_MAGIC if result["valid"] else result["present"],
        "gaps_at_magic": result["gaps_at_target"],
    }


# ---------------------------------------------------------------------------
# 5. Build spin-orbit Hamiltonian matrix in (m_l, sigma) basis
# ---------------------------------------------------------------------------

def build_spin_orbit_hamiltonian(
    n_max: int = 7,
    hw: float = 1.0,
    v_ls: float = 0.0,
    d_ll: float = 0.0,
) -> Dict:
    """
    Build the full single-particle Hamiltonian with spin-orbit coupling.

    The basis is |N, n_r, l, m_l, sigma> where sigma=0 (spin up, +1/2)
    and sigma=1 (spin down, -1/2).

    H = H_HO + H_ll + H_SO
    H_HO = hw*(N + 3/2) * I       (diagonal)
    H_ll = -d_ll * l(l+1) * I     (diagonal, orbital correction)
    H_SO = -v_ls * (l . s)         (block-diagonal in (n_r, l) subspaces)

    l.s = lz*sz + (l+ s- + l- s+)/2

    Parameters
    ----------
    n_max : int
        Number of HO shells.
    hw : float
        Oscillator energy.
    v_ls : float
        Spin-orbit coupling strength.
    d_ll : float
        Orbital correction strength.

    Returns
    -------
    dict with keys:
        'hamiltonian'  : scipy.sparse.csr_matrix
        'states'       : list of (N, n_r, l, m_l, sigma) tuples
        'eigenvalues'  : 1D array (sorted)
        'n_states'     : int
    """
    # Build state list (same ordering as harmonic_shell)
    states: List[Tuple[int, int, int, int, int]] = []
    for N in range(n_max):
        for l in _l_values_in_shell(N):
            n_r = (N - l) // 2
            for m_l in range(-l, l + 1):
                for sigma in (0, 1):
                    states.append((N, n_r, l, m_l, sigma))

    n_states = len(states)

    # Build index lookup: (N, n_r, l, m_l, sigma) -> index
    state_to_idx: Dict[Tuple[int, int, int, int, int], int] = {}
    for i, s in enumerate(states):
        state_to_idx[s] = i

    # Diagonal: HO energies + orbital l(l+1) correction
    diag = np.array([hw * (s[0] + 1.5) - d_ll * s[2] * (s[2] + 1) for s in states])

    # Spin-orbit coupling: l.s in (m_l, sigma) basis
    # Build sparse triplets
    rows: List[int] = []
    cols: List[int] = []
    vals: List[float] = []

    for i, (N, n_r, l, m_l, sigma) in enumerate(states):
        if l == 0:
            continue  # no spin-orbit for s orbitals

        # --- lz * sz term (diagonal) ---
        # sz eigenvalue: +1/2 for sigma=0, -1/2 for sigma=1
        sz = 0.5 - sigma  # sigma=0 -> +0.5, sigma=1 -> -0.5
        lz_sz = m_l * sz
        val = -v_ls * lz_sz
        if abs(val) > 1e-15:
            rows.append(i)
            cols.append(i)
            vals.append(val)

        # --- (l+ s- + l- s+) / 2 term (off-diagonal in spin) ---
        # l+ s-: |m_l, up> -> |m_l+1, down>
        #   matrix element = (1/2) * sqrt((l - m_l)(l + m_l + 1))
        if sigma == 0 and m_l + 1 <= l:
            # This state is |m_l, up>; connects to |m_l+1, down>
            target = (N, n_r, l, m_l + 1, 1)
            j_idx = state_to_idx.get(target)
            if j_idx is not None:
                me = 0.5 * np.sqrt((l - m_l) * (l + m_l + 1))
                val = -v_ls * me
                if abs(val) > 1e-15:
                    rows.append(i)
                    cols.append(j_idx)
                    vals.append(val)
                    # Hermitian conjugate
                    rows.append(j_idx)
                    cols.append(i)
                    vals.append(val)

        # l- s+: |m_l, down> -> |m_l-1, up>
        #   matrix element = (1/2) * sqrt((l + m_l)(l - m_l + 1))
        if sigma == 1 and m_l - 1 >= -l:
            target = (N, n_r, l, m_l - 1, 0)
            j_idx = state_to_idx.get(target)
            if j_idx is not None:
                me = 0.5 * np.sqrt((l + m_l) * (l - m_l + 1))
                val = -v_ls * me
                if abs(val) > 1e-15:
                    # Only add if not already covered by the l+ s- branch above
                    # The l+ s- branch handles (m_l, up) -> (m_l+1, down)
                    # The l- s+ branch handles (m_l, down) -> (m_l-1, up)
                    # These are different pairs, so check for duplication
                    # Actually, the pair (m_l-1, up) <-> (m_l, down) is the
                    # SAME pair as (m', up) <-> (m'+1, down) with m' = m_l-1.
                    # So the l+ s- branch already handled this when processing
                    # state (m_l-1, up). Skip to avoid double-counting.
                    pass

    # Build sparse SO matrix
    H_so = sparse.csr_matrix(
        (vals, (rows, cols)), shape=(n_states, n_states)
    )

    # Full Hamiltonian
    H = sparse.diags(diag) + H_so
    H = sparse.csr_matrix(H)

    # Diagonalize
    eigenvalues = np.sort(np.linalg.eigvalsh(H.toarray()))

    return {
        "hamiltonian": H,
        "states": states,
        "eigenvalues": eigenvalues,
        "n_states": n_states,
    }


# ---------------------------------------------------------------------------
# 6. Level ordering table
# ---------------------------------------------------------------------------

def level_ordering_table(
    n_max: int = 7,
    hw: float = 1.0,
    v_ls: Optional[float] = None,
    d_ll: Optional[float] = None,
    n_closures: Optional[int] = 7,
) -> str:
    """
    Pretty-print the nuclear shell model level ordering.

    Parameters
    ----------
    n_max : int
        Number of HO shells.
    hw : float
        Oscillator energy.
    v_ls : float or None
        Spin-orbit coupling. If None, uses optimal value.
    d_ll : float or None
        Orbital correction. If None, uses optimal value.
    n_closures : int or None
        Number of shell closures to identify.

    Returns
    -------
    str
        Formatted table.
    """
    if v_ls is None or d_ll is None:
        opt = find_optimal_vls(n_max, hw)
        if v_ls is None:
            v_ls = opt["optimal_vls"]
        if d_ll is None:
            d_ll = opt["optimal_d_ll"]
        if v_ls is None:
            return "No optimal v_ls found."

    # Identify magic numbers via ordering verification
    vresult = verify_magic_ordering(_REAL_MAGIC, n_max, hw, v_ls, d_ll)
    levels = vresult["levels"]
    magic = set(vresult["present"])

    lines = []
    lines.append(
        f"Nuclear Shell Model Level Ordering "
        f"(v_ls/hw = {v_ls/hw:.4f}, d_ll/hw = {d_ll/hw:.4f})"
    )
    lines.append("=" * 72)
    lines.append(
        f"{'Label':>10s}  {'N':>3s}  {'l':>3s}  {'j':>5s}  "
        f"{'Energy':>10s}  {'Deg':>4s}  {'Cumul':>6s}  {'Magic':>6s}"
    )
    lines.append("-" * 72)

    cumulative = 0
    for lev in levels:
        cumulative += lev["degeneracy"]
        j_str = f"{int(2*lev['j'])}/2"
        is_magic = "*" if cumulative in magic else ""
        lines.append(
            f"{lev['label']:>10s}  {lev['N']:3d}  {lev['l']:3d}  "
            f"{j_str:>5s}  {lev['energy']:10.4f}  "
            f"{lev['degeneracy']:4d}  {cumulative:6d}  {is_magic:>6s}"
        )
        if cumulative in magic:
            lines.append(f"{'--- Shell closure: ' + str(cumulative) + ' ---':^72s}")

    lines.append("=" * 72)
    lines.append(f"Magic numbers: {sorted(magic)}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 7. Sparsity with spin-orbit
# ---------------------------------------------------------------------------

def sparsity_with_spin_orbit(
    n_max: int = 4,
    hw: float = 1.0,
    v_ls: float = 0.0,
    d_ll: float = 0.0,
) -> Dict:
    """
    Measure how spin-orbit coupling affects Hamiltonian sparsity.

    Parameters
    ----------
    n_max : int
        Number of HO shells.
    hw : float
        Oscillator energy.
    v_ls : float
        Spin-orbit coupling strength.

    Returns
    -------
    dict with sparsity metrics.
    """
    # Without SO: pure diagonal HO Hamiltonian
    so_data_0 = build_spin_orbit_hamiltonian(n_max, hw, v_ls=0.0, d_ll=0.0)
    H_no_so = so_data_0["hamiltonian"]
    n_states = so_data_0["n_states"]
    total_elements = n_states * n_states

    threshold = 1e-15
    nnz_no_so_true = int(np.sum(np.abs(H_no_so.toarray()) > threshold))

    # With SO
    so_data = build_spin_orbit_hamiltonian(n_max, hw, v_ls, d_ll)
    H_with_so = so_data["hamiltonian"]
    nnz_with_so_true = int(np.sum(np.abs(H_with_so.toarray()) > threshold))

    # One-body sparsity
    h1_density_no_so = nnz_no_so_true / total_elements
    h1_density_with_so = nnz_with_so_true / total_elements

    return {
        "n_max": n_max,
        "n_states": n_states,
        "total_elements": total_elements,
        "without_so": {
            "nonzero_elements": int(nnz_no_so_true),
            "density": h1_density_no_so,
            "sparsity": 1.0 - h1_density_no_so,
        },
        "with_so": {
            "nonzero_elements": int(nnz_with_so_true),
            "density": h1_density_with_so,
            "sparsity": 1.0 - h1_density_with_so,
        },
        "so_adds_elements": int(nnz_with_so_true - nnz_no_so_true),
        "note": (
            "Spin-orbit adds off-diagonal (m_l, sigma) couplings within "
            "each (n_r, l) subspace. Two-body ERIs are unchanged because "
            "SO is a one-body operator."
        ),
    }
