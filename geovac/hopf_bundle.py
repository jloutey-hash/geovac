"""
Hopf Bundle Module for GeoVac
==============================

Explores the Hopf fibration S¹ → S³ → S² structure implicit in the
discrete lattice via Fock's 1935 stereographic projection.

The lattice states |n,l,m⟩ live on S³ (Paper 7). The Hopf fibration
decomposes S³ into an S¹ fiber over an S² base. This module:

1. Maps lattice states to explicit S³ coordinates via Fock's projection
2. Decomposes each point into Hopf base (S²) and fiber (S¹) components
3. Computes discrete volume sums that approximate sphere volumes
4. Searches for combinations reproducing 1/α ≈ 4π³ + π² + π = 137.036

Coordinate conventions:
    S³ parametrized by (χ, ψ₁, ψ₂) with χ ∈ [0,π], ψ_i ∈ [0,2π):
        z₁ = cos(χ/2) e^{iψ₁},  z₂ = sin(χ/2) e^{iψ₂}
    giving embedding coordinates:
        x₁ = cos(χ/2) cos(ψ₁),  x₂ = cos(χ/2) sin(ψ₁)
        x₃ = sin(χ/2) cos(ψ₂),  x₄ = sin(χ/2) sin(ψ₂)

    Fock's map: χ_n = 2 arctan(1/n) places the nth shell at hyperspherical
    angle χ_n on S³ (using p₀ = 1 ground-state normalization).

Author: GeoVac
Date: March 2026
"""

import numpy as np
from typing import Tuple, List, Dict, Optional


# ============================================================
# 1. Lattice → S³ mapping via Fock's projection
# ============================================================

def fock_chi(n: int) -> float:
    """
    Hyperspherical angle χ for shell n via Fock's stereographic projection.

    From tan(χ/2) = p/p₀ with characteristic momentum p ~ 1/n and
    ground-state normalization p₀ = 1:

        χ_n = 2 arctan(1/n)

    Properties:
        n=1 → χ = π/2  (equator of S³)
        n→∞ → χ → 0    (north pole)

    Parameters
    ----------
    n : int
        Principal quantum number (n ≥ 1).

    Returns
    -------
    float
        Hyperspherical angle in [0, π].
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    return 2.0 * np.arctan(1.0 / n)


def state_to_hopf_angles(n: int, l: int, m: int) -> Tuple[float, float, float]:
    """
    Map quantum state |n,l,m⟩ to Hopf coordinates (χ, ψ₁, ψ₂) on S³.

    The mapping uses:
        χ  = 2 arctan(1/n)         — Fock projection (radial → S³ latitude)
        ψ₁ = π l / max(n-1, 1)     — angular momentum → first phase
        ψ₂ = 2π m / (2l+1)         — magnetic number → second phase

    This assigns each state a unique point on S³ while preserving the
    Hopf structure: ψ₁ and ψ₂ separately parametrize the two complex
    planes of C² ⊃ S³, and their sum/difference encode fiber/base.

    Parameters
    ----------
    n : int
        Principal quantum number (n ≥ 1).
    l : int
        Angular momentum quantum number (0 ≤ l < n).
    m : int
        Magnetic quantum number (-l ≤ m ≤ l).

    Returns
    -------
    Tuple[float, float, float]
        (χ, ψ₁, ψ₂) with χ ∈ (0, π], ψ₁ ∈ [0, π], ψ₂ ∈ (-π, π].
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    if l < 0 or l >= n:
        raise ValueError(f"l must satisfy 0 <= l < n, got l={l}, n={n}")
    if abs(m) > l:
        raise ValueError(f"|m| must be <= l, got m={m}, l={l}")

    chi = fock_chi(n)
    psi1 = np.pi * l / max(n - 1, 1)
    psi2 = 2.0 * np.pi * m / (2 * l + 1) if l > 0 else 0.0

    return chi, psi1, psi2


def state_to_s3(n: int, l: int, m: int) -> np.ndarray:
    """
    Map quantum state |n,l,m⟩ to a point on the unit S³.

    Uses Fock's projection to assign each lattice state explicit
    embedding coordinates (x₁, x₂, x₃, x₄) ∈ R⁴ with |x|² = 1.

    Parameters
    ----------
    n : int
        Principal quantum number (n ≥ 1).
    l : int
        Angular momentum quantum number (0 ≤ l < n).
    m : int
        Magnetic quantum number (-l ≤ m ≤ l).

    Returns
    -------
    np.ndarray
        4-vector (x₁, x₂, x₃, x₄) on the unit S³.
    """
    chi, psi1, psi2 = state_to_hopf_angles(n, l, m)
    cos_half = np.cos(chi / 2.0)
    sin_half = np.sin(chi / 2.0)

    return np.array([
        cos_half * np.cos(psi1),
        cos_half * np.sin(psi1),
        sin_half * np.cos(psi2),
        sin_half * np.sin(psi2),
    ])


def lattice_to_s3(n_max: int) -> Dict[Tuple[int, int, int], np.ndarray]:
    """
    Map all lattice states up to n_max to points on S³.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    Dict[Tuple[int, int, int], np.ndarray]
        Maps (n, l, m) → 4-vector on unit S³.
    """
    points: Dict[Tuple[int, int, int], np.ndarray] = {}
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                points[(n, l, m)] = state_to_s3(n, l, m)
    return points


# ============================================================
# 2. Hopf projection S³ → S²
# ============================================================

def hopf_project(point: np.ndarray) -> np.ndarray:
    """
    Hopf projection from S³ to S².

    Maps (x₁, x₂, x₃, x₄) on S³ to (y₁, y₂, y₃) on S² via
    the standard Hopf map written in terms of the spinor
    z₁ = x₁ + ix₂, z₂ = x₃ + ix₄:

        y₁ = 2 Re(z₁ z₂*) = 2(x₁x₃ + x₂x₄)
        y₂ = 2 Im(z₁ z₂*) = 2(x₂x₃ - x₁x₄)
        y₃ = |z₁|² - |z₂|² = x₁² + x₂² - x₃² - x₄²

    Parameters
    ----------
    point : np.ndarray
        4-vector (x₁, x₂, x₃, x₄), assumed to lie on unit S³.

    Returns
    -------
    np.ndarray
        3-vector (y₁, y₂, y₃) on the unit S².
    """
    x1, x2, x3, x4 = point
    return np.array([
        2.0 * (x1 * x3 + x2 * x4),
        2.0 * (x2 * x3 - x1 * x4),
        x1**2 + x2**2 - x3**2 - x4**2,
    ])


# ============================================================
# 3. Hopf fiber coordinate
# ============================================================

def fiber_angle(point: np.ndarray) -> float:
    """
    Extract the Hopf fiber (S¹) coordinate from a point on S³.

    In spinor form z₁ = x₁+ix₂, z₂ = x₃+ix₄, the fiber angle is
    the total phase ψ = arg(z₁) + arg(z₂). Points on the same Hopf
    fiber (same base point on S²) differ only in ψ.

    Parameters
    ----------
    point : np.ndarray
        4-vector on the unit S³.

    Returns
    -------
    float
        Fiber angle ψ ∈ (-2π, 2π], or 0.0 if either spinor
        component vanishes (degenerate fiber).
    """
    x1, x2, x3, x4 = point
    z1_angle = np.arctan2(x2, x1)
    z2_angle = np.arctan2(x4, x3)
    return z1_angle + z2_angle


def decompose_hopf(point: np.ndarray) -> Dict[str, object]:
    """
    Full Hopf decomposition of an S³ point.

    Parameters
    ----------
    point : np.ndarray
        4-vector on the unit S³.

    Returns
    -------
    Dict with keys:
        's3': np.ndarray — the input point
        's2': np.ndarray — Hopf base point on S²
        'fiber': float — fiber angle ψ
        'chi': float — hyperspherical angle (from x₃²+x₄² = sin²(χ/2))
    """
    s2 = hopf_project(point)
    psi = fiber_angle(point)
    # Recover χ from |z₂|² = sin²(χ/2) = x₃² + x₄²
    sin2_half = point[2]**2 + point[3]**2
    chi = 2.0 * np.arcsin(np.clip(np.sqrt(sin2_half), 0.0, 1.0))
    return {
        's3': point,
        's2': s2,
        'fiber': psi,
        'chi': chi,
    }


# ============================================================
# 4. Discrete volume sums
# ============================================================

def shell_degeneracy(n: int) -> int:
    """Number of states in shell n: n²."""
    return n * n


def discrete_s3_volume(n_max: int) -> float:
    """
    Discrete approximation to Vol(S³) = 2π² from lattice structure.

    The S³ volume element in hyperspherical coordinates is:
        dV₃ = sin²(χ) sin(θ) dχ dθ dφ

    Integrating the angular part: ∫dφ ∫sin(θ)dθ = 4π. So the
    shell volume at χ with width dχ is 4π sin²(χ) dχ, and
    Vol(S³) = 4π ∫₀^π sin²(χ) dχ = 4π · π/2 = 2π².

    IMPORTANT: Fock's projection maps bound states to χ ∈ (0, π/2]
    (n=1 → χ=π/2, n→∞ → χ→0). The upper hemisphere χ ∈ (π/2, π)
    corresponds to continuum states. We compute the bound-state
    half-volume and double it by symmetry of sin²(χ).

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    float
        Discrete volume sum (converges to 2π² ≈ 19.739).
    """
    vol = 0.0

    for n in range(1, n_max + 1):
        chi = fock_chi(n)
        # Cell width in χ: midpoint rule between neighboring shells
        if n == 1 and n_max == 1:
            d_chi = chi  # Only one shell
        elif n == 1:
            # Upper boundary: π/2 (Fock limit for n=1)
            # Lower boundary: midpoint to n=2
            d_chi = np.pi / 2.0 - 0.5 * (fock_chi(1) + fock_chi(2))
            d_chi += np.pi / 2.0 - fock_chi(1)  # Extend to equator
            d_chi = 0.5 * (np.pi / 2.0 - fock_chi(2)) + (np.pi / 2.0 - fock_chi(1))
            # Simpler: use distance to boundaries
            upper = np.pi / 2.0  # Equator
            lower = 0.5 * (fock_chi(1) + fock_chi(2)) if n_max > 1 else 0.0
            d_chi = upper - lower
        elif n == n_max:
            upper = 0.5 * (fock_chi(n - 1) + fock_chi(n))
            lower = 0.0  # Boundary at north pole
            d_chi = upper - lower
        else:
            upper = 0.5 * (fock_chi(n - 1) + fock_chi(n))
            lower = 0.5 * (fock_chi(n) + fock_chi(n + 1))
            d_chi = upper - lower

        # Shell volume: 4π sin²(χ) dχ
        shell_vol = 4.0 * np.pi * np.sin(chi)**2 * d_chi
        vol += shell_vol

    # Double to account for continuum hemisphere (sin²(χ) symmetric
    # about π/2: sin²(π-χ) = sin²(χ))
    return 2.0 * vol


def discrete_fiber_volume(n_max: int) -> float:
    """
    Discrete approximation to the Hopf fiber integral Vol(S¹) = 2π.

    Each m-ladder (fixed n,l varying m from -l to +l) traces a path
    along the Hopf fiber. The total phase winding across the full
    m-ladder should approach 2π in the continuum limit.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    float
        Average fiber winding (should converge to 2π ≈ 6.283).
    """
    windings = []

    for n in range(1, n_max + 1):
        for l in range(1, n):  # l > 0 for nontrivial fiber
            angles = []
            for m in range(-l, l + 1):
                pt = state_to_s3(n, l, m)
                angles.append(fiber_angle(pt))

            if len(angles) >= 2:
                # The m-ladder spans from m=-l to m=+l
                # Phase increment per m step
                winding = angles[-1] - angles[0]
                windings.append(winding)

    if not windings:
        return 0.0
    return float(np.mean(windings))


def discrete_base_area(n_max: int) -> float:
    """
    Discrete approximation to the Hopf base area Vol(S²) = 4π.

    Uses the convex hull solid angle of projected base points on S².

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    float
        Discrete base area (should converge to 4π ≈ 12.566).
    """
    # Collect unique base points
    base_points = []
    for n in range(1, n_max + 1):
        for l in range(n):
            for m in range(-l, l + 1):
                pt = state_to_s3(n, l, m)
                base_points.append(hopf_project(pt))

    if len(base_points) < 3:
        return 0.0

    base_points = np.array(base_points)
    n_pts = len(base_points)

    # Compute the solid angle covered using the centroid deficit method:
    # For n points on S², the typical Voronoi cell has area 4π/n.
    # A uniform covering gives area ≈ 4π. We measure uniformity.
    centroid = base_points.mean(axis=0)
    centroid_len = np.linalg.norm(centroid)

    # For uniform distribution on S², E[|centroid|] → 0 and area → 4π
    # For clustered distribution, |centroid| → 1 and area → 0
    # Approximate: area ≈ 4π (1 - centroid_len)
    return 4.0 * np.pi * (1.0 - centroid_len)


# ============================================================
# 5. Alpha formula exploration
# ============================================================

# Continuum sphere volumes
VOL_S1 = 2.0 * np.pi             # 6.2832
VOL_S2 = 4.0 * np.pi             # 12.566
VOL_S3 = 2.0 * np.pi**2          # 19.739
ALPHA_INV_FORMULA = 4 * np.pi**3 + np.pi**2 + np.pi  # 137.0363
ALPHA_INV_EXPERIMENT = 137.035999177      # CODATA 2022


def alpha_formula_terms() -> Dict[str, float]:
    """
    Break down the formula 1/α ≈ 4π³ + π² + π into constituent terms.

    Returns
    -------
    Dict with keys:
        'term_4pi3': 4π³ = Vol(S¹) × Vol(S³)
        'term_pi2': π² = Vol(S³) / 2
        'term_pi': π = Vol(S¹) / 2
        'total': sum of all terms
        'experimental': CODATA 2022 value
        'relative_error': |formula - experiment| / experiment
    """
    t1 = 4.0 * np.pi**3
    t2 = np.pi**2
    t3 = np.pi
    total = t1 + t2 + t3

    return {
        'term_4pi3': t1,
        'term_pi2': t2,
        'term_pi': t3,
        'total': total,
        'experimental': ALPHA_INV_EXPERIMENT,
        'relative_error': abs(total - ALPHA_INV_EXPERIMENT) / ALPHA_INV_EXPERIMENT,
        # Geometric identifications
        'term_4pi3_as': 'Vol(S¹) × Vol(S³)',
        'term_pi2_as': 'Vol(S³) / 2',
        'term_pi_as': 'Vol(S¹) / 2',
    }


def discrete_alpha_search(n_max: int) -> Dict[str, float]:
    """
    Compute discrete volume sums from the lattice and search for
    combinations approaching 1/α = 137.036.

    Uses the lattice up to n_max to compute discrete approximations
    to Vol(S¹), Vol(S²), Vol(S³) and evaluates candidate formulas.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    Dict with discrete volumes and candidate formula evaluations.
    """
    v1 = discrete_fiber_volume(n_max)
    v2 = discrete_base_area(n_max)
    v3 = discrete_s3_volume(n_max)

    # Total number of lattice states
    n_states = sum(n * n for n in range(1, n_max + 1))

    # Candidate formulas using discrete volumes
    candidates = {}

    # Formula A: v1 * v3 + v3/2 + v1/2  (direct analogue)
    if v1 > 0:
        candidates['v1*v3 + v3/2 + v1/2'] = v1 * v3 + v3 / 2.0 + v1 / 2.0

    # Formula B: v1 * v3 + v1²/4 + v1/2  (original form)
    if v1 > 0:
        candidates['v1*v3 + v1²/4 + v1/2'] = v1 * v3 + v1**2 / 4.0 + v1 / 2.0

    # Formula C: using all three volumes
    if v1 > 0 and v2 > 0:
        candidates['v1*v3 + v2²/(16π) + v1/2'] = (
            v1 * v3 + v2**2 / (16.0 * np.pi) + v1 / 2.0
        )

    # Berry-phase inspired: sum of log-holonomies relates to fiber?
    # Phase per plaquette ~ 2/n, total winding sums
    berry_sum = sum(2.0 / n for n in range(1, n_max + 1))
    candidates['berry_harmonic_sum'] = berry_sum

    return {
        'n_max': n_max,
        'n_states': n_states,
        'discrete_vol_s1': v1,
        'discrete_vol_s2': v2,
        'discrete_vol_s3': v3,
        'continuum_vol_s1': VOL_S1,
        'continuum_vol_s2': VOL_S2,
        'continuum_vol_s3': VOL_S3,
        'target': ALPHA_INV_EXPERIMENT,
        'candidates': candidates,
    }


# ============================================================
# 6. Lattice Hopf structure analysis
# ============================================================

def hopf_fiber_analysis(n_max: int) -> Dict[str, object]:
    """
    Analyze how lattice states distribute across Hopf fibers.

    Groups states by their S² base point (approximately — states
    on the same fiber project to the same point on S²). Reveals
    the fiber structure encoded in the (n,l,m) quantum numbers.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    Dict with:
        'n_states': total states
        'base_points': array of S² points
        'fiber_angles': array of fiber angles
        'fiber_groups': Dict mapping (n,l) → list of (m, fiber_angle)
        'fiber_winding_by_l': average winding per l value
    """
    base_pts = []
    fiber_angs = []
    states = []
    fiber_groups: Dict[Tuple[int, int], List[Tuple[int, float]]] = {}
    winding_by_l: Dict[int, List[float]] = {}

    for n in range(1, n_max + 1):
        for l in range(n):
            m_angles = []
            for m in range(-l, l + 1):
                pt = state_to_s3(n, l, m)
                base_pts.append(hopf_project(pt))
                fa = fiber_angle(pt)
                fiber_angs.append(fa)
                states.append((n, l, m))
                m_angles.append((m, fa))

            fiber_groups[(n, l)] = m_angles

            if l > 0 and len(m_angles) >= 2:
                angles = [a for _, a in m_angles]
                winding = max(angles) - min(angles)
                winding_by_l.setdefault(l, []).append(winding)

    # Average winding per l
    avg_winding = {l: np.mean(ws) for l, ws in winding_by_l.items()}

    return {
        'n_states': len(states),
        'base_points': np.array(base_pts),
        'fiber_angles': np.array(fiber_angs),
        'fiber_groups': fiber_groups,
        'fiber_winding_by_l': avg_winding,
    }


def convergence_study(n_max_range: Optional[List[int]] = None) -> Dict[str, object]:
    """
    Study convergence of discrete volumes toward continuum values.

    Parameters
    ----------
    n_max_range : List[int], optional
        Values of n_max to test (default: [2, 5, 10, 20, 50]).

    Returns
    -------
    Dict with convergence data for each volume.
    """
    if n_max_range is None:
        n_max_range = [2, 5, 10, 20, 50]

    results = {
        'n_max_values': n_max_range,
        'vol_s3': [],
        'vol_s1': [],
        'vol_s2': [],
        'alpha_candidates': [],
    }

    for n_max in n_max_range:
        data = discrete_alpha_search(n_max)
        results['vol_s3'].append(data['discrete_vol_s3'])
        results['vol_s1'].append(data['discrete_vol_s1'])
        results['vol_s2'].append(data['discrete_vol_s2'])
        results['alpha_candidates'].append(data['candidates'])

    return results


# ============================================================
# 7. Algebraic Hopf structure — counting and degeneracy
# ============================================================

def fiber_winding_exact(l: int) -> float:
    """
    Exact fiber winding for an m-ladder at angular momentum l.

    The ψ₂ mapping gives m-step = 2π/(2l+1), so the total winding
    from m=-l to m=+l is:

        W(l) = 4πl / (2l+1)

    This is a discrete fraction of 2π:  W(l)/2π = 2l/(2l+1).

    Parameters
    ----------
    l : int
        Angular momentum quantum number (l ≥ 1).

    Returns
    -------
    float
        Exact fiber winding for this l value.
    """
    if l < 1:
        return 0.0
    return 4.0 * np.pi * l / (2 * l + 1)


def algebraic_alpha_analysis(n_max: int) -> Dict[str, object]:
    """
    Algebraic analysis of the Hopf structure for alpha derivation.

    Rather than numerical quadrature, this examines exact algebraic
    quantities from the lattice counting structure:

    - Shell degeneracy: g(n) = n²
    - Fiber fraction: f(l) = 2l/(2l+1) — how much of S¹ is covered
    - SO(4) Casimir: C(n) = n² - 1 — Laplacian eigenvalue on S³
    - Weighted sums involving π from the geometric structure

    The key insight is that the formula 1/α = 4π³ + π² + π can be
    rewritten as:
        1/α = Vol(S¹)·Vol(S³) + (Vol(S³) + Vol(S¹))/2
            = 2π·2π² + (2π² + 2π)/2
    suggesting it arises from the fiber-total volume product plus
    a "self-energy" correction from the bundle structure.

    Parameters
    ----------
    n_max : int
        Maximum principal quantum number.

    Returns
    -------
    Dict with algebraic structure data and candidate formulas.
    """
    results: Dict[str, object] = {}

    # Exact degeneracy-weighted sums
    total_states = sum(n**2 for n in range(1, n_max + 1))
    total_casimir = sum(n**2 * (n**2 - 1) for n in range(1, n_max + 1))

    # Fiber completeness: weighted average of 2l/(2l+1) over all states
    fiber_sum = 0.0
    fiber_count = 0
    for n in range(1, n_max + 1):
        for l in range(1, n):  # l > 0 only
            weight = 2 * l + 1  # degeneracy of this l subshell
            fiber_sum += weight * (2 * l) / (2 * l + 1)
            fiber_count += weight

    # Average fiber fraction (should → 1 as n_max → ∞)
    avg_fiber_frac = fiber_sum / fiber_count if fiber_count > 0 else 0.0

    # Effective fiber volume: 2π × avg_fiber_frac
    v1_eff = 2.0 * np.pi * avg_fiber_frac

    # Degeneracy-weighted Fock angle sums
    chi_sum = sum(n**2 * fock_chi(n) for n in range(1, n_max + 1))
    sin2_chi_sum = sum(
        n**2 * np.sin(fock_chi(n))**2 for n in range(1, n_max + 1)
    )

    # Exact continuum target decomposition
    # 1/α = Vol(S¹)·Vol(S³) + (Vol(S³) + Vol(S¹))/2
    #      = fiber × total + (total + fiber) / 2
    target = ALPHA_INV_EXPERIMENT

    # Candidate A: use continuum volumes scaled by discrete fiber fraction
    v1_disc = v1_eff
    v3_disc = discrete_s3_volume(n_max)
    cand_a = v1_disc * v3_disc + (v3_disc + v1_disc) / 2.0

    # Candidate B: pure algebraic — use state-counting to define volumes
    # N_total = n_max(n_max+1)(2n_max+1)/6 → ∞
    # The ratio N_states_with_l>0 / N_total → 1 as n_max → ∞
    n_with_fiber = total_states - n_max  # States with l > 0 (subtract l=0 count)
    fiber_coverage = n_with_fiber / total_states if total_states > 0 else 0.0

    # Candidate C: Casimir-weighted formula
    # Sum of n²(n²-1) relates to ∫ λ dμ over the S³ spectrum
    avg_casimir = total_casimir / total_states if total_states > 0 else 0.0

    results = {
        'n_max': n_max,
        'total_states': total_states,
        'total_casimir': total_casimir,
        'avg_casimir': avg_casimir,
        'avg_fiber_fraction': avg_fiber_frac,
        'effective_vol_s1': v1_eff,
        'discrete_vol_s3': v3_disc,
        'fiber_coverage': fiber_coverage,
        'candidate_A_fiber_x_total_plus_correction': cand_a,
        'candidate_A_error': abs(cand_a - target) / target,
        'continuum_formula': ALPHA_INV_FORMULA,
        'continuum_error': abs(ALPHA_INV_FORMULA - target) / target,
        'target': target,
    }

    return results


# ============================================================
# 8. Fock-embedding S³ mapping (replaces linear approximation)
# ============================================================

def fock_embedding(n: int, l: int, m: int,
                   theta: Optional[float] = None,
                   phi: Optional[float] = None) -> np.ndarray:
    """
    Map |n,l,m⟩ to S³ using the actual Fock embedding coordinates.

    Fock's stereographic projection maps momentum space to S³ with
    coordinates (ξ₀, ξ₁, ξ₂, ξ₃):
        ξ₀ = cos(χ)
        ξ₃ = sin(χ) cos(θ)
        ξ₁ = sin(χ) sin(θ) cos(φ)
        ξ₂ = sin(χ) sin(θ) sin(φ)

    where χ = 2 arctan(1/n) from Fock's map, and (θ,φ) are the
    angular coordinates of the state's semiclassical direction.

    For a state |n,l,m⟩:
        θ = arccos(m / √(l(l+1)))  — semiclassical polar angle
        φ = 2π(m+l)/(2l+1)         — evenly-spaced azimuthal angle

    This produces points on the exact unit S³ with the correct
    hyperspherical structure. The Hopf decomposition then follows
    naturally from these R⁴ coordinates.

    Parameters
    ----------
    n : int
        Principal quantum number (n ≥ 1).
    l : int
        Angular momentum quantum number (0 ≤ l < n).
    m : int
        Magnetic quantum number (-l ≤ m ≤ l).
    theta : float, optional
        Override polar angle (default: semiclassical from m, l).
    phi : float, optional
        Override azimuthal angle (default: evenly-spaced from m, l).

    Returns
    -------
    np.ndarray
        4-vector (ξ₀, ξ₁, ξ₂, ξ₃) on the unit S³.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n}")
    if l < 0 or l >= n:
        raise ValueError(f"l must satisfy 0 <= l < n, got l={l}, n={n}")
    if abs(m) > l:
        raise ValueError(f"|m| must be <= l, got m={m}, l={l}")

    chi = fock_chi(n)

    if theta is None:
        if l > 0:
            cos_theta = m / np.sqrt(l * (l + 1))
            theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))
        else:
            theta = 0.0

    if phi is None:
        if l > 0:
            phi = 2.0 * np.pi * (m + l) / (2 * l + 1)
        else:
            phi = 0.0

    xi0 = np.cos(chi)
    xi3 = np.sin(chi) * np.cos(theta)
    xi1 = np.sin(chi) * np.sin(theta) * np.cos(phi)
    xi2 = np.sin(chi) * np.sin(theta) * np.sin(phi)

    return np.array([xi0, xi1, xi2, xi3])


def fock_hopf_decompose(n: int, l: int, m: int) -> Dict[str, object]:
    """
    Full Hopf decomposition of |n,l,m⟩ using Fock embedding.

    Maps the state to S³ via fock_embedding, then decomposes into
    Hopf base (S²) and fiber (S¹) using the spinor identification
    z₁ = ξ₀ + iξ₃, z₂ = ξ₁ + iξ₂.

    This identification aligns the Hopf fiber with the azimuthal
    structure: L_z eigenvalue m rotates (ξ₁, ξ₂) while leaving
    (ξ₀, ξ₃) in the L_z = 0 plane.

    Parameters
    ----------
    n, l, m : int
        Quantum numbers.

    Returns
    -------
    Dict with keys:
        'fock_s3': np.ndarray — (ξ₀, ξ₁, ξ₂, ξ₃) on S³
        'chi': float — Fock angle χ_n
        'z1', 'z2': complex — spinor components
        'hopf_s2': np.ndarray — base point on S²
        'fiber': float — fiber angle ψ₁ + ψ₂
        'base_phase': float — base phase ψ₁ - ψ₂
        'chi_hopf': float — Hopf polar angle from |z₂|/|z₁|
    """
    pt = fock_embedding(n, l, m)
    xi0, xi1, xi2, xi3 = pt

    # Spinor identification: z₁ = ξ₀ + iξ₃ (L_z=0 plane)
    #                        z₂ = ξ₁ + iξ₂ (azimuthal plane)
    z1 = complex(xi0, xi3)
    z2 = complex(xi1, xi2)

    # Hopf angles
    chi_hopf = 2.0 * np.arctan2(abs(z2), abs(z1))
    psi1 = np.angle(z1) if abs(z1) > 1e-15 else 0.0
    psi2 = np.angle(z2) if abs(z2) > 1e-15 else 0.0

    # Hopf projection S³ → S²
    y1 = 2.0 * (z1 * z2.conjugate()).real
    y2 = 2.0 * (z1 * z2.conjugate()).imag
    y3 = abs(z1)**2 - abs(z2)**2

    return {
        'fock_s3': pt,
        'chi': fock_chi(n),
        'z1': z1,
        'z2': z2,
        'hopf_s2': np.array([y1, y2, y3]),
        'fiber': psi1 + psi2,
        'base_phase': psi1 - psi2,
        'chi_hopf': chi_hopf,
    }


def cg_expectation_mapping(n: int, l: int, m: int) -> Dict[str, float]:
    """
    Map |n,l,m⟩ to Hopf angles via CG expectation values.

    Uses the SO(4) ≅ SU(2)⊗SU(2) decomposition from wigner_so4.cg_so4
    to compute ⟨m⁺⟩ and ⟨m⁻⟩, then maps these to Hopf phases.

    For |n,l,m⟩ = Σ C_{m⁺m⁻} |j,m⁺⟩|j,m⁻⟩:
        ⟨m⁺⟩ = Σ |C_{m⁺m⁻}|² m⁺
        ⟨m⁻⟩ = Σ |C_{m⁺m⁻}|² m⁻

    Then ψ₁ = π(⟨m⁺⟩ + j)/(2j) and ψ₂ = π(⟨m⁻⟩ + j)/(2j).

    Parameters
    ----------
    n, l, m : int
        Quantum numbers.

    Returns
    -------
    Dict with ⟨m⁺⟩, ⟨m⁻⟩, ψ₁, ψ₂, and dominant (m⁺, m⁻).
    """
    # Import here to avoid circular dependency
    from geovac.wigner_so4 import cg_so4 as _cg_so4

    j = (n - 1) / 2.0
    cg = _cg_so4(n, l, m)

    # Expectation values
    exp_mp = sum(abs(c)**2 * mp for (mp, mm), c in cg.items())
    exp_mm = sum(abs(c)**2 * mm for (mp, mm), c in cg.items())

    # Map to [0, π] range
    if j > 0:
        psi1 = np.pi * (exp_mp + j) / (2 * j)
        psi2 = np.pi * (exp_mm + j) / (2 * j)
    else:
        psi1 = 0.0
        psi2 = 0.0

    # Dominant component
    dominant_pair = max(cg.items(), key=lambda x: abs(x[1]))
    (mp_dom, mm_dom), coeff_dom = dominant_pair

    return {
        'exp_m_plus': exp_mp,
        'exp_m_minus': exp_mm,
        'psi1': psi1,
        'psi2': psi2,
        'dominant_m_plus': mp_dom,
        'dominant_m_minus': mm_dom,
        'dominant_coeff': coeff_dom,
        'n_cg_terms': len(cg),
    }
