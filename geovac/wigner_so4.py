"""
SO(4) Wigner D-matrix elements for the hydrogen atom.

This module implements the SO(4) rotation matrices that describe how hydrogen
eigenstates |n,l,m⟩ transform under SO(4) rotations of the momentum-space
three-sphere S³. In the bond-sphere picture (Paper 8), cross-atom matrix
elements are D-matrix elements evaluated at the bond angle γ.

Theory (Paper 8, Sec. III):
    SO(4) ≅ SU(2)⊗SU(2).  For the n-shell, j⁺ = j⁻ = j = (n-1)/2.
    The CG coefficients map |n,l,m⟩ → Σ |j,m⁺⟩|j,m⁻⟩.
    A rotation by angle γ about a single axis acts as
        D^(n)_{l'm', lm}(γ) = Σ C^{lm}_{m⁺m⁻} d^j_{m⁺'m⁺}(γ) d^j_{m⁻'m⁻}(γ) C^{l'm'}_{m⁺'m⁻'}

Functions:
    cg_so4(n, l, m)          -- CG coefficients: |n,l,m⟩ → |j,m⁺⟩⊗|j,m⁻⟩
    wigner_d_su2(j, mp, m, angle) -- SU(2) small d-matrix element d^j_{mp,m}(β)
    wigner_D_so4(n, lp, mp, l, m, gamma) -- SO(4) D-matrix element
    bond_angle(R, p0)        -- cos γ = (p₀² - p_R²)/(p₀² + p_R²)
    d_matrix_block(n, gamma) -- full n×n D-matrix block for the n-shell

References:
    Bander & Itzykson, Rev. Mod. Phys. 38, 330 (1966)
    Shibuya & Wulfman, Proc. R. Soc. A 286, 376 (1965)
    Paper 8: Bond Sphere Geometry (GeoVac, 2026)
"""

from __future__ import annotations

import numpy as np
from typing import Dict, Tuple, List
from functools import lru_cache


# ---------------------------------------------------------------------------
#  CG coefficients:  |n, l, m⟩  →  |j, m⁺⟩ ⊗ |j, m⁻⟩
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1024)
def cg_so4(n: int, l: int, m: int) -> Dict[Tuple[float, float], float]:
    """Return CG coefficients mapping |n,l,m⟩ to SU(2)⊗SU(2) basis.

    For the hydrogen n-shell, SO(4) ≅ SU(2)⊗SU(2) with j = (n-1)/2.
    The angular momentum coupling is:
        |n, l, m⟩ = Σ_{m⁺,m⁻} C^{l,m}_{m⁺,m⁻} |j, m⁺⟩ ⊗ |j, m⁻⟩

    where m⁺ + m⁻ = m, and the CG coefficient is ⟨j,m⁺; j,m⁻ | l,m⟩.

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
    dict
        {(m_plus, m_minus): coefficient} for non-zero coefficients.
        Keys are half-integer floats.
    """
    if n < 1:
        raise ValueError(f"n must be ≥ 1, got {n}")
    if l < 0 or l >= n:
        raise ValueError(f"l must satisfy 0 ≤ l < n={n}, got l={l}")
    if abs(m) > l:
        raise ValueError(f"|m| must be ≤ l={l}, got m={m}")

    j = (n - 1) / 2.0  # half-integer

    # The CG coefficient is ⟨j, m⁺; j, m⁻ | l, m⟩
    # with the constraint m⁺ + m⁻ = m
    # This couples two spin-j representations to total angular momentum l.
    result = {}
    # m⁺ ranges over -j, -j+1, ..., j
    steps = int(2 * j) + 1
    for i_plus in range(steps):
        mp = -j + i_plus
        mm = m - mp  # m⁻ from constraint m⁺ + m⁻ = m
        if abs(mm) > j:
            continue
        coeff = _clebsch_gordan(j, mp, j, mm, l, m)
        if abs(coeff) > 1e-15:
            result[(mp, mm)] = coeff

    return result


# ---------------------------------------------------------------------------
#  SU(2) Wigner small d-matrix
# ---------------------------------------------------------------------------

@lru_cache(maxsize=4096)
def wigner_d_su2(j: float, mp: float, m: float, angle: float) -> float:
    """Compute the SU(2) Wigner small d-matrix element d^j_{mp,m}(β).

    Uses the explicit formula (Varshalovich et al.):
        d^j_{mp,m}(β) = Σ_s  (-1)^{mp-m+s} √((j+mp)!(j-mp)!(j+m)!(j-m)!)
                         / (s!(j+m-s)!(j-mp-s)!(s+mp-m)!)
                         × (cos β/2)^{2j-2s+m-mp} (sin β/2)^{2s+mp-m}

    where the sum runs over all s giving non-negative factorials.
    Sign convention: (-1)^{m'-m+s} gives the standard Condon-Shortley
    phase, matching d^{1/2} = [[cos, -sin],[sin, cos]].

    Parameters
    ----------
    j : float
        Angular momentum quantum number (integer or half-integer ≥ 0).
    mp : float
        Row index (-j ≤ mp ≤ j).
    m : float
        Column index (-j ≤ m ≤ j).
    angle : float
        Rotation angle β in radians.

    Returns
    -------
    float
        The d-matrix element.
    """
    # Convert to ensure we handle half-integers properly
    two_j = int(round(2 * j))
    two_mp = int(round(2 * mp))
    two_m = int(round(2 * m))

    if abs(two_mp) > two_j or abs(two_m) > two_j:
        return 0.0

    cos_half = np.cos(angle / 2.0)
    sin_half = np.sin(angle / 2.0)

    # Integer-valued indices for factorial
    jpm = (two_j + two_m) // 2      # j + m
    jmm = (two_j - two_m) // 2      # j - m
    jpmp = (two_j + two_mp) // 2     # j + mp
    jmmp = (two_j - two_mp) // 2     # j - mp

    prefactor = np.sqrt(
        _factorial(jpm) * _factorial(jmm) *
        _factorial(jpmp) * _factorial(jmmp)
    )

    result = 0.0
    # s ranges: s ≥ 0, j+m-s ≥ 0, j-mp-s ≥ 0, s+mp-m ≥ 0
    s_min = max(0, (two_m - two_mp) // 2)
    s_max = min(jpm, jmmp)

    for s in range(s_min, s_max + 1):
        denom = (
            _factorial(s) *
            _factorial(jpm - s) *
            _factorial(jmmp - s) *
            _factorial(s + (two_mp - two_m) // 2)
        )
        # Exponents for cos and sin
        cos_exp = two_j - 2 * s + (two_m - two_mp) // 2
        sin_exp = 2 * s + (two_mp - two_m) // 2

        mp_minus_m = (two_mp - two_m) // 2  # integer: m' - m
        sign = (-1) ** (mp_minus_m + s)
        term = sign * prefactor / denom

        if cos_exp > 0:
            term *= cos_half ** cos_exp
        elif cos_exp < 0:
            # Should not happen with valid s range
            term = 0.0

        if sin_exp > 0:
            term *= sin_half ** sin_exp
        elif sin_exp < 0:
            term = 0.0

        result += term

    return float(result)


# ---------------------------------------------------------------------------
#  SO(4) Wigner D-matrix element
# ---------------------------------------------------------------------------

def wigner_D_so4(n: int, lp: int, mp: int,
                 l: int, m: int, gamma: float) -> float:
    """Compute SO(4) D-matrix element D^(n)_{l'm', lm}(γ).

    The bond rotation is generated by A (Runge-Lenz), where A = J⁺ - J⁻.
    Since J⁺ and J⁻ commute, exp(iγ A_y) = exp(iγ J⁺_y) exp(-iγ J⁻_y),
    so j⁺ and j⁻ rotate in OPPOSITE directions:

        D^(n)_{l'm', lm}(γ) = Σ_{m⁺,m⁻,m⁺',m⁻'}
            C^{l',m'}_{m⁺',m⁻'} d^j_{m⁺',m⁺}(γ) d^j_{m⁻',m⁻}(-γ) C^{l,m}_{m⁺,m⁻}

    where j = (n-1)/2 and C are CG coefficients from cg_so4().

    Parameters
    ----------
    n : int
        Principal quantum number.
    lp, mp : int
        Bra quantum numbers (l', m').
    l, m : int
        Ket quantum numbers (l, m).
    gamma : float
        Bond angle in radians.

    Returns
    -------
    float
        The D-matrix element.
    """
    if n < 1:
        raise ValueError(f"n must be ≥ 1, got {n}")

    j = (n - 1) / 2.0

    # Get CG decompositions
    cg_ket = cg_so4(n, l, m)     # |n,l,m⟩ → {(m⁺,m⁻): coeff}
    cg_bra = cg_so4(n, lp, mp)   # |n,l',m'⟩ → {(m⁺',m⁻'): coeff}

    result = 0.0
    for (mp_plus, mp_minus), c_bra in cg_bra.items():
        for (m_plus, m_minus), c_ket in cg_ket.items():
            d_plus = wigner_d_su2(j, mp_plus, m_plus, gamma)    # J⁺ rotates by +γ
            d_minus = wigner_d_su2(j, mp_minus, m_minus, -gamma) # J⁻ rotates by -γ
            result += c_bra * d_plus * d_minus * c_ket

    return float(result)


# ---------------------------------------------------------------------------
#  Bond angle
# ---------------------------------------------------------------------------

def bond_angle(R: float, p0: float) -> float:
    """Compute the bond angle γ from internuclear distance R and momentum scale p₀.

    The stereographic mapping gives:
        cos γ = (p₀² - p_R²) / (p₀² + p_R²)

    where p_R = 1/R (momentum conjugate to the bond length).

    At R → ∞: p_R → 0, cos γ → 1, γ → 0 (identity, atoms separate).
    At R → 0: p_R → ∞, cos γ → -1, γ → π (antipodal, united atom).

    Parameters
    ----------
    R : float
        Internuclear distance in bohr. Must be > 0.
    p0 : float
        Momentum scale p₀ = √(2|E|) for the energy shell.

    Returns
    -------
    float
        Bond angle γ in radians, range [0, π].
    """
    if R <= 0:
        raise ValueError(f"R must be > 0, got {R}")
    if p0 <= 0:
        raise ValueError(f"p0 must be > 0, got {p0}")

    p_R = 1.0 / R
    cos_gamma = (p0**2 - p_R**2) / (p0**2 + p_R**2)
    # Clamp for numerical safety
    cos_gamma = max(-1.0, min(1.0, cos_gamma))
    return float(np.arccos(cos_gamma))


# ---------------------------------------------------------------------------
#  Full D-matrix block for an n-shell
# ---------------------------------------------------------------------------

def d_matrix_block(n: int, gamma: float) -> np.ndarray:
    """Compute the full D-matrix block for the hydrogen n-shell.

    Returns an n² × n² matrix with rows/columns ordered by (l, m) with
    l = 0, 1, ..., n-1 and m = -l, -l+1, ..., l for each l.

    Parameters
    ----------
    n : int
        Principal quantum number (n ≥ 1).
    gamma : float
        Bond angle in radians.

    Returns
    -------
    np.ndarray
        Shape (n², n²). The D^(n)(γ) matrix.
    """
    if n < 1:
        raise ValueError(f"n must be ≥ 1, got {n}")

    # Build index list: (l, m) pairs
    states = []
    for l in range(n):
        for m_val in range(-l, l + 1):
            states.append((l, m_val))

    dim = len(states)  # should be n²
    assert dim == n * n, f"Expected {n**2} states, got {dim}"

    D = np.zeros((dim, dim))
    for i, (lp, mp) in enumerate(states):
        for k, (l, m_val) in enumerate(states):
            D[i, k] = wigner_D_so4(n, lp, mp, l, m_val, gamma)

    return D


# ---------------------------------------------------------------------------
#  Internal helpers
# ---------------------------------------------------------------------------

_FACTORIAL_CACHE: List[float] = [1.0]


def _factorial(n: int) -> float:
    """Cached factorial for small n."""
    if n < 0:
        raise ValueError(f"Factorial of negative number {n}")
    while len(_FACTORIAL_CACHE) <= n:
        _FACTORIAL_CACHE.append(_FACTORIAL_CACHE[-1] * len(_FACTORIAL_CACHE))
    return _FACTORIAL_CACHE[n]


@lru_cache(maxsize=4096)
def _clebsch_gordan(j1: float, m1: float, j2: float, m2: float,
                    J: int, M: int) -> float:
    """Compute Clebsch-Gordan coefficient ⟨j1,m1; j2,m2 | J, M⟩.

    Uses the Racah formula. j1, j2 can be half-integer (as floats).

    Parameters
    ----------
    j1, m1 : float
        First angular momentum and projection.
    j2, m2 : float
        Second angular momentum and projection.
    J : int
        Total angular momentum (integer for j1=j2 coupling).
    M : int
        Total projection (integer).

    Returns
    -------
    float
        The CG coefficient.
    """
    # Selection rules
    if abs(m1 + m2 - M) > 1e-10:
        return 0.0
    if J < abs(j1 - j2) - 1e-10 or J > j1 + j2 + 1e-10:
        return 0.0
    if abs(m1) > j1 + 1e-10 or abs(m2) > j2 + 1e-10 or abs(M) > J + 1e-10:
        return 0.0

    # Convert to twice-values for integer arithmetic
    tj1 = int(round(2 * j1))
    tm1 = int(round(2 * m1))
    tj2 = int(round(2 * j2))
    tm2 = int(round(2 * m2))
    tJ = int(round(2 * J))
    tM = int(round(2 * M))

    # Triangle coefficient
    # Δ(j1, j2, J) = √((j1+j2-J)!(j1-j2+J)!(-j1+j2+J)! / (j1+j2+J+1)!)
    a = (tj1 + tj2 - tJ) // 2
    b = (tj1 - tj2 + tJ) // 2
    c = (-tj1 + tj2 + tJ) // 2
    d = (tj1 + tj2 + tJ) // 2 + 1

    if a < 0 or b < 0 or c < 0:
        return 0.0

    # Prefactor
    pf1 = (tj1 + tm1) // 2   # j1 + m1
    pf2 = (tj1 - tm1) // 2   # j1 - m1
    pf3 = (tj2 + tm2) // 2   # j2 + m2
    pf4 = (tj2 - tm2) // 2   # j2 - m2
    pf5 = (tJ + tM) // 2     # J + M
    pf6 = (tJ - tM) // 2     # J - M

    if pf1 < 0 or pf2 < 0 or pf3 < 0 or pf4 < 0 or pf5 < 0 or pf6 < 0:
        return 0.0

    prefactor = np.sqrt(
        (tJ + 1) *
        _factorial(a) * _factorial(b) * _factorial(c) / _factorial(d) *
        _factorial(pf1) * _factorial(pf2) *
        _factorial(pf3) * _factorial(pf4) *
        _factorial(pf5) * _factorial(pf6)
    )

    # Sum over s
    s_min = max(0, (tj2 - tJ + tm1) // 2, (tj1 - tJ - tm2) // 2)
    s_min = max(s_min, 0)
    s_max = min(a, pf1, pf4)

    total = 0.0
    for s in range(s_min, s_max + 1):
        denom = (
            _factorial(s) *
            _factorial(a - s) *
            _factorial(pf1 - s) *
            _factorial(pf4 - s) *
            _factorial(c - pf1 + s) *
            _factorial(b - pf4 + s)
        )
        # Check all factorial args are non-negative
        arg_c = c - pf1 + s
        arg_b = b - pf4 + s
        if arg_c < 0 or arg_b < 0:
            continue
        total += (-1) ** s / denom

    return float(prefactor * total)
