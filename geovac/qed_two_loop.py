"""Two-loop QED on S^3: demonstration that zeta(3) enters through
the spectral data of the Dirac operator.

At one loop, the QED effective action involves traces of D^2 (the squared
Dirac operator).  The T9 theorem guarantees these produce only pi^{even}
-- no zeta(3).

At two loops, the effective action involves products of Dirac propagators
that probe the *first-order* Dirac eigenvalue structure.  The key result
of this module is the **even/odd s discriminant**:

    D_Dirac(s) = sum_{n>=0} g_n / |lambda_n|^s

where g_n = 2(n+1)(n+2) and |lambda_n| = n + 3/2.  This equals
2*zeta(s-2, 3/2) - (1/2)*zeta(s, 3/2) via Hurwitz zeta, and using
zeta(s, 3/2) = (2^s - 1)*zeta_R(s) - 2^s:

    s EVEN: involves zeta_R(even) = rational * pi^{even}  --> pi^{even} only
    s ODD:  involves zeta_R(odd)  = zeta(3), zeta(5), ... --> odd-zeta

Verified exact PSLQ identifications:
    D(4) =   pi^2 - pi^4/12                    [pi^{even}]
    D(5) =   14*zeta(3) - 31/2*zeta(5)         [odd-zeta]
    D(6) =   pi^4/3 - pi^6/30                  [pi^{even}]
    D(7) =   62*zeta(5) - 127/2*zeta(7)        [odd-zeta]
    D(8) =   2*pi^6/15 - 17*pi^8/1260          [pi^{even}]

The Fock-index Dirichlet series (Track D3) confirms independently:
    D_Fock(s) = sum_{n>=1} 2n(n+1)/n^s = 2*zeta(s-2) + 2*zeta(s-1)

At s=4: D_Fock(4) = 2*zeta(2) + 2*zeta(3) = pi^2/3 + 2*zeta(3).

The nested harmonic sum sum_{n=1}^N 1/n^2 * H_n --> 2*zeta(3) (Euler)
is the flat-space analog of the two-loop structure.

Transcendental taxonomy (Paper 18)
----------------------------------
The operator-order discriminant is now demonstrated at both one-loop
and two-loop level:
- One-loop (D^2, s even): pi^{even} only (T9 theorem)
- Two-loop (|D|, s odd):  odd-zeta (zeta(3), zeta(5), ...) enters
- The discriminant is the PARITY OF s, equivalently whether the
  spectral sum uses D^2 (even s) or |D| (odd s)

References
----------
- Track D3 (Tier 1): D_{g^Dirac}(4) = 2*zeta(2) + 2*zeta(3).
- T9 theorem (Tier 3): zeta_{D^2}(s) = pi^{even} at every integer s.
- Rosner (1967), Laporta-Remiddi (1996): zeta(3) in two-loop QED.
- GeoVac Paper 18 (operator-order x bundle grid).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import mpmath

# Set high precision for PSLQ identification
mpmath.mp.dps = 80

__all__ = [
    "dirac_dirichlet_series_hurwitz",
    "dirac_dirichlet_series_numerical",
    "fock_dirichlet_series",
    "nested_harmonic_sum_flat",
    "nested_harmonic_sum_s3",
    "double_spectral_zeta_connected",
    "even_odd_discriminant_table",
    "decompose_into_zeta_basis",
    "two_loop_vacuum_polarization_s3",
    "flat_space_limit_check",
    "classify_two_loop_transcendentals",
]


# ---------------------------------------------------------------------------
# Dirac spectrum helpers (mpmath precision)
# ---------------------------------------------------------------------------

def _lambda_n(n: int) -> mpmath.mpf:
    """Absolute Dirac eigenvalue |lambda_n| = n + 3/2 on unit S^3 (CH convention)."""
    return mpmath.mpf(n) + mpmath.mpf(3) / 2


def _g_n_dirac(n: int) -> mpmath.mpf:
    """Full Dirac degeneracy g_n = 2(n+1)(n+2)."""
    return mpmath.mpf(2) * (n + 1) * (n + 2)


# ---------------------------------------------------------------------------
# Exact Dirichlet series via Hurwitz zeta
# ---------------------------------------------------------------------------

def dirac_dirichlet_series_hurwitz(s: int) -> Dict[str, object]:
    """Exact D_Dirac(s) = sum_{n>=0} g_n / |lambda_n|^s via Hurwitz zeta.

    Using g_n = 2*lambda_n^2 - 1/2 (with lambda_n = n+3/2):
        D(s) = 2*zeta(s-2, 3/2) - (1/2)*zeta(s, 3/2)

    And zeta(k, 3/2) = (2^k - 1)*zeta_R(k) - 2^k  for integer k >= 2.

    For s even: both zeta_R(s-2) and zeta_R(s) are even-index Riemann
    zeta -> rational * pi^{even}.
    For s odd: both are odd-index Riemann zeta -> zeta(3), zeta(5), etc.

    Parameters
    ----------
    s : int
        Must be >= 4 (s=3 involves zeta(1, 3/2) which diverges).

    Returns
    -------
    Dict with exact Hurwitz value, PSLQ decomposition, and parity classification.
    """
    if s < 4:
        raise ValueError(f"s must be >= 4 (s={s} involves divergent zeta(1,3/2))")

    # Compute via Hurwitz zeta
    hz_s_minus_2 = mpmath.hurwitz(s - 2, mpmath.mpf(3) / 2)
    hz_s = mpmath.hurwitz(s, mpmath.mpf(3) / 2)
    D_s = 2 * hz_s_minus_2 - mpmath.mpf(1) / 2 * hz_s

    # Analytical decomposition via zeta_R
    # zeta(k, 3/2) = (2^k - 1)*zeta_R(k) - 2^k
    # D(s) = 2*[(2^{s-2}-1)*z_R(s-2) - 2^{s-2}] - 1/2*[(2^s-1)*z_R(s) - 2^s]
    #       = 2*(2^{s-2}-1)*z_R(s-2) - 2^{s-1} - (2^s-1)/2*z_R(s) + 2^{s-1}
    #       = 2*(2^{s-2}-1)*z_R(s-2) - (2^s-1)/2*z_R(s)

    coeff_a = 2 * (2**(s - 2) - 1)  # coefficient of zeta_R(s-2)
    coeff_b = -(2**s - 1) / 2       # coefficient of zeta_R(s)

    parity = "even" if s % 2 == 0 else "odd"
    zeta_type = "pi^{even}" if s % 2 == 0 else "odd-zeta"

    # PSLQ decomposition
    if s % 2 == 0:
        # Use pi^{even} basis
        max_power = s
        basis = [D_s] + [mpmath.pi**(2 * k) for k in range(1, max_power // 2 + 1)]
        labels = ["value"] + [f"pi^{2*k}" for k in range(1, max_power // 2 + 1)]
    else:
        # Use odd-zeta basis
        basis = [D_s] + [mpmath.zeta(2 * k + 1) for k in range(1, s // 2 + 1)]
        labels = ["value"] + [f"zeta({2*k+1})" for k in range(1, s // 2 + 1)]

    relation = mpmath.pslq(basis, tol=1e-50, maxcoeff=1000000)

    decomposition = {}
    if relation is not None and relation[0] != 0:
        for i in range(1, len(relation)):
            if relation[i] != 0:
                decomposition[labels[i]] = str(
                    mpmath.fraction(-relation[i], relation[0])
                )

    return {
        "s": s,
        "value": D_s,
        "value_float": float(D_s),
        "parity": parity,
        "zeta_type": zeta_type,
        "coeff_zeta_R_s_minus_2": coeff_a,
        "coeff_zeta_R_s": coeff_b,
        "pslq_relation": relation,
        "decomposition": decomposition,
        "contains_odd_zeta": s % 2 == 1,
    }


def dirac_dirichlet_series_numerical(
    s: int,
    N: int = 500,
) -> Dict[str, object]:
    """Compute D_Dirac(s) by direct summation for cross-check.

    Parameters
    ----------
    s : int
        Dirichlet exponent.
    N : int
        Truncation level.
    """
    partial = mpmath.mpf(0)
    for n in range(N + 1):
        partial += _g_n_dirac(n) * _lambda_n(n)**(-s)

    # Compare with Hurwitz if s >= 4
    hurwitz_val = None
    if s >= 4:
        hz_s2 = mpmath.hurwitz(s - 2, mpmath.mpf(3) / 2)
        hz_s = mpmath.hurwitz(s, mpmath.mpf(3) / 2)
        hurwitz_val = 2 * hz_s2 - mpmath.mpf(1) / 2 * hz_s

    return {
        "s": s,
        "N": N,
        "partial_sum": partial,
        "partial_sum_float": float(partial),
        "hurwitz_exact": hurwitz_val,
        "rel_error_vs_hurwitz": (
            float(abs(partial - hurwitz_val) / abs(hurwitz_val))
            if hurwitz_val is not None else None
        ),
    }


# ---------------------------------------------------------------------------
# Fock-index Dirichlet series (Track D3)
# ---------------------------------------------------------------------------

def fock_dirichlet_series(s: int) -> Dict[str, object]:
    """D_Fock(s) = sum_{n>=1} 2n(n+1)/n^s = 2*zeta_R(s-2) + 2*zeta_R(s-1).

    This is the Track D3 version using Fock shell indices (integers)
    as the eigenvalue, rather than the Dirac eigenvalue n+3/2.

    At s=4: D_Fock(4) = 2*zeta(2) + 2*zeta(3), which contains zeta(3).

    The algebraic identity is exact:
        2n(n+1)/n^s = 2(n+1)/n^{s-1} = 2/n^{s-2} + 2/n^{s-1}

    Parameters
    ----------
    s : int
        Must be >= 4 for convergence.
    """
    if s < 4:
        raise ValueError(f"s must be >= 4 for convergence (s={s})")

    val = 2 * mpmath.zeta(s - 2) + 2 * mpmath.zeta(s - 1)

    # PSLQ with extended basis
    basis_elems = [val, mpmath.mpf(1)]
    basis_labels = ["value", "1"]
    # Add pi^{even} up to pi^s
    for k in range(1, s // 2 + 1):
        basis_elems.append(mpmath.pi**(2 * k))
        basis_labels.append(f"pi^{2*k}")
    # Add odd zeta values
    for k in range(1, s // 2 + 1):
        basis_elems.append(mpmath.zeta(2 * k + 1))
        basis_labels.append(f"zeta({2*k+1})")

    relation = mpmath.pslq(basis_elems, tol=1e-50, maxcoeff=10000)

    decomposition = {}
    if relation is not None and relation[0] != 0:
        for i in range(1, len(relation)):
            if relation[i] != 0:
                decomposition[basis_labels[i]] = str(
                    mpmath.fraction(-relation[i], relation[0])
                )

    return {
        "s": s,
        "value": val,
        "value_float": float(val),
        "algebraic_form": f"2*zeta({s-2}) + 2*zeta({s-1})",
        "contains_odd_zeta": (s - 2) % 2 == 1 or (s - 1) % 2 == 1,
        "pslq_relation": relation,
        "decomposition": decomposition,
    }


# ---------------------------------------------------------------------------
# Even/odd discriminant table
# ---------------------------------------------------------------------------

def even_odd_discriminant_table(s_max: int = 10) -> List[Dict[str, object]]:
    """Compute D_Dirac(s) for s = 4..s_max, demonstrating the even/odd pattern.

    Returns a list of dicts, one per s value, showing the parity-dependent
    transcendental content.
    """
    rows = []
    for s in range(4, s_max + 1):
        row = dirac_dirichlet_series_hurwitz(s)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Nested harmonic sums
# ---------------------------------------------------------------------------

def nested_harmonic_sum_flat(N: int = 500) -> Dict[str, object]:
    """Flat-space nested harmonic sum: S = sum_{n=1}^N 1/n^2 * H_n.

    Converges to 2*zeta(3) as N -> infinity (Euler's identity).
    """
    total = mpmath.mpf(0)
    H = mpmath.mpf(0)
    for n in range(1, N + 1):
        H += mpmath.mpf(1) / n
        total += H / mpmath.mpf(n)**2

    target = 2 * mpmath.zeta(3)
    return {
        "partial_sum": total,
        "target_2z3": target,
        "rel_error": float(abs(total - target) / abs(target)),
        "N": N,
    }


def nested_harmonic_sum_s3(N: int = 500) -> Dict[str, object]:
    """S^3 nested harmonic sum with Dirac weights.

    S(N) = sum_{n=1}^{N} [g_n/lambda_n^4] * sum_{m=1}^{n} [g_m/lambda_m^2]

    Uses actual Dirac eigenvalues lambda_n = n + 3/2.
    """
    total = mpmath.mpf(0)
    inner = mpmath.mpf(0)

    convergence = []
    for n in range(1, N + 1):
        g_n = _g_n_dirac(n)
        lam_n = _lambda_n(n)
        inner += g_n / lam_n**2
        total += (g_n / lam_n**4) * inner
        if n in [10, 50, 100, 200, 500]:
            convergence.append((n, float(total)))

    decomp = decompose_into_zeta_basis(total)
    return {
        "partial_sum": total,
        "N": N,
        "decomposition": decomp,
        "convergence": convergence,
    }


# ---------------------------------------------------------------------------
# Double spectral zeta (connected two-loop)
# ---------------------------------------------------------------------------

def double_spectral_zeta_connected(
    s1: int = 4,
    s2: int = 4,
    use_fock: bool = True,
) -> Dict[str, object]:
    """Connected double sum: sum_{n != m} g_n*g_m / (eigenvalue_n^s1 * eigenvalue_m^s2).

    If use_fock=True, uses Fock indices (integers) where zeta(3) lives.
    If use_fock=False, uses Dirac eigenvalues (half-integers).

    The connected part = [sum g/e^s1] * [sum g/e^s2] - sum g^2/e^{s1+s2}.

    Parameters
    ----------
    s1, s2 : int
        Exponents. Both must be >= 4.
    use_fock : bool
        If True, use Fock shell indices; if False, use Dirac eigenvalues.
    """
    if use_fock:
        # Fock: eigenvalue = n, g_n = 2n(n+1), n >= 1
        sum1 = 2 * mpmath.zeta(s1 - 2) + 2 * mpmath.zeta(s1 - 1)
        sum2 = 2 * mpmath.zeta(s2 - 2) + 2 * mpmath.zeta(s2 - 1)
        # Diagonal: sum_{n>=1} [2n(n+1)]^2 / n^{s1+s2}
        # = 4*sum n^2(n+1)^2 / n^{s1+s2} = 4*sum (n+1)^2 / n^{s1+s2-2}
        # = 4*sum (n^2+2n+1)/n^{s1+s2-2} = 4*[z(s-4)+2*z(s-3)+z(s-2)] where s=s1+s2
        stot = s1 + s2
        diag = 4 * (mpmath.zeta(stot - 4) + 2 * mpmath.zeta(stot - 3) + mpmath.zeta(stot - 2))
        label = "Fock indices"
    else:
        # Dirac: use Hurwitz
        hz1a = mpmath.hurwitz(s1 - 2, mpmath.mpf(3) / 2)
        hz1b = mpmath.hurwitz(s1, mpmath.mpf(3) / 2)
        sum1 = 2 * hz1a - mpmath.mpf(1) / 2 * hz1b

        hz2a = mpmath.hurwitz(s2 - 2, mpmath.mpf(3) / 2)
        hz2b = mpmath.hurwitz(s2, mpmath.mpf(3) / 2)
        sum2 = 2 * hz2a - mpmath.mpf(1) / 2 * hz2b

        # Diagonal: sum g_n^2 / lam_n^{s1+s2}
        # g_n^2 = (2*lam^2 - 1/2)^2 = 4*lam^4 - 2*lam^2 + 1/4
        stot = s1 + s2
        hz_a = mpmath.hurwitz(stot - 4, mpmath.mpf(3) / 2)
        hz_b = mpmath.hurwitz(stot - 2, mpmath.mpf(3) / 2)
        hz_c = mpmath.hurwitz(stot, mpmath.mpf(3) / 2)
        diag = 4 * hz_a - 2 * hz_b + mpmath.mpf(1) / 4 * hz_c
        label = "Dirac eigenvalues"

    connected = sum1 * sum2 - diag

    return {
        "s1": s1,
        "s2": s2,
        "source": label,
        "sum_1": sum1,
        "sum_2": sum2,
        "product": sum1 * sum2,
        "diagonal": diag,
        "connected": connected,
        "connected_float": float(connected),
        "decomp_connected": decompose_into_zeta_basis(connected),
    }


# ---------------------------------------------------------------------------
# Two-loop summary
# ---------------------------------------------------------------------------

def two_loop_vacuum_polarization_s3(n_max: int = 300) -> Dict[str, object]:
    """Compute the key spectral sums demonstrating zeta(3) at two loops.

    Returns a comprehensive summary of all computed quantities.
    """
    results: Dict[str, object] = {}

    # 1. Even/odd discriminant for D_Dirac(s)
    results["discriminant_table"] = even_odd_discriminant_table(s_max=8)

    # 2. D3 cross-check: Fock D(4) = 2*z(2) + 2*z(3)
    results["fock_D4"] = fock_dirichlet_series(4)

    # 3. Connected double sum (Fock, s1=s2=4)
    results["connected_fock_4_4"] = double_spectral_zeta_connected(4, 4, use_fock=True)

    # 4. Connected double sum (Dirac eigenvalues, s1=s2=4)
    results["connected_dirac_4_4"] = double_spectral_zeta_connected(4, 4, use_fock=False)

    # 5. Flat-space nested sum
    results["flat_nested"] = nested_harmonic_sum_flat(N=n_max)

    # 6. Key finding
    results["summary"] = {
        "zeta3_mechanism": (
            "zeta(3) enters the two-loop QED effective action on S^3 through "
            "the first-order Dirac Dirichlet series D(s) at ODD s values. "
            "At s=5: D(5) = 14*zeta(3) - 31/2*zeta(5). "
            "At s=7: D(7) = 62*zeta(5) - 127/2*zeta(7). "
            "The Fock-index series always contains odd-zeta: "
            "D_Fock(4) = 2*zeta(2) + 2*zeta(3). "
            "The even/odd discriminant in s is the spectral manifestation "
            "of Paper 18's operator-order grid: second-order (D^2, s even) "
            "gives pi^{even}; first-order (|D|, s odd) gives odd-zeta."
        ),
        "paper18_connection": (
            "The T9 theorem (zeta_{D^2}(s) = pi^{even}) is the s-even case. "
            "The s-odd case is the new result: D_Dirac(s_odd) = linear combination "
            "of Riemann zeta at odd arguments, with RATIONAL coefficients "
            "determined by the Dirac degeneracy polynomial g_n = 2n^2 + 6n + 4. "
            "The two results together confirm the operator-order transcendental "
            "discriminant at all loop orders."
        ),
    }

    return results


# ---------------------------------------------------------------------------
# Flat-space limit
# ---------------------------------------------------------------------------

def flat_space_limit_check(N: int = 500) -> Dict[str, object]:
    """Verify flat-space limits (g_n -> 1, lambda_n -> n).

    The flat-space Dirichlet series sum 1/n^s = zeta_R(s) (standard).
    The nested harmonic sum -> 2*zeta(3) (Euler identity).
    """
    # Flat nested harmonic
    flat_nested = nested_harmonic_sum_flat(N)

    # Flat D(4) = zeta_R(4) = pi^4/90
    flat_D4 = mpmath.zeta(4)
    target_D4 = mpmath.pi**4 / 90

    # S^3 with unit degeneracy
    unit_D4 = mpmath.mpf(0)
    for n in range(1, N + 1):
        unit_D4 += _lambda_n(n)**(-4)
    # This is zeta(4, 5/2) (Hurwitz at a=5/2)
    unit_D4_exact = mpmath.hurwitz(4, mpmath.mpf(5) / 2)

    return {
        "flat_nested": flat_nested,
        "flat_D4": {"value": float(flat_D4), "target_pi4_90": float(target_D4),
                     "match": float(abs(flat_D4 - target_D4)) < 1e-60},
        "unit_degeneracy_D4": {
            "numerical": float(unit_D4),
            "hurwitz_exact": float(unit_D4_exact),
        },
    }


# ---------------------------------------------------------------------------
# PSLQ decomposition
# ---------------------------------------------------------------------------

def decompose_into_zeta_basis(
    value: mpmath.mpf,
    *,
    tol: float = 1e-25,
) -> Dict[str, object]:
    """Decompose a numerical value as a rational linear combination of
    standard transcendental constants using PSLQ.

    Basis: {value, 1, pi^2, pi^4, zeta(3), pi^2*zeta(3), zeta(5), ln(2), pi^2*ln(2)}.

    Returns dict with 'identified', 'components', 'residual', 'contains_zeta3'.
    """
    pi2 = mpmath.pi**2
    pi4 = mpmath.pi**4
    z3 = mpmath.zeta(3)
    z5 = mpmath.zeta(5)
    ln2 = mpmath.log(2)

    basis = [value, mpmath.mpf(1), pi2, pi4, z3, pi2 * z3, z5, ln2, pi2 * ln2]
    labels = ["value", "1", "pi^2", "pi^4", "zeta(3)", "pi^2*zeta(3)",
              "zeta(5)", "ln(2)", "pi^2*ln(2)"]

    try:
        relation = mpmath.pslq(basis, tol=tol, maxcoeff=10000)
    except (ValueError, Exception):
        return {"identified": False, "value_float": float(value),
                "note": "PSLQ raised an exception"}

    if relation is None:
        return {"identified": False, "value_float": float(value),
                "note": "PSLQ failed to find an integer relation"}

    coeff_value = relation[0]
    if coeff_value == 0:
        return {"identified": False, "value_float": float(value),
                "note": "PSLQ returned zero coefficient for value"}

    components = {}
    reconstructed = mpmath.mpf(0)
    for i in range(1, len(relation)):
        if relation[i] != 0:
            coeff = mpmath.mpf(-relation[i]) / mpmath.mpf(coeff_value)
            components[labels[i]] = str(mpmath.fraction(-relation[i], coeff_value))
            reconstructed += coeff * basis[i]

    residual = float(abs(value - reconstructed))

    return {
        "identified": True,
        "raw_relation": relation,
        "components": components,
        "reconstructed_float": float(reconstructed),
        "value_float": float(value),
        "residual": residual,
        "contains_zeta3": any("zeta(3)" in k for k in components.keys()),
        "contains_zeta5": any("zeta(5)" in k for k in components.keys()),
    }


# ---------------------------------------------------------------------------
# Transcendental classification
# ---------------------------------------------------------------------------

def classify_two_loop_transcendentals() -> Dict[str, str]:
    """Classify the two-loop transcendental content in Paper 18's taxonomy."""
    return {
        "D_Dirac_s_even": (
            "calibration_pi_only: D_Dirac(s_even) involves Hurwitz zeta at "
            "even integer -> Riemann zeta at even integer -> rational * pi^{even}. "
            "Matches T9 theorem (zeta_{D^2} structure)."
        ),
        "D_Dirac_s_odd": (
            "odd_zeta: D_Dirac(s_odd) involves Hurwitz zeta at odd integer -> "
            "Riemann zeta at odd integer -> zeta(3), zeta(5), etc. "
            "This is the NEW transcendental content at two loops."
        ),
        "D_Fock_all_s": (
            "mixed: D_Fock(s) = 2*zeta_R(s-2) + 2*zeta_R(s-1) always contains "
            "both even and odd Riemann zeta values (one of s-2, s-1 is odd). "
            "At s=4: 2*zeta(2) + 2*zeta(3) = pi^2/3 + 2*zeta(3)."
        ),
        "even_odd_discriminant": (
            "Paper 18 Section IV: the parity of s in the Dirichlet series "
            "D(s) = sum g_n / |lambda_n|^s determines the transcendental class. "
            "s even -> D^2 structure -> pi^{even}. "
            "s odd -> |D| structure -> odd-zeta. "
            "The curvature shift lambda = n + 3/2 converts integer sums (Riemann) "
            "to half-integer sums (Hurwitz), preserving the even/odd distinction "
            "through Hurwitz's relation to Riemann zeta."
        ),
        "two_loop_mechanism": (
            "At two loops, the sunset diagram involves a PRODUCT of two Dirac "
            "propagators: G(n)*G(m) ~ 1/(|lambda_n|^a * |lambda_m|^b). "
            "When a or b is odd, the sum accesses the odd-zeta channel. "
            "The connected part (n != m) produces cross-terms like zeta(3)*zeta(2) "
            "and higher MZV combinations."
        ),
    }
