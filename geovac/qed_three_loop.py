"""Three-loop QED on S^3: probing depth-3 multiple zeta values.

At one loop, QED on S^3 produces only pi^{even} (T9 theorem, squared Dirac).
At two loops, the first-order Dirac structure exposes zeta(3), zeta(5), and
the vertex topology introduces Catalan G = beta(2) and Dirichlet beta(4).

The depth-k prediction: at k loops, expect depth-k multiple zeta values (MZVs).
At three loops we test whether depth-3 objects like zeta(3,1,1), alternating
Euler sums, or products zeta(3)^2 appear.

Three-loop self-energy topology
-------------------------------
The simplest three-loop self-energy on S^3 is the "iterated sunset":
three internal electron lines connected by two internal photon lines
in a chain topology:

    p --> [n1] --q1--> [n2] --q2--> [n3] --> p

Each vertex satisfies the SO(4) selection rule:
    V1: (n1, n2, q1) — triangle + parity
    V2: (n2, n3, q2) — triangle + parity

The three-loop sum is:

    Sigma^(3) = e^6 * sum_{n1,n2,n3,q1,q2}
        W(n1,n2,q1) * W(n2,n3,q2) * g_{n1} * g_{n2} * g_{n3} * d_{q1}^T * d_{q2}^T
        / (lambda_{n1}^a * lambda_{n2}^b * lambda_{n3}^c * mu_{q1}^p * mu_{q2}^p)

where W is the SO(4) channel count (0, 1, or 2) from qed_vertex.py.

Transcendental taxonomy (Paper 18)
----------------------------------
- Three-loop UNRESTRICTED (all n1,n2,n3): factorizes as D(a)*D(b)*D(c),
  product of three single Dirichlet series. At even exponents, purely pi^{even}.
- Three-loop with vertex restrictions: the parity constraints at each vertex
  create nested parity filtrations that may expose depth-3 multiple Hurwitz zeta.
- Three-loop with CG weights: the SO(4) channel count W creates nested
  min-weighted sums, potentially yielding depth-3 multiple zeta values.

Computational complexity
------------------------
The five-fold sum (n1, n2, n3, q1, q2) is O(N^5). At n_max=15 this is
~750,000 terms (with selection rules cutting ~90%). We use mpmath at 80+ dps
for PSLQ identification.

References
----------
- Laporta, Phys. Lett. B 549 (2002) 115-122 [three-loop QED vertex].
- Bailey & Broadhurst, Math. Comp. 70 (2001) 1719-1736 [MZVs in QFT].
- GeoVac Paper 18 (exchange constant taxonomy).
- GeoVac qed_vertex.py (two-loop infrastructure, SO(4) CG weights).
"""

from __future__ import annotations

import time
from typing import Dict, List, Optional, Tuple

import mpmath

# High precision for PSLQ
mpmath.mp.dps = 80

__all__ = [
    "three_loop_sunset_s3",
    "three_loop_unrestricted",
    "three_loop_vertex_restricted",
    "three_loop_cg_weighted",
    "three_loop_convergence_table",
    "decompose_three_loop",
    "classify_three_loop_transcendentals",
    "three_loop_factorized",
    "three_loop_factorized_convergence",
    "decompose_three_loop_mzv",
    "three_loop_euler_maclaurin_tail",
]


# ---------------------------------------------------------------------------
# Spectrum helpers (mpmath precision, duplicated to avoid cross-module issues)
# ---------------------------------------------------------------------------

def _lambda_n(n: int) -> mpmath.mpf:
    """Absolute Dirac eigenvalue |lambda_n| = n + 3/2 (CH convention)."""
    return mpmath.mpf(n) + mpmath.mpf(3) / 2


def _g_n_dirac(n: int) -> mpmath.mpf:
    """Full Dirac degeneracy g_n = 2(n+1)(n+2)."""
    return mpmath.mpf(2) * (n + 1) * (n + 2)


def _mu_q(q: int) -> mpmath.mpf:
    """Hodge-1 eigenvalue mu_q = q(q+2), q >= 1."""
    return mpmath.mpf(q) * (q + 2)


def _d_q_transverse(q: int) -> mpmath.mpf:
    """Transverse photon degeneracy d_q^T = q(q+2)."""
    return mpmath.mpf(q) * (q + 2)


def _vertex_allowed(n1: int, n2: int, n_gamma: int) -> bool:
    """SO(4) vertex selection rule: triangle + parity."""
    if n_gamma < 1:
        return False
    if n_gamma < abs(n1 - n2):
        return False
    if n_gamma > n1 + n2:
        return False
    if (n1 + n2 + n_gamma) % 2 == 0:
        return False
    return True


def _so4_channel_count(n1: int, n2: int, n_gamma: int) -> int:
    """SO(4) CG channel count (0, 1, or 2). Simplified from qed_vertex.py."""
    if not _vertex_allowed(n1, n2, n_gamma):
        return 0
    from fractions import Fraction
    j1_L = Fraction(n1 + 1, 2)
    j1_R = Fraction(n1, 2)
    j2_L = Fraction(n2, 2)
    j2_R = Fraction(n2 + 1, 2)
    q = n_gamma
    count = 0
    # Component A: ((q+1)/2, (q-1)/2)
    jg_L_A = Fraction(q + 1, 2)
    jg_R_A = Fraction(q - 1, 2)
    if (jg_R_A >= 0
            and abs(j1_L - jg_L_A) <= j2_L <= j1_L + jg_L_A
            and abs(j1_R - jg_R_A) <= j2_R <= j1_R + jg_R_A):
        count += 1
    # Component B: ((q-1)/2, (q+1)/2)
    jg_L_B = Fraction(q - 1, 2)
    jg_R_B = Fraction(q + 1, 2)
    if (jg_L_B >= 0
            and abs(j1_L - jg_L_B) <= j2_L <= j1_L + jg_L_B
            and abs(j1_R - jg_R_B) <= j2_R <= j1_R + jg_R_B):
        count += 1
    return count


# ---------------------------------------------------------------------------
# Dirac Dirichlet series via Hurwitz (for unrestricted comparison)
# ---------------------------------------------------------------------------

def _dirac_D(s: int) -> mpmath.mpf:
    """D_Dirac(s) via Hurwitz zeta."""
    if s < 4:
        raise ValueError(f"Need s >= 4 (got {s})")
    hz_s2 = mpmath.hurwitz(s - 2, mpmath.mpf(3) / 2)
    hz_s = mpmath.hurwitz(s, mpmath.mpf(3) / 2)
    return 2 * hz_s2 - mpmath.mpf(1) / 2 * hz_s


# ---------------------------------------------------------------------------
# Three-loop sums
# ---------------------------------------------------------------------------

def three_loop_sunset_s3(
    n_max: int,
    electron_exponent: int = 4,
    photon_exponent: int = 1,
    use_cg_weights: bool = True,
) -> Dict[str, object]:
    """Three-loop iterated sunset on S^3.

    Topology: three electron lines (n1, n2, n3) connected by two photon
    lines (q1, q2) in a chain:
        [n1] --q1-- [n2] --q2-- [n3]

    Vertex 1: (n1, n2, q1) with SO(4) selection rule
    Vertex 2: (n2, n3, q2) with SO(4) selection rule

    The sum is:
        S = sum W1*W2 * g1*g2*g3 * d_{q1}^T * d_{q2}^T
            / (lam1^a * lam2^a * lam3^a * mu_{q1}^p * mu_{q2}^p)

    where W1 = W(n1,n2,q1), W2 = W(n2,n3,q2) are CG channel counts
    (set to 1 if use_cg_weights=False for uniform vertex restriction).

    Parameters
    ----------
    n_max : int
        Truncation for electron levels (0..n_max) and photon levels (1..2*n_max).
    electron_exponent : int
        Exponent on each electron propagator |lambda_n|^a. Default 4 (s=2 per line).
    photon_exponent : int
        Exponent on each photon propagator mu_q^p. Default 1.
    use_cg_weights : bool
        If True, use SO(4) CG channel count W(n1,n2,q). If False, W=1 for
        all allowed triples (uniform vertex restriction).

    Returns
    -------
    Dict with the sum value, triple counts, and convergence data.
    """
    a = electron_exponent
    p = photon_exponent

    total = mpmath.mpf(0)
    n_quintuples = 0

    # Precompute electron data
    g_arr = [_g_n_dirac(n) for n in range(n_max + 1)]
    lam_arr = [_lambda_n(n) for n in range(n_max + 1)]
    lam_pow_arr = [lam_arr[n] ** a for n in range(n_max + 1)]

    q_max = 2 * n_max

    # Precompute photon data
    d_q_arr = {}
    mu_q_pow_arr = {}
    for q in range(1, q_max + 1):
        d_q_arr[q] = _d_q_transverse(q)
        mu_q_pow_arr[q] = _mu_q(q) ** p

    for n2 in range(n_max + 1):
        g2 = g_arr[n2]
        l2_pow = lam_pow_arr[n2]

        # Build allowed (n1, q1) for vertex 1
        v1_data = []
        for n1 in range(n_max + 1):
            for q1 in range(1, min(q_max, n1 + n2) + 1):
                if not _vertex_allowed(n1, n2, q1):
                    continue
                w1 = _so4_channel_count(n1, n2, q1) if use_cg_weights else 1
                if w1 == 0:
                    continue
                v1_data.append((n1, q1, w1))

        # Build allowed (n3, q2) for vertex 2
        v2_data = []
        for n3 in range(n_max + 1):
            for q2 in range(1, min(q_max, n2 + n3) + 1):
                if not _vertex_allowed(n2, n3, q2):
                    continue
                w2 = _so4_channel_count(n2, n3, q2) if use_cg_weights else 1
                if w2 == 0:
                    continue
                v2_data.append((n3, q2, w2))

        if not v1_data or not v2_data:
            continue

        for n1, q1, w1 in v1_data:
            g1 = g_arr[n1]
            l1_pow = lam_pow_arr[n1]
            dq1 = d_q_arr[q1]
            mq1 = mu_q_pow_arr[q1]
            factor1 = mpmath.mpf(w1) * g1 * dq1 / (l1_pow * mq1)

            for n3, q2, w2 in v2_data:
                g3 = g_arr[n3]
                l3_pow = lam_pow_arr[n3]
                dq2 = d_q_arr[q2]
                mq2 = mu_q_pow_arr[q2]

                contrib = factor1 * mpmath.mpf(w2) * g2 * g3 * dq2 / (l2_pow * l3_pow * mq2)
                total += contrib
                n_quintuples += 1

    return {
        "value": total,
        "value_float": float(total),
        "n_max": n_max,
        "electron_exponent": a,
        "photon_exponent": p,
        "use_cg_weights": use_cg_weights,
        "n_quintuples": n_quintuples,
    }


def three_loop_unrestricted(
    electron_exponent: int = 4,
) -> Dict[str, object]:
    """Three-loop unrestricted: factorizes as D(a)^3 (no vertex restriction).

    S_unr = D(a)^3 where D(a) = sum_n g_n / lambda_n^a.

    At a=4: D(4)^3 = (pi^2 - pi^4/12)^3, purely pi^{even}.

    Parameters
    ----------
    electron_exponent : int
        Exponent on electron propagators.

    Returns
    -------
    Dict with exact Hurwitz value and decomposition.
    """
    a = electron_exponent
    D = _dirac_D(a)
    D_cubed = D ** 3

    return {
        "D_a": D,
        "D_a_float": float(D),
        "D_a_cubed": D_cubed,
        "D_a_cubed_float": float(D_cubed),
        "a": a,
        "transcendental_class": "pi^{even}" if a % 2 == 0 else "odd-zeta",
    }


def three_loop_vertex_restricted(
    n_max: int,
    electron_exponent: int = 4,
    photon_exponent: int = 1,
) -> Dict[str, object]:
    """Three-loop with vertex parity restriction only (W=1 for allowed triples).

    Parameters
    ----------
    n_max : int
        Truncation level.
    electron_exponent : int
        Exponent on electron propagators.
    photon_exponent : int
        Exponent on photon propagators.

    Returns
    -------
    Dict with the sum and decomposition attempt.
    """
    result = three_loop_sunset_s3(
        n_max=n_max,
        electron_exponent=electron_exponent,
        photon_exponent=photon_exponent,
        use_cg_weights=False,
    )
    decomp = decompose_three_loop(result["value"])
    result["decomposition"] = decomp
    return result


def three_loop_cg_weighted(
    n_max: int,
    electron_exponent: int = 4,
    photon_exponent: int = 1,
) -> Dict[str, object]:
    """Three-loop with full SO(4) CG channel count weights.

    Parameters
    ----------
    n_max : int
        Truncation level.
    electron_exponent : int
        Exponent on electron propagators.
    photon_exponent : int
        Exponent on photon propagators.

    Returns
    -------
    Dict with the sum, decomposition, and comparison to unrestricted.
    """
    result = three_loop_sunset_s3(
        n_max=n_max,
        electron_exponent=electron_exponent,
        photon_exponent=photon_exponent,
        use_cg_weights=True,
    )
    decomp = decompose_three_loop(result["value"])
    result["decomposition"] = decomp

    # Compute unrestricted for comparison
    unr = three_loop_unrestricted(electron_exponent=electron_exponent)
    result["unrestricted_value"] = unr["D_a_cubed"]
    result["unrestricted_float"] = unr["D_a_cubed_float"]
    result["difference"] = result["value"] - unr["D_a_cubed"]
    result["difference_float"] = float(result["value"] - unr["D_a_cubed"])

    # Decompose the difference (where new transcendentals would live)
    decomp_diff = decompose_three_loop(result["difference"])
    result["decomposition_difference"] = decomp_diff

    return result


def three_loop_convergence_table(
    n_max_values: Optional[List[int]] = None,
    electron_exponent: int = 4,
    photon_exponent: int = 1,
    use_cg_weights: bool = True,
) -> List[Dict[str, object]]:
    """Compute three-loop sum at increasing n_max to assess convergence.

    Parameters
    ----------
    n_max_values : list of int, optional
        List of n_max values. Default: [3, 5, 7, 10, 12, 15].
    electron_exponent : int
        Exponent on electron propagators.
    photon_exponent : int
        Exponent on photon propagators.
    use_cg_weights : bool
        If True, use CG weights.

    Returns
    -------
    List of dicts, one per n_max value.
    """
    if n_max_values is None:
        n_max_values = [3, 5, 7, 10, 12, 15]

    rows = []
    prev_value = None
    for nm in n_max_values:
        r = three_loop_sunset_s3(
            n_max=nm,
            electron_exponent=electron_exponent,
            photon_exponent=photon_exponent,
            use_cg_weights=use_cg_weights,
        )
        delta = None
        if prev_value is not None:
            delta = float(abs(r["value"] - prev_value))
        prev_value = r["value"]
        rows.append({
            "n_max": nm,
            "value_float": r["value_float"],
            "n_quintuples": r["n_quintuples"],
            "delta_from_previous": delta,
        })
    return rows


# ---------------------------------------------------------------------------
# PSLQ decomposition with extended MZV basis
# ---------------------------------------------------------------------------

def _dirichlet_beta(s: int) -> mpmath.mpf:
    """Dirichlet beta function beta(s) via Hurwitz zeta.

    beta(s) = sum_{n=0}^inf (-1)^n / (2n+1)^s
            = (1/4^s) * [zeta(s, 1/4) - zeta(s, 3/4)]
    """
    return (mpmath.hurwitz(s, mpmath.mpf(1)/4)
            - mpmath.hurwitz(s, mpmath.mpf(3)/4)) / mpmath.mpf(4)**s


def _mzv_3_1_1() -> mpmath.mpf:
    """Multiple zeta value zeta(3,1,1) = sum_{n1>n2>n3>=1} 1/(n1^3 * n2 * n3).

    Known identity: zeta(3,1,1) = zeta(5) - zeta(2)*zeta(3) + (11/2)*zeta(5)/??
    Actually: zeta(3,1,1) is weight 5, depth 3.
    Known: zeta(3,1,1) = (1/12)*pi^2*zeta(3) - (11/2)*zeta(5) + ...
    Let's compute numerically by direct summation.
    """
    # Exact: zeta({3,1,1}) using sum_{n>m>k>=1} 1/(n^3 * m * k)
    # This is very slow. Use the shuffle/stuffle identity.
    # From Zagier's MZV database:
    # zeta(3,1,1) = 2*zeta(5) - zeta(2)*zeta(3)
    # Let's verify by computing a partial sum.
    # Actually from the standard stuffle relations:
    # zeta(3,1,1) = (pi^2/6)*zeta(3) - 3*zeta(5)  [Not sure, compute directly]

    # Direct computation (slow but reliable at moderate N)
    N = 200
    total = mpmath.mpf(0)
    for n1 in range(3, N + 1):
        for n2 in range(2, n1):
            h_n2 = mpmath.mpf(0)
            for n3 in range(1, n2):
                h_n2 += mpmath.mpf(1) / n3
            total += h_n2 / (mpmath.mpf(n1)**3 * n2)
    return total


def _mzv_2_1_1() -> mpmath.mpf:
    """Multiple zeta value zeta(2,1,1) = sum_{n1>n2>n3>=1} 1/(n1^2 * n2 * n3).

    Known: zeta(2,1,1) = zeta(4) = pi^4/90.
    (This is a theorem: all depth-d MZVs of weight d+1 reduce to single zeta.)
    """
    return mpmath.zeta(4)  # pi^4/90


def _li4_half() -> mpmath.mpf:
    """Li_4(1/2) = sum_{n=1}^inf 1/(n^4 * 2^n).

    Known in terms of: (7/4)*zeta(4) - (1/12)*pi^2*log(2)^2 + (1/24)*log(2)^4 - ...
    Actually Li_4(1/2) is a specific constant. Use mpmath.polylog.
    """
    return mpmath.polylog(4, mpmath.mpf(1)/2)


def decompose_three_loop(
    value: mpmath.mpf,
    *,
    tol: float = 1e-25,
    maxcoeff: int = 100000,
) -> Dict[str, object]:
    """Decompose a three-loop sum value into MZV basis using PSLQ.

    The basis is organized by transcendental weight:
    - Weight 0: rational (1)
    - Weight 2: pi^2
    - Weight 4: pi^4, G^2 (Catalan squared, but G has weight 2)
    - Weight 6: pi^6, pi^2*G^2
    - Odd zeta: zeta(3), zeta(5), zeta(7)
    - Products: zeta(3)^2, zeta(3)*zeta(5), pi^2*zeta(3), pi^2*zeta(5)
    - Catalan: G, beta(4), pi^2*G, pi^2*beta(4)
    - Depth-3 MZV: zeta(3,1,1) (computed numerically)
    - Polylog: Li_4(1/2)

    Parameters
    ----------
    value : mpmath.mpf
        The value to decompose.
    tol : float
        PSLQ tolerance.
    maxcoeff : int
        Maximum PSLQ coefficient.

    Returns
    -------
    Dict with identification status, components, and residual.
    """
    pi2 = mpmath.pi**2
    pi4 = mpmath.pi**4
    pi6 = mpmath.pi**6
    z3 = mpmath.zeta(3)
    z5 = mpmath.zeta(5)
    z7 = mpmath.zeta(7)
    G = mpmath.catalan
    beta4 = _dirichlet_beta(4)
    li4_half = _li4_half()
    ln2 = mpmath.log(2)

    # --- Tier 1: basic basis (no depth-3 MZVs) ---
    basis_1 = [value, mpmath.mpf(1), pi2, pi4, pi6,
               z3, z5, z7,
               z3**2, z3 * z5,
               pi2 * z3, pi2 * z5,
               G, beta4,
               pi2 * G, pi2 * beta4,
               G * z3, beta4 * z3]
    labels_1 = ["val", "1", "pi^2", "pi^4", "pi^6",
                "zeta(3)", "zeta(5)", "zeta(7)",
                "zeta(3)^2", "zeta(3)*zeta(5)",
                "pi^2*zeta(3)", "pi^2*zeta(5)",
                "G", "beta(4)",
                "pi^2*G", "pi^2*beta(4)",
                "G*zeta(3)", "beta(4)*zeta(3)"]

    result = _try_pslq(value, basis_1, labels_1, tol, maxcoeff)
    if result["identified"]:
        result["basis_tier"] = 1
        result["note"] = "Identified in standard MZV/Dirichlet basis (no depth-3)"
        return result

    # --- Tier 2: add depth-3 MZVs and Li_4(1/2) ---
    # Compute zeta(3,1,1) numerically
    mzv_311 = _mzv_3_1_1()

    basis_2 = basis_1 + [mzv_311, li4_half, ln2, pi2 * ln2,
                         z3 * G, z3 * beta4, ln2**2, ln2**4]
    labels_2 = labels_1 + ["zeta(3,1,1)", "Li4(1/2)", "ln(2)", "pi^2*ln(2)",
                           "zeta(3)*G", "zeta(3)*beta(4)", "ln(2)^2", "ln(2)^4"]

    result = _try_pslq(value, basis_2, labels_2, tol, maxcoeff)
    if result["identified"]:
        result["basis_tier"] = 2
        result["note"] = "Identified with depth-3 MZV / Li_4(1/2) basis"
        return result

    # --- Tier 3: just try with a smaller maxcoeff and relax tolerance ---
    result = _try_pslq(value, basis_2, labels_2, 1e-15, 1000000)
    if result["identified"]:
        result["basis_tier"] = 3
        result["note"] = "Identified at relaxed tolerance"
        return result

    return {
        "identified": False,
        "value_float": float(value),
        "note": (
            "PSLQ failed with both standard and extended MZV bases. "
            "Likely insufficient convergence (n_max too small) or the "
            "value contains a genuinely new transcendental outside the basis."
        ),
        "basis_tier": None,
    }


def _try_pslq(
    value: mpmath.mpf,
    basis: list,
    labels: list,
    tol: float,
    maxcoeff: int,
) -> Dict[str, object]:
    """Attempt PSLQ identification of value in given basis."""
    try:
        relation = mpmath.pslq(basis, tol=tol, maxcoeff=maxcoeff)
    except Exception as e:
        return {"identified": False, "error": str(e), "value_float": float(value)}

    if relation is None or relation[0] == 0:
        return {"identified": False, "value_float": float(value)}

    components = {}
    reconstructed = mpmath.mpf(0)
    for i in range(1, len(relation)):
        if relation[i] != 0:
            coeff = mpmath.mpf(-relation[i]) / mpmath.mpf(relation[0])
            components[labels[i]] = str(mpmath.fraction(-relation[i], relation[0]))
            reconstructed += coeff * basis[i]

    residual = float(abs(value - reconstructed))

    return {
        "identified": True,
        "raw_relation": relation,
        "components": components,
        "reconstructed_float": float(reconstructed),
        "value_float": float(value),
        "residual": residual,
        "contains_zeta3_squared": "zeta(3)^2" in components,
        "contains_depth3_mzv": "zeta(3,1,1)" in components,
        "contains_catalan": any("G" in k for k in components),
        "contains_dirichlet_beta": any("beta" in k for k in components),
        "new_beyond_two_loop": any(
            k in components
            for k in ["zeta(3)^2", "zeta(3)*zeta(5)", "zeta(7)",
                      "zeta(3,1,1)", "Li4(1/2)",
                      "G*zeta(3)", "beta(4)*zeta(3)",
                      "zeta(3)*G", "zeta(3)*beta(4)"]
        ),
    }


# ---------------------------------------------------------------------------
# Transcendental classification
# ---------------------------------------------------------------------------

def classify_three_loop_transcendentals() -> Dict[str, str]:
    """Classify the expected three-loop transcendental content."""
    return {
        "unrestricted_three_loop": (
            "pi^{even}: D(4)^3 = (pi^2 - pi^4/12)^3, a polynomial in pi^2 "
            "with rational coefficients. No odd-zeta, no Catalan, no MZV."
        ),
        "vertex_restricted_three_loop": (
            "EXPECTED: D_even and D_odd parity filtration at EACH vertex "
            "produces even/odd Hurwitz shifts at quarter-integer lattice. "
            "Two vertices create a NESTED parity filtration: even-even, "
            "even-odd, odd-even, odd-odd subcases. The nested Hurwitz "
            "structure could produce depth-2 Dirichlet L-values (G, beta(4)) "
            "from the vertex topology, as at two loops."
        ),
        "cg_weighted_three_loop": (
            "EXPECTED: The SO(4) CG channel count creates min(n1,n2) and "
            "min(n2,n3) weights at the two vertices, producing a nested "
            "min-weighted sum. At two loops, the single min-weight produced "
            "depth-2 multiple Hurwitz zeta outside standard bases. At three "
            "loops, the nested min-weights are expected to produce depth-3 "
            "multiple Hurwitz zeta values."
        ),
        "depth_k_prediction": (
            "At k loops with vertex restrictions and CG weights, expect "
            "depth-k multiple zeta values. One loop: depth-1 (pi^{even}). "
            "Two loops: depth-2 (unidentified multiple Hurwitz). "
            "Three loops: depth-3 predicted."
        ),
        "flat_space_connection": (
            "In flat-space QED, three-loop corrections involve zeta(3), zeta(5), "
            "pi^6, and Li_4(1/2) (Laporta 2002, Bailey-Broadhurst 2001). "
            "On S^3, the half-integer Dirac spectrum replaces integer sums "
            "with Hurwitz sums, potentially producing Hurwitz analogs of "
            "flat-space MZVs."
        ),
    }


# ---------------------------------------------------------------------------
# Factorized O(N^3) three-loop sum (Track Q-2)
# ---------------------------------------------------------------------------

def three_loop_factorized(
    n_max: int,
    a: int = 4,
    p: int = 1,
    use_cg_weights: bool = True,
) -> Dict[str, object]:
    """Factorized O(N^3) three-loop chain sum on S^3.

    The five-fold sum over (n1, n2, n3, q1, q2) factors as:

        S = sum_{n2=0}^{n_max} g(n2) / lambda(n2)^a * L(n2) * R(n2)

    where L(n2) = sum_{n1,q1} W(n1,n2,q1) * g(n1) * d_T(q1) / (lambda(n1)^a * mu(q1)^p)
      and R(n2) = sum_{n3,q2} W(n2,n3,q2) * g(n3) * d_T(q2) / (lambda(n3)^a * mu(q2)^p)

    Because the vertex weights and spectrum are the same at both vertices,
    L(n2) and R(n2) have the same functional form but with n2 in the "other"
    slot. However, the vertex coupling W(n1, n2, q) is NOT symmetric in
    n1 and n2 in general (the SU(2)_L x SU(2)_R reps differ for the two
    electron lines). So we compute L(n2) and R(n2) separately.

    For the LEFT partial sum, n2 is in position 2 of W(n1, n2, q1):
        L(n2) = sum_{n1} sum_{q1} W(n1, n2, q1) * g(n1) * d_T(q1) / (lam(n1)^a * mu(q1)^p)

    For the RIGHT partial sum, n2 is in position 1 of W(n2, n3, q2):
        R(n2) = sum_{n3} sum_{q2} W(n2, n3, q2) * g(n3) * d_T(q2) / (lam(n3)^a * mu(q2)^p)

    Complexity: O(N^3) instead of O(N^5).

    Parameters
    ----------
    n_max : int
        Truncation for electron modes (0..n_max) and photon modes (1..2*n_max).
    a : int
        Exponent on each electron propagator |lambda_n|^a. Default 4.
    p : int
        Exponent on each photon propagator mu_q^p. Default 1.
    use_cg_weights : bool
        If True, use SO(4) CG channel count W(n1,n2,q). If False, W=1.

    Returns
    -------
    Dict with the factorized sum value and timing.
    """
    t0 = time.perf_counter()

    q_max = 2 * n_max

    # Precompute electron data
    g_arr = [_g_n_dirac(n) for n in range(n_max + 1)]
    lam_pow_arr = [_lambda_n(n) ** a for n in range(n_max + 1)]

    # Precompute photon data
    d_q_cache = {}
    mu_q_pow_cache = {}
    for q in range(1, q_max + 1):
        d_q_cache[q] = _d_q_transverse(q)
        mu_q_pow_cache[q] = _mu_q(q) ** p

    # Compute LEFT partial sums: L(n2) where n2 is in the second slot of W
    # L(n2) = sum_{n1, q1} W(n1, n2, q1) * g(n1) * d_T(q1) / (lam(n1)^a * mu(q1)^p)
    L = [mpmath.mpf(0)] * (n_max + 1)
    for n2 in range(n_max + 1):
        total = mpmath.mpf(0)
        for n1 in range(n_max + 1):
            g1 = g_arr[n1]
            l1_pow = lam_pow_arr[n1]
            q_lo = max(1, abs(n1 - n2))
            q_hi = min(q_max, n1 + n2)
            for q in range(q_lo, q_hi + 1):
                if (n1 + n2 + q) % 2 == 0:
                    continue
                if use_cg_weights:
                    w = _so4_channel_count(n1, n2, q)
                    if w == 0:
                        continue
                else:
                    w = 1
                dq = d_q_cache[q]
                mq_pow = mu_q_pow_cache[q]
                total += mpmath.mpf(w) * g1 * dq / (l1_pow * mq_pow)
        L[n2] = total

    # Compute RIGHT partial sums: R(n2) where n2 is in the first slot of W
    # R(n2) = sum_{n3, q2} W(n2, n3, q2) * g(n3) * d_T(q2) / (lam(n3)^a * mu(q2)^p)
    R = [mpmath.mpf(0)] * (n_max + 1)
    for n2 in range(n_max + 1):
        total = mpmath.mpf(0)
        for n3 in range(n_max + 1):
            g3 = g_arr[n3]
            l3_pow = lam_pow_arr[n3]
            q_lo = max(1, abs(n2 - n3))
            q_hi = min(q_max, n2 + n3)
            for q in range(q_lo, q_hi + 1):
                if (n2 + n3 + q) % 2 == 0:
                    continue
                if use_cg_weights:
                    w = _so4_channel_count(n2, n3, q)
                    if w == 0:
                        continue
                else:
                    w = 1
                dq = d_q_cache[q]
                mq_pow = mu_q_pow_cache[q]
                total += mpmath.mpf(w) * g3 * dq / (l3_pow * mq_pow)
        R[n2] = total

    # Combine: S = sum_{n2} g(n2) / lam(n2)^a * L(n2) * R(n2)
    S = mpmath.mpf(0)
    for n2 in range(n_max + 1):
        S += g_arr[n2] / lam_pow_arr[n2] * L[n2] * R[n2]

    elapsed = time.perf_counter() - t0

    return {
        "value": S,
        "value_float": float(S),
        "n_max": n_max,
        "electron_exponent": a,
        "photon_exponent": p,
        "use_cg_weights": use_cg_weights,
        "time_seconds": elapsed,
    }


def three_loop_factorized_convergence(
    n_max_values: Optional[List[int]] = None,
    a: int = 4,
    p: int = 1,
) -> List[Dict[str, object]]:
    """Convergence table for the factorized three-loop CG-weighted sum.

    Computes the factorized sum at increasing n_max values to assess
    convergence and enable PSLQ identification at high precision.

    Parameters
    ----------
    n_max_values : list of int, optional
        Default: [10, 20, 30, 50, 75, 100].
    a : int
        Electron propagator exponent. Default 4.
    p : int
        Photon propagator exponent. Default 1.

    Returns
    -------
    List of dicts with n_max, value, delta, time_seconds.
    """
    if n_max_values is None:
        n_max_values = [10, 20, 30, 50, 75, 100]

    rows = []
    prev_value = None
    for nm in n_max_values:
        r = three_loop_factorized(n_max=nm, a=a, p=p, use_cg_weights=True)
        delta = None
        if prev_value is not None:
            delta = float(abs(r["value"] - prev_value))
        prev_value = r["value"]
        rows.append({
            "n_max": nm,
            "value": r["value"],
            "value_float": r["value_float"],
            "delta_from_previous": delta,
            "time_seconds": r["time_seconds"],
        })
    return rows


def three_loop_euler_maclaurin_tail(
    n_max: int,
    n_em_terms: int = 10,
    a: int = 4,
    p: int = 1,
) -> mpmath.mpf:
    """Euler-Maclaurin tail correction for the factorized three-loop sum.

    The partial sum increments delta(N) = S(N) - S(N-1) decrease as a
    power law in N. We fit the last few increments to a model
    delta(N) = c * N^alpha and integrate the tail analytically:

        tail = integral_{N}^{inf} c * x^alpha dx = -c * N^{alpha+1} / (alpha+1)

    This works when alpha < -1 (convergent tail). If the increments don't
    decrease fast enough, we return a conservative estimate from the last
    increment alone.

    Parameters
    ----------
    n_max : int
        The base truncation level.
    n_em_terms : int
        Number of increments to fit (minimum 3).
    a : int
        Electron propagator exponent.
    p : int
        Photon propagator exponent.

    Returns
    -------
    mpmath.mpf
        Estimated tail correction S(inf) - S(n_max).
    """
    # Compute a sequence of closely-spaced partial sums
    n_pts = min(n_em_terms, 8)
    s_vals = []

    for i in range(n_pts + 1):
        n_eval = n_max - i
        if n_eval < 3:
            break
        r = three_loop_factorized(n_max=n_eval, a=a, p=p, use_cg_weights=True)
        s_vals.append((n_eval, r["value"]))

    if len(s_vals) < 3:
        return mpmath.mpf(0)

    # Compute increments delta(N) = S(N) - S(N-1)
    deltas = []
    for i in range(len(s_vals) - 1):
        N_i = s_vals[i][0]
        d_i = s_vals[i][1] - s_vals[i + 1][1]
        deltas.append((N_i, d_i))

    if len(deltas) < 2:
        return mpmath.mpf(0)

    # Fit power law: log(delta) = log(c) + alpha * log(N)
    # Use the first and last increments for a simple fit
    N1, d1 = deltas[0]     # delta at n_max
    N2, d2 = deltas[-1]    # delta at n_max - (n_pts - 1)

    if d1 <= 0 or d2 <= 0:
        return mpmath.mpf(0)

    log_d1 = mpmath.log(d1)
    log_d2 = mpmath.log(d2)
    log_N1 = mpmath.log(mpmath.mpf(N1))
    log_N2 = mpmath.log(mpmath.mpf(N2))

    denom_log = log_N1 - log_N2
    if abs(denom_log) < mpmath.mpf(10) ** (-50):
        return mpmath.mpf(0)

    alpha = (log_d1 - log_d2) / denom_log
    log_c = log_d1 - alpha * log_N1
    c = mpmath.exp(log_c)

    # The tail is sum_{k=N+1}^{inf} c * k^alpha ~ integral c * x^alpha dx
    # = c * N^{alpha+1} / (-(alpha+1)) for alpha < -1
    if alpha >= -1:
        # Sum is divergent or very slowly convergent
        # Return a conservative lower bound: just the next increment
        return c * mpmath.mpf(n_max + 1) ** alpha

    tail = c * mpmath.mpf(n_max) ** (alpha + 1) / (-(alpha + 1))
    return tail


def decompose_three_loop_mzv(
    value: mpmath.mpf,
    n_digits: int = 50,
) -> Dict[str, object]:
    """Extended PSLQ decomposition with depth-3 MZV basis.

    The basis includes all expected three-loop transcendentals:
    - Rational: 1
    - Weight 2: pi^2
    - Weight 4: pi^4
    - Weight 6: pi^6
    - Odd zeta: zeta(3), zeta(5), zeta(7)
    - Products: pi^2*zeta(3), pi^2*zeta(5), zeta(3)^2
    - Depth-3 MZV: zeta(3,1,1) (weight 5, depth 3)
    - Polylog: Li_4(1/2)

    Parameters
    ----------
    value : mpmath.mpf
        The value to decompose.
    n_digits : int
        Number of reliable digits in value. Controls PSLQ tolerance.

    Returns
    -------
    Dict with identification status, components, and residual.
    """
    tol = mpmath.mpf(10) ** (-(n_digits - 5))

    pi2 = mpmath.pi ** 2
    pi4 = mpmath.pi ** 4
    pi6 = mpmath.pi ** 6
    z3 = mpmath.zeta(3)
    z5 = mpmath.zeta(5)
    z7 = mpmath.zeta(7)
    G = mpmath.catalan
    beta4 = _dirichlet_beta(4)
    li4_half = _li4_half()
    ln2 = mpmath.log(2)

    # Depth-3 MZV: zeta(3,1,1)
    # Known identity: zeta(3,1,1) = 2*zeta(5) - zeta(2)*zeta(3)
    # (from Hoffman's identity / stuffle relations for weight 5, depth 3)
    # Use the identity rather than slow direct summation for better precision.
    mzv_311 = 2 * z5 - (pi2 / 6) * z3

    pi8 = mpmath.pi ** 8
    pi10 = mpmath.pi ** 10
    pi12 = mpmath.pi ** 12

    # --- Tier 0: pi^{even} only (fast, catches pure polynomial in pi^2) ---
    basis_0 = [value, mpmath.mpf(1), pi2, pi4, pi6, pi8, pi10, pi12]
    labels_0 = ["val", "1", "pi^2", "pi^4", "pi^6", "pi^8", "pi^10", "pi^12"]

    result = _try_pslq(value, basis_0, labels_0, float(tol), 100000)
    if result["identified"]:
        result["basis_tier"] = 0
        result["note"] = "Identified as polynomial in pi^2 (pi^{even} only)"
        return result

    # --- Tier 1: standard MZV/Dirichlet basis ---
    # Keep the basis compact (<=16 elements) for PSLQ reliability
    basis_1 = [
        value, mpmath.mpf(1), pi2, pi4, pi6,
        z3, z5, z7,
        z3 ** 2, z3 * z5,
        pi2 * z3, pi2 * z5,
        G, beta4,
        pi2 * G, pi2 * beta4,
    ]
    labels_1 = [
        "val", "1", "pi^2", "pi^4", "pi^6",
        "zeta(3)", "zeta(5)", "zeta(7)",
        "zeta(3)^2", "zeta(3)*zeta(5)",
        "pi^2*zeta(3)", "pi^2*zeta(5)",
        "G", "beta(4)",
        "pi^2*G", "pi^2*beta(4)",
    ]

    result = _try_pslq(value, basis_1, labels_1, float(tol), 100000)
    if result["identified"]:
        result["basis_tier"] = 1
        result["note"] = "Identified in standard MZV/Dirichlet basis (no depth-3)"
        return result

    # --- Tier 2: add depth-3 MZV, Li_4(1/2), ln(2) products ---
    basis_2 = basis_1 + [
        mzv_311, li4_half, ln2, pi2 * ln2,
        z3 * G, z3 * beta4, ln2 ** 2, ln2 ** 4,
        G ** 2, G * beta4, beta4 ** 2,
    ]
    labels_2 = labels_1 + [
        "zeta(3,1,1)", "Li4(1/2)", "ln(2)", "pi^2*ln(2)",
        "zeta(3)*G", "zeta(3)*beta(4)", "ln(2)^2", "ln(2)^4",
        "G^2", "G*beta(4)", "beta(4)^2",
    ]

    result = _try_pslq(value, basis_2, labels_2, float(tol), 100000)
    if result["identified"]:
        result["basis_tier"] = 2
        result["note"] = "Identified with depth-3 MZV / polylog basis"
        return result

    # --- Tier 3: relax tolerance for lower-precision values ---
    result = _try_pslq(value, basis_2, labels_2, 1e-15, 1000000)
    if result["identified"]:
        result["basis_tier"] = 3
        result["note"] = "Identified at relaxed tolerance (lower confidence)"
        return result

    return {
        "identified": False,
        "value_float": float(value),
        "n_digits_input": n_digits,
        "note": (
            "PSLQ failed with both standard and extended MZV bases at "
            f"{n_digits} digits. Either the value contains a genuinely "
            "new transcendental outside the basis, or more digits are needed."
        ),
        "basis_tier": None,
    }
