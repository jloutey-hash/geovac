"""One-loop electron self-energy on S^3 from spectral mode sums.

This is the first physical-process QED calculation in GeoVac. Existing
modules (qed_vacuum_polarization, qed_two_loop, qed_vertex) compute
vacuum diagrams; this module computes a diagram with an external state.

Physics
-------
The one-loop electron self-energy on S^3 (Feynman gauge) is the spectral
mode sum:

    Sigma(n_ext) = sum_{n_int, q} W(n_ext, n_int, q)
                   * g(n_int) * d_T(q)
                   / (|lambda(n_int)|^{2*s_e} * mu(q)^{s_gamma})

where:
  - n_ext: external electron Dirac level (CH convention, n >= 0)
  - n_int: internal (virtual) electron level (CH convention, n >= 0)
  - q: internal photon level (q >= 1)
  - W(n1, n2, q): SO(4) vertex coupling weight (channel count 0, 1, or 2)
    from ``geovac.qed_vertex.weighted_vertex_coupling``
  - g(n) = 2(n+1)(n+2): full Dirac degeneracy (internal line)
  - d_T(q) = q(q+2): transverse photon degeneracy
  - |lambda(n)| = n + 3/2: absolute Dirac eigenvalue (CH convention)
  - mu(q) = q(q+2): Hodge-1 Laplacian eigenvalue (photon propagator)
  - s_e, s_gamma: propagator exponents (physical: s_e=2, s_gamma=1)

The vertex selection rule (SO(4) triangle + parity) restricts:
  - |n_ext - n_int| <= q <= n_ext + n_int
  - n_ext + n_int + q is odd

At the physical propagator powers (s_e=2, s_gamma=1), the internal
electron propagator enters as 1/|lambda|^4 and the photon as 1/mu_q.

The self-energy depends on the external state n_ext through the vertex
selection rule, which constrains which (n_int, q) pairs can couple to
the external electron. This n_ext dependence is the self-energy's
state-dependent mass renormalization.

Vertex correction (stretch goal)
--------------------------------
The one-loop vertex correction involves TWO vertex insertions and a
double spectral sum over two internal lines. This is O(N^4) in the
cutoff. Implemented as ``vertex_correction_spectral``.

Transcendental taxonomy (Paper 18)
----------------------------------
The self-energy sum mixes:
  - Dirac spectrum (half-integer eigenvalues -> Hurwitz zeta)
  - Hodge-1 spectrum (integer eigenvalues -> products of linear factors)
The vertex selection rule breaks factorization, so the transcendental
content is expected to be richer than single Dirichlet series. At one
loop with the physical propagator powers, the dominant contribution
is pi^{even} with possible Dirichlet beta contamination from the
vertex parity constraint (same mechanism as qed_vertex.py).

References
----------
- Camporesi & Higuchi, J. Geom. Phys. 20 (1996) 1-18 [Dirac on S^3].
- GeoVac qed_vertex.py (SO(4) vertex coupling, selection rules).
- GeoVac qed_two_loop.py (PSLQ decomposition patterns).
- GeoVac Paper 18 (exchange constant taxonomy).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import mpmath

# High precision for PSLQ identification
mpmath.mp.dps = 50

__all__ = [
    "self_energy_spectral",
    "self_energy_table",
    "self_energy_convergence",
    "self_energy_transcendental_class",
    "vertex_correction_spectral",
    "schwinger_convergence",
]


# ---------------------------------------------------------------------------
# Spectrum helpers (mpmath precision, replicating qed_vertex conventions)
# ---------------------------------------------------------------------------

def _lambda_n(n: int) -> mpmath.mpf:
    """Absolute Dirac eigenvalue |lambda_n| = n + 3/2 (CH convention)."""
    return mpmath.mpf(n) + mpmath.mpf(3) / 2


def _g_n_dirac(n: int) -> mpmath.mpf:
    """Full Dirac degeneracy g_n = 2(n+1)(n+2)."""
    return mpmath.mpf(2) * (n + 1) * (n + 2)


def _mu_q(q: int) -> mpmath.mpf:
    """Hodge-1 Laplacian eigenvalue mu_q = q(q+2), q >= 1."""
    return mpmath.mpf(q) * (q + 2)


def _d_q_transverse(q: int) -> mpmath.mpf:
    """Transverse (physical) photon degeneracy d_q^T = q(q+2)."""
    return mpmath.mpf(q) * (q + 2)


def _vertex_allowed(n1: int, n2: int, n_gamma: int) -> bool:
    """Check SO(4) vertex selection rule.

    Requires:
      - Triangle inequality: |n1 - n2| <= n_gamma <= n1 + n2
      - Parity: n1 + n2 + n_gamma is odd (gamma-matrix flips parity)
      - n_gamma >= 1 (photon modes start at q=1)
    """
    if n_gamma < 1:
        return False
    if n_gamma < abs(n1 - n2):
        return False
    if n_gamma > n1 + n2:
        return False
    if (n1 + n2 + n_gamma) % 2 == 0:
        return False
    return True


# ---------------------------------------------------------------------------
# SO(4) channel count (reused from qed_vertex.py logic)
# ---------------------------------------------------------------------------

def _so4_channel_count(n1: int, n2: int, n_gamma: int) -> int:
    """Count the number of SO(4) vector harmonic components that couple
    the Dirac spinor at level n1 to level n2 via a photon at level n_gamma.

    Returns 0, 1, or 2. See qed_vertex.so4_channel_count for full docs.
    """
    if not _vertex_allowed(n1, n2, n_gamma):
        return 0

    from fractions import Fraction

    j1_L = Fraction(n1 + 1, 2)
    j1_R = Fraction(n1, 2)
    j2_L = Fraction(n2, 2)
    j2_R = Fraction(n2 + 1, 2)

    count = 0
    q = n_gamma

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


def _weighted_vertex_coupling(n1: int, n2: int, n_gamma: int) -> mpmath.mpf:
    """Weighted vertex coupling using the SO(4) channel count.

    Returns the number of SO(4) Clebsch-Gordan channels (0, 1, or 2)
    as an mpmath float.
    """
    return mpmath.mpf(_so4_channel_count(n1, n2, n_gamma))


# ---------------------------------------------------------------------------
# One-loop self-energy
# ---------------------------------------------------------------------------

def self_energy_spectral(
    n_ext: int,
    n_max: int,
    s_e: int = 2,
    s_gamma: int = 1,
) -> mpmath.mpf:
    """One-loop electron self-energy on S^3 from spectral mode sum.

    Computes:
        Sigma(n_ext) = sum_{n_int=0}^{n_max} sum_{q=1}^{q_max}
                       W(n_ext, n_int, q) * g(n_int) * d_T(q)
                       / (|lambda(n_int)|^{2*s_e} * mu(q)^{s_gamma})

    where:
      - W is the SO(4) vertex coupling weight (channel count 0, 1, or 2)
      - The sum is over the internal electron (n_int) and photon (q)
      - q_max is bounded by the triangle inequality: q <= n_ext + n_int

    Parameters
    ----------
    n_ext : int
        External electron level (CH convention, n >= 0).
    n_max : int
        Truncation level for internal modes.
    s_e : int
        Half-exponent on internal electron propagator. The actual exponent
        on |lambda| is 2*s_e. Physical value: s_e=2 gives 1/|lambda|^4.
    s_gamma : int
        Exponent on photon propagator. Physical value: s_gamma=1 gives
        1/mu_q.

    Returns
    -------
    mpmath.mpf
        The self-energy spectral sum Sigma(n_ext).
    """
    if n_ext < 0:
        raise ValueError(f"n_ext must be >= 0, got {n_ext}")
    if n_max < 0:
        raise ValueError(f"n_max must be >= 0, got {n_max}")

    s_eff_e = 2 * s_e  # actual exponent on lambda

    total = mpmath.mpf(0)

    for n_int in range(n_max + 1):
        g_int = _g_n_dirac(n_int)
        lam_int = _lambda_n(n_int)
        lam_int_pow = lam_int ** s_eff_e

        # q bounded by triangle inequality: |n_ext - n_int| <= q <= n_ext + n_int
        q_lo = abs(n_ext - n_int)
        q_hi = n_ext + n_int

        for q in range(max(1, q_lo), q_hi + 1):
            if not _vertex_allowed(n_ext, n_int, q):
                continue

            W = _so4_channel_count(n_ext, n_int, q)
            if W == 0:
                continue

            d_T = _d_q_transverse(q)
            mu = _mu_q(q)

            total += mpmath.mpf(W) * g_int * d_T / (lam_int_pow * mu ** s_gamma)

    return total


def self_energy_table(
    n_ext_values: List[int],
    n_max: int,
    s_e: int = 2,
    s_gamma: int = 1,
) -> List[Dict[str, object]]:
    """Compute self-energy for multiple external states.

    Parameters
    ----------
    n_ext_values : List[int]
        List of external electron levels (CH convention).
    n_max : int
        Truncation level for internal modes.
    s_e : int
        Half-exponent on electron propagator.
    s_gamma : int
        Exponent on photon propagator.

    Returns
    -------
    List of dicts with n_ext, Sigma (mpmath), Sigma_float.
    """
    results = []
    for n_ext in n_ext_values:
        sigma = self_energy_spectral(n_ext, n_max, s_e=s_e, s_gamma=s_gamma)
        results.append({
            "n_ext": n_ext,
            "Sigma": sigma,
            "Sigma_float": float(sigma),
        })
    return results


def self_energy_convergence(
    n_ext: int,
    n_max_values: List[int],
    s_e: int = 2,
    s_gamma: int = 1,
) -> List[Dict[str, object]]:
    """Convergence study: Sigma(n_ext) vs n_max.

    Computes the self-energy at increasing cutoffs and reports the
    difference between consecutive values.

    Parameters
    ----------
    n_ext : int
        External electron level (CH convention).
    n_max_values : List[int]
        List of truncation levels (should be sorted ascending).
    s_e : int
        Half-exponent on electron propagator.
    s_gamma : int
        Exponent on photon propagator.

    Returns
    -------
    List of dicts with n_max, Sigma, Sigma_float, delta (change from previous).
    """
    results = []
    prev_sigma = None

    for n_max in n_max_values:
        sigma = self_energy_spectral(n_ext, n_max, s_e=s_e, s_gamma=s_gamma)
        delta = float(abs(sigma - prev_sigma)) if prev_sigma is not None else None
        results.append({
            "n_max": n_max,
            "Sigma": sigma,
            "Sigma_float": float(sigma),
            "delta": delta,
        })
        prev_sigma = sigma

    return results


def self_energy_unrestricted(
    n_max: int,
    s_e: int = 2,
    s_gamma: int = 1,
) -> mpmath.mpf:
    """Unrestricted self-energy: sum over ALL (n_int, q) without vertex rule.

    This is the product D(2*s_e) * D_photon(s_gamma) where:
      D(s) = sum_n g_n / |lambda_n|^s  (Dirac Dirichlet series)
      D_photon(s) = sum_q d_T(q) / mu_q^s  (photon Dirichlet series)

    Factorization: Sigma_unrestricted = D(2*s_e) * D_photon(s_gamma)

    Parameters
    ----------
    n_max : int
        Truncation level.
    s_e : int
        Half-exponent on electron propagator.
    s_gamma : int
        Exponent on photon propagator.

    Returns
    -------
    mpmath.mpf
        The unrestricted sum.
    """
    s_eff_e = 2 * s_e

    D_electron = mpmath.mpf(0)
    for n in range(n_max + 1):
        D_electron += _g_n_dirac(n) / _lambda_n(n) ** s_eff_e

    D_photon = mpmath.mpf(0)
    for q in range(1, 2 * n_max + 1):
        D_photon += _d_q_transverse(q) / _mu_q(q) ** s_gamma

    return D_electron * D_photon


def self_energy_transcendental_class(
    n_ext: int,
    n_max: int,
    s_e: int = 2,
    s_gamma: int = 1,
) -> Dict[str, object]:
    """PSLQ decomposition of Sigma(n_ext) against transcendental basis.

    Basis: {1, pi^2, pi^4, zeta(3), pi^2*zeta(3), zeta(5),
            Catalan_G, beta(4)}

    Parameters
    ----------
    n_ext : int
        External electron level (CH convention).
    n_max : int
        Truncation level.
    s_e : int
        Half-exponent on electron propagator.
    s_gamma : int
        Exponent on photon propagator.

    Returns
    -------
    Dict with the self-energy value, PSLQ decomposition, and
    transcendental classification.
    """
    old_dps = mpmath.mp.dps
    mpmath.mp.dps = 80

    try:
        sigma = self_energy_spectral(n_ext, n_max, s_e=s_e, s_gamma=s_gamma)

        pi2 = mpmath.pi ** 2
        pi4 = mpmath.pi ** 4
        z3 = mpmath.zeta(3)
        z5 = mpmath.zeta(5)
        G = mpmath.catalan
        beta4 = (mpmath.hurwitz(4, mpmath.mpf(1) / 4)
                 - mpmath.hurwitz(4, mpmath.mpf(3) / 4)) / mpmath.power(4, 4)

        basis = [sigma, mpmath.mpf(1), pi2, pi4, z3, pi2 * z3, z5, G, beta4]
        labels = ["value", "1", "pi^2", "pi^4", "zeta(3)", "pi^2*zeta(3)",
                  "zeta(5)", "Catalan_G", "beta(4)"]

        try:
            relation = mpmath.pslq(basis, tol=1e-30, maxcoeff=100000)
        except Exception:
            relation = None

        decomp: Dict[str, object] = {"identified": False}
        if relation is not None and relation[0] != 0:
            components: Dict[str, str] = {}
            reconstructed = mpmath.mpf(0)
            for i in range(1, len(relation)):
                if relation[i] != 0:
                    coeff = mpmath.mpf(-relation[i]) / mpmath.mpf(relation[0])
                    components[labels[i]] = str(
                        mpmath.fraction(-relation[i], relation[0]))
                    reconstructed += coeff * basis[i]
            residual = float(abs(sigma - reconstructed))
            decomp = {
                "identified": True,
                "components": components,
                "residual": residual,
                "contains_pi_even": any(
                    k in components for k in ["pi^2", "pi^4"]),
                "contains_zeta3": any(
                    "zeta(3)" in k for k in components),
                "contains_zeta5": "zeta(5)" in components,
                "contains_catalan": "Catalan_G" in components,
                "contains_beta4": "beta(4)" in components,
            }

        return {
            "n_ext": n_ext,
            "n_max": n_max,
            "s_e": s_e,
            "s_gamma": s_gamma,
            "Sigma": sigma,
            "Sigma_float": float(sigma),
            "decomposition": decomp,
        }
    finally:
        mpmath.mp.dps = old_dps


# ---------------------------------------------------------------------------
# Vertex correction (stretch goal — O(N^4))
# ---------------------------------------------------------------------------

def vertex_correction_spectral(
    n_ext: int,
    n_max: int,
    s_e: int = 2,
    s_gamma: int = 1,
) -> mpmath.mpf:
    """One-loop vertex correction with TWO vertex insertions.

    The vertex correction involves the external electron at level n_ext
    coupling to an internal electron at n_int1 via photon q1, then
    n_int1 coupling back to n_ext via photon q2. The double spectral
    sum is:

        Lambda(n_ext) = sum_{n_int1, q1, q2}
            W(n_ext, n_int1, q1) * W(n_int1, n_ext, q2)
            * g(n_int1) * d_T(q1) * d_T(q2)
            / (|lambda(n_int1)|^{2*s_e} * mu(q1)^{s_gamma} * mu(q2)^{s_gamma})

    This is O(N^3) in n_max (n_int1, q1, q2 summed independently given
    the vertex constraints).

    STRETCH GOAL: This is computationally expensive. n_max <= 10
    recommended.

    Parameters
    ----------
    n_ext : int
        External electron level (CH convention).
    n_max : int
        Truncation level for internal modes.
    s_e : int
        Half-exponent on internal electron propagator.
    s_gamma : int
        Exponent on photon propagator.

    Returns
    -------
    mpmath.mpf
        The vertex correction spectral sum Lambda(n_ext).
    """
    if n_ext < 0:
        raise ValueError(f"n_ext must be >= 0, got {n_ext}")
    if n_max < 0:
        raise ValueError(f"n_max must be >= 0, got {n_max}")

    s_eff_e = 2 * s_e

    total = mpmath.mpf(0)

    for n_int1 in range(n_max + 1):
        g_int1 = _g_n_dirac(n_int1)
        lam_int1 = _lambda_n(n_int1)
        lam_int1_pow = lam_int1 ** s_eff_e

        # First vertex: n_ext -> n_int1 via q1
        q1_lo = abs(n_ext - n_int1)
        q1_hi = n_ext + n_int1

        for q1 in range(max(1, q1_lo), q1_hi + 1):
            if not _vertex_allowed(n_ext, n_int1, q1):
                continue
            W1 = _so4_channel_count(n_ext, n_int1, q1)
            if W1 == 0:
                continue

            d_T1 = _d_q_transverse(q1)
            mu1 = _mu_q(q1)

            # Second vertex: n_int1 -> n_ext via q2
            q2_lo = abs(n_int1 - n_ext)
            q2_hi = n_int1 + n_ext

            for q2 in range(max(1, q2_lo), q2_hi + 1):
                if not _vertex_allowed(n_int1, n_ext, q2):
                    continue
                W2 = _so4_channel_count(n_int1, n_ext, q2)
                if W2 == 0:
                    continue

                d_T2 = _d_q_transverse(q2)
                mu2 = _mu_q(q2)

                total += (mpmath.mpf(W1) * mpmath.mpf(W2)
                          * g_int1 * d_T1 * d_T2
                          / (lam_int1_pow * mu1 ** s_gamma * mu2 ** s_gamma))

    return total


def schwinger_convergence(
    n_max_range: List[int],
    n_ext: int = 0,
    s_e: int = 2,
    s_gamma: int = 1,
) -> List[Dict[str, object]]:
    """Track convergence of the vertex correction toward alpha/(2*pi).

    The Schwinger anomalous magnetic moment result a_e = alpha/(2*pi)
    is the flat-space one-loop vertex correction. On S^3 we compute the
    vertex correction at finite cutoff and track its convergence.

    NOTE: The S^3 spectral sum does NOT directly equal alpha/(2*pi);
    it is a dimensionless spectral sum that should converge to a finite
    value in the large-n_max limit. This function tracks that convergence.

    Parameters
    ----------
    n_max_range : List[int]
        List of cutoff values (should be sorted ascending).
    n_ext : int
        External electron level.
    s_e : int
        Half-exponent on electron propagator.
    s_gamma : int
        Exponent on photon propagator.

    Returns
    -------
    List of dicts with n_max, Lambda, Lambda_float, delta.
    """
    results = []
    prev_val = None

    for n_max in n_max_range:
        val = vertex_correction_spectral(n_ext, n_max, s_e=s_e, s_gamma=s_gamma)
        delta = float(abs(val - prev_val)) if prev_val is not None else None
        results.append({
            "n_max": n_max,
            "Lambda": val,
            "Lambda_float": float(val),
            "delta": delta,
        })
        prev_val = val

    return results
