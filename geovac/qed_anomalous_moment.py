"""Three-point vertex correction and anomalous magnetic moment on S^3.

Extracts the anomalous magnetic moment (g-2)/2 = F_2 from the one-loop
vertex correction on S^3 using m_j-resolved SO(4) Clebsch-Gordan
coefficients.

Physics
-------
The flat-space Schwinger result a_e = alpha/(2*pi) arises from the
one-loop vertex correction diagram. On S^3, the same diagram is a
spectral mode sum over the Dirac and photon towers, but the SO(4)
CG coefficients resolve the m_j dependence that distinguishes the
magnetic form factor F_2 from the charge form factor F_1.

The key object is the magnetic vertex difference:

    B = L(m_j = +1/2) - L(m_j = -1/2)

where L(m_j) is the three-point vertex correction for an external
electron with magnetic quantum number m_j, summed over all internal
electron states and loop photon modes.

The anomalous magnetic moment is:

    F_2 = B / V_tree_magnetic

where V_tree_magnetic is the tree-level magnetic coupling (the
m_j-dependent piece of the tree-level probe-photon vertex).

The SO(4) = SU(2)_L x SU(2)_R CG decomposition uses four CG
coefficients per vertex: two for the state decomposition into
(j_L, m_L) x (j_R, m_R), and two for the coupling with the photon
in each SU(2) factor.

Two-point structural zero
-------------------------
At q_probe=0 (no probe insertion), B = 0 by the Wigner-Eckart theorem:
the two-point vertex correction is m_j-independent since it is a
scalar (rank-0) operator on S^3. This is the SO(4) analog of the
statement that the self-energy renormalizes mass but not spin.

Parity selection
----------------
At q_probe=2 (even parity), B = 0 because the vertex parity rule
requires n_int + n_int + q_probe to be odd, and q_probe=2 forces
this to be even. Only odd q_probe modes have magnetic content.

References
----------
- Schwinger, Phys. Rev. 73 (1948) 416 [a_e = alpha/(2*pi)]
- Camporesi & Higuchi, J. Geom. Phys. 20 (1996) 1-18 [Dirac on S^3]
- GeoVac qed_self_energy.py (shell-summed vertex correction)
- GeoVac qed_vertex.py (SO(4) selection rules)
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Tuple

from sympy import Rational, S
from sympy.physics.wigner import clebsch_gordan

__all__ = [
    "half_integer_range",
    "vertex_amp_polarized",
    "vertex_allowed",
    "get_photon_channels",
    "tree_level_probe_magnetic",
    "tree_level_probe_coupling",
    "compute_vertex_3pt_single_level",
    "compute_anomalous_magnetic_moment",
    "anomalous_moment_convergence",
    "mode_dependent_analysis",
]

_ALPHA_CODATA = 7.2973525693e-3
_SCHWINGER = _ALPHA_CODATA / (2 * math.pi)


# ---------------------------------------------------------------------------
# Half-integer enumeration
# ---------------------------------------------------------------------------

def half_integer_range(j: Rational) -> List[Rational]:
    """Return all m values from -j to +j in integer steps.

    Parameters
    ----------
    j : Rational
        Angular momentum quantum number (integer or half-integer >= 0).

    Returns
    -------
    List of Rational m values from -j to +j.
    """
    result: List[Rational] = []
    m = -j
    while m <= j + Rational(1, 100):
        result.append(m)
        m += 1
    return result


# ---------------------------------------------------------------------------
# SO(4) vertex selection and channel enumeration
# ---------------------------------------------------------------------------

def vertex_allowed(n1: int, n2: int, q: int) -> bool:
    """Check SO(4) vertex selection rule.

    Requires triangle inequality and parity: n1 + n2 + q must be odd.

    Parameters
    ----------
    n1, n2 : int
        Dirac levels (CH convention, n >= 0).
    q : int
        Photon mode (q >= 1).
    """
    if q < 1 or q < abs(n1 - n2) or q > n1 + n2:
        return False
    return (n1 + n2 + q) % 2 == 1


def get_photon_channels(
    n_src: int,
    n_tgt: int,
    q: int,
    j_sL: Rational,
    j_sR: Rational,
    j_tL: Rational,
    j_tR: Rational,
) -> List[Tuple[Rational, Rational]]:
    """Return allowed (jgL, jgR) photon channels for a given vertex.

    The transverse photon at level q has two SO(4) components:
      ((q+1)/2, (q-1)/2) and ((q-1)/2, (q+1)/2).
    Each must satisfy SU(2) triangle inequalities in both L and R.

    Parameters
    ----------
    n_src, n_tgt : int
        Source and target Dirac levels (CH convention).
    q : int
        Photon level.
    j_sL, j_sR : Rational
        SU(2)_L x SU(2)_R quantum numbers of the source state.
    j_tL, j_tR : Rational
        SU(2)_L x SU(2)_R quantum numbers of the target state.

    Returns
    -------
    List of (jgL, jgR) pairs for allowed channels.
    """
    if not vertex_allowed(n_src, n_tgt, q):
        return []
    channels: List[Tuple[Rational, Rational]] = []
    for jgL, jgR in [(Rational(q + 1, 2), Rational(q - 1, 2)),
                      (Rational(q - 1, 2), Rational(q + 1, 2))]:
        if jgL < 0 or jgR < 0:
            continue
        if (abs(j_sL - jgL) <= j_tL <= j_sL + jgL
                and abs(j_sR - jgR) <= j_tR <= j_sR + jgR):
            channels.append((jgL, jgR))
    return channels


# ---------------------------------------------------------------------------
# Polarization-resolved CG vertex amplitude
# ---------------------------------------------------------------------------

def vertex_amp_polarized(
    j_sL: Rational, j_sR: Rational, j_s: Rational, mj_s: Rational,
    j_tL: Rational, j_tR: Rational, j_t: Rational, mj_t: Rational,
    jgL: Rational, jgR: Rational, mgL: Rational, mgR: Rational,
) -> Rational:
    """SO(4) vertex amplitude for a specific photon polarization.

    Computes the product of four CG coefficients:
      C1: source state decomposition  (j_sL, j_sR -> j_s, mj_s)
      C2: target state decomposition  (j_tL, j_tR -> j_t, mj_t)
      C3: left SU(2) coupling         (j_sL, jgL -> j_tL)
      C4: right SU(2) coupling        (j_sR, jgR -> j_tR)

    Summed over all compatible magnetic quantum number partitions
    (mL1, mR1) of the source state.

    Parameters
    ----------
    j_sL, j_sR, j_s, mj_s :
        Source state SU(2)_L x SU(2)_R reps and total j, m_j.
    j_tL, j_tR, j_t, mj_t :
        Target state SU(2)_L x SU(2)_R reps and total j, m_j.
    jgL, jgR, mgL, mgR :
        Photon polarization in SU(2)_L x SU(2)_R.

    Returns
    -------
    sympy Rational (exact)
    """
    total = S.Zero
    for mL1 in half_integer_range(j_sL):
        mR1 = mj_s - mL1
        if abs(mR1) > j_sR:
            continue
        mL2 = mL1 + mgL
        if abs(mL2) > j_tL:
            continue
        mR2 = mR1 + mgR
        if abs(mR2) > j_tR:
            continue
        if mL2 + mR2 != mj_t:
            continue
        c1 = clebsch_gordan(j_sL, j_sR, j_s, mL1, mR1, mj_s)
        if c1 == 0:
            continue
        c2 = clebsch_gordan(j_tL, j_tR, j_t, mL2, mR2, mj_t)
        if c2 == 0:
            continue
        c3 = clebsch_gordan(j_sL, jgL, j_tL, mL1, mgL, mL2)
        c4 = clebsch_gordan(j_sR, jgR, j_tR, mR1, mgR, mR2)
        total += c1 * c2 * c3 * c4
    return total


# ---------------------------------------------------------------------------
# Tree-level probe coupling
# ---------------------------------------------------------------------------

def tree_level_probe_coupling(
    n: int,
    j: Rational,
    mj: Rational,
    q_probe: int = 1,
) -> Rational:
    """Tree-level coupling of |n, j, m_j> to a probe photon of mode q.

    Sums over all probe polarizations and photon channels for the
    n -> n (same state) diagonal matrix element.

    Parameters
    ----------
    n : int
        Dirac level (CH convention).
    j : Rational
        Total angular momentum.
    mj : Rational
        Magnetic quantum number.
    q_probe : int
        Probe photon mode.

    Returns
    -------
    sympy Rational (exact)
    """
    jL = Rational(n + 1, 2)
    jR = Rational(n, 2)

    if not vertex_allowed(n, n, q_probe):
        return S.Zero

    channels = get_photon_channels(n, n, q_probe, jL, jR, jL, jR)
    if not channels:
        return S.Zero

    total = S.Zero
    for jpL, jpR in channels:
        for mpL in half_integer_range(jpL):
            for mpR in half_integer_range(jpR):
                total += vertex_amp_polarized(
                    jL, jR, j, mj,
                    jL, jR, j, mj,
                    jpL, jpR, mpL, mpR)
    return total


def tree_level_probe_magnetic(
    n: int,
    j: Rational,
    q_probe: int = 1,
) -> Rational:
    """Tree-level magnetic coupling: V_tree(+1/2) - V_tree(-1/2).

    The magnetic piece of the tree-level probe vertex is the m_j-odd
    part, which determines the normalization for extracting F_2.

    Parameters
    ----------
    n : int
        Dirac level (CH convention).
    j : Rational
        Total angular momentum.
    q_probe : int
        Probe photon mode.

    Returns
    -------
    sympy Rational (exact)
    """
    v_up = tree_level_probe_coupling(n, j, Rational(1, 2), q_probe)
    v_dn = tree_level_probe_coupling(n, j, Rational(-1, 2), q_probe)
    return v_up - v_dn


# ---------------------------------------------------------------------------
# Three-point vertex correction (per internal level)
# ---------------------------------------------------------------------------

def compute_vertex_3pt_single_level(
    n_ext: int,
    j_ext: Rational,
    mj_ext: Rational,
    n_int: int,
    q_probe: int = 1,
) -> Rational:
    """Contribution from a single internal level to the 3-pt vertex correction.

    The three-point vertex correction at internal level n_int inserts
    a probe photon (q_probe) between two loop vertices:

        ext --[loop photon]--> int --[probe]--> int --[loop photon]--> ext

    The probe insertion is summed over all probe polarizations.
    The loop photon is summed over all allowed modes q_loop.
    The internal electron propagator contributes 1/(lambda^4 * mu_q).

    Parameters
    ----------
    n_ext : int
        External electron level (CH convention).
    j_ext : Rational
        External total angular momentum.
    mj_ext : Rational
        External magnetic quantum number.
    n_int : int
        Internal electron level (CH convention).
    q_probe : int
        Probe photon mode.

    Returns
    -------
    sympy Rational (exact)
    """
    jE_L = Rational(n_ext + 1, 2)
    jE_R = Rational(n_ext, 2)
    jI_L = Rational(n_int + 1, 2)
    jI_R = Rational(n_int, 2)

    lam = Rational(2 * n_int + 3, 2)
    lam4 = lam ** 4

    j_int_min = abs(jI_L - jI_R)
    j_int_max = jI_L + jI_R

    if not vertex_allowed(n_int, n_int, q_probe):
        return S.Zero

    probe_channels = get_photon_channels(
        n_int, n_int, q_probe, jI_L, jI_R, jI_L, jI_R)
    if not probe_channels:
        return S.Zero

    total = S.Zero

    for j_int in half_integer_range(j_int_max):
        if j_int < j_int_min:
            continue
        for mj_int in half_integer_range(j_int):
            for mj_int_prime in half_integer_range(j_int):
                probe_amp = S.Zero
                for jpL, jpR in probe_channels:
                    for mpL in half_integer_range(jpL):
                        for mpR in half_integer_range(jpR):
                            probe_amp += vertex_amp_polarized(
                                jI_L, jI_R, j_int, mj_int,
                                jI_L, jI_R, j_int, mj_int_prime,
                                jpL, jpR, mpL, mpR)
                if probe_amp == 0:
                    continue

                q_lo = max(1, abs(n_ext - n_int))
                q_hi = n_ext + n_int
                for q_loop in range(q_lo, q_hi + 1):
                    if not vertex_allowed(n_ext, n_int, q_loop):
                        continue
                    mu_q = Rational(q_loop * (q_loop + 2))

                    chs1 = get_photon_channels(
                        n_ext, n_int, q_loop, jE_L, jE_R, jI_L, jI_R)
                    chs2 = get_photon_channels(
                        n_int, n_ext, q_loop, jI_L, jI_R, jE_L, jE_R)

                    for jgL1, jgR1 in chs1:
                        for jgL2, jgR2 in chs2:
                            for mgL in half_integer_range(jgL1):
                                for mgR in half_integer_range(jgR1):
                                    if abs(mgL) > jgL2 or abs(mgR) > jgR2:
                                        continue
                                    v1 = vertex_amp_polarized(
                                        jE_L, jE_R, j_ext, mj_ext,
                                        jI_L, jI_R, j_int, mj_int,
                                        jgL1, jgR1, mgL, mgR)
                                    if v1 == 0:
                                        continue
                                    v3 = vertex_amp_polarized(
                                        jI_L, jI_R, j_int, mj_int_prime,
                                        jE_L, jE_R, j_ext, mj_ext,
                                        jgL2, jgR2, mgL, mgR)
                                    if v3 == 0:
                                        continue
                                    total += v1 * probe_amp * v3 / (lam4 * mu_q)

    return total


# ---------------------------------------------------------------------------
# Anomalous magnetic moment extraction
# ---------------------------------------------------------------------------

def compute_anomalous_magnetic_moment(
    n_ext: int = 1,
    n_max: int = 3,
    q_probe: int = 1,
) -> Dict[str, object]:
    """Compute the anomalous magnetic moment F_2 on S^3.

    Sums the three-point vertex correction over internal levels
    0..n_max for both m_j = +1/2 and -1/2, extracts the magnetic
    difference B, and normalizes by the tree-level magnetic coupling
    to get F_2.

    Parameters
    ----------
    n_ext : int
        External electron level (CH convention, n >= 0).
    n_max : int
        Maximum internal level to sum over.
    q_probe : int
        Probe photon mode (default 1; q_probe=2 gives zero by parity).

    Returns
    -------
    Dict with keys:
      - F2: the anomalous magnetic moment (float)
      - F2_over_schwinger: F2 / [alpha/(2*pi)] (float)
      - B_magnetic: the magnetic vertex difference (float)
      - V_tree_magnetic: tree-level magnetic coupling (float)
      - per_level: list of per-n_int contributions
    """
    j_ext = Rational(1, 2)

    v_mag = tree_level_probe_magnetic(n_ext, j_ext, q_probe)
    v_mag_float = float(v_mag)

    per_level: List[Dict[str, object]] = []
    b_cumul = S.Zero

    for n_int in range(n_max + 1):
        contrib_up = compute_vertex_3pt_single_level(
            n_ext, j_ext, Rational(1, 2), n_int, q_probe)
        contrib_dn = compute_vertex_3pt_single_level(
            n_ext, j_ext, Rational(-1, 2), n_int, q_probe)

        b_level = contrib_up - contrib_dn
        b_cumul += b_level

        b_level_f = float(b_level)
        b_cumul_f = float(b_cumul)
        f2 = b_cumul_f / v_mag_float if v_mag_float != 0 else 0.0
        ratio = f2 / _SCHWINGER if _SCHWINGER != 0 else 0.0

        per_level.append({
            "n_int": n_int,
            "B_level": b_level_f,
            "B_cumul": b_cumul_f,
            "F2": f2,
            "F2_over_schwinger": ratio,
        })

    b_total = float(b_cumul)
    f2_total = b_total / v_mag_float if v_mag_float != 0 else 0.0

    return {
        "F2": f2_total,
        "F2_over_schwinger": f2_total / _SCHWINGER if _SCHWINGER != 0 else 0.0,
        "B_magnetic": b_total,
        "V_tree_magnetic": v_mag_float,
        "n_ext": n_ext,
        "n_max": n_max,
        "q_probe": q_probe,
        "per_level": per_level,
    }


def anomalous_moment_convergence(
    n_max_values: List[int],
    n_ext: int = 1,
    q_probe: int = 1,
) -> List[Dict[str, object]]:
    """Convergence study: F_2/Schwinger vs n_max.

    Parameters
    ----------
    n_max_values : List[int]
        Cutoff values (should be sorted ascending).
    n_ext : int
        External electron level.
    q_probe : int
        Probe photon mode.

    Returns
    -------
    List of dicts with n_max, F2, F2_over_schwinger, delta.
    """
    results: List[Dict[str, object]] = []
    prev_f2 = None

    for n_max in n_max_values:
        result = compute_anomalous_magnetic_moment(n_ext, n_max, q_probe)
        f2 = result["F2"]
        delta = abs(f2 - prev_f2) if prev_f2 is not None else None
        results.append({
            "n_max": n_max,
            "F2": f2,
            "F2_over_schwinger": result["F2_over_schwinger"],
            "B_magnetic": result["B_magnetic"],
            "delta": delta,
        })
        prev_f2 = f2

    return results


# ---------------------------------------------------------------------------
# Mode-dependent analysis across n_ext values
# ---------------------------------------------------------------------------

def mode_dependent_analysis(
    n_ext_values: Optional[List[int]] = None,
    n_max: int = 3,
    q_probe: int = 1,
) -> Dict[str, object]:
    """Compute F_2/Schwinger at multiple n_ext values and fit a power law.

    The anomalous magnetic moment depends on which Dirac mode is the
    external electron. Higher modes (larger n_ext) are more suppressed
    because the Dirac eigenvalue lambda = n_ext + 3/2 grows, and the
    curvature correction scales as a power law in 1/lambda.

    This function computes F_2/Schwinger for each n_ext, then fits
    the curvature correction (F_2/Schwinger - 1) to a power law
    in lambda = (2*n_ext + 3)/2:

        F_2/Schwinger = 1 + c1/lambda^2 + c2/lambda^4 + ...

    Parameters
    ----------
    n_ext_values : List[int] or None
        External electron levels (CH convention, n >= 1).
        Defaults to [1, 2, 3].
    n_max : int
        Maximum internal level for the mode sum.
    q_probe : int
        Probe photon mode.

    Returns
    -------
    Dict with keys:
      - per_mode: list of per-n_ext results (n_ext, lambda_ext,
        V_magnetic, B_magnetic, F2, F2_over_schwinger, correction)
      - ratios: list of F2/Schwinger values (same order as n_ext_values)
      - monotone_decreasing: bool, whether F2/Schwinger decreases with n_ext
      - fit_coeffs: dict with c1, c2 (if >= 2 data points)
    """
    if n_ext_values is None:
        n_ext_values = [1, 2, 3]

    per_mode: List[Dict[str, object]] = []
    ratios: List[float] = []

    for n_ext in n_ext_values:
        result = compute_anomalous_magnetic_moment(n_ext, n_max, q_probe)
        lam_ext = (2 * n_ext + 3) / 2.0
        ratio = result["F2_over_schwinger"]
        correction = ratio - 1.0

        per_mode.append({
            "n_ext": n_ext,
            "lambda_ext": lam_ext,
            "V_magnetic": result["V_tree_magnetic"],
            "B_magnetic": result["B_magnetic"],
            "F2": result["F2"],
            "F2_over_schwinger": ratio,
            "correction": correction,
        })
        ratios.append(ratio)

    # Check monotone-decreasing property
    monotone = all(
        ratios[i] > ratios[i + 1]
        for i in range(len(ratios) - 1)
    )

    output: Dict[str, object] = {
        "per_mode": per_mode,
        "ratios": ratios,
        "monotone_decreasing": monotone,
        "n_max": n_max,
        "q_probe": q_probe,
    }

    # Fit curvature expansion if we have enough points
    n_pts = len(per_mode)
    if n_pts >= 2:
        lam_vals = [m["lambda_ext"] for m in per_mode]
        corr_vals = [m["correction"] for m in per_mode]

        # Least-squares fit: correction = c1/lam^2 + c2/lam^4
        # Build design matrix A
        a_rows = [[1.0 / lam ** 2, 1.0 / lam ** 4] for lam in lam_vals]

        if n_pts == 2:
            # Exact solve for 2x2 system
            a00, a01 = a_rows[0]
            a10, a11 = a_rows[1]
            det = a00 * a11 - a01 * a10
            if abs(det) > 1e-30:
                c1 = (a11 * corr_vals[0] - a01 * corr_vals[1]) / det
                c2 = (-a10 * corr_vals[0] + a00 * corr_vals[1]) / det
            else:
                c1 = 0.0
                c2 = 0.0
        else:
            # Least-squares via normal equations (no numpy dependency)
            # A^T A x = A^T b
            ata00 = sum(r[0] ** 2 for r in a_rows)
            ata01 = sum(r[0] * r[1] for r in a_rows)
            ata11 = sum(r[1] ** 2 for r in a_rows)
            atb0 = sum(r[0] * c for r, c in zip(a_rows, corr_vals))
            atb1 = sum(r[1] * c for r, c in zip(a_rows, corr_vals))
            det = ata00 * ata11 - ata01 ** 2
            if abs(det) > 1e-30:
                c1 = (ata11 * atb0 - ata01 * atb1) / det
                c2 = (-ata01 * atb0 + ata00 * atb1) / det
            else:
                c1 = 0.0
                c2 = 0.0

        output["fit_coeffs"] = {"c1": c1, "c2": c2}

    return output
