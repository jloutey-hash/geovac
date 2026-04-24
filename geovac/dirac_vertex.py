"""State-resolved QED vertex correction on the native Dirac graph.

Computes the one-loop vertex correction (the diagram that gives g-2)
by summing over individual Dirac states |n, κ, m_j⟩ on S³, rather
than degenerate shells with g(n) = 2(n+1)(n+2).

The vertex correction for an external state (n_ext, κ_ext, m_j_ext) is:

    Λ(ext) = Σ_{int, q1, q2} C(ext→int, q1) · C(int→ext, q2)
             × 1/(|λ_int|^{2s_e} · μ_{q1}^{s_γ} · μ_{q2}^{s_γ})

where the sum runs over ALL internal Dirac states (n_int, κ_int, m_j_int)
and photon modes q1, q2.

At the SO(4) level (before resolving m_j-dependent CG coefficients),
each state in a shell contributes equally, so the state-resolved sum
must equal the shell-summed version. This is the structural validation.

The m_j-dependent piece — which IS the anomalous magnetic moment —
enters through the angular CG coefficients C(κ_ext, m_j_ext; κ_int, m_j_int; q).
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np
import mpmath

from geovac.dirac_lattice import DiracLattice
from geovac.dirac_matrix_elements import DiracLabel, kappa_to_l

mpmath.mp.dps = 50


def _lambda_n_ch(n_ch: int) -> mpmath.mpf:
    """Camporesi-Higuchi |λ| = n_CH + 3/2."""
    return mpmath.mpf(n_ch) + mpmath.mpf(3) / 2


def _mu_q(q: int) -> mpmath.mpf:
    """Hodge-1 eigenvalue μ_q = q(q+2)."""
    return mpmath.mpf(q) * (q + 2)


def _d_T(q: int) -> mpmath.mpf:
    """Transverse photon degeneracy d_q^T = q(q+2)."""
    return mpmath.mpf(q) * (q + 2)


def _vertex_allowed(n1: int, n2: int, q: int) -> bool:
    """SO(4) vertex selection rule at the shell level."""
    if q < 1:
        return False
    if q < abs(n1 - n2) or q > n1 + n2:
        return False
    if (n1 + n2 + q) % 2 == 0:
        return False
    return True


def _so4_channel_count(n1: int, n2: int, q: int) -> int:
    """Number of SO(4) vector-harmonic components coupling n1->n2 via photon q.

    Checks both the SU(2)_L x SU(2)_R triangle inequalities AND the
    CG parity constraint (n1 + n2 + q must be odd) via _vertex_allowed.
    Returns 0, 1, or 2.
    """
    if not _vertex_allowed(n1, n2, q):
        return 0
    from fractions import Fraction
    j1_L = Fraction(n1 + 1, 2)
    j1_R = Fraction(n1, 2)
    j2_L = Fraction(n2, 2)
    j2_R = Fraction(n2 + 1, 2)
    count = 0
    jg_L_A = Fraction(q + 1, 2)
    jg_R_A = Fraction(q - 1, 2)
    if (jg_R_A >= 0
            and abs(j1_L - jg_L_A) <= j2_L <= j1_L + jg_L_A
            and abs(j1_R - jg_R_A) <= j2_R <= j1_R + jg_R_A):
        count += 1
    jg_L_B = Fraction(q - 1, 2)
    jg_R_B = Fraction(q + 1, 2)
    if (jg_L_B >= 0
            and abs(j1_L - jg_L_B) <= j2_L <= j1_L + jg_L_B
            and abs(j1_R - jg_R_B) <= j2_R <= j1_R + jg_R_B):
        count += 1
    return count


def vertex_correction_state_resolved(
    lat: DiracLattice,
    ext_idx: int,
    q_max: int = 0,
    s_e: int = 2,
    s_gamma: int = 1,
) -> mpmath.mpf:
    """State-resolved one-loop vertex correction for a specific external state.

    Sums over individual internal Dirac states (not degenerate shells).
    At the SO(4) level, all states within a shell contribute equally,
    so summing g(n_int) individual states reproduces the g(n_int) factor
    in the shell-summed formula. The result equals the shell-summed
    vertex correction (which is already per-external-state by SO(4) symmetry).

    Parameters
    ----------
    lat : DiracLattice
        The native Dirac graph.
    ext_idx : int
        Index into lat.labels for the external state.
    q_max : int
        Maximum photon mode. If 0, uses 2*n_max + 1.
    s_e, s_gamma : int
        Propagator exponents.

    Returns
    -------
    mpmath.mpf
        The vertex correction sum.
    """
    ext = lat.labels[ext_idx]
    n_ext_ch = ext.n_fock - 1  # Fock → CH convention
    n_max_ch = lat.n_max - 1
    if q_max <= 0:
        q_max = 2 * n_max_ch + 1

    s_eff = 2 * s_e
    total = mpmath.mpf(0)

    for int_idx, int_lab in enumerate(lat.labels):
        n_int_ch = int_lab.n_fock - 1
        lam_int = _lambda_n_ch(n_int_ch)
        lam_pow = lam_int ** s_eff

        q1_lo = abs(n_ext_ch - n_int_ch)
        q1_hi = min(n_ext_ch + n_int_ch, q_max)

        for q1 in range(max(1, q1_lo), q1_hi + 1):
            W1 = _so4_channel_count(n_ext_ch, n_int_ch, q1)
            if W1 == 0:
                continue

            q2_lo = abs(n_int_ch - n_ext_ch)
            q2_hi = min(n_int_ch + n_ext_ch, q_max)

            for q2 in range(max(1, q2_lo), q2_hi + 1):
                W2 = _so4_channel_count(n_int_ch, n_ext_ch, q2)
                if W2 == 0:
                    continue

                total += (mpmath.mpf(W1) * W2
                          * _d_T(q1) * _d_T(q2)
                          / (lam_pow
                             * _mu_q(q1) ** s_gamma
                             * _mu_q(q2) ** s_gamma))

    return total


def vertex_correction_shell_summed(
    n_ext_ch: int,
    n_max_ch: int,
    s_e: int = 2,
    s_gamma: int = 1,
) -> mpmath.mpf:
    """Shell-summed vertex correction (reproduces qed_self_energy version).

    For cross-validation with the state-resolved version.
    """
    s_eff = 2 * s_e
    total = mpmath.mpf(0)

    for n_int in range(n_max_ch + 1):
        g_int = mpmath.mpf(2) * (n_int + 1) * (n_int + 2)
        lam_int = _lambda_n_ch(n_int)
        lam_pow = lam_int ** s_eff

        q1_lo = abs(n_ext_ch - n_int)
        q1_hi = n_ext_ch + n_int

        for q1 in range(max(1, q1_lo), q1_hi + 1):
            W1 = _so4_channel_count(n_ext_ch, n_int, q1)
            if W1 == 0:
                continue

            q2_lo = abs(n_int - n_ext_ch)
            q2_hi = n_int + n_ext_ch

            for q2 in range(max(1, q2_lo), q2_hi + 1):
                W2 = _so4_channel_count(n_int, n_ext_ch, q2)
                if W2 == 0:
                    continue

                total += (mpmath.mpf(W1) * W2
                          * g_int * _d_T(q1) * _d_T(q2)
                          / (lam_pow
                             * _mu_q(q1) ** s_gamma
                             * _mu_q(q2) ** s_gamma))

    return total


def validate_state_resolved(n_max: int = 3, mode: str = 's3') -> Dict[str, object]:
    """Validate that each state-resolved value equals shell-summed.

    shell_summed gives the per-external-state vertex correction (it
    includes g_int for the internal sum, which is per-state by SO(4)).
    state_resolved sums over individual internal states and should
    reproduce this exactly for every external state.

    Returns a dict with per-state results and max relative error.
    """
    lat = DiracLattice(n_max=n_max, mode=mode)
    n_max_ch = n_max - 1

    results = []
    max_rel_err = 0.0

    unique_n_ext = sorted(set(lab.n_fock for lab in lat.labels))

    for n_fock in unique_n_ext:
        n_ext_ch = n_fock - 1
        shell_val = vertex_correction_shell_summed(n_ext_ch, n_max_ch)

        ext_indices = [i for i, lab in enumerate(lat.labels) if lab.n_fock == n_fock]
        state_vals = []
        for idx in ext_indices:
            val = vertex_correction_state_resolved(lat, idx)
            state_vals.append(float(val))

        per_state_mean = float(np.mean(state_vals)) if state_vals else 0
        per_state_std = float(np.std(state_vals)) if len(state_vals) > 1 else 0

        if abs(shell_val) > 1e-50:
            rel_err = float(abs(per_state_mean - float(shell_val)) / abs(float(shell_val)))
        else:
            rel_err = float(abs(per_state_mean))

        max_rel_err = max(max_rel_err, rel_err)

        results.append({
            "n_fock": n_fock,
            "shell_summed": float(shell_val),
            "per_state_value": per_state_mean,
            "num_states_in_shell": len(ext_indices),
            "per_state_std": per_state_std,
            "rel_error": rel_err,
        })

    return {
        "n_max": n_max,
        "mode": mode,
        "shells": results,
        "max_rel_error": max_rel_err,
        "structural_match": max_rel_err < 1e-10,
    }


def schwinger_convergence_state_resolved(
    n_max_range: List[int],
    n_ext_fock: int = 2,
    s_e: int = 2,
    s_gamma: int = 1,
) -> List[Dict[str, object]]:
    """Track convergence of the state-resolved vertex correction.

    Builds a DiracLattice(mode='s3') at each cutoff and computes the
    vertex correction for a representative external state at n_ext_fock.
    Compares to the shell-summed value and to the Schwinger target α/(2π).

    Parameters
    ----------
    n_max_range : list of int
        Cutoff values (sorted ascending).
    n_ext_fock : int
        External state's Fock quantum number (n_fock ≥ 2 for non-zero).
    s_e, s_gamma : int
        Propagator exponents.
    """
    alpha = 7.2973525693e-3
    schwinger = alpha / (2 * np.pi)
    results = []
    prev_val = None

    for n_max in n_max_range:
        if n_max < n_ext_fock:
            continue

        lat = DiracLattice(n_max=n_max, mode='s3')
        n_ext_ch = n_ext_fock - 1

        shell_val = vertex_correction_shell_summed(n_ext_ch, n_max - 1, s_e=s_e, s_gamma=s_gamma)

        ext_indices = [i for i, lab in enumerate(lat.labels) if lab.n_fock == n_ext_fock]
        state_val = vertex_correction_state_resolved(lat, ext_indices[0], s_e=s_e, s_gamma=s_gamma)

        delta = float(abs(state_val - prev_val)) if prev_val is not None else None

        results.append({
            "n_max": n_max,
            "n_states": lat.num_states,
            "state_resolved": float(state_val),
            "shell_summed": float(shell_val),
            "match_rel_err": float(abs(state_val - shell_val) / abs(shell_val)) if abs(shell_val) > 1e-50 else 0,
            "ratio_to_schwinger": float(state_val) / schwinger if schwinger > 0 else 0,
            "delta": delta,
        })
        prev_val = state_val

    return results
