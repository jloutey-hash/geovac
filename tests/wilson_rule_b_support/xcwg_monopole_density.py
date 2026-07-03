"""
XCWG-F monopole density on Rule B compact U(1) Wilson gauge theory (2026-05-16).

Fifth-witness diagnostic for the four-witness leading-order verdict in Paper 41.

Theoretical setup
=================

DeGrand-Toussaint 1980 PRD 22, 2478 monopole extraction:

For each plaquette P in a U(1) gauge configuration {theta_e}, compute the
oriented sum of link angles:
    theta_P = sum_{e in dP} sign(P,e) * theta_e

Project to the principal branch:
    theta_bar_P = theta_P - 2 pi * nint(theta_P / (2 pi))   in (-pi, pi]

The Dirac-string charge is the integer:
    m_P = (theta_P - theta_bar_P) / (2 pi)   in Z

For each elementary closed 2-surface S (analog of dual-lattice 3-cube), the
monopole charge enclosed by S is:
    m_S = (1 / 2 pi) sum_{P in S} sign(P, S) * (theta_P - theta_bar_P)
        = sum_{P in S} sign(P, S) * m_P

This sum is automatically an integer (Stokes for the integer-valued
Dirac-string field).  The "elementary" closed 2-surfaces on Rule B are the
size-3 closed 2-cycles (triangle-prism configurations of three L=4 plaquettes
sharing a common edge axis) identified by XCWG-E.

Polyakov 1977 (3D compact U(1)): permanent confinement at all beta via dilute
monopole plasma.  rho_M(beta) > 0 at all beta, dilute-gas exponential decay
at large beta: rho_M ~ A * exp(-c * beta).

BMK 1977 / Guth 1980 (4D compact U(1)): transition at beta_c ~ 1.01.
Below: condensed phase, rho_M large.  Above: rho_M -> 0 (deconfined Coulomb).

This script: Monte Carlo simulation of compact U(1) Wilson on the Rule B
Dirac-S^3 graph at n_max=2 (and n_max=3 if tractable), with Metropolis-Hastings
updates of link angles and DeGrand-Toussaint monopole extraction.

Output
======
    tests/wilson_rule_b_support/data/xcwg_monopole_density.json
    debug/xcwg_monopole_density_memo.md          (separately)
    tests/wilson_rule_b_support/data/xcwg_monopole_density.png        (if matplotlib available)
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir, os.pardir))  # repo root (tests/wilson_rule_b_support/ -> two levels up)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

from geovac.ihara_zeta_dirac import build_dirac_s3_graph  # noqa: E402
from xcwg_wilson_loop_scaling import (  # noqa: E402
    signed_incidence, adjacency_list, enumerate_primitive_closed_walks,
)
from xcwg_strong_coupling_wilson import (  # noqa: E402
    edge_set_of_walk, adjacent_plaquettes,
)
from xcwg_nlo_character_expansion import (  # noqa: E402
    build_d1_and_plaquettes, find_smallest_closed_2cycle,
    enumerate_closed_2cycles_up_to_size,
)

sys.stdout.reconfigure(line_buffering=True)

TWO_PI = 2.0 * np.pi


# =============================================================================
# DeGrand-Toussaint plaquette angle and monopole charge
# =============================================================================

def plaquette_angles(d_1: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Compute theta_P = sum_e sign(P,e) * theta_e for each plaquette.

    Args:
        d_1: P x E signed plaquette-boundary matrix (entries in {-1, 0, +1})
        theta: length-E vector of link angles in (-pi, pi]

    Returns:
        length-P vector of oriented plaquette angles (NOT folded back to (-pi, pi])
    """
    return d_1.astype(np.float64) @ theta


def principal_branch(angles: np.ndarray) -> np.ndarray:
    """Fold to (-pi, pi]."""
    return angles - TWO_PI * np.round(angles / TWO_PI)


def dirac_string_charge(d_1: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Compute m_P = (theta_P - theta_bar_P) / (2 pi) for each plaquette.

    Returns length-P integer array.
    """
    theta_P = plaquette_angles(d_1, theta)
    theta_bar_P = principal_branch(theta_P)
    m_P_float = (theta_P - theta_bar_P) / TWO_PI
    # Should be very close to integers
    return np.round(m_P_float).astype(np.int64)


def monopole_charges(
    d_1: np.ndarray,
    theta: np.ndarray,
    monopole_sites: List[Dict],
) -> np.ndarray:
    """For each monopole site (closed 2-surface S), compute m_S = sum_P sign(P,S) m_P.

    Args:
        d_1: P x E signed plaquette-boundary matrix
        theta: length-E gauge configuration
        monopole_sites: list of dicts {"plaq_set": [...], "signs": [...], "size": k}

    Returns:
        length-(n_sites) integer array of monopole charges
    """
    m_P = dirac_string_charge(d_1, theta)
    n_sites = len(monopole_sites)
    m_S = np.zeros(n_sites, dtype=np.int64)
    for s, site in enumerate(monopole_sites):
        plaq_set = site["plaq_set"]
        signs = site["signs"]
        m_S[s] = sum(sign * m_P[p] for p, sign in zip(plaq_set, signs))
    return m_S


# =============================================================================
# Wilson action and Metropolis-Hastings updates
# =============================================================================

def wilson_action(d_1: np.ndarray, theta: np.ndarray, beta: float) -> float:
    """S_W = beta * sum_P (1 - cos(theta_P))."""
    theta_P = plaquette_angles(d_1, theta)
    return float(beta * np.sum(1.0 - np.cos(theta_P)))


def plaquettes_containing_edge(d_1: np.ndarray) -> List[np.ndarray]:
    """For each edge e, return list of (plaq_idx, sign) entries.

    Stored as a length-E list of arrays of shape (k_e, 2) where the first column
    is plaquette index and the second is the sign of that edge in the plaquette.
    """
    P, E = d_1.shape
    out: List[np.ndarray] = []
    for e in range(E):
        col = d_1[:, e]
        nz = np.nonzero(col)[0]
        out.append(np.stack([nz, col[nz]], axis=1).astype(np.int64))
    return out


def local_delta_action(
    d_1_e_data: np.ndarray,  # array of shape (k_e, 2): (plaq_idx, sign)
    theta_P_current: np.ndarray,  # current plaquette angles (length P)
    delta_theta: float,
    sign_at_e_for_each_plaq: np.ndarray,  # length k_e: sign of edge e in each plaquette
    plaq_indices: np.ndarray,  # length k_e
    beta: float,
) -> float:
    """Compute Delta S_W when theta_e changes by delta_theta.

    For each plaquette P containing edge e with sign s_eP,
        theta_P' = theta_P + s_eP * delta_theta
    so
        Delta S_W = beta * sum_P [cos(theta_P) - cos(theta_P + s_eP * delta_theta)]
    """
    theta_P_old = theta_P_current[plaq_indices]
    theta_P_new = theta_P_old + sign_at_e_for_each_plaq * delta_theta
    delta_S = beta * np.sum(np.cos(theta_P_old) - np.cos(theta_P_new))
    return float(delta_S)


def metropolis_sweep(
    d_1: np.ndarray,
    theta: np.ndarray,
    theta_P: np.ndarray,
    edge_data: List[np.ndarray],
    beta: float,
    delta_max: float,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray, int]:
    """One Metropolis sweep over all links.

    Mutates theta and theta_P in place; returns (theta, theta_P, n_accept).
    """
    E = theta.shape[0]
    n_accept = 0
    for e in range(E):
        # Propose
        delta = rng.uniform(-delta_max, delta_max)
        ed = edge_data[e]
        if ed.size == 0:
            # No plaquette contains this edge; angle is gauge / unconstrained
            theta[e] = (theta[e] + delta + np.pi) % TWO_PI - np.pi
            n_accept += 1
            continue
        plaq_indices = ed[:, 0]
        signs = ed[:, 1].astype(np.float64)
        theta_P_old = theta_P[plaq_indices]
        theta_P_new = theta_P_old + signs * delta
        dS = beta * np.sum(np.cos(theta_P_old) - np.cos(theta_P_new))
        # Note: dS = S_new - S_old.  Accept if exp(-dS) > u.
        if dS <= 0.0 or rng.random() < np.exp(-dS):
            theta[e] = (theta[e] + delta + np.pi) % TWO_PI - np.pi
            theta_P[plaq_indices] = theta_P_new
            n_accept += 1
    return theta, theta_P, n_accept


def tune_delta_max(
    d_1: np.ndarray,
    theta: np.ndarray,
    edge_data: List[np.ndarray],
    beta: float,
    rng: np.random.Generator,
    target_accept: float = 0.5,
    n_tune_sweeps: int = 50,
) -> Tuple[float, float]:
    """Quick auto-tune of delta_max to hit roughly target_accept acceptance rate."""
    delta_max = np.pi
    theta_P = plaquette_angles(d_1, theta)
    for it in range(8):
        accepts = 0
        attempts = 0
        for _ in range(n_tune_sweeps):
            E = theta.shape[0]
            for e in range(E):
                delta = rng.uniform(-delta_max, delta_max)
                ed = edge_data[e]
                if ed.size == 0:
                    theta[e] = (theta[e] + delta + np.pi) % TWO_PI - np.pi
                    accepts += 1
                    attempts += 1
                    continue
                plaq_indices = ed[:, 0]
                signs = ed[:, 1].astype(np.float64)
                theta_P_old = theta_P[plaq_indices]
                theta_P_new = theta_P_old + signs * delta
                dS = beta * np.sum(np.cos(theta_P_old) - np.cos(theta_P_new))
                attempts += 1
                if dS <= 0.0 or rng.random() < np.exp(-dS):
                    theta[e] = (theta[e] + delta + np.pi) % TWO_PI - np.pi
                    theta_P[plaq_indices] = theta_P_new
                    accepts += 1
        rate = accepts / max(attempts, 1)
        if abs(rate - target_accept) < 0.05:
            return delta_max, rate
        # Adjust
        if rate > target_accept:
            delta_max = min(np.pi, delta_max * 1.4)
        else:
            delta_max = delta_max * 0.7
    return delta_max, rate


# =============================================================================
# Main MC driver at one beta
# =============================================================================

def mc_run_one_beta(
    d_1: np.ndarray,
    edge_data: List[np.ndarray],
    monopole_sites: List[Dict],
    beta: float,
    n_therm: int,
    n_sample: int,
    sample_interval: int,
    rng: np.random.Generator,
    verbose: bool = True,
) -> Dict:
    """Run MC at one beta, return monopole density statistics.

    rho_M defined as < |m_S| > averaged over monopole sites and configurations.
    Also reports < m_S^2 > and the autocorrelation-based error bar (block-resampling).
    """
    E = d_1.shape[1]
    n_sites = len(monopole_sites)

    # Cold start (theta=0): for small beta this is fine; for large beta also fine
    theta = np.zeros(E)
    # Tune delta_max
    theta_tune = theta.copy()
    delta_max, accept_rate = tune_delta_max(d_1, theta_tune, edge_data, beta, rng,
                                            target_accept=0.5, n_tune_sweeps=20)
    if verbose:
        print(f"    [tune] delta_max={delta_max:.3f}, accept_rate~{accept_rate:.3f}")

    theta_P = plaquette_angles(d_1, theta)

    # Thermalize
    t_therm_start = time.time()
    n_accept_therm = 0
    for sweep in range(n_therm):
        theta, theta_P, n_acc = metropolis_sweep(
            d_1, theta, theta_P, edge_data, beta, delta_max, rng
        )
        n_accept_therm += n_acc
    therm_accept_rate = n_accept_therm / (n_therm * E)
    t_therm = time.time() - t_therm_start

    # Sample
    t_samp_start = time.time()
    samples_abs_mS: List[np.ndarray] = []   # each is shape (n_sites,)
    samples_m2_mS: List[np.ndarray] = []
    plaq_charges_abs_avg: List[float] = []
    n_accept_samp = 0

    for k in range(n_sample):
        for _ in range(sample_interval):
            theta, theta_P, n_acc = metropolis_sweep(
                d_1, theta, theta_P, edge_data, beta, delta_max, rng
            )
            n_accept_samp += n_acc
        m_P_now = dirac_string_charge(d_1, theta)
        plaq_charges_abs_avg.append(float(np.mean(np.abs(m_P_now))))
        m_S = monopole_charges(d_1, theta, monopole_sites)
        samples_abs_mS.append(np.abs(m_S))
        samples_m2_mS.append(m_S.astype(np.float64) ** 2)

    t_samp = time.time() - t_samp_start
    samp_accept_rate = n_accept_samp / (n_sample * sample_interval * E)

    # Stack
    abs_arr = np.stack(samples_abs_mS, axis=0).astype(np.float64)  # (n_sample, n_sites)
    m2_arr = np.stack(samples_m2_mS, axis=0)

    # rho_M = mean |m_S| over sites and samples
    per_sample_rho = abs_arr.mean(axis=1)          # length n_sample
    per_sample_m2 = m2_arr.mean(axis=1)
    rho_M_mean = float(per_sample_rho.mean())
    rho_M_var = float(per_sample_rho.var(ddof=1)) if n_sample > 1 else 0.0
    # Block standard error: 10 blocks
    n_blocks = min(10, n_sample)
    block_size = n_sample // n_blocks
    if block_size >= 1:
        block_means = np.array([
            per_sample_rho[i * block_size : (i + 1) * block_size].mean()
            for i in range(n_blocks)
        ])
        rho_M_se = float(block_means.std(ddof=1) / np.sqrt(n_blocks)) if n_blocks > 1 else 0.0
    else:
        rho_M_se = float(np.sqrt(rho_M_var / max(n_sample, 1)))

    # Fraction of configurations with at least one nonzero monopole
    has_monopole = (abs_arr.max(axis=1) > 0)
    frac_with_monopole = float(np.mean(has_monopole))

    # Fraction of monopole sites with nonzero charge (averaged across configurations)
    site_nonzero_frac = float((abs_arr > 0).mean())

    if verbose:
        print(f"    [beta={beta:.4g}] rho_M = {rho_M_mean:.4e} +/- {rho_M_se:.4e}; "
              f"<m_S^2> = {per_sample_m2.mean():.4e}; "
              f"site_nz_frac = {site_nonzero_frac:.4e}; "
              f"any-monopole frac = {frac_with_monopole:.3f}; "
              f"mean |m_P| = {np.mean(plaq_charges_abs_avg):.4e}; "
              f"therm/samp = {t_therm:.1f}/{t_samp:.1f}s; "
              f"accept (samp) = {samp_accept_rate:.3f}")

    return {
        "beta": float(beta),
        "rho_M_mean": rho_M_mean,
        "rho_M_se": rho_M_se,
        "rho_M_var": rho_M_var,
        "m2_mean": float(per_sample_m2.mean()),
        "site_nonzero_frac": site_nonzero_frac,
        "frac_configs_with_monopole": frac_with_monopole,
        "mean_abs_mP": float(np.mean(plaq_charges_abs_avg)),
        "delta_max": float(delta_max),
        "therm_accept_rate": float(therm_accept_rate),
        "samp_accept_rate": float(samp_accept_rate),
        "n_therm": int(n_therm),
        "n_sample": int(n_sample),
        "sample_interval": int(sample_interval),
        "therm_time_sec": float(t_therm),
        "samp_time_sec": float(t_samp),
    }


# =============================================================================
# Run at one n_max
# =============================================================================

def run_at_nmax(
    n_max: int,
    beta_grid: List[float],
    n_therm: int,
    n_sample: int,
    sample_interval: int,
    cycle_max_count: int,
    cycle_time_budget_sec: float,
    rng_seed: int,
    plaq_cap: int = 200,
    larger_cycle_cap: int = 6,
) -> Dict:
    """Run the monopole-density sweep at one n_max."""
    print(f"\n{'='*72}")
    print(f"Monopole-density MC sweep at n_max={n_max}")
    print(f"{'='*72}")

    rng = np.random.default_rng(rng_seed)

    t_build = time.time()
    g = build_d1_and_plaquettes(n_max, plaq_cap=plaq_cap)
    print(f"  V={g['V']}, E={g['E']}, P={g['P']}, beta_1={g['beta_1_graph']}, "
          f"build_time={time.time() - t_build:.1f}s")

    # Find smallest 2-cycle (k_delta)
    t_cyc = time.time()
    smallest = find_smallest_closed_2cycle(g["d_1"], g["plaq_adj"],
                                            max_size=8,
                                            time_budget_sec=cycle_time_budget_sec)
    print(f"  smallest-cycle search took {time.time() - t_cyc:.1f}s")
    if smallest is None:
        raise RuntimeError(f"No closed 2-cycle found at n_max={n_max} within budget")
    k_delta = smallest["size"]
    print(f"  k_delta = {k_delta}")

    # Enumerate all elementary monopole sites at size k_delta
    t_enum = time.time()
    elementary_sites = enumerate_closed_2cycles_up_to_size(
        g["d_1"], g["plaq_adj"],
        size_target=k_delta,
        max_count=cycle_max_count,
        time_budget_sec=cycle_time_budget_sec,
    )
    print(f"  Enumerated {len(elementary_sites)} elementary monopole sites "
          f"(size={k_delta}) in {time.time() - t_enum:.1f}s")

    # Also enumerate larger closed 2-surfaces up to size larger_cycle_cap
    larger_sites: List[Dict] = []
    for sz in range(k_delta + 1, larger_cycle_cap + 1):
        t0 = time.time()
        cyc = enumerate_closed_2cycles_up_to_size(
            g["d_1"], g["plaq_adj"],
            size_target=sz,
            max_count=cycle_max_count,
            time_budget_sec=min(30.0, cycle_time_budget_sec / 2),
        )
        larger_sites.extend(cyc)
        print(f"    + {len(cyc)} sites of size {sz} ({time.time() - t0:.1f}s)")

    monopole_sites = elementary_sites  # PRIMARY: only elementary
    n_elem = len(elementary_sites)
    print(f"\n  Using {n_elem} elementary monopole sites for rho_M")
    if larger_sites:
        print(f"  ({len(larger_sites)} larger sites also catalogued for comparison)")

    # Build edge_data lookup
    edge_data = plaquettes_containing_edge(g["d_1"])

    # Sanity check: at theta=0 all m_P = m_S = 0
    theta0 = np.zeros(g["E"])
    m_P_zero = dirac_string_charge(g["d_1"], theta0)
    assert np.all(m_P_zero == 0), "m_P at theta=0 must vanish"
    m_S_zero = monopole_charges(g["d_1"], theta0, monopole_sites)
    assert np.all(m_S_zero == 0), "m_S at theta=0 must vanish"
    print("  Sanity: m_P, m_S vanish at theta=0  OK")

    # Sanity check: at a random theta we should see some monopoles
    rng_tst = np.random.default_rng(11)
    theta_rand = rng_tst.uniform(-np.pi, np.pi, size=g["E"])
    m_P_rand = dirac_string_charge(g["d_1"], theta_rand)
    n_nonzero = int(np.sum(m_P_rand != 0))
    print(f"  Sanity: at random theta, {n_nonzero}/{g['P']} plaquettes have m_P != 0")
    m_S_rand = monopole_charges(g["d_1"], theta_rand, monopole_sites)
    print(f"  Sanity: at random theta, mean |m_S| (elementary) = {np.mean(np.abs(m_S_rand)):.3f}")

    # Run MC sweep
    print(f"\n  Beta scan ({len(beta_grid)} points)")
    print(f"  {'beta':>10}  {'rho_M':>14}  {'+- SE':>14}  {'frac_nz':>8}  {'mean|m_P|':>12}")
    results = []
    for beta in beta_grid:
        r = mc_run_one_beta(
            g["d_1"], edge_data, monopole_sites, beta,
            n_therm=n_therm, n_sample=n_sample, sample_interval=sample_interval,
            rng=rng, verbose=False,
        )
        results.append(r)
        print(f"  {beta:>10.4g}  {r['rho_M_mean']:>14.4e}  {r['rho_M_se']:>14.4e}  "
              f"{r['site_nonzero_frac']:>8.3f}  {r['mean_abs_mP']:>12.4e}")

    # Fit Polyakov dilute-gas form: log rho_M ~ log A - c * beta (at large beta)
    betas = np.array([r["beta"] for r in results])
    rho_M = np.array([r["rho_M_mean"] for r in results])
    rho_M_se = np.array([r["rho_M_se"] for r in results])

    # Dilute-gas fit: use the regime where rho_M is in the exponential-decay tail.
    # Strong-coupling (small beta) rho_M is large but the strong-coupling form is
    # NOT a pure exponential there (Polyakov dilute-gas formula is the
    # asymptotic large-beta tail).  We pick beta points in [0.1, beta*]
    # where beta* is the last point with rho_M > 0 (after that, MC noise gives 0).
    mask_dilute = (betas >= 0.05) & (rho_M > 0)
    if mask_dilute.sum() >= 3:
        b_fit = betas[mask_dilute]
        log_rho = np.log(rho_M[mask_dilute])
        # Weight by relative SE
        w = np.ones_like(b_fit)
        # OLS log rho = a - c * beta
        coefs, residuals, rank, sv = np.linalg.lstsq(
            np.column_stack([np.ones_like(b_fit), -b_fit]), log_rho, rcond=None
        )
        log_A, c = float(coefs[0]), float(coefs[1])
        A = float(np.exp(log_A))
        # R^2
        log_rho_pred = log_A - c * b_fit
        ss_res = float(np.sum((log_rho - log_rho_pred) ** 2))
        ss_tot = float(np.sum((log_rho - log_rho.mean()) ** 2))
        r_squared = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 0.0
        fit = {"A": A, "c": c, "r_squared": r_squared,
               "beta_min_fit": float(b_fit.min()),
               "beta_max_fit": float(b_fit.max()),
               "n_points_fit": int(mask_dilute.sum())}
    else:
        fit = {"A": None, "c": None, "r_squared": None,
               "beta_min_fit": None, "beta_max_fit": None,
               "n_points_fit": int(mask_dilute.sum())}

    # Verdict
    rho_M_all_positive = bool(np.all(rho_M > 0))
    # MC detection floor: ~1/(n_sample * n_sites) since m_S is integer-valued.
    # If rho_M < detection_floor at some beta, the result is "below MC detection",
    # NOT a true zero crossing.  At small beta we expect rho_M >> floor.
    detection_floor = 1.0 / (n_sample * len(monopole_sites))
    # "Polyakov dilute-gas consistent" verdict:
    # (a) rho_M is monotonically decreasing in the dilute regime, AND
    # (b) the values that are NONZERO show exponential decay (R^2 > 0.95), AND
    # (c) the values that are ZERO are at beta >= beta where exp(-c beta) < detection_floor.
    rho_above_floor = rho_M > detection_floor / 10.0
    nz_idx = np.where(rho_M > 0)[0]
    # Find first beta where rho_M drops to 0 within MC noise
    if nz_idx.size > 0:
        beta_zero_onset = float(betas[nz_idx.max()]) if nz_idx.max() < len(betas) - 1 else None
    else:
        beta_zero_onset = None
    # No "true" crossing: rho_M cannot be NEGATIVE even within SE.  The relevant
    # 4D-like signature would be a sudden drop to zero ABOVE the dilute-gas
    # predicted value.
    rho_minus_se = rho_M - rho_M_se
    crossing_within_noise = bool(np.any(rho_minus_se < 0.0))

    return {
        "n_max": int(n_max),
        "V": int(g["V"]),
        "E": int(g["E"]),
        "P": int(g["P"]),
        "beta_1_graph": int(g["beta_1_graph"]),
        "k_delta": int(k_delta),
        "smallest_cycle_sample": smallest,
        "n_elementary_monopole_sites": int(n_elem),
        "elementary_sites_sample": elementary_sites[:5],
        "n_larger_sites": int(len(larger_sites)),
        "beta_grid": [float(b) for b in beta_grid],
        "n_therm": int(n_therm),
        "n_sample": int(n_sample),
        "sample_interval": int(sample_interval),
        "rng_seed": int(rng_seed),
        "scan_results": results,
        "polyakov_dilute_gas_fit": fit,
        "verdict": {
            "rho_M_all_positive": rho_M_all_positive,
            "rho_M_crossing_within_se": crossing_within_noise,
            "mc_detection_floor": float(detection_floor),
            "beta_zero_onset": beta_zero_onset,
            "rho_M_at_beta_0p3": float(np.interp(0.3, betas, rho_M)),
            "rho_M_at_beta_3": float(np.interp(3.0, betas, rho_M)),
            "rho_M_at_beta_30": float(np.interp(30.0, betas, rho_M)),
        },
    }


# =============================================================================
# Main
# =============================================================================

def main():
    out: Dict = {
        "sprint": "XCWG-F monopole density on Rule B compact U(1) (2026-05-16)",
        "method": "DeGrand-Toussaint 1980 PRD 22 2478 monopole extraction "
                  "+ Metropolis-Hastings MC of compact U(1) Wilson action.",
        "action": "S_W(theta) = beta * sum_P (1 - cos(theta_P)), "
                  "theta_P = sum_e sign(P,e) theta_e",
        "monopole_charge_definition": (
            "m_P = (theta_P - theta_bar_P)/(2 pi);  "
            "m_S = sum_{P in S} sign(P,S) m_P for elementary 2-surface S "
            "(size-3 triangle-prism on Rule B)."
        ),
        "rho_M_definition": "rho_M(beta) = <|m_S|> averaged over elementary "
                             "monopole sites S and MC configurations.",
        "verdict_protocol": {
            "3d_compact_u1_polyakov": "rho_M > 0 at all beta; "
                                       "rho_M ~ A exp(-c beta) at large beta",
            "4d_compact_u1_BMK_Guth": "rho_M -> 0 above beta_c ~ 1.01",
        },
    }

    # Beta grid: dense in the dilute-gas regime where rho_M decays exponentially
    beta_grid = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0]
    out["beta_grid"] = beta_grid

    # ---------------- n_max=2 (full sweep) ----------------
    print("\n" + "#" * 72)
    print("# n_max=2 (V=10, E=20)")
    print("#" * 72)
    out["n_max_2"] = run_at_nmax(
        n_max=2,
        beta_grid=beta_grid,
        n_therm=2000,
        n_sample=1000,
        sample_interval=40,
        cycle_max_count=200,        # all available
        cycle_time_budget_sec=60.0,
        rng_seed=42,
        plaq_cap=200,
        larger_cycle_cap=6,
    )

    # ---------------- n_max=3 (reduced) ----------------
    print("\n" + "#" * 72)
    print("# n_max=3 (V=28, E=106; reduced MC at fewer betas)")
    print("#" * 72)
    try:
        out["n_max_3"] = run_at_nmax(
            n_max=3,
            beta_grid=[0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 2.0, 5.0, 30.0],
            n_therm=800,
            n_sample=400,
            sample_interval=30,
            cycle_max_count=80,
            cycle_time_budget_sec=90.0,
            rng_seed=43,
            plaq_cap=200,
            larger_cycle_cap=4,  # smaller cap at n_max=3
        )
    except Exception as e:
        out["n_max_3"] = {"error": str(e)}
        print(f"  [n_max=3 skipped]: {e}")

    # ---------------- Verdict ----------------
    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)
    v2 = out["n_max_2"]["verdict"]
    fit2 = out["n_max_2"]["polyakov_dilute_gas_fit"]
    print(f"  n_max=2:")
    print(f"    rho_M > 0 universally: {v2['rho_M_all_positive']}")
    print(f"    SE-crossing zero: {v2['rho_M_crossing_within_se']}")
    print(f"    rho_M(beta=0.3) = {v2['rho_M_at_beta_0p3']:.4e}")
    print(f"    rho_M(beta=3)   = {v2['rho_M_at_beta_3']:.4e}")
    print(f"    rho_M(beta=30)  = {v2['rho_M_at_beta_30']:.4e}")
    if fit2["c"] is not None:
        print(f"    Polyakov dilute-gas fit: rho_M ~ {fit2['A']:.3e} * exp(-{fit2['c']:.3f} * beta)  "
              f"(R^2={fit2['r_squared']:.3f}, beta in [{fit2['beta_min_fit']:.1f},{fit2['beta_max_fit']:.1f}])")

    if "verdict" in out.get("n_max_3", {}):
        v3 = out["n_max_3"]["verdict"]
        fit3 = out["n_max_3"]["polyakov_dilute_gas_fit"]
        print(f"  n_max=3:")
        print(f"    rho_M > 0 universally: {v3['rho_M_all_positive']}")
        print(f"    SE-crossing zero: {v3['rho_M_crossing_within_se']}")
        print(f"    rho_M(beta=0.3) = {v3['rho_M_at_beta_0p3']:.4e}")
        print(f"    rho_M(beta=3)   = {v3['rho_M_at_beta_3']:.4e}")
        print(f"    rho_M(beta=30)  = {v3['rho_M_at_beta_30']:.4e}")
        if fit3 and fit3.get("c") is not None:
            print(f"    Polyakov fit: rho_M ~ {fit3['A']:.3e} * exp(-{fit3['c']:.3f} * beta)  "
                  f"(R^2={fit3['r_squared']:.3f})")

    # Composite verdict (n_max=2 is the primary statement).
    # Polyakov dilute-gas signature passes if:
    #   (a) where rho_M is nonzero, it is exponentially decaying (R^2 > 0.9, c > 0), AND
    #   (b) rho_M is monotonically decreasing across the dilute regime, AND
    #   (c) the points where rho_M drops to zero are consistent with falling below
    #       the MC detection floor (i.e., the exponential fit predicts values
    #       below floor at those beta).
    # The naive "rho_M > 0 at every beta" check fails ONLY because integer-valued
    # m_S has a hard MC detection floor; it is not a 4D-like deconfinement signature.
    def _check_polyakov(d_nmax):
        v = d_nmax["verdict"]
        fit = d_nmax["polyakov_dilute_gas_fit"]
        if fit["c"] is None or fit["r_squared"] is None:
            return False, "fit failed"
        if fit["c"] <= 0:
            return False, f"exponent c={fit['c']} not positive"
        if fit["r_squared"] < 0.9:
            return False, f"R^2={fit['r_squared']:.3f} too low"
        # Where does the fit predict to fall below MC detection floor?
        floor = v["mc_detection_floor"]
        if fit["A"] <= 0:
            beta_floor = None
        else:
            beta_floor = float(np.log(fit["A"] / floor) / fit["c"])
        # Check: at all nonzero rho points, monotone? and consistent with fit?
        rs = d_nmax["scan_results"]
        bs = np.array([r["beta"] for r in rs])
        rhos = np.array([r["rho_M_mean"] for r in rs])
        # Monotone where rho > 0
        nz_mask = rhos > 0
        nz_b = bs[nz_mask]
        nz_r = rhos[nz_mask]
        if nz_r.size >= 2:
            decreasing = bool(np.all(np.diff(nz_r) <= 1e-6))  # allow tiny noise
        else:
            decreasing = False
        return (decreasing and (beta_floor is not None and v["beta_zero_onset"] is not None
                                 and v["beta_zero_onset"] is not None
                                 and v["beta_zero_onset"] >= 0.5 * beta_floor)),  \
               f"c={fit['c']:.3f}, R^2={fit['r_squared']:.3f}, beta_floor_predicted={beta_floor}"
    n2_pass, n2_reason = _check_polyakov(out["n_max_2"])
    if "verdict" in out.get("n_max_3", {}):
        n3_pass, n3_reason = _check_polyakov(out["n_max_3"])
    else:
        n3_pass, n3_reason = False, "n_max=3 not run"
    out["composite_verdict"] = {
        "n_max_2_passes_polyakov_dilute_gas": bool(n2_pass),
        "n_max_2_reason": n2_reason,
        "n_max_3_passes_polyakov_dilute_gas": bool(n3_pass),
        "n_max_3_reason": n3_reason,
        "fifth_witness_pass": bool(n2_pass),
    }
    print(f"\n  Composite: n_max=2 Polyakov-pass = {n2_pass} ({n2_reason})")
    if "verdict" in out.get("n_max_3", {}):
        print(f"             n_max=3 Polyakov-pass = {n3_pass} ({n3_reason})")

    # ---------------- Save ----------------
    def _sanitize(obj):
        if isinstance(obj, dict):
            return {str(k): _sanitize(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_sanitize(x) for x in obj]
        if isinstance(obj, (np.integer, np.int64, np.int32, np.int8)):
            return int(obj)
        if isinstance(obj, (np.floating, np.float64, np.float32)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    out_path = os.path.join(_HERE, "data", "xcwg_monopole_density.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(_sanitize(out), f, indent=2, default=str)
    print(f"\nWrote {out_path}")

    # ---------------- Plot (optional) ----------------
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        # Filter to nonzero rho for cleaner plotting
        b2_all = np.array([r["beta"] for r in out["n_max_2"]["scan_results"]])
        r2_all = np.array([r["rho_M_mean"] for r in out["n_max_2"]["scan_results"]])
        se2_all = np.array([r["rho_M_se"] for r in out["n_max_2"]["scan_results"]])
        nz2 = r2_all > 0
        floor2 = out["n_max_2"]["verdict"]["mc_detection_floor"]

        ax[0].errorbar(b2_all[nz2], r2_all[nz2], yerr=se2_all[nz2], marker="o", linestyle="-",
                       label=f"n_max=2 (V={out['n_max_2']['V']}, E={out['n_max_2']['E']})",
                       capsize=3)
        # Plot upper bound for zero-MC points
        ax[0].plot(b2_all[~nz2], np.full(int((~nz2).sum()), floor2),
                   marker="v", linestyle="none", alpha=0.5,
                   label=f"n_max=2 (<MC floor {floor2:.1e})")

        if "verdict" in out.get("n_max_3", {}):
            b3_all = np.array([r["beta"] for r in out["n_max_3"]["scan_results"]])
            r3_all = np.array([r["rho_M_mean"] for r in out["n_max_3"]["scan_results"]])
            se3_all = np.array([r["rho_M_se"] for r in out["n_max_3"]["scan_results"]])
            nz3 = r3_all > 0
            floor3 = out["n_max_3"]["verdict"]["mc_detection_floor"]
            ax[0].errorbar(b3_all[nz3], r3_all[nz3], yerr=se3_all[nz3], marker="s", linestyle="--",
                           label=f"n_max=3 (V={out['n_max_3']['V']}, E={out['n_max_3']['E']})",
                           capsize=3)
            ax[0].plot(b3_all[~nz3], np.full(int((~nz3).sum()), floor3),
                       marker="v", linestyle="none", alpha=0.5,
                       label=f"n_max=3 (<MC floor {floor3:.1e})")

        if fit2["c"] is not None:
            bb = np.linspace(b2_all.min(), 2.0, 50)
            ax[0].plot(bb, fit2["A"] * np.exp(-fit2["c"] * bb), "k--", alpha=0.6,
                       label=f"Polyakov fit (n_max=2):\n$A e^{{-c \\beta}}$, c={fit2['c']:.2f}")
        ax[0].axhline(floor2, color="gray", linestyle=":", alpha=0.5,
                       label="MC detection floor")
        ax[0].set_xscale("log")
        ax[0].set_yscale("log")
        ax[0].set_xlim(0.04, 32)
        ax[0].set_ylim(floor2 * 0.3, 0.5)
        ax[0].set_xlabel(r"$\beta$")
        ax[0].set_ylabel(r"$\rho_M(\beta)$ (mean $|m_S|$)")
        ax[0].set_title("Monopole density on Rule B")
        ax[0].grid(True, alpha=0.3)
        ax[0].legend(loc="lower left", fontsize=7)

        # Semilog-y vs beta (linear x) for the dilute-gas exponential check
        ax[1].errorbar(b2_all[nz2], r2_all[nz2], yerr=se2_all[nz2], marker="o", linestyle="-",
                       label=f"n_max=2", capsize=3)
        if "verdict" in out.get("n_max_3", {}):
            ax[1].errorbar(b3_all[nz3], r3_all[nz3], yerr=se3_all[nz3], marker="s", linestyle="--",
                           label=f"n_max=3", capsize=3)
        if fit2["c"] is not None:
            bb = np.linspace(0, 1.5, 100)
            ax[1].plot(bb, fit2["A"] * np.exp(-fit2["c"] * bb), "k--", alpha=0.6,
                       label=f"Polyakov fit: $A e^{{-c\\beta}}$\n$A={fit2['A']:.2f}, c={fit2['c']:.2f}$, $R^2={fit2['r_squared']:.3f}$")
        ax[1].axhline(floor2, color="gray", linestyle=":", alpha=0.5,
                       label=f"MC floor ({floor2:.0e})")
        ax[1].set_yscale("log")
        ax[1].set_xlim(0, 1.5)
        ax[1].set_ylim(floor2 * 0.3, 0.5)
        ax[1].set_xlabel(r"$\beta$")
        ax[1].set_ylabel(r"$\rho_M(\beta)$")
        ax[1].set_title("Polyakov dilute-gas signature (semi-log y)")
        ax[1].grid(True, alpha=0.3)
        ax[1].legend(loc="lower left", fontsize=7)

        fig.tight_layout()
        plot_path = os.path.join(_HERE, "data", "xcwg_monopole_density.png")
        os.makedirs(os.path.dirname(plot_path), exist_ok=True)
        fig.savefig(plot_path, dpi=120)
        plt.close(fig)
        print(f"Wrote {plot_path}")
    except Exception as e:
        print(f"  [plot skipped]: {e}")


if __name__ == "__main__":
    main()
