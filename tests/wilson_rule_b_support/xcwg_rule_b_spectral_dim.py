"""
Spectral dimension of the Dirac Rule B graph vs the scalar Fock graph.
=======================================================================

Computes the effective spectral dimension d_s via the heat-kernel return
probability for:

  - Dirac Rule B graph (Paper 29 §RH-C): E1 dipole edges, Δl=±1, parity flip
  - Scalar Fock S³ Coulomb graph (Papers 7, 25): same-l ladders only

at n_max = 3, 4, 5.  Cross-checks: Weyl law and graph-diffusion dimension.

The structural question (May 2026 RG sprint context):
  The scalar Fock graph is structurally a disjoint union of 2D rectangular
  grids P_{n_max-l} × P_{2l+1}, hence d_s = 2 at small t.  Does Rule B's
  cross-l hopping restore the expected S³ 3-dimensionality (d_s ≈ 3), or
  does the graph still effectively look 2D?

Output:
  - tests/wilson_rule_b_support/data/xcwg_rule_b_spectral_dim.json
  - tests/wilson_rule_b_support/data/xcwg_rule_b_spectral_dim.png

CONTEXT: GeoVac, Project_Geometric, 2026-05-15.  Migrated from debug/ to
tests/wilson_rule_b_support/ (2026-07-03 durability migration); consumed by
tests/test_paper41_wilson_witnesses.py.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg
import scipy.sparse.csgraph

# Make geovac importable when run standalone from tests/wilson_rule_b_support/.
_REPO = str(Path(__file__).resolve().parents[2])
if _REPO not in sys.path:
    sys.path.insert(0, _REPO)

from geovac.fock_graph_hodge import FockGraphHodge
from geovac.ihara_zeta_dirac import build_dirac_s3_graph


# ---------------------------------------------------------------------------
# Graph builders
# ---------------------------------------------------------------------------

def build_rule_b_laplacian(n_max: int) -> Tuple[np.ndarray, np.ndarray, int]:
    """Return (Laplacian L=D-A, adjacency A, V) for the Dirac Rule B graph."""
    A, _, deg, _ = build_dirac_s3_graph(n_max, "B")
    V = A.shape[0]
    D = np.diag(deg)
    L = (D - A).astype(np.float64)
    return L, A.astype(np.float64), V


def build_scalar_laplacian(n_max: int) -> Tuple[np.ndarray, np.ndarray, int]:
    """Return (Laplacian L, adjacency A, V) for the scalar Fock S³ graph."""
    h = FockGraphHodge(n_max)
    L = h.node_laplacian_numpy.astype(np.float64)
    # adjacency = D - L
    deg = np.diag(L)
    A = np.diag(deg) - L
    V = h.n_nodes
    return L, A, V


# ---------------------------------------------------------------------------
# Heat-kernel return probability  -->  d_s
# ---------------------------------------------------------------------------

def heat_kernel_return_probability(
    eigvals: np.ndarray, V: int, t_values: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (K(t), P(t)) for an array of times t.

    K(t) = sum_n exp(-lambda_n t)            (heat-kernel trace)
    P(t) = K(t) / V                          (return probability per node)

    Numerically stable: subtract min(lambda) before exp.
    """
    # eigvals shape (V,), t_values shape (T,) -> broadcast to (V, T)
    lam = np.clip(eigvals, 0.0, None)  # protect against tiny negative numerical noise
    # K(t) = sum_n exp(-lam_n t)
    # Use that the constant mode lam = 0 contributes c (number of components)
    # exp(-lam * t) with broadcasting
    arg = -np.outer(lam, t_values)  # (V, T)
    K = np.sum(np.exp(arg), axis=0)  # (T,)
    P = K / V
    return K, P


def extract_spectral_dimension(
    t_values: np.ndarray,
    P_values: np.ndarray,
    c: int,
    V: int,
    t_window: Tuple[float, float] = None,
) -> Dict:
    """Extract d_s from log-log slope of P(t) - P(infty) vs t.

    P(infty) = c / V (number of connected components / number of nodes)
    d_s = -2 d log(P(t) - P(infty)) / d log t

    Parameters
    ----------
    t_values, P_values : np.ndarray
        Time grid and return probability values.
    c : int
        Number of connected components (constant-mode contribution).
    V : int
        Number of vertices.
    t_window : (t_min, t_max) or None
        Range of t for the fit.  If None, uses the regime where
        P - P_infty is in [10*P_infty, 0.5 * P(t_min)].

    Returns
    -------
    dict with keys: d_s_global, d_s_window, t_min, t_max, n_points,
                    fit_residual, log_t, log_dP, P_inf
    """
    P_inf = c / V
    dP = P_values - P_inf
    valid = dP > 1e-14
    if not valid.any():
        return {
            "d_s_global": np.nan,
            "d_s_window": np.nan,
            "n_points": 0,
            "fit_residual": np.nan,
            "P_inf": P_inf,
        }
    log_t_all = np.log(t_values[valid])
    log_dP_all = np.log(dP[valid])

    # Auto-window: pick the regime where log_dP is roughly linear in log_t.
    # Avoid: (a) saturation at small t (P close to 1, dP close to 1 - P_inf
    # which is approximately constant when t is so small that K~V); (b)
    # saturation at large t (dP -> 0 and noise dominates).
    if t_window is None:
        # Use the middle 60% by index of the valid range, but restrict to
        # the regime where dP > 5*P_inf (well above asymptote) AND dP < 0.5
        # (well below transient saturation 1 - P_inf).
        usable = (dP[valid] > 5 * P_inf) & (dP[valid] < 0.5)
        if usable.sum() >= 6:
            t_used = t_values[valid][usable]
            log_t_fit = log_t_all[usable]
            log_dP_fit = log_dP_all[usable]
        else:
            # Fallback: use middle 50% of available points
            n_v = log_t_all.size
            lo = max(1, n_v // 4)
            hi = max(lo + 4, 3 * n_v // 4)
            log_t_fit = log_t_all[lo:hi]
            log_dP_fit = log_dP_all[lo:hi]
            t_used = t_values[valid][lo:hi]
    else:
        t_min, t_max = t_window
        mask_w = (t_values[valid] >= t_min) & (t_values[valid] <= t_max)
        log_t_fit = log_t_all[mask_w]
        log_dP_fit = log_dP_all[mask_w]
        t_used = t_values[valid][mask_w]

    # Linear fit: log(dP) = slope * log(t) + intercept => d_s = -2 * slope
    if log_t_fit.size < 4:
        slope_window = np.nan
        residual_window = np.nan
        t_min_w = np.nan
        t_max_w = np.nan
    else:
        coeffs, residuals_arr, _, _, _ = np.polyfit(
            log_t_fit, log_dP_fit, 1, full=True
        )
        slope_window = coeffs[0]
        # residuals_arr is sum-of-squared-residuals; convert to per-point RMS
        if residuals_arr.size:
            residual_window = float(np.sqrt(residuals_arr[0] / log_t_fit.size))
        else:
            residual_window = 0.0
        t_min_w = float(t_used.min())
        t_max_w = float(t_used.max())
    d_s_window = -2.0 * slope_window if not np.isnan(slope_window) else np.nan

    # Global fit (entire valid range) as comparison
    if log_t_all.size >= 4:
        coeffs_g, _, _, _, _ = np.polyfit(log_t_all, log_dP_all, 1, full=True)
        slope_global = coeffs_g[0]
    else:
        slope_global = np.nan
    d_s_global = -2.0 * slope_global if not np.isnan(slope_global) else np.nan

    return {
        "d_s_global": float(d_s_global),
        "d_s_window": float(d_s_window),
        "t_min": t_min_w,
        "t_max": t_max_w,
        "n_points": int(log_t_fit.size) if not np.isnan(slope_window) else 0,
        "fit_residual": residual_window,
        "P_inf": float(P_inf),
        "log_t_fit": log_t_fit.tolist() if not np.isnan(slope_window) else [],
        "log_dP_fit": log_dP_fit.tolist() if not np.isnan(slope_window) else [],
    }


# ---------------------------------------------------------------------------
# Weyl-law dimension
# ---------------------------------------------------------------------------

def weyl_dimension(eigvals: np.ndarray, V: int, c: int) -> Dict:
    """Extract Weyl dimension from cumulative mode count N(lambda).

    For continuous d-dimensional Laplacian on a compact manifold:
        N(lambda) ~ C * lambda^(d/2)
    so log N vs log lambda has slope d/2 -> d = 2*slope.

    Use the intermediate range, avoiding (a) the gap from zero modes (the
    c lowest eigenvalues are exactly 0 within numerical noise) and (b) the
    high-frequency cutoff at the top of the spectrum.
    """
    # Drop the c constant zero modes
    nz = np.sort(eigvals[eigvals > 1e-8])
    if nz.size < 6:
        return {"d_W": np.nan, "n_points": int(nz.size)}

    # N(lambda) at lambda = nz[k] is (k+1) + c (k+1 nonzero + c zero modes,
    # but the slope is insensitive to the +c offset for large N)
    N = np.arange(c + 1, c + 1 + nz.size, dtype=float)
    log_lambda = np.log(nz)
    log_N = np.log(N)

    # Use middle 50% of indices to avoid edge effects
    n_pts = nz.size
    lo = max(1, n_pts // 4)
    hi = max(lo + 4, 3 * n_pts // 4)

    coeffs, residuals_arr, _, _, _ = np.polyfit(
        log_lambda[lo:hi], log_N[lo:hi], 1, full=True
    )
    slope = coeffs[0]
    if residuals_arr.size:
        residual = float(np.sqrt(residuals_arr[0] / (hi - lo)))
    else:
        residual = 0.0
    d_W = 2.0 * slope
    return {
        "d_W": float(d_W),
        "slope": float(slope),
        "fit_residual": residual,
        "n_points": int(hi - lo),
        "lambda_min_fit": float(nz[lo]),
        "lambda_max_fit": float(nz[hi - 1]),
        "n_eigvals_nonzero": int(nz.size),
    }


# ---------------------------------------------------------------------------
# Diffusion dimension <r^2(t)>
# ---------------------------------------------------------------------------

def diffusion_dimension(
    L: np.ndarray, A: np.ndarray, t_values: np.ndarray
) -> Dict:
    """Extract diffusion dimension from <r^2(t)>_node.

    For a continuous-time random walk with generator L, starting at node v0,
    the probability of being at node u at time t is p(u, v0; t) = (exp(-Lt))[u, v0].
    Mean square graph distance from v0:
        <r^2(t)>_{v0} = sum_u r(u, v0)^2 * p(u, v0; t)
    where r is shortest-path graph distance (number of edges).

    Average over starting nodes v0 to get an averaged <r^2(t)>.

    For a d-dimensional continuous space:
        <r^2(t)> ~ t       (d_diff = 1 in log-log slope, but actually 2/d_s = ?)
    Wait, careful: on a graph with spectral dimension d_s,
        <r^2(t)> ~ t^(2/d_w)
    where d_w is the WALK dimension.  For Euclidean-like graphs d_w = 2 so
    <r^2> ~ t.  For fractal-like graphs d_w > 2.
    The relationship is d_s = 2 * d_f / d_w where d_f is fractal/Hausdorff dim.

    For a "nice" Euclidean-like graph: d_s = d_f and d_w = 2 so
        slope[log <r^2> vs log t] = 1
        d_diff (informal) = 2 / slope = 2

    To get a dimensional readout that matches d_s on Euclidean graphs,
    we use the convention:
        d_diff = 2 / (d log <r^2> / d log t)   when d_w = 2 (random walk)

    This recovers d_diff = 2 for ordinary diffusion (slope = 1) regardless
    of spectral dimension.  Therefore d_diff alone does NOT identify d_s
    unless we also know d_w.  We report the slope as the primary observable.
    """
    V = L.shape[0]
    # Shortest-path distances (BFS, since unweighted)
    A_int = (A > 0.5).astype(int)
    dist = scipy.sparse.csgraph.shortest_path(A_int, directed=False, unweighted=True)
    # Set inf (disconnected pairs) to a large finite value or skip per-row
    finite_mask = np.isfinite(dist)
    # Square the distances
    dist_sq = np.where(finite_mask, dist ** 2, 0.0)

    # Eigendecompose L once (V is small)
    evals, evecs = scipy.linalg.eigh(L)
    evals = np.clip(evals, 0.0, None)

    # exp(-L t) = V_evec diag(exp(-evals t)) V_evec^T
    # p(u, v0; t) = sum_k V_evec[u, k] exp(-evals[k] t) V_evec[v0, k]

    # <r^2(t)>_{v0} = sum_u dist_sq[u, v0] p(u, v0; t)
    # Averaged over v0: <<r^2(t)>>_v0
    #   = (1/V) sum_{v0} sum_u dist_sq[u, v0] p(u, v0; t)
    # = (1/V) sum_{u, v0, k} dist_sq[u, v0] V[u,k] exp(-evals[k]*t) V[v0,k]
    # Define M[k, k'] = sum_{u, v0} dist_sq[u, v0] V[u,k] V[v0, k'] then
    # <<r^2(t)>> = (1/V) sum_k exp(-evals[k]*t) M[k, k] (since k=k' due to symmetric V)
    # Wait, this isn't right since dist_sq is not generally diagonalized by V.

    # Simpler: compute <r^2(t)> for each starting node v0, then average.
    r2_avg = np.zeros_like(t_values)
    for v0 in range(V):
        # Project initial delta state onto eigenbasis
        c_k = evecs[v0, :]  # initial amplitudes
        # p_t[u] = sum_k V[u, k] exp(-evals[k] t) c_k
        # r^2(v0; t) = sum_u dist_sq[u, v0] * p_t[u]
        # Vectorize over t: r^2(v0; t) = sum_u dist_sq[u, v0] *
        #    (V_evec @ diag(exp(-evals*t)) @ c_k)
        # = sum_u dist_sq[u, v0] * sum_k V_evec[u, k] * exp(-evals[k]*t) * c_k[k]
        # = sum_k (sum_u dist_sq[u, v0] * V_evec[u, k]) * c_k[k] * exp(-evals[k]*t)
        # Let w[k] = sum_u dist_sq[u, v0] V_evec[u, k] * c_k[k]
        w_k = (dist_sq[:, v0] @ evecs) * c_k  # (V,)
        # r^2(v0; t) = sum_k w_k * exp(-evals[k] * t)
        contrib = w_k[:, None] * np.exp(-np.outer(evals, t_values))  # (V, T)
        r2_v0 = contrib.sum(axis=0)
        r2_avg += r2_v0
    r2_avg /= V

    # Slope: log <r^2> vs log t in the diffusive regime
    valid = r2_avg > 1e-14
    if valid.sum() < 6:
        return {
            "slope": np.nan,
            "d_diff": np.nan,
            "n_points": 0,
        }
    log_t = np.log(t_values[valid])
    log_r2 = np.log(r2_avg[valid])

    # Use intermediate window: where r^2 has grown well beyond 0 but not yet
    # saturated (to the average squared distance over the graph).
    diam = float(np.where(np.isfinite(dist), dist, 0).max())
    mean_d2 = float(np.where(np.isfinite(dist), dist_sq, 0).mean())
    cap = 0.7 * mean_d2  # below 70% of saturation
    floor = 0.05 * mean_d2  # above 5% of saturation

    in_band = (r2_avg[valid] > floor) & (r2_avg[valid] < cap)
    if in_band.sum() >= 6:
        coeffs, residuals_arr, _, _, _ = np.polyfit(
            log_t[in_band], log_r2[in_band], 1, full=True
        )
        slope = coeffs[0]
        n_pts_used = int(in_band.sum())
        if residuals_arr.size:
            res = float(np.sqrt(residuals_arr[0] / n_pts_used))
        else:
            res = 0.0
    else:
        # Fallback: middle 50% of valid range
        n_v = log_t.size
        lo = max(1, n_v // 4)
        hi = max(lo + 4, 3 * n_v // 4)
        coeffs, residuals_arr, _, _, _ = np.polyfit(
            log_t[lo:hi], log_r2[lo:hi], 1, full=True
        )
        slope = coeffs[0]
        n_pts_used = hi - lo
        if residuals_arr.size:
            res = float(np.sqrt(residuals_arr[0] / n_pts_used))
        else:
            res = 0.0

    # On a "nice" Euclidean-like graph: d_w = 2 (random walk dimension),
    # slope = 1, and d_diff_implied = 2 (NOT d_s).  If slope deviates from
    # 1 it indicates anomalous diffusion (d_w != 2).  We report the slope
    # itself and the saturation value.
    return {
        "slope": float(slope),
        "d_diff_assuming_dw_eq_2": float(2.0 / slope) if slope > 0 else np.nan,
        "fit_residual": res,
        "n_points": int(n_pts_used),
        "graph_diameter": diam,
        "mean_dist_squared": mean_d2,
        "r2_saturation": float(r2_avg[valid][-1]),
        "r2_curve": r2_avg.tolist(),
    }


# ---------------------------------------------------------------------------
# Connected components count
# ---------------------------------------------------------------------------

def count_components(A: np.ndarray) -> int:
    A_int = (A > 0.5).astype(int)
    n_comp, _ = scipy.sparse.csgraph.connected_components(A_int, directed=False)
    return int(n_comp)


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------

def analyze_graph(L: np.ndarray, A: np.ndarray, V: int, name: str,
                   t_values: np.ndarray) -> Dict:
    """Run the full d_s extraction for one graph."""
    c = count_components(A)
    eigvals = np.linalg.eigvalsh(L)
    eigvals = np.clip(eigvals, 0.0, None)

    K, P = heat_kernel_return_probability(eigvals, V, t_values)
    ds_result = extract_spectral_dimension(t_values, P, c, V)
    weyl_result = weyl_dimension(eigvals, V, c)
    diff_result = diffusion_dimension(L, A, t_values)

    return {
        "name": name,
        "V": V,
        "E": int(A.sum() / 2),
        "c_components": c,
        "max_degree": int(np.diag(L).max()),
        "lambda_max": float(eigvals.max()),
        "lambda_1_nonzero": float(np.sort(eigvals[eigvals > 1e-8])[0])
            if (eigvals > 1e-8).sum() > 0 else float("nan"),
        "spectral_dimension": ds_result,
        "weyl_dimension": weyl_result,
        "diffusion": diff_result,
        "eigenvalues": eigvals.tolist(),
        "K_t": K.tolist(),
        "P_t": P.tolist(),
        "t_values": t_values.tolist(),
    }


def make_plots(results: Dict, t_values: np.ndarray, out_path: Path) -> None:
    """Generate diagnostic plots of P(t) and heat-kernel slopes."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    rule_b = results["rule_b"]
    scalar = results["scalar"]
    n_max_values = sorted(rule_b.keys())

    # Row 1: log(P(t) - P_inf) vs log t  for both graphs at each n_max
    for col, n_max in enumerate(n_max_values):
        ax = axes[0, col]
        for r, label, color in [
            (rule_b[n_max], "Rule B (Dirac)", "C0"),
            (scalar[n_max], "Scalar Fock", "C1"),
        ]:
            P = np.asarray(r["P_t"])
            P_inf = r["spectral_dimension"]["P_inf"]
            dP = P - P_inf
            valid = dP > 1e-14
            ax.loglog(
                t_values[valid], dP[valid],
                label=f"{label} (d_s={r['spectral_dimension']['d_s_window']:.2f})",
                color=color, alpha=0.8,
            )
        # Reference slopes
        # d_s = 2 -> slope -1 in dP vs t
        # d_s = 3 -> slope -1.5
        t_ref = t_values[t_values >= 0.5]
        if t_ref.size:
            ax.loglog(t_ref, 0.5 * (t_ref / t_ref[0]) ** (-1.0),
                       "k:", alpha=0.5, label="slope -1 (d_s=2)")
            ax.loglog(t_ref, 0.5 * (t_ref / t_ref[0]) ** (-1.5),
                       "k--", alpha=0.5, label="slope -1.5 (d_s=3)")
        ax.set_xlabel("t")
        ax.set_ylabel("P(t) - P(infty)")
        ax.set_title(f"Heat-kernel return prob (n_max={n_max})")
        ax.legend(fontsize=8)
        ax.grid(True, which="both", alpha=0.3)

    # Row 2 col 0: cumulative mode count N(lambda) for n_max=5
    ax = axes[1, 0]
    n_max = max(n_max_values)
    for r, label, color in [
        (rule_b[n_max], "Rule B", "C0"),
        (scalar[n_max], "Scalar", "C1"),
    ]:
        eigvals = np.array(r["eigenvalues"])
        nz = np.sort(eigvals[eigvals > 1e-8])
        N = np.arange(1, nz.size + 1)
        ax.loglog(nz, N, label=f"{label} (d_W={r['weyl_dimension']['d_W']:.2f})",
                   color=color, alpha=0.8)
    # Reference: lambda^(3/2) and lambda^1
    nz = np.sort(np.array(rule_b[n_max]["eigenvalues"]))
    nz = nz[nz > 1e-8]
    if nz.size:
        ax.loglog(nz, 0.5 * (nz / nz[0]) ** 1.0, "k:", alpha=0.5,
                   label="lambda^1 (d_W=2)")
        ax.loglog(nz, 0.5 * (nz / nz[0]) ** 1.5, "k--", alpha=0.5,
                   label="lambda^(3/2) (d_W=3)")
    ax.set_xlabel("lambda")
    ax.set_ylabel("N(lambda) (cumulative mode count)")
    ax.set_title(f"Weyl law (n_max={n_max})")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    # Row 2 col 1: <r^2(t)> diffusion
    ax = axes[1, 1]
    n_max = max(n_max_values)
    for r, label, color in [
        (rule_b[n_max], "Rule B", "C0"),
        (scalar[n_max], "Scalar", "C1"),
    ]:
        r2 = np.array(r["diffusion"]["r2_curve"])
        sat = r["diffusion"]["mean_dist_squared"]
        slope = r["diffusion"]["slope"]
        valid = r2 > 1e-14
        ax.loglog(t_values[valid], r2[valid],
                   label=f"{label} (slope={slope:.2f}, sat={sat:.2f})",
                   color=color, alpha=0.8)
        ax.axhline(sat, color=color, linestyle=":", alpha=0.3)
    # Reference slope 1 (normal diffusion, d_w=2)
    t_ref = t_values[t_values >= 0.5]
    if t_ref.size:
        ax.loglog(t_ref, 0.3 * (t_ref / t_ref[0]),
                   "k:", alpha=0.5, label="slope 1 (d_w=2)")
    ax.set_xlabel("t")
    ax.set_ylabel("<r^2(t)>")
    ax.set_title(f"Graph-distance diffusion (n_max={n_max})")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)

    # Row 2 col 2: d_s convergence with n_max
    ax = axes[1, 2]
    nm_arr = np.array(n_max_values)
    ds_rb = np.array([rule_b[n]["spectral_dimension"]["d_s_window"]
                       for n in n_max_values])
    ds_sc = np.array([scalar[n]["spectral_dimension"]["d_s_window"]
                       for n in n_max_values])
    dw_rb = np.array([rule_b[n]["weyl_dimension"]["d_W"]
                       for n in n_max_values])
    dw_sc = np.array([scalar[n]["weyl_dimension"]["d_W"]
                       for n in n_max_values])
    ax.plot(nm_arr, ds_rb, "o-", label="Rule B  d_s", color="C0")
    ax.plot(nm_arr, ds_sc, "s-", label="Scalar  d_s", color="C1")
    ax.plot(nm_arr, dw_rb, "o--", label="Rule B  d_W", color="C0", alpha=0.5)
    ax.plot(nm_arr, dw_sc, "s--", label="Scalar  d_W", color="C1", alpha=0.5)
    ax.axhline(3, color="k", linestyle=":", alpha=0.5, label="d=3 (S^3 target)")
    ax.axhline(2, color="gray", linestyle=":", alpha=0.5, label="d=2 (P x P grid)")
    ax.set_xlabel("n_max")
    ax.set_ylabel("dimension")
    ax.set_title("d_s and d_W vs n_max")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Spectral dimension: Dirac Rule B vs Scalar Fock S^3 graph",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] saved {out_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    here = Path(__file__).resolve().parent
    data_path = here / "data" / "xcwg_rule_b_spectral_dim.json"
    plot_path = here / "data" / "xcwg_rule_b_spectral_dim.png"
    data_path.parent.mkdir(parents=True, exist_ok=True)
    plot_path.parent.mkdir(parents=True, exist_ok=True)

    # Time grid: span ~6 decades centered on the relevant scale.  The
    # smallest eigenvalue gap sets the long-t scale; the largest
    # eigenvalue sets the short-t scale.  For n_max = 5 Rule B,
    # lambda_max ~ degree_max ~ 25, so t ~ 1/25 = 0.04 is the
    # short-t cutoff and t ~ 1/lambda_1 ~ 1 is the medium-t scale.
    # Use a wide log-spaced grid.
    t_values = np.logspace(-2.5, 2.5, 200)

    n_max_values = [3, 4, 5]
    results = {
        "metadata": {
            "sprint": "XCWG (spectral dimension of Rule B)",
            "date": "2026-05-15",
            "n_max_values": n_max_values,
            "t_min": float(t_values.min()),
            "t_max": float(t_values.max()),
            "n_time_points": int(t_values.size),
            "method_notes": (
                "d_s extracted from log-log slope of (P(t) - P(infty)) vs t "
                "in the intermediate regime: dP between 5*P_inf and 0.5.  "
                "Constant-mode contribution P(infty) = c/V where c = number "
                "of connected components.  d_s = -2 * d log(P - P_inf)/d log t."
            ),
        },
        "rule_b": {},
        "scalar": {},
    }

    for n_max in n_max_values:
        print(f"\n[n_max={n_max}]")
        # Rule B
        L_b, A_b, V_b = build_rule_b_laplacian(n_max)
        r_b = analyze_graph(L_b, A_b, V_b, f"Rule B n_max={n_max}", t_values)
        ds = r_b["spectral_dimension"]["d_s_window"]
        dW = r_b["weyl_dimension"]["d_W"]
        print(
            f"  Rule B:    V={V_b}, E={r_b['E']}, c={r_b['c_components']}, "
            f"d_s={ds:.3f}, d_W={dW:.3f}"
        )
        results["rule_b"][n_max] = r_b

        # Scalar
        L_s, A_s, V_s = build_scalar_laplacian(n_max)
        r_s = analyze_graph(L_s, A_s, V_s, f"Scalar n_max={n_max}", t_values)
        ds = r_s["spectral_dimension"]["d_s_window"]
        dW = r_s["weyl_dimension"]["d_W"]
        print(
            f"  Scalar:    V={V_s}, E={r_s['E']}, c={r_s['c_components']}, "
            f"d_s={ds:.3f}, d_W={dW:.3f}"
        )
        results["scalar"][n_max] = r_s

    # Strip huge eigenvalue/heat-kernel arrays from final JSON write to
    # keep file size manageable, but keep a summary
    # (full t-grid arrays kept since useful for verification)
    with data_path.open("w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n[data] wrote {data_path}")

    make_plots(results, t_values, plot_path)

    # Print final summary table
    print("\n" + "=" * 70)
    print("SUMMARY  (d_s from heat kernel,  d_W from Weyl law,  slope from <r^2(t)>)")
    print("=" * 70)
    print(f"{'graph':<10} {'n_max':<6} {'V':<5} {'E':<5} {'c':<3} "
          f"{'d_s':<8} {'d_W':<8} {'r2_slope':<10}")
    for n_max in n_max_values:
        for key, label in [("rule_b", "Rule B"), ("scalar", "Scalar")]:
            r = results[key][n_max]
            ds = r["spectral_dimension"]["d_s_window"]
            dW = r["weyl_dimension"]["d_W"]
            slope = r["diffusion"]["slope"]
            print(
                f"{label:<10} {n_max:<6} {r['V']:<5} {r['E']:<5} "
                f"{r['c_components']:<3} {ds:<8.3f} {dW:<8.3f} {slope:<10.3f}"
            )


if __name__ == "__main__":
    main()
