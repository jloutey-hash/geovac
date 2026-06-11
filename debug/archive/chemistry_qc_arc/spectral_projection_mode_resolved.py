#!/usr/bin/env python
"""
Spectral projection mode-resolved analysis: decompose graph-native and
continuum self-energy by internal Fock shell.
==========================================================================

For each n_max in {2, 3, 4}:
  (A) Graph-native: Decompose Sigma by which shell the internal electron
      belongs to.  Sigma_nint[a,b] = sum_{e1,e2} G_gamma[e1,e2]
        * sum_{c in shell(n_int)} V[a,c,e1] * V[c,b,e2]

  (B) Continuum: Decompose the spectral mode sum
        Sigma_cont(n_ext) = sum_{n_int,q} W(...) * g(n_int) * d_T(q)
                            / (|lam(n_int)|^4 * mu(q))
      by restricting n_int to a single shell.

  (C) Extract per-mode ratios rho(n_ext, n_int) = continuum / graph
      and check for patterns (n_ext independence, power law in |lambda|,
      relation to D(s)).

  (D) Save results to JSON and print summary tables.

Uses numpy for n_max >= 4 (sympy eigenvalue bugs).
Uses t = -1/16 (kappa) as the hopping parameter.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import mpmath

# ---------------------------------------------------------------------------
# Add project root to path
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import sympy as sp
from sympy import Rational

from geovac.graph_qed_vertex import (
    build_projection_matrix,
    build_vertex_tensor,
    vertex_tensor_to_matrices,
)
from geovac.graph_qed_photon import (
    build_fock_graph,
    compute_photon_propagator,
)
from geovac.graph_qed_propagator import (
    DiracGraphOperator,
    electron_propagator,
)
from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels


# ---------------------------------------------------------------------------
# Continuum helpers (CH convention: n >= 0)
# ---------------------------------------------------------------------------

def lambda_n_ch(n: int) -> float:
    """Absolute Dirac eigenvalue |lambda_n| = n + 3/2 (CH convention)."""
    return n + 1.5


def g_n_dirac(n: int) -> int:
    """Full Dirac degeneracy g_n = 2(n+1)(n+2)."""
    return 2 * (n + 1) * (n + 2)


def mu_q(q: int) -> int:
    """Hodge-1 Laplacian eigenvalue mu_q = q(q+2), q >= 1."""
    return q * (q + 2)


def d_q_transverse(q: int) -> int:
    """Transverse photon degeneracy d_q^T = q(q+2)."""
    return q * (q + 2)


def vertex_allowed(n1: int, n2: int, q: int) -> bool:
    """SO(4) vertex selection rule: triangle + parity + q >= 1."""
    if q < 1:
        return False
    if q < abs(n1 - n2) or q > n1 + n2:
        return False
    if (n1 + n2 + q) % 2 == 0:
        return False
    return True


def so4_channel_count(n1: int, n2: int, q: int) -> int:
    """Count SO(4) vector harmonic channels (0, 1, or 2)."""
    if not vertex_allowed(n1, n2, q):
        return 0
    j1_L = Fraction(n1 + 1, 2)
    j1_R = Fraction(n1, 2)
    j2_L = Fraction(n2, 2)
    j2_R = Fraction(n2 + 1, 2)
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
# Part A: Graph self-energy decomposed by internal shell
# ---------------------------------------------------------------------------

def graph_self_energy_by_shell(
    n_max: int,
    t: float = -1.0 / 16,
    use_numpy: bool = False,
) -> Dict[str, Any]:
    """Decompose graph-native self-energy by internal electron shell.

    For each n_int shell:
        Sigma_nint[a,b] = sum_{e1,e2} G_gamma[e1,e2]
                          * sum_{c in shell(n_int)} V[a,c,e1] * V[c,b,e2]

    Returns per-shell traces and per-(n_ext, n_int) diagonal blocks.
    """
    print(f"  [graph] Building infrastructure for n_max={n_max} "
          f"({'numpy' if use_numpy else 'sympy'})...")
    t0 = time.time()

    # Build vertex tensor
    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    fock_data = build_fock_graph(n_max)
    E_fock = fock_data.E
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats_sp = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    # Build photon propagator
    # For n_max >= 4, sympy eigenvals() crashes on L1 (complex eigenvalue bug).
    # Build G_gamma via numpy pseudoinverse directly in that case.
    if use_numpy:
        # Build L1 from the Fock graph incidence matrix, then pinv
        L1_np = np.array(fock_data.L1.tolist(), dtype=float)
        G_gamma_np = np.linalg.pinv(L1_np)
    else:
        photon_data = compute_photon_propagator(n_max, exact=True)

    # Convert to numpy if needed
    if use_numpy:
        V_mats_np = [np.array(M.tolist(), dtype=float) for M in V_mats_sp]
        G_gamma = G_gamma_np
    else:
        V_mats = V_mats_sp
        G_gamma = photon_data.G_gamma

    # Identify which Dirac indices belong to each n_fock shell
    shell_indices: Dict[int, List[int]] = {}
    for idx, dl in enumerate(dirac_labels):
        nf = dl.n_fock
        if nf not in shell_indices:
            shell_indices[nf] = []
        shell_indices[nf].append(idx)

    n_fock_shells = sorted(shell_indices.keys())
    print(f"  [graph] Shells: {n_fock_shells}, N_dirac={N_dirac}, "
          f"E_fock={E_fock}, build took {time.time()-t0:.1f}s")

    # Compute shell-decomposed self-energy
    # Sigma_nint[a,b] = sum_{e1,e2} G_gamma[e1,e2]
    #                   * sum_{c in shell(n_int)} V_e1[a,c] * V_e2[c,b]
    # Note: V_e2^T[c,b] = V_e2[b,c], so V_e1 * V_e2^T with restricted c
    #   means: (V_e1[:, shell_c] @ V_e2[:, shell_c].T)

    shell_sigmas = {}  # n_fock -> N_dirac x N_dirac array
    t1 = time.time()

    for nf in n_fock_shells:
        c_indices = shell_indices[nf]

        if use_numpy:
            Sigma_shell = np.zeros((N_dirac, N_dirac))
            for e1 in range(E_fock):
                for e2 in range(E_fock):
                    g_ee = G_gamma[e1, e2]
                    if abs(g_ee) < 1e-15:
                        continue
                    # V_e1[:, c_indices] @ V_e2[:, c_indices].T
                    Ve1_c = V_mats_np[e1][:, c_indices]
                    Ve2_c = V_mats_np[e2][:, c_indices]
                    Sigma_shell += g_ee * (Ve1_c @ Ve2_c.T)
        else:
            Sigma_shell = sp.zeros(N_dirac, N_dirac)
            for e1 in range(E_fock):
                for e2 in range(E_fock):
                    g_ee = G_gamma[e1, e2]
                    if g_ee == 0:
                        continue
                    # Restrict intermediate index to shell
                    for a in range(N_dirac):
                        for b in range(N_dirac):
                            val = sp.S(0)
                            for c in c_indices:
                                val += V_mats[e1][a, c] * V_mats[e2][c, b]
                            if val != 0:
                                Sigma_shell[a, b] += g_ee * val

        shell_sigmas[nf] = Sigma_shell

    print(f"  [graph] Shell decomposition took {time.time()-t1:.1f}s")

    # Extract per-shell traces and per-(n_ext_shell, n_int_shell) blocks
    results = {
        "n_max": n_max,
        "N_dirac": N_dirac,
        "E_fock": E_fock,
        "shells": {},
        "shell_traces": {},
        "gs_block": {},
    }

    # Also compute the per-external-shell diagonal sums
    # For each (n_ext_shell, n_int_shell), sum the diagonal entries
    # of Sigma_nint over external states in n_ext_shell
    n_ext_shell_diag = {}

    for nf_int in n_fock_shells:
        S = shell_sigmas[nf_int]
        if use_numpy:
            tr = float(np.trace(S))
            diag = np.diag(S)
        else:
            tr = float(sp.nsimplify(S.trace(), rational=False))
            diag = [float(S[i, i]) for i in range(N_dirac)]

        results["shell_traces"][str(nf_int)] = tr
        results["shells"][str(nf_int)] = {
            "trace": tr,
            "diag": [float(d) for d in diag],
        }

        # Per external shell diagonal sums
        for nf_ext in n_fock_shells:
            ext_indices = shell_indices[nf_ext]
            diag_sum = sum(float(diag[i]) for i in ext_indices)
            key = f"({nf_ext},{nf_int})"
            n_ext_shell_diag[key] = diag_sum

    results["ext_int_diag_sums"] = n_ext_shell_diag

    # Verify: total trace should match sum of shell traces
    total_trace = sum(results["shell_traces"].values())
    results["total_trace"] = total_trace

    # Ground state (n_fock=1) block of each shell's sigma
    gs_indices = shell_indices.get(1, [])
    for nf_int in n_fock_shells:
        S = shell_sigmas[nf_int]
        if use_numpy:
            gs_block = S[np.ix_(gs_indices, gs_indices)]
            gs_tr = float(np.trace(gs_block))
        else:
            gs_block_entries = []
            gs_tr = 0.0
            for i in gs_indices:
                gs_tr += float(S[i, i])
        results["gs_block"][str(nf_int)] = {
            "trace": gs_tr,
            "n_gs_states": len(gs_indices),
        }

    return results


# ---------------------------------------------------------------------------
# Part B: Continuum self-energy decomposed by internal shell
# ---------------------------------------------------------------------------

def continuum_self_energy_by_shell(
    n_max: int,
) -> Dict[str, Any]:
    """Decompose the continuum spectral self-energy by internal electron shell.

    For each n_ext (CH convention), compute the partial sum over n_int:
        Sigma_cont(n_ext, n_int) = sum_q W(n_ext,n_int,q) * g(n_int) * d_T(q)
                                   / (|lam(n_int)|^4 * mu(q))
    """
    print(f"  [cont] Computing continuum self-energy for n_max={n_max}...")
    t0 = time.time()

    # n_ch ranges: 0 to n_max-1 (since n_fock goes 1..n_max, n_ch = n_fock-1)
    n_ch_max = n_max - 1

    # For each (n_ext_ch, n_int_ch), compute the partial self-energy
    partial = {}  # (n_ext_ch, n_int_ch) -> float

    for n_ext_ch in range(n_ch_max + 1):
        for n_int_ch in range(n_ch_max + 1):
            g_int = g_n_dirac(n_int_ch)
            lam_int = lambda_n_ch(n_int_ch)
            lam_int_pow = lam_int ** 4

            val = 0.0
            q_lo = abs(n_ext_ch - n_int_ch)
            q_hi = n_ext_ch + n_int_ch

            for q in range(max(1, q_lo), q_hi + 1):
                if not vertex_allowed(n_ext_ch, n_int_ch, q):
                    continue
                W = so4_channel_count(n_ext_ch, n_int_ch, q)
                if W == 0:
                    continue
                d_T = d_q_transverse(q)
                mu = mu_q(q)
                val += W * g_int * d_T / (lam_int_pow * mu)

            partial[(n_ext_ch, n_int_ch)] = val

    print(f"  [cont] Computed in {time.time()-t0:.1f}s")

    results = {
        "n_max": n_max,
        "n_ch_max": n_ch_max,
        "partial": {},
    }

    # Store by n_ext and n_int
    for n_ext_ch in range(n_ch_max + 1):
        total_ext = 0.0
        for n_int_ch in range(n_ch_max + 1):
            val = partial[(n_ext_ch, n_int_ch)]
            key = f"({n_ext_ch},{n_int_ch})"
            results["partial"][key] = val
            total_ext += val
        results[f"total_n_ext={n_ext_ch}"] = total_ext

    # Per-n_int totals (sum over n_ext, weighted by g(n_ext))
    for n_int_ch in range(n_ch_max + 1):
        total_int = 0.0
        for n_ext_ch in range(n_ch_max + 1):
            total_int += partial[(n_ext_ch, n_int_ch)]
        results[f"sum_over_next_n_int={n_int_ch}"] = total_int

    return results


# ---------------------------------------------------------------------------
# Part C: Ratio extraction
# ---------------------------------------------------------------------------

def extract_ratios(
    graph_results: Dict,
    cont_results: Dict,
    n_max: int,
) -> Dict[str, Any]:
    """Extract per-mode ratios rho(n_ext, n_int) = continuum / graph.

    The graph self-energy is indexed by n_fock shells (1-indexed);
    the continuum by n_ch shells (0-indexed). n_ch = n_fock - 1.

    For each (n_ext, n_int) pair, compute:
       graph_contrib = sum of Sigma_nint diagonal entries over n_ext shell states
                       / number of states in n_ext shell (average per-state)
       continuum_contrib = Sigma_cont(n_ext_ch, n_int_ch)
       rho = continuum_contrib / graph_contrib  (if graph_contrib != 0)
    """
    n_ch_max = n_max - 1
    ratios = {}
    ratio_table = {}

    for n_ext_ch in range(n_ch_max + 1):
        n_ext_fock = n_ext_ch + 1
        for n_int_ch in range(n_ch_max + 1):
            n_int_fock = n_int_ch + 1

            # Graph: diagonal sum of Sigma_nint over external states
            key_graph = f"({n_ext_fock},{n_int_fock})"
            graph_val = graph_results["ext_int_diag_sums"].get(key_graph, 0.0)

            # Continuum partial
            key_cont = f"({n_ext_ch},{n_int_ch})"
            cont_val = cont_results["partial"].get(key_cont, 0.0)

            # Per-state average for graph (divide by g(n_ext))
            # Actually, the continuum Sigma(n_ext) is already per-state
            # (it's for a single external electron at level n_ext).
            # The graph trace over external shell is sum_{a in ext} Sigma[a,a],
            # which sums over all g(n_ext) Dirac states.
            # To match conventions: Sigma_cont is per-state, graph is total.
            # Ratio: rho = cont_val / (graph_val / g_n_ext)
            #            = cont_val * g_n_ext / graph_val

            # Actually, let's just report both raw numbers and the ratio
            # of per-state averages.
            g_ext = g_n_dirac(n_ext_ch)

            if abs(graph_val) > 1e-15:
                graph_per_state = graph_val / g_ext
                rho = cont_val / graph_per_state
            else:
                graph_per_state = 0.0
                rho = None

            label = f"n_ext_ch={n_ext_ch}, n_int_ch={n_int_ch}"
            ratio_table[label] = {
                "n_ext_ch": n_ext_ch,
                "n_int_ch": n_int_ch,
                "graph_diag_sum": graph_val,
                "graph_per_state": graph_per_state,
                "continuum": cont_val,
                "rho": rho,
                "lam_int": lambda_n_ch(n_int_ch),
                "g_int": g_n_dirac(n_int_ch),
            }

    # Check for n_ext independence of rho
    # For each n_int, collect rho across n_ext
    rho_by_nint = {}
    for label, data in ratio_table.items():
        n_int_ch = data["n_int_ch"]
        rho = data["rho"]
        if rho is not None:
            if n_int_ch not in rho_by_nint:
                rho_by_nint[n_int_ch] = []
            rho_by_nint[n_int_ch].append(rho)

    # Coefficient of variation for each n_int
    rho_cv = {}
    for n_int_ch, rhos in rho_by_nint.items():
        if len(rhos) > 1:
            arr = np.array(rhos)
            mean = np.mean(arr)
            std = np.std(arr)
            cv = std / abs(mean) if abs(mean) > 1e-15 else float('inf')
            rho_cv[str(n_int_ch)] = {
                "mean": float(mean),
                "std": float(std),
                "cv": float(cv),
                "values": [float(r) for r in rhos],
            }
        elif len(rhos) == 1:
            rho_cv[str(n_int_ch)] = {
                "mean": float(rhos[0]),
                "std": 0.0,
                "cv": 0.0,
                "values": [float(rhos[0])],
            }

    # Check for power-law in |lambda_int|
    # Collect (lambda_int, mean_rho) pairs
    lambda_rho_pairs = []
    for n_int_ch_str, data in rho_cv.items():
        n_int_ch = int(n_int_ch_str)
        lam = lambda_n_ch(n_int_ch)
        lambda_rho_pairs.append((lam, data["mean"]))

    power_law_fit = None
    if len(lambda_rho_pairs) >= 2:
        lams = np.array([p[0] for p in lambda_rho_pairs])
        rhos_mean = np.array([p[1] for p in lambda_rho_pairs])
        # Only fit if all rhos are positive (or all negative)
        if np.all(rhos_mean > 0):
            log_lam = np.log(lams)
            log_rho = np.log(rhos_mean)
            if len(log_lam) >= 2:
                coeffs = np.polyfit(log_lam, log_rho, 1)
                power_law_fit = {
                    "exponent": float(coeffs[0]),
                    "prefactor": float(np.exp(coeffs[1])),
                    "log_residual": float(np.std(log_rho - np.polyval(coeffs, log_lam))),
                }

    return {
        "ratio_table": ratio_table,
        "rho_cv_by_nint": rho_cv,
        "lambda_rho_pairs": [(float(l), float(r)) for l, r in lambda_rho_pairs],
        "power_law_fit": power_law_fit,
    }


# ---------------------------------------------------------------------------
# D(s) reference values for comparison
# ---------------------------------------------------------------------------

def dirac_spectral_zeta(s: float, n_max_sum: int = 100) -> float:
    """D(s) = sum_{n=0}^{n_max_sum} g_n * |lambda_n|^{-s}."""
    total = 0.0
    for n in range(n_max_sum + 1):
        total += g_n_dirac(n) * lambda_n_ch(n) ** (-s)
    return total


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def run_analysis() -> Dict[str, Any]:
    """Run the full mode-resolved spectral projection analysis."""

    all_results = {}

    for n_max in [2, 3, 4]:
        print(f"\n{'='*70}")
        print(f"  n_max = {n_max}")
        print(f"{'='*70}")

        use_numpy = (n_max >= 4)

        # Part A: Graph decomposition
        graph_res = graph_self_energy_by_shell(
            n_max, t=-1.0/16, use_numpy=use_numpy
        )

        # Part B: Continuum decomposition
        cont_res = continuum_self_energy_by_shell(n_max)

        # Part C: Ratios
        ratios = extract_ratios(graph_res, cont_res, n_max)

        all_results[f"n_max={n_max}"] = {
            "graph": graph_res,
            "continuum": cont_res,
            "ratios": ratios,
        }

        # Print summary table
        print(f"\n  --- Self-energy shell decomposition (n_max={n_max}) ---")
        print(f"  Graph total trace = {graph_res['total_trace']:.6f}")

        # Shell traces
        print(f"\n  Shell traces (graph):")
        for nf, tr in sorted(graph_res["shell_traces"].items(), key=lambda x: int(x[0])):
            frac = tr / graph_res["total_trace"] if graph_res["total_trace"] != 0 else 0
            print(f"    n_fock={nf}: trace = {tr:.6f}  ({frac*100:.1f}%)")

        # GS block per shell
        print(f"\n  GS (n_fock=1) block trace per internal shell:")
        for nf, data in sorted(graph_res["gs_block"].items(), key=lambda x: int(x[0])):
            print(f"    n_int_fock={nf}: GS trace = {data['trace']:.6f}")

        # Continuum partial sums
        n_ch_max = n_max - 1
        print(f"\n  Continuum self-energy (per n_ext, n_int):")
        for n_ext_ch in range(n_ch_max + 1):
            parts = []
            for n_int_ch in range(n_ch_max + 1):
                key = f"({n_ext_ch},{n_int_ch})"
                val = cont_res["partial"].get(key, 0.0)
                parts.append(f"n_int={n_int_ch}: {val:.6f}")
            total = cont_res.get(f"total_n_ext={n_ext_ch}", 0.0)
            print(f"    n_ext={n_ext_ch}: {', '.join(parts)}  | total={total:.6f}")

        # Ratios
        print(f"\n  Per-mode ratios rho = continuum_per_mode / graph_per_state:")
        header = f"    {'n_ext':>5} {'n_int':>5} {'graph_tot':>12} {'graph/state':>12} {'continuum':>12} {'rho':>12}"
        print(header)
        print(f"    {'-'*len(header)}")
        for label, data in sorted(ratios["ratio_table"].items(),
                                   key=lambda x: (x[1]["n_ext_ch"], x[1]["n_int_ch"])):
            rho_str = f"{data['rho']:.6f}" if data['rho'] is not None else "   N/A"
            print(f"    {data['n_ext_ch']:>5} {data['n_int_ch']:>5} "
                  f"{data['graph_diag_sum']:>12.6f} {data['graph_per_state']:>12.6f} "
                  f"{data['continuum']:>12.6f} {rho_str:>12}")

        # Rho variation across n_ext for each n_int
        print(f"\n  Rho variation across n_ext (per n_int):")
        for n_int_str, data in sorted(ratios["rho_cv_by_nint"].items(),
                                       key=lambda x: int(x[0])):
            print(f"    n_int_ch={n_int_str}: mean={data['mean']:.6f}, "
                  f"CV={data['cv']:.4f}, values={[f'{v:.6f}' for v in data['values']]}")

        # Power law fit
        if ratios["power_law_fit"]:
            pf = ratios["power_law_fit"]
            print(f"\n  Power-law fit rho ~ A * |lambda|^alpha:")
            print(f"    alpha = {pf['exponent']:.4f}, A = {pf['prefactor']:.6f}, "
                  f"log_residual = {pf['log_residual']:.6f}")

    # D(s) reference values for comparison
    print(f"\n{'='*70}")
    print(f"  D(s) reference values (sum to n=100)")
    print(f"{'='*70}")
    for s in [2, 3, 4, 5, 6]:
        ds = dirac_spectral_zeta(s)
        print(f"  D({s}) = {ds:.10f}")

    # Pendant-edge check
    print(f"\n{'='*70}")
    print(f"  Pendant-edge theorem check: Sigma(GS) = 2(n_max-1)/n_max")
    print(f"{'='*70}")
    for nmax_key, data in all_results.items():
        n_max_val = data["graph"]["n_max"]
        total_gs = sum(
            v["trace"]
            for k, v in data["graph"]["gs_block"].items()
        )
        n_gs = list(data["graph"]["gs_block"].values())[0]["n_gs_states"]
        gs_per_state = total_gs / n_gs if n_gs > 0 else 0
        pendant = 2 * (n_max_val - 1) / n_max_val
        print(f"  n_max={n_max_val}: GS trace={total_gs:.6f}, "
              f"per_state={gs_per_state:.6f}, "
              f"pendant={pendant:.6f}, "
              f"match={'YES' if abs(gs_per_state - pendant) < 0.01 else 'NO'}")

    return all_results


# ---------------------------------------------------------------------------
# Save and run
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print("Spectral Projection Mode-Resolved Analysis")
    print("=" * 70)

    results = run_analysis()

    # Prepare JSON-serializable output
    def make_serializable(obj):
        if isinstance(obj, dict):
            return {k: make_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [make_serializable(v) for v in obj]
        elif isinstance(obj, (np.floating, np.integer)):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, float) and (np.isnan(obj) or np.isinf(obj)):
            return str(obj)
        return obj

    output = make_serializable(results)

    out_path = PROJECT_ROOT / "debug" / "data" / "spectral_projection_mode_resolved.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nResults saved to {out_path}")
