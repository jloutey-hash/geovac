"""
Track RH-G: Combinatorial limit toward classical Riemann zeta.

Four investigations on the GeoVac Hopf graphs (scalar S^3 Coulomb,
scalar S^5 Bargmann-Segal, Dirac-S^3 Rules A/B):

  (I)   Primitive-walk length distribution a_L / L vs pi(L)/log(L).
  (II)  Explicit-formula fluctuating term vs RH-zero statistics.
  (III) chi_{-4} twisted Ihara zeta: does it relate to L(s, chi_{-4})?
  (IV)  Continuum limit of Ihara-Bass via N_max -> infinity.

Output: debug/data/riemann_limit_data.json (all numerical data),
then interpreted in debug/riemann_limit_memo.md.
"""

from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp

from geovac.ihara_zeta import (
    _as_integer_adjacency,
    _count_closed_nonbacktracking_walks,
    _count_components,
    _degree_sequence,
    _mobius_to_primitive,
    hashimoto_matrix,
    ihara_zeta_bass,
    is_ramanujan,
)
from geovac.lattice import GeometricLattice
from geovac.nuclear.bargmann_graph import build_bargmann_graph
from geovac.ihara_zeta_dirac import build_dirac_s3_graph

DATA_DIR = Path(__file__).parent / "data"
DATA_DIR.mkdir(exist_ok=True)


# ----------------------------------------------------------------------------
# Graph builders
# ----------------------------------------------------------------------------

def graph_s3(max_n: int):
    """Return (A_connectivity, node_labels) for the scalar S^3 Coulomb graph."""
    lat = GeometricLattice(max_n=max_n)
    A_raw = lat.adjacency
    if hasattr(A_raw, "toarray"):
        A = np.asarray(A_raw.toarray(), dtype=int)
    else:
        A = np.asarray(A_raw, dtype=int)
    # Coerce to 0/1 connectivity
    A = (np.abs(A) > 0).astype(int)
    np.fill_diagonal(A, 0)
    return A, list(lat.states)


def graph_s5(N_max: int):
    """Return (A_connectivity, node_labels) for the Bargmann-Segal S^5 graph."""
    bg = build_bargmann_graph(N_max=N_max)
    A_weighted = np.asarray(bg.adjacency_dense())
    V = A_weighted.shape[0]
    A = (np.abs(A_weighted) > 1e-12).astype(int)
    np.fill_diagonal(A, 0)
    nodes = bg.nodes if hasattr(bg, "nodes") else []
    return A, nodes


def graph_dirac_s3(n_max: int, rule: str):
    """Return (A, node_labels) for the Dirac-S^3 graph, rule 'A' or 'B'."""
    A, nodes, _, _ = build_dirac_s3_graph(n_max=n_max, adjacency_rule=rule)
    A = np.asarray(A, dtype=int)
    return A, nodes


# ----------------------------------------------------------------------------
# Investigation I: Primitive-walk length distribution
# ----------------------------------------------------------------------------

def investigation_I(A: np.ndarray, L_max: int = 20):
    """
    Compute a_L (primitive closed walks of length L) for L = 1..L_max.
    Compare normalization a_L/L to pi(L)/log(L).
    """
    tr = _count_closed_nonbacktracking_walks(A, L_max)
    prim = _mobius_to_primitive(tr)
    # Normalize: a_L / L and compare to pi(L)/log(L)
    data = []
    for L in range(1, L_max + 1):
        a_L = prim.get(L, 0)
        if a_L is None:
            a_L = 0
        # "walk-density" = a_L / L
        density = a_L / L if L > 0 else 0
        # reference: pi(x)/log(x); pi(L) via sympy; for L=1 pi(1)=0
        if L >= 2:
            pi_L = sp.primepi(L)
            ref = float(pi_L) / float(np.log(L))
        else:
            ref = 0
        data.append({
            "L": L,
            "tr_T_L": int(tr[L]),
            "a_L": int(a_L),
            "a_L_over_L": float(density),
            "pi_L_over_log_L": float(ref),
        })
    return data


# ----------------------------------------------------------------------------
# Investigation II: Explicit-formula fluctuating term
# ----------------------------------------------------------------------------

def investigation_II(A: np.ndarray, L_max: int = 20):
    """
    Explicit-formula analog: write down the Ihara-analog of Riemann's
    explicit formula and evaluate the fluctuating term.

    For Ihara zeta, the primitive-walk counting function
        Psi(L) = sum_{L'|L} L' a_{L'}  (= tr(T^L))
    has the "explicit formula"
        Psi(L) = rho(T)^L + conj rho(T)^L + sum_{mu nontriv} mu^L
    Rearranging: the fluctuating term
        Osc(L) := sum_{mu nontriv} mu^L
    should — if anything RH-ish is going on — have moments/statistics
    governed by the distribution of nontrivial T-eigenvalues.

    We compute Osc(L) directly and compare to the GUE pair-correlation
    null (Berry-Tabor Poisson vs RMT GUE) on the normalized argument spec.
    """
    T = hashimoto_matrix(A).astype(float)
    ev = np.linalg.eigvals(T)
    mags = np.abs(ev)
    rho = float(mags.max())

    # Identify trivial eigenvalues: those with |mu| within tol of the Perron value.
    tol = 1e-8
    triv_mask = np.abs(mags - rho) < tol
    nontriv_eigs = ev[~triv_mask]

    # Osc(L) = sum of mu^L over non-trivial eigenvalues.
    # This is real because the eigenvalues come in complex-conjugate pairs
    # (T has real entries).
    data = []
    for L in range(1, L_max + 1):
        Psi_L = float(np.real(np.sum(ev ** L)))
        Osc_L = float(np.real(np.sum(nontriv_eigs ** L)))
        # normalize: Osc_L / sqrt(|V_nontriv|) * rho^{-L}
        N_nt = len(nontriv_eigs)
        scale = (rho ** L) if rho > 0 else 1
        if scale > 0:
            osc_normalized = Osc_L / scale
        else:
            osc_normalized = 0
        data.append({
            "L": L,
            "Psi_L": Psi_L,
            "Osc_L": Osc_L,
            "Osc_normalized": float(osc_normalized),
            "rho_T": rho,
            "n_nontriv": int(N_nt),
        })

    # Basic pair-correlation statistics on the imaginary parts (Im mu),
    # normalized to unit mean spacing. For Ramanujan regular graphs this is
    # classically known to be GUE-like; for irregular graphs it is unknown.
    # We provide raw imag-part spacings for external comparison.
    im_parts = np.sort(np.imag(nontriv_eigs[np.imag(nontriv_eigs) > 1e-9]))
    if len(im_parts) > 1:
        spacings = np.diff(im_parts)
        mean_s = float(spacings.mean())
        normalized_spacings = spacings / mean_s if mean_s > 0 else spacings
        spacing_stats = {
            "n_imag_pos": len(im_parts),
            "mean_spacing": float(mean_s),
            "spacing_cv": float(np.std(spacings) / mean_s) if mean_s > 0 else 0,
            "normalized_spacings": [float(x) for x in normalized_spacings],
        }
    else:
        spacing_stats = {"n_imag_pos": len(im_parts), "note": "too few for statistics"}

    return {
        "series": data,
        "spacing_statistics": spacing_stats,
    }


# ----------------------------------------------------------------------------
# Investigation III: chi_{-4} twisted Ihara zeta  (HEADLINE)
# ----------------------------------------------------------------------------

def chi_minus4(n: int) -> int:
    """Dirichlet character mod 4: chi_{-4}(n) = +1 if n = 1 mod 4,
                                  = -1 if n = 3 mod 4,
                                  =  0 if n even."""
    m = n % 4
    if m == 1:
        return 1
    if m == 3:
        return -1
    return 0


def twisted_ihara_zeta_series(A: np.ndarray, L_max: int = 20):
    """
    The twisted Ihara zeta is defined by the Euler product with a
    character twist on the length:

        zeta_G(s; chi) = prod_{[C] primitive} (1 - chi(L(C)) s^{L(C)})^{-1}

    To verify: via the log-derivative (series in s),
        -s * d/ds log zeta_G(s; chi) = sum_{L>=1} chi(L) * Psi(L) * s^L
                                     = sum_{L>=1} chi(L) * tr(T^L) * s^L
    (Here Psi(L) = tr(T^L) = sum over primitive walks with multiplicity L
     of length d dividing L etc. — same as Riemann's von Mangoldt-style
     Psi.)

    For the Riemann case, log zeta(s, chi) = sum_p (chi(p)/p^s + ...)
    and the corresponding Psi has chi as a weight. The analog here is:
    does the series

        Lambda(s; chi) = sum_{L=1}^{L_max} chi(L) * a_L * s^L

    (with a_L = primitive closed walks of length L) have a meaningful
    analytic structure — specifically, does it mirror the series
    L(s, chi_{-4}) = 1 - 3^{-s} + 5^{-s} - 7^{-s} + ...?

    We compute both and compare.
    """
    tr = _count_closed_nonbacktracking_walks(A, L_max)
    prim = _mobius_to_primitive(tr)

    # Build two versions:
    # (a) twisted trace series: sum_L chi(L) * tr(T^L) * s^L  (series in s)
    # (b) twisted primitive series: sum_L chi(L) * a_L * s^L

    twisted_trace = []
    twisted_primitive = []
    twisted_trace_by_parity = {"odd_only": [], "even_only": []}
    for L in range(1, L_max + 1):
        chi = chi_minus4(L)
        tr_L = int(tr[L])
        a_L = prim.get(L, 0) or 0
        twisted_trace.append({"L": L, "chi_-4(L)": chi, "tr_T_L": tr_L,
                              "chi_tr_L": chi * tr_L})
        twisted_primitive.append({"L": L, "chi_-4(L)": chi, "a_L": int(a_L),
                                  "chi_a_L": chi * int(a_L)})
        # For comparison with L(s, chi_-4) series structure:
        # L(s, chi_-4) has support only on odd integers. The "odd-only"
        # contribution of the graph-side mirrors this.
        if L % 2 == 1:
            twisted_trace_by_parity["odd_only"].append({
                "L": L, "tr_T_L": tr_L, "sign": chi
            })

    # Dirichlet L(s, chi_{-4}) reference coefficients (first 20 odd L):
    # L(s, chi_{-4}) = sum_{n=1}^inf chi_{-4}(n)/n^s
    #               = 1 - 3^{-s} + 5^{-s} - 7^{-s} + 9^{-s} - ...
    L_refs = []
    for n in range(1, L_max + 1):
        c = chi_minus4(n)
        if c != 0:
            L_refs.append({"n": n, "chi(n)": c})

    # The key comparison: on the graph side the "twisted Psi(s)" coefficients
    # are {chi(L) * tr(T^L)}; on the L(s, chi) side they are {chi(L)}.
    # The ratios tr(T^L)/1 = tr(T^L) tell us how the graph weighting deviates
    # from the trivial "all primes have weight 1" Riemann case.
    return {
        "twisted_trace_series": twisted_trace,
        "twisted_primitive_series": twisted_primitive,
        "L_chi_minus4_reference_coefficients": L_refs,
    }


# ----------------------------------------------------------------------------
# Investigation IV: Continuum limit
# ----------------------------------------------------------------------------

def investigation_IV_s5(N_max_list=(2, 3)):
    """
    Continuum-limit exploration: track the normalized Ihara-Bass data
    across a sequence of N_max values. Record
      (V, E, beta_1, rho(T), sqrt(q_max), deviation, #nontriv_zeros)
    and look for patterns in scaling.
    """
    series = []
    for N in N_max_list:
        A, _nodes = graph_s5(N_max=N)
        V = A.shape[0]
        E = int(A.sum()) // 2
        d = _degree_sequence(A)
        q_max = int(d.max()) - 1
        T = hashimoto_matrix(A).astype(float)
        ev = np.linalg.eigvals(T)
        rho = float(np.max(np.abs(ev)))
        # Non-trivial nontriv-max-magnitude
        tol = 1e-8
        mags = np.abs(ev)
        triv_mask = np.abs(mags - rho) < tol
        mu_nt_max = float(mags[~triv_mask].max()) if (~triv_mask).any() else 0
        series.append({
            "N_max": N,
            "V": V, "E": E,
            "q_max": q_max,
            "rho_T": rho,
            "sqrt_q_max": float(np.sqrt(q_max)) if q_max > 0 else 0,
            "mu_nt_max": mu_nt_max,
            "deviation": mu_nt_max - float(np.sqrt(q_max)) if q_max > 0 else 0,
            "n_eigs_total": int(len(ev)),
        })
    return series


# ----------------------------------------------------------------------------
# Master driver
# ----------------------------------------------------------------------------

def run_all():
    output = {"sprint": "RH-G", "description": "Combinatorial limit toward Riemann"}

    # Build all graphs
    print("Building graphs...")
    graphs = {}
    A_s3_3, _ = graph_s3(max_n=3)
    graphs["S3_Coulomb_max_n_3"] = A_s3_3
    A_s5_2, _ = graph_s5(N_max=2)
    graphs["S5_Bargmann_Segal_N_max_2"] = A_s5_2
    A_s5_3, _ = graph_s5(N_max=3)
    graphs["S5_Bargmann_Segal_N_max_3"] = A_s5_3

    # Dirac-S^3 too
    A_dirac_A, _ = graph_dirac_s3(n_max=3, rule="A")
    graphs["Dirac_S3_ruleA_n_max_3"] = A_dirac_A
    A_dirac_B, _ = graph_dirac_s3(n_max=3, rule="B")
    graphs["Dirac_S3_ruleB_n_max_3"] = A_dirac_B

    L_MAX = 20

    # Investigation I
    print("\n=== INVESTIGATION I: primitive-walk length distribution ===")
    inv_I = {}
    for name, A in graphs.items():
        print(f"  [{name}] V={A.shape[0]}, E={int(A.sum())//2}")
        inv_I[name] = investigation_I(A, L_max=L_MAX)
    output["investigation_I"] = inv_I

    # Investigation II
    print("\n=== INVESTIGATION II: explicit-formula fluctuating term ===")
    inv_II = {}
    for name, A in graphs.items():
        print(f"  [{name}]")
        inv_II[name] = investigation_II(A, L_max=L_MAX)
    output["investigation_II"] = inv_II

    # Investigation III (HEADLINE)
    print("\n=== INVESTIGATION III: chi_{-4} twisted Ihara zeta ===")
    inv_III = {}
    for name, A in graphs.items():
        print(f"  [{name}]")
        inv_III[name] = twisted_ihara_zeta_series(A, L_max=L_MAX)
    output["investigation_III"] = inv_III

    # Investigation IV: continuum limit via N_max sweep on S5
    print("\n=== INVESTIGATION IV: continuum limit for S5 ===")
    output["investigation_IV"] = {
        "s5_sweep": investigation_IV_s5(N_max_list=(2, 3)),
        "note": (
            "N_max=4,5 would be natural next sizes but V grows as "
            "(N+1)(N+2)(N+3)/6; the Hashimoto eigensolve becomes expensive. "
            "Reported here for N_max in {2, 3}."
        ),
    }

    # Save
    out_path = DATA_DIR / "riemann_limit_data.json"
    with out_path.open("w") as f:
        json.dump(output, f, indent=2, default=str)
    print(f"\nSaved data to {out_path}")
    return output


if __name__ == "__main__":
    run_all()
