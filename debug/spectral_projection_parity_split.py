#!/usr/bin/env python
"""Spectral projection parity split: decompose graph-native QED quantities
by even/odd Fock shell number of the internal electron.

Goal: look for connections between parity-split projection ratios
(graph-to-continuum) and the Dirichlet beta-function / Catalan's constant G.

The continuum D_even(s) - D_odd(s) = 2^{s-1} * (beta(s) - beta(s-2))
exposes Catalan G = beta(2) and beta(4).  If the graph-native QED
quantities similarly split by internal electron parity, the projection
exchange constant C may factorize by parity sector, revealing connections
to these Dirichlet L-values.

Sections
--------
1. Graph self-energy decomposed by even/odd internal electron n_fock
2. Graph vertex correction / F_2 decomposed by parity
3. Continuum self-energy split by parity of n_int
4. Ratio analysis (graph_even/continuum_even, graph_odd/continuum_odd)
5. PSLQ on parity-split ratios against {1, pi, pi^2, G, beta(4), zeta(3), sqrt(2), sqrt(3)}
6. Output / save to JSON

Author: GeoVac agentic workflow
Date: 2026-04-28
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

# ── mpmath for PSLQ ──────────────────────────────────────────────────
import mpmath
mpmath.mp.dps = 60  # 60 digits internal; PSLQ at 50

# ── GeoVac graph-native QED infrastructure ────────────────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

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
    KAPPA_SCALAR,
)
from geovac.qed_self_energy import (
    _so4_channel_count,
    _lambda_n,
    _g_n_dirac,
    _mu_q,
    _d_q_transverse,
    _vertex_allowed,
)

KAPPA = Rational(-1, 16)   # hopping parameter t = kappa = -1/16

# =====================================================================
# Helpers
# =====================================================================

def _build_graph_ingredients(n_max: int, t, exact: bool):
    """Build all shared graph-QED ingredients for a given n_max.

    Returns dict with V_mats, G_gamma, G_e, dirac_labels, N_dirac, E_fock.
    Uses numpy path for n_max >= 4 (sympy eigenvalue solver has bugs in L1).
    Force exact=False for photon propagator at n_max >= 4.
    """
    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    fock_data = build_fock_graph(n_max)
    entries, N_dirac, V_fock, E_fock = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats_sym = vertex_tensor_to_matrices(entries, N_dirac, E_fock)

    # Photon propagator -- for n_max >= 4, sympy eigenvalue solver has bugs
    # (produces spurious complex eigenvalues), so we bypass compute_photon_propagator
    # and build G_gamma directly via numpy pseudoinverse.
    if n_max >= 4:
        fock_data_for_photon = build_fock_graph(n_max)
        L1_np = np.array(fock_data_for_photon.L1.tolist(), dtype=float)
        G_gamma_direct = np.linalg.pinv(L1_np)
        # Create a minimal PhotonPropagatorData-like object
        class _PhotonStub:
            def __init__(self, G_gamma_numeric):
                self.G_gamma = None
                self.G_gamma_numeric = G_gamma_numeric
        photon_data = _PhotonStub(G_gamma_direct)
    else:
        photon_data = compute_photon_propagator(n_max, exact=exact)

    # Electron propagator
    op = DiracGraphOperator(n_max=n_max, t=t)
    if exact and n_max < 4:
        G_e_sym, _ = electron_propagator(op, exact=True, method='neumann')
        G_e_np = np.array(G_e_sym.tolist(), dtype=float)
    else:
        G_e_sym, _ = electron_propagator(op, exact=False, method='auto')
        G_e_np = G_e_sym if isinstance(G_e_sym, np.ndarray) else np.array(
            G_e_sym.tolist(), dtype=float
        )

    # Convert V_mats to numpy
    V_mats_np = [np.array(vm.tolist(), dtype=float) for vm in V_mats_sym]

    # Photon propagator to numpy
    if n_max < 4 and exact and photon_data.G_gamma is not None:
        G_gamma_np = np.array(photon_data.G_gamma.tolist(), dtype=float)
    else:
        G_gamma_np = photon_data.G_gamma_numeric

    return {
        'V_mats_np': V_mats_np,
        'G_gamma_np': G_gamma_np,
        'G_e_np': G_e_np,
        'dirac_labels': dirac_labels,
        'N_dirac': N_dirac,
        'E_fock': E_fock,
    }


def _parity_index_sets(dirac_labels):
    """Return sets of Dirac state indices with even and odd n_fock."""
    even_idx = [i for i, dl in enumerate(dirac_labels) if dl.n_fock % 2 == 0]
    odd_idx = [i for i, dl in enumerate(dirac_labels) if dl.n_fock % 2 == 1]
    return even_idx, odd_idx


# =====================================================================
# Section 1: Graph self-energy parity decomposition
# =====================================================================

def graph_self_energy_parity(n_max: int, t=KAPPA):
    """Decompose graph self-energy by parity of internal electron n_fock.

    Sigma[a,b] = sum_{e1,e2} G_gamma[e1,e2] * (V_e1 . V_e2^T)[a,b]

    The matrix product (V_e1 . V_e2^T)[a,b] = sum_c V_e1[a,c] * V_e2[b,c]
    sums over ALL intermediate Dirac states c.  We split this sum by
    parity of c's n_fock:

      Sigma_even[a,b] = sum_{e1,e2} G_gamma[e1,e2] *
                        sum_{c: n_fock(c) even} V_e1[a,c] * V_e2[b,c]

      Sigma_odd[a,b]  = sum_{e1,e2} G_gamma[e1,e2] *
                        sum_{c: n_fock(c) odd}  V_e1[a,c] * V_e2[b,c]
    """
    exact = (n_max <= 2)
    ing = _build_graph_ingredients(n_max, t, exact=exact)
    V_mats = ing['V_mats_np']
    G_gamma = ing['G_gamma_np']
    dirac_labels = ing['dirac_labels']
    N = ing['N_dirac']
    E = ing['E_fock']

    even_idx, odd_idx = _parity_index_sets(dirac_labels)

    Sigma_even = np.zeros((N, N))
    Sigma_odd = np.zeros((N, N))

    for e1 in range(E):
        for e2 in range(E):
            g_ee = G_gamma[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            # V_e1[:, even] . V_e2[:, even]^T  → parity-even internal sum
            Ve1_even = V_mats[e1][:, even_idx]  # (N, |even|)
            Ve2_even = V_mats[e2][:, even_idx]
            Sigma_even += g_ee * (Ve1_even @ Ve2_even.T)

            Ve1_odd = V_mats[e1][:, odd_idx]
            Ve2_odd = V_mats[e2][:, odd_idx]
            Sigma_odd += g_ee * (Ve1_odd @ Ve2_odd.T)

    Sigma_full = Sigma_even + Sigma_odd

    return {
        'n_max': n_max,
        'Tr_Sigma_full': float(np.trace(Sigma_full)),
        'Tr_Sigma_even': float(np.trace(Sigma_even)),
        'Tr_Sigma_odd': float(np.trace(Sigma_odd)),
        'ratio_even_odd': (
            float(np.trace(Sigma_even) / np.trace(Sigma_odd))
            if abs(np.trace(Sigma_odd)) > 1e-15 else None
        ),
        'frac_even': float(np.trace(Sigma_even) / np.trace(Sigma_full)),
        'frac_odd': float(np.trace(Sigma_odd) / np.trace(Sigma_full)),
        'N_even': len(even_idx),
        'N_odd': len(odd_idx),
    }


# =====================================================================
# Section 2: Graph vertex correction / F_2 parity decomposition
# =====================================================================

def graph_vertex_parity(n_max: int, t=KAPPA):
    """Decompose graph vertex correction by parity of internal electron.

    Lambda[a,b] = sum_{e1,e2} G_gamma[e1,e2] * (V_e1 . G_e . V_e2^T)[a,b]

    The internal electron enters through G_e.  We split:

      Lambda_even[a,b] = sum_{e1,e2} G_gamma[e1,e2] *
                         V_e1[:, even] . G_e[even, even] . V_e2[:, even]^T

      Lambda_odd[a,b]  = analogous with odd indices

    Additionally there are cross-parity terms
      G_e[even, odd] and G_e[odd, even]
    which arise from the off-diagonal hopping t.  We track these too.
    """
    exact = (n_max <= 2)
    ing = _build_graph_ingredients(n_max, t, exact=exact)
    V_mats = ing['V_mats_np']
    G_gamma = ing['G_gamma_np']
    G_e = ing['G_e_np']
    dirac_labels = ing['dirac_labels']
    N = ing['N_dirac']
    E = ing['E_fock']

    even_idx, odd_idx = _parity_index_sets(dirac_labels)
    ev = np.array(even_idx)
    od = np.array(odd_idx)

    # G_e sub-blocks
    G_e_ee = G_e[np.ix_(ev, ev)]   # even-even
    G_e_oo = G_e[np.ix_(od, od)]   # odd-odd
    G_e_eo = G_e[np.ix_(ev, od)]   # even-odd cross
    G_e_oe = G_e[np.ix_(od, ev)]   # odd-even cross

    Lambda_even = np.zeros((N, N))   # purely even internal
    Lambda_odd = np.zeros((N, N))    # purely odd internal
    Lambda_cross = np.zeros((N, N))  # cross-parity internal

    for e1 in range(E):
        for e2 in range(E):
            g_ee = G_gamma[e1, e2]
            if abs(g_ee) < 1e-15:
                continue

            Ve1_ev = V_mats[e1][:, ev]   # (N, |even|)
            Ve2_ev = V_mats[e2][:, ev]
            Ve1_od = V_mats[e1][:, od]
            Ve2_od = V_mats[e2][:, od]

            # Even internal: V_e1[:, ev] . G_e[ev,ev] . V_e2[:, ev]^T
            Lambda_even += g_ee * (Ve1_ev @ G_e_ee @ Ve2_ev.T)

            # Odd internal: V_e1[:, od] . G_e[od,od] . V_e2[:, od]^T
            Lambda_odd += g_ee * (Ve1_od @ G_e_oo @ Ve2_od.T)

            # Cross terms (even->odd and odd->even)
            Lambda_cross += g_ee * (Ve1_ev @ G_e_eo @ Ve2_od.T)
            Lambda_cross += g_ee * (Ve1_od @ G_e_oe @ Ve2_ev.T)

    Lambda_full = Lambda_even + Lambda_odd + Lambda_cross

    # Bare vertex and F_2 extraction
    V_bare = sum(V_mats)
    tr_norm = np.trace(V_bare @ G_e)

    F2_full = np.trace(Lambda_full) / tr_norm if abs(tr_norm) > 1e-15 else 0.0
    F2_even = np.trace(Lambda_even) / tr_norm if abs(tr_norm) > 1e-15 else 0.0
    F2_odd = np.trace(Lambda_odd) / tr_norm if abs(tr_norm) > 1e-15 else 0.0
    F2_cross = np.trace(Lambda_cross) / tr_norm if abs(tr_norm) > 1e-15 else 0.0

    return {
        'n_max': n_max,
        'Tr_Lambda_full': float(np.trace(Lambda_full)),
        'Tr_Lambda_even': float(np.trace(Lambda_even)),
        'Tr_Lambda_odd': float(np.trace(Lambda_odd)),
        'Tr_Lambda_cross': float(np.trace(Lambda_cross)),
        'F2_full': float(F2_full),
        'F2_even': float(F2_even),
        'F2_odd': float(F2_odd),
        'F2_cross': float(F2_cross),
        'ratio_even_odd_Lambda': (
            float(np.trace(Lambda_even) / np.trace(Lambda_odd))
            if abs(np.trace(Lambda_odd)) > 1e-15 else None
        ),
        'frac_even': float(np.trace(Lambda_even) / np.trace(Lambda_full)),
        'frac_odd': float(np.trace(Lambda_odd) / np.trace(Lambda_full)),
        'frac_cross': float(np.trace(Lambda_cross) / np.trace(Lambda_full)),
    }


# =====================================================================
# Section 3: Continuum self-energy split by parity of n_int
# =====================================================================

def continuum_self_energy_parity(n_ext: int, n_max: int):
    """Split the continuum one-loop self-energy by parity of n_int.

    Sigma(n_ext) = Sigma_even(n_ext) + Sigma_odd(n_ext)

    where Sigma_even sums only over even n_int, Sigma_odd over odd n_int.
    """
    sigma_even = mpmath.mpf(0)
    sigma_odd = mpmath.mpf(0)

    for n_int in range(n_max + 1):
        g_int = _g_n_dirac(n_int)
        lam_int = _lambda_n(n_int)
        lam_int_pow = lam_int ** 4   # s_e=2 → 2*s_e=4

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

            contrib = mpmath.mpf(W) * g_int * d_T / (lam_int_pow * mu)

            if n_int % 2 == 0:
                sigma_even += contrib
            else:
                sigma_odd += contrib

    sigma_full = sigma_even + sigma_odd
    return {
        'n_ext': n_ext,
        'n_max': n_max,
        'Sigma_full': float(sigma_full),
        'Sigma_even': float(sigma_even),
        'Sigma_odd': float(sigma_odd),
        'ratio_even_odd': (
            float(sigma_even / sigma_odd)
            if abs(sigma_odd) > 1e-30 else None
        ),
    }


# =====================================================================
# Section 5: PSLQ on parity-split ratios
# =====================================================================

def pslq_on_ratios(ratios_dict: dict):
    """Run PSLQ at 50 digits on key ratios against a basis of constants.

    Basis: {1, pi, pi^2, G, beta(4), zeta(3), sqrt(2), sqrt(3)}
    where G = Catalan's constant = beta(2).
    """
    mpmath.mp.dps = 60

    # Build basis constants
    pi_val = mpmath.pi
    G_val = mpmath.catalan        # Catalan's constant G = beta(2)
    # beta(4) = pi^4/1536 + G*pi^2/16 - 7*zeta(3)/16  ... compute directly
    # Actually: beta(4) = sum_{k>=0} (-1)^k / (2k+1)^4 = pi^4/768 * ...
    # Safest: compute via Hurwitz
    beta4 = mpmath.mpf(1) / mpmath.mpf(4)**4 * (
        mpmath.hurwitz(4, mpmath.mpf(1)/4) - mpmath.hurwitz(4, mpmath.mpf(3)/4)
    )
    zeta3 = mpmath.zeta(3)
    sqrt2 = mpmath.sqrt(2)
    sqrt3 = mpmath.sqrt(3)

    basis = [
        mpmath.mpf(1),   # 1
        pi_val,           # pi
        pi_val**2,        # pi^2
        G_val,            # G = Catalan
        beta4,            # beta(4)
        zeta3,            # zeta(3)
        sqrt2,            # sqrt(2)
        sqrt3,            # sqrt(3)
    ]
    basis_names = ['1', 'pi', 'pi^2', 'G', 'beta(4)', 'zeta(3)', 'sqrt(2)', 'sqrt(3)']

    results = {}
    for key, val in ratios_dict.items():
        if val is None or not np.isfinite(val):
            results[key] = {'value': val, 'pslq': None}
            continue
        target = mpmath.mpf(val)
        vec = basis + [target]

        try:
            rel = mpmath.pslq(vec, maxcoeff=10000, tol=mpmath.mpf(10)**(-40))
        except Exception:
            rel = None

        if rel is not None:
            # The relation is: sum_i rel[i] * vec[i] = 0
            # So: rel[-1]*target + sum_{i<N} rel[i]*basis[i] = 0
            # => target = -sum_{i<N} (rel[i]/rel[-1]) * basis[i]
            coeff_target = rel[-1]
            if coeff_target != 0:
                expr_parts = []
                for i, name in enumerate(basis_names):
                    c = -rel[i]  # note the sign flip
                    if c != 0:
                        expr_parts.append(f"({c}/{coeff_target})*{name}")
                expr_str = " + ".join(expr_parts) if expr_parts else "0"
            else:
                expr_str = "degenerate (coeff on target = 0)"

            results[key] = {
                'value': float(val),
                'pslq_relation': [int(r) for r in rel],
                'expression': expr_str,
                'residual': float(abs(sum(r * v for r, v in zip(rel, vec)))),
            }
        else:
            results[key] = {
                'value': float(val),
                'pslq': None,
                'note': 'No PSLQ relation found at 50-digit precision',
            }

    return results


# =====================================================================
# Main computation
# =====================================================================

def main():
    t0 = time.time()
    results = {}

    # ------------------------------------------------------------------
    # Section 1: Graph self-energy parity decomposition
    # ------------------------------------------------------------------
    print("=" * 70)
    print("Section 1: Graph self-energy parity decomposition")
    print("=" * 70)

    se_results = {}
    for nmax in [2, 3, 4, 5]:
        print(f"  n_max={nmax} ... ", end="", flush=True)
        t1 = time.time()
        res = graph_self_energy_parity(nmax)
        dt = time.time() - t1
        print(f"done ({dt:.1f}s)")
        print(f"    Tr(Sigma_full)  = {res['Tr_Sigma_full']:.10f}")
        print(f"    Tr(Sigma_even)  = {res['Tr_Sigma_even']:.10f}")
        print(f"    Tr(Sigma_odd)   = {res['Tr_Sigma_odd']:.10f}")
        print(f"    ratio even/odd  = {res['ratio_even_odd']}")
        print(f"    frac even       = {res['frac_even']:.6f}")
        print(f"    frac odd        = {res['frac_odd']:.6f}")
        se_results[str(nmax)] = res
    results['graph_self_energy_parity'] = se_results

    # ------------------------------------------------------------------
    # Section 2: Graph vertex correction / F_2 parity decomposition
    # ------------------------------------------------------------------
    print()
    print("=" * 70)
    print("Section 2: Graph vertex correction / F_2 parity decomposition")
    print("=" * 70)

    vc_results = {}
    for nmax in [2, 3, 4, 5]:
        print(f"  n_max={nmax} ... ", end="", flush=True)
        t1 = time.time()
        res = graph_vertex_parity(nmax)
        dt = time.time() - t1
        print(f"done ({dt:.1f}s)")
        print(f"    F2_full   = {res['F2_full']:.10f}")
        print(f"    F2_even   = {res['F2_even']:.10f}")
        print(f"    F2_odd    = {res['F2_odd']:.10f}")
        print(f"    F2_cross  = {res['F2_cross']:.10f}")
        print(f"    ratio Lambda_even/Lambda_odd = {res['ratio_even_odd_Lambda']}")
        print(f"    frac even   = {res['frac_even']:.6f}")
        print(f"    frac odd    = {res['frac_odd']:.6f}")
        print(f"    frac cross  = {res['frac_cross']:.6f}")
        vc_results[str(nmax)] = res
    results['graph_vertex_parity'] = vc_results

    # ------------------------------------------------------------------
    # Section 3: Continuum self-energy parity split
    # ------------------------------------------------------------------
    print()
    print("=" * 70)
    print("Section 3: Continuum self-energy parity split")
    print("=" * 70)

    cont_results = {}
    for n_ext in [0, 1, 2, 3]:
        for n_max_cont in [20, 50]:
            key = f"n_ext={n_ext}_nmax={n_max_cont}"
            print(f"  {key} ... ", end="", flush=True)
            t1 = time.time()
            res = continuum_self_energy_parity(n_ext, n_max_cont)
            dt = time.time() - t1
            print(f"done ({dt:.1f}s)")
            print(f"    Sigma_full  = {res['Sigma_full']:.10f}")
            print(f"    Sigma_even  = {res['Sigma_even']:.10f}")
            print(f"    Sigma_odd   = {res['Sigma_odd']:.10f}")
            print(f"    ratio even/odd = {res['ratio_even_odd']}")
            cont_results[key] = res
    results['continuum_self_energy_parity'] = cont_results

    # ------------------------------------------------------------------
    # Section 4: Ratio analysis
    # ------------------------------------------------------------------
    print()
    print("=" * 70)
    print("Section 4: Ratio analysis (graph / continuum parity sectors)")
    print("=" * 70)

    ratio_results = {}

    # For each n_max where both graph and continuum exist, compare traces
    # Graph Tr(Sigma) is summed over all external states; continuum is per n_ext.
    # The most meaningful comparison: graph trace vs continuum summed over n_ext.

    for nmax in [2, 3, 4, 5]:
        graph_se = se_results[str(nmax)]
        graph_vc = vc_results[str(nmax)]

        # Continuum at same n_max
        cont_sigma_even_sum = mpmath.mpf(0)
        cont_sigma_odd_sum = mpmath.mpf(0)
        for n_ext in range(nmax + 1):
            for n_int in range(nmax + 1):
                g_int = _g_n_dirac(n_int)
                lam_int = _lambda_n(n_int)
                lam_pow = lam_int ** 4
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
                    contrib = mpmath.mpf(W) * g_int * d_T / (lam_pow * mu)
                    if n_int % 2 == 0:
                        cont_sigma_even_sum += contrib
                    else:
                        cont_sigma_odd_sum += contrib

        cont_sigma_full_sum = float(cont_sigma_even_sum + cont_sigma_odd_sum)
        cont_even_f = float(cont_sigma_even_sum)
        cont_odd_f = float(cont_sigma_odd_sum)

        ratio_even = (
            graph_se['Tr_Sigma_even'] / cont_even_f
            if abs(cont_even_f) > 1e-15 else None
        )
        ratio_odd = (
            graph_se['Tr_Sigma_odd'] / cont_odd_f
            if abs(cont_odd_f) > 1e-15 else None
        )
        ratio_full = (
            graph_se['Tr_Sigma_full'] / cont_sigma_full_sum
            if abs(cont_sigma_full_sum) > 1e-15 else None
        )

        # Ratio of ratios (most interesting for Catalan connection)
        ratio_of_ratios = (
            ratio_even / ratio_odd
            if ratio_even is not None and ratio_odd is not None
            and abs(ratio_odd) > 1e-15 else None
        )

        entry = {
            'n_max': nmax,
            'graph_Tr_Sigma_even': graph_se['Tr_Sigma_even'],
            'graph_Tr_Sigma_odd': graph_se['Tr_Sigma_odd'],
            'cont_Sigma_even_sum': cont_even_f,
            'cont_Sigma_odd_sum': cont_odd_f,
            'ratio_even': ratio_even,
            'ratio_odd': ratio_odd,
            'ratio_full': ratio_full,
            'ratio_of_ratios': ratio_of_ratios,
            'graph_even_odd_ratio': graph_se['ratio_even_odd'],
            'cont_even_odd_ratio': (
                cont_even_f / cont_odd_f
                if abs(cont_odd_f) > 1e-15 else None
            ),
            # F2 parity decomposition
            'F2_even': graph_vc['F2_even'],
            'F2_odd': graph_vc['F2_odd'],
            'F2_cross': graph_vc['F2_cross'],
            'F2_even_odd_ratio': (
                graph_vc['F2_even'] / graph_vc['F2_odd']
                if graph_vc['F2_odd'] is not None
                and abs(graph_vc['F2_odd']) > 1e-15 else None
            ),
        }
        ratio_results[str(nmax)] = entry

        print(f"\n  n_max = {nmax}:")
        print(f"    graph Tr(Sigma_even) = {graph_se['Tr_Sigma_even']:.10f}")
        print(f"    graph Tr(Sigma_odd)  = {graph_se['Tr_Sigma_odd']:.10f}")
        print(f"    cont  Sigma_even_sum = {cont_even_f:.10f}")
        print(f"    cont  Sigma_odd_sum  = {cont_odd_f:.10f}")
        print(f"    ratio_even (graph/cont) = {ratio_even}")
        print(f"    ratio_odd  (graph/cont) = {ratio_odd}")
        print(f"    ratio_full (graph/cont) = {ratio_full}")
        print(f"    ratio_of_ratios         = {ratio_of_ratios}")
        print(f"    graph even/odd          = {graph_se['ratio_even_odd']}")
        print(f"    cont  even/odd          = {entry['cont_even_odd_ratio']}")
        print(f"    F2_even/F2_odd          = {entry['F2_even_odd_ratio']}")

    results['ratio_analysis'] = ratio_results

    # ------------------------------------------------------------------
    # Section 5: PSLQ on parity-split ratios
    # ------------------------------------------------------------------
    print()
    print("=" * 70)
    print("Section 5: PSLQ on parity-split ratios")
    print("=" * 70)

    # Collect key ratios for PSLQ
    pslq_targets = {}
    for nmax_str, entry in ratio_results.items():
        prefix = f"nmax{nmax_str}"
        if entry['ratio_even'] is not None:
            pslq_targets[f"{prefix}_ratio_even"] = entry['ratio_even']
        if entry['ratio_odd'] is not None:
            pslq_targets[f"{prefix}_ratio_odd"] = entry['ratio_odd']
        if entry['ratio_of_ratios'] is not None:
            pslq_targets[f"{prefix}_ratio_of_ratios"] = entry['ratio_of_ratios']
        if entry.get('graph_even_odd_ratio') is not None:
            pslq_targets[f"{prefix}_graph_even_odd"] = entry['graph_even_odd_ratio']
        if entry.get('cont_even_odd_ratio') is not None:
            pslq_targets[f"{prefix}_cont_even_odd"] = entry['cont_even_odd_ratio']
        if entry.get('F2_even_odd_ratio') is not None:
            pslq_targets[f"{prefix}_F2_even_odd"] = entry['F2_even_odd_ratio']

    pslq_results = pslq_on_ratios(pslq_targets)

    for key, val in pslq_results.items():
        print(f"  {key}: ", end="")
        if val.get('pslq_relation') is not None:
            print(f"IDENTIFIED: {val['expression']}  (residual {val['residual']:.2e})")
        elif val.get('pslq') is None and val.get('note'):
            print(f"value={val.get('value', '?'):.10g}  -- {val['note']}")
        else:
            print(f"value={val.get('value', '?')}")

    results['pslq'] = {}
    for k, v in pslq_results.items():
        # Make JSON-serializable
        entry = {}
        for kk, vv in v.items():
            if isinstance(vv, (mpmath.mpf, np.floating)):
                entry[kk] = float(vv)
            elif isinstance(vv, (mpmath.mpc, complex)):
                entry[kk] = str(vv)
            else:
                entry[kk] = vv
        results['pslq'][k] = entry

    # ------------------------------------------------------------------
    # Section 6: Save
    # ------------------------------------------------------------------
    outpath = Path(__file__).parent / "data" / "spectral_projection_parity_split.json"
    outpath.parent.mkdir(parents=True, exist_ok=True)

    # Make everything JSON-serializable
    def _clean(obj):
        if isinstance(obj, dict):
            return {k: _clean(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [_clean(v) for v in obj]
        elif isinstance(obj, (mpmath.mpf, np.floating)):
            return float(obj)
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, (mpmath.mpc, complex)):
            return str(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    with open(outpath, 'w') as f:
        json.dump(_clean(results), f, indent=2)

    dt_total = time.time() - t0
    print(f"\n{'=' * 70}")
    print(f"Total runtime: {dt_total:.1f}s")
    print(f"Results saved to: {outpath}")
    print(f"{'=' * 70}")

    return results


if __name__ == "__main__":
    results = main()
