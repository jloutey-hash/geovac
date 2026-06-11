"""
VP-1: Vector-photon promotion of graph-native QED vertex correction.
=====================================================================

Sprint goal: test whether replacing the SCALAR photon (Fock-edge 1-cochain)
in the graph-native vertex correction with a VECTOR photon ((q, m_q) modes
with Wigner 3j vertex coupling) closes the C × F2 negative result.

Background (from CLAUDE.md):
  - Graph-native scalar QED at n_max=2 gives F2 = 5*sqrt(2)/3 ≈ 2.357.
  - Power-law convergence F2 ~ n^(-0.57) -> 0 in the continuum limit.
  - Projection constant C is rational: C(n_max=3) = 50471424/1779441125.
  - C × F2_scalar diverges from alpha/(2*pi) ≈ 1.16e-3.

Hypothesis: the missing factor is the photon's vector structure
(angular momentum quantum numbers (q, m_q)). With vector photons:

  Lambda_total[a,c] = sum_{b, q, m_q} V(a,b,q,m_q) * G_e^graph[b,b]
                                       * G_gamma_vec(q) * V(c,b,q,m_q)

where:
  - G_e^graph: graph-native Dirac electron propagator (from
    DiracGraphOperator), exact rational.  Uses Camporesi-Higuchi
    eigenvalues |lambda_n| = n + 1/2 in Fock convention.  At t=0 this
    is diag(1/lambda_n) (sign included).
  - V(a, b, q, m_q): Wigner 3j vertex from vector_qed.dirac_vertex_coupling
    (Dirac states a,b ∈ DiracLabel; photon (q, m_q); E1 parity).
  - G_gamma_vec(q) = 1 / [q(q+2)]  (vector Laplacian on S^3, Hodge-1).

F2 extraction (graph-native style):

  F2 = Tr(Lambda_total) / Tr(V_bare * G_e)

where V_bare = sum_{q, m_q} V(:, :, q, m_q) * G_gamma_vec(q)? — this is
ambiguous in the vector case.  We use the EDGE-SUMMED bare vertex
analog: V_bare[a,b] = sum_{q,m_q} V(a,b,q,m_q), the same convention as
graph_qed_self_energy.extract_anomalous_moment.  This gives a
dimensionless F2 directly comparable to the graph-native value.

Selection rules: by construction the Dirac+vector vertex preserves all
8 continuum QED selection rules (vector_qed.py Part 2, Furry recovered
via Dirac spinor phase).

Output:
  - F2_vector(n_max) at n_max ∈ {2, 3, 4}
  - C × F2_vector compared to alpha/(2*pi) = 1.16140973e-3
  - scaling exponent of F2_vector vs F2_scalar's -0.57

Author: VP-1 sprint 2026-05-02
"""
from __future__ import annotations

import json
import time
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import sympy as sp
from sympy import Matrix, Rational, sqrt as sp_sqrt
from sympy.physics.wigner import wigner_3j

# Graph-native infrastructure (electron propagator)
from geovac.graph_qed_propagator import (
    DiracGraphOperator,
    electron_propagator as graph_electron_propagator,
)

# Vector QED infrastructure (Dirac states, vector photon vertex)
from geovac.vector_qed import (
    build_dirac_electron_states,
    build_photon_modes,
    vector_photon_propagator,
    dirac_vertex_coupling,
    build_dirac_vertex_tensor,
)

from geovac.dirac_matrix_elements import DiracLabel, kappa_to_l


ALPHA = 1.0 / 137.035999084
ALPHA_OVER_2PI = ALPHA / (2.0 * np.pi)  # ≈ 1.16140973e-3

# Existing graph-native projection constants (from gn7_continuum_bridge)
PROJECTION_CONSTANT_NMAX3 = Fraction(50471424, 1779441125)
# n_max=2: continuum VP=0, so C is undefined there

# ---------------------------------------------------------------------------
# Vector vertex tensor in EXACT sympy (for n_max=2,3 verification)
# ---------------------------------------------------------------------------

def dirac_vertex_coupling_exact(
    a: DiracLabel,
    b: DiracLabel,
    q: int,
    m_q: int,
) -> sp.Expr:
    """Exact sympy version of vector_qed.dirac_vertex_coupling.

    Returns a sympy expression (typically of the form rational * sqrt(rational) / sqrt(pi)).
    Note: pi enters via the spherical harmonic normalization sqrt(1/(4*pi)).
    """
    # Rule 0: diagonal coupling vanishes (Dirac spinor phase constraint)
    if a == b:
        return sp.S.Zero

    l_a = kappa_to_l(a.kappa)
    l_b = kappa_to_l(b.kappa)

    # Parity selection: l_a + l_b + q must be odd
    if (l_a + l_b + q) % 2 == 0:
        return sp.S.Zero

    # l-triangle
    if q < abs(l_a - l_b) or q > l_a + l_b:
        return sp.S.Zero

    j_a = a.j
    j_b = b.j
    m_j_a = a.m_j
    m_j_b = b.m_j

    # Magnetic conservation
    if -m_j_a + m_q + m_j_b != 0:
        return sp.S.Zero

    # j-triangle
    if q < abs(j_a - j_b) or q > j_a + j_b:
        return sp.S.Zero

    threej = wigner_3j(j_a, q, j_b, -m_j_a, m_q, m_j_b)
    if threej == 0:
        return sp.S.Zero

    prefactor_sq = (2 * j_a + 1) * (2 * q + 1) * (2 * j_b + 1) / (4 * sp.pi)
    prefactor = sp.sqrt(prefactor_sq)

    phase_exp = j_a - m_j_a
    phase = (-1) ** int(phase_exp)

    return phase * prefactor * threej


def build_dirac_vector_vertex_tensor_exact(
    states: List[DiracLabel],
    modes: List[Tuple[int, int]],
) -> List[List[List[sp.Expr]]]:
    """Build V[a, b, k] in exact sympy.

    Returns a 3D nested list V[a][b][k].
    """
    N_e = len(states)
    N_g = len(modes)
    V = [[[sp.S.Zero for _ in range(N_g)] for _ in range(N_e)]
         for _ in range(N_e)]
    for i, a in enumerate(states):
        for j, b in enumerate(states):
            for k, (q, m_q) in enumerate(modes):
                V[i][j][k] = dirac_vertex_coupling_exact(a, b, q, m_q)
    return V


# ---------------------------------------------------------------------------
# Graph-native electron propagator on Dirac graph (atomic mode)
# ---------------------------------------------------------------------------

def build_dirac_electron_G_atomic(
    n_max: int,
    t: sp.Expr = Rational(0),
    exact: bool = True,
):
    """Build G_e on the Dirac graph (l < n, atomic mode), exact sympy.

    Returns
    -------
    (G_e_matrix, dirac_labels, N_dirac)
    """
    op = DiracGraphOperator(n_max=n_max, t=t, mode='atomic')
    G_e, _ = graph_electron_propagator(op, exact=exact)
    return G_e, op.labels, op.N


def graph_dirac_label_index(
    op_labels: List[DiracLabel],
    state: DiracLabel,
) -> int:
    """Find the index of a DiracLabel in the DiracGraphOperator label list.

    The graph propagator labels and vector_qed states should be the same
    set (atomic mode + iter_dirac_labels) but ordering may differ.
    """
    for i, lab in enumerate(op_labels):
        if (lab.n_fock == state.n_fock
                and lab.kappa == state.kappa
                and lab.two_m_j == state.two_m_j):
            return i
    raise KeyError(f"DiracLabel {state} not found in op labels")


# ---------------------------------------------------------------------------
# Vector vertex correction on graph-native pipeline
# ---------------------------------------------------------------------------

def compute_vector_F2_numpy(
    n_max: int,
    q_max: Optional[int] = None,
    exclude_gs_internal: bool = False,
) -> Dict:
    """Compute F2_vector at n_max using numpy (fast, all n_max).

    Lambda_total[a, c] = sum_{b, k} V[a,b,k] * G_e[b,b] * G_gamma(q_k) * V[c,b,k]

    F2 = Tr(Lambda_total) / Tr(V_bare * G_e)
    where V_bare[a,b] = sum_{k} V[a,b,k] (edge-summed vertex,
    matching graph-native extraction in graph_qed_self_energy).

    Parameters
    ----------
    n_max : int
        Electron truncation.
    q_max : int, optional
        Photon truncation (default = n_max - 1, matching vector_qed).
    exclude_gs_internal : bool
        If True, exclude n=1 from internal sums.
    """
    if q_max is None:
        q_max = n_max - 1
    if q_max < 1:
        q_max = 1

    states_vec = build_dirac_electron_states(n_max)
    modes = build_photon_modes(q_max)
    N_e = len(states_vec)
    N_g = len(modes)

    # Vertex tensor (numpy)
    V = build_dirac_vertex_tensor(states_vec, modes)

    # Graph-native electron propagator: build G_e_diag[i] = G_e[i, i] in
    # the order matching states_vec.
    op = DiracGraphOperator(n_max=n_max, t=Rational(0), mode='atomic')
    op_labels = op.labels

    # Graph-native G_e at t=0: diag(1/lambda_n) with sign(chirality)
    # Use the eigenvalues directly (already exact rationals)
    G_e_op = np.zeros(N_e)
    for i, lab in enumerate(states_vec):
        # Find this label in op_labels
        idx = graph_dirac_label_index(op_labels, lab)
        lam_op = float(op._eigenvalues[idx])  # signed eigenvalue
        G_e_op[i] = 1.0 / lam_op

    # Photon propagator
    G_g = np.zeros(N_g)
    for k, (q, m_q) in enumerate(modes):
        G_g[k] = vector_photon_propagator(q)

    # Internal mask
    if exclude_gs_internal:
        internal_mask = np.array([s.n_fock >= 2 for s in states_vec], dtype=bool)
    else:
        internal_mask = np.ones(N_e, dtype=bool)

    # Lambda_total = sum_{b internal, k} V[:,b,k] * G_e[b] * G_g[k] * V[:,b,k]^T
    # Note: this gives Lambda_total[a,c] = sum_{b,k} V[a,b,k] G_e[b] G_g[k] V[c,b,k]
    # which is the symmetric vertex-correction analog.
    Lambda_total = np.zeros((N_e, N_e))
    for b_idx in range(N_e):
        if not internal_mask[b_idx]:
            continue
        ge_b = G_e_op[b_idx]
        for k_idx in range(N_g):
            gg_k = G_g[k_idx]
            weight = ge_b * gg_k
            if abs(weight) < 1e-30:
                continue
            v_col = V[:, b_idx, k_idx]
            Lambda_total += weight * np.outer(v_col, v_col)

    tr_lambda = float(np.trace(Lambda_total))

    # V_bare[a, b] = sum_k V[a, b, k]  (edge-sum convention)
    V_bare = V.sum(axis=2)

    # Tr(V_bare * G_e) = sum_a (V_bare * G_e)[a, a] where G_e is diagonal
    # diag(G_e_op).  V_bare * G_e_diag is (V_bare with each col j scaled by G_e[j])
    # then trace = sum_i V_bare[i,i] * G_e[i].
    tr_norm = float(np.sum(np.diag(V_bare) * G_e_op))

    # Note: V_bare diagonal is zero by Dirac spinor phase constraint
    # (V(a,a,q,m_q)=0 by construction in vector_qed Rule 0).  So this
    # normalization is identically zero — Furry's theorem at the bare
    # vertex level. We need a different normalization for F2.

    # Alternative: use the natural vector-photon normalization where the
    # ANOMALOUS moment F2 is extracted as the ratio of the vertex
    # correction trace to the bare vertex coupling structure summed
    # over allowed channels.  Following continuum QED practice:
    #   F2_extracted = Tr(Lambda_total) / N_dim
    # where N_dim is a normalization choice. The cleanest choice for
    # comparison with the scalar graph-native F2 (which used Tr(V_bare*G_e))
    # is to use the SAME EXTRACTION FORMULA on the vector vertex, but
    # if V_bare diagonal vanishes, fall back to Tr(Lambda_total) / Tr(G_e).

    tr_G_e = float(np.sum(G_e_op))

    # Primary F2: use Tr(Lambda) / Tr(G_e) (analog of bare normalization
    # when diagonal V vanishes — gives a dimensionless number with
    # the right units).
    F2_vector = tr_lambda / tr_G_e if abs(tr_G_e) > 1e-30 else 0.0

    # Secondary normalization options for triangulation
    # (a) by Frobenius norm of V (channel-summed): provides per-channel
    # average vertex strength
    V_frob_sq = float(np.sum(V ** 2))
    F2_vector_frob = tr_lambda / V_frob_sq if V_frob_sq > 1e-30 else 0.0

    # (b) by N_e (matrix dimension): "average diagonal"
    F2_vector_avg_diag = tr_lambda / N_e if N_e > 0 else 0.0

    # GS structural zero check
    gs_indices = [i for i, s in enumerate(states_vec) if s.n_fock == 1]
    if gs_indices:
        gs_block_max = float(np.max(np.abs(Lambda_total[np.ix_(gs_indices, gs_indices)])))
    else:
        gs_block_max = None

    return {
        'n_max': n_max,
        'q_max': q_max,
        'N_electron': N_e,
        'N_photon': N_g,
        'Lambda_trace': tr_lambda,
        'Lambda_frobenius': float(np.linalg.norm(Lambda_total, 'fro')),
        'Lambda_max_entry': float(np.max(np.abs(Lambda_total))),
        'Lambda_GS_block_max': gs_block_max,
        'V_bare_diag_sum': float(np.sum(np.diag(V_bare))),
        'V_bare_diag_is_zero': abs(np.sum(np.diag(V_bare))) < 1e-12,
        'Tr_G_e': tr_G_e,
        'F2_vector_vs_TrGe': F2_vector,
        'F2_vector_vs_VFrob': F2_vector_frob,
        'F2_vector_vs_Ne': F2_vector_avg_diag,
        'V_Frobenius_sq': V_frob_sq,
        'vertex_nonzero': int(np.count_nonzero(V)),
        'vertex_total': V.size,
        'vertex_sparsity': 1.0 - np.count_nonzero(V) / V.size,
        'exclude_gs_internal': exclude_gs_internal,
    }


def compute_vector_F2_exact(n_max: int) -> Dict:
    """Exact sympy version of compute_vector_F2_numpy at small n_max.

    Returns the same dict but with sympy expressions (string-formatted)
    for F2_vector. Only feasible for n_max=2, possibly 3.
    """
    q_max = max(1, n_max - 1)

    states_vec = build_dirac_electron_states(n_max)
    modes = build_photon_modes(q_max)
    N_e = len(states_vec)
    N_g = len(modes)

    V = build_dirac_vector_vertex_tensor_exact(states_vec, modes)

    op = DiracGraphOperator(n_max=n_max, t=Rational(0), mode='atomic')
    op_labels = op.labels

    # G_e_op[i] = 1/lambda_n[op_idx]
    G_e_op = []
    for lab in states_vec:
        idx = graph_dirac_label_index(op_labels, lab)
        lam = op._eigenvalues[idx]
        G_e_op.append(Rational(1) / lam)

    # Photon propagator G_g[k] = 1/[q(q+2)]
    G_g = [Rational(1, q * (q + 2)) for (q, m_q) in modes]

    # Lambda_total[a,c] = sum_{b,k} V[a][b][k] * G_e_op[b] * G_g[k] * V[c][b][k]
    Lambda = sp.zeros(N_e, N_e)
    for b_idx in range(N_e):
        ge_b = G_e_op[b_idx]
        for k_idx in range(N_g):
            gg_k = G_g[k_idx]
            weight = ge_b * gg_k
            for a in range(N_e):
                Vab = V[a][b_idx][k_idx]
                if Vab == 0:
                    continue
                for c in range(N_e):
                    Vcb = V[c][b_idx][k_idx]
                    if Vcb == 0:
                        continue
                    Lambda[a, c] += weight * Vab * Vcb

    # Simplify entries
    Lambda_simp = sp.zeros(N_e, N_e)
    for i in range(N_e):
        for j in range(N_e):
            Lambda_simp[i, j] = sp.simplify(Lambda[i, j])

    tr_lambda = sp.simplify(Lambda_simp.trace())
    tr_lambda_float = float(tr_lambda)

    Tr_G_e = sp.simplify(sum(G_e_op))
    Tr_G_e_float = float(Tr_G_e)

    F2 = sp.simplify(tr_lambda / Tr_G_e) if Tr_G_e != 0 else sp.S.Zero
    F2_float = float(F2)

    return {
        'n_max': n_max,
        'q_max': q_max,
        'N_electron': N_e,
        'N_photon': N_g,
        'Lambda_trace_exact': str(tr_lambda),
        'Lambda_trace_float': tr_lambda_float,
        'Tr_G_e_exact': str(Tr_G_e),
        'Tr_G_e_float': Tr_G_e_float,
        'F2_vector_exact': str(F2),
        'F2_vector_float': F2_float,
    }


# ---------------------------------------------------------------------------
# Comparison and scaling
# ---------------------------------------------------------------------------

def F2_scalar_reference(n_max: int) -> float:
    """Graph-native scalar F2 reference values (from CLAUDE.md)."""
    table = {
        2: 5 * np.sqrt(2) / 3,
        3: 1.873,
        4: 1.589,
        5: 1.396,
        6: 1.253,
    }
    return table.get(n_max)


def projection_constant_C(n_max: int) -> Optional[float]:
    """Projection constant C from gn7 continuum bridge.

    C is well-defined only for n_max>=3 (continuum VP nonzero).
    For n_max=2, returns None.
    """
    if n_max == 2:
        return None
    if n_max == 3:
        return float(PROJECTION_CONSTANT_NMAX3)
    # For larger n_max we could recompute, but it's expensive.
    # Conservatively return None and report only for n_max=3.
    return None


def power_law_fit(xs: List[int], ys: List[float]) -> Tuple[float, float, float]:
    """Fit y = A * x^alpha (log-log) via least squares.

    Returns (A, alpha, R^2).
    """
    if len(xs) < 2:
        return None, None, None
    log_x = np.log(np.array(xs, dtype=float))
    log_y = np.log(np.array(ys, dtype=float))
    coeffs, residuals, *_ = np.polyfit(log_x, log_y, 1, full=True)
    alpha = float(coeffs[0])
    log_A = float(coeffs[1])
    A = np.exp(log_A)
    # R^2
    log_y_pred = log_A + alpha * log_x
    ss_res = np.sum((log_y - log_y_pred) ** 2)
    ss_tot = np.sum((log_y - np.mean(log_y)) ** 2)
    R2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-30 else 1.0
    return A, alpha, R2


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("VP-1: Vector-photon promotion of graph-native QED vertex correction")
    print("=" * 70)
    print(f"alpha = {ALPHA:.10f}")
    print(f"alpha/(2*pi) = {ALPHA_OVER_2PI:.10e}")
    print()

    results = {
        'sprint': 'VP-1',
        'description': (
            'Test whether replacing the scalar (Fock-edge) photon in the '
            'graph-native QED vertex correction with a vector photon '
            '((q, m_q), Wigner 3j vertex) closes the C × F2 negative result.'
        ),
        'alpha': ALPHA,
        'alpha_over_2pi': ALPHA_OVER_2PI,
        'projection_constant_nmax3': str(PROJECTION_CONSTANT_NMAX3),
        'projection_constant_nmax3_float': float(PROJECTION_CONSTANT_NMAX3),
        'numpy_results': {},
        'exact_results': {},
        'scaling_analysis': {},
        'verdict': None,
    }

    # ---- numpy: n_max = 2, 3, 4 ----
    for n_max in [2, 3, 4]:
        print(f"--- n_max = {n_max} (numpy) ---")
        t0 = time.time()
        res = compute_vector_F2_numpy(n_max=n_max)
        elapsed = time.time() - t0
        print(f"  N_e = {res['N_electron']}, N_gamma = {res['N_photon']}, "
              f"q_max = {res['q_max']}")
        print(f"  Vertex sparsity = {res['vertex_sparsity']:.1%} "
              f"({res['vertex_nonzero']}/{res['vertex_total']} nonzero)")
        print(f"  Tr(Lambda) = {res['Lambda_trace']:.6e}")
        print(f"  ||Lambda||_F = {res['Lambda_frobenius']:.6e}")
        print(f"  Lambda(GS,GS) max = {res['Lambda_GS_block_max']:.6e}  "
              f"(should be 0 if vector parity holds)")
        print(f"  V_bare diag sum = {res['V_bare_diag_sum']:.6e}  "
              f"(zero by Furry/spinor phase)")
        print(f"  Tr(G_e) = {res['Tr_G_e']:.6e}")
        print(f"  F2_vector (Tr(Lambda)/Tr(G_e))   = {res['F2_vector_vs_TrGe']:.6e}")
        print(f"  F2_vector (Tr(Lambda)/||V||_F^2) = {res['F2_vector_vs_VFrob']:.6e}")
        print(f"  F2_vector (Tr(Lambda)/N_e)       = {res['F2_vector_vs_Ne']:.6e}")
        print(f"  elapsed: {elapsed:.2f}s")
        print()
        results['numpy_results'][f'n_max_{n_max}'] = res

    # ---- exact: n_max = 2 (and try 3) ----
    for n_max in [2, 3]:
        print(f"--- n_max = {n_max} (exact sympy) ---")
        t0 = time.time()
        try:
            res_exact = compute_vector_F2_exact(n_max=n_max)
            elapsed = time.time() - t0
            print(f"  Tr(Lambda) (exact) = {res_exact['Lambda_trace_exact']}")
            print(f"  Tr(G_e) (exact)    = {res_exact['Tr_G_e_exact']}")
            print(f"  F2_vector (exact)  = {res_exact['F2_vector_exact']}")
            print(f"  F2_vector (float)  = {res_exact['F2_vector_float']:.10e}")
            print(f"  elapsed: {elapsed:.2f}s")
            print()
            results['exact_results'][f'n_max_{n_max}'] = res_exact
        except Exception as e:
            elapsed = time.time() - t0
            print(f"  EXCEPTION ({type(e).__name__}): {e}")
            print(f"  elapsed before exception: {elapsed:.2f}s")
            results['exact_results'][f'n_max_{n_max}'] = {'error': str(e)}

    # ---- comparison: F2_vector vs F2_scalar, and C × F2 ----
    print("=" * 70)
    print("COMPARISON: F2_vector vs F2_scalar; C × F2_vector vs alpha/(2*pi)")
    print("=" * 70)

    # Use the Tr(Lambda)/Tr(G_e) normalization (most natural for matching
    # scalar graph-native extraction Tr(Lambda)/Tr(V_bare*G_e), since
    # V_bare diag = 0 forces a fallback)
    print(f"{'n_max':<7} {'F2_scalar':<12} {'F2_vec/TrGe':<14} {'F2_vec/||V||^2':<16} "
          f"{'C':<14} {'C*F2_vec':<14} {'ratio_vs_a/2pi':<16}")
    print("-" * 110)

    comparison_rows = []
    for n_max in [2, 3, 4]:
        res = results['numpy_results'][f'n_max_{n_max}']
        F2_sc = F2_scalar_reference(n_max)
        F2_vec_TrGe = res['F2_vector_vs_TrGe']
        F2_vec_VFrob = res['F2_vector_vs_VFrob']
        C = projection_constant_C(n_max)
        if C is None:
            C_str = 'N/A'
            CF2_str = 'N/A'
            ratio_str = 'N/A'
            CF2_TrGe = None
            ratio_TrGe = None
        else:
            CF2_TrGe = C * F2_vec_TrGe
            ratio_TrGe = CF2_TrGe / ALPHA_OVER_2PI
            C_str = f"{C:.6e}"
            CF2_str = f"{CF2_TrGe:.6e}"
            ratio_str = f"{ratio_TrGe:.4f}"

        F2_sc_str = f"{F2_sc:.6f}" if F2_sc is not None else "N/A"
        print(f"{n_max:<7} {F2_sc_str:<12} {F2_vec_TrGe:<14.6e} "
              f"{F2_vec_VFrob:<16.6e} {C_str:<14} {CF2_str:<14} {ratio_str:<16}")
        comparison_rows.append({
            'n_max': n_max,
            'F2_scalar_reference': F2_sc,
            'F2_vector_TrGe': F2_vec_TrGe,
            'F2_vector_VFrob': F2_vec_VFrob,
            'F2_vector_avg_diag': res['F2_vector_vs_Ne'],
            'projection_constant_C': C,
            'C_times_F2_vector_TrGe': CF2_TrGe,
            'ratio_to_alpha_over_2pi': ratio_TrGe,
        })

    results['comparison'] = comparison_rows
    print()

    # ---- scaling analysis ----
    print("=" * 70)
    print("SCALING ANALYSIS: F2_vector(n_max) power law")
    print("=" * 70)

    # Use the TrGe normalization
    xs = [2, 3, 4]
    ys_TrGe = [results['numpy_results'][f'n_max_{n}']['F2_vector_vs_TrGe'] for n in xs]
    ys_VFrob = [results['numpy_results'][f'n_max_{n}']['F2_vector_vs_VFrob'] for n in xs]

    print(f"F2_vector (TrGe normalization): {ys_TrGe}")
    print(f"F2_vector (VFrob normalization): {ys_VFrob}")

    # Only fit power-law if all values are positive (else log breaks)
    if all(y > 0 for y in ys_TrGe):
        A_TrGe, alpha_TrGe, R2_TrGe = power_law_fit(xs, ys_TrGe)
        print(f"  TrGe fit:  F2 = {A_TrGe:.4f} * n^{alpha_TrGe:.4f}, R² = {R2_TrGe:.6f}")
        results['scaling_analysis']['F2_vector_TrGe'] = {
            'A': A_TrGe, 'alpha': alpha_TrGe, 'R_squared': R2_TrGe,
        }
    else:
        print(f"  TrGe fit: skipped (some F2 values not positive)")
        results['scaling_analysis']['F2_vector_TrGe'] = {
            'A': None, 'alpha': None, 'R_squared': None,
            'note': 'Not all positive; cannot fit log-log',
        }

    if all(y > 0 for y in ys_VFrob):
        A_VFrob, alpha_VFrob, R2_VFrob = power_law_fit(xs, ys_VFrob)
        print(f"  VFrob fit: F2 = {A_VFrob:.4f} * n^{alpha_VFrob:.4f}, R² = {R2_VFrob:.6f}")
        results['scaling_analysis']['F2_vector_VFrob'] = {
            'A': A_VFrob, 'alpha': alpha_VFrob, 'R_squared': R2_VFrob,
        }
    else:
        print(f"  VFrob fit: skipped (some F2 values not positive)")
        results['scaling_analysis']['F2_vector_VFrob'] = {
            'A': None, 'alpha': None, 'R_squared': None,
            'note': 'Not all positive; cannot fit log-log',
        }

    # F2_scalar reference scaling: -0.573 (CLAUDE.md, n_max=2..6)
    results['scaling_analysis']['F2_scalar_reference_alpha'] = -0.573
    print(f"  F2_scalar reference (CLAUDE.md): F2 ~ n^(-0.573)")

    # ---- verdict ----
    print()
    print("=" * 70)
    print("VERDICT")
    print("=" * 70)

    # Hypothesis: C × F2_vector should converge to alpha/(2*pi).
    # We have only one data point at n_max=3 where C is known.
    if comparison_rows[1]['ratio_to_alpha_over_2pi'] is not None:
        ratio_n3 = comparison_rows[1]['ratio_to_alpha_over_2pi']
        if 0.9 <= ratio_n3 <= 1.1:
            verdict = 'POSITIVE'
            verdict_text = (
                f'C × F2_vector at n_max=3 = {comparison_rows[1]["C_times_F2_vector_TrGe"]:.4e} '
                f'matches alpha/(2*pi) = {ALPHA_OVER_2PI:.4e} within 10% '
                f'(ratio = {ratio_n3:.4f}). Hypothesis CONFIRMED at this truncation.'
            )
        else:
            verdict = 'NEGATIVE'
            verdict_text = (
                f'C × F2_vector at n_max=3 = {comparison_rows[1]["C_times_F2_vector_TrGe"]:.4e} '
                f'is far from alpha/(2*pi) = {ALPHA_OVER_2PI:.4e} (ratio = {ratio_n3:.4f}). '
                f'Vector promotion does NOT close the C × F2 negative result. '
                f'Structural reason: the vector vertex enforces selection rules '
                f'(GS structural zero, m_j conservation, etc.) but does NOT '
                f'self-consistently match the projection constant C, which was '
                f'derived from the SCALAR graph-VP vs continuum-VP ratio. '
                f'The continuum F2 = alpha/(2*pi) requires both the vector '
                f'vertex AND the alpha-dependent coupling at the action level, '
                f'which is the calibration tier (Paper 18) that lives outside '
                f'the graph.'
            )
    else:
        verdict = 'INCONCLUSIVE'
        verdict_text = (
            'C is not defined at n_max=2 (continuum VP zero), and we have '
            'only single-point data at n_max=3. Need n_max=4+ projection '
            'constant for trend analysis.'
        )

    print(verdict)
    print(verdict_text)
    results['verdict'] = verdict
    results['verdict_text'] = verdict_text

    # Save JSON
    output_path = (Path(__file__).parent / 'data'
                   / 'vp1_vector_graph_native.json')
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nSaved: {output_path}")

    return results


if __name__ == "__main__":
    main()
