"""
f2_convergence_nmax56.py

Compute the graph-native anomalous magnetic moment F₂(κ) at n_max=2..6
using numpy for efficiency at large n_max. Also includes:
  - Convergence characterization (power-law, exponential, Richardson)
  - Pendant-edge theorem verification
  - Self-energy spectral decomposition
  - F₂(0) vs F₂(κ) comparison

Scientific goal: characterize F₂(n_max) convergence to determine whether
the graph-native anomalous moment converges to zero, a finite constant,
or diverges as n_max → ∞.

Output: debug/data/f2_convergence_nmax56.json
"""

import json
import sys
import time
from pathlib import Path

# Force UTF-8 stdout on Windows
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
if hasattr(sys.stderr, 'reconfigure'):
    sys.stderr.reconfigure(encoding='utf-8')

import numpy as np
from scipy import optimize

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels, kappa_to_l
from geovac.graph_qed_propagator import DiracGraphOperator, KAPPA_SCALAR
from geovac.graph_qed_photon import build_fock_graph
from geovac.graph_qed_vertex import (
    build_projection_matrix, build_vertex_tensor, vertex_tensor_to_matrices
)

ALPHA = 1.0 / 137.035999084
SCHWINGER = ALPHA / (2 * np.pi)
KAPPA_FLOAT = -1.0 / 16.0


# ---------------------------------------------------------------------------
# Core numerical F₂ computation (numpy only, no sympy eigensolve)
# ---------------------------------------------------------------------------

def compute_f2_numpy(n_max: int, t_val: float = KAPPA_FLOAT):
    """Compute F₂ and self-energy diagnostics purely in numpy.

    Returns a dict with F₂, self-energy trace, eigenvalues, pendant check, etc.
    """
    import sympy as sp
    from sympy import Rational

    t0 = time.time()
    print(f"  n_max={n_max}: building graph objects...", flush=True)

    # 1. Build Fock graph topology (edge list, incidence B, L0, L1)
    fock_data = build_fock_graph(n_max)
    E_fock = fock_data.E
    V_fock = fock_data.V

    # 2. Photon propagator G_gamma = pinv(L1) via numpy
    L1_np = np.array(fock_data.L1.tolist(), dtype=np.float64)
    G_gamma_np = np.linalg.pinv(L1_np)
    print(f"    V_fock={V_fock}, E_fock={E_fock}, "
          f"L1 shape={L1_np.shape}", flush=True)

    # 3. Vertex tensor V[a, b, e] -> list of V_e matrices
    #    build_vertex_tensor uses sympy CG coefficients -> convert to float
    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    entries, N_dirac, V_f, E_f = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    assert E_f == E_fock

    # Convert vertex matrices to numpy
    V_mats_sym = vertex_tensor_to_matrices(entries, N_dirac, E_fock)
    V_mats_np = [np.array(vm.tolist(), dtype=np.float64) for vm in V_mats_sym]

    # Bare vertex: V_bare = sum_e V_e
    V_bare_np = sum(V_mats_np)
    print(f"    N_dirac={N_dirac}, vertex built, "
          f"elapsed={time.time()-t0:.1f}s", flush=True)

    # 4. Electron propagator G_e = D^{-1} via numpy
    t_sym = Rational(-1, 16) if abs(t_val - KAPPA_FLOAT) < 1e-15 else Rational(0)
    if abs(t_val) > 1e-15 and abs(t_val - KAPPA_FLOAT) > 1e-15:
        # General t -- fall back to float-to-rational
        t_sym = sp.nsimplify(t_val, rational=True)
    op = DiracGraphOperator(n_max=n_max, t=t_sym)
    D_np = op.matrix_numpy()
    G_e_np = np.linalg.inv(D_np)
    print(f"    D shape={D_np.shape}, cond={np.linalg.cond(D_np):.2e}, "
          f"elapsed={time.time()-t0:.1f}s", flush=True)

    # 5. Self-energy: Sigma = sum_{e1,e2} G_gamma[e1,e2] * V_{e1} . V_{e2}^T
    Sigma_np = np.zeros((N_dirac, N_dirac))
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_ee = G_gamma_np[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            Sigma_np += g_ee * (V_mats_np[e1] @ V_mats_np[e2].T)
    print(f"    Sigma computed, elapsed={time.time()-t0:.1f}s", flush=True)

    # 6. Vertex correction: Lambda = sum_{e1,e2} G_gamma[e1,e2] * V_{e1} . G_e . V_{e2}^T
    Lambda_np = np.zeros((N_dirac, N_dirac))
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_ee = G_gamma_np[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            Lambda_np += g_ee * (V_mats_np[e1] @ G_e_np @ V_mats_np[e2].T)
    print(f"    Lambda computed, elapsed={time.time()-t0:.1f}s", flush=True)

    # 7. Extract F₂
    tr_lambda = np.trace(Lambda_np)
    tr_norm = np.trace(V_bare_np @ G_e_np)
    F2 = tr_lambda / tr_norm if abs(tr_norm) > 1e-15 else 0.0

    # 8. Self-energy diagnostics
    Sigma_eigs = np.linalg.eigvalsh(Sigma_np)
    Sigma_trace = np.trace(Sigma_np)
    D_trace = np.trace(D_np)

    # Ground state block check (pendant-edge theorem)
    gs_indices = [i for i, lab in enumerate(dirac_labels)
                  if lab.n_fock == 1 and lab.kappa == -1]
    if gs_indices:
        gs_block = Sigma_np[np.ix_(gs_indices, gs_indices)]
        gs_sigma_00 = gs_block[0, 0]
        pendant_predicted = 2.0 * (n_max - 1) / n_max
    else:
        gs_sigma_00 = None
        pendant_predicted = None

    # Number of zero eigenvalues in Sigma
    n_zero_sigma = np.sum(np.abs(Sigma_eigs) < 1e-10)
    n_nonzero_sigma = N_dirac - n_zero_sigma

    elapsed = time.time() - t0
    if gs_sigma_00 is not None:
        print(f"    F2={F2:.8f}, Tr(Sigma)={Sigma_trace:.6f}, "
              f"pendant_pred={pendant_predicted:.8f}, gs_Sigma={gs_sigma_00:.8f}", flush=True)
    else:
        print(f"    F2={F2:.8f}, Tr(Sigma)={Sigma_trace:.6f}", flush=True)
    print(f"    Done. elapsed={elapsed:.1f}s", flush=True)

    return {
        'n_max': n_max,
        't': t_val,
        'N_dirac': N_dirac,
        'E_fock': E_fock,
        'V_fock': V_fock,
        'F2': F2,
        'Tr_Lambda': tr_lambda,
        'Tr_V_bare_G_e': tr_norm,
        'Sigma_trace': Sigma_trace,
        'Sigma_eigenvalues': sorted(Sigma_eigs.tolist()),
        'Sigma_n_zero': int(n_zero_sigma),
        'Sigma_n_nonzero': int(n_nonzero_sigma),
        'Sigma_max_eigval': float(np.max(Sigma_eigs)),
        'Sigma_min_nonzero_eigval': float(np.min(np.abs(Sigma_eigs[np.abs(Sigma_eigs) > 1e-10]))) if n_nonzero_sigma > 0 else None,
        'D_trace': float(D_trace),
        'Tr_Sigma_over_Tr_D': float(Sigma_trace / D_trace) if abs(D_trace) > 1e-15 else None,
        'gs_sigma_00': float(gs_sigma_00) if gs_sigma_00 is not None else None,
        'pendant_predicted': pendant_predicted,
        'pendant_match': bool(abs(gs_sigma_00 - pendant_predicted) < 1e-8) if gs_sigma_00 is not None else None,
        'D_condition': float(np.linalg.cond(D_np)),
        'elapsed_seconds': elapsed,
    }


def compute_f2_at_t0_numpy(n_max: int):
    """Compute F₂ at t=0 (diagonal propagator) using numpy."""
    return compute_f2_numpy(n_max, t_val=0.0)


# ---------------------------------------------------------------------------
# Convergence analysis
# ---------------------------------------------------------------------------

def fit_power_law(n_max_vals, f2_vals):
    """Fit log(F2) = a + b*log(n_max) -> F2 ~ C * n_max^b."""
    log_n = np.log(n_max_vals)
    log_f2 = np.log(f2_vals)
    b, a = np.polyfit(log_n, log_f2, 1)
    C = np.exp(a)
    residuals = log_f2 - (a + b * log_n)
    r_squared = 1 - np.sum(residuals**2) / np.sum((log_f2 - np.mean(log_f2))**2)
    return {'exponent': b, 'prefactor': C, 'R_squared': r_squared}


def fit_exponential(n_max_vals, f2_vals):
    """Fit log(F2) = a + b*n_max -> F2 ~ C * exp(b*n_max)."""
    log_f2 = np.log(f2_vals)
    b, a = np.polyfit(n_max_vals, log_f2, 1)
    C = np.exp(a)
    residuals = log_f2 - (a + b * n_max_vals)
    r_squared = 1 - np.sum(residuals**2) / np.sum((log_f2 - np.mean(log_f2))**2)
    return {'rate': b, 'prefactor': C, 'R_squared': r_squared}


def fit_inverse_power(n_max_vals, f2_vals, power=1):
    """Fit F2 = a + b / n_max^power."""
    x = 1.0 / n_max_vals**power
    b, a = np.polyfit(x, f2_vals, 1)
    residuals = f2_vals - (a + b * x)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((f2_vals - np.mean(f2_vals))**2)
    r_squared = 1 - ss_res / ss_tot
    return {'limit': a, 'coefficient': b, 'power': power, 'R_squared': r_squared}


def richardson_extrapolation(n_vals, f_vals, order=1):
    """Richardson extrapolation assuming F(n) ~ F_inf + c/n^p.

    Uses successive pairs to estimate F_inf and p.
    """
    results = []
    n = np.array(n_vals, dtype=float)
    f = np.array(f_vals, dtype=float)

    # Simple Richardson: for three consecutive values
    # F(n1), F(n2), F(n3) with h_i = 1/n_i
    for i in range(len(n) - 2):
        n1, n2, n3 = n[i], n[i+1], n[i+2]
        f1, f2, f3 = f[i], f[i+1], f[i+2]

        # Estimate p from three-point formula:
        # (f1 - f2) / (f2 - f3) = (h1^p - h2^p) / (h2^p - h3^p)
        # where h_i = 1/n_i
        h1, h2, h3 = 1/n1, 1/n2, 1/n3
        ratio = (f1 - f2) / (f2 - f3) if abs(f2 - f3) > 1e-15 else None

        # Estimate p by solving ratio = (h1^p - h2^p)/(h2^p - h3^p)
        # This is transcendental; use numerical root-finding
        if ratio is not None and ratio > 0:
            def obj(p):
                return (h1**p - h2**p) / (h2**p - h3**p) - ratio
            try:
                p_est = optimize.brentq(obj, 0.01, 10.0)
                # With p estimated, get F_inf
                F_inf = (f2 * h1**p_est - f1 * h2**p_est) / (h1**p_est - h2**p_est)
                results.append({
                    'points': [int(n1), int(n2), int(n3)],
                    'p_estimated': p_est,
                    'F_inf': F_inf,
                })
            except (ValueError, ZeroDivisionError):
                results.append({
                    'points': [int(n1), int(n2), int(n3)],
                    'p_estimated': None,
                    'F_inf': None,
                    'note': 'root-finding failed'
                })
        else:
            results.append({
                'points': [int(n1), int(n2), int(n3)],
                'p_estimated': None,
                'F_inf': None,
                'note': f'invalid ratio={ratio}'
            })

    return results


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def main():
    results = {
        'description': (
            'Graph-native QED F2 convergence study at n_max=2..6. '
            'All computations use numpy for the photon propagator pseudoinverse '
            'and electron propagator inverse. Vertex tensor from CG coefficients '
            'is converted to float.'
        ),
        'schwinger_value': SCHWINGER,
    }

    print("=" * 70, flush=True)
    print("F2 Convergence Sprint: n_max = 2..6", flush=True)
    print("=" * 70, flush=True)
    print(f"Schwinger alpha/(2*pi) = {SCHWINGER:.8e}", flush=True)
    print(flush=True)

    # -----------------------------------------------------------------------
    # Task 1: Compute F₂ at n_max = 2..6
    # -----------------------------------------------------------------------
    f2_table = []
    f2_t0_table = []

    for n_max in range(2, 7):
        print(f"\n--- n_max = {n_max} (t=kappa=-1/16) ---", flush=True)
        row = compute_f2_numpy(n_max, t_val=KAPPA_FLOAT)
        f2_table.append(row)

        print(f"\n--- n_max = {n_max} (t=0) ---", flush=True)
        row_t0 = compute_f2_numpy(n_max, t_val=0.0)
        f2_t0_table.append(row_t0)

    results['f2_at_kappa'] = f2_table
    results['f2_at_t0'] = f2_t0_table

    # Print summary table
    print("\n" + "=" * 70, flush=True)
    print("F2 Convergence Table", flush=True)
    print("-" * 70, flush=True)
    print(f"{'n_max':>5} {'N_dirac':>8} {'E_fock':>7} {'F2(kappa)':>14} "
          f"{'F2(0)':>14} {'ratio':>10} {'pendant?':>10}", flush=True)
    print("-" * 70, flush=True)
    for kappa_row, t0_row in zip(f2_table, f2_t0_table):
        n = kappa_row['n_max']
        f2k = kappa_row['F2']
        f2_0 = t0_row['F2']
        ratio = f2k / f2_0 if abs(f2_0) > 1e-15 else float('nan')
        pendant = "YES" if kappa_row.get('pendant_match') else "NO"
        print(f"{n:>5} {kappa_row['N_dirac']:>8} {kappa_row['E_fock']:>7} "
              f"{f2k:>14.8f} {f2_0:>14.8f} {ratio:>10.6f} {pendant:>10}", flush=True)

    # -----------------------------------------------------------------------
    # Task 2: Convergence characterization
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70, flush=True)
    print("Convergence Analysis", flush=True)
    print("-" * 70, flush=True)

    n_vals = np.array([row['n_max'] for row in f2_table], dtype=float)
    f2_vals = np.array([row['F2'] for row in f2_table])

    # Power-law fit: F2 ~ C * n^b
    pl = fit_power_law(n_vals, f2_vals)
    print(f"Power-law: F2 ~ {pl['prefactor']:.4f} * n^({pl['exponent']:.4f}), "
          f"R^2 = {pl['R_squared']:.6f}", flush=True)

    # Exponential fit: F2 ~ C * exp(b*n)
    ef = fit_exponential(n_vals, f2_vals)
    print(f"Exponential: F2 ~ {ef['prefactor']:.4f} * exp({ef['rate']:.4f}*n), "
          f"R^2 = {ef['R_squared']:.6f}", flush=True)

    # 1/n fit: F2 = a + b/n
    inv1 = fit_inverse_power(n_vals, f2_vals, power=1)
    print(f"1/n fit: F2 -> {inv1['limit']:.6f} + {inv1['coefficient']:.4f}/n, "
          f"R^2 = {inv1['R_squared']:.6f}", flush=True)

    # 1/n^2 fit: F2 = a + b/n^2
    inv2 = fit_inverse_power(n_vals, f2_vals, power=2)
    print(f"1/n^2 fit: F2 -> {inv2['limit']:.6f} + {inv2['coefficient']:.4f}/n^2, "
          f"R^2 = {inv2['R_squared']:.6f}", flush=True)

    # 1/log(n) fit: F2 = a + b/log(n)
    inv_log_x = 1.0 / np.log(n_vals)
    b_log, a_log = np.polyfit(inv_log_x, f2_vals, 1)
    res_log = f2_vals - (a_log + b_log * inv_log_x)
    r2_log = 1 - np.sum(res_log**2) / np.sum((f2_vals - np.mean(f2_vals))**2)
    print(f"1/log(n) fit: F2 -> {a_log:.6f} + {b_log:.4f}/log(n), "
          f"R^2 = {r2_log:.6f}", flush=True)

    # Local power-law exponents (consecutive pairs)
    local_exponents = []
    for i in range(len(n_vals) - 1):
        b_local = (np.log(f2_vals[i+1]) - np.log(f2_vals[i])) / \
                  (np.log(n_vals[i+1]) - np.log(n_vals[i]))
        local_exponents.append({
            'n_pair': [int(n_vals[i]), int(n_vals[i+1])],
            'exponent': b_local,
        })
        print(f"  Local exponent ({int(n_vals[i])},{int(n_vals[i+1])}): {b_local:.4f}", flush=True)

    results['convergence'] = {
        'power_law': pl,
        'exponential': ef,
        'inverse_n': inv1,
        'inverse_n2': inv2,
        'inverse_log_n': {'limit': a_log, 'coefficient': b_log, 'R_squared': r2_log},
        'local_exponents': local_exponents,
    }

    # -----------------------------------------------------------------------
    # Richardson extrapolation
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70, flush=True)
    print("Richardson Extrapolation", flush=True)
    print("-" * 70, flush=True)

    rich = richardson_extrapolation(n_vals.tolist(), f2_vals.tolist())
    for r in rich:
        pts = r['points']
        p = r.get('p_estimated')
        finf = r.get('F_inf')
        note = r.get('note', '')
        print(f"  Points {pts}: p={p:.4f}, F_inf={finf:.8f} {note}" if p is not None
              else f"  Points {pts}: {note}", flush=True)

    results['richardson'] = rich

    # -----------------------------------------------------------------------
    # Task 4: Pendant-edge and self-energy gap analysis
    # -----------------------------------------------------------------------
    print("\n" + "=" * 70, flush=True)
    print("Self-Energy Diagnostics", flush=True)
    print("-" * 70, flush=True)

    se_diagnostics = []
    for row in f2_table:
        n = row['n_max']
        pred = row['pendant_predicted']
        actual = row['gs_sigma_00']
        match = row['pendant_match']
        n_zero = row['Sigma_n_zero']
        n_nonzero = row['Sigma_n_nonzero']
        tr_ratio = row['Tr_Sigma_over_Tr_D']

        se_diagnostics.append({
            'n_max': n,
            'pendant_predicted': pred,
            'gs_sigma_00': actual,
            'pendant_match': match,
            'Sigma_n_zero': n_zero,
            'Sigma_n_nonzero': n_nonzero,
            'Tr_Sigma_over_Tr_D': tr_ratio,
        })

        print(f"  n_max={n}: Sigma(GS)={actual:.8f}, "
              f"pred={pred:.8f}, match={match}, "
              f"zero_eigs={n_zero}, nonzero={n_nonzero}, "
              f"Tr(Sigma)/Tr(D)={tr_ratio:.6f}" if tr_ratio is not None else "", flush=True)

    results['self_energy_diagnostics'] = se_diagnostics

    # -----------------------------------------------------------------------
    # Save results
    # -----------------------------------------------------------------------
    out_path = PROJECT_ROOT / 'debug' / 'data' / 'f2_convergence_nmax56.json'
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Clean up numpy arrays for JSON serialization
    def clean_for_json(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        return obj

    with open(out_path, 'w') as f:
        json.dump(results, f, indent=2, default=clean_for_json)
    print(f"\nResults saved to {out_path}", flush=True)

    return results


if __name__ == '__main__':
    main()
