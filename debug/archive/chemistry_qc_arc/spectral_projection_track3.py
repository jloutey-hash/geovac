"""
Track 3: Cross-diagram projection constants (numpy-only path).

Computes C_VP, C_SE, C_F2 at n_max=2..5 and characterizes their
structural independence. Uses numpy for L1 pseudoinverse (bypasses
sympy eigenvalue bug at n_max>=4).

Output: debug/data/spectral_projection_constants.json
"""

import json
import sys
import time
from pathlib import Path

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

import numpy as np

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from geovac.dirac_matrix_elements import DiracLabel, iter_dirac_labels
from geovac.graph_qed_propagator import DiracGraphOperator
from geovac.graph_qed_photon import build_fock_graph
from geovac.graph_qed_vertex import (
    build_projection_matrix, build_vertex_tensor, vertex_tensor_to_matrices
)
from geovac.qed_self_energy import self_energy_spectral
from geovac.graph_qed_continuum_bridge import continuum_vp_truncated

from sympy import Rational

ALPHA = 1.0 / 137.035999177
SCHWINGER = ALPHA / (2 * np.pi)
KAPPA_FLOAT = -1.0 / 16.0


def compute_graph_qed_numpy(n_max: int, t_val: float = KAPPA_FLOAT):
    """Compute all graph-native QED quantities at given n_max using numpy."""
    t0 = time.time()

    fock_data = build_fock_graph(n_max)
    E_fock = fock_data.E
    V_fock = fock_data.V

    L1_np = np.array(fock_data.L1.tolist(), dtype=np.float64)
    G_gamma_np = np.linalg.pinv(L1_np)

    P, dirac_labels, fock_states = build_projection_matrix(n_max)
    entries, N_dirac, V_f, E_f = build_vertex_tensor(
        n_max, P=P, dirac_labels=dirac_labels, fock_data=fock_data
    )
    V_mats_sym = vertex_tensor_to_matrices(entries, N_dirac, E_fock)
    V_mats_np = [np.array(vm.tolist(), dtype=np.float64) for vm in V_mats_sym]
    V_bare_np = sum(V_mats_np)

    op = DiracGraphOperator(n_max=n_max, t=Rational(-1, 16))
    D_np = op.matrix_numpy()
    G_e_np = np.linalg.inv(D_np)

    # Self-energy
    Sigma_np = np.zeros((N_dirac, N_dirac))
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_ee = G_gamma_np[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            Sigma_np += g_ee * (V_mats_np[e1] @ V_mats_np[e2].T)

    # Vertex correction
    Lambda_np = np.zeros((N_dirac, N_dirac))
    for e1 in range(E_fock):
        for e2 in range(E_fock):
            g_ee = G_gamma_np[e1, e2]
            if abs(g_ee) < 1e-15:
                continue
            Lambda_np += g_ee * (V_mats_np[e1] @ G_e_np @ V_mats_np[e2].T)

    # VP = Tr(V_bare^T . G_e . V_bare . G_e)
    VP_matrix = V_bare_np.T @ G_e_np @ V_bare_np @ G_e_np
    graph_vp_trace = float(np.trace(VP_matrix))

    # F2
    Tr_Lambda = float(np.trace(Lambda_np))
    Tr_V_bare_G_e = float(np.trace(V_bare_np @ G_e_np))
    F2 = Tr_Lambda / Tr_V_bare_G_e if abs(Tr_V_bare_G_e) > 1e-15 else 0.0

    elapsed = time.time() - t0
    return {
        'n_max': n_max,
        'N_dirac': N_dirac,
        'E_fock': E_fock,
        'graph_se_trace': float(np.trace(Sigma_np)),
        'graph_vp_trace': graph_vp_trace,
        'graph_F2': F2,
        'Tr_Lambda': Tr_Lambda,
        'Tr_V_bare_G_e': Tr_V_bare_G_e,
        'elapsed': elapsed,
    }


def compute_continuum_se_trace(n_max_ch: int) -> float:
    """Continuum self-energy trace: sum g(n_ext)*Sigma(n_ext)."""
    total = 0.0
    for n_ext in range(n_max_ch + 1):
        g_ext = 2 * (n_ext + 1) * (n_ext + 2)
        sigma_n = float(self_energy_spectral(n_ext, n_max_ch))
        total += g_ext * sigma_n
    return total


def compute_continuum_vp(n_max_ch: int) -> float:
    """Continuum VP trace (W-weighted spectral sum)."""
    try:
        return float(continuum_vp_truncated(n_max_ch))
    except Exception:
        from geovac.graph_qed_continuum_bridge import _continuum_vp_full_estimate
        return _continuum_vp_full_estimate(n_max_sum=n_max_ch)


def main():
    output_path = Path(__file__).parent / "data" / "spectral_projection_constants.json"
    output_path.parent.mkdir(parents=True, exist_ok=True)

    results = {
        'description': 'Cross-diagram projection constants C_VP, C_SE, C_F2 at n_max=2..5',
        'alpha': ALPHA,
        'schwinger': SCHWINGER,
    }

    all_data = []

    for n_max in [2, 3, 4, 5]:
        print(f"\n{'='*60}")
        print(f"  n_max = {n_max}")
        print(f"{'='*60}")

        n_max_ch = n_max - 1

        # Graph quantities (numpy)
        print(f"  Computing graph QED...", flush=True)
        gd = compute_graph_qed_numpy(n_max)
        print(f"    graph_SE_Tr={gd['graph_se_trace']:.4f}, "
              f"graph_VP_Tr={gd['graph_vp_trace']:.4f}, "
              f"graph_F2={gd['graph_F2']:.6f}  ({gd['elapsed']:.1f}s)")

        # Continuum VP
        print(f"  Computing continuum VP...", flush=True)
        t0 = time.time()
        cont_vp = compute_continuum_vp(n_max_ch)
        print(f"    cont_VP={cont_vp:.6f}  ({time.time()-t0:.1f}s)")

        # Continuum SE trace
        print(f"  Computing continuum SE...", flush=True)
        t0 = time.time()
        cont_se = compute_continuum_se_trace(n_max_ch)
        print(f"    cont_SE_Tr={cont_se:.6f}  ({time.time()-t0:.1f}s)")

        # Projection constants
        C_VP = cont_vp / gd['graph_vp_trace'] if gd['graph_vp_trace'] != 0 else 0
        C_SE = cont_se / gd['graph_se_trace'] if gd['graph_se_trace'] != 0 else 0
        C_F2 = SCHWINGER / gd['graph_F2'] if gd['graph_F2'] != 0 else 0

        print(f"  => C_VP={C_VP:.6f}, C_SE={C_SE:.6f}, C_F2={C_F2:.6e}")

        entry = {
            'n_max': n_max,
            **gd,
            'cont_vp': cont_vp,
            'cont_se_trace': cont_se,
            'C_VP': C_VP,
            'C_SE': C_SE,
            'C_F2': C_F2,
        }
        all_data.append(entry)

    results['data'] = all_data

    # Cross-diagram analysis
    print(f"\n{'='*60}")
    print("  CROSS-DIAGRAM SUMMARY")
    print(f"{'='*60}")
    print(f"  {'n':>3} {'C_VP':>10} {'C_SE':>10} {'C_F2':>12} {'VP/SE':>10} {'SE/F2':>12}")
    print(f"  {'---':>3} {'----------':>10} {'----------':>10} {'------------':>12} {'----------':>10} {'------------':>12}")
    cross = []
    for d in all_data:
        vp_se = d['C_VP'] / d['C_SE'] if d['C_SE'] != 0 else None
        se_f2 = d['C_SE'] / d['C_F2'] if d['C_F2'] != 0 else None
        r = {'n_max': d['n_max'], 'VP/SE': vp_se, 'SE/F2': se_f2}
        cross.append(r)
        vp_se_str = f"{vp_se:.6f}" if vp_se is not None else "N/A"
        se_f2_str = f"{se_f2:.4f}" if se_f2 is not None else "N/A"
        print(f"  {d['n_max']:3d} {d['C_VP']:10.6f} {d['C_SE']:10.6f} {d['C_F2']:12.6e} "
              f"{vp_se_str:>10} {se_f2_str:>12}")
    results['cross'] = cross

    # Power law fits
    print(f"\n  Power law fits (C ~ A * n^beta):")
    for name, key in [('C_VP', 'C_VP'), ('C_SE', 'C_SE'), ('C_F2', 'C_F2')]:
        ns = np.array([d['n_max'] for d in all_data if d[key] > 0])
        cs = np.array([d[key] for d in all_data if d[key] > 0])
        if len(ns) >= 2:
            log_n = np.log(ns)
            log_c = np.log(cs)
            beta, log_A = np.polyfit(log_n, log_c, 1)
            A = np.exp(log_A)
            resid = log_c - (log_A + beta * log_n)
            ss_res = np.sum(resid**2)
            ss_tot = np.sum((log_c - np.mean(log_c))**2)
            R2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
            print(f"    {name}: {A:.4e} * n^{beta:.3f}  (R2={R2:.4f})")
            results[f'{name}_fit'] = {'A': float(A), 'beta': float(beta), 'R2': float(R2)}
        else:
            print(f"    {name}: insufficient data")

    # Graph trace scaling
    print(f"\n  Graph trace scaling:")
    for name, key in [('Tr(Sigma)', 'graph_se_trace'), ('Tr(Pi)', 'graph_vp_trace'), ('F2', 'graph_F2')]:
        ns = np.array([d['n_max'] for d in all_data])
        vs = np.array([d[key] for d in all_data])
        mask = vs > 0
        if np.sum(mask) >= 2:
            log_n = np.log(ns[mask])
            log_v = np.log(vs[mask])
            beta, log_A = np.polyfit(log_n, log_v, 1)
            A = np.exp(log_A)
            print(f"    graph {name} ~ {A:.4f} * n^{beta:.3f}")

    # Key diagnostic: ratio of graph traces
    print(f"\n  Graph internal structure:")
    print(f"  {'n':>3} {'Tr(Pi)/Tr(Sig)':>14} {'Tr(Sig)/F2':>12} {'Tr(Pi)/F2':>10}")
    for d in all_data:
        r1 = d['graph_vp_trace'] / d['graph_se_trace'] if d['graph_se_trace'] != 0 else 0
        r2 = d['graph_se_trace'] / d['graph_F2'] if d['graph_F2'] != 0 else 0
        r3 = d['graph_vp_trace'] / d['graph_F2'] if d['graph_F2'] != 0 else 0
        print(f"  {d['n_max']:3d} {r1:14.6f} {r2:12.4f} {r3:10.4f}")

    with output_path.open('w') as f:
        json.dump(results, f, indent=2, default=lambda x: float(x) if isinstance(x, np.floating) else x)
    print(f"\n  Saved: {output_path}")


if __name__ == "__main__":
    main()
