"""
Track DI Sprint 3C: Graph-Native CI Convergence Analysis.

Uses hybrid h1 (exact diagonal + graph Laplacian off-diagonal)
with analytical Slater V_ee at k=Z. Zero free parameters.

Compares three sequences:
  1. Variational-k (Sprint 3): diagonal h1, optimized k
  2. Fixed k=Z=2 (Sprint 3B): diagonal h1, k=2
  3. Graph-native (this sprint): hybrid h1, k=Z for V_ee, no optimization
"""

import json
import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.casimir_ci import build_graph_native_fci


E_EXACT = -2.903724377

# Sprint 3 variational-k results
VAR_K = [
    {'n_max': 1, 'E': -2.84765625, 'error_pct': 1.9309},
    {'n_max': 2, 'E': -2.85706703, 'error_pct': 1.6068},
    {'n_max': 3, 'E': -2.86012180, 'error_pct': 1.5016},
    {'n_max': 4, 'E': -2.86140461, 'error_pct': 1.4574},
    {'n_max': 5, 'E': -2.86206800, 'error_pct': 1.4346},
    {'n_max': 6, 'E': -2.86246068, 'error_pct': 1.4211},
    {'n_max': 7, 'E': -2.86271605, 'error_pct': 1.4123},
]

# Sprint 3B fixed k=2 results
FIXED_K = [
    {'n_max': 1, 'E': -2.75000000, 'error_pct': 5.2940},
    {'n_max': 2, 'E': -2.83388556, 'error_pct': 2.4051},
    {'n_max': 3, 'E': -2.84352847, 'error_pct': 2.0731},
    {'n_max': 4, 'E': -2.84690325, 'error_pct': 1.9568},
    {'n_max': 5, 'E': -2.84852248, 'error_pct': 1.9011},
    {'n_max': 6, 'E': -2.84944028, 'error_pct': 1.8695},
    {'n_max': 7, 'E': -2.85001882, 'error_pct': 1.8495},
]


def compute_graph_native_sequence(n_max_max: int = 7) -> list:
    """Compute graph-native FCI for n_max=1..n_max_max."""
    print(f"Computing graph-native CI sequence n_max=1..{n_max_max}")
    print("=" * 60)

    results = []
    for n_max in range(1, n_max_max + 1):
        t0 = time.time()
        H = build_graph_native_fci(Z=2, n_max=n_max)
        t_build = time.time() - t0

        E0 = np.linalg.eigvalsh(H)[0]
        error_ha = E0 - E_EXACT
        error_pct = abs(error_ha / E_EXACT) * 100

        entry = {
            'n_max': n_max,
            'n_configs': H.shape[0],
            'E0': float(E0),
            'error_ha': float(error_ha),
            'error_pct': float(error_pct),
            'build_time_s': t_build,
        }
        results.append(entry)
        print(f"  n_max={n_max}: {H.shape[0]:>5} configs, "
              f"E={E0:.8f} Ha, error={error_pct:.4f}% ({error_ha:.6f} Ha) "
              f"[{t_build:.1f}s]")

    return results


def difference_analysis(results: list) -> dict:
    """Successive differences and ratios."""
    energies = [r['E0'] for r in results]
    deltas = [energies[i] - energies[i + 1] for i in range(len(energies) - 1)]
    ratios = []
    for i in range(len(deltas) - 1):
        ratios.append(deltas[i] / deltas[i + 1] if abs(deltas[i + 1]) > 1e-15 else None)
    return {'deltas': deltas, 'ratios': ratios}


def power_law_fit_3param(results: list) -> dict:
    """Fit error(n) = delta_inf + C/n^p."""
    from scipy.optimize import curve_fit
    data = [(r['n_max'], r['error_ha']) for r in results if r['n_max'] >= 2]
    if len(data) < 3:
        return {'status': 'insufficient_data'}
    ns = np.array([d[0] for d in data], dtype=float)
    errors = np.array([d[1] for d in data])

    def model(n, delta_inf, C, p):
        return delta_inf + C / n ** p

    try:
        popt, _ = curve_fit(model, ns, errors, p0=[0.0, errors[0] * ns[0] ** 4, 4.0],
                            maxfev=10000)
        delta_inf, C_fit, p_fit = popt
        fitted = model(ns, *popt)
        rms = np.sqrt(np.mean((errors - fitted) ** 2))
        return {
            'status': 'converged', 'E_inf': float(E_EXACT + delta_inf),
            'delta_inf': float(delta_inf),
            'E_inf_error_pct': abs(delta_inf / E_EXACT) * 100,
            'C': float(C_fit), 'p': float(p_fit), 'rms_residual': float(rms),
        }
    except Exception as e:
        return {'status': f'failed: {e}'}


def power_law_fit_2param(results: list) -> dict:
    """Fit error(n) = C/n^p assuming E_inf = E_exact."""
    data = [(r['n_max'], r['error_ha']) for r in results if r['n_max'] >= 2]
    if len(data) < 2:
        return {'status': 'insufficient_data'}
    ns = np.array([d[0] for d in data], dtype=float)
    errors = np.array([d[1] for d in data])
    coeffs = np.polyfit(np.log(ns), np.log(errors), 1)
    return {
        'status': 'converged', 'C': float(np.exp(coeffs[1])),
        'p': float(-coeffs[0]),
        'rms_residual': float(np.sqrt(np.mean(
            (errors - np.exp(coeffs[1]) / ns ** (-coeffs[0])) ** 2))),
    }


def richardson_extrapolation(results: list, p: float) -> list:
    """Richardson extrapolation."""
    rich = []
    for i in range(1, len(results)):
        n, n_prev = results[i]['n_max'], results[i - 1]['n_max']
        E_n, E_prev = results[i]['E0'], results[i - 1]['E0']
        denom = n ** p - n_prev ** p
        if abs(denom) < 1e-30:
            continue
        E_rich = (n ** p * E_n - n_prev ** p * E_prev) / denom
        error = E_rich - E_EXACT
        rich.append({'n': n, 'E_richardson': float(E_rich),
                     'error_ha': float(error),
                     'error_pct': abs(error / E_EXACT) * 100})
    return rich


def aitken_acceleration(results: list) -> list:
    """Aitken delta-squared."""
    aitken = []
    for i in range(len(results) - 2):
        E0, E1, E2 = results[i]['E0'], results[i + 1]['E0'], results[i + 2]['E0']
        denom = E2 - 2 * E1 + E0
        if abs(denom) < 1e-30:
            continue
        E_ait = E0 - (E1 - E0) ** 2 / denom
        error = E_ait - E_EXACT
        aitken.append({'n_max': results[i]['n_max'], 'E_aitken': float(E_ait),
                       'error_ha': float(error),
                       'error_pct': abs(error / E_EXACT) * 100})
    return aitken


def generate_plots(graph_results: list, diff: dict, rich_p4: list,
                   aitken: list) -> None:
    """Generate three-sequence comparison plots."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("[PLOT] matplotlib not available")
        return

    plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    ns_g = [r['n_max'] for r in graph_results]
    err_g = [r['error_pct'] for r in graph_results]
    ns_v = [r['n_max'] for r in VAR_K]
    err_v = [r['error_pct'] for r in VAR_K]
    ns_f = [r['n_max'] for r in FIXED_K]
    err_f = [r['error_pct'] for r in FIXED_K]

    # Plot 1: Three-sequence error comparison
    ax = axes[0]
    ax.semilogy(ns_v, err_v, 'bo-', markersize=6, label='Variational k')
    ax.semilogy(ns_f, err_f, 'rs-', markersize=6, label='Fixed k=Z=2')
    ax.semilogy(ns_g, err_g, 'g^-', markersize=7, label='Graph-native (hybrid)')
    if rich_p4:
        ax.semilogy([r['n'] for r in rich_p4], [r['error_pct'] for r in rich_p4],
                     'gv--', markersize=5, alpha=0.7, label='Richardson (p=4)')
    ax.set_xlabel('n_max')
    ax.set_ylabel('Error (%)')
    ax.set_title('Three-Sequence Comparison')
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # Plot 2: Log-log convergence (graph-native, n>=2)
    ax = axes[1]
    ns_arr = np.array([r['n_max'] for r in graph_results if r['n_max'] >= 2], dtype=float)
    err_arr = np.array([r['error_ha'] for r in graph_results if r['n_max'] >= 2])
    ax.loglog(ns_arr, err_arr, 'g^-', markersize=8)
    if len(ns_arr) >= 2:
        coeffs = np.polyfit(np.log(ns_arr), np.log(err_arr), 1)
        fit_n = np.linspace(ns_arr[0] * 0.8, ns_arr[-1] * 1.3, 50)
        ax.loglog(fit_n, np.exp(np.polyval(coeffs, np.log(fit_n))),
                  'g--', alpha=0.5, label=f'slope={coeffs[0]:.2f} (p={-coeffs[0]:.2f})')
    ax.set_xlabel('n_max')
    ax.set_ylabel('Error (Ha)')
    ax.set_title('Graph-Native: Power-Law')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 3: Convergence ratios
    ax = axes[2]
    if diff['ratios']:
        rats = [r for r in diff['ratios'] if r is not None]
        ns_rat = [graph_results[i + 1]['n_max'] for i in range(len(rats))]
        ax.plot(ns_rat, rats, 'g^-', markersize=7, label='Graph-native')
    # Variational-k ratios
    vd = [VAR_K[i]['E'] - VAR_K[i + 1]['E'] for i in range(len(VAR_K) - 1)]
    vr = [vd[i] / vd[i + 1] for i in range(len(vd) - 1) if abs(vd[i + 1]) > 1e-15]
    ax.plot([VAR_K[i + 1]['n_max'] for i in range(len(vr))], vr,
            'bo-', markersize=6, label='Variational k')
    ax.axhline(y=1.0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('n_max')
    ax.set_ylabel('delta_n / delta_{n+1}')
    ax.set_title('Convergence Ratios')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(plot_dir, 'track_di_graph_native.png')
    plt.savefig(path, dpi=150)
    print(f"[PLOT] Saved: {path}")
    plt.close()


def main():
    t_start = time.time()

    results = compute_graph_native_sequence(n_max_max=7)

    # Difference analysis
    diff = difference_analysis(results)
    print("\n--- Graph-Native Energy Sequence ---")
    print(f"{'n_max':>5} {'E0 (Ha)':>14} {'Error (Ha)':>12} {'Error (%)':>10} "
          f"{'delta':>12} {'ratio':>8}")
    for i, r in enumerate(results):
        d = f"{diff['deltas'][i]:.6f}" if i < len(diff['deltas']) else ""
        rat = ""
        if i > 0 and i - 1 < len(diff['ratios']):
            rat = f"{diff['ratios'][i-1]:.4f}" if diff['ratios'][i-1] else "---"
        print(f"{r['n_max']:>5} {r['E0']:>14.8f} {r['error_ha']:>12.6f} "
              f"{r['error_pct']:>10.4f} {d:>12} {rat:>8}")

    # Power-law fits
    print("\n--- Power-Law Fit: error(n) = C/n^p (E_inf = E_exact) ---")
    fit_2p = power_law_fit_2param(results)
    if fit_2p['status'] == 'converged':
        print(f"  C = {fit_2p['C']:.6f}, p = {fit_2p['p']:.4f}")

    print("\n--- Power-Law Fit: error(n) = delta_inf + C/n^p (3-param) ---")
    fit_3p = power_law_fit_3param(results)
    if fit_3p['status'] == 'converged':
        print(f"  E_inf = {fit_3p['E_inf']:.8f} (delta = {fit_3p['delta_inf']:.2e})")
        print(f"  E_inf error = {fit_3p['E_inf_error_pct']:.4f}%")
        print(f"  C = {fit_3p['C']:.6f}, p = {fit_3p['p']:.4f}")
        print(f"  RMS = {fit_3p['rms_residual']:.2e} Ha")

    # Richardson
    print("\n--- Richardson Extrapolation ---")
    for p in [2.0, 3.0, 4.0, 5.0]:
        rich = richardson_extrapolation(results, p)
        if rich:
            last = rich[-1]
            print(f"  p={p:.1f}: E(n={last['n']}) = {last['E_richardson']:.8f}, "
                  f"error = {last['error_pct']:.4f}%")

    rich_p4 = richardson_extrapolation(results, 4.0)
    if rich_p4:
        print(f"\n  Richardson (p=4) detail:")
        for r in rich_p4:
            print(f"    n={r['n']}: E = {r['E_richardson']:.8f}, "
                  f"error = {r['error_ha']:.6f} Ha ({r['error_pct']:.4f}%)")

    if fit_3p['status'] == 'converged':
        rich_fit = richardson_extrapolation(results, fit_3p['p'])
        if rich_fit:
            last = rich_fit[-1]
            print(f"\n  Richardson (p={fit_3p['p']:.2f}): "
                  f"E(n={last['n']}) = {last['E_richardson']:.8f}, "
                  f"error = {last['error_pct']:.4f}%")

    # Aitken
    print("\n--- Aitken Delta^2 ---")
    aitken = aitken_acceleration(results)
    for a in aitken:
        print(f"  n_max={a['n_max']}: E = {a['E_aitken']:.8f}, "
              f"error = {a['error_ha']:.6f} Ha ({a['error_pct']:.4f}%)")

    # Three-sequence comparison
    print("\n--- Three-Sequence Comparison ---")
    print(f"{'n_max':>5} | {'Var-k':>10} {'err%':>7} | {'k=2':>10} {'err%':>7} "
          f"| {'Graph':>10} {'err%':>7} | {'Best':>6}")
    print("-" * 78)
    for i in range(min(len(VAR_K), len(FIXED_K), len(results))):
        vk = VAR_K[i]
        fk = FIXED_K[i]
        gr = results[i]
        errs = {'var-k': vk['error_pct'], 'k=2': fk['error_pct'],
                'graph': gr['error_pct']}
        best = min(errs, key=errs.get)
        print(f"{vk['n_max']:>5} | {vk['E']:>10.6f} {vk['error_pct']:>6.2f}% "
              f"| {fk['E']:>10.6f} {fk['error_pct']:>6.2f}% "
              f"| {gr['E0']:>10.6f} {gr['error_pct']:>6.2f}% | {best:>6}")

    # Plots
    generate_plots(results, diff, rich_p4, aitken)

    # Save JSON
    output = {
        'metadata': {
            'track': 'DI Sprint 3C',
            'description': 'Graph-native CI: hybrid h1 + Slater V_ee, zero free parameters',
            'E_exact_He': E_EXACT,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'total_time_s': time.time() - t_start,
        },
        'graph_native_sequence': [{
            'n_max': r['n_max'], 'n_configs': r['n_configs'],
            'E0': r['E0'], 'error_ha': r['error_ha'],
            'error_pct': r['error_pct'], 'build_time_s': r['build_time_s'],
        } for r in results],
        'difference_analysis': {
            'deltas': diff['deltas'],
            'ratios': [r if r is not None else None for r in diff['ratios']],
        },
        'power_law_2param': fit_2p,
        'power_law_3param': fit_3p,
        'richardson_p4': rich_p4,
        'richardson_p_fitted': (
            richardson_extrapolation(results, fit_3p['p'])
            if fit_3p['status'] == 'converged' else []
        ),
        'aitken': aitken,
        'comparison': {
            'variational_k': VAR_K,
            'fixed_k2': FIXED_K,
        },
    }

    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
    os.makedirs(data_dir, exist_ok=True)
    out_path = os.path.join(data_dir, 'track_di_graph_native.json')
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n[DATA] Saved: {out_path}")
    print(f"\nSprint 3C complete. Total time: {time.time() - t_start:.1f}s")


if __name__ == '__main__':
    main()
