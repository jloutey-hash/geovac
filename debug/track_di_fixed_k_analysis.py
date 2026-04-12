"""
Track DI Sprint 3B: Fixed-k=Z=2 Sequence Analysis.

Runs the algebraic CI at fixed k=Z=2 through n_max=7, analyzes
convergence, and compares with the variational-k sequence from Sprint 3.

Key question: does the fixed-k sequence show Schwartz p~4 cusp
convergence, unlike the variational-k sequence (p~1.5, 1.4% floor)?
"""

import json
import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.casimir_ci import build_fci_matrix


E_EXACT = -2.903724377

# Sprint 3 variational-k results for comparison
VAR_K_RESULTS = [
    {'n_max': 1, 'k_var': 1.687500, 'E_var': -2.84765625, 'error_pct': 1.9309},
    {'n_max': 2, 'k_var': 1.762531, 'E_var': -2.85706703, 'error_pct': 1.6068},
    {'n_max': 3, 'k_var': 1.784824, 'E_var': -2.86012180, 'error_pct': 1.5016},
    {'n_max': 4, 'k_var': 1.793911, 'E_var': -2.86140461, 'error_pct': 1.4574},
    {'n_max': 5, 'k_var': 1.798519, 'E_var': -2.86206800, 'error_pct': 1.4346},
    {'n_max': 6, 'k_var': 1.801202, 'E_var': -2.86246068, 'error_pct': 1.4211},
    {'n_max': 7, 'k_var': 1.802917, 'E_var': -2.86271605, 'error_pct': 1.4123},
]


def compute_fixed_k_sequence(n_max_max: int = 7, k: float = 2.0) -> list:
    """Compute FCI energies at fixed k for n_max=1..n_max_max."""
    print(f"Computing fixed k={k} sequence n_max=1..{n_max_max}")
    print("=" * 60)

    results = []
    for n_max in range(1, n_max_max + 1):
        t0 = time.time()
        H = build_fci_matrix(Z=2, n_max=n_max, k_orb=k)
        t_build = time.time() - t0

        eigenvalues = np.linalg.eigvalsh(H)
        E0 = eigenvalues[0]

        error_ha = E0 - E_EXACT
        error_pct = abs(error_ha / E_EXACT) * 100

        entry = {
            'n_max': n_max,
            'n_configs': H.shape[0],
            'k': k,
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
    errors = [r['error_ha'] for r in results]

    deltas = [energies[i] - energies[i + 1] for i in range(len(energies) - 1)]
    ratios = []
    for i in range(len(deltas) - 1):
        if abs(deltas[i + 1]) > 1e-15:
            ratios.append(deltas[i] / deltas[i + 1])
        else:
            ratios.append(None)

    second_diffs = [deltas[i] - deltas[i + 1] for i in range(len(deltas) - 1)]

    return {
        'energies': energies,
        'errors': errors,
        'deltas': deltas,
        'ratios': ratios,
        'second_diffs': second_diffs,
    }


def power_law_fit_2param(results: list) -> dict:
    """Fit error(n) = C/n^p assuming E_inf = E_exact (log-log regression)."""
    data = [(r['n_max'], r['error_ha']) for r in results if r['n_max'] >= 2]
    if len(data) < 2:
        return {'status': 'insufficient_data'}

    ns = np.array([d[0] for d in data], dtype=float)
    errors = np.array([d[1] for d in data])

    log_n = np.log(ns)
    log_err = np.log(errors)
    coeffs = np.polyfit(log_n, log_err, 1)
    p_fit = -coeffs[0]
    C_fit = np.exp(coeffs[1])

    fitted = C_fit / ns ** p_fit
    rms = np.sqrt(np.mean((errors - fitted) ** 2))

    return {'status': 'converged', 'C': float(C_fit), 'p': float(p_fit),
            'rms_residual': float(rms)}


def power_law_fit_3param(results: list) -> dict:
    """Fit error(n) = delta_inf + C/n^p (3-parameter)."""
    from scipy.optimize import curve_fit

    data = [(r['n_max'], r['error_ha']) for r in results if r['n_max'] >= 2]
    if len(data) < 3:
        return {'status': 'insufficient_data'}

    ns = np.array([d[0] for d in data], dtype=float)
    errors = np.array([d[1] for d in data])

    def model(n, delta_inf, C, p):
        return delta_inf + C / n ** p

    p0 = [0.0, errors[0] * ns[0] ** 4, 4.0]
    try:
        popt, _ = curve_fit(model, ns, errors, p0=p0, maxfev=10000)
        delta_inf, C_fit, p_fit = popt
        E_inf = E_EXACT + delta_inf
        fitted = model(ns, *popt)
        rms = np.sqrt(np.mean((errors - fitted) ** 2))

        return {
            'status': 'converged',
            'E_inf': float(E_inf),
            'E_inf_error_pct': abs(delta_inf / E_EXACT) * 100,
            'delta_inf': float(delta_inf),
            'C': float(C_fit),
            'p': float(p_fit),
            'rms_residual': float(rms),
        }
    except Exception as e:
        return {'status': f'failed: {e}'}


def richardson_extrapolation(results: list, p: float) -> list:
    """Richardson extrapolation with assumed power p."""
    rich = []
    for i in range(1, len(results)):
        n = results[i]['n_max']
        n_prev = results[i - 1]['n_max']
        E_n = results[i]['E0']
        E_prev = results[i - 1]['E0']

        np_cur = n ** p
        np_prev = n_prev ** p
        denom = np_cur - np_prev
        if abs(denom) < 1e-30:
            continue

        E_rich = (np_cur * E_n - np_prev * E_prev) / denom
        error = E_rich - E_EXACT
        error_pct = abs(error / E_EXACT) * 100

        rich.append({
            'n': n, 'E_richardson': float(E_rich),
            'error_ha': float(error), 'error_pct': float(error_pct),
        })
    return rich


def aitken_acceleration(results: list) -> list:
    """Aitken delta-squared acceleration."""
    aitken = []
    for i in range(len(results) - 2):
        E0 = results[i]['E0']
        E1 = results[i + 1]['E0']
        E2 = results[i + 2]['E0']
        denom = E2 - 2 * E1 + E0
        if abs(denom) < 1e-30:
            continue
        E_ait = E0 - (E1 - E0) ** 2 / denom
        error = E_ait - E_EXACT
        aitken.append({
            'n_max': results[i]['n_max'],
            'E_aitken': float(E_ait),
            'error_ha': float(error),
            'error_pct': abs(error / E_EXACT) * 100,
        })
    return aitken


def generate_plots(fixed_results: list, diff: dict, rich_p4: list,
                   aitken: list) -> None:
    """Generate comparison plots."""
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

    # --- Plot 1: Error comparison (log scale) ---
    ax = axes[0]
    ns_var = [r['n_max'] for r in VAR_K_RESULTS]
    err_var = [abs((r['E_var'] - E_EXACT) / E_EXACT) * 100 for r in VAR_K_RESULTS]
    ns_fix = [r['n_max'] for r in fixed_results]
    err_fix = [r['error_pct'] for r in fixed_results]

    ax.semilogy(ns_var, err_var, 'bo-', markersize=7, label='Variational k')
    ax.semilogy(ns_fix, err_fix, 'rs-', markersize=7, label='Fixed k=Z=2')
    if rich_p4:
        ax.semilogy([r['n'] for r in rich_p4],
                     [r['error_pct'] for r in rich_p4],
                     'g^--', markersize=5, label='Richardson (p=4)')
    if aitken:
        ax.semilogy([a['n_max'] for a in aitken],
                     [a['error_pct'] for a in aitken],
                     'mv--', markersize=5, label='Aitken')
    ax.set_xlabel('n_max')
    ax.set_ylabel('Error (%)')
    ax.set_title('He CI: Variational vs Fixed k')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # --- Plot 2: Log-log error vs n (fixed k only, n>=2) ---
    ax = axes[1]
    ns_arr = np.array([r['n_max'] for r in fixed_results if r['n_max'] >= 2], dtype=float)
    err_arr = np.array([r['error_ha'] for r in fixed_results if r['n_max'] >= 2])
    ax.loglog(ns_arr, err_arr, 'rs-', markersize=8)
    if len(ns_arr) >= 2:
        log_n = np.log(ns_arr)
        log_err = np.log(err_arr)
        coeffs = np.polyfit(log_n, log_err, 1)
        fit_n = np.linspace(ns_arr[0] * 0.8, ns_arr[-1] * 1.3, 50)
        ax.loglog(fit_n, np.exp(np.polyval(coeffs, np.log(fit_n))),
                  'r--', alpha=0.5, label=f'slope={coeffs[0]:.2f} (p={-coeffs[0]:.2f})')
    ax.set_xlabel('n_max')
    ax.set_ylabel('Error (Ha)')
    ax.set_title('Fixed k=2: Power-Law (n>=2)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # --- Plot 3: Convergence ratios ---
    ax = axes[2]
    if diff['ratios']:
        ns_rat = [fixed_results[i + 1]['n_max'] for i in range(len(diff['ratios']))]
        rats = [r for r in diff['ratios'] if r is not None]
        ns_rat = ns_rat[:len(rats)]
        ax.plot(ns_rat, rats, 'rs-', markersize=7, label='Fixed k=2')
    # Variational-k ratios for comparison
    var_deltas = [VAR_K_RESULTS[i]['E_var'] - VAR_K_RESULTS[i + 1]['E_var']
                  for i in range(len(VAR_K_RESULTS) - 1)]
    var_ratios = [var_deltas[i] / var_deltas[i + 1]
                  for i in range(len(var_deltas) - 1)
                  if abs(var_deltas[i + 1]) > 1e-15]
    ns_vrat = [VAR_K_RESULTS[i + 1]['n_max'] for i in range(len(var_ratios))]
    ax.plot(ns_vrat, var_ratios, 'bo-', markersize=7, label='Variational k')
    ax.axhline(y=1.0, color='gray', linestyle=':', alpha=0.5)
    ax.set_xlabel('n_max')
    ax.set_ylabel('delta_n / delta_{n+1}')
    ax.set_title('Convergence Ratios')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(plot_dir, 'track_di_fixed_k.png')
    plt.savefig(path, dpi=150)
    print(f"[PLOT] Saved: {path}")
    plt.close()


def main():
    t_start = time.time()

    # Step 1: Compute fixed-k sequence
    results = compute_fixed_k_sequence(n_max_max=7, k=2.0)

    # Step 2: Verify against FCI-A
    print("\n--- FCI-A Verification ---")
    fci_a_he_nmax5 = -2.8936  # from paper_fci_atoms.tex
    our_nmax5 = results[4]['E0']
    diff_fcia = abs(our_nmax5 - fci_a_he_nmax5)
    print(f"  FCI-A He n_max=5: {fci_a_he_nmax5:.4f} Ha (0.35% error)")
    print(f"  Our k=2 n_max=5:  {our_nmax5:.8f} Ha ({results[4]['error_pct']:.4f}% error)")
    print(f"  Difference: {diff_fcia:.6f} Ha")
    fcia_match = diff_fcia < 0.01
    print(f"  Match: {'YES' if fcia_match else 'NO'} (threshold 0.01 Ha)")

    # Step 3: Difference analysis
    diff = difference_analysis(results)
    print("\n--- Fixed k=2 Energy Sequence ---")
    print(f"{'n_max':>5} {'E0 (Ha)':>14} {'Error (Ha)':>12} {'Error (%)':>10} "
          f"{'delta':>12} {'ratio':>8}")
    for i, r in enumerate(results):
        d_str = f"{diff['deltas'][i]:.6f}" if i < len(diff['deltas']) else ""
        rat_str = ""
        if i > 0 and i - 1 < len(diff['ratios']):
            rat = diff['ratios'][i - 1]
            rat_str = f"{rat:.4f}" if rat is not None else "---"
        print(f"{r['n_max']:>5} {r['E0']:>14.8f} {r['error_ha']:>12.6f} "
              f"{r['error_pct']:>10.4f} {d_str:>12} {rat_str:>8}")

    # Step 4: Power-law fits
    print("\n--- Power-Law Fit: error(n) = C/n^p (E_inf = E_exact) ---")
    fit_2p = power_law_fit_2param(results)
    if fit_2p['status'] == 'converged':
        print(f"  C = {fit_2p['C']:.6f}, p = {fit_2p['p']:.4f}, "
              f"RMS = {fit_2p['rms_residual']:.2e} Ha")

    print("\n--- Power-Law Fit: error(n) = delta_inf + C/n^p (3-param) ---")
    fit_3p = power_law_fit_3param(results)
    if fit_3p['status'] == 'converged':
        print(f"  E_inf = {fit_3p['E_inf']:.8f} Ha "
              f"(delta from exact = {fit_3p['delta_inf']:.2e})")
        print(f"  E_inf error = {fit_3p['E_inf_error_pct']:.4f}%")
        print(f"  C = {fit_3p['C']:.6f}, p = {fit_3p['p']:.4f}")
        print(f"  RMS = {fit_3p['rms_residual']:.2e} Ha")

    # Step 5: Richardson
    print("\n--- Richardson Extrapolation ---")
    for p_test in [2.0, 3.0, 4.0, 5.0]:
        rich = richardson_extrapolation(results, p_test)
        if rich:
            last = rich[-1]
            print(f"  p={p_test:.1f}: E_rich(n={last['n']}) = "
                  f"{last['E_richardson']:.8f}, error = {last['error_pct']:.4f}%")

    rich_p4 = richardson_extrapolation(results, 4.0)
    if rich_p4:
        print(f"\n  Richardson (p=4) detail:")
        for r in rich_p4:
            print(f"    n={r['n']}: E = {r['E_richardson']:.8f}, "
                  f"error = {r['error_ha']:.6f} Ha ({r['error_pct']:.4f}%)")

    if fit_3p['status'] == 'converged':
        p_fitted = fit_3p['p']
        rich_fitted = richardson_extrapolation(results, p_fitted)
        if rich_fitted:
            last = rich_fitted[-1]
            print(f"\n  Richardson (p={p_fitted:.2f}, fitted): "
                  f"E(n={last['n']}) = {last['E_richardson']:.8f}, "
                  f"error = {last['error_pct']:.4f}%")

    # Aitken
    print("\n--- Aitken Delta^2 ---")
    aitken = aitken_acceleration(results)
    for a in aitken:
        print(f"  n_max={a['n_max']}: E = {a['E_aitken']:.8f}, "
              f"error = {a['error_ha']:.6f} Ha ({a['error_pct']:.4f}%)")

    # Step 6: Comparison table
    print("\n--- Comparison: Variational-k vs Fixed k=2 ---")
    print(f"{'n_max':>5} | {'E(k_var)':>12} {'err(k_var)':>10} | "
          f"{'E(k=2)':>12} {'err(k=2)':>10} | {'k_var':>8} | {'winner':>8}")
    print("-" * 82)
    for i, (vr, fr) in enumerate(zip(VAR_K_RESULTS, results)):
        var_err = vr['error_pct']
        fix_err = fr['error_pct']
        winner = "var-k" if var_err < fix_err else "k=2"
        print(f"{vr['n_max']:>5} | {vr['E_var']:>12.8f} {var_err:>9.4f}% | "
              f"{fr['E0']:>12.8f} {fix_err:>9.4f}% | {vr['k_var']:>8.4f} | {winner:>8}")

    # Plots
    generate_plots(results, diff, rich_p4, aitken)

    # Save JSON
    output = {
        'metadata': {
            'track': 'DI Sprint 3B',
            'description': 'Fixed k=Z=2 algebraic CI convergence',
            'E_exact_He': E_EXACT,
            'k_fixed': 2.0,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'total_time_s': time.time() - t_start,
        },
        'fixed_k_sequence': [{
            'n_max': r['n_max'], 'n_configs': r['n_configs'],
            'k': r['k'], 'E0': r['E0'],
            'error_ha': r['error_ha'], 'error_pct': r['error_pct'],
            'build_time_s': r['build_time_s'],
        } for r in results],
        'difference_analysis': {
            'deltas': diff['deltas'],
            'ratios': [r if r is not None else None for r in diff['ratios']],
            'second_diffs': diff['second_diffs'],
        },
        'power_law_2param': fit_2p,
        'power_law_3param': fit_3p,
        'richardson_p4': rich_p4,
        'richardson_p_fitted': (
            richardson_extrapolation(results, fit_3p['p'])
            if fit_3p['status'] == 'converged' else []
        ),
        'aitken': aitken,
        'fcia_verification': {
            'fcia_nmax5': fci_a_he_nmax5,
            'our_nmax5': float(our_nmax5),
            'difference': float(diff_fcia),
            'match': fcia_match,
        },
        'variational_k_comparison': VAR_K_RESULTS,
    }

    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
    os.makedirs(data_dir, exist_ok=True)
    out_path = os.path.join(data_dir, 'track_di_fixed_k.json')
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n[DATA] Saved: {out_path}")

    total_time = time.time() - t_start
    print(f"\nSprint 3B complete. Total time: {total_time:.1f}s")


if __name__ == '__main__':
    main()
