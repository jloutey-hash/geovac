"""
Track DI Sprint 3: Algebraic CI Convergence Sequence Analysis.

Extends the Casimir CI to n_max=4-7, then analyzes the convergence
structure of the variational energy sequence {E*(n_max)}.

Outputs:
  - debug/data/track_di_sequence.json  (full numerical results)
  - Convergence plots (if matplotlib available)
"""

import json
import sys
import os
import time
import numpy as np

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.casimir_ci import extended_convergence_table


E_EXACT = -2.903724377  # He ground state (Ha)


def compute_sequence(n_max_max: int = 7) -> list:
    """Run extended_convergence_table for n_max=1..n_max_max."""
    print("=" * 70)
    print(f"Phase A: Computing algebraic CI sequence n_max=1..{n_max_max}")
    print("=" * 70)
    return extended_convergence_table(Z=2, n_max_range=range(1, n_max_max + 1))


def difference_analysis(results: list) -> dict:
    """Compute successive differences, ratios, and second differences."""
    energies = [r['E_var'] for r in results]
    errors = [r['E_var'] - E_EXACT for r in results]  # positive (above exact)
    n_values = [r['n_max'] for r in results]

    deltas = []  # energy gain per shell
    for i in range(len(energies) - 1):
        deltas.append(energies[i] - energies[i + 1])  # positive if improving

    ratios = []
    for i in range(len(deltas) - 1):
        if abs(deltas[i + 1]) > 1e-15:
            ratios.append(deltas[i] / deltas[i + 1])
        else:
            ratios.append(None)

    second_diffs = []
    for i in range(len(deltas) - 1):
        second_diffs.append(deltas[i] - deltas[i + 1])

    return {
        'n_values': n_values,
        'energies': energies,
        'errors': errors,
        'deltas': deltas,
        'ratios': ratios,
        'second_diffs': second_diffs,
    }


def power_law_fit(results: list) -> dict:
    """Fit E*(n) = E_inf + C / n^p.

    Uses least-squares on log(error) vs log(n) for data points n>=2
    (n=1 is qualitatively different — no angular correlation).
    """
    # Use n_max >= 2 data points
    data = [(r['n_max'], r['E_var'] - E_EXACT) for r in results if r['n_max'] >= 2]
    if len(data) < 3:
        return {'status': 'insufficient_data', 'n_points': len(data)}

    ns = np.array([d[0] for d in data], dtype=float)
    errors = np.array([d[1] for d in data])

    # Three-parameter fit: E*(n) = E_inf + C/n^p
    # Rearrange: error(n) = E*(n) - E_exact = (E_inf - E_exact) + C/n^p
    # Let delta_inf = E_inf - E_exact, then error(n) = delta_inf + C/n^p
    #
    # Use scipy curve_fit for the 3-parameter fit
    from scipy.optimize import curve_fit

    def model(n, delta_inf, C, p):
        return delta_inf + C / n ** p

    # Initial guess: delta_inf ≈ 0, C ≈ error[0]*ns[0]^4, p ≈ 4
    p0 = [0.0, errors[0] * ns[0] ** 4, 4.0]
    try:
        popt, pcov = curve_fit(model, ns, errors, p0=p0, maxfev=10000)
        delta_inf, C_fit, p_fit = popt
        E_inf = E_EXACT + delta_inf

        # Residuals
        fitted = model(ns, *popt)
        residuals = errors - fitted
        rms_residual = np.sqrt(np.mean(residuals ** 2))

        return {
            'status': 'converged',
            'E_inf': float(E_inf),
            'E_inf_error_pct': abs(delta_inf / E_EXACT) * 100,
            'C': float(C_fit),
            'p': float(p_fit),
            'rms_residual': float(rms_residual),
            'n_points': len(data),
        }
    except Exception as e:
        return {'status': f'failed: {e}', 'n_points': len(data)}


def two_param_power_law_fit(results: list) -> dict:
    """Fit error(n) = C / n^p assuming E_inf = E_exact.

    Uses log-log linear regression: log(error) = log(C) - p*log(n).
    """
    data = [(r['n_max'], r['E_var'] - E_EXACT) for r in results if r['n_max'] >= 2]
    if len(data) < 2:
        return {'status': 'insufficient_data'}

    ns = np.array([d[0] for d in data], dtype=float)
    errors = np.array([d[1] for d in data])

    log_n = np.log(ns)
    log_err = np.log(errors)

    # Linear fit: log(err) = log(C) - p*log(n)
    coeffs = np.polyfit(log_n, log_err, 1)
    p_fit = -coeffs[0]
    C_fit = np.exp(coeffs[1])

    fitted = C_fit / ns ** p_fit
    residuals = errors - fitted
    rms_residual = np.sqrt(np.mean(residuals ** 2))

    return {
        'status': 'converged',
        'C': float(C_fit),
        'p': float(p_fit),
        'rms_residual': float(rms_residual),
    }


def richardson_extrapolation(results: list, p: float) -> list:
    """Apply Richardson extrapolation with assumed power p.

    E_Rich(n) = [n^p * E*(n) - (n-1)^p * E*(n-1)] / [n^p - (n-1)^p]
    """
    rich = []
    for i in range(1, len(results)):
        n = results[i]['n_max']
        n_prev = results[i - 1]['n_max']
        E_n = results[i]['E_var']
        E_prev = results[i - 1]['E_var']

        np_cur = n ** p
        np_prev = n_prev ** p
        denom = np_cur - np_prev
        if abs(denom) < 1e-30:
            continue

        E_rich = (np_cur * E_n - np_prev * E_prev) / denom
        error = E_rich - E_EXACT
        error_pct = abs(error / E_EXACT) * 100

        rich.append({
            'n': n,
            'n_prev': n_prev,
            'E_richardson': float(E_rich),
            'error_ha': float(error),
            'error_pct': float(error_pct),
        })
    return rich


def aitken_acceleration(results: list) -> list:
    """Apply Aitken's delta-squared acceleration.

    E_Aitken(n) = E*(n) - [E*(n+1) - E*(n)]^2 / [E*(n+2) - 2*E*(n+1) + E*(n)]
    """
    aitken = []
    for i in range(len(results) - 2):
        E0 = results[i]['E_var']
        E1 = results[i + 1]['E_var']
        E2 = results[i + 2]['E_var']

        denom = E2 - 2 * E1 + E0
        if abs(denom) < 1e-30:
            continue

        E_ait = E0 - (E1 - E0) ** 2 / denom
        error = E_ait - E_EXACT
        error_pct = abs(error / E_EXACT) * 100

        aitken.append({
            'n_max': results[i]['n_max'],
            'E_aitken': float(E_ait),
            'error_ha': float(error),
            'error_pct': float(error_pct),
        })
    return aitken


def k_star_analysis(results: list) -> dict:
    """Analyze the convergence of the k* sequence."""
    k_vals = [r['k_var'] for r in results]
    n_vals = [r['n_max'] for r in results]

    k_deltas = [k_vals[i + 1] - k_vals[i] for i in range(len(k_vals) - 1)]
    k_ratios = []
    for i in range(len(k_deltas) - 1):
        if abs(k_deltas[i + 1]) > 1e-15:
            k_ratios.append(k_deltas[i] / k_deltas[i + 1])
        else:
            k_ratios.append(None)

    # Estimate k_inf via Richardson on k*
    k_richardson = []
    for i in range(1, len(results)):
        n = results[i]['n_max']
        n_prev = results[i - 1]['n_max']
        k_n = results[i]['k_var']
        k_prev = results[i - 1]['k_var']
        # Assume p=2 for k* convergence (guess)
        for p in [2.0, 4.0]:
            np_c = n ** p
            np_p = n_prev ** p
            denom = np_c - np_p
            if abs(denom) > 1e-30:
                k_rich = (np_c * k_n - np_p * k_prev) / denom
                k_richardson.append({
                    'n': n, 'p': p, 'k_richardson': float(k_rich)
                })

    return {
        'k_values': k_vals,
        'n_values': n_vals,
        'k_deltas': k_deltas,
        'k_ratios': k_ratios,
        'k_richardson': k_richardson,
    }


def generate_plots(results: list, diff_data: dict, rich_p4: list,
                   aitken: list, k_data: dict) -> None:
    """Generate convergence plots."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("[PLOT] matplotlib not available — skipping plots")
        return

    plot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plots')
    os.makedirs(plot_dir, exist_ok=True)

    ns = [r['n_max'] for r in results]
    errors_pct = [r['error_pct'] for r in results]
    errors_ha = [r['E_var'] - E_EXACT for r in results]
    k_vals = [r['k_var'] for r in results]

    # --- Plot 1: Error vs n_max (log scale) ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    ax = axes[0]
    ax.semilogy(ns, errors_pct, 'bo-', markersize=8, label='Variational')
    if rich_p4:
        ax.semilogy([r['n'] for r in rich_p4],
                     [r['error_pct'] for r in rich_p4],
                     'rs-', markersize=6, label='Richardson (p=4)')
    if aitken:
        ax.semilogy([a['n_max'] for a in aitken],
                     [a['error_pct'] for a in aitken],
                     'g^-', markersize=6, label='Aitken Delta^2')
    ax.set_xlabel('n_max')
    ax.set_ylabel('Error (%)')
    ax.set_title('He Algebraic CI Convergence')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # --- Plot 2: Log-log error vs n ---
    ax = axes[1]
    ns_arr = np.array(ns[1:])  # skip n=1
    err_arr = np.array(errors_ha[1:])
    ax.loglog(ns_arr, err_arr, 'bo-', markersize=8)
    # Fit line
    if len(ns_arr) >= 2:
        log_n = np.log(ns_arr)
        log_err = np.log(err_arr)
        coeffs = np.polyfit(log_n, log_err, 1)
        fit_n = np.linspace(ns_arr[0], ns_arr[-1] * 1.5, 50)
        fit_err = np.exp(np.polyval(coeffs, np.log(fit_n)))
        ax.loglog(fit_n, fit_err, 'r--', label=f'slope = {coeffs[0]:.2f}')
    ax.set_xlabel('n_max')
    ax.set_ylabel('Error (Ha)')
    ax.set_title('Power-Law Convergence (n ≥ 2)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # --- Plot 3: k* vs n_max ---
    ax = axes[2]
    ax.plot(ns, k_vals, 'ko-', markersize=8)
    ax.set_xlabel('n_max')
    ax.set_ylabel('k_var')
    ax.set_title('Optimal Exponent k*(n_max)')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(plot_dir, 'track_di_sequence.png')
    plt.savefig(path, dpi=150)
    print(f"[PLOT] Saved: {path}")
    plt.close()


def main():
    """Run the full Sprint 3 analysis."""
    t_start = time.time()

    # Phase A: extend CI
    # Try up to n_max=7; fall back to 6 if too slow
    max_n = 7
    results = compute_sequence(n_max_max=max_n)

    # Phase B: sequence analysis
    print("\n" + "=" * 70)
    print("Phase B: Sequence Analysis")
    print("=" * 70)

    # Step 4-5: Difference analysis
    diff = difference_analysis(results)
    print("\n--- Energy Sequence ---")
    print(f"{'n_max':>5} {'E_var (Ha)':>14} {'Error (Ha)':>12} {'Error (%)':>10} "
          f"{'delta':>12} {'ratio':>8}")
    for i, r in enumerate(results):
        delta_str = f"{diff['deltas'][i]:.6f}" if i < len(diff['deltas']) else ""
        ratio_str = ""
        if i > 0 and i - 1 < len(diff['ratios']):
            rat = diff['ratios'][i - 1]
            ratio_str = f"{rat:.4f}" if rat is not None else "---"
        print(f"{r['n_max']:>5} {r['E_var']:>14.8f} {diff['errors'][i]:>12.6f} "
              f"{r['error_pct']:>10.4f} {delta_str:>12} {ratio_str:>8}")

    # Step 6: Power-law fit
    print("\n--- Power-Law Fit: error(n) = C/n^p (assuming E_inf = E_exact) ---")
    fit_2p = two_param_power_law_fit(results)
    print(f"  Status: {fit_2p['status']}")
    if fit_2p['status'] == 'converged':
        print(f"  C = {fit_2p['C']:.6f}")
        print(f"  p = {fit_2p['p']:.4f}")
        print(f"  RMS residual = {fit_2p['rms_residual']:.2e} Ha")

    print("\n--- Power-Law Fit: error(n) = delta_inf + C/n^p (3-parameter) ---")
    fit_3p = power_law_fit(results)
    print(f"  Status: {fit_3p['status']}")
    if fit_3p['status'] == 'converged':
        print(f"  E_inf = {fit_3p['E_inf']:.8f} Ha "
              f"(exact = {E_EXACT:.8f}, delta = {fit_3p['E_inf'] - E_EXACT:.2e})")
        print(f"  E_inf error = {fit_3p['E_inf_error_pct']:.4f}%")
        print(f"  C = {fit_3p['C']:.6f}")
        print(f"  p = {fit_3p['p']:.4f}")
        print(f"  RMS residual = {fit_3p['rms_residual']:.2e} Ha")

    # Step 7: Richardson extrapolation
    print("\n--- Richardson Extrapolation ---")
    for p_test in [2.0, 3.0, 4.0, 5.0]:
        rich = richardson_extrapolation(results, p_test)
        if rich:
            last = rich[-1]
            print(f"  p={p_test:.1f}: E_rich(n={last['n']}) = {last['E_richardson']:.8f} Ha, "
                  f"error = {last['error_pct']:.4f}%")

    # Use p=4 (Schwartz prediction) for detailed output
    rich_p4 = richardson_extrapolation(results, 4.0)
    if rich_p4:
        print(f"\n  Richardson (p=4) detail:")
        for r in rich_p4:
            print(f"    n={r['n']}: E = {r['E_richardson']:.8f}, "
                  f"error = {r['error_ha']:.6f} Ha ({r['error_pct']:.4f}%)")

    # Also try with fitted p
    if fit_3p['status'] == 'converged':
        p_fitted = fit_3p['p']
        rich_fitted = richardson_extrapolation(results, p_fitted)
        if rich_fitted:
            last = rich_fitted[-1]
            print(f"\n  Richardson (p={p_fitted:.2f}, fitted): "
                  f"E_rich(n={last['n']}) = {last['E_richardson']:.8f}, "
                  f"error = {last['error_pct']:.4f}%")

    # Step 8: Aitken acceleration
    print("\n--- Aitken Delta^2 Acceleration ---")
    aitken = aitken_acceleration(results)
    for a in aitken:
        print(f"  n_max={a['n_max']}: E_aitken = {a['E_aitken']:.8f}, "
              f"error = {a['error_ha']:.6f} Ha ({a['error_pct']:.4f}%)")

    # Step 9: k* sequence
    print("\n--- k* Sequence Analysis ---")
    k_data = k_star_analysis(results)
    print(f"  k* values: {[f'{k:.6f}' for k in k_data['k_values']]}")
    print(f"  k* deltas: {[f'{d:.6f}' for d in k_data['k_deltas']]}")
    if k_data['k_ratios']:
        print(f"  k* ratios: {[f'{r:.4f}' if r else '---' for r in k_data['k_ratios']]}")
    if k_data['k_richardson']:
        for kr in k_data['k_richardson']:
            if kr['n'] == max(r['n_max'] for r in results):
                print(f"  k* Richardson (p={kr['p']:.0f}): {kr['k_richardson']:.6f}")

    # Generate plots
    print("\n--- Generating Plots ---")
    generate_plots(results, diff, rich_p4, aitken, k_data)

    # Save all data
    output = {
        'metadata': {
            'track': 'DI Sprint 3',
            'description': 'Algebraic CI convergence sequence analysis',
            'E_exact_He': E_EXACT,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'total_time_s': time.time() - t_start,
        },
        'sequence': [{
            'n_max': r['n_max'],
            'n_spatial': r['n_spatial'],
            'n_configs': r['n_configs'],
            'k_var': r['k_var'],
            'E_var': r['E_var'],
            'error_pct': r['error_pct'],
            'error_ha': r['error_ha'],
            'build_time_s': r['build_time_s'],
            'opt_time_s': r['opt_time_s'],
        } for r in results],
        'difference_analysis': {
            'deltas': diff['deltas'],
            'ratios': [r if r is not None else None for r in diff['ratios']],
            'second_diffs': diff['second_diffs'],
        },
        'power_law_fit_2param': fit_2p,
        'power_law_fit_3param': fit_3p,
        'richardson_p4': rich_p4,
        'richardson_p_fitted': (
            richardson_extrapolation(results, fit_3p['p'])
            if fit_3p['status'] == 'converged' else []
        ),
        'aitken': aitken,
        'k_star_analysis': {
            'k_values': k_data['k_values'],
            'k_deltas': k_data['k_deltas'],
            'k_ratios': k_data['k_ratios'],
            'k_richardson': k_data['k_richardson'],
        },
    }

    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
    os.makedirs(data_dir, exist_ok=True)
    out_path = os.path.join(data_dir, 'track_di_sequence.json')
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\n[DATA] Saved: {out_path}")

    total_time = time.time() - t_start
    print(f"\n{'=' * 70}")
    print(f"Sprint 3 complete. Total time: {total_time:.1f}s")
    print(f"{'=' * 70}")


if __name__ == '__main__':
    main()
