"""
Overlap diagnostic: compute inter-fiber channel overlap S_avg(R) for BeH2.

This script measures how much the ground adiabatic channels of the two
Be-H bond fibers overlap at the Be center, as a function of R.

If S_avg(R) has a meaningful negative slope in the equilibrium region,
exchange coupling (which scales as ~S^2) could provide the differential
R-dependence needed to fix R_eq. If S_avg is flat, exchange will fail
for the same reason the monopole did.

Output:
  - Table: R, S_avg, E_total, dS/dR
  - Plot: debug/plots/overlap_vs_R.png
  - Data: debug/data/overlap_diagnostic.json
"""

import json
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from geovac.core_screening import CoreScreening
from geovac.ab_initio_pk import AbInitioPK
from geovac.level4_multichannel import solve_level4_h2_multichannel
from geovac.inter_fiber_coupling import compute_overlap_diagnostic

# --- Configuration ---
Z_CENTER = 4      # Be
Z_LIGAND = 1      # H
N_CORE = 2
L_MAX = 2
N_ALPHA = 100
N_RE = 300
R_EQ_EXPT = 2.507  # bohr

R_VALUES = np.array([
    2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 3.0,
    3.3, 3.6, 4.0, 4.5, 5.0, 5.5, 6.0,
])


def main():
    print("=" * 72)
    print("  Inter-Fiber Channel Overlap Diagnostic — BeH2")
    print("=" * 72)

    # Step 1: Solve core + derive PK
    print("\n[1] Solving core (Be 1s^2)...")
    t0 = time.time()
    core = CoreScreening(Z=Z_CENTER, l_max=L_MAX, n_alpha=200)
    core.solve(verbose=False)
    E_core = core.energy
    print(f"    E_core = {E_core:.6f} Ha  ({time.time()-t0:.1f}s)")

    pk = AbInitioPK(core, n_core=N_CORE)
    pk_potentials = [{'C_core': pk.A, 'beta_core': pk.B, 'atom': 'A'}]
    Z_eff = Z_CENTER - N_CORE
    print(f"    Z_eff = {Z_eff}, PK: A={pk.A:.4f}, B={pk.B:.4f}")

    # Step 2: Scan
    print(f"\n[2] Computing overlap at {len(R_VALUES)} R-points...")
    print(f"    {'R':>6s}  {'S_avg':>8s}  {'E_bond':>10s}  {'E_total':>12s}"
          f"  {'time':>6s}")
    print(f"    {'-'*6}  {'-'*8}  {'-'*10}  {'-'*12}  {'-'*6}")

    results = []
    for R in R_VALUES:
        ti = time.time()

        # Solve bond pair
        l4_result = solve_level4_h2_multichannel(
            R=R, Z_A=Z_eff, Z_B=Z_LIGAND,
            l_max=L_MAX, n_alpha=N_ALPHA, n_Re=N_RE,
            verbose=False, pk_potentials=pk_potentials,
        )
        E_bond = l4_result['E_elec']

        # Total energy (block-diagonal, no interbond)
        V_NN = 2.0 * Z_CENTER * Z_LIGAND / R + Z_LIGAND**2 / (2.0 * R)
        zr = Z_CENTER * R
        V_cross = -2.0 * N_CORE * Z_LIGAND * (
            (1.0/R) * (1.0 - (1.0 + zr) * np.exp(-2.0 * zr)))
        E_total = E_core + V_cross + 2.0 * E_bond + V_NN

        # Overlap diagnostic
        ovlp = compute_overlap_diagnostic(
            R=R, Z_A=Z_eff, Z_B=Z_LIGAND,
            l_max=L_MAX, n_alpha=N_ALPHA,
            level4_result=l4_result,
            pk_potentials=pk_potentials,
            n_sample_Re=15,
        )

        dt = time.time() - ti
        print(f"    {R:6.3f}  {ovlp['S_avg']:8.5f}  {2*E_bond:10.4f}"
              f"  {E_total:12.6f}  {dt:6.1f}s")

        results.append({
            'R': float(R),
            'S_avg': ovlp['S_avg'],
            'E_bond_total': float(2 * E_bond),
            'E_total': float(E_total),
            'channel_weights': ovlp['channel_weights'],
            'S_samples': ovlp['S_samples'],
            'Re_samples': ovlp['Re_samples'],
        })

    # Step 3: Compute dS/dR
    R_arr = np.array([r['R'] for r in results])
    S_arr = np.array([r['S_avg'] for r in results])
    E_arr = np.array([r['E_total'] for r in results])

    # Central differences for interior, forward/backward at edges
    dS_dR = np.gradient(S_arr, R_arr)

    print(f"\n[3] Overlap slope analysis")
    print(f"    {'R':>6s}  {'S_avg':>8s}  {'dS/dR':>8s}  {'E_total':>12s}")
    print(f"    {'-'*6}  {'-'*8}  {'-'*8}  {'-'*12}")
    for i in range(len(R_arr)):
        print(f"    {R_arr[i]:6.3f}  {S_arr[i]:8.5f}  {dS_dR[i]:8.5f}"
              f"  {E_arr[i]:12.6f}")

    # Key metrics
    i_eq_region = np.argmin(np.abs(R_arr - R_EQ_EXPT))
    S_short = S_arr[0]
    S_long = S_arr[-1]
    S_range = S_short - S_long
    slope_at_eq = dS_dR[i_eq_region]

    print(f"\n[4] Summary")
    print(f"    S(R={R_arr[0]:.1f}) = {S_short:.5f}")
    print(f"    S(R={R_arr[-1]:.1f}) = {S_long:.5f}")
    print(f"    Range: {S_range:.5f}")
    print(f"    dS/dR near R_eq({R_EQ_EXPT}): {slope_at_eq:.5f}")

    if abs(S_range) < 0.01:
        print(f"\n    ** S_avg is nearly R-independent (range < 0.01).")
        print(f"    ** Exchange will NOT fix R_eq — same failure mode as monopole.")
        exchange_viable = False
    else:
        print(f"\n    ** S_avg has meaningful R-dependence (range = {S_range:.4f}).")
        print(f"    ** Exchange is a viable candidate for R_eq correction.")
        exchange_viable = True

    # Step 5: If viable, fit exchange model E_exch(R) = -K * S_avg(R)^2
    fit_result = None
    if exchange_viable:
        from scipy.interpolate import CubicSpline as CS
        from scipy.optimize import minimize_scalar

        print(f"\n[5] Exchange feasibility fit (spline-interpolated)")

        # Spline the PES and overlap for continuous optimization
        E_spl = CS(R_arr, E_arr)
        S_spl = CS(R_arr, S_arr)

        # Block-diagonal minimum (interpolated)
        res_bd = minimize_scalar(
            lambda r: E_spl(r), bounds=(R_arr[0], R_arr[-1]),
            method='bounded')
        R_eq_bd = res_bd.x

        # Scan K with S^2 model
        best_K_sq = None
        best_err_sq = 999.0
        best_Req_sq = None
        for K in np.linspace(0.01, 15.0, 1500):
            res = minimize_scalar(
                lambda r, _K=K: E_spl(r) - _K * S_spl(r)**2,
                bounds=(R_arr[0], R_arr[-1]), method='bounded')
            err = abs(res.x - R_EQ_EXPT) / R_EQ_EXPT * 100
            if err < best_err_sq:
                best_err_sq = err
                best_K_sq = K
                best_Req_sq = res.x

        # Scan K with linear S model
        best_K_lin = None
        best_err_lin = 999.0
        best_Req_lin = None
        for K in np.linspace(0.01, 15.0, 1500):
            res = minimize_scalar(
                lambda r, _K=K: E_spl(r) - _K * S_spl(r),
                bounds=(R_arr[0], R_arr[-1]), method='bounded')
            err = abs(res.x - R_EQ_EXPT) / R_EQ_EXPT * 100
            if err < best_err_lin:
                best_err_lin = err
                best_K_lin = K
                best_Req_lin = res.x

        print(f"    Block-diagonal R_eq = {R_eq_bd:.4f} bohr"
              f"  ({abs(R_eq_bd - R_EQ_EXPT)/R_EQ_EXPT*100:.1f}% error)")
        print(f"")
        print(f"    Model: E_exch = -K * S(R)^2")
        print(f"    Best K = {best_K_sq:.4f} Ha")
        print(f"    Corrected R_eq = {best_Req_sq:.4f} bohr"
              f"  ({best_err_sq:.1f}% error)")
        print(f"")
        print(f"    Model: E_exch = -K * S(R)")
        print(f"    Best K = {best_K_lin:.4f} Ha")
        print(f"    Corrected R_eq = {best_Req_lin:.4f} bohr"
              f"  ({best_err_lin:.1f}% error)")
        print(f"")
        print(f"    Target: < 15% error")

        # Use the better of the two models
        if best_err_sq <= best_err_lin:
            best_K, best_err, best_Req, model_p = (
                best_K_sq, best_err_sq, best_Req_sq, 2)
        else:
            best_K, best_err, best_Req, model_p = (
                best_K_lin, best_err_lin, best_Req_lin, 1)

        if best_err < 15.0:
            print(f"    ** SUCCESS: Exchange achieves {best_err:.1f}%"
                  f" R_eq error (K={best_K:.4f}, p={model_p})")
        else:
            print(f"    ** Exchange model best: {best_err:.1f}%"
                  f" (K={best_K:.4f}, p={model_p})")

        # Compute corrected PES at grid points for plotting
        E_corrected = E_arr - best_K * S_arr**model_p
        # Also compute on a fine grid for the plot
        R_fine = np.linspace(R_arr[0], R_arr[-1], 200)
        E_corr_fine = E_spl(R_fine) - best_K * S_spl(R_fine)**model_p

        fit_result = {
            'K_best': float(best_K),
            'model_power': model_p,
            'R_eq_corrected': float(best_Req),
            'R_eq_error_pct': float(best_err),
            'R_eq_block_diag': float(R_eq_bd),
            'E_corrected': E_corrected.tolist(),
            'R_fine': R_fine.tolist(),
            'E_corr_fine': E_corr_fine.tolist(),
            'K_sq': float(best_K_sq),
            'err_sq': float(best_err_sq),
            'K_lin': float(best_K_lin),
            'err_lin': float(best_err_lin),
        }

    # Step 6: Channel weight breakdown at a few R values
    print(f"\n[6] Channel weight breakdown (|c_{{l1,l2}}|^2)")
    for idx in [0, len(results)//2, -1]:
        r = results[idx]
        print(f"\n    R = {r['R']:.2f} bohr:")
        cw = r['channel_weights']
        for ch_key in sorted(cw.keys(), key=lambda x: -cw[x]):
            l1, l2 = ch_key
            parity = "+" if (l1+l2) % 2 == 0 else "-"
            print(f"      ({l1},{l2}) [{parity}]: {cw[ch_key]:.5f}")

    # Step 7: Save data
    out_dir = Path(__file__).resolve().parent / 'data'
    out_dir.mkdir(exist_ok=True)
    save_data = {
        'R': R_arr.tolist(),
        'S_avg': S_arr.tolist(),
        'dS_dR': dS_dR.tolist(),
        'E_total': E_arr.tolist(),
        'exchange_viable': exchange_viable,
        'S_range': float(S_range),
        'slope_at_eq': float(slope_at_eq),
    }
    if fit_result:
        save_data['exchange_fit'] = fit_result

    # Convert channel weight keys to strings for JSON
    for r in results:
        r['channel_weights'] = {
            f"({k[0]},{k[1]})": v for k, v in r['channel_weights'].items()
        }
    save_data['per_R'] = results

    data_path = out_dir / 'overlap_diagnostic.json'
    with open(data_path, 'w') as f:
        json.dump(save_data, f, indent=2)
    print(f"\n    Data saved: {data_path}")

    # Step 8: Plot
    plot_dir = Path(__file__).resolve().parent / 'plots'
    plot_dir.mkdir(exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    fig.suptitle('Inter-Fiber Channel Overlap Diagnostic — BeH₂', fontsize=14)

    # (a) S_avg vs R
    ax = axes[0, 0]
    ax.plot(R_arr, S_arr, 'bo-', linewidth=2, markersize=6)
    ax.axvline(R_EQ_EXPT, color='r', linestyle='--', alpha=0.5,
               label=f'R_eq expt = {R_EQ_EXPT}')
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('S_avg(R)')
    ax.set_title('(a) Inter-fiber overlap')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # (b) dS/dR vs R
    ax = axes[0, 1]
    ax.plot(R_arr, dS_dR, 'rs-', linewidth=2, markersize=6)
    ax.axhline(0, color='k', linestyle='-', alpha=0.3)
    ax.axvline(R_EQ_EXPT, color='r', linestyle='--', alpha=0.5)
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('dS/dR')
    ax.set_title('(b) Overlap slope')
    ax.grid(True, alpha=0.3)

    # (c) Block-diagonal PES
    ax = axes[1, 0]
    ax.plot(R_arr, E_arr, 'k-', linewidth=2, label='Block-diagonal')
    if fit_result and 'R_fine' in fit_result:
        ax.plot(fit_result['R_fine'], fit_result['E_corr_fine'],
                'g-', linewidth=2,
                label=f'With exchange (K={fit_result["K_best"]:.3f},'
                      f' p={fit_result["model_power"]})')
    elif fit_result:
        ax.plot(R_arr, fit_result['E_corrected'], 'g--', linewidth=2,
                label=f'With exchange (K={fit_result["K_best"]:.3f})')
    ax.axvline(R_EQ_EXPT, color='r', linestyle='--', alpha=0.5,
               label=f'R_eq expt')
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('E_total (Ha)')
    ax.set_title('(c) PES comparison')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (d) Channel weights vs R (top 4 channels)
    ax = axes[1, 1]
    # Collect channel labels across all R
    all_chs = set()
    for r in results:
        all_chs.update(r['channel_weights'].keys())
    # Sort by average weight
    ch_avg = {}
    for ch in all_chs:
        ch_avg[ch] = np.mean([r['channel_weights'].get(ch, 0.0)
                              for r in results])
    top_chs = sorted(ch_avg.keys(), key=lambda x: -ch_avg[x])[:5]
    colors = ['b', 'r', 'g', 'm', 'c']
    for i, ch in enumerate(top_chs):
        vals = [r['channel_weights'].get(ch, 0.0) for r in results]
        ax.plot(R_arr, vals, f'{colors[i]}o-', linewidth=1.5,
                markersize=4, label=ch)
    ax.set_xlabel('R (bohr)')
    ax.set_ylabel('|c_{l1,l2}|²')
    ax.set_title('(d) Channel weights vs R')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_path = plot_dir / 'overlap_vs_R.png'
    plt.savefig(plot_path, dpi=150)
    print(f"    Plot saved: {plot_path}")
    plt.close()

    print(f"\n{'='*72}")
    print(f"  Done.")
    print(f"{'='*72}")


if __name__ == '__main__':
    main()
