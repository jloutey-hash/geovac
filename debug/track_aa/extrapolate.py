"""Track AA: Convergence extrapolation for D_e(l_max).

Fits D_e(l_max) to three models:
(a) D_e_inf - A/(l_max + 1/2)^2
(b) D_e_inf - A*exp(-B*l_max)
(c) D_e_inf - A/(l_max + 1/2)^4
"""
import json
import numpy as np
from scipy.optimize import curve_fit

D_E_EXACT = 0.17447  # Ha

def model_power2(l, D_inf, A):
    return D_inf - A / (l + 0.5)**2

def model_exp(l, D_inf, A, B):
    return D_inf - A * np.exp(-B * l)

def model_power4(l, D_inf, A):
    return D_inf - A / (l + 0.5)**4

def fit_and_report(lmax_arr, De_arr, label):
    print(f"\n{'='*60}")
    print(f"Extrapolation: {label}")
    print(f"{'='*60}")
    print(f"Data points: l_max = {lmax_arr}, D_e% = {[f'{x:.2f}' for x in De_arr]}")

    results = {}

    # (a) Power-2
    try:
        popt, pcov = curve_fit(model_power2, lmax_arr, De_arr, p0=[100.0, 50.0], maxfev=10000)
        D_inf, A = popt
        resid = np.sqrt(np.mean((model_power2(lmax_arr, *popt) - De_arr)**2))
        results['power2'] = {'D_inf': D_inf, 'A': A, 'rmse': resid}
        print(f"\n(a) D_inf - A/(l+1/2)^2: D_inf = {D_inf:.2f}%, A = {A:.2f}, RMSE = {resid:.3f}%")
    except Exception as e:
        print(f"\n(a) Power-2 fit failed: {e}")

    # (b) Exponential
    try:
        popt, pcov = curve_fit(model_exp, lmax_arr, De_arr, p0=[100.0, 50.0, 0.5], maxfev=10000)
        D_inf, A, B = popt
        resid = np.sqrt(np.mean((model_exp(lmax_arr, *popt) - De_arr)**2))
        results['exp'] = {'D_inf': D_inf, 'A': A, 'B': B, 'rmse': resid}
        print(f"(b) D_inf - A*exp(-B*l): D_inf = {D_inf:.2f}%, A = {A:.2f}, B = {B:.3f}, RMSE = {resid:.3f}%")
    except Exception as e:
        print(f"(b) Exponential fit failed: {e}")

    # (c) Power-4 (Schwartz-like)
    try:
        popt, pcov = curve_fit(model_power4, lmax_arr, De_arr, p0=[100.0, 500.0], maxfev=10000)
        D_inf, A = popt
        resid = np.sqrt(np.mean((model_power4(lmax_arr, *popt) - De_arr)**2))
        results['power4'] = {'D_inf': D_inf, 'A': A, 'rmse': resid}
        print(f"(c) D_inf - A/(l+1/2)^4: D_inf = {D_inf:.2f}%, A = {A:.2f}, RMSE = {resid:.3f}%")
    except Exception as e:
        print(f"(c) Power-4 fit failed: {e}")

    # Best fit
    if results:
        best = min(results.items(), key=lambda x: x[1]['rmse'])
        print(f"\nBest fit: model ({best[0]}), D_inf = {best[1]['D_inf']:.2f}%, RMSE = {best[1]['rmse']:.4f}%")

    return results


if __name__ == '__main__':
    with open('debug/track_aa/convergence_data.json') as f:
        data = json.load(f)

    # Extract D_e% for each method
    for method in ['adiabatic', '2d']:
        lmax_list = []
        de_list = []
        de_cusp_list = []
        for lm in range(1, 7):
            key = f'{method}_lmax{lm}'
            if key in data:
                lmax_list.append(lm)
                de_list.append(data[key]['D_e_pct'])
                de_cusp_list.append(data[key]['D_e_pct_cusp_corrected'])

        if len(lmax_list) >= 3:
            lmax_arr = np.array(lmax_list, dtype=float)
            de_arr = np.array(de_list)
            de_cusp_arr = np.array(de_cusp_list)

            fit_and_report(lmax_arr, de_arr, f'{method} (uncorrected)')
            fit_and_report(lmax_arr, de_cusp_arr, f'{method} (cusp-corrected)')

    # Also fit sigma-only (m_max=0) if we have that data
    # For now, just report the headline numbers
    print("\n\n" + "="*80)
    print("HEADLINE CONVERGENCE TABLE")
    print("="*80)
    print(f"{'l_max':>5} | {'Adiab D_e%':>11} | {'Adiab+cusp':>11} | {'2D D_e%':>11} | {'2D+cusp':>11}")
    print("-"*65)
    for lm in range(1, 7):
        a_key = f'adiabatic_lmax{lm}'
        d_key = f'2d_lmax{lm}'
        a_de = data[a_key]['D_e_pct'] if a_key in data else 0
        a_cusp = data[a_key]['D_e_pct_cusp_corrected'] if a_key in data else 0
        d_de = data[d_key]['D_e_pct'] if d_key in data else 0
        d_cusp = data[d_key]['D_e_pct_cusp_corrected'] if d_key in data else 0
        print(f"{lm:>5} | {a_de:>11.2f} | {a_cusp:>11.2f} | {d_de:>11.2f} | {d_cusp:>11.2f}")
