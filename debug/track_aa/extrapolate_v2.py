"""Track AA: Refined convergence extrapolation.

Key insight: l_max_per_m = {0: l_max, 1: 2} means pi channels are FIXED
at l_max_per_m[1]=2 for all l_max >= 2. Only sigma channels grow.
"""
import json
import numpy as np
from scipy.optimize import curve_fit

D_E_EXACT = 0.17447

with open('debug/track_aa/corrected_summary.json') as f:
    summary = json.load(f)

def model_power2(l, D_inf, A):
    return D_inf - A / (l + 0.5)**2

def model_exp(l, D_inf, A, B):
    return D_inf - A * np.exp(-B * l)

def model_power4(l, D_inf, A):
    return D_inf - A / (l + 0.5)**4

def model_power_n(l, D_inf, A, n):
    return D_inf - A / (l + 0.5)**n

def fit_all(lmax_arr, de_arr, label):
    print(f"\n{'='*60}")
    print(f"{label}")
    print(f"{'='*60}")
    print(f"l_max: {lmax_arr.astype(int).tolist()}")
    print(f"D_e%:  {[f'{x:.2f}' for x in de_arr]}")

    best_rmse = 1e10
    best_name = None
    best_D_inf = None

    for name, func, p0 in [
        ('(l+1/2)^-2', model_power2, [100.0, 50.0]),
        ('exp(-B*l)', model_exp, [100.0, 50.0, 0.5]),
        ('(l+1/2)^-4', model_power4, [100.0, 500.0]),
        ('(l+1/2)^-n', model_power_n, [100.0, 50.0, 3.0]),
    ]:
        try:
            popt, _ = curve_fit(func, lmax_arr, de_arr, p0=p0, maxfev=10000)
            pred = func(lmax_arr, *popt)
            resid = np.sqrt(np.mean((pred - de_arr)**2))
            D_inf = popt[0]
            params = ', '.join([f'{p:.3f}' for p in popt[1:]])
            print(f"  {name}: D_inf = {D_inf:.2f}%, params = [{params}], RMSE = {resid:.4f}%")

            if resid < best_rmse:
                best_rmse = resid
                best_name = name
                best_D_inf = D_inf
        except Exception as e:
            print(f"  {name}: FAILED ({e})")

    if best_name:
        print(f"  ** Best: {best_name}, D_inf = {best_D_inf:.2f}% (RMSE = {best_rmse:.4f}%)")

    return best_D_inf

# 2D solver fits
print("\n\n" + "#"*80)
print("# 2D SOLVER CONVERGENCE")
print("#"*80)

# All points l_max=1-6
lmax_all = np.array([1, 2, 3, 4, 5, 6], dtype=float)
de_all = np.array([summary[f'2d_lmax{lm}']['D_e_pct'] for lm in range(1, 7)])
de_cusp_all = np.array([summary[f'2d_lmax{lm}']['D_e_pct_pes_corrected'] for lm in range(1, 7)])

fit_all(lmax_all, de_all, "2D: all points (uncorrected)")
fit_all(lmax_all, de_cusp_all, "2D: all points (cusp-corrected)")

# Excluding l_max=1
lmax_2p = np.array([2, 3, 4, 5, 6], dtype=float)
de_2p = np.array([summary[f'2d_lmax{lm}']['D_e_pct'] for lm in range(2, 7)])
de_cusp_2p = np.array([summary[f'2d_lmax{lm}']['D_e_pct_pes_corrected'] for lm in range(2, 7)])

fit_all(lmax_2p, de_2p, "2D: l_max=2-6 (uncorrected)")
fit_all(lmax_2p, de_cusp_2p, "2D: l_max=2-6 (cusp-corrected)")

# Even l_max only (where pi channels change)
lmax_even = np.array([2, 4, 6], dtype=float)
de_even = np.array([summary[f'2d_lmax{lm}']['D_e_pct'] for lm in [2, 4, 6]])
de_cusp_even = np.array([summary[f'2d_lmax{lm}']['D_e_pct_pes_corrected'] for lm in [2, 4, 6]])

fit_all(lmax_even, de_even, "2D: even l_max only (uncorrected)")
fit_all(lmax_even, de_cusp_even, "2D: even l_max only (cusp-corrected)")

# Incremental gains
print("\n\nINCREMENTAL GAINS (2D uncorrected)")
print("="*50)
for i in range(1, len(de_all)):
    delta = de_all[i] - de_all[i-1]
    print(f"  l_max {int(lmax_all[i-1])}->{int(lmax_all[i])}: +{delta:.2f} pp")

print("\nINCREMENTAL GAINS (2D cusp-corrected)")
print("="*50)
for i in range(1, len(de_cusp_all)):
    delta = de_cusp_all[i] - de_cusp_all[i-1]
    print(f"  l_max {int(lmax_all[i-1])}->{int(lmax_all[i])}: +{delta:.2f} pp")

# Adiabatic fits (for comparison)
print("\n\n" + "#"*80)
print("# ADIABATIC SOLVER CONVERGENCE")
print("#"*80)

lmax_all_a = np.array([1, 2, 3, 4, 5, 6], dtype=float)
de_all_a = np.array([summary[f'adiabatic_lmax{lm}']['D_e_pct'] for lm in range(1, 7)])
fit_all(lmax_all_a, de_all_a, "Adiabatic: all points (uncorrected)")

# Gap (adiabatic - 2D)
print("\n\nADIABATIC - 2D GAP")
print("="*50)
for lm in range(1, 7):
    gap = summary[f'adiabatic_lmax{lm}']['D_e_pct'] - summary[f'2d_lmax{lm}']['D_e_pct']
    print(f"  l_max={lm}: {gap:.2f} pp")

# Paper 15 comparison
print("\n\nCOMPARISON WITH PAPER 15 TABLE III")
print("="*60)
print(f"{'l_max':>5} {'m_max':>5} | {'Paper15 2D':>10} {'This work':>10} | {'Paper15 Ad':>10} {'This work':>10}")
print("-"*65)
paper15 = {
    (2, 0): (79.5, 90.3),
    (2, 1): (86.2, 96.8),
    (3, 0): (80.2, 91.0),
    (3, 1): (87.1, 97.7),
    (4, 0): (87.0, 98.6),
    (4, 1): (94.1, 105.2),
}
for (lm, mm), (p2d, pad) in sorted(paper15.items()):
    key = f'2d_lmax{lm}'
    if key in summary:
        t2d = summary[key]['D_e_pct']
        tad = summary[f'adiabatic_lmax{lm}']['D_e_pct']
        print(f"{lm:>5} {mm:>5} | {p2d:>10.1f} {t2d:>10.2f} | {pad:>10.1f} {tad:>10.2f}")
    else:
        print(f"{lm:>5} {mm:>5} | {p2d:>10.1f} {'N/A':>10} | {pad:>10.1f} {'N/A':>10}")

print("\nNote: Paper 15 uses n_alpha=200, n_Re=400. This work uses n_alpha=60, n_Re=120.")
print("Small discrepancies expected from grid resolution.")
