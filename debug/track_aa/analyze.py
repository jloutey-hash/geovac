"""Analyze convergence data with correct D_e reference."""
import json
import numpy as np
from scipy.optimize import curve_fit

D_E_EXACT = 0.17447
E_ATOMS = -1.0  # Exact two-atom dissociation limit

with open('debug/track_aa/convergence_data.json') as f:
    data = json.load(f)

# Also apply cusp correction properly
from geovac.cusp_correction import cusp_correction_h2_point

print("CORRECTED CONVERGENCE TABLE")
print("Using E_atoms = -1.0 Ha (exact), D_e_exact = 0.17447 Ha")
print("="*95)
print(f"{'Method':<12} {'l_max':>5} {'N_ch':>5} {'E_min':>10} {'D_e':>8} {'D_e%':>7} {'dE_cusp':>8} {'D_e%+c':>7} {'Note':<10}")
print("="*95)

summary = {}

for method in ['2d', 'adiabatic']:
    for lm in range(1, 7):
        key = f'{method}_lmax{lm}'
        if key not in data:
            continue
        r = data[key]
        E_min = r['E_min']
        R_eq = r['R_eq_approx']

        # Correct D_e using exact dissociation limit
        D_e = E_ATOMS - E_min
        D_e_pct = D_e / D_E_EXACT * 100

        # Cusp correction at R_eq
        dE_cusp = cusp_correction_h2_point(R_eq, l_max=lm)
        E_min_corr = E_min + dE_cusp
        D_e_corr = E_ATOMS - E_min_corr
        D_e_pct_corr = D_e_corr / D_E_EXACT * 100

        note = 'VAR VIOL' if D_e_pct > 100 else ''

        # Count channels
        from geovac.level4_multichannel import _channel_list_extended
        ch = _channel_list_extended(l_max=lm, m_max=1, l_max_per_m={0: lm, 1: min(lm, 2)}, homonuclear=True)
        n_ch = len(ch)

        print(f"{method:<12} {lm:>5} {n_ch:>5} {E_min:>10.6f} {D_e:>8.5f} {D_e_pct:>7.2f} {dE_cusp:>8.5f} {D_e_pct_corr:>7.2f} {note:<10}")

        summary[key] = {
            'method': method, 'l_max': lm, 'n_ch': n_ch,
            'E_min': E_min, 'R_eq': R_eq,
            'D_e': D_e, 'D_e_pct': D_e_pct,
            'dE_cusp_at_Req': dE_cusp,
            'D_e_corrected': D_e_corr,
            'D_e_pct_corrected': D_e_pct_corr,
        }
    print()

# Now do proper PES-level cusp correction for the 2D results
print("\n\nDETAILED PES CUSP CORRECTION (2D solver)")
print("="*80)
for lm in range(1, 7):
    key = f'2d_lmax{lm}'
    if key not in data:
        continue
    r = data[key]
    R_arr = np.array(r['R_grid'])
    E_arr = np.array(r['E_pes'])

    # Apply cusp at each R
    dE = np.array([cusp_correction_h2_point(R, l_max=lm) for R in R_arr])
    E_corr = E_arr + dE

    # D_e from E_atoms = -1.0
    D_e_raw = E_ATOMS - np.min(E_arr)
    D_e_corr = E_ATOMS - np.min(E_corr)

    print(f"l_max={lm}: D_e_raw = {D_e_raw/D_E_EXACT*100:.2f}%, D_e_cusp = {D_e_corr/D_E_EXACT*100:.2f}%")

    summary[key]['D_e_pes_corrected'] = float(D_e_corr)
    summary[key]['D_e_pct_pes_corrected'] = float(D_e_corr / D_E_EXACT * 100)

print("\n\nDETAILED PES CUSP CORRECTION (adiabatic solver)")
print("="*80)
for lm in range(1, 7):
    key = f'adiabatic_lmax{lm}'
    if key not in data:
        continue
    r = data[key]
    R_arr = np.array(r['R_grid'])
    E_arr = np.array(r['E_pes'])

    dE = np.array([cusp_correction_h2_point(R, l_max=lm) for R in R_arr])
    E_corr = E_arr + dE

    D_e_raw = E_ATOMS - np.min(E_arr)
    D_e_corr = E_ATOMS - np.min(E_corr)

    print(f"l_max={lm}: D_e_raw = {D_e_raw/D_E_EXACT*100:.2f}%, D_e_cusp = {D_e_corr/D_E_EXACT*100:.2f}%")

    summary[key]['D_e_pes_corrected'] = float(D_e_corr)
    summary[key]['D_e_pct_pes_corrected'] = float(D_e_corr / D_E_EXACT * 100)


# Convergence extrapolation
print("\n\n" + "="*80)
print("CONVERGENCE EXTRAPOLATION")
print("="*80)

def model_power2(l, D_inf, A):
    return D_inf - A / (l + 0.5)**2

def model_exp(l, D_inf, A, B):
    return D_inf - A * np.exp(-B * l)

def model_power4(l, D_inf, A):
    return D_inf - A / (l + 0.5)**4

for method in ['2d', 'adiabatic']:
    lmax_arr = []
    de_arr = []
    de_cusp_arr = []
    for lm in range(1, 7):
        key = f'{method}_lmax{lm}'
        if key in summary:
            lmax_arr.append(lm)
            de_arr.append(summary[key]['D_e_pct'])
            if 'D_e_pct_pes_corrected' in summary[key]:
                de_cusp_arr.append(summary[key]['D_e_pct_pes_corrected'])

    lmax_arr = np.array(lmax_arr, dtype=float)
    de_arr = np.array(de_arr)

    for label, darr in [('uncorrected', de_arr), ('cusp-corrected', np.array(de_cusp_arr) if de_cusp_arr else None)]:
        if darr is None or len(darr) < 3:
            continue
        print(f"\n--- {method} ({label}) ---")
        print(f"  l_max: {lmax_arr.astype(int).tolist()}")
        print(f"  D_e%:  {[f'{x:.2f}' for x in darr]}")

        for name, func, p0 in [
            ('(l+1/2)^-2', model_power2, [100.0, 50.0]),
            ('exp(-B*l)', model_exp, [100.0, 50.0, 0.5]),
            ('(l+1/2)^-4', model_power4, [100.0, 500.0]),
        ]:
            try:
                popt, _ = curve_fit(func, lmax_arr, darr, p0=p0, maxfev=10000)
                resid = np.sqrt(np.mean((func(lmax_arr, *popt) - darr)**2))
                D_inf = popt[0]
                print(f"  {name}: D_inf = {D_inf:.2f}%, RMSE = {resid:.3f}%")
            except Exception as e:
                print(f"  {name}: FAILED ({e})")


# Final headline table
print("\n\n" + "="*80)
print("HEADLINE CONVERGENCE TABLE")
print("="*80)
print(f"{'l_max':>5} {'N_ch':>5} | {'2D D_e%':>9} {'2D+cusp':>9} | {'Adiab D_e%':>11} {'Adiab+cusp':>11}")
print("-"*65)
for lm in range(1, 7):
    d_key = f'2d_lmax{lm}'
    a_key = f'adiabatic_lmax{lm}'

    from geovac.level4_multichannel import _channel_list_extended
    ch = _channel_list_extended(l_max=lm, m_max=1, l_max_per_m={0: lm, 1: min(lm, 2)}, homonuclear=True)
    n_ch = len(ch)

    d_de = summary.get(d_key, {}).get('D_e_pct', 0)
    d_cusp = summary.get(d_key, {}).get('D_e_pct_pes_corrected', 0)
    a_de = summary.get(a_key, {}).get('D_e_pct', 0)
    a_cusp = summary.get(a_key, {}).get('D_e_pct_pes_corrected', 0)

    print(f"{lm:>5} {n_ch:>5} | {d_de:>9.2f} {d_cusp:>9.2f} | {a_de:>11.2f} {a_cusp:>11.2f}")

# Save corrected summary
with open('debug/track_aa/corrected_summary.json', 'w') as f:
    json.dump({k: v for k, v in summary.items()}, f, indent=2)
print("\nSaved corrected summary to debug/track_aa/corrected_summary.json")
