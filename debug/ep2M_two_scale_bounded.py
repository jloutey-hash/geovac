"""
Track EP-2M: two-scale bounded form for the universal entropy curve.

EP-2h showed log(2)*tanh^2(c*x^gamma) fails (R^2=0.37 vs power law 0.90).
The Taylor expansion forces small-x slope = 2*gamma with prefactor c^2,
incompatible with the 4-decade dynamic range.

Three candidate two-scale forms (each with 3 parameters):

  (P)  S = A * x^gamma                     (pure power law, baseline)
  (T)  S = A * x^gamma / (1 + (A*x^gamma/L)^p)^(1/p)   (rational-saturate)
  (E)  S = L * (1 - exp(-c*x^gamma))^p     (stretched exponential)

where L = log(2) is the saturation ceiling.

Fit each on:
  - EP-2g 22 points (single-center He-like + LiH bond R-sweep)
  - Plus EP-2k HF + H2O R-sweeps (10 more points)

Compare R^2 and max relative residual.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))


def load_data():
    """Load (S_B, w_B/delta_B) from EP-2g, EP-2k."""
    rows = []
    # EP-2g
    g_path = os.path.join(os.path.dirname(__file__), 'data',
                          'ep2g_two_variable_scaling.json')
    with open(g_path) as f:
        g = json.load(f)
    for r in g['rows']:
        rows.append({'S': r['S_B'],
                     'x': r['w_B_dimless'] / r['delta_B'],
                     'src': r['kind']})
    # EP-2k
    k_path = os.path.join(os.path.dirname(__file__), 'data',
                          'ep2k_heavy_atom_bond_sweep.json')
    with open(k_path) as f:
        k = json.load(f)
    for mol in ('HF', 'H2O'):
        for r in k[mol]['rows']:
            if r['delta_B'] < 1e-10:
                continue
            rows.append({'S': r['S_B'],
                         'x': r['w_B_dimless'] / r['delta_B'],
                         'src': f'bond_{mol}'})
    return rows


def fit_power(rows):
    xs = np.array([r['x'] for r in rows])
    ys = np.array([r['S'] for r in rows])
    lx = np.log(xs); ly = np.log(ys)
    g, b = np.polyfit(lx, ly, 1)
    pred = g * lx + b
    r2 = 1 - np.sum((ly - pred) ** 2) / np.sum((ly - ly.mean()) ** 2)
    A = float(np.exp(b))
    pred_S = A * xs ** g
    max_rel = float(np.max(np.abs(pred_S - ys) / ys))
    return {'name': 'power_law', 'A': A, 'gamma': float(g),
            'r2_log': float(r2), 'max_rel': max_rel}


def fit_rational_sat(rows, L=np.log(2.0)):
    """S = A x^g / (1 + (A x^g / L)^p)^(1/p), fit in log-space."""
    xs = np.array([r['x'] for r in rows])
    ys = np.array([r['S'] for r in rows])
    def model(x, A, g, p):
        S0 = A * x ** g
        return S0 / (1 + (S0 / L) ** p) ** (1 / p)
    def log_model(x, A, g, p):
        return np.log(model(x, A, g, p))
    try:
        popt, _ = curve_fit(log_model, xs, np.log(ys),
                            p0=[0.16, 2.2, 2.0],
                            maxfev=20000,
                            bounds=([1e-6, 0.5, 0.5], [100, 10, 20]))
        pred = model(xs, *popt)
        ly = np.log(ys); lp = np.log(pred)
        r2_log = 1 - np.sum((ly - lp) ** 2) / np.sum((ly - ly.mean()) ** 2)
        max_rel = float(np.max(np.abs(pred - ys) / ys))
        return {'name': 'rational_saturate', 'A': float(popt[0]),
                'gamma': float(popt[1]), 'p': float(popt[2]),
                'r2_log': float(r2_log), 'max_rel': max_rel}
    except Exception as e:
        return {'name': 'rational_saturate', 'error': str(e)}


def fit_stretched_exp(rows, L=np.log(2.0)):
    """S = L * (1 - exp(-c*x^g))^p, fit in log-space."""
    xs = np.array([r['x'] for r in rows])
    ys = np.array([r['S'] for r in rows])
    def model(x, c, g, p):
        # clip to avoid log(0) when fit explores 1-exp(0)=0
        return L * (1 - np.exp(-c * x ** g)) ** p
    def log_model(x, c, g, p):
        m = model(x, c, g, p)
        return np.log(np.maximum(m, 1e-300))
    try:
        popt, _ = curve_fit(log_model, xs, np.log(ys),
                            p0=[0.1, 2.0, 1.0],
                            maxfev=20000,
                            bounds=([1e-6, 0.3, 0.1], [100, 10, 20]))
        pred = model(xs, *popt)
        ly = np.log(ys); lp = np.log(pred)
        r2_log = 1 - np.sum((ly - lp) ** 2) / np.sum((ly - ly.mean()) ** 2)
        max_rel = float(np.max(np.abs(pred - ys) / ys))
        return {'name': 'stretched_exp', 'c': float(popt[0]),
                'gamma': float(popt[1]), 'p': float(popt[2]),
                'r2_log': float(r2_log), 'max_rel': max_rel}
    except Exception as e:
        return {'name': 'stretched_exp', 'error': str(e)}


def fit_hill(rows, L=np.log(2.0)):
    """S = L * (x/x0)^g / (1 + (x/x0)^g)  — 2-parameter Hill function."""
    xs = np.array([r['x'] for r in rows])
    ys = np.array([r['S'] for r in rows])
    def model(x, x0, g):
        u = (x / x0) ** g
        return L * u / (1 + u)
    def log_model(x, x0, g):
        return np.log(np.maximum(model(x, x0, g), 1e-300))
    try:
        popt, _ = curve_fit(log_model, xs, np.log(ys),
                            p0=[1.0, 2.2],
                            maxfev=20000,
                            bounds=([1e-3, 0.3], [100, 10]))
        pred = model(xs, *popt)
        ly = np.log(ys); lp = np.log(pred)
        r2_log = 1 - np.sum((ly - lp) ** 2) / np.sum((ly - ly.mean()) ** 2)
        max_rel = float(np.max(np.abs(pred - ys) / ys))
        # Equivalent A in S ~ A x^g at small x.
        A_eff = L / popt[0] ** popt[1]
        return {'name': 'hill', 'x0': float(popt[0]),
                'gamma': float(popt[1]), 'A_eff': float(A_eff),
                'r2_log': float(r2_log), 'max_rel': max_rel}
    except Exception as e:
        return {'name': 'hill', 'error': str(e)}


def main():
    rows = load_data()
    n = len(rows)
    print(f"EP-2M: two-scale fit on {n} (S, x=w/d) points")
    print("=" * 60)

    fits = [fit_power(rows), fit_rational_sat(rows),
            fit_stretched_exp(rows), fit_hill(rows)]
    for f in fits:
        if 'error' in f:
            print(f"\n{f['name']}: ERROR {f['error']}")
            continue
        print(f"\n{f['name']}:")
        for k, v in f.items():
            if k == 'name':
                continue
            print(f"  {k} = {v:.4f}" if isinstance(v, float) else f"  {k} = {v}")

    out = {'rows': rows, 'fits': fits, 'L': float(np.log(2.0))}
    out_path = os.path.join(os.path.dirname(__file__), 'data',
                            'ep2M_two_scale_bounded.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as fh:
        json.dump(out, fh, indent=2)
    print(f"\nSaved -> {out_path}")


if __name__ == '__main__':
    main()
