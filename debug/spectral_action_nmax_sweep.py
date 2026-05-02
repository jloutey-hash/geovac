"""Sweep n_max=3..8 for Gaussian spectral action convergence to K/pi."""
import sys, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))
from geovac.spectral_triple import FockSpectralTriple
from sympy import Rational

K_over_pi = 43.619934066848226
B = 42.0

results = []

for n in range(3, 9):
    print(f"\nn_max={n}:", end=" ", flush=True)
    try:
        st = FockSpectralTriple(n_max=n, kappa=Rational(-1,16),
                                j_type='kramers', adjacency_weights='cg')
    except Exception as e:
        print(f"FAILED to build: {e}")
        break

    D = np.array(st.dirac_operator.tolist(), dtype=float)
    evals = np.linalg.eigvalsh(D)
    evals_sq = evals**2
    dim_H = len(evals)
    max_abs = max(abs(evals))

    # Fine-grained search for K/pi Gaussian match
    # Start with coarse scan
    Ls_coarse = np.linspace(0.5, max(40, 2*max_abs), 20000)
    best_L, best_val, best_err = None, None, 1e10
    for L in Ls_coarse:
        g = float(np.sum(np.exp(-evals_sq / L**2)))
        err = abs(g - K_over_pi)
        if err < best_err:
            best_L, best_val, best_err = L, g, err

    # Refine around best
    Ls_fine = np.linspace(best_L - 0.5, best_L + 0.5, 50000)
    for L in Ls_fine:
        if L <= 0:
            continue
        g = float(np.sum(np.exp(-evals_sq / L**2)))
        err = abs(g - K_over_pi)
        if err < best_err:
            best_L, best_val, best_err = L, g, err

    rel = best_err / K_over_pi

    # Also search for B=42
    best_L2, best_val2, best_err2 = None, None, 1e10
    for L in Ls_coarse:
        g = float(np.sum(np.exp(-evals_sq / L**2)))
        err = abs(g - B)
        if err < best_err2:
            best_L2, best_val2, best_err2 = L, g, err
    Ls_fine2 = np.linspace(best_L2 - 0.5, best_L2 + 0.5, 50000)
    for L in Ls_fine2:
        if L <= 0:
            continue
        g = float(np.sum(np.exp(-evals_sq / L**2)))
        err = abs(g - B)
        if err < best_err2:
            best_L2, best_val2, best_err2 = L, g, err
    rel2 = best_err2 / B

    row = {
        'n_max': n, 'dim_H': dim_H, 'max_abs_lambda': max_abs,
        'Lambda_Kpi': best_L, 'val_Kpi': best_val, 'rel_err_Kpi': rel,
        'Lambda_B': best_L2, 'val_B': best_val2, 'rel_err_B': rel2,
    }
    results.append(row)

    print(f"dim_H={dim_H}, Lambda(K/pi)={best_L:.4f}, "
          f"rel_err={rel:.2e} ({rel*100:.4f}%), "
          f"Lambda/max|lam|={best_L/max_abs:.4f}")

print("\n" + "="*80)
print("Convergence table:")
print(f"{'n_max':>5} {'dim_H':>6} {'Lambda(K/pi)':>12} {'rel_err':>12} "
      f"{'Lambda/|lam|_max':>16} {'Lambda(B)':>10} {'rel_err_B':>12}")
for r in results:
    ratio = r['Lambda_Kpi'] / r['max_abs_lambda']
    print(f"{r['n_max']:>5} {r['dim_H']:>6} {r['Lambda_Kpi']:>12.4f} "
          f"{r['rel_err_Kpi']:>12.2e} {ratio:>16.4f} "
          f"{r['Lambda_B']:>10.4f} {r['rel_err_B']:>12.2e}")

# Check if rel_err is monotonically decreasing
errs = [r['rel_err_Kpi'] for r in results if r['rel_err_Kpi'] < 1]
if len(errs) >= 3:
    print(f"\nrel_err sequence: {[f'{e:.2e}' for e in errs]}")
    ratios = [errs[i+1]/errs[i] for i in range(len(errs)-1)]
    print(f"successive ratios: {[f'{r:.3f}' for r in ratios]}")
    if all(r < 1 for r in ratios):
        print("MONOTONICALLY DECREASING - convergent")
        # Estimate convergence rate
        import math
        log_ratios = [math.log(r) for r in ratios]
        avg_rate = sum(log_ratios)/len(log_ratios)
        print(f"average log(ratio) = {avg_rate:.3f} => convergence factor ~ {math.exp(avg_rate):.3f} per n_max step")
    else:
        print("NOT monotonically decreasing")

# Lambda convergence
lams = [r['Lambda_Kpi'] for r in results if r['rel_err_Kpi'] < 1]
if len(lams) >= 3:
    print(f"\nLambda(K/pi) sequence: {[f'{L:.4f}' for L in lams]}")
    diffs = [lams[i+1] - lams[i] for i in range(len(lams)-1)]
    print(f"successive diffs: {[f'{d:.4f}' for d in diffs]}")
