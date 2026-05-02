"""n_max=5 probe for spectral action convergence trend."""
import sys, numpy as np
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))
from geovac.spectral_triple import FockSpectralTriple
from sympy import Rational

K_over_pi = 43.619934066848226
B = 42.0

for n in [4, 5]:
    print(f"\n{'='*60}")
    print(f"n_max={n}")
    print(f"{'='*60}")

    st = FockSpectralTriple(n_max=n, kappa=Rational(-1,16),
                            j_type='kramers', adjacency_weights='cg')
    D = np.array(st.dirac_operator.tolist(), dtype=float)
    evals = np.linalg.eigvalsh(D)
    evals_sq = evals**2
    dim_H = len(evals)

    print(f"  dim_H = {dim_H}")
    print(f"  |lambda| range: [{min(abs(evals)):.6f}, {max(abs(evals)):.6f}]")

    # Gaussian search for K/pi
    Ls = np.linspace(0.5, 40, 10000)
    best_L, best_val, best_err = None, None, 1e10
    for L in Ls:
        g = float(np.sum(np.exp(-evals_sq / L**2)))
        err = abs(g - K_over_pi)
        if err < best_err:
            best_L, best_val, best_err = L, g, err
    rel = best_err / K_over_pi
    print(f"  Gaussian ~ K/pi at Lambda={best_L:.4f}: value={best_val:.6f}, "
          f"rel_err={rel:.6e} ({rel*100:.4f}%)")

    # Gaussian search for B
    best_L2, best_val2, best_err2 = None, None, 1e10
    for L in Ls:
        g = float(np.sum(np.exp(-evals_sq / L**2)))
        err = abs(g - B)
        if err < best_err2:
            best_L2, best_val2, best_err2 = L, g, err
    rel2 = best_err2 / B
    print(f"  Gaussian ~ B    at Lambda={best_L2:.4f}: value={best_val2:.6f}, "
          f"rel_err={rel2:.6e}")

    # Ratio of Lambdas
    if best_L and best_L2:
        ratio = best_L / best_L2
        print(f"  Lambda(K/pi) / Lambda(B) = {ratio:.6f}")
        print(f"  (K/pi) / B = {K_over_pi/B:.6f}")

    # What is the Lambda relative to CH max?
    ch_max = n - 1 + 1.5  # max CH eigenvalue
    print(f"  CH max eigenvalue: {ch_max}")
    print(f"  Lambda(K/pi) / CH_max = {best_L/ch_max:.6f}")
    print(f"  Lambda(K/pi) / max|lambda_D| = {best_L/max(abs(evals)):.6f}")
