"""
g2_c3_richardson.py -- Richardson extrapolation for c3.

The bottleneck in c3 extraction is the tail correction: B(n) for n > N_max.
Instead of fitting a parametric model B(n) ~ C*n^-p, we use non-parametric
sequence acceleration on the c3(N) sequence, where c3(N) is the c3 estimate
from partial sums through n_int = N.

Method: Richardson extrapolation + Aitken delta-squared acceleration.
If c3(N) = c3 + a/N^q + ..., Richardson eliminates the leading error term.
Multi-level Richardson (Romberg-like) eliminates successive error terms.

Also tries: Wynn epsilon algorithm (the most powerful general-purpose
sequence accelerator for monotone sequences).
"""

import json
import numpy as np
import mpmath
mpmath.mp.dps = 30

# ---------------------------------------------------------------------------
# Constants (exact where possible)
# ---------------------------------------------------------------------------
PI = float(mpmath.pi)
PI2 = PI**2
PI4 = PI**4

ALPHA = 7.2973525693e-3
SCHWINGER = ALPHA / (2 * float(mpmath.pi))

B_VAL = 42.0
F_VAL = PI2 / 6
DELTA = 1.0 / 40
C1 = 0.5
C2_EXACT = (2 - B_VAL * DELTA - F_VAL * DELTA - F_VAL / B_VAL) / 5
X_N1 = 4.0 / 25  # 1/lambda^2 at n_ext=1, lambda=5/2
X3_N1 = X_N1**3


def load_data():
    with open('debug/data/g2_c3_investigation.json') as f:
        data = json.load(f)
    all_B = {int(k): v for k, v in data['all_B'].items()}
    V_mag = data['V_mag']
    return all_B, V_mag


def c3_from_partial_sum(all_B, V_mag, N):
    """c3 estimate using partial sum through n_int = N.

    Note: this is the PARTIAL sum without tail correction,
    so it converges from below (in magnitude) to the true c3.
    """
    cum_B = sum(all_B[k] for k in sorted(all_B.keys()) if k <= N)
    F2_S = cum_B / V_mag / SCHWINGER
    delta = F2_S - 1.0
    residual = delta - C1 * X_N1 - C2_EXACT * X_N1**2
    return residual / X3_N1


def wynn_epsilon(seq):
    """Wynn's epsilon algorithm for sequence acceleration.

    Given a sequence s_0, s_1, ..., s_n, returns accelerated estimates.
    The epsilon table entries at even columns give the accelerated values.
    """
    n = len(seq)
    if n < 3:
        return seq[-1] if seq else None

    # Build epsilon table
    eps = np.zeros((n, n + 1))
    eps[:, 0] = 0.0  # epsilon_{-1}
    eps[:, 1] = np.array(seq, dtype=float)  # epsilon_0 = s_n

    for j in range(2, n + 1):
        for i in range(n - j + 1):
            diff = eps[i+1, j-1] - eps[i, j-1]
            if abs(diff) < 1e-30:
                eps[i, j] = 1e30  # avoid division by zero
            else:
                eps[i, j] = eps[i+1, j-2] + 1.0 / diff

    # Even columns (j=2,4,6,...) give accelerated values
    # Return the best estimate (highest even column, last row)
    best = seq[-1]
    best_col = -1
    for j in range(2, n + 1, 2):
        idx = n - j
        if idx >= 0:
            val = eps[idx, j]
            if abs(val) < 1e20:  # not overflow
                best = val
                best_col = j

    return best, best_col


def richardson_extrapolation(seq_N, seq_vals, p):
    """Richardson extrapolation assuming error ~ 1/N^p.

    Given (N_i, c3(N_i)), eliminates the leading O(1/N^p) error.
    """
    n = len(seq_N)
    if n < 2:
        return seq_vals[-1]

    new_vals = []
    for i in range(n - 1):
        N1, N2 = seq_N[i], seq_N[i+1]
        v1, v2 = seq_vals[i], seq_vals[i+1]
        ratio = (N2 / N1)**p
        new_vals.append((v2 * ratio - v1) / (ratio - 1))

    return new_vals


def multi_level_richardson(seq_N, seq_vals, p_start, p_step=1):
    """Multi-level Richardson: eliminate 1/N^p, 1/N^{p+1}, etc."""
    N_list = list(seq_N)
    vals = list(seq_vals)
    results = [vals[-1]]
    p = p_start

    while len(vals) > 1:
        vals = richardson_extrapolation(N_list, vals, p)
        N_list = N_list[1:]
        results.append(vals[-1] if vals else results[-1])
        p += p_step
        if len(vals) <= 1:
            break

    return results


def main():
    all_B, V_mag = load_data()
    N_max = max(all_B.keys())

    print("=" * 70)
    print("  RICHARDSON EXTRAPOLATION FOR c3")
    print("=" * 70)

    # Build c3(N) sequence for N = 10, 12, 14, ..., 50
    N_vals = list(range(10, N_max + 1, 2))
    c3_vals = [c3_from_partial_sum(all_B, V_mag, N) for N in N_vals]

    print("\n  Raw c3(N) sequence:")
    for N, c3 in zip(N_vals, c3_vals):
        print(f"    N={N:3d}: c3 = {c3:.12e}")

    # Estimate convergence rate from consecutive differences
    print("\n  Consecutive differences (c3(N+2) - c3(N)):")
    diffs = []
    for i in range(len(c3_vals) - 1):
        d = c3_vals[i+1] - c3_vals[i]
        ratio = d / (c3_vals[i] - c3_vals[i-1]) if i > 0 and abs(c3_vals[i] - c3_vals[i-1]) > 1e-30 else 0
        diffs.append(d)
        if i > 0:
            print(f"    N={N_vals[i]:3d}->N={N_vals[i+1]:3d}: d={d:.6e}, ratio={ratio:.4f}")

    # Estimate effective p from ratio: if error ~ 1/N^p, ratio ~ (N/(N+2))^p
    if len(diffs) >= 3:
        ratios = []
        for i in range(1, len(diffs)):
            if abs(diffs[i-1]) > 1e-30:
                r = abs(diffs[i] / diffs[i-1])
                N1 = N_vals[i]
                N2 = N_vals[i+1]
                if r > 0 and r < 1:
                    p_est = -np.log(r) / np.log(N2 / N1)
                    ratios.append((N1, p_est))

        print("\n  Estimated convergence order p:")
        for N, p in ratios[-8:]:
            print(f"    N={N}: p = {p:.3f}")

    # Wynn epsilon on the full sequence
    print("\n  Wynn epsilon acceleration:")

    # Try with different sub-sequences
    for start in [0, 5, 10]:
        sub_vals = c3_vals[start:]
        if len(sub_vals) >= 5:
            result, col = wynn_epsilon(sub_vals)
            print(f"    start={N_vals[start]:3d}, n_pts={len(sub_vals)}: "
                  f"c3 = {result:.12e} (col={col})")

    # Richardson extrapolation with various assumed orders
    print("\n  Multi-level Richardson (p_start, result):")
    N_sub = N_vals[5:]  # from N=20 onwards
    c3_sub = c3_vals[5:]

    for p_start in [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]:
        results = multi_level_richardson(N_sub, c3_sub, p_start)
        final = results[-1]
        print(f"    p={p_start:.1f}: c3 = {final:.12e} (levels={len(results)})")

    # Also try with fractional p increment
    for p_start in [6.0, 6.5, 7.0]:
        results = multi_level_richardson(N_sub, c3_sub, p_start, p_step=2)
        final = results[-1]
        print(f"    p={p_start:.1f} step=2: c3 = {final:.12e}")

    # Best estimate: use Wynn on the most converged sub-sequence
    print("\n" + "=" * 70)
    print("  BEST ESTIMATES")
    print("=" * 70)

    # Wynn on N=20..50
    sub20 = c3_vals[5:]
    wynn_result, wynn_col = wynn_epsilon(sub20)
    print(f"\n  Wynn (N=20..50, {len(sub20)} pts): c3 = {wynn_result:.14e}")

    # Wynn on N=30..50
    sub30 = c3_vals[10:]
    wynn_result30, wynn_col30 = wynn_epsilon(sub30)
    print(f"  Wynn (N=30..50, {len(sub30)} pts): c3 = {wynn_result30:.14e}")

    # Richardson p=6.5 (near estimated convergence order)
    rich_results = multi_level_richardson(N_sub, c3_sub, 6.5)
    print(f"  Richardson (p=6.5, N=20..50):  c3 = {rich_results[-1]:.14e}")

    # Compare to known c3 from tail-corrected extraction
    c3_ref = -5.946063443372761e-07
    print(f"\n  Reference (tail-corrected):     c3 = {c3_ref:.14e}")

    # Spread of estimates
    estimates = [wynn_result, wynn_result30, rich_results[-1]]
    est_mean = np.mean(estimates)
    est_std = np.std(estimates)
    print(f"\n  Mean of best estimates: c3 = {est_mean:.14e}")
    print(f"  Std of best estimates:       +/- {est_std:.4e}")
    print(f"  Significant digits:          ~{-np.log10(abs(est_std/est_mean)):.1f}" if est_mean != 0 else "")

    # PSLQ attempt with best estimate
    print("\n" + "=" * 70)
    print("  PSLQ IDENTIFICATION")
    print("=" * 70)

    c3_best = est_mean

    # Scale c3 for better PSLQ: multiply by known denominators
    # c3 * lambda^6 where lambda = 5/2
    c3_lam6 = c3_best * (5.0/2)**6
    print(f"\n  c3 * lambda^6 = {c3_lam6:.14e}")

    # c3 * 25200 (LCD of c2 denominators)
    c3_scaled = c3_best * 25200
    print(f"  c3 * 25200    = {c3_scaled:.14e}")

    # c3 * 25200 * (5/2)^6
    c3_full = c3_best * 25200 * (5.0/2)**6
    print(f"  c3 * 25200 * lambda^6 = {c3_full:.14e}")

    # Basis elements for PSLQ
    bases = {
        "Minimal T9": {
            "labels": ["1", "pi^2", "pi^4"],
            "vals": [1.0, PI2, PI4],
        },
        "Paper 2 invariants": {
            "labels": ["1", "pi^2", "pi^4", "B*D", "F*D", "F/B",
                       "B*D^2", "F*D^2", "F^2/B", "F^2*D", "B^2*D^2"],
            "vals": [1.0, PI2, PI4, B_VAL*DELTA, F_VAL*DELTA, F_VAL/B_VAL,
                     B_VAL*DELTA**2, F_VAL*DELTA**2, F_VAL**2/B_VAL,
                     F_VAL**2*DELTA, B_VAL**2*DELTA**2],
        },
        "c2-based": {
            "labels": ["1", "c2", "c2^2", "c1*c2", "pi^2", "pi^4", "pi^6"],
            "vals": [1.0, C2_EXACT, C2_EXACT**2, 0.5*C2_EXACT, PI2, PI4, PI**6],
        },
    }

    # Test each scaled version against each basis
    for c3_label, c3_test in [("c3", c3_best),
                               ("c3*lam6", c3_lam6),
                               ("c3*25200", c3_scaled)]:
        print(f"\n  --- Testing {c3_label} = {c3_test:.12e} ---")
        for basis_name, basis in bases.items():
            labels = basis["labels"]
            vals = basis["vals"]
            vec = [mpmath.mpf(str(c3_test))] + [mpmath.mpf(str(v)) for v in vals]
            vec = mpmath.matrix(vec)
            try:
                rel = mpmath.pslq(vec, maxcoeff=10000, tol=mpmath.mpf('1e-6'))
                if rel is not None:
                    check = sum(float(rel[j]) * float(vec[j]) for j in range(len(rel)))
                    print(f"    {basis_name}: {list(rel)}")
                    print(f"      labels: c3_scaled, {', '.join(labels)}")
                    print(f"      check: {check:.3e}")
                    if rel[0] != 0:
                        c3_reconstructed = -sum(int(rel[j]) * float(vec[j])
                                                for j in range(1, len(rel))) / int(rel[0])
                        print(f"      c3_reconstructed = {c3_reconstructed:.12e}")
                else:
                    print(f"    {basis_name}: no relation (tol=1e-6)")
            except Exception as e:
                print(f"    {basis_name}: error: {e}")


if __name__ == "__main__":
    main()
