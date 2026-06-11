"""
g2_c3_asymptotic.py -- Multi-term asymptotic expansion of B(n_int) for c3 extraction.

The single-power-law tail fit has systematic error because the effective exponent
drifts (7.5 -> 7.9). The actual asymptotic form is:

    B(n) = a_0/n^{p0} + a_1/n^{p0+1} + a_2/n^{p0+2} + ...

where p0 is the leading power. From dimensional analysis:
- Denominator: lambda^4 * mu_q ~ n^6
- CG asymptotics with j1=1 (fixed) coupling to j2 ~ n: vertex ~ n^{-1}
- Two vertices + probe + m-sums: net numerator ~ n^{-2}
- Total: B(n) ~ n^{-8}, so p0 = 8 is the natural candidate

Strategy: test candidate p0 values (7, 7.5, 8, 8.5, 9) by computing
n^{p0} * B(n) and fitting a polynomial in 1/n. The correct p0 gives a
clean polynomial (no divergence as n -> inf, no spurious curvature).
Then compute the tail exactly via Hurwitz zeta.
"""

import json
import numpy as np
import mpmath

mpmath.mp.dps = 50

# Constants
ALPHA = 7.2973525693e-3
SCHWINGER = ALPHA / (2 * np.pi)
B_HOPF = 42.0
F_HOPF = float(mpmath.pi**2 / 6)
DELTA_HOPF = 1.0 / 40.0
C1 = 0.5
C2_EXACT = (2 - B_HOPF * DELTA_HOPF - F_HOPF * DELTA_HOPF - F_HOPF / B_HOPF) / 5
X_N1 = 4.0 / 25
X3_N1 = X_N1**3


def load_data():
    with open('debug/data/g2_c3_investigation.json') as f:
        data = json.load(f)
    all_B = {int(k): v for k, v in data['all_B'].items()}
    V_mag = data['V_mag']
    return all_B, V_mag


def analyze_asymptotic_structure(all_B):
    """Analyze the asymptotic structure of B(n) for even and odd n separately."""
    print("=" * 70)
    print("  ASYMPTOTIC STRUCTURE ANALYSIS")
    print("=" * 70)

    even_n = sorted(n for n in all_B if n % 2 == 0 and n >= 6)
    odd_n = sorted(n for n in all_B if n % 2 == 1 and n >= 5)

    # Step 1: Consecutive ratio analysis to estimate leading power
    print("\n  --- Consecutive ratio analysis ---")
    print("  B(n)/B(n-2) ~ (1 - 2/n)^p => effective p = -n/2 * ln(ratio)")

    for label, ns in [("Even", even_n), ("Odd", odd_n)]:
        print(f"\n  {label} subsequence:")
        p_effs = []
        for i in range(1, len(ns)):
            n1, n2 = ns[i-1], ns[i]
            r = abs(all_B[n2] / all_B[n1])
            if r > 0 and r < 1:
                p_eff = -np.log(r) / np.log(float(n2) / n1)
                p_effs.append((n2, p_eff))
                if n2 >= 20:
                    print(f"    n={n2:3d}: p_eff = {p_eff:.6f}")

        if len(p_effs) >= 3:
            ns_eff = [x[0] for x in p_effs[-8:]]
            ps_eff = [x[1] for x in p_effs[-8:]]
            print(f"    Last 8 values: mean={np.mean(ps_eff):.4f}, std={np.std(ps_eff):.4f}")

            # Richardson extrapolation on p_eff sequence
            if len(ps_eff) >= 4:
                # Assume p_eff(n) = p_inf + c/n + ..., so Richardson:
                # p_inf ≈ (n2*p2 - n1*p1)/(n2 - n1)
                rich_ps = []
                for j in range(len(ps_eff) - 1):
                    n1_r, n2_r = ns_eff[j], ns_eff[j+1]
                    p1_r, p2_r = ps_eff[j], ps_eff[j+1]
                    p_rich = (n2_r * p2_r - n1_r * p1_r) / (n2_r - n1_r)
                    rich_ps.append(p_rich)
                print(f"    Richardson extrapolated: {[f'{p:.4f}' for p in rich_ps[-4:]]}")

    # Step 2: Test candidate leading powers
    print("\n\n  --- Testing candidate leading powers p0 ---")
    print("  For each p0, compute f(n) = n^p0 * B(n), fit polynomial in 1/n")

    for label, ns in [("Even", even_n), ("Odd", odd_n)]:
        Bs = [all_B[n] for n in ns]

        print(f"\n  {label} subsequence (n = {ns[0]}..{ns[-1]}, {len(ns)} points):")

        best_p0 = None
        best_rms = 1e10

        for p0 in [7.0, 7.5, 7.75, 8.0, 8.25, 8.5, 9.0]:
            # Compute f(n) = n^p0 * B(n)
            fn = [float(n)**p0 * B for n, B in zip(ns, Bs)]
            inv_n = [1.0 / n for n in ns]

            # Use only n >= 10 for fitting
            mask = [n >= 10 for n in ns]
            fn_fit = [f for f, m in zip(fn, mask) if m]
            inv_n_fit = [x for x, m in zip(inv_n, mask) if m]

            if len(fn_fit) < 5:
                continue

            # Fit polynomial of degree 3 in 1/n
            coeffs = np.polyfit(inv_n_fit, fn_fit, 3)
            fn_pred = np.polyval(coeffs, inv_n_fit)
            rms = np.sqrt(np.mean((np.array(fn_fit) - fn_pred)**2))
            rel_rms = rms / np.mean(np.abs(fn_fit))

            # Check if f(n) converges (a0 = coeffs[-1])
            a0 = coeffs[-1]
            a0_sign = "+" if a0 > 0 else "-"

            print(f"    p0={p0:5.2f}: a0={a0:+12.6e}, rms={rel_rms:.3e}, "
                  f"coeffs=[{coeffs[0]:.4e}, {coeffs[1]:.4e}, {coeffs[2]:.4e}, {a0:.4e}]")

            if rel_rms < best_rms:
                best_rms = rel_rms
                best_p0 = p0

        print(f"    Best fit: p0 = {best_p0:.2f} (relative RMS = {best_rms:.3e})")

    # Step 3: Detailed analysis at p0 = 8 (physics prediction)
    print("\n\n  --- Detailed analysis at p0 = 8 ---")

    for label, ns in [("Even", even_n), ("Odd", odd_n)]:
        Bs = [all_B[n] for n in ns]

        # Compute f(n) = n^8 * B(n) for n >= 10
        mask = [n >= 10 for n in ns]
        ns_fit = [n for n, m in zip(ns, mask) if m]
        fn = [float(n)**8 * B for n, B, m in zip(ns, Bs, mask) if m]
        inv_n = [1.0 / n for n in ns_fit]

        print(f"\n  {label} sub-sequence, n^8 * B(n):")
        for n, f in zip(ns_fit, fn):
            print(f"    n={n:3d}: n^8*B(n) = {f:+14.6e}")

        # Fit degree-4 polynomial in 1/n
        coeffs4 = np.polyfit(inv_n, fn, 4)
        fn_pred4 = np.polyval(coeffs4, inv_n)
        rms4 = np.sqrt(np.mean((np.array(fn) - fn_pred4)**2)) / np.mean(np.abs(fn))

        print(f"    Polynomial fit (degree 4): a0={coeffs4[4]:.8e}")
        print(f"      a1={coeffs4[3]:.8e}, a2={coeffs4[2]:.8e}, "
              f"a3={coeffs4[1]:.8e}, a4={coeffs4[0]:.8e}")
        print(f"      relative RMS = {rms4:.3e}")

        # Compare with degree-3 and degree-5
        for deg in [2, 3, 5, 6]:
            c = np.polyfit(inv_n, fn, deg)
            fp = np.polyval(c, inv_n)
            rms = np.sqrt(np.mean((np.array(fn) - fp)**2)) / np.mean(np.abs(fn))
            print(f"      degree {deg}: a0={c[-1]:.8e}, rms={rms:.3e}")

    # Step 4: Also test p0 = 8 on the COMBINED (non-parity-separated) sequence
    print("\n\n  --- Combined sequence test (no parity separation) ---")
    ns_all = sorted(n for n in all_B if n >= 10)
    fn_all = [float(n)**8 * all_B[n] for n in ns_all]
    inv_n_all = [1.0/n for n in ns_all]

    # Even/odd alternation: fit with even/odd indicator
    # f(n) = (a0_e + a1_e/n + ...) if n even, (a0_o + a1_o/n + ...) if n odd
    # Combined: f(n) = a0 + a1/n + ... + (-1)^n * (b0 + b1/n + ...)

    for n, f in zip(ns_all[-10:], fn_all[-10:]):
        print(f"    n={n:3d}: n^8*B(n) = {f:+14.6e}")


def extract_c3_multiterm(all_B, V_mag, p0=8, n_fit_min=10, n_poly_deg=4):
    """Extract c3 using multi-term asymptotic tail correction.

    Fits B(n) = sum_k a_k / n^{p0+k} for even and odd subsequences,
    then computes tail analytically via Hurwitz zeta.
    """
    print("=" * 70)
    print(f"  c3 EXTRACTION WITH MULTI-TERM TAIL (p0={p0})")
    print("=" * 70)

    n_max = max(all_B.keys())
    cum_B = sum(all_B[k] for k in sorted(all_B.keys()) if k <= n_max)

    # Separate even and odd
    even_n = sorted(n for n in all_B if n % 2 == 0 and n >= n_fit_min)
    odd_n = sorted(n for n in all_B if n % 2 == 1 and n >= n_fit_min)

    tails = {}

    for label, ns in [("even", even_n), ("odd", odd_n)]:
        Bs = [all_B[n] for n in ns]
        fn = [float(n)**p0 * B for n, B in zip(ns, Bs)]
        inv_n = [1.0 / n for n in ns]

        coeffs = np.polyfit(inv_n, fn, n_poly_deg)
        # coeffs[0] = highest power, coeffs[-1] = constant (a0)
        # f(n) = a0 + a1/n + a2/n^2 + ... = sum_{k=0}^{deg} coeffs[deg-k] * (1/n)^k
        # So B(n) = sum_k a_k / n^{p0+k} where a_k = coeffs[deg-k]

        a = [coeffs[n_poly_deg - k] for k in range(n_poly_deg + 1)]

        fn_pred = np.polyval(coeffs, inv_n)
        rms = np.sqrt(np.mean((np.array(fn) - fn_pred)**2))

        # Compute tail: sum_{n > n_max, parity} sum_k a_k / n^{p0+k}
        # For even n > n_max: n = n_max+2, n_max+4, ... = 2*(n_max/2 + 1), 2*(n_max/2 + 2), ...
        # sum_{m=M}^inf (2m)^{-(p0+k)} = 2^{-(p0+k)} * zeta(p0+k, M)

        if label == "even":
            M = n_max // 2 + 1  # first even n after n_max is 2*M
        else:
            M_odd = (n_max + 1) // 2  # n_max = 50, first odd > 50 is 51 = 2*26-1
            # odd n: n = 2m-1, m = 1,2,3,...; n > n_max => m > (n_max+1)/2
            M = (n_max + 2) // 2  # first m such that 2m-1 > n_max

        tail = 0.0
        tail_terms = []
        for k in range(n_poly_deg + 1):
            s = p0 + k
            if label == "even":
                # sum_{m=M}^inf (2m)^{-s} = 2^{-s} * zeta(s, M)
                term = a[k] * float(2**(-s) * mpmath.zeta(s, M))
            else:
                # sum_{m=M}^inf (2m-1)^{-s} = sum_{m=M}^inf (2m-1)^{-s}
                # = 2^{-s} * [zeta(s, M-1/2)] (Hurwitz at half-integer)
                term = a[k] * float(2**(-s) * mpmath.zeta(s, mpmath.mpf(M) - mpmath.mpf(1)/2))
            tail += term
            tail_terms.append(term)

        tails[label] = {
            'coeffs': a,
            'tail': tail,
            'tail_terms': tail_terms,
            'rms': rms,
        }

        print(f"\n  {label.upper()} subsequence:")
        print(f"    Asymptotic coefficients a_k (B(n) = sum a_k/n^{{{p0}+k}}):")
        for k, ak in enumerate(a):
            print(f"      a_{k} = {ak:+14.8e}")
        print(f"    Fit RMS = {rms:.3e}")
        print(f"    Tail (n>{n_max}) = {tail:.8e}")
        print(f"    Tail breakdown: {[f'{t:.4e}' for t in tail_terms]}")

    tail_total = tails["even"]["tail"] + tails["odd"]["tail"]
    total_B = cum_B + tail_total

    F2_S = total_B / V_mag / SCHWINGER
    delta = F2_S - 1.0
    c1_term = C1 * X_N1
    c2_term = C2_EXACT * X_N1**2
    residual = delta - c1_term - c2_term
    c3 = residual / X3_N1

    print(f"\n  --- c3 extraction ---")
    print(f"  cum_B (n=0..{n_max}) = {cum_B:.15e}")
    print(f"  tail (n>{n_max})     = {tail_total:.8e}")
    print(f"  total_B              = {total_B:.15e}")
    print(f"  F2/S                 = {F2_S:.15f}")
    print(f"  delta                = {delta:.15f}")
    print(f"  c1*x                 = {c1_term:.15f}")
    print(f"  c2*x^2               = {c2_term:.15f}")
    print(f"  residual             = {residual:.6e}")
    print(f"  c3                   = {c3:.12e}")

    # Uncertainty: compare different polynomial degrees
    print(f"\n  --- Stability under polynomial degree ---")
    c3_vals = []
    for deg in range(2, 7):
        c3_d = extract_c3_quiet(all_B, V_mag, p0, n_fit_min, deg)
        c3_vals.append(c3_d)
        print(f"    deg={deg}: c3 = {c3_d:.12e}")

    c3_mean = np.mean(c3_vals)
    c3_std = np.std(c3_vals)
    print(f"  Mean: {c3_mean:.12e}")
    print(f"  Std:  {c3_std:.4e}")
    print(f"  Significant digits: ~{-np.log10(abs(c3_std/c3_mean)):.1f}" if c3_mean != 0 else "")

    # Stability under p0
    print(f"\n  --- Stability under leading power p0 ---")
    c3_p0s = []
    for p0_test in [7.5, 7.75, 8.0, 8.25, 8.5]:
        c3_p = extract_c3_quiet(all_B, V_mag, p0_test, n_fit_min, 4)
        c3_p0s.append((p0_test, c3_p))
        print(f"    p0={p0_test:.2f}: c3 = {c3_p:.12e}")

    # Stability under fit range
    print(f"\n  --- Stability under fit range (n_fit_min) ---")
    for nfm in [8, 10, 15, 20, 25]:
        c3_f = extract_c3_quiet(all_B, V_mag, p0, nfm, 4)
        print(f"    n_fit_min={nfm:2d}: c3 = {c3_f:.12e}")

    return c3


def extract_c3_quiet(all_B, V_mag, p0, n_fit_min, n_poly_deg):
    """Quiet version of extract_c3_multiterm — returns just c3."""
    n_max = max(all_B.keys())
    cum_B = sum(all_B[k] for k in sorted(all_B.keys()) if k <= n_max)

    tail_total = 0.0
    for label in ["even", "odd"]:
        ns = sorted(n for n in all_B if n % 2 == (0 if label == "even" else 1) and n >= n_fit_min)
        if len(ns) < n_poly_deg + 2:
            continue
        Bs = [all_B[n] for n in ns]
        fn = [float(n)**p0 * B for n, B in zip(ns, Bs)]
        inv_n = [1.0 / n for n in ns]
        coeffs = np.polyfit(inv_n, fn, n_poly_deg)
        a = [coeffs[n_poly_deg - k] for k in range(n_poly_deg + 1)]

        if label == "even":
            M = n_max // 2 + 1
        else:
            M = (n_max + 2) // 2

        for k in range(n_poly_deg + 1):
            s = p0 + k
            if label == "even":
                tail_total += a[k] * float(2**(-s) * mpmath.zeta(s, M))
            else:
                tail_total += a[k] * float(2**(-s) * mpmath.zeta(s, mpmath.mpf(M) - mpmath.mpf(1)/2))

    total_B = cum_B + tail_total
    F2_S = total_B / V_mag / SCHWINGER
    delta = F2_S - 1.0
    residual = delta - C1 * X_N1 - C2_EXACT * X_N1**2
    return residual / X3_N1


def pslq_scan(c3_best):
    """PSLQ scan of c3 against T9-compatible bases."""
    print("\n" + "=" * 70)
    print("  PSLQ IDENTIFICATION OF c3")
    print("=" * 70)

    pi = float(mpmath.pi)
    pi2 = pi**2
    pi4 = pi**4
    pi6 = pi**6
    B, F, D = 42.0, pi2/6, 1.0/40
    c2 = C2_EXACT

    bases = {
        "Minimal": {
            "labels": ["1", "pi^2", "pi^4"],
            "vals": [1.0, pi2, pi4],
        },
        "T9 invariants": {
            "labels": ["1", "pi^2", "pi^4", "B*D", "F*D", "F/B"],
            "vals": [1.0, pi2, pi4, B*D, F*D, F/B],
        },
        "c2-related": {
            "labels": ["1", "pi^2", "pi^4", "c2", "c2^2"],
            "vals": [1.0, pi2, pi4, c2, c2**2],
        },
        "Full T9": {
            "labels": ["1", "pi^2", "pi^4", "B*D", "F*D", "F/B", "D^2",
                       "B*D^2", "F*D^2", "F^2/B", "(BD)^2", "(FD)^2", "BFD^2"],
            "vals": [1.0, pi2, pi4, B*D, F*D, F/B, D**2,
                     B*D**2, F*D**2, F**2/B, (B*D)**2, (F*D)**2, B*F*D**2],
        },
        "Extended pi": {
            "labels": ["1", "pi^2", "pi^4", "pi^6"],
            "vals": [1.0, pi2, pi4, pi6],
        },
        "Odd zeta check": {
            "labels": ["1", "pi^2", "zeta(3)", "pi^4"],
            "vals": [1.0, pi2, float(mpmath.zeta(3)), pi4],
        },
    }

    # Scale c3 for better PSLQ matching
    scalings = {
        "c3": c3_best,
        "c3*5": c3_best * 5,
        "c3*25200": c3_best * 25200,
        "c3*lambda^6": c3_best * (5.0/2)**6,
        "c3*25200*lambda^6": c3_best * 25200 * (5.0/2)**6,
    }

    for scale_name, c3_scaled in scalings.items():
        print(f"\n  --- {scale_name} = {c3_scaled:.12e} ---")
        for basis_name, basis in bases.items():
            labels = basis["labels"]
            vals = basis["vals"]
            vec = [mpmath.mpf(str(c3_scaled))] + [mpmath.mpf(str(v)) for v in vals]
            vec = mpmath.matrix(vec)
            try:
                for tol_exp in [-4, -6, -8]:
                    rel = mpmath.pslq(vec, maxcoeff=10000, tol=mpmath.mpf(10)**tol_exp)
                    if rel is not None:
                        check = sum(float(rel[j]) * float(vec[j]) for j in range(len(rel)))
                        max_coeff = max(abs(int(r)) for r in rel)
                        if max_coeff < 5000:
                            print(f"    {basis_name} (tol=1e{tol_exp}): {list(rel)}")
                            print(f"      labels: {scale_name}, {', '.join(labels)}")
                            print(f"      check: {check:.3e}, max_coeff: {max_coeff}")
                            if int(rel[0]) != 0:
                                c3_recon = -sum(int(rel[j]) * float(vec[j])
                                               for j in range(1, len(rel))) / int(rel[0])
                                print(f"      c3_reconstructed = {c3_recon:.12e}")
                            break
                else:
                    print(f"    {basis_name}: no relation (all tol)")
            except Exception as e:
                print(f"    {basis_name}: error: {e}")


def main():
    all_B, V_mag = load_data()

    # Part 1: Analyze asymptotic structure
    analyze_asymptotic_structure(all_B)

    # Part 2: Extract c3 with multi-term tail
    c3_best = extract_c3_multiterm(all_B, V_mag, p0=8.0)

    # Part 3: PSLQ scan
    pslq_scan(c3_best)


if __name__ == "__main__":
    main()
