"""Sprint MR-B: M2 transcendental signature in the spectral action convergence rate.

Tests two distinct convergence rates on the GeoVac S^3 Dirac spectral triple:

(A) Truncated -> continuum tail rate
    r(n_max, t) = K_continuum(t) - K_{n_max}(t) = sum_{n>n_max} g_n exp(-(n+3/2)^2 t)
    For each fixed t (i.e. fixed Lambda^2 = 1/t), fit r ~ C * n^beta * exp(-alpha*n^2)
    PSLQ alpha and C/n_max^beta against M2 ring sqrt(pi)*Q + pi^2 * Q.

(B) Modular residual rate (the cleaner M2 test)
    K(t) = (sqrt(pi)/2) * t^{-3/2} - (sqrt(pi)/4) * t^{-1/2} + epsilon(t)
    For small t, epsilon(t) ~ C * exp(-A/t) where A is predicted to be exactly pi^2
    (Jacobi theta modular transformation gives the exp(-pi^2/t) suppression).
    Fit log|epsilon(t)| ~ -A/t + B*log(t) + C
    PSLQ A against pi^2 and other M2 candidates.

References:
- CLAUDE.md sec 2 "Spectral action on S^3 Dirac" (v2.26.1, May 2026)
- Paper 35 / debug/spectral_action_sd_exactness.json
- Paper 18 sec III.7 master Mellin engine
- L2 sprint memo 2026-05-04 (4/pi M1 signature established)
"""
import json
from pathlib import Path
import numpy as np
from mpmath import (
    mp, mpf, exp as mpexp, log as mplog, fsum, nstr,
    pi as mp_pi, sqrt as mpsqrt, pslq, mpc
)

mp.dps = 100  # 100 decimal places for PSLQ

OUTDIR = Path(__file__).resolve().parent / "data"
OUTDIR.mkdir(parents=True, exist_ok=True)
OUTPATH = OUTDIR / "mr_b_spectral_action_rate.json"

# =============================================================================
# Building blocks
# =============================================================================

def heat_kernel_truncated(t, n_max):
    """K_{n_max}(t) = sum_{n=0}^{n_max} 2(n+1)(n+2) exp(-(n+3/2)^2 t)."""
    t = mpf(t)
    terms = []
    for n in range(n_max + 1):
        g = mpf(2 * (n + 1) * (n + 2))
        lam_sq = (mpf(n) + mpf('1.5'))**2
        terms.append(g * mpexp(-lam_sq * t))
    return fsum(terms)

def heat_kernel_continuum(t, n_terms_safe=400):
    """K(t) summed to a safe n_terms; for our t values (>= 0.05) this is exact to mp.dps."""
    t = mpf(t)
    terms = []
    for n in range(n_terms_safe):
        g = mpf(2 * (n + 1) * (n + 2))
        lam_sq = (mpf(n) + mpf('1.5'))**2
        terms.append(g * mpexp(-lam_sq * t))
    return fsum(terms)

def two_term_sd(t):
    """Two-term SD asymptotic: (sqrt(pi)/2)*t^{-3/2} - (sqrt(pi)/4)*t^{-1/2}."""
    t = mpf(t)
    sp = mpsqrt(mp_pi)
    return sp / 2 / t**mpf('1.5') - sp / 4 / t**mpf('0.5')

def modular_residual(t, n_terms_safe=400):
    """epsilon(t) = K(t) - 2-term SD."""
    K = heat_kernel_continuum(t, n_terms_safe)
    return K - two_term_sd(t)

# =============================================================================
# PSLQ helpers
# =============================================================================

def m2_basis():
    """Build the M2 transcendental basis: rationals, sqrt(pi), pi, pi^2, products."""
    return [
        ('1', mpf(1)),
        ('sqrt(pi)', mpsqrt(mp_pi)),
        ('pi', mp_pi),
        ('pi^2', mp_pi**2),
        ('pi^(3/2)', mp_pi**mpf('1.5')),
        ('1/pi', 1/mp_pi),
        ('1/pi^2', 1/mp_pi**2),
        ('1/sqrt(pi)', 1/mpsqrt(mp_pi)),
        ('pi^3', mp_pi**3),
        ('pi^4', mp_pi**4),
    ]

def try_pslq(target, max_coeff=10**6, verbose=False):
    """Try to find an integer relation between target and the M2 basis."""
    basis = m2_basis()
    vec = [target] + [val for _, val in basis]
    names = ['target'] + [name for name, _ in basis]
    try:
        rel = pslq(vec, maxcoeff=max_coeff, tol=mpf(10)**(-50))
        if rel is None:
            return None
        # Build a human-readable form: target = -sum(rel[i+1]*basis[i])/rel[0]
        if rel[0] == 0:
            return None
        c0 = rel[0]
        contributions = []
        for name, coef in zip(names[1:], rel[1:]):
            if coef != 0:
                contributions.append(f"({-coef}/{c0}) * {name}")
        return {
            'relation': rel,
            'names': names,
            'reconstruction': ' + '.join(contributions) if contributions else '0',
            'c0': int(c0),
            'coeffs': [int(c) for c in rel],
        }
    except Exception as e:
        return {'error': str(e)}

# =============================================================================
# PART A: truncated -> continuum tail rate
# =============================================================================

def part_a_results():
    """For several t values (i.e. several Lambda), compute r(n_max, t) for n_max in 1..30."""
    out = {}

    # Test t values:
    # t = 1            -> Lambda^2 = 1
    # t = 1/pi         -> Lambda^2 = pi  (a clean M2 ring choice)
    # t = 1/Lambda_inf^2 -> Lambda = Lambda_inf (the alpha-conjecture cubic root)
    Lambda_inf_sq = mpf('3.7102454679060528505')**2
    t_values = [
        ('t=1', mpf(1)),
        ('t=1/pi', 1 / mp_pi),
        ('t=1/Lambda_inf^2', 1 / Lambda_inf_sq),
        ('t=4/pi^2', 4 / mp_pi**2),
        ('t=1/4', mpf('0.25')),
    ]

    n_max_list = list(range(1, 31))
    for label, t in t_values:
        K_cont = heat_kernel_continuum(t, n_terms_safe=400)
        residuals = []
        log_residuals = []
        for n_max in n_max_list:
            K_trunc = heat_kernel_truncated(t, n_max)
            r = K_cont - K_trunc
            residuals.append(float(r))
            if r > 0:
                log_residuals.append(float(mplog(r)))
            else:
                log_residuals.append(None)

        # Asymptotic prediction: r(n_max) ~ g_{n_max+1} * exp(-(n_max+5/2)^2 * t)
        #                                = 2(n_max+2)(n_max+3) * exp(-(n_max+5/2)^2 * t)
        # log r(n_max) ~ log(2(n_max+2)(n_max+3)) - (n_max+5/2)^2 * t
        #
        # Fit: log r = -t * (n + 5/2)^2 + log[2(n+2)(n+3)] + correction
        # The leading -t*(n+5/2)^2 is exact (single-term tail dominates).
        # Extract t from the slope of log r vs (n+5/2)^2.

        # Use n_max in 5..25 for the fit (avoid early-n curvature, late-n underflow)
        valid = [(n, lr) for n, lr in zip(n_max_list, log_residuals)
                 if lr is not None and -700 < lr < 0]
        valid_fit = [(n, lr) for n, lr in valid if 5 <= n <= 25]

        if len(valid_fit) >= 5:
            ns = np.array([n for n, _ in valid_fit])
            lrs = np.array([lr for _, lr in valid_fit])
            # Subtract the predicted log degeneracy 2(n+2)(n+3) and fit log[r/g] vs (n+5/2)^2
            log_g = np.log(2 * (ns + 2) * (ns + 3))
            x = (ns + 2.5)**2
            y = lrs - log_g  # should be ~ -t * x + small constant
            # Linear fit
            A = np.vstack([x, np.ones_like(x)]).T
            slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
            alpha_fit = -slope
            const_fit = float(np.exp(intercept))

            # Verify alpha_fit ≈ t
            t_float = float(t)
            rel_err = (alpha_fit - t_float) / t_float
        else:
            alpha_fit = None
            const_fit = None
            rel_err = None

        # PSLQ check on alpha_fit (the rate constant)
        pslq_alpha = None
        if alpha_fit is not None and len(valid_fit) >= 8:
            # High-precision rate extraction: r(n)/r(n-1) at large n approaches
            # exp(-t * [2n + 5]) * (2(n+2)(n+3)) / (2(n+1)(n+2))
            # log[r(n)/r(n-1)] ~ -t * (2n+5) + log[(n+3)/(n+1)]
            # We want to extract t very precisely.
            # Use the highest-quality data points.
            best_two = sorted([(n, lr) for n, lr in valid if 15 <= n <= 25], key=lambda x: -x[1])
            if len(best_two) >= 2:
                # Use successive differences with high precision arithmetic
                diffs = []
                for i in range(1, len(valid)):
                    n0, lr0 = valid[i-1]
                    n1, lr1 = valid[i]
                    if n1 - n0 == 1 and lr1 < lr0:
                        # log[r(n1)/r(n0)] = -t*(2*n0+5+1) + log[(n0+3)/(n0+1)]
                        # Wait: g(n+1)/g(n) = (n+2)(n+3) / ((n+1)(n+2)) = (n+3)/(n+1)
                        # exp diff: (n+5/2+1)^2 - (n+5/2)^2 = 2(n+5/2) + 1 = 2n+6
                        # So log r(n+1) - log r(n) = -t*(2n+6) + log[(n+3)/(n+1)]
                        # => -t = [(lr1 - lr0) - log((n0+3)/(n0+1))] / (2*n0+6)
                        n0_mp = mpf(n0)
                        d = mpf(lr1 - lr0) - mplog(mpf(n0+3) / mpf(n0+1))
                        t_from_pair = -d / (2*n0_mp + 6)
                        diffs.append((n0, t_from_pair))
                if len(diffs) >= 3:
                    # Take the average of the latest 5 (most asymptotic)
                    latest = diffs[-min(8, len(diffs)):]
                    t_estimates = [v for _, v in latest]
                    t_extracted = fsum(t_estimates) / len(t_estimates)
                    pslq_alpha = try_pslq(t_extracted)
                    pslq_alpha_t_extracted = nstr(t_extracted, 30)
                    # And the final-pair high-precision estimate
                    if len(diffs) >= 1:
                        final = diffs[-1][1]
                        pslq_alpha_final = try_pslq(final)
                else:
                    pslq_alpha_t_extracted = None
                    pslq_alpha_final = None
            else:
                pslq_alpha_t_extracted = None
                pslq_alpha_final = None
        else:
            pslq_alpha_t_extracted = None
            pslq_alpha_final = None

        out[label] = {
            't_value_50dps': nstr(t, 50),
            't_value_float': float(t),
            'K_continuum_30dps': nstr(K_cont, 30),
            'n_max_list': n_max_list,
            'residuals_log_float': [float(lr) if lr is not None else None
                                    for lr in log_residuals],
            'residuals_float': residuals,
            'fit_alpha_via_lstsq': alpha_fit,
            'fit_alpha_rel_err_vs_t': rel_err,
            'fit_const': const_fit,
            't_extracted_high_precision_50dps': pslq_alpha_t_extracted,
            'pslq_t_extracted': pslq_alpha,
            'pslq_t_final_pair': pslq_alpha_final if len(valid_fit) >= 8 else None,
            'asymptotic_prediction': 'r(n_max) ~ 2(n_max+2)(n_max+3) * exp(-t*(n_max+5/2)^2)',
            'note': 'rate constant alpha = t exactly (single-term tail dominates)',
        }

    return out

# =============================================================================
# PART B: modular residual epsilon(t) = K(t) - 2-term SD ~ C * exp(-pi^2/t)
# =============================================================================

def part_b_results():
    """Extract the exponent A in epsilon(t) ~ C * exp(-A/t) and PSLQ against pi^2."""
    out = {}

    # The continuum heat kernel via direct sum is unstable for very small t
    # because we need huge n_terms. Use n_terms_safe scaled appropriately.
    # For t = 0.1, dominant n ~ 1/(0.1) ~ 10, need n_terms >> 10. Use 500.
    # For t = 0.05, dominant n ~ 20, need n_terms ~ 200. Use 1000.
    # For t = 0.02, dominant n ~ 50, need n_terms ~ 500. Use 2000.
    # For t = 0.01, dominant n ~ 100, need n_terms ~ 1000. Use 4000.
    # The exponent of suppression is exp(-(n+3/2)^2 * t), so n_terms = 5/sqrt(t) suffices.

    t_values = [mpf(s) for s in
                ['0.5', '0.4', '0.3', '0.25', '0.2', '0.15', '0.12',
                 '0.1', '0.08', '0.06', '0.05', '0.04', '0.03']]
    samples = []
    for t in t_values:
        # Choose n_terms to push to 100-digit precision: need (n+3/2)^2 * t > 250 (since e^-250 < 10^-108)
        # i.e. n_terms > sqrt(250/t)
        n_terms = max(400, int(2 * float(mpsqrt(mpf(250) / t))))
        K = heat_kernel_continuum(t, n_terms_safe=n_terms)
        eps = K - two_term_sd(t)
        log_eps = mplog(abs(eps))
        samples.append({
            't': nstr(t, 20),
            't_float': float(t),
            'K_30dps': nstr(K, 30),
            'two_term_30dps': nstr(two_term_sd(t), 30),
            'epsilon_30dps': nstr(eps, 30),
            'epsilon_sign': int(1 if eps > 0 else -1),
            'log_abs_epsilon': nstr(log_eps, 30),
            'log_abs_epsilon_float': float(log_eps),
            'one_over_t_float': float(1 / t),
            'n_terms_used': n_terms,
        })

    # Fit log|eps(t)| ~ -A/t + B*log(t) + C
    # Use small-t data points for the asymptotic regime
    fit_data = [(float(s['t']), float(s['log_abs_epsilon']))
                for s in samples if float(s['t']) <= 0.1]

    if len(fit_data) >= 4:
        ts_f = np.array([d[0] for d in fit_data])
        ys_f = np.array([d[1] for d in fit_data])
        # Design: y = -A * (1/t) + B * log(t) + C
        X = np.vstack([1/ts_f, np.log(ts_f), np.ones_like(ts_f)]).T
        coef, *_ = np.linalg.lstsq(X, ys_f, rcond=None)
        A_fit, B_fit, C_fit = coef
        A_fit_neg = -A_fit  # because design has +1/t * (-A); we want A = -coef[0]
        # Predicted: A = pi^2 (from Jacobi modular transformation)
        pi2 = float(mp_pi)**2
        rel_err_pi2 = (A_fit_neg - pi2) / pi2
    else:
        A_fit, B_fit, C_fit, A_fit_neg, rel_err_pi2 = (None,)*5

    # High-precision extraction of A using two pairs of points
    # log|eps(t1)| - log|eps(t2)| ~ -A*(1/t1 - 1/t2) + B*log(t1/t2)
    # If we use ratios (t1, t2) and (t3, t4) we can solve for A and B.
    # Simpler: fit at very small t where the log(t) term is sub-leading.
    pi2_mp = mp_pi**2
    A_estimates = []
    for i in range(len(samples) - 1):
        s1, s2 = samples[i], samples[i+1]
        if float(s1['t']) > 0.15:
            continue
        t1, t2 = mpf(s1['t']), mpf(s2['t'])
        le1 = mpf(s1['log_abs_epsilon'])
        le2 = mpf(s2['log_abs_epsilon'])
        # log r(t2) - log r(t1) = -A*(1/t2 - 1/t1) + B*log(t2/t1)
        # Use 3-point method to solve for A, B simultaneously when we have 3 points
        # Or just take simple slope and accept B-leakage
        slope = (le2 - le1) / (1/t2 - 1/t1)  # this is approximately -A
        A_est = -slope
        A_estimates.append({
            't1': nstr(t1, 10), 't2': nstr(t2, 10),
            'A_estimate': nstr(A_est, 30),
            'A_estimate_float': float(A_est),
            'rel_err_vs_pi^2': float((A_est - pi2_mp) / pi2_mp),
        })

    # Best (smallest-t) estimate
    A_best = mpf(A_estimates[-1]['A_estimate']) if A_estimates else None

    # Three-point extraction: solve for A and B simultaneously
    A_3pt = None
    B_3pt = None
    if len(samples) >= 3:
        # Use the three smallest-t samples
        small_t_samples = sorted(samples, key=lambda s: float(s['t']))[:6]
        if len(small_t_samples) >= 3:
            # log eps = -A/t + B*log(t) + C. Use any 3 points.
            for triple_idx in [(0,1,2), (0,1,3), (0,2,4), (1,2,3)]:
                if max(triple_idx) >= len(small_t_samples):
                    continue
                pts = [small_t_samples[i] for i in triple_idx]
                t1, t2, t3 = [mpf(p['t']) for p in pts]
                le1, le2, le3 = [mpf(p['log_abs_epsilon']) for p in pts]
                # linear system:
                # le_i = -A/t_i + B*log(t_i) + C
                M = mp.matrix(3, 3)
                rhs = mp.matrix(3, 1)
                for i, (ti, lei) in enumerate(zip([t1,t2,t3], [le1,le2,le3])):
                    M[i, 0] = -1/ti
                    M[i, 1] = mplog(ti)
                    M[i, 2] = 1
                    rhs[i, 0] = lei
                try:
                    sol = mp.lu_solve(M, rhs)
                    A_cand, B_cand, _ = sol[0,0], sol[1,0], sol[2,0]
                    if A_cand > 0:
                        A_3pt = A_cand
                        B_3pt = B_cand
                        used_triple = triple_idx
                        break
                except Exception:
                    continue

    # PSLQ: A = pi^2 ?
    pslq_A = None
    if A_best is not None:
        pslq_A = try_pslq(A_best)
    pslq_A_3pt = None
    if A_3pt is not None:
        pslq_A_3pt = try_pslq(A_3pt)

    # Subleading fit: refit with B fixed at theoretical values to get clean A
    # After Jacobi inversion, the modular tail of K(t) is:
    # K(t) = (sqrt(pi)/2) t^{-3/2} - (sqrt(pi)/4) t^{-1/2} + ... + (modular terms)
    # The modular term involves theta_2 / theta_3 corrections at zero.
    # For the unit S^3 Dirac, the spectrum is |lambda_n| = n + 3/2, the half-integer
    # shift means the spectral sum involves (1/2)-shifted theta. The modular dual
    # is governed by theta(0; -1/tau) ~ exp(-pi^2/(t)) at small t, with prefactor
    # picking up powers of t from the modular weight.

    out['samples'] = samples
    out['fit_lstsq'] = {
        'A_pred_pi^2': float(mp_pi**2),
        'A_fit': A_fit_neg,
        'B_fit': B_fit,
        'C_fit': C_fit,
        'rel_err_pi^2': rel_err_pi2,
    }
    out['pairwise_A_estimates'] = A_estimates
    if A_best is not None:
        out['A_best_estimate_50dps'] = nstr(A_best, 50)
        out['A_best_estimate_float'] = float(A_best)
        out['A_best_rel_err_pi^2'] = float((A_best - mp_pi**2) / mp_pi**2)
        out['pslq_A_best'] = pslq_A
    if A_3pt is not None:
        out['A_3point_estimate_50dps'] = nstr(A_3pt, 50)
        out['B_3point_estimate'] = nstr(B_3pt, 30)
        out['A_3point_float'] = float(A_3pt)
        out['A_3point_rel_err_pi^2'] = float((A_3pt - mp_pi**2) / mp_pi**2)
        out['pslq_A_3point'] = pslq_A_3pt

    return out

# =============================================================================
# Main
# =============================================================================

def main():
    print("="*72)
    print("Sprint MR-B: spectral action convergence rate")
    print("="*72)
    print()

    print("PART A: truncated -> continuum tail rate at fixed Lambda")
    print("-"*72)
    a = part_a_results()
    for label, result in a.items():
        print(f"\n  {label}: t = {result['t_value_float']}")
        if result.get('fit_alpha_via_lstsq') is not None:
            print(f"    fit alpha (lstsq)         = {result['fit_alpha_via_lstsq']:.10f}")
            print(f"    rel err vs t              = {result['fit_alpha_rel_err_vs_t']:.2e}")
        if result.get('t_extracted_high_precision_50dps'):
            print(f"    t_extracted (high prec)   = {result['t_extracted_high_precision_50dps']}")
        if result.get('pslq_t_extracted'):
            ps = result['pslq_t_extracted']
            if isinstance(ps, dict) and 'reconstruction' in ps:
                print(f"    PSLQ:                       {ps['reconstruction']}")

    print()
    print("PART B: modular residual epsilon(t) = K(t) - 2-term SD")
    print("-"*72)
    b = part_b_results()
    print()
    print("  Samples (t, log|epsilon|):")
    for s in b['samples']:
        print(f"    t = {s['t_float']:.4f}    log|eps| = {s['log_abs_epsilon_float']:+10.4f}    "
              f"sign = {s['epsilon_sign']:+d}    n_terms_used = {s['n_terms_used']}")

    print()
    print("  Pairwise A estimates (smallest t pairs are best):")
    for est in b['pairwise_A_estimates']:
        print(f"    t1={est['t1']:>10}  t2={est['t2']:>10}  "
              f"A = {est['A_estimate_float']:.10f}    "
              f"rel_err vs pi^2 = {est['rel_err_vs_pi^2']:+.4e}")

    if 'A_best_estimate_50dps' in b:
        print()
        print(f"  Best (smallest-t) A estimate (50 dps):")
        print(f"    A = {b['A_best_estimate_50dps']}")
        print(f"    pi^2 = {nstr(mp_pi**2, 50)}")
        print(f"    rel err = {b['A_best_rel_err_pi^2']:+.4e}")

    if 'A_3point_estimate_50dps' in b:
        print()
        print(f"  3-point fit (corrects for B*log(t) subleading):")
        print(f"    A = {b['A_3point_estimate_50dps']}")
        print(f"    rel err vs pi^2 = {b['A_3point_rel_err_pi^2']:+.4e}")
        print(f"    B = {b['B_3point_estimate']}")

    if 'pslq_A_3point' in b and b['pslq_A_3point']:
        ps = b['pslq_A_3point']
        if isinstance(ps, dict) and 'reconstruction' in ps:
            print(f"    PSLQ A: {ps['reconstruction']}")
            print(f"    coeffs: {ps.get('coeffs')}")

    if 'pslq_A_best' in b and b['pslq_A_best']:
        ps = b['pslq_A_best']
        if isinstance(ps, dict) and 'reconstruction' in ps:
            print(f"    PSLQ A_best: {ps['reconstruction']}")

    # Save everything
    payload = {
        'description': 'Sprint MR-B: spectral action convergence rate vs M2 ring',
        'mp_dps': mp.dps,
        'part_A': a,
        'part_B': b,
    }
    OUTPATH.write_text(json.dumps(payload, indent=2, default=str))
    print(f"\nSaved: {OUTPATH}")

if __name__ == '__main__':
    main()
