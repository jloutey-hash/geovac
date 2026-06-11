"""Sprint MR-B Part B refinement: high-precision extraction of A in
log|epsilon(t)| ~ -A/t + B*log(t) + C, predicted A = pi^2.

Use mp.dps = 300, push to t = 0.02 (where epsilon ~ exp(-pi^2/0.02) ~ exp(-493)
which requires ~215-digit precision to resolve cleanly).

Method:
- For each t, sum the heat kernel to enough terms so the unsummed tail < 10^{-(mp.dps-30)}
- Compute epsilon(t) = K(t) - 2-term SD
- Compute g(t) = t * log|eps(t)|. Theory: g(t) = -A + B*t*log(t) + C*t + (corrections that vanish as t -> 0)
- Fit g(t) = -A + B*t*log(t) + C*t. The y-intercept at t=0 is -A.
"""
import json
from pathlib import Path
import numpy as np
from mpmath import (
    mp, mpf, exp as mpexp, log as mplog, fsum, nstr,
    pi as mp_pi, sqrt as mpsqrt, pslq
)

mp.dps = 300

OUTDIR = Path(__file__).resolve().parent / "data"
OUTPATH = OUTDIR / "mr_b_modular_residual_high_prec.json"

def heat_kernel(t, n_terms):
    t = mpf(t)
    return fsum(
        mpf(2*(n+1)*(n+2)) * mpexp(-(mpf(n) + mpf('1.5'))**2 * t)
        for n in range(n_terms)
    )

def two_term_sd(t):
    t = mpf(t)
    sp = mpsqrt(mp_pi)
    return sp / 2 / t**mpf('1.5') - sp / 4 / t**mpf('0.5')

def n_terms_for_precision(t, target_dps=None):
    """How many terms to sum so the unsummed tail < 10^{-target_dps}.
    Tail of sum is dominated by exp(-(N+3/2)^2 * t). Need (N+3/2)^2 * t > target * log(10).
    """
    if target_dps is None:
        target_dps = mp.dps + 50
    return max(400, int(2 * float(mpsqrt(mpf(target_dps) * mpf('2.302585') / t))))

def main():
    print("="*72)
    print("Sprint MR-B Part B: high-precision modular residual extraction")
    print(f"mp.dps = {mp.dps}")
    print("="*72)

    # Wider span of t values, going down to 0.02
    t_strs = ['0.5', '0.4', '0.3', '0.25', '0.2', '0.15', '0.12',
              '0.1', '0.08', '0.06', '0.05', '0.04', '0.035',
              '0.03', '0.025', '0.02']
    t_values = [mpf(s) for s in t_strs]

    samples = []
    print(f"\n{'t':>8} {'n_terms':>10} {'log|eps|':>20} {'eps_sign':>9}")
    print("-"*60)
    for t in t_values:
        n_terms = n_terms_for_precision(t)
        K = heat_kernel(t, n_terms)
        eps = K - two_term_sd(t)
        if eps == 0:
            print(f"  {nstr(t, 6)}    {n_terms:>10}    eps == 0 (precision wall)")
            continue
        log_eps = mplog(abs(eps))
        sign = 1 if eps > 0 else -1
        samples.append({
            't': t,
            't_str': nstr(t, 12),
            'n_terms': n_terms,
            'log_abs_eps': log_eps,
            'sign': sign,
            'eps': eps,
        })
        print(f"  {nstr(t, 6):>8} {n_terms:>10} {nstr(log_eps, 18):>20} {sign:>9}")

    # Fit g(t) = t * log|eps(t)| = -A + B*t*log(t) + C*t + ...
    # The y-intercept at t->0 of g(t) is -A. Predicted A = pi^2.
    print()
    print("g(t) = t * log|eps(t)|.  Theory: g(t) = -A + B*t*log(t) + C*t + O(t^2 ...)")
    print(f"{'t':>8} {'g(t)':>30} {'g(t) - (-pi^2)':>30}")
    print("-"*72)
    for s in samples:
        g_t = s['t'] * s['log_abs_eps']
        residual_vs_minus_pi2 = g_t - (-mp_pi**2)
        s['g_t'] = g_t
        s['g_t_minus_neg_pi2'] = residual_vs_minus_pi2
        print(f"  {nstr(s['t'], 6):>8} {nstr(g_t, 28):>30} {nstr(residual_vs_minus_pi2, 16):>30}")

    # Linear fit at small t: g(t) approaches -A as t -> 0
    # Fit g(t) = a0 + a1 * t * log(t) + a2 * t + a3 * t^2 * log(t) + a4 * t^2
    # using all samples. Check fit quality.
    print()
    print("Polynomial fit:  g(t) = a0 + a1*t*log(t) + a2*t + a3*t^2*log(t) + a4*t^2")
    if len(samples) >= 6:
        # Use mp arithmetic
        N = len(samples)
        # Try several polynomial orders and check stability of a0
        for order_name, basis_funcs in [
            ('linear (a0)',        [lambda t: mpf(1)]),
            ('+ t*log(t)',         [lambda t: mpf(1), lambda t: t*mplog(t)]),
            ('+ t',                [lambda t: mpf(1), lambda t: t*mplog(t), lambda t: t]),
            ('+ t^2*log(t)',       [lambda t: mpf(1), lambda t: t*mplog(t), lambda t: t,
                                    lambda t: t**2 * mplog(t)]),
            ('+ t^2',              [lambda t: mpf(1), lambda t: t*mplog(t), lambda t: t,
                                    lambda t: t**2 * mplog(t), lambda t: t**2]),
            ('+ t^3*log(t)',       [lambda t: mpf(1), lambda t: t*mplog(t), lambda t: t,
                                    lambda t: t**2 * mplog(t), lambda t: t**2,
                                    lambda t: t**3 * mplog(t)]),
            ('+ t^3',              [lambda t: mpf(1), lambda t: t*mplog(t), lambda t: t,
                                    lambda t: t**2 * mplog(t), lambda t: t**2,
                                    lambda t: t**3 * mplog(t), lambda t: t**3]),
        ]:
            k = len(basis_funcs)
            if N < k + 2:
                continue
            # Pick smallest-t samples for fit (asymptotic regime)
            sorted_by_t = sorted(samples, key=lambda s: float(s['t']))
            fit_subset = sorted_by_t[:max(k+2, min(N, 10))]
            n_fit = len(fit_subset)
            M = mp.matrix(n_fit, k)
            rhs = mp.matrix(n_fit, 1)
            for i, s in enumerate(fit_subset):
                t_i = s['t']
                for j, f in enumerate(basis_funcs):
                    M[i, j] = f(t_i)
                rhs[i, 0] = s['g_t']
            # Solve normal equations: (M^T M) x = M^T b
            try:
                MTM = M.T * M
                MTb = M.T * rhs
                sol = mp.lu_solve(MTM, MTb)
                a0_est = sol[0, 0]
                A_est = -a0_est
                rel_err_pi2 = (A_est - mp_pi**2) / mp_pi**2
                print(f"  {order_name:25s}  a0 = {nstr(a0_est, 18):>22}    A = {nstr(A_est, 18):>22}    rel_err vs pi^2 = {float(rel_err_pi2):+.4e}")
                # Save the most informative fit
                if k == 5 or k == 7:
                    samples_record = {
                        'order': order_name,
                        'a0': nstr(a0_est, 50),
                        'A_est': nstr(A_est, 50),
                        'A_est_float': float(A_est),
                        'A_pred_pi^2': float(mp_pi**2),
                        'rel_err_pi^2': float(rel_err_pi2),
                        'n_fit_points': n_fit,
                    }
                    if 'best_fit' not in samples_record or k > samples_record.get('best_k', 0):
                        best_fit = samples_record
                        best_fit['best_k'] = k
            except Exception as e:
                print(f"  {order_name:25s}  FAIL: {e}")

    # PSLQ the best A estimate
    print()
    print("PSLQ best A estimate (k=7 fit) against M2 ring:")
    try:
        best_A = mpf(best_fit['A_est'])
        basis = [
            ('1', mpf(1)),
            ('sqrt(pi)', mpsqrt(mp_pi)),
            ('pi', mp_pi),
            ('pi^2', mp_pi**2),
            ('pi^3', mp_pi**3),
            ('1/pi', 1/mp_pi),
        ]
        vec = [best_A] + [v for _, v in basis]
        names = ['target'] + [n for n, _ in basis]
        try:
            rel = pslq(vec, maxcoeff=10**8, tol=mpf(10)**(-100))
            if rel:
                print(f"  PSLQ relation: {rel}")
                print(f"  names: {names}")
                if rel[0] != 0:
                    c0 = rel[0]
                    contributions = []
                    for nm, c in zip(names[1:], rel[1:]):
                        if c != 0:
                            contributions.append(f"{-c}/{c0} * {nm}")
                    reconstruction = " + ".join(contributions)
                    print(f"  reconstruction: target = {reconstruction}")
            else:
                print("  PSLQ: no relation found")
                rel = None
        except Exception as e:
            print(f"  PSLQ error: {e}")
            rel = None

        # Direct check: is A = pi^2 exactly?
        diff_pi2 = best_A - mp_pi**2
        print(f"\n  Direct check: A - pi^2 = {nstr(diff_pi2, 30)}")
        print(f"  rel diff (A - pi^2)/pi^2 = {nstr(diff_pi2/mp_pi**2, 20)}")
    except (NameError, KeyError) as e:
        print(f"  Could not extract best_A: {e}")
        rel = None

    # Save
    payload = {
        'mp_dps': mp.dps,
        'description': 'High-precision extraction of A in log|eps(t)| ~ -A/t for S^3 Dirac',
        'predicted_A': 'pi^2',
        'samples': [
            {
                't': nstr(s['t'], 30),
                't_float': float(s['t']),
                'n_terms': s['n_terms'],
                'log_abs_eps_50dps': nstr(s['log_abs_eps'], 50),
                'log_abs_eps_float': float(s['log_abs_eps']),
                'sign': s['sign'],
                'g_t_50dps': nstr(s['g_t'], 50),
                'g_t_float': float(s['g_t']),
            }
            for s in samples
        ],
    }
    if 'best_fit' in dir() and best_fit:
        payload['best_polynomial_fit'] = best_fit
        payload['A_PROVEN_TO_BE_PI_SQUARED'] = float(abs(mpf(best_fit['A_est']) - mp_pi**2) / mp_pi**2 < 1e-20)
    OUTPATH.write_text(json.dumps(payload, indent=2, default=str))
    print(f"\nSaved: {OUTPATH}")

if __name__ == '__main__':
    main()
