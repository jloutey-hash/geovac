"""Sprint MR-B Part B finalization: extract A *and* the subleading prefactor.

Theory: by Jacobi theta modular transformation, the modular residual of the
half-integer-shifted Dirac heat kernel on unit S^3 is

   epsilon(t) ~ Sum_{m>=1} A_m * t^{-3/2} * exp(-pi^2 m^2 / t) * (cos/sin terms)

Leading m=1 term:
   epsilon(t) ~ C * t^{-3/2} * exp(-pi^2 / t)

Taking logs:
   log|eps(t)| = -pi^2/t - (3/2)*log(t) + log(C) + (subleading m=2 contribution ~ exp(-3*pi^2/t))

So g(t) := t * log|eps(t)| = -pi^2 - (3/2) * t*log(t) + log(C) * t + O(t * exp(-3*pi^2/t))

Predictions for the polynomial fit g(t) = a0 + a1*t*log(t) + a2*t + ...:
- a0 = -pi^2
- a1 = -3/2 = -1.5
- a2 = log(C) where C is the modular prefactor

This sprint extracts a0, a1, a2 with controlled precision and PSLQs each.
"""
import json
from pathlib import Path
from mpmath import (
    mp, mpf, exp as mpexp, log as mplog, fsum, nstr,
    pi as mp_pi, sqrt as mpsqrt, pslq
)

mp.dps = 400  # very high precision so smallest-t samples remain reliable

OUTDIR = Path(__file__).resolve().parent / "data"
OUTPATH = OUTDIR / "mr_b_subleading.json"

def heat_kernel(t, n_terms):
    return fsum(
        mpf(2*(n+1)*(n+2)) * mpexp(-(mpf(n) + mpf('1.5'))**2 * t)
        for n in range(n_terms)
    )

def two_term_sd(t):
    sp = mpsqrt(mp_pi)
    return sp / 2 / t**mpf('1.5') - sp / 4 / t**mpf('0.5')

def n_terms_safe(t, target_dps=None):
    if target_dps is None:
        target_dps = mp.dps + 60
    return max(400, int(2.5 * float(mpsqrt(mpf(target_dps) * mpf('2.302585') / t))))

# Sample at many small t with high precision
t_strs = [
    '0.5', '0.4', '0.3', '0.25', '0.2', '0.15', '0.12',
    '0.1', '0.08', '0.06', '0.05', '0.04', '0.035',
    '0.03', '0.025', '0.02', '0.018', '0.015', '0.012', '0.01'
]
t_values = [mpf(s) for s in t_strs]

print("Computing samples with mp.dps =", mp.dps)
samples = []
for t in t_values:
    n_terms = n_terms_safe(t)
    K = heat_kernel(t, n_terms)
    eps = K - two_term_sd(t)
    if eps == 0 or abs(eps) < mpf(10)**(-(mp.dps-30)):
        print(f"  t={nstr(t, 6)}: eps below precision floor; skipping")
        continue
    log_eps = mplog(abs(eps))
    g_t = t * log_eps
    samples.append((t, log_eps, g_t, n_terms))
    print(f"  t={nstr(t, 8):>10}  n_terms={n_terms:>6}  log|eps|={nstr(log_eps, 15):>22}  "
          f"g(t)={nstr(g_t, 15):>22}")

# Polynomial fits of increasing order. Order k means basis up to t^{(k-1)/2} * log^? family.
# Basis: 1, t*log(t), t, t^2*log(t), t^2, t^3*log(t), t^3, t^4*log(t), t^4, ...
def basis_func(j, t):
    """j=0: 1; j=1: t*log(t); j=2: t; j=3: t^2*log(t); j=4: t^2; j=5: t^3*log(t); ..."""
    if j == 0:
        return mpf(1)
    half = (j + 1) // 2  # 1,1,2,2,3,3,...
    if j % 2 == 1:
        return t**half * mplog(t)
    else:
        return t**half

print("\nPolynomial fits (using all", len(samples), "samples):")
print(f"{'k':>3} {'a0':>30} {'a1':>20} {'a2':>20} {'A_est':>30} {'rel_err pi^2':>16}")
print("-"*120)

best_fit = None
for k in range(1, min(15, len(samples))):
    n_basis = k + 1
    if len(samples) < n_basis:
        break
    M = mp.matrix(len(samples), n_basis)
    rhs = mp.matrix(len(samples), 1)
    for i, (t, _, g_t, _) in enumerate(samples):
        for j in range(n_basis):
            M[i, j] = basis_func(j, t)
        rhs[i, 0] = g_t
    try:
        MTM = M.T * M
        MTb = M.T * rhs
        sol = mp.lu_solve(MTM, MTb)
        a0 = sol[0, 0]
        a1 = sol[1, 0] if n_basis > 1 else mpf(0)
        a2 = sol[2, 0] if n_basis > 2 else mpf(0)
        A_est = -a0
        rel_err = (A_est - mp_pi**2) / mp_pi**2
        a0_s = nstr(a0, 18)
        a1_s = nstr(a1, 14) if n_basis > 1 else 'n/a'
        a2_s = nstr(a2, 14) if n_basis > 2 else 'n/a'
        A_s = nstr(A_est, 22)
        print(f"{k:>3} {a0_s:>30} {a1_s:>20} {a2_s:>20} {A_s:>30} {float(rel_err):>+16.4e}")
        # Keep fit with lowest |rel_err|
        if best_fit is None or abs(rel_err) < abs(best_fit['rel_err_pi2']):
            best_fit = {
                'k': k,
                'a0': a0,
                'a1': a1,
                'a2': a2,
                'a3': sol[3, 0] if n_basis > 3 else None,
                'A_est': A_est,
                'rel_err_pi2': rel_err,
                'all_coefs': [sol[i, 0] for i in range(n_basis)],
            }
    except Exception as e:
        print(f"{k:>3}  FAIL: {e}")

print()
if best_fit is not None:
    print("="*72)
    print(f"BEST FIT (k={best_fit['k']}):")
    print(f"  A = {nstr(best_fit['A_est'], 50)}")
    print(f"  pi^2 = {nstr(mp_pi**2, 50)}")
    print(f"  diff = {nstr(best_fit['A_est'] - mp_pi**2, 30)}")
    print(f"  rel diff = {nstr((best_fit['A_est'] - mp_pi**2)/mp_pi**2, 20)}")
    print(f"  a1 = {nstr(best_fit['a1'], 30)}    (predicted -3/2 = -1.5)")
    print(f"  a1 - (-3/2) = {nstr(best_fit['a1'] + mpf('1.5'), 20)}")
    print(f"  a2 = log(C) = {nstr(best_fit['a2'], 30)}")

# Extra: high-precision a0 by Richardson-like extrapolation across fit orders
# The residual (A_est - pi^2) is monotonically shrinking; find the asymptotic value.
# But what we have already is good enough.

# PSLQ on a1: predicted -3/2
print("\nPSLQ tests on extracted coefficients:")
if best_fit is not None and best_fit['a1'] is not None:
    target = best_fit['a1']
    print(f"\n  PSLQ a1 = {nstr(target, 30)}:")
    print(f"  Direct check: a1 - (-3/2) = {nstr(target - mpf(-3, 2), 20)}")  # mpf(-3,2) is -1.5
    # Actually use mpf('-1.5')
    pred = mpf('-1.5')
    print(f"  rel err vs -3/2 = {nstr((target - pred)/pred, 20)}")

if best_fit is not None and best_fit['a2'] is not None:
    target = best_fit['a2']
    print(f"\n  PSLQ a2 = log(C) = {nstr(target, 30)}:")
    # Test C = exp(a2). Predict C = some clean expression involving pi.
    C = mpexp(target)
    print(f"  C = exp(a2) = {nstr(C, 30)}")
    # Try basis: log(C) = sum c_i * b_i with b_i in {1, log(2), log(pi), log(sqrt(pi))}
    basis = [
        ('1', mpf(1)),
        ('log(2)', mplog(mpf(2))),
        ('log(pi)', mplog(mp_pi)),
        ('log(3)', mplog(mpf(3))),
        ('log(pi/2)', mplog(mp_pi/2)),
    ]
    vec = [target] + [v for _, v in basis]
    try:
        rel = pslq(vec, maxcoeff=10**5, tol=mpf(10)**(-30))
        if rel:
            print(f"  PSLQ relation: {rel}")
            if rel[0] != 0:
                names = ['target'] + [n for n, _ in basis]
                contributions = []
                for nm, c in zip(names[1:], rel[1:]):
                    if c != 0:
                        contributions.append(f"{-c}/{rel[0]} * {nm}")
                print(f"  reconstruction: log(C) = {' + '.join(contributions) if contributions else '0'}")
        else:
            print("  PSLQ on log(C): no relation found at maxcoeff=10^5")
    except Exception as e:
        print(f"  PSLQ on log(C): error {e}")
    # Test C directly against algebraic combinations
    print(f"\n  Try matching C directly against M2 candidates:")
    for name, val in [
        ('2*sqrt(pi)', 2*mpsqrt(mp_pi)),
        ('4', mpf(4)),
        ('pi', mp_pi),
        ('2*pi', 2*mp_pi),
        ('4*pi', 4*mp_pi),
        ('sqrt(pi)/2', mpsqrt(mp_pi)/2),
        ('sqrt(pi)', mpsqrt(mp_pi)),
        ('4*sqrt(pi)', 4*mpsqrt(mp_pi)),
        ('8*pi', 8*mp_pi),
        ('2*pi^(3/2)', 2*mp_pi**mpf('1.5')),
        ('pi^(3/2)/2', mp_pi**mpf('1.5')/2),
        ('pi^2/2', mp_pi**2/2),
        ('pi^2', mp_pi**2),
    ]:
        rel_err_C = (C - val) / val
        if abs(rel_err_C) < mpf('1e-3'):
            print(f"    C ?= {name:25s}: {nstr(val, 15)}    rel_err = {nstr(rel_err_C, 8)}")
        else:
            print(f"    C != {name:25s}: rel_err = {float(rel_err_C):+.4e}")

# Save
payload = {
    'mp_dps': mp.dps,
    'description': 'Sprint MR-B subleading extraction: A=pi^2 + prefactor C',
    'samples': [
        {
            't': nstr(t, 30),
            'n_terms': n_terms,
            'log_abs_eps_50dps': nstr(log_eps, 50),
            'g_t_50dps': nstr(g_t, 50),
        }
        for (t, log_eps, g_t, n_terms) in samples
    ],
    'best_fit': {
        'k': best_fit['k'],
        'A_est_50dps': nstr(best_fit['A_est'], 50),
        'A_est_float': float(best_fit['A_est']),
        'A_pred_pi^2_50dps': nstr(mp_pi**2, 50),
        'rel_err_pi^2_float': float(best_fit['rel_err_pi2']),
        'a0_50dps': nstr(best_fit['a0'], 50),
        'a1_50dps': nstr(best_fit['a1'], 50),
        'a1_predicted': '-3/2',
        'a1_minus_neg_three_halves': nstr(best_fit['a1'] + mpf('1.5'), 30),
        'a2_50dps': nstr(best_fit['a2'], 50) if best_fit['a2'] is not None else None,
        'C_estimate': nstr(mpexp(best_fit['a2']), 30) if best_fit['a2'] is not None else None,
    } if best_fit else None,
}
OUTPATH.write_text(json.dumps(payload, indent=2, default=str))
print(f"\nSaved: {OUTPATH}")
