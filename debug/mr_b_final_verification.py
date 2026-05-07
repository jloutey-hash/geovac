"""Sprint MR-B final verification.

Verified by analytical Jacobi modular transformation:
  K(t) = -d/dt[theta_2(t)] - (1/4) * theta_2(t)
  theta_2(t) = sum_{m in Z+1/2} exp(-m^2 t) = sqrt(pi/t) * [1 - 2 exp(-pi^2/t) + 2 exp(-4 pi^2/t) - ...]

Differentiating and combining:
  K(t) = (sqrt(pi)/2) t^{-3/2} - (sqrt(pi)/4) t^{-1/2}              [2-term SD]
         + 2 sqrt(pi) pi^2 t^{-5/2} exp(-pi^2/t)
         - sqrt(pi) t^{-3/2} exp(-pi^2/t)
         + (1/2) sqrt(pi) t^{-1/2} exp(-pi^2/t)
         + O(exp(-4 pi^2/t))

So:
  epsilon(t) = K(t) - 2-term SD
             = exp(-pi^2/t) * sqrt(pi) * [2 pi^2 t^{-5/2} - t^{-3/2} + (1/2) t^{-1/2}] + O(exp(-4 pi^2/t))

Leading at small t: epsilon(t) ~ 2 pi^{5/2} * t^{-5/2} * exp(-pi^2/t)

This script verifies all three predictions:
  (i)   exponent = pi^2                           -- M2 signature in exponent
  (ii)  power of t = -5/2                         -- rational, structural
  (iii) prefactor C = 2 pi^{5/2}                  -- M2 signature in prefactor
"""
import json
from pathlib import Path
from mpmath import (
    mp, mpf, exp as mpexp, log as mplog, fsum, nstr,
    pi as mp_pi, sqrt as mpsqrt
)

mp.dps = 400

OUTDIR = Path(__file__).resolve().parent / "data"
OUTPATH = OUTDIR / "mr_b_final_verification.json"

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

# Predicted closed form for the modular residual:
def epsilon_predicted_leading(t):
    """Leading m=1 contribution: exp(-pi^2/t) * sqrt(pi) * [2 pi^2 t^{-5/2} - t^{-3/2} + (1/2) t^{-1/2}]"""
    return mpexp(-mp_pi**2 / t) * mpsqrt(mp_pi) * (
        2 * mp_pi**2 * t**mpf('-2.5')
        - t**mpf('-1.5')
        + mpf('0.5') * t**mpf('-0.5')
    )

def epsilon_predicted_two_terms(t):
    """Two modular terms: m=1 and m=2 (m=2 gives exp(-4 pi^2 / t))."""
    sp = mpsqrt(mp_pi)
    leading = mpexp(-mp_pi**2 / t) * sp * (
        2 * mp_pi**2 * t**mpf('-2.5')
        - t**mpf('-1.5')
        + mpf('0.5') * t**mpf('-0.5')
    )
    # m=2 term: theta_2 has +2 exp(-4 pi^2 / t), so similar derivation gives
    # 2 * (factor of 2 from -2 in bracket * -1 from sign at m=2: -2 * (-1)^2 = +2)
    # K_modular_m=2 contribution: -(d/dt)[sqrt(pi/t) * 2 (-1)^2 exp(-4pi^2/t)] - (1/4)[sqrt(pi/t)*2*exp(-4pi^2/t)]
    # but sign is positive at m=2 (theta_4 has alternating; at m=2 sign is +)
    # Actually theta_4(s) = 1 - 2*e^-s + 2*e^-4s - ... so theta_2(t) = sqrt(pi/t)[1 - 2e^-pi^2/t + 2e^-4pi^2/t - ...]
    # m=2 contributes: theta_2_m2(t) = sqrt(pi/t) * 2 * exp(-4 pi^2 / t)
    # K_modular_m2 = -d/dt[theta_2_m2] - (1/4) theta_2_m2
    #              = -d/dt[sqrt(pi/t) * 2 exp(-4 pi^2 / t)] - (1/4) sqrt(pi/t) * 2 exp(-4 pi^2 / t)
    # d/dt[sqrt(pi/t) * 2 exp(-4 pi^2 / t)]
    # = 2[(-1/2)sqrt(pi)*t^{-3/2} * exp(-4pi^2/t) + sqrt(pi/t) * (4 pi^2 / t^2) * exp(-4 pi^2 / t)]
    # = exp(-4 pi^2 / t) * [- sqrt(pi) t^{-3/2} + 8 pi^{5/2} t^{-5/2}]
    # K_modular_m2 = -[-sqrt(pi) t^{-3/2} + 8 pi^{5/2} t^{-5/2}] e^{-4pi^2/t} - (1/2) sqrt(pi) t^{-1/2} e^{-4pi^2/t}
    #              = e^{-4 pi^2 / t} [sqrt(pi) t^{-3/2} - 8 pi^{5/2} t^{-5/2} - (1/2) sqrt(pi) t^{-1/2}]
    # Wait, but theta_4 has alternating signs. Let me reread.
    # theta_4(s) = sum_{n=-inf}^{inf} (-1)^n e^{-n^2 s} = 1 - 2e^{-s} + 2e^{-4s} - 2e^{-9s} + ...
    # So m=1 contribution = -2 e^{-s} (in theta_4); m=2 contribution = +2 e^{-4s}.
    # In theta_2(t) = sqrt(pi/t) theta_4(pi^2/t):
    #   m=0: sqrt(pi/t) * 1
    #   m=1: sqrt(pi/t) * (-2) exp(-pi^2 / t)         <- this is what we used above
    #   m=2: sqrt(pi/t) * (+2) exp(-4 pi^2 / t)
    # So m=2 gives +2 sqrt(pi/t) exp(-4 pi^2 / t).
    # The sign is FLIPPED relative to m=1 (which had -2). But the structural form is the same with pi^2 -> 4*pi^2.
    # So m=2 contribution to epsilon(t) is:
    # eps_m2(t) = -2 (sign factor: m=1 had -2, m=2 has +2, ratio=-1) * (replace pi^2 with 4 pi^2 in exponent and prefactor)
    # i.e. eps_m2 = -1 * exp(-4 pi^2 / t) * sqrt(pi) * [2 (4 pi^2) t^{-5/2} - t^{-3/2} + (1/2) t^{-1/2}]
    #             = -exp(-4 pi^2 / t) * sqrt(pi) * [8 pi^2 t^{-5/2} - t^{-3/2} + (1/2) t^{-1/2}]
    # Hmm let me re-derive carefully without the shortcut.
    # epsilon = K - 2-term SD; with theta_2(t) = sqrt(pi/t)[1 - 2e^-x + 2e^-4x - ...] where x=pi^2/t
    # The "1" gives the SD; the rest gives epsilon.
    # epsilon = -d/dt[sqrt(pi/t)*(-2 e^-x + 2 e^-4x - ...)] - (1/4)*sqrt(pi/t)*(-2e^-x + 2e^-4x - ...)
    # Define f_k(t) = sqrt(pi/t) * sign_k * 2 * exp(-k^2 pi^2 / t), where sign_k = (-1)^k for k>=1. So sign_1=-1, sign_2=+1, sign_3=-1.
    # eps = sum_{k>=1} [-df_k/dt - (1/4) f_k]
    # df_k/dt = sign_k * 2 * [(-1/2) sqrt(pi) t^{-3/2} exp(-k^2 pi^2/t) + sqrt(pi/t) * (k^2 pi^2 / t^2) * exp(-k^2 pi^2/t)]
    #         = sign_k * 2 * exp(-k^2 pi^2 / t) * sqrt(pi) * [-(1/2) t^{-3/2} + k^2 pi^2 t^{-5/2}]
    # -df_k/dt = sign_k * 2 * exp(-k^2 pi^2 / t) * sqrt(pi) * [(1/2) t^{-3/2} - k^2 pi^2 t^{-5/2}]
    #          = sign_k * exp(-k^2 pi^2 / t) * sqrt(pi) * [t^{-3/2} - 2 k^2 pi^2 t^{-5/2}]
    # -(1/4) f_k = -(1/4) sign_k * 2 * sqrt(pi/t) * exp(-k^2 pi^2 / t) = -(1/2) sign_k * sqrt(pi) t^{-1/2} exp(-k^2 pi^2 / t)
    # Sum:
    # contribution_k = sign_k * sqrt(pi) * exp(-k^2 pi^2 / t) * [t^{-3/2} - 2 k^2 pi^2 t^{-5/2} - (1/2) t^{-1/2}]
    # Wait, let me re-examine my m=1 derivation. I had:
    #   eps_m1 = exp(-pi^2/t) * sqrt(pi) * [2 pi^2 t^{-5/2} - t^{-3/2} + (1/2) t^{-1/2}]
    # That has signs: +2 pi^2 t^{-5/2}, -t^{-3/2}, +(1/2) t^{-1/2}.
    # sign_1 = -1 in my formula above. Then sign_k * [t^{-3/2} - 2 k^2 pi^2 t^{-5/2} - (1/2) t^{-1/2}]
    #        = -1 * [t^{-3/2} - 2 pi^2 t^{-5/2} - (1/2) t^{-1/2}]
    #        = -t^{-3/2} + 2 pi^2 t^{-5/2} + (1/2) t^{-1/2}
    # That matches! So my formula is consistent.
    # m=2: sign_2 = +1, factor [t^{-3/2} - 8 pi^2 t^{-5/2} - (1/2) t^{-1/2}]
    #     = t^{-3/2} - 8 pi^2 t^{-5/2} - (1/2) t^{-1/2}
    # multiplied by sqrt(pi) * exp(-4 pi^2 / t)
    m2 = mpexp(-4 * mp_pi**2 / t) * sp * (
        t**mpf('-1.5')
        - 8 * mp_pi**2 * t**mpf('-2.5')
        - mpf('0.5') * t**mpf('-0.5')
    )
    return leading + m2

# ======== Verification ========

print("="*72)
print(f"Sprint MR-B Final Verification (mp.dps = {mp.dps})")
print("="*72)
print()
print("Predicted closed form:")
print("  epsilon(t) = exp(-pi^2/t) * sqrt(pi) * [2 pi^2 t^{-5/2} - t^{-3/2} + (1/2) t^{-1/2}]")
print("             + (-1) * exp(-4 pi^2/t) * sqrt(pi) * [8 pi^2 t^{-5/2} - t^{-3/2} + (1/2) t^{-1/2}]")
print("             + O(exp(-9 pi^2/t))")
print()

# Tabulate measured epsilon(t) vs predicted at multiple t values
t_values = [mpf(s) for s in
            ['0.5', '0.3', '0.2', '0.15', '0.1', '0.08', '0.06', '0.04',
             '0.03', '0.025', '0.02', '0.018', '0.015', '0.012']]
results = []
print(f"{'t':>10} {'measured eps':>30} {'predicted eps (m=1+m=2)':>30} {'ratio':>20} {'rel err':>14}")
print("-"*110)
for t in t_values:
    n_terms = n_terms_safe(t)
    K = heat_kernel(t, n_terms)
    eps_measured = K - two_term_sd(t)
    eps_predicted_full = epsilon_predicted_two_terms(t)
    eps_predicted_leading = epsilon_predicted_leading(t)
    if abs(eps_predicted_full) > 0:
        ratio = eps_measured / eps_predicted_full
        rel_err = (eps_measured - eps_predicted_full) / eps_predicted_full
    else:
        ratio = mpf('nan')
        rel_err = mpf('nan')
    results.append({
        't': nstr(t, 12),
        't_float': float(t),
        'eps_measured_50dps': nstr(eps_measured, 50),
        'eps_predicted_m1_plus_m2_50dps': nstr(eps_predicted_full, 50),
        'eps_predicted_m1_only_50dps': nstr(eps_predicted_leading, 50),
        'ratio_measured_over_predicted': nstr(ratio, 30),
        'ratio_minus_one': nstr(ratio - 1, 20),
        'rel_err_m1_only': nstr((eps_measured - eps_predicted_leading) / eps_predicted_leading, 20),
        'rel_err_m1_plus_m2': nstr(rel_err, 20),
    })
    print(f"  {nstr(t, 6):>8}  {nstr(eps_measured, 12):>28}  "
          f"{nstr(eps_predicted_full, 12):>28}  {nstr(ratio, 16):>20}  {float(rel_err):>+14.4e}")

# The quality of the m=1+m=2 prediction tells us what the next correction is.
# If the m=1+m=2 prediction matches to ~10^{-16} at t=0.5 but to ~10^{-???} at t=0.012,
# the error scales as exp(-9 pi^2 / t) (m=3 contribution).

# Pin down the leading prefactor C = 2 pi^{5/2} via direct extraction:
# At small t, epsilon(t) / [t^{-5/2} * exp(-pi^2/t)] -> 2 pi^{5/2}
print()
print("Direct extraction of leading prefactor:")
print(f"   C := epsilon(t) / [t^{{-5/2}} * exp(-pi^2/t)]   should converge to  2 * pi^{{5/2}}")
print(f"{'t':>10} {'C(t)':>30} {'C - 2pi^(5/2)':>30}")
print("-"*80)
target_C = 2 * mp_pi**mpf('2.5')
print(f"{'predicted':>10} {nstr(target_C, 28):>30} {'0':>30}")
for t in t_values:
    n_terms = n_terms_safe(t)
    K = heat_kernel(t, n_terms)
    eps = K - two_term_sd(t)
    # extract C: eps = C * t^{-5/2} * exp(-pi^2/t) * (1 + corrections)
    C_t = eps / (t**mpf('-2.5') * mpexp(-mp_pi**2 / t))
    diff = C_t - target_C
    print(f"  {nstr(t, 6):>8} {nstr(C_t, 28):>30} {nstr(diff, 16):>30}")

# The C - 2pi^(5/2) values should equal -1*pi^(1/2)*t * pi^... = -sqrt(pi) * t + (1/2) sqrt(pi) * t^2
# from the subleading -t^{-3/2} and +(1/2) t^{-1/2} terms (relative to t^{-5/2}).
# That is, eps/[t^{-5/2} exp] = sqrt(pi) [2 pi^2 - t + (1/2) t^2] + [m=2 correction]
# So C(t) - 2 pi^{5/2} = sqrt(pi) * (-t + (1/2) t^2) + ...
# At t=0.5: prediction = sqrt(pi)*(-0.5 + 0.125) = -0.375*sqrt(pi) = -0.6647...

print()
print("Check the subleading prediction:  C(t) - 2 pi^{5/2} approx sqrt(pi) * (-t + (1/2) t^2):")
print(f"{'t':>10} {'measured C - target':>30} {'sqrt(pi)*(-t + t^2/2)':>30} {'diff':>30}")
print("-"*100)
for t in t_values:
    n_terms = n_terms_safe(t)
    K = heat_kernel(t, n_terms)
    eps = K - two_term_sd(t)
    C_t = eps / (t**mpf('-2.5') * mpexp(-mp_pi**2 / t))
    measured = C_t - target_C
    predicted_sub = mpsqrt(mp_pi) * (-t + (t**2) / 2)
    diff = measured - predicted_sub
    print(f"  {nstr(t, 6):>8} {nstr(measured, 16):>30} {nstr(predicted_sub, 16):>30} {nstr(diff, 16):>30}")

# Save
payload = {
    'mp_dps': mp.dps,
    'description': 'Sprint MR-B final verification: epsilon(t) closed form vs measured',
    'closed_form': {
        'leading_m1': 'exp(-pi^2/t) * sqrt(pi) * [2 pi^2 t^{-5/2} - t^{-3/2} + (1/2) t^{-1/2}]',
        'subleading_m2': '-exp(-4 pi^2/t) * sqrt(pi) * [8 pi^2 t^{-5/2} - t^{-3/2} + (1/2) t^{-1/2}]',
        'prefactor_C_leading': '2 pi^{5/2}',
        'prefactor_C_value_50dps': nstr(target_C, 50),
        'prefactor_C_value_float': float(target_C),
    },
    'comparison_table': results,
}
OUTPATH.write_text(json.dumps(payload, indent=2, default=str))
print(f"\nSaved: {OUTPATH}")
