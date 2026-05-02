"""Verify the SD two-term exactness and compute a_2 directly.

The key observation: at Lambda ~ 3.71, the two-term SD approximation
  K/pi = (sqrt(pi)/2)*L^3 + (-sqrt(pi)/4)*L
gives K/pi to machine precision.

Question: is a_2 genuinely zero, or is it just that a_2/L is negligible?
"""
import numpy as np
from scipy.optimize import brentq
from mpmath import mp, mpf, exp as mpexp, fsum, nstr, pi as mp_pi, sqrt as mpsqrt, erfc

mp.dps = 80

K_over_pi = mpf('43.619934066848226')

def continuum_heat_kernel(t, n_terms=50000):
    """K(t) = Tr exp(-t D^2) = sum 2(n+1)(n+2) exp(-(n+3/2)^2 t)."""
    terms = []
    for n in range(n_terms):
        g = 2*(n+1)*(n+2)
        lam_sq = (n + mpf('3')/2)**2
        terms.append(g * mpexp(-lam_sq * t))
    return fsum(terms)

# Compute SD coefficients by direct extraction
# K(t) = a_0 t^{-3/2} + a_1 t^{-1/2} + a_2 t^{1/2} + a_3 t^{3/2} + ...
# a_0 = lim_{t->0} t^{3/2} K(t)
# a_1 = lim_{t->0} [K(t) - a_0 t^{-3/2}] * t^{1/2}
# a_2 = lim_{t->0} [K(t) - a_0 t^{-3/2} - a_1 t^{-1/2}] * t^{-1/2}

print("SD coefficient extraction from continuum S^3 Dirac heat kernel")
print("="*70)

a0_exact = mpsqrt(mp_pi) / 2
a1_exact = -mpsqrt(mp_pi) / 4
print(f"Theoretical a_0 = sqrt(pi)/2 = {nstr(a0_exact, 20)}")
print(f"Theoretical a_1 = -sqrt(pi)/4 = {nstr(a1_exact, 20)}")

# Extract a_2 at various small t
print(f"\nExtracting a_2 from [K(t) - a_0/t^(3/2) - a_1/t^(1/2)] / t^(1/2):")
for t_val in [mpf('1e-2'), mpf('5e-3'), mpf('1e-3'), mpf('5e-4'), mpf('1e-4')]:
    K = continuum_heat_kernel(t_val, n_terms=50000)
    remainder = K - a0_exact / t_val**mpf('1.5') - a1_exact / t_val**mpf('0.5')
    a2_est = remainder / t_val**mpf('0.5')
    print(f"  t = {nstr(t_val, 4)}: a_2 = {nstr(a2_est, 15)}")

# Use the smallest t for the best estimate
t_small = mpf('1e-4')
K_small = continuum_heat_kernel(t_small, n_terms=50000)
remainder = K_small - a0_exact / t_small**mpf('1.5') - a1_exact / t_small**mpf('0.5')
a2_est = remainder / t_small**mpf('0.5')
print(f"\nBest estimate: a_2 = {nstr(a2_est, 20)}")
print(f"  = {nstr(a2_est/mpsqrt(mp_pi), 20)} * sqrt(pi)")

# Check against known formula: a_2 for D^2 on unit S^3
# Gilkey formula with E = R/4 = 3/2, d=3, rank(V)=2
# a_2_curv = (1/180)(5R^2 - 2|Ric|^2 + 2|Rm|^2) = (180-24+24)/180 = 1
# a_2_E = (1/2)E^2 - (1/6)RE = 9/8 - 3/2 = -3/8
# a_2_Omega: spin connection curvature on S^3
# For constant curvature: Omega_ij = (1/4)R_{abij} gamma^a gamma^b
# |Omega|^2 = (1/12)Omega_ij Omega^ij = ...
# On S^3: R_{abcd} = g_ac g_bd - g_ad g_bc, so
# Omega_ij = (1/4)(g_ai g_bj - g_aj g_bi) gamma^a gamma^b = (1/4)(gamma_i gamma_j - gamma_j gamma_i) = (1/2)sigma_ij
# (1/12)Omega_ij Omega^ij = (1/12)*(1/4)*R_{abij}R^{cdij}*tr(gamma^a gamma^b gamma_c gamma_d)
# This is getting complicated. Let me just use the numerical value.

# Predicted: a_2 = (4pi)^{-3/2} * 2 * integral * 2*pi^2
# where integral = a_2_curv + a_2_E + a_2_Omega
# = 1 + (-3/8) + a_2_Omega

print(f"\n--- SD coefficient table ---")
print(f"a_0 = {nstr(a0_exact, 15)} = sqrt(pi)/2")
print(f"a_1 = {nstr(a1_exact, 15)} = -sqrt(pi)/4")
print(f"a_2 = {nstr(a2_est, 15)}")

# Now compute the spectral action and check two-term vs three-term
print(f"\n--- Spectral action at various Lambda ---")
L_Kpi = mpf('3.710245467906053')

for L_val in [mpf('3.5'), mpf('3.71'), L_Kpi, mpf('4.0'), mpf('5.0')]:
    exact = continuum_heat_kernel(1/L_val**2, n_terms=50000)
    two_term = a0_exact * L_val**3 + a1_exact * L_val
    three_term = two_term + a2_est / L_val

    err_2 = exact - two_term
    err_3 = exact - three_term

    print(f"\n  Lambda = {nstr(L_val, 6)}:")
    print(f"    Exact     = {nstr(exact, 20)}")
    print(f"    2-term SD = {nstr(two_term, 20)}")
    print(f"    3-term SD = {nstr(three_term, 20)}")
    print(f"    |exact - 2-term| = {nstr(abs(err_2), 6)}")
    print(f"    |exact - 3-term| = {nstr(abs(err_3), 6)}")

# Key structural question: the cubic
print(f"\n--- Cubic equation ---")
print(f"K/pi = a_0 L^3 + a_1 L + a_2/L + ...")
print(f"     = (sqrt(pi)/4)(2L^3 - L) + a_2/L + ...")
print(f"")
print(f"If a_2/L is negligible:")
print(f"  K/pi = (sqrt(pi)/4)(2L^3 - L)")
print(f"  4*K/(pi*sqrt(pi)) = 2L^3 - L")
target = 4 * K_over_pi / mpsqrt(mp_pi)
print(f"  4*K/pi^(3/2) = {nstr(target, 15)}")
L_cubic = L_Kpi
cubic_val = 2*L_cubic**3 - L_cubic
print(f"  2L^3 - L at L=3.710 = {nstr(cubic_val, 15)}")
print(f"  Ratio = {nstr(cubic_val/target, 15)}")

# Is a_2 / L really negligible?
a2_correction = a2_est / L_Kpi
print(f"\n  a_2/L at L=3.710 = {nstr(a2_correction, 10)}")
print(f"  vs K/pi = {nstr(K_over_pi, 10)}")
print(f"  fractional = {nstr(a2_correction/K_over_pi, 6)}")
