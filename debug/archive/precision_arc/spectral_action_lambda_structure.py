"""Deep structural analysis of Lambda_inf = 3.71025 for K/pi Gaussian match.

Tests whether Lambda_inf has a structural origin in S^3 Dirac spectral data,
or is just the IVT-guaranteed root of a monotone function.
"""
import numpy as np
from scipy.optimize import brentq
from mpmath import mp, mpf, exp as mpexp, fsum, nstr, pi as mp_pi, zeta, sqrt as mpsqrt

mp.dps = 50

K_over_pi_exact = mpf('43.619934066848226')
B = mpf('42')
F = mp_pi**2 / 6
Delta = mpf('1') / 40

def continuum_gaussian_action(L, n_terms=10000):
    L2 = L * L
    terms = []
    for n in range(n_terms):
        g = 2 * (n+1) * (n+2)
        lam_sq = (n + mpf('1.5'))**2
        terms.append(g * mpexp(-lam_sq / L2))
    return fsum(terms)

def continuum_gaussian_np(L, n_terms=10000):
    ns = np.arange(n_terms)
    gs = 2.0 * (ns+1) * (ns+2)
    lam_sqs = (ns + 1.5)**2
    return float(np.sum(gs * np.exp(-lam_sqs / L**2)))

# Precise continuum Lambda for K/pi
f_Kpi = lambda L: continuum_gaussian_np(L) - float(K_over_pi_exact)
L_Kpi = brentq(f_Kpi, 1.0, 50.0, xtol=1e-15, rtol=1e-15)

print(f"Lambda_inf = {L_Kpi:.15f}")
print(f"Lambda_inf^2 = {L_Kpi**2:.15f}")
L2 = L_Kpi**2

# Test 1: Is Lambda related to the SD coefficients?
# On unit S^3, the SD coefficients are:
# a_0 = sqrt(pi)  (Dirac on S^3, Volume/(4*pi)^{3/2} * 4)
# a_1 = sqrt(pi)  (curvature correction)
# a_2 = sqrt(pi)/8
# These are from the continuum heat kernel Tr exp(-tD^2) ~ sum a_k t^{k-3/2}
print("\n--- Test 1: SD coefficient relation ---")
a0 = float(mpsqrt(mp_pi))
a1 = float(mpsqrt(mp_pi))
a2 = float(mpsqrt(mp_pi) / 8)
print(f"  a_0 = sqrt(pi) = {a0:.10f}")
print(f"  a_1 = sqrt(pi) = {a1:.10f}")
print(f"  a_2 = sqrt(pi)/8 = {a2:.10f}")
# CC spectral action: Tr f(D^2/L^2) ~ a_0 L^3 f_3 + a_1 L f_1 + a_2/L f_{-1} + ...
# For Gaussian f(u) = e^{-u}: f_3 = Gamma(2) = 1, f_1 = Gamma(1) = 1
# So leading: a_0 * L^3 + a_1 * L = sqrt(pi)*(L^3 + L)
leading = a0 * L_Kpi**3 + a1 * L_Kpi
print(f"  a_0*L^3 + a_1*L = sqrt(pi)*(L^3 + L) = {leading:.6f}")
print(f"  vs K/pi = {float(K_over_pi_exact):.6f}")
print(f"  ratio = {leading/float(K_over_pi_exact):.6f}")
# Check: does sqrt(pi)*(L^3 + L) = K/pi?
# If so, L^3 + L = K/(pi*sqrt(pi)) = K/pi^(3/2)
target_cubic = float(K_over_pi_exact) / a0
print(f"  L^3 + L target = {target_cubic:.10f}")
actual_cubic = L_Kpi**3 + L_Kpi
print(f"  L^3 + L actual = {actual_cubic:.10f}")
print(f"  ratio = {actual_cubic/target_cubic:.10f}")

# Test 2: Is Lambda the solution of L^3 + L = K/pi^{3/2}?
# L^3 + L = 43.62/sqrt(pi) = 43.62/1.7725 = 24.61
print(f"\n--- Test 2: Cubic equation L^3 + L = K/pi^{3/2} ---")
from scipy.optimize import fsolve
K_pi32 = float(K_over_pi_exact) / a0
cubic_root = fsolve(lambda L: L**3 + L - K_pi32, 3.0)[0]
print(f"  K/pi^(3/2) = {K_pi32:.10f}")
print(f"  Cubic root = {cubic_root:.10f}")
print(f"  Lambda_inf = {L_Kpi:.10f}")
print(f"  Diff = {cubic_root - L_Kpi:.6e}")
print(f"  Rel = {(cubic_root - L_Kpi)/L_Kpi:.6e}")

# Test 3: Higher-order SD corrections
# Actually compute: does the SD expansion at Lambda_inf converge to K/pi?
print(f"\n--- Test 3: SD expansion convergence at Lambda = {L_Kpi:.6f} ---")
# Heat kernel: K(t) = Tr exp(-t D^2)
# For Dirac on unit S^3: K(t) = sum_{n=0}^inf 2(n+1)(n+2) exp(-(n+3/2)^2 t)
# SD expansion: K(t) ~ a_0 t^{-3/2} + a_1 t^{-1/2} + a_2 t^{1/2} + a_3 t^{3/2} + ...
# Spectral action: Tr f(D^2/L^2) = integral_0^inf K(t/L^2) f_tilde(t) dt
# For Gaussian f(u) = e^{-u}: Tr exp(-D^2/L^2) = K(1/L^2)
# = a_0 (1/L^2)^{-3/2} + a_1 (1/L^2)^{-1/2} + a_2 (1/L^2)^{1/2} + ...
# = a_0 L^3 + a_1 L + a_2/L + a_3/L^3 + ...
# Need to know a_3 etc. Fit from exact heat kernel at small t.

# Actually let me just compute the SD terms directly
# K(t) = Tr exp(-t D^2) exact at many t-values, then fit
from numpy.polynomial import polynomial as P

t_vals = np.logspace(-3, -1, 100)
K_vals = np.array([continuum_gaussian_np(1.0/np.sqrt(t)) for t in t_vals])

# K(t) = a_0 t^{-3/2} + a_1 t^{-1/2} + a_2 t^{1/2} + a_3 t^{3/2}
# Multiply by t^{3/2}: t^{3/2} K(t) = a_0 + a_1 t + a_2 t^2 + a_3 t^3
y = K_vals * t_vals**1.5
# Fit polynomial in t
coeffs = np.polyfit(t_vals, y, 6)  # up to t^6 = a_6 t^{7.5-1.5} = a_6 t^6
# Coefficients are in descending order
coeffs_ascending = coeffs[::-1]
print(f"  SD coefficients (from heat kernel fit):")
for k in range(min(7, len(coeffs_ascending))):
    print(f"    a_{k} = {coeffs_ascending[k]:.10f}")

# Now compute SD approximation at Lambda_inf
L = L_Kpi
SD_terms = []
for k in range(7):
    power = 3 - 2*k  # 3, 1, -1, -3, -5, -7, -9
    term = coeffs_ascending[k] * L**power
    SD_terms.append(term)
    cumsum = sum(SD_terms)
    print(f"    SD order {k}: term = {term:+.6e}, cumsum = {cumsum:.6f} "
          f"(vs K/pi = {float(K_over_pi_exact):.6f}, err = {(cumsum-float(K_over_pi_exact))/float(K_over_pi_exact):.2e})")

# Test 4: Is Lambda meaningful as a physical scale on S^3?
print(f"\n--- Test 4: Lambda as physical scale ---")
print(f"  Lowest eigenvalue: |lambda_0| = 3/2 = 1.5")
print(f"  Lambda/|lambda_0| = {L_Kpi/1.5:.6f}")
print(f"  Lambda^2/|lambda_0|^2 = {L_Kpi**2/2.25:.6f}")
print(f"  Number of modes with |lambda| < Lambda:")
n_below = 0
for n in range(1000):
    lam = n + 1.5
    if lam < L_Kpi:
        n_below += 2*(n+1)*(n+2)
print(f"    {n_below} modes (n = 0, 1, 2 contributes 4+12+24 = 40 = Delta^-1)")
print(f"    Note: Lambda = 3.71 cuts between n=2 (|lambda|=3.5) and n=3 (|lambda|=4.5)")
print(f"    This is EXACTLY at n_max = 3 boundary!")
n_included = sum(2*(n+1)*(n+2) for n in range(3))  # n=0,1,2
print(f"    Modes at n <= 2: {n_included} = Delta^-1 = 40")
print(f"    Modes at n <= 3: {n_included + 2*4*5} = {n_included + 40} = 80")

# Test 5: Lambda from implicit equation?
print(f"\n--- Test 5: What if K/pi = Tr exp(-D^2/L^2) defines L? ---")
print(f"  K = pi * Tr exp(-D^2/L^2) at L = {L_Kpi:.10f}")
print(f"  This is NOT a derivation: it's inverting the spectral action")
print(f"  The IVT guarantees a solution for any target in (0, inf)")
print(f"  The question: is L = {L_Kpi:.6f} DERIVABLE from S^3 data alone?")

# Final structural comparison
print(f"\n--- Final verdict ---")
print(f"  Lambda^2 = {L2:.10f}")
# Check: B-2 = 40 = Delta^-1, and Lambda^2 ~ 13.77
print(f"  Lambda^2 / B = {L2/42:.10f}")
print(f"  Lambda^2 / (B-2) = {L2/40:.10f}")
print(f"  Lambda^2 / Delta^-1 = {L2/40:.10f}")
print(f"  B * Delta * Lambda^2 = {42*L2/40:.10f}")
print(f"  Lambda^2 + B = {L2+42:.6f} vs K/pi = {float(K_over_pi_exact):.6f}")
print(f"  Lambda^2 * pi = {L2*np.pi:.6f} vs K = {float(K_over_pi_exact)*np.pi:.6f}")
print(f"  Lambda^2 * e = {L2*np.e:.6f}")
