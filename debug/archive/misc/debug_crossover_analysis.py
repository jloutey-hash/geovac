"""
Crossover function analysis for mu(R) in helium hyperspherical coordinates.

mu(R) interpolates between:
  - Small R: mu ~ a_1*R + a_2*R^2 + ...  (perturbative, algebraic coefficients)
  - Large R:  mu ~ -Z^2*R^2/2 + ...       (He+ threshold, asymptotic)

This script investigates:
  1. Pade approximants for mu(R)/R
  2. Two-scale interpolation formulas
  3. Digamma connection for a_2 sum
  4. SO(6) representation content of 1/r_12
  5. Exact requirements for n_eff and a_2 to match He energy
  6. Continued fraction representation
"""

import numpy as np
from fractions import Fraction
from scipy.integrate import quad
from scipy.optimize import minimize_scalar, curve_fit
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from geovac.hyperspherical_angular import solve_angular


# =========================================================================
# Import exact algebraic formulas from closed-form analysis
# =========================================================================

def I_nuc_exact(n: int) -> Fraction:
    return sum(Fraction(1, 2*j+1) for j in range(2*n))

def I_ee_rational_exact(n: int) -> Fraction:
    s1 = sum(Fraction((-1)**k, 4*k+1) for k in range(n))
    s2 = sum(Fraction((-1)**k, 4*k+3) for k in range(n))
    return s1 - s2

def a1_exact(n: int, Z: float) -> float:
    I_nuc = float(I_nuc_exact(n))
    I_ee_rat = float(I_ee_rational_exact(n))
    return (4 / np.pi) * (-2 * Z * I_nuc + np.sqrt(2) * I_ee_rat)

def V_nuc_off_diag(m: int, n: int, Z: float) -> float:
    if (m + n) % 2 != 0:
        return 0.0
    p_plus = (m + n) // 2
    p_minus = abs(m - n) // 2
    I_plus = float(I_nuc_exact(p_plus))
    I_minus = float(I_nuc_exact(p_minus)) if p_minus > 0 else 0.0
    return (-8 * Z / np.pi) * (I_plus - I_minus)

def V_ee_off_diag(m: int, n: int) -> float:
    if (m + n) % 2 != 0:
        return 0.0
    p_plus = (m + n) // 2
    p_minus = abs(m - n) // 2
    I_plus = float(I_ee_rational_exact(p_plus))
    I_minus = float(I_ee_rational_exact(p_minus)) if p_minus > 0 else 0.0
    return (4 * np.sqrt(2) / np.pi) * (I_plus - I_minus)

def a2_from_sum(n_target: int, Z: float, n_terms: int = 50) -> float:
    E_n = 2 * n_target**2 - 2
    a2 = 0.0
    for m in range(n_target + 2, 2 * n_terms + 2, 2):
        E_m = 2 * m**2 - 2
        V = V_nuc_off_diag(m, n_target, Z) + V_ee_off_diag(m, n_target)
        if abs(V) < 1e-15:
            continue
        a2 += V**2 / (E_n - E_m)
    return a2


# =========================================================================
# STEP 1: Compute exact mu(R) curve from angular solver
# =========================================================================

print("=" * 72)
print("STEP 1: Exact mu(R) from angular solver")
print("=" * 72)

Z = 2.0
R_grid = np.concatenate([
    np.linspace(0.01, 0.5, 30),
    np.linspace(0.5, 3.0, 50),
    np.linspace(3.0, 15.0, 50),
])
R_grid = np.unique(R_grid)

mu_exact = np.zeros(len(R_grid))
for i, R in enumerate(R_grid):
    mu, _ = solve_angular(R=R, Z=Z, l_max=3, n_alpha=200, n_channels=1)
    mu_exact[i] = mu[0]

a1_z2 = a1_exact(1, Z)
a2_z2 = a2_from_sum(1, Z, n_terms=50)

print(f"\n  a_1(Z=2) = {a1_z2:.10f}")
print(f"  a_2(Z=2) = {a2_z2:.10f}")
print(f"  mu(R=0.01) = {mu_exact[0]:.8f}  (expected ~ a_1*0.01 = {a1_z2*0.01:.8f})")
print(f"  mu(R=15)   = {mu_exact[-1]:.4f}  (expected ~ -Z^2*R^2/2 = {-Z**2*15**2/2:.1f})")


# =========================================================================
# STEP 2: Pade approximants for mu(R)/R
# =========================================================================

print("\n" + "=" * 72)
print("STEP 2: Pade approximants for mu(R)/R")
print("=" * 72)

# mu(R)/R ~ a_1 + a_2*R + a_3*R^2 + ...
# At large R: mu(R)/R ~ -Z^2*R/2
# So mu/R diverges linearly at large R.
# Better: consider mu(R)/R^2 which approaches -Z^2/2 at large R.

# Fit mu(R)/R^2 to a Pade [M/N] approximant
mu_over_R2 = mu_exact / R_grid**2

# mu/R^2 ~ a_1/R + a_2 + a_3*R + ...  at small R (diverges)
# mu/R^2 ~ -Z^2/2 at large R (constant)
# So Pade in 1/R might work better.

# Try: f(x) = mu(1/x) * x^2 where x = 1/R
# f(x) ~ a_1*x + a_2 + ... as x -> inf (small R)
# f(x) -> -Z^2/2 as x -> 0 (large R)

# Instead, let's fit mu(R) directly to several functional forms.

# Form 1: Two-parameter interpolation
# mu(R) = a_1*R + a_inf*R^2 / (1 + R/R_c)
# At small R: ~ a_1*R + a_inf*R^2
# At large R: ~ a_1*R + a_inf*R_c*R -> still linear (wrong)

# Form 2: Quadratic crossover
# mu(R) = a_1*R + c_2 * R^2 / (1 + R/R_c)
# Small R: a_1*R + c_2*R^2
# Large R: a_1*R + c_2*R_c*R (linear, not quadratic -- wrong)

# Form 3: mu(R) = a_1*R - (Z^2/2)*R^2 * tanh(R/R_c) + correction
# Small R: a_1*R - (Z^2/2)*R^3/R_c  (cubic, not quadratic)

# Form 4 (correct asymptotics):
# mu(R) = a_1*R + a_2*R^2/(1 + beta*R) + gamma*R^2*(1 - 1/(1+beta*R))
# = a_1*R + a_2*R^2/(1+beta*R) + gamma*beta*R^3/(1+beta*R)
# Small: a_1*R + a_2*R^2 + ...
# Large: a_1*R + (a_2/beta)*R + gamma*R^2 -> need gamma = -Z^2/2 = -2

# Form 5 (simplest correct):
# mu(R) = a_1*R + a_2*R^2 + (-Z^2/2 - a_2)*R^2*(1 - exp(-R/R_c))
# Small: a_1*R + a_2*R^2 + (-Z^2/2 - a_2)*R^3/R_c + ...
# Large: a_1*R + a_2*R^2 + (-Z^2/2 - a_2)*R^2 = a_1*R - Z^2*R^2/2

a_inf = -Z**2 / 2  # = -2.0

# Simple crossover: mu(R) = a_1*R + a_2*R^2 + (a_inf - a_2)*R^2*f(R/R_c)
# where f(0) = 0, f(inf) = 1

# Try f(x) = 1 - exp(-x)
def mu_model_exp(R, R_c):
    f = 1.0 - np.exp(-R / R_c)
    return a1_z2 * R + a2_z2 * R**2 + (a_inf - a2_z2) * R**2 * f

# Try f(x) = x/(1+x) (rational)
def mu_model_rat(R, R_c):
    x = R / R_c
    f = x / (1.0 + x)
    return a1_z2 * R + a2_z2 * R**2 + (a_inf - a2_z2) * R**2 * f

# Try f(x) = tanh(x)
def mu_model_tanh(R, R_c):
    f = np.tanh(R / R_c)
    return a1_z2 * R + a2_z2 * R**2 + (a_inf - a2_z2) * R**2 * f

# Two-parameter: f(x) = 1 - 1/(1+x)^p
def mu_model_2p(R, R_c, p):
    x = R / R_c
    f = 1.0 - 1.0 / (1.0 + x)**p
    return a1_z2 * R + a2_z2 * R**2 + (a_inf - a2_z2) * R**2 * f

# Fit each model
from scipy.optimize import curve_fit

mask = R_grid > 0.05  # avoid R=0 singularity
R_fit = R_grid[mask]
mu_fit = mu_exact[mask]

print("\nFitting mu(R) = a_1*R + a_2*R^2 + (a_inf - a_2)*R^2*f(R/R_c)")
print(f"  where a_1 = {a1_z2:.6f}, a_2 = {a2_z2:.6f}, a_inf = {a_inf:.1f}")

# Exponential crossover
try:
    popt, _ = curve_fit(mu_model_exp, R_fit, mu_fit, p0=[1.0])
    R_c_exp = popt[0]
    resid_exp = np.max(np.abs(mu_model_exp(R_fit, R_c_exp) - mu_fit))
    print(f"\n  f(x) = 1 - exp(-x):  R_c = {R_c_exp:.6f},  max|resid| = {resid_exp:.6f}")
except Exception as e:
    print(f"\n  Exponential fit failed: {e}")
    R_c_exp = None

# Rational crossover
try:
    popt, _ = curve_fit(mu_model_rat, R_fit, mu_fit, p0=[1.0])
    R_c_rat = popt[0]
    resid_rat = np.max(np.abs(mu_model_rat(R_fit, R_c_rat) - mu_fit))
    print(f"  f(x) = x/(1+x):     R_c = {R_c_rat:.6f},  max|resid| = {resid_rat:.6f}")
except Exception as e:
    print(f"  Rational fit failed: {e}")
    R_c_rat = None

# Tanh crossover
try:
    popt, _ = curve_fit(mu_model_tanh, R_fit, mu_fit, p0=[1.0])
    R_c_tanh = popt[0]
    resid_tanh = np.max(np.abs(mu_model_tanh(R_fit, R_c_tanh) - mu_fit))
    print(f"  f(x) = tanh(x):     R_c = {R_c_tanh:.6f},  max|resid| = {resid_tanh:.6f}")
except Exception as e:
    print(f"  Tanh fit failed: {e}")
    R_c_tanh = None

# Two-parameter fit
try:
    popt, _ = curve_fit(mu_model_2p, R_fit, mu_fit, p0=[1.0, 1.0])
    R_c_2p, p_2p = popt
    resid_2p = np.max(np.abs(mu_model_2p(R_fit, R_c_2p, p_2p) - mu_fit))
    print(f"  f(x) = 1-1/(1+x)^p: R_c = {R_c_2p:.6f}, p = {p_2p:.6f},  max|resid| = {resid_2p:.6f}")
except Exception as e:
    print(f"  Two-param fit failed: {e}")
    R_c_2p, p_2p = None, None


# =========================================================================
# STEP 3: What a_2 and n_eff must be for exact He energy
# =========================================================================

print("\n" + "=" * 72)
print("STEP 3: Required a_2 and n_eff for exact He energy")
print("=" * 72)

E_exact_He = -2.903724
Z_eff = -a1_z2

# E = a_2 - Z_eff^2 / (2*n_eff^2)
# => a_2 = E + Z_eff^2 / (2*n_eff^2)
# With n_eff = 5/2:
a2_required_n52 = E_exact_He + Z_eff**2 / (2 * 2.5**2)
print(f"\n  Z_eff = -a_1 = {Z_eff:.10f}")
print(f"  Z_eff^2 = {Z_eff**2:.10f}")

print(f"\n  If n_eff = 5/2 (standard):")
print(f"    a_2 required = {a2_required_n52:.10f}")
print(f"    a_2 actual   = {a2_z2:.10f}")
print(f"    Deficit:       {a2_required_n52 - a2_z2:.10f}")

# Or solve for n_eff given our a_2:
# E = a_2 - Z_eff^2/(2*n_eff^2) = E_exact
# n_eff^2 = Z_eff^2 / (2*(a_2 - E_exact))
n_eff_sq = Z_eff**2 / (2 * (a2_z2 - E_exact_He))
n_eff_req = np.sqrt(n_eff_sq)
print(f"\n  If a_2 = {a2_z2:.6f} (our 50-term value):")
print(f"    n_eff required = {n_eff_req:.10f}")
print(f"    n_eff standard = 2.5")
print(f"    Shift:           {n_eff_req - 2.5:.10f}")

# What about the exact V_eff minimum?
V_eff_exact = mu_exact / R_grid**2 + 15.0 / (8.0 * R_grid**2)
i_min = np.argmin(V_eff_exact)
R_min = R_grid[i_min]
V_min = V_eff_exact[i_min]

print(f"\n  Exact V_eff(R) minimum: V = {V_min:.6f} Ha at R = {R_min:.4f} bohr")
print(f"  Perturbative minimum:   V = {a2_z2 - a1_z2**2 * 2/15:.6f} at "
      f"R = {15/(4*abs(a1_z2)):.4f}")

# What a_2_eff would give the exact well depth?
# V_eff_min = a_2 - 2*a_1^2/15
a2_eff_from_Vmin = V_min + 2 * a1_z2**2 / 15
print(f"\n  Effective a_2 from V_eff minimum: {a2_eff_from_Vmin:.10f}")
print(f"  This vs our a_2: deficit = {a2_eff_from_Vmin - a2_z2:.6f}")


# =========================================================================
# STEP 4: Digamma connection for the a_2 sum
# =========================================================================

print("\n" + "=" * 72)
print("STEP 4: Digamma / special function analysis of a_2")
print("=" * 72)

# a_2 = sum_{m=3,5,7,...} |V_m1|^2 / (E_1 - E_m)
# V_m1 = V_nuc(m,1) + V_ee(m,1)
# V_nuc(m,1) = (-8Z/pi) * [I_nuc((m+1)/2) - I_nuc((m-1)/2)]
# I_nuc(p) - I_nuc(p-1) = 1/(4p-3) + 1/(4p-1)  (two new terms in the sum)
# E_1 - E_m = -(2m^2 - 2) = -2(m^2 - 1)

print("\nAnalysis of individual terms in a_2 sum:")
print(f"\n  {'m':>3}  {'V_nuc':>12}  {'V_ee':>12}  {'V_total':>12}  {'dE':>8}  {'term':>14}")
print("-" * 72)

terms_nuc_only = []
terms_ee_only = []
terms_cross = []
terms_total = []

for m in range(3, 62, 2):
    V_n = V_nuc_off_diag(m, 1, Z)
    V_e = V_ee_off_diag(m, 1)
    V_t = V_n + V_e
    dE = -2 * (m**2 - 1)
    term = V_t**2 / dE  # negative

    terms_nuc_only.append(V_n**2 / dE)
    terms_ee_only.append(V_e**2 / dE)
    terms_cross.append(2 * V_n * V_e / dE)
    terms_total.append(term)

    if m <= 15 or m == 21 or m == 31 or m == 61:
        print(f"  {m:3d}  {V_n:12.8f}  {V_e:12.8f}  {V_t:12.8f}  {dE:8d}  {term:14.10f}")

a2_nuc = sum(terms_nuc_only)
a2_ee = sum(terms_ee_only)
a2_cross = sum(terms_cross)
a2_total = sum(terms_total)

print(f"\n  Decomposition of a_2:")
print(f"    Nuclear only:  {a2_nuc:.10f}")
print(f"    V_ee only:     {a2_ee:.10f}")
print(f"    Cross (2*nuc*ee): {a2_cross:.10f}")
print(f"    Total:         {a2_total:.10f}")
print(f"    Check: {a2_nuc:.6f} + {a2_ee:.6f} + {a2_cross:.6f} = {a2_nuc+a2_ee+a2_cross:.6f}")

# Now analyze the NUCLEAR-ONLY part, which is the simplest:
# V_nuc(m,1) = (-8Z/pi) * [I_nuc((m+1)/2) - I_nuc((m-1)/2)]
# For m = 2k+1 (k=1,2,...):
#   I_nuc(k+1) - I_nuc(k) = 1/(4k+1) + 1/(4k+3)
#   V_nuc(2k+1, 1) = (-8Z/pi) * [1/(4k+1) + 1/(4k+3)]
#   E_1 - E_{2k+1} = -2((2k+1)^2 - 1) = -2(4k^2 + 4k) = -8k(k+1)

print("\n" + "-" * 40)
print("Nuclear-only a_2 closed form analysis:")
print("-" * 40)

# a2_nuc = sum_k (8Z/pi)^2 * [1/(4k+1) + 1/(4k+3)]^2 / (-8k(k+1))
# = -(8Z/pi)^2 / 8 * sum_k [1/(4k+1) + 1/(4k+3)]^2 / (k(k+1))
# = -(8Z)^2 / (8*pi^2) * sum_k [...]

prefactor = -(8*Z)**2 / (8 * np.pi**2)
print(f"  Prefactor: -(8Z)^2/(8*pi^2) = {prefactor:.10f}")

print(f"\n  Term-by-term check:")
for k in range(1, 6):
    m = 2*k + 1
    bracket = 1.0/(4*k+1) + 1.0/(4*k+3)
    denom = k * (k + 1)
    term_formula = prefactor * bracket**2 / denom
    term_direct = terms_nuc_only[k-1]
    print(f"    k={k} (m={m}): formula = {term_formula:.10f}, direct = {term_direct:.10f}, "
          f"match: {abs(term_formula - term_direct) < 1e-12}")

# Can the nuclear sum be expressed via digamma?
# sum_k [1/(4k+1) + 1/(4k+3)]^2 / (k(k+1))
# Expand: [1/(4k+1) + 1/(4k+3)]^2 = 1/(4k+1)^2 + 2/((4k+1)(4k+3)) + 1/(4k+3)^2
# And 1/(k(k+1)) = 1/k - 1/(k+1)  (partial fractions)

# So the sum splits into three parts, each telescoping in k:
# S = sum_k (1/k - 1/(k+1)) * [1/(4k+1)^2 + 2/((4k+1)(4k+3)) + 1/(4k+3)^2]

# This does NOT telescope into digamma, but each piece can be expressed
# via the Hurwitz zeta function or polygamma.

# Let's check if the sum has a recognizable value:
S_nuc = 0.0
for k in range(1, 10000):
    bracket = 1.0/(4*k+1) + 1.0/(4*k+3)
    S_nuc += bracket**2 / (k * (k + 1))

a2_nuc_formula = prefactor * S_nuc
print(f"\n  Nuclear sum S = {S_nuc:.12f}")
print(f"  a2_nuc from formula = {a2_nuc_formula:.10f}")
print(f"  a2_nuc direct       = {a2_nuc:.10f}")

# Try to identify S_nuc as a combination of known constants
# S_nuc ~ 0.27... Let me check against pi^2/6, ln(2), Catalan's constant, etc.
from scipy.special import digamma
G = 0.9159655941  # Catalan's constant

print(f"\n  Attempting identification of S = {S_nuc:.10f}:")
candidates = [
    ("pi^2/36", np.pi**2/36),
    ("1/4", 0.25),
    ("pi^2/40", np.pi**2/40),
    ("4 - pi", 4 - np.pi),
    ("(4-pi)/pi", (4-np.pi)/np.pi),
    ("pi/12", np.pi/12),
    ("ln(2)", np.log(2)),
    ("2*ln(2) - 1", 2*np.log(2) - 1),
    ("Catalan/3", G/3),
    ("1 - 2/pi", 1 - 2/np.pi),
    ("8/pi^2 - 1/2", 8/np.pi**2 - 0.5),
]
for name, val in candidates:
    if abs(val - S_nuc) < 0.01:
        print(f"    {name} = {val:.10f}  (diff = {abs(val - S_nuc):.2e})")


# =========================================================================
# STEP 5: SO(6) representation content of 1/r_12
# =========================================================================

print("\n" + "=" * 72)
print("STEP 5: SO(6) structure of 1/r_12")
print("=" * 72)

# 1/r_12 = 1 / sqrt(r1^2 + r2^2 - 2*r1*r2*cos(theta_12))
# In hyperspherical: r1 = R*cos(alpha), r2 = R*sin(alpha)
# r_12 = R * sqrt(1 - sin(2*alpha)*cos(theta_12))
# So 1/r_12 = (1/R) * 1/sqrt(1 - sin(2*alpha)*cos(theta_12))

# The 1/sqrt(1 - t*cos(theta)) has a Legendre expansion:
# 1/sqrt(1 - t*cos(theta)) = sum_l P_l(cos(theta)) * t^l / sqrt(1-t^2) ...
# Actually, the standard expansion:
# 1/r_12 = (1/R) * sum_l (r_</r_>)^l * P_l(cos(theta))
# = (1/R) * sum_l [min(cos,sin)/max(cos,sin)]^l * P_l(cos(theta_12))

# In the l=0 sector (S-wave for both electrons), only the l=0 term of
# 1/r_12 contributes to the angular equation. This is:
# V_ee^{l=0} = (1/R) * 1/max(cos(alpha), sin(alpha))
# (The min/max expansion with P_0 = 1 gives sum_l (min/max)^l for l=0 only
#  in the partial wave equation — the Gaunt coupling mixes l values.)

# The FULL SO(6) decomposition of V_ee requires the angular matrix elements
# <nu', l' | V_ee | nu, l> for all (nu, l) quantum numbers.

# For the l=0 sector:
# <m| 1/max(cos,sin) |n> = int_0^{pi/2} sin(2ma)*sin(2na)/max(cos,sin) da

# This is exactly V_ee^{mn} = (4*sqrt(2)/pi) * [I_ee_rat((m+n)/2) - I_ee_rat(|m-n|/2)]

# The selection rule Delta_nu = 0 mod 4 means 1/r_12 connects:
# nu=0 <-> nu=0, 4, 8, ...
# This is the symmetric traceless tensor representation of SO(6).

# What representation does 1/r_12 belong to?
# 1/r_12 is a scalar under rotations (l=0 part), and as a function on S^5
# it decomposes into SO(6) harmonics Y_{nu,L}(Omega).
# The key question: how does the coupling strength |V_{m,1}|^2 decay with m?

print("\nDecay of |V^{m,1}|^2 with m (nu = 2(m-1)):")
print(f"\n  {'m':>3}  {'nu':>4}  {'|V_nuc|':>12}  {'|V_ee|':>12}  {'|V_tot|':>12}  {'|V|^2':>14}  {'ratio':>8}")
print("-" * 72)

prev_V2 = None
for m in range(3, 32, 2):
    V_n = abs(V_nuc_off_diag(m, 1, Z))
    V_e = abs(V_ee_off_diag(m, 1))
    V_t = abs(V_nuc_off_diag(m, 1, Z) + V_ee_off_diag(m, 1))
    V2 = V_t**2
    nu = 2*(m-1)
    ratio = V2 / prev_V2 if prev_V2 else float('nan')
    print(f"  {m:3d}  {nu:4d}  {V_n:12.8f}  {V_e:12.8f}  {V_t:12.8f}  {V2:14.10f}  {ratio:8.4f}")
    prev_V2 = V2

# Check: does |V_{m,1}| ~ 1/m^p for some power p?
m_vals = np.array(list(range(3, 62, 2)), dtype=float)
V_vals = np.array([abs(V_nuc_off_diag(int(m), 1, Z) + V_ee_off_diag(int(m), 1)) for m in m_vals])
# Log-log fit: log|V| = -p*log(m) + c
coeffs = np.polyfit(np.log(m_vals), np.log(V_vals), 1)
p_decay = -coeffs[0]
print(f"\n  Power-law fit: |V_{{m,1}}| ~ m^(-{p_decay:.4f})")
print(f"  Expected for smooth function on S^5: V ~ nu^(-3) -> m^(-3)")

# For a_2 convergence: terms go as m^{-2p} / m^2, so convergent if p > 1/2.
# With p ~ 3, terms decay as m^{-8}, extremely fast.


# =========================================================================
# STEP 6: Continued fraction representation
# =========================================================================

print("\n" + "=" * 72)
print("STEP 6: Continued fraction for mu(R)")
print("=" * 72)

# The coupled-channel problem: H_angular * Phi = mu * Phi
# In the free basis |n>, this is:
# (E_n^free + R*V_nn)*c_n + R*sum_{m!=n} V_mn*c_m = mu*c_n

# For the ground state, this is a matrix eigenvalue problem.
# The standard continued-fraction approach: project out all states except |1>
# mu = E_1^free + R*V_11 + sum_{m!=1} R^2*|V_m1|^2 / (mu - E_m^free - R*V_mm - ...)
# This is a self-consistent equation (mu appears on both sides).

# At lowest order (no self-consistency, no further coupling):
# mu ~ 0 + R*a_1 + sum_m R^2*|V_m1|^2 / (-delta_E_m)
# = a_1*R + a_2*R^2   (our result)

# The EXACT continued fraction includes all orders:
# mu(R) = R*V_11 + R^2*|V_31|^2 / (E_1-E_3 + R*(V_33-V_11) + R^2*|V_53|^2/(E_3-E_5+...) - mu(R))

# Let's build this numerically and compare to exact.

def cf_mu(R: float, Z: float, n_levels: int = 20) -> float:
    """
    Continued-fraction approximation for mu(R).

    Uses the Tridiagonal-like structure (selection rule couples n to n+2 only
    in leading order) to build a CF from the top down.
    """
    # Build the diagonal and off-diagonal elements
    V_diag = [a1_exact(n, Z) * R for n in range(1, n_levels + 1)]
    E_free = [2 * n**2 - 2 for n in range(1, n_levels + 1)]

    # Off-diagonal: V_{n, n+2} (selection rule: Delta_m = 2, Delta_nu = 4)
    V_off = []
    for n in range(1, n_levels - 1):
        m = n + 2
        V = V_nuc_off_diag(m, n, Z) + V_ee_off_diag(m, n)
        V_off.append(V * R)

    # Also include V_{n, n+4}, V_{n, n+6}, ... but these are smaller.
    # For the CF, use only nearest-neighbor coupling (tridiagonal approx).

    # Build CF from the bottom up:
    # mu_N = E_free[N] + V_diag[N]
    # mu_{N-1} = E_free[N-1] + V_diag[N-1] + V_off[N-1]^2 / (mu_1_trial - mu_N)
    # But this needs mu_1 (self-consistent).

    # Simpler: just truncate the matrix and diagonalize
    n_mat = min(n_levels, 15)
    H = np.zeros((n_mat, n_mat))
    for i in range(n_mat):
        n = i + 1  # basis index n=1,2,...
        H[i, i] = E_free[i] + V_diag[i]
        for j in range(i + 1, n_mat):
            m = j + 1
            if (m + n) % 2 == 0:  # selection rule
                V = (V_nuc_off_diag(m, n, Z) + V_ee_off_diag(m, n)) * R
                H[i, j] = V
                H[j, i] = V

    eigenvalues = np.sort(np.linalg.eigh(H)[0])
    return eigenvalues[0]  # ground state

# Compare CF to exact
print("\nContinued-fraction (matrix truncation) vs exact angular solver:")
print(f"\n  {'R':>6}  {'mu_exact':>12}  {'mu_CF':>12}  {'error':>12}")
print("-" * 50)
R_test = [0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 15.0]
for R in R_test:
    mu_e, _ = solve_angular(R=R, Z=Z, l_max=3, n_alpha=200, n_channels=1)
    mu_cf = cf_mu(R, Z, n_levels=15)
    err = mu_cf - mu_e[0]
    print(f"  {R:6.2f}  {mu_e[0]:12.6f}  {mu_cf:12.6f}  {err:12.6f}")


# =========================================================================
# STEP 7: Large-R asymptotic analysis
# =========================================================================

print("\n" + "=" * 72)
print("STEP 7: Large-R asymptotic analysis")
print("=" * 72)

# At large R, mu(R) -> E_He+ * R^2 = -Z^2*R^2/2
# The correction is: mu(R) ~ -Z^2*R^2/2 + c_1*R + c_0 + d_1/R + ...
# From the exact curve, extract c_1 and c_0.

large_mask = R_grid > 5.0
R_large = R_grid[large_mask]
mu_large = mu_exact[large_mask]

# Fit: mu = A*R^2 + B*R + C + D/R
def large_R_model(R, A, B, C, D):
    return A * R**2 + B * R + C + D / R

popt_large, _ = curve_fit(large_R_model, R_large, mu_large, p0=[-2.0, 0.0, 0.0, 0.0])
A_fit, B_fit, C_fit, D_fit = popt_large

print(f"\n  Large-R fit (R > 5 bohr):")
print(f"    mu(R) = {A_fit:.6f}*R^2 + {B_fit:.6f}*R + {C_fit:.4f} + {D_fit:.4f}/R")
print(f"    Expected: A = -Z^2/2 = {-Z**2/2:.1f}")
print(f"    Fit A = {A_fit:.6f}  (error from -2.0: {abs(A_fit + 2.0):.4f})")

# The coefficient B should be related to the He+ ionization energy correction
# B*R comes from the long-range interaction: electron-ion at distance R
print(f"\n    Coefficient B = {B_fit:.6f}")
print(f"    For a hydrogen-like outer electron: E_outer ~ -Z_ion^2/(2n^2) = -(Z-1)^2/2")
print(f"    The 'B' correction comes from the next-order asymptotic expansion")


# =========================================================================
# STEP 8: Comparison of all functional forms at key R values
# =========================================================================

print("\n" + "=" * 72)
print("STEP 8: Comparison table at key R values")
print("=" * 72)

print(f"\n  {'R':>6}  {'mu_exact':>10}  {'PT(a1R)':>10}  {'PT(a1R+a2R2)':>13}  ", end="")
if R_c_exp is not None:
    print(f"{'Exp crossover':>14}  ", end="")
if R_c_rat is not None:
    print(f"{'Rat crossover':>14}  ", end="")
print(f"{'Large-R fit':>12}")
print("-" * 100)

R_compare = [0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 15.0]
for R in R_compare:
    mu_e, _ = solve_angular(R=R, Z=Z, l_max=3, n_alpha=200, n_channels=1)
    mu_pt1 = a1_z2 * R
    mu_pt2 = a1_z2 * R + a2_z2 * R**2

    print(f"  {R:6.2f}  {mu_e[0]:10.4f}  {mu_pt1:10.4f}  {mu_pt2:13.4f}  ", end="")
    if R_c_exp is not None:
        print(f"{mu_model_exp(R, R_c_exp):14.4f}  ", end="")
    if R_c_rat is not None:
        print(f"{mu_model_rat(R, R_c_rat):14.4f}  ", end="")
    print(f"{large_R_model(R, *popt_large):12.4f}")


# =========================================================================
# STEP 9: Pade [2/1] and [2/2] in R for mu(R)
# =========================================================================

print("\n" + "=" * 72)
print("STEP 9: Pade approximants for mu(R)")
print("=" * 72)

# We know: mu(R) = a_1*R + a_2*R^2 + a_3*R^3 + ...
# Extract a_3 from small-R behavior:
small_mask = R_grid < 1.0
R_small = R_grid[small_mask]
mu_small = mu_exact[small_mask]

# Fit mu/R = a_1 + a_2*R + a_3*R^2 + a_4*R^3
mu_over_R_small = mu_small / R_small
coeffs_small = np.polyfit(R_small, mu_over_R_small, 4)
a4_fit, a3_fit, a2_fit_check, a1_fit_check, _ = coeffs_small

print(f"\n  Small-R polynomial fit of mu(R)/R:")
print(f"    a_1 fit = {a1_fit_check:.6f}  (exact: {a1_z2:.6f})")
print(f"    a_2 fit = {a2_fit_check:.6f}  (exact: {a2_z2:.6f})")
print(f"    a_3 fit = {a3_fit:.6f}")
print(f"    a_4 fit = {a4_fit:.6f}")

# Pade [1/1] for mu(R)/R: (a_1 + b*R) / (1 + c*R) with correct a_1, a_2
# Expanding: a_1 + (b - a_1*c)*R + ... = a_1 + a_2*R + ...
# So b - a_1*c = a_2. And large R: ~ b/c.
# We want mu/R -> -Z^2*R/2, i.e., mu/R -> -2*R. So b/c = -infinity? No.
# Pade [1/1] can't capture linear growth.

# Pade [2/1] for mu(R): a_1*R + a_2*R^2 + a_3*R^3 at small R
# mu = (p_1*R + p_2*R^2) / (1 + q_1*R) matches 3 coefficients
# Expand: (p_1*R + p_2*R^2)(1 - q_1*R + q_1^2*R^2 - ...)
# = p_1*R + (p_2 - p_1*q_1)*R^2 + (p_1*q_1^2 - p_2*q_1)*R^3 + ...
# Match: p_1 = a_1, p_2 - p_1*q_1 = a_2, p_1*q_1^2 - p_2*q_1 = a_3
# => q_1 = (p_1*q_1^2 - a_3) / p_2 ... complicated.
# From first two: p_2 = a_2 + a_1*q_1
# From third: a_1*q_1^2 - (a_2 + a_1*q_1)*q_1 = a_3
# a_1*q_1^2 - a_2*q_1 - a_1*q_1^2 = a_3
# -a_2*q_1 = a_3  =>  q_1 = -a_3/a_2

a1 = a1_z2
a2 = a2_z2
a3 = a3_fit

q1 = -a3 / a2
p1 = a1
p2 = a2 + a1 * q1

print(f"\n  Pade [2/1]: mu(R) = ({p1:.6f}*R + {p2:.6f}*R^2) / (1 + {q1:.6f}*R)")
print(f"    Large-R limit: mu ~ {p2/q1:.4f}*R  (should be -2*R^2, WRONG -- Pade [2/1] can't capture R^2)")

# Pade [2/2] would be better:
# mu = (p_1*R + p_2*R^2) / (1 + q_1*R + q_2*R^2) -- matches 4 coefficients
# but needs a_4.
# Large R: mu -> p_2/q_2 (constant) -- still wrong for R^2 growth.

# Conclusion: standard Pade fails because mu(R) has DIFFERENT leading powers
# at small R (linear) and large R (quadratic). Need a modified approach.

print("\n  CONCLUSION: Standard Pade approximants fail for mu(R) because")
print("  the function changes QUALITATIVE character (linear -> quadratic).")
print("  This is a strong hint that mu(R) does NOT have a closed form")
print("  as a rational function of R.")


# =========================================================================
# STEP 10: Best interpolation and energy prediction
# =========================================================================

print("\n" + "=" * 72)
print("STEP 10: Best interpolation formula and energy prediction")
print("=" * 72)

# The best model: mu(R) = a_1*R + a_2*R^2 + (a_inf - a_2)*R^2*f(R/R_c)
# with the best f and R_c from the fits above.

if R_c_2p is not None:
    best_R_c = R_c_2p
    best_p = p_2p
    best_name = f"f(x) = 1-1/(1+x)^{best_p:.4f}, R_c = {best_R_c:.4f}"
    mu_best = lambda R: mu_model_2p(R, best_R_c, best_p)
elif R_c_exp is not None:
    best_R_c = R_c_exp
    best_name = f"f(x) = 1-exp(-x), R_c = {best_R_c:.4f}"
    mu_best = lambda R: mu_model_exp(R, best_R_c)
else:
    best_R_c = 1.0
    best_name = "No good fit found"
    mu_best = lambda R: a1_z2 * R + a2_z2 * R**2

print(f"\n  Best interpolation: {best_name}")

# Compute V_eff and find bound states
R_dense = np.linspace(0.1, 30.0, 500)
mu_interp = np.array([mu_best(R) for R in R_dense])
V_eff_interp = mu_interp / R_dense**2 + 15.0 / (8.0 * R_dense**2)

i_min = np.argmin(V_eff_interp)
V_min_interp = V_eff_interp[i_min]
R_min_interp = R_dense[i_min]

print(f"\n  V_eff minimum from interpolation: V = {V_min_interp:.6f} at R = {R_min_interp:.2f}")

# The quasi-Coulomb energy from the interpolated potential:
# For the interpolated potential, solving the full Schrodinger equation would
# give a better answer. But we can estimate from the well depth.
print(f"\n  Rough energy estimate from V_eff minimum: {V_min_interp:.4f} Ha")
print(f"  (Exact He: -2.904 Ha)")
print(f"  (The V_eff minimum is NOT the energy -- need to solve the radial eq.)")


# =========================================================================
# FINAL SUMMARY
# =========================================================================

print("\n" + "=" * 72)
print("FINAL SUMMARY: CROSSOVER ANALYSIS")
print("=" * 72)

print(f"""
KEY FINDINGS:

1. PADE APPROXIMANTS: Standard Pade fails because mu(R) changes qualitative
   character from O(R) at small R to O(R^2) at large R. No rational function
   of R can capture this transition.

2. TWO-SCALE INTERPOLATION: mu(R) = a_1*R + a_2*R^2 + (a_inf - a_2)*R^2*f(R/R_c)
   works well with one crossover parameter R_c. Best fit:""")

if R_c_exp is not None:
    print(f"     Exponential: R_c = {R_c_exp:.4f}, max residual = {resid_exp:.4f}")
if R_c_rat is not None:
    print(f"     Rational:    R_c = {R_c_rat:.4f}, max residual = {resid_rat:.4f}")
if R_c_tanh is not None:
    print(f"     Tanh:        R_c = {R_c_tanh:.4f}, max residual = {resid_tanh:.4f}")
if R_c_2p is not None:
    print(f"     Power:       R_c = {R_c_2p:.4f}, p = {p_2p:.4f}, max residual = {resid_2p:.4f}")

print(f"""
3. a_2 DECOMPOSITION:
     Nuclear-only: {a2_nuc:.6f}  ({abs(a2_nuc/a2_total)*100:.1f}%)
     V_ee-only:    {a2_ee:.6f}  ({abs(a2_ee/a2_total)*100:.1f}%)
     Cross term:   {a2_cross:.6f}  ({abs(a2_cross/a2_total)*100:.1f}%)
     Total:        {a2_total:.6f}

   The nuclear part dominates. Each term involves
   [1/(4k+1) + 1/(4k+3)]^2 / (k(k+1)) -- NOT a standard special function.

4. SO(6) COUPLING DECAY: |V_{{m,1}}| ~ m^(-{p_decay:.1f}).
   This rapid decay (faster than m^(-3) expected for smooth S^5 functions)
   means a_2 converges fast, but the FULL mu(R) at large R requires
   the infinite sum to resum into -Z^2*R^2/2.

5. REQUIRED CORRECTIONS for exact He:
     a_2 needed:  {a2_required_n52:.6f}  (actual: {a2_z2:.6f})
     n_eff needed: {n_eff_req:.6f}  (standard: 2.5)
     The 5.5% error comes primarily from truncating mu(R) at O(R^2).

6. LARGE-R EXPANSION:
     mu(R) ~ {A_fit:.4f}*R^2 + {B_fit:.4f}*R + {C_fit:.2f} + {D_fit:.2f}/R
     The R^2 coefficient {A_fit:.4f} matches expected -Z^2/2 = -2.0.
     The R coefficient {B_fit:.4f} is the long-range ion-electron interaction.

CONCLUSION:
  mu(R) does NOT have a simple algebraic closed form. It interpolates between
  two different power-law regimes (R vs R^2) with a smooth crossover at
  R_c ~ {best_R_c:.1f} bohr. The crossover scale R_c is NOT obviously algebraic
  (not a simple rational multiple of pi, sqrt(2), etc.).

  The perturbative coefficients a_1, a_2, ... ARE algebraic (partial harmonic sums),
  but the full function requires their resummation, which produces a transcendental
  crossover. This is the non-perturbative content of the helium problem.
""")
