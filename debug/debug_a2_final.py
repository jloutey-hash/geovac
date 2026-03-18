"""
Final computation of a_2 using the exact sympy results already obtained.

From the sympy computation:
  <3|V_nuc|1> = 24/35  (I_nuc)
  <3|V_ee|1>  = -2*sqrt(2)/35  (I_ee)
  <5|V_nuc|1> = 40/99  (I_nuc)
  <5|V_ee|1>  = 2*sqrt(2)/99  (I_ee)

Full matrix element: V_{mn} = (4/pi) * (-Z * I_nuc + I_ee)

a_2(n=1) = |V_{31}|^2 / (E_1 - E_3) + |V_{51}|^2 / (E_1 - E_5) + ...
"""

import numpy as np
import sympy as sp
from sympy import pi, sqrt, Rational, simplify

Z_s = sp.Symbol('Z', positive=True)

print("=" * 72)
print("EXACT a_2(nu=0) from sympy off-diagonal elements")
print("=" * 72)

# Free energies: E_n = 2n^2 - 2
E = lambda n: 2 * n**2 - 2

# Matrix elements V_{m,1} = (4/pi) * (-Z * I_nuc(m,1) + I_ee(m,1))
# From sympy:
V_31_expr = Rational(4, 1) / pi * (-Z_s * Rational(24, 35) + (-2 * sqrt(2) / 35))
V_51_expr = Rational(4, 1) / pi * (-Z_s * Rational(40, 99) + (2 * sqrt(2) / 99))

V_31_expr = simplify(V_31_expr)
V_51_expr = simplify(V_51_expr)

print(f"\n  V_31(Z) = {V_31_expr}")
print(f"  V_51(Z) = {V_51_expr}")

V_31_z2 = float(V_31_expr.subs(Z_s, 2))
V_51_z2 = float(V_51_expr.subs(Z_s, 2))
print(f"\n  V_31(Z=2) = {V_31_z2:.10f}")
print(f"  V_51(Z=2) = {V_51_z2:.10f}")

# a_2 terms
# E_1 = 0, E_3 = 16, E_5 = 48
term_31 = V_31_expr**2 / (E(1) - E(3))  # 0 - 16 = -16
term_51 = V_51_expr**2 / (E(1) - E(5))  # 0 - 48 = -48

a2_two_terms = simplify(term_31 + term_51)
a2_two_z2 = float(a2_two_terms.subs(Z_s, 2))

print(f"\n  a_2 contribution from m=3: {float(term_31.subs(Z_s, 2)):.10f}")
print(f"  a_2 contribution from m=5: {float(term_51.subs(Z_s, 2)):.10f}")
print(f"  a_2 (two terms):           {a2_two_z2:.10f}")

# Expand in Z
a2_expanded = sp.expand(a2_two_terms)
poly = sp.Poly(a2_expanded, Z_s)
print(f"\n  Z-structure of a_2 (two-term approximation):")
for power, coeff in sorted(poly.as_dict().items(), reverse=True):
    c = simplify(coeff)
    print(f"    Z^{power[0]}: {c} = {float(c):.10f}")

# Numerical comparison with 20-term sum from debug_algebraic_coefficients.py
a2_20term = -0.24387012  # from previous run
print(f"\n  Comparison:")
print(f"    a_2 (2 terms, exact):  {a2_two_z2:.10f}")
print(f"    a_2 (20 terms, FD):    {a2_20term:.10f}")
print(f"    Missing terms:         {a2_20term - a2_two_z2:.10f}")
print(f"    Fraction captured:     {a2_two_z2 / a2_20term * 100:.1f}%")

# =============================================================================
# Now put it all together: the quasi-Coulomb hyperradial equation
# =============================================================================

print("\n" + "=" * 72)
print("QUASI-COULOMB HYPERRADIAL EQUATION")
print("=" * 72)

# At nu=0: mu(R) ~ a_1*R + a_2*R^2
# V_eff(R) = mu(R)/R^2 + 15/(8R^2)
#           = a_1/R + a_2 + 15/(8R^2)
# = (15/8)/R^2 + a_1/R + a_2

# This is: -1/2 F'' + [(15/8)/R^2 + a_1/R + a_2] F = E F
# With l_eff(l_eff+1)/2 = 15/8 => l_eff = 3/2 (exact)

# Substituting E' = E - a_2:
# -1/2 F'' + [(15/8)/R^2 + a_1/R] F = E' F
# This is the radial hydrogen equation with Z_eff = -a_1, l = 3/2
# E'_N = -Z_eff^2 / (2*n_eff^2), n_eff = N + l_eff + 1 = N + 5/2

a1_exact = sp.Rational(8, 1) * (-4 * Z_s + sqrt(2)) / (3 * pi)

print(f"\n  a_1 = {a1_exact}")
print(f"  a_2 = {a2_two_terms}  (leading two terms)")

Z_eff = -a1_exact
n_eff = sp.Rational(5, 2)  # N=0 ground state

# E_gs = a_2 - Z_eff^2 / (2 * (5/2)^2) = a_2 - 2*Z_eff^2/25
E_gs = a2_two_terms - Z_eff**2 / (2 * n_eff**2)
E_gs_simplified = simplify(E_gs)

print(f"\n  Exact ground state energy formula (N=0):")
print(f"  E = a_2 - Z_eff^2 / (2*(5/2)^2)")
print(f"    = a_2 - 2*a_1^2/25")

E_gs_z2 = float(E_gs_simplified.subs(Z_s, 2))
print(f"\n  E_gs(Z=2) = {E_gs_z2:.6f} Ha")
print(f"  Exact He:   -2.903724 Ha")
print(f"  Error:       {abs(E_gs_z2 + 2.903724) / 2.903724 * 100:.2f}%")

# With 20-term a_2
a2_full = -0.24387012
a1_z2 = float(a1_exact.subs(Z_s, 2))
Z_eff_z2 = -a1_z2
E_gs_full = a2_full - Z_eff_z2**2 / (2 * 2.5**2)
print(f"\n  With 20-term a_2:")
print(f"  E_gs(Z=2) = {E_gs_full:.6f} Ha")
print(f"  Error:       {abs(E_gs_full + 2.903724) / 2.903724 * 100:.2f}%")

# =============================================================================
# Test: what about the exact (numeric) V_eff minimum?
# =============================================================================

print("\n" + "=" * 72)
print("EXACT V_eff ANALYSIS (numerical)")
print("=" * 72)

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from geovac.hyperspherical_angular import solve_angular

R_grid = np.linspace(0.1, 15.0, 300)
mu_exact = np.zeros(len(R_grid))
for i, R in enumerate(R_grid):
    mu, _ = solve_angular(R=R, Z=2.0, l_max=3, n_alpha=200, n_channels=1)
    mu_exact[i] = mu[0]

V_eff_exact = mu_exact / R_grid**2 + 15.0 / (8.0 * R_grid**2)
V_eff_pt1 = a1_z2 / R_grid + 15.0 / (8.0 * R_grid**2)
V_eff_pt2 = a1_z2 / R_grid + a2_full + 15.0 / (8.0 * R_grid**2)

# Find minimum of exact V_eff
i_min = np.argmin(V_eff_exact)
R_min = R_grid[i_min]
V_min = V_eff_exact[i_min]
print(f"\n  Exact V_eff minimum: V = {V_min:.6f} Ha at R = {R_min:.2f} bohr")

# Compare perturbative V_eff at the minimum
V_pt1_at_min = a1_z2 / R_min + 15.0 / (8.0 * R_min**2)
V_pt2_at_min = a1_z2 / R_min + a2_full + 15.0 / (8.0 * R_min**2)
print(f"  PT1 V_eff at R_min:   {V_pt1_at_min:.6f} Ha (error: {abs(V_pt1_at_min - V_min):.4f})")
print(f"  PT2 V_eff at R_min:   {V_pt2_at_min:.6f} Ha (error: {abs(V_pt2_at_min - V_min):.4f})")

# The perturbative minimum
# V_pt = 15/(8R^2) + a1/R + a2
# dV/dR = -15/(4R^3) - a1/R^2 = 0
# a1*R = -15/4 => R_min_pt = -15/(4*a1) = 15/(4*|a1|)
R_min_pt = 15.0 / (4.0 * abs(a1_z2))
V_min_pt = a1_z2 / R_min_pt + a2_full + 15.0 / (8.0 * R_min_pt**2)
print(f"\n  Perturbative V_eff minimum:")
print(f"    R_min_pt = 15/(4*|a_1|) = {R_min_pt:.4f} bohr")
print(f"    V_min_pt = {V_min_pt:.6f} Ha")

# For comparison: V_eff at large R
print(f"\n  Asymptotic V_eff(R=15): exact = {V_eff_exact[-1]:.4f}, "
      f"He+ threshold = -2.0000")

# =============================================================================
# WHY THE PERTURBATIVE APPROACH MISSES: mu(R) grows as ~R^2
# =============================================================================

print("\n" + "=" * 72)
print("WHY PERTURBATION THEORY FAILS AT LARGE R")
print("=" * 72)

# At large R, mu(R) -> -Z^2 * R^2 / 2 (quadratic!)
# The perturbation series mu = a_1*R + a_2*R^2 + ... must have
# a_2 -> -Z^2/2 = -2 for He. But our perturbative a_2 is only -0.244.
# Higher order terms contribute the rest of the -R^2 coefficient.

# Fit mu(R) at large R to extract the true quadratic coefficient
large_mask = R_grid > 8.0
coeffs_large = np.polyfit(R_grid[large_mask], mu_exact[large_mask], 2)
print(f"\n  Large-R fit: mu(R) = {coeffs_large[0]:.4f}*R^2 + {coeffs_large[1]:.4f}*R + {coeffs_large[2]:.4f}")
print(f"  Expected: mu(R) -> -Z^2*R^2/2 = -2*R^2  (He+ threshold)")
print(f"  Our perturbative a_2 = {a2_full:.4f}  (should converge to -2.0)")
print(f"  This means higher-order terms contribute {-2.0 - a2_full:.4f} to the R^2 coefficient")

# The exact V_eff(R) = mu(R)/R^2 + 15/(8R^2) at large R:
# ~ -Z^2/2 + (correction)/R + 15/(8R^2)
# This is the Coulomb field of the He+ ion on the outer electron.

print(f"""
PHYSICAL INTERPRETATION:
  - At small R: mu(R) ~ a_1*R + a_2*R^2 (perturbative, algebraic)
  - At large R: mu(R) ~ -Z^2*R^2/2 + ... (hydrogenic, asymptotic)

  The perturbation series captures the small-R regime where both electrons
  are close to the nucleus (correlation-dominated). The large-R regime
  (one electron bound, one far away) requires non-perturbative physics.

  The 5.5% error comes from using the small-R Coulomb form in a region
  where the potential is qualitatively different.
""")

# =============================================================================
# FINAL: Table of all algebraic results
# =============================================================================

print("=" * 72)
print("COMPLETE TABLE OF ALGEBRAIC RESULTS")
print("=" * 72)

print("""
EXACT CLOSED FORMS
------------------
Free eigenvalues:  mu_free(nu) = nu*(nu+4)/2  [SO(6) Casimir]

First-order:  a_1(n, Z) = (4/pi) * (-2Z * I_nuc(n) + I_ee(n))

  n=1 (nu=0): a_1 = 8*(sqrt(2) - 4Z)/(3*pi)
  n=2 (nu=2): a_1 = 128*(2*sqrt(2) - 11Z)/(105*pi)
  n=3 (nu=4): a_1 = 8*(1091*sqrt(2) - 6508Z)/(3465*pi)
  n=4 (nu=6): a_1 = 512*(218*sqrt(2) - 1423Z)/(45045*pi)

  General: I_nuc(n) in Q, I_ee(n) in Q*sqrt(2)

Selection rule:  <m|V|n> = 0  unless m+n even  (Dnu = 0 mod 4)

Off-diagonal:
  V_31 = 8*(-12Z - sqrt(2))/(35*pi)
  V_51 = 8*(-20Z + sqrt(2))/(99*pi)

ALGEBRAIC BUT NOT SIMPLE CLOSED FORM
-------------------------------------
a_2(nu=0, Z) = sum_{m=3,5,7,...} |V_{m1}|^2 / (E_1 - E_m)
             ~ -0.234 (2 terms) vs -0.244 (20 terms)
             Quadratic in Z: c_2*Z^2 + c_1*Z*sqrt(2) + c_0
             Converges, but no simple closed form.

QUASI-COULOMB HYPERRADIAL EQUATION
-----------------------------------
V_eff(R) = (15/8)/R^2 + a_1/R + a_2   [truncated at O(R^0)]

l_eff = 3/2  (exact)
Z_eff = -a_1 = 8*(4Z - sqrt(2))/(3*pi)
n_eff = N + 5/2   (N = 0, 1, 2, ...)

E_N = a_2 - Z_eff^2 / (2*(N + 5/2)^2)

Ground state (N=0, Z=2):
  E = -2.744 Ha  (5.5% error vs exact -2.904 Ha)

Accuracy limited by: mu(R) ~ -Z^2*R^2/2 at large R,
not captured by finite-order perturbation theory.
""")
