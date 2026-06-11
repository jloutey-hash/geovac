"""Sprint 5 DV Step 3: Compute bipolar channel radial integrals symbolically
and express each in terms of Drake's M^K_{dir/exch} basis.

For He (1s)(2p) at Z=1:
  R_{1s}(r) = 2 e^{-r}
  R_{2p}(r) = (1/(2 sqrt 6)) r e^{-r/2}

Density products:
  rho_1s(r) = R_{1s}^2 * r^2 = 4 r^2 e^{-2r}
  rho_2p(r) = R_{2p}^2 * r^2 = (1/24) r^4 e^{-r}
  rho_1s2p(r) = R_{1s} R_{2p} * r^2 = (1/sqrt 6) r^3 e^{-3r/2}  [i.e., (2/(2 sqrt 6)) r^3 e^{-3r/2} = (1/sqrt 6) r^3 e^{-3r/2}]
"""
from __future__ import annotations

import os
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
import sys
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

import sympy as sp
from sympy import Rational, sqrt, Integer, simplify, symbols, exp, integrate, oo, log
from sympy.physics.wigner import wigner_3j, wigner_9j, clebsch_gordan

from geovac.breit_integrals import breit_ss_radial


# Symbolic radial integration variables
r1, r2, r_outer = symbols('r_1 r_2 r_outer', positive=True)

# Hydrogenic wavefunctions at Z=1 (radial part only; r^2 for volume element added separately)
R1s = 2 * exp(-r1)
R1s_r2 = 2 * exp(-r2)
R2p_r1 = r1 * exp(-r1 / 2) / (2 * sqrt(6))
R2p_r2 = r2 * exp(-r2 / 2) / (2 * sqrt(6))


def compute_bipolar_integral(p1, p2, Q):
    """Compute ∫∫ [rho_1(r_1) * r_1^2] * [rho_2(r_2) * r_2^2] * (r_1^p1 r_2^p2 / r_>^Q) dr_1 dr_2

    for DIRECT: rho_1 = R_1s^2, rho_2 = R_2p^2
    for EXCHANGE: rho_1 = R_1s R_2p, rho_2 = R_2p R_1s

    Returns direct and exchange versions.
    """
    # Direct densities (with r^2 volume element)
    n_1s_1 = R1s**2 * r1**2  # = 4 r_1^2 e^{-2 r_1}
    n_2p_2 = (R2p_r2)**2 * r2**2  # = (1/24) r_2^4 e^{-r_2}
    # Exchange densities
    n_1s2p_1 = R1s * R2p_r1 * r1**2  # = (1/sqrt 6) r_1^3 e^{-3 r_1 / 2}
    n_2p1s_2 = R2p_r2 * R1s_r2 * r2**2  # same form as n_1s2p_1 in variable r2

    # Kernel: r_1^p1 r_2^p2 / r_>^Q, split by region r_1 < r_2 vs r_1 > r_2
    # Region I: r_1 < r_2, r_> = r_2
    # Region II: r_1 > r_2, r_> = r_1
    # Direct integral = ∫_0^∞ dr_2 [n_2p_2 / r_2^Q_if_r1_<_r_2] ∫_0^{r_2} dr_1 [n_1s_1 * r_1^{p1} r_2^{p2}]
    #                 + ∫_0^∞ dr_2 [n_2p_2 * r_2^{p2}] ∫_{r_2}^∞ dr_1 [n_1s_1 * r_1^{p1} / r_1^Q]

    integrand_direct_I = n_1s_1 * n_2p_2 * r1**p1 * r2**p2 / r2**Q  # r_1 < r_2
    integrand_direct_II = n_1s_1 * n_2p_2 * r1**p1 * r2**p2 / r1**Q  # r_1 > r_2

    I1 = integrate(integrate(integrand_direct_I, (r1, 0, r2)), (r2, 0, oo))
    I2 = integrate(integrate(integrand_direct_II, (r1, r2, oo)), (r2, 0, oo))
    direct_val = simplify(I1 + I2)

    integrand_exch_I = n_1s2p_1 * n_2p1s_2 * r1**p1 * r2**p2 / r2**Q
    integrand_exch_II = n_1s2p_1 * n_2p1s_2 * r1**p1 * r2**p2 / r1**Q

    I1e = integrate(integrate(integrand_exch_I, (r1, 0, r2)), (r2, 0, oo))
    I2e = integrate(integrate(integrand_exch_II, (r1, r2, oo)), (r2, 0, oo))
    exch_val = simplify(I1e + I2e)

    return direct_val, exch_val


# Drake M^K integrals (for comparison)
M0_dir = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 0, Z=1)
M0_exch = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 0, Z=1)
M1_dir = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 1, Z=1)
M1_exch = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 1, Z=1)
M2_dir = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, 2, Z=1)
M2_exch = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, 2, Z=1)

print("Drake M^K integrals (numerical floats):")
print(f"  M^0_dir  = {float(M0_dir):+.6e}")
print(f"  M^0_exch = {float(M0_exch):+.6e}")
print(f"  M^1_dir  = {float(M1_dir):+.6e}")
print(f"  M^1_exch = {float(M1_exch):+.6e}")
print(f"  M^2_dir  = {float(M2_dir):+.6e}")
print(f"  M^2_exch = {float(M2_exch):+.6e}")

# Verify the bipolar (0, 2, 2) kernel r_2^2/r_>^5 on direct
print("\n--- Bipolar channel: K=2, (k1=0, k2=2), kernel r_1^0 r_2^2 / r_>^5 ---")
I_d_02_2, I_e_02_2 = compute_bipolar_integral(0, 2, 5)
print(f"  Direct:   {simplify(I_d_02_2)}")
print(f"           = {float(I_d_02_2):+.6e}")
print(f"  Exchange: {simplify(I_e_02_2)}")
print(f"           = {float(I_e_02_2):+.6e}")

print("\n--- Bipolar channel: K=2, (k1=1, k2=1), kernel r_1 r_2 / r_>^5 ---")
I_d_11_2, I_e_11_2 = compute_bipolar_integral(1, 1, 5)
print(f"  Direct:   {simplify(I_d_11_2)}")
print(f"           = {float(I_d_11_2):+.6e}")
print(f"  Exchange: {simplify(I_e_11_2)}")
print(f"           = {float(I_e_11_2):+.6e}")

print("\n--- Bipolar channel: K=2, (k1=2, k2=0), kernel r_1^2 / r_>^5 ---")
I_d_20_2, I_e_20_2 = compute_bipolar_integral(2, 0, 5)
print(f"  Direct:   {simplify(I_d_20_2)}")
print(f"           = {float(I_d_20_2):+.6e}")
print(f"  Exchange: {simplify(I_e_20_2)}")
print(f"           = {float(I_e_20_2):+.6e}")

# Check the identity: is (0,2,2) direct + (2,0,2) direct = M^2_dir (with Drake kernel r_<^2/r_>^5)?
print("\n--- Drake kernel: r_<^2 / r_>^5 (Drake's M^2) ---")
# r_<^2 / r_>^5 in regions:
#   r_1 < r_2: r_1^2 / r_2^5
#   r_1 > r_2: r_2^2 / r_1^5
# We need kernel_p1=2, p2=0, Q=5 if r_1 < r_2 (r_1^2 / r_2^5 matches p1=2, p2=0, divide by r_2^5)
# and kernel_p1=0, p2=2, Q=5 if r_1 > r_2 (r_2^2 / r_1^5 matches p1=0, p2=2, divide by r_1^5)
# This IS NOT a single (p1, p2, Q) form -- it's a piecewise function.
# Let's compute it directly.

def compute_Drake_M_K(K_Drake):
    """Drake's M^K = ∫∫ [dens]_1 [dens]_2 (r_<^K / r_>^{K+3}) r_1^2 r_2^2 dr_1 dr_2"""
    n_1s_1 = R1s**2 * r1**2
    n_2p_2 = (R2p_r2)**2 * r2**2
    n_1s2p_1 = R1s * R2p_r1 * r1**2
    n_2p1s_2 = R2p_r2 * R1s_r2 * r2**2

    # Region I: r_1 < r_2: r_< = r_1, r_> = r_2
    int_dir_I = n_1s_1 * n_2p_2 * r1**K_Drake / r2**(K_Drake + 3)
    int_dir_II = n_1s_1 * n_2p_2 * r2**K_Drake / r1**(K_Drake + 3)
    M_dir = simplify(
        integrate(integrate(int_dir_I, (r1, 0, r2)), (r2, 0, oo)) +
        integrate(integrate(int_dir_II, (r1, r2, oo)), (r2, 0, oo))
    )
    int_exch_I = n_1s2p_1 * n_2p1s_2 * r1**K_Drake / r2**(K_Drake + 3)
    int_exch_II = n_1s2p_1 * n_2p1s_2 * r2**K_Drake / r1**(K_Drake + 3)
    M_exch = simplify(
        integrate(integrate(int_exch_I, (r1, 0, r2)), (r2, 0, oo)) +
        integrate(integrate(int_exch_II, (r1, r2, oo)), (r2, 0, oo))
    )
    return M_dir, M_exch


for K_d in (0, 1, 2):
    mdir, mexch = compute_Drake_M_K(K_d)
    print(f"Drake M^{K_d}_dir (computed)  = {simplify(mdir - float('inf') * 0)}")
    print(f"  float: {float(mdir):+.6e}")
    print(f"Drake M^{K_d}_exch (computed) = {simplify(mexch)}")
    print(f"  float: {float(mexch):+.6e}")
