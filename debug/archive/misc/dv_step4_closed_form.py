"""Sprint 5 DV Step 4: Compute the bipolar radial integrals in closed form
by piecewise direct symbolic integration.

For hydrogenic orbitals at Z=1:
  rho_1s(r) = R_1s^2 * r^2 = 4 r^2 e^{-2r}
  rho_2p(r) = R_2p^2 * r^2 = (1/24) r^4 e^{-r}
  rho_1s2p(r) = R_1s R_2p * r^2 = (1/sqrt 6) r^3 e^{-3r/2}

For channel (k_1, k_2, K=2), the radial kernel is r_1^{k1} r_2^{k2} / r_>^{2K+1} = r_1^{k1} r_2^{k2} / r_>^5.

Direct integral: I_dir(k1, k2, K) = ∫∫ rho_1s(r_1) rho_2p(r_2) r_1^{k1} r_2^{k2} / r_>^5 dr_1 dr_2
Exchange integral: I_exch(k1, k2, K) = ∫∫ rho_1s2p(r_1) rho_1s2p(r_2) r_1^{k1} r_2^{k2} / r_>^5 dr_1 dr_2

We split r_1 < r_2 vs r_1 > r_2 and use sympy to evaluate both pieces exactly.

For each bipolar channel, we then express I in terms of Drake's basis
{M^0_dir, M^1_dir, M^2_dir, M^0_exch, M^1_exch, M^2_exch}.
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
from sympy import Rational, sqrt, Integer, simplify, symbols, exp, integrate, oo, log, nsimplify, Matrix

from geovac.breit_integrals import breit_ss_radial

r1, r2 = symbols('r_1 r_2', positive=True)


def bipolar_int_pieces(p1, p2, Q, dens1, dens2):
    """Compute ∫∫ dens1(r_1) dens2(r_2) r_1^p1 r_2^p2 / r_>^Q dr_1 dr_2
    in the two regions r_1 < r_2 and r_1 > r_2.

    dens1, dens2 are sympy expressions in r1 and r2 respectively (already include r^2).
    """
    # Region I: r_1 < r_2 (so r_> = r_2)
    integrand_I = dens1 * dens2 * r1**p1 * r2**p2 / r2**Q
    val_I = integrate(integrand_I, (r1, 0, r2))
    val_I = integrate(val_I, (r2, 0, oo))
    # Region II: r_1 > r_2 (so r_> = r_1)
    integrand_II = dens1 * dens2 * r1**p1 * r2**p2 / r1**Q
    val_II = integrate(integrand_II, (r1, r2, oo))
    val_II = integrate(val_II, (r2, 0, oo))
    return simplify(val_I + val_II)


# Direct densities (with r^2 volume element included)
dens_1s_sq = 4 * r1**2 * exp(-2 * r1)
dens_2p_sq = Rational(1, 24) * r2**4 * exp(-r2)
# Exchange densities
dens_1s2p_r1 = r1**3 * exp(-3 * r1 / 2) / sqrt(6)
dens_2p1s_r2 = r2**3 * exp(-3 * r2 / 2) / sqrt(6)

# ============================================================
# STEP A: Drake M^K basis, computed symbolically from the piecewise r_<, r_>
# kernel (should match breit_ss_radial at Z=1)
# ============================================================

def drake_M_K(K_d):
    """Drake's M^K = ∫∫ [dens] (r_<^K / r_>^{K+3}) r_1^2 r_2^2 dr_1 dr_2"""
    # r_< = min(r_1, r_2), r_> = max
    # Region I: r_1 < r_2: r_< = r_1
    # Region II: r_1 > r_2: r_< = r_2
    # Both regions have r_>^(K+3) in denominator; kernel = r_<^K = r_(min)^K

    # Direct:
    I_I = integrate(integrate(
        dens_1s_sq * dens_2p_sq * r1**K_d / r2**(K_d + 3),
        (r1, 0, r2)), (r2, 0, oo))
    I_II = integrate(integrate(
        dens_1s_sq * dens_2p_sq * r2**K_d / r1**(K_d + 3),
        (r1, r2, oo)), (r2, 0, oo))
    M_d = simplify(I_I + I_II)
    # Exchange:
    E_I = integrate(integrate(
        dens_1s2p_r1 * dens_2p1s_r2 * r1**K_d / r2**(K_d + 3),
        (r1, 0, r2)), (r2, 0, oo))
    E_II = integrate(integrate(
        dens_1s2p_r1 * dens_2p1s_r2 * r2**K_d / r1**(K_d + 3),
        (r1, r2, oo)), (r2, 0, oo))
    M_e = simplify(E_I + E_II)
    return M_d, M_e


print("Computing Drake M^K_dir, M^K_exch (k=0, 1, 2) from symbolic integration...")
M_drake = {}
for K_d in (0, 1, 2):
    md, me = drake_M_K(K_d)
    M_drake[(K_d, 'dir')] = md
    M_drake[(K_d, 'exch')] = me
    # Cross-check against production module:
    md_prod = breit_ss_radial(1, 0, 2, 1, 1, 0, 2, 1, K_d, Z=1)
    me_prod = breit_ss_radial(1, 0, 2, 1, 2, 1, 1, 0, K_d, Z=1)
    diff_d = simplify(md - md_prod)
    diff_e = simplify(me - me_prod)
    status_d = "MATCH" if diff_d == 0 else f"DIFFER ({float(diff_d):+.2e})"
    status_e = "MATCH" if diff_e == 0 else f"DIFFER ({float(diff_e):+.2e})"
    print(f"  M^{K_d}_dir  = {md}")
    print(f"    vs production: {status_d}")
    print(f"  M^{K_d}_exch = {me}")
    print(f"    vs production: {status_e}")

# ============================================================
# STEP B: Bipolar channel integrals for K=2
# ============================================================
print("\n\n=== Bipolar channel integrals for K=2 ===")
bipolar_K2 = {}
for k1, k2 in [(0, 2), (1, 1), (2, 0)]:
    I_d = bipolar_int_pieces(k1, k2, 5, dens_1s_sq, dens_2p_sq)
    I_e = bipolar_int_pieces(k1, k2, 5, dens_1s2p_r1, dens_2p1s_r2)
    bipolar_K2[(k1, k2, 'dir')] = I_d
    bipolar_K2[(k1, k2, 'exch')] = I_e
    print(f"\n  (k1={k1}, k2={k2}, K=2), kernel r_1^{k1} r_2^{k2} / r_>^5")
    print(f"    Direct:   {simplify(I_d)}")
    print(f"    Exchange: {simplify(I_e)}")

# ============================================================
# STEP C: Express each bipolar integral in the Drake basis
# ============================================================
# The bipolar integrals I_d(k1, k2, K=2), I_e(k1, k2, K=2) should be LINEAR
# combinations of Drake's M^0, M^1, M^2 (direct or exchange respectively).
#
# Specifically, we expect (for the SAME region and SAME density product):
#   I(k1, k2, K=2)_dir = rational combination of M^K_dir for various K
# where the rationals come from the POLYNOMIAL IDENTITY relating
# r_1^{k1} r_2^{k2} / r_>^{2K+1} to r_<^Kd / r_>^{Kd+3} for various Kd.
#
# SPECIFICALLY:
#   Bipolar (k1, k2, K) with Q=2K+1=5:
#     Region I (r_1 < r_2):  r_1^{k1} r_2^{k2} / r_2^5 = r_1^{k1} r_2^{k2-5}
#     Region II (r_1 > r_2): r_1^{k1} r_2^{k2} / r_1^5 = r_1^{k1-5} r_2^{k2}
#
#   Drake M^K (Q=K+3):
#     Region I: r_1^K / r_2^{K+3}
#     Region II: r_2^K / r_1^{K+3}
#
# For Drake's M^K (K=0,1,2) we have exponents:
#   Region I: (K, -K-3)   ∈ {(0,-3), (1,-4), (2,-5)}
#   Region II: (-K-3, K) ∈ {(-3,0), (-4,1), (-5,2)}
#
# For bipolar (k1, k2, K=2) (Q=5):
#   Region I: (k1, k2-5)
#   Region II: (k1-5, k2)
#
# For the bipolar kernel to match a Drake kernel in Region I:
#   Bipolar(0,2): (0, -3) — MATCHES Drake M^0 Region I!
#   Bipolar(1,1): (1, -4) — MATCHES Drake M^1 Region I!
#   Bipolar(2,0): (2, -5) — MATCHES Drake M^2 Region I!
#
# And in Region II:
#   Bipolar(0,2): (-5, 2) — MATCHES Drake M^2 Region II!
#   Bipolar(1,1): (-4, 1) — MATCHES Drake M^1 Region II!
#   Bipolar(2,0): (-3, 0) — MATCHES Drake M^0 Region II!
#
# So the bipolar channels are HYBRIDS: region-I and region-II correspond
# to DIFFERENT Drake M^K.
#
# Bipolar(0, 2): Region I = M^0 Region I, Region II = M^2 Region II
# Bipolar(1, 1): Region I = M^1 Region I, Region II = M^1 Region II
# Bipolar(2, 0): Region I = M^2 Region I, Region II = M^0 Region II

# So we need SEPARATE region integrals for each Drake M^K:
def drake_M_K_piecewise(K_d):
    """Return (I_d_regionI, I_d_regionII, I_e_regionI, I_e_regionII)."""
    I_I = integrate(integrate(
        dens_1s_sq * dens_2p_sq * r1**K_d / r2**(K_d + 3),
        (r1, 0, r2)), (r2, 0, oo))
    I_II = integrate(integrate(
        dens_1s_sq * dens_2p_sq * r2**K_d / r1**(K_d + 3),
        (r1, r2, oo)), (r2, 0, oo))
    E_I = integrate(integrate(
        dens_1s2p_r1 * dens_2p1s_r2 * r1**K_d / r2**(K_d + 3),
        (r1, 0, r2)), (r2, 0, oo))
    E_II = integrate(integrate(
        dens_1s2p_r1 * dens_2p1s_r2 * r2**K_d / r1**(K_d + 3),
        (r1, r2, oo)), (r2, 0, oo))
    return simplify(I_I), simplify(I_II), simplify(E_I), simplify(E_II)


print("\n\n=== Drake M^K piecewise (region I and II contributions) ===")
pieces = {}
for K_d in (0, 1, 2):
    I_I, I_II, E_I, E_II = drake_M_K_piecewise(K_d)
    pieces[(K_d, 'dir', 'I')] = I_I
    pieces[(K_d, 'dir', 'II')] = I_II
    pieces[(K_d, 'exch', 'I')] = E_I
    pieces[(K_d, 'exch', 'II')] = E_II
    print(f"\n  M^{K_d}_dir Region I:  {I_I}")
    print(f"                       = {float(I_I):+.6e}")
    print(f"  M^{K_d}_dir Region II: {I_II}")
    print(f"                       = {float(I_II):+.6e}")
    print(f"  M^{K_d}_exch Region I:  {E_I}")
    print(f"                       = {float(E_I):+.6e}")
    print(f"  M^{K_d}_exch Region II: {E_II}")
    print(f"                       = {float(E_II):+.6e}")

# ============================================================
# STEP D: Verify the bipolar = piecewise Drake identity
# ============================================================
print("\n\n=== Bipolar K=2 channel integrals in Drake basis ===")

# Bipolar(0, 2) Region I = Drake M^0 Region I
# Bipolar(0, 2) Region II = Drake M^2 Region II
bipolar_02_dir_check = simplify(
    bipolar_K2[(0, 2, 'dir')] - (pieces[(0, 'dir', 'I')] + pieces[(2, 'dir', 'II')])
)
bipolar_02_exch_check = simplify(
    bipolar_K2[(0, 2, 'exch')] - (pieces[(0, 'exch', 'I')] + pieces[(2, 'exch', 'II')])
)
print(f"Bipolar(0,2) dir  - [M^0_dir_I + M^2_dir_II]  = {bipolar_02_dir_check}")
print(f"Bipolar(0,2) exch - [M^0_exch_I + M^2_exch_II] = {bipolar_02_exch_check}")

# Bipolar(1,1) = M^1 Region I + M^1 Region II = Full M^1
bipolar_11_dir_check = simplify(
    bipolar_K2[(1, 1, 'dir')] - (pieces[(1, 'dir', 'I')] + pieces[(1, 'dir', 'II')])
)
bipolar_11_exch_check = simplify(
    bipolar_K2[(1, 1, 'exch')] - (pieces[(1, 'exch', 'I')] + pieces[(1, 'exch', 'II')])
)
print(f"Bipolar(1,1) dir  - M^1_dir   = {bipolar_11_dir_check}")
print(f"Bipolar(1,1) exch - M^1_exch  = {bipolar_11_exch_check}")

# Bipolar(2,0) = M^2 Region I + M^0 Region II
bipolar_20_dir_check = simplify(
    bipolar_K2[(2, 0, 'dir')] - (pieces[(2, 'dir', 'I')] + pieces[(0, 'dir', 'II')])
)
bipolar_20_exch_check = simplify(
    bipolar_K2[(2, 0, 'exch')] - (pieces[(2, 'exch', 'I')] + pieces[(0, 'exch', 'II')])
)
print(f"Bipolar(2,0) dir  - [M^2_dir_I + M^0_dir_II]   = {bipolar_20_dir_check}")
print(f"Bipolar(2,0) exch - [M^2_exch_I + M^0_exch_II] = {bipolar_20_exch_check}")

# ============================================================
# STEP E: check density symmetry I <-> II for direct (1s)^2 (2p)^2 vs exchange (1s*2p)^2
# ============================================================
# For the EXCHANGE density (rho_1s2p is the same function of r_1 as of r_2), we have
#     rho(r_1) rho(r_2) symmetric under r_1 <-> r_2
# which means M^K_exch has Region I = Region II (up to r-swap relabeling of variables).
# So M^K_exch = 2 * Region I = 2 * Region II.
print("\n--- Check exchange density r_1<->r_2 symmetry ---")
for K_d in (0, 1, 2):
    ratio = simplify(pieces[(K_d, 'exch', 'II')] / pieces[(K_d, 'exch', 'I')])
    print(f"  M^{K_d}_exch Region II / Region I = {ratio}")
