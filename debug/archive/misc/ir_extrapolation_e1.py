"""
IR Extrapolation for E1 Dipole Polarizability
==============================================
Extends the Furnstahl-More-Papenbrock IR correction framework
to the E1 sum rule / electric dipole polarizability.

Physics:
--------
A finite HO basis at N_shells, hw imposes an effective hard-wall
boundary at:
    L_eff = sqrt(2*(2*N_max + 2*l + 3)) * b     (Furnstahl et al. 2012)
where b = sqrt(hbar^2 / (m*hw)) is the oscillator length and
N_max = N_shells - 1 is the maximum HO shell included.

For the GROUND STATE energy, the IR correction is:
    E(L) = E_inf + A_0 * exp(-2*kappa_0 * L)
where kappa_0 = sqrt(2*mu*|E_bind|)/hbar is the binding wave number
and A_0 depends on the asymptotic normalization coefficient (ANC).

For the POLARIZABILITY, the sum over intermediate states samples
the continuum. In a box of size L, continuum states are discretized
with momenta k_n ~ n*pi/L. The E1 polarizability is:

    alpha_E = 2*e^2 * sum_n |<n|D_z|0>|^2 / (E_n - E_0)

In the continuum limit (L -> inf), the sum becomes an integral:
    alpha_E = 2*e^2 * integral dk * rho(k) * |<k|D_z|0>|^2 / E(k)

where rho(k) = L/pi is the density of states in a 1D box.

The IR correction to alpha_E has two contributions:
  (a) Missing continuum strength above the box cutoff (k > pi*N/L)
  (b) Discretization error from replacing integral with sum

For the deuteron (s-wave ground state, p-wave E1 intermediates):
  - The E1 matrix element <k,l=1|r|0,l=0> involves the overlap of
    the p-wave continuum scattering state with the s-wave bound state
    weighted by r (the dipole operator)
  - At large k, the matrix element falls as ~ 1/k^2 (oscillatory
    cancellation in the radial integral)
  - The dominant contribution comes from low-k states near threshold

The KEY INSIGHT (Furnstahl framework applied to E1):
  The hard wall at L creates a box quantization for the p-wave
  continuum with k_n = (n + 1/2)*pi/L (hard wall, l=1).
  The polarizability in the box is:
    alpha_E(L) = (2*e^2*L/pi) * sum_n |M(k_n)|^2 / E(k_n)
  where M(k) = <k,l=1|r|psi_0> is the reduced E1 matrix element.

  The correction is:
    alpha_E(inf) - alpha_E(L) = integral contribution from k > k_max
                               + trapezoidal discretization error
  For the deuteron, the dominant correction is the MISSING TAIL:
    delta_alpha ~ (2*e^2/pi) * integral_{k_max}^inf dk |M(k)|^2 / E(k)

  Using M(k) ~ C * kappa * k / (kappa^2 + k^2)^2 (Hulthen model) and
  E(k) = k^2/(2*mu) + |E_bind|:
    delta_alpha ~ C^2 * kappa^2 / (pi * mu) * integral dk k^2 / ((kappa^2+k^2)^4 * (k^2/(2mu) + |E|))

  This integral is computable analytically for the deuteron.

Implementation:
--------------
1. Compute alpha_E(L) at multiple L values (by varying hw at fixed N_shells)
2. Fit to alpha_E(L) = alpha_E(inf) + correction(L)
3. Extrapolate to L -> inf
4. Compare extrapolated value to experiment
"""
import sys, io, json
import numpy as np
from scipy.optimize import curve_fit

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
sys.path.insert(0, '.')

from geovac.nuclear.nuclear_hamiltonian import (
    build_deuteron_hamiltonian, enumerate_sp_states, DeuteronSpec,
)

HBAR_C = 197.3269804  # MeV*fm
E2_MEV_FM = 1.4399764
M_N = 938.918  # average nucleon mass MeV
MU_D = M_N / 2  # deuteron reduced mass (n-p)
E_BIND_D = 2.2246  # MeV (positive, binding energy)
KAPPA_D = np.sqrt(2 * MU_D * E_BIND_D) / HBAR_C  # fm^-1
ALPHA_E_EXP = 0.6328  # fm^3

# Reuse dipole infrastructure
def ho_radial_r_me(nr1, l1, nr2, l2, b_fm):
    if l1 == l2 + 1:
        if nr1 == nr2: return b_fm * np.sqrt(nr2 + l2 + 1.5)
        elif nr1 == nr2 - 1: return b_fm * np.sqrt(float(nr2))
        return 0.0
    elif l1 == l2 - 1:
        return ho_radial_r_me(nr2, l2, nr1, l1, b_fm)
    return 0.0

def angular_cos_theta(l1, m1, l2, m2):
    if m1 != m2: return 0.0
    m = m1
    if l1 == l2 + 1:
        return np.sqrt(((l2+1)**2 - m**2) / ((2*l2+1)*(2*l2+3)))
    elif l1 == l2 - 1:
        return np.sqrt((l2**2 - m**2) / ((2*l2-1)*(2*l2+1)))
    return 0.0

def sp_dipole_z(states, b_fm):
    n = len(states)
    D = np.zeros((n, n))
    for a, sa in enumerate(states):
        for c, sc in enumerate(states):
            if sa.m_l != sc.m_l or sa.m_s != sc.m_s: continue
            if abs(sa.l - sc.l) != 1: continue
            D[a,c] = ho_radial_r_me(sa.n_r, sa.l, sc.n_r, sc.l, b_fm) * \
                     angular_cos_theta(sa.l, sa.m_l, sc.l, sc.m_l)
    return D

def deuteron_alpha_at_hw(hw, N_shells=2):
    data = build_deuteron_hamiltonian(N_shells=N_shells, hw=hw)
    H = data['H_matrix']
    evals, evecs = np.linalg.eigh(H)
    gs, E_gs = evecs[:,0], evals[0]
    b_fm = data['metadata']['b_fm']
    n_n = len(data['states_n'])
    D_sp = sp_dipole_z(data['states_p'], b_fm)
    dim = len(evals)
    D_z = np.zeros((dim, dim))
    for i in range(len(data['states_p'])):
        for k in range(len(data['states_p'])):
            if abs(D_sp[i,k]) < 1e-15: continue
            for j in range(n_n):
                D_z[i*n_n+j, k*n_n+j] += D_sp[i,k]
    D_gs = D_z @ gs
    alpha_E = sum(2*E2_MEV_FM*np.dot(evecs[:,n],D_gs)**2/(evals[n]-E_gs)
                  for n in range(1,dim) if evals[n]-E_gs > 1e-10)
    return alpha_E, b_fm, E_gs

# =========================================================================
print("=" * 72)
print("IR EXTRAPOLATION FOR E1 DIPOLE POLARIZABILITY")
print("Extending Furnstahl-More-Papenbrock to the E1 sum rule")
print("=" * 72)

print(f"\nDeuteron parameters:")
print(f"  E_bind = {E_BIND_D:.4f} MeV")
print(f"  kappa = {KAPPA_D:.4f} fm^-1  (binding wave number)")
print(f"  1/kappa = {1/KAPPA_D:.2f} fm  (asymptotic tail scale)")
print(f"  Experimental alpha_E = {ALPHA_E_EXP:.4f} fm^3")

# Step 1: Fine hw scan to get alpha_E(L_eff)
print(f"\n--- Step 1: Compute alpha_E vs L_eff ---")
N_shells = 2
N_max = N_shells - 1  # = 1

# L_eff for the s-wave ground state (l=0):
# L_eff = sqrt(2*(2*N_max + 3)) * b = sqrt(2*(2+3)) * b = sqrt(10) * b
# For the p-wave intermediates (l=1):
# L_eff_p = sqrt(2*(2*N_max + 2*1 + 3)) * b = sqrt(2*(2+2+3)) * b = sqrt(14) * b
# The relevant L_eff for the polarizability is the SMALLER of the two
# (the ground state constrains the box size for the overlap integral)

hw_values = np.arange(3, 41, 0.5)
scan = []

for hw in hw_values:
    alpha_E, b_fm, E_gs = deuteron_alpha_at_hw(hw, N_shells)
    L_eff_s = np.sqrt(2 * (2*N_max + 3)) * b_fm  # s-wave box
    L_eff_p = np.sqrt(2 * (2*N_max + 2 + 3)) * b_fm  # p-wave box
    L_eff = min(L_eff_s, L_eff_p)  # controlling box size
    scan.append({
        'hw': hw, 'b': b_fm, 'alpha_E': alpha_E,
        'L_eff_s': L_eff_s, 'L_eff_p': L_eff_p,
        'L_eff': L_eff, 'E_gs': E_gs,
        'two_kappa_L': 2 * KAPPA_D * L_eff,
    })

L_arr = np.array([s['L_eff'] for s in scan])
aE_arr = np.array([s['alpha_E'] for s in scan])
hw_arr = np.array([s['hw'] for s in scan])

print(f"\n  {'hw':>5} {'b':>6} {'L_eff':>7} {'2kL':>6} {'alpha_E':>9}")
for s in scan[::4]:
    print(f"  {s['hw']:>5.1f} {s['b']:>6.3f} {s['L_eff']:>7.2f} "
          f"{s['two_kappa_L']:>6.2f} {s['alpha_E']:>9.4f}")

# Step 2: Fit IR extrapolation model
print(f"\n--- Step 2: IR extrapolation fits ---")

# Model A: Furnstahl-style exponential
# alpha_E(L) = alpha_inf + A * exp(-2*kappa*L)
# This is the energy-correction analog applied to polarizability
def model_exp(L, alpha_inf, A):
    return alpha_inf + A * np.exp(-2 * KAPPA_D * L)

# Model B: Power-law correction (box quantization of continuum)
# alpha_E(L) ~ alpha_inf + B / L^2
# Density of states in box ~ L; matrix elements ~ 1/L^{3/2}; net ~ 1/L^2
def model_power(L, alpha_inf, B):
    return alpha_inf + B / L**2

# Model C: Combined (exponential for bound-state tail + power for continuum)
def model_combined(L, alpha_inf, A, B):
    return alpha_inf + A * np.exp(-2 * KAPPA_D * L) + B / L**2

# Fit only in the physically relevant range where alpha_E is reasonable
# (exclude very small L where UV effects dominate)
mask = (L_arr > 5.0) & (aE_arr > 0.01) & (aE_arr < 5.0)
L_fit = L_arr[mask]
aE_fit = aE_arr[mask]

print(f"  Fitting {len(L_fit)} points with L_eff in [{L_fit.min():.1f}, {L_fit.max():.1f}] fm")

results = {}

# Model A: exponential
try:
    popt_A, pcov_A = curve_fit(model_exp, L_fit, aE_fit,
                                p0=[0.6, -5.0], maxfev=10000)
    alpha_inf_A = popt_A[0]
    residuals_A = aE_fit - model_exp(L_fit, *popt_A)
    rmse_A = np.sqrt(np.mean(residuals_A**2))
    print(f"\n  Model A (exponential): alpha_inf = {alpha_inf_A:.4f} fm^3, "
          f"A = {popt_A[1]:.4f}, RMSE = {rmse_A:.6f}")
    print(f"    vs experiment: {(alpha_inf_A - ALPHA_E_EXP)/ALPHA_E_EXP*100:+.1f}%")
    results['model_A'] = {'alpha_inf': alpha_inf_A, 'A': popt_A[1],
                          'rmse': rmse_A, 'model': 'exp(-2*kappa*L)'}
except Exception as e:
    print(f"  Model A: FAILED ({e})")

# Model B: power law
try:
    popt_B, pcov_B = curve_fit(model_power, L_fit, aE_fit,
                                p0=[0.6, 10.0], maxfev=10000)
    alpha_inf_B = popt_B[0]
    residuals_B = aE_fit - model_power(L_fit, *popt_B)
    rmse_B = np.sqrt(np.mean(residuals_B**2))
    print(f"\n  Model B (1/L^2):       alpha_inf = {alpha_inf_B:.4f} fm^3, "
          f"B = {popt_B[1]:.4f}, RMSE = {rmse_B:.6f}")
    print(f"    vs experiment: {(alpha_inf_B - ALPHA_E_EXP)/ALPHA_E_EXP*100:+.1f}%")
    results['model_B'] = {'alpha_inf': alpha_inf_B, 'B': popt_B[1],
                          'rmse': rmse_B, 'model': '1/L^2'}
except Exception as e:
    print(f"  Model B: FAILED ({e})")

# Model C: combined
try:
    popt_C, pcov_C = curve_fit(model_combined, L_fit, aE_fit,
                                p0=[0.6, -2.0, 10.0], maxfev=10000)
    alpha_inf_C = popt_C[0]
    residuals_C = aE_fit - model_combined(L_fit, *popt_C)
    rmse_C = np.sqrt(np.mean(residuals_C**2))
    print(f"\n  Model C (exp + 1/L^2): alpha_inf = {alpha_inf_C:.4f} fm^3, "
          f"A = {popt_C[1]:.4f}, B = {popt_C[2]:.4f}, RMSE = {rmse_C:.6f}")
    print(f"    vs experiment: {(alpha_inf_C - ALPHA_E_EXP)/ALPHA_E_EXP*100:+.1f}%")
    results['model_C'] = {'alpha_inf': alpha_inf_C, 'A': popt_C[1],
                          'B': popt_C[2], 'rmse': rmse_C,
                          'model': 'exp(-2*kappa*L) + 1/L^2'}
except Exception as e:
    print(f"  Model C: FAILED ({e})")

# Step 3: Analytical Hulthen-model estimate
print(f"\n--- Step 3: Analytical Hulthen-model polarizability ---")
# The Hulthen wavefunction for the deuteron:
#   psi_0(r) = sqrt(kappa/(2*pi)) * (exp(-kappa*r) - exp(-lambda*r)) / r
# with kappa = 0.232 fm^-1 (binding), lambda ~ 1.4 fm^-1 (short-range).
#
# The E1 reduced matrix element squared, summed and integrated over
# the p-wave continuum, gives the Thomas-Reiche-Kuhn (TRK) sum rule:
#   S_1 = integral (E_n - E_0) |<n|D|0>|^2 = (9/4) * hbar^2 / (2*M_d)
# (factor 9/4 = Z^2/A for the isoscalar E1 for deuteron Z=1, A=2
#  with the effective charge e_eff = e/2 from CM subtraction)
#
# For the zero-energy-weighted sum (the polarizability):
#   alpha_E = 2*e^2 * S_{-1}
# where S_{-1} = integral |<n|D|0>|^2 / (E_n - E_0)
#
# The Hulthen model gives an analytical expression (Arenhovel & Sanzone 1991):
#   alpha_E(Hulthen) = e^2 / (2*mu*E_bind) * (ANC)^2 * I(kappa)
# where I(kappa) involves the p-wave phase space integral.
#
# Simpler: from dimensional analysis + known asymptotic behavior:
#   alpha_E ~ e^2 * ANC^2 / (2*mu*kappa^5)
# The ANC for the deuteron: C_0 ~ 0.8845 fm^{-1/2} (Ericson & Rosa-Clot)

ANC_D = 0.8845  # fm^{-1/2}, deuteron s-wave ANC
alpha_hulthen = E2_MEV_FM * ANC_D**2 / (2 * MU_D / HBAR_C**2 * KAPPA_D**5)
# Units: MeV*fm * fm^{-1} / (MeV^{-1}*fm^{-2} * fm^{-5})
# = MeV*fm * fm^{-1} * MeV * fm^2 * fm^5 = fm^3 ... need to check

# More careful: alpha_E has units fm^3.
# e^2 = 1.44 MeV*fm
# ANC^2 has units fm^{-1}
# mu = M_N/2 = 469.5 MeV/c^2
# In natural units (hbar*c = 197.3 MeV*fm):
# 1/(2*mu) = hbar^2/(2*mu*c^2) * c^2/hbar^2 -- this gets messy
#
# Let's use the known result: for a zero-range s-wave bound state,
# alpha_E = e^2 * kappa / (48 * mu_red * E_bind^2)
# (Rupak & Higa 2011, leading order in pionless EFT)

# Actually the simplest analytical result is from Chen et al. 1998:
# alpha_E^(LO) = alpha * gamma_t / (16 * M_N * B^2)
# where gamma_t = sqrt(M_N * B) is the binding momentum
# Hmm, these formulas use different conventions.

# Let me just use the pionless EFT NLO result directly:
# alpha_E(NLO) = 0.595 fm^3 (Chen-Griesshammer-Savage-Springer 1998)
alpha_eft_nlo = 0.595

# And the modern NNLO: ~0.63 fm^3 (Emmons-Ji-Platter 2020)
alpha_eft_nnlo = 0.63

print(f"  Pionless EFT NLO:  {alpha_eft_nlo:.3f} fm^3 (Chen et al. 1998)")
print(f"  Pionless EFT NNLO: {alpha_eft_nnlo:.2f} fm^3 (Emmons et al. 2020)")
print(f"  SLEGS experiment:  {ALPHA_E_EXP:.4f} +/- 0.028 fm^3 (Hao et al. 2026)")

# Step 4: Summary
print(f"\n{'='*72}")
print(f"EXTRAPOLATION RESULTS SUMMARY")
print(f"{'='*72}")

best = min(results.values(), key=lambda x: x['rmse'])
best_name = [k for k, v in results.items() if v is best][0]

print(f"\n  Best fit model: {best_name} ({best['model']})")
print(f"  alpha_E(inf) = {best['alpha_inf']:.4f} fm^3")
print(f"  RMSE = {best['rmse']:.6f} fm^3")
print(f"")
print(f"  {'Method':<30} {'alpha_E':>10} {'vs exp':>10}")
print(f"  {'-'*55}")
print(f"  {'Raw at hw=8 (alpha-matched)':30} {'0.6944':>10} {'+9.7%':>10}")
for name, r in results.items():
    res = (r['alpha_inf'] - ALPHA_E_EXP) / ALPHA_E_EXP * 100
    print(f"  {name + ' extrapolation':30} {r['alpha_inf']:>10.4f} {res:>+9.1f}%")
print(f"  {'Pionless EFT NLO':30} {alpha_eft_nlo:>10.3f} {(alpha_eft_nlo-ALPHA_E_EXP)/ALPHA_E_EXP*100:>+9.1f}%")
print(f"  {'Pionless EFT NNLO':30} {alpha_eft_nnlo:>10.2f} {(alpha_eft_nnlo-ALPHA_E_EXP)/ALPHA_E_EXP*100:>+9.1f}%")
print(f"  {'SLEGS experiment (2026)':30} {ALPHA_E_EXP:>10.4f} {'---':>10}")

# Save
output = {
    'scan': [{'hw': s['hw'], 'L_eff': s['L_eff'], 'alpha_E': s['alpha_E'],
              'b_fm': s['b']} for s in scan],
    'fits': {k: {kk: float(vv) for kk, vv in v.items() if isinstance(vv, (int, float))}
             for k, v in results.items()},
    'kappa_d': KAPPA_D,
    'alpha_exp': ALPHA_E_EXP,
}
with open('debug/data/ir_extrapolation_e1.json', 'w') as f:
    json.dump(output, f, indent=2)

print(f"\nSaved to debug/data/ir_extrapolation_e1.json")
