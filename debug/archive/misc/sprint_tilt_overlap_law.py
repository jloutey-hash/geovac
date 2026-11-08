"""
Tilt-model test: is the balanced-solver residual TILT eps'(R) proportional to the
block-block bond-overlap slope S'(R)?

LiH balanced bond block = Li-valence(Z_eff=1.0) + H(Z=1): BOTH centers hydrogen-
like Z=1, so the leading bond overlap is the parameter-free 1s-1s Slater overlap
  S(R)   = e^{-R} (1 + R + R^2/3)
  S'(R)  = -e^{-R} (R/3)(1 + R)
Hypothesis A: eps(R)  proportional to  S(R)      => eps'(R) = c S'(R)
Hypothesis B: eps(R)  proportional to  S(R)^2    => eps'(R) = c (2 S S')

eps'(R) = E_c'(R) - E_true'(R) is offset-free.  E_true'(R) from experiment:
harmonic k(R-R_e) (tight window) and spectroscopic Morse (wider), cross-checked.
No reference PES needed (pyscf unavailable on Windows).
"""
from __future__ import annotations
import json
import numpy as np

R_E, E_EXACT, D_E = 3.015, -8.0706, 0.0924
CM1_TO_HA = 1.0 / 219474.6313702
MU_ME = (7.016004*1.007825)/(7.016004+1.007825) * 1822.888486
K_TRUE = MU_ME * (1405.65*CM1_TO_HA)**2          # 0.0659 Ha/bohr^2
A_MORSE = np.sqrt(K_TRUE/(2*D_E))                 # Morse range param

# ---- overlap (1s-1s, Z=1, both centers) ----
def S(R):   return np.exp(-R)*(1 + R + R**2/3)
def Sp(R):  return -np.exp(-R)*(R/3)*(1 + R)

# ---- true well slope (experiment only) ----
def Etrue_p_harm(R):  return K_TRUE*(R - R_E)
def Etrue_p_morse(R):
    x = np.exp(-A_MORSE*(R - R_E))
    return 2*D_E*A_MORSE*(1 - x)*x

# ---- computed balanced LiH n=2 curve (live fine grid, current solver) ----
d = json.load(open('debug/data/pk_amplification_finegrid.json'))
Rg = np.array(d['grid']); Eg = np.array(d['E'])         # E on raw -15.21 scale (offset irrelevant for slopes)

# smooth model of computed PES; analytic derivative
cE = np.polyfit(Rg - R_E, Eg, 5)
pE = np.poly1d(cE); dpE = pE.deriv(1)
def Ec_p(R):  return dpE(R - R_E)

# ---- evaluate on an interior window (avoid grid edges) ----
win = (Rg >= 2.85) & (Rg <= 3.60)
Rw = Rg[win]
epsp_harm  = Ec_p(Rw) - Etrue_p_harm(Rw)
epsp_morse = Ec_p(Rw) - Etrue_p_morse(Rw)
Spw  = Sp(Rw)
SSpw = 2*S(Rw)*Sp(Rw)

def regress(y, x, label):
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    slope, intc = coef
    yhat = A @ coef
    ss_res = np.sum((y - yhat)**2); ss_tot = np.sum((y - y.mean())**2)
    r2 = 1 - ss_res/ss_tot
    print(f"  {label:34s}: c = {slope:+.4f} Ha, intercept = {intc:+.5f}, R^2 = {r2:.4f}")
    return slope, intc, r2

print(f"[ref] k_true={K_TRUE:.4f} Ha/bohr^2, a_Morse={A_MORSE:.3f} bohr^-1")
print(f"[overlap @R_e]  S={S(R_E):.4f}  S'={Sp(R_E):.4f}")
print(f"\nWindow R in [{Rw.min():.2f},{Rw.max():.2f}], {len(Rw)} pts")
print("\n-- Hypothesis A: eps'(R) = c * S'(R)  (eps ~ S) --")
cA_h = regress(epsp_harm,  Spw,  "harmonic E_true'")
cA_m = regress(epsp_morse, Spw,  "Morse E_true'")
print("\n-- Hypothesis B: eps'(R) = c * 2 S S'  (eps ~ S^2) --")
cB_h = regress(epsp_harm,  SSpw, "harmonic E_true'")
cB_m = regress(epsp_morse, SSpw, "Morse E_true'")

# constancy of the pointwise ratio eps'/S' across R (Hypothesis A)
ratio = epsp_harm / Spw
print("\n-- pointwise eps'(R)/S'(R) (Hyp A, harmonic) across the window --")
for R, rr in zip(Rw, ratio):
    print(f"   R={R:.3f}  eps'={Ec_p(R)-Etrue_p_harm(R):+.4f}  S'={Sp(R):+.4f}  ratio={rr:+.4f}")
print(f"   ratio mean={ratio.mean():+.4f} Ha, std={ratio.std():.4f} "
      f"(rel spread {ratio.std()/abs(ratio.mean())*100:.1f}%)")

# reference-free anchor at R_e
epsp_Re = Ec_p(R_E)  # = eps'(R_e) exactly since E_true'(R_e)=0
print(f"\n[anchor] eps'(R_e=3.015) = E_c'(R_e) = {epsp_Re:+.4f} Ha/bohr (reference-free)")
print(f"         predicted by Hyp A: c*S'(R_e) = {cA_h[0]*Sp(R_E):+.4f} "
      f"(c from window fit)")

out = {
    'k_true': K_TRUE, 'a_morse': A_MORSE,
    'S_Re': float(S(R_E)), 'Sp_Re': float(Sp(R_E)),
    'epsp_Re': float(epsp_Re),
    'hypA_harmonic': {'c': float(cA_h[0]), 'r2': float(cA_h[2])},
    'hypA_morse': {'c': float(cA_m[0]), 'r2': float(cA_m[2])},
    'hypB_harmonic': {'c': float(cB_h[0]), 'r2': float(cB_h[2])},
    'hypB_morse': {'c': float(cB_m[0]), 'r2': float(cB_m[2])},
    'ratio_mean': float(ratio.mean()), 'ratio_relspread_pct': float(ratio.std()/abs(ratio.mean())*100),
}
json.dump(out, open('debug/data/tilt_overlap_law.json', 'w'), indent=2)
print("\n[saved] debug/data/tilt_overlap_law.json")
