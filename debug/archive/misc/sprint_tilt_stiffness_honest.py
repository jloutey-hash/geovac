"""
Honest stiffness check: evaluate the WELL curvature at each curve's OWN minimum
(the omega_e-relevant point), not at R_true (which sits on the inner wall and
inflates the curvature).  Reconcile with Paper 17's reported composed omega_e.
"""
from __future__ import annotations
import json
import numpy as np

R_TRUE = 3.015
CM1_TO_HA = 1.0 / 219474.6313702
MU_ME = (7.016004*1.007825)/(7.016004+1.007825) * 1822.888486
K_TRUE = MU_ME * (1405.65*CM1_TO_HA)**2
W_EXP = 1405.65

def omega_e(curv):                       # Ha/bohr^2 -> cm^-1
    return np.sqrt(max(curv, 0)/MU_ME) / CM1_TO_HA

def analyze(label, R, E, w_lo, w_hi):
    """Fit near the minimum, report R_eq, curvature AT the min, and omega_e."""
    c = np.polyfit(R, E, 4); p = np.poly1d(c); dp = p.deriv(1); ddp = p.deriv(2)
    roots = dp.r; roots = roots[np.isreal(roots)].real
    mins = [r for r in roots if ddp(r) > 0 and R.min() < r < R.max()]
    r_eq = min(mins, key=lambda r: p(r)) if mins else R[np.argmin(E)]
    curv_at_min = float(ddp(r_eq))
    curv_at_true = float(ddp(R_TRUE))
    print(f"{label:22s} R_eq={r_eq:.3f}  "
          f"curv@min={curv_at_min:.4f} (w_e={omega_e(curv_at_min):.0f} cm^-1, "
          f"{(omega_e(curv_at_min)/W_EXP-1)*100:+.1f}%)  |  "
          f"curv@R_true={curv_at_true:.4f} (w_e={omega_e(curv_at_true):.0f})")
    return r_eq, curv_at_min, curv_at_true

print(f"experiment: w_e = {W_EXP:.0f} cm^-1,  k_true = {K_TRUE:.4f} Ha/bohr^2\n")

# balanced n=2 (live fine grid)
d = json.load(open('debug/data/pk_amplification_finegrid.json'))
analyze('balanced n=2', np.array(d['grid']), np.array(d['E']), 2.9, 3.6)

# composed l_dependent (from algebraic_pk report), narrow window around each min
CMP_R = np.array([2.500, 2.700, 2.844, 2.989, 3.133, 3.278, 3.422, 3.567, 3.711,
                  3.856, 4.000, 4.500])
def cwin(E_full):
    full_R = np.array([2.000,2.250,2.500,2.700,2.844,2.989,3.133,3.278,3.422,
                       3.567,3.711,3.856,4.000,4.500,5.125,5.750,6.375,7.000])
    idx = [np.where(full_R==r)[0][0] for r in CMP_R]
    return np.array(E_full)[idx]
CMP_L2 = [-8.126168,-8.157216,-8.179985,-8.190854,-8.197908,-8.201961,-8.202484,-8.202702,-8.200522,-8.197152,-8.192090,-8.185661,-8.179401,-8.152712,-8.116517,-8.079466,-8.044413,-8.011630]
CMP_L3 = [-8.133047,-8.166645,-8.193048,-8.207655,-8.217559,-8.224602,-8.228381,-8.231704,-8.232734,-8.232534,-8.230662,-8.227533,-8.224221,-8.207323,-8.180851,-8.151514,-8.121998,-8.093078]
CMP_L4 = [-8.146910,-8.185534,-8.217897,-8.237515,-8.251063,-8.261677,-8.268755,-8.275213,-8.279222,-8.281740,-8.282369,-8.281533,-8.280334,-8.269630,-8.248972,-8.223839,-8.197442,-8.170719]
analyze('composed l=2', CMP_R, cwin(CMP_L2), 3.0, 3.4)
analyze('composed l=3', CMP_R, cwin(CMP_L3), 3.2, 3.7)
analyze('composed l=4', CMP_R, cwin(CMP_L4), 3.4, 4.0)

print("\nPaper 17 reports composed LiH (l-dep PK) w_e = 1471 cm^-1 (+4.6%).")
print("=> reconciles with curv@min, NOT curv@R_true. The '~2x too stiff at R_true'")
print("   is an inner-wall evaluation-point artifact; at its own minimum the well")
print("   stiffness is only mildly high. Corrects the earlier over-statement.")
