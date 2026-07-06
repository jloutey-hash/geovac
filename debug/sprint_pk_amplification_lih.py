"""
Sprint: chemistry-error projection decomposition -- LiH tilt/curvature autopsy
==============================================================================
Diagnostic (no solver run; uses current-convention CACHED curves only).

Question: why does the balanced solver's ENERGY converge (1.7%->0.2%, n_max 2->3)
while its GEOMETRY (R_eq) drifts WORSE (7.0%->8.8%)?

Decomposition. For a computed PES E_c(R) with residual eps(R)=E_c(R)-E_true(R):
  - The true minimum is at R_true=3.015 (experiment), where E_true'(R_true)=0.
  - Therefore eps'(R_true) = E_c'(R_true) EXACTLY (reference-free).
  - The computed minimum shifts by  dR = -E_c'(R_true) / E_c''(R_true).
So the geometry error is governed by the residual TILT eps'(R_true) and the well
curvature E_c''(R_true) -- NOT by the residual MAGNITUDE eps(R_true) (= energy error).

We compare:
  eps(R_true)   -- energy error at the true geometry   (the "energy accuracy")
  eps'(R_true)  -- residual tilt at the true geometry   (drives R_eq)
  E_c''(R_true) -- computed well stiffness              (vs experimental k = mu*w_e^2)
across n_max (balanced) and l_max (composed/PK).

Reference well (LiH X^1Sigma+): experimental spectroscopic constants only.
"""
from __future__ import annotations
import json
import numpy as np
from pathlib import Path

# ---------------------------------------------------------------- reference
R_TRUE = 3.015            # bohr, experimental R_e
E_EXACT = -8.0706         # Ha, non-rel Coulomb FCI/CBS total
D_E = 0.0924              # Ha, well depth

# LiH true well curvature from spectroscopy: k = mu * w_e^2  (atomic units)
CM1_TO_HA = 1.0 / 219474.6313702
W_E_CM1 = 1405.65         # cm^-1 (Huber-Herzberg)
U_TO_ME = 1822.888486     # electron masses per u
M_LI7, M_H1 = 7.016004, 1.007825
MU_U = (M_LI7 * M_H1) / (M_LI7 + M_H1)
MU_ME = MU_U * U_TO_ME
W_E_HA = W_E_CM1 * CM1_TO_HA
K_TRUE = MU_ME * W_E_HA ** 2      # Ha / bohr^2
print(f"[ref] mu = {MU_U:.6f} u = {MU_ME:.2f} m_e ; w_e = {W_E_HA:.6e} Ha ;"
      f" k_true = mu*w_e^2 = {K_TRUE:.5f} Ha/bohr^2")

# ---------------------------------------------------------------- cached curves
# BALANCED coupled (PK-free), current Track-CD convention: E_fci are full BO totals
# comparable to E_EXACT.  Source: debug/data/balanced_coupled_lih_nmax{2,3}_analytical.json
BAL_N2_R = np.array([2.0, 2.5, 3.0, 3.015, 3.5, 4.0, 5.0])
BAL_N2_E = np.array([-7.824551, -7.894006, -7.928962, -7.929426,
                     -7.927324, -7.895162, -7.771569])
BAL_N3_R = np.array([2.9, 3.015, 3.1, 3.3, 3.5])
BAL_N3_E = np.array([-8.0493805415, -8.0545240769, -8.0570773414,
                     -8.0591647662, -8.0562794055])

# COMPOSED / PK (l_dependent PK), total BO energy convention (E_composed = E_core +
# V_cross + E_elec + V_NN_bare).  Source: debug/archive/misc/algebraic_pk_lmax_report.md
CMP_R = np.array([2.000, 2.250, 2.500, 2.700, 2.844, 2.989, 3.133, 3.278, 3.422,
                  3.567, 3.711, 3.856, 4.000, 4.500, 5.125, 5.750, 6.375, 7.000])
CMP_L2_E = np.array([-8.126168, -8.157216, -8.179985, -8.190854, -8.197908,
                     -8.201961, -8.202484, -8.202702, -8.200522, -8.197152,
                     -8.192090, -8.185661, -8.179401, -8.152712, -8.116517,
                     -8.079466, -8.044413, -8.011630])
CMP_L3_E = np.array([-8.133047, -8.166645, -8.193048, -8.207655, -8.217559,
                     -8.224602, -8.228381, -8.231704, -8.232734, -8.232534,
                     -8.230662, -8.227533, -8.224221, -8.207323, -8.180851,
                     -8.151514, -8.121998, -8.093078])
CMP_L4_E = np.array([-8.146910, -8.185534, -8.217897, -8.237515, -8.251063,
                     -8.261677, -8.268755, -8.275213, -8.279222, -8.281740,
                     -8.282369, -8.281533, -8.280334, -8.269630, -8.248972,
                     -8.223839, -8.197442, -8.170719])


def autopsy(label: str, R: np.ndarray, E: np.ndarray, window: float = 0.9):
    """Local polynomial autopsy of a PES about R_TRUE."""
    # window of points around the minimum for a stable local fit
    m = np.abs(R - R_TRUE) <= window
    Rw, Ew = R[m], E[m]
    order = min(4, len(Rw) - 1)
    x = Rw - R_TRUE
    c = np.polyfit(x, Ew, order)          # highest power first, in x=R-R_TRUE
    p = np.poly1d(c)
    dp = p.deriv(1)
    ddp = p.deriv(2)
    # computed minimum (root of derivative nearest R_TRUE, within data range)
    roots = dp.r
    roots = roots[np.isreal(roots)].real
    roots = roots[(roots > (Rw.min() - R_TRUE - 0.2)) & (roots < (Rw.max() - R_TRUE + 0.2))]
    if len(roots):
        # pick the one that is a minimum (ddp>0) and lowest energy
        mins = [r for r in roots if ddp(r) > 0]
        r_eq_x = min(mins, key=lambda r: p(r)) if mins else roots[np.argmin([p(r) for r in roots])]
    else:
        r_eq_x = x[np.argmin(Ew)]
    R_eq = r_eq_x + R_TRUE

    E_at_true = float(p(0.0))
    slope_true = float(dp(0.0))            # = eps'(R_true)  (since E_true'(R_true)=0)
    curv_true = float(ddp(0.0))            # = E_c''(R_true)
    eps_true = E_at_true - E_EXACT         # energy error at true geometry
    dR = R_eq - R_TRUE
    dR_from_tilt = -slope_true / curv_true if curv_true else float('nan')

    print(f"\n=== {label} ===  (order {order} fit, {len(Rw)} pts)")
    print(f"  energy error  eps(R_true)      = {eps_true:+.4f} Ha "
          f"({abs(eps_true)/abs(E_EXACT)*100:.2f}% of |E|)")
    print(f"  residual tilt eps'(R_true)     = {slope_true:+.4f} Ha/bohr")
    print(f"  computed curv E_c''(R_true)    = {curv_true:+.4f} Ha/bohr^2"
          f"   (k_true={K_TRUE:.4f}, ratio {curv_true/K_TRUE:.2f}x)")
    print(f"  R_eq                           = {R_eq:.3f} bohr "
          f"(error {abs(dR)/R_TRUE*100:.1f}%)")
    print(f"  dR observed = {dR:+.3f} ;  -eps'/E_c'' = {dR_from_tilt:+.3f}  (consistency)")
    print(f"  hypothetical dR if well had true stiffness: "
          f"-eps'/k_true = {-slope_true/K_TRUE:+.3f} bohr")
    return {
        'label': label, 'n_pts_fit': int(len(Rw)), 'order': int(order),
        'eps_true_Ha': eps_true, 'eps_pct': abs(eps_true)/abs(E_EXACT)*100,
        'tilt_Ha_per_bohr': slope_true, 'curv_Ha_per_bohr2': curv_true,
        'curv_over_ktrue': curv_true/K_TRUE, 'R_eq': R_eq,
        'R_eq_err_pct': abs(dR)/R_TRUE*100, 'dR': dR,
        'dR_from_tilt': dR_from_tilt, 'dR_if_true_stiffness': -slope_true/K_TRUE,
    }


rows = []
print("\n" + "=" * 70)
print("BALANCED coupled (PK-free): energy converges, does the TILT?")
print("=" * 70)
rows.append(autopsy('balanced n_max=2', BAL_N2_R, BAL_N2_E, window=1.2))
rows.append(autopsy('balanced n_max=3', BAL_N3_R, BAL_N3_E, window=0.6))

print("\n" + "=" * 70)
print("COMPOSED / PK (l_dependent): energy AND geometry both drift")
print("=" * 70)
rows.append(autopsy('composed l_max=2', CMP_R, CMP_L2_E, window=1.0))
rows.append(autopsy('composed l_max=3', CMP_R, CMP_L3_E, window=1.0))
rows.append(autopsy('composed l_max=4', CMP_R, CMP_L4_E, window=1.0))

out = {
    'reference': {'R_true': R_TRUE, 'E_exact': E_EXACT, 'D_e': D_E,
                  'k_true_Ha_bohr2': K_TRUE, 'w_e_cm1': W_E_CM1, 'mu_u': MU_U},
    'rows': rows,
}
Path('debug/data').mkdir(parents=True, exist_ok=True)
with open('debug/data/pk_amplification_lih.json', 'w') as f:
    json.dump(out, f, indent=2)
print("\n[saved] debug/data/pk_amplification_lih.json")

# ------ headline trend: does the tilt heal with basis while energy does? ------
print("\n" + "=" * 70)
print("HEADLINE TREND")
print("=" * 70)
b2, b3 = rows[0], rows[1]
print(f"  balanced energy error : {b2['eps_pct']:.2f}% -> {b3['eps_pct']:.2f}%"
      f"  (factor {b2['eps_pct']/b3['eps_pct']:.1f} better)")
print(f"  balanced tilt eps'    : {b2['tilt_Ha_per_bohr']:+.4f} -> "
      f"{b3['tilt_Ha_per_bohr']:+.4f} Ha/bohr  (does NOT heal)")
print(f"  balanced R_eq error   : {b2['R_eq_err_pct']:.1f}% -> {b3['R_eq_err_pct']:.1f}%"
      f"  (drifts WORSE)")
print(f"  computed stiffness    : {b2['curv_over_ktrue']:.2f}x / "
      f"{b3['curv_over_ktrue']:.2f}x  the true k (well too stiff both)")
